//! Experimental single-pass partition (issue #323).
//!
//! An alternative to DADA2's divisive loop, after UNOISE3 (Edgar 2016,
//! doi:10.1101/081257): visit the uniques once, most abundant first. Each one
//! either becomes a new centre or joins the existing centre expected to produce
//! the most reads of it. Nothing is reassigned afterwards.
//!
//! The birth test is DADA2's, with two choices settled on #323:
//!
//! - The expected reads of a raw are `λ × the centre's own reads`, not the
//!   cluster's running total, so a raw's fate does not depend on which raws
//!   joined the cluster before it.
//! - A raw's best centre is the one with the most expected reads, the quantity
//!   the divisive shuffle assigns by.
//!
//! The MIN_FOLD, MIN_HAMMING, MIN_ABUNDANCE and singleton gates apply as in
//! the divisive bud, and the abundance p-value is Bonferroni-corrected by the
//! number of uniques as there.
//!
//! Comparisons are `b_compare_parallel`'s, run from each new centre exactly as
//! the divisive loop runs them after a split, so the screen, alignment and λ
//! are identical between the two partitions. Raws more abundant than a new
//! centre are already decided, and greedy mode skips them.

use std::time::Instant;

use rayon::prelude::*;

use crate::cluster::{
    BatchTiming, CompCost, CompareTiming, NOT_DEFERRED, align_deferred, deferred_pairs,
};
use crate::containers::{B, Bi, BirthType, Comparison, Raw, Sub};
use crate::nwalign::{
    AlignBuffers, AlignParams, GAP_GLYPH, SubOrDefer, sub_new_or_defer, sub_new_with_buf,
};
use crate::pval::{compute_lambda, get_pA};

/// The centre a raw would join: most expected reads seen so far.
#[derive(Clone)]
pub struct Best {
    /// `λ × centre reads`.
    pub e: f64,
    pub cluster: usize,
    pub comp: Comparison,
}

/// Visiting order: decreasing reads, then input index, so ties go to the
/// first-seen raw as in `B::assign_center`.
pub fn order(b: &B) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..b.raws.len()).collect();
    idx.sort_by(|&x, &y| b.raws[y].reads.cmp(&b.raws[x].reads).then(x.cmp(&y)));
    idx
}

/// Fold cluster `k`'s stored comparisons into `best`. A tie keeps the earlier
/// cluster.
///
/// `b_compare_parallel` stores a comparison only if `λ × total reads` beats
/// the raw's best expected reads so far, which every comparison that could
/// raise `best` does, so nothing that matters here is filtered out.
pub fn update_best(best: &mut [Option<Best>], b: &B, k: usize) {
    let center = b.clusters[k]
        .center
        .expect("update_best: cluster has no centre");
    let center_reads = b.raws[center].reads as f64;
    for comp in &b.clusters[k].comp {
        let e = comp.lambda * center_reads;
        let slot = &mut best[comp.index as usize];
        if slot.as_ref().is_none_or(|cur| e > cur.e) {
            *slot = Some(Best {
                e,
                cluster: k,
                comp: comp.clone(),
            });
        }
    }
}

/// Greedy locks for new cluster `k`, as `b_p_update` sets them once for a newly
/// born cluster: its centre, and every raw whose best centre is now `k` and
/// which `k`'s centre alone is expected to produce more reads of than were
/// observed. A locked raw is skipped by later comparisons, as under the
/// divisive loop. Cluster 0's locks come from the initial `b_p_update`.
pub fn lock_new_members(b: &mut B, best: &[Option<Best>], k: usize) {
    let center = b.clusters[k]
        .center
        .expect("lock_new_members: cluster has no centre");
    b.raws[center].lock = true;
    for ci in 0..b.clusters[k].comp.len() {
        let j = b.clusters[k].comp[ci].index as usize;
        if let Some(bst) = &best[j]
            && bst.cluster == k
            && bst.e > b.raws[j].reads as f64
        {
            b.raws[j].lock = true;
        }
    }
}

/// Gates and thresholds for the birth test.
#[derive(Clone, Copy)]
pub struct BirthTest {
    pub omega_a: f64,
    pub omega_p: f64,
    pub min_fold: f64,
    pub min_hamming: u32,
    pub min_abund: u32,
    pub detect_singletons: bool,
}

/// A raw that passed the birth test.
pub struct Birth {
    pub kind: BirthType,
    pub pval: f64,
}

/// Does raw `j`, whose best centre is `best`, become a new centre?
pub fn birth_test(b: &B, j: usize, best: &Best, t: &BirthTest) -> Option<Birth> {
    let raw = &b.raws[j];
    let center = b.clusters[best.cluster]
        .center
        .expect("birth_test: cluster has no centre");
    let center_reads = b.raws[center].reads;
    if raw.reads < t.min_abund || best.comp.hamming < t.min_hamming {
        return None;
    }
    if t.min_fold > 1.0 && (raw.reads as f64) < t.min_fold * best.e {
        return None;
    }
    let p = get_pA(
        raw.reads,
        raw.prior,
        best.comp.lambda,
        best.comp.hamming,
        center_reads,
        t.detect_singletons,
    );
    let p_a = p * b.raws.len() as f64;
    if p_a < t.omega_a {
        return Some(Birth {
            kind: BirthType::Abundance,
            pval: p_a,
        });
    }
    if raw.prior && p < t.omega_p {
        return Some(Birth {
            kind: BirthType::Prior,
            pval: p,
        });
    }
    None
}

/// Create a cluster centred on raw `j`, born from `best.cluster`. Membership is
/// assigned later by [`finalize`]; returns the new cluster's index.
pub fn add_center(b: &mut B, j: usize, best: &Best, birth: &Birth) -> usize {
    let mut bi = Bi::new(b.raws.len() as u32);
    bi.birth_type = birth.kind.clone();
    bi.birth_from = best.cluster as u32;
    bi.birth_pval = birth.pval;
    bi.birth_e = best.e;
    bi.birth_fold = b.raws[j].reads as f64 / best.e.max(f64::MIN_POSITIVE);
    bi.birth_comp = best.comp.clone();
    bi.center = Some(j);
    bi.seq = b.raws[j].seq.clone();
    b.add_cluster(bi)
}

/// Create a cluster centred on raw `j` under the UNOISE rule. The rule has no
/// p-value, fold or parent, so the birth fields keep their defaults apart from
/// the type.
pub fn add_unoise_center(b: &mut B, j: usize) -> usize {
    let mut bi = Bi::new(b.raws.len() as u32);
    bi.birth_type = BirthType::Abundance;
    bi.center = Some(j);
    bi.seq = b.raws[j].seq.clone();
    b.add_cluster(bi)
}

/// Rebuild every cluster's membership from the pass's assignments: the centre
/// first, then members in input order. Each member's `comp` becomes its
/// comparison with its own centre, which the post-run p-value pass and the
/// outputs read.
pub fn finalize(b: &mut B, assign: &[usize], comps: &[Option<Comparison>]) {
    for bi in &mut b.clusters {
        bi.raws.clear();
        bi.reads = 0;
    }
    let centers: Vec<usize> = b
        .clusters
        .iter()
        .map(|bi| bi.center.expect("finalize: cluster has no centre"))
        .collect();
    for (k, &c) in centers.iter().enumerate() {
        b.bi_add_raw(k, c);
        // Set directly: `b_compare_parallel` records the centre's comparison
        // with itself in `self_` always, but stores it only if it passes the
        // store filter.
        b.raws[c].comp = Comparison {
            i: k as u32,
            index: c as u32,
            lambda: b.clusters[k].self_,
            hamming: 0,
        };
    }
    for (j, &k) in assign.iter().enumerate() {
        if centers[k] == j {
            continue;
        }
        b.bi_add_raw(k, j);
        // `None` keeps the raw's comparison with cluster 0's centre, from the
        // initial compare: the UNOISE rule's unassigned raws.
        if let Some(comp) = &comps[j] {
            debug_assert_eq!(comp.i as usize, k);
            b.raws[j].comp = comp.clone();
        }
    }
}

// ---------------------------------------------------------------------------
// UNOISE birth rule (`--birth-rule unoise`)
// ---------------------------------------------------------------------------

/// UNOISE's α (Edgar 2016; vsearch `--unoise_alpha` default).
pub const UNOISE_ALPHA: f64 = 2.0;

/// Maximum skew a member at distance `d` may have: β(d) = 1/2^(αd+1).
pub fn unoise_beta(d: u32) -> f64 {
    0.5f64.powf(UNOISE_ALPHA * d as f64 + 1.0)
}

/// Gap columns inside the alignment: reference bases opposite a query gap,
/// and query bases between two aligned reference bases. Leading and trailing
/// gaps on either sequence are not counted, as under vsearch's `--iddef 2`.
pub fn internal_gap_columns(sub: &Sub) -> u32 {
    internal_columns(sub).1
}

/// The alignment's columns between its first and last aligned reference
/// bases, and how many of those are gaps. End gaps on either sequence are
/// outside this span, as under vsearch's `--iddef 2`.
pub fn internal_columns(sub: &Sub) -> (u32, u32) {
    let mapped: Vec<(usize, u16)> = sub
        .map
        .iter()
        .enumerate()
        .filter(|&(_, &q)| q != GAP_GLYPH)
        .map(|(i, &q)| (i, q))
        .collect();
    let (Some(&(first, _)), Some(&(last, _))) = (mapped.first(), mapped.last()) else {
        return (0, 0);
    };
    // Deletions: unmapped reference bases between the first and last mapped.
    let deletions = (first..=last).filter(|&i| sub.map[i] == GAP_GLYPH).count() as u32;
    // Insertions: query bases skipped between consecutive mapped positions.
    let insertions: u32 = mapped
        .windows(2)
        .map(|w| (w[1].1 - w[0].1 - 1) as u32)
        .sum();
    (
        (last - first + 1) as u32 + insertions,
        deletions + insertions,
    )
}

/// UNOISE's distance d: substitutions plus internal gap columns.
pub fn unoise_d(sub: &Sub) -> u32 {
    sub.nsubs() as u32 + internal_gap_columns(sub)
}

/// `--min-abund` under `--birth-rule unoise` when the flag is omitted: vsearch
/// `--cluster_unoise`'s `--minsize` default. The UNOISE rule has no abundance
/// test, so this is its only birth gate; at 1, every unique no centre absorbs
/// becomes a centre, and each new centre is aligned against every unique.
pub const UNOISE_MIN_ABUND_DEFAULT: u32 = 8;

/// Identity at which UNOISE3 counts a read to a ZOTU: USEARCH `-otutab`'s
/// default, used with ZOTUs as with 97% OTUs.
pub const UNOISE_COUNT_ID: f64 = 0.97;

/// A raw's standing against the centres compared with it so far.
#[derive(Clone)]
pub struct UnoiseBest {
    /// Smallest d among centres the raw qualifies for (skew ≤ β(d)). Decides
    /// births only.
    pub qualified: Option<(u32, usize, Comparison)>,
    /// The centre with the highest identity of at least [`UNOISE_COUNT_ID`],
    /// at any skew: (identity, cluster, comparison). Decides which ASV the
    /// raw's reads count toward, as UNOISE3's `-otutab` does.
    pub nearest: Option<(f64, usize, Comparison)>,
}

/// Identity as vsearch `--iddef 2` defines it: matches over the alignment's
/// columns between its first and last aligned reference bases.
pub fn unoise_identity(sub: &Sub) -> f64 {
    let (cols, gaps) = internal_columns(sub);
    if cols == 0 {
        return 0.0;
    }
    let matches = cols - gaps - sub.nsubs() as u32;
    matches as f64 / cols as f64
}

/// One pair's result in [`unoise_compare`]: `(d, identity, lambda, nsubs)`,
/// `None` when the screen rejected it.
type Aligned = Option<(u32, f64, f64, u32)>;

/// Align new centre `k` against every raw and record two things for each.
///
/// - Births: the qualifying centre with the smallest d (skew ≤ β(d)). Only a
///   raw with at most half the centre's reads can qualify, since β(d) ≤ 1/2.
/// - Counting: the centre with the highest identity of at least
///   [`UNOISE_COUNT_ID`], at any skew, as UNOISE3's `-otutab` maps every read
///   to its closest ZOTU.
///
/// Ties keep the earlier, more abundant centre in both, as `-otutab` breaks
/// ties by ZOTU file order. Also sets the centre's `self_`. Counts pairs into
/// `b.nalign` and screened-out pairs into `b.nshroud`, as `b_compare` does.
///
/// Returns the same [`CompareTiming`] as [`b_compare_parallel`], so the
/// compare phase is attributed under `--metrics-json` (#328). With `measure`,
/// the worker time is split into screen, DP, traceback and post-alignment.
/// `stored` counts the pairs that updated a raw's best centre.
///
/// Under `--batch-align` the pairs bound for the vectorized kernel are
/// aligned in a second pass, in same-length batches (#345). Unlike
/// `b_compare`, which batches only in its parallel path, this runs at any
/// thread count, since there is no separate serial path.
pub fn unoise_compare(
    b: &mut B,
    k: usize,
    best: &mut [UnoiseBest],
    err_mat: &[f64],
    ncol: usize,
    params: &AlignParams,
    measure: bool,
) -> CompareTiming {
    let t_setup = Instant::now();
    let center = b.clusters[k]
        .center
        .expect("unoise_compare: cluster has no centre");
    let center_reads = b.raws[center].reads;
    let use_quals = b.use_quals;
    let raws = &b.raws;
    let make_buf = || {
        if measure {
            AlignBuffers::measuring()
        } else {
            AlignBuffers::new()
        }
    };
    let score = |raw: &Raw, sub: &Sub| {
        let lambda = compute_lambda(raw, Some(sub), err_mat, ncol, use_quals);
        (
            unoise_d(sub),
            unoise_identity(sub),
            lambda,
            sub.nsubs() as u32,
        )
    };
    // `--batch-align` (#345): the screen pass leaves pairs bound for the
    // vectorized kernel unaligned, and `align_deferred` aligns them 16 at a
    // time and fills in their results, as in `b_compare_parallel`.
    let batch = params.batch_align && !params.screen_audit;
    let compare_one = |buf: &mut AlignBuffers, j: usize| {
        let t0 = measure.then(Instant::now);
        let raw = &raws[j];
        let (aligned, deferred) = if batch {
            match sub_new_or_defer(&raws[center], raw, params, buf, None) {
                SubOrDefer::Defer => (None, true),
                SubOrDefer::Done(sub) => (sub.map(|sub| score(raw, &sub)), false),
            }
        } else {
            let sub = sub_new_with_buf(&raws[center], raw, params, buf);
            (sub.map(|sub| score(raw, &sub)), false)
        };
        // The screen resets the per-pair timers, so a pair it rejected, or
        // deferred, reports only its screen time.
        let cost = match t0 {
            Some(t) => CompCost {
                total: t.elapsed().as_nanos() as u64,
                screen: buf.last_screen_nanos,
                dp: buf.last_dp_nanos,
                post: buf.last_post_nanos,
                traceback: buf.last_traceback_nanos,
            },
            None => CompCost::default(),
        };
        let len = if deferred {
            raw.seq.len() as u32
        } else {
            NOT_DEFERRED
        };
        ((aligned, cost), len)
    };
    let setup = t_setup.elapsed();
    let t_map = Instant::now();
    let mut results: Vec<(Aligned, CompCost)> = Vec::new();
    let mut defer_len: Vec<u32> = Vec::new();
    (0..raws.len())
        .into_par_iter()
        .map_init(make_buf, compare_one)
        .unzip_into_vecs(&mut results, &mut defer_len);
    let screen_pass = t_map.elapsed();
    let mut batch_timing = batch.then(|| BatchTiming {
        screen_pass,
        ..BatchTiming::default()
    });
    let batch_cost = match batch_timing.as_mut() {
        Some(bt) => {
            let pending = deferred_pairs(&defer_len);
            let scan = t_map.elapsed() - screen_pass;
            let cost = align_deferred(
                &mut results,
                pending,
                raws,
                center,
                params,
                measure,
                bt,
                |raw, sub, slot: &mut (Aligned, CompCost)| slot.0 = Some(score(raw, sub)),
            );
            bt.prep += scan;
            cost
        }
        None => CompCost::default(),
    };
    drop(defer_len);
    let map = t_map.elapsed();
    let t_serial = Instant::now();
    let mut cost = batch_cost;
    let (mut considered, mut aligned_n, mut stored) = (0u64, 0u64, 0u64);
    for (j, &(aligned, item)) in results.iter().enumerate() {
        cost = cost.add(item);
        considered += 1;
        let Some((d, id, lambda, nsubs)) = aligned else {
            b.nshroud += 1;
            continue;
        };
        aligned_n += 1;
        if j == center {
            b.clusters[k].self_ = lambda;
            continue;
        }
        let comp = Comparison {
            i: k as u32,
            index: j as u32,
            lambda,
            hamming: nsubs,
        };
        if id >= UNOISE_COUNT_ID {
            let slot = &mut best[j].nearest;
            if slot.as_ref().is_none_or(|&(cur, _, _)| id > cur) {
                *slot = Some((id, k, comp.clone()));
                stored += 1;
            }
        }
        let skew = b.raws[j].reads as f64 / center_reads as f64;
        if skew > unoise_beta(d) {
            continue;
        }
        let slot = &mut best[j].qualified;
        if slot.as_ref().is_none_or(|&(cur, _, _)| d < cur) {
            *slot = Some((d, k, comp));
            stored += 1;
        }
    }
    let serial = t_serial.elapsed();
    let t_free = Instant::now();
    drop(results);
    let free = t_free.elapsed();
    b.nalign += considered;
    CompareTiming {
        map,
        serial,
        free,
        setup,
        stored,
        batch: batch_timing,
        screened: considered,
        aligned: aligned_n,
        ..cost.timing()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nwalign::al2subs;

    fn enc(s: &str) -> Vec<u8> {
        s.bytes()
            .map(|c| match c {
                b'A' => 1,
                b'C' => 2,
                b'G' => 3,
                b'T' => 4,
                _ => b'-',
            })
            .collect()
    }

    fn d(al0: &str, al1: &str) -> (usize, u32) {
        let sub = al2subs(&enc(al0), &enc(al1));
        (sub.nsubs(), internal_gap_columns(&sub))
    }

    #[test]
    fn beta_matches_the_preprint() {
        assert_eq!(unoise_beta(1), 1.0 / 8.0);
        assert_eq!(unoise_beta(2), 1.0 / 32.0);
        assert_eq!(unoise_beta(3), 1.0 / 128.0);
    }

    #[test]
    fn internal_gaps_are_counted_by_column() {
        assert_eq!(d("ACGTACGTAC", "ACGTACGTAC"), (0, 0));
        assert_eq!(d("ACGTACGTAC", "ACGAACGTAC"), (1, 0));
        // A 3-base deletion is three columns, not one event.
        assert_eq!(d("ACGTACGTAC", "ACG---GTAC"), (0, 3));
        // A 2-base insertion in the query.
        assert_eq!(d("ACGT--ACGTAC", "ACGTGGACGTAC"), (0, 2));
        // Substitution plus internal gap.
        assert_eq!(d("ACGTACGTAC", "ACG-ACGTTC"), (1, 1));
    }

    #[test]
    fn identity_is_over_internal_columns() {
        let id = |a: &str, b: &str| unoise_identity(&al2subs(&enc(a), &enc(b)));
        assert_eq!(id("ACGTACGTAC", "ACGAACGTAC"), 0.9);
        // End gaps are outside the span: 6 columns, all matches.
        assert_eq!(id("ACGTACGTAC", "--GTACGT--"), 1.0);
        // A 2-base insertion adds two columns: 10 of 12 match.
        assert_eq!(id("ACGT--ACGTAC", "ACGTGGACGTAC"), 10.0 / 12.0);
    }

    #[test]
    fn end_gaps_are_not_counted() {
        // Leading and trailing gaps on either sequence.
        assert_eq!(d("ACGTACGTAC", "--GTACGT--"), (0, 0));
        assert_eq!(d("--GTACGTAC--", "ACGTACGTACGT"), (0, 0));
        // An internal gap still counts beside end gaps.
        assert_eq!(d("ACGTACGTAC", "--GT-CGT--"), (0, 1));
    }
}
