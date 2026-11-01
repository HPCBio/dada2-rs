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

use crate::containers::{B, Bi, BirthType, Comparison};
use crate::pval::get_pA;

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

/// Rebuild every cluster's membership from the pass's assignments: the centre
/// first, then members in input order. Each member's `comp` becomes its
/// comparison with its own centre, which the post-run p-value pass and the
/// outputs read.
pub fn finalize(b: &mut B, assign: &[usize], best: &[Option<Best>]) {
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
        if let Some(bst) = &best[j] {
            debug_assert_eq!(bst.cluster, k);
            b.raws[j].comp = bst.comp.clone();
        }
    }
}
