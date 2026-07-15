//! Cluster operations: compare, shuffle, bud.
//!
//! Ports `cluster.cpp`, excluding R/Rcpp and RcppParallel wrappers.
//! Parallel comparisons use Rayon in place of RcppParallel.

use std::sync::OnceLock;

use rayon::prelude::*;

use crate::containers::{B, Bi, BirthType, Comparison};
use crate::nwalign::{AlignBuffers, AlignParams, sub_new_with_buf};
use crate::pval::compute_lambda;

/// Maximum chunk size for the parallel raw-compare loop in `b_compare_parallel`
/// (passed to rayon's `with_max_len`). Default `32`. Overridable for tuning via
/// the `DADA2_RS_PAR_GRAIN` env var (must be > 0; invalid values fall back to
/// the default). Read once per process and cached. Undocumented in `--help`:
/// this is a tuning knob, not user-facing config.
fn par_max_len() -> usize {
    static VALUE: OnceLock<usize> = OnceLock::new();
    *VALUE.get_or_init(|| {
        std::env::var("DADA2_RS_PAR_GRAIN")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or(32)
    })
}

// ---------------------------------------------------------------------------
// b_compare  (serial)
// ---------------------------------------------------------------------------

/// Align every Raw to the center of cluster `i`, compute lambda under `err_mat`,
/// and store comparisons that could attract a Raw to this cluster.
///
/// Serial version. Equivalent to C++ `b_compare`.
pub fn b_compare(
    b: &mut B,
    i: usize,
    err_mat: &[f64],
    ncol: usize,
    params: &AlignParams,
    greedy: bool,
    verbose: bool,
) {
    let center_idx = b.clusters[i]
        .center
        .expect("b_compare: cluster has no center");
    let center_reads = b.raws[center_idx].reads;

    if verbose {
        eprint!("C{i}LU:");
    }

    let mut buf = AlignBuffers::new();
    for index in 0..b.raws.len() {
        let skip = greedy && (b.raws[index].reads > center_reads || b.raws[index].lock);

        let sub = if skip {
            None
        } else {
            let s = sub_new_with_buf(&b.raws[center_idx], &b.raws[index], params, &mut buf);
            b.nalign += 1;
            if s.is_none() {
                b.nshroud += 1;
            }
            s
        };

        let lambda = compute_lambda(&b.raws[index], sub.as_ref(), err_mat, ncol, b.use_quals);

        if index == center_idx {
            b.clusters[i].self_ = lambda;
        }

        let total_reads = b.reads as f64;
        if lambda * total_reads > b.raws[index].e_minmax {
            let new_e = lambda * center_reads as f64;
            if new_e > b.raws[index].e_minmax {
                b.raws[index].e_minmax = new_e;
            }
            let update_raw = i == 0 || index == center_idx;
            let comp = Comparison {
                i: i as u32,
                index: index as u32,
                lambda,
                hamming: sub.as_ref().map_or(0, |s| s.nsubs() as u32),
            };
            b.clusters[i].comp.push(comp.clone());
            if update_raw {
                b.raws[index].comp = comp;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// b_compare_parallel
// ---------------------------------------------------------------------------

/// Parallel version of `b_compare` using Rayon.
/// Equivalent to C++ `b_compare_parallel`.
/// Returns `(map, serial, busy)` durations: `map` is the parallel-pass wall
/// time, `serial` the post-processing store-loop wall time, and `busy` the
/// summed per-item compute time across all worker threads. `busy / (map ×
/// nthreads)` is the map's parallel efficiency — if well below 1, threads idle
/// inside the parallel region (tail load-imbalance); near 1 with low OS-level
/// utilization points to memory-bandwidth stalls instead. Returned (not
/// global-accumulated) so it stays correct under nested per-sample parallelism.
pub fn b_compare_parallel(
    b: &mut B,
    i: usize,
    err_mat: &[f64],
    ncol: usize,
    params: &AlignParams,
    greedy: bool,
    measure: bool,
) -> (
    std::time::Duration,
    std::time::Duration,
    std::time::Duration,
) {
    let center_idx = b.clusters[i]
        .center
        .expect("b_compare_parallel: cluster has no center");
    let center_reads = b.raws[center_idx].reads;
    let nraw = b.raws.len();
    let use_quals = b.use_quals;
    let t_map = std::time::Instant::now();

    // Read-only parallel pass over raws.
    //
    // Each result carries:
    //   lambda  — error-model probability
    //   hamming — substitution count (u32::MAX = kmer-shrouded / no alignment)
    //   skipped — true when greedy mode skipped this raw entirely
    let raws = b.raws.as_slice();
    // Load balancing: raws are abundance-sorted so per-task cost is skewed
    // (high-abundance raws trigger full alignment; low-abundance get kmer-
    // screened or greedy-skipped). Limiting the maximum task size gives
    // rayon's work-stealing more splits to rebalance across workers, at a
    // small fixed per-task overhead. `map_init` caches AlignBuffers per
    // worker, so buffer reuse still holds across many small tasks.
    //
    // 32 was chosen empirically: ~7% faster than the default on an 8-core
    // box with F3D0 (nraw≈2000); larger thread counts on skewed workloads
    // benefit more. Smaller values (16, 8) were not meaningfully better at
    // 8 threads, but may help at higher thread counts — override via the
    // `DADA2_RS_PAR_GRAIN` env var to tune for your workload.
    // Per-item compute time (4th tuple field, nanos) is summed after collect to
    // get total worker-busy time without cross-thread atomic contention.
    let comps: Vec<(f64, u32, bool, u64)> = (0..nraw)
        .into_par_iter()
        .with_max_len(par_max_len())
        .map_init(AlignBuffers::new, |buf, index| {
            // Per-item timing only under `measure` (verbose) — keeps the hot
            // alignment loop allocation/Instant-free in production runs.
            let t0 = measure.then(std::time::Instant::now);
            let raw = &raws[index];
            let (lambda, hamming, skipped) = if greedy && (raw.reads > center_reads || raw.lock) {
                let lambda = compute_lambda(raw, None, err_mat, ncol, use_quals);
                (lambda, u32::MAX, true)
            } else {
                let sub = sub_new_with_buf(&raws[center_idx], raw, params, buf);
                let lambda = compute_lambda(raw, sub.as_ref(), err_mat, ncol, use_quals);
                let hamming = sub.as_ref().map_or(u32::MAX, |s| s.nsubs() as u32);
                (lambda, hamming, false)
            };
            let nanos = t0.map_or(0, |t| t.elapsed().as_nanos() as u64);
            (lambda, hamming, skipped, nanos)
        })
        .collect();
    let map_dur = t_map.elapsed();
    let busy_dur = std::time::Duration::from_nanos(comps.iter().map(|c| c.3).sum());

    // Serial post-processing: selectively store comparisons.
    let t_serial = std::time::Instant::now();
    let total_reads = b.reads as f64;
    for (index, (lambda, hamming, skipped, _busy)) in comps.into_iter().enumerate() {
        // Match serial b_compare counting: only count non-skipped raws.
        if !skipped {
            b.nalign += 1;
            if hamming == u32::MAX {
                b.nshroud += 1;
            }
        }

        if index == center_idx {
            b.clusters[i].self_ = lambda;
        }

        if lambda * total_reads > b.raws[index].e_minmax {
            let new_e = lambda * center_reads as f64;
            if new_e > b.raws[index].e_minmax {
                b.raws[index].e_minmax = new_e;
            }
            let update_raw = i == 0 || index == center_idx;
            let comp = Comparison {
                i: i as u32,
                index: index as u32,
                lambda,
                hamming: if hamming == u32::MAX { 0 } else { hamming },
            };
            b.clusters[i].comp.push(comp.clone());
            if update_raw {
                b.raws[index].comp = comp;
            }
        }
    }
    (map_dur, t_serial.elapsed(), busy_dur)
}

// ---------------------------------------------------------------------------
// b_shuffle2
// ---------------------------------------------------------------------------

/// Diagnostics from one [`b_shuffle2`] call, used to quantify the rescan
/// redundancy that gates any future incremental best-cluster tracking (see the
/// "reduce work, don't parallelize" conclusion in docs/results.md). Purely
/// observational — collecting it changes no partition behaviour.
#[derive(Clone, Copy, Default)]
pub struct ShuffleStats {
    /// Raws relocated to a different cluster (summed over the loop's iterations
    /// for `b_shuffle_converge`; a single call's moves for `b_shuffle2`).
    pub moves: usize,
    /// Comparisons scanned to (re)build the best-cluster map. For the serial
    /// scan this is Σ over clusters of comp-vec length per call; for the
    /// incremental driver it is the realised work (one build + the per-iteration
    /// recomputes), so the reduction is directly measurable.
    pub comps_scanned: usize,
    /// Cluster count at this call (context for the scan cost).
    pub nclusters: usize,
    /// Move-pass iterations that relocated no raws (pure convergence checks).
    pub zero_move_calls: usize,
    /// Move-pass iterations run.
    pub calls: usize,
}

impl ShuffleStats {
    /// Whether any raw moved — the loop-termination signal callers use.
    pub fn shuffled(&self) -> bool {
        self.moves > 0
    }
}

/// Move each Raw to the cluster that maximises its expected read count.
/// The center of a cluster may not be reassigned.
/// Returns a [`ShuffleStats`] whose `shuffled()` is true iff any Raws moved.
/// Equivalent to C++ `b_shuffle2`.
pub fn b_shuffle2(b: &mut B) -> ShuffleStats {
    let nraw = b.raws.len();

    // Initialise best-E and best-comparison trackers from cluster 0.
    // During initialisation b_compare is always run on cluster 0 first, so
    // its comp vec contains an entry for every raw (in raw-index order).
    let mut emax: Vec<f64> = vec![f64::NEG_INFINITY; nraw];
    let mut compmax: Vec<Comparison> = vec![Comparison::default(); nraw];

    let mut comps_scanned = b.clusters[0].comp.len();
    let c0_reads = b.clusters[0].reads as f64;
    for comp in &b.clusters[0].comp {
        let idx = comp.index as usize;
        emax[idx] = comp.lambda * c0_reads;
        compmax[idx] = comp.clone();
    }

    // Scan remaining clusters for better matches.
    for ci in 1..b.clusters.len() {
        comps_scanned += b.clusters[ci].comp.len();
        let ci_reads = b.clusters[ci].reads as f64;
        for comp in &b.clusters[ci].comp {
            let idx = comp.index as usize;
            let e = comp.lambda * ci_reads;
            if e > emax[idx] {
                emax[idx] = e;
                compmax[idx] = comp.clone();
            }
        }
    }

    // Move raws to their best cluster.
    // Iterate backwards because bi_pop_raw uses swap_remove.
    let mut moves = 0usize;
    for ci in 0..b.clusters.len() {
        let mut r = b.clusters[ci].raws.len();
        while r > 0 {
            r -= 1;
            let raw_idx = b.clusters[ci].raws[r];
            let best_ci = compmax[raw_idx].i as usize;
            if best_ci != ci {
                if b.clusters[ci].center == Some(raw_idx) {
                    // Centers may not leave their cluster.
                    continue;
                }
                b.bi_pop_raw(ci, r);
                b.bi_add_raw(best_ci, raw_idx);
                b.raws[raw_idx].comp = compmax[raw_idx].clone();
                moves += 1;
            }
        }
    }
    ShuffleStats {
        moves,
        comps_scanned,
        nclusters: b.clusters.len(),
        zero_move_calls: usize::from(moves == 0),
        calls: 1,
    }
}

// ---------------------------------------------------------------------------
// b_shuffle_converge — incremental shuffle-to-convergence
// ---------------------------------------------------------------------------

/// A pointer to one stored comparison: `clusters[ci].comp[off]`. The candidate
/// index holds these 8-byte references instead of copying each comparison's
/// `lambda`/`hamming` (issue #85) — the values live once, in the per-cluster
/// `comp` vecs (which are also the shuffle's contiguous initial-build source).
/// Valid for the whole run because a cluster's `comp` vec is written once, by
/// `b_compare` at the cluster's creation, and never mutated afterwards.
#[derive(Clone, Copy)]
pub struct CandRef {
    pub ci: u32,
    pub off: u32,
}

/// Per-raw candidate index, maintained across the whole `run_dada`: each
/// cluster's stored comparisons are appended once, when the cluster is created
/// (see [`index_add_cluster`]). This lets a shuffle recompute a raw's best
/// cluster from just its own candidates, without rescanning every cluster — and
/// crucially it is built incrementally (O(new comps) per bud) rather than
/// rebuilt per shuffle loop, which is what makes the incremental driver a net
/// win rather than a wash. Entries are [`CandRef`] pointers into the per-cluster
/// `comp` vecs, not copies, so the comparison data is stored once (issue #85).
pub type CandIndex = Vec<Vec<CandRef>>;

/// Append cluster `ci`'s stored comparisons to the per-raw candidate index as
/// pointers. Must be called once per cluster, immediately after `b_compare`/
/// `b_compare_parallel` populates `clusters[ci].comp`, and in ascending `ci`
/// order (cluster 0 first, then each bud's new cluster). Ascending order means
/// each raw's candidate list stays cluster-ascending, so a strict-`>` scan keeps
/// the lowest cluster index on ties — matching the serial scan's tie-break.
pub fn index_add_cluster(index: &mut CandIndex, b: &B, ci: usize) {
    for (off, comp) in b.clusters[ci].comp.iter().enumerate() {
        index[comp.index as usize].push(CandRef {
            ci: ci as u32,
            off: off as u32,
        });
    }
}

/// Raw's best cluster over its candidate list at the clusters' current reads,
/// dereferencing each [`CandRef`] into the owning cluster's `comp` vec.
/// Ascending-ci order + strict `>` reproduces the serial lowest-ci tie-break.
fn best_from_cands(cands: &[CandRef], raw: usize, clusters: &[Bi]) -> Comparison {
    let mut best_e = f64::NEG_INFINITY;
    let mut best = Comparison::default();
    for c in cands {
        let comp = &clusters[c.ci as usize].comp[c.off as usize];
        let e = comp.lambda * clusters[c.ci as usize].reads as f64;
        if e > best_e {
            best_e = e;
            best = Comparison {
                i: c.ci,
                index: raw as u32,
                lambda: comp.lambda,
                hamming: comp.hamming,
            };
        }
    }
    best
}

/// Incremental equivalent of looping [`b_shuffle2`] to convergence (or
/// `max_shuffle`). Rebuilds `compmax` once from the persistent candidate
/// `index` (one pass over the comps), then after each move pass only recomputes
/// raws whose candidate clusters' read counts changed — instead of rescanning
/// every cluster every iteration. Byte-identical to the serial loop: `compmax`
/// is kept exactly equal to what a full rebuild at the current reads would
/// produce (same max, same lowest-ci tie-break), so the move decisions are
/// identical. `comps_scanned` reports the realised scan work (initial build +
/// per-iteration recomputes) so the reduction against the serial baseline is
/// directly measurable.
pub fn b_shuffle_converge(b: &mut B, index: &CandIndex, max_shuffle: usize) -> ShuffleStats {
    let nraw = b.raws.len();

    // Initial build: every raw's true best at the current reads. Done the
    // serial way — a contiguous, cache-friendly scan of the per-cluster comp
    // vecs (cluster-major), NOT the raw-major inverted index. This is the bulk
    // of the work each loop, and sequential access here is far cheaper per
    // comparison than the index's scattered reads (which is what made a
    // fully-index-based build a net wall loss despite fewer comparisons). The
    // index is reserved for the reconcile, where the touched-raw volume is
    // small. Byte-identical to b_shuffle2's build: ascending ci + strict `>`
    // keeps the lowest-ci max.
    let mut compmax = vec![Comparison::default(); nraw];
    let mut emax = vec![f64::NEG_INFINITY; nraw];
    let mut comps_scanned = 0usize;
    for bi in &b.clusters {
        let ci_reads = bi.reads as f64;
        for comp in &bi.comp {
            let idx = comp.index as usize;
            let e = comp.lambda * ci_reads;
            if e > emax[idx] {
                emax[idx] = e;
                compmax[idx] = comp.clone();
            }
        }
        comps_scanned += bi.comp.len();
    }
    // emax was only needed to build compmax; the reconcile recomputes affected
    // raws from the index, so it is not carried forward.

    // Reads the map is currently consistent with (for dirty detection).
    let mut reads_used: Vec<u32> = b.clusters.iter().map(|c| c.reads).collect();

    // Reused scratch for the affected-raw set (avoid per-iteration nraw alloc).
    let mut in_affected = vec![false; nraw];
    let mut affected: Vec<usize> = Vec::new();

    let mut total_moves = 0usize;
    let mut nshuffle = 0usize;
    let mut zero_move_calls = 0usize;
    loop {
        // Move pass — identical to b_shuffle2's, using the current compmax.
        let mut moves = 0usize;
        for ci in 0..b.clusters.len() {
            let mut r = b.clusters[ci].raws.len();
            while r > 0 {
                r -= 1;
                let raw_idx = b.clusters[ci].raws[r];
                let best_ci = compmax[raw_idx].i as usize;
                if best_ci != ci {
                    if b.clusters[ci].center == Some(raw_idx) {
                        continue;
                    }
                    b.bi_pop_raw(ci, r);
                    b.bi_add_raw(best_ci, raw_idx);
                    b.raws[raw_idx].comp = compmax[raw_idx].clone();
                    moves += 1;
                }
            }
        }
        total_moves += moves;
        if moves == 0 {
            zero_move_calls += 1;
        }
        nshuffle += 1;
        if moves == 0 || nshuffle >= max_shuffle {
            break;
        }

        // Reconcile: any raw with a comp in a cluster whose reads changed may
        // have a new best. Collect those raws (dedup) and recompute each from
        // its candidates — provably correct at the current reads regardless of
        // increase/decrease, avoiding stale-max ordering hazards.
        for (ci, ru) in reads_used.iter_mut().enumerate() {
            if b.clusters[ci].reads != *ru {
                for comp in &b.clusters[ci].comp {
                    let raw = comp.index as usize;
                    if !in_affected[raw] {
                        in_affected[raw] = true;
                        affected.push(raw);
                    }
                }
                *ru = b.clusters[ci].reads;
            }
        }
        for &raw in &affected {
            compmax[raw] = best_from_cands(&index[raw], raw, &b.clusters);
            comps_scanned += index[raw].len();
            in_affected[raw] = false;
        }
        affected.clear();
    }

    ShuffleStats {
        moves: total_moves,
        comps_scanned,
        nclusters: b.clusters.len(),
        zero_move_calls,
        calls: nshuffle,
    }
}

// ---------------------------------------------------------------------------
// b_bud
// ---------------------------------------------------------------------------

/// Find the Raw with the smallest abundance p-value. If it passes the
/// Bonferroni-corrected significance threshold, pop it into a new cluster.
///
/// Returns `Some(new_cluster_idx)` on a successful bud, `None` otherwise.
/// Equivalent to C++ `b_bud`.
pub fn b_bud(
    b: &mut B,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
    verbose: bool,
    raws_scanned: &mut u64,
) -> Option<usize> {
    let nraw = b.raws.len() as f64;
    let init_center = b.clusters[0]
        .center
        .expect("b_bud: cluster 0 has no center");

    // Redundancy accounting (observational only): the full min-p scan visits
    // every non-center raw across all clusters each call, yet at most one buds.
    // Recording the scan cost bounds the headroom an incremental min-p
    // structure (e.g. a p-ordered queue) could reclaim. Knowable up front, so
    // it stays correct across the function's several early returns.
    *raws_scanned = b
        .clusters
        .iter()
        .map(|bi| bi.raws.len().saturating_sub(1) as u64)
        .sum();

    // (cluster_idx, position_r, raw_index) for non-prior and prior minimums.
    let mut mini: Option<(usize, usize, usize)> = None;
    let mut mini_prior: Option<(usize, usize, usize)> = None;
    let mut min_p = b.raws[init_center].p;
    let mut min_p_prior = b.raws[init_center].p;
    let mut min_reads = b.raws[init_center].reads;
    let mut min_reads_prior = b.raws[init_center].reads;

    for ci in 0..b.clusters.len() {
        // r=1: skip position 0, which is the center of the cluster.
        for r in 1..b.clusters[ci].raws.len() {
            let raw_idx = b.clusters[ci].raws[r];
            let raw = &b.raws[raw_idx];

            if raw.reads < min_abund {
                continue;
            }
            let hamming = raw.comp.hamming;
            let lambda = raw.comp.lambda;

            if hamming < min_hamming {
                continue;
            }
            let fold_ok = min_fold <= 1.0
                || raw.reads as f64 >= min_fold * lambda * b.clusters[ci].reads as f64;
            if !fold_ok {
                continue;
            }

            // Non-prior minimum p-value.
            if raw.p < min_p || (raw.p == min_p && raw.reads > min_reads) {
                mini = Some((ci, r, raw_idx));
                min_p = raw.p;
                min_reads = raw.reads;
            }
            // Prior-sequence minimum p-value.
            if raw.prior
                && (raw.p < min_p_prior || (raw.p == min_p_prior && raw.reads > min_reads_prior))
            {
                mini_prior = Some((ci, r, raw_idx));
                min_p_prior = raw.p;
                min_reads_prior = raw.reads;
            }
        }
    }

    let p_a = min_p * nraw;
    let p_p = min_p_prior;

    // Abundance-based bud.
    if p_a < b.omega_a
        && let Some((ci, r, raw_idx)) = mini
    {
        // Capture pre-pop state.
        let expected = b.raws[raw_idx].comp.lambda * b.clusters[ci].reads as f64;
        let birth_comp = b.raws[raw_idx].comp.clone();
        let birth_fold = b.raws[raw_idx].reads as f64 / expected.max(f64::MIN_POSITIVE);
        let nraw_total = b.raws.len() as u32;

        b.bi_pop_raw(ci, r);

        let mut new_bi = Bi::new(nraw_total);
        new_bi.birth_type = BirthType::Abundance;
        new_bi.birth_from = ci as u32;
        new_bi.birth_pval = p_a;
        new_bi.birth_fold = birth_fold;
        new_bi.birth_e = expected;
        new_bi.birth_comp = birth_comp;

        let new_ci = b.add_cluster(new_bi);
        b.bi_add_raw(new_ci, raw_idx);
        b.assign_center(new_ci);

        if verbose {
            eprint!(", Division (naive): Raw {raw_idx} from Bi {ci}, pA={p_a:.2e}");
        }
        return Some(new_ci);
    }

    // Prior-based bud.
    if p_p < b.omega_p
        && let Some((ci, r, raw_idx)) = mini_prior
    {
        let expected = b.raws[raw_idx].comp.lambda * b.clusters[ci].reads as f64;
        let birth_comp = b.raws[raw_idx].comp.clone();
        let birth_fold = b.raws[raw_idx].reads as f64 / expected.max(f64::MIN_POSITIVE);
        let nraw_total = b.raws.len() as u32;

        b.bi_pop_raw(ci, r);

        let mut new_bi = Bi::new(nraw_total);
        new_bi.birth_type = BirthType::Prior;
        new_bi.birth_pval = p_p;
        new_bi.birth_fold = birth_fold;
        new_bi.birth_e = expected;
        new_bi.birth_comp = birth_comp;

        let new_ci = b.add_cluster(new_bi);
        b.bi_add_raw(new_ci, raw_idx);
        b.assign_center(new_ci);

        if verbose {
            eprint!(", Division (prior): Raw {raw_idx} from Bi {ci}, pP={p_p:.2e}");
        }
        return Some(new_ci);
    }

    if verbose {
        let (raw_idx_str, reads, ci_str) = match mini {
            Some((ci, r, _)) => {
                let raw_idx = b.clusters[ci].raws[r];
                (raw_idx.to_string(), b.raws[raw_idx].reads, ci.to_string())
            }
            None => (
                init_center.to_string(),
                b.raws[init_center].reads,
                String::from("0"),
            ),
        };
        eprint!(
            ", No Division. Minimum pA={p_a:.2e} (Raw {raw_idx_str} w/ {reads} reads in Bi {ci_str})."
        );
    }
    None
}

/// Incremental equivalent of [`b_bud`] (issue #85).
///
/// Instead of rescanning every non-center raw across all clusters, this combines
/// the per-cluster cached minima (`Bi::bud_min` / `bud_min_prior`) that
/// `b_p_update` maintains, seeded — exactly as the scan is — with cluster 0's
/// center. The combine is O(nclusters); the cache refresh rides along on the
/// p-update pass at no extra memory traffic. The abundance/prior significance
/// tests and the pop-into-new-cluster tail are byte-identical to [`b_bud`].
///
/// Correctness rests on the invariant (verified in issue #85) that between
/// consecutive bud calls a candidate raw's p or eligibility changes *iff* its
/// cluster was flagged `update_e` — the exact set `b_p_update` reprices and
/// re-caches. A debug-only cross-check asserts the combined selection equals the
/// serial scan every call.
///
/// `min_fold`/`min_hamming`/`min_abund` must match the values passed to
/// `b_p_update` (they define which members were cached as candidates).
/// `combine_len` receives the O(nclusters) combine cost, mirroring `b_bud`'s
/// `raws_scanned` out-param for the redundancy diagnostic.
#[cfg_attr(not(debug_assertions), allow(unused_variables))]
pub fn b_bud_incremental(
    b: &mut B,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
    verbose: bool,
    combine_len: &mut u64,
) -> Option<usize> {
    let nraw = b.raws.len() as f64;
    let init_center = b.clusters[0]
        .center
        .expect("b_bud_incremental: cluster 0 has no center");

    *combine_len = b.clusters.len() as u64;

    // Seed with cluster 0's center, then combine per-cluster cached minima in
    // ascending cluster order with the same strict replacement as the scan
    // (keeps the lowest cluster index on a full (p, reads) tie; each cluster's
    // cache already resolved the lowest position within it).
    let mut mini: Option<(usize, usize, usize)> = None;
    let mut mini_prior: Option<(usize, usize, usize)> = None;
    let mut min_p = b.raws[init_center].p;
    let mut min_p_prior = b.raws[init_center].p;
    let mut min_reads = b.raws[init_center].reads;
    let mut min_reads_prior = b.raws[init_center].reads;

    for ci in 0..b.clusters.len() {
        if let Some(c) = b.clusters[ci].bud_min
            && (c.p < min_p || (c.p == min_p && c.reads > min_reads))
        {
            mini = Some((ci, c.r, c.raw_idx));
            min_p = c.p;
            min_reads = c.reads;
        }
        if let Some(c) = b.clusters[ci].bud_min_prior
            && (c.p < min_p_prior || (c.p == min_p_prior && c.reads > min_reads_prior))
        {
            mini_prior = Some((ci, c.r, c.raw_idx));
            min_p_prior = c.p;
            min_reads_prior = c.reads;
        }
    }

    // Debug guardrail: the incremental combine must match a full serial scan
    // exactly (same winners, same p/reads). Catches any cache-staleness bug the
    // update_e invariant would otherwise hide. Compiled out of release.
    #[cfg(debug_assertions)]
    {
        let (s_mini, s_min_p, s_min_reads, s_mini_prior, s_min_p_prior, s_min_reads_prior) =
            crate::pval::b_bud_scan_select(b, min_fold, min_hamming, min_abund);
        debug_assert_eq!(mini, s_mini, "b_bud_incremental: abundance mini mismatch");
        debug_assert_eq!(min_p.to_bits(), s_min_p.to_bits(), "min_p mismatch");
        debug_assert_eq!(min_reads, s_min_reads, "min_reads mismatch");
        debug_assert_eq!(mini_prior, s_mini_prior, "prior mini mismatch");
        debug_assert_eq!(
            min_p_prior.to_bits(),
            s_min_p_prior.to_bits(),
            "min_p_prior"
        );
        debug_assert_eq!(min_reads_prior, s_min_reads_prior, "min_reads_prior");
    }

    let p_a = min_p * nraw;
    let p_p = min_p_prior;

    // Abundance-based bud (tail identical to b_bud).
    if p_a < b.omega_a
        && let Some((ci, r, raw_idx)) = mini
    {
        let expected = b.raws[raw_idx].comp.lambda * b.clusters[ci].reads as f64;
        let birth_comp = b.raws[raw_idx].comp.clone();
        let birth_fold = b.raws[raw_idx].reads as f64 / expected.max(f64::MIN_POSITIVE);
        let nraw_total = b.raws.len() as u32;

        b.bi_pop_raw(ci, r);

        let mut new_bi = Bi::new(nraw_total);
        new_bi.birth_type = BirthType::Abundance;
        new_bi.birth_from = ci as u32;
        new_bi.birth_pval = p_a;
        new_bi.birth_fold = birth_fold;
        new_bi.birth_e = expected;
        new_bi.birth_comp = birth_comp;

        let new_ci = b.add_cluster(new_bi);
        b.bi_add_raw(new_ci, raw_idx);
        b.assign_center(new_ci);

        if verbose {
            eprint!(", Division (naive): Raw {raw_idx} from Bi {ci}, pA={p_a:.2e}");
        }
        return Some(new_ci);
    }

    // Prior-based bud (tail identical to b_bud).
    if p_p < b.omega_p
        && let Some((ci, r, raw_idx)) = mini_prior
    {
        let expected = b.raws[raw_idx].comp.lambda * b.clusters[ci].reads as f64;
        let birth_comp = b.raws[raw_idx].comp.clone();
        let birth_fold = b.raws[raw_idx].reads as f64 / expected.max(f64::MIN_POSITIVE);
        let nraw_total = b.raws.len() as u32;

        b.bi_pop_raw(ci, r);

        let mut new_bi = Bi::new(nraw_total);
        new_bi.birth_type = BirthType::Prior;
        new_bi.birth_pval = p_p;
        new_bi.birth_fold = birth_fold;
        new_bi.birth_e = expected;
        new_bi.birth_comp = birth_comp;

        let new_ci = b.add_cluster(new_bi);
        b.bi_add_raw(new_ci, raw_idx);
        b.assign_center(new_ci);

        if verbose {
            eprint!(", Division (prior): Raw {raw_idx} from Bi {ci}, pP={p_p:.2e}");
        }
        return Some(new_ci);
    }

    if verbose {
        let (raw_idx_str, reads, ci_str) = match mini {
            Some((ci, r, _)) => {
                let raw_idx = b.clusters[ci].raws[r];
                (raw_idx.to_string(), b.raws[raw_idx].reads, ci.to_string())
            }
            None => (
                init_center.to_string(),
                b.raws[init_center].reads,
                String::from("0"),
            ),
        };
        eprint!(
            ", No Division. Minimum pA={p_a:.2e} (Raw {raw_idx_str} w/ {reads} reads in Bi {ci_str})."
        );
    }
    None
}
