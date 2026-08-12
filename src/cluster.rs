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
    /// Portion of `comps_scanned` spent on the initial full build — the
    /// cluster-major rescan of every comp vec, paid fresh on every call (i.e.
    /// once per bud round). Any scheme that avoids rebuilding from scratch is
    /// bidding for this portion, so it bounds such a change's payoff (#124).
    pub comps_build: usize,
    /// Portion of `comps_scanned` spent recomputing affected raws from the
    /// candidate index after a move pass — raws whose best cluster may
    /// genuinely have changed, so this part does not go away by caching.
    pub comps_reconcile: usize,
    /// Candidates the build's bound skipped without evaluating (#124).
    ///
    /// `comps_build + comps_build_pruned` is the volume the old unbounded
    /// cluster-major scan would have evaluated, so the ratio is the prune rate
    /// realised on this workload. Note `comps_build` charges the hot clusters'
    /// comparisons twice when a raw's bounded pass re-touches them, which is
    /// real work — so the two do not sum to exactly the stored comp count.
    pub comps_build_pruned: usize,
    /// Wall time in the initial full build.
    ///
    /// Paired with `comps_build` this gives ns/comp for the sequential,
    /// cluster-major access pattern. Comparing that against the reconcile's
    /// ns/comp is what turns the comp-count split into a time split. Measured
    /// ~2x apart, so comp counts alone overstate the build's share: the
    /// contiguous scan is markedly cheaper per comparison than the scattered.
    pub build_time: std::time::Duration,
    /// Wall time in the affected-raw reconcile, summed over iterations.
    ///
    /// This walks the raw-major inverted index, so its reads are scattered and
    /// cost ~2x the build's per comparison. Any change that shifts work here
    /// from the build can be a net loss despite scanning fewer comparisons —
    /// price it in ns/comp, not comparisons (#87 was sized this way).
    pub reconcile_time: std::time::Duration,
    /// Cluster count at this call (context for the scan cost).
    pub nclusters: usize,
    /// Move-pass iterations that relocated no raws (pure convergence checks).
    ///
    /// NOTE: for `b_shuffle_converge` this is exactly the number of calls —
    /// every converge call ends in one zero-move iteration. It is not a
    /// measure of wasted work: a zero-move iteration breaks before the
    /// reconcile and scans nothing.
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
        // The serial scan is all build, by construction — it has no reconcile.
        comps_build: comps_scanned,
        comps_build_pruned: 0,
        comps_reconcile: 0,
        // Not instrumented on the serial path: it is the historical baseline,
        // superseded by the incremental driver.
        build_time: std::time::Duration::ZERO,
        reconcile_time: std::time::Duration::ZERO,
        nclusters: b.clusters.len(),
        zero_move_calls: usize::from(moves == 0),
        calls: 1,
    }
}

// ---------------------------------------------------------------------------
// b_shuffle_converge — incremental shuffle-to-convergence
// ---------------------------------------------------------------------------

/// Number of highest-abundance clusters the shuffle build scans cluster-major
/// before switching to the bounded raw-major pass (#124).
///
/// The modelled optimum is flat between 8 and 32 and degrades above it (at 128
/// pass 1 starts scanning clusters the bound would have pruned anyway), so 32
/// is chosen for margin rather than as a tuned value. Override with
/// `DADA2RS_SHUFFLE_HOT_CLUSTERS` to sweep it; it changes only the bound's
/// tightness, never the resulting partition.
const HOT_CLUSTERS: usize = 32;

fn hot_cluster_count() -> usize {
    std::env::var("DADA2RS_SHUFFLE_HOT_CLUSTERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(HOT_CLUSTERS)
}

/// One candidate cluster for a raw: (cluster index, lambda, hamming) — the
/// fields the shuffle scan needs from a stored `Comparison`.
#[derive(Clone, Copy)]
pub struct Cand {
    pub ci: u32,
    pub lambda: f64,
    pub hamming: u32,
}

/// Per-raw candidate index, maintained across the whole `run_dada`: each
/// cluster's stored comparisons are appended once, when the cluster is created
/// (see [`index_add_cluster`]). This lets a shuffle recompute a raw's best
/// cluster from just its own candidates, without rescanning every cluster — and
/// crucially it is built incrementally (O(new comps) per bud) rather than
/// rebuilt per shuffle loop, which is what makes the incremental driver a net
/// win rather than a wash.
/// Each raw's candidates are held **lambda-descending** (ties broken by
/// ascending `ci`). `lambda` is fixed by `b_compare` and never changes
/// afterwards — only the per-cluster scalar `reads` moves — so this order is
/// static and can be established once per candidate. It is what makes the
/// bounded build in [`b_shuffle_converge`] possible: a scan in this order can
/// stop as soon as `lambda × reads_bound` falls below the best `e` found so
/// far, because no later candidate can beat it.
///
/// Note this is *no longer* ascending-`ci` order, so scans over it cannot rely
/// on strict `>` for the lowest-`ci` tie-break and must compare `ci` explicitly
/// (see [`best_from_cands`]).
pub type CandIndex = Vec<Vec<Cand>>;

/// Insert cluster `ci`'s stored comparisons into the per-raw candidate index.
/// Must be called once per cluster, immediately after `b_compare`/
/// `b_compare_parallel` populates `clusters[ci].comp`, and in ascending `ci`
/// order (cluster 0 first, then each bud's new cluster).
///
/// Each candidate is placed so the raw's list stays lambda-descending. A full
/// re-sort (or a flat CSR rebuild) per bud would cost O(all comps) — exactly
/// the full-scan cost the bounded build exists to avoid — so instead each new
/// candidate is inserted at its position, O(list length) per candidate and only
/// for the raws this cluster actually compared against.
///
/// Ties: `ci` is the largest cluster index so far, and the insertion point is
/// the first entry with a *strictly smaller* lambda, so an equal-lambda
/// candidate lands after the ones already present — preserving ascending `ci`
/// within a run of equal lambdas.
pub fn index_add_cluster(index: &mut CandIndex, b: &B, ci: usize) {
    for comp in &b.clusters[ci].comp {
        let cands = &mut index[comp.index as usize];
        let at = cands.partition_point(|c| c.lambda >= comp.lambda);
        cands.insert(
            at,
            Cand {
                ci: ci as u32,
                lambda: comp.lambda,
                hamming: comp.hamming,
            },
        );
    }
}

/// Raw's best cluster over its candidate list at the clusters' current reads.
///
/// The list is lambda-descending, not `ci`-ascending, so the lowest-`ci`
/// tie-break is applied explicitly rather than falling out of strict `>`. The
/// resulting argmax is order-independent and identical to the serial scan's.
fn best_from_cands(cands: &[Cand], raw: usize, clusters: &[Bi]) -> Comparison {
    let mut best_e = f64::NEG_INFINITY;
    let mut best = Comparison::default();
    for c in cands {
        let e = c.lambda * clusters[c.ci as usize].reads as f64;
        if e > best_e || (e == best_e && c.ci < best.i) {
            best_e = e;
            best = Comparison {
                i: c.ci,
                index: raw as u32,
                lambda: c.lambda,
                hamming: c.hamming,
            };
        }
    }
    best
}

/// Build every raw's best cluster (`compmax`) and its expected read count
/// at the clusters' current `reads`, using the two-level bounded scan (#124).
/// Returns `(compmax, comps_examined, comps_pruned)`. `emax` is scratch — only
/// needed to resolve the max while scanning — so it is not returned; the
/// reconcile recomputes affected raws from the index instead.
///
/// `n_hot` is the number of highest-abundance clusters scanned cluster-major
/// before the bounded raw-major pass. It is a pure **performance** parameter:
/// the result is identical for every value of `n_hot`, which is what the
/// `bounded_build_matches_full_scan` test pins down. `n_hot >= nclusters`
/// degenerates to exactly the old unbounded cluster-major scan.
fn build_compmax(b: &B, index: &CandIndex, n_hot: usize) -> (Vec<Comparison>, usize, usize) {
    let nraw = b.raws.len();
    let mut compmax = vec![Comparison::default(); nraw];
    let mut emax = vec![f64::NEG_INFINITY; nraw];
    let mut comps_scanned = 0usize;

    let nclusters = b.clusters.len();
    let n_hot = n_hot.min(nclusters);
    // Sorted by reads descending; `ci` breaks ties so the choice of hot set is
    // deterministic run to run (it does not affect the result, only the bound).
    let mut hot: Vec<u32> = (0..nclusters as u32).collect();
    hot.sort_unstable_by_key(|&c| (std::cmp::Reverse(b.clusters[c as usize].reads), c));
    let reads_hot_min = if n_hot < nclusters {
        b.clusters[hot[n_hot] as usize].reads as f64
    } else {
        // Every cluster is hot: pass 2 has nothing left to find.
        0.0
    };
    let mut is_hot = vec![false; nclusters];
    for &c in &hot[..n_hot] {
        is_hot[c as usize] = true;
    }

    // Pass 1: hot clusters, cluster-major (contiguous).
    for &ci in &hot[..n_hot] {
        let bi = &b.clusters[ci as usize];
        let ci_reads = bi.reads as f64;
        for comp in &bi.comp {
            let idx = comp.index as usize;
            let e = comp.lambda * ci_reads;
            if e > emax[idx] || (e == emax[idx] && ci < compmax[idx].i) {
                emax[idx] = e;
                compmax[idx] = comp.clone();
            }
        }
        comps_scanned += bi.comp.len();
    }

    // Pass 2: cold candidates, raw-major, with the exact bound.
    //
    // Skipped entirely when every cluster was hot: pass 1 then already
    // evaluated every candidate exactly, so this pass could only re-touch them.
    // That is the common case for small pools (per-sample `dada`, where a B has
    // far fewer than HOT_CLUSTERS clusters), which therefore keep exactly the
    // old build's cost rather than paying a pass that can find nothing.
    let mut comps_pruned = 0usize;
    if n_hot < nclusters {
        for (raw, cands) in index.iter().enumerate() {
            let mut best_e = emax[raw];
            let mut examined = 0usize;
            for c in cands {
                // `cands` is lambda-descending, so once the bound fails here it
                // fails for every remaining candidate.
                //
                // Strictly `<`, not `<=`. A candidate whose bound merely *equals*
                // `best_e` can still tie it exactly, and a tie is won by the
                // lower `ci` — which, because the list is ordered by lambda and
                // not by `ci`, may be this later candidate. Exiting on `<=`
                // prunes it away and silently returns the wrong cluster (see
                // `bounded_build_keeps_lowest_ci_on_cross_lambda_ties`). The
                // pruning cost of this is nil: it only scans on through an
                // exact-equality run.
                if c.lambda * reads_hot_min < best_e {
                    break;
                }
                examined += 1;
                if is_hot[c.ci as usize] {
                    continue; // already evaluated exactly in pass 1
                }
                let e = c.lambda * b.clusters[c.ci as usize].reads as f64;
                if e > best_e || (e == best_e && c.ci < compmax[raw].i) {
                    best_e = e;
                    compmax[raw] = Comparison {
                        i: c.ci,
                        index: raw as u32,
                        lambda: c.lambda,
                        hamming: c.hamming,
                    };
                }
            }
            emax[raw] = best_e;
            comps_scanned += examined;
            comps_pruned += cands.len() - examined;
        }
    }
    (compmax, comps_scanned, comps_pruned)
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

    // Initial build: every raw's true best at the current reads.
    //
    // Two-level bounded scan (#124). The full cluster-major scan this replaces
    // evaluated every stored comparison — ~9.6e6 per build on a pooled 16S run
    // — to discover ~550 moves. Since `lambda` is fixed and only `reads`
    // changes, a raw's candidate list can be held lambda-descending once and
    // scanned with an exact early exit: if `lambda × reads_bound <= best_e`, no
    // later candidate can win.
    //
    // A single global `max_reads` bound is useless on real data, because
    // cluster abundances are power-law skewed and one huge cluster inflates the
    // bound for every raw (modelled in dev/shuffle-model: it regresses to
    // 0.83x). But the huge clusters are *few*, so:
    //
    //   pass 1 — the top HOT_CLUSTERS by reads, scanned cluster-major, the
    //            contiguous access pattern the old build used, seeding the map;
    //   pass 2 — everything else, raw-major over the candidate index, bounded
    //            by `reads_hot_min`, the smallest read count in pass 1's set,
    //            which under skew is orders of magnitude below max_reads.
    //
    // Exactness: both passes take the max with an explicit lowest-`ci`
    // tie-break, which is order-independent, so the result equals the serial
    // scan's `compmax` exactly — same max, same tie-break — and the move
    // decisions are unchanged.
    let t_build = std::time::Instant::now();
    let (mut compmax, mut comps_scanned, comps_pruned) =
        build_compmax(b, index, hot_cluster_count());
    let comps_build = comps_scanned;
    // Includes the compmax/emax allocation above: a scheme that reused the map
    // across buds would keep that buffer alive anyway, so it belongs here.
    let build_time = t_build.elapsed();
    let mut reconcile_time = std::time::Duration::ZERO;
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
        let t_rec = std::time::Instant::now();
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
        reconcile_time += t_rec.elapsed();
    }

    ShuffleStats {
        moves: total_moves,
        comps_scanned,
        comps_build,
        comps_build_pruned: comps_pruned,
        comps_reconcile: comps_scanned - comps_build,
        build_time,
        reconcile_time,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::containers::Raw;

    /// Synthetic partition with a skewed cluster-abundance distribution and a
    /// wide lambda spectrum — the regime the bounded build is designed for, and
    /// the one where a wrong bound would actually change the argmax.
    ///
    /// Returns the partition plus its candidate index (built exactly as
    /// `run_dada` builds it: `index_add_cluster` once per cluster, ascending).
    fn synthetic(nraw: usize, nclusters: usize) -> (B, CandIndex) {
        // Deterministic LCG — the test must not depend on rand.
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let raws: Vec<Raw> = (0..nraw)
            .map(|i| Raw::new(vec![1, 2, 3, 4], None, 1 + (i as u32 % 7), false))
            .collect();
        let mut b = B::new(raws, 1e-40, 1e-4, false);

        b.clusters.clear();
        for ci in 0..nclusters {
            let mut bi = Bi::new(nraw as u32);
            bi.i = ci as u32;
            // Power-law abundances: cluster 0 dwarfs the tail, which is what
            // makes a single global max_reads bound useless.
            // A band of clusters shares an identical read count so that, paired
            // with the duplicated lambdas below, raws hit exact `e` ties across
            // clusters. Without this the lowest-ci tie-break is never exercised
            // and a wrong tie-break would pass the test silently.
            bi.reads = if (10..20).contains(&ci) {
                50_000
            } else {
                1 + (1_000_000.0 / ((ci + 1) as f64)) as u32
            };
            for raw in 0..nraw {
                // Cluster 0 holds a comparison for every raw (b_compare runs
                // unscreened on it); the rest are sparse.
                if ci == 0 || next() % 8 == 0 {
                    // Lambda spanning many orders of magnitude, plus exact
                    // duplicates (the `% 5` bucket) so the lowest-ci tie-break
                    // is genuinely exercised rather than assumed unreachable.
                    // Within the equal-reads band, force a shared lambda often
                    // enough that ties in `e = lambda * reads` actually arise
                    // across those clusters.
                    let lambda = if (10..20).contains(&ci) && raw % 3 == 0 {
                        1e-3
                    } else if next() % 5 == 0 {
                        1e-3
                    } else {
                        ((next() % 1000) as f64 / -50.0).exp()
                    };
                    bi.comp.push(Comparison {
                        i: ci as u32,
                        index: raw as u32,
                        lambda,
                        hamming: (next() % 8) as u32,
                    });
                }
            }
            b.clusters.push(bi);
        }

        let mut index: CandIndex = vec![Vec::new(); nraw];
        for ci in 0..nclusters {
            index_add_cluster(&mut index, &b, ci);
        }
        (b, index)
    }

    /// The unbounded cluster-major scan the bounded build replaces: ascending
    /// `ci` with strict `>`, which keeps the lowest `ci` on ties. This is the
    /// reference semantics, taken straight from `b_shuffle2`.
    fn full_scan(b: &B) -> Vec<Comparison> {
        let mut compmax = vec![Comparison::default(); b.raws.len()];
        let mut emax = vec![f64::NEG_INFINITY; b.raws.len()];
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
        }
        compmax
    }

    /// The bounded build must be *exact*: identical to the full scan for every
    /// `n_hot`, including the lowest-`ci` tie-break. `n_hot` may only change how
    /// much work is done, never the answer — if it can change the answer, the
    /// ASV-concordance guardrail is live and the optimization is invalid.
    #[test]
    fn bounded_build_matches_full_scan() {
        let (b, index) = synthetic(2000, 60);
        let want = full_scan(&b);

        for n_hot in [0usize, 1, 2, 8, 32, 59, 60, 128] {
            let (got, examined, pruned) = build_compmax(&b, &index, n_hot);
            assert_eq!(got.len(), want.len());
            for (raw, (g, w)) in got.iter().zip(want.iter()).enumerate() {
                assert_eq!(
                    (g.i, g.index, g.lambda, g.hamming),
                    (w.i, w.index, w.lambda, w.hamming),
                    "n_hot={n_hot}: raw {raw} got cluster {}, want {}",
                    g.i,
                    w.i,
                );
            }
            // n_hot >= nclusters degenerates to the full scan, so it prunes
            // nothing; below that the bound must actually be doing work, or the
            // whole change is cost without benefit.
            if n_hot < b.clusters.len() {
                assert!(pruned > 0, "n_hot={n_hot}: bound pruned nothing");
            }
            assert!(examined > 0);
        }
    }

    /// The candidate index must stay lambda-descending as clusters are added,
    /// since the build's early exit is only sound on a sorted list. Ties must
    /// keep ascending `ci`, which is what preserves the lowest-`ci` tie-break.
    #[test]
    fn candidate_index_is_lambda_descending() {
        let (_b, index) = synthetic(500, 40);
        for (raw, cands) in index.iter().enumerate() {
            for w in cands.windows(2) {
                assert!(
                    w[0].lambda > w[1].lambda || (w[0].lambda == w[1].lambda && w[0].ci < w[1].ci),
                    "raw {raw}: candidates out of order ({}, ci {}) then ({}, ci {})",
                    w[0].lambda,
                    w[0].ci,
                    w[1].lambda,
                    w[1].ci,
                );
            }
        }
    }

    /// Cross-lambda ties: two clusters can reach the *same* `e` from different
    /// `lambda × reads` pairs. The candidate list is ordered by lambda, so the
    /// tied candidate with the lower `ci` — the one the full scan keeps — can
    /// sit *after* the higher-`ci` one. The bound must not prune it away, which
    /// is why the early exit tests `<` and not `<=`.
    ///
    /// Regression test for a real bug: with a `<=` exit this returned cluster 3
    /// where the full scan returns cluster 1.
    #[test]
    fn bounded_build_keeps_lowest_ci_on_cross_lambda_ties() {
        let raws = vec![Raw::new(vec![1, 2, 3, 4], None, 1, false)];
        let mut b = B::new(raws, 1e-40, 1e-4, false);
        b.clusters.clear();

        // e = lambda * reads:  c0 -> ~0, c2 -> ~0, and c1/c3/c4 all tie at 2.0.
        // The full scan keeps c1 (lowest ci). The tie is reached from three
        // different lambda/reads pairs, which pins both passes' tie-breaks:
        //   - c3 has the largest lambda, so it sorts FIRST in the
        //     lambda-descending index -> exercises pass 2's tie-break;
        //   - c4 has the largest reads, so pass 1 (which walks clusters in
        //     reads-descending order) visits it FIRST -> exercises pass 1's.
        for (ci, reads, lambda) in [
            (0u32, 1u32, 1e-9),
            (1, 200, 0.01),
            (2, 1, 1e-9),
            (3, 100, 0.02),
            (4, 400, 0.005),
        ] {
            let mut bi = Bi::new(1);
            bi.i = ci;
            bi.reads = reads;
            bi.comp.push(Comparison {
                i: ci,
                index: 0,
                lambda,
                hamming: 0,
            });
            b.clusters.push(bi);
        }
        assert_eq!(0.01 * 200.0, 0.02 * 100.0, "fixture must tie exactly");
        assert_eq!(0.01 * 200.0, 0.005 * 400.0, "fixture must tie exactly");

        let mut index: CandIndex = vec![Vec::new(); 1];
        for ci in 0..b.clusters.len() {
            index_add_cluster(&mut index, &b, ci);
        }
        assert_eq!(index[0][0].ci, 3, "c3 must sort first (larger lambda)");

        let want = full_scan(&b);
        assert_eq!(want[0].i, 1, "reference: lowest ci wins the tie");
        for n_hot in [0usize, 1, 2, 3, 4, 5] {
            let (got, _, _) = build_compmax(&b, &index, n_hot);
            assert_eq!(got[0].i, 1, "n_hot={n_hot}: tie must resolve to cluster 1");
        }
    }
}
