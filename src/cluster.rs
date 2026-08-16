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

/// Verify, every reconcile, that `compmax`/`emax` equal what a full rescan
/// would produce (#136).
///
/// Always on under `debug_assertions`; opt-in in release via
/// `DADA2RS_RECONCILE_VERIFY=1`. The release path exists because **the test
/// fixtures cannot reach the cases that matter**: exact cross-cluster ties in
/// `lambda × reads` never occur on them, and neither does an incumbent whose
/// fall changes the argmax — so mutations to the tie-break clause and to the
/// necessity marking both pass every fixture-based test, including this
/// invariant. Only production-scale data (thousands of clusters, millions of
/// candidate lists) exercises them, and a debug build is far too slow to run
/// there.
///
/// It is O(nraw x candidates) per reconcile — precisely the work being avoided
/// — so it is for validation runs, never timing runs.
fn reconcile_verify() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| {
        cfg!(debug_assertions) || std::env::var_os("DADA2RS_RECONCILE_VERIFY").is_some()
    })
}

/// Force the shuffle's reconcile to fully rescan every affected raw's candidate
/// list, disabling the #136 incremental update.
///
/// Same purpose as [`shuffle_no_prune`]: both arms of an A/B from one binary.
fn reconcile_full() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var_os("DADA2RS_RECONCILE_FULL").is_some())
}

/// Disable the #132 dirty-cluster pruning in the shuffle's move pass, forcing
/// the full scan.
///
/// Exists so both arms of an A/B run from **one binary**: the pruned and
/// unpruned paths must produce byte-identical partitions, and building two
/// binaries to check that has twice produced a run where the intended arm was
/// not actually the one measured. `DADA2RS_SHUFFLE_NO_PRUNE=1` selects the old
/// behaviour. Read once per process; undocumented in `--help` (a diagnostic,
/// not user-facing config).
fn shuffle_no_prune() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| std::env::var_os("DADA2RS_SHUFFLE_NO_PRUNE").is_some())
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

/// Per-comparison cost, in nanoseconds, split into the two halves that matter
/// for #127. `screen` is the k-mer prefilter, paid on every comparison;
/// `align` is the DP alignment and its post-processing, paid only by pairs the
/// screen lets through. `total - screen - align` is `compute_lambda` plus
/// per-item overhead.
#[derive(Clone, Copy, Default)]
struct CompCost {
    total: u64,
    screen: u64,
    dp: u64,
    post: u64,
}

/// Timing breakdown returned by [`b_compare_parallel`].
///
/// `busy / (map × nthreads)` is the map's parallel efficiency — if well below
/// 1, threads idle inside the parallel region (tail load-imbalance); near 1
/// with low OS-level utilization points to memory-bandwidth stalls instead.
/// Returned (not global-accumulated) so it stays correct under nested
/// per-sample parallelism.
#[derive(Clone, Copy, Default)]
pub struct CompareTiming {
    /// Parallel-pass wall time.
    pub map: std::time::Duration,
    /// Post-processing store-loop wall time (serial).
    pub serial: std::time::Duration,
    /// Summed per-item compute time across all worker threads.
    pub busy: std::time::Duration,
    /// Portion of `busy` in the k-mer screen (zero unless `measure`).
    ///
    /// The screen runs on 100% of comparisons while the aligner runs only on
    /// the 12-25% that pass it, so this ratio decides where `b_compare`'s cost
    /// actually lives — and therefore whether the lever is a better index or a
    /// better aligner (#127).
    pub screen: std::time::Duration,
    /// Portion of `busy` in the DP alignment kernel proper.
    pub dp: std::time::Duration,
    /// Portion of `busy` in `al2subs` + quality mapping, i.e. the
    /// post-alignment work paid only by pairs the screen passed.
    pub post: std::time::Duration,
    /// Comparisons that reached the k-mer screen (i.e. were not greedy-skipped).
    pub screened: u64,
    /// Comparisons that passed the screen and were aligned.
    pub aligned: u64,
}

/// Parallel version of `b_compare` using Rayon.
/// Equivalent to C++ `b_compare_parallel`.
///
/// Under `measure` the returned [`CompareTiming`] also carries the screen/align
/// split; the extra timer pair is verbose-only because it is not free relative
/// to a sub-microsecond k-mer screen.
pub fn b_compare_parallel(
    b: &mut B,
    i: usize,
    err_mat: &[f64],
    ncol: usize,
    params: &AlignParams,
    greedy: bool,
    measure: bool,
) -> CompareTiming {
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
    let comps: Vec<(f64, u32, bool, CompCost)> = (0..nraw)
        .into_par_iter()
        .with_max_len(par_max_len())
        .map_init(
            || {
                if measure {
                    AlignBuffers::measuring()
                } else {
                    AlignBuffers::new()
                }
            },
            |buf, index| {
                // Per-item timing only under `measure` (verbose) — keeps the hot
                // alignment loop allocation/Instant-free in production runs.
                let t0 = measure.then(std::time::Instant::now);
                let raw = &raws[index];
                let (lambda, hamming, skipped) = if greedy && (raw.reads > center_reads || raw.lock)
                {
                    let lambda = compute_lambda(raw, None, err_mat, ncol, use_quals);
                    (lambda, u32::MAX, true)
                } else {
                    let sub = sub_new_with_buf(&raws[center_idx], raw, params, buf);
                    let lambda = compute_lambda(raw, sub.as_ref(), err_mat, ncol, use_quals);
                    let hamming = sub.as_ref().map_or(u32::MAX, |s| s.nsubs() as u32);
                    (lambda, hamming, false)
                };
                // A greedy skip runs neither half, so leave both at zero rather
                // than reading the buffer's stale values from the previous item.
                let (screen, dp, post) = if skipped || !measure {
                    (0, 0, 0)
                } else {
                    (
                        buf.last_screen_nanos,
                        buf.last_dp_nanos,
                        buf.last_post_nanos,
                    )
                };
                let nanos = t0.map_or(0, |t| t.elapsed().as_nanos() as u64);
                (
                    lambda,
                    hamming,
                    skipped,
                    CompCost {
                        total: nanos,
                        screen,
                        dp,
                        post,
                    },
                )
            },
        )
        .collect();
    let map_dur = t_map.elapsed();
    let cost = CompCost {
        total: comps.iter().map(|c| c.3.total).sum(),
        screen: comps.iter().map(|c| c.3.screen).sum(),
        dp: comps.iter().map(|c| c.3.dp).sum(),
        post: comps.iter().map(|c| c.3.post).sum(),
    };
    let busy_dur = std::time::Duration::from_nanos(cost.total);
    // Denominators for the ns/comp figures: the screen runs on every
    // non-greedy-skipped comparison, the aligner only on those it passes.
    let screened = comps.iter().filter(|c| !c.2).count() as u64;
    let aligned = comps.iter().filter(|c| !c.2 && c.1 != u32::MAX).count() as u64;

    // Serial post-processing: selectively store comparisons.
    let t_serial = std::time::Instant::now();
    let total_reads = b.reads as f64;
    for (index, (lambda, hamming, skipped, _cost)) in comps.into_iter().enumerate() {
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
    CompareTiming {
        map: map_dur,
        serial: t_serial.elapsed(),
        busy: busy_dur,
        screen: std::time::Duration::from_nanos(cost.screen),
        dp: std::time::Duration::from_nanos(cost.dp),
        post: std::time::Duration::from_nanos(cost.post),
        screened,
        aligned,
    }
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
    /// Raw count at this call. With `nclusters` this normalises every other
    /// counter here, which is what makes the numbers comparable across
    /// platforms (a PacBio pool has far fewer, far longer raws than a MiSeq
    /// one, so absolute comparison counts alone say little).
    pub nraw: usize,
    /// Wall time in the move passes.
    ///
    /// The third phase of a shuffle, and the one that went unmeasured until
    /// #124's build arm closed: `shuffle - build - reconcile` was a consistent
    /// 15-19% of shuffle time with nothing accounting for it. A move pass walks
    /// every cluster's membership list — O(nraw) per iteration, with a
    /// scattered `compmax[raw]` read each — to discover a few hundred moves.
    pub move_time: std::time::Duration,
    /// Raws visited by the move passes, summed over iterations. Against
    /// `moves`, this is the move pass's redundancy ratio — the same
    /// scanned-per-outcome figure that sized the build.
    pub move_raws_scanned: usize,
    /// Raws collected into the affected set and recomputed by the reconcile,
    /// summed over iterations. `comps_reconcile / reconcile_affected` is the
    /// mean candidate-list length actually walked.
    pub reconcile_affected: usize,
    /// Raws a *full* (unpruned) move pass would have visited, against
    /// `move_raws_scanned` for what the dirty-cluster pass actually visited
    /// (issue #132). Their ratio is the realised prune.
    ///
    /// Kept because the prune is workload-dependent — it was 64-67% on MiSeq
    /// and 57% on PacBio — so a future change that erodes it should be visible
    /// rather than inferred from wall time.
    pub move_raws_unpruned: usize,
    /// Dirty clusters summed over the pruned passes, for the mean per pass.
    pub move_dirty_clusters: usize,
    /// Move passes that followed a reconcile and were therefore pruned. The
    /// remainder follow a build and must scan every cluster.
    pub move_passes_prunable: usize,
    /// Total move passes, pruned or not.
    pub move_passes: usize,
    /// Affected raws whose full candidate rescan is *provably necessary*,
    /// because their current best cluster's reads **decreased** (issue #136).
    ///
    /// The reconcile currently rescans every affected raw's whole candidate
    /// list. Most of that is redundant — only 0.006-0.034% of recomputes change
    /// anything — and #124 called the redundancy unreachable, on the grounds
    /// that finding which raws changed costs the scattered access being avoided.
    ///
    /// That is true of *avoiding the touch*, but not of doing it more cheaply.
    /// A raw's rescan is only needed when its current best shrank: if that
    /// cluster's reads held or rose, its score still beats every *unchanged*
    /// candidate (theirs are unchanged and were already below the old max), so
    /// only the changed candidates need testing — and those are reachable by
    /// walking the changed clusters' own comp vectors, which is sequential and
    /// already happening to collect the affected set.
    ///
    /// So this counter, against `reconcile_affected`, is the fraction of the
    /// reconcile that a cheaper scheme could not avoid. Measurement only.
    pub reconcile_rescan_raws: usize,
    /// Comparisons those provably-necessary rescans would walk, against
    /// `comps_reconcile` for what the current reconcile walks.
    pub reconcile_comps_rescan: usize,
    /// Reconcile time in the *collect* half: walking the changed clusters' own
    /// comp vectors to gather the affected raws. Sequential, cluster-major, and
    /// unavoidable — a cheaper scheme still has to do this.
    pub reconcile_collect_time: std::time::Duration,
    /// Reconcile time in the *rescan* half: recomputing each affected raw's
    /// best over its full candidate list. Scattered, raw-major, and the part
    /// `reconcile_rescan_raws` says is mostly unnecessary.
    pub reconcile_rescan_time: std::time::Duration,
    /// Of those recomputes, how many changed the raw's best cluster.
    ///
    /// This is the sizing number for any reconcile optimization: the rest of
    /// the recomputes returned the value already in `compmax`, so a scheme that
    /// could predict them is bidding for that fraction — exactly the way
    /// `comps_build` bounded #87's payoff. A low ratio means the reconcile is
    /// mostly re-deriving what it already knew.
    ///
    /// NOTE the incremental path (#136) counts **promotions**, not net changes:
    /// a raw can be promoted more than once within one reconcile as successive
    /// changed clusters are visited, and could in principle return to where it
    /// started. Measured drift against the full-rescan path is ~0.02% (55,289
    /// vs 55,279 on MiSeq R1), far below anything this figure is used to decide.
    /// The dirty-cluster marking is unaffected — it keys on the raw's *holder*,
    /// which does not move during a reconcile, so repeated marking is
    /// idempotent.
    pub reconcile_changed: usize,
    /// Work done by the **first** reconcile of this call — the one immediately
    /// following the post-bud move pass ([#139](https://github.com/HPCBio/dada2-rs/issues/139),
    /// reviving #87's projection).
    ///
    /// This is the direct measurement of what carrying `compmax` across buds
    /// would *relocate*. #87 removes the per-bud build but forces the
    /// newly-budded state to be reconciled instead, and this first reconcile is
    /// that reconcile. Averaged reconcile volume understates it, because a bud
    /// both adds a cluster and steals raws from many others, perturbing more
    /// comp volume than a mid-loop iteration does.
    ///
    /// **The model has changed since #87 was closed.** It priced the relocated
    /// work at the *full-rescan* reconcile's scattered rate, and #136 replaced
    /// that with an incremental update: most affected raws are now settled by a
    /// sequential cluster-major test, and only those whose own best cluster
    /// fell pay a candidate rescan. So the relocated cost now has two parts,
    /// counted separately here, and #87's break-even must be re-derived rather
    /// than re-read.
    ///
    /// Pairs walked by the incremental test (sequential, cluster-major).
    pub pairs_first_reconcile: usize,
    /// (raw, changed-cluster) pairs walked by the incremental test across *all*
    /// reconciles — the denominator that turns `reconcile_collect_time` into a
    /// per-pair rate, so the projection can be priced at this run's own cost.
    pub pairs_reconcile: usize,
    /// Candidate comparisons rescanned by the first reconcile (scattered,
    /// raw-major) — the part that survives #136.
    pub comps_first_reconcile: usize,
    /// Calls that reached at least one reconcile, the denominator for both.
    /// Below `calls`, because a converge call that moves nothing breaks before
    /// reconciling.
    pub first_reconcile_calls: usize,
    /// Times the incremental reconcile's replacement decided on an **exact
    /// tie** in `lambda × reads`, resolved to the lower `ci` (#136).
    ///
    /// Worth counting because a passing correctness run cannot distinguish
    /// "the tie rule is right" from "ties never happen". They never happen on
    /// any fixture in this repo — which is why dropping the tie clause passes
    /// every fixture-based test — so this says whether the case is live on real
    /// data, and therefore whether #124's latent-bug class is reachable here.
    pub reconcile_tie_breaks: usize,
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
        comps_reconcile: 0,
        // Not instrumented on the serial path: it is the historical baseline,
        // superseded by the incremental driver.
        build_time: std::time::Duration::ZERO,
        reconcile_time: std::time::Duration::ZERO,
        nclusters: b.clusters.len(),
        nraw,
        // Not instrumented on the serial path (superseded by the incremental
        // driver); it has no reconcile and its move pass is not split out.
        move_time: std::time::Duration::ZERO,
        move_raws_scanned: 0,
        move_raws_unpruned: 0,
        move_dirty_clusters: 0,
        move_passes_prunable: 0,
        move_passes: 0,
        reconcile_affected: 0,
        reconcile_rescan_raws: 0,
        reconcile_tie_breaks: 0,
        pairs_first_reconcile: 0,
        pairs_reconcile: 0,
        comps_first_reconcile: 0,
        first_reconcile_calls: 0,
        reconcile_comps_rescan: 0,
        reconcile_collect_time: std::time::Duration::ZERO,
        reconcile_rescan_time: std::time::Duration::ZERO,
        reconcile_changed: 0,
        zero_move_calls: usize::from(moves == 0),
        calls: 1,
    }
}

// ---------------------------------------------------------------------------
// b_shuffle_converge — incremental shuffle-to-convergence
// ---------------------------------------------------------------------------

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
pub type CandIndex = Vec<Vec<Cand>>;

/// Append cluster `ci`'s stored comparisons to the per-raw candidate index.
/// Must be called once per cluster, immediately after `b_compare`/
/// `b_compare_parallel` populates `clusters[ci].comp`, and in ascending `ci`
/// order (cluster 0 first, then each bud's new cluster). Ascending order means
/// each raw's candidate list stays cluster-ascending, so a strict-`>` scan keeps
/// the lowest cluster index on ties — matching the serial scan's tie-break.
pub fn index_add_cluster(index: &mut CandIndex, b: &B, ci: usize) {
    for comp in &b.clusters[ci].comp {
        index[comp.index as usize].push(Cand {
            ci: ci as u32,
            lambda: comp.lambda,
            hamming: comp.hamming,
        });
    }
}

/// One cluster's move pass: relocate every member whose `compmax` names a
/// different cluster. Shared by the full and dirty-cluster scans so the two
/// cannot drift apart (issue #132).
///
/// Iterates backwards because `bi_pop_raw` uses `swap_remove`: the vacated slot
/// is filled from the end, which a descending walk has already passed.
#[inline]
fn move_pass_cluster(b: &mut B, ci: usize, compmax: &[Comparison], scanned: &mut usize) -> usize {
    *scanned += b.clusters[ci].raws.len();
    let mut moves = 0usize;
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
    moves
}

/// Does candidate `(e, ci)` displace incumbent `(e_cur, ci_cur)`?
///
/// `best_from_cands_scored` is argmax by `(score, -ci)`: it scans ascending `ci` with
/// strict `>`, so an exact tie resolves to the **lowest** `ci`. Reproducing that
/// pairwise needs the tie clause — strict `>` alone is wrong whenever the
/// incumbent has the higher `ci`, which the incremental reconcile (#136) can
/// easily produce since it never scans in `ci` order from scratch.
///
/// Exact float equality is deliberate. Cross-cluster exact ties in
/// `lambda × reads` are the case where #124's pruning arm hid a latent bug that
/// survived both benchmarking and ASV concordance.
#[inline]
fn beats(e: f64, e_cur: f64, ci: u32, ci_cur: u32) -> bool {
    e > e_cur || (e == e_cur && ci < ci_cur)
}

/// Note that `raw`'s best cluster changed: count it, and mark the cluster that
/// currently *holds* the raw dirty so the next move pass visits it (#132).
///
/// Its current cluster, not the new best — the move pass finds movers by
/// walking membership, so it has to look where the raw is now.
#[inline]
#[allow(clippy::too_many_arguments)]
fn record_change(
    raw: usize,
    new_ci: u32,
    compmax: &[Comparison],
    b: &B,
    dirty_cluster: &mut [bool],
    dirty_list: &mut Vec<u32>,
    reconcile_changed: &mut usize,
) {
    if new_ci == compmax[raw].i {
        return;
    }
    *reconcile_changed += 1;
    let cur = b.raw_cluster[raw];
    if cur != u32::MAX && !dirty_cluster[cur as usize] {
        dirty_cluster[cur as usize] = true;
        dirty_list.push(cur);
    }
}

/// Raw's best cluster over its candidate list at the clusters' current reads,
/// with the winning score.
///
/// Ascending-`ci` order + strict `>` reproduces the serial lowest-`ci`
/// tie-break. The score is returned because the incremental reconcile (#136)
/// keeps `emax` in step with `compmax`, and recomputing it from the returned
/// `Comparison` would re-read the cluster's reads for nothing.
fn best_from_cands_scored(cands: &[Cand], raw: usize, clusters: &[Bi]) -> (Comparison, f64) {
    let mut best_e = f64::NEG_INFINITY;
    let mut best = Comparison::default();
    for c in cands {
        let e = c.lambda * clusters[c.ci as usize].reads as f64;
        if e > best_e {
            best_e = e;
            best = Comparison {
                i: c.ci,
                index: raw as u32,
                lambda: c.lambda,
                hamming: c.hamming,
            };
        }
    }
    (best, best_e)
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

    // #132: `b.raw_cluster` must agree with actual membership, or the move
    // pass prunes against the wrong clusters and silently drops moves. The
    // map is maintained across the whole run by `bi_add_raw`, far from here,
    // so the invariant is asserted rather than assumed.
    //
    // Debug-only: it is an O(nraw) walk, which is the cost the pruning exists
    // to avoid. An end-to-end equivalence test does *not* reliably catch a
    // stale map — on a small fixture most raws sit in cluster 0, which gets
    // marked dirty anyway, so the scan happens to stay complete.
    #[cfg(debug_assertions)]
    {
        for (ci, c) in b.clusters.iter().enumerate() {
            for &raw in &c.raws {
                debug_assert_eq!(
                    b.raw_cluster[raw], ci as u32,
                    "raw_cluster[{raw}] says {} but raw is in cluster {ci}",
                    b.raw_cluster[raw]
                );
            }
        }
    }

    // Initial build: every raw's true best at the current reads. Done the
    // serial way — a contiguous, cache-friendly scan of the per-cluster comp
    // vecs (cluster-major), NOT the raw-major inverted index. This is the bulk
    // of the work each loop, and sequential access here is far cheaper per
    // comparison than the index's scattered reads (which is what made a
    // fully-index-based build a net wall loss despite fewer comparisons). The
    // index is reserved for the reconcile, where the touched-raw volume is
    // small. Byte-identical to b_shuffle2's build: ascending ci + strict `>`
    // keeps the lowest-ci max.
    let t_build = std::time::Instant::now();
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
    let comps_build = comps_scanned;
    // Includes the compmax/emax allocation above: a scheme that reused the map
    // across buds would keep that buffer alive anyway, so it belongs here.
    let build_time = t_build.elapsed();
    let mut reconcile_time = std::time::Duration::ZERO;
    // `emax` is carried forward (#136). It used to be dropped here, because the
    // reconcile recomputed every affected raw from its candidate list and so
    // never needed the incumbent's score. The incremental reconcile does need
    // it: testing a changed candidate against the current best is exactly a
    // comparison of scores. +8 bytes/raw, against a pooled peak measured in GB.

    // Reads the map is currently consistent with (for dirty detection).
    let mut reads_used: Vec<u32> = b.clusters.iter().map(|c| c.reads).collect();

    // Reused scratch for the affected-raw set (avoid per-iteration nraw alloc).
    let mut in_affected = vec![false; nraw];
    let mut affected: Vec<usize> = Vec::new();

    // --- #132: dirty-cluster move pass ---
    // After a reconcile, only raws whose `compmax` changed can move, so the
    // move pass need visit only the clusters *holding* those raws. Tracking
    // dirty clusters (rather than dirty raws) keeps the walk cluster-major and
    // sequential — the access pattern is the whole point, since #124 showed a
    // scattered pass at 12.4-14.1 ns/comp loses to this 2.6-3.1 ns/raw one
    // unless it prunes below ~35%.
    //
    // `b.raw_cluster` supplies the raw -> cluster mapping and is maintained
    // globally by `bi_add_raw`.
    let mut dirty_cluster = vec![false; b.clusters.len()];
    let mut dirty_list: Vec<u32> = Vec::new();
    // False only on the first pass, which follows the build and must scan
    // everything (`compmax` was just rebuilt wholesale). Every later pass is
    // preceded by a reconcile, which sets this.
    let mut after_reconcile = false;
    let mut move_raws_unpruned = 0usize;
    let mut move_dirty_clusters = 0usize;
    let mut move_passes_prunable = 0usize;
    let mut move_passes = 0usize;

    let mut total_moves = 0usize;
    let mut nshuffle = 0usize;
    let mut zero_move_calls = 0usize;
    let mut move_time = std::time::Duration::ZERO;
    let mut move_raws_scanned = 0usize;
    let mut reconcile_affected = 0usize;
    // Clusters whose reads *decreased* in the current reconcile. A raw whose
    // best sits in one of these is the only case a cheaper reconcile could not
    // settle without a full candidate rescan — see `reconcile_rescan_raws`.
    let mut reads_fell = vec![false; b.clusters.len()];
    let mut reads_fell_list: Vec<u32> = Vec::new();
    let mut reconcile_rescan_raws = 0usize;
    let mut reconcile_tie_breaks = 0usize;
    let mut pairs_first_reconcile = 0usize;
    let mut comps_first_reconcile = 0usize;
    let mut first_reconcile_calls = 0usize;
    let mut pairs_reconcile = 0usize;
    let mut reconcile_comps_rescan = 0usize;
    let mut reconcile_collect_time = std::time::Duration::ZERO;
    let mut reconcile_rescan_time = std::time::Duration::ZERO;
    let mut reconcile_changed = 0usize;
    loop {
        // Move pass — identical to b_shuffle2's, using the current compmax.
        // Timed separately (#124): build + reconcile left 15-19% of shuffle
        // time unaccounted for, and this is where it goes.
        let t_move = std::time::Instant::now();
        let mut moves = 0usize;
        move_passes += 1;
        move_raws_unpruned += b.clusters.iter().map(|c| c.raws.len()).sum::<usize>();
        if after_reconcile && !shuffle_no_prune() {
            move_passes_prunable += 1;
            move_dirty_clusters += dirty_list.len();
            // Ascending cluster order, matching the full scan. The move outcome
            // is order-independent (`compmax` is fixed for the pass and each raw
            // goes to `compmax[raw].i` regardless of when it is visited), but
            // keeping the order identical means the *sequence* of pops and adds
            // is too — so `swap_remove` shuffles membership vectors the same
            // way, and anything downstream that reads them positionally sees
            // what it saw before.
            dirty_list.sort_unstable();
            #[allow(clippy::needless_range_loop)]
            // Indexed rather than iterated: `move_pass_cluster` takes `&mut b`,
            // and `dirty_list` is a local so the borrows do not overlap, but an
            // iterator over it would hold a borrow across the call.
            for k in 0..dirty_list.len() {
                let ci = dirty_list[k] as usize;
                moves += move_pass_cluster(b, ci, &compmax, &mut move_raws_scanned);
            }
        } else {
            for ci in 0..b.clusters.len() {
                moves += move_pass_cluster(b, ci, &compmax, &mut move_raws_scanned);
            }
        }
        for &ci in &dirty_list {
            dirty_cluster[ci as usize] = false;
        }
        dirty_list.clear();
        dirty_cluster.resize(b.clusters.len(), false);
        move_time += t_move.elapsed();

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
        let t_collect = std::time::Instant::now();
        // #139: `nshuffle == 1` means this is the first reconcile of the call —
        // the one settling the state the bud just created, and therefore the
        // work #87 would relocate.
        let first_reconcile = nshuffle == 1;
        let comps_before_rec = comps_scanned;
        let mut pairs_this_reconcile = 0usize;
        reads_fell.resize(b.clusters.len(), false);
        let incremental = !reconcile_full();

        // --- Pass A1: necessity marking (#136) ---
        // A raw needs a full candidate rescan only if its current best cluster's
        // reads FELL. If they held or rose, that cluster still beats every
        // *unchanged* candidate — theirs are unchanged and were already below
        // the old max — so only the changed candidates need testing, and those
        // are reachable cluster-major in A2.
        //
        // This must run before A2 and separately from it: it reads the
        // *entry-time* incumbent, and A2 mutates `compmax` as it goes.
        for (ci, ru) in reads_used.iter().enumerate() {
            if b.clusters[ci].reads < *ru && !reads_fell[ci] {
                reads_fell[ci] = true;
                reads_fell_list.push(ci as u32);
            }
        }
        for &ci in &reads_fell_list {
            for comp in &b.clusters[ci as usize].comp {
                let raw = comp.index as usize;
                if compmax[raw].i == ci && !in_affected[raw] {
                    in_affected[raw] = true;
                    affected.push(raw);
                }
            }
        }
        let n_rescan = affected.len();
        reconcile_rescan_raws += n_rescan;
        for &raw in &affected {
            reconcile_comps_rescan += index[raw].len();
        }

        // --- Pass A2: incremental test (#136) ---
        // Every changed cluster, ascending `ci`, skipping rescan-marked raws.
        //
        // The replacement rule reproduces `best_from_cands` exactly, which is
        // argmax by (score, -ci): strict `>` alone is NOT sufficient, because an
        // incumbent with a *higher* ci must lose an exact tie to a lower-ci
        // candidate. Exact float equality is deliberate — cross-cluster exact
        // ties are the case where #124's pruning arm hid a latent bug.
        // Indexed rather than zipped: the loop body needs `&mut` access to
        // `compmax`/`emax`/`dirty_*` alongside `&b`, so holding iterators over
        // `b.clusters` and `reads_used` across it would over-borrow.
        #[allow(clippy::needless_range_loop)]
        for ci in 0..b.clusters.len() {
            let new_reads = b.clusters[ci].reads;
            if new_reads == reads_used[ci] {
                continue;
            }
            let reads_f = new_reads as f64;
            if incremental {
                pairs_this_reconcile += b.clusters[ci].comp.len();
                for comp in &b.clusters[ci].comp {
                    let raw = comp.index as usize;
                    if in_affected[raw] {
                        continue; // queued for full rescan
                    }
                    let e = comp.lambda * reads_f;
                    if compmax[raw].i == ci as u32 {
                        // Same cluster, re-priced at the new reads. Its own
                        // score is stale, so this is an assignment, not a test.
                        emax[raw] = e;
                    } else if beats(e, emax[raw], ci as u32, compmax[raw].i) {
                        if e == emax[raw] {
                            reconcile_tie_breaks += 1;
                        }
                        emax[raw] = e;
                        record_change(
                            raw,
                            ci as u32,
                            &compmax,
                            b,
                            &mut dirty_cluster,
                            &mut dirty_list,
                            &mut reconcile_changed,
                        );
                        compmax[raw] = comp.clone();
                    }
                }
            } else {
                // Baseline path (DADA2RS_RECONCILE_FULL): collect every raw in a
                // changed cluster and rescan them all, as before #136.
                for comp in &b.clusters[ci].comp {
                    let raw = comp.index as usize;
                    if !in_affected[raw] {
                        in_affected[raw] = true;
                        affected.push(raw);
                    }
                }
            }
            reads_used[ci] = new_reads;
        }
        reconcile_collect_time += t_collect.elapsed();
        reconcile_affected += affected.len();

        // --- Pass B: full rescan, for the marked raws only ---
        let t_rescan = std::time::Instant::now();
        for &raw in &affected {
            let (new_best, new_e) = best_from_cands_scored(&index[raw], raw, &b.clusters);
            if new_best.i != compmax[raw].i {
                record_change(
                    raw,
                    new_best.i,
                    &compmax,
                    b,
                    &mut dirty_cluster,
                    &mut dirty_list,
                    &mut reconcile_changed,
                );
            }
            compmax[raw] = new_best;
            emax[raw] = new_e;
            comps_scanned += index[raw].len();
            in_affected[raw] = false;
        }
        reconcile_rescan_time += t_rescan.elapsed();
        pairs_reconcile += pairs_this_reconcile;
        if first_reconcile {
            first_reconcile_calls = 1;
            pairs_first_reconcile = pairs_this_reconcile;
            comps_first_reconcile = comps_scanned - comps_before_rec;
        }

        // #136: `compmax`/`emax` must equal what a full rescan would produce.
        // The incremental path reaches that answer by a different route —
        // testing only changed candidates against a carried incumbent — so the
        // equality is asserted, not assumed. See `reconcile_verify` for why
        // this must be runnable in release: the fixtures cannot reach the cases
        // that matter.
        if reconcile_verify() {
            for raw in 0..nraw {
                let (want, want_e) = best_from_cands_scored(&index[raw], raw, &b.clusters);
                assert_eq!(
                    compmax[raw].i, want.i,
                    "reconcile: compmax[{raw}] = cluster {} but a full rescan says {}",
                    compmax[raw].i, want.i
                );
                assert_eq!(
                    emax[raw], want_e,
                    "reconcile: emax[{raw}] = {} but a full rescan says {want_e}",
                    emax[raw]
                );
            }
        }

        affected.clear();
        for &ci in &reads_fell_list {
            reads_fell[ci as usize] = false;
        }
        reads_fell_list.clear();
        // The next move pass follows a reconcile, so its dirty set is known —
        // *even when nothing changed*, in which case the set is empty and the
        // pass provably moves nothing. Setting this only on a change would
        // credit a no-op pass with a full scan and understate the prune.
        after_reconcile = true;
        reconcile_time += t_rec.elapsed();
    }

    ShuffleStats {
        moves: total_moves,
        comps_scanned,
        comps_build,
        comps_reconcile: comps_scanned - comps_build,
        build_time,
        reconcile_time,
        nclusters: b.clusters.len(),
        nraw,
        move_time,
        move_raws_scanned,
        move_raws_unpruned,
        move_dirty_clusters,
        move_passes_prunable,
        move_passes,
        reconcile_affected,
        reconcile_rescan_raws,
        reconcile_tie_breaks,
        pairs_first_reconcile,
        pairs_reconcile,
        comps_first_reconcile,
        first_reconcile_calls,
        reconcile_comps_rescan,
        reconcile_collect_time,
        reconcile_rescan_time,
        reconcile_changed,
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
mod reconcile_rule_tests {
    use super::*;

    /// Reference semantics: exactly what `best_from_cands_scored` does — ascending
    /// `ci`, strict `>`.
    fn reference_argmax(cands: &[(f64, u32)]) -> u32 {
        let mut sorted = cands.to_vec();
        sorted.sort_by_key(|&(_, ci)| ci);
        let mut best_e = f64::NEG_INFINITY;
        let mut best_ci = u32::MAX;
        for &(e, ci) in &sorted {
            if e > best_e {
                best_e = e;
                best_ci = ci;
            }
        }
        best_ci
    }

    /// Fold `beats` over the candidates in an arbitrary order and require the
    /// same winner as the reference.
    ///
    /// The arbitrary order is the point: the incremental reconcile visits only
    /// *changed* clusters, so it meets candidates in an order unrelated to `ci`,
    /// and an incumbent with a higher `ci` than a later candidate is routine.
    fn fold_beats(cands: &[(f64, u32)]) -> u32 {
        let mut e_cur = f64::NEG_INFINITY;
        let mut ci_cur = u32::MAX;
        for &(e, ci) in cands {
            if beats(e, e_cur, ci, ci_cur) {
                e_cur = e;
                ci_cur = ci;
            }
        }
        ci_cur
    }

    /// **Exact ties, in both orders.** This is the case no fixture in this repo
    /// reaches — `lambda × reads` collisions across clusters simply do not occur
    /// on them — so dropping the tie clause passes every end-to-end test,
    /// including the full-rescan invariant. It is caught only here.
    #[test]
    fn exact_ties_resolve_to_lowest_ci() {
        // Incumbent has the HIGHER ci and is met first: the tie must flip it.
        assert_eq!(fold_beats(&[(10.0, 7), (10.0, 3)]), 3);
        assert_eq!(reference_argmax(&[(10.0, 7), (10.0, 3)]), 3);
        // Incumbent has the LOWER ci: the tie must not flip it.
        assert_eq!(fold_beats(&[(10.0, 3), (10.0, 7)]), 3);
        assert_eq!(reference_argmax(&[(10.0, 3), (10.0, 7)]), 3);
        // Three-way tie, worst order.
        assert_eq!(fold_beats(&[(5.0, 9), (5.0, 4), (5.0, 6)]), 4);
        assert_eq!(reference_argmax(&[(5.0, 9), (5.0, 4), (5.0, 6)]), 4);
    }

    /// A tie produced the way production produces one: different `lambda` and
    /// `reads` whose product is bit-identical.
    #[test]
    fn ties_from_distinct_lambda_reads_products() {
        let a = (0.25_f64 * 400.0, 11u32); // 100.0
        let b = (0.5_f64 * 200.0, 4u32); // 100.0
        assert_eq!(
            a.0, b.0,
            "test precondition: products must be bit-identical"
        );
        assert_eq!(fold_beats(&[a, b]), 4);
        assert_eq!(fold_beats(&[b, a]), 4);
    }

    /// Order-independence on non-tied inputs, including negatives and zero.
    #[test]
    fn fold_matches_reference_in_any_order() {
        let cands = [(3.0, 5), (7.5, 2), (7.5, 9), (0.0, 1), (-2.0, 0), (7.5, 4)];
        let want = reference_argmax(&cands);
        assert_eq!(
            want, 2,
            "reference: highest score 7.5, lowest ci among ties"
        );
        let mut perm: Vec<_> = cands.to_vec();
        for rot in 0..cands.len() {
            perm.rotate_left(rot.min(1));
            assert_eq!(fold_beats(&perm), want, "rotation {rot} disagreed");
        }
        assert_eq!(fold_beats(&[(-2.0, 0)]), 0, "single candidate");
    }
}
