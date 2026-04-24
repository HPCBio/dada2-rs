//! Cluster operations: compare, shuffle, bud.
//!
//! Ports `cluster.cpp`, excluding R/Rcpp and RcppParallel wrappers.
//! Parallel comparisons use Rayon in place of RcppParallel.

use rayon::prelude::*;

use crate::containers::{B, Bi, BirthType, Comparison};
use crate::nwalign::{sub_new_with_buf, AlignBuffers, AlignParams};
use crate::pval::compute_lambda;

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
    let center_idx = b.clusters[i].center.expect("b_compare: cluster has no center");
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
pub fn b_compare_parallel(
    b: &mut B,
    i: usize,
    err_mat: &[f64],
    ncol: usize,
    params: &AlignParams,
    greedy: bool,
) {
    let center_idx = b.clusters[i].center.expect("b_compare_parallel: cluster has no center");
    let center_reads = b.raws[center_idx].reads;
    let nraw = b.raws.len();
    let use_quals = b.use_quals;

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
    // benefit more. Smaller values (16, 8) were not meaningfully better.
    const PAR_MAX_LEN: usize = 32;
    let comps: Vec<(f64, u32, bool)> = (0..nraw)
        .into_par_iter()
        .with_max_len(PAR_MAX_LEN)
        .map_init(AlignBuffers::new, |buf, index| {
            let raw = &raws[index];
            if greedy && (raw.reads > center_reads || raw.lock) {
                let lambda = compute_lambda(raw, None, err_mat, ncol, use_quals);
                (lambda, u32::MAX, true)
            } else {
                let sub = sub_new_with_buf(&raws[center_idx], raw, params, buf);
                let lambda = compute_lambda(raw, sub.as_ref(), err_mat, ncol, use_quals);
                let hamming = sub.as_ref().map_or(u32::MAX, |s| s.nsubs() as u32);
                (lambda, hamming, false)
            }
        })
        .collect();

    // Serial post-processing: selectively store comparisons.
    let total_reads = b.reads as f64;
    for (index, (lambda, hamming, skipped)) in comps.into_iter().enumerate() {
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
}

// ---------------------------------------------------------------------------
// b_shuffle2
// ---------------------------------------------------------------------------

/// Move each Raw to the cluster that maximises its expected read count.
/// The center of a cluster may not be reassigned.
/// Returns `true` if any Raws were moved.
/// Equivalent to C++ `b_shuffle2`.
pub fn b_shuffle2(b: &mut B) -> bool {
    let nraw = b.raws.len();

    // Initialise best-E and best-comparison trackers from cluster 0.
    // During initialisation b_compare is always run on cluster 0 first, so
    // its comp vec contains an entry for every raw (in raw-index order).
    let mut emax: Vec<f64> = vec![f64::NEG_INFINITY; nraw];
    let mut compmax: Vec<Comparison> = vec![Comparison::default(); nraw];

    let c0_reads = b.clusters[0].reads as f64;
    for comp in &b.clusters[0].comp {
        let idx = comp.index as usize;
        emax[idx] = comp.lambda * c0_reads;
        compmax[idx] = comp.clone();
    }

    // Scan remaining clusters for better matches.
    for ci in 1..b.clusters.len() {
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
    let mut shuffled = false;
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
                shuffled = true;
            }
        }
    }
    shuffled
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
) -> Option<usize> {
    let nraw = b.raws.len() as f64;
    let init_center = b.clusters[0].center.expect("b_bud: cluster 0 has no center");

    // (cluster_idx, position_r, raw_index) for non-prior and prior minimums.
    let mut mini:       Option<(usize, usize, usize)> = None;
    let mut mini_prior: Option<(usize, usize, usize)> = None;
    let mut min_p       = b.raws[init_center].p;
    let mut min_p_prior = b.raws[init_center].p;
    let mut min_reads       = b.raws[init_center].reads;
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
            let lambda  = raw.comp.lambda;

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
                && (raw.p < min_p_prior
                    || (raw.p == min_p_prior && raw.reads > min_reads_prior))
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
    if p_a < b.omega_a {
        if let Some((ci, r, raw_idx)) = mini {
            // Capture pre-pop state.
            let expected   = b.raws[raw_idx].comp.lambda * b.clusters[ci].reads as f64;
            let birth_comp = b.raws[raw_idx].comp.clone();
            let birth_fold = b.raws[raw_idx].reads as f64 / expected.max(f64::MIN_POSITIVE);
            let nraw_total = b.raws.len() as u32;

            b.bi_pop_raw(ci, r);

            let mut new_bi = Bi::new(nraw_total);
            new_bi.birth_type = BirthType::Abundance;
            new_bi.birth_from = ci as u32;
            new_bi.birth_pval = p_a;
            new_bi.birth_fold = birth_fold;
            new_bi.birth_e    = expected;
            new_bi.birth_comp = birth_comp;

            let new_ci = b.add_cluster(new_bi);
            b.bi_add_raw(new_ci, raw_idx);
            b.assign_center(new_ci);

            if verbose {
                eprintln!(
                    "\n\t[cluster] Division (naive): Raw {raw_idx} from Bi {ci}, pA={p_a:.2e}"
                );
            }
            return Some(new_ci);
        }
    }

    // Prior-based bud.
    if p_p < b.omega_p {
        if let Some((ci, r, raw_idx)) = mini_prior {
            let expected   = b.raws[raw_idx].comp.lambda * b.clusters[ci].reads as f64;
            let birth_comp = b.raws[raw_idx].comp.clone();
            let birth_fold = b.raws[raw_idx].reads as f64 / expected.max(f64::MIN_POSITIVE);
            let nraw_total = b.raws.len() as u32;

            b.bi_pop_raw(ci, r);

            let mut new_bi = Bi::new(nraw_total);
            new_bi.birth_type = BirthType::Prior;
            new_bi.birth_pval = p_p;
            new_bi.birth_fold = birth_fold;
            new_bi.birth_e    = expected;
            new_bi.birth_comp = birth_comp;

            let new_ci = b.add_cluster(new_bi);
            b.bi_add_raw(new_ci, raw_idx);
            b.assign_center(new_ci);

            if verbose {
                eprintln!(
                    "\n\t[cluster] Division (prior): Raw {raw_idx} from Bi {ci}, pP={p_p:.2e}"
                );
            }
            return Some(new_ci);
        }
    }

    if verbose {
        let (raw_idx_str, ci_str) = match mini {
            Some((ci, r, _)) => {
                let raw_idx = b.clusters[ci].raws[r];
                (raw_idx.to_string(), ci.to_string())
            }
            None => (init_center.to_string(), String::from("0")),
        };
        eprintln!(
            "\n\t[cluster] No Division. Minimum pA={p_a:.2e} (Raw {raw_idx_str} in Bi {ci_str})."
        );
    }
    None
}
