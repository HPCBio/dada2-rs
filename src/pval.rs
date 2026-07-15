//! Abundance p-value and error-model lambda calculations.
// Function names like calc_pA / get_pA intentionally match the C++ source.
#![allow(non_snake_case)]
//!
//! Ported from `pval.cpp`.
//!
//! The Poisson CDF previously called via `Rcpp::ppois` is replaced by
//! `statrs::distribution::Poisson`, which uses the same regularised
//! incomplete gamma function that R's `ppois` calls internally.

use statrs::distribution::{DiscreteCDF, Poisson};

use crate::containers::{B, Raw, Sub};

/// Minimum value of the conditioning normaliser below which the second-order
/// Taylor approximation `E - E²/2` is used instead of `1 - exp(-E)`.
/// Matches C++ `TAIL_APPROX_CUTOFF` in `dada.h`.
const TAIL_APPROX_CUTOFF: f64 = 1e-7;

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

/// Update abundance p-values for every Raw in the partition.
///
/// For each cluster whose `update_e` flag is set, recomputes `raw.p` for all
/// member Raws and clears the flag. When `greedy` is true, also locks Raws
/// whose expected read count from the cluster center alone already exceeds
/// their observed count (preventing them from budding a new cluster).
///
/// Returns the number of Raws whose `p` was recomputed this call — i.e. the
/// members of every cluster with `update_e` set. This is the per-round p-churn
/// that gated the b_bud incremental follow-up (issue #85).
///
/// While repricing each dirty cluster's members, this also refreshes that
/// cluster's cached best budding candidate (`Bi::bud_min` / `bud_min_prior`),
/// applying the same eligibility filters and tie-break as `b_bud`'s scan. The
/// member loop already touches exactly these Raws with `bi_reads` in hand, so
/// the cache is maintained at no extra memory traffic — and `b_bud` then
/// combines the per-cluster minima in O(nclusters) instead of rescanning every
/// Raw (issue #85). `min_fold`/`min_hamming`/`min_abund` are `b_bud`'s
/// candidate filters, threaded through so the cache matches its scan exactly.
///
/// Equivalent to C++ `b_p_update`.
pub fn b_p_update(
    b: &mut B,
    greedy: bool,
    detect_singletons: bool,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
) -> u64 {
    use crate::containers::BudCand;
    let mut raws_repriced = 0u64;
    for i in 0..b.clusters.len() {
        if b.clusters[i].update_e {
            // Clone indices to avoid holding a shared borrow on b.clusters
            // while mutating b.raws.
            let members: Vec<usize> = b.clusters[i].raws.clone();
            let bi_reads = b.clusters[i].reads;
            // Refresh this cluster's cached bud candidate in the same pass.
            // Position 0 is skipped, mirroring b_bud's `for r in 1..len`.
            let mut bud_min: Option<BudCand> = None;
            let mut bud_min_prior: Option<BudCand> = None;
            for (r, &raw_idx) in members.iter().enumerate() {
                let p = get_pA(
                    b.raws[raw_idx].reads,
                    b.raws[raw_idx].prior,
                    b.raws[raw_idx].comp.lambda,
                    b.raws[raw_idx].comp.hamming,
                    bi_reads,
                    detect_singletons,
                );
                b.raws[raw_idx].p = p;
                raws_repriced += 1;

                if r == 0 {
                    continue; // center slot: never a bud candidate
                }
                let raw = &b.raws[raw_idx];
                if raw.reads < min_abund || raw.comp.hamming < min_hamming {
                    continue;
                }
                let fold_ok = min_fold <= 1.0
                    || raw.reads as f64 >= min_fold * raw.comp.lambda * bi_reads as f64;
                if !fold_ok {
                    continue;
                }
                let cand = BudCand {
                    p,
                    reads: raw.reads,
                    r,
                    raw_idx,
                };
                // Ascending r with strict replacement keeps the first-in-scan
                // (min p, then max reads, then lowest r) — identical to b_bud.
                if bud_min.is_none_or(|m| p < m.p || (p == m.p && cand.reads > m.reads)) {
                    bud_min = Some(cand);
                }
                if raw.prior
                    && bud_min_prior.is_none_or(|m| p < m.p || (p == m.p && cand.reads > m.reads))
                {
                    bud_min_prior = Some(cand);
                }
            }
            b.clusters[i].bud_min = bud_min;
            b.clusters[i].bud_min_prior = bud_min_prior;
            b.clusters[i].update_e = false;
        }

        if greedy && b.clusters[i].check_locks {
            let center_idx = b.clusters[i].center;
            let center_reads = center_idx.map_or(0, |ci| b.raws[ci].reads);
            let members: Vec<usize> = b.clusters[i].raws.clone();
            for raw_idx in members {
                // Lock if the center alone expects more reads than observed.
                let e_from_center = center_reads as f64 * b.raws[raw_idx].comp.lambda;
                if e_from_center > b.raws[raw_idx].reads as f64 {
                    b.raws[raw_idx].lock = true;
                }
                // Always lock the center itself.
                if Some(raw_idx) == center_idx {
                    b.raws[raw_idx].lock = true;
                }
            }
            b.clusters[i].check_locks = false;
        }
    }
    raws_repriced
}

/// Non-mutating serial equivalent of `b_bud`'s candidate selection: scans every
/// non-position-0 raw across all clusters and returns the abundance and prior
/// minima it would pick, as `(mini, min_p, min_reads, mini_prior, min_p_prior,
/// min_reads_prior)` where each `mini*` is `Option<(ci, r, raw_idx)>`. Used only
/// by `b_bud`'s debug cross-check against the incremental cache (issue #85).
#[cfg(debug_assertions)]
#[allow(clippy::type_complexity)]
pub fn b_bud_scan_select(
    b: &B,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
) -> (
    Option<(usize, usize, usize)>,
    f64,
    u32,
    Option<(usize, usize, usize)>,
    f64,
    u32,
) {
    let init_center = b.clusters[0]
        .center
        .expect("b_bud_scan_select: cluster 0 has no center");
    let mut mini: Option<(usize, usize, usize)> = None;
    let mut mini_prior: Option<(usize, usize, usize)> = None;
    let mut min_p = b.raws[init_center].p;
    let mut min_p_prior = b.raws[init_center].p;
    let mut min_reads = b.raws[init_center].reads;
    let mut min_reads_prior = b.raws[init_center].reads;

    for ci in 0..b.clusters.len() {
        for r in 1..b.clusters[ci].raws.len() {
            let raw_idx = b.clusters[ci].raws[r];
            let raw = &b.raws[raw_idx];
            if raw.reads < min_abund || raw.comp.hamming < min_hamming {
                continue;
            }
            let fold_ok = min_fold <= 1.0
                || raw.reads as f64 >= min_fold * raw.comp.lambda * b.clusters[ci].reads as f64;
            if !fold_ok {
                continue;
            }
            if raw.p < min_p || (raw.p == min_p && raw.reads > min_reads) {
                mini = Some((ci, r, raw_idx));
                min_p = raw.p;
                min_reads = raw.reads;
            }
            if raw.prior
                && (raw.p < min_p_prior || (raw.p == min_p_prior && raw.reads > min_reads_prior))
            {
                mini_prior = Some((ci, r, raw_idx));
                min_p_prior = raw.p;
                min_reads_prior = raw.reads;
            }
        }
    }
    (
        mini,
        min_p,
        min_reads,
        mini_prior,
        min_p_prior,
        min_reads_prior,
    )
}

/// Abundance p-value: P(X ≥ `reads` | Poisson(λ = `e_reads`)).
///
/// When `prior` is false the p-value is conditioned on the sequence being
/// present (i.e. at least one read), normalising by `1 - exp(-e_reads)`.
/// A second-order Taylor expansion replaces the normaliser when it would
/// underflow below `TAIL_APPROX_CUTOFF`.
///
/// # Panics
/// Panics in debug builds if `reads == 0` (callers must guard this case).
/// Equivalent to C++ `calc_pA`.
pub fn calc_pA(reads: u32, e_reads: f64, prior: bool) -> f64 {
    debug_assert!(reads > 0, "calc_pA: reads must be > 0");

    // R uses ppois(reads-1, e_reads, lower.tail=FALSE), which under the hood
    // calls the regularised lower incomplete gamma function. statrs's `sf`
    // does the same nominally — sf(x) = gamma_lr(x+1, λ) — but its
    // gamma_lr implementation underflows to 0 for very small λ, while R's
    // pgamma stays accurate down to ~1e-300. For tiny λ we fall back to the
    // direct upper-tail Poisson series in log space:
    //   P(X >= k | λ) = e^{-λ} λ^k/k! * Σ_{j>=0} λ^j / ∏_{i=1..j}(k+i)
    // which is dominated by the leading e^{-λ} λ^k/k! for small λ.
    let pois = match Poisson::new(e_reads) {
        Ok(p) => p,
        Err(_) => return 1.0, // degenerate: e_reads <= 0
    };
    let mut pval = pois.sf((reads - 1) as u64);
    if pval == 0.0 && e_reads > 0.0 {
        pval = poisson_upper_tail_direct(reads, e_reads);
    }

    if prior {
        return pval;
    }

    // Condition on the sequence being present (at least one read observed).
    let norm = 1.0 - (-e_reads).exp();
    let norm = if norm < TAIL_APPROX_CUTOFF {
        // 2nd-order Taylor: 1 - e^{-E} ≈ E - E²/2 for small E.
        e_reads - 0.5 * e_reads * e_reads
    } else {
        norm
    };
    pval / norm
}

/// Compute lambda: the probability under the error model that `raw`'s
/// sequence arose from the reference sequence described by `sub`.
///
/// `err_mat` is a flat row-major matrix of shape 16 × `ncol` where rows
/// index the 16 nucleotide transitions (ref_nt × 4 + query_nt, with A=0,
/// C=1, G=2, T=3) and columns index the rounded quality score.
///
/// Returns `0.0` when `sub` is `None` (sequence was outside the k-mer
/// distance threshold and was not aligned).
/// Equivalent to C++ `compute_lambda_ts`.
pub fn compute_lambda(
    raw: &Raw,
    sub: Option<&Sub>,
    err_mat: &[f64],
    ncol: usize,
    use_quals: bool,
) -> f64 {
    let sub = match sub {
        Some(s) => s,
        None => return 0.0,
    };

    let len = raw.len();

    // Initialise per-position transition and quality-index vectors.
    // Default: diagonal match (nti * 4 + nti for both ref and query).
    let mut tvec = vec![0usize; len];
    let mut qind = vec![0usize; len];

    for pos in 0..len {
        let nti = raw.seq[pos]
            .checked_sub(1)
            .filter(|&n| n < 4)
            .expect("non-ACGT nucleotide in compute_lambda") as usize;
        tvec[pos] = nti * 4 + nti;
        qind[pos] = if use_quals {
            raw.qual.as_ref().map_or(0, |q| q[pos] as usize)
        } else {
            0
        };
    }

    // Override transition at each substitution position.
    for s in 0..sub.nsubs() {
        let pos0 = sub.pos[s] as usize;
        let pos1 = sub.map[pos0] as usize;

        debug_assert!(
            pos0 < sub.len0 as usize,
            "sub pos0 {pos0} >= len0 {}",
            sub.len0
        );
        debug_assert!(pos1 < len, "sub pos1 {pos1} >= raw len {len}");

        let nti0 = sub.nt0[s].saturating_sub(1) as usize;
        let nti1 = sub.nt1[s].saturating_sub(1) as usize;
        tvec[pos1] = nti0 * 4 + nti1;
    }

    // Lambda = product of error probabilities across all positions.
    let lambda: f64 = (0..len)
        .map(|pos| err_mat[tvec[pos] * ncol + qind[pos]])
        .product();

    debug_assert!(
        (0.0..=1.0).contains(&lambda),
        "lambda {lambda} outside [0,1]"
    );
    lambda
}

/// Direct upper-tail Poisson computation in log space, robust to extreme
/// underflow.
///
/// Used as a fallback when `statrs::Poisson::sf` returns exactly 0 due to
/// `gamma_lr` losing precision for very small λ. The leading term
/// `e^{-λ} λ^reads / reads!` dominates for small λ; the series correction
/// `1 + λ/(k+1) + λ²/((k+1)(k+2)) + ...` converges quickly.
fn poisson_upper_tail_direct(reads: u32, lambda: f64) -> f64 {
    debug_assert!(reads > 0);
    debug_assert!(lambda > 0.0);

    let k = reads as f64;
    let log_lambda = lambda.ln();
    let log_kfact: f64 = (1..=reads).map(|i| (i as f64).ln()).sum();
    let log_leading = -lambda + k * log_lambda - log_kfact;

    let mut series_sum: f64 = 1.0;
    let mut term: f64 = 1.0;
    for j in 1..10000u32 {
        term *= lambda / (k + j as f64);
        if !term.is_finite() || term <= 0.0 {
            break;
        }
        let new_sum = series_sum + term;
        if new_sum == series_sum {
            break;
        }
        series_sum = new_sum;
    }

    (log_leading + series_sum.ln()).exp()
}

/// Self-production probability: the probability that a sequence is produced
/// from itself given a compact 4×4 match-probability matrix.
///
/// Returns the product of `err[nti][nti]` (match diagonal) over all positions.
/// Equivalent to C++ `get_self`.
#[allow(dead_code)]
pub fn get_self(seq: &[u8], err: &[[f64; 4]; 4]) -> f64 {
    seq.iter().fold(1.0, |acc, &nt| {
        let nti = (nt as usize).saturating_sub(1).min(3);
        acc * err[nti][nti]
    })
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Per-raw p-value calculation, handling all sentinel cases before calling
/// `calc_pA`. Factored out to allow use from `b_p_update` without holding
/// simultaneous borrows on both `B.raws` and `B.clusters`.
/// Equivalent to C++ `get_pA`.
fn get_pA(
    reads: u32,
    prior: bool,
    lambda: f64,
    hamming: u32,
    bi_reads: u32,
    detect_singletons: bool,
) -> f64 {
    if reads == 1 && !prior && !detect_singletons {
        // Singleton: no abundance p-value is applied.
        return 1.0;
    }
    if hamming == 0 {
        // Cluster center (or exact match): always valid.
        return 1.0;
    }
    if lambda == 0.0 {
        // Zero expected reads: reject unconditionally.
        return 0.0;
    }
    let e_reads = lambda * bi_reads as f64;
    calc_pA(reads, e_reads, prior || detect_singletons)
}
