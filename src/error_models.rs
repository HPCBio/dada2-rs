//! Error model estimation.
//!
//! Ports the error model functions from R's `errorModels.R`:
//! - `loessErrfun`          → [`loess_errfun`]
//! - `makeBinnedQualErrfun` → [`binned_qual_errfun`]
//! - `PacBioErrfun`         → [`pacbio_errfun`]
//! - `noqualErrfun`         → [`noqual_errfun`]
//! - `inflateErr`           → [`inflate_err`]
//! - `accumulateTrans`      → [`accumulate_trans`]
//! - `getBadBases`          → [`get_bad_bases`]
//! - `isBadBaseFP`          → [`is_bad_base_fp`]
//!
//! ## Transition matrix layout
//!
//! All 16 × `nq` matrices are flat, row-major `Vec<T>`:
//! - Row `nti * 4 + ntj` is the transition from nucleotide `nti` to `ntj`
//!   (A=0, C=1, G=2, T=3).
//! - Columns index quality scores, with column 0 corresponding to Q0.

use std::collections::{HashMap, HashSet};

use statrs::distribution::{DiscreteCDF, Poisson};

const MAX_ERROR_RATE: f64 = 0.25;
const MIN_ERROR_RATE: f64 = 1e-7;

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Gaussian elimination with partial pivoting.
///
/// Solves `A * x = b` where `a` is the row-major `n×n` matrix stored in a
/// flat slice and `b` is the right-hand-side vector.
/// Both slices are modified in place.  Returns `Some(x)` on success, `None`
/// if the system is (numerically) singular.
fn solve_linear(a: &mut [f64], b: &mut [f64], n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        // Find the row with the largest absolute value in this column (pivot).
        let pivot_row = (col..n).max_by(|&r1, &r2| {
            a[r1 * n + col]
                .abs()
                .partial_cmp(&a[r2 * n + col].abs())
                .unwrap()
        })?;

        if pivot_row != col {
            for k in 0..n {
                a.swap(col * n + k, pivot_row * n + k);
            }
            b.swap(col, pivot_row);
        }

        let pivot = a[col * n + col];
        if pivot.abs() < 1e-12 {
            return None;
        }

        for row in (col + 1)..n {
            let f = a[row * n + col] / pivot;
            a[row * n + col] = 0.0;
            for k in (col + 1)..n {
                let v = a[col * n + k] * f;
                a[row * n + k] -= v;
            }
            let bv = b[col] * f;
            b[row] -= bv;
        }
    }

    // Back substitution.
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i * n + j] * x[j];
        }
        let diag = a[i * n + i];
        if diag.abs() < 1e-12 {
            return None;
        }
        x[i] = s / diag;
    }
    Some(x)
}

/// Locally-weighted polynomial regression (LOESS).
///
/// Mirrors R's `loess(y ~ x, data, weights=w)` with `span = 0.75` and
/// `degree = 2`.  Observations with non-finite `y` or zero weight are
/// excluded from fitting.  Predictions outside the valid data range are
/// extrapolated flat from the nearest finite prediction.
///
/// Returns a `Vec<Option<f64>>` aligned to `xs`; `None` at a position means
/// the local fit could not be computed there.
fn loess_predict(
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    span: f64,
    degree: usize,
) -> Vec<Option<f64>> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert_eq!(xs.len(), weights.len());

    let n = xs.len();

    // Collect indices of valid observations (finite y, positive weight).
    let valid: Vec<usize> = (0..n)
        .filter(|&i| ys[i].is_finite() && weights[i] > 0.0)
        .collect();
    let nv = valid.len();

    // Clamp polynomial degree to what the data allows.
    let eff_degree = degree.min(nv.saturating_sub(1));
    if nv <= eff_degree {
        return vec![None; n];
    }

    let n_local = ((span * nv as f64).ceil() as usize)
        .max(eff_degree + 1)
        .min(nv);
    let p = eff_degree + 1; // number of polynomial coefficients

    // R's `predict.loess` returns NA for x outside the fitted-data range;
    // `loessErrfun` then flat-fills via `pred[minrli]` / `pred[maxrli]`. We
    // match that: out-of-range x → None here so `extrapolate_flat` handles
    // the boundary fill, instead of a polynomial extrapolation that drifts.
    let (x_min, x_max) = valid
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &i| {
            (lo.min(xs[i]), hi.max(xs[i]))
        });

    (0..n)
        .map(|pred_idx| {
            let x0 = xs[pred_idx];

            if x0 < x_min || x0 > x_max {
                return None;
            }

            // Sort valid observations by distance to x0.
            let mut dists: Vec<(usize, f64)> =
                valid.iter().map(|&i| (i, (xs[i] - x0).abs())).collect();
            dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            let max_dist = dists[n_local - 1].1;

            // Combined weight = tricubic(distance) × user-supplied weight.
            let ws: Vec<(usize, f64)> = dists[..n_local]
                .iter()
                .map(|&(i, d)| {
                    let u = if max_dist > 0.0 { d / max_dist } else { 0.0 };
                    let tri = if u < 1.0 {
                        (1.0 - u.powi(3)).powi(3)
                    } else {
                        0.0
                    };
                    (i, tri * weights[i])
                })
                .filter(|&(_, w)| w > 0.0)
                .collect();

            if ws.len() < p {
                return None;
            }

            // Build normal equations X^T W X a = X^T W y.
            let mut xtx = vec![0.0f64; p * p];
            let mut xty = vec![0.0f64; p];

            for &(i, w) in &ws {
                let xi = xs[i];
                let yi = ys[i];

                // row = [1, xi, xi^2, …]
                let mut row = vec![1.0f64; p];
                let mut xpow = xi;
                for j in 1..p {
                    row[j] = xpow;
                    xpow *= xi;
                }

                for j in 0..p {
                    xty[j] += w * row[j] * yi;
                    for k in 0..p {
                        xtx[j * p + k] += w * row[j] * row[k];
                    }
                }
            }

            let coeffs = solve_linear(&mut xtx, &mut xty, p)?;

            // Evaluate polynomial at x0.
            let mut pred = 0.0f64;
            let mut xpow = 1.0f64;
            for j in 0..p {
                pred += coeffs[j] * xpow;
                xpow *= x0;
            }
            Some(pred)
        })
        .collect()
}

/// Apply flat extrapolation at the ends of a LOESS prediction vector.
///
/// Fills `None` entries at the low end with the first finite value and at the
/// high end with the last finite value.  Returns the filled vector.
fn extrapolate_flat(raw: Vec<Option<f64>>, n: usize) -> Vec<f64> {
    let valid: Vec<(usize, f64)> = raw
        .iter()
        .enumerate()
        .filter_map(|(i, p)| p.map(|v| (i, v)))
        .collect();

    if valid.is_empty() {
        return vec![f64::NAN; n];
    }

    let (min_i, min_v) = *valid.first().unwrap();
    let (max_i, max_v) = *valid.last().unwrap();

    let mut out = vec![f64::NAN; n];
    for &(i, v) in &valid {
        out[i] = v;
    }
    for i in 0..min_i {
        out[i] = min_v;
    }
    for i in (max_i + 1)..n {
        out[i] = max_v;
    }
    out
}

/// Clamp a rate to the canonical `[MIN_ERROR_RATE, MAX_ERROR_RATE]` interval.
#[inline]
fn clamp_rate(r: f64) -> f64 {
    r.clamp(MIN_ERROR_RATE, MAX_ERROR_RATE)
}

/// Expand the 12 off-diagonal estimated rates into the full 16 × `nq` error matrix.
///
/// The 12 off-diagonal rows (0-indexed) are ordered by iterating source
/// nucleotide A→C→G→T and, for each, by destination nucleotide A→C→G→T
/// (skipping the diagonal):
///
/// ```text
/// 0: A→C,  1: A→G,  2: A→T,
/// 3: C→A,  4: C→G,  5: C→T,
/// 6: G→A,  7: G→C,  8: G→T,
/// 9: T→A, 10: T→C, 11: T→G
/// ```
///
/// Diagonal (self-transition) entries are inserted as
/// `1 − sum(off-diagonals for that source nucleotide)`, mirroring the R code.
fn expand_err_matrix(off_diag: &[f64], nq: usize) -> Vec<f64> {
    debug_assert_eq!(off_diag.len(), 12 * nq);

    let mut err = vec![0.0f64; 16 * nq];

    // (off_row_start, [dest_full_rows], diag_full_row)
    //   full_row = nti * 4 + ntj
    let groups: [(usize, [usize; 3], usize); 4] = [
        (0, [1, 2, 3], 0),     // A: A→C, A→G, A→T; diag=A→A (row 0)
        (3, [4, 6, 7], 5),     // C: C→A, C→G, C→T; diag=C→C (row 5)
        (6, [8, 9, 11], 10),   // G: G→A, G→C, G→T; diag=G→G (row 10)
        (9, [12, 13, 14], 15), // T: T→A, T→C, T→G; diag=T→T (row 15)
    ];

    for (off_start, full_rows, diag_row) in groups {
        for q in 0..nq {
            let vals = [
                off_diag[(off_start) * nq + q],
                off_diag[(off_start + 1) * nq + q],
                off_diag[(off_start + 2) * nq + q],
            ];
            for (k, &fr) in full_rows.iter().enumerate() {
                err[fr * nq + q] = vals[k];
            }
            err[diag_row * nq + q] = 1.0 - vals[0] - vals[1] - vals[2];
        }
    }

    err
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Estimate error rates using LOESS smoothing of observed transition counts.
///
/// Mirrors R's `loessErrfun(trans)`.  For each of the 12 off-diagonal
/// transition types, fits `log10((errs + 1) / tot) ~ quality` (weighted by
/// `tot`) with `span = 0.75` and `degree = 2`, then converts predictions back
/// with `10^pred`.  Predictions are extrapolated flat outside the data range.
/// All rates are clamped to `[1e-7, 0.25]`.  Diagonal entries are filled as
/// `1 − sum(off-diagonals)`.
///
/// # Arguments
/// - `trans`: flat 16 × `nq` row-major matrix of transition counts
/// - `qual_scores`: quality score value for each column (length `nq`)
///
/// # Returns
/// Flat 16 × `nq` row-major error rate matrix.
pub fn loess_errfun(trans: &[u32], qual_scores: &[f64]) -> Vec<f64> {
    let nq = qual_scores.len();
    assert_eq!(trans.len(), 16 * nq, "trans length must be 16 * nq");

    let mut off_diag = vec![0.0f64; 12 * nq];
    let mut off_row = 0usize;

    for nti in 0..4usize {
        // Total counts across all destinations for this source nucleotide.
        let tot: Vec<f64> = (0..nq)
            .map(|q| {
                (0..4usize)
                    .map(|ntj| trans[(nti * 4 + ntj) * nq + q] as f64)
                    .sum()
            })
            .collect();

        for ntj in 0..4usize {
            if nti == ntj {
                continue;
            }

            let errs: Vec<f64> = (0..nq)
                .map(|q| trans[(nti * 4 + ntj) * nq + q] as f64)
                .collect();

            // rlogp = log10((errs+1)/tot); NA where tot==0.
            let rlogp: Vec<f64> = (0..nq)
                .map(|q| {
                    if tot[q] == 0.0 {
                        f64::NAN
                    } else {
                        ((errs[q] + 1.0) / tot[q]).log10()
                    }
                })
                .collect();

            let raw_pred = loess_predict(qual_scores, &rlogp, &tot, 0.75, 2);
            let filled = extrapolate_flat(raw_pred, nq);

            for q in 0..nq {
                let rate = if filled[q].is_finite() {
                    10.0f64.powf(filled[q])
                } else {
                    MIN_ERROR_RATE
                };
                off_diag[off_row * nq + q] = clamp_rate(rate);
            }
            off_row += 1;
        }
    }

    expand_err_matrix(&off_diag, nq)
}

/// Estimate error rates ignoring quality scores.
///
/// Mirrors R's `noqualErrfun(trans, pseudocount=1)`.  Aggregates all
/// transition counts across quality scores, adds `pseudocount` to each row,
/// then estimates rates as `count / total_from_that_nucleotide`.  The same
/// rate is broadcast to every quality-score column.
///
/// # Arguments
/// - `trans`: flat 16 × `nq` row-major matrix of transition counts
/// - `nq`: number of quality-score columns
/// - `pseudocount`: added to each row sum before computing rates (R default: 1)
///
/// # Returns
/// Flat 16 × `nq` row-major error rate matrix (all columns identical).
pub fn noqual_errfun(trans: &[u32], nq: usize, pseudocount: f64) -> Vec<f64> {
    assert_eq!(trans.len(), 16 * nq, "trans length must be 16 * nq");

    // obs[r] = sum over quality columns + pseudocount
    let obs: Vec<f64> = (0..16)
        .map(|r| (0..nq).map(|q| trans[r * nq + q] as f64).sum::<f64>() + pseudocount)
        .collect();

    let mut off_diag = vec![0.0f64; 12 * nq];
    let mut off_row = 0usize;

    for nti in 0..4usize {
        let tot_init: f64 = (0..4usize).map(|ntj| obs[nti * 4 + ntj]).sum();

        for ntj in 0..4usize {
            if nti == ntj {
                continue;
            }

            let rate = if tot_init > 0.0 {
                clamp_rate(obs[nti * 4 + ntj] / tot_init)
            } else {
                MIN_ERROR_RATE
            };

            // Same rate broadcast to every quality column.
            for q in 0..nq {
                off_diag[off_row * nq + q] = rate;
            }
            off_row += 1;
        }
    }

    expand_err_matrix(&off_diag, nq)
}

/// Estimate error rates using piecewise linear interpolation between binned quality scores.
///
/// Mirrors the function returned by R's `makeBinnedQualErrfun(binnedQ)(trans)`.
/// For each off-diagonal transition, the observed rate at each anchor quality
/// score is linearly interpolated between adjacent anchors.  Rates outside the
/// anchor range are extrapolated flat.
///
/// # Arguments
/// - `trans`: flat 16 × `nq` row-major transition count matrix
/// - `qual_scores`: quality score values — **must be consecutive integers 0, 1, 2, …**
/// - `binned_quals`: anchor quality score values (sorted ascending)
///
/// # Errors
/// Returns `Err` if quality scores are not 0-based consecutive integers, or if
/// the observed data range falls outside the anchor range.
pub fn binned_qual_errfun(
    trans: &[u32],
    qual_scores: &[f64],
    binned_quals: &[f64],
) -> Result<Vec<f64>, String> {
    let nq = qual_scores.len();
    assert_eq!(trans.len(), 16 * nq, "trans length must be 16 * nq");

    // Verify that qual_scores are 0, 1, 2, …
    for (i, &q) in qual_scores.iter().enumerate() {
        if (q - i as f64).abs() > 1e-9 {
            return Err(format!(
                "Unexpected Q score series: expected {i} but got {q}"
            ));
        }
    }

    // Determine min/max quality scores that have any observations.
    let col_totals: Vec<f64> = (0..nq)
        .map(|q| (0..16usize).map(|r| trans[r * nq + q] as f64).sum())
        .collect();

    let observed_qs: Vec<f64> = qual_scores
        .iter()
        .enumerate()
        .filter(|&(q, _)| col_totals[q] > 0.0)
        .map(|(_, &v)| v)
        .collect();

    if observed_qs.is_empty() {
        return Err("No observed transitions".to_string());
    }

    let qmax = observed_qs
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let qmin = observed_qs.iter().cloned().fold(f64::INFINITY, f64::min);
    let bmax = binned_quals
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let bmin = binned_quals.iter().cloned().fold(f64::INFINITY, f64::min);

    if qmax > bmax {
        return Err(format!(
            "Input data contains quality score {qmax} higher than max binned value {bmax}"
        ));
    }
    if qmin < bmin {
        return Err(format!(
            "Input data contains quality score {qmin} lower than min binned value {bmin}"
        ));
    }

    let mut off_diag = vec![0.0f64; 12 * nq];
    let mut off_row = 0usize;

    for nti in 0..4usize {
        let tot: Vec<f64> = (0..nq)
            .map(|q| {
                (0..4usize)
                    .map(|ntj| trans[(nti * 4 + ntj) * nq + q] as f64)
                    .sum()
            })
            .collect();

        for ntj in 0..4usize {
            if nti == ntj {
                continue;
            }

            let errs: Vec<f64> = (0..nq)
                .map(|q| trans[(nti * 4 + ntj) * nq + q] as f64)
                .collect();

            // Observed rate per quality score (NaN where tot == 0).
            let p: Vec<f64> = (0..nq)
                .map(|q| {
                    if tot[q] > 0.0 {
                        errs[q] / tot[q]
                    } else {
                        f64::NAN
                    }
                })
                .collect();

            // Piecewise linear interpolation between adjacent anchor pairs.
            let mut pred = vec![f64::NAN; nq];
            for i in 0..binned_quals.len().saturating_sub(1) {
                let lo_q = binned_quals[i] as usize;
                let hi_q = binned_quals[i + 1] as usize;

                if lo_q >= nq || hi_q >= nq {
                    continue;
                }

                let lo_p = p[lo_q];
                let hi_p = p[hi_q];

                if lo_p.is_finite() && hi_p.is_finite() {
                    let steps = hi_q - lo_q + 1;
                    for k in 0..steps {
                        let t = if steps > 1 {
                            k as f64 / (steps - 1) as f64
                        } else {
                            0.0
                        };
                        pred[lo_q + k] = lo_p * (1.0 - t) + hi_p * t;
                    }
                }
            }

            // Flat extrapolation outside the anchor range.
            let raw_opt: Vec<Option<f64>> = pred
                .iter()
                .map(|&v| if v.is_finite() { Some(v) } else { None })
                .collect();
            let filled = extrapolate_flat(raw_opt, nq);

            for q in 0..nq {
                off_diag[off_row * nq + q] = if filled[q].is_finite() {
                    clamp_rate(filled[q])
                } else {
                    MIN_ERROR_RATE
                };
            }
            off_row += 1;
        }
    }

    Ok(expand_err_matrix(&off_diag, nq))
}

/// Estimate error rates for PacBio CCS data.
///
/// Mirrors R's `PacBioErrfun(trans)`.  If the last quality-score column is
/// Q93, it applies LOESS fitting to quality scores 0–92 and a maximum-
/// likelihood estimate `(count + 1) / (group_total + 4)` for Q93.  Otherwise
/// falls back to [`loess_errfun`].
///
/// # Arguments
/// - `trans`: flat 16 × `nq` row-major transition count matrix
/// - `qual_scores`: quality score values for each column
pub fn pacbio_errfun(trans: &[u32], qual_scores: &[f64]) -> Vec<f64> {
    let nq = qual_scores.len();
    assert_eq!(trans.len(), 16 * nq, "trans length must be 16 * nq");

    if nq > 0 && (qual_scores[nq - 1] - 93.0).abs() < 1e-9 {
        // Fit LOESS on columns 0..(nq-1).
        let sub_nq = nq - 1;
        let trans_sub: Vec<u32> = (0..16)
            .flat_map(|r| (0..sub_nq).map(move |q| trans[r * nq + q]))
            .collect();
        let qs_sub = &qual_scores[..sub_nq];
        let err_sub = loess_errfun(&trans_sub, qs_sub); // 16 × sub_nq

        // MLE for Q93: (count + 1) / (group_total + 4).
        let q93 = nq - 1;
        let mut err93 = vec![0.0f64; 16];
        for nti in 0..4usize {
            let tot93: f64 = (0..4usize)
                .map(|ntj| trans[(nti * 4 + ntj) * nq + q93] as f64)
                .sum();
            for ntj in 0..4usize {
                let c = trans[(nti * 4 + ntj) * nq + q93] as f64;
                err93[nti * 4 + ntj] = (c + 1.0) / (tot93 + 4.0);
            }
        }

        // Combine: err_sub (16 × sub_nq) followed by err93 column.
        let mut full_err = vec![0.0f64; 16 * nq];
        for r in 0..16 {
            for q in 0..sub_nq {
                full_err[r * nq + q] = err_sub[r * sub_nq + q];
            }
            full_err[r * nq + q93] = err93[r];
        }
        full_err
    } else {
        loess_errfun(trans, qual_scores)
    }
}

/// Inflate error rates by a multiplicative factor with saturation.
///
/// Mirrors R's `inflateErr(err, inflation, inflateSelfTransitions=FALSE)`.
/// Applies `new = rate * k / (1 + (k−1) * rate)` to prevent rates exceeding
/// 1.0.  By default only off-diagonal (substitution) rates are inflated.
///
/// # Arguments
/// - `err`: flat 16 × `nq` error rate matrix
/// - `nq`: number of quality-score columns
/// - `inflation`: inflation factor
/// - `inflate_self`: if `true`, also inflate diagonal (self-transition) rates
pub fn inflate_err(err: &[f64], nq: usize, inflation: f64, inflate_self: bool) -> Vec<f64> {
    assert_eq!(err.len(), 16 * nq, "err length must be 16 * nq");

    let mut out = err.to_vec();
    for nti in 0..4usize {
        for ntj in 0..4usize {
            if nti == ntj && !inflate_self {
                continue;
            }
            let r = nti * 4 + ntj;
            for q in 0..nq {
                let rate = err[r * nq + q];
                out[r * nq + q] = (rate * inflation) / (1.0 + (inflation - 1.0) * rate);
            }
        }
    }
    out
}

/// Accumulate transition count matrices with potentially different numbers of quality columns.
///
/// Mirrors R's `accumulateTrans(trans)`.  Sums corresponding cells across all
/// input matrices.  The output has `max_nq` columns where `max_nq` is the
/// maximum `nq` seen in the inputs; columns not present in a shorter matrix
/// are implicitly zero.
///
/// # Arguments
/// - `matrices`: slice of `(flat_trans_slice, nq)` pairs
///
/// # Returns
/// `(accumulated_trans, max_nq)`
pub fn accumulate_trans(matrices: &[(&[u32], usize)]) -> (Vec<u32>, usize) {
    if matrices.is_empty() {
        return (vec![0u32; 16], 1);
    }

    let max_nq = matrices.iter().map(|(_, nq)| *nq).max().unwrap_or(1);
    let mut result = vec![0u32; 16 * max_nq];

    for &(mat, nq) in matrices {
        debug_assert_eq!(mat.len(), 16 * nq);
        for r in 0..16 {
            for q in 0..nq {
                result[r * max_nq + q] += mat[r * nq + q];
            }
        }
    }

    (result, max_nq)
}

/// Identify bad base positions that disproportionately drive cluster births.
///
/// Mirrors R's `getBadBases(clust, birth_subs, omegaB, minOccurence)`.
///
/// A position is flagged as "bad" if it appears in the birth substitutions of
/// more 1-Hamming-distance clusters than expected by chance (Poisson null with
/// rate `n_ones / seqlen`), after Bonferroni correction for `seqlen` tests,
/// and it occurs at least `min_occurrence` times.
///
/// # Arguments
/// - `seqlen`: length of the reference sequence
/// - `birth_subs_pos_1ham`: birth-substitution positions (0-indexed) drawn
///   from clusters with `birth_ham == 1`; one entry per cluster (since ham==1
///   implies exactly one birth substitution)
/// - `omega_b`: Bonferroni-corrected p-value threshold (R default: `1e-20`)
/// - `min_occurrence`: minimum count per position to be considered (R default: 4)
///
/// # Returns
/// Sorted vector of 0-indexed bad base positions.
pub fn get_bad_bases(
    seqlen: usize,
    birth_subs_pos_1ham: &[u16],
    omega_b: f64,
    min_occurrence: usize,
) -> Vec<u16> {
    if birth_subs_pos_1ham.is_empty() || seqlen == 0 {
        return Vec::new();
    }

    // Count occurrences of each position.
    let mut counts: HashMap<u16, usize> = HashMap::new();
    for &pos in birth_subs_pos_1ham {
        *counts.entry(pos).or_insert(0) += 1;
    }

    // Poisson null: lambda = n_ones / seqlen.
    // For ham==1, each cluster contributes exactly one position, so
    // n_ones == birth_subs_pos_1ham.len().
    let n_ones = birth_subs_pos_1ham.len();
    let lambda = n_ones as f64 / seqlen as f64;

    let pois = match Poisson::new(lambda) {
        Ok(p) => p,
        Err(_) => return Vec::new(),
    };

    let mut bad: Vec<u16> = counts
        .iter()
        .filter(|&(&_pos, &count)| {
            if count < min_occurrence {
                return false;
            }
            // R: ppois(count, lambda, lower.tail=FALSE) = P(X > count)
            //   = 1 - P(X <= count) = 1 - cdf(count)
            // Bonferroni: multiply by seqlen.
            let pval_corrected = (1.0 - pois.cdf(count as u64)) * seqlen as f64;
            pval_corrected < omega_b
        })
        .map(|(&pos, _)| pos)
        .collect();

    bad.sort_unstable();
    bad
}

/// Flag clusters whose birth substitutions fall predominantly on bad base positions.
///
/// Mirrors R's `isBadBaseFP(clust, birth_subs, minFraction, omegaB, minOccurence)`.
///
/// # Arguments
/// - `birth_subs_by_cluster`: for each cluster, a slice of its birth-substitution
///   positions (0-indexed)
/// - `bad_bases`: sorted vector of bad base positions from [`get_bad_bases`]
/// - `min_fraction`: minimum fraction of birth subs on bad bases to flag a cluster
///   (R default: 0.51)
///
/// # Returns
/// `Vec<bool>` of length `birth_subs_by_cluster.len()` — `true` when a cluster
/// is likely a false positive caused by bad bases.
pub fn is_bad_base_fp(
    birth_subs_by_cluster: &[&[u16]],
    bad_bases: &[u16],
    min_fraction: f64,
) -> Vec<bool> {
    let bad_set: HashSet<u16> = bad_bases.iter().copied().collect();

    birth_subs_by_cluster
        .iter()
        .map(|positions| {
            if positions.is_empty() {
                return false;
            }
            let bad_count = positions.iter().filter(|&&p| bad_set.contains(&p)).count();
            (bad_count as f64 / positions.len() as f64) >= min_fraction
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_uniform_trans(nq: usize, diag_count: u32, off_count: u32) -> Vec<u32> {
        let mut trans = vec![0u32; 16 * nq];
        for nti in 0..4 {
            for ntj in 0..4 {
                let count = if nti == ntj { diag_count } else { off_count };
                for q in 0..nq {
                    trans[(nti * 4 + ntj) * nq + q] = count;
                }
            }
        }
        trans
    }

    #[test]
    fn test_noqual_errfun_rows_sum_to_one() {
        let nq = 5;
        let trans = make_uniform_trans(nq, 100, 1);
        let err = noqual_errfun(&trans, nq, 1.0);
        assert_eq!(err.len(), 16 * nq);

        for nti in 0..4 {
            let row_sum: f64 = (0..4).map(|ntj| err[(nti * 4 + ntj) * nq]).sum();
            assert!(
                (row_sum - 1.0).abs() < 1e-6,
                "nti={nti} row sum = {row_sum}"
            );
        }
    }

    #[test]
    fn test_loess_errfun_shape() {
        let nq = 41;
        let qs: Vec<f64> = (0..nq).map(|i| i as f64).collect();
        let trans = make_uniform_trans(nq, 10_000, 10);
        let err = loess_errfun(&trans, &qs);
        assert_eq!(err.len(), 16 * nq);
        // All rates should be in [0, 1].
        for &r in &err {
            assert!(r >= 0.0 && r <= 1.0, "rate {r} out of [0,1]");
        }
    }

    #[test]
    fn test_inflate_err_increases_off_diagonal() {
        let nq = 5;
        let trans = make_uniform_trans(nq, 1000, 1);
        let qs: Vec<f64> = (0..nq).map(|i| i as f64).collect();
        let err = loess_errfun(&trans, &qs);
        let inflated = inflate_err(&err, nq, 2.0, false);

        for nti in 0..4 {
            for ntj in 0..4 {
                if nti != ntj {
                    for q in 0..nq {
                        assert!(
                            inflated[(nti * 4 + ntj) * nq + q] >= err[(nti * 4 + ntj) * nq + q],
                            "inflated rate should be >= original"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_accumulate_trans_sums_correctly() {
        let nq1 = 3;
        let nq2 = 5;
        let t1 = vec![1u32; 16 * nq1];
        let t2 = vec![2u32; 16 * nq2];
        let (acc, max_nq) = accumulate_trans(&[(&t1, nq1), (&t2, nq2)]);
        assert_eq!(max_nq, nq2);
        assert_eq!(acc.len(), 16 * nq2);
        // Columns 0..nq1 should have 1+2=3; columns nq1..nq2 should have 0+2=2.
        for r in 0..16 {
            for q in 0..nq1 {
                assert_eq!(acc[r * nq2 + q], 3, "r={r} q={q}");
            }
            for q in nq1..nq2 {
                assert_eq!(acc[r * nq2 + q], 2, "r={r} q={q}");
            }
        }
    }

    #[test]
    fn test_expand_err_matrix_diagonal() {
        let nq = 1;
        // All off-diagonals = 0.1 → diagonal should = 1 - 0.3 = 0.7
        let off = vec![0.1f64; 12 * nq];
        let err = expand_err_matrix(&off, nq);
        for nti in 0..4 {
            let diag_row = nti * 4 + nti;
            assert!((err[diag_row * nq] - 0.7).abs() < 1e-10);
        }
    }

    #[test]
    fn test_get_bad_bases_empty_input() {
        let bad = get_bad_bases(100, &[], 1e-10, 4);
        assert!(bad.is_empty());
    }
}
