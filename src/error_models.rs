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
use std::io;

use statrs::distribution::{DiscreteCDF, Poisson};

use crate::loess::{extrapolate_flat, loess_predict};
// The LOESS smoother, its configuration, and the default clamp bounds now live
// in [`crate::loess`]; re-export the public types so downstream consumers
// (`learn_errors`, the CLI) keep using `error_models::{...}` paths unchanged.
pub use crate::loess::{DEFAULT_MAX_ERROR_RATE, DEFAULT_MIN_ERROR_RATE, LoessConfig, LoessSurface};

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

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
/// Off-diagonal rates are clamped to `[config.min_error_rate, config.max_error_rate]`
/// (the `default` preset uses `[1e-7, 0.25]`; the `r-dada2` preset uses `[0, 1]`).
/// Diagonal entries are filled as `1 − sum(off-diagonals)` regardless of preset.
///
/// # Arguments
/// - `trans`: flat 16 × `nq` row-major matrix of transition counts
/// - `qual_scores`: quality score value for each column (length `nq`)
/// - `config`: surface + clamp bounds — see [`LoessConfig`]
///
/// # Returns
/// Flat 16 × `nq` row-major error rate matrix.
pub fn loess_errfun(trans: &[u32], qual_scores: &[f64], config: &LoessConfig) -> Vec<f64> {
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

            let raw_pred = loess_predict(qual_scores, &rlogp, &tot, 0.75, 2, config.surface);
            let filled = extrapolate_flat(raw_pred, nq);

            for q in 0..nq {
                let rate = if filled[q].is_finite() {
                    10.0f64.powf(filled[q])
                } else {
                    config.min_error_rate
                };
                off_diag[off_row * nq + q] = config.clamp(rate);
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
pub fn noqual_errfun(trans: &[u32], nq: usize, pseudocount: f64, config: &LoessConfig) -> Vec<f64> {
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
                config.clamp(obs[nti * 4 + ntj] / tot_init)
            } else {
                config.min_error_rate
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
    config: &LoessConfig,
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
                    config.clamp(filled[q])
                } else {
                    config.min_error_rate
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
pub fn pacbio_errfun(trans: &[u32], qual_scores: &[f64], config: &LoessConfig) -> Vec<f64> {
    let nq = qual_scores.len();
    assert_eq!(trans.len(), 16 * nq, "trans length must be 16 * nq");

    if nq > 0 && (qual_scores[nq - 1] - 93.0).abs() < 1e-9 {
        // Fit LOESS on columns 0..(nq-1).
        let sub_nq = nq - 1;
        let trans_sub: Vec<u32> = (0..16)
            .flat_map(|r| (0..sub_nq).map(move |q| trans[r * nq + q]))
            .collect();
        let qs_sub = &qual_scores[..sub_nq];
        let err_sub = loess_errfun(&trans_sub, qs_sub, config); // 16 × sub_nq

        // MLE for Q93: (count + 1) / (group_total + 4).
        let q93 = nq - 1;
        let mut err93 = [0.0f64; 16];
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
        loess_errfun(trans, qual_scores, config)
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
// External (user-supplied) error function
// ---------------------------------------------------------------------------

/// Row labels written to the trans/err TSVs exchanged with the external script.
///
/// Order matches the standard 16-row layout used elsewhere: outer = ref nt
/// (A,C,G,T), inner = query nt (A,C,G,T).  Same as R DADA2's
/// `paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))`.
const EXTERNAL_ROW_LABELS: [&str; 16] = [
    "A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A",
    "T2C", "T2G", "T2T",
];

/// RAII temp-dir guard used by [`external_errfun`].  Cleaned up on drop.
struct TempDir(std::path::PathBuf);

use std::sync::atomic::{AtomicU64, Ordering};
static TEMP_DIR_SEQ: AtomicU64 = AtomicU64::new(0);

impl TempDir {
    fn new(prefix: &str) -> io::Result<Self> {
        let seq = TEMP_DIR_SEQ.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("{prefix}-{}-{seq}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        Ok(Self(dir))
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Run a user-supplied command to fit error rates from a transition matrix.
///
/// `cmd` is whitespace-split into argv; two file paths — input trans TSV and
/// output err TSV — are appended as the final two arguments.  Both files use
/// R's `read.table(..., row.names=1, header=TRUE, check.names=FALSE)` layout:
/// the header line is a leading tab followed by `nq` integer column labels
/// (`0`, `1`, …, `nq-1`); each of the 16 data rows starts with a row label
/// (`A2A`, `A2C`, …, `T2T`) followed by `nq` whitespace-separated values.
///
/// The script must write the err matrix in the same shape (R's
/// `write.table(err, args[2], sep="\t", quote=FALSE, col.names=NA)` produces
/// it directly).
///
/// On success returns a flat row-major `Vec<f64>` of length `16 * nq`.
/// All values must be finite and in `[0, 1]`; row labels in the err TSV are
/// not strictly checked but their count must be exactly 16.
pub fn external_errfun(trans: &[u32], nq: usize, cmd: &str) -> Result<Vec<f64>, String> {
    if trans.len() != 16 * nq {
        return Err(format!(
            "external_errfun: trans length {} != 16 * nq ({})",
            trans.len(),
            16 * nq,
        ));
    }
    let argv: Vec<&str> = cmd.split_whitespace().collect();
    if argv.is_empty() {
        return Err("external_errfun: --errfun-cmd is empty".into());
    }

    let tmp = TempDir::new("dada2-rs-errfun")
        .map_err(|e| format!("external_errfun: failed to create temp dir: {e}"))?;
    let trans_path = tmp.0.join("trans.tsv");
    let err_path = tmp.0.join("err.tsv");

    write_trans_tsv(&trans_path, trans, nq)
        .map_err(|e| format!("external_errfun: write {}: {e}", trans_path.display()))?;

    let output = std::process::Command::new(argv[0])
        .args(&argv[1..])
        .arg(&trans_path)
        .arg(&err_path)
        .output()
        .map_err(|e| format!("external_errfun: failed to spawn '{cmd}': {e}"))?;

    if !output.status.success() {
        return Err(format!(
            "external_errfun: command '{cmd}' exited with {}\n--- stderr ---\n{}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim_end(),
        ));
    }

    let err_text = std::fs::read_to_string(&err_path)
        .map_err(|e| format!("external_errfun: read {}: {e}", err_path.display()))?;
    parse_err_tsv(&err_text, nq)
}

fn write_trans_tsv(path: &std::path::Path, trans: &[u32], nq: usize) -> io::Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    // Header: leading tab, then 0, 1, …, nq-1
    for q in 0..nq {
        write!(w, "\t{q}")?;
    }
    writeln!(w)?;
    for r in 0..16 {
        write!(w, "{}", EXTERNAL_ROW_LABELS[r])?;
        for q in 0..nq {
            write!(w, "\t{}", trans[r * nq + q])?;
        }
        writeln!(w)?;
    }
    Ok(())
}

fn parse_err_tsv(text: &str, nq: usize) -> Result<Vec<f64>, String> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| "external_errfun: err TSV is empty".to_string())?;
    let header_cols: Vec<&str> = header.split('\t').collect();
    if header_cols.len() != nq + 1 {
        return Err(format!(
            "external_errfun: err TSV header has {} fields, expected {} (1 row-name + {} quality scores)",
            header_cols.len(),
            nq + 1,
            nq,
        ));
    }

    let mut err = vec![0.0f64; 16 * nq];
    let mut row = 0usize;
    for line in lines {
        if line.is_empty() {
            continue;
        }
        if row >= 16 {
            return Err("external_errfun: err TSV has more than 16 data rows".into());
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != nq + 1 {
            return Err(format!(
                "external_errfun: err TSV row {} has {} fields, expected {}",
                row + 1,
                cols.len(),
                nq + 1,
            ));
        }
        for q in 0..nq {
            let s = cols[q + 1];
            let v: f64 = s.parse().map_err(|e| {
                format!(
                    "external_errfun: err TSV row {} column {} ({:?}): {e}",
                    row + 1,
                    q,
                    s,
                )
            })?;
            if !v.is_finite() || !(0.0..=1.0).contains(&v) {
                return Err(format!(
                    "external_errfun: err[{}][{}] = {v} is not a valid probability in [0, 1]",
                    EXTERNAL_ROW_LABELS[row], q,
                ));
            }
            err[row * nq + q] = v;
        }
        row += 1;
    }
    if row != 16 {
        return Err(format!(
            "external_errfun: err TSV has {row} data rows, expected 16",
        ));
    }
    Ok(err)
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
        let err = noqual_errfun(&trans, nq, 1.0, &LoessConfig::default());
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
    fn loess_config_default_and_r_dada2_presets() {
        let d = LoessConfig::default();
        assert!(matches!(d.surface, LoessSurface::Direct));
        assert_eq!(d.max_error_rate, DEFAULT_MAX_ERROR_RATE);
        assert_eq!(d.min_error_rate, DEFAULT_MIN_ERROR_RATE);

        let r = LoessConfig::r_dada2();
        assert!(matches!(r.surface, LoessSurface::Interpolate { cell }
                         if (cell - 0.2).abs() < 1e-12));
        // Both presets clamp to [1e-7, 0.25] — R DADA2's loessErrfun does too
        // (errorModels.R:53-56).  The presets differ only in the surface.
        assert_eq!(r.max_error_rate, DEFAULT_MAX_ERROR_RATE);
        assert_eq!(r.min_error_rate, DEFAULT_MIN_ERROR_RATE);

        for cfg in [&d, &r] {
            assert_eq!(cfg.clamp(0.5), 0.25); // upper-clamps
            assert_eq!(cfg.clamp(1e-10), 1e-7); // lower-clamps
            assert_eq!(cfg.clamp(0.1), 0.1); // in-range pass-through
        }
    }

    #[test]
    fn test_loess_errfun_shape() {
        let nq = 41;
        let qs: Vec<f64> = (0..nq).map(|i| i as f64).collect();
        let trans = make_uniform_trans(nq, 10_000, 10);
        let err = loess_errfun(&trans, &qs, &LoessConfig::default());
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
        let err = loess_errfun(&trans, &qs, &LoessConfig::default());
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

    /// Exercises `external_errfun` via a tiny POSIX-shell script: read trans,
    /// emit a constant 0.01 error rate everywhere. Verifies the round-trip
    /// (write trans TSV, spawn process, read back err TSV, validate shape).
    /// Skipped on non-unix targets where `/bin/sh` semantics may differ.
    #[cfg(unix)]
    #[test]
    fn test_external_errfun_roundtrip() {
        let nq = 4;
        let trans: Vec<u32> = (0..16 * nq).map(|i| i as u32).collect();

        let dir = TempDir::new("dada2-rs-extfun-test").unwrap();
        let script = dir.0.join("constant_err.sh");
        // Awk script: print header from trans, then 16 rows of constant 0.01
        let body = "#!/bin/sh\n\
                    awk -F'\\t' 'BEGIN{OFS=\"\\t\"}\n\
                    NR==1{print; nq=NF-1; next}\n\
                    {printf \"%s\", $1; for(i=1;i<=nq;i++) printf \"\\t0.01\"; print \"\"}' \"$1\" > \"$2\"\n";
        std::fs::write(&script, body).unwrap();
        // Make executable
        use std::os::unix::fs::PermissionsExt;
        let mut perms = std::fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&script, perms).unwrap();

        // Invoke via `/bin/sh <script>` rather than exec'ing the script file
        // directly: tests run in parallel, and exec'ing a just-written file can
        // race with another test's fork() and fail with ETXTBSY ("Text file
        // busy") on Linux. Running it through sh only *reads* the script.
        let cmd = format!("/bin/sh {}", script.to_string_lossy());
        let err = external_errfun(&trans, nq, &cmd).expect("external errfun should succeed");
        assert_eq!(err.len(), 16 * nq);
        for v in &err {
            assert!((v - 0.01).abs() < 1e-12, "value {v} != 0.01");
        }
        // Hold dir alive until after the spawn so the script isn't cleaned up early.
        drop(dir);
    }

    /// A script that exits with non-zero status should surface a useful error.
    #[cfg(unix)]
    #[test]
    fn test_external_errfun_propagates_nonzero_exit() {
        let dir = TempDir::new("dada2-rs-extfun-test").unwrap();
        let script = dir.0.join("fail.sh");
        std::fs::write(&script, "#!/bin/sh\necho 'oops' >&2\nexit 7\n").unwrap();
        use std::os::unix::fs::PermissionsExt;
        let mut perms = std::fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&script, perms).unwrap();

        let trans = vec![0u32; 16 * 4];
        // Run via `/bin/sh <script>` to avoid the parallel-test ETXTBSY race
        // (see test_external_errfun_roundtrip).
        let result = external_errfun(&trans, 4, &format!("/bin/sh {}", script.to_string_lossy()));
        let err_msg = result.expect_err("script exits 7 — should fail");
        assert!(err_msg.contains("exited with"), "got: {err_msg}");
        assert!(
            err_msg.contains("oops"),
            "stderr should be propagated, got: {err_msg}"
        );
        drop(dir);
    }

    /// A script returning a value outside [0,1] should be rejected.
    #[cfg(unix)]
    #[test]
    fn test_external_errfun_rejects_invalid_probability() {
        let dir = TempDir::new("dada2-rs-extfun-test").unwrap();
        let script = dir.0.join("bad_prob.sh");
        // Header line + 16 rows where the first cell is 1.5
        std::fs::write(
            &script,
            "#!/bin/sh\nawk -F'\\t' 'BEGIN{OFS=\"\\t\"}\n\
             NR==1{print; nq=NF-1; next}\n\
             {printf \"%s\\t1.5\", $1; for(i=2;i<=nq;i++) printf \"\\t0.01\"; print \"\"}' \"$1\" > \"$2\"\n",
        )
        .unwrap();
        use std::os::unix::fs::PermissionsExt;
        let mut perms = std::fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&script, perms).unwrap();

        let trans = vec![0u32; 16 * 4];
        // Run via `/bin/sh <script>` to avoid the parallel-test ETXTBSY race
        // (see test_external_errfun_roundtrip).
        let result = external_errfun(&trans, 4, &format!("/bin/sh {}", script.to_string_lossy()));
        let err_msg = result.expect_err("1.5 is not a valid probability");
        assert!(
            err_msg.contains("not a valid probability"),
            "got: {err_msg}"
        );
        drop(dir);
    }
}
