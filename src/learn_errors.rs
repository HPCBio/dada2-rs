//! Error model learning — iterative self-consistency estimation.
//!
//! Mirrors R's `learnErrors()`: given a set of dereplicated samples (supplied
//! as the JSON output from the `subsample` subcommand), this module runs the
//! DADA2 algorithm iteratively, accumulating transition counts and re-fitting
//! the error model until convergence.
//!
//! ## Output
//! [`LearnErrorsResult`] contains three flat row-major matrices (16 × `nq`):
//! - `trans`   — accumulated transition *counts* from the final iteration.
//! - `err_in`  — error *rates* fed into the last DADA run.
//! - `err_out` — error *rates* estimated from `trans` via the chosen error function.

use std::io;
use std::path::PathBuf;

use serde::Deserialize;

use crate::containers::Raw;
use crate::dada::{DadaParams, RawInput, dada_uniques};
use crate::error_models::{
    accumulate_trans, binned_qual_errfun, loess_errfun, noqual_errfun, pacbio_errfun,
};
use crate::misc::nt_encode;
use crate::nwalign::{raw_align, AlignParams};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const MIN_ERR: f64 = 1e-7;
const MAX_ERR: f64 = 0.25;
/// Maximum difference between consecutive error-model iterations considered converged.
const CONVERGENCE_TOL: f64 = 1e-7;

// ---------------------------------------------------------------------------
// Error function selection
// ---------------------------------------------------------------------------

/// Which error model fitting function to apply to the accumulated transition matrix.
#[derive(Clone, Debug)]
pub enum ErrFun {
    /// Locally-weighted polynomial regression (default for Illumina).
    Loess,
    /// Quality-score-free: one rate per transition type, broadcast across all Q.
    Noqual {
        pseudocount: f64,
    },
    /// Piecewise linear interpolation between anchor quality bins.
    BinnedQual {
        bins: Vec<f64>,
    },
    /// PacBio-specific model.
    PacBio,
}

impl ErrFun {
    /// Apply the error function to a transition count matrix and return error rates (16 × nq).
    pub fn apply(&self, trans: &[u32], nq: usize) -> Result<Vec<f64>, String> {
        let qual_scores: Vec<f64> = (0..nq).map(|q| q as f64).collect();
        match self {
            ErrFun::Loess => Ok(loess_errfun(trans, &qual_scores)),
            ErrFun::Noqual { pseudocount } => {
                Ok(noqual_errfun(trans, nq, *pseudocount))
            }
            ErrFun::BinnedQual { bins } => {
                binned_qual_errfun(trans, &qual_scores, bins)
            }
            ErrFun::PacBio => Ok(pacbio_errfun(trans, &qual_scores)),
        }
    }
}

// ---------------------------------------------------------------------------
// JSON input types (mirror the subsample / derep output)
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct UniqueEntry {
    sequence: String,
    count: u64,
    mean_quality: Vec<f64>,
}

#[derive(Deserialize)]
struct DerepJson {
    uniques: Vec<UniqueEntry>,
}

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

/// Output of [`learn_errors`].
pub struct LearnErrorsResult {
    /// Flat row-major 16 × `nq` accumulated transition count matrix.
    pub trans: Vec<u32>,
    /// Flat row-major 16 × `nq` error-rate matrix used as input to the last DADA run.
    pub err_in: Vec<f64>,
    /// Flat row-major 16 × `nq` error-rate matrix estimated from `trans`.
    pub err_out: Vec<f64>,
    /// Number of quality-score columns.
    pub nq: usize,
    /// Whether the model converged within `max_consist` iterations.
    pub converged: bool,
    /// Number of iterations completed.
    pub iterations: usize,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Load and parse a single JSON (or gzip-compressed JSON) derep file into a `Vec<RawInput>`.
fn load_json(path: &PathBuf) -> io::Result<Vec<RawInput>> {
    let derep: DerepJson = crate::misc::read_json_file(path)?;

    let inputs = derep
        .uniques
        .into_iter()
        .map(|u| RawInput {
            seq: u.sequence,
            abundance: u.count as u32,
            prior: false,
            quals: Some(u.mean_quality),
        })
        .collect();

    Ok(inputs)
}

/// Determine `nq` (number of quality columns) from the maximum rounded quality
/// score seen across all input sets.  Falls back to 41 when no quality data is
/// present.
fn detect_nq(all_inputs: &[Vec<RawInput>]) -> usize {
    let max_q = all_inputs
        .iter()
        .flat_map(|s| s.iter())
        .filter_map(|r| r.quals.as_deref())
        .flat_map(|qs| qs.iter().copied())
        .map(|q| q.round() as usize)
        .max()
        .unwrap_or(40);
    max_q + 1
}

/// Build a Phred-score-based initial error model (16 × `nq`).
///
/// Off-diagonal rates: 10^(-q/10) / 3 (split evenly across three substitution types).
/// Diagonal rates: 1 - 10^(-q/10).
/// All rates are clamped to `[MIN_ERR, MAX_ERR]`.
fn init_err_model(nq: usize) -> Vec<f64> {
    let mut err = vec![0.0f64; 16 * nq];
    for q in 0..nq {
        let prob_error = 10.0f64.powf(-(q as f64) / 10.0);
        let prob_mismatch = (prob_error / 3.0).clamp(MIN_ERR, MAX_ERR);
        let prob_match = (1.0 - prob_error).clamp(MIN_ERR, MAX_ERR);

        for nti in 0..4usize {
            for ntj in 0..4usize {
                let row = nti * 4 + ntj;
                err[row * nq + q] = if nti == ntj { prob_match } else { prob_mismatch };
            }
        }
    }
    err
}

/// Build the 16 × `nq` transition count matrix from a single DADA result.
///
/// For every input raw that was assigned to a cluster, aligns the raw against
/// its cluster center and counts (ref_nt → query_nt) transitions at each
/// quality-score column.
fn build_trans_mat(
    inputs: &[RawInput],
    result: &crate::dada::DadaResult,
    align_params: &AlignParams,
    nq: usize,
) -> Vec<u32> {
    let mut trans = vec![0u32; 16 * nq];

    // Build Raw stubs for cluster centers (no quality scores needed).
    let center_raws: Vec<Raw> = result
        .clusters
        .iter()
        .map(|c| Raw::new(c.sequence.clone(), None, 0, false))
        .collect();

    for (i, inp) in inputs.iter().enumerate() {
        let ci = match result.map[i] {
            Some(ci) => ci,
            None => continue,
        };

        let quals = match &inp.quals {
            Some(q) => q.as_slice(),
            None => continue, // no quality data — cannot build trans
        };

        let seq: Vec<u8> = inp.seq.bytes().map(nt_encode).collect();
        let raw_query = Raw::new(seq, Some(quals), inp.abundance, false);
        let raw_center = &center_raws[ci];

        // Align center (ref = al[0]) against the raw (query = al[1]).
        let al = match raw_align(raw_center, &raw_query, align_params) {
            Some(a) => a,
            None => continue,
        };

        let al_ref = &al[0];
        let al_qry = &al[1];
        let reads = inp.abundance;
        let mut qpos = 0usize; // position in original query sequence

        for alpos in 0..al_ref.len() {
            let nt0 = al_ref[alpos]; // reference nucleotide (encoded)
            let nt1 = al_qry[alpos]; // query nucleotide (encoded)

            // Advance query position for any non-gap query character.
            let qry_is_nt = matches!(nt1, 1..=4);

            if matches!(nt0, 1..=4) && qry_is_nt {
                // Both reference and query are ACGT — accumulate.
                let q = (quals[qpos].round() as usize).min(nq - 1);
                let row = (nt0 as usize - 1) * 4 + (nt1 as usize - 1);
                trans[row * nq + q] = trans[row * nq + q].saturating_add(reads);
            }

            // Advance query position for any non-gap query character (incl. N).
            if nt1 != b'-' {
                qpos += 1;
            }
        }
    }

    trans
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Learn an error model from pre-subsampled derep JSON files.
///
/// # Arguments
/// - `json_files`    — paths to JSON files produced by the `subsample` subcommand
/// - `errfun`        — which error model fitting function to use
/// - `dada_params`   — DADA2 algorithm parameters (error matrix is overwritten each iteration)
/// - `align_params`  — NW alignment parameters
/// - `max_consist`   — maximum self-consistency iterations (R default: 10)
/// - `verbose`       — print progress to stderr
///
/// # Returns
/// [`LearnErrorsResult`] with transition counts and error-rate matrices, or an I/O error.
pub fn learn_errors(
    json_files: &[PathBuf],
    errfun: &ErrFun,
    mut dada_params: DadaParams,
    align_params: &AlignParams,
    max_consist: usize,
    verbose: bool,
) -> io::Result<LearnErrorsResult> {
    // ---- Load all samples ----
    if json_files.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "No JSON input files provided"));
    }

    let all_inputs: Vec<Vec<RawInput>> = json_files
        .iter()
        .map(load_json)
        .collect::<io::Result<_>>()?;

    let nq = detect_nq(&all_inputs);

    if verbose {
        eprintln!(
            "[learn_errors] {} sample(s), nq={}, max_consist={}",
            all_inputs.len(),
            nq,
            max_consist,
        );
    }

    // ---- Initialise error model ----
    let mut err = init_err_model(nq);
    dada_params.err_ncol = nq;

    let mut converged = false;
    let mut iterations = 0usize;
    let mut trans_final = vec![0u32; 16 * nq];
    let mut err_in_final = err.clone();
    let mut err_out_final = err.clone();

    for iter in 0..max_consist {
        iterations = iter + 1;

        // ---- Run DADA on each sample ----
        let mut sample_trans_pairs: Vec<(Vec<u32>, usize)> = Vec::new();

        for (si, inputs) in all_inputs.iter().enumerate() {
            dada_params.err_mat = err.clone();

            match dada_uniques(inputs, &dada_params) {
                Ok(result) => {
                    let t = build_trans_mat(inputs, &result, align_params, nq);
                    sample_trans_pairs.push((t, nq));

                    if verbose {
                        let nclusters = result.clusters.len();
                        eprintln!(
                            "[learn_errors] iter={} sample={}: {} cluster(s)",
                            iter + 1,
                            si + 1,
                            nclusters,
                        );
                    }
                }
                Err(e) => {
                    if verbose {
                        eprintln!(
                            "[learn_errors] iter={} sample={}: dada_uniques failed: {}",
                            iter + 1,
                            si + 1,
                            e,
                        );
                    }
                }
            }
        }

        if sample_trans_pairs.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "All DADA runs failed; cannot estimate error model",
            ));
        }

        // ---- Accumulate transition matrices ----
        let refs: Vec<(&[u32], usize)> = sample_trans_pairs
            .iter()
            .map(|(t, nq)| (t.as_slice(), *nq))
            .collect();
        let (acc_trans, _) = accumulate_trans(&refs);

        // ---- Estimate new error rates ----
        let new_err = errfun.apply(&acc_trans, nq).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("errfun failed: {e}"))
        })?;

        // ---- Check convergence ----
        let max_delta = err
            .iter()
            .zip(new_err.iter())
            .map(|(&a, &b)| (a - b).abs())
            .fold(0.0f64, f64::max);

        if verbose {
            eprintln!(
                "[learn_errors] iter={}: max |err_in - err_out| = {:.2e}",
                iter + 1,
                max_delta,
            );
        }

        err_in_final = err.clone();
        trans_final = acc_trans;
        err_out_final = new_err.clone();

        if max_delta < CONVERGENCE_TOL {
            converged = true;
            if verbose {
                eprintln!("[learn_errors] converged after {} iteration(s)", iter + 1);
            }
            break;
        }

        err = new_err;
    }

    if !converged && verbose {
        eprintln!(
            "[learn_errors] did not converge within {} iteration(s)",
            max_consist
        );
    }

    Ok(LearnErrorsResult {
        trans: trans_final,
        err_in: err_in_final,
        err_out: err_out_final,
        nq,
        converged,
        iterations,
    })
}
