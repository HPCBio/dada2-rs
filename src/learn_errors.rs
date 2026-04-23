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

use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::containers::Raw;
use crate::dada::{DadaParams, RawInput, dada_uniques};
use crate::derep::dereplicate;
use crate::error_models::{
    accumulate_trans, binned_qual_errfun, loess_errfun, noqual_errfun, pacbio_errfun,
};
use crate::misc::nt_encode;
use crate::nwalign::{raw_align_with_buf, AlignBuffers, AlignParams};

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
// FASTQ loading with subsampling
// ---------------------------------------------------------------------------

fn is_gz(path: &PathBuf) -> bool {
    path.extension().and_then(|e| e.to_str()) == Some("gz")
}

/// Dereplicate FASTQ files and return one `Vec<RawInput>` per sample.
///
/// Files are processed in order (or shuffled when `randomize` is true) and
/// processing stops once the cumulative base count reaches `nbases`.
pub fn load_fastq_samples(
    paths: &[PathBuf],
    nbases: u64,
    randomize: bool,
    seed: Option<u64>,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    verbose: bool,
) -> io::Result<Vec<Vec<RawInput>>> {
    let mut ordered: Vec<&PathBuf> = paths.iter().collect();
    if randomize {
        use rand::seq::SliceRandom as _;
        use rand::SeedableRng as _;
        if let Some(s) = seed {
            ordered.shuffle(&mut rand::rngs::SmallRng::seed_from_u64(s));
        } else {
            ordered.shuffle(&mut rand::thread_rng());
        }
    }

    let mut all_inputs = Vec::new();
    let mut total_bases: u64 = 0;

    for path in ordered {
        let derep = if is_gz(path) {
            dereplicate(MultiGzDecoder::new(File::open(path)?), phred_offset, pool, verbose)?
        } else {
            dereplicate(File::open(path)?, phred_offset, pool, verbose)?
        };

        let file_bases: u64 = derep
            .uniques
            .iter()
            .map(|(seq, count)| seq.len() as u64 * count)
            .sum();

        let inputs: Vec<RawInput> = derep
            .uniques
            .into_iter()
            .enumerate()
            .map(|(i, (seq, count))| RawInput {
                seq: String::from_utf8(seq)
                    .unwrap_or_default(),
                abundance: count as u32,
                prior: false,
                quals: Some(derep.quals[i].clone()),
            })
            .collect();

        if verbose {
            eprintln!(
                "[learn-errors] loaded {} unique(s) from {} ({} bases)",
                inputs.len(),
                path.display(),
                file_bases,
            );
        }

        all_inputs.push(inputs);
        total_bases += file_bases;

        if total_bases >= nbases {
            if verbose {
                eprintln!(
                    "[learn-errors] reached {} bases after {} file(s); stopping subsampling",
                    total_bases,
                    all_inputs.len(),
                );
            }
            break;
        }
    }

    Ok(all_inputs)
}

// ---------------------------------------------------------------------------
// JSON-backed sample loading (for errors-from-sample)
// ---------------------------------------------------------------------------

/// Wire format for one unique sequence in a derep/sample JSON file.
#[derive(Deserialize)]
struct UniqueEntryJson {
    sequence: String,
    count: u64,
    mean_quality: Vec<f64>,
}

/// Top-level structure of a derep/sample JSON file.
#[derive(Deserialize)]
struct SampleJson {
    uniques: Vec<UniqueEntryJson>,
}

/// Load pre-computed derep JSON files (from `sample` or `derep`) as samples
/// for error learning.
///
/// Each file is expected to have the same structure written by the `derep` and
/// `sample` subcommands: a top-level object with an `"uniques"` array of
/// `{sequence, count, mean_quality}` entries.
pub fn load_derep_samples(paths: &[PathBuf]) -> io::Result<Vec<Vec<RawInput>>> {
    paths
        .iter()
        .map(|path| {
            let file = File::open(path)?;
            let sample: SampleJson = serde_json::from_reader(file).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("failed to parse {}: {e}", path.display()),
                )
            })?;

            let inputs = sample
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
        })
        .collect()
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

    let mut buf = AlignBuffers::new();
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
        let al = match raw_align_with_buf(raw_center, &raw_query, align_params, &mut buf) {
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

// ---------------------------------------------------------------------------
// Per-iteration diagnostics
// ---------------------------------------------------------------------------

/// Cluster-count summary for one sample in one iteration.
#[derive(Serialize)]
struct SampleIterDiag {
    sample: usize,
    n_clusters: usize,
    total_reads: u32,
    n_initial: usize,
    n_abundance: usize,
    n_prior: usize,
    n_singleton: usize,
    nalign: u32,
    nshroud: u32,
}

/// Written to `<diag_dir>/iter_NNN.json` after each full iteration.
#[derive(Serialize)]
struct IterDiag {
    iter: usize,
    converged: bool,
    max_delta: f64,
    samples: Vec<SampleIterDiag>,
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Learn an error model from pre-loaded dereplicated samples.
///
/// # Arguments
/// - `all_inputs`    — one `Vec<RawInput>` per sample (from `load_fastq_samples`)
/// - `errfun`        — which error model fitting function to use
/// - `dada_params`   — DADA2 algorithm parameters (error matrix is overwritten each iteration)
/// - `align_params`  — NW alignment parameters
/// - `max_consist`   — maximum self-consistency iterations (R default: 10)
/// - `verbose`       — print progress to stderr
/// - `diag_dir`      — if `Some`, write per-iteration cluster diagnostics as
///                     `iter_001.json`, `iter_002.json`, … into this directory
///
/// # Returns
/// [`LearnErrorsResult`] with transition counts and error-rate matrices, or an I/O error.
pub fn learn_errors(
    all_inputs: Vec<Vec<RawInput>>,
    errfun: &ErrFun,
    mut dada_params: DadaParams,
    align_params: &AlignParams,
    max_consist: usize,
    verbose: bool,
    diag_dir: Option<&Path>,
) -> io::Result<LearnErrorsResult> {
    if all_inputs.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "No input samples provided"));
    }

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
        // Set the error matrix once per iteration (shared across all samples).
        dada_params.err_mat = err.clone();

        // Parallel pass: each sample runs dada_uniques independently.
        // Rayon's work-stealing distributes threads across samples; if only
        // one sample exists all threads serve its inner b_compare_parallel.
        let collect_diags = diag_dir.is_some();
        let sample_results: Vec<(usize, Result<(Vec<u32>, crate::dada::DadaResult), String>)> =
            all_inputs
                .par_iter()
                .enumerate()
                .map(|(si, inputs)| {
                    let outcome = dada_uniques(inputs, &dada_params).map(|result| {
                        let t = build_trans_mat(inputs, &result, align_params, nq);
                        (t, result)
                    });
                    (si, outcome)
                })
                .collect();

        // Serial post-processing: logging and diagnostics.
        let mut sample_trans_pairs: Vec<(Vec<u32>, usize)> = Vec::new();
        let mut sample_diags: Vec<SampleIterDiag> = Vec::new();

        for (si, outcome) in sample_results {
            match outcome {
                Ok((t, result)) => {
                    sample_trans_pairs.push((t, nq));

                    if verbose {
                        eprintln!(
                            "[learn_errors] iter={} sample={}: {} cluster(s)",
                            iter + 1,
                            si + 1,
                            result.clusters.len(),
                        );
                    }

                    if collect_diags {
                        use crate::containers::BirthType;
                        let total_reads = result.clusters.iter().map(|c| c.reads).sum();
                        let mut diag = SampleIterDiag {
                            sample: si + 1,
                            n_clusters: result.clusters.len(),
                            total_reads,
                            n_initial:   0,
                            n_abundance: 0,
                            n_prior:     0,
                            n_singleton: 0,
                            nalign:  result.nalign,
                            nshroud: result.nshroud,
                        };
                        for c in &result.clusters {
                            match c.birth_type {
                                BirthType::Initial   => diag.n_initial   += 1,
                                BirthType::Abundance => diag.n_abundance += 1,
                                BirthType::Prior     => diag.n_prior     += 1,
                                BirthType::Singleton => diag.n_singleton += 1,
                            }
                        }
                        sample_diags.push(diag);
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

        let iter_converged = max_delta < CONVERGENCE_TOL;

        if verbose {
            eprintln!(
                "[learn_errors] iter={}: max |err_in - err_out| = {:.2e}",
                iter + 1,
                max_delta,
            );
        }

        // ---- Write per-iteration diagnostics ----
        if let Some(dir) = diag_dir {
            let diag = IterDiag {
                iter: iter + 1,
                converged: iter_converged,
                max_delta,
                samples: sample_diags,
            };
            let path = dir.join(format!("iter_{:03}.json", iter + 1));
            let json = serde_json::to_string_pretty(&diag)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            std::fs::write(&path, json)?;
            if verbose {
                eprintln!("[learn_errors] diagnostics written to {}", path.display());
            }
        }

        err_in_final = err.clone();
        trans_final = acc_trans;
        err_out_final = new_err.clone();

        if iter_converged {
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
