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

use crate::cluster_trace::TraceParams;
use crate::containers::Raw;
use crate::dada::{DadaParams, RawInput, dada_uniques_cached};
use crate::derep::dereplicate;
use crate::error_models::{
    accumulate_trans, binned_qual_errfun, external_errfun, loess_errfun, noqual_errfun,
    pacbio_errfun,
};
use crate::misc::nt_encode;
use crate::nwalign::{AlignBuffers, AlignParams, raw_align_with_buf};

// ---------------------------------------------------------------------------
// Public parameter types
// ---------------------------------------------------------------------------

/// Diagnostic and trace output options for [`learn_errors`].
pub struct LearnDiagOptions<'a> {
    pub verbose: bool,
    pub diag_dir: Option<&'a Path>,
    pub cluster_trace_dir: Option<&'a Path>,
    pub trace_params: TraceParams,
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

// TODO: check why these aren't used
#[allow(dead_code)]
const MIN_ERR: f64 = 1e-7;
#[allow(dead_code)]
const MAX_ERR: f64 = 0.25;
/// Maximum difference between consecutive error-model iterations considered converged.
const CONVERGENCE_TOL: f64 = 1e-7;
/// Bail with a "stalled" stop reason if `max_delta` fails to set a new
/// best-so-far for this many consecutive iterations. Catches the oscillation
/// pattern where err_out cycles between a few values without settling — most
/// commonly seen with too-loose k-mer screening (e.g. k=4 on Illumina).
const STALL_PATIENCE: usize = 3;
/// Relative improvement required to count as a new best-so-far. Without
/// this floor, tiny floating-point drifts would forever look like progress.
const STALL_REL_IMPROVE: f64 = 0.99;

// ---------------------------------------------------------------------------
// Error function selection
// ---------------------------------------------------------------------------

/// Which error model fitting function to apply to the accumulated transition matrix.
#[derive(Clone, Debug)]
pub enum ErrFun {
    /// Locally-weighted polynomial regression (default for Illumina).
    Loess,
    /// Quality-score-free: one rate per transition type, broadcast across all Q.
    Noqual { pseudocount: f64 },
    /// Piecewise linear interpolation between anchor quality bins.
    BinnedQual { bins: Vec<f64> },
    /// PacBio-specific model.
    PacBio,
    /// User-supplied command that reads a trans TSV and writes an err TSV.
    /// See [`external_errfun`] for the wire format.
    External { command: String },
}

impl ErrFun {
    /// Apply the error function to a transition count matrix and return error rates (16 × nq).
    pub fn apply(&self, trans: &[u32], nq: usize) -> Result<Vec<f64>, String> {
        let qual_scores: Vec<f64> = (0..nq).map(|q| q as f64).collect();
        match self {
            ErrFun::Loess => Ok(loess_errfun(trans, &qual_scores)),
            ErrFun::Noqual { pseudocount } => Ok(noqual_errfun(trans, nq, *pseudocount)),
            ErrFun::BinnedQual { bins } => binned_qual_errfun(trans, &qual_scores, bins),
            ErrFun::PacBio => Ok(pacbio_errfun(trans, &qual_scores)),
            ErrFun::External { command } => external_errfun(trans, nq, command),
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
    /// Why the iteration loop stopped.
    pub stop_reason: StopReason,
}

/// Parameters captured from a `learn-errors` / `errors-from-sample` run, embedded
/// in the JSON output for provenance and downstream consistency checks.
///
/// All fields here either feed `dada_uniques` directly or determine how the
/// transition counts that produced this err model were generated. Embedding them
/// lets a downstream `dada` invocation either (a) warn when its own CLI flags
/// disagree with what learned the model, or (b) opt-in inherit them so the same
/// parameter vector is used end-to-end.
///
/// `Deserialize` is permissive (`#[serde(default)]` everywhere) so error-model
/// JSON files written before this field existed still load.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct LearnedErrParams {
    /// Error-estimation function name as a stable string token: `"loess"`,
    /// `"noqual"`, `"binned-qual"`, or `"pacbio"`.
    #[serde(default)]
    pub errfun: String,
    /// Pseudocount for `noqual` errfun (None for other variants).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub errfun_pseudocount: Option<f64>,
    /// Quality-score bin edges for `binned-qual` errfun (None for others).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub errfun_bins: Option<Vec<f64>>,
    /// User-supplied command for `external` errfun (None for others).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub errfun_cmd: Option<String>,
    /// Maximum self-consistency iterations the model was allowed to run.
    #[serde(default)]
    pub max_consist: usize,

    // ---- DadaParams: significance / fold thresholds ----
    #[serde(default)]
    pub omega_a: f64,
    #[serde(default)]
    pub omega_c: f64,
    #[serde(default)]
    pub omega_p: f64,
    #[serde(default)]
    pub min_fold: f64,
    #[serde(default)]
    pub min_hamming: u32,
    #[serde(default)]
    pub min_abund: u32,
    #[serde(default)]
    pub detect_singletons: bool,
    #[serde(default)]
    pub use_quals: bool,
    #[serde(default)]
    pub greedy: bool,

    // ---- AlignParams ----
    #[serde(default)]
    pub match_score: i32,
    #[serde(default)]
    pub mismatch: i32,
    #[serde(default)]
    pub gap_p: i32,
    #[serde(default)]
    pub homo_gap_p: i32,
    #[serde(default)]
    pub use_kmers: bool,
    #[serde(default)]
    pub kdist_cutoff: f64,
    #[serde(default)]
    pub kmer_size: usize,
    #[serde(default)]
    pub band: i32,
    #[serde(default)]
    pub vectorized: bool,
    #[serde(default)]
    pub gapless: bool,
}

/// Outcome of the self-consistency iteration in `learn_errors`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StopReason {
    /// `max_delta < CONVERGENCE_TOL`: the model is stable.
    Converged,
    /// `max_delta` stopped improving for `STALL_PATIENCE` iterations.
    /// Typically signals oscillation (e.g. too-loose k-mer screening).
    Stalled,
    /// Reached `max_consist` without converging or stalling.
    MaxConsistReached,
}

// ---------------------------------------------------------------------------
// FASTQ loading with subsampling
// ---------------------------------------------------------------------------

fn is_gz(path: &Path) -> bool {
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
        use rand::SeedableRng as _;
        use rand::seq::SliceRandom as _;
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
            dereplicate(
                MultiGzDecoder::new(File::open(path)?),
                phred_offset,
                pool,
                verbose,
            )?
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
                seq: String::from_utf8(seq).unwrap_or_default(),
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
    #[serde(default)]
    sort_order: Option<String>,
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
            let sample: SampleJson = crate::misc::read_tagged_json(path, &["derep", "sample"])
                .map_err(|e| {
                    io::Error::new(e.kind(), format!("failed to parse {}: {e}", path.display()))
                })?;

            let already_sorted = sample.sort_order.as_deref() == Some("abundance_desc");
            let mut inputs: Vec<RawInput> = sample
                .uniques
                .into_iter()
                .map(|u| RawInput {
                    seq: u.sequence,
                    abundance: u.count as u32,
                    prior: false,
                    quals: Some(u.mean_quality),
                })
                .collect();
            // Defensive sort: derep JSONs produced by older versions (or by
            // other tools) may not be abundance-sorted. The DADA2 algorithm
            // assumes the most-abundant raw is at index 0 (issue #4).
            // `sort_by` is stable so ties preserve the file's original order.
            if !already_sorted {
                inputs.sort_by(|a, b| b.abundance.cmp(&a.abundance));
            }
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

/// Run R's init pass and return the err_mat to feed into iter-1.
///
/// Matches R `dada.R` `nconsist=0` step:
///   1. `erri <- matrix(1, nrow=16, ncol=...)` — every entry 1.0.
///   2. Run `dada_uniques` with `MAX_CLUST=1` on every sample (no buds; just
///      aligns each raw against cluster 0 to populate transitions).
///   3. Accumulate trans across samples and fit `errorEstimationFunction`.
///   4. Force diagonal rows (A2A, C2C, G2G, T2T) to 1.0 at every q.
///
/// Issue #4 traced the iter-1 over-budding to a different starting err_mat
/// shape than R's; the c3f4e88 shortcut (match=1.0, mismatch=MAX_ERR/3 at
/// every q) closed most of the trans-count gap but left R's per-Q
/// calibration on the table. Running the actual init pass produces a
/// Q-resolved off-diagonal shape derived from real data instead.
///
/// Side effect: populates `raw_cache` so iter 1's encoding + k-mer build is
/// reused.
fn run_init_pass(
    all_inputs: &[Vec<RawInput>],
    raw_cache: &mut [Option<Vec<Raw>>],
    errfun: &ErrFun,
    dada_params: &mut DadaParams,
    nq: usize,
    verbose: bool,
) -> io::Result<Vec<f64>> {
    let align_params = dada_params.align;
    // R: erri = matrix(1, nrow=16, ncol=...). Every entry, not just diagonal.
    let saved_max_clust = dada_params.max_clust;
    dada_params.err_mat = vec![1.0f64; 16 * nq];
    dada_params.err_ncol = nq;
    dada_params.max_clust = 1;

    let sample_results: Vec<(usize, Result<Vec<u32>, String>)> = all_inputs
        .par_iter()
        .zip(raw_cache.par_iter_mut())
        .enumerate()
        .map(|(si, (inputs, cache_slot))| {
            let cached = cache_slot.take();
            let outcome =
                dada_uniques_cached(inputs, cached, dada_params).map(|(result, reused_raws)| {
                    *cache_slot = Some(reused_raws);
                    build_trans_mat(inputs, &result, &align_params, nq)
                });
            (si, outcome)
        })
        .collect();

    dada_params.max_clust = saved_max_clust;

    let mut sample_trans: Vec<Vec<u32>> = Vec::new();
    for (si, outcome) in sample_results {
        match outcome {
            Ok(t) => sample_trans.push(t),
            Err(e) if verbose => {
                eprintln!(
                    "[learn_errors] init pass sample={}: dada_uniques failed: {}",
                    si + 1,
                    e
                );
            }
            Err(_) => {}
        }
    }

    if sample_trans.is_empty() {
        return Err(io::Error::other(
            "All DADA runs failed during init pass; cannot estimate error model",
        ));
    }

    let refs: Vec<(&[u32], usize)> = sample_trans.iter().map(|t| (t.as_slice(), nq)).collect();
    let (acc_trans, _) = accumulate_trans(&refs);

    let mut new_err = errfun
        .apply(&acc_trans, nq)
        .map_err(|e| io::Error::other(format!("init-pass errfun failed: {e}")))?;

    // R: err[c(1,6,11,16),] <- 1.0 — force self-transitions to 1 at every q.
    // Row indices are 1-based in R; here they are 0,5,10,15 (A2A, C2C, G2G, T2T).
    for &row in &[0usize, 5, 10, 15] {
        for q in 0..nq {
            new_err[row * nq + q] = 1.0;
        }
    }

    if verbose {
        eprintln!("[learn_errors] init pass complete (R nconsist=0)");
    }

    Ok(new_err)
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
    // Build Raw stubs for cluster centers (no quality scores needed).
    let center_raws: Vec<Raw> = result
        .clusters
        .iter()
        .map(|c| Raw::new(c.sequence.clone(), None, 0, false))
        .collect();

    // Parallel reduction: each rayon task builds a private 16*nq trans
    // matrix for a subset of inputs, then we sum them. Each task gets its
    // own AlignBuffers via map_init so buffer reuse holds across the
    // task's items without locking.
    inputs
        .par_iter()
        .enumerate()
        .fold(
            || (vec![0u32; 16 * nq], AlignBuffers::new()),
            |(mut trans, mut buf), (i, inp)| {
                let ci = match result.map[i] {
                    Some(ci) => ci,
                    None => return (trans, buf),
                };
                let quals = match &inp.quals {
                    Some(q) => q.as_slice(),
                    None => return (trans, buf),
                };

                let seq: Vec<u8> = inp.seq.bytes().map(nt_encode).collect();
                let raw_query = Raw::new(seq, Some(quals), inp.abundance, false);
                let raw_center = &center_raws[ci];

                // Align center (ref = al0) against the raw (query = al1).
                if raw_align_with_buf(raw_center, &raw_query, align_params, &mut buf).is_none() {
                    return (trans, buf);
                }
                let al_ref: &[u8] = &buf.al0;
                let al_qry: &[u8] = &buf.al1;
                let reads = inp.abundance;
                let mut qpos = 0usize;

                for alpos in 0..al_ref.len() {
                    let nt0 = al_ref[alpos];
                    let nt1 = al_qry[alpos];
                    let qry_is_nt = matches!(nt1, 1..=4);
                    if matches!(nt0, 1..=4) && qry_is_nt {
                        let q = (quals[qpos].round() as usize).min(nq - 1);
                        let row = (nt0 as usize - 1) * 4 + (nt1 as usize - 1);
                        trans[row * nq + q] = trans[row * nq + q].saturating_add(reads);
                    }
                    if nt1 != b'-' {
                        qpos += 1;
                    }
                }
                (trans, buf)
            },
        )
        .map(|(t, _buf)| t)
        .reduce(
            || vec![0u32; 16 * nq],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b.iter()) {
                    *x = x.saturating_add(*y);
                }
                a
            },
        )
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
/// - `dada_params`   — DADA2 algorithm parameters (error matrix is overwritten each iteration;
///   alignment parameters are taken from `dada_params.align`)
/// - `max_consist`   — maximum self-consistency iterations (R default: 10)
/// - `diag`          — diagnostic/trace output options (verbosity, directories, trace config)
///
/// # Returns
/// [`LearnErrorsResult`] with transition counts and error-rate matrices, or an I/O error.
pub fn learn_errors(
    all_inputs: Vec<Vec<RawInput>>,
    errfun: &ErrFun,
    mut dada_params: DadaParams,
    max_consist: usize,
    diag: LearnDiagOptions<'_>,
) -> io::Result<LearnErrorsResult> {
    let LearnDiagOptions {
        verbose,
        diag_dir,
        cluster_trace_dir,
        trace_params,
    } = diag;
    if all_inputs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No input samples provided",
        ));
    }

    type SampleResults = Vec<(usize, Result<(Vec<u32>, crate::dada::DadaResult), String>)>;

    let nq = detect_nq(&all_inputs);

    if verbose {
        eprintln!(
            "[learn_errors] {} sample(s), nq={}, max_consist={}",
            all_inputs.len(),
            nq,
            max_consist,
        );
    }

    dada_params.err_ncol = nq;

    // Per-sample cache of `Vec<Raw>` (encoded sequences + k-mer vectors),
    // reused across outer iterations so we don't re-encode inputs or
    // recompute k-mers on every self-consistency pass. Populated by the
    // init pass below; reused across the main loop.
    let mut raw_cache: Vec<Option<Vec<Raw>>> = (0..all_inputs.len()).map(|_| None).collect();

    // ---- Init pass (R nconsist=0) ----
    // Run dada with all-1.0 err_mat and MAX_CLUST=1 on every sample, fit
    // errfun on the resulting trans, and force diagonal self-transitions
    // to 1.0. This produces iter-1's err_in with the same shape and
    // per-Q calibration that R uses.
    let mut err = run_init_pass(
        &all_inputs,
        &mut raw_cache,
        errfun,
        &mut dada_params,
        nq,
        verbose,
    )?;

    let mut converged = false;
    let mut iterations = 0usize;
    let mut trans_final = vec![0u32; 16 * nq];
    let mut err_in_final = err.clone();
    let mut err_out_final = err.clone();
    let mut stop_reason = StopReason::MaxConsistReached;
    let mut best_max_delta = f64::INFINITY;
    let mut iters_since_best = 0usize;

    for iter in 0..max_consist {
        iterations = iter + 1;

        // ---- Run DADA on each sample ----
        // Set the error matrix once per iteration (shared across all samples).
        dada_params.err_mat = err.clone();

        // Parallel pass: each sample runs dada_uniques independently.
        // Rayon's work-stealing distributes threads across samples; if only
        // one sample exists all threads serve its inner b_compare_parallel.
        //
        // `raw_cache` is zipped in so each sample's Raw Vec flows back out of
        // DADA and into the next iteration's call without reallocation.
        let collect_diags = diag_dir.is_some();
        let sample_results: SampleResults = all_inputs
            .par_iter()
            .zip(raw_cache.par_iter_mut())
            .enumerate()
            .map(|(si, (inputs, cache_slot))| {
                let cached = cache_slot.take();
                let outcome = dada_uniques_cached(inputs, cached, &dada_params).map(
                    |(result, reused_raws)| {
                        *cache_slot = Some(reused_raws);
                        let t = build_trans_mat(inputs, &result, &dada_params.align, nq);
                        (t, result)
                    },
                );
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
                            n_initial: 0,
                            n_abundance: 0,
                            n_prior: 0,
                            n_singleton: 0,
                            nalign: result.nalign,
                            nshroud: result.nshroud,
                        };
                        for c in &result.clusters {
                            match c.birth_type {
                                BirthType::Initial => diag.n_initial += 1,
                                BirthType::Abundance => diag.n_abundance += 1,
                                BirthType::Prior => diag.n_prior += 1,
                                BirthType::Singleton => diag.n_singleton += 1,
                            }
                        }
                        sample_diags.push(diag);
                    }

                    if let Some(dir) = cluster_trace_dir {
                        let path = dir.join(format!(
                            "cluster_iter_{:03}_sample_{:03}.json",
                            iter + 1,
                            si + 1,
                        ));
                        if let Err(e) = crate::cluster_trace::write_trace(
                            &path,
                            &format!("sample_{:03}", si + 1),
                            Some(iter + 1),
                            &all_inputs[si],
                            &result,
                            Some(&dada_params.err_mat),
                            nq,
                            trace_params,
                            true, // compact: trace files can be large; minify
                        ) {
                            eprintln!(
                                "[learn_errors] WARN: failed to write cluster trace {}: {e}",
                                path.display(),
                            );
                        } else if verbose {
                            eprintln!("[learn_errors] cluster trace written to {}", path.display());
                        }
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
            return Err(io::Error::other(
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
        let new_err = errfun
            .apply(&acc_trans, nq)
            .map_err(|e| io::Error::other(format!("errfun failed: {e}")))?;

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
            let json = serde_json::to_string_pretty(&diag).map_err(io::Error::other)?;
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
            stop_reason = StopReason::Converged;
            if verbose {
                eprintln!("[learn_errors] converged after {} iteration(s)", iter + 1);
            }
            break;
        }

        // ---- Stall detection: bail if max_delta hasn't set a new
        // best-so-far for STALL_PATIENCE consecutive iterations.  This
        // catches oscillation (e.g. k-mer screen too loose → err_out
        // cycles between a small set of values). Without this we'd run
        // through max_consist iterations producing identical-cluster
        // but slightly-different err_out tables. ----
        if max_delta < best_max_delta * STALL_REL_IMPROVE {
            best_max_delta = max_delta;
            iters_since_best = 0;
        } else {
            iters_since_best += 1;
            if iters_since_best >= STALL_PATIENCE {
                stop_reason = StopReason::Stalled;
                if verbose {
                    eprintln!(
                        "[learn_errors] stalled at iter {}: max_delta has not improved \
                         on best-so-far ({:.2e}) for {} iteration(s); stopping early",
                        iter + 1,
                        best_max_delta,
                        STALL_PATIENCE,
                    );
                }
                break;
            }
        }

        err = new_err;
    }

    if !converged && verbose && stop_reason == StopReason::MaxConsistReached {
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
        stop_reason,
    })
}
