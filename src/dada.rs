//! Core DADA2 algorithm — Rcpp-free entry points.
//!
//! Ports the logic from `Rmain.cpp`, stripping all R/Rcpp bindings:
//!
//! - `dada_uniques`: validates input, constructs `Raw` objects, runs DADA,
//!   computes final p-values, and returns a `DadaResult`.
//! - `run_dada`: the inner algorithm loop (initial compare → bud/shuffle
//!   iterations → p-value updates).
//!
//! ## Removed from the C++ original
//! - `SSE`/`X64` SIMD dispatch — LLVM auto-vectorises scalar loops.
//! - `Rcpp::checkUserInterrupt()` — no R event loop.
//! - `b_make_*` output formatters — those produced R data frames; callers
//!   should build their own output from `DadaResult`.
//! - `final_consensus` — kept in `DadaParams` for future use but currently
//!   has no effect (mirrors C++ where it is passed through but unused in the
//!   loop itself).

use crate::cluster::{b_bud, b_compare, b_compare_parallel, b_shuffle2};
use crate::containers::{B, BirthType, Raw};
use crate::kmers::{raw_assign_kmers, KMER_SIZE_MAX, KMER_SIZE_MIN};
use crate::misc::nt_encode;
use crate::nwalign::AlignParams;
use crate::pval::{b_p_update, calc_pA};

/// Maximum shuffle iterations before giving up on convergence.
/// Matches C++ `MAX_SHUFFLE`.
const MAX_SHUFFLE: usize = 10;

/// Maximum accepted sequence length (buffer guard from C++ `SEQLEN`).
const SEQLEN: usize = 9999;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// All tuning parameters for the DADA2 algorithm.
pub struct DadaParams {
    /// Alignment parameters (method selection, scoring, band).
    pub align: AlignParams,
    /// Flat row-major error rate matrix with shape 16 × `err_ncol`.
    /// Row `r * 16 + q` holds the error probability for transition `r` at
    /// quality score `q`, where transitions are indexed as ref_nt*4 + query_nt.
    pub err_mat: Vec<f64>,
    /// Number of quality-score columns in `err_mat`.
    pub err_ncol: usize,
    /// Significance threshold for abundance-based cluster splitting.
    pub omega_a: f64,
    /// Significance threshold for prior-sequence splitting.
    pub omega_p: f64,
    /// Per-raw p-value threshold below which a raw is not corrected to its
    /// cluster center (maps to `NA` in the read-to-cluster assignment).
    pub omega_c: f64,
    /// Apply singleton detection (detect_singletons in C++).
    pub detect_singletons: bool,
    /// Maximum number of clusters. `0` means unlimited (use `nraw`).
    pub max_clust: usize,
    /// Minimum fold-enrichment above expected for a raw to bud a new cluster.
    pub min_fold: f64,
    /// Minimum Hamming distance for a raw to bud a new cluster.
    pub min_hamming: u32,
    /// Minimum read abundance for a raw to bud a new cluster.
    pub min_abund: u32,
    /// Whether quality scores are available and should be used.
    pub use_quals: bool,
    /// Reserved for future use (matches C++ `final_consensus` parameter).
    pub final_consensus: bool,
    /// Use Rayon for parallel comparisons.
    pub multithread: bool,
    /// Write progress to stderr.
    pub verbose: bool,
    /// Greedy mode: lock Raws whose expected abundance already exceeds observed.
    pub greedy: bool,
}

/// A single unique sequence with its abundance, optional quality profile,
/// and prior flag.
pub struct RawInput {
    /// ASCII nucleotide sequence (A/C/G/T/N, upper or lower case).
    pub seq: String,
    /// Number of reads with this exact sequence.
    pub abundance: u32,
    /// When `true`, this sequence is presumed genuine regardless of p-value.
    pub prior: bool,
    /// Per-position Phred quality scores. Must have the same length as `seq`
    /// if provided. Supply `None` when quality data is unavailable.
    pub quals: Option<Vec<f64>>,
}

/// Per-cluster summary produced by `dada_uniques`.
pub struct ClusterSummary {
    /// Integer-encoded representative (center) sequence.
    pub sequence: Vec<u8>,
    /// Total reads assigned to this cluster.
    pub reads: u32,
    /// Indices (into the input `RawInput` slice) of member Raws.
    pub members: Vec<usize>,
    pub birth_type: BirthType,
    /// Index of the parent cluster that this one was split from.
    pub birth_from: u32,
    /// Bonferroni-corrected p-value that triggered this cluster's creation.
    pub birth_pval: f64,
    /// Fold-enrichment above expectation at birth.
    pub birth_fold: f64,
    /// Expected read count at the time of birth.
    pub birth_e: f64,
}

/// Output of `dada_uniques`.
pub struct DadaResult {
    /// One entry per cluster, in cluster order (cluster 0 is the initial
    /// catch-all; clusters 1+ are buds).
    pub clusters: Vec<ClusterSummary>,
    /// For each input Raw (in input order), the index of the cluster it maps
    /// to.  `None` means the Raw's final p-value fell below `omega_c` and it
    /// was not corrected to any center.
    pub map: Vec<Option<usize>>,
    /// Final abundance p-value for each input Raw (in input order).
    pub pvals: Vec<f64>,
    /// Total pairwise alignments performed.
    pub nalign: u32,
    /// Comparisons screened out by k-mer distance.
    pub nshroud: u32,
}

// ---------------------------------------------------------------------------
// dada_uniques
// ---------------------------------------------------------------------------

/// Validate inputs, construct `Raw` objects, run DADA, compute final
/// p-values, and return a `DadaResult`.
///
/// Equivalent to the logic in C++ `dada_uniques`, minus the Rcpp layer and
/// the R-specific output formatters.
pub fn dada_uniques(inputs: &[RawInput], params: &DadaParams) -> Result<DadaResult, String> {
    // ---- Input validation ----
    let nraw = inputs.len();
    if nraw == 0 {
        return Err("Zero input sequences.".into());
    }
    let maxlen = inputs.iter().map(|r| r.seq.len()).max().unwrap_or(0);
    let minlen = inputs.iter().map(|r| r.seq.len()).min().unwrap_or(0);

    if maxlen >= SEQLEN {
        return Err(format!(
            "Input sequences exceed the maximum allowed length ({SEQLEN})."
        ));
    }
    let k = params.align.kmer_size;
    if !(KMER_SIZE_MIN..=KMER_SIZE_MAX).contains(&k) {
        return Err(format!(
            "kmer_size {k} out of supported range ({KMER_SIZE_MIN}..={KMER_SIZE_MAX})."
        ));
    }
    if minlen <= k {
        return Err(format!(
            "All input sequences must be longer than the k-mer size ({k})."
        ));
    }
    if params.err_mat.len() != 16 * params.err_ncol {
        return Err(format!(
            "Error matrix length {} does not match 16 × {} = {}.",
            params.err_mat.len(),
            params.err_ncol,
            16 * params.err_ncol
        ));
    }
    let has_quals = inputs.iter().any(|r| r.quals.is_some());
    if has_quals {
        for (i, inp) in inputs.iter().enumerate() {
            match &inp.quals {
                Some(q) if q.len() != inp.seq.len() => {
                    return Err(format!(
                        "Sequence {i}: quality length {} does not match sequence length {}.",
                        q.len(),
                        inp.seq.len()
                    ));
                }
                _ => {}
            }
        }
    }

    // ---- Construct Raw objects ----
    let mut raws: Vec<Raw> = inputs
        .iter()
        .enumerate()
        .map(|(i, inp)| {
            let seq: Vec<u8> = inp.seq.bytes().map(nt_encode).collect();
            let qual = if has_quals { inp.quals.as_deref() } else { None };
            let mut raw = Raw::new(seq, qual, inp.abundance, inp.prior);
            raw.index = i as u32;
            raw
        })
        .collect();

    // ---- Assign k-mers ----
    if params.align.use_kmers {
        for raw in &mut raws {
            raw_assign_kmers(raw, k);
        }
    }

    // ---- Run core algorithm ----
    let mut b = run_dada(raws, params);

    // ---- Final per-raw p-value pass ----
    // Determines raw->correct, which controls the read-to-cluster map.
    let mut pvals = vec![0.0f64; nraw];
    for ci in 0..b.clusters.len() {
        let members: Vec<usize> = b.clusters[ci].raws.clone();
        let center_idx = b.clusters[ci].center;
        let ci_reads = b.clusters[ci].reads;
        for raw_idx in members {
            let is_center = Some(raw_idx) == center_idx;
            let (p, correct) = if is_center {
                (1.0, true)
            } else {
                let lambda = b.raws[raw_idx].comp.lambda;
                let p = calc_pA(b.raws[raw_idx].reads, lambda * ci_reads as f64, true);
                let correct = p >= params.omega_c;
                (p, correct)
            };
            b.raws[raw_idx].p = p;
            b.raws[raw_idx].correct = correct;
            pvals[b.raws[raw_idx].index as usize] = p;
        }
    }

    // ---- Build map ----
    let mut map: Vec<Option<usize>> = vec![None; nraw];
    for ci in 0..b.clusters.len() {
        for &raw_idx in &b.clusters[ci].raws {
            if b.raws[raw_idx].correct {
                map[b.raws[raw_idx].index as usize] = Some(ci);
            }
        }
    }

    // ---- Build cluster summaries ----
    let clusters = b
        .clusters
        .iter()
        .map(|bi| ClusterSummary {
            sequence: bi.seq.clone(),
            reads: bi.reads,
            members: bi.raws.clone(),
            birth_type: bi.birth_type.clone(),
            birth_from: bi.birth_from,
            birth_pval: bi.birth_pval,
            birth_fold: bi.birth_fold,
            birth_e: bi.birth_e,
        })
        .collect();

    Ok(DadaResult {
        clusters,
        map,
        pvals,
        nalign: b.nalign,
        nshroud: b.nshroud,
    })
}

// ---------------------------------------------------------------------------
// run_dada
// ---------------------------------------------------------------------------

/// Core DADA2 algorithm loop.
///
/// 1. Initialises a single cluster containing all Raws.
/// 2. Compares all Raws to cluster 0 (no k-mer screen: `kdist_cutoff = 1.0`).
/// 3. Computes initial abundance p-values.
/// 4. Iterates: bud → compare new cluster → shuffle to convergence →
///    update p-values — until no significant bud is found or `max_clust` is
///    reached.
///
/// Returns the final partition `B`.  Callers are responsible for any
/// post-processing (final p-values, map construction, output formatting).
///
/// Equivalent to C++ `run_dada`.
pub fn run_dada(raws: Vec<Raw>, params: &DadaParams) -> B {
    let mut bb = B::new(raws, params.omega_a, params.omega_p, params.use_quals);

    // Initial compare: no k-mer distance screen so that cluster 0 accumulates
    // comparisons for every Raw (required by b_shuffle2).
    let init_params = AlignParams { kdist_cutoff: 1.0, ..params.align };

    if params.multithread {
        b_compare_parallel(
            &mut bb, 0,
            &params.err_mat, params.err_ncol,
            &init_params,
            params.greedy,
        );
    } else {
        b_compare(
            &mut bb, 0,
            &params.err_mat, params.err_ncol,
            &init_params,
            params.greedy,
            params.verbose,
        );
    }
    b_p_update(&mut bb, params.greedy, params.detect_singletons);

    let max_clust = if params.max_clust == 0 {
        bb.raws.len()
    } else {
        params.max_clust
    };

    while bb.clusters.len() < max_clust {
        let newi = match b_bud(
            &mut bb,
            params.min_fold,
            params.min_hamming,
            params.min_abund,
            params.verbose,
        ) {
            Some(i) => i,
            None => break,
        };

        if params.verbose {
            eprint!("[dada] New Cluster C{newi}:");
        }

        if params.multithread {
            b_compare_parallel(
                &mut bb, newi,
                &params.err_mat, params.err_ncol,
                &params.align,
                params.greedy,
            );
        } else {
            b_compare(
                &mut bb, newi,
                &params.err_mat, params.err_ncol,
                &params.align,
                params.greedy,
                params.verbose,
            );
        }

        // Shuffle until stable or MAX_SHUFFLE reached.
        let mut nshuffle = 0usize;
        loop {
            let shuffled = b_shuffle2(&mut bb);
            if params.verbose {
                eprint!("S");
            }
            nshuffle += 1;
            if !shuffled || nshuffle >= MAX_SHUFFLE {
                break;
            }
        }
        if params.verbose && nshuffle >= MAX_SHUFFLE {
            eprintln!(" Warning: Reached maximum ({MAX_SHUFFLE}) shuffles.");
        }

        b_p_update(&mut bb, params.greedy, params.detect_singletons);
    }

    if params.verbose {
        eprintln!(
            "[dada] ALIGN: {} aligns, {} shrouded ({} raw).",
            bb.nalign,
            bb.nshroud,
            bb.raws.len()
        );
    }

    bb
}
