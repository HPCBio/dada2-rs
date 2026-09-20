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

use rayon::prelude::*;

use crate::cluster::{
    CandIndex, ShuffleCarry, b_bud_incremental, b_compare, b_compare_parallel, b_shuffle_converge,
    index_add_cluster,
};
use crate::containers::{B, BirthType, Raw, Sub};
use crate::error::{
    BirthSubRecord, ClusterStats, birth_sub_records, cluster_quality, cluster_stats,
    transition_counts,
};
use crate::kmers::{KMER_SIZE_MAX, KMER_SIZE_MIN, assign_kmer_order, raw_assign_kmers};
use crate::metrics::{MeasureLevel, RawCounters, RunMetrics, RunShape};
use crate::minimizers::{self, MINIMIZER_K_MAX, MINIMIZER_K_MIN, MINIMIZER_W_MAX, MINIMIZER_W_MIN};
use crate::misc::nt_encode;
use crate::nwalign::ScreenBackend;
use crate::nwalign::{AlignBuffers, AlignParams, sub_new_with_buf};
use crate::progress::ProgressRecord;
use crate::pval::{b_p_update, calc_pA};
use crate::rec_print;
use std::sync::OnceLock;
use std::time::Duration;

/// Maximum shuffle iterations before giving up on convergence.
/// Matches C++ `MAX_SHUFFLE`.
const MAX_SHUFFLE: usize = 10;

/// Maximum accepted sequence length (buffer guard from C++ `SEQLEN`).
const SEQLEN: usize = 9999;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// All tuning parameters for the DADA2 algorithm.
///
/// `Clone` so `dada_uniques_cached` can hand `run_dada` a copy carrying an
/// error matrix extended to cover the data's quality range (issue #102) without
/// forcing every caller to build one.
#[derive(Clone)]
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
    #[allow(dead_code)]
    pub final_consensus: bool,
    /// Use Rayon for parallel comparisons.
    pub multithread: bool,
    /// Write progress to stderr.
    pub verbose: bool,
    /// Sample label prefixed to each bud-round progress record.
    ///
    /// Set only when samples are denoised concurrently. With a single
    /// `run_dada` in flight there is nothing to disambiguate, so the records
    /// stay byte-identical to the R output they were ported from (#172).
    pub progress_tag: Option<String>,
    /// How much instrumentation to collect, independent of `verbose`.
    ///
    /// `verbose` decides what is *printed*; this decides what is *measured*.
    /// They are separate because [`MeasureLevel::Attribution`] costs 2-4
    /// `Instant::now()` calls per comparison, so `--metrics-json` can collect
    /// the free tier on a production run without paying for it (issue #162).
    pub measure: MeasureLevel,
    /// Greedy mode: lock Raws whose expected abundance already exceeds observed.
    pub greedy: bool,
    /// Compute auxiliary outputs (R DADA2 parity: `$clustering`, `$birth_subs`,
    /// `$subqual`, `$clusterquals`).
    ///
    /// When `true`, `dada_uniques` runs an extra final-subs alignment pass to
    /// produce per-cluster substitution stats (n0/n1/nunq/birth_qave/post-hoc
    /// p-value), per-cluster mean quality at each position, per-cluster birth
    /// substitution records, and a 16 × `err_ncol` transition-by-quality
    /// matrix. Defaults to `false` — the pass costs roughly one alignment per
    /// raw against its cluster center.
    pub aux_outputs: bool,
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
    /// Per-position *integer* Phred quality SUM over the reads in this unique
    /// (deferred division, issue #23). Must have the same length as `seq` when
    /// present; `None` when quality data is unavailable.
    ///
    /// We store the integer sum rather than the f64 mean for compactness: 4
    /// bytes/position instead of 8, which roughly halves the dominant resident
    /// field across all samples held in memory. The mean Phred the rest of the
    /// pipeline consumes is recovered on demand via [`RawInput::mean_quals`] as
    /// `sum / abundance`.
    ///
    /// This is bit-identical to storing the f64 mean: the sum is exact, and
    /// `abundance` is the per-position read count for every covered position
    /// because dereplication groups by exact full-sequence identity (all reads
    /// in a unique are byte-identical → identical length → full coverage). If
    /// approximate/length-collapsing dereplication is ever introduced, this
    /// representation must carry per-position counts too.
    pub quals: Option<Vec<u32>>,
}

/// Per-cluster summary produced by `dada_uniques`.
pub struct ClusterSummary {
    /// Integer-encoded representative (center) sequence.
    pub sequence: Vec<u8>,
    /// Total reads assigned to this cluster.
    pub reads: u32,
    /// Indices (into the input `RawInput` slice) of member Raws.
    pub members: Vec<usize>,
    /// Hamming distance to the cluster center for each entry in `members`
    /// (parallel slice; same length as `members`).
    pub member_hammings: Vec<u32>,
    /// Per-member λ (transition-probability product against the center).
    pub member_lambdas: Vec<f64>,
    /// Per-member final abundance p-value.
    pub member_pvals: Vec<f64>,
    pub birth_type: BirthType,
    /// Index of the parent cluster that this one was split from.
    pub birth_from: u32,
    /// Bonferroni-corrected p-value that triggered this cluster's creation.
    pub birth_pval: f64,
    /// Fold-enrichment above expectation at birth.
    pub birth_fold: f64,
    /// Expected read count at the time of birth.
    pub birth_e: f64,
    /// Hamming distance from this cluster's center to its birth-parent's
    /// center at the time of budding (0 for the initial cluster).
    pub birth_hamming: u32,
}

/// R DADA2 parity outputs computed when `DadaParams::aux_outputs` is set.
///
/// Mirrors the additional fields R's `dada()` returns alongside its main
/// clustering result: `$clustering`, `$birth_subs`, `$subqual`, `$clusterquals`.
pub struct DadaAux {
    /// Per-cluster summary stats (R `$clustering`): n0, n1, nunq, birth_qave,
    /// post-hoc abundance p-value.
    pub cluster_stats: Vec<ClusterStats>,
    /// Per-cluster read-weighted mean quality at each reference position
    /// (R `$clusterquals`). Outer length = nclust, inner length = `cluster_quality_maxlen`.
    /// Positions outside the cluster center or with no covering reads are NaN.
    pub cluster_quality: Vec<Vec<f64>>,
    /// Maximum reference length used to size each `cluster_quality` row.
    pub cluster_quality_maxlen: usize,
    /// Per-substitution records from each cluster's birth alignment
    /// (R `$birth_subs`).
    pub birth_subs: Vec<BirthSubRecord>,
    /// Flat row-major 16 × `transitions_ncol` transition-by-quality count
    /// matrix (R `$subqual`). `result[t*ncol + q]` = reads with transition `t`
    /// (ref_nt*4 + query_nt) at quality `q`.
    pub transitions: Vec<u32>,
    /// Number of quality columns in `transitions` (1 when no quals were used).
    pub transitions_ncol: usize,
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
    #[allow(dead_code)]
    pub pvals: Vec<f64>,
    /// Total pairwise alignments performed.
    pub nalign: u64,
    /// Comparisons screened out by k-mer distance.
    pub nshroud: u64,
    /// Auxiliary R-DADA2-parity outputs. `Some` only when
    /// `DadaParams::aux_outputs` was true.
    pub aux: Option<DadaAux>,
    /// Bonferroni divisor used in the abundance p-value test
    /// (`p_a = min_p * nraw`); equals the unique-input count. Constant across
    /// all births in this run. Lets a consumer recover the abundance-scale
    /// p-value of a `Prior`-born cluster as `birth_pval * nraw`.
    pub nraw: u32,
    /// `omega_a` / `omega_p` thresholds in force for this run, echoed so the
    /// trace/consumer can compare against the exact values used rather than a
    /// presumed default.
    pub omega_a: f64,
    pub omega_p: f64,
    /// Instrumentation from this invocation, when `DadaParams::measure` asked
    /// for any. Reaches the user through `--metrics-json`; nothing serializes
    /// it into a subcommand's own output (issue #162).
    pub metrics: Option<RunMetrics>,
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
    let (result, _raws) = dada_uniques_cached(inputs, None, params)?;
    Ok(result)
}

/// Variant of [`dada_uniques`] that accepts a pre-built `Vec<Raw>` (with
/// k-mer vectors already populated) via `cached` to skip per-iteration
/// setup. Returns the `Vec<Raw>` alongside the result so it can be fed back
/// into the next call.
///
/// Pass `None` on the first call; pass `Some(raws)` returned from a prior
/// call on subsequent iterations. The caller is responsible for ensuring
/// `inputs` hasn't changed between calls (we do not re-validate the seq
/// bytes when a cache is reused — only the mutable iteration state is
/// reset).
///
/// Used by `learn_errors` to avoid re-encoding sequences and rebuilding
/// k-mer vectors on every self-consistency iteration.
pub fn dada_uniques_cached(
    inputs: &[RawInput],
    cached: Option<Vec<Raw>>,
    params: &DadaParams,
) -> Result<(DadaResult, Vec<Raw>), String> {
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
    if params.align.screen_backend == ScreenBackend::Minimizer {
        let (mk, mw) = (params.align.minimizer_k, params.align.minimizer_w);
        if !(MINIMIZER_K_MIN..=MINIMIZER_K_MAX).contains(&mk) {
            return Err(format!(
                "minimizer_k {mk} out of supported range ({MINIMIZER_K_MIN}..={MINIMIZER_K_MAX})."
            ));
        }
        if !(MINIMIZER_W_MIN..=MINIMIZER_W_MAX).contains(&mw) {
            return Err(format!(
                "minimizer_w {mw} out of supported range ({MINIMIZER_W_MIN}..={MINIMIZER_W_MAX})."
            ));
        }
        // Not an error: a sequence shorter than one full window gets an empty
        // sketch, and `screen_dist_minimizer` deliberately fails *open* on it
        // (aligns rather than screens out). That is correct but silent, so say
        // so once here rather than letting an unexplained loss of screening
        // efficiency look like a performance bug.
        if minlen < mk + mw - 1 && params.verbose {
            eprintln!(
                "[dada] warning: shortest input is {minlen} nt, below the minimizer \
                 window ({} nt for k={mk}, w={mw}); those sequences bypass the screen \
                 and are always aligned",
                mk + mw - 1
            );
        }
    }
    // Reject non-ACGT here rather than in the inner loop (issue #101). `N` in
    // particular is a user error with a known remedy — DADA2's workflow requires
    // `maxN=0` — so it deserves the same treatment as the other input problems
    // above, not a panic from a rayon worker inside `compute_lambda` that names
    // neither the sample nor the sequence. `U` is accepted because `nt_encode`
    // folds it to `T`; everything else (including `N` and `-`) has no transition
    // row in the error matrix and cannot be scored.
    for (i, inp) in inputs.iter().enumerate() {
        if let Some(off) = inp.seq.bytes().position(|b| {
            !matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
            )
        }) {
            let ch = inp.seq.as_bytes()[off] as char;
            return Err(format!(
                "Sequence {i} contains a non-ACGT base {ch:?} at position {off} \
                 (0-based). DADA2 cannot assign an error rate to it. Filter the \
                 reads first — `filter-and-trim --max-n 0` removes reads \
                 containing N, as R DADA2's workflow requires."
            ));
        }
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

    // ---- Extend the error model if the data has higher quality than it covers ----
    // Mirrors R DADA2 (`dada.R:302-312`), which repeats the last column up to the
    // maximum observed quality rather than failing: a model learned on one run,
    // or on a `--nbases` subsample that missed the top quality bin, would
    // otherwise abort on data R processes fine. Without this the index runs off
    // the end of `err_mat` inside `compute_lambda` (issue #102).
    //
    // Two deliberate departures from R:
    //  - the warning is NOT verbose-gated. R extends silently by default;
    //    extrapolated error rates shift low-abundance calls, so a run that does
    //    it should say so.
    //  - `qmax` uses `round`, matching what our own indexing does
    //    (`Raw::from_qual_sums`) and what `learn_errors::detect_nq` sizes the
    //    matrix with. R uses `ceiling` (`dada.R:290`), but adopting that here
    //    would fire on the STANDARD workflow: a mean of 39.4 rounds to 39, so
    //    detect_nq produces nq=40, while ceiling gives qmax=40 >= 40 and would
    //    extend a matrix that already covers the data. The trigger has to match
    //    the index that can actually overflow.
    let extended_err: Option<(Vec<f64>, usize)> = if has_quals && params.err_ncol > 0 {
        let qmax = inputs
            .iter()
            .filter_map(|inp| {
                let q = inp.quals.as_ref()?;
                let c = inp.abundance.max(1) as f64;
                q.iter().map(|&s| (s as f64 / c).round() as usize).max()
            })
            .max()
            .unwrap_or(0);
        if qmax >= params.err_ncol {
            let old_ncol = params.err_ncol;
            let new_ncol = qmax + 1;
            let mut ext = vec![0.0f64; 16 * new_ncol];
            for t in 0..16 {
                let row = &params.err_mat[t * old_ncol..(t + 1) * old_ncol];
                ext[t * new_ncol..t * new_ncol + old_ncol].copy_from_slice(row);
                // Repeat the last learned column across the new ones.
                let last = row[old_ncol - 1];
                for q in old_ncol..new_ncol {
                    ext[t * new_ncol + q] = last;
                }
            }
            eprintln!(
                "dada2-rs: warning: input has quality up to Q{qmax} but the error \
                 model covers Q0-Q{}. Extending the model by repeating its last \
                 column (Q{}) for Q{}-Q{qmax}. Rates for the extended columns are \
                 extrapolated, not learned; re-run learn-errors on this data to \
                 avoid it.",
                old_ncol - 1,
                old_ncol - 1,
                old_ncol,
            );
            Some((ext, new_ncol))
        } else {
            None
        }
    } else {
        None
    };
    // Owned params only when the matrix had to grow, so the common path keeps
    // borrowing the caller's and copies nothing.
    let owned_params: Option<DadaParams> = extended_err.map(|(err_mat, err_ncol)| DadaParams {
        err_mat,
        err_ncol,
        ..params.clone()
    });
    let params: &DadaParams = owned_params.as_ref().unwrap_or(params);

    // ---- Build or reset Raw objects ----
    let raws: Vec<Raw> = match cached {
        Some(mut raws) if raws.len() == nraw => {
            // Reuse path: reset per-iteration mutable state only. seq/qual/
            // kmer vectors persist across iterations.
            for raw in &mut raws {
                raw.reset_for_iteration();
            }
            raws
        }
        _ => {
            // Fresh build: encode sequences and populate k-mer vectors.
            let mut raws: Vec<Raw> = inputs
                .iter()
                .enumerate()
                .map(|(i, inp)| {
                    let seq: Vec<u8> = inp.seq.bytes().map(nt_encode).collect();
                    let qual = if has_quals {
                        inp.quals.as_deref()
                    } else {
                        None
                    };
                    let mut raw = Raw::from_qual_sums(seq, qual, inp.abundance, inp.prior);
                    raw.index = i as u32;
                    raw
                })
                .collect();

            if params.align.use_kmers {
                match params.align.screen_backend {
                    ScreenBackend::Kmer => {
                        for raw in &mut raws {
                            raw_assign_kmers(raw, k);
                        }
                    }
                    ScreenBackend::Minimizer => {
                        let (mk, mw) = (params.align.minimizer_k, params.align.minimizer_w);
                        for raw in &mut raws {
                            // `kord` is still needed: the gapless fast path in
                            // `raw_align_dp` keys off the k-order vector, not off
                            // the screen. `kmer8` is deliberately left `None` —
                            // building both screens would pay twice for one.
                            raw.minimizers = Some(minimizers::sketch(&raw.seq, mk, mw));
                            raw.kord = Some(assign_kmer_order(&raw.seq, k));
                            // ...except under audit, which compares the two
                            // screens pair-by-pair and therefore needs both.
                            if params.align.screen_audit {
                                raw_assign_kmers(raw, k);
                            }
                        }
                    }
                }
            }

            // Resident-footprint accounting for memory profiling (#32). Fires
            // once per fresh Raw build (not on learn-errors' cached-reuse
            // iterations). Per-Raw resident = seq + qual + the k-mer screen
            // vectors (kmer8 = 4^k bytes; kord = (len-k+1)×2 bytes). The u16
            // k-mer frequency vector is no longer stored (#32). In pooled mode
            // this is the whole resident set; in pseudo/per-sample it is one
            // sample's share — multiply by `--sample-jobs` for the peak.
            if params.verbose {
                let nr = raws.len();
                let (mut kmer_b, mut seq_b) = (0usize, 0usize);
                for r in &raws {
                    kmer_b += r.kmer8.as_ref().map_or(0, |v| v.resident_bytes())
                        + r.minimizers.as_ref().map_or(0, |m| m.resident_bytes())
                        + r.kord.as_ref().map_or(0, |v| v.len() * 2);
                    seq_b += r.seq.len() + r.qual.as_ref().map_or(0, |q| q.len());
                }
                let mb = |b: usize| b as f64 / (1024.0 * 1024.0);
                let screen_repr = match params.align.screen_backend {
                    ScreenBackend::Minimizer => format!(
                        "minimizer sketch k={}/w={}",
                        params.align.minimizer_k, params.align.minimizer_w
                    ),
                    ScreenBackend::Kmer if k >= crate::kmers::SPARSE_KMER_MIN => {
                        "sparse #43".to_string()
                    }
                    ScreenBackend::Kmer => "dense".to_string(),
                };
                eprintln!(
                    "[dada] resident Raw footprint: {nr} raws, seq+qual {:.1} MB, \
                     screen vectors {:.1} MB ({:.0} B/raw) [kord k={k}; screen {screen_repr}; \
                     u16 k-mer freq not stored, #32]",
                    mb(seq_b),
                    mb(kmer_b),
                    if nr > 0 {
                        kmer_b as f64 / nr as f64
                    } else {
                        0.0
                    },
                );

                // K-mer complexity / diversity diagnostics (#43). Two signals:
                //  - per-Raw FILL = distinct k-mers / positional max (len-k+1).
                //    ~100% for high-complexity amplicon reads; a low value flags
                //    low-complexity/repetitive sequence content.
                //  - pooled DIVERSITY = union of distinct k-mers across all
                //    uniques vs the 4^k space, plus mean sharing (Σ positional /
                //    union). High sharing = tight homologous amplicon; low
                //    sharing + high occupancy flags diverse/contaminated or
                //    over-amplified data (inflated PCR error/chimera k-mers, not
                //    biological diversity). Cost: one 4^k-bit presence bitmap
                //    (≤8 KB through k8) + a single pass over the screens.
                if nr > 0 && raws.iter().any(|r| r.minimizers.is_some()) {
                    // The minimizer counterpart of the kmer8 fill/diversity
                    // block below. Two signals worth watching:
                    //  - DENSITY = sketch entries / (len - k + 1). Winnowing
                    //    predicts ~2/(w+1); well above that means the sketch is
                    //    not compressing (short reads, or low-complexity content
                    //    defeating the window), which costs memory without
                    //    buying specificity.
                    //  - POOLED SHARING = Σ per-raw entries / distinct
                    //    minimizers across the pool. High sharing means the
                    //    sketch is landing on conserved regions and the screen
                    //    will be permissive; near-1 means every raw is picking
                    //    its own minimizers and the screen will be aggressive.
                    let mut entries_sum = 0usize;
                    let mut positional_sum = 0usize;
                    let mut union: std::collections::HashSet<u64> =
                        std::collections::HashSet::new();
                    let mut unusable = 0usize;
                    let mw = params.align.minimizer_w;
                    let mk = params.align.minimizer_k;
                    for r in &raws {
                        if let Some(m) = &r.minimizers {
                            entries_sum += m.total() as usize;
                            positional_sum += r.len().saturating_sub(mk - 1);
                            union.extend(m.hashes());
                            if m.is_empty() {
                                unusable += 1;
                            }
                        }
                    }
                    let density = entries_sum as f64 / positional_sum.max(1) as f64;
                    eprintln!(
                        "[dada] minimizer sketch: mean {:.0} entries/raw, density {:.3} \
                         (winnowing predicts {:.3} at w={mw}), {unusable} raws unsketchable \
                         (screen bypassed)",
                        entries_sum as f64 / nr as f64,
                        density,
                        2.0 / (mw as f64 + 1.0),
                    );
                    eprintln!(
                        "[dada] minimizer pooled diversity: {} distinct minimizers, \
                         mean sharing {:.1}× across {nr} uniques",
                        union.len(),
                        entries_sum as f64 / union.len().max(1) as f64,
                    );
                }

                if nr > 0 && raws.iter().any(|r| r.kmer8.is_some()) {
                    let nk = crate::kmers::n_kmers(k);
                    let mut bitmap = vec![0u64; nk.div_ceil(64)];
                    let mut distinct_sum = 0usize; // Σ per-raw distinct k-mers
                    let mut positional_sum = 0usize; // Σ (len - k + 1)
                    for r in &raws {
                        if let Some(screen) = &r.kmer8 {
                            distinct_sum += screen.distinct_kmers();
                            positional_sum += r.len().saturating_sub(k - 1);
                            screen.for_each_present_index(|idx| {
                                bitmap[idx >> 6] |= 1u64 << (idx & 63);
                            });
                        }
                    }
                    let union: usize = bitmap.iter().map(|w| w.count_ones() as usize).sum();
                    let fill_pct = 100.0 * distinct_sum as f64 / positional_sum.max(1) as f64;
                    let dense_pct = 100.0 * distinct_sum as f64 / (nr as f64 * nk as f64);
                    eprintln!(
                        "[dada] kmer8 fill: mean {:.0} / {:.0} distinct k-mers/raw \
                         ({fill_pct:.1}% of positional max), {dense_pct:.1}% of dense 4^k [k={k}]",
                        distinct_sum as f64 / nr as f64,
                        positional_sum as f64 / nr as f64,
                    );
                    eprintln!(
                        "[dada] kmer8 pooled diversity: {union} distinct k-mers \
                         ({:.1}% of 4^k space), mean sharing {:.0}× across {nr} uniques",
                        100.0 * union as f64 / nk as f64,
                        positional_sum as f64 / union.max(1) as f64,
                    );
                }
            }
            raws
        }
    };

    // ---- Run core algorithm ----
    let mut b = run_dada(raws, params);

    if params.align.screen_audit {
        // Counters are process-global, so with `--sample-jobs > 1` this is a
        // cumulative snapshot across every sample denoised so far rather than
        // this sample's own figures. Intended for a single pooled run.
        eprintln!("{}", crate::minimizers::audit::summary().report());
        let (hits, tot) = (
            crate::nwalign::GAPLESS_HITS.load(std::sync::atomic::Ordering::Relaxed),
            crate::nwalign::GAPLESS_TOTAL.load(std::sync::atomic::Ordering::Relaxed),
        );
        eprintln!(
            "[screen-audit] gapless shortcut: {hits} / {tot} aligned pairs ({:.2}%)",
            100.0 * hits as f64 / tot.max(1) as f64
        );
    }

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
        .map(|bi| {
            let members = bi.raws.clone();
            let mut member_hammings = Vec::with_capacity(members.len());
            let mut member_lambdas = Vec::with_capacity(members.len());
            let mut member_pvals = Vec::with_capacity(members.len());
            for &raw_idx in &members {
                member_hammings.push(b.raws[raw_idx].comp.hamming);
                member_lambdas.push(b.raws[raw_idx].comp.lambda);
                member_pvals.push(b.raws[raw_idx].p);
            }
            ClusterSummary {
                sequence: bi.seq.clone(),
                reads: bi.reads,
                members,
                member_hammings,
                member_lambdas,
                member_pvals,
                birth_type: bi.birth_type.clone(),
                birth_from: bi.birth_from,
                birth_pval: bi.birth_pval,
                birth_fold: bi.birth_fold,
                birth_e: bi.birth_e,
                birth_hamming: bi.birth_comp.hamming,
            }
        })
        .collect();

    // ---- Aux outputs (R DADA2 parity: $clustering, $birth_subs, $subqual,
    //      $clusterquals) ----
    let aux = if params.aux_outputs {
        Some(compute_aux(&b, params, has_quals))
    } else {
        None
    };

    // The abundance p-value's Bonferroni divisor is the unique count `b.raws`
    // used by the bud scan; it must equal the input count (no raws are dropped,
    // only reassigned between clusters). Locked here so the trace can record it.
    debug_assert_eq!(b.raws.len(), inputs.len(), "nraw divisor != input count");
    let result = DadaResult {
        clusters,
        map,
        pvals,
        nalign: b.nalign,
        nshroud: b.nshroud,
        aux,
        nraw: b.raws.len() as u32,
        omega_a: params.omega_a,
        omega_p: params.omega_p,
        metrics: b.metrics.take(),
    };

    // Reclaim Raws for the caller to pass back on the next iteration.
    Ok((result, std::mem::take(&mut b.raws)))
}

// ---------------------------------------------------------------------------
// Aux-output computation
// ---------------------------------------------------------------------------

/// Compute the R-DADA2-parity outputs (`DadaAux`).
///
/// Re-aligns every Raw against its cluster center (final-subs pass) and each
/// cluster center against its parent (birth-subs pass), both with the k-mer
/// screen disabled (`use_kmers=false, kdist_cutoff=1.0`) so every comparison
/// produces a Sub. Mirrors the `FinalSubsParallel` block in C++ `Rmain.cpp`.
fn compute_aux(b: &B, params: &DadaParams, has_quals: bool) -> DadaAux {
    // Final-subs alignment params: no kmer screen.
    let final_align = AlignParams {
        use_kmers: false,
        kdist_cutoff: 1.0,
        ..params.align
    };
    // Birth-subs alignment params: keep kmer use, no kdist screen.
    let birth_align = AlignParams {
        kdist_cutoff: 1.0,
        ..params.align
    };

    let final_subs = compute_final_subs(b, &final_align);
    let birth_subs = compute_birth_subs(b, &birth_align);

    let cluster_stats_v = cluster_stats(b, &final_subs, &birth_subs, has_quals);
    let maxlen = b.raws.iter().map(|r| r.seq.len()).max().unwrap_or(0);
    let cluster_quality_v = cluster_quality(b, &final_subs, has_quals, maxlen);
    let birth_records = birth_sub_records(&birth_subs, has_quals);
    let ncol = if has_quals { params.err_ncol } else { 1 };
    let transitions = transition_counts(b, &final_subs, has_quals, ncol);

    DadaAux {
        cluster_stats: cluster_stats_v,
        cluster_quality: cluster_quality_v,
        cluster_quality_maxlen: maxlen,
        birth_subs: birth_records,
        transitions,
        transitions_ncol: ncol,
    }
}

/// For each Raw in `b`, align it against its cluster's center and store the
/// resulting `Sub` indexed by `raw.index`. Raws not assigned to a cluster
/// (none in current usage) get `None`. Parallel via Rayon.
fn compute_final_subs(b: &B, align: &AlignParams) -> Vec<Option<Sub>> {
    // (cluster_idx, raw_idx) work items.
    let pairs: Vec<(usize, usize)> = b
        .clusters
        .iter()
        .enumerate()
        .flat_map(|(ci, bi)| bi.raws.iter().map(move |&ri| (ci, ri)))
        .collect();

    let computed: Vec<(u32, Option<Sub>)> = pairs
        .par_iter()
        .map_init(AlignBuffers::new, |buf, &(ci, raw_idx)| {
            let center_idx = match b.clusters[ci].center {
                Some(c) => c,
                None => return (b.raws[raw_idx].index, None),
            };
            let sub = sub_new_with_buf(&b.raws[center_idx], &b.raws[raw_idx], align, buf);
            (b.raws[raw_idx].index, sub)
        })
        .collect();

    let mut out: Vec<Option<Sub>> = (0..b.raws.len()).map(|_| None).collect();
    for (idx, sub) in computed {
        out[idx as usize] = sub;
    }
    out
}

/// For each cluster `i ≥ 1`, align its center against its birth parent's
/// center.  Cluster 0 (and any cluster missing a center) gets `None`.
/// Parallel via Rayon.
fn compute_birth_subs(b: &B, align: &AlignParams) -> Vec<Option<Sub>> {
    (0..b.clusters.len())
        .into_par_iter()
        .map_init(AlignBuffers::new, |buf, ci| {
            if ci == 0 {
                return None;
            }
            let parent_ci = b.clusters[ci].birth_from as usize;
            let parent_center = b.clusters[parent_ci].center?;
            let center = b.clusters[ci].center?;
            sub_new_with_buf(&b.raws[parent_center], &b.raws[center], align, buf)
        })
        .collect()
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
/// Interval, in seconds, between `--verbose` progress lines from the bud loop.
/// Default `30`; `0` disables them. Overridable via `DADA2RS_PROGRESS_SECS`.
///
/// Every other figure `run_dada` prints is a **total**, summed over the whole
/// loop and reported once at the end. That hides any change in the phase mix
/// over the run — and OS-level core usage is observably not flat, ramping
/// before it plateaus, which means the end-of-run means (effective cores, map
/// parallel efficiency) may describe no part of the run. These lines report
/// *deltas since the previous line* so the shape is visible, not just the mean.
///
/// Time-based rather than every N clusters: per-cluster cost is not constant,
/// so a cluster-stride would sample unevenly in exactly the dimension under
/// examination. The cluster index is printed on every line so the two can still
/// be cross-referenced.
fn progress_secs() -> u64 {
    static VALUE: OnceLock<u64> = OnceLock::new();
    *VALUE.get_or_init(|| {
        std::env::var("DADA2RS_PROGRESS_SECS")
            .ok()
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(30)
    })
}

/// Accumulator values at the previous progress line, so the next one can report
/// deltas rather than running totals.
#[derive(Clone, Copy)]
struct ProgressMark {
    at: Duration,
    clusters: usize,
    map: Duration,
    store: Duration,
    busy: Duration,
    shuffle: Duration,
    bud: Duration,
    pupdate: Duration,
    screened: u64,
    aligned: u64,
}

pub fn run_dada(raws: Vec<Raw>, params: &DadaParams) -> B {
    use std::time::{Duration, Instant};
    let mut bb = B::new(raws, params.omega_a, params.omega_p, params.use_quals);

    // Cumulative phase timers. Only `b_compare_parallel` is multithreaded;
    // shuffle/bud/p_update are serial, so their share quantifies the Amdahl
    // serial fraction that caps thread utilization (printed under verbose).
    let (mut t_compare, mut t_shuffle, mut t_bud, mut t_pupdate) = (
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
    );
    // Split of `compare` into the parallel alignment map vs. the serial store,
    // plus summed worker-busy time to derive the map's parallel efficiency.
    let (mut t_cmp_map, mut t_cmp_serial, mut t_cmp_busy) =
        (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    // #143: the rest of `b_compare`, which map+store left unattributed — the
    // serial reduction over the map's results, the free of that result vector,
    // and the pre-map setup. Plus `index_add_cluster`, which sits just outside
    // the compare timer and so lands in run_dada's own residual.
    let (mut t_cmp_agg, mut t_cmp_free, mut t_cmp_setup) =
        (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let mut t_index_add = Duration::ZERO;
    // Store-loop denominators: it scans every raw but pushes only some, so the
    // two rates have to be separated before anything is designed against it.
    let (mut n_cmp_scanned, mut n_cmp_stored) = (0u64, 0u64);
    // Split of the map's worker-busy time into the k-mer screen (paid on every
    // comparison) vs. the DP alignment (paid only by pairs the screen passes),
    // with the matching denominators. This ratio is what #127 exists to
    // measure: it decides whether `b_compare`'s lever is a cheaper screen or a
    // cheaper aligner.
    let (mut t_cmp_screen, mut t_cmp_dp, mut t_cmp_post) =
        (Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let (mut n_cmp_screened, mut n_cmp_aligned) = (0u64, 0u64);
    // Shuffle rescan-redundancy accounting (verbose-only diagnostics).
    let (mut shuf_calls, mut shuf_moves, mut shuf_zero_move_calls) = (0u64, 0u64, 0u64);
    let mut shuf_comps_scanned = 0u64;
    // Split of the scan work: the per-call full build vs the affected-raw
    // reconcile. The two have different access patterns, so this split (with
    // the timings below) is the starting point for sizing any change to the
    // shuffle — see #124.
    let (mut shuf_comps_build, mut shuf_comps_reconcile) = (0u64, 0u64);
    // #139: converge calls that actually ran a build (all of them without the
    // carry; typically one with it).
    let mut shuf_builds = 0u64;
    // b_shuffle_converge invocations (one per bud round) = number of full
    // builds paid; the per-build average is the useful unit here.
    let mut shuf_converge_calls = 0u64;
    // Wall time in each scan phase. Divided by the comp counts above these give
    // ns/comp for the sequential (build) vs scattered (reconcile) access
    // patterns. That ratio is the thing to model before betting on any
    // access-pattern change here: comparison counts alone overstated one such
    // bet ~2x (#87), because a scattered comp costs ~2x a sequential one.
    let (mut t_shuf_build, mut t_shuf_reconcile) =
        (std::time::Duration::ZERO, std::time::Duration::ZERO);
    // The move pass — the third phase, and the one build+reconcile left
    // unaccounted for (#124). Plus the reconcile's internals: how many raws it
    // recomputed and how many of those recomputes actually changed the raw's
    // best cluster, which is what bounds any future reconcile optimization.
    let mut t_shuf_move = std::time::Duration::ZERO;
    let (mut shuf_move_raws, mut shuf_rec_affected, mut shuf_rec_changed) = (0u64, 0u64, 0u64);
    // Reconcile rescan-necessity projection (#136).
    let (mut shuf_rec_rescan, mut shuf_rec_rescan_comps) = (0u64, 0u64);
    let mut shuf_rec_ties: u64 = 0;
    // #139 (reviving #87's projection): what carrying compmax across buds would
    // relocate, split by access pattern.
    let (mut shuf_first_rec_pairs, mut shuf_first_rec_comps) = (0u64, 0u64);
    let mut shuf_rec_pairs: u64 = 0;
    let mut shuf_first_rec_calls: u64 = 0;
    let mut t_rec_collect = std::time::Duration::ZERO;
    let mut t_rec_rescan = std::time::Duration::ZERO;
    // #132 dirty-cluster move-pass diagnostics.
    let (mut shuf_move_unpruned, mut shuf_move_dirty) = (0u64, 0u64);
    let (mut shuf_move_prunable, mut shuf_move_passes) = (0u64, 0u64);
    // b_bud scan-redundancy accounting (verbose-only diagnostics).
    let (mut bud_calls, mut bud_success, mut bud_raws_scanned) = (0u64, 0u64, 0u64);
    // p-update churn: raws whose p was recomputed per round (see issue #85).
    // Only the in-loop p_update rounds count — the pre-loop call reprices the
    // whole partition once and is not part of the per-bud churn a p-ordered
    // structure would face.
    let mut pupd_rounds = 0u64;
    let mut pupd_stats = crate::pval::PUpdateStats::default();

    // One progress record spans the whole bud loop. A DADA2 progress line is
    // assembled as [New Cluster CN:][CNLU:][S per shuffle][the NEXT round's
    // Division fragment], so it is closed at the top of the following round
    // rather than the bottom of its own. Buffering it and emitting one line
    // keeps that text byte-identical while making the write atomic (#172).
    let mut rec = ProgressRecord::for_verbose(params.verbose, params.progress_tag.clone());

    // Initial compare: no k-mer distance screen so that cluster 0 accumulates
    // comparisons for every Raw (required by b_shuffle2).
    let init_params = AlignParams {
        kdist_cutoff: 1.0,
        ..params.align
    };

    let t = Instant::now();
    if params.multithread {
        let ct = b_compare_parallel(
            &mut bb,
            0,
            &params.err_mat,
            params.err_ncol,
            &init_params,
            params.greedy,
            params.measure.attribution(),
        );
        t_cmp_map += ct.map;
        t_cmp_serial += ct.serial;
        t_cmp_busy += ct.busy;
        t_cmp_screen += ct.screen;
        t_cmp_dp += ct.dp;
        t_cmp_post += ct.post;
        n_cmp_screened += ct.screened;
        n_cmp_aligned += ct.aligned;
        t_cmp_agg += ct.agg;
        t_cmp_free += ct.free;
        t_cmp_setup += ct.setup;
        n_cmp_scanned += bb.raws.len() as u64;
        n_cmp_stored += ct.stored;
    } else {
        // `b_compare`'s record drives the `C{i}LU:` progress print, which is a
        // verbose concern, not a measurement one.
        b_compare(
            &mut bb,
            0,
            &params.err_mat,
            params.err_ncol,
            &init_params,
            params.greedy,
            &mut rec,
        );
    }
    t_compare += t.elapsed();
    // Persistent per-raw candidate index for the incremental shuffle driver.
    // Appended once per cluster, in ascending cluster order, right after each
    // compare populates that cluster's comps.
    let mut cand_index: CandIndex = vec![Vec::new(); bb.raws.len()];
    let t = Instant::now();
    index_add_cluster(&mut cand_index, &bb, 0);
    t_index_add += t.elapsed();
    let t = Instant::now();
    b_p_update(
        &mut bb,
        params.greedy,
        params.detect_singletons,
        params.min_fold,
        params.min_hamming,
        params.min_abund,
    );
    t_pupdate += t.elapsed();

    let max_clust = if params.max_clust == 0 {
        bb.raws.len()
    } else {
        params.max_clust
    };

    let mut shuffle_carry = ShuffleCarry::new();

    // Progress-line state (see `progress_secs`). `t_loop` is the bud loop's own
    // clock, so the reported `t=` excludes the setup above and lines up with the
    // phase totals rather than with process wall time.
    let t_loop = Instant::now();
    let progress_every = if params.verbose { progress_secs() } else { 0 };
    let mut mark = ProgressMark {
        at: Duration::ZERO,
        clusters: bb.clusters.len(),
        map: Duration::ZERO,
        store: Duration::ZERO,
        busy: Duration::ZERO,
        shuffle: Duration::ZERO,
        bud: Duration::ZERO,
        pupdate: Duration::ZERO,
        screened: 0,
        aligned: 0,
    };
    let nthreads = rayon::current_num_threads().max(1);

    while bb.clusters.len() < max_clust {
        let t = Instant::now();
        let mut bud_scanned = 0u64;
        let bud = b_bud_incremental(
            &mut bb,
            params.min_fold,
            params.min_hamming,
            params.min_abund,
            &mut rec,
            &mut bud_scanned,
        );
        t_bud += t.elapsed();
        bud_calls += 1;
        bud_raws_scanned += bud_scanned;
        if bud.is_some() {
            bud_success += 1;
        }
        let newi = match bud {
            Some(i) => i,
            None => break,
        };

        // Close the previous round's line (which this round's bud just ended
        // with its Division fragment), then open this one.
        rec.flush();
        rec_print!(rec, "New Cluster C{newi}:");

        let t = Instant::now();
        if params.multithread {
            let ct = b_compare_parallel(
                &mut bb,
                newi,
                &params.err_mat,
                params.err_ncol,
                &params.align,
                params.greedy,
                params.measure.attribution(),
            );
            t_cmp_map += ct.map;
            t_cmp_serial += ct.serial;
            t_cmp_busy += ct.busy;
            t_cmp_screen += ct.screen;
            t_cmp_dp += ct.dp;
            t_cmp_post += ct.post;
            n_cmp_screened += ct.screened;
            n_cmp_aligned += ct.aligned;
            t_cmp_agg += ct.agg;
            t_cmp_free += ct.free;
            t_cmp_setup += ct.setup;
            n_cmp_scanned += bb.raws.len() as u64;
            n_cmp_stored += ct.stored;
        } else {
            b_compare(
                &mut bb,
                newi,
                &params.err_mat,
                params.err_ncol,
                &params.align,
                params.greedy,
                &mut rec,
            );
        }
        t_compare += t.elapsed();
        // Append the new cluster's comps to the persistent candidate index
        // (ascending cluster order preserved: newi is the largest index so far).
        let t = Instant::now();
        index_add_cluster(&mut cand_index, &bb, newi);
        t_index_add += t.elapsed();

        // Shuffle until stable or MAX_SHUFFLE reached — incremental driver.
        // Redundancy accounting: comps_scanned is now the realised scan work
        // (one build + per-iteration recomputes), so comparing it to the serial
        // baseline's counts shows the reduction directly.
        let t = Instant::now();
        // #139: one carry, threaded through every bud round. When the carry is
        // off, `b_shuffle_converge` resets it internally and rebuilds, so this
        // is behaviour-neutral by default.
        let st = b_shuffle_converge(&mut bb, &cand_index, MAX_SHUFFLE, &mut shuffle_carry);
        shuf_converge_calls += 1;
        shuf_calls += st.calls as u64;
        shuf_moves += st.moves as u64;
        shuf_comps_scanned += st.comps_scanned as u64;
        shuf_comps_build += st.comps_build as u64;
        if st.comps_build > 0 {
            shuf_builds += 1;
        }
        shuf_comps_reconcile += st.comps_reconcile as u64;
        t_shuf_build += st.build_time;
        t_shuf_reconcile += st.reconcile_time;
        t_shuf_move += st.move_time;
        shuf_move_raws += st.move_raws_scanned as u64;
        shuf_move_unpruned += st.move_raws_unpruned as u64;
        shuf_move_dirty += st.move_dirty_clusters as u64;
        shuf_move_prunable += st.move_passes_prunable as u64;
        shuf_move_passes += st.move_passes as u64;
        shuf_rec_affected += st.reconcile_affected as u64;
        shuf_rec_rescan += st.reconcile_rescan_raws as u64;
        shuf_rec_ties += st.reconcile_tie_breaks as u64;
        shuf_first_rec_pairs += st.pairs_first_reconcile as u64;
        shuf_rec_pairs += st.pairs_reconcile as u64;
        shuf_first_rec_comps += st.comps_first_reconcile as u64;
        shuf_first_rec_calls += st.first_reconcile_calls as u64;
        shuf_rec_rescan_comps += st.reconcile_comps_rescan as u64;
        t_rec_collect += st.reconcile_collect_time;
        t_rec_rescan += st.reconcile_rescan_time;
        shuf_rec_changed += st.reconcile_changed as u64;
        shuf_zero_move_calls += st.zero_move_calls as u64;
        rec_print!(rec, "{}", "S".repeat(st.calls));
        t_shuffle += t.elapsed();
        if params.verbose && st.calls >= MAX_SHUFFLE {
            eprintln!("Warning: Reached maximum ({MAX_SHUFFLE}) shuffles.");
        }

        let t = Instant::now();
        let repriced = b_p_update(
            &mut bb,
            params.greedy,
            params.detect_singletons,
            params.min_fold,
            params.min_hamming,
            params.min_abund,
        );
        t_pupdate += t.elapsed();
        pupd_rounds += 1;
        pupd_stats += repriced;

        if progress_every > 0 {
            let now = t_loop.elapsed();
            if (now - mark.at).as_secs() >= progress_every {
                let wall = (now - mark.at).as_secs_f64();
                let map = (t_cmp_map - mark.map).as_secs_f64();
                let store = (t_cmp_serial - mark.store).as_secs_f64();
                let busy = (t_cmp_busy - mark.busy).as_secs_f64();
                let shuffle = (t_shuffle - mark.shuffle).as_secs_f64();
                let other = ((t_bud - mark.bud) + (t_pupdate - mark.pupdate)).as_secs_f64();
                let screened = n_cmp_screened - mark.screened;
                let aligned = n_cmp_aligned - mark.aligned;
                // Effective cores over this window: worker-busy time plus the
                // serial time that ran single-threaded, divided by the window.
                // Same construction as the end-of-run figure, so the lines and
                // the total are directly comparable.
                let eff = if wall > 0.0 {
                    (busy + (wall - map).max(0.0)) / wall
                } else {
                    0.0
                };
                let map_eff = if map > 0.0 {
                    busy / (map * nthreads as f64)
                } else {
                    0.0
                };
                let align_frac = if screened > 0 {
                    aligned as f64 / screened as f64 * 100.0
                } else {
                    0.0
                };
                eprintln!(
                    "\n[dada] progress t={:.0}s cluster {} (+{} in {:.0}s): \
                     map={:.1}s store={:.1}s shuffle={:.1}s bud+pupd={:.1}s  \
                     eff cores {:.1}/{}  map eff {:.0}%  align {:.2}%",
                    now.as_secs_f64(),
                    bb.clusters.len(),
                    bb.clusters.len() - mark.clusters,
                    wall,
                    map,
                    store,
                    shuffle,
                    other,
                    eff,
                    nthreads,
                    map_eff * 100.0,
                    align_frac,
                );
                mark = ProgressMark {
                    at: now,
                    clusters: bb.clusters.len(),
                    map: t_cmp_map,
                    store: t_cmp_serial,
                    busy: t_cmp_busy,
                    shuffle: t_shuffle,
                    bud: t_bud,
                    pupdate: t_pupdate,
                    screened: n_cmp_screened,
                    aligned: n_cmp_aligned,
                };
            }
        }
    }

    // Close the final line: the loop exits after a bud that found no
    // division, and that `, No Division.` fragment is the last thing in the
    // record.
    rec.flush();

    // Advisory: did the probe's sample turn out to represent the run?
    //
    // The probe reads the first few clusters, but greedy skipping is NOT
    // stationary -- it grows as clusters accumulate and raws lock, so the early
    // window screens a larger fraction of raws than the run does. On pooled
    // PacBio the window screened 93.7% against a 45.7% run average, overstating
    // the saving 2.05x, always in the direction of building the index.
    //
    // This costs nothing to check: `screened` is already summed for the phase
    // split, and the cluster count is to hand. Re-run the SAME decision function
    // on the rate the run actually exhibited and say so if the verdict flips.
    // There is no threshold here -- the test is whether the answer changes, not
    // whether some ratio crossed a line.
    if let Some(d) = bb.minimizer_index_decision
        && d.forced.is_none()
        && d.clusters > 0
        && !bb.clusters.is_empty()
    {
        let actual = n_cmp_screened / bb.clusters.len() as u64;
        let hindsight = crate::minimizers::decide_from_probe(
            crate::minimizers::ProbeTimings {
                scatter_ns: d.scatter_ns,
                merge_ns: d.merge_ns,
                array_ns: d.array_ns,
            },
            actual as usize,
            d.threads,
            d.entries,
            d.distinct,
            None,
        );
        if hindsight.use_index != d.use_index {
            eprintln!(
                "[dada] warning: the minimizer index choice may have been wrong. The probe \
                 measured {} screened comparisons per cluster over its {} sampled cluster(s), \
                 but the run averaged {} over {} clusters; at that rate the saving is \
                 {:.2} ms against a {:.2} ms scatter, which favours {}. The screen is exact \
                 either way -- this is a speed advisory, not a correctness one. Override with \
                 DADA2RS_MINIMIZER_INDEX={}.",
                d.ncomp,
                d.clusters,
                actual,
                bb.clusters.len(),
                hindsight.saving_ns / 1e6,
                hindsight.scatter_ns / 1e6,
                if hindsight.use_index {
                    "the index"
                } else {
                    "the per-pair merge-join"
                },
                u8::from(hindsight.use_index),
            );
        }
    }

    if params.verbose {
        eprintln!(
            // The leading newline used to terminate the bud loop's dangling
            // partial line. The progress record now ends itself, so keeping it
            // would leave a blank line behind (#172).
            "ALIGN: {} aligns, {} shrouded ({} raw).",
            bb.nalign,
            bb.nshroud,
            bb.raws.len()
        );
        // The index choice is MEASURED on the first cluster, not predicted, and it
        // reverses between workloads (see `minimizers::decide_from_probe`). Report
        // the measurement: a run that declines the index would otherwise look like
        // an unexplained `setup` of 0.00s in the attribution below.
        if let Some(d) = bb.minimizer_index_decision {
            let why = match d.forced {
                Some(true) => " [FORCED on by DADA2RS_MINIMIZER_INDEX=1]",
                Some(false) => " [FORCED off by DADA2RS_MINIMIZER_INDEX=0]",
                None => "",
            };
            if d.forced == Some(false) {
                eprintln!(
                    "[dada] minimizer index: not built (per-pair merge-join) \
                     [FORCED off by DADA2RS_MINIMIZER_INDEX=0]"
                );
            } else {
                eprintln!(
                    "[dada] minimizer index: {} — {} probed cluster(s), mean: scatter {:.2} ms \
                 vs saving {:.2} ms (= ({:.0} - {:.1} ns/pair) x {} screened comps / \
                 {} threads); {} distinct minimizers, {} postings, mean posting {:.0}{}",
                    if d.use_index {
                        "USED (scatter per cluster)"
                    } else {
                        "declined (per-pair merge-join)"
                    },
                    d.clusters,
                    d.scatter_ns / 1e6,
                    d.saving_ns / 1e6,
                    d.merge_ns,
                    d.array_ns,
                    d.ncomp,
                    d.threads,
                    d.distinct,
                    d.entries,
                    d.sharing,
                    why,
                );
            }
        }
        eprintln!(
            "[dada] phase times (serial except compare-map): compare={:.2}s (map={:.2}s parallel, store={:.2}s serial)  shuffle={:.2}s  bud={:.2}s  p_update={:.2}s",
            t_compare.as_secs_f64(),
            t_cmp_map.as_secs_f64(),
            t_cmp_serial.as_secs_f64(),
            t_shuffle.as_secs_f64(),
            t_bud.as_secs_f64(),
            t_pupdate.as_secs_f64(),
        );
        // The optimisation attribution -- compare attribution and split, map
        // parallel efficiency, the shuffle phases and scan split, bud
        // redundancy, p-update churn, and the #132/#136/#139 projections --
        // is no longer printed. It was ~130 lines of ns/comp tables written
        // for specific issues, which buried the handful of lines that answer
        // "how big is this run and is it configured sanely" (#162).
        //
        // Every one of those quantities is in `--metrics-json`, read from the
        // same accumulators, so nothing was lost in the move -- verified by
        // `dev/check_metrics_superset.py` against a 30-sample ITS2 sweep and
        // a pooled run of the same data.
        eprintln!(
            "[dada] phase attribution, compare/shuffle splits and optimisation \
             counters: pass --metrics-json <path> (add --metrics-attribution \
             for the per-comparison timings)."
        );
    }

    // ---- Machine-readable metrics (#162) -------------------------------
    //
    // Built from the SAME accumulators the prose above printed, so the two
    // cannot drift. Gated on `measure`, not `verbose`: `--metrics-json` asks
    // for numbers without asking for output.
    if params.measure.enabled() {
        let nthreads = rayon::current_num_threads().max(1);

        // Recomputed here rather than threaded out of the Raw build: it is one
        // O(nraw) walk over the same raws, so the numbers are identical, and
        // the build runs in a scope that does not reach this one.
        let footprint = {
            let (mut kmer_b, mut seq_b) = (0usize, 0usize);
            for r in &bb.raws {
                kmer_b += r.kmer8.as_ref().map_or(0, |v| v.resident_bytes())
                    + r.minimizers.as_ref().map_or(0, |m| m.resident_bytes())
                    + r.kord.as_ref().map_or(0, |v| v.len() * 2);
                seq_b += r.seq.len() + r.qual.as_ref().map_or(0, |q| q.len());
            }
            let nr = bb.raws.len();
            crate::metrics::FootprintMetrics {
                nraw: nr,
                seq_qual_bytes: seq_b as u64,
                screen_vector_bytes: kmer_b as u64,
                screen_bytes_per_raw: (nr > 0).then(|| kmer_b as f64 / nr as f64),
                screen_repr: match params.align.screen_backend {
                    ScreenBackend::Minimizer => format!(
                        "minimizer sketch k={}/w={}",
                        params.align.minimizer_k, params.align.minimizer_w
                    ),
                    ScreenBackend::Kmer
                        if params.align.kmer_size >= crate::kmers::SPARSE_KMER_MIN =>
                    {
                        "sparse #43".to_string()
                    }
                    ScreenBackend::Kmer => "dense".to_string(),
                },
            }
        };

        let counters = RawCounters {
            t_compare,
            t_shuffle,
            t_bud,
            t_pupdate,
            t_index_add,
            t_loop: t_loop.elapsed(),
            t_cmp_map,
            t_cmp_serial,
            t_cmp_busy,
            t_cmp_agg,
            t_cmp_free,
            t_cmp_setup,
            t_cmp_screen,
            t_cmp_dp,
            t_cmp_post,
            n_cmp_scanned,
            n_cmp_stored,
            n_cmp_screened,
            n_cmp_aligned,
            t_shuf_build,
            t_shuf_reconcile,
            t_shuf_move,
            shuf_comps_scanned,
            shuf_comps_build,
            shuf_comps_reconcile,
            shuf_calls,
            shuf_converge_calls,
            shuf_builds,
            shuf_moves,
            shuf_zero_move_calls,
            shuf_move_raws,
            t_rec_collect,
            t_rec_rescan,
            shuf_rec_affected,
            shuf_rec_changed,
            shuf_rec_rescan,
            shuf_rec_rescan_comps,
            shuf_rec_ties,
            shuf_rec_pairs,
            shuf_move_unpruned,
            shuf_move_dirty,
            shuf_move_prunable,
            shuf_move_passes,
            shuf_first_rec_pairs,
            shuf_first_rec_comps,
            shuf_first_rec_calls,
            bud_calls,
            bud_success,
            bud_raws_scanned,
            pupd_rounds,
            pupd: crate::metrics::PUpdateMetrics {
                repriced: pupd_stats.repriced,
                exit_singleton: pupd_stats.exit_singleton,
                exit_center: pupd_stats.exit_center,
                exit_zero_lambda: pupd_stats.exit_zero_lambda,
                full_calc: pupd_stats.full_calc,
                lock_scanned: pupd_stats.lock_scanned,
                dirty_clusters: pupd_stats.dirty_clusters,
                clean_clusters: pupd_stats.clean_clusters,
                // Filled from `pupd_rounds` in `finish`.
                rounds: 0,
            },
            footprint: Some(footprint),
        };

        let shape = RunShape {
            nraw: bb.raws.len(),
            reads: bb.reads,
            nclusters: bb.clusters.len(),
            nalign: bb.nalign,
            nshroud: bb.nshroud,
            threads: nthreads,
            multithread: params.multithread,
            // The bare backend name, not the `backend_repr` log line: this is a
            // field, not a sentence.
            align_backend: format!("{:?}", params.align.backend).to_lowercase(),
            screen_backend: format!("{:?}", params.align.screen_backend).to_lowercase(),
            kmer_size: params.align.kmer_size,
            kdist_cutoff: params.align.kdist_cutoff,
            band: params.align.band,
            minimizer_k: params.align.minimizer_k,
            minimizer_w: params.align.minimizer_w,
            gates: crate::gates::report(),
        };

        let mut m: RunMetrics = counters.finish(shape, nthreads, params.measure);
        m.index = bb.minimizer_index_decision.map(|d| {
            crate::metrics::IndexDecision {
                use_index: d.use_index,
                forced: d.forced,
                probed_clusters: d.clusters,
                scatter_ns: d.scatter_ns,
                saving_ns: d.saving_ns,
                merge_ns: d.merge_ns,
                array_ns: d.array_ns,
                screened_comps: d.ncomp as u64,
                threads: d.threads,
                distinct_minimizers: d.distinct as u64,
                postings: d.entries as u64,
                mean_posting: d.sharing,
                // Filled by the caller's hindsight check above when it fired;
                // false here means the probe and the run agreed.
                hindsight_disagrees: false,
            }
        });
        bb.metrics = Some(m);
    }

    bb
}
