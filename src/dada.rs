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
use crate::kmers::{KMER_SIZE_MAX, KMER_SIZE_MIN, raw_assign_kmers};
use crate::misc::nt_encode;
use crate::nwalign::{AlignBuffers, AlignParams, sub_new_with_buf};
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
    pub nalign: u32,
    /// Comparisons screened out by k-mer distance.
    pub nshroud: u32,
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
                for raw in &mut raws {
                    raw_assign_kmers(raw, k);
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
                        + r.kord.as_ref().map_or(0, |v| v.len() * 2);
                    seq_b += r.seq.len() + r.qual.as_ref().map_or(0, |q| q.len());
                }
                let mb = |b: usize| b as f64 / (1024.0 * 1024.0);
                let screen_repr = if k >= crate::kmers::SPARSE_KMER_MIN {
                    "sparse #43"
                } else {
                    "dense"
                };
                eprintln!(
                    "[dada] resident Raw footprint: {nr} raws, seq+qual {:.1} MB, \
                     k-mer vectors {:.1} MB ({:.0} B/raw) [k={k}; kmer8 {screen_repr}; u16 k-mer freq not stored, #32]",
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
/// Mirrors `cluster::reconcile_full` for the report, which needs to know which
/// path produced the numbers.
fn reconcile_full_env() -> bool {
    std::env::var_os("DADA2RS_RECONCILE_FULL").is_some()
}

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
    let mut shuf_nraw = 0usize;
    // b_bud scan-redundancy accounting (verbose-only diagnostics).
    let (mut bud_calls, mut bud_success, mut bud_raws_scanned) = (0u64, 0u64, 0u64);
    // p-update churn: raws whose p was recomputed per round (see issue #85).
    // Only the in-loop p_update rounds count — the pre-loop call reprices the
    // whole partition once and is not part of the per-bud churn a p-ordered
    // structure would face.
    let (mut pupd_rounds, mut pupd_raws_repriced) = (0u64, 0u64);

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
            params.verbose,
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
            0,
            &params.err_mat,
            params.err_ncol,
            &init_params,
            params.greedy,
            params.verbose,
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

    while bb.clusters.len() < max_clust {
        let t = Instant::now();
        let mut bud_scanned = 0u64;
        let bud = b_bud_incremental(
            &mut bb,
            params.min_fold,
            params.min_hamming,
            params.min_abund,
            params.verbose,
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

        if params.verbose {
            eprint!("\nNew Cluster C{newi}:");
        }

        let t = Instant::now();
        if params.multithread {
            let ct = b_compare_parallel(
                &mut bb,
                newi,
                &params.err_mat,
                params.err_ncol,
                &params.align,
                params.greedy,
                params.verbose,
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
                params.verbose,
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
        shuf_nraw = st.nraw;
        shuf_zero_move_calls += st.zero_move_calls as u64;
        if params.verbose {
            eprint!("{}", "S".repeat(st.calls));
        }
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
        pupd_raws_repriced += repriced;
    }

    if params.verbose {
        eprintln!(
            "\nALIGN: {} aligns, {} shrouded ({} raw).",
            bb.nalign,
            bb.nshroud,
            bb.raws.len()
        );
        // Parallel efficiency of the map: worker-busy time / (map wall ×
        // threads). Near 1.0 → threads compute the whole map wall (low OS-level
        // utilization then implies memory-bandwidth stalls); well below 1.0 →
        // threads idle inside the parallel region (tail load-imbalance).
        let nthreads = rayon::current_num_threads().max(1);
        let map_eff = if t_cmp_map.as_secs_f64() > 0.0 {
            t_cmp_busy.as_secs_f64() / (t_cmp_map.as_secs_f64() * nthreads as f64)
        } else {
            0.0
        };
        eprintln!(
            "[dada] phase times (serial except compare-map): compare={:.2}s (map={:.2}s parallel, store={:.2}s serial)  shuffle={:.2}s  bud={:.2}s  p_update={:.2}s",
            t_compare.as_secs_f64(),
            t_cmp_map.as_secs_f64(),
            t_cmp_serial.as_secs_f64(),
            t_shuffle.as_secs_f64(),
            t_bud.as_secs_f64(),
            t_pupdate.as_secs_f64(),
        );
        // #143: attribute the rest of `b_compare`. Before this, map+store
        // covered only part of the compare timer and the remainder was
        // reported as nothing at all. Each line below is serial, so it caps
        // node occupancy exactly the way the shuffle does.
        {
            let cmp = t_compare.as_secs_f64();
            let pct = |x: f64| if cmp > 0.0 { 100.0 * x / cmp } else { 0.0 };
            let per = |x: f64, n: u64| {
                if n > 0 { x * 1e9 / n as f64 } else { 0.0 }
            };
            let residual = cmp
                - (t_cmp_map.as_secs_f64()
                    + t_cmp_serial.as_secs_f64()
                    + t_cmp_agg.as_secs_f64()
                    + t_cmp_free.as_secs_f64()
                    + t_cmp_setup.as_secs_f64());
            eprintln!(
                "[dada] compare attribution (of {:.2}s over {} raw-visits, {} stored):",
                cmp, n_cmp_scanned, n_cmp_stored,
            );
            eprintln!(
                "[dada]   map          {:8.2}s ({:4.1}%)  parallel",
                t_cmp_map.as_secs_f64(),
                pct(t_cmp_map.as_secs_f64()),
            );
            eprintln!(
                "[dada]   reduction    {:8.2}s ({:4.1}%)  serial  {:6.1} ns/raw  (summed costs + denominators)",
                t_cmp_agg.as_secs_f64(),
                pct(t_cmp_agg.as_secs_f64()),
                per(t_cmp_agg.as_secs_f64(), n_cmp_scanned),
            );
            eprintln!(
                "[dada]   store        {:8.2}s ({:4.1}%)  serial  {:6.1} ns/raw  {:6.1} ns/stored",
                t_cmp_serial.as_secs_f64(),
                pct(t_cmp_serial.as_secs_f64()),
                per(t_cmp_serial.as_secs_f64(), n_cmp_scanned),
                per(t_cmp_serial.as_secs_f64(), n_cmp_stored),
            );
            eprintln!(
                "[dada]   free         {:8.2}s ({:4.1}%)  serial  {:6.1} ns/raw  (result vector)",
                t_cmp_free.as_secs_f64(),
                pct(t_cmp_free.as_secs_f64()),
                per(t_cmp_free.as_secs_f64(), n_cmp_scanned),
            );
            eprintln!(
                "[dada]   setup        {:8.2}s ({:4.1}%)  serial",
                t_cmp_setup.as_secs_f64(),
                pct(t_cmp_setup.as_secs_f64()),
            );
            eprintln!(
                "[dada]   unattributed {:8.2}s ({:4.1}%)",
                residual,
                pct(residual),
            );
            eprintln!(
                "[dada]   (outside compare) index_add_cluster {:.2}s  {:6.1} ns/raw",
                t_index_add.as_secs_f64(),
                per(t_index_add.as_secs_f64(), n_cmp_scanned),
            );
        }
        eprintln!(
            "[dada] map parallel efficiency: {:.0}% (busy={:.0}s / map={:.0}s × {} threads)",
            100.0 * map_eff,
            t_cmp_busy.as_secs_f64(),
            t_cmp_map.as_secs_f64(),
            nthreads,
        );
        // Screen-vs-align split of the map's worker-busy time (#127).
        //
        // Read it against the shroud rate: the screen's share is what the
        // aligner never gets to avoid. A screen-dominated profile points at
        // replacing the prefilter (an index that skips comparisons outright);
        // an align-dominated one points back at the aligner, which is already
        // well-explored. `other` is compute_lambda plus per-item overhead,
        // including this instrumentation's own timer cost.
        {
            let busy = t_cmp_busy.as_secs_f64();
            let screen = t_cmp_screen.as_secs_f64();
            let dp = t_cmp_dp.as_secs_f64();
            let post = t_cmp_post.as_secs_f64();
            let align = dp + post;
            let pct = |x: f64| if busy > 0.0 { 100.0 * x / busy } else { 0.0 };
            let per = |x: f64, n: u64| {
                if n > 0 { x * 1e9 / n as f64 } else { 0.0 }
            };
            eprintln!(
                "[dada] compare split (of {:.2}s busy over {} comparisons):",
                busy, n_cmp_screened,
            );
            eprintln!(
                "[dada]   kmer screen  {:8.2}s ({:4.1}%)  {:>13} comps  {:6.0} ns/comp  (every comparison)",
                screen,
                pct(screen),
                n_cmp_screened,
                per(screen, n_cmp_screened),
            );
            eprintln!(
                "[dada]   align total  {:8.2}s ({:4.1}%)  {:>13} comps  {:6.0} ns/comp  ({:.1}% passed the screen)",
                align,
                pct(align),
                n_cmp_aligned,
                per(align, n_cmp_aligned),
                if n_cmp_screened > 0 {
                    100.0 * n_cmp_aligned as f64 / n_cmp_screened as f64
                } else {
                    0.0
                },
            );
            // DP-vs-post split (#127): a high align total is ambiguous between a
            // slow kernel and expensive post-processing, and the two point at
            // completely different fixes.
            eprintln!(
                "[dada]     dp kernel  {:8.2}s ({:4.1}%)  {:>13} comps  {:6.0} ns/comp",
                dp,
                pct(dp),
                n_cmp_aligned,
                per(dp, n_cmp_aligned),
            );
            eprintln!(
                "[dada]     al2subs    {:8.2}s ({:4.1}%)  {:>13} comps  {:6.0} ns/comp  (+ qual mapping)",
                post,
                pct(post),
                n_cmp_aligned,
                per(post, n_cmp_aligned),
            );
            eprintln!(
                "[dada]   other        {:8.2}s ({:4.1}%)  (compute_lambda + per-item overhead)",
                busy - screen - align,
                pct(busy - screen - align),
            );
        }
        // Shuffle rescan redundancy: how much of the full-rescan work each
        // b_shuffle2 call actually translated into moves. High scanned/move and
        // a high zero-move-call fraction bound the headroom for an incremental
        // best-cluster tracker (the "reduce work" lever in docs/results.md).
        let scanned_per_move = if shuf_moves > 0 {
            shuf_comps_scanned as f64 / shuf_moves as f64
        } else {
            f64::INFINITY
        };
        let zero_pct = if shuf_calls > 0 {
            100.0 * shuf_zero_move_calls as f64 / shuf_calls as f64
        } else {
            0.0
        };
        eprintln!(
            "[dada] shuffle redundancy: {} calls, {} moves, {} comps scanned ({:.0} scanned/move), {} zero-move calls ({:.0}%)",
            shuf_calls,
            shuf_moves,
            shuf_comps_scanned,
            scanned_per_move,
            shuf_zero_move_calls,
            zero_pct,
        );
        // Where that scan work actually goes. Zero-move iterations break before
        // the reconcile, so they scan nothing — the zero-move percentage above
        // is NOT a measure of wasted scanning. The removable work is the full
        // build, re-paid once per bud round; the reconcile is only incurred by
        // raws whose best cluster may actually have changed. `comps/build` also
        // tracks cluster structure, which is what made #87's payoff flip sign
        // between amplicons.
        let build_pct = if shuf_comps_scanned > 0 {
            100.0 * shuf_comps_build as f64 / shuf_comps_scanned as f64
        } else {
            0.0
        };
        // #139: with the carry on, most converge calls do no build at all, so
        // dividing by the call count would report a per-build cost for builds
        // that never happened. Divide by the calls that actually built.
        let per_build = if shuf_builds > 0 {
            shuf_comps_build as f64 / shuf_builds as f64
        } else {
            0.0
        };
        eprintln!(
            "[dada] shuffle scan split: build={} ({:.0}% of scanned) over {} builds ({:.0} comps/build), reconcile={} ({:.0}%)",
            shuf_comps_build,
            build_pct,
            shuf_builds,
            per_build,
            shuf_comps_reconcile,
            100.0 - build_pct,
        );
        // Time split, and the ns/comp each access pattern actually costs. The
        // build is a contiguous cluster-major walk; the reconcile walks the
        // raw-major inverted index (scattered). Measured at ~2x apart, so any
        // scheme that shifts work between the two must be priced in ns/comp
        // rather than comparisons — that gap is what made #87 a much smaller
        // win than its comp counts implied.
        let ns_per = |d: std::time::Duration, n: u64| {
            if n > 0 {
                d.as_secs_f64() * 1e9 / n as f64
            } else {
                0.0
            }
        };
        // Full phase accounting for the shuffle (#124). build + reconcile
        // alone left 15-19% of shuffle time unexplained; the move pass is that
        // remainder. `other` should now be near zero — if it is not, there is
        // still an unmeasured phase and any optimization here is being sized
        // against an incomplete denominator.
        //
        // The three columns to compare across platforms: each phase's share of
        // shuffle, its ns per unit of work, and its redundancy (work done per
        // outcome produced). A phase with a high share AND high redundancy is
        // the only kind worth attacking.
        let shuf_s = t_shuffle.as_secs_f64();
        let other = shuf_s
            - t_shuf_build.as_secs_f64()
            - t_shuf_reconcile.as_secs_f64()
            - t_shuf_move.as_secs_f64();
        let pct = |d: f64| {
            if shuf_s > 0.0 {
                100.0 * d / shuf_s
            } else {
                0.0
            }
        };
        eprintln!(
            "[dada] shuffle phases ({:.2}s total over {} converge calls, nraw={}, nclusters={}):",
            shuf_s,
            shuf_converge_calls,
            shuf_nraw,
            bb.clusters.len(),
        );
        eprintln!(
            "[dada]   build     {:>8.2}s ({:>4.1}%)  {:>13} comps  {:>6.2} ns/comp",
            t_shuf_build.as_secs_f64(),
            pct(t_shuf_build.as_secs_f64()),
            shuf_comps_build,
            ns_per(t_shuf_build, shuf_comps_build),
        );
        eprintln!(
            "[dada]   reconcile {:>8.2}s ({:>4.1}%)  {:>13} comps  {:>6.2} ns/comp  over {} raws, {} ({:.1}%) actually changed cluster",
            t_shuf_reconcile.as_secs_f64(),
            pct(t_shuf_reconcile.as_secs_f64()),
            shuf_comps_reconcile,
            ns_per(t_shuf_reconcile, shuf_comps_reconcile),
            shuf_rec_affected,
            shuf_rec_changed,
            if shuf_rec_affected > 0 {
                100.0 * shuf_rec_changed as f64 / shuf_rec_affected as f64
            } else {
                0.0
            },
        );
        eprintln!(
            "[dada]   move      {:>8.2}s ({:>4.1}%)  {:>13} raws   {:>6.2} ns/raw   for {} moves ({:.0} raws scanned/move)",
            t_shuf_move.as_secs_f64(),
            pct(t_shuf_move.as_secs_f64()),
            shuf_move_raws,
            ns_per(t_shuf_move, shuf_move_raws),
            shuf_moves,
            if shuf_moves > 0 {
                shuf_move_raws as f64 / shuf_moves as f64
            } else {
                f64::INFINITY
            },
        );
        // Reconcile rescan-necessity (#136, measurement only).
        //
        // #124 called the reconcile's redundancy unreachable: 99.97% of
        // recomputes return the value already in `compmax`, and finding which
        // raws changed costs the scattered access being avoided. That rules out
        // *skipping* the touch, not doing it more cheaply. A raw only needs its
        // full candidate rescan when its current best cluster shrank; otherwise
        // that cluster still beats every unchanged candidate and only the
        // changed ones need testing, which a sequential walk of the changed
        // clusters' comps already reaches.
        //
        // This reports the fraction that a cheaper scheme could NOT avoid. If
        // it is most of them, the redundancy really is unreachable and #124's
        // verdict stands as written.
        if shuf_rec_affected > 0 {
            // Mode-aware: the two paths make `affected` mean different things.
            // Under the baseline full rescan it is every raw in a changed
            // cluster, and this counter *projects* how many genuinely need
            // rescanning. Under the incremental path only those raws are
            // collected at all, so a percentage there would measure the counter
            // against itself and always read 100%.
            if reconcile_full_env() {
                eprintln!(
                    "[dada]   reconcile rescan-necessity (#136): {} of {} affected raws must rescan \
                     ({:.1}%), {} of {} comps ({:.1}%); collect {:.2}s + rescan {:.2}s, \
                     so ~{:.2}s of {:.2}s reconcile is avoidable",
                    shuf_rec_rescan,
                    shuf_rec_affected,
                    100.0 * shuf_rec_rescan as f64 / shuf_rec_affected as f64,
                    shuf_rec_rescan_comps,
                    shuf_comps_reconcile,
                    if shuf_comps_reconcile > 0 {
                        100.0 * shuf_rec_rescan_comps as f64 / shuf_comps_reconcile as f64
                    } else {
                        0.0
                    },
                    t_rec_collect.as_secs_f64(),
                    t_rec_rescan.as_secs_f64(),
                    // The rescan half scaled by the share of comps a cheaper scheme
                    // would not have to walk. Optimistic: it credits the survivors
                    // with the same per-comp rate, and the collect half is unchanged.
                    t_rec_rescan.as_secs_f64()
                        * if shuf_comps_reconcile > 0 {
                            1.0 - shuf_rec_rescan_comps as f64 / shuf_comps_reconcile as f64
                        } else {
                            0.0
                        },
                    t_shuf_reconcile.as_secs_f64(),
                );
            } else {
                eprintln!(
                    "[dada]   reconcile incremental (#136): {} raws fully rescanned over {} \
                     comps; collect {:.2}s + rescan {:.2}s of {:.2}s reconcile; \
                     {} exact-tie tie-breaks",
                    shuf_rec_rescan,
                    shuf_rec_rescan_comps,
                    t_rec_collect.as_secs_f64(),
                    t_rec_rescan.as_secs_f64(),
                    t_shuf_reconcile.as_secs_f64(),
                    shuf_rec_ties,
                );
            }
        }
        // #139: would carrying `compmax` across buds pay? Reviving #87's
        // projection, with its cost model re-derived — #87 priced the relocated
        // work at the full-rescan reconcile's scattered rate, and #136 replaced
        // that with an incremental update.
        //
        // Skipping the per-bud build saves `comps_build x build_ns`. It costs
        // the first reconcile instead, which post-#136 is a sequential
        // cluster-major walk plus a much smaller scattered rescan. Both are
        // measured here rather than assumed, and priced at this run's own rates.
        //
        // Caveat: the pair and comp counts are specific to the first reconcile,
        // but the ns rates are averaged over all of them. A post-bud reconcile
        // touches more clusters than a mid-loop one, so if its per-unit cost
        // differs the estimate drifts. Directionally it is the volume, not the
        // rate, that separates the regimes.
        //
        // With the carry on (the default) the projection is moot — the build is
        // already skipped, and `shuf_comps_build` counts only the single initial
        // build that seeds the map. Report what was realised instead of
        // projecting what could be, so the two arms are never confused.
        if crate::cluster::shuffle_carry() {
            eprintln!(
                "[dada]   #87 carry (#139): ACTIVE -- {} of {} converge calls built \
                 ({} comps, {:.1}s); {} entered on a carried map",
                shuf_builds,
                shuf_converge_calls,
                shuf_comps_build,
                t_shuf_build.as_secs_f64(),
                shuf_converge_calls.saturating_sub(shuf_builds),
            );
        } else if shuf_first_rec_calls > 0 && shuf_comps_build > 0 {
            let build_ns = t_shuf_build.as_secs_f64() * 1e9 / shuf_comps_build as f64;
            let collect_ns = if shuf_rec_pairs > 0 {
                t_rec_collect.as_secs_f64() * 1e9 / shuf_rec_pairs as f64
            } else {
                0.0
            };
            let rescan_ns = if shuf_comps_reconcile > 0 {
                t_rec_rescan.as_secs_f64() * 1e9 / shuf_comps_reconcile as f64
            } else {
                0.0
            };
            let saved = shuf_comps_build as f64 * build_ns;
            let cost =
                shuf_first_rec_pairs as f64 * collect_ns + shuf_first_rec_comps as f64 * rescan_ns;
            eprintln!(
                "[dada]   #87 projection (#139): skipping the per-bud build saves {:.1}s \
                 ({} comps at {:.2} ns), costs {:.1}s relocated ({} pairs at {:.2} ns + \
                 {} comps at {:.2} ns) over {} calls => {} by {:.2}x",
                saved / 1e9,
                shuf_comps_build,
                build_ns,
                cost / 1e9,
                shuf_first_rec_pairs,
                collect_ns,
                shuf_first_rec_comps,
                rescan_ns,
                shuf_first_rec_calls,
                if cost < saved { "WIN" } else { "LOSS" },
                if cost > 0.0 {
                    saved / cost
                } else {
                    f64::INFINITY
                },
            );
        }
        // #132: how much the dirty-cluster pruning actually bought. The prune
        // is workload-dependent (64-67% MiSeq, 57% PacBio), so it is reported
        // rather than assumed — an erosion should be visible here, not inferred
        // from wall time.
        if shuf_move_passes > 0 {
            let pruned = shuf_move_unpruned.saturating_sub(shuf_move_raws);
            eprintln!(
                "[dada]   move pruning (#132): {} of {} raws skipped ({:.1}%), \
                 {} of {} passes pruned ({:.0}%), mean {:.1} dirty clusters/pass of {}",
                pruned,
                shuf_move_unpruned,
                if shuf_move_unpruned > 0 {
                    100.0 * pruned as f64 / shuf_move_unpruned as f64
                } else {
                    0.0
                },
                shuf_move_prunable,
                shuf_move_passes,
                100.0 * shuf_move_prunable as f64 / shuf_move_passes as f64,
                if shuf_move_prunable > 0 {
                    shuf_move_dirty as f64 / shuf_move_prunable as f64
                } else {
                    0.0
                },
                bb.clusters.len(),
            );
        }
        eprintln!(
            "[dada]   other     {:>8.2}s ({:>4.1}%)  (loop bookkeeping; large = an unmeasured phase)",
            other,
            pct(other),
        );
        eprintln!(
            "[dada] shuffle scan time: build={:.2}s ({:.2} ns/comp)  reconcile={:.2}s ({:.2} ns/comp)  build={:.0}% of shuffle time",
            t_shuf_build.as_secs_f64(),
            ns_per(t_shuf_build, shuf_comps_build),
            t_shuf_reconcile.as_secs_f64(),
            ns_per(t_shuf_reconcile, shuf_comps_reconcile),
            if t_shuffle.as_secs_f64() > 0.0 {
                100.0 * t_shuf_build.as_secs_f64() / t_shuffle.as_secs_f64()
            } else {
                0.0
            },
        );
        // b_bud combine cost: with the incremental candidate cache (#85) each
        // bud combines per-cluster minima in O(nclusters) instead of rescanning
        // every raw. This reports the combine volume — compare to the historical
        // ~nraw scanned/bud to see the reduction.
        let combine_per_bud = if bud_calls > 0 {
            bud_raws_scanned as f64 / bud_calls as f64
        } else {
            f64::INFINITY
        };
        eprintln!(
            "[dada] bud redundancy: {} calls ({} budded), {} clusters combined ({:.0} combined/bud, incremental cache)",
            bud_calls, bud_success, bud_raws_scanned, combine_per_bud,
        );
        // p-update churn between bud calls (issue #85 gate): raws whose p is
        // recomputed per in-loop round. A p-ordered incremental budding
        // structure must re-key each of these. Compare repriced/round to nraw:
        // if it approaches the scanned/bud figure above, re-keying costs about
        // as much as the scan it would replace — i.e. no obvious net win.
        let repriced_per_round = if pupd_rounds > 0 {
            pupd_raws_repriced as f64 / pupd_rounds as f64
        } else {
            0.0
        };
        eprintln!(
            "[dada] p-update churn: {} rounds, {} raws repriced ({:.0} repriced/round, nraw={})",
            pupd_rounds,
            pupd_raws_repriced,
            repriced_per_round,
            bb.raws.len(),
        );
    }

    bb
}
