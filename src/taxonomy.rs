//! Naive Bayes k-mer taxonomic classifier.
//!
//! Ports `taxonomy.cpp`, excluding all R/Rcpp and RcppParallel wrappers.
//!
//! ## Algorithm overview
//! 1. Build a log-probability table: for each genus × k-mer, log P(kmer | genus),
//!    smoothed by a cross-genus kmer prior (Laplace-like).
//! 2. For each query sequence, compute the sorted k-mer index array and find the
//!    genus that maximises the sum of log probabilities (ties broken uniformly at
//!    random — reservoir sampling).
//! 3. Bootstrap confidence: resample 1/8 of the query's k-mers `NBOOT` times and
//!    record how often each bootstrap hit agrees with the full assignment at every
//!    taxonomic level.
//!
//! Sequences shorter than 50 bp receive no assignment (`None`).
//!
//! ## Differences from the C++ original
//! - Random numbers are generated on-the-fly via `rand::thread_rng()` instead of
//!   being pre-allocated by `Rcpp::runif`.
//! - `RcppParallel::parallelFor` is replaced by Rayon `par_iter`.
//! - All indexing is 0-based; callers should add 1 if they need R-style output.
//! - `Rcpp::checkUserInterrupt()` is removed (no R event loop).

use rayon::prelude::*;

/// Number of bootstrap replicates.  Matches C++ `NBOOT`.
pub const NBOOT: usize = 100;

/// K-mer size used for classification.  Matches the hard-coded `k=8` in C++.
pub const TAX_K: usize = 8;

/// Minimum query sequence length for an assignment to be attempted.
const MIN_SEQ_LEN: usize = 50;

// ---------------------------------------------------------------------------
// Output type
// ---------------------------------------------------------------------------

/// Output of [`assign_taxonomy`].
pub struct TaxonomyResult {
    /// For each query sequence (in input order), the 0-indexed best-genus index,
    /// or `None` if the sequence was too short to classify.
    pub assignments: Vec<Option<usize>>,
    /// Bootstrap agreement counts with shape `[nseq][nlevel]`.
    /// `boot_counts[i][l]` is the number of `NBOOT` bootstrap replicates whose
    /// genus assignment matched the full assignment at taxonomic level `l`.
    pub boot_counts: Vec<Vec<u32>>,
    /// Per-bootstrap genus assignments with shape `[nseq][NBOOT]`.
    /// `boot_taxa[i][b]` is the 0-indexed genus chosen in replicate `b` for
    /// sequence `i`, or `None` if the sequence was too short.
    pub boot_taxa: Vec<Vec<Option<usize>>>,
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Compute the k-mer index for an ASCII window of length `k`.
/// Returns `None` if any base is not A/C/G/T (case-sensitive).
/// Equivalent to C++ `tax_kmer`.
fn tax_kmer(window: &[u8], k: usize) -> Option<usize> {
    let mut kmer = 0usize;
    for &b in &window[..k] {
        let nti = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        kmer = 4 * kmer + nti;
    }
    Some(kmer)
}

/// Fill a presence/absence bit-vector (`kvec`) for all valid k-mers in `seq`.
/// `kvec` must have length `4^k`.
/// Equivalent to C++ `tax_kvec`.
fn tax_kvec(seq: &[u8], k: usize, kvec: &mut [u8]) {
    kvec.fill(0);
    let klen = seq.len().saturating_sub(k - 1);
    for i in 0..klen {
        if let Some(km) = tax_kmer(&seq[i..], k) {
            kvec[km] = 1;
        }
    }
}

/// Collect sorted k-mer indices for all valid k-mers in `seq`.
/// Equivalent to C++ `tax_karray`.
fn tax_karray(seq: &[u8], k: usize) -> Vec<usize> {
    let klen = seq.len().saturating_sub(k - 1);
    let mut arr: Vec<usize> = (0..klen)
        .filter_map(|i| tax_kmer(&seq[i..], k))
        .collect();
    arr.sort_unstable();
    arr
}

/// Find the genus that maximises `Σ log P(kmer | genus)` over the query's
/// sorted k-mer array.
///
/// Ties are broken by uniform reservoir sampling (each tied genus has equal
/// probability of being selected).  This matches the C++ behaviour exactly.
///
/// Returns `(best_genus_idx, max_logp)`.
/// Equivalent to C++ `get_best_genus`.
fn get_best_genus(
    karray: &[usize],
    n_kmers: usize,
    ngenus: usize,
    lgk: &[f32],
) -> (usize, f32) {
    let mut max_g = 0usize;
    let mut max_logp = f32::NEG_INFINITY;
    let mut nmax = 0u32;

    for g in 0..ngenus {
        let lgk_v = &lgk[g * n_kmers..];
        let mut logp = 0.0f32;
        let mut early_exit = false;

        for &km in karray {
            logp += lgk_v[km];
            // Early-exit: this genus can't beat the current best.
            if logp < max_logp {
                early_exit = true;
                break;
            }
        }
        if early_exit {
            continue;
        }

        if max_logp > 0.0 || logp > max_logp {
            // New maximum.
            max_logp = logp;
            max_g = g;
            nmax = 1;
        } else if logp == max_logp {
            // Tied: keep with probability 1/nmax (reservoir sampling).
            nmax += 1;
            if rand::random::<f64>() < 1.0 / nmax as f64 {
                max_g = g;
            }
        }
    }
    (max_g, max_logp)
}

/// Classify a single query sequence (and optionally its reverse complement).
///
/// Returns `None` if the sequence is shorter than `MIN_SEQ_LEN`.
/// Otherwise returns `(best_genus, karray)` where `karray` is the sorted
/// k-mer index array of the winning orientation (forward or RC).
fn classify_seq(
    seq: &[u8],
    rc: Option<&[u8]>,
    k: usize,
    n_kmers: usize,
    ngenus: usize,
    lgk: &[f32],
) -> Option<(usize, Vec<usize>)> {
    if seq.len() < MIN_SEQ_LEN {
        return None;
    }

    let karray = tax_karray(seq, k);
    let (mut best_g, mut best_logp) = get_best_genus(&karray, n_kmers, ngenus, lgk);
    let mut best_karray = karray;

    if let Some(rc_seq) = rc {
        let karray_rc = tax_karray(rc_seq, k);
        let (g_rc, logp_rc) = get_best_genus(&karray_rc, n_kmers, ngenus, lgk);
        if logp_rc > best_logp {
            best_g = g_rc;
            best_logp = logp_rc;
            best_karray = karray_rc;
        }
        let _ = best_logp; // suppress unused warning after RC path
    }

    Some((best_g, best_karray))
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Assign taxonomy to each query sequence using a naive Bayes k-mer classifier.
///
/// # Parameters
/// - `seqs`: query sequences in ASCII (A/C/G/T).
/// - `rcs`: reverse complements of each query; supply `&[]` if `try_rc` is false.
///   When `try_rc` is true this slice must have the same length as `seqs`.
/// - `refs`: reference sequences in ASCII.
/// - `ref_to_genus`: 0-indexed genus for each reference (same length as `refs`).
/// - `genus_tax`: flat row-major matrix of shape `[ngenus][nlevel]` mapping each
///   genus to its taxon IDs at each level.  Used only for bootstrap agreement
///   counting; pass `&[]` with `nlevel = 0` to skip.
/// - `nlevel`: number of taxonomic levels (columns of `genus_tax`).
/// - `try_rc`: if true, also classify each sequence's reverse complement and
///   keep whichever orientation scores higher.
/// - `verbose`: print progress to stderr.
///
/// Equivalent to C++ `C_assign_taxonomy2`.
pub fn assign_taxonomy(
    seqs: &[&[u8]],
    rcs: &[&[u8]],
    refs: &[&[u8]],
    ref_to_genus: &[usize],
    genus_tax: &[usize],
    nlevel: usize,
    try_rc: bool,
    verbose: bool,
) -> Result<TaxonomyResult, String> {
    // ---- Validate ----
    let nseq = seqs.len();
    if nseq == 0 {
        return Err("No sequences provided to classify.".into());
    }
    let nref = refs.len();
    if nref != ref_to_genus.len() {
        return Err(format!(
            "Length mismatch: {} references but {} ref_to_genus entries.",
            nref,
            ref_to_genus.len()
        ));
    }
    let ngenus = if nlevel > 0 { genus_tax.len() / nlevel } else { 0 };
    if nlevel > 0 && genus_tax.len() != ngenus * nlevel {
        return Err(format!(
            "genus_tax length {} is not divisible by nlevel {}.",
            genus_tax.len(),
            nlevel
        ));
    }
    if try_rc && rcs.len() != nseq {
        return Err(format!(
            "try_rc=true but rcs has {} entries (expected {}).",
            rcs.len(),
            nseq
        ));
    }
    for (i, &g) in ref_to_genus.iter().enumerate() {
        if g >= ngenus && ngenus > 0 {
            return Err(format!(
                "ref_to_genus[{i}] = {g} is out of range (ngenus = {ngenus})."
            ));
        }
    }

    let k = TAX_K;
    let n_kmers = 1usize << (2 * k);

    // ---- Build genus k-mer counts (M_g + 1) and cross-genus kmer prior ----
    // genus_num_plus1[g] = number of refs in genus g, plus 1.
    let mut genus_count = vec![1.0f32; ngenus]; // starts at 1 (the "+1")
    for &g in ref_to_genus {
        genus_count[g] += 1.0;
    }

    // kmer_prior[km] = (number of genera in which km appears) + 0.5, / (1 + nref)
    // lgk[g * n_kmers + km] starts as raw count of refs-in-genus-g that contain km.
    let mut kmer_prior = vec![0.0f32; n_kmers];
    let mut lgk = vec![0.0f32; ngenus * n_kmers];
    let mut ref_kv = vec![0u8; n_kmers];

    for (i, seq) in refs.iter().enumerate() {
        tax_kvec(seq, k, &mut ref_kv);
        let g = ref_to_genus[i];
        let lgk_v = &mut lgk[g * n_kmers..(g + 1) * n_kmers];
        for km in 0..n_kmers {
            if ref_kv[km] != 0 {
                lgk_v[km] += 1.0;
                kmer_prior[km] += 1.0;
            }
        }
    }

    // Finalise kmer prior.
    for km in 0..n_kmers {
        kmer_prior[km] = (kmer_prior[km] + 0.5) / (1.0 + nref as f32);
    }

    // Convert counts to log probabilities: log((count + prior) / genus_num_plus1).
    for g in 0..ngenus {
        let lgk_v = &mut lgk[g * n_kmers..(g + 1) * n_kmers];
        let denom = genus_count[g];
        for km in 0..n_kmers {
            lgk_v[km] = ((lgk_v[km] + kmer_prior[km]) / denom).ln();
        }
    }

    if verbose {
        eprintln!("Finished processing reference fasta.");
    }

    // ---- Classify each query sequence in parallel ----
    // Each element: Option<(best_genus, karray)>.
    let classified: Vec<Option<(usize, Vec<usize>)>> = (0..nseq)
        .into_par_iter()
        .map(|j| {
            let rc = if try_rc { Some(rcs[j]) } else { None };
            classify_seq(seqs[j], rc, k, n_kmers, ngenus, &lgk)
        })
        .collect();

    // ---- Bootstrap in parallel ----
    let boot_results: Vec<(Vec<u32>, Vec<Option<usize>>)> = (0..nseq)
        .into_par_iter()
        .map(|j| {
            let mut boot_counts = vec![0u32; nlevel];
            let mut boot_taxa: Vec<Option<usize>> = vec![None; NBOOT];

            if let Some((best_g, ref karray)) = classified[j] {
                let arraylen = karray.len();
                let sample_size = arraylen / 8;

                for b in 0..NBOOT {
                    // Sample sample_size k-mers uniformly at random from karray.
                    let bootarray: Vec<usize> = (0..sample_size)
                        .map(|_| {
                            let idx = (arraylen as f64 * rand::random::<f64>()) as usize;
                            karray[idx.min(arraylen - 1)]
                        })
                        .collect();

                    let (boot_g, _) = get_best_genus(&bootarray, n_kmers, ngenus, &lgk);
                    boot_taxa[b] = Some(boot_g);

                    // Count levels where boot_g and best_g agree, stopping at first mismatch.
                    if nlevel > 0 {
                        for l in 0..nlevel {
                            let bg = genus_tax[boot_g * nlevel + l];
                            let mg = genus_tax[best_g * nlevel + l];
                            if bg == mg {
                                boot_counts[l] += 1;
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
            (boot_counts, boot_taxa)
        })
        .collect();

    // ---- Assemble result ----
    let assignments: Vec<Option<usize>> = classified.iter().map(|c| c.as_ref().map(|x| x.0)).collect();
    let boot_counts: Vec<Vec<u32>> = boot_results.iter().map(|(bc, _)| bc.clone()).collect();
    let boot_taxa: Vec<Vec<Option<usize>>> = boot_results.into_iter().map(|(_, bt)| bt).collect();

    Ok(TaxonomyResult { assignments, boot_counts, boot_taxa })
}
