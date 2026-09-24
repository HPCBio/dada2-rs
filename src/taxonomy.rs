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
//! - Random numbers are generated on-the-fly per sequence via a local `SmallRng`
//!   keyed on the sequence, so results do not depend on input order (issue #187)
//!   instead of being pre-allocated by `Rcpp::runif`.
//! - Runs are reproducible by default: the seed has a fixed value unless the
//!   caller overrides it (see [`assign_taxonomy`]).
//! - `RcppParallel::parallelFor` is replaced by Rayon `par_iter`.
//! - All indexing is 0-based; callers should add 1 if they need R-style output.
//! - `Rcpp::checkUserInterrupt()` is removed (no R event loop).

use std::collections::{BTreeSet, HashMap};

use rand::SeedableRng;
use rand::distributions::{Distribution, Standard};
use rand::rngs::SmallRng;

use crate::sequence_table::md5_seed;
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
    #[allow(dead_code)]
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
    let mut arr: Vec<usize> = (0..klen).filter_map(|i| tax_kmer(&seq[i..], k)).collect();
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
    rng: &mut SmallRng,
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
            let u: f64 = Standard.sample(rng);
            if u < 1.0 / nmax as f64 {
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
    rng: &mut SmallRng,
) -> Option<(usize, Vec<usize>)> {
    if seq.len() < MIN_SEQ_LEN {
        return None;
    }

    let karray = tax_karray(seq, k);
    let (mut best_g, mut best_logp) = get_best_genus(&karray, n_kmers, ngenus, lgk, rng);
    let mut best_karray = karray;

    if let Some(rc_seq) = rc {
        let karray_rc = tax_karray(rc_seq, k);
        let (g_rc, logp_rc) = get_best_genus(&karray_rc, n_kmers, ngenus, lgk, rng);
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
// RNG construction
// ---------------------------------------------------------------------------

/// Build a per-sequence `SmallRng`.
///
/// The stream is derived from the sequence, so a sequence's result depends only
/// on the sequence, the reference and the seed -- not on its position in the
/// input, how many others were submitted, or how Rayon scheduled them.
///
/// Keying on the position instead (`seed ^ index`) was reproducible only for a
/// fixed input in a fixed order: shuffling 3994 queries moved 7.6% of the
/// assignments (issue #187). That is the same symptom as R DADA2's
/// `assignTaxonomy` (benjjneb/dada2#1115), from a different cause -- theirs is
/// one C-side stream advancing across sequences.
///
/// There is no unseeded mode. Sampling from entropy would only make a run
/// unreproducible; a caller who wants to measure the bootstrap's sensitivity
/// varies the seed instead, which gives the same spread and can be repeated.
#[inline]
fn make_rng(seed: u64, seq: &[u8]) -> SmallRng {
    SmallRng::seed_from_u64(seed ^ md5_seed(seq))
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
/// - `seed`: RNG seed. Each sequence's stream is derived from the sequence
///   itself, so output is identical regardless of input order, input set, or
///   Rayon thread scheduling.
/// - `verbose`: print progress to stderr.
///
/// Equivalent to C++ `C_assign_taxonomy2`.
pub fn assign_taxonomy(
    seqs: &[&[u8]],
    rcs: &[&[u8]],
    ref_db: &TaxonomyRef<'_>,
    opts: TaxonomyOptions,
) -> Result<TaxonomyResult, String> {
    let TaxonomyRef {
        refs,
        ref_to_genus,
        genus_tax,
        nlevel,
    } = *ref_db;
    let TaxonomyOptions {
        try_rc,
        seed,
        verbose,
    } = opts;
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
    let ngenus = genus_tax.len().checked_div(nlevel).unwrap_or(0);
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
    for km in kmer_prior.iter_mut().take(n_kmers) {
        *km = (*km + 0.5) / (1.0 + nref as f32);
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
    // Each sequence gets its own RNG, derived from the sequence itself.
    let classified: Vec<Option<(usize, Vec<usize>)>> = (0..nseq)
        .into_par_iter()
        .map(|j| {
            let mut rng = make_rng(seed, seqs[j]);
            let rc = if try_rc { Some(rcs[j]) } else { None };
            classify_seq(seqs[j], rc, k, n_kmers, ngenus, &lgk, &mut rng)
        })
        .collect();

    // ---- Bootstrap in parallel ----
    // Use a different RNG stream from classification by XOR-ing with a constant.
    let boot_results: Vec<(Vec<u32>, Vec<Option<usize>>)> = (0..nseq)
        .into_par_iter()
        .map(|j| {
            let mut rng = make_rng(seed ^ 0xdead_beef_cafe_0000, seqs[j]);
            let mut boot_counts = vec![0u32; nlevel];
            let mut boot_taxa: Vec<Option<usize>> = vec![None; NBOOT];

            if let Some((best_g, ref karray)) = classified[j] {
                let arraylen = karray.len();
                let sample_size = arraylen / 8;

                for b in boot_taxa.iter_mut().take(NBOOT) {
                    // Sample sample_size k-mers uniformly at random from karray.
                    let bootarray: Vec<usize> = (0..sample_size)
                        .map(|_| {
                            let u: f64 = Standard.sample(&mut rng);
                            let idx = (arraylen as f64 * u) as usize;
                            karray[idx.min(arraylen - 1)]
                        })
                        .collect();

                    let (boot_g, _) = get_best_genus(&bootarray, n_kmers, ngenus, &lgk, &mut rng);
                    *b = Some(boot_g);

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
    let assignments: Vec<Option<usize>> =
        classified.iter().map(|c| c.as_ref().map(|x| x.0)).collect();
    let boot_counts: Vec<Vec<u32>> = boot_results.iter().map(|(bc, _)| bc.clone()).collect();
    let boot_taxa: Vec<Vec<Option<usize>>> = boot_results.into_iter().map(|(_, bt)| bt).collect();

    Ok(TaxonomyResult {
        assignments,
        boot_counts,
        boot_taxa,
    })
}

// ---------------------------------------------------------------------------
// assignSpecies — exact-match species assignment
// ---------------------------------------------------------------------------

/// Reference database for [`assign_taxonomy`].
#[derive(Clone, Copy)]
pub struct TaxonomyRef<'a> {
    pub refs: &'a [&'a [u8]],
    pub ref_to_genus: &'a [usize],
    pub genus_tax: &'a [usize],
    pub nlevel: usize,
}

/// Runtime options for [`assign_taxonomy`].
#[derive(Clone, Copy)]
pub struct TaxonomyOptions {
    pub try_rc: bool,
    pub seed: u64,
    pub verbose: bool,
}

/// Reference database for [`assign_species`].
#[derive(Clone, Copy)]
pub struct SpeciesRef<'a> {
    pub ref_seqs: &'a [&'a [u8]],
    pub ref_genus: &'a [&'a str],
    pub ref_species: &'a [&'a str],
}

/// Runtime options for [`assign_species`].
#[derive(Clone, Copy)]
pub struct SpeciesOptions {
    pub max_species: usize,
    pub try_rc: bool,
    pub verbose: bool,
}

/// Per-query result from [`assign_species`].
pub struct SpeciesHit {
    /// Unambiguous genus (after Escherichia/Shigella merging); `None` if
    /// the query matched references with conflicting genera or had no hit.
    pub genus: Option<String>,
    /// Species assignment; `None` if no hit, or if the number of distinct
    /// matching species exceeds `max_species`.
    pub species: Option<String>,
}

fn rc_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            _ => b'N',
        })
        .collect()
}

fn escsh(s: &str) -> &str {
    if s.contains("Escherichia") || s.contains("Shigella") {
        "Escherichia/Shigella"
    } else {
        s
    }
}

/// Assign species to each query by exact sequence matching.
///
/// Reference FASTA is expected in the format `>SeqID genus species`; callers
/// pre-parse genus and species strings and pass them as parallel slices.
///
/// `max_species`:
/// - `1` — unambiguous only (R's default `allowMultiple=FALSE`)
/// - `0` — unlimited (R's `allowMultiple=TRUE`)
/// - `N > 1` — up to N distinct species joined with `"/"`
///
/// Equivalent to R's `assignSpecies`.
/// For each query, the set of reference indices whose sequence *contains* it.
///
/// Mirrors R's `vcountPDict(PDict(seqs), sread(refs)) > 0`. `PDict` requires
/// equal-width patterns, so R batches queries by length; the same grouping is
/// what makes a rolling hash usable here.
///
/// Rabin-Karp per length class: hash each `len`-wide window of a reference once
/// and look it up, rather than searching each (query, reference) pair
/// separately -- 352k references and 517 Mbases make the pairwise form
/// hopeless. Hash collisions are resolved by comparing the bytes, so the
/// result is exact.
///
/// `n_fwd` is the number of forward queries; patterns at or above it are
/// reverse complements and fold back onto `pattern_index - n_fwd`.
fn find_containing_refs(
    patterns: &[&[u8]],
    ref_seqs: &[&[u8]],
    n_fwd: usize,
) -> Vec<BTreeSet<usize>> {
    const BASE: u64 = 0x100_0000_01b3; // FNV prime, arbitrary but odd

    // Length class -> (hash -> pattern indices), plus the class's leading
    // coefficient for the rolling update.
    let mut by_len: HashMap<usize, HashMap<u64, Vec<usize>>> = HashMap::new();
    for (pi, p) in patterns.iter().enumerate() {
        if p.is_empty() {
            continue;
        }
        let h = p.iter().fold(0u64, |acc, &b| {
            acc.wrapping_mul(BASE).wrapping_add(b as u64)
        });
        by_len
            .entry(p.len())
            .or_default()
            .entry(h)
            .or_default()
            .push(pi);
    }
    /// One query-length class: the width, `BASE^(width-1)` for the rolling
    /// update, and that width's hash → pattern-index table.
    type LenClass<'a> = (usize, u64, &'a HashMap<u64, Vec<usize>>);

    let classes: Vec<LenClass<'_>> = by_len
        .iter()
        .map(|(&len, map)| {
            // BASE^(len-1), the weight of the byte leaving the window.
            let high = (0..len.saturating_sub(1)).fold(1u64, |acc, _| acc.wrapping_mul(BASE));
            (len, high, map)
        })
        .collect();

    ref_seqs
        .par_iter()
        .enumerate()
        .fold(
            || vec![BTreeSet::new(); n_fwd],
            |mut acc, (ri, &r)| {
                for &(len, high, map) in &classes {
                    if r.len() < len {
                        continue;
                    }
                    let mut h = r[..len]
                        .iter()
                        .fold(0u64, |a, &b| a.wrapping_mul(BASE).wrapping_add(b as u64));
                    let mut start = 0usize;
                    loop {
                        if let Some(cands) = map.get(&h) {
                            for &pi in cands {
                                if patterns[pi] == &r[start..start + len] {
                                    acc[if pi >= n_fwd { pi - n_fwd } else { pi }].insert(ri);
                                }
                            }
                        }
                        if start + len >= r.len() {
                            break;
                        }
                        h = h
                            .wrapping_sub((r[start] as u64).wrapping_mul(high))
                            .wrapping_mul(BASE)
                            .wrapping_add(r[start + len] as u64);
                        start += 1;
                    }
                }
                acc
            },
        )
        .reduce(
            || vec![BTreeSet::new(); n_fwd],
            |mut a, b| {
                for (dst, src) in a.iter_mut().zip(b) {
                    dst.extend(src);
                }
                a
            },
        )
}

pub fn assign_species(
    seqs: &[&[u8]],
    ref_db: &SpeciesRef<'_>,
    opts: SpeciesOptions,
) -> Vec<SpeciesHit> {
    let SpeciesRef {
        ref_seqs,
        ref_genus,
        ref_species,
    } = *ref_db;
    let SpeciesOptions {
        max_species,
        try_rc,
        verbose,
    } = opts;
    assert_eq!(ref_seqs.len(), ref_genus.len());
    assert_eq!(ref_seqs.len(), ref_species.len());

    // A query matches a reference when it occurs *within* it, not when the two
    // are equal: R's `vcountPDict(PDict(seqs), sread(refs)) > 0` (taxonomy.R).
    // The distinction is the whole feature -- species references are full-length
    // 16S and queries are sub-region amplicons, so equality never holds (#200).
    let mut patterns: Vec<&[u8]> = seqs.to_vec();
    let rcs: Vec<Vec<u8>> = if try_rc {
        seqs.iter().map(|q| rc_seq(q)).collect()
    } else {
        Vec::new()
    };
    // Query i's reverse complement is pattern `seqs.len() + i`. R reverse-
    // complements the *reference* instead; searching for RC(query) in the
    // reference is the same predicate.
    patterns.extend(rcs.iter().map(|v| v.as_slice()));

    let hits_per_query = find_containing_refs(&patterns, ref_seqs, seqs.len());

    let mut results = Vec::with_capacity(seqs.len());
    let mut n_assigned = 0usize;

    for (qi, &query) in seqs.iter().enumerate() {
        let _ = query;
        let hit_indices: Vec<usize> = hits_per_query[qi].iter().copied().collect();

        if hit_indices.is_empty() {
            results.push(SpeciesHit {
                genus: None,
                species: None,
            });
            continue;
        }

        // Genus: unique values after E/Shigella merging; unambiguous only.
        let mut genera: Vec<&str> = hit_indices.iter().map(|&i| escsh(ref_genus[i])).collect();
        genera.sort_unstable();
        genera.dedup();
        let genus = if genera.len() == 1 {
            Some(genera[0].to_string())
        } else {
            None
        };

        // Species: unique values; accept if count ≤ max_species (0 = unlimited).
        let mut spp: Vec<&str> = hit_indices.iter().map(|&i| ref_species[i]).collect();
        spp.sort_unstable();
        spp.dedup();
        let species = if max_species == 0 || spp.len() <= max_species {
            Some(spp.join("/"))
        } else {
            None
        };

        if species.is_some() {
            n_assigned += 1;
        }
        results.push(SpeciesHit { genus, species });
    }

    if verbose {
        eprintln!(
            "{} out of {} were assigned to the species level.",
            n_assigned,
            seqs.len()
        );
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Queries, references, reference→genus, genus→taxon-IDs.
    type Fixture = (Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<usize>, Vec<usize>);

    /// Build a reference set of `ngenus` genera, each with several near-identical
    /// members, plus queries that sit between two genera so the bootstrap has
    /// something to disagree about. Without that ambiguity the test could not
    /// fail even with a broken RNG key.
    fn fixture() -> Fixture {
        let base: Vec<u8> = (0..400u32)
            .map(|i| b"ACGT"[((i * 7 + i / 3) % 4) as usize])
            .collect();
        let mut refs = Vec::new();
        let mut ref_to_genus = Vec::new();
        for g in 0..4usize {
            for m in 0..3usize {
                let mut r = base.clone();
                for p in 0..40usize {
                    let idx = (g * 97 + m * 13 + p * 9) % r.len();
                    r[idx] = b"ACGT"[(g + p + m) % 4];
                }
                refs.push(r);
                ref_to_genus.push(g);
            }
        }
        // Queries: blends of genus 0 and genus 1 members, so the bootstrap
        // draw decides between them.
        let mut queries = Vec::new();
        for q in 0..12usize {
            let mut s = refs[0].clone();
            for p in 0..(20 + q) {
                let idx = (q * 31 + p * 17) % s.len();
                s[idx] = refs[3][idx];
            }
            queries.push(s);
        }
        let genus_tax: Vec<usize> = (0..4usize).flat_map(|g| [g / 2, g]).collect();
        (queries, refs, ref_to_genus, genus_tax)
    }

    fn classify(
        order: &[usize],
        queries: &[Vec<u8>],
        refs: &[Vec<u8>],
        r2g: &[usize],
        gt: &[usize],
    ) -> Vec<(Vec<u8>, Vec<Option<usize>>)> {
        let seqs: Vec<&[u8]> = order.iter().map(|&i| queries[i].as_slice()).collect();
        let refv: Vec<&[u8]> = refs.iter().map(|r| r.as_slice()).collect();
        let res = assign_taxonomy(
            &seqs,
            &[],
            &TaxonomyRef {
                refs: &refv,
                ref_to_genus: r2g,
                genus_tax: gt,
                nlevel: 2,
            },
            TaxonomyOptions {
                try_rc: false,
                seed: 42,
                verbose: false,
            },
        )
        .expect("classification failed");
        order
            .iter()
            .enumerate()
            .map(|(pos, &i)| (queries[i].clone(), res.boot_taxa[pos].clone()))
            .collect()
    }

    /// A query matches a reference it is *contained in*, not only one it equals
    /// (issue #200; R does `vcountPDict(...) > 0`).
    ///
    /// The flanks are what make this able to fail: with whole-sequence equality
    /// the query is never found, which is what shipped. Species references are
    /// full-length 16S and queries are sub-region amplicons, so containment is
    /// the ordinary case and equality is the degenerate one.
    #[test]
    fn species_match_is_containment_not_equality() {
        let query: Vec<u8> = b"ACGTACGTACGTTTGGCCAA".to_vec();
        let mut reference = b"TTTTTTGGGG".to_vec();
        reference.extend_from_slice(&query);
        reference.extend_from_slice(b"CCCCAAAATT");
        assert!(
            reference.len() > query.len(),
            "flanks must make this a proper substring"
        );

        let refs: Vec<&[u8]> = vec![reference.as_slice()];
        let hits = assign_species(
            &[query.as_slice()],
            &SpeciesRef {
                ref_seqs: &refs,
                ref_genus: &["Blautia"],
                ref_species: &["coccoides"],
            },
            SpeciesOptions {
                max_species: 1,
                try_rc: false,
                verbose: false,
            },
        );
        assert_eq!(hits[0].genus.as_deref(), Some("Blautia"));
        assert_eq!(hits[0].species.as_deref(), Some("coccoides"));
    }

    /// Equality is just containment with empty flanks, and must keep working.
    #[test]
    fn species_match_still_accepts_an_exact_reference() {
        let query: Vec<u8> = b"ACGTACGTACGTTTGGCCAA".to_vec();
        let refs: Vec<&[u8]> = vec![query.as_slice()];
        let hits = assign_species(
            &[query.as_slice()],
            &SpeciesRef {
                ref_seqs: &refs,
                ref_genus: &["Blautia"],
                ref_species: &["coccoides"],
            },
            SpeciesOptions {
                max_species: 1,
                try_rc: false,
                verbose: false,
            },
        );
        assert_eq!(hits[0].species.as_deref(), Some("coccoides"));
    }

    /// `--try-rc` has to use the same predicate: R reverse-complements the
    /// reference, we reverse-complement the query, and containment makes those
    /// the same test.
    #[test]
    fn species_try_rc_also_matches_by_containment() {
        let query: Vec<u8> = b"ACGTACGTACGTTTGGCCAA".to_vec();
        let rc = rc_seq(&query);
        let mut reference = b"TTTTTTGGGG".to_vec();
        reference.extend_from_slice(&rc);
        reference.extend_from_slice(b"CCCCAAAATT");

        let refs: Vec<&[u8]> = vec![reference.as_slice()];
        let opts = |try_rc| SpeciesOptions {
            max_species: 1,
            try_rc,
            verbose: false,
        };
        let sref = SpeciesRef {
            ref_seqs: &refs,
            ref_genus: &["Blautia"],
            ref_species: &["coccoides"],
        };
        assert_eq!(
            assign_species(&[query.as_slice()], &sref, opts(false))[0].species,
            None,
            "forward-only must not find the reverse-complemented reference"
        );
        assert_eq!(
            assign_species(&[query.as_slice()], &sref, opts(true))[0]
                .species
                .as_deref(),
            Some("coccoides")
        );
    }

    /// A seeded run must give each sequence the same answer no matter where it
    /// sits in the input (issue #187; the R analogue is benjjneb/dada2#1115).
    ///
    /// Compares `boot_taxa` rather than the final call: a changed RNG stream
    /// always perturbs the per-replicate draws, while a final assignment only
    /// moves when a replicate crosses the confidence threshold. Asserting on the
    /// final call would make this test far less able to fail.
    #[test]
    fn seeded_assignment_is_independent_of_input_order() {
        let (queries, refs, r2g, gt) = fixture();
        let n = queries.len();
        let forward: Vec<usize> = (0..n).collect();
        let reversed: Vec<usize> = (0..n).rev().collect();
        let rotated: Vec<usize> = (0..n).map(|i| (i + 5) % n).collect();

        let a = classify(&forward, &queries, &refs, &r2g, &gt);
        for order in [reversed, rotated] {
            let b = classify(&order, &queries, &refs, &r2g, &gt);
            for (seq, boots) in &b {
                let (_, expect) = a.iter().find(|(s, _)| s == seq).expect("sequence missing");
                assert_eq!(
                    boots, expect,
                    "bootstrap draws changed when the input order changed; \
                     the per-sequence RNG is keyed on position, not on the sequence"
                );
            }
        }
    }

    /// The stream must also not depend on which *other* sequences were submitted.
    #[test]
    fn seeded_assignment_is_independent_of_the_input_set() {
        let (queries, refs, r2g, gt) = fixture();
        let all: Vec<usize> = (0..queries.len()).collect();
        let subset: Vec<usize> = vec![7, 2, 9];

        let a = classify(&all, &queries, &refs, &r2g, &gt);
        let b = classify(&subset, &queries, &refs, &r2g, &gt);
        for (seq, boots) in &b {
            let (_, expect) = a.iter().find(|(s, _)| s == seq).expect("sequence missing");
            assert_eq!(
                boots, expect,
                "bootstrap draws depend on the rest of the input set"
            );
        }
    }
}
