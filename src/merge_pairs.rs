//! Paired-end read merging — port of DADA2's `mergePairs` from `paired.R`.
//!
//! ## Workflow
//!
//! For each sample (a matched set of forward/reverse FASTQ + dada JSON files):
//!
//! 1. The forward and reverse FASTQ files are re-dereplicated to recover the
//!    read → unique-index mapping.
//! 2. The dada JSON files supply the unique-index → ASV-index map (`map` field,
//!    only present when `dada --show-map` was used) and the ASV sequences.
//! 3. For every read, the two maps are composed to give (fwd_asv, rev_asv).
//!    Reads where either direction is unassigned (map entry = `null`) are
//!    silently dropped.
//! 4. Distinct (fwd_asv, rev_asv) pairs are counted, then each is attempted:
//!    the forward ASV sequence is aligned (ends-free Needleman-Wunsch) against
//!    the reverse-complement of the reverse ASV sequence.  If the overlap is
//!    long enough, has few enough mismatches, and no indels, the merge is
//!    accepted and the merged amplicon sequence is assembled.
//!
//! ## Unique-index ordering guarantee
//!
//! DADA2's dereplication assigns unique indices in "first-seen" order across
//! the read stream.  This project's `dereplicate()` function uses a
//! chunk-parallel fold/reduce pattern that preserves the same first-seen
//! ordering regardless of thread count (rayon's reduce respects element order
//! for contiguous partitions).  Re-dereplicating the same FASTQ file therefore
//! yields identical unique indices — no need to save a separate derep JSON.
//!
//! **Caveat**: the dada JSON **must** have been produced with `--show-map`;
//! without that flag the `map` field is absent and merging is impossible.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::path::Path;

use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};

use crate::derep::dereplicate;
use crate::misc::{intstr, nt_decode};
use crate::nwalign::align_endsfree;

// ---------------------------------------------------------------------------
// Alignment constants (same as core DADA2 algorithm)
// ---------------------------------------------------------------------------

const MATCH_SCORE: i32 = 5;
const MISMATCH: i32 = -4;
const GAP_P: i32 = -8;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Tuning parameters for paired-end merging.
pub struct MergeParams {
    /// Minimum overlap length (nmatch + nmismatch + nindel ≥ min_overlap).
    pub min_overlap: u32,
    /// Maximum mismatches in the overlap region.
    pub max_mismatch: u32,
    /// When true, include rejected merges in the output (with `accept = false`).
    pub return_rejects: bool,
    /// When true, concatenate fwd + N-spacer + RC(rev) instead of merging.
    pub just_concatenate: bool,
    /// Length of the N spacer used when `just_concatenate` is true.
    pub concat_nnn_len: usize,
    /// When true, trim portions of fwd/rev that overhang past the other read.
    pub trim_overhang: bool,
    /// Phred quality-score offset for FASTQ re-dereplication.
    pub phred_offset: u8,
    /// Print per-sample progress to stderr.
    pub verbose: bool,
}

// ---------------------------------------------------------------------------
// dada JSON deserialization (only the fields we need)
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct AsvJson {
    sequence: String,
}

#[derive(Deserialize)]
struct DadaJsonInput {
    asvs: Vec<AsvJson>,
    /// unique-index → ASV-index mapping; only present when `--show-map` was used.
    map: Option<Vec<Option<usize>>>,
}

// ---------------------------------------------------------------------------
// Output structures
// ---------------------------------------------------------------------------

/// One accepted (or rejected) merged amplicon sequence.
#[derive(Serialize)]
pub struct MergedPair {
    /// Merged amplicon sequence (empty string when `accept = false`).
    pub sequence: String,
    /// Number of read-pairs that produced this merge.
    pub abundance: u64,
    /// 0-based index of the forward ASV.
    pub forward: usize,
    /// 0-based index of the reverse ASV.
    pub reverse: usize,
    /// Matching positions in the overlap region.
    pub nmatch: u32,
    /// Mismatching positions in the overlap region.
    pub nmismatch: u32,
    /// Indel positions in the overlap region.
    pub nindel: u32,
    /// Whether this merge met all acceptance criteria.
    pub accept: bool,
}

/// Merging results for one sample.
#[derive(Serialize)]
pub struct SampleMergeResult {
    /// Sample name (derived from the forward dada JSON file stem).
    pub sample: String,
    /// Total read-pairs where both directions were assigned to an ASV.
    pub total_pairs: u64,
    /// Read-pairs that produced an accepted merge.
    pub accepted_pairs: u64,
    /// Number of distinct merged sequences.
    pub num_merged: usize,
    /// Merged (and optionally rejected) pairs, sorted by abundance descending.
    pub merged: Vec<MergedPair>,
}

// ---------------------------------------------------------------------------
// Core helpers
// ---------------------------------------------------------------------------

/// Reverse-complement an ASCII DNA sequence (A/C/G/T/N, case-insensitive).
fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            _ => b'N',
        })
        .map(|b| b as char)
        .collect()
}

/// Analyse the overlap region in a ends-free NW alignment of fwd vs RC(rev).
///
/// Returns `(nmatch, nmismatch, nindel, ov_left, ov_right)` where `ov_left`
/// and `ov_right` are the first/last alignment column where **both** strands
/// have a non-gap base.  Returns `None` if there is no such overlap.
fn analyze_overlap(al0: &[u8], al1: &[u8]) -> Option<(u32, u32, u32, usize, usize)> {
    let n = al0.len();

    // First and last columns where both strands have a base.
    let left = (0..n).find(|&i| al0[i] != b'-' && al1[i] != b'-')?;
    let right = (0..n).rev().find(|&i| al0[i] != b'-' && al1[i] != b'-')?;

    if left > right {
        return None;
    }

    let mut nmatch = 0u32;
    let mut nmismatch = 0u32;
    let mut nindel = 0u32;

    for i in left..=right {
        match (al0[i] == b'-', al1[i] == b'-') {
            (true, true) => {} // shouldn't occur in a valid alignment
            (false, false) => {
                if al0[i] == al1[i] {
                    nmatch += 1;
                } else {
                    nmismatch += 1;
                }
            }
            _ => nindel += 1,
        }
    }

    Some((nmatch, nmismatch, nindel, left, right))
}

/// Assemble the merged amplicon sequence from the alignment.
///
/// `prefer_fwd = true` (R's default `prefer = 1`) uses the forward strand in
/// the overlap region.
///
/// Without `trim_overhang`:
/// - fwd bases before the overlap are included (fwd prefix).
/// - RC(rev) bases after the overlap are included (rev suffix).
/// - Any fwd bases *after* the overlap (fwd right-overhang) and any RC(rev)
///   bases *before* the overlap (rcrev left-overhang) are also included — this
///   matches R's default behaviour.
///
/// With `trim_overhang`:
/// - fwd right-overhang and rcrev left-overhang are omitted.
fn build_merged(
    al0: &[u8],
    al1: &[u8],
    ov_left: usize,
    ov_right: usize,
    trim_overhang: bool,
    prefer_fwd: bool,
) -> String {
    let n = al0.len();
    let mut result: Vec<u8> = Vec::with_capacity(n);

    // --- Region before the overlap ---
    for i in 0..ov_left {
        match (al0[i] == b'-', al1[i] == b'-') {
            // fwd has a base, rcrev has a gap → fwd prefix (always include)
            (false, true) => result.push(nt_decode(al0[i])),
            // rcrev has a base, fwd has a gap → rcrev left-overhang
            (true, false) => {
                if !trim_overhang {
                    result.push(nt_decode(al1[i]));
                }
            }
            // Both gap or both base outside the overlap shouldn't happen in a
            // well-formed ends-free alignment, but handle gracefully.
            _ => {}
        }
    }

    // --- Overlap region ---
    for i in ov_left..=ov_right {
        match (al0[i] == b'-', al1[i] == b'-') {
            (false, false) => {
                // Both have a base: use the preferred strand.
                result.push(nt_decode(if prefer_fwd { al0[i] } else { al1[i] }));
            }
            // Gap in fwd within the overlap (indel): use rcrev base.
            (true, false) => result.push(nt_decode(al1[i])),
            // Gap in rcrev within the overlap (indel): use fwd base.
            (false, true) => result.push(nt_decode(al0[i])),
            (true, true) => {}
        }
    }

    // --- Region after the overlap ---
    for i in (ov_right + 1)..n {
        match (al0[i] == b'-', al1[i] == b'-') {
            // rcrev has a base, fwd has a gap → rcrev suffix (always include)
            (true, false) => result.push(nt_decode(al1[i])),
            // fwd has a base, rcrev has a gap → fwd right-overhang
            (false, true) => {
                if !trim_overhang {
                    result.push(nt_decode(al0[i]));
                }
            }
            _ => {}
        }
    }

    String::from_utf8_lossy(&result).into_owned()
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Process one sample: re-dereplicate FASTQs, load dada JSONs, merge pairs.
///
/// The four paths must correspond to the same biological sample.  Files are
/// opened, processed, and closed within this call; nothing is held across
/// samples.
pub fn merge_sample(
    sample_name: &str,
    fwd_dada_path: &Path,
    rev_dada_path: &Path,
    fwd_fastq_path: &Path,
    rev_fastq_path: &Path,
    params: &MergeParams,
    pool: &rayon::ThreadPool,
) -> io::Result<SampleMergeResult> {
    // ---- Load dada JSONs (plain or gzip-compressed) ----
    let fwd_dada: DadaJsonInput = crate::misc::read_json_file(fwd_dada_path)?;
    let rev_dada: DadaJsonInput = crate::misc::read_json_file(rev_dada_path)?;

    let fwd_map = fwd_dada.map.as_ref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: 'map' field is absent — re-run `dada` with `--show-map`",
                fwd_dada_path.display()
            ),
        )
    })?;
    let rev_map = rev_dada.map.as_ref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: 'map' field is absent — re-run `dada` with `--show-map`",
                rev_dada_path.display()
            ),
        )
    })?;

    // ---- Re-dereplicate FASTQs ----
    // The ordering of unique sequences is deterministic (first-seen globally)
    // regardless of thread count, so the indices here match those in the dada
    // JSON map.
    let is_gz = |p: &Path| {
        p.file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.ends_with(".gz"))
            .unwrap_or(false)
    };

    let fwd_derep = if is_gz(fwd_fastq_path) {
        dereplicate(
            MultiGzDecoder::new(File::open(fwd_fastq_path)?),
            params.phred_offset,
            pool,
            params.verbose,
        )?
    } else {
        dereplicate(
            File::open(fwd_fastq_path)?,
            params.phred_offset,
            pool,
            params.verbose,
        )?
    };

    let rev_derep = if is_gz(rev_fastq_path) {
        dereplicate(
            MultiGzDecoder::new(File::open(rev_fastq_path)?),
            params.phred_offset,
            pool,
            params.verbose,
        )?
    } else {
        dereplicate(
            File::open(rev_fastq_path)?,
            params.phred_offset,
            pool,
            params.verbose,
        )?
    };

    // ---- Validate pairing ----
    let n_reads = fwd_derep.map.len();
    if rev_derep.map.len() != n_reads {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Sample '{}': forward FASTQ has {} reads but reverse has {}",
                sample_name,
                n_reads,
                rev_derep.map.len()
            ),
        ));
    }

    // Sanity-check that the dada map sizes are plausible.
    if fwd_map.len() != fwd_derep.uniques.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Sample '{}': forward dada map length ({}) ≠ forward unique count ({}); \
                 check that the same FASTQ was used for both `dada` and `merge-pairs`",
                sample_name,
                fwd_map.len(),
                fwd_derep.uniques.len()
            ),
        ));
    }
    if rev_map.len() != rev_derep.uniques.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Sample '{}': reverse dada map length ({}) ≠ reverse unique count ({}); \
                 check that the same FASTQ was used for both `dada` and `merge-pairs`",
                sample_name,
                rev_map.len(),
                rev_derep.uniques.len()
            ),
        ));
    }

    // ---- Count (fwd_asv, rev_asv) pairs ----
    let mut pair_counts: HashMap<(usize, usize), u64> = HashMap::new();
    let mut total_pairs: u64 = 0;

    for i in 0..n_reads {
        let fu = fwd_derep.map[i];
        let ru = rev_derep.map[i];

        // fwd_map[fu] is the ASV index for the fu-th forward unique (None = unassigned).
        let fa = match fwd_map.get(fu).and_then(|x| *x) {
            Some(a) => a,
            None => continue,
        };
        let ra = match rev_map.get(ru).and_then(|x| *x) {
            Some(a) => a,
            None => continue,
        };

        *pair_counts.entry((fa, ra)).or_insert(0) += 1;
        total_pairs += 1;
    }

    // ---- Attempt merge for each distinct pair ----
    let mut merged: Vec<MergedPair> = Vec::with_capacity(pair_counts.len());
    let mut accepted_pairs: u64 = 0;

    for ((fi, ri), count) in &pair_counts {
        let fwd_seq = &fwd_dada.asvs[*fi].sequence;
        let rev_seq = &rev_dada.asvs[*ri].sequence;
        let rc_rev = reverse_complement(rev_seq);

        // Just-concatenate mode: no alignment required.
        if params.just_concatenate {
            let spacer = "N".repeat(params.concat_nnn_len);
            merged.push(MergedPair {
                sequence: format!("{fwd_seq}{spacer}{rc_rev}"),
                abundance: *count,
                forward: *fi,
                reverse: *ri,
                nmatch: 0,
                nmismatch: 0,
                nindel: 0,
                accept: true,
            });
            accepted_pairs += count;
            continue;
        }

        // Encode sequences for the NW aligner (1=A, 2=C, 3=G, 4=T, 5=N).
        let fwd_enc = intstr(fwd_seq.as_bytes());
        let rev_enc = intstr(rc_rev.as_bytes());

        // Ends-free NW alignment (band = -1 → unbanded).
        let al = align_endsfree(&fwd_enc, &rev_enc, MATCH_SCORE, MISMATCH, GAP_P, -1);

        let ov = analyze_overlap(&al[0], &al[1]);

        let (nmatch, nmismatch, nindel, ov_left, ov_right) = match ov {
            Some(v) => v,
            None => {
                // No overlap at all.
                if params.return_rejects {
                    merged.push(MergedPair {
                        sequence: String::new(),
                        abundance: *count,
                        forward: *fi,
                        reverse: *ri,
                        nmatch: 0,
                        nmismatch: 0,
                        nindel: 0,
                        accept: false,
                    });
                }
                continue;
            }
        };

        let overlap_len = nmatch + nmismatch + nindel;
        let accept = overlap_len >= params.min_overlap
            && nmismatch <= params.max_mismatch
            && nindel == 0;

        if !accept && !params.return_rejects {
            continue;
        }

        let sequence = if accept {
            accepted_pairs += count;
            build_merged(&al[0], &al[1], ov_left, ov_right, params.trim_overhang, true)
        } else {
            String::new()
        };

        merged.push(MergedPair {
            sequence,
            abundance: *count,
            forward: *fi,
            reverse: *ri,
            nmatch,
            nmismatch,
            nindel,
            accept,
        });
    }

    // Sort by abundance descending for readability.
    merged.sort_unstable_by(|a, b| b.abundance.cmp(&a.abundance));

    let num_merged = merged.iter().filter(|m| m.accept).count();

    Ok(SampleMergeResult {
        sample: sample_name.to_string(),
        total_pairs,
        accepted_pairs,
        num_merged,
        merged,
    })
}
