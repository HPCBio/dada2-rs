//! Paired-end read merging — port of DADA2's `mergePairs` from `paired.R`.
//!
//! ## Workflow
//!
//! For each sample (a matched set of forward/reverse FASTQ + dada JSON files):
//!
//! 1. The forward and reverse FASTQ files are re-dereplicated to recover the
//!    read → unique-index mapping.
//! 2. The dada JSON files supply the unique-index → ASV-index map (`map`
//!    field, always emitted by `dada` / `dada-pooled`) and the ASV sequences.
//! 3. For every read, the two maps are composed to give (fwd_asv, rev_asv).
//!    Reads where either direction is unassigned (map entry = `null`) are
//!    silently dropped.
//! 4. Distinct (fwd_asv, rev_asv) pairs are counted, then each is attempted:
//!    the forward ASV sequence is aligned (ends-free, unbanded
//!    Needleman-Wunsch, with R `mergePairs`'s own scores — see
//!    [`merge_scores`]) against the reverse-complement of the reverse ASV
//!    sequence.  If the overlap has at least `min_overlap` matching bases and
//!    at most `max_mismatch` mismatches and indels together, the merge is
//!    accepted and the merged amplicon sequence is assembled.
//!
//! ## Unique-index ordering guarantee
//!
//! `dereplicate()` returns uniques sorted by abundance descending (stable;
//! ties keep first-seen order) — matching R `derepFastq` ordering. The
//! chunk-parallel fold/reduce that builds the unique set runs in
//! deterministic order, so re-dereplicating the same FASTQ file yields
//! identical unique indices, regardless of thread count.  No need to save
//! a separate derep JSON: the dada JSON's `map` field references unique
//! indices that re-derepping the source FASTQ will reproduce exactly.
//!
//! **Caveat**: dada JSON files saved with the pre-sort first-seen ordering
//! cannot be merged with a post-sort dereplicate — re-derep the FASTQ
//! through the current pipeline first.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::path::Path;

use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};

use crate::derep::dereplicate;
use crate::misc::WithPath;
use crate::misc::{intstr, nt_decode};
use crate::nwalign::{AlignBuffers, align_endsfree_with_buf};

// ---------------------------------------------------------------------------
// Alignment scores
// ---------------------------------------------------------------------------

/// `(match, mismatch, gap)` for aligning a forward ASV against RC(reverse).
///
/// Not the denoising scores (5 / −4 / −8): R's `mergePairs` swaps in its own
/// (`R/paired.R`, "prioritize zero-mismatch merges"). At `maxMismatch == 0`
/// any mismatch or gap costs −64, so a short perfect overlap always beats a
/// longer imperfect one at another offset; under 5 / −4 / −8 the longer one
/// can win and the pair is then rejected (#272).
fn merge_scores(max_mismatch: u32) -> (i32, i32, i32) {
    if max_mismatch == 0 {
        (1, -64, -64)
    } else {
        (1, -8, -8)
    }
}

/// Align a forward ASV against RC(reverse) ends-free and unbanded, as R's
/// `mergePairs` does, and count the overlap (R `C_eval_pair`).
///
/// Returns `(nmatch, nmismatch, nindel, ov_left, ov_right)`, or `None` when
/// the reads do not overlap.
fn align_pair(
    fwd_seq: &str,
    rc_rev: &str,
    max_mismatch: u32,
    buf: &mut AlignBuffers,
) -> Option<(u32, u32, u32, usize, usize)> {
    let (match_score, mismatch, gap_p) = merge_scores(max_mismatch);
    align_endsfree_with_buf(
        &intstr(fwd_seq.as_bytes()),
        &intstr(rc_rev.as_bytes()),
        match_score,
        mismatch,
        gap_p,
        -1,
        buf,
    );
    analyze_overlap(&buf.al0, &buf.al1)
}

/// R's acceptance rule: at least `minOverlap` matching bases, and at most
/// `maxMismatch` mismatches and indels together.
fn accepts(nmatch: u32, nmismatch: u32, nindel: u32, params: &MergeParams) -> bool {
    nmatch >= params.min_overlap && nmismatch + nindel <= params.max_mismatch
}

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Tuning parameters for paired-end merging.
pub struct MergeParams {
    /// Minimum matching bases in the overlap (R `minOverlap`: `nmatch >= min_overlap`).
    pub min_overlap: u32,
    /// Maximum mismatches plus indels in the overlap (R `maxMismatch`).
    pub max_mismatch: u32,
    /// When true, include rejected merges in the output (with `accept = false`).
    pub return_rejects: bool,
    /// When true, concatenate fwd + N-spacer + RC(rev) instead of merging.
    pub just_concatenate: bool,
    /// When true, pairs that fail to merge (no overlap, or failing the overlap
    /// criteria) are rescued by concatenating fwd + N-spacer + RC(rev) and
    /// accepted, rather than being dropped. Useful for amplicons such as ITS
    /// whose reads may not overlap. Takes precedence over `return_rejects`.
    pub rescue_unmerged: bool,
    /// Length of the N spacer used when `just_concatenate` is true.
    pub concat_nnn_len: usize,
    /// When true, trim portions of fwd/rev that overhang past the other read.
    pub trim_overhang: bool,
    /// Phred quality-score offset for FASTQ re-dereplication.
    pub phred_offset: u8,
    /// When true, verify that the fwd and rev dada JSONs carry the same
    /// `sample` field, that it equals the resolved sample name, and that both
    /// FASTQ filenames contain the sample name as a substring.
    pub check_sample_ids: bool,
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
    /// Sample identifier; absent in dada JSONs produced before --sample-name.
    sample: Option<String>,
    /// File name (no directory) `dada` read: a FASTQ or a derep JSON.
    input_file: Option<String>,
    /// The FASTQ behind `input_file`; absent before #111.
    source_fastq: Option<String>,
    asvs: Vec<AsvJson>,
    /// unique-index → ASV-index mapping; absent in dada JSONs produced
    /// before the map became part of the default output.
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
    /// True when the sequence was produced by concatenation rather than an
    /// overlap merge — i.e. via `--just-concatenate`, or via `--rescue-unmerged`
    /// for a pair that failed the overlap criteria.
    pub concatenated: bool,
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
/// Count matches, mismatches and indels in the overlap, as R's `C_eval_pair`
/// (`evaluate.cpp`) does.
///
/// The overlap runs from the column where **both** strands have begun (the
/// later of their first bases) to the column where the first of them ends.
/// A gap inside that span is an indel, including one at its edge: when an
/// aligned base faces a gap at the first or last overlap column, R counts it,
/// so a pair that is otherwise a perfect overlap fails `maxMismatch = 0`.
/// Starting at the first column where both have a base would skip it (#272).
///
/// Returns `(nmatch, nmismatch, nindel, ov_left, ov_right)`, or `None` when
/// the strands do not overlap.
fn analyze_overlap(al0: &[u8], al1: &[u8]) -> Option<(u32, u32, u32, usize, usize)> {
    let first = |al: &[u8]| al.iter().position(|&c| c != b'-');
    let last = |al: &[u8]| al.iter().rposition(|&c| c != b'-');
    let left = first(al0)?.max(first(al1)?);
    let right = last(al0)?.min(last(al1)?);
    if left > right {
        return None;
    }

    let mut nmatch = 0u32;
    let mut nmismatch = 0u32;
    let mut nindel = 0u32;

    for i in left..=right {
        match (al0[i] == b'-', al1[i] == b'-') {
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
            (true, false) if !trim_overhang => {
                result.push(nt_decode(al1[i]));
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
            (false, true) if !trim_overhang => {
                result.push(nt_decode(al0[i]));
            }
            _ => {}
        }
    }

    String::from_utf8_lossy(&result).into_owned()
}

/// Default-on provenance check: warn (to stderr) when the FASTQ a dada JSON
/// records having been computed from does not match the FASTQ now being passed
/// for that orientation. A mismatch usually means the four positional file
/// lists have drifted out of alignment (e.g. a glob expanded to a different
/// set), which would silently merge the wrong samples. This only warns.
///
/// Compares against `source_fastq`. An older dada JSON lacks it, and its
/// `input_file` is used only when that is not a derep JSON, which can never
/// match a FASTQ name (#111); otherwise the check is skipped.
fn warn_on_input_mismatch(label: &str, dada: &DadaJsonInput, fastq_path: &Path, dada_path: &Path) {
    let recorded = dada.source_fastq.as_deref().or_else(|| {
        dada.input_file
            .as_deref()
            .filter(|f| !f.ends_with(".json") && !f.ends_with(".json.gz"))
    });
    let Some(recorded) = recorded else { return };
    let passed = fastq_path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default();
    if recorded != passed {
        eprintln!(
            "[merge-pairs] warning: {label} dada '{}' was computed from '{recorded}', \
             but the {label} FASTQ passed is '{passed}' — check that the file lists line up",
            dada_path.display(),
        );
    }
}

// ---------------------------------------------------------------------------
// Sample-ID sanity check
// ---------------------------------------------------------------------------

fn check_sample_ids(
    sample_name: &str,
    fwd_sample: Option<&str>,
    rev_sample: Option<&str>,
    fwd_dada_path: &Path,
    rev_dada_path: &Path,
    fwd_fastq_path: &Path,
    rev_fastq_path: &Path,
) -> io::Result<()> {
    let mismatch = |msg: String| io::Error::new(io::ErrorKind::InvalidData, msg);

    match (fwd_sample, rev_sample) {
        (Some(f), Some(r)) if f != r => {
            return Err(mismatch(format!(
                "sample-id check: forward dada '{}' has sample '{f}' but reverse dada '{}' has sample '{r}'",
                fwd_dada_path.display(),
                rev_dada_path.display(),
            )));
        }
        (Some(f), _) if f != sample_name => {
            return Err(mismatch(format!(
                "sample-id check: forward dada '{}' has sample '{f}' but resolved sample name is '{sample_name}'",
                fwd_dada_path.display(),
            )));
        }
        (_, Some(r)) if r != sample_name => {
            return Err(mismatch(format!(
                "sample-id check: reverse dada '{}' has sample '{r}' but resolved sample name is '{sample_name}'",
                rev_dada_path.display(),
            )));
        }
        _ => {}
    }

    for (label, path) in [("forward", fwd_fastq_path), ("reverse", rev_fastq_path)] {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default();
        if !name.contains(sample_name) {
            return Err(mismatch(format!(
                "sample-id check: {label} FASTQ '{}' does not contain sample name '{sample_name}'",
                path.display(),
            )));
        }
    }

    Ok(())
}

/// `((fwd_asv, rev_asv), reads)`.
type PairCount = ((usize, usize), u64);

/// Count reads per `(fwd_asv, rev_asv)` pair, skipping reads unassigned in
/// either direction. Returns the pairs in first-occurrence order over reads,
/// as R's `unique(pairdf)` does, so abundance ties sort in R's order rather
/// than `HashMap` order (#240), plus the number of reads counted.
fn count_pairs(
    fwd_reads: &[usize],
    rev_reads: &[usize],
    fwd_map: &[Option<usize>],
    rev_map: &[Option<usize>],
) -> (Vec<PairCount>, u64) {
    let mut index: HashMap<(usize, usize), usize> = HashMap::new();
    let mut pairs: Vec<PairCount> = Vec::new();
    let mut total: u64 = 0;
    for (&fu, &ru) in fwd_reads.iter().zip(rev_reads) {
        // fwd_map[fu] is the ASV index for the fu-th forward unique (None = unassigned).
        let (Some(fa), Some(ra)) = (
            fwd_map.get(fu).copied().flatten(),
            rev_map.get(ru).copied().flatten(),
        ) else {
            continue;
        };
        let k = *index.entry((fa, ra)).or_insert_with(|| {
            pairs.push(((fa, ra), 0));
            pairs.len() - 1
        });
        pairs[k].1 += 1;
        total += 1;
    }
    (pairs, total)
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
    // Accept "dada" (independent), "dada-pooled", and "dada-pseudo" — same schema.
    let fwd_dada: DadaJsonInput =
        crate::misc::read_tagged_json(fwd_dada_path, &["dada", "dada-pooled", "dada-pseudo"])
            .with_path(fwd_dada_path)?;
    let rev_dada: DadaJsonInput =
        crate::misc::read_tagged_json(rev_dada_path, &["dada", "dada-pooled", "dada-pseudo"])
            .with_path(rev_dada_path)?;

    // Provenance warning (always on): does each dada JSON's recorded source
    // FASTQ match the FASTQ now being passed for that orientation?
    warn_on_input_mismatch("forward", &fwd_dada, fwd_fastq_path, fwd_dada_path);
    warn_on_input_mismatch("reverse", &rev_dada, rev_fastq_path, rev_dada_path);

    if params.check_sample_ids {
        check_sample_ids(
            sample_name,
            fwd_dada.sample.as_deref(),
            rev_dada.sample.as_deref(),
            fwd_dada_path,
            rev_dada_path,
            fwd_fastq_path,
            rev_fastq_path,
        )?;
    }

    let fwd_map = fwd_dada.map.as_ref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: 'map' field is absent — re-run `dada` with the current dada2-rs",
                fwd_dada_path.display()
            ),
        )
    })?;
    let rev_map = rev_dada.map.as_ref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: 'map' field is absent — re-run `dada` with the current dada2-rs",
                rev_dada_path.display()
            ),
        )
    })?;

    // ---- Re-dereplicate FASTQs ----
    // The ordering of unique sequences is deterministic (sorted by abundance
    // descending; ties preserve first-seen order) regardless of thread
    // count, so the indices here match those in the dada JSON map.
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
    let (pair_counts, total_pairs) = count_pairs(&fwd_derep.map, &rev_derep.map, fwd_map, rev_map);

    // ---- Attempt merge for each distinct pair ----
    let mut merged: Vec<MergedPair> = Vec::with_capacity(pair_counts.len());
    let mut accepted_pairs: u64 = 0;
    let mut align_buf = AlignBuffers::new();
    let spacer = "N".repeat(params.concat_nnn_len);

    // Build a concatenated (accepted) pair: fwd + N-spacer + RC(rev). Used by
    // both `--just-concatenate` and `--rescue-unmerged`.
    let make_concat = |fwd_seq: &str, rc_rev: &str, fi: usize, ri: usize, count: u64| MergedPair {
        sequence: format!("{fwd_seq}{spacer}{rc_rev}"),
        abundance: count,
        forward: fi,
        reverse: ri,
        nmatch: 0,
        nmismatch: 0,
        nindel: 0,
        accept: true,
        concatenated: true,
    };

    for ((fi, ri), count) in &pair_counts {
        let fwd_seq = &fwd_dada.asvs[*fi].sequence;
        let rev_seq = &rev_dada.asvs[*ri].sequence;
        let rc_rev = reverse_complement(rev_seq);

        // Just-concatenate mode: no alignment required.
        if params.just_concatenate {
            merged.push(make_concat(fwd_seq, &rc_rev, *fi, *ri, *count));
            accepted_pairs += count;
            continue;
        }

        let ov = align_pair(fwd_seq, &rc_rev, params.max_mismatch, &mut align_buf);

        let (nmatch, nmismatch, nindel, ov_left, ov_right) = match ov {
            Some(v) => v,
            None => {
                // No overlap at all.
                if params.rescue_unmerged {
                    merged.push(make_concat(fwd_seq, &rc_rev, *fi, *ri, *count));
                    accepted_pairs += count;
                } else if params.return_rejects {
                    merged.push(MergedPair {
                        sequence: String::new(),
                        abundance: *count,
                        forward: *fi,
                        reverse: *ri,
                        nmatch: 0,
                        nmismatch: 0,
                        nindel: 0,
                        accept: false,
                        concatenated: false,
                    });
                }
                continue;
            }
        };

        let accept = accepts(nmatch, nmismatch, nindel, params);

        // Rescue pairs that overlapped but failed the acceptance criteria by
        // concatenating them (takes precedence over return_rejects).
        if !accept && params.rescue_unmerged {
            merged.push(make_concat(fwd_seq, &rc_rev, *fi, *ri, *count));
            accepted_pairs += count;
            continue;
        }

        if !accept && !params.return_rejects {
            continue;
        }

        let sequence = if accept {
            accepted_pairs += count;
            build_merged(
                &align_buf.al0,
                &align_buf.al1,
                ov_left,
                ov_right,
                params.trim_overhang,
                true,
            )
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
            concatenated: false,
        });
    }

    // Abundance descending. Stable, like R's `order()`, so ties keep
    // first-occurrence order (#240).
    merged.sort_by_key(|a| std::cmp::Reverse(a.abundance));

    let num_merged = merged.iter().filter(|m| m.accept).count();

    Ok(SampleMergeResult {
        sample: sample_name.to_string(),
        total_pairs,
        accepted_pairs,
        num_merged,
        merged,
    })
}

#[cfg(test)]
mod count_pairs_tests {
    use super::count_pairs;

    /// Fifty pairs of one read each: every abundance ties, so the order is
    /// entirely the tie-break. It must be first occurrence over reads, as R's
    /// `unique(pairdf)` gives; `HashMap` order would match by chance with
    /// negligible probability (#240).
    #[test]
    fn tied_pairs_keep_first_occurrence_order() {
        let n = 50;
        // Read i pairs forward unique i with reverse unique n-1-i.
        let fwd: Vec<usize> = (0..n).collect();
        let rev: Vec<usize> = (0..n).rev().collect();
        let map: Vec<Option<usize>> = (0..n).map(Some).collect();
        let (pairs, total) = count_pairs(&fwd, &rev, &map, &map);
        let want: Vec<((usize, usize), u64)> = (0..n).map(|i| ((i, n - 1 - i), 1)).collect();
        assert_eq!(pairs, want);
        assert_eq!(total, n as u64);
    }

    #[test]
    fn repeats_accumulate_and_unassigned_reads_are_skipped() {
        let fwd = [0, 1, 0, 2, 0];
        let rev = [0, 1, 0, 0, 1];
        let fmap = [Some(0), Some(1), None];
        let rmap = [Some(0), Some(1)];
        let (pairs, total) = count_pairs(&fwd, &rev, &fmap, &rmap);
        assert_eq!(pairs, vec![((0, 0), 2), ((1, 1), 1), ((0, 1), 1)]);
        assert_eq!(total, 4);
    }
}

#[cfg(test)]
mod merge_scoring_tests {
    use super::*;

    /// An amplicon whose forward and reverse reads share a 16-base perfect
    /// overlap, plus a 60-column overlap at another offset with 4 mismatches
    /// (a 44 bp near-repeat), the shape length-variable ITS2 produces (#272).
    /// Returns (amplicon, forward read, RC(reverse) read).
    fn repeat_amplicon() -> (String, String, String) {
        let n = 150; // forward read length
        let mut state: u64 = 0x2720_2720;
        let mut base = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 33) as usize % 4]
        };
        let mut a: Vec<u8> = (0..n - 16).map(|_| base()).collect();
        let head = a[n - 60..n - 44].to_vec();
        a.extend_from_slice(&head); // a[n-16..n] = a[n-60..n-44]
        let period = a[n - 44..n].to_vec();
        a.extend_from_slice(&period); // a[n..n+44] = a[n-44..n] ...
        for i in [n + 5, n + 15, n + 25, n + 35] {
            a[i] = if a[i] == b'A' { b'C' } else { b'A' }; // ... with 4 differences
        }
        a.extend((0..40).map(|_| base()));
        let amplicon = String::from_utf8(a).unwrap();
        let fwd = amplicon[..n].to_string();
        let rc_rev = amplicon[n - 16..].to_string();
        (amplicon, fwd, rc_rev)
    }

    fn params(max_mismatch: u32) -> MergeParams {
        MergeParams {
            min_overlap: 12,
            max_mismatch,
            return_rejects: false,
            rescue_unmerged: false,
            trim_overhang: false,
            just_concatenate: false,
            concat_nnn_len: 10,
            phred_offset: 33,
            check_sample_ids: true,
            verbose: false,
        }
    }

    /// With the denoising scores (5 / −4 / −8) the longer imperfect overlap
    /// wins, so the pair would be rejected: the pre-#272 behaviour.
    #[test]
    fn denoising_scores_prefer_the_long_imperfect_overlap() {
        let (_, fwd, rc_rev) = repeat_amplicon();
        let mut buf = AlignBuffers::new();
        align_endsfree_with_buf(
            &intstr(fwd.as_bytes()),
            &intstr(rc_rev.as_bytes()),
            5,
            -4,
            -8,
            -1,
            &mut buf,
        );
        let (nmatch, nmismatch, nindel, _, _) = analyze_overlap(&buf.al0, &buf.al1).unwrap();
        assert_eq!((nmatch, nmismatch, nindel), (56, 4, 0));
        assert!(!accepts(nmatch, nmismatch, nindel, &params(0)));
    }

    /// R's `mergePairs` scores find the 16-base perfect overlap, accept it, and
    /// rebuild the amplicon exactly.
    #[test]
    fn merge_scores_find_the_perfect_overlap_and_merge() {
        let (amplicon, fwd, rc_rev) = repeat_amplicon();
        let mut buf = AlignBuffers::new();
        let (nmatch, nmismatch, nindel, left, right) =
            align_pair(&fwd, &rc_rev, 0, &mut buf).unwrap();
        assert_eq!((nmatch, nmismatch, nindel), (16, 0, 0));
        assert!(accepts(nmatch, nmismatch, nindel, &params(0)));
        let merged = build_merged(&buf.al0, &buf.al1, left, right, false, true);
        assert_eq!(merged, amplicon);
    }

    /// R's `C_eval_pair` starts the overlap where both strands have begun, so
    /// a base facing a gap at that first column is an indel. This is the shape
    /// of a real MiSeq V4 pair: the reverse read's first base is gapped against
    /// the forward read (a gap and a mismatch both cost −64), and R rejects it.
    /// Starting at the first column where both have a base would hide it.
    #[test]
    fn overlap_starts_where_both_strands_have_begun() {
        let fwd = b"CCCCC-GGACT----";
        let rcrev = b"-----TGGACTAAAA";
        let (nmatch, nmismatch, nindel, left, right) = analyze_overlap(fwd, rcrev).unwrap();
        assert_eq!((left, right), (5, 10));
        assert_eq!((nmatch, nmismatch, nindel), (5, 0, 1));
        assert!(!accepts(nmatch, nmismatch, nindel, &params(0)));
    }

    /// R's acceptance rule: `minOverlap` counts matching bases only, and indels
    /// count against `maxMismatch` rather than being refused outright.
    #[test]
    fn acceptance_follows_r() {
        assert!(!accepts(11, 0, 0, &params(0)), "11 matches < minOverlap 12");
        assert!(accepts(12, 0, 0, &params(0)));
        assert!(!accepts(11, 1, 0, &params(1)), "a mismatch is not a match");
        assert!(
            accepts(20, 0, 1, &params(1)),
            "one indel within maxMismatch 1"
        );
        assert!(!accepts(20, 1, 1, &params(1)), "mismatch + indel = 2 > 1");
    }
}
