//! Remove primer sequences from FASTQ reads, mirroring R's `removePrimers()`.
//!
//! Primers are matched using a sliding-window Hamming search with IUPAC
//! ambiguity support.  Orientation detection (orient=TRUE) reverse-complements
//! reads that match primers only in the RC direction before trimming.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use noodles::bgzf;
use noodles::fastq;
use rayon::prelude::*;

use crate::misc::WithPath;

const CHUNK_SIZE: usize = 4096;

// ---------------------------------------------------------------------------
// IUPAC bitmask table
// ---------------------------------------------------------------------------

// Standard IUPAC nucleotide bitmask encoding: A=1, C=2, G=4, T/U=8.
// Ambiguity codes are bitwise OR of their constituent nucleotides.
// This convention matches the one used by the `bio` crate
// (https://docs.rs/bio, src/alphabets/dna.rs).
const fn make_iupac_table() -> [u8; 128] {
    let mut t = [0u8; 128];
    t[b'A' as usize] = 0b0001;
    t[b'a' as usize] = 0b0001;
    t[b'C' as usize] = 0b0010;
    t[b'c' as usize] = 0b0010;
    t[b'G' as usize] = 0b0100;
    t[b'g' as usize] = 0b0100;
    t[b'T' as usize] = 0b1000;
    t[b't' as usize] = 0b1000;
    t[b'U' as usize] = 0b1000;
    t[b'u' as usize] = 0b1000;
    // Ambiguity codes
    t[b'R' as usize] = 0b0101;
    t[b'r' as usize] = 0b0101; // A|G
    t[b'Y' as usize] = 0b1010;
    t[b'y' as usize] = 0b1010; // C|T
    t[b'S' as usize] = 0b0110;
    t[b's' as usize] = 0b0110; // G|C
    t[b'W' as usize] = 0b1001;
    t[b'w' as usize] = 0b1001; // A|T
    t[b'K' as usize] = 0b1100;
    t[b'k' as usize] = 0b1100; // G|T
    t[b'M' as usize] = 0b0011;
    t[b'm' as usize] = 0b0011; // A|C
    t[b'B' as usize] = 0b1110;
    t[b'b' as usize] = 0b1110; // C|G|T
    t[b'D' as usize] = 0b1101;
    t[b'd' as usize] = 0b1101; // A|G|T
    t[b'H' as usize] = 0b1011;
    t[b'h' as usize] = 0b1011; // A|C|T
    t[b'V' as usize] = 0b0111;
    t[b'v' as usize] = 0b0111; // A|C|G
    t[b'N' as usize] = 0b1111;
    t[b'n' as usize] = 0b1111; // A|C|G|T
    t
}

const IUPAC_BITS: [u8; 128] = make_iupac_table();

#[inline]
fn iupac_matches(pattern_base: u8, seq_base: u8) -> bool {
    let p = pattern_base as usize;
    let s = seq_base as usize;
    p < 128 && s < 128 && IUPAC_BITS[p] & IUPAC_BITS[s] != 0
}

// ---------------------------------------------------------------------------
// Primer search
// ---------------------------------------------------------------------------

/// Count Hamming mismatches between `primer` and a same-length window of
/// `read`, using IUPAC matching.
fn count_mismatches(primer: &[u8], window: &[u8]) -> usize {
    primer
        .iter()
        .zip(window.iter())
        .filter(|&(&p, &s)| !iupac_matches(p, s))
        .count()
}

/// Return the `(start, end)` (end exclusive) of the **first** window in `read`
/// where `primer` matches with ≤ `max_mismatch` mismatches.
///
/// Mirrors R's `vmatchPattern(fixed=FALSE, max.mismatch=N)` boundary behaviour:
/// the primer may overhang either read boundary by up to `max_mismatch` bases,
/// with each out-of-bounds base counted as one mismatch.
fn find_first_match(primer: &[u8], read: &[u8], max_mismatch: usize) -> Option<(usize, usize)> {
    let plen = primer.len();
    let rlen = read.len();
    if plen == 0 || rlen == 0 {
        return None;
    }
    let i_min = -(max_mismatch as isize);
    let i_max = rlen as isize - 1 + max_mismatch as isize;
    for raw_i in i_min..=i_max {
        let oob_left = if raw_i < 0 { (-raw_i) as usize } else { 0 };
        let right_end = raw_i + plen as isize;
        let oob_right = if right_end > rlen as isize {
            (right_end - rlen as isize) as usize
        } else {
            0
        };
        if oob_left + oob_right > max_mismatch {
            continue;
        }
        let read_start = raw_i.max(0) as usize;
        let read_end = right_end.min(rlen as isize) as usize;
        if read_start >= read_end {
            continue;
        }
        let inner_mm = count_mismatches(
            &primer[oob_left..plen - oob_right],
            &read[read_start..read_end],
        );
        if oob_left + inner_mm + oob_right <= max_mismatch {
            return Some((read_start, read_end));
        }
    }
    None
}

/// Return the `(start, end)` (end exclusive) of the **last** window in `read`
/// where `primer` matches with ≤ `max_mismatch` mismatches.
///
/// Mirrors R's behaviour for reverse primers: when multiple matches exist, use
/// the rightmost one so the longest possible internal fragment is retained.
/// Also supports boundary overhang (see [`find_first_match`]).
fn find_last_match(primer: &[u8], read: &[u8], max_mismatch: usize) -> Option<(usize, usize)> {
    let plen = primer.len();
    let rlen = read.len();
    if plen == 0 || rlen == 0 {
        return None;
    }
    let i_min = -(max_mismatch as isize);
    let i_max = rlen as isize - 1 + max_mismatch as isize;
    for raw_i in (i_min..=i_max).rev() {
        let oob_left = if raw_i < 0 { (-raw_i) as usize } else { 0 };
        let right_end = raw_i + plen as isize;
        let oob_right = if right_end > rlen as isize {
            (right_end - rlen as isize) as usize
        } else {
            0
        };
        if oob_left + oob_right > max_mismatch {
            continue;
        }
        let read_start = raw_i.max(0) as usize;
        let read_end = right_end.min(rlen as isize) as usize;
        if read_start >= read_end {
            continue;
        }
        let inner_mm = count_mismatches(
            &primer[oob_left..plen - oob_right],
            &read[read_start..read_end],
        );
        if oob_left + inner_mm + oob_right <= max_mismatch {
            return Some((read_start, read_end));
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Semi-global alignment for indel-tolerant primer matching
//
// Mirrors R Biostrings `matchPattern(..., with.indels = TRUE)`, which uses
// Levenshtein edit distance (mismatches and indels each cost 1).
//
// The alignment is semi-global: `primer` is the query (fully consumed, no free
// end gaps) and `read` is the target (free to start / end anywhere).
//
// dp[i][j] = minimum edit distance to align primer[0..i] ending at read[..j].
//
// Initialisation:
//   dp[0][j] = 0   — primer may start at any position in the read
//   dp[i][0] = i   — i deletions needed when no read bases are available
// ---------------------------------------------------------------------------

fn semi_global_dp(primer: &[u8], read: &[u8]) -> Vec<Vec<usize>> {
    let m = primer.len();
    let n = read.len();
    let mut dp = vec![vec![0usize; n + 1]; m + 1];
    for (i, row) in dp.iter_mut().enumerate().skip(1) {
        row[0] = i;
    }
    #[allow(clippy::needless_range_loop)]
    for i in 1..=m {
        for j in 1..=n {
            let sub = if iupac_matches(primer[i - 1], read[j - 1]) {
                0
            } else {
                1
            };
            dp[i][j] = (dp[i - 1][j - 1] + sub)
                .min(dp[i - 1][j] + 1) // deletion from read (gap in read)
                .min(dp[i][j - 1] + 1); // insertion into read (gap in primer)
        }
    }
    dp
}

/// Trace back through `dp` from `(primer.len(), j_end)` and return the read
/// index where the match begins.
fn backtrack_start(dp: &[Vec<usize>], primer: &[u8], read: &[u8], j_end: usize) -> usize {
    let mut i = primer.len();
    let mut j = j_end;
    while i > 0 {
        let d = dp[i][j];
        if j > 0 {
            let sub = if iupac_matches(primer[i - 1], read[j - 1]) {
                0
            } else {
                1
            };
            if dp[i - 1][j - 1] + sub == d {
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if dp[i - 1][j] + 1 == d {
            i -= 1; // gap in read
        } else {
            j -= 1; // gap in primer
        }
    }
    j
}

/// Indel-tolerant version of [`find_first_match`].
///
/// Returns the `(start, end)` pair with the **smallest start position** among
/// all matches with edit distance ≤ `max_mismatch`, breaking ties by smallest
/// end.  Matches R's first-IRanges-element selection for the forward primer.
fn find_first_match_indels(
    primer: &[u8],
    read: &[u8],
    max_mismatch: usize,
) -> Option<(usize, usize)> {
    if primer.is_empty() || read.is_empty() {
        return None;
    }
    let dp = semi_global_dp(primer, read);
    let m = primer.len();
    let n = read.len();
    (1..=n)
        .filter(|&j| dp[m][j] <= max_mismatch)
        .map(|j| (backtrack_start(&dp, primer, read, j), j))
        .min_by_key(|&(start, end)| (start, end))
}

/// Indel-tolerant version of [`find_last_match`].
///
/// Returns the `(start, end)` pair with the **largest start position** among
/// all matches with edit distance ≤ `max_mismatch`, breaking ties by largest
/// end.  Matches R's last-IRanges-element selection for the reverse primer.
fn find_last_match_indels(
    primer: &[u8],
    read: &[u8],
    max_mismatch: usize,
) -> Option<(usize, usize)> {
    if primer.is_empty() || read.is_empty() {
        return None;
    }
    let dp = semi_global_dp(primer, read);
    let m = primer.len();
    let n = read.len();
    (1..=n)
        .filter(|&j| dp[m][j] <= max_mismatch)
        .map(|j| (backtrack_start(&dp, primer, read, j), j))
        .max_by_key(|&(start, end)| (start, end))
}

// ---------------------------------------------------------------------------
// Reverse complement
// ---------------------------------------------------------------------------

/// RC for read sequences (ACGT + N only). Non-ACGT bases become N.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

/// RC for primer sequences with full IUPAC ambiguity support.
///
/// Mirrors R DADA2's `rc()` helper.  Output is upper-case.
/// Self-complementary codes: S (G|C), W (A|T), N.
pub fn iupac_reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            b'R' | b'r' => b'Y', // A|G -> C|T
            b'Y' | b'y' => b'R', // C|T -> A|G
            b'S' | b's' => b'S', // G|C -> G|C
            b'W' | b'w' => b'W', // A|T -> A|T
            b'K' | b'k' => b'M', // G|T -> A|C
            b'M' | b'm' => b'K', // A|C -> G|T
            b'B' | b'b' => b'V', // C|G|T -> A|C|G
            b'V' | b'v' => b'B', // A|C|G -> C|G|T
            b'D' | b'd' => b'H', // A|G|T -> A|C|T
            b'H' | b'h' => b'D', // A|C|T -> A|G|T
            b'N' | b'n' => b'N',
            _ => b'N',
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

pub struct RemovePrimersParams {
    pub primer_fwd: Vec<u8>,
    pub primer_rev: Option<Vec<u8>>,
    pub max_mismatch: usize,
    pub allow_indels: bool,
    pub trim_fwd: bool,
    pub trim_rev: bool,
    pub orient: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct PrimerStats {
    pub reads_in: u64,
    pub reads_out: u64,
    pub reads_reoriented: u64,
}

// ---------------------------------------------------------------------------
// Per-read processing
// ---------------------------------------------------------------------------

/// Apply primer detection and trimming to one read.
///
/// Returns `Some((trimmed_seq, trimmed_qual, reoriented))` when the read
/// passes, where `reoriented` is `true` if the read was reverse-complemented
/// before trimming.  Returns `None` when the read should be discarded.
fn process_read(
    seq: &[u8],
    qual: &[u8],
    params: &RemovePrimersParams,
) -> Option<(Vec<u8>, Vec<u8>, bool)> {
    let mm = params.max_mismatch;
    let find_fwd = |p: &[u8], r: &[u8]| {
        if params.allow_indels {
            find_first_match_indels(p, r, mm)
        } else {
            find_first_match(p, r, mm)
        }
    };
    let find_rev = |p: &[u8], r: &[u8]| {
        if params.allow_indels {
            find_last_match_indels(p, r, mm)
        } else {
            find_last_match(p, r, mm)
        }
    };

    let fwd_hit = find_fwd(&params.primer_fwd, seq);

    // Resolve working orientation: forward if fwd primer found, or RC if orient
    // is enabled and the fwd primer matches in the reverse complement.
    let (work_seq, work_qual, fwd_hit, reoriented) = match (fwd_hit, params.orient) {
        (Some(hit), _) => (seq.to_vec(), qual.to_vec(), hit, false),
        (None, true) => {
            let rc_seq = reverse_complement(seq);
            let rc_qual: Vec<u8> = qual.iter().copied().rev().collect();
            let rc_hit = find_fwd(&params.primer_fwd, &rc_seq)?;
            (rc_seq, rc_qual, rc_hit, true)
        }
        (None, false) => return None,
    };

    // Reverse primer: require last match in the (possibly flipped) read.
    let rev_hit = params
        .primer_rev
        .as_ref()
        .map(|rev| find_rev(rev, &work_seq));

    if let Some(None) = rev_hit {
        return None; // rev primer provided but not found
    }
    let rev_hit = rev_hit.flatten();

    // Trim boundaries (0-based, exclusive end).
    // Mirrors R: first = end(match.fwd) + 1 (1-based) → fwd_hit.1 (0-based)
    //            last  = start(match.rev) - 1 (1-based) → rev_hit.0 (0-based)
    let trim_start = if params.trim_fwd { fwd_hit.1 } else { 0 };
    let trim_end = match (params.trim_rev, rev_hit) {
        (true, Some((rev_start, _))) => rev_start,
        _ => work_seq.len(),
    };

    if trim_end <= trim_start {
        return None;
    }

    Some((
        work_seq[trim_start..trim_end].to_vec(),
        work_qual[trim_start..trim_end].to_vec(),
        reoriented,
    ))
}

// ---------------------------------------------------------------------------
// I/O helpers (mirrors filter_trim.rs)
// ---------------------------------------------------------------------------

fn is_gz(path: &Path) -> bool {
    path.extension().and_then(|e| e.to_str()) == Some("gz")
}

type FastqReader = fastq::io::Reader<Box<dyn BufRead>>;

fn open_reader(path: &Path) -> io::Result<FastqReader> {
    let file = File::open(path).with_path(path)?;
    let inner: Box<dyn BufRead> = if is_gz(path) {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    Ok(fastq::io::Reader::new(inner))
}

fn open_writer(path: &Path, compress: bool, threads: usize) -> io::Result<Box<dyn Write>> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let file = File::create(path).with_path(path)?;
    if compress {
        if threads > 1 {
            let w = bgzf::multithreaded_writer::Builder::default()
                .set_worker_count(NonZeroUsize::new(threads).unwrap())
                .build_from_writer(file);
            Ok(Box::new(w))
        } else {
            Ok(Box::new(GzEncoder::new(file, Compression::default())))
        }
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn write_record(
    out: &mut dyn Write,
    name: &[u8],
    desc: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    out.write_all(b"@")?;
    out.write_all(name)?;
    if !desc.is_empty() {
        out.write_all(b" ")?;
        out.write_all(desc)?;
    }
    out.write_all(b"\n")?;
    out.write_all(seq)?;
    out.write_all(b"\n+\n")?;
    out.write_all(qual)?;
    out.write_all(b"\n")?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn remove_primers(
    input: &Path,
    output: &Path,
    params: &RemovePrimersParams,
    compress: bool,
    threads: usize,
    verbose: bool,
) -> io::Result<PrimerStats> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    let mut reader = open_reader(input)?;
    let mut writer = open_writer(output, compress, threads)?;

    let mut reads_in: u64 = 0;
    let mut reads_out: u64 = 0;
    let mut reads_reoriented: u64 = 0;
    let mut record = fastq::Record::default();

    type ReadRecord = (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);
    // Batch of (name, description, seq, qual) owned buffers.
    let mut batch: Vec<ReadRecord> = Vec::with_capacity(CHUNK_SIZE);

    loop {
        // Fill one chunk from the reader (single-threaded I/O).
        batch.clear();
        loop {
            match reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    batch.push((
                        record.name().to_vec(),
                        record.description().to_vec(),
                        record.sequence().to_vec(),
                        record.quality_scores().to_vec(),
                    ));
                    if batch.len() == CHUNK_SIZE {
                        break;
                    }
                }
                Err(e) => return Err(e),
            }
        }
        if batch.is_empty() {
            break;
        }
        reads_in += batch.len() as u64;

        type ProcessResult = Option<(Vec<u8>, Vec<u8>, bool)>;
        // Process all reads in the chunk in parallel using the local pool.
        let results: Vec<ProcessResult> = pool.install(|| {
            batch
                .par_iter()
                .map(|(_, _, seq, qual)| process_read(seq, qual, params))
                .collect()
        });

        // Write passing reads in original order (single-threaded I/O).
        for ((name, desc, _, _), result) in batch.iter().zip(results.iter()) {
            if let Some((seq_out, qual_out, reoriented)) = result {
                write_record(&mut *writer, name, desc, seq_out, qual_out)?;
                reads_out += 1;
                if *reoriented {
                    reads_reoriented += 1;
                }
            }
        }
    }

    writer.flush()?;

    if reads_out == 0 {
        drop(writer);
        let _ = std::fs::remove_file(output);
        if verbose {
            eprintln!(
                "[remove-primers] {} — no reads passed primer detection, output not written",
                input.display()
            );
        }
    } else if verbose {
        let pct = reads_out as f64 * 100.0 / reads_in as f64;
        eprintln!(
            "[remove-primers] {} → {} reads in, {} out ({:.1}%), {} reoriented",
            input.display(),
            reads_in,
            reads_out,
            pct,
            reads_reoriented,
        );
    }

    Ok(PrimerStats {
        reads_in,
        reads_out,
        reads_reoriented,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_exact() {
        assert!(iupac_matches(b'A', b'A'));
        assert!(iupac_matches(b'C', b'C'));
        assert!(!iupac_matches(b'A', b'C'));
    }

    #[test]
    fn iupac_ambig() {
        assert!(iupac_matches(b'R', b'A')); // R = A|G
        assert!(iupac_matches(b'R', b'G'));
        assert!(!iupac_matches(b'R', b'C'));
        assert!(iupac_matches(b'Y', b'C')); // Y = C|T
        assert!(iupac_matches(b'Y', b'T'));
        assert!(!iupac_matches(b'Y', b'A'));
        assert!(iupac_matches(b'N', b'T')); // N matches anything
    }

    #[test]
    fn find_first_left_boundary_overhang() {
        // Primer "AACGT" (5bp), read "CGTAAA" (6bp), max_mismatch=2.
        // raw_i=-2: oob_left=2, primer[2..5]="CGT" vs read[0..3]="CGT" (0 mm).
        // Total = 2+0+0 = 2 ≤ 2.  First valid start → (0, 3).
        let primer = b"AACGT";
        let read = b"CGTAAA";
        assert_eq!(find_first_match(primer, read, 2), Some((0, 3)));
        // With max_mismatch=1 the overhang budget is only 1, so raw_i=-2 is skipped.
        // raw_i=-1: oob_left=1, primer[1..5]="ACGT" vs read[0..4]="CGTA" — mismatches C≠A,G=G,T=T,A≠A? "CGTA" vs "ACGT": C≠A, G=G, T=T, A≠A = 2 inner mm → total 1+2=3 > 1. Skip.
        // No earlier match → None (or check deeper into the read).
        // Just verify the 2-budget case works.
        assert_eq!(find_first_match(primer, read, 0), None);
    }

    #[test]
    fn find_last_right_boundary_overhang() {
        // Primer "CGTAA" (5bp), read "AAACGT" (6bp), max_mismatch=2.
        // raw_i=4: right_end=9, oob_right=3 > 2. Skip.
        // raw_i=3: right_end=8, oob_right=2. primer[0..3]="CGT" vs read[3..6]="CGT" = 0 mm.
        // Total = 0+0+2 = 2 ≤ 2. Last (rightmost) → (3, 6).
        let primer = b"CGTAA";
        let read = b"AAACGT";
        assert_eq!(find_last_match(primer, read, 2), Some((3, 6)));
    }

    #[test]
    fn find_first_exact() {
        let read = b"AAACGTCGT";
        let primer = b"CGT";
        assert_eq!(find_first_match(primer, read, 0), Some((3, 6)));
    }

    #[test]
    fn find_last_exact() {
        let read = b"AAACGTCGT";
        let primer = b"CGT";
        assert_eq!(find_last_match(primer, read, 0), Some((6, 9)));
    }

    #[test]
    fn find_with_mismatch() {
        let read = b"AAACGACGT"; // CGT → CGA at pos 3, then CGT at pos 6
        let primer = b"CGT";
        assert_eq!(find_first_match(primer, read, 0), Some((6, 9)));
        assert_eq!(find_first_match(primer, read, 1), Some((3, 6)));
    }

    #[test]
    fn find_iupac_primer() {
        // R matches A or G
        let primer = b"CGR"; // matches CGA and CGG and CGG
        let read = b"AAACGACGT";
        assert_eq!(find_first_match(primer, read, 0), Some((3, 6))); // CGA matches CGR
    }

    #[test]
    fn process_read_basic_trim() {
        // fwd primer "AAA" at pos 0..3, rev primer "TTT" at pos 6..9
        // expected trimmed region: pos 3..6 = "CCC"
        let seq = b"AAACCCTTT";
        let qual = b"IIIIIIIII";
        let params = RemovePrimersParams {
            primer_fwd: b"AAA".to_vec(),
            primer_rev: Some(b"TTT".to_vec()),
            max_mismatch: 0,
            allow_indels: false,
            trim_fwd: true,
            trim_rev: true,
            orient: false,
        };
        let result = process_read(seq, qual, &params);
        assert_eq!(
            result.as_ref().map(|(s, _, _)| s.as_slice()),
            Some(b"CCC".as_slice())
        );
        assert_eq!(
            result.as_ref().map(|(_, q, _)| q.as_slice()),
            Some(b"III".as_slice())
        );
        assert_eq!(result.as_ref().map(|&(_, _, r)| r), Some(false));
    }

    #[test]
    fn process_read_orient_flip() {
        // fwd primer "AAA" only matches in RC of this read
        // read: RC("AAACCC") = "GGGTTT"
        let seq = b"GGGTTT";
        let qual = b"ABCDEF";
        let params = RemovePrimersParams {
            primer_fwd: b"AAA".to_vec(),
            primer_rev: None,
            max_mismatch: 0,
            allow_indels: false,
            trim_fwd: true,
            trim_rev: false,
            orient: true,
        };
        let result = process_read(seq, qual, &params);
        // RC of GGGTTT = AAACCC; after trimming fwd primer (AAA), seq = CCC
        assert_eq!(
            result.as_ref().map(|(s, _, _)| s.as_slice()),
            Some(b"CCC".as_slice())
        );
        assert_eq!(result.as_ref().map(|&(_, _, r)| r), Some(true));
    }

    #[test]
    fn process_read_no_primer_discarded() {
        let seq = b"TTTGGGCCC";
        let qual = b"IIIIIIIII";
        let params = RemovePrimersParams {
            primer_fwd: b"AAA".to_vec(),
            primer_rev: None,
            max_mismatch: 0,
            allow_indels: false,
            trim_fwd: true,
            trim_rev: false,
            orient: false,
        };
        assert!(process_read(seq, qual, &params).is_none());
    }

    // -------------------------------------------------------------------------
    // iupac_reverse_complement tests
    // -------------------------------------------------------------------------

    #[test]
    fn iupac_rc_acgt() {
        assert_eq!(iupac_reverse_complement(b"ACGT"), b"ACGT");
    }

    #[test]
    fn iupac_rc_ambig() {
        // Single-base complement pairs
        assert_eq!(iupac_reverse_complement(b"R"), b"Y"); // A|G -> C|T
        assert_eq!(iupac_reverse_complement(b"Y"), b"R"); // C|T -> A|G
        assert_eq!(iupac_reverse_complement(b"K"), b"M"); // G|T -> A|C
        assert_eq!(iupac_reverse_complement(b"M"), b"K"); // A|C -> G|T
        assert_eq!(iupac_reverse_complement(b"B"), b"V"); // C|G|T -> A|C|G
        assert_eq!(iupac_reverse_complement(b"V"), b"B"); // A|C|G -> C|G|T
        assert_eq!(iupac_reverse_complement(b"D"), b"H"); // A|G|T -> A|C|T
        assert_eq!(iupac_reverse_complement(b"H"), b"D"); // A|C|T -> A|G|T
        // Self-complementary codes
        assert_eq!(iupac_reverse_complement(b"S"), b"S"); // G|C
        assert_eq!(iupac_reverse_complement(b"W"), b"W"); // A|T
        assert_eq!(iupac_reverse_complement(b"N"), b"N");
        // Two-base: reverse + complement
        assert_eq!(iupac_reverse_complement(b"RR"), b"YY");
        assert_eq!(iupac_reverse_complement(b"RY"), b"RY"); // palindrome: rev(YR) -> complement -> RY
    }

    #[test]
    fn iupac_rc_lowercase() {
        assert_eq!(iupac_reverse_complement(b"acgt"), b"ACGT");
    }

    // -------------------------------------------------------------------------
    // Indel mode tests
    // -------------------------------------------------------------------------

    #[test]
    fn indel_exact_match_same_as_hamming() {
        // With no indels needed, both modes should agree.
        let primer = b"ACGT";
        let read = b"NNACGTNNN";
        assert_eq!(
            find_first_match(primer, read, 0),
            find_first_match_indels(primer, read, 0),
        );
    }

    #[test]
    fn indel_single_deletion_in_read() {
        // Primer: ACGT  Read: AGTAAA — primer with C deleted (edit dist 1).
        // With strict Hamming (max=0) no 4-base window is an exact match.
        // With indels (max=1) the primer aligns via one deletion.
        let primer = b"ACGT";
        let read = b"AGTAAA";
        assert_eq!(find_first_match(primer, read, 0), None);
        let hit = find_first_match_indels(primer, read, 1);
        assert!(hit.is_some(), "indel mode should find 1-deletion match");
        assert_eq!(hit.unwrap().0, 0); // match starts at read position 0
    }

    #[test]
    fn indel_single_insertion_in_read() {
        // Primer: ACT  Read: ACGTAA — read has an extra G (1 insertion, edit dist 1).
        // Hamming (max=0): no exact 3-base window.
        // Indel (max=1): ACT aligns ACGT with one insertion.
        let primer = b"ACT";
        let read = b"ACGTAA";
        assert_eq!(find_first_match(primer, read, 0), None);
        let hit = find_first_match_indels(primer, read, 1);
        assert!(hit.is_some(), "indel mode should find 1-insertion match");
    }

    #[test]
    fn indel_last_match() {
        // Two candidate positions; last_match_indels should return the rightmost start.
        let primer = b"CGT";
        let read = b"ACGTNNNCGT";
        let first = find_first_match_indels(primer, read, 0).unwrap();
        let last = find_last_match_indels(primer, read, 0).unwrap();
        assert!(first.0 <= last.0);
        assert_eq!(last, (7, 10));
    }

    #[test]
    fn process_read_indel_mode_trims() {
        // fwd primer "ACGT" with 1-deletion tolerance; primer appears as "AGT" at start
        let seq = b"AGTCCCNNN";
        let qual = b"IIIIIIIII";
        let params = RemovePrimersParams {
            primer_fwd: b"ACGT".to_vec(),
            primer_rev: None,
            max_mismatch: 1,
            allow_indels: true,
            trim_fwd: true,
            trim_rev: false,
            orient: false,
        };
        let result = process_read(seq, qual, &params);
        assert!(
            result.is_some(),
            "read should pass with 1-deletion indel tolerance"
        );
        // Trimmed region starts after the matched primer region
        let (trimmed_seq, _, _) = result.unwrap();
        assert!(
            !trimmed_seq.starts_with(b"AGT"),
            "primer region should be trimmed"
        );
    }
}
