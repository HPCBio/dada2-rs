//! Filter and trim FASTQ reads, mirroring R's `filterAndTrim` / `fastqFilter` /
//! `fastqPairedFilter` functions.
//!
//! Filtering is applied in this order (matching R):
//!   1. maxLen (on the raw, untrimmed read)
//!   2. trimLeft / trimRight
//!   3. truncQ  (truncate at first base with Phred ≤ truncQ)
//!   4. truncLen (discard if too short; truncate to length)
//!   5. minLen
//!   6. maxN
//!   7. minQ
//!   8. maxEE
//!   9. rm_phix
//!  10. rm_lowcomplex

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use noodles::fasta;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use noodles::fastq;

use crate::filter::match_ref;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Per-direction filtering parameters.
///
/// For paired data create two copies — one for forward, one for reverse.
pub struct FilterParams {
    /// Truncate reads at the first Phred score ≤ this value (default 2).
    pub trunc_q: u8,
    /// Truncate reads to this length; discard if shorter (0 = disabled).
    pub trunc_len: usize,
    /// Remove this many bases from the 5′ end.
    pub trim_left: usize,
    /// Remove this many bases from the 3′ end.
    pub trim_right: usize,
    /// Discard reads longer than this before trimming (0 = no limit).
    pub max_len: usize,
    /// Discard reads shorter than this after all trimming (default 20).
    pub min_len: usize,
    /// Discard reads with more than this many N bases (default 0).
    pub max_n: usize,
    /// Discard reads with any Phred score below this value (0 = disabled).
    pub min_q: u8,
    /// Discard reads with expected errors above this value (Inf = disabled).
    pub max_ee: f64,
    /// Pre-loaded phiX genome sequence (forward strand only); set to `Some` to
    /// enable phiX filtering.  The reverse complement is derived at runtime.
    pub phix_genome: Option<Vec<u8>>,
    /// Discard reads with 2-mer Shannon richness below this value (0.0 = disabled).
    pub rm_lowcomplex: f64,
    /// Phred quality offset (33 for Sanger/Illumina 1.8+).
    pub phred_offset: u8,
}

/// Per-sample filter statistics.
#[derive(Debug, Clone, Copy)]
pub struct SampleStats {
    pub reads_in: u64,
    pub reads_out: u64,
}

/// I/O options shared by all filter entry points.
#[derive(Debug, Clone, Copy)]
pub struct WriteOptions {
    pub compress: bool,
    pub verbose: bool,
}

/// Input/output path pair for a paired-end filter run.
pub struct PairedFiles<'a> {
    pub fwd_in: &'a Path,
    pub rev_in: &'a Path,
    pub fwd_out: &'a Path,
    pub rev_out: &'a Path,
}

// ---------------------------------------------------------------------------
// PhiX helpers
// ---------------------------------------------------------------------------

/// Read the first sequence from a FASTA file using noodles.
pub fn read_fasta_first_seq(path: &Path) -> io::Result<Vec<u8>> {
    let file = File::open(path)?;
    let mut reader = fasta::io::Reader::new(BufReader::new(file));
    let record = reader
        .records()
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "FASTA file is empty"))??;
    Ok(record.sequence().as_ref().to_vec())
}

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

fn phix_genomes(fwd: Vec<u8>) -> (Vec<u8>, Vec<u8>) {
    let rc = reverse_complement(&fwd);
    (fwd, rc)
}

fn is_phix(seq: &[u8], phix_fwd: &[u8], phix_rc: &[u8]) -> bool {
    let seqs = [seq];
    let hits_f = match_ref(&seqs, phix_fwd, 16, true);
    let hits_r = match_ref(&seqs, phix_rc, 16, true);
    // Mirror R: (hits >= minMatches) | (hits.rc >= minMatches), minMatches=2
    hits_f[0] >= 2 || hits_r[0] >= 2
}

// ---------------------------------------------------------------------------
// Low-complexity helper
// ---------------------------------------------------------------------------

/// Shannon kmer richness using 2-mers over {A, C, G, T}.
///
/// Returns exp(Shannon entropy), i.e. the effective number of distinct 2-mers,
/// matching R's `seqComplexity(seq, kmerSize=2)`.  Maximum value is 16.
fn seq_complexity_2mer(seq: &[u8]) -> f64 {
    #[inline]
    fn nt_idx(b: u8) -> Option<usize> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    let mut counts = [0u32; 16];
    let mut total = 0u32;

    for w in seq.windows(2) {
        if let (Some(i), Some(j)) = (nt_idx(w[0]), nt_idx(w[1])) {
            counts[i * 4 + j] += 1;
            total += 1;
        }
    }

    if total == 0 {
        return 1.0;
    }

    let h: f64 = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / total as f64;
            -p * p.ln()
        })
        .sum();

    h.exp()
}

// ---------------------------------------------------------------------------
// Core per-read filter
// ---------------------------------------------------------------------------

/// Apply all trim/filter operations to a single read.
///
/// Returns `None` when the read should be discarded; otherwise returns
/// `Some((trimmed_seq, trimmed_qual))`.
fn filter_read(seq: &[u8], qual: &[u8], p: &FilterParams) -> Option<(Vec<u8>, Vec<u8>)> {
    // 1. maxLen (raw, before any trimming)
    if p.max_len > 0 && seq.len() > p.max_len {
        return None;
    }

    // 2. trimLeft / trimRight boundaries
    let start = p.trim_left;
    if start >= seq.len() {
        return None;
    }
    let mut end = seq.len(); // exclusive
    if p.trim_right > 0 {
        // Need at least trim_right + 1 bases after trimLeft
        if end.saturating_sub(start) <= p.trim_right {
            return None;
        }
        end -= p.trim_right;
    }

    let seq = &seq[start..end];
    let qual = &qual[start..end];

    // 3. truncQ: truncate at first position where Phred ≤ trunc_q
    //    Phred = byte - phred_offset, so condition is byte ≤ trunc_q + phred_offset.
    let cutoff = p.trunc_q.saturating_add(p.phred_offset);
    let trunc_at = qual.iter().position(|&q| q <= cutoff).unwrap_or(qual.len());
    let seq = &seq[..trunc_at];
    let qual = &qual[..trunc_at];

    // 4. truncLen: discard if shorter than required; then trim to length.
    if p.trunc_len > 0 && seq.len() < p.trunc_len {
        return None;
    }
    let (seq, qual) = if p.trunc_len > 0 && seq.len() > p.trunc_len {
        (&seq[..p.trunc_len], &qual[..p.trunc_len])
    } else {
        (seq, qual)
    };

    // 5. minLen
    if seq.len() < p.min_len {
        return None;
    }

    // 6. maxN
    let n_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
    if n_count > p.max_n {
        return None;
    }

    // 7. minQ
    if p.min_q > 0 {
        let min_phred = qual
            .iter()
            .map(|&q| q.saturating_sub(p.phred_offset))
            .min()
            .unwrap_or(0);
        if min_phred < p.min_q {
            return None;
        }
    }

    // 8. maxEE  EE = Σ 10^(−Q/10)
    if p.max_ee.is_finite() {
        let ee: f64 = qual
            .iter()
            .map(|&q| {
                let phred = q as f64 - p.phred_offset as f64;
                10_f64.powf(-phred / 10.0)
            })
            .sum();
        if ee > p.max_ee {
            return None;
        }
    }

    Some((seq.to_vec(), qual.to_vec()))
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

fn is_gz(path: &Path) -> bool {
    path.extension().and_then(|e| e.to_str()) == Some("gz")
}

type FastqReader = fastq::io::Reader<Box<dyn BufRead>>;

fn open_reader(path: &Path) -> io::Result<FastqReader> {
    let file = File::open(path)?;
    let inner: Box<dyn BufRead> = if is_gz(path) {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    Ok(fastq::io::Reader::new(inner))
}

fn open_writer(path: &Path, compress: bool) -> io::Result<Box<dyn Write>> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let file = File::create(path)?;
    if compress {
        Ok(Box::new(GzEncoder::new(file, Compression::default())))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// Write one FASTQ record (raw bytes).
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
// Public entry points
// ---------------------------------------------------------------------------

/// Filter and trim a single-end FASTQ file.
pub fn filter_single(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    opts: WriteOptions,
) -> io::Result<SampleStats> {
    let WriteOptions { compress, verbose } = opts;
    let phix = params.phix_genome.clone().map(phix_genomes);

    let mut reader = open_reader(input)?;
    let mut writer = open_writer(output, compress)?;

    let mut reads_in: u64 = 0;
    let mut reads_out: u64 = 0;
    let mut record = fastq::Record::default();

    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e),
        }

        reads_in += 1;
        let seq = record.sequence();
        let qual = record.quality_scores();

        let Some((seq_out, qual_out)) = filter_read(seq, qual, params) else {
            continue;
        };

        // rm_phix (on trimmed sequence)
        if let Some((ref fwd, ref rc)) = phix {
            if is_phix(&seq_out, fwd, rc) {
                continue;
            }
        }

        // rm_lowcomplex (on trimmed sequence)
        if params.rm_lowcomplex > 0.0 && seq_complexity_2mer(&seq_out) < params.rm_lowcomplex {
            continue;
        }

        write_record(
            &mut *writer,
            record.name(),
            record.description(),
            &seq_out,
            &qual_out,
        )?;
        reads_out += 1;
    }

    writer.flush()?;

    if reads_out == 0 {
        // Release the writer (flushes/finalises gzip if applicable) then delete.
        drop(writer);
        let _ = std::fs::remove_file(output);
        if verbose {
            eprintln!(
                "[filter-and-trim] {} — all reads filtered out, output not written",
                input.display()
            );
        }
    } else if verbose {
        let pct = reads_out as f64 * 100.0 / reads_in as f64;
        eprintln!(
            "[filter-and-trim] {} → {} reads in, {} out ({:.1}%)",
            input.display(),
            reads_in,
            reads_out,
            pct
        );
    }

    Ok(SampleStats {
        reads_in,
        reads_out,
    })
}

/// Filter and trim a paired-end FASTQ file pair.
///
/// A pair is written to the outputs only when **both** reads pass all filters.
pub fn filter_paired(
    files: &PairedFiles<'_>,
    params_fwd: &FilterParams,
    params_rev: &FilterParams,
    opts: WriteOptions,
) -> io::Result<SampleStats> {
    let PairedFiles {
        fwd_in,
        rev_in,
        fwd_out,
        rev_out,
    } = files;
    let WriteOptions { compress, verbose } = opts;
    let phix = params_fwd
        .phix_genome
        .clone()
        .or_else(|| params_rev.phix_genome.clone())
        .map(phix_genomes);

    let mut reader_fwd = open_reader(fwd_in)?;
    let mut reader_rev = open_reader(rev_in)?;
    let mut writer_fwd = open_writer(fwd_out, compress)?;
    let mut writer_rev = open_writer(rev_out, compress)?;

    let mut reads_in: u64 = 0;
    let mut reads_out: u64 = 0;
    let mut rec_fwd = fastq::Record::default();
    let mut rec_rev = fastq::Record::default();

    loop {
        let n_fwd = reader_fwd.read_record(&mut rec_fwd)?;
        let n_rev = reader_rev.read_record(&mut rec_rev)?;

        match (n_fwd, n_rev) {
            (0, 0) => break,
            (0, _) | (_, 0) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "Paired FASTQ files have different record counts: {} vs {}",
                        fwd_in.display(),
                        rev_in.display()
                    ),
                ));
            }
            _ => {}
        }

        reads_in += 1;

        // Apply per-read filters independently.
        let filt_fwd = filter_read(rec_fwd.sequence(), rec_fwd.quality_scores(), params_fwd);
        let filt_rev = filter_read(rec_rev.sequence(), rec_rev.quality_scores(), params_rev);

        // Both must pass.
        let (Some((seq_f, qual_f)), Some((seq_r, qual_r))) = (filt_fwd, filt_rev) else {
            continue;
        };

        // rm_phix: discard pair if either read hits phiX.
        if let Some((ref phix_fwd, ref phix_rc)) = phix {
            if is_phix(&seq_f, phix_fwd, phix_rc) || is_phix(&seq_r, phix_fwd, phix_rc) {
                continue;
            }
        }

        // rm_lowcomplex: discard pair if either read is too simple.
        if params_fwd.rm_lowcomplex > 0.0 && seq_complexity_2mer(&seq_f) < params_fwd.rm_lowcomplex
        {
            continue;
        }
        if params_rev.rm_lowcomplex > 0.0 && seq_complexity_2mer(&seq_r) < params_rev.rm_lowcomplex
        {
            continue;
        }

        write_record(
            &mut *writer_fwd,
            rec_fwd.name(),
            rec_fwd.description(),
            &seq_f,
            &qual_f,
        )?;
        write_record(
            &mut *writer_rev,
            rec_rev.name(),
            rec_rev.description(),
            &seq_r,
            &qual_r,
        )?;
        reads_out += 1;
    }

    writer_fwd.flush()?;
    writer_rev.flush()?;

    if reads_out == 0 {
        drop(writer_fwd);
        drop(writer_rev);
        let _ = std::fs::remove_file(fwd_out);
        let _ = std::fs::remove_file(rev_out);
        if verbose {
            eprintln!(
                "[filter-and-trim] {}/{} — all pairs filtered out, outputs not written",
                fwd_in.display(),
                rev_in.display()
            );
        }
    } else if verbose {
        let pct = reads_out as f64 * 100.0 / reads_in as f64;
        eprintln!(
            "[filter-and-trim] {}/{} → {} pairs in, {} out ({:.1}%)",
            fwd_in.display(),
            rev_in.display(),
            reads_in,
            reads_out,
            pct
        );
    }

    Ok(SampleStats {
        reads_in,
        reads_out,
    })
}
