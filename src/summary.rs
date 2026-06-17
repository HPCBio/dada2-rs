use std::io::{self, BufReader};

use noodles::fastq;
use rayon::prelude::*;

/// Highest quality score we track per position. 0..=MAX_QUAL inclusive.
/// Covers Illumina (≤41), Sanger (≤40), and PacBio HiFi (≤93) ranges.
const MAX_QUAL: usize = 93;

/// Configuration for the optional per-read sequence-complexity histogram, a port
/// of DADA2's `seqComplexity`/`plotComplexity` (Benjamin Callahan; see
/// dada2/R/filter.R and plot-methods.R). Complexity is the effective number of
/// k-mers in a read — `exp(Shannon entropy)` over its overlapping k-mer counts,
/// ranging `[1, 4^kmer_size]`. Values are binned over `[0, 4^kmer_size]`.
#[derive(Clone, Copy)]
pub struct ComplexityConfig {
    pub kmer_size: u8,
    pub bins: usize,
}

impl ComplexityConfig {
    /// Maximum possible complexity, `4^kmer_size` (all k-mers equally frequent).
    fn max_complexity(&self) -> f64 {
        4f64.powi(self.kmer_size as i32)
    }
}

/// Effective number of k-mers in `seq` = `exp(-Σ p·ln p)` over the read's
/// overlapping k-mer counts. K-mers containing non-ACGT bases are skipped
/// (matching Biostrings `oligonucleotideFrequency`, which ignores ambiguous
/// k-mers). A read with no valid k-mer yields 1.0 (R's `sindex` of an all-zero
/// frequency vector is `exp(0) = 1`).
fn seq_complexity(seq: &[u8], k: usize) -> f64 {
    if seq.len() < k || k == 0 {
        return 1.0;
    }
    let mut counts = vec![0u32; 1usize << (2 * k)];
    let mut total = 0u64;
    'window: for w in seq.windows(k) {
        let mut idx = 0usize;
        for &b in w {
            let code = match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => continue 'window,
            };
            idx = (idx << 2) | code;
        }
        counts[idx] += 1;
        total += 1;
    }
    if total == 0 {
        return 1.0;
    }
    let t = total as f64;
    let mut h = 0.0;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / t;
            h -= p * p.ln();
        }
    }
    h.exp()
}

pub struct QualitySummary {
    pub total_reads: u64,
    sums: Vec<f64>,
    counts: Vec<u64>,
    /// Per-cycle quality distribution: `hist[pos][q]` is the count of reads
    /// with quality `q` at zero-based cycle `pos`.
    hist: Vec<[u64; MAX_QUAL + 1]>,
    /// Optional per-read complexity histogram. `None` unless requested.
    complexity: Option<ComplexityConfig>,
    /// `complexity_hist[bin]` = number of reads whose effective k-mer count falls
    /// in `bin` (equal-width bins over `[0, 4^kmer_size]`). Empty when disabled.
    complexity_hist: Vec<u64>,
}

impl QualitySummary {
    pub fn with_complexity(complexity: Option<ComplexityConfig>) -> Self {
        let complexity_hist = match complexity {
            Some(cfg) => vec![0; cfg.bins],
            None => Vec::new(),
        };
        Self {
            total_reads: 0,
            sums: Vec::new(),
            counts: Vec::new(),
            hist: Vec::new(),
            complexity,
            complexity_hist,
        }
    }

    fn add_record(&mut self, sequence: &[u8], quality: &[u8], phred_offset: u8) {
        if quality.len() > self.sums.len() {
            self.sums.resize(quality.len(), 0.0);
            self.counts.resize(quality.len(), 0);
            self.hist.resize(quality.len(), [0; MAX_QUAL + 1]);
        }
        for (i, &q) in quality.iter().enumerate() {
            let q_phred = (q as i16) - (phred_offset as i16);
            self.sums[i] += q_phred as f64;
            self.counts[i] += 1;
            let idx = q_phred.clamp(0, MAX_QUAL as i16) as usize;
            self.hist[i][idx] += 1;
        }
        if let Some(cfg) = self.complexity {
            let si = seq_complexity(sequence, cfg.kmer_size as usize);
            let frac = si / cfg.max_complexity();
            let bin = ((frac * cfg.bins as f64) as usize).min(cfg.bins - 1);
            self.complexity_hist[bin] += 1;
        }
        self.total_reads += 1;
    }

    fn merge(mut self, other: QualitySummary) -> QualitySummary {
        let len = self.sums.len().max(other.sums.len());
        self.sums.resize(len, 0.0);
        self.counts.resize(len, 0);
        self.hist.resize(len, [0; MAX_QUAL + 1]);
        for (i, (s, c)) in other.sums.iter().zip(other.counts.iter()).enumerate() {
            self.sums[i] += s;
            self.counts[i] += c;
        }
        for (i, row) in other.hist.iter().enumerate() {
            for (q, &n) in row.iter().enumerate() {
                self.hist[i][q] += n;
            }
        }
        if self.complexity_hist.len() < other.complexity_hist.len() {
            self.complexity_hist.resize(other.complexity_hist.len(), 0);
        }
        for (i, &n) in other.complexity_hist.iter().enumerate() {
            self.complexity_hist[i] += n;
        }
        if self.complexity.is_none() {
            self.complexity = other.complexity;
        }
        self.total_reads += other.total_reads;
        self
    }

    pub fn mean_quality_per_position(&self) -> Vec<f64> {
        self.sums
            .iter()
            .zip(self.counts.iter())
            .map(|(sum, &count)| if count > 0 { sum / count as f64 } else { 0.0 })
            .collect()
    }

    /// Per-position read coverage (reads with a base at each cycle).
    pub fn reads_per_position(&self) -> &[u64] {
        &self.counts
    }

    /// Per-position quality histogram trimmed to the highest quality observed
    /// across any position. Returns `(max_quality, hist[pos][0..=max_quality])`.
    pub fn quality_histogram(&self) -> (usize, Vec<Vec<u64>>) {
        let mut max_q = 0usize;
        for row in &self.hist {
            for (q, &n) in row.iter().enumerate() {
                if n > 0 && q > max_q {
                    max_q = q;
                }
            }
        }
        let trimmed = self.hist.iter().map(|row| row[..=max_q].to_vec()).collect();
        (max_q, trimmed)
    }

    /// The per-read complexity histogram, if it was requested.
    /// Returns `(kmer_size, bins, counts)` where `counts[bin]` covers the
    /// effective-k-mer-count range `[bin, bin+1) · 4^kmer_size / bins`.
    pub fn complexity_histogram(&self) -> Option<(u8, usize, &[u64])> {
        self.complexity
            .map(|cfg| (cfg.kmer_size, cfg.bins, self.complexity_hist.as_slice()))
    }
}

/// Records per thread per chunk — total chunk size scales with thread count.
const RECORDS_PER_THREAD: usize = 10_000;

pub fn process<R: io::Read>(
    reader: R,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    complexity: Option<ComplexityConfig>,
) -> io::Result<QualitySummary> {
    let chunk_size = RECORDS_PER_THREAD * pool.current_num_threads();
    let buf = BufReader::new(reader);
    let mut fastq_reader = fastq::io::Reader::new(buf);
    let mut overall = QualitySummary::with_complexity(complexity);

    // Only retain the sequence when complexity is requested — otherwise it's
    // dead weight in the chunk buffer.
    let want_seq = complexity.is_some();

    loop {
        // Read a chunk sequentially — the reader is a stream and cannot be shared.
        let mut chunk: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(chunk_size);
        let mut record = fastq::Record::default();
        let mut error: Option<io::Error> = None;

        for _ in 0..chunk_size {
            match fastq_reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    let seq = if want_seq {
                        record.sequence().to_vec()
                    } else {
                        Vec::new()
                    };
                    chunk.push((seq, record.quality_scores().to_vec()));
                }
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }

        if let Some(e) = error {
            return Err(e);
        }
        if chunk.is_empty() {
            break;
        }

        let done = chunk.len() < chunk_size;

        // Process the chunk in parallel within the configured thread pool.
        let partial = pool.install(|| {
            chunk
                .par_iter()
                .fold(
                    || QualitySummary::with_complexity(complexity),
                    |mut acc, (sequence, quality)| {
                        acc.add_record(sequence, quality, phred_offset);
                        acc
                    },
                )
                .reduce(
                    || QualitySummary::with_complexity(complexity),
                    QualitySummary::merge,
                )
        });

        overall = overall.merge(partial);

        if done {
            break;
        }
    }

    Ok(overall)
}
