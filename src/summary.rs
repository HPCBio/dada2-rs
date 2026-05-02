use std::io::{self, BufReader};

use noodles::fastq;
use rayon::prelude::*;

/// Highest quality score we track per position. 0..=MAX_QUAL inclusive.
/// Covers Illumina (≤41), Sanger (≤40), and PacBio HiFi (≤93) ranges.
const MAX_QUAL: usize = 93;

pub struct QualitySummary {
    pub total_reads: u64,
    sums: Vec<f64>,
    counts: Vec<u64>,
    /// Per-cycle quality distribution: `hist[pos][q]` is the count of reads
    /// with quality `q` at zero-based cycle `pos`.
    hist: Vec<[u64; MAX_QUAL + 1]>,
}

impl QualitySummary {
    pub fn new() -> Self {
        Self {
            total_reads: 0,
            sums: Vec::new(),
            counts: Vec::new(),
            hist: Vec::new(),
        }
    }

    fn add_record(&mut self, quality: &[u8], phred_offset: u8) {
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
}

/// Records per thread per chunk — total chunk size scales with thread count.
const RECORDS_PER_THREAD: usize = 10_000;

pub fn process<R: io::Read>(
    reader: R,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
) -> io::Result<QualitySummary> {
    let chunk_size = RECORDS_PER_THREAD * pool.current_num_threads();
    let buf = BufReader::new(reader);
    let mut fastq_reader = fastq::io::Reader::new(buf);
    let mut overall = QualitySummary::new();

    loop {
        // Read a chunk sequentially — the reader is a stream and cannot be shared.
        let mut chunk: Vec<Vec<u8>> = Vec::with_capacity(chunk_size);
        let mut record = fastq::Record::default();
        let mut error: Option<io::Error> = None;

        for _ in 0..chunk_size {
            match fastq_reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => chunk.push(record.quality_scores().to_vec()),
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
                .fold(QualitySummary::new, |mut acc, quality| {
                    acc.add_record(quality, phred_offset);
                    acc
                })
                .reduce(QualitySummary::new, QualitySummary::merge)
        });

        overall = overall.merge(partial);

        if done {
            break;
        }
    }

    Ok(overall)
}
