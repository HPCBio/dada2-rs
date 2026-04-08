use std::io::{self, BufReader};

use noodles::fastq;
use rayon::prelude::*;

pub struct QualitySummary {
    pub total_reads: u64,
    sums: Vec<f64>,
    counts: Vec<u64>,
}

impl QualitySummary {
    pub fn new() -> Self {
        Self {
            total_reads: 0,
            sums: Vec::new(),
            counts: Vec::new(),
        }
    }

    fn add_record(&mut self, quality: &[u8], phred_offset: u8) {
        if quality.len() > self.sums.len() {
            self.sums.resize(quality.len(), 0.0);
            self.counts.resize(quality.len(), 0);
        }
        for (i, &q) in quality.iter().enumerate() {
            self.sums[i] += ((q as i16) - (phred_offset as i16)) as f64;
            self.counts[i] += 1;
        }
        self.total_reads += 1;
    }

    fn merge(mut self, other: QualitySummary) -> QualitySummary {
        let len = self.sums.len().max(other.sums.len());
        self.sums.resize(len, 0.0);
        self.counts.resize(len, 0);
        for (i, (s, c)) in other.sums.iter().zip(other.counts.iter()).enumerate() {
            self.sums[i] += s;
            self.counts[i] += c;
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
                Err(e) => { error = Some(e); break; }
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
