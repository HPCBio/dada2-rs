use std::collections::HashMap;
use std::io::{self, BufReader};

use noodles::fastq;
use rayon::prelude::*;

/// Mirrors the R dada2 `derep` class.
///
/// - `uniques`: unique sequences in first-seen order, each paired with its read count.
/// - `quals`:   mean Phred quality score per position for each unique sequence;
///              `quals[i][j]` is the mean score at position `j` for unique `i`.
/// - `map`:     for each input read (in order), the index into `uniques` of the
///              unique sequence it maps to.
pub struct Derep {
    pub uniques: Vec<(Vec<u8>, u64)>,
    pub quals: Vec<Vec<f64>>,
    pub map: Vec<usize>,
}

/// Per-thread accumulator.  Can be merged in order to preserve read ordering.
struct PartialDerep {
    seq_order: Vec<Vec<u8>>,
    seq_to_idx: HashMap<Vec<u8>, usize>,
    counts: Vec<u64>,
    qual_sums: Vec<Vec<f64>>,
    qual_cnts: Vec<Vec<u64>>,
    map: Vec<usize>,
}

impl PartialDerep {
    fn new() -> Self {
        Self {
            seq_order: Vec::new(),
            seq_to_idx: HashMap::new(),
            counts: Vec::new(),
            qual_sums: Vec::new(),
            qual_cnts: Vec::new(),
            map: Vec::new(),
        }
    }

    fn add_record(&mut self, seq: Vec<u8>, qual: &[u8], phred_offset: u8) {
        let idx = match self.seq_to_idx.get(&seq) {
            Some(&i) => i,
            None => {
                let i = self.seq_order.len();
                self.seq_to_idx.insert(seq.clone(), i);
                self.seq_order.push(seq.clone());
                self.counts.push(0);
                self.qual_sums.push(Vec::new());
                self.qual_cnts.push(Vec::new());
                i
            }
        };

        self.counts[idx] += 1;
        self.map.push(idx);

        let sums = &mut self.qual_sums[idx];
        let cnts = &mut self.qual_cnts[idx];
        if qual.len() > sums.len() {
            sums.resize(qual.len(), 0.0);
            cnts.resize(qual.len(), 0);
        }
        for (j, &q) in qual.iter().enumerate() {
            sums[j] += (q as i16 - phred_offset as i16) as f64;
            cnts[j] += 1;
        }
    }

    /// Merge `other` (which covers reads that come *after* `self`) into `self`.
    /// `other`'s local indices are remapped into `self`'s index space so that
    /// the combined `map` remains in correct read order.
    fn merge(mut self, other: PartialDerep) -> PartialDerep {
        let mut remap = vec![0usize; other.seq_order.len()];

        for (i, seq) in other.seq_order.iter().enumerate() {
            let idx = if let Some(&j) = self.seq_to_idx.get(seq) {
                // Sequence already known — accumulate quality and count.
                let sums = &mut self.qual_sums[j];
                let cnts = &mut self.qual_cnts[j];
                let osums = &other.qual_sums[i];
                let ocnts = &other.qual_cnts[i];
                if osums.len() > sums.len() {
                    sums.resize(osums.len(), 0.0);
                    cnts.resize(ocnts.len(), 0);
                }
                for k in 0..osums.len() {
                    sums[k] += osums[k];
                    cnts[k] += ocnts[k];
                }
                self.counts[j] += other.counts[i];
                j
            } else {
                let j = self.seq_order.len();
                self.seq_to_idx.insert(seq.clone(), j);
                self.seq_order.push(seq.clone());
                self.counts.push(other.counts[i]);
                self.qual_sums.push(other.qual_sums[i].clone());
                self.qual_cnts.push(other.qual_cnts[i].clone());
                j
            };
            remap[i] = idx;
        }

        for &local_idx in &other.map {
            self.map.push(remap[local_idx]);
        }

        self
    }

    fn into_derep(self) -> Derep {
        let quals = self
            .qual_sums
            .iter()
            .zip(self.qual_cnts.iter())
            .map(|(sums, cnts)| {
                sums.iter()
                    .zip(cnts.iter())
                    .map(|(&s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
                    .collect()
            })
            .collect();

        let uniques = self
            .seq_order
            .into_iter()
            .zip(self.counts)
            .collect();

        Derep { uniques, quals, map: self.map }
    }
}

/// Records assigned to each thread per chunk — total chunk size scales with thread count.
const RECORDS_PER_THREAD: usize = 10_000;

pub fn dereplicate<R: io::Read>(
    reader: R,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    verbose: bool,
) -> io::Result<Derep> {
    let chunk_size = RECORDS_PER_THREAD * pool.current_num_threads();
    let buf = BufReader::new(reader);
    let mut fastq_reader = fastq::io::Reader::new(buf);
    let mut overall = PartialDerep::new();

    loop {
        // Read a chunk sequentially — the reader is a stream and cannot be shared.
        let mut chunk: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(chunk_size);
        let mut record = fastq::Record::default();
        let mut error: Option<io::Error> = None;

        for _ in 0..chunk_size {
            match fastq_reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => chunk.push((
                    record.sequence().to_vec(),
                    record.quality_scores().to_vec(),
                )),
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

        // Dereplicate the chunk in parallel, then merge order-preserving into overall.
        let partial = pool.install(|| {
            chunk
                .par_iter()
                .fold(PartialDerep::new, |mut acc, (seq, qual)| {
                    acc.add_record(seq.clone(), qual, phred_offset);
                    acc
                })
                .reduce(PartialDerep::new, PartialDerep::merge)
        });

        overall = overall.merge(partial);

        if done {
            break;
        }
    }

    let derep = overall.into_derep();
    if verbose {
        eprintln!(
            "[derep] {} raw sequences -> {} unique sequences",
            derep.map.len(),
            derep.uniques.len()
        );
    }
    Ok(derep)
}
