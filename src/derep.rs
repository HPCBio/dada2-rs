use std::collections::HashMap;
use std::io::{self, BufReader};

use noodles::fastq;
use rayon::prelude::*;

/// Mirrors the R dada2 `derep` class.
///
/// - `uniques`: unique sequences sorted by read count descending (stable;
///   ties preserve first-seen order). Matches R `derepFastq` ordering.
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
        let quals: Vec<Vec<f64>> = self
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

        let n = self.seq_order.len();

        // Sort uniques by abundance descending (stable: first-seen order
        // within count ties).  Matches R `derepFastq` ordering, which the
        // downstream DADA2 algorithm assumes when traversing raws — the
        // most-abundant raw lands at index 0 so the cluster-0 center is
        // both the most abundant and the lowest-indexed eligible raw.
        // Issue #4 traced part of the over-budding to our previous
        // first-seen ordering disagreeing with R.
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| self.counts[b].cmp(&self.counts[a]));

        // Apply the permutation to uniques and quals; remap each `map`
        // entry from old → new index.
        let mut new_seq_order: Vec<Vec<u8>> = Vec::with_capacity(n);
        let mut new_counts: Vec<u64> = Vec::with_capacity(n);
        let mut new_quals: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut old_to_new: Vec<usize> = vec![0; n];
        let seq_order_owned = self.seq_order;
        let mut seq_iter: Vec<Option<Vec<u8>>> = seq_order_owned.into_iter().map(Some).collect();
        let mut quals_iter: Vec<Option<Vec<f64>>> = quals.into_iter().map(Some).collect();
        for (new_idx, &old_idx) in order.iter().enumerate() {
            old_to_new[old_idx] = new_idx;
            new_seq_order.push(seq_iter[old_idx].take().unwrap());
            new_counts.push(self.counts[old_idx]);
            new_quals.push(quals_iter[old_idx].take().unwrap());
        }
        let map: Vec<usize> = self.map.into_iter().map(|i| old_to_new[i]).collect();

        let uniques = new_seq_order.into_iter().zip(new_counts).collect();

        Derep { uniques, quals: new_quals, map }
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
