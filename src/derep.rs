use std::collections::HashMap;
use std::io::{self, BufReader};

use noodles::fastq;
use rayon::prelude::*;

/// Per-position Phred quality sums are stored as `u32` (deferred division,
/// issue #23: keep the integer sum and divide by abundance on demand).
///
/// **Capacity cap.** A sum overflows `u32` only when a *single* unique sequence
/// accumulates more than `u32::MAX / Q_max` reads. With the Phred ceiling (~93)
/// that is roughly **46 million reads of one byte-identical sequence** — far
/// beyond any realistic dataset (per-sample totals rarely exceed ~1M paired-end;
/// 46M reads of one *exact* sequence would imply ~1B+ reads pooled into highly
/// repetitive data). If this is ever hit, the fix is to widen the qual-sum
/// storage (`RawInput.quals` / `Derep.quals`) to `u64`. Until then we fail
/// loudly via [`checked_qual_sum`] rather than silently saturate.
pub const QUAL_SUM_MAX: u32 = u32::MAX;

/// Convert an accumulated per-position Phred sum (`f64`, but integer-valued) to
/// the stored `u32`, panicking with a clear, actionable message on overflow
/// instead of silently saturating — `f64 as u32` clamps to `u32::MAX`, which
/// would corrupt the recovered mean quality without any signal. See
/// [`QUAL_SUM_MAX`] for when this can happen (essentially never).
#[inline]
pub fn checked_qual_sum(sum: f64) -> u32 {
    assert!(
        sum <= QUAL_SUM_MAX as f64,
        "per-position Phred quality sum ({sum:.0}) overflows u32: a single unique \
         sequence has more than ~46M reads. dada2-rs stores quality as integer \
         sums per position (issue #23); to support data this deep, widen \
         RawInput.quals / Derep.quals to u64. Please open an issue reporting your \
         per-sample read counts."
    );
    // Negatives (malformed quals below the Phred offset) pass the cap check and,
    // as before, saturate to 0 here — that pathology is handled upstream, not by
    // this overflow guard.
    sum.round() as u32
}

/// Mirrors the R dada2 `derep` class.
///
/// - `uniques`: unique sequences sorted by read count descending, ties broken
///   by sequence (lexical, ascending). Matches R `derepFastq` ordering, where
///   `qtables2` builds uniques in lexical order and the final stable
///   `order(decreasing=TRUE)` leaves equal-count uniques in that lexical order.
/// - `quals`:   per-position *integer* Phred quality SUM over the reads in each
///              unique (deferred division, issue #23); `quals[i][j]` is the
///              summed score at position `j` for unique `i`. Consumers recover
///              the mean on demand as `sum / count` (the unique's abundance) —
///              the per-position count equals the unique's read count because
///              dereplication groups by exact full-sequence identity, so every
///              read covers every position. Storing the sum is bit-identical to
///              the old f64 mean while using 4 bytes/position instead of 8.
/// - `map`:     for each input read (in order), the index into `uniques` of the
///              unique sequence it maps to. **Built only when requested** (see
///              `build_map` in [`dereplicate`]); empty otherwise. It scales with
///              the number of *reads* (8 bytes each), so it dominates resident
///              memory for deep samples and is opt-in (issue #36). The read count
///              is always available via `n_reads` regardless.
/// - `n_reads`: total input reads (= Σ unique counts). Always populated, so read
///              counts don't require the per-read `map`.
pub struct Derep {
    pub uniques: Vec<(Vec<u8>, u64)>,
    pub quals: Vec<Vec<u32>>,
    pub map: Vec<usize>,
    pub n_reads: u64,
}

/// Per-thread accumulator.  Can be merged in order to preserve read ordering.
struct PartialDerep {
    seq_order: Vec<Vec<u8>>,
    seq_to_idx: HashMap<Vec<u8>, usize>,
    counts: Vec<u64>,
    qual_sums: Vec<Vec<f64>>,
    qual_cnts: Vec<Vec<u64>>,
    map: Vec<usize>,
    /// When false, the per-read `map` is not accumulated (issue #36) — uniques,
    /// quals, and counts are unaffected.
    build_map: bool,
}

impl PartialDerep {
    fn new(build_map: bool) -> Self {
        Self {
            seq_order: Vec::new(),
            seq_to_idx: HashMap::new(),
            counts: Vec::new(),
            qual_sums: Vec::new(),
            qual_cnts: Vec::new(),
            map: Vec::new(),
            build_map,
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
        if self.build_map {
            self.map.push(idx);
        }

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

        if self.build_map {
            for &local_idx in &other.map {
                self.map.push(remap[local_idx]);
            }
        }

        self
    }

    fn into_derep(self) -> Derep {
        // Keep the integer per-position sum (deferred division, issue #23):
        // consumers divide by the unique's count on demand. `qual_cnts` is
        // redundant under exact-match dereplication (every read covers every
        // position), so it equals the unique's count and we don't store it.
        let quals: Vec<Vec<u32>> = self
            .qual_sums
            .iter()
            .map(|sums| sums.iter().map(|&s| checked_qual_sum(s)).collect())
            .collect();

        let n = self.seq_order.len();

        // Sort uniques by abundance descending, ties broken by sequence
        // (lexical, ascending) to match R `derepFastq`: its `qtables2` builds
        // uniques in lexical order, then a stable `order(decreasing=TRUE)`
        // leaves equal-count uniques in that lexical order. The downstream
        // DADA2 algorithm assumes this ordering when traversing raws — the
        // most-abundant raw lands at index 0 so the cluster-0 center is both
        // the most abundant and the lowest-indexed eligible raw. Issue #4
        // traced part of the over-budding to our previous first-seen ordering
        // disagreeing with R; the lexical tie-break closes the remaining gap
        // for equal-abundance uniques. (Total order on distinct sequences, so
        // sort stability is irrelevant.)
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            self.counts[b]
                .cmp(&self.counts[a])
                .then_with(|| self.seq_order[a].cmp(&self.seq_order[b]))
        });

        // Apply the permutation to uniques and quals; remap each `map`
        // entry from old → new index.
        let mut new_seq_order: Vec<Vec<u8>> = Vec::with_capacity(n);
        let mut new_counts: Vec<u64> = Vec::with_capacity(n);
        let mut new_quals: Vec<Vec<u32>> = Vec::with_capacity(n);
        let mut old_to_new: Vec<usize> = vec![0; n];
        let seq_order_owned = self.seq_order;
        let mut seq_iter: Vec<Option<Vec<u8>>> = seq_order_owned.into_iter().map(Some).collect();
        let mut quals_iter: Vec<Option<Vec<u32>>> = quals.into_iter().map(Some).collect();
        for (new_idx, &old_idx) in order.iter().enumerate() {
            old_to_new[old_idx] = new_idx;
            new_seq_order.push(seq_iter[old_idx].take().unwrap());
            new_counts.push(self.counts[old_idx]);
            new_quals.push(quals_iter[old_idx].take().unwrap());
        }
        // Read count is the sum of per-unique counts — always exact, and
        // independent of whether the per-read map was built (issue #36).
        let n_reads: u64 = new_counts.iter().sum();
        let map: Vec<usize> = if self.build_map {
            self.map.into_iter().map(|i| old_to_new[i]).collect()
        } else {
            Vec::new()
        };

        let uniques = new_seq_order.into_iter().zip(new_counts).collect();

        Derep {
            uniques,
            quals: new_quals,
            map,
            n_reads,
        }
    }
}

/// Records assigned to each thread per chunk — total chunk size scales with thread count.
const RECORDS_PER_THREAD: usize = 10_000;

/// Dereplicate a FASTQ stream into unique sequences with per-position quality
/// sums and read counts.
///
/// `build_map`: when `true`, also record the per-read → unique-index `map`
/// (needed only for read-level traceback, e.g. the `derep --show-map` output;
/// `merge-pairs` builds its own internally). When `false`, the `map` is skipped
/// entirely — it scales with read count (8 bytes/read) and dominates resident
/// memory for deep pooled runs (issue #36). Uniques/quals/counts are identical
/// either way; the read count is always available via `Derep::n_reads`.
pub fn dereplicate<R: io::Read>(
    reader: R,
    phred_offset: u8,
    build_map: bool,
    pool: &rayon::ThreadPool,
    verbose: bool,
) -> io::Result<Derep> {
    let chunk_size = RECORDS_PER_THREAD * pool.current_num_threads();
    let buf = BufReader::new(reader);
    let mut fastq_reader = fastq::io::Reader::new(buf);
    let mut overall = PartialDerep::new(build_map);

    loop {
        // Read a chunk sequentially — the reader is a stream and cannot be shared.
        let mut chunk: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(chunk_size);
        let mut record = fastq::Record::default();
        let mut error: Option<io::Error> = None;

        for _ in 0..chunk_size {
            match fastq_reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => chunk.push((record.sequence().to_vec(), record.quality_scores().to_vec())),
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

        // Dereplicate the chunk in parallel, then merge order-preserving into overall.
        let partial = pool.install(|| {
            chunk
                .par_iter()
                .fold(
                    || PartialDerep::new(build_map),
                    |mut acc, (seq, qual)| {
                        acc.add_record(seq.clone(), qual, phred_offset);
                        acc
                    },
                )
                .reduce(|| PartialDerep::new(build_map), PartialDerep::merge)
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
            derep.n_reads,
            derep.uniques.len()
        );
    }
    Ok(derep)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_qual_sum_passes_through_in_range() {
        assert_eq!(checked_qual_sum(0.0), 0);
        assert_eq!(checked_qual_sum(37.4), 37); // rounds like the f64-mean path did
        assert_eq!(checked_qual_sum(QUAL_SUM_MAX as f64), QUAL_SUM_MAX);
        // Negative (malformed quals below the Phred offset) saturates to 0 as before.
        assert_eq!(checked_qual_sum(-5.0), 0);
    }

    #[test]
    #[should_panic(expected = "overflows u32")]
    fn checked_qual_sum_panics_on_overflow() {
        // ~46M reads of one unique at Q40 would land here; we fail loudly rather
        // than silently saturate to u32::MAX.
        checked_qual_sum(QUAL_SUM_MAX as f64 + 1.0);
    }
}
