//! Dereplication: one sample (R `derepFastq`) and a pool of samples (R
//! `combineDereps2`, `R/multiSample.R`). Ported from DADA2 by Benjamin
//! Callahan.

use std::collections::HashMap;
use std::io::{self, BufReader};
use std::sync::OnceLock;

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
///              unique sequence it maps to.
pub struct Derep {
    pub uniques: Vec<(Vec<u8>, u64)>,
    pub quals: Vec<Vec<u32>>,
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
        let map: Vec<usize> = self.map.into_iter().map(|i| old_to_new[i]).collect();

        let uniques = new_seq_order.into_iter().zip(new_counts).collect();

        Derep {
            uniques,
            quals: new_quals,
            map,
        }
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

/// How [`DerepPool::finish`] breaks abundance ties (#260). Selected by
/// `DADA2RS_POOL_TIEBREAK`.
///
/// - `first-seen` (default, or unset): by first appearance across samples in
///   input order, as R's `combineDereps2` does with its stable `order()`.
///   Pooled results depend on the order samples are given, as R's do. On the
///   cluster A/B this matched R's pooled ASV set exactly on MiSeq and ITS2
///   (F and R) and cut PacBio's residual from 17 ASVs to 10.
/// - `lexical`: by sequence, as `derepFastq` orders a single sample. The
///   pre-#260 behaviour, kept as a result-changing arm for comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoolTiebreak {
    Lexical,
    FirstSeen,
}

impl PoolTiebreak {
    /// Parse a `DADA2RS_POOL_TIEBREAK` value.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "" | "first-seen" => Ok(Self::FirstSeen),
            "lexical" => Ok(Self::Lexical),
            v => Err(format!(
                "DADA2RS_POOL_TIEBREAK={v:?} is not recognised; expected first-seen or lexical"
            )),
        }
    }

    /// The value as it would be written in the environment.
    pub fn label(self) -> &'static str {
        match self {
            Self::Lexical => "lexical",
            Self::FirstSeen => "first-seen",
        }
    }
}

/// The resolved tie-break, read once per process. An unparseable value is
/// fatal: a mistyped arm silently running the default would compare the
/// default with itself (see `gates`).
pub fn pool_tiebreak() -> PoolTiebreak {
    static VALUE: OnceLock<PoolTiebreak> = OnceLock::new();
    *VALUE.get_or_init(|| match std::env::var("DADA2RS_POOL_TIEBREAK") {
        Ok(v) => PoolTiebreak::parse(&v).unwrap_or_else(|e| panic!("{e}")),
        Err(_) => PoolTiebreak::FirstSeen,
    })
}

/// Folds per-sample dereplications into one pooled unique table for
/// `dada-pooled` (R `combineDereps2`). Samples are added one at a time and can
/// be dropped after [`DerepPool::add`], so only the accumulator and the sample
/// being added are resident (#41).
#[derive(Default)]
pub struct DerepPool {
    seq_to_merged: HashMap<Vec<u8>, usize>,
    /// Merged uniques in first-seen order until [`DerepPool::finish`].
    seqs: Vec<Vec<u8>>,
    /// Per-position Phred sums, added as integers: per-sample quals are already
    /// `u32` sums (#23), so `u32` halves the largest merge intermediate (#39).
    qual_sums: Vec<Vec<u32>>,
    abundance: Vec<u32>,
    local_to_merged: Vec<Vec<usize>>,
    sample_counts: Vec<Vec<u32>>,
}

/// The pooled unique table, most abundant first.
pub struct PooledDerep {
    /// Merged unique sequences.
    pub seqs: Vec<Vec<u8>>,
    /// Per-position Phred sums for each merged unique, across all samples.
    pub qual_sums: Vec<Vec<u32>>,
    /// Total reads of each merged unique, across all samples.
    pub abundance: Vec<u32>,
    /// Per sample, in the order added: each local unique's merged index.
    pub local_to_merged: Vec<Vec<usize>>,
    /// Per sample: each local unique's read count.
    pub sample_counts: Vec<Vec<u32>>,
}

impl DerepPool {
    pub fn new() -> Self {
        Self::default()
    }

    /// Fold one sample into the pool. Samples must be added in input order:
    /// `PooledDerep::local_to_merged` is indexed by it.
    pub fn add(&mut self, derep: &Derep) {
        let mut local_map: Vec<usize> = Vec::with_capacity(derep.uniques.len());
        let mut counts: Vec<u32> = Vec::with_capacity(derep.uniques.len());
        for ((seq, count), qual) in derep.uniques.iter().zip(derep.quals.iter()) {
            let count_u32 = *count as u32;
            let mu = match self.seq_to_merged.get(seq) {
                Some(&j) => {
                    self.abundance[j] += count_u32;
                    // `qual` is already this unique's per-position Phred
                    // sum; accumulate sums across samples.
                    for (p, &q) in qual.iter().enumerate() {
                        self.qual_sums[j][p] = self.qual_sums[j][p].checked_add(q).expect(
                            "merged per-position Phred sum overflows u32 \
                                     (pooled depth extreme); widen DerepPool::qual_sums to u64",
                        );
                    }
                    j
                }
                None => {
                    let j = self.seqs.len();
                    self.seq_to_merged.insert(seq.clone(), j);
                    self.seqs.push(seq.clone());
                    self.qual_sums.push(qual.clone());
                    self.abundance.push(count_u32);
                    j
                }
            };
            local_map.push(mu);
            counts.push(count_u32);
        }
        self.local_to_merged.push(local_map);
        self.sample_counts.push(counts);
    }

    /// Merged uniques so far.
    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    /// Reads pooled so far.
    pub fn total_reads(&self) -> u32 {
        self.abundance.iter().sum()
    }

    /// Order the pool by descending abundance and hand it over (#219).
    ///
    /// The DADA2 loop assumes the most abundant unique is at index 0: `b_bud`'s
    /// scan is `for r in 1..` ("r=0 is the center"), which holds only because
    /// cluster 0 starts with every raw in input order and `assign_center` picks
    /// the most abundant. Neither R nor we move the centre into place. A pool
    /// left in first-seen order made index 0 permanently unbuddable: a
    /// 66,937-read organism on the pinned 95-sample PacBio run, which R calls
    /// as its own ASV. Order also decides saturated births (pA = 0, 880 of
    /// 2,490 pooled MiSeq births), where the tie-break is reads, then position.
    ///
    /// Ties are broken by [`pool_tiebreak`]: in first-seen order by default,
    /// as `combineDereps2`'s stable `order()` leaves them (#260).
    pub fn finish(self) -> PooledDerep {
        self.finish_with(pool_tiebreak())
    }

    /// [`DerepPool::finish`] with an explicit tie-break.
    pub fn finish_with(self, tiebreak: PoolTiebreak) -> PooledDerep {
        let DerepPool {
            seq_to_merged,
            seqs,
            qual_sums,
            abundance,
            mut local_to_merged,
            sample_counts,
        } = self;
        drop(seq_to_merged);
        let mut perm: Vec<usize> = (0..seqs.len()).collect();
        // `sort_by` is stable, so `FirstSeen` keeps the pool's first-seen order
        // among equal abundances.
        perm.sort_by(|&a, &b| {
            let by_abundance = abundance[b].cmp(&abundance[a]);
            match tiebreak {
                PoolTiebreak::Lexical => by_abundance.then_with(|| seqs[a].cmp(&seqs[b])),
                PoolTiebreak::FirstSeen => by_abundance,
            }
        });
        let mut old_to_new = vec![0usize; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx;
        }
        let mut seq_slots: Vec<Option<Vec<u8>>> = seqs.into_iter().map(Some).collect();
        let mut qual_slots: Vec<Option<Vec<u32>>> = qual_sums.into_iter().map(Some).collect();
        let seqs = perm
            .iter()
            .map(|&i| seq_slots[i].take().expect("permutation is a bijection"))
            .collect();
        let qual_sums = perm
            .iter()
            .map(|&i| qual_slots[i].take().expect("permutation is a bijection"))
            .collect();
        let abundance = perm.iter().map(|&i| abundance[i]).collect();
        // Every per-sample map indexes the pool, so it has to follow the
        // permutation or the split-back attributes reads to the wrong sequence.
        for local in local_to_merged.iter_mut() {
            for mu in local.iter_mut() {
                *mu = old_to_new[*mu];
            }
        }
        PooledDerep {
            seqs,
            qual_sums,
            abundance,
            local_to_merged,
            sample_counts,
        }
    }
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

    /// A sample whose uniques carry a constant per-position quality sum, so
    /// merged sums are easy to predict.
    fn sample(uniques: &[(&str, u64, u32)]) -> Derep {
        Derep {
            uniques: uniques
                .iter()
                .map(|&(s, c, _)| (s.as_bytes().to_vec(), c))
                .collect(),
            quals: uniques.iter().map(|&(s, _, q)| vec![q; s.len()]).collect(),
            map: Vec::new(),
        }
    }

    fn pool(samples: &[Derep]) -> PooledDerep {
        pool_with(samples, PoolTiebreak::FirstSeen)
    }

    fn pool_with(samples: &[Derep], tiebreak: PoolTiebreak) -> PooledDerep {
        let mut pool = DerepPool::new();
        for s in samples {
            pool.add(s);
        }
        pool.finish_with(tiebreak)
    }

    /// #219: a unique that first appears late but is most abundant overall
    /// must still land at index 0, where the DADA2 loop expects the centre.
    #[test]
    fn pooled_most_abundant_is_at_index_zero() {
        let p = pool(&[
            sample(&[("AAAA", 5, 10), ("CCCC", 1, 10)]),
            sample(&[("GGGG", 4, 10), ("CCCC", 9, 10)]),
        ]);
        assert_eq!(p.seqs[0], b"CCCC");
        assert_eq!(p.abundance, vec![10, 5, 4]);
        assert_eq!(p.qual_sums[0], vec![20; 4]);
    }

    /// Each sample's map still points at the sequence it named, after the
    /// reorder, and its counts are the sample's own, not the pooled total.
    #[test]
    fn pooled_sample_maps_follow_the_reorder() {
        let samples = [
            sample(&[("AAAA", 5, 10), ("CCCC", 1, 10)]),
            sample(&[("GGGG", 4, 10), ("CCCC", 9, 10), ("TTTT", 2, 10)]),
        ];
        let p = pool(&samples);
        for (s, derep) in samples.iter().enumerate() {
            for (lu, (seq, count)) in derep.uniques.iter().enumerate() {
                let mu = p.local_to_merged[s][lu];
                assert_eq!(&p.seqs[mu], seq);
                assert_eq!(p.sample_counts[s][lu] as u64, *count);
            }
        }
    }

    /// Abundances and quality sums do not depend on the order samples are added
    /// in, and neither does the order when no abundances tie. (Ties follow
    /// input order under `first-seen`; see the test below.)
    #[test]
    fn pooled_table_is_independent_of_fold_order() {
        let a = sample(&[("AAAA", 5, 10), ("CCCC", 1, 7)]);
        let b = sample(&[("GGGG", 4, 3), ("CCCC", 9, 11), ("TTTT", 2, 5)]);
        let ab = pool(&[a, b]);
        let a = sample(&[("AAAA", 5, 10), ("CCCC", 1, 7)]);
        let b = sample(&[("GGGG", 4, 3), ("CCCC", 9, 11), ("TTTT", 2, 5)]);
        let ba = pool(&[b, a]);
        assert_eq!(ab.seqs, ba.seqs);
        assert_eq!(ab.abundance, ba.abundance);
        assert_eq!(ab.qual_sums, ba.qual_sums);
    }

    /// The `lexical` arm: ties by sequence, whatever the input order. Default
    /// before #260.
    #[test]
    fn pooled_lexical_arm_breaks_ties_by_sequence() {
        let p = pool_with(
            &[
                sample(&[("TTTT", 3, 10), ("GGGG", 2, 10)]),
                sample(&[("AAAA", 3, 10), ("CCCC", 2, 10)]),
            ],
            PoolTiebreak::Lexical,
        );
        let order: Vec<&[u8]> = p.seqs.iter().map(|s| s.as_slice()).collect();
        assert_eq!(order, [&b"AAAA"[..], b"TTTT", b"CCCC", b"GGGG"]);
    }

    /// Unset means `first-seen` (#260): the default `pool()` keeps input order
    /// among ties.
    #[test]
    fn pooled_default_breaks_ties_in_first_seen_order() {
        let p = pool(&[
            sample(&[("TTTT", 3, 10), ("GGGG", 2, 10)]),
            sample(&[("AAAA", 3, 10), ("CCCC", 2, 10)]),
        ]);
        let order: Vec<&[u8]> = p.seqs.iter().map(|s| s.as_slice()).collect();
        assert_eq!(order, [&b"TTTT"[..], b"AAAA", b"GGGG", b"CCCC"]);
    }

    /// R's `combineDereps2` rule: ties keep first appearance across samples,
    /// in input order, so reversing the samples reverses the tied pair.
    #[test]
    fn pooled_first_seen_ties_follow_input_order() {
        let a = || sample(&[("TTTT", 3, 10), ("GGGG", 2, 10)]);
        let b = || sample(&[("AAAA", 3, 10), ("CCCC", 2, 10)]);
        let seqs = |p: PooledDerep| -> Vec<Vec<u8>> { p.seqs };
        let ab = seqs(pool_with(&[a(), b()], PoolTiebreak::FirstSeen));
        let ba = seqs(pool_with(&[b(), a()], PoolTiebreak::FirstSeen));
        let s = |v: &[&str]| -> Vec<Vec<u8>> { v.iter().map(|x| x.as_bytes().to_vec()).collect() };
        assert_eq!(ab, s(&["TTTT", "AAAA", "GGGG", "CCCC"]));
        assert_eq!(ba, s(&["AAAA", "TTTT", "CCCC", "GGGG"]));
    }

    /// The tie-break never overrides abundance: index 0 is the most abundant
    /// under either rule (#219).
    #[test]
    fn pooled_first_seen_keeps_most_abundant_first() {
        let p = pool_with(
            &[
                sample(&[("AAAA", 5, 10), ("CCCC", 1, 10)]),
                sample(&[("GGGG", 4, 10), ("CCCC", 9, 10)]),
            ],
            PoolTiebreak::FirstSeen,
        );
        assert_eq!(p.seqs[0], b"CCCC");
        assert_eq!(p.abundance, vec![10, 5, 4]);
    }

    #[test]
    fn pool_tiebreak_parses_and_rejects() {
        assert_eq!(PoolTiebreak::parse(""), Ok(PoolTiebreak::FirstSeen));
        assert_eq!(PoolTiebreak::parse("lexical"), Ok(PoolTiebreak::Lexical));
        assert_eq!(
            PoolTiebreak::parse("first-seen"),
            Ok(PoolTiebreak::FirstSeen)
        );
        assert!(PoolTiebreak::parse("firstseen").is_err());
    }

    #[test]
    #[should_panic(expected = "overflows u32")]
    fn pooled_qual_sum_overflow_panics() {
        pool(&[sample(&[("AAAA", 1, u32::MAX)]), sample(&[("AAAA", 1, 1)])]);
    }
}
