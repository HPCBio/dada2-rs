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

/// Per-position EE histogram is binned on a log10 scale over
/// `[10^EE_LOG_MIN, 10^EE_LOG_MAX]`. EE below the low edge clamps to bin 0,
/// above the high edge to the last bin. The range comfortably brackets the
/// standard maxEE cutoffs (2/3/5/7) with fine low-end resolution.
const EE_LOG_MIN: f64 = -4.0;
const EE_LOG_MAX: f64 = 3.0;

/// Configuration for the optional per-position expected-error (EE) metrics. EE
/// at a position is the read's *cumulative* expected error up to that base:
/// `Σ 10^(-Q_i/10)` for `i ≤ pos` — the same quantity `filter-and-trim`
/// thresholds against `maxEE`. We aggregate the per-read cumulative EE across
/// reads at each position: exact mean/min/max, plus quantiles from a per-position
/// log-binned histogram.
#[derive(Clone, Copy)]
pub struct ExpectedErrorConfig {
    pub bins: usize,
}

/// What extra per-read metrics `summary` should compute, beyond the per-position
/// quality aggregates it always produces.
#[derive(Clone, Copy, Default)]
pub struct SummaryConfig {
    pub complexity: Option<ComplexityConfig>,
    pub expected_error: Option<ExpectedErrorConfig>,
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

/// Whether a quality-value distribution with `distinct_len` distinct levels is
/// judged "binned": 2..=`binned_threshold` distinct levels. Shared by the
/// per-sample `summary` and the run-level `summary-merge` so both use one rule.
pub fn judge_binned(distinct_len: usize, binned_threshold: usize) -> bool {
    distinct_len >= 2 && distinct_len <= binned_threshold
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
    /// Optional per-position cumulative-EE metrics. `None` unless requested.
    expected_error: Option<ExpectedErrorConfig>,
    /// Per-position accumulators over per-read cumulative EE. Empty when disabled.
    ee_sum: Vec<f64>,
    ee_min: Vec<f64>,
    ee_max: Vec<f64>,
    /// `ee_hist[pos][bin]` = reads whose cumulative EE at `pos` falls in log-bin
    /// `bin` (see `EE_LOG_MIN`/`EE_LOG_MAX`). Used for per-position quantiles.
    ee_hist: Vec<Vec<u64>>,
}

impl QualitySummary {
    pub fn with_config(config: SummaryConfig) -> Self {
        let complexity_hist = match config.complexity {
            Some(cfg) => vec![0; cfg.bins],
            None => Vec::new(),
        };
        Self {
            total_reads: 0,
            sums: Vec::new(),
            counts: Vec::new(),
            hist: Vec::new(),
            complexity: config.complexity,
            complexity_hist,
            expected_error: config.expected_error,
            ee_sum: Vec::new(),
            ee_min: Vec::new(),
            ee_max: Vec::new(),
            ee_hist: Vec::new(),
        }
    }

    /// Log-bin index for a cumulative EE value, clamped to `[0, bins)`.
    fn ee_bin(ee: f64, bins: usize) -> usize {
        if ee <= 0.0 {
            return 0;
        }
        let frac = (ee.log10() - EE_LOG_MIN) / (EE_LOG_MAX - EE_LOG_MIN);
        ((frac * bins as f64).floor() as isize).clamp(0, bins as isize - 1) as usize
    }

    /// Representative EE value at the center of log-bin `bin`.
    fn ee_bin_center(bin: usize, bins: usize) -> f64 {
        let log = EE_LOG_MIN + (bin as f64 + 0.5) / bins as f64 * (EE_LOG_MAX - EE_LOG_MIN);
        10f64.powf(log)
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
        if let Some(cfg) = self.expected_error {
            let n = quality.len();
            if n > self.ee_sum.len() {
                self.ee_sum.resize(n, 0.0);
                self.ee_min.resize(n, f64::INFINITY);
                self.ee_max.resize(n, 0.0);
                self.ee_hist.resize(n, vec![0; cfg.bins]);
            }
            // Cumulative expected error: Σ 10^(-Q/10) up to and including pos.
            let mut ee = 0.0;
            for (i, &q) in quality.iter().enumerate() {
                let q_phred = ((q as i16) - (phred_offset as i16)).max(0) as f64;
                ee += 10f64.powf(-q_phred / 10.0);
                self.ee_sum[i] += ee;
                if ee < self.ee_min[i] {
                    self.ee_min[i] = ee;
                }
                if ee > self.ee_max[i] {
                    self.ee_max[i] = ee;
                }
                self.ee_hist[i][Self::ee_bin(ee, cfg.bins)] += 1;
            }
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
        // Expected-error accumulators.
        let ee_len = self.ee_sum.len().max(other.ee_sum.len());
        if ee_len > 0 {
            let bins = self
                .expected_error
                .or(other.expected_error)
                .map(|c| c.bins)
                .unwrap_or(0);
            self.ee_sum.resize(ee_len, 0.0);
            self.ee_min.resize(ee_len, f64::INFINITY);
            self.ee_max.resize(ee_len, 0.0);
            self.ee_hist.resize(ee_len, vec![0; bins]);
            for i in 0..other.ee_sum.len() {
                self.ee_sum[i] += other.ee_sum[i];
                if other.ee_min[i] < self.ee_min[i] {
                    self.ee_min[i] = other.ee_min[i];
                }
                if other.ee_max[i] > self.ee_max[i] {
                    self.ee_max[i] = other.ee_max[i];
                }
                for (b, &n) in other.ee_hist[i].iter().enumerate() {
                    self.ee_hist[i][b] += n;
                }
            }
        }
        if self.expected_error.is_none() {
            self.expected_error = other.expected_error;
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

    /// Distinct quality values observed anywhere in the file, ascending, each
    /// with its total read-base count summed over all positions. Derived from
    /// the per-position histogram — no extra pass over the reads.
    pub fn quality_value_counts(&self) -> Vec<(u8, u64)> {
        let mut totals = [0u64; MAX_QUAL + 1];
        for row in &self.hist {
            for (q, &n) in row.iter().enumerate() {
                totals[q] += n;
            }
        }
        totals
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n > 0)
            .map(|(q, &n)| (q as u8, n))
            .collect()
    }

    /// Summary of the quality-value distribution used to judge whether the
    /// data is binned. `binned_threshold` is the maximum number of distinct Q
    /// values for the data to be judged binned (real binned schemes have 2–8
    /// levels; continuous Illumina/HiFi data has dozens). The full distinct set
    /// and counts are always returned so borderline cases can be judged by eye.
    pub fn quality_metrics(&self, binned_threshold: usize) -> QualityMetrics {
        let counts = self.quality_value_counts();
        let distinct: Vec<u8> = counts.iter().map(|&(q, _)| q).collect();
        let value_counts: Vec<u64> = counts.iter().map(|&(_, n)| n).collect();
        let min_quality = distinct.first().copied().unwrap_or(0);
        let max_quality = distinct.last().copied().unwrap_or(0);
        let is_binned = judge_binned(distinct.len(), binned_threshold);
        let binned_scores = is_binned.then(|| distinct.clone());
        QualityMetrics {
            distinct_quality_values: distinct,
            quality_value_counts: value_counts,
            min_quality,
            max_quality,
            is_binned,
            binned_scores,
        }
    }

    /// The per-read complexity histogram, if it was requested.
    /// Returns `(kmer_size, bins, counts)` where `counts[bin]` covers the
    /// effective-k-mer-count range `[bin, bin+1) · 4^kmer_size / bins`.
    pub fn complexity_histogram(&self) -> Option<(u8, usize, &[u64])> {
        self.complexity
            .map(|cfg| (cfg.kmer_size, cfg.bins, self.complexity_hist.as_slice()))
    }

    /// Per-position cumulative expected-error metrics, if requested. Each vector
    /// is indexed by zero-based position. `mean`/`min`/`max` are exact; the
    /// quantiles (`q25`, `median`, `q75`) are read off the per-position log-binned
    /// histogram (bin-center resolution).
    pub fn expected_error_metrics(&self) -> Option<ExpectedErrorMetrics> {
        let cfg = self.expected_error?;
        let npos = self.ee_sum.len();
        let mut mean = vec![0.0; npos];
        let mut min = vec![0.0; npos];
        let mut max = vec![0.0; npos];
        let mut q25 = vec![0.0; npos];
        let mut median = vec![0.0; npos];
        let mut q75 = vec![0.0; npos];
        for i in 0..npos {
            let count = self.counts[i];
            if count == 0 {
                continue;
            }
            mean[i] = self.ee_sum[i] / count as f64;
            min[i] = self.ee_min[i];
            max[i] = self.ee_max[i];
            q25[i] = self.ee_quantile(i, 0.25, cfg.bins);
            median[i] = self.ee_quantile(i, 0.50, cfg.bins);
            q75[i] = self.ee_quantile(i, 0.75, cfg.bins);
        }
        Some(ExpectedErrorMetrics {
            mean,
            min,
            max,
            q25,
            median,
            q75,
        })
    }

    /// Quantile `q` of the cumulative-EE distribution at position `i`, taken from
    /// the log-binned histogram (returns the center EE of the crossing bin).
    fn ee_quantile(&self, i: usize, q: f64, bins: usize) -> f64 {
        let row = &self.ee_hist[i];
        let total: u64 = row.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let target = q * total as f64;
        let mut cum = 0u64;
        for (bin, &n) in row.iter().enumerate() {
            cum += n;
            if cum as f64 >= target {
                return Self::ee_bin_center(bin, bins);
            }
        }
        Self::ee_bin_center(bins - 1, bins)
    }
}

/// Summary of the quality-value distribution (see `quality_metrics`), used to
/// answer whether a FASTQ's quality scores are binned and to what levels.
pub struct QualityMetrics {
    /// Distinct Q values observed anywhere in the file, ascending.
    pub distinct_quality_values: Vec<u8>,
    /// Total read-base count for each distinct value (parallel to above).
    pub quality_value_counts: Vec<u64>,
    pub min_quality: u8,
    pub max_quality: u8,
    /// True when the number of distinct Q values is small enough to be binned.
    pub is_binned: bool,
    /// The bin levels (the distinct set) when `is_binned`, else `None`.
    pub binned_scores: Option<Vec<u8>>,
}

/// Per-position cumulative expected-error metrics (see `expected_error_metrics`).
pub struct ExpectedErrorMetrics {
    pub mean: Vec<f64>,
    pub min: Vec<f64>,
    pub max: Vec<f64>,
    pub q25: Vec<f64>,
    pub median: Vec<f64>,
    pub q75: Vec<f64>,
}

/// Records per thread per chunk — total chunk size scales with thread count.
const RECORDS_PER_THREAD: usize = 10_000;

pub fn process<R: io::Read>(
    reader: R,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    config: SummaryConfig,
) -> io::Result<QualitySummary> {
    let chunk_size = RECORDS_PER_THREAD * pool.current_num_threads();
    let buf = BufReader::new(reader);
    let mut fastq_reader = fastq::io::Reader::new(buf);
    let mut overall = QualitySummary::with_config(config);

    // Only retain the sequence when complexity is requested — otherwise it's
    // dead weight in the chunk buffer.
    let want_seq = config.complexity.is_some();

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
                    || QualitySummary::with_config(config),
                    |mut acc, (sequence, quality)| {
                        acc.add_record(sequence, quality, phred_offset);
                        acc
                    },
                )
                .reduce(
                    || QualitySummary::with_config(config),
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

#[cfg(test)]
mod tests {
    use super::*;

    // Phred+33: 'I' = Q40 (err 1e-4), '+' = Q10 (err 0.1).
    #[test]
    fn cumulative_ee_is_running_sum_of_per_base_error() {
        let cfg = SummaryConfig {
            complexity: None,
            expected_error: Some(ExpectedErrorConfig { bins: 200 }),
        };
        let mut s = QualitySummary::with_config(cfg);
        // Two reads so mean/min/max differ from any single read.
        s.add_record(b"ACGT", b"IIII", 33); // EE each base 1e-4
        s.add_record(b"ACGT", b"++++", 33); // EE each base 1e-1
        let m = s.expected_error_metrics().unwrap();
        // Position 0: per-base EE only.
        assert!((m.min[0] - 1e-4).abs() < 1e-12);
        assert!((m.max[0] - 1e-1).abs() < 1e-12);
        // Mean of cumulative EE at last position = mean of the two read totals.
        let expected_mean = (4.0 * 1e-4 + 4.0 * 1e-1) / 2.0;
        assert!((m.mean[3] - expected_mean).abs() < 1e-12);
        // Cumulative is monotonic non-decreasing within the max curve.
        assert!(m.max[3] >= m.max[0]);
        assert!((m.max[3] - 0.4).abs() < 1e-12);
    }

    /// Phred+33 char for a given Q value.
    fn qchar(q: u8) -> u8 {
        q + 33
    }

    #[test]
    fn judge_binned_rule() {
        assert!(!judge_binned(0, 8)); // no data
        assert!(!judge_binned(1, 8)); // constant quality — not binned
        assert!(judge_binned(2, 8));
        assert!(judge_binned(8, 8));
        assert!(!judge_binned(9, 8)); // above threshold — continuous
    }

    #[test]
    fn binned_data_detected() {
        let mut s = QualitySummary::with_config(SummaryConfig::default());
        // Only four distinct quality levels, NovaSeq-style {2,12,23,37}.
        let quals = [qchar(2), qchar(12), qchar(23), qchar(37)];
        s.add_record(b"ACGT", &quals, 33);
        s.add_record(b"ACGT", &[qchar(37), qchar(37), qchar(2), qchar(12)], 33);
        let m = s.quality_metrics(8);
        assert!(m.is_binned);
        assert_eq!(m.distinct_quality_values, vec![2, 12, 23, 37]);
        assert_eq!(m.binned_scores, Some(vec![2, 12, 23, 37]));
        assert_eq!(m.min_quality, 2);
        assert_eq!(m.max_quality, 37);
    }

    #[test]
    fn continuous_data_not_binned() {
        let mut s = QualitySummary::with_config(SummaryConfig::default());
        // Ten distinct quality levels — above the default threshold of 8.
        let quals: Vec<u8> = (20u8..30).map(qchar).collect();
        s.add_record(&vec![b'A'; quals.len()], &quals, 33);
        let m = s.quality_metrics(8);
        assert!(!m.is_binned);
        assert_eq!(m.binned_scores, None);
        assert_eq!(m.distinct_quality_values.len(), 10);
    }

    #[test]
    fn quality_value_counts_sum_matches_total_bases() {
        let mut s = QualitySummary::with_config(SummaryConfig::default());
        s.add_record(b"ACGT", &[qchar(10), qchar(20), qchar(20), qchar(30)], 33);
        s.add_record(b"AC", &[qchar(20), qchar(30)], 33);
        let total_bases: u64 = s.quality_value_counts().iter().map(|&(_, n)| n).sum();
        assert_eq!(total_bases, 6);
    }

    #[test]
    fn complexity_of_a_homopolymer_is_one() {
        // A single 2-mer ("AA") repeated -> one distinct k-mer -> exp(0) = 1.
        assert!((seq_complexity(b"AAAAAA", 2) - 1.0).abs() < 1e-12);
        // Equal counts of all four bases' 2-mers approach higher diversity.
        assert!(seq_complexity(b"ACGTACGTACGT", 2) > 3.0);
        // Read shorter than k, or all-N, falls back to 1.0 (matches R sindex).
        assert_eq!(seq_complexity(b"A", 2), 1.0);
        assert_eq!(seq_complexity(b"NNNN", 2), 1.0);
    }
}
