//! Quality summaries: `summary`, `summary-merge`.

use std::{fs::File, io};

use flate2::read::MultiGzDecoder;

use dada2_rs::{cli, misc, summary};

use misc::{Tagged, read_tagged_json};
use serde::Serialize;
use summary::{ComplexityConfig, ExpectedErrorConfig, SummaryConfig, judge_binned, process};

use super::common::{check_input_paths, fastq_stem};
use misc::WithPath;

pub(crate) fn run_summary(args: cli::SummaryArgs) -> io::Result<()> {
    let cli::SummaryArgs {
        input,
        sample_name,
        phred_offset,
        threads,
        output,
        compact,
        complexity,
        complexity_kmer_size,
        complexity_bins,
        expected_error,
        ee_bins,
        binned_threshold,
        report,
    } = args;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    let summary_cfg = SummaryConfig {
        complexity: complexity.then_some(ComplexityConfig {
            kmer_size: complexity_kmer_size,
            bins: complexity_bins,
        }),
        expected_error: expected_error.then_some(ExpectedErrorConfig { bins: ee_bins }),
    };

    let summary = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
        process(
            MultiGzDecoder::new(File::open(&input).with_path(&input)?),
            phred_offset,
            &pool,
            summary_cfg,
        )?
    } else {
        process(
            File::open(&input).with_path(&input)?,
            phred_offset,
            &pool,
            summary_cfg,
        )?
    };

    /// `complexity_histogram[bin]` = number of reads whose effective
    /// k-mer count (`exp(Shannon entropy)`, range `[1, 4^kmer_size]`)
    /// falls in equal-width `bin` over `[0, 4^kmer_size]`.
    #[derive(Serialize)]
    struct ComplexityOutput {
        kmer_size: u8,
        bins: usize,
        histogram: Vec<u64>,
    }

    /// Per-position cumulative expected error aggregated across reads.
    /// Each array is indexed by zero-based position; `mean`/`min`/`max`
    /// are exact, the quantiles come from a log-binned histogram.
    #[derive(Serialize)]
    struct ExpectedErrorOutput {
        mean: Vec<f64>,
        min: Vec<f64>,
        q25: Vec<f64>,
        median: Vec<f64>,
        q75: Vec<f64>,
        max: Vec<f64>,
    }

    /// Quality-value distribution metrics: whether the FASTQ's quality
    /// scores are binned and, if so, to what levels. Always present.
    #[derive(Serialize)]
    struct QualityMetricsOutput {
        distinct_quality_values: Vec<u8>,
        quality_value_counts: Vec<u64>,
        min_quality: u8,
        max_quality: u8,
        is_binned: bool,
        #[serde(skip_serializing_if = "Option::is_none")]
        binned_scores: Option<Vec<u8>>,
    }

    #[derive(Serialize)]
    struct SummaryOutput {
        sample: String,
        total_reads: u64,
        mean_quality_per_position: Vec<f64>,
        reads_per_position: Vec<u64>,
        max_quality: usize,
        /// `quality_histogram[pos][q]` = count of reads with quality `q`
        /// at zero-based cycle `pos`. Inner length is `max_quality + 1`.
        quality_histogram: Vec<Vec<u64>>,
        /// Per-read sequence-complexity histogram; present only when
        /// `--complexity` was requested.
        #[serde(skip_serializing_if = "Option::is_none")]
        complexity: Option<ComplexityOutput>,
        /// Per-position cumulative expected error; present only when
        /// `--expected-error` was requested.
        #[serde(skip_serializing_if = "Option::is_none")]
        expected_error: Option<ExpectedErrorOutput>,
        /// Quality-value distribution + binning verdict. Always present.
        quality_metrics: QualityMetricsOutput,
    }

    let sample = sample_name.unwrap_or_else(|| fastq_stem(&input));
    let (max_quality, quality_histogram) = summary.quality_histogram();
    let complexity_out = summary
        .complexity_histogram()
        .map(|(kmer_size, bins, hist)| ComplexityOutput {
            kmer_size,
            bins,
            histogram: hist.to_vec(),
        });
    let expected_error_out = summary
        .expected_error_metrics()
        .map(|m| ExpectedErrorOutput {
            mean: m.mean,
            min: m.min,
            q25: m.q25,
            median: m.median,
            q75: m.q75,
            max: m.max,
        });
    let metrics = summary.quality_metrics(binned_threshold);

    if report {
        let coverage = summary.reads_per_position();
        let read_len = coverage.len();
        let binned = match &metrics.binned_scores {
            Some(levels) => {
                let levels: Vec<String> = levels.iter().map(|q| q.to_string()).collect();
                format!("yes ({} levels: {})", levels.len(), levels.join(", "))
            }
            None => "no".to_string(),
        };
        eprintln!("=== FASTQ quality metrics: {sample} ===");
        eprintln!("total sequences:   {}", summary.total_reads);
        eprintln!("read length:       {read_len} cycles");
        eprintln!(
            "quality range:     Q{}–Q{} ({} distinct values)",
            metrics.min_quality,
            metrics.max_quality,
            metrics.distinct_quality_values.len()
        );
        eprintln!("binned:            {binned}");
        if let Some(ee) = &expected_error_out
            && let (Some(&min), Some(&max)) = (ee.min.last(), ee.max.last())
        {
            eprintln!("expected error:    {min:.3}–{max:.3} (cumulative, final position)");
        }
    }

    let quality_metrics = QualityMetricsOutput {
        distinct_quality_values: metrics.distinct_quality_values,
        quality_value_counts: metrics.quality_value_counts,
        min_quality: metrics.min_quality,
        max_quality: metrics.max_quality,
        is_binned: metrics.is_binned,
        binned_scores: metrics.binned_scores,
    };

    let out = SummaryOutput {
        sample,
        total_reads: summary.total_reads,
        mean_quality_per_position: summary.mean_quality_per_position(),
        reads_per_position: summary.reads_per_position().to_vec(),
        max_quality,
        quality_histogram,
        complexity: complexity_out,
        expected_error: expected_error_out,
        quality_metrics,
    };

    let tagged = Tagged::new("summary", out);
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_summary_merge(args: cli::SummaryMergeArgs) -> io::Result<()> {
    let cli::SummaryMergeArgs {
        inputs,
        output,
        compact,
        binned_threshold,
        report,
        expected_bins,
    } = args;
    check_input_paths("input", &inputs)?;
    use std::collections::BTreeMap;

    /// Only the fields we need out of a `summary` JSON.
    #[derive(serde::Deserialize)]
    struct SummaryJson {
        #[serde(default)]
        sample: Option<String>,
        #[serde(default)]
        total_reads: u64,
        #[serde(default)]
        quality_metrics: Option<QualityMetricsJson>,
    }
    #[derive(serde::Deserialize)]
    struct QualityMetricsJson {
        distinct_quality_values: Vec<u8>,
        quality_value_counts: Vec<u64>,
    }

    if inputs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "summary-merge: no input summary JSON files given",
        ));
    }

    // Union of Q -> summed base count across all samples.
    let mut union: BTreeMap<u8, u64> = BTreeMap::new();
    let mut total_reads: u64 = 0;
    let mut n_samples: u64 = 0;
    // Per-sample distinct sets, to flag samples missing run-level bins.
    let mut per_sample: Vec<(String, Vec<u8>)> = Vec::new();

    for path in &inputs {
        let parsed: SummaryJson = read_tagged_json(path, &["summary"]).with_path(path)?;
        let qm = match parsed.quality_metrics {
            Some(qm) => qm,
            None => {
                eprintln!(
                    "warning: {} has no quality_metrics block; skipping",
                    path.display()
                );
                continue;
            }
        };
        if qm.distinct_quality_values.len() != qm.quality_value_counts.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: quality_metrics distinct/counts length mismatch",
                    path.display()
                ),
            ));
        }
        let sample = parsed.sample.unwrap_or_else(|| fastq_stem(path));
        for (&q, &c) in qm
            .distinct_quality_values
            .iter()
            .zip(qm.quality_value_counts.iter())
        {
            *union.entry(q).or_insert(0) += c;
        }
        per_sample.push((sample, qm.distinct_quality_values));
        total_reads += parsed.total_reads;
        n_samples += 1;
    }

    if union.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "summary-merge: no usable quality_metrics found in any input",
        ));
    }

    let distinct: Vec<u8> = union.keys().copied().collect();
    let counts: Vec<u64> = union.values().copied().collect();
    let run_set: std::collections::BTreeSet<u8> = distinct.iter().copied().collect();
    let min_quality = *distinct.first().unwrap();
    let max_quality = *distinct.last().unwrap();
    let is_binned = judge_binned(distinct.len(), binned_threshold);
    let binned_scores = is_binned.then(|| distinct.clone());

    // Per-sample flags: which run-level bins this sample never emitted.
    #[derive(Serialize)]
    struct SampleFlag {
        sample: String,
        distinct_quality_values: Vec<u8>,
        missing_from_run: Vec<u8>,
    }
    let samples: Vec<SampleFlag> = per_sample
        .into_iter()
        .map(|(sample, s_distinct)| {
            let s_set: std::collections::BTreeSet<u8> = s_distinct.iter().copied().collect();
            let missing_from_run = run_set.difference(&s_set).copied().collect::<Vec<u8>>();
            SampleFlag {
                sample,
                distinct_quality_values: s_distinct,
                missing_from_run,
            }
        })
        .collect();

    // Optional validation against a declared bin scheme.
    #[derive(Serialize)]
    struct ExpectedValidation {
        consistent: bool,
        observed_not_in_expected: Vec<u8>,
        expected_not_observed: Vec<u8>,
    }
    let expected_validation = (!expected_bins.is_empty()).then(|| {
        let exp: std::collections::BTreeSet<u8> = expected_bins.iter().copied().collect();
        let observed_not_in_expected = run_set.difference(&exp).copied().collect::<Vec<u8>>();
        let expected_not_observed = exp.difference(&run_set).copied().collect::<Vec<u8>>();
        ExpectedValidation {
            consistent: observed_not_in_expected.is_empty(),
            observed_not_in_expected,
            expected_not_observed,
        }
    });

    if report {
        let binned = match &binned_scores {
            Some(levels) => {
                let levels: Vec<String> = levels.iter().map(|q| q.to_string()).collect();
                format!("yes ({} levels: {})", levels.len(), levels.join(", "))
            }
            None => "no".to_string(),
        };
        eprintln!("=== Run-level FASTQ quality metrics ===");
        eprintln!("samples:           {n_samples}");
        eprintln!("total sequences:   {total_reads}");
        eprintln!(
            "quality range:     Q{min_quality}–Q{max_quality} ({} distinct values)",
            distinct.len()
        );
        eprintln!("binned:            {binned}");
        if let Some(v) = &expected_validation {
            if v.consistent {
                eprintln!("expected scheme:   CONSISTENT");
            } else {
                let bad: Vec<String> = v
                    .observed_not_in_expected
                    .iter()
                    .map(|q| format!("Q{q}"))
                    .collect();
                eprintln!(
                    "expected scheme:   VIOLATION (observed {} not in declared scheme)",
                    bad.join(", ")
                );
            }
            if !v.expected_not_observed.is_empty() {
                let missing: Vec<String> = v
                    .expected_not_observed
                    .iter()
                    .map(|q| format!("Q{q}"))
                    .collect();
                eprintln!(
                    "                   declared but unobserved: {}",
                    missing.join(", ")
                );
            }
        }
        let odd: Vec<&SampleFlag> = samples
            .iter()
            .filter(|s| !s.missing_from_run.is_empty())
            .collect();
        if !odd.is_empty() {
            eprintln!("samples missing run-level bins:");
            for s in odd {
                let m: Vec<String> = s.missing_from_run.iter().map(|q| format!("Q{q}")).collect();
                eprintln!("  {}: missing {}", s.sample, m.join(", "));
            }
        }
    }

    #[derive(Serialize)]
    struct QualityMetricsOut {
        distinct_quality_values: Vec<u8>,
        quality_value_counts: Vec<u64>,
        min_quality: u8,
        max_quality: u8,
        is_binned: bool,
        #[serde(skip_serializing_if = "Option::is_none")]
        binned_scores: Option<Vec<u8>>,
    }
    #[derive(Serialize)]
    struct MergeOut {
        n_samples: u64,
        total_reads: u64,
        quality_metrics: QualityMetricsOut,
        #[serde(skip_serializing_if = "Vec::is_empty")]
        expected_bins: Vec<u8>,
        #[serde(skip_serializing_if = "Option::is_none")]
        expected_validation: Option<ExpectedValidation>,
        samples: Vec<SampleFlag>,
    }

    let out = MergeOut {
        n_samples,
        total_reads,
        quality_metrics: QualityMetricsOut {
            distinct_quality_values: distinct,
            quality_value_counts: counts,
            min_quality,
            max_quality,
            is_binned,
            binned_scores,
        },
        expected_bins,
        expected_validation,
        samples,
    };

    let tagged = Tagged::new("summary-merge", out);
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}
