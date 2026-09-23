# `summary`

Compute per-position quality metrics from a FASTQ file. Output is a JSON object
on stdout (or `--output`); optional metrics are off by default.

This is the dada2-rs counterpart of R DADA2's `plotQualityProfile()`, the survey
step near the start of the DADA2 SOP tutorial: look at raw read quality before
choosing truncation and filtering parameters. `summary` computes the statistics
and emits JSON; `scripts/plot_quality_profile.R` turns that JSON into the
familiar figure.

```bash
dada2-rs summary reads.fastq.gz --report -o summary.json
```

## Input

**`<INPUT>`** — a single FASTQ file, uncompressed or gzipped.

**`--sample-name`** — sample identifier written to the output JSON's `sample`
field. Defaults to the filename stem of the input FASTQ.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Metrics

All of these are optional; the base run reports per-position quality only.

**`--complexity`** — compute a per-read sequence-complexity histogram, matching
DADA2's `seqComplexity()` / `plotComplexity()`: the effective number of k-mers
per read, `exp(Shannon entropy)` over its k-mer counts.

**`--complexity-kmer-size`** (default 2) — k-mer size for the complexity
calculation, matching the DADA2 default. Only used with `--complexity`.

**`--complexity-bins`** (default 100) — number of histogram bins, spanning
`[0, 4^kmer_size]`. Matches `plotComplexity()`'s default. Only used with
`--complexity`.

**`--expected-error`** — compute per-position cumulative expected-error (EE)
metrics: `Σ 10^(-Q/10)` along each read, aggregated across reads into
mean/median/min/max/quartiles per position. This is the quantity
`filter-and-trim` thresholds with `--max-ee`, so it is the direct way to judge
`maxEE` and truncation choices before committing to them.

**`--ee-bins`** (default 200) — number of log-spaced histogram bins backing the
EE quantiles. Only used with `--expected-error`.

**`--binned-threshold`** (default 8) — the maximum number of distinct quality
values for the data to be judged *binned*. NovaSeq and NextSeq instruments
collapse Phred scores to a handful of levels; continuous Illumina and PacBio
HiFi data have dozens of distinct values, so the two regimes separate cleanly
well below this threshold.

Binned quality data changes how the error model should be fit — see
[Binned quality scores](../findings/binned-quality.md) and
`learn-errors --errfun binned-qual`.

## Performance

**`--threads`** (default 1) — threads for parallel processing. Does not affect
results.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON instead of pretty-printed.

## Diagnostics

**`--report`** — also print a human-readable metrics report to stderr: total
sequences, quality and EE ranges, and whether quality scores are binned along
with the observed levels. Stdout stays clean JSON, so this is safe to use in a
pipe.

## See also

- `summary-merge` — union per-sample summaries into a run-level report
- `scripts/plot_quality_profile.R` — reproduce DADA2's `plotQualityProfile()`
  figure from one or more `summary` JSONs
- [Binned quality scores](../findings/binned-quality.md)
- [`filter-and-trim`](filter-and-trim.md) — the step these metrics inform
