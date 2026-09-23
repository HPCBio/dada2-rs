# Plotting scripts

Four R scripts that render dada2-rs JSON into the figures DADA2 users know.
All take one or more JSON files and write a PDF; see
[conventions](index.md#conventions) for the shared flags.

## `plot_quality_profile.R`

The **first thing to run on a new dataset** — R's `plotQualityProfile()`, and
the step the DADA2 SOP opens with. A heatmap of cycle × quality with per-cycle
mean, quartile and cumulative-read lines, from which you choose truncation
length and `maxEE`.

```bash
dada2-rs summary reads.fastq.gz -o sample.json
Rscript scripts/plot_quality_profile.R --out=quality.pdf sample.json
```

```text
Rscript plot_quality_profile.R [--aggregate] [--out=plot.pdf]
                               [--width=8] [--height=5]
                               summary1.json [summary2.json ...]
```

Defaults to `quality_profile.pdf`. The per-cycle statistics are a direct port
and match R's Mean and Q50 **exactly** on the fixtures.

## `plot_complexity.R`

R's `plotComplexity()`: a per-file histogram of each read's effective
oligonucleotide number — `exp(Shannon entropy)` over its k-mer counts, in
`[1, 4^kmerSize]`. Low-complexity reads are the homopolymer and repeat junk you
may want filtered.

```bash
dada2-rs summary reads.fastq.gz --complexity -o sample.json
Rscript scripts/plot_complexity.R --out=complexity.pdf sample.json
```

Requires `summary --complexity`; without it the input carries no histogram.

One difference from R worth knowing: dada2-rs streams **all** reads and emits a
pre-binned histogram, where R's `plotComplexity()` subsamples via `FastqSampler`.
The script therefore draws the bins it is given rather than recomputing
`geom_histogram()` over raw values. The underlying per-read statistic is
bit-identical to `dada2:::seqComplexity()` on the fixtures.

## `plot_expected_error.R`

Per-position **cumulative** expected error — the running `Σ 10^(-Q/10)` that
`filter-and-trim` thresholds against `maxEE`. Where the curves cross the
reference lines at EE = 2, 3, 5 and 7, reads start being discarded, so this
reads directly as "what will `--max-ee` cost me at each truncation length".

```bash
dada2-rs summary reads.fastq.gz --expected-error -o sample.json
Rscript scripts/plot_expected_error.R --out=ee.pdf sample.json
```

```text
Rscript plot_expected_error.R [--out=plot.pdf] [--linear]
                              [--width=8] [--height=5]
                              summary1.json [summary2.json ...]
```

Log10 y-axis by default; `--linear` switches it. Requires
`summary --expected-error`.

This one has no DADA2 counterpart. It was inspired by Remi Maglione's
`Qual_vs_MaxEE` plot, but is driven by the **true per-read cumulative-EE
distribution** that `summary` aggregates — exact mean, min, max and quartiles —
rather than EE derived from the per-position mean-quality curve, which
understates the tail.

## `plot_errors.R`

R's `plotErrors()`: a 4×4 panel, one per nucleotide transition, of the fitted
error model.

- **points** — observed rate from the transition counts
- **black line** — estimated rate (`err_out`, the fit)
- **dashed line** — nominal rate fed into the final run (`err_in`)
- **red line** — theoretical Phred rate, `(1/3) × 10^(-Q/10)`

```bash
Rscript scripts/plot_errors.R err.json errors.pdf
Rscript scripts/plot_errors.R --help
```

### Beyond `plotErrors()`: weighting by observation mass

The addition worth knowing about. Points can be scaled and coloured by **how
much of the data actually sits at each quality score**, via `--mass=size`,
`colour`, `both` (the default) or `none`.

This matters on [binned quality scores](../findings/binned-quality.md) — NovaSeq,
Revio, i100 — where nearly all the mass sits on three or four Q values and the
visually prominent low-mass points carry almost no weight in the fit. On a
binned dataset the unweighted plot invites you to worry about points the model
barely sees.

Other options: `--out`, `--width`, `--height`, `--title`, and `--nti` / `--ntj`
to restrict which from/to nucleotides are drawn. `--help` lists the full set.
