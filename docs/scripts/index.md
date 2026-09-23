# Helper scripts

`scripts/` holds small R and Python helpers that sit either side of the
dada2-rs pipeline: they turn its JSON into the figures DADA2 users expect, carry
an R-fitted error model in, or roll per-sample counts up into a table. None of
them is required to run a pipeline, and none of them does inference.

They are documented here and only here — a subcommand page links to a script
rather than explaining it, so there is one place to keep current.

| Script | Language | What it is for |
|---|---|---|
| [`plot_quality_profile.R`](plotting.md#plot_quality_profiler) | R | quality heatmap, R's `plotQualityProfile()` |
| [`plot_complexity.R`](plotting.md#plot_complexityr) | R | read-complexity histogram, R's `plotComplexity()` |
| [`plot_expected_error.R`](plotting.md#plot_expected_errorr) | R | cumulative expected error along the read |
| [`plot_errors.R`](plotting.md#plot_errorsr) | R | error-model panel, R's `plotErrors()` |
| [`learnerrors_to_dada2rs.R`](pipeline-helpers.md#learnerrors_to_dada2rsr) | R | R `learnErrors()` `.rds` → dada2-rs error model |
| [`track_reads.py`](pipeline-helpers.md#track_readspy) | Python | per-sample read counts across every step |
| [`examples/import_mia.R`](pipeline-helpers.md#examplesimport_miar) | R | template for loading the final tables into `mia` |

## What is not here

`dev/` holds maintainer tooling — the benchmark harness, the concordance
runners, release and build-config scripts. It is deliberately undocumented on
this site: **everything in `scripts/` is documented, nothing in `dev/` is.**
That split is what keeps this page from drifting into a list of things nobody
outside the repository runs.

The install recipe had encoded half of that rule for a while before it was
written down: it globs `scripts/*.R` and `scripts/*.py`, so the two shell
scripts that used to sit alongside them were never installed anyway. Moving
them to `dev/maintenance/` made the boundary explicit rather than incidental.

## Running them

From a checkout, invoke them directly — that is the form used throughout these
pages:

```bash
Rscript scripts/plot_errors.R err.json errors.pdf
```

`make install` / `just install` also copies them onto `PATH`, namespaced as
`dada2-rs-<name>` with underscores becoming hyphens, so the same call becomes:

```bash
dada2-rs-plot-errors err.json errors.pdf
```

| Script | Installed as |
|---|---|
| `plot_quality_profile.R` | `dada2-rs-plot-quality-profile` |
| `plot_complexity.R` | `dada2-rs-plot-complexity` |
| `plot_expected_error.R` | `dada2-rs-plot-expected-error` |
| `plot_errors.R` | `dada2-rs-plot-errors` |
| `learnerrors_to_dada2rs.R` | `dada2-rs-learnerrors-to-dada2rs` |
| `track_reads.py` | `dada2-rs-track-reads` |

`cargo install dada2-rs` installs **only the binary** — it has no way to place
these — so a crates.io install needs a checkout for the helpers. See
[Installation](../installation.md#installing-onto-path).

`examples/import_mia.R` is not installed; it is a template to copy.

## Requirements

| Script | Needs |
|---|---|
| all four plotting scripts | `jsonlite`, `ggplot2` |
| `plot_errors.R` | also `optparse` |
| `learnerrors_to_dada2rs.R` | `jsonlite` |
| `track_reads.py` | Python 3, standard library only |
| `examples/import_mia.R` | `mia`, `SummarizedExperiment` (Bioconductor) |

```r
install.packages(c("jsonlite", "ggplot2", "optparse"))
```

## Conventions

The four plotting scripts share an interface: they take one or more `summary`
or `learn-errors` JSON files and write a PDF.

```bash
Rscript scripts/<script>.R [--out=plot.pdf] [--width=8] [--height=5] input.json ...
```

`--aggregate` (on `plot_quality_profile.R` and `plot_complexity.R`) pools every
input into a single panel instead of faceting per file — the equivalent of
passing a vector of files to the R function and setting `aggregate = TRUE`.

## Provenance

The plotting statistics are ports of methods from the DADA2 R package by
Benjamin Callahan — `plotQualityProfile()`, `plotComplexity()`,
`seqComplexity()` and `plotErrors()`, all in `dada2/R/plot-methods.R` and
`dada2/R/filter.R`. Each script's header cites the original and records how
closely it was validated against it. Credit for the methods is his.
