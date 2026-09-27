# `learn-errors`

Fit an error model from FASTQ or derep/sample JSON files.

```bash
dada2-rs learn-errors derep/*.json.gz --threads 24 -o err.json
```

Reads one or more FASTQ files (dereplicated on the fly) or pre-computed
derep/sample JSON files (`.json` / `.json.gz`, as written by `derep` and
`sample`), accumulates them up to `--nbases` total bases, then iteratively runs
the DADA2 algorithm and re-fits the chosen error model until self-consistency.

Output is a JSON object with three flat 16 × nq matrices:

| Field | Meaning |
|---|---|
| `trans` | accumulated transition counts |
| `err_in` | error rates used in the final DADA run |
| `err_out` | error rates estimated from `trans` |

`dada` consumes `err_out` by default.

!!! warning "Cross-sample diversity: `--nbases` accumulates whole samples"
    Files are taken **whole** — in the supplied order, or shuffled with
    `--randomize` — and each contributes all of its bases until the running
    total reaches `--nbases`. Reads are never subsampled *within* a file.

    With modern deep runs a single sample can supply the entire budget (400k
    reads × 250 bp = 100M bases = the default), so the error model may be
    learned from the diversity of just one or a few samples. `--randomize` only
    shuffles which samples are drawn first; it does not guarantee representation
    across the run.

    Raise `--nbases`, pass a hand-picked set of inputs, or pre-`sample` each
    file to spread learning across samples. Tracked in
    [issue #68](https://github.com/HPCBio/dada2-rs/issues/68).

## Input

**`<INPUT>...`** — FASTQ (`.fastq` / `.fastq.gz` / `.fq` / `.fq.gz`) or
derep/sample JSON (`.json` / `.json.gz`) files.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Error model

### Sampling

**`--nbases`** (default 1e8) — stop after accumulating at least this many total
bases. See the caveat above; the model is not necessarily converged at the
default, and moving `--nbases` can shift the fit more than a `--kdist-cutoff`
change does. See
[learn-errors --nbases convergence](../findings/learn-errors-nbases-convergence.md).

**`--randomize`** — process input files in random order. Shuffles sample order
only.

**`--seed`** — RNG seed for reproducible `--randomize`.

### Fitting function

**`--errfun`** (default `loess`) — one of `loess`, `noqual`, `binned-qual`,
`pacbio`, `external`.

`loess` vs R DADA2: the native Rust LOESS now matches R's `loess()` to
double-precision round-off on **both** surfaces — ~1e-14 or better in every
regime `dev/loess-oracle` tests, including the sparse-quality case
(issues #213, #215, #217). Since #205 the default surface is `interpolate`,
which is what R's `loessErrfun` uses, so a stock `learn-errors` run now agrees
with R's error model without any extra flags.

Use `--loess-surface direct` for the historical behaviour. Direct is equally
exact against `loess(surface = "direct")` — but R DADA2 never calls that, so
it differs from R by the whole direct-vs-interpolate gap: on a 95-sample
pooled PacBio run that was the largest remaining term in the error model
(median 1.088e-03 between surfaces, versus 2.599e-04 for the k-mer screen).

For bit-for-bit parity with R's `loessErrfun`:

```bash
dada2-rs learn-errors ... \
  --errfun external \
  --errfun-cmd "Rscript examples/external_errfun/loess_reference.R"
```

See [issue #14](https://github.com/HPCBio/dada2-rs/issues/14) for the full
decomposition, and [the LOESS findings
page](../findings/loess-error-model-correctness.md) for why the error model is
worth this much attention.

To skip fitting entirely and use a model R has **already** learned, convert it
with [`learnerrors_to_dada2rs.R`](../scripts/pipeline-helpers.md#learnerrors_to_dada2rsr)
and pass the result to `dada --error-model`.

Both routes, the wire contract for writing your own script, and the caveats are
in [Using an external error
model](../walkthroughs/external-error-models.md).

Use `binned-qual` on NovaSeq/NextSeq-style data where Phred scores are collapsed
to a handful of levels — `summary --report` will tell you whether that is the
case. See [Binned quality scores](../findings/binned-quality.md).

**`--pseudocount`** (default 1.0) — added to each transition total. Only used
with `--errfun noqual`.

**`--binned-quals`** — comma-separated anchor quality values for piecewise-linear
interpolation, e.g. `0,10,20,30,40`. Only used with `--errfun binned-qual`.

!!! tip "Observed bins can be narrower than the instrument's documented ones"
    `summary --report` tells you which quality values are *present in your data*,
    which is not necessarily the scheme the vendor documents. NovaSeq is
    documented as binning to **2, 11, 25, 37**, but the Q2 tail is routinely
    removed by trimming or `--trunc-q`, so a trimmed dataset commonly shows only
    **11, 25, 37**.

    Either set works. An anchor with no observations behind it contributes
    nothing — the interpolation segment that would need it is skipped and the
    values below the lowest *observed* anchor are flat-filled, so on one NovaSeq
    soil dataset `2,11,25,37` and `11,25,37` produce **bit-identical** error
    matrices.

    What matters is that the anchors **bracket** the observed range: a quality
    score outside them is an error, not a warning, and the same is true in R.

**`--errfun-cmd`** — command to invoke for `--errfun external`. Whitespace-split
into argv; the trans-input and err-output file paths are appended as the final
two arguments. Both files use R's
`read.table(..., row.names = 1, header = TRUE, check.names = FALSE)` layout.
See [Using an external error model](../walkthroughs/external-error-models.md)
for the contract and the shipped reference scripts, and
[`plot_errors.R`](../scripts/plotting.md#plot_errorsr) for visualising the
result.

### LOESS knobs

**`--loess-surface`** (default `interpolate`) — `interpolate` or `direct`.
`interpolate` fits the local polynomial at kd-tree vertices and blends between
them with cubic Hermite, which is what R's `loess()` does by default and
therefore what `loessErrfun` produces. `direct` evaluates the polynomial at
every query point.

**`--loess-cell`** (default 0.2) — maximum fraction of observations per
kd-tree cell, R's `loess.control(cell=)`. Interpolate surface only.

**`--loess-max-rate`** (default 0.25) / **`--loess-min-rate`** (default 1e-7) —
the clamp applied to fitted off-diagonal rates, matching R's post-fit step
(`errorModels.R:53-56`). Unaffected by the surface.

!!! note "`--loess-preset` is deprecated"
    It was a bundle over the four knobs above, with `default` selecting
    `direct` and `r-dada2` selecting `interpolate`. Since #205 `interpolate`
    is the default, so `--loess-preset r-dada2` is redundant and
    `--loess-preset default` is a confusing name for the non-default surface.
    The flag still works and maps to `--loess-surface`, but warns. Use
    `--loess-surface direct` in place of `--loess-preset default`.

**`--loess-surface`** — `direct` evaluates the local polynomial at every query
point (matches R `loess(surface = "direct")`). `interpolate` builds a 1-D
kd-tree partition, fits at each vertex, and blends with cubic Hermite at
queries (matches R's default `loess()`). Applies to `--errfun loess` and
`--errfun pacbio` only; ignored by `noqual`, `binned-qual` and `external`.

**`--loess-cell`** — maximum fraction of observations allowed per kd-tree cell
before it is subdivided. Only used with `--loess-surface interpolate`. Mirrors
R's `loess.control(cell = ...)`; R's default is 0.2.

**`--loess-max-rate`** / **`--loess-min-rate`** — upper and lower clamps applied
to off-diagonal error rates after fitting. Apply to `loess`, `pacbio`, `noqual`
and `binned-qual`; ignored by `external`. Both presets default to 0.25 and 1e-7,
matching R DADA2. Set to `1.0` / `0.0` respectively to disable.

**`--max-consist`** (default 10) — maximum self-consistency iterations, R's
`MAX_CONSIST`.

## Denoising

These mirror R's `setDadaOpt()` surface and apply to the `dada` runs performed
inside each self-consistency iteration; see the
[parameters page](../parameters.md) for the equivalency table and the
[`dada` page](dada.md#denoising) for what each one does.

One default differs from `dada`: **`--omega-c` defaults to 0 here**, matching R
DADA2's `learnErrors()`, which hard-codes `OMEGA_C = 0` in its internal `dada()`
calls and so overrides the standard `dada()` default of 1e-40. Error inference
should not merge by abundance p-value. Pass `--omega-c 1e-40` to use the
standard value instead.

## Alignment

Identical in meaning to [`dada`'s alignment flags](dada.md#alignment):
`--band` (default 16), `--gap-p` (−8), `--homo-gap-p` (falls back to `--gap-p`),
`--match` (5), `--mismatch` (−4), `--align-backend`.

## Screening

Identical in meaning to [`dada`'s screening flags](dada.md#screening):
`--kdist-cutoff` (default 0.42), `--kmer-size` (default 5), `--no-kmer-screen`.

Note that the cutoff used here and the one used for denoising can be set
independently, and that decoupling them is the safe way to speed up denoising —
see [KDIST cutoff decoupling](../findings/kdist-cutoff-decoupling.md).

## Performance

**`--threads`** (default 1) — threads for parallel sample processing. Does not
affect results.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON instead of pretty-printed.

## Diagnostics

**`--verbose`** — per-iteration progress to stderr.

**`--diag-dir`** — directory for per-iteration cluster diagnostics
(`iter_001.json`, …). Each file holds cluster counts and a birth-type breakdown
per sample for that iteration. Created if absent.

**`--cluster-trace-dir`** — directory for full per-iteration cluster traces
(`cluster_iter_NNN_sample_NNN.json`). Each file describes the full cluster
structure for one sample at one iteration: cluster centers, members with their
hamming distance, λ, expected reads, and abundance p-value, plus the `err`
matrix used for that iteration. See `examples/cluster_trace/` for plotting
scripts.

**`--trace-no-members`** — omit the per-cluster `members` array from trace
files, emitting only centers and birth metadata. Reduces trace size ~10×.

**`--trace-min-abund`** (default 1) — only include trace members at or above
this abundance.

## Experimental

Identical in meaning to [`dada`'s experimental flags](dada.md#experimental):
`--screen-backend`, `--minimizer-k`, `--minimizer-w`, `--screen-audit`,
`--wfa-max-edits`.

## See also

- `errors-from-sample` — fit a model from a single pre-sampled input
- [Binned quality scores](../findings/binned-quality.md)
- [learn-errors --nbases convergence](../findings/learn-errors-nbases-convergence.md)
