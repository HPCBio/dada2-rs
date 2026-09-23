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

`loess` vs R DADA2: the native Rust LOESS is algorithmically equivalent to R's
`loess(surface = "direct")` — bit-exact to machine precision on real data,
validated against `examples/external_errfun/loess_reference_direct.R`. R
DADA2's `loessErrfun`, however, calls `loess(...)` with R's default
`surface = "interpolate"`, which fits the local polynomial at kd-tree vertices
and interpolates between them. The two surfaces disagree by ~1e-3 absolute /
~4% relative at low-Q edges; the downstream impact on inference is minimal
(~1 read per sample on a 362-sample benchmark).

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
once and pass it to `dada` with `--error-model`:

```bash
Rscript scripts/learnerrors_to_dada2rs.R errF.rds errF.json
dada2-rs dada sample.derep.json.gz --error-model errF.json -o sample.dada.json
```

The input `.rds` may be the list `learnErrors()` returns or a bare 16-row error
matrix. This is the closest route to R's output, since nothing about the model
is re-derived.

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

**`--errfun-cmd`** — command to invoke for `--errfun external`. Whitespace-split
into argv; the trans-input and err-output file paths are appended as the final
two arguments. Both files use R's
`read.table(..., row.names = 1, header = TRUE, check.names = FALSE)` layout.
See [Using an external error model](../walkthroughs/external-error-models.md)
for the contract and the shipped reference scripts.

### LOESS knobs

**`--loess-preset`** (default `default`) — resolves a bundle of related knobs
(`--loess-surface`, `--loess-cell`, `--loess-max-rate`, `--loess-min-rate`).
Any of those flags passed explicitly overrides the preset for that knob.

| Preset | Surface | Cell | Max rate | Min rate |
|---|---|---|---|---|
| `default` | `direct` | — | 0.25 | 1e-7 |
| `r-dada2` | `interpolate` | 0.2 | 0.25 | 1e-7 |

`r-dada2` mirrors R DADA2's `loessErrfun`: R's default `loess()` surface plus
the same `[1e-7, 0.25]` clamp R applies after the fit
(`errorModels.R:53-56`). Both presets clamp to the same range; they differ only
in the fitting surface.

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
