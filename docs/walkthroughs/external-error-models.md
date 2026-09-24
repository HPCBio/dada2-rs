# Using an external error model

`dada2-rs` does not have to fit the error model itself. You can hand the fit to
an R or Python script, or skip fitting entirely and use a model R DADA2 has
already learned. Both routes are supported, and both exist because the project
needed them to establish R concordance in the first place — see
[the LOESS findings page](../findings/loess-error-model-correctness.md).

This is a **fidelity** facility, not a correctness one. The built-in error
functions are the recommended path; nothing here produces a *more accurate*
model, it produces *R's* model. Reach for it when agreeing with a particular R
run is the requirement.

Three reasons this comes up:

- **R parity.** You are reproducing a published table, or diffing against an R
  pipeline during a migration, and need output as close to that run as possible.
- **Binned quality data.** The stock `loessErrfun` is known to misbehave on
  NovaSeq-style data, and the R community has proposed several replacements. You
  want to try one without waiting for it to be ported.
- **A new chemistry.** You are prototyping a model and do not want to rebuild
  `dada2-rs` to test each iteration.

## Which route

| You have | Use | What happens |
|---|---|---|
| A model R already fitted (`.rds`) | `scripts/learnerrors_to_dada2rs.R` | Converted once to JSON, passed to `dada --error-model`. Nothing is re-fitted. |
| A fitting function in R or Python | `learn-errors --errfun external` | `dada2-rs` runs your script once per self-consistency iteration. |
| Neither | the built-in `--errfun` values | `loess`, `noqual`, `binned-qual`, `pacbio` |

The first route is the closest to R's output, because the model is R's. The
second re-runs the self-consistency loop under `dada2-rs`, so the *loop* is ours
and only the fit is yours.

## Route 1 — reuse R DADA2's fitted model

```bash
# In R, once
# saveRDS(errF, "errF.rds")   # the learnErrors() result, or getErrors(errF)

Rscript scripts/learnerrors_to_dada2rs.R errF.rds errF.json

dada2-rs dada sample.derep.json.gz --error-model errF.json -o sample.dada.json
```

The `.rds` may be either the list `learnErrors()` returns (its `$err_out` is
used) or a bare 16-row error-rate matrix. Row order must be `A2A, A2C, A2G,
A2T, C2A, …, T2T`. The converter writes both `err_in` and `err_out` to the same
matrix, so `--use-err-in` makes no difference downstream.

## Route 2 — supply the fitting function

```bash
dada2-rs learn-errors *.fastq.gz \
    --errfun external \
    --errfun-cmd "Rscript examples/external_errfun/loess_reference.R" \
    -o err.json
```

`dada2-rs` writes a transition matrix to a temporary TSV, runs
`<your command> <trans-tsv> <err-tsv>`, and reads the rates back — **once per
self-consistency iteration**, so a slow script is paid for repeatedly.

The command string is whitespace-split into argv with no shell interpolation;
wrap it in a shell script if you need quoting.

### Scripts that ship with the repo

In `examples/external_errfun/`:

| File | Needs | What it is |
|---|---|---|
| `loess_reference.R` | base R | Verbatim port of `dada2:::loessErrfun`. The bit-parity reference. |
| `loess_reference_direct.R` | base R | Same, with `surface = "direct"` — what our native LOESS matches to machine precision. |
| `loess_modified.R` | R + dplyr, magrittr | `span = 2` with `log10(tot)` weights (Salazar/Ruscheweyh, dada2#938) plus a Q40 floor for monotonicity (Holland-Moritz, dada2#791). Aimed at binned data. |
| `pacbio_reference.R` | base R | The PacBio errfun as an external reference. |
| `noqual.py` | Python 3 stdlib | `noqualErrfun` equivalent; demonstrates the contract without R. |

### Writing your own

Your script must read a 16 × `nq` matrix of non-negative integers from
`argv[1]`, write a 16 × `nq` matrix of rates in `[0, 1]` to `argv[2]`, and use
the same layout for both — a leading tab, then integer column labels `0 … nq-1`,
then 16 rows labelled `A2A … T2T`. That is exactly R's
`read.table(..., row.names = 1, header = TRUE, check.names = FALSE)` and
`write.table(err, path, sep = "\t", quote = FALSE, col.names = NA)`.

`dada2-rs` validates the dimensions and the value range, then enforces the
diagonal self-transition probabilities itself — so you may either reconstruct the
diagonal or leave those rows alone.

## Evaluating error models on binned data

This is the use case the external route was most useful for. Binned-quality
chemistries collapse Phred scores to a handful of levels, stock `loessErrfun`
has no guard for that, and several replacements have been proposed in
[dada2#1307](https://github.com/benjjneb/dada2/issues/1307) without a controlled
comparison between them.

Run `summary --report` first to find out whether your data is binned at all.
Then the arms worth comparing on **one fixed input** are the built-in
`binned-qual`, `loess` and `pacbio`, against `loess_modified.R` via the external
route. Holding the reads fixed is what isolates the errfun term — the published
comparisons mostly vary binning *and* errfun together, which is why they cannot
settle the question. That sweep is [#98](https://github.com/HPCBio/dada2-rs/issues/98).

What the project has measured so far is on the
[binned quality findings pages](../findings/binned-quality.md), and the short
version is that the answer is dataset-dependent and not yet predictable.

!!! warning "Score on ASVs, not on the curves"
    Every arm here changes the error rates visibly. That is not the question.
    Compare the final chimera-filtered tables with `dev/compare_asvs.py` using
    set identity and abundance — a curve that looks better is not evidence, and
    `n_asv` alone is a trap.

## Two things that will bite you

**Pin the input set.** An error model is sensitive to *which sequences it was
fitted from*, to roughly the same degree it is sensitive to how much data it
saw — see [`--nbases` and error-model
convergence](../findings/learn-errors-nbases-convergence.md). If you compare two
errfuns on different inputs, or let `--nbases` truncate a differently-ordered
file list, you are measuring two things at once. Use the same derep JSONs for
every arm.

**Provenance is recorded but not verified.** The `--errfun-cmd` string is written
into the output JSON's `params` block alongside `errfun: "external"`, so a run is
self-describing. The *contents* of the script are not hashed — editing it in
place leaves no trace in the output.

## See also

- [`learn-errors`](../commands/learn-errors.md) — all the errfun flags
- [`dada`](../commands/dada.md) — `--error-model`
- [The LOESS error model](../findings/loess-error-model-correctness.md) — why
  these routes exist and what they showed
- [Binned quality scores](../findings/binned-quality.md) — what has been measured
