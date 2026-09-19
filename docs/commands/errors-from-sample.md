# `errors-from-sample`

Learn an error model from pre-computed sample JSON files.

```bash
dada2-rs sample reads/*.fastq.gz -o sampled/ --gzip
dada2-rs errors-from-sample sampled/*.json.gz --threads 24 -o err.json
```

Reads JSON files produced by [`sample`](sample.md) (or [`derep`](derep.md)) and
iteratively runs the DADA2 algorithm, re-fitting the chosen error model until
self-consistency. Output is the same three flat 16 × nq matrices as
[`learn-errors`](learn-errors.md): `trans`, `err_in`, `err_out`.

**This is the two-step form of `learn-errors`.** The difference matters: `sample`
controls *what* is accumulated, so pre-sampling each file and then fitting here
spreads learning across samples rather than letting one deep sample fill the
whole budget. If you are reaching for `learn-errors --nbases` and worrying about
which samples it actually saw, this pair is the answer.

## Input

**`<INPUT>...`** — sample or derep JSON files (`.json` / `.json.gz`).

There is no `--nbases` / `--randomize` / `--phred-offset` here; those belong to
the [`sample`](sample.md) step that produced the input.

## Error model, Denoising, Alignment, Screening, Diagnostics, Experimental

Identical in meaning to [`learn-errors`](learn-errors.md). In particular
`--omega-c` defaults to 0 here too, matching R DADA2's `learnErrors()`.

## Performance

**`--threads`** (default 1) — threads for parallel sample processing.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON.

## See also

- [`learn-errors`](learn-errors.md) — the one-step form
- [`sample`](sample.md)
- [learn-errors --nbases convergence](../findings/learn-errors-nbases-convergence.md)
