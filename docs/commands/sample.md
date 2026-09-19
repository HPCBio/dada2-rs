# `sample`

Dereplicate and subsample FASTQ files, writing one JSON file per sample.

```bash
dada2-rs sample reads/*.fastq.gz -o sampled/ --nbases 1e8 --gzip
```

Processes input FASTQ files in order (or shuffled with `--randomize`),
dereplicating each and writing a JSON file to `--output-dir`. Processing stops
once the cumulative base count reaches `--nbases`.

Each output file uses the same format as [`derep`](derep.md) and can be passed
directly to [`errors-from-sample`](errors-from-sample.md). Pre-sampling is the
way to spread error learning across many samples rather than letting one deep
sample dominate — see the caveat on
[`learn-errors`](learn-errors.md).

## Input

**`<INPUT>...`** — FASTQ files (`.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`).

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Filtering

**`--nbases`** (default 1e8) — stop after accumulating at least this many total
bases across input files. Files are taken whole.

**`--randomize`** — process input files in random order instead of the supplied
order. This shuffles sample order only; reads are never subsampled within a
file.

**`--seed`** — RNG seed for reproducible `--randomize`.

## Performance

**`--threads`** (default 1) — threads for dereplication.

## Output

**`--output-dir` / `-o`** — directory for per-sample JSON files, created if
absent.

**`--pretty`** — pretty-print the JSON; the default is compact (~34% smaller).

**`--gzip`** — write `{sample}.json.gz`, read back transparently downstream.

## Diagnostics

**`--verbose`** — progress to stderr.
