# `dada-pooled`

Denoise multiple samples with full pooling — R DADA2's `dada(..., pool = TRUE)`.

```bash
dada2-rs dada-pooled derep/*.json.gz \
  --error-model err.json -o dada/ --threads 24 --gzip
```

Per-sample uniques are merged into one combined table (abundances summed,
qualities abundance-weighted-averaged), DADA2 runs **once** on the merged table,
and one JSON file per sample is written to the output directory containing only
the ASVs present in that sample.

Pooling is the most sensitive mode for rare variants, and the most expensive.
Peak memory is driven by the merged unique table rather than any single sample.

## Input

**`<INPUT>...`** — FASTQ or derep/sample JSON files, one per sample.

**`--sample-names`** — comma-separated, one per input. Defaults to filename
stems.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

**`--prior`** — FASTA of prior sequences; exact matches against the merged table
are exempt from the abundance p-value filter and split via `--omega-p`.

## Error model

`--error-model` (required), `--use-err-in`, `--inherit-err-params` — identical
in meaning to [`dada`'s error-model flags](dada.md#error-model).

## Denoising, Alignment, Screening

Identical in meaning to `dada`: see [Denoising](dada.md#denoising),
[Alignment](dada.md#alignment) and [Screening](dada.md#screening).

`--kmer-size` deserves particular attention here. Pooling produces a large,
diverse unique table, which is exactly the regime where the k-mer screen runs on
everything — so on PacBio HiFi, leaving it at 5 is costly. See
[K-mer screen size](../findings/kmer-size-screening.md).

## Performance

**`--threads`** (default 1) — threads for dereplication and the DADA2 comparison
map. Pooled runs scale to high thread counts, but the shuffle phase has a serial
floor that grows as a share of runtime as threads rise; see
[Thread scaling & memory placement](../findings/thread-scaling-and-placement.md).

## Output

**`--output-dir` / `-o`** (required) — directory for per-sample
`{sample}.json`, created if absent.

**`--compact`** — minified JSON.

**`--gzip`** — write `{sample}.json.gz`.

## Diagnostics

**`--verbose`** — progress to stderr.

**`--failed-uniques`** — TSV of uniques that failed to denoise. Because pooled
denoising runs once on the merged unique table, "failed" is a *global* property
(`map == null` on the merged index); for each failed merged unique a row is
emitted per sample it appears in, carrying that sample's read count.

**`--pooled-record`** — write a self-contained pooled record (merged uniques
with pooled abundance, the global map, and global ASVs) for
`kdist-calibrate --from-dada-pooled`. Off by default. Give a path **outside**
`--output-dir` so it does not join the per-sample `*.json.gz` glob; gzip follows
the path's `.gz` extension.

**`--cluster-trace`**, **`--trace-no-members`**, **`--trace-min-abund`** — as in
[`dada`](dada.md#diagnostics).

## Experimental

Identical in meaning to [`dada`'s experimental flags](dada.md#experimental).
Note that on pooled PacBio the minimizer index choice reversed relative to
Illumina, which is why the backend now probes rather than predicts; see
[Minimizers as the screen](../findings/minimizer-screening.md).

## See also

- [Pooled denoising algorithm](../algorithm-dada-pooled.md)
- [`dada-pseudo`](dada-pseudo.md) — cheaper, nearly as sensitive
