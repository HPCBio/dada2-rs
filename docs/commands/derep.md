# `derep`

Dereplicate a FASTQ file into unique sequences with read counts.

```bash
dada2-rs derep reads.fastq.gz -o sample.derep.json.gz
```

Produces the equivalent of R DADA2's `derep` class: a set of unique sequences
with read counts, per-unique integer Phred quality sums (`qual_sum`; mean =
sum / count), and a read-to-unique mapping.

Pre-dereplicating is worth doing whenever you will run `dada` more than once on
the same data — it avoids re-reading the FASTQ on every parameter change.

## Input

**`<INPUT>`** — a single FASTQ file, uncompressed or gzipped.

**`--sample-name`** — sample identifier written to the output JSON's `sample`
field. Downstream subcommands (`dada`, `dada-pooled`) pick this up as their
default when their own sample-name flag is omitted. Defaults to the filename
stem.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Performance

**`--threads`** (default 1) — threads for parallel processing.

## Output

**`--output` / `-o`** — write JSON here instead of stdout. A path ending in
`.gz` is gzip-compressed, and is read back transparently by every subcommand
that accepts derep JSON.

**`--show-map`** — include the per-read mapping (read index → unique index).

**`--pretty`** — pretty-print the JSON. The default is compact (minified),
which is ~34% smaller on disk.

## Diagnostics

**`--verbose`** — progress to stderr.

## See also

- [`sample`](sample.md) — dereplicate *and* subsample, one JSON per sample
- [`dada`](dada.md)
