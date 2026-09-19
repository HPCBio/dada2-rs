# `filter-and-trim`

Filter and trim a single sample's FASTQ reads. Mirrors R's `filterAndTrim()`
for one sample.

```bash
# single-end
dada2-rs filter-and-trim --fwd raw.fastq.gz --filt filt.fastq.gz --trunc-len 240 --max-ee 2

# paired-end
dada2-rs filter-and-trim \
  --fwd raw_R1.fastq.gz --filt filt_R1.fastq.gz \
  --rev raw_R2.fastq.gz --filt-rev filt_R2.fastq.gz \
  --trunc-len 240 200 --max-ee 2 2
```

Several parameters accept **paired values**: give one value to apply it to both
directions, or two space-separated values (forward first, then reverse). These
are `--trunc-q`, `--trunc-len`, `--trim-left`, `--trim-right`, `--max-len`,
`--min-len`, `--max-ee` and `--rm-lowcomplex`.

Use [`summary --expected-error`](summary.md#metrics) first to choose
`--max-ee` and `--trunc-len` from the data rather than by habit.

## Input

**`--fwd`** (required) — forward (R1) input FASTQ.

**`--rev`** — reverse (R2) input FASTQ. Supplying it enables paired-end mode,
and `--filt-rev` then becomes required.

**`--sample-name`** — sample identifier for the output JSON. Defaults to the
`--fwd` filename stem.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Trimming

Applied in this order, before the filtering rules below.

**`--trunc-q`** (default 2) — truncate each read at the first Phred score at or
below this value.

**`--trunc-len`** (default 0 = off) — truncate reads to this many bases,
discarding reads that are shorter.

**`--trim-left`** / **`--trim-right`** (default 0) — remove this many bases from
the 5′ / 3′ end.

## Filtering

**`--max-len`** (default 0 = no limit) — discard reads longer than this, applied
*before* trimming.

**`--min-len`** (default 20) — discard reads shorter than this after all
trimming.

**`--max-n`** (default 0) — discard reads with more than this many `N` bases.
The default of 0 discards any read containing an `N`, matching R; DADA2's core
cannot handle ambiguous bases.

**`--min-q`** (default 0 = off) — discard reads containing any Phred score below
this value.

**`--max-ee`** — discard reads whose total expected errors, `Σ 10^(-Q/10)`,
exceed this. Omit for no EE filtering. This is usually the most useful quality
filter, and the most defensible one to tune.

**`--phix-genome`** — FASTA of the phiX genome; matching reads are removed. Omit
to skip phiX filtering.

**`--rm-lowcomplex`** (default 0 = off) — discard reads whose 2-mer Shannon
richness falls below this value.

## Performance

**`--threads`** (default 1) — threads for bgzf output compression. Values above
1 switch the output to bgzf (blocked gzip), which is valid gzip but seekable.

## Output

**`--filt`** (required) — forward (R1) output FASTQ.

**`--filt-rev`** — reverse (R2) output FASTQ; required when `--rev` is given.

**`--compress`** (default true) — gzip-compress output files.

**`--output` / `-o`** — write the JSON summary here instead of stdout.

**`--compact`** — minified JSON.

## Diagnostics

**`--verbose`** — progress to stderr.

## See also

- [`remove-primers`](remove-primers.md) — primer removal, with the same
  filtering and trimming flags available afterwards
- [`summary`](summary.md) — choose thresholds from the data
