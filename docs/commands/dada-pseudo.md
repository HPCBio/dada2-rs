# `dada-pseudo`

Denoise multiple samples with pseudo-pooling — R DADA2's
`dada(..., pool = "pseudo")`.

```bash
dada2-rs dada-pseudo derep/*.json.gz \
  --error-model err.json -o dada/ --threads 24 --gzip
```

Two per-sample rounds. Round 1 denoises each sample independently with no
priors. The ASVs from round 1 are pooled into a sequence table and a prior set
is selected using R DADA2's `PSEUDO_PREVALENCE` / `PSEUDO_ABUNDANCE` rule.
Round 2 re-runs each sample with those priors flagged, routed through
`--omega-p`. One `{sample}.json` per sample is written to `--output-dir`.

This recovers most of pooling's sensitivity to rare variants at per-sample cost
and memory.

## Input

**`<INPUT>...`** — FASTQ or derep/sample JSON files, one per sample.

**`--sample-names`** — comma-separated, one per input. Defaults to filename
stems.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

There is no `--prior`: pseudo-pooling derives its own prior set from round 1.

## Error model

`--error-model` (required), `--use-err-in`, `--inherit-err-params` — as in
[`dada`](dada.md#error-model).

**`--reestimate-err-between-rounds`** — re-fit the error model from round 1 and
use it for round 2, instead of using the supplied model for both rounds.

This exists because R implements pseudo-pooling as two turns of the
self-consistency loop rather than two `dada()` calls, and re-fits
`err <- errorEstimationFunction(trans)` at the end of every turn with no
`selfConsist` guard (`dada.R:371-378`). A hand-rolled two-round run in R
therefore does not reproduce `pool = "pseudo"` exactly, and this flag is how
dada2-rs can emulate R's behaviour. Treat it as the leading explanation for
round-2 differences against R, not a confirmed one.

## Pseudo-pooling

**`--pseudo-prevalence`** (default 2) — number of samples an ASV must appear in
to become a round-2 prior. R's `PSEUDO_PREVALENCE`; equivalent to
`seq-table-to-fasta --prevalence`.

**`--pseudo-min-abundance`** — total abundance across samples for an ASV to
become a prior. R's `PSEUDO_ABUNDANCE`, whose default is `Inf` (i.e. off).

The two rules are OR'd, matching R:
`colSums(st>0) >= PSEUDO_PREVALENCE | colSums(st) >= PSEUDO_ABUNDANCE`.

**`--priors-out`** — write the selected round-2 priors to this FASTA. Useful for
checking what round 1 actually promoted.

## Denoising, Alignment, Screening

Identical in meaning to `dada`: see [Denoising](dada.md#denoising),
[Alignment](dada.md#alignment) and [Screening](dada.md#screening).

## Performance

**`--threads`** (default 1) — threads for dereplication and DADA2 comparisons.

**`--sample-jobs`** — samples to denoise concurrently, each on its own
`threads / sample-jobs` sub-pool. Defaults to `round(threads / 4)`.

**`--cache-samples`** — keep every sample's uniques in memory across both
rounds. **Off by default, and it should normally stay off.** By default
`dada-pseudo` *streams*: each sample is dropped after round 1 and re-read
(re-dereplicated) in round 2, bounding peak memory to `--sample-jobs` samples in
flight rather than all samples at once. Streaming is both faster and lighter on
large runs — re-dereplication turns out to be cheaper than carrying the cache —
so the cache is pure overhead except in narrow cases.

## Output

**`--output-dir` / `-o`** (required), **`--compact`**, **`--gzip`** — as in
[`dada-pooled`](dada-pooled.md#output).

## Diagnostics

**`--verbose`**, **`--failed-uniques`** — as in [`dada`](dada.md#diagnostics).

## Experimental

Identical in meaning to [`dada`'s experimental flags](dada.md#experimental).

## See also

- [Pseudo-pooled denoising algorithm](../algorithm-dada-pseudo.md)
- [Pseudo-pooling priors vs error model](../findings/pseudo-pooling-priors-vs-error-model.md)
