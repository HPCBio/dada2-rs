# `seq-table-to-fasta`

Convert a [`make-sequence-table`](make-sequence-table.md) JSON to FASTA, one
record per sequence, using the sequence ID as the header.

```bash
dada2-rs seq-table-to-fasta seqtab.json -o asvs.fasta
```

To extract pseudo-pooling priors by hand — mirroring R DADA2's `pool="pseudo"`
selection rule — pass `--prevalence` and/or `--min-abundance`. Note that
[`dada-pseudo`](dada-pseudo.md) applies this rule internally, so you only need
this for inspection or for a manual two-round workflow.

## Input

**`<INPUT>`** — sequence table JSON.

## Filtering

**`--prevalence`** — keep only sequences present in at least this many samples
(R's `PSEUDO_PREVALENCE`, default 2 in R). Omit to disable.

**`--min-abundance`** — keep only sequences whose total abundance across samples
is at least this value (R's `PSEUDO_ABUNDANCE`). Omit to disable, which is
equivalent to R's default of `Inf`.

The two rules are OR'd.

## Output

**`--output` / `-o`** — write the FASTA here instead of stdout.
