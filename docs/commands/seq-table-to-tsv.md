# `seq-table-to-tsv`

Convert a sequence table JSON to a tab-delimited count table.

```bash
dada2-rs seq-table-to-tsv seqtab.nochim.json -o seqtab.tsv
```

Reads JSON from [`make-sequence-table`](make-sequence-table.md) or
[`remove-bimera-denovo`](remove-bimera-denovo.md) and writes a TSV with sequence
IDs as rows and sample names as columns.

## Input

**`<INPUT>`** — sequence table JSON.

## Filtering

Both rules mirror R DADA2's pseudo-pooling prior selection and are OR'd:
`colSums(st>0) >= PSEUDO_PREVALENCE | colSums(st) >= PSEUDO_ABUNDANCE`.

**`--prevalence`** — keep only sequences present in at least this many samples
(R's `PSEUDO_PREVALENCE`). Omit to disable.

**`--min-abundance`** — keep only sequences whose total abundance is at least
this value (R's `PSEUDO_ABUNDANCE`). Omit to disable.

## Output

**`--output` / `-o`** — write the TSV here instead of stdout.
