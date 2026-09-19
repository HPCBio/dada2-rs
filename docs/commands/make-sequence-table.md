# `make-sequence-table`

Build a sample-by-sequence feature table. Mirrors R's `makeSequenceTable()`.

```bash
dada2-rs make-sequence-table dada/*.json -o seqtab.json
```

Reads JSON produced by [`dada`](dada.md) (one file per sample) or
[`merge-pairs`](merge-pairs.md) (one file containing multiple samples) and
assembles a flat count matrix, samples × sequences.

## Input

**`<INPUT>...`** — `dada` or `merge-pairs` JSON files.

**`--sample-names`** — sample name for each input file. Only applies to
single-sample `dada` files; `merge-pairs` files carry sample names internally.
If provided, the length must match the number of inputs.

## Filtering

**`--min-len`** / **`--max-len`** — discard ASV sequences shorter / longer than
this length, inclusive. The usual reason to set these is removing off-target
amplicons, which show up as a length mode well away from the expected amplicon
size.

## Output

**`--order-by`** (default `abundance`) — column order: `abundance` (decreasing
total), `nsamples` (decreasing number of samples present in), or `none`
(first-seen order).

**`--hash`** (default `md5`) — hash algorithm for sequence identifiers; `md5` or
`sha1`.

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON.

## See also

- [`remove-bimera-denovo`](remove-bimera-denovo.md) — the usual next step
- [`seq-table-to-tsv`](seq-table-to-tsv.md)
