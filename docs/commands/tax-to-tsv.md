# `tax-to-tsv`

Convert an [`assign-taxonomy`](assign-taxonomy.md) or
[`assign-species`](assign-species.md) JSON to a TSV table.

```bash
dada2-rs tax-to-tsv tax.json -o tax.tsv
```

Emits one row per assignment, with the sequence ID first, followed by one column
per taxonomic level in the order they appear in the input JSON. Unassigned
levels are written as `NA`, matching R DADA2 output.

## Input

**`<INPUT>`** — JSON from `assign-taxonomy` or `assign-species`.

## Output

**`--na-string`** (default `NA`) — string written for unassigned (null) levels.

**`--output` / `-o`** — write the TSV here instead of stdout.
