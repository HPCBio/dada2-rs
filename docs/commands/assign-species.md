# `assign-species`

Fill in the Species column of an [`assign-taxonomy`](assign-taxonomy.md) JSON by
exact match. Mirrors R DADA2's `addSpecies()`.

```bash
dada2-rs assign-species tax.json --ref-fasta silva_species_assignment.fa.gz -o tax.species.json
```

Reads an `assign-taxonomy` JSON, runs exact-match species assignment against
`--ref-fasta`, and writes a JSON file with the same shape. The `Species` level
is appended, or replaced if already present.

Species is only filled when the reference's genus matches the query's assigned
Genus level (where one is present), using R's `matchGenera` rules: exact match,
a `"Genus "` prefix, or the split-genus forms `Genus/…` and `…/Genus`.

The reference FASTA must use the DADA2 species-assignment format, where each
header holds three whitespace-delimited fields — accession, genus, species:

```
>AY123456 Staphylococcus aureus
```

## Input

**`<INPUT>`** — taxonomy JSON from `assign-taxonomy`.

**`--ref-fasta`** (required) — reference FASTA with `>ID genus species` headers.

## Classification

**`--allow-multiple`** (default 1) — maximum distinct species returned per
query; `0` means unlimited. The default of 1 returns only unambiguous
assignments, matching R's `allowMultiple = FALSE`.

**`--try-rc`** — also try the reverse complement of each query.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON.

## Diagnostics

**`--verbose`** — progress to stderr.
