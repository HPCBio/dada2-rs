# `assign-taxonomy`

Assign taxonomy using a Naive Bayes k-mer classifier. Mirrors R's
`assignTaxonomy()`.

```bash
dada2-rs assign-taxonomy asvs.fasta \
  --ref-fasta silva_nr99_train_set.fa.gz --threads 24 -o tax.json
```

The query may be a FASTA file or a
[`make-sequence-table`](make-sequence-table.md) JSON. The reference FASTA must
use DADA2-formatted headers, where the description is a semicolon-separated
taxonomy string:

```
>Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;
```

Output is a JSON object with a `levels` array and an `assignments` array — one
entry per query — containing the sequence, its assigned taxonomy (null where
confidence is below `--min-boot`), and optionally the raw bootstrap counts.

## Input

**`<INPUT>`** — query sequences: FASTA (`.fa` / `.fa.gz` / `.fasta`) or a
sequence-table JSON.

**`--ref-fasta`** (required) — reference FASTA with semicolon-delimited taxonomy
headers.

## Classification

**`--min-boot`** (default 50) — minimum bootstrap confidence, 0–100, to assign a
taxonomic level. Levels below it are reported as null rather than guessed.

**`--try-rc`** — also classify the reverse complement of each query and keep the
better-scoring orientation. Worth setting when read orientation is not
guaranteed upstream.

**`--tax-levels`** (default
`Kingdom,Phylum,Class,Order,Family,Genus,Species`) — comma-separated level
names, applied in order. Must match the depth of the reference's taxonomy
strings.

**`--seed`** — RNG seed for reproducible bootstrap sampling. Set it if you need
byte-identical reruns.

## Performance

**`--threads`** (default 1) — threads for parallel query classification.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--output-bootstraps`** — include raw bootstrap counts in the output.

**`--compact`** — minified JSON.

## Diagnostics

**`--verbose`** — progress to stderr.

## See also

- [`assign-species`](assign-species.md) — exact-match species assignment
- [`tax-to-tsv`](tax-to-tsv.md)
