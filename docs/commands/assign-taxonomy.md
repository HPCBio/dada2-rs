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

**`--seed`** — RNG seed for the bootstrap. Has a fixed default, so runs are
reproducible without setting it.

!!! tip "Reproducibility, and how to measure the bootstrap's spread"
    Each sequence's bootstrap stream is derived from the **sequence itself**, so
    its classification depends only on the sequence, the reference and the seed.
    Reruns, a different `--threads`, a reordered input and a different set of
    companion sequences all give the same answer for the same sequence.

    There is no unseeded mode, because sampling from entropy only makes a run
    unreproducible. **To see how sensitive a call is, vary the seed** — that
    gives the same spread and can be repeated:

    ```bash
    for s in 1 2 3 4 5; do
      dada2-rs assign-taxonomy asvs.fa --ref-fasta silva.fa.gz --seed "$s" \
        -o "tax.$s.json"
    done
    ```

    Calls that move across seeds are the ones sitting near `--min-boot`. On one
    3,994-query test that was ~7.6% of assignments, all at Family and Genus.
    That spread is inherent to the bootstrap, not a defect — but it is worth
    knowing which of your calls live in it.

    R DADA2's `assignTaxonomy` cannot offer the same guarantee: its RNG advances
    across sequences C-side and `set.seed()` does not reach it, so its output
    depends on input order ([benjjneb/dada2#1115](https://github.com/benjjneb/dada2/issues/1115)).
    Comparisons against R therefore have to be statistical rather than
    exact-match.

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
