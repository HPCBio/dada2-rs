# `reference-eval`

Evaluate a query ASV set against a reference (truth) set by Needleman-Wunsch
alignment.

```bash
dada2-rs reference-eval --asvs run_pooled.json.gz --reference truth.fasta \
  --per-asv per_asv.tsv --per-ref per_ref.tsv --threads 24 -o eval.json
```

Classifies each ASV as TP / near / FP by ends-free global alignment, reports
reference recovery (FN) and precision/recall, and — for a `dada` /
`dada-pooled` JSON input — joins the classification to each ASV's p-value.

This is a pure diagnostic: it only reads existing fields and never changes
inference. Output rows are keyed by ASV sequence, and reference header
annotations pass through verbatim so downstream R or Python analyses can extend
the table.

!!! note "Scope"
    v1 (issue #91) handles a single query set against a single reference. There
    is no TN bucket, no cross-sample logic, and no abundance modeling.

A truth set is only meaningful for the `(mock, region, primer)` combination it
was derived from — an Illumina subregion collapses alleles that a full-length
PacBio run resolves — so recovery fractions are not comparable across platforms
unless both truth sets were derived for their own reads.

## Input

**`--asvs`** (required) — query ASVs: a FASTA, or a `dada` / `dada-pooled` JSON.
The format is auto-detected.

**`--reference`** (required) — reference/truth FASTA, e.g. mock-community
alleles.

**`--non-chimeric`** — the non-chimeric survivor set, as a
`remove-bimera-denovo` / `make-sequence-table` JSON or a FASTA. ASVs whose
sequence is absent are flagged chimeric, and the summary then reports a
chimera × reference-class 2×2: does the FP tail consist of chimeras, and does
chimera removal drop any TP allele?

**`--nraw`** — Bonferroni divisor from a pooled cluster trace, so a `Prior` ASV's
`birth_pval` is reported on the abundance scale as `p_a = birth_pval * nraw`.
Omit to leave `p_a` blank.

## Evaluation

**`--max-diffs`** (default 0) — edit distance (mismatches plus internal indels)
at or below which an ASV counts as a true positive. 0 means an exact,
alignment-internal match.

**`--near-diffs`** (default 3) — upper edit-distance bound for the report-only
"near" bucket, i.e. `max_diffs < edit <= near_diffs`. Must be at least
`--max-diffs`.

## Alignment

**`--match`** (5), **`--mismatch`** (−4), **`--gap-p`** (−8) — NW scoring.

**`--band`** (default 32) — band radius, deliberately over-provisioned. Internal
indels need band coverage; end-length differences are free under ends-free
alignment.

## Screening

**`--kdist-screen`** — optional permissive k-mer prefilter cutoff applied before
NW. Omit to disable, which is recommended for small reference sets: too tight a
screen can drop an ASV's true reference and misclassify it as a false positive.

**`--kmer-size`** (default 5) — k-mer size for that prefilter.

## Performance

**`--threads`** (default 1) — threads for alignment.

## Output

**`--output` / `-o`** — write the summary JSON here instead of stdout.

**`--per-asv`** — write the per-ASV classification table (TSV) here.

**`--per-ref`** — write the per-reference recovery table (TSV) here.

**`--compact`** — minified summary JSON.
