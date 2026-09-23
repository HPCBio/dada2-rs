# Pipeline helpers

Three scripts that are not plots: one carries an R error model in, one rolls
per-sample counts up into a table, and one is a template for getting the final
tables into Bioconductor.

## `learnerrors_to_dada2rs.R`

Convert a model R DADA2 has **already fitted** into a dada2-rs error model, so
inference runs here on R's numbers. This is the closest route to R's output,
because nothing about the model is re-derived.

```bash
# In R:  saveRDS(errF, "errF.rds")   # the learnErrors() result

Rscript scripts/learnerrors_to_dada2rs.R errF.rds errF.json
dada2-rs dada sample.derep.json.gz --error-model errF.json -o sample.dada.json
```

The `.rds` may be either the list `learnErrors()` returns — its `$err_out` is
used — or a bare 16-row error-rate matrix, e.g. from
`saveRDS(getErrors(errF), ...)`. Row order must be `A2A, A2C, A2G, A2T, C2A, …,
T2T`.

The output sets `err_in` and `err_out` to the same matrix, so `--use-err-in`
makes no difference downstream.

The other direction — having dada2-rs call an R or Python *fitting function* per
self-consistency iteration — is `--errfun external`; see
[Using an external error model](../walkthroughs/external-error-models.md), which
covers both routes and when to pick which.

## `track_reads.py`

The SOP's "track reads through the pipeline" table: one row per sample, one
column per step, so a sample that loses its reads somewhere is visible at a
glance. Standard library only.

```bash
python3 scripts/track_reads.py \
    -f filtered/*.json \
    -d seqtab_R1.json seqtab_R2.json \
    -m merged.json \
    -s seqtab.nochim.json \
    -o track.tsv
```

| Flag | Input |
|---|---|
| `-f`, `--filter-and-trim` | one or more `filter-and-trim` JSONs, one per sample |
| `-d`, `--dada` | one (single-end) or two (paired: R1 R2) `make-sequence-table` JSONs built from the per-sample `dada` outputs |
| `-m`, `--merge-pairs` | `merge-pairs` JSON; omit for single-end |
| `-l`, `--length-filtered` | optional sequence table filtered by ASV length |
| `-s`, `--seqtab` | final sequence table, after `remove-bimera-denovo` |
| `-o`, `--output` | write TSV here instead of stdout |

Every input is optional, so the table can be built for whichever steps a run
actually has. Samples missing from a step get `0` in that column rather than
being dropped.

## `examples/import_mia.R`

A **template, not a tool** — it lives in `examples/` because it hardcodes its
input filenames and is meant to be copied and edited, not run with arguments.

It loads the two TSV exports into a `TreeSummarizedExperiment` for the
Bioconductor [`mia`](https://bioconductor.org/packages/mia/) ecosystem:

```bash
dada2-rs seq-table-to-tsv seqtab.nochim.json -o seqtab.tsv
dada2-rs tax-to-tsv taxonomy.json -o taxonomy.tsv
Rscript examples/import_mia.R
```

It intersects the two tables on sequence ID and fails loudly if they share
none — which is the failure you want, since a silent empty intersection
produces an object that looks valid and holds nothing. Sample metadata goes in
as `colData`; the template marks the line.
