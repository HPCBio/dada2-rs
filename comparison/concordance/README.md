# Concordance guardrail — dada2-rs vs R DADA2

A CI sanity check that short `dada2-rs` runs stay consistent with R DADA2 on
small fixtures. It compares the **post-chimera ASVs and their counts** against a
**static R reference CSV** committed to the repo — CI never runs R.

This is a *guardrail*, not an exact-match test: different alignment kernels,
tie-breaks, and floating point mean the two will never be bit-identical. The
comparison is threshold-based and meant to catch **regressions or divergence**.

## Pieces

| File | Role |
|---|---|
| `compare_to_reference.py` | Diffs a dada2-rs seqtab JSON vs the R reference CSV — recall, precision, log-count correlation. Pure stdlib. Warn-only by default; `--gate` to fail. |
| `run_illumina.sh` | dada2-rs paired-end pipeline → `seqtab.nochim.json`. |
| `run_pacbio.sh` | dada2-rs single-end PacBio pipeline → `seqtab.nochim.json` (k=5, to match R). |
| `write_reference.R` | Run **once** in R to produce the reference CSV. Defines the CSV schema. |
| `reference/` | Committed reference CSVs (you generate these). |
| `data/pacbio/` | Committed subsampled PacBio fixture (you add this). |
| `../../.github/workflows/concordance.yml` | The workflow: PRs + branch pushes + manual on main. |

Both toy fixtures are committed under `comparison/concordance/data/` (the repo-root
`data/` is gitignored — it holds the large local-only benchmark datasets — so CI
fixtures live here):

- **Illumina** (`data/illumina/`, ~0.9 MB): 2 samples of the MiSeq SOP tutorial
  data — the exact dataset used in the standard DADA2 tutorial (full set:
  <https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip>). The 2-sample
  subset here ships with the dada2 R package as its built-in test data.
- **PacBio** (`data/pacbio/`, ~2 MB): 2 samples × 1500 reads subsampled from a
  downsampled Sequel IIe 16S set (see `data/pacbio/README.md`). Denoises to ~55
  ASVs in seconds.

Both run in CI today, and only their R reference CSVs generated to enable
the comparison.

## Simple reference CSV schema

Long format, one row per (ASV, sample) with count > 0, from `seqtab.nochim`:

```
sequence,sample,count
ACGT...,sam1,142
ACGT...,sam2,98
```

`compare_to_reference.py` upper-cases sequences and compares on **total count per
ASV** (summed across samples), which is robust to per-sample assignment noise.

## What you need to generate (one-time, needs R + dada2)

The workflow runs the dada2-rs pipelines today but **skips the comparison** until
the reference CSVs (and the PacBio fixture) are committed. To enable it:

1. **Illumina reference** — on a machine with R + dada2 (fixture already committed):
   ```bash
   Rscript write_reference.R illumina data/illumina reference/illumina_seqtab_nochim.csv
   ```
   Commit `reference/illumina_seqtab_nochim.csv`.

2. **PacBio reference** — fixture already committed under `data/pacbio/`:
   ```bash
   Rscript write_reference.R pacbio data/pacbio reference/pacbio_seqtab_nochim.csv \
       AGRGTTYGATYMTGGCTCAG RGYTACCTTGTTACGACTT
   ```
   Commit `reference/pacbio_seqtab_nochim.csv`.

> **Parameters must match.** `write_reference.R` and the `run_*.sh` scripts use
> the same truncation/filter/error settings and `pool=FALSE`. If you change one,
> change both and regenerate the reference.

## Tuning and turning on the gatekeeper

The workflow is **warn-only**: it prints metrics (and a job summary) but stays
green. Once a few runs show stable numbers, edit the thresholds in
`concordance.yml` (`MIN_RECALL` / `MIN_PRECISION` / `MIN_COUNT_CORR`) and add
`--gate` to the compare steps to make divergence fail the build.

## Local use

```bash
cargo build --release
comparison/concordance/run_illumina.sh \
    ./target/release/dada2-rs comparison/concordance/data/illumina /tmp/ill 4
comparison/concordance/compare_to_reference.py \
    --rs /tmp/ill/seqtab.nochim.json \
    --reference comparison/concordance/reference/illumina_seqtab_nochim.csv \
    --min-abundance 2
```
