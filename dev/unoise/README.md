# Single-pass and UNOISE arm tools (#323)

Development scripts for the experimental `--partition single-pass` and
`--birth-rule unoise` (issue #323). They are not installed and not part of the
CLI. The flags are documented in `docs/commands/dada.md`.

| Script | Purpose |
| --- | --- |
| `run_arms.sh` | Run the arm grid on one dataset: learn once, then denoise, merge, remove chimeras and optionally score each arm against a truth set |
| `summarize_arms.py` | One row per arm: ASV counts, set identity and Bray-Curtis against a baseline arm, dada time, and `reference-eval` scores |
| `phase_share.py` | Phase 0: the share of denoising time that the single pass could remove, from `--metrics-json` |

Each script's header documents its inputs, outputs and environment.

## The arms

| Arm | Flags | `_indel` variant |
| --- | --- | --- |
| `divisive` | none (`dada` as released) | `--indel-rates` |
| `sp` | `--partition single-pass` | `--indel-rates` |
| `sp_reassign` | `--partition single-pass --reassign` | `--indel-rates` |
| `unoise` | `--partition single-pass --birth-rule unoise` | none: d counts gaps by definition |

Compare `sp_indel` with `unoise` to test the learned error model against the
fixed rule with indels counted in both (`summarize_arms.py --pair
sp_indel:unoise`). The indel λ counts indel events and the UNOISE rule counts
gap columns, so the two still count indels differently; see the command page.

One `learn-errors` run serves every arm. Learning with `--indel-rates-out`
leaves the substitution model byte-identical, so the arms with and without
indels share it.

`run_arms.sh` checks that each arm's output JSON records the partition, birth
rule, `--reassign`, indel λ and `--min-abund` it was asked for, and stops if
not. Arms that silently ran alike would look like a null result.

## `--min-abund`

UNOISE needs an abundance cutoff (vsearch `--cluster_unoise` defaults to
`--minsize 8`); DADA2 defaults to 1. Run each dataset twice, sharing one
`LEARN_DIR`:

- `MIN_ABUND` unset: each arm at its own default, `unoise` at
  `UNOISE_MIN_ABUND` (8). This compares the methods as each would be used.
- `MIN_ABUND=8`: every arm at 8. This isolates the partition and birth rule
  from the cutoff.

## Running a dataset

The manifest is tab-separated with no header: `sample  fwd  [rev]`. Inputs must
already be trimmed and filtered, for example the `filtered/` directory of a
`dev/concordance/run_illumina.sh` or `run_pacbio.sh` run. Row order selects the
training set when `--nbases` truncates.

```bash
BIN=target/release/dada2-rs

# MiSeq V4
THREADS=24 bash dev/unoise/run_arms.sh $BIN out/miseq miseq.tsv
THREADS=24 MIN_ABUND=8 LEARN_DIR=out/miseq/learn \
  bash dev/unoise/run_arms.sh $BIN out/miseq_ma8 miseq.tsv

# NovaSeq, binned quality scores (bins from `dada2-rs summary --report`)
LEARN_EXTRA="--errfun binned-qual --binned-quals 2,11,25,37" \
  THREADS=24 bash dev/unoise/run_arms.sh $BIN out/novaseq novaseq.tsv

# PacBio HiFi, MSA-1003 (k=7, the PacBio recommendation)
LEARN_EXTRA="--errfun pacbio" EXTRA="--band 32 --kmer-size 7" \
  REFERENCE=msa1003_truth.fa THREADS=24 \
  bash dev/unoise/run_arms.sh $BIN out/msa1003 msa1003.tsv

# ITS2: no fixed-length truncation in the prep (cutadapt only)
THREADS=24 bash dev/unoise/run_arms.sh $BIN out/its2 its2.tsv

python3 dev/unoise/summarize_arms.py out/msa1003 --pair sp_indel:unoise --tsv out/msa1003.tsv
```

A run stopped part-way resumes: finished arms (those with `<arm>/.done`) are
skipped unless `FORCE=1`, and learning is skipped when `LEARN_DIR` already holds
the error models and indel rates.
