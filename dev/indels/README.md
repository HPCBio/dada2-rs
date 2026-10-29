# Indel-aware denoising tools (#291)

Development scripts for the experimental indel-aware λ. They are not installed
and not part of the CLI; user-facing pieces live in `scripts/plot_indel_errors.R`,
`scripts/flag_artifacts.py` (the post-denoising flagger) and
`examples/external_indel_fit/`. Truth sets come from `dev/insilico_pcr.py`,
which stays in `dev/` because it serves every mock (#166).

| Script | Purpose |
| --- | --- |
| `indel_spikein.py` | Spike known indel variants (and homopolymer controls) into a real run |
| `indel_spikein_eval.py` | Score `dada` output on spiked samples: recall by size, position, type, reads moved |
| `indel_lambda_dryrun.py` | Offline test of an indel-aware λ from `--indel-events` tables, before it was built |
| `indel_pair_concordance.py` | Separate molecular (PCR) from sequencing indels by read-pair concordance |
| `indel_fit_proto.py` | Held-out comparison of indel-rate fits (pooling, LOESS, GLMs) |
| `indel_gam_proto.R` | The same for mgcv GAMs: hierarchical, surface, monotone, centre vs flanking quality |
| `run_indel_ab.sh` | Full-run A/B under `dada-pooled`: substitution-only λ against `--indel-rates` |
| `indel_ab_compare.py` | Judge each indel-arm-only ASV by how it tracks its nearest default ASV across samples |
| `hp_variant_counts.py` | Homopolymer slippage rate per run per parent, from exact read matches |
| `hp_structure_features.py` | Structure covariates per run: Rfam consensus pairing (`cmalign`), local folding (`RNAplfold`, DNA parameters), inverted repeats |
| `hp_structure_proto.R` | Held-out test of structure and sequence context as slippage covariates |
| `align_model_subs.py` | Where an alignment model moves substitutions: realigns a fixed partition under ends-free linear and global affine scoring and classifies the differences (#309; needs Biopython) |

Each script's header documents its inputs, outputs and usage.
