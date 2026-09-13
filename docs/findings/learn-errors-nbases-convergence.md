# `--nbases`: when does the error model stop moving?

**Verdict: the default `--nbases 1e8` is inside the noise floor of *which samples
you happened to draw*.** On a full MiSeq run, raising `--nbases` below ~500 Mbp
moves the fitted model by about as much as simply reshuffling the sample order at
a fixed budget. Only above 500 Mbp does added depth move the model *less* than
sample choice does. The practical reading: **`--nbases` is not a convergence knob
until it is large enough to stop being a sampling knob.**

This confirms and quantifies the aside recorded during the
[KDIST cutoff work](kdist-cutoff-decoupling.md) — that the error model was "not
converged at nbases=1e8" — and connects it to
[issue #68](https://github.com/HPCBio/dada2-rs/issues/68) (cross-sample diversity
in `learn-errors`).

## Why a naive ladder cannot answer this

`learn-errors` accumulates **whole files**. It never subsamples reads *within* a
file, so each input contributes all of its bases until the running total reaches
`--nbases` (this is documented in the flag's own help). A plain ladder therefore
varies **how many samples were read** at the same time as **how much data**, and
the two cannot be separated after the fact.

On the 20-sample MiSeq SOP subset the confound is total: at ~1.7 Mbp per sample,
`--nbases 2e6` is roughly *one* sample and `2e7` is roughly twelve. A "depth"
ladder there is really a sample-count ladder.

So the experiment (`dev/run_nbases_ladder.sh`) runs two arms:

- **Arm A — depth.** Every sample is pre-cut to the *same* read depth, then
  `--nbases` is laddered. A rung now adds depth uniformly instead of adding
  samples.
- **Arm B — composition.** `--nbases` is held *fixed* and only the sample draw
  varies (`--randomize`, five seeds). This is the **noise channel**: whatever
  spread it shows is what a model moves by for reasons that have nothing to do
  with how much data you gave it.

Arm B is the part that makes Arm A interpretable, and it is the piece a
convergence scan usually omits.

## Results

Full MiSeq run, forward reads, `loess` errfun. Distances are max/RMS difference
in `log10(rate)` over off-diagonal entries, aligned by quality column.

### Arm A — depth, sample composition held fixed

| step | max fold | rms Δlog10 |
|---|---|---|
| 20 Mb → 100 Mb | 1.43× | 0.042 |
| 100 Mb → 500 Mb | 1.41× | 0.034 |
| **500 Mb → 1 Gb** | **1.11×** | 0.014 |
| 1 Gb → 5 Gb | 1.00× | 0.000 |

The last rung is free: the run holds ~780 M transitions, so 1 Gb and 5 Gb fit the
same data and produce the identical model. The real top of the ladder is 1 Gb.

### Arm B — sample draw, depth held fixed at `--nbases 1e8`

| step | max fold |
|---|---|
| seed1 → seed2 | 1.25× |
| seed2 → seed3 | 1.26× |
| seed3 → seed4 | 1.30× |
| seed4 → seed5 | 1.41× |

Five models fitted on the *same amount of data* differ by **1.25–1.43×**.

### Putting them together

| regime | depth effect | composition noise | reading |
|---|---|---|---|
| ≤ 100 Mb | 1.43× | 1.25–1.43× | indistinguishable — `--nbases` is a sampling knob |
| 100 → 500 Mb | 1.41× | 1.25–1.43× | still inside the noise |
| **500 Mb → 1 Gb** | **1.11×** | 1.25–1.43× | **depth effect finally drops below composition** |

**The default `1e8` sits 1.38× from the full-data model, while merely changing the
sample draw moves things 1.25–1.43×.** At that budget the two are the same size,
so a model fitted at the default is as much an artefact of *which samples were
read first* as of the data volume.

## What this dictates

- **Do not treat `--nbases 1e8` as converged.** On a run of this size the model is
  still ~1.4× from its full-data limit in its worst entry, and that gap is
  indistinguishable from sampling noise.
- **Report the composition channel whenever `--nbases` is varied.** A `--nbases`
  result without a fixed-budget reshuffle arm cannot distinguish "more data
  helped" from "different samples were read". This is the same discipline as the
  timing control channel in [measuring on a NUMA node](measuring-on-numa.md), and
  it is what turned this scan from ambiguous into decisive.
- **Prefer spreading learning across samples over raising the budget.** The
  `--nbases` help already recommends pre-`sample`ing each file; these numbers say
  why. Equalising per-sample depth (Arm A) is what made the depth signal legible
  at all.
- **How much does this matter downstream?** Less than it looks. The
  [minimizer screen work](minimizer-screening.md) found a **90.9%** relative
  change in `err_out` propagating to only 0.24% read-level L1 and one churned
  ASV. The error model is markedly more robust to perturbation than denoising is,
  so a 1.4× model difference is unlikely to move a table much — but it has not
  been measured directly here, and "unlikely" is not "measured".

## Caveats

- **One run, one platform.** Forward reads of a single MiSeq run. Whether ~500 Mbp
  is the general stabilisation point or specific to this run's diversity is
  untested; a more diverse run should need more.
- **The ladder saturates.** 1 Gb and 5 Gb are the same model because the data ran
  out at ~780 M transitions. The true convergence point could be above the data
  available here.
- **Arm A and Arm B use slightly different populations** — Arm A on equal-depth
  samples, Arm B on the full filtered set — so the comparison between arms is
  between comparable transition counts (~1e8) rather than identical inputs.
- **`nq` differs between rungs** (40 at 20 Mb, 41 above). `dev/compare_error_models.py`
  aligns by quality column and truncates; comparing without that manufactures
  enormous spurious differences that land exactly on the `nq` boundary. That
  mistake produced an apparent 100× effect in an earlier pass of this analysis.

## Reproducing

```bash
# One read orientation only -- a mixed R1/R2 set is rejected, because
# learn-errors would fit ONE model across two different error profiles.
mkdir r1 && ln -s /path/to/raw/*_R1_*.fastq.gz r1/

bash dev/run_nbases_ladder.sh ./target/release/dada2-rs r1 out 32
#   EQUAL_DEPTH=<reads>   per-sample depth for arm A (set near the shallowest sample)
#   LADDER="..."          arm A rungs
#   FIXED_NBASES=1e8      arm B budget
#   SEEDS="1 2 3 4 5"     arm B draws
```
