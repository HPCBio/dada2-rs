# Choosing the error-model training set

**There is no measured best practice.** How many bases to train on, whether
`--randomize` helps, and which samples to train on are all settled by advice
from the R DADA2 maintainers, not by measurement. A survey of the benjjneb/dada2
issue tracker, tutorials and workflows (October 2026, about 100 threads) found:

- a few small comparisons run by users, none of them against a truth set;
- no peer-reviewed study of training-set choice.

This page records that advice and labels each item **advice** or **evidence**,
so later results can replace it. That measurement is
[issue #105](https://github.com/HPCBio/dada2-rs/issues/105). Until it reports,
treat the defaults below as provisional.

## What is measured

| Question | Evidence | Source |
|---|---|---|
| How many bases | Small effect above 1e8 on the error plots; ASV sets 98.5% shared at 1e8 vs 5e8 | User comparisons, no truth set |
| How many bases (dada2-rs) | Below ~500 Mbp, changing the budget moves the model about as much as changing which samples are read | [nbases convergence](findings/learn-errors-nbases-convergence.md), one MiSeq run, likely heterogeneous |
| `--randomize` | Five shuffled sample orders at 1e8 differ by 1.25–1.43× | Same finding (Arm B); nothing in R found |
| Low-diversity or same-environment samples | +18% reads kept for sediment samples with a sediment-only model | User report, read retention only |
| Mock or control samples only | A better-looking error plot on PacBio full-length 16S | User report, plots only |
| One model per run vs pooled runs | None found | — |

## How training samples are used

R and dada2-rs behave the same way:

- **Samples are taken whole.** Files are read in order, each adds all of its
  bases, and reading stops once the total passes `--nbases`. Reads are never
  subsampled within a file, so one deep sample can use the whole budget (see the
  warning on [learn-errors](commands/learn-errors.md)).
- **`--randomize` only shuffles the file order.** It does not spread reads across
  samples.
- **Training samples are denoised one at a time.** Each sample is denoised
  separately and the transition counts are summed before the error function is
  fitted. Ben Callahan
  ([dada2#2054](https://github.com/benjjneb/dada2/issues/2054)): "The samples
  themselves are not pooled. The error rates learned from each sample
  (individually) are averaged across samples."
- **The fit starts high.** The first iteration treats the most abundant
  sequence as correct and everything else as error, so rates approach their
  final value from above.

## Maintainer advice

### How many bases (advice)

- 1e8 was chosen to match the old default of one million 100 nt reads; "it
  really is the number of bases that is critical"
  ([dada2#443](https://github.com/benjjneb/dada2/issues/443)). The
  [big-data workflow](https://benjjneb.github.io/dada2/bigdata.html) calls 1e8
  "more than adequate".
- More data gives a negligible gain and costs only time
  ([dada2#672](https://github.com/benjjneb/dada2/issues/672),
  [dada2#938](https://github.com/benjjneb/dada2/issues/938),
  [dada2#1740](https://github.com/benjjneb/dada2/issues/1740)). If the error plot
  fits poorly on a large dataset, try raising `nbases`
  ([tutorial 1.8](https://benjjneb.github.io/dada2/tutorial_1_8.html)).
- Downsample a single very deep sample before training. With several samples,
  the R maintainer prefers keeping samples whole and going over the budget
  ([dada2#2054](https://github.com/benjjneb/dada2/issues/2054),
  [dada2#240](https://github.com/benjjneb/dada2/issues/240)). This differs from
  the suggestion on [learn-errors](commands/learn-errors.md) to pre-`sample` each
  file, and neither has been measured.
- No platform-specific budget advice was found. The PacBio workflow uses the
  default 1e8.

### Which samples (advice)

- Train on lower-diversity samples. In very diverse samples (soil, sediment) the
  fit converges above the true rates, which is conservative but loses rare
  variants. Add lower-diversity samples prepared the same way
  ([dada2#851](https://github.com/benjjneb/dada2/issues/851),
  [dada2#2151](https://github.com/benjjneb/dada2/issues/2151),
  [dada2#2208](https://github.com/benjjneb/dada2/issues/2208)).
- But "probably a little diversity is better", so the nucleotides at each read
  position vary ([dada2#2080](https://github.com/benjjneb/dada2/issues/2080)).
  Choosing samples for low diversity can also make the training set
  unrepresentative of the run, the concern in #105.
- Training only on mock or high-replication controls is "even preferred",
  provided they went through the same PCR and sequencing
  ([dada2#2032](https://github.com/benjjneb/dada2/issues/2032)). For
  low-duplication PacBio data, train on samples with more than 50% duplication
  ([dada2#2135](https://github.com/benjjneb/dada2/issues/2135)).
- Choose whole samples at random, never single reads from many samples: pooling
  reads from different communities counts real variants as errors
  ([dada2#100](https://github.com/benjjneb/dada2/issues/100)).

### Runs and batches (advice)

- Learn one model per sequencing run, or per batch whose PCR or sequencing may
  differ. Mixing runs gives rates "too high for one, too low for the other"
  ([dada2#110](https://github.com/benjjneb/dada2/issues/110),
  [dada2#240](https://github.com/benjjneb/dada2/issues/240); repeated in more
  than ten threads).
- Pool runs only if their models agree; otherwise use the higher model
  ([dada2#257](https://github.com/benjjneb/dada2/issues/257)).
- Different sample types or treatments within one run can share a model
  ([dada2#851](https://github.com/benjjneb/dada2/issues/851),
  [dada2#1079](https://github.com/benjjneb/dada2/issues/1079)). One user's run
  that mixed environments did better with a model per environment
  ([dada2#1317](https://github.com/benjjneb/dada2/issues/1317)).

### Signs of a poor fit (advice)

- "Terminated before convergence": in almost all cases, proceed
  ([FAQ](https://benjjneb.github.io/dada2/faq.html)).
- Overestimating rates is less harmful than underestimating them
  ([dada2#722](https://github.com/benjjneb/dada2/issues/722)).
- "Very few reads" or an empty matrix: drop near-empty samples, and remove
  poly-G reads on two-colour instruments
  ([dada2#469](https://github.com/benjjneb/dada2/issues/469),
  [dada2#2147](https://github.com/benjjneb/dada2/issues/2147)). On binned-quality
  data, see [binned quality](findings/binned-quality.md).

## Evidence in detail

- **Budget sweep on NovaSeq**
  ([dada2#791](https://github.com/benjjneb/dada2/issues/791#issuecomment-554588665)):
  1e8, 1e9 and 1e10 bases (1, 9 and 81 samples) gave similar models, except that
  the Q40 G→T rate fell from 3.3e-4 to 4.3e-6 at 1e10. No ASV comparison.
- **1e8 vs 5e8**
  ([dada2#1780](https://github.com/benjjneb/dada2/issues/1780)): 2881 vs 2871
  ASVs, 98.5% shared, described as "virtually indistinguishable".
- **Model per environment**
  ([dada2#1317](https://github.com/benjjneb/dada2/issues/1317)): one MiSeq-type
  run of sediment, seawater and fish samples. A sediment-only model kept about
  18% more sediment reads than a model trained on 15 random samples from the
  whole run; water gained about 5%, fish samples were unchanged.
- **Controls only**
  ([dada2#2032](https://github.com/benjjneb/dada2/issues/2032)): PacBio
  full-length 16S, training on Zymo controls gave a better-looking error plot
  than training on low-duplication saliva samples.
- **Single-sample batch**
  ([dada2#1830](https://github.com/benjjneb/dada2/issues/1830)): a model learned
  from one late sample made it an ordination outlier; processing it with the main
  batch removed the separation.

## Not measured

- Whether `--randomize` produces a better model, or only a different one.
- Low-diversity vs diverse vs same-environment training samples, against a truth
  set or on held-out samples (#105).
- How far error rates differ between runs, and what a pooled-run model costs.
- Any platform-specific budget. dada2-rs's concordance runners use 1e8 for
  Illumina and 2e8 for PacBio as provisional values
  ([issue #225](https://github.com/HPCBio/dada2-rs/issues/225)).
- Whether a deep sample should be downsampled or kept whole when it is one of
  several training samples ([issue #68](https://github.com/HPCBio/dada2-rs/issues/68)).
