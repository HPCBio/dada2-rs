# Binned quality on Illumina — the ITS2 arm, and two dead mechanisms

**Verdict: the errfun mismatch is real on fungal ITS2 but roughly a third as
consequential as on bacterial 16S from the same sequencing run — and the two
mechanisms we had proposed to explain *when* it matters are both refuted by this
comparison.**

This page carries the negative result. It is the most useful page in the series
for deciding what to do next, and the least satisfying.

## Why this comparison exists

The [NovaSeq 6000 16S result](binned-quality-illumina-novaseq.md) showed a large
errfun effect. The obvious question is whether that generalises past bacterial
16S. The same BioProject sequenced a **fungal ITS2** arm on the same instrument,
which is the cheapest available test.

It is also the arm whose metadata is wrong in the *opposite* direction from the
16S arm — the region label `(ITS2)` is right while the primer names
`ITS1FI2-ITS2` are not; the real pair is **gITS7 + ITS4**. See
[Reading the prep before the result](reading-the-prep.md). Reads were trimmed
with `dev/cutadapt_ITS.sh` before anything here.

!!! warning "This is *not* a controlled amplicon comparison"
    It is tempting to describe the two arms as "same run, same prep, only the
    amplicon differs". That is wrong and the error matters. The two arms are
    **two separate PCRs** with different primer chemistry and, in all
    likelihood, different cycle counts. Nothing about the library preparation is
    documented: cycle number, template quality, amplification consistency across
    samples, and index strategy are all unknown, and this deposit has already
    concealed two separate prep traps.

    So "amplicon" here is confounded with PCR conditions and with whatever else
    differs between the two reactions. The comparison is worth making — it is
    the best available — but a difference between the arms **cannot be
    attributed to the amplicon alone**.

## Setup

30 samples, 2,503,123 read pairs, 825,214 merged R1 uniques — **identical derep
inputs in both arms**, `dada-pooled`, k=5, band 16, NW backend. Only `--errfun`
differs. All four `learn-errors` runs converge (loess 7/6 iterations,
`binned-qual` 4/4).

## The result

| endpoint | `loess` | `binned-qual` | Δ | Jaccard | L1 (union) | L1 (shared) | reads |
|---|---|---|---|---|---|---|---|
| R1 | 3,617 | 3,394 | +6.6% | 0.915 | 6.7% | 6.0% | +3.7% |
| R2 | 3,577 | 3,461 | +3.4% | 0.905 | 10.2% | 7.6% | +5.7% |
| **non-chimeric** | **3,181** | **3,020** | **+5.3%** | **0.928** | **11.1%** | **10.7%** | **+11.8%** |

Against 16S on the same instrument (Jaccard **0.706**, L1 **29.0%**, reads
**+17.0%**), every disagreement axis is substantially smaller. The direction is
unchanged — `loess` fits lower error rates, splits more, and produces more ASVs —
and the read-retention gap persists, which makes it the most robust signature in
the series: on both amplicons `binned-qual` declines to place reads that `loess`
places.

Two structural differences are worth recording because they are large:

- **Merging is nearly free.** 97.1% / 94.7%, against 83.8% / 81.7% on 16S. The
  ITS2 product is ~271 bp on 2×250, leaving ~230 nt of overlap where V3–V4
  leaves ~57. Merged lengths run 152–443 bp (median 271) — that spread is real
  fungal ITS length variation, and its presence is itself confirmation the
  primer trimming was correct.
- **Chimera removal is almost a no-op**, removing 0.1–0.9% of uniques across the
  depth ladder against roughly 40% on 16S. Read this cautiously: chimera
  detection is *harder* on ITS, which has no conserved anchors and highly
  variable length, so a low call rate is not evidence of few chimeras.

## Mechanism 1 — mass concentration — predicts the wrong sign

The series overview previously explained the PacBio/Illumina split by **mass
concentration**: an errfun mismatch can only express itself where observation
mass actually sits, so platforms concentrating mass in one quality bin should be
insensitive. ITS2 tests it directly and breaks it.

| | 16S R1 | **ITS2 R1** | PacBio HiFi |
|---|---|---|---|
| Off-bin mass | 25.4% | **40.9%** | ~8% |
| Share in dominant bin (Q37) | 65.3% | **49.6%** | ~91% |
| `loess`/`binned` error ratio at that bin | 0.685 | **0.330** | 0.92 |
| Mass-weighted mean error, `loess` vs `binned` | 31.2% low | **38.6% low** | 17.8% low |
| **Resulting inference disagreement (Jaccard)** | **0.706** | **0.928** | ~0.99 |

Every input to the model says ITS2 should be **more** sensitive, not less. Its
quality mass is more spread, and the two error models disagree with each other
roughly *twice as hard* — `loess` underestimates the error rate at the dominant
bin by a factor of three, against 1.5× on 16S. The inference nonetheless
disagrees far less.

**Mass concentration is not a sufficient predictor of errfun sensitivity.** It
may still be necessary — PacBio's null is consistent with it — but it cannot
carry the explanatory weight the overview gave it, and any prediction derived
from it alone (including the standing prediction that i100 will behave like
PacBio) is now unsupported.

## Mechanism 2 — coverage per variant — predicts the wrong direction

The replacement hypothesis: the error model only changes calls sitting near the
OMEGA_A abundance boundary, so what should matter is how well-supported each
decision is. 16S carries a **median 19 reads per ASV** against ITS2's **60**, and
its 22,498 ASVs against ITS2's 3,181 at comparable depth. That predicts the
effect should **appear** in ITS2 once it is rarefied to 16S-like coverage.

Rarefying is also the only direction available — you cannot make a dataset
better-covered by discarding reads — and predicting an effect to emerge in the
near-null arm is the falsifiable direction, so ITS2 was the right arm to
subsample.

Method: nested subsampling at 50 / 25 / 12.5% (`seqkit sample` at fixed seed;
`dev/subsample_fastq.py` is the dependency-free equivalent), with the
**full-depth error models reused at every rung** via `--error-model` +
`--inherit-err-params`. Refitting `learn-errors` per rung would have confounded
the mechanism with the fit. No parameter-disagreement warnings fired, so the
full-depth run is a valid 100% rung.

| frac | `loess` | `binned` | ΔASV | Jaccard | L1u | L1s | read Δ | **median coverage** |
|---|---|---|---|---|---|---|---|---|
| 1.000 | 3,181 | 3,020 | +5.3% | 0.928 | 11.1% | 10.7% | +11.8% | 60 |
| 0.500 | 2,562 | 2,447 | +4.7% | 0.930 | 11.1% | 10.4% | +10.2% | 45 |
| 0.250 | 1,961 | 1,833 | +7.0% | 0.912 | 9.4% | 8.5% | +7.3% | 34 |
| 0.125 | 1,433 | 1,338 | +7.1% | 0.904 | 8.6% | 7.3% | +5.3% | 29 |

An 8× depth cut drags median coverage from 60 to 29 — the same order as 16S's 19
— and Jaccard moves only 0.928 → 0.904. It is not heading for 0.706.

**Worse, L1 and read retention move the wrong way entirely**: divergence
*shrinks* as depth falls (11.1% → 8.6%; +11.8% → +5.3%). The R1 endpoint agrees
more sharply (L1u 6.7% → 3.8%). A hypothesis that predicts growth and observes
decay is refuted, not merely unsupported.

Sanity checks: read counts halve cleanly (1,251,572 → 618,297 → 308,848) and the
top-ASV share is invariant across rungs (11.9% → 12.2%), so rarefaction removed
evidence without reshaping the community.

### What the direction does tell us

Starving the data makes the two arms agree **more**. The natural reading is that
divergence tracks **how many denoising decisions are made at all**, not how
well-supported each one is: with less data fewer uniques clear the abundance
threshold, both arms make the same conservative calls, and the models have
nowhere to disagree.

That is a real observation, but it does not rescue a predictor. Rarefaction cuts
ITS2's decision count too, and the effect still did not grow toward 16S.

## A third candidate, unresolved

The remaining idea is structural: **how densely variants are packed in sequence
space.** Soil 16S is full of near-identical taxa a substitution or two apart —
exactly where the error rate decides the call — while ITS2 is a fast-evolving
spacer whose real taxa are many substitutions apart. Rarefaction cannot change
that, which would explain a flat ladder.

Measured as nearest-neighbour Hamming distance to a more abundant ASV:

| | median NN distance | at distance 1 | ≤ 2 | ≤ 5 |
|---|---|---|---|---|
| 16S | **7** | 18.5% | 29.7% | 46.2% |
| ITS2 | **83** | 18.0% | 23.8% | 30.2% |

**This is not yet evidence.** The median difference is 12× and in the predicted
direction, but the **distance-1 fraction is identical** (18.5% vs 18.0%) and
distance 1–2 is where the error model actually decides. The profile also
compares the top 1,500 ASVs of each table, which is 6.7% of the 16S table
against 47% of the ITS2 table — different slices of very different
distributions. A real test needs **matched abundance strata**, and until then
this belongs in the hypothesis column.

## What this dictates

1. **We cannot currently predict which datasets are errfun-sensitive.** Two
   mechanisms proposed, two refuted, one unresolved and confounded. This is the
   honest state and it makes the practical recommendation *stronger*, not
   weaker: nobody gets to skip the check on the grounds that their data "looks
   concentrated."
2. **On binned Illumina, report the errfun in methods, and test both arms where
   the result matters.** This survives both refutations because it never
   depended on either mechanism.
3. **The effect is not a property of "Illumina" or of binning alone.** Two
   amplicons off one instrument differ threefold in sensitivity. Whatever the
   governing variable is, it is not the sequencer.
4. **Retract mass concentration as the series' explanation.** It remains a
   plausible *necessary* condition and it is consistent with the PacBio null,
   but it does not predict sensitivity, and the i100 prediction derived from it
   should be treated as an open question rather than an expectation.
5. **Do not re-run the coverage ladder.** It is answered. Rarefying reduces
   divergence on this dataset; that avenue is closed.
6. **Read-retention divergence is the most robust axis.** It is the one effect
   that persists across both amplicons and every rung of the ladder. If a single
   diagnostic is ever built for errfun mismatch, this is the quantity to watch —
   not ASV count, which is nearly invariant on both arms.

## Scope

No truth set, so this measures **disagreement between two models, not the
accuracy of either**. NovaSeq 6000 specifically. Library preparation is
undocumented; absolute rates should not be read as typical for soil fungi or
compared across datasets. And per the warning above, the 16S/ITS2 contrast is
confounded with PCR conditions and cannot be attributed to amplicon identity
alone.
