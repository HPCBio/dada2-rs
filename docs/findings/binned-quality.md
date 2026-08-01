# Binned quality scores — overview

Modern sequencers increasingly emit **binned** quality scores: instead of a
continuous Phred range, every base carries one of a handful of discrete values.
PacBio Revio emits 7 levels; Illumina i100 emits `{12, 24, 38}`; NovaSeq 6000
emits `{2, 11, 25, 37}`. DADA2's error model is fit *per quality value*, so
collapsing the quality axis changes what the model can express — and therefore
what gets denoised.

dada2-rs provides `--errfun binned-qual` for this case. These pages record what
we measured when we varied the error model on binned input.

## The headline: the answer is platform-dependent

This is the single most important thing to carry away, and it is why these pages
are split by platform rather than merged into one verdict.

| Platform | Effect of choosing the wrong errfun on binned input |
|---|---|
| **PacBio HiFi** (SequelIIe, Revio) | Small, and largely concealed by chimera removal. Post-chimera ASV sets agree to within 0.7%; the mismatched arm ends up **31 fewer** ASVs (−0.13%). |
| **Illumina NovaSeq 6000** | Richness barely moves (**+3.4%** ASVs) but composition and abundance do: Jaccard **0.706**, **29.0%** L1 divergence, and one arm carries **17% more reads** into the final table. An *abundance* error, not a richness error. |

Do not generalise either result to the other platform. And note **which axis
moves**: on NovaSeq an errfun mismatch would pass any review that checks only ASV
counts, while corrupting every abundance-weighted downstream analysis.

!!! warning "The Illumina result was retracted once before"
    Its first version reported a −26% ASV difference that proved to be
    library-prep artifact — retained degenerate primers behind heterogeneity
    spacers at both ends, inflating the table 3–4× and *reversing the direction*
    of the effect. The current figures are from a clean re-run on trimmed reads.
    The most transferable lesson in this whole series is the methodological one:
    **characterise primers and spacers before attributing anything to the error
    model.** See [Reading the prep before the result](reading-the-prep.md).

## Why the platforms differ

The mechanism is **mass concentration**, and it is worth understanding because it
predicts where the next platform will land.

An error model fit on binned data is wrong mainly *between* the bins, where the
LOESS fit interpolates across sparse support. Whether that matters depends on how
much observation mass actually sits in the interpolated region:

- **PacBio HiFi** concentrates ~91% of observation mass in the top bin. The badly
  interpolated region holds ~8% of the mass, so a fit that is 16× wrong there
  barely reaches the likelihood computation. The mismatch mostly *cannot express
  itself*.
- **Illumina** quality declines along the read, so mass is spread across bins.
  On the NovaSeq run measured here, **24–25% of mass is off-bin**, the dominant
  bin holds only ~63–65% of mass, and the error reaches that dominant bin
  directly (top-bin ratio **0.67–0.69**, versus 0.92 on PacBio). The fitting
  error lands where the reads actually are. Conversely the worst troughs — 86–113×
  — carry under 1% of the mass and are nearly irrelevant.

So the rule of thumb is: **the errfun choice matters in proportion to how much
quality mass sits away from a single dominant bin.**

A second-order contributor on Illumina: with few bins there are few anchors for
the fit. With only 3 populated columns, stock LOESS fits an exact quadratic
through three points, and the [#95] graceful-degradation guard (≤ 5 populated
columns) can fire — putting the fallback in genuinely degenerate territory rather
than merely mis-shaped.

## Reading these pages

Both pages report **pre- and post-chimera endpoints**, deliberately. Post-chimera
answers "does this change my results?"; pre-chimera is where the mechanism is
visible. Several times in this investigation a real pre-chimera difference
vanished downstream — collapsing to the endpoint alone would have reported null
results and hidden real effects that resurface under a different platform.

- [**PacBio (SequelIIe & Revio)**](binned-quality-error-model.md) — binning
  itself is safe and never cost reference recovery; the errfun mismatch is real
  at the model level but washes out downstream. Includes the truth-anchored mock
  community arms.
- [**Illumina NovaSeq 6000**](binned-quality-illumina-novaseq.md) — the errfun
  mismatch reaches the final table intact, as an abundance and read-retention
  effect rather than a richness one. Measures disagreement, not accuracy: there
  is no truth set on this dataset.
- [**Reading the prep before the result**](reading-the-prep.md) — the
  methodological companion to the NovaSeq page: how untrimmed degenerate primers
  produced a confident, reproducible, sign-reversed wrong answer, and the checks
  that catch it.

## A caveat that applies to both

Neither page tests **NovaSeq X**, **i100 at inference**, or **ITS/non-16S
amplicons** end to end. The NovaSeq result in particular is from a **NovaSeq
6000**; the newer NovaSeq X has different quality behaviour and may well shift
toward the PacBio picture. See each page's *Path forward*.

[#95]: https://github.com/HPCBio/dada2-rs/issues/95
