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
| **Illumina NovaSeq 6000** | **Large, and it survives chimera removal.** Final tables differ by **−16.7%** in ASV count with post-chimera Jaccard **0.58** and ~27% abundance churn. |

Do not generalise either result to the other platform. A reader who takes the
PacBio null as reassurance will make a materially wrong table on NovaSeq.

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
  On the NovaSeq run measured here, **15–16% of mass is off-bin** and the error
  reaches the *dominant* bin (top-bin ratio 0.83–0.85, versus 0.92 on PacBio).
  The fitting error lands where the reads actually are.

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
  mismatch reaches the final table intact. Measures disagreement, not accuracy:
  there is no truth set on this dataset.

## A caveat that applies to both

Neither page tests **NovaSeq X**, **i100 at inference**, or **ITS/non-16S
amplicons** end to end. The NovaSeq result in particular is from a **NovaSeq
6000**; the newer NovaSeq X has different quality behaviour and may well shift
toward the PacBio picture. See each page's *Path forward*.

[#95]: https://github.com/HPCBio/dada2-rs/issues/95
