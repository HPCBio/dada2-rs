# Binned quality scores — overview

Modern sequencers increasingly emit **binned** quality scores: instead of a
continuous Phred range, every base carries one of a handful of discrete values.
PacBio Revio emits 7 levels; Illumina i100 emits `{12, 24, 38}`; NovaSeq 6000
emits `{2, 11, 25, 37}`. DADA2's error model is fit *per quality value*, so
collapsing the quality axis changes what the model can express — and therefore
what gets denoised.

dada2-rs provides `--errfun binned-qual` for this case. These pages record what
we measured when we varied the error model on binned input.

## The headline: the answer is dataset-dependent, and we cannot yet predict it

This is the single most important thing to carry away, and it is why these pages
are split rather than merged into one verdict.

| Dataset | Effect of choosing the wrong errfun on binned input |
|---|---|
| **PacBio HiFi** (SequelIIe, Revio) | Small, and largely concealed by chimera removal. Post-chimera ASV sets agree to within 0.7%; the mismatched arm ends up **31 fewer** ASVs (−0.13%). |
| **NovaSeq 6000 — soil 16S V3–V4** | Richness barely moves (**+3.4%** ASVs) but composition and abundance do: Jaccard **0.706**, **22.2%** L1 abundance divergence on shared ASVs, and one arm carries **17% more reads** into the final table. An *abundance* error, not a richness error. |
| **NovaSeq 6000 — soil ITS2** | Same instrument, same deposit, **roughly a third the effect**: Jaccard **0.928**, L1 **10.7%** shared, **+11.8%** reads. |

Do not generalise any of these to the others. Note **which axis moves**: on the
16S arm an errfun mismatch would pass any review that checks only ASV counts,
while corrupting every abundance-weighted downstream analysis.

The third row is the awkward one, and it is why this section no longer says
"platform-dependent". **Two amplicons off one instrument differ threefold in
sensitivity**, so the sequencer is not the unit of generalisation — and, as of
this writing, neither is anything else we have been able to measure.

!!! warning "The Illumina result was retracted once before"
    Its first version reported a −26% ASV difference that proved to be
    library-prep artifact — retained degenerate primers behind heterogeneity
    spacers at both ends, inflating the table 3–4× and *reversing the direction*
    of the effect. The current figures are from a clean re-run on trimmed reads.
    The most transferable lesson in this whole series is the methodological one:
    **characterise primers and spacers before attributing anything to the error
    model.** See [Reading the prep before the result](reading-the-prep.md).

## Why the datasets differ — two refuted explanations

!!! danger "Mass concentration was this section's answer. It is refuted."
    Earlier versions of this page explained the split via **mass concentration**
    and offered it as a predictor of where the next platform would land. The
    ITS2 arm refutes it, and a depth ladder refutes the replacement hypothesis
    too. Both are documented below because a ruled-out mechanism is worth as
    much as a confirmed one — it stops the idea being re-litigated.

**Mass concentration** was the proposal that an errfun mismatch can only express
itself where observation mass actually sits, so a platform concentrating mass in
one quality bin is insensitive by construction. It fits PacBio (~91% of mass in
the top bin, mismatch washes out) and it fits soil 16S (mass spread, dominant bin
only 63–65%, mismatch lands where the reads are).

It fails on ITS2 **in the wrong direction**. That arm is *less* concentrated
(40.9% off-bin, dominant bin 49.6%) and its two error models disagree roughly
twice as hard (`loess`/`binned` ratio at the dominant bin **0.330** against 16S's
0.685) — yet its inference diverges about a third as much. Every input to the
model predicts higher sensitivity; the outcome is lower.

**Coverage per variant** was the replacement: the error model only flips calls
near the OMEGA_A abundance boundary, so what should matter is how well-supported
each decision is. 16S carries a median 19 reads per ASV against ITS2's 60. That
predicts the effect should appear in ITS2 once rarefied.

It fails **in the wrong direction too**. Rarefying ITS2 8×, to median coverage
29, moves Jaccard only 0.928 → 0.904, and abundance divergence *shrinks* rather
than grows (L1 11.1% → 8.6%; read gap +11.8% → +5.3%). Less data makes the two
arms agree **more**, apparently because fewer uniques clear the abundance
threshold at all, so fewer decisions get made for the models to disagree about.

Details and the full ladder are on the [ITS2 page](binned-quality-illumina-its2.md).

**Where that leaves prediction.** Mass concentration may still be *necessary* —
nothing here contradicts the PacBio null — but it is not *sufficient*, and no
prediction should rest on it alone. In particular the standing expectation that
**i100 will behave like PacBio** because it is 98.9% single-bin is now an open
question rather than a forecast. A third hypothesis, sequence-space density, is
under investigation in [#118] and is not yet evidence.

A second-order contributor on Illumina, unaffected by the above: with few bins
there are few anchors for the fit. With only 3 populated columns, stock LOESS
fits an exact quadratic through three points, and the [#95] graceful-degradation
guard (≤ 5 populated columns) can fire — putting the fallback in genuinely
degenerate territory rather than merely mis-shaped.

## What to do in the absence of a mechanism

The practical guidance never depended on the mechanism, and losing it makes the
advice **stronger** rather than weaker:

- **Report the errfun in methods** on any binned dataset.
- **Test both arms where the result matters.** No dataset property we can
  currently measure licenses skipping this — "my data looks concentrated" is
  exactly the reasoning ITS2 falsified.
- **Watch read retention, not ASV count.** The read-retention gap is the one
  axis that persisted across both amplicons and every rung of the depth ladder;
  ASV counts were nearly invariant on both arms and would have reported a null.

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
- [**Illumina NovaSeq 6000 — soil 16S**](binned-quality-illumina-novaseq.md) —
  the errfun mismatch reaches the final table intact, as an abundance and
  read-retention effect rather than a richness one. Measures disagreement, not
  accuracy: there is no truth set on this dataset.
- [**Illumina NovaSeq 6000 — soil ITS2**](binned-quality-illumina-its2.md) — the
  same instrument at a third the sensitivity, plus the depth ladder. This is the
  page that refutes both proposed mechanisms; read it before trusting any
  prediction about a new dataset.
- [**Reading the prep before the result**](reading-the-prep.md) — the
  methodological companion to the NovaSeq pages: how untrimmed degenerate primers
  produced a confident, reproducible, sign-reversed wrong answer, and the checks
  that catch it.

## A caveat that applies to all of them

No page tests **NovaSeq X** or **i100 at inference** end to end. The Illumina
results are from a **NovaSeq 6000**; NovaSeq X has different quality behaviour.
i100 was previously expected to resemble PacBio on mass-concentration grounds —
that expectation no longer has support and needs measuring.

A caveat specific to the two Illumina pages: they come from **one deposit whose
library preparation is undocumented**, and the 16S and ITS2 arms are **separate
PCRs**, so amplicon identity is confounded with PCR conditions. That deposit has
already concealed two distinct prep traps (see
[Reading the prep](reading-the-prep.md)). Treat cross-arm differences as real but
not cleanly attributable.

[#95]: https://github.com/HPCBio/dada2-rs/issues/95
[#118]: https://github.com/HPCBio/dada2-rs/issues/118
