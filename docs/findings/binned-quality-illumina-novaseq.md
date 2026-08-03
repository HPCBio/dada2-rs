# Binned quality on Illumina — NovaSeq 6000

**Verdict: on binned NovaSeq 6000 data, choosing the wrong `--errfun` is an
*abundance* error, not a richness error.** The two models find nearly the same
number of ASVs (+3.4%) but disagree substantially on *which* ones (Jaccard 0.706)
and heavily on *how much of each* (22.2% L1 on shared ASVs, 29.0% over the
union), and one arm carries 17% more reads into the final table than the other.

!!! warning "This page was rewritten after a full retraction"
    An earlier version reported a −26.1% ASV difference. **That was library-prep
    artifact**, as was a subsequent −16.7% partial correction. All figures here
    are from a clean re-run on `cutadapt`-trimmed reads. The failure is worth
    reading on its own account — see [Reading the prep before the
    result](reading-the-prep.md), which is the methodological companion to this
    page.

## Scope and provenance

Public data, **BioProject [PRJNA1504839]** — soil bacterial and fungal
communities in *Quercus liaotungensis* and *Pinus tabuliformis* forests in
Lanzhou; Illumina NovaSeq 6000, topsoil and subsoil. This page covers the
**bacterial 16S** arm only; 30 samples, 2,445,276 read pairs.

Qualifications:

- **The amplicon is V3–V4 (341F/805R), not V4.** The NCBI BioProject metadata
  states `Design: 515F-806R (V4)` and is **wrong** — established by direct primer
  search. Reads were trimmed with `cutadapt` against the actual primers
  (`dev/cutadapt_v3v4.sh`, `--discard-untrimmed`) before this analysis. See
  [Reading the prep](reading-the-prep.md).
- **Documented sampling, undocumented library prep.** PCR cycle counts, template
  quality, and index strategy are not published. Absolute error and chimera rates
  should not be read as typical for soil or compared across datasets. The
  **within-dataset A/B survives this** — everything but the errfun is held fixed.
- **No truth set.** This measures *disagreement between two models*, not the
  accuracy of either. Deciding which arm is right needs a mock community on
  binned Illumina — the same dependency blocking [#44].
- **NovaSeq 6000 specifically.** Illumina has since released the **NovaSeq X**,
  whose quality behaviour differs. If its mass is less spread across bins, the
  effect measured here should shrink toward the PacBio picture. **This result
  should not be assumed to hold on NovaSeq X.**

## Setup

Two arms over **identical** derep inputs — 1,225,523 (R1) / 1,304,044 (R2) merged
uniques from 2,445,276 reads, k=5, band 16, NW backend, `dada-pooled`.
`learn-errors` used the same pinned samples in both arms (nq=38, bins
`{2, 11, 25, 37}`).

| Arm | `--errfun` | Correct for binned input? |
|---|---|---|
| **default** | `loess` | no — the mismatch under test |
| **binned** | `binned-qual` | yes |

**Only the errfun differs.** Same reads, same derep, same denoising parameters,
same chimera settings.

## The error model, before any inference

All four `learn-errors` runs converge (7, 5, 9, 8 iterations), so this compares
two well-defined models. Off-bin mass is **25.4% (R1) / 23.9% (R2)** across 27 of
38 populated quality columns — far more spread than PacBio's ~8%. (As on PacBio,
off-bin values arise from derep quality-averaging across the reads collapsed into
a unique; they are not a binning-detection failure.)

Against `binned-qual`, the `loess` fit is:

| | R1 | R2 | PacBio HiFi (for contrast) |
|---|---|---|---|
| Mass-weighted mean error rate | **31.2% low** | **29.0% low** | 17.8% low |
| Ratio at the **dominant** bin (Q37) | **0.685** | **0.668** | 0.92 |
| Share of mass in that bin | 65.3% | 62.8% | ~91% |
| Worst trough | 86× at Q32 | 113× at Q32 | 16× at Q35 |
| Mass in the trough | 0.69% | 0.75% | ~8.4% |

The second row was long treated as decisive. On PacBio the `loess` error is small
where the data is dense, so the mismatch has nowhere to express itself; here it
lands squarely on the dominant bin, with `loess` underestimating the error rate
by ~32% exactly where two thirds of all observations sit. The 86–113× troughs
are attention-grabbing and nearly irrelevant — they carry under 1% of the mass.

!!! danger "This is a description of the model, not a working explanation"
    The **mass-concentration** account above is consistent with everything on
    this page, and it is nonetheless **not a predictor of errfun sensitivity**.
    The [ITS2 arm](binned-quality-illumina-its2.md) of this same deposit is *less*
    concentrated and its two error models disagree *twice as hard*, yet its
    inference diverges about a third as much. A replacement hypothesis —
    coverage per variant — was refuted by a depth ladder in the wrong direction.

    Read the model-level numbers here as a characterisation of what the two fits
    do, not as the reason the effect is large. The reason is unknown.

## Inference — four endpoints

| endpoint | `loess` | `binned-qual` | Δ | Jaccard | L1 (union) | L1 (shared) | reads |
|---|---|---|---|---|---|---|---|
| R1 | 10,541 | 9,701 | +8.7% | 0.817 | 27.8% | 19.1% | +11.9% |
| R2 | 11,851 | 11,283 | +5.0% | 0.869 | 10.4% | 6.8% | +3.1% |
| merged (uniques) | 37,355 | 36,315 | +2.9% | 0.741 | 28.1% | 22.2% | +18.1% |
| **non-chimeric** | **22,498** | **21,754** | **+3.4%** | **0.706** | **29.0%** | **22.2%** | **+17.0%** |

!!! note "Two L1 variants, because they answer different questions"
    L1 is the summed absolute difference in **relative** abundance per ASV, each
    table normalised to 1. It ranges 0–2, and **half of it is the total variation
    distance** — the fraction of community mass that would have to move to turn
    one table into the other.

    **Union L1** sums over all ASVs in either arm, scoring arm-specific ones as
    zero in the other. It is the honest "how different are these two tables"
    number, but it partly re-counts the composition difference already in
    Jaccard. **Shared L1** restricts to ASVs both arms found, renormalised within
    that subset, isolating *pure abundance* disagreement on the ASVs they agree
    exist. Shared mass is 96.7% / 95.4% of each arm, so the two are measuring
    nearly the same reads, differently attributed.

    Earlier versions of this page reported shared-only L1. Quote the variant you
    mean.

Length distributions are identical in both arms (median 426 bp, range 143–447),
so none of this is geometry.

> The merged row counts merged-pair uniques, not a denoised table — singletons
> appear there (12,034 / 11,884) where the denoised rows have none. It is not
> comparable *in kind* to the rows around it, though both arms are counted
> identically so the comparison across it is valid.

**Richness barely moves.** At +3.4% post-chimera the ASV-count difference is
close to a non-event, and lands near the ~0.3% predicted by the insert-only
collapse performed during the retraction. Anyone reading this A/B through `n_asv`
alone would conclude the errfun does not matter.

**Composition and abundance move a lot.** Jaccard 0.706 means roughly 30% of the
union of the two tables is arm-specific — 3.3% of `loess` reads and 4.7% of
`binned-qual` reads sit in ASVs the other arm never produced. Restricting to the
ASVs *both* arms found, they still disagree by **22.2% L1** on relative
abundance — an 11.1% total variation distance, on sequences whose existence is
not in dispute. That is a scale that changes ordination, differential-abundance
testing, and any abundance-weighted diversity metric.

For comparison, the post-hoc insert-only collapse during the retraction estimated
16.3% shared L1. The clean re-run puts it at **22.2%** on the same statistic — so
the abundance effect is somewhat *larger* than the collapse suggested, even as the
ASV-count effect vanished.

**Read retention differs more than anything else.** Of 2,445,276 input reads, the
`loess` arm carries **97.9%** into the R1 table against `binned-qual`'s **87.5%** —
a **10.4-point** gap, and 17.0% more reads in the final table. The binned arm is
not merely partitioning the same reads differently; it is **declining to place
reads the `loess` arm places**. This is a distinct axis from ASV count and is the
largest single difference between the arms.

### Direction

`loess` fits lower error rates, so variants are harder to explain away as
sequencing error, so it splits more and produces **more** ASVs. That is what
happens here at every endpoint, and it matches PacBio.

This is worth stating explicitly because the **pre-trimming analysis showed the
opposite** — `loess` produced consistently *fewer* ASVs, an inversion recorded at
the time as an unexplained anomaly. It was prep artifact. The direction now
agrees with the model across both platforms, and there is no anomaly left to
explain.

### Chimera removal does not rescue it

| | `loess` | `binned-qual` |
|---|---|---|
| Merge rate | 83.8% | 81.7% |
| Uniques removed as chimeric | 39.8% | 40.1% |
| Reads retained through chimera removal | 95.0% | 95.9% |
| Final ASVs | 22,498 | 21,754 |

The arms are pruned at **nearly identical rates**, so chimera removal is not
treating one arm harshly. It removes a comparable fraction from each and leaves
the disagreement standing: Jaccard moves only 0.741 → 0.706 across the step. The
sequences the arms disagree about are largely not the ones being called chimeric.

## What this dictates

1. **On binned Illumina, report the errfun in methods.** It is a
   result-changing parameter. The PacBio null does not transfer — but note *what*
   it changes: richness is nearly stable, so an errfun mismatch will pass a review
   that only checks ASV counts while corrupting every abundance-weighted
   downstream result.
2. **The [#98] guard should fire here, and its message should name binned
   Illumina.** This is the first measured case where warning the user prevents a
   real error rather than an invisible one. Detection is available from
   `count == 1` uniques, whose `qual_sum` is the raw unaveraged quality vector.
   Per the project's stated preference this stays a **strong warning, not an
   error**, so users can experiment deliberately.
3. **Do not predict the effect from mass concentration.** *(Revised.)* This page
   previously recommended exactly that, on the grounds that the fraction of mass
   in the dominant bin separates the PacBio null from the effect here. The
   [ITS2 arm](binned-quality-illumina-its2.md) refutes it: less concentrated,
   error models disagreeing twice as hard, a third the effect. The derived
   prediction that **i100 should behave like PacBio** (98.9% single-bin) is
   therefore an open question, not an expectation. "Illumina" is still the wrong
   unit of generalisation — but so, on current evidence, is every other unit we
   have tried.
4. **Do not treat `loess` non-convergence as the signal.** An earlier version of
   this page proposed exactly that, on the basis that 3 of 4 runs hit
   `max_consist`. After correct trimming all four converge, so non-convergence
   tracked prep quality, not errfun suitability. Retracted.
5. **Verify the prep before trusting any absolute count on unfamiliar data.**
   See [Reading the prep before the result](reading-the-prep.md) — this analysis
   was fully corrupted, and briefly *sign-reversed*, by an untrimmed library that
   raised no warnings anywhere in the pipeline.

## Path forward

1. **The ITS2 arm of the same BioProject — DONE**, and it changed the series.
   The effect is about a third the size, and it refuted the mechanism this page
   used to assert. See [Illumina NovaSeq 6000 — soil
   ITS2](binned-quality-illumina-its2.md). Note the arms are **separate PCRs**,
   so the difference is not cleanly attributable to the amplicon.
2. **NovaSeq X.** The most consequential gap. Until measured, this verdict is
   scoped to the 6000.
3. **i100 at inference.** i100 was characterised at the model level only —
   **98.9% of mass at Q38**, *more* concentrated than PacBio. This was previously
   listed as a clean falsifiable prediction that i100 would show the PacBio null;
   with mass concentration refuted as a predictor, it is now simply **unmeasured**
   and worth running on its own account.
4. **The fixed-reads errfun sweep** ([#98]) — `binned-qual` vs `loess` vs
   `pacbio` vs R's `loessErrfun_mod4` on one binned input. Isolates the errfun
   term exactly.
5. **A truth-anchored binned Illumina dataset.** Everything here measures
   disagreement. Deciding *which* arm is right needs a mock community on binned
   Illumina — the same dependency blocking [#44].
6. **A read-overlap QC screen ahead of the pipeline.** Tooling worth having
   regardless of the binned-quality question; see
   [Reading the prep](reading-the-prep.md).

[PRJNA1504839]: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1504839
[#98]: https://github.com/HPCBio/dada2-rs/issues/98
[#44]: https://github.com/HPCBio/dada2-rs/issues/44
