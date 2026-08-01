# Binned quality on Illumina — NovaSeq 6000

**Verdict:** on binned Illumina NovaSeq 6000 reads, the choice of `--errfun`
**changes the published table**. `loess` and `binned-qual` on identical inputs
produce final non-chimeric tables that differ by **−16.7% in ASV count**, share
only **58%** of their ASVs (Jaccard 0.577), and disagree by **~27% in abundance
even on the ASVs they share**. Chimera removal does **not** absorb the
difference.

Those figures are after correcting for a library-prep artifact in this dataset
(retained heterogeneity spacers + forward primer, which alone inflate ASV counts
by 39–46%). Uncorrected the gap looks larger still — −26.1%, Jaccard 0.43 — but
that is not all errfun. See [Prep artifacts](#prep-artifacts-read-this-before-the-numbers).

This is the opposite of the [PacBio result](binned-quality-error-model.md), where
the same mismatch was small and concealed downstream. See
[the overview](binned-quality.md) for why the platforms diverge.

!!! warning "This measures disagreement, not accuracy"
    There is no truth set for this dataset, and the `loess` arm **failed to
    converge**. Nothing here establishes which arm is closer to the true
    community — only that the two arms cannot both be right, and that the choice
    is not absorbed downstream. Do not read the ASV counts as one arm
    "over-calling."

## Scope and provenance

Public data, **BioProject [PRJNA1504839]** — soil bacterial and fungal
communities in *Quercus liaotungensis* and *Pinus tabuliformis* forests in
Lanzhou; Illumina NovaSeq 6000, 16S rRNA V4 (bacteria) and ITS2 (fungi), topsoil
and subsoil. This page covers the **bacterial 16S** arm only.

Qualifications:

- **Documented sampling, undocumented library prep.** The sampling design is
  published, but PCR cycle counts, template quality, and index strategy are not.
  Absolute error rates and chimera rates should not be read as typical for soil,
  and should not be compared across datasets. The **within-dataset A/B survives
  this** — everything but the errfun is held fixed.
- **The NCBI design annotation is wrong.** BioProject metadata states
  `Design: 515F-806R (V4)`. The data are **V3–V4 (341F/805R, ~444 bp)**, shown
  by direct primer search: 341F is present in **99.96%** of merged reads at the
  5' end, and the 515F site sits ~150 bp *internally*. If these were V4, 515F
  would be at position 0 and 341F absent. See
  [Prep artifacts](#prep-artifacts-read-this-before-the-numbers).
- **NovaSeq 6000 specifically.** Illumina has since released the **NovaSeq X**,
  whose quality behaviour differs. If its mass is less spread across bins, the
  effect measured here should shrink toward the PacBio picture. **This result
  should not be assumed to hold on NovaSeq X.**

## Setup

Two arms over **identical** derep inputs — 30 samples, 898,604 merged uniques,
1,750,945 reads, k=5, band 16, NW backend, `dada-pooled`. `learn-errors` used the
same 8 pinned samples in both arms (nq=38, bins `{2, 11, 25, 37}`).

| Arm | `--errfun` | Correct for binned input? |
|---|---|---|
| **default** | `loess` | no — the mismatch under test |
| **binned** | `binned-qual` | yes |

**Only the errfun differs.** Same reads, same derep, same denoising parameters,
same chimera settings.

## The error model, before any inference

Off-bin mass is **15.0% (R1) / 16.3% (R2)**, with 26 of 38 quality columns
populated — far more spread than PacBio's ~8%. (As on PacBio, off-bin values come
from derep quality-averaging across the reads collapsed into a unique; they are
not a binning-detection failure.)

Against `binned-qual`, the `loess` fit is:

| | R1 | R2 | PacBio HiFi (for contrast) |
|---|---|---|---|
| Mass-weighted mean error rate | **23.4% low** | **26.6% low** | 17.8% low |
| Ratio at the **top** (dominant) bin | **0.85** | **0.83** | 0.92 |
| Worst trough | — | **25× at Q32** | 16× at Q35 |
| Mass in the trough region | 13.9% | 15.1% | ~8.4% |

The decisive row is the second. On PacBio the error is small where the data is
dense; here it reaches the dominant bin directly.

### Non-convergence

| Arm | R1 | R2 |
|---|---|---|
| `loess` | did not converge (`max_consist` 10) | **did not converge**, oscillating non-monotonically (1.05e-4 → 2.86e-4, *rising*) |
| `binned-qual` | did not converge, but residual 10–20× smaller | **converged to exactly 0 in 6 iterations** |

Every PacBio run in this investigation converged in 4–5 iterations under both
errfuns. The `loess` arm's final iterate is therefore not a well-behaved model,
which is part of why its *direction* of effect resists explanation (below).

## Prep artifacts — read this before the numbers

Two properties of this library, both discovered from the merge output, affect
every absolute count on this page. Both are shared by the two arms, so the A/B
survives; neither should be ignored when reading magnitudes.

**1. The amplicon is V3–V4, not V4.** Primer search across merged reads:

| motif | found in | median position |
|---|---|---|
| **341F** `CCTACGGGNGGCWGCAG` | **99.96%** | position 2 |
| **515F** `GTGYCAGCMGCCGCGGTAA` | 98.08% | **position 152** (internal) |
| 806R / 805R (rev-comp) | 0.00% | — |

The geometry is fully self-consistent: V3–V4 is ~464 bp with primers; less the
21 bp reverse primer gives 443 bp, the observed modal merged length. Reads are
uniformly 240+240 (`R1len + R2len = 480` for every merge), so 480 − 443 = **37 bp
expected overlap**, matching the observed median of 35.

**2. Heterogeneity spacers and the forward primer are retained.** 341F appears at
**5 distinct offsets (0–4 nt)** — Fadrosh-style spacers that boost base diversity
on patterned flow cells. The reverse primer *was* trimmed; the forward was not.
Consequently the same biological sequence appears at several offsets as several
distinct ASVs. Stripping spacer + 341F collapses **39.3%** (`loess`) / **46.1%**
(`binned-qual`) of post-chimera ASVs.

!!! danger "`filter-and-trim` cannot fix this — use cutadapt"
    `filter-and-trim` removes leading bases by **fixed offset only**
    (`--trim-left`); it does not match primer sequence. With the primer at 5
    different offsets there is no correct value: too small leaves spacer bases,
    too large eats real bases, and both split one biological sequence into
    several ASVs. Primer removal for designs like this must happen upstream with
    `cutadapt` (or R's `removePrimers()`), with `filter-and-trim` used for
    quality filtering only. Tracked in [#113].

!!! warning "All headline figures on this page are the corrected ones"
    Comparisons below are computed on **spacer- and primer-stripped** sequences.
    This is a *post-hoc* collapse, not a re-run: properly, the reads should be
    re-trimmed and re-denoised, which would not give identical results. Treat the
    corrected numbers as a close approximation, not an exact re-analysis.

    | post-chimera | as-is | **stripped** |
    |---|---|---|
    | ASV Δ | −26.1% | **−16.7%** |
    | Jaccard | 0.433 | **0.577** |
    | binned-only reads | 31.3% | **16.8%** |
    | shared-ASV L1 | 29.4% | **27.2%** |

    About a third of the apparent errfun effect was prep artifact. The remainder
    is large and real.

## Inference — four endpoints

| endpoint | `loess` ASVs | `binned-qual` ASVs | Δ | Jaccard | binned-only reads | shared-ASV L1 |
|---|---|---|---|---|---|---|
| R1 | 8,186 | 9,120 | −10.2% | 0.735 | 13.4% | 19.9% |
| R2 | 11,890 | 14,574 | −18.4% | 0.716 | 13.3% | 14.7% |
| merged | 47,763 | 62,892 | −24.1% | 0.545 | 25.2% | 31.0% |
| non-chimeric (as-is) | 26,097 | 35,314 | −26.1% | 0.433 | 31.3% | 29.4% |
| **non-chimeric (corrected)** | **15,849** | **19,019** | **−16.7%** | **0.577** | **16.8%** | **27.2%** |

The divergence **grows monotonically down the pipeline**. It is not a
low-abundance tail: at R1 and R2 the arm-specific ASVs contain **zero
singletons**, with medians of 48–96 reads.

The corrected row strips spacer + forward primer (see
[Prep artifacts](#prep-artifacts-read-this-before-the-numbers)); **the bolded row
is the one to quote.** The R1, R2 and merged rows are as-is and are inflated by
the same artifact — the correction was computed only at the post-chimera
endpoint, so the intermediate magnitudes are upper bounds, though the monotonic
trend holds regardless.

### The merge rate is amplicon geometry, not library quality

Both arms merge at essentially the same rate — **44.09%** (`loess`) vs **42.95%**
(`binned-qual`), same per-sample spread (min 35.3%, max 51.7%). So the loss is
**upstream of the error model** and shared by both arms.

With only ~37 bp of nominal overlap on a 443 bp product, quality trimming eats
directly into the margin and pairs fall under `minOverlap`'s floor of 12. **A
~44% merge rate is unremarkable for V3–V4 on 2×250.**

Library background is **disconfirmed**: 99.97% of merges have overlap < 50 bp in
a tight unimodal distribution, merged lengths are unimodal at 443–449 bp with a
hard floor at 277 bp and **0.00% below 200 bp**, and **100.0%** of accepted
merges carry **zero mismatches** in the overlap. Primer dimer and mis-amplification
would produce short fragments and a multimodal profile. Neither is present.

!!! tip "Screen read overlap *before* the pipeline, not after"
    Everything in this section was recovered retrospectively from merge output.
    A per-sample overlap/fragment-length heatmap against the **expected amplicon
    size**, run as QC ahead of denoising, would have surfaced all of it up front:
    the wrong primer pair, the retained spacers, and the thin overlap budget. It
    applies to every paired run, independent of the binned-quality question.

But merging **amplifies** the divergence: Jaccard drops 0.716 → 0.545. A merged
sequence is determined by *both* its R1 and R2 ASV, so the two independent
disagreements compose rather than average.

> Note the merged row counts merged-pair uniques, not a denoised table —
> singletons appear there (3,911 / 6,874) where R1 and R2 had none. It is not
> comparable *in kind* to the rows above it, though both arms are counted
> identically so the trend across it is real.

### Chimera removal does not rescue it

| | `loess` | `binned-qual` |
|---|---|---|
| ASVs removed as chimeric | 45.4% | 43.8% |
| Reads retained | 78.4% | 78.7% |
| Final ASVs | 26,097 | 35,314 |

The two arms are pruned at **nearly the same rate** — chimera removal is not
treating one arm harshly. It removes a comparable fraction from each and leaves
the disagreement standing, because the sequences the arms disagree about are
largely *not* the ones it calls chimeric. On as-is sequences post-chimera Jaccard
(0.433) is even *lower* than pre-chimera (0.545); on corrected sequences the
post-chimera figure is 0.577.

On PacBio, chimera removal reduced a +0.21% pre-chimera effect to −0.13%. Here
the difference persists at −16.7% (corrected) after removal.

## Two observations without explanations

Recorded because they are real and reproducible, not because they are understood.

**1. The direction is unexplained.** `loess` fits *lower* error rates, which
should make variants harder to explain as sequencing error and therefore produce
*more* ASVs — the direction seen on PacBio. Instead `loess` produces
**consistently fewer** ASVs at every endpoint. A mass-weighted average is the
wrong statistic for predicting this: the abundance p-value depends on the error
rate at the specific transition and quality of the differing position, not on the
average. Combined with the `loess` arm's non-convergence, there is no basis here
for a directional rule.

**2. `loess` under-splits *and* retains more reads.**

| | reads in table | % of 1,750,945 |
|---|---|---|
| R1 `loess` | 1,703,269 | 97.28% |
| R1 `binned-qual` | 1,605,329 | 91.68% |
| R2 `loess` | 1,702,811 | 97.25% |
| R2 `binned-qual` | 1,684,125 | 96.18% |

The binned arm makes more ASVs *and* discards more reads (5.6 points more at R1).
It is therefore not merely partitioning the same reads more finely — it is
dropping reads the `loess` arm keeps. This is a distinct axis from the ASV count.

## What this dictates

0. **Check the amplicon and primer handling before trusting any absolute count.**
   This dataset's NCBI annotation names the wrong primer pair, and retained
   spacers plus forward primer inflate ASV counts by 39–46%. Both were invisible
   until the merge output was examined, and `--trim-left` cannot remove a
   variable-offset primer ([#113]). See
   [Prep artifacts](#prep-artifacts-read-this-before-the-numbers).
1. **On binned Illumina, the errfun is a result-changing parameter, not a
   diagnostic detail.** It must be chosen deliberately and reported in methods.
   The PacBio null does not transfer.
2. **The [#98] guard should fire here, and its message should name binned
   Illumina.** This is the first measured case where warning the user prevents a
   real error rather than an invisible one. Detection is available from
   `count == 1` uniques, whose `qual_sum` is the raw unaveraged quality vector.
   Per the project's stated preference this stays a **strong warning, not an
   error**, so users can still experiment deliberately.
3. **Treat `loess` non-convergence as a reportable condition.** Both `loess` arms
   hit `max_consist` with R2 oscillating upward, while `binned-qual` R2 converged
   cleanly. Non-convergence is currently recorded in the artifact
   (`stop_reason`) but is easy to miss; it deserves surfacing, and it is an
   argument for the guard independent of which arm one considers correct.
4. **A community report that `binned-qual` inflates ASV counts is *consistent
   with* this run** — `binned-qual` yields 10–26% more ASVs than `loess` at every
   endpoint here, whereas on PacBio binning fixed reads either cut ASVs 4× or left
   them within 0.3%. This raises Illumina from one of three hypotheses to the
   likely setting for that report. It is **not confirmation**: this A/B is not the
   clean test, which is the fixed-reads errfun sweep in [#98].

## Path forward

1. **The ITS2 arm of the same BioProject.** Same study, same sequencer, same
   library prep — so a difference there is amplicon-driven, not study-to-study
   variation. This is the cheapest test of whether the result generalises past
   bacterial 16S.
2. **NovaSeq X.** The most consequential gap. If its quality mass is less spread,
   the effect should shrink toward PacBio's. Until measured, this page's verdict
   is scoped to the 6000.
3. **i100 at inference.** i100 was characterised at the model level only —
   **98.9% of mass at Q38**, *more* concentrated than PacBio, against the
   prediction that Illumina would spread. i100 may therefore behave like HiFi and
   not like NovaSeq, which would make "Illumina" the wrong unit of generalisation
   entirely — the right one being mass concentration. Worth running to inference.
4. **The fixed-reads errfun sweep** ([#98]) — `binned-qual` vs `loess` vs
   `pacbio` vs R's `loessErrfun_mod4` on one binned input. Isolates the errfun
   term exactly and is the clean test of the ASV-inflation report.
5. **A truth-anchored binned Illumina dataset.** Everything here measures
   disagreement. Deciding *which* arm is right needs a mock community on binned
   Illumina — the same dependency blocking [#44].
6. **A read-overlap QC screen ahead of the pipeline.** The merge rate here was
   diagnosed retrospectively, from merge output, and was initially misread as a
   possible library-quality problem before the length distribution showed it was
   amplicon geometry. A per-sample overlap/fragment-length heatmap against the
   expected amplicon size answers that up front, and would also have caught the
   V3-V4/V4 metadata discrepancy before denoising rather than after. Worth having
   as tooling regardless of the binned-quality question, since it applies to
   every paired run.

[PRJNA1504839]: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1504839
[#98]: https://github.com/HPCBio/dada2-rs/issues/98
[#44]: https://github.com/HPCBio/dada2-rs/issues/44
[#113]: https://github.com/HPCBio/dada2-rs/issues/113
