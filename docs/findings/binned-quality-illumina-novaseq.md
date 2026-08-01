# Binned quality on Illumina — NovaSeq 6000

**Verdict:** on binned Illumina NovaSeq 6000 reads, the choice of `--errfun`
**changes the published table**. `loess` and `binned-qual` on identical inputs
produce final non-chimeric tables that differ by **−26.1% in ASV count**, share
fewer than half their ASVs (**Jaccard 0.43**), and disagree by **~29% in
abundance even on the ASVs they share**. Chimera removal does **not** absorb the
difference — it widens it.

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
and subsoil. This page covers the **16S V4** arm only.

Two important qualifications:

- **Documented sampling, undocumented library prep.** The sampling design is
  published, but PCR cycle counts, template quality, and index strategy are not.
  Absolute error rates and chimera rates should not be read as typical for soil,
  and should not be compared across datasets. The **within-dataset A/B survives
  this** — everything but the errfun is held fixed. The low merge rate suggests
  these libraries may carry substantial background; see
  [the merge section](#the-merge-compounds-it).
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

## Inference — four endpoints

| endpoint | `loess` ASVs | `binned-qual` ASVs | Δ | Jaccard | binned-only reads | shared-ASV L1 |
|---|---|---|---|---|---|---|
| R1 | 8,186 | 9,120 | −10.2% | 0.735 | 13.4% | 19.9% |
| R2 | 11,890 | 14,574 | −18.4% | 0.716 | 13.3% | 14.7% |
| merged | 47,763 | 62,892 | −24.1% | 0.545 | 25.2% | 31.0% |
| **non-chimeric** | **26,097** | **35,314** | **−26.1%** | **0.433** | **31.3%** | **29.4%** |

The divergence **grows monotonically down the pipeline**. It is not a
low-abundance tail: at R1 and R2 the arm-specific ASVs contain **zero
singletons**, with medians of 48–96 reads.

### The merge compounds it

Both arms merge at essentially the same rate — **44.09%** (`loess`) vs **42.95%**
(`binned-qual`), same per-sample spread (min 35.3%, max 51.7%). So the ~44% merge
loss is **upstream of the error model**, shared by both arms, and is a property
of the reads rather than the errfun.

!!! warning "That rate is low for V4, and sample quality is not ruled out"
    16S V4 on 2×250 should overlap generously; losing over half the pairs is not
    normal. Because these are SRA data with **undocumented library prep**, we
    cannot distinguish sequencing quality from library background — primer dimer,
    mis-amplification, or off-target product would all present as pairs that fail
    to overlap at the expected insert size.

    A **read-overlap QC step before the pipeline proper** is the right way to
    catch this: a per-sample heatmap of observed overlap/fragment length, where
    the dense region should coincide with the expected amplicon size. Samples
    whose mass sits away from that size carry background rather than target. This
    run had no such screen available, so **the possibility that these libraries
    are simply poor cannot be excluded.**

    This does **not** undermine the errfun A/B — both arms received identical
    reads, and the merge rates differ by 1.1 points. But it does mean the
    *absolute* ASV counts, chimera fractions, and read-retention figures on this
    page should not be read as characteristic of soil, or of NovaSeq 6000, or of
    binned data generally.

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
largely *not* the ones it calls chimeric. Post-chimera Jaccard (0.433) is *lower*
than pre-chimera (0.545).

On PacBio, chimera removal reduced a +0.21% pre-chimera effect to −0.13%. Here it
takes a −24.1% difference to −26.1%.

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
   16S V4.
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
6. **A read-overlap QC screen ahead of the pipeline.** The 44% merge rate should
   have been caught and characterised *before* denoising, not noticed afterwards.
   A per-sample overlap/fragment-length heatmap against the expected amplicon
   size separates "these libraries carry background" from "this dataset is hard,"
   and would let a run like this be qualified up front. Worth having as tooling
   regardless of the binned-quality question, since it applies to every paired
   run — and it is the missing control that keeps this page's absolute numbers
   uninterpretable.

[PRJNA1504839]: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1504839
[#98]: https://github.com/HPCBio/dada2-rs/issues/98
[#44]: https://github.com/HPCBio/dada2-rs/issues/44
