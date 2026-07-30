# Binned quality scores & the `binned-qual` error model

**Verdict (PacBio, mock community):** collapsing PacBio quality scores to 7 binned
levels and fitting with `--errfun binned-qual` **costs nothing in reference
recovery and sharply reduces the residual error-variant tail**. On the *same reads*,
binning cut ASVs 4× and non-chimeric near-variants 12×, while recovering the
**identical** set of 43/52 truth alleles. The effect is **denoising, not filtering**.

Scope: this is PacBio HiFi (SequelIIe + Revio Kinnex 16S) on a single mock
community (ATCC MSA-1003), where truth is known. It is **not** yet evidence for
Illumina binned chemistries (i100, NovaSeq), nor for real communities where
"absorbed a near-variant" cannot be distinguished from "lost a real low-abundance
strain." See *Path forward*.

## Setup

Three arms over the same MSA-1003 libraries and barcodes (96 samples), all with
k=7, band 32, `OMEGA_A` 1e-40, NW backend, `dada-pooled`:

| Arm | Sequencer | Quality | `--errfun` | nq |
|---|---|---|---|---|
| **A** | SequelIIe | full range (Q2–Q93) | `pacbio` | 94 |
| **B** | Revio | natively binned, 7 levels | `binned-qual` | 41 |
| **C** | SequelIIe, re-binned | 7 levels, from arm A's reads | `binned-qual` | 41 |

Arm C is the controlled one: `dev/bin_fastq_quals.py` re-bins arm A's raw FASTQs to
the exact level set Revio emits — `{3, 10, 17, 22, 27, 35, 40}` — so the *only*
varying term is the quality representation. Binning must be applied to raw FASTQs
**before derep**: dereplication averages qualities across the reads collapsed into a
unique, so `bin→derep` and `derep→bin` disagree, and native binning is `bin→derep`.

Evaluation is truth-anchored via `reference-eval` against the 52-allele MSA-1003
set, joined to `remove-bimera-denovo` output (`consensus`,
`--min-fold-parent-over-abundance 3.5`) so every count is post-chimera-removal.

## The error model is fit almost entirely from one quality value

In the natively binned Revio arm, **98.9% of all error-model observation mass sits
at Q40 alone**; the 7 nominal bins together hold 99.06%. Re-binned SequelIIe is
more extreme still (99.74% at Q40), because ~99.5% of its mass sits at Q≥50 across
~54 distinct values that all collapse into the top bin.

A wrinkle worth knowing before reading a `trans` matrix: the binned Revio model has
**32 populated Q columns despite the FASTQ carrying only 7 distinct values**. That
is derep quality-averaging partially *de-binning* the data, and it accounts for the
0.94% of off-bin mass. It is not a failure of binned-quality detection.

## Evidence — arm C vs arm A (same reads, only the quality representation differs)

| Metric | A: unbinned (`pacbio`) | C: binned (`binned-qual`) |
|---|---|---|
| ASVs | 1464 | **369** (−4×) |
| Truth alleles recovered | 43/52 | **43/52 — identical set** |
| Non-chimeric "near" ASVs (1–3 diffs from truth) | 1148 | **93** (−12.3×) |
| ...per 100k reads | 54.9 | **4.5** |
| Non-chimeric FPs | 54 | **21** |

Recovery is identical allele-for-allele, not merely equal in count.

**It is denoising, not filtering.** Re-prepping under binned qualities changes maxEE
and drops 1.7% of uniques (arm C's input is a strict *subset* of arm A's: 351497 of
357653, zero arm-C-only). That is not the cause: of arm A's 1148 non-chimeric near
ASVs, 1060 are no longer ASVs in arm C — and **100.0% of those are still present as
input uniques** (median abundance 66, max 373). The error model absorbed them into
their parents; the filter-dropped uniques contribute nothing.

## Cross-sequencer arm (B vs A) — consistent, and now de-confounded

Comparing native Revio to SequelIIe gave the same directional result (9.7 vs 54.7
non-chimeric near-variants per 100k reads) but confounded quality representation
with sequencer chemistry and 3× depth. Arm C resolves it: **the effect is the
quality representation.** Re-binned SequelIIe (4.5/100k) is in fact cleaner than
native Revio (9.7/100k).

Two artifacts of that cross-sequencer comparison are worth recording because they
are easy to misread:

- **Chimera calling is unaffected by the error model.** Controlled for ASV set, on
  the 626 shared ASVs both arms call **114 chimeric, with zero disagreements**. The
  raw rate gap (45.8% Revio vs 15.0% SequelIIe) is entirely carried by run-specific
  ASVs (69.8% vs 12.5% chimeric) — deeper sequencing of the *same* libraries
  surfaces more rare library chimeras.
- **`--min-fold-parent-over-abundance` is a near-dead lever here.** Sweeping
  1.5 → 32 moves the chimeric fraction under 1% on both arms (Revio 46.4%→45.7%,
  SequelIIe 15.0%→14.8%).

## What this dictates

1. **`binned-qual` is the right default for natively binned PacBio data.** It does
   not degrade recovery even when the fit is driven by essentially one quality
   value, and it materially reduces the spurious-variant tail.
2. **Never evaluate this on raw ASV counts.** Chimera removal does most of the work
   and removes the excess preferentially — the binned Revio arm's extra ASVs were
   69.8% chimeric. Use post-chimera non-chimeric near-variants per read plus truth
   recovery. `n_asv` is a trap here twice over.
3. **A community report that `binned-qual` inflates ASV counts is not reproduced
   here** — on fixed reads, binning *cut* ASVs 4×. That report may have been
   pre-chimera counts, a different errfun implementation, or Illumina data. The
   clean test is an errfun sweep on fixed reads (below), which is a different
   experiment from any arm run here.
4. **Chimera removal can delete a real allele, independent of the error model.** In
   the Revio arm, `NC_004461.1:1722239-1723890(-)` (*S. epidermidis* ATCC 12228) is
   an exact truth match at 691 reads (rank 161/1345, p=1.8e-62) yet is flagged
   chimeric at **every** min-fold from 1.5 to 32. No threshold rescues it; raising
   min-fold makes chimera calls harder, so its putative parents clear even a 32×
   bar. This is the expected multi-allele 16S failure mode — alleles of one genome
   decomposing into each other — and it is a chimera-removal false positive, not a
   denoising defect.

## Path forward

The result is strong but narrow: one mock community, one platform, truth known. Two
gaps, in priority order.

1. **A real community, full-quality Revio run** — the nodule-based 16S dataset
   behind the trimera question. Being full-quality, it supports the same controlled
   `pacbio` vs `binned-qual` A/B as arm C. Without a truth set the endpoint must
   change — no TP/near/FP — so use ASV set-identity, abundance-weighted read
   retention, and the pre/post-chimera split. This is the test that distinguishes
   "absorbs noise" from "under-splits genuine low-abundance variation": if binning
   removes ASVs that are *abundant and non-chimeric* in a real community, that is
   under-splitting, not denoising. The trimera question makes it doubly
   informative, since higher-order chimeras are exactly what a coarser error model
   might reclassify.
2. **An errfun sweep on fixed reads** — `binned-qual` vs `loess` vs `pacbio` vs R's
   `loessErrfun_mod4` (via `--errfun external`), all on the same derep inputs. This
   isolates the errfun term exactly and is the only clean test of the ASV-inflation
   report. Cheaper than a new sequencing arm: it needs only `learn-errors` + `dada`
   re-runs.

Illumina binned chemistries (i100 `{12,24,38}`, NovaSeq ~4 bins) remain untested and
should not inherit this verdict — see also the binned-quality wildcard in
[Band size & platform defaults](band-size-platform-defaults.md).
