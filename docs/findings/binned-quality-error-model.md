# Binned quality scores & the `binned-qual` error model

**Verdict (PacBio):** `--errfun binned-qual` on 7-level binned quality scores is
**safe to use — it never cost reference recovery, and what it changes depends on
the chemistry of the reads being binned, not on binning itself.**

- On **SequelIIe** reads (mock community), binning the *same reads* cut ASVs 4× and
  non-chimeric near-variants 12× while recovering the **identical** set of 43/52
  truth alleles. The effect was **denoising, not filtering**.
- On **Revio** reads (a real community, 237 samples), binning is close to a
  **no-op**: ASV set Jaccard 0.97, 99.9% of reads in shared ASVs, ASV count within
  0.3%. The learned error models agree to **2.1%** mass-weighted.

The unifying explanation: binning absorbs a residual error-variant tail *when one
exists*. SequelIIe's full-quality model carries such a tail; Revio's does not, so
there is nothing left to absorb. Since Revio is the platform that actually emits
binned data, the practical reading is that **working with binned qualities loses
essentially nothing**.

Scope: PacBio HiFi only, on two datasets — one mock community (ATCC MSA-1003,
SequelIIe + Revio, truth known) and one real community (a 237-sample PacBio Kinnex
16S library sequenced on a Revio in 2024). It is **not** evidence for Illumina
binned chemistries (i100, NovaSeq). See *Path forward*.

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

**Off-bin mass is much larger on real data: 7.7%, versus 0.94% on the mock** (Q40
alone holds 90.8% rather than 98.9%). Off-bin values arise from averaging qualities
across the reads collapsed into a unique, so they scale with how *abundant* the
dominant uniques are. A mock community is close to evenly covered; the real dataset
is dominated by a single taxon whose top ASV carries ~3.25M reads, producing far
heavier averaging. Expect real communities with dominant taxa to de-bin
substantially more than mocks.

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

## Real-community arm (Revio) — binning is close to a no-op

The mock arms above all bin *SequelIIe* reads. Repeating the controlled A/B on
**Revio** reads from a real community gives a very different magnitude — a
237-sample PacBio Kinnex 16S library sequenced on a Revio in 2024, denoised with
`dada-pseudo` (streaming), k=7, band 32. Both arms used the same reads and the same
pinned `learn-errors` sample, so only the quality representation differs.

**The learned error models agree to 2.1%.** Comparing each bin against the full
model averaged over that bin's *source* Q range, weighted by observation mass:
every bin falls within 0.06 log10, and the overall mass-weighted mean error rate is
5.097e-4 (full) vs 4.988e-4 (binned). Per-transition, all 12 substitution types land
in ratio 0.70–1.25 (median |log10 ratio| 0.047), so the aggregate agreement is not
masking compensating per-transition errors.

> **Comparison trap.** A naive same-Q comparison shows binned Q40 = 4.6e-4 against
> full-model Q40 = 1.29e-3 — 2.8× apart, and an artifact. In the full model Q40 is a
> *mediocre* quality (HiFi mass sits at Q83–93); in the binned model Q40 is the top
> bin absorbing everything ≥40. A bin must be compared against its source Q range,
> mass-weighted.

**Inference is nearly unchanged:**

| Metric | full (`pacbio`) | binned (`binned-qual`) |
|---|---|---|
| Distinct ASVs | 45,455 | 45,598 (+0.3%) |
| Median ASVs/sample | 690 | 694 |
| Total reads | 23,256,894 | 23,256,894 |
| ASV set Jaccard | — | **0.9685** |
| Reads in shared ASVs | 99.92% | 99.70% |

Nothing like the 4× collapse seen on SequelIIe. The reconciliation is that binning
absorbs a residual error-variant tail only where one exists: SequelIIe's
full-quality arm carried 54.7 non-chimeric near-variants per 100k reads, native
Revio only 9.7. Revio's full-quality model is already clean, so binning has little
left to do.

This also **retires the under-splitting concern** raised when the mock result stood
alone. If binning were over-absorbing genuine low-abundance variation, a diverse
real community would lose ASVs. Instead the binned arm has *more* run-specific ASVs
than the full arm (800 vs 657), at a higher median abundance (33 vs 22).

## What this dictates

1. **`binned-qual` is the right default for natively binned PacBio data.** It does
   not degrade recovery even when the fit is driven by essentially one quality
   value, and on Revio it reproduces full-quality inference almost exactly.
2. **Expect the magnitude to be chemistry-dependent, and do not quote the 4×.**
   Binning cleans up a residual error-variant tail only where one exists. The large
   SequelIIe effect says more about that platform's full-quality error model than
   about binning. On Revio, expect a near-no-op.
3. **Never evaluate this on raw ASV counts.** Chimera removal does most of the work
   and removes the excess preferentially — the binned Revio mock arm's extra ASVs
   were 69.8% chimeric. Use post-chimera non-chimeric near-variants per read plus
   truth recovery. `n_asv` is a trap here twice over.
4. **A community report that `binned-qual` inflates ASV counts is not reproduced
   here** — binning fixed reads either cut ASVs 4× (SequelIIe) or left them within
   0.3% (Revio). That report may have been pre-chimera counts, a different errfun
   implementation, or Illumina data. The clean test is an errfun sweep on fixed
   reads (below), which is a different experiment from any arm run here.
5. **Chimera removal can delete a real allele, independent of the error model.** In
   the Revio mock arm, `NC_004461.1:1722239-1723890(-)` (*S. epidermidis* ATCC 12228) is
   an exact truth match at 691 reads (rank 161/1345, p=1.8e-62) yet is flagged
   chimeric at **every** min-fold from 1.5 to 32. No threshold rescues it; raising
   min-fold makes chimera calls harder, so its putative parents clear even a 32×
   bar. This is the expected multi-allele 16S failure mode — alleles of one genome
   decomposing into each other — and it is a chimera-removal false positive, not a
   denoising defect.

## Path forward

Two datasets, one platform family, one of them without truth. Remaining gaps, in
priority order.

1. **Chimera removal on the binned real-community arm.** The full-quality arm ran
   45,455 → 23,022 ASVs (49.4% chimeric, 95.12% of reads retained); the binned arm's
   post-chimera table is still outstanding. Given the two ASV sets agree at Jaccard
   0.97 this is expected to be uneventful, but it closes the comparison and tests
   whether a coarser error model shifts chimera calling on real data — the mock said
   it does not (114/114 agreement on shared ASVs).
2. **An errfun sweep on fixed reads** ([#98]) — `binned-qual` vs `loess` vs
   `pacbio` vs R's `loessErrfun_mod4` (via `--errfun external`), all on the same
   derep inputs. Note the arms above vary binning *and* errfun together; a sweep on
   one fixed binned input isolates the errfun term exactly, and is the only clean
   test of the ASV-inflation report. Cheaper than a new sequencing arm: it needs
   only `learn-errors` + `dada` re-runs.

[#98]: https://github.com/HPCBio/dada2-rs/issues/98

Illumina binned chemistries (i100 `{12,24,38}`, NovaSeq ~4 bins) remain untested and
should not inherit this verdict — see also the binned-quality wildcard in
[Band size & platform defaults](band-size-platform-defaults.md).
