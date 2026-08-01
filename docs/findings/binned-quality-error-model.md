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

**Chimera removal is likewise unaffected:**

| | full (`pacbio`) | binned (`binned-qual`) |
|---|---|---|
| ASVs post-chimera | 45,455 → 23,022 (49.4% chim) | 45,598 → 23,174 (49.2%) |
| Reads retained | 95.12% | 95.15% |
| Non-chimeric ASV Jaccard | — | 0.9701 |

On the 44,798 ASVs present in *both* pre-chimera sets, the two arms agree on
**99.833%** of chimera calls — 21,975 agree chimeric, 22,748 agree clean, 75
disagreements total. That reproduces the mock's 114/114 agreement at roughly 600×
the scale. Across both datasets, **chimera calling is independent of the error
model**, so a coarser error model does not shift what gets removed.

Note ~49% of ASVs are chimeric here while carrying under 5% of reads — the same
large-but-low-mass tail seen in the mock. A counterintuitive corollary: chimeric
ASVs have a *higher* median abundance (28) than retained ones (22), because chimera
formation requires abundant parents. Abundance thresholds are therefore a poor
proxy for chimera filtering on data like this.

## Choosing the wrong errfun for binned input — and how chimera removal hides it

A user who reaches for `--errfun pacbio` on binned PacBio data does not get an
error, a warning, or the PacBio model. `pacbio_errfun` tests whether the last
quality column is Q93; on binned input it is not, so the function silently falls
back to plain LOESS. R's `PacBioErrfun` at least emits a `message()`; dada2-rs
prints nothing. The run's artifacts are indistinguishable from a legitimate
`pacbio` run.

Tested end to end on the real dataset, same reads and same pinned `learn-errors`
sample, varying only the errfun. **The full chain is worth reading in order,
because the endpoint alone is misleading:**

1. **The fitted model is wrong, measurably.** Mass-weighted mean error rate
   4.10e-4 against 4.99e-4 for `binned-qual` — **17.8% low**, and 0.80× the
   full-quality reference fit from the same sample.
2. **The error is concentrated in a spurious trough.** Across Q31–Q39 LOESS
   underestimates by up to **16×** (Q35: 8.10e-5 vs 1.30e-3). That region is not
   empty — it holds ~8.4% of observation mass, the off-bin values created by derep
   quality-averaging. There are real anchors at Q27/Q35/Q40 with sparse support
   between them, and the local fit dives into the gap.
3. **Inference shifts in the predicted direction.** A low error model makes
   variants harder to explain as errors, so the mismatched arm **over-splits:
   +96 ASVs (+0.21%)**, 45,694 vs 45,598. Of the 479 arm-specific ASVs, **98%
   appear in exactly one sample** — the rare, near-threshold population a slightly
   low model would tip. They are spread across the run (median 2 per sample; 61
   samples have none) rather than confined to a few pathological samples.
4. **Chimera removal then conceals it.** 59% of those extra ASVs are chimeric.
   Post-chimera the mismatched arm has **31 *fewer*** ASVs (−0.13%) — the effect
   reverses sign and rounds away. All three arms land within 0.7% of each other,
   with read retention identical to two decimals.

| arm | pre-chimera | post-chimera | chimeric | reads retained |
|---|---|---|---|---|
| full quality, `pacbio` (reference) | 45,455 | 23,022 | 49.4% | 95.12% |
| binned, `binned-qual` (correct) | 45,598 | 23,174 | 49.2% | 95.15% |
| binned, `pacbio` (**mismatched**) | 45,694 | 23,143 | 49.4% | 95.16% |

**Read step 4 as concealment, not correction.** Chimera removal discarding those
ASVs does not establish that they were chimeric, nor that the error model handled
them correctly — two errors happened to cancel. The cancellation is contingent on
things this finding does not control:

- It requires chimera removal to be run at all, at settings that catch them.
- It is chemistry-dependent. ~91% of HiFi observation mass sits in the top bin,
  where the mismatched model is only 8% low, so the badly-wrong interpolated region
  barely reaches the likelihood computation. **That offset is a property of HiFi,
  not of binning** (see *Path forward*).
- Nothing informs the user any of it happened.

Why the mass concentration matters is worth stating plainly: the mismatch is large
where the data is sparse and small where the data is dense, so on this platform it
mostly cannot express itself. On a platform whose quality mass is spread across its
bins, the same fitting error would land where the reads actually are.

## What this dictates

1. **`binned-qual` is the right default for natively binned PacBio data.** It does
   not degrade recovery even when the fit is driven by essentially one quality
   value, and on Revio it reproduces full-quality inference almost exactly.
2. **Expect the magnitude to be chemistry-dependent, and do not quote the 4×.**
   Binning cleans up a residual error-variant tail only where one exists. The large
   SequelIIe effect says more about that platform's full-quality error model than
   about binning. On Revio, expect a near-no-op.
3. **Report both endpoints — pre- and post-chimera — not just the one that
   matters.** Post-chimera answers "does this change my results?"; pre-chimera is
   where the mechanism is *visible*. Three times in this investigation a real
   pre-chimera difference vanished downstream (the mock ASV counts, an R
   comparison, and the errfun mismatch above). Collapsing to the endpoint would
   have reported three null results and hidden three real effects, each of which
   could resurface under a different downstream configuration or platform. Raw
   `n_asv` on its own remains a trap in both directions.
4. **Do not use `--errfun pacbio` on binned input, and do not rely on chimera
   removal to cover it if you do.** The fallback is silent, the resulting model is
   ~18% low with a 16× trough, and the fact that the consequence washes out here is
   a HiFi-specific accident (see the section above). A warning belongs in the tool —
   tracked in [#98], with the detection available exactly from `count == 1`
   uniques, whose `qual_sum` is the raw unaveraged quality vector.
5. **A community report that `binned-qual` inflates ASV counts is not reproduced
   here** — binning fixed reads either cut ASVs 4× (SequelIIe) or left them within
   0.3% (Revio). That report may have been pre-chimera counts, a different errfun
   implementation, or Illumina data. The clean test is an errfun sweep on fixed
   reads (below), which is a different experiment from any arm run here.
6. **Chimera removal can delete a real allele, independent of the error model.** In
   the Revio mock arm, `NC_004461.1:1722239-1723890(-)` (*S. epidermidis* ATCC 12228) is
   an exact truth match at 691 reads (rank 161/1345, p=1.8e-62) yet is flagged
   chimeric at **every** min-fold from 1.5 to 32. No threshold rescues it; raising
   min-fold makes chimera calls harder, so its putative parents clear even a 32×
   bar. This is the expected multi-allele 16S failure mode — alleles of one genome
   decomposing into each other — and it is a chimera-removal false positive, not a
   denoising defect.

## Path forward

Two datasets, one platform family, one of them without truth. Remaining gaps:

1. **Binned Illumina — where the offset above should not hold.** The reason the
   errfun mismatch washes out on HiFi is that ~91% of observation mass sits in one
   bin, so a fit that is badly wrong between bins lands where almost no data is.
   Illumina does not have that structure: quality declines along the read, so i100's
   `{12, 24, 38}` and NovaSeq's ~4 levels each carry substantial mass, much more of
   it in the interpolated region. Two things compound it — with only 3 anchors stock
   LOESS fits an exact quadratic through three points (the oscillation case in
   [#98]), and the #95 graceful-degradation guard triggers at ≤ 5 populated columns,
   so on Illumina it may actually fire, putting the fallback in genuinely degenerate
   territory rather than merely mis-shaped. **Prediction: the mismatched-errfun
   penalty is larger on binned Illumina than on binned PacBio.** The PacBio null
   should not be read as reassurance here.
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
