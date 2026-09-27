# OMEGA_A is well calibrated on PacBio HiFi

**Verdict:** `OMEGA_A = 1e-40`, inherited from DADA2's Illumina calibration and
applied platform-blind, is **not too conservative on PacBio HiFi**. Both arms of
a truth-set probe on a 96-sample pooled ATCC MSA-1003 run agree, and they fail
in opposite directions, which is what makes the result usable:

- **True-positive arm** — 43 of 52 truth alleles recovered. The one additional
  call that relaxing `OMEGA_A` would rescue needs roughly a **25-order**
  relaxation to clear the threshold. The remaining 8 are not threshold misses at
  all: they are resolution limits or PCR/detection dropout, including a
  confirmed *Bifidobacterium* primer bias.
- **False-positive arm** — 263 FPs, and **all of them sit far below `omega_a`**,
  median `p_a` around **1e-140**. They are the chimera tail, not marginal calls
  admitted by a permissive threshold.

So the threshold is not the lever it was hypothesised to be. Tightening it would
not remove the false positives, because they are nowhere near it; loosening it
would not recover the missing alleles, because they are nowhere near it either.

## Why this needed a truth set

This question was **parked for exactly the right reason**, and the reason is
worth preserving. The PacBio churn that motivated it — ASVs that move when the
k-mer distance cutoff changes — concentrates in the low-abundance band, roughly
4 to 76 reads, which is precisely the `OMEGA_A` abundance-p-value borderline.
That correlation is suggestive and it is *not* evidence.

A concordance A/B against a reference run measures **whether an ASV moved**, not
**whether it was correct**. Every arm of a cutoff sweep can disagree with every
other and none of them is thereby wrong. Asking whether `OMEGA_A` is
mis-calibrated is a question about correctness, so it was unanswerable until
there was something to be correct *against*.

The gating resource has since landed: the ATCC MSA-1003 HiFi 16S truth set,
allele-level (52 alleles, not 20 strains), plus the `reference-eval` subcommand
that scores ASVs against truth with a p-value join
([issue 91](https://github.com/HPCBio/dada2-rs/issues/91)).

## Why both arms were needed

A TP arm alone would have been unfalsifiable in the direction that matters. "We
recover 43 of 52" is compatible with a well-calibrated threshold *and* with an
over-conservative one; the difference is whether the 9 misses are near the
boundary. Joining each miss to its abundance p-value is what turns the count
into a verdict — and it showed a ~25-order gap, not a near miss.

Symmetrically, an FP arm alone would have said only that we make 263 false
calls, with no way to tell an over-permissive threshold from a chimera problem.
The median `p_a` of 1e-140 settles it: these clear the threshold by 100 orders
of magnitude. No setting of `OMEGA_A` that keeps the true positives would
exclude them.

The two arms are independent and could have disagreed. They did not.

## What this dictates

- **Do not tune `OMEGA_A` for HiFi.** The hypothesis that DADA2's
  Illumina-calibrated default is wrong for long reads is, on this mock, refuted.
  Treat a proposal to change it as needing new evidence, not as an open question.
- **The low-abundance churn needs a different explanation.** The 4–76 read band
  where cutoff changes churn ASVs is still real; `OMEGA_A` is simply not the
  mechanism. See
  [KDIST cutoff decoupling](kdist-cutoff-decoupling.md).
- **The 263 FPs are a chimera question.** They are the tail that de-novo chimera
  removal does not catch, which points at
  [issue 55](https://github.com/HPCBio/dada2-rs/issues/55) (higher-order and
  trimera screening), not at a denoising threshold.
- **`OMEGA_C` is not covered by this**, and cannot be by the same instrument.
  The probe addressed `OMEGA_A`, which decides whether an ASV is *born* and so
  shows up in truth-set scoring. `OMEGA_C` decides read *attribution* after
  clustering, leaving the per-orientation ASV set unchanged — an ASV-level TP/FP
  arm is largely blind to it. The pair was parked together and only one half has
  been answered; the other is
  [issue 228](https://github.com/HPCBio/dada2-rs/issues/228).
- **This verdict is mock-specific and platform-specific.** It is one community,
  one region, one primer pair, on PacBio HiFi. An Illumina truth set on the
  *same* mock is needed before any cross-platform statement, because an Illumina
  subregion collapses alleles that HiFi resolves — see
  [issue 166](https://github.com/HPCBio/dada2-rs/issues/166).

## Provenance

Recorded from the probe reported on
[issue 93](https://github.com/HPCBio/dada2-rs/issues/93) and
[issue 166](https://github.com/HPCBio/dada2-rs/issues/166); tooling from
[issue 91](https://github.com/HPCBio/dada2-rs/issues/91) (`reference-eval`).
Truth set and its constraints:
[issue 166](https://github.com/HPCBio/dada2-rs/issues/166), the umbrella for
mock-community truth-set work.
