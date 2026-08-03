# Findings

This section is a running lab-notebook of **decided outcomes**: experiments we
have run, what they showed, and — most importantly — **what each result dictates
for how the algorithm and the roadmap should proceed**.

It is deliberately separate from [Performance → Results](../results.md), which is
the head-to-head benchmark *scoreboard* (tables of wall time and memory vs R
DADA2). A finding here is not "how fast are we"; it is "we tested a hypothesis,
here is the evidence, and here is the path it opens or closes."

## What belongs here

- **A hypothesis and a verdict.** Each page states what was tested, the A/B setup,
  the numbers, and a plain conclusion.
- **The consequence.** Every finding ends with what it means going forward — a
  recommended default, a lever to pursue, or a direction ruled out.
- **Negative results are first-class.** A carefully-disproven idea is as valuable
  as a win because it stops us re-litigating it. Pages say "do not retry" where
  that applies.
- **Correctness framing.** Where a change touches inference, ASV concordance is
  reported (via `dev/compare_asvs.py`) using set-identity + abundance, not raw
  `n_asv` counts — see the [concordance
  guardrail](../benchmarking.md#5-concordance-validation-tooling).

## Index

- [KDIST cutoff decoupling](kdist-cutoff-decoupling.md) — the learn-errors and
  dada k-mer-distance cutoffs can be set independently; the dada-stage speedup is
  safe to take, while pushing the same cutoff into learn-errors perturbs the error
  model and churns real-abundance ASVs.
- [K-mer screen size](kmer-size-screening.md) — `--kmer-size` has ~no effect on
  the final chimera-filtered table on either platform; it is a speed/memory knob,
  not an accuracy knob (k=5 Illumina, k=7 PacBio for speed).
- [Pseudo-pooling: priors, not a re-fitted error model](pseudo-pooling-priors-vs-error-model.md)
  — our priors-only `dada-pseudo` beats emulating R's between-round error-model
  re-fit on every axis (+3,118 reads recovered, 709 fewer ASVs, 72× fewer ASVs
  unexplained by priors); R emulation stays opt-in, and confirming R's actual
  behaviour is still the open gate.
- [**Binned quality scores**](binned-quality.md) — a series, because the answer is
  **dataset-dependent** and we cannot yet predict which datasets are sensitive:
    - [PacBio (SequelIIe & Revio)](binned-quality-error-model.md) — binning never
      cost reference recovery (identical 43/52 alleles); what it *changes* depends
      on chemistry (4× ASV cut on SequelIIe, near-no-op on Revio). Choosing the
      wrong errfun is real at the model level but washes out downstream.
    - [Illumina NovaSeq 6000 — soil 16S](binned-quality-illumina-novaseq.md) —
      the errfun mismatch reaches the final table, but as an **abundance** error
      rather than a richness one: ASV counts agree to +3.4% while Jaccard is
      0.706, abundance L1 is 22.2% on shared ASVs, and one arm retains 17% more
      reads. It would pass a review that checks only ASV counts.
    - [Illumina NovaSeq 6000 — soil ITS2](binned-quality-illumina-its2.md) — the
      same instrument at **a third the sensitivity**, and the page that refutes
      both proposed mechanisms: mass concentration predicts the wrong sign, and a
      depth ladder shows coverage-per-variant predicts the wrong direction.
      A negative result, and the most decision-relevant page in the series.
    - [Reading the prep before the result](reading-the-prep.md) — the NovaSeq
      analysis was first published, then retracted twice: untrimmed degenerate
      primers behind heterogeneity spacers inflated the table 3–4× and *reversed
      the direction* of the effect, with no warning anywhere in the pipeline.
      The checks that catch it, for any dataset you did not prepare yourself.
- [Band size & platform-aware defaults](band-size-platform-defaults.md) — the
  16/32 Illumina/HiFi band default is vindicated; the two platforms fail
  band-tightening through opposite mechanisms, so a single global band would be
  wrong.
