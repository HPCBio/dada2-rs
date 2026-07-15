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
- [Band size & platform-aware defaults](band-size-platform-defaults.md) — the
  16/32 Illumina/HiFi band default is vindicated; the two platforms fail
  band-tightening through opposite mechanisms, so a single global band would be
  wrong.
