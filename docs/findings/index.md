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
- [Carrying `compmax` across buds](shuffle-compmax-carry.md) — the remaining
  serial lever in pooled `dada` is worth −7.5% wall on 16S and a +10.5%
  *regression* on ITS2 from the same NovaSeq run; comparison counts overstate
  the prize ~2× because the scattered reconcile costs ~2× per comparison. Any
  implementation must be gated and off by default. #87 closed on this evidence;
  the phase itself is reconsidered as a whole in #124.
- [Redesigning `b_shuffle`](shuffle-build-scan.md) — **negative result, and the
  phase is closed.** Two build-scan redesigns were falsified: data layout, then
  an exact bounded scan that pruned 56% of comparisons and *still* cost +4.8%
  wall because the survivors ran at the scattered rate (6.83 → 14.09 ns/comp).
  Completing the phase accounting then removed the reason to keep looking —
  `b_shuffle` is 5.6–27% of `run_dada`, not the 46–55% of pooled wall that
  motivated #124, and no remaining phase is worth more than 3%. The lesson that
  generalises beyond the arithmetic: re-measure an optimisation's *premise*, not
  only its design. Two of its three closures were later reopened and built —
  the move pass in #132 (−4.3 to −5.6% of `run_dada`) and the reconcile in #136
  (**−5.2 to −28.3%**, the top of that range on a diverse soil pool) — both by
  finding a **cluster-major route** to information the page had priced at the
  scattered rate. "Unreachable" turned out to be a
  claim about a route, not a quantity. Then the third fell too: #139 deleted the
  build scan rather than redesigning it, carrying `compmax` across buds for
  **−28.8 to −32.0% of `run_dada`** on soil 16S and reversing #87's +10.5% ITS2
  regression into a −8.1 to −13.7% win — because #136 had changed the price of
  the work #87 was costed against. A closure can expire when a *dependency*
  changes, not only when the idea improves. Also the page on why exactness is worth enforcing on a change
  you throw away — it caught a latent tie-break bug and is what makes the verdict
  final rather than ambiguous.
- [Inside `b_compare`: screen vs align](compare-screen-vs-align.md) — the phase
  that is 63–88% of `run_dada`, split by time for the first time. **Both
  platforms are align-dominated**, PacBio most of all, refuting the prediction
  that k=7's 16× larger k-mer vectors would make HiFi screen-bound: the screen
  does cost 3.1× more per comparison there, but the alignment it guards costs
  7.7× more. The k-mer screen closes as a perf target — it returns 3–7× on what
  it costs, and a perfect free replacement caps out at 15–23%. What is left is
  the banded NW **DP kernel at 30–63% of `run_dada`**, the largest single share
  measured in this project. It is execution-bound to ~48 threads and
  **memory-bandwidth-bound above that**, which the rolling-`d16` change exploits
  for 3.5–6% — a result that took three attempts to measure, see
  [Measuring on a NUMA node](measuring-on-numa.md). Ends where it did not expect
  to: the largest available win was **running wider** (48 threads beats 24 by
  28–41%, output-identical), and a 24/48/96/128 sweep on both platforms then
  found a **flat serial floor** — which flips the ranking. At 48 threads the serial block is
  over half of `run_dada` and `b_shuffle` is two-thirds of it, so
  [#124's closure](shuffle-build-scan.md) turns out to have been conditional on
  thread count. Includes a section on which of these numbers travel to other
  hardware and which do not.
- [Inside `b_compare`: the serial remainder](compare-serial-fold.md) — the 329 s
  (48% of `run_dada`) that was neither the parallel map nor the store. Attributed
  first, on the issue's own insistence that an optimisation not be designed
  against a share-of-wall figure before re-measuring what it is made of — and the
  attribution made the fix obvious: **six serial passes over the map's result
  vector**, folded into the store loop that already walked it. Seven passes become
  one, for **−24 to −30% of `run_dada`** on two soil pools, byte-identical. The
  cache-residency explanation for why those passes were expensive is **falsified**
  here — ITS2 reduces at the same ~29 ns/raw-visit as soil 16S with a fraction of
  the working set, so it was six passes, not where the data lived. Leaves the
  serial **store** as the whole remaining block, at a suspiciously flat
  18.4–19.6 ns/raw-visit across a 5× range in pool size.
- [Measuring on a NUMA node](measuring-on-numa.md) — **a methodology result that
  reversed a verdict.** The benchmark node has two NUMA domains and nothing was
  ever pinned, so page placement re-rolled every run and replicates of the *same
  binary* disagreed by up to 25% on the serial scattered phases. A real 3.5–6%
  DP improvement was twice measured as "flat, below the noise floor" and twice
  recommended for closure; with placement fixed it is unambiguous, replicates
  agreeing to 0.4 points. Also: the obvious fix (bind to one domain) is the
  wrong one — it slows the parallel map 21% — and `--interleave=all` is both
  more reproducible and closer to production. Includes which earlier results
  this puts in question.
- [Band size & platform-aware defaults](band-size-platform-defaults.md) — the
  16/32 Illumina/HiFi band default is vindicated; the two platforms fail
  band-tightening through opposite mechanisms, so a single global band would be
  wrong.
