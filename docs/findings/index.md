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
- [Minimizers as the pre-alignment screen](minimizer-screening.md) — **negative
  result so far, and a page about being wrong twice.** On the 20-sample MiSeq SOP
  a winnowed-minimizer screen does **not** reproduce the k-mer screen's table: it
  fragments clusters into spurious low-abundance ASVs, one at Hamming-1 from a
  3,066-read parent. Calibrating its cutoff cuts the disagreement from 17 churned
  ASVs to 2 (L1 1.10% → 0.34%) but not to zero — against a control whose noise
  floor is **exactly** zero. Falsifies two of its own design claims: that a
  cutoff transfers between metrics sharing an algebraic form (it does not — 0.42
  passes 27.6% of pairs on one metric and 9.0% on the other, matching at ≈0.64),
  and that k=11 was a sane default (k=8 is). Also: the earlier revision of the
  page concluded the *opposite* from 2-sample fixtures, so it is now the
  reference case for **why a concordance fixture is only evidence about the part
  of the distribution it contains** — an 8-ASV table has no low-abundance tail
  for a screen to fragment. Independently: PacBio at k=5 measures at **100.00%**
  pass, a literal no-op; and ASV *ordering* is nondeterministic run-to-run on the
  existing default backend.
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
- [Inside `b_compare`: the store scan](compare-store-scan.md) — the serial store
  was never doing serial work; it was walking a 160-byte `Raw` to read 8 bytes of
  it. The tell was a **per-unit cost that did not move with working-set size**
  (18.4–19.6 ns/raw-visit across a 5× range in pool size), which is a cache-line
  signature, not a residency one. A microbenchmark put the strided read at
  **83–87% of the scan**; hoisting `e_minmax` into a dense array parallel to
  `raw_cluster` cut the store **71%** and `run_dada` **20.1–21.3%** on soil 16S
  (−15.2% on ITS2), byte-identical, with the untouched phases flat. Effective
  cores go 31.4 → 38.4 of 64. Then two follow-on levers were built and
  **neither merged**: removing the map's per-call allocation collapses `sys` by
  94% and returns *nothing* — wall +2.5%, `user` +16%, total CPU up — because
  that `sys` time was never on the critical path; and map load imbalance, the
  suspected cause of the 86–90% parallel efficiency, measures at ±0.1% on both
  pool shapes, so the residual is unmeasured work rather than idle threads.
  Ends with the three instrument errors behind those nulls — a hypothesis killed
  by arithmetic instead of node time, a skew model that could not exhibit skew,
  and the right estimator applied to the one arm whose variance was intrinsic —
  and the rule they share: **a null is evidence only once the instrument has
  been shown capable of returning something else.**
- [Thread scaling & memory placement](thread-scaling-and-placement.md) — **two
  results with opposite characters.** How many threads to use is strongly
  data-dependent: the best count differs *two-fold between two reads of the same
  pool*, and is predicted by the map's screen/align split, which the verbose log
  already prints (50.8% screen scales to 96 threads, 67.7% peaks at 64, 83–86%
  peaks at 48 and degrades past it). Where to put the memory is **not**
  data-dependent — binding one job to one NUMA domain is worth **25–29% on every
  arm**, across a 35-point gap in screen share, and gives 2.6–2.9× throughput
  when two jobs each take a domain. Also the page that reconciles this with
  [Measuring on a NUMA node](measuring-on-numa.md), which concluded the opposite:
  that page's own phase table **already showed the serial phases 23–37% faster
  bound**, matching what we just measured; the map penalty at 64 threads
  outweighed it and the total said "binding loses". Ends with four instrument
  errors, including one that inverted the sign of a reported trend.
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
