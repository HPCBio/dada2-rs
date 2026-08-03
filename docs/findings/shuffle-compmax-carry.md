# Carrying `compmax` across buds: a workload-dependent optimization

**Verdict:** carrying the best-cluster map (`compmax`) across bud rounds — the
remaining serial lever in pooled `dada` after [#86][86] — is **conditionally
worth it, and the condition is the workload, not the platform label**. On a
NovaSeq soil 16S run it projects to **−6.4% to −7.5% wall** on R1 and only
**−1.7% to −2.3%** on R2; on ITS2 from the *same instrument and run* it projects
to a **+10.4% to +10.5% wall penalty**. Because the sign of the result flips
between two amplicons of one sequencing run, any implementation must be
**optional and gated, off by default**. Given a modest and uncertain upside
against a nontrivial correctness risk, this is filed as **low priority** — but
explicitly **not ruled out**, since prior optimizations in this codebase have
outperformed their projections.

The error model (binned vs. default quality) has **no meaningful effect** on this
result — the deciding statistic shifts by at most 0.02 across all four arms.

## Background

Pooled `dada` alternates two phases. `b_compare` is parallel; `b_shuffle` is
serial, and after [#86][86] made the shuffle incremental it became the dominant
serial cost — the Amdahl ceiling on pooled runs. Within `b_shuffle_converge`
the work splits in two:

- **build** — construct `compmax` (every raw's best cluster) by walking the
  per-cluster comp vecs. A contiguous, cluster-major scan. Paid **fresh on every
  call**, i.e. once per bud round.
- **reconcile** — after a move pass, recompute the raws whose best cluster may
  have changed, via the raw-major inverted index. **Scattered** reads.

[Issue #87][87] proposes carrying `compmax` across buds to remove the repeated
build. The question is what that actually saves.

## The cost model

The proposal does not *delete* the build's work so much as **relocate** it. A
bud adds a cluster and steals raws from many others, changing their `reads` and
so marking their comps affected — pushing that work through the reconcile path
instead. So the change saves `build_time` and pays

```
f × comps_build × reconcile_ns_per_comp
```

where `f` is the post-bud reconcile volume as a fraction of the build volume.
It is a net win only when

```
f  <  build_ns_per_comp ÷ reconcile_ns_per_comp
```

This matters because the two phases have **different access patterns**, and in
this codebase access pattern rather than op count is the binding cost (four
prior scattered-access negative results). Comparison *counts* cannot answer the
question; only ns/comp can.

Instrumentation for this lives in `ShuffleStats` ([#122][122], [#123][123]) and
prints under `--verbose`.

## Setup

One NovaSeq 6000 soil run, four arms — two amplicons × two error models — full
pooled mode, release build, per read direction (R1/R2), 8 logs total:

| Arm | Amplicon | Error model |
| --- | --- | --- |
| `full-pooling` | 16S | default |
| `full-pooling-binned` | 16S | binned quality |
| `full-pooling-ITS` | ITS2 | default |
| `full-pooling-ITS-binned` | ITS2 | binned quality |

Extraction: `dev/compare_shuffle_87.py <run_dirs...>`.

## Results

| Run | Read | ratio | build % *comps* | build % *time* | `f` | break-even | verdict | proj % wall |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16S default | R1 | 1.96× | 71% | 41% | 0.34 | 0.51 | **WIN** | **−7.5%** |
| 16S default | R2 | 1.85× | 66% | 36% | 0.49 | 0.54 | WIN | −2.3% |
| 16S binned | R1 | 1.98× | 70% | 40% | 0.36 | 0.50 | **WIN** | **−6.4%** |
| 16S binned | R2 | 1.87× | 65% | 37% | 0.50 | 0.53 | WIN | −1.7% |
| ITS2 default | R1 | 2.83× | 61% | 27% | 0.65 | 0.35 | **LOSS** | **+10.5%** |
| ITS2 default | R2 | 2.64× | 59% | 27% | 0.69 | 0.38 | LOSS | +10.5% |
| ITS2 binned | R1 | 2.64× | 60% | 27% | 0.67 | 0.38 | **LOSS** | +10.4% |
| ITS2 binned | R2 | 2.51× | 59% | 27% | 0.71 | 0.40 | LOSS | +10.4% |

### 1. Comparison counts overstate the prize by ~2×

Reconcile comparisons cost **1.85–2.83× more than build comparisons** (build
4.3–6.3 ns/comp; reconcile 11.4–13.4 ns/comp). The build is 59–71% of
comparisons but only **27–41% of time**. Sizing this work from comp counts
alone would have overstated #87's ceiling roughly twofold — a direct
re-confirmation of the access-pattern cost model.

### 2. The verdict flips between amplicons of the same run

16S `f` is 0.34–0.50; ITS2 `f` is 0.65–0.71. **No overlap.** Both terms of the
inequality move the same direction: ITS2 also has a *lower* break-even
(0.35–0.40 vs 0.50–0.54) because its scattered reconcile is relatively more
expensive. Fewer clusters means each bud perturbs a larger share of total comp
volume, so `f` rises as cluster count falls. `comps/build` separates the two
regimes ~4× (16S 8.8–10.0e6, ITS2 2.3–2.8e6) and is the natural gate.

This is the load-bearing result: **16S and ITS2 here are the same instrument,
the same run, the same pipeline.** The sign of the optimization is set by the
workload's cluster structure, not by anything a user would think to configure.

### 3. The error model does not matter

Largest binned-vs-default shift in `f` is **0.02** (16S R1 0.34→0.36; ITS2 R2
0.69→0.71), well inside the R1/R2 spread within either arm. Verdicts, ratios and
projections are unchanged. This finding needs **no error-model caveat**.

### 4. Averaged reconcile volume is not a valid proxy for `f`

`f` was first inferred from the *mean* reconcile volume, giving 0.16–0.24 (16S)
and 0.49–0.52 (ITS2). Measuring the **post-bud** reconcile directly ([#123][123])
gave 0.34–0.50 and 0.65–0.71 — the proxy understated `f` by ~35–110%, because a
bud perturbs more than a mid-loop iteration does. The proxy would have called
16S a comfortable win (−12 to −14% wall) rather than the marginal one it is.
**Measure the post-bud reconcile; do not average.**

## What this dictates

- **Do not implement #87 as an unconditional optimization.** It is a ~10% wall
  *regression* on ITS2-like workloads.
- **If implemented, gate it on cluster structure** (`comps/build` or cluster
  count) and ship it **off by default**. The 4× separation gives a clean
  threshold, but the gate must be validated on more than these two workloads
  before it could ever default on.
- **Expect less than the projection.** The model assumes the build is removed
  entirely and the only new cost is the relocated reconcile. A real
  implementation also pays bookkeeping to keep `compmax` valid across buds and
  must still scan the newly-budded cluster's comps. That unmodelled overhead
  plausibly erases the R2 margin outright (`f` 0.49–0.50 vs break-even
  0.53–0.54 — a gap of 0.03–0.05).
- **Correctness is the real cost.** [#86][86] achieved byte-identical output by
  keeping `compmax` exactly equal to a full rebuild at current reads. Carrying it
  across buds is a strictly harder invariant with more stale-max hazards. Any
  attempt must clear the ASV-level concordance guardrail, not just wall time.
- **Priority: low, but not closed.** Realistic upside is ~5–7% wall on
  16S-R1-like workloads and ~0 on R2 — materially smaller than #86's −16%, and
  workload-conditional. Kept open because optimizations here have surprised us
  before (the k-mer memory work in [#32][32]/[#43][43] most notably), so a future
  evaluation is warranted rather than a closed door.

## Reproducing

Run pooled `dada` with `--verbose` and collect:

```
[dada] shuffle scan split: build=… (…% of scanned) over … builds …
[dada] shuffle scan time: build=…s (… ns/comp)  reconcile=…s (… ns/comp)  …
[dada] #87 projection: post-bud reconcile=… over … buds, f=… vs break-even … → WIN|LOSS …
```

Then `dev/compare_shuffle_87.py <run_dirs...> [--csv out.csv]` for the table.

[86]: https://github.com/HPCBio/dada2-rs/pull/86
[87]: https://github.com/HPCBio/dada2-rs/issues/87
[122]: https://github.com/HPCBio/dada2-rs/pull/122
[123]: https://github.com/HPCBio/dada2-rs/pull/123
[32]: https://github.com/HPCBio/dada2-rs/issues/32
[43]: https://github.com/HPCBio/dada2-rs/issues/43
