# `b_shuffle`: pruning the build scan, and the full phase accounting

**Verdict: NEGATIVE, and the phase is now closed.** Two redesigns of `b_shuffle`'s
build scan were falsified — data layout, then exact bound pruning, the latter
built to a byte-identical implementation that **pruned 56% of the comparisons and
still cost +4.8% wall**. Completing the phase accounting afterwards then removed
the reason to keep looking: **`b_shuffle` is 5.6–27% of `run_dada`, not the
46–55% of pooled wall that motivated [#124][124]**, and the remaining phases have
either no exploitable headroom or a prize under 3%. #124 closes as fully
assessed; the real target is `b_compare`, which is 70–92% of `run_dada`.

Three results, in the order they were established:

1. **Layout rework is a dead end** — the build is already within 1.8× of the same
   arithmetic with no gather at all.
2. **Bounded pruning is structurally unable to win** — the break-even is set by
   access pattern, not by how well the bound prunes.
3. **The phase is smaller than believed** — and the premise in #124's title was
   wrong by 2–8×.

A fourth thing is worth extracting and may outlast the rest: **the exactness
discipline paid for itself even though the change was thrown away.** See "What
the guardrail bought".

## Background

After [#86][86], `b_shuffle` was believed to be 46–55% of pooled wall and the
Amdahl ceiling on pooled runs. It splits into a **build** (cluster-major rebuild
of `compmax`, paid fresh each bud round), a **reconcile** (raw-major, scattered),
and a **move pass** — which, as it turned out, nobody had ever timed.

[#87][87] tried to remove the build by carrying `compmax` across buds and was
closed: it relocated work into the more expensive scattered path, and its sign
flipped between two amplicons of one run (`shuffle-compmax-carry.md`).

That left the phase itself as the target. #124's starting observation is that
during a shuffle **`lambda_ij` is fixed** — only the scalar `reads_j` changes —
so the phase is `argmax_j (lambda_ij × reads_j)` with one factor static,
evaluated billions of times to discover a few hundred moves per build.

## Result 1: layout rework is a dead end

Modelled in `dev/shuffle-model/build_model.rs`, at nraw 285k / 400 clusters:

| Variant | ns/comp | vs current |
| --- | ---: | ---: |
| A current (AoS 24 B comps, `emax` + `compmax` 24 B) | 1.65 | 1.00× |
| B fused 16 B per-raw map | 1.84 | 0.89× |
| C fused map + SoA 12 B comp streams | 1.78 | 0.93× |
| D `emax` + `best_ci` u32, SoA comps | 1.54 | 1.07× |
| E control: sequential stream, no gather at all | 0.91 | 1.81× |

The build is already within **1.8×** of the same arithmetic with no gather at
all, and the best real variant captures 7% of that. Two of three compaction
designs are *negative*. Repeating at nraw 1e6 and 3e6 does not change the
ranking. **Do not spend effort on the build's data layout.**

## Result 2: bounded pruning — modelled a win, measured a loss

The design: hold each raw's candidate list lambda-descending (static, since
lambda never changes), scan the top-`T` clusters by reads cluster-major to seed
the map, then scan the rest raw-major with the exact early exit
`lambda × reads_T < best_e`. A single global `max_reads` bound was modelled and
rejected first — power-law cluster abundance makes it useless (0.83×).

Modelled at **8.9–23.5% of candidates examined, 1.06–4.24×, never a loss**
across a 24-cell sweep. Measured on a 384-sample pooled MiSeq run (`dada_fwd` =
R1, `dada_rev` = R2), via `dev/run_shuffle_ab.sh`:

| | base R1 | branch R1 | base R2 | branch R2 |
| --- | ---: | ---: | ---: | ---: |
| candidates examined | 3.68e9 (100%) | 2.56e9 (**44.4%**) | 3.39e9 (100%) | 2.28e9 (**43.0%**) |
| build ns/comp | 6.83 | **14.09** | 6.86 | **14.46** |
| build time | 25.10s | **36.07s** | 23.27s | **32.95s** |
| reconcile ns/comp | 12.08 | 12.37 | 13.30 | 13.17 |
| shuffle total | 43.22s | 55.53s | 52.43s | 62.61s |

Wall: `dada_fwd` +6.2%, `dada_rev` +5.4%, **pipeline total +4.8%**. Effective
cores −4.7% (the phase is serial, so a slower build directly lowers utilization).
Peak RSS unchanged (+0.1–0.2% on the denoise steps).

**The prune worked — better than the local fixtures suggested — and it did not
matter.** The bounded pass runs at 14.09 ns/comp against the sequential build's
6.83. Note it lands right on top of the *reconcile's* 12.08–13.30 ns/comp: pass 2
is the reconcile's access pattern, and inherits its cost.

Break-even therefore required

```
examined × 12.4 ns  <  total × 6.83 ns    =>   examined < ~35% of candidates
```

Achieved: 43–44%. And even a *perfect* pass 2 with zero overhead, running at the
pure scattered rate, would have cost 2.56e9 × 12.08 ns = 30.9s against the
baseline's 25.10s — **still a loss**. The gap between 14.09 and the 12.4 floor is
a per-raw index chase paid on every raw of every build even when the bound prunes
immediately; removing it entirely would not have changed the verdict.

### Why the model was wrong

The model priced the bounded pass at ~1.7 ns/comp. That was measured on Apple
silicon, where the per-raw map is cache-resident and the gather is nearly free.
The target's scattered raw-major access costs **11.4–13.4 ns/comp** — a number
[#87][87] had already measured and which this very page quoted while presenting
a model that contradicted it.

**The transferable lesson: model an access-pattern change with this repo's own
measured ns/comp constants, never with a laptop's.** The modelled *examined
fraction* was hardware-independent and still came in optimistic (9–24% modelled,
43–44% measured, because real pools reach thousands of clusters where `T=32`
covers a fraction of a percent). But the fatal error was the per-comparison cost,
and it was avoidable from data already in this repository. This is the **fifth**
scattered-access negative result here, and the first where the wrong constant was
sitting in an adjacent findings page the whole time.

## Result 3: the full phase accounting

Closing the build scan pointed #124 at the reconcile, and that prompted a
question nobody had answered: arithmetic on the A/B logs showed **15–19% of
shuffle time was unattributed** in every arm. The missing phase was the **move
pass**, which had never been timed. Verbose-only, behaviour-neutral
instrumentation was added (commit [`302c60f`][302c60f]) and both benchmark
platforms re-run.

362-sample pooled MiSeq (R1, R2) and the 95-sample PacBio HiFi set:

| | MiSeq R1 | MiSeq R2 | PacBio |
| --- | ---: | ---: | ---: |
| merged uniques | 272,574 | 296,893 | 547,273 |
| clusters | 2,564 | 2,007 | 2,844 |
| `run_dada` | 307.9s | 206.6s | 953.4s |
| `compare` (map parallel) | 254.7s | 144.7s | **877.0s** |
| **shuffle** | **46.85s** | **55.47s** | **53.20s** |
| **shuffle as % of `run_dada`** | **15.2%** | **26.8%** | **5.6%** |
| build | 27.28s (58.2%) | 25.49s (46.0%) | 19.76s (37.1%) |
| build ns/comp | 7.29 | 7.47 | **3.53** |
| reconcile | 10.37s (22.1%) | 21.19s (38.2%) | 21.58s (40.6%) |
| reconcile ns/comp | 12.79 | 14.13 | 12.50 |
| reconcile recomputes | 162.5M | 274.8M | 514.5M |
| …that changed cluster | 55,279 (0.034%) | 94,941 (0.035%) | 30,280 (**0.006%**) |
| move | 8.66s (18.5%) | 8.58s (15.5%) | 11.21s (21.1%) |
| move ns/raw | 4.37 | 4.52 | 3.07 |
| raws scanned per move | 1,973 | 1,865 | 2,857 |
| other (bookkeeping) | 0.53s (1.1%) | 0.22s (0.4%) | 0.65s (1.2%) |

The residual is ~1%, so this is the whole phase.

### The premise was wrong

**`b_shuffle` is 15.2% / 26.8% / 5.6% of `run_dada`** — not 46–55% of pooled
wall. #124's Amdahl argument (kill shuffle, serial fraction 62% → 25%) does not
survive contact with these runs. `compare` is 83% / 70% / 92% of `run_dada`.
Deleting `b_shuffle` outright would buy 15% / 27% / 6%, and on the platform where
pooled runs hurt most — PacBio, at 953s — it would buy the least.

This is the result that closes #124. The two falsified redesigns each cost real
effort; the number that made them worth attempting had never been re-measured.

### The constants generalise

Reconcile costs **12.50–14.13 ns/comp** across two platforms and a 143× range in
raw count. The ~12.4 ns/comp scattered constant is a property of this codebase's
access pattern, not of one dataset, and is safe to design against.

The build, meanwhile, gets **cheaper** with scale — 3.53 ns/comp on PacBio's
547k raws against ~7.4 on MiSeq's 273–297k. Independent confirmation of Result 1:
the cluster-major build is at its floor and improves as the pool grows.

### Reconcile redundancy is enormous and not exploitable

Only **0.006–0.034%** of reconcile recomputes actually change the raw's cluster —
and the fraction falls as the pool grows, so the waste scales *against* us. That
is far more nominal headroom than the build ever had.

It is not reachable. The reconcile already runs at the scattered rate, so there
is no cheap baseline to lose against (the asymmetry that sank Result 2 does not
apply) — but identifying which 0.03% of raws changed requires touching them, at
~12.4 ns/comp, which is the cost being avoided. A **perfect oracle** costing
nothing would save 3.4% / 10.3% / 2.3% of `run_dada`; any real mechanism gets a
fraction of that.

### The move pass: the only favourable cost shape, and still too small

Now visible for the first time: 15–21% of shuffle, scanning 1,865–2,857 raws for
every move it produces. It runs once per shuffle call over all `nraw`.

It has the one mechanism here whose cost model does not immediately sink it. A
raw can only move if its `compmax` changed, and `compmax` changes only in the
build or the reconcile. After a build every raw is dirty, so those passes are
irreducible — but the rest need only the reconcile-changed set, which is tiny in
absolute terms (55k, 95k, 30k for an entire run). Builds are 2,563 of 7,266 calls
(R1), 2,006 of 6,397 (R2), and 2,843 of 6,679 (PacBio), so **57–69% of move
passes are reducible to near-zero work**, and unlike every prior attempt the
replacement set is small enough that the scattered rate does not matter.

Sized: ~5.6s / 5.9s / 6.4s ⇒ **1.8% / 2.9% / 0.7% of `run_dada`**.

**Not built.** A sub-3% prize does not justify a change that must preserve
byte-identical move ordering (ascending `ci`) — precisely the invariant where
Result 2 hid a latent bug. Recorded here so it does not have to be rediscovered,
and so that anyone reopening it starts from the sizing rather than from the
mechanism.

## What the guardrail bought

The change was reverted, and the exactness work was still worth it.

The implementation was held to byte-identical output and cleared it: **2514
ASVs, churn = 0** on the real pooled run, plus three unit tests, each
mutation-checked to confirm it fails when its invariant is broken.

That discipline caught a **real latent correctness bug in the design itself**.
Two clusters can reach the same `e` from different `lambda × reads` pairs; since
the candidate list is ordered by lambda, the tied candidate with the lower `ci` —
the one the full scan keeps — can sit *later* in the list. Exiting the scan on
`<=` pruned it away and silently returned the wrong cluster. It needs exact
float equality, so it would very likely have survived a benchmark and a
concordance run, then surfaced months later as an unexplained ASV diff.

The payoff is that **the performance verdict above is unambiguous**. Because the
partition is provably identical, +4.8% wall is a statement about access patterns
and nothing else — not a slower build that also happened to cluster differently.
A perf A/B whose two arms might disagree on output cannot be interpreted at all,
and would have left this question open instead of closed. Establishing identity
first is what makes a negative result *final* rather than merely discouraging.

## What this dictates

- **`b_shuffle` is done.** [#124][124] closes as fully assessed. Build scan
  falsified twice (layout, bounded pruning); reconcile has huge redundancy and no
  mechanism to reach it; the move pass has a mechanism but a <3% prize. Do not
  reopen without both a fundamentally different mechanism *and* a fresh
  measurement of the phase's share.
- **`b_compare` is the target.** 70–92% of `run_dada` on every arm measured, and
  877s of PacBio's 953s. Any further pooled-`dada` performance work starts there.
- **Re-measure the premise before costing an optimisation, not just the design.**
  Both #124 redesigns were correctly executed against a share-of-wall figure that
  had gone stale. A cost model gates *how* you build; re-measuring the phase gates
  *whether* to.
- **Design constants for anything that follows:** sequential build ~3.5–7.5
  ns/comp, scattered raw-major ~12.4–14.1 ns/comp, sequential move scan ~3–4.5
  ns/raw. Anything scattered must prune below ~35% of candidates merely to break
  even.
- **Keep the instrumentation and the harnesses.** The `shuffle phases` verbose
  table is behaviour-neutral and is what produced Result 3;
  `dev/shuffle-model/` and `dev/run_shuffle_ab.sh` both survive the revert.

## Reproducing

```sh
# model (recalibrate its ns/comp to the target before trusting a cell)
cd dev/shuffle-model && rustc -O build_model.rs -o build_model && ./build_model

# the A/B that produced Result 2
bash dev/run_shuffle_ab.sh illumina /path/to/miseq_raw <base-bin> <branch-bin> out 32

# the phase table (Result 3) — any pooled run, no A/B needed
dada2-rs dada-pooled ... --verbose 2>&1 | grep -A4 "shuffle phases"
```

The bounded-pruning implementation is preserved in git history at
[`53d7add`][53d7add] and its revert at [`ff40014`][ff40014].

[86]: https://github.com/HPCBio/dada2-rs/pull/86
[87]: https://github.com/HPCBio/dada2-rs/issues/87
[124]: https://github.com/HPCBio/dada2-rs/issues/124
[53d7add]: https://github.com/HPCBio/dada2-rs/commit/53d7add
[ff40014]: https://github.com/HPCBio/dada2-rs/commit/ff40014
[302c60f]: https://github.com/HPCBio/dada2-rs/commit/302c60f
