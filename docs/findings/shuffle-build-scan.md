# Pruning the `b_shuffle` build scan

**Verdict: NEGATIVE — bounded pruning of the shuffle build is falsified, and
`b_shuffle`'s build phase is now believed to be at its floor.** A two-level
bounded scan was modelled at 1.06–4.24× and built to a byte-identical
implementation, then measured on a 384-sample pooled MiSeq run: it **pruned 56%
of the comparisons and still lost**, because each surviving comparison cost
**2.06× more** (6.83 → 14.09 ns/comp). Net **+4.8% wall**. Reverted.

Two things are worth extracting, and the second matters more than the verdict:

1. The break-even is **structural, not a tuning problem**. The bounded scan's
   access pattern costs ~12.4 ns/comp — a rate this repo had *already measured*
   — so beating a 6.83 ns/comp sequential build required pruning to **under
   ~35%** of candidates. At the achievable 44%, no implementation of this design
   could have won. Tuning the bound was never going to close it.
2. **The exactness discipline paid for itself even though the change was
   thrown away.** See "What the guardrail bought" below.

Layout rework of the same phase is separately falsified (Result 1), so together
these close the build scan as a target and move [#124][124] to the reconcile.

## Background

After [#86][86], `b_shuffle` is 46–55% of pooled wall and the Amdahl ceiling on
pooled runs. It splits into a **build** (cluster-major rebuild of `compmax`,
paid fresh each bud round) and a **reconcile** (raw-major, scattered).
[#87][87] tried to remove the build by carrying `compmax` across buds and was
closed: it relocated work into the more expensive scattered path, and its sign
flipped between two amplicons of one run (`shuffle-compmax-carry.md`).

That left the phase itself as the target. #124's starting observation is that
during a shuffle **`lambda_ij` is fixed** — only the scalar `reads_j` changes —
so the phase is `argmax_j (lambda_ij × reads_j)` with one factor static,
evaluated billions of times to discover a few hundred moves per build.

## Result 1: layout rework is a dead end (holds)

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

- **The build scan is closed.** Layout is falsified (Result 1); bounded pruning
  is falsified (Result 2). Both bets have now been paid for. Do not reopen
  without a fundamentally different mechanism — and price it in ns/comp against
  the measured constants before writing code.
- **Any future bet here must beat ~6.8 ns/comp sequential.** That is the number
  to design against. Anything that touches raws in a scattered order starts at
  ~12.4 and must prune below ~35% of candidates merely to break even.
- **[#124][124] should re-point at the reconcile.** It is now the larger half of
  the phase (24–40% of scanned comparisons, 10.2–20.2s of a 55–63s shuffle) and
  is untouched. It is already raw-major and already paying the scattered rate, so
  unlike the build it has no cheap-access baseline to lose against — the
  asymmetry that sank this attempt does not apply there.
- **Keep the harnesses.** `dev/shuffle-model/` and `dev/run_shuffle_ab.sh` both
  survive the revert and are the measurement path for anything that follows.

## Reproducing

```sh
# model (recalibrate its ns/comp to the target before trusting a cell)
cd dev/shuffle-model && rustc -O build_model.rs -o build_model && ./build_model

# the A/B that produced the table above
bash dev/run_shuffle_ab.sh illumina /path/to/miseq_raw <base-bin> <branch-bin> out 32
```

The implementation is preserved in git history at [`53d7add`][53d7add] and its
revert at [`ff40014`][ff40014].

[86]: https://github.com/HPCBio/dada2-rs/pull/86
[87]: https://github.com/HPCBio/dada2-rs/issues/87
[124]: https://github.com/HPCBio/dada2-rs/issues/124
[53d7add]: https://github.com/HPCBio/dada2-rs/commit/53d7add
[ff40014]: https://github.com/HPCBio/dada2-rs/commit/ff40014
