# Measuring on a NUMA node: how a 5% win looked like nothing for two months

*Issue [#127](https://github.com/HPCBio/dada2-rs/issues/127). A methodology
result, not an algorithm one — but it changes the status of at least one earlier
verdict and possibly two.*

## The short version

The benchmark node has **two NUMA domains of 64 physical cores each**, with
remote memory at distance 32 against local 10. Nothing was ever bound to a
domain. So page placement was decided by first touch — in practice by whichever
threads of the 128-wide parallel map got there first — and it **re-rolled on
every run**.

The result was a run-to-run variance of **25% on the single-threaded scattered
phases**, against changes worth 3–6%. Under those conditions a real improvement
is indistinguishable from nothing, and we twice concluded it *was* nothing.
Fixing the placement — by any policy — removes it.

| control phase (untouched by the change under test) | default placement | fixed placement |
|---|---|---|
| `b_shuffle` | 53.6–66.9 s (**25%**) | 65.3–65.5 s (**0.3%**) |
| `store` | 17.6–21.6 s (**23%**) | 20.2–20.2 s (**0.0%**) |
| `run_dada` | 301.9–324.7 s (7.5%) | 372.2–372.7 s (**0.13%**) |

(Fixed placement here is `--interleave=all`; see
[below](#fixing-placement-is-a-measurement-tool-not-a-tuning-recommendation)
for why that rather than binding.)

With placement fixed, `main` reproduces itself to **0.04–0.45%** on the DP
kernel and the parallel map. Signal-to-noise on the same question went from
under 1× to roughly 10–100×.

## What it cost

The [rolling-`d16` DP change](compare-screen-vs-align.md) removes ~5/6 of the
DP kernel's memory traffic. Measured under default placement it looked like this:

| measurement | verdict at the time |
|---|---|
| 24-thread production A/B, both platforms | "flat — MiSeq +1.9%, PacBio −2.8%, both inside a 2–6% noise floor" |
| 128-thread PacBio pair | "−5.9%, but the untouched controls moved 12–19%, so unattributable" |
| 128-thread A/B/A/B | "**null** — the two paired estimates disagree by 6 points; close the branch" |

Measured with placement fixed (both arms bound to one domain), the same change
on the same data:

| threads | rep1 | rep2 | screen (untouched control) |
|---|---|---|---|
| 24 | −3.71% | −3.39% | −1.1% / −0.9% |
| 64 | −4.72% | −4.73% | −1.1% / −0.5% |
| 128 (SMT) | −5.80% | −6.18% | −1.1% / −0.1% |

Replicates of each arm agree to **0.4 percentage points**, the effect scales
monotonically with thread count as a memory-pressure mechanism should, and it
is confirmed independently under `numactl --interleave=all` (−4.3%). The change
was always real. What changed was the rig.

Note that the A/B/A/B design was not the problem. It was built to cancel
*monotonic drift* — and the noise here is not drift, it is placement, which
re-rolls per run and cancels nothing. A better design cannot rescue a
measurement whose variance exceeds its signal; only reducing the variance can.

## Fixing placement is a measurement tool, not a tuning recommendation

Three policies at 64 threads, `main` only, two replicates each. The middle
column is what `dev/numa_pin.sh` now does by default; the right column is what
production does today.

| | `--cpunodebind=0 --membind=0` | `--interleave=all` | default (first touch) |
|---|---|---|---|
| DP kernel ns/comp | 149,747 | **119,776** | 120,878 |
| parallel map | 286.8 s | **235.4 s** | 242.9 s |
| `b_shuffle` (serial) | **50.4 s** | 65.3 s | 55.2 s |
| `store` (serial) | **12.8 s** | 20.2 s | 18.8 s |
| `run_dada` (mean) | 391.4 s | **372.5 s** | 372.2 s |

**Binding is the outlier, and it is the obvious fix that turns out to be wrong.**
It forces all 64 threads to draw bandwidth from one node's memory controllers,
costing 24% on the DP kernel and 21% on the parallel map. Interleaving and
default placement are within 1% of each other on the kernel — first touch on a
64-thread map already spreads pages much as interleaving does — so **the default
is not leaving throughput on the table.**

The phases want opposite policies. The parallel map wants pages spread across
both domains for bandwidth; the single-threaded scattered `b_shuffle` and
`store` want them local, and pay 24–58% for interleaving. Net, on `run_dada`,
spreading wins.

What separates the three is **reproducibility**, and that is the whole point:

| replicate-to-replicate spread | `run_dada` | `b_shuffle` | `store` |
|---|---|---|---|
| `--interleave=all` | **0.13%** | **0.31%** | **0.0%** |
| `--cpunodebind`/`--membind` | 0.90% | 5.6% | 1.6% |
| default (first touch) | 5.8% | 12.9% | 22.9% |

`--interleave=all` is the most reproducible of the three — more so than binding
— which is why it is the default in `dev/numa_pin.sh`. Two further reasons to
prefer it for benchmarking:

- **It does not distort what is being measured.** Binding makes the parallel map
  21% slower, inflating its share of `run_dada` and shifting exactly the phase
  ratios these benchmarks report. Interleaving stays within 3% of default
  placement.
- **It cannot oversubscribe.** `--interleave=all` sets memory policy only, never
  CPU affinity. Binding restricts cores too, so a thread count sized for the
  whole node silently becomes SMT contention — which has already happened here
  once, in a 128-thread run into a 64-core domain, caught only by the
  `cpu allocation` line in `--verbose`.

**For production, this is a predictability knob, not a speed one**: interleaved
and default placement tie on mean `run_dada` (372.5 s vs 372.2 s). What
interleaving buys is that the run takes the same time twice.

A useful side result: bound-64 (149,747 ns/comp) lands almost exactly on the
*default-placement* 96-thread value (144,845), and interleaved-64 (119,776) on
the 48-thread value (116,199). **Binding 64 threads to one domain reproduces the
memory pressure of 96 threads across two** — independent confirmation that the
DP kernel's degradation at high thread counts is bandwidth rather than anything
about scheduling.

## None of this was visible without the phase instrumentation

Worth drawing out, because it is an argument for instrumenting phases you have
already decided *not* to optimise.

The variance was detectable only because every run reports timings for phases
the change under test could not possibly have affected. For a DP-kernel change
those controls are the k-mer screen (same parallel map, same data, different
code path) and `b_shuffle` / `store` (serial, different phase entirely). A
control that moves *more* than the signal is the cheapest possible refutation,
and it fired twice here before anything wrong was published.

Those channels exist as a side effect of two efforts that concluded elsewhere:

- The `shuffle phases` table came from [#124](shuffle-build-scan.md), which
  closed `b_shuffle` as a perf target. The instrumentation was kept anyway, and
  became the control channel for a `b_compare` experiment with no connection to
  `b_shuffle`.
- The `compare split` came from [#127](compare-screen-vs-align.md) and supplies
  the screen-vs-DP separation the A/Bs are read through.
- The `cpu allocation` line was added late in #127 to make archived runs
  self-describing about their allocation. It immediately caught a 128-thread run
  binding into a 64-core domain — an SMT oversubscription that would otherwise
  have looked like a regression in the change under test.

The general form: **instrument phases you are not optimising**, because their
only job in a later experiment is to stay still, and a still control is what
converts "the number moved" into "the change did it".

## Which earlier results this puts in question

Every A/B this project has run on that node carries the same uncontrolled term.
The exposure is not uniform — it scales with how much scattered, latency-bound
memory access a change touches:

- **[#124's bounded pruning arm](shuffle-build-scan.md) (+4.8% wall, reverted).**
  The most exposed result we have. It is a `b_shuffle` change — the exact phase
  that varies 25% under default placement — and its verdict sits near that
  floor. The
  *mechanism* finding behind it stands on its own (the bounded pass provably
  runs at the scattered rate, 6.83 → 14.09 ns/comp, and break-even needed a
  prune fraction the design could not reach), so the conclusion is probably
  safe. But the headline number should not be quoted to two significant figures,
  and if that arm is ever revisited it needs re-measuring with placement fixed.
- **The `b_shuffle` phase splits** (build/reconcile/move) and the ns/comp
  constants derived from them. Ratios within a run are less exposed than
  absolutes, since all phases in a run share one placement roll.
- **Anything quoted at "±2–6% noise floor"** in
  [the `b_compare` findings](compare-screen-vs-align.md). That floor was
  measured, not assumed — but it was measured under uncontrolled placement, so
  it describes the rig rather than the machine.

Less exposed: the `b_compare` *ratios* (screen vs align, DP vs `al2subs`), which
come from phases measured simultaneously in the same run, and every
count-based quantity (comparisons, screen pass rates, cells per alignment).

## What we do about it

`dev/numa_pin.sh` is sourced by `dev/run_shuffle_ab.sh` and sets a `numactl`
prefix applied to **both arms** of an A/B. It degrades to a no-op without
`numactl`, on single-domain machines, and where `numactl` exists but the policy
cannot actually be applied (containers, restricted cpusets).

```bash
NUMA_POLICY=interleave   # default: numactl --interleave=all
NUMA_POLICY=bind         # --cpunodebind=N --membind=N; lowest latency for the
                         #   serial scattered phases, but slows the parallel
                         #   map 21%, and refuses to oversubscribe
NUMA_POLICY=none         # default first-touch — for comparing against results
                         #   measured before this existed
NUMA_NODE=0              # which domain, for NUMA_POLICY=bind
```

Under `bind` it **refuses when the thread count exceeds the physical cores in
the target domain**: binding restricts CPUs as well as memory, so a thread count
sized for the whole node silently becomes SMT contention that reads as a
regression in whatever is under test. That is not hypothetical — a 128-thread
run into a 64-core domain did exactly this, and only the `cpu allocation` line
in `--verbose` caught it. `interleave` has no such failure mode, which is part
of why it is the default.

## What this dictates

1. **Pin both arms of any A/B on a multi-domain machine.** It is one command and
   it moves the noise floor by more than an order of magnitude.
2. **Always report an untouched control.** Every conclusion here rests on a
   phase the change could not have affected — the k-mer screen for a DP change,
   `b_shuffle`/`store` for anything in `b_compare`. A control that moves more
   than the signal is the cheapest possible refutation, and it caught two of
   this issue's wrong turns before they were published.
3. **"Below the noise floor" is a claim about the rig.** It licenses "we cannot
   measure this here", never "there is no effect". The distinction is not
   pedantic: it is the difference between closing a branch and fixing a
   benchmark.
4. **Replicate the same arm, not just the two arms.** `main` against `main` is
   what makes the floor visible; without it there is no scale against which a
   difference means anything.
