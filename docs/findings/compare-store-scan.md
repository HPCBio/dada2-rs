# Inside `b_compare`: the store is a cache-line problem

**Verdict: confirmed and built.** The serial store in `b_compare` was not doing
serial work — it was walking a 160-byte struct to read 8 bytes of it. Moving
that one field into a dense array cuts the store by **71%** and `run_dada` by
**20.1–21.3%** on NovaSeq soil 16S, byte-identical.

## What was left to explain

[The serial remainder](compare-serial-fold.md) folded six reduction passes into
the store loop for −24 to −30% of `run_dada`, and ended by naming what remained:
the store itself, at a "suspiciously flat 18.4–19.6 ns/raw-visit across a 5×
range in pool size."

That flatness was the clue. A phase whose cost per unit does not move when the
working set moves five-fold is not bound by residency. It is paying a fixed
per-access price — which is what a cache miss on every access looks like.

## The commit path is negligible; it is all scan

Before measuring anything on a cluster, the arithmetic settles what the store
*is*. On soil 16S R1 the store walks 11,888,798,623 raw-visits and commits
12,082,545 of them: **0.102%**. Whatever the store costs, 99.9% of it is the
scan over the result vector and the per-raw state it touches, not the work done
on a hit. So the target is the scan, and the commit path can be ignored.

This is the cost-model-first discipline the project keeps re-learning: the
question "what is this phase made of" is answerable on paper here, and answering
it first avoided building against the wrong 0.1%.

## The microbenchmark: `examples/store_bandwidth`

Seven arms decompose the scan into the streams it actually performs, at the real
pool size (1,225,523 raws) and the real store rate (0.00102). Reads go through
`black_box` — the first draft reported 574 GB/s because the single-array arms
auto-vectorised while the paired arms could not, an artifact of the instrument,
not a result.

| arm | what it isolates |
|---|---|
| `both48` / `both16` | result vector + strided `Raw` read, at 48 B and 16 B elements |
| `dense48` / `dense16` | result vector + a *dense* `f64` array |
| `comps48` / `comps16` | the result vector alone |
| `stride` | the strided `Raw` read alone |

The answer was unambiguous and reproduced on two machines and under both the
`release` and `release-native` profiles: **the strided read is 83–87% of the
scan.** The result vector, which is a clean sequential stream, is nearly free.

`Raw` is 160 bytes. `e_minmax` is an `f64`. Reading it once per raw pulls a
64-byte line to use 8 bytes, and the raws are visited in an order that defeats
the prefetcher. Held instead in a `Vec<f64>` on `B` — the same parallel-array
pattern `B::raw_cluster` already uses — the identical reads become contiguous.

## Result

Two replicates per arm, NovaSeq soil 16S and soil ITS2, both reads, exclusive
node (compute-5-6, EPYC 7713), 64 threads, `numactl --interleave=all`.

| soil 16S | `main` | dense array | change |
|---|---|---|---|
| store R1 | 227.2 / 223.7 s | 63.8 / 63.6 s | **−71.7%** |
| store R2 | 275.4 / 277.9 s | 80.2 / 81.2 s | **−70.8%** |
| store ns/raw-visit R1 | 19.1 / 18.8 | 5.4 / 5.3 | 18.9 → 5.4 |
| store ns/raw-visit R2 | 18.7 / 18.9 | 5.5 / 5.5 | 18.8 → 5.5 |
| `run_dada` R1 | 834.2 / 820.5 s | 654.1 / 667.4 s | **−20.1%** |
| `run_dada` R2 | 926.4 / 960.5 s | 735.4 / 749.1 s | **−21.3%** |

Soil ITS2, the smaller and more variable pool, returns −69% of store time and
**−15.2%** of `run_dada` on both reads. All 300 per-sample outputs across the
six 16S runs are byte-identical once the embedded build hash is excluded, and
240/240 on ITS2.

The controls are flat. `b_shuffle` sits at 112.5–115.3 s across every 16S arm
and `reconcile` at 40.2–41.2 ns/comp — the untouched phases did not move, which
is what rules out a whole-run effect masquerading as a phase win.

The observable the user noticed first was core utilisation. Effective cores —
`(map thread-busy + serial wall) / run_dada` — rise from **31.4 to 38.4** on R1
and **29.4 to 36.5** on R2, of 64. Amdahl is the whole story of this phase, and
the store was the largest remaining serial block.

## What generalises

The synthetic under-predicted. `store_bandwidth` put the strided read at 83–87%
of the scan, which brackets a −6 to −17% `run_dada` win; the real deficit
multiplied out to −20%. Reported here because the direction of the miss matters:
the microbenchmark isolated the right *mechanism* while understating its
*magnitude*, so its ranking of levers was trustworthy even though its arithmetic
was not. Plan against a synthetic's ordering, not its numbers.

The reusable shape is the diagnostic, not the fix. **A per-unit cost that does
not move with working-set size is a cache-line problem, not a residency
problem** — and the fix is a layout change, not a smaller working set. The same
signature is worth looking for anywhere the profile shows a flat ns/unit across
pools that differ several-fold. `Raw` is 160 bytes and other phases read single
fields from it.

## What this leaves

The store is now 12–15% of compare, down from 33–37%, and 4.9–5.5 ns/raw-visit
is close to the sequential-stream floor the synthetic measured. There is no
large serial block left in `b_compare`: the parallel map is 66% of the phase and
bandwidth-bound above 48 threads, which is
[a different problem](compare-screen-vs-align.md) with a different lever.

Two further levers were built and measured, and **neither is merged.** They are
the rest of this page, because the reasons are more useful than the result.

## The allocation levers: a real cost that does not pay

`b_compare_parallel` allocates its `nraw`-long result vector on every call and
frees it — 45.1 MB on soil ITS2, 58.8 MB on 16S, over thousands of calls. Two
ways to stop paying for that were built:

- **Lever A** — collect at 16 B instead of 48 B in production runs, dropping
  `CompCost` from the element (`--verbose` still needs the 48 B form for its
  attribution). This works by *accident*: it puts the vector under a threshold.
- **Lever D** — hold one buffer on `B` and fill it with `collect_into_vec`,
  removing the allocation at either width.

Both collapse `sys` time — 70% for A, **94% for D** (674.6 → 38.6 s on ITS2 R2).
The mechanism is glibc's dynamic `mmap` threshold, which caps at 32 MB: a vector
above the cap is `mmap`ped and `munmap`ped every call rather than recycled from
the arena, so the kernel re-faults and re-zeroes every page each time.
`examples/parallel_overhead` confirms it directly — at 16 B the vector is
15.0–19.6 MB, under the cap, and fresh-vs-reuse shows no consistent difference;
at 48 B it is 45.1–58.8 MB, over the cap, and the gap is **~5 ms per call at
both pool sizes**.

That arithmetic is consistent end to end: 4,017 calls × 5 ms ≈ 20 s against a
131 s map, and 674 core-seconds of `sys` is ~10.5 s of wall across 64 threads.

**And removing it returns nothing.**

| soil ITS2 R2, non-verbose | baseline | D (scratch buffer) | A (16 B) |
|---|---|---|---|
| real | 258.3 s | 264.7 s (**+2.5%**) | 257.7 s (−1.5%) |
| user | 7,381 s | 8,592 s (**+16%**) | 8,536 s (+11%) |
| sys | 674.6 s | 38.6 s (−94%) | 202 s (−70%) |

Both levers convert a large `sys` saving into an equal-or-larger `user` cost and
a wall time that does not improve. Total CPU goes *up*: 8,056 core-seconds
becomes 8,631 (D) or 8,738 (A).

The leading account for the `user` rise — untested, and stated as a hypothesis —
is that a freshly-`mmap`ped page is zeroed by the kernel immediately before a
worker writes it, so those lines are cache-warm and the store is nearly free,
while a reused buffer's lines are cold and dirty and every store pays a
read-for-ownership plus a writeback. Same work, moved out of kernel zeroing and
into user-mode stalls, and *more* of it, since the kernel zeroes with
non-temporal stores and a `collect` does not.

`parallel_overhead` cannot settle that, and says so: its synthetic reports reuse
as **faster**, the opposite of production. Its work function is pure ALU while
production's map streams 1.7 GB of k-mer vectors, so it has nothing to contend
with. Same failure mode as `screen_bandwidth` — a stand-in reproduces ordering,
not magnitudes, when it stands in for real work.

**Consequence: do not retry either lever on the strength of the `sys` number.**
The `sys` time is real, well-understood, and not on the critical path. Anything
that revisits this needs to explain the `user` rise first, with hardware
counters rather than a timer.

## Load imbalance in the map: falsified

The same investigation was pointed at a second target. Production reports 86–90%
map parallel efficiency, and the residual looked like threads idling at the tail
of each collect waiting on stragglers — plausible, since raws are
abundance-sorted and a few percent pay a full alignment while the rest are
screened out.

It is not happening. `join_uniform` and `join_skewed` do the **same total work**
per call and differ only in its distribution:

| | ITS2 shape (0.8% heavy) | 16S shape (4.5% heavy) |
|---|---|---|
| skewed vs work-matched uniform | +0.1% | −0.1% |

Zero, on both parameterisations, at production's measured 14× aligner/screen
cost ratio. The heavy region is ~235 tasks of 29,360 spread over 64 threads, so
`with_max_len(32)` gives work-stealing far more splits than it needs.

**So the 86–90% figure is not idle time**, and an earlier reading of it here as
"1,166 core-seconds lost to imbalance" was wrong. `busy` is the sum of timers
taken *inside* the map closure: it never measured rayon's per-task dispatch, the
collect's stores into the destination, or the timer calls themselves. The gap is
mostly unmeasured work. It is bounded at ≤14% of the map and separating it needs
a different instrument.

**Consequence: the load-balancing knobs are closed.** `DADA2RS_PAR_GRAIN` is
already doing its job, and there is no parallel-dampening problem to chase.

## Three instrument errors, and what each cost

This page's negative results all came from measurement mistakes caught late, so
they are recorded rather than quietly fixed:

1. **A hypothesis killed by accounting, not by a run.** Rayon spin-wait was the
   first explanation for the `user` rise. The whole ITS2 R2 run has only 838
   core-seconds of CPU outside the map's summed busy time, while lever D's
   `user` rise alone is 1,211 — there was no room for it. Doing that arithmetic
   before booking node time is the cheapest step in this entire page.
2. **A skew model that could not exhibit the thing it measured.** The first
   imbalance arm marked every 32nd index heavy while the task grain was 32, so
   every task held exactly one heavy item and mean task cost equalled max task
   cost *by construction*. It reported "no imbalance" and would have done so
   whatever the truth was. A skew uniform at the granularity of a task is not a
   skew.
3. **The right estimator applied to the wrong arm.** Reporting the per-round
   minimum defends against a noisy neighbour, and is correct for five of the six
   arms. `collect_fresh`'s spread is *intrinsic* — glibc allocator state, i.e.
   whether a round recycled rather than remapped — so the minimum systematically
   selects the rounds where the effect under test did not occur. On `min` the
   allocation cost appeared to shrink as the vector grew (2.5 ms at 45.1 MB,
   1.3 ms at 58.8 MB); on `median` it is 5.04 and 5.10 ms, essentially identical.
   The tell was in the output the whole time: that arm runs at `med/min`
   1.06–1.10 while every other arm sits at 1.00–1.01.

The generalisation, and the reason this section exists: **an arm that cannot
produce the effect, and an estimator that selects against it, both return a
clean-looking null.** A null is only evidence once the instrument has been shown
capable of returning something else.
