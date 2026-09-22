# Threading the serial steps: the unit of parallelism is the sample

**The largest end-to-end speedups in this project are not from the denoiser.**
They come from running per-sample work concurrently where R DADA2 runs it
serially — and from one non-obvious correction: a single sample cannot feed a
24-thread pool, so the thing to run in parallel is *samples*, not the inner loop.

Four results, in decreasing order of how much they change what you should do:

1. **A single sample starves the pool.** Denoising one sample on all threads
   plateaued at ~6.8× speedup on 16 threads while burning **+72% CPU** — workers
   spin-waiting on a comparison map too small to divide. Running J samples on
   sub-pools of `threads/J` fixed it.
2. **Per-sample scaling does not predict aggregate throughput.** The single-sample
   thread sweep implied ~8 threads per sample; the samples-in-flight sweep put the
   plateau at **~4**. The default was retuned on that, not on the first sweep.
3. **The headline ratios are statements about the baseline.** `remove-primers` is
   284× faster than R's — because R's `removePrimers` is serial and ours fuses
   filtering into the same pass. That is a real wall-clock win and it is not a
   claim about our algorithm.
4. **A third of our measured wall time was once a harness artifact** — and it ran
   *against* us. See [the measurement trap](#the-measurement-trap) below.

## Where the time actually went

From the [benchmark scoreboard](../results.md#per-step-comparisons) — PacBio
pooled, k=7, dada2-rs vs R-single:

| Step | dada2-rs (s) | effective cores | R (s) | speedup |
|---|---:|---:|---:|---:|
| remove_primers | 19.0 | 20.9 | 5410.2 | **284.3×** |
| learn | 25.2 | 20.6 | 159.7 | 6.3× |
| dada | 698.8 | 18.8 | 6226.2 | 8.9× |
| make_table | 0.3 | 1.0 | 0.1 | **0.2×** |
| remove_bimera | 17.8 | 23.6 | 19.8 | 1.1× |
| **TOTAL** | **761.1** | 19.0 | 11885.4 | **15.6×** |

On Illumina the same shape appears in `merge-pairs`, at **34–59×** across pooling
modes. Both are steps where R processes samples one at a time.

**Read the bottom two rows as carefully as the top one.** `make_table` runs at 1.0
effective cores and is *slower than R*; `remove_bimera` sits at parity despite
using 23.6 cores. Threading bought nothing in either place — the first because the
work is a serial table assembly, the second because R already parallelises it.
The wins are concentrated exactly where R is serial and the work is independent
per sample, and nowhere else.

## A single sample cannot feed the pool

`dada-pseudo` originally denoised samples serially, each `dada_uniques` running on
the full thread pool. The thread sweep showed the failure clearly: speedup
plateaued at **~6.8× on 16 threads while CPU rose 72%**. By contrast `--pool
false`, which runs one sample per process on one thread each, reached **85%
parallel efficiency**. The pool was not short of work overall; it was short of
work *per sample*.

The fix is `for_each_sample_concurrent`: J workers, each owning a rayon sub-pool
of about `threads/J`, pulling samples from an atomic counter and calling
`sub_pool.install(|| dada_uniques(...))` so each per-sample map is pinned to a
right-sized pool. It now backs four call sites:

| Where | J | Why that J |
|---|---|---|
| `dada-pseudo` rounds | `threads/4` | per-sample denoising; both rounds fan, priors marked serially at the existing barrier |
| multi-input `dada` | `threads/4` | single-pass load → denoise → write; also bounds memory to J samples in flight |
| `dada-pooled` load | up to `threads`, ~1 thread each | derep is I/O and hashing, not DP — across-sample concurrency fills cores better than dividing one derep |
| `merge-pairs` | rayon `par_iter` over samples | samples are independent; `collect` preserves input order |

**J is a property of the phase, not a global setting.** Denoising wants a few
threads per sample; the pooled derep front wants one. Both were measured, not
assumed.

## Per-sample scaling does not predict aggregate throughput

This is the part worth carrying to other work. The single-sample thread sweep
suggested ~8 threads per sample was the right size, and the first default was set
to `round(threads/8)`. The samples-in-flight sweep then disagreed: the wall-time
curve plateaus at **~4 threads per sample**, because aggregate throughput prefers
more samples in flight than per-sample scaling does.

At 24 threads that is J=3 (148 s) against the J=6 plateau (135 s) — **~9% left on
the table by tuning on the wrong curve**. The default became `round(threads/4)`,
worth ~6% end-to-end on pseudo. J=8 was marginally better still but within
run-to-run noise of J=6 and cost more memory, so it was not taken.

The general form: *the efficiency of one unit of work under N threads tells you
very little about the throughput of many units sharing N threads.* The later
[thread-scaling study](thread-scaling-and-placement.md) is the developed version
of this, and finds the best thread count differs two-fold between two reads of the
same pool.

## Correctness: byte-identical, and tested for it

Every one of these changes is output-identical regardless of how many jobs run,
and each shipped with a determinism regression test comparing `--threads 1`
against `--threads N`:

- `merge-pairs` — `collect` preserves input order.
- `dada-pooled` load — results reassembled by input index before inference.
- `dada` / `dada-pseudo` — `dada_uniques` is deterministic and files are written
  by sample name; round-1 prior selection is a set union, so it does not depend on
  completion order.

A threading change that alters output is a bug, not a tradeoff. Nested rayon is
safe here because a sub-pool `install` runs inline and work-steals rather than
oversubscribing.

## The measurement trap

The benchmark harness once ran `filter-and-trim`, `remove-primers`,
`remove-bimera-denovo` and `merge-pairs` single-threaded while giving R
`multithread=N`. On a 24-thread MiSeq run, `filter` and `remove_bimera` together
accounted for **~1/3 of dada2-rs's wall time purely from that asymmetry** — an
unfairness pointing *at* us, which is the direction that does not get caught by
wishful thinking.

The lesson is the one the [NUMA work](measuring-on-numa.md) later made
unavoidable: **a speedup ratio is a statement about two configurations, not about
two implementations.** Quote it with what both sides were allowed to do. The
open form of this question — whether R gains as much from NUMA binding as we do
— is [#155](https://github.com/HPCBio/dada2-rs/issues/155), and until it is
answered the binding results are ours alone, not a comparative claim.

## The instrument came first

Worth noting the order of events. The effective-cores column — CPU time over wall
time, per step — was added to the harness *before* the threading work, together
with `--thread-sweep`. Nothing in the tables above is visible without it: a step
can be 3× faster than R and still be leaving 90% of the machine idle.

That same column is what later exposed pooled `dada` running at 9.9–12.9 of 24
effective cores, which started the Amdahl analysis and, from there,
[the entire `b_shuffle` and `b_compare` arc](compare-screen-vs-align.md). The
cheapest instrument in this project found the most expensive problem in it.

## What this dictates

- **`--sample-jobs` defaults to `threads/4`** on `dada` (multi-input) and
  `dada-pseudo`, and to 1 at ≤4 threads, which is the old serial path. Raise it
  for many small samples; lower it if J samples in flight strains memory, since
  it also bounds the resident set.
- **Do not quote an end-to-end speedup without the configuration.** The 284× is
  real and is mostly a fact about R's serial `removePrimers`; `cutadapt` closes
  much of that gap for anyone who wants it closed.
- **Check effective cores before optimising anything.** A low value means the
  machine is idle and the fix is scheduling; a high value with a poor ratio means
  the fix is in the algorithm. These call for opposite work, and the column tells
  you which one you have.
