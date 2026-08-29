# Tuning for your data

**How many threads should I give `dada-pooled`, and does pinning memory help?**

This page answers both from measurement rather than from the usual advice, which
is "use all your cores." On one of our benchmark pools that advice costs 5% of
wall time, and on another it leaves 45% on the table. The difference is a
property of your data that `--verbose` already prints.

!!! warning "This guide is incomplete by construction"
    Every number here comes from four benchmark arms on one machine. Your
    amplicon, your instrument, and your hardware all move the answer. Treat the
    **rules** as portable and the **numbers** as illustrative — and if you have
    the node time, spend twenty minutes measuring your own data instead of
    trusting ours. [Diagnostics](diagnostics.md) has the tooling.

## The one number that predicts the rest

`dada-pooled --verbose` reports how the denoiser's inner loop divides:

```
[dada]   kmer screen  11078.36s (50.8%)    11605515622 comps    955 ns/comp  (every comparison)
[dada]   align total   9088.68s (41.7%)      522407633 comps  17398 ns/comp  (4.5% passed the screen)
```

Those two shares behave completely differently:

- The **k-mer screen** runs on every comparison and streams large per-sequence
  vectors. It is limited by how fast memory can be reached.
- The **aligner** runs only on pairs the screen lets through, and does dense
  arithmetic. It is limited by how many cores you have.

So a pool where the screen dominates is a **memory-bound** workload, and a pool
where alignment takes a large share is a **compute-bound** one. That single
distinction drives everything below.

**What moves it, in our data:**

| | screen share | screen pass rate |
|---|---|---|
| Soil 16S (V3–V4), read 1 | 50.8% | 4.5% |
| Soil 16S, read 2 | 67.7% | 2.4% |
| Soil ITS2, read 1 | 85.8% | 0.8% |
| Soil ITS2, read 2 | 83.2% | 0.9% |

Two things are worth noticing. **Amplicon and diversity dominate**: 16S is a
conserved gene at moderate diversity, so many sequences resemble each other and
clear the screen; ITS2 is hypervariable at high diversity, so few do. These two
pools came from the *same sequencing run and the same samples*, so that gap is
biology and amplicon design, not instrument.

But **read position matters too**, and by a surprising amount. 16S read 1 and
read 2 are the same library and the same organisms; quality decay along read 2
changes the unique-sequence structure enough to move the screen share by 17
points — and, as below, to change the best thread count.

## Threads

Measured on a 128-core EPYC 7713 (two NUMA domains), two replicates per point,
`run_dada` seconds:

| threads | 16S R1 | 16S R2 | ITS2 R1 | ITS2 R2 |
|---|---|---|---|---|
| 24 | 1072.5 | 1076.3 | 221.1 | 308.5 |
| 48 | 713.8 | 767.8 | **194.0** | **254.9** |
| 64 | 638.9 | **711.5** | 198.2 | 258.2 |
| 96 | **594.3** | 734.9 | 203.8 | 263.2 |

**The compute-bound pool keeps rewarding cores; the memory-bound one stops at 48
and then gets slower.** Past its knee, adding threads is not merely wasted — it
is actively harmful, because more threads contend for the same memory.

### Rule of thumb

| your screen share | suggested threads |
|---|---|
| below ~55% (align-heavy) | scale freely; we saw gains to 96 |
| ~55–75% | expect a knee around 64 |
| above ~75% (screen-heavy) | ~48; more will cost you |

Below 24 threads we saw no knee on any pool: **if you have 24 cores or fewer,
use them all.** Most of this page only starts to matter above that.

### If you have a queue rather than one urgent sample

Wall time is not the only cost. Parallel efficiency relative to 24 threads:

| | 48 | 64 | 96 |
|---|---|---|---|
| 16S R1 | 75% | 63% | 45% |
| ITS2 R1 | 57% | 42% | 27% |

Every arm we measured buys wall time at roughly **four times the core-seconds**.
So the two questions have opposite answers:

- **One sample, and you are waiting on it** → use the thread count from the table.
- **Many samples through a shared machine** → use *fewer* threads per job and run
  more jobs. Throughput is maximised near the low end of this range, not the high.

## Memory placement (NUMA)

!!! tip "This is the one recommendation that does *not* depend on your data"
    Unlike thread count, the NUMA effect held at nearly the same size on both
    pools — 25% on the compute-bound one and 27–29% on the memory-bound one —
    despite a 35-point gap in screen share. We expected the compute-bound pool to
    gain much less, because arithmetic can hide memory latency; it did not.
    (16S read 2 is still outstanding.) If you have a multi-domain machine, this
    is worth trying whatever your amplicon is.

On a multi-socket or multi-die machine, memory is faster from some cores than
others. By default the kernel scatters a job's pages across all of it, so a
sizeable fraction of every access pays the remote penalty.

Binding one job to a single NUMA domain, at 48 threads (`run_dada` seconds):

| configuration | ITS2 R1 | ITS2 R2 | 16S R1 |
|---|---|---|---|
| one job, interleaved (default) | 197.8 | 259.7 | 708.1 |
| **one job, bound to one domain** | **143.9** | **185.2** | **528.1** |
| two jobs, interleaved | 195.2 / 193.8 | 258.7 / 258.8 | 656.0 / 651.8 |
| **two jobs, one domain each** | **139.0 / 138.1** | **179.3 / 178.1** | **527.8 / 540.4** |
| *binding gain* | *−27%* | *−29%* | *−25%* |

Throughput, two bound jobs against one default job: **2.6× (16S) to 2.9× (ITS2)**.

```bash
# one job, pinned to a domain
numactl --cpunodebind=0 --membind=0 dada2-rs dada-pooled --threads 48 ...

# two jobs, one per domain — ~2.9x the throughput of a single default job
numactl --cpunodebind=0 --membind=0 dada2-rs dada-pooled --threads 48 ... &
numactl --cpunodebind=1 --membind=1 dada2-rs dada-pooled --threads 48 ... &
wait
```

Three things this measurement showed that are worth knowing:

1. **Almost all of the gain is locality, not packing.** A single bound job
   captures ~92% of what two bound jobs get. Pinning is worth doing even if you
   only ever run one sample at a time.
2. **The gain appears in every phase**, from the threaded inner loop to the
   single-threaded bookkeeping — on 16S: map −27%, store −22%, shuffle −21%,
   p-update −31%. Phases with completely different characters improving by
   similar amounts is the signature of memory *latency*, not bandwidth.
3. **Two unpinned jobs do not slow each other down** (0% on ITS2, and 16S was
   7.6% *faster* per job than running one alone — unexplained). Either way the
   machine was never short of bandwidth. It was waiting on distant memory.

**The thread count must fit inside one domain.** Our node has 64 cores per
domain, so 48 works and 96 cannot be bound this way. `numactl --hardware` shows
your layout. If your threads do not fit a domain, use `--interleave=all`
instead — binding a domain-sized-plus thread count onto one domain makes things
worse, which is how our own benchmarking guidance came to recommend interleaving
in the first place.

## What we have not tested

- **Binned-quality data at scale** beyond the pools above (see
  [Binned quality scores](findings/binned-quality.md)).
- **Thread scaling under NUMA binding** — the thread table was measured
  interleaved, so the knees may shift once memory is pinned.
- **Anything but x86-64 servers.** No Apple silicon, no ARM server parts, no
  single-domain workstations, where none of the NUMA section applies.
- **`dada` and `dada-pseudo`** — everything here is `dada-pooled`, which is the
  mode with the largest single-process footprint.

If you measure your own data, the tooling is in
[Diagnostics](diagnostics.md#-verbose-progress-lines-the-shape-of-a-run-not-just-its-mean)
and `dev/run_concurrency_test.sh`. Results that contradict this page are welcome
— it is a summary of four arms, not a law.
