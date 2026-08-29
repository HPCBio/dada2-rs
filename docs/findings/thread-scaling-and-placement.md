# Thread scaling and memory placement

**Two results with opposite characters.** How many threads to use is strongly
data-dependent, and predictable from a number the verbose log already prints.
Where to put the memory is *not* data-dependent, is worth 25–29% on every arm
measured, and reverses a recommendation this project made from its own earlier
data.

## Why the sweep happened

The [progress lines](../diagnostics.md#-verbose-progress-lines-the-shape-of-a-run-not-just-its-mean)
added in #150 showed pooled `run_dada` is not a uniform workload: occupancy
climbs from ~22 to ~42 effective cores of 64 across a single run, so the
end-of-run mean of ~29 describes no part of it. The ramp is the **serial fraction
shrinking** — `b_shuffle` plus `p_update` fall from 15.7 s to 5.7 s per 30-second
window — not the workload changing.

That made every thread-count figure we had suspect, because all of them were
run-level means over a run spanning nearly a factor of two.

## Threads: data-dependent, and predictable

Two replicates per point, `run_dada` seconds, 128-core EPYC 7713:

| threads | 16S R1 | 16S R2 | ITS2 R1 | ITS2 R2 |
|---|---|---|---|---|
| 24 | 1072.5 | 1076.3 | 221.1 | 308.5 |
| 48 | 713.8 | 767.8 | **194.0** | **254.9** |
| 64 | 638.9 | **711.5** | 198.2 | 258.2 |
| 96 | **594.3** | 734.9 | 203.8 | 263.2 |

**The best thread count differs two-fold between two reads of the same pool.**
And it is predicted by the map's screen/align split, which `--verbose` reports:

| arm | k-mer screen (% of map busy) | best threads |
|---|---|---|
| 16S R1 | 50.8% | 96 (still scaling) |
| 16S R2 | 67.7% | 64 |
| ITS2 R2 | 83.2% | 48 |
| ITS2 R1 | 85.8% | 48 |

The screen streams per-sequence k-mer vectors and is memory-bound; the aligner
does dense DP and is core-bound. A screen-dominated pool saturates early. The
ordering is monotonic across four arms and **predicts within a pool as well as
across pools** — 16S's two reads differ by 17 points of screen share and have
different optima — so it tracks arithmetic intensity, not the dataset's name.

Two axes move it, and both are visible here. **Amplicon and diversity dominate**:
these pools came from the *same run and the same samples*, so the 16S/ITS2 gap is
biology, not instrument. **Read position matters too**: same library, same
organisms, 17 points of difference from quality decay along read 2.

### Threads buy latency, not throughput

Parallel efficiency relative to 24 threads: 75%/63%/45% (16S R1 at 48/64/96) and
57%/42%/27% (ITS2 R1). **Every arm buys wall time at roughly 4× the
core-seconds.** So the two questions have opposite answers — one urgent sample
wants the table above; a queue of samples wants fewer threads and more jobs.

## Placement: not data-dependent, and worth 25–29%

Binding one job to a single NUMA domain at 48 threads:

| | ITS2 R1 | ITS2 R2 | 16S R1 | 16S R2 |
|---|---|---|---|---|
| one job, `--interleave=all` | 197.8 | 259.7 | 708.1 | 785.5 |
| **one job, `--cpunodebind=0 --membind=0`** | **143.9** | **185.2** | **528.1** | **575.0** |
| two jobs, interleaved | 195.2 / 193.8 | 258.7 / 258.8 | 656.0 / 651.8 | 826.0 / 834.5 |
| **two jobs, one domain each** | **139.0 / 138.1** | **179.3 / 178.1** | **527.8 / 540.4** | **530.4 / 559.4** |
| *binding gain* | *−27%* | *−29%* | *−25%* | *−27%* |

A 25–29% band across a **35-point gap in screen share**. We predicted the
compute-bound pool would gain much less, since arithmetic can hide memory
latency. It did not.

Three things the arms settle:

1. **It is locality, not packing.** A single bound job captures ~92% of what two
   bound jobs get. The `sbind` arm exists solely to separate these, because
   `split` confounds them.
2. **It is latency, not bandwidth.** Every phase gains by a similar amount
   despite completely different characters — 16S: map −27%, store −22%, shuffle
   −21%, p-update −31%. And two *unpinned* jobs barely affect each other (−8% to
   +6% across arms), which a bandwidth shortage would not permit.
3. **Unpinned performance is a lottery.** That −8%-to-+6% spread is the arm where
   placement is left to the kernel; where pages land decides the run.

Two bound jobs give **2.6–2.9× the throughput** of one default job, and better
per-job latency at the same time. All 600 outputs across all arms byte-identical.

## Reconciling this with "binding is the obvious fix that turns out to be wrong"

[Measuring on a NUMA node](measuring-on-numa.md) concluded the opposite, and
`dev/numa_pin.sh` still recommends interleaving on its authority. Both are right,
and **the earlier page already contained the mechanism** — nobody read it that
way, this author included:

| PacBio, 64 threads | bound | interleaved |
|---|---|---|
| parallel map | 286.8 s | **235.4 s** |
| `b_shuffle` (serial) | **50.4 s** | 65.3 s |
| `store` (serial) | **12.8 s** | 20.2 s |
| `run_dada` | 391.4 s | **372.5 s** |

The serial phases were already **23% and 37% faster bound** — within a point or
two of what we just measured on entirely different pools. The map penalty simply
outweighed them, so the run-level total said "binding loses" and the phase table
underneath went unexamined.

What changes with 48 threads is the map penalty. Sixty-four threads bound to a
64-core domain draws all its bandwidth from one node's controllers; 48 leaves
headroom. The serial gain is unchanged, the map penalty shrinks, and the sign
flips.

**Caveat, and it is a real one:** the old measurement is PacBio at 64 threads and
the new one is NovaSeq at 48, so platform and thread count are confounded. The
reconciliation is *consistent* with thread-count-relative-to-domain-size, and the
matching serial numbers are good evidence for it, but it is not isolated. The
clean test is one pool bound at 48 **and** 64 — worth doing before the default in
`dev/numa_pin.sh` changes.

## What this dictates

- **`dev/numa_pin.sh` is annotated, not changed.** It exists for
  *reproducibility*, and whether binding replicates as tightly as interleaving is
  unmeasured. But "interleave beats bind" is no longer safe as a general claim.
- **Published absolute numbers are conservative.** Everything in
  [Results](../results.md) was gathered interleaved at sub-domain thread counts.
  A/B *deltas* are unaffected — both arms always shared the policy.
- **The benchmark harness now applies placement symmetrically** (`--numa`), since
  pinning one stack and not the other would inflate every speedup by ~25%.
  Whether R DADA2 gains as much is [#155](https://github.com/HPCBio/dada2-rs/issues/155).
- **The thread table above was measured interleaved**, so the knees may shift once
  memory is pinned. Not yet retested.
- User-facing guidance: [Tuning for your data](../tuning-for-your-data.md).

## Instrument errors, and what each cost

Four, all of which produced confident and wrong readings:

1. **Comparing fixed wall-clock windows across thread counts.** At the same `t=`,
   a 96-thread run and a 48-thread run are at different cluster positions in a
   workload whose serial fraction collapses over the run. This reported the
   penalty as largest *early*; aligned by cluster it is largest *late* — **the
   sign of the trend inverted**. Now encoded in `dev/analyze_thread_sweep.py`
   rather than remembered.
2. **Locating a knee from one side of it.** From 48/64/96 alone, ITS2 degraded
   monotonically and was reported as "already past the knee below 48". Adding 24
   showed 48 *is* the peak, beating 24 by 12–17%. A sweep that does not bracket
   the peak cannot find it.
3. **A predictor that is itself thread-dependent.** Screen share reads 80.0% at
   24 threads and 85.8% at 48 for the same arm, because the bandwidth-bound half
   degrades faster under contention. Adding 24 to ITS2 silently made its number
   incomparable to 16S's — inside the table presenting them as a calibration.
   Now pinned to a reference thread count.
4. **The right estimator on the wrong arm.** Reporting per-round minima defends
   against a noisy neighbour and was correct for five of six benchmark arms; for
   the one whose variance is *intrinsic* it selects against the effect under
   test. See [the store scan](compare-store-scan.md).

The pattern across all four: **each produced a plausible, quotable number, and
nothing in the output said the instrument could not answer the question.**
