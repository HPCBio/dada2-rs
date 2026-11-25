# Batched alignment: 16 pairs per SIMD vector

*Issue [#337](https://github.com/HPCBio/dada2-rs/issues/337). EPYC 7713,
48 threads, `release-native` build unless stated otherwise, three replicates
per arm. Output byte-identical in every run.*

**Verdict:** aligning the pairs that pass the k-mer screen 16 at a time, one
pair per SIMD lane, never changes output and never measurably costs time.
`--batch-align` is now on by default.

**Results:** pooled runs at 48 threads, median of three:

- −16% on MiSeq.
- −28% on PacBio at k = 7, and −42% on a build without native CPU flags.
- −63% on PacBio at k = 5.
- −53% on MiSeq `dada-pseudo`.
- No change on NovaSeq, where the k-mer screen dominates.

**What sets the gain:** batching makes alignment about 2–3× faster end to end,
so a run gains in proportion to alignment's share of its comparison time.

The page follows the work in order: the kernel, then wiring it into `b_compare`,
then three regressions found on the cluster and their fixes, then the results.

## Starting point

The banded Needleman–Wunsch kernel is the largest single cost in pooled `dada`:
44%, 30% and 63% of `run_dada` on MiSeq R1, MiSeq R2 and PacBio
([screen vs align](compare-screen-vs-align.md)). It ran at 1.2–1.8 ns per DP
cell, well below what an `i16` SIMD kernel should reach.

The kernel vectorizes *within* one pair, along anti-diagonals. At band 16 an
anti-diagonal holds about 17 cells, so the per-diagonal overhead dominates.

The idea came from an outside view. A design discussion with no access to the
code surveyed alternative backends: difference recurrences, adaptive bands,
bit-parallel edit distance, and anchored DP. It ranked *inter-sequence*
batching first, because it leaves every alignment unchanged.

That fits how `b_compare` works: it aligns one cluster centre against every
unique. So many pairs share one sequence, and pairs of equal length share the
band's whole control flow. Earlier tuning had only looked inside the per-pair
kernel.

## The kernel

`src/nwalign_batch.rs` runs the per-pair algorithm (`align_vectorized_with_buf`,
R's `nwalign_vectorized2`) on 16 pairs of identical shape at once. DP storage is
lane-interleaved.

One layout choice decided the result. The first sequence is stored *reversed*,
so the anti-diagonal fill becomes one flat loop. Without it the 16-lane kernel
was only 1.17× faster than the per-pair one on a laptop; with it, 2.9×.

Output is byte-identical by construction. That is checked by fuzzing (4,500 random
batches, every lane against the per-pair kernel), by every real pair in the benchmark, and by a parity test on the
full `b_compare` state.

**Kernel benchmark** (`examples/batch_nw_bench.rs`): full pools, 100 centres,
NUMA-pinned. 16-lane speedup per DP cell over the per-pair kernel:

| Build | Pool | 1 thread | 24 threads | 48 threads |
|---|---|---|---|---|
| `release` | MiSeq (362) | 4.94× | 4.47× | 4.53× |
| `release` | NovaSeq 16S | 4.98× | 5.05× | 5.57× |
| `release` | PacBio 95 | 3.04× | 3.18× | 3.08× |
| `release-native` | MiSeq (362) | 4.14× | 3.58× | 3.70× |
| `release-native` | NovaSeq 16S | 4.14× | 3.93× | 4.35× |
| `release-native` | PacBio 95 | 3.11× | 3.39× | 2.98× |

- **The first cluster run "collapsed" at 16+ threads.** The cause was the
  harness, not the kernel: each batch moved its lane vectors out, which meant
  one allocation per pair. A control arm that allocated on purpose collapsed
  the same way (0.11–0.18×), which confirmed it.
- **32 lanes lost at 48 threads,** where production runs, though they were
  3–7% faster at low thread counts. So 16 lanes it is.
- **Traceback takes 25–40% of batched time.** Two attempts to cut it:
  - Skipping the pointer-matrix zero-fill saved 6–7% on Illumina.
  - Interleaving the 16 lanes' tracebacks was 23–39% slower, and was reverted.
- **PacBio gains less** (about 3× against 4–5.5×). Its pointer matrix is about
  1.6 MB per batch, which exceeds Zen 3's 512 KB L2.

## Wiring it into `b_compare`

Under `--batch-align`, `b_compare_parallel` runs in two passes:

1. **Screen pass.** It handles every unique as before, except that a pair bound
   for the vectorized kernel is *deferred*.
2. **Align pass.** The deferred pairs are sorted by length and aligned 16 at a
   time. Groups smaller than 8 use the per-pair kernel. λ and the stored
   comparison are then computed exactly as before.

Some pairs never reach the vectorized kernel, so they keep their existing
paths: the gapless shortcut, WFA, affine gaps, `--homo-gap-p`, and reads over
3,500 bp. `--screen-audit` turns batching off. A single-threaded run takes the
serial `b_compare` and never batches.

Grouping by length fills batches well. 99.4% of PacBio pairs land in full
16-lane batches, and the mean fill is 15.9–16.0 lanes on every pool.

On a laptop (6 threads, subsets), runs were 2.4–3.2× faster:

| Run | Per-pair | Batched |
|---|---|---|
| MiSeq `learn-errors` | 4.94 s | 1.67 s |
| MiSeq `dada-pooled` | 1.72 s | 0.69 s |
| NovaSeq `dada-pooled` | 99.5 s | 41.8 s |
| PacBio `learn-errors` | 89.0 s | 33.0 s |
| PacBio `dada-pooled` | 52.4 s | 16.2 s |

## Three regressions on the cluster

The first full-pool A/B, on MiSeq 362, gave identical outputs but disappointing
times. `learn-errors` was **35% slower** (11.7 → 15.7 s) and `dada-pooled` only
3% faster. The laptop had shown none of this. Each fix below was found by
measurement on the node.

**Per-cluster buffers.** The screen pass built a fresh batch buffer per rayon
job, hundreds per cluster comparison (about 260 at 10 threads on a laptop), and
each one zero-filled 0.15–1.6 MB. At
48 threads that churned the allocator. Thread-local scratch fixed it:

- `learn-errors` back to level with per-pair.
- `dada-pooled`'s `dada` phase from 71.8 to 59.2 s.

**A lost overlap.** In the per-pair path each unique is screened and aligned
together, so alignment compute hides the screen's memory stalls. Split into two
passes, the screen pass streams k-mer vectors alone, and it stops scaling past
about 24 threads. That is limited by memory bandwidth and has no cheap fix:

- Aligning each rayon job's deferred pairs in place would restore the overlap.
  It was rejected because per-job pools are too small to fill PacBio batches
  (about 120 pairs over about 60 lengths).
- Collecting the deferred pairs inside the screen pass with a rayon `fold` made
  that pass 2.7× slower at 48 threads (22 → 60 s). It was reverted.

**Join-point stalls.** After the screen pass, a "prep" step finds the deferred
pairs and sorts them, and it ran as two small parallel steps. A first guess was
that the cost came from re-reading the 15 MB results array. Collecting the
lengths into a compact vector instead cut prep only from 8.1 to 6.8 s, so that
guess was mostly wrong.

The real cause was rayon's join points. `learn-errors` and `dada-pseudo` run
samples concurrently, each with its own `b_compare` inside, so the scheduler
always has other work queued. A thread waiting at a join picks up another
sample's work, and its own sample's sequential chain stalls.

Per-pass `learn-errors` metrics ([#338](https://github.com/HPCBio/dada2-rs/issues/338))
measured it: prep took 218 thread-seconds over 30,749 calls, about 7 ms per call
for roughly 2,500 uniques. Running prep serially fixed it:

- MiSeq `learn-errors` from 13.8 to 9.6 s.
- MiSeq `dada-pooled` from 75.3 to 70.5 s.

A tuning gate then compared serial and parallel prep on every platform. Only
NovaSeq `dada-pooled`, with 2.3M uniques, preferred parallel, by 1–2%. Serial
prep was better by 8–10% on NovaSeq and PacBio `learn-errors`, and by 30% on
MiSeq. Prep is now always serial and the gate is gone. A size threshold is left
as an option in [#342](https://github.com/HPCBio/dada2-rs/issues/342).

## Results

Wall time, median of three, 48 threads. In every row the error models matched
by `cmp` and the denoised outputs by `diff -r`.

| Run | Denoising wall, per-pair | batched | change | `learn-errors` per-pair | batched | change |
|---|---|---|---|---|---|---|
| MiSeq 362, `dada-pooled` | 84.0 s | 70.5 s | −16% | 13.9 s | 9.6 s | −31% |
| MiSeq 362, `dada-pseudo` | 35.6 s | 16.7 s | −53% | 11.7 s | 8.5 s | −28% |
| PacBio 95 k = 7, `dada-pooled` | 394.9 s | 284.9 s | −28% | 14.4 s | 11.7 s | −19% |
| PacBio 95 k = 7, `release` build | 550.7 s | 319.8 s | −42% | 18.8 s | 13.4 s | −29% |
| PacBio 95 k = 7, single pass (`pvalue`) | 400.1 s | 281.5 s | −30% | – | – | – |
| PacBio 95 k = 5, `dada-pooled` | 2,056 s | 770 s | −63% | 38.5 s | 17.7 s | −54% |
| NovaSeq 16S, `dada-pooled` | 2,006 s | 2,033 s | +1% | 103.7 s | 104.3 s | tie |

- **The two MiSeq rows ran on different days.** Their `learn-errors` numbers
  are not comparable with each other: compare the two arms within a row, not
  across rows.
- **`release` build.** The default build, which the bioconda and crates.io
  binaries use, gains the most. Without AVX2 the per-pair kernel slows down more
  than the batched one: compare took 1.5× longer in the per-pair arm, against
  1.18× in the batched arm.
- **Single pass.** `--partition single-pass` with the p-value rule goes through
  the same compare, so it gains the same.
- **Not yet covered.** The UNOISE rule (`--birth-rule unoise`) has its own,
  per-pair compare and does not batch yet
  ([#345](https://github.com/HPCBio/dada2-rs/issues/345)).
- **Memory** was unchanged in every run.

## Why the gains differ

How much of the pairs the k-mer screen passes doesn't predict the gain by
itself, because a 1,450 bp alignment costs far more than a 250 bp one. What does
predict it is **alignment's share of the per-pair run's compare time**. That
share is estimated below as the per-pair map time minus the batched screen pass.

| `dada-pooled` | Aligned / screened | Alignment share of compare | Change |
|---|---|---|---|
| NovaSeq 16S, k = 5 (2.3M uniques) | 1.7% | 11% | +1% |
| MiSeq 362, k = 5 | 25% | 50% | −16% |
| PacBio 95, k = 7 | 12.2% | 61% | −28% |
| PacBio 95, k = 5 | 98.9% | ~98% | −63% |

Share of wall × (1 − 1/alignment speedup) reproduces the measured changes.

The pools were chosen to cover both ends. NovaSeq 16S is screen-heavy and
alignment-light; PacBio at k = 5 is the opposite, since the screen passes nearly
every pair. Other pools should fall between them, by their alignment share. The
NovaSeq ITS2 pool was not run.

On NovaSeq the ceiling is about 9% of wall, and serial prep (64–72 s per run)
takes back what the align pass saves. There, the k-mer screen (1,006 s) and the
store (303 s) are the costs to attack, not alignment.

## Consequences

- **`--batch-align` is on by default** on `dada`, `dada-pooled`, `dada-pseudo`,
  `learn-errors` and `errors-from-sample`. `--batch-align false` restores
  per-pair alignment. A test runs every mode both ways at 4 threads and
  requires identical output.
- **Prep stays serial.** Do not parallelise small steps inside `b_compare`
  without measuring under nested sample parallelism. The laptop cannot show the
  stall.
- **Ruled out:**
  - 32 lanes.
  - Interleaving the lanes' tracebacks.
  - Aligning deferred pairs inside each screen-pass job.
  - Collecting deferred pairs with a `fold`.
  - Parallel prep in general.
- **Next, if anything:** [#342](https://github.com/HPCBio/dada2-rs/issues/342)
  lists the remaining candidates:
  - overlapping the screen and align passes;
  - packing traceback pointers into 2 bits, for PacBio;
  - producing substitutions during traceback.

  The UNOISE port is [#345](https://github.com/HPCBio/dada2-rs/issues/345).
  AVX2 dispatch for the default build
  ([#336](https://github.com/HPCBio/dada2-rs/issues/336)) would be worth about
  12% more on batched PacBio.
