# Benchmark results

Head-to-head performance vs R DADA2 on representative datasets. Numbers come from
the benchmark harness (see [Tooling & metrics](benchmarking.md)) run on our
cluster — the datasets are large, so these are run manually and the tables below
are regenerated from each run's `summary.csv` with
[`bench_table.py`](#regenerating-these-tables).

## Methodology

- **Function-vs-function.** Each mode compares the dada2-rs subcommand against
  its direct R analog, one process each: `dada-pooled` ↔ `pool=TRUE`, multi-input
  `dada` ↔ `pool=FALSE`, `dada-pseudo` ↔ `pool="pseudo"`.
- **Wall** is the fair end-to-end time (R as a single process). **Speedup** =
  R wall ÷ dada2-rs wall.
- **Peak RSS** is the per-process high-water mark. The ratio is reported as
  **dada2-rs ÷ R**, so it reads directly as a fraction of the comparator's
  memory: **<1× means dada2-rs uses *less*** (e.g. 0.7× = 70% of R's RSS, a 30%
  reduction) and **>1× means more** (2.0× = twice R's). This fraction-of-the-
  comparator convention is the one commonly used in memory benchmarks to surface
  a reduction directly, and it handles both wins and increases in one column — at
  the cost of running in the *opposite* direction to Speedup (where higher is
  better). **Bold** marks the rows where dada2-rs wins the axis.
- **Build:** `release-native` unless noted. See
  [build target](benchmarking.md#6-build-target-matters).
- **Correctness** (ASV concordance) is validated separately — see
  [concordance tooling](benchmarking.md#5-concordance-validation-tooling).

!!! note "Populate from a cluster run"
    We generate the results using `dev/benchmark/bench_pooled.py` 
    and summarize with:
    ```bash
    python3 dev/benchmark/bench_table.py \
        "pooled=bench_true/summary.csv" \
        "per-sample=bench_false/summary.csv" \
        "pseudo=bench_pseudo/summary.csv"
    ```
    then paste the table below.

## PacBio HiFi

* `dada2-rs v0.1.1-9a0c5da3` (v0.1.1 prerelease, `release-native` build)
* Data is from Hergenrother 2024, 93 total samples, PacBio Sequel IIe.
* Default settings using FASTQ input for `dada` in both tools
* 24 threads, 96GB memory allocated on an exclusive node run

### Overall workflow

PacBio is benchmarked from `--kmer-size 5` (matched to R's fixed
`KMER_SIZE = 5` — isolates kernel + threading speed) to `--kmer-size 7`
(dada2-rs recommended default for PacBio — adds the k-mer screen-effectiveness gain). See the [PacBio notes](benchmarking.md#7-pacbio-vs-illumina-specifics).

*Overall workflow (remove-primers FASTQ -> removing chimeras)*

| Run | Pooling mode | k-mer | sample jobs | dada2-rs wall (s) | R wall (s) | Speedup (R÷rs) | dada2-rs peak (MB) | R peak (MB) | Peak RSS (rs÷R) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PacBio_pooled_k5    | full | 5 | NA | 4336.8 | 11885.4 | **2.7×** | 29469 | 40062 | **0.7×** |
| PacBio_pooled_k6    | full | 6 | NA | 1052.0 | 11885.4* | **11.3x** | 34365 | 40062* | **0.9x** |
| PacBio_pooled_k7    | full | 7 | NA | 761.1 | 11885.4* | **15.6x** | 53781 | 40062* | *1.3x* |
| PacBio_pseudo_k5    | pseudo | 5 | 6 | 2080.3 | 9299.1 | **4.5×** | 6173 | 3081 | *2.0×* |
| PacBio_pseudo_k6    | pseudo | 6 | 6 | 496.5 | 9299.1* | **18.7x** | 6877 | 3081* | *2.2x* |
| PacBio_pseudo_k7    | pseudo | 7 | 6 | 381.2 | 9299.1* | **24.4x** | 10907 | 3081* | *3.5x* |
| PacBio_nopool_k5    | none | 5 | 6 | 859.0 | 7219.1 | **8.4×** | 5873 | 2783 | *2.1×* |
| PacBio_nopool_k6    | none | 6 | 6 | 223.1 | 7219.1* | **32.4x** | 6206 | 2783* | *2.2x* |
| PacBio_nopool_k7    | none | 7 | 6 | 180.3 | 7219.1* | **40.0x**| 8640 | 2783* | *3.1x* |
| PacBio_pseudo_k5_j1 | pseudo | 5 | 1 | 2207.9 | 9299.1* | **4.2x** | 2880 | 3081* | **0.9x** |
| PacBio_pseudo_k6_j1 | pseudo | 6 | 1 | 642.0 | 9299.1* | **14.5x** | 2801 | 3081* | **0.9x** |
| PacBio_pseudo_k7_j1 | pseudo | 7 | 1 | 555.0 | 9299.1* | **16.8x** | 3717 | 3081* | *1.2x* |
| PacBio_nopool_k5_j1 | none | 5 | 1 | 914.5 | 7219.1* | **7.9x** | 2280 | 2783* | **0.8x** |
| PacBio_nopool_k6_j1 | none | 6 | 1 | 288.4 | 7219.1* | **25.0x** | 2158 | 2783* | **0.8x** |
| PacBio_nopool_k7_j1 | none | 7 | 1 | 262.7 | 7219.1* | **27.5x** | 2733 | 2783* | *1.0x* |

`*` - values from k=5 run. 

### Preliminary observation and results

**Recommendation:** use `--kmer-size 7` for PacBio runs, especially if your compute has a decent memory footprint. If you need to reduce memory and try to retain throughput, use `--kmer-size 7` and `--sample-jobs 1` for PacBio runs

#### k-mer size

- R DADA2 *currently* has a fixed (hard-coded) k-mer size of 5 for screening. This appears to have the effect of reducing k-mer screening efficiency almost to zero (a no-op), with almost every unique sequence proceeding to NW alignment. 
- Increasing the k-mer length marginally has a significant effect, dramatically reducing walltime. 
- Preliminary followup (to be added) indicates this has essentially **no impact** on the ASVs derived post-chimera removal as well as their counts, and may actually increase sensitivity slightly with a few more ASVs recovered in the `dada2-rs` runs.
- We are also evaluating the overall effect increasing k-mer sizes have on clustering with the different pooling results. Initial k=8 results confirms this increases memory usage dramatically (>100GB), so we don't anticipate this to have much practical use. 

#### ASVs

- Pre-chimera results differ slightly from R DADA2 (with some ASVs found unique to both runs). However as noted above these are essentially removed with chimera removal, suggesting their lower abundance

#### Memory vs walltime

- `dada2-rs` currently uses more memory by default (last column) under default settings, particular with higher k-mer settings. This is at the `dada-pooled`/`dada-pseudo`/`dada` stage (see below)
- `dada-pseudo` and `dada` (per-sample) memory overhead can be alleviated to run each sample serially by setting the number of concurrent jobs to one (`--sample-jobs 1`), with a small walltime penalty. 
- We're looking into potential paths to further reduce this footprint

### Per step comparisons

**PacBio pooled, k=7**  — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| remove_primers | 19.0 | 20.9 | 31 | 5410.2 | **284.3×** |
| learn | 25.2 | 20.6 | 2791 | 159.7 | **6.3×** |
| dada | 698.8 | 18.8 | 53781 | 6226.2 | **8.9×** |
| make_table | 0.3 | 1.0 | 111 | 0.1 | 0.2× |
| remove_bimera | 17.8 | 23.6 | 103 | 19.8 | 1.1× |
| TOTAL | 761.1 | 19.0 | 53781 | 11885.4 | **15.6×** |

**PacBio pseudo, k=7**  — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| remove_primers | 19.0 | 20.8 | 29 | 5380.1 | **283.6×** |
| learn | 25.1 | 20.6 | 2809 | 155.2 | **6.2×** |
| dada | 329.5 | 21.7 | 10907 | 3688.6 | **11.2×** |
| make_table | 0.2 | 1.0 | 75 | 0.1 | 0.3× |
| remove_bimera | 7.4 | 23.5 | 85 | 8.3 | 1.1× |
| TOTAL | 381.2 | 21.6 | 10907 | 9299.1 | **24.4×** |

**PacBio, no pooling, k=7** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| remove_primers | 19.1 | 20.8 | 29 | 5445.1 | **285.7×** |
| learn | 25.1 | 20.5 | 2717 | 159.9 | **6.4×** |
| dada | 131.9 | 21.0 | 8640 | 1542.7 | **11.7×** |
| make_table | 0.1 | 1.0 | 54 | 0.0 | 0.3× |
| remove_bimera | 4.0 | 23.3 | 69 | 4.5 | 1.1× |
| TOTAL | 180.3 | 21.0 | 8640 | 7219.1 | **40.0×** |

Note that most of the walltime improvement is actually dominated by the `remove-primers` step (which combines two steps from DADA2: `removePrimers` and `filterAndTrim`, and parallelizes both steps; R `removePrimers` processes data serially). In many cases alternative tools like `cutadapt` can help close this gap and may prove more flexible for others. However significant improvements are also apparent for both the `learn-errors` and `dada` steps for each pooling mode.

### Streaming vs cached (`--cache-samples`)

`dada-pseudo` streams samples by default (re-reading each sample per round) rather
than holding every sample's dereplicated state resident. `--cache-samples` opts
into the older all-in-memory behavior. This characterizes that tradeoff on PacBio
(k=5, node-exclusive single rep), comparing streaming (default) as baseline
against cached and with R. `--sample-jobs 1` is used to make R comparisons more apples-to-apples. Numbers are filtered per stack with`compare_bench.py --stack <stack>`; only the `dada` step moves materially (other
steps run on identical inputs).

| Stack | mode | dada wall (s) | dada peak | Δ wall | Δ peak |
|---|---|---:|---:|---:|---:|
| dada2-rs | streaming (default) | 2092.4 | 2.34 GB | — | — |
| dada2-rs | cached | 2033.5 | 8.18 GB | **−2.8%** | **+249%** |
| R-split | streaming (default) | 3737.8 | 2.86 GB | — | — |
| R-split | cached | 3454.6 | 14.86 GB | **−7.6%** | **+420%** |
| R-single | streaming (default) | 3735.1 | 3.02 GB† | — | — |
| R-single | cached | 3400.4 | 15.69 GB† | **−9.0%** | **+420%** |

`†` R-single runs the whole pipeline as one process, so peak RSS is only resolvable
process-wide (not per step); it reflects the `dada` phase, which dominates.

**Verdict:** caching buys a small walltime gain on the `dada` step at a large peak-memory
cost, and the pattern holds across all three stacks. Streaming is the right
default; `--cache-samples` stays an opt-in for the walltime-bound, memory-rich case.

Two cross-stack notes: 
- streaming-vs-streaming: `dada2-rs` is leaner and faster on `dada` (2.34 GB / 2092 s vs 2.86 GB / 3738 s)
- cached-vs-cached, dada2-rs is ~1.8× leaner (8.18 GB vs 14.86 GB) — R's resident-derep cache is not magically compact.

See [issue #22](https://github.com/HPCBio/dada2-rs/issues/22) for the keep-but-default-off
decision.

### Pre/post A/B: dropping the resident 16-bit k-mer vector ([#32](https://github.com/HPCBio/dada2-rs/issues/32))

A dada2-rs-vs-itself check isolating the [#32](https://github.com/HPCBio/dada2-rs/issues/32)
change (PR #34): **`pre`** is `main` immediately before it, **`post`** drops the
resident `4^k` u16 k-mer frequency vector from each `Raw` (recomputed only on the
near-never `kmer_dist8` overflow). Full 93-sample PacBio dataset, node-exclusive,
no R. **Output is byte-identical** — `compare_asvs.py` pre-vs-post on the final
chimera-filtered table: 2082 ASVs vs 2082, churn 0.

Only the `dada` step moves materially (`remove_primers` runs on identical inputs;
`learn` sheds the same per-`Raw` vector and so also drops ~−35%):

| mode | k | dada peak (pre → post) | peak Δ | dada wall Δ | dada cores Δ |
|---|---|---|---:|---:|---:|
| pooled (`--pool true`) | 5 | 28.72 → 27.59 GB | −3.9% | −0.2% (flat) | flat |
| pooled (`--pool true`) | 7 | 52.33 → 35.69 GB | **−31.8%** | **−5.6%** | +3.4% |
| pseudo cached (`--cache-samples`) | 7 | 14.29 → 10.27 GB | **−28.1%** | **−4.5%** | +2.7% |

**The memory win scales with `4^k`.** The dropped vector is `4^k × 2` bytes/unique
— 2 KB at k5, 32 KB at k7. Pooled sheds −1.1 GB at k5 vs −16.6 GB at k7, a 14.7×
ratio ≈ the theoretical 16×, so the reduction is provably exactly the removed
vector and nothing else. The `4^k` k-mer arrays are the dominant resident term in
the all-samples modes — pooled `dada` is strongly k-related (k7 ~52 GB vs k5 ~29
GB), not k-independent.

**The wall bonus is memory-bandwidth relief**, and it appears only where the
resident set was large enough to saturate bandwidth. At k7, pooled and cached had
depressed cores (18.2, 19.8) that recover (+3.4%, +2.7%) once ~16 GB / ~4 GB of
k-mer data is shed, cutting `dada` wall −5.6% / −4.5%. At k5 pooled, cores were
already 22.8/24 (below the bandwidth ceiling) and the wall is flat — nothing to
recover. The win improving *only* where cores were depressed (not where they were
saturated) rules out a build-time CPU saving, which would show flat at both k.
This partially addresses [#33](https://github.com/HPCBio/dada2-rs/issues/33).
Streaming pseudo is not shown: its small per-sample working set holds little
resident k-mer data, so there is no u16 term to shed (≈flat).

## Illumina MiSeq (F1000, 384 samples)

The majority of recent memory work focused on PacBio data. How does this affect Illumina data processing? 

* `dada2-rs v0.1.1-9a0c5da3` (v0.1.1 prerelease, `release-native` build)
* Data is the F1000 'MiSeqSOP' full data set: 384 samples, MiSeq v2 run.
* Default settings using FASTQ input for `dada` in both tools
* 24 threads, 96GB memory allocated on an exclusive node run

*Overall workflow (filter and trim FASTQ -> removing chimeras)*

| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (rs÷R) |
|---|---:|---:|---:|---:|---:|---:|
| MiSeqSOP_pooled_k5 | 417.8 | 1284.6  | **3.1×** | 3653 | 6105  | **0.6×** |
| MiSeqSOP_pooled_k6 | 414.0 | 1284.6* | **3.1x** | 3666 | 6105* | **0.6x** |
| MiSeqSOP_pooled_k7 | 416.1 | 1284.6* | **3.1x** | 3660 | 6105* | **0.6x** |
| MiSeqSOP_pseudo_k5 | 165.6 | 1224.1  | **7.4×** | 1483 | 2001  | **0.7×** |
| MiSeqSOP_pseudo_k6 | 165.5 | 1224.1* | **7.4×** | 1517 | 2001* | **0.8x** |
| MiSeqSOP_pseudo_k7 | 164.8 | 1224.1* | **7.4×** | 1424 | 2001* | **0.7x** |
| MiSeqSOP_nopool_k5 | 93.2  | 791.3   | **8.5×** | 1554 | 1673  | **0.9×** |
| MiSeqSOP_nopool_k6 | 94.7  | 791.3   | **8.4x** | 1456 | 1673* | **0.9×** |
| MiSeqSOP_nopool_k7 | 94.2  | 791.3   | **8.4x** | 1414 | 1673* | **0.9×** |

`*` - values from k=5 run

Tests prior to memory updates were from 2.4x (pooled, k=5) to 5.9x (no pooling, k=5). This round suggests a modest improvement: 3.1x to 8.5x, all at k=5. 

Prior data also suggested a small improvement for pooled samples using a higher k-mer, but in this latest run it seems to disappear, replaced with an overall improvement in performance and memory. 

### Per step comparisons

**MiSeqSOP pooled, k=5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 11.4 | 20.3 | 21 | 30.8 | **2.7×** |
| learn_fwd | 18.9 | 21.5 | 642 | 83.0 | **4.4×** |
| learn_rev | 23.0 | 20.7 | 770 | 107.8 | **4.7×** |
| dada_fwd | 196.1 | 12.9 | 3653 | 436.1 | **2.2×** |
| dada_rev | 157.8 | 9.9 | 3009 | 355.1 | **2.3×** |
| merge | 7.2 | 15.6 | 1652 | 245.8 | **34.1×** |
| make_table | 0.6 | 1.0 | 337 | 0.3 | 0.5× |
| remove_bimera | 2.8 | 23.1 | 47 | 2.6 | 0.9× |
| TOTAL | 417.8 | 12.9 | 3653 | 1284.6 | **3.1×** |

**MiSeqSOP, pseudo, k5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 11.3 | 20.4 | 21 | 31.7 | **2.8×** |
| learn_fwd | 18.7 | 21.6 | 648 | 83.2 | **4.4×** |
| learn_rev | 22.4 | 21.1 | 772 | 109.4 | **4.9×** |
| dada_fwd | 65.4 | 20.8 | 686 | 433.9 | **6.6×** |
| dada_rev | 43.1 | 19.9 | 568 | 332.9 | **7.7×** |
| merge | 4.1 | 17.1 | 1483 | 208.9 | **51.5×** |
| make_table | 0.2 | 1.0 | 144 | 0.1 | 0.6× |
| remove_bimera | 0.3 | 20.6 | 21 | 0.2 | 0.8× |
| TOTAL | 165.6 | 20.5 | 1483 | 1224.1 | **7.4×** |

**MiSeqSOP, no pooling, k5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 11.2 | 20.6 | 21 | 31.1 | **2.8×** |
| learn_fwd | 19.0 | 21.5 | 635 | 80.8 | **4.3×** |
| learn_rev | 22.4 | 21.0 | 769 | 107.6 | **4.8×** |
| dada_fwd | 22.2 | 19.9 | 573 | 196.4 | **8.8×** |
| dada_rev | 14.8 | 19.0 | 467 | 155.5 | **10.5×** |
| merge | 3.4 | 15.5 | 1554 | 198.5 | **59.1×** |
| make_table | 0.1 | 1.0 | 62 | 0.1 | 1.0× |
| remove_bimera | 0.1 | 17.8 | 21 | 0.1 | 0.8× |
| TOTAL | 93.2 | 20.2 | 1554 | 791.3 | **8.5×** |

Notably the big improvement is with `merge-pairs`, largely due to threading, though overall steps are generally faster.  One interesting point: the full pooling run suggests that the two `dada` steps are not fully utilizing all cores, which may be a further point to assess. 

### Pooled `dada` core under-utilization: root cause (Amdahl serial fraction)

The depressed core counts in the pooled `dada` steps (12.9 for `dada_fwd`, 9.9
for `dada_rev`, of 24 available) are **not** a parallel-efficiency problem — the
compare step's parallel alignment map is already 96–97% efficient at 24 cores.
The low averages are entirely the *serial* phases of `run_dada`
(`b_shuffle2` / `b_bud` / `b_p_update` / the compare store-loop) running at ~1
core and dragging the average down. Verbose phase timings from the pooled k5 run
(272,580 fwd / 296,898 rev raws):

| Phase | `dada_fwd` | `dada_rev` | Parallel? |
|---|---:|---:|---|
| compare-map | 98.9s | 58.1s | ✅ 24 cores @ 96–97% |
| compare-store | 12.7s | 11.1s | ❌ serial |
| **shuffle** (`b_shuffle2`) | **49.6s** | **57.8s** | ❌ serial |
| **bud** (`b_bud`) | **21.5s** | **18.9s** | ❌ serial |
| p_update (`b_p_update`) | 2.9s | 4.1s | ❌ serial |

Weighting these serial/parallel splits reproduces the reported 12.9 / 9.9
average core counts almost exactly, confirming the diagnosis. In `dada_rev`,
`shuffle` alone (57.8s) exceeds the entire parallel map (58.1s).

This is **pooled-specific**: `b_shuffle2` and `b_bud` both scan all
clusters × all raws every round. In no-pool / pseudo each sample is small so
these scans are trivial and cores stay ~20; pooling merges all 384 samples into
one `B` (~273k–297k raws, hundreds of clusters), so the O(nraw × nclusters)
serial scans balloon. Overall throughput is unaffected (still 3.1× vs R), but
there is significant headroom.

#### Attempted fix: parallelizing the serial scans (NEGATIVE RESULT — do not retry)

The obvious idea — thread-parallelize `b_shuffle2`'s per-raw best-cluster
selection (~50–58s serial) and `b_bud`'s min-p-value scan (~19–22s serial) —
was implemented byte-identical and A/B'd on this exact 384-sample pool (four
arms, one run each, all confirmed distinct commits via the bench `version.txt`):

| Metric (Δ% vs main) | `b_bud` only | `b_shuffle2` only | both |
|---|---:|---:|---:|
| learn_fwd wall | +1.7% | **+97.5%** | +102.0% |
| learn_rev wall | −0.1% | **+137.5%** | +121.4% |
| dada_fwd wall | +0.6% | **+62.7%** | +67.9% |
| dada_rev wall | +1.7% | **+80.7%** | +80.0% |
| dada_fwd cores | +3.2% | −12.1% | −13.1% |

Both changes were **net losses**:

- **`b_shuffle2` (atomic scatter-reduce) is catastrophic** — it doubles
  `learn-errors` and inflates pooled `dada` 60–80%, while *lowering* dada cores.
  It flattens all comparisons into a fresh allocation every call, then does
  random-access atomic CAS scatter (cache-hostile, contended on cluster 0's
  all-raw entries). Called thousands of times, this costs far more than the
  cheap sequential serial scan it replaced.
- **`b_bud` (lightweight per-cluster reduce) is a wash** — learn stays flat (no
  fan-out contention), dada cores tick up ~3%, but cpu rises ~4–5% and wall is
  unchanged (within run-to-run noise). The ~20s bud phase parallelized, but
  allocation/dispatch overhead ate the savings. No wall benefit.

**Why threading can't help here.** `learn-errors` runs the *same* dada algorithm
yet stays near-saturated (~21 cores) because it fans out across 384 independent
per-sample jobs (`all_inputs.par_iter()` in `learn_errors.rs`), so one sample's
serial phase overlaps another's parallel phase. Pooled `dada` is a single
monolithic run: its serial scans are a global barrier with no other work to fill
the cores, *and* the scans are cheap sequential/bandwidth-bound memory traffic,
not parallelizable compute. Spreading bandwidth-bound work across cores that
share one memory bus adds synchronization for no gain.

The remaining lever is therefore **algorithmic — reduce the scans' total work**,
not thread them: incremental/cached best-cluster tracking in `b_shuffle2` (only
re-score raws whose candidate clusters changed this round, vs rescanning all
clusters × all raws every shuffle) and an incremental structure for `b_bud`'s
minimum. The ~12-core pooled utilization is the inherent shape of a
single-population run, not reclaimable waste. This lever is pursued next, and
it works.

#### Incremental shuffle-to-convergence (WIN — merged approach)

First, the redundancy was quantified. Instrumenting the serial shuffle on the
384-sample pool (verbose `[dada] shuffle redundancy` line) showed the waste is
enormous: **10.3B / 10.8B comparisons scanned** (fwd/rev) to perform ~1.0M moves
— **~10,000 comps scanned per move** — and **31–35% of shuffle calls move zero
raws** (the terminal convergence-check of each bud's shuffle loop rescans
everything to move nothing). Within a shuffle loop the comparisons are *fixed*;
only cluster read counts change as raws move.

`b_shuffle_converge` exploits this: it maintains a persistent per-raw candidate
index (`CandIndex`, appended once per cluster — O(new comps) per bud, never
rebuilt), rebuilds the best-cluster map once per loop, then after each move pass
only recomputes raws whose candidate clusters' reads changed. It is
byte-identical to looping `b_shuffle2` (same max, same lowest-ci tie-break), so
moves are identical — confirmed by an order-independent ASV+count map match
against `main` at scale (the `seqtab.json` *row order* differs only by
per-process `HashMap` nondeterminism in `sequence_table.rs`, which varies across
any two runs including main-vs-main; concordance is checked on the ASV+count map,
not raw bytes).

The **cache-locality subtlety was decisive**. A first cut rebuilt the map each
loop by reading the raw-major inverted index — scattered, cache-missing access.
It cut comparisons by the predicted ~55% (10.3B → 4.5B) yet made the shuffle
phase *slower* (48→59s fwd): 4.5B scattered reads cost more wall than the serial
scan's 10.3B *sequential* reads — the same bandwidth/cache lesson as the
threading attempt. The fix (the merged version) does the per-loop initial build
the serial way — a contiguous cache-friendly scan of the per-cluster comp vecs
(the bulk of the work) — and uses the inverted index *only* for the reconcile,
where the touched-raw volume is small.

Result on the 384-sample pool (vs `main`, one run each, distinct commits):

| Metric (Δ% vs main) | scattered-index cut | **cache-friendly (merged)** |
|---|---:|---:|
| dada_fwd shuffle phase | 48→59s (+22%) | **48→32s (−33%)** |
| dada_rev shuffle phase | 56→60s (+6%) | **56→42s (−25%)** |
| dada_fwd wall | +4.2% | **−10.4%** |
| dada_rev wall | +3.5% | **−8.7%** |
| learn_fwd / learn_rev wall | ≈0 | **−2.9% / −2.6%** |
| **total pipeline wall** | +3.5% | **−8.4%** (427→391s) |
| dada peak RSS | +5% | +4.6–5.5% |

`learn-errors` speeds up too because it runs the same dada/shuffle internally.
The one cost is dada RSS +~5% (the index double-stores each comparison) — but
the *pipeline* peak RSS is unchanged, since the peak lives in `merge`, not
`dada`. Follow-ups (issue #85): apply the same incremental idea to `b_bud`
(done — below), and replace per-cluster `comp` storage *with* the index to erase
the duplication (deferred).

#### Incremental `b_bud` — per-cluster candidate cache (WIN)

With shuffle addressed, `b_bud` became the next-largest serial phase: it rescanned
every non-center raw across all clusters each call to find the minimum-p budding
candidate — **271k / 297k raws scanned per bud** (fwd/rev), ~19–22s serial total.

The redundancy gate here flipped the design. A full shuffle between buds reprices
only **4.4% / 8.6%** of raws (verbose `[dada] p-update churn`), which *looked*
like green light for a p-ordered heap — but cost-modeling killed it: churn ×
O(log nraw) *scattered* heap ops approaches (rev: exceeds) the cache-friendly
*sequential* scan it would replace. Same bandwidth/cache lesson, a third time.

The design that survives folds the candidate tracking into the pass that already
touches exactly those churned raws. `b_p_update`, while repricing each dirty
cluster's members, now also caches that cluster's best abundance/prior candidate
(`Bi::bud_min`, using `b_bud`'s exact filters and tie-break). `b_bud_incremental`
then combines the per-cluster minima in **O(nclusters)** instead of O(nraw),
seeded identically with cluster 0's center. Cache maintenance rides the existing
p-update memory traffic at ~zero extra cost — no heap, no log-n, no scatter.

Correctness rests on a verified invariant: between consecutive buds a candidate's
p or eligibility changes *iff* its cluster was flagged `update_e` (`b_compare`
writes `raw.comp` only for the initial cluster/centers; every move and reads
change goes through `bi_add_raw`/`bi_pop_raw`, which flag it). A debug-only
cross-check asserts the combine equals a full serial scan on *every* bud — it
passed across the entire integration suite in debug.

Result on the 384-sample pool (vs the shuffle-win parent, one run each, distinct
commits `f1f6c71`→`44fcca6`):

| Metric | parent (serial bud) | **incremental bud** | Δ |
|---|---:|---:|---:|
| `bud` phase (fwd / rev) | ~20s / ~20s | **0.07s / 0.05s** | ~eliminated |
| dada_fwd wall | 182.5s | **166.4s** | **−8.8%** |
| dada_rev wall | 148.4s | **131.2s** | **−11.6%** |
| dada_fwd cores (of 24) | 13.6 | **15.0** | +1.5 |
| dada_rev cores (of 24) | 10.2 | **11.5** | +1.3 |
| **total pipeline wall** | 391.4s | **359.0s** | **−8.3%** |
| pipeline peak RSS | 1623 MB | 1614 MB | ≈0 |

Byte-identical: all 725 sequences, 362 samples, and 197,242 pre- / 75,849
post-chimera nonzero counts match exactly (order-independent). No RSS cost (the
cache is ~100 KB). The win **stacks** on the shuffle win — cumulatively `main`→now
is ~427→359s (**≈−16%**). With `bud` at ~0, the shuffle phase is now the dominant
remaining serial fraction and the next ceiling on core utilization.

## Regenerating these tables

After a run, distill its `summary.csv` to Markdown:

```bash
# scorecard across modes (one row per run)
python3 dev/benchmark/bench_table.py \
    "pooled=bench_true/summary.csv" \
    "per-sample=bench_false/summary.csv" \
    "pseudo=bench_pseudo/summary.csv"

# per-step breakdown for one run
python3 dev/benchmark/bench_table.py --per-step bench_pseudo/summary.csv
```

Paste the output into the tables above.
