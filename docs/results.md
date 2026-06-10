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
- **Peak RSS** is the per-process high-water mark; the ratio is R ÷ dada2-rs
  (>1× means dada2-rs uses less memory, <1× means more).
- **Build:** `release-native` unless noted. See
  [build target](benchmarking.md#6-build-target-matters).
- **Correctness** (ASV concordance) is validated separately — see
  [concordance tooling](benchmarking.md#5-concordance-validation-tooling).

!!! note "Populate from a cluster run"
    We generate the results using `comparison/benchmark/bench_pooled.py` 
    and summarize with:
    ```bash
    python3 comparison/benchmark/bench_table.py \
        "pooled=bench_true/summary.csv" \
        "per-sample=bench_false/summary.csv" \
        "pseudo=bench_pseudo/summary.csv"
    ```
    then paste the table below.

## PacBio HiFi

* `dada2-rs v0.1.1-9a0c5da3` (v0.1.1 prerelease)
* Data is from Hergenrother 2024, 93 total samples, PacBio Sequel IIe.
* Default settings using FASTQ input for `dada` in both tools

### Overall workflow

PacBio is benchmarked from `--kmer-size 5` (matched to R's fixed
`KMER_SIZE = 5` — isolates kernel + threading speed) to `--kmer-size 7`
(dada2-rs recommended default for PacBio — adds the k-mer screen-effectiveness gain). See the [PacBio notes](benchmarking.md#7-pacbio-vs-illumina-specifics).

** Results in progress, more to be added **

| Run | Pooling mode | k-mer | sample jobs | dada2-rs wall (s) | R wall (s) | Speedup (rs vs R) | dada2-rs peak (MB) | R peak (MB) | Peak RSS (rs vs R) |
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
against cached. Numbers are filtered per stack with
`compare_bench.py --stack <stack>`; only the `dada` step moves materially (other
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

### Pre/post-update A/B (memory + speed work)

A dada2-rs-vs-itself regression check (PacBio, node-exclusive, no cache flag).

**`pre`** is `main` before recent memory/speed improvements, **`post`** is current `main`. This
spans the whole update window — a *composite* of the memory-representation change
(#23, `qual_sum:[u32]` deferred division) and the `learn-errors`/`dada` speedups (issue
#3 checklist: skip per-iteration setup, parallel `build_trans_mat`, SIMD `kmer_dist8`),
not an isolated change. **R is not a variable here, so no R numbers apply.** 

All rows are the same PacBio dataset, varying mode and k. `learn` wall falls ~−13 to −15% in every row below; the `dada` step is broken out:

| mode | k | dada wall Δ | dada CPU Δ | dada peak (pre → post) | peak Δ |
|---|---|---:|---:|---:|---:|
| pooled (`--pool true`) | 5 | **−14.9%** | −15.5% | 36.14 → 28.61 GB | **−20.8%** |
| pooled (`--pool true`) | 7 | **−14.4%** | −15.2% | 36.04 → 28.69 GB | **−20.4%** |
| pseudo (streaming default) | 5 | **−14.9%** | −15.1% | 6.00 → 6.30 GB | +4.9% |
| pseudo (cached, `--cache-samples`) | 7 | **−16.8%** | −17.4% | 18.27 → 13.12 GB | **−28.2%** |

**Verdict:** the **speed win is mode-independent**
— `dada` wall and CPU fall together (~−15 to −17%) at flat cores (−0.3 to −1.0%) in every mode, i.e. genuinely less work, not better parallel packing. 
- The **memory win is resident-bound**: the modes that hold all samples resident shed a large fraction of the high-water mark (pooled −20%, ~7.5 GB; cached −28%, ~5.2 GB), while the streaming pseudo path keeps only a small working set, so there is no RSS reduction to capture (a small ~300 MB uptick, single rep — read as flat, not a regression). 
- The peak ordering tracks residency exactly — pooled 36 GB > cached 18 GB > streaming 6 GB. Pooled `dada` peak is ~36 GB at both k values: at that scale the resident all-samples data dominates RSS, the 4^k screen array is washed out, and the memory delta is representation-driven, not k-related.

## Illumina MiSeq (F1000, 384 samples)

The majority of recent memory work focused on PacBio data. How does this affect Illumina data processing? 

* `dada2-rs v0.1.1-9a0c5da3` (v0.1.1 prerelease)
* Data is the F1000 'MiSeqSOP' full data set: 384 samples, MiSeq v2 run.
* Default settings using FASTQ input for `dada` in both tools

* 24 threads using release-native `dada2-rs`, `v0.1.1-a20fee47`.

*Overall workflow (filter and trim FASTQ -> removing chimeras)*

| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (R÷rs) |
|---|---:|---:|---:|---:|---:|---:|
| MiSeqSOP_pooled_k5 | 544.5 | 1361.4 | **2.5×** | 4535 | 6082 | 1.3× |
| MiSeqSOP_pooled_k6 | 488.9 | 1361.4* | **2.8x** | 4529 | 6082* | 1.3x |
| MiSeqSOP_pooled_k7 | 487.0 | 1361.4* | **2.8x** | 4521 | 6082* | 1.3x |
| MiSeqSOP_pseudo_k5 | 245.2 | 1158.1 | **4.7×** | 2208 | 1841 | 0.8× |
| MiSeqSOP_pseudo_k6 | 246.5 | 1158.1* | **4.7x** | 2229 | 1841* | 0.8x |
| MiSeqSOP_pseudo_k7 | 244.8 | 1158.1* | **4.7x** | 2242 | 1841* | 0.8x |
| MiSeqSOP_nopool_k5 | 132.3 | 778.7 | **5.9×** | 1457 | 1659 | 1.1× |
| MiSeqSOP_nopool_k6 | 151.1 | 778.7* | **5.2x** | 1542 | 1659* | 1.1x |
| MiSeqSOP_nopool_k7 | 151.7 | 778.7* | **5.2x** | 1469 | 1659* | 1.1x |

`*` - values from k=5 run

Moderate k-mer difference in pooled runs only, but an overall improvement in performance (2.5-5.9x). We only recommend k=5 for most runs, maybe k=6 for pooled runs.

### Per step comparisons

**MiSeqSOP pooled, k=5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 13.7 | 20.9 | 23 | 40.9 | **3.0×** |
| learn_fwd | 38.3 | 22.3 | 718 | 89.6 | **2.3×** |
| learn_rev | 38.8 | 22.0 | 852 | 111.9 | **2.9×** |
| dada_fwd | 251.2 | 15.8 | 4535 | 467.7 | **1.9×** |
| dada_rev | 187.9 | 12.3 | 3622 | 368.3 | **2.0×** |
| merge | 8.6 | 18.7 | 1569 | 256.3 | **29.7×** |
| make_table | 0.5 | 1.0 | 329 | 0.3 | 0.6× |
| remove_bimera | 5.5 | 21.9 | 47 | 4.2 | 0.8× |
| TOTAL | 544.5 | 15.7 | 4535 | 1361.4 | **2.5×** |

**MiSeqSOP, pseudo, k5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 13.2 | 21.0 | 21 | 36.4 | **2.8×** |
| learn_fwd | 32.8 | 22.4 | 721 | 80.7 | **2.5×** |
| learn_rev | 33.4 | 21.9 | 848 | 100.0 | **3.0×** |
| dada_fwd | 96.7 | 20.1 | 2208 | 413.1 | **4.3×** |
| dada_rev | 63.5 | 18.9 | 1551 | 303.4 | **4.8×** |
| merge | 4.9 | 17.3 | 1583 | 203.9 | **41.6×** |
| make_table | 0.2 | 1.0 | 139 | 0.1 | 0.7× |
| remove_bimera | 0.4 | 17.3 | 21 | 0.3 | 0.9× |
| TOTAL | 245.2 | 20.3 | 2208 | 1158.1 | **4.7×** |

**MiSeqSOP, no pooling, k5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 13.3 | 20.9 | 21 | 35.5 | **2.7×** |
| learn_fwd | 33.1 | 22.3 | 724 | 80.7 | **2.4×** |
| learn_rev | 33.1 | 22.0 | 852 | 108.8 | **3.3×** |
| dada_fwd | 30.2 | 20.9 | 594 | 197.0 | **6.5×** |
| dada_rev | 18.8 | 20.2 | 470 | 142.4 | **7.6×** |
| merge | 3.4 | 18.2 | 1457 | 193.7 | **56.6×** |
| make_table | 0.1 | 1.0 | 60 | 0.1 | 1.2× |
| remove_bimera | 0.3 | 11.6 | 21 | 0.2 | 0.6× |
| TOTAL | 132.3 | 21.3 | 1457 | 778.7 | **5.9×** |

Notably the big improvement is with `merge-pairs`, largely due to threading, though overall steps are all faster.

## Regenerating these tables

After a run, distill its `summary.csv` to Markdown:

```bash
# scorecard across modes (one row per run)
python3 comparison/benchmark/bench_table.py \
    "pooled=bench_true/summary.csv" \
    "per-sample=bench_false/summary.csv" \
    "pseudo=bench_pseudo/summary.csv"

# per-step breakdown for one run
python3 comparison/benchmark/bench_table.py --per-step bench_pseudo/summary.csv
```

Paste the output into the tables above.
