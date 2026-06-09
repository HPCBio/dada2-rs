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

Data is from Hergenrother 2024, 93 total samples, PacBio Sequel IIe.

PacBio is benchmarked from `--kmer-size 5` (matched to R's fixed
`KMER_SIZE = 5` — isolates kernel + threading speed) to `--kmer-size 7`
(dada2-rs recommended default for PacBio — adds the k-mer screen-effectiveness gain). See the [PacBio notes](benchmarking.md#7-pacbio-vs-illumina-specifics).

| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (R÷rs) |
|---|---:|---:|---:|---:|---:|---:|
| PacBio_pooled_k5 | 6549.1 | 15162.8 | **2.3×** | 36931 | 40024 | **1.1×** |
| PacBio_pooled_k6 | 1442.8 | 15162.8* | **10.5x** | 42067 | 40024* | **1.0x** |
| PacBio_pooled_k7 | 1065.0 | 15162.8* | **14.2x** | 61082 | 40024* | **0.7x** |
| PacBio_pseudo_k5 | 3466.6 | 10880.4 | **3.1×** | 13842 | 3069 | **0.2×** |
| PacBio_pseudo_k6 | 784.5 | 10880.4* | **13.9x** | 15110 | 3069* | **0.2x** |
| PacBio_pseudo_k7 | 609.8 | 10880.4* | **17.8x** | 19219 | 3069* | **0.2x** |
| PacBio_nopool_k5 | 1535.8 | 8129.6 | **5.3×** | 5838 | 2808 | **0.5×** |
| PacBio_nopool_k6 | 378.1 | 8129.6* | **21.5x** | 6356 | 2808* | **0.4x** |
| PacBio_nopool_k7 | 283.7 | 8129.6* | **28.7x** | 8873 | 2808* | **0.3x** |

`*` - values from k=5 run

Note that R DADA2 currently has a fixed (hard-coded) k-mer size of 5 for screening. This seems to have the effect of reducing k-mer screening efficiency almost to a no-op, with almost every unique sequence proceeding to NW alignment.

Our current implementation currently uses more memory (last column), an issue we are chasing down. Two of the pooling modes have a low memory setting which we are also benchmarking which should alleviate this with some performance tradeoff.

### Preliminary observation and results

**Recommendation:** use `--kmer-size 7` for PacBio runs, especially if your compute has a decent memory footprint. 

- Increasing the k-mer length marginally has a significant effect, reducing walltime dramatically. 
- Preliminary followup (to be added) indicates this has essentially **no impact** on the ASVs derived post-chimera removal as well as their counts (to be added below), and may actually increase sensitivity slightly, with a few more ASVs recovered in the `dada2-rs` runs
- Pre-chimera results differ slightly from R DADA2 (with some ASVs found unique to both runs). However as noted above these are essentially removed with chimera removal, suggesting their lower abundance
- `dada2-rs` utilizes more memory than expected for these steps, which we are in the process of evaluating

We are also evaluating the general effect that increasing k-mer size has on clustering with the different pooling results, including basic exploratory testing using higher k-mers. Initial k=8 results confirms this increases memory usage dramatically (>100GB), so we don't anticipate this will have practical use.

### Per step comparisons

**PacBio pooled, k=7** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| remove_primers | 23.0 | 20.8 | 31 | 5356.4 | **232.6×** |
| learn | 29.4 | 20.4 | 2850 | 205.9 | **7.0×** |
| dada | 984.4 | 20.1 | 61082 | 9494.9 | **9.6×** |
| make_table | 0.3 | 1.0 | 111 | 0.1 | 0.3× |
| remove_bimera | 27.8 | 23.7 | 73 | 30.5 | 1.1× |
| TOTAL | 1065.0 | 20.2 | 61082 | 15162.8 | **14.2×** |

**PacBio pseudo, k=7** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| remove_primers | 23.5 | 20.0 | 29 | 5434.3 | **230.9×** |
| learn | 29.7 | 20.2 | 2873 | 215.2 | **7.2×** |
| dada | 542.2 | 21.3 | 19219 | 5135.7 | **9.5×** |
| make_table | 0.2 | 1.0 | 70 | 0.1 | 0.3× |
| remove_bimera | 14.0 | 22.6 | 86 | 12.7 | 0.9× |
| TOTAL | 609.8 | 21.2 | 19219 | 10880.4 | **17.8×** |

**PacBio, no pooling, k=7** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| remove_primers | 27.6 | 20.4 | 29 | 5373.3 | **195.0×** |
| learn | 33.9 | 20.7 | 2902 | 239.4 | **7.1×** |
| dada | 214.4 | 21.6 | 8873 | 2426.9 | **11.3×** |
| make_table | 0.1 | 1.0 | 53 | 0.0 | 0.3× |
| remove_bimera | 7.7 | 22.0 | 70 | 8.1 | 1.1× |
| TOTAL | 283.7 | 21.4 | 8873 | 8129.6 | **28.7×** |

Note that most of the walltime improvement is actually dominated by the `remove-primers` step (which combines two steps from DADA2: `removePrimers` and filterAndTrim). However significant improvements are also apparent for both the `learn-errors` and `dada` steps for each pooling mode.

### Streaming vs cached (`--cache-samples`)

`dada-pseudo` streams samples by default (re-reading each sample per round) rather
than holding every sample's dereplicated state resident. `--cache-samples` opts
into the old all-in-memory behavior. This characterizes that tradeoff on PacBio
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
Two cross-stack notes: streaming-vs-streaming, dada2-rs is both leaner and faster on
`dada` (2.34 GB / 2092 s vs 2.86 GB / 3738 s); cached-vs-cached, dada2-rs is ~1.8×
leaner (8.18 GB vs 14.86 GB) — R's resident-derep cache is not magically compact.
See [issue #22](https://github.com/HPCBio/dada2-rs/issues/22) for the keep-but-default-off
decision.

### Pre/post-update A/B (memory + speed work)

A dada2-rs-vs-itself regression check (PacBio, node-exclusive, no cache flag):
**`pre`** is `main` before the memory/speed work, **`post`** is current `main`. This
spans the whole update window — it is a *composite* of the memory-representation change
(#23, `qual_sum:[u32]` deferred division) and the `learn-errors`/`dada` speedups (issue
#3 checklist: skip per-iteration setup, parallel `build_trans_mat`, SIMD `kmer_dist8`),
not an isolated change. R is not a variable here, so no R numbers apply. All rows are the
same PacBio dataset, varying mode and k. `learn` wall falls ~−13 to −15% in every row
below; the `dada` step is broken out:

| mode | k | dada wall Δ | dada CPU Δ | dada peak (pre → post) | peak Δ |
|---|---|---:|---:|---:|---:|
| pooled (`--pool true`) | 5 | **−14.9%** | −15.5% | 36.14 → 28.61 GB | **−20.8%** |
| pooled (`--pool true`) | 7 | **−14.4%** | −15.2% | 36.04 → 28.69 GB | **−20.4%** |
| pseudo (streaming default) | 5 | **−14.9%** | −15.1% | 6.00 → 6.30 GB | +4.9% |
| pseudo (cached, `--cache-samples`) | 7 | **−16.8%** | −17.4% | 18.27 → 13.12 GB | **−28.2%** |

**Verdict:** the **speed win is mode-independent** — `dada` wall and CPU fall together
(~−15 to −17%) at flat cores (−0.3 to −1.0%) in every mode, i.e. genuinely less work, not
better parallel packing. The **memory win is resident-bound**: the modes that hold all
samples resident shed a large fraction of the high-water mark (pooled −20%, ~7.5 GB;
cached −28%, ~5.2 GB), while the streaming pseudo path keeps only a small working set, so
there is no RSS reduction to capture (a small ~300 MB uptick, single rep — read as flat,
not a regression). The peak ordering tracks residency exactly — pooled 36 GB > cached
18 GB > streaming 6 GB. Pooled `dada` peak is ~36 GB at both k values: at that scale the
resident all-samples data dominates RSS, the 4^k screen array is washed out, and the
memory delta is representation-driven, not k-related.

## Illumina MiSeq (F1000, 384 samples)

Run using 24 threads using release-native `dada2-rs`, `v0.1.1-a20fee47`.

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
