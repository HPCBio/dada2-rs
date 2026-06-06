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

## Illumina MiSeq (SOP)

!!! note "Populate from a cluster run"
    Generate with:
    ```bash
    python3 comparison/benchmark/bench_table.py \
        "pooled=bench_true/summary.csv" \
        "per-sample=bench_false/summary.csv" \
        "pseudo=bench_pseudo/summary.csv"
    ```
    then paste the table below.

Run using 24 threads using release-native `dada2-rs`, `v0.1.1-a20fee47`.

*Overall workflow (filter and trim FASTQ -> removing chimeras)*

| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (R÷rs) |
|---|---:|---:|---:|---:|---:|---:|
| MiSeqSOP_pooled_k5 | 544.5 | 1361.4 | 2.5× | 4535 | 6082 | 1.3× |
| MiSeqSOP_pooled_k6 | 488.9 | 1361.4* | 2.8x | 4529 | 6082* | 1.3x |
| MiSeqSOP_pooled_k7 | 487.0 | 1361.4* | 2.8x | 4521 | 6082* | 1.3x |
| MiSeqSOP_pseudo_k5 | 245.2 | 1158.1 | 4.7× | 2208 | 1841 | 0.8× |
| MiSeqSOP_pseudo_k6 | 246.5 | 1158.1* | 4.7x | 2229 | 1841* | 0.8x |
| MiSeqSOP_pseudo_k7 | 244.8 | 1158.1* | 4.7x | 2242 | 1841* | 0.8x |
| MiSeqSOP_nopool_k5 | 132.3 | 778.7 | 5.9× | 1457 | 1659 | 1.1× |
| MiSeqSOP_nopool_k6 | 151.1 | 778.7* | 5.2x | 1542 | 1659* | 1.1x |
| MiSeqSOP_nopool_k7 | 151.7 | 778.7* | 5.2x | 1469 | 1659* | 1.1x |

`*` - values from k=5 run

Moderate k-mer difference in pooled runs only; we only recommend k=5 or for pooled runs k=6.

### Per step comparisons

### Per step comparison at k=5, full pooling:

**MiSeqSOP pooled, k=5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 13.7 | 20.9 | 23 | 40.9 | 3.0× |
| learn_fwd | 38.3 | 22.3 | 718 | 89.6 | 2.3× |
| learn_rev | 38.8 | 22.0 | 852 | 111.9 | 2.9× |
| dada_fwd | 251.2 | 15.8 | 4535 | 467.7 | 1.9× |
| dada_rev | 187.9 | 12.3 | 3622 | 368.3 | 2.0× |
| merge | 8.6 | 18.7 | 1569 | 256.3 | 29.7× |
| make_table | 0.5 | 1.0 | 329 | 0.3 | 0.6× |
| remove_bimera | 5.5 | 21.9 | 47 | 4.2 | 0.8× |
| TOTAL | 544.5 | 15.7 | 4535 | 1361.4 | 2.5× |

**MiSeqSOP, pseudo, k5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 13.2 | 21.0 | 21 | 36.4 | 2.8× |
| learn_fwd | 32.8 | 22.4 | 721 | 80.7 | 2.5× |
| learn_rev | 33.4 | 21.9 | 848 | 100.0 | 3.0× |
| dada_fwd | 96.7 | 20.1 | 2208 | 413.1 | 4.3× |
| dada_rev | 63.5 | 18.9 | 1551 | 303.4 | 4.8× |
| merge | 4.9 | 17.3 | 1583 | 203.9 | 41.6× |
| make_table | 0.2 | 1.0 | 139 | 0.1 | 0.7× |
| remove_bimera | 0.4 | 17.3 | 21 | 0.3 | 0.9× |
| TOTAL | 245.2 | 20.3 | 2208 | 1158.1 | 4.7× |

**MiSeqSOP, no pooling, k5** — dada2-rs vs R (R-single end-to-end wall)

| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |
|---|---:|---:|---:|---:|---:|
| filter | 13.3 | 20.9 | 21 | 35.5 | 2.7× |
| learn_fwd | 33.1 | 22.3 | 724 | 80.7 | 2.4× |
| learn_rev | 33.1 | 22.0 | 852 | 108.8 | 3.3× |
| dada_fwd | 30.2 | 20.9 | 594 | 197.0 | 6.5× |
| dada_rev | 18.8 | 20.2 | 470 | 142.4 | 7.6× |
| merge | 3.4 | 18.2 | 1457 | 193.7 | 56.6× |
| make_table | 0.1 | 1.0 | 60 | 0.1 | 1.2× |
| remove_bimera | 0.3 | 11.6 | 21 | 0.2 | 0.6× |
| TOTAL | 132.3 | 21.3 | 1457 | 778.7 | 5.9× |

## PacBio HiFi

PacBio is benchmarked at both `--kmer-size 5` (matched to R's fixed
`KMER_SIZE = 5` — isolates kernel + threading speed) and `--kmer-size 7`
(dada2-rs default — adds the k-mer screen-effectiveness gain). See the
[PacBio notes](benchmarking.md#7-pacbio-vs-illumina-specifics).

| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (R÷rs) |
|---|---:|---:|---:|---:|---:|---:|
| pooled (k=5) | | | | | | |
| pseudo (k=5) | | | | | | |
| pseudo (k=7) | | | | | | |

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
