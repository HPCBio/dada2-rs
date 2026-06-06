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

| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (R÷rs) |
|---|---:|---:|---:|---:|---:|---:|
| MiSeqSOP_pooled_k5 | 544.5 | 1361.4 | 2.5× | 4535 | 6082 | 1.3× |
| MiSeqSOP_pooled_k6 | 488.9 | — | — | 4529 | — | — |
| MiSeqSOP_pooled_k7 | 487.0 | — | — | 4521 | — | — |
| MiSeqSOP_pseudo_k5 | 245.2 | 1158.1 | 4.7× | 2208 | 1841 | 0.8× |
| MiSeqSOP_pseudo_k6 | 246.5 | — | — | 2229 | — | — |
| MiSeqSOP_pseudo_k7 | 244.8 | — | — | 2242 | — | — |
| MiSeqSOP_nopool_k5 | 132.3 | 778.7 | 5.9× | 1457 | 1659 | 1.1× |
| MiSeqSOP_nopool_k6 | 151.1 | — | — | 1542 | — | — |
| MiSeqSOP_nopool_k7 | 151.7 | — | — | 1469 | — | — |

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

After a cluster run, distill its `summary.csv` to Markdown:

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
