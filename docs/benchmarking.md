# Performance & Benchmarking

This document is the reference for the tooling and log output used to evaluate
`dada2-rs` performance and to benchmark it head-to-head against R DADA2. It is
written in plain Markdown so it can be dropped into MkDocs or Sphinx (via
`myst-parser`) for a ReadTheDocs site with no changes.

- **Where things live:** `comparison/benchmark/` (the head-to-head harness),
  `comparison/` (validation/error-model scripts), `scripts/` (plotting and
  read-tracking helpers).
- **What you get:** per-step wall time, CPU time, **effective cores**, and peak
  RSS for both stacks; scaling sweeps; and ASV/error-matrix concordance checks.

---

## 1. The benchmark harness (`comparison/benchmark/`)

A single Python driver runs the **full pipeline** for each stack, times every
step, captures per-process peak memory, and prints a symmetric side-by-side
table plus a CSV.

| File | Role |
|---|---|
| `bench_pooled.py` | Orchestrator (Python). Runs dada2-rs and/or R, captures timing/RSS/CPU, writes `summary.csv`. |
| `run_dada2_pooled.R` | R reference pipeline as **one process** (fair end-to-end wall + one overall RSS). |
| `bench_step.R` | R reference pipeline as **one process per step** (gives per-step RSS). |

### 1.1 Function-vs-function design

Each `--pool` mode invokes the dada2-rs subcommand that is the direct analog of
the R function — one process each, so the comparison is *function vs function*,
not *harness vs harness*:

| `--pool` | dada2-rs subcommand | R analog |
|---|---|---|
| `true` | `dada-pooled` | `dada(..., pool = TRUE)` |
| `false` | `dada` (multi-input, `--sample-jobs`) | `dada(..., pool = FALSE, multithread = N)` |
| `pseudo` | `dada-pseudo` | `dada(..., pool = "pseudo")` |

Pipeline stages (each a timed row in the output):

- **Illumina (paired):** `filter` → `learn_fwd` / `learn_rev` → `dada_fwd` /
  `dada_rev` → `merge` → `make_table` → `remove_bimera`
- **PacBio (single-end):** `remove_primers` (trim + orient + filter in one) →
  `learn` → `dada` → `make_table` → `remove_bimera`

### 1.2 How metrics are captured

All timing/memory comes from `os.wait4()` on each child process — no
`/usr/bin/time`, no `/proc` polling, so it is portable across the cluster and
macOS:

| Metric | Source | Meaning |
|---|---|---|
| `wall_s` | wall clock around the process | elapsed time |
| `cpu_s` | `ru_utime + ru_stime` | total CPU-seconds consumed |
| `cores` | `cpu_s / wall_s` | **effective cores used** (≈ threads = fully utilized) |
| `peak_rss` | `ru_maxrss` (kB on Linux, bytes on macOS → normalized) | peak resident set high-water mark |

For **fanned** steps (`filter`, `remove_primers`, which run one subprocess per
sample concurrently), the row aggregates: `wall_s` = wall clock of the whole
batch, `cpu_s` = sum across children, `peak_rss` = max of any single child.

> **`cores` is the key diagnostic.** It tells you whether a step actually uses
> the threads you gave it. `cores ≈ threads` → saturated; `cores ≪ threads` →
> under-utilized. See §3 for how to read it together with `cpu_s`.

### 1.3 Key options

| Option | Purpose |
|---|---|
| `platform` (positional) | `illumina` or `pacbio` |
| `input` (positional) | directory of raw FASTQ files |
| `--dada2rs PATH` | **required** — path to the binary (no auto-discovery; build target matters, see §6) |
| `--threads N` | total threads |
| `--pool {true,false,pseudo}` | denoising mode (see §1.1) |
| `--sample-jobs N` | samples denoised concurrently for `pseudo`/`false` (default `round(threads/4)`) |
| `--low-memory` | `pseudo` only: stream samples instead of caching (bounds peak RSS) |
| `--run-r` | also run the R pipeline (off by default → dada2-rs only) |
| `--r-mode {split,single,both}` | how to run R (see §1.5) |
| `--no-run-rust` | skip the dada2-rs pipeline |
| `--nbases N` | bases to subsample for error learning |
| `--thread-sweep N,N,...` | thread-scaling study (§2.1) |
| `--sample-jobs-sweep N,N,...` | samples-in-flight study (§2.2) |

Platform-specific filter/denoise knobs: `--trunc-len`, `--max-ee`, `--trunc-q`,
`--max-n` (Illumina); `--min-len`, `--max-len`, `--band`, `--homo-gap`,
`--kmer-size`, `--primer-fwd`, `--primer-rev`, `--max-mismatch` (PacBio). Run
`python3 comparison/benchmark/bench_pooled.py --help` for the full list.

### 1.4 Output

The console prints a per-step table for each stack and a `HEADLINE` block
(end-to-end speedup, pooled/per-sample denoise speedup, peak-RSS ratio). The
machine-readable form is `<outdir>/summary.csv`:

```csv
stack,step,wall_s,cpu_s,cores,maxrss_kb
dada2-rs,dada_fwd,80.53,1562.65,19.40,2219336
...
R-split,dada_fwd,447.67,2164.32,4.84,1470632
R-single,TOTAL,1208.38,5417.63,4.48,2021736
```

`maxrss_kb` is the exact value; the printed table rounds RSS to MiB. Sweeps
write `sweep.csv` / `sweep_jobs.csv` (§2).

### 1.5 R run modes (`--r-mode`)

R is captured two ways because they trade off differently:

- **`single`** — the whole R pipeline as one process. Gives a **fair end-to-end
  wall** and one overall peak RSS. Per-step wall comes from R's `system.time()`;
  per-step RSS/cores are not available (one accumulating process).
- **`split`** — one `Rscript` per step (state passed via `.rds`). Gives
  **per-step RSS**, but each step reloads R + dada2 (~150–200 MB, several
  seconds), so per-step *wall* is inflated.
- **`both`** (default) reports both; the `HEADLINE` uses `single` for the fair
  end-to-end wall and `split` only where per-step RSS is needed.

> **Caveat — not perfectly apples-to-apples on RSS:** dada2-rs's per-step RSS is
> per separate process (memory released between steps); R-single's is one
> accumulating process. The numbers are directionally comparable, not identical
> in methodology.

### 1.6 Examples

```bash
# dada2-rs only, Illumina, pooled
python3 comparison/benchmark/bench_pooled.py illumina /data/MiSeqSOP \
    --dada2rs target/release-native/dada2-rs --threads 24 --pool true

# head-to-head vs R, pseudo-pooling
python3 comparison/benchmark/bench_pooled.py illumina /data/MiSeqSOP \
    --dada2rs target/release-native/dada2-rs --threads 24 --pool pseudo --run-r

# PacBio HiFi (raw, primered reads), supply primers
python3 comparison/benchmark/bench_pooled.py pacbio /data/HiFi \
    --dada2rs target/release-native/dada2-rs --threads 24 --pool pseudo \
    --primer-fwd AGRGTTYGATYMTGGCTCAG --primer-rev RGYTACCTTGTTACGACTT --run-r

# pseudo in streaming (low-memory) mode
python3 comparison/benchmark/bench_pooled.py illumina /data/MiSeqSOP \
    --dada2rs target/release-native/dada2-rs --threads 24 --pool pseudo --low-memory
```

### 1.7 Distilling results to Markdown

`bench_table.py` turns one or more `summary.csv` files into a Markdown table —
useful since the large-dataset runs happen on the cluster and the tables are
pasted into the docs (see [Benchmark results](results.md)).

```bash
# Scorecard — one row per run, end-to-end head-to-head:
python3 comparison/benchmark/bench_table.py \
    "pooled=bench_true/summary.csv" \
    "per-sample=bench_false/summary.csv" \
    "pseudo=bench_pseudo/summary.csv"

# Per-step breakdown for a single run:
python3 comparison/benchmark/bench_table.py --per-step bench_pseudo/summary.csv
```

`LABEL=path` names each run; R columns compare against the `R-single` rows (fair
end-to-end wall + overall RSS), and rust-only runs simply omit them.

---

## 2. Scaling studies (dada2-rs only; skip R)

Both sweeps **prepare inputs once** (filter + learn) and then re-run **only the
denoise step** at each setting — isolating the scaling behavior of the inference.

### 2.1 Thread sweep — `--thread-sweep 1,2,4,8,16,24`

Varies `--threads`, pinning `--sample-jobs 1` so it measures **single-sample
thread scaling**. Output (`sweep.csv`): `threads, wall_s, cpu_s, cores, speedup,
efficiency`. `efficiency = speedup / (threads / baseline_threads)`.

Reads the Amdahl story: efficiency rolling off and `cpu_s` *climbing* with
threads = the per-sample alignment map is too small to feed the threads, so
rayon workers **spin-wait** (burning CPU without cutting wall).

### 2.2 Sample-jobs sweep — `--sample-jobs-sweep 1,2,3,4,6,8`

At **fixed `--threads`**, varies `--sample-jobs` (samples denoised concurrently,
each on a `threads / jobs` sub-pool). Applies to `--pool pseudo` or `false`.
Output (`sweep_jobs.csv`): `sample_jobs, threads_per_job, wall_s, cpu_s, cores,
maxrss_kb, speedup`, with the fastest row marked.

How to read it:

- **Minimum `wall_s`** is the optimum number of samples-in-flight.
- **`cpu_s` rising again** at high `jobs` = sub-pools too small, per-map overhead
  returns (the mirror image of the spin in §2.1).
- **`peak_rss` rising with `jobs`** = more concurrent working sets — the
  speed/memory trade for picking `jobs`.

The default `--sample-jobs = round(threads/4)` came from this sweep (the
wall-time curve plateaus around ~4 threads/sample).

---

## 3. Interpreting the metrics

**`cores` (CPU/wall) finds the problem; `cpu_s` confirms the fix.** A step that
is busy but not faster shows up as high `cores` *and* high `cpu_s` *and* poor
speedup — that signature is **spin-waiting**, not useful work. When you fix it
(e.g. right-sizing sub-pools), the wall drops **and `cpu_s` drops** (reclaimed
spin) — so `cores` may look unchanged because both its numerator and denominator
shrink. **Trust `wall_s` and `cpu_s`, not `cores` alone, when judging a fix.**

Rules of thumb:

| Observation | Interpretation |
|---|---|
| `cores ≈ threads` | step is saturating the machine |
| `cores ≪ threads`, low `cpu_s` | genuinely idle (load imbalance, too few samples, serial backbone) |
| `cores` high, `cpu_s` high, poor speedup | spin-waiting on an under-fed parallel region |
| `cpu_s` roughly equal between stacks, wall differs | same work, different parallel efficiency |

**Peak RSS caveats:** it is the kernel high-water mark of a *single* process. For
fanned/concurrent steps the reported value is the max of one worker, so it
*understates* the true concurrent system memory when many samples run at once.

---

## 4. Built-in instrumentation (the binary's own logs)

Run any `dada` / `dada-pooled` / `dada-pseudo` with `--verbose` for a per-run
breakdown on stderr (also captured in the harness `*.log` files):

```
ALIGN: 123456 aligns, 654321 shrouded (1000 raw).
[dada] phase times (serial except compare-map): compare=12.3s (map=11.0s parallel, store=1.3s serial)  shuffle=0.4s  bud=0.2s  p_update=0.3s
[dada] map parallel efficiency: 74% (busy=320s / map=11.0s × 24 threads)
```

- **`ALIGN`** — alignments performed vs *shrouded* (skipped by the k-mer
  pre-screen). On long reads at small `k` the shroud count is ~0 (screen is a
  no-op); larger `k` engages it.
- **`phase times`** — only `compare-map` is parallel. `shuffle`, `bud`,
  `p_update`, and the `store` half of compare are **serial**. If those dominate,
  you are Amdahl-bound on the greedy cluster-building backbone, not on the map.
- **`map parallel efficiency`** — how well the parallel region itself uses
  threads; well below 100% means in-region load imbalance.

**Error-model parameter warning.** `dada*` records the alignment params in the
error-model JSON and warns when the denoise params disagree:

```
[dada] warning: 1 dada parameter(s) differ from error model errors_pacbio.json;
pass --inherit-err-params to adopt the err model's values:
  homo_gap_p = -1 (err model: -8)
```

This means the error model was *learned* with different params than you are
*denoising* with. Either pass the same params to `learn-errors` (preferred — so
the model matches), or `--inherit-err-params` to adopt the model's values.
Common on PacBio if `--homo-gap-p` is set for `dada` but not for `learn-errors`.

---

## 5. Concordance & validation tooling

Performance only matters if the results are right. These check correctness
against R, independent of timing:

| Tool | What it checks |
|---|---|
| `scripts/compare_asvs.py` | ASV table A/B diff (sequences present, per-sample abundances) between two runs |
| `comparison/run_rust_errors.sh` | runs filter → derep → `learn-errors` on one sample with R-matched params |
| `comparison/compare_errors.R` | compares the learned error matrices against R's `learnErrors()` |
| `comparison/run_kmer_sweep.sh` | sweeps `--kmer-size`, reporting shroud %, ASV count, wall, RSS (issue #15) |
| `comparison/plot_cluster_diag.R` | plots per-iteration cluster diagnostics from `--diag-dir` output |
| `scripts/plot_errors.R` | plots an error-model JSON |
| `scripts/track_reads.py` | per-stage read-count tracking (the DADA2 "track" table) |
| `scripts/summarize_learn_errors.py` | summarizes a learned error model |

Typical correctness loop for a benchmark run:

1. Run the harness for the stack(s).
2. Export both final sequence tables to TSV (`seq-table-to-tsv`) or compare the
   per-sample `dada` JSONs with `compare_asvs.py`.
3. Confirm ASV concordance before trusting any speed claim.

---

## 6. Build target matters

`--dada2rs` is **required** — the harness does not auto-discover, because the
build target changes the numbers:

- `target/release` — portable, reproducible across nodes. Use for numbers others
  should be able to reproduce.
- `target/release-native` — built with `-C target-cpu=native` (the
  `release-native` profile). Best-case single-machine performance (AVX2/AVX-512
  on x86, full NEON on Apple Silicon), but **not portable** across
  microarchitectures.

For a head-to-head where R's C was compiled with the cluster's native gcc, prefer
`release-native` so the alignment kernel isn't handicapped to SSE2.

---

## 7. PacBio vs Illumina specifics

- **PacBio input is raw, primered reads.** `remove-primers` trims primers,
  orients (`--orient`), and applies the length/quality filters in one pass
  (mirrors R `removePrimers()` + `filterAndTrim()`). Supply `--primer-fwd` /
  `--primer-rev` (5'→3'; the reverse primer is reverse-complemented internally).
- **Alignment params must be consistent.** Pass `--band`, `--homo-gap-p`, and
  `--kmer-size` to **both** `learn-errors` and the denoise step (the harness does
  this). Mismatches trigger the §4 warning and subtly change results.
- **k-mer size.** On ~1.5 kb HiFi reads `k=5` makes the pre-screen a no-op (≈
  every pair is fully aligned → slower); `k=7` is the dada2-rs default. R fixes
  `KMER_SIZE = 5`, so a matched-`k` (k=5) run isolates kernel+threading speed,
  while a `k=7` run additionally shows the screen-effectiveness gain. Per issue
  #15 the final chimera-filtered ASVs are ~unchanged across `k`.

---

## 8. Quick reference — recommended workflow

1. **Build** the binary you intend to measure (`release` or `release-native`).
2. **Single mode, both stacks:** `--pool {true|false|pseudo} --run-r` →
   read the `HEADLINE` and `summary.csv`.
3. **Diagnose under-utilization:** look at the `cores` column; for a deeper view
   run one denoise with `--verbose` and read the `[dada] phase times` /
   `map parallel efficiency` lines.
4. **Quantify scaling:** `--thread-sweep` (thread scaling) and
   `--sample-jobs-sweep` (samples-in-flight, with `peak_rss`).
5. **Validate results:** `compare_asvs.py` / `compare_errors.R` before reporting.
6. **Memory-constrained / huge sample sets:** add `--low-memory` (pseudo) and
   compare peak RSS vs wall.
