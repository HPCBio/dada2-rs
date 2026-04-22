# dada2-rs

This is an experimental implementation of DADA2 in Rust, using Claude Sonnet 4.6 for the bulk of the work. 

This work will be following [rewrites.bio](https://rewrites.bio), including full credit to the original work by Ben Callahan and other DADA2 contributors. 

## Plans

The short term plan: 

1) Port key code over to Rust: dereplication, error modeling, denoising, merging, chimera removal
2) Decouple any underlying C++ code from R and reimplement key steps as subcommands, R classes as Rust structs/classes, etc.
3) Allow intermediate outputs (in JSON) that can be evaluated for debugging purposes or for plotting in R, Python, etc.
4) Add basic regression tests that follow those within the original dada2 repository and expect results similar to those expected

We do anticipate porting taxonomic classification at a later point

## Building

Requires a recent stable Rust toolchain ([rustup.rs](https://rustup.rs)).

```bash
cargo build --release
# binary at target/release/dada2-rs
```

## Subcommands

| Subcommand | Description |
|---|---|
| `filter-and-trim` | Filter and trim FASTQ reads (mirrors `filterAndTrim`) |
| `derep` | Dereplicate a FASTQ file |
| `sample` | Dereplicate and subsample FASTQ files, one JSON per sample |
| `errors-from-sample` | Learn error model from derep JSON files |
| `learn-errors` | Learn error model directly from FASTQ files |
| `dada` | Denoise a sample using a learned error model |
| `merge-pairs` | Merge denoised forward and reverse reads |
| `make-sequence-table` | Build a sample × sequence count table |
| `remove-bimera-denovo` | Remove chimeric sequences |
| `summary` | Per-position quality metrics from a FASTQ |

Run `dada2-rs <subcommand> --help` for full parameter documentation.

## Typical MiSeq workflow

The steps below follow the [DADA2 MiSeq SOP](http://benjjneb.github.io/dada2/tutorial.html).
All intermediate outputs are JSON and can be inspected or plotted independently.

### 1. Filter and trim

```bash
dada2-rs filter-and-trim \
  --fwd  raw/sample_R1.fastq.gz --filt filtered/sample_R1.fastq.gz \
  --rev  raw/sample_R2.fastq.gz --filt-rev filtered/sample_R2.fastq.gz \
  --trunc-len 240 160 \
  --max-n 0 --max-ee 2 2 --trunc-q 2 \
  --compress --verbose
```

### 2. Learn the error model

**Option A — from pre-computed derep JSON files** (via `sample` subcommand):

```bash
# Dereplicate and subsample forward reads across all samples
dada2-rs sample filtered/*_R1.fastq.gz \
  --output-dir sample_json/ --nbases 100000000 --verbose

# Learn error model
dada2-rs errors-from-sample sample_json/*.json \
  --errfun loess -o errors_fwd.json --verbose
```

**Option B — directly from FASTQ** (single command):

```bash
dada2-rs learn-errors filtered/*_R1.fastq.gz \
  --nbases 100000000 --errfun loess \
  -o errors_fwd.json --verbose
```

#### Visualise cluster diagnostics during error learning

Pass `--diag-dir` to emit a `iter_NNN.json` file for each self-consistency
iteration, then plot with the bundled R script:

```bash
dada2-rs errors-from-sample sample_json/*.json \
  --errfun loess --diag-dir diag_fwd/ -o errors_fwd.json --verbose

Rscript comparison/plot_cluster_diag.R diag_fwd/ diag_fwd/cluster_diag.pdf
```

The plot shows cluster counts, birth-type breakdown, convergence trace, and
alignment work across iterations.

### 3. Denoise each sample

```bash
dada2-rs dada filtered/sample_R1.fastq.gz \
  --error-model errors_fwd.json --show-map \
  -o dada/sample_R1.json --verbose
```

Repeat for reverse reads using the reverse error model.

### 4. Merge paired reads

```bash
dada2-rs merge-pairs \
  --fwd-dada dada/fwd/*.json \
  --rev-dada dada/rev/*.json \
  --fwd-fastq filtered/fwd/*.fastq.gz \
  --rev-fastq filtered/rev/*.fastq.gz \
  -o merged.json --verbose
```

### 5. Build sequence table and remove chimeras

```bash
dada2-rs make-sequence-table merged.json -o seqtab.json
dada2-rs remove-bimera-denovo seqtab.json --method consensus -o seqtab_nochim.json
```

### Visualise the error model

```bash
Rscript plot_errors.R errors_fwd.json errors_fwd.pdf
```

## Comparing with R DADA2

The `comparison/` directory contains scripts for validating output against
the reference R implementation on the MiSeq SOP dataset:

```bash
# Run Rust pipeline (filter → derep → error model)
bash comparison/run_rust_errors.sh

# Compare error matrices against R's learnErrors()
Rscript comparison/compare_errors.R \
  comparison/comparison_out/F3D0_S188_L001_R1_001_filtered.fastq.gz \
  comparison/comparison_out/F3D0_S188_L001_R1_001_errors_rust.json
```

Requires R packages: `dada2`, `jsonlite`, `ggplot2`, `gridExtra`.

## AI Assistance Disclosure
This tool was written with the assistance of AI coding agents, specificall Claude Code, using Sonnet 4.6. All commits using AI are noted.

Correctness is validated by comparing output against DADA2 v1.36 on a suite of real sequencing datasets - not by manual code review alone. 
AI generated the implementation; humans defined the validation criteria, made some key coding updates, and verified results.
