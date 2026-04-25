# dada2-rs

An experimental implementation of DADA2 in Rust, using Claude Code (specically Sonnet 4.6 and Opus 4.6/4.7) for the bulk of the work. 

## Current state

The implementation is roughly 4.5x faster than the R DADA2 implementation (not including R loading overhead), and uses about 1/8 the memory. See also the Project Background for how this started and how we intend to move forward on this implementation.

## Plans

The short term plan: 

- [ ] Port key code over to Rust: dereplication, error modeling, denoising, merging, chimera removal
- [ ] Decouple any underlying C++ code from R and reimplement key steps as subcommands, R classes as Rust structs/classes, etc.
- [ ] Intermediate outputs (in JSON) that can be evaluated for debugging purposes or for plotting in R, Python, etc.
- [ ] Add basic regression tests for each stage that follow those within the original DADA2 repository

Long-term plans:

- [ ] Assess and optimize critical steps such as `learn-errors` and `dada`, in particular alignment steps within `dada` which are currently slower than the DADA2 implementation
- [ ] Porting and optimizing taxonomic classification
- [ ] Optimizing chimera detection
- [ ] Add some functionality to R or Python to allow custom error model generation

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

## Project background

The original purpose of the project was exploratory: walk through the steps in the DADA2 workflow to understand the underlying implementation for each key command, initially to replicate results from the R DADA2 workflow but to also explore potential paths for improving the implementation, as we use this quite extensively in our own research. 

This is now at an interesting stage, as the current implementation up to the 'learn errors' error model step, closely matches results from DADA2 but runs about 5x faster with 8x less memory. 

### Learning

This was also meant to be a learning opportunity on several fronts. I have been learning Rust over the last few years in my (vanishingly small) spare time, and actually planned a general port of DADA2 a few years back. I had programmed various languages over the years (Perl, Python, R, C++, and a little Java), but apart from R and Python I'm pretty, um, rusty, 

### AI

I have noticed a dramatic improvement in coding-based AI tools and agents, also noted by many others in the field. However there has been a lot of discussion on the best approaches to reimplementation strategies involving these approaches, and many controversies along the way. 

I've long been involved in open-source development and open-science projects, and I'm also the leader of a bioinformatics core group. I understand the community responses to this as well as the potential benefits to this community, and I believe community standards are needed that ensure we follow some general guidelines that ensure some degree of consistency and support, provenance that ensure the original implementors retain clear credit for their work, and where community members can play an active role.

### Guidelines

AI is having a clear, fundamental, and disruptive impact, and simply ignoring this is to one's detriment. Some of us are also mentors, and as such we have an obligation to understand this wildly changing landscape and help prepare students, postdocs, and scientist on how to best use AI for their future career path. However many controveries exist over the use of these tools, in some cases arguably crossing moral and ethical boundaries. Standards are sorely needed.

Thankfully, within the bioinformatics community these are starting to coalesce, for example [rewrites.bio](https://rewrites.bio). Therefore, this project will follow [rewrites.bio](https://rewrites.bio) guidelines as closely as possible to: 

* Ensure the original DADA2 developers and contributors are acknowledged,
* Follows the original implementation's details,
* Include tests and benchmarks for the work,
* Utilize consistent libraries where possible,
* Release as open source and follow the original licensing
* Should interest arise: develop a community that can contribute.

#### Caveat

One point where this implementation will vary from the rewrites.bio standards: due to some key implementation details (conversion of R/C++ to Rust including error models), results will vary slightly. However we will strive to reproduce results as closely as possible, including the future possibility to add bridging code for Rust/R (and maybe Rust/Python) for custom error model analysis to more closely emulate what we see from the original implementation. 

We also do not want to prevent additional outcomes or functionality that may come from the work in this project by being constrained to emulating the original code. For example, one interesting side benefit for PacBio HiFi reads has come from exposing the k-mer size as an option: alternative k-mer lengths appear to improve performance for PacBio denoising; this is something that needs to be explored more but could result in a substantial improvement in processing. 