# dada2-rs

📖 **Documentation:** [dada2-rs.readthedocs.io](https://dada2-rs.readthedocs.io) — installation, Illumina/PacBio walkthroughs, and the performance/benchmarking reference.

An experimental implementation of DADA2 in Rust, using Claude Code (specically Sonnet 4.6 and Opus 4.6/4.7/4.8) for the bulk of the work. 

## Implementations

Rust ports of:

  | Step | DADA2 (R) | dada2-rs |
  |---|---|---|
  | Filter/trimming FASTQ | `filterAndTrim` | `filter-and-trim` |
  | Filter/trimming FASTQ (PacBio) | `removePrimers`,`filterAndTrim` | `remove-primers` (one step) |
  | Dereplication | `derepFastq` | `derep` |
  | Error models | `learnErrors` | `learn-errors` |
  | Denoising | `dada` | `dada` |
  | Denoising, pooled | `dada(pool = TRUE)` | `dada-pooled` |
  | Denoising, pseudo-pooled | `dada(pool = "pseudo")` | `dada-pseudo` |
  | Merging | `mergePairs` | `merge-pairs` |
  | Chimera removal | `removeBimeraDenovo` | `remove-bimera-denovo` |
  | RDP taxonomic classifier | `assignTaxonomy` + `assignSpecies` | `assign-taxonomy` + `assign-species` |
  | Merging sequence tables | `mergeSequenceTables` | `make-sequence-table` (accepts multiple inputs) |
  | Making sequence tables from multiple inputs | `makeSequenceTable` | `make-sequence-table` |

Denoising reproduces R DADA2's ASV tables cell for cell on the CI concordance
fixtures; see [benchmarking](docs/benchmarking.md) and the
[changelog](CHANGELOG.md).

- Current error models:
  - `loessErrfun`
  - `PacBioErrfun`
  - `noqualErrfun`
  - `makeBinnedErrfun`
  - Custom error models in R and/or Python - *Experimental, tested*
- Other functionality:
  - Helper functions to convert sequence tables to TSV or FASTA (`seq-table-to-fasta`, `seq-table-to-tsv`)
  - Helper function to convert taxonomic output to TSV (`tax-to-tsv`)
  - Basic scripts and examples for comparing results between runs to trace
    differences
  - Experimental sub-sampling function for input FASTQ (`sample`) and related error model function (`errors-from-sample`) that mirrors `learn-errors` (can be used for bootstrapping)
  - Intermediate outputs (in JSON) - can be evaluated for debugging purposes or for plotting in R, Python, etc.
  - Experimental [wavefront-alignment (WFA) backend](docs/parameters.md#experimental-wfa-alignment-backend) (`--align-backend wfa2`), a developer build gated behind `--features wfa` (not in the published crate)
  - Per-position quality summaries, with binned-quality detection (`summary`, `summary-merge`)
  - Dedicated Docker builds available

## Building & installing

Requires a recent stable Rust toolchain ([rustup.rs](https://rustup.rs)).

```bash
cargo build --release
# binary at target/release/dada2-rs
```

Common tasks are wrapped by a `justfile` (canonical) and an equivalent
portable `Makefile` — use whichever you have; targets/recipes match:

```bash
just build            # or: make build      → release binary
just check            # or: make check      → fmt-check + clippy + tests

# install the binary + helper R/Python scripts onto PATH; scripts are
# namespaced as dada2-rs-<name> (e.g. dada2-rs-plot-errors). Use the
# env-var form for overrides — it works for both just and make:
PREFIX=/usr/local make install
```

For the native (`-C target-cpu=native`) build and Docker, see
**[Installation](docs/installation.md)**.

## Subcommands

Grouped as in `dada2-rs --help`; each subcommand has a page under
[`docs/commands/`](docs/commands/index.md).

| Group | Subcommands |
|---|---|
| Workflow | `remove-primers`, `filter-and-trim`, `derep`, `learn-errors`, `dada`, `dada-pooled`, `dada-pseudo`, `merge-pairs`, `make-sequence-table`, `remove-bimera-denovo`, `assign-taxonomy`, `assign-species` |
| Quality control | `summary`, `summary-merge` |
| Export | `seq-table-to-tsv`, `seq-table-to-fasta`, `tax-to-tsv` |
| Diagnostics | `sample`, `errors-from-sample`, `chimera-diagnostics`, `kdist-calibrate`, `reference-eval` |

## Usage & walkthroughs

End-to-end examples live in the documentation:

- **[Illumina MiSeq walkthrough](docs/walkthroughs/walkthrough-illumina.md)** (paired-end)
- **[PacBio HiFi walkthrough](docs/walkthroughs/walkthrough-pacbio.md)** (single-end, primer removal + PacBio-tuned params)
- Binned quality scores: **[Illumina](docs/walkthroughs/walkthrough-illumina-binned.md)** and **[PacBio](docs/walkthroughs/walkthrough-pacbio-binned.md)**

Run `dada2-rs <subcommand> --help` for the full parameter set of any step.

## Benchmarks & comparison with R DADA2

Performance is benchmarked head-to-head against R DADA2 with the harness in
`dev/benchmark/`. See the docs for the tooling, metrics, and results:

- **[Performance — tooling & metrics](docs/benchmarking.md)** (the harness,
  `cores`/`cpu_s`/peak-RSS, scaling sweeps, built-in logs, concordance checks)
- **[Benchmark results](docs/results.md)** (head-to-head scorecards by platform
  and pooling mode)

## AI Assistance Disclosure

See **[the project overview](docs/about.md)** for the project's origins and goals, which are directly relevant here.

This tool was written with the assistance of AI coding agents, specifically Claude Code, using Sonnet and Opus LLMs. All commits using AI assistance are openly noted.

Correctness is validated by comparing output against R DADA2 on a suite of real sequencing datasets - not by manual code review alone. Most of that testing used DADA2 v1.36.0 (R v4.5.1) on the cluster; local comparisons use DADA2 v1.40.0 (R v4.6.0).
AI generated the implementation; humans defined the validation criteria, made some key coding updates, and verified results. 

## Citation

`dada2-rs` is a reimplementation; if you use it, please cite the original work that
describes the algorithm (see [`CITATION.cff`](CITATION.cff)):

* Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. DADA2:
  High-resolution sample inference from Illumina amplicon data. *Nature Methods*.
  2016;13:581-583. doi:[10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
* Rosen MJ, Callahan BJ, Fisher DS, Holmes SP. Denoising PCR-amplified metagenome
  data. *BMC Bioinformatics*. 2012;13:283.
  doi:[10.1186/1471-2105-13-283](https://doi.org/10.1186/1471-2105-13-283)
