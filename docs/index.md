# dada2-rs

An experimental implementation of [DADA2](https://benjjneb.github.io/dada2/) in
Rust. The goal is to reproduce the DADA2 amplicon workflow closely while
exploring paths to better performance — currently several times faster than R
DADA2 with a fraction of the memory.

!!! note "Provenance"
    `dada2-rs` is a reimplementation. If you use it, please cite the original
    DADA2 work (see `CITATION.cff` in the repository). This project follows the
    [rewrites.bio](https://rewrites.bio) guidelines as closely as possible.

## What's implemented

Rust ports of the core DADA2 functions:

| Step | DADA2 (R) | dada2-rs subcommand |
|---|---|---|
| Filter / trim | `filterAndTrim` | `filter-and-trim` |
| Primer removal | `removePrimers` | `remove-primers` |
| Dereplication | `derepFastq` | `derep` |
| Error models | `learnErrors` | `learn-errors` |
| Denoising | `dada` | `dada`, `dada-pooled`, `dada-pseudo` |
| Merging | `mergePairs` | `merge-pairs` |
| Sequence table | `makeSequenceTable` | `make-sequence-table` |
| Chimera removal | `removeBimeraDenovo` | `remove-bimera-denovo` |
| Taxonomy | `assignTaxonomy` / `assignSpecies` | `assign-taxonomy` / `assign-species` |

The three denoising subcommands mirror R's `pool` argument as distinct,
explicit commands:

- **`dada`** — per-sample (independent) denoising; accepts one or more inputs.
- **`dada-pooled`** — full pooling (R `pool = TRUE`).
- **`dada-pseudo`** — pseudo-pooling (R `pool = "pseudo"`).

## Building

Requires a recent stable Rust toolchain ([rustup.rs](https://rustup.rs)).

```bash
cargo build --release
# binary at target/release/dada2-rs

# best-case, machine-specific build (AVX2/AVX-512 / full NEON):
RUSTFLAGS="-C target-cpu=native" cargo build --profile release-native
# binary at target/release-native/dada2-rs
```

Run `dada2-rs <subcommand> --help` for full parameter documentation.

## Where to go next

- **[Installation](installation.md)** — build from source (standard and native
  release) or use a container image.
- **[Illumina MiSeq walkthrough](walkthrough-illumina.md)** — paired-end,
  end-to-end.
- **[PacBio HiFi walkthrough](walkthrough-pacbio.md)** — single-end, with primer
  removal and PacBio-tuned parameters.
- **[Performance & Benchmarking](benchmarking.md)** — the tooling, log output,
  and metrics for evaluating performance and running head-to-head comparisons
  against R DADA2.
- **[Citation](citation.md)** — how to cite the original DADA2 work.

!!! info "Documentation status"
    This site is being built out incrementally. Some background and reference
    material still lives in the
    [project README](https://github.com/HPCBio/dada2-rs#readme) and will migrate
    here over time.
