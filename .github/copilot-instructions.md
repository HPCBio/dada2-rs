# Copilot Instructions for `dada2-rs`

## Build, test, and lint commands

- Use the same commands as CI unless there is a specific reason not to:
  - `cargo build --release`
  - `cargo test --release`
  - `cargo fmt --all -- --check`
  - `cargo clippy`
- For a quick local check, `cargo test` is also documented in `CONTRIBUTING.md`, but CI uses the `--release` variants for build and test.
- To run a single test, pass the test name substring, for example:
  - `cargo test --release test_external_errfun_roundtrip`
  - `cargo test --release test_vectorized_vs_endsfree_one_gap`
- There is one ignored benchmark-style test in `src/kmers.rs`; run it explicitly with:
  - `cargo test --release -- --ignored bench_kmer_dist8 --nocapture`

## High-level architecture

- `src/main.rs` is the only binary entry point. It dispatches all CLI subcommands defined in `src/cli.rs`.
- The main amplicon-processing pipeline is staged through JSON outputs:
  1. `filter-and-trim` in `src/filter_trim.rs` filters FASTQ reads.
  2. `derep` in `src/derep.rs` collapses reads into unique sequences, mean qualities, and an optional read-to-unique map.
  3. `learn-errors` / `errors-from-sample` in `src/learn_errors.rs` iteratively fit error matrices by repeatedly running the core DADA algorithm and re-estimating transition rates with `src/error_models.rs`.
  4. `dada` in `src/dada.rs` denoises one sample using a learned error model.
  5. `merge-pairs` in `src/merge_pairs.rs` re-dereplicates FASTQ files, composes those mappings with `dada` output maps, and merges forward/reverse ASVs.
  6. `make-sequence-table` in `src/sequence_table.rs` converts `dada` or `merge-pairs` JSON into a sample-by-sequence count matrix.
  7. `remove-bimera-denovo` in `src/remove_bimera.rs` filters chimeras from a sequence table.
- The DADA core is split across several low-level modules:
  - `src/containers.rs` holds the Rust equivalents of the original `Raw`, cluster, and partition structures.
  - `src/cluster.rs`, `src/pval.rs`, `src/nwalign.rs`, `src/kmers.rs`, and `src/error.rs` implement comparison, budding, p-value, alignment, k-mer screening, and auxiliary error outputs used by `src/dada.rs`.
- Taxonomy is a separate path from denoising:
  - `assign-taxonomy` and `assign-species` are wired in `src/main.rs`.
  - `src/taxonomy.rs` implements the naive-Bayes k-mer classifier and bootstrap logic.
- `comparison/` contains scripts used to compare Rust outputs against the reference R DADA2 implementation, especially around error-model learning.

## Key conventions specific to this repository

- Prioritize algorithmic fidelity to the R/C++ DADA2 reference implementation. Many functions are documented with the exact R or C++ counterpart they mirror; preserve that correspondence when changing behavior.
- All pipeline JSON outputs are wrapped with `dada2_rs_command` and `dada2_rs_version` via `misc::Tagged`. Downstream commands validate those tags before deserializing, so keep this contract intact when adding or changing output formats.
- Dereplication order is not arbitrary: unique sequences are sorted by abundance descending, with stable first-seen tie handling. `merge-pairs` depends on that deterministic indexing when it re-dereplicates FASTQ inputs and joins them back to `dada` maps.
- `merge-pairs` requires `dada` JSON produced with `--show-map`; without the `map` field, merging is intentionally rejected.
- Error-model JSON embeds the parameters used to learn the model (`LearnedErrParams`). Keep those fields aligned with the actual inference parameters so `dada --inherit-err-params` and mismatch warnings remain meaningful.
- Readers in `misc.rs` transparently accept `.gz` FASTQ, FASTA, and JSON files based on filename extension. Reuse those helpers instead of adding parallel I/O paths.
- When comments are needed, explain the algorithmic reason or a fidelity constraint rather than restating the code.
