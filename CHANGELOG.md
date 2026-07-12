# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) (pre-1.0, so
minor versions may carry breaking changes).

## [0.2.0] - Unreleased

### Added
- `summary` subcommand: per-position quality metrics, cumulative expected-error
  metrics, and an optional per-read sequence-complexity histogram (#8).
- `sample` and `errors-from-sample` subcommands for subsampling input FASTQs and
  learning an error model from the subsample (usable for bootstrapping).
- `--failed-uniques` diagnostic on `dada` / `dada-pseudo` / `dada-pooled`:
  emits a TSV of uniques that failed to seed a partition (#60).
- `chimera-diagnostics` for higher-order (trimera) screening, with
  nearest-parent and flank gates (#54).
- Experimental **WFA alignment backend** (`--align-backend wfa2`) with a
  `--wfa-max-edits` edit-budget cap (#49, #51). ASV-equivalent to Needleman-
  Wunsch on tested Illumina and PacBio HiFi data, but not byte-identical.
- `dada-pooled --gzip` writes the per-sample JSON files gzip-compressed as
  `{sample}.json.gz` (#70). `dada` and `dada-pseudo` gain the same `--gzip`
  flag, for consistency across the per-sample denoisers.
- `merge-pairs` records input-file provenance and warns on mismatch (#10).
- `just` / `make` task runners for build, install, test, and docs (#46).

### Fixed
- `learn-errors` and `errors-from-sample` now gzip-compress their output when
  the `-o` path ends in `.gz`, for consistency with the other JSON-emitting
  subcommands (previously the file was written uncompressed) (#71).
- `kdist-calibrate --from-dada` now resolves derep inputs whose filenames were
  renamed by the pipeline (e.g. `{sample}.derep.R1.json.gz`): it falls back from
  the exact `{sample}.json[.gz]` match to a `{sample}.*.json[.gz]` scan,
  disambiguating an ambiguous prefix by the derep JSON's own `sample` field
  (#72).
- `kdist-calibrate --from-dada` now reads gzip-compressed `dada` output JSONs
  (`*.json.gz`, e.g. from `dada-pooled --gzip`); previously the gzip bytes were
  parsed as raw JSON and failed with `expected value, line 1 column 1`.
- `make-sequence-table` now reads gzip-compressed `dada` / `merge-pairs` inputs
  (`*.json.gz`, e.g. from `dada --gzip` / `dada-pooled --gzip`); previously the
  gzip bytes were parsed as raw JSON and failed with `expected value, line 1
  column 1`.
- Single-file JSON `-o` outputs now honor a `.gz` output path, routing through
  the shared gzip-aware writer instead of a raw write (matching `learn-errors`,
  #71): `summary`, `summary-merge`, `merge-pairs`, `remove-primers`,
  `filter-and-trim`, `make-sequence-table`, `remove-bimera-denovo`,
  `assign-taxonomy`, `assign-species`, and `dada`'s single-sample output.

### Changed
- The experimental WFA backend is now gated behind an **off-by-default `wfa`
  Cargo feature** (#63). Default builds — and the published crate — are
  Needleman-Wunsch only; selecting `--align-backend wfa2` without the feature
  errors rather than silently falling back to NW. Building WFA requires a source
  checkout (`cargo build --features wfa`) because it depends on a git crate that
  cannot ship on crates.io.
- Denoising/error-learning verbose logs now echo the active alignment backend
  (#51).
- Trimmed the published crate: development, benchmarking, concordance data
  (multi-MB PacBio FASTQs), examples, notes, and CI/infra files are excluded.

### Performance
- Sparse k-mer-8 screen, gated to k ≥ 8, cutting resident memory in the
  high-k regime (#43).
- `dada-pooled` streams dereplication into merging instead of holding all
  dereps in memory, and uses an integer merge-quality accumulator (#39, #41).
- Various k-mer-vector memory reductions in the alignment screen (#32).

### CI
- R DADA2 concordance guardrail for Illumina and PacBio, with abundance/recall/
  precision gates (#35).

## [0.1.0] - 2026-06-03

- Initial crates.io release: Rust ports of the core DADA2 pipeline steps
  (`filter-and-trim`, `derep`, `learn-errors`, `dada`, `merge-pairs`,
  `remove-bimera-denovo`, `assign-taxonomy` / `assign-species`, sequence-table
  helpers).

[0.2.0]: https://github.com/HPCBio/dada2-rs/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/HPCBio/dada2-rs/releases/tag/v0.1.0
