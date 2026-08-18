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
- The `nalign` and `nshroud` counters no longer overflow on large pooled runs
  (#143). Both count pairwise comparisons (`nraw × nclusters`) and were `u32`,
  so any pool exceeding ~4.3 billion comparisons wrapped silently: a 1.23 M-
  unique NovaSeq soil 16S pool reported 3,015,581,030 alignments against an
  actual 11,605,515,622. The wrapped values reached the `ALIGN:` verbose line
  *and* the `nalign` / `nshroud` fields of `dada`, `dada-pseudo` and
  `dada-pooled` output JSON, so affected artifacts carry wrong values in those
  two fields. Nothing reads either counter — no clustering, error-model,
  p-value or table output was affected, and no re-run is needed to correct an
  archived result. Smaller pools (under ~4.3e9 comparisons) were never
  affected.
- LOESS error models no longer collapse to a uniform `1e-7` floor on input with
  five or fewer populated quality columns (#95). The tricube kernel zeroes the
  farthest neighbour in the local window, which left a `span = 0.75`,
  `degree = 2` fit short of the coefficients it needed; the local fit now falls
  back to the highest degree the surviving points can identify. Affected
  `--errfun loess` (and the `pacbio` loess path) on binned-quality data, where a
  handful of distinct Q values is normal. Fits with six or more populated
  columns are byte-identical to before. `loess_errfun` and `pacbio_errfun` now
  return `Result` and report an error rather than emitting an unusable uniform
  matrix when every transition fails to fit.
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
- `Raw::e_minmax` moves off `Raw` into a dense `Vec<f64>` on `B`, parallel to
  `B::raw_cluster` (#147). `b_compare`'s serial store reads this value once per
  raw per call; `Raw` is 160 bytes, so as a field each read pulled a full
  64-byte cache line to use 8 bytes of it, and that strided read measured as
  **83–87% of the store scan**. Contiguous, the store falls from 18.8–19.1 to
  4.9–5.5 ns per raw-visit — **−71% of store time** — worth **−20.1 to −21.3%
  of `run_dada`** on NovaSeq soil 16S and −15.2% on ITS2, byte-identical on
  both pools across two replicates. Effective cores rise 31.4 → 38.4 (R1) and
  29.4 → 36.5 (R2) of 64.
- `b_compare`'s serial reduction is folded into the store loop it already
  shared a result vector with (#143): seven passes over that vector become
  one. Worth **−24 to −30% of `run_dada`** on two NovaSeq soil pools
  (16S and ITS2, both reads), byte-identical. The reduction had been six
  separate passes costing 337 s on soil 16S — 31% of `run_dada` — on a phase
  whose parallel map is already bandwidth-saturated and cannot absorb it.
- `b_shuffle` carries `compmax` across bud rounds instead of rebuilding it
  per bud (#139), deleting the build scan: 9,700 rebuilds collapse to 1.
  Measured **−45 to −48% of `run_dada`** on soil 16S and −15.4 to −15.9% on
  ITS2 (64 threads, post-#143), byte-identical across arms.
  `DADA2RS_SHUFFLE_NO_CARRY=1` restores the per-bud rebuild.
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
