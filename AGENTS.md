# Agent Guide — dada2-rs

Guidance for AI agents (Claude Code, Copilot, etc.) working in this repository.

## What this project is

A Rust port of the [DADA2](https://benjjneb.github.io/dada2/) amplicon sequence denoising algorithm. **Algorithmic fidelity to the R reference is the primary goal.** When behaviour differs from R DADA2, that is a bug unless explicitly documented.

The R source lives at `https://github.com/benjjneb/dada2` (key files: `R/dada.R`, `R/errorModels.R`, `R/chimeras.R`, `R/taxonomy.R`).

## Build and test

```bash
cargo build --release          # binary → target/release/dada2-rs
cargo test                     # all unit + integration tests
cargo clippy -- -D warnings    # must be clean before committing
cargo fmt --all -- --check     # must be clean before committing
```

The pre-commit hook at `.githooks/pre-commit` runs both checks automatically. Activate it once after cloning:

```bash
git config core.hooksPath .githooks
```

## Repository layout

| Path | Purpose |
|---|---|
| `src/main.rs` | CLI dispatch — one handler block per subcommand |
| `src/cli.rs` | Clap structs for every subcommand |
| `src/error_models.rs` | Error rate estimation (loess, noqual, binned, external) |
| `src/dada.rs` | Core denoising algorithm |
| `src/learn_errors.rs` | Error learning (drives `learn-errors` / `errors-from-sample`) |
| `src/merge_pairs.rs` | Paired-end merging |
| `src/filter_trim.rs` | Per-read filtering and trimming |
| `src/summary.rs` | Per-position quality summaries (with full quality histograms) |
| `src/taxonomy.rs` | RDP-style k-mer taxonomic classification |
| `src/remove_bimera.rs` | Chimera detection and removal |
| `src/sequence_table.rs` | Sample × ASV count table construction |
| `scripts/` | Helper scripts: `track_reads.py`, `plot_quality_profile.R` |
| `testdata/` | Small FASTQ fixtures for regression tests |
| `build.rs` | Computes `DADA2_RS_VERSION_FULL` (semver + git SHA on non-tag builds) |

## JSON output conventions

All subcommand outputs are flat JSON objects tagged with `dada2_rs_command` and `dada2_rs_version` at the top level — **not** a `{"type": ..., "data": ...}` envelope. This is produced by the `Tagged<T>` wrapper with `#[serde(flatten)]` in `src/main.rs`. Scripts and tests must not assume an envelope.

## Versioning

`DADA2_RS_VERSION_FULL` is computed in `build.rs`:
- Tagged release (`vX.Y.Z`) → bare `X.Y.Z`
- Any other commit → `X.Y.Z-<8-char SHA>`
- Docker / CI builds inject it via the `DADA2_RS_VERSION_FULL` environment variable so that `.git` is not required inside the container

## Code conventions

- Follow standard Rust idioms. Run `cargo fmt` before committing.
- Comments explain *why*, not *what*. Avoid restating what the code already says.
- Do not add error handling for cases that cannot happen inside Rust-controlled code. Validate only at system boundaries (CLI input, external file I/O, external process output).
- Do not introduce abstractions speculatively. Three similar blocks are preferable to a premature helper that obscures intent.
- Algorithmic code that mirrors R or C++ DADA2 should cite the corresponding R function in a short doc comment.

## Testing guidance

- The external-process tests in `error_models::tests` create temp directories using a global `AtomicU64` counter (not timestamps) to guarantee uniqueness across parallel test threads. Follow the same pattern when writing tests that involve temp files.
- When adding a new subcommand or changing output shape, add or update the corresponding integration test in `src/main.rs` or the relevant module's `#[cfg(test)]` block.
- Run `cargo test` with the default parallelism; tests must not share mutable global state.

## Docker

```bash
# Local build (version injected at build time):
docker build \
  --build-arg DADA2_RS_VERSION_FULL=$(grep -m1 '^version = ' Cargo.toml | cut -d'"' -f2)-local \
  -t dada2-rs:dev .
```

CI builds (`.github/workflows/docker.yml`) compute the version from `GITHUB_REF` and inject it automatically.

## What to avoid

- Do not change numerical outputs of the denoising, error-model, or chimera-removal code without first verifying that the change matches R DADA2's output on the MiSeq SOP test data.
- Do not remove the `map` field from `DadaOutput` / `DadaPooled` output — it is always emitted and consumed by downstream tools.
- Do not add `--show-map` or similar flags to re-enable suppressed output; the map is unconditional.
- Do not alter the flat JSON output shape; downstream scripts (`track_reads.py`, `plot_quality_profile.R`) depend on it.
