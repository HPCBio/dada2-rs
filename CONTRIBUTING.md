# Contributing to dada2-rs

Thank you for your interest in contributing! This project is an experimental Rust port of the [DADA2](https://benjjneb.github.io/dada2/) amplicon sequence denoising algorithm.

## Prerequisites

- A recent stable Rust toolchain ([rustup.rs](https://rustup.rs))
- Familiarity with the DADA2 algorithm is helpful; the original paper and the R source at [benjjneb/dada2](https://github.com/benjjneb/dada2) are the primary references

## Getting started

```bash
git clone https://github.com/HPCBio/dada2-rs
cd dada2-rs
cargo build --release
cargo test
```

## How to contribute

1. **Open an issue first** for any non-trivial change so we can discuss the approach before you invest time implementing it.
2. Fork the repository and create a branch from `main`.
3. Make your changes, keeping commits focused and the commit messages descriptive.
4. Run `cargo test` and `cargo clippy` before submitting.
5. Open a pull request against `main` with a clear description of what changed and why.

## Algorithmic fidelity

The primary goal of this project is a faithful port of R DADA2. When in doubt, behaviour should match the R reference implementation. If you find a divergence, please open an issue with a minimal reproducible example — per-raw p-value traces (via `DADA2_TRACE_PVAL`) are often useful for narrowing down where the outputs diverge.

## Code style

Follow standard Rust idioms (`cargo fmt`, `cargo clippy`). Comments should explain *why* something is done, not *what* — especially where the code mirrors a non-obvious algorithmic choice from the R or C++ source.

## Versioning and releases

This project uses [Semantic Versioning](https://semver.org/).

- **PATCH** releases are for backward-compatible fixes
- **MINOR** releases are for backward-compatible features
- **MAJOR** releases are for backward-incompatible changes, especially CLI,
  output-format, or workflow changes that users would need to adapt to

Version bumps are manual. When preparing a release:

1. Update the version in `Cargo.toml`
2. Update the version in `CITATION.cff`
3. Run `cargo test --release`, `cargo fmt --all -- --check`, and `cargo clippy`
4. Commit the release changes
5. Push a matching Git tag in the form `vX.Y.Z`

The CI workflow checks that `Cargo.toml` and `CITATION.cff` stay in sync, and
the release workflow verifies that the pushed tag matches those files before it
publishes a GitHub release.

## Reporting bugs

Please open a GitHub issue with:
- The command you ran
- The dada2-rs version (`dada2-rs --version`)
- A minimal input that reproduces the problem (a small FASTQ excerpt is ideal)
- Expected vs. observed output

## License

By contributing you agree that your contributions will be licensed under the GNU Lesser General Public License v3 (LGPL-3.0-only), the same license as the original DADA2 project.
