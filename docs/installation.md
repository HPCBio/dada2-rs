# Installation

## Prerequisites

A recent stable Rust toolchain. The easiest way to get one is
[rustup.rs](https://rustup.rs):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

## Build from source

### Standard release

Portable and reproducible across machines — use this when others should be able
to rebuild and match your binary.

```bash
cargo build --release
# binary at target/release/dada2-rs
```

### Native release (recommended for production runs)

Built with `-C target-cpu=native` via the `release-native` profile, so LLVM uses
every instruction-set extension the build host supports (AVX2/AVX-512 on x86-64,
the full NEON feature set on Apple Silicon). This noticeably speeds up the
alignment kernel.

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --profile release-native
# binary at target/release-native/dada2-rs
```

!!! warning "Portability trade-off"
    A `release-native` binary only runs on CPUs at or above the build host's
    microarchitecture. Build it on (or for) the machine you will run it on. For
    a portable artifact, use the standard `release` build.

## Task runners (`just` / `make`)

Common build, test, and install tasks are wrapped by a `justfile` and an
equivalent `Makefile`. The `justfile` is the canonical source; the `Makefile`
is a portable mirror for hosts that have `make` but not `just`. Recipe and
target names match (a CI check enforces this), so use whichever you have:

```bash
just build         # or: make build         → release binary
just build-native  # or: make build-native  → native-CPU build
just check         # or: make check         → fmt-check + clippy + tests
just --list        # or: make help          → list available tasks
```

### Installing onto PATH

`install` copies the binary plus the user-facing helper scripts (the R plotting
scripts and a couple of Python utilities) into `$PREFIX/bin`. Helper scripts are
made executable and namespaced as `dada2-rs-<name>` to avoid PATH collisions —
for example `scripts/plot_errors.R` installs as `dada2-rs-plot-errors`.

```bash
PREFIX=/usr/local make install        # or: PREFIX=/usr/local just install
```

Overrides use the **env-var form** (`PREFIX=...`, `DESTDIR=...`), which behaves
identically for both `just` and `make`. `DESTDIR` supports staged installs for
packaging:

```bash
DESTDIR=/tmp/stage PREFIX=/usr make install   # installs into /tmp/stage/usr/bin
```

Development/benchmarking scripts live under `dev/` and are intentionally **not**
installed. To remove an installation, use `make uninstall` (same `PREFIX`).

## Docker

Dedicated container images are published for the project — see the `Dockerfile`
and the Docker build workflow in the repository for the current image
coordinates and tags.

## Verify

```bash
target/release/dada2-rs --help
target/release/dada2-rs --version
```

You should see the list of [subcommands](index.md#whats-implemented). Run
`dada2-rs <subcommand> --help` for the full parameter set of any step.

## Next steps

- [Illumina MiSeq walkthrough](walkthrough-illumina.md)
- [PacBio HiFi walkthrough](walkthrough-pacbio.md)
- [Performance & Benchmarking](benchmarking.md)
