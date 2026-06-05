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
