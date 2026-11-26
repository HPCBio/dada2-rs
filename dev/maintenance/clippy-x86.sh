#!/usr/bin/env bash
# clippy-x86.sh — lint the crate as x86-64 from a non-x86 host (#336).
#
# The AVX2 builds of the NW kernels and their `*_variants_match` tests are
# `#[cfg(target_arch = "x86_64")]`, so an arm64 build never compiles them and a
# mistake there first shows up on CI's ubuntu runner. This runs clippy for an
# x86-64 target instead. It checks and lints only: nothing is linked or run
# (Rosetta does not expose AVX2, so the variant tests could not run anyway).
#
# Needs rustup with the target's std installed. Homebrew's `cargo`/`rustc` ship
# no x86 std, so the rustup toolchain's bin is put first on PATH.
#
#   just clippy-x86        (or: make clippy-x86)
#   FEATURES=wfa just clippy-x86
set -euo pipefail
cd "$(dirname "$0")/../.."

if [ "$(uname -m)" = x86_64 ]; then
    echo "clippy-x86: host is already x86-64; \`just clippy\` covers it." >&2
    exit 0
fi

case "$(uname -s)" in
    Darwin) target=x86_64-apple-darwin ;;
    *) target=x86_64-unknown-linux-gnu ;;
esac
toolchain=${RUSTUP_TOOLCHAIN:-stable}

command -v rustup >/dev/null || {
    echo "clippy-x86: needs rustup (https://rustup.rs)" >&2
    exit 1
}
sysroot=$(rustup run "$toolchain" rustc --print sysroot)
if [ ! -d "$sysroot/lib/rustlib/$target" ]; then
    echo "clippy-x86: the $toolchain toolchain has no $target std. Install it with:" >&2
    echo "  rustup target add --toolchain $toolchain $target" >&2
    exit 1
fi

feat=()
[ -n "${FEATURES:-}" ] && feat=(--features "$FEATURES")

# A separate target dir: the x86 artifacts would otherwise evict the host ones.
PATH="$sysroot/bin:$PATH" CARGO_TARGET_DIR=target/clippy-x86 \
    cargo clippy --target "$target" --all-targets ${feat[@]+"${feat[@]}"} -- -D warnings
echo "clippy-x86: $target clean"
