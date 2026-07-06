# justfile — CANONICAL source of project build/dev/deploy tasks for dada2-rs.
#
# This file is authoritative. The Makefile is a hand-maintained portable mirror
# (deployment hosts are guaranteed to have `make`, not `just`). When you change
# a task here, mirror it in the Makefile — `just check-build-sync` (run in CI /
# the pre-commit hook) fails if the two drift apart at the recipe/target level.
#
#   just build                          release binary
#   just build-native                   native-CPU build (non-portable, fastest)
#   PREFIX=/opt/dada2 just install      binary + helper scripts onto PATH
#   just check                          fmt-check + clippy + tests
#
# Note the env-var form for overrides (PREFIX=... DESTDIR=...): it works
# identically for both `just` and `make`, unlike make's `make install PREFIX=...`
# trailing-arg form which `just` does not accept.
#
# Helper R/Python scripts install namespaced as `dada2-rs-<name>` (e.g.
# scripts/plot_errors.R -> dada2-rs-plot-errors) to avoid PATH collisions.

bin      := "dada2-rs"
prefix   := env_var_or_default("PREFIX", "/usr/local")
destdir  := env_var_or_default("DESTDIR", "")

# list recipes (default)
default:
    @just --list

# release binary
build:
    cargo build --release

# native-CPU build (cluster/in-house; non-portable binary)
build-native:
    RUSTFLAGS="-C target-cpu=native" cargo build --profile release-native

# cross-compile Linux ARM64 (requires `cross`)
build-arm64:
    cross build --release --target aarch64-unknown-linux-gnu

# copy binary + helper scripts onto PATH (scripts made executable)
install: build
    #!/usr/bin/env bash
    set -euo pipefail
    bindir="{{destdir}}{{prefix}}/bin"
    install -d "$bindir"
    install -m 755 target/release/{{bin}} "$bindir/{{bin}}"
    for s in scripts/*.R scripts/*.py; do
        name=$(basename "$s" | sed -E 's/\.(R|py)$//; s/_/-/g')
        echo "install -m 755 $s $bindir/{{bin}}-$name"
        install -m 755 "$s" "$bindir/{{bin}}-$name"
    done

# remove installed binary + helper scripts
uninstall:
    #!/usr/bin/env bash
    set -euo pipefail
    bindir="{{destdir}}{{prefix}}/bin"
    rm -f "$bindir/{{bin}}"
    for s in scripts/*.R scripts/*.py; do
        name=$(basename "$s" | sed -E 's/\.(R|py)$//; s/_/-/g')
        echo "rm -f $bindir/{{bin}}-$name"
        rm -f "$bindir/{{bin}}-$name"
    done

# fmt-check + clippy (warnings as errors) + tests
check: fmt-check clippy test

fmt:
    cargo fmt

fmt-check:
    cargo fmt --check

clippy:
    cargo clippy --all-targets -- -D warnings

test:
    cargo test

# live-reload mkdocs site
docs-serve:
    mkdocs serve

# build static mkdocs site
docs-build:
    mkdocs build

clean:
    cargo clean

# publish an NW-only crate to crates.io (strips the WFA git dep; see issue #63).
# Pass --dry-run to verify without uploading: `just publish-crate --dry-run`.
publish-crate *ARGS:
    ./scripts/publish-crate.sh {{ARGS}}

# fail if justfile recipes and Makefile targets have drifted apart
check-build-sync:
    ./scripts/check-build-sync.sh
