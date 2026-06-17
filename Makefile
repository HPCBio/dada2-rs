# Makefile — portable MIRROR of the justfile (the canonical task source).
#
# The justfile is authoritative; this file exists so deployment hosts that have
# `make` but not `just` get the same tasks. Keep the two in step: `make
# check-build-sync` (CI / pre-commit) fails if recipe and target names drift.
#
# Common entry points:
#   make build                       release binary
#   make build-native                native-CPU build (non-portable, fastest)
#   make install PREFIX=/opt/dada2   binary + helper scripts onto PATH
#   make check                       fmt-check + clippy + tests
#
# Helper R/Python scripts are installed namespaced as `dada2-rs-<name>` (e.g.
# scripts/plot_errors.R -> dada2-rs-plot-errors) to avoid PATH collisions.

PREFIX  ?= /usr/local
DESTDIR ?=
BINDIR   = $(DESTDIR)$(PREFIX)/bin

BIN     = dada2-rs
SCRIPTS = $(wildcard scripts/*.R) $(wildcard scripts/*.py)

.PHONY: all build build-native build-arm64 install uninstall \
        check fmt fmt-check clippy test clean docs-serve docs-build \
        check-build-sync help

all: build

## build: release binary -> target/release/$(BIN)
build:
	cargo build --release

## build-native: native-CPU build (cluster/in-house; non-portable binary)
build-native:
	RUSTFLAGS="-C target-cpu=native" cargo build --profile release-native

## build-arm64: cross-compile Linux ARM64 (requires `cross`)
build-arm64:
	cross build --release --target aarch64-unknown-linux-gnu

## install: copy binary + helper scripts onto PATH (scripts made executable)
install: build
	install -d "$(BINDIR)"
	install -m 755 target/release/$(BIN) "$(BINDIR)/$(BIN)"
	@for s in $(SCRIPTS); do \
		name=$$(basename "$$s" | sed -E 's/\.(R|py)$$//; s/_/-/g'); \
		dest="$(BINDIR)/$(BIN)-$$name"; \
		echo "install -m 755 $$s $$dest"; \
		install -m 755 "$$s" "$$dest"; \
	done

## uninstall: remove installed binary + helper scripts
uninstall:
	rm -f "$(BINDIR)/$(BIN)"
	@for s in $(SCRIPTS); do \
		name=$$(basename "$$s" | sed -E 's/\.(R|py)$$//; s/_/-/g'); \
		dest="$(BINDIR)/$(BIN)-$$name"; \
		echo "rm -f $$dest"; \
		rm -f "$$dest"; \
	done

## check: fmt-check + clippy (warnings as errors) + tests
check: fmt-check clippy test

fmt:
	cargo fmt

fmt-check:
	cargo fmt --check

clippy:
	cargo clippy --all-targets -- -D warnings

test:
	cargo test

## docs-serve: live-reload mkdocs site
docs-serve:
	mkdocs serve

## docs-build: build static mkdocs site
docs-build:
	mkdocs build

clean:
	cargo clean

## check-build-sync: fail if justfile recipes and Makefile targets have drifted
check-build-sync:
	./scripts/check-build-sync.sh

## help: list documented targets
help:
	@grep -E '^## ' $(MAKEFILE_LIST) | sed -E 's/^## /  /'
