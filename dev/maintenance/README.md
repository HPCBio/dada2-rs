# Maintenance scripts

Tooling a **maintainer** runs against the repository itself — releases, and
build-config guards. Nothing here is part of a user's analysis.

| Script | What it does |
|---|---|
| `check-build-sync.sh` | fails if `justfile` recipes and `Makefile` targets have drifted apart; runs in CI and from `.githooks/pre-commit` |
| `publish-crate.sh` | publishes an NW-only crate to crates.io, temporarily stripping the `wfa` git dependency that crates.io rejects (issue #63) |

Both `cd` to the repository root, so they can be invoked from anywhere; prefer
`just check-build-sync` / `just publish-crate` (or the `make` equivalents).

## Release checklist

crates.io freezes the README with each version, and a published version can be
yanked but never replaced, so everything a user reads must be right before step
5.

1. **Release PR** from a `release/vX.Y.Z` branch:
   - bump `version` in `Cargo.toml` and `CITATION.cff`, and refresh `Cargo.lock`
     (`cargo update -p dada2-rs --offline`);
   - rename `## [Unreleased]` in `CHANGELOG.md` to `## [X.Y.Z] - YYYY-MM-DD`,
     mark changes to results, flags or output **BREAKING**, and update the
     compare links at the bottom.
2. **Docs review**, in the same PR:
   - `README.md` against `dada2-rs --help` (subcommands and groups), the
     changelog's headline items, and every relative link (pages under `docs/`
     move);
   - `docs/index.md` and `docs/installation.md` for the same.
3. **Dry run:** `just publish-crate --dry-run` packages and verifies with no
   warnings and leaves the working tree clean. Check the file count against the
   last release; growth usually means a new directory needs an `exclude`.
4. **Tag:** after the PR merges and CI passes on the merge commit,
   `git tag -a vX.Y.Z <merge-sha> -m vX.Y.Z` and push it. The tag runs
   `release.yml` (binaries) and `docker.yml` (images); confirm both pass and the
   GitHub Release has its tarballs.
5. **Publish:** `just publish-crate` (needs `cargo login`), then check the
   crates.io page renders the README.
6. **Reopen:** add an empty `## [Unreleased]` section and its
   `compare/vX.Y.Z...HEAD` link to `CHANGELOG.md` on `main`, so new entries do
   not land in the released section (v0.2.0 was tagged without this).

## Why this is not in `scripts/`

`scripts/` is for helpers a **user** runs as part of an analysis — the plotting
scripts, the error-model converter, read tracking — and those are documented on
ReadTheDocs. Keeping the two apart makes that boundary checkable rather than a
convention: **everything in `scripts/` is documented, nothing in `dev/` is**
(issue #193).
