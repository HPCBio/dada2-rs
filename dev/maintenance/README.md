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
6. **Bioconda:** open the version-bump PR by hand (see below).
7. **Reopen:** add an empty `## [Unreleased]` section and its
   `compare/vX.Y.Z...HEAD` link to `CHANGELOG.md` on `main`, so new entries do
   not land in the released section (v0.2.0 was tagged without this).

## Bioconda version bump

The recipe is `recipes/dada2-rs/` in
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) (#325), built
from the crates.io tarball. Bioconda's autobump bot has no crates.io support,
so it never sees a new release and nothing opens on our side; every bump is a
manual PR until #330 moves the source to a GitHub Release asset. Bump PRs the
bot does open @mention the `recipe-maintainers` listed in `meta.yaml`.

1. In a clone of the `HPCBio/bioconda-recipes` fork, sync `master` with
   upstream and branch (`bump-dada2-rs-X.Y.Z`).
2. In `recipes/dada2-rs/meta.yaml`, set `version` and `sha256`, and reset
   `build: number` to 0. The sha256 is the `checksum` field of
   `https://crates.io/api/v1/crates/dada2-rs/X.Y.Z`.
3. **First bump after 0.5.0 only:** delete `build-tags-new-build-dir.patch`,
   `patch` from the build requirements, and the `patch -p1` line in
   `build.sh`. The 0.5.0 crate predates that fix (PR #329).
4. `bioconda-utils lint --packages dada2-rs` from the clone's root. To install
   it: `conda create -n bioconda -c conda-forge -c bioconda bioconda-utils`.
5. Push and open the PR. Bioconda's CI builds linux-64, linux-aarch64 and
   osx-arm64; the test asserts the exact `--version` string, so a wrong
   version or a leftover patch fails there.

## Why this is not in `scripts/`

`scripts/` is for helpers a **user** runs as part of an analysis — the plotting
scripts, the error-model converter, read tracking — and those are documented on
ReadTheDocs. Keeping the two apart makes that boundary checkable rather than a
convention: **everything in `scripts/` is documented, nothing in `dev/` is**
(issue #193).
