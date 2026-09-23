# Maintenance scripts

Tooling a **maintainer** runs against the repository itself — releases, and
build-config guards. Nothing here is part of a user's analysis.

| Script | What it does |
|---|---|
| `check-build-sync.sh` | fails if `justfile` recipes and `Makefile` targets have drifted apart; runs in CI and from `.githooks/pre-commit` |
| `publish-crate.sh` | publishes an NW-only crate to crates.io, temporarily stripping the `wfa` git dependency that crates.io rejects (issue #63) |

Both `cd` to the repository root, so they can be invoked from anywhere; prefer
`just check-build-sync` / `just publish-crate` (or the `make` equivalents).

## Why this is not in `scripts/`

`scripts/` is for helpers a **user** runs as part of an analysis — the plotting
scripts, the error-model converter, read tracking — and those are documented on
ReadTheDocs. Keeping the two apart makes that boundary checkable rather than a
convention: **everything in `scripts/` is documented, nothing in `dev/` is**
(issue #193).
