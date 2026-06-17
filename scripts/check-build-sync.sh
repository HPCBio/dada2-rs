#!/usr/bin/env bash
# check-build-sync.sh — guard against the justfile (canonical) and Makefile
# (portable mirror) drifting apart at the recipe/target level.
#
# It compares the SET OF NAMES, not bodies: the realistic failure mode for thin
# task wrappers is adding a recipe to one file and forgetting the other. Body
# divergence in one-liners is low-risk and not worth a just->make translator.
#
# Names that legitimately differ between the tools are normalized away below
# (just's `default` lister vs make's `all`/`help`).
set -euo pipefail

cd "$(dirname "$0")/.."

# Names present in one tool but not meaningfully comparable in the other.
IGNORE='^(default|all|help)$'

just_names=$(just --summary | tr ' ' '\n' | grep -vE "$IGNORE" | sort -u)

# Makefile targets: lines like `name:` (not variable assignments, not patterns).
make_names=$(grep -E '^[a-zA-Z0-9_-]+:' Makefile \
    | sed -E 's/:.*//' \
    | grep -vE "$IGNORE" | sort -u)

if [ "$just_names" != "$make_names" ]; then
    echo "ERROR: justfile and Makefile have drifted apart." >&2
    echo "--- only in justfile ---" >&2
    comm -23 <(echo "$just_names") <(echo "$make_names") >&2
    echo "--- only in Makefile ---" >&2
    comm -13 <(echo "$just_names") <(echo "$make_names") >&2
    exit 1
fi

echo "OK: justfile recipes and Makefile targets are in sync."
