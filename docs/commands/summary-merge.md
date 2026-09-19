# `summary-merge`

Merge per-sample [`summary`](summary.md) JSONs into a run-level quality and
binning report.

```bash
dada2-rs summary summaries/*.json --report -o run_summary.json
```

Unions the observed quality-value distribution across all samples, so rare bins
missed by any single low-count sample — a seldom-called Q2, say — are recovered
at the run level. This matters because binned-quality detection on one shallow
sample can easily miss a level that the run as a whole does use.

## Input

**`<INPUTS>...`** — per-sample `summary` JSON files. Gzip is fine, and `-` reads
stdin.

## Metrics

**`--binned-threshold`** (default 8) — the maximum number of distinct quality
values for the run to be judged *binned*. Matches the `summary` default.

**`--expected-bins`** — declared bin levels to validate the run against, e.g.
`2,12,24,40`. When set, observed ⊆ expected is reported as "consistent"; any
observed value outside the declared set is a violation, which means either a
wrong declaration or a silently changed instrument scheme. Declared-but-unobserved
bins are reported as informational, not as errors.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON.

## Diagnostics

**`--report`** — also print a human-readable run-level report to stderr. Stdout
stays clean JSON.

## See also

- [`summary`](summary.md)
- [Binned quality scores](../findings/binned-quality.md)
