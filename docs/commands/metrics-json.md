# Run metrics (`--metrics-json`)

`dada`, `dada-pooled` and `dada-pseudo` can write a structured record of what a
run did and how long each part took:

```bash
dada2-rs dada-pooled derep/*.json.gz --error-model err.json \
  -o dada/ --threads 24 --metrics-json run_metrics.json
```

**This is the only place the attribution detail lives.** `--verbose` used to
print it as ~130 lines of `ns/comp` tables; as of #162 it prints the run's shape
and results and points here instead. The compare attribution and split, map
parallel efficiency, the shuffle phases and scan split, bud redundancy,
p-update churn and the optimisation projections are all here and nowhere else.

## The two measurement levels

Measurement is **separate from verbosity**. `--verbose` decides what is
*printed*; these flags decide what is *measured*. They are separate because the
two cost wildly different amounts.

| Level | How to get it | Cost |
|---|---|---|
| `phases` | `--metrics-json <path>` | Effectively free |
| `attribution` | add `--metrics-attribution` | **Slows the run** |

**`phases` is free** and safe to leave on for production runs. It carries the
phase wall times, the store-loop counters, the `screened` / `aligned`
denominators, the resident footprint and the minimizer index decision. All of
these are either a handful of `Instant::now()` calls per bud round, or counters
that the store loop already folds in.

**`attribution` is not free.** It adds 2–4 `Instant::now()` calls to *every
comparison*, and comparisons reach 1.4e10 on a diverse pool. It buys `busy`, the
map's parallel efficiency, and the screen / DP / `al2subs` split. Never enable it
on a run whose wall time you intend to report.

`--metrics-attribution` requires `--metrics-json`, so the timers can never run
with nowhere to write.

!!! note "Absent is not zero"
    Fields that a level did not measure are **omitted**, not written as `0.0`.
    `compare.split` is absent at the `phases` level, and `compare.attribution`
    is absent on a single-threaded run, because the serial `b_compare` path
    returns no timing at all. Check `measure_level` before concluding something
    was free.

Measurement never changes results: the ASV output is byte-identical at every
level.

## Shape

```json
{
  "schema_version": 1,
  "measure_level": "phases",
  "wall_seconds": 812.4,
  "pipeline": { "derep": 41.2, "merge": 18.9, "dada": 701.1, "output": 51.2 },
  "runs": [
    {
      "sample": "__pooled__",
      "run":     { "nraw": 2.4e6, "threads": 24, "screen_backend": "kmer", "...": null },
      "phases":  { "compare": 640.2, "shuffle": 48.1, "bud": 9.7, "...": null },
      "compare": { "total": 640.2, "attribution": { "...": null }, "screened": 0, "aligned": 0 },
      "shuffle": { "build": 12.1, "reconcile": 31.8, "move_pass": 2.2, "...": null },
      "footprint": { "nraw": 2400000, "screen_repr": "dense", "...": null },
      "index":   { "use_index": true, "...": null }
    }
  ]
}
```

One entry in `runs` per `run_dada` invocation:

- **`dada`** — one per input sample.
- **`dada-pooled`** — exactly one, labelled `__pooled__`, because pooling
  denoises the merged table once. The per-sample output files are slices of that
  single run, not separate runs.
- **`dada-pseudo`** — one per sample for round 2, tagged `"round": 2`. Round 1
  is not collected.

`schema_version` bumps only when a field changes meaning or disappears; new
fields are added without a bump, so a consumer that ignores unknown keys keeps
working.

## Reading it

```bash
# where did a pooled run's time go?
jq '.runs[0].phases' run_metrics.json

# screen vs align, the question behind thread-count choice
jq '.runs[0].compare.split' run_metrics.json    # needs --metrics-attribution

# did the minimizer index get built, and was the probe right?
jq '.runs[0].index | {use_index, hindsight_disagrees}' run_metrics.json

# slowest sample in a per-sample run
jq -r '.runs | sort_by(-.phases.compare)[0] | "\(.sample) \(.phases.compare)s"' run_metrics.json
```

## Gotchas

- **`phases.loop_total` is not a supertotal.** The initial pre-loop compare is
  counted in `phases.compare` but runs before the loop timer starts, so
  `compare > loop_total` is normal.
- **`compare.attribution.unattributed`** is `total` minus the named parts. A
  large value means a phase is missing from the split, not that the work was
  free.
- **A non-empty `run.gates`** means `DADA2RS_*` tuning gates were active, so
  this run is not stock and its timings are not comparable with one that is.
- **Rates are `null`, not `0.0`, when their denominator is zero** — so "not
  measured" stays distinguishable from "measured as free".

## See also

- [`dada`](dada.md) · [`dada-pooled`](dada-pooled.md) · [`dada-pseudo`](dada-pseudo.md)
- [Inside b_compare — screen vs align](../findings/compare-screen-vs-align.md)
- [Thread scaling & memory placement](../findings/thread-scaling-and-placement.md)
