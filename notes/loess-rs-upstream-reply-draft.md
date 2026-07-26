# Draft reply — thisisamirv/loess-rs#1 ("User-defined weights?")

**Status: NOT POSTED.** Draft only, pending review.

Context: `loess-rs` 0.9.0 (now under `thisisamirv/loess-project`) shipped
`custom_weights()` in response to our issue #1. Probe results that back the
numbers below are summarised in the loess-rs feasibility issue on this repo.

Note: the maintainer plans to retire the `loess-rs` repo, so the reply should
probably be posted on `loess-project` if #1 disappears first.

---

Thanks for shipping this — I've now tested `loess-rs` 0.9.0 against R's
`stats::loess` on real data (16S error-model fits, 12 transitions × 40 quality
bins, weighted by observation counts).

`custom_weights` is spot on: with `surface_mode("direct")`, `iterations(0)`,
`degree("quadratic")`, `fraction(0.75)`, your fitted values match R's
`loess(rlogp ~ q, weights = tot)` to **~2e-15 relative** — machine epsilon.
Exactly what I needed, thank you.

Two observations, one of which may be a bug:

1. **Interpolation surface at the boundaries.** `surface_mode("interpolation")`
   with `cell(0.2)` matches R's `surface="interpolate"` to machine epsilon
   across the interior, but diverges at the first and last ~3 x-values — up to
   0.24 in log10 terms. R's kernel (`ehg124`/`ehg128` in `loessf.f`) places
   kd-tree vertices at cell medians and blends with a cubic Hermite that uses
   the vertex-fit derivatives; the interior agreement suggests you're doing
   something very close, so this looks like boundary-vertex handling
   specifically. Happy to send the reproducing dataset.

2. **Degenerate neighborhoods at small n.** With `fraction(0.75)` and
   `degree("quadratic")` you produce good fits down to n=4 (R manages n=3). At
   n=3 the fit fails. Not a complaint — it's better than my own implementation,
   which fails at n<=5 — but worth documenting the supported floor, since
   binned-quality sequencing data legitimately gives you 3–4 distinct x-values.

On tests: my earlier offer stands and is now concrete. I have R `stats::loess`
reference fixtures with per-observation weights across both surfaces. If a
`tests/r_reference/` module with tabulated expected values would be useful, say
the word and I'll open a PR.

Also — issue #1 here is still open, and you mentioned retiring this repo; worth
closing it with a pointer to loess-project.
