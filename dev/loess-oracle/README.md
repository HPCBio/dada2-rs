# loess-oracle

A comparison harness for the LOESS smoother behind `loess_errfun`. It fits one
set of real `(q, rlogp, tot)` triples with several implementations and reports
how far each lands from R's `stats::loess`, measured on the **clamped error
rate** that actually reaches inference rather than on raw fit residuals.

Built for issue [#84](https://github.com/HPCBio/dada2-rs/issues/84) (evaluating
the upstream [`loess-rs`](https://github.com/thisisamirv/loess-project) crate);
it is what turned up [#95](https://github.com/HPCBio/dada2-rs/issues/95). It is
**not** part of the dada2-rs build — own workspace, `publish = false`, and the
`loess-rs` dependency never reaches the shipped crate.

Distinct from `dev/compare_loess.R`, which compares whole `err_out` matrices
from a finished `learn-errors` run. This harness goes one level down, to the
individual weighted fits, and can synthesise sparse/binned inputs.

## Running

```sh
./run.sh <learn_errors.json> [outdir]     # everything, tables to stdout
```

Needs `cargo`, `python3`, and `Rscript` (base `stats` only — no dada2 needed).
The JSON must be a native dada2-rs `learn-errors` output, since it has to carry
a `trans` block.

Piecewise, if you want one case:

```sh
python3 extract_triples.py errs.json out/triples.tsv          # full grid
python3 extract_triples.py errs.json out/t.tsv --bins 12,24,38 # binned shape
python3 extract_triples.py errs.json out/t.tsv --nbins 5 --anchors-only
Rscript r_truth.R out/triples.tsv out/r_pred.tsv
cargo run -q -- out/triples.tsv out/rust_pred.tsv
python3 compare.py out/r_pred.tsv out/rust_pred.tsv --floor --by-q --signed
```

## Arms

| arm | what |
|---|---|
| `ours_direct`, `ours_interp` | dada2-rs `loess_predict` — the real code via path dep, not a copy |
| `lrs_direct`, `lrs_interp` | `loess-rs` 0.9.0 configured to match R (`iterations(0)`, no padding) |
| `lrs_robust_mad` | + 4 bisquare robustness iterations, MAD scale |
| `lrs_reflect`, `lrs_extend` | + non-flat boundary padding |

Ground truth is `loess(rlogp ~ q, weights = tot)` at `span = 0.75`,
`degree = 2`, both of R's surfaces.

## Reading the output — two traps

**`NA` is not the 1e-7 floor.** Our `loess_predict` evaluates at every grid
point and `loess_errfun` writes failures out as `min_error_rate`, so a floored
cell in an `ours_*` arm is real product behaviour. `loess-rs` 0.9.0 has **no
public predict-at-new-x** — `fit(&x, &y)` returns values only at the supplied
points — so it simply cannot speak to an unpopulated quality column. Those cells
are `NA` and skipped. Watch the `cells` column: on the binned shape the
`lrs_*` arms cover 192/480 because 28 of the 40 columns have no observations.
That coverage gap is itself a blocker for adopting the crate, since
`loess_errfun` must emit a rate for every column.

**R is not a usable oracle below ~6 populated columns.** With 3 anchors R fits
an exact quadratic and oscillates — on this data it returns rates that *rise*
with quality and hit the 0.25 ceiling — while emitting warnings that `dada()`
suppresses. A large `ours_direct` divergence in the binned case means R is
wrong, not us. See the decision recorded on #95.

## Findings (F3D0 MiSeqSOP, `loess-rs` 0.9.0, 2026-07)

Full 40-column grid, 28 populated, 12 transitions — `|log10(rate/rate_R)|`:

| arm | vs R | median | p95 | max | >2x |
|---|---|---|---|---|---|
| `ours_direct` | direct | 4.4e-14 | 1.0e-12 | 7.7e-12 | 0/480 |
| `lrs_direct` | direct | 2.3e-15 | 7.1e-15 | 2.5e-14 | 0/480 |
| `ours_interp` | interpolate | 1.0e-12 | 7.3e-04 | 1.2e-03 | 0/480 |
| `lrs_interp` | interpolate | 5.1e-15 | 2.3e-01 | 2.4e-01 | 0/480 |
| `lrs_robust_mad` | direct | 3.5e-02 | 1.8e-01 | 4.1e-01 | 8/480 |
| `lrs_reflect` | direct | 4.4e-02 | 2.1e-01 | 3.2e-01 | 13/480 |
| `lrs_extend` | direct | 2.3e-02 | 1.3e-01 | 2.5e-01 | 0/480 |

- **`custom_weights` is exactly R's `weights=`** — `w_ij = w_j x K(d_ij/h)`, and
  it reproduces R's weighted direct fit to machine epsilon, cleaner than our own
  port.
- **`lrs_interp` is not R's `ehg128`.** Machine-epsilon across the interior
  (q=15–36), diverging up to 0.24 log10 (~1.7x) at the boundary vertices
  (q=12–14, 37–39). Our port is ~200x closer there — and that is the surface the
  `r-dada2` preset uses.
- **Robustness and boundary padding move real amounts**, concentrated at
  q=31–37 where counts are highest, and robustness *raises* high-Q rates. Tracked
  in [#96](https://github.com/HPCBio/dada2-rs/issues/96).

Sparsity sweep (anchors only, `--nbins N --anchors-only`), max `|log10|` vs R:

| n_valid | `ours_direct` | `lrs_direct` | `ours_interp` | `lrs_interp` |
|---|---|---|---|---|
| 3 | 3.2e-15 | 1.5e-15 | 2.5e-01 | 1.5e-01 |
| 4 | 2.7e-15 | 4.0e-15 | 7.3e-01 | 2.8e-04 |
| 5 | 8.8e-15 | 4.9e-15 | 5.4e-01 | 2.8e-03 |
| 6 | 3.7e-12 | 2.4e-14 | 7.1e-01 | 1.9e-02 |
| 8 | 2.2e-12 | 2.3e-14 | 2.1e-01 | 4.5e-02 |

Three things to take from this:

- **Both direct implementations handle n=3.** An earlier write-up of this
  evaluation claimed `loess-rs` failed at n=3; that was an artefact of a
  `idx.len() < 4` guard in the throwaway probe this harness replaced — our code,
  not theirs. It matches R to ~4e-16 at n=3.
- **Our `ours_direct` only got here via #95.** Before that fix every one of
  these cells was pinned at the 1e-7 floor for n <= 5.
- **`ours_interp` is the weak arm at small n** (p95 up to 0.73 log10), where
  `lrs_interp` is much tighter. The kd-tree vertex partition has little to work
  with at 4–6 points. Not currently a product concern — nothing ships a
  sparse-input interpolate path — but it is the one place `loess-rs` is
  genuinely better on a surface we care about, and worth revisiting alongside
  binned-quality work.

Also worth noting for #96: `lrs_robust_mad` is unstable at small n — at n=7 it
reached 5.0 log10 divergence and floored 4/84 cells. Robustness iterations need
enough points to estimate a scale from.

## Status

Parked with #84. Re-run it if `loess-rs` cuts a 1.0, adds predict-at-new-x, or
fixes the interpolate boundary — and as a starting point for the #96 arms.
