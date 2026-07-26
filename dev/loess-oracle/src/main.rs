//! Fit one `(q, rlogp, tot)` triple set with several LOESS implementations and
//! emit their predictions on a common grid, for comparison against R's
//! `stats::loess` (see `r_truth.R`).
//!
//! Arms:
//!   - `ours_direct`  / `ours_interp`  — dada2-rs `loess_predict`, the real code
//!   - `lrs_direct`   / `lrs_interp`   — upstream `loess-rs` 0.9.0, R-matching config
//!   - `lrs_robust_mad`                — + bisquare robustness w/ MAD scale
//!   - `lrs_reflect`  / `lrs_extend`   — + non-flat boundary padding
//!
//! Usage: `loess-oracle <triples.tsv> <out.tsv>`
//!
//! Predictions are reported both as the raw log10 fit (`pred`) and as the error
//! rate that actually reaches inference (`rate` = `10^pred`, clamped to
//! `[1e-7, 0.25]` exactly as `loess_errfun` does). Compare on `rate`: two fits
//! can differ visibly in log space and be identical after clamping, and vice
//! versa.
//!
//! # `NA` vs the 1e-7 floor — do not conflate them
//!
//! `loess_predict` evaluates at *every* point of the quality grid, returning
//! `None` where the local fit fails; `loess_errfun` then writes those out as
//! `min_error_rate`. So for the `ours_*` arms a floored cell is real product
//! behaviour and is reported as `rate = 1e-7`.
//!
//! `loess-rs` 0.9.0 has **no public predict-at-new-x**: `fit(&x, &y)` returns
//! fitted values only at the supplied points (the `predict` in
//! `engine/executor.rs` is gated behind its `dev` feature). It therefore cannot
//! produce a value at an unpopulated quality column at all. Those cells are
//! written as `NA` and skipped by `compare.py` — reporting them as floored would
//! blame the crate for a fit failure that never happened, and would make it look
//! catastrophically worse than it is on binned input. The `cells` column in the
//! comparison output shows how much of the grid each arm could actually speak
//! to. This gap is itself a real blocker for adopting the crate; see README.

use dada2_rs::loess::{extrapolate_flat, loess_predict, LoessSurface};
use loess_rs::prelude::Loess;
use std::collections::BTreeMap;
use std::io::Write;

/// Matches `loess::DEFAULT_MIN_ERROR_RATE` / `DEFAULT_MAX_ERROR_RATE`.
const MIN_RATE: f64 = 1e-7;
const MAX_RATE: f64 = 0.25;

/// DADA2's `loessErrfun` span and degree — not configurable here on purpose.
/// The point of the harness is to compare *implementations* at the settings we
/// actually ship, not to sweep the parameters.
const SPAN: f64 = 0.75;
const DEGREE: usize = 2;
/// R's `loess.control(cell = 0.2)` default, used by the `r-dada2` preset.
const CELL: f64 = 0.2;

struct Series {
    q: Vec<f64>,
    tot: Vec<f64>,
    rlogp: Vec<f64>,
}

fn read_triples(path: &str) -> BTreeMap<String, Series> {
    let txt = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {path}: {e} (run extract_triples.py first)"));
    let mut out: BTreeMap<String, Series> = BTreeMap::new();
    for line in txt.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let e = out.entry(f[0].to_string()).or_insert(Series {
            q: vec![],
            tot: vec![],
            rlogp: vec![],
        });
        e.q.push(f[1].parse().expect("column 2 (q) must be numeric"));
        e.tot
            .push(f[3].parse().expect("column 4 (tot) must be numeric"));
        // Non-numeric / "nan" means tot == 0: no observations at this Q.
        e.rlogp.push(f[4].parse::<f64>().unwrap_or(f64::NAN));
    }
    out
}

fn ours_fit(s: &Series, surface: LoessSurface) -> Vec<f64> {
    let raw = loess_predict(&s.q, &s.rlogp, &s.tot, SPAN, DEGREE, surface);
    extrapolate_flat(raw, s.q.len())
}

/// Fit with `loess-rs` on the valid (`tot > 0`, finite `rlogp`) subset, then
/// scatter the fitted values back onto the full q grid and flat-fill the ends.
///
/// The subsetting is what makes this a fair comparison, not a convenience:
/// `loess_predict` excludes zero-weight observations *before* choosing the
/// neighborhood, and R's `loess()` drops the same rows via `na.omit`. Handing
/// `loess-rs` the full grid instead would let NaN rows consume neighborhood
/// slots and set the bandwidth, which neither of the other two implementations
/// does.
fn lrs_fit(s: &Series, cfg: impl Fn(Loess<f64>) -> Loess<f64>) -> Vec<f64> {
    let n = s.q.len();
    let idx: Vec<usize> = (0..n)
        .filter(|&i| s.rlogp[i].is_finite() && s.tot[i] > 0.0)
        .collect();
    if idx.is_empty() {
        return vec![f64::NAN; n];
    }
    let x: Vec<f64> = idx.iter().map(|&i| s.q[i]).collect();
    let y: Vec<f64> = idx.iter().map(|&i| s.rlogp[i]).collect();
    let w: Vec<f64> = idx.iter().map(|&i| s.tot[i]).collect();

    let builder = cfg(Loess::new()
        .fraction(SPAN)
        .degree("quadratic")
        .custom_weights(w));

    let fitted = match builder.build().and_then(|m| m.fit(&x, &y)) {
        Ok(r) => r.y,
        Err(e) => {
            // Expected at small n — loess-rs fails at n_valid <= 3. Report it
            // as NaN rather than aborting, so one degenerate transition does
            // not lose the whole run.
            eprintln!("  loess-rs declined a fit (n_valid={}): {e:?}", idx.len());
            return vec![f64::NAN; n];
        }
    };

    let mut raw: Vec<Option<f64>> = vec![None; n];
    for (k, &i) in idx.iter().enumerate() {
        raw[i] = Some(fitted[k]);
    }
    // Flat-fill outside the fitted range, mirroring R's minrli/maxrli and our
    // extrapolate_flat. Interior gaps stay NaN: the crate has no way to predict
    // there, which is a coverage limitation and not a failed fit.
    extrapolate_flat(raw, n)
}

/// `loess_errfun`'s back-transform: `10^pred`, clamped.
///
/// `floors_failures` distinguishes the two kinds of missing value. For the
/// `ours_*` arms a non-finite fit is what `loess_errfun` really turns into
/// `min_error_rate` (the silent-failure mode of issue #95), so report the floor.
/// For `loess-rs` a non-finite entry means the crate cannot predict at that
/// point at all, so report `NA` and let `compare.py` skip the cell.
fn format_cell(p: f64, floors_failures: bool) -> (String, String) {
    if p.is_finite() {
        let r = (10f64.powf(p)).clamp(MIN_RATE, MAX_RATE);
        (format!("{p}"), format!("{r}"))
    } else if floors_failures {
        (String::from("NA"), format!("{MIN_RATE:e}"))
    } else {
        (String::from("NA"), String::from("NA"))
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: loess-oracle <triples.tsv> <out.tsv>");
        std::process::exit(2);
    }
    let (inf, ouf) = (&args[1], &args[2]);

    let data = read_triples(inf);
    let mut out = std::fs::File::create(ouf).expect("cannot create output file");
    writeln!(out, "pair\tq\tarm\tpred\trate").unwrap();

    for (pair, s) in &data {
        // (name, predictions, does a missing value mean the shipped errfun
        // would floor it?)
        let arms: Vec<(&str, Vec<f64>, bool)> = vec![
            ("ours_direct", ours_fit(s, LoessSurface::Direct), true),
            (
                "ours_interp",
                ours_fit(s, LoessSurface::Interpolate { cell: CELL }),
                true,
            ),
            // R-matching config: no robustness iterations (R's family =
            // "gaussian") and no boundary padding.
            (
                "lrs_direct",
                lrs_fit(s, |b| {
                    b.surface_mode("direct")
                        .iterations(0)
                        .boundary_policy("noboundary")
                }),
                false,
            ),
            (
                "lrs_interp",
                lrs_fit(s, |b| {
                    b.surface_mode("interpolation")
                        .cell(CELL)
                        .iterations(0)
                        .boundary_policy("noboundary")
                }),
                false,
            ),
            // The features that motivated the evaluation (issue #96).
            (
                "lrs_robust_mad",
                lrs_fit(s, |b| {
                    b.surface_mode("direct")
                        .iterations(4)
                        .robustness_method("bisquare")
                        .scaling_method("mad")
                        .boundary_policy("noboundary")
                }),
                false,
            ),
            (
                "lrs_reflect",
                lrs_fit(s, |b| {
                    b.surface_mode("direct")
                        .iterations(0)
                        .boundary_policy("reflect")
                }),
                false,
            ),
            (
                "lrs_extend",
                lrs_fit(s, |b| {
                    b.surface_mode("direct")
                        .iterations(0)
                        .boundary_policy("extend")
                }),
                false,
            ),
        ];
        for (arm, pred, floors_failures) in arms {
            for (k, &p) in pred.iter().enumerate() {
                let (pred_s, rate_s) = format_cell(p, floors_failures);
                writeln!(out, "{}\t{}\t{}\t{}\t{}", pair, s.q[k], arm, pred_s, rate_s).unwrap();
            }
        }
    }
    eprintln!("wrote {ouf}");
}
