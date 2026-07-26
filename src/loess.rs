//! Locally-weighted polynomial regression (LOESS) smoother.
//!
//! A self-contained 1-D LOESS implementation mirroring R's `loess(y ~ x,
//! weights=w)` with `span = 0.75`, `degree = 2`.  Extracted from the error
//! model layer ([`crate::error_models`]) so the pure smoother — surface
//! selection, kd-tree partitioning, local weighted least-squares, and
//! flat extrapolation — lives apart from the DADA2 errfun family that
//! consumes it.
//!
//! Ports the relevant pieces of R's `stats` LOESS kernel
//! (`loessc.c`, `loessf.f`: `ehg124`, `ehg128`).

/// Default clamp bounds applied to off-diagonal rates returned by
/// [`crate::error_models::loess_errfun`] (and the loess-fallback paths in
/// `pacbio_errfun` / `binned_qual_errfun`).  Match R DADA2's `loessErrfun`
/// clamp at errorModels.R:53-56 (`MAX_ERROR_RATE <- 0.25;
/// MIN_ERROR_RATE <- 1e-7`, applied to the off-diagonal `est` matrix before
/// the diagonals are computed).  Both the `default` and `r-dada2` presets use
/// these.
pub const DEFAULT_MAX_ERROR_RATE: f64 = 0.25;
pub const DEFAULT_MIN_ERROR_RATE: f64 = 1e-7;

/// Gaussian elimination with partial pivoting.
///
/// Solves `A * x = b` where `a` is the row-major `n×n` matrix stored in a
/// flat slice and `b` is the right-hand-side vector.
/// Both slices are modified in place.  Returns `Some(x)` on success, `None`
/// if the system is (numerically) singular.
fn solve_linear(a: &mut [f64], b: &mut [f64], n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        // Find the row with the largest absolute value in this column (pivot).
        let pivot_row = (col..n).max_by(|&r1, &r2| {
            a[r1 * n + col]
                .abs()
                .partial_cmp(&a[r2 * n + col].abs())
                .unwrap()
        })?;

        if pivot_row != col {
            for k in 0..n {
                a.swap(col * n + k, pivot_row * n + k);
            }
            b.swap(col, pivot_row);
        }

        let pivot = a[col * n + col];
        if pivot.abs() < 1e-12 {
            return None;
        }

        for row in (col + 1)..n {
            let f = a[row * n + col] / pivot;
            a[row * n + col] = 0.0;
            for k in (col + 1)..n {
                let v = a[col * n + k] * f;
                a[row * n + k] -= v;
            }
            let bv = b[col] * f;
            b[row] -= bv;
        }
    }

    // Back substitution.
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i * n + j] * x[j];
        }
        let diag = a[i * n + i];
        if diag.abs() < 1e-12 {
            return None;
        }
        x[i] = s / diag;
    }
    Some(x)
}

/// Choice of fitting surface for [`loess_predict`].
///
/// Mirrors R's `loess(...)` `surface` parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LoessSurface {
    /// Evaluate the local polynomial directly at every query point.
    /// Bit-equivalent to R `loess(surface = "direct")`.
    Direct,
    /// Build a kd-tree partition of the data, fit local polynomials at each
    /// vertex, then blend neighboring vertex polynomials with cubic Hermite
    /// at query points.  Mirrors R `loess(surface = "interpolate")`, which
    /// is R DADA2's default `loessErrfun` path.
    ///
    /// `cell` controls the maximum number of observations per kd-tree cell;
    /// cells with more than `floor(cell * span * nv)` points are subdivided.
    /// R's default is `0.2` (see `loess.control`).
    Interpolate { cell: f64 },
}

/// Resolved configuration for the LOESS errfun and the loess-fallback path
/// inside [`crate::error_models::pacbio_errfun`] /
/// [`crate::error_models::binned_qual_errfun`].
///
/// The CLI layer builds this from a preset (`default` or `r-dada2`) plus any
/// explicit overrides.  Two presets are exposed:
///
/// | Knob | `default` | `r-dada2` |
/// |---|---|---|
/// | `surface` | `Direct` | `Interpolate { cell: 0.2 }` |
/// | `max_error_rate` | `0.25` | `0.25` |
/// | `min_error_rate` | `1e-7` | `1e-7` |
///
/// Both presets clamp off-diagonals to `[1e-7, 0.25]` — R DADA2's
/// `loessErrfun` (errorModels.R:53-56) applies the same clamp after the
/// `loess()` fit.  The presets differ only in the fitting surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LoessConfig {
    pub surface: LoessSurface,
    /// Upper clamp on off-diagonal rate output. `1.0` disables the upper clamp.
    pub max_error_rate: f64,
    /// Lower clamp on off-diagonal rate output. `0.0` disables the lower clamp.
    pub min_error_rate: f64,
}

impl Default for LoessConfig {
    /// The dada2-rs default: direct surface + R-style `[1e-7, 0.25]` clamp.
    fn default() -> Self {
        Self {
            surface: LoessSurface::Direct,
            max_error_rate: DEFAULT_MAX_ERROR_RATE,
            min_error_rate: DEFAULT_MIN_ERROR_RATE,
        }
    }
}

impl LoessConfig {
    /// `r-dada2` preset: surface=`Interpolate { cell: 0.2 }`, clamp at
    /// `[1e-7, 0.25]`.  Mirrors R DADA2's `loessErrfun` — both the R
    /// default `loess()` surface and the post-fit clamp applied at
    /// errorModels.R:53-56 (`MAX_ERROR_RATE <- 0.25; MIN_ERROR_RATE <- 1e-7`).
    pub fn r_dada2() -> Self {
        Self {
            surface: LoessSurface::Interpolate { cell: 0.2 },
            max_error_rate: DEFAULT_MAX_ERROR_RATE,
            min_error_rate: DEFAULT_MIN_ERROR_RATE,
        }
    }

    pub(crate) fn clamp(&self, r: f64) -> f64 {
        r.clamp(self.min_error_rate, self.max_error_rate)
    }
}

/// Fit a local weighted polynomial at `x0` and return the coefficients
/// in raw basis `[1, x, x², …]`.
///
/// Shared by [`loess_predict`]'s Direct (one call per query) and Interpolate
/// (one call per kd-tree vertex) paths.  Returns `None` only if the
/// neighborhood is empty after weighting, or if even a local constant fit is
/// numerically singular.
///
/// # Degree fallback
///
/// The tricube kernel assigns weight exactly zero to the farthest included
/// neighbor (`u == 1.0`), so a neighborhood of `n_local` points yields at most
/// `n_local - 1` positively-weighted observations.  When that is fewer than the
/// `p` coefficients a degree-`p-1` polynomial needs, we fit the highest degree
/// the surviving points *can* support instead of giving up.
///
/// This matters only in the degenerate regime.  With `span = 0.75` and
/// `degree = 2` (`p = 3`) the fallback engages at `nv <= 5`
/// (`floor(0.75 * 5) - 1 = 2 < 3`); from `nv >= 6` there are always at least
/// `p` surviving points and the fit is unchanged.  Previously this path
/// returned `None`, which propagated through [`extrapolate_flat`] to `NaN` and
/// then to `min_error_rate` for every cell — a silently floored error model on
/// binned-quality data, where 3–4 distinct Q values is normal.  See issue #95.
fn fit_local_at(
    x0: f64,
    valid: &[usize],
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    n_local: usize,
    p: usize,
) -> Option<Vec<f64>> {
    let mut dists: Vec<(usize, f64)> = valid.iter().map(|&i| (i, (xs[i] - x0).abs())).collect();
    dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let max_dist = dists[n_local - 1].1;

    let ws: Vec<(usize, f64)> = dists[..n_local]
        .iter()
        .map(|&(i, d)| {
            let u = if max_dist > 0.0 { d / max_dist } else { 0.0 };
            let tri = if u < 1.0 {
                (1.0 - u.powi(3)).powi(3)
            } else {
                0.0
            };
            (i, tri * weights[i])
        })
        .filter(|&(_, w)| w > 0.0)
        .collect();

    // Distinct x values bound the degree we can identify: fitting a degree-d
    // polynomial needs d+1 distinct abscissae regardless of how many
    // observations share them.
    let mut distinct_xs: Vec<f64> = ws.iter().map(|&(i, _)| xs[i]).collect();
    distinct_xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    distinct_xs.dedup();

    // Highest `p` the surviving neighborhood can support.  Equals `p` whenever
    // the neighborhood is non-degenerate, so the common path is unchanged.
    let mut p_eff = p.min(ws.len()).min(distinct_xs.len());

    while p_eff >= 1 {
        // Normal equations X^T W X a = X^T W y in raw basis [1, x, x², …].
        let mut xtx = vec![0.0f64; p_eff * p_eff];
        let mut xty = vec![0.0f64; p_eff];

        for &(i, w) in &ws {
            let xi = xs[i];
            let yi = ys[i];

            let mut row = vec![1.0f64; p_eff];
            let mut xpow = xi;
            for j in row.iter_mut().take(p_eff).skip(1) {
                *j = xpow;
                xpow *= xi;
            }

            for j in 0..p_eff {
                xty[j] += w * row[j] * yi;
                for k in 0..p_eff {
                    xtx[j * p_eff + k] += w * row[j] * row[k];
                }
            }
        }

        if let Some(coeffs) = solve_linear(&mut xtx, &mut xty, p_eff) {
            return Some(coeffs);
        }
        // Numerically rank-deficient at this degree — drop one and retry.
        p_eff -= 1;
    }

    None
}

/// Evaluate a polynomial in raw basis at `x` via Horner-equivalent loop.
#[inline]
fn eval_poly(coeffs: &[f64], x: f64) -> f64 {
    let mut acc = 0.0;
    let mut xpow = 1.0;
    for &c in coeffs {
        acc += c * xpow;
        xpow *= x;
    }
    acc
}

/// Build the unique sorted vertex positions of a 1-D kd-tree partition of
/// `sorted_valid_xs`.  Mirrors R's `ehg124` (`stats/src/loessf.f`):
/// cells with more than `threshold` observations are subdivided at the
/// **median data point** (rank `m = floor((l+u)/2)` within the cell),
/// and the new vertex is the median point's x value `x[pi(m)]`.
/// Left cell takes indices `l..=m`, right cell `m+1..=u`.
///
/// Returns at minimum `[x_min, x_max]`.  Vertices are sorted ascending and
/// deduplicated; consecutive pairs form the leaf cells.
fn build_kd_vertices_1d(sorted_valid_xs: &[f64], threshold: usize) -> Vec<f64> {
    let n = sorted_valid_xs.len();
    let x_min = sorted_valid_xs[0];
    let x_max = sorted_valid_xs[n - 1];
    let mut vertices = vec![x_min, x_max];

    // Stack of inclusive index ranges into `sorted_valid_xs`.
    let mut stack: Vec<(usize, usize)> = vec![(0, n - 1)];
    while let Some((l, u)) = stack.pop() {
        let count = u - l + 1;
        if count <= threshold {
            continue;
        }
        // R's `ehg124`: m = floor((l + u) / 2), vertex = x[pi(m)].
        // (l, u are 1-indexed in Fortran; here zero-indexed but the
        // arithmetic is identical.)
        let m = (l + u) / 2;
        if m == l || m == u {
            // No room to subdivide further while keeping both halves nonempty.
            continue;
        }
        let vertex_x = sorted_valid_xs[m];
        vertices.push(vertex_x);
        // Left: l..=m, right: m+1..=u (the median point belongs to the left).
        stack.push((l, m));
        stack.push((m + 1, u));
    }

    vertices.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vertices.dedup();
    vertices
}

/// Evaluate the value and first derivative of a polynomial in raw basis
/// `[c_0, c_1, c_2, …]` (so `f(x) = Σ c_j x^j`) at `x`.
#[inline]
fn eval_poly_and_deriv(coeffs: &[f64], x: f64) -> (f64, f64) {
    // value = Σ c_j x^j ;  derivative = Σ j·c_j x^(j-1)
    let mut val = 0.0;
    let mut xpow = 1.0;
    for &c in coeffs {
        val += c * xpow;
        xpow *= x;
    }
    let mut der = 0.0;
    let mut xpow_dm1 = 1.0; // x^(j-1) for j = 1
    for (j, &c) in coeffs.iter().enumerate().skip(1) {
        der += (j as f64) * c * xpow_dm1;
        xpow_dm1 *= x;
    }
    (val, der)
}

/// Locally-weighted polynomial regression (LOESS).
///
/// Mirrors R's `loess(y ~ x, data, weights=w)` with `span = 0.75` and
/// `degree = 2`.  Observations with non-finite `y` or zero weight are
/// excluded from fitting.  Predictions outside the valid data range are
/// returned as `None` so [`extrapolate_flat`] can fill them with the nearest
/// finite prediction.
///
/// `surface` selects between direct per-query fits and kd-tree-vertex fits
/// with smoothstep blending — see [`LoessSurface`].
///
/// Returns a `Vec<Option<f64>>` aligned to `xs`; `None` at a position means
/// the local fit could not be computed there.
pub(crate) fn loess_predict(
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    span: f64,
    degree: usize,
    surface: LoessSurface,
) -> Vec<Option<f64>> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert_eq!(xs.len(), weights.len());

    let n = xs.len();

    let valid: Vec<usize> = (0..n)
        .filter(|&i| ys[i].is_finite() && weights[i] > 0.0)
        .collect();
    let nv = valid.len();

    let eff_degree = degree.min(nv.saturating_sub(1));
    if nv <= eff_degree {
        return vec![None; n];
    }

    // R's `loess` (with surface="direct") uses `floor(span * n)` for the
    // neighborhood size; see `simpleLoess` → C kernel in `loessc.c`. dada2-rs
    // previously used `ceil`, which agrees when `span * nv` is integer but
    // differs by 1 otherwise — enough to nudge the local fit at nontrivial
    // numbers of observations. See issue #14 checklist item 1.
    let n_local = ((span * nv as f64).floor() as usize)
        .max(eff_degree + 1)
        .min(nv);
    let p = eff_degree + 1;

    let (x_min, x_max) = valid
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &i| {
            (lo.min(xs[i]), hi.max(xs[i]))
        });

    match surface {
        LoessSurface::Direct => (0..n)
            .map(|pred_idx| {
                let x0 = xs[pred_idx];
                if x0 < x_min || x0 > x_max {
                    return None;
                }
                let coeffs = fit_local_at(x0, &valid, xs, ys, weights, n_local, p)?;
                Some(eval_poly(&coeffs, x0))
            })
            .collect(),

        LoessSurface::Interpolate { cell } => {
            // Partition the valid x-range into kd-tree cells whose populations
            // are bounded by `floor(cell * span * nv)`.  Fit the local
            // polynomial at every unique vertex; queries are blended between
            // their enclosing cell's two vertex polynomials with a cubic
            // smoothstep.
            let threshold = (cell * span * nv as f64).floor() as usize;
            let threshold = threshold.max(p); // sanity floor

            let mut sorted_xs: Vec<f64> = valid.iter().map(|&i| xs[i]).collect();
            sorted_xs.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let vertices = build_kd_vertices_1d(&sorted_xs, threshold);

            // Fit at each vertex.  Store coefficients per vertex; `None` if
            // the local fit was rank-deficient.
            let vertex_coeffs: Vec<Option<Vec<f64>>> = vertices
                .iter()
                .map(|&v| fit_local_at(v, &valid, xs, ys, weights, n_local, p))
                .collect();

            (0..n)
                .map(|pred_idx| {
                    let x0 = xs[pred_idx];
                    if x0 < x_min || x0 > x_max {
                        return None;
                    }

                    // Find the cell [vertices[k], vertices[k+1]] containing x0.
                    // upper_bound: smallest k with vertices[k] > x0.  Handles
                    // x0 == x_min (k=1) and x0 == x_max (k=vertices.len()-1).
                    let upper = vertices.partition_point(|&v| v <= x0);
                    let (lo_idx, hi_idx) = if upper == 0 {
                        (0, 1.min(vertices.len() - 1))
                    } else if upper >= vertices.len() {
                        (vertices.len() - 2, vertices.len() - 1)
                    } else {
                        (upper - 1, upper)
                    };

                    let a = vertices[lo_idx];
                    let b = vertices[hi_idx];
                    let c_a = vertex_coeffs[lo_idx].as_ref()?;
                    let c_b = vertex_coeffs[hi_idx].as_ref()?;

                    let (f_a, fp_a) = eval_poly_and_deriv(c_a, a);
                    let (f_b, fp_b) = eval_poly_and_deriv(c_b, b);

                    if a == b {
                        return Some(f_a);
                    }
                    // Cubic Hermite interpolation between vertex-fits, per
                    // R `ehg128` (stats/src/loessf.f):
                    //   φ₀(h) = (1-h)²(1+2h)   value at left vertex
                    //   φ₁(h) = h²(3-2h)       value at right vertex
                    //   ψ₀(h) = h(1-h)²        derivative at left vertex
                    //   ψ₁(h) = h²(h-1)        derivative at right vertex
                    //   f(x) = φ₀·f(a) + φ₁·f(b) + (ψ₀·f'(a) + ψ₁·f'(b))·(b-a)
                    let h = (x0 - a) / (b - a);
                    let omh = 1.0 - h;
                    let phi_0 = omh * omh * (1.0 + 2.0 * h);
                    let phi_1 = h * h * (3.0 - 2.0 * h);
                    let psi_0 = h * omh * omh;
                    let psi_1 = h * h * (h - 1.0);
                    Some(phi_0 * f_a + phi_1 * f_b + (psi_0 * fp_a + psi_1 * fp_b) * (b - a))
                })
                .collect()
        }
    }
}

/// Apply flat extrapolation at the ends of a LOESS prediction vector.
///
/// Fills `None` entries at the low end with the first finite value and at the
/// high end with the last finite value.  Returns the filled vector.
pub(crate) fn extrapolate_flat(raw: Vec<Option<f64>>, n: usize) -> Vec<f64> {
    let valid: Vec<(usize, f64)> = raw
        .iter()
        .enumerate()
        .filter_map(|(i, p)| p.map(|v| (i, v)))
        .collect();

    if valid.is_empty() {
        return vec![f64::NAN; n];
    }

    let (min_i, min_v) = *valid.first().unwrap();
    let (max_i, max_v) = *valid.last().unwrap();

    let mut out = vec![f64::NAN; n];
    for &(i, v) in &valid {
        out[i] = v;
    }
    for i in out.iter_mut().take(min_i) {
        *i = min_v;
    }
    for i in out.iter_mut().take(n).skip(max_i + 1) {
        *i = max_v;
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Sparse quality anchors with the shape binned-quality data produces:
    /// a handful of distinct Q values, monotonically decreasing log10 rate,
    /// weighted by wildly unequal observation counts.
    fn sparse_case(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let xs: Vec<f64> = (0..n)
            .map(|i| 12.0 + 27.0 * i as f64 / (n - 1) as f64)
            .collect();
        let ys: Vec<f64> = (0..n)
            .map(|i| -1.5 - 3.0 * i as f64 / (n - 1) as f64)
            .collect();
        let ws: Vec<f64> = (0..n).map(|i| 100.0 * (i + 1) as f64).collect();
        (xs, ys, ws)
    }

    /// Regression test for issue #95: with `span = 0.75`, `degree = 2`, the
    /// tricube kernel zeroes the farthest included neighbor, leaving fewer
    /// than the 3 coefficients a quadratic needs whenever `nv <= 5`.  That used
    /// to return `None` at every query point, which `extrapolate_flat` turned
    /// into `NaN` and `loess_errfun` then floored at `min_error_rate` — a
    /// silently useless error model on binned-quality input.
    #[test]
    fn sparse_neighborhoods_still_fit() {
        for n in 3..=8 {
            let (xs, ys, ws) = sparse_case(n);
            for surface in [
                LoessSurface::Direct,
                LoessSurface::Interpolate { cell: 0.2 },
            ] {
                let pred = loess_predict(&xs, &ys, &ws, 0.75, 2, surface);
                assert!(
                    pred.iter().all(|p| p.is_some_and(|v| v.is_finite())),
                    "n={n} surface={surface:?} produced a non-finite prediction: {pred:?}"
                );
            }
        }
    }

    /// The fit at n <= 5 should not merely be finite — it should track the
    /// data.  These points are exactly collinear in log space, so any sane
    /// local fit reproduces them.
    #[test]
    fn sparse_fit_reproduces_collinear_data() {
        for n in 3..=6 {
            let (xs, ys, ws) = sparse_case(n);
            let pred = loess_predict(&xs, &ys, &ws, 0.75, 2, LoessSurface::Direct);
            for (i, p) in pred.iter().enumerate() {
                let v = p.expect("fit should succeed");
                assert!(
                    (v - ys[i]).abs() < 1e-9,
                    "n={n} i={i}: fit {v} != data {}",
                    ys[i]
                );
            }
        }
    }

    /// The degree fallback must engage *only* in the degenerate regime.  With
    /// enough observations the neighborhood supports a full quadratic, and a
    /// quadratic fit reproduces quadratic data exactly — which a linear
    /// fallback would not.
    #[test]
    fn dense_neighborhoods_keep_full_degree() {
        let xs: Vec<f64> = (0..40).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 - 0.3 * x + 0.01 * x * x).collect();
        let ws = vec![1.0f64; 40];

        let pred = loess_predict(&xs, &ys, &ws, 0.75, 2, LoessSurface::Direct);
        for (i, p) in pred.iter().enumerate() {
            let v = p.expect("fit should succeed");
            assert!(
                (v - ys[i]).abs() < 1e-8,
                "i={i}: quadratic fit {v} != {} — degree was silently reduced",
                ys[i]
            );
        }
    }

    /// Repeated x values do not add identifiability: three observations at two
    /// distinct abscissae can only support a line, and the solve must not be
    /// attempted at quadratic degree.
    #[test]
    fn repeated_abscissae_fall_back_to_supportable_degree() {
        let xs = vec![10.0, 10.0, 20.0, 20.0];
        let ys = vec![-2.0, -2.0, -4.0, -4.0];
        let ws = vec![50.0, 50.0, 80.0, 80.0];

        let pred = loess_predict(&xs, &ys, &ws, 0.75, 2, LoessSurface::Direct);
        for (i, p) in pred.iter().enumerate() {
            let v = p.expect("fit should succeed despite duplicate x");
            assert!(v.is_finite(), "i={i} produced {v}");
        }
    }

    /// A neighborhood with no positively-weighted observations still has no
    /// fit to give — `None` remains the answer, and `extrapolate_flat` turns a
    /// fully-empty prediction vector into NaN rather than a bogus value.
    #[test]
    fn no_valid_observations_still_yields_none() {
        let xs = vec![1.0, 2.0, 3.0];
        let ys = vec![f64::NAN, f64::NAN, f64::NAN];
        let ws = vec![1.0, 1.0, 1.0];

        let pred = loess_predict(&xs, &ys, &ws, 0.75, 2, LoessSurface::Direct);
        assert!(pred.iter().all(|p| p.is_none()));
        assert!(extrapolate_flat(pred, 3).iter().all(|v| v.is_nan()));
    }
}
