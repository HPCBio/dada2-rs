//! Abundance p-value and error-model lambda calculations.
// Function names like calc_pA / get_pA intentionally match the C++ source.
#![allow(non_snake_case)]
//!
//! Ported from `pval.cpp`.
//!
//! The Poisson CDF previously called via `Rcpp::ppois` is replaced by
//! `statrs::distribution::Poisson`, which uses the same regularised
//! incomplete gamma function that R's `ppois` calls internally.

use statrs::distribution::{DiscreteCDF, Poisson};

use crate::containers::{B, Raw, Sub};
use std::sync::OnceLock;

/// Minimum value of the conditioning normaliser below which the second-order
/// Taylor approximation `E - E²/2` is used instead of `1 - exp(-E)`.
/// Matches C++ `TAIL_APPROX_CUTOFF` in `dada.h`.
const TAIL_APPROX_CUTOFF: f64 = 1e-7;

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

/// Attribution of `b_p_update`'s cost, gathered by counting rather than timing
/// (issue #154).
///
/// The phase is 10-12% of `run_dada` on soil ITS2 and has never been optimised,
/// but its headline figure — ~56 ns per repricing — is an average over four
/// paths of very different cost. `get_pA` returns early and free for
/// singletons, for cluster centres, and for zero lambda; only the remainder
/// reaches [`calc_pA`], which evaluates a regularised incomplete gamma. Whether
/// the phase is worth parallelising, and what the payoff would be, depends
/// entirely on the mix — an average over a cheap majority and an expensive
/// minority is not a cost model.
///
/// **Counted, not timed, on purpose.** Two `Instant::now()` calls cost ~40 ns
/// against a ~56 ns budget, so per-item timing here would measure the
/// instrument rather than the phase. Counters are near-free; the per-path
/// prices come from `examples/pval_cost.rs` instead.
#[derive(Clone, Copy, Default, Debug)]
pub struct PUpdateStats {
    /// Raws whose `p` was recomputed (members of every cluster with `update_e`).
    pub repriced: u64,
    /// Returned `1.0` immediately: singleton, not a prior, singleton detection off.
    pub exit_singleton: u64,
    /// Returned `1.0` immediately: `hamming == 0`, i.e. a cluster centre or exact match.
    pub exit_center: u64,
    /// Returned `0.0` immediately: `lambda == 0`, outside the k-mer screen.
    pub exit_zero_lambda: u64,
    /// Reached `calc_pA` — the Poisson upper tail. The expensive path.
    pub full_calc: u64,
    /// Raws walked by the greedy lock pass, which is a second traversal of the
    /// same members and is not counted in `repriced`.
    pub lock_scanned: u64,
    /// Clusters found dirty (`update_e` set) and repriced.
    pub dirty_clusters: u64,
    /// Clusters skipped because `update_e` was clear.
    pub clean_clusters: u64,
}

impl std::ops::AddAssign for PUpdateStats {
    fn add_assign(&mut self, o: Self) {
        self.repriced += o.repriced;
        self.exit_singleton += o.exit_singleton;
        self.exit_center += o.exit_center;
        self.exit_zero_lambda += o.exit_zero_lambda;
        self.full_calc += o.full_calc;
        self.lock_scanned += o.lock_scanned;
        self.dirty_clusters += o.dirty_clusters;
        self.clean_clusters += o.clean_clusters;
    }
}

/// Distance, in loop iterations, at which `b_p_update` prefetches the `Raw` it
/// will read next (issue #154). Default `16`; `0` disables. Overridable via
/// `DADA2RS_PUPDATE_PREFETCH` so both arms of an A/B come from one binary.
///
/// 16 is where the benefit saturates: on soil ITS2 the sweep reads 23.6 ns per
/// repricing at 4, 19.2 at 8, then 17.0 at both 16 and 32 against a 27.7 ns
/// baseline. 16 buys the same as 32 with half the outstanding prefetches.
pub const PUPDATE_PREFETCH_DEFAULT: usize = 16;

pub fn prefetch_distance() -> usize {
    static VALUE: OnceLock<usize> = OnceLock::new();
    *VALUE.get_or_init(|| {
        std::env::var("DADA2RS_PUPDATE_PREFETCH")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(PUPDATE_PREFETCH_DEFAULT)
    })
}

/// Ask the memory system for the two fields `b_p_update` is about to read from
/// `raw`, `prefetch_distance()` iterations before it reads them.
///
/// **Why this and not a layout change.** The phase is memory-bound on a random
/// gather: over 91% of repricings take an early exit doing no arithmetic, and
/// the cost is a cache miss on a 160-byte `Raw` read for ~20 bytes of it.
/// Three layout fixes were priced (`examples/pval_layout.rs`) and all are worse
/// than this one:
///
/// | candidate | soil ITS2 | soil 16S |
/// |---|---|---|
/// | today | 27.7 ns | 29.0 ns |
/// | `repr(C, align(64))` | 37.0 (worse) | 40.4 (worse) |
/// | `p`/`comp` in dense arrays | 41.7 (worse) | 44.8 (worse) |
/// | every hot field packed to 32 B | **11.7** | 33.3 (worse) |
/// | **prefetch, d=16** | **17.0** | **20.2** |
///
/// The packed layout is the only competitive one and it wins by fitting a 32 MB
/// per-CCD L3 (26.4 MB on ITS2) — so it *loses* on 16S (39.2 MB), needs `reads`
/// and `prior` moved out of `Raw`, and its advantage swung 37% between runs
/// because an array sitting on the cache-capacity boundary depends on having
/// that cache to itself. Prefetching wins on both pools, needs no refactor,
/// varied 0% at the plateau, and does not decay as pools grow.
///
/// Addresses come from `addr_of!` on the real fields rather than hardcoded
/// offsets: `Raw` is `repr(Rust)`, so its layout is the compiler's to choose and
/// a hardcoded offset would silently prefetch the wrong line after any field
/// reordering.
///
/// **Byte-identical by construction**: a prefetch has no architectural effect.
#[inline(always)]
fn prefetch_raw(raw: &Raw) {
    #[cfg(target_arch = "x86_64")]
    {
        use std::arch::x86_64::{_MM_HINT_T0, _mm_prefetch};
        // SAFETY: both pointers are derived from a live `&Raw`, so they are
        // valid for reads. `_mm_prefetch` is a hint with no architectural
        // effect and cannot fault regardless.
        unsafe {
            _mm_prefetch(
                std::ptr::addr_of!(raw.comp.lambda) as *const i8,
                _MM_HINT_T0,
            );
            _mm_prefetch(std::ptr::addr_of!(raw.reads) as *const i8, _MM_HINT_T0);
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = raw;
}

/// Update abundance p-values for every Raw in the partition.
///
/// For each cluster whose `update_e` flag is set, recomputes `raw.p` for all
/// member Raws and clears the flag. When `greedy` is true, also locks Raws
/// whose expected read count from the cluster center alone already exceeds
/// their observed count (preventing them from budding a new cluster).
///
/// Returns the number of Raws whose `p` was recomputed this call — i.e. the
/// members of every cluster with `update_e` set. This is the per-round p-churn
/// that gated the b_bud incremental follow-up (issue #85).
///
/// While repricing each dirty cluster's members, this also refreshes that
/// cluster's cached best budding candidate (`Bi::bud_min` / `bud_min_prior`),
/// applying the same eligibility filters and tie-break as `b_bud`'s scan. The
/// member loop already touches exactly these Raws with `bi_reads` in hand, so
/// the cache is maintained at no extra memory traffic — and `b_bud` then
/// combines the per-cluster minima in O(nclusters) instead of rescanning every
/// Raw (issue #85). `min_fold`/`min_hamming`/`min_abund` are `b_bud`'s
/// candidate filters, threaded through so the cache matches its scan exactly.
///
/// Equivalent to C++ `b_p_update`.
pub fn b_p_update(
    b: &mut B,
    greedy: bool,
    detect_singletons: bool,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
) -> PUpdateStats {
    use crate::containers::BudCand;
    let mut st = PUpdateStats::default();
    for i in 0..b.clusters.len() {
        if !b.clusters[i].update_e {
            st.clean_clusters += 1;
        }
        if b.clusters[i].update_e {
            st.dirty_clusters += 1;
            // Clone indices to avoid holding a shared borrow on b.clusters
            // while mutating b.raws.
            let members: Vec<usize> = b.clusters[i].raws.clone();
            let bi_reads = b.clusters[i].reads;
            // Refresh this cluster's cached bud candidate in the same pass.
            // Position 0 is skipped, mirroring b_bud's `for r in 1..len`.
            let mut bud_min: Option<BudCand> = None;
            let mut bud_min_prior: Option<BudCand> = None;
            let pf = prefetch_distance();
            for (r, &raw_idx) in members.iter().enumerate() {
                // The member list is fully known before the loop runs, so the
                // address needed `pf` iterations from now is available now.
                if pf > 0
                    && let Some(&ahead) = members.get(r + pf)
                {
                    prefetch_raw(&b.raws[ahead]);
                }
                let p = get_pA_counted(
                    b.raws[raw_idx].reads,
                    b.raws[raw_idx].prior,
                    b.raws[raw_idx].comp.lambda,
                    b.raws[raw_idx].comp.hamming,
                    bi_reads,
                    detect_singletons,
                    Some(&mut st),
                );
                b.raws[raw_idx].p = p;
                st.repriced += 1;

                if r == 0 {
                    continue; // center slot: never a bud candidate
                }
                let raw = &b.raws[raw_idx];
                if raw.reads < min_abund || raw.comp.hamming < min_hamming {
                    continue;
                }
                let fold_ok = min_fold <= 1.0
                    || raw.reads as f64 >= min_fold * raw.comp.lambda * bi_reads as f64;
                if !fold_ok {
                    continue;
                }
                let cand = BudCand {
                    p,
                    reads: raw.reads,
                    r,
                    raw_idx,
                };
                // Ascending r with strict replacement keeps the first-in-scan
                // (min p, then max reads, then lowest r) — identical to b_bud.
                if bud_min.is_none_or(|m| p < m.p || (p == m.p && cand.reads > m.reads)) {
                    bud_min = Some(cand);
                }
                if raw.prior
                    && bud_min_prior.is_none_or(|m| p < m.p || (p == m.p && cand.reads > m.reads))
                {
                    bud_min_prior = Some(cand);
                }
            }
            b.clusters[i].bud_min = bud_min;
            b.clusters[i].bud_min_prior = bud_min_prior;
            b.clusters[i].update_e = false;
        }

        if greedy && b.clusters[i].check_locks {
            let center_idx = b.clusters[i].center;
            let center_reads = center_idx.map_or(0, |ci| b.raws[ci].reads);
            let members: Vec<usize> = b.clusters[i].raws.clone();
            for raw_idx in members {
                st.lock_scanned += 1;
                // Lock if the center alone expects more reads than observed.
                let e_from_center = center_reads as f64 * b.raws[raw_idx].comp.lambda;
                if e_from_center > b.raws[raw_idx].reads as f64 {
                    b.raws[raw_idx].lock = true;
                }
                // Always lock the center itself.
                if Some(raw_idx) == center_idx {
                    b.raws[raw_idx].lock = true;
                }
            }
            b.clusters[i].check_locks = false;
        }
    }
    st
}

/// At loop termination, cross-check the incremental `bud_min` cache against a
/// full serial scan. If the scan finds a candidate the cache did not offer,
/// the loop stopped early and every raw behind that candidate is unreachable.
///
/// `b_bud_incremental` already asserts this equivalence -- but only under
/// `#[cfg(debug_assertions)]`, so it has never run on a production workload.
/// A 547k-unique pooled PacBio run cannot be done in a debug build, which is
/// precisely where a stale cache would matter (issue #219: a 66,937-read
/// sequence at hamming 311 with an underflowed lambda, whose `get_pA` value is
/// 0.0, never appears among 2,818 divisions).
///
/// O(raws) once at termination -- free there, prohibitive per round.
/// Non-fatal: it reports, it does not abort a long run.
pub fn report_missed_bud(
    b: &B,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
    verbose: bool,
) -> Option<(usize, f64, u32)> {
    let (mini, min_p, min_reads, _, _, _) = b_bud_scan_select(b, min_fold, min_hamming, min_abund);
    let nraw = b.raws.len() as f64;
    let p_a = min_p * nraw;
    if p_a >= b.omega_a {
        return None; // the scan agrees: nothing left to bud
    }
    let (ci, _, raw_idx) = mini?;
    if verbose {
        let r = &b.raws[raw_idx];
        eprintln!(
            "[dada] WARNING: bud loop stopped while a full scan still finds a candidate \
             (issue #219): raw {raw_idx} in cluster {ci}, {} reads, hamming {}, lambda {:.3e}, \
             p {:.3e}, pA {:.3e} < omega_a {:.3e}",
            r.reads, r.comp.hamming, r.comp.lambda, min_p, p_a, b.omega_a
        );
    }
    Some((raw_idx, min_p, min_reads))
}

/// Non-mutating serial equivalent of `b_bud`'s candidate selection: scans every
/// non-position-0 raw across all clusters and returns the abundance and prior
/// minima it would pick, as `(mini, min_p, min_reads, mini_prior, min_p_prior,
/// min_reads_prior)` where each `mini*` is `Option<(ci, r, raw_idx)>`. Used only
/// by `b_bud`'s debug cross-check against the incremental cache (issue #85).
#[allow(clippy::type_complexity)]
pub fn b_bud_scan_select(
    b: &B,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
) -> (
    Option<(usize, usize, usize)>,
    f64,
    u32,
    Option<(usize, usize, usize)>,
    f64,
    u32,
) {
    let init_center = b.clusters[0]
        .center
        .expect("b_bud_scan_select: cluster 0 has no center");
    let mut mini: Option<(usize, usize, usize)> = None;
    let mut mini_prior: Option<(usize, usize, usize)> = None;
    let mut min_p = b.raws[init_center].p;
    let mut min_p_prior = b.raws[init_center].p;
    let mut min_reads = b.raws[init_center].reads;
    let mut min_reads_prior = b.raws[init_center].reads;

    for ci in 0..b.clusters.len() {
        for r in 1..b.clusters[ci].raws.len() {
            let raw_idx = b.clusters[ci].raws[r];
            let raw = &b.raws[raw_idx];
            if raw.reads < min_abund || raw.comp.hamming < min_hamming {
                continue;
            }
            let fold_ok = min_fold <= 1.0
                || raw.reads as f64 >= min_fold * raw.comp.lambda * b.clusters[ci].reads as f64;
            if !fold_ok {
                continue;
            }
            if raw.p < min_p || (raw.p == min_p && raw.reads > min_reads) {
                mini = Some((ci, r, raw_idx));
                min_p = raw.p;
                min_reads = raw.reads;
            }
            if raw.prior
                && (raw.p < min_p_prior || (raw.p == min_p_prior && raw.reads > min_reads_prior))
            {
                mini_prior = Some((ci, r, raw_idx));
                min_p_prior = raw.p;
                min_reads_prior = raw.reads;
            }
        }
    }
    (
        mini,
        min_p,
        min_reads,
        mini_prior,
        min_p_prior,
        min_reads_prior,
    )
}

/// Abundance p-value: P(X ≥ `reads` | Poisson(λ = `e_reads`)).
///
/// When `prior` is false the p-value is conditioned on the sequence being
/// present (i.e. at least one read), normalising by `1 - exp(-e_reads)`.
/// A second-order Taylor expansion replaces the normaliser when it would
/// underflow below `TAIL_APPROX_CUTOFF`.
///
/// # Panics
/// Panics in debug builds if `reads == 0` (callers must guard this case).
/// Equivalent to C++ `calc_pA`.
pub fn calc_pA(reads: u32, e_reads: f64, prior: bool) -> f64 {
    debug_assert!(reads > 0, "calc_pA: reads must be > 0");

    // R uses ppois(reads-1, e_reads, lower.tail=FALSE), which under the hood
    // calls the regularised lower incomplete gamma function. statrs's `sf`
    // does the same nominally — sf(x) = gamma_lr(x+1, λ) — but its
    // gamma_lr implementation underflows to 0 for very small λ, while R's
    // pgamma stays accurate down to ~1e-300. For tiny λ we fall back to the
    // direct upper-tail Poisson series in log space:
    //   P(X >= k | λ) = e^{-λ} λ^k/k! * Σ_{j>=0} λ^j / ∏_{i=1..j}(k+i)
    // which is dominated by the leading e^{-λ} λ^k/k! for small λ.
    // Zero expected reads. R evaluates `ppois(reads-1, 0, lower.tail=FALSE)`,
    // which is P(X > reads-1 | lambda = 0) = 0 for every reads >= 1 -- verified
    // against R directly. We returned 1.0 here, the opposite end of the range,
    // and that inversion is not cosmetic: the post-loop omega_c pass
    // (`dada.rs`) sets `correct = p >= omega_c`, so 1.0 CORRECTS the raw into
    // its cluster where R's 0.0 leaves it unplaced.
    //
    // On the pinned 95-sample pooled PacBio run that folded a 66,937-read
    // sequence -- hamming 311 from the centre, lambda underflowed to exactly
    // zero -- into cluster 0, while R left those reads out entirely. Same
    // mechanism as issue #204, where we place every read and R leaves
    // ~1.6-1.9% unplaced.
    //
    // Only `e_reads == 0` is reachable: lambda and reads are both
    // non-negative. A non-finite e_reads would be a bug upstream, and R would
    // return NaN there rather than a probability, so it is not special-cased
    // into a silent answer.
    if e_reads <= 0.0 {
        return 0.0;
    }
    let pois = match Poisson::new(e_reads) {
        Ok(p) => p,
        Err(_) => return 0.0, // unreachable given the guard above
    };
    let mut pval = pois.sf((reads - 1) as u64);
    if pval == 0.0 && e_reads > 0.0 {
        pval = poisson_upper_tail_direct(reads, e_reads);
    }

    if prior {
        return pval;
    }

    // Condition on the sequence being present (at least one read observed).
    let norm = 1.0 - (-e_reads).exp();
    let norm = if norm < TAIL_APPROX_CUTOFF {
        // 2nd-order Taylor: 1 - e^{-E} ≈ E - E²/2 for small E.
        e_reads - 0.5 * e_reads * e_reads
    } else {
        norm
    };
    pval / norm
}

/// Compute lambda: the probability under the error model that `raw`'s
/// sequence arose from the reference sequence described by `sub`.
///
/// `err_mat` is a flat row-major matrix of shape 16 × `ncol` where rows
/// index the 16 nucleotide transitions (ref_nt × 4 + query_nt, with A=0,
/// C=1, G=2, T=3) and columns index the rounded quality score.
///
/// Returns `0.0` when `sub` is `None` (sequence was outside the k-mer
/// distance threshold and was not aligned).
/// Equivalent to C++ `compute_lambda_ts`.
pub fn compute_lambda(
    raw: &Raw,
    sub: Option<&Sub>,
    err_mat: &[f64],
    ncol: usize,
    use_quals: bool,
) -> f64 {
    let sub = match sub {
        Some(s) => s,
        None => return 0.0,
    };

    let len = raw.len();

    // Initialise per-position transition and quality-index vectors.
    // Default: diagonal match (nti * 4 + nti for both ref and query).
    let mut tvec = vec![0usize; len];
    let mut qind = vec![0usize; len];

    for pos in 0..len {
        let nti = raw.seq[pos]
            .checked_sub(1)
            .filter(|&n| n < 4)
            // Unreachable: dada_uniques_cached rejects non-ACGT input up front
            // (issue #101), which is where a user-facing error belongs — this is
            // inside a rayon worker and knows neither the sample nor the sequence.
            // Kept as a hard invariant rather than silently scoring an unknown
            // base as a match, which would corrupt lambda instead of failing.
            .expect(
                "non-ACGT nucleotide reached compute_lambda; input validation in \
                 dada_uniques_cached should have rejected it (bug)",
            ) as usize;
        tvec[pos] = nti * 4 + nti;
        qind[pos] = if use_quals {
            raw.qual.as_ref().map_or(0, |q| q[pos] as usize)
        } else {
            0
        };
    }

    // Override transition at each substitution position.
    for s in 0..sub.nsubs() {
        let pos0 = sub.pos[s] as usize;
        let pos1 = sub.map[pos0] as usize;

        debug_assert!(
            pos0 < sub.len0 as usize,
            "sub pos0 {pos0} >= len0 {}",
            sub.len0
        );
        debug_assert!(pos1 < len, "sub pos1 {pos1} >= raw len {len}");

        let nti0 = sub.nt0[s].saturating_sub(1) as usize;
        let nti1 = sub.nt1[s].saturating_sub(1) as usize;
        tvec[pos1] = nti0 * 4 + nti1;
    }

    // Lambda = product of error probabilities across all positions.
    let lambda: f64 = (0..len)
        .map(|pos| err_mat[tvec[pos] * ncol + qind[pos]])
        .product();

    debug_assert!(
        (0.0..=1.0).contains(&lambda),
        "lambda {lambda} outside [0,1]"
    );
    lambda
}

/// Direct upper-tail Poisson computation in log space, robust to extreme
/// underflow.
///
/// Used as a fallback when `statrs::Poisson::sf` returns exactly 0 due to
/// `gamma_lr` losing precision for very small λ. The leading term
/// `e^{-λ} λ^reads / reads!` dominates for small λ; the series correction
/// `1 + λ/(k+1) + λ²/((k+1)(k+2)) + ...` converges quickly.
fn poisson_upper_tail_direct(reads: u32, lambda: f64) -> f64 {
    debug_assert!(reads > 0);
    debug_assert!(lambda > 0.0);

    let k = reads as f64;
    let log_lambda = lambda.ln();
    let log_kfact: f64 = (1..=reads).map(|i| (i as f64).ln()).sum();
    let log_leading = -lambda + k * log_lambda - log_kfact;

    let mut series_sum: f64 = 1.0;
    let mut term: f64 = 1.0;
    for j in 1..10000u32 {
        term *= lambda / (k + j as f64);
        if !term.is_finite() || term <= 0.0 {
            break;
        }
        let new_sum = series_sum + term;
        if new_sum == series_sum {
            break;
        }
        series_sum = new_sum;
    }

    (log_leading + series_sum.ln()).exp()
}

/// Self-production probability: the probability that a sequence is produced
/// from itself given a compact 4×4 match-probability matrix.
///
/// Returns the product of `err[nti][nti]` (match diagonal) over all positions.
/// Equivalent to C++ `get_self`.
#[allow(dead_code)]
pub fn get_self(seq: &[u8], err: &[[f64; 4]; 4]) -> f64 {
    seq.iter().fold(1.0, |acc, &nt| {
        let nti = (nt as usize).saturating_sub(1).min(3);
        acc * err[nti][nti]
    })
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Per-raw p-value calculation, handling all sentinel cases before calling
/// `calc_pA`. Factored out to allow use from `b_p_update` without holding
/// simultaneous borrows on both `B.raws` and `B.clusters`.
/// Equivalent to C++ `get_pA`.
/// Abundance p-value for one Raw, with optional path attribution (#154).
///
/// Was `get_pA`; the un-counted wrapper had no callers once `b_p_update`
/// started attributing, so it folded away rather than sitting as dead code.
///
/// `stats` is `None` on any path that does not care, so the counting folds
/// away; `b_p_update` passes `Some` because the mix of early exits versus full
/// evaluations is the thing that decides whether this phase is worth
/// parallelising. Equivalent to C++ `get_pA`.
fn get_pA_counted(
    reads: u32,
    prior: bool,
    lambda: f64,
    hamming: u32,
    bi_reads: u32,
    detect_singletons: bool,
    mut stats: Option<&mut PUpdateStats>,
) -> f64 {
    macro_rules! bump {
        ($field:ident) => {
            if let Some(s) = stats.as_deref_mut() {
                s.$field += 1;
            }
        };
    }
    if reads == 1 && !prior && !detect_singletons {
        // Singleton: no abundance p-value is applied.
        bump!(exit_singleton);
        return 1.0;
    }
    if hamming == 0 {
        // Cluster center (or exact match): always valid.
        bump!(exit_center);
        return 1.0;
    }
    if lambda == 0.0 {
        // Zero expected reads: reject unconditionally.
        bump!(exit_zero_lambda);
        return 0.0;
    }
    bump!(full_calc);
    let e_reads = lambda * bi_reads as f64;
    calc_pA(reads, e_reads, prior || detect_singletons)
}

#[cfg(test)]
mod tests {
    use crate::containers::{B, Comparison, Raw};

    /// Build a two-raw pool: an abundant centre and an abundant, very divergent
    /// member. Mirrors the shape that stranded a 66,937-read PacBio sequence
    /// (issue #219).
    fn pool_with_divergent_member() -> B {
        let a = b"ACGTACGTAGCTAGCTAAGGCCTTAGCTAGCTACGTACGTTTGACTGACAGCTTAAGGCCA".to_vec();
        let mut c = a.clone();
        for i in (0..c.len()).step_by(3) {
            c[i] = if c[i] == b'A' { b'T' } else { b'A' };
        }
        let raws = vec![
            Raw::new(a, None, 200_000, false),
            Raw::new(c, None, 66_937, false),
        ];
        let mut b = B::new(raws, 1e-40, 1e-4, false);
        b.clusters[0].center = Some(0);
        b
    }

    /// A non-singleton whose lambda underflows to exactly zero must get
    /// `p = 0.0` -- maximally significant, so `b_bud` can promote it. Returning
    /// 1.0 makes it permanently unbuddable, which is what stranded a
    /// 66,937-read / hamming-311 sequence inside cluster 0 on the 95-sample
    /// PacBio run while R called it as its own ASV with 106,853 reads.
    ///
    /// The bug is invisible on singletons: the singleton branch returns 1.0
    /// anyway, and 1,741 of the 1,742 raws in that class on the real run were
    /// singletons.
    #[test]
    fn zero_lambda_non_singleton_is_maximally_significant() {
        let mut b = pool_with_divergent_member();
        // The state the trace recorded: real alignment, underflowed lambda.
        b.raws[1].comp = Comparison {
            i: 0,
            index: 1,
            lambda: 0.0,
            hamming: 311,
        };
        b.clusters[0].update_e = true;
        super::b_p_update(&mut b, false, false, 1.0, 1, 1);
        assert_eq!(
            b.raws[1].p, 0.0,
            "a zero-lambda non-singleton must be maximally significant, not p=1.0"
        );
    }
}

#[cfg(test)]
mod diag_tests {
    /// Zero expected reads must be MAXIMALLY significant, matching R's
    /// `ppois(reads-1, 0, lower.tail=FALSE) == 0` (checked against R itself).
    /// We returned 1.0, which flipped the post-loop omega_c verdict from
    /// "leave unplaced" to "correct into this cluster" and folded a
    /// 66,937-read PacBio sequence into a cluster 311 mismatches away
    /// (#219, and the mechanism behind #204).
    #[test]
    fn calc_pa_with_zero_e_reads_is_zero_like_r() {
        for reads in [1u32, 2, 66_937] {
            for prior in [true, false] {
                assert_eq!(
                    super::calc_pA(reads, 0.0, prior),
                    0.0,
                    "reads={reads} prior={prior}"
                );
            }
        }
    }

    /// The guard must not disturb ordinary values.
    #[test]
    fn calc_pa_unchanged_for_positive_e_reads() {
        let p = super::calc_pA(5, 2.0, true);
        assert!(p > 0.0 && p < 1.0, "expected a real probability, got {p}");
    }

    use crate::containers::{B, Comparison, Raw};

    /// The diagnostic must actually fire on the state it was written for --
    /// a null from an instrument is worth nothing until the instrument is
    /// shown capable of returning something else. This reproduces the #219
    /// raw: 66,937 reads, hamming 311, lambda underflowed to 0, cached p 1.0.
    /// The cross-check must fire when the cache would miss a candidate. Built
    /// by hand: a raw whose `p` makes it a valid bud target, with the
    /// per-cluster `bud_min` cache left empty as a stale cache would leave it.
    #[test]
    fn report_missed_bud_fires_when_the_cache_is_stale() {
        let a = b"ACGTACGTAGCTAGCTAAGGCCTTAGCTAGCTACGTACGTTTGACTGACAGCTTAAGGCCA".to_vec();
        let z = b"TTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCC".to_vec();
        let raws = vec![
            Raw::new(a, None, 200_000, false),
            Raw::new(z, None, 66_937, false),
        ];
        let mut b = B::new(raws, 1e-40, 1e-4, false);
        b.clusters[0].center = Some(0);
        b.raws[1].comp = Comparison {
            i: 0,
            index: 1,
            lambda: 0.0,
            hamming: 311,
        };
        b.raws[1].p = 0.0; // what get_pA gives for lambda == 0
        // The scan seeds its minimum from cluster 0's centre, which get_pA
        // forces to 1.0 via the `hamming == 0` branch. Without this the seed
        // is Raw::new's initial 0.0 and nothing can beat it.
        b.raws[0].p = 1.0;
        b.clusters[0].bud_min = None; // the stale cache: offers nothing
        let hit = super::report_missed_bud(&b, 1.0, 1, 1, false);
        assert_eq!(
            hit.map(|h| h.0),
            Some(1),
            "the full scan must still find raw 1"
        );
    }

    /// ... and must stay quiet when there is genuinely nothing to bud.
    #[test]
    fn report_missed_bud_quiet_when_nothing_qualifies() {
        let a = b"ACGTACGTAGCTAGCTAAGGCCTTAGCTAGCTACGTACGTTTGACTGACAGCTTAAGGCCA".to_vec();
        let raws = vec![Raw::new(a, None, 200_000, false)];
        let mut b = B::new(raws, 1e-40, 1e-4, false);
        b.clusters[0].center = Some(0);
        assert!(super::report_missed_bud(&b, 1.0, 1, 1, false).is_none());
    }
}
