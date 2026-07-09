//! Chimera / bimera detection.
//!
//! Ports `chimera.cpp`, excluding Rcpp wrappers and RcppParallel.
//! Parallel processing uses Rayon.
//!
//! A *bimera* is a sequence that can be explained as a chimeric join of two
//! more-abundant "parent" sequences: the left portion of the query matches one
//! parent and the right portion matches another.

use rayon::prelude::*;

use crate::nwalign::{
    AlignBackend, AlignBuffers, VectorizedAlignScores, align_vectorized_with_buf,
};
use crate::wfa::{align_wfa_endsfree_with_buf, wfa_cost_cap};

/// Align `sq` against `par` (ends-free) into `buf`, dispatching on `backend`.
/// `Nw` uses the vectorized NW; `Wfa2` the experimental WFA backend. WFA takes
/// `i32` scores, so the `i16` `VectorizedAlignScores` are widened.
#[inline]
fn bimera_align(
    sq: &[u8],
    par: &[u8],
    scores: &VectorizedAlignScores,
    backend: AlignBackend,
    max_steps: i32,
    buf: &mut AlignBuffers,
) {
    match backend {
        AlignBackend::Nw => align_vectorized_with_buf(sq, par, scores, buf),
        AlignBackend::Wfa2 => align_wfa_endsfree_with_buf(
            sq,
            par,
            scores.match_score as i32,
            scores.mismatch as i32,
            scores.gap_p as i32,
            scores.band,
            max_steps,
            buf,
        ),
    }
}

// ---------------------------------------------------------------------------
// Private alignment helpers
// ---------------------------------------------------------------------------

/// Hamming distance between two aligned sequences, ignoring leading and
/// trailing end-gaps.
///
/// End-gaps are positions where either strand is entirely gap-only up to the
/// first non-gap position (left) or after the last non-gap position (right).
/// Equivalent to C++ `get_ham_endsfree`.
fn get_ham_endsfree(s0: &[u8], s1: &[u8]) -> usize {
    let len = s0.len();

    // Find start of internal region.
    let mut i = 0usize;
    let mut gap0 = true;
    let mut gap1 = true;
    while (gap0 || gap1) && i < len {
        gap0 = gap0 && s0[i] == b'-';
        gap1 = gap1 && s1[i] == b'-';
        if gap0 || gap1 {
            i += 1;
        }
    }

    // Find end of internal region.
    let mut j = len as isize - 1;
    gap0 = true;
    gap1 = true;
    while (gap0 || gap1) && j >= i as isize {
        gap0 = gap0 && s0[j as usize] == b'-';
        gap1 = gap1 && s1[j as usize] == b'-';
        if gap0 || gap1 {
            j -= 1;
        }
    }

    if j < i as isize {
        return 0;
    }

    (i..=j as usize).filter(|&p| s0[p] != s1[p]).count()
}

/// Compute left and right overlap lengths between query `al[0]` and parent
/// `al[1]` in a pairwise alignment.
///
/// Returns `(left, right, left_oo, right_oo)`:
/// - `left`:     exact matches from the left before the first mismatch.
/// - `right`:    exact matches from the right before the first mismatch.
/// - `left_oo`:  `left` extended by one mismatch (one-off).
/// - `right_oo`: `right` extended by one mismatch (one-off).
///
/// `left_oo` / `right_oo` are `0` when `allow_one_off` is false.
/// Equivalent to C++ `get_lr`.
fn get_lr(
    s0: &[u8],
    s1: &[u8],
    allow_one_off: bool,
    max_shift: usize,
) -> (usize, usize, usize, usize) {
    let len = s0.len();

    // ---- Left overlap ----
    let mut pos = 0usize;
    // Skip leading gaps in query.
    while pos < len && s0[pos] == b'-' {
        pos += 1;
    }
    // Credit ends-free parent gaps up to max_shift.
    let mut left = 0usize;
    while pos < len && s1[pos] == b'-' && left < max_shift {
        pos += 1;
        left += 1;
    }
    // Count matching positions.
    while pos < len && s0[pos] == s1[pos] {
        pos += 1;
        left += 1;
    }

    let left_oo = if allow_one_off {
        // Step one past the mismatch and count further matches.
        let mut loo = left;
        pos += 1; // step over mismatch
        if pos < len && s0[pos] != b'-' {
            loo += 1; // credit the mismatch position itself if not a gap
        }
        while pos < len && s0[pos] == s1[pos] {
            pos += 1;
            loo += 1;
        }
        loo
    } else {
        0
    };

    // ---- Right overlap ----
    let mut pos = len as isize - 1;
    // Skip trailing gaps in query.
    while pos >= 0 && s0[pos as usize] == b'-' {
        pos -= 1;
    }
    // Credit ends-free parent gaps up to max_shift.
    let mut right = 0usize;
    while pos >= 0 && s1[pos as usize] == b'-' && right < max_shift {
        pos -= 1;
        right += 1;
    }
    // Count matching positions.
    while pos >= 0 && s0[pos as usize] == s1[pos as usize] {
        pos -= 1;
        right += 1;
    }

    let right_oo = if allow_one_off {
        let mut roo = right;
        pos -= 1; // step over mismatch
        if pos >= 0 && s0[pos as usize] != b'-' {
            roo += 1;
        }
        while pos >= 0 && s0[pos as usize] == s1[pos as usize] {
            pos -= 1;
            roo += 1;
        }
        roo
    } else {
        0
    };

    (left, right, left_oo, right_oo)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Alignment and one-off parameters shared by bimera detection functions.
#[derive(Clone, Copy, Debug)]
pub struct BimeraAlignParams {
    pub allow_one_off: bool,
    pub min_one_off_par_dist: usize,
    pub match_score: i16,
    pub mismatch: i16,
    pub gap_p: i16,
    pub max_shift: i32,
    /// Pairwise-alignment backend (issue #49). `Nw` uses the vectorized NW;
    /// `Wfa2` routes through the experimental WFA backend.
    pub backend: AlignBackend,
    /// WFA edit-budget cap (issue #51), in edit operations; `0` = unbounded.
    /// Only meaningful for the `Wfa2` backend. See [`AlignParams::wfa_max_edits`].
    pub wfa_max_edits: i32,
}

/// Determine whether `sq` is a bimera of any pair from `parents`.
///
/// A sequence is a bimera when the left portion of `sq` exactly matches one
/// parent and the right portion exactly matches another, with the two portions
/// together covering the entire length of `sq`.
///
/// `allow_one_off`:        also accept junctions with a single mismatch.
/// `min_one_off_par_dist`: a parent is only considered for one-off detection
///                         if it has at least this many mismatches with `sq`
///                         (prevents self-folding being flagged as chimeric).
///
/// Alignment uses the vectorized banded NW (`end_gap_p = 0`, ends-free).
/// Equivalent to C++ `C_is_bimera`.
#[allow(dead_code)]
pub fn is_bimera(sq: &[u8], parents: &[&[u8]], params: &BimeraAlignParams) -> bool {
    let mut buf = AlignBuffers::new();
    is_bimera_with_buf(sq, parents, params, &mut buf)
}

/// Buffer-reusing variant of [`is_bimera`]. See [`AlignBuffers`].
pub fn is_bimera_with_buf(
    sq: &[u8],
    parents: &[&[u8]],
    params: &BimeraAlignParams,
    buf: &mut AlignBuffers,
) -> bool {
    let BimeraAlignParams {
        allow_one_off,
        min_one_off_par_dist,
        match_score,
        mismatch,
        gap_p,
        max_shift,
        backend,
        wfa_max_edits,
    } = *params;
    let align_scores = VectorizedAlignScores {
        match_score,
        mismatch,
        gap_p,
        end_gap_p: 0,
        band: max_shift,
    };
    let max_steps = wfa_cost_cap(wfa_max_edits, gap_p as i32);
    let sqlen = sq.len();
    let mut max_left = 0usize;
    let mut max_right = 0usize;
    let mut oo_max_left = 0usize;
    let mut oo_max_right = 0usize;
    let mut oo_max_left_oo = 0usize;
    let mut oo_max_right_oo = 0usize;

    for &par in parents {
        bimera_align(sq, par, &align_scores, backend, max_steps, buf);
        let (al0, al1) = buf.alignment();
        let (left, right, left_oo, right_oo) = get_lr(al0, al1, allow_one_off, max_shift as usize);

        // Skip identity / pure-shift / internal-indel parents.
        if left + right >= sqlen {
            continue;
        }

        if left > max_left {
            max_left = left;
        }
        if right > max_right {
            max_right = right;
        }

        if allow_one_off && get_ham_endsfree(al0, al1) >= min_one_off_par_dist {
            if left > oo_max_left {
                oo_max_left = left;
            }
            if right > oo_max_right {
                oo_max_right = right;
            }
            if left_oo > oo_max_left_oo {
                oo_max_left_oo = left_oo;
            }
            if right_oo > oo_max_right_oo {
                oo_max_right_oo = right_oo;
            }
        }

        if max_left + max_right >= sqlen {
            return true;
        }
        if allow_one_off
            && (oo_max_left + oo_max_right_oo >= sqlen || oo_max_left_oo + oo_max_right >= sqlen)
        {
            return true;
        }
    }
    false
}

/// Per-sequence bimera result (one entry per column of the count matrix).
pub struct BimeraFlags {
    /// Number of samples in which a bimeric model was found for this sequence.
    pub nflag: u32,
    /// Number of samples in which this sequence was present (abundance > 0).
    pub nsam: u32,
}

/// Detect bimeras across a multi-sample count table.
///
/// `mat` is a flat column-major count matrix of shape `nrow × ncol`:
/// `mat[i + j * nrow]` is the abundance of sequence `j` in sample `i`.
/// `seqs[j]` is the ASCII sequence for column `j`.
///
/// For each sequence `j`, a potential parent `k` in sample `i` is any
/// sequence where `mat[i + k*nrow] > min_fold * mat[i + j*nrow]` and
/// `mat[i + k*nrow] >= min_abund`.
///
/// Alignments between each (query, parent) pair are cached across samples to
/// avoid redundant computation.
///
/// Runs in parallel over sequences (columns) using Rayon.
/// Equivalent to C++ `C_table_bimera2`.
pub fn table_bimera2(
    mat: &[u32],
    nrow: usize,
    ncol: usize,
    seqs: &[&[u8]],
    min_fold: f64,
    min_abund: u32,
    params: &BimeraAlignParams,
) -> Vec<BimeraFlags> {
    let BimeraAlignParams {
        allow_one_off,
        min_one_off_par_dist,
        match_score,
        mismatch,
        gap_p,
        max_shift,
        backend,
        wfa_max_edits,
    } = *params;
    let align_scores = VectorizedAlignScores {
        match_score,
        mismatch,
        gap_p,
        end_gap_p: 0,
        band: max_shift,
    };
    let max_steps = wfa_cost_cap(wfa_max_edits, gap_p as i32);
    assert_eq!(mat.len(), nrow * ncol, "mat length must be nrow * ncol");
    assert_eq!(seqs.len(), ncol, "seqs length must equal ncol");

    type Cache = Vec<Option<(usize, usize, usize, usize, bool)>>;

    (0..ncol)
        .into_par_iter()
        .map_init(AlignBuffers::new, |buf, j| {
            let sqlen = seqs[j].len();
            let mut nsam = 0u32;
            let mut nflag = 0u32;

            // Cache computed (left, right, left_oo, right_oo, allowed) per parent k.
            // None = not yet computed.
            let mut cache: Cache = vec![None; ncol];

            for i in 0..nrow {
                let j_abund = mat[i + j * nrow];
                if j_abund == 0 {
                    continue;
                }
                nsam += 1;

                let mut max_left = 0usize;
                let mut max_right = 0usize;
                let mut oo_max_left = 0usize;
                let mut oo_max_right = 0usize;
                let mut oo_max_left_oo = 0usize;
                let mut oo_max_right_oo = 0usize;

                for k in 0..ncol {
                    let k_abund = mat[i + k * nrow];
                    if k_abund as f64 <= min_fold * j_abund as f64 || k_abund < min_abund {
                        continue;
                    }

                    // Compute alignment if not cached for this (j, k) pair.
                    if cache[k].is_none() {
                        bimera_align(seqs[j], seqs[k], &align_scores, backend, max_steps, buf);
                        let (al0, al1) = buf.alignment();
                        let (left, right, left_oo, right_oo) =
                            get_lr(al0, al1, allow_one_off, max_shift as usize);
                        let allowed =
                            allow_one_off && get_ham_endsfree(al0, al1) >= min_one_off_par_dist;

                        // Invalidate identity/pure-shift/internal-indel parents.
                        let (l, r, loo, roo) = if left + right < sqlen {
                            (left, right, left_oo, right_oo)
                        } else {
                            (0, 0, 0, 0)
                        };
                        cache[k] = Some((l, r, loo, roo, allowed));
                    }

                    let (l, r, loo, roo, allowed) = cache[k].unwrap();
                    if l > max_left {
                        max_left = l;
                    }
                    if r > max_right {
                        max_right = r;
                    }
                    if allow_one_off && allowed {
                        if l > oo_max_left {
                            oo_max_left = l;
                        }
                        if r > oo_max_right {
                            oo_max_right = r;
                        }
                        if loo > oo_max_left_oo {
                            oo_max_left_oo = loo;
                        }
                        if roo > oo_max_right_oo {
                            oo_max_right_oo = roo;
                        }
                    }
                } // for k

                if max_left + max_right >= sqlen
                    || (allow_one_off
                        && (oo_max_left + oo_max_right_oo >= sqlen
                            || oo_max_left_oo + oo_max_right >= sqlen))
                {
                    nflag += 1;
                }
            } // for i

            BimeraFlags { nflag, nsam }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Coverage diagnostics (higher-order chimera screen)
// ---------------------------------------------------------------------------

/// Per-sequence bimera *coverage* diagnostic.
///
/// The boolean `is_bimera` decision collapses the question "does one junction
/// of two more-abundant parents cover the whole read?" into a yes/no and throws
/// away the underlying coverage. That coverage is exactly the signal for
/// higher-order chimeras (trimeras): a read whose best single-junction model
/// *almost* spans the full length, leaving a small internal gap, is a candidate
/// chimera of three or more parents.
///
/// Coverage here is measured under the **strict** (no one-off) model so the gap
/// boundaries are not blurred by junction mismatches, and under the **pooled**
/// abundance model (abundances summed across samples) for parent selection.
/// All coordinates are in query-base units.
pub struct BimeraDiagnostic {
    /// Strict single-junction decision: `max_left + max_right >= sqlen`.
    pub is_bimera: bool,
    pub sqlen: usize,
    /// Longest strict left overlap with any valid parent.
    pub max_left: usize,
    /// Longest strict right overlap with any valid parent.
    pub max_right: usize,
    /// Parent column index achieving `max_left` (`None` if no valid parent).
    pub best_left_parent: Option<usize>,
    /// Parent column index achieving `max_right`.
    pub best_right_parent: Option<usize>,
    /// Residual uncovered query interval `[gap_start, gap_end)`. Empty when the
    /// read is a clean bimera (`gap_start == gap_end == sqlen`).
    pub gap_start: usize,
    pub gap_end: usize,
    /// Third-parent confirmatory test: the parent (other than the two end
    /// parents) with the fewest mismatches across the residual gap, and that
    /// mismatch count. `gap_mismatches == 0` with a `third_parent` means three
    /// parents fully reconstruct the read — a confirmed trimera. `None` /
    /// large counts suggest the gap is novel biological sequence instead.
    pub third_parent: Option<usize>,
    pub gap_mismatches: usize,
    /// Minimum ends-free Hamming distance to any single candidate parent. The
    /// key disambiguator: a few-SNP *variant* of one abundant parent has a small
    /// value here and is NOT a chimera, even though it leaves a coverage gap; a
    /// genuine mosaic is close to no single parent (large value). See
    /// [`get_ham_endsfree`]. `usize::MAX` when there are no candidate parents.
    pub nearest_parent_dist: usize,
    /// Fewest mismatches across the residual gap achieved by either *end* parent
    /// (`best_left_parent` / `best_right_parent`). The baseline the third parent
    /// must beat: if `gap_mismatches` is not meaningfully below this, the gap is
    /// explained just as well by an end parent (i.e. a variant), not a third
    /// source. `usize::MAX` when not computed (bimera / no parents).
    pub gap_end_parent_mismatches: usize,
}

/// Count mismatches of parent `al_p` against query `al_q` across the query
/// interval `[q_start, q_end)`, walking alignment columns and mapping non-gap
/// query columns to query base indices. Parent gaps (deletions relative to the
/// query) within the region count as mismatches; parent insertions (query gap)
/// are ignored as they do not advance the query coordinate.
fn gap_mismatch_count(al_q: &[u8], al_p: &[u8], q_start: usize, q_end: usize) -> usize {
    let mut qi = 0usize;
    let mut mm = 0usize;
    for c in 0..al_q.len() {
        if al_q[c] != b'-' {
            if qi >= q_start && qi < q_end && al_p[c] != al_q[c] {
                mm += 1;
            }
            qi += 1;
            if qi >= q_end {
                break;
            }
        }
    }
    mm
}

/// Compute [`BimeraDiagnostic`]s for every sequence under the pooled model.
///
/// `pooled[j]` is the summed abundance of sequence `j` across all samples.
/// Parent selection mirrors the pooled `is_bimera` path: a sequence `k` is a
/// candidate parent of `j` when `pooled[k] > min_fold * pooled[j]` and
/// `pooled[k] >= min_abund`.
///
/// Runs two alignment passes per query: pass 1 finds the best left/right
/// coverage (and end parents); pass 2 — only when the read is not already a
/// strict bimera — re-aligns the parents to score the residual gap for a third
/// parent. Parallel over sequences via Rayon.
pub fn pooled_diagnostics(
    pooled: &[u32],
    seqs: &[&[u8]],
    min_fold: f64,
    min_abund: u32,
    params: &BimeraAlignParams,
) -> Vec<BimeraDiagnostic> {
    let BimeraAlignParams {
        match_score,
        mismatch,
        gap_p,
        max_shift,
        backend,
        wfa_max_edits,
        ..
    } = *params;
    let align_scores = VectorizedAlignScores {
        match_score,
        mismatch,
        gap_p,
        end_gap_p: 0,
        band: max_shift,
    };
    let max_steps = wfa_cost_cap(wfa_max_edits, gap_p as i32);
    let ncol = seqs.len();
    assert_eq!(pooled.len(), ncol, "pooled length must equal seqs length");

    (0..ncol)
        .into_par_iter()
        .map_init(AlignBuffers::new, |buf, j| {
            let sq = seqs[j];
            let sqlen = sq.len();
            let abund = pooled[j];

            let mut d = BimeraDiagnostic {
                is_bimera: false,
                sqlen,
                max_left: 0,
                max_right: 0,
                best_left_parent: None,
                best_right_parent: None,
                gap_start: sqlen,
                gap_end: sqlen,
                third_parent: None,
                gap_mismatches: 0,
                nearest_parent_dist: usize::MAX,
                gap_end_parent_mismatches: usize::MAX,
            };

            if abund == 0 {
                return d;
            }
            let parents: Vec<usize> = (0..ncol)
                .filter(|&k| {
                    k != j && pooled[k] as f64 > min_fold * abund as f64 && pooled[k] >= min_abund
                })
                .collect();
            if parents.len() < 2 {
                return d;
            }

            // Pass 1: strict left/right coverage, the parents achieving it, and
            // the nearest single parent (ends-free Hamming) — the disambiguator
            // between a few-SNP variant and a genuine mosaic.
            let mut nearest = usize::MAX;
            for &k in &parents {
                bimera_align(sq, seqs[k], &align_scores, backend, max_steps, buf);
                let (al0, al1) = buf.alignment();
                let dist = get_ham_endsfree(al0, al1);
                if dist < nearest {
                    nearest = dist;
                }
                let (left, right, _, _) = get_lr(al0, al1, false, max_shift as usize);
                // Skip identity / pure-shift / internal-indel parents.
                if left + right >= sqlen {
                    continue;
                }
                if left > d.max_left {
                    d.max_left = left;
                    d.best_left_parent = Some(k);
                }
                if right > d.max_right {
                    d.max_right = right;
                    d.best_right_parent = Some(k);
                }
            }
            d.nearest_parent_dist = nearest;

            d.is_bimera = d.max_left + d.max_right >= sqlen;
            if d.is_bimera {
                return d;
            }
            d.gap_start = d.max_left.min(sqlen);
            d.gap_end = sqlen.saturating_sub(d.max_right).max(d.gap_start);

            // Pass 2: score the residual gap. The two end parents give the
            // baseline a third source must beat; the best non-end parent is the
            // third-parent candidate.
            let mut best_mm = usize::MAX;
            let mut best_k = None;
            let mut end_mm = usize::MAX;
            for &k in &parents {
                bimera_align(sq, seqs[k], &align_scores, backend, max_steps, buf);
                let (al0, al1) = buf.alignment();
                let mm = gap_mismatch_count(al0, al1, d.gap_start, d.gap_end);
                if Some(k) == d.best_left_parent || Some(k) == d.best_right_parent {
                    if mm < end_mm {
                        end_mm = mm;
                    }
                } else if mm < best_mm {
                    best_mm = mm;
                    best_k = Some(k);
                }
            }
            d.third_parent = best_k;
            d.gap_mismatches = if best_mm == usize::MAX { 0 } else { best_mm };
            d.gap_end_parent_mismatches = end_mm;
            d
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> BimeraAlignParams {
        BimeraAlignParams {
            allow_one_off: false,
            min_one_off_par_dist: 4,
            match_score: 5,
            mismatch: -4,
            gap_p: -8,
            max_shift: 16,
            backend: AlignBackend::Nw,
            wfa_max_edits: 0,
        }
    }

    /// A query built as L+M+R, with one high-abundance parent each contributing
    /// L (left), R (right), and M (middle), must be flagged as NOT a strict
    /// bimera but with a residual gap that a distinct third parent fills cleanly.
    #[test]
    fn pooled_diagnostics_detects_trimera() {
        const L: &[u8] = b"ACGTACGTACGTACGTACGT"; // 20
        const M: &[u8] = b"TTTTGGGGCCCCAAAATTTT"; // 20
        const R: &[u8] = b"GATTACAGATTACAGATTAC"; // 20

        let query: Vec<u8> = [L, M, R].concat(); // sqlen 60
        let par_a: Vec<u8> = [
            L,
            &b"AAAAAAAAAAAAAAAAAAAA"[..],
            &b"CCCCCCCCCCCCCCCCCCCC"[..],
        ]
        .concat();
        let par_b: Vec<u8> = [
            &b"GCGCGCGCGCGCGCGCGCGC"[..],
            &b"GCGCGCGCGCGCGCGCGCGC"[..],
            R,
        ]
        .concat();
        let par_c: Vec<u8> = [
            &b"CACACACACACACACACACA"[..],
            M,
            &b"TGTGTGTGTGTGTGTGTGTG"[..],
        ]
        .concat();

        let seqs: Vec<&[u8]> = vec![&query, &par_a, &par_b, &par_c];
        // query low abundance, parents high.
        let pooled = vec![2u32, 100, 100, 100];

        let d = pooled_diagnostics(&pooled, &seqs, 1.5, 2, &params());
        let q = &d[0];

        assert!(
            !q.is_bimera,
            "single junction should not cover the full read"
        );
        assert_eq!(q.max_left, 20, "left parent covers L");
        assert_eq!(q.max_right, 20, "right parent covers R");
        assert_eq!(q.best_left_parent, Some(1));
        assert_eq!(q.best_right_parent, Some(2));
        assert_eq!((q.gap_start, q.gap_end), (20, 40), "gap is the middle M");
        assert_eq!(q.third_parent, Some(3), "parent C fills the gap");
        // C matches the 20-base middle up to minor ends-free boundary slack.
        assert!(
            q.gap_mismatches <= 2,
            "third parent should fill the gap nearly cleanly, got {}",
            q.gap_mismatches
        );
        // A genuine mosaic is close to no single parent...
        assert!(
            q.nearest_parent_dist >= 20,
            "mosaic should not be near any single parent, got {}",
            q.nearest_parent_dist
        );
        // ...and the third parent must beat the end parents in the gap.
        assert!(
            q.gap_mismatches < q.gap_end_parent_mismatches,
            "third parent ({}) must beat end parents ({}) in the gap",
            q.gap_mismatches,
            q.gap_end_parent_mismatches
        );
    }

    /// A few-SNP variant of a single abundant parent leaves a coverage gap but
    /// is NOT a chimera: it sits very close to one parent, and no third parent
    /// explains the gap better than that parent does. (The confound found on the
    /// real nodA data.)
    #[test]
    fn pooled_diagnostics_rejects_snp_variant() {
        // 60 bp parent with low internal self-similarity.
        let par_a: Vec<u8> =
            b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT".to_vec();
        assert_eq!(par_a.len(), 60);
        // Query = par_a with 3 internal substitutions at 20, 30, 40.
        let mut query = par_a.clone();
        for &p in &[20usize, 30, 40] {
            query[p] = if query[p] == b'A' { b'C' } else { b'A' };
        }
        let par_x: Vec<u8> =
            b"TTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGGGGGGGG".to_vec();
        let par_y: Vec<u8> =
            b"GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA".to_vec();

        let seqs: Vec<&[u8]> = vec![&query, &par_a, &par_x, &par_y];
        let pooled = vec![2u32, 100, 100, 100];

        let d = pooled_diagnostics(&pooled, &seqs, 1.5, 2, &params());
        let q = &d[0];

        assert!(!q.is_bimera);
        assert_eq!(
            q.nearest_parent_dist, 3,
            "variant is 3 SNPs from its single parent"
        );
        // The gap is explained at least as well by an end parent as by any
        // third parent, so it must not pass a "third beats end" trimera gate.
        assert!(
            q.gap_mismatches >= q.gap_end_parent_mismatches,
            "no third parent should beat the end parent for a SNP variant \
             (third={}, end={})",
            q.gap_mismatches,
            q.gap_end_parent_mismatches
        );
    }

    /// A genuine bimera (L from one parent, R from another) is flagged
    /// is_bimera with no residual gap.
    #[test]
    fn pooled_diagnostics_flags_plain_bimera() {
        const L: &[u8] = b"ACGTACGTACGTACGTACGT";
        const R: &[u8] = b"GATTACAGATTACAGATTAC";

        let query: Vec<u8> = [L, R].concat(); // 40
        let par_a: Vec<u8> = [L, &b"CCCCCCCCCCCCCCCCCCCC"[..]].concat();
        let par_b: Vec<u8> = [&b"GGGGGGGGGGGGGGGGGGGG"[..], R].concat();

        let seqs: Vec<&[u8]> = vec![&query, &par_a, &par_b];
        let pooled = vec![2u32, 100, 100];

        let d = pooled_diagnostics(&pooled, &seqs, 1.5, 2, &params());
        let q = &d[0];

        assert!(q.is_bimera);
        assert_eq!(q.gap_start, q.gap_end, "no residual gap for a clean bimera");
    }
}
