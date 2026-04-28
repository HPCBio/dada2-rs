//! Chimera / bimera detection.
//!
//! Ports `chimera.cpp`, excluding Rcpp wrappers and RcppParallel.
//! Parallel processing uses Rayon.
//!
//! A *bimera* is a sequence that can be explained as a chimeric join of two
//! more-abundant "parent" sequences: the left portion of the query matches one
//! parent and the right portion matches another.

use rayon::prelude::*;

use crate::nwalign::{AlignBuffers, align_vectorized_with_buf};

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
pub fn is_bimera(
    sq: &[u8],
    parents: &[&[u8]],
    allow_one_off: bool,
    min_one_off_par_dist: usize,
    match_score: i16,
    mismatch: i16,
    gap_p: i16,
    max_shift: i32,
) -> bool {
    let mut buf = AlignBuffers::new();
    is_bimera_with_buf(
        sq,
        parents,
        allow_one_off,
        min_one_off_par_dist,
        match_score,
        mismatch,
        gap_p,
        max_shift,
        &mut buf,
    )
}

/// Buffer-reusing variant of [`is_bimera`]. See [`AlignBuffers`].
pub fn is_bimera_with_buf(
    sq: &[u8],
    parents: &[&[u8]],
    allow_one_off: bool,
    min_one_off_par_dist: usize,
    match_score: i16,
    mismatch: i16,
    gap_p: i16,
    max_shift: i32,
    buf: &mut AlignBuffers,
) -> bool {
    let sqlen = sq.len();
    let mut max_left = 0usize;
    let mut max_right = 0usize;
    let mut oo_max_left = 0usize;
    let mut oo_max_right = 0usize;
    let mut oo_max_left_oo = 0usize;
    let mut oo_max_right_oo = 0usize;

    for &par in parents {
        align_vectorized_with_buf(sq, par, match_score, mismatch, gap_p, 0, max_shift, buf);
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
    allow_one_off: bool,
    min_one_off_par_dist: usize,
    match_score: i16,
    mismatch: i16,
    gap_p: i16,
    max_shift: i32,
) -> Vec<BimeraFlags> {
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
                        align_vectorized_with_buf(
                            seqs[j],
                            seqs[k],
                            match_score,
                            mismatch,
                            gap_p,
                            0,
                            max_shift,
                            buf,
                        );
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

                if max_left + max_right >= sqlen {
                    nflag += 1;
                } else if allow_one_off
                    && (oo_max_left + oo_max_right_oo >= sqlen
                        || oo_max_left_oo + oo_max_right >= sqlen)
                {
                    nflag += 1;
                }
            } // for i

            BimeraFlags { nflag, nsam }
        })
        .collect()
}
