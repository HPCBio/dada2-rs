//! Needleman-Wunsch alignment and substitution compression.
//!
//! Ports `nwalign_endsfree.cpp` and `nwalign_vectorized.cpp`, excluding all
//! R/Rcpp export wrappers.
//!
//! ## Alignment representation
//! All functions return `[Vec<u8>; 2]` (a pair of integer-encoded, gap-annotated
//! sequences), replacing the C++ `char **al` heap-allocated pair.
//! Gaps are encoded as `b'-'` (byte 45).
//!
//! ## p-matrix encoding
//! Traceback pointer values: `1` = diagonal, `2` = left (gap in s1), `3` = up (gap in s2).

use crate::containers::{Raw, Sub};
use crate::kmers::{assign_kmer, kmer_dist, kord_dist};

/// Sentinel used in `Sub::map` to indicate that a reference position aligns
/// to a gap in the query.  Matches C++ `GAP_GLYPH = 9999`.
pub const GAP_GLYPH: u16 = 9999;

/// Score sentinel for out-of-band DP cells.
const BAND_SENTINEL: i32 = -9999;

// ---------------------------------------------------------------------------
// AlignParams
// ---------------------------------------------------------------------------

/// Parameters controlling alignment method selection in `raw_align`.
#[derive(Clone, Copy)]
pub struct AlignParams {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_p: i32,
    pub homo_gap_p: i32,
    pub use_kmers: bool,
    pub kdist_cutoff: f64,
    /// K-mer size used for the pre-alignment screen and for building the
    /// k-mer / k-order vectors on each `Raw`. Must match the `k` used when
    /// `raw_assign_kmers` populated those vectors (otherwise the distance
    /// indices are garbage).  Valid range: 3..=8.
    pub kmer_size: usize,
    /// Band radius. Negative means unbanded.
    pub band: i32,
    pub vectorized: bool,
    pub gapless: bool,
}

/// Scoring parameters for [`align_vectorized_with_buf`].
///
/// Uses `i16` to match the SIMD-friendly DP tables; `end_gap_p = 0` gives
/// ends-free behaviour, `end_gap_p = gap_p` gives standard NW edge costs.
#[derive(Clone, Copy, Debug)]
pub struct VectorizedAlignScores {
    pub match_score: i16,
    pub mismatch: i16,
    pub gap_p: i16,
    pub end_gap_p: i16,
    pub band: i32,
}

// ---------------------------------------------------------------------------
// AlignBuffers
// ---------------------------------------------------------------------------

/// Reusable scratch buffers for alignment. Pass the same instance across many
/// alignments in a tight loop to avoid re-allocating the DP/traceback matrices.
///
/// One instance is not safe to share across threads; give each worker its own
/// (see `rayon::iter::ParallelIterator::map_init`).
///
/// `al0`/`al1` hold the traceback output of the most recent alignment; the
/// `_with_buf` functions write into them and callers read them in place to
/// avoid per-alignment `Vec<u8>` allocations (~2× per `raw_align`).
#[derive(Default)]
pub struct AlignBuffers {
    // Scalar DP (align_endsfree, align_endsfree_homo, align_standard).
    d32: Vec<i32>,
    p32: Vec<u8>,
    // Vectorized DP (align_vectorized).
    d16: Vec<i16>,
    // Traceback pointers hold only 0..=3, so `u8` halves this matrix's
    // memory traffic vs `i16`. The DP kernel is bandwidth-bound on the
    // streamed `d`/`p` writes, so the narrower store is a direct win.
    p8: Vec<u8>,
    diag_buf: Vec<i16>,
    // Homopolymer masks (align_endsfree_homo).
    homo1: Vec<bool>,
    homo2: Vec<bool>,
    // Traceback output. See struct doc.
    pub al0: Vec<u8>,
    pub al1: Vec<u8>,
}

impl AlignBuffers {
    pub fn new() -> Self {
        Self::default()
    }

    /// Borrow the most recent alignment output as a pair of slices.
    #[inline]
    pub fn alignment(&self) -> (&[u8], &[u8]) {
        (&self.al0, &self.al1)
    }
}

/// Grow `v` to length `n` (no-op if already ≥ n) and fill the first `n`
/// elements with `val`. Reuses the existing allocation when capacity allows.
#[inline]
fn reset_buf<T: Copy>(v: &mut Vec<T>, n: usize, val: T) {
    if v.len() < n {
        v.clear();
        v.resize(n, val);
    } else {
        v[..n].fill(val);
    }
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Return true if byte `b` encodes a nucleotide (A/C/G/T/N = 1..=5).
#[inline]
fn is_nt(b: u8) -> bool {
    matches!(b, 1..=5)
}

/// Compute per-position homopolymer mask into `mask`: true at positions inside
/// a run of three or more identical nucleotides. Resizes `mask` to `seq.len()`.
/// Equivalent to C++ `homo1`/`homo2`.
fn homopolymer_mask_into(seq: &[u8], mask: &mut Vec<bool>) {
    let n = seq.len();
    reset_buf(mask, n, false);
    let mut run_start = 0;
    while run_start < n {
        let mut run_end = run_start + 1;
        while run_end < n && seq[run_end] == seq[run_start] {
            run_end += 1;
        }
        if run_end - run_start >= 3 {
            for k in mask.iter_mut().take(run_end).skip(run_start) {
                *k = true;
            }
        }
        run_start = run_end;
    }
}

/// Trace back through the pointer matrix `p` and fill `al0`/`al1` with the
/// alignment pair. Clears the buffers first so any prior alignment is
/// overwritten; reuses the existing allocation.
/// Shared by `align_endsfree`, `align_endsfree_homo`, and `align_standard`.
#[allow(clippy::too_many_arguments)]
fn traceback_into(
    p: &[u8],
    ncol: usize,
    s1: &[u8],
    s2: &[u8],
    len1: usize,
    len2: usize,
    al0: &mut Vec<u8>,
    al1: &mut Vec<u8>,
) {
    al0.clear();
    al1.clear();
    al0.reserve(len1 + len2);
    al1.reserve(len1 + len2);
    let mut i = len1;
    let mut j = len2;
    while i > 0 || j > 0 {
        match p[i * ncol + j] {
            1 => {
                al0.push(s1[i - 1]);
                al1.push(s2[j - 1]);
                i -= 1;
                j -= 1;
            }
            2 => {
                al0.push(b'-');
                al1.push(s2[j - 1]);
                j -= 1;
            }
            3 => {
                al0.push(s1[i - 1]);
                al1.push(b'-');
                i -= 1;
            }
            _ => panic!("nwalign traceback: invalid pointer value at ({i},{j})"),
        }
    }
    al0.reverse();
    al1.reverse();
}

/// Compute (lband, rband) adjusted for length difference.
fn band_adjust(len1: usize, len2: usize, band: i32) -> (i32, i32) {
    if len2 > len1 {
        (band, band + (len2 - len1) as i32)
    } else if len1 > len2 {
        (band + (len1 - len2) as i32, band)
    } else {
        (band, band)
    }
}

/// Fill band-boundary sentinels into a flat DP matrix.
fn fill_band_sentinels(
    d: &mut [i32],
    ncol: usize,
    len1: usize,
    len2: usize,
    lband: i32,
    rband: i32,
    band: i32,
) {
    if band >= 0 && (band < len1 as i32 || band < len2 as i32) {
        for i in 0..=len1 {
            let li = i as i32 - lband - 1;
            if li >= 0 {
                d[i * ncol + li as usize] = BAND_SENTINEL;
            }
            let ri = i as i32 + rband + 1;
            if ri <= len2 as i32 {
                d[i * ncol + ri as usize] = BAND_SENTINEL;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Standard NW alignment (ends-free)
// ---------------------------------------------------------------------------

/// Banded end-gap-free Needleman-Wunsch alignment.
///
/// End gaps (at the beginning/end of either sequence) are free (score 0).
/// Interior gaps have penalty `gap_p` (should be negative).
/// `band < 0` disables banding.
/// Equivalent to C++ `nwalign_endsfree`.
#[allow(dead_code)]
pub fn align_endsfree(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
) -> [Vec<u8>; 2] {
    let mut buf = AlignBuffers::new();
    align_endsfree_with_buf(s1, s2, match_score, mismatch, gap_p, band, &mut buf);
    [std::mem::take(&mut buf.al0), std::mem::take(&mut buf.al1)]
}

/// Buffer-reusing variant of [`align_endsfree`]. Fills `buf.al0`/`buf.al1`
/// with the alignment pair; read them via `buf.alignment()` or the fields.
pub fn align_endsfree_with_buf(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
    buf: &mut AlignBuffers,
) {
    let len1 = s1.len();
    let len2 = s2.len();
    let ncol = len2 + 1;
    let nrow = len1 + 1;

    reset_buf(&mut buf.d32, nrow * ncol, 0i32);
    reset_buf(&mut buf.p32, nrow * ncol, 0u8);
    {
        let d = &mut buf.d32[..nrow * ncol];
        let p = &mut buf.p32[..nrow * ncol];

        // Initialise edges (ends-free: score 0).
        for slot in p.iter_mut().step_by(ncol).take(len1 + 1) {
            *slot = 3;
        }
        p[..=len2].fill(2);

        let (lband, rband) = band_adjust(len1, len2, band);
        fill_band_sentinels(d, ncol, len1, len2, lband, rband, band);

        for i in 1..=len1 {
            let l = if band >= 0 {
                (i as i32 - lband).max(1) as usize
            } else {
                1
            };
            let r = if band >= 0 {
                (i as i32 + rband).min(len2 as i32) as usize
            } else {
                len2
            };
            for j in l..=r {
                let left = if i == len1 {
                    d[i * ncol + j - 1]
                } else {
                    d[i * ncol + j - 1] + gap_p
                };
                let up = if j == len2 {
                    d[(i - 1) * ncol + j]
                } else {
                    d[(i - 1) * ncol + j] + gap_p
                };
                let diag = d[(i - 1) * ncol + j - 1]
                    + if s1[i - 1] == s2[j - 1] {
                        match_score
                    } else {
                        mismatch
                    };

                if up >= diag && up >= left {
                    d[i * ncol + j] = up;
                    p[i * ncol + j] = 3;
                } else if left >= diag {
                    d[i * ncol + j] = left;
                    p[i * ncol + j] = 2;
                } else {
                    d[i * ncol + j] = diag;
                    p[i * ncol + j] = 1;
                }
            }
        }
    }
    traceback_into(
        &buf.p32[..nrow * ncol],
        ncol,
        s1,
        s2,
        len1,
        len2,
        &mut buf.al0,
        &mut buf.al1,
    );
}

// ---------------------------------------------------------------------------
// Homopolymer-aware end-gap-free NW
// ---------------------------------------------------------------------------

/// Banded end-gap-free NW with position-specific homopolymer gap penalties.
///
/// Gaps inside homopolymer runs (length ≥ 3) use `homo_gap_p` instead of
/// `gap_p`.  Equivalent to C++ `nwalign_endsfree_homo`.
#[allow(dead_code)]
pub fn align_endsfree_homo(s1: &[u8], s2: &[u8], params: &AlignParams) -> [Vec<u8>; 2] {
    let mut buf = AlignBuffers::new();
    align_endsfree_homo_with_buf(s1, s2, params, &mut buf);
    [std::mem::take(&mut buf.al0), std::mem::take(&mut buf.al1)]
}

/// Buffer-reusing variant of [`align_endsfree_homo`]. Fills `buf.al0`/`buf.al1`.
pub fn align_endsfree_homo_with_buf(
    s1: &[u8],
    s2: &[u8],
    params: &AlignParams,
    buf: &mut AlignBuffers,
) {
    let AlignParams {
        match_score,
        mismatch,
        gap_p,
        homo_gap_p,
        band,
        ..
    } = *params;
    let len1 = s1.len();
    let len2 = s2.len();
    let ncol = len2 + 1;
    let nrow = len1 + 1;

    homopolymer_mask_into(s1, &mut buf.homo1);
    homopolymer_mask_into(s2, &mut buf.homo2);
    reset_buf(&mut buf.d32, nrow * ncol, 0i32);
    reset_buf(&mut buf.p32, nrow * ncol, 0u8);
    {
        let homo1 = &buf.homo1[..len1];
        let homo2 = &buf.homo2[..len2];
        let d = &mut buf.d32[..nrow * ncol];
        let p = &mut buf.p32[..nrow * ncol];

        for slot in p.iter_mut().step_by(ncol).take(len1 + 1) {
            *slot = 3;
        }
        p[..=len2].fill(2);

        let (lband, rband) = band_adjust(len1, len2, band);
        fill_band_sentinels(d, ncol, len1, len2, lband, rband, band);

        for i in 1..=len1 {
            let l = if band >= 0 {
                (i as i32 - lband).max(1) as usize
            } else {
                1
            };
            let r = if band >= 0 {
                (i as i32 + rband).min(len2 as i32) as usize
            } else {
                len2
            };
            for j in l..=r {
                let left = if i == len1 {
                    d[i * ncol + j - 1]
                } else if homo2[j - 1] {
                    d[i * ncol + j - 1] + homo_gap_p
                } else {
                    d[i * ncol + j - 1] + gap_p
                };
                let up = if j == len2 {
                    d[(i - 1) * ncol + j]
                } else if homo1[i - 1] {
                    d[(i - 1) * ncol + j] + homo_gap_p
                } else {
                    d[(i - 1) * ncol + j] + gap_p
                };
                let diag = d[(i - 1) * ncol + j - 1]
                    + if s1[i - 1] == s2[j - 1] {
                        match_score
                    } else {
                        mismatch
                    };

                if up >= diag && up >= left {
                    d[i * ncol + j] = up;
                    p[i * ncol + j] = 3;
                } else if left >= diag {
                    d[i * ncol + j] = left;
                    p[i * ncol + j] = 2;
                } else {
                    d[i * ncol + j] = diag;
                    p[i * ncol + j] = 1;
                }
            }
        }
    }
    traceback_into(
        &buf.p32[..nrow * ncol],
        ncol,
        s1,
        s2,
        len1,
        len2,
        &mut buf.al0,
        &mut buf.al1,
    );
}

// ---------------------------------------------------------------------------
// Standard (non-ends-free) NW — not used in core DADA2, included for parity
// ---------------------------------------------------------------------------

/// Standard banded Needleman-Wunsch (edge gaps are penalised).
/// Not used in the core DADA2 algorithm.  Equivalent to C++ `nwalign`.
#[allow(dead_code)]
pub fn align_standard(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
) -> [Vec<u8>; 2] {
    let mut buf = AlignBuffers::new();
    align_standard_with_buf(s1, s2, match_score, mismatch, gap_p, band, &mut buf);
    [std::mem::take(&mut buf.al0), std::mem::take(&mut buf.al1)]
}

/// Buffer-reusing variant of [`align_standard`]. Fills `buf.al0`/`buf.al1`.
#[allow(dead_code)]
pub fn align_standard_with_buf(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
    buf: &mut AlignBuffers,
) {
    let len1 = s1.len();
    let len2 = s2.len();
    let ncol = len2 + 1;
    let nrow = len1 + 1;

    reset_buf(&mut buf.d32, nrow * ncol, 0i32);
    reset_buf(&mut buf.p32, nrow * ncol, 0u8);
    {
        let d = &mut buf.d32[..nrow * ncol];
        let p = &mut buf.p32[..nrow * ncol];

        for i in 1..=len1 {
            d[i * ncol] = d[(i - 1) * ncol] + gap_p;
            p[i * ncol] = 3;
        }
        for j in 1..=len2 {
            d[j] = d[j - 1] + gap_p;
            p[j] = 2;
        }

        let (lband, rband) = band_adjust(len1, len2, band);
        fill_band_sentinels(d, ncol, len1, len2, lband, rband, band);

        for i in 1..=len1 {
            let l = if band >= 0 {
                (i as i32 - lband).max(1) as usize
            } else {
                1
            };
            let r = if band >= 0 {
                (i as i32 + rband).min(len2 as i32) as usize
            } else {
                len2
            };
            for j in l..=r {
                let left = d[i * ncol + j - 1] + gap_p;
                let up = d[(i - 1) * ncol + j] + gap_p;
                let diag = d[(i - 1) * ncol + j - 1]
                    + if s1[i - 1] == s2[j - 1] {
                        match_score
                    } else {
                        mismatch
                    };

                if up >= diag && up >= left {
                    d[i * ncol + j] = up;
                    p[i * ncol + j] = 3;
                } else if left >= diag {
                    d[i * ncol + j] = left;
                    p[i * ncol + j] = 2;
                } else {
                    d[i * ncol + j] = diag;
                    p[i * ncol + j] = 1;
                }
            }
        }
    }
    traceback_into(
        &buf.p32[..nrow * ncol],
        ncol,
        s1,
        s2,
        len1,
        len2,
        &mut buf.al0,
        &mut buf.al1,
    );
}

// ---------------------------------------------------------------------------
// Gapless alignment
// ---------------------------------------------------------------------------

/// Position-by-position alignment without gaps.
/// Shorter sequence is padded with gaps on the right.
/// Equivalent to C++ `nwalign_gapless`.
#[allow(dead_code)]
pub fn align_gapless(s1: &[u8], s2: &[u8]) -> [Vec<u8>; 2] {
    let mut buf = AlignBuffers::new();
    align_gapless_with_buf(s1, s2, &mut buf);
    [std::mem::take(&mut buf.al0), std::mem::take(&mut buf.al1)]
}

/// Buffer-reusing variant of [`align_gapless`]. Fills `buf.al0`/`buf.al1`.
pub fn align_gapless_with_buf(s1: &[u8], s2: &[u8], buf: &mut AlignBuffers) {
    let len = s1.len().max(s2.len());
    buf.al0.clear();
    buf.al1.clear();
    buf.al0.reserve(len);
    buf.al1.reserve(len);
    for i in 0..len {
        buf.al0.push(if i < s1.len() { s1[i] } else { b'-' });
        buf.al1.push(if i < s2.len() { s2[i] } else { b'-' });
    }
}

// ---------------------------------------------------------------------------
// Vectorized (diagonal-banded) NW  — port of nwalign_vectorized2
// ---------------------------------------------------------------------------

/// DP inner loop with `up ≥ left ≥ diag` tie-breaking precedence.
/// Equivalent to C++ `dploop_vec`.
#[inline]
#[allow(clippy::too_many_arguments)]
fn dploop(
    d: &mut [i16],
    p: &mut [u8],
    d_prev: &[i16],
    diag_buf: &[i16],
    left_off: usize,
    up_off: usize,
    out_off: usize,
    col_min: usize,
    n: usize,
    gap_p: i16,
    swap: bool,
) {
    // Bind exact-length subslices up front so the per-cell accesses below
    // carry no bounds checks: the single range-slice panics at most once
    // (before the loop), then the zipped iterators are checkless. This is
    // what lets LLVM auto-vectorize — computed-offset indexing inside the
    // loop (`d_prev[left_off + k]`) inserts a panic branch per access that
    // blocks the vectorizer.
    let left = &d_prev[left_off..left_off + n];
    let up = &d_prev[up_off..up_off + n];
    let diag = &diag_buf[col_min..col_min + n];
    let d_out = &mut d[out_off..out_off + n];
    let p_out = &mut p[out_off..out_off + n];

    // `saturating_add` instead of plain `+` so that at-or-beyond i16 range
    // (long reads, ~>4kbp at default scoring) scores pin to the bounds
    // deterministically rather than wrapping; the outer length guard in
    // raw_align still routes truly long inputs to the i32 path. The LLVM
    // auto-vectorizer lowers saturating i16 adds to `paddsw`/`psubsw`.
    //
    // The `swap` precedence is loop-invariant, so it selects one of two
    // branch-free loop bodies rather than being tested per cell (matches the
    // C++ `dploop_vec` / `dploop_vec_swap` split). Diag wins only when it is
    // strictly greater, preserving the original `entry >= diag` tie-break.
    if swap {
        // precedence: left ≥ up ≥ diag
        for ((((o, pp), &l), &u), &dg) in d_out
            .iter_mut()
            .zip(p_out.iter_mut())
            .zip(left)
            .zip(up)
            .zip(diag)
        {
            let l = l.saturating_add(gap_p);
            let u = u.saturating_add(gap_p);
            let (mut entry, mut pentry) = if l >= u { (l, 2u8) } else { (u, 3u8) };
            if dg > entry {
                entry = dg;
                pentry = 1u8;
            }
            *o = entry;
            *pp = pentry;
        }
    } else {
        // precedence: up ≥ left ≥ diag
        for ((((o, pp), &l), &u), &dg) in d_out
            .iter_mut()
            .zip(p_out.iter_mut())
            .zip(left)
            .zip(up)
            .zip(diag)
        {
            let l = l.saturating_add(gap_p);
            let u = u.saturating_add(gap_p);
            let (mut entry, mut pentry) = if u >= l { (u, 3u8) } else { (l, 2u8) };
            if dg > entry {
                entry = dg;
                pentry = 1u8;
            }
            *o = entry;
            *pp = pentry;
        }
    }
}

/// Cache-friendly diagonal-banded Needleman-Wunsch using `i16` DP tables.
///
/// Processes anti-diagonals (i+j = constant) left to right, giving sequential
/// memory access patterns that auto-vectorise with LLVM.  `end_gap_p = 0`
/// gives ends-free behaviour; `end_gap_p = gap_p` gives standard NW edge costs.
/// Equivalent to C++ `nwalign_vectorized2`.
#[allow(dead_code)]
pub fn align_vectorized(
    s1_in: &[u8],
    s2_in: &[u8],
    scores: &VectorizedAlignScores,
) -> [Vec<u8>; 2] {
    let mut buf = AlignBuffers::new();
    align_vectorized_with_buf(s1_in, s2_in, scores, &mut buf);
    [std::mem::take(&mut buf.al0), std::mem::take(&mut buf.al1)]
}

/// Buffer-reusing variant of [`align_vectorized`]. Fills `buf.al0`/`buf.al1`
/// in the order of the original `s1_in`/`s2_in` inputs (internal swap is
/// undone before returning).
pub fn align_vectorized_with_buf(
    s1_in: &[u8],
    s2_in: &[u8],
    scores: &VectorizedAlignScores,
    buf: &mut AlignBuffers,
) {
    let VectorizedAlignScores {
        match_score,
        mismatch,
        gap_p,
        end_gap_p,
        band,
    } = *scores;
    // Ensure s1 is the shorter sequence; record whether we swapped.
    let swap = s1_in.len() > s2_in.len();
    let (s1, s2) = if swap { (s2_in, s1_in) } else { (s1_in, s2_in) };
    let len1 = s1.len();
    let len2 = s2.len();

    let band = if band < 0 { len2 as i32 } else { band };
    let band_usize = band as usize;

    // Compressed matrix dimensions (diagonal layout).
    // Column index for original cell (i,j): (2*start_col + j - i) / 2
    let start_col = 1 + band_usize.min(len1).div_ceil(2);
    let ncol = 2 + start_col + (band_usize + len2 - len1).min(len2) / 2;
    let nrow = len1 + len2 + 1;

    reset_buf(&mut buf.d16, ncol * nrow, 0i16);
    reset_buf(&mut buf.p8, ncol * nrow, 0u8);
    reset_buf(&mut buf.diag_buf, ncol, 0i16);
    let d = &mut buf.d16[..ncol * nrow];
    let p = &mut buf.p8[..ncol * nrow];
    let diag_buf = &mut buf.diag_buf[..ncol];

    // Sentinel fill: columns 0,1 and ncol-2,ncol-1 in every row act as hard
    // band boundaries.  fill_val is chosen so fill_val + gap_p doesn't overflow.
    let min_score = mismatch.min(gap_p).min(match_score).min(0);
    let fill_val = i16::MIN.wrapping_sub(min_score);
    for row in 0..nrow {
        d[row * ncol] = fill_val;
        d[row * ncol + 1] = fill_val;
        d[row * ncol + ncol - 2] = fill_val;
        d[row * ncol + ncol - 1] = fill_val;
    }

    // Starting cell (0,0) in compressed coordinates.
    d[start_col] = 0;
    p[start_col] = 0;

    // Fill "left column" (gaps in s2 at the start of s1) — ends-free edge.
    {
        let mut row = 1usize;
        let mut col = start_col - 1;
        let mut d_free = end_gap_p;
        let limit = 1 + band_usize.min(len1);
        while row < limit {
            d[row * ncol + col] = d_free;
            p[row * ncol + col] = 3;
            if row.is_multiple_of(2) {
                col = col.saturating_sub(1);
            }
            row += 1;
            d_free = d_free.saturating_add(end_gap_p);
        }
    }

    // Fill "top row" (gaps in s1 at the start of s2) — ends-free edge.
    {
        let mut row = 1usize;
        let mut col = start_col;
        let mut d_free = end_gap_p;
        let limit = 1 + (band_usize + len2 - len1).min(len2);
        while row < limit {
            d[row * ncol + col] = d_free;
            p[row * ncol + col] = 2;
            if row % 2 == 1 {
                col += 1;
            }
            row += 1;
            d_free = d_free.saturating_add(end_gap_p);
        }
    }

    // Main DP: iterate over anti-diagonals (row = i + j).
    let mut row = 2usize;
    let mut col_min = start_col;
    let mut col_max = start_col;
    let mut i_max = 0i64; // 0-indexed into s1 (decrements along anti-diag)
    let mut j_min = 0usize; // 0-indexed into s2 (increments along anti-diag)
    let mut even = true;
    let mut recalc_left = false;
    let mut recalc_right = false;

    while row <= len1 + len2 {
        let n = col_max - col_min + 1;

        // --- Fill diag_buf for this anti-diagonal ---
        // Cell (i,j) in the original NW matrix uses s1[i-1] vs s2[j-1].
        // Here i_max / j_min are 0-indexed, so s1[i_max] vs s2[j_min].
        //
        // The active band can extend one cell past the valid sequence range at
        // the lower-right corner (equal-length sequences, banded mode).  The
        // C++ reference reads s2[len2] in that case (null-terminator UB).
        // We guard with explicit bounds checks and use fill_val so the sentinel
        // columns in d suppress any influence on the traceback.
        {
            let base = (row - 2) * ncol + col_min;
            let d_prev2 = &d[base..base + n];
            let diag_out = &mut diag_buf[col_min..col_min + n];

            // The in-bounds cells of this anti-diagonal form a contiguous range
            // [k_lo, k_hi): `si = i_max - k` in [0, len1) and `sj = j_min + k`
            // in [0, len2). Hoisting that bounds test out of the per-cell loop
            // (computing the range once) lets the hot middle run guard-free with
            // sequential `s1`/`s2` walks — matching R's branchless scalar
            // diag-fill — while staying memory-safe (no null-terminator UB).
            // Cells outside the range take `fill_val`; sentinel `d` columns then
            // suppress any influence on the traceback.
            let k_lo = (i_max - len1 as i64 + 1).clamp(0, n as i64) as usize;
            let k_hi = (i_max + 1)
                .min(len2 as i64 - j_min as i64)
                .clamp(0, n as i64) as usize;
            let k_hi = k_hi.max(k_lo);

            // Boundary cells (si ≥ len1 below k_lo; si < 0 or sj ≥ len2 at/above
            // k_hi): out of range.
            for (dst, &prev) in diag_out[..k_lo].iter_mut().zip(&d_prev2[..k_lo]) {
                *dst = prev.saturating_add(fill_val);
            }
            for (dst, &prev) in diag_out[k_hi..].iter_mut().zip(&d_prev2[k_hi..]) {
                *dst = prev.saturating_add(fill_val);
            }

            // In-bounds middle, guard-free. `s2` is read forward; `s1` in reverse
            // (si decreases as k increases). The slice bounds are provably valid
            // because every k in [k_lo, k_hi) maps to an in-range si/sj.
            if k_lo < k_hi {
                let s1_mid =
                    &s1[(i_max - (k_hi as i64 - 1)) as usize..=(i_max - k_lo as i64) as usize];
                let s2_mid = &s2[j_min + k_lo..j_min + k_hi];
                for (((dst, &prev), &a), &b) in diag_out[k_lo..k_hi]
                    .iter_mut()
                    .zip(&d_prev2[k_lo..k_hi])
                    .zip(s1_mid.iter().rev())
                    .zip(s2_mid.iter())
                {
                    let score = if a == b { match_score } else { mismatch };
                    *dst = prev.saturating_add(score);
                }
            }
        }

        // --- Compute d/p for this row using the previous row ---
        // left  = d[(row-1)*ncol + col_min - even]
        // up    = d[(row-1)*ncol + col_min + 1 - even]
        // out   = d[row*ncol + col_min]
        let even_off = if even { 1 } else { 0 }; // even=true → subtract 1
        let prev_base = (row - 1) * ncol;
        let left_off = prev_base + col_min - even_off;
        let up_off = prev_base + col_min + 1 - even_off;
        let _out_off = row * ncol + col_min;

        // Split d into previous-row (read) and current-row (write) slices.
        let (d_prev_part, d_cur_part) = d.split_at_mut(row * ncol);
        let d_prev = &d_prev_part[prev_base..];
        let d_cur = &mut d_cur_part[..ncol];

        let (p_prev_part, p_cur_part) = p.split_at_mut(row * ncol);
        let _ = p_prev_part; // unused; p_cur accessed via index below
        let p_cur = &mut p_cur_part[..ncol];

        dploop(
            d_cur,
            p_cur,
            d_prev,
            diag_buf,
            left_off - prev_base, // relative to d_prev
            up_off - prev_base,
            col_min, // relative to d_cur / p_cur
            col_min,
            n,
            gap_p,
            swap,
        );

        // --- Band transition: widen active columns at wedge boundaries ---
        // C++ decrements j_min unconditionally (unsigned underflow 0 → SIZE_MAX).
        // The next main-update iteration increments it back to 0.  Using
        // wrapping_sub matches this; the intervening diag_buf fill safely
        // handles the out-of-range index via the bounds guard above.
        if row == band_usize.min(len1) {
            col_min = col_min.saturating_sub(1);
            i_max += 1;
            j_min = j_min.wrapping_sub(1);
        }
        if row == (band_usize + len2 - len1).min(len2) {
            col_max += 1;
        }

        // --- End-gap recalculation (when end_gap_p > gap_p, e.g. ends-free) ---
        // left_off and up_off are absolute indices into d (= prev_base + relative).
        if end_gap_p > gap_p {
            // Left boundary: gap in s2 extending to end of s1
            if recalc_left {
                let d_free = d[left_off].saturating_add(end_gap_p);
                let cur = d[row * ncol + col_min];
                let pcur = p[row * ncol + col_min];
                if d_free > cur {
                    d[row * ncol + col_min] = d_free;
                    p[row * ncol + col_min] = 2;
                } else if !(d_free != cur || !swap && pcur != 1 || swap && pcur == 2) {
                    p[row * ncol + col_min] = 2;
                }
            }
            if i_max == len1 as i64 - 1 {
                recalc_left = true;
            }

            // Right boundary: gap in s1 extending to end of s2
            if recalc_right {
                let d_free = d[up_off + col_max - col_min].saturating_add(end_gap_p);
                let cur = d[row * ncol + col_max];
                let pcur = p[row * ncol + col_max];
                if d_free > cur {
                    d[row * ncol + col_max] = d_free;
                    p[row * ncol + col_max] = 3;
                } else if !(d_free != cur || !swap && pcur == 3 || swap && pcur != 1) {
                    p[row * ncol + col_max] = 3;
                }
            }
            let j_max_1idx = row.div_ceil(2) + col_max - start_col;
            if j_max_1idx == len2 {
                recalc_right = true;
            }
        }

        // --- Update col_min, col_max, i_max, j_min for next anti-diagonal ---
        // j_min uses wrapping arithmetic because the band transition above may
        // set it to usize::MAX (matching C++ unsigned underflow); the increment
        // here wraps it back to 0.
        let band_mod2 = band % 2;
        if (row as i32) < band && (row as i32) < len1 as i32 {
            // Upper triangle for s1
            if even {
                col_min = col_min.saturating_sub(1);
            }
            i_max += 1;
        } else if i_max < (len1 as i64) - 1 {
            // Banded area
            if band_mod2 == 0 {
                if even {
                    j_min = j_min.wrapping_add(1);
                } else {
                    i_max += 1;
                }
            } else if even {
                col_min = col_min.saturating_sub(1);
                i_max += 1;
            } else {
                col_min += 1;
                j_min = j_min.wrapping_add(1);
            }
        } else {
            // Lower triangle for s1
            if !even {
                col_min += 1;
            }
            j_min = j_min.wrapping_add(1);
        }

        let top_limit = (band_usize + len2 - len1).min(len2);
        if row < top_limit {
            if !even {
                col_max += 1;
            }
        } else if row.div_ceil(2) + col_max - start_col < len2 {
            let full_band = band_usize + len2 - len1;
            if full_band.is_multiple_of(2) {
                if even {
                    col_max = col_max.saturating_sub(1);
                } else {
                    col_max += 1;
                }
            }
            // no action for odd full_band
        } else if even {
            col_max = col_max.saturating_sub(1);
        }

        row += 1;
        even = !even;
    }

    // --- Traceback through compressed p matrix ---
    // Reborrow via disjoint fields: p8 (shared) + al0/al1 (mutable each) are
    // three distinct fields of `buf`, so NLL permits holding them simultaneously.
    let p_ro = &buf.p8[..ncol * nrow];
    let al0 = &mut buf.al0;
    let al1 = &mut buf.al1;
    al0.clear();
    al1.clear();
    al0.reserve(len1 + len2);
    al1.reserve(len1 + len2);

    let mut i = len1;
    let mut j = len2;
    while i > 0 || j > 0 {
        // Compressed column: (2*start_col + j - i) / 2  (C-style truncating division).
        // j - i can be odd, which is why C++ just truncates — the correct column
        // for (i,j) in the anti-diagonal layout is floor((2*start_col + j - i) / 2).
        let col_signed = 2 * start_col as i64 + j as i64 - i as i64;
        debug_assert!(
            col_signed >= 0,
            "vectorized traceback: col_signed={col_signed} < 0 at i={i} j={j}"
        );
        let col = (col_signed / 2) as usize;
        match p_ro[(i + j) * ncol + col] {
            1 => {
                al0.push(s1[i - 1]);
                al1.push(s2[j - 1]);
                i -= 1;
                j -= 1;
            }
            2 => {
                al0.push(b'-');
                al1.push(s2[j - 1]);
                j -= 1;
            }
            3 => {
                al0.push(s1[i - 1]);
                al1.push(b'-');
                i -= 1;
            }
            v => panic!("vectorized traceback: invalid pointer {v} at i={i} j={j}"),
        }
    }
    al0.reverse();
    al1.reverse();

    // Restore original input ordering if we swapped the DP inputs.
    if swap {
        std::mem::swap(&mut buf.al0, &mut buf.al1);
    }
}

// ---------------------------------------------------------------------------
// WFA backend (experimental) — wavefront alignment via wfa2lib-rs
// ---------------------------------------------------------------------------
//
// Wraps the pure-Rust WFA aligner (HPCBio fork of COMBINE-lab/wfa2lib-rs) so it
// can stand in for `align_endsfree`. WFA *minimises a penalty* (cost) where the
// match cost is ≤ 0, whereas DADA2 *maximises a score* (match = +5, mismatch =
// −4, gap = −8). The two are equivalent under sign inversion (score → cost):
//
//     match_      = -match_score   (≤ 0, e.g. -5)
//     mismatch    = -mismatch      (> 0, e.g.  4)
//     gap_extension = -gap_p       (> 0, e.g.  8)   gap_opening = 0  (linear gap)
//
// wfa2lib-rs's `new_affine` applies the Eizenga adjustment internally when
// `match_ < 0`, so a non-zero match score reproduces the same optimum as the
// scalar DP. Sequences are passed as the raw 1..=5 nt encoding: WFA matches by
// byte equality and its EOS sentinels (`b'!'`=33, `b'?'`=63) never collide with
// 1..=5, so no ASCII round-trip is needed.

use wfa2lib_rs::aligner::{AffineAligner, AlignmentScope};
use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

/// Convert a WFA CIGAR (`M`/`X`/`I`/`D` ops, pattern=s1, text=s2) into the
/// gap-annotated `[al0, al1]` pair the rest of the module consumes.
///
/// Op semantics (from `Cigar::check_alignment`): `M`/`X` consume one base from
/// each strand; `I` advances the text only (gap in pattern → `al0`); `D`
/// advances the pattern only (gap in text → `al1`).
#[allow(dead_code)]
fn cigar_to_alignment_into(ops: &[u8], s1: &[u8], s2: &[u8], al0: &mut Vec<u8>, al1: &mut Vec<u8>) {
    al0.clear();
    al1.clear();
    al0.reserve(ops.len());
    al1.reserve(ops.len());
    let mut p = 0usize; // pattern (s1) position
    let mut t = 0usize; // text (s2) position
    for &op in ops {
        match op {
            b'M' | b'X' => {
                al0.push(s1[p]);
                al1.push(s2[t]);
                p += 1;
                t += 1;
            }
            b'I' => {
                al0.push(b'-');
                al1.push(s2[t]);
                t += 1;
            }
            b'D' => {
                al0.push(s1[p]);
                al1.push(b'-');
                p += 1;
            }
            v => panic!("WFA cigar_to_alignment: unknown op {v} ({})", v as char),
        }
    }
}

/// Ends-free Needleman-Wunsch via the WFA backend. Fills `buf.al0`/`buf.al1`
/// with the same alignment-pair representation as [`align_endsfree`].
///
/// `match_score`/`mismatch`/`gap_p` use the DADA2 score convention (match > 0,
/// the rest < 0); they are converted to WFA penalties internally.
///
/// NOTE: constructs a fresh `AffineAligner` per call. The hot path will want a
/// cached aligner in `AlignBuffers`; this is correctness-first.
#[allow(dead_code)]
pub fn align_wfa_endsfree_with_buf(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    buf: &mut AlignBuffers,
) {
    let penalties = WavefrontPenalties::new_affine(AffinePenalties {
        match_: -match_score,
        mismatch: -mismatch,
        gap_opening: 0,
        gap_extension: -gap_p,
    });
    let mut aligner = AffineAligner::new(penalties);
    // Compute the full CIGAR, not just the score (default is ComputeScore, which
    // leaves the cigar empty). The field is pub; a public setter would be a good
    // upstream addition.
    aligner.alignment_scope = AlignmentScope::ComputeAlignment;
    // Full ends-free: leading/trailing gaps on either strand are free.
    aligner.set_alignment_free_ends(
        s1.len() as i32,
        s1.len() as i32,
        s2.len() as i32,
        s2.len() as i32,
    );
    aligner.align_endsfree(s1, s2);
    let cigar = aligner.cigar();
    cigar_to_alignment_into(cigar.operations_slice(), s1, s2, &mut buf.al0, &mut buf.al1);
}

// ---------------------------------------------------------------------------
// Substitution compression
// ---------------------------------------------------------------------------

/// Convert an alignment pair into a `Sub` (compressed substitution record).
///
/// Records substitutions of `al1` relative to `al0`, ignoring positions
/// where either strand has an N (encoded 5).  `Sub::q0`/`Sub::q1` are left
/// empty; fill them via `sub_new` if quality scores are needed.
/// Equivalent to C++ `al2subs`.
pub fn al2subs(al0: &[u8], al1: &[u8]) -> Sub {
    let alen = al0.len();
    debug_assert_eq!(al0.len(), al1.len());

    // First pass: count reference length and substitution count.
    let mut len0 = 0u32;
    let mut nsubs = 0usize;
    for i in 0..alen {
        let nt0 = is_nt(al0[i]);
        let nt1 = is_nt(al1[i]);
        if nt0 {
            len0 += 1;
        }
        if nt0 && nt1 && al0[i] != al1[i] && al0[i] != 5 && al1[i] != 5 {
            nsubs += 1;
        }
    }

    let mut map = vec![GAP_GLYPH; len0 as usize];
    let mut pos = Vec::with_capacity(nsubs);
    let mut nt0_vec = Vec::with_capacity(nsubs);
    let mut nt1_vec = Vec::with_capacity(nsubs);

    // Second pass: fill map and substitution arrays.
    let mut i0: i64 = -1;
    let mut i1: i64 = -1;
    for i in 0..alen {
        let nt0 = is_nt(al0[i]);
        let nt1 = is_nt(al1[i]);
        if nt0 {
            i0 += 1;
        }
        if nt1 {
            i1 += 1;
        }

        if nt0 {
            map[i0 as usize] = if nt1 { i1 as u16 } else { GAP_GLYPH };
        }
        if nt0 && nt1 && al0[i] != al1[i] && al0[i] != 5 && al1[i] != 5 {
            pos.push(i0 as u16);
            nt0_vec.push(al0[i]);
            nt1_vec.push(al1[i]);
        }
    }

    Sub {
        len0,
        map,
        pos,
        nt0: nt0_vec,
        nt1: nt1_vec,
        q0: Vec::new(),
        q1: Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Dispatcher and sub_new
// ---------------------------------------------------------------------------

/// Select and run the appropriate alignment for two `Raw` objects.
///
/// Returns `None` if k-mer screening determines the sequences are too
/// dissimilar to be worth aligning (i.e. they will produce a NULL Sub).
/// Equivalent to C++ `raw_align`.
#[allow(dead_code)]
pub fn raw_align(raw1: &Raw, raw2: &Raw, p: &AlignParams) -> Option<[Vec<u8>; 2]> {
    let mut buf = AlignBuffers::new();
    raw_align_with_buf(raw1, raw2, p, &mut buf)?;
    Some([std::mem::take(&mut buf.al0), std::mem::take(&mut buf.al1)])
}

/// Buffer-reusing variant of [`raw_align`]. On `Some(())`, the alignment is
/// in `buf.al0`/`buf.al1` (read via `buf.alignment()`).
pub fn raw_align_with_buf(
    raw1: &Raw,
    raw2: &Raw,
    p: &AlignParams,
    buf: &mut AlignBuffers,
) -> Option<()> {
    // --- K-mer screening ---
    let mut kdist = 0.0f64;
    let mut kodist = -1.0f64; // sentinel: different from kdist when use_kmers=false

    if p.use_kmers {
        let k = p.kmer_size;
        // Prefer 8-bit kmer distance; fall back to 16-bit on overflow.
        kdist = match (&raw1.kmer8, &raw2.kmer8) {
            (Some(k1), Some(k2)) => {
                let d8 = k1.dist8(k2, raw1.len(), raw2.len(), k);
                if d8 < 0.0 {
                    // Overflow (a k-mer occurs ≥255× in both seqs): fall back to
                    // the exact 16-bit distance. The u16 vectors are not kept
                    // resident (issue #32), so recompute them from sequence here
                    // — this path is essentially never hit for amplicon data.
                    let v1 = assign_kmer(&raw1.seq, k);
                    let v2 = assign_kmer(&raw2.seq, k);
                    kmer_dist(&v1, raw1.len(), &v2, raw2.len(), k)
                } else {
                    d8
                }
            }
            // No 8-bit vectors built (e.g. cluster-center raws): no k-mer screen,
            // exactly as before — previously both kmer8 and the u16 kmer were
            // absent together, yielding kdist = 0.0.
            _ => 0.0,
        };

        if p.gapless
            && let (Some(o1), Some(o2)) = (&raw1.kord, &raw2.kord)
        {
            kodist = kord_dist(o1, raw1.len(), o2, raw2.len(), k);
        }
    }

    if p.use_kmers && kdist > p.kdist_cutoff {
        return None; // Outside k-mer distance threshold → NULL alignment.
    }

    // --- Method selection ---
    if p.band == 0 || (p.gapless && (kodist - kdist).abs() < f64::EPSILON) {
        align_gapless_with_buf(&raw1.seq, &raw2.seq, buf);
        return Some(());
    }
    // Homopolymer gapping: when a distinct homopolymer gap penalty is in
    // effect, the alignment must use the homopolymer-aware scalar DP. The
    // vectorized aligner has no homopolymer-gap support, so it is disabled in
    // this mode — mirroring R DADA2, which forces `VECTORIZED_ALIGNMENT <-
    // FALSE` whenever `HOMOPOLYMER_GAP_PENALTY != GAP_PENALTY` (dada.R:229-230).
    // Without this, a `--homo-gap-p` setting would be silently ignored by the
    // vectorized path and diverge from R.
    let use_homo = p.homo_gap_p != p.gap_p && p.homo_gap_p <= 0;
    // Long-read guard: align_vectorized uses i16 DP tables. With the default
    // DADA2 scoring (match=5, mismatch=-4, gap_p=-8) cumulative scores can
    // approach ±8·N, so we must fall back to the i32 path before overflow
    // can distort the optimum. 3500 bp leaves ~10% headroom in i16 range.
    const VECTORIZED_MAX_LEN: usize = 3500;
    let too_long = raw1.len() > VECTORIZED_MAX_LEN || raw2.len() > VECTORIZED_MAX_LEN;
    if p.vectorized && !too_long && !use_homo {
        align_vectorized_with_buf(
            &raw1.seq,
            &raw2.seq,
            &VectorizedAlignScores {
                match_score: p.match_score as i16,
                mismatch: p.mismatch as i16,
                gap_p: p.gap_p as i16,
                end_gap_p: 0,
                band: p.band,
            },
            buf,
        );
        return Some(());
    }
    if use_homo {
        align_endsfree_homo_with_buf(&raw1.seq, &raw2.seq, p, buf);
        return Some(());
    }
    align_endsfree_with_buf(
        &raw1.seq,
        &raw2.seq,
        p.match_score,
        p.mismatch,
        p.gap_p,
        p.band,
        buf,
    );
    Some(())
}

/// Align two `Raw` objects and return the compressed substitution record,
/// with quality scores filled in when both Raws carry quality data.
///
/// Returns `None` when the k-mer screen rejects the pair (equivalent to a
/// NULL Sub in the C++ code).
/// Equivalent to C++ `sub_new`.
#[allow(dead_code)]
pub fn sub_new(raw0: &Raw, raw1: &Raw, params: &AlignParams) -> Option<Sub> {
    let mut buf = AlignBuffers::new();
    sub_new_with_buf(raw0, raw1, params, &mut buf)
}

/// Buffer-reusing variant of [`sub_new`]. See [`AlignBuffers`].
pub fn sub_new_with_buf(
    raw0: &Raw,
    raw1: &Raw,
    params: &AlignParams,
    buf: &mut AlignBuffers,
) -> Option<Sub> {
    raw_align_with_buf(raw0, raw1, params, buf)?;
    let mut sub = al2subs(&buf.al0, &buf.al1);

    if let (Some(q0), Some(q1)) = (&raw0.qual, &raw1.qual) {
        sub.q0 = sub.pos.iter().map(|&pos| q0[pos as usize]).collect();
        sub.q1 = sub
            .pos
            .iter()
            .map(|&pos| q1[sub.map[pos as usize] as usize])
            .collect();
    }
    Some(sub)
}

#[cfg(test)]
mod bench_align {
    use super::*;

    /// Isolated kernel micro-benchmark: time `align_vectorized_with_buf` on one
    /// fixed sequence pair, many iterations, with a reused buffer. Strips away
    /// every pipeline confound (k-mer screen, greedy, threading, counting) so
    /// the per-alignment number is directly comparable to R's `nwalign(vec=TRUE)`.
    ///
    /// Reads a 2-record FASTA from the `BENCH_PAIR` env var (default
    /// `/tmp/bench_pair.fasta`). Run explicitly:
    ///   BENCH_PAIR=/tmp/bench_pair.fasta cargo test --release \
    ///     -- --ignored bench_align_vectorized --nocapture
    #[test]
    #[ignore]
    fn bench_align_vectorized() {
        let path = std::env::var("BENCH_PAIR").unwrap_or_else(|_| "/tmp/bench_pair.fasta".into());
        let text = std::fs::read_to_string(&path).expect("read BENCH_PAIR fasta");
        let seqs: Vec<Vec<u8>> = text
            .split('>')
            .filter(|r| !r.trim().is_empty())
            .map(|rec| {
                let body: String = rec.lines().skip(1).collect();
                body.bytes()
                    .filter_map(|b| match b {
                        b'A' | b'a' => Some(1u8),
                        b'C' | b'c' => Some(2),
                        b'G' | b'g' => Some(3),
                        b'T' | b't' => Some(4),
                        b'N' | b'n' => Some(5),
                        _ => None,
                    })
                    .collect()
            })
            .collect();
        assert!(seqs.len() >= 2, "need >=2 records in {path}");
        let (s1, s2) = (&seqs[0], &seqs[1]);

        let scores = VectorizedAlignScores {
            match_score: 5,
            mismatch: -4,
            gap_p: -8,
            end_gap_p: 0,
            band: 32,
        };
        let mut buf = AlignBuffers::new();

        // Warmup.
        for _ in 0..200 {
            align_vectorized_with_buf(s1, s2, &scores, &mut buf);
        }

        let iters = 20_000usize;
        let t0 = std::time::Instant::now();
        for _ in 0..iters {
            align_vectorized_with_buf(s1, s2, &scores, &mut buf);
            std::hint::black_box(buf.alignment());
        }
        let dt = t0.elapsed();

        let us = dt.as_secs_f64() * 1e6 / iters as f64;
        let aps = iters as f64 / dt.as_secs_f64();
        println!(
            "  align_vectorized: len1={} len2={} band=32  {:.2} us/align  {:.0} aligns/s  ({iters} iters)",
            s1.len(),
            s2.len(),
            us,
            aps,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn score_alignment(al: &[Vec<u8>; 2], match_score: i32, mismatch: i32, gap_p: i32) -> i32 {
        let al0 = &al[0];
        let al1 = &al[1];
        let n = al0.len();
        let mut score = 0i32;
        let mut i0 = false; // was previous al0[k] a gap?
        let mut i1 = false;
        for k in 0..n {
            let g0 = al0[k] == b'-';
            let g1 = al1[k] == b'-';
            if !g0 && !g1 {
                score += if al0[k] == al1[k] {
                    match_score
                } else {
                    mismatch
                };
                i0 = false;
                i1 = false;
            } else if g0 {
                // gap in s1 (end-gap free if at start or end)
                let _ = i1;
                i1 = true;
                i0 = false;
                let _ = i0; // end-gap: skip penalty (ends-free)
                // Only penalise interior gaps
                let at_start = al0[..k].iter().all(|&b| b == b'-');
                let at_end = al0[k + 1..].iter().all(|&b| b == b'-');
                if !at_start && !at_end {
                    score += gap_p;
                }
            } else {
                i0 = true;
                i1 = false;
                let at_start = al1[..k].iter().all(|&b| b == b'-');
                let at_end = al1[k + 1..].iter().all(|&b| b == b'-');
                if !at_start && !at_end {
                    score += gap_p;
                }
            }
        }
        score
    }

    /// Encode a DNA string (A/C/G/T) to u8 nt-index (1-4).
    fn encode(seq: &str) -> Vec<u8> {
        seq.bytes()
            .map(|b| match b {
                b'A' | b'a' => 1,
                b'C' | b'c' => 2,
                b'G' | b'g' => 3,
                b'T' | b't' => 4,
                _ => 5,
            })
            .collect()
    }

    fn decode_al(al: &[Vec<u8>; 2]) -> (String, String) {
        let to_str = |v: &Vec<u8>| {
            v.iter()
                .map(|&b| match b {
                    1 => 'A',
                    2 => 'C',
                    3 => 'G',
                    4 => 'T',
                    b'-' => '-',
                    _ => 'N',
                })
                .collect::<String>()
        };
        (to_str(&al[0]), to_str(&al[1]))
    }

    const MATCH: i32 = 5;
    const MM: i32 = -4;
    const GAP: i32 = -8;
    const BAND: i32 = 16;

    fn check_endsfree_score(s1: &[u8], s2: &[u8], expected: i32, label: &str) {
        let al = align_endsfree(s1, s2, MATCH, MM, GAP, BAND);
        let got = score_alignment(&al, MATCH, MM, GAP);
        assert_eq!(got, expected, "{label}: endsfree score mismatch");
    }

    /// Assert that align_vectorized produces the same optimal score as align_endsfree.
    fn compare_alignments(s1: &[u8], s2: &[u8], label: &str) {
        let ef = align_endsfree(s1, s2, MATCH, MM, GAP, BAND);
        let ve = align_vectorized(
            s1,
            s2,
            &VectorizedAlignScores {
                match_score: MATCH as i16,
                mismatch: MM as i16,
                gap_p: GAP as i16,
                end_gap_p: 0,
                band: BAND,
            },
        );

        let score_ef = score_alignment(&ef, MATCH, MM, GAP);
        let score_ve = score_alignment(&ve, MATCH, MM, GAP);

        if score_ef != score_ve {
            let (ef0, ef1) = decode_al(&ef);
            let (ve0, ve1) = decode_al(&ve);
            panic!(
                "{label}: score mismatch: endsfree={score_ef} vectorized={score_ve}\n  EF: {ef0}\n      {ef1}\n  VE: {ve0}\n      {ve1}"
            );
        }
    }

    /// Assert the WFA backend produces the same optimal score as `align_endsfree`.
    fn compare_wfa(s1: &[u8], s2: &[u8], label: &str) {
        let ef = align_endsfree(s1, s2, MATCH, MM, GAP, BAND);
        let mut buf = AlignBuffers::new();
        align_wfa_endsfree_with_buf(s1, s2, MATCH, MM, GAP, &mut buf);
        let wfa = [buf.al0.clone(), buf.al1.clone()];

        let score_ef = score_alignment(&ef, MATCH, MM, GAP);
        let score_wfa = score_alignment(&wfa, MATCH, MM, GAP);
        if score_ef != score_wfa {
            let (e0, e1) = decode_al(&ef);
            let (w0, w1) = decode_al(&wfa);
            panic!(
                "{label}: WFA score mismatch: endsfree={score_ef} wfa={score_wfa}\n  EF: {e0}\n      {e1}\n  WFA: {w0}\n       {w1}"
            );
        }
    }

    #[test]
    fn test_wfa_vs_endsfree_identical() {
        let s = encode("ACGTACGTACGT");
        compare_wfa(&s, &s, "identical-short");
    }

    #[test]
    fn test_wfa_vs_endsfree_one_sub() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTTCGTACGT");
        compare_wfa(&s1, &s2, "one-sub");
    }

    #[test]
    fn test_wfa_vs_endsfree_one_gap() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGACGTACGT");
        compare_wfa(&s1, &s2, "one-gap");
    }

    #[test]
    fn test_wfa_vs_endsfree_different_lengths() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTACGTACGTAC");
        compare_wfa(&s1, &s2, "diff-len-short");
    }

    /// Find and print the first few low-edit pairs where WFA disagrees with
    /// `align_endsfree`, with both alignments, to diagnose the cause.
    ///   cargo test --release --bins wfa_diagnose -- --ignored --nocapture
    #[test]
    #[ignore]
    fn wfa_diagnose() {
        let nts = [1u8, 2, 3, 4];
        let mut st: u64 = 0x1234_5678_9ABC_DEF0;
        let mut rng = |st: &mut u64, m: usize| {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*st >> 33) as usize) % m
        };
        let mut shown = 0;
        for _ in 0..200_000 {
            if shown >= 5 {
                break;
            }
            let len = 40 + rng(&mut st, 60);
            let s1: Vec<u8> = (0..len).map(|_| nts[rng(&mut st, 4)]).collect();
            let mut s2 = s1.clone();
            let nedits = 1 + rng(&mut st, 1);
            for _ in 0..nedits {
                let p = rng(&mut st, s2.len());
                match rng(&mut st, 3) {
                    0 => s2[p] = nts[rng(&mut st, 4)],
                    1 => s2.insert(p, nts[rng(&mut st, 4)]),
                    _ => {
                        s2.remove(p);
                    }
                }
            }
            let ef = align_endsfree(&s1, &s2, 5, -4, -8, -1);
            let mut buf = AlignBuffers::new();
            align_wfa_endsfree_with_buf(&s1, &s2, 5, -4, -8, &mut buf);
            let wfa = [buf.al0.clone(), buf.al1.clone()];
            let sef = score_alignment(&ef, 5, -4, -8);
            let swfa = score_alignment(&wfa, 5, -4, -8);
            if sef != swfa {
                let (e0, e1) = decode_al(&ef);
                let (w0, w1) = decode_al(&wfa);
                println!(
                    "--- disagreement (edits={nedits}) ef={sef} wfa={swfa} len1={} len2={} ---",
                    s1.len(),
                    s2.len()
                );
                println!("EF : {e0}");
                println!("     {e1}");
                println!("WFA: {w0}");
                println!("     {w1}");
                shown += 1;
            }
        }
    }

    /// Randomized parity stress test: WFA must match `align_endsfree`'s optimal
    /// score across many random pairs (varied lengths, indels). Run explicitly:
    ///   cargo test --release -- --ignored sweep_wfa_parity --nocapture
    #[test]
    #[ignore]
    fn sweep_wfa_parity() {
        let nts = [1u8, 2, 3, 4];
        let mut st: u64 = 0x1234_5678_9ABC_DEF0;
        let mut rng = |st: &mut u64, m: usize| {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*st >> 33) as usize) % m
        };
        // Bucket failures by number of edits to test the hypothesis that
        // divergence (not a logic bug) drives disagreement. DADA2 only aligns
        // k-mer-screened, very-similar pairs, so the low-edit buckets are the
        // ones that matter.
        let mut fails_by_edits = [0u32; 9];
        let mut total_by_edits = [0u32; 9];
        for _ in 0..10_000 {
            let len = 40 + rng(&mut st, 210); // 40..250, amplicon-ish
            let s1: Vec<u8> = (0..len).map(|_| nts[rng(&mut st, 4)]).collect();
            let mut s2 = s1.clone();
            let nedits = rng(&mut st, 9); // 0..8
            for _ in 0..nedits {
                if s2.is_empty() {
                    break;
                }
                let p = rng(&mut st, s2.len());
                match rng(&mut st, 3) {
                    0 => s2[p] = nts[rng(&mut st, 4)],
                    1 => s2.insert(p, nts[rng(&mut st, 4)]),
                    _ => {
                        s2.remove(p);
                    }
                }
            }
            if s2.len() < 3 {
                continue;
            }
            let ef = align_endsfree(&s1, &s2, 5, -4, -8, -1);
            let mut buf = AlignBuffers::new();
            align_wfa_endsfree_with_buf(&s1, &s2, 5, -4, -8, &mut buf);
            let wfa = [buf.al0.clone(), buf.al1.clone()];
            total_by_edits[nedits] += 1;
            if score_alignment(&ef, 5, -4, -8) != score_alignment(&wfa, 5, -4, -8) {
                fails_by_edits[nedits] += 1;
            }
        }
        for e in 0..9 {
            println!(
                "  edits={e}: {}/{} disagree",
                fails_by_edits[e], total_by_edits[e]
            );
        }
        // INFORMATIONAL (not asserted): documents the known divergence between
        // the WFA backend and `align_endsfree`. WFA-via-Eizenga ends-free does
        // not reproduce DADA2's exact score-maximizing optimum at gap/homopolymer
        // boundaries and sequence ends (free end-gap crediting differs), so it is
        // NOT yet a drop-in for align_endsfree. Tracked for the fork's ends-free
        // handling. See wfa_diagnose for concrete counterexamples.
        let low_edit_fails: u32 = fails_by_edits[..=4].iter().sum();
        let total_fails: u32 = fails_by_edits.iter().sum();
        println!(
            "  TOTAL: {total_fails} disagree ({low_edit_fails} at <=4 edits) \
             — known WFA ends-free divergence, see doc comment"
        );
    }

    #[test]
    fn test_align_endsfree_identical_short() {
        let s = encode("ACGTACGTACGT");
        check_endsfree_score(&s, &s, 60, "identical-short"); // 12 matches * 5
    }

    #[test]
    fn test_align_endsfree_one_sub() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTTCGTACGT"); // one sub at pos 4
        check_endsfree_score(&s1, &s2, 11 * 5 + (-4), "one-sub");
    }

    #[test]
    fn test_align_endsfree_one_gap() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGACGTACGT"); // one deletion (11 matches, 1 interior gap)
        check_endsfree_score(&s1, &s2, 11 * 5 + (-8), "one-gap");
    }

    #[test]
    fn test_vectorized_vs_endsfree_identical() {
        let s = encode("ACGTACGTACGT");
        compare_alignments(&s, &s, "identical-short");
    }

    #[test]
    fn test_vectorized_vs_endsfree_one_sub() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTTCGTACGT");
        compare_alignments(&s1, &s2, "one-sub");
    }

    #[test]
    fn test_vectorized_vs_endsfree_one_gap() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGACGTACGT");
        compare_alignments(&s1, &s2, "one-gap");
    }

    #[test]
    fn test_vectorized_vs_endsfree_equal_length_240() {
        let nts: [u8; 4] = [1, 2, 3, 4];
        let mut state: u64 = 99991;
        let next_nt = |st: &mut u64| -> u8 {
            *st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            nts[((*st >> 33) as usize) % 4]
        };
        let s1: Vec<u8> = (0..240).map(|_| next_nt(&mut state)).collect();
        compare_alignments(&s1, &s1.clone(), "identical-240");

        let mut s2 = s1.clone();
        s2[239] ^= 3;
        compare_alignments(&s1, &s2, "last-mismatch-240");
    }

    #[test]
    fn test_vectorized_vs_endsfree_divergent_240() {
        let nts: [u8; 4] = [1, 2, 3, 4];
        let mut state: u64 = 12345;
        let next_nt = |st: &mut u64| -> u8 {
            *st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            nts[((*st >> 33) as usize) % 4]
        };
        let s1: Vec<u8> = (0..240).map(|_| next_nt(&mut state)).collect();
        let mut s2 = s1.clone();
        for i in (0..240).step_by(50) {
            s2[i] = nts[(s2[i] as usize) % 4 ^ 1];
        }
        compare_alignments(&s1, &s2, "divergent-240");
    }

    #[test]
    fn test_vectorized_vs_endsfree_different_lengths() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTACGTACGTAC"); // s2 longer
        compare_alignments(&s1, &s2, "diff-len-short");

        let nts: [u8; 4] = [1, 2, 3, 4];
        let mut state: u64 = 77777;
        let next_nt = |st: &mut u64| -> u8 {
            *st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            nts[((*st >> 33) as usize) % 4]
        };
        let s1: Vec<u8> = (0..230).map(|_| next_nt(&mut state)).collect();
        let s2: Vec<u8> = (0..240).map(|_| next_nt(&mut state)).collect();
        compare_alignments(&s1, &s2, "diff-len-230-vs-240");
    }

    /// Randomized parity stress test for the diag-fill band-corner range
    /// arithmetic: `align_vectorized` must match `align_endsfree`'s optimal
    /// score across many random pairs (varied lengths, indels, bands). Run
    /// explicitly (it does ~30k alignments):
    ///   cargo test --release -- --ignored sweep_vectorized_parity --nocapture
    #[test]
    #[ignore]
    fn sweep_vectorized_parity() {
        let nts = [1u8, 2, 3, 4];
        let mut st: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rng = |st: &mut u64, m: usize| {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*st >> 33) as usize) % m
        };
        let mut fails = 0;
        for _ in 0..10_000 {
            let len = 5 + rng(&mut st, 60);
            let s1: Vec<u8> = (0..len).map(|_| nts[rng(&mut st, 4)]).collect();
            let mut s2 = s1.clone();
            for _ in 0..rng(&mut st, 5) {
                if s2.is_empty() {
                    break;
                }
                let p = rng(&mut st, s2.len());
                match rng(&mut st, 3) {
                    0 => s2[p] = nts[rng(&mut st, 4)],       // substitution
                    1 => s2.insert(p, nts[rng(&mut st, 4)]), // insertion
                    _ => {
                        s2.remove(p);
                    } // deletion
                }
            }
            if s2.len() < 3 {
                continue;
            }
            for &band in &[8i32, 16, 32] {
                let ef = align_endsfree(&s1, &s2, 5, -4, -8, band);
                let ve = align_vectorized(
                    &s1,
                    &s2,
                    &VectorizedAlignScores {
                        match_score: 5,
                        mismatch: -4,
                        gap_p: -8,
                        end_gap_p: 0,
                        band,
                    },
                );
                if score_alignment(&ef, 5, -4, -8) != score_alignment(&ve, 5, -4, -8) {
                    fails += 1;
                }
            }
        }
        assert_eq!(fails, 0, "vectorized/endsfree score parity must hold");
    }

    /// A distinct homopolymer gap penalty must route `raw_align` through the
    /// homopolymer-aware scalar DP, not the vectorized aligner (which has no
    /// homopolymer-gap support). Mirrors R DADA2 forcing VECTORIZED_ALIGNMENT
    /// off when HOMOPOLYMER_GAP_PENALTY != GAP_PENALTY (dada.R:229-230).
    #[test]
    fn homo_gap_penalty_routes_to_homopolymer_aligner() {
        // Unequal lengths force a gapped alignment (the gapless fast-path
        // requires equal lengths). This pair is chosen so the homopolymer-aware
        // and vectorized aligners place the gap *differently*: the cheap
        // homopolymer gap (-1) is preferred inside the A-run, while the uniform
        // gap (-8) vectorized path resolves it elsewhere — so the test
        // genuinely distinguishes which aligner ran.
        let s1 = encode("AGAAAAGGGGTTTAAAAAATTTTTCCCC");
        let s2 = encode("AAAAAGGGGTTTAAAAAATTTTTCCCC");

        let params = AlignParams {
            match_score: 5,
            mismatch: -4,
            gap_p: -8,
            homo_gap_p: -1,
            use_kmers: true,
            kdist_cutoff: 0.42,
            kmer_size: 5,
            band: 16,
            vectorized: true, // must be overridden by the homopolymer branch
            gapless: true,
        };

        // Reference: the homopolymer-aware aligner called directly.
        let mut rbuf = AlignBuffers::new();
        align_endsfree_homo_with_buf(&s1, &s2, &params, &mut rbuf);
        let want = (rbuf.al0.clone(), rbuf.al1.clone());

        // Sanity: the vectorized aligner gives a *different* alignment here, so
        // this test genuinely distinguishes the dispatch (the pre-fix path used
        // vectorized and would silently drop the homopolymer penalty).
        let vec = align_vectorized(
            &s1,
            &s2,
            &VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                end_gap_p: 0,
                band: 16,
            },
        );
        assert_ne!(
            (vec[0].clone(), vec[1].clone()),
            want,
            "test precondition: vectorized and homopolymer alignments must differ"
        );

        // Dispatched path with vectorized=true AND a homopolymer penalty set.
        let mut r1 = Raw::new(s1, None, 10, false);
        let mut r2 = Raw::new(s2, None, 5, false);
        crate::kmers::raw_assign_kmers(&mut r1, params.kmer_size);
        crate::kmers::raw_assign_kmers(&mut r2, params.kmer_size);

        let mut buf = AlignBuffers::new();
        raw_align_with_buf(&r1, &r2, &params, &mut buf).expect("alignment produced");
        assert_eq!(
            (buf.al0.clone(), buf.al1.clone()),
            want,
            "homopolymer gap penalty must use the homopolymer-aware aligner, not vectorized"
        );
    }
}
