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
use crate::kmers::{kmer_dist, kmer_dist8, kord_dist, KMER_SIZE};

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
    /// Band radius. Negative means unbanded.
    pub band: i32,
    pub vectorized: bool,
    pub gapless: bool,
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Return true if byte `b` encodes a nucleotide (A/C/G/T/N = 1..=5).
#[inline]
fn is_nt(b: u8) -> bool {
    matches!(b, 1..=5)
}

/// Compute per-position homopolymer mask: true at positions inside a run of
/// three or more identical nucleotides.  Equivalent to C++ `homo1`/`homo2`.
fn homopolymer_mask(seq: &[u8]) -> Vec<bool> {
    let n = seq.len();
    let mut mask = vec![false; n];
    let mut run_start = 0;
    while run_start < n {
        let mut run_end = run_start + 1;
        while run_end < n && seq[run_end] == seq[run_start] {
            run_end += 1;
        }
        if run_end - run_start >= 3 {
            for k in run_start..run_end {
                mask[k] = true;
            }
        }
        run_start = run_end;
    }
    mask
}

/// Trace back through the pointer matrix `p` and build the alignment pair.
/// Shared by `align_endsfree`, `align_endsfree_homo`, and `align_standard`.
fn traceback(p: &[u8], ncol: usize, s1: &[u8], s2: &[u8], len1: usize, len2: usize) -> [Vec<u8>; 2] {
    let mut al0 = Vec::with_capacity(len1 + len2);
    let mut al1 = Vec::with_capacity(len1 + len2);
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
    [al0, al1]
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
fn fill_band_sentinels(d: &mut [i32], ncol: usize, len1: usize, len2: usize, lband: i32, rband: i32, band: i32) {
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
pub fn align_endsfree(s1: &[u8], s2: &[u8], match_score: i32, mismatch: i32, gap_p: i32, band: i32) -> [Vec<u8>; 2] {
    let len1 = s1.len();
    let len2 = s2.len();
    let ncol = len2 + 1;
    let nrow = len1 + 1;

    let mut d = vec![0i32; nrow * ncol];
    let mut p = vec![0u8; nrow * ncol];

    // Initialise edges (ends-free: score 0).
    for i in 0..=len1 { p[i * ncol] = 3; }
    for j in 0..=len2 { p[j] = 2; }

    let (lband, rband) = band_adjust(len1, len2, band);
    fill_band_sentinels(&mut d, ncol, len1, len2, lband, rband, band);

    for i in 1..=len1 {
        let l = if band >= 0 { (i as i32 - lband).max(1) as usize } else { 1 };
        let r = if band >= 0 { (i as i32 + rband).min(len2 as i32) as usize } else { len2 };
        for j in l..=r {
            let left = if i == len1 { d[i * ncol + j - 1] } else { d[i * ncol + j - 1] + gap_p };
            let up   = if j == len2 { d[(i - 1) * ncol + j] } else { d[(i - 1) * ncol + j] + gap_p };
            let diag = d[(i - 1) * ncol + j - 1]
                + if s1[i - 1] == s2[j - 1] { match_score } else { mismatch };

            if up >= diag && up >= left {
                d[i * ncol + j] = up;   p[i * ncol + j] = 3;
            } else if left >= diag {
                d[i * ncol + j] = left; p[i * ncol + j] = 2;
            } else {
                d[i * ncol + j] = diag; p[i * ncol + j] = 1;
            }
        }
    }
    traceback(&p, ncol, s1, s2, len1, len2)
}

// ---------------------------------------------------------------------------
// Homopolymer-aware end-gap-free NW
// ---------------------------------------------------------------------------

/// Banded end-gap-free NW with position-specific homopolymer gap penalties.
///
/// Gaps inside homopolymer runs (length ≥ 3) use `homo_gap_p` instead of
/// `gap_p`.  Equivalent to C++ `nwalign_endsfree_homo`.
pub fn align_endsfree_homo(
    s1: &[u8], s2: &[u8],
    match_score: i32, mismatch: i32,
    gap_p: i32, homo_gap_p: i32,
    band: i32,
) -> [Vec<u8>; 2] {
    let len1 = s1.len();
    let len2 = s2.len();
    let ncol = len2 + 1;
    let nrow = len1 + 1;

    let homo1 = homopolymer_mask(s1);
    let homo2 = homopolymer_mask(s2);

    let mut d = vec![0i32; nrow * ncol];
    let mut p = vec![0u8; nrow * ncol];

    for i in 0..=len1 { p[i * ncol] = 3; }
    for j in 0..=len2 { p[j] = 2; }

    let (lband, rband) = band_adjust(len1, len2, band);
    fill_band_sentinels(&mut d, ncol, len1, len2, lband, rband, band);

    for i in 1..=len1 {
        let l = if band >= 0 { (i as i32 - lband).max(1) as usize } else { 1 };
        let r = if band >= 0 { (i as i32 + rband).min(len2 as i32) as usize } else { len2 };
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
                + if s1[i - 1] == s2[j - 1] { match_score } else { mismatch };

            if up >= diag && up >= left {
                d[i * ncol + j] = up;   p[i * ncol + j] = 3;
            } else if left >= diag {
                d[i * ncol + j] = left; p[i * ncol + j] = 2;
            } else {
                d[i * ncol + j] = diag; p[i * ncol + j] = 1;
            }
        }
    }
    traceback(&p, ncol, s1, s2, len1, len2)
}

// ---------------------------------------------------------------------------
// Standard (non-ends-free) NW — not used in core DADA2, included for parity
// ---------------------------------------------------------------------------

/// Standard banded Needleman-Wunsch (edge gaps are penalised).
/// Not used in the core DADA2 algorithm.  Equivalent to C++ `nwalign`.
#[allow(dead_code)]
pub fn align_standard(s1: &[u8], s2: &[u8], match_score: i32, mismatch: i32, gap_p: i32, band: i32) -> [Vec<u8>; 2] {
    let len1 = s1.len();
    let len2 = s2.len();
    let ncol = len2 + 1;
    let nrow = len1 + 1;

    let mut d = vec![0i32; nrow * ncol];
    let mut p = vec![0u8; nrow * ncol];

    for i in 1..=len1 { d[i * ncol] = d[(i - 1) * ncol] + gap_p; p[i * ncol] = 3; }
    for j in 1..=len2 { d[j] = d[j - 1] + gap_p; p[j] = 2; }

    let (lband, rband) = band_adjust(len1, len2, band);
    fill_band_sentinels(&mut d, ncol, len1, len2, lband, rband, band);

    for i in 1..=len1 {
        let l = if band >= 0 { (i as i32 - lband).max(1) as usize } else { 1 };
        let r = if band >= 0 { (i as i32 + rband).min(len2 as i32) as usize } else { len2 };
        for j in l..=r {
            let left = d[i * ncol + j - 1] + gap_p;
            let up   = d[(i - 1) * ncol + j] + gap_p;
            let diag = d[(i - 1) * ncol + j - 1]
                + if s1[i - 1] == s2[j - 1] { match_score } else { mismatch };

            if up >= diag && up >= left {
                d[i * ncol + j] = up;   p[i * ncol + j] = 3;
            } else if left >= diag {
                d[i * ncol + j] = left; p[i * ncol + j] = 2;
            } else {
                d[i * ncol + j] = diag; p[i * ncol + j] = 1;
            }
        }
    }
    traceback(&p, ncol, s1, s2, len1, len2)
}

// ---------------------------------------------------------------------------
// Gapless alignment
// ---------------------------------------------------------------------------

/// Position-by-position alignment without gaps.
/// Shorter sequence is padded with gaps on the right.
/// Equivalent to C++ `nwalign_gapless`.
pub fn align_gapless(s1: &[u8], s2: &[u8]) -> [Vec<u8>; 2] {
    let len = s1.len().max(s2.len());
    let al0 = (0..len).map(|i| if i < s1.len() { s1[i] } else { b'-' }).collect();
    let al1 = (0..len).map(|i| if i < s2.len() { s2[i] } else { b'-' }).collect();
    [al0, al1]
}

// ---------------------------------------------------------------------------
// Vectorized (diagonal-banded) NW  — port of nwalign_vectorized2
// ---------------------------------------------------------------------------

/// DP inner loop with `up ≥ left ≥ diag` tie-breaking precedence.
/// Equivalent to C++ `dploop_vec`.
#[inline]
fn dploop(
    d: &mut [i16], p: &mut [i16],
    d_prev: &[i16],
    diag_buf: &[i16],
    left_off: usize, up_off: usize, out_off: usize,
    col_min: usize, n: usize,
    gap_p: i16, swap: bool,
) {
    for k in 0..n {
        let left  = d_prev[left_off + k] + gap_p;
        let diag  = diag_buf[col_min + k];
        let up    = d_prev[up_off + k] + gap_p;

        let (entry, pentry) = if swap {
            // left ≥ up ≥ diag
            if left >= up { (left, 2i16) } else { (up, 3i16) }
        } else {
            // up ≥ left ≥ diag
            if up >= left { (up, 3i16) } else { (left, 2i16) }
        };
        let (entry, pentry) = if entry >= diag { (entry, pentry) } else { (diag, 1i16) };

        d[out_off + k] = entry;
        p[out_off + k] = pentry;
    }
}

/// Cache-friendly diagonal-banded Needleman-Wunsch using `i16` DP tables.
///
/// Processes anti-diagonals (i+j = constant) left to right, giving sequential
/// memory access patterns that auto-vectorise with LLVM.  `end_gap_p = 0`
/// gives ends-free behaviour; `end_gap_p = gap_p` gives standard NW edge costs.
/// Equivalent to C++ `nwalign_vectorized2`.
pub fn align_vectorized(
    s1_in: &[u8], s2_in: &[u8],
    match_score: i16, mismatch: i16,
    gap_p: i16, end_gap_p: i16,
    band: i32,
) -> [Vec<u8>; 2] {
    // Ensure s1 is the shorter sequence; record whether we swapped.
    let swap = s1_in.len() > s2_in.len();
    let (s1, s2) = if swap { (s2_in, s1_in) } else { (s1_in, s2_in) };
    let len1 = s1.len();
    let len2 = s2.len();

    let band = if band < 0 { len2 as i32 } else { band };
    let band_usize = band as usize;

    // Compressed matrix dimensions (diagonal layout).
    // Column index for original cell (i,j): (2*start_col + j - i) / 2
    let start_col = 1 + (1 + band_usize.min(len1)) / 2;
    let ncol = 2 + start_col + (band_usize + len2 - len1).min(len2) / 2;
    let nrow = len1 + len2 + 1;

    let mut d = vec![0i16; ncol * nrow];
    let mut p = vec![0i16; ncol * nrow];
    let mut diag_buf = vec![0i16; ncol];

    // Sentinel fill: columns 0,1 and ncol-2,ncol-1 in every row act as hard
    // band boundaries.  fill_val is chosen so fill_val + gap_p doesn't overflow.
    let min_score = mismatch.min(gap_p).min(match_score).min(0);
    let fill_val = i16::MIN.wrapping_sub(min_score);
    for row in 0..nrow {
        d[row * ncol]         = fill_val;
        d[row * ncol + 1]     = fill_val;
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
            if row % 2 == 0 { col = col.saturating_sub(1); }
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
            if row % 2 == 1 { col += 1; }
            row += 1;
            d_free = d_free.saturating_add(end_gap_p);
        }
    }

    // Main DP: iterate over anti-diagonals (row = i + j).
    let mut row       = 2usize;
    let mut col_min   = start_col;
    let mut col_max   = start_col;
    let mut i_max     = 0i64; // 0-indexed into s1 (decrements along anti-diag)
    let mut j_min     = 0usize; // 0-indexed into s2 (increments along anti-diag)
    let mut even      = true;
    let mut recalc_left  = false;
    let mut recalc_right = false;

    while row <= len1 + len2 {
        let n = col_max - col_min + 1;

        // --- Fill diag_buf for this anti-diagonal ---
        // Cell (i,j) in the original NW matrix uses s1[i-1] vs s2[j-1].
        // Here i_max / j_min are 0-indexed, so s1[i_max] vs s2[j_min].
        {
            let mut si = i_max;
            let mut sj = j_min;
            for k in 0..n {
                debug_assert!(si >= 0 && (si as usize) < len1,
                    "diag_buf si={si} out of range (len1={len1})");
                debug_assert!(sj < len2,
                    "diag_buf sj={sj} out of range (len2={len2})");
                let score = if s1[si as usize] == s2[sj] { match_score } else { mismatch };
                diag_buf[col_min + k] = d[(row - 2) * ncol + col_min + k] + score;
                si -= 1;
                sj += 1;
            }
        }

        // --- Compute d/p for this row using the previous row ---
        // left  = d[(row-1)*ncol + col_min - even]
        // up    = d[(row-1)*ncol + col_min + 1 - even]
        // out   = d[row*ncol + col_min]
        let even_off = if even { 1 } else { 0 }; // even=true → subtract 1
        let prev_base = (row - 1) * ncol;
        let left_off  = prev_base + col_min - even_off;
        let up_off    = prev_base + col_min + 1 - even_off;
        let _out_off  = row * ncol + col_min;

        // Split d into previous-row (read) and current-row (write) slices.
        let (d_prev_part, d_cur_part) = d.split_at_mut(row * ncol);
        let d_prev = &d_prev_part[prev_base..];
        let d_cur  = &mut d_cur_part[..ncol];

        let (p_prev_part, p_cur_part) = p.split_at_mut(row * ncol);
        let _ = p_prev_part; // unused; p_cur accessed via index below
        let p_cur = &mut p_cur_part[..ncol];

        dploop(
            d_cur, p_cur,
            d_prev,
            &diag_buf,
            left_off - prev_base,   // relative to d_prev
            up_off - prev_base,
            col_min,                // relative to d_cur / p_cur
            col_min,
            n,
            gap_p, swap,
        );

        // --- Band transition: widen active columns at wedge boundaries ---
        if row == band_usize.min(len1) {
            col_min = col_min.saturating_sub(1);
            i_max += 1;
            if j_min > 0 { j_min -= 1; }
        }
        if row == (band_usize + len2 - len1).min(len2) {
            col_max += 1;
        }

        // --- End-gap recalculation (when end_gap_p > gap_p, e.g. ends-free) ---
        // left_off and up_off are absolute indices into d (= prev_base + relative).
        if end_gap_p > gap_p {
            // Left boundary: gap in s2 extending to end of s1
            if recalc_left {
                let d_free = d[left_off] + end_gap_p;
                let cur    = d[row * ncol + col_min];
                let pcur   = p[row * ncol + col_min];
                if d_free > cur {
                    d[row * ncol + col_min] = d_free;
                    p[row * ncol + col_min] = 2;
                } else if !swap && d_free == cur && pcur == 1 {
                    p[row * ncol + col_min] = 2;
                } else if swap && d_free == cur && pcur != 2 {
                    p[row * ncol + col_min] = 2;
                }
            }
            if i_max == len1 as i64 - 1 { recalc_left = true; }

            // Right boundary: gap in s1 extending to end of s2
            if recalc_right {
                let d_free = d[up_off + col_max - col_min] + end_gap_p;
                let cur    = d[row * ncol + col_max];
                let pcur   = p[row * ncol + col_max];
                if d_free > cur {
                    d[row * ncol + col_max] = d_free;
                    p[row * ncol + col_max] = 3;
                } else if !swap && d_free == cur && pcur != 3 {
                    p[row * ncol + col_max] = 3;
                } else if swap && d_free == cur && pcur == 1 {
                    p[row * ncol + col_max] = 3;
                }
            }
            let j_max_1idx = (row + 1) / 2 + col_max - start_col;
            if j_max_1idx == len2 { recalc_right = true; }
        }

        // --- Update col_min, col_max, i_max, j_min for next anti-diagonal ---
        let band_mod2 = band % 2;
        if (row as i32) < band && (row as i32) < len1 as i32 {
            // Upper triangle for s1
            if even { col_min = col_min.saturating_sub(1); }
            i_max += 1;
        } else if i_max < (len1 as i64) - 1 {
            // Banded area
            if band_mod2 == 0 {
                if even  { j_min += 1; }
                else     { i_max += 1; }
            } else {
                if even  { col_min = col_min.saturating_sub(1); i_max += 1; }
                else     { col_min += 1; j_min += 1; }
            }
        } else {
            // Lower triangle for s1
            if !even { col_min += 1; }
            j_min += 1;
        }

        let top_limit = (band_usize + len2 - len1).min(len2);
        if row < top_limit {
            if !even { col_max += 1; }
        } else if (row + 1) / 2 + col_max - start_col < len2 {
            let full_band = band_usize + len2 - len1;
            if full_band % 2 == 0 {
                if even { col_max = col_max.saturating_sub(1); }
                else    { col_max += 1; }
            }
            // no action for odd full_band
        } else {
            if even { col_max = col_max.saturating_sub(1); }
        }

        row += 1;
        even = !even;
    }

    // --- Traceback through compressed p matrix ---
    let mut al0 = Vec::with_capacity(len1 + len2);
    let mut al1 = Vec::with_capacity(len1 + len2);
    let mut i = len1;
    let mut j = len2;

    while i > 0 || j > 0 {
        // Compressed column: (2*start_col + j - i) / 2
        let col_signed = 2 * start_col as i64 + j as i64 - i as i64;
        debug_assert!(col_signed >= 0 && col_signed % 2 == 0,
            "vectorized traceback: bad col_signed={col_signed} at i={i} j={j}");
        let col = (col_signed / 2) as usize;
        match p[(i + j) * ncol + col] {
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

    // Restore original sequence ordering if we swapped.
    if swap { [al1, al0] } else { [al0, al1] }
}

// ---------------------------------------------------------------------------
// Substitution compression
// ---------------------------------------------------------------------------

/// Convert an alignment pair into a `Sub` (compressed substitution record).
///
/// Records substitutions of `al[1]` relative to `al[0]`, ignoring positions
/// where either strand has an N (encoded 5).  `Sub::q0`/`Sub::q1` are left
/// empty; fill them via `sub_new` if quality scores are needed.
/// Equivalent to C++ `al2subs`.
pub fn al2subs(al: &[Vec<u8>; 2]) -> Sub {
    let al0 = &al[0];
    let al1 = &al[1];
    let alen = al0.len();

    // First pass: count reference length and substitution count.
    let mut len0 = 0u32;
    let mut nsubs = 0usize;
    for i in 0..alen {
        let nt0 = is_nt(al0[i]);
        let nt1 = is_nt(al1[i]);
        if nt0 { len0 += 1; }
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
        if nt0 { i0 += 1; }
        if nt1 { i1 += 1; }

        if nt0 {
            map[i0 as usize] = if nt1 { i1 as u16 } else { GAP_GLYPH };
        }
        if nt0 && nt1 && al0[i] != al1[i] && al0[i] != 5 && al1[i] != 5 {
            pos.push(i0 as u16);
            nt0_vec.push(al0[i]);
            nt1_vec.push(al1[i]);
        }
    }

    Sub { len0, map, pos, nt0: nt0_vec, nt1: nt1_vec, q0: Vec::new(), q1: Vec::new() }
}

// ---------------------------------------------------------------------------
// Dispatcher and sub_new
// ---------------------------------------------------------------------------

/// Select and run the appropriate alignment for two `Raw` objects.
///
/// Returns `None` if k-mer screening determines the sequences are too
/// dissimilar to be worth aligning (i.e. they will produce a NULL Sub).
/// Equivalent to C++ `raw_align`.
pub fn raw_align(raw1: &Raw, raw2: &Raw, p: &AlignParams) -> Option<[Vec<u8>; 2]> {
    // --- K-mer screening ---
    let mut kdist = 0.0f64;
    let mut kodist = -1.0f64; // sentinel: different from kdist when use_kmers=false

    if p.use_kmers {
        // Prefer 8-bit kmer distance; fall back to 16-bit on overflow.
        kdist = match (&raw1.kmer8, &raw2.kmer8) {
            (Some(k1), Some(k2)) => {
                let d8 = kmer_dist8(k1, raw1.len(), k2, raw2.len(), KMER_SIZE);
                if d8 < 0.0 {
                    // Overflow: use 16-bit vectors.
                    match (&raw1.kmer, &raw2.kmer) {
                        (Some(k1), Some(k2)) => kmer_dist(k1, raw1.len(), k2, raw2.len(), KMER_SIZE),
                        _ => 0.0,
                    }
                } else { d8 }
            }
            _ => match (&raw1.kmer, &raw2.kmer) {
                (Some(k1), Some(k2)) => kmer_dist(k1, raw1.len(), k2, raw2.len(), KMER_SIZE),
                _ => 0.0,
            },
        };

        if p.gapless {
            if let (Some(o1), Some(o2)) = (&raw1.kord, &raw2.kord) {
                kodist = kord_dist(o1, raw1.len(), o2, raw2.len(), KMER_SIZE);
            }
        }
    }

    if p.use_kmers && kdist > p.kdist_cutoff {
        return None; // Outside k-mer distance threshold → NULL alignment.
    }

    // --- Method selection ---
    if p.band == 0 || (p.gapless && (kodist - kdist).abs() < f64::EPSILON) {
        return Some(align_gapless(&raw1.seq, &raw2.seq));
    }
    if p.vectorized {
        return Some(align_vectorized(
            &raw1.seq, &raw2.seq,
            p.match_score as i16, p.mismatch as i16, p.gap_p as i16,
            0, // end_gap_p = 0 for ends-free
            p.band,
        ));
    }
    if p.homo_gap_p != p.gap_p && p.homo_gap_p <= 0 {
        return Some(align_endsfree_homo(
            &raw1.seq, &raw2.seq,
            p.match_score, p.mismatch, p.gap_p, p.homo_gap_p, p.band,
        ));
    }
    Some(align_endsfree(&raw1.seq, &raw2.seq, p.match_score, p.mismatch, p.gap_p, p.band))
}

/// Align two `Raw` objects and return the compressed substitution record,
/// with quality scores filled in when both Raws carry quality data.
///
/// Returns `None` when the k-mer screen rejects the pair (equivalent to a
/// NULL Sub in the C++ code).
/// Equivalent to C++ `sub_new`.
pub fn sub_new(raw0: &Raw, raw1: &Raw, params: &AlignParams) -> Option<Sub> {
    let al = raw_align(raw0, raw1, params)?;
    let mut sub = al2subs(&al);

    if let (Some(q0), Some(q1)) = (&raw0.qual, &raw1.qual) {
        sub.q0 = sub.pos.iter().map(|&pos| q0[pos as usize]).collect();
        sub.q1 = sub.pos.iter()
            .map(|&pos| q1[sub.map[pos as usize] as usize])
            .collect();
    }
    Some(sub)
}
