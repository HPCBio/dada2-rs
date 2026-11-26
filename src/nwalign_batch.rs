//! Inter-sequence batched banded NW (experimental, `--batch-align`).
//!
//! [`align_vectorized_batch`] runs the algorithm of
//! [`crate::nwalign::align_vectorized_with_buf`] (R `nwalign_vectorized2`) on
//! `L` pairs at once. The band walk (`col_min`, `col_max`, `i_max`, `j_min`,
//! edge fills) depends only on `(len1, len2, band)`, so pairs of identical
//! shape run it in lockstep; DP storage is lane-interleaved (`cell * L + lane`)
//! and every per-cell operation becomes `L` contiguous lanes. Output per lane
//! is byte-identical to the single-pair kernel (see the tests).
//!
//! Used by `cluster::b_compare_parallel` under `--batch-align`;
//! `examples/batch_nw_bench.rs` measures it against the single-pair kernel.

use crate::containers::Sub;
use crate::nwalign::{GAP_GLYPH, VectorizedAlignScores, align_vectorized_with_buf, reset_pointers};

/// What the batched kernel produces after the fill.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Out {
    /// Nothing: the fill alone, for timing.
    Fill,
    /// The two aligned strings per lane, in `BatchBuffers::al`.
    Align,
    /// The substitution record per lane, in `BatchBuffers::subs` (#342).
    Sub,
}

/// Reusable scratch for [`align_vectorized_batch`].
#[derive(Default)]
pub struct BatchBuffers {
    d: Vec<i16>,
    p: Vec<u8>,
    diag: Vec<i16>,
    s1: Vec<u8>,
    s2: Vec<u8>,
    /// Alignment of each lane, in the order of its input pair.
    pub al: Vec<(Vec<u8>, Vec<u8>)>,
    /// [`align_vectorized_batch_subs`]: each lane's substitutions of its
    /// second sequence against its first, as `al2subs` on `al` gives them,
    /// without the quality fields.
    pub subs: Vec<Sub>,
    /// Time the traceback into `last_traceback_nanos` (metrics attribution).
    pub measure: bool,
    /// Nanos in the last batch's traceback, all lanes, when `measure`.
    pub last_traceback_nanos: u64,
}

impl BatchBuffers {
    pub fn new() -> Self {
        Self::default()
    }
}

#[inline]
fn reset_buf<T: Copy>(v: &mut Vec<T>, n: usize, val: T) {
    if v.len() < n {
        v.clear();
        v.resize(n, val);
    } else {
        v[..n].fill(val);
    }
}

/// `nwalign::dploop` over `n * L` interleaved cells.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn dploop(
    d: &mut [i16],
    p: &mut [u8],
    d_prev: &[i16],
    diag: &[i16],
    left_off: usize,
    up_off: usize,
    diag_off: usize,
    m: usize,
    gap_p: i16,
    swap: bool,
) {
    let left = &d_prev[left_off..left_off + m];
    let up = &d_prev[up_off..up_off + m];
    let diag = &diag[diag_off..diag_off + m];
    let d_out = &mut d[..m];
    let p_out = &mut p[..m];
    if swap {
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

/// Align `pairs.len()` pairs, all with the same `(a.len(), b.len())`, `L` at
/// a time, passing each alignment to `sink` in input order. A short final
/// group is padded by repeating its first pair. Lane buffers are reused, so
/// the steady state allocates nothing; moving them out per batch is the
/// suspected cause of the bench collapsing at 16+ threads on the cluster
/// (`L16alloc` arm in `examples/batch_nw_bench.rs` is the control).
pub fn align_batch_all<const L: usize>(
    pairs: &[(&[u8], &[u8])],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
    mut sink: impl FnMut(&[u8], &[u8]),
) {
    for chunk in pairs.chunks(L) {
        let mut lanes = [chunk[0]; L];
        lanes[..chunk.len()].copy_from_slice(chunk);
        align_vectorized_batch::<L>(&lanes, scores, buf);
        for (al0, al1) in &buf.al[..chunk.len()] {
            sink(al0, al1);
        }
    }
}

/// Batched [`align_vectorized_with_buf`]: `L` pairs of identical shape.
///
/// # Panics
/// If the pairs do not all share `(a.len(), b.len())`.
pub fn align_vectorized_batch<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
) {
    batch_impl::<L>(pairs, scores, buf, Out::Align);
}

/// [`align_vectorized_batch`], but the traceback writes each lane's [`Sub`]
/// (reference = the pair's first sequence) into `buf.subs` instead of the
/// aligned strings, so they need neither building nor a second pass in
/// `al2subs` (#342). `buf.al` is left untouched. The `q0`/`q1` fields are
/// empty; the caller adds qualities.
pub fn align_vectorized_batch_subs<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
) {
    batch_impl::<L>(pairs, scores, buf, Out::Sub);
}

/// The DP fill of [`align_vectorized_batch`] without the traceback, for timing
/// the two apart. Leaves `buf.al` untouched.
pub fn align_vectorized_batch_fill_only<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
) {
    batch_impl::<L>(pairs, scores, buf, Out::Fill);
}

/// Runs the AVX2 build when [`crate::simd::avx2`] allows it (issue #336); the
/// output is identical either way.
fn batch_impl<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
    out: Out,
) {
    #[cfg(target_arch = "x86_64")]
    if crate::simd::avx2() {
        // SAFETY: `simd::avx2` is true only when the CPU has AVX2.
        unsafe { batch_avx2::<L>(pairs, scores, buf, out) };
        return;
    }
    batch_body::<L>(pairs, scores, buf, out);
}

/// [`align_vectorized_batch`] compiled for the baseline target, whatever the
/// host supports. For the variant tests and benches.
pub fn align_vectorized_batch_baseline<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
) {
    batch_body::<L>(pairs, scores, buf, Out::Align);
}

/// [`align_vectorized_batch`] compiled with AVX2.
///
/// # Safety
/// The CPU must support AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn align_vectorized_batch_avx2<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
) {
    batch_body::<L>(pairs, scores, buf, Out::Align);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn batch_avx2<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
    out: Out,
) {
    batch_body::<L>(pairs, scores, buf, out);
}

#[inline(always)]
fn batch_body<const L: usize>(
    pairs: &[(&[u8], &[u8]); L],
    scores: &VectorizedAlignScores,
    buf: &mut BatchBuffers,
    out: Out,
) {
    let VectorizedAlignScores {
        match_score,
        mismatch,
        gap_p,
        begin_gap_p,
        end_gap_p,
        band,
    } = *scores;
    let (la, lb) = (pairs[0].0.len(), pairs[0].1.len());
    assert!(
        pairs.iter().all(|(a, b)| a.len() == la && b.len() == lb),
        "align_vectorized_batch: pairs must share one shape"
    );
    let swap = la > lb;
    let (len1, len2) = if swap { (lb, la) } else { (la, lb) };

    let band = if band < 0 { len2 as i32 } else { band };
    let band_usize = band as usize;
    let start_col = 1 + band_usize.min(len1).div_ceil(2);
    let ncol = 2 + start_col + (band_usize + len2 - len1).min(len2) / 2;
    let nrow = len1 + len2 + 1;

    // Interleaved sequences: s[pos * L + lane], s1 reversed.
    buf.s1.clear();
    buf.s1.resize(len1 * L, 0);
    buf.s2.clear();
    buf.s2.resize(len2 * L, 0);
    for (lane, &(a, b)) in pairs.iter().enumerate() {
        let (s1, s2) = if swap { (b, a) } else { (a, b) };
        for (pos, &c) in s1.iter().rev().enumerate() {
            buf.s1[pos * L + lane] = c;
        }
        for (pos, &c) in s2.iter().enumerate() {
            buf.s2[pos * L + lane] = c;
        }
    }
    let s1 = &buf.s1[..];
    let s2 = &buf.s2[..];

    const DROWS: usize = 3;
    let w = ncol * L;
    reset_buf(&mut buf.d, w * DROWS, 0i16);
    reset_pointers(&mut buf.p, w * nrow);
    reset_buf(&mut buf.diag, w, 0i16);
    let p = &mut buf.p[..w * nrow];
    let diag_buf = &mut buf.diag[..w];

    let min_score = mismatch.min(gap_p).min(match_score).min(0);
    let fill_val = i16::MIN.wrapping_sub(min_score);

    let (d_a, d_rest) = buf.d[..w * DROWS].split_at_mut(w);
    let (d_b, d_c) = d_rest.split_at_mut(w);
    let mut slots: [&mut [i16]; DROWS] = [d_a, d_b, d_c];

    // Cell `c` of a row, all lanes.
    let lanes = |c: usize| c * L..c * L + L;

    let left_fill_limit = 1 + band_usize.min(len1);
    let top_fill_limit = 1 + (band_usize + len2 - len1).min(len2);
    let left_col = |row: usize| (start_col - 1).saturating_sub((row - 1) / 2);
    let top_col = |row: usize| start_col + row / 2;

    for (row, slot) in slots.iter_mut().enumerate().take(2) {
        for c in [0, 1, ncol - 2, ncol - 1] {
            slot[lanes(c)].fill(fill_val);
        }
        if row == 0 {
            slot[lanes(start_col)].fill(0);
        } else {
            if row < left_fill_limit {
                slot[lanes(left_col(row))].fill(begin_gap_p);
                p[row * w..][lanes(left_col(row))].fill(3);
            }
            if row < top_fill_limit {
                slot[lanes(top_col(row))].fill(begin_gap_p);
                p[row * w..][lanes(top_col(row))].fill(2);
            }
        }
    }
    p[lanes(start_col)].fill(0);
    let mut row_free = begin_gap_p;

    let mut row = 2usize;
    let mut col_min = start_col;
    let mut col_max = start_col;
    let mut i_max = 0i64;
    let mut j_min = 0usize;
    let mut even = true;
    let mut recalc_left = false;
    let mut recalc_right = false;

    while row <= len1 + len2 {
        let n = col_max - col_min + 1;

        row_free = row_free.saturating_add(begin_gap_p);
        {
            let cur = &mut slots[2];
            for c in [0, 1, ncol - 2, ncol - 1] {
                cur[lanes(c)].fill(fill_val);
            }
            if row < left_fill_limit {
                cur[lanes(left_col(row))].fill(row_free);
                p[row * w..][lanes(left_col(row))].fill(3);
            }
            if row < top_fill_limit {
                cur[lanes(top_col(row))].fill(row_free);
                p[row * w..][lanes(top_col(row))].fill(2);
            }
        }

        {
            let d_prev2 = &slots[0][col_min * L..(col_min + n) * L];
            let diag_out = &mut diag_buf[col_min * L..(col_min + n) * L];

            let k_lo = (i_max - len1 as i64 + 1).clamp(0, n as i64) as usize;
            let k_hi = (i_max + 1)
                .min(len2 as i64 - j_min as i64)
                .clamp(0, n as i64) as usize;
            let k_hi = k_hi.max(k_lo);

            for (dst, &prev) in diag_out[..k_lo * L].iter_mut().zip(&d_prev2[..k_lo * L]) {
                *dst = prev.saturating_add(fill_val);
            }
            for (dst, &prev) in diag_out[k_hi * L..].iter_mut().zip(&d_prev2[k_hi * L..]) {
                *dst = prev.saturating_add(fill_val);
            }
            // s1 is stored reversed, so both sequences advance with k and the
            // in-band middle is one flat run of (k_hi - k_lo) * L cells.
            if k_lo < k_hi {
                let m = (k_hi - k_lo) * L;
                let r0 = (len1 as i64 - 1 - i_max + k_lo as i64) as usize;
                let sj0 = j_min.wrapping_add(k_lo);
                let a = &s1[r0 * L..r0 * L + m];
                let b = &s2[sj0 * L..sj0 * L + m];
                let prev = &d_prev2[k_lo * L..k_lo * L + m];
                let dst = &mut diag_out[k_lo * L..k_lo * L + m];
                for (((o, &pv), &x), &y) in dst.iter_mut().zip(prev).zip(a).zip(b) {
                    let score = if x == y { match_score } else { mismatch };
                    *o = pv.saturating_add(score);
                }
            }
        }

        let even_off = if even { 1 } else { 0 };
        let left_off = col_min - even_off;
        let up_off = col_min + 1 - even_off;

        {
            let [_, d_prev, d_cur] = &mut slots;
            let base = row * w + col_min * L;
            dploop(
                &mut d_cur[col_min * L..],
                &mut p[base..base + n * L],
                d_prev,
                diag_buf,
                left_off * L,
                up_off * L,
                col_min * L,
                n * L,
                gap_p,
                swap,
            );
        }

        if row == band_usize.min(len1) {
            col_min = col_min.saturating_sub(1);
            i_max += 1;
            j_min = j_min.wrapping_sub(1);
        }
        if row == (band_usize + len2 - len1).min(len2) {
            col_max += 1;
        }

        if end_gap_p > gap_p {
            if recalc_left {
                for lane in 0..L {
                    let d_free = slots[1][left_off * L + lane].saturating_add(end_gap_p);
                    let ci = col_min * L + lane;
                    let cur = slots[2][ci];
                    let pi = row * w + ci;
                    let pcur = p[pi];
                    if d_free > cur {
                        slots[2][ci] = d_free;
                        p[pi] = 2;
                    } else if !(d_free != cur || !swap && pcur != 1 || swap && pcur == 2) {
                        p[pi] = 2;
                    }
                }
            }
            if i_max == len1 as i64 - 1 {
                recalc_left = true;
            }

            if recalc_right {
                let src = up_off + col_max - col_min;
                for lane in 0..L {
                    let d_free = slots[1][src * L + lane].saturating_add(end_gap_p);
                    let ci = col_max * L + lane;
                    let cur = slots[2][ci];
                    let pi = row * w + ci;
                    let pcur = p[pi];
                    if d_free > cur {
                        slots[2][ci] = d_free;
                        p[pi] = 3;
                    } else if !(d_free != cur || !swap && pcur == 3 || swap && pcur != 1) {
                        p[pi] = 3;
                    }
                }
            }
            let j_max_1idx = row.div_ceil(2) + col_max - start_col;
            if j_max_1idx == len2 {
                recalc_right = true;
            }
        }

        let band_mod2 = band % 2;
        if (row as i32) < band && (row as i32) < len1 as i32 {
            if even {
                col_min = col_min.saturating_sub(1);
            }
            i_max += 1;
        } else if i_max < (len1 as i64) - 1 {
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
        } else if even {
            col_max = col_max.saturating_sub(1);
        }

        slots.rotate_left(1);
        row += 1;
        even = !even;
    }

    if out == Out::Fill {
        return;
    }
    // Per-lane traceback.
    let t_tb = buf.measure.then(std::time::Instant::now);
    let p_ro = &buf.p[..w * nrow];
    if out == Out::Sub {
        buf.subs.resize_with(L, Sub::default);
        for (lane, (&(a, b), sub)) in pairs.iter().zip(buf.subs.iter_mut()).enumerate() {
            let ptr = |i: usize, j: usize| {
                let col = ((2 * start_col as i64 + j as i64 - i as i64) / 2) as usize;
                p_ro[((i + j) * ncol + col) * L + lane]
            };
            traceback_sub(ptr, a, b, swap, sub);
        }
        if let Some(t) = t_tb {
            buf.last_traceback_nanos = t.elapsed().as_nanos() as u64;
        }
        return;
    }
    buf.al.resize_with(L, Default::default);
    for (lane, &(a, b)) in pairs.iter().enumerate() {
        let (s1, s2) = if swap { (b, a) } else { (a, b) };
        let (al0, al1) = &mut buf.al[lane];
        al0.clear();
        al1.clear();
        let mut i = len1;
        let mut j = len2;
        while i > 0 || j > 0 {
            let col = ((2 * start_col as i64 + j as i64 - i as i64) / 2) as usize;
            match p_ro[((i + j) * ncol + col) * L + lane] {
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
                v => panic!("batched traceback: invalid pointer {v} at i={i} j={j}"),
            }
        }
        al0.reverse();
        al1.reverse();
        if swap {
            std::mem::swap(al0, al1);
        }
    }
    if let Some(t) = t_tb {
        buf.last_traceback_nanos = t.elapsed().as_nanos() as u64;
    }
}

/// One lane's traceback, written as the [`Sub`] of `b` against `a` that
/// `al2subs` would build from the aligned strings (#342).
///
/// `ptr(i, j)` is the pointer at DP cell `(i, j)` of the kernel's `(s1, s2)`,
/// which is `(b, a)` when `swap`. The walk runs from the end, so each diagonal
/// step knows both positions and bases: it sets `map` and records a
/// substitution where the bases differ and neither is N (5). A gap step leaves
/// the reference position at `GAP_GLYPH`. Substitutions arrive last first and
/// are reversed once at the end. Assumes sequences hold only 1..=5, as raws do,
/// so that `len0` is the reference length.
#[inline(always)]
fn traceback_sub(ptr: impl Fn(usize, usize) -> u8, a: &[u8], b: &[u8], swap: bool, sub: &mut Sub) {
    let (len1, len2) = if swap {
        (b.len(), a.len())
    } else {
        (a.len(), b.len())
    };
    sub.len0 = a.len() as u32;
    sub.map.clear();
    sub.map.resize(a.len(), GAP_GLYPH);
    sub.pos.clear();
    sub.nt0.clear();
    sub.nt1.clear();
    sub.q0.clear();
    sub.q1.clear();
    let (mut i, mut j) = (len1, len2);
    while i > 0 || j > 0 {
        match ptr(i, j) {
            1 => {
                i -= 1;
                j -= 1;
                // Reference (a) and query (b) positions of this column.
                let (ia, ib) = if swap { (j, i) } else { (i, j) };
                let (ca, cb) = (a[ia], b[ib]);
                sub.map[ia] = ib as u16;
                if ca != cb && ca != 5 && cb != 5 {
                    sub.pos.push(ia as u16);
                    sub.nt0.push(ca);
                    sub.nt1.push(cb);
                }
            }
            2 => j -= 1,
            3 => i -= 1,
            v => panic!("batched traceback: invalid pointer {v} at i={i} j={j}"),
        }
    }
    sub.pos.reverse();
    sub.nt0.reverse();
    sub.nt1.reverse();
}

/// Single-pair reference, for the A/B and the identity tests.
pub fn align_single(
    a: &[u8],
    b: &[u8],
    scores: &VectorizedAlignScores,
    buf: &mut crate::nwalign::AlignBuffers,
) -> (Vec<u8>, Vec<u8>) {
    align_vectorized_with_buf(a, b, scores, buf);
    (buf.al0.clone(), buf.al1.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nwalign::AlignBuffers;

    fn rng(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    }

    /// `base` mutated with substitutions and indels, then forced to length `len`.
    fn variant(base: &[u8], len: usize, st: &mut u64) -> Vec<u8> {
        let mut v = Vec::with_capacity(len + 8);
        for &c in base {
            match rng(st) % 40 {
                0 => {}
                1 => {
                    v.push(c);
                    v.push(1 + (rng(st) % 4) as u8);
                }
                2 | 3 => v.push(1 + (rng(st) % 4) as u8),
                _ => v.push(c),
            }
        }
        v.resize(len, 1 + (rng(st) % 4) as u8);
        v
    }

    fn check_identity<const L: usize>(cases: usize, seed: u64) {
        let mut st = seed;
        let mut single = AlignBuffers::new();
        let mut batch = BatchBuffers::new();
        for case in 0..cases {
            let la = 20 + (rng(&mut st) % 300) as usize;
            let lb = la.saturating_sub(12) + (rng(&mut st) % 25) as usize;
            let band = [-1, 1, 2, 3, 7, 8, 16, 17, 32][(rng(&mut st) % 9) as usize];
            let (begin, end) =
                [(0i16, 0i16), (-8, -8), (0, -8), (-8, 0)][(rng(&mut st) % 4) as usize];
            let scores = VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                begin_gap_p: begin,
                end_gap_p: end,
                band,
            };
            let center: Vec<u8> = (0..la).map(|_| 1 + (rng(&mut st) % 4) as u8).collect();
            let raws: Vec<Vec<u8>> = (0..L).map(|_| variant(&center, lb, &mut st)).collect();
            // Alternate which side is the center so both `swap` paths are hit.
            let pairs: [(&[u8], &[u8]); L] = std::array::from_fn(|k| {
                if case % 2 == 0 {
                    (center.as_slice(), raws[k].as_slice())
                } else {
                    (raws[k].as_slice(), center.as_slice())
                }
            });
            align_vectorized_batch::<L>(&pairs, &scores, &mut batch);
            for (lane, &(a, b)) in pairs.iter().enumerate() {
                let want = align_single(a, b, &scores, &mut single);
                assert_eq!(
                    batch.al[lane], want,
                    "case {case} lane {lane}: la={la} lb={lb} band={band} ends=({begin},{end})"
                );
            }
        }
    }

    /// The baseline and AVX2 builds of the batched kernel agree (#336),
    /// including where `i16` scores saturate (`match_score` 40 on long pairs).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn batch_variants_match() {
        if !crate::simd::host_has_avx2("batch_variants_match") {
            return;
        }
        const L: usize = 16;
        let mut st = 0x336u64;
        let (mut base, mut avx) = (BatchBuffers::new(), BatchBuffers::new());
        for case in 0..400 {
            let long = case % 10 == 0;
            let la = if long { 1500 } else { 20 } + (rng(&mut st) % 300) as usize;
            let lb = la.saturating_sub(12) + (rng(&mut st) % 25) as usize;
            let band = [-1, 1, 2, 3, 7, 8, 16, 17, 32][(rng(&mut st) % 9) as usize];
            let (begin, end) =
                [(0i16, 0i16), (-8, -8), (0, -8), (-8, 0)][(rng(&mut st) % 4) as usize];
            let scores = VectorizedAlignScores {
                match_score: if long { 40 } else { 5 },
                mismatch: -4,
                gap_p: -8,
                begin_gap_p: begin,
                end_gap_p: end,
                band,
            };
            let center: Vec<u8> = (0..la).map(|_| 1 + (rng(&mut st) % 4) as u8).collect();
            let raws: Vec<Vec<u8>> = (0..L).map(|_| variant(&center, lb, &mut st)).collect();
            let pairs: [(&[u8], &[u8]); L] = std::array::from_fn(|k| {
                if case % 2 == 0 {
                    (center.as_slice(), raws[k].as_slice())
                } else {
                    (raws[k].as_slice(), center.as_slice())
                }
            });
            align_vectorized_batch_baseline::<L>(&pairs, &scores, &mut base);
            // SAFETY: `host_has_avx2` checked the CPU above.
            unsafe { align_vectorized_batch_avx2::<L>(&pairs, &scores, &mut avx) };
            assert!(
                base.al[..L] == avx.al[..L],
                "case {case}: la={la} lb={lb} band={band} ends=({begin},{end})"
            );
        }
    }

    /// The fused traceback's `Sub` equals `al2subs` on the same lane's aligned
    /// strings (#342): both orientations, every end model, bands from 1 to
    /// unbanded, and Ns (5) in both sequences, which never count as
    /// substitutions.
    #[test]
    fn batch_subs_match_al2subs() {
        const L: usize = 16;
        let mut st = 0x342u64;
        let (mut al, mut fused) = (BatchBuffers::new(), BatchBuffers::new());
        let mut nsubs = 0;
        let mut nn = 0;
        for case in 0..600 {
            let la = 20 + (rng(&mut st) % 300) as usize;
            let lb = la.saturating_sub(12) + (rng(&mut st) % 25) as usize;
            let band = [-1, 1, 2, 3, 7, 8, 16, 17, 32][(rng(&mut st) % 9) as usize];
            let (begin, end) =
                [(0i16, 0i16), (-8, -8), (0, -8), (-8, 0)][(rng(&mut st) % 4) as usize];
            let scores = VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                begin_gap_p: begin,
                end_gap_p: end,
                band,
            };
            let mut center: Vec<u8> = (0..la).map(|_| 1 + (rng(&mut st) % 4) as u8).collect();
            let mut raws: Vec<Vec<u8>> = (0..L).map(|_| variant(&center, lb, &mut st)).collect();
            // Sprinkle Ns into both sides.
            for seq in std::iter::once(&mut center).chain(raws.iter_mut()) {
                for c in seq.iter_mut() {
                    if rng(&mut st).is_multiple_of(50) {
                        *c = 5;
                    }
                }
            }
            let pairs: [(&[u8], &[u8]); L] = std::array::from_fn(|k| {
                if case % 2 == 0 {
                    (center.as_slice(), raws[k].as_slice())
                } else {
                    (raws[k].as_slice(), center.as_slice())
                }
            });
            align_vectorized_batch::<L>(&pairs, &scores, &mut al);
            align_vectorized_batch_subs::<L>(&pairs, &scores, &mut fused);
            for lane in 0..L {
                let (al0, al1) = &al.al[lane];
                let want = crate::nwalign::al2subs(al0, al1);
                assert_eq!(
                    fused.subs[lane], want,
                    "case {case} lane {lane}: la={la} lb={lb} band={band} ends=({begin},{end})"
                );
                nsubs += want.nsubs();
                nn += al0
                    .iter()
                    .zip(al1.iter())
                    .filter(|&(&x, &y)| x == 5 || y == 5)
                    .count();
            }
        }
        // The cases must exercise substitutions and N columns, not only identity.
        assert!(
            nsubs > 10_000 && nn > 1_000,
            "nsubs {nsubs}, N columns {nn}"
        );
    }

    #[test]
    fn batch_matches_single_l8() {
        check_identity::<8>(3000, 0x9E3779B97F4A7C15);
    }

    #[test]
    fn batch_matches_single_l16() {
        check_identity::<16>(1500, 0xD1B54A32D192ED03);
    }

    /// Release builds no longer clear the pointer matrix (`reset_pointers`).
    /// Fill it with valid-looking pointers before every call: a traceback that
    /// read an unwritten cell would follow one and change the alignment. Only
    /// meaningful without debug assertions (`cargo test --release`), since
    /// debug builds still zero the matrix.
    #[test]
    fn stale_pointers_do_not_leak() {
        let mut st = 0xA076_1D64_78BD_642Fu64;
        let mut single = AlignBuffers::new();
        let mut batch = BatchBuffers::new();
        for case in 0..1500 {
            let la = 20 + (rng(&mut st) % 300) as usize;
            let lb = la.saturating_sub(12) + (rng(&mut st) % 25) as usize;
            let band = [-1, 1, 2, 3, 7, 8, 16, 17, 32][(rng(&mut st) % 9) as usize];
            let (begin, end) =
                [(0i16, 0i16), (-8, -8), (0, -8), (-8, 0)][(rng(&mut st) % 4) as usize];
            let scores = VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                begin_gap_p: begin,
                end_gap_p: end,
                band,
            };
            let center: Vec<u8> = (0..la).map(|_| 1 + (rng(&mut st) % 4) as u8).collect();
            let raws: Vec<Vec<u8>> = (0..16).map(|_| variant(&center, lb, &mut st)).collect();
            let pairs: [(&[u8], &[u8]); 16] = std::array::from_fn(|k| {
                if case % 2 == 0 {
                    (center.as_slice(), raws[k].as_slice())
                } else {
                    (raws[k].as_slice(), center.as_slice())
                }
            });
            let n = 64 * (la + lb + 2) * 16;
            batch.p.resize(n.max(batch.p.len()), 0);
            for x in batch.p.iter_mut() {
                *x = 1 + (rng(&mut st) % 3) as u8;
            }
            align_vectorized_batch::<16>(&pairs, &scores, &mut batch);
            for (lane, &(a, b)) in pairs.iter().enumerate() {
                single.poison_pointers(n / 16, &mut st);
                let got = align_single(a, b, &scores, &mut single);
                let want = align_single(a, b, &scores, &mut AlignBuffers::new());
                assert_eq!(got, want, "single, case {case} lane {lane}");
                assert_eq!(batch.al[lane], want, "batch, case {case} lane {lane}");
            }
        }
    }

    #[test]
    fn batch_matches_single_l1() {
        check_identity::<1>(2000, 0x2545F4914F6CDD1D);
    }
}
