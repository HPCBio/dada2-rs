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
use crate::minimizers;
// The experimental WFA backend lives in [`crate::wfa`]; the ends-free dispatch
// below routes to it when selected. `align_wfa_endsfree_with_buf` and
// `wfa_cost_cap` exist in every build (a stub / pure arithmetic respectively);
// the `DADA2RS_ALIGN_BACKEND=wfa` override switch is compiled only with the
// `wfa` feature.
#[cfg(feature = "wfa")]
use crate::wfa::USE_WFA_BACKEND;
use crate::wfa::{align_wfa_endsfree_with_buf, wfa_cost_cap};

/// Sentinel used in `Sub::map` to indicate that a reference position aligns
/// to a gap in the query.  Matches C++ `GAP_GLYPH = 9999`.
pub const GAP_GLYPH: u16 = 9999;

/// Score sentinel for out-of-band DP cells.
const BAND_SENTINEL: i32 = -9999;

// ---------------------------------------------------------------------------
// AlignParams
// ---------------------------------------------------------------------------

/// Pairwise-alignment backend (issue #49).
///
/// `Nw` is the default scalar/vectorized Needleman-Wunsch path. `Wfa2` routes
/// the ends-free path through the experimental WFA backend (wfa2lib-rs). WFA is
/// ASV-equivalent to NW on tested Illumina and PacBio HiFi data but its
/// alignments are not byte-identical (see `sweep_wfa_parity`).
#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    clap::ValueEnum,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum AlignBackend {
    /// Needleman-Wunsch (scalar/vectorized), the default.
    #[default]
    Nw,
    /// Experimental wavefront alignment via wfa2lib-rs.
    Wfa2,
}

/// Which pre-alignment screen decides whether a pair is worth aligning.
///
/// Both backends feed the same `kdist_cutoff` gate and neither defines the
/// clusters — the screen only avoids alignments (see
/// `docs/findings/kmer-size-screening.md`). They differ in *which* k-mers are
/// consulted: `Kmer` compares full `4^k` frequency vectors, `Minimizer` compares
/// a winnowed sketch. Experimental, opt-in, and default `Kmer`, on the same
/// footing as [`AlignBackend::Wfa2`]: this deviates from the R/C++ reference,
/// whose screen is the ESPRIT frequency vector.
#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    clap::ValueEnum,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum ScreenBackend {
    /// ESPRIT-style `4^k` k-mer frequency vectors, the default and the R
    /// reference's behaviour.
    #[default]
    Kmer,
    /// Experimental winnowed-minimizer sketch (see [`crate::minimizers`]).
    Minimizer,
}

/// Parameters controlling alignment method selection in `raw_align`.
#[derive(Clone, Copy)]
pub struct AlignParams {
    /// Which pairwise-alignment backend to use for the ends-free path.
    pub backend: AlignBackend,
    /// WFA edit-budget cap (issue #51), in edit operations. WFA aborts once the
    /// alignment needs more than this many edits and the pair falls back to the
    /// banded NW path (byte-identical to the NW backend for that pair). Bounds
    /// WFA's O(n·s) cost on divergent non-error-copy pairs that survive the
    /// k-mer screen; real error-copies sit far under any sane budget given the
    /// ~99.9% denoising-identity assumption. `0` = unbounded. Ignored for the NW
    /// backend. Converted to a WFA cost via [`wfa_cost_cap`].
    pub wfa_max_edits: i32,
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_p: i32,
    pub homo_gap_p: i32,
    pub use_kmers: bool,
    pub kdist_cutoff: f64,
    /// Which pre-alignment screen to use (experimental; default
    /// [`ScreenBackend::Kmer`]).
    pub screen_backend: ScreenBackend,
    /// K-mer size for the minimizer sketch. Ignored unless `screen_backend` is
    /// [`ScreenBackend::Minimizer`]. Independent of `kmer_size`: the sketch's
    /// memory does not scale with it, so it is chosen purely for discrimination.
    pub minimizer_k: usize,
    /// Winnowing window, in k-mers, for the minimizer sketch. Ignored unless
    /// `screen_backend` is [`ScreenBackend::Minimizer`].
    pub minimizer_w: usize,
    /// Diagnostic: evaluate BOTH screens on every comparison and align the
    /// union, recording the pair-level agreement in [`minimizers::audit`].
    /// Requires `screen_backend == Minimizer` and both screens built on each
    /// `Raw`. Changes the work done (and therefore the timings) but not the
    /// backend's own screen verdicts, so the ASVs are the minimizer backend's.
    pub screen_audit: bool,
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
    /// Screen-audit state: whether the last comparison passed each screen.
    /// Written by `raw_align_with_buf`, consumed by `sub_new_with_buf`, which is
    /// the first point at which the substitution count exists.
    pub audit_kmer_pass: bool,
    pub audit_mini_pass: bool,
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
    /// When true, [`raw_align_with_buf`] and [`sub_new_with_buf`] record the
    /// k-mer-screen and alignment times of their most recent call into
    /// `last_screen_nanos` / `last_align_nanos` (issue #127).
    ///
    /// Off in production runs: the timer pair costs ~40-50 ns against a k-mer
    /// screen that can itself be under a microsecond, so leaving it always-on
    /// would distort the very ratio it exists to measure. Enabled only under
    /// `--verbose`, where the reported numbers carry that overhead explicitly.
    pub measure: bool,
    /// Nanos in the k-mer screen of the last call. Paid on *every* comparison.
    pub last_screen_nanos: u64,
    /// Nanos in the DP alignment kernel proper (traceback included, `al2subs`
    /// excluded) of the last call.
    ///
    /// Split out from `last_align_nanos` because the first MiSeq measurement put
    /// `align+subs` at ~15.4 us per aligned pair — around 1.9 ns per DP cell for
    /// a ~250 bp band-16 alignment, far above what a vectorized i16 kernel
    /// should cost. Without this split the cause is ambiguous between a slow
    /// kernel and expensive post-processing, so the number cannot be acted on.
    pub last_dp_nanos: u64,
    /// Nanos in `al2subs` + quality mapping of the last call — the
    /// post-alignment work in [`sub_new_with_buf`].
    pub last_post_nanos: u64,
    /// Nanos in DP alignment + `al2subs` + quality mapping of the last call,
    /// i.e. `last_dp_nanos + last_post_nanos`.
    /// Zero when the screen shrouded the pair, which is the whole point of the
    /// split: the screen's cost is unavoidable, the aligner's is what it buys.
    pub last_align_nanos: u64,
}

impl AlignBuffers {
    pub fn new() -> Self {
        Self::default()
    }

    /// Enable per-call screen/align timing (see [`AlignBuffers::measure`]).
    pub fn measuring() -> Self {
        Self {
            measure: true,
            ..Self::default()
        }
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

    // The `d` score matrix is a **3-row rolling buffer**, not the full
    // `ncol × nrow` matrix (issue #127). The DP recurrence for anti-diagonal
    // `row` reads only rows `row-1` (left/up) and `row-2` (diag), and the
    // traceback below reads `p` alone and never touches `d` — so every `d` row
    // older than two is write-only. Keeping the full matrix cost 4 of the 6
    // bytes of memory traffic per DP cell (a zeroing pass plus the write, over
    // `i16` scores) for data that is never read back; on a 1.5 kbp HiFi pair at
    // band 32 that is a 199 KB per-alignment allocation retired down to 210 B.
    //
    // `p` stays full-size: the traceback genuinely needs every pointer.
    //
    // Dropping the `d` zeroing pass is safe because no cell is read before it
    // is written. Within a row, `dploop` writes exactly `[col_min, col_max]`;
    // the only cells outside that range which the next two rows read are
    // `col_min-1` and `col_max+1`, and each of those is *always* either an
    // edge-fill cell (while the ends-free edges are still active) or one of the
    // four sentinel columns `{0, 1, ncol-2, ncol-1}`. Both are rewritten for
    // every row below, which the rolling buffer requires in any case since the
    // slot carries stale values from three rows back.
    const DROWS: usize = 3;
    reset_buf(&mut buf.d16, ncol * DROWS, 0i16);
    reset_buf(&mut buf.p8, ncol * nrow, 0u8);
    reset_buf(&mut buf.diag_buf, ncol, 0i16);
    let p = &mut buf.p8[..ncol * nrow];
    let diag_buf = &mut buf.diag_buf[..ncol];

    // Sentinel value for the hard band-boundary columns. Chosen so that
    // fill_val + gap_p doesn't overflow.
    let min_score = mismatch.min(gap_p).min(match_score).min(0);
    let fill_val = i16::MIN.wrapping_sub(min_score);

    // Rolling slots, held in DP order: `[diag (row-2), prev (row-1), cur]`.
    // Rotating left by one after each anti-diagonal re-labels them for the next
    // row (old `prev` → `diag`, old `cur` → `prev`, old `diag` → the slot about
    // to be overwritten), which keeps the three borrows disjoint by position.
    let (d_a, d_rest) = buf.d16[..ncol * DROWS].split_at_mut(ncol);
    let (d_b, d_c) = d_rest.split_at_mut(ncol);
    let mut slots: [&mut [i16]; DROWS] = [d_a, d_b, d_c];

    // Per-row initialisation, applied to `cur` as each anti-diagonal is
    // reached rather than in an up-front pass over the whole matrix.
    // `row_free` is the accumulated ends-free edge score at this row; it is
    // carried incrementally so that saturation matches the original repeated
    // `saturating_add` exactly.
    let left_fill_limit = 1 + band_usize.min(len1);
    let top_fill_limit = 1 + (band_usize + len2 - len1).min(len2);
    // Column of the left/top edge cell at `row`, matching the original
    // decrement-after-even / increment-after-odd walks.
    let left_col = |row: usize| (start_col - 1).saturating_sub((row - 1) / 2);
    let top_col = |row: usize| start_col + row / 2;

    // Rows 0 and 1 are consumed by the first main-loop iteration (row 2), so
    // they are initialised here; rows ≥ 2 are initialised inside the loop.
    for (row, slot) in slots.iter_mut().enumerate().take(2) {
        slot[0] = fill_val;
        slot[1] = fill_val;
        slot[ncol - 2] = fill_val;
        slot[ncol - 1] = fill_val;
        if row == 0 {
            // Starting cell (0,0) in compressed coordinates.
            slot[start_col] = 0;
        } else {
            if row < left_fill_limit {
                slot[left_col(row)] = end_gap_p;
                p[row * ncol + left_col(row)] = 3;
            }
            if row < top_fill_limit {
                slot[top_col(row)] = end_gap_p;
                p[row * ncol + top_col(row)] = 2;
            }
        }
    }
    p[start_col] = 0;
    // Edge score for row 1; advanced once per row below.
    let mut row_free = end_gap_p;

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

        // --- Initialise this row's slot ---
        // The slot still holds row-3's values, so the sentinels and the
        // ends-free edge cells are rewritten here. This replaces the original
        // up-front passes over the full matrix; ordering within the row is
        // preserved (sentinels first, then edge fills, then the DP), which
        // matters where an edge cell lands on a sentinel column.
        row_free = row_free.saturating_add(end_gap_p);
        {
            let cur = &mut slots[2];
            cur[0] = fill_val;
            cur[1] = fill_val;
            cur[ncol - 2] = fill_val;
            cur[ncol - 1] = fill_val;
            if row < left_fill_limit {
                cur[left_col(row)] = row_free;
                p[row * ncol + left_col(row)] = 3;
            }
            if row < top_fill_limit {
                cur[top_col(row)] = row_free;
                p[row * ncol + top_col(row)] = 2;
            }
        }

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
            let d_prev2 = &slots[0][col_min..col_min + n];
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
        // Offsets are now *relative to the row*, since each rolling slot is
        // exactly one row wide.
        let even_off = if even { 1 } else { 0 }; // even=true → subtract 1
        let left_off = col_min - even_off;
        let up_off = col_min + 1 - even_off;

        let p_cur = &mut p[row * ncol..row * ncol + ncol];

        // `slots` is held as [diag(row-2), prev(row-1), cur]; destructuring
        // gives the disjoint borrows `dploop` needs.
        let [_, d_prev, d_cur] = &mut slots;
        dploop(
            d_cur, p_cur, d_prev, diag_buf, left_off, up_off, col_min, col_min, n, gap_p, swap,
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
        // `left_off`/`up_off` index the previous row (`slots[1]`); `col_min` and
        // `col_max` may already have been advanced by the band transition
        // above, matching the original absolute-index arithmetic.
        if end_gap_p > gap_p {
            // Left boundary: gap in s2 extending to end of s1
            if recalc_left {
                let d_free = slots[1][left_off].saturating_add(end_gap_p);
                let cur = slots[2][col_min];
                let pcur = p[row * ncol + col_min];
                if d_free > cur {
                    slots[2][col_min] = d_free;
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
                let d_free = slots[1][up_off + col_max - col_min].saturating_add(end_gap_p);
                let cur = slots[2][col_max];
                let pcur = p[row * ncol + col_max];
                if d_free > cur {
                    slots[2][col_max] = d_free;
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

        // Re-label the rolling slots for the next anti-diagonal: old `prev`
        // becomes `diag`, old `cur` becomes `prev`, and old `diag` becomes the
        // slot that row+1 overwrites.
        slots.rotate_left(1);

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

/// One-line, human-readable descriptor of the alignment backend in effect, for
/// `--verbose` logs so archived runs are self-labeling (issue #51). The WFA
/// edit-budget cap is shown only for the WFA backend (it is inert under NW).
pub fn backend_repr(p: &AlignParams) -> String {
    match p.backend {
        AlignBackend::Nw => "alignment backend: nw (Needleman-Wunsch)".to_string(),
        AlignBackend::Wfa2 => {
            let cap = if p.wfa_max_edits > 0 {
                format!("wfa-max-edits={}", p.wfa_max_edits)
            } else {
                "wfa-max-edits=0 (unbounded)".to_string()
            };
            format!("alignment backend: wfa2 (experimental WFA; {cap})")
        }
    }
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
    raw_align_with_screen(raw1, raw2, p, buf, None)
}

/// [`raw_align_with_buf`] with the screen distance optionally supplied by the
/// caller.
///
/// `screen` short-circuits the pairwise screen computation only — the value is
/// used exactly as the computed one would have been, for the cutoff gate and for
/// the DP's banding decisions alike. Its purpose is
/// [`crate::minimizers::MinimizerIndex`], which derives the same distance for
/// every raw against one cluster center in a single scatter pass, making the
/// per-pair merge-join redundant. Passing a value that differs from what the
/// pairwise path would compute changes the run's output.
pub fn raw_align_with_screen(
    raw1: &Raw,
    raw2: &Raw,
    p: &AlignParams,
    buf: &mut AlignBuffers,
    screen: Option<f64>,
) -> Option<()> {
    // --- K-mer screening ---
    // Timed under `buf.measure` (#127). The screen and the alignment are the
    // two halves of a comparison's cost, and only the second is avoidable —
    // the screen runs on all of them, the aligner only on those it passes.
    let t_screen = buf.measure.then(std::time::Instant::now);
    if buf.measure {
        buf.last_screen_nanos = 0;
        buf.last_dp_nanos = 0;
        buf.last_post_nanos = 0;
        buf.last_align_nanos = 0;
    }
    let mut kdist = 0.0f64;
    let mut kodist = -1.0f64; // sentinel: different from kdist when use_kmers=false

    if p.use_kmers {
        let k = p.kmer_size;
        kdist = match (screen, p.screen_backend) {
            (Some(d), _) => d,
            (None, backend) => match backend {
                ScreenBackend::Kmer => {
                    // Prefer 8-bit kmer distance; fall back to 16-bit on overflow.
                    match (&raw1.kmer8, &raw2.kmer8) {
                        (Some(k1), Some(k2)) => {
                            let d8 = k1.dist8(k2, raw1.len(), raw2.len(), k);
                            if d8 < 0.0 {
                                // Overflow (a k-mer occurs ≥255× in both seqs): fall
                                // back to the exact 16-bit distance. The u16 vectors
                                // are not kept resident (issue #32), so recompute
                                // them from sequence here — this path is essentially
                                // never hit for amplicon data.
                                let v1 = assign_kmer(&raw1.seq, k);
                                let v2 = assign_kmer(&raw2.seq, k);
                                kmer_dist(&v1, raw1.len(), &v2, raw2.len(), k)
                            } else {
                                d8
                            }
                        }
                        // No 8-bit vectors built (e.g. cluster-center raws): no
                        // k-mer screen, exactly as before — previously both kmer8
                        // and the u16 kmer were absent together, yielding
                        // kdist = 0.0.
                        _ => 0.0,
                    }
                }
                ScreenBackend::Minimizer => screen_dist_minimizer(raw1, raw2),
            },
        };

        if p.screen_audit {
            // Evaluate the *other* screen too and stash both verdicts. `kdist`
            // itself is left as the active backend's value so the alignment path
            // (banding, the gapless shortcut) is exactly what that backend would
            // have taken — the audit changes which pairs are aligned, never how.
            let other = match p.screen_backend {
                ScreenBackend::Minimizer => match (&raw1.kmer8, &raw2.kmer8) {
                    (Some(k1), Some(k2)) => {
                        let d8 = k1.dist8(k2, raw1.len(), raw2.len(), k);
                        if d8 < 0.0 { 0.0 } else { d8 }
                    }
                    _ => 0.0,
                },
                ScreenBackend::Kmer => screen_dist_minimizer(raw1, raw2),
            };
            let (kmer_d, mini_d) = match p.screen_backend {
                ScreenBackend::Minimizer => (other, kdist),
                ScreenBackend::Kmer => (kdist, other),
            };
            buf.audit_kmer_pass = kmer_d <= p.kdist_cutoff;
            buf.audit_mini_pass = mini_d <= p.kdist_cutoff;
        }

        if p.gapless
            && let (Some(o1), Some(o2)) = (&raw1.kord, &raw2.kord)
        {
            kodist = kord_dist(o1, raw1.len(), o2, raw2.len(), k);
        }
    }

    if let Some(t) = t_screen {
        buf.last_screen_nanos = t.elapsed().as_nanos() as u64;
    }

    if p.use_kmers && kdist > p.kdist_cutoff {
        if p.screen_audit && (buf.audit_kmer_pass || buf.audit_mini_pass) {
            // The active backend shrouded this pair but the other backend would
            // have aligned it. Under audit we align anyway, so the disagreement
            // can be bucketed by the substitution count the aligner finds —
            // otherwise a disagreement is just a count, with no way to tell a
            // lost error copy from a correctly-rejected stranger.
        } else {
            if p.screen_audit {
                minimizers::audit::record(buf.audit_kmer_pass, buf.audit_mini_pass, None);
            }
            return None; // Outside k-mer distance threshold → NULL alignment.
        }
    }

    let t_align = buf.measure.then(std::time::Instant::now);
    let r = raw_align_dp(raw1, raw2, p, buf, kdist, kodist);
    if let Some(t) = t_align {
        buf.last_dp_nanos = t.elapsed().as_nanos() as u64;
        buf.last_align_nanos = buf.last_dp_nanos;
    }
    r
}

/// Gapless-shortcut hit counters, active only under `--screen-audit`.
pub static GAPLESS_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static GAPLESS_TOTAL: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Minimizer-sketch screen distance for a pair (experimental,
/// [`ScreenBackend::Minimizer`]).
///
/// Thin wrapper over [`minimizers::screen_dist`], which owns the fail-open rule
/// so this path and the index path in `b_compare_parallel` cannot drift apart.
#[inline]
fn screen_dist_minimizer(raw1: &Raw, raw2: &Raw) -> f64 {
    let shared = match (&raw1.minimizers, &raw2.minimizers) {
        (Some(m1), Some(m2)) => minimizers::shared_count(m1, m2),
        _ => 0,
    };
    minimizers::screen_dist(shared, raw1.minimizers.as_ref(), raw2.minimizers.as_ref())
}

/// The alignment half of [`raw_align_with_buf`], split out so the k-mer screen
/// and the DP can be timed separately (#127). Behaviour is unchanged: this is
/// the body that followed the screen's early return.
fn raw_align_dp(
    raw1: &Raw,
    raw2: &Raw,
    p: &AlignParams,
    buf: &mut AlignBuffers,
    kdist: f64,
    kodist: f64,
) -> Option<()> {
    // --- Method selection ---
    let take_gapless = p.band == 0 || (p.gapless && (kodist - kdist).abs() < f64::EPSILON);
    if p.screen_audit {
        // The gapless shortcut fires on `kodist == kdist`: the positional and
        // compositional k-mer distances agreeing means no shifts, hence no
        // indels. That predicate is only meaningful when `kdist` comes from the
        // SAME k-mer space as `kodist` -- under the minimizer backend `kdist` is
        // a sketch distance, so the equality essentially never holds and this
        // path silently stops firing. Counted so that claim is measured, not
        // assumed.
        GAPLESS_TOTAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if take_gapless {
            GAPLESS_HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
    if take_gapless {
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
    // Experimental WFA backend (issue #49), selected via `p.backend` (CLI
    // `--align-backend wfa2`) or the undocumented `DADA2RS_ALIGN_BACKEND=wfa`
    // override (handy for sweeps). Replaces only the vectorized/ends-free path;
    // the gapless fast-path (above) and the homopolymer-aware path (below) are
    // left untouched. NOTE: WFA ends-free is not byte-identical to align_endsfree
    // (see sweep_wfa_parity / wfa_endsfree_known_divergence), though it is
    // ASV-equivalent on tested data. Root cause is upstream WFA2-lib #102
    // (suboptimal endsfree alignment when match score != 0): DADA2's match=+5
    // means a free end-gap changes the number of scored columns, which WFA's
    // cost-model pruning can't see, so it under-credits free end-gaps. The edit
    // cap does NOT fix this (these are low-edit pairs that finish under budget);
    // it only bounds the divergent-pair cost. Tracked for the fork's ends-free
    // handling — see project_wfa_dependency_wfa2lib_rs.
    #[cfg(feature = "wfa")]
    if (p.backend == AlignBackend::Wfa2 || *USE_WFA_BACKEND) && !use_homo {
        align_wfa_endsfree_with_buf(
            &raw1.seq,
            &raw2.seq,
            p.match_score,
            p.mismatch,
            p.gap_p,
            p.band,
            wfa_cost_cap(p.wfa_max_edits, p.gap_p),
            buf,
        );
        return Some(());
    }
    // NW-only build: WFA code is not compiled (issue #63). Selecting the WFA
    // backend aborts via the stub rather than silently aligning with NW.
    #[cfg(not(feature = "wfa"))]
    if p.backend == AlignBackend::Wfa2 && !use_homo {
        align_wfa_endsfree_with_buf(
            &raw1.seq,
            &raw2.seq,
            p.match_score,
            p.mismatch,
            p.gap_p,
            p.band,
            wfa_cost_cap(p.wfa_max_edits, p.gap_p),
            buf,
        );
        return Some(());
    }
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
    sub_new_with_screen(raw0, raw1, params, buf, None)
}

/// [`sub_new_with_buf`] with the screen distance optionally supplied by the
/// caller; see [`raw_align_with_screen`].
pub fn sub_new_with_screen(
    raw0: &Raw,
    raw1: &Raw,
    params: &AlignParams,
    buf: &mut AlignBuffers,
    screen: Option<f64>,
) -> Option<Sub> {
    raw_align_with_screen(raw0, raw1, params, buf, screen)?;
    let audit = params.screen_audit && params.use_kmers;
    // Post-alignment work is charged to the alignment half (#127): it is paid
    // only by pairs the screen let through, so it belongs to what the aligner
    // costs, not to the screen's unavoidable baseline.
    let t_post = buf.measure.then(std::time::Instant::now);
    let mut sub = al2subs(&buf.al0, &buf.al1);

    if let (Some(q0), Some(q1)) = (&raw0.qual, &raw1.qual) {
        sub.q0 = sub.pos.iter().map(|&pos| q0[pos as usize]).collect();
        sub.q1 = sub
            .pos
            .iter()
            .map(|&pos| q1[sub.map[pos as usize] as usize])
            .collect();
    }
    if let Some(t) = t_post {
        buf.last_post_nanos = t.elapsed().as_nanos() as u64;
        buf.last_align_nanos += buf.last_post_nanos;
    }
    if audit {
        minimizers::audit::record(buf.audit_kmer_pass, buf.audit_mini_pass, Some(sub.nsubs()));
    }
    // Under audit the gate is the union of both screens, so a pair the ACTIVE
    // backend shrouded may have been aligned purely to classify it. Returning
    // its Sub would let the audit change the run's output, which would defeat
    // the point of auditing that backend.
    if audit && !screen_pass_for_backend(params, buf) {
        return None;
    }
    Some(sub)
}

/// Whether the active backend's own screen passed the last comparison, for
/// discarding audit-only alignments. Reads the verdicts stashed by
/// `raw_align_with_buf`.
#[inline]
fn screen_pass_for_backend(p: &AlignParams, buf: &AlignBuffers) -> bool {
    match p.screen_backend {
        ScreenBackend::Kmer => buf.audit_kmer_pass,
        ScreenBackend::Minimizer => buf.audit_mini_pass,
    }
}

#[cfg(test)]
mod bench_align {
    use super::*;

    /// Long-read WFA-vs-NW kernel benchmark (issue #51). Times a single fixed
    /// near-identical ~1.5 kb pair (representative of k-mer-screened
    /// learn-errors pairs) across alignment configs to localize the PacBio
    /// slowdown. Run:
    ///   cargo test --release -- --ignored bench_wfa_long --nocapture
    #[test]
    #[ignore]
    #[cfg(feature = "wfa")]
    fn bench_wfa_long() {
        use std::time::Instant;
        use wfa2lib_rs::aligner::{AffineAligner, AlignmentScope};
        use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

        // Deterministic ~1450 bp pair, near-identical: a handful of subs + indels,
        // like the similar sequences the aligner actually sees after k-mer screen.
        let len: usize = 1450;
        let nts = [1u8, 2, 3, 4];
        let mut st: u64 = 0x00C0_FFEE_1234_5678;
        let mut rng = |st: &mut u64, m: usize| {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*st >> 33) as usize) % m
        };
        let s1: Vec<u8> = (0..len).map(|_| nts[rng(&mut st, 4)]).collect();
        let mut s2 = s1.clone();
        for _ in 0..6 {
            let p = rng(&mut st, s2.len());
            s2[p] = nts[rng(&mut st, 4)];
        }
        for _ in 0..2 {
            let p = rng(&mut st, s2.len());
            s2.remove(p);
        }
        let band = 32i32;
        let iters = 3000usize;

        let bench = |label: &str, mut f: Box<dyn FnMut()>| {
            for _ in 0..100 {
                f();
            } // warmup
            let t0 = Instant::now();
            for _ in 0..iters {
                f();
            }
            let us = t0.elapsed().as_secs_f64() * 1e6 / iters as f64;
            println!("  {label:<34} {us:8.2} us/align");
            us
        };

        // Affine penalties matching the NW default scheme (Eizenga-adjusted).
        let penalties = || {
            WavefrontPenalties::new_affine(AffinePenalties {
                match_: -5,
                mismatch: 4,
                gap_opening: 0,
                gap_extension: 8,
            })
        };

        println!(
            "\nlen1={} len2={} band={band} iters={iters}",
            s1.len(),
            s2.len()
        );

        // --- NW baseline (vectorized, ends-free) ---
        {
            let scores = VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                end_gap_p: 0,
                band,
            };
            let mut buf = AlignBuffers::new();
            let (s1, s2) = (s1.clone(), s2.clone());
            bench(
                "NW vectorized (baseline)",
                Box::new(move || align_vectorized_with_buf(&s1, &s2, &scores, &mut buf)),
            );
        }

        // --- WFA, full ends-free, ComputeAlignment (current production path) ---
        {
            let mut a = AffineAligner::new(penalties());
            a.alignment_scope = AlignmentScope::ComputeAlignment;
            let (s1, s2) = (s1.clone(), s2.clone());
            bench(
                "WFA full ends-free + CIGAR",
                Box::new(move || {
                    a.set_alignment_free_ends(
                        s1.len() as i32,
                        s1.len() as i32,
                        s2.len() as i32,
                        s2.len() as i32,
                    );
                    a.align_endsfree(&s1, &s2);
                    std::hint::black_box(a.cigar().operations_slice());
                }),
            );
        }

        // --- WFA, bounded ends-free (~band), ComputeAlignment (hypothesis #1) ---
        {
            let mut a = AffineAligner::new(penalties());
            a.alignment_scope = AlignmentScope::ComputeAlignment;
            let (s1, s2) = (s1.clone(), s2.clone());
            bench(
                "WFA bounded ends-free + CIGAR",
                Box::new(move || {
                    a.set_alignment_free_ends(band, band, band, band);
                    a.align_endsfree(&s1, &s2);
                    std::hint::black_box(a.cigar().operations_slice());
                }),
            );
        }

        // --- WFA, full ends-free, score-only (isolates traceback cost, #2) ---
        {
            let mut a = AffineAligner::new(penalties());
            a.alignment_scope = AlignmentScope::ComputeScore;
            let (s1, s2) = (s1.clone(), s2.clone());
            bench(
                "WFA full ends-free, score-only",
                Box::new(move || {
                    a.set_alignment_free_ends(
                        s1.len() as i32,
                        s1.len() as i32,
                        s2.len() as i32,
                        s2.len() as i32,
                    );
                    std::hint::black_box(a.align_endsfree(&s1, &s2));
                }),
            );
        }

        // --- WFA end2end + CIGAR (no free-ends at all; lower bound on WFA cost) ---
        {
            let mut a = AffineAligner::new(penalties());
            a.alignment_scope = AlignmentScope::ComputeAlignment;
            let (s1, s2) = (s1.clone(), s2.clone());
            bench(
                "WFA end2end + CIGAR",
                Box::new(move || {
                    a.align_end2end(&s1, &s2);
                    std::hint::black_box(a.cigar().operations_slice());
                }),
            );
        }

        // --- BiWFA + CIGAR (O(s) memory, long-read oriented, #3) ---
        {
            let mut a = AffineAligner::new(penalties());
            a.alignment_scope = AlignmentScope::ComputeAlignment;
            let (s1, s2) = (s1.clone(), s2.clone());
            bench(
                "WFA BiWFA + CIGAR",
                Box::new(move || {
                    a.align_biwfa(&s1, &s2);
                    std::hint::black_box(a.cigar().operations_slice());
                }),
            );
        }
        println!();
    }

    /// WFA-vs-NW across divergence levels (issue #51). WFA is O(n·s); NW is
    /// banded O(n·band). This sweeps edit distance on a ~1.5 kb pair to find
    /// where WFA's cost crosses NW's — the suspected source of the pipeline
    /// slowdown (divergent-but-k-mer-screened pairs).
    ///   cargo test --release -- --ignored bench_wfa_divergence --nocapture
    #[test]
    #[ignore]
    #[cfg(feature = "wfa")]
    fn bench_wfa_divergence() {
        use std::time::Instant;
        use wfa2lib_rs::aligner::{AffineAligner, AlignmentScope};
        use wfa2lib_rs::heuristic::HeuristicStrategy;
        use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

        let len: usize = 1450;
        let nts = [1u8, 2, 3, 4];
        let mut st: u64 = 0xABCD_1234;
        let mut rng = |st: &mut u64, m: usize| {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*st >> 33) as usize) % m
        };
        let s1: Vec<u8> = (0..len).map(|_| nts[rng(&mut st, 4)]).collect();
        let band = 32i32;
        let iters = 1000usize;
        let penalties = || {
            WavefrontPenalties::new_affine(AffinePenalties {
                match_: -5,
                mismatch: 4,
                gap_opening: 0,
                gap_extension: 8,
            })
        };

        println!("\nlen={len} band={band} iters={iters}");
        println!(
            "  {:>6}  {:>10}  {:>12}  {:>8}  {:>12}  {:>10}",
            "edits", "NW us", "WFA us", "WFA/NW", "WFAband us", "band/NW"
        );
        for &nedits in &[0usize, 5, 10, 25, 50, 100, 200] {
            let mut s2 = s1.clone();
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
            // NW vectorized
            let scores = VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                end_gap_p: 0,
                band,
            };
            let mut buf = AlignBuffers::new();
            for _ in 0..50 {
                align_vectorized_with_buf(&s1, &s2, &scores, &mut buf);
            }
            let t0 = Instant::now();
            for _ in 0..iters {
                align_vectorized_with_buf(&s1, &s2, &scores, &mut buf);
            }
            let nw_us = t0.elapsed().as_secs_f64() * 1e6 / iters as f64;

            // WFA full ends-free + CIGAR (current production config)
            let mut a = AffineAligner::new(penalties());
            a.alignment_scope = AlignmentScope::ComputeAlignment;
            let setf = |a: &mut AffineAligner, s1: &[u8], s2: &[u8]| {
                a.set_alignment_free_ends(
                    s1.len() as i32,
                    s1.len() as i32,
                    s2.len() as i32,
                    s2.len() as i32,
                );
                a.align_endsfree(s1, s2);
            };
            for _ in 0..50 {
                setf(&mut a, &s1, &s2);
            }
            let t0 = Instant::now();
            for _ in 0..iters {
                setf(&mut a, &s1, &s2);
                std::hint::black_box(a.cigar().operations_slice());
            }
            let wfa_us = t0.elapsed().as_secs_f64() * 1e6 / iters as f64;

            // WFA banded (BandedStatic ±band) + CIGAR — NW-band equivalent (#51).
            let mut ab = AffineAligner::new(penalties());
            ab.alignment_scope = AlignmentScope::ComputeAlignment;
            ab.set_heuristic(HeuristicStrategy::BandedStatic {
                min_k: -band,
                max_k: band,
            });
            for _ in 0..50 {
                setf(&mut ab, &s1, &s2);
            }
            let t0 = Instant::now();
            for _ in 0..iters {
                setf(&mut ab, &s1, &s2);
                std::hint::black_box(ab.cigar().operations_slice());
            }
            let wfab_us = t0.elapsed().as_secs_f64() * 1e6 / iters as f64;

            println!(
                "  {nedits:>6}  {nw_us:>10.2}  {wfa_us:>12.2}  {:>7.2}x  {wfab_us:>12.2}  {:>9.2}x",
                wfa_us / nw_us,
                wfab_us / nw_us
            );
        }
        println!();
    }

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
    #[cfg(feature = "wfa")]
    fn compare_wfa(s1: &[u8], s2: &[u8], label: &str) {
        let ef = align_endsfree(s1, s2, MATCH, MM, GAP, BAND);
        let mut buf = AlignBuffers::new();
        align_wfa_endsfree_with_buf(s1, s2, MATCH, MM, GAP, BAND, i32::MAX, &mut buf);
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
    #[cfg(feature = "wfa")]
    fn test_wfa_vs_endsfree_identical() {
        let s = encode("ACGTACGTACGT");
        compare_wfa(&s, &s, "identical-short");
    }

    #[test]
    #[cfg(feature = "wfa")]
    fn test_wfa_vs_endsfree_one_sub() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTTCGTACGT");
        compare_wfa(&s1, &s2, "one-sub");
    }

    #[test]
    #[cfg(feature = "wfa")]
    fn test_wfa_vs_endsfree_one_gap() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGACGTACGT");
        compare_wfa(&s1, &s2, "one-gap");
    }

    #[test]
    #[cfg(feature = "wfa")]
    fn test_wfa_vs_endsfree_different_lengths() {
        let s1 = encode("ACGTACGTACGT");
        let s2 = encode("ACGTACGTACGTAC");
        compare_wfa(&s1, &s2, "diff-len-short");
    }

    // ---- WFA edit-budget cap (issue #51) -------------------------------------

    /// WFA *cost* (Eizenga convention: matches free, interior mismatches and
    /// gaps penalised) of an ends-free alignment — the unit `--wfa-max-edits` is
    /// converted into via [`wfa_cost_cap`]. Lets a test assert a pair genuinely
    /// exceeds (or sits under) a cap, so the NW-fallback checks aren't vacuous.
    #[cfg(feature = "wfa")]
    fn wfa_edit_cost(al: &[Vec<u8>; 2]) -> i32 {
        let (al0, al1) = (&al[0], &al[1]);
        let mut cost = 0;
        for k in 0..al0.len() {
            let g0 = al0[k] == b'-';
            let g1 = al1[k] == b'-';
            if !g0 && !g1 {
                if al0[k] != al1[k] {
                    cost += MM.abs();
                }
            } else {
                let row = if g0 { al0 } else { al1 };
                let interior = !row[..k].iter().all(|&b| b == b'-')
                    && !row[k + 1..].iter().all(|&b| b == b'-');
                if interior {
                    cost += GAP.abs();
                }
            }
        }
        cost
    }

    #[cfg(feature = "wfa")]
    fn nw_vectorized(s1: &[u8], s2: &[u8]) -> [Vec<u8>; 2] {
        align_vectorized(
            s1,
            s2,
            &VectorizedAlignScores {
                match_score: MATCH as i16,
                mismatch: MM as i16,
                gap_p: GAP as i16,
                end_gap_p: 0,
                band: BAND,
            },
        )
    }

    /// When a pair's alignment exceeds the edit budget, WFA aborts and the pair
    /// must fall back to the banded NW path — byte-identical to the NW backend.
    #[test]
    #[cfg(feature = "wfa")]
    fn wfa_cap_triggers_nw_fallback_byte_identical() {
        // ~6 substitutions over 20 bp: cost 24 > the 2-edit (cost-16) budget.
        let s1 = encode("ACGTACGTACGTACGTACGT");
        let s2 = encode("AGGTAAGTACCTACCTAGGT");
        let cap = wfa_cost_cap(2, GAP); // 2 edits -> cost 16
        let mut buf = AlignBuffers::new();
        align_wfa_endsfree_with_buf(&s1, &s2, MATCH, MM, GAP, BAND, cap, &mut buf);

        // Non-vacuity: the (uncapped) optimal really does exceed the budget, so
        // the fallback path is the one under test.
        let nw = nw_vectorized(&s1, &s2);
        assert!(
            wfa_edit_cost(&nw) > cap,
            "test pair must exceed the cap (cost {} <= cap {cap})",
            wfa_edit_cost(&nw),
        );

        assert_eq!(buf.al0, nw[0], "capped WFA al0 must equal NW fallback");
        assert_eq!(buf.al1, nw[1], "capped WFA al1 must equal NW fallback");
    }

    /// A pair within budget must be untouched by the cap: a generous cap yields
    /// exactly the uncapped WFA alignment (i.e. the fallback did NOT fire).
    #[test]
    #[cfg(feature = "wfa")]
    fn wfa_cap_within_budget_matches_uncapped() {
        let s1 = encode("ACGTACGTACGTACGT");
        let s2 = encode("ACGTACGTTCGTACGT"); // 1 substitution
        let cap = wfa_cost_cap(50, GAP); // 50 edits -> cost 400, far above

        let mut capped = AlignBuffers::new();
        let mut uncapped = AlignBuffers::new();
        align_wfa_endsfree_with_buf(&s1, &s2, MATCH, MM, GAP, BAND, cap, &mut capped);
        align_wfa_endsfree_with_buf(&s1, &s2, MATCH, MM, GAP, BAND, i32::MAX, &mut uncapped);

        assert_eq!(
            capped.al0, uncapped.al0,
            "in-budget al0 must match uncapped"
        );
        assert_eq!(
            capped.al1, uncapped.al1,
            "in-budget al1 must match uncapped"
        );

        // And it really is under budget, so the match is the WFA path, not NW.
        let cost = wfa_edit_cost(&[uncapped.al0.clone(), uncapped.al1.clone()]);
        assert!(
            cost <= cap,
            "pair should be within budget (cost {cost} > cap {cap})"
        );
    }

    /// `wfa_cost_cap` converts an edit budget into a WFA cost (`edits·|gap_p|`),
    /// with `<= 0` meaning unbounded and no overflow on saturation.
    #[test]
    fn wfa_cost_cap_converts_edits_to_cost() {
        assert_eq!(wfa_cost_cap(50, -8), 400);
        assert_eq!(wfa_cost_cap(40, -4), 160);
        assert_eq!(wfa_cost_cap(1, -8), 8);
        assert_eq!(wfa_cost_cap(3, 8), 24); // positive gap_p: |gap_p| is used
        assert_eq!(wfa_cost_cap(0, -8), i32::MAX); // 0 = unbounded
        assert_eq!(wfa_cost_cap(-5, -8), i32::MAX); // negative = unbounded
        assert_eq!(wfa_cost_cap(i32::MAX, -8), i32::MAX); // saturating, no panic
    }

    /// Isolation experiment: compare WFA *global* (`align_end2end`) against the
    /// scalar *global* `align_standard` (both linear gap, penalized ends). If
    /// these agree where ends-free diverges, the divergence is purely in
    /// ends-free free-end-gap crediting, not the interior gap model — meaning
    /// affine gap scoring would not address it.
    ///   cargo test --release --bins wfa_global_isolation -- --ignored --nocapture
    #[test]
    #[ignore]
    #[cfg(feature = "wfa")]
    fn wfa_global_isolation() {
        use wfa2lib_rs::aligner::{AffineAligner, AlignmentScope};
        use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

        fn wfa_global(s1: &[u8], s2: &[u8]) -> [Vec<u8>; 2] {
            let penalties = WavefrontPenalties::new_affine(AffinePenalties {
                match_: -5,
                mismatch: 4,
                gap_opening: 0,
                gap_extension: 8,
            });
            let mut a = AffineAligner::new(penalties);
            a.alignment_scope = AlignmentScope::ComputeAlignment;
            a.align_end2end(s1, s2);
            let mut al0 = Vec::new();
            let mut al1 = Vec::new();
            crate::wfa::cigar_to_alignment_into(
                a.cigar().operations_slice(),
                s1,
                s2,
                &mut al0,
                &mut al1,
            );
            [al0, al1]
        }

        let nts = [1u8, 2, 3, 4];
        let mut st: u64 = 0xDEAD_BEEF_0BAD_F00D;
        let mut rng = |st: &mut u64, m: usize| {
            *st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*st >> 33) as usize) % m
        };
        let mut global_fails = 0u32;
        let mut endsfree_fails = 0u32;
        let n = 10_000;
        for _ in 0..n {
            let len = 40 + rng(&mut st, 60);
            let s1: Vec<u8> = (0..len).map(|_| nts[rng(&mut st, 4)]).collect();
            let mut s2 = s1.clone();
            for _ in 0..rng(&mut st, 6) {
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
            // Global: WFA end2end vs scalar align_standard (band=-1, unbanded).
            let g_scalar = align_standard(&s1, &s2, 5, -4, -8, -1);
            let g_wfa = wfa_global(&s1, &s2);
            if score_alignment(&g_scalar, 5, -4, -8) != score_alignment(&g_wfa, 5, -4, -8) {
                global_fails += 1;
            }
            // Ends-free: WFA endsfree vs scalar align_endsfree, for comparison.
            let e_scalar = align_endsfree(&s1, &s2, 5, -4, -8, -1);
            let mut buf = AlignBuffers::new();
            align_wfa_endsfree_with_buf(&s1, &s2, 5, -4, -8, -1, i32::MAX, &mut buf);
            let e_wfa = [buf.al0.clone(), buf.al1.clone()];
            if score_alignment(&e_scalar, 5, -4, -8) != score_alignment(&e_wfa, 5, -4, -8) {
                endsfree_fails += 1;
            }
        }
        println!("  WFA global   vs align_standard : {global_fails}/{n} disagree");
        println!("  WFA endsfree vs align_endsfree : {endsfree_fails}/{n} disagree");
    }

    /// Find and print the first few low-edit pairs where WFA disagrees with
    /// `align_endsfree`, with both alignments, to diagnose the cause.
    ///   cargo test --release --bins wfa_diagnose -- --ignored --nocapture
    #[test]
    #[ignore]
    #[cfg(feature = "wfa")]
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
            align_wfa_endsfree_with_buf(&s1, &s2, 5, -4, -8, -1, i32::MAX, &mut buf);
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

    /// Deterministic characterization of the known WFA ends-free divergence
    /// (upstream WFA2-lib #102: suboptimal endsfree alignment when match != 0).
    ///
    /// `s2` is exactly `s1` with the leading base dropped, so the optimal
    /// ends-free alignment credits a FREE leading end-gap and scores every one
    /// of the 51 shared columns as a match (51 * 5 = 255). WFA cannot see that
    /// the free end-gap changes the scored-column count (its cost model treats
    /// match as 0), so it instead places the gap one column in — a penalized
    /// internal gap — and scores 255 - 8 = 247. This is NOT a tie-break or a
    /// logic bug on our side: it is the upstream paradigm mismatch, and the edit
    /// cap does NOT mask it (a 1-edit pair finishes far under any budget, so
    /// there is no NW fallback). It is the mechanism behind WFA's low-abundance
    /// ASV over-calls (jaccard ~0.999, not 1.000) seen at scale (issue #51).
    ///
    /// Asserted as a REGRESSION GUARD on the documented behavior: if a future
    /// WFA/fork change makes WFA optimal here, this test will fail loudly and
    /// should be updated to assert equality (i.e. the divergence is fixed).
    #[test]
    #[cfg(feature = "wfa")]
    fn wfa_endsfree_known_divergence() {
        // 52 nt; s2 = s1[1..] (leading base removed) — a single deletion.
        let s1 = encode("AACAGCGCAAACCAACTCGCTAGCTAGCAAAATCTTGTGTTTCTGCCTAGCG");
        let s2 = encode("ACAGCGCAAACCAACTCGCTAGCTAGCAAAATCTTGTGTTTCTGCCTAGCG");
        assert_eq!(s2.len(), s1.len() - 1);

        // align_endsfree finds the true optimum: free leading end-gap, every
        // shared column a match.
        let ef = align_endsfree(&s1, &s2, 5, -4, -8, -1);
        let sef = score_alignment(&ef, 5, -4, -8);
        assert_eq!(
            sef, 255,
            "align_endsfree should credit the free leading gap"
        );

        // WFA (uncapped) is strictly suboptimal here — the #102 symptom.
        let mut buf = AlignBuffers::new();
        align_wfa_endsfree_with_buf(&s1, &s2, 5, -4, -8, -1, i32::MAX, &mut buf);
        let wfa = [buf.al0.clone(), buf.al1.clone()];
        let swfa = score_alignment(&wfa, 5, -4, -8);
        assert!(
            swfa < sef,
            "WFA ends-free is expected to be suboptimal here (#102); \
             got wfa={swfa} vs endsfree={sef}. If WFA now matches, the upstream \
             divergence may be fixed — update this guard to assert equality."
        );
        assert_eq!(
            swfa, 247,
            "documented WFA score: one penalized internal gap"
        );
    }

    /// Randomized parity stress test: WFA must match `align_endsfree`'s optimal
    /// score across many random pairs (varied lengths, indels). Run explicitly:
    ///   cargo test --release -- --ignored sweep_wfa_parity --nocapture
    #[test]
    #[ignore]
    #[cfg(feature = "wfa")]
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
            align_wfa_endsfree_with_buf(&s1, &s2, 5, -4, -8, -1, i32::MAX, &mut buf);
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

    // -----------------------------------------------------------------------
    // Rolling-buffer differential test (issue #127)
    // -----------------------------------------------------------------------
    // Verbatim copy of the pre-#127 full-matrix `align_vectorized_with_buf`,
    // kept as a byte-identity oracle for the 3-row rolling `d16` rewrite. The
    // rolling version must produce the *same alignment strings*, not merely the
    // same optimal score (which `sweep_vectorized_parity` already checks against
    // `align_endsfree`) — a score-equal but differently-placed gap would change
    // `al2subs` output and therefore the error model.
    #[allow(clippy::needless_range_loop)]
    fn align_vectorized_reference(
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

    /// A/B the rolling-buffer DP against the retained full-matrix reference on
    /// synthetic pairs at both platforms' geometries.
    ///
    ///     cargo test --release -- --ignored bench_rolling_d16 --nocapture
    ///
    /// NOTE: this machine is not the benchmark node. The change is a
    /// memory-traffic reduction, so a laptop with large caches and high
    /// per-core bandwidth will *understate* the win — especially for the
    /// long-read geometry, whose 199 KB full matrix is the case that does not
    /// fit L2 on the target. Treat the ratio as a lower bound; see
    /// docs/findings/compare-screen-vs-align.md.
    #[test]
    #[ignore]
    fn bench_rolling_d16() {
        let mut state = 0x243F6A88_85A308D3u64;
        for &(len, band, label) in &[
            (240usize, 16i32, "MiSeq R1  (240bp, band 16)"),
            (160, 16, "MiSeq R2  (160bp, band 16)"),
            (1455, 32, "PacBio    (1455bp, band 32)"),
        ] {
            let s1: Vec<u8> = (0..len)
                .map(|_| 1 + (rng_next(&mut state) % 4) as u8)
                .collect();
            let mut s2 = s1.clone();
            for k in (0..len).step_by(50) {
                s2[k] = 1 + (s2[k] % 4);
            }
            let scores = VectorizedAlignScores {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                end_gap_p: 0,
                band,
            };
            let iters = if len > 1000 { 2_000 } else { 20_000 };
            let cells = (2 * len + 1) * (band as usize + 1);

            let mut buf = AlignBuffers::new();
            let mut time = |f: &mut dyn FnMut(&mut AlignBuffers)| {
                for _ in 0..iters / 10 {
                    f(&mut buf);
                }
                let t = std::time::Instant::now();
                for _ in 0..iters {
                    f(&mut buf);
                }
                t.elapsed().as_secs_f64() / iters as f64 * 1e9
            };
            let new_ns =
                time(&mut |b: &mut AlignBuffers| align_vectorized_with_buf(&s1, &s2, &scores, b));
            let old_ns =
                time(&mut |b: &mut AlignBuffers| align_vectorized_reference(&s1, &s2, &scores, b));
            println!(
                "  {label}: full-matrix {old_ns:>8.0} ns ({:.2} ns/cell)  ->  rolling {new_ns:>8.0} ns ({:.2} ns/cell)   {:.2}x",
                old_ns / cells as f64,
                new_ns / cells as f64,
                old_ns / new_ns
            );
        }
    }

    /// The same A/B run on N threads at once. The rolling buffer is a
    /// *memory-traffic* change, so its effect should appear only when enough
    /// cores contend for bandwidth — the condition on the 24-thread benchmark
    /// node, which a single-threaded bench cannot reproduce at all.
    ///
    ///     cargo test --release -- --ignored bench_rolling_d16_threaded --nocapture
    #[test]
    #[ignore]
    fn bench_rolling_d16_threaded() {
        for &(len, band, label) in &[
            (240usize, 16i32, "MiSeq R1  (240bp, band 16)"),
            (1455, 32, "PacBio    (1455bp, band 32)"),
        ] {
            println!("  {label}");
            // Ladder defaults to 1..=8 but should be taken up to the node's
            // full core count — bandwidth contention is the whole hypothesis,
            // and it does not appear until the cores are saturated:
            //   DADA2RS_BENCH_THREADS=1,2,4,8,12,24 cargo test --release ...
            let ladder: Vec<usize> = match std::env::var("DADA2RS_BENCH_THREADS") {
                Ok(v) => v
                    .split(',')
                    .filter_map(|x| x.trim().parse().ok())
                    .filter(|&n: &usize| n > 0)
                    .collect(),
                Err(_) => vec![1, 2, 4, 8],
            };
            for &nthreads in &ladder {
                let mut out = Vec::new();
                for reference in [true, false] {
                    let elapsed: f64 = std::thread::scope(|sc| {
                        let hs: Vec<_> = (0..nthreads)
                            .map(|t| {
                                sc.spawn(move || {
                                    let mut st = 0x243F6A88_85A308D3u64 ^ (t as u64);
                                    let s1: Vec<u8> = (0..len)
                                        .map(|_| 1 + (rng_next(&mut st) % 4) as u8)
                                        .collect();
                                    let mut s2 = s1.clone();
                                    for k in (0..len).step_by(50) {
                                        s2[k] = 1 + (s2[k] % 4);
                                    }
                                    let scores = VectorizedAlignScores {
                                        match_score: 5,
                                        mismatch: -4,
                                        gap_p: -8,
                                        end_gap_p: 0,
                                        band,
                                    };
                                    let iters = if len > 1000 { 1_000 } else { 10_000 };
                                    let mut buf = AlignBuffers::new();
                                    for _ in 0..iters / 10 {
                                        align_vectorized_with_buf(&s1, &s2, &scores, &mut buf);
                                    }
                                    let t0 = std::time::Instant::now();
                                    for _ in 0..iters {
                                        if reference {
                                            align_vectorized_reference(&s1, &s2, &scores, &mut buf);
                                        } else {
                                            align_vectorized_with_buf(&s1, &s2, &scores, &mut buf);
                                        }
                                    }
                                    t0.elapsed().as_secs_f64() / iters as f64 * 1e9
                                })
                            })
                            .collect();
                        hs.into_iter().map(|h| h.join().unwrap()).sum::<f64>() / nthreads as f64
                    });
                    out.push(elapsed);
                }
                println!(
                    "    {nthreads:>2} thread(s): full-matrix {:>8.0} ns  ->  rolling {:>8.0} ns   {:.2}x",
                    out[0],
                    out[1],
                    out[0] / out[1]
                );
            }
        }
    }

    /// Deterministic xorshift so failures are reproducible without a dep.
    fn rng_next(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    }

    /// Fuzz the rolling-buffer DP against the full-matrix reference and require
    /// byte-identical `al0`/`al1`.
    #[test]
    fn rolling_d16_matches_full_matrix() {
        let mut state = 0x2545F491_4F6CDD1Du64;
        let mut checked = 0usize;
        // Length pairs cover equal, shorter-first, longer-first and the
        // band-degenerate cases; bands cover both parities and unbanded (-1).
        for &(l1, l2) in &[
            (1usize, 1usize),
            (1, 5),
            (5, 1),
            (2, 3),
            (10, 10),
            (10, 13),
            (13, 10),
            (17, 16),
            (16, 17),
            (33, 33),
            (40, 40),
            (40, 55),
            (55, 40),
            (64, 64),
            (100, 97),
            (97, 100),
            (128, 128),
            (240, 240),
            (160, 160),
            (240, 160),
            (160, 240),
            (300, 299),
        ] {
            for &band in &[-1i32, 1, 2, 3, 4, 8, 15, 16, 17, 32, 64] {
                for trial in 0..12 {
                    // Trial 0 makes the sequences identical (the common case in
                    // production); later trials diverge progressively.
                    let s1: Vec<u8> = (0..l1)
                        .map(|_| 1 + (rng_next(&mut state) % 4) as u8)
                        .collect();
                    let s2: Vec<u8> = if trial == 0 && l1 == l2 {
                        s1.clone()
                    } else {
                        let mut v = if l2 <= l1 {
                            s1[..l2.min(l1)].to_vec()
                        } else {
                            s1.clone()
                        };
                        v.resize(l2, 1);
                        // Mutate an increasing fraction of positions.
                        for slot in v.iter_mut() {
                            if rng_next(&mut state) % 16 < trial as u64 {
                                *slot = 1 + (rng_next(&mut state) % 4) as u8;
                            }
                        }
                        v
                    };
                    for &end_gap_p in &[0i16, -8] {
                        let scores = VectorizedAlignScores {
                            match_score: 5,
                            mismatch: -4,
                            gap_p: -8,
                            end_gap_p,
                            band,
                        };
                        let mut a = AlignBuffers::new();
                        let mut b = AlignBuffers::new();
                        align_vectorized_with_buf(&s1, &s2, &scores, &mut a);
                        align_vectorized_reference(&s1, &s2, &scores, &mut b);
                        assert_eq!(
                            (&a.al0, &a.al1),
                            (&b.al0, &b.al1),
                            "rolling d16 diverged from the full-matrix reference: \
                             l1={l1} l2={l2} band={band} end_gap_p={end_gap_p} trial={trial}\n\
                             s1={s1:?}\ns2={s2:?}"
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 5000, "fuzz coverage too thin: {checked}");
    }

    /// Buffers are reused across calls in production (`map_init` hands each
    /// rayon worker one `AlignBuffers`). A rolling slot therefore carries stale
    /// values from a *previous, longer* alignment — the case a fresh-buffer
    /// test cannot reach. Replay a descending-length sequence through a single
    /// reused buffer and require the same output as fresh buffers.
    #[test]
    fn rolling_d16_survives_buffer_reuse() {
        let mut state = 0x9E3779B9_7F4A7C15u64;
        let scores = VectorizedAlignScores {
            match_score: 5,
            mismatch: -4,
            gap_p: -8,
            end_gap_p: 0,
            band: 16,
        };
        let mut reused = AlignBuffers::new();
        for &len in &[400usize, 240, 239, 120, 33, 32, 7, 250, 3, 260] {
            let s1: Vec<u8> = (0..len)
                .map(|_| 1 + (rng_next(&mut state) % 4) as u8)
                .collect();
            let mut s2 = s1.clone();
            for k in (0..len).step_by(7) {
                s2[k] = 1 + (s2[k] % 4);
            }
            let mut fresh = AlignBuffers::new();
            align_vectorized_with_buf(&s1, &s2, &scores, &mut fresh);
            align_vectorized_with_buf(&s1, &s2, &scores, &mut reused);
            assert_eq!(
                (&reused.al0, &reused.al1),
                (&fresh.al0, &fresh.al1),
                "reused buffer diverged from a fresh one at len={len}"
            );
        }
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
            backend: AlignBackend::Nw,
            wfa_max_edits: 0,
            match_score: 5,
            mismatch: -4,
            gap_p: -8,
            homo_gap_p: -1,
            use_kmers: true,
            kdist_cutoff: 0.42,
            screen_backend: ScreenBackend::Kmer,
            minimizer_k: crate::minimizers::MINIMIZER_K,
            minimizer_w: crate::minimizers::MINIMIZER_W,
            screen_audit: false,
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
