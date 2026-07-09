//! WFA backend (experimental) — wavefront alignment via wfa2lib-rs.
//!
//! Wraps the pure-Rust WFA aligner (HPCBio fork of COMBINE-lab/wfa2lib-rs) so it
//! can stand in for [`crate::nwalign::align_endsfree`]. WFA *minimises a penalty*
//! (cost) where the match cost is ≤ 0, whereas DADA2 *maximises a score*
//! (match = +5, mismatch = −4, gap = −8). The two are equivalent under sign
//! inversion (score → cost):
//!
//! ```text
//!     match_        = -match_score   (≤ 0, e.g. -5)
//!     mismatch      = -mismatch      (> 0, e.g.  4)
//!     gap_extension = -gap_p         (> 0, e.g.  8)   gap_opening = 0  (linear gap)
//! ```
//!
//! wfa2lib-rs's `new_affine` applies the Eizenga adjustment internally when
//! `match_ < 0`, so a non-zero match score reproduces the same optimum as the
//! scalar DP. Sequences are passed as the raw 1..=5 nt encoding: WFA matches by
//! byte equality and its EOS sentinels (`b'!'`=33, `b'?'`=63) never collide with
//! 1..=5, so no ASCII round-trip is needed.
//!
//! The experimental backend is compiled only under `--features wfa` (it pulls in
//! the git `wfa2lib-rs` dependency, which cannot ship on crates.io). In a default
//! build [`align_wfa_endsfree_with_buf`] is a stub that aborts with a clear
//! message; [`wfa_cost_cap`] is always available since it is pure arithmetic.

use crate::nwalign::AlignBuffers;
#[cfg(feature = "wfa")]
use crate::nwalign::{VectorizedAlignScores, align_vectorized_with_buf};

#[cfg(feature = "wfa")]
use std::cell::RefCell;
#[cfg(feature = "wfa")]
use std::sync::LazyLock;
#[cfg(feature = "wfa")]
use wfa2lib_rs::aligner::{AffineAligner, AlignStatus, AlignmentScope};
#[cfg(feature = "wfa")]
use wfa2lib_rs::heuristic::HeuristicStrategy;
#[cfg(feature = "wfa")]
use wfa2lib_rs::penalties::{AffinePenalties, WavefrontPenalties};

/// Opt-in switch for the experimental WFA alignment backend (issue #49).
/// Set `DADA2RS_ALIGN_BACKEND=wfa` to route the ends-free path through WFA.
/// Read once at first alignment.
#[cfg(feature = "wfa")]
pub(crate) static USE_WFA_BACKEND: LazyLock<bool> = LazyLock::new(|| {
    std::env::var("DADA2RS_ALIGN_BACKEND")
        .map(|v| v.eq_ignore_ascii_case("wfa"))
        .unwrap_or(false)
});

/// Experimental WFA edit-budget cap (issue #51). When set to a positive value,
/// WFA aborts as soon as the alignment *cost* (WFA score, in penalty units)
/// exceeds this bound, and the pair falls back to the banded NW path. This
/// exploits the DADA2 denoising assumption that real error-copies are extremely
/// similar (~99.9% identity): clearly-divergent pairs that survive the k-mer
/// screen otherwise pay WFA's full O(n·s) extend cost, which dominates the
/// PacBio slowdown. The fallback keeps results byte-identical to NW for exactly
/// those capped pairs while WFA's fast path serves the similar majority.
///
/// `0` / unset disables the cap (`i32::MAX`). The value is a WFA *cost*, not an
/// edit count: with default scoring an indel costs `-gap_p` (8) and a mismatch
/// `-mismatch` (4), so e.g. a 40-edit budget ≈ 320.
#[cfg(feature = "wfa")]
static WFA_MAX_STEPS: LazyLock<i32> = LazyLock::new(|| {
    std::env::var("DADA2RS_WFA_MAX_STEPS")
        .ok()
        .and_then(|v| v.parse::<i32>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(i32::MAX)
});

/// Convert an edit-operation budget into a WFA *cost* cap (the unit
/// `set_max_alignment_steps` compares against). The worst case for `E` edits is
/// all of them being the most expensive op — a gap, cost `|gap_p|` — so an
/// alignment needing ≤ `E` edits always has cost ≤ `E·|gap_p|`. Capping cost at
/// that bound therefore never truncates an alignment within budget, regardless
/// of the scoring scheme. `max_edits <= 0` means unbounded.
#[inline]
pub fn wfa_cost_cap(max_edits: i32, gap_p: i32) -> i32 {
    if max_edits > 0 {
        max_edits.saturating_mul(gap_p.abs())
    } else {
        i32::MAX
    }
}

/// Convert a WFA CIGAR (`M`/`X`/`I`/`D` ops, pattern=s1, text=s2) into the
/// gap-annotated `[al0, al1]` pair the rest of the aligner consumes.
///
/// Op semantics (from `Cigar::check_alignment`): `M`/`X` consume one base from
/// each strand; `I` advances the text only (gap in pattern → `al0`); `D`
/// advances the pattern only (gap in text → `al1`).
#[cfg(feature = "wfa")]
pub(crate) fn cigar_to_alignment_into(
    ops: &[u8],
    s1: &[u8],
    s2: &[u8],
    al0: &mut Vec<u8>,
    al1: &mut Vec<u8>,
) {
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
/// with the same alignment-pair representation as
/// [`crate::nwalign::align_endsfree`].
///
/// `match_score`/`mismatch`/`gap_p` use the DADA2 score convention (match > 0,
/// the rest < 0); they are converted to WFA penalties internally. `band` is the
/// DADA2 banding radius (matches the NW path): a `BandedStatic ±band` heuristic
/// is applied so WFA explores the same diagonal band NW does, which both matches
/// NW's band-limited semantics and bounds WFA's O(n·s) cost on divergent pairs
/// (issue #51). `band < 0` means unbanded.
///
/// A per-thread `AffineAligner` is cached and reused across calls so wfa2lib-rs
/// reaches its zero-alloc steady state; it is rebuilt only when the scoring or
/// band changes (both constant within a run).
#[cfg(feature = "wfa")]
#[allow(clippy::too_many_arguments)]
pub fn align_wfa_endsfree_with_buf(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
    max_steps: i32,
    buf: &mut AlignBuffers,
) {
    // A `thread_local` (rather than a field on `AlignBuffers`) keeps the aligner
    // off any struct that crosses threads: `WavefrontAligner` holds raw pointers
    // (`!Send`), and `AlignBuffers` is carried in a rayon `fold` accumulator that
    // must be `Send`. The aligner never moves between threads, so TLS is the
    // right home.
    //
    // Cache key is (match, mismatch, gap_p, band) — all constant within a run.
    type WfaCache = ((i32, i32, i32, i32), AffineAligner);
    thread_local! {
        static WFA: RefCell<Option<WfaCache>> = const { RefCell::new(None) };
    }
    WFA.with(|cell| {
        let mut slot = cell.borrow_mut();
        let key = (match_score, mismatch, gap_p, band);
        if slot.as_ref().map(|(k, _)| *k) != Some(key) {
            let penalties = WavefrontPenalties::new_affine(AffinePenalties {
                match_: -match_score,
                mismatch: -mismatch,
                gap_opening: 0,
                gap_extension: -gap_p,
            });
            let mut aligner = AffineAligner::new(penalties);
            // Compute the full CIGAR, not just the score (default ComputeScore
            // leaves the cigar empty).
            aligner.alignment_scope = AlignmentScope::ComputeAlignment;
            // Match NW's diagonal band so WFA stays band-limited and its O(n·s)
            // cost can't blow up on divergent-but-screened pairs (issue #51).
            if band >= 0 {
                aligner.set_heuristic(HeuristicStrategy::BandedStatic {
                    min_k: -band,
                    max_k: band,
                });
            }
            *slot = Some((key, aligner));
        }
        let aligner = &mut slot.as_mut().unwrap().1;
        // Edit-budget cap (issue #51): abort WFA once cost exceeds the bound so
        // divergent pairs don't pay the full O(n·s) extend. Capped pairs fall
        // back to NW below. `max_steps` is the caller's cost cap (from
        // `wfa_cost_cap`); the `DADA2RS_WFA_MAX_STEPS` env var overrides it for
        // ad-hoc sweeps. i32::MAX = unbounded.
        let cap = if *WFA_MAX_STEPS != i32::MAX {
            *WFA_MAX_STEPS
        } else {
            max_steps
        };
        aligner.set_max_alignment_steps(cap);
        // Full ends-free: leading/trailing gaps on either strand are free.
        aligner.set_alignment_free_ends(
            s1.len() as i32,
            s1.len() as i32,
            s2.len() as i32,
            s2.len() as i32,
        );
        aligner.align_endsfree(s1, s2);
        // On MaxStepsReached the CIGAR is left empty; fall back to the banded NW
        // ends-free path so the pair still gets a valid alignment (byte-identical
        // to the pure-NW backend for exactly these capped, divergent pairs).
        if aligner.status() == AlignStatus::MaxStepsReached {
            align_vectorized_with_buf(
                s1,
                s2,
                &VectorizedAlignScores {
                    match_score: match_score as i16,
                    mismatch: mismatch as i16,
                    gap_p: gap_p as i16,
                    end_gap_p: 0,
                    band,
                },
                buf,
            );
            return;
        }
        let ops = aligner.cigar().operations_slice();
        cigar_to_alignment_into(ops, s1, s2, &mut buf.al0, &mut buf.al1);
    });
}

/// NW-only build stub for the WFA backend (issue #63). The experimental WFA
/// backend is compiled only under `--features wfa` (it pulls in the git
/// `wfa2lib-rs` dependency, which cannot ship on crates.io). Selecting
/// `--align-backend wfa2` in a default build reaches this stub and aborts with a
/// clear message rather than silently falling back to NW.
#[cfg(not(feature = "wfa"))]
#[allow(clippy::too_many_arguments)]
pub fn align_wfa_endsfree_with_buf(
    _s1: &[u8],
    _s2: &[u8],
    _match_score: i32,
    _mismatch: i32,
    _gap_p: i32,
    _band: i32,
    _max_steps: i32,
    _buf: &mut AlignBuffers,
) {
    panic!(
        "the WFA alignment backend (--align-backend wfa2) is not available in \
         this build: dada2-rs was compiled without the experimental `wfa` \
         feature. Rebuild from source with `cargo build --features wfa` to \
         enable it (see issue #63)."
    );
}
