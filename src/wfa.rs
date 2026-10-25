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

use crate::nwalign::{AlignBuffers, EndModel};
#[cfg(feature = "wfa")]
use crate::nwalign::{VectorizedAlignScores, align_vectorized_with_buf, band_adjust};

#[cfg(feature = "wfa")]
use std::cell::RefCell;
#[cfg(feature = "wfa")]
use std::sync::LazyLock;
#[cfg(feature = "wfa")]
use wfa2lib_rs::aligner::{
    Affine2pAligner, AffineAligner, AlignStatus, AlignmentScope, Compat, WavefrontAligner,
};
#[cfg(feature = "wfa")]
use wfa2lib_rs::heuristic::HeuristicStrategy;
#[cfg(feature = "wfa")]
use wfa2lib_rs::penalties::{
    Affine2pPenalties, AffinePenalties, LinearPenalties, WavefrontPenalties,
};

/// `DADA2RS_WFA_LINEAR=1` runs WFA with gap-linear penalties (one wavefront)
/// instead of affine with gap opening 0 (three). Same scores and, on the
/// fixtures, byte-identical error models; about 2x faster uncapped with half
/// the peak memory. An A/B switch until confirmed at scale (#51).
#[cfg(feature = "wfa")]
static WFA_LINEAR: LazyLock<bool> =
    LazyLock::new(|| std::env::var("DADA2RS_WFA_LINEAR").is_ok_and(|v| v == "1"));

/// Alignment span for the alignment-model experiment (issue #309); `None`
/// follows `--end-model`. `s1` is always the cluster centre and `s2` the read.
#[cfg(feature = "wfa")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum WfaSpan {
    /// Global: every end gap is charged.
    EndToEnd,
    /// Free end gaps on both sequences (DADA2's model).
    EndsFree,
    /// Free trailing gaps on both sequences (`--end-model free-3p`).
    Free3p,
    /// Free end gaps on the centre only: the read is aligned end to end,
    /// anywhere within the centre (semi-global / glocal).
    CentreFree,
    /// Free end gaps on the read only: the centre is aligned end to end.
    ReadFree,
    /// Extension: anchored at the start, trimmed to the best-scoring prefix.
    Extension,
}

/// Experimental alignment-model settings (issue #309), from developer-only
/// environment overrides read once per run. All unset is today's model:
/// span from `--end-model`, linear gaps (affine with gap opening 0).
///
/// - `DADA2RS_WFA_SPAN`: `e2e`, `endsfree`, `centre-free`, `read-free`,
///   `extension`
/// - `DADA2RS_WFA_GAP_OPEN=N`: gap-affine opening cost (WFA penalty units)
/// - `DADA2RS_WFA_AFFINE2P=o1,e1,o2,e2`: two-piece gap-affine
#[cfg(feature = "wfa")]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct WfaExperiment {
    pub span: Option<WfaSpan>,
    pub gap_open: i32,
    pub affine2p: Option<[i32; 4]>,
}

#[cfg(feature = "wfa")]
impl WfaExperiment {
    fn from_env() -> Self {
        let span = std::env::var("DADA2RS_WFA_SPAN")
            .ok()
            .map(|v| match v.as_str() {
                "e2e" => WfaSpan::EndToEnd,
                "endsfree" => WfaSpan::EndsFree,
                "free-3p" => WfaSpan::Free3p,
                "centre-free" => WfaSpan::CentreFree,
                "read-free" => WfaSpan::ReadFree,
                "extension" => WfaSpan::Extension,
                other => panic!("DADA2RS_WFA_SPAN: unknown span `{other}`"),
            });
        let gap_open = std::env::var("DADA2RS_WFA_GAP_OPEN")
            .ok()
            .map(|v| {
                v.parse::<i32>()
                    .expect("DADA2RS_WFA_GAP_OPEN: integer >= 0")
            })
            .unwrap_or(0);
        let affine2p = std::env::var("DADA2RS_WFA_AFFINE2P").ok().map(|v| {
            let x: Vec<i32> = v
                .split(',')
                .map(|n| n.trim().parse().expect("DADA2RS_WFA_AFFINE2P: o1,e1,o2,e2"))
                .collect();
            assert_eq!(x.len(), 4, "DADA2RS_WFA_AFFINE2P: o1,e1,o2,e2");
            [x[0], x[1], x[2], x[3]]
        });
        let e = WfaExperiment {
            span,
            gap_open,
            affine2p,
        };
        if e != WfaExperiment::default() {
            assert!(
                !*WFA_LINEAR || (e.gap_open == 0 && e.affine2p.is_none()),
                "DADA2RS_WFA_LINEAR cannot use affine gap penalties"
            );
            eprintln!("[wfa] alignment-model experiment (#309): {e:?}");
        }
        e
    }
}

#[cfg(feature = "wfa")]
static WFA_EXPERIMENT: LazyLock<WfaExperiment> = LazyLock::new(WfaExperiment::from_env);

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
    // Extension can stop before either sequence ends; the unaligned tails are
    // end overhangs, as unscored as any other terminal gap.
    for &b in &s1[p..] {
        al0.push(b);
        al1.push(b'-');
    }
    for &b in &s2[t..] {
        al0.push(b'-');
        al1.push(b);
    }
}

/// Ends-free Needleman-Wunsch via the WFA backend. Fills `buf.al0`/`buf.al1`
/// with the same alignment-pair representation as
/// [`crate::nwalign::align_endsfree`].
///
/// `match_score`/`mismatch`/`gap_p` use the DADA2 score convention (match > 0,
/// the rest < 0); they are converted to WFA penalties internally. `band` is the
/// DADA2 banding radius (matches the NW path): a `BandedStatic` heuristic over
/// NW's own diagonal range (`band_adjust`, widened by the length difference) is
/// applied per pair, so WFA explores the same band NW does, which both matches
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
    ends: EndModel,
    max_steps: i32,
    buf: &mut AlignBuffers,
) {
    align_wfa_with_experiment(
        s1,
        s2,
        match_score,
        mismatch,
        gap_p,
        band,
        ends,
        max_steps,
        &WFA_EXPERIMENT,
        buf,
    );
}

/// [`align_wfa_endsfree_with_buf`] with explicit experiment settings (#309),
/// so tests need no environment.
#[cfg(feature = "wfa")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn align_wfa_with_experiment(
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
    ends: EndModel,
    max_steps: i32,
    exp: &WfaExperiment,
    buf: &mut AlignBuffers,
) {
    // A `thread_local` (rather than a field on `AlignBuffers`) keeps the aligner
    // off any struct that crosses threads: `WavefrontAligner` holds raw pointers
    // (`!Send`), and `AlignBuffers` is carried in a rayon `fold` accumulator that
    // must be `Send`. The aligner never moves between threads, so TLS is the
    // right home.
    //
    // Cache key: the scoring, band and experiment penalties, all constant
    // within a run.
    type Key = (i32, i32, i32, i32, i32, Option<[i32; 4]>);
    thread_local! {
        static WFA: RefCell<Option<(Key, AffineAligner)>> = const { RefCell::new(None) };
        static WFA_2P: RefCell<Option<(Key, Affine2pAligner)>> = const { RefCell::new(None) };
        static WFA_LIN: RefCell<Option<(Key, WavefrontAligner<1>)>> =
            const { RefCell::new(None) };
    }
    let key = (
        match_score,
        mismatch,
        gap_p,
        band,
        exp.gap_open,
        exp.affine2p,
    );
    if let Some([o1, e1, o2, e2]) = exp.affine2p {
        WFA_2P.with(|cell| {
            let mut slot = cell.borrow_mut();
            if slot.as_ref().map(|(k, _)| *k) != Some(key) {
                let penalties = WavefrontPenalties::new_affine2p(Affine2pPenalties {
                    match_: -match_score,
                    mismatch: -mismatch,
                    gap_opening1: o1,
                    gap_extension1: e1,
                    gap_opening2: o2,
                    gap_extension2: e2,
                });
                *slot = Some((key, configure(Affine2pAligner::new(penalties))));
            }
            let aligner = &mut slot.as_mut().unwrap().1;
            align_with(
                aligner,
                s1,
                s2,
                match_score,
                mismatch,
                gap_p,
                band,
                ends,
                max_steps,
                exp.span,
                buf,
            );
        });
        return;
    }
    if *WFA_LINEAR {
        WFA_LIN.with(|cell| {
            let mut slot = cell.borrow_mut();
            if slot.as_ref().map(|(k, _)| *k) != Some(key) {
                let penalties = WavefrontPenalties::new_linear(LinearPenalties {
                    match_: -match_score,
                    mismatch: -mismatch,
                    indel: -gap_p,
                });
                *slot = Some((key, configure(WavefrontAligner::<1>::new(penalties))));
            }
            let aligner = &mut slot.as_mut().unwrap().1;
            align_with(
                aligner,
                s1,
                s2,
                match_score,
                mismatch,
                gap_p,
                band,
                ends,
                max_steps,
                exp.span,
                buf,
            );
        });
        return;
    }
    WFA.with(|cell| {
        let mut slot = cell.borrow_mut();
        if slot.as_ref().map(|(k, _)| *k) != Some(key) {
            let penalties = WavefrontPenalties::new_affine(AffinePenalties {
                match_: -match_score,
                mismatch: -mismatch,
                gap_opening: exp.gap_open,
                gap_extension: -gap_p,
            });
            *slot = Some((key, configure(AffineAligner::new(penalties))));
        }
        let aligner = &mut slot.as_mut().unwrap().1;
        align_with(
            aligner,
            s1,
            s2,
            match_score,
            mismatch,
            gap_p,
            band,
            ends,
            max_steps,
            exp.span,
            buf,
        );
    });
}

/// One-time setup shared by the affine and linear aligners.
#[cfg(feature = "wfa")]
fn configure<const N: usize>(mut aligner: WavefrontAligner<N>) -> WavefrontAligner<N> {
    // Compute the full CIGAR, not just the score (default ComputeScore
    // leaves the cigar empty).
    aligner.alignment_scope = AlignmentScope::ComputeAlignment;
    // Fix the C WFA2-lib bugs the fork reproduces by default
    // (HPCBio/wfa2lib-rs#32), notably ends-free with a match reward
    // stopping at the first free end rather than the best one
    // (Marco-Sola-Lab/WFA2-lib#102, HPCBio/wfa2lib-rs#37 and #38).
    // dada2-rs needs NW's optimum, not C-identical output.
    aligner.set_compat(Compat::Corrected);
    aligner
}

/// Per-pair band, cap and free ends, then align (or fall back to NW).
#[cfg(feature = "wfa")]
#[allow(clippy::too_many_arguments)]
fn align_with<const N: usize>(
    aligner: &mut WavefrontAligner<N>,
    s1: &[u8],
    s2: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_p: i32,
    band: i32,
    ends: EndModel,
    max_steps: i32,
    span: Option<WfaSpan>,
    buf: &mut AlignBuffers,
) {
    // Match NW's diagonal band so WFA stays band-limited and its O(n·s)
    // cost can't blow up on divergent-but-screened pairs (issue #51). The
    // band is per pair: NW widens one side by the length difference, and a
    // symmetric ±band cut WFA off from leading gaps NW can take. A WFA
    // diagonal is text offset minus pattern offset, NW's `j - i`.
    if band >= 0 {
        let (lband, rband) = band_adjust(s1.len(), s2.len(), band);
        aligner.set_heuristic(HeuristicStrategy::BandedStatic {
            min_k: -lband,
            max_k: rband,
        });
    }
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
    // Ends-free, with each free end bounded by what the banded NW path can
    // reach: `band + |len1 - len2|` (see `band_adjust`). Whole-length free
    // ends let WFA (C and the port alike) stop on the no-overlap alignment
    // (all end gaps, score 0) for divergent pairs once begin-free gaps are
    // honored (HPCBio/wfa2lib-rs#17); on PacBio this collapsed the
    // unbounded arm's ASVs (#51). Unbanded keeps whole-length free ends.
    let (p_free, t_free) = if band >= 0 {
        let fe = band + (s1.len() as i32 - s2.len() as i32).abs();
        (fe.min(s1.len() as i32), fe.min(s2.len() as i32))
    } else {
        (s1.len() as i32, s2.len() as i32)
    };
    // Penalised ends (issue #289): global is end-to-end; free-3p frees
    // trailing gaps only. With a free end the match reward changes which
    // alignment is best, so free-3p relies on the corrected mode and the
    // real reward (HPCBio/wfa2lib-rs#31); global ranks alike either way.
    let span = span.unwrap_or(match ends {
        EndModel::Global => WfaSpan::EndToEnd,
        EndModel::EndsFree => WfaSpan::EndsFree,
        EndModel::Free3p => WfaSpan::Free3p,
    });
    match span {
        WfaSpan::EndToEnd => {
            aligner.set_alignment_end_to_end();
            aligner.align_end2end(s1, s2);
        }
        WfaSpan::EndsFree => {
            aligner.set_alignment_free_ends(p_free, p_free, t_free, t_free);
            aligner.align_endsfree(s1, s2);
        }
        WfaSpan::Free3p => {
            aligner.set_alignment_free_ends(0, p_free, 0, t_free);
            aligner.align_endsfree(s1, s2);
        }
        WfaSpan::CentreFree => {
            aligner.set_alignment_free_ends(p_free, p_free, 0, 0);
            aligner.align_endsfree(s1, s2);
        }
        WfaSpan::ReadFree => {
            aligner.set_alignment_free_ends(0, 0, t_free, t_free);
            aligner.align_endsfree(s1, s2);
        }
        WfaSpan::Extension => {
            aligner.set_alignment_extension();
            aligner.align_endsfree(s1, s2);
        }
    }
    // On MaxStepsReached the CIGAR is left empty; fall back to the banded NW
    // path under the same end model, so the pair still gets a valid alignment
    // (byte-identical to the pure-NW backend for exactly these capped pairs).
    if aligner.status() == AlignStatus::MaxStepsReached {
        let (begin_free, end_free) = ends.free_ends();
        let edge = |free: bool| if free { 0 } else { gap_p as i16 };
        align_vectorized_with_buf(
            s1,
            s2,
            &VectorizedAlignScores {
                match_score: match_score as i16,
                mismatch: mismatch as i16,
                gap_p: gap_p as i16,
                begin_gap_p: edge(begin_free),
                end_gap_p: edge(end_free),
                band,
            },
            buf,
        );
        return;
    }
    let ops = aligner.cigar().operations_slice();
    cigar_to_alignment_into(ops, s1, s2, &mut buf.al0, &mut buf.al1);
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
    _ends: EndModel,
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

#[cfg(all(test, feature = "wfa"))]
mod experiment_tests {
    use super::*;

    /// A deterministic pseudo-random sequence in the 1..=4 encoding.
    fn seq(seed: u64, len: usize) -> Vec<u8> {
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                1 + ((s >> 33) % 4) as u8
            })
            .collect()
    }

    fn align(s1: &[u8], s2: &[u8], exp: WfaExperiment) -> (Vec<u8>, Vec<u8>) {
        let mut buf = AlignBuffers::new();
        align_wfa_with_experiment(
            s1,
            s2,
            5,
            -4,
            -8,
            32,
            EndModel::EndsFree,
            i32::MAX,
            &exp,
            &mut buf,
        );
        (buf.al0.clone(), buf.al1.clone())
    }

    fn ungapped(al: &[u8]) -> Vec<u8> {
        al.iter().copied().filter(|&b| b != b'-').collect()
    }

    /// A read missing the centre's first 120 bases: centre-free and
    /// end-to-end with affine gaps both leave one leading gap (#309).
    #[test]
    fn truncated_read_aligns_as_one_leading_gap() {
        let centre = seq(7, 400);
        let read = centre[120..].to_vec();
        for exp in [
            WfaExperiment {
                span: Some(WfaSpan::CentreFree),
                ..Default::default()
            },
            WfaExperiment {
                span: Some(WfaSpan::EndToEnd),
                gap_open: 8,
                ..Default::default()
            },
            WfaExperiment {
                span: Some(WfaSpan::EndToEnd),
                affine2p: Some([4, 8, 24, 2]),
                ..Default::default()
            },
        ] {
            let (al0, al1) = align(&centre, &read, exp);
            assert_eq!(ungapped(&al0), centre, "{exp:?}");
            assert!(
                al1[..120].iter().all(|&b| b == b'-'),
                "{exp:?}: gap not leading"
            );
            assert!(
                al1[120..].iter().all(|&b| b != b'-'),
                "{exp:?}: interior gaps"
            );
        }
    }

    /// Extension may stop before the sequences end; the tails are padded as
    /// end overhangs, so both rows still spell the inputs.
    #[test]
    fn extension_pads_unaligned_tails() {
        let centre = seq(11, 300);
        let mut read = centre[..200].to_vec();
        read.extend(seq(99, 100));
        let exp = WfaExperiment {
            span: Some(WfaSpan::Extension),
            ..Default::default()
        };
        let (al0, al1) = align(&centre, &read, exp);
        assert_eq!(ungapped(&al0), centre);
        assert_eq!(ungapped(&al1), read);
        assert_eq!(al0.len(), al1.len());
    }
}
