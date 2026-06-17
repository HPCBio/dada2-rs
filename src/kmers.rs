//! K-mer frequency and distance functions.
//!
//! Ported from `kmers.cpp`. The original C++ contains SSE2 SIMD paths
//! (`kmer_dist_SSEi`, `kmer_dist_SSEi_8`, `kord_dist_SSEi`). The scalar
//! loops here are written in a form that LLVM auto-vectorises when compiled
//! with `-C target-cpu=native`, providing equivalent throughput without
//! manual intrinsics.
//!
//! Nucleotide encoding: A=1, C=2, G=3, T=4 (matching `misc.rs`).
//! K-mer indices are computed in base-4 with A→0, C→1, G→2, T→3,
//! i.e. `index = sum over positions of (nt - 1) * 4^(k-1-i)`.

use crate::containers::Raw;
use std::cmp::Ordering;

/// Default k-mer size, matching the C++ `KMER_SIZE` constant in `dada.h`.
/// At runtime this is overridable via `AlignParams::kmer_size`; the const is
/// kept only so callers without an `AlignParams` (e.g. tests, defaults) can
/// reference it.
#[allow(dead_code)]
pub const KMER_SIZE: usize = 5;

/// Valid k-mer size range. `assign_kmer` requires k ≥ 3 (smaller k produces
/// uselessly coarse vectors) and k ≤ 8 (4^k ≤ 65536 fits in our u16 vectors
/// and matches the C++ assumption).
pub const KMER_SIZE_MIN: usize = 3;
pub const KMER_SIZE_MAX: usize = 8;

/// Smallest k for which the resident 8-bit k-mer screen is stored *sparsely*
/// (issue #43). The sparse `(index, count)` representation costs a fixed
/// ~6 KB/raw (≈seqlen entries), independent of k, whereas the dense `4^k` array
/// grows as `4^k`. The PacBio 93-sample pooled A/B (dense vs sparse, per k)
/// pinned the crossover precisely at where `4^k` passes ~6 KB:
///   - k5 (1 KB) / k6 (4 KB): dense is *smaller* than sparse, and its contiguous
///     `sum-of-min` sweep is cache-resident + auto-vectorised, so dense wins both
///     memory and wall (k6 sparse regressed +20.8% RSS *and* +20.1% dada wall).
///   - k7 (16 KB): sparse wins RSS (−22.6%) but 16 KB still fits cache, so the
///     dense sweep beats the branchy merge-join → +25.7% dada wall. A tradeoff,
///     and k7 is the production default — so we keep it dense.
///   - k8 (64 KB): dense spills cache; merge-join wins *both* (−68.3% RSS,
///     −17.3% wall) and makes k8 feasible at all (42.8 GB dense → 13.6 GB).
///
/// So sparse is gated to k≥8, where it is a pure win with no regression at any
/// smaller k; k≤7 keeps the byte-identical dense path.
pub const SPARSE_KMER_MIN: usize = 8;

/// Number of possible k-mers for a given k: 4^k = 1 << (2*k).
#[inline]
pub fn n_kmers(k: usize) -> usize {
    1 << (2 * k)
}

/// Encode a window of integer-encoded nucleotides to a k-mer index.
/// Returns `None` if any nucleotide value is not 1..=4 (A/C/G/T).
#[inline]
fn encode_kmer(window: &[u8]) -> Option<usize> {
    let mut kmer = 0usize;
    for &nt in window {
        let nti = nt as usize;
        if !(1..=4).contains(&nti) {
            return None;
        }
        kmer = 4 * kmer + (nti - 1);
    }
    Some(kmer)
}

// ---------------------------------------------------------------------------
// Assignment functions
// ---------------------------------------------------------------------------

/// Build a 16-bit k-mer frequency vector from an integer-encoded sequence.
///
/// `k` must be 3..=8 and strictly less than `seq.len()`.
/// Non-ACGT positions (values outside 1..=4) are silently skipped, matching
/// the C++ behaviour where `kmer = 999999` causes the k-mer to be ignored.
/// Counts saturate at `u16::MAX` via `saturating_add`.
/// Equivalent to C++ `assign_kmer`.
pub fn assign_kmer(seq: &[u8], k: usize) -> Vec<u16> {
    debug_assert!((3..=8).contains(&k), "k must be 3..=8");
    debug_assert!(seq.len() > k, "sequence must be longer than k");
    let mut kvec = vec![0u16; n_kmers(k)];
    for window in seq.windows(k) {
        if let Some(idx) = encode_kmer(window) {
            kvec[idx] = kvec[idx].saturating_add(1);
        }
    }
    kvec
}

/// Build an 8-bit k-mer frequency vector, saturating at 255.
///
/// Computed by downconverting the 16-bit vector, matching the C++
/// implementation that first builds a `uint16_t` vector then truncates.
/// Equivalent to C++ `assign_kmer8`.
pub fn assign_kmer8(seq: &[u8], k: usize) -> Vec<u8> {
    assign_kmer(seq, k)
        .into_iter()
        .map(|v| v.min(255) as u8)
        .collect()
}

/// Build a *sparse* 8-bit k-mer screen: sorted `(kmer_index, count)` pairs for
/// the k-mers actually present in `seq`, counts saturating at 255 (issue #43).
///
/// The index fits `u16` for all supported k (`4^8 = 65536` ≤ `u16::MAX + 1`,
/// the largest index being `4^k - 1 = 65535` at k=8). Entries are emitted in
/// ascending index order (a dense scan), which the merge-join in
/// [`kmer_dist8_sparse`] relies on. Built via the dense `u16` vector then
/// filtered; the dense intermediate is transient (freed per call, bounded by
/// thread count), while the resident result is ~`seqlen` entries instead of
/// `4^k` bytes — the whole point at large k.
pub fn assign_kmer8_sparse(seq: &[u8], k: usize) -> Vec<(u16, u8)> {
    assign_kmer(seq, k)
        .into_iter()
        .enumerate()
        .filter(|&(_, c)| c > 0)
        .map(|(idx, c)| (idx as u16, c.min(255) as u8))
        .collect()
}

/// Build the ordered k-mer vector: the k-mer index at each sequence position.
///
/// `kord[i]` is the index of the k-mer starting at position `i`. Positions
/// containing non-ACGT nucleotides produce index `0`, matching the C++
/// behaviour where the slot is initialised to 0 and never overwritten for
/// invalid k-mers.
/// Equivalent to C++ `assign_kmer_order`.
pub fn assign_kmer_order(seq: &[u8], k: usize) -> Vec<u16> {
    debug_assert!((1..=8).contains(&k), "k must be 1..=8");
    debug_assert!(seq.len() > k, "sequence must be longer than k");
    seq.windows(k)
        .map(|w| encode_kmer(w).unwrap_or(0) as u16)
        .collect()
}

/// Resident 8-bit k-mer screen representation (issue #43).
///
/// `Dense` is the classic `4^k`-byte frequency array (used for k < [`SPARSE_KMER_MIN`]);
/// `Sparse` is sorted `(index, count)` pairs for present k-mers only (k ≥
/// `SPARSE_KMER_MIN`). Both yield identical [`kmer_dist8`] results — `Sparse`
/// merely skips the buckets where the dense `min` would be zero.
#[derive(Clone)]
pub enum KmerScreen {
    Dense(Vec<u8>),
    Sparse(Vec<(u16, u8)>),
}

impl KmerScreen {
    /// 8-bit k-mer distance against another screen of the same representation.
    ///
    /// Returns `-1.0` on saturation overflow (a k-mer at ≥255 in *both* seqs),
    /// matching [`kmer_dist8`]. Both raws in a run share `k`, hence the same
    /// variant; mixed variants are a programming error.
    #[inline]
    pub fn dist8(&self, other: &KmerScreen, len1: usize, len2: usize, k: usize) -> f64 {
        match (self, other) {
            (KmerScreen::Dense(a), KmerScreen::Dense(b)) => kmer_dist8(a, len1, b, len2, k),
            (KmerScreen::Sparse(a), KmerScreen::Sparse(b)) => {
                kmer_dist8_sparse(a, len1, b, len2, k)
            }
            _ => unreachable!("mixed dense/sparse k-mer screens (k must be uniform per run)"),
        }
    }

    /// Resident byte footprint of this screen (for `--verbose` accounting).
    #[inline]
    pub fn resident_bytes(&self) -> usize {
        match self {
            KmerScreen::Dense(v) => v.len(),
            KmerScreen::Sparse(v) => v.len() * std::mem::size_of::<(u16, u8)>(),
        }
    }

    /// Number of *distinct* k-mers present in this sequence (sparse entries, or
    /// dense nonzero buckets). Bounded by `len - k + 1`; below that signals
    /// low-complexity/repetitive content. For `--verbose` diagnostics (#43).
    #[inline]
    pub fn distinct_kmers(&self) -> usize {
        match self {
            KmerScreen::Dense(v) => v.iter().filter(|&&c| c > 0).count(),
            KmerScreen::Sparse(v) => v.len(),
        }
    }

    /// Invoke `f` with each present k-mer index, for tallying the pooled union
    /// of distinct k-mers across all uniques (`--verbose` diagnostics, #43).
    #[inline]
    pub fn for_each_present_index<F: FnMut(usize)>(&self, mut f: F) {
        match self {
            KmerScreen::Dense(v) => {
                for (i, &c) in v.iter().enumerate() {
                    if c > 0 {
                        f(i);
                    }
                }
            }
            KmerScreen::Sparse(v) => {
                for &(idx, _) in v {
                    f(idx as usize);
                }
            }
        }
    }
}

/// Populate the resident k-mer screen fields (`kmer8`, `kord`) on a `Raw`.
///
/// The exact 16-bit frequency vector is intentionally NOT stored (issue #32):
/// it dominated pooled RSS at k7 and is only consulted on the `kmer_dist8`
/// overflow fallback, where it is recomputed from `seq` via [`assign_kmer`].
/// The screen is stored sparsely for k ≥ [`SPARSE_KMER_MIN`] (issue #43).
pub fn raw_assign_kmers(raw: &mut Raw, k: usize) {
    let screen = if k >= SPARSE_KMER_MIN {
        KmerScreen::Sparse(assign_kmer8_sparse(&raw.seq, k))
    } else {
        KmerScreen::Dense(assign_kmer8(&raw.seq, k))
    };
    raw.kmer8 = Some(screen);
    raw.kord = Some(assign_kmer_order(&raw.seq, k));
}

// ---------------------------------------------------------------------------
// Distance functions
// ---------------------------------------------------------------------------

/// K-mer frequency-vector distance between two sequences.
///
/// `dist = 1 - dotsum / (min(len1, len2) - k + 1)`
/// where `dotsum = Σ min(kv1[i], kv2[i])`.
///
/// Equivalent to C++ `kmer_dist` and `kmer_dist_SSEi`.
pub fn kmer_dist(kv1: &[u16], len1: usize, kv2: &[u16], len2: usize, k: usize) -> f64 {
    let dotsum: u32 = kv1
        .iter()
        .zip(kv2.iter())
        .map(|(&a, &b)| a.min(b) as u32)
        .sum();
    let scale = (len1.min(len2) - k + 1) as f64;
    1.0 - dotsum as f64 / scale
}

/// K-mer frequency-vector distance using 8-bit vectors.
///
/// Returns `-1.0` if any element of the element-wise minimum equals 255,
/// indicating saturation overflow and an unreliable result.
/// Equivalent to C++ `kmer_dist_SSEi_8`.
///
/// LLVM auto-vectorises this to NEON `umin.16b` on aarch64 (loop processes
/// 32 B/iter with 4 parallel u32 accumulators) and equivalent SSE/AVX on
/// x86-64. Measured at ~33–39 GB/s (memory-bandwidth-bound) for k≥5 on
/// Apple Silicon — see the `bench` module in this file. No manual SIMD is
/// needed.
pub fn kmer_dist8(kv1: &[u8], len1: usize, kv2: &[u8], len2: usize, k: usize) -> f64 {
    let mut dotsum = 0u32;
    let mut overflow = false;
    for (&a, &b) in kv1.iter().zip(kv2.iter()) {
        let m = a.min(b);
        if m == 255 {
            overflow = true;
        }
        dotsum += m as u32;
    }
    if overflow {
        return -1.0;
    }
    let scale = (len1.min(len2) - k + 1) as f64;
    1.0 - dotsum as f64 / scale
}

/// Sparse counterpart of [`kmer_dist8`] over sorted `(index, count)` pairs
/// (issue #43).
///
/// The dense distance is `1 - Σ min(a[i], b[i]) / scale`; `min` is zero at any
/// bucket where either count is zero, so only k-mers present in *both* sequences
/// contribute. A merge-join over the two ascending-index vectors visits exactly
/// those, giving a result bit-identical to the dense path while touching
/// ~`seqlen` entries instead of `4^k`. Overflow (`min == 255`) is detected on the
/// matching buckets, the only place a saturated `min` can arise.
pub fn kmer_dist8_sparse(
    kv1: &[(u16, u8)],
    len1: usize,
    kv2: &[(u16, u8)],
    len2: usize,
    k: usize,
) -> f64 {
    let mut dotsum = 0u32;
    let (mut i, mut j) = (0usize, 0usize);
    while i < kv1.len() && j < kv2.len() {
        let (ai, ac) = kv1[i];
        let (bi, bc) = kv2[j];
        match ai.cmp(&bi) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let m = ac.min(bc);
                if m == 255 {
                    return -1.0;
                }
                dotsum += m as u32;
                i += 1;
                j += 1;
            }
        }
    }
    let scale = (len1.min(len2) - k + 1) as f64;
    1.0 - dotsum as f64 / scale
}

/// Ordered k-mer distance (fraction of positionally mismatched k-mers).
///
/// Returns `-1.0` if the sequence lengths differ.
/// Equivalent to C++ `kord_dist` and `kord_dist_SSEi`.
pub fn kord_dist(kord1: &[u16], len1: usize, kord2: &[u16], len2: usize, k: usize) -> f64 {
    if len1 != len2 {
        return -1.0;
    }
    let klen = match len1.checked_sub(k.saturating_sub(1)) {
        Some(l) if l > 0 => l,
        _ => return -1.0,
    };
    let dotsum: u32 = kord1[..klen]
        .iter()
        .zip(kord2[..klen].iter())
        .map(|(a, b)| (a == b) as u32)
        .sum();
    1.0 - dotsum as f64 / klen as f64
}

#[cfg(test)]
mod sparse_tests {
    use super::*;

    // Deterministic pseudo-random ACGT (1..=4) sequence.
    fn make_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((s >> 33) % 4) as u8 + 1
            })
            .collect()
    }

    /// Sparse builder must round-trip to the dense vector exactly.
    #[test]
    fn sparse_reconstructs_dense() {
        for &k in &[6usize, 7, 8] {
            for seed in 0..8 {
                let seq = make_seq(300, seed);
                let dense = assign_kmer8(&seq, k);
                let sparse = assign_kmer8_sparse(&seq, k);
                // ascending, unique indices
                assert!(sparse.windows(2).all(|w| w[0].0 < w[1].0));
                let mut reconstructed = vec![0u8; n_kmers(k)];
                for &(idx, c) in &sparse {
                    reconstructed[idx as usize] = c;
                }
                assert_eq!(dense, reconstructed, "k={k} seed={seed}");
            }
        }
    }

    /// Sparse merge-join distance must be bit-identical to the dense sweep.
    #[test]
    fn sparse_dist_matches_dense() {
        for &k in &[6usize, 7, 8] {
            for seed in 0..16 {
                let a = make_seq(250, seed);
                let b = make_seq(250, seed + 100);
                let (da, db) = (assign_kmer8(&a, k), assign_kmer8(&b, k));
                let (sa, sb) = (assign_kmer8_sparse(&a, k), assign_kmer8_sparse(&b, k));
                let dense = kmer_dist8(&da, a.len(), &db, b.len(), k);
                let sparse = kmer_dist8_sparse(&sa, a.len(), &sb, b.len(), k);
                assert_eq!(dense.to_bits(), sparse.to_bits(), "k={k} seed={seed}");
                // identical sequences → distance 0
                let sparse_self = kmer_dist8_sparse(&sa, a.len(), &sa, a.len(), k);
                assert_eq!(
                    sparse_self.to_bits(),
                    0.0f64.to_bits(),
                    "self k={k} seed={seed}"
                );
            }
        }
    }

    /// Overflow (a k-mer at ≥255 in both seqs) must return -1.0 in both paths.
    #[test]
    fn sparse_dist_overflow_matches_dense() {
        let k = 6;
        // a 300-A homopolymer: the all-A k-mer occurs 295× (>255) → saturates
        let seq = vec![1u8; 300];
        let dense = assign_kmer8(&seq, k);
        let sparse = assign_kmer8_sparse(&seq, k);
        let d = kmer_dist8(&dense, seq.len(), &dense, seq.len(), k);
        let s = kmer_dist8_sparse(&sparse, seq.len(), &sparse, seq.len(), k);
        assert_eq!(d, -1.0);
        assert_eq!(s, -1.0);
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    use std::time::Instant;

    fn make_kv(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (s >> 33) as u8 & 0x07 // small values to avoid saturation
            })
            .collect()
    }

    #[test]
    #[ignore] // run explicitly: cargo test --release -- --ignored bench_kmer_dist8 --nocapture
    fn bench_kmer_dist8() {
        for &k in &[3usize, 5, 7, 8] {
            let n = 1 << (2 * k);
            let a = make_kv(n, 11);
            let b = make_kv(n, 22);
            let iters = 2_000_000 / n.max(1); // amortise call count across sizes
            let reps = 100usize;

            // Warmup
            let mut acc = 0.0f64;
            for _ in 0..(iters / 10).max(1) {
                acc += kmer_dist8(&a, 250, &b, 250, k);
            }

            let t0 = Instant::now();
            for _ in 0..reps {
                for _ in 0..iters {
                    acc += kmer_dist8(&a, 250, &b, 250, k);
                }
            }
            let dt = t0.elapsed();
            std::hint::black_box(acc);

            let total_calls = (iters * reps) as f64;
            let ns_per_call = dt.as_nanos() as f64 / total_calls;
            let bytes_per_call = 2.0 * n as f64;
            let gbs = bytes_per_call * total_calls / dt.as_secs_f64() / 1e9;
            println!(
                "  k={k:>1} (n={n:>5}): {ns_per_call:>7.1} ns/call, {gbs:>6.1} GB/s, {:>10.0} calls/s, {} iters × {} reps",
                total_calls / dt.as_secs_f64(),
                iters,
                reps,
            );
        }
    }
}
