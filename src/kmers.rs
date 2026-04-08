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

/// Default k-mer size, matching the C++ `KMER_SIZE` constant in `dada.h`.
pub const KMER_SIZE: usize = 5;

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
        if nti < 1 || nti > 4 {
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

/// Populate all k-mer fields (`kmer`, `kmer8`, `kord`) on a `Raw`.
pub fn raw_assign_kmers(raw: &mut Raw, k: usize) {
    raw.kmer = Some(assign_kmer(&raw.seq, k));
    raw.kmer8 = Some(assign_kmer8(&raw.seq, k));
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
