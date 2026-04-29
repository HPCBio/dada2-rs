//! Alignment evaluation and k-mer distance utilities.
//!
//! Ports `evaluate.cpp`, excluding R/Rcpp wrappers.
//! The NW alignment wrappers (`C_nwalign`) are omitted — callers can use
//! the functions in `nwalign` directly.

use crate::kmers::{assign_kmer, assign_kmer_order, kmer_dist, kord_dist};
use crate::misc::nt_encode;

// ---------------------------------------------------------------------------
// Alignment evaluation
// ---------------------------------------------------------------------------

/// Count matches, mismatches, and internal indels in an aligned pair.
///
/// End gaps (leading/trailing gaps in either strand) are excluded; only the
/// internal region of the alignment is scored.
///
/// Returns `(matches, mismatches, indels)`.
/// Equivalent to C++ `C_eval_pair`.
#[allow(dead_code)]
pub fn eval_pair(al0: &[u8], al1: &[u8]) -> Result<(u32, u32, u32), String> {
    if al0.len() != al1.len() {
        return Err(format!(
            "Aligned strings have different lengths ({} vs {}).",
            al0.len(),
            al1.len()
        ));
    }
    let len = al0.len();
    if len == 0 {
        return Ok((0, 0, 0));
    }

    // Find the start of the internal region (skip leading end-gaps).
    let mut start = 0usize;
    let mut s0gap = true;
    let mut s1gap = true;
    while (s0gap || s1gap) && start < len {
        s0gap = s0gap && al0[start] == b'-';
        s1gap = s1gap && al1[start] == b'-';
        if s0gap || s1gap {
            start += 1;
        }
    }

    // Find the end of the internal region (skip trailing end-gaps).
    let mut end = len;
    s0gap = true;
    s1gap = true;
    while (s0gap || s1gap) && end > start {
        end -= 1;
        s0gap = s0gap && al0[end] == b'-';
        s1gap = s1gap && al1[end] == b'-';
    }
    // `end` is now the last internal position (inclusive).

    let mut matches = 0u32;
    let mut mismatches = 0u32;
    let mut indels = 0u32;
    for i in start..=end {
        if al0[i] == b'-' || al1[i] == b'-' {
            indels += 1;
        } else if al0[i] == al1[i] {
            matches += 1;
        } else {
            mismatches += 1;
        }
    }
    Ok((matches, mismatches, indels))
}

/// Compute the consensus of two aligned sequences.
///
/// Positions where the sequences agree take that base.  A gap in one strand
/// is filled by the other.  Mismatches are resolved in favour of sequence
/// `prefer` (1 = `al0` wins, 2 = `al1` wins; anything else → `b'N'`).
///
/// When `trim_overhang` is true, positions where `al0` leads with gaps
/// (i.e. `al1` overhangs on the left) or `al1` trails with gaps (right
/// overhang) are set to gaps before compaction.
///
/// Remaining gaps are removed from the output.
/// Equivalent to C++ `C_pair_consensus`.
#[allow(dead_code)]
pub fn pair_consensus(
    al0: &[u8],
    al1: &[u8],
    prefer: u8,
    trim_overhang: bool,
) -> Result<Vec<u8>, String> {
    if al0.len() != al1.len() {
        return Err(format!(
            "Aligned strings have different lengths ({} vs {}).",
            al0.len(),
            al1.len()
        ));
    }
    let len = al0.len();
    let mut out: Vec<u8> = (0..len)
        .map(|i| {
            if al0[i] == al1[i] || al1[i] == b'-' {
                al0[i]
            } else if al0[i] == b'-' {
                al1[i]
            } else {
                match prefer {
                    1 => al0[i],
                    2 => al1[i],
                    _ => b'N',
                }
            }
        })
        .collect();

    if trim_overhang {
        // Blank positions where al0 has a leading gap (al1 left-overhangs).
        for i in 0..len {
            if al0[i] != b'-' {
                break;
            }
            out[i] = b'-';
        }
        // Blank positions where al1 has a trailing gap (al1 right-overhangs).
        for i in (0..len).rev() {
            if al1[i] != b'-' {
                break;
            }
            out[i] = b'-';
        }
    }

    // Remove gaps.
    Ok(out.into_iter().filter(|&b| b != b'-').collect())
}

// ---------------------------------------------------------------------------
// Sequence checking
// ---------------------------------------------------------------------------

/// Return `true` if every byte of `seq` is one of A/C/G/T (upper-case only).
/// Equivalent to C++ `C_isACGT`.
#[allow(dead_code)]
pub fn is_acgt(seq: &[u8]) -> bool {
    seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

// ---------------------------------------------------------------------------
// K-mer distances on ASCII sequence pairs
// ---------------------------------------------------------------------------

/// Compute k-mer frequency-vector distances for each (s1[i], s2[i]) pair.
///
/// Inputs are ASCII sequences; integer-encoding is handled internally.
/// Equivalent to the R-exported `kmer_dist` in `evaluate.cpp`.
#[allow(dead_code)]
pub fn kmer_dist_pairs(s1: &[&[u8]], s2: &[&[u8]], k: usize) -> Result<Vec<f64>, String> {
    if s1.len() != s2.len() {
        return Err(format!(
            "Mismatched pair counts ({} vs {}).",
            s1.len(),
            s2.len()
        ));
    }
    let n_kmers = 1usize << (2 * k);
    s1.iter()
        .zip(s2.iter())
        .map(|(a, b)| {
            let enc_a: Vec<u8> = a.iter().map(|&c| nt_encode(c)).collect();
            let enc_b: Vec<u8> = b.iter().map(|&c| nt_encode(c)).collect();
            let kv_a = assign_kmer(&enc_a, k);
            let kv_b = assign_kmer(&enc_b, k);
            if kv_a.len() != n_kmers || kv_b.len() != n_kmers {
                return Err(format!("k-mer vector length mismatch (k={k})."));
            }
            Ok(kmer_dist(&kv_a, enc_a.len(), &kv_b, enc_b.len(), k))
        })
        .collect()
}

/// Compute ordered k-mer distances for each (s1[i], s2[i]) pair.
///
/// Equivalent to the R-exported `kord_dist` in `evaluate.cpp` (scalar path,
/// no SSE).
#[allow(dead_code)]
pub fn kord_dist_pairs(s1: &[&[u8]], s2: &[&[u8]], k: usize) -> Result<Vec<f64>, String> {
    if s1.len() != s2.len() {
        return Err(format!(
            "Mismatched pair counts ({} vs {}).",
            s1.len(),
            s2.len()
        ));
    }
    s1.iter()
        .zip(s2.iter())
        .map(|(a, b)| {
            let enc_a: Vec<u8> = a.iter().map(|&c| nt_encode(c)).collect();
            let enc_b: Vec<u8> = b.iter().map(|&c| nt_encode(c)).collect();
            let ko_a = assign_kmer_order(&enc_a, k);
            let ko_b = assign_kmer_order(&enc_b, k);
            Ok(kord_dist(&ko_a, enc_a.len(), &ko_b, enc_b.len(), k))
        })
        .collect()
}

/// Count positional k-mer matches for each (s1[i], s2[i]) pair.
///
/// A positional match is a position `j` (within the shorter k-mer order
/// array) where `kord1[j] == kord2[j]`.
/// Equivalent to the R-exported `kmer_matches` in `evaluate.cpp`.
#[allow(dead_code)]
pub fn kmer_matches_pairs(s1: &[&[u8]], s2: &[&[u8]], k: usize) -> Result<Vec<u32>, String> {
    if s1.len() != s2.len() {
        return Err(format!(
            "Mismatched pair counts ({} vs {}).",
            s1.len(),
            s2.len()
        ));
    }
    s1.iter()
        .zip(s2.iter())
        .map(|(a, b)| {
            let enc_a: Vec<u8> = a.iter().map(|&c| nt_encode(c)).collect();
            let enc_b: Vec<u8> = b.iter().map(|&c| nt_encode(c)).collect();
            let ko_a = assign_kmer_order(&enc_a, k);
            let ko_b = assign_kmer_order(&enc_b, k);
            let klen1 = enc_a.len().saturating_sub(k - 1);
            let klen2 = enc_b.len().saturating_sub(k - 1);
            let min_len = klen1.min(klen2);
            let matches = ko_a[..min_len]
                .iter()
                .zip(ko_b[..min_len].iter())
                .filter(|(x, y)| x == y)
                .count() as u32;
            Ok(matches)
        })
        .collect()
}

/// Count shared k-mers (sum of min(count_a, count_b) over all k-mers) for
/// each (s1[i], s2[i]) pair.
///
/// This is the numerator used in `kmer_dist`.
/// Equivalent to the R-exported `kdist_matches` in `evaluate.cpp`.
#[allow(dead_code)]
pub fn kdist_matches_pairs(s1: &[&[u8]], s2: &[&[u8]], k: usize) -> Result<Vec<u32>, String> {
    if s1.len() != s2.len() {
        return Err(format!(
            "Mismatched pair counts ({} vs {}).",
            s1.len(),
            s2.len()
        ));
    }
    s1.iter()
        .zip(s2.iter())
        .map(|(a, b)| {
            let enc_a: Vec<u8> = a.iter().map(|&c| nt_encode(c)).collect();
            let enc_b: Vec<u8> = b.iter().map(|&c| nt_encode(c)).collect();
            let kv_a = assign_kmer(&enc_a, k);
            let kv_b = assign_kmer(&enc_b, k);
            let dotsum: u32 = kv_a
                .iter()
                .zip(kv_b.iter())
                .map(|(&x, &y)| x.min(y) as u32)
                .sum();
            Ok(dotsum)
        })
        .collect()
}
