//! Filtering and utility functions.
//!
//! Ports `filter.cpp`, excluding Rcpp export wrappers.
//! Both functions are pure-logic utilities with no dependency on DADA2
//! data structures.

use std::collections::HashSet;

// ---------------------------------------------------------------------------
// match_ref
// ---------------------------------------------------------------------------

/// For each sequence in `seqs`, count how many length-`word_size` words appear
/// in the set of words derived from `reference` (treated circularly: the
/// reference is extended by its first `word_size` characters before indexing).
///
/// When `non_overlapping` is true, each match advances the search position by
/// `word_size` instead of 1, preventing the same reference word from being
/// counted at adjacent offsets.
///
/// Equivalent to C++ `C_matchRef`.
pub fn match_ref(
    seqs: &[&[u8]],
    reference: &[u8],
    word_size: usize,
    non_overlapping: bool,
) -> Vec<u32> {
    let reflen = reference.len();

    // Build the circular-reference word set.
    // C++: ref.append(ref, 0, word_size)  →  extend by first word_size bytes.
    let mut extended = reference.to_vec();
    if word_size <= reflen {
        extended.extend_from_slice(&reference[..word_size]);
    }

    let mut phash: HashSet<&[u8]> = HashSet::with_capacity(reflen);
    for i in 0..reflen {
        if i + word_size <= extended.len() {
            phash.insert(&extended[i..i + word_size]);
        }
    }

    seqs.iter()
        .map(|seq| {
            let slen = seq.len();
            if slen < word_size {
                return 0;
            }
            let mut count = 0u32;
            let mut j = 0usize;
            while j + word_size <= slen {
                if phash.contains(&seq[j..j + word_size]) {
                    count += 1;
                    j += if non_overlapping { word_size } else { 1 };
                } else {
                    j += 1;
                }
            }
            count
        })
        .collect()
}

// ---------------------------------------------------------------------------
// matrix_ee
// ---------------------------------------------------------------------------

/// Compute the expected number of errors for each row of `matrix`.
///
/// Each element is a Phred quality score (non-negative integer). A value of
/// `-1` signals a missing entry (analogous to R's `NA_INTEGER`) and terminates
/// processing of that row early.
///
/// Expected errors = Σ 10^(−Q/10) for each non-missing Q in the row.
///
/// Equivalent to C++ `C_matrixEE`.
#[allow(dead_code)]
pub fn matrix_ee(matrix: &[Vec<i32>]) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .take_while(|&&q| q >= 0)
                .map(|&q| 10_f64.powf(-q as f64 / 10.0))
                .sum()
        })
        .collect()
}
