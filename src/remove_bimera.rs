//! De novo bimera removal from a sequence table.
//!
//! Mirrors R's `removeBimeraDenovo` with `method = "consensus"` (default),
//! `"pooled"`, or `"per-sample"`.
//!
//! Input/output: [`SequenceTable`] (the JSON produced by `make-sequence-table`).

use rayon::prelude::*;

use crate::chimera::{is_bimera_with_buf, table_bimera2};
use crate::nwalign::AlignBuffers;
use crate::sequence_table::SequenceTable;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

pub enum Method {
    Consensus,
    Pooled,
    PerSample,
}

pub struct BimeraParams {
    pub min_fold_parent_over_abundance: f64,
    pub min_parent_abundance: u32,
    pub allow_one_off: bool,
    pub min_one_off_parent_distance: usize,
    pub max_shift: i32,
    /// Consensus only: fraction of samples in which a sequence must be flagged.
    pub min_sample_fraction: f64,
    /// Consensus only: number of unflagged samples to ignore.
    pub ignore_n_negatives: u32,
    pub match_score: i16,
    pub mismatch: i16,
    pub gap_p: i16,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Remove bimeric sequences from `table` and return the filtered table.
pub fn remove_bimera_denovo(
    table: SequenceTable,
    method: &Method,
    params: &BimeraParams,
    verbose: bool,
) -> SequenceTable {
    let nrow = table.samples.len();
    let ncol = table.sequences.len();

    if ncol == 0 || nrow == 0 {
        return table;
    }

    // Build column-major u32 matrix: mat[i + j*nrow] = counts[i][j]
    let mat: Vec<u32> = {
        let counts = &table.counts;
        (0..ncol)
            .flat_map(|j| (0..nrow).map(move |i| counts[i][j] as u32))
            .collect()
    };

    let seq_bytes: Vec<Vec<u8>> = table
        .sequences
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let seqs_ref: Vec<&[u8]> = seq_bytes.iter().map(|v| v.as_slice()).collect();

    match method {
        Method::Consensus => {
            let col_flags = consensus_flags(&mat, nrow, ncol, &seqs_ref, params);
            if verbose {
                eprintln!(
                    "[remove-bimera-denovo] Identified {} bimera(s) out of {} sequence(s).",
                    col_flags.iter().filter(|&&b| b).count(),
                    ncol
                );
            }
            drop_columns(table, &col_flags)
        }
        Method::Pooled => {
            let col_flags = pooled_flags(&mat, nrow, ncol, &seqs_ref, params);
            if verbose {
                eprintln!(
                    "[remove-bimera-denovo] Identified {} bimera(s) out of {} sequence(s).",
                    col_flags.iter().filter(|&&b| b).count(),
                    ncol
                );
            }
            drop_columns(table, &col_flags)
        }
        Method::PerSample => {
            // cell_bim[i][j] = true if seq j is bimeric in sample i
            let cell_bim = per_sample_cell_flags(&mat, nrow, ncol, &seqs_ref, params);
            let n_flagged: usize = cell_bim
                .iter()
                .flat_map(|r| r.iter())
                .filter(|&&b| b)
                .count();
            if verbose {
                eprintln!(
                    "[remove-bimera-denovo] Zeroed {} cell(s) across {} sample(s).",
                    n_flagged, nrow
                );
            }
            zero_and_drop(table, &cell_bim, nrow, ncol)
        }
    }
}

// ---------------------------------------------------------------------------
// Method implementations
// ---------------------------------------------------------------------------

/// Consensus: table_bimera2 + fraction vote per sequence.
fn consensus_flags(
    mat: &[u32],
    nrow: usize,
    ncol: usize,
    seqs: &[&[u8]],
    p: &BimeraParams,
) -> Vec<bool> {
    let flags = table_bimera2(
        mat,
        nrow,
        ncol,
        seqs,
        p.min_fold_parent_over_abundance,
        p.min_parent_abundance,
        p.allow_one_off,
        p.min_one_off_parent_distance,
        p.match_score,
        p.mismatch,
        p.gap_p,
        p.max_shift,
    );

    flags
        .iter()
        .map(|f| {
            // Mirror R: nflag >= nsam || (nflag > 0 && nflag >= (nsam - ignoreN) * minFrac)
            f.nflag >= f.nsam
                || (f.nflag > 0
                    && f.nflag as f64
                        >= f.nsam.saturating_sub(p.ignore_n_negatives) as f64
                            * p.min_sample_fraction)
        })
        .collect()
}

/// Pooled: sum abundances across all samples, then run is_bimera on each seq.
fn pooled_flags(
    mat: &[u32],
    nrow: usize,
    ncol: usize,
    seqs: &[&[u8]],
    p: &BimeraParams,
) -> Vec<bool> {
    let pooled: Vec<u32> = (0..ncol)
        .map(|j| (0..nrow).map(|i| mat[i + j * nrow]).sum())
        .collect();

    (0..ncol)
        .into_par_iter()
        .map_init(AlignBuffers::new, |buf, j| {
            let abund = pooled[j];
            if abund == 0 {
                return false;
            }
            let parents: Vec<&[u8]> = (0..ncol)
                .filter(|&k| {
                    k != j
                        && pooled[k] as f64 > p.min_fold_parent_over_abundance * abund as f64
                        && pooled[k] >= p.min_parent_abundance
                })
                .map(|k| seqs[k])
                .collect();

            if parents.len() < 2 {
                return false;
            }

            is_bimera_with_buf(
                seqs[j],
                &parents,
                p.allow_one_off,
                p.min_one_off_parent_distance,
                p.match_score,
                p.mismatch,
                p.gap_p,
                p.max_shift,
                buf,
            )
        })
        .collect()
}

/// Per-sample: returns a cell-level bimera flag matrix (nrow × ncol).
fn per_sample_cell_flags(
    mat: &[u32],
    nrow: usize,
    ncol: usize,
    seqs: &[&[u8]],
    p: &BimeraParams,
) -> Vec<Vec<bool>> {
    let mut buf = AlignBuffers::new();
    (0..nrow)
        .map(|i| {
            (0..ncol)
                .map(|j| {
                    let abund = mat[i + j * nrow];
                    if abund == 0 {
                        return false;
                    }
                    let parents: Vec<&[u8]> = (0..ncol)
                        .filter(|&k| {
                            let k_abund = mat[i + k * nrow];
                            k != j
                                && k_abund as f64 > p.min_fold_parent_over_abundance * abund as f64
                                && k_abund >= p.min_parent_abundance
                        })
                        .map(|k| seqs[k])
                        .collect();

                    if parents.len() < 2 {
                        return false;
                    }

                    is_bimera_with_buf(
                        seqs[j],
                        &parents,
                        p.allow_one_off,
                        p.min_one_off_parent_distance,
                        p.match_score,
                        p.mismatch,
                        p.gap_p,
                        p.max_shift,
                        &mut buf,
                    )
                })
                .collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Table filtering helpers
// ---------------------------------------------------------------------------

/// Drop columns flagged as bimeric (consensus / pooled).
fn drop_columns(table: SequenceTable, bimeric: &[bool]) -> SequenceTable {
    let keep = |j: usize| !bimeric[j];
    let ncol = table.sequences.len();

    let sequences = table
        .sequences
        .into_iter()
        .enumerate()
        .filter_map(|(j, s)| keep(j).then_some(s))
        .collect();
    let sequence_ids = table
        .sequence_ids
        .into_iter()
        .enumerate()
        .filter_map(|(j, s)| keep(j).then_some(s))
        .collect();
    let counts = table
        .counts
        .into_iter()
        .map(|row| {
            row.into_iter()
                .enumerate()
                .filter_map(|(j, c)| keep(j).then_some(c))
                .collect()
        })
        .collect();

    SequenceTable {
        samples: table.samples,
        sequences,
        sequence_ids,
        counts,
    }
}

/// Zero flagged cells then drop columns that are entirely zero (per-sample).
fn zero_and_drop(
    table: SequenceTable,
    cell_bim: &[Vec<bool>],
    nrow: usize,
    ncol: usize,
) -> SequenceTable {
    // Build zeroed count matrix.
    let mut counts: Vec<Vec<u64>> = table.counts.clone();
    for i in 0..nrow {
        for j in 0..ncol {
            if cell_bim[i][j] {
                counts[i][j] = 0;
            }
        }
    }

    // Determine which columns are all-zero.
    let keep: Vec<bool> = (0..ncol)
        .map(|j| (0..nrow).any(|i| counts[i][j] > 0))
        .collect();

    let sequences = table
        .sequences
        .into_iter()
        .enumerate()
        .filter_map(|(j, s)| keep[j].then_some(s))
        .collect();
    let sequence_ids = table
        .sequence_ids
        .into_iter()
        .enumerate()
        .filter_map(|(j, s)| keep[j].then_some(s))
        .collect();
    let counts = counts
        .into_iter()
        .map(|row| {
            row.into_iter()
                .enumerate()
                .filter_map(|(j, c)| keep[j].then_some(c))
                .collect()
        })
        .collect();

    SequenceTable {
        samples: table.samples,
        sequences,
        sequence_ids,
        counts,
    }
}
