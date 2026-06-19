//! Higher-order chimera (trimera) diagnostics over a sequence table.
//!
//! `remove-bimera-denovo` answers a binary question per sequence — is it a
//! chimera of two more-abundant parents? — and discards the coverage that led
//! to that answer. Some lower-abundance reads (notably on long amplicons such
//! as full-length 16S or `nodA`, and in low-biomass samples) survive bimera
//! removal yet look chimeric, consistent with being chimeras of *three or more*
//! parents. The bimera coverage signal is exactly what flags these: a read
//! whose best single-junction model nearly spans its full length leaves a small
//! internal gap that a third parent can fill.
//!
//! This module runs [`pooled_diagnostics`] over a [`SequenceTable`] and emits a
//! per-sequence TSV. It does not filter or modify the table.

use std::io::{self, Write};

use crate::chimera::{BimeraAlignParams, BimeraDiagnostic, pooled_diagnostics};
use crate::remove_bimera::BimeraParams;
use crate::sequence_table::SequenceTable;

/// One TSV row of diagnostics for a sequence.
pub struct DiagnosticRow<'a> {
    pub sequence_id: &'a str,
    pub length: usize,
    pub total_abundance: u32,
    pub n_samples: u32,
    pub is_bimera: bool,
    pub max_left: usize,
    pub max_right: usize,
    pub cover: usize,
    pub cover_frac: f64,
    pub gap_len: usize,
    pub gap_start: usize,
    pub gap_end: usize,
    pub best_left_parent: Option<&'a str>,
    pub best_right_parent: Option<&'a str>,
    pub third_parent: Option<&'a str>,
    pub gap_mismatches: usize,
    /// Heuristic screen: survived the strict bimera test but a single junction
    /// covers at least `min_cover_frac` of the read, leaving a real gap.
    pub trimera_suspect: bool,
}

/// Run the pooled bimera coverage diagnostic over `table`.
///
/// `min_cover_frac` sets the `trimera_suspect` threshold on `cover / length`.
/// Returns one [`DiagnosticRow`] per sequence, in table column order.
pub fn run_diagnostics<'a>(
    table: &'a SequenceTable,
    params: &BimeraParams,
    min_cover_frac: f64,
) -> Vec<DiagnosticRow<'a>> {
    let ncol = table.sequences.len();
    let nrow = table.samples.len();
    if ncol == 0 || nrow == 0 {
        return Vec::new();
    }

    let pooled: Vec<u32> = (0..ncol)
        .map(|j| (0..nrow).map(|i| table.counts[i][j] as u32).sum())
        .collect();
    let n_samples: Vec<u32> = (0..ncol)
        .map(|j| (0..nrow).filter(|&i| table.counts[i][j] > 0).count() as u32)
        .collect();

    let seq_bytes: Vec<&[u8]> = table.sequences.iter().map(|s| s.as_bytes()).collect();

    let align_params = BimeraAlignParams {
        allow_one_off: params.allow_one_off,
        min_one_off_par_dist: params.min_one_off_parent_distance,
        match_score: params.match_score,
        mismatch: params.mismatch,
        gap_p: params.gap_p,
        max_shift: params.max_shift,
        backend: params.backend,
        wfa_max_edits: params.wfa_max_edits,
    };

    let diags = pooled_diagnostics(
        &pooled,
        &seq_bytes,
        params.min_fold_parent_over_abundance,
        params.min_parent_abundance,
        &align_params,
    );

    let id = |k: Option<usize>| k.map(|k| table.sequence_ids[k].as_str());

    diags
        .iter()
        .enumerate()
        .map(|(j, d): (usize, &BimeraDiagnostic)| {
            let cover = (d.max_left + d.max_right).min(d.sqlen);
            let cover_frac = if d.sqlen > 0 {
                cover as f64 / d.sqlen as f64
            } else {
                0.0
            };
            let gap_len = d.gap_end.saturating_sub(d.gap_start);
            let trimera_suspect = !d.is_bimera && gap_len > 0 && cover_frac >= min_cover_frac;
            DiagnosticRow {
                sequence_id: table.sequence_ids[j].as_str(),
                length: d.sqlen,
                total_abundance: pooled[j],
                n_samples: n_samples[j],
                is_bimera: d.is_bimera,
                max_left: d.max_left,
                max_right: d.max_right,
                cover,
                cover_frac,
                gap_len,
                gap_start: d.gap_start,
                gap_end: d.gap_end,
                best_left_parent: id(d.best_left_parent),
                best_right_parent: id(d.best_right_parent),
                third_parent: id(d.third_parent),
                gap_mismatches: d.gap_mismatches,
                trimera_suspect,
            }
        })
        .collect()
}

const HEADER: &str = "sequence_id\tlength\ttotal_abundance\tn_samples\tis_bimera\t\
max_left\tmax_right\tcover\tcover_frac\tgap_len\tgap_start\tgap_end\t\
best_left_parent\tbest_right_parent\tthird_parent\tgap_mismatches\ttrimera_suspect";

/// Write diagnostics rows as TSV (header + one row per sequence) to `w`.
pub fn write_tsv<W: Write>(rows: &[DiagnosticRow<'_>], w: &mut W) -> io::Result<()> {
    writeln!(w, "{HEADER}")?;
    let na = "NA";
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.sequence_id,
            r.length,
            r.total_abundance,
            r.n_samples,
            r.is_bimera,
            r.max_left,
            r.max_right,
            r.cover,
            r.cover_frac,
            r.gap_len,
            r.gap_start,
            r.gap_end,
            r.best_left_parent.unwrap_or(na),
            r.best_right_parent.unwrap_or(na),
            r.third_parent.unwrap_or(na),
            r.gap_mismatches,
            r.trimera_suspect,
        )?;
    }
    Ok(())
}
