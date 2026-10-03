//! Helpers shared by more than one command group.

use std::{io, path::Path};

use dada2_rs::sequence_table;

use sequence_table::SequenceTable;

/// Default WFA edit-budget cap (issue #51), in edit operations, used when the
/// experimental `--align-backend wfa2` is selected without an explicit
/// `--wfa-max-edits`. See [`nwalign::AlignParams::wfa_max_edits`].
pub(super) const WFA_MAX_EDITS_DEFAULT: i32 = 50;

/// Verify that every path given for `flag` exists on disk.
///
/// A shell glob that matches nothing is handed to us verbatim (bash and zsh
/// only drop it under `nullglob`/`failglob`), so "no files here" arrives as one
/// entry whose name still contains `*`, `?`, or `[`. Say so explicitly, rather
/// than letting a downstream count comparison blame some other flag.
pub(super) fn check_input_paths(flag: &str, paths: &[std::path::PathBuf]) -> io::Result<()> {
    let missing: Vec<&std::path::PathBuf> = paths.iter().filter(|p| !p.exists()).collect();
    if missing.is_empty() {
        return Ok(());
    }
    let unexpanded = missing
        .iter()
        .any(|p| p.to_string_lossy().contains(['*', '?', '[']));
    let list: Vec<String> = missing
        .iter()
        .take(5)
        .map(|p| p.display().to_string())
        .collect();
    let more = if missing.len() > 5 {
        format!(" (and {} more)", missing.len() - 5)
    } else {
        String::new()
    };
    let hint = if unexpanded {
        "; the pattern looks like a shell glob that matched no files, \
         so it was passed through literally — check the directory path"
    } else {
        ""
    };
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!(
            "{flag}: {} of {} path(s) do not exist: {}{more}{hint}",
            missing.len(),
            paths.len(),
            list.join(", ")
        ),
    ))
}

/// Emit a one-line note when homopolymer gapping is active. Because
/// `homo_gap_p != gap_p` forces the slow scalar aligner — the vectorized SIMD
/// path can't do homopolymer gaps, exactly why R sets `VECTORIZED_ALIGNMENT
/// <- FALSE` (dada.R:229-230) — this is an easy performance gotcha to miss
/// (e.g. a stray `--homo-gap-p -1` on HiFi). Surfaced under --verbose.
pub(super) fn note_homopolymer_gapping(verbose: bool, gap_p: i32, homo_gap_p: i32) {
    if verbose && homo_gap_p != gap_p {
        eprintln!(
            "[align] note: homopolymer gapping active (homo_gap_p={homo_gap_p} != \
             gap_p={gap_p}) — vectorized aligner disabled; alignments use the slower scalar DP."
        );
    }
}

/// Sequence-table column filter mirroring R DADA2's pseudo-pooling prior selection:
///   keep[j] = (n_samples_present[j] >= prevalence) || (total_abundance[j] >= min_abundance)
/// When both thresholds are `None` every column is kept.
pub(super) fn select_sequences(
    table: &SequenceTable,
    prevalence: Option<u32>,
    min_abundance: Option<u64>,
) -> Vec<usize> {
    let nseq = table.sequences.len();
    if prevalence.is_none() && min_abundance.is_none() {
        return (0..nseq).collect();
    }
    (0..nseq)
        .filter(|&j| {
            let by_prev = prevalence.is_some_and(|p| {
                let n_present = table.counts.iter().filter(|row| row[j] > 0).count() as u32;
                n_present >= p
            });
            let by_abund = min_abundance.is_some_and(|m| {
                let total: u64 = table.counts.iter().map(|row| row[j]).sum();
                total >= m
            });
            by_prev || by_abund
        })
        .collect()
}

/// Original file name without the directory path, for provenance. Recorded in
/// dada/derep JSON so downstream commands (e.g. `merge-pairs`) can verify they
/// were handed the file the JSON was actually computed from.
pub(super) fn file_basename(path: &Path) -> String {
    path.file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown")
        .to_string()
}

pub(super) fn fastq_stem(path: &Path) -> String {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown");

    for suffix in &[".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if let Some(stem) = name.strip_suffix(suffix) {
            return stem.to_string();
        }
    }
    // Fallback: use whatever Path::file_stem gives.
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string()
}
