//! Failed-to-denoise unique diagnostic (issue #60).
//!
//! DADA2 reports which input uniques fell out of denoising via the per-unique
//! `map`: a unique whose final abundance p-value is below `omega_c` and that was
//! not corrected to any cluster center gets a `null` (`None`) map entry — Ben
//! Callahan's guidance on benjjneb/dada2#1899 is that these NAs are how you trace
//! reads that failed to denoise. This module turns that signal into a direct
//! artifact: a tidy long-format TSV of the dropped uniques and the reads they
//! cost, so users don't have to hand-join the dada `map` to the derep uniques.
//!
//! Format: a header line `sequence<TAB>sample<TAB>reads`, then one row per failed
//! unique per sample it appears in. For per-sample modes (`dada`, `dada-pseudo`)
//! each row's `reads` is the unique's in-sample abundance; for `dada-pooled` the
//! "failed" decision is global (made once on the merged unique table) and a row
//! is emitted per sample the failed merged unique appears in.

use std::io::{self, Write};
use std::path::Path;

/// One failed-to-denoise unique, scoped to a single sample.
pub struct Row {
    pub sequence: String,
    pub sample: String,
    pub reads: u32,
}

/// Write `rows` as the tidy long-format TSV (with header) to `path`, creating
/// parent directories as needed. Rows are sorted by sample, then descending
/// reads, then sequence, so output is deterministic regardless of the order in
/// which concurrently denoised samples contributed them.
pub fn write_tsv(path: &Path, mut rows: Vec<Row>) -> io::Result<usize> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    rows.sort_by(|a, b| {
        a.sample
            .cmp(&b.sample)
            .then(b.reads.cmp(&a.reads))
            .then(a.sequence.cmp(&b.sequence))
    });
    let mut w = io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "sequence\tsample\treads")?;
    for r in &rows {
        writeln!(w, "{}\t{}\t{}", r.sequence, r.sample, r.reads)?;
    }
    w.flush()?;
    Ok(rows.len())
}
