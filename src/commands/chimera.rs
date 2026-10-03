//! Chimeras: `remove-bimera-denovo`, `chimera-diagnostics`.

use std::io;

use dada2_rs::{chimera_diagnostics, cli, misc, remove_bimera, sequence_table};

use misc::{Tagged, read_tagged_json};
use remove_bimera::{BimeraParams, Method, remove_bimera_denovo};
use sequence_table::SequenceTable;

use super::common::WFA_MAX_EDITS_DEFAULT;
use misc::WithPath;

pub(crate) fn run_remove_bimera_denovo(args: cli::RemoveBimeraDenovoArgs) -> io::Result<()> {
    let cli::RemoveBimeraDenovoArgs {
        input,
        method,
        min_fold_parent_over_abundance,
        min_parent_abundance,
        allow_one_off,
        min_one_off_parent_distance,
        max_shift,
        match_score,
        mismatch,
        gap_p,
        align_backend,
        wfa_max_edits,
        min_sample_fraction,
        ignore_n_negatives,
        threads,
        verbose,
        output,
        compact,
    } = args;
    let table: SequenceTable =
        read_tagged_json(&input, &["make-sequence-table", "remove-bimera-denovo"])
            .with_path(&input)?;

    let method = match method.as_str() {
        "pooled" => Method::Pooled,
        "per-sample" => Method::PerSample,
        _ => Method::Consensus,
    };
    // R's defaults for these two differ by method (issue #211).
    let (default_fold, default_abund) = method.default_parent_abundance();
    let params = BimeraParams {
        min_fold_parent_over_abundance: min_fold_parent_over_abundance.unwrap_or(default_fold),
        min_parent_abundance: min_parent_abundance.unwrap_or(default_abund),
        allow_one_off,
        min_one_off_parent_distance,
        max_shift,
        min_sample_fraction,
        ignore_n_negatives,
        match_score,
        mismatch,
        gap_p,
        backend: align_backend.unwrap_or_default(),
        wfa_max_edits: wfa_max_edits.unwrap_or(WFA_MAX_EDITS_DEFAULT),
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;
    let filtered = pool.install(|| remove_bimera_denovo(table, &method, &params, verbose));

    let tagged = Tagged::new("remove-bimera-denovo", filtered);
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;
    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_chimera_diagnostics(args: cli::ChimeraDiagnosticsArgs) -> io::Result<()> {
    let cli::ChimeraDiagnosticsArgs {
        input,
        min_fold_parent_over_abundance,
        min_parent_abundance,
        max_shift,
        match_score,
        mismatch,
        gap_p,
        align_backend,
        wfa_max_edits,
        trimera_min_parent_dist,
        trimera_min_gap,
        trimera_max_gap_error,
        trimera_min_flank,
        threads,
        output,
    } = args;
    let table: SequenceTable =
        read_tagged_json(&input, &["make-sequence-table", "remove-bimera-denovo"])
            .with_path(&input)?;

    // One-off parameters are unused by the diagnostic (coverage is
    // measured strictly); fill with the standard defaults.
    let params = BimeraParams {
        min_fold_parent_over_abundance,
        min_parent_abundance,
        allow_one_off: false,
        min_one_off_parent_distance: 4,
        max_shift,
        min_sample_fraction: 0.9,
        ignore_n_negatives: 1,
        match_score,
        mismatch,
        gap_p,
        backend: align_backend.unwrap_or_default(),
        wfa_max_edits: wfa_max_edits.unwrap_or(WFA_MAX_EDITS_DEFAULT),
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;
    let crit = chimera_diagnostics::TrimeraCriteria {
        min_parent_dist: trimera_min_parent_dist,
        min_gap_len: trimera_min_gap,
        max_gap_error_frac: trimera_max_gap_error,
        min_flank: trimera_min_flank,
    };
    let rows = pool.install(|| chimera_diagnostics::run_diagnostics(&table, &params, crit));

    let mut out: Box<dyn io::Write> = match output {
        Some(path) => Box::new(io::BufWriter::new(std::fs::File::create(&path)?)),
        None => Box::new(io::BufWriter::new(io::stdout())),
    };
    chimera_diagnostics::write_tsv(&rows, &mut out)?;
    out.flush()?;
    Ok(())
}
