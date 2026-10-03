//! Paired-end merging: `merge-pairs`.

use std::io;

use rayon::prelude::*;

use dada2_rs::{cli, merge_pairs, misc};

use misc::Tagged;
use serde::Serialize;

use super::common::check_input_paths;

pub(crate) fn run_merge_pairs(args: cli::MergePairsArgs) -> io::Result<()> {
    let cli::MergePairsArgs {
        fwd_dada,
        rev_dada,
        fwd_fastq,
        rev_fastq,
        min_overlap,
        max_mismatch,
        return_rejects,
        just_concatenate,
        rescue_unmerged,
        concat_nnn_len,
        trim_overhang,
        sample_names,
        check_sample_ids,
        phred_offset,
        threads,
        output,
        compact,
        verbose,
    } = args;
    // ---- Validate that every input path actually exists ----
    // An unmatched shell glob is passed through literally (bash/zsh
    // without `nullglob`/`failglob`), so a wrong directory shows up as
    // a single bogus entry rather than zero. Catching that here keeps
    // the length mismatch below from blaming the wrong flag.
    for (flag, paths) in [
        ("--fwd-dada", &fwd_dada),
        ("--rev-dada", &rev_dada),
        ("--fwd-fastq", &fwd_fastq),
        ("--rev-fastq", &rev_fastq),
    ] {
        check_input_paths(flag, paths)?;
    }

    // ---- Validate that all four lists have the same length ----
    let n = fwd_dada.len();
    for (flag, len) in [
        ("--rev-dada", rev_dada.len()),
        ("--fwd-fastq", fwd_fastq.len()),
        ("--rev-fastq", rev_fastq.len()),
    ] {
        if len != n {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{flag} has {len} entries but --fwd-dada has {n}; \
                     all four file lists must have the same length"
                ),
            ));
        }
    }

    // ---- Resolve sample names ----
    let names: Vec<String> = match sample_names {
        Some(names) => {
            if names.len() != n {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--sample-names has {} entries but {} sample(s) were given",
                        names.len(),
                        n
                    ),
                ));
            }
            names
        }
        None => fwd_dada
            .iter()
            .map(|p| {
                // Strip .json/.json.gz (and any preceding fastq-style extensions) from the stem.
                let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("unknown");
                // Strip trailing .json.gz or .json then apply the FASTQ-stem logic.
                let without_json = name
                    .strip_suffix(".json.gz")
                    .or_else(|| name.strip_suffix(".json"))
                    .unwrap_or(name);
                for suffix in &[".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
                    if let Some(s) = without_json.strip_suffix(suffix) {
                        return s.to_string();
                    }
                }
                without_json.to_string()
            })
            .collect(),
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    let params = merge_pairs::MergeParams {
        min_overlap,
        max_mismatch,
        return_rejects,
        just_concatenate,
        rescue_unmerged,
        concat_nnn_len,
        trim_overhang,
        phred_offset,
        check_sample_ids,
        verbose,
    };

    // Samples are independent, so merge them in parallel across the pool
    // (each sample is mostly serial internally, so across-sample fan-out
    // is what saturates cores). `collect` preserves input order, so the
    // output is identical to a serial run. Nested rayon (merge_sample ->
    // dereplicate also installs on this pool) runs inline + work-steals.
    let results: Vec<merge_pairs::SampleMergeResult> = pool.install(|| {
        (0..n)
            .into_par_iter()
            .map(|i| {
                if verbose {
                    eprintln!("[merge-pairs] sample '{}' ({}/{})", names[i], i + 1, n);
                }
                let result = merge_pairs::merge_sample(
                    &names[i],
                    &fwd_dada[i],
                    &rev_dada[i],
                    &fwd_fastq[i],
                    &rev_fastq[i],
                    &params,
                    &pool,
                )?;
                if verbose {
                    eprintln!(
                        "[merge-pairs] '{}': {}/{} read-pairs accepted → {} merged sequence(s)",
                        names[i], result.accepted_pairs, result.total_pairs, result.num_merged,
                    );
                }
                Ok(result)
            })
            .collect::<io::Result<Vec<_>>>()
    })?;

    #[derive(Serialize)]
    struct MergePairsOutput {
        samples: Vec<merge_pairs::SampleMergeResult>,
    }
    let tagged = Tagged::new("merge-pairs", MergePairsOutput { samples: results });
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
