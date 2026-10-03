//! Evaluation: `reference-eval`.

use std::io;

use dada2_rs::{cli, reference_eval};

pub(crate) fn run_reference_eval(args: cli::ReferenceEvalArgs) -> io::Result<()> {
    let cli::ReferenceEvalArgs {
        asvs,
        reference,
        max_diffs,
        near_diffs,
        kdist_screen,
        kmer_size,
        match_score,
        mismatch,
        gap_p,
        band,
        nraw,
        non_chimeric,
        per_asv,
        per_ref,
        output,
        threads,
        compact,
    } = args;
    reference_eval::run(&reference_eval::Params {
        asvs,
        reference,
        max_diffs,
        near_diffs,
        kdist_screen,
        kmer_size,
        match_score,
        mismatch,
        gap_p,
        band,
        nraw,
        non_chimeric,
        per_asv,
        per_ref,
        out: output,
        threads,
        compact,
    })?;
    Ok(())
}
