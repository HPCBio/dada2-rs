#![allow(clippy::doc_overindented_list_items)]
use std::{io, process::ExitCode};

// The module tree lives in the `dada2_rs` library crate (src/lib.rs); bring the
// modules this binary uses into scope so the existing `foo::Bar` paths resolve.
use dada2_rs::{cli, misc};

use clap::FromArgMatches;
use cli::{Cli, Commands};
use misc::DADA2_RS_VERSION;

mod commands;

/// Print a failure and exit non-zero, instead of letting the std runtime print
/// `Error: {:?}` of the `io::Error` — which dumps the struct
/// (`Custom { kind: Other, error: "..." }`) and buries the message.
///
/// Format is stable and meant to be both readable and parseable:
///
/// ```text
/// dada2-rs: error[<ErrorKind>]: <message>
/// ```
///
/// The bracketed token is `io::ErrorKind`'s debug name (`NotFound`,
/// `InvalidInput`, `InvalidData`, `Other`, …), so a consumer can split on
/// `"]: "` for the message or match the kind without parsing prose. The message
/// is the error's `Display`, which already carries the offending path where
/// `WithPath` was applied. Emitted on stderr; stdout stays clean for data.
fn die(e: &io::Error) -> ExitCode {
    eprintln!("dada2-rs: error[{:?}]: {e}", e.kind());
    ExitCode::FAILURE
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => die(&e),
    }
}

fn run() -> io::Result<()> {
    let cli = Cli::from_arg_matches(&cli::command().get_matches()).unwrap_or_else(|e| e.exit());

    // Before anything else, and irrespective of `--verbose`: a `DADA2RS_*`
    // variable that is set but not recognised means the run is not the one that
    // was asked for, and the entire failure mode is that nobody thinks to check
    // (#145). Two full soil-pool A/B runs were lost to exactly this.
    dada2_rs::gates::warn_unrecognised();
    // Separately: a gate that adds validation work on top of whatever arm is
    // running makes every timing in the log meaningless (#154).
    dada2_rs::gates::warn_timing_invalidating();
    // And a gate that changes results rather than timings (#157).
    dada2_rs::gates::warn_result_changing();

    let command = match cli.command {
        Some(c) => c,
        None => {
            eprintln!("dada2-rs {DADA2_RS_VERSION}");
            eprintln!();
            cli::command().print_help()?;
            return Ok(());
        }
    };

    match command {
        Commands::Summary(args) => commands::qc::run_summary(args)?,
        Commands::SummaryMerge(args) => commands::qc::run_summary_merge(args)?,
        Commands::Derep(args) => commands::derep::run_derep(args)?,
        Commands::Dada(args) => commands::dada::run_dada(args)?,
        Commands::DadaPooled(args) => commands::dada::run_dada_pooled(args)?,
        Commands::DadaPseudo(args) => commands::dada::run_dada_pseudo(args)?,
        Commands::MergePairs(args) => commands::merge::run_merge_pairs(args)?,
        Commands::RemovePrimers(args) => commands::prep::run_remove_primers(args)?,
        Commands::FilterAndTrim(args) => commands::prep::run_filter_and_trim(args)?,
        Commands::MakeSequenceTable(args) => commands::seqtable::run_make_sequence_table(args)?,
        Commands::RemoveBimeraDenovo(args) => commands::chimera::run_remove_bimera_denovo(args)?,
        Commands::ChimeraDiagnostics(args) => commands::chimera::run_chimera_diagnostics(args)?,
        Commands::SeqTableToTsv(args) => commands::seqtable::run_seq_table_to_tsv(args)?,
        Commands::SeqTableToFasta(args) => commands::seqtable::run_seq_table_to_fasta(args)?,
        Commands::TaxToTsv(args) => commands::taxonomy::run_tax_to_tsv(args)?,
        Commands::Sample(args) => commands::derep::run_sample(args)?,
        Commands::ErrorsFromSample(args) => commands::errors::run_errors_from_sample(args)?,
        Commands::AssignTaxonomy(args) => commands::taxonomy::run_assign_taxonomy(args)?,
        Commands::AssignSpecies(args) => commands::taxonomy::run_assign_species(args)?,
        Commands::LearnErrors(args) => commands::errors::run_learn_errors(args)?,
        Commands::KdistCalibrate(args) => commands::errors::run_kdist_calibrate(args)?,
        Commands::ReferenceEval(args) => commands::eval::run_reference_eval(args)?,
    }

    Ok(())
}
