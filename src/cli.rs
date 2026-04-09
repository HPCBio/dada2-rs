use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(about = "DADA2 toolkit")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Compute per-position quality metrics from a FASTQ file
    Summary {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for parallel processing
        #[arg(long, default_value_t = 1)]
        threads: usize,
    },

    /// Dereplicate sequences from a FASTQ file
    ///
    /// Produces the equivalent of the R dada2 `derep` class: a set of unique
    /// sequences with read counts, per-unique mean quality profiles, and a
    /// read-to-unique mapping.
    Derep {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for parallel processing
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Include the per-read mapping (read index → unique index) in the output
        #[arg(long)]
        show_map: bool,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print progress information to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Subsample FASTQ files from a directory up to a maximum base count
    ///
    /// Reads FASTQ files from a directory one by one (optionally in random
    /// order), dereplicated each, and writes one JSON file per input to the
    /// output directory.  Processing stops as soon as the cumulative base
    /// count meets or exceeds `--nbases`.
    ///
    /// Output filenames are derived from the input FASTQ stems
    /// (e.g. `sample1.fastq.gz` → `sample1.json`).
    Subsample {
        /// Directory containing input FASTQ files (.fastq, .fastq.gz, .fq, .fq.gz)
        input_dir: PathBuf,

        /// Directory where per-sample JSON derep files will be written
        output_dir: PathBuf,

        /// Stop after accumulating at least this many total bases across all
        /// processed samples (mirrors R's learnErrors `nbases` parameter)
        #[arg(long, default_value_t = 100_000_000)]
        nbases: u64,

        /// Process files in random order instead of sorted order
        #[arg(long)]
        randomize: bool,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for parallel processing
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Print per-file progress to stderr
        #[arg(long)]
        verbose: bool,
    },
}
