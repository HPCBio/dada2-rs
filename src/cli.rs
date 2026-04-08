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

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print progress information to stderr
        #[arg(long)]
        verbose: bool,
    },
}
