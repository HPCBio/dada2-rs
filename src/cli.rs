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

    /// Denoise a FASTQ file using the DADA2 algorithm
    ///
    /// Dereplicate the input FASTQ in memory, then run DADA2 sample inference
    /// using the supplied error model.  Outputs a JSON object describing the
    /// inferred ASVs.
    ///
    /// By default `err_out` from the error model file is used as the error
    /// matrix.  Pass `--use-err-in` to use `err_in` instead.
    Dada {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

        /// JSON error model file produced by the `learn-errors` subcommand
        #[arg(long)]
        error_model: PathBuf,

        /// Use `err_in` from the error model instead of `err_out`
        #[arg(long)]
        use_err_in: bool,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads (used for both dereplication and DADA2 comparisons)
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long, default_value_t = 1e-40)]
        omega_a: f64,

        /// Significance threshold for reads not corrected to any center (omega_c)
        #[arg(long, default_value_t = 0.0)]
        omega_c: f64,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long, default_value_t = 1e-4)]
        omega_p: f64,

        /// Minimum fold-enrichment above expected for cluster splitting
        #[arg(long, default_value_t = 1.0)]
        min_fold: f64,

        /// Minimum Hamming distance required for cluster splitting
        #[arg(long, default_value_t = 1)]
        min_hamming: u32,

        /// Minimum read abundance required for cluster splitting
        #[arg(long, default_value_t = 1)]
        min_abund: u32,

        /// Use singleton detection
        #[arg(long)]
        detect_singletons: bool,

        /// Include a per-unique-to-ASV index map in the output
        #[arg(long)]
        show_map: bool,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Learn an error model from subsampled derep JSON files
    ///
    /// Reads one or more JSON files produced by the `subsample` subcommand,
    /// iteratively runs the DADA2 algorithm and re-fits the chosen error
    /// model until self-consistency.
    ///
    /// Output is a JSON object with three flat 16 × nq matrices:
    ///   `trans`   — accumulated transition counts,
    ///   `err_in`  — error rates used in the final DADA run,
    ///   `err_out` — error rates estimated from `trans`.
    LearnErrors {
        /// One or more JSON derep files produced by the `subsample` subcommand
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Error model fitting function to use
        ///
        /// Allowed values: loess (default), noqual, binned-qual, pacbio
        #[arg(long, default_value = "loess")]
        errfun: String,

        /// Pseudocount added to each transition total (only used with --errfun noqual)
        #[arg(long, default_value_t = 1.0)]
        pseudocount: f64,

        /// Anchor quality-score bins for piecewise-linear interpolation
        ///
        /// Comma-separated list of quality score values, e.g. "0,10,20,30,40".
        /// Only used with --errfun binned-qual.
        #[arg(long, value_delimiter = ',')]
        binned_quals: Option<Vec<f64>>,

        /// Maximum self-consistency iterations (mirrors R's MAX_CONSIST)
        #[arg(long, default_value_t = 10)]
        max_consist: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long, default_value_t = 1e-40)]
        omega_a: f64,

        /// Significance threshold for omega_c (reads not corrected to any center)
        #[arg(long, default_value_t = 0.0)]
        omega_c: f64,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long, default_value_t = 1e-4)]
        omega_p: f64,

        /// Minimum fold-enrichment above expected for cluster splitting
        #[arg(long, default_value_t = 1.0)]
        min_fold: f64,

        /// Minimum Hamming distance required for cluster splitting
        #[arg(long, default_value_t = 1)]
        min_hamming: u32,

        /// Minimum read abundance required for cluster splitting
        #[arg(long, default_value_t = 1)]
        min_abund: u32,

        /// Use singleton detection (detect singletons as genuine)
        #[arg(long)]
        detect_singletons: bool,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print per-iteration progress to stderr
        #[arg(long)]
        verbose: bool,
    },
}
