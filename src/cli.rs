use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand};

use crate::misc::DADA2_RS_VERSION;
use crate::nwalign::{AlignBackend, ScreenBackend};

/// Help headings shared across subcommands. Keep this list closed: a new
/// heading means a new category for every subcommand that has such a flag,
/// so prefer reusing one of these over inventing a variant.
const H_INPUT: &str = "Input";
const H_OUTPUT: &str = "Output";
const H_ERRMODEL: &str = "Error model";
const H_DENOISE: &str = "Denoising";
const H_ALIGN: &str = "Alignment";
const H_SCREEN: &str = "Screening";
const H_METRICS: &str = "Metrics";
const H_PSEUDO: &str = "Pseudo-pooling";
const H_FILTER: &str = "Filtering";
const H_TRIM: &str = "Trimming";
const H_PRIMER: &str = "Primers";
const H_MERGE: &str = "Merging";
const H_CHIMERA: &str = "Chimera";
const H_EVAL: &str = "Evaluation";
const H_REGIME: &str = "Pooling regime";
const H_TAX: &str = "Classification";
const H_DIAG: &str = "Diagnostics";
const H_PERF: &str = "Performance";
const H_EXP: &str = "Experimental";

/// Points at the subcommand's full parameter reference on ReadTheDocs. Detailed
/// prose lives there, not in `--help` (issue #168).
macro_rules! docs_link {
    ($page:literal) => {
        concat!(
            "Full parameter reference: https://dada2-rs.readthedocs.io/en/latest/commands/",
            $page,
            "/"
        )
    };
}

#[derive(Parser)]
#[command(
    about = "DADA2 toolkit",
    long_about = "DADA2 toolkit\n\n\
                  For subcommands that take a single JSON input file, pass `-` \
                  to read from stdin (gzip auto-detected from the leading magic \
                  bytes). Output flags such as -o/--output remain explicit; \
                  omit them to write to stdout.",
    version = DADA2_RS_VERSION,
    disable_version_flag = true,
)]
pub struct Cli {
    /// Print the dada2-rs version and exit.
    #[arg(short = 'v', long = "version", action = ArgAction::Version)]
    version: Option<bool>,

    #[command(subcommand)]
    pub command: Option<Commands>,
}

/// Top-level `--help` groups, in display order (issue #83). clap 4 has no
/// subcommand groups, so [`command`] renders this list itself. Every visible
/// subcommand must appear exactly once; a test enforces it.
const COMMAND_GROUPS: &[(&str, &[&str])] = &[
    (
        "Workflow",
        &[
            "remove-primers",
            "filter-and-trim",
            "derep",
            "learn-errors",
            "dada",
            "dada-pooled",
            "dada-pseudo",
            "merge-pairs",
            "make-sequence-table",
            "remove-bimera-denovo",
            "assign-taxonomy",
            "assign-species",
        ],
    ),
    ("Quality control", &["summary", "summary-merge"]),
    (
        "Export",
        &["seq-table-to-tsv", "seq-table-to-fasta", "tax-to-tsv"],
    ),
    (
        "Diagnostics",
        &[
            "sample",
            "errors-from-sample",
            "chimera-diagnostics",
            "kdist-calibrate",
            "reference-eval",
        ],
    ),
];

/// The top-level command with its subcommands listed under [`COMMAND_GROUPS`].
/// Parse and print help through this, not `Cli::command()`.
pub fn command() -> clap::Command {
    use clap::CommandFactory;

    let cmd = Cli::command();
    let styles = cmd.get_styles();
    let (header, literal) = (styles.get_header(), styles.get_literal());
    let width = COMMAND_GROUPS
        .iter()
        .flat_map(|(_, names)| names.iter().map(|n| n.len()))
        .max()
        .unwrap_or(0);

    let mut groups = String::new();
    for (title, names) in COMMAND_GROUPS {
        groups.push_str(&format!("{header}{title}:{header:#}\n"));
        for name in *names {
            let about = cmd
                .find_subcommand(name)
                .and_then(|s| s.get_about())
                .map(|a| a.to_string())
                .unwrap_or_default();
            groups.push_str(&format!("  {literal}{name:width$}{literal:#}  {about}\n"));
        }
        groups.push('\n');
    }

    let template = format!(
        "{{before-help}}{{about-with-newline}}\n{{usage-heading}} {{usage}}\n\n\
         {groups}{header}Options:{header:#}\n{{options}}{{after-help}}"
    );
    cmd.help_template(template)
        .after_help("Run `dada2-rs help <COMMAND>` for a subcommand's flags.")
}

#[derive(Subcommand)]
pub enum Commands {
    /// Compute per-position quality metrics from a FASTQ file
    #[command(after_help = docs_link!("summary"))]
    Summary(SummaryArgs),

    /// Merge per-sample `summary` JSONs into a run-level quality/binning report
    #[command(after_help = docs_link!("summary-merge"))]
    SummaryMerge(SummaryMergeArgs),

    /// Dereplicate sequences from a FASTQ file
    #[command(after_help = docs_link!("derep"))]
    Derep(DerepArgs),

    /// Denoise one or more samples independently (R DADA2 `pool=FALSE`)
    #[command(after_help = docs_link!("dada"))]
    Dada(DadaArgs),

    /// Denoise multiple samples with full pooling (R DADA2 `pool=TRUE`)
    #[command(after_help = docs_link!("dada-pooled"))]
    DadaPooled(DadaPooledArgs),

    /// Denoise multiple samples with pseudo-pooling (R DADA2 `pool="pseudo"`)
    #[command(after_help = docs_link!("dada-pseudo"))]
    DadaPseudo(DadaPseudoArgs),

    /// Merge denoised forward and reverse reads into full-length amplicons
    #[command(after_help = docs_link!("merge-pairs"))]
    MergePairs(MergePairsArgs),

    /// Remove primer sequences from a FASTQ file
    #[command(after_help = docs_link!("remove-primers"))]
    RemovePrimers(RemovePrimersArgs),

    /// Filter and trim a single sample's FASTQ reads
    #[command(after_help = docs_link!("filter-and-trim"))]
    FilterAndTrim(FilterAndTrimArgs),

    /// Build a sample-by-sequence feature table
    #[command(after_help = docs_link!("make-sequence-table"))]
    MakeSequenceTable(MakeSequenceTableArgs),

    /// Remove bimeric sequences from a sequence table
    #[command(after_help = docs_link!("remove-bimera-denovo"))]
    RemoveBimeraDenovo(RemoveBimeraDenovoArgs),

    /// Screen a sequence table for higher-order chimeras (trimeras)
    #[command(after_help = docs_link!("chimera-diagnostics"))]
    ChimeraDiagnostics(ChimeraDiagnosticsArgs),

    /// Convert a sequence table JSON to a tab-delimited count table
    #[command(after_help = docs_link!("seq-table-to-tsv"))]
    SeqTableToTsv(SeqTableToTsvArgs),

    /// Assign taxonomy to sequences using a Naive Bayes k-mer classifier
    #[command(after_help = docs_link!("assign-taxonomy"))]
    AssignTaxonomy(AssignTaxonomyArgs),

    /// Fill in the Species column of an assign-taxonomy JSON by exact match
    #[command(after_help = docs_link!("assign-species"))]
    AssignSpecies(AssignSpeciesArgs),

    /// Convert an assign-taxonomy or assign-species JSON to a TSV table
    #[command(after_help = docs_link!("tax-to-tsv"))]
    TaxToTsv(TaxToTsvArgs),

    /// Convert a make-sequence-table JSON file to FASTA
    #[command(after_help = docs_link!("seq-table-to-fasta"))]
    SeqTableToFasta(SeqTableToFastaArgs),

    /// Dereplicate and subsample FASTQ files, writing one JSON file per sample
    #[command(after_help = docs_link!("sample"))]
    Sample(SampleArgs),

    /// Learn an error model from pre-computed sample JSON files
    #[command(after_help = docs_link!("errors-from-sample"))]
    ErrorsFromSample(ErrorsFromSampleArgs),

    /// Learn an error model from FASTQ or derep/sample JSON files
    #[command(
        after_help = concat!(
            "CAVEAT: --nbases accumulates whole samples, so one deep sample can fill\n\
             the budget and the model may reflect only a few samples' diversity.\n\n",
            docs_link!("learn-errors"),
        ),
    )]
    LearnErrors(LearnErrorsArgs),

    /// Calibrate the k-mer screen: emit kdist vs true alignment divergence
    #[command(after_help = docs_link!("kdist-calibrate"))]
    KdistCalibrate(KdistCalibrateArgs),

    /// Evaluate a query ASV set against a reference (truth) set by NW alignment
    #[command(after_help = docs_link!("reference-eval"))]
    ReferenceEval(ReferenceEvalArgs),
}

/// Denoising, alignment and screening flags shared by `dada`, `dada-pooled` and
/// `dada-pseudo`. Every flag is optional: omitted ones resolve from the error
/// model or the built-in default (`resolve_dada_params`).
#[derive(Args)]
pub struct DadaDenoiseArgs {
    /// Significance threshold for abundance-based cluster splitting (R OMEGA_A)
    #[arg(long, help_heading = H_DENOISE)]
    pub omega_a: Option<f64>,

    /// Significance threshold for reads not corrected to any center (R OMEGA_C)
    #[arg(long, help_heading = H_DENOISE)]
    pub omega_c: Option<f64>,

    /// Significance threshold for prior-sequence splitting (R OMEGA_P)
    #[arg(long, help_heading = H_DENOISE)]
    pub omega_p: Option<f64>,

    /// Minimum fold-enrichment above expected for splitting (R MIN_FOLD)
    #[arg(long, help_heading = H_DENOISE)]
    pub min_fold: Option<f64>,

    /// Minimum Hamming distance required for splitting (R MIN_HAMMING)
    #[arg(long, help_heading = H_DENOISE)]
    pub min_hamming: Option<u32>,

    /// Minimum read abundance required for splitting (R MIN_ABUNDANCE)
    #[arg(long, help_heading = H_DENOISE)]
    pub min_abund: Option<u32>,

    /// Detect singletons as genuine (R DETECT_SINGLETONS) [omit to inherit]
    #[arg(long, help_heading = H_DENOISE)]
    pub detect_singletons: Option<bool>,

    /// Maximum number of clusters to infer, 0 = unlimited (R MAX_CLUST)
    #[arg(long, help_heading = H_DENOISE)]
    pub max_clust: Option<usize>,

    /// Use greedy clustering (R GREEDY) [omit to inherit]
    #[arg(long, help_heading = H_DENOISE)]
    pub greedy: Option<bool>,

    /// Use quality scores in the error model (R USE_QUALS) [omit to inherit]
    #[arg(long, help_heading = H_DENOISE)]
    pub use_quals: Option<bool>,

    /// Band radius (R BAND_SIZE): 16 Illumina, 32 PacBio HiFi, -1 unbanded
    #[arg(long, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub band: Option<i32>,

    /// Gap penalty for the Needleman-Wunsch alignment (R GAP_PENALTY)
    #[arg(long, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub gap_p: Option<i32>,

    /// Homopolymer-run gap penalty (R HOMOPOLYMER_GAP_PENALTY) [default: --gap-p]
    #[arg(long, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub homo_gap_p: Option<i32>,

    /// Match score for the Needleman-Wunsch alignment (R MATCH)
    #[arg(long = "match", allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub match_score: Option<i32>,

    /// Mismatch score for the Needleman-Wunsch alignment (R MISMATCH)
    #[arg(long, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub mismatch: Option<i32>,

    /// Pairwise aligner; `wfa2` requires a build with `--features wfa`
    #[arg(long, value_enum, help_heading = H_ALIGN)]
    pub align_backend: Option<AlignBackend>,

    /// K-mer distance cutoff for the pre-alignment screen (R KDIST_CUTOFF)
    #[arg(long, help_heading = H_SCREEN)]
    pub kdist_cutoff: Option<f64>,

    /// K-mer size for the screen, 3-8 (R KMER_SIZE); use 6-7 on PacBio HiFi
    #[arg(long, help_heading = H_SCREEN)]
    pub kmer_size: Option<usize>,

    /// Disable the k-mer screen and align every pair (much slower)
    #[arg(long, help_heading = H_SCREEN)]
    pub no_kmer_screen: Option<bool>,
}

/// Denoising, alignment and screening flags shared by `learn-errors` and
/// `errors-from-sample`, with R `learnErrors` defaults.
#[derive(Args)]
pub struct LearnDenoiseArgs {
    /// Significance threshold for abundance-based cluster splitting (R OMEGA_A)
    #[arg(long, default_value = "1e-40", help_heading = H_DENOISE)]
    pub omega_a: f64,

    /// Threshold for reads not corrected to any center; R learnErrors uses 0
    #[arg(long, default_value = "0", help_heading = H_DENOISE)]
    pub omega_c: f64,

    /// Significance threshold for prior-sequence splitting (R OMEGA_P)
    #[arg(long, default_value = "1e-4", help_heading = H_DENOISE)]
    pub omega_p: f64,

    /// Minimum fold-enrichment above expected for splitting (R MIN_FOLD)
    #[arg(long, default_value_t = 1.0, help_heading = H_DENOISE)]
    pub min_fold: f64,

    /// Minimum Hamming distance required for splitting (R MIN_HAMMING)
    #[arg(long, default_value_t = 1, help_heading = H_DENOISE)]
    pub min_hamming: u32,

    /// Minimum read abundance required for splitting (R MIN_ABUNDANCE)
    #[arg(long, default_value_t = 1, help_heading = H_DENOISE)]
    pub min_abund: u32,

    /// Detect singletons as genuine (R DETECT_SINGLETONS)
    #[arg(long, help_heading = H_DENOISE)]
    pub detect_singletons: bool,

    /// Maximum number of clusters to infer, 0 = unlimited (R MAX_CLUST)
    #[arg(long, default_value_t = 0, help_heading = H_DENOISE)]
    pub max_clust: usize,

    /// Use greedy clustering (R GREEDY) [omit for the default, true]
    #[arg(long, help_heading = H_DENOISE)]
    pub greedy: Option<bool>,

    /// Use quality scores in the error model (R USE_QUALS) [omit for true]
    #[arg(long, help_heading = H_DENOISE)]
    pub use_quals: Option<bool>,

    /// Band radius (R BAND_SIZE): 16 Illumina, 32 PacBio HiFi, -1 unbanded
    #[arg(long, default_value_t = 16, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub band: i32,

    /// Gap penalty for the Needleman-Wunsch alignment (R GAP_PENALTY) [default: -8]
    #[arg(long, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub gap_p: Option<i32>,

    /// Homopolymer-run gap penalty (R HOMOPOLYMER_GAP_PENALTY) [default: --gap-p]
    #[arg(long, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub homo_gap_p: Option<i32>,

    /// Match score for the Needleman-Wunsch alignment (R MATCH)
    #[arg(long = "match", default_value_t = 5, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub match_score: i32,

    /// Mismatch score for the Needleman-Wunsch alignment (R MISMATCH)
    #[arg(long, default_value_t = -4, allow_hyphen_values = true, help_heading = H_ALIGN)]
    pub mismatch: i32,

    /// Pairwise aligner; `wfa2` requires a build with `--features wfa`
    #[arg(long, value_enum, help_heading = H_ALIGN)]
    pub align_backend: Option<AlignBackend>,

    /// K-mer distance cutoff for the pre-alignment screen (R KDIST_CUTOFF)
    #[arg(long, default_value_t = 0.42, help_heading = H_SCREEN)]
    pub kdist_cutoff: f64,

    /// K-mer size for the screen, 3-8 (R KMER_SIZE); use 6-7 on PacBio HiFi
    #[arg(long, default_value_t = 5, help_heading = H_SCREEN)]
    pub kmer_size: usize,

    /// Disable the k-mer screen and align every pair (much slower)
    #[arg(long, help_heading = H_SCREEN)]
    pub no_kmer_screen: bool,
}

/// Experimental screen and aligner knobs shared by every denoising subcommand.
#[derive(Args)]
pub struct ExperimentalArgs {
    /// EXPERIMENTAL: pre-alignment screen; `minimizer` needs a tuned cutoff
    #[arg(long, value_enum, help_heading = H_EXP)]
    pub screen_backend: Option<ScreenBackend>,

    /// EXPERIMENTAL: k-mer size for the minimizer sketch, 5-31 [default: 8]
    #[arg(long, help_heading = H_EXP)]
    pub minimizer_k: Option<usize>,

    /// EXPERIMENTAL: minimizer winnowing window in k-mers, 1-64 [default: 5]
    #[arg(long, help_heading = H_EXP)]
    pub minimizer_w: Option<usize>,

    /// EXPERIMENTAL: run both screens and report disagreements (much slower)
    #[arg(long, default_value_t = false, help_heading = H_EXP)]
    pub screen_audit: bool,

    /// EXPERIMENTAL: WFA edit-budget cap, 0 = unbounded [default: 50]
    #[arg(long, help_heading = H_EXP)]
    pub wfa_max_edits: Option<i32>,
}

/// Error-function and self-consistency flags shared by `learn-errors` and
/// `errors-from-sample`.
#[derive(Args)]
pub struct ErrModelFitArgs {
    /// Error-model fitting function
    #[arg(long, default_value = "loess", help_heading = H_ERRMODEL,
          value_parser = ["loess", "noqual", "binned-qual", "pacbio", "external"])]
    pub errfun: String,

    /// Pseudocount added to each transition total (--errfun noqual only)
    #[arg(long, default_value_t = 1.0, help_heading = H_ERRMODEL)]
    pub pseudocount: f64,

    /// Anchor quality bins, e.g. "0,10,20,30,40" (--errfun binned-qual only)
    #[arg(long, value_delimiter = ',', help_heading = H_ERRMODEL)]
    pub binned_quals: Option<Vec<f64>>,

    /// Command to invoke for --errfun external; input/output paths are appended
    #[arg(long, help_heading = H_ERRMODEL)]
    pub errfun_cmd: Option<String>,

    /// DEPRECATED: use --loess-surface
    #[arg(long, help_heading = H_ERRMODEL, hide = true,
          value_parser = ["default", "r-dada2"])]
    pub loess_preset: Option<String>,

    /// LOESS fitting surface [default: interpolate, as R's loess() uses]
    #[arg(long, value_parser = ["direct", "interpolate"], help_heading = H_ERRMODEL)]
    pub loess_surface: Option<String>,

    /// Max fraction of observations per kd-tree cell (interpolate surface only)
    #[arg(long, help_heading = H_ERRMODEL)]
    pub loess_cell: Option<f64>,

    /// Upper clamp on fitted off-diagonal error rates; 1.0 disables
    #[arg(long, help_heading = H_ERRMODEL)]
    pub loess_max_rate: Option<f64>,

    /// Lower clamp on fitted off-diagonal error rates; 0.0 disables
    #[arg(long, help_heading = H_ERRMODEL)]
    pub loess_min_rate: Option<f64>,

    /// Maximum self-consistency iterations (R MAX_CONSIST)
    #[arg(long, default_value_t = 10, help_heading = H_ERRMODEL)]
    pub max_consist: usize,
}

/// Per-iteration diagnostics shared by `learn-errors` and `errors-from-sample`.
#[derive(Args)]
pub struct LearnDiagArgs {
    /// Directory for per-iteration cluster diagnostics (iter_NNN.json)
    #[arg(long, help_heading = H_DIAG)]
    pub diag_dir: Option<PathBuf>,

    /// Directory for full per-iteration cluster traces
    #[arg(long, help_heading = H_DIAG)]
    pub cluster_trace_dir: Option<PathBuf>,

    /// Omit the per-cluster `members` array from trace files (~10x smaller)
    #[arg(long, help_heading = H_DIAG)]
    pub trace_no_members: bool,

    /// Only trace members with abundance >= this value
    #[arg(long, default_value_t = 1, help_heading = H_DIAG)]
    pub trace_min_abund: u32,
}

/// Arguments for [`Commands::Summary`].
#[derive(Args)]
pub struct SummaryArgs {
    /// Input FASTQ file (uncompressed or gzipped)
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Sample identifier for the output JSON [default: input filename stem]
    #[arg(long, help_heading = H_INPUT)]
    pub sample_name: Option<String>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Also compute the per-read sequence-complexity histogram
    #[arg(long, help_heading = H_METRICS)]
    pub complexity: bool,

    /// K-mer size for --complexity
    #[arg(long, default_value_t = 2, help_heading = H_METRICS)]
    pub complexity_kmer_size: u8,

    /// Histogram bins for --complexity, spanning `[0, 4^kmer_size]`
    #[arg(long, default_value_t = 100, help_heading = H_METRICS)]
    pub complexity_bins: usize,

    /// Also compute per-position cumulative expected-error (EE) metrics
    #[arg(long, help_heading = H_METRICS)]
    pub expected_error: bool,

    /// Log-spaced histogram bins backing the EE quantiles
    #[arg(long, default_value_t = 200, help_heading = H_METRICS)]
    pub ee_bins: usize,

    /// Max distinct quality values before the data is judged "binned"
    #[arg(long, default_value_t = 8, help_heading = H_METRICS)]
    pub binned_threshold: usize,

    /// Number of threads for parallel processing
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Also print a human-readable metrics report to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub report: bool,
}

/// Arguments for [`Commands::SummaryMerge`].
#[derive(Args)]
pub struct SummaryMergeArgs {
    /// Per-sample `summary` JSON files (gzip ok; `-` reads stdin)
    #[arg(help_heading = H_INPUT)]
    pub inputs: Vec<PathBuf>,

    /// Max distinct quality values before the run is judged "binned"
    #[arg(long, default_value_t = 8, help_heading = H_METRICS)]
    pub binned_threshold: usize,

    /// Declared bin levels to validate the run against, e.g. `2,12,24,40`
    #[arg(long, value_delimiter = ',', help_heading = H_METRICS)]
    pub expected_bins: Vec<u8>,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Also print a human-readable run-level report to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub report: bool,
}

/// Arguments for [`Commands::Derep`].
#[derive(Args)]
pub struct DerepArgs {
    /// Input FASTQ file (uncompressed or gzipped)
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Sample identifier for the output JSON [default: input filename stem]
    #[arg(long, help_heading = H_INPUT)]
    pub sample_name: Option<String>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Number of threads for parallel processing
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON here instead of stdout; a `.gz` path is gzip-compressed
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Include the per-read mapping (read index -> unique index)
    #[arg(long, help_heading = H_OUTPUT)]
    pub show_map: bool,

    /// Pretty-print the JSON; the default is compact (~34% smaller)
    #[arg(long, help_heading = H_OUTPUT)]
    pub pretty: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::Dada`].
#[derive(Args)]
pub struct DadaArgs {
    /// Input FASTQ or derep/sample JSON files; >1 input requires --output-dir
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    /// Sample identifier for the output JSON [default: input filename stem]
    #[arg(long, help_heading = H_INPUT)]
    pub sample_name: Option<String>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// FASTA of prior sequences; exact matches skip the abundance p-value filter
    #[arg(long, help_heading = H_INPUT)]
    pub prior: Option<PathBuf>,

    /// JSON error model file produced by the `learn-errors` subcommand
    #[arg(long, help_heading = H_ERRMODEL)]
    pub error_model: PathBuf,

    /// Use `err_in` from the error model instead of `err_out`
    #[arg(long, help_heading = H_ERRMODEL)]
    pub use_err_in: bool,

    /// Inherit unspecified algorithm parameters from the model's `params` block
    #[arg(long, help_heading = H_ERRMODEL)]
    pub inherit_err_params: bool,

    #[command(flatten)]
    pub denoise: DadaDenoiseArgs,

    /// Number of threads for dereplication and DADA2 comparisons
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Samples to denoise concurrently, multi-input only [default: threads/4]
    #[arg(long, help_heading = H_PERF)]
    pub sample_jobs: Option<usize>,

    /// Write JSON to this file instead of stdout (single input only)
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Directory for per-sample `{sample}.json` (required for >1 input)
    #[arg(long, help_heading = H_OUTPUT)]
    pub output_dir: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Gzip the per-sample JSON files (`{sample}.json.gz`)
    #[arg(long, help_heading = H_OUTPUT)]
    pub gzip: bool,

    /// Write machine-readable run metrics (JSON) to this file
    #[arg(long, help_heading = H_DIAG)]
    pub metrics_json: Option<PathBuf>,

    /// Add per-comparison timings to --metrics-json; slows the run measurably
    #[arg(long, requires = "metrics_json", help_heading = H_DIAG)]
    pub metrics_attribution: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    /// Emit R-parity per-cluster diagnostics in the output JSON
    #[arg(long, help_heading = H_DIAG)]
    pub aux_outputs: bool,

    /// Write a full cluster trace (clusters.json) to this file
    #[arg(long, help_heading = H_DIAG)]
    pub cluster_trace: Option<PathBuf>,

    /// Omit the per-cluster `members` array from the trace
    #[arg(long, help_heading = H_DIAG)]
    pub trace_no_members: bool,

    /// Only trace members with abundance >= this value
    #[arg(long, default_value_t = 1, help_heading = H_DIAG)]
    pub trace_min_abund: u32,

    /// Write a TSV of uniques that failed to denoise to this file
    #[arg(long, help_heading = H_DIAG)]
    pub failed_uniques: Option<PathBuf>,

    #[command(flatten)]
    pub experimental: ExperimentalArgs,
}

/// Arguments for [`Commands::DadaPooled`].
#[derive(Args)]
pub struct DadaPooledArgs {
    /// Input FASTQ or derep/sample JSON files, one per sample
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    /// Comma-separated sample names, one per input [default: input filename stems]
    #[arg(long, value_delimiter = ',', help_heading = H_INPUT)]
    pub sample_names: Option<Vec<String>>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// FASTA of prior sequences; exact matches skip the abundance p-value filter
    #[arg(long, help_heading = H_INPUT)]
    pub prior: Option<PathBuf>,

    /// JSON error model file produced by the `learn-errors` subcommand
    #[arg(long, help_heading = H_ERRMODEL)]
    pub error_model: PathBuf,

    /// Use `err_in` from the error model instead of `err_out`
    #[arg(long, help_heading = H_ERRMODEL)]
    pub use_err_in: bool,

    /// Inherit unspecified algorithm parameters from the model's `params` block
    #[arg(long, help_heading = H_ERRMODEL)]
    pub inherit_err_params: bool,

    #[command(flatten)]
    pub denoise: DadaDenoiseArgs,

    /// Number of threads for dereplication and DADA2 comparisons
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Output directory for per-sample `{sample}.json` (created if absent)
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output_dir: PathBuf,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Gzip the per-sample JSON files (`{sample}.json.gz`)
    #[arg(long, help_heading = H_OUTPUT)]
    pub gzip: bool,

    /// Write machine-readable run metrics (JSON) to this file
    #[arg(long, help_heading = H_DIAG)]
    pub metrics_json: Option<PathBuf>,

    /// Add per-comparison timings to --metrics-json; slows the run measurably
    #[arg(long, requires = "metrics_json", help_heading = H_DIAG)]
    pub metrics_attribution: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    /// Write a TSV of uniques that failed to denoise to this file
    #[arg(long, help_heading = H_DIAG)]
    pub failed_uniques: Option<PathBuf>,

    /// Write a self-contained pooled record here, for `kdist-calibrate`
    #[arg(long, help_heading = H_DIAG)]
    pub pooled_record: Option<PathBuf>,

    /// Write a full cluster trace (clusters.json) to this file
    #[arg(long, help_heading = H_DIAG)]
    pub cluster_trace: Option<PathBuf>,

    /// Omit the per-cluster `members` array from the trace
    #[arg(long, help_heading = H_DIAG)]
    pub trace_no_members: bool,

    /// Only trace members with abundance >= this value
    #[arg(long, default_value_t = 1, help_heading = H_DIAG)]
    pub trace_min_abund: u32,

    #[command(flatten)]
    pub experimental: ExperimentalArgs,
}

/// Arguments for [`Commands::DadaPseudo`].
#[derive(Args)]
pub struct DadaPseudoArgs {
    /// Input FASTQ or derep/sample JSON files, one per sample
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    /// Comma-separated sample names, one per input [default: input filename stems]
    #[arg(long, value_delimiter = ',', help_heading = H_INPUT)]
    pub sample_names: Option<Vec<String>>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// JSON error model file produced by the `learn-errors` subcommand
    #[arg(long, help_heading = H_ERRMODEL)]
    pub error_model: PathBuf,

    /// Use `err_in` from the error model instead of `err_out`
    #[arg(long, help_heading = H_ERRMODEL)]
    pub use_err_in: bool,

    /// Inherit unspecified algorithm parameters from the model's `params` block
    #[arg(long, help_heading = H_ERRMODEL)]
    pub inherit_err_params: bool,

    /// Re-fit the error model from round 1 and use it for round 2
    #[arg(long, help_heading = H_ERRMODEL)]
    pub reestimate_err_between_rounds: bool,

    /// Samples an ASV must appear in to become a prior (R PSEUDO_PREVALENCE)
    #[arg(long, default_value_t = 2, help_heading = H_PSEUDO)]
    pub pseudo_prevalence: u32,

    /// Total abundance for an ASV to become a prior (R PSEUDO_ABUNDANCE)
    #[arg(long, help_heading = H_PSEUDO)]
    pub pseudo_min_abundance: Option<u64>,

    /// Write the selected round-2 priors to this FASTA
    #[arg(long, help_heading = H_PSEUDO)]
    pub priors_out: Option<PathBuf>,

    #[command(flatten)]
    pub denoise: DadaDenoiseArgs,

    /// Number of threads for dereplication and DADA2 comparisons
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Samples to denoise concurrently [default: threads/4]
    #[arg(long, help_heading = H_PERF)]
    pub sample_jobs: Option<usize>,

    /// Keep every sample's uniques in memory across both rounds
    #[arg(long, help_heading = H_PERF)]
    pub cache_samples: bool,

    /// Output directory for per-sample `{sample}.json` (created if absent)
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output_dir: PathBuf,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Gzip the per-sample JSON files (`{sample}.json.gz`)
    #[arg(long, help_heading = H_OUTPUT)]
    pub gzip: bool,

    /// Write machine-readable run metrics (JSON) to this file
    #[arg(long, help_heading = H_DIAG)]
    pub metrics_json: Option<PathBuf>,

    /// Add per-comparison timings to --metrics-json; slows the run measurably
    #[arg(long, requires = "metrics_json", help_heading = H_DIAG)]
    pub metrics_attribution: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    /// Write a TSV of uniques that failed to denoise to this file
    #[arg(long, help_heading = H_DIAG)]
    pub failed_uniques: Option<PathBuf>,

    #[command(flatten)]
    pub experimental: ExperimentalArgs,
}

/// Arguments for [`Commands::MergePairs`].
#[derive(Args)]
pub struct MergePairsArgs {
    /// Forward dada JSON files
    #[arg(long, required = true, num_args = 1.., help_heading = H_INPUT)]
    pub fwd_dada: Vec<PathBuf>,

    /// Reverse dada JSON files
    #[arg(long, required = true, num_args = 1.., help_heading = H_INPUT)]
    pub rev_dada: Vec<PathBuf>,

    /// Forward FASTQ files, re-dereplicated to recover the read->unique map
    #[arg(long, required = true, num_args = 1.., help_heading = H_INPUT)]
    pub fwd_fastq: Vec<PathBuf>,

    /// Reverse FASTQ files, re-dereplicated to recover the read->unique map
    #[arg(long, required = true, num_args = 1.., help_heading = H_INPUT)]
    pub rev_fastq: Vec<PathBuf>,

    /// Comma-separated sample names, one per input set [default: --fwd-dada stems]
    #[arg(long, value_delimiter = ',', help_heading = H_INPUT)]
    pub sample_names: Option<Vec<String>>,

    /// Phred offset for FASTQ re-dereplication
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Minimum overlap between the forward and RC(reverse) ASVs
    #[arg(long, default_value_t = 12, help_heading = H_MERGE)]
    pub min_overlap: u32,

    /// Maximum mismatches allowed in the overlap region
    #[arg(long, default_value_t = 0, help_heading = H_MERGE)]
    pub max_mismatch: u32,

    /// Concatenate with an N spacer instead of merging
    #[arg(long, help_heading = H_MERGE)]
    pub just_concatenate: bool,

    /// Concatenate pairs that fail to merge instead of dropping them
    #[arg(long, help_heading = H_MERGE)]
    pub rescue_unmerged: bool,

    /// Number of N characters in the concatenation spacer
    #[arg(long, default_value_t = 10, help_heading = H_MERGE)]
    pub concat_nnn_len: usize,

    /// Trim read overhangs past the overlap
    #[arg(long, help_heading = H_MERGE)]
    pub trim_overhang: bool,

    /// Number of threads, used within each sample for dereplication
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Include rejected merges (with `accept: false`) in the output
    #[arg(long, help_heading = H_OUTPUT)]
    pub return_rejects: bool,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Verify that the dada JSONs and FASTQ names agree on the sample
    #[arg(long, help_heading = H_DIAG)]
    pub check_sample_ids: bool,

    /// Print per-sample progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::RemovePrimers`].
#[derive(Args)]
pub struct RemovePrimersArgs {
    /// Input FASTQ file (uncompressed or gzipped)
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Sample identifier for the output JSON [default: input filename stem]
    #[arg(long, help_heading = H_INPUT)]
    pub sample_name: Option<String>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Forward primer, 5'->3' (IUPAC ambiguity codes accepted)
    #[arg(long, help_heading = H_PRIMER)]
    pub primer_fwd: String,

    /// Reverse primer, 5'->3'; omit to skip reverse primer detection
    #[arg(long, help_heading = H_PRIMER)]
    pub primer_rev: Option<String>,

    /// Reverse-complement --primer-rev before matching
    #[arg(long, default_value_t = true, help_heading = H_PRIMER)]
    pub rc_primer_rev: bool,

    /// Maximum mismatches allowed when matching each primer
    #[arg(long, default_value_t = 2, help_heading = H_PRIMER)]
    pub max_mismatch: usize,

    /// Also allow indels when matching primers (edit distance; slower)
    #[arg(long, help_heading = H_PRIMER)]
    pub allow_indels: bool,

    /// Trim the forward primer from the 5' end of each read
    #[arg(long, default_value_t = true, help_heading = H_PRIMER)]
    pub trim_fwd: bool,

    /// Trim the reverse primer from the 3' end of each read
    #[arg(long, default_value_t = true, help_heading = H_PRIMER)]
    pub trim_rev: bool,

    /// Flip reads that match primers only in the reverse complement
    #[arg(long, default_value_t = true, help_heading = H_PRIMER)]
    pub orient: bool,

    /// Truncate reads at the first Phred score <= this value
    #[arg(long, help_heading = H_TRIM)]
    pub trunc_q: Option<u8>,

    /// Truncate reads to this many bases; discard if shorter
    #[arg(long, help_heading = H_TRIM)]
    pub trunc_len: Option<usize>,

    /// Remove this many bases from the 5' end of the primer-trimmed read
    #[arg(long, help_heading = H_TRIM)]
    pub trim_left: Option<usize>,

    /// Remove this many bases from the 3' end of the primer-trimmed read
    #[arg(long, help_heading = H_TRIM)]
    pub trim_right: Option<usize>,

    /// Discard reads longer than this before quality trimming
    #[arg(long, help_heading = H_FILTER)]
    pub max_len: Option<usize>,

    /// Discard reads shorter than this after all trimming
    #[arg(long, help_heading = H_FILTER)]
    pub min_len: Option<usize>,

    /// Discard reads with more than this many N bases
    #[arg(long, help_heading = H_FILTER)]
    pub max_n: Option<usize>,

    /// Discard reads with any Phred score below this value
    #[arg(long, help_heading = H_FILTER)]
    pub min_q: Option<u8>,

    /// Discard reads with expected errors above this threshold
    #[arg(long, help_heading = H_FILTER)]
    pub max_ee: Option<f64>,

    /// FASTA of the phiX genome; matching reads are removed
    #[arg(long, help_heading = H_FILTER)]
    pub phix_genome: Option<PathBuf>,

    /// Discard reads with 2-mer Shannon richness below this value
    #[arg(long, help_heading = H_FILTER)]
    pub rm_lowcomplex: Option<f64>,

    /// Threads for primer matching and bgzf output; >1 enables bgzf
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Output FASTQ file
    #[arg(long, short = 'f', help_heading = H_OUTPUT)]
    pub fout: PathBuf,

    /// Gzip-compress the output FASTQ file
    #[arg(long, default_value_t = true, help_heading = H_OUTPUT)]
    pub compress: bool,

    /// Write JSON stats here instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::FilterAndTrim`].
#[derive(Args)]
pub struct FilterAndTrimArgs {
    /// Forward (R1) input FASTQ file
    #[arg(long, required = true, help_heading = H_INPUT)]
    pub fwd: PathBuf,

    /// Reverse (R2) input FASTQ file; enables paired-end mode
    #[arg(long, help_heading = H_INPUT)]
    pub rev: Option<PathBuf>,

    /// Sample identifier for the output JSON [default: --fwd filename stem]
    #[arg(long, help_heading = H_INPUT)]
    pub sample_name: Option<String>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Truncate reads at the first Phred score <= this value
    #[arg(long, default_value = "2", num_args = 1..=2, help_heading = H_TRIM)]
    pub trunc_q: Vec<u8>,

    /// Truncate reads to this many bases, discarding shorter ones; 0 = off
    #[arg(long, default_value = "0", num_args = 1..=2, help_heading = H_TRIM)]
    pub trunc_len: Vec<usize>,

    /// Remove this many bases from the 5' end
    #[arg(long, default_value = "0", num_args = 1..=2, help_heading = H_TRIM)]
    pub trim_left: Vec<usize>,

    /// Remove this many bases from the 3' end
    #[arg(long, default_value = "0", num_args = 1..=2, help_heading = H_TRIM)]
    pub trim_right: Vec<usize>,

    /// Discard reads longer than this before trimming; 0 = no limit
    #[arg(long, default_value = "0", num_args = 1..=2, help_heading = H_FILTER)]
    pub max_len: Vec<usize>,

    /// Discard reads shorter than this after all trimming
    #[arg(long, default_value = "20", num_args = 1..=2, help_heading = H_FILTER)]
    pub min_len: Vec<usize>,

    /// Discard reads with more than this many Ns; 0 = discard any N
    #[arg(long, default_value_t = 0, help_heading = H_FILTER)]
    pub max_n: usize,

    /// Discard reads with any Phred score below this value; 0 = off
    #[arg(long, default_value_t = 0, help_heading = H_FILTER)]
    pub min_q: u8,

    /// Discard reads with expected errors above this; omit for no EE filter
    #[arg(long, num_args = 1..=2, help_heading = H_FILTER)]
    pub max_ee: Vec<f64>,

    /// FASTA of the phiX genome; matching reads are removed
    #[arg(long, help_heading = H_FILTER)]
    pub phix_genome: Option<PathBuf>,

    /// Discard reads with 2-mer Shannon richness below this; 0 = off
    #[arg(long, default_value = "0", num_args = 1..=2, help_heading = H_FILTER)]
    pub rm_lowcomplex: Vec<f64>,

    /// Threads for bgzf output compression; >1 enables bgzf
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Forward (R1) output FASTQ file
    #[arg(long, required = true, help_heading = H_OUTPUT)]
    pub filt: PathBuf,

    /// Reverse (R2) output FASTQ file; required when --rev is given
    #[arg(long, help_heading = H_OUTPUT)]
    pub filt_rev: Option<PathBuf>,

    /// Gzip-compress output files
    #[arg(long, default_value_t = true, help_heading = H_OUTPUT)]
    pub compress: bool,

    /// Write JSON summary here instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::MakeSequenceTable`].
#[derive(Args)]
pub struct MakeSequenceTableArgs {
    /// JSON files from `dada` (one per sample) or `merge-pairs` (multi-sample)
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    /// Comma-separated sample names, one per input; single-sample `dada` files only
    #[arg(long, value_delimiter = ',', help_heading = H_INPUT)]
    pub sample_names: Vec<String>,

    /// Discard ASVs shorter than this length (inclusive)
    #[arg(long, help_heading = H_FILTER)]
    pub min_len: Option<usize>,

    /// Discard ASVs longer than this length (inclusive)
    #[arg(long, help_heading = H_FILTER)]
    pub max_len: Option<usize>,

    /// Column order for sequences
    #[arg(long, default_value = "abundance", help_heading = H_OUTPUT,
          value_parser = ["abundance", "nsamples", "none"])]
    pub order_by: String,

    /// Hash algorithm used to generate sequence identifiers
    #[arg(long, default_value = "md5", help_heading = H_OUTPUT,
          value_parser = ["md5", "sha1"])]
    pub hash: String,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,
}

/// Arguments for [`Commands::RemoveBimeraDenovo`].
#[derive(Args)]
pub struct RemoveBimeraDenovoArgs {
    /// Sequence table JSON produced by `make-sequence-table`
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Bimera detection method
    #[arg(long, default_value = "consensus", help_heading = H_CHIMERA,
          value_parser = ["consensus", "pooled", "per-sample"])]
    pub method: String,

    /// Parent fold-abundance minimum [default: 1.5 consensus, 2 pooled/per-sample]
    #[arg(long, help_heading = H_CHIMERA)]
    pub min_fold_parent_over_abundance: Option<f64>,

    /// Parent abundance minimum [default: 2 consensus, 8 pooled/per-sample]
    #[arg(long, help_heading = H_CHIMERA)]
    pub min_parent_abundance: Option<u32>,

    /// Also flag sequences one mismatch/indel away from an exact bimera
    #[arg(long, default_value_t = false, help_heading = H_CHIMERA)]
    pub allow_one_off: bool,

    /// Minimum mismatches to parent required for one-off detection
    #[arg(long, default_value_t = 4, help_heading = H_CHIMERA)]
    pub min_one_off_parent_distance: usize,

    /// (consensus) Fraction of samples a sequence must be flagged in
    #[arg(long, default_value_t = 0.9, help_heading = H_CHIMERA)]
    pub min_sample_fraction: f64,

    /// (consensus) Unflagged samples to ignore in the fraction vote
    #[arg(long, default_value_t = 1, help_heading = H_CHIMERA)]
    pub ignore_n_negatives: u32,

    /// Maximum shift in the ends-free alignment to potential parents
    #[arg(long, default_value_t = 16, help_heading = H_ALIGN)]
    pub max_shift: i32,

    /// Match score for the parent alignment (R MATCH)
    #[arg(long = "match", default_value_t = 5, help_heading = H_ALIGN)]
    pub match_score: i16,

    /// Mismatch penalty for the parent alignment (R MISMATCH)
    #[arg(long, allow_hyphen_values = true, default_value_t = -4, help_heading = H_ALIGN)]
    pub mismatch: i16,

    /// Gap penalty for the parent alignment (R GAP_PENALTY)
    #[arg(long, allow_hyphen_values = true, default_value_t = -8, help_heading = H_ALIGN)]
    pub gap_p: i16,

    /// Pairwise aligner; `wfa2` requires a build with `--features wfa`
    #[arg(long, value_enum, help_heading = H_ALIGN)]
    pub align_backend: Option<AlignBackend>,

    /// Number of threads for parallel bimera detection
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    /// EXPERIMENTAL: WFA edit-budget cap, 0 = unbounded [default: 50]
    #[arg(long, help_heading = H_EXP)]
    pub wfa_max_edits: Option<i32>,
}

/// Arguments for [`Commands::ChimeraDiagnostics`].
#[derive(Args)]
pub struct ChimeraDiagnosticsArgs {
    /// Sequence table JSON from `make-sequence-table` or `remove-bimera-denovo`
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Minimum fold-difference in abundance for a sequence to be a parent
    #[arg(long, default_value_t = 1.5, help_heading = H_CHIMERA)]
    pub min_fold_parent_over_abundance: f64,

    /// Minimum abundance for a sequence to be a parent
    #[arg(long, default_value_t = 2, help_heading = H_CHIMERA)]
    pub min_parent_abundance: u32,

    /// Minimum distance to the nearest single parent to flag a suspect
    #[arg(long, default_value_t = 15, help_heading = H_CHIMERA)]
    pub trimera_min_parent_dist: usize,

    /// Minimum residual gap (bp) for a credible third segment
    #[arg(long, default_value_t = 20, help_heading = H_CHIMERA)]
    pub trimera_min_gap: usize,

    /// Maximum third-parent mismatch fraction across the gap
    #[arg(long, default_value_t = 0.10, help_heading = H_CHIMERA)]
    pub trimera_max_gap_error: f64,

    /// Minimum length (bp) of each end flank
    #[arg(long, default_value_t = 30, help_heading = H_CHIMERA)]
    pub trimera_min_flank: usize,

    /// Maximum shift in the ends-free alignment to potential parents
    #[arg(long, default_value_t = 16, help_heading = H_ALIGN)]
    pub max_shift: i32,

    /// Match score for the parent alignment (R MATCH)
    #[arg(long = "match", default_value_t = 5, help_heading = H_ALIGN)]
    pub match_score: i16,

    /// Mismatch penalty for the parent alignment (R MISMATCH)
    #[arg(long, allow_hyphen_values = true, default_value_t = -4, help_heading = H_ALIGN)]
    pub mismatch: i16,

    /// Gap penalty for the parent alignment (R GAP_PENALTY)
    #[arg(long, allow_hyphen_values = true, default_value_t = -8, help_heading = H_ALIGN)]
    pub gap_p: i16,

    /// Pairwise aligner; `wfa2` requires a build with `--features wfa`
    #[arg(long, value_enum, help_heading = H_ALIGN)]
    pub align_backend: Option<AlignBackend>,

    /// Number of threads for parallel diagnostics
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write TSV to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// EXPERIMENTAL: WFA edit-budget cap, 0 = unbounded [default: 50]
    #[arg(long, help_heading = H_EXP)]
    pub wfa_max_edits: Option<i32>,
}

/// Arguments for [`Commands::SeqTableToTsv`].
#[derive(Args)]
pub struct SeqTableToTsvArgs {
    /// Sequence table JSON from `make-sequence-table` or `remove-bimera-denovo`
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Keep only sequences present in >= this many samples (R PSEUDO_PREVALENCE)
    #[arg(long, help_heading = H_FILTER)]
    pub prevalence: Option<u32>,

    /// Keep only sequences with total abundance >= this (R PSEUDO_ABUNDANCE)
    #[arg(long, help_heading = H_FILTER)]
    pub min_abundance: Option<u64>,

    /// Write TSV to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,
}

/// Arguments for [`Commands::AssignTaxonomy`].
#[derive(Args)]
pub struct AssignTaxonomyArgs {
    /// Query sequences: FASTA (.fa/.fa.gz/.fasta) or sequence-table JSON
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Reference FASTA with semicolon-delimited taxonomy strings as headers
    #[arg(long, help_heading = H_INPUT)]
    pub ref_fasta: PathBuf,

    /// Minimum bootstrap confidence to assign a level, 0-100
    #[arg(long, default_value_t = 50, help_heading = H_TAX)]
    pub min_boot: u32,

    /// Also classify the reverse complement and keep the better orientation
    #[arg(long, help_heading = H_TAX)]
    pub try_rc: bool,

    /// Comma-separated names for taxonomic levels, applied in order
    #[arg(
        long,
        default_value = "Kingdom,Phylum,Class,Order,Family,Genus,Species",
        value_delimiter = ',',
        help_heading = H_TAX
    )]
    pub tax_levels: Vec<String>,

    /// RNG seed for reproducible bootstrap sampling
    #[arg(long, default_value_t = 0x9E37_79B9_7F4A_7C15, help_heading = H_TAX)]
    pub seed: u64,

    /// Number of threads for parallel query classification
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Include raw bootstrap counts in the output
    #[arg(long, help_heading = H_OUTPUT)]
    pub output_bootstraps: bool,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::AssignSpecies`].
#[derive(Args)]
pub struct AssignSpeciesArgs {
    /// Taxonomy JSON produced by `assign-taxonomy`
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Reference FASTA with ">ID genus species" headers
    #[arg(long, help_heading = H_INPUT)]
    pub ref_fasta: PathBuf,

    /// Max distinct species per query; 1 = unambiguous only, 0 = unlimited
    #[arg(long, default_value_t = 1, help_heading = H_TAX)]
    pub allow_multiple: usize,

    /// Also try the reverse complement of each query
    #[arg(long, help_heading = H_TAX)]
    pub try_rc: bool,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::TaxToTsv`].
#[derive(Args)]
pub struct TaxToTsvArgs {
    /// JSON file produced by `assign-taxonomy` or `assign-species`
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// String written for unassigned (null) taxonomy levels
    #[arg(long, default_value = "NA", help_heading = H_OUTPUT)]
    pub na_string: String,

    /// Write TSV to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,
}

/// Arguments for [`Commands::SeqTableToFasta`].
#[derive(Args)]
pub struct SeqTableToFastaArgs {
    /// JSON file produced by the `make-sequence-table` subcommand
    #[arg(help_heading = H_INPUT)]
    pub input: PathBuf,

    /// Keep only sequences present in >= this many samples (R PSEUDO_PREVALENCE)
    #[arg(long, help_heading = H_FILTER)]
    pub prevalence: Option<u32>,

    /// Keep only sequences with total abundance >= this (R PSEUDO_ABUNDANCE)
    #[arg(long, help_heading = H_FILTER)]
    pub min_abundance: Option<u64>,

    /// Write FASTA to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,
}

/// Arguments for [`Commands::Sample`].
#[derive(Args)]
pub struct SampleArgs {
    /// FASTQ files (.fastq, .fastq.gz, .fq, .fq.gz) to process
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Stop after this many total bases; whole samples are taken in order
    #[arg(long, default_value_t = 100_000_000, help_heading = H_FILTER)]
    pub nbases: u64,

    /// Process input files in random order (shuffles sample order only)
    #[arg(long, help_heading = H_FILTER)]
    pub randomize: bool,

    /// RNG seed for reproducible --randomize
    #[arg(long, help_heading = H_FILTER)]
    pub seed: Option<u64>,

    /// Number of threads for dereplication
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Directory for per-sample JSON files (created if absent)
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output_dir: PathBuf,

    /// Pretty-print the JSON; the default is compact (~34% smaller)
    #[arg(long, help_heading = H_OUTPUT)]
    pub pretty: bool,

    /// Gzip each per-sample JSON (`{sample}.json.gz`)
    #[arg(long, help_heading = H_OUTPUT)]
    pub gzip: bool,

    /// Print progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,
}

/// Arguments for [`Commands::ErrorsFromSample`].
#[derive(Args)]
pub struct ErrorsFromSampleArgs {
    /// Sample or derep JSON files (`.json` / `.json.gz`)
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    #[command(flatten)]
    pub fit: ErrModelFitArgs,

    #[command(flatten)]
    pub denoise: LearnDenoiseArgs,

    /// Number of threads for parallel sample processing
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print per-iteration progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    #[command(flatten)]
    pub diag: LearnDiagArgs,

    #[command(flatten)]
    pub experimental: ExperimentalArgs,
}

/// Arguments for [`Commands::LearnErrors`].
#[derive(Args)]
pub struct LearnErrorsArgs {
    /// FASTQ or derep/sample JSON files to learn from
    #[arg(required = true, help_heading = H_INPUT)]
    pub input: Vec<PathBuf>,

    /// Phred offset: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7)
    #[arg(long, default_value_t = 33, help_heading = H_INPUT)]
    pub phred_offset: u8,

    /// Stop after this many total bases; whole samples are taken in order
    #[arg(long, default_value_t = 100_000_000, help_heading = H_ERRMODEL)]
    pub nbases: u64,

    /// Process input files in random order (shuffles sample order only)
    #[arg(long, help_heading = H_ERRMODEL)]
    pub randomize: bool,

    /// RNG seed for reproducible --randomize
    #[arg(long, help_heading = H_ERRMODEL)]
    pub seed: Option<u64>,

    #[command(flatten)]
    pub fit: ErrModelFitArgs,

    #[command(flatten)]
    pub denoise: LearnDenoiseArgs,

    /// Number of threads for parallel sample processing
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write JSON to this file instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Emit compact (minified) JSON instead of pretty-printed
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,

    /// Print per-iteration progress to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    #[command(flatten)]
    pub diag: LearnDiagArgs,

    #[command(flatten)]
    pub experimental: ExperimentalArgs,
}

/// Arguments for [`Commands::KdistCalibrate`].
#[derive(Args)]
pub struct KdistCalibrateArgs {
    /// Derep JSON files, or `dada` output with --from-dada[-pooled]
    #[arg(required = true, help_heading = H_INPUT)]
    pub inputs: Vec<PathBuf>,

    /// With --from-dada: directory holding the derep JSONs that fed `dada`
    #[arg(long, help_heading = H_INPUT)]
    pub derep_dir: Option<PathBuf>,

    /// Pair within each sample instead of pooling all uniques into one set
    #[arg(long, help_heading = H_REGIME)]
    pub per_sample: bool,

    /// Link each unique to its nearest more-abundant neighbour; pair with --per-sample
    #[arg(long, help_heading = H_REGIME)]
    pub nearest_parent: bool,

    /// Post-inference: read `dada` output and label uniques by cluster role
    #[arg(long, help_heading = H_REGIME)]
    pub from_dada: bool,

    /// Post-inference: read a `dada-pooled` `_pooled.json[.gz]` record
    #[arg(long, conflicts_with = "from_dada", help_heading = H_REGIME)]
    pub from_dada_pooled: bool,

    /// Derive the minimizer cutoff matching the k-mer pass rate, then stop
    #[arg(long, help_heading = H_REGIME)]
    pub derive_cutoff: bool,

    /// With --derive-cutoff: sample pairs uniformly instead of abundance-weighted
    #[arg(long, help_heading = H_REGIME)]
    pub derive_uniform_pairs: bool,

    /// Max pairs per population; no effect under --nearest-parent
    #[arg(long, default_value_t = 200_000, conflicts_with = "nearest_parent", help_heading = H_REGIME)]
    pub max_pairs: usize,

    /// Subsample each sample to at most this many uniques; 0 = keep all
    #[arg(long, default_value_t = 0, help_heading = H_REGIME)]
    pub max_uniques: usize,

    /// RNG seed for reproducible subsampling
    #[arg(long, default_value_t = 0x9E37_79B9_7F4A_7C15, help_heading = H_REGIME)]
    pub seed: u64,

    /// K-mer size (R default 5; PacBio full-length wants 7)
    #[arg(long, default_value_t = 5, help_heading = H_SCREEN)]
    pub k: usize,

    /// Screen cutoff used for the `screened_in` flag and leakage summary
    #[arg(long, default_value_t = 0.42, help_heading = H_SCREEN)]
    pub cutoff: f64,

    /// Divergence above which a screened-in pair counts as leaked
    #[arg(long, default_value_t = 5.0, help_heading = H_SCREEN)]
    pub leak_pct: f64,

    /// Alignment band radius; negative = unbanded, and unbanded is correct here
    #[arg(long, default_value_t = -1, allow_negative_numbers = true, help_heading = H_ALIGN)]
    pub band: i32,

    /// Threads for the parallel alignment
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write CSV here instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Print per-population progress and the leakage summary to stderr
    #[arg(long, help_heading = H_DIAG)]
    pub verbose: bool,

    /// EXPERIMENTAL: which screen to calibrate; cutoffs do NOT transfer
    #[arg(long, value_enum, default_value_t = ScreenBackend::Kmer, help_heading = H_EXP)]
    pub screen_backend: ScreenBackend,

    /// EXPERIMENTAL (with `--screen-backend minimizer`): sketch k-mer size
    #[arg(long, default_value_t = crate::minimizers::MINIMIZER_K, help_heading = H_EXP)]
    pub minimizer_k: usize,

    /// EXPERIMENTAL (with `--screen-backend minimizer`): winnowing window
    #[arg(long, default_value_t = crate::minimizers::MINIMIZER_W, help_heading = H_EXP)]
    pub minimizer_w: usize,
}

/// Arguments for [`Commands::ReferenceEval`].
#[derive(Args)]
pub struct ReferenceEvalArgs {
    /// Query ASVs: a FASTA, or a `dada`/`dada-pooled` JSON (auto-detected)
    #[arg(long, help_heading = H_INPUT)]
    pub asvs: PathBuf,

    /// Reference/truth FASTA, e.g. mock-community alleles
    #[arg(long, help_heading = H_INPUT)]
    pub reference: PathBuf,

    /// Non-chimeric survivor set (JSON or FASTA) for the chimera 2x2
    #[arg(long, help_heading = H_INPUT)]
    pub non_chimeric: Option<PathBuf>,

    /// Bonferroni divisor `nraw`, so a Prior ASV reports p_a on the abundance scale
    #[arg(long, help_heading = H_INPUT)]
    pub nraw: Option<f64>,

    /// Edit distance at/below which an ASV is a true positive; 0 = exact
    #[arg(long, default_value_t = 0, help_heading = H_EVAL)]
    pub max_diffs: u32,

    /// Upper edit-distance bound for the report-only "near" bucket
    #[arg(long, default_value_t = 3, help_heading = H_EVAL)]
    pub near_diffs: u32,

    /// Match score for the NW alignment
    #[arg(long = "match", allow_hyphen_values = true, default_value_t = 5, help_heading = H_ALIGN)]
    pub match_score: i32,

    /// Mismatch score for the NW alignment
    #[arg(long, allow_hyphen_values = true, default_value_t = -4, help_heading = H_ALIGN)]
    pub mismatch: i32,

    /// Gap penalty for the NW alignment
    #[arg(long, allow_hyphen_values = true, default_value_t = -8, help_heading = H_ALIGN)]
    pub gap_p: i32,

    /// Alignment band radius; over-provisioned by default
    #[arg(long, allow_hyphen_values = true, default_value_t = 32, help_heading = H_ALIGN)]
    pub band: i32,

    /// Permissive k-mer prefilter cutoff before NW; omit to disable
    #[arg(long, help_heading = H_SCREEN)]
    pub kdist_screen: Option<f64>,

    /// K-mer size for the optional prefilter
    #[arg(long, default_value_t = 5, help_heading = H_SCREEN)]
    pub kmer_size: usize,

    /// Number of threads for alignment
    #[arg(long, default_value_t = 1, help_heading = H_PERF)]
    pub threads: usize,

    /// Write the summary JSON here instead of stdout
    #[arg(long, short = 'o', help_heading = H_OUTPUT)]
    pub output: Option<PathBuf>,

    /// Write the per-ASV classification table (TSV) here
    #[arg(long, help_heading = H_OUTPUT)]
    pub per_asv: Option<PathBuf>,

    /// Write the per-reference recovery table (TSV) here
    #[arg(long, help_heading = H_OUTPUT)]
    pub per_ref: Option<PathBuf>,

    /// Emit compact (minified) summary JSON
    #[arg(long, help_heading = H_OUTPUT)]
    pub compact: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    /// A subcommand missing from `COMMAND_GROUPS` would vanish from `--help`.
    #[test]
    fn every_subcommand_is_grouped_once() {
        // A debug build of the full command overflows the 2 MB test-thread stack.
        std::thread::Builder::new()
            .stack_size(16 << 20)
            .spawn(check_grouping)
            .unwrap()
            .join()
            .unwrap();
    }

    fn check_grouping() {
        let mut grouped: Vec<&str> = COMMAND_GROUPS
            .iter()
            .flat_map(|(_, names)| names.iter().copied())
            .collect();
        grouped.sort_unstable();
        let mut actual: Vec<String> = Cli::command()
            .get_subcommands()
            .filter(|s| !s.is_hide_set() && s.get_name() != "help")
            .map(|s| s.get_name().to_string())
            .collect();
        actual.sort_unstable();
        assert_eq!(grouped, actual);
    }
}
