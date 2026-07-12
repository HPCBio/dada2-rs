use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};

use crate::misc::DADA2_RS_VERSION;
use crate::nwalign::AlignBackend;

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

#[derive(Subcommand)]
pub enum Commands {
    /// Compute per-position quality metrics from a FASTQ file
    #[command(display_order = 1)]
    Summary {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

        /// Sample identifier included in the output JSON's `sample` field.
        /// Defaults to the filename stem of the input FASTQ.
        #[arg(long)]
        sample_name: Option<String>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for parallel processing
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Also compute a per-read sequence-complexity histogram (DADA2's
        /// `seqComplexity`/`plotComplexity`): the effective number of k-mers
        /// per read, `exp(Shannon entropy)` over its k-mer counts.
        #[arg(long)]
        complexity: bool,

        /// K-mer size for the complexity calculation (DADA2 default 2). Only
        /// used when `--complexity` is set.
        #[arg(long, default_value_t = 2)]
        complexity_kmer_size: u8,

        /// Number of histogram bins for complexity, spanning `[0, 4^kmer_size]`
        /// (DADA2 `plotComplexity` default 100). Only used with `--complexity`.
        #[arg(long, default_value_t = 100)]
        complexity_bins: usize,

        /// Also compute per-position cumulative expected-error (EE) metrics:
        /// `Σ 10^(-Q/10)` along each read, aggregated across reads into
        /// mean/median/min/max/quartiles per position. Useful for judging
        /// `filter-and-trim` `maxEE`/truncation choices.
        #[arg(long)]
        expected_error: bool,

        /// Number of log-spaced histogram bins backing the EE quantiles. Only
        /// used when `--expected-error` is set.
        #[arg(long, default_value_t = 200)]
        ee_bins: usize,

        /// Maximum number of distinct quality values for the data to be judged
        /// "binned" (NovaSeq/NextSeq collapse Phred to a handful of levels).
        /// Continuous Illumina/HiFi data has dozens of distinct values.
        #[arg(long, default_value_t = 8)]
        binned_threshold: usize,

        /// Also print a human-readable metrics report to stderr (total
        /// sequences, quality/EE ranges, whether quality scores are binned and
        /// their levels). Stdout stays clean JSON for piping.
        #[arg(long)]
        report: bool,
    },

    /// Merge per-sample `summary` JSONs into a run-level quality/binning report
    ///
    /// Unions the observed quality-value distribution across all samples so rare
    /// bins missed by any single low-count sample (e.g. a seldom-called Q2) are
    /// recovered at the run level. Optionally validates the run against a
    /// user-declared bin scheme via `--expected-bins`.
    SummaryMerge {
        /// Per-sample `summary` JSON files (gzip ok; `-` reads stdin).
        inputs: Vec<PathBuf>,

        /// Write JSON output to this file instead of stdout.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed.
        #[arg(long)]
        compact: bool,

        /// Maximum number of distinct quality values for the run to be judged
        /// "binned". Matches the `summary` default.
        #[arg(long, default_value_t = 8)]
        binned_threshold: usize,

        /// Also print a human-readable run-level report to stderr. Stdout stays
        /// clean JSON for piping.
        #[arg(long)]
        report: bool,

        /// Declared bin levels to validate the run against, e.g. `2,12,24,40`.
        /// When set, observed ⊆ expected is "consistent"; any observed value
        /// outside the declared set is a violation (a wrong or silently changed
        /// scheme). Declared-but-unobserved bins are reported as informational.
        #[arg(long, value_delimiter = ',')]
        expected_bins: Vec<u8>,
    },

    /// Dereplicate sequences from a FASTQ file
    ///
    /// Produces the equivalent of the R dada2 `derep` class: a set of unique
    /// sequences with read counts, per-unique integer Phred quality sums
    /// (`qual_sum`; mean = sum / count), and a read-to-unique mapping.
    #[command(display_order = 4)]
    Derep {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

        /// Sample identifier embedded in the output JSON's `sample` field.
        /// Downstream subcommands (`dada`, `dada-pooled`) pick this up as a
        /// default when their own sample-name flag is omitted.
        /// Defaults to the filename stem of the input FASTQ.
        #[arg(long)]
        sample_name: Option<String>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for parallel processing
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Include the per-read mapping (read index → unique index) in the output
        #[arg(long)]
        show_map: bool,

        /// Write JSON output to this file instead of stdout. When the path ends
        /// in `.gz` the output is gzip-compressed (read back transparently).
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Pretty-print the JSON. Default output is compact (minified), which is
        /// ~34% smaller on disk; pass this for human-readable output.
        #[arg(long)]
        pretty: bool,

        /// Print progress information to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Denoise a sample using the DADA2 algorithm
    ///
    /// Accepts either a FASTQ file (dereplicated in memory) or a JSON file
    /// produced by the `derep` or `sample` subcommand (`.json` / `.json.gz`).
    /// Pre-dereplicated input avoids re-reading the FASTQ when iterating on
    /// parameters.  Outputs a JSON object describing the inferred ASVs.
    ///
    /// By default `err_out` from the error model file is used as the error
    /// matrix.  Pass `--use-err-in` to use `err_in` instead.
    #[command(display_order = 8)]
    Dada {
        /// One or more input files: FASTQ (uncompressed or gzipped) or a
        /// derep/sample JSON. With a single input the result is written to
        /// `--output`/`-o` (or stdout). With more than one input the samples are
        /// processed independently/serially (NOT pooled) and one `{sample}.json`
        /// per sample is written to `--output-dir` (which is then required).
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// JSON error model file produced by the `learn-errors` subcommand
        #[arg(long)]
        error_model: PathBuf,

        /// Use `err_in` from the error model instead of `err_out`
        #[arg(long)]
        use_err_in: bool,

        /// Sample identifier included in the output JSON.
        /// Defaults to the filename stem of the input FASTQ.
        #[arg(long)]
        sample_name: Option<String>,

        /// FASTA file of prior sequences (uncompressed or gzip-compressed).
        ///
        /// Each sequence in the file that matches a dereplicated unique exactly
        /// is flagged as a prior, making it immune to the abundance p-value
        /// filter. Prior-based splitting uses --omega-p instead of --omega-a.
        #[arg(long)]
        prior: Option<PathBuf>,

        /// Inherit any unspecified algorithm parameters (omega_*, min_*,
        /// detect_singletons, band, homo_gap_p, kdist_cutoff, kmer_size,
        /// no_kmer_screen) from the error model JSON's `params` block. Any
        /// flag passed explicitly on the CLI still wins. Without this flag,
        /// the built-in CLI defaults apply and a warning is emitted for each
        /// CLI value that disagrees with the err model's value.
        #[arg(long)]
        inherit_err_params: bool,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads (used for both dereplication and DADA2 comparisons)
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// (Multi-input only) Number of samples to denoise concurrently, each on
        /// its own `threads / sample-jobs` sub-pool. A single sample's comparison
        /// map is often too small to feed many threads, so fanning samples across
        /// smaller sub-pools keeps every core fed (~4 threads/sample is the sweet
        /// spot) and bounds memory to this many samples in flight. Defaults to
        /// round(threads / 4) (1 at <=4 threads = serial). Ignored for a single
        /// input. Dial it down for very large/complex samples if memory is tight.
        #[arg(long)]
        sample_jobs: Option<usize>,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long)]
        omega_a: Option<f64>,

        /// Significance threshold for reads not corrected to any center (omega_c)
        #[arg(long)]
        omega_c: Option<f64>,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long)]
        omega_p: Option<f64>,

        /// Minimum fold-enrichment above expected for cluster splitting
        #[arg(long)]
        min_fold: Option<f64>,

        /// Minimum Hamming distance required for cluster splitting
        #[arg(long)]
        min_hamming: Option<u32>,

        /// Minimum read abundance required for cluster splitting
        #[arg(long)]
        min_abund: Option<u32>,

        /// Use singleton detection (tri-state: omit to inherit / use default,
        /// `true` or `false` to set explicitly).
        #[arg(long)]
        detect_singletons: Option<bool>,

        /// Alignment band radius, matching R's `BAND_SIZE` parameter.
        /// 16 = Illumina default. 32 = recommended for PacBio HiFi 16S amplicons
        /// (per the DADA2 LRAS manuscript). -1 = unbanded (O(n²), rarely needed).
        #[arg(long, allow_hyphen_values = true)]
        band: Option<i32>,

        /// Homopolymer-run gap penalty. Matches R's `HOMOPOLYMER_GAP_PENALTY`;
        /// PacBio pipelines typically set this closer to 0 (e.g. `-1`) because
        /// homopolymer indels are the dominant error mode. Defaults to --gap-p
        /// when unset (R's HOMOPOLYMER_GAP_PENALTY = NULL).
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// Gap penalty for the Needleman-Wunsch alignment (R's GAP_PENALTY).
        #[arg(long, allow_hyphen_values = true)]
        gap_p: Option<i32>,

        /// Match score for the Needleman-Wunsch alignment (R's MATCH).
        #[arg(long = "match", allow_hyphen_values = true)]
        match_score: Option<i32>,

        /// Mismatch score for the Needleman-Wunsch alignment (R's MISMATCH).
        #[arg(long, allow_hyphen_values = true)]
        mismatch: Option<i32>,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend (wfa2lib-rs) — ASV-equivalent on tested
        /// Illumina and PacBio HiFi data, but alignments are not byte-identical.
        /// `wfa2` requires a build with `--features wfa`; a default build (and the
        /// published crate) errors when it is selected. See issue #63.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. WFA aborts a pair once it needs more than this many
        /// edits and falls back to the NW path for that pair (NW-identical there).
        /// The budget is an absolute edit count, NOT a fraction of read length:
        /// denoising only aligns near-identical reads (~99.9% identity), so real
        /// error-copies stay a few edits apart regardless of read length, while
        /// divergent non-error-copy pairs that slip past the k-mer screen are
        /// bounded. Ignored for the `nw` backend. 0 = unbounded. [default: 50]
        /// (Internally an edit budget E maps to a WFA cost of E·|gap_p|, e.g.
        /// 50·8 = 400 with default scoring; the DADA2RS_WFA_MAX_STEPS env
        /// override is specified in those raw cost units, not edits.)
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// Maximum number of clusters to infer (R's MAX_CLUST). 0 = unlimited.
        #[arg(long)]
        max_clust: Option<usize>,

        /// Use greedy clustering (R's GREEDY). Tri-state: omit to inherit / use
        /// default.
        #[arg(long)]
        greedy: Option<bool>,

        /// Use quality scores in the error model (R's USE_QUALS). Tri-state:
        /// omit to inherit / use default.
        #[arg(long)]
        use_quals: Option<bool>,

        /// K-mer distance cutoff for the pre-alignment screen. Pairs with
        /// k-mer distance above this threshold are not aligned (matches R's
        /// `KDIST_CUTOFF`). Lower values screen more aggressively (faster,
        /// more false negatives); raise for divergent sequences.
        #[arg(long)]
        kdist_cutoff: Option<f64>,

        /// K-mer size used for the pre-alignment screen and for the Raw
        /// k-mer vectors (matches R's `KMER_SIZE`). 5 is the DADA2 default,
        /// tuned for 16S/ITS-length amplicons. Valid range: 3..=8 (8 is the
        /// hard ceiling — k-mer indices must fit in u16). For PacBio HiFi do NOT
        /// use k=5: on ~1.4 kb reads the screen is a no-op there (~every pair is
        /// aligned, ~4–5× slower at scale). Use k=7 for speed or k=6 to cap
        /// memory; both give effectively identical ASVs. Memory scales as 4^k
        /// per Raw (k=5 → 1KB, k=6 → 4KB, k=7 → 16KB, k=8 → 64KB).
        #[arg(long)]
        kmer_size: Option<usize>,

        /// Disable the k-mer pre-alignment screen (every pair is aligned).
        /// Much slower; use only when the screen is wrongly filtering valid
        /// comparisons. Tri-state: omit to inherit / use default.
        #[arg(long)]
        no_kmer_screen: Option<bool>,

        /// Emit R-DADA2-parity per-cluster diagnostics in the output JSON:
        /// `cluster_stats` (n0/n1/nunq/birth_qave/post-hoc pval),
        /// `cluster_quality` (mean quality at each reference position),
        /// `birth_subs` (the substitutions that drove each cluster split),
        /// and `transitions` (16 × nq transition-by-quality matrix).
        ///
        /// Adds one alignment per Raw against its cluster center; only enable
        /// when you need the extra diagnostics.
        #[arg(long)]
        aux_outputs: bool,

        /// Write a single full cluster trace (clusters.json) to this file.
        ///
        /// Describes the final cluster structure: cluster centers, members
        /// with hamming/λ/pval, birth metadata. Useful for ASV-calling QC
        /// and one-off plots. See examples/cluster_trace/.
        #[arg(long)]
        cluster_trace: Option<PathBuf>,

        /// Skip the per-cluster `members` array in the trace; emit only
        /// cluster centers and birth metadata.
        #[arg(long)]
        trace_no_members: bool,

        /// In the trace, only include members with abundance >= this value.
        #[arg(long, default_value_t = 1)]
        trace_min_abund: u32,

        /// Write JSON output to this file instead of stdout (single input only)
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output directory for per-sample JSON files (required when more than
        /// one input is given; created if absent). One `{sample}.json` per
        /// sample, same convention as `dada-pooled`.
        #[arg(long)]
        output_dir: Option<PathBuf>,

        /// Write a TSV of uniques that failed to denoise (final p-value < omega_c,
        /// `map == null`) to this file. Tidy long format with a header:
        /// `sequence<TAB>sample<TAB>reads`, one row per failed unique per sample
        /// it appears in. Answers "what failed to denoise and how many reads did
        /// it cost?" without hand-joining the `map` to the derep uniques.
        #[arg(long)]
        failed_uniques: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Gzip-compress the per-sample JSON files (written as `{sample}.json.gz`)
        #[arg(long)]
        gzip: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Denoise multiple samples with full pooling (R DADA2 `pool=TRUE`)
    ///
    /// Each input may be a FASTQ file (uncompressed or gzipped) or a JSON file
    /// produced by the `derep` or `sample` subcommand (`.json` / `.json.gz`),
    /// independently per sample.  Per-sample uniques are merged into one
    /// combined table (abundances summed, qualities abundance-weighted-averaged),
    /// DADA2 is run once on the merged table, and one JSON file per sample is
    /// written into the output directory containing only the ASVs present in
    /// that sample.
    #[command(display_order = 9)]
    DadaPooled {
        /// One or more input files — FASTQ (.fastq/.fastq.gz) or derep/sample JSON
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// JSON error model file produced by the `learn-errors` subcommand
        #[arg(long)]
        error_model: PathBuf,

        /// Use `err_in` from the error model instead of `err_out`
        #[arg(long)]
        use_err_in: bool,

        /// FASTA file of prior sequences (uncompressed or gzip-compressed).
        /// Sequences in this file are flagged as priors in the merged unique
        /// table, exempt from the abundance p-value filter.
        #[arg(long)]
        prior: Option<PathBuf>,

        /// Inherit any unspecified algorithm parameters from the error model
        /// JSON's `params` block. See the `dada` subcommand for full semantics.
        #[arg(long)]
        inherit_err_params: bool,

        /// Sample names, one per input FASTQ. Defaults to filename stems
        /// (e.g. `sample1.fastq.gz` → `sample1`).
        #[arg(long, value_delimiter = ',')]
        sample_names: Option<Vec<String>>,

        /// Output directory for per-sample JSON files (created if absent)
        #[arg(long, short = 'o')]
        output_dir: PathBuf,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for dereplication and DADA2 comparisons
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long)]
        omega_a: Option<f64>,

        /// Significance threshold for reads not corrected to any center (omega_c)
        #[arg(long)]
        omega_c: Option<f64>,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long)]
        omega_p: Option<f64>,

        /// Minimum fold-enrichment above expected for cluster splitting
        #[arg(long)]
        min_fold: Option<f64>,

        /// Minimum Hamming distance required for cluster splitting
        #[arg(long)]
        min_hamming: Option<u32>,

        /// Minimum read abundance required for cluster splitting
        #[arg(long)]
        min_abund: Option<u32>,

        /// Use singleton detection (omit to inherit / use default).
        #[arg(long)]
        detect_singletons: Option<bool>,

        /// Alignment band radius (matches R's `BAND_SIZE`).
        #[arg(long, allow_hyphen_values = true)]
        band: Option<i32>,

        /// Homopolymer-run gap penalty (matches R's `HOMOPOLYMER_GAP_PENALTY`).
        /// Defaults to --gap-p when unset (R's HOMOPOLYMER_GAP_PENALTY = NULL).
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// Gap penalty for the Needleman-Wunsch alignment (R's GAP_PENALTY).
        #[arg(long, allow_hyphen_values = true)]
        gap_p: Option<i32>,

        /// Match score for the Needleman-Wunsch alignment (R's MATCH).
        #[arg(long = "match", allow_hyphen_values = true)]
        match_score: Option<i32>,

        /// Mismatch score for the Needleman-Wunsch alignment (R's MISMATCH).
        #[arg(long, allow_hyphen_values = true)]
        mismatch: Option<i32>,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend (wfa2lib-rs) — ASV-equivalent on tested
        /// Illumina and PacBio HiFi data, but alignments are not byte-identical.
        /// `wfa2` requires a build with `--features wfa`; a default build (and the
        /// published crate) errors when it is selected. See issue #63.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. WFA aborts a pair once it needs more than this many
        /// edits and falls back to the NW path for that pair (NW-identical there).
        /// The budget is an absolute edit count, NOT a fraction of read length:
        /// denoising only aligns near-identical reads (~99.9% identity), so real
        /// error-copies stay a few edits apart regardless of read length, while
        /// divergent non-error-copy pairs that slip past the k-mer screen are
        /// bounded. Ignored for the `nw` backend. 0 = unbounded. [default: 50]
        /// (Internally an edit budget E maps to a WFA cost of E·|gap_p|, e.g.
        /// 50·8 = 400 with default scoring; the DADA2RS_WFA_MAX_STEPS env
        /// override is specified in those raw cost units, not edits.)
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// Maximum number of clusters to infer (R's MAX_CLUST). 0 = unlimited.
        #[arg(long)]
        max_clust: Option<usize>,

        /// Use greedy clustering (R's GREEDY). Tri-state: omit to inherit / use
        /// default.
        #[arg(long)]
        greedy: Option<bool>,

        /// Use quality scores in the error model (R's USE_QUALS). Tri-state:
        /// omit to inherit / use default.
        #[arg(long)]
        use_quals: Option<bool>,

        /// K-mer distance cutoff for the pre-alignment screen.
        #[arg(long)]
        kdist_cutoff: Option<f64>,

        /// K-mer size used for the pre-alignment screen.
        #[arg(long)]
        kmer_size: Option<usize>,

        /// Disable the k-mer pre-alignment screen.
        #[arg(long)]
        no_kmer_screen: Option<bool>,

        /// Write a TSV of uniques that failed to denoise to this file. Tidy long
        /// format with a header: `sequence<TAB>sample<TAB>reads`. Because pooled
        /// denoising runs once on the merged unique table, "failed" is a global
        /// property (`result.map == null` on the merged index); for each failed
        /// merged unique a row is emitted per sample it appears in, carrying that
        /// sample's read count.
        #[arg(long)]
        failed_uniques: Option<PathBuf>,

        /// Write a self-contained pooled record (merged uniques with pooled
        /// abundance + the global map + global ASVs) to this path, for
        /// `kdist-calibrate --from-dada-pooled`. Off by default. Give an explicit
        /// path OUTSIDE `--output-dir` so it doesn't join the per-sample
        /// `*.json.gz` glob; gzip follows the path's `.gz` extension.
        #[arg(long)]
        pooled_record: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Gzip-compress the per-sample JSON files (written as `{sample}.json.gz`)
        #[arg(long)]
        gzip: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Denoise multiple samples with pseudo-pooling (R DADA2 `pool="pseudo"`)
    ///
    /// Two per-sample rounds. Round 1 denoises each sample independently (no
    /// priors). The ASVs from round 1 are pooled into a sequence table and a
    /// prior set is selected using R DADA2's PSEUDO_PREVALENCE / PSEUDO_ABUNDANCE
    /// rule (`--pseudo-prevalence` / `--pseudo-min-abundance`). Round 2 re-runs
    /// each sample with those priors flagged (routed through `--omega-p`). One
    /// `{sample}.json` per sample is written to `--output-dir`.
    #[command(display_order = 10)]
    DadaPseudo {
        /// One or more input files — FASTQ (.fastq/.fastq.gz) or derep/sample JSON
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// JSON error model file produced by the `learn-errors` subcommand
        #[arg(long)]
        error_model: PathBuf,

        /// Use `err_in` from the error model instead of `err_out`
        #[arg(long)]
        use_err_in: bool,

        /// Inherit any unspecified algorithm parameters from the error model
        /// JSON's `params` block. See the `dada` subcommand for full semantics.
        #[arg(long)]
        inherit_err_params: bool,

        /// Sample names, one per input FASTQ. Defaults to filename stems.
        #[arg(long, value_delimiter = ',')]
        sample_names: Option<Vec<String>>,

        /// Output directory for per-sample JSON files (created if absent)
        #[arg(long, short = 'o')]
        output_dir: PathBuf,

        /// Minimum number of samples a sequence must appear in to become a
        /// round-2 prior (R DADA2 PSEUDO_PREVALENCE). Equivalent to
        /// seq-table-to-fasta's --prevalence.
        #[arg(long, default_value_t = 2)]
        pseudo_prevalence: u32,

        /// Minimum total abundance across samples to become a prior (R DADA2
        /// PSEUDO_ABUNDANCE). Equivalent to seq-table-to-fasta's --min-abundance.
        #[arg(long)]
        pseudo_min_abundance: Option<u64>,

        /// Optional FASTA dump of the selected round-2 priors.
        #[arg(long)]
        priors_out: Option<PathBuf>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for dereplication and DADA2 comparisons
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Number of samples to denoise concurrently, each on its own
        /// `threads / sample-jobs` sub-pool. A single sample's comparison map is
        /// often too small to feed many threads, so fanning samples across
        /// smaller sub-pools keeps every core fed (~4 threads/sample is the sweet
        /// spot — the sample-jobs sweep's wall-time curve plateaus there).
        /// Defaults to round(threads / 4) (1 at <=4 threads, i.e. the original
        /// serial behavior). Trades a little peak memory (this many concurrent
        /// working sets) for much better thread utilization; dial it down for
        /// PacBio (k=7, larger per-sample state) if memory is tight.
        #[arg(long)]
        sample_jobs: Option<usize>,

        /// Cache every sample's uniques in memory across both pseudo-pooling
        /// rounds. By default dada-pseudo STREAMS: each sample is dropped after
        /// round 1 and re-read (re-dereplicated) in round 2, bounding peak
        /// memory to `--sample-jobs` samples in flight rather than all samples
        /// at once. Streaming is both faster and lighter on large runs (the
        /// retained all-samples cache is pure overhead — re-dereplication is
        /// cheaper than carrying it), so it is the default; pass --cache-samples
        /// to force the old all-in-memory behavior. Output is identical.
        #[arg(long)]
        cache_samples: bool,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long)]
        omega_a: Option<f64>,

        /// Significance threshold for reads not corrected to any center (omega_c)
        #[arg(long)]
        omega_c: Option<f64>,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long)]
        omega_p: Option<f64>,

        /// Minimum fold-enrichment above expected for cluster splitting
        #[arg(long)]
        min_fold: Option<f64>,

        /// Minimum Hamming distance required for cluster splitting
        #[arg(long)]
        min_hamming: Option<u32>,

        /// Minimum read abundance required for cluster splitting
        #[arg(long)]
        min_abund: Option<u32>,

        /// Use singleton detection (omit to inherit / use default).
        #[arg(long)]
        detect_singletons: Option<bool>,

        /// Alignment band radius (matches R's `BAND_SIZE`).
        #[arg(long, allow_hyphen_values = true)]
        band: Option<i32>,

        /// Homopolymer-run gap penalty (matches R's `HOMOPOLYMER_GAP_PENALTY`).
        /// Defaults to --gap-p when unset (R's HOMOPOLYMER_GAP_PENALTY = NULL).
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// Gap penalty for the Needleman-Wunsch alignment (R's GAP_PENALTY).
        #[arg(long, allow_hyphen_values = true)]
        gap_p: Option<i32>,

        /// Match score for the Needleman-Wunsch alignment (R's MATCH).
        #[arg(long = "match", allow_hyphen_values = true)]
        match_score: Option<i32>,

        /// Mismatch score for the Needleman-Wunsch alignment (R's MISMATCH).
        #[arg(long, allow_hyphen_values = true)]
        mismatch: Option<i32>,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend (wfa2lib-rs) — ASV-equivalent on tested
        /// Illumina and PacBio HiFi data, but alignments are not byte-identical.
        /// `wfa2` requires a build with `--features wfa`; a default build (and the
        /// published crate) errors when it is selected. See issue #63.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. WFA aborts a pair once it needs more than this many
        /// edits and falls back to the NW path for that pair (NW-identical there).
        /// The budget is an absolute edit count, NOT a fraction of read length:
        /// denoising only aligns near-identical reads (~99.9% identity), so real
        /// error-copies stay a few edits apart regardless of read length, while
        /// divergent non-error-copy pairs that slip past the k-mer screen are
        /// bounded. Ignored for the `nw` backend. 0 = unbounded. [default: 50]
        /// (Internally an edit budget E maps to a WFA cost of E·|gap_p|, e.g.
        /// 50·8 = 400 with default scoring; the DADA2RS_WFA_MAX_STEPS env
        /// override is specified in those raw cost units, not edits.)
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// Maximum number of clusters to infer (R's MAX_CLUST). 0 = unlimited.
        #[arg(long)]
        max_clust: Option<usize>,

        /// Use greedy clustering (R's GREEDY). Tri-state: omit to inherit / use
        /// default.
        #[arg(long)]
        greedy: Option<bool>,

        /// Use quality scores in the error model (R's USE_QUALS). Tri-state:
        /// omit to inherit / use default.
        #[arg(long)]
        use_quals: Option<bool>,

        /// K-mer distance cutoff for the pre-alignment screen.
        #[arg(long)]
        kdist_cutoff: Option<f64>,

        /// K-mer size used for the pre-alignment screen.
        #[arg(long)]
        kmer_size: Option<usize>,

        /// Disable the k-mer pre-alignment screen.
        #[arg(long)]
        no_kmer_screen: Option<bool>,

        /// Write a TSV of uniques that failed to denoise (per-sample `map == null`)
        /// to this file. Tidy long format with a header:
        /// `sequence<TAB>sample<TAB>reads`, one row per failed unique per sample.
        /// Failures are decided in the per-sample round-2 pass, so this is a
        /// per-sample signal (same semantics as `dada`).
        #[arg(long)]
        failed_uniques: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Gzip-compress the per-sample JSON files (written as `{sample}.json.gz`)
        #[arg(long)]
        gzip: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Merge denoised forward and reverse reads into full-length amplicons
    ///
    /// For each sample, the forward and reverse FASTQ files are re-dereplicated
    /// to reconstruct the read → unique mapping, which is composed with the
    /// unique → ASV mapping from the dada JSON files to count every
    /// (forward ASV, reverse ASV) pair.  Each distinct pair is then aligned
    /// (ends-free Needleman-Wunsch of the forward ASV against the
    /// reverse-complement of the reverse ASV) and accepted or rejected based on
    /// overlap length, mismatches, and indels.
    ///
    ///
    /// Files are matched by position: the first `--fwd-dada` corresponds to the
    /// first `--rev-dada`, `--fwd-fastq`, and `--rev-fastq`.  For hundreds of
    /// samples use shell globbing, e.g.:
    ///
    ///   dada2-rs merge-pairs \
    ///     --fwd-dada fwd_dada/*.json \
    ///     --rev-dada rev_dada/*.json \
    ///     --fwd-fastq fwd_fastq/*.fastq.gz \
    ///     --rev-fastq rev_fastq/*.fastq.gz
    #[command(display_order = 11)]
    MergePairs {
        /// Forward dada JSON files
        #[arg(long, required = true, num_args = 1..)]
        fwd_dada: Vec<PathBuf>,

        /// Reverse dada JSON files
        #[arg(long, required = true, num_args = 1..)]
        rev_dada: Vec<PathBuf>,

        /// Forward FASTQ files — re-dereplicated to recover read→unique mapping
        #[arg(long, required = true, num_args = 1..)]
        fwd_fastq: Vec<PathBuf>,

        /// Reverse FASTQ files — re-dereplicated to recover read→unique mapping
        #[arg(long, required = true, num_args = 1..)]
        rev_fastq: Vec<PathBuf>,

        /// Minimum overlap length between forward and RC(reverse) ASVs
        #[arg(long, default_value_t = 12)]
        min_overlap: u32,

        /// Maximum mismatches allowed in the overlap region
        #[arg(long, default_value_t = 0)]
        max_mismatch: u32,

        /// Include rejected merges (with `accept: false`) in the output
        #[arg(long)]
        return_rejects: bool,

        /// Concatenate forward and RC(reverse) with an N spacer instead of merging
        #[arg(long)]
        just_concatenate: bool,

        /// Rescue pairs that fail to merge by concatenating them, as in
        /// --just-concatenate, instead of dropping them. Useful for
        /// variable-length amplicons (e.g. ITS) whose reads may not overlap.
        /// Rescued reads are marked `concatenated: true` and this takes
        /// precedence over --return-rejects.
        #[arg(long)]
        rescue_unmerged: bool,

        /// Number of N characters in the concatenation spacer
        #[arg(long, default_value_t = 10)]
        concat_nnn_len: usize,

        /// Trim overhanging portions of forward/reverse reads past the overlap
        #[arg(long)]
        trim_overhang: bool,

        /// Override sample names (defaults to stems of --fwd-dada files)
        #[arg(long, num_args = 1..)]
        sample_names: Option<Vec<String>>,

        /// Verify that the fwd and rev dada JSONs carry the same `sample`
        /// field, that it matches the resolved sample name, and that both
        /// FASTQ filenames contain the sample name as a substring.
        #[arg(long)]
        check_sample_ids: bool,

        /// Phred quality-score offset for FASTQ re-dereplication
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads (used within each sample for dereplication)
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print per-sample progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Remove primer sequences from a FASTQ file
    ///
    /// Mirrors R's `removePrimers()`.  Detects and trims forward (and
    /// optionally reverse) primers from each read using mismatch-tolerant
    /// IUPAC-aware matching.  Reads lacking a primer match are discarded.
    ///
    /// With `--orient` (default), reads that match primers only in the
    /// reverse-complement direction are flipped before trimming.
    ///
    /// Outputs a trimmed FASTQ file; JSON stats (reads_in / reads_out) go to
    /// stdout or the file given by `-o`.
    #[command(display_order = 2)]
    RemovePrimers {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

        /// Output FASTQ file
        #[arg(long, short = 'f')]
        fout: PathBuf,

        /// Sample identifier included in the output JSON's `sample` field.
        /// Defaults to the filename stem of the input FASTQ.
        #[arg(long)]
        sample_name: Option<String>,

        /// Forward primer sequence (IUPAC ambiguity codes accepted,
        /// e.g. AGRGTTYGATYMTGGCTCAG)
        #[arg(long)]
        primer_fwd: String,

        /// Reverse primer sequence in its 5'→3' (catalog / synthesis) direction
        /// (IUPAC ambiguity codes accepted).  It is automatically
        /// reverse-complemented before matching (see `--rc-primer-rev`).
        /// Omit to skip reverse primer detection.
        #[arg(long)]
        primer_rev: Option<String>,

        /// Automatically reverse-complement `--primer-rev` before matching.
        /// Primers are conventionally specified 5'→3'; the reverse primer must
        /// be RC'd to match the orientation it appears in reads.  Pass
        /// `--rc-primer-rev false` only when supplying `--primer-rev` already
        /// as it appears in the read (i.e. already reverse-complemented).
        #[arg(long, default_value_t = true)]
        rc_primer_rev: bool,

        /// Maximum mismatches allowed when matching each primer
        #[arg(long, default_value_t = 2)]
        max_mismatch: usize,

        /// Allow insertions and deletions (indels) when matching primers, in
        /// addition to mismatches.  Uses Levenshtein edit distance where each
        /// mismatch or indel counts as 1 toward `--max-mismatch`.
        /// Significantly slower than the default mismatch-only mode.
        #[arg(long)]
        allow_indels: bool,

        /// Trim the forward primer from the 5′ end of each read
        #[arg(long, default_value_t = true)]
        trim_fwd: bool,

        /// Trim the reverse primer from the 3′ end of each read
        #[arg(long, default_value_t = true)]
        trim_rev: bool,

        /// Detect and correct read orientation: reads that match primers only
        /// in the reverse complement are flipped before trimming
        #[arg(long, default_value_t = true)]
        orient: bool,

        /// Gzip-compress the output FASTQ file
        #[arg(long, default_value_t = true)]
        compress: bool,

        /// Number of threads for parallel primer matching and bgzf output compression.
        /// Values > 1 enable bgzf (blocked gzip) output, which is valid gzip but seekable.
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Truncate reads at first Phred score ≤ this value (omit to disable).
        /// Applied after primer trimming.
        #[arg(long)]
        trunc_q: Option<u8>,

        /// Truncate reads to this many bases; discard if shorter (omit to disable).
        /// Applied after primer trimming.
        #[arg(long)]
        trunc_len: Option<usize>,

        /// Remove this many bases from the 5′ end of the primer-trimmed read.
        #[arg(long)]
        trim_left: Option<usize>,

        /// Remove this many bases from the 3′ end of the primer-trimmed read.
        #[arg(long)]
        trim_right: Option<usize>,

        /// Discard reads longer than this before quality trimming (omit = no limit).
        #[arg(long)]
        max_len: Option<usize>,

        /// Discard reads shorter than this after all trimming (omit = no minimum).
        #[arg(long)]
        min_len: Option<usize>,

        /// Discard reads with more than this many N bases.
        #[arg(long)]
        max_n: Option<usize>,

        /// Discard reads with any Phred score below this value (omit to disable).
        #[arg(long)]
        min_q: Option<u8>,

        /// Discard reads with expected errors above this threshold (omit to disable).
        #[arg(long)]
        max_ee: Option<f64>,

        /// Path to a FASTA file containing the phiX genome; reads matching it are removed.
        /// Omit to skip phiX filtering.
        #[arg(long)]
        phix_genome: Option<PathBuf>,

        /// Discard reads with 2-mer Shannon richness below this value (omit to disable).
        #[arg(long)]
        rm_lowcomplex: Option<f64>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3–1.7).
        /// Only relevant when quality-based filter options are used.
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Write JSON stats (reads_in / reads_out) to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Filter and trim a single sample's FASTQ reads
    ///
    /// Mirrors R's `filterAndTrim` function for a single sample.  Pass the
    /// forward (R1) input/output file pair; for paired-end data also supply
    /// `--rev` / `--filt-rev`.
    ///
    /// For parameters that accept paired values (`--trunc-len`, `--trim-left`,
    /// etc.) provide either one value (applied to both directions) or two
    /// space-separated values (first for forward, second for reverse).
    #[command(display_order = 3)]
    FilterAndTrim {
        /// Forward (R1) input FASTQ file
        #[arg(long, required = true)]
        fwd: PathBuf,

        /// Forward (R1) output FASTQ file
        #[arg(long, required = true)]
        filt: PathBuf,

        /// Reverse (R2) input FASTQ file (enables paired-end mode)
        #[arg(long)]
        rev: Option<PathBuf>,

        /// Reverse (R2) output FASTQ file (required when --rev is given)
        #[arg(long)]
        filt_rev: Option<PathBuf>,

        /// Sample identifier included in the output JSON.
        /// Defaults to the filename stem of --fwd.
        #[arg(long)]
        sample_name: Option<String>,

        /// Gzip-compress output files
        #[arg(long, default_value_t = true)]
        compress: bool,

        /// Number of threads for bgzf output compression.
        /// Values > 1 enable bgzf (blocked gzip) output, which is valid gzip but seekable.
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Truncate reads at first Phred score ≤ this value.
        /// One value (both directions) or two (fwd rev).
        #[arg(long, default_value = "2", num_args = 1..=2)]
        trunc_q: Vec<u8>,

        /// Truncate reads to this many bases; discard if shorter (0 = disabled).
        /// One value or two (fwd rev).
        #[arg(long, default_value = "0", num_args = 1..=2)]
        trunc_len: Vec<usize>,

        /// Remove this many bases from the 5′ end.
        /// One value or two (fwd rev).
        #[arg(long, default_value = "0", num_args = 1..=2)]
        trim_left: Vec<usize>,

        /// Remove this many bases from the 3′ end.
        /// One value or two (fwd rev).
        #[arg(long, default_value = "0", num_args = 1..=2)]
        trim_right: Vec<usize>,

        /// Discard reads longer than this before trimming (0 = no limit).
        /// One value or two (fwd rev).
        #[arg(long, default_value = "0", num_args = 1..=2)]
        max_len: Vec<usize>,

        /// Discard reads shorter than this after all trimming.
        /// One value or two (fwd rev).
        #[arg(long, default_value = "20", num_args = 1..=2)]
        min_len: Vec<usize>,

        /// Discard reads with more than this many N bases (0 = discard any N).
        #[arg(long, default_value_t = 0)]
        max_n: usize,

        /// Discard reads with any Phred score below this value (0 = disabled).
        #[arg(long, default_value_t = 0)]
        min_q: u8,

        /// Discard reads with expected errors above this threshold.
        /// One value or two (fwd rev). Omit for no EE filtering.
        #[arg(long, num_args = 1..=2)]
        max_ee: Vec<f64>,

        /// Path to a FASTA file containing the phiX genome; reads matching it are removed.
        /// Omit to skip phiX filtering.
        #[arg(long)]
        phix_genome: Option<PathBuf>,

        /// Discard reads with 2-mer Shannon richness below this value (0 = disabled).
        /// One value or two (fwd rev).
        #[arg(long, default_value = "0", num_args = 1..=2)]
        rm_lowcomplex: Vec<f64>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Write JSON summary to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Build a sample-by-sequence feature table
    ///
    /// Reads one or more JSON files produced by the `dada` or `merge-pairs`
    /// subcommands and assembles a flat count matrix (samples × sequences).
    #[command(display_order = 12)]
    MakeSequenceTable {
        /// One or more JSON files from `dada` (one file per sample) or
        /// `merge-pairs` (one file containing multiple samples).
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Sample name for each input file.
        ///
        /// Only applies to single-sample `dada` files; merge-pairs files carry
        /// sample names internally.  If provided, length must match --input.
        #[arg(long, num_args = 1..)]
        sample_names: Vec<String>,

        /// Order sequences (columns) by decreasing total abundance, number of
        /// samples present in, or leave in first-seen order.
        #[arg(long, default_value = "abundance",
              value_parser = ["abundance", "nsamples", "none"])]
        order_by: String,

        /// Discard ASV sequences shorter than this length (inclusive).
        /// Useful for removing off-target amplicons.
        #[arg(long)]
        min_len: Option<usize>,

        /// Discard ASV sequences longer than this length (inclusive).
        /// Useful for removing off-target amplicons.
        #[arg(long)]
        max_len: Option<usize>,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Hash algorithm used to generate sequence identifiers.
        #[arg(long, default_value = "md5",
              value_parser = ["md5", "sha1"])]
        hash: String,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,
    },

    /// Remove bimeric sequences from a sequence table
    ///
    /// Reads a JSON file produced by `make-sequence-table` and removes sequences
    /// identified as bimeras (chimeras of two more-abundant parents).
    /// Mirrors R's `removeBimeraDenovo`.
    #[command(display_order = 13)]
    RemoveBimeraDenovo {
        /// Sequence table JSON produced by `make-sequence-table`
        input: PathBuf,

        /// Bimera detection method
        #[arg(long, default_value = "consensus",
              value_parser = ["consensus", "pooled", "per-sample"])]
        method: String,

        /// Minimum fold-difference in abundance for a sequence to be a parent
        #[arg(long, default_value_t = 1.5)]
        min_fold_parent_over_abundance: f64,

        /// Minimum abundance for a sequence to be a parent
        #[arg(long, default_value_t = 2)]
        min_parent_abundance: u32,

        /// Also flag sequences one mismatch/indel away from an exact bimera
        #[arg(long, default_value_t = false)]
        allow_one_off: bool,

        /// Minimum mismatches to parent required for one-off bimera detection
        #[arg(long, default_value_t = 4)]
        min_one_off_parent_distance: usize,

        /// Maximum shift in ends-free alignment to potential parents
        #[arg(long, default_value_t = 16)]
        max_shift: i32,

        /// Match score for the parent alignment (mirrors R's `MATCH`)
        #[arg(long = "match", default_value_t = 5)]
        match_score: i16,

        /// Mismatch penalty for the parent alignment (mirrors R's `MISMATCH`)
        #[arg(long, allow_hyphen_values = true, default_value_t = -4)]
        mismatch: i16,

        /// Gap penalty for the parent alignment (mirrors R's `GAP_PENALTY`).
        /// R's `removeBimeraDenovo` honors the global setDadaOpt gap penalty.
        #[arg(long, allow_hyphen_values = true, default_value_t = -8)]
        gap_p: i16,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend (wfa2lib-rs) — ASV-equivalent on tested
        /// Illumina and PacBio HiFi data, but alignments are not byte-identical.
        /// `wfa2` requires a build with `--features wfa`; a default build (and the
        /// published crate) errors when it is selected. See issue #63.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. WFA aborts a pair once it needs more than this many
        /// edits and falls back to the NW path for that pair (NW-identical there).
        /// The budget is an absolute edit count, NOT a fraction of read length:
        /// denoising only aligns near-identical reads (~99.9% identity), so real
        /// error-copies stay a few edits apart regardless of read length, while
        /// divergent non-error-copy pairs that slip past the k-mer screen are
        /// bounded. Ignored for the `nw` backend. 0 = unbounded. [default: 50]
        /// (Internally an edit budget E maps to a WFA cost of E·|gap_p|, e.g.
        /// 50·8 = 400 with default scoring; the DADA2RS_WFA_MAX_STEPS env
        /// override is specified in those raw cost units, not edits.)
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// (consensus) Fraction of samples a sequence must be flagged in
        #[arg(long, default_value_t = 0.9)]
        min_sample_fraction: f64,

        /// (consensus) Number of unflagged samples to ignore in fraction vote
        #[arg(long, default_value_t = 1)]
        ignore_n_negatives: u32,

        /// Number of threads for parallel bimera detection
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,
    },

    /// Screen a sequence table for higher-order chimeras (trimeras)
    ///
    /// Reads JSON from `make-sequence-table` or `remove-bimera-denovo` and emits
    /// a per-sequence TSV of bimera *coverage* metrics. Unlike the boolean
    /// `remove-bimera-denovo` decision, this retains how much of each read a
    /// single two-parent junction explains: a read that survives bimera removal
    /// yet is nearly covered (`cover_frac` high) leaves a small internal gap a
    /// third parent can fill — a trimera suspect. Useful on long amplicons
    /// (full-length 16S, nodA) and low-biomass samples where complex chimeras
    /// are more likely. Coverage uses the pooled (across-sample) abundance model.
    #[command(display_order = 14)]
    ChimeraDiagnostics {
        /// Sequence table JSON produced by `make-sequence-table` or `remove-bimera-denovo`
        input: PathBuf,

        /// Minimum fold-difference in abundance for a sequence to be a parent
        #[arg(long, default_value_t = 1.5)]
        min_fold_parent_over_abundance: f64,

        /// Minimum abundance for a sequence to be a parent
        #[arg(long, default_value_t = 2)]
        min_parent_abundance: u32,

        /// Maximum shift in ends-free alignment to potential parents
        #[arg(long, default_value_t = 16)]
        max_shift: i32,

        /// Match score for the parent alignment (mirrors R's `MATCH`)
        #[arg(long = "match", default_value_t = 5)]
        match_score: i16,

        /// Mismatch penalty for the parent alignment (mirrors R's `MISMATCH`)
        #[arg(long, allow_hyphen_values = true, default_value_t = -4)]
        mismatch: i16,

        /// Gap penalty for the parent alignment (mirrors R's `GAP_PENALTY`)
        #[arg(long, allow_hyphen_values = true, default_value_t = -8)]
        gap_p: i16,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. 0 = unbounded. [default: 50]
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// Minimum distance to the nearest single parent to flag a trimera
        /// suspect. Rejects few-SNP variants of one abundant parent, which leave
        /// a coverage gap but are not chimeras.
        #[arg(long, default_value_t = 15)]
        trimera_min_parent_dist: usize,

        /// Minimum residual gap length (bp) for a credible third segment.
        /// Rejects one-off bimeras (gap of ~1 base).
        #[arg(long, default_value_t = 20)]
        trimera_min_gap: usize,

        /// Maximum third-parent mismatch fraction across the gap (clean fit)
        #[arg(long, default_value_t = 0.10)]
        trimera_max_gap_error: f64,

        /// Minimum length (bp) of each end flank. A 3-segment mosaic needs two
        /// substantial flanks; rejects tiny-flank divergent singletons.
        #[arg(long, default_value_t = 30)]
        trimera_min_flank: usize,

        /// Number of threads for parallel diagnostics
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Write TSV output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },

    /// Convert a sequence table JSON to a tab-delimited count table
    ///
    /// Reads JSON produced by `make-sequence-table` or `remove-bimera-denovo`
    /// and writes a TSV with sequence IDs as rows and sample names as columns.
    ///
    /// Pass `--prevalence` and/or `--min-abundance` to filter rows the same
    /// way R DADA2's pseudo-pooling selects priors
    /// (`colSums(st>0) >= PSEUDO_PREVALENCE | colSums(st) >= PSEUDO_ABUNDANCE`).
    #[command(display_order = 16)]
    SeqTableToTsv {
        /// Sequence table JSON produced by `make-sequence-table` or `remove-bimera-denovo`
        input: PathBuf,

        /// Keep only sequences present in at least this many samples
        /// (mirrors R DADA2's `PSEUDO_PREVALENCE`).  Omit to disable.
        #[arg(long)]
        prevalence: Option<u32>,

        /// Keep only sequences whose total abundance is at least this value
        /// (mirrors R DADA2's `PSEUDO_ABUNDANCE`).  Omit to disable.
        #[arg(long)]
        min_abundance: Option<u64>,

        /// Write TSV output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },

    /// Assign taxonomy to sequences using a Naive Bayes k-mer classifier
    ///
    /// Mirrors R's `assignTaxonomy`.  The query input may be a FASTA file or
    /// a sequence-table JSON produced by `make-sequence-table`.  The reference
    /// FASTA must use DADA2-formatted headers where the description is a
    /// semicolon-separated taxonomy string, e.g.
    ///
    ///   >Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;
    ///
    /// Output is a JSON object with a `levels` array and an `assignments`
    /// array — one entry per query — containing the sequence, its assigned
    /// taxonomy (null where confidence is below `--min-boot`), and optionally
    /// the raw bootstrap counts.
    #[command(display_order = 14)]
    AssignTaxonomy {
        /// Query sequences: FASTA (.fa/.fa.gz/.fasta) or sequence-table JSON
        input: PathBuf,

        /// Reference FASTA with semicolon-delimited taxonomy strings as headers
        #[arg(long)]
        ref_fasta: PathBuf,

        /// Minimum bootstrap confidence to assign a taxonomic level (0–100)
        #[arg(long, default_value_t = 50)]
        min_boot: u32,

        /// Also classify the reverse complement of each query and keep the
        /// better-scoring orientation
        #[arg(long)]
        try_rc: bool,

        /// Include raw bootstrap counts in the output
        #[arg(long)]
        output_bootstraps: bool,

        /// Comma-separated names for taxonomic levels (applied in order)
        #[arg(
            long,
            default_value = "Kingdom,Phylum,Class,Order,Family,Genus,Species",
            value_delimiter = ','
        )]
        tax_levels: Vec<String>,

        /// RNG seed for reproducible bootstrap sampling
        #[arg(long)]
        seed: Option<u64>,

        /// Number of threads for parallel query classification
        #[arg(long, default_value_t = 1)]
        threads: usize,

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

    /// Fill in the Species column of an assign-taxonomy JSON by exact match
    ///
    /// Mirrors R DADA2's `addSpecies()`.  Reads a JSON file produced by
    /// `assign-taxonomy`, runs exact-match species assignment against
    /// `--ref-fasta`, and writes a JSON file with the same shape.  The
    /// "Species" level is appended (or replaced if already present), and is
    /// only filled when the species reference's genus matches the query's
    /// assigned Genus level (when present), using R's `matchGenera` rules
    /// (exact, "Genus " prefix, or `Genus/…`/`…/Genus` split-genus forms).
    ///
    /// The reference FASTA must use the DADA2 species-assignment format
    /// where each header contains three whitespace-delimited fields:
    /// accession, genus, species, e.g.
    ///
    ///   >AY123456 Staphylococcus aureus
    #[command(display_order = 15)]
    AssignSpecies {
        /// Taxonomy JSON produced by `assign-taxonomy`
        input: PathBuf,

        /// Reference FASTA with ">ID genus species" headers
        #[arg(long)]
        ref_fasta: PathBuf,

        /// Maximum distinct species to return per query (0 = unlimited).
        /// With the default of 1 only unambiguous assignments are returned,
        /// matching R's `allowMultiple=FALSE`.
        #[arg(long, default_value_t = 1)]
        allow_multiple: usize,

        /// Also try the reverse complement of each query
        #[arg(long)]
        try_rc: bool,

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

    /// Convert an assign-taxonomy or assign-species JSON to a TSV table
    ///
    /// Emits one row per assignment with the sequence ID first, followed by
    /// one column per taxonomic level in the order they appear in the input
    /// JSON.  Unassigned levels are written as `NA` (matching R DADA2 output).
    #[command(display_order = 18)]
    TaxToTsv {
        /// JSON file produced by `assign-taxonomy` or `assign-species`
        input: PathBuf,

        /// String written for unassigned (null) taxonomy levels
        #[arg(long, default_value = "NA")]
        na_string: String,

        /// Write TSV output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },

    /// Convert a make-sequence-table JSON file to FASTA
    ///
    /// Writes one record per sequence using the sequence ID as the header.
    ///
    /// To extract pseudo-pooling priors (mirroring R DADA2's
    /// `pool="pseudo"` selection rule
    /// `colSums(st>0) >= PSEUDO_PREVALENCE | colSums(st) >= PSEUDO_ABUNDANCE`),
    /// pass `--prevalence` and/or `--min-abundance`.
    #[command(display_order = 17)]
    SeqTableToFasta {
        /// JSON file produced by the `make-sequence-table` subcommand
        input: PathBuf,

        /// Keep only sequences present in at least this many samples
        /// (mirrors R DADA2's `PSEUDO_PREVALENCE`, default 2 in R).
        /// Omit to disable the prevalence rule.
        #[arg(long)]
        prevalence: Option<u32>,

        /// Keep only sequences whose total abundance across samples is at
        /// least this value (mirrors R DADA2's `PSEUDO_ABUNDANCE`).  Omit to
        /// disable the abundance rule (equivalent to R's default of `Inf`).
        #[arg(long)]
        min_abundance: Option<u64>,

        /// Write FASTA output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },

    /// Dereplicate and subsample FASTQ files, writing one JSON file per sample
    ///
    /// Processes input FASTQ files in order (or shuffled when `--randomize` is
    /// set), dereplicating each one and writing a JSON file to `--output-dir`.
    /// Processing stops once the cumulative base count reaches `--nbases`.
    ///
    /// Each output file uses the same format as the `derep` subcommand and can
    /// be passed directly to `errors-from-sample`.
    #[command(display_order = 5)]
    Sample {
        /// One or more FASTQ files (.fastq, .fastq.gz, .fq, .fq.gz) to process
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Directory to write per-sample JSON files into (created if absent)
        #[arg(long, short = 'o')]
        output_dir: PathBuf,

        /// Stop after accumulating at least this many total bases across input files
        #[arg(long, default_value_t = 100_000_000)]
        nbases: u64,

        /// Process input files in random order instead of the supplied order
        #[arg(long)]
        randomize: bool,

        /// RNG seed for reproducible randomization (only used with --randomize)
        #[arg(long)]
        seed: Option<u64>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Number of threads for dereplication
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Pretty-print the JSON. Default output is compact (minified), which is
        /// ~34% smaller on disk; pass this for human-readable output.
        #[arg(long)]
        pretty: bool,

        /// Gzip each per-sample JSON (writes `{sample}.json.gz`). Read back
        /// transparently by downstream subcommands.
        #[arg(long)]
        gzip: bool,

        /// Print progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Learn an error model from pre-computed sample JSON files
    ///
    /// Reads one or more JSON files produced by the `sample` subcommand (or
    /// `derep`) and iteratively runs the DADA2 algorithm, re-fitting the
    /// chosen error model until self-consistency — a clean reimplementation of
    /// R's `learnErrors()`.
    ///
    /// Output is a JSON object with three flat 16 × nq matrices:
    ///   `trans`   — accumulated transition counts,
    ///   `err_in`  — error rates used in the final DADA run,
    ///   `err_out` — error rates estimated from `trans`.
    #[command(display_order = 7)]
    ErrorsFromSample {
        /// One or more derep JSON files produced by `sample` or `derep`
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Error model fitting function to use.
        ///
        /// Allowed values: loess (default), noqual, binned-qual, pacbio, external.
        ///
        /// Note on `loess` vs R DADA2: the native Rust loess is
        /// algorithmically equivalent to R's `loess(surface = "direct")` —
        /// bit-exact to machine precision on real data (validated against
        /// `examples/external_errfun/loess_reference_direct.R`). R DADA2's
        /// `loessErrfun`, however, calls `loess(...)` with R's default
        /// `surface = "interpolate"`, which fits the local polynomial at
        /// kd-tree vertices and interpolates between them. The two surfaces
        /// disagree by ~1e-3 absolute / ~4% relative at low-Q edges; the
        /// downstream ASV impact on dada2 inference is minimal (~1 read per
        /// sample on a 362-sample benchmark).
        ///
        /// For bit-for-bit parity with R DADA2's `loessErrfun`, use:
        ///   --errfun external --errfun-cmd "Rscript examples/external_errfun/loess_reference.R"
        /// See issue #14 for the full decomposition.
        #[arg(long, default_value = "loess",
              value_parser = ["loess", "noqual", "binned-qual", "pacbio", "external"])]
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

        /// External command to invoke when --errfun external is used.
        ///
        /// Whitespace-split into argv; the trans-input and err-output file
        /// paths are appended as the final two arguments. Both files use
        /// R's `read.table(..., row.names=1, header=TRUE, check.names=FALSE)`
        /// layout. See examples/external_errfun/ for reference scripts.
        #[arg(long)]
        errfun_cmd: Option<String>,

        /// LOESS configuration preset.  Resolves a bundle of related knobs
        /// (`--loess-surface`, `--loess-cell`, `--loess-max-rate`,
        /// `--loess-min-rate`); any of those flags passed explicitly overrides
        /// the preset's value for that knob.
        ///
        /// - `default`: surface=direct, max-rate=0.25, min-rate=1e-7 — the
        ///   historical dada2-rs behavior.
        /// - `r-dada2`: surface=interpolate, cell=0.2, max-rate=0.25,
        ///   min-rate=1e-7.  Mirrors R DADA2's `loessErrfun` — R's default
        ///   `loess()` surface plus the same `[1e-7, 0.25]` clamp R DADA2
        ///   applies after the fit (errorModels.R:53-56).
        ///
        /// Both presets clamp to the same range; they differ only in the
        /// fitting surface.
        #[arg(long, default_value = "default",
              value_parser = ["default", "r-dada2"])]
        loess_preset: String,

        /// LOESS fitting surface (overrides preset).
        /// Only applies to `--errfun loess` and `--errfun pacbio`; ignored by
        /// `noqual`, `binned-qual`, and `external`.
        ///
        /// `direct` evaluates the local polynomial at every query point
        /// (matches R `loess(surface = "direct")`).  `interpolate` builds a
        /// 1-D kd-tree partition, fits at each vertex, and blends with cubic
        /// Hermite at queries (matches R's default `loess()`).
        #[arg(long, value_parser = ["direct", "interpolate"])]
        loess_surface: Option<String>,

        /// Maximum fraction of observations allowed per kd-tree cell before
        /// it is subdivided.  Only used with `--loess-surface interpolate`
        /// (i.e. `loess` and `pacbio` errfuns only).
        /// Mirrors R `loess.control(cell = ...)`; R's default is 0.2.
        #[arg(long)]
        loess_cell: Option<f64>,

        /// Upper clamp applied to off-diagonal error rates after fitting.
        /// Applies to `loess`, `pacbio`, `noqual`, and `binned-qual`; ignored
        /// by `external`.  Both presets default to 0.25 (matching R DADA2).
        /// Set to `1.0` to disable.
        #[arg(long)]
        loess_max_rate: Option<f64>,

        /// Lower clamp applied to off-diagonal error rates after fitting.
        /// Applies to `loess`, `pacbio`, `noqual`, and `binned-qual`; ignored
        /// by `external`.  Both presets default to 1e-7 (matching R DADA2).
        /// Set to `0.0` to disable.
        #[arg(long)]
        loess_min_rate: Option<f64>,

        /// Maximum self-consistency iterations (mirrors R's MAX_CONSIST)
        #[arg(long, default_value_t = 10)]
        max_consist: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long, default_value = "1e-40")]
        omega_a: f64,

        /// Significance threshold for omega_c (reads not corrected to any center).
        /// Defaults to 0, matching R DADA2's `learnErrors()` (which hard-codes
        /// OMEGA_C=0 in its internal dada() calls, overriding the standard
        /// `dada()` default of 1e-40). Pass `--omega-c 1e-40` to use the
        /// standard `dada()` value instead.
        #[arg(long, default_value = "0")]
        omega_c: f64,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long, default_value = "1e-4")]
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

        /// Alignment band radius, matching R's `BAND_SIZE` parameter.
        /// 16 = Illumina default. 32 = recommended for PacBio HiFi 16S amplicons
        /// (per the DADA2 LRAS manuscript). -1 = unbanded (O(n²), rarely needed).
        #[arg(long, default_value_t = 16, allow_hyphen_values = true)]
        band: i32,

        /// Homopolymer-run gap penalty, matching R's `HOMOPOLYMER_GAP_PENALTY`.
        /// Defaults to --gap-p when unset (R's HOMOPOLYMER_GAP_PENALTY = NULL),
        /// i.e. -8 for both Illumina and PacBio HiFi. Lower values (closer to 0)
        /// can help with older CLR or Nanopore data where homopolymer indels
        /// dominate.
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// Gap penalty for the Needleman-Wunsch alignment (R's GAP_PENALTY).
        /// Defaults to -8 when unset.
        #[arg(long, allow_hyphen_values = true)]
        gap_p: Option<i32>,

        /// Match score for the Needleman-Wunsch alignment (R's MATCH).
        #[arg(long = "match", default_value_t = 5, allow_hyphen_values = true)]
        match_score: i32,

        /// Mismatch score for the Needleman-Wunsch alignment (R's MISMATCH).
        #[arg(long, default_value_t = -4, allow_hyphen_values = true)]
        mismatch: i32,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend (wfa2lib-rs) — ASV-equivalent on tested
        /// Illumina and PacBio HiFi data, but alignments are not byte-identical.
        /// `wfa2` requires a build with `--features wfa`; a default build (and the
        /// published crate) errors when it is selected. See issue #63.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. WFA aborts a pair once it needs more than this many
        /// edits and falls back to the NW path for that pair (NW-identical there).
        /// The budget is an absolute edit count, NOT a fraction of read length:
        /// denoising only aligns near-identical reads (~99.9% identity), so real
        /// error-copies stay a few edits apart regardless of read length, while
        /// divergent non-error-copy pairs that slip past the k-mer screen are
        /// bounded. Ignored for the `nw` backend. 0 = unbounded. [default: 50]
        /// (Internally an edit budget E maps to a WFA cost of E·|gap_p|, e.g.
        /// 50·8 = 400 with default scoring; the DADA2RS_WFA_MAX_STEPS env
        /// override is specified in those raw cost units, not edits.)
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// Maximum number of clusters to infer (R's MAX_CLUST). 0 = unlimited.
        #[arg(long, default_value_t = 0)]
        max_clust: usize,

        /// Use greedy clustering (R's GREEDY). Tri-state: omit for default
        /// (true), or set explicitly.
        #[arg(long)]
        greedy: Option<bool>,

        /// Use quality scores in the error model (R's USE_QUALS). Tri-state:
        /// omit for default (true), or set explicitly.
        #[arg(long)]
        use_quals: Option<bool>,

        /// K-mer distance cutoff for the pre-alignment screen. Pairs with
        /// k-mer distance above this threshold are not aligned (matches R's
        /// `KDIST_CUTOFF`). Lower values screen more aggressively (faster,
        /// more false negatives); raise for divergent sequences.
        #[arg(long, default_value_t = 0.42)]
        kdist_cutoff: f64,

        /// K-mer size used for the pre-alignment screen and for the Raw
        /// k-mer vectors (matches R's `KMER_SIZE`). 5 is the DADA2 default,
        /// tuned for 16S/ITS-length amplicons. Valid range: 3..=8 (8 is the
        /// hard ceiling — k-mer indices must fit in u16). For PacBio HiFi do NOT
        /// use k=5: on ~1.4 kb reads the screen is a no-op there (~every pair is
        /// aligned, ~4–5× slower at scale). Use k=7 for speed or k=6 to cap
        /// memory; both give effectively identical ASVs. Memory scales as 4^k
        /// per Raw (k=5 → 1KB, k=6 → 4KB, k=7 → 16KB, k=8 → 64KB).
        #[arg(long, default_value_t = 5)]
        kmer_size: usize,

        /// Disable the k-mer pre-alignment screen (every pair is aligned).
        /// Much slower; use only when the screen is wrongly filtering valid
        /// comparisons.
        #[arg(long)]
        no_kmer_screen: bool,

        /// Number of threads for parallel sample processing
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Directory to write per-iteration cluster diagnostics (iter_001.json, …)
        ///
        /// Each file contains cluster counts and birth-type breakdown per sample
        /// for that iteration. The directory is created if it does not exist.
        #[arg(long)]
        diag_dir: Option<PathBuf>,

        /// Directory to write full per-iteration cluster traces
        /// (cluster_iter_NNN_sample_NNN.json).
        ///
        /// Each file describes the full cluster structure for one sample at one
        /// iteration: cluster centers, members with their hamming distance, λ,
        /// expected reads, and abundance p-value, plus the err matrix used for
        /// that iteration. See examples/cluster_trace/ for plotting scripts.
        #[arg(long)]
        cluster_trace_dir: Option<PathBuf>,

        /// Skip the per-cluster `members` array in trace files; emit only
        /// cluster centers and birth metadata. Reduces trace size ~10×.
        #[arg(long)]
        trace_no_members: bool,

        /// In trace files, only include members with abundance >= this value.
        /// Default 1 (include all).
        #[arg(long, default_value_t = 1)]
        trace_min_abund: u32,

        /// Print per-iteration progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Learn an error model from FASTQ or derep/sample JSON files
    ///
    /// Reads one or more FASTQ files (dereplicated on the fly) or pre-computed
    /// derep/sample JSON files (`.json` / `.json.gz`, as written by the `derep`
    /// and `sample` subcommands), subsamples them up to `--nbases` total bases,
    /// then iteratively runs the DADA2 algorithm and re-fits the chosen error
    /// model until self-consistency.
    ///
    /// Output is a JSON object with three flat 16 × nq matrices:
    ///   `trans`   — accumulated transition counts,
    ///   `err_in`  — error rates used in the final DADA run,
    ///   `err_out` — error rates estimated from `trans`.
    ///
    /// CAVEAT (cross-sample diversity): accumulation is *sample-level*. Files
    /// are taken whole — in supplied order, or shuffled with `--randomize` — and
    /// each contributes all of its bases until the running total reaches
    /// `--nbases`. Reads are never subsampled *within* a file. With modern deep
    /// runs a single sample can supply the entire `--nbases` budget (e.g. 400k
    /// reads × 250 bp = 100M bases = the default budget), so the error model may
    /// be learned from the diversity of just one or a few samples. `--randomize`
    /// only shuffles which samples are drawn first; it does not guarantee
    /// representation across the run. Raise `--nbases`, pass a hand-picked set
    /// of inputs, or pre-`sample` each file to spread learning across samples.
    /// (Tracking: issue #68 — cross-sample diversity in learn-errors.)
    #[command(display_order = 6)]
    LearnErrors {
        /// FASTQ (.fastq/.fastq.gz/.fq/.fq.gz) or derep/sample JSON
        /// (.json/.json.gz) files to learn from
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Stop after accumulating at least this many total bases across input
        /// files. NOTE: accumulation is sample-level (whole files); a single
        /// deep sample can fill the entire budget — see the command help.
        #[arg(long, default_value_t = 100_000_000)]
        nbases: u64,

        /// Process input files in random order instead of the supplied order.
        /// NOTE: shuffles sample ORDER only — does not subsample reads within a
        /// sample, so it does not guarantee cross-sample diversity.
        #[arg(long)]
        randomize: bool,

        /// RNG seed for reproducible randomization (only used with --randomize)
        #[arg(long)]
        seed: Option<u64>,

        /// Phred quality score offset (33 for Sanger/Illumina 1.8+)
        #[arg(long, default_value_t = 33)]
        phred_offset: u8,

        /// Error model fitting function to use.
        ///
        /// Allowed values: loess (default), noqual, binned-qual, pacbio, external.
        ///
        /// Note on `loess` vs R DADA2: the native Rust loess is
        /// algorithmically equivalent to R's `loess(surface = "direct")` —
        /// bit-exact to machine precision on real data (validated against
        /// `examples/external_errfun/loess_reference_direct.R`). R DADA2's
        /// `loessErrfun`, however, calls `loess(...)` with R's default
        /// `surface = "interpolate"`, which fits the local polynomial at
        /// kd-tree vertices and interpolates between them. The two surfaces
        /// disagree by ~1e-3 absolute / ~4% relative at low-Q edges; the
        /// downstream ASV impact on dada2 inference is minimal (~1 read per
        /// sample on a 362-sample benchmark).
        ///
        /// For bit-for-bit parity with R DADA2's `loessErrfun`, use:
        ///   --errfun external --errfun-cmd "Rscript examples/external_errfun/loess_reference.R"
        /// See issue #14 for the full decomposition.
        #[arg(long, default_value = "loess",
              value_parser = ["loess", "noqual", "binned-qual", "pacbio", "external"])]
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

        /// External command to invoke when --errfun external is used.
        ///
        /// Whitespace-split into argv; the trans-input and err-output file
        /// paths are appended as the final two arguments. Both files use
        /// R's `read.table(..., row.names=1, header=TRUE, check.names=FALSE)`
        /// layout. See examples/external_errfun/ for reference scripts.
        #[arg(long)]
        errfun_cmd: Option<String>,

        /// LOESS configuration preset.  Resolves a bundle of related knobs
        /// (`--loess-surface`, `--loess-cell`, `--loess-max-rate`,
        /// `--loess-min-rate`); any of those flags passed explicitly overrides
        /// the preset's value for that knob.
        ///
        /// - `default`: surface=direct, max-rate=0.25, min-rate=1e-7 — the
        ///   historical dada2-rs behavior.
        /// - `r-dada2`: surface=interpolate, cell=0.2, max-rate=0.25,
        ///   min-rate=1e-7.  Mirrors R DADA2's `loessErrfun` — R's default
        ///   `loess()` surface plus the same `[1e-7, 0.25]` clamp R DADA2
        ///   applies after the fit (errorModels.R:53-56).
        ///
        /// Both presets clamp to the same range; they differ only in the
        /// fitting surface.
        #[arg(long, default_value = "default",
              value_parser = ["default", "r-dada2"])]
        loess_preset: String,

        /// LOESS fitting surface (overrides preset).
        /// Only applies to `--errfun loess` and `--errfun pacbio`; ignored by
        /// `noqual`, `binned-qual`, and `external`.
        ///
        /// `direct` evaluates the local polynomial at every query point
        /// (matches R `loess(surface = "direct")`).  `interpolate` builds a
        /// 1-D kd-tree partition, fits at each vertex, and blends with cubic
        /// Hermite at queries (matches R's default `loess()`).
        #[arg(long, value_parser = ["direct", "interpolate"])]
        loess_surface: Option<String>,

        /// Maximum fraction of observations allowed per kd-tree cell before
        /// it is subdivided.  Only used with `--loess-surface interpolate`
        /// (i.e. `loess` and `pacbio` errfuns only).
        /// Mirrors R `loess.control(cell = ...)`; R's default is 0.2.
        #[arg(long)]
        loess_cell: Option<f64>,

        /// Upper clamp applied to off-diagonal error rates after fitting.
        /// Applies to `loess`, `pacbio`, `noqual`, and `binned-qual`; ignored
        /// by `external`.  Both presets default to 0.25 (matching R DADA2).
        /// Set to `1.0` to disable.
        #[arg(long)]
        loess_max_rate: Option<f64>,

        /// Lower clamp applied to off-diagonal error rates after fitting.
        /// Applies to `loess`, `pacbio`, `noqual`, and `binned-qual`; ignored
        /// by `external`.  Both presets default to 1e-7 (matching R DADA2).
        /// Set to `0.0` to disable.
        #[arg(long)]
        loess_min_rate: Option<f64>,

        /// Maximum self-consistency iterations (mirrors R's MAX_CONSIST)
        #[arg(long, default_value_t = 10)]
        max_consist: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long, default_value = "1e-40")]
        omega_a: f64,

        /// Significance threshold for omega_c (reads not corrected to any center).
        /// Defaults to 0, matching R DADA2's `learnErrors()` (which hard-codes
        /// OMEGA_C=0 in its internal dada() calls, overriding the standard
        /// `dada()` default of 1e-40). Pass `--omega-c 1e-40` to use the
        /// standard `dada()` value instead.
        #[arg(long, default_value = "0")]
        omega_c: f64,

        /// Significance threshold for prior-sequence splitting (omega_p)
        #[arg(long, default_value = "1e-4")]
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

        /// Alignment band radius, matching R's `BAND_SIZE` parameter.
        /// 16 = Illumina default. 32 = recommended for PacBio HiFi 16S amplicons
        /// (per the DADA2 LRAS manuscript). -1 = unbanded (O(n²), rarely needed).
        #[arg(long, default_value_t = 16, allow_hyphen_values = true)]
        band: i32,

        /// Homopolymer-run gap penalty, matching R's `HOMOPOLYMER_GAP_PENALTY`.
        /// Defaults to --gap-p when unset (R's HOMOPOLYMER_GAP_PENALTY = NULL),
        /// i.e. -8 for both Illumina and PacBio HiFi. Lower values (closer to 0)
        /// can help with older CLR or Nanopore data where homopolymer indels
        /// dominate.
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// Gap penalty for the Needleman-Wunsch alignment (R's GAP_PENALTY).
        /// Defaults to -8 when unset.
        #[arg(long, allow_hyphen_values = true)]
        gap_p: Option<i32>,

        /// Match score for the Needleman-Wunsch alignment (R's MATCH).
        #[arg(long = "match", default_value_t = 5, allow_hyphen_values = true)]
        match_score: i32,

        /// Mismatch score for the Needleman-Wunsch alignment (R's MISMATCH).
        #[arg(long, default_value_t = -4, allow_hyphen_values = true)]
        mismatch: i32,

        /// Pairwise alignment backend. `nw` (default) is Needleman-Wunsch; `wfa2`
        /// is the experimental WFA backend (wfa2lib-rs) — ASV-equivalent on tested
        /// Illumina and PacBio HiFi data, but alignments are not byte-identical.
        /// `wfa2` requires a build with `--features wfa`; a default build (and the
        /// published crate) errors when it is selected. See issue #63.
        #[arg(long, value_enum)]
        align_backend: Option<AlignBackend>,

        /// EXPERIMENTAL (with `--align-backend wfa2`): WFA edit-budget cap, in
        /// edit operations. WFA aborts a pair once it needs more than this many
        /// edits and falls back to the NW path for that pair (NW-identical there).
        /// The budget is an absolute edit count, NOT a fraction of read length:
        /// denoising only aligns near-identical reads (~99.9% identity), so real
        /// error-copies stay a few edits apart regardless of read length, while
        /// divergent non-error-copy pairs that slip past the k-mer screen are
        /// bounded. Ignored for the `nw` backend. 0 = unbounded. [default: 50]
        /// (Internally an edit budget E maps to a WFA cost of E·|gap_p|, e.g.
        /// 50·8 = 400 with default scoring; the DADA2RS_WFA_MAX_STEPS env
        /// override is specified in those raw cost units, not edits.)
        #[arg(long)]
        wfa_max_edits: Option<i32>,

        /// Maximum number of clusters to infer (R's MAX_CLUST). 0 = unlimited.
        #[arg(long, default_value_t = 0)]
        max_clust: usize,

        /// Use greedy clustering (R's GREEDY). Tri-state: omit for default
        /// (true), or set explicitly.
        #[arg(long)]
        greedy: Option<bool>,

        /// Use quality scores in the error model (R's USE_QUALS). Tri-state:
        /// omit for default (true), or set explicitly.
        #[arg(long)]
        use_quals: Option<bool>,

        /// K-mer distance cutoff for the pre-alignment screen. Pairs with
        /// k-mer distance above this threshold are not aligned (matches R's
        /// `KDIST_CUTOFF`). Lower values screen more aggressively (faster,
        /// more false negatives); raise for divergent sequences.
        #[arg(long, default_value_t = 0.42)]
        kdist_cutoff: f64,

        /// K-mer size used for the pre-alignment screen and for the Raw
        /// k-mer vectors (matches R's `KMER_SIZE`). 5 is the DADA2 default,
        /// tuned for 16S/ITS-length amplicons. Valid range: 3..=8 (8 is the
        /// hard ceiling — k-mer indices must fit in u16). For PacBio HiFi do NOT
        /// use k=5: on ~1.4 kb reads the screen is a no-op there (~every pair is
        /// aligned, ~4–5× slower at scale). Use k=7 for speed or k=6 to cap
        /// memory; both give effectively identical ASVs. Memory scales as 4^k
        /// per Raw (k=5 → 1KB, k=6 → 4KB, k=7 → 16KB, k=8 → 64KB).
        #[arg(long, default_value_t = 5)]
        kmer_size: usize,

        /// Disable the k-mer pre-alignment screen (every pair is aligned).
        /// Much slower; use only when the screen is wrongly filtering valid
        /// comparisons.
        #[arg(long)]
        no_kmer_screen: bool,

        /// Number of threads for parallel sample processing
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// Write JSON output to this file instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

        /// Directory to write per-iteration cluster diagnostics (iter_001.json, …)
        ///
        /// Each file contains cluster counts and birth-type breakdown per sample
        /// for that iteration. The directory is created if it does not exist.
        #[arg(long)]
        diag_dir: Option<PathBuf>,

        /// Directory to write full per-iteration cluster traces
        /// (cluster_iter_NNN_sample_NNN.json).
        ///
        /// Each file describes the full cluster structure for one sample at one
        /// iteration: cluster centers, members with their hamming distance, λ,
        /// expected reads, and abundance p-value, plus the err matrix used for
        /// that iteration. See examples/cluster_trace/ for plotting scripts.
        #[arg(long)]
        cluster_trace_dir: Option<PathBuf>,

        /// Skip the per-cluster `members` array in trace files; emit only
        /// cluster centers and birth metadata. Reduces trace size ~10×.
        #[arg(long)]
        trace_no_members: bool,

        /// In trace files, only include members with abundance >= this value.
        /// Default 1 (include all).
        #[arg(long, default_value_t = 1)]
        trace_min_abund: u32,

        /// Print per-iteration progress to stderr
        #[arg(long)]
        verbose: bool,
    },

    /// Calibrate the k-mer screen: emit kdist vs true alignment divergence
    ///
    /// For sampled pairs of unique sequences, reports the k-mer distance
    /// (`KDIST_CUTOFF` screen metric) alongside the true UNBANDED ends-free
    /// alignment divergence, so the 0.42 cutoff (nominally ~10%, calibrated on
    /// Illumina 16S) can be checked per dataset / platform / k / pooling regime.
    /// Outputs CSV: sample,kdist,edits,core_len,pct_div,screened_in,ab_i,ab_j.
    ///
    /// PICKING A REGIME. The default pools all input uniques into one set
    /// (full-pool geometry baseline). For the abundance-aware modes below,
    /// prefer a biologically meaningful population instead:
    ///
    ///   * `--per-sample` for the independent (per-sample denoising) regime.
    ///   * `--from-dada-pooled` for a true pooled run — pass the
    ///     `_pooled.json[.gz]` record from `dada-pooled`, NOT the raw derep
    ///     union, so the pool is scored once with pooled abundances and no
    ///     cross-sample double-counting.
    ///
    /// EXAMPLES
    ///
    /// All-vs-all geometry baseline, per sample (fast: pairs capped by --max-pairs):
    ///
    ///   dada2-rs kdist-calibrate derep/*.json.gz --k 5 --per-sample --threads 24 -o kdist.allvall.csv
    ///
    /// Abundance mode (nearest more-abundant parent): ALWAYS pair with --per-sample. --nearest-parent scans ALL uniques and ignores --max-pairs, so pooling the raw dereps is both O(n^2) huge and cross-sample noise (a unique's "parent" may live in another sample):
    ///
    ///   dada2-rs kdist-calibrate derep/*.json.gz --k 5 --nearest-parent --per-sample --threads 24 -o kdist.abundance.csv
    ///
    /// Pooled regime, post-inference: score the merged pool as one population:
    ///
    ///   dada2-rs kdist-calibrate run_pooled.json.gz --k 5 --from-dada-pooled --threads 24 -o kdist.pooled.csv
    KdistCalibrate {
        /// One or more derep JSON files (`.json` / `.json.gz`)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// k-mer size (DADA2 default 5; PacBio full-length wants 7)
        #[arg(long, default_value_t = 5)]
        k: usize,

        /// Screen cutoff used for the `screened_in` flag / leakage summary
        #[arg(long, default_value_t = 0.42)]
        cutoff: f64,

        /// Divergence above which a screened-in pair is "leaked" (too far to be
        /// an error copy). Crude — the true ceiling is abundance-dependent.
        #[arg(long, default_value_t = 5.0)]
        leak_pct: f64,

        /// Alignment band radius; negative = unbanded (the correct, default
        /// choice — a band truncates the divergence of distant pairs).
        #[arg(long, default_value_t = -1, allow_negative_numbers = true)]
        band: i32,

        /// Max pairs computed PER population (random-subsample above this to
        /// bound the O(n^2) cost). Applies to all-pairs mode and the center
        /// pairs of --from-dada[-pooled]; it has NO effect under
        /// --nearest-parent (which scans every unique uncapped — bound that with
        /// --max-uniques instead), so passing both is rejected.
        #[arg(long, default_value_t = 200_000, conflicts_with = "nearest_parent")]
        max_pairs: usize,

        /// Randomly subsample each sample to at most this many uniques before
        /// pairing (0 = keep all). Unlike --max-pairs this DOES apply under
        /// --nearest-parent, so it is the lever for bounding that mode's uncapped
        /// O(n^2) scan.
        #[arg(long, default_value_t = 0)]
        max_uniques: usize,

        /// Compute pairs WITHIN each sample (per-sample / independent regime)
        /// instead of pooling all uniques into one set (full-pool regime)
        #[arg(long)]
        per_sample: bool,

        /// Abundance-aware mode: instead of random pairs, link each unique to
        /// its nearest MORE-abundant neighbour (candidate error-copy parent) and
        /// report the screen's headroom above the real error-copy distances.
        /// Output columns change to sample,ab,parent_ab,ab_ratio,kdist,edits,
        /// core_len,pct_div,screened_in.
        ///
        /// CAVEAT: this scans EVERY unique against its more-abundant prefix and
        /// does NOT honour --max-pairs, so it is O(n^2) in the population size.
        /// Pair it with --per-sample (small per-sample n, and each parent link is
        /// a genuine within-sample error copy); for a pooled run use
        /// --from-dada-pooled. Do NOT run it over a pooled raw-derep union.
        #[arg(long)]
        nearest_parent: bool,

        /// Post-inference mode: treat the positional inputs as `dada` output
        /// JSONs (not derep JSONs) and label every input unique by what
        /// denoising did to it — center (survived as an ASV), member (absorbed
        /// as an error copy), or failed (shed by the abundance test). Requires
        /// --derep-dir. Output columns change to sample,class,cluster,ab,
        /// center_ab,ab_ratio,birth_type,birth_pval,kdist,edits,core_len,
        /// pct_div,band_req,screened_in.
        #[arg(long)]
        from_dada: bool,

        /// Pooled post-inference mode: treat the positional inputs as the
        /// `_pooled.json[.gz]` record(s) written by `dada-pooled` and screen the
        /// merged unique table against the single global partition. Self-contained
        /// (no --derep-dir), so the pool is assessed as one population with pooled
        /// abundances — not the re-aggregated per-sample splits. Same output
        /// columns as --from-dada (sample = `__pooled__`).
        #[arg(long, conflicts_with = "from_dada")]
        from_dada_pooled: bool,

        /// With --from-dada: directory holding the derep JSONs that fed `dada`.
        /// Matched to each output by sample name: an exact `{sample}.json[.gz]`
        /// first, else a `{sample}.*.json[.gz]` file (e.g. a pipeline-renamed
        /// `{sample}.derep.R1.json.gz`); ambiguous prefixes are resolved by the
        /// derep JSON's own `sample` field.
        #[arg(long)]
        derep_dir: Option<PathBuf>,

        /// Threads for the parallel alignment
        #[arg(long, default_value_t = 1)]
        threads: usize,

        /// RNG seed for subsampling (reproducible)
        #[arg(long, default_value_t = 0x9E37_79B9_7F4A_7C15)]
        seed: u64,

        /// Write CSV here instead of stdout
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Print per-population progress + leakage summary to stderr
        #[arg(long)]
        verbose: bool,
    },
}
