use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};

use crate::misc::DADA2_RS_VERSION;

#[derive(Parser)]
#[command(
    about = "DADA2 toolkit",
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
    Summary {
        /// Input FASTQ file (uncompressed or gzipped)
        input: PathBuf,

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

    /// Denoise a sample using the DADA2 algorithm
    ///
    /// Accepts either a FASTQ file (dereplicated in memory) or a JSON file
    /// produced by the `derep` or `sample` subcommand (`.json` / `.json.gz`).
    /// Pre-dereplicated input avoids re-reading the FASTQ when iterating on
    /// parameters.  Outputs a JSON object describing the inferred ASVs.
    ///
    /// By default `err_out` from the error model file is used as the error
    /// matrix.  Pass `--use-err-in` to use `err_in` instead.
    Dada {
        /// Input file: FASTQ (uncompressed or gzipped) or a derep/sample JSON.
        input: PathBuf,

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
        /// homopolymer indels are the dominant error mode.
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// K-mer distance cutoff for the pre-alignment screen. Pairs with
        /// k-mer distance above this threshold are not aligned (matches R's
        /// `KDIST_CUTOFF`). Lower values screen more aggressively (faster,
        /// more false negatives); raise for divergent sequences.
        #[arg(long)]
        kdist_cutoff: Option<f64>,

        /// K-mer size used for the pre-alignment screen and for the Raw
        /// k-mer vectors (matches R's `KMER_SIZE`). 5 is the DADA2 default,
        /// tuned for 16S/ITS-length amplicons. Valid range: 3..=8.
        /// Memory scales as 4^k per Raw (k=5 → 1KB, k=7 → 16KB).
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

    /// Denoise multiple samples with full pooling (R DADA2 `pool=TRUE`)
    ///
    /// Each input may be a FASTQ file (uncompressed or gzipped) or a JSON file
    /// produced by the `derep` or `sample` subcommand (`.json` / `.json.gz`),
    /// independently per sample.  Per-sample uniques are merged into one
    /// combined table (abundances summed, qualities abundance-weighted-averaged),
    /// DADA2 is run once on the merged table, and one JSON file per sample is
    /// written into the output directory containing only the ASVs present in
    /// that sample.
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
        #[arg(long, allow_hyphen_values = true)]
        homo_gap_p: Option<i32>,

        /// K-mer distance cutoff for the pre-alignment screen.
        #[arg(long)]
        kdist_cutoff: Option<f64>,

        /// K-mer size used for the pre-alignment screen.
        #[arg(long)]
        kmer_size: Option<usize>,

        /// Disable the k-mer pre-alignment screen.
        #[arg(long)]
        no_kmer_screen: Option<bool>,

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

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

    /// Filter and trim a single sample's FASTQ reads
    ///
    /// Mirrors R's `filterAndTrim` function for a single sample.  Pass the
    /// forward (R1) input/output file pair; for paired-end data also supply
    /// `--rev` / `--filt-rev`.
    ///
    /// For parameters that accept paired values (`--trunc-len`, `--trim-left`,
    /// etc.) provide either one value (applied to both directions) or two
    /// space-separated values (first for forward, second for reverse).
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

    /// Convert a sequence table JSON to a tab-delimited count table
    ///
    /// Reads JSON produced by `make-sequence-table` or `remove-bimera-denovo`
    /// and writes a TSV with sequence IDs as rows and sample names as columns.
    ///
    /// Pass `--prevalence` and/or `--min-abundance` to filter rows the same
    /// way R DADA2's pseudo-pooling selects priors
    /// (`colSums(st>0) >= PSEUDO_PREVALENCE | colSums(st) >= PSEUDO_ABUNDANCE`).
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

    /// Convert a make-sequence-table JSON file to FASTA
    ///
    /// Writes one record per sequence using the sequence ID as the header.
    ///
    /// To extract pseudo-pooling priors (mirroring R DADA2's
    /// `pool="pseudo"` selection rule
    /// `colSums(st>0) >= PSEUDO_PREVALENCE | colSums(st) >= PSEUDO_ABUNDANCE`),
    /// pass `--prevalence` and/or `--min-abundance`.
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

        /// Output compact (minified) JSON instead of pretty-printed
        #[arg(long)]
        compact: bool,

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
    ErrorsFromSample {
        /// One or more derep JSON files produced by `sample` or `derep`
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Error model fitting function to use
        ///
        /// Allowed values: loess (default), noqual, binned-qual, pacbio, external
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

        /// External command to invoke when --errfun external is used.
        ///
        /// Whitespace-split into argv; the trans-input and err-output file
        /// paths are appended as the final two arguments. Both files use
        /// R's `read.table(..., row.names=1, header=TRUE, check.names=FALSE)`
        /// layout. See examples/external_errfun/ for reference scripts.
        #[arg(long)]
        errfun_cmd: Option<String>,

        /// Maximum self-consistency iterations (mirrors R's MAX_CONSIST)
        #[arg(long, default_value_t = 10)]
        max_consist: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long, default_value = "1e-40")]
        omega_a: f64,

        /// Significance threshold for omega_c (reads not corrected to any center)
        #[arg(long, default_value = "1e-40")]
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
        /// -8 = default for both Illumina and PacBio HiFi. Lower values (closer
        /// to 0) can help with older CLR or Nanopore data where homopolymer
        /// indels dominate.
        #[arg(long, default_value_t = -8, allow_hyphen_values = true)]
        homo_gap_p: i32,

        /// K-mer distance cutoff for the pre-alignment screen. Pairs with
        /// k-mer distance above this threshold are not aligned (matches R's
        /// `KDIST_CUTOFF`). Lower values screen more aggressively (faster,
        /// more false negatives); raise for divergent sequences.
        #[arg(long, default_value_t = 0.42)]
        kdist_cutoff: f64,

        /// K-mer size used for the pre-alignment screen and for the Raw
        /// k-mer vectors (matches R's `KMER_SIZE`). 5 is the DADA2 default,
        /// tuned for 16S/ITS-length amplicons. Valid range: 3..=8.
        /// Memory scales as 4^k per Raw (k=5 → 1KB, k=7 → 16KB).
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

    /// Learn an error model from FASTQ files
    ///
    /// Reads one or more FASTQ files, dereplicates and subsamples them on the
    /// fly up to `--nbases` total bases, then iteratively runs the DADA2
    /// algorithm and re-fits the chosen error model until self-consistency.
    ///
    /// Output is a JSON object with three flat 16 × nq matrices:
    ///   `trans`   — accumulated transition counts,
    ///   `err_in`  — error rates used in the final DADA run,
    ///   `err_out` — error rates estimated from `trans`.
    LearnErrors {
        /// One or more FASTQ files (.fastq, .fastq.gz, .fq, .fq.gz) to learn from
        #[arg(required = true)]
        input: Vec<PathBuf>,

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

        /// Error model fitting function to use
        ///
        /// Allowed values: loess (default), noqual, binned-qual, pacbio, external
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

        /// External command to invoke when --errfun external is used.
        ///
        /// Whitespace-split into argv; the trans-input and err-output file
        /// paths are appended as the final two arguments. Both files use
        /// R's `read.table(..., row.names=1, header=TRUE, check.names=FALSE)`
        /// layout. See examples/external_errfun/ for reference scripts.
        #[arg(long)]
        errfun_cmd: Option<String>,

        /// Maximum self-consistency iterations (mirrors R's MAX_CONSIST)
        #[arg(long, default_value_t = 10)]
        max_consist: usize,

        /// Significance threshold for abundance-based cluster splitting (omega_a)
        #[arg(long, default_value = "1e-40")]
        omega_a: f64,

        /// Significance threshold for omega_c (reads not corrected to any center)
        #[arg(long, default_value = "1e-40")]
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
        /// -8 = default for both Illumina and PacBio HiFi. Lower values (closer
        /// to 0) can help with older CLR or Nanopore data where homopolymer
        /// indels dominate.
        #[arg(long, default_value_t = -8, allow_hyphen_values = true)]
        homo_gap_p: i32,

        /// K-mer distance cutoff for the pre-alignment screen. Pairs with
        /// k-mer distance above this threshold are not aligned (matches R's
        /// `KDIST_CUTOFF`). Lower values screen more aggressively (faster,
        /// more false negatives); raise for divergent sequences.
        #[arg(long, default_value_t = 0.42)]
        kdist_cutoff: f64,

        /// K-mer size used for the pre-alignment screen and for the Raw
        /// k-mer vectors (matches R's `KMER_SIZE`). 5 is the DADA2 default,
        /// tuned for 16S/ITS-length amplicons. Valid range: 3..=8.
        /// Memory scales as 4^k per Raw (k=5 → 1KB, k=7 → 16KB).
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
}
