#![allow(clippy::doc_overindented_list_items)]
use std::{fs::File, io, path::Path};

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rand::seq::SliceRandom as _;
use rayon::prelude::*;

mod chimera;
mod cli;
mod cluster;
mod cluster_trace;
mod containers;
mod dada;
mod derep;
mod error;
mod error_models;
mod evaluate;
mod filter;
mod filter_trim;
mod kmers;
mod learn_errors;
mod merge_pairs;
mod misc;
mod nwalign;
mod pval;
mod remove_bimera;
mod remove_primers;
mod sequence_table;
mod summary;
mod taxonomy;

use clap::CommandFactory;
use cli::{Cli, Commands};
use containers::BirthType;
use derep::dereplicate;
use filter_trim::{
    FilterParams, PairedFiles, WriteOptions, filter_paired, filter_single, read_fasta_first_seq,
};
use learn_errors::{
    ErrFun, LearnDiagOptions, LearnedErrParams, learn_errors, load_derep_samples,
    load_fastq_samples,
};
use misc::{DADA2_RS_VERSION, Tagged, read_fasta_records, read_tagged_json};
use nwalign::AlignParams;
use remove_bimera::{BimeraParams, Method, remove_bimera_denovo};
use remove_primers::{RemovePrimersParams, iupac_reverse_complement, remove_primers};
use sequence_table::{HashAlgo, OrderBy, SequenceTable, make_sequence_table};
use serde::Serialize;
use summary::process;
use taxonomy::{
    SpeciesHit, SpeciesOptions, SpeciesRef, TaxonomyOptions, TaxonomyRef, assign_species,
    assign_taxonomy,
};

use crate::error_models::{LoessConfig, LoessSurface};
use crate::misc::WithPath;

/// Resolve a [`LoessConfig`] from CLI inputs: preset + per-knob overrides.
/// `--loess-surface`, `--loess-cell`, `--loess-max-rate`, and `--loess-min-rate`
/// each override the preset's value for that knob if supplied.  `--loess-cell`
/// is ignored unless the resolved surface is `Interpolate`.
fn resolve_loess_config(
    preset: &str,
    surface: Option<&str>,
    cell: Option<f64>,
    max_rate: Option<f64>,
    min_rate: Option<f64>,
) -> LoessConfig {
    let base = match preset {
        "r-dada2" => LoessConfig::r_dada2(),
        _ => LoessConfig::default(),
    };
    let surface = match surface {
        Some("interpolate") => {
            let c = cell.unwrap_or(match base.surface {
                LoessSurface::Interpolate { cell } => cell,
                LoessSurface::Direct => 0.2,
            });
            LoessSurface::Interpolate { cell: c }
        }
        Some("direct") => LoessSurface::Direct,
        _ => match base.surface {
            LoessSurface::Interpolate { cell: base_cell } => LoessSurface::Interpolate {
                cell: cell.unwrap_or(base_cell),
            },
            LoessSurface::Direct => LoessSurface::Direct,
        },
    };
    LoessConfig {
        surface,
        max_error_rate: max_rate.unwrap_or(base.max_error_rate),
        min_error_rate: min_rate.unwrap_or(base.min_error_rate),
    }
}

/// Build a [`LearnedErrParams`] snapshot from the resolved errfun + dada/align
/// params, for embedding in the err-model JSON. Captures everything dada cares
/// about so a downstream invocation can validate or inherit.
fn build_learned_err_params(
    errfun: &ErrFun,
    max_consist: usize,
    dp: &dada::DadaParams,
    ap: &AlignParams,
) -> LearnedErrParams {
    let (errfun_name, errfun_pseudocount, errfun_bins, errfun_cmd, loess_cfg) = match errfun {
        ErrFun::Loess { config } => ("loess", None, None, None, Some(config)),
        ErrFun::Noqual {
            pseudocount,
            config,
        } => ("noqual", Some(*pseudocount), None, None, Some(config)),
        ErrFun::BinnedQual { bins, config } => {
            ("binned-qual", None, Some(bins.clone()), None, Some(config))
        }
        ErrFun::PacBio { config } => ("pacbio", None, None, None, Some(config)),
        ErrFun::External { command } => ("external", None, None, Some(command.clone()), None),
    };
    LearnedErrParams {
        errfun: errfun_name.to_string(),
        errfun_pseudocount,
        errfun_bins,
        errfun_cmd,
        loess: loess_cfg.map(Into::into),
        max_consist,
        omega_a: dp.omega_a,
        // Deliberately not embedded: learn-time `omega_c` (R default 0)
        // differs from dada-time (R default 1e-40) and must not transfer.
        omega_c: None,
        omega_p: dp.omega_p,
        min_fold: dp.min_fold,
        min_hamming: dp.min_hamming,
        min_abund: dp.min_abund,
        detect_singletons: dp.detect_singletons,
        use_quals: dp.use_quals,
        greedy: dp.greedy,
        match_score: ap.match_score,
        mismatch: ap.mismatch,
        gap_p: ap.gap_p,
        homo_gap_p: ap.homo_gap_p,
        use_kmers: ap.use_kmers,
        kdist_cutoff: ap.kdist_cutoff,
        kmer_size: ap.kmer_size,
        band: ap.band,
        vectorized: ap.vectorized,
        gapless: ap.gapless,
    }
}

#[derive(Serialize, Copy, Clone)]
struct DadaRunParams {
    omega_a: f64,
    omega_c: f64,
    omega_p: f64,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
    detect_singletons: bool,
    band: i32,
    homo_gap_p: i32,
    kdist_cutoff: f64,
    kmer_size: usize,
    use_kmers: bool,
    /// Denoising mode: false = independent per-sample (`dada`), true = full
    /// pooling across samples (`dada-pooled`, R DADA2 pool=TRUE).
    pool: bool,
    /// Number of unique input sequences flagged as priors (0 when no `--prior`).
    n_prior: usize,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    let command = match cli.command {
        Some(c) => c,
        None => {
            eprintln!("dada2-rs {DADA2_RS_VERSION}");
            eprintln!();
            Cli::command().print_help()?;
            return Ok(());
        }
    };

    match command {
        Commands::Summary {
            input,
            sample_name,
            phred_offset,
            threads,
            output,
            compact,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            let summary = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                process(
                    MultiGzDecoder::new(File::open(&input).with_path(&input)?),
                    phred_offset,
                    &pool,
                )?
            } else {
                process(File::open(&input).with_path(&input)?, phred_offset, &pool)?
            };

            #[derive(Serialize)]
            struct SummaryOutput {
                sample: String,
                total_reads: u64,
                mean_quality_per_position: Vec<f64>,
                reads_per_position: Vec<u64>,
                max_quality: usize,
                /// `quality_histogram[pos][q]` = count of reads with quality `q`
                /// at zero-based cycle `pos`. Inner length is `max_quality + 1`.
                quality_histogram: Vec<Vec<u64>>,
            }

            let sample = sample_name.unwrap_or_else(|| fastq_stem(&input));
            let (max_quality, quality_histogram) = summary.quality_histogram();
            let out = SummaryOutput {
                sample,
                total_reads: summary.total_reads,
                mean_quality_per_position: summary.mean_quality_per_position(),
                reads_per_position: summary.reads_per_position().to_vec(),
                max_quality,
                quality_histogram,
            };

            let tagged = Tagged::new("summary", out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::Derep {
            input,
            sample_name,
            phred_offset,
            threads,
            output,
            show_map,
            compact,
            verbose,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            let derep = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                dereplicate(
                    MultiGzDecoder::new(File::open(&input).with_path(&input)?),
                    phred_offset,
                    &pool,
                    verbose,
                )?
            } else {
                dereplicate(
                    File::open(&input).with_path(&input)?,
                    phred_offset,
                    &pool,
                    verbose,
                )?
            };

            #[derive(Serialize)]
            struct UniqueEntry<'a> {
                sequence: &'a str,
                count: u64,
                mean_quality: &'a [f64],
            }

            #[derive(Serialize)]
            struct DerepOutput<'a> {
                sample: &'a str,
                total_reads: usize,
                unique_sequences: usize,
                /// "abundance_desc" — produced by `dereplicate()`; lets dada /
                /// dada-pooled skip the defensive abundance sort on reload.
                sort_order: &'static str,
                uniques: Vec<UniqueEntry<'a>>,
                #[serde(skip_serializing_if = "Option::is_none")]
                map: Option<&'a [usize]>,
            }

            let mut uniq_entries = Vec::with_capacity(derep.uniques.len());
            for (i, (seq, count)) in derep.uniques.iter().enumerate() {
                let sequence = std::str::from_utf8(seq)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                uniq_entries.push(UniqueEntry {
                    sequence,
                    count: *count,
                    mean_quality: &derep.quals[i],
                });
            }

            let sample = sample_name.unwrap_or_else(|| fastq_stem(&input));
            let derep_out = DerepOutput {
                sample: &sample,
                total_reads: derep.map.len(),
                unique_sequences: derep.uniques.len(),
                sort_order: "abundance_desc",
                uniques: uniq_entries,
                map: if show_map { Some(&derep.map) } else { None },
            };

            let tagged = Tagged::new("derep", derep_out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::Dada {
            input,
            error_model,
            use_err_in,
            sample_name,
            prior,
            inherit_err_params,
            phred_offset,
            threads,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            band,
            homo_gap_p,
            kdist_cutoff,
            kmer_size,
            no_kmer_screen,
            aux_outputs,
            cluster_trace,
            trace_no_members,
            trace_min_abund,
            output,
            output_dir,
            compact,
            verbose,
        } => {
            // ---- Load uniques from FASTQ or a derep/sample JSON ----
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            // ---- Multi-input: independent per-sample (serial, NOT pooled) ----
            if input.len() > 1 {
                // Reject single-sample-only flags.
                let mut bad: Vec<&str> = Vec::new();
                if output.is_some() {
                    bad.push("--output/-o");
                }
                if sample_name.is_some() {
                    bad.push("--sample-name");
                }
                if aux_outputs {
                    bad.push("--aux-outputs");
                }
                if cluster_trace.is_some() {
                    bad.push("--cluster-trace");
                }
                if trace_no_members {
                    bad.push("--trace-no-members");
                }
                // `trace_min_abund` has a default of 1; only flag if non-default.
                if trace_min_abund != 1 {
                    bad.push("--trace-min-abund");
                }
                if !bad.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "{} {} single-sample-only and cannot be combined with multiple inputs; use --output-dir for per-sample output",
                            bad.join(", "),
                            if bad.len() == 1 { "is" } else { "are" },
                        ),
                    ));
                }
                let output_dir = output_dir.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "more than one input given: --output-dir is required (one {sample}.json per sample is written there)",
                    )
                })?;
                std::fs::create_dir_all(&output_dir)?;

                // Load the optional prior set once.
                let prior_set: Option<std::collections::HashSet<String>> = match prior {
                    Some(ref prior_path) => Some(
                        read_fasta_records(prior_path)
                            .with_path(prior_path)?
                            .into_iter()
                            .map(|(_, seq)| String::from_utf8_lossy(&seq).to_ascii_uppercase())
                            .collect(),
                    ),
                    None => None,
                };

                // Resolve parameters once (shared across samples).
                let resolved = resolve_dada_params(
                    &error_model,
                    use_err_in,
                    inherit_err_params,
                    threads,
                    false, // aux_outputs (rejected above for multi-input)
                    false, // pool
                    verbose,
                    omega_a,
                    omega_c,
                    omega_p,
                    min_fold,
                    min_hamming,
                    min_abund,
                    detect_singletons,
                    band,
                    homo_gap_p,
                    kdist_cutoff,
                    kmer_size,
                    no_kmer_screen,
                )?;

                for path in &input {
                    let (derep, json_sample) =
                        load_derep_for_dada(path, phred_offset, &pool, verbose)?;
                    let mut raw_inputs: Vec<dada::RawInput> = derep
                        .uniques
                        .into_iter()
                        .zip(derep.quals)
                        .map(|((seq, count), quals)| {
                            let sequence = String::from_utf8(seq).unwrap_or_else(|e| {
                                String::from_utf8_lossy(e.as_bytes()).into_owned()
                            });
                            dada::RawInput {
                                seq: sequence,
                                abundance: count as u32,
                                prior: false,
                                quals: Some(quals),
                            }
                        })
                        .collect();
                    if raw_inputs.is_empty() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("{}: no uniques found", path.display()),
                        ));
                    }
                    if let Some(ref set) = prior_set {
                        let n_marked = mark_priors(&mut raw_inputs, set);
                        if verbose {
                            eprintln!(
                                "[dada] {} of {} unique(s) marked as prior in {}",
                                n_marked,
                                raw_inputs.len(),
                                path.display(),
                            );
                        }
                    }
                    let sample = json_sample.unwrap_or_else(|| fastq_stem(path));
                    let json = denoise_and_serialize(
                        "dada",
                        &sample,
                        &raw_inputs,
                        &resolved.params,
                        &resolved.run,
                        &pool,
                        compact,
                        verbose,
                    )?;
                    let out_path = output_dir.join(format!("{sample}.json"));
                    std::fs::write(&out_path, &json)?;
                    if verbose {
                        eprintln!("[dada] wrote {}", out_path.display());
                    }
                }
                return Ok(());
            }

            // ---- Single input: preserve existing behavior byte-for-byte ----
            if output_dir.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--output-dir is only valid with multiple inputs; for a single input use --output/-o or stdout",
                ));
            }
            let input = &input[0];

            let (derep, json_sample) = load_derep_for_dada(input, phred_offset, &pool, verbose)?;

            let mut raw_inputs: Vec<dada::RawInput> = derep
                .uniques
                .into_iter()
                .zip(derep.quals)
                .map(|((seq, count), quals)| {
                    let sequence = String::from_utf8(seq)
                        .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
                    dada::RawInput {
                        seq: sequence,
                        abundance: count as u32,
                        prior: false,
                        quals: Some(quals),
                    }
                })
                .collect();

            if raw_inputs.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("{}: no uniques found", input.display()),
                ));
            }

            // ---- Mark prior sequences ----
            if let Some(ref prior_path) = prior {
                let prior_seqs: std::collections::HashSet<String> = read_fasta_records(prior_path)
                    .with_path(prior_path)?
                    .into_iter()
                    .map(|(_, seq)| String::from_utf8_lossy(&seq).to_ascii_uppercase())
                    .collect();
                let n_marked = mark_priors(&mut raw_inputs, &prior_seqs);
                if verbose {
                    eprintln!(
                        "[dada] {} of {} unique(s) marked as prior from {}",
                        n_marked,
                        raw_inputs.len(),
                        prior_path.display(),
                    );
                }
            }

            // ---- Load error model + resolve parameters (shared helper) ----
            let resolved = resolve_dada_params(
                &error_model,
                use_err_in,
                inherit_err_params,
                threads,
                aux_outputs,
                false, // pool
                verbose,
                omega_a,
                omega_c,
                omega_p,
                min_fold,
                min_hamming,
                min_abund,
                detect_singletons,
                band,
                homo_gap_p,
                kdist_cutoff,
                kmer_size,
                no_kmer_screen,
            )?;
            let dada_params = resolved.params;
            let run_params = resolved.run;
            let nq = resolved.nq;

            // ---- Run DADA2 ----
            let result = pool
                .install(|| dada::dada_uniques(&raw_inputs, &dada_params))
                .map_err(io::Error::other)?;

            if verbose {
                eprintln!(
                    "[dada] {} ASV(s) from {} unique input(s); {} aligns, {} shrouded",
                    result.clusters.len(),
                    raw_inputs.len(),
                    result.nalign,
                    result.nshroud,
                );
            }

            // ---- Optional cluster trace ----
            if let Some(ref trace_path) = cluster_trace {
                if let Some(parent) = trace_path.parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent)?;
                }
                let sample_name = input
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("dada")
                    .to_string();
                let trace_params = cluster_trace::TraceParams {
                    no_members: trace_no_members,
                    min_abund: trace_min_abund,
                };
                cluster_trace::write_trace(
                    trace_path,
                    &sample_name,
                    None, // no iteration: this is the final dada run
                    &raw_inputs,
                    &result,
                    Some(&dada_params.err_mat),
                    nq,
                    trace_params,
                    compact,
                )?;
                if verbose {
                    eprintln!("[dada] cluster trace written to {}", trace_path.display());
                }
            }

            // ---- Serialize output ----
            #[derive(Serialize)]
            struct AsvEntry {
                sequence: String,
                abundance: u32,
                birth_type: String,
                birth_pval: f64,
                birth_fold: f64,
                birth_e: f64,
            }

            #[derive(Serialize)]
            struct DadaStats {
                nalign: u32,
                nshroud: u32,
            }

            #[derive(Serialize)]
            struct ClusterStatJson {
                sequence: String,
                abundance: u32,
                n0: u32,
                n1: u32,
                nunq: u32,
                pval: f64,
                #[serde(skip_serializing_if = "Option::is_none")]
                birth_from: Option<usize>,
                birth_pval: f64,
                birth_fold: f64,
                birth_ham: u32,
                birth_e: f64,
                #[serde(skip_serializing_if = "Option::is_none")]
                birth_qave: Option<f64>,
            }

            #[derive(Serialize)]
            struct BirthSubJson {
                cluster: usize,
                pos: u16,
                nt0: char,
                nt1: char,
                #[serde(skip_serializing_if = "Option::is_none")]
                qual: Option<u8>,
            }

            #[derive(Serialize)]
            struct AuxJson {
                cluster_stats: Vec<ClusterStatJson>,
                cluster_quality: Vec<Vec<f64>>,
                cluster_quality_maxlen: usize,
                birth_subs: Vec<BirthSubJson>,
                transitions: Vec<u32>,
                transitions_ncol: usize,
            }

            #[derive(Serialize)]
            struct DadaOutput {
                sample: String,
                num_asvs: usize,
                total_reads: u32,
                asvs: Vec<AsvEntry>,
                stats: DadaStats,
                params: DadaRunParams,
                map: Vec<Option<usize>>,
                #[serde(skip_serializing_if = "Option::is_none")]
                aux: Option<AuxJson>,
            }

            let sample = sample_name
                .or(json_sample)
                .unwrap_or_else(|| fastq_stem(input));
            let total_reads: u32 = result.clusters.iter().map(|c| c.reads).sum();

            let asvs: Vec<AsvEntry> = result
                .clusters
                .iter()
                .map(|c| {
                    let sequence: String = c
                        .sequence
                        .iter()
                        .map(|&b| misc::nt_decode(b) as char)
                        .collect();
                    let birth_type = match &c.birth_type {
                        BirthType::Initial => "Initial",
                        BirthType::Abundance => "Abundance",
                        BirthType::Prior => "Prior",
                        BirthType::Singleton => "Singleton",
                    }
                    .to_string();
                    AsvEntry {
                        sequence,
                        abundance: c.reads,
                        birth_type,
                        birth_pval: c.birth_pval,
                        birth_fold: c.birth_fold,
                        birth_e: c.birth_e,
                    }
                })
                .collect();

            let aux_json = result.aux.as_ref().map(|a| {
                let cluster_stats_j = a
                    .cluster_stats
                    .iter()
                    .map(|c| {
                        let sequence: String = c
                            .sequence
                            .iter()
                            .map(|&b| misc::nt_decode(b) as char)
                            .collect();
                        ClusterStatJson {
                            sequence,
                            abundance: c.abundance,
                            n0: c.n0,
                            n1: c.n1,
                            nunq: c.nunq,
                            pval: c.pval,
                            birth_from: c.birth_from,
                            birth_pval: c.birth_pval,
                            birth_fold: c.birth_fold,
                            birth_ham: c.birth_ham,
                            birth_e: c.birth_e,
                            birth_qave: c.birth_qave,
                        }
                    })
                    .collect();
                let birth_subs_j = a
                    .birth_subs
                    .iter()
                    .map(|r| BirthSubJson {
                        cluster: r.cluster,
                        pos: r.pos,
                        nt0: r.nt0 as char,
                        nt1: r.nt1 as char,
                        qual: r.qual,
                    })
                    .collect();
                AuxJson {
                    cluster_stats: cluster_stats_j,
                    cluster_quality: a.cluster_quality.clone(),
                    cluster_quality_maxlen: a.cluster_quality_maxlen,
                    birth_subs: birth_subs_j,
                    transitions: a.transitions.clone(),
                    transitions_ncol: a.transitions_ncol,
                }
            });

            let out = DadaOutput {
                sample,
                num_asvs: asvs.len(),
                total_reads,
                asvs,
                stats: DadaStats {
                    nalign: result.nalign,
                    nshroud: result.nshroud,
                },
                params: DadaRunParams {
                    n_prior: raw_inputs.iter().filter(|r| r.prior).count(),
                    ..run_params
                },
                map: result.map,
                aux: aux_json,
            };

            let tagged = Tagged::new("dada", out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::DadaPooled {
            input,
            error_model,
            use_err_in,
            prior,
            inherit_err_params,
            sample_names,
            output_dir,
            phred_offset,
            threads,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            band,
            homo_gap_p,
            kdist_cutoff,
            kmer_size,
            no_kmer_screen,
            compact,
            verbose,
        } => {
            use std::collections::{HashMap, HashSet};

            let n_samples = input.len();
            if let Some(ref names) = sample_names
                && names.len() != n_samples
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--sample-names has {} entries but {} input file(s) given",
                        names.len(),
                        n_samples
                    ),
                ));
            }

            std::fs::create_dir_all(&output_dir)?;

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            // ---- Per-sample dereplication (or load from derep/sample JSON) ----
            // Phase timers (printed under --verbose): derep and merge are
            // serial, so their share quantifies the single-threaded front of a
            // pooled run vs the parallel `run_dada`.
            let t_derep = std::time::Instant::now();
            let mut dereps: Vec<derep::Derep> = Vec::with_capacity(n_samples);
            let mut json_samples: Vec<Option<String>> = Vec::with_capacity(n_samples);
            for path in &input {
                let (d, name) = load_derep_for_dada(path, phred_offset, &pool, verbose)?;
                dereps.push(d);
                json_samples.push(name);
            }
            let t_derep = t_derep.elapsed();

            // Resolve sample names: CLI override > JSON-embedded > filename stem.
            let sample_names: Vec<String> = match sample_names {
                Some(names) => names,
                None => input
                    .iter()
                    .zip(json_samples)
                    .map(|(p, js)| js.unwrap_or_else(|| fastq_stem(p)))
                    .collect(),
            };

            // ---- Merge across samples (abundance-weighted quality average) ----
            let t_merge = std::time::Instant::now();
            let mut seq_to_merged: HashMap<Vec<u8>, usize> = HashMap::new();
            let mut merged_seqs: Vec<Vec<u8>> = Vec::new();
            let mut merged_qual_sum: Vec<Vec<f64>> = Vec::new();
            let mut merged_total: Vec<u32> = Vec::new();
            let mut local_to_merged: Vec<Vec<usize>> = Vec::with_capacity(n_samples);

            for derep in &dereps {
                let mut local_map: Vec<usize> = Vec::with_capacity(derep.uniques.len());
                for ((seq, count), qual) in derep.uniques.iter().zip(derep.quals.iter()) {
                    let count_u32 = *count as u32;
                    let mu = match seq_to_merged.get(seq) {
                        Some(&i) => {
                            merged_total[i] += count_u32;
                            for (p, &q) in qual.iter().enumerate() {
                                merged_qual_sum[i][p] += q * count_u32 as f64;
                            }
                            i
                        }
                        None => {
                            let i = merged_seqs.len();
                            seq_to_merged.insert(seq.clone(), i);
                            merged_seqs.push(seq.clone());
                            merged_qual_sum
                                .push(qual.iter().map(|&q| q * count_u32 as f64).collect());
                            merged_total.push(count_u32);
                            i
                        }
                    };
                    local_map.push(mu);
                }
                local_to_merged.push(local_map);
            }

            let n_merged = merged_seqs.len();
            if verbose {
                eprintln!(
                    "[dada-pooled] {} sample(s) → {} merged unique(s), {} total reads",
                    n_samples,
                    n_merged,
                    merged_total.iter().sum::<u32>()
                );
            }

            // ---- Build merged RawInput list ----
            let mut raw_inputs: Vec<dada::RawInput> = (0..n_merged)
                .map(|i| {
                    let total = merged_total[i] as f64;
                    let mean_qual: Vec<f64> =
                        merged_qual_sum[i].iter().map(|&s| s / total).collect();
                    let sequence: String = String::from_utf8(merged_seqs[i].clone())
                        .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
                    dada::RawInput {
                        seq: sequence,
                        abundance: merged_total[i],
                        prior: false,
                        quals: Some(mean_qual),
                    }
                })
                .collect();

            if raw_inputs.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "All input FASTQ files contain no reads",
                ));
            }

            // ---- Mark prior sequences ----
            if let Some(ref prior_path) = prior {
                let prior_seqs: HashSet<String> = read_fasta_records(prior_path)
                    .with_path(prior_path)?
                    .into_iter()
                    .map(|(_, seq)| String::from_utf8_lossy(&seq).to_ascii_uppercase())
                    .collect();
                let n_marked = mark_priors(&mut raw_inputs, &prior_seqs);
                if verbose {
                    eprintln!(
                        "[dada-pooled] {} of {} merged unique(s) marked as prior from {}",
                        n_marked,
                        raw_inputs.len(),
                        prior_path.display(),
                    );
                }
            }

            // ---- Load error model + resolve parameters (shared helper) ----
            let resolved = resolve_dada_params(
                &error_model,
                use_err_in,
                inherit_err_params,
                threads,
                false, // aux_outputs
                true,  // pool
                verbose,
                omega_a,
                omega_c,
                omega_p,
                min_fold,
                min_hamming,
                min_abund,
                detect_singletons,
                band,
                homo_gap_p,
                kdist_cutoff,
                kmer_size,
                no_kmer_screen,
            )?;
            let dada_params = resolved.params;
            let mut run_params = resolved.run;
            run_params.n_prior = raw_inputs.iter().filter(|r| r.prior).count();

            // ---- Run DADA once on the merged table ----
            let t_merge = t_merge.elapsed();
            let t_dada = std::time::Instant::now();
            let result = pool
                .install(|| dada::dada_uniques(&raw_inputs, &dada_params))
                .map_err(io::Error::other)?;
            let t_dada = t_dada.elapsed();

            if verbose {
                eprintln!(
                    "[dada-pooled] {} ASV(s) from {} merged unique(s); {} aligns, {} shrouded",
                    result.clusters.len(),
                    raw_inputs.len(),
                    result.nalign,
                    result.nshroud,
                );
            }

            // ---- Per-sample output ----
            let t_output = std::time::Instant::now();
            #[derive(Serialize)]
            struct AsvEntry {
                sequence: String,
                abundance: u32,
                birth_type: String,
                birth_pval: f64,
                birth_fold: f64,
                birth_e: f64,
            }
            #[derive(Serialize)]
            struct DadaStats {
                nalign: u32,
                nshroud: u32,
            }
            #[derive(Serialize)]
            struct DadaOutput {
                sample: String,
                num_asvs: usize,
                total_reads: u32,
                asvs: Vec<AsvEntry>,
                stats: DadaStats,
                params: DadaRunParams,
                map: Vec<Option<usize>>,
            }

            for (s, sample_name) in sample_names.iter().enumerate() {
                // Sum per-cluster reads for this sample by walking its local uniques.
                let mut cluster_reads: Vec<u32> = vec![0u32; result.clusters.len()];
                for (lu, &mu) in local_to_merged[s].iter().enumerate() {
                    if let Some(c) = result.map[mu] {
                        cluster_reads[c] += dereps[s].uniques[lu].1 as u32;
                    }
                }

                // Filter to clusters present in this sample; renumber globally → locally.
                let mut global_to_local: Vec<Option<usize>> = vec![None; result.clusters.len()];
                let mut asvs: Vec<AsvEntry> = Vec::new();
                for (c, cluster) in result.clusters.iter().enumerate() {
                    if cluster_reads[c] == 0 {
                        continue;
                    }
                    global_to_local[c] = Some(asvs.len());
                    let sequence: String = cluster
                        .sequence
                        .iter()
                        .map(|&b| misc::nt_decode(b) as char)
                        .collect();
                    let birth_type = match &cluster.birth_type {
                        BirthType::Initial => "Initial",
                        BirthType::Abundance => "Abundance",
                        BirthType::Prior => "Prior",
                        BirthType::Singleton => "Singleton",
                    }
                    .to_string();
                    asvs.push(AsvEntry {
                        sequence,
                        abundance: cluster_reads[c],
                        birth_type,
                        birth_pval: cluster.birth_pval,
                        birth_fold: cluster.birth_fold,
                        birth_e: cluster.birth_e,
                    });
                }

                let total_reads: u32 = cluster_reads.iter().sum();

                // Per-sample local-unique → local-cluster map (mirrors single-sample dada).
                let map: Vec<Option<usize>> = (0..dereps[s].uniques.len())
                    .map(|lu| {
                        let mu = local_to_merged[s][lu];
                        result.map[mu].and_then(|c| global_to_local[c])
                    })
                    .collect();

                let n_asvs = asvs.len();
                let out = DadaOutput {
                    sample: sample_name.clone(),
                    num_asvs: n_asvs,
                    total_reads,
                    asvs,
                    stats: DadaStats {
                        nalign: result.nalign,
                        nshroud: result.nshroud,
                    },
                    params: run_params,
                    map,
                };

                let tagged = Tagged::new("dada-pooled", out);
                let json = if compact {
                    serde_json::to_string(&tagged)
                } else {
                    serde_json::to_string_pretty(&tagged)
                }
                .map_err(io::Error::other)?;

                let path = output_dir.join(format!("{sample_name}.json"));
                std::fs::write(&path, &json)?;
                if verbose {
                    eprintln!(
                        "[dada-pooled] wrote {} ({} ASV(s), {} reads)",
                        path.display(),
                        n_asvs,
                        total_reads
                    );
                }
            }
            if verbose {
                let t_output = t_output.elapsed();
                let total = (t_derep + t_merge + t_dada + t_output)
                    .as_secs_f64()
                    .max(1e-9);
                let pct = |d: std::time::Duration| 100.0 * d.as_secs_f64() / total;
                eprintln!(
                    "[dada-pooled] phase wall times: derep={:.1}s ({:.0}%, serial)  merge={:.1}s ({:.0}%, serial)  run_dada={:.1}s ({:.0}%, parallel)  output={:.1}s ({:.0}%, serial)",
                    t_derep.as_secs_f64(),
                    pct(t_derep),
                    t_merge.as_secs_f64(),
                    pct(t_merge),
                    t_dada.as_secs_f64(),
                    pct(t_dada),
                    t_output.as_secs_f64(),
                    pct(t_output),
                );
            }
        }

        Commands::DadaPseudo {
            input,
            error_model,
            use_err_in,
            inherit_err_params,
            sample_names,
            output_dir,
            pseudo_prevalence,
            pseudo_min_abundance,
            priors_out,
            phred_offset,
            threads,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            band,
            homo_gap_p,
            kdist_cutoff,
            kmer_size,
            no_kmer_screen,
            compact,
            verbose,
        } => {
            use std::collections::{HashMap, HashSet};

            let n_samples = input.len();
            if let Some(ref names) = sample_names
                && names.len() != n_samples
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--sample-names has {} entries but {} input file(s) given",
                        names.len(),
                        n_samples
                    ),
                ));
            }

            std::fs::create_dir_all(&output_dir)?;

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            // ---- Load + cache each sample's uniques once ----
            let mut sample_raws: Vec<Vec<dada::RawInput>> = Vec::with_capacity(n_samples);
            let mut json_samples: Vec<Option<String>> = Vec::with_capacity(n_samples);
            for path in &input {
                let (derep, json_sample) = load_derep_for_dada(path, phred_offset, &pool, verbose)?;
                let raws: Vec<dada::RawInput> = derep
                    .uniques
                    .into_iter()
                    .zip(derep.quals)
                    .map(|((seq, count), quals)| {
                        let sequence = String::from_utf8(seq)
                            .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
                        dada::RawInput {
                            seq: sequence,
                            abundance: count as u32,
                            prior: false,
                            quals: Some(quals),
                        }
                    })
                    .collect();
                if raws.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("{}: no uniques found", path.display()),
                    ));
                }
                sample_raws.push(raws);
                json_samples.push(json_sample);
            }

            // Resolve sample names: CLI override > JSON-embedded > filename stem.
            let sample_names: Vec<String> = match sample_names {
                Some(names) => names,
                None => input
                    .iter()
                    .zip(json_samples)
                    .map(|(p, js)| js.unwrap_or_else(|| fastq_stem(p)))
                    .collect(),
            };

            // Resolve parameters once (shared across both rounds and all samples).
            let resolved = resolve_dada_params(
                &error_model,
                use_err_in,
                inherit_err_params,
                threads,
                false, // aux_outputs
                false, // pool (pseudo is per-sample, not pooled)
                verbose,
                omega_a,
                omega_c,
                omega_p,
                min_fold,
                min_hamming,
                min_abund,
                detect_singletons,
                band,
                homo_gap_p,
                kdist_cutoff,
                kmer_size,
                no_kmer_screen,
            )?;

            // ---- Round 1: denoise each sample independently (no priors) ----
            if verbose {
                eprintln!("[dada-pseudo] round 1: {n_samples} sample(s), no priors");
            }
            // Collect per-sample round-1 ASV (sequence, abundance) pairs.
            let mut round1_asvs: Vec<Vec<(String, u32)>> = Vec::with_capacity(n_samples);
            for (s, raws) in sample_raws.iter().enumerate() {
                let result = pool
                    .install(|| dada::dada_uniques(raws, &resolved.params))
                    .map_err(io::Error::other)?;
                let asvs: Vec<(String, u32)> = result
                    .clusters
                    .iter()
                    .map(|c| {
                        let sequence: String = c
                            .sequence
                            .iter()
                            .map(|&b| misc::nt_decode(b) as char)
                            .collect();
                        (sequence, c.reads)
                    })
                    .collect();
                if verbose {
                    eprintln!(
                        "[dada-pseudo]   round 1 {}: {} ASV(s)",
                        sample_names[s],
                        asvs.len(),
                    );
                }
                round1_asvs.push(asvs);
            }

            // ---- Build a sequence table from round-1 ASVs (samples × sequences) ----
            let mut seq_index: HashMap<String, usize> = HashMap::new();
            let mut sequences: Vec<String> = Vec::new();
            for asvs in &round1_asvs {
                for (seq, _) in asvs {
                    if !seq_index.contains_key(seq) {
                        seq_index.insert(seq.clone(), sequences.len());
                        sequences.push(seq.clone());
                    }
                }
            }
            let nseq = sequences.len();
            let mut counts: Vec<Vec<u64>> = vec![vec![0u64; nseq]; n_samples];
            for (s, asvs) in round1_asvs.iter().enumerate() {
                for (seq, abund) in asvs {
                    let j = seq_index[seq];
                    counts[s][j] += *abund as u64;
                }
            }
            let table = SequenceTable {
                samples: sample_names.clone(),
                sequence_ids: sequences.iter().map(|s| HashAlgo::Sha1.digest(s)).collect(),
                sequences,
                counts,
            };

            // ---- Select round-2 priors via R's PSEUDO rule ----
            let selected = select_sequences(&table, Some(pseudo_prevalence), pseudo_min_abundance);
            let prior_set: HashSet<String> = selected
                .iter()
                .map(|&j| table.sequences[j].to_ascii_uppercase())
                .collect();

            if verbose {
                eprintln!(
                    "[dada-pseudo] selected {} prior sequence(s) from {} round-1 ASV(s) (prevalence>={}{})",
                    prior_set.len(),
                    table.sequences.len(),
                    pseudo_prevalence,
                    match pseudo_min_abundance {
                        Some(m) => format!(" or total-abundance>={m}"),
                        None => String::new(),
                    },
                );
            }

            // ---- Optional priors FASTA dump ----
            if let Some(ref priors_path) = priors_out {
                if let Some(parent) = priors_path.parent()
                    && !parent.as_os_str().is_empty()
                {
                    std::fs::create_dir_all(parent)?;
                }
                let mut fasta = String::new();
                for (i, &j) in selected.iter().enumerate() {
                    fasta.push_str(&format!(">prior{}\n{}\n", i + 1, table.sequences[j]));
                }
                std::fs::write(priors_path, fasta)?;
                if verbose {
                    eprintln!("[dada-pseudo] wrote priors to {}", priors_path.display());
                }
            }

            if verbose && (n_samples < 2 || prior_set.is_empty()) {
                eprintln!(
                    "[dada-pseudo] note: {} — round 2 is equivalent to round 1 (no priors applied)",
                    if n_samples < 2 {
                        "fewer than 2 samples".to_string()
                    } else {
                        "no sequences met the prior-selection threshold".to_string()
                    },
                );
            }

            // ---- Round 2: re-denoise each sample with the priors flagged ----
            if verbose {
                eprintln!("[dada-pseudo] round 2: denoising with priors");
            }
            for (s, sample_name) in sample_names.iter().enumerate() {
                let raws = &mut sample_raws[s];
                let n_marked = mark_priors(raws, &prior_set);
                if verbose {
                    eprintln!(
                        "[dada-pseudo]   round 2 {sample_name}: {n_marked} of {} unique(s) flagged as prior",
                        raws.len(),
                    );
                }
                let json = denoise_and_serialize(
                    "dada-pseudo",
                    sample_name,
                    raws,
                    &resolved.params,
                    &resolved.run,
                    &pool,
                    compact,
                    verbose,
                )?;
                let out_path = output_dir.join(format!("{sample_name}.json"));
                std::fs::write(&out_path, &json)?;
                if verbose {
                    eprintln!("[dada-pseudo] wrote {}", out_path.display());
                }
            }
        }

        Commands::MergePairs {
            fwd_dada,
            rev_dada,
            fwd_fastq,
            rev_fastq,
            min_overlap,
            max_mismatch,
            return_rejects,
            just_concatenate,
            concat_nnn_len,
            trim_overhang,
            sample_names,
            check_sample_ids,
            phred_offset,
            threads,
            output,
            compact,
            verbose,
        } => {
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
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::RemovePrimers {
            input,
            fout,
            sample_name,
            primer_fwd,
            primer_rev,
            rc_primer_rev,
            max_mismatch,
            allow_indels,
            trim_fwd,
            trim_rev,
            orient,
            compress,
            threads,
            trunc_q,
            trunc_len,
            trim_left,
            trim_right,
            max_len,
            min_len,
            max_n,
            min_q,
            max_ee,
            phix_genome,
            rm_lowcomplex,
            phred_offset,
            output,
            compact,
            verbose,
        } => {
            if allow_indels && verbose {
                eprintln!("[remove-primers] indel mode enabled — expect ~4× slower matching");
            }
            let primer_rev_bytes = primer_rev.map(|s| {
                let b = s.into_bytes();
                if rc_primer_rev {
                    iupac_reverse_complement(&b)
                } else {
                    b
                }
            });
            let filter_params = if trunc_q.is_some()
                || trunc_len.is_some()
                || trim_left.is_some()
                || trim_right.is_some()
                || max_len.is_some()
                || min_len.is_some()
                || max_n.is_some()
                || min_q.is_some()
                || max_ee.is_some()
                || phix_genome.is_some()
                || rm_lowcomplex.is_some()
            {
                let phix_seq: Option<Vec<u8>> = phix_genome
                    .as_deref()
                    .map(read_fasta_first_seq)
                    .transpose()?;
                Some(FilterParams {
                    trunc_q: trunc_q.unwrap_or(0),
                    trunc_len: trunc_len.unwrap_or(0),
                    trim_left: trim_left.unwrap_or(0),
                    trim_right: trim_right.unwrap_or(0),
                    max_len: max_len.unwrap_or(0),
                    min_len: min_len.unwrap_or(0),
                    max_n: max_n.unwrap_or(usize::MAX),
                    min_q: min_q.unwrap_or(0),
                    max_ee: max_ee.unwrap_or(f64::INFINITY),
                    phix_genome: phix_seq,
                    rm_lowcomplex: rm_lowcomplex.unwrap_or(0.0),
                    phred_offset,
                })
            } else {
                None
            };
            // Validate filter params before processing.
            if let Some(ref fp) = filter_params
                && fp.max_len > 0
                && fp.min_len > fp.max_len
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--min-len ({}) is greater than --max-len ({}); no read can satisfy both",
                        fp.min_len, fp.max_len
                    ),
                ));
            }
            let params = RemovePrimersParams {
                primer_fwd: primer_fwd.into_bytes(),
                primer_rev: primer_rev_bytes,
                max_mismatch,
                allow_indels,
                trim_fwd,
                trim_rev,
                orient,
                filter_params,
            };
            let sample = sample_name.unwrap_or_else(|| fastq_stem(&input));
            let stats = remove_primers(&input, &fout, &params, compress, threads, verbose)?;

            #[derive(Serialize)]
            struct RemovePrimersOutput {
                sample: String,
                reads_in: u64,
                reads_out: u64,
                reads_reoriented: u64,
                #[serde(skip_serializing_if = "is_zero")]
                reads_filter_fail: u64,
            }
            fn is_zero(v: &u64) -> bool {
                *v == 0
            }
            let tagged = Tagged::new(
                "remove-primers",
                RemovePrimersOutput {
                    sample,
                    reads_in: stats.reads_in,
                    reads_out: stats.reads_out,
                    reads_reoriented: stats.reads_reoriented,
                    reads_filter_fail: stats.reads_filter_fail,
                },
            );
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::FilterAndTrim {
            fwd,
            filt,
            rev,
            filt_rev,
            sample_name,
            compress,
            threads,
            trunc_q,
            trunc_len,
            trim_left,
            trim_right,
            max_len,
            min_len,
            max_n,
            min_q,
            max_ee,
            phix_genome,
            rm_lowcomplex,
            phred_offset,
            output,
            compact,
            verbose,
        } => {
            // ---- Validate paired-end files ----
            if rev.is_some() {
                filt_rev.as_ref().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--filt-rev is required when --rev is provided",
                    )
                })?;
            }

            // ---- Helper: expand a 1-or-2-element Vec into (fwd_val, rev_val) ----
            macro_rules! pair {
                ($v:expr, $default:expr) => {{
                    let v = &$v;
                    if v.is_empty() {
                        ($default, $default)
                    } else if v.len() == 1 {
                        (v[0], v[0])
                    } else {
                        (v[0], v[1])
                    }
                }};
            }

            let (trunc_q_f, trunc_q_r) = pair!(trunc_q, 2u8);
            let (trunc_len_f, trunc_len_r) = pair!(trunc_len, 0usize);
            let (trim_left_f, trim_left_r) = pair!(trim_left, 0usize);
            let (trim_right_f, trim_right_r) = pair!(trim_right, 0usize);
            let (max_len_f, max_len_r) = pair!(max_len, 0usize);
            let (min_len_f, min_len_r) = pair!(min_len, 20usize);
            let (max_ee_f, max_ee_r) = if max_ee.is_empty() {
                (f64::INFINITY, f64::INFINITY)
            } else if max_ee.len() == 1 {
                (max_ee[0], max_ee[0])
            } else {
                (max_ee[0], max_ee[1])
            };
            let (rm_lowcomplex_f, rm_lowcomplex_r) = pair!(rm_lowcomplex, 0.0f64);

            let phix_seq: Option<Vec<u8>> = phix_genome
                .as_deref()
                .map(read_fasta_first_seq)
                .transpose()?;

            let make_params = |tq, tl, trl, trr, ml, mnl, ee, rlc| FilterParams {
                trunc_q: tq,
                trunc_len: tl,
                trim_left: trl,
                trim_right: trr,
                max_len: ml,
                min_len: mnl,
                max_n,
                min_q,
                max_ee: ee,
                phix_genome: phix_seq.clone(),
                rm_lowcomplex: rlc,
                phred_offset,
            };

            let params_fwd = make_params(
                trunc_q_f,
                trunc_len_f,
                trim_left_f,
                trim_right_f,
                max_len_f,
                min_len_f,
                max_ee_f,
                rm_lowcomplex_f,
            );
            let params_rev = make_params(
                trunc_q_r,
                trunc_len_r,
                trim_left_r,
                trim_right_r,
                max_len_r,
                min_len_r,
                max_ee_r,
                rm_lowcomplex_r,
            );

            let sample = sample_name.unwrap_or_else(|| fastq_stem(&fwd));
            let opts = WriteOptions {
                compress,
                threads,
                verbose,
            };

            let stats = if let (Some(rev_in), Some(rev_out)) = (rev, filt_rev) {
                filter_paired(
                    &PairedFiles {
                        fwd_in: &fwd,
                        rev_in: &rev_in,
                        fwd_out: &filt,
                        rev_out: &rev_out,
                    },
                    &params_fwd,
                    &params_rev,
                    opts,
                )?
            } else {
                filter_single(&fwd, &filt, &params_fwd, opts)?
            };

            #[derive(Serialize)]
            struct FilterAndTrimOutput {
                sample: String,
                reads_in: u64,
                reads_out: u64,
            }
            let tagged = Tagged::new(
                "filter-and-trim",
                FilterAndTrimOutput {
                    sample,
                    reads_in: stats.reads_in,
                    reads_out: stats.reads_out,
                },
            );
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::MakeSequenceTable {
            input,
            sample_names,
            order_by,
            min_len,
            max_len,
            hash,
            output,
            compact,
        } => {
            if !sample_names.is_empty() && sample_names.len() != input.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--sample-names has {} entries but {} input file(s) were given",
                        sample_names.len(),
                        input.len()
                    ),
                ));
            }
            let order = match order_by.as_str() {
                "abundance" => OrderBy::Abundance,
                "nsamples" => OrderBy::NSamples,
                _ => OrderBy::None,
            };
            let names_opt = if sample_names.is_empty() {
                None
            } else {
                Some(sample_names.as_slice())
            };
            let hash_algo = if hash == "sha1" {
                HashAlgo::Sha1
            } else {
                HashAlgo::Md5
            };
            let paths: Vec<&Path> = input.iter().map(|p| p.as_path()).collect();
            let mut table = make_sequence_table(&paths, names_opt, order, hash_algo)?;

            if min_len.is_some() || max_len.is_some() {
                let keep: Vec<usize> = table
                    .sequences
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| {
                        min_len.is_none_or(|mn| s.len() >= mn)
                            && max_len.is_none_or(|mx| s.len() <= mx)
                    })
                    .map(|(j, _)| j)
                    .collect();
                table.sequences = keep.iter().map(|&j| table.sequences[j].clone()).collect();
                table.sequence_ids = keep
                    .iter()
                    .map(|&j| table.sequence_ids[j].clone())
                    .collect();
                for row in &mut table.counts {
                    *row = keep.iter().map(|&j| row[j]).collect();
                }
            }

            let tagged = Tagged::new("make-sequence-table", table);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;
            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::RemoveBimeraDenovo {
            input,
            method,
            min_fold_parent_over_abundance,
            min_parent_abundance,
            allow_one_off,
            min_one_off_parent_distance,
            max_shift,
            min_sample_fraction,
            ignore_n_negatives,
            threads,
            verbose,
            output,
            compact,
        } => {
            let table: SequenceTable =
                read_tagged_json(&input, &["make-sequence-table", "remove-bimera-denovo"])
                    .with_path(&input)?;

            let method = match method.as_str() {
                "pooled" => Method::Pooled,
                "per-sample" => Method::PerSample,
                _ => Method::Consensus,
            };
            let params = BimeraParams {
                min_fold_parent_over_abundance,
                min_parent_abundance,
                allow_one_off,
                min_one_off_parent_distance,
                max_shift,
                min_sample_fraction,
                ignore_n_negatives,
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
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
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::SeqTableToTsv {
            input,
            prevalence,
            min_abundance,
            output,
        } => {
            let table: SequenceTable =
                read_tagged_json(&input, &["make-sequence-table", "remove-bimera-denovo"])
                    .with_path(&input)?;
            let keep = select_sequences(&table, prevalence, min_abundance);

            let mut out: Box<dyn io::Write> = match output {
                Some(ref path) => Box::new(io::BufWriter::new(std::fs::File::create(path)?)),
                None => Box::new(io::BufWriter::new(std::io::stdout())),
            };

            // Header: sequence_id <TAB> sample1 <TAB> sample2 ...
            write!(out, "sequence_id")?;
            for sample in &table.samples {
                write!(out, "\t{sample}")?;
            }
            writeln!(out)?;

            // One row per kept sequence: id <TAB> count_per_sample...
            for &j in &keep {
                write!(out, "{}", table.sequence_ids[j])?;
                for sample_counts in &table.counts {
                    write!(out, "\t{}", sample_counts[j])?;
                }
                writeln!(out)?;
            }
            out.flush()?;
        }

        Commands::SeqTableToFasta {
            input,
            prevalence,
            min_abundance,
            output,
        } => {
            let table: SequenceTable =
                read_tagged_json(&input, &["make-sequence-table", "remove-bimera-denovo"])
                    .with_path(&input)?;

            if table.sequences.len() != table.sequence_ids.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sequence_ids and sequences lengths differ",
                ));
            }

            let keep = select_sequences(&table, prevalence, min_abundance);

            let mut out: Box<dyn io::Write> = match output {
                Some(ref path) => Box::new(io::BufWriter::new(std::fs::File::create(path)?)),
                None => Box::new(io::BufWriter::new(std::io::stdout())),
            };

            for &j in &keep {
                writeln!(out, ">{}\n{}", table.sequence_ids[j], table.sequences[j])?;
            }
            out.flush()?;
        }

        Commands::TaxToTsv {
            input,
            na_string,
            output,
        } => {
            #[derive(serde::Deserialize)]
            struct TaxAssignment {
                sequence_id: String,
                taxonomy: Vec<Option<String>>,
            }
            #[derive(serde::Deserialize)]
            struct TaxJson {
                levels: Vec<String>,
                assignments: Vec<TaxAssignment>,
            }

            let tax: TaxJson = read_tagged_json(&input, &["assign-taxonomy", "assign-species"])
                .with_path(&input)?;

            let mut out: Box<dyn io::Write> = match output {
                Some(ref path) => Box::new(io::BufWriter::new(std::fs::File::create(path)?)),
                None => Box::new(io::BufWriter::new(std::io::stdout())),
            };

            // Header: sequence_id <TAB> level1 <TAB> level2 ...
            write!(out, "sequence_id")?;
            for level in &tax.levels {
                write!(out, "\t{level}")?;
            }
            writeln!(out)?;

            for a in &tax.assignments {
                write!(out, "{}", a.sequence_id)?;
                for l in 0..tax.levels.len() {
                    let cell = a
                        .taxonomy
                        .get(l)
                        .and_then(|x| x.as_deref())
                        .unwrap_or(na_string.as_str());
                    write!(out, "\t{cell}")?;
                }
                writeln!(out)?;
            }
            out.flush()?;
        }

        Commands::Sample {
            input,
            output_dir,
            nbases,
            randomize,
            seed,
            phred_offset,
            threads,
            compact,
            verbose,
        } => {
            std::fs::create_dir_all(&output_dir)?;

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            // Optionally shuffle file order.
            let mut ordered: Vec<&std::path::PathBuf> = input.iter().collect();
            if randomize {
                use rand::SeedableRng as _;
                if let Some(s) = seed {
                    ordered.shuffle(&mut rand::rngs::SmallRng::seed_from_u64(s));
                } else {
                    ordered.shuffle(&mut rand::thread_rng());
                }
            }

            #[derive(Serialize)]
            struct UniqueEntry<'a> {
                sequence: &'a str,
                count: u64,
                mean_quality: &'a [f64],
            }
            #[derive(Serialize)]
            struct DerepOutput<'a> {
                sample: &'a str,
                total_reads: usize,
                unique_sequences: usize,
                sort_order: &'static str,
                uniques: Vec<UniqueEntry<'a>>,
            }
            #[derive(Serialize)]
            struct SampleSummary {
                samples_processed: usize,
                total_bases: u64,
                total_reads: u64,
                output_files: Vec<String>,
            }

            let mut total_bases: u64 = 0;
            let mut total_reads: u64 = 0;
            let mut output_files: Vec<String> = Vec::new();

            for path in &ordered {
                let is_gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
                let derep = if is_gz {
                    dereplicate(
                        MultiGzDecoder::new(File::open(path).with_path(path)?),
                        phred_offset,
                        &pool,
                        verbose,
                    )?
                } else {
                    dereplicate(
                        File::open(path).with_path(path)?,
                        phred_offset,
                        &pool,
                        verbose,
                    )?
                };

                let file_bases: u64 = derep
                    .uniques
                    .iter()
                    .map(|(seq, count)| seq.len() as u64 * count)
                    .sum();
                let file_reads: u64 = derep.map.len() as u64;

                // Build a stem for the output filename, stripping up to two extensions.
                let stem = {
                    let p = path.as_path();
                    let s1 = p.file_stem().unwrap_or_default();
                    let s1_path = std::path::Path::new(s1);
                    if s1_path.extension().is_some() {
                        s1_path
                            .file_stem()
                            .unwrap_or(s1)
                            .to_string_lossy()
                            .into_owned()
                    } else {
                        s1.to_string_lossy().into_owned()
                    }
                };
                let out_path = output_dir.join(format!("{stem}.json"));

                let mut uniq_entries = Vec::with_capacity(derep.uniques.len());
                for (i, (seq, count)) in derep.uniques.iter().enumerate() {
                    let sequence = std::str::from_utf8(seq)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                    uniq_entries.push(UniqueEntry {
                        sequence,
                        count: *count,
                        mean_quality: &derep.quals[i],
                    });
                }
                let sample_out = DerepOutput {
                    sample: &stem,
                    total_reads: derep.map.len(),
                    unique_sequences: uniq_entries.len(),
                    sort_order: "abundance_desc",
                    uniques: uniq_entries,
                };

                let unique_count = sample_out.unique_sequences;
                let tagged = Tagged::new("sample", sample_out);
                let json = if compact {
                    serde_json::to_string(&tagged)
                } else {
                    serde_json::to_string_pretty(&tagged)
                }
                .map_err(io::Error::other)?;

                std::fs::write(&out_path, &json)?;
                output_files.push(out_path.display().to_string());
                total_bases += file_bases;
                total_reads += file_reads;

                if verbose {
                    eprintln!(
                        "[sample] wrote {} ({} unique(s), {} bases)",
                        out_path.display(),
                        unique_count,
                        file_bases,
                    );
                }

                if total_bases >= nbases {
                    if verbose {
                        eprintln!(
                            "[sample] reached {} bases after {} file(s); stopping",
                            total_bases,
                            output_files.len(),
                        );
                    }
                    break;
                }
            }

            let summary = SampleSummary {
                samples_processed: output_files.len(),
                total_bases,
                total_reads,
                output_files,
            };
            let tagged = Tagged::new("sample", summary);
            let summary_json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;
            println!("{summary_json}");
        }

        Commands::ErrorsFromSample {
            input,
            errfun,
            pseudocount,
            binned_quals,
            errfun_cmd,
            loess_preset,
            loess_surface,
            loess_cell,
            loess_max_rate,
            loess_min_rate,
            max_consist,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            band,
            homo_gap_p,
            kdist_cutoff,
            kmer_size,
            no_kmer_screen,
            threads,
            output,
            compact,
            diag_dir,
            cluster_trace_dir,
            trace_no_members,
            trace_min_abund,
            verbose,
        } => {
            let loess_config = resolve_loess_config(
                &loess_preset,
                loess_surface.as_deref(),
                loess_cell,
                loess_max_rate,
                loess_min_rate,
            );
            let err_fun = match errfun.as_str() {
                "loess" => ErrFun::Loess {
                    config: loess_config,
                },
                "noqual" => ErrFun::Noqual {
                    pseudocount,
                    config: loess_config,
                },
                "binned-qual" => {
                    let bins = binned_quals.ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--binned-quals is required when --errfun binned-qual is used",
                        )
                    })?;
                    ErrFun::BinnedQual {
                        bins,
                        config: loess_config,
                    }
                }
                "pacbio" => ErrFun::PacBio {
                    config: loess_config,
                },
                "external" => {
                    let command = errfun_cmd.clone().ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--errfun-cmd is required when --errfun external is used",
                        )
                    })?;
                    if command.trim().is_empty() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--errfun-cmd cannot be empty",
                        ));
                    }
                    ErrFun::External { command }
                }
                other => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown errfun '{other}'; expected one of: loess, noqual, binned-qual, pacbio, external"
                        ),
                    ));
                }
            };

            let align_params = AlignParams {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                homo_gap_p,
                use_kmers: !no_kmer_screen,
                kdist_cutoff,
                kmer_size,
                band,
                vectorized: true,
                gapless: true,
            };

            let dada_params = dada::DadaParams {
                align: align_params,
                err_mat: Vec::new(),
                err_ncol: 0,
                omega_a,
                omega_c,
                omega_p,
                detect_singletons,
                max_clust: 0,
                min_fold,
                min_hamming,
                min_abund,
                use_quals: true,
                final_consensus: false,
                multithread: threads > 1,
                verbose,
                greedy: true,
                aux_outputs: false,
            };

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            let all_inputs = load_derep_samples(&input)?;

            if verbose {
                eprintln!(
                    "[errors-from-sample] loaded {} sample(s) from JSON",
                    all_inputs.len()
                );
            }

            if let Some(ref dir) = diag_dir {
                std::fs::create_dir_all(dir)?;
            }

            let params_snapshot =
                build_learned_err_params(&err_fun, max_consist, &dada_params, &align_params);

            if let Some(ref dir) = cluster_trace_dir {
                std::fs::create_dir_all(dir)?;
            }
            let trace_params = cluster_trace::TraceParams {
                no_members: trace_no_members,
                min_abund: trace_min_abund,
            };

            let result = pool.install(|| {
                learn_errors(
                    all_inputs,
                    &err_fun,
                    dada_params,
                    max_consist,
                    LearnDiagOptions {
                        verbose,
                        diag_dir: diag_dir.as_deref(),
                        cluster_trace_dir: cluster_trace_dir.as_deref(),
                        trace_params,
                    },
                )
            })?;

            #[derive(Serialize)]
            struct LearnErrorsOutput {
                nq: usize,
                converged: bool,
                stop_reason: learn_errors::StopReason,
                iterations: usize,
                params: LearnedErrParams,
                trans: Vec<Vec<u32>>,
                err_in: Vec<Vec<f64>>,
                err_out: Vec<Vec<f64>>,
            }

            fn flat_to_rows_u32(flat: &[u32], nq: usize) -> Vec<Vec<u32>> {
                (0..16)
                    .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
                    .collect()
            }
            fn flat_to_rows_f64(flat: &[f64], nq: usize) -> Vec<Vec<f64>> {
                (0..16)
                    .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
                    .collect()
            }

            let out = LearnErrorsOutput {
                nq: result.nq,
                converged: result.converged,
                stop_reason: result.stop_reason,
                iterations: result.iterations,
                params: params_snapshot,
                trans: flat_to_rows_u32(&result.trans, result.nq),
                err_in: flat_to_rows_f64(&result.err_in, result.nq),
                err_out: flat_to_rows_f64(&result.err_out, result.nq),
            };

            let tagged = Tagged::new("errors-from-sample", out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::AssignTaxonomy {
            input,
            ref_fasta,
            min_boot,
            try_rc,
            output_bootstraps,
            tax_levels,
            seed,
            threads,
            output,
            compact,
            verbose,
        } => {
            const MIN_REF_LEN: usize = 20;
            const DADA2_UNSPEC: &str = "_DADA2_UNSPECIFIED";

            // ---- Read queries ----
            let queries = read_query_sequences(&input).with_path(&input)?;
            let query_seqs: Vec<&[u8]> = queries.iter().map(|(_, s)| s.as_slice()).collect();
            let rcs: Vec<Vec<u8>> = if try_rc {
                query_seqs.iter().map(|&s| rc_bytes(s)).collect()
            } else {
                vec![]
            };
            let rc_refs: Vec<&[u8]> = rcs.iter().map(|s| s.as_slice()).collect();

            // ---- Read and parse reference FASTA ----
            let raw_refs = read_fasta_records(&ref_fasta).with_path(&ref_fasta)?;
            let raw_refs: Vec<(String, Vec<u8>)> = raw_refs
                .into_iter()
                .filter(|(_, seq)| seq.len() >= MIN_REF_LEN)
                .collect();

            if raw_refs.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "No reference sequences passed the minimum length filter.",
                ));
            }
            if !raw_refs[0].0.contains(';') {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Reference header does not look like a taxonomy string (no ';'). \
                     Use --ref-fasta with a DADA2-formatted taxonomy reference.",
                ));
            }

            // Parse each header into semicolon-delimited fields, finding max depth.
            let tax_fields: Vec<Vec<String>> = raw_refs
                .iter()
                .map(|(hdr, _)| {
                    hdr.split(';')
                        .filter(|s| !s.is_empty())
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>()
                })
                .collect();
            let max_depth = tax_fields.iter().map(|f| f.len()).max().unwrap_or(0);

            // Pad shorter strings with _DADA2_UNSPECIFIED.
            let tax_padded: Vec<Vec<String>> = tax_fields
                .into_iter()
                .map(|mut f| {
                    while f.len() < max_depth {
                        f.push(DADA2_UNSPEC.to_string());
                    }
                    f
                })
                .collect();

            // Build unique taxonomy strings and ref→genus mapping.
            let full_strings: Vec<String> = tax_padded.iter().map(|f| f.join(";")).collect();
            let mut genus_uniq: Vec<String> = {
                let mut seen = std::collections::HashSet::new();
                let mut v = Vec::new();
                for s in &full_strings {
                    if seen.insert(s.clone()) {
                        v.push(s.clone());
                    }
                }
                v
            };
            genus_uniq.sort_unstable(); // stable ordering for reproducibility
            let genus_to_idx: std::collections::HashMap<&str, usize> = genus_uniq
                .iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i))
                .collect();
            let ref_to_genus: Vec<usize> = full_strings
                .iter()
                .map(|s| genus_to_idx[s.as_str()])
                .collect();

            // Split each unique genus string into level fields.
            let genus_fields: Vec<Vec<String>> = genus_uniq
                .iter()
                .map(|s| {
                    let mut f: Vec<String> = s
                        .split(';')
                        .filter(|x| !x.is_empty())
                        .map(|x| x.to_string())
                        .collect();
                    while f.len() < max_depth {
                        f.push(DADA2_UNSPEC.to_string());
                    }
                    f
                })
                .collect();

            // Build integer-factor matrix [ngenus × nlevel] (1-based, sorted alpha).
            let ngenus = genus_uniq.len();
            let nlevel = max_depth;
            let mut genus_tax = vec![0usize; ngenus * nlevel];
            for l in 0..nlevel {
                let mut level_vals: Vec<String> =
                    genus_fields.iter().map(|f| f[l].clone()).collect();
                level_vals.sort_unstable();
                level_vals.dedup();
                let level_map: std::collections::HashMap<&str, usize> = level_vals
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.as_str(), i + 1))
                    .collect();
                for g in 0..ngenus {
                    genus_tax[g * nlevel + l] = level_map[genus_fields[g][l].as_str()];
                }
            }

            let ref_seqs: Vec<&[u8]> = raw_refs.iter().map(|(_, s)| s.as_slice()).collect();

            if verbose {
                eprintln!(
                    "[assign-taxonomy] {} queries, {} references, {} unique taxa, {} levels",
                    queries.len(),
                    ref_seqs.len(),
                    ngenus,
                    nlevel,
                );
            }

            // ---- Run classifier ----
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;

            let result = pool
                .install(|| {
                    assign_taxonomy(
                        &query_seqs,
                        &rc_refs,
                        &TaxonomyRef {
                            refs: &ref_seqs,
                            ref_to_genus: &ref_to_genus,
                            genus_tax: &genus_tax,
                            nlevel,
                        },
                        TaxonomyOptions {
                            try_rc,
                            seed,
                            verbose,
                        },
                    )
                })
                .map_err(io::Error::other)?;

            // ---- Assemble output ----
            #[derive(Serialize)]
            struct TaxAssignment {
                sequence_id: String,
                sequence: String,
                taxonomy: Vec<Option<String>>,
                #[serde(skip_serializing_if = "Option::is_none")]
                bootstrap: Option<Vec<u32>>,
            }
            #[derive(Serialize)]
            struct AssignTaxOutput {
                levels: Vec<String>,
                assignments: Vec<TaxAssignment>,
            }

            let out_levels: Vec<String> = tax_levels
                .iter()
                .take(nlevel)
                .cloned()
                .chain((tax_levels.len()..nlevel).map(|i| format!("Level{}", i + 1)))
                .collect();

            let assignments: Vec<TaxAssignment> = queries
                .iter()
                .enumerate()
                .map(|(i, (id, seq))| {
                    let (taxonomy, bootstrap) = if let Some(g) = result.assignments[i] {
                        let fields: Vec<&str> =
                            genus_fields[g].iter().map(|s| s.as_str()).collect();
                        let boot = &result.boot_counts[i];
                        let mut tax = Vec::with_capacity(nlevel);
                        let mut passed = true;
                        for l in 0..nlevel {
                            let b = boot.get(l).copied().unwrap_or(0);
                            if passed && b >= min_boot {
                                let s = fields.get(l).copied().unwrap_or(DADA2_UNSPEC);
                                tax.push(if s == DADA2_UNSPEC {
                                    None
                                } else {
                                    Some(s.to_string())
                                });
                            } else {
                                passed = false;
                                tax.push(None);
                            }
                        }
                        let boot_out = if output_bootstraps {
                            Some(boot.clone())
                        } else {
                            None
                        };
                        (tax, boot_out)
                    } else {
                        let boot_out = if output_bootstraps {
                            Some(vec![0u32; nlevel])
                        } else {
                            None
                        };
                        (vec![None; nlevel], boot_out)
                    };
                    TaxAssignment {
                        sequence_id: id.clone(),
                        sequence: String::from_utf8_lossy(seq).into_owned(),
                        taxonomy,
                        bootstrap,
                    }
                })
                .collect();

            let out = AssignTaxOutput {
                levels: out_levels,
                assignments,
            };
            let tagged = Tagged::new("assign-taxonomy", out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::AssignSpecies {
            input,
            ref_fasta,
            allow_multiple,
            try_rc,
            output,
            compact,
            verbose,
        } => {
            const MIN_REF_LEN: usize = 20;

            // ---- Read input taxonomy JSON ----
            #[derive(serde::Deserialize)]
            struct TaxAssignmentIn {
                sequence_id: String,
                sequence: String,
                taxonomy: Vec<Option<String>>,
                #[serde(default)]
                bootstrap: Option<Vec<u32>>,
            }
            #[derive(serde::Deserialize)]
            struct AssignTaxIn {
                levels: Vec<String>,
                assignments: Vec<TaxAssignmentIn>,
            }
            let tax_in: AssignTaxIn =
                read_tagged_json(&input, &["assign-taxonomy"]).with_path(&input)?;

            let query_seqs_owned: Vec<Vec<u8>> = tax_in
                .assignments
                .iter()
                .map(|a| a.sequence.as_bytes().to_vec())
                .collect();
            let query_seqs: Vec<&[u8]> = query_seqs_owned.iter().map(|s| s.as_slice()).collect();

            // ---- Read and parse species reference FASTA ----
            let raw_refs = read_fasta_records(&ref_fasta).with_path(&ref_fasta)?;
            let mut ref_seqs_owned: Vec<Vec<u8>> = Vec::new();
            let mut ref_genus_owned: Vec<String> = Vec::new();
            let mut ref_species_owned: Vec<String> = Vec::new();

            for (header, seq) in raw_refs {
                if seq.len() < MIN_REF_LEN {
                    continue;
                }
                let mut fields = header.split_whitespace();
                let _id = fields.next(); // accession, ignored
                let genus = fields.next().unwrap_or("").to_string();
                let species = fields.next().unwrap_or("").to_string();
                if genus.is_empty() || species.is_empty() {
                    continue;
                }
                ref_seqs_owned.push(seq);
                ref_genus_owned.push(genus);
                ref_species_owned.push(species);
            }

            if ref_seqs_owned.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "No valid reference sequences found. Expected '>ID genus species' headers.",
                ));
            }

            if verbose {
                eprintln!(
                    "[assign-species] {} queries, {} references",
                    query_seqs.len(),
                    ref_seqs_owned.len(),
                );
            }

            let ref_seqs: Vec<&[u8]> = ref_seqs_owned.iter().map(|s| s.as_slice()).collect();
            let ref_genus: Vec<&str> = ref_genus_owned.iter().map(|s| s.as_str()).collect();
            let ref_species: Vec<&str> = ref_species_owned.iter().map(|s| s.as_str()).collect();

            let hits: Vec<SpeciesHit> = assign_species(
                &query_seqs,
                &SpeciesRef {
                    ref_seqs: &ref_seqs,
                    ref_genus: &ref_genus,
                    ref_species: &ref_species,
                },
                SpeciesOptions {
                    max_species: allow_multiple,
                    try_rc,
                    verbose,
                },
            );

            // ---- Build new levels: drop existing "Species", append new "Species" ----
            let genus_idx = tax_in.levels.iter().position(|l| l == "Genus");
            let species_idx = tax_in.levels.iter().position(|l| l == "Species");
            let new_levels: Vec<String> = tax_in
                .levels
                .iter()
                .filter(|l| *l != "Species")
                .cloned()
                .chain(std::iter::once("Species".to_string()))
                .collect();

            // ---- Combine taxonomy + species hit per assignment ----
            #[derive(Serialize)]
            struct TaxAssignmentOut {
                sequence_id: String,
                sequence: String,
                taxonomy: Vec<Option<String>>,
                #[serde(skip_serializing_if = "Option::is_none")]
                bootstrap: Option<Vec<u32>>,
            }

            let assignments_out: Vec<TaxAssignmentOut> = tax_in
                .assignments
                .into_iter()
                .zip(hits)
                .map(|(a, hit)| {
                    // Genus matching: when input has a Genus level, only fill species
                    // if the species hit's genus matches the assigned genus
                    // (mirrors R's matchGenera rules).
                    let species = match (genus_idx, hit.genus.as_deref()) {
                        (Some(gi), Some(hit_gen)) => {
                            match a.taxonomy.get(gi).and_then(|x| x.as_deref()) {
                                Some(assigned_gen) if match_genera(assigned_gen, hit_gen) => {
                                    hit.species
                                }
                                _ => None,
                            }
                        }
                        _ => hit.species,
                    };

                    // Drop the old Species column from taxonomy + bootstrap, then
                    // append the species call (no bootstrap entry — exact-match).
                    let mut new_tax: Vec<Option<String>> = a
                        .taxonomy
                        .into_iter()
                        .enumerate()
                        .filter(|(i, _)| Some(*i) != species_idx)
                        .map(|(_, t)| t)
                        .collect();
                    new_tax.push(species);

                    let new_boot = a.bootstrap.map(|mut b| {
                        if let Some(si) = species_idx
                            && si < b.len()
                        {
                            b.remove(si);
                        }
                        b
                    });

                    TaxAssignmentOut {
                        sequence_id: a.sequence_id,
                        sequence: a.sequence,
                        taxonomy: new_tax,
                        bootstrap: new_boot,
                    }
                })
                .collect();

            #[derive(Serialize)]
            struct AssignSpeciesOutput {
                levels: Vec<String>,
                assignments: Vec<TaxAssignmentOut>,
            }

            let out = AssignSpeciesOutput {
                levels: new_levels,
                assignments: assignments_out,
            };
            let tagged = Tagged::new("assign-species", out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::LearnErrors {
            input,
            nbases,
            randomize,
            seed,
            phred_offset,
            errfun,
            pseudocount,
            binned_quals,
            errfun_cmd,
            loess_preset,
            loess_surface,
            loess_cell,
            loess_max_rate,
            loess_min_rate,
            max_consist,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            band,
            homo_gap_p,
            kdist_cutoff,
            kmer_size,
            no_kmer_screen,
            threads,
            output,
            compact,
            diag_dir,
            cluster_trace_dir,
            trace_no_members,
            trace_min_abund,
            verbose,
        } => {
            let loess_config = resolve_loess_config(
                &loess_preset,
                loess_surface.as_deref(),
                loess_cell,
                loess_max_rate,
                loess_min_rate,
            );
            let err_fun = match errfun.as_str() {
                "loess" => ErrFun::Loess {
                    config: loess_config,
                },
                "noqual" => ErrFun::Noqual {
                    pseudocount,
                    config: loess_config,
                },
                "binned-qual" => {
                    let bins = binned_quals.ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--binned-quals is required when --errfun binned-qual is used",
                        )
                    })?;
                    ErrFun::BinnedQual {
                        bins,
                        config: loess_config,
                    }
                }
                "pacbio" => ErrFun::PacBio {
                    config: loess_config,
                },
                "external" => {
                    let command = errfun_cmd.clone().ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--errfun-cmd is required when --errfun external is used",
                        )
                    })?;
                    if command.trim().is_empty() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--errfun-cmd cannot be empty",
                        ));
                    }
                    ErrFun::External { command }
                }
                other => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown errfun '{other}'; expected one of: loess, noqual, binned-qual, pacbio, external"
                        ),
                    ));
                }
            };

            let align_params = AlignParams {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                homo_gap_p,
                use_kmers: !no_kmer_screen,
                kdist_cutoff,
                kmer_size,
                band,
                vectorized: true,
                gapless: true,
            };

            let dada_params = dada::DadaParams {
                align: align_params,
                err_mat: Vec::new(), // overwritten each iteration
                err_ncol: 0,         // overwritten each iteration
                omega_a,
                omega_c,
                omega_p,
                detect_singletons,
                max_clust: 0,
                min_fold,
                min_hamming,
                min_abund,
                use_quals: true,
                final_consensus: false,
                multithread: threads > 1,
                verbose,
                greedy: true,
                aux_outputs: false,
            };

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)?;
            let all_inputs = load_fastq_samples(
                &input,
                nbases,
                randomize,
                seed,
                phred_offset,
                &pool,
                verbose,
            )?;

            if let Some(ref dir) = diag_dir {
                std::fs::create_dir_all(dir)?;
            }

            let params_snapshot =
                build_learned_err_params(&err_fun, max_consist, &dada_params, &align_params);

            if let Some(ref dir) = cluster_trace_dir {
                std::fs::create_dir_all(dir)?;
            }
            let trace_params = cluster_trace::TraceParams {
                no_members: trace_no_members,
                min_abund: trace_min_abund,
            };

            let result = pool.install(|| {
                learn_errors(
                    all_inputs,
                    &err_fun,
                    dada_params,
                    max_consist,
                    LearnDiagOptions {
                        verbose,
                        diag_dir: diag_dir.as_deref(),
                        cluster_trace_dir: cluster_trace_dir.as_deref(),
                        trace_params,
                    },
                )
            })?;

            // Serialize: represent the three matrices as Vec<Vec<T>> (16 rows × nq cols).
            #[derive(Serialize)]
            struct LearnErrorsOutput {
                nq: usize,
                converged: bool,
                stop_reason: learn_errors::StopReason,
                iterations: usize,
                /// Provenance: parameters used for the dada_uniques runs that
                /// produced this err model. Embedded so a downstream `dada`
                /// invocation can validate or inherit them. See
                /// `LearnedErrParams` for field details.
                params: LearnedErrParams,
                /// Transition counts: 16 rows (ref_nt*4+query_nt), nq columns.
                trans: Vec<Vec<u32>>,
                /// Error rates fed into the final DADA run: 16 × nq.
                err_in: Vec<Vec<f64>>,
                /// Error rates estimated from `trans`: 16 × nq.
                err_out: Vec<Vec<f64>>,
            }

            fn flat_to_rows_u32(flat: &[u32], nq: usize) -> Vec<Vec<u32>> {
                (0..16)
                    .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
                    .collect()
            }
            fn flat_to_rows_f64(flat: &[f64], nq: usize) -> Vec<Vec<f64>> {
                (0..16)
                    .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
                    .collect()
            }

            let out = LearnErrorsOutput {
                nq: result.nq,
                converged: result.converged,
                stop_reason: result.stop_reason,
                iterations: result.iterations,
                params: params_snapshot,
                trans: flat_to_rows_u32(&result.trans, result.nq),
                err_in: flat_to_rows_f64(&result.err_in, result.nq),
                err_out: flat_to_rows_f64(&result.err_out, result.nq),
            };

            let tagged = Tagged::new("learn-errors", out);
            let json = if compact {
                serde_json::to_string(&tagged)
            } else {
                serde_json::to_string_pretty(&tagged)
            }
            .map_err(io::Error::other)?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }
    }

    Ok(())
}

/// Collect all FASTQ files (`.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`) from a directory.
/// Returns paths in arbitrary order; the caller is responsible for sorting or shuffling.
/// Read query sequences from a FASTA file or a sequence-table JSON.
///
/// Returns `(sequence_id, sequence_bytes)` pairs.  The input format is
/// detected by file extension: `.json` (or `.json.gz`) triggers JSON parse;
/// anything else is treated as FASTA.
fn read_query_sequences(path: &Path) -> io::Result<Vec<(String, Vec<u8>)>> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    let is_json = ext == "json"
        || path
            .file_stem()
            .and_then(|s| Path::new(s).extension())
            .and_then(|e| e.to_str())
            == Some("json");

    if is_json {
        #[derive(serde::Deserialize)]
        struct SeqTable {
            sequences: Vec<String>,
            sequence_ids: Vec<String>,
        }
        let table: SeqTable =
            read_tagged_json(path, &["make-sequence-table", "remove-bimera-denovo"])
                .with_path(path)?;
        if table.sequences.len() != table.sequence_ids.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sequence_ids and sequences differ in length",
            ));
        }
        Ok(table
            .sequence_ids
            .into_iter()
            .zip(table.sequences)
            .map(|(id, seq)| (id, seq.into_bytes()))
            .collect())
    } else {
        read_fasta_records(path).with_path(path)
    }
}

fn rc_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            _ => b'N',
        })
        .collect()
}

/// Derive a base stem from a FASTQ path by stripping recognised extensions.
///
/// `sample1.fastq.gz` → `"sample1"`,  `sample2.fq` → `"sample2"`.
/// Mirrors R DADA2's `matchGenera()`.  Returns `true` when the genus assigned
/// by the taxonomy classifier (`gen_tax`) matches the genus of a species
/// reference hit (`gen_binom`).  Tolerates split genus names like
/// `Escherichia/Shigella` and the "Candidatus X" prefix form.
fn match_genera(gen_tax: &str, gen_binom: &str) -> bool {
    if gen_tax == gen_binom {
        return true;
    }
    // gen_tax starts with "<gen_binom> " — e.g. "Candidatus Saccharimonas" vs "Candidatus".
    if gen_tax.len() > gen_binom.len()
        && gen_tax.starts_with(gen_binom)
        && gen_tax.as_bytes()[gen_binom.len()] == b' '
    {
        return true;
    }
    // gen_binom is a "/"-split genus that contains gen_tax at either end.
    if gen_binom.starts_with(gen_tax)
        && gen_binom.len() > gen_tax.len()
        && gen_binom.as_bytes()[gen_tax.len()] == b'/'
    {
        return true;
    }
    if gen_binom.ends_with(gen_tax)
        && gen_binom.len() > gen_tax.len()
        && gen_binom.as_bytes()[gen_binom.len() - gen_tax.len() - 1] == b'/'
    {
        return true;
    }
    false
}

/// `true` when `path` looks like a JSON file (`.json` or `.json.gz`).
fn is_json_path(path: &Path) -> bool {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    ext == "json"
        || path
            .file_stem()
            .and_then(|s| Path::new(s).extension())
            .and_then(|e| e.to_str())
            == Some("json")
}

/// Error-model JSON shape shared by `dada`, `dada-pooled`, and `dada-pseudo`.
#[derive(serde::Deserialize)]
struct ErrorModelJson {
    nq: usize,
    err_in: Vec<Vec<f64>>,
    err_out: Vec<Vec<f64>>,
    #[serde(default)]
    params: Option<LearnedErrParams>,
}

/// Bundle of resolved DADA parameters plus the scalars needed to serialize a
/// `DadaRunParams` provenance block. Produced by [`resolve_dada_params`] so the
/// `dada`, `dada-pooled`, and `dada-pseudo` handlers all share one code path.
struct ResolvedDada {
    params: dada::DadaParams,
    nq: usize,
    run: DadaRunParams,
}

/// Load the error model and resolve every DADA parameter via the three-tier
/// precedence (CLI explicit > inherited from err-model `params` > built-in
/// default). Emits the same warnings the inline handlers used to.
///
/// `aux_outputs` and `pool` are handler-specific and passed in. `n_prior` is
/// filled in later by the caller (priors are marked after this point), so it is
/// left at 0 here.
#[allow(clippy::too_many_arguments)]
fn resolve_dada_params(
    error_model: &Path,
    use_err_in: bool,
    inherit_err_params: bool,
    threads: usize,
    aux_outputs: bool,
    pool: bool,
    verbose: bool,
    omega_a: Option<f64>,
    omega_c: Option<f64>,
    omega_p: Option<f64>,
    min_fold: Option<f64>,
    min_hamming: Option<u32>,
    min_abund: Option<u32>,
    detect_singletons: Option<bool>,
    band: Option<i32>,
    homo_gap_p: Option<i32>,
    kdist_cutoff: Option<f64>,
    kmer_size: Option<usize>,
    no_kmer_screen: Option<bool>,
) -> io::Result<ResolvedDada> {
    let em: ErrorModelJson = read_tagged_json(error_model, &["learn-errors", "errors-from-sample"])
        .with_path(error_model)?;
    let nq = em.nq;
    let rows = if use_err_in { &em.err_in } else { &em.err_out };
    if rows.len() != 16 || rows.iter().any(|r| r.len() != nq) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Error model matrix must be 16 × {nq}, got {} rows",
                rows.len()
            ),
        ));
    }
    let err_mat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();

    if inherit_err_params && em.params.is_none() {
        eprintln!(
            "[dada] warning: --inherit-err-params requested but error model {} has no `params` block (likely produced by a pre-provenance version); falling back to built-in defaults",
            error_model.display(),
        );
    }

    let p = em.params.as_ref();
    macro_rules! resolve {
        ($cli:expr, $em_field:ident, $default:expr) => {{
            match ($cli, inherit_err_params, p) {
                (Some(v), _, _) => v,
                (None, true, Some(em_params)) => em_params.$em_field,
                _ => $default,
            }
        }};
    }
    let omega_a = resolve!(omega_a, omega_a, 1e-40);
    // `omega_c` is intentionally not inherited from the err model:
    // learn-errors uses 0 (R DADA2 convention), dada uses 1e-40.
    let omega_c = omega_c.unwrap_or(1e-40);
    let omega_p = resolve!(omega_p, omega_p, 1e-4);
    let min_fold = resolve!(min_fold, min_fold, 1.0);
    let min_hamming = resolve!(min_hamming, min_hamming, 1);
    let min_abund = resolve!(min_abund, min_abund, 1);
    let detect_singletons = resolve!(detect_singletons, detect_singletons, false);
    let band = resolve!(band, band, 16);
    let homo_gap_p = resolve!(homo_gap_p, homo_gap_p, -8);
    let kdist_cutoff = resolve!(kdist_cutoff, kdist_cutoff, 0.42);
    let kmer_size = resolve!(kmer_size, kmer_size, 5);
    let use_kmers = match (no_kmer_screen, inherit_err_params, p) {
        (Some(no), _, _) => !no,
        (None, true, Some(em_params)) => em_params.use_kmers,
        _ => true,
    };

    // ---- Consistency warnings (only when NOT inheriting) ----
    if !inherit_err_params && let Some(em_params) = p {
        let mut mismatches: Vec<String> = Vec::new();
        macro_rules! check {
            ($name:literal, $cli_val:expr, $em_val:expr) => {
                if $cli_val != $em_val {
                    mismatches.push(format!(
                        "  {} = {:?} (err model: {:?})",
                        $name, $cli_val, $em_val
                    ));
                }
            };
        }
        check!("omega_a", omega_a, em_params.omega_a);
        check!("omega_p", omega_p, em_params.omega_p);
        check!("min_fold", min_fold, em_params.min_fold);
        check!("min_hamming", min_hamming, em_params.min_hamming);
        check!("min_abund", min_abund, em_params.min_abund);
        check!(
            "detect_singletons",
            detect_singletons,
            em_params.detect_singletons
        );
        check!("band", band, em_params.band);
        check!("homo_gap_p", homo_gap_p, em_params.homo_gap_p);
        check!("kdist_cutoff", kdist_cutoff, em_params.kdist_cutoff);
        check!("kmer_size", kmer_size, em_params.kmer_size);
        check!("use_kmers", use_kmers, em_params.use_kmers);
        if !mismatches.is_empty() {
            eprintln!(
                "[dada] warning: {} dada parameter(s) differ from error model {}; pass --inherit-err-params to adopt the err model's values:",
                mismatches.len(),
                error_model.display(),
            );
            for line in &mismatches {
                eprintln!("{line}");
            }
        }
    }

    let align_params = AlignParams {
        match_score: 5,
        mismatch: -4,
        gap_p: -8,
        homo_gap_p,
        use_kmers,
        kdist_cutoff,
        kmer_size,
        band,
        vectorized: true,
        gapless: true,
    };

    let params = dada::DadaParams {
        align: align_params,
        err_mat,
        err_ncol: nq,
        omega_a,
        omega_c,
        omega_p,
        detect_singletons,
        max_clust: 0,
        min_fold,
        min_hamming,
        min_abund,
        use_quals: true,
        final_consensus: false,
        multithread: threads > 1,
        verbose,
        greedy: true,
        aux_outputs,
    };

    let run = DadaRunParams {
        omega_a,
        omega_c,
        omega_p,
        min_fold,
        min_hamming,
        min_abund,
        detect_singletons,
        band,
        homo_gap_p,
        kdist_cutoff,
        kmer_size,
        use_kmers,
        pool,
        n_prior: 0,
    };

    Ok(ResolvedDada { params, nq, run })
}

/// Mark each RawInput whose (uppercased) sequence is in `prior_set` as a prior.
/// Returns the number marked.
fn mark_priors(
    raws: &mut [dada::RawInput],
    prior_set: &std::collections::HashSet<String>,
) -> usize {
    let mut n = 0usize;
    for inp in raws.iter_mut() {
        if prior_set.contains(&inp.seq.to_ascii_uppercase()) {
            inp.prior = true;
            n += 1;
        }
    }
    n
}

/// Run `dada_uniques` on one sample's uniques and serialize the standard
/// (non-pooled) `DadaOutput` JSON, returning the pretty/compact string. Used by
/// the multi-input `dada` path and by `dada-pseudo`. The single-input `dada`
/// path keeps its own inline serialization because it also emits aux outputs
/// and cluster traces.
#[allow(clippy::too_many_arguments)]
fn denoise_and_serialize(
    tag: &'static str,
    sample: &str,
    raw_inputs: &[dada::RawInput],
    params: &dada::DadaParams,
    run_params: &DadaRunParams,
    pool: &rayon::ThreadPool,
    compact: bool,
    verbose: bool,
) -> io::Result<String> {
    #[derive(Serialize)]
    struct AsvEntry {
        sequence: String,
        abundance: u32,
        birth_type: String,
        birth_pval: f64,
        birth_fold: f64,
        birth_e: f64,
    }
    #[derive(Serialize)]
    struct DadaStats {
        nalign: u32,
        nshroud: u32,
    }
    #[derive(Serialize)]
    struct DadaOutput {
        sample: String,
        num_asvs: usize,
        total_reads: u32,
        asvs: Vec<AsvEntry>,
        stats: DadaStats,
        params: DadaRunParams,
        map: Vec<Option<usize>>,
    }

    let result = pool
        .install(|| dada::dada_uniques(raw_inputs, params))
        .map_err(io::Error::other)?;

    if verbose {
        eprintln!(
            "[{tag}] {} ASV(s) from {} unique input(s); {} aligns, {} shrouded",
            result.clusters.len(),
            raw_inputs.len(),
            result.nalign,
            result.nshroud,
        );
    }

    let total_reads: u32 = result.clusters.iter().map(|c| c.reads).sum();
    let asvs: Vec<AsvEntry> = result
        .clusters
        .iter()
        .map(|c| {
            let sequence: String = c
                .sequence
                .iter()
                .map(|&b| misc::nt_decode(b) as char)
                .collect();
            let birth_type = match &c.birth_type {
                BirthType::Initial => "Initial",
                BirthType::Abundance => "Abundance",
                BirthType::Prior => "Prior",
                BirthType::Singleton => "Singleton",
            }
            .to_string();
            AsvEntry {
                sequence,
                abundance: c.reads,
                birth_type,
                birth_pval: c.birth_pval,
                birth_fold: c.birth_fold,
                birth_e: c.birth_e,
            }
        })
        .collect();

    let mut run_params = *run_params;
    run_params.n_prior = raw_inputs.iter().filter(|r| r.prior).count();

    let out = DadaOutput {
        sample: sample.to_string(),
        num_asvs: asvs.len(),
        total_reads,
        asvs,
        stats: DadaStats {
            nalign: result.nalign,
            nshroud: result.nshroud,
        },
        params: run_params,
        map: result.map,
    };

    let tagged = Tagged::new(tag, out);
    if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)
}

/// Build a [`derep::Derep`] for `dada` / `dada-pooled` from either a FASTQ file
/// (uncompressed or gzipped) or a derep/sample JSON file.
///
/// JSON inputs are defensively sorted by abundance descending — DADA2 assumes
/// the most-abundant unique is at index 0.  The `map` (read → unique) field is
/// only populated from the FASTQ path; JSON inputs leave it empty since neither
/// `dada` nor `dada-pooled` consult it.
///
/// Returns the dereplicated table plus the JSON's embedded `sample` field
/// when present; FASTQ inputs always return `None` for the name.
fn load_derep_for_dada(
    path: &Path,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    verbose: bool,
) -> io::Result<(derep::Derep, Option<String>)> {
    if is_json_path(path) {
        #[derive(serde::Deserialize)]
        struct UniqueEntryJson {
            sequence: String,
            count: u64,
            mean_quality: Vec<f64>,
        }
        #[derive(serde::Deserialize)]
        struct SampleJson {
            #[serde(default)]
            sample: Option<String>,
            #[serde(default)]
            sort_order: Option<String>,
            uniques: Vec<UniqueEntryJson>,
        }
        let parsed: SampleJson = read_tagged_json(path, &["derep", "sample"]).with_path(path)?;
        let sample_name = parsed.sample;
        let mut entries = parsed.uniques;
        // Skip the defensive sort when the producer has declared the order.
        // Older JSONs without `sort_order` get sorted, matching prior behaviour.
        if parsed.sort_order.as_deref() != Some("abundance_desc") {
            entries.sort_by_key(|a| std::cmp::Reverse(a.count));
        }
        let mut uniques: Vec<(Vec<u8>, u64)> = Vec::with_capacity(entries.len());
        let mut quals: Vec<Vec<f64>> = Vec::with_capacity(entries.len());
        for u in entries {
            if !u.mean_quality.is_empty() && u.mean_quality.len() != u.sequence.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}: mean_quality length {} != sequence length {}",
                        path.display(),
                        u.mean_quality.len(),
                        u.sequence.len(),
                    ),
                ));
            }
            uniques.push((u.sequence.into_bytes(), u.count));
            quals.push(u.mean_quality);
        }
        Ok((
            derep::Derep {
                uniques,
                quals,
                map: Vec::new(),
            },
            sample_name,
        ))
    } else if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        let derep = dereplicate(
            MultiGzDecoder::new(File::open(path).with_path(path)?),
            phred_offset,
            pool,
            verbose,
        )?;
        Ok((derep, None))
    } else {
        let derep = dereplicate(
            File::open(path).with_path(path)?,
            phred_offset,
            pool,
            verbose,
        )?;
        Ok((derep, None))
    }
}

/// Sequence-table column filter mirroring R DADA2's pseudo-pooling prior selection:
///   keep[j] = (n_samples_present[j] >= prevalence) || (total_abundance[j] >= min_abundance)
/// When both thresholds are `None` every column is kept.
fn select_sequences(
    table: &SequenceTable,
    prevalence: Option<u32>,
    min_abundance: Option<u64>,
) -> Vec<usize> {
    let nseq = table.sequences.len();
    if prevalence.is_none() && min_abundance.is_none() {
        return (0..nseq).collect();
    }
    (0..nseq)
        .filter(|&j| {
            let by_prev = prevalence.is_some_and(|p| {
                let n_present = table.counts.iter().filter(|row| row[j] > 0).count() as u32;
                n_present >= p
            });
            let by_abund = min_abundance.is_some_and(|m| {
                let total: u64 = table.counts.iter().map(|row| row[j]).sum();
                total >= m
            });
            by_prev || by_abund
        })
        .collect()
}

fn fastq_stem(path: &Path) -> String {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown");

    for suffix in &[".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if let Some(stem) = name.strip_suffix(suffix) {
            return stem.to_string();
        }
    }
    // Fallback: use whatever Path::file_stem gives.
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn match_genera_rules() {
        // Exact match.
        assert!(match_genera("Lactobacillus", "Lactobacillus"));

        // Split-genus reference: "/"-joined name on either side of gen_tax.
        assert!(match_genera("Escherichia", "Escherichia/Shigella"));
        assert!(match_genera("Shigella", "Escherichia/Shigella"));
        assert!(!match_genera("Salmonella", "Escherichia/Shigella"));

        // "Candidatus X" form: gen_tax has the prefix word that gen_binom is.
        assert!(match_genera("Candidatus Saccharimonas", "Candidatus"));

        // Mismatches.
        assert!(!match_genera("Lactobacillus", "Streptococcus"));
        // Substring without separator must not match.
        assert!(!match_genera("Lacto", "Lactobacillus"));
        assert!(!match_genera("Lactobacillus", "Lacto"));
    }
}
