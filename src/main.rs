use std::{fs::File, io, path::Path};

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rand::seq::SliceRandom as _;

mod chimera;
mod remove_bimera;
mod sequence_table;
mod cli;
mod merge_pairs;
mod cluster;
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
mod misc;
mod nwalign;
mod pval;
mod summary;
mod taxonomy;

use cli::{Cli, Commands};
use containers::BirthType;
use derep::dereplicate;
use filter_trim::{FilterParams, filter_single, filter_paired, read_fasta_first_seq};
use remove_bimera::{BimeraParams, Method, remove_bimera_denovo};
use sequence_table::{HashAlgo, OrderBy, SequenceTable, make_sequence_table};
use learn_errors::{ErrFun, learn_errors, load_derep_samples, load_fastq_samples};
use nwalign::AlignParams;
use serde::Serialize;
use summary::process;

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Summary { input, phred_offset, threads, output, compact } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            let summary = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                process(MultiGzDecoder::new(File::open(&input)?), phred_offset, &pool)?
            } else {
                process(File::open(&input)?, phred_offset, &pool)?
            };

            #[derive(Serialize)]
            struct SummaryOutput {
                total_reads: u64,
                mean_quality_per_position: Vec<f64>,
            }

            let out = SummaryOutput {
                total_reads: summary.total_reads,
                mean_quality_per_position: summary.mean_quality_per_position(),
            };

            let json = if compact {
                serde_json::to_string(&out)
            } else {
                serde_json::to_string_pretty(&out)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::Derep { input, phred_offset, threads, output, show_map, compact, verbose } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            let derep = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                dereplicate(MultiGzDecoder::new(File::open(&input)?), phred_offset, &pool, verbose)?
            } else {
                dereplicate(File::open(&input)?, phred_offset, &pool, verbose)?
            };

            #[derive(Serialize)]
            struct UniqueEntry<'a> {
                sequence: &'a str,
                count: u64,
                mean_quality: &'a [f64],
            }

            #[derive(Serialize)]
            struct DerepOutput<'a> {
                total_reads: usize,
                unique_sequences: usize,
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

            let derep_out = DerepOutput {
                total_reads: derep.map.len(),
                unique_sequences: derep.uniques.len(),
                uniques: uniq_entries,
                map: if show_map { Some(&derep.map) } else { None },
            };

            let json = if compact {
                serde_json::to_string(&derep_out)
            } else {
                serde_json::to_string_pretty(&derep_out)
            }.map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::Dada {
            input,
            error_model,
            use_err_in,
            phred_offset,
            threads,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            show_map,
            output,
            compact,
            verbose,
        } => {
            // ---- Dereplicate FASTQ in memory ----
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            let derep = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                dereplicate(MultiGzDecoder::new(File::open(&input)?), phred_offset, &pool, verbose)?
            } else {
                dereplicate(File::open(&input)?, phred_offset, &pool, verbose)?
            };

            let raw_inputs: Vec<dada::RawInput> = derep
                .uniques
                .into_iter()
                .zip(derep.quals.into_iter())
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
                    "FASTQ contains no reads",
                ));
            }

            // ---- Load error model JSON ----
            // The learn-errors output stores err_in and err_out as Vec<Vec<f64>>
            // with 16 rows and nq columns.
            #[derive(serde::Deserialize)]
            struct ErrorModelJson {
                nq: usize,
                err_in: Vec<Vec<f64>>,
                err_out: Vec<Vec<f64>>,
            }

            let em: ErrorModelJson = misc::read_json_file(&error_model)?;

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
            // Flatten row-major: row r, column q → index r*nq + q
            let err_mat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();

            // ---- Build algorithm parameters ----
            let align_params = AlignParams {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                homo_gap_p: -8,
                use_kmers: true,
                kdist_cutoff: 0.42,
                band: -1,
                vectorized: false,
                gapless: false,
            };

            let dada_params = dada::DadaParams {
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
            };

            // ---- Run DADA2 ----
            let result = dada::dada_uniques(&raw_inputs, &dada_params)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            if verbose {
                eprintln!(
                    "[dada] {} ASV(s) from {} unique input(s); {} aligns, {} shrouded",
                    result.clusters.len(),
                    raw_inputs.len(),
                    result.nalign,
                    result.nshroud,
                );
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
            struct DadaOutput {
                num_asvs: usize,
                total_reads: u32,
                asvs: Vec<AsvEntry>,
                stats: DadaStats,
                #[serde(skip_serializing_if = "Option::is_none")]
                map: Option<Vec<Option<usize>>>,
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

            let out = DadaOutput {
                num_asvs: asvs.len(),
                total_reads,
                asvs,
                stats: DadaStats {
                    nalign: result.nalign,
                    nshroud: result.nshroud,
                },
                map: if show_map { Some(result.map) } else { None },
            };

            let json = if compact {
                serde_json::to_string(&out)
            } else {
                serde_json::to_string_pretty(&out)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
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
            phred_offset,
            threads,
            output,
            compact,
            verbose,
        } => {
            // ---- Validate that all four lists have the same length ----
            let n = fwd_dada.len();
            for (flag, len) in [
                ("--rev-dada",   rev_dada.len()),
                ("--fwd-fastq",  fwd_fastq.len()),
                ("--rev-fastq",  rev_fastq.len()),
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
                        let name = p
                            .file_name()
                            .and_then(|n| n.to_str())
                            .unwrap_or("unknown");
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
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            let params = merge_pairs::MergeParams {
                min_overlap,
                max_mismatch,
                return_rejects,
                just_concatenate,
                concat_nnn_len,
                trim_overhang,
                phred_offset,
                verbose,
            };

            let mut results: Vec<merge_pairs::SampleMergeResult> = Vec::with_capacity(n);

            for i in 0..n {
                if verbose {
                    eprintln!(
                        "[merge-pairs] sample '{}' ({}/{})",
                        names[i],
                        i + 1,
                        n
                    );
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
                        names[i],
                        result.accepted_pairs,
                        result.total_pairs,
                        result.num_merged,
                    );
                }

                results.push(result);
            }

            let json = if compact {
                serde_json::to_string(&results)
            } else {
                serde_json::to_string_pretty(&results)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

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
            compress,
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
            threads,
            output,
            compact,
            verbose,
        } => {
            // ---- Validate file-list lengths ----
            let n = fwd.len();
            if filt.len() != n {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("--filt has {} entries but --fwd has {n}", filt.len()),
                ));
            }
            let paired = rev.is_some();
            let (rev_files, filt_rev_files) = if paired {
                let rv = rev.unwrap();
                let fr = filt_rev.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--filt-rev is required when --rev is provided",
                    )
                })?;
                if rv.len() != n {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("--rev has {} entries but --fwd has {n}", rv.len()),
                    ));
                }
                if fr.len() != n {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("--filt-rev has {} entries but --fwd has {n}", fr.len()),
                    ));
                }
                (rv, fr)
            } else {
                (vec![], vec![])
            };

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

            let (trunc_q_f, trunc_q_r)       = pair!(trunc_q, 2u8);
            let (trunc_len_f, trunc_len_r)   = pair!(trunc_len, 0usize);
            let (trim_left_f, trim_left_r)   = pair!(trim_left, 0usize);
            let (trim_right_f, trim_right_r) = pair!(trim_right, 0usize);
            let (max_len_f, max_len_r)       = pair!(max_len, 0usize);
            let (min_len_f, min_len_r)       = pair!(min_len, 20usize);
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
                trunc_q_f, trunc_len_f, trim_left_f, trim_right_f,
                max_len_f, min_len_f, max_ee_f, rm_lowcomplex_f,
            );
            let params_rev = make_params(
                trunc_q_r, trunc_len_r, trim_left_r, trim_right_r,
                max_len_r, min_len_r, max_ee_r, rm_lowcomplex_r,
            );

            // ---- Run (optionally parallel across samples) ----
            #[derive(Serialize)]
            struct SampleResult {
                sample: String,
                reads_in: u64,
                reads_out: u64,
            }

            let results: io::Result<Vec<SampleResult>> = if threads > 1 {
                use rayon::prelude::*;

                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                // Build index list, process in parallel.
                let indices: Vec<usize> = (0..n).collect();
                let par_results: Vec<io::Result<SampleResult>> = pool.install(|| {
                    indices.par_iter().map(|&i| {
                        let sample = fastq_stem(&fwd[i]);
                        let stats = if paired {
                            filter_paired(
                                &fwd[i], &rev_files[i],
                                &filt[i], &filt_rev_files[i],
                                &params_fwd, &params_rev,
                                compress, verbose,
                            )?
                        } else {
                            filter_single(&fwd[i], &filt[i], &params_fwd, compress, verbose)?
                        };
                        Ok(SampleResult { sample, reads_in: stats.reads_in, reads_out: stats.reads_out })
                    }).collect()
                });

                par_results.into_iter().collect()
            } else {
                let mut out = Vec::with_capacity(n);
                for i in 0..n {
                    let sample = fastq_stem(&fwd[i]);
                    let stats = if paired {
                        filter_paired(
                            &fwd[i], &rev_files[i],
                            &filt[i], &filt_rev_files[i],
                            &params_fwd, &params_rev,
                            compress, verbose,
                        )?
                    } else {
                        filter_single(&fwd[i], &filt[i], &params_fwd, compress, verbose)?
                    };
                    out.push(SampleResult { sample, reads_in: stats.reads_in, reads_out: stats.reads_out });
                }
                Ok(out)
            };

            let results = results?;

            let json = if compact {
                serde_json::to_string(&results)
            } else {
                serde_json::to_string_pretty(&results)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::MakeSequenceTable {
            input,
            sample_names,
            order_by,
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
                "nsamples"  => OrderBy::NSamples,
                _           => OrderBy::None,
            };
            let names_opt = if sample_names.is_empty() { None } else { Some(sample_names.as_slice()) };
            let hash_algo = if hash == "sha1" { HashAlgo::Sha1 } else { HashAlgo::Md5 };
            let paths: Vec<&Path> = input.iter().map(|p| p.as_path()).collect();
            let table = make_sequence_table(&paths, names_opt, order, hash_algo)?;
            let json = if compact {
                serde_json::to_string(&table)
            } else {
                serde_json::to_string_pretty(&table)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
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
            let bytes = std::fs::read(&input)?;
            let table: SequenceTable = serde_json::from_slice(&bytes)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let method = match method.as_str() {
                "pooled"     => Method::Pooled,
                "per-sample" => Method::PerSample,
                _            => Method::Consensus,
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
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            let filtered = pool.install(|| remove_bimera_denovo(table, &method, &params, verbose));

            let json = if compact {
                serde_json::to_string(&filtered)
            } else {
                serde_json::to_string_pretty(&filtered)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            match output {
                Some(path) => std::fs::write(&path, &json)?,
                None => println!("{json}"),
            }
        }

        Commands::SeqTableToTsv { input, output } => {
            let bytes = std::fs::read(&input)?;
            let table: SequenceTable = serde_json::from_slice(&bytes)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

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

            // One row per sequence: id <TAB> count_per_sample...
            for (j, id) in table.sequence_ids.iter().enumerate() {
                write!(out, "{id}")?;
                for sample_counts in &table.counts {
                    write!(out, "\t{}", sample_counts[j])?;
                }
                writeln!(out)?;
            }
            out.flush()?;
        }

        Commands::SeqTableToFasta { input, output } => {
            #[derive(serde::Deserialize)]
            struct SeqTable {
                sequences: Vec<String>,
                sequence_ids: Vec<String>,
            }

            let bytes = std::fs::read(&input)?;
            let table: SeqTable = serde_json::from_slice(&bytes)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            if table.sequences.len() != table.sequence_ids.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sequence_ids and sequences lengths differ",
                ));
            }

            let mut out: Box<dyn io::Write> = match output {
                Some(ref path) => Box::new(io::BufWriter::new(std::fs::File::create(path)?)),
                None => Box::new(io::BufWriter::new(std::io::stdout())),
            };

            for (id, seq) in table.sequence_ids.iter().zip(table.sequences.iter()) {
                writeln!(out, ">{id}\n{seq}")?;
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
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

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
                total_reads: usize,
                unique_sequences: usize,
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
                    dereplicate(MultiGzDecoder::new(File::open(path)?), phred_offset, &pool, verbose)?
                } else {
                    dereplicate(File::open(path)?, phred_offset, &pool, verbose)?
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
                        s1_path.file_stem().unwrap_or(s1).to_string_lossy().into_owned()
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
                    total_reads: derep.map.len(),
                    unique_sequences: uniq_entries.len(),
                    uniques: uniq_entries,
                };

                let json = if compact {
                    serde_json::to_string(&sample_out)
                } else {
                    serde_json::to_string_pretty(&sample_out)
                }
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                std::fs::write(&out_path, &json)?;
                output_files.push(out_path.display().to_string());
                total_bases += file_bases;
                total_reads += file_reads;

                if verbose {
                    eprintln!(
                        "[sample] wrote {} ({} unique(s), {} bases)",
                        out_path.display(),
                        sample_out.unique_sequences,
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
            let summary_json = if compact {
                serde_json::to_string(&summary)
            } else {
                serde_json::to_string_pretty(&summary)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            println!("{summary_json}");
        }

        Commands::ErrorsFromSample {
            input,
            errfun,
            pseudocount,
            binned_quals,
            max_consist,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            threads,
            output,
            compact,
            diag_dir,
            verbose,
        } => {
            let err_fun = match errfun.as_str() {
                "loess" => ErrFun::Loess,
                "noqual" => ErrFun::Noqual { pseudocount },
                "binned-qual" => {
                    let bins = binned_quals.ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--binned-quals is required when --errfun binned-qual is used",
                        )
                    })?;
                    ErrFun::BinnedQual { bins }
                }
                "pacbio" => ErrFun::PacBio,
                other => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown errfun '{other}'; expected one of: loess, noqual, binned-qual, pacbio"
                        ),
                    ));
                }
            };

            let align_params = AlignParams {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                homo_gap_p: -8,
                use_kmers: true,
                kdist_cutoff: 0.42,
                band: -1,
                vectorized: false,
                gapless: false,
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
            };

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

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

            let result = pool.install(|| {
                learn_errors(all_inputs, &err_fun, dada_params, &align_params, max_consist, verbose, diag_dir.as_deref())
            })?;

            #[derive(Serialize)]
            struct LearnErrorsOutput {
                nq: usize,
                converged: bool,
                iterations: usize,
                trans: Vec<Vec<u32>>,
                err_in: Vec<Vec<f64>>,
                err_out: Vec<Vec<f64>>,
            }

            fn flat_to_rows_u32(flat: &[u32], nq: usize) -> Vec<Vec<u32>> {
                (0..16).map(|r| flat[r * nq..(r + 1) * nq].to_vec()).collect()
            }
            fn flat_to_rows_f64(flat: &[f64], nq: usize) -> Vec<Vec<f64>> {
                (0..16).map(|r| flat[r * nq..(r + 1) * nq].to_vec()).collect()
            }

            let out = LearnErrorsOutput {
                nq: result.nq,
                converged: result.converged,
                iterations: result.iterations,
                trans: flat_to_rows_u32(&result.trans, result.nq),
                err_in: flat_to_rows_f64(&result.err_in, result.nq),
                err_out: flat_to_rows_f64(&result.err_out, result.nq),
            };

            let json = if compact {
                serde_json::to_string(&out)
            } else {
                serde_json::to_string_pretty(&out)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

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
            max_consist,
            omega_a,
            omega_c,
            omega_p,
            min_fold,
            min_hamming,
            min_abund,
            detect_singletons,
            threads,
            output,
            compact,
            diag_dir,
            verbose,
        } => {
            let err_fun = match errfun.as_str() {
                "loess" => ErrFun::Loess,
                "noqual" => ErrFun::Noqual { pseudocount },
                "binned-qual" => {
                    let bins = binned_quals.ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--binned-quals is required when --errfun binned-qual is used",
                        )
                    })?;
                    ErrFun::BinnedQual { bins }
                }
                "pacbio" => ErrFun::PacBio,
                other => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown errfun '{other}'; expected one of: loess, noqual, binned-qual, pacbio"
                        ),
                    ));
                }
            };

            let align_params = AlignParams {
                match_score: 5,
                mismatch: -4,
                gap_p: -8,
                homo_gap_p: -8,
                use_kmers: true,
                kdist_cutoff: 0.42,
                band: -1,
                vectorized: false,
                gapless: false,
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
            };

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            let all_inputs = load_fastq_samples(&input, nbases, randomize, seed, phred_offset, &pool, verbose)?;

            if let Some(ref dir) = diag_dir {
                std::fs::create_dir_all(dir)?;
            }

            let result = pool.install(|| {
                learn_errors(all_inputs, &err_fun, dada_params, &align_params, max_consist, verbose, diag_dir.as_deref())
            })?;

            // Serialize: represent the three matrices as Vec<Vec<T>> (16 rows × nq cols).
            #[derive(Serialize)]
            struct LearnErrorsOutput {
                nq: usize,
                converged: bool,
                iterations: usize,
                /// Transition counts: 16 rows (ref_nt*4+query_nt), nq columns.
                trans: Vec<Vec<u32>>,
                /// Error rates fed into the final DADA run: 16 × nq.
                err_in: Vec<Vec<f64>>,
                /// Error rates estimated from `trans`: 16 × nq.
                err_out: Vec<Vec<f64>>,
            }

            fn flat_to_rows_u32(flat: &[u32], nq: usize) -> Vec<Vec<u32>> {
                (0..16).map(|r| flat[r * nq..(r + 1) * nq].to_vec()).collect()
            }
            fn flat_to_rows_f64(flat: &[f64], nq: usize) -> Vec<Vec<f64>> {
                (0..16).map(|r| flat[r * nq..(r + 1) * nq].to_vec()).collect()
            }

            let out = LearnErrorsOutput {
                nq: result.nq,
                converged: result.converged,
                iterations: result.iterations,
                trans: flat_to_rows_u32(&result.trans, result.nq),
                err_in: flat_to_rows_f64(&result.err_in, result.nq),
                err_out: flat_to_rows_f64(&result.err_out, result.nq),
            };

            let json = if compact {
                serde_json::to_string(&out)
            } else {
                serde_json::to_string_pretty(&out)
            }
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

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

/// Derive a base stem from a FASTQ path by stripping recognised extensions.
///
/// `sample1.fastq.gz` → `"sample1"`,  `sample2.fq` → `"sample2"`.
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
