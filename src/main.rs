use std::{fs::File, io, path::Path};

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rand::seq::SliceRandom as _;

mod chimera;
mod cli;
mod cluster;
mod containers;
mod dada;
mod derep;
mod error;
mod error_models;
mod evaluate;
mod filter;
mod kmers;
mod learn_errors;
mod misc;
mod nwalign;
mod pval;
mod summary;
mod taxonomy;

use cli::{Cli, Commands};
use derep::dereplicate;
use learn_errors::{ErrFun, learn_errors};
use nwalign::AlignParams;
use serde::Serialize;
use summary::process;

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Summary { input, phred_offset, threads } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            let summary = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                process(MultiGzDecoder::new(File::open(&input)?), phred_offset, &pool)?
            } else {
                process(File::open(&input)?, phred_offset, &pool)?
            };

            println!("Total reads: {}", summary.total_reads);
            println!("Mean quality per position (Phred):");
            println!("{:<8} {:>10}", "Position", "MeanQ");
            for (pos, mean) in summary.mean_quality_per_position().iter().enumerate() {
                println!("{:<8} {:>10.2}", pos + 1, mean);
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

        Commands::Subsample { input_dir, output_dir, nbases, randomize, phred_offset, threads, verbose } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            // Discover FASTQ files in the input directory.
            let mut fastq_files = collect_fastq_files(&input_dir)?;
            if fastq_files.is_empty() {
                eprintln!("No FASTQ files found in {}", input_dir.display());
                return Ok(());
            }

            if randomize {
                fastq_files.shuffle(&mut rand::thread_rng());
            } else {
                fastq_files.sort();
            }

            std::fs::create_dir_all(&output_dir)?;

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

            let mut total_bases: u64 = 0;
            let mut total_reads: u64 = 0;
            let mut files_written: usize = 0;

            if verbose {
                eprintln!(
                    "[subsample] {} FASTQ file(s) found; target {} bases",
                    fastq_files.len(),
                    nbases
                );
            }

            for path in &fastq_files {
                let derep = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
                    dereplicate(MultiGzDecoder::new(File::open(path)?), phred_offset, &pool, verbose)?
                } else {
                    dereplicate(File::open(path)?, phred_offset, &pool, verbose)?
                };

                // Bases for this file = sum(count * sequence_length) over uniques.
                let file_bases: u64 = derep.uniques.iter()
                    .map(|(seq, count)| seq.len() as u64 * count)
                    .sum();
                let file_reads = derep.map.len() as u64;

                // Serialize to JSON and write.
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
                let out_data = DerepOutput {
                    total_reads: derep.map.len(),
                    unique_sequences: derep.uniques.len(),
                    uniques: uniq_entries,
                };
                let json = serde_json::to_string_pretty(&out_data)
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                let stem = fastq_stem(path);
                let out_path = output_dir.join(format!("{stem}.json"));
                std::fs::write(&out_path, json)?;

                total_bases += file_bases;
                total_reads += file_reads;
                files_written += 1;

                if verbose {
                    eprintln!(
                        "[subsample] {} -> {} ({} reads, {} bases; {}/{} bases total)",
                        path.display(),
                        out_path.display(),
                        file_reads,
                        file_bases,
                        total_bases,
                        nbases,
                    );
                }

                // Mirror R's `if(NBASES > nbases) { break }` — strict greater-than
                // so a file that lands exactly on the target is still included.
                if total_bases > nbases {
                    if verbose {
                        eprintln!("[subsample] base target reached; stopping.");
                    }
                    break;
                }
            }

            eprintln!(
                "[subsample] done: {} total bases from {} total reads across {} file(s)",
                total_bases,
                total_reads,
                files_written,
            );
        }

        Commands::LearnErrors {
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
            output,
            compact,
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
                multithread: false,
                verbose,
                greedy: true,
            };

            let result = learn_errors(&input, &err_fun, dada_params, &align_params, max_consist, verbose)?;

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
fn collect_fastq_files(dir: &Path) -> io::Result<Vec<std::path::PathBuf>> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file() && is_fastq(&path) {
            files.push(path);
        }
    }
    Ok(files)
}

/// Return true when `path` has a recognised FASTQ extension.
fn is_fastq(path: &Path) -> bool {
    let name = match path.file_name().and_then(|n| n.to_str()) {
        Some(n) => n,
        None => return false,
    };
    name.ends_with(".fastq")
        || name.ends_with(".fastq.gz")
        || name.ends_with(".fq")
        || name.ends_with(".fq.gz")
}

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
