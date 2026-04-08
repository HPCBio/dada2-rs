use std::{fs::File, io};

use clap::Parser;
use flate2::read::MultiGzDecoder;

mod cli;
mod containers;
mod derep;
mod kmers;
mod misc;
mod pval;
mod summary;

use cli::{Cli, Commands};
use derep::dereplicate;
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

        Commands::Derep { input, phred_offset, threads, show_map, compact } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

            let derep = if input.extension().and_then(|e| e.to_str()) == Some("gz") {
                dereplicate(MultiGzDecoder::new(File::open(&input)?), phred_offset, &pool)?
            } else {
                dereplicate(File::open(&input)?, phred_offset, &pool)?
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

            let output = DerepOutput {
                total_reads: derep.map.len(),
                unique_sequences: derep.uniques.len(),
                uniques: uniq_entries,
                map: if show_map { Some(&derep.map) } else { None },
            };

            let json = if compact {
                serde_json::to_string(&output)
            } else {
                serde_json::to_string_pretty(&output)
            }.map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            println!("{json}");
        }
    }

    Ok(())
}
