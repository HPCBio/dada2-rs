use std::{fs::File, io};

use clap::Parser;
use flate2::read::MultiGzDecoder;

mod cli;
mod summary;

use cli::{Cli, Commands};
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
    }

    Ok(())
}
