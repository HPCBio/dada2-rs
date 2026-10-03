//! Dereplication: `derep`, `sample`.

use std::{fs::File, io};

use flate2::read::MultiGzDecoder;
use rand::seq::SliceRandom as _;

use dada2_rs::{cli, derep, misc};

use derep::dereplicate;
use misc::Tagged;
use serde::Serialize;

use super::common::{check_input_paths, fastq_stem, file_basename};
use misc::WithPath;

pub(crate) fn run_derep(args: cli::DerepArgs) -> io::Result<()> {
    let cli::DerepArgs {
        input,
        sample_name,
        phred_offset,
        threads,
        output,
        show_map,
        pretty,
        verbose,
    } = args;
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
        /// Per-position integer Phred SUM; consumers divide by `count` on demand.
        qual_sum: &'a [u32],
    }

    #[derive(Serialize)]
    struct DerepOutput<'a> {
        sample: &'a str,
        /// Original input file name (no directory) for provenance.
        input_file: String,
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
        let sequence =
            std::str::from_utf8(seq).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        uniq_entries.push(UniqueEntry {
            sequence,
            count: *count,
            qual_sum: &derep.quals[i],
        });
    }

    let sample = sample_name.unwrap_or_else(|| fastq_stem(&input));
    let derep_out = DerepOutput {
        sample: &sample,
        input_file: file_basename(&input),
        total_reads: derep.map.len(),
        unique_sequences: derep.uniques.len(),
        sort_order: "abundance_desc",
        uniques: uniq_entries,
        map: if show_map { Some(&derep.map) } else { None },
    };

    let tagged = Tagged::new("derep", derep_out);
    let json = if pretty {
        serde_json::to_string_pretty(&tagged)
    } else {
        serde_json::to_string(&tagged)
    }
    .map_err(io::Error::other)?;

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_sample(args: cli::SampleArgs) -> io::Result<()> {
    let cli::SampleArgs {
        input,
        output_dir,
        nbases,
        randomize,
        seed,
        phred_offset,
        threads,
        pretty,
        gzip,
        verbose,
    } = args;
    check_input_paths("input", &input)?;
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
        /// Per-position integer Phred SUM; consumers divide by `count` on demand.
        qual_sum: &'a [u32],
    }
    #[derive(Serialize)]
    struct DerepOutput<'a> {
        sample: &'a str,
        /// Original input file name (no directory) for provenance.
        input_file: String,
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
        let ext = if gzip { "json.gz" } else { "json" };
        let out_path = output_dir.join(format!("{stem}.{ext}"));

        let mut uniq_entries = Vec::with_capacity(derep.uniques.len());
        for (i, (seq, count)) in derep.uniques.iter().enumerate() {
            let sequence = std::str::from_utf8(seq)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            uniq_entries.push(UniqueEntry {
                sequence,
                count: *count,
                qual_sum: &derep.quals[i],
            });
        }
        let sample_out = DerepOutput {
            sample: &stem,
            input_file: file_basename(path),
            total_reads: derep.map.len(),
            unique_sequences: uniq_entries.len(),
            sort_order: "abundance_desc",
            uniques: uniq_entries,
        };

        let unique_count = sample_out.unique_sequences;
        let tagged = Tagged::new("sample", sample_out);
        let json = if pretty {
            serde_json::to_string_pretty(&tagged)
        } else {
            serde_json::to_string(&tagged)
        }
        .map_err(io::Error::other)?;

        misc::write_maybe_gz(&out_path, json.as_bytes())?;
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
    let summary_json = if pretty {
        serde_json::to_string_pretty(&tagged)
    } else {
        serde_json::to_string(&tagged)
    }
    .map_err(io::Error::other)?;
    println!("{summary_json}");
    Ok(())
}
