//! Sequence tables: `make-sequence-table`, `seq-table-to-tsv`, `seq-table-to-fasta`.

use std::{io, path::Path};

use dada2_rs::{cli, misc, sequence_table};

use misc::{Tagged, read_tagged_json};
use sequence_table::{HashAlgo, OrderBy, SequenceTable, make_sequence_table};

use super::common::{check_input_paths, select_sequences};
use misc::WithPath;

pub(crate) fn run_make_sequence_table(args: cli::MakeSequenceTableArgs) -> io::Result<()> {
    let cli::MakeSequenceTableArgs {
        input,
        sample_names,
        order_by,
        min_len,
        max_len,
        hash,
        output,
        compact,
    } = args;
    check_input_paths("input", &input)?;
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
                min_len.is_none_or(|mn| s.len() >= mn) && max_len.is_none_or(|mx| s.len() <= mx)
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
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_seq_table_to_tsv(args: cli::SeqTableToTsvArgs) -> io::Result<()> {
    let cli::SeqTableToTsvArgs {
        input,
        prevalence,
        min_abundance,
        output,
    } = args;
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
    Ok(())
}

pub(crate) fn run_seq_table_to_fasta(args: cli::SeqTableToFastaArgs) -> io::Result<()> {
    let cli::SeqTableToFastaArgs {
        input,
        prevalence,
        min_abundance,
        output,
    } = args;
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
    Ok(())
}
