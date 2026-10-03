//! Taxonomy: `assign-taxonomy`, `assign-species`, `tax-to-tsv`.

use std::{io, path::Path};

use dada2_rs::{cli, misc, taxonomy};

use misc::{Tagged, read_fasta_records, read_tagged_json};
use serde::Serialize;
use taxonomy::{
    SpeciesHit, SpeciesOptions, SpeciesRef, TaxonomyOptions, TaxonomyRef, assign_species,
    assign_taxonomy,
};

use misc::WithPath;

pub(crate) fn run_assign_taxonomy(args: cli::AssignTaxonomyArgs) -> io::Result<()> {
    let cli::AssignTaxonomyArgs {
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
    } = args;
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
        let mut level_vals: Vec<String> = genus_fields.iter().map(|f| f[l].clone()).collect();
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
                let fields: Vec<&str> = genus_fields[g].iter().map(|s| s.as_str()).collect();
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
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_assign_species(args: cli::AssignSpeciesArgs) -> io::Result<()> {
    let cli::AssignSpeciesArgs {
        input,
        ref_fasta,
        allow_multiple,
        try_rc,
        output,
        compact,
        verbose,
    } = args;
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
    let tax_in: AssignTaxIn = read_tagged_json(&input, &["assign-taxonomy"]).with_path(&input)?;

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
                (Some(gi), Some(hit_gen)) => match a.taxonomy.get(gi).and_then(|x| x.as_deref()) {
                    Some(assigned_gen) if match_genera(assigned_gen, hit_gen) => hit.species,
                    _ => None,
                },
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
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_tax_to_tsv(args: cli::TaxToTsvArgs) -> io::Result<()> {
    let cli::TaxToTsvArgs {
        input,
        na_string,
        output,
    } = args;
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

    let tax: TaxJson =
        read_tagged_json(&input, &["assign-taxonomy", "assign-species"]).with_path(&input)?;

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
