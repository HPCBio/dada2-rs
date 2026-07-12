//! Build a sample-by-sequence feature table from `dada` or `merge-pairs` JSON output.
//!
//! Mirrors R's `makeSequenceTable`.  Each input file is either:
//!   - A single-sample `dada` JSON object (`{ "asvs": [...] }`), or
//!   - A multi-sample `merge-pairs` JSON array (`[{ "sample": "...", "merged": [...] }, ...]`).
//!
//! Output is a flat matrix:
//! ```json
//! { "samples": [...], "sequences": [...], "counts": [[...], ...] }
//! ```
//! Rows = samples, columns = sequences (ordered by decreasing total abundance by default).

use std::collections::HashMap;
use std::io;
use std::path::Path;

use md5::{Digest as _, Md5};
use serde::{Deserialize, Serialize};
pub use serde_json::Value;
use sha1::Sha1;

// ---------------------------------------------------------------------------
// Deserialisation helpers for the two supported input formats
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct DadaAsv {
    sequence: String,
    abundance: u64,
}

#[derive(Deserialize)]
struct DadaOutput {
    sample: Option<String>,
    asvs: Vec<DadaAsv>,
}

#[derive(Deserialize)]
struct MergeEntry {
    sequence: String,
    abundance: u64,
    accept: bool,
}

#[derive(Deserialize)]
struct SampleMergeResult {
    sample: String,
    merged: Vec<MergeEntry>,
}

// ---------------------------------------------------------------------------
// Internal representation: one sample's sequence → count map
// ---------------------------------------------------------------------------

struct SampleCounts {
    name: String,
    counts: HashMap<String, u64>,
}

// ---------------------------------------------------------------------------
// Input parsing
// ---------------------------------------------------------------------------

/// Parse one JSON file.  Returns one `SampleCounts` per sample found.
///
/// Dispatches on the `dada2_rs_command` tag:
/// * `"dada"` / `"dada-pooled"` / `"dada-pseudo"` → single-sample dada output (same schema)
/// * `"merge-pairs"` → multi-sample merge-pairs output (a `samples` array)
fn parse_file(path: &Path, sample_name: Option<&str>) -> io::Result<Vec<SampleCounts>> {
    let bytes = crate::misc::read_all_maybe_gz(path)?;
    let value: Value = serde_json::from_slice(&bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let cmd = value
        .get("dada2_rs_command")
        .and_then(|v| v.as_str())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: missing dada2_rs_command tag (expected 'dada', 'dada-pooled', 'dada-pseudo', or 'merge-pairs')",
                    path.display()
                ),
            )
        })?
        .to_owned();

    match cmd.as_str() {
        // "dada-pooled" and "dada-pseudo" emit the same per-sample DadaOutput
        // schema as "dada"; only the tag differs (it records pooling mode).
        "dada" | "dada-pooled" | "dada-pseudo" => {
            let out: DadaOutput = serde_json::from_value(value)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let name = sample_name
                .map(str::to_owned)
                .or(out.sample)
                .unwrap_or_else(|| file_stem(path));
            let counts = out
                .asvs
                .into_iter()
                .map(|a| (a.sequence, a.abundance))
                .collect();
            Ok(vec![SampleCounts { name, counts }])
        }
        "merge-pairs" => {
            #[derive(Deserialize)]
            struct MergePairsFile {
                samples: Vec<SampleMergeResult>,
            }
            let out: MergePairsFile = serde_json::from_value(value)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(out
                .samples
                .into_iter()
                .map(|r| SampleCounts {
                    name: r.sample,
                    counts: r
                        .merged
                        .into_iter()
                        .filter(|m| m.accept && !m.sequence.is_empty())
                        .map(|m| (m.sequence, m.abundance))
                        .collect(),
                })
                .collect())
        }
        other => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: dada2_rs_command={other:?}, expected 'dada', 'dada-pooled', or 'merge-pairs'",
                path.display()
            ),
        )),
    }
}

fn file_stem(path: &Path) -> String {
    path.file_stem()
        .and_then(|s| {
            // Strip a second extension if the stem itself ends in e.g. ".json"
            // so "sample1.dada.json" → "sample1.dada", which is fine.
            s.to_str()
        })
        .unwrap_or("unknown")
        .to_owned()
}

#[cfg(test)]
mod parse_tests {
    use super::parse_file;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    /// A gzipped `dada` output JSON (`*.json.gz`, e.g. from `dada --gzip` /
    /// `dada-pooled --gzip`) must be gunzipped before parsing — otherwise the
    /// gzip magic bytes hit serde as raw JSON and fail with "expected value,
    /// line 1 column 1". Regression for the make-sequence-table read path.
    #[test]
    fn parse_file_reads_gzipped_dada_output() {
        let dir = std::env::temp_dir().join(format!("dada2rs_seqtab_gz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("F3D0.json.gz");
        let mut gz = GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            Compression::default(),
        );
        write!(
            gz,
            r#"{{"dada2_rs_command":"dada","sample":"F3D0",
                 "asvs":[{{"sequence":"ACGTACGT","abundance":7}}]}}"#
        )
        .unwrap();
        gz.finish().unwrap();

        let got = parse_file(&path, None).expect("gzipped dada output should parse");
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].name, "F3D0");
        assert_eq!(got[0].counts.get("ACGTACGT"), Some(&7));

        let _ = std::fs::remove_dir_all(&dir);
    }
}

// ---------------------------------------------------------------------------
// Column ordering
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub enum HashAlgo {
    Md5,
    Sha1,
}

impl HashAlgo {
    pub fn digest(self, seq: &str) -> String {
        match self {
            HashAlgo::Md5 => format!("{:x}", Md5::digest(seq.as_bytes())),
            HashAlgo::Sha1 => format!("{:x}", Sha1::digest(seq.as_bytes())),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum OrderBy {
    Abundance,
    NSamples,
    None,
}

// ---------------------------------------------------------------------------
// Output type
// ---------------------------------------------------------------------------

#[derive(Serialize, Deserialize)]
pub struct SequenceTable {
    pub samples: Vec<String>,
    pub sequences: Vec<String>,
    /// Hash-based unique identifier for each sequence (parallel to `sequences`).
    pub sequence_ids: Vec<String>,
    /// counts[i][j] = count for samples[i], sequences[j]
    pub counts: Vec<Vec<u64>>,
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Build a `SequenceTable` from one or more input JSON files.
///
/// `sample_names` is an optional per-file name override; it only applies to
/// single-sample dada files (merge-pairs files carry names internally).
/// If provided its length must equal `inputs.len()`.
pub fn make_sequence_table(
    inputs: &[&Path],
    sample_names: Option<&[String]>,
    order_by: OrderBy,
    hash_algo: HashAlgo,
) -> io::Result<SequenceTable> {
    let mut all_samples: Vec<SampleCounts> = Vec::new();

    for (i, path) in inputs.iter().enumerate() {
        let name_override = sample_names.and_then(|ns| ns.get(i)).map(String::as_str);
        let mut parsed = parse_file(path, name_override)?;
        all_samples.append(&mut parsed);
    }

    // Collect all unique sequences in first-seen order.
    let mut seq_index: HashMap<String, usize> = HashMap::new();
    let mut sequences: Vec<String> = Vec::new();
    for s in &all_samples {
        for seq in s.counts.keys() {
            if !seq_index.contains_key(seq) {
                seq_index.insert(seq.clone(), sequences.len());
                sequences.push(seq.clone());
            }
        }
    }

    let nseq = sequences.len();
    let nsamp = all_samples.len();

    // Build the count matrix (samples × sequences).
    let mut counts: Vec<Vec<u64>> = vec![vec![0u64; nseq]; nsamp];
    for (i, s) in all_samples.iter().enumerate() {
        for (seq, &cnt) in &s.counts {
            let j = seq_index[seq];
            counts[i][j] += cnt;
        }
    }

    // Determine column order.
    let col_order: Vec<usize> = match order_by {
        OrderBy::None => (0..nseq).collect(),
        OrderBy::Abundance => {
            let totals: Vec<u64> = (0..nseq)
                .map(|j| counts.iter().map(|row| row[j]).sum())
                .collect();
            let mut order: Vec<usize> = (0..nseq).collect();
            order.sort_by(|&a, &b| totals[b].cmp(&totals[a]));
            order
        }
        OrderBy::NSamples => {
            let present: Vec<usize> = (0..nseq)
                .map(|j| counts.iter().filter(|row| row[j] > 0).count())
                .collect();
            let mut order: Vec<usize> = (0..nseq).collect();
            order.sort_by(|&a, &b| present[b].cmp(&present[a]));
            order
        }
    };

    // Reorder columns.
    let sequences: Vec<String> = col_order.iter().map(|&j| sequences[j].clone()).collect();
    let counts: Vec<Vec<u64>> = counts
        .into_iter()
        .map(|row| col_order.iter().map(|&j| row[j]).collect())
        .collect();

    let samples: Vec<String> = all_samples.into_iter().map(|s| s.name).collect();
    let sequence_ids: Vec<String> = sequences.iter().map(|s| hash_algo.digest(s)).collect();

    Ok(SequenceTable {
        samples,
        sequences,
        sequence_ids,
        counts,
    })
}
