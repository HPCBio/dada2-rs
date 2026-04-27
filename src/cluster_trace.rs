//! Per-iteration cluster trace dump.
//!
//! Writes a tagged JSON document describing the full structure of one DADA
//! result — cluster centers, members with their hamming/λ/pval, birth
//! metadata — so downstream tools (R/Python plotting scripts) can inspect
//! how clustering evolves across self-consistency iterations or the final
//! ASV partition produced by `dada`.
//!
//! Sequences are deduplicated at the top-level `sequences` array; centers
//! and members reference them by integer index to keep file sizes
//! manageable on samples with many similar uniques.

use std::io;
use std::path::Path;

use serde::Serialize;

use crate::containers::BirthType;
use crate::dada::{DadaResult, RawInput};
use crate::misc::Tagged;

/// Knobs controlling what to include in the trace.  Defaults match the CLI:
/// every member written, no per-cluster sub matrices.
#[derive(Clone, Copy, Debug)]
pub struct TraceParams {
    /// Skip the `members` array entirely — emit only cluster centers and
    /// birth metadata.  Smallest files; sufficient for centers-vs-iter plots.
    pub no_members: bool,
    /// When `no_members` is false, only include members whose abundance is
    /// at least this value.  `1` means "include everything".
    pub min_abund: u32,
}

impl Default for TraceParams {
    fn default() -> Self {
        Self {
            no_members: false,
            min_abund: 1,
        }
    }
}

#[derive(Serialize)]
struct MemberJson {
    raw_seq_id: usize,
    abundance: u32,
    hamming: u32,
    lambda: f64,
    e_reads: f64,
    pval: f64,
}

#[derive(Serialize)]
struct ClusterJson {
    id: usize,
    center_seq_id: usize,
    abundance: u32,
    n_members: usize,
    birth_type: &'static str,
    birth_from: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    birth_pval: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    birth_fold: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    birth_e: Option<f64>,
    birth_hamming: u32,
    /// Number of members with abundance < `trace_min_abund` that were
    /// present in the cluster but excluded from the `members` array.
    /// Always 0 when `no_members` is false and `min_abund` is 1.
    members_excluded: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    members: Option<Vec<MemberJson>>,
}

#[derive(Serialize)]
struct ClusterTrace {
    sample: String,
    /// Iteration number (1-based) for `learn-errors` runs; `None` for the
    /// final-only `dada` run.
    #[serde(skip_serializing_if = "Option::is_none")]
    iteration: Option<usize>,
    nclust: usize,
    total_reads: u32,
    /// Knob values that produced this file; useful when consumers want
    /// to know whether members are complete.
    trace_no_members: bool,
    trace_min_abund: u32,
    nq: usize,
    /// Error matrix used as input to the dada run that produced these
    /// clusters (16 rows × `nq` cols, row-major).  `None` when not
    /// available (e.g. `dada` subcommand passes its err model directly).
    #[serde(skip_serializing_if = "Option::is_none")]
    err_in: Option<Vec<Vec<f64>>>,
    /// Deduplicated sequence pool referenced by `clusters[*].center_seq_id`
    /// and `clusters[*].members[*].raw_seq_id`.
    sequences: Vec<String>,
    clusters: Vec<ClusterJson>,
}

fn birth_label(b: &BirthType) -> &'static str {
    match b {
        BirthType::Initial => "Initial",
        BirthType::Abundance => "Abundance",
        BirthType::Prior => "Prior",
        BirthType::Singleton => "Singleton",
    }
}

/// Build the JSON document and write it to `path`.
///
/// `inputs` is the slice of `RawInput` originally passed to `dada_uniques`;
/// it provides the ASCII sequences and per-raw abundances that are
/// referenced by index.  `err_in` (if supplied) is a flat row-major
/// 16 × `nq` matrix.
#[allow(clippy::too_many_arguments)]
pub fn write_trace(
    path: &Path,
    sample: &str,
    iteration: Option<usize>,
    inputs: &[RawInput],
    result: &DadaResult,
    err_in: Option<&[f64]>,
    nq: usize,
    params: TraceParams,
    compact: bool,
) -> io::Result<()> {
    use std::collections::HashMap;

    // Deduplicate sequences: centers added first, then any member raw
    // sequences not already present.  String-keyed because center sequences
    // come from Vec<u8> (need decoding) while RawInput.seq is already String.
    let mut sequences: Vec<String> = Vec::new();
    let mut intern_table: HashMap<String, usize> = HashMap::new();

    let intern =
        |s: &str, table: &mut HashMap<String, usize>, pool: &mut Vec<String>| -> usize {
            if let Some(&id) = table.get(s) {
                return id;
            }
            let id = pool.len();
            pool.push(s.to_string());
            table.insert(s.to_string(), id);
            id
        };

    let total_reads: u32 = result.clusters.iter().map(|c| c.reads).sum();

    let mut cluster_jsons: Vec<ClusterJson> = Vec::with_capacity(result.clusters.len());
    for (ci, c) in result.clusters.iter().enumerate() {
        let center_str = decode_seq(&c.sequence);
        let center_seq_id = intern(&center_str, &mut intern_table, &mut sequences);

        let mut members_excluded = 0usize;
        let members_field = if params.no_members {
            None
        } else {
            let mut out: Vec<MemberJson> = Vec::with_capacity(c.members.len());
            for k in 0..c.members.len() {
                let raw_idx = c.members[k];
                let abund = inputs[raw_idx].abundance;
                if abund < params.min_abund {
                    members_excluded += 1;
                    continue;
                }
                let seq_id = intern(&inputs[raw_idx].seq, &mut intern_table, &mut sequences);
                out.push(MemberJson {
                    raw_seq_id: seq_id,
                    abundance: abund,
                    hamming: c.member_hammings[k],
                    lambda: c.member_lambdas[k],
                    e_reads: c.member_lambdas[k] * c.reads as f64,
                    pval: c.member_pvals[k],
                });
            }
            Some(out)
        };

        let (birth_pval, birth_fold, birth_e) = match c.birth_type {
            BirthType::Initial => (None, None, None),
            _ => (Some(c.birth_pval), Some(c.birth_fold), Some(c.birth_e)),
        };

        cluster_jsons.push(ClusterJson {
            id: ci,
            center_seq_id,
            abundance: c.reads,
            n_members: c.members.len(),
            birth_type: birth_label(&c.birth_type),
            birth_from: c.birth_from,
            birth_pval,
            birth_fold,
            birth_e,
            birth_hamming: c.birth_hamming,
            members_excluded,
            members: members_field,
        });
    }

    let err_in_rows = err_in.map(|flat| {
        (0..16)
            .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
            .collect::<Vec<Vec<f64>>>()
    });

    let trace = ClusterTrace {
        sample: sample.to_string(),
        iteration,
        nclust: result.clusters.len(),
        total_reads,
        trace_no_members: params.no_members,
        trace_min_abund: params.min_abund,
        nq,
        err_in: err_in_rows,
        sequences,
        clusters: cluster_jsons,
    };

    let tagged = Tagged::new("cluster-trace", trace);
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(|e| io::Error::other(e))?;
    std::fs::write(path, json)?;
    Ok(())
}

/// Decode an integer-encoded sequence (A=1,C=2,G=3,T=4,N=5) back to ASCII.
fn decode_seq(seq: &[u8]) -> String {
    seq.iter()
        .map(|&b| match b {
            1 => 'A',
            2 => 'C',
            3 => 'G',
            4 => 'T',
            5 => 'N',
            _ => '?',
        })
        .collect()
}
