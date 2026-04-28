//! Post-processing output computation.
//!
//! Ports `error.cpp`, excluding Rcpp data-frame/matrix constructors.
//! Each function returns plain Rust types; callers convert to their
//! preferred output format (JSON, TSV, etc.).
//!
//! All functions that accept `subs` expect a slice indexed by `raw.index`
//! (0-indexed input order), matching the layout produced by the final
//! parallel sub-computation pass in `dada_uniques`.

use std::collections::HashMap;

use crate::containers::{B, Sub};
use crate::misc::nt_decode;
use crate::nwalign::GAP_GLYPH;
use crate::pval::calc_pA;

// ---------------------------------------------------------------------------
// Cluster statistics (mirrors b_make_clustering_df)
// ---------------------------------------------------------------------------

/// Per-cluster summary statistics.
pub struct ClusterStats {
    /// Integer-encoded representative sequence.
    pub sequence: Vec<u8>,
    /// Total reads from `correct=true` members.
    pub abundance: u32,
    /// Reads from `correct=true` members with 0 substitutions vs the center.
    pub n0: u32,
    /// Reads from `correct=true` members with 1 substitution vs the center.
    pub n1: u32,
    /// Number of `correct=true` unique sequences in this cluster.
    pub nunq: u32,
    /// Post-hoc abundance p-value for this cluster's center.
    pub pval: f64,
    /// 0-indexed index of the parent cluster (`None` for cluster 0).
    pub birth_from: Option<usize>,
    pub birth_pval: f64,
    pub birth_fold: f64,
    pub birth_ham: u32,
    pub birth_e: f64,
    /// Mean quality score of substitutions that drove cluster birth.
    /// `None` when quality data is unavailable or there are no birth subs.
    pub birth_qave: Option<f64>,
}

/// Compute per-cluster summary statistics.
///
/// `subs[raw.index]` is the `Sub` between `bi.center` and `raw` for every
/// raw in the partition (computed with `use_kmers=false`).
/// `birth_subs[i]` is the `Sub` between the parent center and cluster `i`'s
/// center (`None` for cluster 0).
/// Equivalent to C++ `b_make_clustering_df` (minus the R data-frame wrapper)
/// plus the post-hoc p-value logic.
pub fn cluster_stats(
    b: &B,
    subs: &[Option<Sub>],
    birth_subs: &[Option<Sub>],
    has_quals: bool,
) -> Vec<ClusterStats> {
    let pvals = post_hoc_pvals(b);
    let nclust = b.clusters.len();

    (0..nclust)
        .map(|ci| {
            let bi = &b.clusters[ci];
            let sequence = bi.seq.clone();

            let mut abundance = 0u32;
            let mut n0 = 0u32;
            let mut n1 = 0u32;
            let mut nunq = 0u32;

            for &raw_idx in &bi.raws {
                let raw = &b.raws[raw_idx];
                if !raw.correct {
                    continue;
                }
                abundance += raw.reads;
                nunq += 1;
                if let Some(sub) = subs.get(raw.index as usize).and_then(|s| s.as_ref()) {
                    if sub.nsubs() == 0 {
                        n0 += raw.reads;
                    } else if sub.nsubs() == 1 {
                        n1 += raw.reads;
                    }
                }
            }

            let (birth_from, birth_pval, birth_fold, birth_ham, birth_e, birth_qave) = if ci == 0 {
                (None, 0.0, 0.0, 0u32, 0.0, None)
            } else {
                let bq = if has_quals {
                    birth_subs
                        .get(ci)
                        .and_then(|s| s.as_ref())
                        .filter(|s| !s.q1.is_empty())
                        .map(|s| {
                            let sum: f64 = s.q1.iter().map(|&q| q as f64).sum();
                            sum / s.nsubs() as f64
                        })
                } else {
                    None
                };
                (
                    Some(bi.birth_from as usize),
                    bi.birth_pval,
                    bi.birth_fold,
                    bi.birth_comp.hamming,
                    bi.birth_e,
                    bq,
                )
            };

            ClusterStats {
                sequence,
                abundance,
                n0,
                n1,
                nunq,
                pval: pvals[ci],
                birth_from,
                birth_pval,
                birth_fold,
                birth_ham,
                birth_e,
                birth_qave,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Post-hoc p-values
// ---------------------------------------------------------------------------

/// Compute the post-hoc abundance p-value for each cluster center.
///
/// For each cluster `j`, the expected read count is the sum of
/// `lambda * bi_reads` over all other clusters `i` whose comparison list
/// includes the center of cluster `j`.  The p-value is then
/// `P(X ≥ center_reads | Poisson(tot_e))` (unconditional).
///
/// Equivalent to the post-hoc pval block in C++ `b_make_clustering_df`.
fn post_hoc_pvals(b: &B) -> Vec<f64> {
    let nclust = b.clusters.len();

    // Map from raw.index → cluster index, for centers only.
    let center_of: HashMap<u32, usize> = b
        .clusters
        .iter()
        .enumerate()
        .filter_map(|(ci, bi)| bi.center.map(|c| (b.raws[c].index, ci)))
        .collect();

    // Sum expected contributions from other clusters.
    let mut tot_e = vec![0.0f64; nclust];
    for (ci, cluster) in b.clusters.iter().enumerate() {
        for comp in &cluster.comp {
            if let Some(&j) = center_of.get(&comp.index) {
                if j != ci {
                    tot_e[j] += comp.lambda * cluster.reads as f64;
                }
            }
        }
    }

    (0..nclust)
        .map(|ci| {
            if let Some(center_idx) = b.clusters[ci].center {
                calc_pA(b.raws[center_idx].reads, tot_e[ci], true)
            } else {
                1.0
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Transition counts (mirrors b_make_transition_by_quality_matrix)
// ---------------------------------------------------------------------------

/// Accumulate observed nucleotide transition counts, stratified by quality.
///
/// Returns a flat row-major matrix of shape `16 × ncol` stored as `Vec<u32>`,
/// where `result[t_ij * ncol + qual]` holds the count of transitions of type
/// `t_ij = ref_nti * 4 + query_nti` at quality score `qual`.
///
/// When `has_quals` is false, `ncol` is forced to 1 and all counts go into
/// column 0.
///
/// Only `correct=true` members of each cluster are counted.
/// Equivalent to C++ `b_make_transition_by_quality_matrix`.
pub fn transition_counts(b: &B, subs: &[Option<Sub>], has_quals: bool, ncol: usize) -> Vec<u32> {
    let ncol = if has_quals { ncol } else { 1 };
    let mut mat = vec![0u32; 16 * ncol];

    for ci in 0..b.clusters.len() {
        let center_idx = match b.clusters[ci].center {
            Some(c) => c,
            None => continue,
        };
        let center_seq = &b.raws[center_idx].seq;
        let center_len = center_seq.len();

        for &raw_idx in &b.clusters[ci].raws {
            let raw = &b.raws[raw_idx];
            if !raw.correct {
                continue;
            }
            let sub = match subs.get(raw.index as usize).and_then(|s| s.as_ref()) {
                Some(s) => s,
                None => continue,
            };

            for pos0 in 0..center_len {
                let pos1 = sub.map[pos0] as usize;
                if sub.map[pos0] == GAP_GLYPH {
                    continue;
                }
                let nti0 = center_seq[pos0].saturating_sub(1) as usize;
                let nti1 = raw.seq[pos1].saturating_sub(1) as usize;
                let qual = if has_quals {
                    raw.qual.as_ref().map_or(0, |q| q[pos1] as usize)
                } else {
                    0
                };
                let t_ij = 4 * nti0 + nti1;
                if t_ij < 16 && qual < ncol {
                    mat[t_ij * ncol + qual] += raw.reads;
                }
            }
        }
    }
    mat
}

// ---------------------------------------------------------------------------
// Cluster quality matrix (mirrors b_make_cluster_quality_matrix)
// ---------------------------------------------------------------------------

/// Compute read-weighted mean quality per position for each cluster.
///
/// Returns a `Vec` of length `nclust`, each element a `Vec<f64>` of length
/// `maxlen`.  Position `p` of cluster `ci` holds the mean quality at reference
/// position `p`, weighted by read counts.  Positions beyond the center length
/// or not covered by any correct read are `f64::NAN`.
///
/// Only `correct=true` members contribute.
/// Equivalent to C++ `b_make_cluster_quality_matrix`.
pub fn cluster_quality(
    b: &B,
    subs: &[Option<Sub>],
    has_quals: bool,
    maxlen: usize,
) -> Vec<Vec<f64>> {
    let nclust = b.clusters.len();
    let mut result: Vec<Vec<f64>> = vec![vec![f64::NAN; maxlen]; nclust];
    if !has_quals {
        return result;
    }

    for ci in 0..nclust {
        let center_idx = match b.clusters[ci].center {
            Some(c) => c,
            None => continue,
        };
        let center_len = b.raws[center_idx].seq.len().min(maxlen);
        let mut nreads = vec![0u32; center_len];
        let mut qual_sum = vec![0.0f64; center_len];

        for &raw_idx in &b.clusters[ci].raws {
            let raw = &b.raws[raw_idx];
            if !raw.correct {
                continue;
            }
            let sub = match subs.get(raw.index as usize).and_then(|s| s.as_ref()) {
                Some(s) => s,
                None => continue,
            };
            let qual = match &raw.qual {
                Some(q) => q,
                None => continue,
            };

            for pos0 in 0..center_len {
                if sub.map[pos0] == GAP_GLYPH {
                    continue;
                }
                let pos1 = sub.map[pos0] as usize;
                nreads[pos0] += raw.reads;
                qual_sum[pos0] += qual[pos1] as f64 * raw.reads as f64;
            }
        }

        for pos0 in 0..center_len {
            result[ci][pos0] = if nreads[pos0] > 0 {
                qual_sum[pos0] / nreads[pos0] as f64
            } else {
                f64::NAN
            };
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Birth substitution records (mirrors b_make_birth_subs_df)
// ---------------------------------------------------------------------------

/// One substitution from a cluster's birth alignment.
pub struct BirthSubRecord {
    /// 0-indexed position in the reference (parent center) sequence.
    pub pos: u16,
    /// Reference nucleotide (ASCII).
    pub nt0: u8,
    /// Query nucleotide (ASCII).
    pub nt1: u8,
    /// Quality score at this position in the query (if available).
    pub qual: Option<u8>,
    /// 0-indexed index of the cluster this substitution belongs to.
    pub cluster: usize,
}

/// Collect per-substitution records from each cluster's birth alignment.
///
/// Equivalent to C++ `b_make_birth_subs_df` (minus the R wrapper).
pub fn birth_sub_records(
    b: &B,
    birth_subs: &[Option<Sub>],
    has_quals: bool,
) -> Vec<BirthSubRecord> {
    let mut records = Vec::new();
    for (ci, maybe_sub) in birth_subs.iter().enumerate() {
        let sub = match maybe_sub {
            Some(s) => s,
            None => continue,
        };
        for s in 0..sub.nsubs() {
            records.push(BirthSubRecord {
                pos: sub.pos[s],
                nt0: nt_decode(sub.nt0[s]),
                nt1: nt_decode(sub.nt1[s]),
                qual: if has_quals {
                    sub.q1.get(s).copied()
                } else {
                    None
                },
                cluster: ci,
            });
        }
    }
    records
}
