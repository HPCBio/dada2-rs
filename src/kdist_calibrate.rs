//! k-mer-distance vs true-divergence calibration (hidden `kdist-calibrate` subcommand).
//!
//! DADA2's k-mer screen skips alignment for pairs with k-mer distance above
//! `KDIST_CUTOFF = 0.42` (nominally ~10% nucleotide divergence, calibrated on
//! Illumina 16S). The k-mer distance traces to ESPRIT (Sun et al. 2009, whose
//! reference implementation is gone); the 0.42/10% calibration to DADA2
//! (Callahan et al. 2016). This re-derives the relationship empirically on real
//! data: for sampled unique-sequence pairs it emits the k-mer distance
//! (`kmer_dist8`, our port of the ESPRIT metric) alongside the true UNBANDED
//! `align_endsfree` divergence, so the constant can be checked per dataset /
//! platform / k / pooling regime.
//!
//! Alignment is unbanded by default (`--band -1`): the curve must measure the
//! true divergence of distant pairs, which a band would truncate. It is the
//! expensive part, so it is parallelised across `--threads`.
//!
//! POOLING: with multiple derep JSONs, `--per-sample` computes pairs WITHIN each
//! sample (the per-sample / independent regime); the default pools all uniques
//! into one set and computes pairs across the union (the full-pool regime).
//! Pseudo's screen population is per-sample (priors change the partition, not
//! which pairs are screened), so model it with `--per-sample`.

use crate::kmers::{assign_kmer8, kmer_dist8};
use crate::misc::WithPath as _;
use crate::nwalign::{AlignBuffers, align_endsfree_with_buf};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const GAP: u8 = b'-';

/// Parameters for [`run`] (mirrors the CLI flags).
pub struct Params {
    pub k: usize,
    pub cutoff: f64,
    pub leak_pct: f64,
    pub band: i32,
    pub max_pairs: usize,
    pub max_uniques: usize,
    pub per_sample: bool,
    pub nearest_parent: bool,
    /// Post-inference mode: positional inputs are `dada` output JSONs; each is
    /// paired with its derep JSON (from `derep_dir`) so every input unique can be
    /// labelled by what denoising did to it (center / member / failed).
    pub from_dada: bool,
    /// Directory holding the derep JSONs that fed `dada` (matched by `sample`
    /// name → `{sample}.json` / `.json.gz`). Required with `from_dada`.
    pub derep_dir: Option<PathBuf>,
    pub threads: usize,
    pub seed: u64,
    pub output: Option<PathBuf>,
    pub verbose: bool,
}

fn encode(seq: &str) -> Vec<u8> {
    seq.bytes()
        .map(|b| match b {
            b'A' | b'a' => 1,
            b'C' | b'c' => 2,
            b'G' | b'g' => 3,
            b'T' | b't' => 4,
            _ => 5, // N etc. — never a valid k-mer, never matches
        })
        .collect()
}

/// Stats of an ends-free alignment. Returns (edits, core_len, band_req):
/// - edits/core_len: substitution+indel columns over the aligned core (terminal
///   gap overhang trimmed — that's length difference, not divergence).
/// - band_req: the MINIMUM diagonal band that would reproduce this alignment =
///   max over the path of |#seq1-bases − #seq2-bases consumed so far|. A banded
///   aligner with band < band_req cannot find this alignment (it gets truncated
///   to a worse score), so band_req is the cost of correctly aligning the pair.
fn aln_divergence(a: &[u8], b: &[u8]) -> (usize, usize, usize) {
    let n = a.len();
    let mut lo = 0;
    while lo < n && (a[lo] == GAP || b[lo] == GAP) {
        lo += 1;
    }
    let mut hi = n;
    while hi > lo && (a[hi - 1] == GAP || b[hi - 1] == GAP) {
        hi -= 1;
    }
    // band_req over the full path (the band applies to the whole DP matrix):
    // a gap in seq2 advances only seq1 (offset +1); a gap in seq1, only seq2
    // (offset -1); a match/mismatch advances both (no change).
    let (mut off, mut band_req) = (0i32, 0i32);
    for k in 0..n {
        if b[k] == GAP {
            off += 1;
        } else if a[k] == GAP {
            off -= 1;
        }
        band_req = band_req.max(off.abs());
    }
    let mut edits = 0;
    for k in lo..hi {
        if a[k] == GAP || b[k] == GAP || a[k] != b[k] {
            edits += 1;
        }
    }
    (edits, hi - lo, band_req as usize)
}

/// One sample's uniques: encoded sequences, abundances, and k-mer screens.
struct Sample {
    name: String,
    enc: Vec<Vec<u8>>,
    counts: Vec<u64>,
    kmers: Vec<Vec<u8>>,
}

/// Load a derep JSON (`uniques[].sequence` + `count`), gzip-transparent. Only
/// derep JSONs are accepted (the screen operates on per-sample uniques).
fn load_derep(path: &Path, k: usize, max_uniques: usize, seed: u64) -> io::Result<Sample> {
    let f = File::open(path).with_path(path)?;
    let mut txt = String::new();
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        MultiGzDecoder::new(f).read_to_string(&mut txt)?;
    } else {
        io::BufReader::new(f).read_to_string(&mut txt)?;
    }
    let v: serde_json::Value =
        serde_json::from_str(&txt).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let uniques = v.get("uniques").and_then(|u| u.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: not a derep JSON (no `uniques`)", path.display()),
        )
    })?;
    let name = v
        .get("sample")
        .and_then(|s| s.as_str())
        .map(str::to_string)
        .unwrap_or_else(|| {
            path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("?")
                .to_string()
        });
    let mut seqs: Vec<(Vec<u8>, u64)> = uniques
        .iter()
        .filter_map(|e| {
            let s = e.get("sequence")?.as_str()?;
            let c = e.get("count").and_then(|c| c.as_u64()).unwrap_or(1);
            Some((encode(s), c))
        })
        .collect();
    // Optional per-sample subsample of uniques to bound the O(n^2) pair count.
    // Random (not abundance-top) keeps the divergence distribution unbiased.
    if max_uniques > 0 && seqs.len() > max_uniques {
        let mut st = seed ^ (seqs.len() as u64).wrapping_mul(0x9E37_79B9);
        // partial Fisher–Yates: move `max_uniques` random picks to the front.
        for i in 0..max_uniques {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            let j = i + ((st >> 33) as usize) % (seqs.len() - i);
            seqs.swap(i, j);
        }
        seqs.truncate(max_uniques);
    }
    let (enc, counts): (Vec<_>, Vec<_>) = seqs.into_iter().unzip();
    let kmers = enc.iter().map(|e| assign_kmer8(e, k)).collect();
    Ok(Sample {
        name,
        enc,
        counts,
        kmers,
    })
}

/// A `dada` output paired with its derep input, ready for post-inference
/// classification. Input uniques are in denoising input order (abundance-desc),
/// so `map[i]` is the cluster index Raw `i` was corrected to (`None` = failed
/// the abundance test, i.e. did not survive denoising). Centers are indexed by
/// cluster id, carrying the birth metadata that lets us trace *why* an ASV
/// exists (`Prior` = born from a pseudo-pool prior; `birth_pval` = how close
/// that birth was to OMEGA_A).
struct DadaSample {
    name: String,
    // input uniques (denoising input order, aligned with `map`)
    enc: Vec<Vec<u8>>,
    counts: Vec<u64>,
    kmers: Vec<Vec<u8>>,
    map: Vec<Option<usize>>,
    // cluster centers, indexed by cluster id
    c_enc: Vec<Vec<u8>>,
    c_kmers: Vec<Vec<u8>>,
    c_ab: Vec<u64>,
    c_birth: Vec<String>,
    c_birth_pval: Vec<f64>,
}

/// Read a derep JSON's uniques in denoising input order (abundance-descending,
/// matching `load_derep_for_dada`): no subsample — indices must line up with the
/// dada `map`. Returns (encoded seqs, counts).
fn load_derep_aligned(path: &Path) -> io::Result<(Vec<Vec<u8>>, Vec<u64>)> {
    let f = File::open(path).with_path(path)?;
    let mut txt = String::new();
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        MultiGzDecoder::new(f).read_to_string(&mut txt)?;
    } else {
        io::BufReader::new(f).read_to_string(&mut txt)?;
    }
    let v: serde_json::Value =
        serde_json::from_str(&txt).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let uniques = v.get("uniques").and_then(|u| u.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: not a derep JSON (no `uniques`)", path.display()),
        )
    })?;
    let mut seqs: Vec<(Vec<u8>, u64)> = uniques
        .iter()
        .filter_map(|e| {
            let s = e.get("sequence")?.as_str()?;
            let c = e.get("count").and_then(|c| c.as_u64()).unwrap_or(1);
            Some((encode(s), c))
        })
        .collect();
    // Mirror load_derep_for_dada: only re-sort when the producer didn't declare
    // abundance_desc; the sort is stable so ties keep file order (= input order).
    if v.get("sort_order").and_then(|s| s.as_str()) != Some("abundance_desc") {
        seqs.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    }
    Ok(seqs.into_iter().unzip())
}

/// Load a `dada` output JSON and pair it with its derep input from `derep_dir`.
fn load_dada(dada_path: &Path, derep_dir: &Path, k: usize) -> io::Result<DadaSample> {
    let f = File::open(dada_path).with_path(dada_path)?;
    let v: serde_json::Value = serde_json::from_reader(io::BufReader::new(f))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let name = v
        .get("sample")
        .and_then(|s| s.as_str())
        .map(str::to_string)
        .unwrap_or_else(|| {
            dada_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("?")
                .to_string()
        });
    let asvs = v.get("asvs").and_then(|a| a.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: not a dada output JSON (no `asvs`)",
                dada_path.display()
            ),
        )
    })?;
    let (mut c_enc, mut c_ab, mut c_birth, mut c_birth_pval) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for a in asvs {
        let seq = a.get("sequence").and_then(|s| s.as_str()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "asv entry missing `sequence`")
        })?;
        c_enc.push(encode(seq));
        c_ab.push(a.get("abundance").and_then(|x| x.as_u64()).unwrap_or(0));
        c_birth.push(
            a.get("birth_type")
                .and_then(|s| s.as_str())
                .unwrap_or("?")
                .to_string(),
        );
        c_birth_pval.push(
            a.get("birth_pval")
                .and_then(|x| x.as_f64())
                .unwrap_or(f64::NAN),
        );
    }
    let map: Vec<Option<usize>> = v
        .get("map")
        .and_then(|m| m.as_array())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "dada output missing `map`"))?
        .iter()
        .map(|e| e.as_u64().map(|u| u as usize))
        .collect();

    // Locate and load the matching derep input (sample.json[.gz]).
    let derep_path = [
        derep_dir.join(format!("{name}.json")),
        derep_dir.join(format!("{name}.json.gz")),
    ]
    .into_iter()
    .find(|p| p.exists())
    .ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "{}: no derep JSON for sample `{name}` in {}",
                dada_path.display(),
                derep_dir.display()
            ),
        )
    })?;
    let (enc, counts) = load_derep_aligned(&derep_path)?;
    if enc.len() != map.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{name}: derep has {} uniques but dada map has {} entries — \
                 mismatched inputs (regenerate from the same derep)",
                enc.len(),
                map.len()
            ),
        ));
    }
    let kmers = enc.iter().map(|e| assign_kmer8(e, k)).collect();
    let c_kmers = c_enc.iter().map(|e| assign_kmer8(e, k)).collect();
    Ok(DadaSample {
        name,
        enc,
        counts,
        kmers,
        map,
        c_enc,
        c_kmers,
        c_ab,
        c_birth,
        c_birth_pval,
    })
}

/// Build the (i, j) pair list for a population of `n` uniques: enumerate all if
/// `n*(n-1)/2 <= max_pairs`, else draw `max_pairs` random pairs (with possible
/// repeats — fine for a calibration scatter).
fn pairs_for(n: usize, max_pairs: usize, seed: u64) -> Vec<(usize, usize)> {
    let total = n.saturating_mul(n.saturating_sub(1)) / 2;
    if n < 2 {
        return Vec::new();
    }
    if total <= max_pairs {
        let mut v = Vec::with_capacity(total);
        for i in 0..n {
            for j in (i + 1)..n {
                v.push((i, j));
            }
        }
        return v;
    }
    let mut st = seed;
    let mut rnd = |m: usize| {
        st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((st >> 33) as usize) % m
    };
    (0..max_pairs)
        .map(|_| {
            let i = rnd(n);
            let mut j = rnd(n);
            if i == j {
                j = (j + 1) % n;
            }
            (i.min(j), i.max(j))
        })
        .collect()
}

pub fn run(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    if inputs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no input derep JSON(s) given",
        ));
    }
    if p.from_dada {
        return run_from_dada(inputs, p);
    }
    let loaded: Vec<Sample> = inputs
        .iter()
        .map(|path| load_derep(path, p.k, p.max_uniques, p.seed))
        .collect::<io::Result<_>>()?;

    // Form populations: one per sample (per-sample) or a single merged pool.
    let pops: Vec<Sample> = if p.per_sample {
        loaded
    } else {
        let mut pool = Sample {
            name: "pool".into(),
            enc: Vec::new(),
            counts: Vec::new(),
            kmers: Vec::new(),
        };
        for s in loaded {
            pool.enc.extend(s.enc);
            pool.counts.extend(s.counts);
            pool.kmers.extend(s.kmers);
        }
        vec![pool]
    };

    let mut w: Box<dyn Write> = match &p.output {
        Some(path) => Box::new(BufWriter::new(File::create(path).with_path(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads)
        .build()
        .map_err(io::Error::other)?;

    if p.nearest_parent {
        nearest_parent_mode(&mut *w, &pool, &pops, p)?;
    } else {
        pairs_mode(&mut *w, &pool, &pops, p)?;
    }
    w.flush()?;
    Ok(())
}

/// All-pairs (random-subsampled) mode: kdist vs divergence over sampled pairs.
fn pairs_mode(
    w: &mut dyn Write,
    pool: &rayon::ThreadPool,
    pops: &[Sample],
    p: &Params,
) -> io::Result<()> {
    writeln!(
        w,
        "sample,kdist,edits,core_len,pct_div,band_req,screened_in,ab_i,ab_j"
    )?;
    let (mut tot, mut scr, mut leak) = (0u64, 0u64, 0u64);
    let mut band_all: Vec<usize> = Vec::new();
    let mut band_scr: Vec<usize> = Vec::new();
    for (pi, s) in pops.iter().enumerate() {
        let n = s.enc.len();
        let seed = p.seed.wrapping_add(pi as u64).wrapping_mul(0x100000001B3);
        let pairs = pairs_for(n, p.max_pairs, seed);
        if p.verbose {
            let total = n.saturating_mul(n.saturating_sub(1)) / 2;
            eprintln!(
                "[kdist] {} : {n} uniques, {total} pairs -> {} computed (k={}, band={}, {} threads)",
                s.name,
                pairs.len(),
                p.k,
                p.band,
                p.threads,
            );
        }
        let rows: Vec<(f64, usize, usize, f64, usize, bool)> = pool.install(|| {
            pairs
                .par_iter()
                .map_init(AlignBuffers::new, |buf, &(i, j)| {
                    let kd = kmer_dist8(
                        &s.kmers[i],
                        s.enc[i].len(),
                        &s.kmers[j],
                        s.enc[j].len(),
                        p.k,
                    );
                    align_endsfree_with_buf(&s.enc[i], &s.enc[j], 5, -4, -8, p.band, buf);
                    let (edits, core, band_req) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (kd, edits, core, pct, band_req, kd < p.cutoff)
                })
                .collect()
        });
        for (idx, &(kd, edits, core, pct, band_req, sin)) in rows.iter().enumerate() {
            let (i, j) = pairs[idx];
            tot += 1;
            band_all.push(band_req);
            if sin {
                scr += 1;
                band_scr.push(band_req);
                if pct > p.leak_pct {
                    leak += 1;
                }
            }
            writeln!(
                w,
                "{},{kd:.4},{edits},{core},{pct:.3},{band_req},{},{},{}",
                s.name, sin as u8, s.counts[i], s.counts[j]
            )?;
        }
    }
    if p.verbose && tot > 0 {
        eprintln!(
            "[kdist] {tot} pairs: screened-in (kdist<{}) {scr} ({:.1}%); of those {leak} are >{}% divergent (leakage)",
            p.cutoff,
            100.0 * scr as f64 / tot as f64,
            p.leak_pct,
        );
        eprintln!("[kdist] {}", band_fit("all pairs", &band_all));
        eprintln!("[kdist] {}", band_fit("screened-in", &band_scr));
    }
    Ok(())
}

/// Post-inference driver: load each (dada output, derep input) pair and emit the
/// labelled pairwise comparisons.
fn run_from_dada(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    let derep_dir = p.derep_dir.as_deref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "--from-dada requires --derep-dir <DIR> (where the matching derep JSONs live)",
        )
    })?;
    let samples: Vec<DadaSample> = inputs
        .iter()
        .map(|path| load_dada(path, derep_dir, p.k))
        .collect::<io::Result<_>>()?;

    let mut w: Box<dyn Write> = match &p.output {
        Some(path) => Box::new(BufWriter::new(File::create(path).with_path(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads)
        .build()
        .map_err(io::Error::other)?;
    from_dada_mode(&mut *w, &pool, &samples, p)?;
    w.flush()?;
    Ok(())
}

/// A row of work: align input unique `i` (or center `i` for center pairs) against
/// a partner center `c`, under a class label.
enum Job {
    /// Member `i` absorbed by its own center `c` (an error copy denoising fixed).
    Member { i: usize, c: usize },
    /// Failed unique `i` (map==None) vs its nearest center `c` by k-mer distance.
    Failed { i: usize, c: usize },
    /// Two surviving ASV centers `a`,`b` (the inter-ASV resolution floor).
    CenterPair { a: usize, b: usize },
}

/// Post-inference classification mode. For every input unique, denoising decided
/// one of three fates; we align it against the relevant center and tag the row so
/// the divergence/k-mer-distance can be read *conditioned on what dada did*:
///   - member  → its own center: real error copies, the within-cluster cloud.
///   - failed  → nearest center: shed by the abundance test but not assigned.
///   - center  → other centers: survivors vs each other (resolution floor).
///
/// `birth_type`/`birth_pval` ride along on the partner center so prior-born ASVs
/// (pseudo-pool) and births near OMEGA_A are visible in the same table.
fn from_dada_mode(
    w: &mut dyn Write,
    pool: &rayon::ThreadPool,
    samples: &[DadaSample],
    p: &Params,
) -> io::Result<()> {
    writeln!(
        w,
        "sample,class,cluster,ab,center_ab,ab_ratio,birth_type,birth_pval,\
         kdist,edits,core_len,pct_div,band_req,screened_in"
    )?;
    for s in samples {
        let ncenters = s.c_enc.len();
        // Build the job list from each unique's fate.
        let mut jobs: Vec<Job> = Vec::with_capacity(s.enc.len());
        let (mut n_center, mut n_member, mut n_failed) = (0u64, 0u64, 0u64);
        for (i, m) in s.map.iter().enumerate() {
            match m {
                None => {
                    // nearest center by k-mer distance (cheap pre-pass, no align)
                    if ncenters > 0 {
                        let c = (0..ncenters)
                            .min_by(|&a, &b| {
                                let ka = kmer_dist8(
                                    &s.kmers[i],
                                    s.enc[i].len(),
                                    &s.c_kmers[a],
                                    s.c_enc[a].len(),
                                    p.k,
                                );
                                let kb = kmer_dist8(
                                    &s.kmers[i],
                                    s.enc[i].len(),
                                    &s.c_kmers[b],
                                    s.c_enc[b].len(),
                                    p.k,
                                );
                                ka.partial_cmp(&kb).unwrap()
                            })
                            .unwrap();
                        jobs.push(Job::Failed { i, c });
                    }
                    n_failed += 1;
                }
                Some(c) => {
                    let c = *c;
                    if s.enc[i] == s.c_enc[c] {
                        n_center += 1; // the representative itself — no self-align
                    } else {
                        jobs.push(Job::Member { i, c });
                        n_member += 1;
                    }
                }
            }
        }
        // Inter-center pairs (resolution floor): enumerate / subsample.
        let cpairs = pairs_for(ncenters, p.max_pairs, p.seed ^ 0xC0FFEE);
        for (a, b) in cpairs {
            jobs.push(Job::CenterPair { a, b });
        }
        if p.verbose {
            eprintln!(
                "[kdist] {} : {} uniques ({n_center} centers, {n_member} members, {n_failed} failed), \
                 {ncenters} ASVs, {} jobs (k={}, band={}, {} threads)",
                s.name,
                s.enc.len(),
                jobs.len(),
                p.k,
                p.band,
                p.threads,
            );
        }
        // Align every job in parallel.
        let rows: Vec<(usize, f64, usize, usize, f64, usize)> = pool.install(|| {
            jobs.par_iter()
                .map_init(AlignBuffers::new, |buf, job| {
                    let (ei, ej, ki, kj, c) = match job {
                        Job::Member { i, c } => {
                            (&s.enc[*i], &s.c_enc[*c], &s.kmers[*i], &s.c_kmers[*c], *c)
                        }
                        Job::Failed { i, c } => {
                            (&s.enc[*i], &s.c_enc[*c], &s.kmers[*i], &s.c_kmers[*c], *c)
                        }
                        Job::CenterPair { a, b } => (
                            &s.c_enc[*a],
                            &s.c_enc[*b],
                            &s.c_kmers[*a],
                            &s.c_kmers[*b],
                            *b,
                        ),
                    };
                    let kd = kmer_dist8(ki, ei.len(), kj, ej.len(), p.k);
                    align_endsfree_with_buf(ei, ej, 5, -4, -8, p.band, buf);
                    let (edits, core, band_req) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (c, kd, edits, core, pct, band_req)
                })
                .collect()
        });
        // Emit, pairing each row back with its job for labels/abundances.
        // Track the failed class by abundance: singletons can never seed an ASV
        // under the default (≥2 reads, unless --detect-singletons), so a failed
        // singleton fails for that reason, NOT for being distant — split it out.
        let (mut f_singleton, mut f_singleton_in, mut f_multi, mut f_multi_in) =
            (0u64, 0u64, 0u64, 0u64);
        for (job, &(c, kd, edits, core, pct, band_req)) in jobs.iter().zip(&rows) {
            let (class, ab, center_ab) = match job {
                Job::Member { i, .. } => ("member", s.counts[*i], s.c_ab[c]),
                Job::Failed { i, .. } => ("failed", s.counts[*i], s.c_ab[c]),
                Job::CenterPair { a, .. } => ("center_pair", s.c_ab[*a], s.c_ab[c]),
            };
            if let Job::Failed { .. } = job {
                let within = (kd < p.cutoff) as u64;
                if ab <= 1 {
                    f_singleton += 1;
                    f_singleton_in += within;
                } else {
                    f_multi += 1;
                    f_multi_in += within;
                }
            }
            let ratio = center_ab as f64 / ab.max(1) as f64;
            writeln!(
                w,
                "{},{class},{c},{ab},{center_ab},{ratio:.2},{},{:.3e},{kd:.4},{edits},{core},{pct:.3},{band_req},{}",
                s.name,
                s.c_birth[c],
                s.c_birth_pval[c],
                (kd < p.cutoff) as u8,
            )?;
        }
        if p.verbose {
            let f_total = f_singleton + f_multi;
            if f_total > 0 {
                eprintln!(
                    "[kdist] {} : {f_total} failed | singletons {f_singleton} ({f_singleton_in} within cutoff) \
                     | multi-read {f_multi} ({f_multi_in} within cutoff) — failed singletons are the \
                     --detect-singletons tradeoff, not distance",
                    s.name,
                );
            }
            let priors = s.c_birth.iter().filter(|b| b.as_str() == "Prior").count();
            if priors > 0 {
                eprintln!(
                    "[kdist] {} : {priors}/{ncenters} ASVs born from priors (pseudo); \
                     filter the table on class=center_pair,birth_type=Prior to see their nearest survivor",
                    s.name,
                );
            }
        }
    }
    Ok(())
}

/// Candidate band sizes for the band-fit summary (DADA2 default is 16).
const BAND_SWEEP: [usize; 7] = [2, 4, 8, 16, 32, 64, 128];

/// For each candidate band B, the fraction of alignments whose true path fits
/// within B (band_req <= B) — i.e. that a banded aligner at B would compute
/// correctly. The complement is distorted/effectively-rejected by that band.
fn band_fit(label: &str, band_reqs: &[usize]) -> String {
    let n = band_reqs.len();
    if n == 0 {
        return format!("{label} band-fit: (none)");
    }
    let parts: Vec<String> = BAND_SWEEP
        .iter()
        .map(|&b| {
            let f = band_reqs.iter().filter(|&&r| r <= b).count();
            format!("≤{b}:{:.1}%", 100.0 * f as f64 / n as f64)
        })
        .collect();
    let mx = band_reqs.iter().copied().max().unwrap_or(0);
    format!("{label} band-fit ({n}, max_req {mx}): {}", parts.join(" "))
}

/// Divergence below which a nearest-parent link is treated as a clear
/// error-copy candidate when computing the screen "headroom".
const CLEAR_ERROR_COPY_PCT: f64 = 3.0;

/// Abundance-aware mode: for each unique, find its nearest MORE-abundant
/// neighbour (its candidate error-copy "parent", mirroring DADA2's greedy
/// center-based comparison) by k-mer distance, then align that one pair for the
/// true divergence. The distribution of parent-link kdist is the empirical
/// error-copy distance ceiling; `cutoff − ceiling` is the screen's headroom.
fn nearest_parent_mode(
    w: &mut dyn Write,
    pool: &rayon::ThreadPool,
    pops: &[Sample],
    p: &Params,
) -> io::Result<()> {
    writeln!(
        w,
        "sample,ab,parent_ab,ab_ratio,kdist,edits,core_len,pct_div,band_req,screened_in"
    )?;
    for s in pops {
        let n = s.enc.len();
        // abundance-descending order (stable by index for ties)
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| s.counts[b].cmp(&s.counts[a]).then(a.cmp(&b)));
        if p.verbose {
            eprintln!(
                "[kdist] {} : {n} uniques, nearest more-abundant parent each (k={}, band={}, {} threads)",
                s.name, p.k, p.band, p.threads,
            );
        }
        // For each non-top unique (position r), scan the more-abundant prefix
        // order[0..r] for the min-kdist parent, then align that pair.
        let rows: Vec<(usize, usize, f64, usize, usize, f64, usize)> = pool.install(|| {
            (1..n)
                .into_par_iter()
                .map_init(AlignBuffers::new, |buf, r| {
                    let i = order[r];
                    let (mut best_kd, mut parent) = (f64::INFINITY, order[0]);
                    for &c in &order[0..r] {
                        let kd = kmer_dist8(
                            &s.kmers[i],
                            s.enc[i].len(),
                            &s.kmers[c],
                            s.enc[c].len(),
                            p.k,
                        );
                        if kd < best_kd {
                            best_kd = kd;
                            parent = c;
                        }
                    }
                    align_endsfree_with_buf(&s.enc[i], &s.enc[parent], 5, -4, -8, p.band, buf);
                    let (edits, core, band_req) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (i, parent, best_kd, edits, core, pct, band_req)
                })
                .collect()
        });
        // Headroom: among clear error-copy links (<= CLEAR_ERROR_COPY_PCT
        // divergence) the max kdist is the ceiling the cutoff must cover.
        let (mut linked, mut total) = (0u64, 0u64);
        let mut ceiling = 0.0f64;
        let mut kds: Vec<f64> = Vec::with_capacity(rows.len());
        let mut band_ec: Vec<usize> = Vec::new(); // band_req of clear error-copy links
        for &(i, parent, kd, edits, core, pct, band_req) in &rows {
            total += 1;
            if kd < p.cutoff {
                linked += 1;
            }
            if pct <= CLEAR_ERROR_COPY_PCT {
                ceiling = ceiling.max(kd);
                band_ec.push(band_req);
            }
            kds.push(kd);
            let ratio = s.counts[parent] as f64 / s.counts[i].max(1) as f64;
            writeln!(
                w,
                "{},{},{},{ratio:.2},{kd:.4},{edits},{core},{pct:.3},{band_req},{}",
                s.name,
                s.counts[i],
                s.counts[parent],
                (kd < p.cutoff) as u8
            )?;
        }
        if p.verbose && total > 0 {
            kds.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p90 = kds[(kds.len() * 9 / 10).min(kds.len() - 1)];
            eprintln!(
                "[kdist] {} : {total} children | nearest-parent kdist median {:.3} p90 {:.3} | \
                 {linked} ({:.1}%) within cutoff {} | clear-error-copy ceiling {:.3} -> headroom {:.3}",
                s.name,
                kds[kds.len() / 2],
                p90,
                100.0 * linked as f64 / total as f64,
                p.cutoff,
                ceiling,
                p.cutoff - ceiling,
            );
            eprintln!(
                "[kdist] {} : {}",
                s.name,
                band_fit("clear-error-copy", &band_ec)
            );
        }
    }
    Ok(())
}
