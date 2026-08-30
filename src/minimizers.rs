//! Minimizer-based pre-alignment screening (experimental).
//!
//! An alternative to the `4^k` frequency-vector screen in [`crate::kmers`].
//! Both answer the same question — *is this pair similar enough to be worth
//! aligning?* — and both feed the same `kdist_cutoff` gate in
//! `raw_align_with_buf`. Neither is the cluster definition; see
//! `docs/findings/kmer-size-screening.md`.
//!
//! # Why minimizers
//!
//! The ESPRIT-era screen ([Sun et al. 2009]) compares full `4^k` count vectors.
//! That has two structural costs this module attacks:
//!
//! 1. **Memory is `4^k` per raw, independent of sequence length** — 16 KB at
//!    k=7, which is what forced the sparse representation in issue #43 and what
//!    gates raising the k cap (#44).
//! 2. **Specificity does not scale with read length.** At k=5 a 1.5 kb HiFi read
//!    has ~1,500 k-mers drawn from a 1,024-slot space, so *every* pair shares
//!    nearly everything and the screen is close to a no-op — the measured reason
//!    PacBio needs k=7 (`docs/findings/kmer-size-screening.md`).
//!
//! A winnowed minimizer sketch is `O(len / w)` entries regardless of `k`, so `k`
//! becomes free to raise for specificity while memory *falls*. That decouples
//! the two knobs the frequency vector welds together.
//!
//! # The metric is deliberately the same shape
//!
//! [`minimizer_dist`] computes
//!
//! ```text
//! dist = 1 - Σ min(count_a[m], count_b[m]) / min(|sketch_a|, |sketch_b|)
//! ```
//!
//! which is [`crate::kmers::kmer_dist8`]'s formula evaluated on a winnowed
//! subsample of the k-mer space rather than on all of it. This is a choice, not
//! a coincidence: it keeps `KDIST_CUTOFF` interpretable across both backends, so
//! an A/B changes *which k-mers are consulted* and nothing else about the gate.
//!
//! [Sun et al. 2009]: https://doi.org/10.1093/nar/gkp285

use std::collections::HashMap;

/// Default k-mer size for the minimizer sketch.
///
/// Deliberately much larger than the frequency screen's `KMER_SIZE = 5`. The
/// sketch's memory does not depend on `k`, so `k` is chosen purely for
/// discrimination: it must be large enough that unrelated amplicons share
/// almost no minimizers, and small enough that a single substitution — which
/// destroys the `k` k-mers spanning it — removes only a small fraction of the
/// sketch. At k=11 a substitution invalidates 11 of ~240 k-mer positions in a
/// 250 bp read, i.e. ~2-3 of ~48 minimizers.
pub const MINIMIZER_K: usize = 11;

/// Default winnowing window, in k-mers.
///
/// The sketch retains ~`2/(w+1)` of positions (Schleimer et al. 2003's density
/// bound for winnowing), so `w = 5` keeps roughly a third — dense enough that a
/// substitution cannot wipe out a whole read's shared minimizers, sparse enough
/// to be a real memory reduction.
pub const MINIMIZER_W: usize = 5;

/// Valid ranges. `k` is bounded above by 31 because two bits per base must fit
/// a `u64`; below by 7 because shorter k-mers are not discriminating enough to
/// be worth winnowing (use the frequency screen instead).
pub const MINIMIZER_K_MIN: usize = 7;
pub const MINIMIZER_K_MAX: usize = 31;
pub const MINIMIZER_W_MIN: usize = 1;
pub const MINIMIZER_W_MAX: usize = 64;

/// SplitMix64 finalizer.
///
/// Winnowing on the raw 2-bit packing would be pathologically biased: `A` packs
/// to `00`, so a poly-A run is the numerically smallest k-mer in any window it
/// touches and would be selected everywhere, concentrating sketches on exactly
/// the low-complexity content that carries no discriminating signal. Avalanching
/// the packed value first makes minimizer selection effectively uniform over
/// k-mer identity.
#[inline]
fn hash64(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

/// A winnowed minimizer sketch: sorted `(hash, count)` pairs.
///
/// Sorted-and-counted rather than a plain set so [`minimizer_dist`] can be the
/// same `Σ min(count)` merge-join as [`crate::kmers::kmer_dist8_sparse`], and so
/// a repeated minimizer (a real signal in repetitive amplicons) is not silently
/// collapsed. Counts saturate at 255, matching the 8-bit frequency screen.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MinimizerSketch {
    entries: Vec<(u64, u8)>,
    /// Σ counts — the multiset size, and the denominator in [`minimizer_dist`].
    total: u32,
}

impl MinimizerSketch {
    /// Number of *distinct* minimizers retained.
    #[inline]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Multiset size (Σ counts), the normaliser in [`minimizer_dist`].
    #[inline]
    pub fn total(&self) -> u32 {
        self.total
    }

    /// Resident byte footprint, for `--verbose` accounting alongside
    /// [`crate::kmers::KmerScreen::resident_bytes`].
    #[inline]
    pub fn resident_bytes(&self) -> usize {
        self.entries.len() * std::mem::size_of::<(u64, u8)>()
    }

    /// The distinct minimizer hashes, ascending.
    #[inline]
    pub fn hashes(&self) -> impl Iterator<Item = u64> + '_ {
        self.entries.iter().map(|&(h, _)| h)
    }

    #[inline]
    fn entries(&self) -> &[(u64, u8)] {
        &self.entries
    }
}

/// Build the winnowed minimizer sketch of an integer-encoded sequence
/// (A=1, C=2, G=3, T=4, as produced by `misc::nt_encode`).
///
/// Forward strand only. DADA2 compares reads that are already in a consistent
/// orientation, so canonical (strand-symmetric) minimizers would merge k-mers
/// the error model treats as distinct.
///
/// Non-ACGT positions break the k-mer run: any window containing one is skipped,
/// mirroring the frequency screen's silent rejection of invalid k-mers.
///
/// Windows are winnowed with the standard monotonic deque, taking the **leftmost**
/// minimum on ties (Schleimer et al. 2003). Leftmost-on-ties is what makes the
/// sketch stable under local edits: it guarantees adjacent windows agree on the
/// same occurrence when one exists, so a substitution perturbs only the windows
/// it actually spans.
pub fn sketch(seq: &[u8], k: usize, w: usize) -> MinimizerSketch {
    debug_assert!((MINIMIZER_K_MIN..=MINIMIZER_K_MAX).contains(&k));
    debug_assert!((MINIMIZER_W_MIN..=MINIMIZER_W_MAX).contains(&w));

    // (hash, position) of each k-mer, with invalid k-mers omitted. Positions are
    // kept so the deque can evict by window, and so ties resolve leftmost.
    let mut selected: Vec<u64> = Vec::new();
    // Monotonic deque of (hash, pos), increasing in hash.
    let mut deque: std::collections::VecDeque<(u64, usize)> = std::collections::VecDeque::new();

    let mut packed: u64 = 0;
    let mut valid_run = 0usize; // consecutive valid bases ending at the current one
    let mask: u64 = if k >= 32 { !0 } else { (1u64 << (2 * k)) - 1 };
    let mut last_emitted: Option<u64> = None;

    for (pos, &nt) in seq.iter().enumerate() {
        match nt {
            1..=4 => {
                packed = ((packed << 2) | (nt as u64 - 1)) & mask;
                valid_run += 1;
            }
            _ => {
                // Invalid base: the k-mer run restarts, and so does winnowing —
                // a window straddling the break has no meaningful minimum.
                valid_run = 0;
                deque.clear();
                last_emitted = None;
                continue;
            }
        }
        if valid_run < k {
            continue;
        }
        let kmer_pos = pos + 1 - k; // start offset of this k-mer
        let h = hash64(packed);

        // Evict k-mers that have fallen out of the window.
        while let Some(&(_, p)) = deque.front() {
            if p + w <= kmer_pos {
                deque.pop_front();
            } else {
                break;
            }
        }
        // Maintain increasing-hash order; `<` (not `<=`) keeps the earlier
        // occurrence on a tie, which is the leftmost-minimum rule.
        while let Some(&(hb, _)) = deque.back() {
            if hb > h {
                deque.pop_back();
            } else {
                break;
            }
        }
        deque.push_back((h, kmer_pos));

        // A full window has been seen once w k-mers are in flight.
        if kmer_pos + 1 >= w {
            let &(mh, _) = deque.front().expect("deque non-empty after push");
            // Emit only when the window's minimum *hash* changes.
            //
            // Deduplicating on hash rather than on (hash, position) is a
            // deliberate divergence from the mapping-oriented convention
            // (minimap2 et al. key on position, because they need to know
            // *where* a seed matched). We only need to know *whether* content is
            // shared, and per-occurrence emission actively hurts that: inside a
            // homopolymer every window has the same minimal k-mer at a new
            // offset, so the run emits at density 1.0 against winnowing's
            // ~2/(w+1). That inflates its count, and `Σ min(count)` would then
            // read two sequences as similar merely because both contain a
            // homopolymer — a false-similarity mode on exactly the
            // low-complexity content that carries no discriminating signal.
            // Only *consecutive* duplicates collapse, so a genuine repeat at two
            // separated loci is still counted twice.
            if last_emitted != Some(mh) {
                selected.push(mh);
                last_emitted = Some(mh);
            }
        }
    }

    build_entries(selected)
}

/// Collapse selected hashes into the sorted `(hash, count)` representation.
fn build_entries(mut selected: Vec<u64>) -> MinimizerSketch {
    selected.sort_unstable();
    let mut entries: Vec<(u64, u8)> = Vec::with_capacity(selected.len());
    let mut total: u32 = 0;
    for h in selected {
        match entries.last_mut() {
            Some((last_h, c)) if *last_h == h => {
                *c = c.saturating_add(1);
            }
            _ => entries.push((h, 1)),
        }
        total = total.saturating_add(1);
    }
    // `total` must agree with Σ counts even after saturation, or the distance
    // denominator and numerator would be on different scales.
    let summed: u32 = entries.iter().map(|&(_, c)| c as u32).sum();
    MinimizerSketch {
        entries,
        total: summed.min(total),
    }
}

/// Minimizer-sketch distance, the drop-in counterpart of
/// [`crate::kmers::kmer_dist8`].
///
/// `dist = 1 - Σ min(count_a, count_b) / min(total_a, total_b)`, so the result
/// is in `[0, 1]`: 0 for identical sketches, 1 when nothing is shared. An empty
/// sketch on either side yields 1.0 (nothing is known to be shared), which is
/// above any usable cutoff and therefore screens the pair out — the safe
/// direction is *not* taken here, see [`sketch_is_usable`].
pub fn minimizer_dist(a: &MinimizerSketch, b: &MinimizerSketch) -> f64 {
    let scale = a.total.min(b.total);
    if scale == 0 {
        return 1.0;
    }
    let (ea, eb) = (a.entries(), b.entries());
    let mut shared: u32 = 0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < ea.len() && j < eb.len() {
        let (ha, ca) = ea[i];
        let (hb, cb) = eb[j];
        match ha.cmp(&hb) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += ca.min(cb) as u32;
                i += 1;
                j += 1;
            }
        }
    }
    1.0 - shared as f64 / scale as f64
}

/// Whether a sequence is long enough for its sketch to mean anything.
///
/// A sequence shorter than `k + w - 1` produces no full window and therefore an
/// empty sketch, which [`minimizer_dist`] reports as maximally distant. Screening
/// such a pair *out* would be a silent false negative — it would drop the pair
/// from consideration entirely on the basis of no evidence. Callers must fall
/// back to aligning (i.e. treat the screen as passing) when this is false.
#[inline]
pub fn sketch_is_usable(seq_len: usize, k: usize, w: usize) -> bool {
    seq_len >= k + w - 1
}

// ---------------------------------------------------------------------------
// Inverted index
// ---------------------------------------------------------------------------

/// Minimizer → raw-index posting lists over a fixed set of `Raw`s.
///
/// `b_compare` is called once per cluster and scans **all** raws against that
/// cluster's center, so the raws are the stable side of the comparison and the
/// center is the query. Building the index once per `dada` run and querying it
/// with each center replaces `nraw` independent merge-joins with a single
/// scatter over the center's posting lists.
///
/// This is an **exact** acceleration of [`minimizer_dist`], not an
/// approximation: [`MinimizerIndex::shared_counts`] recovers the same
/// `Σ min(count_a, count_b)` numerator that the pairwise merge-join computes,
/// so index-on and index-off produce bit-identical distances. Only
/// minimizer-vs-k-mer is a behavioural change; keeping those two axes separable
/// is the point.
pub struct MinimizerIndex {
    postings: HashMap<u64, Vec<(u32, u8)>>,
    nraw: usize,
}

impl MinimizerIndex {
    /// Build over the sketches of every raw, in `Raw::index` order.
    ///
    /// Raws with no usable sketch contribute no postings; they are recovered by
    /// [`MinimizerIndex::query_into`]'s unusable-sketch fallback.
    pub fn build(sketches: &[Option<MinimizerSketch>]) -> Self {
        // One posting entry per (raw, distinct minimizer) incidence.
        let incidences: usize = sketches
            .iter()
            .map(|s| s.as_ref().map_or(0, |s| s.len()))
            .sum();
        let mut postings: HashMap<u64, Vec<(u32, u8)>> =
            HashMap::with_capacity(incidences / 8 + 16);
        for (idx, s) in sketches.iter().enumerate() {
            let Some(s) = s else { continue };
            for &(h, c) in s.entries() {
                postings.entry(h).or_default().push((idx as u32, c));
            }
        }
        MinimizerIndex {
            postings,
            nraw: sketches.len(),
        }
    }

    /// Number of `(raw, minimizer)` incidences stored.
    pub fn n_postings(&self) -> usize {
        self.postings.values().map(|v| v.len()).sum()
    }

    /// Number of distinct minimizers in the pool.
    pub fn n_keys(&self) -> usize {
        self.postings.len()
    }

    /// Resident byte footprint (postings + key overhead), for `--verbose`.
    pub fn resident_bytes(&self) -> usize {
        self.n_postings() * std::mem::size_of::<(u32, u8)>()
            + self.n_keys() * (std::mem::size_of::<u64>() + std::mem::size_of::<Vec<(u32, u8)>>())
    }

    /// Scatter `query`'s minimizers over the postings, writing each raw's shared
    /// count into `out[raw_index]`.
    ///
    /// `out` is resized to `nraw` and fully cleared, so the caller may reuse one
    /// buffer across clusters. Clearing is `O(nraw)` sequential — cheaper than
    /// tracking touched indices, and `b_compare` walks all `nraw` regardless.
    pub fn shared_counts(&self, query: &MinimizerSketch, out: &mut Vec<u32>) {
        out.clear();
        out.resize(self.nraw, 0);
        for &(h, qc) in query.entries() {
            let Some(list) = self.postings.get(&h) else {
                continue;
            };
            for &(idx, c) in list {
                // Σ min(count_query, count_raw) — the same numerator the
                // pairwise merge-join accumulates.
                out[idx as usize] += qc.min(c) as u32;
            }
        }
    }
}

/// Distance from a precomputed shared count, matching [`minimizer_dist`].
///
/// The pair of `(shared, query_total, raw_total)` is exactly what
/// [`MinimizerIndex::shared_counts`] plus the two sketches supply, so this is
/// the index-side evaluation of the same formula.
#[inline]
pub fn dist_from_shared(shared: u32, total_a: u32, total_b: u32) -> f64 {
    let scale = total_a.min(total_b);
    if scale == 0 {
        return 1.0;
    }
    1.0 - shared as f64 / scale as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode(s: &str) -> Vec<u8> {
        s.bytes()
            .map(|b| match b {
                b'A' => 1,
                b'C' => 2,
                b'G' => 3,
                b'T' => 4,
                _ => 5,
            })
            .collect()
    }

    fn make_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((s >> 33) % 4) as u8 + 1
            })
            .collect()
    }

    #[test]
    fn identical_sequences_have_zero_distance() {
        for seed in 0..8 {
            let seq = make_seq(250, seed);
            let s = sketch(&seq, MINIMIZER_K, MINIMIZER_W);
            assert!(!s.is_empty(), "seed {seed} produced an empty sketch");
            assert_eq!(minimizer_dist(&s, &s), 0.0);
        }
    }

    #[test]
    fn distance_is_bounded_and_symmetric() {
        for seed in 0..16 {
            let a = sketch(&make_seq(250, seed), MINIMIZER_K, MINIMIZER_W);
            let b = sketch(&make_seq(250, seed + 100), MINIMIZER_K, MINIMIZER_W);
            let d = minimizer_dist(&a, &b);
            assert!((0.0..=1.0).contains(&d), "d={d} out of range");
            assert_eq!(d, minimizer_dist(&b, &a));
        }
    }

    /// The screen's whole purpose: a one-substitution error copy must land far
    /// below the cutoff, while an unrelated sequence lands far above it.
    #[test]
    fn separates_error_copies_from_unrelated() {
        const CUTOFF: f64 = 0.42;
        for seed in 0..16 {
            let seq = make_seq(250, seed);
            let mut mutated = seq.clone();
            mutated[125] = (mutated[125] % 4) + 1; // single substitution
            let (a, m) = (
                sketch(&seq, MINIMIZER_K, MINIMIZER_W),
                sketch(&mutated, MINIMIZER_K, MINIMIZER_W),
            );
            let d_err = minimizer_dist(&a, &m);
            assert!(
                d_err < CUTOFF,
                "seed {seed}: 1-sub error copy screened OUT at d={d_err}"
            );

            let u = sketch(&make_seq(250, seed + 1000), MINIMIZER_K, MINIMIZER_W);
            let d_unrel = minimizer_dist(&a, &u);
            assert!(
                d_unrel > CUTOFF,
                "seed {seed}: unrelated pair passed at d={d_unrel}"
            );
        }
    }

    /// Poly-A must not dominate selection — the reason for the hash avalanche.
    #[test]
    fn homopolymer_does_not_capture_every_window() {
        // 40 A's embedded in random sequence; without hashing, the all-A k-mer
        // would be the minimum of every window it touches.
        let mut seq = make_seq(250, 7);
        seq[100..140].fill(1);
        let s = sketch(&seq, MINIMIZER_K, MINIMIZER_W);
        let poly_a = hash64(0); // 2-bit packing of AAAA... is all zero bits
        let poly_count = s
            .entries()
            .iter()
            .find(|&&(h, _)| h == poly_a)
            .map_or(0, |&(_, c)| c as usize);
        assert!(
            poly_count * 4 < s.total() as usize,
            "poly-A took {poly_count} of {} sketch mass",
            s.total()
        );
    }

    #[test]
    fn invalid_bases_break_the_run() {
        let mut seq = make_seq(120, 3);
        let clean = sketch(&seq, MINIMIZER_K, MINIMIZER_W);
        seq[60] = 5; // N
        let broken = sketch(&seq, MINIMIZER_K, MINIMIZER_W);
        // The N invalidates the windows spanning it, so the sketch must shrink
        // but not vanish.
        assert!(broken.total() < clean.total());
        assert!(broken.total() > 0);
    }

    #[test]
    fn short_sequences_are_flagged_unusable() {
        let k = MINIMIZER_K;
        let w = MINIMIZER_W;
        assert!(!sketch_is_usable(k + w - 2, k, w));
        assert!(sketch_is_usable(k + w - 1, k, w));
        let short = make_seq(k + w - 2, 1);
        assert!(sketch(&short, k, w).is_empty());
    }

    /// The index must reproduce the pairwise merge-join exactly — this is what
    /// lets the findings doc attribute any output change to the *metric* rather
    /// than to the acceleration.
    #[test]
    fn index_reproduces_pairwise_distance_exactly() {
        let seqs: Vec<Vec<u8>> = (0..64u64)
            .map(|s| make_seq(200 + (s as usize % 40), s))
            .collect();
        let sketches: Vec<Option<MinimizerSketch>> = seqs
            .iter()
            .map(|s| Some(sketch(s, MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let index = MinimizerIndex::build(&sketches);
        let mut counts = Vec::new();
        for q in &sketches {
            let q = q.as_ref().unwrap();
            index.shared_counts(q, &mut counts);
            for (i, s) in sketches.iter().enumerate() {
                let s = s.as_ref().unwrap();
                let via_index = dist_from_shared(counts[i], q.total(), s.total());
                let pairwise = minimizer_dist(q, s);
                assert_eq!(
                    via_index.to_bits(),
                    pairwise.to_bits(),
                    "index/pairwise disagree at raw {i}"
                );
            }
        }
    }

    #[test]
    fn sketch_is_smaller_than_the_dense_frequency_screen() {
        // The memory claim, asserted rather than left to prose: at k=7 the dense
        // frequency screen is 16 KB/raw regardless of length.
        let seq = make_seq(250, 42);
        let s = sketch(&seq, MINIMIZER_K, MINIMIZER_W);
        assert!(s.resident_bytes() < crate::kmers::n_kmers(7));
    }

    #[test]
    fn sketch_handles_all_valid_encodings() {
        let seq = encode("ACGTACGTACGTACGTACGTACGTACGTACGT");
        let s = sketch(&seq, MINIMIZER_K, MINIMIZER_W);
        assert!(!s.is_empty());
        assert_eq!(minimizer_dist(&s, &s), 0.0);
    }
}
