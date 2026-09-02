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
//! # The metric shares a formula, but NOT an operating point
//!
//! [`minimizer_dist`] computes
//!
//! ```text
//! dist = 1 - Σ min(count_a[m], count_b[m]) / min(total_a, total_b)
//! ```
//!
//! which is [`crate::kmers::kmer_dist8`]'s formula evaluated on a winnowed
//! subsample of the k-mer space rather than on all of it.
//!
//! **The shared formula does not make `KDIST_CUTOFF` transfer, and assuming it
//! did was a real error on this branch.** A cutoff is a property of the distance
//! *distribution* a metric induces, not of its algebra. Measured on a MiSeq SOP
//! sample, at the same `0.42`:
//!
//! | screen | pairs passed |
//! |---|---|
//! | k-mer, k=5 | 27.6% |
//! | minimizer, k=8 | 9.0% |
//!
//! The matching operating point is near **0.64**. Running the minimizer screen
//! at the k-mer screen's cutoff over-screens by ~3x, which starves raws of any
//! cluster match and fragments them into spurious low-abundance ASVs. See
//! `docs/findings/minimizer-screening.md`; a per-backend calibrated cutoff is a
//! prerequisite for promoting this backend.
//!
//! [Sun et al. 2009]: https://doi.org/10.1093/nar/gkp285

use std::collections::HashMap;

/// Default k-mer size for the minimizer sketch.
///
/// Larger than the frequency screen's `KMER_SIZE = 5`, because the sketch's
/// memory does not depend on `k` — so `k` is free to be chosen for
/// discrimination alone.
///
/// **11 was a back-of-envelope choice and measurement contradicted it.** The
/// argument was that a substitution destroys only the `k` k-mers spanning it,
/// so ~2-3 of ~48 minimizers in a 250 bp read. True, but it ignored the
/// aggregate pass rate: at k=11 the screen shrouds pairs the k-mer screen
/// aligns at 6-10 substitutions (14-26 per sample on the MiSeq SOP), and k=8
/// drives that to zero everywhere tested while still passing only ~9% of pairs.
///
/// Now **8**, the value the sweeps settled on. 11 was a back-of-envelope choice:
/// a substitution destroys only the `k` k-mers spanning it, so ~2-3 of ~48
/// minimizers in a 250 bp read. True, but it ignored the aggregate pass rate --
/// at k=11 the screen shrouds pairs the k-mer screen aligns at 6-10 substitutions
/// (14-26 per sample on the MiSeq SOP), and k=8 drives that to zero everywhere
/// tested. k=8 also beat k=9 on ASV set agreement at every cutoff on the
/// 362-sample MiSeq run.
pub const MINIMIZER_K: usize = 8;

/// Default screen cutoff for the minimizer backend, replacing `KDIST_CUTOFF`'s
/// 0.42 when the cutoff is not given explicitly.
///
/// A cutoff is a property of the distance *distribution* a metric induces, not of
/// its algebra, so 0.42 does not carry over: it passes ~28% of pairs on the
/// frequency vector and ~9% on the sketch, over-screening ~3x and fragmenting
/// clusters. Sweeps put the useful region at **0.62-0.65** on Illumina, with a
/// sharp count-L1 minimum at 0.63 on pooled ITS2 and the read-retention crossing
/// at ~0.636.
///
/// This is a *starting point*, not a substitute for calibration. The right value
/// tracks the k-mer screen's pass rate, which spans 0.70%-26.8% across measured
/// workloads: PacBio HiFi wants ~0.50 and low-diversity MiSeq ~0.70-0.80. Derive
/// it with `kdist-calibrate --screen-backend minimizer` and confirm by sweep; see
/// `docs/findings/minimizer-screening.md`.
pub const MINIMIZER_KDIST_CUTOFF: f64 = 0.63;

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
    /// Distinct minimizer hashes, ascending.
    ///
    /// Struct-of-arrays rather than `Vec<(u64, u8)>`: that tuple aligns to 8 and
    /// so occupies **16** bytes per entry, of which 7 are padding — 44% of a
    /// structure whose entire justification is being smaller than a `4^k` array.
    /// Split, an entry costs 9 bytes. The merge-join reads `hashes` densely and
    /// touches `counts` only on a match, which is the better access pattern
    /// anyway.
    hashes: Vec<u64>,
    counts: Vec<u8>,
    /// Σ counts — the multiset size, and the denominator in [`minimizer_dist`].
    total: u32,
}

impl MinimizerSketch {
    /// Number of *distinct* minimizers retained.
    #[inline]
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
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
        self.hashes.len() * std::mem::size_of::<u64>() + self.counts.len()
    }

    /// The distinct minimizer hashes, ascending.
    #[inline]
    pub fn hashes(&self) -> impl Iterator<Item = u64> + '_ {
        self.hashes.iter().copied()
    }

    /// `(hash, count)` pairs, ascending by hash.
    #[inline]
    fn entries(&self) -> impl Iterator<Item = (u64, u8)> + '_ {
        self.hashes.iter().copied().zip(self.counts.iter().copied())
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
    let mut hashes: Vec<u64> = Vec::with_capacity(selected.len());
    let mut counts: Vec<u8> = Vec::with_capacity(selected.len());
    for h in selected {
        match (hashes.last(), counts.last_mut()) {
            (Some(&last_h), Some(c)) if last_h == h => *c = c.saturating_add(1),
            _ => {
                hashes.push(h);
                counts.push(1);
            }
        }
    }
    hashes.shrink_to_fit();
    counts.shrink_to_fit();
    // Σ counts *after* saturation, so the distance's numerator and denominator
    // are on the same scale even in the (amplicon-implausible) saturating case.
    let total: u32 = counts.iter().map(|&c| c as u32).sum();
    MinimizerSketch {
        hashes,
        counts,
        total,
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
    dist_from_shared(shared_count(a, b), a.total, b.total)
}

/// `Σ min(count_a[m], count_b[m])` over the two sketches — the numerator of
/// [`minimizer_dist`], and the same quantity [`MinimizerIndex::shared_counts`]
/// accumulates by scatter. Factored out so the pairwise and index paths compute
/// one definition of "shared".
pub fn shared_count(a: &MinimizerSketch, b: &MinimizerSketch) -> u32 {
    let mut shared: u32 = 0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.hashes.len() && j < b.hashes.len() {
        match a.hashes[i].cmp(&b.hashes[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += a.counts[i].min(b.counts[j]) as u32;
                i += 1;
                j += 1;
            }
        }
    }
    shared
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

/// Default ceiling on [`IndexDecision::score`], above which the inverted index
/// costs more than it saves and the per-pair merge-join is used instead.
///
/// Three measured configurations bracket the crossover; this is not a calibrated
/// constant, and `DADA2RS_MINIMIZER_INDEX_MAX` overrides it.
///
/// | configuration | score | index verdict |
/// |---|---|---|
/// | pooled ITS2 | 0.073 | wins by 13.9% |
/// | per-sample PacBio | 0.185 | wins by 4.4% |
/// | pooled PacBio | 0.571 | **loses by 31.6%** |
///
/// The default sits between the last win and the first loss, biased toward the
/// low side because **the two errors are not symmetric**: choosing the index
/// wrongly cost +26% against the k-mer screen on pooled PacBio, while declining
/// it wrongly costs 2-4% on the configurations where it wins.
pub const MINIMIZER_INDEX_MAX_SCORE: f64 = 0.30;

/// Whether to build the inverted index for a given pool, and the evidence for it.
#[derive(Debug, Clone, Copy)]
pub struct IndexDecision {
    pub use_index: bool,
    /// Σ sketch entries over all raws — equals `n_postings` if the index is built.
    pub entries: usize,
    /// Distinct minimizers in the pool — equals `n_keys` if the index is built.
    pub distinct: usize,
    /// Mean posting-list length, `entries / distinct`.
    pub sharing: f64,
    /// `sharing × threads / nraw`. See [`decide_index`].
    pub score: f64,
    pub threshold: f64,
    /// `Some(true|false)` when `DADA2RS_MINIMIZER_INDEX` forced the outcome.
    pub forced: Option<bool>,
}

/// `HashSet<u64>` over values that are *already* avalanched.
///
/// The sketch stores SplitMix64 output, so re-hashing it with SipHash buys
/// nothing and costs ~20 ns per probe — which matters at 262 M inserts on a
/// pooled PacBio run. Only `write_u64` is reachable because the key is a `u64`.
#[derive(Default, Clone, Copy)]
struct IdentityHasher(u64);

impl std::hash::Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!("IdentityHasher is only used for u64 keys")
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = v;
    }
}

#[derive(Default, Clone, Copy)]
struct IdentityBuildHasher;

impl std::hash::BuildHasher for IdentityBuildHasher {
    type Hasher = IdentityHasher;
    fn build_hasher(&self) -> IdentityHasher {
        IdentityHasher(0)
    }
}

/// Should this pool use the inverted index, or the per-pair merge-join?
///
/// # Why there is a choice at all
///
/// The index makes the screen `O(1)` per pair, but pays for it with a **serial**
/// per-cluster scatter over the postings. Whether that trade is good is not a
/// property of the backend — it is a property of the workload, and it reverses:
///
/// | configuration | `setup` | map saved | setup/saved | verdict |
/// |---|---|---|---|---|
/// | pooled ITS2 | 16.14s | 56.14s | 0.29 | index wins 13.9% |
/// | per-sample PacBio | 33.44s | 76.26s | 0.44 | index wins 4.4% |
/// | pooled PacBio | **162.69s** | 43.72s | **3.72** | index **loses** 31.6% |
///
/// On pooled PacBio the indexed arm is slower than the k-mer screen it replaces.
/// Leaving that to an env var means the right answer depends on a flag no user
/// would know to set, so it is decided here.
///
/// # The score
///
/// Serial scatter work is `nclusters × entries × sharing`; the parallel saving it
/// buys is `ncomps × entries × c / threads`, and `ncomps ≈ nclusters × nraw`.
/// The ratio drops `nclusters` and `entries`, leaving
///
/// ```text
/// score = sharing × threads / nraw          where sharing = entries / distinct
/// ```
///
/// which ranks the three configurations above in the correct order (0.073,
/// 0.185, 0.571) against their measured ratios (0.29, 0.44, 3.72). **Three
/// points determine a bracket, not a threshold** — the crossover lies somewhere
/// in 0.19..0.57 and [`MINIMIZER_INDEX_MAX_SCORE`] picks a biased-low value
/// inside it.
///
/// `threads` must be the threads that will actually run *this pool's* compare
/// map — `rayon::current_num_threads()` inside the sub-pool, which is 4 under
/// per-sample concurrency and 48 for a single pooled run. Passing the global
/// count for a per-sample run inflates the score ~12x and would decline the
/// index exactly where it wins.
pub fn decide_index(
    sketches: &[Option<MinimizerSketch>],
    threads: usize,
    forced: Option<bool>,
) -> IndexDecision {
    let nraw = sketches.len();
    let entries: usize = sketches
        .iter()
        .map(|s| s.as_ref().map_or(0, |s| s.len()))
        .sum();

    // Distinct minimizers, counted without building the postings: a set of u64
    // is ~40k entries where the index would allocate 262M postings, so a
    // declined index costs a hash pass rather than a gigabyte.
    let mut keys: std::collections::HashSet<u64, IdentityBuildHasher> =
        std::collections::HashSet::with_capacity_and_hasher(1024, IdentityBuildHasher);
    for s in sketches.iter().flatten() {
        for h in s.hashes() {
            keys.insert(h);
        }
    }
    let distinct = keys.len();

    let threshold = std::env::var("DADA2RS_MINIMIZER_INDEX_MAX")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite() && *v > 0.0)
        .unwrap_or(MINIMIZER_INDEX_MAX_SCORE);

    let sharing = if distinct == 0 {
        0.0
    } else {
        entries as f64 / distinct as f64
    };
    let score = if nraw == 0 {
        0.0
    } else {
        sharing * threads.max(1) as f64 / nraw as f64
    };

    IndexDecision {
        use_index: forced.unwrap_or(score <= threshold),
        entries,
        distinct,
        sharing,
        score,
        threshold,
        forced,
    }
}

/// `DADA2RS_MINIMIZER_INDEX` as an **explicit override**, not a fallback: when it
/// names a mode, [`decide_index`] still computes its score for the log but the
/// override decides. Unset (or `auto`) leaves the choice to the rule.
///
/// | value | effect |
/// |---|---|
/// | `0`, `false`, `off`, `no` | force the per-pair merge-join |
/// | `1`, `true`, `on`, `yes` | force the inverted index |
/// | `auto`, unset | let [`decide_index`] choose |
///
/// Anything else **warns and falls back to auto**. It used to be that any value
/// other than `0` meant "on", so `DADA2RS_MINIMIZER_INDEX=true` silently changed
/// meaning when the rule landed — and on this branch a silently-ineffective A/B
/// arm has produced retracted results more than once. An unrecognised value is
/// far more likely to be a typo in a sweep script than an intent to use the
/// default, so it says so on stderr rather than quietly doing something else.
pub fn index_env_override() -> Option<bool> {
    let raw = std::env::var("DADA2RS_MINIMIZER_INDEX").ok()?;
    match raw.trim().to_ascii_lowercase().as_str() {
        "0" | "false" | "off" | "no" => Some(false),
        "1" | "true" | "on" | "yes" => Some(true),
        "" | "auto" => None,
        other => {
            eprintln!(
                "[dada] warning: DADA2RS_MINIMIZER_INDEX={other:?} is not recognised \
                 (expected 0/1, false/true, off/on, no/yes, or auto); \
                 leaving the index choice to the selection rule"
            );
            None
        }
    }
}

/// Result of deriving a minimizer cutoff by matching the k-mer screen's pass rate.
#[derive(Debug, Clone, Copy)]
pub struct MatchedPass {
    pub n_pairs: usize,
    /// Fraction of sampled pairs the k-mer screen passes at its own cutoff.
    pub kmer_pass: f64,
    /// The derived minimizer cutoff (raw quantile).
    pub cutoff: f64,
    /// The same, rounded to the 0.01 grid every sweep on this branch used.
    pub cutoff_rounded: f64,
    /// Minimizer pass rate at `cutoff_rounded` — should sit near `kmer_pass`;
    /// it will not match exactly, because rounding moves it.
    pub mini_pass: f64,
}

/// Derive a minimizer cutoff that passes the same fraction of pairs as the k-mer
/// screen does at `kmer_cutoff`.
///
/// # Why matched pass rate, and why it is cheap
///
/// The right cutoff has to be derived per dataset, because the k-mer screen's own
/// pass rate spans **0.70%-26.8%** across the workloads measured here — a fixed
/// minimizer cutoff cannot track that. Matching the pass rate predicted the sweep
/// optimum on 5 of 5 Illumina workloads.
///
/// Stated as a quantile it costs nothing. If the k-mer screen passes `q` of
/// pairs, the matching minimizer cutoff is simply **the `q`-quantile of the
/// minimizer distance distribution** — so this needs the two distance
/// distributions over a sample of pairs and *no alignment at all*. Alignment is
/// what makes `kdist-calibrate` expensive (every sampled pair aligned unbanded,
/// ~2.2 M DP cells on 1.5 kb reads), and it is only needed for the true-divergence
/// columns, which this rule never consults.
///
/// # What it does not do
///
/// It reproduces the k-mer screen's *selectivity*, which is the safe target, not
/// the cheapest cutoff that still agrees with it. On PacBio HiFi the ASV table is
/// **identical from 0.45 to 0.60**, and matched-pass picks 0.50 — inside that
/// plateau but at its expensive end, worth ~19 points of wall clock against 0.45.
/// Finding the cheap edge of a plateau needs the ASV table, so a sweep still
/// beats this; it is a default, not an answer.
pub fn matched_pass_cutoff(kmer_d: &[f64], mini_d: &[f64], kmer_cutoff: f64) -> MatchedPass {
    let n = kmer_d.len().min(mini_d.len());
    if n == 0 {
        return MatchedPass {
            n_pairs: 0,
            kmer_pass: 0.0,
            cutoff: MINIMIZER_KDIST_CUTOFF,
            cutoff_rounded: MINIMIZER_KDIST_CUTOFF,
            mini_pass: 0.0,
        };
    }
    let passed = kmer_d[..n].iter().filter(|&&d| d <= kmer_cutoff).count();
    let kmer_pass = passed as f64 / n as f64;

    let mut sorted: Vec<f64> = mini_d[..n].to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // The q-quantile, clamped so a pass rate of 0 or 1 still indexes in range.
    let idx = ((kmer_pass * n as f64).ceil() as usize).clamp(1, n) - 1;
    let cutoff = sorted[idx];
    let cutoff_rounded = (cutoff * 100.0).round() / 100.0;
    let mini_pass = sorted.iter().filter(|&&d| d <= cutoff_rounded).count() as f64 / n as f64;

    MatchedPass {
        n_pairs: n,
        kmer_pass,
        cutoff,
        cutoff_rounded,
        mini_pass,
    }
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
            for (h, c) in s.entries() {
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
    /// `counts` is owned by the caller (`B`) and reused across clusters, grown to
    /// `nraw` on first use and zeroed each call.
    ///
    /// # Where this phase's time actually goes
    ///
    /// This runs in `b_compare_parallel`'s **serial** setup, and on pooled ITS2 it
    /// measured 23.7% of the compare phase. Two things it is *not*:
    ///
    /// - **Not the zeroing.** A 3.3 MB memset is ~165 us, so ~0.6 s across 3,414
    ///   clusters. A variant that tracked touched indices to skip it added a
    ///   branch and a push per increment below and measured **20% slower** setup
    ///   on PacBio (5.97 -> 7.70 s). Removed.
    /// - **Not measurably the allocation either**, at least at PacBio's scale.
    ///   Reusing the buffer instead of allocating per cluster is 6.06 vs 5.97 s
    ///   there — inside noise. The reuse is kept because at pooled-ITS2's 3.3 MB
    ///   the allocation crosses glibc's mmap threshold (an mmap plus ~825
    ///   first-touch page faults plus munmap per cluster) whereas PacBio's 300 KB
    ///   does not, so the two workloads are not expected to behave alike. **That
    ///   is a mechanism argument, not a measurement** — it is untested at 3.3 MB.
    ///
    /// What remains is **the scatter itself**, which is the bulk of the phase on
    /// both workloads measured: ~9 s of 15.2 on pooled ITS2, ~5.7 s of 6.0 on
    /// PacBio. It is `O(Σ|posting(h)|)` over the centre's minimizers and
    /// single-threaded, so it is the Amdahl term worth attacking next — by
    /// parallelising it, or by asking whether the index earns its keep at high
    /// thread counts at all, since the pairwise merge-join it replaces is `O(sketch)`
    /// per pair but fully parallel (`DADA2RS_MINIMIZER_INDEX=0` measures that).
    /// See `docs/findings/minimizer-screening.md`.
    pub fn shared_counts(&self, query: &MinimizerSketch, counts: &mut Vec<u32>) {
        if counts.len() < self.nraw {
            counts.resize(self.nraw, 0);
        }
        counts[..self.nraw].fill(0);
        for (h, qc) in query.entries() {
            let Some(list) = self.postings.get(&h) else {
                continue;
            };
            for &(idx, c) in list {
                // Σ min(count_query, count_raw) — the same numerator the
                // pairwise merge-join accumulates.
                counts[idx as usize] += qc.min(c) as u32;
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

/// The screen verdict for a pair, with the fail-open rule applied.
///
/// **This is the single definition of what the minimizer screen decides**, used
/// by both the pairwise path and the index path so the two cannot drift. It
/// returns `0.0` — "as close as possible", i.e. *screen passes* — whenever there
/// is no evidence to judge on: either sketch missing, or either sketch empty.
///
/// Emptiness rather than sequence length is the right test. A sequence can be
/// long enough for a window (`sketch_is_usable`) and still sketch to nothing if
/// it is dense with non-ACGT, and in that case [`minimizer_dist`] would report
/// 1.0 — maximal distance — and screen the pair out on the basis of no evidence
/// at all. Failing open costs an alignment; failing closed silently loses ASVs.
///
/// `shared` is ignored on the fail-open paths, so the index may pass whatever it
/// accumulated.
#[inline]
pub fn screen_dist(shared: u32, a: Option<&MinimizerSketch>, b: Option<&MinimizerSketch>) -> f64 {
    match (a, b) {
        (Some(x), Some(y)) if x.total() > 0 && y.total() > 0 => {
            dist_from_shared(shared, x.total(), y.total())
        }
        _ => 0.0,
    }
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
            .find(|&(h, _)| h == poly_a)
            .map_or(0, |(_, c)| c as usize);
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

    /// The screen must fail OPEN on absent or empty sketches, on both paths.
    /// Failing closed there would drop pairs on no evidence at all.
    #[test]
    fn screen_fails_open_without_evidence() {
        let full = sketch(&make_seq(250, 5), MINIMIZER_K, MINIMIZER_W);
        let empty = sketch(&make_seq(4, 6), MINIMIZER_K, MINIMIZER_W);
        assert!(empty.is_empty());
        // Absent on either side.
        assert_eq!(screen_dist(0, None, Some(&full)), 0.0);
        assert_eq!(screen_dist(0, Some(&full), None), 0.0);
        // Present but empty on either side -- minimizer_dist would say 1.0 here.
        assert_eq!(minimizer_dist(&full, &empty), 1.0);
        assert_eq!(screen_dist(0, Some(&full), Some(&empty)), 0.0);
        assert_eq!(screen_dist(0, Some(&empty), Some(&full)), 0.0);
        // Two usable sketches still go through the real formula.
        assert_eq!(
            screen_dist(0, Some(&full), Some(&full)),
            dist_from_shared(0, full.total(), full.total())
        );
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

    /// The derived cutoff is the k-mer pass-rate quantile, so if the k-mer screen
    /// passes 30% of pairs the cutoff must sit at the 30th percentile of the
    /// minimizer distances.
    #[test]
    fn matched_pass_is_the_kmer_pass_quantile() {
        // k-mer distances 0.00..0.99; at cutoff 0.30 exactly 31 of 100 pass.
        let kd: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        // Minimizer distances 0.00..0.99 as well, so the quantile is readable.
        let md: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let r = matched_pass_cutoff(&kd, &md, 0.30);
        assert_eq!(r.n_pairs, 100);
        assert!((r.kmer_pass - 0.31).abs() < 1e-9, "{}", r.kmer_pass);
        assert!(
            (r.cutoff_rounded - 0.30).abs() < 1e-9,
            "{}",
            r.cutoff_rounded
        );
    }

    /// A cutoff must be derived, not transferred: the same pass rate maps to a
    /// very different number when the minimizer distribution is shifted, which is
    /// exactly why sharing `kmer_dist8`'s algebra transferred nothing.
    #[test]
    fn matched_pass_tracks_the_minimizer_distribution_not_the_kmer_one() {
        let kd: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let shifted: Vec<f64> = (0..100).map(|i| 0.40 + i as f64 / 250.0).collect();
        let a = matched_pass_cutoff(&kd, &kd, 0.30).cutoff_rounded;
        let b = matched_pass_cutoff(&kd, &shifted, 0.30).cutoff_rounded;
        assert!(b > a + 0.15, "{a} vs {b}");
    }

    /// Degenerate pass rates must still index in range rather than panic.
    #[test]
    fn matched_pass_handles_all_and_nothing() {
        let md: Vec<f64> = (0..50).map(|i| i as f64 / 50.0).collect();
        let none = vec![1.0; 50]; // k-mer screen passes nothing at 0.42
        let all = vec![0.0; 50]; // ... or everything
        let r0 = matched_pass_cutoff(&none, &md, 0.42);
        assert_eq!(r0.kmer_pass, 0.0);
        assert!(r0.cutoff.is_finite());
        let r1 = matched_pass_cutoff(&all, &md, 0.42);
        assert_eq!(r1.kmer_pass, 1.0);
        assert!((r1.cutoff - 0.98).abs() < 1e-9, "{}", r1.cutoff);

        let empty = matched_pass_cutoff(&[], &[], 0.42);
        assert_eq!(empty.n_pairs, 0);
        assert_eq!(empty.cutoff, MINIMIZER_KDIST_CUTOFF);
    }

    /// The score must reproduce the ORDER of the three measured configurations,
    /// and the default threshold must land the index where it actually won.
    ///
    /// | configuration | sharing | threads | nraw | measured setup/saved |
    /// |---|---|---|---|---|
    /// | pooled ITS2 | 1258 | 48 | 825,214 | 0.29 — index wins |
    /// | per-sample PacBio | 440 | 4 | 9,500 | 0.44 — index wins |
    /// | pooled PacBio | 6509 | 48 | 547,273 | 3.72 — index LOSES |
    #[test]
    fn score_ranks_the_measured_configurations() {
        let score =
            |sharing: f64, threads: usize, nraw: usize| sharing * threads as f64 / nraw as f64;
        let its2 = score(1258.0, 48, 825_214);
        let pacbio_per_sample = score(440.0, 4, 9_500);
        let pacbio_pooled = score(6509.0, 48, 547_273);

        assert!(its2 < pacbio_per_sample, "{its2} !< {pacbio_per_sample}");
        assert!(
            pacbio_per_sample < pacbio_pooled,
            "{pacbio_per_sample} !< {pacbio_pooled}"
        );
        // The two the index won stay under the default; the one it lost does not.
        assert!(its2 <= MINIMIZER_INDEX_MAX_SCORE);
        assert!(pacbio_per_sample <= MINIMIZER_INDEX_MAX_SCORE);
        assert!(pacbio_pooled > MINIMIZER_INDEX_MAX_SCORE);
    }

    /// A pool of identical sequences shares every minimizer, so sharing == nraw
    /// and the score is ~threads — the worst case, and it must decline.
    #[test]
    fn identical_pool_declines_the_index() {
        let seq = make_seq(600, 42);
        let sketches: Vec<Option<MinimizerSketch>> = (0..400)
            .map(|_| Some(sketch(&seq, MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let d = decide_index(&sketches, 48, None);
        assert!(!d.use_index, "score {} should decline", d.score);
        assert!(d.score > d.threshold);
    }

    /// A pool of unrelated sequences shares almost nothing, so posting lists are
    /// short and the index is the right call.
    ///
    /// Judged at 4 threads, not 48: `nraw` is the denominator, and 400 uniques is
    /// a pool that would never be given 48 threads. At 48 this scores 0.306 and
    /// declines — correctly, since a serial scatter cannot pay for itself across
    /// 400 raws no matter how diverse they are.
    #[test]
    fn diverse_pool_uses_the_index() {
        let sketches: Vec<Option<MinimizerSketch>> = (0..400)
            .map(|i| Some(sketch(&make_seq(600, 1000 + i), MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let d = decide_index(&sketches, 4, None);
        assert!(d.use_index, "score {} should accept", d.score);
        assert!(d.distinct > 0 && d.entries > 0);
    }

    /// Scale-free version of the same claim: at identical `threads` and `nraw`,
    /// sharing is the only thing that moves, and a diverse pool must score far
    /// below a degenerate one. This is the property the rule actually relies on.
    #[test]
    fn diversity_dominates_the_score() {
        let n = 400;
        let seq = make_seq(600, 77);
        let identical: Vec<Option<MinimizerSketch>> = (0..n)
            .map(|_| Some(sketch(&seq, MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let diverse: Vec<Option<MinimizerSketch>> = (0..n)
            .map(|i| Some(sketch(&make_seq(600, 5000 + i), MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let a = decide_index(&identical, 48, None).score;
        let b = decide_index(&diverse, 48, None).score;
        assert!(b * 20.0 < a, "diverse {b} vs identical {a}");
    }

    /// `entries` and `distinct` must equal what the index would report, or the
    /// score is computed from different quantities than the ones it models.
    #[test]
    fn decision_counts_match_the_built_index() {
        let sketches: Vec<Option<MinimizerSketch>> = (0..64)
            .map(|i| Some(sketch(&make_seq(400, 7000 + i), MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let d = decide_index(&sketches, 8, None);
        let idx = MinimizerIndex::build(&sketches);
        assert_eq!(d.entries, idx.n_postings());
        assert_eq!(d.distinct, idx.n_keys());
    }

    /// The env override must win in both directions, whatever the score says.
    #[test]
    fn forced_override_beats_the_score() {
        let seq = make_seq(600, 11);
        let sketches: Vec<Option<MinimizerSketch>> = (0..400)
            .map(|_| Some(sketch(&seq, MINIMIZER_K, MINIMIZER_W)))
            .collect();
        assert!(decide_index(&sketches, 48, Some(true)).use_index);
        assert!(!decide_index(&sketches, 48, Some(false)).use_index);

        let diverse: Vec<Option<MinimizerSketch>> = (0..400)
            .map(|i| Some(sketch(&make_seq(600, 2000 + i), MINIMIZER_K, MINIMIZER_W)))
            .collect();
        assert!(!decide_index(&diverse, 48, Some(false)).use_index);
        assert!(decide_index(&diverse, 48, Some(true)).use_index);
    }

    /// Threads is the pool that runs THIS compare. Reading the global count for a
    /// per-sample run inflates the score ~12x and would decline the index exactly
    /// where it measured a win, so the score must be sensitive to it.
    #[test]
    fn score_scales_with_threads() {
        let sketches: Vec<Option<MinimizerSketch>> = (0..200)
            .map(|i| Some(sketch(&make_seq(500, 3000 + i), MINIMIZER_K, MINIMIZER_W)))
            .collect();
        let a = decide_index(&sketches, 4, None).score;
        let b = decide_index(&sketches, 48, None).score;
        assert!((b / a - 12.0).abs() < 1e-9, "{a} -> {b}");
    }

    /// Empty and sketchless pools must not divide by zero.
    #[test]
    fn degenerate_pools_are_safe() {
        let d = decide_index(&[], 8, None);
        assert_eq!(d.score, 0.0);
        assert!(d.use_index, "an empty pool is trivially under threshold");

        let none: Vec<Option<MinimizerSketch>> = vec![None, None];
        let d = decide_index(&none, 8, None);
        assert_eq!(d.entries, 0);
        assert_eq!(d.distinct, 0);
        assert_eq!(d.sharing, 0.0);
        assert_eq!(d.score, 0.0);
    }
}

// ---------------------------------------------------------------------------
// Screen audit
// ---------------------------------------------------------------------------

/// Pair-level agreement between the two screen backends, for validating the
/// minimizer screen against the k-mer screen it would replace.
///
/// This is the instrument for the question that decides whether the backend is
/// safe: on a real pool, does the minimizer screen pass a **superset**, a
/// **subset**, or a **different set** of the pairs that `kdist` passes? Wall
/// clock and ASV counts both answer that only indirectly — ASV counts in
/// particular are a known trap here (`docs/findings/kmer-size-screening.md`),
/// because a screen change can churn intermediate clusters while the final table
/// stays the same size.
///
/// The counts alone would still be ambiguous, because a *disagreement* is only
/// bad if the pair mattered. So under audit the gate is the **union** of the two
/// screens — every pair either backend would align is aligned — and each
/// disagreement is bucketed by the substitution count the alignment actually
/// found. A pair the minimizer screen rejects at 1 substitution is a candidate
/// lost error-copy; the same rejection at 40 substitutions is the screen doing
/// its job better.
///
/// Global rather than threaded through the return types because this is a
/// diagnostic path only: it already pays two screens and aligns the union, so
/// relaxed atomic contention is not what makes it slow.
pub mod audit {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

    /// Substitution-count buckets: 0, 1, 2, 3, 4, 5, 6-10, >10.
    pub const NBUCKETS: usize = 8;

    #[inline]
    fn bucket(nsubs: usize) -> usize {
        match nsubs {
            0..=5 => nsubs,
            6..=10 => 6,
            _ => 7,
        }
    }

    static COMPARISONS: AtomicU64 = AtomicU64::new(0);
    static BOTH_PASS: AtomicU64 = AtomicU64::new(0);
    static KMER_ONLY: AtomicU64 = AtomicU64::new(0);
    static MINI_ONLY: AtomicU64 = AtomicU64::new(0);
    static NEITHER: AtomicU64 = AtomicU64::new(0);

    // Substitution histograms for the two disagreement classes and, for scale,
    // the agreeing-pass class.
    static H_KMER_ONLY: [AtomicU64; NBUCKETS] = [const { AtomicU64::new(0) }; NBUCKETS];
    static H_MINI_ONLY: [AtomicU64; NBUCKETS] = [const { AtomicU64::new(0) }; NBUCKETS];
    static H_BOTH: [AtomicU64; NBUCKETS] = [const { AtomicU64::new(0) }; NBUCKETS];

    /// Record one comparison's two verdicts and, when the union gate aligned it,
    /// the substitution count found. `nsubs` is `None` only when neither screen
    /// passed, in which case no alignment was performed.
    pub fn record(kmer_pass: bool, mini_pass: bool, nsubs: Option<usize>) {
        COMPARISONS.fetch_add(1, Relaxed);
        let (class, hist) = match (kmer_pass, mini_pass) {
            (true, true) => (&BOTH_PASS, Some(&H_BOTH)),
            (true, false) => (&KMER_ONLY, Some(&H_KMER_ONLY)),
            (false, true) => (&MINI_ONLY, Some(&H_MINI_ONLY)),
            (false, false) => (&NEITHER, None),
        };
        class.fetch_add(1, Relaxed);
        if let (Some(h), Some(n)) = (hist, nsubs) {
            h[bucket(n)].fetch_add(1, Relaxed);
        }
    }

    /// Snapshot of the counters, for reporting.
    #[derive(Clone, Copy, Debug, Default)]
    pub struct Summary {
        pub comparisons: u64,
        pub both_pass: u64,
        pub kmer_only: u64,
        pub mini_only: u64,
        pub neither: u64,
        pub h_kmer_only: [u64; NBUCKETS],
        pub h_mini_only: [u64; NBUCKETS],
        pub h_both: [u64; NBUCKETS],
    }

    fn load_hist(h: &[AtomicU64; NBUCKETS]) -> [u64; NBUCKETS] {
        std::array::from_fn(|i| h[i].load(Relaxed))
    }

    pub fn summary() -> Summary {
        Summary {
            comparisons: COMPARISONS.load(Relaxed),
            both_pass: BOTH_PASS.load(Relaxed),
            kmer_only: KMER_ONLY.load(Relaxed),
            mini_only: MINI_ONLY.load(Relaxed),
            neither: NEITHER.load(Relaxed),
            h_kmer_only: load_hist(&H_KMER_ONLY),
            h_mini_only: load_hist(&H_MINI_ONLY),
            h_both: load_hist(&H_BOTH),
        }
    }

    pub fn reset() {
        for c in [&COMPARISONS, &BOTH_PASS, &KMER_ONLY, &MINI_ONLY, &NEITHER] {
            c.store(0, Relaxed);
        }
        for h in [&H_KMER_ONLY, &H_MINI_ONLY, &H_BOTH] {
            for b in h {
                b.store(0, Relaxed);
            }
        }
    }

    /// Bucket labels, aligned with [`bucket`].
    pub const BUCKET_LABELS: [&str; NBUCKETS] = ["0", "1", "2", "3", "4", "5", "6-10", ">10"];

    impl Summary {
        /// Human-readable report, one line per fact worth reading.
        pub fn report(&self) -> String {
            let n = self.comparisons.max(1) as f64;
            let pct = |v: u64| 100.0 * v as f64 / n;
            let kmer_pass = self.both_pass + self.kmer_only;
            let mini_pass = self.both_pass + self.mini_only;
            let hist = |h: &[u64; NBUCKETS]| {
                BUCKET_LABELS
                    .iter()
                    .zip(h.iter())
                    .map(|(l, v)| format!("{l}:{v}"))
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            let mut s = String::new();
            s.push_str(&format!(
                "[screen-audit] {} comparisons; kmer passed {} ({:.2}%), minimizer passed {} ({:.2}%)\n",
                self.comparisons,
                kmer_pass,
                pct(kmer_pass),
                mini_pass,
                pct(mini_pass),
            ));
            s.push_str(&format!(
                "[screen-audit] agreement: both {} ({:.2}%), neither {} ({:.2}%), \
                 kmer-only {} ({:.4}%), minimizer-only {} ({:.4}%)\n",
                self.both_pass,
                pct(self.both_pass),
                self.neither,
                pct(self.neither),
                self.kmer_only,
                pct(self.kmer_only),
                self.mini_only,
                pct(self.mini_only),
            ));
            // The decisive line: a kmer-only pair at low substitution count is a
            // pair the minimizer screen would have shrouded but that the aligner
            // says is a plausible error copy.
            s.push_str(&format!(
                "[screen-audit] nsubs | kmer-only (minimizer would MISS these): {}\n",
                hist(&self.h_kmer_only)
            ));
            s.push_str(&format!(
                "[screen-audit] nsubs | minimizer-only (extra alignments): {}\n",
                hist(&self.h_mini_only)
            ));
            s.push_str(&format!(
                "[screen-audit] nsubs | both passed (for scale): {}",
                hist(&self.h_both)
            ));
            s
        }

        /// Substitution counts at or below which a shrouded pair could still
        /// have been a real error copy. Pairs the minimizer screen misses in
        /// this range are the ones that can cost ASVs.
        pub fn kmer_only_close_pairs(&self) -> u64 {
            self.h_kmer_only[..=3].iter().sum()
        }
    }
}
