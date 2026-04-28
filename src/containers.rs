/// Default initial cluster buffer size (mirrors C++ RAWBUF / CLUSTBUF).
/// Vec growth is automatic in Rust so this only affects the initial allocation.
const INIT_RAWS_CAPACITY: usize = 50;
const INIT_CLUSTERS_CAPACITY: usize = 50;

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

/// Brief summary of a comparison between a cluster center and a Raw.
/// Equivalent to the C++ `Comparison` struct.
#[derive(Debug, Clone, Default)]
pub struct Comparison {
    pub i: u32,
    pub index: u32,
    pub lambda: f64,
    pub hamming: u32,
}

// ---------------------------------------------------------------------------
// Sub
// ---------------------------------------------------------------------------

/// Compressed substitutions between two aligned sequences.
/// All positions are 0-indexed in the alignment.
/// Equivalent to the C++ `Sub` struct.
#[derive(Debug, Clone)]
pub struct Sub {
    /// Length of the reference sequence.
    pub len0: u32,
    /// Map from each reference-sequence position to its column in the alignment.
    pub map: Vec<u16>,
    /// Reference-sequence positions of each substitution.
    pub pos: Vec<u16>,
    /// Reference nucleotide at each substitution (integer-encoded).
    pub nt0: Vec<u8>,
    /// Query nucleotide at each substitution (integer-encoded).
    pub nt1: Vec<u8>,
    /// Quality score at each substitution position in the reference.
    pub q0: Vec<u8>,
    /// Quality score at each substitution position in the query.
    pub q1: Vec<u8>,
}

impl Sub {
    pub fn nsubs(&self) -> usize {
        self.pos.len()
    }
}

// ---------------------------------------------------------------------------
// BirthType
// ---------------------------------------------------------------------------

/// How a cluster was created. Equivalent to the C++ `birth_type` char field.
#[derive(Debug, Clone, Default)]
pub enum BirthType {
    /// First cluster, created at initialization ("I").
    #[default]
    Initial,
    /// Born from an abundance p-value test ("A").
    Abundance,
    /// Born from a prior-sequence p-value test ("P").
    Prior,
    /// Born from a singleton p-value test ("S").
    Singleton,
}

// ---------------------------------------------------------------------------
// Raw
// ---------------------------------------------------------------------------

/// A unique sequence variant with its quality and abundance data.
/// Equivalent to the C++ `Raw` struct.
pub struct Raw {
    /// Sequence encoded as A=1, C=2, G=3, T=4 (mirrors C++ char* seq).
    pub seq: Vec<u8>,
    /// Per-position quality scores rounded to u8; `None` when qualities are
    /// unavailable. Stored as u8 to match the C++ memory-saving approach.
    pub qual: Option<Vec<u8>>,
    /// Prior reasons to expect this sequence to be genuine.
    pub prior: bool,
    /// 16-bit k-mer frequency vector; populated by kmers module.
    pub kmer: Option<Vec<u16>>,
    /// 8-bit compressed k-mer frequency vector; populated by kmers module.
    pub kmer8: Option<Vec<u8>>,
    /// K-mers in the order they appear along the sequence; populated by kmers module.
    pub kord: Option<Vec<u16>>,
    /// Number of reads of this unique sequence.
    pub reads: u32,
    /// Index of this Raw in `B.raws`.
    pub index: u32,
    /// Abundance p-value relative to the current cluster.
    pub p: f64,
    /// Sentinel value used during min/max expected-abundance calculations.
    pub e_minmax: f64,
    /// Most recent comparison result against a cluster center.
    pub comp: Comparison,
    /// When true, this Raw is locked to its current cluster in greedy mode.
    pub lock: bool,
    /// When true, this Raw will be error-corrected to its cluster center.
    pub correct: bool,
}

impl Raw {
    /// Construct a Raw from an already-encoded sequence and optional quality scores.
    ///
    /// `seq` must already use the integer encoding (A=1, C=2, G=3, T=4).
    /// `qual` values are rounded to `u8`, matching the C++ implementation.
    pub fn new(seq: Vec<u8>, qual: Option<&[f64]>, reads: u32, prior: bool) -> Self {
        Raw {
            qual: qual.map(|q| q.iter().map(|&v| v.round() as u8).collect()),
            seq,
            prior,
            kmer: None,
            kmer8: None,
            kord: None,
            reads,
            index: 0,
            p: 0.0,
            e_minmax: -999.0,
            comp: Comparison::default(),
            lock: false,
            correct: true,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Reset per-iteration mutable state so this Raw can be fed back into a
    /// fresh DADA run without re-encoding the sequence or recomputing k-mer
    /// vectors. Leaves `seq`, `qual`, `kmer*`, `kord`, `reads`, `prior`
    /// intact — those are fixed for the life of the input. `index` is
    /// reassigned by `B::new`.
    pub fn reset_for_iteration(&mut self) {
        self.p = 0.0;
        self.e_minmax = -999.0;
        self.comp = Comparison::default();
        self.lock = false;
        self.correct = true;
    }
}

// ---------------------------------------------------------------------------
// Bi
// ---------------------------------------------------------------------------

/// A single cluster in the partition.
/// Equivalent to the C++ `Bi` struct.
pub struct Bi {
    /// Integer-encoded representative sequence for this cluster.
    pub seq: Vec<u8>,
    /// Index of the representative Raw in `B.raws`; `None` until assigned.
    pub center: Option<usize>,
    /// Indices of member Raws in `B.raws`.
    pub raws: Vec<usize>,
    /// Sum of reads across all member Raws.
    pub reads: u32,
    /// Index of this cluster in `B.clusters`.
    pub i: u32,
    /// Recalculate expected abundances on next pass when true.
    pub update_e: bool,
    /// Raws should be reassigned between clusters when true.
    pub shuffle: bool,
    /// Check whether Raws should be locked to this cluster when true.
    pub check_locks: bool,
    /// Self-production genotype error probability.
    pub self_: f64,
    /// Total number of Raws across the entire partition (used as p-value denominator).
    pub totraw: u32,
    pub birth_type: BirthType,
    /// Index of the cluster from which this one was split.
    pub birth_from: u32,
    /// Bonferroni-corrected p-value that triggered this cluster's creation.
    pub birth_pval: f64,
    /// Fold-enrichment above expectation at birth.
    pub birth_fold: f64,
    /// Expected read count at the time of birth.
    pub birth_e: f64,
    pub birth_comp: Comparison,
    /// Comparisons with all Raws that could potentially join this cluster.
    pub comp: Vec<Comparison>,
}

impl Bi {
    pub fn new(totraw: u32) -> Self {
        Bi {
            seq: Vec::new(),
            center: None,
            raws: Vec::with_capacity(INIT_RAWS_CAPACITY),
            reads: 0,
            i: 0,
            update_e: true,
            shuffle: true,
            check_locks: true,
            self_: 0.0,
            totraw,
            birth_type: BirthType::default(),
            birth_from: 0,
            birth_pval: 0.0,
            birth_fold: 1.0,
            birth_e: 0.0,
            birth_comp: Comparison::default(),
            comp: Vec::new(),
        }
    }

    pub fn nraw(&self) -> usize {
        self.raws.len()
    }
}

// ---------------------------------------------------------------------------
// B
// ---------------------------------------------------------------------------

/// The full partition of reads into clusters.
/// Equivalent to the C++ `B` struct.
pub struct B {
    /// All unique sequences, in their original input order.
    pub raws: Vec<Raw>,
    /// All clusters.
    pub clusters: Vec<Bi>,
    /// Total read count across all Raws.
    pub reads: u32,
    /// Number of pairwise alignments performed.
    pub nalign: u32,
    /// Number of comparisons screened out by k-mer distance ("shrouded").
    pub nshroud: u32,
    /// Significance threshold for abundance-based cluster splitting.
    pub omega_a: f64,
    /// Significance threshold for singleton detection.
    pub omega_p: f64,
    pub use_quals: bool,
    /// Sorted unique lambda values (for CDF precomputation).
    pub lams: Vec<f64>,
    /// CDF values corresponding to `lams`.
    pub cdf: Vec<f64>,
}

impl B {
    /// Create a new partition and initialize it with a single cluster containing
    /// all Raws.
    ///
    /// `raws` must carry already-encoded sequences. Raw indices are assigned here.
    pub fn new(mut raws: Vec<Raw>, omega_a: f64, omega_p: f64, use_quals: bool) -> Self {
        let reads = raws.iter().map(|r| r.reads).sum();
        for (i, raw) in raws.iter_mut().enumerate() {
            raw.index = i as u32;
        }
        let mut b = B {
            raws,
            clusters: Vec::with_capacity(INIT_CLUSTERS_CAPACITY),
            reads,
            nalign: 0,
            nshroud: 0,
            omega_a,
            omega_p,
            use_quals,
            lams: Vec::new(),
            cdf: Vec::new(),
        };
        b.init();
        b
    }

    /// Reset to a single cluster containing all Raws.
    /// Equivalent to C++ `b_init`.
    pub fn init(&mut self) {
        self.clusters.clear();
        let nraw = self.raws.len() as u32;

        let mut bi = Bi::new(nraw);
        bi.birth_type = BirthType::Initial;
        bi.birth_e = self.reads as f64;
        self.clusters.push(bi);
        self.nalign = 0;
        self.nshroud = 0;

        // Populate the single cluster. Collect (index, reads) first to avoid
        // holding a borrow on self.raws while mutating self.clusters.
        let members: Vec<(usize, u32)> = self
            .raws
            .iter()
            .map(|r| (r.index as usize, r.reads))
            .collect();
        for (idx, reads) in members {
            self.clusters[0].raws.push(idx);
            self.clusters[0].reads += reads;
        }
        self.clusters[0].update_e = true;

        self.census(0);
        self.assign_center(0);
    }

    /// Add a new cluster. Sets its `i` field and returns its index.
    /// Equivalent to C++ `b_add_bi`.
    pub fn add_cluster(&mut self, mut bi: Bi) -> usize {
        let idx = self.clusters.len();
        bi.i = idx as u32;
        self.clusters.push(bi);
        idx
    }

    /// Add Raw `raw_idx` to cluster `bi_idx`. Updates cluster reads and flags.
    /// Equivalent to C++ `bi_add_raw`.
    pub fn bi_add_raw(&mut self, bi_idx: usize, raw_idx: usize) {
        let reads = self.raws[raw_idx].reads;
        self.clusters[bi_idx].raws.push(raw_idx);
        self.clusters[bi_idx].reads += reads;
        self.clusters[bi_idx].update_e = true;
    }

    /// Remove the Raw at position `r` within cluster `bi_idx` using swap-remove,
    /// preserving O(1) removal at the cost of reordering. Returns the removed
    /// Raw's index in `B.raws`.
    ///
    /// Mirrors the C++ `bi_pop_raw` behaviour: the popped slot is filled by the
    /// last element, so caller must not rely on positional stability.
    /// Equivalent to C++ `bi_pop_raw`.
    pub fn bi_pop_raw(&mut self, bi_idx: usize, r: usize) -> usize {
        assert!(
            r < self.clusters[bi_idx].raws.len(),
            "bi_pop_raw: index {r} out of range (nraw={})",
            self.clusters[bi_idx].raws.len()
        );
        let raw_idx = self.clusters[bi_idx].raws.swap_remove(r);
        let reads = self.raws[raw_idx].reads;
        self.clusters[bi_idx].reads -= reads;
        self.clusters[bi_idx].update_e = true;
        raw_idx
    }

    /// Recompute the read count for cluster `i` from its member Raws.
    /// Equivalent to C++ `bi_census`.
    pub fn census(&mut self, i: usize) {
        let reads: u32 = {
            let raws = &self.raws;
            self.clusters[i].raws.iter().map(|&ri| raws[ri].reads).sum()
        };
        self.clusters[i].reads = reads;
    }

    /// Set the most-abundant Raw as the center of cluster `i` and copy its
    /// encoded sequence into `Bi.seq`.
    /// Equivalent to C++ `bi_assign_center`.
    pub fn assign_center(&mut self, i: usize) {
        let ci: Option<usize> = {
            let raws = &self.raws;
            self.clusters[i]
                .raws
                .iter()
                .copied()
                .max_by_key(|&ri| raws[ri].reads)
        };
        if let Some(ci) = ci {
            self.clusters[i].center = Some(ci);
            self.clusters[i].seq = self.raws[ci].seq.clone();
        }
    }
}
