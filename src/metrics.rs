//! Machine-readable run metrics (`--metrics-json`).
//!
//! The same quantities `--verbose` prints as prose, as a structured document.
//! This exists so tooling — `dev/run_screen_sweep*.sh`, the benchmark harness —
//! stops parsing `^\[dada\]` lines out of stderr. Once a consumer reads this
//! instead, the prose is free to change shape without breaking it (issue #162).
//!
//! Two rules hold this together:
//!
//! * **Collection is gated by `DadaParams::measure`, not `verbose`.** Printing
//!   and instrumentation are separate switches: `--metrics-json` turns
//!   measurement on without printing anything, and a future quiet tier can
//!   print without paying for timers. See [`crate::dada::run_dada`].
//! * **Nothing here is derived from the prose.** Both the prose and this
//!   document read the same accumulators, so they cannot drift apart.
//!
//! Durations are seconds (`f64`); rates are nanoseconds per unit, and are
//! `null` when their denominator is zero rather than `0.0`, so a consumer can
//! tell "not measured" from "measured as free".

use std::time::Duration;

use serde::Serialize;

/// Bumped when a field changes meaning or disappears. Additive fields do not
/// bump it: a consumer that ignores unknown keys keeps working.
pub const METRICS_SCHEMA_VERSION: u32 = 1;

/// How much instrumentation a run pays for.
///
/// Orthogonal to stderr verbosity: `--metrics-json` measures without printing,
/// and a quiet run can still print. The split exists because the two classes
/// have very different costs (issue #162).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum MeasureLevel {
    /// No instrumentation. Production speed.
    #[default]
    Off,
    /// Free tier: phase wall times, the store-loop counters and screened /
    /// aligned denominators (folded into the store loop by #143, so they are
    /// ALU on already-loaded values), the footprint walk, and the index probe.
    /// A few `Instant::now()` calls per bud round.
    Phases,
    /// Adds the per-comparison timers: `AlignBuffers::measuring()` and the
    /// per-item `t0`. These produce `busy`, map parallel efficiency and the
    /// screen / dp / al2subs split — at the price of 2-4 `Instant::now()` calls
    /// on *every* comparison, and comparisons reach 1.4e10 on a diverse pool.
    /// Never enable this on a run whose wall time you intend to report.
    Attribution,
}

impl MeasureLevel {
    /// Whether the per-comparison timers run. This is the value that reaches
    /// `b_compare_parallel`'s `measure` argument.
    pub fn attribution(self) -> bool {
        self >= MeasureLevel::Attribution
    }

    /// Whether anything is collected at all.
    pub fn enabled(self) -> bool {
        self > MeasureLevel::Off
    }
}

fn secs(d: Duration) -> f64 {
    d.as_secs_f64()
}

/// Nanoseconds per unit, or `None` when the denominator is zero.
fn per(d: Duration, n: u64) -> Option<f64> {
    (n > 0).then(|| d.as_secs_f64() * 1e9 / n as f64)
}

/// What the run was configured to do — the parameters that change the numbers
/// below, so an archived metrics file is self-describing.
#[derive(Debug, Clone, Default, Serialize)]
pub struct RunShape {
    pub nraw: usize,
    pub reads: u32,
    pub nclusters: usize,
    /// Pairwise alignments performed, and comparisons the screen shrouded.
    /// `u64` is required: these reach 1.4e10 on a diverse pool.
    pub nalign: u64,
    pub nshroud: u64,
    pub threads: usize,
    pub multithread: bool,
    pub align_backend: String,
    pub screen_backend: String,
    pub kmer_size: usize,
    pub kdist_cutoff: f64,
    pub band: i32,
    /// Only meaningful under the minimizer screen.
    pub minimizer_k: usize,
    pub minimizer_w: usize,
    /// Active `DADA2RS_*` tuning gates, as `gates::describe` renders them. A
    /// non-empty list means this run is not stock and its timings are not
    /// comparable with one that is.
    pub gates: Vec<String>,
}

/// Wall time in each phase of the bud loop. Serial except `compare.map`.
#[derive(Debug, Clone, Default, Serialize)]
pub struct PhaseTimes {
    pub compare: f64,
    pub shuffle: f64,
    pub bud: f64,
    pub p_update: f64,
    pub index_add: f64,
    /// The bud loop only. **Not** a supertotal of the phases above: the initial
    /// pre-loop compare is counted in `compare` but runs before this timer
    /// starts, so `compare > loop_total` is normal and not an inconsistency.
    pub loop_total: f64,
}

/// Where `b_compare` spent its time, and on how many comparisons.
#[derive(Debug, Clone, Default, Serialize)]
pub struct CompareMetrics {
    /// Wall time in the compare phase. Always real, whichever path ran.
    pub total: f64,
    /// Phase breakdown. `None` on a single-threaded run: the serial
    /// `b_compare` path returns no [`crate::cluster::CompareTiming`] at all, so
    /// its internals are *unmeasured*, not zero. Emitting zeros here would read
    /// as "the map was free".
    #[serde(skip_serializing_if = "Option::is_none")]
    pub attribution: Option<CompareAttribution>,
    /// Outside the compare timer, reported here because it tracks the same work.
    pub index_add_cluster: f64,
    /// Comparisons that reached the screen (were not greedy-skipped). Counted
    /// in the store-loop fold regardless of level, so this is always present —
    /// and the minimizer index probe reads it, which is why it must never
    /// become measurement-gated (issue #162).
    pub screened: u64,
    /// Comparisons that passed the screen and were aligned. Always present, as
    /// above.
    pub aligned: u64,
    /// Split of summed worker-busy time inside the parallel map. `None` unless
    /// the run measured at [`MeasureLevel::Attribution`] — absent here means
    /// "not measured", never "measured as zero".
    #[serde(skip_serializing_if = "Option::is_none")]
    pub split: Option<CompareSplit>,
    /// Summed per-item compute across all workers. Attribution level only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub busy: Option<f64>,
    /// `busy / (map × threads)`. Near 1.0 with low CPU utilization implies
    /// memory-bandwidth stalls; well below 1.0 implies threads idling on tail
    /// load-imbalance.
    pub map_parallel_efficiency: Option<f64>,
}

/// The phase breakdown of `b_compare_parallel`. `map` is the parallel pass;
/// everything else is serial. `unattributed` is `total` minus the named parts —
/// a large value means a phase is missing from this split, not that the work
/// was free.
#[derive(Debug, Clone, Default, Serialize)]
pub struct CompareAttribution {
    pub map: f64,
    pub reduction: f64,
    pub store: f64,
    pub free: f64,
    pub setup: f64,
    pub unattributed: f64,
    /// Raw-visits the store loop scanned.
    pub scanned: u64,
    /// Comparisons the store loop retained.
    pub stored: u64,
    pub ns_per_scanned: Option<f64>,
    pub ns_per_stored: Option<f64>,
}

/// The screen-vs-align question (#127): the screen runs on every comparison,
/// the aligner only on those that pass it, so this ratio decides whether the
/// lever is a cheaper screen or a cheaper aligner.
#[derive(Debug, Clone, Default, Serialize)]
pub struct CompareSplit {
    pub screen: f64,
    pub dp_kernel: f64,
    pub al2subs: f64,
    /// `busy` minus the three named parts.
    pub other: f64,
    pub screen_ns_per_screened: Option<f64>,
    pub dp_ns_per_aligned: Option<f64>,
    pub al2subs_ns_per_aligned: Option<f64>,
}

/// The three scan phases of `b_shuffle`, which have different access patterns:
/// `build` is sequential, `reconcile` is scattered, and a scattered comparison
/// costs roughly twice a sequential one — so comparison counts alone mis-size
/// any change here by about 2x (#87).
#[derive(Debug, Clone, Default, Serialize)]
pub struct ShuffleMetrics {
    pub build: f64,
    pub reconcile: f64,
    pub move_pass: f64,
    pub comps_build: u64,
    pub comps_reconcile: u64,
    pub build_ns_per_comp: Option<f64>,
    pub reconcile_ns_per_comp: Option<f64>,
    pub calls: u64,
    pub converge_calls: u64,
    pub builds: u64,
    pub moves: u64,
    pub zero_move_calls: u64,
    pub move_raws: u64,
}

/// `b_bud` scan redundancy: the bud pass rescans raws every round to find the
/// next cluster centre, and this is how much of that scanning repeats.
#[derive(Debug, Clone, Default, Serialize)]
pub struct BudMetrics {
    pub calls: u64,
    /// Calls that actually produced a new cluster.
    pub successes: u64,
    pub raws_scanned: u64,
    pub raws_per_call: Option<f64>,
}

/// p-update churn (#85): how many raws get repriced per round, and by which
/// exit path. `full_calc` is the only expensive one — the Poisson upper tail —
/// so a high ratio of cheap exits means a p-ordered structure would buy little.
#[derive(Debug, Clone, Default, Serialize)]
pub struct PUpdateMetrics {
    /// In-loop rounds only. The pre-loop call reprices the whole partition once
    /// and is not part of the per-bud churn.
    pub rounds: u64,
    pub repriced: u64,
    pub exit_singleton: u64,
    pub exit_center: u64,
    pub exit_zero_lambda: u64,
    pub full_calc: u64,
    pub lock_scanned: u64,
    pub dirty_clusters: u64,
    pub clean_clusters: u64,
}

/// Reconcile internals (#136): how much of the scattered reconcile pass is
/// recomputation that changes nothing, which bounds any future optimization.
#[derive(Debug, Clone, Default, Serialize)]
pub struct ReconcileMetrics {
    pub collect: f64,
    pub rescan: f64,
    /// Raws the reconcile recomputed.
    pub affected: u64,
    /// Of those, the ones whose best cluster actually changed.
    pub changed: u64,
    /// Raws that genuinely required a rescan.
    pub rescan_required: u64,
    pub rescan_comps: u64,
    /// Ties, where the incremental path cannot shortcut.
    pub ties: u64,
    pub pairs: u64,
}

/// Move-pass dirty-set pruning (#132).
#[derive(Debug, Clone, Default, Serialize)]
pub struct MovePruneMetrics {
    pub unpruned: u64,
    pub dirty: u64,
    pub prunable: u64,
    pub passes: u64,
}

/// Carrying `compmax` across buds (#139, reviving #87's projection): what the
/// carry relocates, split by access pattern.
#[derive(Debug, Clone, Default, Serialize)]
pub struct CarryMetrics {
    pub first_reconcile_pairs: u64,
    pub first_reconcile_comps: u64,
    pub first_reconcile_calls: u64,
}

/// Resident bytes per Raw (#32). In pooled mode this is the whole resident set;
/// in pseudo or per-sample mode it is one sample's share, so multiply by
/// `--sample-jobs` for the peak.
#[derive(Debug, Clone, Default, Serialize)]
pub struct FootprintMetrics {
    pub nraw: usize,
    pub seq_qual_bytes: u64,
    pub screen_vector_bytes: u64,
    pub screen_bytes_per_raw: Option<f64>,
    /// How the screen vectors are stored: `dense`, `sparse #43`, or
    /// `minimizer sketch k=../w=..`.
    pub screen_repr: String,
}

/// The minimizer index choice, which is **measured on the first clusters, not
/// predicted**, and reverses between workloads. Recorded because a run that
/// declines the index otherwise shows an unexplained `setup` near zero.
#[derive(Debug, Clone, Default, Serialize)]
pub struct IndexDecision {
    pub use_index: bool,
    /// `Some` when `DADA2RS_MINIMIZER_INDEX` overrode the probe.
    pub forced: Option<bool>,
    pub probed_clusters: usize,
    pub scatter_ns: f64,
    pub saving_ns: f64,
    pub merge_ns: f64,
    pub array_ns: f64,
    pub screened_comps: u64,
    pub threads: usize,
    pub distinct_minimizers: u64,
    pub postings: u64,
    pub mean_posting: f64,
    /// Set when the run's actual screened-comps-per-cluster disagreed with the
    /// probe's sample badly enough that the other choice would have been
    /// faster. A speed advisory: the screen is exact either way.
    pub hindsight_disagrees: bool,
}

/// Wall time of the pipeline stages around `run_dada`, which the denoiser
/// itself never sees. Present only for the fields the running subcommand
/// actually has.
#[derive(Debug, Clone, Default, Serialize)]
pub struct PipelineTimes {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub derep: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub merge: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dada: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<f64>,
}

/// Everything `run_dada` measured about one invocation.
#[derive(Debug, Clone, Default, Serialize)]
pub struct RunMetrics {
    pub run: RunShape,
    pub phases: PhaseTimes,
    pub compare: CompareMetrics,
    pub shuffle: ShuffleMetrics,
    pub reconcile: ReconcileMetrics,
    pub move_pruning: MovePruneMetrics,
    pub carry_87: CarryMetrics,
    pub bud: BudMetrics,
    pub p_update: PUpdateMetrics,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub footprint: Option<FootprintMetrics>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub index: Option<IndexDecision>,
}

/// The document written to `--metrics-json`.
///
/// One `run_dada` invocation per entry: `dada-pooled` denoises once and emits
/// one, `dada` over several inputs emits one per sample, and `dada-pseudo`
/// emits one per sample per round.
#[derive(Debug, Clone, Default, Serialize)]
pub struct MetricsDocument {
    pub schema_version: u32,
    /// What was instrumented. A consumer that finds `compare.split` absent
    /// should check this before concluding the split was zero.
    pub measure_level: MeasureLevel,
    /// Wall time of the whole subcommand, including everything outside
    /// `run_dada`.
    pub wall_seconds: f64,
    #[serde(skip_serializing_if = "PipelineTimes::is_empty")]
    pub pipeline: PipelineTimes,
    pub runs: Vec<LabelledRun>,
}

/// A `run_dada` invocation, tagged with which sample and round produced it.
#[derive(Debug, Clone, Serialize)]
pub struct LabelledRun {
    /// Sample name, or `__pooled__` for a pooled run's single invocation.
    pub sample: String,
    /// `dada-pseudo` round (1 or 2); absent for the other subcommands.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub round: Option<u8>,
    #[serde(flatten)]
    pub metrics: RunMetrics,
}

impl PipelineTimes {
    fn is_empty(&self) -> bool {
        self.derep.is_none() && self.merge.is_none() && self.dada.is_none() && self.output.is_none()
    }
}

impl MetricsDocument {
    pub fn new(wall: Duration, measure_level: MeasureLevel) -> Self {
        Self {
            schema_version: METRICS_SCHEMA_VERSION,
            measure_level,
            wall_seconds: secs(wall),
            pipeline: PipelineTimes::default(),
            runs: Vec::new(),
        }
    }

    pub fn push(&mut self, sample: impl Into<String>, round: Option<u8>, metrics: RunMetrics) {
        self.runs.push(LabelledRun {
            sample: sample.into(),
            round,
            metrics,
        });
    }
}

/// The accumulators `run_dada` fills as it goes.
///
/// Kept as raw `Duration`s and counts rather than formatted numbers, so the
/// prose printer and [`RunMetrics`] read exactly the same values and cannot
/// disagree. [`RawCounters::finish`] does every derivation in one place.
#[derive(Debug, Clone, Default)]
pub struct RawCounters {
    pub t_compare: Duration,
    pub t_shuffle: Duration,
    pub t_bud: Duration,
    pub t_pupdate: Duration,
    pub t_index_add: Duration,
    pub t_loop: Duration,

    pub t_cmp_map: Duration,
    pub t_cmp_serial: Duration,
    pub t_cmp_busy: Duration,
    pub t_cmp_agg: Duration,
    pub t_cmp_free: Duration,
    pub t_cmp_setup: Duration,
    pub t_cmp_screen: Duration,
    pub t_cmp_dp: Duration,
    pub t_cmp_post: Duration,
    pub n_cmp_scanned: u64,
    pub n_cmp_stored: u64,
    pub n_cmp_screened: u64,
    pub n_cmp_aligned: u64,

    pub t_shuf_build: Duration,
    pub t_shuf_reconcile: Duration,
    pub t_shuf_move: Duration,
    pub shuf_comps_build: u64,
    pub shuf_comps_reconcile: u64,
    pub shuf_calls: u64,
    pub shuf_converge_calls: u64,
    pub shuf_builds: u64,
    pub shuf_moves: u64,
    pub shuf_zero_move_calls: u64,
    pub shuf_move_raws: u64,

    pub t_rec_collect: Duration,
    pub t_rec_rescan: Duration,
    pub shuf_rec_affected: u64,
    pub shuf_rec_changed: u64,
    pub shuf_rec_rescan: u64,
    pub shuf_rec_rescan_comps: u64,
    pub shuf_rec_ties: u64,
    pub shuf_rec_pairs: u64,

    pub shuf_move_unpruned: u64,
    pub shuf_move_dirty: u64,
    pub shuf_move_prunable: u64,
    pub shuf_move_passes: u64,

    pub shuf_first_rec_pairs: u64,
    pub shuf_first_rec_comps: u64,
    pub shuf_first_rec_calls: u64,

    pub bud_calls: u64,
    pub bud_success: u64,
    pub bud_raws_scanned: u64,

    pub pupd_rounds: u64,
    pub pupd: PUpdateMetrics,

    pub footprint: Option<FootprintMetrics>,
}

impl RawCounters {
    /// Derive the serializable document. `threads` is the rayon pool width the
    /// map actually ran on, which is not always `--threads`.
    pub fn finish(&self, run: RunShape, threads: usize, level: MeasureLevel) -> RunMetrics {
        let multithread = run.multithread;
        let busy = self.t_cmp_busy;
        let cmp_total = self.t_compare;
        let named = self.t_cmp_map
            + self.t_cmp_serial
            + self.t_cmp_agg
            + self.t_cmp_free
            + self.t_cmp_setup;
        let other_busy = busy
            .saturating_sub(self.t_cmp_screen)
            .saturating_sub(self.t_cmp_dp)
            .saturating_sub(self.t_cmp_post);

        RunMetrics {
            run,
            phases: PhaseTimes {
                compare: secs(self.t_compare),
                shuffle: secs(self.t_shuffle),
                bud: secs(self.t_bud),
                p_update: secs(self.t_pupdate),
                index_add: secs(self.t_index_add),
                loop_total: secs(self.t_loop),
            },
            compare: CompareMetrics {
                total: secs(cmp_total),
                attribution: multithread.then(|| CompareAttribution {
                    map: secs(self.t_cmp_map),
                    reduction: secs(self.t_cmp_agg),
                    store: secs(self.t_cmp_serial),
                    free: secs(self.t_cmp_free),
                    setup: secs(self.t_cmp_setup),
                    unattributed: secs(cmp_total.saturating_sub(named)),
                    scanned: self.n_cmp_scanned,
                    stored: self.n_cmp_stored,
                    ns_per_scanned: per(self.t_cmp_serial, self.n_cmp_scanned),
                    ns_per_stored: per(self.t_cmp_serial, self.n_cmp_stored),
                }),
                index_add_cluster: secs(self.t_index_add),
                screened: self.n_cmp_screened,
                aligned: self.n_cmp_aligned,
                // Everything below is produced by the per-comparison timers, so
                // it is absent — not zero — at the cheap level.
                busy: level.attribution().then(|| secs(busy)),
                map_parallel_efficiency: (level.attribution()
                    && self.t_cmp_map.as_secs_f64() > 0.0
                    && threads > 0)
                    .then(|| busy.as_secs_f64() / (self.t_cmp_map.as_secs_f64() * threads as f64)),
                split: level.attribution().then(|| CompareSplit {
                    screen: secs(self.t_cmp_screen),
                    dp_kernel: secs(self.t_cmp_dp),
                    al2subs: secs(self.t_cmp_post),
                    other: secs(other_busy),
                    screen_ns_per_screened: per(self.t_cmp_screen, self.n_cmp_screened),
                    dp_ns_per_aligned: per(self.t_cmp_dp, self.n_cmp_aligned),
                    al2subs_ns_per_aligned: per(self.t_cmp_post, self.n_cmp_aligned),
                }),
            },
            shuffle: ShuffleMetrics {
                build: secs(self.t_shuf_build),
                reconcile: secs(self.t_shuf_reconcile),
                move_pass: secs(self.t_shuf_move),
                comps_build: self.shuf_comps_build,
                comps_reconcile: self.shuf_comps_reconcile,
                build_ns_per_comp: per(self.t_shuf_build, self.shuf_comps_build),
                reconcile_ns_per_comp: per(self.t_shuf_reconcile, self.shuf_comps_reconcile),
                calls: self.shuf_calls,
                converge_calls: self.shuf_converge_calls,
                builds: self.shuf_builds,
                moves: self.shuf_moves,
                zero_move_calls: self.shuf_zero_move_calls,
                move_raws: self.shuf_move_raws,
            },
            reconcile: ReconcileMetrics {
                collect: secs(self.t_rec_collect),
                rescan: secs(self.t_rec_rescan),
                affected: self.shuf_rec_affected,
                changed: self.shuf_rec_changed,
                rescan_required: self.shuf_rec_rescan,
                rescan_comps: self.shuf_rec_rescan_comps,
                ties: self.shuf_rec_ties,
                pairs: self.shuf_rec_pairs,
            },
            move_pruning: MovePruneMetrics {
                unpruned: self.shuf_move_unpruned,
                dirty: self.shuf_move_dirty,
                prunable: self.shuf_move_prunable,
                passes: self.shuf_move_passes,
            },
            carry_87: CarryMetrics {
                first_reconcile_pairs: self.shuf_first_rec_pairs,
                first_reconcile_comps: self.shuf_first_rec_comps,
                first_reconcile_calls: self.shuf_first_rec_calls,
            },
            bud: BudMetrics {
                calls: self.bud_calls,
                successes: self.bud_success,
                raws_scanned: self.bud_raws_scanned,
                raws_per_call: (self.bud_calls > 0)
                    .then(|| self.bud_raws_scanned as f64 / self.bud_calls as f64),
            },
            p_update: PUpdateMetrics {
                rounds: self.pupd_rounds,
                ..self.pupd.clone()
            },
            footprint: self.footprint.clone(),
            index: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn full() -> MeasureLevel {
        MeasureLevel::Attribution
    }

    #[test]
    fn rates_are_null_not_zero_when_unmeasured() {
        let c = RawCounters::default();
        let m = c.finish(RunShape::default(), 4, full());
        // A zero denominator must be distinguishable from a measured zero.
        assert!(m.compare.attribution.is_none());
        assert!(
            m.compare
                .split
                .as_ref()
                .unwrap()
                .dp_ns_per_aligned
                .is_none()
        );
        assert_eq!(m.compare.screened, 0);
        // map == 0 means the map never ran, so efficiency is unknown, not 0%.
        assert!(m.compare.map_parallel_efficiency.is_none());
    }

    #[test]
    fn rates_are_computed_when_denominators_are_present() {
        let c = RawCounters {
            t_cmp_serial: Duration::from_nanos(1000),
            n_cmp_scanned: 10,
            t_cmp_map: Duration::from_secs(1),
            t_cmp_busy: Duration::from_secs(2),
            ..RawCounters::default()
        };
        let m = c.finish(
            RunShape {
                multithread: true,
                ..RunShape::default()
            },
            4,
            full(),
        );
        assert_eq!(m.compare.attribution.unwrap().ns_per_scanned, Some(100.0));
        // busy / (map * threads) = 2 / (1 * 4)
        assert_eq!(m.compare.map_parallel_efficiency, Some(0.5));
    }

    #[test]
    fn unattributed_is_the_residual_and_never_negative() {
        let c = RawCounters {
            t_compare: Duration::from_secs(10),
            t_cmp_map: Duration::from_secs(6),
            t_cmp_serial: Duration::from_secs(1),
            ..RawCounters::default()
        };
        let mt = || RunShape {
            multithread: true,
            ..RunShape::default()
        };
        let m = c.finish(mt(), 1, full());
        assert_eq!(m.compare.attribution.unwrap().unattributed, 3.0);

        // Named parts exceeding the total (nested timers) must clamp, not wrap.
        let c = RawCounters {
            t_compare: Duration::from_secs(1),
            t_cmp_map: Duration::from_secs(5),
            ..RawCounters::default()
        };
        let m = c.finish(mt(), 1, full());
        assert_eq!(m.compare.attribution.unwrap().unattributed, 0.0);
    }

    /// The whole point of the cheap level: an absent split must not be
    /// readable as a measured zero, or a consumer will average it in.
    #[test]
    fn cheap_level_omits_attribution_entirely() {
        let c = RawCounters {
            t_cmp_map: Duration::from_secs(4),
            t_cmp_busy: Duration::from_secs(8),
            t_cmp_screen: Duration::from_secs(3),
            n_cmp_screened: 100,
            ..RawCounters::default()
        };

        let mt = RunShape {
            multithread: true,
            ..RunShape::default()
        };
        let cheap = c.finish(mt.clone(), 2, MeasureLevel::Phases);
        assert!(cheap.compare.split.is_none());
        assert!(cheap.compare.busy.is_none());
        assert!(cheap.compare.map_parallel_efficiency.is_none());
        // ...but the free tier survives, including the denominators the index
        // probe depends on.
        assert_eq!(cheap.compare.attribution.as_ref().unwrap().map, 4.0);
        assert_eq!(cheap.compare.screened, 100);

        let j = serde_json::to_value(&cheap).unwrap();
        assert!(j["compare"].get("split").is_none());
        assert!(j["compare"].get("busy").is_none());

        // Same counters at the full level do surface it.
        let full = c.finish(mt, 2, MeasureLevel::Attribution);
        assert_eq!(full.compare.busy, Some(8.0));
        assert_eq!(full.compare.split.as_ref().unwrap().screen, 3.0);
    }

    /// The serial `b_compare` path returns no timing, so its breakdown must be
    /// absent rather than a pile of zeros a consumer would average in.
    #[test]
    fn single_threaded_run_omits_the_compare_breakdown() {
        let c = RawCounters {
            t_compare: Duration::from_secs(9),
            ..RawCounters::default()
        };
        let serial = c.finish(RunShape::default(), 1, MeasureLevel::Attribution);
        assert!(serial.compare.attribution.is_none());
        // The phase total is still real and must survive.
        assert_eq!(serial.compare.total, 9.0);
        let j = serde_json::to_value(&serial).unwrap();
        assert!(j["compare"].get("attribution").is_none());
        assert_eq!(j["compare"]["total"], 9.0);
    }

    #[test]
    fn measure_level_orders_by_cost() {
        assert!(MeasureLevel::Off < MeasureLevel::Phases);
        assert!(MeasureLevel::Phases < MeasureLevel::Attribution);
        assert!(!MeasureLevel::Off.enabled());
        assert!(MeasureLevel::Phases.enabled());
        assert!(!MeasureLevel::Phases.attribution());
        assert!(MeasureLevel::Attribution.attribution());
        assert_eq!(MeasureLevel::default(), MeasureLevel::Off);
    }

    #[test]
    fn document_omits_empty_pipeline_but_keeps_runs() {
        let doc = MetricsDocument::new(Duration::from_secs(3), MeasureLevel::Phases);
        let j = serde_json::to_value(&doc).unwrap();
        assert_eq!(j["schema_version"], METRICS_SCHEMA_VERSION);
        assert_eq!(j["wall_seconds"], 3.0);
        assert_eq!(j["measure_level"], "phases");
        assert!(j.get("pipeline").is_none());
        assert!(j["runs"].as_array().unwrap().is_empty());
    }

    #[test]
    fn labelled_run_flattens_metrics_and_drops_absent_round() {
        let mut doc = MetricsDocument::new(Duration::ZERO, MeasureLevel::Attribution);
        doc.push("s1", None, RunMetrics::default());
        doc.push("s2", Some(2), RunMetrics::default());
        let j = serde_json::to_value(&doc).unwrap();
        let runs = j["runs"].as_array().unwrap();
        assert_eq!(runs[0]["sample"], "s1");
        assert!(runs[0].get("round").is_none());
        // flattened, not nested under `metrics`
        assert!(runs[0].get("phases").is_some());
        assert!(runs[0].get("metrics").is_none());
        assert_eq!(runs[1]["round"], 2);
    }
}
