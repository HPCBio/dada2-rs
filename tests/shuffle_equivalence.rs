//! Shuffle optimisations that must not change results.
//!
//! Each submodule pins one work reduction in the pooled shuffle against the
//! code path it replaced: both arms run from one binary, selected by an env
//! gate, and their full per-sample output must be byte-identical.

mod common;

use std::path::Path;
use std::process::Command;

use common::{BIN, fixture, learn_errors};

/// Run `dada-pooled` on the two committed fixtures into `out`, with `env` set.
fn pooled(out: &Path, threads: &str, env: &[(&str, &str)]) {
    let errs = learn_errors(out, threads);
    let res = Command::new(BIN)
        .args(["dada-pooled", "--threads", threads, "--error-model"])
        .arg(&errs)
        .arg("--output-dir")
        .arg(out)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .envs(env.iter().copied())
        .output()
        .expect("dada-pooled");
    assert!(
        res.status.success(),
        "dada-pooled failed ({env:?}): {}",
        String::from_utf8_lossy(&res.stderr)
    );
}

mod move_prune {
    //! Equivalence tests for the dirty-cluster move pass (issue #132).
    //!
    //! The shuffle's move pass visits only the clusters holding raws whose best
    //! cluster changed in the preceding reconcile, instead of every cluster every
    //! iteration. That is a pure work reduction: the partition it produces must be
    //! **identical** to the full scan's, because the two differ only in which
    //! clusters are examined, never in what happens to a raw once examined.
    //!
    //! This is the invariant #124's pruning arm nearly broke — a byte-identical
    //! partition depends on ascending-`ci` order, strict `>`, and a lowest-`ci`
    //! tie-break, and a subtle violation there survived both benchmarking and ASV
    //! concordance before exact-equality testing caught it. So these tests compare
    //! the full per-sample output, not ASV counts.
    //!
    //! `DADA2RS_SHUFFLE_NO_PRUNE=1` forces the full scan, so both arms run from one
    //! binary — which also removes the failure mode where an A/B is built from the
    //! wrong checkout and silently measures the same code twice.

    use std::path::Path;

    use crate::common::{assert_same_outputs, scratch};

    fn run_pooled(out: &Path, prune: bool, threads: &str) {
        let env: &[_] = if prune {
            &[]
        } else {
            &[("DADA2RS_SHUFFLE_NO_PRUNE", "1")]
        };
        super::pooled(out, threads, env);
    }

    /// The pruned and unpruned move passes must produce identical output.
    #[test]
    fn dirty_cluster_prune_matches_full_scan() {
        let tmp = scratch("prune");
        let (pruned, full) = (tmp.join("pruned"), tmp.join("full"));
        std::fs::create_dir_all(&pruned).unwrap();
        std::fs::create_dir_all(&full).unwrap();

        run_pooled(&pruned, true, "1");
        run_pooled(&full, false, "1");
        assert_same_outputs(&pruned, &full, "pruned vs full move pass, single-threaded");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
    /// stored in, which is upstream of the move pass — so this exercises the prune
    /// against a different (still deterministic) input ordering.
    #[test]
    fn dirty_cluster_prune_matches_full_scan_threaded() {
        let tmp = scratch("prune_mt");
        let (pruned, full) = (tmp.join("pruned"), tmp.join("full"));
        std::fs::create_dir_all(&pruned).unwrap();
        std::fs::create_dir_all(&full).unwrap();

        run_pooled(&pruned, true, "4");
        run_pooled(&full, false, "4");
        assert_same_outputs(&pruned, &full, "pruned vs full move pass, 4 threads");

        let _ = std::fs::remove_dir_all(&tmp);
    }
}

mod reconcile {
    //! Equivalence tests for the incremental reconcile (issue #136).
    //!
    //! The reconcile updates each affected raw's best cluster incrementally --
    //! testing only the *changed* candidates against a carried incumbent -- instead
    //! of rescanning every affected raw's whole candidate list. It reaches the same
    //! answer by a different route, so the equality is not structural the way the
    //! move-pass prune's was: it depends on the replacement rule reproducing
    //! `best_from_cands` exactly, including the lowest-`ci` tie-break on exact
    //! float equality.
    //!
    //! That tie-break is exactly where #124's pruning arm hid a latent bug, which
    //! survived both benchmarking and ASV concordance. So these tests compare full
    //! per-sample output, not ASV counts -- and the stronger check is the
    //! `debug_assertions` invariant in `b_shuffle_converge`, which compares against
    //! a real full rescan every iteration.
    //!
    //! `DADA2RS_RECONCILE_FULL=1` forces the old full-rescan path, so both arms run
    //! from one binary -- which also removes the failure mode where an A/B is built
    //! from the wrong checkout and silently measures the same code twice.

    use std::path::Path;

    use crate::common::{assert_same_outputs, scratch};

    fn run_pooled(out: &Path, incremental: bool, threads: &str) {
        let env: &[_] = if incremental {
            &[]
        } else {
            &[("DADA2RS_RECONCILE_FULL", "1")]
        };
        super::pooled(out, threads, env);
    }

    /// The incremental and full reconcilees must produce identical output.
    #[test]
    fn incremental_reconcile_matches_full_rescan() {
        let tmp = scratch("recon");
        let (incremental, full) = (tmp.join("incremental"), tmp.join("full"));
        std::fs::create_dir_all(&incremental).unwrap();
        std::fs::create_dir_all(&full).unwrap();

        run_pooled(&incremental, true, "1");
        run_pooled(&full, false, "1");
        assert_same_outputs(
            &incremental,
            &full,
            "incremental vs full reconcile, single-threaded",
        );

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
    /// stored in, which is upstream of the reconcile -- so this exercises the
    /// incremental rule against a different (still deterministic) candidate order,
    /// where a tie-break error is more likely to surface.
    #[test]
    fn incremental_reconcile_matches_full_rescan_threaded() {
        let tmp = scratch("recon_mt");
        let (incremental, full) = (tmp.join("incremental"), tmp.join("full"));
        std::fs::create_dir_all(&incremental).unwrap();
        std::fs::create_dir_all(&full).unwrap();

        run_pooled(&incremental, true, "4");
        run_pooled(&full, false, "4");
        assert_same_outputs(
            &incremental,
            &full,
            "incremental vs full reconcile, 4 threads",
        );

        let _ = std::fs::remove_dir_all(&tmp);
    }
}

mod carry {
    //! Equivalence tests for carrying `compmax` across bud rounds (issue #139).
    //!
    //! The shuffle normally rebuilds its `compmax`/`emax` map from scratch on every
    //! `b_shuffle_converge` call — one full cluster-major pass over every comp in
    //! the pool, per bud. The carry (on by default; `DADA2RS_SHUFFLE_NO_CARRY=1`
    //! disables it) keeps the map alive across buds instead, letting the first reconcile of the next call repair the parts the
    //! bud invalidated (the parent cluster's fallen reads, the new cluster's
    //! never-folded comps).
    //!
    //! That is a pure work reduction *if and only if* the carried map ends up
    //! exactly where a rebuild would have put it: same argmax, same lowest-`ci`
    //! tie-break. The partition must therefore be **byte-identical** between arms.
    //!
    //! Why this needs more than a fixture A/B: #136's four seeded mutations all
    //! survived end-to-end equivalence on these same fixtures, because a small pool
    //! never reaches the states that separate the routes. So the real guard is
    //! `DADA2RS_RECONCILE_VERIFY=1`, which asserts after *every* reconcile that
    //! `compmax`/`emax` equal a full rescan — turned on here, and the only check
    //! that can see a carry-specific divergence the moment it happens rather than
    //! if it survives to the output. The carry makes that assertion strictly
    //! stronger, because the state it validates now spans bud rounds.

    use std::path::Path;

    use crate::common::{assert_same_outputs, scratch};

    fn run_pooled(out: &Path, carry: bool, threads: &str, verify: bool) {
        let mut env = vec![];
        if !carry {
            env.push(("DADA2RS_SHUFFLE_NO_CARRY", "1"));
        }
        if verify {
            env.push(("DADA2RS_RECONCILE_VERIFY", "1"));
        }
        super::pooled(out, threads, &env);
    }

    /// The carried and rebuilt maps must produce identical output.
    #[test]
    fn carried_compmax_matches_rebuild() {
        let tmp = scratch("carry");
        let (carried, rebuilt) = (tmp.join("carried"), tmp.join("rebuilt"));
        std::fs::create_dir_all(&carried).unwrap();
        std::fs::create_dir_all(&rebuilt).unwrap();

        run_pooled(&carried, true, "1", false);
        run_pooled(&rebuilt, false, "1", false);
        assert_same_outputs(&carried, &rebuilt, "carried vs rebuilt, single-threaded");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
    /// stored in, which is upstream of the map — so this exercises the carry against
    /// a different (still deterministic) input ordering.
    #[test]
    fn carried_compmax_matches_rebuild_threaded() {
        let tmp = scratch("carry_mt");
        let (carried, rebuilt) = (tmp.join("carried"), tmp.join("rebuilt"));
        std::fs::create_dir_all(&carried).unwrap();
        std::fs::create_dir_all(&rebuilt).unwrap();

        run_pooled(&carried, true, "4", false);
        run_pooled(&rebuilt, false, "4", false);
        assert_same_outputs(&carried, &rebuilt, "carried vs rebuilt, 4 threads");

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// The carry under the full-rescan invariant.
    ///
    /// This is the test that can actually fail for a carry-specific reason. The two
    /// above compare *outputs*, which a divergence only reaches if it survives to
    /// the partition; this asserts the carried `compmax`/`emax` equal a full rescan
    /// after every single reconcile, including the first one of each bud round —
    /// the one where the carry's freshly-relocated work lands.
    #[test]
    fn carried_compmax_survives_full_rescan_verify() {
        let tmp = scratch("carry_vfy");
        let out = tmp.join("carried");
        std::fs::create_dir_all(&out).unwrap();

        run_pooled(&out, true, "1", true);

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// The carry must also compose with the #132 move-pass prune disabled.
    ///
    /// The two interact: the carry decides what the reconcile has to repair, and the
    /// prune decides which clusters the *next* move pass visits based on what that
    /// reconcile marked dirty. A carry bug that under-marks would be masked by the
    /// unpruned scan, so both combinations need to agree with the baseline.
    #[test]
    fn carried_compmax_matches_rebuild_unpruned() {
        let tmp = scratch("carry_np");
        let (carried, rebuilt) = (tmp.join("carried"), tmp.join("rebuilt"));
        std::fs::create_dir_all(&carried).unwrap();
        std::fs::create_dir_all(&rebuilt).unwrap();

        let no_prune = ("DADA2RS_SHUFFLE_NO_PRUNE", "1");
        super::pooled(&carried, "1", &[no_prune]);
        super::pooled(
            &rebuilt,
            "1",
            &[no_prune, ("DADA2RS_SHUFFLE_NO_CARRY", "1")],
        );
        assert_same_outputs(&carried, &rebuilt, "carried vs rebuilt, unpruned move pass");

        let _ = std::fs::remove_dir_all(&tmp);
    }
}
