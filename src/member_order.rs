//! Experimental member-list ordering for the reproducibility-floor experiment
//! (issue #157). **Changes results; not R-equivalent.**
//!
//! `b_bud` breaks ties on position within a cluster's member list, and that
//! order is an artifact of insertion and `swap_remove`. Where the abundance
//! p-value underflows to zero, position is the operative comparator (see
//! `docs/findings/r-parity-floor-and-ceiling.md`). This module lets one binary
//! run the same input under different, deliberate member orders, so the spread
//! of results measures that floor.
//!
//! Selected by `DADA2RS_MEMBER_ORDER`:
//!
//! - `insertion` (default, or unset): no-op. Byte-identical to not having this
//!   module.
//! - `sorted`: members ascending by raw index (derep order), centre first.
//! - `shuffle:<seed>`: a seeded random order, centre first.
//!
//! A cluster is reordered only when its membership changed since the last bud
//! round (`update_e`), so order stays arbitrary but stable, as it is in any
//! real implementation. Those are exactly the clusters `b_p_update` reprices
//! and whose cached bud candidate (stored by position) it rebuilds, so the
//! cache stays consistent without extra invalidation.
//!
//! The centre is moved to position 0 in every non-default arm: `b_bud` never
//! considers position 0, and a centre elsewhere would reproduce #219/#239
//! instead of measuring tie sensitivity.

use std::sync::OnceLock;

use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;

use crate::containers::B;

/// The member-order arm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MemberOrder {
    Insertion,
    Sorted,
    Shuffle(u64),
}

impl MemberOrder {
    /// Parse a `DADA2RS_MEMBER_ORDER` value.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "" | "insertion" => Ok(Self::Insertion),
            "sorted" => Ok(Self::Sorted),
            v => match v.strip_prefix("shuffle:").map(str::parse::<u64>) {
                Some(Ok(seed)) => Ok(Self::Shuffle(seed)),
                _ => Err(format!(
                    "DADA2RS_MEMBER_ORDER={v:?} is not recognised; expected \
                     insertion, sorted or shuffle:<seed>"
                )),
            },
        }
    }

    /// The value as it would be written in the environment.
    pub fn label(self) -> String {
        match self {
            Self::Insertion => "insertion".into(),
            Self::Sorted => "sorted".into(),
            Self::Shuffle(seed) => format!("shuffle:{seed}"),
        }
    }
}

/// The resolved arm, read once per process.
///
/// An unparseable value is fatal rather than falling back to `insertion`: a
/// mistyped arm silently running the baseline is the failure `gates` exists to
/// prevent, and here it would make an experiment compare the baseline with
/// itself.
pub fn member_order() -> MemberOrder {
    static VALUE: OnceLock<MemberOrder> = OnceLock::new();
    *VALUE.get_or_init(|| match std::env::var("DADA2RS_MEMBER_ORDER") {
        Ok(v) => MemberOrder::parse(&v).unwrap_or_else(|e| panic!("{e}")),
        Err(_) => MemberOrder::Insertion,
    })
}

/// Applies the arm to a partition, one bud round at a time.
pub struct MemberOrderer {
    arm: MemberOrder,
    rng: Option<SmallRng>,
}

impl MemberOrderer {
    pub fn new(arm: MemberOrder) -> Self {
        let rng = match arm {
            MemberOrder::Shuffle(seed) => Some(SmallRng::seed_from_u64(seed)),
            _ => None,
        };
        Self { arm, rng }
    }

    /// Reorder the member list of every cluster flagged `update_e`. Call
    /// immediately before `b_p_update`. Returns the number of clusters
    /// reordered (always 0 for `insertion`).
    pub fn apply(&mut self, b: &mut B) -> usize {
        if self.arm == MemberOrder::Insertion {
            return 0;
        }
        let mut n = 0;
        // Ascending cluster order, so one RNG stream gives a deterministic run.
        for bi in b.clusters.iter_mut().filter(|bi| bi.update_e) {
            if let Some(center) = bi.center
                && let Some(pos) = bi.raws.iter().position(|&r| r == center)
            {
                bi.raws.swap(0, pos);
            }
            let rest = bi.raws.get_mut(1..).unwrap_or_default();
            match self.arm {
                MemberOrder::Sorted => rest.sort_unstable(),
                MemberOrder::Shuffle(_) => {
                    rest.shuffle(self.rng.as_mut().expect("shuffle arm has an rng"))
                }
                MemberOrder::Insertion => unreachable!(),
            }
            n += 1;
        }
        n
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::containers::Raw;

    /// A one-cluster partition of `n` raws, centred off position 0 so pinning
    /// has something to do.
    fn partition(n: usize, center: usize) -> B {
        let raws = (0..n)
            .map(|i| {
                Raw::new(
                    vec![1, 2, 3, (i % 4) as u8 + 1],
                    None,
                    (n - i) as u32,
                    false,
                )
            })
            .collect();
        let mut b = B::new(raws, 1e-40, 1e-4, false);
        b.clusters[0].center = Some(center);
        b
    }

    #[test]
    fn non_default_arms_pin_the_centre_and_keep_the_members() {
        for arm in [MemberOrder::Sorted, MemberOrder::Shuffle(3)] {
            let mut b = partition(40, 17);
            b.clusters[0].update_e = true;
            assert_eq!(MemberOrderer::new(arm).apply(&mut b), 1);
            let raws = &b.clusters[0].raws;
            assert_eq!(raws[0], 17, "{arm:?}: centre must be at position 0");
            let mut got = raws.clone();
            got.sort_unstable();
            assert_eq!(got, (0..40).collect::<Vec<_>>(), "{arm:?}: members changed");
            if arm == MemberOrder::Sorted {
                assert!(raws[1..].windows(2).all(|w| w[0] < w[1]));
            }
        }
    }

    #[test]
    fn shuffle_is_seeded_and_insertion_is_a_no_op() {
        let order = |arm| {
            let mut b = partition(40, 0);
            b.clusters[0].update_e = true;
            MemberOrderer::new(arm).apply(&mut b);
            b.clusters[0].raws.clone()
        };
        assert_eq!(
            order(MemberOrder::Shuffle(3)),
            order(MemberOrder::Shuffle(3))
        );
        assert_ne!(
            order(MemberOrder::Shuffle(3)),
            order(MemberOrder::Shuffle(4))
        );
        assert_eq!(order(MemberOrder::Insertion), (0..40).collect::<Vec<_>>());
    }

    #[test]
    fn clean_clusters_are_left_alone() {
        let mut b = partition(40, 17);
        b.clusters[0].update_e = false;
        let before = b.clusters[0].raws.clone();
        assert_eq!(MemberOrderer::new(MemberOrder::Shuffle(1)).apply(&mut b), 0);
        assert_eq!(b.clusters[0].raws, before);
    }

    #[test]
    fn parse_accepts_the_three_arms_and_rejects_the_rest() {
        assert_eq!(MemberOrder::parse(""), Ok(MemberOrder::Insertion));
        assert_eq!(MemberOrder::parse("insertion"), Ok(MemberOrder::Insertion));
        assert_eq!(MemberOrder::parse("sorted"), Ok(MemberOrder::Sorted));
        assert_eq!(MemberOrder::parse("shuffle:7"), Ok(MemberOrder::Shuffle(7)));
        for bad in ["sort", "shuffle", "shuffle:", "shuffle:x", "random:1"] {
            assert!(
                MemberOrder::parse(bad).is_err(),
                "{bad:?} should be rejected"
            );
        }
        assert_eq!(MemberOrder::Shuffle(7).label(), "shuffle:7");
    }
}
