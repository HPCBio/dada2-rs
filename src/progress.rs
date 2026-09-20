//! One bud round's progress line, assembled across phases and emitted once.
//!
//! DADA2's divisive progress output is a *partial* line built from several
//! writes as a bud round proceeds — the new cluster's id, the shuffle's `S`
//! per iteration, and finally the division that closed the round:
//!
//! ```text
//! New Cluster C5:C5LU:SS, Division (naive): Raw 42 from Bi 3, pA=1.2e-45
//! ```
//!
//! That is faithful to the R original, which writes the same fragments with
//! `Rprintf` (`Rmain.cpp:317`, `cluster.cpp:40,329,344,348`). It is also safe
//! there, because R's `dada()` denoises one sample per C++ call.
//!
//! dada2-rs runs up to `--sample-jobs` samples **concurrently in one process**,
//! and then the fragments collide: one sample's dangling partial line and
//! another's `eprintln!` land on the same output line. A filter anchored on
//! `^\[dada\]` then discards the collided line *and the message glued to it*,
//! which is how a 30-sample ITS2 run lost 26 of its 30 `resident Raw footprint`
//! lines with no trace (issue #172).
//!
//! Buffering the fragments and emitting the finished record with a single
//! `eprintln!` fixes that without changing a character of the text.

/// Accumulates one bud round's progress fragments, emitting them as one line.
///
/// Disabled (`None`) when not verbose, so the fragments are never formatted on
/// a production run. Dropping a record without calling [`Self::flush`] discards
/// it, which is what should happen when a round ends early.
#[derive(Debug, Default)]
pub struct ProgressRecord {
    buf: Option<String>,
    /// Prefixed to each record to say which sample produced it. Set only when
    /// samples run concurrently — with one `run_dada` in flight there is
    /// nothing to disambiguate, and the text then stays byte-identical to R's.
    tag: Option<String>,
}

impl ProgressRecord {
    /// A disabled record. Every push is a no-op and nothing is ever printed.
    pub fn disabled() -> Self {
        Self {
            buf: None,
            tag: None,
        }
    }

    /// An enabled record, optionally tagged with the sample that owns it.
    pub fn new(tag: Option<String>) -> Self {
        Self {
            buf: Some(String::new()),
            tag,
        }
    }

    /// `disabled()` or `new(tag)`, whichever `verbose` calls for.
    pub fn for_verbose(verbose: bool, tag: Option<String>) -> Self {
        if verbose {
            Self::new(tag)
        } else {
            Self::disabled()
        }
    }

    /// Whether anything is being collected. Callers use this to skip
    /// formatting work that would otherwise be thrown away.
    pub fn enabled(&self) -> bool {
        self.buf.is_some()
    }

    /// Append a fragment. Accepts `format_args!`, so call sites read exactly
    /// as the `eprint!` they replace.
    pub fn push(&mut self, args: std::fmt::Arguments<'_>) {
        if let Some(b) = self.buf.as_mut() {
            use std::fmt::Write;
            // Writing to a String is infallible.
            let _ = b.write_fmt(args);
        }
    }

    /// Emit the record as one atomic line and reset for the next round.
    ///
    /// A record holding nothing prints nothing: a round that produced no
    /// fragments should not leave a blank line behind.
    pub fn flush(&mut self) {
        let Some(b) = self.buf.as_mut() else { return };
        if b.is_empty() {
            return;
        }
        match &self.tag {
            Some(t) => eprintln!("[{t}] {b}"),
            None => eprintln!("{b}"),
        }
        b.clear();
    }
}

/// Append to a [`ProgressRecord`], mirroring `eprint!`'s call shape.
#[macro_export]
macro_rules! rec_print {
    ($rec:expr, $($arg:tt)*) => {
        $rec.push(format_args!($($arg)*))
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disabled_collects_nothing() {
        let mut r = ProgressRecord::disabled();
        assert!(!r.enabled());
        rec_print!(r, "New Cluster C{}:", 5);
        // Nothing buffered, so a flush cannot print.
        assert!(r.buf.is_none());
        r.flush();
    }

    #[test]
    fn fragments_accumulate_in_order() {
        let mut r = ProgressRecord::new(None);
        assert!(r.enabled());
        rec_print!(r, "New Cluster C{}:", 5);
        rec_print!(r, "C{}LU:", 5);
        rec_print!(r, "{}", "S".repeat(2));
        rec_print!(
            r,
            ", Division (naive): Raw {} from Bi {}, pA={:.2e}",
            42,
            3,
            1.2e-45
        );
        assert_eq!(
            r.buf.as_deref().unwrap(),
            "New Cluster C5:C5LU:SS, Division (naive): Raw 42 from Bi 3, pA=1.20e-45"
        );
    }

    /// The buffer must reset, or round 2 would repeat round 1's text.
    #[test]
    fn flush_clears_the_buffer() {
        let mut r = ProgressRecord::new(None);
        rec_print!(r, "first");
        r.flush();
        assert_eq!(r.buf.as_deref().unwrap(), "");
        rec_print!(r, "second");
        assert_eq!(r.buf.as_deref().unwrap(), "second");
    }

    /// An untagged record must be byte-identical to the R-derived text, so a
    /// serial or pooled run stays directly comparable with the original.
    #[test]
    fn tag_is_absent_unless_asked_for() {
        let mut untagged = ProgressRecord::new(None);
        rec_print!(untagged, "New Cluster C1:");
        assert!(!untagged.buf.as_deref().unwrap().starts_with('['));

        let mut tagged = ProgressRecord::new(Some("sam1F".to_string()));
        rec_print!(tagged, "New Cluster C1:");
        // The tag is applied at flush, not stored in the buffer.
        assert_eq!(tagged.buf.as_deref().unwrap(), "New Cluster C1:");
        assert_eq!(tagged.tag.as_deref(), Some("sam1F"));
    }

    #[test]
    fn for_verbose_selects_the_right_mode() {
        assert!(!ProgressRecord::for_verbose(false, None).enabled());
        assert!(ProgressRecord::for_verbose(true, None).enabled());
        assert!(
            ProgressRecord::for_verbose(false, Some("s".into()))
                .buf
                .is_none()
        );
    }
}
