//! Read preparation: `remove-primers`, `filter-and-trim`.

use std::io;

use dada2_rs::{cli, filter_trim, misc, remove_primers};

use filter_trim::{
    FilterParams, PairedFiles, WriteOptions, filter_paired, filter_single, read_fasta_first_seq,
};
use misc::Tagged;
use remove_primers::{RemovePrimersParams, iupac_reverse_complement, remove_primers};
use serde::Serialize;

use super::common::fastq_stem;

pub(crate) fn run_remove_primers(args: cli::RemovePrimersArgs) -> io::Result<()> {
    let cli::RemovePrimersArgs {
        input,
        fout,
        sample_name,
        primer_fwd,
        primer_rev,
        rc_primer_rev,
        max_mismatch,
        allow_indels,
        trim_fwd,
        trim_rev,
        orient,
        compress,
        threads,
        trunc_q,
        trunc_len,
        trim_left,
        trim_right,
        max_len,
        min_len,
        max_n,
        min_q,
        max_ee,
        phix_genome,
        rm_lowcomplex,
        phred_offset,
        output,
        compact,
        verbose,
    } = args;
    if allow_indels && verbose {
        eprintln!("[remove-primers] indel mode enabled — expect ~4× slower matching");
    }
    let primer_rev_bytes = primer_rev.map(|s| {
        let b = s.into_bytes();
        if rc_primer_rev {
            iupac_reverse_complement(&b)
        } else {
            b
        }
    });
    let filter_params = if trunc_q.is_some()
        || trunc_len.is_some()
        || trim_left.is_some()
        || trim_right.is_some()
        || max_len.is_some()
        || min_len.is_some()
        || max_n.is_some()
        || min_q.is_some()
        || max_ee.is_some()
        || phix_genome.is_some()
        || rm_lowcomplex.is_some()
    {
        let phix_seq: Option<Vec<u8>> = phix_genome
            .as_deref()
            .map(read_fasta_first_seq)
            .transpose()?;
        Some(FilterParams {
            trunc_q: trunc_q.unwrap_or(0),
            trunc_len: trunc_len.unwrap_or(0),
            trim_left: trim_left.unwrap_or(0),
            trim_right: trim_right.unwrap_or(0),
            max_len: max_len.unwrap_or(0),
            min_len: min_len.unwrap_or(0),
            max_n: max_n.unwrap_or(usize::MAX),
            min_q: min_q.unwrap_or(0),
            max_ee: max_ee.unwrap_or(f64::INFINITY),
            phix_genome: phix_seq,
            rm_lowcomplex: rm_lowcomplex.unwrap_or(0.0),
            phred_offset,
        })
    } else {
        None
    };
    // Validate filter params before processing.
    if let Some(ref fp) = filter_params
        && fp.max_len > 0
        && fp.min_len > fp.max_len
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "--min-len ({}) is greater than --max-len ({}); no read can satisfy both",
                fp.min_len, fp.max_len
            ),
        ));
    }
    let params = RemovePrimersParams {
        primer_fwd: primer_fwd.into_bytes(),
        primer_rev: primer_rev_bytes,
        max_mismatch,
        allow_indels,
        trim_fwd,
        trim_rev,
        orient,
        filter_params,
    };
    let sample = sample_name.unwrap_or_else(|| fastq_stem(&input));
    let stats = remove_primers(&input, &fout, &params, compress, threads, verbose)?;

    #[derive(Serialize)]
    struct RemovePrimersOutput {
        sample: String,
        reads_in: u64,
        reads_out: u64,
        reads_reoriented: u64,
        #[serde(skip_serializing_if = "is_zero")]
        reads_filter_fail: u64,
    }
    fn is_zero(v: &u64) -> bool {
        *v == 0
    }
    let tagged = Tagged::new(
        "remove-primers",
        RemovePrimersOutput {
            sample,
            reads_in: stats.reads_in,
            reads_out: stats.reads_out,
            reads_reoriented: stats.reads_reoriented,
            reads_filter_fail: stats.reads_filter_fail,
        },
    );
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_filter_and_trim(args: cli::FilterAndTrimArgs) -> io::Result<()> {
    let cli::FilterAndTrimArgs {
        fwd,
        filt,
        rev,
        filt_rev,
        sample_name,
        compress,
        threads,
        trunc_q,
        trunc_len,
        trim_left,
        trim_right,
        max_len,
        min_len,
        max_n,
        min_q,
        max_ee,
        phix_genome,
        rm_lowcomplex,
        phred_offset,
        output,
        compact,
        verbose,
    } = args;
    // ---- Validate paired-end files ----
    if rev.is_some() {
        filt_rev.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "--filt-rev is required when --rev is provided",
            )
        })?;
    }

    // ---- Helper: expand a 1-or-2-element Vec into (fwd_val, rev_val) ----
    macro_rules! pair {
        ($v:expr, $default:expr) => {{
            let v = &$v;
            if v.is_empty() {
                ($default, $default)
            } else if v.len() == 1 {
                (v[0], v[0])
            } else {
                (v[0], v[1])
            }
        }};
    }

    let (trunc_q_f, trunc_q_r) = pair!(trunc_q, 2u8);
    let (trunc_len_f, trunc_len_r) = pair!(trunc_len, 0usize);
    let (trim_left_f, trim_left_r) = pair!(trim_left, 0usize);
    let (trim_right_f, trim_right_r) = pair!(trim_right, 0usize);
    let (max_len_f, max_len_r) = pair!(max_len, 0usize);
    let (min_len_f, min_len_r) = pair!(min_len, 20usize);
    let (max_ee_f, max_ee_r) = if max_ee.is_empty() {
        (f64::INFINITY, f64::INFINITY)
    } else if max_ee.len() == 1 {
        (max_ee[0], max_ee[0])
    } else {
        (max_ee[0], max_ee[1])
    };
    let (rm_lowcomplex_f, rm_lowcomplex_r) = pair!(rm_lowcomplex, 0.0f64);

    let phix_seq: Option<Vec<u8>> = phix_genome
        .as_deref()
        .map(read_fasta_first_seq)
        .transpose()?;

    let make_params = |tq, tl, trl, trr, ml, mnl, ee, rlc| FilterParams {
        trunc_q: tq,
        trunc_len: tl,
        trim_left: trl,
        trim_right: trr,
        max_len: ml,
        min_len: mnl,
        max_n,
        min_q,
        max_ee: ee,
        phix_genome: phix_seq.clone(),
        rm_lowcomplex: rlc,
        phred_offset,
    };

    let params_fwd = make_params(
        trunc_q_f,
        trunc_len_f,
        trim_left_f,
        trim_right_f,
        max_len_f,
        min_len_f,
        max_ee_f,
        rm_lowcomplex_f,
    );
    let params_rev = make_params(
        trunc_q_r,
        trunc_len_r,
        trim_left_r,
        trim_right_r,
        max_len_r,
        min_len_r,
        max_ee_r,
        rm_lowcomplex_r,
    );

    let sample = sample_name.unwrap_or_else(|| fastq_stem(&fwd));
    let opts = WriteOptions {
        compress,
        threads,
        verbose,
    };

    let stats = if let (Some(rev_in), Some(rev_out)) = (rev, filt_rev) {
        filter_paired(
            &PairedFiles {
                fwd_in: &fwd,
                rev_in: &rev_in,
                fwd_out: &filt,
                rev_out: &rev_out,
            },
            &params_fwd,
            &params_rev,
            opts,
        )?
    } else {
        filter_single(&fwd, &filt, &params_fwd, opts)?
    };

    #[derive(Serialize)]
    struct FilterAndTrimOutput {
        sample: String,
        reads_in: u64,
        reads_out: u64,
    }
    let tagged = Tagged::new(
        "filter-and-trim",
        FilterAndTrimOutput {
            sample,
            reads_in: stats.reads_in,
            reads_out: stats.reads_out,
        },
    );
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}
