use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use serde::Serialize;
use serde::de::DeserializeOwned;

/// Write `bytes` to `path`, gzip-compressing when the path ends in `.gz`
/// (otherwise written verbatim). Parent directories are created as needed.
///
/// Plain gzip (not bgzf): these JSON artifacts are read back whole via
/// [`read_all_maybe_gz`], so bgzf's blocked/seekable layout and multithreaded
/// compression buy nothing here — plain gzip is smaller and simpler, and the
/// output is read transparently by the existing `.gz`-sniffing reader.
pub fn write_maybe_gz(path: &Path, bytes: &[u8]) -> io::Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    let file = File::create(path)?;
    let is_gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    if is_gz {
        let mut w = GzEncoder::new(file, Compression::default());
        w.write_all(bytes)?;
        w.finish()?;
    } else {
        let mut w = BufWriter::new(file);
        w.write_all(bytes)?;
        w.flush()?;
    }
    Ok(())
}

/// Returns `true` when `path` is the stdin sentinel `-`.
fn is_stdin(path: &Path) -> bool {
    path.as_os_str() == "-"
}

/// A trait for call paths to add the path causing the error
pub trait WithPath<T> {
    fn with_path(self, path: &Path) -> io::Result<T>;
}

impl<T> WithPath<T> for io::Result<T> {
    fn with_path(self, path: &Path) -> io::Result<T> {
        self.map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
    }
}

/// Read the full contents of `path` into a `Vec<u8>`.  Treats `-` as stdin and
/// transparently decompresses gzip — by extension when reading a real file,
/// by sniffing the magic bytes (`1f 8b`) when reading stdin.
pub fn read_all_maybe_gz(path: &Path) -> io::Result<Vec<u8>> {
    if is_stdin(path) {
        let mut raw = Vec::new();
        io::stdin().lock().read_to_end(&mut raw)?;
        if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
            let mut out = Vec::new();
            MultiGzDecoder::new(raw.as_slice()).read_to_end(&mut out)?;
            Ok(out)
        } else {
            Ok(raw)
        }
    } else {
        let file = File::open(path)?;
        let is_gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
        let mut out = Vec::new();
        if is_gz {
            MultiGzDecoder::new(file).read_to_end(&mut out)?;
        } else {
            BufReader::new(file).read_to_end(&mut out)?;
        }
        Ok(out)
    }
}

/// Process peak resident-set size (high-water mark) in kibibytes, via
/// `getrusage(RUSAGE_SELF).ru_maxrss`.
///
/// `ru_maxrss` is monotonic (a high-water mark, never decreasing), so logging it
/// at successive phase boundaries reveals *which phase* drove the peak: the phase
/// after which it jumps owns that memory. Linux reports `ru_maxrss` in kB; macOS
/// in bytes — normalized to kB here, matching the benchmark harness.
pub fn peak_rss_kb() -> u64 {
    // SAFETY: getrusage with a zeroed rusage out-param is always sound.
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut ru) } != 0 {
        return 0;
    }
    let maxrss = ru.ru_maxrss as u64;
    if cfg!(target_os = "macos") {
        maxrss / 1024 // bytes -> kB
    } else {
        maxrss // already kB on Linux
    }
}

/// Format a path for use in error messages.  `-` becomes `<stdin>`.
fn display_input(path: &Path) -> String {
    if is_stdin(path) {
        "<stdin>".to_string()
    } else {
        path.display().to_string()
    }
}

/// Crate version, embedded in every JSON output as `dada2_rs_version`.
///
/// On a tagged release build (HEAD == `v<CARGO_PKG_VERSION>` or
/// `<CARGO_PKG_VERSION>`) this is the bare semver string; otherwise it is
/// `<CARGO_PKG_VERSION>-<short-sha>`. See `build.rs`.
pub const DADA2_RS_VERSION: &str = env!("DADA2_RS_VERSION_FULL");

/// Wraps a serializable output with the dada2-rs command name and version.
/// The two tag fields are emitted at the top of the resulting JSON object,
/// followed by the inner struct's fields (via `#[serde(flatten)]`).
#[derive(Serialize)]
pub struct Tagged<T: Serialize> {
    pub dada2_rs_command: &'static str,
    pub dada2_rs_version: &'static str,
    #[serde(flatten)]
    pub inner: T,
}

impl<T: Serialize> Tagged<T> {
    pub fn new(command: &'static str, inner: T) -> Self {
        Self {
            dada2_rs_command: command,
            dada2_rs_version: DADA2_RS_VERSION,
            inner,
        }
    }
}

/// Read a tagged JSON file and validate its `dada2_rs_command` is one of
/// `expected`.  Returns the inner payload on success.
///
/// The tag is checked first, against a `serde_json::Value` parse, so the
/// caller gets a clear "wrong command" error instead of a confusing
/// "missing field X" error when the file came from the wrong subcommand.
///
/// Errors with `InvalidData` if the tag is missing or mismatched.
/// Transparently decompresses gzip — by `.gz` extension for real files, or by
/// gzip-magic detection when reading stdin (`path == "-"`).
pub fn read_tagged_json<T: DeserializeOwned>(path: &Path, expected: &[&str]) -> io::Result<T> {
    let bytes = read_all_maybe_gz(path)?;
    let value: serde_json::Value = serde_json::from_slice(&bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let label = display_input(path);
    let cmd = value.get("dada2_rs_command").and_then(|v| v.as_str());
    match cmd {
        Some(c) if expected.contains(&c) => {}
        Some(c) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{label}: dada2_rs_command={c:?}, expected one of {expected:?}"),
            ));
        }
        None => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{label}: missing dada2_rs_command tag (expected one of {expected:?})"),
            ));
        }
    }

    serde_json::from_value(value)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{label}: {e}")))
}

/// Open a JSON file and deserialize it, transparently decompressing gzip when
/// the path ends with `.gz` (e.g. `foo.json.gz`).
#[allow(dead_code)]
pub fn read_json_file<T: DeserializeOwned>(path: &Path) -> io::Result<T> {
    let file = File::open(path)?;
    let is_gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    if is_gz {
        serde_json::from_reader(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        serde_json::from_reader(BufReader::new(file))
    }
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Read all records from a FASTA file, returning `(header, sequence)` pairs.
///
/// The header is the full line after `>`, trimmed of leading/trailing
/// whitespace.  The function transparently decompresses gzip when the path
/// ends with `.gz`.
pub fn read_fasta_records(path: &Path) -> io::Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(path)?;
    let is_gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    let reader: Box<dyn io::Read> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut cur_header: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    for line_result in BufReader::new(reader).lines() {
        let line = line_result?;
        let trimmed = line.trim_end();
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(header) = cur_header.take() {
                records.push((header, std::mem::take(&mut cur_seq)));
            }
            cur_header = Some(rest.trim().to_string());
        } else if !trimmed.is_empty() {
            cur_seq.extend_from_slice(trimmed.as_bytes());
        }
    }
    if let Some(header) = cur_header {
        records.push((header, cur_seq));
    }
    Ok(records)
}

/// Nucleotide integer encoding used throughout dada2: A=1, C=2, G=3, T=4, N=5.
/// Gap characters (b'-') pass through unchanged in both directions.
pub const NT_A: u8 = 1;
pub const NT_C: u8 = 2;
pub const NT_G: u8 = 3;
pub const NT_T: u8 = 4;
pub const NT_N: u8 = 5;

/// Encode one ASCII nucleotide byte to its integer representation.
/// Returns 0 for unrecognized characters.
pub fn nt_encode(b: u8) -> u8 {
    match b {
        b'A' | b'a' => NT_A,
        b'C' | b'c' => NT_C,
        b'G' | b'g' => NT_G,
        b'T' | b't' | b'U' | b'u' => NT_T,
        b'N' | b'n' => NT_N,
        b'-' => b'-',
        _ => 0,
    }
}

/// Decode one integer-encoded nucleotide back to its ASCII representation.
/// Returns b'?' for unrecognized values.
pub fn nt_decode(b: u8) -> u8 {
    match b {
        NT_A => b'A',
        NT_C => b'C',
        NT_G => b'G',
        NT_T => b'T',
        NT_N => b'N',
        b'-' => b'-',
        _ => b'?',
    }
}

/// Encode an ASCII nucleotide slice, returning a new `Vec<u8>`.
/// Equivalent to C++ `intstr`.
pub fn intstr(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| nt_encode(b)).collect()
}

/// Decode an integer-encoded slice back to ASCII, returning a new `Vec<u8>`.
/// Equivalent to C++ `ntstr`.
#[allow(dead_code)]
pub fn ntstr(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| nt_decode(b)).collect()
}

/// Print an alignment of two integer-encoded sequences to stderr.
/// Equivalent to C++ `align_print`.
#[allow(dead_code)]
pub fn align_print(al0: &[u8], al1: &[u8]) {
    assert_eq!(
        al0.len(),
        al1.len(),
        "alignment strands must have equal length"
    );
    eprintln!("{}", String::from_utf8_lossy(&ntstr(al0)));
    let mid: String = al0
        .iter()
        .zip(al1.iter())
        .map(|(a, b)| if a == b { '|' } else { ' ' })
        .collect();
    eprintln!("{mid}");
    eprintln!("{}", String::from_utf8_lossy(&ntstr(al1)));
}

/// Print a 4×4 error rate matrix to stderr.
/// Equivalent to C++ `err_print`.
#[allow(dead_code)]
pub fn err_print(err: &[[f64; 4]; 4]) {
    for (i, row) in err.iter().enumerate() {
        if i == 0 {
            eprint!("{{");
        } else {
            eprint!(" ");
        }
        eprint!("{{");
        for (j, val) in row.iter().enumerate() {
            eprint!("{val:.6}");
            if j < 3 {
                eprint!(", ");
            }
        }
        if i < 3 {
            eprintln!("}},");
        } else {
            eprintln!("}}}}");
        }
    }
}

// ---------------------------------------------------------------------------
// CPU allocation topology (issue #127)
// ---------------------------------------------------------------------------

/// Parse a Linux CPU list (`Cpus_allowed_list` / sysfs style) such as
/// `"0-11,128-139"` or `"0,2,4-7"` into individual CPU ids.
///
/// Returns `None` if any element fails to parse, so a malformed file degrades
/// to "unknown" rather than to a confidently wrong core count.
pub fn parse_cpu_list(s: &str) -> Option<Vec<usize>> {
    let mut out = Vec::new();
    for part in s.trim().split(',').filter(|p| !p.is_empty()) {
        match part.split_once('-') {
            Some((lo, hi)) => {
                let (lo, hi) = (lo.trim().parse().ok()?, hi.trim().parse::<usize>().ok()?);
                if hi < lo {
                    return None;
                }
                out.extend(lo..=hi);
            }
            None => out.push(part.trim().parse().ok()?),
        }
    }
    Some(out)
}

/// Logical CPUs this process may run on, and how many *physical* cores they
/// map to.
///
/// The second number is the one that matters and the one no environment
/// variable reports: a SLURM allocation of "24" can be 24 cores or 12 cores
/// with both SMT siblings of each, and `$SLURM_NTASKS`, `$SLURM_CPUS_ON_NODE`
/// and `$SLURM_TASKS_PER_NODE` all read `24` either way. Cores are counted as
/// distinct `(physical_package_id, core_id)` pairs from sysfs.
///
/// Returns `(logical, physical)`; `physical` is `None` when sysfs topology is
/// unavailable (containers, non-Linux, restricted mounts).
#[cfg(target_os = "linux")]
pub fn cpu_allocation() -> (Option<usize>, Option<usize>) {
    use std::collections::HashSet;

    let status = match std::fs::read_to_string("/proc/self/status") {
        Ok(s) => s,
        Err(_) => return (None, None),
    };
    let list = status
        .lines()
        .find_map(|l| l.strip_prefix("Cpus_allowed_list:"))
        .and_then(parse_cpu_list);
    let Some(cpus) = list else {
        return (None, None);
    };

    let mut cores = HashSet::new();
    for &c in &cpus {
        let base = format!("/sys/devices/system/cpu/cpu{c}/topology");
        let core = std::fs::read_to_string(format!("{base}/core_id"));
        let pkg = std::fs::read_to_string(format!("{base}/physical_package_id"));
        match (core, pkg) {
            (Ok(core), Ok(pkg)) => {
                cores.insert((pkg.trim().to_string(), core.trim().to_string()));
            }
            // Topology unreadable for even one CPU: report logical only rather
            // than an undercount that would read as a false SMT warning.
            _ => return (Some(cpus.len()), None),
        }
    }
    (Some(cpus.len()), Some(cores.len()))
}

/// Non-Linux fallback. The logical CPU count is available everywhere, but
/// there is no portable core-topology source, so physical cores are reported
/// as unknown rather than guessed — an undercount would fire a spurious SMT
/// warning, which is worse than saying nothing.
#[cfg(not(target_os = "linux"))]
pub fn cpu_allocation() -> (Option<usize>, Option<usize>) {
    (
        std::thread::available_parallelism().ok().map(|n| n.get()),
        None,
    )
}

/// One-line `--verbose` description of the CPU allocation the run actually got,
/// so archived logs are self-describing about it (issue #127).
///
/// Motivation: pooled benchmarks were run for months at `--threads 24` on
/// allocations that were 12 physical cores, which halved throughput invisibly.
/// Nothing in the environment revealed it — and the in-process
/// `map parallel efficiency` statistic cannot either, being
/// `busy / (map × nthreads)`: threads each running at half speed on shared SMT
/// siblings still report as fully efficient.
///
/// Returns the description plus an optional warning when `threads` exceeds the
/// physical cores available.
pub fn cpu_allocation_repr(threads: usize) -> (String, Option<String>) {
    let (logical, physical) = cpu_allocation();
    let host = std::fs::read_to_string("/proc/sys/kernel/hostname")
        .ok()
        .map(|h| h.trim().to_string())
        .filter(|h| !h.is_empty());
    format_allocation(
        threads,
        logical,
        physical,
        host.as_deref(),
        scheduler_job_id().as_deref(),
        scheduler_allocation(),
    )
}

/// CPUs the batch scheduler says it allocated on this node, with the scheduler's
/// name — or `None` when not running under a recognised one (a laptop, a
/// workstation, an unrecognised scheduler).
///
/// Covers the common HPC schedulers rather than SLURM alone, since `dada2-rs`
/// is also run on PBS/LSF/SGE sites and on ordinary machines.
fn scheduler_allocation() -> Option<(&'static str, usize)> {
    // SLURM. SLURM_JOB_CPUS_PER_NODE is authoritative but can be written in
    // run-length form ("48", "48(x2)", "24,48"); take the first count, which
    // is this node's. Fall back to cpus-per-task × ntasks.
    if std::env::var_os("SLURM_JOB_ID").is_some() {
        if let Ok(v) = std::env::var("SLURM_JOB_CPUS_PER_NODE") {
            let head: String = v.chars().take_while(|c| c.is_ascii_digit()).collect();
            if let Ok(n) = head.parse() {
                return Some(("SLURM", n));
            }
        }
        if let Some(cpt) = std::env::var("SLURM_CPUS_PER_TASK")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
        {
            let ntasks = std::env::var("SLURM_NTASKS")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(1);
            return Some(("SLURM", cpt * ntasks));
        }
    }
    // PBS Pro / OpenPBS / TORQUE: NCPUS is set from the select statement.
    if std::env::var_os("PBS_JOBID").is_some() {
        for k in ["NCPUS", "PBS_NP"] {
            if let Some(n) = std::env::var(k).ok().and_then(|v| v.parse().ok()) {
                return Some(("PBS", n));
            }
        }
    }
    // LSF.
    if std::env::var_os("LSB_JOBID").is_some()
        && let Some(n) = std::env::var("LSB_DJOB_NUMPROC")
            .ok()
            .and_then(|v| v.parse().ok())
    {
        return Some(("LSF", n));
    }
    // Grid Engine (SGE/UGE/OGE).
    if std::env::var_os("JOB_ID").is_some()
        && let Some(n) = std::env::var("NSLOTS").ok().and_then(|v| v.parse().ok())
    {
        return Some(("Grid Engine", n));
    }
    None
}

/// The scheduler's job id, for tying an archived log to an allocation.
fn scheduler_job_id() -> Option<String> {
    ["SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID", "JOB_ID"]
        .iter()
        .find_map(|k| std::env::var(k).ok())
        .filter(|v| !v.is_empty())
}

/// Formatting half of [`cpu_allocation_repr`], split out so the Linux-only
/// output can be tested from any platform.
fn format_allocation(
    threads: usize,
    logical: Option<usize>,
    physical: Option<usize>,
    host: Option<&str>,
    job_id: Option<&str>,
    allocation: Option<(&str, usize)>,
) -> (String, Option<String>) {
    let mut desc = format!("{threads} threads");
    match (logical, physical) {
        (Some(l), Some(p)) => desc.push_str(&format!(" on {l} logical CPUs / {p} physical cores")),
        (Some(l), None) => desc.push_str(&format!(" on {l} logical CPUs (core topology unknown)")),
        _ => desc.push_str(" (allocation topology unavailable)"),
    }
    if let Some((sched, n)) = allocation {
        desc.push_str(&format!("; {sched} allocated {n} CPUs"));
    }
    if let Some(h) = host {
        desc.push_str(&format!(" [host {h}"));
        if let Some(job) = job_id {
            desc.push_str(&format!(", job {job}"));
        }
        desc.push(']');
    }

    // Two distinct hazards. Both messages are written to be true off a batch
    // scheduler as well as on one — a laptop can oversubscribe SMT siblings
    // just as a compute node can — so scheduler-specific advice is appended
    // only when a scheduler was actually detected.
    let sched_hint = |sched: &str| match sched {
        "SLURM" => " Under SLURM, request cores rather than tasks: \
                    --ntasks=1 --cpus-per-task=N --threads-per-core=1, and set --threads \
                    from $SLURM_CPUS_PER_TASK (not $SLURM_NTASKS, which cannot \
                    distinguish the two)."
            .to_string(),
        "PBS" => " Under PBS, select cores explicitly (e.g. -l select=1:ncpus=N) and set \
                  --threads from $NCPUS."
            .to_string(),
        other => format!(" Set --threads from the CPU count {other} allocated, not the node's."),
    };

    let warn = match (physical, allocation) {
        // More threads than cores this process may use: SMT oversubscription.
        (Some(p), alloc) if threads > p => {
            let mut w = format!(
                "WARNING: {threads} threads share {p} physical cores (SMT siblings), so \
                 throughput may not exceed {p} cores' worth."
            );
            if let Some((sched, _)) = alloc {
                w.push_str(&sched_hint(sched));
            }
            Some(w)
        }
        // The scheduler granted fewer CPUs than the process can actually reach.
        // Cores are only enforced by affinity or a cpuset cgroup when the site
        // configures it or the step is launched through the scheduler's
        // launcher; a batch script invoking the binary directly frequently runs
        // unconfined. Reported, not corrected: running wider is usually faster
        // for the job and worse for the node, so it is the operator's call.
        (Some(p), Some((sched, n))) if p > n => Some(format!(
            "NOTE: {sched} allocated {n} CPUs but this process is not confined to them \
             ({p} physical cores are reachable). Threads may land on cores allocated to \
             other jobs, and timings are not attributable to the allocation."
        )),
        _ => None,
    };
    (desc, warn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_cpu_list_ranges_and_singletons() {
        // The two real allocations from issue #127: 12 cores presented as 24
        // logical CPUs, and the corrected 24-core request.
        assert_eq!(parse_cpu_list("0-11,128-139").unwrap().len(), 24);
        assert_eq!(parse_cpu_list("0-23,128-151").unwrap().len(), 48);
        assert_eq!(parse_cpu_list("0-31").unwrap().len(), 32);
        assert_eq!(parse_cpu_list("0,2,4-7").unwrap(), vec![0, 2, 4, 5, 6, 7]);
        assert_eq!(parse_cpu_list("  3  ").unwrap(), vec![3]);
        assert_eq!(parse_cpu_list("").unwrap(), Vec::<usize>::new());
    }

    #[test]
    fn parse_cpu_list_rejects_malformed() {
        // Must return None, never a partial list: an undercount would surface
        // as a spurious SMT warning.
        assert!(parse_cpu_list("0-").is_none());
        assert!(parse_cpu_list("7-3").is_none());
        assert!(parse_cpu_list("abc").is_none());
        assert!(parse_cpu_list("0-11,x").is_none());
    }

    /// The real allocations from issue #127, so the Linux output is pinned
    /// from any platform.
    #[test]
    fn format_allocation_matches_the_issue_127_cases() {
        // Interactive `srun --pty` under -n 24: bound to 12 cores.
        let (desc, warn) = format_allocation(
            24,
            Some(24),
            Some(12),
            Some("compute-5-6"),
            Some("2352266"),
            Some(("SLURM", 24)),
        );
        assert_eq!(
            desc,
            "24 threads on 24 logical CPUs / 12 physical cores; SLURM allocated 24 CPUs \
             [host compute-5-6, job 2352266]"
        );
        let w = warn.unwrap();
        assert!(w.contains("24 threads share 12 physical cores"), "{w}");
        assert!(w.contains("$SLURM_CPUS_PER_TASK"), "{w}");

        // Interactive under -c 24 --threads-per-core=1: 24 cores, both siblings.
        let (desc, warn) = format_allocation(
            24,
            Some(48),
            Some(24),
            Some("compute-5-6"),
            None,
            Some(("SLURM", 24)),
        );
        assert!(desc.starts_with("24 threads on 48 logical CPUs / 24 physical cores"));
        assert!(warn.is_none(), "threads <= cores must not warn: {warn:?}");

        // The batch case: unconfined, so the whole node is reachable. This is
        // what the production runs actually look like, and it must be flagged
        // as *unattributable* rather than reported as a 128-core allocation.
        let (desc, warn) = format_allocation(
            48,
            Some(256),
            Some(128),
            Some("compute-5-6"),
            Some("2352276"),
            Some(("SLURM", 48)),
        );
        assert_eq!(
            desc,
            "48 threads on 256 logical CPUs / 128 physical cores; SLURM allocated 48 CPUs \
             [host compute-5-6, job 2352276]"
        );
        let w = warn.expect("unconfined job must be flagged");
        assert!(w.contains("SLURM allocated 48 CPUs"), "{w}");
        assert!(w.contains("not confined"), "{w}");
    }

    /// Off a batch scheduler the output must stay generic: no allocation
    /// clause, no scheduler-specific advice, and no unattributable-allocation
    /// NOTE (there is no allocation to be unconfined from).
    #[test]
    fn format_allocation_is_generic_without_a_scheduler() {
        // Laptop, sane thread count.
        let (desc, warn) = format_allocation(8, Some(16), Some(8), None, None, None);
        assert_eq!(desc, "8 threads on 16 logical CPUs / 8 physical cores");
        assert!(warn.is_none(), "{warn:?}");

        // Laptop oversubscribing SMT: warn, but with no scheduler advice.
        let (_, warn) = format_allocation(16, Some(16), Some(8), None, None, None);
        let w = warn.expect("must warn when threads exceed cores");
        assert!(w.contains("16 threads share 8 physical cores"), "{w}");
        assert!(
            !w.contains("SLURM"),
            "leaked SLURM advice off-scheduler: {w}"
        );
        assert!(!w.contains("--cpus-per-task"), "leaked SLURM advice: {w}");

        // A wide machine with a small thread count is not "unconfined".
        let (_, warn) = format_allocation(8, Some(256), Some(128), None, None, None);
        assert!(warn.is_none(), "{warn:?}");
    }

    #[test]
    fn format_allocation_handles_other_schedulers() {
        let (desc, warn) = format_allocation(
            64,
            Some(64),
            Some(32),
            Some("nid001"),
            Some("12345.pbs"),
            Some(("PBS", 64)),
        );
        assert!(desc.contains("PBS allocated 64 CPUs"), "{desc}");
        assert!(desc.contains("[host nid001, job 12345.pbs]"), "{desc}");
        let w = warn.unwrap();
        assert!(w.contains("$NCPUS"), "{w}");
        assert!(!w.contains("SLURM"), "{w}");

        // An unrecognised scheduler still gets a usable, non-SLURM hint.
        let (_, warn) = format_allocation(8, Some(8), Some(4), None, None, Some(("LSF", 8)));
        let w = warn.unwrap();
        assert!(w.contains("LSF"), "{w}");
        assert!(!w.contains("--cpus-per-task"), "{w}");
    }

    #[test]
    fn cpu_allocation_repr_is_always_printable() {
        // Must degrade gracefully wherever topology is unreadable (macOS,
        // containers) rather than panicking or omitting the thread count.
        let (desc, _) = cpu_allocation_repr(8);
        assert!(desc.starts_with("8 threads"), "unexpected: {desc}");
    }

    #[test]
    fn cpu_allocation_repr_warns_only_when_oversubscribed() {
        let (_, warn) = cpu_allocation_repr(1);
        if let (_, Some(p)) = cpu_allocation() {
            assert!(p >= 1);
            // One thread can never oversubscribe.
            assert!(warn.is_none(), "spurious warning: {warn:?}");
            let (_, warn_many) = cpu_allocation_repr(p + 1);
            assert!(warn_many.is_some(), "missing warning at {} threads", p + 1);
            assert!(warn_many.unwrap().contains("SLURM_CPUS_PER_TASK"));
        }
    }
}
