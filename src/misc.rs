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

#[cfg(not(target_os = "linux"))]
pub fn cpu_allocation() -> (Option<usize>, Option<usize>) {
    (None, None)
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
        std::env::var("SLURM_JOB_ID").ok().as_deref(),
    )
}

/// Formatting half of [`cpu_allocation_repr`], split out so the Linux-only
/// output can be tested from any platform.
fn format_allocation(
    threads: usize,
    logical: Option<usize>,
    physical: Option<usize>,
    host: Option<&str>,
    slurm_job: Option<&str>,
) -> (String, Option<String>) {
    let mut desc = format!("{threads} threads");
    match (logical, physical) {
        (Some(l), Some(p)) => desc.push_str(&format!(" on {l} logical CPUs / {p} physical cores")),
        (Some(l), None) => desc.push_str(&format!(" on {l} logical CPUs (core topology unknown)")),
        _ => desc.push_str(" (allocation topology unavailable)"),
    }
    if let Some(h) = host {
        desc.push_str(&format!(" [host {h}"));
        if let Some(job) = slurm_job {
            desc.push_str(&format!(", SLURM job {job}"));
        }
        desc.push(']');
    }

    let warn = physical.filter(|&p| threads > p).map(|p| {
        format!(
            "WARNING: {threads} threads share {p} physical cores (SMT siblings), so \
             throughput may not exceed {p} cores' worth. Under SLURM request cores, \
             not tasks: --ntasks=1 --cpus-per-task=N --threads-per-core=1, and set \
             --threads from $SLURM_CPUS_PER_TASK (not $SLURM_NTASKS, which cannot \
             distinguish the two)."
        )
    });
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

    /// The two real allocations from issue #127, so the Linux output is
    /// pinned from any platform.
    #[test]
    fn format_allocation_matches_the_issue_127_cases() {
        // -n 24 on the 128-core node: Cpus_allowed_list 0-11,128-139
        let (desc, warn) =
            format_allocation(24, Some(24), Some(12), Some("compute-5-6"), Some("2352266"));
        assert_eq!(
            desc,
            "24 threads on 24 logical CPUs / 12 physical cores [host compute-5-6, SLURM job 2352266]"
        );
        let w = warn.expect("must warn: 24 threads on 12 cores");
        assert!(w.contains("24 threads share 12 physical cores"), "{w}");

        // -c 24 --threads-per-core=1: Cpus_allowed_list 0-23,128-151
        let (desc, warn) = format_allocation(24, Some(48), Some(24), Some("compute-5-6"), None);
        assert_eq!(
            desc,
            "24 threads on 48 logical CPUs / 24 physical cores [host compute-5-6]"
        );
        assert!(
            warn.is_none(),
            "must not warn when threads <= cores: {warn:?}"
        );
    }

    #[test]
    fn format_allocation_degrades_without_topology() {
        let (desc, warn) = format_allocation(8, None, None, None, None);
        assert_eq!(desc, "8 threads (allocation topology unavailable)");
        assert!(warn.is_none());
        let (desc, warn) = format_allocation(8, Some(16), None, Some("h"), None);
        assert_eq!(
            desc,
            "8 threads on 16 logical CPUs (core topology unknown) [host h]"
        );
        assert!(warn.is_none(), "unknown core count must not warn");
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
