use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use serde::Serialize;
use serde::de::DeserializeOwned;

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
/// Transparently decompresses gzip when the path ends in `.gz`.
pub fn read_tagged_json<T: DeserializeOwned>(path: &Path, expected: &[&str]) -> io::Result<T> {
    let file = File::open(path)?;
    let is_gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    let value: serde_json::Value = if is_gz {
        serde_json::from_reader(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        serde_json::from_reader(BufReader::new(file))
    }
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let cmd = value.get("dada2_rs_command").and_then(|v| v.as_str());
    match cmd {
        Some(c) if expected.contains(&c) => {}
        Some(c) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: dada2_rs_command={c:?}, expected one of {expected:?}",
                    path.display()
                ),
            ));
        }
        None => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: missing dada2_rs_command tag (expected one of {expected:?})",
                    path.display()
                ),
            ));
        }
    }

    serde_json::from_value(value).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: {e}", path.display()),
        )
    })
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
