use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use serde::de::DeserializeOwned;

/// Open a JSON file and deserialize it, transparently decompressing gzip when
/// the path ends with `.gz` (e.g. `foo.json.gz`).
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
pub fn ntstr(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| nt_decode(b)).collect()
}

/// Print an alignment of two integer-encoded sequences to stderr.
/// Equivalent to C++ `align_print`.
pub fn align_print(al0: &[u8], al1: &[u8]) {
    assert_eq!(al0.len(), al1.len(), "alignment strands must have equal length");
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
pub fn err_print(err: &[[f64; 4]; 4]) {
    for (i, row) in err.iter().enumerate() {
        if i == 0 { eprint!("{{"); } else { eprint!(" "); }
        eprint!("{{");
        for (j, val) in row.iter().enumerate() {
            eprint!("{val:.6}");
            if j < 3 { eprint!(", "); }
        }
        if i < 3 { eprintln!("}},"); } else { eprintln!("}}}}"); }
    }
}
