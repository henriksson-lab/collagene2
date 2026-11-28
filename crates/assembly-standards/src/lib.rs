//! Declarative standard grammars and deterministic part-sequence linting.
//! No cleavage/assembly simulator, host inference, network or process execution.
pub mod assessment;
pub mod grammar;
pub mod profile;

use std::fmt;
pub const MAX_PROFILE_BYTES: usize = 128 * 1024;
pub const MAX_REPORT_BYTES: usize = 1024 * 1024;
pub const ENGINE_VERSION: &str = "assembly-standard-assessment-1";
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    UnsupportedProfile,
    UnsupportedVersion,
    Invalid(&'static str),
    CatalogMismatch,
    Limit(&'static str),
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedProfile => f.write_str("unsupported assembly-standard profile"),
            Self::UnsupportedVersion => f.write_str("unsupported assembly-standard schema version"),
            Self::Invalid(s) => write!(f, "invalid assembly-standard input: {s}"),
            Self::CatalogMismatch => f.write_str("assembly-standard catalog pin/geometry mismatch"),
            Self::Limit(s) => write!(f, "assembly-standard resource bound: {s}"),
        }
    }
}
impl std::error::Error for Error {}

pub(crate) fn require(ok: bool, message: &'static str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(Error::Invalid(message))
    }
}
pub(crate) fn text(s: &str, max: usize) -> bool {
    !s.trim().is_empty() && s.len() <= max && !s.chars().any(char::is_control)
}
pub(crate) fn concrete(s: &str, max: usize) -> bool {
    !s.is_empty() && s.len() <= max && s.bytes().all(|b| b"ACGT".contains(&b))
}
pub(crate) fn mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
pub(crate) fn bounded_json<T: serde::Serialize>(value: &T, max: usize) -> Result<Vec<u8>> {
    struct Sink {
        bytes: Vec<u8>,
        max: usize,
    }
    impl std::io::Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() > self.max.saturating_sub(self.bytes.len()) {
                return Err(std::io::Error::other("assembly-standard output budget"));
            }
            self.bytes.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut sink = Sink {
        bytes: Vec::new(),
        max,
    };
    serde_json::to_writer(&mut sink, value).map_err(|_| Error::Limit("serialized bytes"))?;
    Ok(sink.bytes)
}
