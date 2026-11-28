//! Stream existing serde JSON bytes into an identity without an unbounded buffer.
use collagene_model::Error;
use serde::Serialize;
use std::io::Write;

const MAX_BYTES: usize = 32 << 20;

/// Measure a domain's declared byte ceiling before hashing or cloning its DTO.
pub(crate) fn preflight(value: &impl Serialize, maximum: usize) -> Result<(), Error> {
    struct Count(usize);
    impl Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_sub(bytes.len())
                .ok_or_else(|| std::io::Error::other("bounded JSON admission"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(Count(maximum.min(MAX_BYTES)), value).map_err(|_| Error::TooLarge)
}

struct Hash {
    hash: blake3::Hasher,
    bytes: usize,
}

impl Write for Hash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded JSON identity"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

pub(crate) fn hash(value: &impl Serialize) -> Result<String, Error> {
    let mut writer = Hash {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut writer, value).map_err(|_| Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
