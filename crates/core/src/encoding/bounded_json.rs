//! Stream existing serde JSON bytes into an identity without an unbounded buffer.
use collagene_model::Error;
use serde::Serialize;
use std::io::Write;

const MAX_BYTES: usize = 32 << 20;

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
