use crate::{Error, Result};
use cap_fs_ext::{FollowSymlinks, OpenOptionsFollowExt};
use cap_std::{
    ambient_authority,
    fs::{Dir, OpenOptions},
};
use std::{
    io::{Read, Write},
    path::Path,
};
use subtle::ConstantTimeEq;

/// No Debug/Serialize implementation: secrets must not appear in logs or data.
pub struct AccessKey([u8; 32]);
impl AccessKey {
    /// Operator-selected file containing exactly 64 hex characters plus an optional
    /// final newline. Actual entropy comes from key generation, not this parser.
    pub fn read(path: &Path) -> Result<Self> {
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let name = path.file_name().ok_or(Error::Configuration("key file"))?;
        let directory = Dir::open_ambient_dir(parent, ambient_authority())?;
        let mut options = OpenOptions::new();
        options.read(true).follow(FollowSymlinks::No);
        let mut file = directory.open_with(name, &options)?.into_std();
        let metadata = file.metadata()?;
        if !metadata.is_file() || metadata.len() > 66 {
            return Err(Error::Configuration("key file"));
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            if metadata.mode() & 0o077 != 0 {
                return Err(Error::Configuration("key file must be owner-only"));
            }
        }
        let mut bytes = Vec::new();
        Read::by_ref(&mut file).take(67).read_to_end(&mut bytes)?;
        if bytes.ends_with(b"\r\n") {
            bytes.truncate(bytes.len() - 2);
        } else if bytes.ends_with(b"\n") {
            bytes.pop();
        }
        Ok(Self(decode_hex(&bytes).ok_or(Error::Configuration(
            "key file must contain 64 hex characters",
        ))?))
    }
    /// Explicit administration action. Never prints the key, overwrites a file,
    /// or generates a key as a side effect of starting the server.
    pub fn generate_file(path: &Path) -> Result<()> {
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let name = path.file_name().ok_or(Error::Configuration("key file"))?;
        let directory = Dir::open_ambient_dir(parent, ambient_authority())?;
        let mut random = [0u8; 32];
        getrandom::fill(&mut random).map_err(|_| Error::Entropy)?;
        let mut options = OpenOptions::new();
        options
            .write(true)
            .create_new(true)
            .follow(FollowSymlinks::No);
        #[cfg(unix)]
        {
            use cap_std::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        let mut file = directory.open_with(name, &options)?.into_std();
        let alphabet = b"0123456789abcdef";
        let mut text = [0u8; 65];
        for (i, b) in random.iter().enumerate() {
            text[i * 2] = alphabet[usize::from(b >> 4)];
            text[i * 2 + 1] = alphabet[usize::from(b & 15)];
        }
        text[64] = b'\n';
        file.write_all(&text)?;
        file.sync_all()?;
        Ok(())
    }
    pub(crate) fn same_secret(&self, other: &Self) -> bool {
        bool::from(self.0.ct_eq(&other.0))
    }
    pub(crate) fn authenticates(&self, value: &[u8]) -> bool {
        let Some(encoded) = value.strip_prefix(b"Bearer ") else {
            return false;
        };
        let Some(candidate) = decode_hex(encoded) else {
            return false;
        };
        bool::from(self.0.ct_eq(&candidate))
    }
}
fn decode_hex(encoded: &[u8]) -> Option<[u8; 32]> {
    if encoded.len() != 64 {
        return None;
    }
    let mut bytes = [0u8; 32];
    fn digit(value: u8) -> Option<u8> {
        match value {
            b'0'..=b'9' => Some(value - b'0'),
            b'a'..=b'f' => Some(value - b'a' + 10),
            b'A'..=b'F' => Some(value - b'A' + 10),
            _ => None,
        }
    }
    for (i, pair) in encoded.chunks_exact(2).enumerate() {
        bytes[i] = digit(pair[0])? * 16 + digit(pair[1])?;
    }
    Some(bytes)
}
