//! Bounded, checksummed current binary envelopes.
use super::schema::{CANDIDATE_SCHEMA_VERSION, WORKER_SCHEMA_VERSION};
use crate::formats::sequence::MAX_IMPORT_BYTES;
use collagene_model::Error;
use serde::{Deserialize, Serialize};
// Current selectors whose accepted payloads can contain complete workspace archives.
// Stateless and session dispatch must admit the same envelope before routing.
pub(super) const WORKSPACE_REQUEST_SELECTORS: &[u32] = &[
    19, 20, 42, 43, 69, 72, 82, 83, 84, 85, 96, 100, 103, 104, 105, 106, 107, 122, 126, 128, 129,
    132, 133, 140, 145, 146, 150, 153, 154, 156, 157, 158, 159, 160, 161, 163, 164, 165, 166, 167,
    168, 169, 171, 172, 173, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188,
    189, 190, 191, 197, 212,
];

fn current_version(magic: &[u8; 4]) -> u16 {
    match magic {
        b"CGQ\0" | b"CGR\0" => WORKER_SCHEMA_VERSION,
        b"CGP\0" => CANDIDATE_SCHEMA_VERSION,
        _ => 1,
    }
}
pub(super) const HEADER: usize = 42;

pub const MAX_WORKSPACE_MESSAGE_BYTES: usize = 128 * 1024 * 1024;
struct BoundedPayload {
    bytes: Vec<u8>,
    limit: usize,
}
impl postcard::ser_flavors::Flavor for BoundedPayload {
    type Output = Vec<u8>;
    fn try_push(&mut self, byte: u8) -> Result<(), postcard::Error> {
        if self.bytes.len() >= self.limit {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes.push(byte);
        Ok(())
    }
    fn try_extend(&mut self, bytes: &[u8]) -> Result<(), postcard::Error> {
        if bytes.len() > self.limit.saturating_sub(self.bytes.len()) {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }
    fn finalize(self) -> Result<Self::Output, postcard::Error> {
        Ok(self.bytes)
    }
}
pub(crate) fn encode<T: Serialize>(magic: &[u8; 4], value: &T) -> Result<Vec<u8>, Error> {
    encode_with_limit(magic, value, MAX_IMPORT_BYTES)
}
pub(super) fn encode_with_limit<T: Serialize>(
    magic: &[u8; 4],
    value: &T,
    maximum: usize,
) -> Result<Vec<u8>, Error> {
    let payload = postcard::serialize_with_flavor(
        value,
        BoundedPayload {
            bytes: Vec::new(),
            limit: maximum,
        },
    )
    .map_err(|error| {
        if matches!(error, postcard::Error::SerializeBufferFull) {
            Error::TooLarge
        } else {
            Error::Format(error.to_string())
        }
    })?;
    let mut bytes = Vec::with_capacity(HEADER + payload.len());
    bytes.extend_from_slice(magic);
    bytes.extend_from_slice(&current_version(magic).to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    Ok(bytes)
}
fn envelope_info(
    magic: &[u8; 4],
    bytes: &[u8],
    frame_length: usize,
    maximum: usize,
) -> Result<(u16, usize), Error> {
    if bytes.len() < HEADER || &bytes[..4] != magic {
        return Err(Error::Format(
            "not a supported Collagene binary file/message".into(),
        ));
    }
    let version = u16::from_le_bytes([bytes[4], bytes[5]]);
    if version != current_version(magic) {
        return Err(Error::Unsupported(format!(
            "binary schema version {version}"
        )));
    }
    let len = u32::from_le_bytes(bytes[6..10].try_into().unwrap()) as usize;
    if len > maximum {
        return Err(Error::TooLarge);
    }
    if bytes.len() > frame_length || HEADER.checked_add(len) != Some(frame_length) {
        return Err(Error::Format("binary length mismatch".into()));
    }
    Ok((version, len))
}
pub(crate) fn decode<T: for<'de> Deserialize<'de>>(
    magic: &[u8; 4],
    bytes: &[u8],
) -> Result<T, Error> {
    decode_with_limit(magic, bytes, MAX_IMPORT_BYTES)
}
pub(super) fn decode_with_limit<T: for<'de> Deserialize<'de>>(
    magic: &[u8; 4],
    bytes: &[u8],
    maximum: usize,
) -> Result<T, Error> {
    envelope_info(magic, bytes, bytes.len(), maximum)?;
    if blake3::hash(&bytes[HEADER..]).as_bytes() != &bytes[10..42] {
        return Err(Error::Format("binary checksum mismatch".into()));
    }
    let outcome = collagene_model::decode::with_native_frame_budget(maximum, || {
        postcard::take_from_bytes(&bytes[HEADER..])
    });
    let (result, rest) = outcome.value.map_err(|e| {
        if outcome.resource_exhausted {
            Error::TooLarge
        } else {
            Error::Format(e.to_string())
        }
    })?;
    if !rest.is_empty() {
        return Err(Error::Format("unexpected trailing payload".into()));
    }
    Ok(result)
}
// Postcard structs encode their fields consecutively; the borrowed prefix probes only
// id/revision/variant and cannot allocate the untrusted large variant's payload.
pub(super) fn worker_limit(
    magic: &[u8; 4],
    bytes: &[u8],
    frame_length: usize,
    allowed_large_variants: &[u32],
) -> Result<usize, Error> {
    let (_, len) = envelope_info(magic, bytes, frame_length, MAX_WORKSPACE_MESSAGE_BYTES)?;
    let ((_, revision, variant), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
            .map_err(|e| Error::Format(e.to_string()))?;
    if revision.len() > 128 {
        return Err(Error::TooLarge);
    }
    if len <= MAX_IMPORT_BYTES {
        return Ok(MAX_IMPORT_BYTES);
    }
    if !allowed_large_variants.iter().any(|&tag| tag == variant) {
        return Err(Error::TooLarge);
    }
    Ok(MAX_WORKSPACE_MESSAGE_BYTES)
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Molecule, Topology};

    #[test]
    fn native_model_resource_rejections_preserve_too_large_and_malformed_is_format() {
        let mut molecule = Molecule::dna("source", "A", Topology::Linear).unwrap();
        molecule
            .metadata
            .insert("source-evidence".into(), "x".repeat(1024 * 1024 + 1));
        // Construct a current, checksummed hostile frame directly, bypassing the
        // public saver whose borrowed admission correctly refuses this model.
        let oversized = encode(b"CG2\0", &molecule).unwrap();
        assert!(matches!(
            decode::<Molecule>(b"CG2\0", &oversized),
            Err(Error::TooLarge)
        ));
        let malformed = encode(b"CG2\0", &0u8).unwrap();
        assert!(matches!(
            decode::<Molecule>(b"CG2\0", &malformed),
            Err(Error::Format(_))
        ));
        let valid = Molecule::dna("valid", "A", Topology::Linear).unwrap();
        assert_eq!(
            decode::<Molecule>(b"CG2\0", &encode(b"CG2\0", &valid).unwrap()).unwrap(),
            valid
        );
    }
}
