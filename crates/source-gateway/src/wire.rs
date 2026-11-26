//! Shared native/WASM DTOs. No provider token, signed URL or raw HTTP response is
//! returned by catalog discovery. Request Debug deliberately omits credentials.
use collagene_core::formats::ImportReport;
use serde::{Deserialize, Serialize};
use std::fmt;
pub const VERSION: u16 = 1;
pub const MAX_REQUEST: usize = 16 * 1024;
pub const MAX_RESPONSE: usize = 64 * 1024 * 1024;
const HEADER: usize = 42;
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub id: u64,
    /// Fresh cryptographically random client-session nonce; not an account ID.
    /// Together with the server's access check it scopes cancellation/catalogs.
    pub client: [u8; 32],
    pub operation: Operation,
}
impl fmt::Debug for Request {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SourceGatewayRequest")
            .field("id", &self.id)
            .finish_non_exhaustive()
    }
}
#[derive(Serialize, Deserialize)]
pub enum Operation {
    Import {
        source: String,
        addgene_token: Option<String>,
    },
    Catalog {
        source: String,
        addgene_token: Option<String>,
    },
    Select {
        handle: [u8; 32],
        response_blake3: String,
        sequence_id: u64,
        addgene_token: Option<String>,
    },
    Cancel {
        request_id: u64,
    },
    DisposeCatalog {
        handle: [u8; 32],
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Response {
    pub id: u64,
    pub result: Result<Output, Failure>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Output {
    Imported(ImportReport),
    Catalog(Catalog),
    Canceled { found: bool },
    CatalogDisposed { found: bool },
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Catalog {
    pub handle: [u8; 32],
    pub response_blake3: String,
    pub plasmid_id: u64,
    pub name: String,
    pub choices: Vec<Choice>,
    pub expires_in_seconds: u64,
}
impl fmt::Debug for Catalog {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SourceCatalog")
            .field("plasmid_id", &self.plasmid_id)
            .field("name", &self.name)
            .field("choices", &self.choices)
            .field("response_blake3", &self.response_blake3)
            .field("expires_in_seconds", &self.expires_in_seconds)
            .finish_non_exhaustive()
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Choice {
    pub sequence_id: u64,
    pub description: String,
    pub origin: Origin,
    pub completeness: Completeness,
    pub length: Option<usize>,
    pub catalog_sequence_blake3: Option<String>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Origin {
    Depositor,
    Addgene,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Completeness {
    Full,
    Partial,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ErrorCode {
    InvalidRequest,
    Forbidden,
    Busy,
    Timeout,
    Canceled,
    ExpiredCatalog,
    Provider,
    TooLarge,
    Internal,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Failure {
    pub code: ErrorCode,
    pub message: String,
    pub retry_after_seconds: Option<u64>,
}
impl Failure {
    pub fn new(code: ErrorCode, message: &'static str) -> Self {
        Self {
            code,
            message: message.into(),
            retry_after_seconds: None,
        }
    }
}
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum WireError {
    #[error("Invalid bounded source-gateway envelope")]
    Invalid,
    #[error("Source-gateway payload exceeds its supported size")]
    TooLarge,
}
struct Bounded {
    bytes: Vec<u8>,
    limit: usize,
}
impl postcard::ser_flavors::Flavor for Bounded {
    type Output = Vec<u8>;
    fn try_push(&mut self, byte: u8) -> postcard::Result<()> {
        self.try_extend(&[byte])
    }
    fn try_extend(&mut self, bytes: &[u8]) -> postcard::Result<()> {
        if bytes.len() > self.limit.saturating_sub(self.bytes.len()) {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes
            .try_reserve(bytes.len())
            .map_err(|_| postcard::Error::SerializeBufferFull)?;
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }
    fn finalize(self) -> postcard::Result<Vec<u8>> {
        Ok(self.bytes)
    }
}
fn encode<T: Serialize>(magic: &[u8; 4], value: &T, max: usize) -> Result<Vec<u8>, WireError> {
    let body: Vec<u8> = postcard::serialize_with_flavor(
        value,
        Bounded {
            bytes: Vec::new(),
            limit: max - HEADER,
        },
    )
    .map_err(|_| WireError::TooLarge)?;
    let mut out = Vec::new();
    out.try_reserve_exact(HEADER + body.len())
        .map_err(|_| WireError::TooLarge)?;
    out.extend_from_slice(magic);
    out.extend(VERSION.to_le_bytes());
    out.extend((body.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&body).as_bytes());
    out.extend(body);
    Ok(out)
}
fn decode<'a, T: Deserialize<'a>>(
    magic: &[u8; 4],
    bytes: &'a [u8],
    max: usize,
) -> Result<T, WireError> {
    if bytes.len() > max {
        return Err(WireError::TooLarge);
    }
    if bytes.len() < HEADER || &bytes[..4] != magic || bytes[4..6] != VERSION.to_le_bytes() {
        return Err(WireError::Invalid);
    }
    let length =
        u32::from_le_bytes(bytes[6..10].try_into().map_err(|_| WireError::Invalid)?) as usize;
    if length != bytes.len() - HEADER
        || blake3::hash(&bytes[HEADER..]).as_bytes() != &bytes[10..HEADER]
    {
        return Err(WireError::Invalid);
    }
    let (value, rest) =
        postcard::take_from_bytes(&bytes[HEADER..]).map_err(|_| WireError::Invalid)?;
    if !rest.is_empty() {
        return Err(WireError::Invalid);
    }
    Ok(value)
}
pub fn encode_request(value: &Request) -> Result<Vec<u8>, WireError> {
    encode(b"CGSQ", value, MAX_REQUEST)
}
pub fn decode_request(bytes: &[u8]) -> Result<Request, WireError> {
    decode(b"CGSQ", bytes, MAX_REQUEST)
}
pub fn encode_response(value: &Response) -> Result<Vec<u8>, WireError> {
    encode(b"CGSR", value, MAX_RESPONSE)
}
pub fn decode_response(bytes: &[u8]) -> Result<Response, WireError> {
    decode(b"CGSR", bytes, MAX_RESPONSE)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn versioned_envelope_integrity_bounds_and_secret_free_debug() {
        let request = Request {
            id: 42,
            client: [7; 32],
            operation: Operation::Catalog {
                source: "addgene:12345".into(),
                addgene_token: Some("PRIVATE_SOURCE_TOKEN".into()),
            },
        };
        assert!(!format!("{request:?}").contains("PRIVATE_SOURCE_TOKEN"));
        let bytes = encode_request(&request).unwrap();
        let decoded = decode_request(&bytes).unwrap();
        assert_eq!(decoded.id, 42);
        for position in [0, 4, 6, 10, bytes.len() - 1] {
            let mut bad = bytes.clone();
            bad[position] ^= 1;
            assert!(decode_request(&bad).is_err());
        }
        let mut trailing = bytes.clone();
        trailing.push(0);
        assert!(decode_request(&trailing).is_err());
        let request = Request {
            id: 1,
            client: [7; 32],
            operation: Operation::Import {
                source: "X".repeat(MAX_REQUEST),
                addgene_token: None,
            },
        };
        assert_eq!(encode_request(&request).unwrap_err(), WireError::TooLarge);
    }
}
