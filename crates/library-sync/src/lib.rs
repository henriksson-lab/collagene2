//! A separate reliable-mutation namespace; the ordinary v1 library API is intact.
//! No actor, lab authority, native checkpoint, filesystem path or SQL is on wire.
#[cfg(not(target_arch = "wasm32"))]
pub mod native;
pub use collagene_library_protocol as library;
use library::{
    Digest, Document, DocumentId, Folder, FolderId, LabContext, Name, NativeSelection, Payload,
    PayloadKind, RevisionId, RevisionSummary, Text, WorkspaceSharing,
};
use serde::{Deserialize, Serialize};
pub const SCHEMA_VERSION: u16 = 1;
pub const HEADER_BYTES: usize = 42;
pub const MAX_MESSAGE_BYTES: usize = library::MAX_MESSAGE_BYTES;
pub const MAX_RESPONSE_BYTES: usize = 256 << 10;
pub const MEDIA_TYPE: &str = "application/vnd.collagene.library-sync+postcard";
/// Separate host route; no reinterpretation of ordinary /api/library v1 bytes.
pub const ENDPOINT: &str = "/api/library/sync";
pub const INTENT_DOMAIN: &[u8] = b"collagene-library-mutation-v1\0";
const REQUEST_MAGIC: &[u8; 4] = b"CSQ\0";
const RESPONSE_MAGIC: &[u8; 4] = b"CSR\0";
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("invalid sync field: {0}")]
    Invalid(&'static str),
    #[error("sync payload exceeds bounds")]
    TooLarge,
    #[error("unsupported sync schema {0}")]
    Version(u16),
    #[error("malformed sync frame")]
    Malformed,
    #[error("sync checksum mismatch")]
    Integrity,
    #[error("sync reply does not match the pinned mutation")]
    Mismatch,
    #[error("library field validation: {0}")]
    Library(#[from] library::Error),
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "Text<128>", into = "Text<128>")]
pub struct MutationId(Text<128>);
impl MutationId {
    pub fn new(s: impl AsRef<str>) -> Result<Self, Error> {
        Self::try_from(Text::new(s)?)
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
}
impl TryFrom<Text<128>> for MutationId {
    type Error = Error;
    fn try_from(s: Text<128>) -> Result<Self, Error> {
        let Some(rest) = s.strip_prefix("mut-") else {
            return Err(Error::Invalid("mutation ID"));
        };
        if rest.len() != 32
            || !rest
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(Error::Invalid("mutation ID"));
        }
        Ok(Self(s))
    }
}
impl From<MutationId> for Text<128> {
    fn from(id: MutationId) -> Self {
        id.0
    }
}
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum MutationOperation{
    CreateFolder{parent:Option<FolderId>,name:Name},
    RenameFolder{id:FolderId,expected_generation:u64,name:Name},
    MoveFolder{id:FolderId,expected_generation:u64,parent:Option<FolderId>},
    DeleteEmptyFolder{id:FolderId,expected_generation:u64},
    ImportDocument{folder:Option<FolderId>,name:Name,payload:Payload},
    SaveRevision{id:DocumentId,base_revision:RevisionId,payload:Payload},
    RenameDocument{id:DocumentId,expected_generation:u64,name:Name},
    MoveDocument{id:DocumentId,expected_generation:u64,folder:Option<FolderId>},
}}
dto! {pub struct Mutation{pub id:MutationId,pub operation:MutationOperation}}
// Field/variant order is a shared fixed v1 contract with native mutations.rs.
// Payload contents are bound through their verified digest, never a huge JSON
// integer-array serialization. Mutation nonce/private checkpoints are excluded.
#[derive(Serialize)]
struct Descriptor<'a> {
    kind: PayloadKind,
    archive_hash: &'a Digest,
    archive_size: u64,
    selection: &'a Option<NativeSelection>,
    sharing: &'a Option<WorkspaceSharing>,
}
impl<'a> From<&'a Payload> for Descriptor<'a> {
    fn from(p: &'a Payload) -> Self {
        Self {
            kind: p.kind,
            archive_hash: &p.archive_hash,
            archive_size: p.bytes.len() as u64,
            selection: &p.selection,
            sharing: &p.sharing,
        }
    }
}
#[derive(Serialize)]
enum IntentOperation<'a> {
    CreateFolder {
        parent: &'a Option<FolderId>,
        name: &'a Name,
    },
    RenameFolder {
        id: &'a FolderId,
        expected_generation: u64,
        name: &'a Name,
    },
    MoveFolder {
        id: &'a FolderId,
        expected_generation: u64,
        parent: &'a Option<FolderId>,
    },
    DeleteEmptyFolder {
        id: &'a FolderId,
        expected_generation: u64,
    },
    ImportDocument {
        folder: &'a Option<FolderId>,
        name: &'a Name,
        payload: Descriptor<'a>,
    },
    SaveRevision {
        id: &'a DocumentId,
        base_revision: &'a RevisionId,
        payload: Descriptor<'a>,
    },
    RenameDocument {
        id: &'a DocumentId,
        expected_generation: u64,
        name: &'a Name,
    },
    MoveDocument {
        id: &'a DocumentId,
        expected_generation: u64,
        folder: &'a Option<FolderId>,
    },
}
impl MutationOperation {
    pub fn payload(&self) -> Option<&Payload> {
        match self {
            Self::ImportDocument { payload, .. } | Self::SaveRevision { payload, .. } => {
                Some(payload)
            }
            _ => None,
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::CreateFolder { name, .. } => library::validate_name(name)?,
            Self::RenameFolder {
                expected_generation,
                name,
                ..
            }
            | Self::RenameDocument {
                expected_generation,
                name,
                ..
            } => {
                if *expected_generation == 0 {
                    return Err(Error::Invalid("metadata generation"));
                }
                library::validate_name(name)?;
            }
            Self::MoveFolder {
                id,
                expected_generation,
                parent,
            } => {
                if *expected_generation == 0 || parent.as_ref() == Some(id) {
                    return Err(Error::Invalid("folder move/generation"));
                }
            }
            Self::DeleteEmptyFolder {
                expected_generation,
                ..
            }
            | Self::MoveDocument {
                expected_generation,
                ..
            } => {
                if *expected_generation == 0 {
                    return Err(Error::Invalid("metadata generation"));
                }
            }
            Self::ImportDocument { name, payload, .. } => {
                library::validate_name(name)?;
                payload.validate()?;
            }
            Self::SaveRevision { payload, .. } => payload.validate()?,
        }
        Ok(())
    }
    pub fn intent_digest(&self) -> Result<Digest, Error> {
        use MutationOperation as O;
        let intent = match self {
            O::CreateFolder { parent, name } => IntentOperation::CreateFolder { parent, name },
            O::RenameFolder {
                id,
                expected_generation,
                name,
            } => IntentOperation::RenameFolder {
                id,
                expected_generation: *expected_generation,
                name,
            },
            O::MoveFolder {
                id,
                expected_generation,
                parent,
            } => IntentOperation::MoveFolder {
                id,
                expected_generation: *expected_generation,
                parent,
            },
            O::DeleteEmptyFolder {
                id,
                expected_generation,
            } => IntentOperation::DeleteEmptyFolder {
                id,
                expected_generation: *expected_generation,
            },
            O::ImportDocument {
                folder,
                name,
                payload,
            } => IntentOperation::ImportDocument {
                folder,
                name,
                payload: payload.into(),
            },
            O::SaveRevision {
                id,
                base_revision,
                payload,
            } => IntentOperation::SaveRevision {
                id,
                base_revision,
                payload: payload.into(),
            },
            O::RenameDocument {
                id,
                expected_generation,
                name,
            } => IntentOperation::RenameDocument {
                id,
                expected_generation: *expected_generation,
                name,
            },
            O::MoveDocument {
                id,
                expected_generation,
                folder,
            } => IntentOperation::MoveDocument {
                id,
                expected_generation: *expected_generation,
                folder,
            },
        };
        // Descriptor strings are bounded by the reused DTO constructors. No DNA
        // cloning, payload hashing, native decoding or scientific replay occurs.
        let json = serde_json::to_vec(&intent).map_err(|_| Error::Malformed)?;
        if json.len() > 16 * 1024 {
            return Err(Error::TooLarge);
        }
        let mut hash = blake3::Hasher::new();
        hash.update(INTENT_DOMAIN);
        hash.update(&json);
        Digest::new(hash.finalize().to_hex().as_str()).map_err(Into::into)
    }
}
impl Mutation {
    pub fn validate(&self) -> Result<(), Error> {
        self.operation.validate()
    }
    pub fn intent_digest(&self) -> Result<Digest, Error> {
        self.operation.intent_digest()
    }
}
dto! {pub enum Applied{
    CreatedFolder(Folder),RenamedFolder(Folder),MovedFolder(Folder),DeletedEmptyFolder(FolderId),
    Imported(Document),Saved{document:Document,revision:RevisionSummary},RenamedDocument(Document),MovedDocument(Document),
}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RejectionCode {
    Invalid,
    NotFound,
    AlreadyExists,
    NotEmpty,
    TooLarge,
    Conflict,
}
dto! {pub enum MutationConflict{Folder{current:Folder},Document{current:Document}}}
dto! {pub struct Rejection{pub code:RejectionCode,pub message:Text<1024>,pub conflict:Option<MutationConflict>}}
dto! {pub enum MutationOutcome{Applied(Applied),Rejected(Rejection)}}
dto! {pub struct Receipt{pub id:MutationId,pub intent:Digest,pub catalog_generation:u64,pub outcome:MutationOutcome}}
impl Applied {
    fn validate(&self) -> Result<(), Error> {
        match self {
            Self::CreatedFolder(f) | Self::RenamedFolder(f) | Self::MovedFolder(f) => {
                f.validate()?
            }
            Self::Imported(d) | Self::RenamedDocument(d) | Self::MovedDocument(d) => {
                d.validate()?
            }
            Self::Saved { document, revision } => {
                document.validate()?;
                revision.validate()?;
                if document.head != *revision {
                    return Err(Error::Invalid("saved head"));
                }
            }
            Self::DeletedEmptyFolder(_) => {}
        }
        Ok(())
    }
}
impl Rejection {
    fn validate(&self) -> Result<(), Error> {
        if self.message.trim().is_empty()
            || self
                .message
                .chars()
                .any(|c| c.is_control() && !matches!(c, '\n' | '\t'))
            || (self.code == RejectionCode::Conflict) != self.conflict.is_some()
        {
            return Err(Error::Invalid("durable rejection"));
        }
        match &self.conflict {
            Some(MutationConflict::Folder { current }) => current.validate()?,
            Some(MutationConflict::Document { current }) => current.validate()?,
            None => {}
        }
        Ok(())
    }
}
impl Receipt {
    pub fn validate(&self) -> Result<(), Error> {
        match &self.outcome {
            MutationOutcome::Applied(a) => a.validate(),
            MutationOutcome::Rejected(r) => r.validate(),
        }
    }
    /// Exact historical receipt; a replay never substitutes a newer head or
    /// recomputes an earlier conflict against today's catalog.
    pub fn validate_for(&self, mutation: &Mutation) -> Result<(), Error> {
        self.validate()?;
        if self.id != mutation.id || self.intent != mutation.intent_digest()? {
            return Err(Error::Mismatch);
        }
        use Applied as A;
        use MutationOperation as O;
        let matches = match (&mutation.operation, &self.outcome) {
            (O::CreateFolder { parent, name }, MutationOutcome::Applied(A::CreatedFolder(f))) => {
                f.parent == *parent && f.name == *name
            }
            (
                O::RenameFolder {
                    id,
                    expected_generation,
                    name,
                },
                MutationOutcome::Applied(A::RenamedFolder(f)),
            ) => f.id == *id && f.name == *name && f.generation > *expected_generation,
            (
                O::MoveFolder {
                    id,
                    expected_generation,
                    parent,
                },
                MutationOutcome::Applied(A::MovedFolder(f)),
            ) => f.id == *id && f.parent == *parent && f.generation > *expected_generation,
            (
                O::DeleteEmptyFolder { id, .. },
                MutationOutcome::Applied(A::DeletedEmptyFolder(f)),
            ) => id == f,
            (
                O::ImportDocument {
                    folder,
                    name,
                    payload,
                },
                MutationOutcome::Applied(A::Imported(d)),
            ) => {
                d.folder == *folder
                    && d.name == *name
                    && d.head.parent.is_none()
                    && payload_matches(&d.head, payload)
            }
            (
                O::SaveRevision {
                    id,
                    base_revision,
                    payload,
                },
                MutationOutcome::Applied(A::Saved { document, revision }),
            ) => {
                document.id == *id
                    && revision.parent.as_ref() == Some(base_revision)
                    && payload_matches(revision, payload)
            }
            (
                O::RenameDocument {
                    id,
                    expected_generation,
                    name,
                },
                MutationOutcome::Applied(A::RenamedDocument(d)),
            ) => d.id == *id && d.name == *name && d.generation > *expected_generation,
            (
                O::MoveDocument {
                    id,
                    expected_generation,
                    folder,
                },
                MutationOutcome::Applied(A::MovedDocument(d)),
            ) => d.id == *id && d.folder == *folder && d.generation > *expected_generation,
            (_, MutationOutcome::Rejected(Rejection { conflict: None, .. })) => true,
            (
                O::RenameFolder {
                    id,
                    expected_generation,
                    ..
                }
                | O::MoveFolder {
                    id,
                    expected_generation,
                    ..
                }
                | O::DeleteEmptyFolder {
                    id,
                    expected_generation,
                },
                MutationOutcome::Rejected(Rejection {
                    conflict: Some(MutationConflict::Folder { current }),
                    ..
                }),
            ) => current.id == *id && current.generation != *expected_generation,
            (
                O::SaveRevision {
                    id, base_revision, ..
                },
                MutationOutcome::Rejected(Rejection {
                    conflict: Some(MutationConflict::Document { current }),
                    ..
                }),
            ) => current.id == *id && current.head.id != *base_revision,
            (
                O::RenameDocument {
                    id,
                    expected_generation,
                    ..
                }
                | O::MoveDocument {
                    id,
                    expected_generation,
                    ..
                },
                MutationOutcome::Rejected(Rejection {
                    conflict: Some(MutationConflict::Document { current }),
                    ..
                }),
            ) => current.id == *id && current.generation != *expected_generation,
            _ => false,
        };
        if !matches {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}
fn payload_matches(r: &RevisionSummary, p: &Payload) -> bool {
    r.archive_hash == p.archive_hash
        && r.archive_size == p.bytes.len() as u64
        && r.kind == p.kind
        && r.selection == p.selection
        && r.sharing == p.sharing
}
dto! {pub enum Request{Apply{mutation:Mutation},Lookup{mutation_id:MutationId,expected_intent:Digest}}}
dto! {pub struct RequestEnvelope{pub request_id:u64,pub expected_lab:library::LabId,pub request:Request}}
dto! {pub enum ResponseOutcome{Receipt(Receipt),Unknown{mutation_id:MutationId,intent:Digest},IntentMismatch{mutation_id:MutationId},Failure(library::Failure)}}
dto! {pub struct ResponseEnvelope{pub request_id:u64,pub context:Option<LabContext>,pub outcome:ResponseOutcome}}
impl RequestEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("correlation"));
        }
        if let Request::Apply { mutation } = &self.request {
            mutation.validate()?;
        }
        Ok(())
    }
}
impl ResponseEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("correlation"));
        }
        if let Some(context) = &self.context {
            library::validate_name(&context.name)?;
        } else if !matches!(
            &self.outcome,
            ResponseOutcome::Failure(library::Failure {
                code: library::ErrorCode::Unauthenticated | library::ErrorCode::Forbidden,
                ..
            })
        ) {
            return Err(Error::Invalid("host lab context"));
        }
        match &self.outcome {
            ResponseOutcome::Receipt(r) => r.validate()?,
            ResponseOutcome::Failure(f) => {
                if f.message.trim().is_empty()
                    || f.message
                        .chars()
                        .any(|c| c.is_control() && !matches!(c, '\n' | '\t'))
                    || f.conflict.is_some()
                    || f.code == library::ErrorCode::Conflict
                {
                    return Err(Error::Invalid("non-durable boundary failure"));
                }
            }
            _ => {}
        }
        Ok(())
    }
    pub fn validate_for(&self, request: &RequestEnvelope) -> Result<(), Error> {
        self.validate()?;
        if self.request_id != request.request_id {
            return Err(Error::Mismatch);
        }
        if self
            .context
            .as_ref()
            .is_some_and(|context| context.lab_id != request.expected_lab)
        {
            return Err(Error::Mismatch);
        }
        let matches = match (&request.request, &self.outcome) {
            (Request::Apply { mutation }, ResponseOutcome::Receipt(receipt)) => {
                receipt.validate_for(mutation)?;
                true
            }
            (
                Request::Lookup {
                    mutation_id,
                    expected_intent,
                },
                ResponseOutcome::Receipt(r),
            ) => &r.id == mutation_id && &r.intent == expected_intent,
            (
                Request::Lookup {
                    mutation_id,
                    expected_intent,
                },
                ResponseOutcome::Unknown {
                    mutation_id: actual,
                    intent,
                },
            ) => actual == mutation_id && intent == expected_intent,
            (Request::Apply { mutation }, ResponseOutcome::IntentMismatch { mutation_id }) => {
                &mutation.id == mutation_id
            }
            (
                Request::Lookup { mutation_id, .. },
                ResponseOutcome::IntentMismatch {
                    mutation_id: actual,
                },
            ) => actual == mutation_id,
            (_, ResponseOutcome::Failure(_)) => true,
            _ => false,
        };
        if !matches {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
    pub fn validate_for_in_context(
        &self,
        request: &RequestEnvelope,
        lab: &library::LabId,
    ) -> Result<(), Error> {
        self.validate_for(request)?;
        if self.context.as_ref().is_none_or(|c| &c.lab_id != lab) {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}

struct Buffer {
    bytes: Vec<u8>,
    limit: usize,
}
impl postcard::ser_flavors::Flavor for Buffer {
    type Output = Vec<u8>;
    fn try_push(&mut self, b: u8) -> Result<(), postcard::Error> {
        if self.bytes.len() >= self.limit {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes.push(b);
        Ok(())
    }
    fn try_extend(&mut self, b: &[u8]) -> Result<(), postcard::Error> {
        if b.len() > self.limit.saturating_sub(self.bytes.len()) {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes.extend_from_slice(b);
        Ok(())
    }
    fn finalize(self) -> Result<Self::Output, postcard::Error> {
        Ok(self.bytes)
    }
}
fn encode<T: Serialize>(magic: &[u8; 4], v: &T, limit: usize) -> Result<Vec<u8>, Error> {
    let payload = postcard::serialize_with_flavor(
        v,
        Buffer {
            bytes: Vec::new(),
            limit,
        },
    )
    .map_err(|e| {
        if matches!(e, postcard::Error::SerializeBufferFull) {
            Error::TooLarge
        } else {
            Error::Malformed
        }
    })?;
    let mut bytes = Vec::with_capacity(HEADER_BYTES + payload.len());
    bytes.extend_from_slice(magic);
    bytes.extend_from_slice(&SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    Ok(bytes)
}
fn checked_payload<'a>(magic: &[u8; 4], bytes: &'a [u8], limit: usize) -> Result<&'a [u8], Error> {
    if bytes.len() < HEADER_BYTES || &bytes[..4] != magic {
        return Err(Error::Malformed);
    }
    let version = u16::from_le_bytes([bytes[4], bytes[5]]);
    if version != SCHEMA_VERSION {
        return Err(Error::Version(version));
    }
    let length =
        u32::from_le_bytes(bytes[6..10].try_into().map_err(|_| Error::Malformed)?) as usize;
    if length > limit {
        return Err(Error::TooLarge);
    }
    if length.checked_add(HEADER_BYTES) != Some(bytes.len()) {
        return Err(Error::Malformed);
    }
    let raw = &bytes[HEADER_BYTES..];
    if blake3::hash(raw).as_bytes() != &bytes[10..42] {
        return Err(Error::Integrity);
    }
    Ok(raw)
}
fn decode<T: for<'de> Deserialize<'de>>(
    magic: &[u8; 4],
    bytes: &[u8],
    limit: usize,
) -> Result<T, Error> {
    let raw = checked_payload(magic, bytes, limit)?;
    let (value, tail) = postcard::take_from_bytes(raw).map_err(|_| Error::Malformed)?;
    if !tail.is_empty() {
        return Err(Error::Malformed);
    }
    Ok(value)
}
/// Borrow just the exchange ID and required database pin before allocating a
/// typed operation/archive. This pin is a precondition, never authorization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RequestPrefix<'a> {
    pub request_id: u64,
    pub expected_lab: &'a str,
}
pub fn request_prefix(bytes: &[u8]) -> Result<RequestPrefix<'_>, Error> {
    let raw = checked_payload(REQUEST_MAGIC, bytes, MAX_MESSAGE_BYTES)?;
    let ((request_id, expected_lab), _) =
        postcard::take_from_bytes::<(u64, &str)>(raw).map_err(|_| Error::Malformed)?;
    if request_id == 0
        || expected_lab.len() != 36
        || !expected_lab.starts_with("lab-")
        || !expected_lab[4..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(Error::Invalid("request lab pin"));
    }
    Ok(RequestPrefix {
        request_id,
        expected_lab,
    })
}
pub fn encode_request(r: &RequestEnvelope) -> Result<Vec<u8>, Error> {
    r.validate()?;
    encode(REQUEST_MAGIC, r, MAX_MESSAGE_BYTES)
}
pub fn decode_request(bytes: &[u8]) -> Result<RequestEnvelope, Error> {
    let value: RequestEnvelope = decode(REQUEST_MAGIC, bytes, MAX_MESSAGE_BYTES)?;
    value.validate()?;
    Ok(value)
}
pub fn encode_response(r: &ResponseEnvelope) -> Result<Vec<u8>, Error> {
    r.validate()?;
    encode(RESPONSE_MAGIC, r, MAX_RESPONSE_BYTES)
}
pub fn decode_response(bytes: &[u8]) -> Result<ResponseEnvelope, Error> {
    let value: ResponseEnvelope = decode(RESPONSE_MAGIC, bytes, MAX_RESPONSE_BYTES)?;
    value.validate()?;
    Ok(value)
}
pub fn decode_response_for(
    bytes: &[u8],
    request: &RequestEnvelope,
) -> Result<ResponseEnvelope, Error> {
    let r = decode_response(bytes)?;
    r.validate_for(request)?;
    Ok(r)
}
