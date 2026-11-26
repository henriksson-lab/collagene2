//! Raw checked IPC for one trusted desktop command: collagene_rpc(bytes).
//! The host supplies window access and the owned storage root out of band.
//! There are no paths, identities, credentials or network commands in this DTO.
pub use collagene_library_protocol as library;
use library::{Blob, Digest, LabContext, Text};
use serde::{Deserialize, Serialize};
pub const SCHEMA_VERSION: u16 = 1;
pub const HEADER_BYTES: usize = 42;
pub const MAX_EXPORT_BYTES: usize = 64 << 20;
pub const MAX_ARCHIVE_BYTES: usize = 64 << 20;
pub const MAX_REQUEST_BYTES: usize = 129 << 20;
pub const MAX_RESPONSE_BYTES: usize = 65 << 20;
pub const MAX_CONTROL_BYTES: usize = 1024;
pub const MAX_PREFERENCE_BYTES: usize = 1024;
const REQUEST_MAGIC: &[u8; 4] = b"DTQ\0";
const RESPONSE_MAGIC: &[u8; 4] = b"DTR\0";

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("invalid desktop field: {0}")]
    Invalid(&'static str),
    #[error("desktop IPC bound exceeded")]
    TooLarge,
    #[error("unsupported desktop IPC schema {0}")]
    Version(u16),
    #[error("malformed desktop IPC frame")]
    Malformed,
    #[error("desktop IPC checksum mismatch")]
    Integrity,
    #[error("desktop response changed operation/request pins")]
    Mismatch,
    #[error(transparent)]
    Library(#[from] library::Error),
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "Text<128>", into = "Text<128>")]
pub struct OperationId(Text<128>);
fn operation_id(value: &str) -> bool {
    value.len() == 35
        && value.starts_with("op-")
        && value[3..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl OperationId {
    pub fn new(value: impl AsRef<str>) -> Result<Self, Error> {
        Self::try_from(Text::new(value)?)
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
}
impl TryFrom<Text<128>> for OperationId {
    type Error = Error;
    fn try_from(value: Text<128>) -> Result<Self, Error> {
        if !operation_id(&value) {
            return Err(Error::Invalid("operation ID"));
        }
        Ok(Self(value))
    }
}
impl From<OperationId> for Text<128> {
    fn from(value: OperationId) -> Self {
        value.0
    }
}
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WorkingSlot {
    Current,
    Previous,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OperationState {
    Unknown,
    Queued,
    Running,
    CancelRequested,
    Committing,
    Completed,
    Cancelled,
    Failed,
}
dto! {pub enum Command{
    Bootstrap,
    Library(library::Request),
    WorkingMetadata,
    ReadWorking{slot:WorkingSlot,expected_version:Option<Digest>},
    SaveCurrent{expected_version:Option<Digest>,archive:Blob<MAX_ARCHIVE_BYTES>},
    ReplaceWorking{expected_version:Option<Digest>,current:Blob<MAX_ARCHIVE_BYTES>,previous:Blob<MAX_ARCHIVE_BYTES>},
    RecoverPrevious{expected_version:Digest},
    ReadPreferences,
    SavePreferences{expected_version:Option<Digest>,bytes:Blob<MAX_PREFERENCE_BYTES>},
    Cancel{operation_id:OperationId},
    OperationStatus{operation_id:OperationId},
    ExportFile{suggested_name:Text<256>,bytes:Blob<MAX_EXPORT_BYTES>},
}}
dto! {pub struct RequestEnvelope{pub request_id:u64,pub operation_id:OperationId,pub command:Command}}
dto! {pub struct SlotMetadata{pub archive_hash:Digest,pub bytes:u64}}
dto! {pub struct WorkingMetadata{pub version:Option<Digest>,pub current:Option<SlotMetadata>,pub previous:Option<SlotMetadata>}}
dto! {pub struct Preferences{pub version:Option<Digest>,pub bytes:Blob<MAX_PREFERENCE_BYTES>}}
dto! {pub struct Bootstrap{pub library_context:LabContext,pub working:WorkingMetadata,pub preferences:Preferences}}
dto! {pub enum Success{
    Bootstrap(Bootstrap),Library(library::Success),WorkingMetadata(WorkingMetadata),
    WorkingArchive{slot:WorkingSlot,version:Digest,archive:Blob<MAX_ARCHIVE_BYTES>},
    Preferences(Preferences),Cancelled{operation_id:OperationId,status:OperationState},Status{operation_id:OperationId,status:OperationState},
    Exported{name:Text<256>,hash:Digest,bytes:u64},
}}
dto! {pub enum Conflict{Working(WorkingMetadata),Preferences(Preferences)}}
dto! {pub enum Failure{Invalid{message:Text<1024>},Conflict{current:Conflict},Cancelled,Unavailable,Library(library::Failure)}}
dto! {pub enum Outcome{Success(Success),Failure(Failure)}}
dto! {pub struct ResponseEnvelope{pub request_id:u64,pub operation_id:OperationId,pub outcome:Outcome}}
/// A display basename, never a path or filesystem capability.
pub fn validate_export_name(name: &str) -> Result<(), Error> {
    if name.trim().is_empty()
        || name.len() > 256
        || name
            .chars()
            .any(|c| c.is_control() || matches!(c, '/' | '\\' | ':'))
    {
        return Err(Error::Invalid("export basename"));
    }
    Ok(())
}
/// Stable, portable name owned by this operation. Suggested Unicode characters
/// are conservatively mapped; the original request still binds exact intent.
pub fn export_file_name(operation: &OperationId, suggested: &str) -> Result<Text<256>, Error> {
    validate_export_name(suggested)?;
    let mut safe: String = suggested
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_') {
                c
            } else {
                '_'
            }
        })
        .collect();
    if safe.len() > 180 {
        // Keep exchange-format extensions even for a long construct name.
        // Sanitization above makes byte indexes ASCII character boundaries.
        let extension = safe
            .rfind('.')
            .filter(|&dot| {
                dot > 0
                    && safe.len() - dot <= 17
                    && safe[dot + 1..].bytes().all(|b| b.is_ascii_alphanumeric())
            })
            .map(|dot| safe[dot..].to_owned())
            .unwrap_or_default();
        safe.truncate(180 - extension.len());
        safe.push_str(&extension);
    }
    Text::new(format!(
        "{}-{}",
        operation.as_str(),
        safe.trim_end_matches('.')
    ))
    .map_err(Error::from)
}
fn archive_shape(bytes: &[u8]) -> Result<(), Error> {
    if bytes.is_empty() {
        return Err(Error::Invalid("empty native working archive"));
    }
    if bytes.len() > MAX_ARCHIVE_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn message(s: &str) -> Result<(), Error> {
    if s.trim().is_empty()
        || s.chars()
            .any(|c| c.is_control() && !matches!(c, '\n' | '\t'))
    {
        return Err(Error::Invalid("public error message"));
    }
    Ok(())
}
impl RequestEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("request ID"));
        }
        match &self.command {
            Command::Library(r) => r.validate()?,
            Command::ExportFile { suggested_name, .. } => validate_export_name(suggested_name)?,
            Command::SaveCurrent { archive, .. } => archive_shape(archive)?,
            Command::ReplaceWorking {
                current, previous, ..
            } => {
                archive_shape(current)?;
                archive_shape(previous)?;
            }
            _ => {}
        }
        Ok(())
    }
}
impl SlotMetadata {
    pub fn validate(&self) -> Result<(), Error> {
        if self.bytes == 0 || self.bytes > MAX_ARCHIVE_BYTES as u64 {
            return Err(Error::Invalid("working slot size"));
        }
        Ok(())
    }
}
impl WorkingMetadata {
    pub fn validate(&self) -> Result<(), Error> {
        if self.version.is_none() && (self.current.is_some() || self.previous.is_some()) {
            return Err(Error::Invalid("unversioned working slots"));
        }
        for slot in [&self.current, &self.previous].into_iter().flatten() {
            slot.validate()?;
        }
        Ok(())
    }
}
impl Preferences {
    pub fn validate(&self) -> Result<(), Error> {
        if self.version.is_none() && !self.bytes.is_empty() {
            return Err(Error::Invalid("unversioned preferences"));
        }
        Ok(())
    }
}
impl Bootstrap {
    pub fn validate(&self) -> Result<(), Error> {
        library::validate_name(&self.library_context.name)?;
        self.working.validate()?;
        self.preferences.validate()
    }
}
impl Success {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Exported { name, bytes, .. } => {
                validate_export_name(name)?;
                if *bytes > MAX_EXPORT_BYTES as u64 {
                    return Err(Error::TooLarge);
                }
            }
            Self::Bootstrap(b) => b.validate()?,
            Self::Library(s) => s.validate()?,
            Self::WorkingMetadata(m) => m.validate()?,
            Self::WorkingArchive { archive, .. } => archive_shape(archive)?,
            Self::Preferences(p) => p.validate()?,
            Self::Cancelled { status, .. }
                if !matches!(
                    status,
                    OperationState::Cancelled
                        | OperationState::CancelRequested
                        | OperationState::Committing
                        | OperationState::Completed
                        | OperationState::Failed
                ) =>
            {
                return Err(Error::Invalid("cancellation acknowledgment"));
            }
            _ => {}
        }
        Ok(())
    }
}
impl Failure {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Invalid { message: s } => message(s)?,
            Self::Conflict {
                current: Conflict::Working(m),
            } => m.validate()?,
            Self::Conflict {
                current: Conflict::Preferences(p),
            } => p.validate()?,
            Self::Library(f) => {
                message(&f.message)?;
                if (f.code == library::ErrorCode::Conflict) != f.conflict.is_some() {
                    return Err(Error::Invalid("library conflict"));
                }
                match &f.conflict {
                    Some(library::Conflict::Folder { current }) => current.validate()?,
                    Some(library::Conflict::Document { current }) => current.validate()?,
                    _ => {}
                }
            }
            _ => {}
        }
        Ok(())
    }
}
fn slot_matches(m: &Option<SlotMetadata>, bytes: &[u8]) -> bool {
    m.as_ref()
        .is_some_and(|m| m.bytes == bytes.len() as u64 && m.archive_hash == Digest::of(bytes))
}
impl ResponseEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("request ID"));
        }
        match &self.outcome {
            Outcome::Success(s) => s.validate(),
            Outcome::Failure(f) => f.validate(),
        }
    }
    pub fn validate_for(&self, request: &RequestEnvelope) -> Result<(), Error> {
        self.validate()?;
        if self.request_id != request.request_id || self.operation_id != request.operation_id {
            return Err(Error::Mismatch);
        }
        let matches = match (&request.command, &self.outcome) {
            (
                Command::ExportFile {
                    suggested_name,
                    bytes,
                },
                Outcome::Success(Success::Exported {
                    name,
                    hash,
                    bytes: size,
                }),
            ) => {
                *name == export_file_name(&request.operation_id, suggested_name)?
                    && *hash == Digest::of(bytes)
                    && *size == bytes.len() as u64
            }
            (Command::Bootstrap, Outcome::Success(Success::Bootstrap(_))) => true,
            (Command::Library(r), Outcome::Success(Success::Library(s))) => library_matches(r, s),
            (Command::WorkingMetadata, Outcome::Success(Success::WorkingMetadata(_))) => true,
            (
                Command::ReadWorking {
                    slot,
                    expected_version,
                },
                Outcome::Success(Success::WorkingArchive {
                    slot: actual,
                    version,
                    ..
                }),
            ) => slot == actual && expected_version.as_ref().is_none_or(|v| v == version),
            (
                Command::SaveCurrent { archive, .. },
                Outcome::Success(Success::WorkingMetadata(m)),
            ) => m.version.is_some() && slot_matches(&m.current, archive),
            (
                Command::ReplaceWorking {
                    current, previous, ..
                },
                Outcome::Success(Success::WorkingMetadata(m)),
            ) => {
                m.version.is_some()
                    && slot_matches(&m.current, current)
                    && slot_matches(&m.previous, previous)
            }
            (Command::RecoverPrevious { .. }, Outcome::Success(Success::WorkingMetadata(m))) => {
                m.version.is_some() && m.current.is_some()
            }
            (Command::ReadPreferences, Outcome::Success(Success::Preferences(_))) => true,
            (Command::SavePreferences { bytes, .. }, Outcome::Success(Success::Preferences(p))) => {
                p.version.is_some() && p.bytes == *bytes
            }
            (
                Command::Cancel { operation_id },
                Outcome::Success(Success::Cancelled {
                    operation_id: actual,
                    ..
                }),
            ) => operation_id == actual,
            (
                Command::OperationStatus { operation_id },
                Outcome::Success(Success::Status {
                    operation_id: actual,
                    ..
                }),
            ) => operation_id == actual,
            (Command::Library(r), Outcome::Failure(Failure::Library(f))) => {
                library_failure_matches(r, f)
            }
            (
                Command::ReadWorking {
                    expected_version, ..
                }
                | Command::SaveCurrent {
                    expected_version, ..
                }
                | Command::ReplaceWorking {
                    expected_version, ..
                },
                Outcome::Failure(Failure::Conflict {
                    current: Conflict::Working(m),
                }),
            ) => &m.version != expected_version,
            (
                Command::RecoverPrevious { expected_version },
                Outcome::Failure(Failure::Conflict {
                    current: Conflict::Working(m),
                }),
            ) => m.version.as_ref() != Some(expected_version),
            (
                Command::SavePreferences {
                    expected_version, ..
                },
                Outcome::Failure(Failure::Conflict {
                    current: Conflict::Preferences(p),
                }),
            ) => &p.version != expected_version,
            (
                _,
                Outcome::Failure(
                    Failure::Invalid { .. } | Failure::Cancelled | Failure::Unavailable,
                ),
            ) => true,
            _ => false,
        };
        if !matches {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}
fn payload_matches(r: &library::RevisionSummary, p: &library::Payload) -> bool {
    r.archive_hash == p.archive_hash
        && r.archive_size == p.bytes.len() as u64
        && r.kind == p.kind
        && r.selection == p.selection
        && r.sharing == p.sharing
}
fn page_matches(r: &library::PageRequest, generation: u64, n: usize, next: Option<u32>) -> bool {
    n <= usize::from(r.limit)
        && r.expected_catalog_generation
            .is_none_or(|g| g == generation)
        && next.is_none_or(|offset| offset == r.offset.saturating_add(n as u32) && n > 0)
}
fn library_matches(r: &library::Request, s: &library::Success) -> bool {
    use library::{Request as R, Success as S};
    match (r, s) {
        (R::Browse { parent, page }, S::Browsed(p)) => {
            page_matches(page, p.catalog_generation, p.entries.len(), p.next_offset)
                && p.entries.iter().all(|e| match e {
                    library::Entry::Folder(f) => f.parent == *parent,
                    library::Entry::Document(d) => d.folder == *parent,
                })
        }
        (R::CreateFolder { parent, name }, S::CreatedFolder(f)) => {
            f.parent == *parent && f.name == *name
        }
        (
            R::RenameFolder {
                id,
                expected_generation,
                name,
            },
            S::RenamedFolder(f),
        ) => f.id == *id && f.name == *name && f.generation > *expected_generation,
        (
            R::MoveFolder {
                id,
                expected_generation,
                parent,
            },
            S::MovedFolder(f),
        ) => f.id == *id && f.parent == *parent && f.generation > *expected_generation,
        (R::DeleteEmptyFolder { id, .. }, S::DeletedEmptyFolder(actual)) => id == actual,
        (
            R::ImportDocument {
                folder,
                name,
                payload,
            },
            S::Imported(d),
        ) => {
            d.folder == *folder
                && d.name == *name
                && d.head.parent.is_none()
                && payload_matches(&d.head, payload)
        }
        (
            R::OpenDocument {
                id,
                revision: selected,
            },
            S::Opened {
                document, revision, ..
            },
        ) => {
            document.id == *id
                && selected
                    .as_ref()
                    .map_or(document.head.id == revision.id, |id| id == &revision.id)
        }
        (
            R::SaveRevision {
                id,
                base_revision,
                payload,
            },
            S::Saved { document, revision },
        ) => {
            document.id == *id
                && revision.parent.as_ref() == Some(base_revision)
                && payload_matches(revision, payload)
        }
        (
            R::RenameDocument {
                id,
                expected_generation,
                name,
            },
            S::RenamedDocument(d),
        ) => d.id == *id && d.name == *name && d.generation > *expected_generation,
        (
            R::MoveDocument {
                id,
                expected_generation,
                folder,
            },
            S::MovedDocument(d),
        ) => d.id == *id && d.folder == *folder && d.generation > *expected_generation,
        (R::Search { page, .. }, S::Searched(p)) => {
            page_matches(page, p.catalog_generation, p.entries.len(), p.next_offset)
        }
        (R::ListRevisions { id, page }, S::Revisions(p)) => {
            p.document == *id
                && page_matches(page, p.catalog_generation, p.revisions.len(), p.next_offset)
        }
        _ => false,
    }
}
fn library_failure_matches(request: &library::Request, failure: &library::Failure) -> bool {
    use library::{Conflict as C, Request as R};
    match (&failure.conflict, request) {
        (None, _) => true,
        (
            Some(C::Catalog { actual_generation }),
            R::Browse { page, .. } | R::Search { page, .. } | R::ListRevisions { page, .. },
        ) => page
            .expected_catalog_generation
            .is_some_and(|expected| expected != *actual_generation),
        (
            Some(C::Folder { current }),
            R::RenameFolder {
                id,
                expected_generation,
                ..
            }
            | R::MoveFolder {
                id,
                expected_generation,
                ..
            }
            | R::DeleteEmptyFolder {
                id,
                expected_generation,
            },
        ) => current.id == *id && current.generation != *expected_generation,
        (
            Some(C::Document { current }),
            R::SaveRevision {
                id, base_revision, ..
            },
        ) => current.id == *id && current.head.id != *base_revision,
        (
            Some(C::Document { current }),
            R::RenameDocument {
                id,
                expected_generation,
                ..
            }
            | R::MoveDocument {
                id,
                expected_generation,
                ..
            },
        ) => current.id == *id && current.generation != *expected_generation,
        _ => false,
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
fn encode<T: Serialize>(magic: &[u8; 4], value: &T, limit: usize) -> Result<Vec<u8>, Error> {
    let raw = postcard::serialize_with_flavor(
        value,
        Buffer {
            bytes: vec![],
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
    let mut bytes = Vec::with_capacity(HEADER_BYTES + raw.len());
    bytes.extend_from_slice(magic);
    bytes.extend_from_slice(&SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&(raw.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&raw).as_bytes());
    bytes.extend_from_slice(&raw);
    Ok(bytes)
}
fn payload<'a>(magic: &[u8; 4], bytes: &'a [u8], limit: usize) -> Result<&'a [u8], Error> {
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
    Ok(&bytes[HEADER_BYTES..])
}
fn checked<'a>(magic: &[u8; 4], bytes: &'a [u8], limit: usize) -> Result<&'a [u8], Error> {
    let raw = payload(magic, bytes, limit)?;
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
    let raw = checked(magic, bytes, limit)?;
    let (value, tail) = postcard::take_from_bytes(raw).map_err(|_| Error::Malformed)?;
    if !tail.is_empty() {
        return Err(Error::Malformed);
    }
    Ok(value)
}
pub fn encode_request(r: &RequestEnvelope) -> Result<Vec<u8>, Error> {
    r.validate()?;
    encode(REQUEST_MAGIC, r, MAX_REQUEST_BYTES)
}
pub fn decode_request(bytes: &[u8]) -> Result<RequestEnvelope, Error> {
    let value: RequestEnvelope = decode(REQUEST_MAGIC, bytes, MAX_REQUEST_BYTES)?;
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
    let value = decode_response(bytes)?;
    value.validate_for(request)?;
    Ok(value)
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Control {
    Cancel { operation_id: OperationId },
    OperationStatus { operation_id: OperationId },
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ControlRequest {
    pub request_id: u64,
    pub operation_id: OperationId,
    pub control: Control,
}
/// Untrusted lane-selection hint only. It performs no checksum/content
/// validation and grants no authority; the chosen blocking dispatcher must call
/// decode_request before access. Borrowing the prefix avoids copying or hashing
/// a 64 MiB export on the UI/async executor.
pub fn export_request_hint(bytes: &[u8]) -> bool {
    if bytes.len() < HEADER_BYTES || bytes.len() > MAX_REQUEST_BYTES || &bytes[..4] != REQUEST_MAGIC
    {
        return false;
    }
    postcard::take_from_bytes::<(u64, &str, u32)>(&bytes[HEADER_BYTES..])
        .is_ok_and(|((id, op, tag), _)| id != 0 && operation_id(op) && tag == 11)
}
/// Fast path only for small controls. None defers checksum/full operation
/// validation to decode_request; it does not validate a normal archive job.
pub fn decode_control(bytes: &[u8]) -> Result<Option<ControlRequest>, Error> {
    let raw = payload(REQUEST_MAGIC, bytes, MAX_REQUEST_BYTES)?;
    let ((request_id, id, tag), _) =
        postcard::take_from_bytes::<(u64, &str, u32)>(raw).map_err(|_| Error::Malformed)?;
    if request_id == 0 || !operation_id(id) {
        return Err(Error::Invalid("control prefix"));
    }
    if tag < 9 || tag == 11 {
        return Ok(None);
    }
    if tag > 10 {
        return Err(Error::Malformed);
    }
    if bytes.len() > MAX_CONTROL_BYTES {
        return Err(Error::TooLarge);
    }
    let request = decode_request(bytes)?;
    let control = match request.command {
        Command::Cancel { operation_id } => Control::Cancel { operation_id },
        Command::OperationStatus { operation_id } => Control::OperationStatus { operation_id },
        _ => return Err(Error::Malformed),
    };
    Ok(Some(ControlRequest {
        request_id: request.request_id,
        operation_id: request.operation_id,
        control,
    }))
}
