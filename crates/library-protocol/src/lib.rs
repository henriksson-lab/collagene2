//! Shared-library wire boundary, independent of scientific models and storage.
//! Authentication and the lab namespace come exclusively from a trusted host.
//! A checksum proves byte integrity, never identity, permission or physical DNA.
mod bounded;
pub use bounded::{Blob, Rows, Text};
use serde::{Deserialize, Serialize};

pub const SCHEMA_VERSION: u16 = 1;
pub const HEADER_BYTES: usize = 42;
pub const MAX_DOCUMENT_BYTES: usize = 32 << 20;
pub const MAX_ARCHIVE_BYTES: usize = 64 << 20;
pub const MAX_MESSAGE_BYTES: usize = 65 << 20;
pub const MAX_PAGE_ROWS: usize = 200;
pub const MAX_OFFSET: u32 = 100_000;
pub const MEDIA_TYPE: &str = "application/vnd.collagene.library+postcard";
/// A target precondition, never a source of membership or namespace authority.
pub const EXPECTED_LAB_HEADER: &str = "x-collagene-expected-lab";
const REQUEST_MAGIC: &[u8; 4] = b"CLQ\0";
const RESPONSE_MAGIC: &[u8; 4] = b"CLR\0";

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("library message exceeds its size bound")]
    TooLarge,
    #[error("invalid library field: {0}")]
    Invalid(&'static str),
    #[error("malformed library envelope or payload")]
    Malformed,
    #[error("unsupported library schema version {0}")]
    Version(u16),
    #[error("library payload checksum mismatch")]
    Integrity,
    #[error("library response does not match its request")]
    ResponseMismatch,
}
macro_rules! id {
    ($name:ident,$prefix:literal) => {
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
        #[serde(try_from = "Text<128>", into = "Text<128>")]
        pub struct $name(Text<128>);
        impl $name {
            pub fn new(s: impl AsRef<str>) -> Result<Self, Error> {
                Self::try_from(Text::new(s)?)
            }
            pub fn as_str(&self) -> &str {
                self.0.as_str()
            }
        }
        impl TryFrom<Text<128>> for $name {
            type Error = Error;
            fn try_from(s: Text<128>) -> Result<Self, Error> {
                let Some(h) = s.strip_prefix($prefix) else {
                    return Err(Error::Invalid("opaque identifier prefix"));
                };
                if h.len() != 32 || !lower_hex(h) {
                    return Err(Error::Invalid("opaque identifier"));
                }
                Ok(Self(s))
            }
        }
        impl From<$name> for Text<128> {
            fn from(s: $name) -> Self {
                s.0
            }
        }
    };
}
id!(LabId, "lab-");
id!(FolderId, "fld-");
id!(DocumentId, "doc-");
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "Text<64>", into = "Text<64>")]
pub struct Digest(Text<64>);
impl Digest {
    pub fn new(s: impl AsRef<str>) -> Result<Self, Error> {
        Self::try_from(Text::new(s)?)
    }
    pub fn of(bytes: &[u8]) -> Self {
        Self(Text::new(blake3::hash(bytes).to_hex().to_string()).expect("fixed hash"))
    }
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}
impl TryFrom<Text<64>> for Digest {
    type Error = Error;
    fn try_from(s: Text<64>) -> Result<Self, Error> {
        if s.len() != 64 || !lower_hex(&s) {
            return Err(Error::Invalid("digest"));
        }
        Ok(Self(s))
    }
}
impl From<Digest> for Text<64> {
    fn from(s: Digest) -> Self {
        s.0
    }
}
/// Server-derived immutable history identity; distinct from archive/content hash.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct RevisionId(Digest);
impl RevisionId {
    pub fn new(s: impl AsRef<str>) -> Result<Self, Error> {
        Ok(Self(Digest::new(s)?))
    }
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}
pub type Name = Text<1024>;
pub type NativeId = Text<128>;
fn lower_hex(s: &str) -> bool {
    s.bytes()
        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn native_id(s: &str) -> bool {
    !s.is_empty() && !s.chars().any(|c| c.is_control() || c.is_whitespace())
}
pub fn validate_name(name: &str) -> Result<(), Error> {
    if name.trim().is_empty()
        || name.len() > 1024
        || name.chars().count() > 256
        || name.chars().any(char::is_control)
    {
        return Err(Error::Invalid("display name"));
    }
    Ok(())
}
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PayloadKind {
    NativeDocument,
    WorkspaceArchive,
}
impl PayloadKind {
    pub fn maximum_bytes(self) -> usize {
        match self {
            Self::NativeDocument => MAX_DOCUMENT_BYTES,
            Self::WorkspaceArchive => MAX_ARCHIVE_BYTES,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WorkspaceSharing {
    DocumentClosure,
    ExplicitWholeWorkspace,
}
dto! {pub struct NativeSelection {pub document_id:NativeId,pub revision:Digest,pub head_event:Option<NativeId>}}
dto! {pub struct Payload {
    pub kind:PayloadKind,pub bytes:Blob<MAX_ARCHIVE_BYTES>,pub archive_hash:Digest,
    pub selection:Option<NativeSelection>,pub sharing:Option<WorkspaceSharing>,
}}
fn validate_payload_descriptor(
    kind: PayloadKind,
    size: usize,
    selection: &Option<NativeSelection>,
    sharing: &Option<WorkspaceSharing>,
) -> Result<(), Error> {
    if size == 0 || size > kind.maximum_bytes() {
        return Err(Error::TooLarge);
    }
    match kind {
        PayloadKind::NativeDocument if selection.is_some() || sharing.is_some() => {
            return Err(Error::Invalid("native document sharing/selection"));
        }
        PayloadKind::WorkspaceArchive if selection.is_none() || sharing.is_none() => {
            return Err(Error::Invalid("workspace sharing/selection required"));
        }
        _ => {}
    }
    if let Some(s) = selection {
        if !native_id(&s.document_id) || s.head_event.as_ref().is_some_and(|id| !native_id(id)) {
            return Err(Error::Invalid("native selection identity"));
        }
    }
    Ok(())
}
impl Payload {
    /// Byte/shape integrity only. Host MUST separately decode native science/archive
    /// formats and derive/approve the selected sharing closure before persistence.
    pub fn validate(&self) -> Result<(), Error> {
        validate_payload_descriptor(self.kind, self.bytes.len(), &self.selection, &self.sharing)?;
        if self.archive_hash != Digest::of(&self.bytes) {
            return Err(Error::Integrity);
        }
        Ok(())
    }
}
dto! {pub struct Folder {pub id:FolderId,pub parent:Option<FolderId>,pub name:Name,pub generation:u64}}
dto! {pub struct RevisionSummary {
    pub id:RevisionId,pub parent:Option<RevisionId>,pub number:u64,pub archive_hash:Digest,pub archive_size:u64,
    pub kind:PayloadKind,pub selection:Option<NativeSelection>,pub sharing:Option<WorkspaceSharing>,pub created_at:u64,
}}
dto! {pub struct Document {pub id:DocumentId,pub folder:Option<FolderId>,pub name:Name,pub generation:u64,pub head:RevisionSummary}}
dto! {pub enum Entry {Folder(Folder),Document(Document)}}
impl Folder {
    pub fn validate(&self) -> Result<(), Error> {
        validate_name(&self.name)?;
        if self.generation == 0 || self.parent.as_ref() == Some(&self.id) {
            return Err(Error::Invalid("folder generation/parent"));
        }
        Ok(())
    }
}
impl RevisionSummary {
    pub fn validate(&self) -> Result<(), Error> {
        if self.number == 0
            || self.parent.as_ref() == Some(&self.id)
            || (self.number == 1) != self.parent.is_none()
        {
            return Err(Error::Invalid("revision ancestry"));
        }
        let size = usize::try_from(self.archive_size).map_err(|_| Error::TooLarge)?;
        validate_payload_descriptor(self.kind, size, &self.selection, &self.sharing)
    }
}
impl Document {
    pub fn validate(&self) -> Result<(), Error> {
        validate_name(&self.name)?;
        if self.generation == 0 {
            return Err(Error::Invalid("document generation"));
        }
        self.head.validate()
    }
}
impl Entry {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Folder(f) => f.validate(),
            Self::Document(d) => d.validate(),
        }
    }
    pub fn name(&self) -> &str {
        match self {
            Self::Folder(f) => &f.name,
            Self::Document(d) => &d.name,
        }
    }
    pub fn id(&self) -> &str {
        match self {
            Self::Folder(f) => f.id.as_str(),
            Self::Document(d) => d.id.as_str(),
        }
    }
}
dto! {pub struct PageRequest {pub offset:u32,pub limit:u16,pub expected_catalog_generation:Option<u64>}}
impl Default for PageRequest {
    fn default() -> Self {
        Self {
            offset: 0,
            limit: 100,
            expected_catalog_generation: None,
        }
    }
}
impl PageRequest {
    pub fn validate(&self) -> Result<(), Error> {
        if self.limit == 0
            || usize::from(self.limit) > MAX_PAGE_ROWS
            || self.offset > MAX_OFFSET
            || (self.offset > 0 && self.expected_catalog_generation.is_none())
        {
            return Err(Error::Invalid("pagination"));
        }
        Ok(())
    }
}
dto! {pub enum SearchScope {AllFolders,Folder{folder:Option<FolderId>,recursive:bool}}}
dto! {pub enum Request {
    Browse{parent:Option<FolderId>,page:PageRequest},
    CreateFolder{parent:Option<FolderId>,name:Name},
    RenameFolder{id:FolderId,expected_generation:u64,name:Name},
    MoveFolder{id:FolderId,expected_generation:u64,parent:Option<FolderId>},
    DeleteEmptyFolder{id:FolderId,expected_generation:u64},
    ImportDocument{folder:Option<FolderId>,name:Name,payload:Payload},
    OpenDocument{id:DocumentId,revision:Option<RevisionId>},
    SaveRevision{id:DocumentId,base_revision:RevisionId,payload:Payload},
    RenameDocument{id:DocumentId,expected_generation:u64,name:Name},
    MoveDocument{id:DocumentId,expected_generation:u64,folder:Option<FolderId>},
    Search{query:Text<256>,scope:SearchScope,page:PageRequest},
    ListRevisions{id:DocumentId,page:PageRequest},
}}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Browse { page, .. } | Self::ListRevisions { page, .. } => page.validate()?,
            Self::CreateFolder { name, .. } => validate_name(name)?,
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
                    return Err(Error::Invalid("metadata CAS generation"));
                }
                validate_name(name)?;
            }
            Self::MoveFolder {
                id,
                expected_generation,
                parent,
            } => {
                if *expected_generation == 0 || parent.as_ref() == Some(id) {
                    return Err(Error::Invalid("folder move/CAS"));
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
                    return Err(Error::Invalid("metadata CAS generation"));
                }
            }
            Self::ImportDocument { name, payload, .. } => {
                validate_name(name)?;
                payload.validate()?;
            }
            Self::SaveRevision { payload, .. } => payload.validate()?,
            Self::Search { query, page, .. } => {
                if query.trim().is_empty() || query.chars().any(char::is_control) {
                    return Err(Error::Invalid("search query"));
                }
                page.validate()?;
            }
            Self::OpenDocument { .. } => {}
        }
        Ok(())
    }
}
dto! {pub struct Page {pub catalog_generation:u64,pub entries:Rows<Entry,MAX_PAGE_ROWS>,pub next_offset:Option<u32>}}
dto! {pub struct RevisionPage {pub document:DocumentId,pub catalog_generation:u64,pub revisions:Rows<RevisionSummary,MAX_PAGE_ROWS>,pub next_offset:Option<u32>}}
dto! {pub enum Success {
    Browsed(Page),CreatedFolder(Folder),RenamedFolder(Folder),MovedFolder(Folder),DeletedEmptyFolder(FolderId),
    Imported(Document),Opened{document:Document,revision:RevisionSummary,payload:Payload},Saved{document:Document,revision:RevisionSummary},
    RenamedDocument(Document),MovedDocument(Document),Searched(Page),Revisions(RevisionPage),
}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ErrorCode {
    Unauthenticated,
    Forbidden,
    NotFound,
    AlreadyExists,
    NotEmpty,
    Invalid,
    TooLarge,
    Conflict,
    Unsupported,
    Unavailable,
}
dto! {pub enum Conflict {Catalog{actual_generation:u64},Folder{current:Folder},Document{current:Document}}}
dto! {pub struct Failure {pub code:ErrorCode,pub message:Text<1024>,pub conflict:Option<Conflict>}}
dto! {pub enum Outcome {Success(Success),Failure(Failure)}}
dto! {pub struct RequestEnvelope {pub request_id:u64,pub request:Request}}
dto! {pub struct LabContext {pub lab_id:LabId,pub name:Name}}
dto! {pub struct ResponseEnvelope {pub request_id:u64,pub context:Option<LabContext>,pub outcome:Outcome}}

fn page_fields(_generation: u64, next: Option<u32>) -> Result<(), Error> {
    if next.is_some_and(|n| n == 0 || n > MAX_OFFSET) {
        return Err(Error::Invalid("page generation/continuation"));
    }
    Ok(())
}
impl Page {
    pub fn validate(&self) -> Result<(), Error> {
        page_fields(self.catalog_generation, self.next_offset)?;
        for e in self.entries.iter() {
            e.validate()?;
        }
        let mut ids = std::collections::BTreeSet::new();
        if self.entries.iter().any(|entry| !ids.insert(entry.id())) {
            return Err(Error::Invalid("duplicate page identity"));
        }
        if self
            .entries
            .windows(2)
            .any(|w| (w[0].name(), w[0].id()) >= (w[1].name(), w[1].id()))
        {
            return Err(Error::Invalid("page order or duplicate entries"));
        }
        Ok(())
    }
}
impl Success {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Browsed(p) | Self::Searched(p) => p.validate(),
            Self::CreatedFolder(f) | Self::RenamedFolder(f) | Self::MovedFolder(f) => f.validate(),
            Self::Imported(d) | Self::RenamedDocument(d) | Self::MovedDocument(d) => d.validate(),
            Self::DeletedEmptyFolder(_) => Ok(()),
            Self::Opened {
                document,
                revision,
                payload,
            } => {
                document.validate()?;
                revision.validate()?;
                payload.validate()?;
                if revision.number > document.head.number
                    || (revision.number == document.head.number && revision.id != document.head.id)
                    || revision.archive_hash != payload.archive_hash
                    || revision.archive_size != payload.bytes.len() as u64
                    || revision.kind != payload.kind
                    || revision.selection != payload.selection
                    || revision.sharing != payload.sharing
                {
                    return Err(Error::Invalid("opened revision payload binding"));
                }
                Ok(())
            }
            Self::Saved { document, revision } => {
                document.validate()?;
                revision.validate()?;
                if document.head != *revision {
                    return Err(Error::Invalid("saved head binding"));
                }
                Ok(())
            }
            Self::Revisions(p) => {
                page_fields(p.catalog_generation, p.next_offset)?;
                for r in p.revisions.iter() {
                    r.validate()?;
                }
                if p.revisions.windows(2).any(|w| {
                    w[0].number.checked_sub(1) != Some(w[1].number)
                        || w[0].parent.as_ref() != Some(&w[1].id)
                }) {
                    return Err(Error::Invalid("revision page order"));
                }
                Ok(())
            }
        }
    }
}
impl RequestEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("request correlation"));
        }
        self.request.validate()
    }
}
impl ResponseEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("response correlation"));
        }
        if let Some(context) = &self.context {
            validate_name(&context.name)?;
        }
        if self.context.is_none()
            && !matches!(
                &self.outcome,
                Outcome::Failure(Failure {
                    code: ErrorCode::Unauthenticated | ErrorCode::Forbidden,
                    ..
                })
            )
        {
            return Err(Error::Invalid("host lab context required"));
        }
        match &self.outcome {
            Outcome::Success(s) => s.validate(),
            Outcome::Failure(f) => {
                if f.message.trim().is_empty()
                    || f.message
                        .chars()
                        .any(|c| c.is_control() && !matches!(c, '\n' | '\t'))
                    || (f.code == ErrorCode::Conflict) != f.conflict.is_some()
                {
                    return Err(Error::Invalid("failure message/conflict"));
                }
                match &f.conflict {
                    Some(Conflict::Folder { current }) => current.validate(),
                    Some(Conflict::Document { current }) => current.validate(),
                    _ => Ok(()),
                }
            }
        }
    }
    pub fn validate_for_in_context(
        &self,
        request: &RequestEnvelope,
        expected: &LabId,
    ) -> Result<(), Error> {
        self.validate_for(request)?;
        if self.context.as_ref().is_none_or(|c| &c.lab_id != expected) {
            return Err(Error::ResponseMismatch);
        }
        Ok(())
    }
    /// Check server answer correlation, operation, selected identities and CAS
    /// ancestry. This does not authenticate a server or validate native science.
    pub fn validate_for(&self, request: &RequestEnvelope) -> Result<(), Error> {
        request.validate()?;
        self.validate()?;
        if self.request_id != request.request_id {
            return Err(Error::ResponseMismatch);
        }
        let Outcome::Success(s) = &self.outcome else {
            if let Outcome::Failure(Failure {
                conflict: Some(conflict),
                ..
            }) = &self.outcome
            {
                let matches = match (conflict, &request.request) {
                    (
                        Conflict::Catalog { actual_generation },
                        Request::Browse { page, .. }
                        | Request::Search { page, .. }
                        | Request::ListRevisions { page, .. },
                    ) => page
                        .expected_catalog_generation
                        .is_some_and(|expected| expected != *actual_generation),
                    (
                        Conflict::Folder { current },
                        Request::RenameFolder {
                            id,
                            expected_generation,
                            ..
                        }
                        | Request::MoveFolder {
                            id,
                            expected_generation,
                            ..
                        }
                        | Request::DeleteEmptyFolder {
                            id,
                            expected_generation,
                        },
                    ) => current.id == *id && current.generation != *expected_generation,
                    (
                        Conflict::Document { current },
                        Request::SaveRevision {
                            id, base_revision, ..
                        },
                    ) => current.id == *id && current.head.id != *base_revision,
                    (
                        Conflict::Document { current },
                        Request::RenameDocument {
                            id,
                            expected_generation,
                            ..
                        }
                        | Request::MoveDocument {
                            id,
                            expected_generation,
                            ..
                        },
                    ) => current.id == *id && current.generation != *expected_generation,
                    _ => false,
                };
                if !matches {
                    return Err(Error::ResponseMismatch);
                }
            }
            return Ok(());
        };
        let ok = match (&request.request, s) {
            (Request::Browse { parent, page }, Success::Browsed(p)) => {
                page_matches(page, p.catalog_generation, p.entries.len(), p.next_offset)
                    && p.entries.iter().all(|e| match e {
                        Entry::Folder(f) => f.parent == *parent,
                        Entry::Document(d) => d.folder == *parent,
                    })
            }
            (Request::CreateFolder { parent, name }, Success::CreatedFolder(f)) => {
                f.parent == *parent && f.name == *name
            }
            (
                Request::RenameFolder {
                    id,
                    expected_generation,
                    name,
                },
                Success::RenamedFolder(f),
            ) => f.id == *id && f.name == *name && f.generation > *expected_generation,
            (
                Request::MoveFolder {
                    id,
                    expected_generation,
                    parent,
                },
                Success::MovedFolder(f),
            ) => f.id == *id && f.parent == *parent && f.generation > *expected_generation,
            (Request::DeleteEmptyFolder { id, .. }, Success::DeletedEmptyFolder(actual)) => {
                id == actual
            }
            (
                Request::ImportDocument {
                    folder,
                    name,
                    payload,
                },
                Success::Imported(d),
            ) => {
                d.folder == *folder
                    && d.name == *name
                    && d.head.parent.is_none()
                    && revision_payload(&d.head, payload)
            }
            (
                Request::OpenDocument {
                    id,
                    revision: selected,
                },
                Success::Opened {
                    document, revision, ..
                },
            ) => {
                document.id == *id
                    && selected
                        .as_ref()
                        .map_or(document.head.id == revision.id, |id| id == &revision.id)
            }
            (
                Request::SaveRevision {
                    id,
                    base_revision,
                    payload,
                },
                Success::Saved { document, revision },
            ) => {
                document.id == *id
                    && revision.parent.as_ref() == Some(base_revision)
                    && revision_payload(revision, payload)
            }
            (
                Request::RenameDocument {
                    id,
                    expected_generation,
                    name,
                },
                Success::RenamedDocument(d),
            ) => d.id == *id && d.name == *name && d.generation > *expected_generation,
            (
                Request::MoveDocument {
                    id,
                    expected_generation,
                    folder,
                },
                Success::MovedDocument(d),
            ) => d.id == *id && d.folder == *folder && d.generation > *expected_generation,
            (Request::Search { page, .. }, Success::Searched(p)) => {
                page_matches(page, p.catalog_generation, p.entries.len(), p.next_offset)
            }
            (Request::ListRevisions { id, page }, Success::Revisions(p)) => {
                p.document == *id
                    && page_matches(page, p.catalog_generation, p.revisions.len(), p.next_offset)
            }
            _ => false,
        };
        if !ok {
            return Err(Error::ResponseMismatch);
        }
        Ok(())
    }
}
fn revision_payload(r: &RevisionSummary, p: &Payload) -> bool {
    r.archive_hash == p.archive_hash
        && r.archive_size == p.bytes.len() as u64
        && r.kind == p.kind
        && r.selection == p.selection
        && r.sharing == p.sharing
}
fn page_matches(request: &PageRequest, generation: u64, count: usize, next: Option<u32>) -> bool {
    count <= usize::from(request.limit)
        && request
            .expected_catalog_generation
            .is_none_or(|g| g == generation)
        && next.is_none_or(|n| n == request.offset.saturating_add(count as u32) && count > 0)
}

struct BoundedPayload(Vec<u8>);
impl postcard::ser_flavors::Flavor for BoundedPayload {
    type Output = Vec<u8>;
    fn try_push(&mut self, b: u8) -> Result<(), postcard::Error> {
        if self.0.len() >= MAX_MESSAGE_BYTES {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.0.push(b);
        Ok(())
    }
    fn try_extend(&mut self, b: &[u8]) -> Result<(), postcard::Error> {
        if b.len() > MAX_MESSAGE_BYTES.saturating_sub(self.0.len()) {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.0.extend_from_slice(b);
        Ok(())
    }
    fn finalize(self) -> Result<Self::Output, postcard::Error> {
        Ok(self.0)
    }
}
fn encode<T: Serialize>(magic: &[u8; 4], value: &T) -> Result<Vec<u8>, Error> {
    let payload =
        postcard::serialize_with_flavor(value, BoundedPayload(Vec::new())).map_err(|e| {
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
fn decode<T: for<'de> Deserialize<'de>>(magic: &[u8; 4], bytes: &[u8]) -> Result<T, Error> {
    if bytes.len() < HEADER_BYTES || &bytes[..4] != magic {
        return Err(Error::Malformed);
    }
    let version = u16::from_le_bytes([bytes[4], bytes[5]]);
    if version != SCHEMA_VERSION {
        return Err(Error::Version(version));
    }
    let size = u32::from_le_bytes(bytes[6..10].try_into().map_err(|_| Error::Malformed)?) as usize;
    if size > MAX_MESSAGE_BYTES {
        return Err(Error::TooLarge);
    }
    if size.checked_add(HEADER_BYTES) != Some(bytes.len()) {
        return Err(Error::Malformed);
    }
    let raw = &bytes[HEADER_BYTES..];
    if blake3::hash(raw).as_bytes() != &bytes[10..42] {
        return Err(Error::Integrity);
    }
    let (value, tail) = postcard::take_from_bytes(raw).map_err(|_| Error::Malformed)?;
    if !tail.is_empty() {
        return Err(Error::Malformed);
    }
    Ok(value)
}
pub fn encode_request(request: &RequestEnvelope) -> Result<Vec<u8>, Error> {
    request.validate()?;
    encode(REQUEST_MAGIC, request)
}
pub fn decode_request(bytes: &[u8]) -> Result<RequestEnvelope, Error> {
    let value: RequestEnvelope = decode(REQUEST_MAGIC, bytes)?;
    value.validate()?;
    Ok(value)
}
pub fn encode_response(response: &ResponseEnvelope) -> Result<Vec<u8>, Error> {
    response.validate()?;
    encode(RESPONSE_MAGIC, response)
}
pub fn decode_response(bytes: &[u8]) -> Result<ResponseEnvelope, Error> {
    let value: ResponseEnvelope = decode(RESPONSE_MAGIC, bytes)?;
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
