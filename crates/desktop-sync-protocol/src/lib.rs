//! Metadata-only desktop remote commands. The native host supplies repository,
//! window and credential authority. This protocol never transfers biological archives.
pub use collagene_desktop_protocol::OperationId;
pub use collagene_library_protocol as library;
use library::{Digest, DocumentId, FolderId, LabId, Rows, Text};
use serde::{Deserialize, Serialize};
use std::fmt;

pub const SCHEMA_VERSION: u16 = 4;
mod publication;
mod removal;
pub use publication::*;
pub use removal::*;
pub const HEADER_BYTES: usize = 42;
/// Limits include the complete frame header.
pub const MAX_REQUEST_BYTES: usize = 64 << 10;
pub const MAX_RESPONSE_BYTES: usize = 3 << 20;
pub const MAX_SECRET_BYTES: usize = 2048;
pub const MAX_PAGE_ROWS: usize = 200;
const REQUEST_MAGIC: &[u8; 4] = b"DSQ\0";
const RESPONSE_MAGIC: &[u8; 4] = b"DSR\0";

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("invalid desktop sync field: {0}")]
    Invalid(&'static str),
    #[error("desktop sync frame bound exceeded")]
    TooLarge,
    #[error("unsupported desktop sync schema {0}")]
    Version(u16),
    #[error("malformed desktop sync frame")]
    Malformed,
    #[error("desktop sync checksum mismatch")]
    Integrity,
    #[error("desktop sync reply binding mismatch")]
    Mismatch,
    #[error(transparent)]
    Library(#[from] library::Error),
}
macro_rules! dto {($($t:tt)*) => {#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)] #[serde(deny_unknown_fields)] $($t)*};}
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
                if s.len() != $prefix.len() + 32
                    || !s.starts_with($prefix)
                    || !s[$prefix.len()..]
                        .bytes()
                        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                {
                    return Err(Error::Invalid(stringify!($name)));
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
id!(RemoteId, "rem-");
id!(PushJobId, "mut-");
id!(TreeJobId, "tree-");
id!(ReviewId, "review-");

/// Intentionally has no Display/AsRef implementation and never appears in a reply.
/// A host must not log request frames or serialized commands containing this value.
#[derive(Clone, PartialEq, Eq, Serialize)]
#[serde(transparent)]
pub struct Secret(Text<MAX_SECRET_BYTES>);
impl Secret {
    pub fn new(s: impl AsRef<str>) -> Result<Self, Error> {
        let s = s.as_ref();
        if s.is_empty() || !s.bytes().all(|b| (0x21..=0x7e).contains(&b)) {
            return Err(Error::Invalid("credential"));
        }
        Ok(Self(Text::new(s)?))
    }
    /// Only the native credential provider should read this value.
    pub fn expose_secret(&self) -> &str {
        self.0.as_str()
    }
}
impl fmt::Debug for Secret {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("Secret([REDACTED])")
    }
}
impl<'de> Deserialize<'de> for Secret {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let text = Text::<MAX_SECRET_BYTES>::deserialize(d)?;
        Self::new(text.as_str()).map_err(serde::de::Error::custom)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CredentialRetention {
    Session,
    Remember,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PullMode {
    CleanFastForward,
    ExplicitCheckout,
}
dto! {pub struct Remote {pub id:RemoteId,pub name:Text<256>,pub endpoint:Text<2048>,pub expected_lab:LabId}}
dto! {pub struct RemoteContext {pub remote:Remote,pub lab:library::LabContext}}
dto! {pub struct LocalPageRequest {pub offset:u16,pub limit:u16}}
impl Default for LocalPageRequest {
    fn default() -> Self {
        Self {
            offset: 0,
            limit: 100,
        }
    }
}
dto! {pub enum TreeRoot {LabRoot,Folder{id:FolderId,parent:Option<FolderId>}}}
dto! {pub struct TreeLimits {pub max_nodes:u16,pub max_depth:u16,pub max_pages:u16,pub page_size:u16}}
impl Default for TreeLimits {
    fn default() -> Self {
        Self {
            max_nodes: 1000,
            max_depth: 64,
            max_pages: 2000,
            page_size: 100,
        }
    }
}
dto! {pub enum Command {
    Remotes,
    AddRemote{name:Text<256>,endpoint:Text<2048>,expected_lab:LabId},
    SetCredential{remote:RemoteId,secret:Secret,retention:CredentialRetention},
    DeleteCredential{remote:RemoteId},
    CredentialStatus{remote:RemoteId},
    ClearSession{remote:Option<RemoteId>},
    Browse{remote:RemoteId,parent:Option<FolderId>,page:library::PageRequest},
    Search{remote:RemoteId,query:Text<256>,scope:library::SearchScope,page:library::PageRequest},
    ListRevisions{remote:RemoteId,document:DocumentId,page:library::PageRequest},
    CloneDocument{remote:RemoteId,document:DocumentId,folder:Option<FolderId>},
    QueuePush{local_document:DocumentId,expected_local_revision:library::RevisionId,remote:RemoteId},
    ResumePush{job:PushJobId},
    Pending{page:LocalPageRequest},
    CancelUnattempted{job:PushJobId},
    DiscardRejected{job:PushJobId},
    Pull{local_document:DocumentId,expected_local_revision:library::RevisionId,remote:RemoteId,mode:PullMode},
    BeginTreeClone{remote:RemoteId,root:TreeRoot,destination:Option<FolderId>,limits:TreeLimits},
    ResumeTreeClone{job:TreeJobId,max_steps:u16},
    RetryTreeClone{job:TreeJobId},
    TreeJobs{page:LocalPageRequest},
    QueueTreeFolderMetadata{job:TreeJobId,local_folder:FolderId},
    QueueTreeDocumentMetadata{job:TreeJobId,local_document:DocumentId},
    ListBookmarks{page:LocalPageRequest,expected_checkpoint:Option<Digest>},
    PrepareMetadataPublication{target:PublicationTarget,remote:RemoteId},
    QueueMetadataPublication{review:MetadataReview},
    InspectPushChoices{local_document:DocumentId,expected_local_revision:library::RevisionId,remote:RemoteId},
    PrepareReviewedPush{local_document:DocumentId,expected_local_revision:library::RevisionId,remote:RemoteId,choices:Rows<PushChoice,128>},
    QueueReviewedPush{local_document:DocumentId,expected_local_revision:library::RevisionId,remote:RemoteId,review:ReviewId},
    DiscardPushReview{review:ReviewId},
    ListDeletedFolderPublications{remote:RemoteId,page:LocalPageRequest,expected_checkpoint:Option<Digest>},
    PrepareRemoteRemoval{remote:RemoteId},
    ConfirmRemoteRemoval{remote:RemoteId,review:ReviewId},
    DiscardRemoteRemovalReview{review:ReviewId},
    CheckRemoteRemovalOutcome{remote:RemoteId},
    Discussion{remote:RemoteId,request:library::discussion::Request},
}}
dto! {pub struct RequestEnvelope {pub request_id:u64,pub operation_id:OperationId,pub command:Command}}
dto! {pub enum LocalTarget {Document(DocumentId),Folder(FolderId),DeletedFolder(FolderId)}}
dto! {pub struct Rejection {pub code:library::ErrorCode,pub message:Text<1024>,pub conflict:Option<library::Conflict>}}
dto! {pub struct PendingPush {pub job:PushJobId,pub remote:RemoteId,pub target:LocalTarget,pub local_revision:Option<library::RevisionId>,pub attempted:bool,pub rejection:Option<Rejection>}}
dto! {pub struct PendingPage {pub items:Rows<PendingPush,MAX_PAGE_ROWS>,pub next_offset:Option<u16>}}
// Documentary local history references, never material or remote authority.
dto! {pub enum Bookmark {
    Snapshot{local_document:DocumentId,revision:library::RevisionId},
    NativeBranch{document:Text<128>,revision:Digest,event:Option<Text<128>>,undo_count:u32,redo_count:u32,incoming_revision:Digest,incoming_event:Option<Text<128>>},
}}
dto! {pub struct BookmarkPage {pub items:Rows<Bookmark,MAX_PAGE_ROWS>,pub next_offset:Option<u16>,pub checkpoint:Option<Digest>}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TreePhase {
    Discovering,
    Scanning,
    Importing,
    Complete,
    Invalidated,
    Blocked,
}
dto! {pub enum TreeIssue {
    RemoteUnavailable,CatalogChanged{actual:u64},RemoteRejected{code:library::ErrorCode},
    InvalidGraph,Limit,AlreadyBound,LocalParentChanged{folder:FolderId},LocalMutationRejected,LocalStateChanged,
}}
dto! {pub struct TreeStatus {
    pub job:TreeJobId,pub remote:RemoteId,pub root:TreeRoot,pub destination:Option<FolderId>,
    pub phase:TreePhase,pub catalog_generation:Option<u64>,pub scanned_nodes:u16,
    pub folders_created:u16,pub documents_created:u16,pub issue:Option<TreeIssue>,
}}
dto! {pub struct TreePage {pub items:Rows<TreeStatus,MAX_PAGE_ROWS>,pub next_offset:Option<u16>}}
dto! {pub enum PushOutcome {Applied,Rejected(Rejection)}}
dto! {pub struct PushReceipt {pub job:PushJobId,pub intent:Digest,pub catalog_generation:u64,pub outcome:PushOutcome}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CredentialState {
    Missing,
    Session,
    Remembered,
}
dto! {pub enum Success {
    Remotes(Rows<Remote,128>),AddedRemote(Remote),
    CredentialSet{remote:RemoteId,retention:CredentialRetention},CredentialDeleted{remote:RemoteId},
    CredentialStatus{remote:RemoteId,state:CredentialState},SessionsCleared{remote:Option<RemoteId>},
    Browsed(library::Page),Searched(library::Page),Revisions(library::RevisionPage),
    Cloned{source_document:DocumentId,source_revision:library::RevisionId,local:library::Document},
    Queued(PendingPush),Resumed(PushReceipt),Pending(PendingPage),
    CancelledUnattempted{job:PushJobId},DiscardedRejected{job:PushJobId},
    Pulled{local:library::Document},TreeStarted(TreeStatus),TreeUpdated(TreeStatus),TreeJobs(TreePage),
    Bookmarks(BookmarkPage),
    MetadataPrepared(MetadataReview),PushChoices(PushChoices),PushPrepared(PushPublicationReview),PushReviewDiscarded{review:ReviewId},DeletedFolderPublications(DeletedFolderPage),
    Discussion(library::discussion::ResponseEnvelope),
    RemoteRemovalPrepared{review:RemoteRemovalReview},RemoteRemoved(RemovalReceipt),RemoteRemovalReviewDiscarded{review:ReviewId},RemoteRemovalOutcome{remote:RemoteId,receipt:Option<RemovalReceipt>},
}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ErrorCode {
    Invalid,
    Unavailable,
    NotFound,
    AlreadyBound,
    Pending,
    NoChanges,
    UnboundFolder,
    NativeRootChanged,
    HistoryConflict,
    Cancelled,
    TooLarge,
    BookmarkChanged,
    ReviewChanged,
    CredentialForgetFailed,
    RemovalChangedAfterForget,
}
dto! {pub enum Failure {
    Error{code:ErrorCode,message:Text<1024>},
    LocalConflict{current:library::Document},
    Remote(library::Failure),
}}
dto! {pub enum Outcome {Success(Success),Failure(Failure)}}
dto! {pub struct ResponseEnvelope {pub request_id:u64,pub operation_id:OperationId,pub context:Option<RemoteContext>,pub outcome:Outcome}}

fn name(s: &str) -> Result<(), Error> {
    library::validate_name(s)?;
    Ok(())
}
fn endpoint(s: &str) -> Result<(), Error> {
    let origin = s
        .strip_prefix("https://")
        .or_else(|| s.strip_prefix("http://"))
        .ok_or(Error::Invalid("remote endpoint"))?;
    let authority = origin.strip_suffix('/').unwrap_or(origin);
    if authority.is_empty()
        || !authority
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b".-:[]".contains(&b))
    {
        return Err(Error::Invalid("remote endpoint origin"));
    }
    // The native adapter additionally parses/canonicalizes the origin and enforces
    // TLS / literal-loopback policy. A locator never grants network authority.
    Ok(())
}

impl Remote {
    pub fn validate(&self) -> Result<(), Error> {
        name(&self.name)?;
        endpoint(&self.endpoint)
    }
}
impl RemoteContext {
    pub fn validate(&self) -> Result<(), Error> {
        self.remote.validate()?;
        name(&self.lab.name)?;
        if self.lab.lab_id != self.remote.expected_lab {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}
impl LocalPageRequest {
    pub fn validate(&self) -> Result<(), Error> {
        if self.limit == 0 || usize::from(self.limit) > MAX_PAGE_ROWS {
            return Err(Error::Invalid("local page"));
        }
        Ok(())
    }
}
impl TreeLimits {
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_nodes == 0
            || self.max_nodes > 2000
            || self.max_depth == 0
            || self.max_depth > 128
            || self.max_pages == 0
            || self.max_pages > 10000
            || self.page_size == 0
            || usize::from(self.page_size) > MAX_PAGE_ROWS
        {
            return Err(Error::Invalid("tree limits"));
        }
        Ok(())
    }
}
impl Command {
    pub fn remote_id(&self) -> Option<&RemoteId> {
        match self {
            Self::SetCredential { remote, .. }
            | Self::DeleteCredential { remote }
            | Self::CredentialStatus { remote }
            | Self::ClearSession {
                remote: Some(remote),
            }
            | Self::Browse { remote, .. }
            | Self::Search { remote, .. }
            | Self::ListRevisions { remote, .. }
            | Self::Discussion { remote, .. }
            | Self::CloneDocument { remote, .. }
            | Self::QueuePush { remote, .. }
            | Self::PrepareRemoteRemoval { remote }
            | Self::ConfirmRemoteRemoval { remote, .. }
            | Self::CheckRemoteRemovalOutcome { remote }
            | Self::PrepareMetadataPublication { remote, .. }
            | Self::ListDeletedFolderPublications { remote, .. }
            | Self::InspectPushChoices { remote, .. }
            | Self::PrepareReviewedPush { remote, .. }
            | Self::QueueReviewedPush { remote, .. }
            | Self::Pull { remote, .. }
            | Self::BeginTreeClone { remote, .. } => Some(remote),
            Self::QueueMetadataPublication { review } => Some(&review.description.remote.id),
            _ => None,
        }
    }
    pub fn read_request(&self) -> Option<library::Request> {
        match self {
            Self::Browse { parent, page, .. } => Some(library::Request::Browse {
                parent: parent.clone(),
                page: page.clone(),
            }),
            Self::Search {
                query, scope, page, ..
            } => Some(library::Request::Search {
                query: query.clone(),
                scope: scope.clone(),
                page: page.clone(),
            }),
            Self::ListRevisions { document, page, .. } => Some(library::Request::ListRevisions {
                id: document.clone(),
                page: page.clone(),
            }),
            _ => None,
        }
    }
}
impl RequestEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("request ID"));
        }
        if let Some(read) = self.command.read_request() {
            read.validate()?;
        }
        match &self.command {
            Command::Discussion { request, .. } => request.validate()?,
            Command::ListDeletedFolderPublications {
                page,
                expected_checkpoint,
                ..
            } => {
                page.validate()?;
                if page.offset > 0 && expected_checkpoint.is_none() {
                    return Err(Error::Invalid("deleted folder continuation"));
                }
            }
            Command::QueueMetadataPublication { review } => review.description.validate()?,
            Command::PrepareReviewedPush { choices, .. } => validate_choices(choices)?,
            Command::AddRemote {
                name: n,
                endpoint: e,
                ..
            } => {
                name(n)?;
                endpoint(e)?;
            }
            Command::Pending { page } | Command::TreeJobs { page } => page.validate()?,
            Command::ListBookmarks {
                page,
                expected_checkpoint,
            } => {
                page.validate()?;
                if page.offset > 0 && expected_checkpoint.is_none() {
                    return Err(Error::Invalid("bookmark continuation checkpoint"));
                }
            }
            Command::BeginTreeClone { root, limits, .. } => {
                limits.validate()?;
                if let TreeRoot::Folder {
                    id,
                    parent: Some(parent),
                } = root
                {
                    if id == parent {
                        return Err(Error::Invalid("tree root parent"));
                    }
                }
            }
            Command::ResumeTreeClone { max_steps, .. } if *max_steps == 0 || *max_steps > 1000 => {
                return Err(Error::Invalid("tree steps"));
            }
            _ => {}
        }
        Ok(())
    }
}
fn public_message(s: &str) -> Result<(), Error> {
    if s.trim().is_empty()
        || s.chars()
            .any(|c| c.is_control() && !matches!(c, '\n' | '\t'))
    {
        return Err(Error::Invalid("public error message"));
    }
    Ok(())
}
fn rejection(r: &Rejection) -> Result<(), Error> {
    // This is a historical receipt summary, not a new library response. Its
    // authenticated response/context binding was checked when it entered the
    // native outbox; do not invent a lab context to validate its metadata.
    public_message(&r.message)?;
    if (r.code == library::ErrorCode::Conflict) != r.conflict.is_some() {
        return Err(Error::Invalid("rejection code/conflict"));
    }
    match &r.conflict {
        Some(library::Conflict::Folder { current }) => current.validate()?,
        Some(library::Conflict::Document { current }) => current.validate()?,
        _ => {}
    }
    Ok(())
}

impl PendingPush {
    pub fn validate(&self) -> Result<(), Error> {
        if let Some(r) = &self.rejection {
            if !self.attempted {
                return Err(Error::Invalid("unattempted rejection"));
            }
            rejection(r)?;
        }
        Ok(())
    }
}
impl TreeStatus {
    pub fn validate(&self) -> Result<(), Error> {
        if self.scanned_nodes > 2000
            || u32::from(self.folders_created) + u32::from(self.documents_created)
                > u32::from(self.scanned_nodes)
        {
            return Err(Error::Invalid("tree progress"));
        }
        if self.phase == TreePhase::Complete
            && (self.issue.is_some()
                || self.folders_created + self.documents_created != self.scanned_nodes)
        {
            return Err(Error::Invalid("complete tree progress"));
        }
        Ok(())
    }
}
impl ResponseEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("request ID"));
        }
        if let Some(c) = &self.context {
            c.validate()?;
        }
        match &self.outcome {
            Outcome::Failure(Failure::Error { message, .. }) => public_message(message)?,
            Outcome::Failure(Failure::LocalConflict { current }) => current.validate()?,
            Outcome::Failure(Failure::Remote(f)) => {
                library::ResponseEnvelope {
                    request_id: self.request_id,
                    context: self.context.as_ref().map(|c| c.lab.clone()),
                    outcome: library::Outcome::Failure(f.clone()),
                }
                .validate()?;
            }
            Outcome::Success(s) => match s {
                Success::Remotes(rs) => {
                    let mut ids = std::collections::BTreeSet::new();
                    for r in rs.iter() {
                        r.validate()?;
                        if !ids.insert(&r.id) {
                            return Err(Error::Invalid("duplicate remote"));
                        }
                    }
                }
                Success::AddedRemote(r) => r.validate()?,
                Success::MetadataPrepared(r) => r.description.validate()?,
                Success::RemoteRemovalPrepared { review } => review.description.validate()?,
                Success::RemoteRemoved(receipt) => receipt.validate()?,
                Success::RemoteRemovalOutcome { remote, receipt } => {
                    if let Some(receipt) = receipt {
                        receipt.validate()?;
                        if &receipt.remote.id != remote {
                            return Err(Error::Mismatch);
                        }
                    }
                }
                Success::DeletedFolderPublications(p) => {
                    for f in p.items.iter() {
                        f.remote_folder.validate()?;
                    }
                }
                Success::PushPrepared(r) => {
                    r.configuration.validate()?;
                    if r.remote != r.configuration.id {
                        return Err(Error::Mismatch);
                    }
                    if r.entries.iter().any(|entry| {
                        entry.retained_objects > 10000
                            || entry.omitted_objects > 10000
                            || entry.omitted_sequence_edges > 200000
                    }) {
                        return Err(Error::TooLarge);
                    }
                }
                Success::PushChoices(r) => {
                    let mut seen = std::collections::BTreeSet::new();
                    for c in r.requirements.iter() {
                        if c.component_iri.is_empty()
                            || !seen.insert((&c.original_hash, &c.component_iri))
                            || c.required_sequences
                                .iter()
                                .any(|s| !c.available_sequences.iter().any(|a| a == s))
                        {
                            return Err(Error::Invalid("push choices"));
                        }
                    }
                }
                Success::Cloned { local, .. } | Success::Pulled { local } => local.validate()?,
                Success::Queued(p) => p.validate()?,
                Success::Pending(p) => {
                    let mut ids = std::collections::BTreeSet::new();
                    for item in p.items.iter() {
                        item.validate()?;
                        if !ids.insert(&item.job) {
                            return Err(Error::Invalid("duplicate pending job"));
                        }
                    }
                }
                Success::Resumed(p) => {
                    if let PushOutcome::Rejected(r) = &p.outcome {
                        rejection(r)?
                    }
                }
                Success::TreeStarted(t) | Success::TreeUpdated(t) => t.validate()?,
                Success::TreeJobs(p) => {
                    let mut ids = std::collections::BTreeSet::new();
                    for item in p.items.iter() {
                        item.validate()?;
                        if !ids.insert(&item.job) {
                            return Err(Error::Invalid("duplicate tree job"));
                        }
                    }
                }
                Success::Browsed(p) | Success::Searched(p) => p.validate()?,
                Success::Revisions(p) => library::Success::Revisions(p.clone()).validate()?,
                Success::Discussion(reply) => reply.validate()?,
                Success::Bookmarks(p) => {
                    for item in p.items.iter() {
                        if let Bookmark::NativeBranch {
                            document,
                            event,
                            incoming_event,
                            ..
                        } = item
                        {
                            public_message(document)?;
                            for event in [event, incoming_event].into_iter().flatten() {
                                public_message(event)?;
                            }
                        }
                    }
                }
                _ => {}
            },
        }
        Ok(())
    }
    pub fn validate_in_remote(
        &self,
        request: &RequestEnvelope,
        expected: &Remote,
    ) -> Result<(), Error> {
        self.validate_for(request)?;
        expected.validate()?;
        if self.context.as_ref().is_none_or(|c| &c.remote != expected) {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
    pub fn validate_for(&self, r: &RequestEnvelope) -> Result<(), Error> {
        r.validate()?;
        self.validate()?;
        if self.request_id != r.request_id || self.operation_id != r.operation_id {
            return Err(Error::Mismatch);
        }
        if let (Some(expected), Some(c)) = (r.command.remote_id(), &self.context) {
            if &c.remote.id != expected {
                return Err(Error::Mismatch);
            }
        }
        if let Command::Discussion { remote, request } = &r.command {
            return match &self.outcome {
                Outcome::Success(Success::Discussion(reply)) => {
                    let context = self.context.as_ref().ok_or(Error::Mismatch)?;
                    if &context.remote.id != remote || context.lab != reply.context {
                        return Err(Error::Mismatch);
                    }
                    reply.validate_for(
                        &library::discussion::RequestEnvelope {
                            request_id: r.request_id,
                            request: request.clone(),
                        },
                        &context.remote.expected_lab,
                    )?;
                    Ok(())
                }
                Outcome::Failure(Failure::Error { .. }) => Ok(()),
                _ => Err(Error::Mismatch),
            };
        }
        if let Some(read) = r.command.read_request() {
            let outcome = match &self.outcome {
                Outcome::Success(Success::Browsed(p)) => {
                    library::Outcome::Success(library::Success::Browsed(p.clone()))
                }
                Outcome::Success(Success::Searched(p)) => {
                    library::Outcome::Success(library::Success::Searched(p.clone()))
                }
                Outcome::Success(Success::Revisions(p)) => {
                    library::Outcome::Success(library::Success::Revisions(p.clone()))
                }
                Outcome::Failure(Failure::Remote(f)) => library::Outcome::Failure(f.clone()),
                Outcome::Failure(Failure::Error { .. }) => return Ok(()),
                _ => return Err(Error::Mismatch),
            };
            let c = self.context.as_ref().ok_or(Error::Mismatch)?;
            let req = library::RequestEnvelope {
                request_id: r.request_id,
                request: read,
            };
            library::ResponseEnvelope {
                request_id: r.request_id,
                context: Some(c.lab.clone()),
                outcome,
            }
            .validate_for_in_context(&req, &c.remote.expected_lab)?;
            return Ok(());
        }
        let matches = match (&r.command, &self.outcome) {
            (
                Command::ListDeletedFolderPublications {
                    remote,
                    page,
                    expected_checkpoint,
                },
                Outcome::Success(Success::DeletedFolderPublications(v)),
            ) => {
                remote == &v.remote
                    && v.items.len() <= usize::from(page.limit)
                    && (page.offset == 0 || expected_checkpoint == &v.checkpoint)
                    && v.next_offset.is_none_or(|next| next > page.offset)
            }
            (
                Command::PrepareRemoteRemoval { remote },
                Outcome::Success(Success::RemoteRemovalPrepared { review }),
            ) => remote == &review.description.remote.id,
            (
                Command::ConfirmRemoteRemoval { remote, .. },
                Outcome::Success(Success::RemoteRemoved(receipt)),
            ) => remote == &receipt.remote.id,
            (
                Command::DiscardRemoteRemovalReview { review },
                Outcome::Success(Success::RemoteRemovalReviewDiscarded { review: actual }),
            ) => review == actual,
            (
                Command::CheckRemoteRemovalOutcome { remote },
                Outcome::Success(Success::RemoteRemovalOutcome {
                    remote: actual,
                    receipt,
                }),
            ) => remote == actual && receipt.as_ref().is_none_or(|r| &r.remote.id == remote),
            (
                Command::PrepareMetadataPublication { target, remote },
                Outcome::Success(Success::MetadataPrepared(v)),
            ) => target == &v.description.target && remote == &v.description.remote.id,
            (
                Command::QueueMetadataPublication { review },
                Outcome::Success(Success::Queued(p)),
            ) => {
                review.description.remote.id == p.remote
                    && match (&review.description.target, &p.target) {
                        (PublicationTarget::Folder(a), LocalTarget::Folder(b))
                        | (PublicationTarget::DeletedFolder(a), LocalTarget::DeletedFolder(b)) => {
                            a == b
                        }
                        (PublicationTarget::DocumentMetadata(a), LocalTarget::Document(b)) => {
                            a == b
                        }
                        _ => false,
                    }
            }
            (
                Command::InspectPushChoices {
                    local_document,
                    expected_local_revision,
                    remote,
                },
                Outcome::Success(Success::PushChoices(v)),
            ) => {
                local_document == &v.local_document
                    && expected_local_revision == &v.local_revision
                    && remote == &v.remote
            }
            (
                Command::PrepareReviewedPush {
                    local_document,
                    expected_local_revision,
                    remote,
                    ..
                },
                Outcome::Success(Success::PushPrepared(v)),
            ) => {
                local_document == &v.local_document
                    && expected_local_revision == &v.local_revision
                    && remote == &v.remote
            }
            (
                Command::QueueReviewedPush {
                    local_document,
                    expected_local_revision,
                    remote,
                    ..
                },
                Outcome::Success(Success::Queued(v)),
            ) => {
                remote == &v.remote
                    && v.target == LocalTarget::Document(local_document.clone())
                    && v.local_revision.as_ref() == Some(expected_local_revision)
            }
            (
                Command::DiscardPushReview { review },
                Outcome::Success(Success::PushReviewDiscarded { review: actual }),
            ) => review == actual,
            (Command::Remotes, Outcome::Success(Success::Remotes(_))) => true,
            (
                Command::AddRemote {
                    name,
                    endpoint,
                    expected_lab,
                },
                Outcome::Success(Success::AddedRemote(a)),
            ) => &a.name == name && &a.endpoint == endpoint && &a.expected_lab == expected_lab,
            (
                Command::SetCredential {
                    remote, retention, ..
                },
                Outcome::Success(Success::CredentialSet {
                    remote: actual,
                    retention: actual_retention,
                }),
            ) => remote == actual && retention == actual_retention,
            (
                Command::DeleteCredential { remote },
                Outcome::Success(Success::CredentialDeleted { remote: actual }),
            ) => remote == actual,
            (
                Command::CredentialStatus { remote },
                Outcome::Success(Success::CredentialStatus { remote: actual, .. }),
            ) => remote == actual,
            (
                Command::ClearSession { remote },
                Outcome::Success(Success::SessionsCleared { remote: actual }),
            ) => remote == actual,
            (
                Command::CloneDocument {
                    document, folder, ..
                },
                Outcome::Success(Success::Cloned {
                    source_document,
                    local,
                    ..
                }),
            ) => document == source_document && folder == &local.folder && self.context.is_some(),
            (
                Command::QueuePush {
                    local_document,
                    expected_local_revision,
                    remote,
                },
                Outcome::Success(Success::Queued(p)),
            ) => {
                p.target == LocalTarget::Document(local_document.clone())
                    && &p.remote == remote
                    && p.local_revision.as_ref() == Some(expected_local_revision)
                    && !p.attempted
                    && p.rejection.is_none()
            }
            (
                Command::QueueTreeFolderMetadata { local_folder, .. },
                Outcome::Success(Success::Queued(p)),
            ) => {
                p.target == LocalTarget::Folder(local_folder.clone())
                    && p.local_revision.is_none()
                    && !p.attempted
                    && p.rejection.is_none()
            }
            (
                Command::QueueTreeDocumentMetadata { local_document, .. },
                Outcome::Success(Success::Queued(p)),
            ) => {
                p.target == LocalTarget::Document(local_document.clone())
                    && p.local_revision.is_none()
                    && !p.attempted
                    && p.rejection.is_none()
            }
            (Command::ResumePush { job }, Outcome::Success(Success::Resumed(p))) => job == &p.job,
            (Command::Pending { page }, Outcome::Success(Success::Pending(p))) => {
                local_page(page, p.items.len(), p.next_offset)
            }
            (
                Command::CancelUnattempted { job },
                Outcome::Success(Success::CancelledUnattempted { job: actual }),
            )
            | (
                Command::DiscardRejected { job },
                Outcome::Success(Success::DiscardedRejected { job: actual }),
            ) => job == actual,
            (Command::Pull { local_document, .. }, Outcome::Success(Success::Pulled { local })) => {
                local_document == &local.id && self.context.is_some()
            }
            (
                Command::BeginTreeClone {
                    remote,
                    root,
                    destination,
                    ..
                },
                Outcome::Success(Success::TreeStarted(t)),
            ) => {
                remote == &t.remote
                    && root == &t.root
                    && destination == &t.destination
                    && t.phase
                        == if matches!(root, TreeRoot::LabRoot) {
                            TreePhase::Scanning
                        } else {
                            TreePhase::Discovering
                        }
            }
            (
                Command::ResumeTreeClone { job, .. } | Command::RetryTreeClone { job },
                Outcome::Success(Success::TreeUpdated(t)),
            ) => job == &t.job,
            (Command::TreeJobs { page }, Outcome::Success(Success::TreeJobs(p))) => {
                local_page(page, p.items.len(), p.next_offset)
            }
            (
                Command::ListBookmarks {
                    page,
                    expected_checkpoint,
                },
                Outcome::Success(Success::Bookmarks(p)),
            ) => {
                self.context.is_none()
                    && local_page(page, p.items.len(), p.next_offset)
                    && expected_checkpoint
                        .as_ref()
                        .is_none_or(|expected| p.checkpoint.as_ref() == Some(expected))
            }
            (
                Command::QueuePush {
                    local_document,
                    expected_local_revision,
                    ..
                }
                | Command::Pull {
                    local_document,
                    expected_local_revision,
                    ..
                },
                Outcome::Failure(Failure::LocalConflict { current }),
            ) => local_document == &current.id && expected_local_revision != &current.head.id,
            (_, Outcome::Failure(Failure::Error { .. })) => true,
            (_, Outcome::Failure(Failure::Remote(_))) => self.context.is_some(),
            _ => false,
        };
        if !matches {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}
fn local_page(request: &LocalPageRequest, n: usize, next: Option<u16>) -> bool {
    n <= usize::from(request.limit)
        && next.is_none_or(|v| n > 0 && usize::from(v) == usize::from(request.offset) + n)
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
            limit: limit - HEADER_BYTES,
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
    if bytes.len() > limit || length > limit - HEADER_BYTES {
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

/// Validate the exact current request/response correlation before encoding.
pub fn encode_response_for_frame(r: &ResponseEnvelope, request: &[u8]) -> Result<Vec<u8>, Error> {
    let q = decode_request(request)?;
    r.validate_for(&q)?;
    encode_response(r)
}
