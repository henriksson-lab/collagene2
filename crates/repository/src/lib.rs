//! Offline owned-file repository. Network access is injected only into explicit
//! clone/pull/resume calls. Credentials never enter this database's checkpoints
//! or biological archives. Native graph identity is distinct from remote history.
#![cfg(not(target_arch = "wasm32"))]
pub mod folder_publication;
pub mod http;
pub mod publication_review;
mod remote;
use collagene_core::{
    protocol,
    workspace::{self, Workspace},
};
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sharing::{
    self as sharing, Selection,
    import::{HeadBookmark, ImportOptions},
};
use collagene_library_sync as sync;
pub use remote::{Channel, RemoteTransport, SyncResult, TransportError};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, io::Write, path::Path};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Storage(#[from] store::Error),
    #[error(transparent)]
    Transaction(#[from] durable::MutationError),
    #[error("invalid native archive or exact selection")]
    Native,
    #[error("SBOL source evidence is outside the supported publication scope")]
    SbolSource,
    #[error("invalid repository state or field: {0}")]
    Invalid(&'static str),
    #[error("local mutation cancelled before commit")]
    Cancelled,
    #[error("repository state changed; reload before retrying")]
    Stale,
    #[error("local mutation rejected: {0}")]
    LocalRejected(String),
    #[error("remote operation failed: {0}")]
    Remote(String),
    #[error("this local item has an unacknowledged push to that remote")]
    Pending,
    #[error("there is no native or metadata change to publish")]
    NoChanges,
    #[error("the parent folder must be explicitly bound/pushed first")]
    UnboundFolder,
    #[error("the remote document is already cloned; pull it explicitly")]
    AlreadyBound,
    #[error(
        "the remote replaced the native root identity; clone separately to preserve local work"
    )]
    NativeRootChanged,
    #[error(transparent)]
    Merge(#[from] sharing::import::ImportError),
}
pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RemoteConfig {
    pub id: String,
    pub name: String,
    pub endpoint: String,
    pub lab: wire::LabId,
}
impl RemoteConfig {
    fn validate(&self) -> Result<()> {
        let url = url::Url::parse(&self.endpoint).map_err(|_| Error::Invalid("remote endpoint"))?;
        if self.id.len() != 36
            || !self.id.starts_with("rem-")
            || !self.id[4..]
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
            || self.name.trim().is_empty()
            || self.name.len() > 256
            || self.name.chars().any(char::is_control)
            || self.endpoint.len() > 2048
            || !matches!(url.scheme(), "http" | "https")
            || url.host_str().is_none()
            || !url.username().is_empty()
            || url.password().is_some()
            || url.query().is_some()
            || url.fragment().is_some()
        {
            return Err(Error::Invalid(
                "remote locator (credentials/query/fragment forbidden)",
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub enum NativePin {
    Document { revision: String },
    Workspace(Selection),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DocumentBinding {
    pub remote: wire::Document,
    pub common_native: NativePin,
    pub common_local_revision: store::RevisionId,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FolderBinding {
    pub remote: wire::Folder,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SnapshotBookmark {
    pub document: store::DocumentId,
    pub revision: store::RevisionId,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
enum Target {
    Document(String),
    Folder(String),
    DeleteFolder(String),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Pending {
    id: sync::MutationId,
    remote: RemoteConfig,
    target: Target,
    chunks: usize,
    frame_hash: wire::Digest,
    native: Option<NativePin>,
    local_revision: Option<store::RevisionId>,
    attempted: bool,
    rejection: Option<sync::Receipt>,
}
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct State {
    epoch: u64,
    remotes: BTreeMap<String, RemoteConfig>,
    documents: BTreeMap<String, BTreeMap<String, DocumentBinding>>,
    folders: BTreeMap<String, BTreeMap<String, FolderBinding>>,
    pending: BTreeMap<String, Pending>,
    bookmarks: Vec<HeadBookmark>,
    #[serde(default)]
    snapshot_bookmarks: Vec<SnapshotBookmark>,
}
/// Token-free native outbox projection. Frames, archive bytes and credentials
/// remain private implementation details and are never exposed by this view.
#[derive(Debug, Clone)]
pub enum OutboxTarget {
    Document(store::DocumentId),
    Folder(store::FolderId),
    DeletedFolder(store::FolderId),
}
#[derive(Debug, Clone)]
pub struct OutboxMetadata {
    pub id: String,
    pub remote_id: String,
    pub target: OutboxTarget,
    pub local_revision: Option<store::RevisionId>,
    pub attempted: bool,
    pub rejection: Option<sync::Receipt>,
}
pub struct Repository {
    db: store::Library,
    access: store::LabAccess,
    state: State,
    checkpoint: Option<store::Digest>,
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > 4 << 20)
        {
            return Err(std::io::Error::other("repository metadata limit"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn state_bytes(state: &State) -> Result<Vec<u8>> {
    if state.remotes.len() > 128
        || state.documents.len() > 100_000
        || state.folders.len() > 100_000
        || state.pending.len() > 1000
        || state.bookmarks.len() > 1000
        || state.snapshot_bookmarks.len() > 1000
    {
        return Err(Error::Invalid("repository metadata bounds"));
    }
    let mut out = Bounded(vec![]);
    serde_json::to_writer(&mut out, state)
        .map_err(|_| Error::Invalid("repository metadata bounds"))?;
    Ok(out.0)
}
fn state_write(state: &State, expected: Option<store::Digest>) -> Result<durable::CheckpointWrite> {
    Ok(durable::CheckpointWrite {
        key: durable::CheckpointKey::repository("state"),
        expected,
        value: Some(durable::CheckpointValue {
            version: 1,
            bytes: state_bytes(state)?,
        }),
    })
}
fn mutation_error(e: Error) -> durable::MutationError {
    durable::MutationError::Aborted(e.to_string())
}
fn convert<T: Serialize, U: serde::de::DeserializeOwned>(value: &T) -> Result<U> {
    let mut out = Bounded(vec![]);
    serde_json::to_writer(&mut out, value).map_err(|_| Error::Invalid("metadata conversion"))?;
    serde_json::from_slice(&out.0).map_err(|_| Error::Invalid("metadata conversion"))
}
fn selected(payload: &store::Payload) -> Result<NativePin> {
    payload.validate()?;
    match payload.kind {
        store::PayloadKind::NativeDocument => {
            let molecule = protocol::load_document(&payload.bytes).map_err(|_| Error::Native)?;
            Ok(NativePin::Document {
                revision: molecule.revision(),
            })
        }
        store::PayloadKind::WorkspaceArchive => {
            let frame = workspace::load_workspace(&payload.bytes).map_err(|_| Error::Native)?;
            let pin = payload.selection.as_ref().ok_or(Error::Native)?;
            let d = frame
                .documents
                .iter()
                .find(|d| d.id == pin.document_id)
                .ok_or(Error::Native)?;
            if d.head != pin.revision.as_str() || d.head_event != pin.head_event {
                return Err(Error::Native);
            }
            Ok(NativePin::Workspace(Selection {
                document_id: d.id.clone(),
                revision: d.head.clone(),
                head_event: d.head_event.clone(),
            }))
        }
    }
}
/// Private full workspaces stay local; sharing describes archive scope, never
/// permission to publish it. Queueing separately derives a precise closure.
pub fn private_workspace(frame: &Workspace, pin: &Selection) -> Result<store::Payload> {
    let bytes = workspace::save_workspace(frame).map_err(|_| Error::Native)?;
    let payload = store::Payload {
        kind: store::PayloadKind::WorkspaceArchive,
        archive_hash: store::Digest::parse(wire::Digest::of(&bytes).as_str())?,
        bytes,
        selection: Some(store::NativeSelection {
            document_id: pin.document_id.clone(),
            revision: store::Digest::parse(pin.revision.clone())?,
            head_event: pin.head_event.clone(),
        }),
        sharing: Some(store::WorkspaceSharing::ExplicitWholeWorkspace),
    };
    selected(&payload)?;
    Ok(payload)
}
fn wire_payload(payload: store::Payload) -> Result<wire::Payload> {
    Ok(wire::Payload {
        kind: match payload.kind {
            store::PayloadKind::NativeDocument => wire::PayloadKind::NativeDocument,
            store::PayloadKind::WorkspaceArchive => wire::PayloadKind::WorkspaceArchive,
        },
        bytes: wire::Blob::new(payload.bytes).map_err(|_| Error::Native)?,
        archive_hash: wire::Digest::new(payload.archive_hash.as_str())
            .map_err(|_| Error::Native)?,
        selection: payload.selection.as_ref().map(convert).transpose()?,
        sharing: payload.sharing.as_ref().map(convert).transpose()?,
    })
}
fn stored_payload(payload: wire::Payload) -> Result<store::Payload> {
    Ok(store::Payload {
        kind: match payload.kind {
            wire::PayloadKind::NativeDocument => store::PayloadKind::NativeDocument,
            wire::PayloadKind::WorkspaceArchive => store::PayloadKind::WorkspaceArchive,
        },
        bytes: payload.bytes.into_vec(),
        archive_hash: store::Digest::parse(payload.archive_hash.as_str())?,
        selection: payload.selection.as_ref().map(convert).transpose()?,
        sharing: payload.sharing.as_ref().map(convert).transpose()?,
    })
}
fn publication_error(error: sharing::Error) -> Error {
    match error {
        sharing::Error::SbolSource(_) => Error::SbolSource,
        _ => Error::Native,
    }
}
fn publish(payload: store::Payload) -> Result<(wire::Payload, NativePin)> {
    let pin = selected(&payload)?;
    match &pin {
        NativePin::Document { .. } => Ok((wire_payload(payload)?, pin)),
        NativePin::Workspace(selection) => {
            let frame = workspace::load_workspace(&payload.bytes).map_err(|_| Error::Native)?;
            let closure =
                sharing::document_closure(&frame, selection).map_err(publication_error)?;
            Ok((
                wire::Payload {
                    kind: wire::PayloadKind::WorkspaceArchive,
                    archive_hash: wire::Digest::of(&closure.archive),
                    bytes: wire::Blob::new(closure.archive).map_err(|_| Error::Native)?,
                    selection: Some(wire::NativeSelection {
                        document_id: wire::NativeId::new(&selection.document_id)
                            .map_err(|_| Error::Native)?,
                        revision: wire::Digest::new(&selection.revision)
                            .map_err(|_| Error::Native)?,
                        head_event: selection
                            .head_event
                            .as_ref()
                            .map(|id| wire::NativeId::new(id).map_err(|_| Error::Native))
                            .transpose()?,
                    }),
                    sharing: Some(wire::WorkspaceSharing::DocumentClosure),
                },
                pin,
            ))
        }
    }
}
impl Repository {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let mut db = store::Library::open(path)?;
        let access = db.local_repository_access("offline-plasmids", "Local plasmid repository")?;
        let mut result = Self {
            db,
            access,
            state: State::default(),
            checkpoint: None,
        };
        result.reload()?;
        if result.checkpoint.is_none() {
            result.commit(State::default(), vec![])?;
        }
        Ok(result)
    }
    pub fn reload(&mut self) -> Result<()> {
        match self
            .db
            .read_checkpoint(&self.access, &durable::CheckpointKey::repository("state"))?
        {
            Some(c) => {
                if c.value.version != 1 || c.value.bytes.len() > 4 << 20 {
                    return Err(Error::Invalid("state version/size"));
                }
                let state: State = serde_json::from_slice(&c.value.bytes)
                    .map_err(|_| Error::Invalid("state encoding"))?;
                state_bytes(&state)?;
                for r in state.remotes.values() {
                    r.validate()?;
                }
                self.state = state;
                self.checkpoint = Some(c.digest);
            }
            None => {
                self.state = State::default();
                self.checkpoint = None;
            }
        }
        Ok(())
    }
    fn commit(&mut self, next: State, mut extra: Vec<durable::CheckpointWrite>) -> Result<()> {
        extra.push(state_write(&next, self.checkpoint.clone())?);
        self.db.commit_checkpoints(&self.access, &extra)?;
        self.reload()
    }
    fn local_outcome_with_id<F, G>(
        &mut self,
        operation: durable::MutationOperation,
        supplied_id: Option<durable::MutationId>,
        change: F,
        gate: G,
    ) -> Result<durable::MutationOutcome>
    where
        F: FnOnce(&mut State, &durable::Applied) -> Result<()>,
        G: FnOnce() -> bool,
    {
        let reliable = supplied_id.is_some();
        let mutation = durable::Mutation {
            id: match supplied_id {
                Some(id) => id,
                None => self.db.new_mutation_id()?,
            },
            operation,
        };
        if reliable {
            let intent = mutation.intent_digest()?;
            if let Some(receipt) = self
                .db
                .lookup_reliable(&self.access, &mutation.id, &intent)?
            {
                self.reload()?;
                return Ok(receipt.outcome);
            }
        }
        use durable::MutationOperation as O;
        if let O::ImportDocument { payload, .. } | O::SaveRevision { payload, .. } =
            &mutation.operation
        {
            selected(payload)?;
        }
        let mut next = self.state.clone();
        let expected = self.checkpoint.clone();
        let mut cancelled = false;
        let result = self
            .db
            .apply_reliable_with(&self.access, &mutation, |receipt| {
                let writes = if let durable::MutationOutcome::Applied(applied) = &receipt.outcome {
                    change(&mut next, applied).map_err(mutation_error)?;
                    next.epoch = next
                        .epoch
                        .checked_add(1)
                        .ok_or(durable::MutationError::Invalid("repository epoch"))?;
                    vec![state_write(&next, expected).map_err(mutation_error)?]
                } else {
                    vec![]
                };
                if !gate() {
                    cancelled = true;
                    return Err(durable::MutationError::Aborted(
                        "local operation cancelled before commit".into(),
                    ));
                }
                Ok(writes)
            });
        if cancelled {
            return Err(Error::Cancelled);
        }
        let applied = result?;
        self.reload()?;
        Ok(applied.receipt.outcome)
    }
    fn local_outcome<F, G>(
        &mut self,
        operation: durable::MutationOperation,
        change: F,
        gate: G,
    ) -> Result<durable::MutationOutcome>
    where
        F: FnOnce(&mut State, &durable::Applied) -> Result<()>,
        G: FnOnce() -> bool,
    {
        self.local_outcome_with_id(operation, None, change, gate)
    }
    /// Stable native mutation nonce for an IPC retry. The local repository scope
    /// and full intent bind the durable receipt; a matching replay reruns neither
    /// the mutation nor the checkpoint/gate callback, even after process restart.
    pub fn mutate_local_reliable_with_gate<F: FnOnce() -> bool>(
        &mut self,
        mutation: durable::Mutation,
        gate: F,
    ) -> Result<durable::MutationOutcome> {
        self.local_outcome_with_id(mutation.operation, Some(mutation.id), |_, _| Ok(()), gate)
    }
    fn local<F>(
        &mut self,
        operation: durable::MutationOperation,
        change: F,
    ) -> Result<durable::Applied>
    where
        F: FnOnce(&mut State, &durable::Applied) -> Result<()>,
    {
        match self.local_outcome(operation, change, || true)? {
            durable::MutationOutcome::Applied(a) => Ok(a),
            durable::MutationOutcome::Rejected(r) => Err(Error::LocalRejected(r.message)),
        }
    }
    /// Structured rejection for trusted native adapters, preserving atomic
    /// repository checkpoint semantics and the existing mutate_local wrapper.
    pub fn mutate_local_outcome(
        &mut self,
        operation: durable::MutationOperation,
    ) -> Result<durable::MutationOutcome> {
        self.mutate_local_outcome_with_gate(operation, || true)
    }
    /// Gate runs inside the transaction immediately before its prepared
    /// checkpoint/receipt commits. False rolls back the complete local effect.
    pub fn mutate_local_outcome_with_gate<F: FnOnce() -> bool>(
        &mut self,
        operation: durable::MutationOperation,
        gate: F,
    ) -> Result<durable::MutationOutcome> {
        self.local_outcome(operation, |_, _| Ok(()), gate)
    }
    /// Actual offline folder/document CRUD, including rename/move/empty delete.
    /// Every successful local mutation commits its repository checkpoint atomically.
    pub fn mutate_local(
        &mut self,
        operation: durable::MutationOperation,
    ) -> Result<durable::Applied> {
        self.local(operation, |_, _| Ok(()))
    }
    pub fn browse(
        &mut self,
        parent: Option<&store::FolderId>,
        offset: u32,
        limit: u16,
        generation: Option<u64>,
    ) -> Result<store::Page> {
        Ok(self
            .db
            .browse(&self.access, parent, offset, limit, generation)?)
    }
    pub fn document(&self, id: &store::DocumentId) -> Result<store::Document> {
        Ok(self.db.document(&self.access, id)?)
    }
    pub fn folder(&self, id: &store::FolderId) -> Result<store::Folder> {
        Ok(self.db.folder(&self.access, id)?)
    }
    pub fn open_document(
        &self,
        id: &store::DocumentId,
        revision: Option<&store::RevisionId>,
    ) -> Result<store::Payload> {
        let head = self.db.document(&self.access, id)?;
        let payload = self
            .db
            .payload(&self.access, id, revision.unwrap_or(&head.head.id))?;
        selected(&payload)?;
        Ok(payload)
    }
    pub fn add_remote(
        &mut self,
        name: &str,
        endpoint: &str,
        lab: wire::LabId,
    ) -> Result<RemoteConfig> {
        let id = self.db.new_mutation_id()?.as_str().replace("mut-", "rem-");
        let remote = RemoteConfig {
            id: id.clone(),
            name: name.into(),
            endpoint: endpoint.into(),
            lab,
        };
        remote.validate()?;
        let mut next = self.state.clone();
        next.remotes.insert(id, remote.clone());
        self.commit(next, vec![])?;
        Ok(remote)
    }
    /// Context is obtained from the owned native database, never from IPC data.
    pub fn local_context(&self) -> Result<wire::LabContext> {
        let lab = self.db.lab(&self.access)?;
        Ok(wire::LabContext {
            lab_id: wire::LabId::new(lab.id.as_str())
                .map_err(|_| Error::Invalid("local context"))?,
            name: wire::Text::new(lab.name).map_err(|_| Error::Invalid("local context"))?,
        })
    }
    /// Internal read/write snapshot pin for a trusted native service. This is
    /// synchronization metadata, not a credential or biological document hash.
    pub fn checkpoint_version(&self) -> Option<store::Digest> {
        self.checkpoint.clone()
    }
    pub fn outbox_metadata(&self) -> Result<Vec<OutboxMetadata>> {
        self.state
            .pending
            .values()
            .map(|pending| {
                let target = match &pending.target {
                    Target::Document(id) => {
                        OutboxTarget::Document(store::DocumentId::parse(id.clone())?)
                    }
                    Target::Folder(id) => OutboxTarget::Folder(store::FolderId::parse(id.clone())?),
                    Target::DeleteFolder(id) => {
                        OutboxTarget::DeletedFolder(store::FolderId::parse(id.clone())?)
                    }
                };
                Ok(OutboxMetadata {
                    id: pending.id.as_str().into(),
                    remote_id: pending.remote.id.clone(),
                    target,
                    local_revision: pending.local_revision.clone(),
                    attempted: pending.attempted,
                    rejection: pending.rejection.clone(),
                })
            })
            .collect()
    }
    pub fn remotes(&self) -> impl Iterator<Item = &RemoteConfig> {
        self.state.remotes.values()
    }
    pub fn document_binding(
        &self,
        id: &store::DocumentId,
        remote: &str,
    ) -> Option<&DocumentBinding> {
        self.state.documents.get(id.as_str())?.get(remote)
    }
    pub fn folder_binding(&self, id: &store::FolderId, remote: &str) -> Option<&FolderBinding> {
        self.state.folders.get(id.as_str())?.get(remote)
    }
    pub fn pending_ids(&self) -> Vec<String> {
        self.state.pending.keys().cloned().collect()
    }
    pub fn snapshot_bookmarks(&self) -> &[SnapshotBookmark] {
        &self.state.snapshot_bookmarks
    }
    pub fn bookmarks(&self) -> &[HeadBookmark] {
        &self.state.bookmarks
    }
}

/// Persistent, snapshot-pinned recursive folder checkout jobs.
pub mod folders;
