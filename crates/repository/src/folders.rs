//! Explicit snapshot subtree checkout. Jobs are private native checkpoints,
//! not plasmid attachments or remote operations. Completed local heads are never
//! overwritten when a scan, immutable fetch or interrupted job is resumed.
use super::*;
use std::collections::{BTreeSet, VecDeque};
const MAX_JOBS: usize = 128;
const MAX_JOB_BYTES: usize = 8 << 20;
const MAX_NODES: usize = 2000;
const MAX_PAGES: usize = 10_000;

#[derive(Debug, thiserror::Error)]
pub enum FolderError {
    #[error(transparent)]
    Repository(#[from] Error),
    #[error("invalid tree job: {0}")]
    Invalid(&'static str),
    #[error("remote catalog changed to generation {actual}")]
    CatalogChanged { actual: u64 },
    #[error("remote request unavailable; tree job remains incomplete")]
    Unavailable,
    #[error("remote tree request rejected: {0:?}")]
    Rejected(wire::ErrorCode),
    #[error("remote tree reply is not bound to its request and lab")]
    Protocol,
    #[error("tree scan bounds exceeded")]
    Limit,
    #[error("a remote tree node is already bound locally; no existing local head was replaced")]
    AlreadyBound,
}
impl From<store::Error> for FolderError {
    fn from(e: store::Error) -> Self {
        Self::Repository(Error::Storage(e))
    }
}
impl From<durable::MutationError> for FolderError {
    fn from(e: durable::MutationError) -> Self {
        Self::Repository(Error::Transaction(e))
    }
}
pub type TreeResult<T> = std::result::Result<T, FolderError>;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum TreeRoot {
    LabRoot,
    Folder {
        id: wire::FolderId,
        parent: Option<wire::FolderId>,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TreeLimits {
    pub max_nodes: usize,
    pub max_depth: usize,
    pub max_pages: usize,
    pub page_size: u16,
}
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
impl TreeLimits {
    fn validate(&self) -> TreeResult<()> {
        if self.max_nodes == 0
            || self.max_nodes > MAX_NODES
            || self.max_depth == 0
            || self.max_depth > store::MAX_FOLDER_DEPTH
            || self.max_pages == 0
            || self.max_pages > MAX_PAGES
            || self.page_size == 0
            || usize::from(self.page_size) > wire::MAX_PAGE_ROWS
        {
            return Err(FolderError::Invalid("limits"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TreePhase {
    Discovering,
    Scanning,
    Importing,
    Complete,
    Invalidated,
    Blocked,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum TreeIssue {
    RemoteUnavailable,
    CatalogChanged { actual: u64 },
    RemoteRejected { code: wire::ErrorCode },
    InvalidGraph,
    Limit,
    AlreadyBound,
    LocalParentChanged { folder: store::FolderId },
    LocalMutationRejected,
    LocalStateChanged,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum TreeNode {
    Folder {
        folder: wire::Folder,
        depth: usize,
    },
    Document {
        document: wire::Document,
        depth: usize,
    },
}
impl TreeNode {
    fn id(&self) -> &str {
        match self {
            Self::Folder { folder, .. } => folder.id.as_str(),
            Self::Document { document, .. } => document.id.as_str(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Cursor {
    parent: Option<wire::FolderId>,
    offset: u32,
    depth: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TreeJob {
    pub version: u16,
    pub id: String,
    pub remote: RemoteConfig,
    pub root: TreeRoot,
    pub destination: Option<store::FolderId>,
    pub limits: TreeLimits,
    pub phase: TreePhase,
    pub catalog_generation: Option<u64>,
    pub pages: usize,
    pub nodes: Vec<TreeNode>,
    pub folders: BTreeMap<String, store::FolderId>,
    pub documents: BTreeMap<String, store::DocumentId>,
    pub local_folder_generations: BTreeMap<String, u64>,
    pub issue: Option<TreeIssue>,
    discovery: Option<Cursor>,
    queue: VecDeque<Cursor>,
    scan_complete: bool,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TreeProgress {
    pub id: String,
    pub phase: TreePhase,
    pub catalog_generation: Option<u64>,
    pub scanned_nodes: usize,
    pub folders_created: usize,
    pub documents_created: usize,
    pub issue: Option<TreeIssue>,
}
impl TreeJob {
    pub fn progress(&self) -> TreeProgress {
        TreeProgress {
            id: self.id.clone(),
            phase: self.phase,
            catalog_generation: self.catalog_generation,
            scanned_nodes: self.nodes.len(),
            folders_created: self.folders.len(),
            documents_created: self.documents.len(),
            issue: self.issue.clone(),
        }
    }
    fn add_node(&mut self, node: TreeNode) -> TreeResult<()> {
        if self.nodes.len() >= self.limits.max_nodes {
            return Err(FolderError::Limit);
        }
        if self.nodes.iter().any(|old| old.id() == node.id()) {
            return Err(FolderError::Protocol);
        }
        let depth = match &node {
            TreeNode::Folder { depth, .. } | TreeNode::Document { depth, .. } => *depth,
        };
        if depth > self.limits.max_depth {
            return Err(FolderError::Limit);
        }
        if let TreeNode::Folder { folder, .. } = &node {
            self.queue.push_back(Cursor {
                parent: Some(folder.id.clone()),
                offset: 0,
                depth,
            });
        }
        self.nodes.push(node);
        Ok(())
    }
    fn scan_step<T: RemoteTransport>(&mut self, transport: &mut T) -> TreeResult<()> {
        if let Some(cursor) = self.discovery.clone() {
            let page = browse(self, &cursor, transport)?;
            self.catalog_generation = Some(page.catalog_generation);
            self.pages += 1;
            let TreeRoot::Folder { id, .. } = &self.root else {
                return Err(FolderError::Invalid("root discovery"));
            };
            if let Some(found) = page.entries.iter().find_map(|e| match e {
                wire::Entry::Folder(f) if &f.id == id => Some(f.clone()),
                _ => None,
            }) {
                self.discovery = None;
                self.add_node(TreeNode::Folder {
                    folder: found,
                    depth: 1,
                })?;
                self.phase = TreePhase::Scanning;
            } else if let Some(offset) = page.next_offset {
                self.discovery = Some(Cursor { offset, ..cursor });
            } else {
                return Err(FolderError::Rejected(wire::ErrorCode::NotFound));
            }
            return Ok(());
        }
        let Some(cursor) = self.queue.front().cloned() else {
            self.scan_complete = true;
            self.phase = TreePhase::Importing;
            return Ok(());
        };
        let page = browse(self, &cursor, transport)?;
        self.catalog_generation = Some(page.catalog_generation);
        self.pages += 1;
        for entry in page.entries.into_vec() {
            match entry {
                wire::Entry::Folder(folder) => {
                    if folder.parent != cursor.parent {
                        return Err(FolderError::Protocol);
                    }
                    self.add_node(TreeNode::Folder {
                        folder,
                        depth: cursor.depth + 1,
                    })?;
                }
                wire::Entry::Document(document) => {
                    if document.folder != cursor.parent {
                        return Err(FolderError::Protocol);
                    }
                    self.add_node(TreeNode::Document {
                        document,
                        depth: cursor.depth + 1,
                    })?;
                }
            }
        }
        if let Some(offset) = page.next_offset {
            self.queue.front_mut().ok_or(FolderError::Protocol)?.offset = offset;
        } else {
            self.queue.pop_front();
        }
        if self.queue.is_empty() {
            self.scan_complete = true;
            self.phase = TreePhase::Importing;
        }
        Ok(())
    }
}
fn job_key(id: &str) -> TreeResult<durable::CheckpointKey> {
    if id.len() != 37
        || !id.starts_with("tree-")
        || !id[5..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(FolderError::Invalid("job ID"));
    }
    Ok(durable::CheckpointKey::repository(format!("job-{id}")))
}
struct JobBuffer(Vec<u8>);
impl Write for JobBuffer {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_JOB_BYTES.saturating_sub(self.0.len()) {
            return Err(std::io::Error::other("tree job limit"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encoded(job: &TreeJob) -> TreeResult<Vec<u8>> {
    job.limits.validate()?;
    job_key(&job.id)?;
    if job.version != 1
        || job.nodes.len() > job.limits.max_nodes
        || job.pages > job.limits.max_pages
        || job.queue.len() > job.limits.max_nodes
        || job.folders.len() + job.documents.len() > job.nodes.len()
        || job.local_folder_generations.len() > MAX_NODES + 1
    {
        return Err(FolderError::Invalid("persisted job bounds"));
    }
    let mut b = JobBuffer(Vec::new());
    serde_json::to_writer(&mut b, job).map_err(|_| FolderError::Limit)?;
    Ok(b.0)
}
fn job_write(
    job: &TreeJob,
    expected: Option<store::Digest>,
) -> TreeResult<durable::CheckpointWrite> {
    Ok(durable::CheckpointWrite {
        key: job_key(&job.id)?,
        expected,
        value: Some(durable::CheckpointValue {
            version: 1,
            bytes: encoded(job)?,
        }),
    })
}
fn exchange<T: RemoteTransport>(
    remote: &RemoteConfig,
    request: wire::Request,
    transport: &mut T,
) -> TreeResult<wire::Success> {
    let request = wire::RequestEnvelope {
        request_id: 1,
        request,
    };
    let bytes = wire::encode_request(&request).map_err(|_| FolderError::Protocol)?;
    let reply = transport
        .exchange(remote, Channel::Library, &bytes)
        .map_err(|_| FolderError::Unavailable)?;
    let reply = wire::decode_response_for(&reply, &request).map_err(|_| FolderError::Protocol)?;
    reply
        .validate_for_in_context(&request, &remote.lab)
        .map_err(|_| FolderError::Protocol)?;
    match reply.outcome {
        wire::Outcome::Success(s) => Ok(s),
        wire::Outcome::Failure(wire::Failure {
            conflict: Some(wire::Conflict::Catalog { actual_generation }),
            ..
        }) => Err(FolderError::CatalogChanged {
            actual: actual_generation,
        }),
        wire::Outcome::Failure(f) => Err(FolderError::Rejected(f.code)),
    }
}
fn browse<T: RemoteTransport>(
    job: &TreeJob,
    cursor: &Cursor,
    transport: &mut T,
) -> TreeResult<wire::Page> {
    if job.pages >= job.limits.max_pages {
        return Err(FolderError::Limit);
    }
    let result = exchange(
        &job.remote,
        wire::Request::Browse {
            parent: cursor.parent.clone(),
            page: wire::PageRequest {
                offset: cursor.offset,
                limit: job.limits.page_size,
                expected_catalog_generation: job.catalog_generation,
            },
        },
        transport,
    )?;
    let wire::Success::Browsed(page) = result else {
        return Err(FolderError::Protocol);
    };
    Ok(page)
}
fn issue(error: &FolderError) -> TreeIssue {
    match error {
        FolderError::Unavailable => TreeIssue::RemoteUnavailable,
        FolderError::CatalogChanged { actual } => TreeIssue::CatalogChanged { actual: *actual },
        FolderError::Rejected(code) => TreeIssue::RemoteRejected { code: *code },
        FolderError::Limit => TreeIssue::Limit,
        FolderError::AlreadyBound => TreeIssue::AlreadyBound,
        FolderError::Repository(Error::Stale)
        | FolderError::Repository(Error::Transaction(
            durable::MutationError::CheckpointConflict { .. },
        )) => TreeIssue::LocalStateChanged,
        _ => TreeIssue::InvalidGraph,
    }
}
fn exact_payload<T: RemoteTransport>(
    job: &TreeJob,
    pin: &wire::Document,
    transport: &mut T,
) -> TreeResult<store::Payload> {
    let reply = exchange(
        &job.remote,
        wire::Request::OpenDocument {
            id: pin.id.clone(),
            revision: Some(pin.head.id.clone()),
        },
        transport,
    )?;
    let wire::Success::Opened {
        document,
        revision,
        payload,
    } = reply
    else {
        return Err(FolderError::Protocol);
    };
    if document.id != pin.id || revision != pin.head {
        return Err(FolderError::Protocol);
    }
    let stored = stored_payload(payload)?;
    let native = selected(&stored)?;
    if let NativePin::Workspace(selection) = native {
        let frame = workspace::load_workspace(&stored.bytes)
            .map_err(|_| FolderError::Repository(Error::Native))?;
        let closure = sharing::inspect_document_closure(&frame, &selection)
            .map_err(|_| FolderError::Repository(Error::Native))?;
        // Never import unrelated native records under a selected-plasmid label.
        // Do not silently replace pinned bytes with a different current head.
        if frame
            .documents
            .iter()
            .map(|d| d.id.clone())
            .collect::<BTreeSet<_>>()
            != closure.documents
            || frame.revisions.keys().cloned().collect::<BTreeSet<_>>() != closure.revisions
            || frame
                .history
                .iter()
                .map(|e| e.id.clone())
                .collect::<BTreeSet<_>>()
                != closure.events
            || frame
                .candidates
                .iter()
                .map(|c| c.plan_revision.clone())
                .collect::<BTreeSet<_>>()
                != closure.candidates
            || frame.attachments.keys().cloned().collect::<BTreeSet<_>>() != closure.attachments
            || frame.documents.iter().any(|d| {
                closure
                    .document_heads
                    .get(&d.id)
                    .is_none_or(|h| h.revision != d.head || h.head_event != d.head_event)
            })
        {
            return Err(FolderError::Repository(Error::Native));
        }
    }
    Ok(stored)
}
impl Repository {
    fn load_tree_job(&mut self, id: &str) -> TreeResult<(TreeJob, store::Digest)> {
        let checkpoint = self
            .db
            .read_checkpoint(&self.access, &job_key(id)?)?
            .ok_or(FolderError::Invalid("unknown job"))?;
        if checkpoint.value.version != 1 || checkpoint.value.bytes.len() > MAX_JOB_BYTES {
            return Err(FolderError::Invalid("job checkpoint"));
        }
        let job: TreeJob = serde_json::from_slice(&checkpoint.value.bytes)
            .map_err(|_| FolderError::Invalid("job encoding"))?;
        if job.id != id {
            return Err(FolderError::Invalid("job key binding"));
        }
        encoded(&job)?;
        Ok((job, checkpoint.digest))
    }
    fn save_tree_job(&mut self, job: &TreeJob, expected: store::Digest) -> TreeResult<()> {
        self.db
            .commit_checkpoints(&self.access, &[job_write(job, Some(expected))?])?;
        Ok(())
    }
    pub fn tree_clone_job(&mut self, id: &str) -> TreeResult<TreeJob> {
        self.load_tree_job(id).map(|v| v.0)
    }
    pub fn tree_clone_jobs(&mut self) -> TreeResult<Vec<TreeJob>> {
        let value = self.db.read_checkpoint(
            &self.access,
            &durable::CheckpointKey::repository("tree-jobs"),
        )?;
        let ids: Vec<String> = match value {
            None => vec![],
            Some(c) => {
                if c.value.version != 1 || c.value.bytes.len() > 16384 {
                    return Err(FolderError::Invalid("job index"));
                }
                serde_json::from_slice(&c.value.bytes)
                    .map_err(|_| FolderError::Invalid("job index"))?
            }
        };
        if ids.len() > MAX_JOBS {
            return Err(FolderError::Limit);
        }
        ids.iter().map(|id| self.tree_clone_job(id)).collect()
    }
    /// Persist the selected scope before network access. Nothing is recursively
    /// uploaded, deleted or merged into an existing local native document.
    pub fn begin_tree_clone(
        &mut self,
        remote_id: &str,
        root: TreeRoot,
        destination: Option<store::FolderId>,
        limits: TreeLimits,
    ) -> TreeResult<String> {
        limits.validate()?;
        if matches!(&root,TreeRoot::Folder{id,parent:Some(parent)}if id==parent) {
            return Err(FolderError::Invalid(
                "selected folder cannot be its own parent",
            ));
        }
        self.reload()?;
        let remote = self
            .state
            .remotes
            .get(remote_id)
            .cloned()
            .ok_or(FolderError::Invalid("unknown remote"))?;
        let mut generations = BTreeMap::new();
        if let Some(id) = &destination {
            let folder = self.db.folder(&self.access, id)?;
            generations.insert(id.as_str().into(), folder.generation);
        }
        let (discovery, queue, phase) = match &root {
            TreeRoot::LabRoot => (
                None,
                VecDeque::from([Cursor {
                    parent: None,
                    offset: 0,
                    depth: 0,
                }]),
                TreePhase::Scanning,
            ),
            TreeRoot::Folder { parent, .. } => (
                Some(Cursor {
                    parent: parent.clone(),
                    offset: 0,
                    depth: 0,
                }),
                VecDeque::new(),
                TreePhase::Discovering,
            ),
        };
        let id = self.db.new_mutation_id()?.as_str().replace("mut-", "tree-");
        let job = TreeJob {
            version: 1,
            id: id.clone(),
            remote,
            root,
            destination,
            limits,
            phase,
            catalog_generation: None,
            pages: 0,
            nodes: vec![],
            folders: BTreeMap::new(),
            documents: BTreeMap::new(),
            local_folder_generations: generations,
            issue: None,
            discovery,
            queue,
            scan_complete: false,
        };
        let key = durable::CheckpointKey::repository("tree-jobs");
        let old = self.db.read_checkpoint(&self.access, &key)?;
        let mut ids: Vec<String> = match &old {
            None => vec![],
            Some(c) => {
                if c.value.version != 1 || c.value.bytes.len() > 16384 {
                    return Err(FolderError::Invalid("job index"));
                }
                serde_json::from_slice(&c.value.bytes)
                    .map_err(|_| FolderError::Invalid("job index"))?
            }
        };
        if ids.len() >= MAX_JOBS {
            return Err(FolderError::Limit);
        }
        ids.push(id.clone());
        self.db.commit_checkpoints(
            &self.access,
            &[
                // A retiring/reconfigured remote must invalidate a previously
                // validated tree creation, even when the job index was unchanged.
                state_write(&self.state, self.checkpoint.clone())?,
                job_write(&job, None)?,
                durable::CheckpointWrite {
                    key,
                    expected: old.map(|c| c.digest),
                    value: Some(durable::CheckpointValue {
                        version: 1,
                        bytes: serde_json::to_vec(&ids)
                            .map_err(|_| FolderError::Invalid("job index"))?,
                    }),
                },
            ],
        )?;
        Ok(id)
    }
    pub fn restart_tree_scan(&mut self, id: &str) -> TreeResult<TreeProgress> {
        let (mut job, digest) = self.load_tree_job(id)?;
        if !job.folders.is_empty() || !job.documents.is_empty() {
            return Err(FolderError::Invalid(
                "cannot discard a partially imported snapshot",
            ));
        }
        job.nodes.clear();
        job.catalog_generation = None;
        job.pages = 0;
        job.scan_complete = false;
        job.issue = None;
        match &job.root {
            TreeRoot::LabRoot => {
                job.phase = TreePhase::Scanning;
                job.discovery = None;
                job.queue = VecDeque::from([Cursor {
                    parent: None,
                    offset: 0,
                    depth: 0,
                }]);
            }
            TreeRoot::Folder { parent, .. } => {
                job.phase = TreePhase::Discovering;
                job.discovery = Some(Cursor {
                    parent: parent.clone(),
                    offset: 0,
                    depth: 0,
                });
                job.queue.clear();
            }
        }
        self.save_tree_job(&job, digest)?;
        Ok(job.progress())
    }
    /// Explicitly retry a blocked import after the caller resolves its conflict.
    /// No applied node is forgotten and no immutable remote pin is refreshed.
    pub fn retry_tree_clone(&mut self, id: &str) -> TreeResult<TreeProgress> {
        let (mut job, digest) = self.load_tree_job(id)?;
        if job.phase != TreePhase::Blocked || !job.scan_complete {
            return Err(FolderError::Invalid(
                "retry requires a blocked completed scan",
            ));
        }
        job.phase = TreePhase::Importing;
        job.issue = None;
        self.save_tree_job(&job, digest)?;
        Ok(job.progress())
    }
    /// Explicit acknowledgment of a locally renamed/moved parent, after reviewing
    /// its current generation. Descendant imports keep that private metadata.
    pub fn accept_tree_parent_change(
        &mut self,
        id: &str,
        folder: &store::FolderId,
        expected_current_generation: u64,
    ) -> TreeResult<()> {
        let (mut job, digest) = self.load_tree_job(id)?;
        if !job.local_folder_generations.contains_key(folder.as_str()) {
            return Err(FolderError::Invalid("folder not in clone job"));
        }
        if self.db.folder(&self.access, folder)?.generation != expected_current_generation {
            return Err(FolderError::Repository(Error::Stale));
        }
        job.local_folder_generations
            .insert(folder.as_str().into(), expected_current_generation);
        self.save_tree_job(&job, digest)
    }
    fn tree_parent(&self, job: &TreeJob, node: &TreeNode) -> TreeResult<Option<store::FolderId>> {
        let parent = match node {
            TreeNode::Folder { folder, .. } => {
                if matches!(&job.root,TreeRoot::Folder{id,..}if id==&folder.id) {
                    return Ok(job.destination.clone());
                }
                &folder.parent
            }
            TreeNode::Document { document, .. } => &document.folder,
        };
        parent
            .as_ref()
            .map(|id| {
                job.folders
                    .get(id.as_str())
                    .cloned()
                    .ok_or(FolderError::Invalid("unmapped/cyclic parent"))
            })
            .transpose()
            .map(|p| p.or_else(|| job.destination.clone()))
    }
    fn import_tree_node<T: RemoteTransport>(
        &mut self,
        job: &mut TreeJob,
        digest: store::Digest,
        node: TreeNode,
        transport: &mut T,
    ) -> TreeResult<()> {
        self.reload()?;
        let parent = self.tree_parent(job, &node)?;
        if let Some(id) = &parent {
            let current = self.db.folder(&self.access, id)?;
            if job.local_folder_generations.get(id.as_str()) != Some(&current.generation) {
                job.phase = TreePhase::Blocked;
                job.issue = Some(TreeIssue::LocalParentChanged { folder: id.clone() });
                return self.save_tree_job(job, digest);
            }
        }
        let rid = job.remote.id.clone();
        let operation = match &node {
            TreeNode::Folder { folder, .. } => {
                if self
                    .state
                    .folders
                    .values()
                    .filter_map(|m| m.get(&rid))
                    .any(|b| b.remote.id == folder.id)
                {
                    return Err(FolderError::AlreadyBound);
                }
                durable::MutationOperation::CreateFolder {
                    parent,
                    name: folder.name.to_string(),
                }
            }
            TreeNode::Document { document, .. } => {
                if self
                    .state
                    .documents
                    .values()
                    .filter_map(|m| m.get(&rid))
                    .any(|b| b.remote.id == document.id)
                {
                    return Err(FolderError::AlreadyBound);
                }
                let payload = exact_payload(job, document, transport)?;
                durable::MutationOperation::ImportDocument {
                    folder: parent,
                    name: document.name.to_string(),
                    payload,
                }
            }
        };
        let native = match &operation {
            durable::MutationOperation::ImportDocument { payload, .. } => Some(selected(payload)?),
            _ => None,
        };
        let mutation = durable::Mutation {
            id: self.db.new_mutation_id()?,
            operation,
        };
        let mut next_state = self.state.clone();
        let mut next_job = job.clone();
        let expected_state = self.checkpoint.clone();
        let result = self
            .db
            .apply_reliable_with(&self.access, &mutation, |receipt| {
                match (&node, &receipt.outcome) {
                    (
                        TreeNode::Folder { folder, .. },
                        durable::MutationOutcome::Applied(durable::Applied::CreatedFolder(local)),
                    ) => {
                        next_job
                            .folders
                            .insert(folder.id.as_str().into(), local.id.clone());
                        next_job
                            .local_folder_generations
                            .insert(local.id.as_str().into(), local.generation);
                        next_state
                            .folders
                            .entry(local.id.as_str().into())
                            .or_default()
                            .insert(
                                rid.clone(),
                                FolderBinding {
                                    remote: folder.clone(),
                                },
                            );
                    }
                    (
                        TreeNode::Document { document, .. },
                        durable::MutationOutcome::Applied(durable::Applied::Imported(local)),
                    ) => {
                        next_job
                            .documents
                            .insert(document.id.as_str().into(), local.id.clone());
                        next_state
                            .documents
                            .entry(local.id.as_str().into())
                            .or_default()
                            .insert(
                                rid.clone(),
                                DocumentBinding {
                                    remote: document.clone(),
                                    common_native: native.clone().ok_or_else(|| {
                                        durable::MutationError::Invalid("tree native pin")
                                    })?,
                                    common_local_revision: local.head.id.clone(),
                                },
                            );
                    }
                    (_, durable::MutationOutcome::Rejected(_)) => {
                        next_job.phase = TreePhase::Blocked;
                        next_job.issue = Some(TreeIssue::LocalMutationRejected);
                    }
                    _ => return Err(durable::MutationError::Invalid("tree local receipt")),
                }
                if matches!(receipt.outcome, durable::MutationOutcome::Applied(_)) {
                    next_state.epoch = next_state
                        .epoch
                        .checked_add(1)
                        .ok_or_else(|| durable::MutationError::Invalid("repository epoch"))?;
                }
                Ok(vec![
                    state_write(&next_state, expected_state).map_err(mutation_error)?,
                    job_write(&next_job, Some(digest))
                        .map_err(|_| durable::MutationError::Invalid("tree checkpoint"))?,
                ])
            })?;
        let _ = result;
        *job = next_job;
        self.reload()?;
        Ok(())
    }
    /// Each step performs one bounded scan page, local folder creation or pinned
    /// immutable document import. A dropped/failed transport keeps progress.
    pub fn resume_tree_clone<T: RemoteTransport>(
        &mut self,
        id: &str,
        max_steps: usize,
        transport: &mut T,
    ) -> TreeResult<TreeProgress> {
        if max_steps == 0 || max_steps > 1000 {
            return Err(FolderError::Invalid("step budget"));
        }
        self.reload()?;
        for _ in 0..max_steps {
            let (mut job, digest) = self.load_tree_job(id)?;
            if self.state.remotes.get(&job.remote.id) != Some(&job.remote) {
                return Err(FolderError::Repository(Error::Stale));
            }
            if matches!(
                job.phase,
                TreePhase::Complete | TreePhase::Invalidated | TreePhase::Blocked
            ) {
                return Ok(job.progress());
            }
            job.issue = None;
            if matches!(job.phase, TreePhase::Discovering | TreePhase::Scanning) {
                if let Err(error) = job.scan_step(transport) {
                    job.issue = Some(issue(&error));
                    job.phase = match error {
                        FolderError::Unavailable => job.phase,
                        FolderError::CatalogChanged { .. }
                        | FolderError::Rejected(wire::ErrorCode::NotFound) => {
                            TreePhase::Invalidated
                        }
                        _ => TreePhase::Blocked,
                    };
                    self.save_tree_job(&job, digest)?;
                    return Ok(job.progress());
                }
                self.save_tree_job(&job, digest)?;
                continue;
            }
            let node=job.nodes.iter().find(|n|matches!(n,TreeNode::Folder{folder,..}if!job.folders.contains_key(folder.id.as_str()))).or_else(||job.nodes.iter().find(|n|matches!(n,TreeNode::Document{document,..}if!job.documents.contains_key(document.id.as_str())))).cloned();
            let Some(node) = node else {
                job.phase = TreePhase::Complete;
                self.save_tree_job(&job, digest)?;
                return Ok(job.progress());
            };
            if let Err(error) = self.import_tree_node(&mut job, digest.clone(), node, transport) {
                job.issue = Some(issue(&error));
                if !matches!(error, FolderError::Unavailable) {
                    job.phase = TreePhase::Blocked;
                }
                self.save_tree_job(&job, digest)?;
                return Ok(job.progress());
            }
        }
        Ok(self.tree_clone_job(id)?.progress())
    }
    /// Queue one explicit document metadata change from a cloned tree. A virtual
    /// remote root mounted in a local folder remains the remote root; the local
    /// destination does not need (and must not invent) a remote folder binding.
    pub fn queue_tree_document_metadata(
        &mut self,
        job_id: &str,
        local_id: &store::DocumentId,
    ) -> TreeResult<String> {
        let (job, _) = self.load_tree_job(job_id)?;
        self.reload()?;
        if !job.documents.values().any(|id| id == local_id) {
            return Err(FolderError::Invalid("document outside clone"));
        }
        let local = self.db.document(&self.access, local_id)?;
        let binding = self
            .document_binding(local_id, &job.remote.id)
            .ok_or(FolderError::Invalid("unbound tree document"))?;
        // Rename is independent of placement. In particular, an unbound local
        // destination must not prevent a harmless rename or imply a remote move.
        let operation = if local.name != binding.remote.name.as_str() {
            sync::MutationOperation::RenameDocument {
                id: binding.remote.id.clone(),
                expected_generation: binding.remote.generation,
                name: wire::Text::new(&local.name)
                    .map_err(|_| FolderError::Invalid("document name"))?,
            }
        } else {
            let folder = if matches!(job.root, TreeRoot::LabRoot) && local.folder == job.destination
            {
                None
            } else {
                local
                    .folder
                    .as_ref()
                    .map(|id| {
                        self.folder_binding(id, &job.remote.id)
                            .map(|b| b.remote.id.clone())
                            .ok_or(FolderError::Repository(Error::UnboundFolder))
                    })
                    .transpose()?
            };
            if folder == binding.remote.folder {
                return Err(FolderError::Repository(Error::NoChanges));
            }
            sync::MutationOperation::MoveDocument {
                id: binding.remote.id.clone(),
                expected_generation: binding.remote.generation,
                folder,
            }
        };
        Ok(self.enqueue(
            job.remote,
            Target::Document(local_id.as_str().into()),
            operation,
            None,
            None,
        )?)
    }
    /// Mounted subtree placement is local bookkeeping, not an implicit remote
    /// move. Actual rename/move intent uses the acknowledged remote generation.
    pub fn queue_tree_folder_metadata(
        &mut self,
        job_id: &str,
        local_id: &store::FolderId,
    ) -> TreeResult<String> {
        let (job, _) = self.load_tree_job(job_id)?;
        self.reload()?;
        let local = self.db.folder(&self.access, local_id)?;
        let binding = self
            .folder_binding(local_id, &job.remote.id)
            .ok_or(FolderError::Invalid("unbound tree folder"))?;
        if !job.folders.values().any(|id| id == local_id) {
            return Err(FolderError::Invalid("folder outside clone"));
        }
        let mounted_root = matches!(&job.root,TreeRoot::Folder{id,..}if job.folders.get(id.as_str())==Some(local_id));
        let parent = if mounted_root && local.parent == job.destination {
            binding.remote.parent.clone()
        } else if matches!(job.root, TreeRoot::LabRoot) && local.parent == job.destination {
            None
        } else {
            local
                .parent
                .as_ref()
                .map(|id| {
                    self.folder_binding(id, &job.remote.id)
                        .map(|b| b.remote.id.clone())
                        .ok_or(FolderError::Repository(Error::UnboundFolder))
                })
                .transpose()?
        };
        let operation = if local.name != binding.remote.name.as_str() {
            sync::MutationOperation::RenameFolder {
                id: binding.remote.id.clone(),
                expected_generation: binding.remote.generation,
                name: wire::Text::new(&local.name)
                    .map_err(|_| FolderError::Invalid("folder name"))?,
            }
        } else if parent != binding.remote.parent {
            sync::MutationOperation::MoveFolder {
                id: binding.remote.id.clone(),
                expected_generation: binding.remote.generation,
                parent,
            }
        } else {
            return Err(FolderError::Repository(Error::NoChanges));
        };
        Ok(self.enqueue(
            job.remote,
            Target::Folder(local_id.as_str().into()),
            operation,
            None,
            None,
        )?)
    }
}
