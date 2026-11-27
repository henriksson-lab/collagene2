use super::*;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    Library,
    Sync,
    Discussion,
}
#[derive(Debug, thiserror::Error)]
pub enum TransportError {
    #[error("remote unavailable or response interrupted")]
    Unavailable,
}
/// Credentials/HTTP/Tauri IPC belong to this injected adapter, never to RemoteConfig
/// or native archives. Calls are synchronous so hosts can run them off the UI thread.
pub trait RemoteTransport {
    fn exchange(
        &mut self,
        remote: &RemoteConfig,
        channel: Channel,
        request: &[u8],
    ) -> std::result::Result<Vec<u8>, TransportError>;
}
#[derive(Debug, Clone)]
pub enum SyncResult {
    Applied(sync::Receipt),
    Rejected(sync::Receipt),
}
fn chunk_key(id: &str, index: usize) -> durable::CheckpointKey {
    durable::CheckpointKey::repository(format!("outbox-{id}-{index}"))
}
impl Repository {
    fn remote(&self, id: &str) -> Result<RemoteConfig> {
        self.state
            .remotes
            .get(id)
            .cloned()
            .ok_or(Error::Invalid("unknown remote"))
    }
    fn request<T: RemoteTransport>(
        &self,
        remote: &RemoteConfig,
        request: wire::Request,
        transport: &mut T,
    ) -> Result<wire::Success> {
        let request = wire::RequestEnvelope {
            request_id: 1,
            request,
        };
        let bytes =
            wire::encode_request(&request).map_err(|_| Error::Invalid("library request"))?;
        let reply = transport
            .exchange(remote, Channel::Library, &bytes)
            .map_err(|e| Error::Remote(e.to_string()))?;
        let reply = wire::decode_response_for(&reply, &request)
            .map_err(|_| Error::Remote("Unbound library reply".into()))?;
        reply
            .validate_for_in_context(&request, &remote.lab)
            .map_err(|_| Error::Remote("Reply from a different lab".into()))?;
        match reply.outcome {
            wire::Outcome::Success(s) => Ok(s),
            wire::Outcome::Failure(e) => Err(Error::Remote(e.message.into_string())),
        }
    }
    fn folder_destination(
        &self,
        parent: Option<&store::FolderId>,
        remote: &str,
    ) -> Result<Option<wire::FolderId>> {
        parent
            .map(|id| {
                self.folder_binding(id, remote)
                    .map(|b| b.remote.id.clone())
                    .ok_or(Error::UnboundFolder)
            })
            .transpose()
    }
    fn busy(&self, target: &Target, remote: &str) -> bool {
        self.state.pending.values().any(|p| {
            p.remote.id == remote
                && match (&p.target, target) {
                    (
                        Target::Folder(a) | Target::DeleteFolder(a),
                        Target::Folder(b) | Target::DeleteFolder(b),
                    ) => a == b,
                    (Target::Document(a), Target::Document(b)) => a == b,
                    _ => false,
                }
        })
    }
    pub(crate) fn enqueue(
        &mut self,
        remote: RemoteConfig,
        target: Target,
        operation: sync::MutationOperation,
        native: Option<NativePin>,
        local_revision: Option<store::RevisionId>,
    ) -> Result<String> {
        if self.busy(&target, &remote.id) {
            return Err(Error::Pending);
        }
        let id = sync::MutationId::new(self.db.new_mutation_id()?.as_str())
            .map_err(|_| Error::Invalid("mutation ID"))?;
        let mutation = sync::Mutation {
            id: id.clone(),
            operation,
        };
        let request = sync::RequestEnvelope {
            request_id: 1,
            expected_lab: remote.lab.clone(),
            request: sync::Request::Apply { mutation },
        };
        let frame = sync::encode_request(&request).map_err(|_| Error::Invalid("outbox bounds"))?;
        let chunks = frame
            .chunks(durable::MAX_CHECKPOINT_BYTES)
            .collect::<Vec<_>>();
        if chunks.len() > 2 {
            return Err(Error::Invalid("outbox chunks"));
        }
        let mut writes = vec![];
        for (index, bytes) in chunks.iter().enumerate() {
            writes.push(durable::CheckpointWrite {
                key: chunk_key(id.as_str(), index),
                expected: None,
                value: Some(durable::CheckpointValue {
                    version: 1,
                    bytes: bytes.to_vec(),
                }),
            });
        }
        let pending = Pending {
            id: id.clone(),
            remote,
            target,
            chunks: chunks.len(),
            frame_hash: wire::Digest::of(&frame),
            native,
            local_revision,
            attempted: false,
            rejection: None,
        };
        let mut next = self.state.clone();
        next.pending.insert(id.as_str().into(), pending);
        self.commit(next, writes)?;
        Ok(id.as_str().into())
    }
    /// Queue immutable publication bytes before touching the transport. Private
    /// unrelated documents, local UI and repository checkpoints are not sent.
    pub fn queue_push(&mut self, id: &store::DocumentId, remote_id: &str) -> Result<String> {
        self.queue_push_inner(id, remote_id, None)
    }
    /// Bind a native UI push to the saved local head it actually displayed.
    /// The service captures its checkpoint internally; clients supply no authority.
    pub fn queue_push_pinned(
        &mut self,
        id: &store::DocumentId,
        expected_local: &store::RevisionId,
        remote_id: &str,
        expected_checkpoint: &Option<store::Digest>,
    ) -> Result<String> {
        self.reload()?;
        if &self.checkpoint != expected_checkpoint {
            return Err(Error::Stale);
        }
        self.queue_push_inner(id, remote_id, Some(expected_local))
    }
    fn queue_push_inner(
        &mut self,
        id: &store::DocumentId,
        remote_id: &str,
        expected_local: Option<&store::RevisionId>,
    ) -> Result<String> {
        let remote = self.remote(remote_id)?;
        let doc = self.db.document(&self.access, id)?;
        if expected_local.is_some_and(|expected| expected != &doc.head.id) {
            return Err(Error::Stale);
        }
        let payload = self.db.payload(&self.access, id, &doc.head.id)?;
        let (payload, native) = publish(payload)?;
        self.queue_prepared_push(doc, remote, payload, native, false)
    }
    /// Internal prepared-publication path. The reviewed adapter has already
    /// rechecked its opaque source/remote pins and derived exact public scope.
    pub(super) fn queue_prepared_push(
        &mut self,
        doc: store::Document,
        remote: RemoteConfig,
        payload: wire::Payload,
        native: NativePin,
        documentary_update: bool,
    ) -> Result<String> {
        payload.validate().map_err(|_| Error::Native)?;
        if payload.kind == wire::PayloadKind::WorkspaceArchive
            && payload.sharing != Some(wire::WorkspaceSharing::DocumentClosure)
        {
            return Err(Error::Native);
        }
        let id = &doc.id;
        let remote_id = &remote.id;
        let operation = if let Some(binding) = self.document_binding(id, remote_id) {
            if binding.common_native == native
                && (!documentary_update || binding.remote.head.archive_hash == payload.archive_hash)
            {
                return Err(Error::NoChanges);
            }
            sync::MutationOperation::SaveRevision {
                id: binding.remote.id.clone(),
                base_revision: binding.remote.head.id.clone(),
                payload,
            }
        } else {
            sync::MutationOperation::ImportDocument {
                folder: self.folder_destination(doc.folder.as_ref(), remote_id)?,
                name: wire::Text::new(&doc.name).map_err(|_| Error::Invalid("document name"))?,
                payload,
            }
        };
        self.enqueue(
            remote,
            Target::Document(id.as_str().into()),
            operation,
            Some(native),
            Some(doc.head.id),
        )
    }
    /// Explicit one-step metadata push; rename precedes move when both changed.
    pub fn queue_document_metadata(
        &mut self,
        id: &store::DocumentId,
        remote_id: &str,
    ) -> Result<String> {
        let remote = self.remote(remote_id)?;
        let local = self.db.document(&self.access, id)?;
        let b = self
            .document_binding(id, remote_id)
            .ok_or(Error::Invalid("unbound document"))?;
        let folder = self.folder_destination(local.folder.as_ref(), remote_id)?;
        let operation = if local.name != b.remote.name.as_str() {
            sync::MutationOperation::RenameDocument {
                id: b.remote.id.clone(),
                expected_generation: b.remote.generation,
                name: wire::Text::new(&local.name).map_err(|_| Error::Invalid("document name"))?,
            }
        } else if folder != b.remote.folder {
            sync::MutationOperation::MoveDocument {
                id: b.remote.id.clone(),
                expected_generation: b.remote.generation,
                folder,
            }
        } else {
            return Err(Error::NoChanges);
        };
        self.enqueue(
            remote,
            Target::Document(id.as_str().into()),
            operation,
            None,
            None,
        )
    }
    pub fn queue_folder(&mut self, id: &store::FolderId, remote_id: &str) -> Result<String> {
        let remote = self.remote(remote_id)?;
        let local = self.db.folder(&self.access, id)?;
        let parent = self.folder_destination(local.parent.as_ref(), remote_id)?;
        let operation = if let Some(b) = self.folder_binding(id, remote_id) {
            if local.name != b.remote.name.as_str() {
                sync::MutationOperation::RenameFolder {
                    id: b.remote.id.clone(),
                    expected_generation: b.remote.generation,
                    name: wire::Text::new(&local.name)
                        .map_err(|_| Error::Invalid("folder name"))?,
                }
            } else if parent != b.remote.parent {
                sync::MutationOperation::MoveFolder {
                    id: b.remote.id.clone(),
                    expected_generation: b.remote.generation,
                    parent,
                }
            } else {
                return Err(Error::NoChanges);
            }
        } else {
            sync::MutationOperation::CreateFolder {
                parent,
                name: wire::Text::new(&local.name).map_err(|_| Error::Invalid("folder name"))?,
            }
        };
        self.enqueue(
            remote,
            Target::Folder(id.as_str().into()),
            operation,
            None,
            None,
        )
    }
    pub fn queue_deleted_folder(
        &mut self,
        id: &store::FolderId,
        remote_id: &str,
    ) -> Result<String> {
        match self.db.folder(&self.access, id) {
            Err(store::Error::NotFound) => {}
            Ok(_) => return Err(Error::Invalid("delete local empty folder explicitly first")),
            Err(e) => return Err(e.into()),
        }
        let remote = self.remote(remote_id)?;
        let b = self
            .folder_binding(id, remote_id)
            .ok_or(Error::Invalid("unbound folder"))?;
        let op = sync::MutationOperation::DeleteEmptyFolder {
            id: b.remote.id.clone(),
            expected_generation: b.remote.generation,
        };
        self.enqueue(
            remote,
            Target::DeleteFolder(id.as_str().into()),
            op,
            None,
            None,
        )
    }
    fn pending_mutation(
        &mut self,
        p: &Pending,
    ) -> Result<(sync::Mutation, Vec<(durable::CheckpointKey, store::Digest)>)> {
        if p.chunks == 0 || p.chunks > 2 {
            return Err(Error::Invalid("outbox chunks"));
        }
        let mut frame = vec![];
        let mut checkpoints = vec![];
        for index in 0..p.chunks {
            let c = self
                .db
                .read_checkpoint(&self.access, &chunk_key(p.id.as_str(), index))?
                .ok_or(Error::Invalid("missing outbox chunk"))?;
            if c.value.version != 1
                || frame
                    .len()
                    .checked_add(c.value.bytes.len())
                    .is_none_or(|n| n > wire::MAX_MESSAGE_BYTES + 42)
            {
                return Err(Error::Invalid("outbox frame size/version"));
            }
            frame.extend_from_slice(&c.value.bytes);
            checkpoints.push((c.key, c.digest));
        }
        if wire::Digest::of(&frame) != p.frame_hash {
            return Err(Error::Invalid("outbox integrity"));
        }
        let request = sync::decode_request(&frame).map_err(|_| Error::Invalid("outbox framing"))?;
        if request.expected_lab != p.remote.lab {
            return Err(Error::Invalid("outbox expected lab changed"));
        }
        let sync::Request::Apply { mutation } = request.request else {
            return Err(Error::Invalid("outbox operation"));
        };
        if mutation.id != p.id {
            return Err(Error::Invalid("outbox mutation ID"));
        }
        Ok((mutation, checkpoints))
    }
    fn sync_request<T: RemoteTransport>(
        &self,
        remote: &RemoteConfig,
        request: sync::Request,
        transport: &mut T,
    ) -> Result<sync::ResponseOutcome> {
        let request = sync::RequestEnvelope {
            request_id: 1,
            expected_lab: remote.lab.clone(),
            request,
        };
        let bytes = sync::encode_request(&request).map_err(|_| Error::Invalid("sync request"))?;
        let response = transport
            .exchange(remote, Channel::Sync, &bytes)
            .map_err(|e| Error::Remote(e.to_string()))?;
        let response = sync::decode_response_for(&response, &request)
            .map_err(|_| Error::Remote("Unbound sync response".into()))?;
        response
            .validate_for_in_context(&request, &remote.lab)
            .map_err(|_| Error::Remote("Sync reply from a different lab".into()))?;
        Ok(response.outcome)
    }
    /// Lookup first after every restart or interrupted reply. Unknown resends the
    /// exact immutable ID+intent; an existing receipt is never executed twice.
    pub fn resume_push<T: RemoteTransport>(
        &mut self,
        id: &str,
        transport: &mut T,
    ) -> Result<SyncResult> {
        let pending = self
            .state
            .pending
            .get(id)
            .cloned()
            .ok_or(Error::Invalid("unknown pending mutation"))?;
        if self.remote(&pending.remote.id)? != pending.remote {
            return Err(Error::Stale);
        }
        if let Some(receipt) = &pending.rejection {
            return Ok(SyncResult::Rejected(receipt.clone()));
        }
        let (mutation, chunks) = self.pending_mutation(&pending)?;
        let mut attempted = self.state.clone();
        attempted.pending.get_mut(id).unwrap().attempted = true;
        self.commit(attempted, vec![])?;
        let lookup = self.sync_request(
            &pending.remote,
            sync::Request::Lookup {
                mutation_id: mutation.id.clone(),
                expected_intent: mutation
                    .intent_digest()
                    .map_err(|_| Error::Invalid("mutation intent"))?,
            },
            transport,
        )?;
        let result = match lookup {
            sync::ResponseOutcome::Unknown { .. } => self.sync_request(
                &pending.remote,
                sync::Request::Apply {
                    mutation: mutation.clone(),
                },
                transport,
            )?,
            other => other,
        };
        let receipt = match result {
            sync::ResponseOutcome::Receipt(r) => r,
            sync::ResponseOutcome::IntentMismatch { .. } => {
                return Err(Error::Remote(
                    "Mutation identity belongs to another intent".into(),
                ));
            }
            sync::ResponseOutcome::Failure(f) => {
                return Err(Error::Remote(f.message.into_string()));
            }
            _ => {
                return Err(Error::Remote(
                    "Mutation outcome is unresolved; queued bytes retained".into(),
                ));
            }
        };
        receipt
            .validate_for(&mutation)
            .map_err(|_| Error::Remote("Receipt does not match queued mutation".into()))?;
        let mut next = self.state.clone();
        if matches!(receipt.outcome, sync::MutationOutcome::Rejected(_)) {
            next.pending.get_mut(id).ok_or(Error::Stale)?.rejection = Some(receipt.clone());
            self.commit(next, vec![])?;
            return Ok(SyncResult::Rejected(receipt));
        }
        if !next.pending.contains_key(id) {
            return Err(Error::Stale);
        }
        let sync::MutationOutcome::Applied(applied) = &receipt.outcome else {
            unreachable!()
        };
        match (&pending.target, applied) {
            (
                Target::Document(local),
                sync::Applied::Imported(remote)
                | sync::Applied::Saved {
                    document: remote, ..
                },
            ) => {
                let native = pending
                    .native
                    .clone()
                    .ok_or(Error::Invalid("outbox native pin"))?;
                let local_revision = pending
                    .local_revision
                    .clone()
                    .ok_or(Error::Invalid("outbox source revision"))?;
                next.documents.entry(local.clone()).or_default().insert(
                    pending.remote.id.clone(),
                    DocumentBinding {
                        remote: remote.clone(),
                        common_native: native,
                        common_local_revision: local_revision,
                    },
                );
            }
            (
                Target::Document(local),
                sync::Applied::RenamedDocument(remote) | sync::Applied::MovedDocument(remote),
            ) => {
                next.documents
                    .get_mut(local)
                    .and_then(|m| m.get_mut(&pending.remote.id))
                    .ok_or(Error::Stale)?
                    .remote = remote.clone();
            }
            (
                Target::Folder(local),
                sync::Applied::CreatedFolder(remote)
                | sync::Applied::RenamedFolder(remote)
                | sync::Applied::MovedFolder(remote),
            ) => {
                next.folders.entry(local.clone()).or_default().insert(
                    pending.remote.id.clone(),
                    FolderBinding {
                        remote: remote.clone(),
                    },
                );
            }
            (Target::DeleteFolder(local), sync::Applied::DeletedEmptyFolder(_)) => {
                next.folders
                    .get_mut(local)
                    .ok_or(Error::Stale)?
                    .remove(&pending.remote.id);
            }
            _ => {
                return Err(Error::Remote(
                    "Receipt operation differs from outbox target".into(),
                ));
            }
        }
        next.pending.remove(id);
        let writes = chunks
            .into_iter()
            .map(|(key, digest)| durable::CheckpointWrite {
                key,
                expected: Some(digest),
                value: None,
            })
            .collect();
        self.commit(next, writes)?;
        Ok(SyncResult::Applied(receipt))
    }
    /// Before the first transport attempt, cancellation cannot hide a remote
    /// effect. Once any attempt starts, resolve the durable outcome instead.
    pub fn cancel_unattempted(&mut self, id: &str) -> Result<()> {
        let pending = self
            .state
            .pending
            .get(id)
            .cloned()
            .ok_or(Error::Invalid("unknown outbox"))?;
        if pending.attempted {
            return Err(Error::Invalid("resolve attempted remote outcome first"));
        }
        let (_, chunks) = self.pending_mutation(&pending)?;
        let writes = chunks
            .into_iter()
            .map(|(key, digest)| durable::CheckpointWrite {
                key,
                expected: Some(digest),
                value: None,
            })
            .collect();
        let mut next = self.state.clone();
        next.pending.remove(id);
        self.commit(next, writes)
    }
    /// Only an acknowledged rejection may be discarded. An interrupted create
    /// cannot be forgotten then reissued under a new ID and duplicate its effect.
    pub fn discard_rejected(&mut self, id: &str) -> Result<()> {
        let p = self
            .state
            .pending
            .get(id)
            .cloned()
            .ok_or(Error::Invalid("unknown outbox"))?;
        if p.rejection.is_none() {
            return Err(Error::Invalid(
                "resolve unknown remote outcome before discarding",
            ));
        }
        let (_, chunks) = self.pending_mutation(&p)?;
        let writes = chunks
            .into_iter()
            .map(|(key, digest)| durable::CheckpointWrite {
                key,
                expected: Some(digest),
                value: None,
            })
            .collect();
        let mut next = self.state.clone();
        next.pending.remove(id);
        self.commit(next, writes)
    }
    fn remote_document<T: RemoteTransport>(
        &self,
        remote: &RemoteConfig,
        id: wire::DocumentId,
        transport: &mut T,
    ) -> Result<(wire::Document, wire::Payload)> {
        match self.request(
            remote,
            wire::Request::OpenDocument { id, revision: None },
            transport,
        )? {
            wire::Success::Opened {
                document, payload, ..
            } => Ok((document, payload)),
            _ => Err(Error::Remote("Unexpected document response".into())),
        }
    }
    pub fn clone_from<T: RemoteTransport>(
        &mut self,
        remote_id: &str,
        id: wire::DocumentId,
        folder: Option<store::FolderId>,
        transport: &mut T,
    ) -> Result<store::Document> {
        if self
            .state
            .documents
            .values()
            .filter_map(|m| m.get(remote_id))
            .any(|b| b.remote.id == id)
        {
            return Err(Error::AlreadyBound);
        }
        let remote = self.remote(remote_id)?;
        let (document, payload) = self.remote_document(&remote, id, transport)?;
        let payload = stored_payload(payload)?;
        let native = selected(&payload)?;
        let source = document.clone();
        let rid = remote_id.to_owned();
        let result = self.local(
            durable::MutationOperation::ImportDocument {
                folder,
                name: document.name.into_string(),
                payload,
            },
            move |state, applied| {
                let durable::Applied::Imported(local) = applied else {
                    return Err(Error::Invalid("local clone receipt"));
                };
                state
                    .documents
                    .entry(local.id.as_str().into())
                    .or_default()
                    .insert(
                        rid,
                        DocumentBinding {
                            remote: source,
                            common_native: native,
                            common_local_revision: local.head.id.clone(),
                        },
                    );
                Ok(())
            },
        )?;
        match result {
            durable::Applied::Imported(d) => Ok(d),
            _ => Err(Error::Invalid("clone receipt")),
        }
    }
    pub fn pull<T: RemoteTransport>(
        &mut self,
        id: &store::DocumentId,
        remote_id: &str,
        options: &ImportOptions,
        transport: &mut T,
    ) -> Result<store::Document> {
        self.pull_inner(id, remote_id, options, transport, None)
    }
    /// Check the displayed saved head before any remote request. A later local
    /// mutation is still rejected by the existing immutable-head/checkpoint CAS.
    pub fn pull_pinned<T: RemoteTransport>(
        &mut self,
        id: &store::DocumentId,
        expected_local: &store::RevisionId,
        remote_id: &str,
        options: &ImportOptions,
        transport: &mut T,
        expected_checkpoint: &Option<store::Digest>,
    ) -> Result<store::Document> {
        self.reload()?;
        if &self.checkpoint != expected_checkpoint {
            return Err(Error::Stale);
        }
        self.pull_inner(id, remote_id, options, transport, Some(expected_local))
    }
    fn pull_inner<T: RemoteTransport>(
        &mut self,
        id: &store::DocumentId,
        remote_id: &str,
        options: &ImportOptions,
        transport: &mut T,
        expected_local: Option<&store::RevisionId>,
    ) -> Result<store::Document> {
        if self.busy(&Target::Document(id.as_str().into()), remote_id) {
            return Err(Error::Pending);
        }
        let remote = self.remote(remote_id)?;
        let binding = self
            .document_binding(id, remote_id)
            .cloned()
            .ok_or(Error::Invalid("unbound document"))?;
        let current = self.db.document(&self.access, id)?;
        if expected_local.is_some_and(|expected| expected != &current.head.id) {
            return Err(Error::Stale);
        }
        let local = self.db.payload(&self.access, id, &current.head.id)?;
        let (remote_doc, incoming) = self.remote_document(&remote, binding.remote.id, transport)?;
        let incoming = stored_payload(incoming)?;
        let incoming_pin = selected(&incoming)?;
        let (local_frame, incoming_frame, selection) =
            match (selected(&local)?, incoming_pin.clone()) {
                (NativePin::Workspace(local_pin), NativePin::Workspace(remote_pin)) => {
                    if local_pin.document_id != remote_pin.document_id {
                        return Err(Error::NativeRootChanged);
                    }
                    let local_frame =
                        workspace::load_workspace(&local.bytes).map_err(|_| Error::Native)?;
                    let remote_frame =
                        workspace::load_workspace(&incoming.bytes).map_err(|_| Error::Native)?;
                    let narrow = sharing::project_document_closure(&remote_frame, &remote_pin)
                        .map_err(publication_error)?
                        .workspace;
                    (local_frame, narrow, remote_pin)
                }
                // Bare native records have no shared event graph to prove a merge.
                _ => return Err(Error::NativeRootChanged),
            };
        let merged = sharing::import::merge_document_closure(
            &local_frame,
            &incoming_frame,
            &selection,
            options,
        )?;
        let payload = private_workspace(&merged.workspace, &selection)?;
        let bookmark = merged.displaced_selected;
        let rid = remote_id.to_owned();
        let applied = self.local(
            durable::MutationOperation::SaveRevision {
                id: id.clone(),
                base_revision: current.head.id,
                payload,
            },
            move |state, applied| {
                let durable::Applied::Saved { document, .. } = applied else {
                    return Err(Error::Invalid("local pull receipt"));
                };
                if let Some(bookmark) = bookmark {
                    state.bookmarks.push(bookmark);
                }
                state
                    .documents
                    .entry(document.id.as_str().into())
                    .or_default()
                    .insert(
                        rid,
                        DocumentBinding {
                            remote: remote_doc,
                            common_native: incoming_pin,
                            common_local_revision: document.head.id.clone(),
                        },
                    );
                Ok(())
            },
        )?;
        match applied {
            durable::Applied::Saved { document, .. } => Ok(document),
            _ => Err(Error::Invalid("pull receipt")),
        }
    }
    /// Eventless .cg2 records cannot prove native event ancestry. An explicit
    /// snapshot checkout instead pins the current local library revision and
    /// atomically bookmarks that entire immutable body before changing the head.
    pub fn pull_snapshot<T: RemoteTransport>(
        &mut self,
        id: &store::DocumentId,
        remote_id: &str,
        expected_local: &store::RevisionId,
        transport: &mut T,
    ) -> Result<store::Document> {
        if self.busy(&Target::Document(id.as_str().into()), remote_id) {
            return Err(Error::Pending);
        }
        let current = self.db.document(&self.access, id)?;
        if &current.head.id != expected_local {
            return Err(Error::Stale);
        }
        let local = self.db.payload(&self.access, id, &current.head.id)?;
        if !matches!(selected(&local)?, NativePin::Document { .. }) {
            return Err(Error::Invalid(
                "use graph-preserving pull for workspace archives",
            ));
        }
        let remote = self.remote(remote_id)?;
        let binding = self
            .document_binding(id, remote_id)
            .cloned()
            .ok_or(Error::Invalid("unbound document"))?;
        let (remote_doc, incoming) = self.remote_document(&remote, binding.remote.id, transport)?;
        let incoming = stored_payload(incoming)?;
        let native = selected(&incoming)?;
        if !matches!(native, NativePin::Document { .. }) {
            return Err(Error::NativeRootChanged);
        }
        let bookmark = SnapshotBookmark {
            document: id.clone(),
            revision: current.head.id.clone(),
        };
        let rid = remote_id.to_owned();
        let applied = self.local(
            durable::MutationOperation::SaveRevision {
                id: id.clone(),
                base_revision: current.head.id,
                payload: incoming,
            },
            move |state, applied| {
                let durable::Applied::Saved { document, .. } = applied else {
                    return Err(Error::Invalid("local snapshot checkout receipt"));
                };
                state.snapshot_bookmarks.push(bookmark);
                state
                    .documents
                    .entry(document.id.as_str().into())
                    .or_default()
                    .insert(
                        rid,
                        DocumentBinding {
                            remote: remote_doc,
                            common_native: native,
                            common_local_revision: document.head.id.clone(),
                        },
                    );
                Ok(())
            },
        )?;
        match applied {
            durable::Applied::Saved { document, .. } => Ok(document),
            _ => Err(Error::Invalid("snapshot checkout receipt")),
        }
    }
    /// Clone only the selected folder metadata. Children are selected/cloned
    /// explicitly; this is deliberately not a recursive background sync engine.
    pub fn clone_folder<T: RemoteTransport>(
        &mut self,
        remote_id: &str,
        id: wire::FolderId,
        remote_parent: Option<wire::FolderId>,
        local_parent: Option<store::FolderId>,
        transport: &mut T,
    ) -> Result<store::Folder> {
        let remote = self.remote(remote_id)?;
        if self
            .state
            .folders
            .values()
            .filter_map(|m| m.get(remote_id))
            .any(|b| b.remote.id == id)
        {
            return Err(Error::AlreadyBound);
        }
        let mut page = wire::PageRequest::default();
        let found = loop {
            let result = self.request(
                &remote,
                wire::Request::Browse {
                    parent: remote_parent.clone(),
                    page: page.clone(),
                },
                transport,
            )?;
            let wire::Success::Browsed(result) = result else {
                return Err(Error::Remote("Unexpected folder listing".into()));
            };
            if let Some(folder) = result.entries.iter().find_map(|entry| match entry {
                wire::Entry::Folder(f) if f.id == id => Some(f.clone()),
                _ => None,
            }) {
                break folder;
            }
            let Some(offset) = result.next_offset else {
                return Err(Error::Remote(
                    "Selected remote folder is unavailable in this parent".into(),
                ));
            };
            page.offset = offset;
            page.expected_catalog_generation = Some(result.catalog_generation);
        };
        let rid = remote_id.to_owned();
        let applied = self.local(
            durable::MutationOperation::CreateFolder {
                parent: local_parent,
                name: found.name.as_str().into(),
            },
            move |state, applied| {
                let durable::Applied::CreatedFolder(local) = applied else {
                    return Err(Error::Invalid("folder clone receipt"));
                };
                state
                    .folders
                    .entry(local.id.as_str().into())
                    .or_default()
                    .insert(rid, FolderBinding { remote: found });
                Ok(())
            },
        )?;
        match applied {
            durable::Applied::CreatedFolder(f) => Ok(f),
            _ => Err(Error::Invalid("folder clone receipt")),
        }
    }
}
