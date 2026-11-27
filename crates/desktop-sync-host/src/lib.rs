//! Saved-catalog desktop remote service. Shell/window trust is supplied outside
//! this data protocol. No command accepts a path, archive upload or arbitrary
//! remote mutation, and local work remains independent of this service.
#![cfg(not(target_arch = "wasm32"))]
mod adapter;
mod projection;
mod publication;
mod removal;
pub use adapter::{ProductionAdapter, SyncAdapter};
use collagene_desktop_credentials as credentials;
use collagene_desktop_sync_protocol as wire;
pub use collagene_desktop_sync_protocol::{MAX_REQUEST_BYTES, MAX_RESPONSE_BYTES};
use collagene_lab_library as store;
use collagene_library_protocol as library;
use collagene_library_sharing::{
    self as sharing,
    import::{ImportOptions, SelectedPolicy},
};
use collagene_repository::{self as repository, RemoteConfig, Repository, folders};
use projection::convert;
#[derive(Debug, thiserror::Error)]
pub enum HostError {
    #[error("invalid desktop sync frame")]
    Wire,
    #[error("credential adapter does not belong to this native repository")]
    Scope,
    #[error("desktop sync service unavailable")]
    Unavailable,
}
#[derive(Debug)]
struct Fault {
    failure: wire::Failure,
    context: Option<wire::RemoteContext>,
}
type Result<T> = std::result::Result<T, Fault>;
impl Fault {
    fn error(code: wire::ErrorCode, message: &'static str) -> Self {
        Self {
            failure: wire::Failure::Error {
                code,
                message: library::Text::new(message).expect("bounded static native message"),
            },
            context: None,
        }
    }
    fn invalid(message: &'static str) -> Self {
        Self::error(wire::ErrorCode::Invalid, message)
    }
    fn too_large() -> Self {
        Self::error(
            wire::ErrorCode::TooLarge,
            "Native metadata exceeds this bounded page.",
        )
    }
    fn unavailable() -> Self {
        Self::error(
            wire::ErrorCode::Unavailable,
            "The native operation could not complete; existing pending outcomes are retained.",
        )
    }
}
fn repository_error(error: repository::Error) -> Fault {
    use repository::Error as E;
    match error {
        E::Pending => Fault::error(
            wire::ErrorCode::Pending,
            "Resolve the existing pending operation first.",
        ),
        E::NoChanges => Fault::error(
            wire::ErrorCode::NoChanges,
            "No saved native changes to publish.",
        ),
        E::UnboundFolder => Fault::error(
            wire::ErrorCode::UnboundFolder,
            "Bind or explicitly push the parent folder first.",
        ),
        E::AlreadyBound => Fault::error(
            wire::ErrorCode::AlreadyBound,
            "This remote item is already bound; pull it explicitly.",
        ),
        E::NativeRootChanged => Fault::error(
            wire::ErrorCode::NativeRootChanged,
            "The remote native root changed; clone separately to preserve local work.",
        ),
        E::Cancelled => Fault::error(
            wire::ErrorCode::Cancelled,
            "Operation cancelled before acceptance.",
        ),
        E::Merge(_) => Fault::error(
            wire::ErrorCode::HistoryConflict,
            "Saved native histories diverge; review an explicit checkout.",
        ),
        E::SbolSource => Fault::invalid(
            "SBOL source evidence cannot be shared in its current scope. Keep the complete native archive and review the original design.",
        ),
        E::Storage(store::Error::NotFound) => Fault::error(
            wire::ErrorCode::NotFound,
            "Item not found in this native repository.",
        ),
        E::Invalid(_) | E::Native => {
            Fault::invalid("The saved native state or selected operation is invalid.")
        }
        _ => Fault::unavailable(),
    }
}
fn folder_error(error: folders::FolderError) -> Fault {
    match error {
        folders::FolderError::Repository(e) => repository_error(e),
        folders::FolderError::Limit => Fault::too_large(),
        folders::FolderError::AlreadyBound => Fault::error(
            wire::ErrorCode::AlreadyBound,
            "An existing local tree binding was preserved.",
        ),
        folders::FolderError::Invalid(_) => {
            Fault::invalid("Invalid saved tree job or mount operation.")
        }
        _ => Fault::unavailable(),
    }
}
fn credential_error(error: credentials::Error) -> Fault {
    match error {
        credentials::Error::InvalidScope | credentials::Error::InvalidSecret => {
            Fault::invalid("Invalid native credential scope or token.")
        }
        credentials::Error::Missing => Fault::error(
            wire::ErrorCode::Unavailable,
            "No credential is available for this exact native repository and remote.",
        ),
        credentials::Error::Locked => Fault::error(
            wire::ErrorCode::Unavailable,
            "The native credential vault is locked; no session fallback was installed.",
        ),
        credentials::Error::Denied => Fault::error(
            wire::ErrorCode::Unavailable,
            "Native credential access was denied.",
        ),
        _ => Fault::unavailable(),
    }
}
fn accepted(gate: &mut impl FnMut() -> bool) -> Result<()> {
    if gate() {
        Ok(())
    } else {
        Err(Fault::error(
            wire::ErrorCode::Cancelled,
            "Operation cancelled before acceptance.",
        ))
    }
}
fn document(id: &library::DocumentId) -> Result<store::DocumentId> {
    store::DocumentId::parse(id.as_str())
        .map_err(|_| Fault::invalid("Invalid local document identity."))
}
fn folder(id: &library::FolderId) -> Result<store::FolderId> {
    store::FolderId::parse(id.as_str())
        .map_err(|_| Fault::invalid("Invalid local folder identity."))
}
fn revision(id: &library::RevisionId) -> Result<store::RevisionId> {
    store::RevisionId::parse(id.as_str())
        .map_err(|_| Fault::invalid("Invalid local library revision."))
}
/// The adapter must already be native-configured; a DSQ frame never selects a
/// repository/path or grants window/credential authority to a caller.
pub struct SyncHost<A> {
    repository: Repository,
    adapter: A,
    local_lab: library::LabId,
    reviews: publication::Reviews,
    removals: removal::Reviews,
}
impl SyncHost<ProductionAdapter> {
    pub fn production(
        repository: Repository,
        options: repository::http::HttpOptions,
    ) -> std::result::Result<Self, HostError> {
        let local = repository
            .local_context()
            .map_err(|_| HostError::Unavailable)?
            .lab_id;
        let adapter = repository::http::HttpRemoteTransport::new(
            credentials::Credentials::new(local),
            options,
        )
        .map_err(|_| HostError::Unavailable)?;
        Self::new(repository, adapter)
    }
}
impl<A: SyncAdapter> SyncHost<A> {
    pub fn new(repository: Repository, mut adapter: A) -> std::result::Result<Self, HostError> {
        let local = repository
            .local_context()
            .map_err(|_| HostError::Unavailable)?
            .lab_id;
        if adapter.local_repository_lab() != local {
            return Err(HostError::Scope);
        }
        Ok(Self {
            repository,
            adapter,
            local_lab: local,
            reviews: publication::Reviews::new()?,
            removals: removal::Reviews::new()?,
        })
    }
    pub fn adapter_mut(&mut self) -> &mut A {
        &mut self.adapter
    }
    pub fn repository(&self) -> &Repository {
        &self.repository
    }
    pub fn into_parts(self) -> (Repository, A) {
        (self.repository, self.adapter)
    }
    fn remote(&self, id: &wire::RemoteId) -> Result<RemoteConfig> {
        self.repository
            .remotes()
            .find(|r| r.id == id.as_str())
            .cloned()
            .ok_or_else(|| {
                Fault::error(
                    wire::ErrorCode::NotFound,
                    "Configured remote not found in this repository.",
                )
            })
    }
    fn pending(&self, id: &str) -> Result<wire::PendingPush> {
        let items = self
            .repository
            .outbox_metadata()
            .map_err(repository_error)?;
        projection::pending(items.iter().find(|p| p.id == id).ok_or_else(|| {
            Fault::error(wire::ErrorCode::NotFound, "Pending operation not found.")
        })?)
    }
    fn saved(
        &self,
        id: &store::DocumentId,
        expected: &store::RevisionId,
    ) -> Result<store::Document> {
        let current = self.repository.document(id).map_err(repository_error)?;
        if &current.head.id != expected {
            return Err(Fault {
                failure: wire::Failure::LocalConflict {
                    current: convert(&current)?,
                },
                context: None,
            });
        }
        Ok(current)
    }
    fn read_remote(
        &mut self,
        remote: &RemoteConfig,
        request: library::Request,
        request_id: u64,
    ) -> Result<(wire::RemoteContext, library::Success)> {
        if !matches!(
            request,
            library::Request::Browse { .. }
                | library::Request::Search { .. }
                | library::Request::ListRevisions { .. }
        ) {
            return Err(Fault::invalid(
                "Only native metadata reads are permitted here.",
            ));
        }
        let request = library::RequestEnvelope {
            request_id,
            request,
        };
        let bytes = library::encode_request(&request)
            .map_err(|_| Fault::invalid("Invalid remote metadata request."))?;
        let response = self
            .adapter
            .exchange(remote, repository::Channel::Library, &bytes)
            .map_err(|_| Fault::unavailable())?;
        let response =
            library::decode_response_for(&response, &request).map_err(|_| Fault::unavailable())?;
        response
            .validate_for_in_context(&request, &remote.lab)
            .map_err(|_| Fault::unavailable())?;
        let lab = response.context.ok_or_else(Fault::unavailable)?;
        let context = wire::RemoteContext {
            remote: projection::remote(remote)?,
            lab,
        };
        match response.outcome {
            library::Outcome::Success(s) => Ok((context, s)),
            library::Outcome::Failure(f) => Err(Fault {
                failure: wire::Failure::Remote(library::Failure {
                    code: f.code,
                    message: library::Text::new(
                        "The remote metadata request was rejected; refresh its current state.",
                    )
                    .expect("static public message"),
                    conflict: f.conflict,
                }),
                context: Some(context),
            }),
        }
    }
    fn context(&mut self, remote: &RemoteConfig, request_id: u64) -> Result<wire::RemoteContext> {
        self.read_remote(
            remote,
            library::Request::Browse {
                parent: None,
                page: library::PageRequest {
                    offset: 0,
                    limit: 1,
                    expected_catalog_generation: None,
                },
            },
            request_id,
        )
        .map(|v| v.0)
    }
    fn tree_status(&mut self, id: &str) -> Result<wire::TreeStatus> {
        let job = self.repository.tree_clone_job(id).map_err(folder_error)?;
        projection::tree(&job)
    }
    fn execute(
        &mut self,
        request: &wire::RequestEnvelope,
        gate: &mut impl FnMut() -> bool,
        context: &mut Option<wire::RemoteContext>,
    ) -> Result<wire::Success> {
        self.repository.reload().map_err(repository_error)?;
        if self
            .repository
            .local_context()
            .map_err(repository_error)?
            .lab_id
            != self.local_lab
        {
            return Err(Fault::unavailable());
        }
        let checkpoint = self.repository.checkpoint_version();
        if matches!(
            &request.command,
            wire::Command::SetCredential { .. }
                | wire::Command::DeleteCredential { .. }
                | wire::Command::ClearSession { .. }
        ) {
            self.reviews.clear();
            self.removals.clear();
        }
        if let Some(success) = self.execute_removal(&request.command, gate)? {
            return Ok(success);
        }
        if let Some(success) = self.execute_publication(&request.command, gate)? {
            return Ok(success);
        }
        use wire::{Command as C, Success as S};
        Ok(match &request.command {
            C::PrepareRemoteRemoval { .. }
            | C::ConfirmRemoteRemoval { .. }
            | C::DiscardRemoteRemovalReview { .. }
            | C::CheckRemoteRemovalOutcome { .. }
            | C::PrepareMetadataPublication { .. }
            | C::QueueMetadataPublication { .. }
            | C::InspectPushChoices { .. }
            | C::PrepareReviewedPush { .. }
            | C::QueueReviewedPush { .. }
            | C::DiscardPushReview { .. }
            | C::ListDeletedFolderPublications { .. } => {
                return Err(Fault::invalid("Publication command was not dispatched."));
            }
            C::Remotes => S::Remotes(
                library::Rows::new(
                    self.repository
                        .remotes()
                        .map(projection::remote)
                        .collect::<Result<Vec<_>>>()?,
                )
                .map_err(|_| Fault::too_large())?,
            ),
            C::AddRemote {
                name,
                endpoint,
                expected_lab,
            } => {
                accepted(gate)?;
                let remote = self
                    .repository
                    .add_remote(name.as_str(), endpoint.as_str(), expected_lab.clone())
                    .map_err(repository_error)?;
                S::AddedRemote(projection::remote(&remote)?)
            }
            C::SetCredential {
                remote,
                secret,
                retention,
            } => {
                let config = self.remote(remote)?;
                let secret = credentials::Secret::new(secret.expose_secret().as_bytes().to_vec())
                    .map_err(credential_error)?;
                let persistence = match retention {
                    wire::CredentialRetention::Session => credentials::Persistence::Session,
                    wire::CredentialRetention::Remember => credentials::Persistence::Remember,
                };
                accepted(gate)?;
                self.adapter
                    .store_credential(&config, secret, persistence)
                    .map_err(credential_error)?;
                S::CredentialSet {
                    remote: remote.clone(),
                    retention: *retention,
                }
            }
            C::DeleteCredential { remote } => {
                let config = self.remote(remote)?;
                accepted(gate)?;
                self.adapter
                    .delete_credential(&config)
                    .map_err(credential_error)?;
                S::CredentialDeleted {
                    remote: remote.clone(),
                }
            }
            C::CredentialStatus { remote } => {
                let config = self.remote(remote)?;
                let state = match self
                    .adapter
                    .credential_status(&config)
                    .map_err(credential_error)?
                {
                    credentials::CredentialStatus::Missing => wire::CredentialState::Missing,
                    credentials::CredentialStatus::Session => wire::CredentialState::Session,
                    credentials::CredentialStatus::Remembered => wire::CredentialState::Remembered,
                };
                S::CredentialStatus {
                    remote: remote.clone(),
                    state,
                }
            }
            C::ClearSession { remote } => {
                let config = remote.as_ref().map(|id| self.remote(id)).transpose()?;
                accepted(gate)?;
                self.adapter
                    .clear_session(config.as_ref())
                    .map_err(credential_error)?;
                S::SessionsCleared {
                    remote: remote.clone(),
                }
            }
            C::Discussion {
                remote,
                request: operation,
            } => {
                let config = self.remote(remote)?;
                let request = library::discussion::RequestEnvelope {
                    request_id: request.request_id,
                    request: operation.clone(),
                };
                let bytes = library::discussion::encode_request(&request)
                    .map_err(|_| Fault::invalid("Invalid revision discussion request."))?;
                // Every remote discussion action is explicit. It is never queued,
                // replayed or attached to a local scientific document.
                accepted(gate)?;
                let bytes = self
                    .adapter
                    .exchange(&config, repository::Channel::Discussion, &bytes)
                    .map_err(|_| Fault::unavailable())?;
                let reply = library::discussion::decode_response(&bytes)
                    .map_err(|_| Fault::unavailable())?;
                reply
                    .validate_for(&request, &config.lab)
                    .map_err(|_| Fault::unavailable())?;
                *context = Some(wire::RemoteContext {
                    remote: projection::remote(&config)?,
                    lab: reply.context.clone(),
                });
                S::Discussion(reply)
            }
            C::Browse { remote, .. }
            | C::Search { remote, .. }
            | C::ListRevisions { remote, .. } => {
                let config = self.remote(remote)?;
                let read = request
                    .command
                    .read_request()
                    .ok_or_else(|| Fault::invalid("Invalid metadata read."))?;
                let (proof, result) = self.read_remote(&config, read, request.request_id)?;
                *context = Some(proof);
                match result {
                    library::Success::Browsed(p) => S::Browsed(p),
                    library::Success::Searched(p) => S::Searched(p),
                    library::Success::Revisions(p) => S::Revisions(p),
                    _ => return Err(Fault::unavailable()),
                }
            }
            C::CloneDocument {
                remote,
                document: source,
                folder: destination,
            } => {
                let config = self.remote(remote)?;
                *context = Some(self.context(&config, request.request_id)?);
                let destination = destination.as_ref().map(folder).transpose()?;
                accepted(gate)?;
                let local = self
                    .repository
                    .clone_from(
                        remote.as_str(),
                        source.clone(),
                        destination,
                        &mut self.adapter,
                    )
                    .map_err(repository_error)?;
                let binding = self
                    .repository
                    .document_binding(&local.id, remote.as_str())
                    .ok_or_else(Fault::unavailable)?;
                S::Cloned {
                    source_document: binding.remote.id.clone(),
                    source_revision: binding.remote.head.id.clone(),
                    local: convert(&local)?,
                }
            }
            C::QueuePush {
                local_document,
                expected_local_revision,
                remote,
            } => {
                self.remote(remote)?;
                let id = document(local_document)?;
                let expected = revision(expected_local_revision)?;
                self.saved(&id, &expected)?;
                accepted(gate)?;
                let queued = self
                    .repository
                    .queue_push_pinned(&id, &expected, remote.as_str(), &checkpoint)
                    .map_err(repository_error)?;
                S::Queued(self.pending(&queued)?)
            }
            C::ResumePush { job } => {
                self.pending(job.as_str())?;
                accepted(gate)?;
                let result = self
                    .repository
                    .resume_push(job.as_str(), &mut self.adapter)
                    .map_err(repository_error)?;
                let receipt = match result {
                    repository::SyncResult::Applied(r) | repository::SyncResult::Rejected(r) => r,
                };
                S::Resumed(projection::receipt(&receipt)?)
            }
            C::Pending { page } => {
                let items = self
                    .repository
                    .outbox_metadata()
                    .map_err(repository_error)?
                    .iter()
                    .map(projection::pending)
                    .collect::<Result<Vec<_>>>()?;
                let (items, next_offset) = projection::page(items, page)?;
                S::Pending(wire::PendingPage {
                    items: library::Rows::new(items).map_err(|_| Fault::too_large())?,
                    next_offset,
                })
            }
            C::ListBookmarks {
                page,
                expected_checkpoint,
            } => {
                let checkpoint = self
                    .repository
                    .checkpoint_version()
                    .as_ref()
                    .map(convert)
                    .transpose()?;
                if expected_checkpoint
                    .as_ref()
                    .is_some_and(|expected| checkpoint.as_ref() != Some(expected))
                {
                    return Err(Fault::error(
                        wire::ErrorCode::BookmarkChanged,
                        "The local repository changed; refresh retained branches before continuing.",
                    ));
                }
                let mut items = Vec::new();
                for bookmark in self.repository.bookmarks() {
                    let old = &bookmark.document;
                    items.push(wire::Bookmark::NativeBranch {
                        document: library::Text::new(&old.id).map_err(|_| Fault::too_large())?,
                        revision: library::Digest::new(&old.head)
                            .map_err(|_| Fault::invalid("Invalid retained native revision."))?,
                        event: old
                            .head_event
                            .as_ref()
                            .map(library::Text::new)
                            .transpose()
                            .map_err(|_| Fault::too_large())?,
                        undo_count: u32::try_from(old.undo.len())
                            .map_err(|_| Fault::too_large())?,
                        redo_count: u32::try_from(old.redo.len())
                            .map_err(|_| Fault::too_large())?,
                        incoming_revision: library::Digest::new(
                            &bookmark.imported_selection.revision,
                        )
                        .map_err(|_| Fault::invalid("Invalid imported native revision."))?,
                        incoming_event: bookmark
                            .imported_selection
                            .head_event
                            .as_ref()
                            .map(library::Text::new)
                            .transpose()
                            .map_err(|_| Fault::too_large())?,
                    });
                }
                for bookmark in self.repository.snapshot_bookmarks() {
                    items.push(wire::Bookmark::Snapshot {
                        local_document: convert(&bookmark.document)?,
                        revision: convert(&bookmark.revision)?,
                    });
                }
                let (items, next_offset) = projection::page(items, page)?;
                S::Bookmarks(wire::BookmarkPage {
                    items: library::Rows::new(items).map_err(|_| Fault::too_large())?,
                    next_offset,
                    checkpoint,
                })
            }
            C::CancelUnattempted { job } => {
                self.pending(job.as_str())?;
                accepted(gate)?;
                self.repository
                    .cancel_unattempted(job.as_str())
                    .map_err(repository_error)?;
                S::CancelledUnattempted { job: job.clone() }
            }
            C::DiscardRejected { job } => {
                self.pending(job.as_str())?;
                accepted(gate)?;
                self.repository
                    .discard_rejected(job.as_str())
                    .map_err(repository_error)?;
                S::DiscardedRejected { job: job.clone() }
            }
            C::Pull {
                local_document,
                expected_local_revision,
                remote,
                mode,
            } => {
                let config = self.remote(remote)?;
                let id = document(local_document)?;
                let expected = revision(expected_local_revision)?;
                self.saved(&id, &expected)?;
                let current = self
                    .repository
                    .open_document(&id, Some(&expected))
                    .map_err(repository_error)?;
                let binding = self
                    .repository
                    .document_binding(&id, remote.as_str())
                    .cloned()
                    .ok_or_else(|| {
                        Fault::invalid("The saved document is not bound to this remote.")
                    })?;
                let local = match &current.selection {
                    Some(pin) => {
                        let native = sharing::Selection {
                            document_id: pin.document_id.clone(),
                            revision: pin.revision.as_str().into(),
                            head_event: pin.head_event.clone(),
                        };
                        if *mode == wire::PullMode::CleanFastForward
                            && binding.common_native
                                != repository::NativePin::Workspace(native.clone())
                        {
                            return Err(Fault::error(
                                wire::ErrorCode::HistoryConflict,
                                "Private native changes require an explicit checkout decision.",
                            ));
                        }
                        let options = ImportOptions {
                            selected_policy: match mode {
                                wire::PullMode::CleanFastForward => SelectedPolicy::FastForward {
                                    expected_local: native,
                                },
                                wire::PullMode::ExplicitCheckout => SelectedPolicy::Checkout {
                                    expected_local: native,
                                },
                            },
                            activate_selected: true,
                        };
                        *context = Some(self.context(&config, request.request_id)?);
                        accepted(gate)?;
                        self.repository
                            .pull_pinned(
                                &id,
                                &expected,
                                remote.as_str(),
                                &options,
                                &mut self.adapter,
                                &checkpoint,
                            )
                            .map_err(repository_error)?
                    }
                    None => {
                        if *mode != wire::PullMode::ExplicitCheckout {
                            return Err(Fault::error(
                                wire::ErrorCode::HistoryConflict,
                                "Eventless snapshots require an explicit revision-pinned checkout.",
                            ));
                        }
                        *context = Some(self.context(&config, request.request_id)?);
                        accepted(gate)?;
                        self.repository
                            .pull_snapshot(&id, remote.as_str(), &expected, &mut self.adapter)
                            .map_err(repository_error)?
                    }
                };
                S::Pulled {
                    local: convert(&local)?,
                }
            }
            C::BeginTreeClone {
                remote,
                root,
                destination,
                limits,
            } => {
                self.remote(remote)?;
                let root = match root {
                    wire::TreeRoot::LabRoot => folders::TreeRoot::LabRoot,
                    wire::TreeRoot::Folder { id, parent } => folders::TreeRoot::Folder {
                        id: id.clone(),
                        parent: parent.clone(),
                    },
                };
                let destination = destination.as_ref().map(folder).transpose()?;
                let limits = folders::TreeLimits {
                    max_nodes: usize::from(limits.max_nodes),
                    max_depth: usize::from(limits.max_depth),
                    max_pages: usize::from(limits.max_pages),
                    page_size: limits.page_size,
                };
                accepted(gate)?;
                let job = self
                    .repository
                    .begin_tree_clone(remote.as_str(), root, destination, limits)
                    .map_err(folder_error)?;
                S::TreeStarted(self.tree_status(&job)?)
            }
            C::ResumeTreeClone { job, max_steps } => {
                self.tree_status(job.as_str())?;
                accepted(gate)?;
                self.repository
                    .resume_tree_clone(
                        job.as_str(),
                        usize::from(*max_steps).min(4),
                        &mut self.adapter,
                    )
                    .map_err(folder_error)?;
                S::TreeUpdated(self.tree_status(job.as_str())?)
            }
            C::RetryTreeClone { job } => {
                self.tree_status(job.as_str())?;
                accepted(gate)?;
                self.repository
                    .retry_tree_clone(job.as_str())
                    .map_err(folder_error)?;
                S::TreeUpdated(self.tree_status(job.as_str())?)
            }
            C::TreeJobs { page } => {
                let jobs = self.repository.tree_clone_jobs().map_err(folder_error)?;
                let items = jobs
                    .iter()
                    .map(projection::tree)
                    .collect::<Result<Vec<_>>>()?;
                let (items, next_offset) = projection::page(items, page)?;
                S::TreeJobs(wire::TreePage {
                    items: library::Rows::new(items).map_err(|_| Fault::too_large())?,
                    next_offset,
                })
            }
            C::QueueTreeFolderMetadata { job, local_folder } => {
                let id = folder(local_folder)?;
                self.tree_status(job.as_str())?;
                accepted(gate)?;
                let queued = self
                    .repository
                    .queue_tree_folder_metadata(job.as_str(), &id)
                    .map_err(folder_error)?;
                S::Queued(self.pending(&queued)?)
            }
            C::QueueTreeDocumentMetadata {
                job,
                local_document,
            } => {
                let id = document(local_document)?;
                self.tree_status(job.as_str())?;
                accepted(gate)?;
                let queued = self
                    .repository
                    .queue_tree_document_metadata(job.as_str(), &id)
                    .map_err(folder_error)?;
                S::Queued(self.pending(&queued)?)
            }
        })
    }
    /// Only the native shell calls this after its own window/IPC authorization.
    /// accept must mark the operation irrevocably accepted before returning true;
    /// later cancellation cannot erase uncertain network or committed outcomes.
    pub fn dispatch_raw(
        &mut self,
        bytes: &[u8],
        accept: &mut impl FnMut() -> bool,
    ) -> std::result::Result<Vec<u8>, HostError> {
        let request = wire::decode_request(bytes).map_err(|_| HostError::Wire)?;
        let mut context = None;
        let outcome = match self.execute(&request, accept, &mut context) {
            Ok(success) => wire::Outcome::Success(success),
            Err(fault) => {
                if context.is_none() {
                    context = fault.context;
                }
                wire::Outcome::Failure(fault.failure)
            }
        };
        let response = wire::ResponseEnvelope {
            request_id: request.request_id,
            operation_id: request.operation_id.clone(),
            context,
            outcome,
        };
        response
            .validate_for(&request)
            .map_err(|_| HostError::Wire)?;
        wire::encode_response_for_frame(&response, bytes).map_err(|_| HostError::Unavailable)
    }
}
