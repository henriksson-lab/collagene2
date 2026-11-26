//! Local library v1 bridge. Writes go through Repository's atomic mutation and
//! private checkpoint path; the second native connection is read-only in use.
use crate::WindowAccess;
use crate::{DesktopHost, HostError, HostResult, convert};
use collagene_desktop_protocol::OperationId;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sync as sync;
fn folder(id: &wire::FolderId) -> HostResult<store::FolderId> {
    store::FolderId::parse(id.as_str()).map_err(|_| HostError::Wire)
}
fn document(id: &wire::DocumentId) -> HostResult<store::DocumentId> {
    store::DocumentId::parse(id.as_str()).map_err(|_| HostError::Wire)
}
fn revision(id: &wire::RevisionId) -> HostResult<store::RevisionId> {
    store::RevisionId::parse(id.as_str()).map_err(|_| HostError::Wire)
}
fn payload(payload: store::Payload) -> HostResult<wire::Payload> {
    Ok(wire::Payload {
        kind: match payload.kind {
            store::PayloadKind::NativeDocument => wire::PayloadKind::NativeDocument,
            store::PayloadKind::WorkspaceArchive => wire::PayloadKind::WorkspaceArchive,
        },
        bytes: wire::Blob::new(payload.bytes).map_err(|_| HostError::Unavailable)?,
        archive_hash: wire::Digest::new(payload.archive_hash.as_str())
            .map_err(|_| HostError::Unavailable)?,
        selection: payload.selection.as_ref().map(convert).transpose()?,
        sharing: payload.sharing.as_ref().map(convert).transpose()?,
    })
}
fn operation(request: &wire::Request) -> Option<sync::MutationOperation> {
    use sync::MutationOperation as O;
    use wire::Request as R;
    Some(match request {
        R::CreateFolder { parent, name } => O::CreateFolder {
            parent: parent.clone(),
            name: name.clone(),
        },
        R::RenameFolder {
            id,
            expected_generation,
            name,
        } => O::RenameFolder {
            id: id.clone(),
            expected_generation: *expected_generation,
            name: name.clone(),
        },
        R::MoveFolder {
            id,
            expected_generation,
            parent,
        } => O::MoveFolder {
            id: id.clone(),
            expected_generation: *expected_generation,
            parent: parent.clone(),
        },
        R::DeleteEmptyFolder {
            id,
            expected_generation,
        } => O::DeleteEmptyFolder {
            id: id.clone(),
            expected_generation: *expected_generation,
        },
        R::ImportDocument {
            folder,
            name,
            payload,
        } => O::ImportDocument {
            folder: folder.clone(),
            name: name.clone(),
            payload: payload.clone(),
        },
        R::SaveRevision {
            id,
            base_revision,
            payload,
        } => O::SaveRevision {
            id: id.clone(),
            base_revision: base_revision.clone(),
            payload: payload.clone(),
        },
        R::RenameDocument {
            id,
            expected_generation,
            name,
        } => O::RenameDocument {
            id: id.clone(),
            expected_generation: *expected_generation,
            name: name.clone(),
        },
        R::MoveDocument {
            id,
            expected_generation,
            folder,
        } => O::MoveDocument {
            id: id.clone(),
            expected_generation: *expected_generation,
            folder: folder.clone(),
        },
        _ => return None,
    })
}
fn rejection(rejected: durable::Rejection) -> HostResult<wire::Failure> {
    let code = match rejected.code {
        durable::RejectionCode::Invalid => wire::ErrorCode::Invalid,
        durable::RejectionCode::NotFound => wire::ErrorCode::NotFound,
        durable::RejectionCode::AlreadyExists => wire::ErrorCode::AlreadyExists,
        durable::RejectionCode::NotEmpty => wire::ErrorCode::NotEmpty,
        durable::RejectionCode::TooLarge => wire::ErrorCode::TooLarge,
        durable::RejectionCode::Conflict => wire::ErrorCode::Conflict,
    };
    let conflict = rejected.conflict.as_ref().map(convert).transpose()?;
    // Domain messages are native bounded public descriptions; do not propagate
    // a database, filesystem or parser error into the raw IPC channel.
    Ok(wire::Failure {
        code,
        message: wire::Text::new(rejected.message).map_err(|_| HostError::Unavailable)?,
        conflict,
    })
}
impl DesktopHost {
    pub(crate) fn local_library(
        &mut self,
        access: &WindowAccess,
        op: &OperationId,
        request: &wire::Request,
        accepted: &std::sync::atomic::AtomicBool,
    ) -> HostResult<wire::Success> {
        if let Some(operation) = operation(request) {
            self.repository
                .reload()
                .map_err(|_| HostError::Unavailable)?;
            let id = sync::MutationId::new(op.as_str().replacen("op-", "mut-", 1))
                .map_err(|_| HostError::Wire)?;
            let mutation = sync::native::to_storage_mutation(&sync::Mutation { id, operation })
                .map_err(|_| HostError::Wire)?;
            let control = self.control.clone();
            let outcome = self
                .repository
                .mutate_local_reliable_with_gate(mutation, || {
                    crate::commit_gate(&control, access, op, accepted).is_ok()
                })
                .map_err(|error| match error {
                    collagene_repository::Error::Cancelled => HostError::Cancelled,
                    collagene_repository::Error::Transaction(
                        durable::MutationError::IntentMismatch,
                    ) => HostError::Invalid(
                        "Operation ID was already used for another local mutation.",
                    ),
                    collagene_repository::Error::Native => {
                        HostError::Invalid("Invalid private native document/workspace.")
                    }
                    _ => HostError::Unavailable,
                })?;
            // Either a newly committed receipt or the exact durable replay has
            // established its outcome, even if the window is revoked afterward.
            accepted.store(true, std::sync::atomic::Ordering::Release);
            return match outcome {
                durable::MutationOutcome::Applied(applied) => convert(&applied),
                durable::MutationOutcome::Rejected(rejected) => {
                    Err(HostError::Library(rejection(rejected)?))
                }
            };
        }
        use wire::{Request as R, Success as S};
        let result = match request {
            R::Browse { parent, page } => self
                .reads
                .browse(
                    &self.local_access,
                    parent.as_ref().map(folder).transpose()?.as_ref(),
                    page.offset,
                    page.limit,
                    page.expected_catalog_generation,
                )
                .and_then(|page| {
                    convert(&page)
                        .map(S::Browsed)
                        .map_err(|_| store::Error::Corrupt("desktop metadata"))
                }),
            R::Search { query, scope, page } => {
                let scope = match scope {
                    wire::SearchScope::AllFolders => store::SearchScope::AllFolders,
                    wire::SearchScope::Folder {
                        folder: parent,
                        recursive,
                    } => store::SearchScope::Folder {
                        folder: parent.as_ref().map(folder).transpose()?,
                        recursive: *recursive,
                    },
                };
                self.reads
                    .search(
                        &self.local_access,
                        query.as_str(),
                        &scope,
                        page.offset,
                        page.limit,
                        page.expected_catalog_generation,
                    )
                    .and_then(|page| {
                        convert(&page)
                            .map(S::Searched)
                            .map_err(|_| store::Error::Corrupt("desktop metadata"))
                    })
            }
            R::ListRevisions { id, page } => self
                .reads
                .revision_page(
                    &self.local_access,
                    &document(id)?,
                    page.offset,
                    page.limit,
                    page.expected_catalog_generation,
                )
                .and_then(|page| {
                    convert(&page)
                        .map(S::Revisions)
                        .map_err(|_| store::Error::Corrupt("desktop metadata"))
                }),
            R::OpenDocument {
                id,
                revision: selected,
            } => {
                let local = document(id)?;
                let result = (|| -> HostResult<S> {
                    let current = self.repository.document(&local).map_err(|e| match e {
                        collagene_repository::Error::Storage(store::Error::NotFound) => {
                            HostError::Library(public_error(
                                wire::ErrorCode::NotFound,
                                "Item not found in this local library.",
                                None,
                            ))
                        }
                        _ => HostError::Unavailable,
                    })?;
                    let selected = selected
                        .as_ref()
                        .map(revision)
                        .transpose()?
                        .unwrap_or_else(|| current.head.id.clone());
                    let metadata = self
                        .reads
                        .revision(&self.local_access, &local, &selected)
                        .map_err(storage_error)?;
                    let body = self
                        .repository
                        .open_document(&local, Some(&selected))
                        .map_err(|_| HostError::Unavailable)?;
                    Ok(S::Opened {
                        document: convert(&current)?,
                        revision: convert(&metadata)?,
                        payload: payload(body)?,
                    })
                })();
                return result;
            }
            _ => return Err(HostError::Wire),
        };
        result.map_err(storage_error)
    }
}
fn public_error(
    code: wire::ErrorCode,
    message: &str,
    conflict: Option<wire::Conflict>,
) -> wire::Failure {
    wire::Failure {
        code,
        message: wire::Text::new(message).expect("bounded static desktop message"),
        conflict,
    }
}
fn storage_error(error: store::Error) -> HostError {
    HostError::Library(match error {
        store::Error::NotFound => public_error(
            wire::ErrorCode::NotFound,
            "Item not found in this local library.",
            None,
        ),
        store::Error::CatalogConflict { actual, .. } => public_error(
            wire::ErrorCode::Conflict,
            "The local catalog changed; restart this listing.",
            Some(wire::Conflict::Catalog {
                actual_generation: actual,
            }),
        ),
        store::Error::Invalid(_) => public_error(
            wire::ErrorCode::Invalid,
            "Invalid local library operation.",
            None,
        ),
        store::Error::Limit(_) => public_error(
            wire::ErrorCode::TooLarge,
            "The local library limit was reached.",
            None,
        ),
        _ => public_error(
            wire::ErrorCode::Unavailable,
            "The local library operation could not complete.",
            None,
        ),
    })
}
