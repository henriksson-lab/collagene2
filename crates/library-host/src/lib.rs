//! Native validation/dispatch only. No listener, login, actor claim or membership
//! inference. The caller supplies a freshly authenticated server-side LabAccess.
#![cfg(not(target_arch = "wasm32"))]
use collagene_core::{protocol, workspace};
use collagene_lab_library as store;
use collagene_library_protocol as wire;
use serde::{Serialize, de::DeserializeOwned};
use std::{collections::BTreeSet, io::Write};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid library framing or fields")]
    Wire,
    #[error("invalid native document, archive or selected native event")]
    Native,
    #[error("trusted publication scope does not match this lab, archive and selection")]
    PublicationScope,
    #[error("host library boundary unavailable")]
    Unavailable,
}
/// Exact retained entities, not proof that they are an appropriate closure.
/// The trusted host derives the permitted closure under its publication rules;
/// it MUST NOT merely echo a client-provided inventory into `document_closure`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArchiveInventory {
    pub documents: BTreeSet<String>,
    pub revisions: BTreeSet<String>,
    pub events: BTreeSet<String>,
    pub candidates: BTreeSet<String>,
    pub attachments: BTreeSet<String>,
}
/// Read-only native inspection; bytes/hash and private scope receipt are bound.
/// This contains no scientific/material authority and implements no Deserialize.
#[derive(Debug)]
pub struct PublicationReview {
    archive_hash: wire::Digest,
    selection: Option<wire::NativeSelection>,
    sharing: Option<wire::WorkspaceSharing>,
    inventory: Option<ArchiveInventory>,
}
impl PublicationReview {
    pub fn inventory(&self) -> Option<&ArchiveInventory> {
        self.inventory.as_ref()
    }
    pub fn selection(&self) -> Option<&wire::NativeSelection> {
        self.selection.as_ref()
    }
    pub fn archive_hash(&self) -> &wire::Digest {
        &self.archive_hash
    }
}
/// Non-wire publication evidence, bound to an already authenticated lab. Scope
/// is separate from membership and from the request's untrusted sharing enum.
/// Ordinary Save can derive this automatically from server-owned closure rules.
#[derive(Debug)]
pub struct PublicationScope {
    lab: store::LabId,
    archive_hash: wire::Digest,
    selection: wire::NativeSelection,
    sharing: wire::WorkspaceSharing,
    inventory: ArchiveInventory,
}
impl PublicationScope {
    /// Bind an independently derived permitted document/dependency closure to
    /// an exact validated archive. Extra documents/history/attachments fail.
    pub fn document_closure(
        access: &store::LabAccess,
        review: &PublicationReview,
        permitted: &ArchiveInventory,
    ) -> Result<Self, Error> {
        if review.sharing != Some(wire::WorkspaceSharing::DocumentClosure)
            || review.inventory.as_ref() != Some(permitted)
        {
            return Err(Error::PublicationScope);
        }
        Self::bind(access, review, wire::WorkspaceSharing::DocumentClosure)
    }
    /// Trusted host calls only after receiving explicit whole-workspace sharing
    /// intent through its authenticated interaction, never from the enum alone.
    pub fn explicit_whole_workspace(
        access: &store::LabAccess,
        review: &PublicationReview,
    ) -> Result<Self, Error> {
        if review.sharing != Some(wire::WorkspaceSharing::ExplicitWholeWorkspace) {
            return Err(Error::PublicationScope);
        }
        Self::bind(
            access,
            review,
            wire::WorkspaceSharing::ExplicitWholeWorkspace,
        )
    }
    fn bind(
        access: &store::LabAccess,
        review: &PublicationReview,
        sharing: wire::WorkspaceSharing,
    ) -> Result<Self, Error> {
        Ok(Self {
            lab: access.lab_id().clone(),
            archive_hash: review.archive_hash.clone(),
            selection: review.selection.clone().ok_or(Error::PublicationScope)?,
            sharing,
            inventory: review.inventory.clone().ok_or(Error::PublicationScope)?,
        })
    }
    fn validate(&self, access: &store::LabAccess, review: &PublicationReview) -> Result<(), Error> {
        if &self.lab != access.lab_id()
            || self.archive_hash != review.archive_hash
            || Some(&self.selection) != review.selection.as_ref()
            || Some(self.sharing) != review.sharing
            || Some(&self.inventory) != review.inventory.as_ref()
        {
            return Err(Error::PublicationScope);
        }
        Ok(())
    }
}
/// Existing native loaders enforce native versions, checksums, expanded archive
/// bounds, molecule/history/run structure. No scientific report is trusted here.
pub fn review_payload(payload: &wire::Payload) -> Result<PublicationReview, Error> {
    payload.validate().map_err(|_| Error::Wire)?;
    let inventory = match payload.kind {
        wire::PayloadKind::NativeDocument => {
            protocol::load_document(payload.bytes.as_slice()).map_err(|_| Error::Native)?;
            None
        }
        wire::PayloadKind::WorkspaceArchive => {
            let frame =
                workspace::load_workspace(payload.bytes.as_slice()).map_err(|_| Error::Native)?;
            Some(workspace_inventory(payload, &frame)?)
        }
    };
    Ok(PublicationReview {
        archive_hash: payload.archive_hash.clone(),
        selection: payload.selection.clone(),
        sharing: payload.sharing,
        inventory,
    })
}
fn workspace_inventory(
    payload: &wire::Payload,
    frame: &workspace::Workspace,
) -> Result<ArchiveInventory, Error> {
    let selected = payload.selection.as_ref().ok_or(Error::Native)?;
    let document = frame
        .documents
        .iter()
        .find(|d| d.id == selected.document_id.as_str())
        .ok_or(Error::Native)?;
    if document.head != selected.revision.as_str()
        || document.head_event.as_deref()
            != selected.head_event.as_ref().map(wire::NativeId::as_str)
    {
        return Err(Error::Native);
    }
    if frame
        .revisions
        .get(&document.head)
        .ok_or(Error::Native)?
        .revision()
        != selected.revision.as_str()
    {
        return Err(Error::Native);
    }
    Ok(ArchiveInventory {
        documents: frame.documents.iter().map(|d| d.id.clone()).collect(),
        revisions: frame.revisions.keys().cloned().collect(),
        events: frame.history.iter().map(|e| e.id.clone()).collect(),
        candidates: frame
            .candidates
            .iter()
            .map(|c| c.plan_revision.clone())
            .collect(),
        attachments: frame.attachments.keys().cloned().collect(),
    })
}
/// Automatic ordinary-Save scope. Native documents need no archive scope;
/// workspace input must already be the exact independently derived closure.
/// Whole-workspace sharing deliberately requires a separate trusted interaction.
/// No request-supplied list or sharing enum can authorize extra local records.
pub fn document_publication_scope(
    access: &store::LabAccess,
    payload: &wire::Payload,
) -> Result<Option<PublicationScope>, Error> {
    payload.validate().map_err(|_| Error::Wire)?;
    if payload.kind == wire::PayloadKind::NativeDocument {
        review_payload(payload)?;
        return Ok(None);
    }
    if payload.sharing != Some(wire::WorkspaceSharing::DocumentClosure) {
        return Err(Error::PublicationScope);
    }
    let frame = workspace::load_workspace(payload.bytes.as_slice()).map_err(|_| Error::Native)?;
    let actual = workspace_inventory(payload, &frame)?;
    let selected = payload.selection.as_ref().ok_or(Error::Native)?;
    let selection = collagene_library_sharing::Selection {
        document_id: selected.document_id.as_str().into(),
        revision: selected.revision.as_str().into(),
        head_event: selected.head_event.as_ref().map(|id| id.as_str().into()),
    };
    let derived = collagene_library_sharing::inspect_document_closure(&frame, &selection)
        .map_err(|_| Error::PublicationScope)?;
    // IDs alone cannot bind dependency heads: a source document may now have an
    // unrelated unpublished head while a run needs an older exact source event.
    if frame.documents.len() != derived.document_heads.len()
        || frame.documents.iter().any(|doc| {
            derived
                .document_heads
                .get(&doc.id)
                .is_none_or(|head| head.revision != doc.head || head.head_event != doc.head_event)
        })
    {
        return Err(Error::PublicationScope);
    }
    let permitted = ArchiveInventory {
        documents: derived.documents,
        revisions: derived.revisions,
        events: derived.events,
        candidates: derived.candidates,
        attachments: derived.attachments,
    };
    let review = PublicationReview {
        archive_hash: payload.archive_hash.clone(),
        selection: payload.selection.clone(),
        sharing: payload.sharing,
        inventory: Some(actual),
    };
    PublicationScope::document_closure(access, &review, &permitted).map(Some)
}

fn publication(
    access: &store::LabAccess,
    payload: &wire::Payload,
    scope: Option<&PublicationScope>,
) -> Result<(), Error> {
    let review = review_payload(payload)?;
    if review.inventory.is_some() {
        scope
            .ok_or(Error::PublicationScope)?
            .validate(access, &review)?;
    }
    Ok(())
}

// Structural metadata only. Binary archive bytes take the direct ownership-aware
// path below, never a JSON conversion. Independent DTOs are checked again by the
// wire validator before any response is emitted.
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > 4 << 20)
        {
            return Err(std::io::Error::other("bounded library metadata"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn convert<T: Serialize, U: DeserializeOwned>(value: &T) -> Result<U, Error> {
    let mut bytes = Bounded(Vec::new());
    serde_json::to_writer(&mut bytes, value).map_err(|_| Error::Unavailable)?;
    serde_json::from_slice(&bytes.0).map_err(|_| Error::Unavailable)
}
fn folder(id: &wire::FolderId) -> store::FolderId {
    store::FolderId::parse(id.as_str()).expect("validated matching identifier contract")
}
fn document(id: &wire::DocumentId) -> store::DocumentId {
    store::DocumentId::parse(id.as_str()).expect("validated matching identifier contract")
}
fn revision(id: &wire::RevisionId) -> store::RevisionId {
    store::RevisionId::parse(id.as_str()).expect("validated matching revision contract")
}
fn payload_to_store(payload: &wire::Payload) -> Result<store::Payload, Error> {
    Ok(store::Payload {
        kind: match payload.kind {
            wire::PayloadKind::NativeDocument => store::PayloadKind::NativeDocument,
            wire::PayloadKind::WorkspaceArchive => store::PayloadKind::WorkspaceArchive,
        },
        bytes: payload.bytes.as_slice().to_vec(),
        archive_hash: store::Digest::parse(payload.archive_hash.as_str())
            .map_err(|_| Error::Wire)?,
        selection: payload.selection.as_ref().map(convert).transpose()?,
        sharing: payload.sharing.as_ref().map(convert).transpose()?,
    })
}
fn payload_to_wire(payload: store::Payload) -> Result<wire::Payload, Error> {
    let result = wire::Payload {
        kind: match payload.kind {
            store::PayloadKind::NativeDocument => wire::PayloadKind::NativeDocument,
            store::PayloadKind::WorkspaceArchive => wire::PayloadKind::WorkspaceArchive,
        },
        bytes: wire::Blob::new(payload.bytes).map_err(|_| Error::Unavailable)?,
        archive_hash: wire::Digest::new(payload.archive_hash.as_str())
            .map_err(|_| Error::Unavailable)?,
        selection: payload.selection.as_ref().map(convert).transpose()?,
        sharing: payload.sharing.as_ref().map(convert).transpose()?,
    };
    result.validate().map_err(|_| Error::Unavailable)?;
    Ok(result)
}
fn failure(
    code: wire::ErrorCode,
    message: &str,
    conflict: Option<wire::Conflict>,
) -> wire::Failure {
    wire::Failure {
        code,
        message: wire::Text::new(message).expect("static bounded public message"),
        conflict,
    }
}
fn boundary_failure(error: Error) -> wire::Failure {
    match error {
        Error::Wire | Error::Native => failure(
            wire::ErrorCode::Invalid,
            "Invalid native library payload or selected document/event.",
            None,
        ),
        Error::PublicationScope => failure(
            wire::ErrorCode::Forbidden,
            "Workspace publication scope is missing or does not match this lab and exact archive.",
            None,
        ),
        Error::Unavailable => failure(
            wire::ErrorCode::Unavailable,
            "The library service could not complete this operation.",
            None,
        ),
    }
}
fn storage_failure(
    db: &store::Library,
    access: &store::LabAccess,
    request: &wire::Request,
    error: store::Error,
) -> wire::Failure {
    use store::Error as E;
    use wire::{Conflict, ErrorCode as C, Request as R};
    match error {
        E::NotFound => failure(C::NotFound, "Item not found in this lab.", None),
        E::AlreadyExists => failure(
            C::AlreadyExists,
            "That name already exists in this folder.",
            None,
        ),
        E::NotEmpty => failure(C::NotEmpty, "Only an empty folder can be deleted.", None),
        E::Invalid(_) => failure(C::Invalid, "Invalid library operation.", None),
        E::Limit(_) => failure(
            C::TooLarge,
            "A library storage or size limit was reached.",
            None,
        ),
        E::CatalogConflict { actual, .. } => failure(
            C::Conflict,
            "The catalog changed; restart this listing.",
            Some(Conflict::Catalog {
                actual_generation: actual,
            }),
        ),
        E::Conflict { .. } | E::MetadataConflict { .. } => {
            let conflict = match request {
                R::RenameFolder { id, .. }
                | R::MoveFolder { id, .. }
                | R::DeleteEmptyFolder { id, .. } => db
                    .folder(access, &folder(id))
                    .ok()
                    .and_then(|current| convert(&current).ok())
                    .map(|current| Conflict::Folder { current }),
                R::SaveRevision { id, .. }
                | R::RenameDocument { id, .. }
                | R::MoveDocument { id, .. } => db
                    .document(access, &document(id))
                    .ok()
                    .and_then(|current| convert(&current).ok())
                    .map(|current| Conflict::Document { current }),
                _ => None,
            };
            if conflict.is_some() {
                failure(
                    C::Conflict,
                    "The item changed. Keep your edits and review the current revision.",
                    conflict,
                )
            } else {
                failure(
                    C::Unavailable,
                    "The current conflict state is unavailable; keep your edits.",
                    None,
                )
            }
        }
        E::Corrupt(_) | E::UnsupportedDatabase | E::Database(_) => failure(
            C::Unavailable,
            "The library service could not complete this operation.",
            None,
        ),
    }
}

enum DispatchError {
    Storage(store::Error),
    Boundary(Error),
}
impl From<store::Error> for DispatchError {
    fn from(e: store::Error) -> Self {
        Self::Storage(e)
    }
}
impl From<Error> for DispatchError {
    fn from(e: Error) -> Self {
        Self::Boundary(e)
    }
}
fn execute(
    db: &mut store::Library,
    access: &store::LabAccess,
    request: &wire::Request,
    scope: Option<&PublicationScope>,
) -> Result<wire::Success, DispatchError> {
    use wire::{Request as R, Success as S};
    Ok(match request {
        R::Browse { parent, page } => S::Browsed(convert(&db.browse(
            access,
            parent.as_ref().map(folder).as_ref(),
            page.offset,
            page.limit,
            page.expected_catalog_generation,
        )?)?),
        R::CreateFolder { parent, name } => S::CreatedFolder(convert(&db.create_folder(
            access,
            parent.as_ref().map(folder).as_ref(),
            name.as_str(),
        )?)?),
        R::RenameFolder {
            id,
            expected_generation,
            name,
        } => S::RenamedFolder(convert(&db.rename_folder(
            access,
            &folder(id),
            *expected_generation,
            name.as_str(),
        )?)?),
        R::MoveFolder {
            id,
            expected_generation,
            parent,
        } => S::MovedFolder(convert(&db.move_folder(
            access,
            &folder(id),
            *expected_generation,
            parent.as_ref().map(folder).as_ref(),
        )?)?),
        R::DeleteEmptyFolder {
            id,
            expected_generation,
        } => {
            db.delete_empty_folder(access, &folder(id), *expected_generation)?;
            S::DeletedEmptyFolder(id.clone())
        }
        R::ImportDocument {
            folder: destination,
            name,
            payload,
        } => {
            publication(access, payload, scope)?;
            S::Imported(convert(&db.create_document(
                access,
                destination.as_ref().map(folder).as_ref(),
                name.as_str(),
                &payload_to_store(payload)?,
            )?)?)
        }
        R::SaveRevision {
            id,
            base_revision,
            payload,
        } => {
            publication(access, payload, scope)?;
            let current = db.save(
                access,
                &document(id),
                &revision(base_revision),
                &payload_to_store(payload)?,
            )?;
            let revision = convert(&current.head)?;
            S::Saved {
                document: convert(&current)?,
                revision,
            }
        }
        R::OpenDocument {
            id,
            revision: selected,
        } => {
            let current = db.document(access, &document(id))?;
            let pin = selected
                .as_ref()
                .map(revision)
                .unwrap_or_else(|| current.head.id.clone());
            let metadata = db.revision(access, &document(id), &pin)?;
            let payload = payload_to_wire(db.payload(access, &document(id), &pin)?)?;
            // Stored bytes never bypass current loader/pin checks. Reads need no
            // new publication scope: sharing was checked at the immutable write.
            review_payload(&payload)?;
            S::Opened {
                document: convert(&current)?,
                revision: convert(&metadata)?,
                payload,
            }
        }
        R::RenameDocument {
            id,
            expected_generation,
            name,
        } => S::RenamedDocument(convert(&db.rename_document(
            access,
            &document(id),
            *expected_generation,
            name.as_str(),
        )?)?),
        R::MoveDocument {
            id,
            expected_generation,
            folder: destination,
        } => S::MovedDocument(convert(&db.move_document(
            access,
            &document(id),
            *expected_generation,
            destination.as_ref().map(folder).as_ref(),
        )?)?),
        R::Search { query, scope, page } => {
            let native = match scope {
                wire::SearchScope::AllFolders => store::SearchScope::AllFolders,
                wire::SearchScope::Folder {
                    folder: parent,
                    recursive,
                } => store::SearchScope::Folder {
                    folder: parent.as_ref().map(folder),
                    recursive: *recursive,
                },
            };
            S::Searched(convert(&db.search(
                access,
                query.as_str(),
                &native,
                page.offset,
                page.limit,
                page.expected_catalog_generation,
            )?)?)
        }
        R::ListRevisions { id, page } => S::Revisions(convert(&db.revision_page(
            access,
            &document(id),
            page.offset,
            page.limit,
            page.expected_catalog_generation,
        )?)?),
    })
}
/// Dispatch an authenticated native request. Caller MUST recheck membership on
/// each invocation. There is deliberately no actor/lab assertion parameter.
pub fn dispatch(
    db: &mut store::Library,
    access: &store::LabAccess,
    request: &wire::RequestEnvelope,
    scope: Option<&PublicationScope>,
) -> Result<wire::ResponseEnvelope, Error> {
    request.validate().map_err(|_| Error::Wire)?;
    let lab = db.lab(access).map_err(|_| Error::Unavailable)?;
    let context = wire::LabContext {
        lab_id: wire::LabId::new(lab.id.as_str()).map_err(|_| Error::Unavailable)?,
        name: wire::Text::new(lab.name).map_err(|_| Error::Unavailable)?,
    };
    let outcome = match execute(db, access, &request.request, scope) {
        Ok(success) => wire::Outcome::Success(success),
        Err(DispatchError::Boundary(e)) => wire::Outcome::Failure(boundary_failure(e)),
        Err(DispatchError::Storage(e)) => {
            wire::Outcome::Failure(storage_failure(db, access, &request.request, e))
        }
    };
    let response = wire::ResponseEnvelope {
        request_id: request.request_id,
        context: Some(context),
        outcome,
    };
    response
        .validate_for(request)
        .map_err(|_| Error::Unavailable)?;
    Ok(response)
}
/// Decode only bounded checked frames; malformed requests return a transport
/// error because their correlation ID cannot be trusted or fabricated.
pub fn dispatch_frame(
    db: &mut store::Library,
    access: &store::LabAccess,
    bytes: &[u8],
    scope: Option<&PublicationScope>,
) -> Result<Vec<u8>, Error> {
    let request = wire::decode_request(bytes).map_err(|_| Error::Wire)?;
    let response = dispatch(db, access, &request, scope)?;
    wire::encode_response(&response).map_err(|_| Error::Unavailable)
}
