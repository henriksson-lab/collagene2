//! Explicit documentary publication of one immutable saved native catalog head.
//! A review is local process state, never deserializable scientific authority.
use super::*;
use sharing::sbol_subset::{PublicationChoices, PublicationReview, ReviewedPublication};

#[derive(Debug, Clone)]
struct Pins {
    document: store::Document,
    payload_hash: store::Digest,
    native: NativePin,
    remote: RemoteConfig,
    checkpoint: Option<store::Digest>,
    binding: wire::Digest,
    destination: Option<wire::FolderId>,
}
/// Source-bound review. Only readonly counts/selected references are exposed;
/// omitted private SBOL identifiers and bytes never become UI response fields.
#[derive(Debug, Clone)]
pub struct PushReview {
    pins: Pins,
    publication: PublicationReview,
}
#[derive(Debug, Clone)]
pub struct ApprovedPushReview {
    pins: Pins,
    publication: ReviewedPublication,
}
impl PushReview {
    pub fn local_document(&self) -> &store::DocumentId {
        &self.pins.document.id
    }
    pub fn local_revision(&self) -> &store::RevisionId {
        &self.pins.document.head.id
    }
    pub fn remote(&self) -> &RemoteConfig {
        &self.pins.remote
    }
    pub fn checkpoint(&self) -> &Option<store::Digest> {
        &self.pins.checkpoint
    }
    pub fn native(&self) -> &NativePin {
        &self.pins.native
    }
    pub fn review(&self) -> &PublicationReview {
        &self.publication
    }
    /// The host calls this only for an explicit acknowledgement of the retained
    /// review object; a client-supplied hash cannot construct this token.
    pub fn approve(&self) -> ApprovedPushReview {
        ApprovedPushReview {
            pins: self.pins.clone(),
            publication: self.publication.approve(),
        }
    }
}
fn binding_hash(binding: Option<&DocumentBinding>) -> Result<wire::Digest> {
    Ok(wire::Digest::of(
        &serde_json::to_vec(&binding).map_err(|_| Error::Invalid("binding encoding"))?,
    ))
}
fn same_document(a: &store::Document, b: &store::Document) -> Result<bool> {
    Ok(serde_json::to_vec(a).map_err(|_| Error::Native)?
        == serde_json::to_vec(b).map_err(|_| Error::Native)?)
}
fn destination(
    repo: &Repository,
    doc: &store::Document,
    remote: &str,
) -> Result<Option<wire::FolderId>> {
    // A mounted/bound document's local placement is not an implicit remote
    // move. Only first publication needs an explicitly bound destination.
    if repo.document_binding(&doc.id, remote).is_some() {
        return Ok(None);
    }
    doc.folder
        .as_ref()
        .map(|folder| {
            repo.folder_binding(folder, remote)
                .map(|binding| binding.remote.id.clone())
                .ok_or(Error::UnboundFolder)
        })
        .transpose()
}
fn selected_workspace(payload: &store::Payload) -> Result<(Workspace, NativePin, Selection)> {
    let native = selected(payload)?;
    let NativePin::Workspace(selection) = &native else {
        return Err(Error::Invalid(
            "documentary review requires a saved native workspace archive",
        ));
    };
    let workspace = workspace::load_workspace(&payload.bytes).map_err(|_| Error::Native)?;
    Ok((workspace, native.clone(), selection.clone()))
}
impl Repository {
    /// Discover only structurally reachable alternatives in the selected saved
    /// archive. No first sequence is selected and no public bytes are queued.
    pub fn push_review_choices(
        &mut self,
        id: &store::DocumentId,
        expected_local: &store::RevisionId,
        remote_id: &str,
    ) -> Result<Vec<sharing::sbol_subset::ChoiceRequirement>> {
        self.reload()?;
        if !self.state.remotes.contains_key(remote_id) {
            return Err(Error::Invalid("unknown remote"));
        }
        let document = self.db.document(&self.access, id)?;
        if &document.head.id != expected_local {
            return Err(Error::Stale);
        }
        let payload = self.db.payload(&self.access, id, expected_local)?;
        let (workspace, _, selection) = selected_workspace(&payload)?;
        sharing::sbol_subset::choice_requirements(&workspace, &selection).map_err(publication_error)
    }
    /// Reload canonical state before deriving a bounded review from the exact
    /// saved local head. No transport, outbox or local archive is changed.
    pub fn prepare_push_review(
        &mut self,
        id: &store::DocumentId,
        expected_local: &store::RevisionId,
        remote_id: &str,
        choices: &PublicationChoices,
    ) -> Result<PushReview> {
        self.reload()?;
        let remote = self
            .state
            .remotes
            .get(remote_id)
            .cloned()
            .ok_or(Error::Invalid("unknown remote"))?;
        let document = self.db.document(&self.access, id)?;
        if &document.head.id != expected_local {
            return Err(Error::Stale);
        }
        let payload = self.db.payload(&self.access, id, expected_local)?;
        let (workspace, native, selection) = selected_workspace(&payload)?;
        let publication =
            sharing::sbol_subset::prepare_documentary_subset(&workspace, &selection, choices)
                .map_err(publication_error)?;
        let pins = Pins {
            payload_hash: payload.archive_hash,
            native,
            checkpoint: self.checkpoint.clone(),
            binding: binding_hash(self.document_binding(id, remote_id))?,
            destination: destination(self, &document, remote_id)?,
            document,
            remote,
        };
        Ok(PushReview { pins, publication })
    }
    /// Re-derive the approved subset from the still-identical private revision,
    /// then put exact validated public bytes in the ordinary reliable outbox.
    /// Once queued, retries/restart use those immutable bytes, not this token.
    pub fn queue_reviewed_push(&mut self, approved: &ApprovedPushReview) -> Result<String> {
        self.reload()?;
        let pins = &approved.pins;
        if self.checkpoint != pins.checkpoint
            || self.state.remotes.get(&pins.remote.id) != Some(&pins.remote)
        {
            return Err(Error::Stale);
        }
        let document = self.db.document(&self.access, &pins.document.id)?;
        if !same_document(&document, &pins.document)?
            || binding_hash(self.document_binding(&document.id, &pins.remote.id))? != pins.binding
            || destination(self, &document, &pins.remote.id)? != pins.destination
        {
            return Err(Error::Stale);
        }
        let payload = self
            .db
            .payload(&self.access, &document.id, &document.head.id)?;
        if payload.archive_hash != pins.payload_hash {
            return Err(Error::Stale);
        }
        let (workspace, native, selection) = selected_workspace(&payload)?;
        if native != pins.native {
            return Err(Error::Stale);
        }
        let projected = sharing::sbol_subset::project_document_closure_reviewed(
            &workspace,
            &selection,
            &approved.publication,
        )
        .map_err(publication_error)?;
        // Independently validate the result's ordinary supported public closure;
        // the acknowledgement never bypasses native identity/closure checks.
        sharing::inspect_document_closure(&projected.workspace, &selection)
            .map_err(publication_error)?;
        let bytes = workspace::save_workspace(&projected.workspace).map_err(|_| Error::Native)?;
        let payload = wire::Payload {
            kind: wire::PayloadKind::WorkspaceArchive,
            archive_hash: wire::Digest::of(&bytes),
            bytes: wire::Blob::new(bytes).map_err(|_| Error::Native)?,
            selection: Some(wire::NativeSelection {
                document_id: wire::NativeId::new(&selection.document_id)
                    .map_err(|_| Error::Native)?,
                revision: wire::Digest::new(&selection.revision).map_err(|_| Error::Native)?,
                head_event: selection
                    .head_event
                    .as_ref()
                    .map(wire::NativeId::new)
                    .transpose()
                    .map_err(|_| Error::Native)?,
            }),
            sharing: Some(wire::WorkspaceSharing::DocumentClosure),
        };
        self.queue_prepared_push(document, pins.remote.clone(), payload, native, true)
    }
}
