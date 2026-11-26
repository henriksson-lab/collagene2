//! Native document capture and library import. The caller owns the workspace
//! transaction and must recheck its epoch/open guard before applying a result.
use crate::{
    library_bindings::{self, Binding, Bindings},
    library_view::{Captured, LocalDocument, Opened},
};
use collagene_core::{
    history::DocumentHistory,
    workspace::{Workspace, WorkspaceAttachment},
};
use collagene_library_protocol as wire;
use collagene_library_sharing::{
    Selection,
    import::{ImportOptions, SelectedPolicy},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const BRANCH_KIND: &str = "lab-library-branches-v1";
type Attachments = BTreeMap<String, WorkspaceAttachment>;

pub fn bindings(attachments: &Attachments) -> Result<Bindings, String> {
    let mut values = attachments
        .values()
        .filter(|attachment| attachment.kind == library_bindings::ATTACHMENT_KIND);
    let Some(value) = values.next() else {
        return Ok(Bindings::default());
    };
    if values.next().is_some() || value.media_type != "application/json" {
        return Err("Workspace has ambiguous or unsupported lab-library link metadata.".into());
    }
    Bindings::decode(&value.bytes)
}

/// Commit only after the native import or a validated remote save receipt.
pub fn record_binding(attachments: &mut Attachments, binding: Binding) -> Result<(), String> {
    let mut records = bindings(attachments)?;
    records.record(binding)?;
    replace_attachment(
        attachments,
        library_bindings::ATTACHMENT_KIND,
        "Lab library links",
        records.encode()?,
    )
}

fn replace_attachment(
    attachments: &mut Attachments,
    kind: &str,
    name: &str,
    bytes: Vec<u8>,
) -> Result<(), String> {
    let old: Vec<_> = attachments
        .iter()
        .filter(|(_, value)| value.kind == kind)
        .map(|(id, _)| id.clone())
        .collect();
    if old.len() > 1 {
        return Err("Workspace has duplicate library bookkeeping attachments.".into());
    }
    if old.is_empty() && attachments.len() >= collagene_core::workspace::MAX_ATTACHMENTS {
        return Err(
            "Workspace attachment limit reached; library bookkeeping was not changed.".into(),
        );
    }
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if attachments.get(&hash).is_some_and(|value| {
        value.kind != kind || value.bytes != bytes || value.media_type != "application/json"
    }) {
        return Err("Library bookkeeping conflicts with another immutable attachment.".into());
    }
    let next = WorkspaceAttachment {
        kind: kind.into(),
        name: name.into(),
        media_type: "application/json".into(),
        bytes,
    };
    for id in old {
        attachments.remove(&id);
    }
    attachments.insert(hash, next);
    Ok(())
}

pub fn current(frame: &Workspace, epoch: u64) -> Result<LocalDocument, String> {
    let document = frame
        .documents
        .iter()
        .find(|document| document.id == frame.active_document)
        .ok_or("No active plasmid to share")?;
    let molecule = frame
        .head(&document.id)
        .map_err(|error| error.to_string())?;
    Ok(LocalDocument {
        id: document.id.clone(),
        name: molecule.name.clone(),
        revision: wire::Digest::new(&document.head).map_err(|error| error.to_string())?,
        head_event: document.head_event.clone(),
        workspace_epoch: epoch,
    })
}

#[derive(Clone, Default)]
pub struct CaptureRequest {
    pub choices: collagene_library_sharing::sbol_subset::PublicationChoices,
    pub approved: Option<collagene_library_sharing::sbol_subset::ReviewedPublication>,
    pub expected: Option<LocalDocument>,
}
#[derive(Clone)]
pub struct CaptureReview {
    pub local: LocalDocument,
    pub review: collagene_library_sharing::sbol_subset::PublicationReview,
}
pub enum CaptureOutcome {
    Ready(Captured),
    ReviewRequired(CaptureReview),
    ChoicesRequired {
        local: LocalDocument,
        requirements: Vec<collagene_library_sharing::sbol_subset::ChoiceRequirement>,
    },
}

/// Save the complete private snapshot to the embedded native catalog. This is
/// selected only by the trusted desktop IPC transport, never by a server URL.
/// Public publication choices/acknowledgements cannot silently change its scope.
pub async fn capture_private_local(
    frame: Workspace,
    epoch: u64,
    request: CaptureRequest,
) -> Result<CaptureOutcome, String> {
    bindings(&frame.attachments)?;
    let local = current(&frame, epoch)?;
    if request.approved.is_some() || !request.choices.is_empty() {
        return Err("Publication review does not apply to a private local catalog save. Prepare a fresh local save.".into());
    }
    if request
        .expected
        .as_ref()
        .is_some_and(|expected| expected != &local)
    {
        return Err("The local plasmid changed; prepare a fresh local save.".into());
    }
    let prepared = crate::library_capture::prepare_private(frame)?;
    let (frame, selection) = prepared.into_parts();
    let bytes = crate::app::encode_workspace(frame, u32::MAX - 302).await?;
    let payload = crate::library_capture::private_payload(&selection, bytes)?;
    Ok(CaptureOutcome::Ready(Captured { local, payload }))
}

pub async fn capture(
    frame: Workspace,
    epoch: u64,
    request: CaptureRequest,
) -> Result<CaptureOutcome, String> {
    // Malformed saved links cannot silently turn a revision push into a new import.
    bindings(&frame.attachments)?;
    let local = current(&frame, epoch)?;
    if request.approved.is_some() && request.expected.is_none() {
        return Err(
            "Reviewed publication requires its exact local document and workspace epoch pin."
                .into(),
        );
    }
    if request
        .expected
        .as_ref()
        .is_some_and(|expected| expected != &local)
    {
        return Err(
            "The selected plasmid changed; prepare a new publication review. Nothing was sent."
                .into(),
        );
    }
    let selection = Selection {
        document_id: local.id.clone(),
        revision: local.revision.as_str().into(),
        head_event: local.head_event.clone(),
    };
    let projected = if let Some(approved) = &request.approved {
        collagene_library_sharing::sbol_subset::project_document_closure_reviewed(
            &frame, &selection, approved,
        )
        .map_err(|error| format!("Publication review no longer matches this workspace: {error}"))?
    } else {
        match collagene_library_sharing::project_document_closure(&frame, &selection) {
            Ok(projected) => projected,
            Err(_) => {
                let requirements =
                    collagene_library_sharing::sbol_subset::choice_requirements(&frame, &selection)
                        .map_err(|error| {
                            format!("Cannot inspect reachable SBOL sequence choices: {error}")
                        })?;
                let missing = requirements.iter().any(|requirement| {
                    requirement.required_sequences.is_empty()
                        && request
                            .choices
                            .get(&requirement.original_hash)
                            .and_then(|choices| choices.get(&requirement.component_iri))
                            .is_none()
                });
                if missing {
                    return Ok(CaptureOutcome::ChoicesRequired {
                        local,
                        requirements,
                    });
                }
                let review = collagene_library_sharing::sbol_subset::prepare_documentary_subset(
                    &frame,
                    &selection,
                    &request.choices,
                )
                .map_err(|error| format!("Cannot prepare this plasmid's publication: {error}"))?;
                // No archive is encoded and no request is sent until the user
                // explicitly approves this exact, source-bound documentary review.
                return Ok(CaptureOutcome::ReviewRequired(CaptureReview {
                    local,
                    review,
                }));
            }
        }
    };
    let bytes = crate::app::encode_workspace(projected.workspace, u32::MAX - 301).await?;
    let payload = wire::Payload {
        kind: wire::PayloadKind::WorkspaceArchive,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).map_err(|error| error.to_string())?,
        selection: Some(wire::NativeSelection {
            document_id: wire::Text::new(&selection.document_id)
                .map_err(|error| error.to_string())?,
            revision: local.revision.clone(),
            head_event: selection
                .head_event
                .map(wire::Text::new)
                .transpose()
                .map_err(|error| error.to_string())?,
        }),
        sharing: Some(wire::WorkspaceSharing::DocumentClosure),
    };
    payload.validate().map_err(|error| error.to_string())?;
    Ok(CaptureOutcome::Ready(Captured { local, payload }))
}

pub struct Decoded {
    pub frame: Workspace,
    pub selection: Selection,
}

pub async fn decode(opened: &Opened) -> Result<Decoded, String> {
    opened
        .payload
        .validate()
        .map_err(|error| error.to_string())?;
    if !opened.guard.is_current() {
        return Err("Library session changed before import.".into());
    }
    let (frame, selection) = match opened.payload.kind {
        wire::PayloadKind::WorkspaceArchive => {
            let pin = opened
                .payload
                .selection
                .as_ref()
                .ok_or("Library archive has no exact document selection")?;
            let selection = Selection {
                document_id: pin.document_id.to_string(),
                revision: pin.revision.as_str().into(),
                head_event: pin
                    .head_event
                    .as_ref()
                    .map(|event| event.as_str().to_owned()),
            };
            let frame = crate::app::decode_workspace(
                opened.payload.bytes.as_slice().to_vec(),
                u32::MAX - 302,
            )
            .await?;
            // Even an explicitly shared whole archive contributes only this plasmid's
            // declared dependency closure to the local workspace.
            let projected = if opened.guard.is_private_local()
                && opened.payload.sharing == Some(wire::WorkspaceSharing::ExplicitWholeWorkspace)
            {
                collagene_library_sharing::project_private_document_closure(&frame, &selection)
            } else {
                // Historical native DocumentClosure records and all public
                // transports retain the original strict publication rule.
                collagene_library_sharing::project_document_closure(&frame, &selection)
            }
            .map_err(|error| format!("Cannot import this plasmid's native history: {error}"))?;
            (projected.workspace, selection)
        }
        wire::PayloadKind::NativeDocument => {
            let molecule = collagene_core::protocol::load_document(&opened.payload.bytes)
                .map_err(|error| error.to_string())?;
            let mut identity = blake3::Hasher::new();
            for part in [
                opened.remote.server.as_str(),
                opened.remote.lab.as_str(),
                opened.document.id.as_str(),
            ] {
                identity.update(&(part.len() as u64).to_le_bytes());
                identity.update(part.as_bytes());
            }
            let id = format!("library-{}", identity.finalize().to_hex());
            let history = DocumentHistory::new(&id, molecule).map_err(|error| error.to_string())?;
            let selection = Selection {
                document_id: id.clone(),
                revision: history.head_revision().into(),
                head_event: Some(history.head_event().into()),
            };
            let mut frame = Workspace {
                name: "Library plasmid".into(),
                active_document: id,
                documents: vec![],
                revisions: BTreeMap::new(),
                history: vec![],
                candidates: vec![],
                active_candidate: None,
                attachments: BTreeMap::new(),
            };
            history
                .contribute(&mut frame)
                .map_err(|error| error.to_string())?;
            frame.validate().map_err(|error| error.to_string())?;
            (frame, selection)
        }
    };
    if !opened.guard.is_current() {
        return Err("Library session changed during import; local work was not changed.".into());
    }
    Ok(Decoded { frame, selection })
}

pub fn merge(local: &Workspace, decoded: &Decoded, opened: &Opened) -> Result<Workspace, String> {
    let records = bindings(&local.attachments)?;
    let mut options = ImportOptions::default();
    if let Some(binding) = records.get(&opened.remote, &decoded.selection.document_id) {
        if binding.remote_document == opened.document.id {
            if let Some(document) = local
                .documents
                .iter()
                .find(|document| document.id == binding.local_document)
            {
                // Pull only fast-forwards a still-clean local copy. A dirty or
                // unrelated head goes through strict conflict handling instead.
                if document.head == binding.local_revision.as_str()
                    && document.head_event == binding.local_head_event
                {
                    let expected_local = Selection {
                        document_id: document.id.clone(),
                        revision: document.head.clone(),
                        head_event: document.head_event.clone(),
                    };
                    options.selected_policy = if opened.revision.id != opened.document.head.id {
                        // Choosing a historical revision is an explicit checkout of
                        // a clean copy. Keep its previous navigation as a bookmark;
                        // dirty local branches still take the strict conflict path.
                        SelectedPolicy::Checkout { expected_local }
                    } else {
                        SelectedPolicy::FastForward { expected_local }
                    };
                }
            }
        }
    }
    let merger = if opened.guard.is_private_local()
        && opened.payload.sharing == Some(wire::WorkspaceSharing::ExplicitWholeWorkspace)
    {
        collagene_library_sharing::import::merge_private_document_closure
    } else {
        collagene_library_sharing::import::merge_document_closure
    };
    let report = merger(local, &decoded.frame, &decoded.selection, &options)
        .map_err(|error| format!("Library import did not change local work: {error}"))?;
    let mut frame = report.workspace;
    if let Some(bookmark) = report.displaced_selected {
        retain_branch(&mut frame.attachments, bookmark)?;
    }
    frame.validate().map_err(|error| error.to_string())?;
    Ok(frame)
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Branches {
    version: u16,
    records: Vec<collagene_library_sharing::import::HeadBookmark>,
}

fn retain_branch(
    attachments: &mut Attachments,
    bookmark: collagene_library_sharing::import::HeadBookmark,
) -> Result<(), String> {
    let mut existing = attachments
        .values()
        .filter(|attachment| attachment.kind == BRANCH_KIND);
    let mut records = if let Some(attachment) = existing.next() {
        if existing.next().is_some()
            || attachment.media_type != "application/json"
            || attachment.bytes.len() > 4 * 1024 * 1024
        {
            return Err("Unsupported local branch bookmarks; import was not committed.".into());
        }
        let value: Branches = serde_json::from_slice(&attachment.bytes)
            .map_err(|_| "Malformed local branch bookmarks; import was not committed.")?;
        if value.version != 1 || value.records.len() > 128 {
            return Err("Unsupported local branch bookmarks; import was not committed.".into());
        }
        value
    } else {
        Branches {
            version: 1,
            records: vec![],
        }
    };
    if !records.records.iter().any(|old| {
        old.document.id == bookmark.document.id
            && old.document.head_event == bookmark.document.head_event
            && old.document.head == bookmark.document.head
    }) {
        if records.records.len() >= 128 {
            return Err("Local branch bookmark limit reached; export or review branches before importing another update.".into());
        }
        records.records.push(bookmark);
    }
    let bytes = serde_json::to_vec(&records).map_err(|error| error.to_string())?;
    if bytes.len() > 4 * 1024 * 1024 {
        return Err("Local branch bookmark size limit reached; import was not committed.".into());
    }
    replace_attachment(
        attachments,
        BRANCH_KIND,
        "Retained local library branches",
        bytes,
    )
}

#[cfg(test)]
mod documentary_review_tests {
    use super::*;
    use collagene_core::{Molecule, Topology};
    use collagene_design_exchange as sbol;
    use serde_json::json;
    use std::{
        future::Future,
        pin::pin,
        task::{Context, Poll, Waker},
    };

    fn fixture() -> Workspace {
        let raw=serde_json::to_vec(&json!({"@context":{"sbol":sbol::SBOL},"@graph":[
            {"@id":"https://example.test/c","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.test"},"sbol:type":[{"@id":sbol::DNA},{"@id":sbol::LINEAR}],"sbol:hasSequence":{"@id":"https://example.test/s"}},
            {"@id":"https://example.test/s","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.test"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"AGTCCGATGCTACGTAGGACTTCA"},
            {"@id":"https://example.test/private-activity","@type":"https://example.test/Activity","https://example.test/note":"PRIVATE-NEVER-PUBLISHED"}
        ]})).unwrap();
        let stock = sbol::import(&raw)
            .unwrap()
            .project(
                &sbol::Selection {
                    component: "https://example.test/c".into(),
                    sequence: "https://example.test/s".into(),
                    topology: Some(Topology::Linear),
                },
                sbol::ProjectionPolicy::AllowOmissions,
            )
            .unwrap()
            .molecule;
        let mut workspace = Workspace::from_documents(
            "Private local workspace",
            vec![
                stock,
                Molecule::dna(
                    "Private unrelated plasmid",
                    "CCCGGGTTTAAA",
                    Topology::Circular,
                )
                .unwrap(),
            ],
            0,
        )
        .unwrap();
        workspace.attachments.insert(
            blake3::hash(&raw).to_hex().to_string(),
            WorkspaceAttachment {
                kind: collagene_library_sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
                name: "original-design.sbol.jsonld".into(),
                media_type: "application/ld+json".into(),
                bytes: raw,
            },
        );
        workspace
    }
    #[test]
    fn unapproved_subset_capture_returns_review_before_archive_worker_or_transport() {
        let workspace = fixture();
        let original = collagene_core::workspace::save_workspace(&workspace).unwrap();
        let mut future = pin!(capture(workspace.clone(), 7, CaptureRequest::default()));
        let Poll::Ready(Ok(CaptureOutcome::ReviewRequired(review))) = future
            .as_mut()
            .poll(&mut Context::from_waker(Waker::noop()))
        else {
            panic!("Unapproved publication must return review before awaiting archive encoding");
        };
        assert_eq!(review.local, current(&workspace, 7).unwrap());
        assert_eq!(review.review.entries().len(), 1);
        let entry = &review.review.entries()[0];
        assert_eq!(entry.omissions.unrelated_objects, 1);
        assert_ne!(entry.original_hash, entry.subset_hash);
        assert_eq!(review.review.inventory().documents.len(), 1);
        assert_eq!(
            original,
            collagene_core::workspace::save_workspace(&workspace).unwrap()
        );
        // A valid backend acknowledgement cannot authorize a different UI epoch.
        let request = CaptureRequest {
            approved: Some(review.review.approve()),
            expected: Some(review.local),
            ..Default::default()
        };
        let mut stale = pin!(capture(workspace, 8, request));
        assert!(matches!(
            stale.as_mut().poll(&mut Context::from_waker(Waker::noop())),
            Poll::Ready(Err(_))
        ));
    }
}
