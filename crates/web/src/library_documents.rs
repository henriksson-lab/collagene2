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

pub async fn capture(frame: Workspace, epoch: u64) -> Result<Captured, String> {
    // Malformed saved links cannot silently turn a revision push into a new import.
    bindings(&frame.attachments)?;
    let local = current(&frame, epoch)?;
    let selection = Selection {
        document_id: local.id.clone(),
        revision: local.revision.as_str().into(),
        head_event: local.head_event.clone(),
    };
    let projected = collagene_library_sharing::project_document_closure(&frame, &selection)
        .map_err(|error| format!("Cannot share this plasmid's history: {error}"))?;
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
    Ok(Captured { local, payload })
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
            let projected = collagene_library_sharing::project_document_closure(&frame, &selection)
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
    let report = collagene_library_sharing::import::merge_document_closure(
        local,
        &decoded.frame,
        &decoded.selection,
        &options,
    )
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
