//! Separate private native persistence from explicit public document publication.
//! This module performs no transport, hostname detection or scientific inference.
use collagene_core::workspace::Workspace;
use collagene_library_protocol as wire;
use collagene_library_sharing::Selection;

pub struct PrivateCapture {
    workspace: Workspace,
    selection: Selection,
}
impl PrivateCapture {
    pub fn into_parts(self) -> (Workspace, Selection) {
        (self.workspace, self.selection)
    }
}
/// The complete private workspace stays local, including original exchange
/// evidence, unrelated local branches and native candidate/navigation state.
pub fn prepare_private(workspace: Workspace) -> Result<PrivateCapture, String> {
    workspace.validate().map_err(|error| error.to_string())?;
    let document = workspace
        .documents
        .iter()
        .find(|document| document.id == workspace.active_document)
        .ok_or("No active local plasmid to save")?;
    let selection = Selection {
        document_id: document.id.clone(),
        revision: document.head.clone(),
        head_event: document.head_event.clone(),
    };
    Ok(PrivateCapture {
        workspace,
        selection,
    })
}
/// The native archive worker owns encoding. Its full bounded bytes and selected
/// pin are subsequently validated by the native repository before any commit.
pub fn private_payload(selection: &Selection, bytes: Vec<u8>) -> Result<wire::Payload, String> {
    let payload = wire::Payload {
        kind: wire::PayloadKind::WorkspaceArchive,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).map_err(|error| error.to_string())?,
        selection: Some(wire::NativeSelection {
            document_id: wire::Text::new(&selection.document_id)
                .map_err(|error| error.to_string())?,
            revision: wire::Digest::new(&selection.revision).map_err(|error| error.to_string())?,
            head_event: selection
                .head_event
                .as_ref()
                .map(wire::Text::new)
                .transpose()
                .map_err(|error| error.to_string())?,
        }),
        sharing: Some(wire::WorkspaceSharing::ExplicitWholeWorkspace),
    };
    payload.validate().map_err(|error| error.to_string())?;
    Ok(payload)
}
/// This existing browser pane publishes selected document closures only. Reject
/// private local snapshots before serializing or constructing any HTTP request.
/// A future whole-workspace publication UI needs a separate explicit boundary.
pub fn ensure_public_document_request(request: &wire::Request) -> Result<(), String> {
    let payload = match request {
        wire::Request::ImportDocument { payload, .. }
        | wire::Request::SaveRevision { payload, .. } => payload,
        _ => return Ok(()),
    };
    if payload.kind == wire::PayloadKind::WorkspaceArchive
        && payload.sharing != Some(wire::WorkspaceSharing::DocumentClosure)
    {
        return Err("Private local workspace snapshots cannot be sent by the shared plasmid library. Prepare an explicit selected-document publication instead.".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::workspace::{self, WorkspaceAttachment};
    use collagene_model::{Molecule, Topology};
    pub(super) fn fixture() -> Workspace {
        let raw=br#"{"@context":{"sbol":"http://sbols.org/v3#"},"@graph":[{"@id":"https://private.invalid/original","@type":"https://private.invalid/Activity","https://private.invalid/note":"LOCAL-ORIGINAL-MUST-SURVIVE"}]}"#.to_vec();
        let mut molecule = Molecule::dna(
            "Offline plasmid",
            "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
            Topology::Linear,
        )
        .unwrap();
        let hash = blake3::hash(&raw).to_hex().to_string();
        molecule
            .metadata
            .insert("sbol:source_hash".into(), hash.clone());
        molecule.metadata.insert(
            "sbol:component_iri".into(),
            "https://private.invalid/c".into(),
        );
        molecule.metadata.insert(
            "sbol:sequence_iri".into(),
            "https://private.invalid/s".into(),
        );
        molecule.metadata.insert(
            "sbol:exchange_profile".into(),
            "collagene-sbol-3.1-structural-jsonld-1".into(),
        );
        let mut workspace = Workspace::from_documents(
            "Private local work",
            vec![
                molecule,
                Molecule::dna("Unrelated local branch", "CCCGGGTTTAAA", Topology::Circular)
                    .unwrap(),
            ],
            0,
        )
        .unwrap();
        workspace.attachments.insert(
            hash,
            WorkspaceAttachment {
                kind: collagene_library_sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
                name: "Exact private original.jsonld".into(),
                media_type: "application/ld+json".into(),
                bytes: raw,
            },
        );
        workspace
    }
    #[test]
    fn private_capture_preserves_exact_full_workspace_and_selected_native_pin() {
        let w = fixture();
        let before = serde_json::to_value(&w).unwrap();
        assert!(
            collagene_library_sharing::project_document_closure(
                &w,
                &Selection {
                    document_id: w.active_document.clone(),
                    revision: w.documents[0].head.clone(),
                    head_event: w.documents[0].head_event.clone()
                }
            )
            .is_err()
        );
        let (frame, pin) = prepare_private(w).unwrap().into_parts();
        let payload = private_payload(&pin, workspace::save_workspace(&frame).unwrap()).unwrap();
        assert_eq!(
            payload.sharing,
            Some(wire::WorkspaceSharing::ExplicitWholeWorkspace)
        );
        let restored = workspace::load_workspace(payload.bytes.as_slice()).unwrap();
        assert_eq!(serde_json::to_value(&restored).unwrap(), before);
        let selected = payload.selection.as_ref().unwrap();
        assert_eq!(selected.document_id.as_str(), pin.document_id);
        assert_eq!(selected.revision.as_str(), pin.revision);
        assert_eq!(
            selected.head_event.as_ref().map(|s| s.as_str()),
            pin.head_event.as_deref()
        );
    }
    #[test]
    fn public_request_gate_refuses_private_import_and_save_before_transport() {
        let (w, pin) = prepare_private(fixture()).unwrap().into_parts();
        let payload = private_payload(&pin, workspace::save_workspace(&w).unwrap()).unwrap();
        let import = wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("Local plasmid").unwrap(),
            payload: payload.clone(),
        };
        let save = wire::Request::SaveRevision {
            id: wire::DocumentId::new(format!("doc-{:032x}", 1)).unwrap(),
            base_revision: wire::RevisionId::new("11".repeat(32)).unwrap(),
            payload: payload.clone(),
        };
        assert!(ensure_public_document_request(&import).is_err());
        assert!(ensure_public_document_request(&save).is_err());
        let mut public = payload;
        public.sharing = Some(wire::WorkspaceSharing::DocumentClosure);
        assert!(
            ensure_public_document_request(&wire::Request::ImportDocument {
                folder: None,
                name: wire::Text::new("Selected publication").unwrap(),
                payload: public
            })
            .is_ok()
        );
        assert!(
            ensure_public_document_request(&wire::Request::Browse {
                parent: None,
                page: wire::PageRequest::default()
            })
            .is_ok()
        );
    }
}
