//! Documentary intended targets retain actual native source ancestry without widening it.
#[path = "../../core/tests/common/actual_strand.rs"]
mod fixture;
use collagene_core::{
    assembly_intent,
    history::DocumentHistory,
    run_material,
    workspace::{self, RevisionEvent, Workspace, WorkspaceDocument},
};
use collagene_library_sharing::{self as sharing, Selection};
fn adopt(w: &mut Workspace, target: collagene_core::Molecule, parents: Vec<String>) -> Selection {
    let id = "independent-actual-strand-design".to_string();
    let revision = target.revision();
    let event = RevisionEvent::new(
        id.clone(),
        revision.clone(),
        parents,
        "Independently declared actual source strand",
    )
    .unwrap();
    w.revisions.insert(revision.clone(), target);
    w.documents.push(WorkspaceDocument {
        id: id.clone(),
        head: revision.clone(),
        undo: vec![],
        redo: vec![],
        head_event: Some(event.id.clone()),
        undo_events: vec![],
        redo_events: vec![],
    });
    w.history.push(event.clone());
    w.active_document = id.clone();
    w.validate().unwrap();
    Selection {
        document_id: id,
        revision,
        head_event: Some(event.id),
    }
}
#[test]
fn actual_material_target_requires_existing_exact_parents_and_replays_their_owned_run() {
    let (mut w, selected, input, r) = fixture::fixture();
    let unrelated = DocumentHistory::new(
        "private-unrelated",
        collagene_core::Molecule::dna(
            "Unrelated private plasmid",
            "CCCCAAAATTTTGGGG",
            collagene_core::Topology::Circular,
        )
        .unwrap(),
    )
    .unwrap();
    let private = unrelated.id.clone();
    unrelated.contribute(&mut w).unwrap();
    let report = assembly_intent::build(&[input], &r, Some(&w)).unwrap();
    assert_eq!(report.target.sequence.as_str(), "AATTCGTAAG");
    let mut missing = w.clone();
    let orphan = adopt(&mut missing, report.target.clone(), vec![]);
    assert!(sharing::project_document_closure(&missing, &orphan).is_err());
    let selection = adopt(
        &mut w,
        report.target,
        vec![selected.event().event_id.clone()],
    );
    let before = workspace::save_workspace(&w).unwrap();
    for projection in [
        sharing::project_document_closure(&w, &selection).unwrap(),
        sharing::project_private_document_closure(&w, &selection).unwrap(),
    ] {
        assert!(
            !projection
                .workspace
                .documents
                .iter()
                .any(|d| d.id == private)
        );
        let bytes = workspace::save_workspace(&projection.workspace).unwrap();
        let loaded = workspace::load_workspace(&bytes).unwrap();
        let origins = run_material::documentary_origins(&loaded, &[selected.clone()], &[])
            .unwrap()
            .materials;
        origins[0]
            .replay(&loaded.revisions[&selected.event().revision])
            .unwrap();
        assert_eq!(
            loaded
                .head(&selection.document_id)
                .unwrap()
                .sequence
                .as_str(),
            "AATTCGTAAG"
        );
        assert!(
            loaded
                .history
                .iter()
                .any(|e| e.id == selected.event().event_id)
        );
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn forged_unkeyed_descriptor_cannot_publish_private_event_or_claim_another_run() {
    let (mut w, selected, input, r) = fixture::fixture();
    let private = DocumentHistory::new(
        "private-foreign-event",
        collagene_core::Molecule::dna(
            "Private sequence",
            "CCCCAAAATTTTGGGG",
            collagene_core::Topology::Circular,
        )
        .unwrap(),
    )
    .unwrap();
    let foreign = collagene_core::run_lineage::EventBinding {
        document_id: private.id.clone(),
        event_id: private.head_event().into(),
        revision: private.head_revision().into(),
    };
    private.contribute(&mut w).unwrap();
    let report = assembly_intent::build(&[input], &r, Some(&w)).unwrap();
    let original = workspace::save_workspace(&w).unwrap();
    for use_private in [true, false] {
        let mut forged = report.target.clone();
        let mut descriptors = assembly_intent::material_source_descriptors(&forged).unwrap();
        let descriptor = &mut descriptors[0];
        let collagene_core::run_material::MaterialSelection::Run(choice) =
            &mut descriptor.selection
        else {
            panic!("actual run fixture");
        };
        if use_private {
            descriptor.event = foreign.clone();
            choice.event = foreign.clone();
        } else {
            choice.run_id = "0".repeat(64);
        }
        let json = serde_json::to_string(&descriptors).unwrap();
        forged
            .metadata
            .insert(assembly_intent::MATERIAL_METADATA.into(), json.clone());
        forged.metadata.insert(
            "assembly_material_sources:revision".into(),
            blake3::hash(json.as_bytes()).to_hex().to_string(),
        );
        // Shape/hash alone accepts this documentary claim, intentionally.
        assert!(assembly_intent::material_source_descriptors(&forged).is_ok());
        let mut copied = w.clone();
        let selection = adopt(&mut copied, forged, vec![selected.event().event_id.clone()]);
        assert!(sharing::project_document_closure(&copied, &selection).is_err());
        assert!(sharing::project_private_document_closure(&copied, &selection).is_err());
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), original);
}
