#[path = "../../core/tests/common/domestication_oe.rs"]
mod fixture;
use collagene_core::{
    Molecule, Topology, domestication_preparation as preparation,
    run_lineage::{self, EventBinding, InputEventBinding},
    verification::Status,
    workspace::{self, Workspace, WorkspaceAttachment},
};
use collagene_library_sharing::{self as sharing, Selection};
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn current_designer_review_is_omitted_while_actual_native_run_and_sources_replay() {
    let (input, request) = fixture::fixture(Topology::Linear);
    let prepared = preparation::materialize(&input, &request).unwrap();
    assert!(prepared.proposal_ready, "{:?}", prepared.proposal_findings);
    let source = Molecule::dna(
        "Unrelated private construct",
        "TTTTGGGGAAAA",
        Topology::Circular,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Native run plus private designer review",
        vec![input.molecule.clone(), prepared.target.clone(), source],
        1,
    )
    .unwrap();
    let source_event = event(&w, 0);
    let target_event = event(&w, 1);
    let captured = run_lineage::record_run(
        &w,
        &prepared.plan,
        target_event.clone(),
        vec![InputEventBinding {
            input_id: input.id,
            event: source_event.clone(),
        }],
        false,
    )
    .unwrap();
    let product = captured.receipt.products.last().unwrap();
    let selection = Selection {
        document_id: product.document_id.clone(),
        revision: product.revision.clone(),
        head_event: Some(product.event_id.clone()),
    };
    let run_id = captured.receipt.run_id;
    let mut w = captured.workspace;
    let bytes=serde_json::to_vec(&serde_json::json!({"version":preparation::REVIEW_EVIDENCE_VERSION,"source":source_event,"target_event":target_event,"check":prepared.check,"prepared":prepared})).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: preparation::REVIEW_EVIDENCE_KIND.into(),
            name: "Full current domestication design review".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    w.validate().unwrap();
    let before = workspace::save_workspace(&w).unwrap();
    let projected = sharing::project_document_closure(&w, &selection).unwrap();
    assert!(!projected.inventory.attachments.contains(&hash));
    assert!(projected.inventory.omitted.attachments >= 1);
    assert!(
        projected
            .workspace
            .attachments
            .values()
            .any(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
    );
    assert!(
        !projected
            .workspace
            .revisions
            .values()
            .any(|m| m.name == "Unrelated private construct")
    );
    assert!(
        projected
            .workspace
            .history
            .iter()
            .any(|e| e.id == source_event.event_id)
    );
    assert!(
        projected
            .workspace
            .history
            .iter()
            .any(|e| e.id == target_event.event_id)
    );
    let replay = run_lineage::replay_run(&projected.workspace, &run_id).unwrap();
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(replay.products_match);
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    for unsupported in ["domestication-intent-v1", "domestication-intent-v3"] {
        let mut bad = w.clone();
        bad.attachments.get_mut(&hash).unwrap().kind = unsupported.into();
        assert!(sharing::project_document_closure(&bad, &selection).is_err());
    }
}
