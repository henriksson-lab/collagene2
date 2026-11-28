#[path = "../../core/tests/common/protocol_layout.rs"]
mod fixture;
use collagene_core::{documents::workspace, output::layout};
use collagene_library_sharing::{self as sharing, Selection};
#[test]
fn ordinary_plasmid_omits_plate_but_explicit_plate_scope_preserves_exact_dependencies() {
    let (w, r) = fixture::fixture();
    let captured = layout::capture(&w, &r).unwrap();
    let before = workspace::save_workspace(&captured.workspace).unwrap();
    let d = &captured.workspace.documents[0];
    let selection = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let ordinary = sharing::project_document_closure(&captured.workspace, &selection).unwrap();
    assert!(
        !ordinary
            .inventory
            .attachments
            .contains(&captured.attachment_hash)
    );
    assert!(ordinary.inventory.omitted.attachments >= 2);
    let projected =
        sharing::project_protocol_layout_closure(&captured.workspace, &captured.attachment_hash)
            .unwrap();
    assert!(
        projected
            .review
            .inventory
            .attachments
            .contains(&captured.attachment_hash)
    );
    assert_eq!(projected.workspace.documents.len(), 1);
    assert!(
        !projected
            .workspace
            .revisions
            .values()
            .any(|m| m.name == "Private unrelated")
    );
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&projected.workspace).unwrap())
            .unwrap();
    let saved = layout::load(&reopened, &captured.attachment_hash).unwrap();
    assert_eq!(
        serde_json::to_value(saved.request).unwrap(),
        serde_json::to_value(r).unwrap()
    );
    assert_eq!(
        workspace::save_workspace(&captured.workspace).unwrap(),
        before
    );
}
#[test]
fn missing_or_corrupt_plate_dependencies_refuse_without_guessing_other_records() {
    let (w, r) = fixture::fixture();
    let captured = layout::capture(&w, &r).unwrap();
    let mut missing = captured.workspace;
    let batch = match &r.assignments[0].source {
        layout::Source::Planned { batch_hash, .. } => batch_hash,
        _ => unreachable!(),
    };
    missing.attachments.remove(batch);
    assert!(sharing::project_protocol_layout_closure(&missing, &captured.attachment_hash).is_err());
    assert!(sharing::inspect_protocol_layout_closure(&w, &"0".repeat(64)).is_err());
}
#[test]
fn released_layout_retains_actual_run_and_replays_after_projection() {
    let (w, r) = fixture::released_fixture();
    let captured = layout::capture(&w, &r).unwrap();
    let run_id = match &r.assignments[0].source {
        layout::Source::Released { run_id, .. } => run_id,
        _ => unreachable!(),
    };
    let projected =
        sharing::project_protocol_layout_closure(&captured.workspace, &captured.attachment_hash)
            .unwrap();
    layout::load(&projected.workspace, &captured.attachment_hash).unwrap();
    let replay =
        collagene_core::material::lineage::replay_run(&projected.workspace, run_id).unwrap();
    assert!(replay.products_match);
    assert!(replay.engine_matches);
    assert_eq!(
        replay.current_report.status,
        collagene_core::verification::evidence::Status::ConsistentUnderAssumptions
    );
    assert!(
        !projected
            .workspace
            .revisions
            .values()
            .any(|m| m.name == "Private unrelated")
    );
}
