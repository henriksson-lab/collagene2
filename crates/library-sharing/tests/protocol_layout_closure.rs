#[path = "../../core/tests/common/protocol_layout.rs"]
mod fixture;
#[path = "../../core/tests/common/paired_library.rs"]
mod paired;
use collagene_core::{documents::workspace, output::layout};
use collagene_library_sharing::{self as sharing, Selection};

#[test]
fn paired_design_plate_publication_keeps_selected_and_unselected_bottom_events() {
    let (frame, request, unselected_bottom) = paired::fixture();
    let captured = layout::capture(&frame, &request).unwrap();
    let dependencies = layout::source_bindings(&captured.workspace, &captured.layout).unwrap();
    assert!(
        dependencies
            .events
            .iter()
            .any(|event| event.event_id == unselected_bottom)
    );
    let projected =
        sharing::project_protocol_layout_closure(&captured.workspace, &captured.attachment_hash)
            .unwrap();
    assert_eq!(projected.workspace.documents.len(), 3);
    assert!(
        projected
            .workspace
            .history
            .iter()
            .any(|event| event.id == unselected_bottom)
    );
    assert!(
        !projected
            .workspace
            .revisions
            .values()
            .any(|m| m.name == "Private unrelated")
    );
    let restarted =
        workspace::load_workspace(&workspace::save_workspace(&projected.workspace).unwrap())
            .unwrap();
    let retained = layout::load(&restarted, &captured.attachment_hash).unwrap();
    assert_eq!(
        serde_json::to_value(retained.request).unwrap(),
        serde_json::to_value(&request).unwrap()
    );
    for format in [layout::ExportFormat::Html, layout::ExportFormat::Pdf] {
        let exported = layout::export(&restarted, &captured.attachment_hash, format).unwrap();
        assert!(!exported.bytes.is_empty());
    }
    let batch_hash = match &request.assignments[0].source {
        layout::Source::Planned { batch_hash, .. } => batch_hash,
        _ => unreachable!(),
    };
    let record =
        collagene_core::planning::combinatorial::batch::load(&restarted, batch_hash).unwrap();
    assert_eq!(record.rows[0].input_bindings.len(), 2);
    let mut missing = restarted;
    missing
        .history
        .retain(|event| event.id != unselected_bottom);
    assert!(collagene_core::planning::combinatorial::batch::load(&missing, batch_hash).is_err());
}
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
