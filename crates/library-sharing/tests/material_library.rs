#[path = "../../core/tests/common/material_library.rs"]
mod material_library;
use collagene_core::{documents::workspace, material::lineage as run_lineage, output::layout};
use collagene_library_sharing::{self as sharing, Selection};
#[test]
fn released_actual_bottom_closure_retains_native_origin_and_omits_unrelated_private_catalog() {
    let (commit, pin, run) = material_library::released();
    let selection = Selection {
        document_id: pin.document_id.clone(),
        revision: pin.revision.clone(),
        head_event: Some(pin.event_id.clone()),
    };
    let projected = sharing::project_document_closure(&commit.workspace, &selection).unwrap();
    assert!(
        !projected
            .workspace
            .documents
            .iter()
            .any(|d| projected.workspace.revisions[&d.head].name == "Unrelated private head")
    );
    assert!(
        !projected
            .workspace
            .attachments
            .values()
            .any(|a| a.kind == collagene_core::planning::combinatorial::batch::KIND)
    );
    let restarted =
        workspace::load_workspace(&workspace::save_workspace(&projected.workspace).unwrap())
            .unwrap();
    let replay = run_lineage::replay_run(&restarted, &run).unwrap();
    assert!(replay.products_match);
    assert_eq!(
        replay
            .current_report
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        material_library::EXPECTED
    );
}
#[test]
fn explicit_layout_closes_indivisible_material_binding_and_native_source_dependencies() {
    let (saved, declaration) = material_library::captured();
    let layout = layout::capture(
        &saved.workspace,
        &layout::Request {
            name: "Explicit native library plate".into(),
            plate: layout::Plate::Wells96,
            assignments: vec![layout::Assignment {
                id: "bottom-A1".into(),
                well: layout::Well { row: 0, column: 0 },
                label: "Native bottom".into(),
                source: layout::Source::Planned {
                    batch_hash: saved.attachment_hash.clone(),
                    construct_id: "native-circle".into(),
                    declaration_revision: declaration,
                },
            }],
        },
    )
    .unwrap();
    let projected =
        sharing::project_protocol_layout_closure(&layout.workspace, &layout.attachment_hash)
            .unwrap();
    let record = collagene_core::planning::combinatorial::batch::load(
        &projected.workspace,
        &saved.attachment_hash,
    )
    .unwrap();
    record
        .validate_source_closure(&projected.workspace)
        .unwrap();
    let replay = collagene_core::planning::combinatorial::batch::replay(
        &projected.workspace,
        &saved.attachment_hash,
        0,
    )
    .unwrap();
    assert_eq!(
        replay.current.status,
        Some(collagene_core::verification::evidence::Status::ConsistentUnderAssumptions)
    );
    assert!(
        !projected
            .workspace
            .documents
            .iter()
            .any(|d| projected.workspace.revisions[&d.head].name == "Unrelated private head")
    );
}
