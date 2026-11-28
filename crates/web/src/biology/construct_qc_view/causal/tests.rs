use super::{
    model::Capture,
    selection::{self, Form},
};
use collagene_core::{
    documents::{history::DocumentHistory, reference_edit::Operation, workspace::Workspace},
    material::lineage::EventBinding,
};
use collagene_model::{Feature, Location, Molecule, Topology};

fn archive(history: &DocumentHistory) -> Workspace {
    let mut workspace =
        Workspace::from_documents("causal capture", vec![history.molecule.clone()], 0).unwrap();
    workspace.documents.clear();
    workspace.history.clear();
    workspace.revisions.clear();
    workspace.active_document = history.id.clone();
    history.contribute(&mut workspace).unwrap();
    workspace.validate().unwrap();
    workspace
}
fn fixture() -> (DocumentHistory, Workspace, EventBinding, Form) {
    let mut molecule = Molecule::dna("source", "ATGAAACCCTAA", Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 12),
        qualifiers: vec![],
        needs_review: false,
    });
    let mut history = DocumentHistory::new("coding", molecule).unwrap();
    let source = history.binding();
    history
        .commit(
            vec![Operation::Rename {
                name: "product".into(),
            }],
            "rename",
        )
        .unwrap();
    let form = Form {
        source_event: source.event_id,
        source_feature: "0".into(),
        product_feature: "0".into(),
        rationale: "Explicitly compare the retained source and product coding roles.".into(),
        genetic_code: "1".into(),
    };
    let workspace = archive(&history);
    let product = history.binding();
    (history, workspace, product, form)
}
#[test]
fn same_content_new_event_invalidates_captured_workspace_even_with_historical_product_pin() {
    let (mut history, workspace, product, form) = fixture();
    let capture = Capture::prepare(
        &workspace,
        &selection::events(&workspace),
        product.clone(),
        9,
        form.clone(),
    )
    .unwrap();
    capture.guard(&workspace, &product, 9, &form).unwrap();
    let original = history.molecule.revision();
    history
        .commit(
            vec![Operation::Rename {
                name: "temporary".into(),
            }],
            "rename",
        )
        .unwrap();
    history
        .commit(
            vec![Operation::Rename {
                name: "product".into(),
            }],
            "rename back",
        )
        .unwrap();
    assert_eq!(history.molecule.revision(), original);
    assert_ne!(history.binding().event_id, product.event_id);
    let updated = archive(&history);
    assert!(capture.guard(&updated, &product, 9, &form).is_err());
}
#[test]
fn changed_cds_declaration_or_code_cannot_reuse_a_worker_capture() {
    let (_, workspace, product, form) = fixture();
    let capture = Capture::prepare(
        &workspace,
        &selection::events(&workspace),
        product.clone(),
        1,
        form.clone(),
    )
    .unwrap();
    let mut changed = form.clone();
    changed.genetic_code = "2".into();
    assert!(capture.guard(&workspace, &product, 1, &changed).is_err());
    let mut changed = form.clone();
    changed.rationale.push_str(" Changed pairing premise.");
    assert!(capture.guard(&workspace, &product, 1, &changed).is_err());
    assert!(capture.guard(&workspace, &product, 2, &form).is_err());
}
#[test]
fn stale_or_undeclared_annotation_choices_never_become_native_requests() {
    let (_, workspace, product, mut form) = fixture();
    let choices = selection::events(&workspace);
    form.rationale.clear();
    assert!(form.request(&workspace, &choices, &product).is_err());
    form.rationale = "Explicit role declaration".into();
    form.source_feature = "5".into();
    assert!(form.request(&workspace, &choices, &product).is_err());
    form.source_feature = "0".into();
    form.source_event = "a descriptive label".into();
    assert!(form.request(&workspace, &choices, &product).is_err());
}

#[test]
fn source_and_cds_menus_keep_selected_identity_with_at_most_one_hundred_rows() {
    let choices = (0..150)
        .map(|index| selection::EventChoice {
            binding: EventBinding {
                document_id: format!("document-{index}"),
                event_id: format!("event-{index:064x}"),
                revision: "a".repeat(64),
            },
            label: format!("Matching source {index}"),
        })
        .collect::<Vec<_>>();
    let selected = &choices[149].binding.event_id;
    let rows = selection::event_menu(&choices, "Matching", selected);
    assert_eq!(rows.len(), 100);
    assert_eq!(rows[0].binding.event_id, *selected);
    assert_eq!(
        rows.iter()
            .filter(|row| row.binding.event_id == *selected)
            .count(),
        1
    );
    let rows = selection::event_menu(&choices, "document-128", "");
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].binding.document_id, "document-128");
    let rows = selection::event_menu(&choices, &choices[128].binding.event_id, selected);
    assert_eq!(rows.len(), 2, "full exact event IDs remain searchable");
    let cds = (0..150)
        .map(|index| selection::CdsChoice {
            index,
            label: format!("{} · matching CDS", index + 1),
        })
        .collect::<Vec<_>>();
    let rows = selection::cds_menu(&cds, "matching", "149");
    assert_eq!(rows.len(), 100);
    assert_eq!(rows[0].index, 149);
    assert_eq!(selection::cds_menu(&cds, "150 ·", "")[0].index, 149);
}
#[test]
fn navigation_refuses_empty_or_outside_product_coordinates_without_truncation() {
    assert!(selection::navigation_ranges_valid(&[(0, 3), (9, 12)], 12));
    for ranges in [
        vec![],
        vec![(3, 3)],
        vec![(4, 3)],
        vec![(0, 13)],
        vec![(0, 3), (11, 13)],
    ] {
        assert!(!selection::navigation_ranges_valid(&ranges, 12));
    }
}
#[test]
fn accepted_report_uses_native_request_binding_and_cannot_ignore_changed_limits() {
    let (_, workspace, product, form) = fixture();
    let capture =
        Capture::prepare(&workspace, &selection::events(&workspace), product, 1, form).unwrap();
    let report =
        collagene_core::construct_qc::causal::assess(&workspace, &capture.request).unwrap();
    capture.validate_report(&workspace, &report).unwrap();
    let mut changed = capture.clone();
    changed.request.limits.max_work = 1;
    assert!(changed.validate_report(&workspace, &report).is_err());
    changed = capture;
    changed.request.pairs[0]
        .source_cds
        .rationale
        .push_str(" Different premise.");
    assert!(changed.validate_report(&workspace, &report).is_err());
}
