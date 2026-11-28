use super::model::{self, Candidate, Capture, Draft, PcrDraft};
use collagene_core::{
    chemistry::{methylation::MethylationState, tracked_duplex::checked_selection},
    cloning::screening as native,
    documents::workspace::Workspace,
};
use collagene_model::{Molecule, Topology};

fn fixture() -> (Workspace, crate::sources::identity::SourceInfo, Draft) {
    let workspace = Workspace::from_documents(
        "screening",
        vec![
            Molecule::dna("intended source", "AAAAGAATTCTTTT", Topology::Linear).unwrap(),
            Molecule::dna("parent source", "CCCCCCCCCCCCCCCCCCCC", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let source = model::info(&workspace, &workspace.history[0].id).unwrap();
    let parent = workspace
        .history
        .iter()
        .find(|e| e.document_id != source.binding.document_id)
        .and_then(|e| model::info(&workspace, &e.id))
        .unwrap();
    let mut draft = Draft::default();
    draft.add_component(&source).unwrap();
    draft.add_component(&parent).unwrap();
    draft.add_decoy().unwrap();
    draft.subjects[1].label = "Intended plus parent carryover".into();
    draft.subjects[1].role = native::Role::Decoy(native::DecoyKind::ParentCarryover);
    draft.subjects[1].components = draft.components.iter().map(|c| c.id.clone()).collect();
    for component in &mut draft.components {
        component.declare_generation = true;
        component.generation.dam = MethylationState::Unmodified;
        component.generation.dcm = MethylationState::Unmodified;
        component.generation.cpg = MethylationState::Unmodified;
    }
    draft.candidates.push(Candidate::Digest {
        label: "EcoRI".into(),
        enzymes: vec![checked_selection("EcoRI").unwrap().enzymes[0].id.clone()],
    });
    (workspace, source, draft)
}

#[test]
fn explicit_mixture_and_canonical_original_inputs_survive_capture_without_inferred_decoys() {
    let (workspace, source, draft) = fixture();
    let capture = Capture::prepare(&workspace, &source, 7, draft.clone()).unwrap();
    assert_eq!(capture.request.components.len(), 2);
    assert_eq!(
        capture.request.subjects[1].components,
        draft.subjects[1].components
    );
    assert_eq!(
        capture.request.subjects[1].role,
        native::Role::Decoy(native::DecoyKind::ParentCarryover)
    );
    for component in &capture.request.components {
        assert_eq!(component.revision, component.molecule.revision());
        assert!(matches!(
            component.origin,
            collagene_core::material::derivation::SourceOrigin::OriginalReference
        ));
    }
    let report = native::assess(&capture.request).unwrap();
    capture.validate_report(&report).unwrap();
    assert_eq!(report.tests[0].signatures[1].sizes_twice_bp, [14, 40]);
    let export: serde_json::Value =
        serde_json::from_slice(&capture.export_json(&report).unwrap()).unwrap();
    assert_eq!(export["format"], "collagene-screening-ui-evidence-1");
    assert_eq!(
        export["component_events"][0]["event"]["event_id"],
        source.binding.event_id
    );
}

#[test]
fn workspace_event_epoch_and_every_draft_setting_guard_response_export_and_navigation() {
    let (workspace, source, draft) = fixture();
    let capture = Capture::prepare(&workspace, &source, 7, draft.clone()).unwrap();
    capture.guard(&workspace, &source, 7, &draft).unwrap();
    let mut changed = workspace.clone();
    changed.name = "renamed workspace".into();
    assert!(capture.guard(&changed, &source, 7, &draft).is_err());
    assert!(capture.guard(&workspace, &source, 8, &draft).is_err());
    let mut other_event = source.clone();
    other_event.binding.event_id = workspace
        .history
        .iter()
        .find(|e| e.id != source.binding.event_id)
        .unwrap()
        .id
        .clone();
    assert!(capture.guard(&workspace, &other_event, 7, &draft).is_err());
    let mut changed = draft.clone();
    changed.tolerance = "2".into();
    assert!(capture.guard(&workspace, &source, 7, &changed).is_err());
    let mut changed = draft.clone();
    changed.components[1].generation.dam = MethylationState::Unknown;
    assert!(capture.guard(&workspace, &source, 7, &changed).is_err());
    let mut changed = draft.clone();
    changed.subjects[1].components.pop();
    assert!(capture.guard(&workspace, &source, 7, &changed).is_err());
    let mut changed = draft.clone();
    changed.candidates.push(Candidate::Pcr(PcrDraft {
        name: "PCR".into(),
        forward: "ACGT".into(),
        forward_tail: "GG".into(),
        reverse: "CGTA".into(),
        reverse_tail: String::new(),
        bottom: false,
    }));
    assert!(capture.guard(&workspace, &source, 7, &changed).is_err());
}

#[test]
fn only_the_exact_active_original_reference_gets_native_recognition_navigation() {
    let (workspace, source, draft) = fixture();
    let capture = Capture::prepare(&workspace, &source, 7, draft).unwrap();
    let report = native::assess(&capture.request).unwrap();
    assert_eq!(
        capture.recognition_ranges(&report, 0, 0, &source).unwrap(),
        [(4, 10)]
    );
    let other = model::info(&workspace, &capture.draft.components[1].binding.event_id).unwrap();
    assert!(capture.recognition_ranges(&report, 0, 0, &other).is_err());
    assert!(
        capture
            .recognition_ranges(&report, usize::MAX, 0, &source)
            .is_err()
    );
    let mut changed = report.clone();
    changed.cells[0].cleavage.as_mut().unwrap().sites[0].recognition_start = source.length + 1;
    assert!(capture.recognition_ranges(&changed, 0, 0, &source).is_err());
}

#[test]
fn deletion_prunes_only_explicit_component_links_and_never_manufactures_replacement_sequence() {
    let (_, _, mut draft) = fixture();
    let removed = draft.components[1].id.clone();
    let intended_id = draft.components[0].id.clone();
    draft.remove_component(&removed);
    assert_eq!(draft.components.len(), 1);
    assert_eq!(draft.subjects[1].components, [intended_id]);
    assert_eq!(
        draft.subjects[1].role,
        native::Role::Decoy(native::DecoyKind::ParentCarryover)
    );
}

#[test]
fn event_menu_is_bounded_and_search_uses_exact_identity_even_with_repeated_names() {
    let molecules = (0..110)
        .map(|_| Molecule::dna("same name", "ACGT", Topology::Linear).unwrap())
        .collect();
    let workspace = Workspace::from_documents("many explicit events", molecules, 0).unwrap();
    assert_eq!(model::event_choices(&workspace, "", "").len(), 100);
    let event = &workspace.history[0];
    let filtered = model::event_choices(&workspace, &event.id, "");
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].binding.event_id, event.id);
    let selected = model::event_choices(&workspace, "does not match any name", &event.id);
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].binding.event_id, event.id);
}

#[test]
fn invalid_or_unassociated_sources_and_numeric_inputs_fail_before_dispatch() {
    let (workspace, source, mut draft) = fixture();
    draft.components[0].binding.revision = "0".repeat(64);
    assert!(Capture::prepare(&workspace, &source, 0, draft).is_err());
    for value in ["", "-1", "1.5", " 1", "999999999999999999999999999999"] {
        assert!(model::integer(value, "limit").is_err());
    }
    assert_eq!(model::integer("2", "tolerance").unwrap(), 2);
}
