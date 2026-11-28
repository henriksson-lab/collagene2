use collagene_core::documents::{
    history::{
        DocumentHistory,
        transition::{AnalysisOrigin, AnalysisProducer, Transition},
    },
    reference_edit::Operation,
    workspace::{RevisionEvent, Workspace},
};
use collagene_model::{Molecule, Topology};

fn source() -> DocumentHistory {
    DocumentHistory::new(
        "source",
        Molecule::dna("original", "AACCGGTT", Topology::Circular).unwrap(),
    )
    .unwrap()
}
fn archive(history: &DocumentHistory) -> Workspace {
    let mut workspace =
        Workspace::from_documents("typed edits", vec![history.molecule.clone()], 0).unwrap();
    workspace.documents.clear();
    workspace.history.clear();
    workspace.revisions.clear();
    workspace.active_document = history.id.clone();
    history.contribute(&mut workspace).unwrap();
    workspace.validate().unwrap();
    workspace
}

#[test]
fn ordered_native_actions_replay_the_exact_saved_product_after_reopen() {
    let mut history = source();
    let before = history.binding();
    history
        .commit(
            vec![
                Operation::ReplaceSequence {
                    start: 2,
                    end: 4,
                    replacement: "TTA".into(),
                },
                Operation::RotateOrigin { new_origin: 3 },
                Operation::Rename {
                    name: "edited".into(),
                },
            ],
            "an arbitrary descriptive label",
        )
        .unwrap();
    let bytes = collagene_core::documents::workspace::save_workspace(&archive(&history)).unwrap();
    let workspace = collagene_core::documents::workspace::load_workspace(&bytes).unwrap();
    let event = workspace
        .history
        .iter()
        .find(|event| event.id == history.head_event())
        .unwrap();
    let Transition::ReferenceEdit {
        source,
        actions,
        supplemental_sources,
    } = &event.transition
    else {
        panic!("native edit receipt missing")
    };
    assert_eq!(source, &before);
    assert!(supplemental_sources.is_empty());
    let original = history.revision_molecule(&source.event_id).unwrap();
    assert_eq!(
        Operation::apply_all(&original, actions).unwrap(),
        history.molecule
    );
    let mut reordered = actions.clone();
    reordered.swap(0, 1);
    assert_ne!(
        Operation::apply_all(&original, &reordered).unwrap(),
        history.molecule
    );
    let changed = RevisionEvent::new(
        event.document_id.clone(),
        event.revision.clone(),
        event.parents.clone(),
        event.operation.clone(),
        Transition::ReferenceEdit {
            source: source.clone(),
            actions: reordered,
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    assert_ne!(changed.id, event.id);
}

#[test]
fn mutation_outside_history_cannot_be_laundered_into_a_native_edit() {
    let mut history = source();
    let head = history.head_event().to_owned();
    history.molecule.name = "untracked mutation".into();
    assert!(
        history
            .commit(
                vec![Operation::Rename {
                    name: "next".into()
                }],
                "edit"
            )
            .is_err()
    );
    assert_eq!(history.head_event(), head);
    assert_eq!(history.event_count(), 1);
}

#[test]
fn computational_preview_is_explicitly_distinct_from_an_edit_recipe() {
    let source = source();
    let product = Molecule::dna("predicted", "TTTT", Topology::Linear).unwrap();
    let origin = AnalysisOrigin::preview(AnalysisProducer::AssemblyIntent, product.revision());
    let derived = source
        .derive_event(
            source.head_event(),
            "prediction",
            product,
            "ReplaceSequence",
            origin.clone(),
        )
        .unwrap();
    let mut workspace = archive(&source);
    derived.contribute(&mut workspace).unwrap();
    workspace.validate().unwrap();
    let event = workspace
        .history
        .iter()
        .find(|event| event.id == derived.head_event())
        .unwrap();
    assert_eq!(event.transition, Transition::Analysis { origin });
    assert_eq!(event.parents, [source.head_event()]);
}

#[test]
fn rehashed_edit_with_a_false_source_revision_is_rejected_by_workspace_ownership() {
    let mut history = source();
    history
        .commit(
            vec![Operation::Rename {
                name: "current".into(),
            }],
            "rename",
        )
        .unwrap();
    let mut workspace = archive(&history);
    let event = workspace
        .history
        .iter_mut()
        .find(|event| event.id == history.head_event())
        .unwrap();
    let Transition::ReferenceEdit {
        source,
        actions,
        supplemental_sources,
    } = &event.transition
    else {
        panic!("native edit receipt missing")
    };
    let mut false_source = source.clone();
    false_source.revision = history.molecule.revision();
    let forged = RevisionEvent::new(
        event.document_id.clone(),
        event.revision.clone(),
        event.parents.clone(),
        event.operation.clone(),
        Transition::ReferenceEdit {
            source: false_source,
            actions: actions.clone(),
            supplemental_sources: supplemental_sources.clone(),
        },
    )
    .unwrap();
    let forged_id = forged.id.clone();
    *event = forged;
    workspace
        .documents
        .iter_mut()
        .find(|document| document.id == history.id)
        .unwrap()
        .head_event = Some(forged_id);
    assert!(workspace.validate().is_err());
}

#[test]
fn supplemental_parent_admission_precedes_action_execution_and_deduplication() {
    let mut history = source();
    let before = history.binding();
    let parents = vec![before.event_id.clone(); 129];
    let error = history
        .commit_with_source_events(
            vec![Operation::ReplaceSequence {
                start: 9,
                end: 10,
                replacement: "A".into(),
            }],
            "oversized parent request",
            &parents,
        )
        .unwrap_err();
    assert!(matches!(error, collagene_model::Error::TooLarge));
    assert_eq!(history.binding(), before);
    assert_eq!(history.event_count(), 1);
}

#[test]
fn complete_transition_entry_is_admitted_before_event_hash_or_history_retention() {
    let history = source();
    let source = history.binding();
    let empty_action = Operation::ReplaceSequence {
        start: 0,
        end: 0,
        replacement: String::new(),
    };
    let action_envelope = serde_json::to_vec(&vec![empty_action]).unwrap().len();
    // The action array fits its own exact limit; the required source/transition
    // envelope does not fit a native archive entry. No edit execution is needed
    // to reject this unexportable event at its constructor.
    let replacement =
        "A".repeat(collagene_core::documents::reference_edit::MAX_EDIT_BYTES - action_envelope);
    let error = RevisionEvent::new(
        source.document_id.clone(),
        source.revision.clone(),
        vec![source.event_id.clone()],
        "bounded event",
        Transition::ReferenceEdit {
            source,
            actions: vec![Operation::ReplaceSequence {
                start: 0,
                end: 0,
                replacement,
            }],
            supplemental_sources: vec![],
        },
    )
    .unwrap_err();
    assert!(matches!(error, collagene_model::Error::TooLarge));
    assert_eq!(history.event_count(), 1);
}
