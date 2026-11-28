use collagene_core::documents::history::transition::{AnalysisOrigin, AnalysisProducer};
use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::documents::{
    history::DocumentHistory,
    workspace::{self, Workspace},
};
use collagene_core::sequence::reference_transfer::{self, Part, RequestCopy, Selection};
use collagene_model::{Molecule, Topology};
fn source(id: &str) -> DocumentHistory {
    DocumentHistory::new(
        id,
        Molecule::dna("identical retained read", "ACGT", Topology::Linear).unwrap(),
    )
    .unwrap()
}
#[test]
fn multi_read_product_keeps_every_exact_event_and_source_navigation() {
    let mut a = source("a");
    let b = source("b");
    let old = a.head_event().to_string();
    let edit_end = a.molecule.sequence.len();
    a.commit(
        vec![
            Edit::ReplaceSequence {
                start: 0,
                end: edit_end,
                replacement: "AGGT".into(),
            },
            Edit::Rename {
                name: "newer read state".into(),
            },
            Edit::SetTopology {
                topology: Topology::Linear,
            },
        ],
        "edit",
    )
    .unwrap();
    let before = (
        a.head_event().to_string(),
        a.navigation_event_ids(),
        b.head_event().to_string(),
        b.navigation_event_ids(),
    );
    let consensus = Molecule::dna("inferred consensus", "ACGTACGT", Topology::Linear).unwrap();
    let origin = AnalysisOrigin::preview(AnalysisProducer::ReadContig, consensus.revision());
    let product = DocumentHistory::derive_from_events(
        &[(&a, &old), (&b, b.head_event()), (&a, &old)],
        "consensus",
        consensus,
        "Read-supported computational assembly; not physical synthesis",
        origin,
    )
    .unwrap();
    assert_eq!(
        product.summaries()[0].parents,
        vec![old.clone(), b.head_event().to_string()]
    );
    assert_eq!(product.event_count(), 1);
    assert!(!product.can_undo());
    assert_eq!(
        before,
        (
            a.head_event().to_string(),
            a.navigation_event_ids(),
            b.head_event().to_string(),
            b.navigation_event_ids()
        )
    );
    let mut w = Workspace::from_documents("test", vec![a.molecule.clone()], 0).unwrap();
    w.documents.clear();
    w.revisions.clear();
    w.history.clear();
    a.contribute(&mut w).unwrap();
    b.contribute(&mut w).unwrap();
    product.contribute(&mut w).unwrap();
    w.active_document = product.id.clone();
    let loaded = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    let restored = DocumentHistory::restore_all(&loaded).unwrap();
    let p = restored.iter().find(|d| d.id == "consensus").unwrap();
    assert_eq!(
        p.summaries()[0].parents,
        vec![old, b.head_event().to_string()]
    );
}
#[test]
fn bad_or_excessive_source_bindings_reject_without_mutating_any_source() {
    let a = source("a");
    let b = source("b");
    let m = Molecule::dna("product", "ACGTACGT", Topology::Linear).unwrap();
    let before = (a.head_event().to_string(), b.head_event().to_string());
    let origin = AnalysisOrigin::preview(AnalysisProducer::ReadContig, m.revision());
    assert!(
        DocumentHistory::derive_from_events(&[], "product", m.clone(), "assembly", origin.clone())
            .is_err()
    );
    assert!(
        DocumentHistory::derive_from_events(
            &[(&a, a.head_event()), (&b, "missing")],
            "product",
            m.clone(),
            "assembly",
            origin.clone(),
        )
        .is_err()
    );
    assert!(
        DocumentHistory::derive_from_events(
            &[(&a, a.head_event()), (&b, b.head_event())],
            "b",
            m.clone(),
            "assembly",
            origin.clone(),
        )
        .is_err()
    );
    let too_many = vec![(&a, a.head_event()); 129];
    assert!(
        DocumentHistory::derive_from_events(&too_many, "product", m, "assembly", origin.clone())
            .is_err()
    );
    assert_eq!(
        before,
        (a.head_event().to_string(), b.head_event().to_string())
    );
}

#[test]
fn in_place_reference_transfer_keeps_both_events_and_undo_redo_after_reload() {
    let mut destination = source("destination");
    let mut original = source("original");
    let original_event = original.head_event().to_owned();
    let edit_end = original.molecule.sequence.len();
    original
        .commit(
            vec![
                Edit::ReplaceSequence {
                    start: 0,
                    end: edit_end,
                    replacement: "AAAA".into(),
                },
                Edit::Rename {
                    name: "later original".into(),
                },
                Edit::SetTopology {
                    topology: Topology::Linear,
                },
            ],
            "Later edit",
        )
        .unwrap();
    let original_head = original.head_event().to_owned();
    let before = destination.head_event().to_owned();
    let copied = original.revision_molecule(&original_event).unwrap();
    let selected = Selection {
        parts: vec![Part {
            start: 0,
            end: 4,
            reverse: false,
        }],
        reverse: false,
    };
    let fragment_revision = reference_transfer::copy_revision(&RequestCopy {
        source: copied.clone(),
        selection: selected.clone(),
    })
    .unwrap();
    let source_event = collagene_core::material::lineage::EventBinding {
        document_id: original.id.clone(),
        event_id: original_event.clone(),
        revision: copied.revision(),
    };
    destination
        .commit_with_source_events(
            vec![
                Edit::PasteReference {
                    source: copied,
                    source_event,
                    selection: selected,
                    fragment_revision,
                    start: 2,
                    end: 2,
                },
                Edit::Rename {
                    name: "reference paste".into(),
                },
            ],
            "Reference copy, not a cloning reaction",
            &[original_event.clone(), original_event.clone()],
        )
        .unwrap();
    let after = destination.head_event().to_owned();
    assert_eq!(
        destination
            .summaries()
            .iter()
            .find(|event| event.current)
            .unwrap()
            .parents,
        vec![before.clone(), original_event]
    );
    let mut frame =
        Workspace::from_documents("transfer", vec![destination.molecule.clone()], 0).unwrap();
    frame.documents.clear();
    frame.revisions.clear();
    frame.history.clear();
    destination.contribute(&mut frame).unwrap();
    original.contribute(&mut frame).unwrap();
    frame.active_document = destination.id.clone();
    let bytes = workspace::save_workspace(&frame).unwrap();
    let loaded = workspace::load_workspace(&bytes).unwrap();
    let mut histories = DocumentHistory::restore_all(&loaded).unwrap();
    let restored = histories
        .iter_mut()
        .find(|doc| doc.id == "destination")
        .unwrap();
    assert!(restored.undo().unwrap());
    assert_eq!(restored.head_event(), before);
    assert_eq!(restored.molecule.sequence.as_str(), "ACGT");
    assert!(restored.redo().unwrap());
    assert_eq!(restored.head_event(), after);
    assert_eq!(restored.molecule.sequence.as_str(), "ACACGTGT");
    assert_eq!(original.head_event(), original_head);
}
#[test]
fn source_only_ancestry_is_reversible_and_bad_parent_fields_fail_atomically() {
    let mut destination = source("destination");
    let original = source("original");
    let before = destination.head_event().to_owned();
    let molecule = destination.molecule.clone();
    destination
        .commit_with_source_events(
            vec![Edit::Rename {
                name: molecule.name.clone(),
            }],
            "Recorded ancillary reference ancestry",
            &[original.head_event().into()],
        )
        .unwrap();
    assert_ne!(destination.head_event(), before);
    assert_eq!(destination.molecule, molecule);
    assert!(destination.undo().unwrap());
    assert_eq!(destination.head_event(), before);
    let count = destination.event_count();
    let navigation = destination.navigation_event_ids();
    assert!(
        destination
            .commit_with_source_events(
                vec![Edit::Rename {
                    name: molecule.name.clone()
                }],
                "invalid",
                &["not-an-event".into()]
            )
            .is_err()
    );
    assert_eq!(destination.head_event(), before);
    assert_eq!(destination.event_count(), count);
    assert_eq!(destination.navigation_event_ids(), navigation);
}
