use collagene_core::documents::{
    history::DocumentHistory,
    workspace::{self, Workspace},
};
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
    a.commit(
        Molecule::dna("newer read state", "AGGT", Topology::Linear).unwrap(),
        "edit",
    )
    .unwrap();
    let before = (
        a.head_event().to_string(),
        a.navigation_event_ids(),
        b.head_event().to_string(),
        b.navigation_event_ids(),
    );
    let product = DocumentHistory::derive_from_events(
        &[(&a, &old), (&b, b.head_event()), (&a, &old)],
        "consensus",
        Molecule::dna("inferred consensus", "ACGTACGT", Topology::Linear).unwrap(),
        "Read-supported computational assembly; not physical synthesis",
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
    assert!(DocumentHistory::derive_from_events(&[], "product", m.clone(), "assembly").is_err());
    assert!(
        DocumentHistory::derive_from_events(
            &[(&a, a.head_event()), (&b, "missing")],
            "product",
            m.clone(),
            "assembly"
        )
        .is_err()
    );
    assert!(
        DocumentHistory::derive_from_events(
            &[(&a, a.head_event()), (&b, b.head_event())],
            "b",
            m.clone(),
            "assembly"
        )
        .is_err()
    );
    let too_many = vec![(&a, a.head_event()); 129];
    assert!(DocumentHistory::derive_from_events(&too_many, "product", m, "assembly").is_err());
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
    original
        .commit(
            Molecule::dna("later original", "AAAA", Topology::Linear).unwrap(),
            "Later edit",
        )
        .unwrap();
    let original_head = original.head_event().to_owned();
    let before = destination.head_event().to_owned();
    destination
        .commit_with_source_events(
            Molecule::dna("reference paste", "ACACGTGT", Topology::Linear).unwrap(),
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
            molecule.clone(),
            "Recorded reference transfer",
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
            .commit_with_source_events(molecule, "invalid", &["not-an-event".into()])
            .is_err()
    );
    assert_eq!(destination.head_event(), before);
    assert_eq!(destination.event_count(), count);
    assert_eq!(destination.navigation_event_ids(), navigation);
}
