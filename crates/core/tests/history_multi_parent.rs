use collagene_core::{
    Molecule, Topology,
    history::DocumentHistory,
    workspace::{self, Workspace},
};
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
