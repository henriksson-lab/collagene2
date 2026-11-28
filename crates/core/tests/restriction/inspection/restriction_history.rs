//! Source-role, native-history and real product integration fixtures for the MVP path.
use collagene_core::{
    chemistry,
    cloning::restriction::{history::*, inspection::InspectionOptions},
    documents::workspace::{RevisionEvent, Workspace},
    formats::sequence as formats,
    material::lineage::EventBinding as SourceBinding,
    sequence::analysis,
};
use collagene_model::{Molecule, Sequence, Topology};
const GB: &str = "LOCUS       ZZ123456                  12 bp    DNA     circular SYN 01-JAN-2026\nDEFINITION  Independently authored digest workflow.\nACCESSION   ZZ123456\nVERSION     ZZ123456.1\nFEATURES             Location/Qualifiers\n     misc_feature    join(11..12,1..2)\n                     /label=\"origin crossing\"\nORIGIN\n        1 tttgaattca aa\n//\n";
fn original() -> Molecule {
    formats::import_text(GB).unwrap().molecules.remove(0)
}
fn enzyme() -> Vec<analysis::Enzyme> {
    vec![
        analysis::starter_enzymes()
            .into_iter()
            .find(|e| e.name == "EcoRI")
            .unwrap(),
    ]
}
fn frame() -> Workspace {
    Workspace::from_documents("MVP fixture", vec![original()], 0).unwrap()
}
fn binding(w: &Workspace, index: usize) -> SourceBinding {
    let d = &w.documents[index];
    SourceBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn genbank_to_fresh_circular_digest_retains_annotations_original_and_exact_source_parent() {
    let w = frame();
    let b = binding(&w, 0);
    let before = serde_json::to_value(&w).unwrap();
    let p = prepare_install(&w, &b, &enzyme(), &InspectionOptions::default(), 0).unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    assert_eq!(p.workspace.documents.len(), 2);
    assert_eq!(p.workspace.active_document, p.document_id);
    let product = p.workspace.revisions.get(&p.product_revision).unwrap();
    assert_eq!(product.sequence.as_str(), "AATTCAAATTTG");
    assert_eq!(product.topology, Topology::Linear);
    assert!(product.original_genbank.is_none());
    let feature = product
        .features
        .iter()
        .find(|f| f.label() == "origin crossing")
        .unwrap();
    assert_eq!(feature.location.extract(&product.sequence).unwrap(), "AATT");
    let positions = feature
        .location
        .spans()
        .into_iter()
        .flat_map(|(a, b, _)| a..b)
        .collect::<Vec<_>>();
    assert_eq!(positions, vec![6, 7, 8, 9]);
    let event = p
        .workspace
        .history
        .iter()
        .find(|e| e.id == p.event_id)
        .unwrap();
    assert_eq!(event.parents, vec![b.event_id.clone()]);
    assert_eq!(
        p.workspace.revisions[&b.revision]
            .original_genbank
            .as_deref(),
        Some(GB)
    );
    assert_eq!(
        p.workspace.revisions[&b.revision].sequence.as_str(),
        "TTTGAATTCAAA"
    );
}
#[test]
fn identical_content_roles_and_new_same_content_events_never_reuse_old_binding() {
    let mut w =
        Workspace::from_documents("equal source roles", vec![original(), original()], 0).unwrap();
    let a = binding(&w, 0);
    let b = binding(&w, 1);
    assert_eq!(a.revision, b.revision);
    assert_ne!(source_role(&a).unwrap(), source_role(&b).unwrap());
    assert!(prepare_install(&w, &b, &enzyme(), &InspectionOptions::default(), 0).is_err());
    let old = installation_revision(&w, &a, &enzyme(), &InspectionOptions::default(), 0).unwrap();
    let new = RevisionEvent::new(
        &a.document_id,
        &a.revision,
        vec![a.event_id.clone()],
        "New event with unchanged content",
    )
    .unwrap();
    w.documents[0].head_event = Some(new.id.clone());
    w.history.push(new);
    w.validate().unwrap();
    assert!(prepare_install(&w, &a, &enzyme(), &InspectionOptions::default(), 0).is_err());
    let current = binding(&w, 0);
    assert_ne!(
        old,
        installation_revision(&w, &current, &enzyme(), &InspectionOptions::default(), 0).unwrap()
    );
    let p = prepare_install(&w, &current, &enzyme(), &InspectionOptions::default(), 0).unwrap();
    assert_eq!(
        p.workspace
            .history
            .iter()
            .find(|e| e.id == p.event_id)
            .unwrap()
            .parents,
        vec![current.event_id]
    );
}
#[test]
fn native_roundtrip_keeps_recipe_as_unverified_evidence_and_prepared_binding_is_checked() {
    let w = frame();
    let b = binding(&w, 0);
    let mut p = install(&w, &b, &enzyme(), &InspectionOptions::default(), 0).unwrap();
    let loaded = collagene_core::documents::workspace::load_workspace(&p.archive).unwrap();
    p.validate_archive(&loaded).unwrap();
    let product = &loaded.revisions[&p.product_revision];
    assert!(chemistry::has_unmapped_physical_evidence(product));
    assert!(chemistry::Duplex::from_molecule(product).is_err());
    assert!(
        product
            .metadata
            .keys()
            .any(|k| k.starts_with("restriction_inspection:"))
    );
    assert!(!product.metadata.contains_key("verified"));
    let mut restored =
        collagene_core::documents::history::DocumentHistory::restore_all(&loaded).unwrap();
    let child = restored.iter_mut().find(|d| d.id == p.document_id).unwrap();
    assert_eq!(child.head_event(), p.event_id);
    assert_eq!(child.summaries()[0].parents, vec![b.event_id.clone()]);
    p.source = binding(&loaded, 1);
    assert!(p.validate_archive(&loaded).is_err());
}
#[test]
fn capacity_failure_is_atomic_and_never_drops_original_documents() {
    let w = Workspace::from_documents(
        "at capacity",
        vec![Molecule::dna("tiny", "GAATTC", Topology::Circular).unwrap(); 1000],
        0,
    )
    .unwrap();
    let b = binding(&w, 0);
    let fingerprint = collagene_core::material::lineage::workspace_fingerprint(&w).unwrap();
    assert!(prepare_install(&w, &b, &enzyme(), &InspectionOptions::default(), 0).is_err());
    assert_eq!(w.documents.len(), 1000);
    assert_eq!(
        fingerprint,
        collagene_core::material::lineage::workspace_fingerprint(&w).unwrap()
    );
}
#[test]
fn no_product_is_installed_from_possible_sites_or_truncated_inspection() {
    let mut m = original();
    m.sequence = Sequence::dna("TTTGAANTCAAA").unwrap();
    m.original_genbank = None;
    m.features.clear();
    let w = Workspace::from_documents("uncertain", vec![m], 0).unwrap();
    let b = binding(&w, 0);
    assert!(prepare_install(&w, &b, &enzyme(), &InspectionOptions::default(), 0).is_err());
    let w = frame();
    let b = binding(&w, 0);
    let mut o = InspectionOptions::default();
    o.max_comparisons = 1;
    assert!(prepare_install(&w, &b, &enzyme(), &o, 0).is_err());
    assert_eq!(w.documents.len(), 1);
}
