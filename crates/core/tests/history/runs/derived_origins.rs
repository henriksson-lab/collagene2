//! Actual native target-free installations exercise mixed material ancestry.
use collagene_core::{
    analysis::Primer,
    catalog_digest,
    material_derivation::{self as native, Derivation, PcrSelection, SourceOrigin},
    run_lineage::EventBinding,
    run_material::{self, MaterialOrigin, MaterialOriginRef, MaterialSelection},
    verification,
    workspace::{self, RevisionEvent, Workspace, WorkspaceDocument},
};
use collagene_model::{Error, Molecule, Topology};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn original() -> Workspace {
    Workspace::from_documents(
        "Original native reference",
        vec![Molecule::dna("Stock", STOCK, Topology::Linear).unwrap()],
        0,
    )
    .unwrap()
}
fn pcr(w: &Workspace) -> native::PreparedInstallation {
    let source = w.head(&w.documents[0].id).unwrap();
    let forward = Primer::new("Forward", "AGTCCGAT", "").unwrap();
    let reverse = Primer::new("Reverse", "TGACGATC", "").unwrap();
    let mut previews = verification::simulate_pcr_documents(
        source,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(previews.len(), 1);
    let preview = previews.remove(0);
    native::prepare_install(
        w,
        &native::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: event(w, 0),
                origin: SourceOrigin::OriginalReference,
            }],
            derivation: Derivation::ExactPcr {
                forward,
                reverse,
                selection: PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: preview.molecule.revision(),
                    geometry: preview.geometry,
                },
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            name: "Actual target-free PCR".into(),
            save_primers: false,
        },
    )
    .unwrap()
}
fn digest_request(source: EventBinding, origin: MaterialOrigin, index: usize) -> native::Request {
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: source,
            origin: SourceOrigin::Saved(origin),
        }],
        derivation: Derivation::CatalogDigest {
            selection: catalog_digest::select(&c, &[id]).unwrap(),
            options: Default::default(),
            fragment_index: index,
        },
        name: "Actual target-free digest".into(),
        save_primers: false,
    }
}
#[test]
fn exact_derived_event_is_offered_and_equal_content_document_alias_cannot_select_it() {
    let installed = pcr(&original());
    let product = installed.receipt.product.clone();
    let inspection = run_material::inspect_material(&installed.workspace, &product, &None).unwrap();
    assert_eq!(
        inspection.origins,
        vec![MaterialSelection::Derived(product.clone())]
    );
    let chosen = run_material::inspect_material(
        &installed.workspace,
        &product,
        &Some(MaterialSelection::Derived(product.clone())),
    )
    .unwrap();
    let origin = chosen.selected.unwrap();
    origin
        .validate_for_workspace(&installed.workspace, &product)
        .unwrap();
    let mut w = installed.workspace;
    let alias_id = "same-content-different-document".to_string();
    let alias_event = RevisionEvent::new(
        &alias_id,
        &product.revision,
        vec![],
        "Imported equal-content reference",
    )
    .unwrap();
    let alias = EventBinding {
        document_id: alias_id.clone(),
        event_id: alias_event.id.clone(),
        revision: product.revision.clone(),
    };
    w.documents.push(WorkspaceDocument {
        id: alias_id,
        head: product.revision.clone(),
        undo: vec![],
        redo: vec![],
        head_event: Some(alias_event.id.clone()),
        undo_events: vec![],
        redo_events: vec![],
    });
    w.history.push(alias_event);
    w.validate().unwrap();
    assert!(
        run_material::inspect_material(&w, &alias, &None)
            .unwrap()
            .origins
            .is_empty()
    );
    assert!(origin.validate_for_workspace(&w, &alias).is_err());
    assert!(
        run_material::inspect_material(&w, &alias, &Some(MaterialSelection::Derived(product)))
            .is_err()
    );
}
#[test]
fn nested_pcr_then_digest_keeps_exact_sources_and_actual_staggered_state_after_archive_reload() {
    let first = pcr(&original());
    let origin =
        MaterialOrigin::Derived(native::capsule(&first.workspace, &first.receipt.product).unwrap());
    let second = native::prepare_install(
        &first.workspace,
        &digest_request(first.receipt.product.clone(), origin, 1),
    )
    .unwrap();
    let bytes = workspace::save_workspace(&second.workspace).unwrap();
    let restored = workspace::load_workspace(&bytes).unwrap();
    let product = &second.receipt.product;
    let source = restored.revisions[&product.revision].clone();
    let origin = MaterialOrigin::Derived(native::capsule(&restored, product).unwrap());
    origin.validate_for_workspace(&restored, product).unwrap();
    let fresh = origin.replay(&source).unwrap();
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert!(!fresh.incomplete());
    let before = source.clone();
    let forward = Primer::new("New forward", "AATTCCGT", "").unwrap();
    let reverse = Primer::new("New reverse", "TGACGATC", "").unwrap();
    let previews = collagene_core::material_pcr::simulate(
        &source,
        &origin,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(previews.len(), 1);
    assert_eq!(
        previews[0].molecule.sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    assert_eq!(source, before);
    let MaterialOrigin::Derived(capsule) = &origin else {
        panic!("wrong origin")
    };
    let record = native::load_record(&capsule.record_bytes).unwrap();
    assert!(
        native::source_bindings(&record)
            .unwrap()
            .contains(&first.receipt.product)
    );
    assert_eq!(
        native::required_derivation_products(&record).unwrap(),
        vec![first.receipt.product]
    );
}
#[test]
fn one_shared_ledger_limits_all_roots_and_refuses_oversized_opaque_package_before_decode() {
    let first = pcr(&original());
    let capsule = native::capsule(&first.workspace, &first.receipt.product).unwrap();
    let origin = MaterialOrigin::Derived(capsule.clone());
    assert_eq!(run_material::MAX_AGGREGATE_ORIGIN_BYTES, 32 * 1024 * 1024);
    run_material::material_origins_preflight(std::iter::repeat_n(&origin, 1000)).unwrap();
    assert!(matches!(
        run_material::material_origins_preflight(std::iter::repeat_n(&origin, 1001)),
        Err(Error::TooLarge)
    ));
    let mut oversized = capsule;
    oversized.record_bytes = vec![0; native::MAX_RECORD_BYTES + 1];
    assert!(matches!(
        run_material::material_capsules_preflight([MaterialOriginRef::Derived(&oversized)]),
        Err(Error::TooLarge)
    ));
}
#[test]
fn opaque_derived_nesting_cannot_reset_the_native_eight_generation_bound() {
    let first = pcr(&original());
    let mut w = first.workspace;
    let mut product = first.receipt.product;
    // First transformation atdepth1 plus seven actual catalog transformations.
    // After the first cleavage the chosen fragment has no remaining EcoRI site;
    // subsequent transforms retain those actual unchanged free ends explicitly.
    for layer in 2..=8 {
        let origin = MaterialOrigin::Derived(native::capsule(&w, &product).unwrap());
        let next = native::prepare_install(
            &w,
            &digest_request(product, origin, if layer == 2 { 1 } else { 0 }),
        )
        .unwrap();
        w = next.workspace;
        product = next.receipt.product;
    }
    let origin = MaterialOrigin::Derived(native::capsule(&w, &product).unwrap());
    let request = digest_request(product, origin, 0);
    let before = workspace::save_workspace(&w).unwrap();
    assert!(native::prepare_install(&w, &request).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
