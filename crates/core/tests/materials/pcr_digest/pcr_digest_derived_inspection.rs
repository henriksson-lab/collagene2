//! Actual retained unequal strands survive a target-free derivation and restart.
#[path = "../../common/material_library.rs"]
mod material_library;
use collagene_core::{
    chemistry::{FivePrime, Strand},
    cloning::{pcr::digest as inspection, restriction::catalog as catalog_digest},
    documents::workspace,
    material::{derivation, run::MaterialOrigin},
    sequence::analysis::Primer,
};
fn selection(name: &str) -> catalog_digest::CatalogSelection {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == name)
        .unwrap()
        .id
        .clone();
    catalog_digest::select(&catalog, &[id]).unwrap()
}
fn fixture() -> (workspace::Workspace, MaterialOrigin, inspection::Request) {
    let (w, event, origin) =
        material_library::native_source("AGTCCGATGCTAC", "GTAGGTAGCATCGGACT", 13);
    let prepared = derivation::prepare_install(
        &w,
        &derivation::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: event,
                origin: derivation::SourceOrigin::Saved(origin),
            }],
            derivation: derivation::Derivation::CatalogDigest {
                conditions: None,
                selection: selection("SmaI"),
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Retained unequal native fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let w = workspace::load_workspace(&workspace::save_workspace(&prepared.workspace).unwrap())
        .unwrap();
    let origin =
        MaterialOrigin::Derived(derivation::capsule(&w, &prepared.receipt.product).unwrap());
    let request = inspection::Request {
        source_role: "actual-bottom-template".into(),
        template_strand: Strand::Bottom,
        forward: Primer::new("Bottom F", "GTAGGTAG", "CCGAATTC").unwrap(),
        reverse: Primer::new("Bottom R", "AGTCCGAT", "GGGAATTC").unwrap(),
        selection: selection("EcoRI"),
        options: Default::default(),
    };
    (w, origin, request)
}
#[test]
fn derived_bottom_pcr_digest_uses_extra_actual_calls_after_restart() {
    let (w, origin, request) = fixture();
    let source = &w.revisions[&origin.event().revision];
    let before = workspace::save_workspace(&w).unwrap();
    let fresh = origin.replay(source).unwrap();
    assert_eq!(fresh.tracked().duplex.summary().unwrap().top_length, 13);
    assert_eq!(fresh.tracked().duplex.summary().unwrap().bottom_length, 17);
    let report = inspection::inspect_from_material(source, &origin, &request, &w).unwrap();
    report
        .validate_binding(
            &inspection::inspection_revision_from_material(source, &origin, &request, &w).unwrap(),
        )
        .unwrap();
    assert_eq!(report.preview.pcr_geometries.len(), 1);
    assert_eq!(
        report.preview.pcr_geometries[0].template_strand,
        Strand::Bottom
    );
    assert_eq!(
        report.preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCGTAGGTAGCATCGGACTGAATTCCC"
    );
    let digest = report.preview.inspection.as_ref().unwrap();
    assert_eq!(digest.fragments.len(), 3);
    assert_eq!(
        digest.fragments[0]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        FivePrime::Hydroxyl
    );
    assert_eq!(
        digest.fragments[1]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        FivePrime::Phosphate
    );
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    assert!(inspection::inspect(source, &request).is_err());
}
#[test]
fn derived_source_owner_alias_primers_and_competing_geometries_are_not_selected_implicitly() {
    let (w, origin, request) = fixture();
    let source = &w.revisions[&origin.event().revision];
    let report = inspection::inspect_from_material(source, &origin, &request, &w).unwrap();
    let foreign =
        workspace::Workspace::from_documents("Equal DNA is not ownership", vec![source.clone()], 0)
            .unwrap();
    assert!(inspection::inspect_from_material(source, &origin, &request, &foreign).is_err());
    let MaterialOrigin::Derived(mut alias) = origin.clone() else {
        panic!("Derived")
    };
    alias.product.event_id.push_str("-alias");
    assert!(
        inspection::preflight_from_material(source, &MaterialOrigin::Derived(alias), &request, &w)
            .is_err()
    );
    let mut changed = request.clone();
    changed.forward.tail = "AAGAATTC".into();
    assert!(
        report
            .validate_binding(
                &inspection::inspection_revision_from_material(source, &origin, &changed, &w)
                    .unwrap()
            )
            .is_err()
    );
    let mut ambiguous = request.clone();
    ambiguous.forward = Primer::new("F", "G", "").unwrap();
    ambiguous.reverse = Primer::new("R", "C", "").unwrap();
    // Whole-frame geometry remains native: ambiguity or out-of-frame extension
    // is explicit; no desired target is supplied to select a result.
    match inspection::inspect_from_material(source, &origin, &ambiguous, &w) {
        Ok(p) => {
            assert_ne!(p.preview.pcr_geometries.len(), 1);
            assert!(p.preview.inspection.is_none());
        }
        Err(e) => assert!(matches!(e, collagene_model::Error::Unsupported(_))),
    }
}
