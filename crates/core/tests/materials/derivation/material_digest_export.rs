#[allow(dead_code)]
#[path = "../../common/material_library.rs"]
mod fixture;
use collagene_core::{
    cloning::restriction::{
        catalog,
        setups::{CleavageMode, Recipe},
    },
    formats::sequence::{self, TextFormat},
    material::derivation::{self as d, digest_export},
    protocol::{
        execution::process_request,
        request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
};

fn request(
    source: collagene_core::material::lineage::EventBinding,
    origin: d::SourceOrigin,
) -> digest_export::Request {
    let catalog = catalog::builtin_catalog().unwrap();
    let eco = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    digest_export::Request {
        recipe: Recipe {
            source,
            origin,
            selection: catalog::select(&catalog, &[eco.id.clone()]).unwrap(),
            options: Default::default(),
            cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
            conditions: None,
        },
        fragment_index: 1,
        name: "Detached right fragment".into(),
    }
}

#[test]
fn saved_staggered_fragment_export_uses_actual_strands_and_discloses_detachment() {
    let (w, source, origin) =
        fixture::native_source("AATTGCTACGGAATTCTTACG", "AATTCGTAAGAATTCCGTAGC", 17);
    let molecule = &w.revisions[&source.revision];
    let req = request(source, d::SourceOrigin::Saved(origin));
    let product = req.materialize(molecule).unwrap();
    assert_eq!(product.molecule.sequence.as_str(), "AATTCTTACG");
    assert_eq!(product.duplex.top_length, 10);
    assert_eq!(product.duplex.bottom_length, 10);
    assert_eq!(product.duplex.bottom_start, 4);
    assert_eq!(product.molecule.name, req.name);
    assert!(
        product
            .findings
            .iter()
            .any(|s| s.contains("both strand sequences"))
    );
    let gb = sequence::export_with_report(&product.molecule, TextFormat::GenBank).unwrap();
    assert!(
        gb.warnings
            .iter()
            .any(|s| s.contains("reusable physical material"))
    );
    assert!(gb.text.contains("source:material_digest_export:losses"));
    let imported = sequence::import_text(&gb.text).unwrap();
    assert_eq!(imported.molecules[0].sequence, product.molecule.sequence);
    let recognition = d::digest::inspect_digest(
        "unattested reimport",
        &imported.molecules[0],
        &d::SourceOrigin::OriginalReference,
        &req.recipe.selection,
        &req.recipe.options,
    )
    .unwrap();
    assert!(recognition.source_material.is_none());
    assert!(!recognition.inspection.products_determined);
    let detached_frame = collagene_core::documents::workspace::Workspace::from_documents(
        "Detached reference",
        vec![imported.molecules[0].clone()],
        0,
    )
    .unwrap();
    let detached_request = request(
        fixture::event(&detached_frame, 0),
        d::SourceOrigin::OriginalReference,
    );
    assert!(
        detached_request
            .materialize(&imported.molecules[0])
            .is_err()
    );
}

#[test]
fn derived_origin_and_worker_roundtrip_keep_exact_export_binding() {
    let (w, source, origin) =
        fixture::native_source("AATTGCTACGGAATTCTTACG", "AATTCGTAAGAATTCCGTAGC", 17);
    let mut req = request(source.clone(), d::SourceOrigin::Saved(origin));
    let native = d::Request {
        sources: vec![d::BoundSource {
            role: d::SourceRole::Template,
            event: source,
            origin: req.recipe.origin.clone(),
        }],
        derivation: d::Derivation::CatalogDigest {
            conditions: None,
            selection: req.recipe.selection.clone(),
            options: req.recipe.options.clone(),
            fragment_index: 1,
        },
        name: "Retained fragment".into(),
        save_primers: false,
    };
    let installed = d::prepare_install(&w, &native).unwrap();
    req.recipe.source = installed.receipt.product.clone();
    req.recipe.origin =
        d::SourceOrigin::Saved(collagene_core::material::run::MaterialOrigin::Derived(
            d::capsule(&installed.workspace, &req.recipe.source).unwrap(),
        ));
    req.fragment_index = 0;
    let source = installed.workspace.revisions[&req.recipe.source.revision].clone();
    let revision = req.revision(&source).unwrap();
    let bytes = request::Request {
        id: 43,
        revision: revision.clone(),
        job: AnalysisJob::MaterializeMaterialDigestFragment {
            source: source.clone(),
            request: req.clone(),
        },
    }
    .encode()
    .unwrap();
    let ((id, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!((id, tag), (43, 197));
    let mut superseded = bytes.clone();
    superseded[4..6].copy_from_slice(&121u16.to_le_bytes());
    assert!(process_request(&superseded).is_err());
    let wire = process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    assert_eq!(tag, 187);
    let response = Response::decode(&wire).unwrap();
    assert_eq!(response.revision, revision);
    let AnalysisResult::MaterialDigestProduct(product) = response.result else {
        panic!("unexpected export response")
    };
    product.validate_bindings_for(&req, &source).unwrap();
    assert_eq!(product.molecule.sequence.as_str(), "AATTCTTACG");
    req.name = "Different request".into();
    assert!(product.validate_bindings_for(&req, &source).is_err());
}

#[test]
fn export_refuses_wrong_event_limits_and_changed_worker_identity() {
    let (w, source, origin) =
        fixture::native_source("AATTGCTACGGAATTCTTACG", "AATTCGTAAGAATTCCGTAGC", 17);
    let molecule = &w.revisions[&source.revision];
    let mut req = request(source, d::SourceOrigin::Saved(origin));
    req.recipe.source.event_id = format!("event-{}", "0".repeat(64));
    assert!(req.revision(molecule).is_err());
    let d::SourceOrigin::Saved(origin) = &req.recipe.origin else {
        unreachable!()
    };
    req.recipe.source = origin.event().clone();
    req.recipe.options.max_comparisons = 1;
    assert!(req.materialize(molecule).is_err());
    req.recipe.options = Default::default();
    let bytes = request::Request {
        id: 1,
        revision: "0".repeat(64),
        job: AnalysisJob::MaterializeMaterialDigestFragment {
            source: molecule.clone(),
            request: req,
        },
    }
    .encode()
    .unwrap();
    assert!(process_request(&bytes).is_err());
}

#[test]
fn ambiguous_reference_retains_recognition_without_inventing_exportable_strands() {
    use collagene_core::{cloning::restriction::setups, documents::workspace::Workspace};
    use collagene_model::{Molecule, Topology};
    let molecule = Molecule::dna(
        "Unresolved EcoRI reference",
        "TTTGAANTCAAA",
        Topology::Linear,
    )
    .unwrap();
    let frame =
        Workspace::from_documents("Ambiguous reference", vec![molecule.clone()], 0).unwrap();
    let req = request(
        fixture::event(&frame, 0),
        d::SourceOrigin::OriginalReference,
    );
    let inspected = req.recipe.inspect(&molecule).unwrap();
    assert_eq!(inspected.inspection.sites.len(), 1);
    assert!(!inspected.inspection.sites[0].certain);
    assert!(inspected.source_material.is_none());
    assert!(!inspected.inspection.products_determined);
    assert!(inspected.inspection.fragments.is_empty());
    assert!(req.materialize(&molecule).is_err());
    let (saved, hash) = setups::capture(
        &catalog::builtin_catalog().unwrap(),
        &frame,
        &req.recipe,
        "Unresolved recognition",
        &inspected,
    )
    .unwrap();
    assert!(
        setups::load(&saved, &hash)
            .unwrap()
            .report
            .source_material
            .is_none()
    );
    let revision = req.recipe.inspection_revision(&molecule).unwrap();
    let bytes = request::Request {
        id: 7,
        revision: revision.clone(),
        job: AnalysisJob::InspectMaterialDigest {
            source_role: inspected.inspection.source_role.clone(),
            source: molecule,
            origin: req.recipe.origin,
            selection: req.recipe.selection,
            options: req.recipe.options,
        },
    }
    .encode()
    .unwrap();
    let result = Response::decode(&process_request(&bytes).unwrap()).unwrap();
    let AnalysisResult::MaterialDigestInspection(report) = result.result else {
        panic!("expected recognition evidence")
    };
    assert_eq!(report.context_revision, revision);
    assert!(report.source_material.is_none());
    assert!(report.incomplete);
}

#[test]
fn reaction_choices_bind_detached_provenance_without_changing_assumed_cleavage() {
    use collagene_core::cloning::restriction::conditions::{ChosenFormulation, ReactionConditions};
    let (w, source, origin) =
        fixture::native_source("AATTGCTACGGAATTCTTACG", "AATTCGTAAGAATTCCGTAGC", 17);
    let molecule = &w.revisions[&source.revision];
    let mut req = request(source, d::SourceOrigin::Saved(origin));
    let reference = req.materialize(molecule).unwrap();
    let scientific_revision = req.recipe.inspection_revision(molecule).unwrap();
    req.recipe.conditions = Some(ReactionConditions {
        schema_version: 1,
        formulations: vec![ChosenFormulation {
            enzyme_id: "ecori".into(),
            formulation_id: Some("neb-r3101".into()),
        }],
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
        incubation_minutes: 30,
        minimum_activity_percent: 75,
        activity_context: Some("NEB unit-assay relative activity".into()),
        temperature_context: Some("vendor recommended incubation".into()),
    });
    let chosen = req.materialize(molecule).unwrap();
    assert_ne!(reference.context_revision, chosen.context_revision);
    assert_eq!(chosen.inspection_revision, scientific_revision);
    assert_eq!(chosen.molecule.sequence, reference.molecule.sequence);
    assert_eq!(
        serde_json::to_value(&chosen.duplex).unwrap(),
        serde_json::to_value(&reference.duplex).unwrap()
    );
    let expected = serde_json::to_string(req.recipe.conditions.as_ref().unwrap()).unwrap();
    assert_eq!(
        chosen.molecule.metadata["source:restriction_conditions:options"],
        expected
    );
    for format in [TextFormat::GenBank, TextFormat::Embl] {
        let output = sequence::export_with_report(&chosen.molecule, format).unwrap();
        assert!(output.text.contains("restriction_conditions:options="));
        assert!(output.text.contains("unverified producer data"));
        let imported = sequence::import_text(&output.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(imported.sequence, chosen.molecule.sequence);
        assert!(
            !imported
                .metadata
                .contains_key("restriction_conditions:verified")
        );
    }
    let mut changed = req.clone();
    changed
        .recipe
        .conditions
        .as_mut()
        .unwrap()
        .incubation_minutes = 31;
    assert_eq!(
        changed.recipe.inspection_revision(molecule).unwrap(),
        scientific_revision
    );
    assert_ne!(changed.revision(molecule).unwrap(), chosen.context_revision);
    assert!(chosen.validate_bindings_for(&changed, molecule).is_err());
    let mut tampered = chosen.clone();
    tampered
        .molecule
        .metadata
        .remove("source:restriction_conditions:options");
    assert!(tampered.validate_bindings_for(&req, molecule).is_err());
}
