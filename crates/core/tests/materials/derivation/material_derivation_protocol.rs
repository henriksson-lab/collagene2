//! Native material selectors, source identity and current encoding integrity.
use collagene_core::{
    analysis::Primer,
    catalog_digest, material_derivation as native,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    run_lineage::EventBinding,
    run_material::MaterialOrigin,
    transfer,
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
use collagene_model::{Molecule, Topology};
const TEMPLATE: &str = "ACGACCCGGGATCGATCCGTA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn tag(v: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(v).unwrap())
        .unwrap()
        .0
}
fn call(
    job: AnalysisJob,
    revision: String,
    job_tag: u32,
    result_tag: u32,
) -> protocol::response::Response {
    assert_eq!(tag(&job), job_tag);
    let bytes = protocol::request::Request::encode(&Request {
        id: 177,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&90u16.to_le_bytes());
    assert_eq!(&bytes[6..], &old[6..]);
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let reply = protocol::response::Response::decode(&bytes).unwrap();
    assert_eq!(reply.revision, revision);
    assert_eq!(tag(&reply.result), result_tag, "{:?}", reply.result);
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&90u16.to_le_bytes());
    assert_eq!(&bytes[6..], &old[6..]);
    assert!(protocol::response::Response::decode(&old).is_err());
    reply
}
fn installed() -> (Workspace, EventBinding) {
    let source = Molecule::dna("Actual template", TEMPLATE, Topology::Linear).unwrap();
    let w = Workspace::from_documents("PCR source", vec![source.clone()], 0).unwrap();
    let forward = Primer::new("F", "ACGAC", "").unwrap();
    let reverse = Primer::new("R", "TACGG", "").unwrap();
    let mut previews = verification::simulate_pcr_documents(
        &source,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(previews.len(), 1);
    let p = previews.remove(0);
    let request = native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: event(&w, 0),
            origin: native::SourceOrigin::OriginalReference,
        }],
        derivation: native::Derivation::ExactPcr {
            forward,
            reverse,
            selection: native::PcrSelection {
                selected_index: 0,
                expected_preview_revision: p.molecule.revision(),
                geometry: p.geometry,
            },
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        name: "Reusable amplicon".into(),
        save_primers: true,
    };
    let revision = native::revision(&w, &request).unwrap();
    let reply = call(
        AnalysisJob::InstallMaterialDerivation {
            workspace: w,
            request,
        },
        revision.clone(),
        171,
        161,
    );
    let AnalysisResult::MaterialDerivationInstalled(value) = reply.result else {
        unreachable!()
    };
    value.validate_binding(&revision).unwrap();
    let w = workspace::load_workspace(&value.archive).unwrap();
    value.receipt.validate_archive(&w).unwrap();
    assert_eq!(
        w.revisions[&value.receipt.product.revision]
            .sequence
            .as_str(),
        TEMPLATE
    );
    assert!(value.receipt.primer_attachment.is_some());
    (w, value.receipt.product)
}
#[test]
fn current_install_and_digest_jobs_bind_actual_material() {
    let (w, product) = installed();
    let source = w.revisions[&product.revision].clone();
    let origin = native::SourceOrigin::Saved(MaterialOrigin::Derived(
        native::capsule(&w, &product).unwrap(),
    ));
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "SmaI")
        .unwrap()
        .id
        .clone();
    let selection = catalog_digest::select(&catalog, &[enzyme]).unwrap();
    let options = Default::default();
    let role = collagene_core::restriction_history::source_role(&product).unwrap();
    let context =
        native::digest_inspection_revision(&role, &source, &origin, &selection, &options).unwrap();
    let response = call(
        AnalysisJob::InspectMaterialDigest {
            source_role: role,
            source: source.clone(),
            origin,
            selection,
            options,
        },
        context.clone(),
        172,
        162,
    );
    let AnalysisResult::MaterialDigestInspection(report) = response.result else {
        unreachable!()
    };
    report.validate_binding(&context).unwrap();
    assert!(report.inspection.products_determined);
    assert_eq!(report.inspection.fragments.len(), 2);
    assert!(!report.incomplete);
    assert_eq!(source, w.revisions[&product.revision]);
}
#[test]
fn reconstruction_operation_and_evidence_are_explicit_and_replayable() {
    let (w, product) = installed();
    let source = w.revisions[&product.revision].clone();
    let target = Molecule::dna("Independent expected DNA", TEMPLATE, Topology::Linear).unwrap();
    let origin = native::capsule(&w, &product).unwrap();
    let plan = DesignPlan {
        name: "Retrieve native material against independent DNA".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "amplicon".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::ReconstructMaterial {
            template: Source::Input("amplicon".into()),
            origin,
        }],
    };
    assert_eq!(tag(&plan.steps[0]), 26);
    assert_eq!(protocol::schema::CANDIDATE_SCHEMA_VERSION, 31);
    assert_eq!(workspace::PLAN_SCHEMA_VERSION, 29);
    let bytes = transfer::save_candidate(&transfer::Candidate::new(plan.clone())).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&(protocol::schema::CANDIDATE_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(transfer::load_candidate(&old).is_err());
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert!(
        report
            .engine_version
            .eq(collagene_core::verification::ENGINE_VERSION)
    );
    assert_eq!(tag(&report.physical_processing_reports[0].evidence), 6);
}
#[test]
fn dropping_only_duplex_summary_cannot_launder_native_material_into_a_plain_source() {
    let (w, product) = installed();
    let mut source = w.revisions[&product.revision].clone();
    source.metadata.remove("duplex_model");
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &source
    ));
    assert!(
        native::preview_pcr(
            &source,
            &native::SourceOrigin::OriginalReference,
            &Primer::new("F", "ACGAC", "").unwrap(),
            &Primer::new("R", "TACGG", "").unwrap(),
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    assert!(
        collagene_core::snapgene_export::export_snapgene(
            &source,
            collagene_core::snapgene_export::Policy::AllowOmissions
        )
        .is_err()
    );
    let origin = MaterialOrigin::Derived(native::capsule(&w, &product).unwrap());
    let mut reference = collagene_core::material_pcr::simulate(
        &w.revisions[&product.revision],
        &origin,
        &Primer::new("F", "ACGAC", "").unwrap(),
        &Primer::new("R", "TACGG", "").unwrap(),
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap()
    .remove(0)
    .molecule;
    reference.metadata.remove("duplex_model");
    assert!(
        reference
            .metadata
            .contains_key(collagene_core::material_pcr::ORIGIN_METADATA)
    );
    assert!(
        native::preview_pcr(
            &reference,
            &native::SourceOrigin::OriginalReference,
            &Primer::new("F", "ACGAC", "").unwrap(),
            &Primer::new("R", "TACGG", "").unwrap(),
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
}

#[test]
fn ordinary_preparation_compiler_owns_reconstruction_new_pcr_and_physical_assembly() {
    use collagene_core::assembly_recipe::{self as a, FragmentSpec, Preparation, RunPreparation};
    let (w, product) = installed();
    let amp = w.revisions[&product.revision].clone();
    let other = Molecule::dna("Second original", "CGTATTGGAACG", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent assembly intent",
        "ACGACCCGGGATCGATCCGTATTGGAACG",
        Topology::Linear,
    )
    .unwrap();
    let inputs = vec![
        PlanInput {
            id: "amplicon".into(),
            revision: amp.revision(),
            molecule: amp,
        },
        PlanInput {
            id: "other".into(),
            revision: other.revision(),
            molecule: other,
        },
    ];
    let recipe = a::Recipe {
        name: "Fresh preparation from reusable material".into(),
        fragments: vec![
            FragmentSpec {
                generation: None,
                role: "amplicon-role".into(),
                input_id: "amplicon".into(),
                reverse: false,
                preparation: Preparation::FromDerived {
                    origin: native::capsule(&w, &product).unwrap(),
                    preparation: RunPreparation::ExactPcr {
                        forward: Primer::new("New F", "ACGAC", "").unwrap(),
                        reverse: Primer::new("New R", "TACGG", "").unwrap(),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    },
                },
            },
            FragmentSpec {
                generation: None,
                role: "other-role".into(),
                input_id: "other".into(),
                reverse: false,
                preparation: Preparation::ExactPcr {
                    forward: Primer::new("Other F", "CGTA", "").unwrap(),
                    reverse: Primer::new("Other R", "CGTT", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
            },
        ],
        overlaps: vec![4],
        circular: false,
        options: Default::default(),
    };
    let policy = collagene_core::homology_assembly::Options {
        overlaps: vec![4],
        topology: Topology::Linear,
        min_paired_bases: 4,
        ..Default::default()
    };
    let preview = a::preview_homology(&inputs, &target, &recipe, &policy).unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "ACGACCCGGGATCGATCCGTATTGGAACG"
    );
    assert!(matches!(
        preview.compiled.plan.steps[0],
        Operation::ReconstructMaterial { .. }
    ));
    let orders =
        collagene_core::plan_oligos::build(&preview.compiled.plan, &target, &Default::default())
            .unwrap();
    assert_eq!(orders.oligos.records.len(), 4);
    let verification::PhysicalProcessingEvidence::DerivedOrigin(evidence) =
        &preview.verification.physical_processing_reports[0].evidence
    else {
        unreachable!()
    };
    let mut excerpt = evidence.clone();
    excerpt.material.origin_findings = vec!["A".repeat(2048)];
    excerpt.material.origin_findings_total = 1;
    excerpt.material.origin_findings_truncated = true;
    excerpt.validate_documentary().unwrap();
    excerpt.material.origin_findings_total = 2;
    excerpt.material.origin_findings_truncated = false;
    assert!(excerpt.validate_documentary().is_err());
}
