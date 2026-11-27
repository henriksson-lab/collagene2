//! Current worker material-origin reconstruction, identity and output fixtures.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest, cloning_outputs,
    history::DocumentHistory,
    material_derivation as native,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin, MaterialSelection},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const CUT: &str = "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const AMPLICON: &str = "CCGTAGGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn digest_fixture() -> (Workspace, native::Capsule, EventBinding) {
    let source = Molecule::dna("Original stock", STOCK, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent original PCR target", STOCK, Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Mixed origin literal",
        vec![source.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Original actual PCR".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source)],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("Original F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Original R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved full PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selected = run_material::Selection {
        run_id: saved.receipt.run.run_id.clone(),
        step_index: 0,
        event: saved.receipt.construct.clone(),
    };
    let origin = run_material::capsule(&saved.workspace, &selected).unwrap();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let prepared = native::prepare_install(
        &saved.workspace,
        &native::Request {
            source: selected.event,
            origin: native::SourceOrigin::Saved(MaterialOrigin::Run(origin)),
            derivation: native::Derivation::CatalogDigest {
                selection: catalog_digest::select(&catalog, &[enzyme]).unwrap(),
                options: Default::default(),
                fragment_index: 1,
            },
            name: "Opened EcoRI fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let pin = prepared.receipt.product;
    let cap = native::capsule(&prepared.workspace, &pin).unwrap();
    (prepared.workspace, cap, pin)
}
fn prepared_plan(
    w: &mut Workspace,
    cap: &native::Capsule,
    pin: &EventBinding,
) -> (DesignPlan, EventBinding) {
    let target =
        Molecule::dna("Independent second PCR target", AMPLICON, Topology::Linear).unwrap();
    DocumentHistory::new("independent-second-target", target.clone())
        .unwrap()
        .contribute(w)
        .unwrap();
    let target_pin = event(w, w.documents.len() - 1);
    let source = w.revisions[&pin.revision].clone();
    assert_eq!(source.sequence.as_str(), CUT);
    let plan = DesignPlan {
        name: "PCR from exact saved digest fragment".into(),
        target_revision: target.revision(),
        inputs: vec![input("fragment", source)],
        steps: vec![
            Operation::ReconstructMaterial {
                template: Source::Input("fragment".into()),
                origin: cap.clone(),
            },
            Operation::Pcr {
                template: Source::Step(0),
                forward: Primer::new("New F", "CCGTAGGA", "").unwrap(),
                reverse: Primer::new("New R", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        ],
    };
    (plan, target_pin)
}

use collagene_core::{
    material_pcr, pcr_outputs,
    protocol::{self, AnalysisJob, AnalysisResult},
    transfer,
};
fn wire<T: serde::Serialize>(magic: &[u8; 4], version: u16, payload: &T) -> Vec<u8> {
    let payload = postcard::to_allocvec(payload).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&version.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    bytes
}
fn tag<T: serde::Serialize>(value: &T) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(value).unwrap())
        .unwrap()
        .0
}
fn checked_roundtrip(
    job: AnalysisJob,
    revision: String,
    job_tag: u32,
    result_tag: u32,
) -> protocol::Response {
    let request = protocol::Request {
        id: job_tag,
        revision: revision.clone(),
        job,
    };
    assert_eq!(tag(&request.job), job_tag);
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        protocol::WORKER_SCHEMA_VERSION
    );
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(protocol::process_request(&old).is_err());
    let encoded = protocol::process_request(&bytes).unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::decode_response(&encoded).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, revision);
    assert_eq!(tag(&response.result), result_tag);
    assert!(
        !matches!(response.result, AnalysisResult::Error(_)),
        "{:?}",
        response.result
    );
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &encoded[6..]);
    assert!(protocol::decode_response(&old).is_err());
    response
}
#[test]
fn current_worker_operations_reconstruct_literal_material_and_bind_outputs() {
    let (w, cap, pin) = digest_fixture();
    let before = serde_json::to_value(&w).unwrap();
    let selection = Some(MaterialSelection::Derived(pin.clone()));
    let revision = run_material::material_inspection_revision(&w, &pin, &selection).unwrap();
    let inspected = checked_roundtrip(
        AnalysisJob::InspectMaterialOrigin {
            workspace: w.clone(),
            event: pin.clone(),
            selection,
        },
        revision.clone(),
        168,
        158,
    );
    let AnalysisResult::MaterialOriginInspection(report) = inspected.result else {
        panic!("material origin inspection")
    };
    report.validate_binding(&revision).unwrap();
    assert!(matches!(report.selected, Some(MaterialOrigin::Derived(_))));
    let request = native::load_record(&cap.record_bytes).unwrap().request;
    let revision = native::revision(&w, &request).unwrap();
    let installed = checked_roundtrip(
        AnalysisJob::InstallMaterialDerivation {
            workspace: w.clone(),
            request,
        },
        revision.clone(),
        174,
        164,
    );
    let AnalysisResult::MaterialDerivationInstalled(installed) = installed.result else {
        panic!("digest installation")
    };
    installed.validate_binding(&revision).unwrap();
    assert_eq!(
        workspace::load_workspace(&installed.archive)
            .unwrap()
            .revisions[&installed.receipt.product.revision]
            .sequence
            .as_str(),
        CUT
    );
    let source = w.revisions[&pin.revision].clone();
    let origin = MaterialOrigin::Derived(cap);
    let forward = Primer::new("Native-material F", "CCGTAGGA", "").unwrap();
    let reverse = Primer::new("Native-material R", "TGACGATC", "").unwrap();
    let revision = material_pcr::revision(
        &source,
        &origin,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    let previews = checked_roundtrip(
        AnalysisJob::MaterialPcrDocuments {
            molecule: source.clone(),
            origin: origin.clone(),
            forward: forward.clone(),
            reverse: reverse.clone(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        revision,
        169,
        159,
    );
    let AnalysisResult::MaterialPcrDocuments(previews) = previews.result else {
        panic!("material PCR")
    };
    assert_eq!(previews.len(), 1);
    assert_eq!(previews[0].molecule.sequence.as_str(), AMPLICON);
    let evidence: material_pcr::OriginEvidence =
        serde_json::from_str(&previews[0].molecule.metadata[material_pcr::ORIGIN_METADATA])
            .unwrap();
    assert_eq!(evidence.selection, origin.selection());
    assert!(!evidence.origin_incomplete);
    let request = pcr_outputs::Request {
        source,
        forward,
        reverse,
        selected_index: 0,
        expected_product_revision: previews[0].molecule.revision(),
        options: pcr_outputs::Options {
            name: "Explicit native-material PCR output".into(),
            save_product: true,
            save_primers: true,
        },
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let revision = material_pcr::outputs_revision(&origin, &request).unwrap();
    let outputs = checked_roundtrip(
        AnalysisJob::PrepareMaterialPcrOutputs { origin, request },
        revision.clone(),
        170,
        160,
    );
    let AnalysisResult::MaterialPcrOutputsPrepared(outputs) = outputs.result else {
        panic!("material PCR outputs")
    };
    outputs.validate_binding(&revision).unwrap();
    assert_eq!(outputs.product.unwrap().sequence.as_str(), AMPLICON);
    assert_eq!(outputs.primers.unwrap().records.len(), 2);
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
}
#[test]
fn current_material_origin_encoding_preserves_operation_and_exact_evidence() {
    assert_eq!(protocol::CANDIDATE_SCHEMA_VERSION, 29);
    assert_eq!(workspace::PLAN_SCHEMA_VERSION, 27);
    let (mut w, cap, pin) = digest_fixture();
    let (plan, target) = prepared_plan(&mut w, &cap, &pin);
    assert_eq!(tag(&plan.steps[0]), 27);
    let candidate = transfer::save_candidate(&transfer::Candidate::new(plan.clone())).unwrap();
    assert_eq!(
        u16::from_le_bytes([candidate[4], candidate[5]]),
        protocol::CANDIDATE_SCHEMA_VERSION
    );
    let mut old = candidate.clone();
    old[4..6].copy_from_slice(&20u16.to_le_bytes());
    assert_eq!(&old[6..], &candidate[6..]);
    assert!(transfer::load_candidate(&old).is_err());
    let report = verification::verify(&plan, &w.revisions[&target.revision]);
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
    let request = protocol::Request {
        id: 27,
        revision: target.revision.clone(),
        job: AnalysisJob::Verify {
            plan: plan.clone(),
            target: w.revisions[&target.revision].clone(),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert!(protocol::process_request(&wire(b"CGQ\0", 89, &request)).is_err());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert!(protocol::decode_response(&wire(b"CGR\0", 89, &response)).is_err());
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        target,
        vec![InputEventBinding {
            input_id: "fragment".into(),
            event: pin,
        }],
        &cloning_outputs::Options {
            name: "Recorded material flow".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let mut record =
        collagene_core::run_lineage::load_run(&saved.workspace, &saved.receipt.run.run_id).unwrap();
    assert_eq!(record.plan_schema, workspace::PLAN_SCHEMA_VERSION);
    // An unsupported internal plan representation is refused even inside a current run envelope.
    record.plan_schema = workspace::PLAN_SCHEMA_VERSION - 1;
    assert!(collagene_core::run_lineage::save_run_record(&record).is_err());
}
#[test]
fn source_step_reconstruction_is_refused_before_material_can_be_recovered() {
    let (mut w, cap, pin) = digest_fixture();
    let (mut plan, target) = prepared_plan(&mut w, &cap, &pin);
    plan.steps.insert(
        1,
        Operation::ReconstructMaterial {
            template: Source::Step(0),
            origin: cap,
        },
    );
    assert!(run_material::plan_preflight(&plan).is_err());
    let report = verification::verify(&plan, &w.revisions[&target.revision]);
    assert_eq!(report.status, Status::Invalid);
    assert!(report.products.is_empty());
}
#[test]
fn current_encoding_retains_documentary_producer_stamps_without_granting_authority() {
    let (_, cap, _) = digest_fixture();
    let native::SourceOrigin::Saved(MaterialOrigin::Run(origin)) =
        native::load_record(&cap.record_bytes)
            .unwrap()
            .request
            .origin
    else {
        panic!("the digest must retain its actual saved PCR run source")
    };
    let plan = transfer::load_candidate(&origin.candidate).unwrap().plan;
    assert!(
        plan.steps
            .iter()
            .all(|s| matches!(s, Operation::Pcr { .. }))
    );
    let report = verification::verify(&plan, &origin.target);
    assert!(
        report
            .engine_version
            .eq(collagene_core::verification::ENGINE_VERSION)
    );
    let response = protocol::Response {
        id: 1,
        revision: origin.target.revision(),
        result: AnalysisResult::Verification(report),
    };
    protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &response)).unwrap();
    let mut historical = response;
    let AnalysisResult::Verification(report) = &mut historical.result else {
        unreachable!()
    };
    report.engine_version = report
        .engine_version
        .replace(collagene_core::verification::ENGINE_VERSION, "/science-24");
    protocol::decode_response(&wire(
        b"CGR\0",
        protocol::WORKER_SCHEMA_VERSION,
        &historical,
    ))
    .unwrap();
    assert!(protocol::decode_response(&wire(b"CGR\0", 89, &historical)).is_err());
}
