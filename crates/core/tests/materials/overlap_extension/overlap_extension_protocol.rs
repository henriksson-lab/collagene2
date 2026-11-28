//! Real reciprocal extension stays physical across binary plans and closed origins.
use collagene_core::{
    cloning::overlap::extension::Options,
    documents::{
        transfer::{self, Candidate},
        workspace::Workspace,
    },
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::Selection,
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EXPECTED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
fn fixture() -> (DesignPlan, Molecule) {
    let original = Molecule::dna("Original unchanged stock", ORIGINAL, Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent GGAC to TT requirement",
        EXPECTED,
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Two actual mutagenic PCRs and reciprocal extension".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: original.revision(),
            molecule: original,
        }],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("original".into()),
                forward: Primer::new("Left outer", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Left mutagenic", "TACGTAGC", "GGTCTGAAAA").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("original".into()),
                forward: Primer::new("Right mutagenic", "TTCAGACC", "GCTACGTATT").unwrap(),
                reverse: Primer::new("Right outer", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::OverlapExtension {
                left: OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                },
                right: OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                },
                options: Options {
                    overlap: 18,
                    min_paired_bases: 8,
                    ..Default::default()
                },
            },
            Operation::Pcr {
                template: Source::Step(2),
                forward: Primer::new("Outer forward", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Outer reverse", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        ],
    };
    (plan, target)
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&payload).as_bytes());
    out.extend(payload);
    out
}
fn run(job: AnalysisJob, revision: &str) -> AnalysisResult {
    let request = Request {
        id: 81,
        revision: revision.into(),
        job,
    };
    let reply = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    assert_eq!(
        &reply[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    protocol::response::Response::decode(&reply).unwrap().result
}
#[test]
fn real_binary_pcr_extension_outer_pcr_matches_independent_target_and_exports() {
    let (plan, target) = fixture();
    let original = serde_json::to_value(&plan.inputs[0].molecule).unwrap();
    let AnalysisResult::Verification(report) = run(
        AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
        &target.revision(),
    ) else {
        panic!("verification");
    };
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
    assert_eq!(report.reconstructed_count(), 4);
    let extension = report
        .physical_processing_reports
        .iter()
        .find(|step| step.step_index == 2)
        .expect("retained reciprocal-extension evidence");
    let collagene_core::verification::evidence::PhysicalProcessingEvidence::OverlapExtension(
        geometry,
    ) = &extension.evidence
    else {
        panic!("extension geometry");
    };
    geometry.validate_documentary().unwrap();
    assert_eq!(geometry.options.overlap, 18);
    assert_eq!(geometry.homology.len(), 18);
    assert_eq!(geometry.fill.duplex.top.as_str(), EXPECTED);
    assert!(geometry.terminal_pairings.complete);
    assert!(geometry.terminal_pairings.unique);
    assert_eq!(
        report
            .product(0)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        "AGTCCGATGCTACGTATTTTCAGACC"
    );
    assert_eq!(
        report
            .product(1)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        "GCTACGTATTTTCAGACCTGTAGATCGTCA"
    );
    assert_eq!(
        report
            .product(2)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        EXPECTED
    );
    assert_eq!(
        report
            .product(3)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        EXPECTED
    );
    assert_eq!(
        serde_json::to_value(&plan.inputs[0].molecule).unwrap(),
        original
    );
    let AnalysisResult::ProtocolHtml(html) = run(
        AnalysisJob::ExportProtocolHtml {
            plan: plan.clone(),
            target: target.clone(),
        },
        &target.revision(),
    ) else {
        panic!("HTML");
    };
    assert!(html.to_ascii_lowercase().contains("extension"));
    let AnalysisResult::ProtocolPdf(pdf) = run(
        AnalysisJob::ExportProtocolPdf {
            plan: plan.clone(),
            target: target.clone(),
        },
        &target.revision(),
    ) else {
        panic!("PDF");
    };
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    let mut wrong = target;
    wrong.sequence = plan.inputs[0].molecule.sequence.clone();
    let mut wrong_plan = plan;
    wrong_plan.target_revision = wrong.revision();
    let AnalysisResult::Verification(report) = run(
        AnalysisJob::Verify {
            plan: wrong_plan,
            target: wrong.clone(),
        },
        &wrong.revision(),
    ) else {
        panic!("verification");
    };
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(
        report
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        EXPECTED
    );
}
#[test]
fn extension_candidates_reject_unsupported_internal_versions() {
    let (plan, target) = fixture();
    let encoded = postcard::to_allocvec(&plan.steps[2]).unwrap();
    let (tag, _) = postcard::take_from_bytes::<u32>(&encoded).unwrap();
    assert_eq!(tag, 22);
    let candidate = Candidate::new(plan.clone());
    let bytes = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&wire(b"CGP\0", 18, &candidate)).is_err());
    assert!(transfer::load_candidate(&bytes).is_ok());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 17, &candidate)).is_err());
    let request = Request {
        id: 1,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: plan.clone(),
            target,
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 80, &request)).is_err());
    let mut plain = plan;
    plain.steps.truncate(1);
    let candidate = Candidate::new(plain);
    assert!(transfer::load_candidate(&wire(b"CGP\0", 1, &candidate)).is_err());
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn closed_extension_origin_cannot_hide_new_semantics_in_old_worker_envelope() {
    let (plan, target) = fixture();
    let workspace = Workspace::from_documents(
        "Actual extension origin",
        vec![plan.inputs[0].molecule.clone(), target.clone()],
        1,
    )
    .unwrap();
    let commit = cloning_outputs::finalize(
        &workspace,
        &plan,
        event(&workspace, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(&workspace, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved extension".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let save = Request {
        id: 3,
        revision: commit
            .workspace
            .head(&commit.workspace.active_document)
            .unwrap()
            .revision(),
        job: AnalysisJob::SaveWorkspace {
            workspace: commit.workspace.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 80, &save)).is_err());
    let saved_reply =
        protocol::execution::process_request(&protocol::request::Request::encode(&save).unwrap())
            .unwrap();
    let AnalysisResult::WorkspaceSaved(archive) =
        protocol::response::Response::decode(&saved_reply)
            .unwrap()
            .result
    else {
        panic!("archive");
    };
    let loaded = collagene_core::documents::workspace::load_workspace(&archive).unwrap();
    assert_eq!(loaded.candidates.len(), commit.workspace.candidates.len());
    let load = Request {
        id: 4,
        revision: blake3::hash(&archive).to_hex().to_string(),
        job: AnalysisJob::LoadWorkspace { bytes: archive },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 80, &load)).is_err());
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&load).unwrap())
            .is_ok()
    );
    let record = run_lineage::load_run(&commit.workspace, &commit.receipt.run.run_id).unwrap();
    assert_eq!(
        record.plan_schema,
        collagene_core::documents::workspace::PLAN_SCHEMA_VERSION
    );
    assert!(run_lineage::save_run_record(&record).is_ok());
    let mut downgraded = record;
    downgraded.plan_schema = 15;
    assert!(run_lineage::save_run_record(&downgraded).is_err());
    let selection = Selection {
        run_id: commit.receipt.run.run_id,
        step_index: 3,
        event: commit.receipt.construct.clone(),
    };
    let origin = collagene_core::material::capsule::capsule(&commit.workspace, &selection).unwrap();
    let saved = commit.workspace.revisions[&selection.event.revision].clone();
    let continuation = DesignPlan {
        name: "Freshly reconstructed saved extension".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "saved".into(),
            revision: saved.revision(),
            molecule: saved,
        }],
        steps: vec![Operation::ReconstructRunProduct {
            template: Source::Input("saved".into()),
            origin,
        }],
    };
    let request = Request {
        id: 2,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: continuation,
            target,
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 80, &request)).is_err());
    let reply = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    let response = protocol::response::Response::decode(&reply).unwrap();
    let AnalysisResult::Verification(report) = &response.result else {
        panic!("verification");
    };
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        EXPECTED
    );
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 80, &response)).is_err());
}

#[test]
fn competing_terminal_pairings_remain_incomplete_after_outer_pcr() {
    let left = Molecule::dna("Left stock", "AGTCCGATGCTACGTAAA", Topology::Linear).unwrap();
    let right = Molecule::dna("Right stock", "AAATTCAGACCTGTAGATCGTCA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent selected candidate",
        "AGTCCGATGCTACGTAAATTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "An outer PCR does not resolve competing terminal pairings".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left.revision(),
                molecule: left,
            },
            PlanInput {
                id: "right".into(),
                revision: right.revision(),
                molecule: right,
            },
        ],
        steps: vec![
            Operation::DeclareSourceTermini {
                template: Source::Input("left".into()),
                left: collagene_core::chemistry::Terminus::PCR,
                right: collagene_core::chemistry::Terminus::PCR,
            },
            Operation::DeclareSourceTermini {
                template: Source::Input("right".into()),
                left: collagene_core::chemistry::Terminus::PCR,
                right: collagene_core::chemistry::Terminus::PCR,
            },
            Operation::OverlapExtension {
                left: OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                },
                right: OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                },
                options: Options {
                    overlap: 3,
                    min_paired_bases: 2,
                    ..Default::default()
                },
            },
            Operation::Pcr {
                template: Source::Step(2),
                forward: Primer::new("Outer forward", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Outer reverse", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        ],
    };
    let report = collagene_core::verification::engine::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.reconstructed_count(), 4, "{:?}", report.findings);
    assert_eq!(
        report
            .product(3)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        target.sequence.as_str()
    );
    let collagene_core::verification::evidence::PhysicalProcessingEvidence::OverlapExtension(
        geometry,
    ) = &report.physical_processing_reports[0].evidence
    else {
        panic!("extension geometry");
    };
    assert_eq!(geometry.terminal_pairings.overlaps, vec![2, 3]);
    assert!(geometry.terminal_pairings.complete);
    assert!(!geometry.terminal_pairings.unique);
}
