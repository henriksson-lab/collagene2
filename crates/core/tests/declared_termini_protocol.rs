//! Explicit source assumptions cannot hide under old candidate/worker/archive versions.
use collagene_core::{
    Molecule, Topology,
    chemistry::Terminus,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    transfer::{self, Candidate},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
fn fixture() -> (DesignPlan, Molecule) {
    let source = Molecule::dna(
        "original intact linear reference",
        "ACGTTGCA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "independent terminal-state target",
        "ACGTTGCA",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "explicit original-source termini".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::DeclareSourceTermini {
            template: Source::Input("s".into()),
            left: Terminus::CLEAVED,
            right: Terminus::CLEAVED,
        }],
    };
    (plan, target)
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
#[test]
fn new_operation_requires_candidate_ten_and_worker_thirty_four_while_old_candidates_remain_readable()
 {
    let (plan, target) = fixture();
    let candidate = Candidate::new(plan.clone());
    let encoded = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    transfer::load_candidate(&encoded).unwrap();
    assert!(transfer::load_candidate(&wire(b"CGP\0", 9, &candidate)).is_err());
    let mut old = plan.clone();
    old.steps = vec![Operation::Overlap {
        fragments: vec![Source::Input("s".into())],
        overlaps: vec![],
        circular: false,
    }];
    assert!(transfer::load_candidate(&wire(b"CGP\0", 9, &Candidate::new(old))).is_err());
    let request = Request {
        id: 160,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 33, &request)).is_err());
    let response = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::Verification(report) = response.result else {
        panic!("verification expected")
    };
    assert_eq!(
        report.engine_version,
        collagene_core::verification::ENGINE_VERSION
    );
    assert!(
        protocol::decode_response(&wire(
            b"CGR\0",
            33,
            &Response {
                id: 160,
                revision: target.revision(),
                result: AnalysisResult::Verification(report)
            }
        ))
        .is_err()
    );
    assert_eq!(postcard::to_allocvec(&plan.steps[0]).unwrap()[0], 16);
}
#[test]
fn inline_and_binary_archives_reject_new_semantics_in_an_old_worker_envelope() {
    let (plan, target) = fixture();
    let mut frame = Workspace::from_documents(
        "declared source",
        vec![plan.inputs[0].molecule.clone(), target.clone()],
        1,
    )
    .unwrap();
    frame.candidates.push(Candidate::new(plan));
    let archive = workspace::save_workspace(&frame).unwrap();
    workspace::load_workspace(&archive).unwrap();
    for job in [
        AnalysisJob::SaveWorkspace { workspace: frame },
        AnalysisJob::LoadWorkspace {
            bytes: archive.clone(),
        },
    ] {
        let revision = match &job {
            AnalysisJob::LoadWorkspace { .. } => blake3::hash(&archive).to_hex().to_string(),
            _ => target.revision(),
        };
        let request = Request {
            id: 160,
            revision,
            job,
        };
        assert!(protocol::process_request(&wire(b"CGQ\0", 33, &request)).is_err());
    }
}
#[test]
fn archived_plan_eight_cannot_be_relabeled_seven() {
    use std::io::{Cursor, Read, Write};
    let (plan, target) = fixture();
    let mut frame = Workspace::from_documents(
        "declared source",
        vec![plan.inputs[0].molecule.clone(), target],
        1,
    )
    .unwrap();
    frame.candidates.push(Candidate::new(plan));
    let bytes = workspace::save_workspace(&frame).unwrap();
    let mut input = zip::ZipArchive::new(Cursor::new(bytes)).unwrap();
    let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for index in 0..input.len() {
        let mut entry = input.by_index(index).unwrap();
        let name = entry.name().to_string();
        let mut bytes = vec![];
        entry.read_to_end(&mut bytes).unwrap();
        if name == "manifest.json" {
            let mut m: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(m["plan_schema"], workspace::PLAN_SCHEMA_VERSION);
            m["plan_schema"] = 7.into();
            bytes = serde_json::to_vec(&m).unwrap();
        }
        output
            .start_file(name, zip::write::SimpleFileOptions::default())
            .unwrap();
        output.write_all(&bytes).unwrap();
    }
    assert!(workspace::load_workspace(&output.finish().unwrap().into_inner()).is_err());
}

#[test]
fn declared_preparation_cannot_hide_in_an_old_existing_recipe_job() {
    use collagene_core::assembly_recipe::{FragmentSpec, Preparation, Recipe};
    let (plan, target) = fixture();
    let recipe = Recipe {
        name: "declared overlap source".into(),
        fragments: vec![FragmentSpec {
            generation: None,
            role: "arm".into(),
            input_id: "s".into(),
            preparation: Preparation::DeclaredDuplex {
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED,
            },
            reverse: false,
        }],
        overlaps: vec![],
        circular: false,
        options: Default::default(),
    };
    let revision =
        collagene_core::assembly_recipe::preview_revision(&plan.inputs, &target, &recipe).unwrap();
    let request = Request {
        id: 160,
        revision,
        job: AnalysisJob::PreviewOverlapRecipe {
            inputs: plan.inputs,
            target,
            recipe,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 33, &request)).is_err());
    let response = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    protocol::decode_response(&response).unwrap();
}
