//! Tracked operations never hide in historical plan/workspace/worker schemas.
use collagene_core::{
    Molecule, Topology,
    methylation::MethylationState,
    modified_duplex::GenerationAssumptions,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    transfer::{self, Candidate},
    verification::*,
    workspace::Workspace,
};
fn fixture() -> (Molecule, DesignPlan) {
    let source = Molecule::dna("source", "AAAAGAATTCCCC", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AAAAG", Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "checked digestion".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("s".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::Unmodified,
                    dcm: MethylationState::Unmodified,
                    cpg: MethylationState::Unmodified,
                },
            },
            Operation::CheckedDigest {
                template: Source::Step(0),
                enzymes: vec!["EcoRI".into()],
                profile_revision: collagene_core::modified_duplex::profile_revision(),
                fragment_index: 0,
            },
        ],
    };
    (target, plan)
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
fn appended_operations_require_candidate_six_and_worker_seventeen_in_nested_plans() {
    let (target, plan) = fixture();
    let candidate = Candidate::new(plan.clone());
    let bytes = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(&bytes[4..6], &7u16.to_le_bytes());
    transfer::load_candidate(&bytes).unwrap();
    for version in 1..6 {
        assert!(transfer::load_candidate(&wire(b"CGP\0", version, &candidate)).is_err());
    }
    for job in [
        AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
        AnalysisJob::ExportProtocolPdf {
            plan: plan.clone(),
            target: target.clone(),
        },
    ] {
        let request = Request {
            id: 501,
            revision: target.revision(),
            job,
        };
        for version in 1..17 {
            assert!(protocol::process_request(&wire(b"CGQ\0", version, &request)).is_err());
        }
        let bytes = protocol::encode_request(&request).unwrap();
        assert_eq!(&bytes[4..6], &23u16.to_le_bytes());
        protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    }
}
#[test]
fn inline_and_binary_workspace_plan_graphs_cannot_downcast_tracked_operations() {
    let (target, plan) = fixture();
    let mut workspace = Workspace::from_documents(
        "tracked",
        vec![plan.inputs[0].molecule.clone(), target.clone()],
        1,
    )
    .unwrap();
    workspace.candidates.push(Candidate::new(plan));
    let archive = collagene_core::workspace::save_workspace(&workspace).unwrap();
    let jobs = [
        AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
        AnalysisJob::LoadWorkspace {
            bytes: archive.clone(),
        },
    ];
    for job in jobs {
        let revision = match &job {
            AnalysisJob::LoadWorkspace { bytes } => blake3::hash(bytes).to_hex().to_string(),
            _ => target.revision(),
        };
        let request = Request {
            id: 502,
            revision,
            job,
        };
        assert!(protocol::process_request(&wire(b"CGQ\0", 16, &request)).is_err());
    }
    for result in [
        AnalysisResult::WorkspaceLoaded(workspace),
        AnalysisResult::WorkspaceSaved(archive),
    ] {
        let response = Response {
            id: 502,
            revision: target.revision(),
            result,
        };
        assert!(protocol::decode_response(&wire(b"CGR\0", 16, &response)).is_err());
        protocol::decode_response(&wire(b"CGR\0", 17, &response)).unwrap();
    }
}
#[test]
fn science_nine_reports_require_new_envelopes_while_eight_and_seven_stay_historical() {
    let (target, plan) = fixture();
    let report = verify(&plan, &target);
    assert!(report.engine_version.ends_with("/science-9"));
    let response = Response {
        id: 503,
        revision: target.revision(),
        result: AnalysisResult::Verification(report.clone()),
    };
    assert!(protocol::decode_response(&wire(b"CGR\0", 16, &response)).is_err());
    assert!(protocol::decode_response(&wire(b"CGR\0", 21, &response)).is_err());
    protocol::decode_response(&wire(b"CGR\0", 22, &response)).unwrap();
    let mut historical_eight = report.clone();
    historical_eight.engine_version = "0.1.0/science-8".into();
    let historical = Response {
        id: 503,
        revision: target.revision(),
        result: AnalysisResult::Verification(historical_eight),
    };
    protocol::decode_response(&wire(b"CGR\0", 17, &historical)).unwrap();

    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    let response = Response {
        id: 503,
        revision: target.revision(),
        result: AnalysisResult::ProtocolPdf(pdf),
    };
    assert!(protocol::decode_response(&wire(b"CGR\0", 16, &response)).is_err());
    let mut old = report;
    old.engine_version = "0.1.0/science-7".into();
    let response = Response {
        id: 503,
        revision: target.revision(),
        result: AnalysisResult::Verification(old),
    };
    let decoded = protocol::decode_response(&wire(b"CGR\0", 16, &response)).unwrap();
    let AnalysisResult::Verification(old) = decoded.result else {
        panic!()
    };
    assert_eq!(old.engine_version, "0.1.0/science-7");
    assert_ne!(old.engine_version, ENGINE_VERSION);
}
#[test]
fn native_archive_plan_four_floor_is_checked_even_when_its_candidate_hashes_are_valid() {
    use std::io::{Cursor, Read, Write};
    let (target, plan) = fixture();
    let mut workspace =
        Workspace::from_documents("tracked", vec![plan.inputs[0].molecule.clone(), target], 1)
            .unwrap();
    workspace.candidates.push(Candidate::new(plan));
    let archive = collagene_core::workspace::save_workspace(&workspace).unwrap();
    for schema in 1..4 {
        let mut input = zip::ZipArchive::new(Cursor::new(&archive)).unwrap();
        let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for index in 0..input.len() {
            let mut entry = input.by_index(index).unwrap();
            let name = entry.name().to_owned();
            let mut data = Vec::new();
            entry.read_to_end(&mut data).unwrap();
            if name == "manifest.json" {
                let mut manifest: serde_json::Value = serde_json::from_slice(&data).unwrap();
                assert_eq!(manifest["plan_schema"], 5);
                manifest["plan_schema"] = serde_json::json!(schema);
                data = serde_json::to_vec(&manifest).unwrap();
            }
            output
                .start_file(name, zip::write::SimpleFileOptions::default())
                .unwrap();
            output.write_all(&data).unwrap();
        }
        let mislabeled = output.finish().unwrap().into_inner();
        assert!(collagene_core::workspace::load_workspace(&mislabeled).is_err());
    }
}
#[test]
fn retained_runs_and_their_workspaces_cannot_hide_new_semantics_under_old_outer_variants() {
    use collagene_core::run_lineage::*;
    let (target, plan) = fixture();
    let workspace = Workspace::from_documents(
        "retained",
        vec![plan.inputs[0].molecule.clone(), target.clone()],
        1,
    )
    .unwrap();
    let binding = |index: usize| EventBinding {
        document_id: workspace.documents[index].id.clone(),
        event_id: workspace.documents[index].head_event.clone().unwrap(),
        revision: workspace.documents[index].head.clone(),
    };
    let target_binding = binding(1);
    let input_bindings = vec![InputEventBinding {
        input_id: "s".into(),
        event: binding(0),
    }];
    let request = Request {
        id: 504,
        revision: workspace_fingerprint(&workspace).unwrap(),
        job: AnalysisJob::RecordCloningRun {
            workspace: workspace.clone(),
            plan: plan.clone(),
            target_binding: target_binding.clone(),
            input_bindings: input_bindings.clone(),
            allow_incomplete: false,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 16, &request)).is_err());
    let commit = record_run(&workspace, &plan, target_binding, input_bindings, false).unwrap();
    let record =
        load_run_record(&commit.workspace.attachments[&commit.receipt.attachment_hash].bytes)
            .unwrap();
    assert_eq!(record.plan_schema, 5);
    let mut mislabeled = serde_json::to_value(&record).unwrap();
    mislabeled["plan_schema"] = serde_json::json!(3);
    assert!(load_run_record(&serde_json::to_vec(&mislabeled).unwrap()).is_err());
    let replay = replay_run(&commit.workspace, &commit.receipt.run_id).unwrap();
    assert!(replay.products_match);
    let mut historical_reply = replay.clone();
    historical_reply.current_report.engine_version = "0.1.0/science-7".into();
    historical_reply.engine_matches =
        historical_reply.recorded_engine == historical_reply.current_report.engine_version;
    historical_reply.validate().unwrap();
    assert!(
        historical_reply
            .validate_for_workspace(&commit.workspace)
            .is_err()
    );
    let response = Response {
        id: 504,
        revision: workspace_fingerprint(&workspace).unwrap(),
        result: AnalysisResult::RunCommit(commit.clone()),
    };
    assert!(protocol::decode_response(&wire(b"CGR\0", 16, &response)).is_err());
    let response = Response {
        id: 504,
        revision: workspace_fingerprint(&commit.workspace).unwrap(),
        result: AnalysisResult::RunReplay(replay),
    };
    assert!(protocol::decode_response(&wire(b"CGR\0", 16, &response)).is_err());
    for job in [
        AnalysisJob::SaveWorkspace {
            workspace: commit.workspace.clone(),
        },
        AnalysisJob::ReplayCloningRun {
            workspace: commit.workspace.clone(),
            run_id: commit.receipt.run_id,
        },
    ] {
        let request = Request {
            id: 504,
            revision: workspace_fingerprint(&commit.workspace).unwrap(),
            job,
        };
        assert!(protocol::process_request(&wire(b"CGQ\0", 16, &request)).is_err());
    }
}
#[test]
fn genuine_schema_three_run_identity_remains_readable_with_its_recorded_engine() {
    use collagene_core::run_lineage::*;
    let source = Molecule::dna("old source", "AACCGG", Topology::Linear).unwrap();
    let target = source.clone();
    let workspace = Workspace::from_documents("history", vec![source.clone()], 0).unwrap();
    let event = EventBinding {
        document_id: workspace.documents[0].id.clone(),
        event_id: workspace.documents[0].head_event.clone().unwrap(),
        revision: source.revision(),
    };
    let plan = DesignPlan {
        name: "old model".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("s".into())],
            overlaps: vec![],
            circular: false,
        }],
    };
    let commit = record_run(
        &workspace,
        &plan,
        event.clone(),
        vec![InputEventBinding {
            input_id: "s".into(),
            event,
        }],
        false,
    )
    .unwrap();
    let mut record =
        load_run_record(&commit.workspace.attachments[&commit.receipt.attachment_hash].bytes)
            .unwrap();
    record.plan_schema = 3;
    record.assessment.engine_version = "0.1.0/science-7".into();
    record.run_id = blake3::hash(
        &serde_json::to_vec(&(
            "cloning-run-v1",
            record.schema_version,
            record.plan_schema,
            &record.plan_revision,
            &record.plan,
            &record.target,
            &record.inputs,
            record
                .products
                .iter()
                .map(|p| (&p.step_index, &p.revision))
                .collect::<Vec<_>>(),
            &record.assessment,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let historical = load_run_record(&serde_json::to_vec(&record).unwrap()).unwrap();
    assert_eq!(historical.plan_schema, 3);
    assert_eq!(historical.assessment.engine_version, "0.1.0/science-7");
    assert_ne!(historical.assessment.engine_version, ENGINE_VERSION);
}
