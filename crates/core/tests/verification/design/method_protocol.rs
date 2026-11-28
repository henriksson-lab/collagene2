use collagene_core::{
    analysis::Primer,
    overlap_profiles::ProfileOptions,
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    tailing::TailingBase,
    transfer::Candidate,
    verification::*,
};
use collagene_model::{Molecule, Topology};
fn envelope(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    b
}
fn source() -> Molecule {
    Molecule::dna("source", "GATATCAAAA", Topology::Circular).unwrap()
}
fn plan(m: &Molecule) -> DesignPlan {
    DesignPlan {
        name: "tail candidate".into(),
        target_revision: m.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: m.revision(),
            molecule: m.clone(),
        }],
        steps: vec![
            Operation::Digest {
                template: Source::Input("s".into()),
                enzymes: vec!["EcoRV".into()],
                fragment_index: 0,
            },
            Operation::Tail3Prime {
                template: Source::Step(0),
                base: TailingBase::A,
            },
        ],
    }
}
#[test]
fn schema_eleven_tail_and_profile_jobs_bind_every_source_and_option() {
    let m = source();
    let req = Request {
        id: 5,
        revision: m.revision(),
        job: AnalysisJob::TailPreview {
            molecule: m.clone(),
            recipe: TailSubstrate::Digest {
                enzymes: vec!["EcoRV".into()],
                fragment_index: 0,
                blunt_before: false,
            },
            base: TailingBase::A,
        },
    };
    let bytes = protocol::request::Request::encode(&req).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(&bytes).unwrap(),
    )
    .unwrap();
    let AnalysisResult::TailPreview(p) = response.result else {
        panic!()
    };
    assert_eq!(p.source_revision, m.revision());
    assert_eq!(p.tailing.duplex.top.as_str(), "ATCAAAAGATA");
    for v in 1..11 {
        assert!(protocol::execution::process_request(&envelope(b"CGQ\0", v, &req)).is_err());
    }
    let mols = vec![
        Molecule::dna("a", "AAAACC", Topology::Linear).unwrap(),
        Molecule::dna("b", "AACCGG", Topology::Linear).unwrap(),
    ];
    let reverse = vec![false, true];
    let options = ProfileOptions {
        scan_repeated_overlaps: false,
        ..Default::default()
    };
    let binding =
        protocol::bindings::overlap_profile_revision(&mols, &reverse, &[2], false, &options)
            .unwrap();
    let req = Request {
        id: 6,
        revision: binding.clone(),
        job: AnalysisJob::CheckOverlapProfile {
            molecules: mols,
            reverse,
            overlaps: vec![2],
            circular: false,
            options,
        },
    };
    let r = protocol::response::Response::decode(
        &protocol::execution::process_request(&protocol::request::Request::encode(&req).unwrap())
            .unwrap(),
    )
    .unwrap();
    assert_eq!(r.revision, binding);
    let AnalysisResult::OverlapProfile(report) = r.result else {
        panic!()
    };
    assert!(report.geometry_valid);
    assert_eq!(report.orientations, vec![false, true]);
    let mut changed = req;
    changed.revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&changed).unwrap()
        )
        .is_err()
    );
}
#[test]
fn newer_operations_cannot_hide_in_old_candidates_verify_pdf_or_workspace_messages() {
    let m = source();
    let plan = plan(&m);
    let c = Candidate::new(plan.clone());
    let bytes = collagene_core::transfer::save_candidate(&c).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    for v in 1..5 {
        assert!(collagene_core::transfer::load_candidate(&envelope(b"CGP\0", v, &c)).is_err());
    }
    let mut workspace =
        collagene_core::workspace::Workspace::from_documents("bundle", vec![m.clone()], 0).unwrap();
    workspace.candidates.push(c);
    for job in [
        AnalysisJob::Verify {
            plan: plan.clone(),
            target: m.clone(),
        },
        AnalysisJob::ExportProtocolPdf {
            plan,
            target: m.clone(),
        },
        AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
    ] {
        let req = Request {
            id: 7,
            revision: m.revision(),
            job,
        };
        for v in 1..11 {
            assert!(
                protocol::execution::process_request(&envelope(b"CGQ\0", v, &req)).is_err(),
                "accepted worker version{v}"
            );
        }
    }
    let response = Response {
        id: 7,
        revision: m.revision(),
        result: AnalysisResult::WorkspaceLoaded(workspace),
    };
    for v in 1..11 {
        assert!(protocol::response::Response::decode(&envelope(b"CGR\0", v, &response)).is_err());
    }
}
#[test]
fn tail_and_profile_responses_reject_stale_or_forged_evidence() {
    let m = source();
    let preview = simulate_tailing(
        &m,
        &TailSubstrate::Digest {
            enzymes: vec!["EcoRV".into()],
            fragment_index: 0,
            blunt_before: false,
        },
        TailingBase::A,
    )
    .unwrap();
    let mut response = Response {
        id: 8,
        revision: m.revision(),
        result: AnalysisResult::TailPreview(preview),
    };
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_ok()
    );
    response.revision = "0".repeat(64);
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_err()
    );
    if let AnalysisResult::TailPreview(ref mut p) = response.result {
        response.revision = m.revision();
        p.tailing.additions.clear();
    }
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_err()
    );
}
#[test]
fn pcr_recipe_unknown_fields_never_become_end_chemistry_declarations() {
    let recipe = TailSubstrate::ExactPcr {
        forward: Primer::new("f", "ACGA", "").unwrap(),
        reverse: Primer::new("r", "GTTA", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let mut value = serde_json::to_value(recipe).unwrap();
    value["ExactPcr"]["phosphorylated"] = serde_json::json!(true);
    assert!(serde_json::from_value::<TailSubstrate>(value).is_err());
}
