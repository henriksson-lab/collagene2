//! Current binary admission of declared vector graphs and IVT cut alternatives.
use collagene_core::{
    chemistry::tracked_duplex::GenerationAssumptions,
    cloning::restriction::{catalog, ivt_linearization},
    construct_qc::{Verdict, vector_class},
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
    },
    rna_chemistry::{Provenance, ivt},
};
use collagene_model::{Molecule, Topology};

fn raw<T: serde::Serialize>(magic: &[u8; 4], value: &T) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    bytes
}
fn execute(request: &Request) -> Response {
    let encoded = request.encode().unwrap();
    assert_eq!(&encoded[4..6], &WORKER_SCHEMA_VERSION.to_le_bytes());
    Response::decode(&protocol::execution::process_request(&encoded).unwrap()).unwrap()
}
fn picker(source: &Molecule) -> ivt_linearization::Request {
    ivt_linearization::Request {
        ivt: ivt::Request {
            source_revision: source.revision(),
            profile: ivt::profiles().unwrap().remove(0).reference,
            transcript: ivt::Transcript {
                start: 17,
                length: 7,
                coding_strand: ivt::CodingStrand::Forward,
            },
            cap: ivt::CapChoice::NotSpecified,
            initiator: ivt::InitiatorRequirement::NoAdditionalRequirement,
            linearization: None,
            poly_a: None,
            provenance: Provenance {
                statement: "Literal source interval, no performed cut".into(),
                reference: None,
            },
        },
        source: ivt_linearization::SourceInput::DeclaredOriginal {
            generation: GenerationAssumptions::default(),
        },
        selection: catalog::select(
            &catalog::builtin_catalog().unwrap(),
            &["smai".into(), "ecori".into()],
        )
        .unwrap(),
        conditions: None,
        policy: ivt_linearization::Policy {
            allowed_end_kinds: vec![ivt::EndKind::Blunt],
            require_known_conditions: false,
            require_resolved_ivt_checks: false,
            require_complete_origin: true,
        },
        limits: ivt_linearization::Limits {
            max_comparisons: 1,
            max_sites_per_enzyme: 64,
            max_report_bytes: ivt_linearization::MAX_REPORT_BYTES,
        },
    }
}

#[test]
fn vector_discovery_and_empty_declared_graph_roundtrip_without_invented_role_evidence() {
    let discovery = Request {
        id: 11,
        revision: vector_class::catalogue_revision().unwrap(),
        job: AnalysisJob::ListVectorClassProfiles,
    };
    let response = execute(&discovery);
    let AnalysisResult::VectorClassProfiles(profiles) = response.result else {
        panic!("profile discovery")
    };
    assert_eq!(profiles, vector_class::catalogue().unwrap());
    let source = Molecule::dna(
        "No annotated role premises",
        "ACGTACGTACGT",
        Topology::Linear,
    )
    .unwrap();
    let declaration = vector_class::Request {
        source_revision: source.revision(),
        profile: profiles[0].reference.clone(),
        annotations: vec![],
        references: vec![],
        bindings: vec![],
        boundaries: None,
        spacing: vec![],
        literal_parts: vec![],
        literal_part_limits: Default::default(),
    };
    let request = Request {
        id: 12,
        revision: vector_class::revision(&source, &declaration).unwrap(),
        job: AnalysisJob::AssessVectorClass {
            molecule: source.clone(),
            request: declaration.clone(),
        },
    };
    let encoded = request.encode().unwrap();
    let mut session = protocol::session::AnalysisSession::default();
    for response in [
        execute(&request),
        Response::decode(&session.process_request(&encoded).unwrap()).unwrap(),
    ] {
        assert_eq!(response.revision, request.revision);
        let AnalysisResult::VectorClass(report) = &response.result else {
            panic!("vector report")
        };
        report.validate_for(&source, &declaration).unwrap();
        assert!(report.bindings.is_empty());
        assert!(
            report
                .findings
                .iter()
                .any(|f| f.verdict == Verdict::Indeterminate)
        );
        let mut forged = response;
        let AnalysisResult::VectorClass(report) = &mut forged.result else {
            unreachable!()
        };
        report.source_revision = "0".repeat(64);
        assert!(Response::decode(&raw(b"CGR\0", &forged)).is_err());
    }
    let mut stale = request.clone();
    stale.revision = "0".repeat(64);
    assert!(stale.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &stale)).is_err());
    let mut profile = request;
    let AnalysisJob::AssessVectorClass { request, .. } = &mut profile.job else {
        unreachable!()
    };
    request.profile.identity = "0".repeat(64);
    assert!(profile.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &profile)).is_err());
    let mut forged = Response {
        id: 11,
        revision: discovery.revision,
        result: AnalysisResult::VectorClassProfiles(profiles),
    };
    let AnalysisResult::VectorClassProfiles(profiles) = &mut forged.result else {
        unreachable!()
    };
    profiles.clear();
    assert!(Response::decode(&raw(b"CGR\0", &forged)).is_err());
}

#[test]
fn ivt_picker_binary_keeps_every_limited_alternative_and_refuses_stale_or_forged_evidence() {
    let source = Molecule::dna(
        "Bounded IVT alternatives",
        "TAATACGACTCACTATAAGTTCCCGGGTT",
        Topology::Linear,
    )
    .unwrap();
    let declaration = picker(&source);
    let request = Request {
        id: 13,
        revision: ivt_linearization::revision(&source, &declaration).unwrap(),
        job: AnalysisJob::AssessIvtLinearization {
            molecule: source.clone(),
            request: declaration.clone(),
        },
    };
    let encoded = request.encode().unwrap();
    let mut session = protocol::session::AnalysisSession::default();
    for response in [
        execute(&request),
        Response::decode(&session.process_request(&encoded).unwrap()).unwrap(),
    ] {
        let AnalysisResult::IvtLinearization(report) = &response.result else {
            panic!("picker report")
        };
        report.validate_for(&source, &declaration).unwrap();
        assert_eq!(report.candidates.len(), 2);
        assert_eq!(report.comparisons, 1);
        for candidate in &report.candidates {
            assert_eq!(
                candidate.state,
                ivt_linearization::CandidateState::Unresolved
            );
            assert!(!candidate.cleavage.complete);
            assert!(
                !candidate
                    .issues
                    .iter()
                    .any(|i| i.kind == ivt_linearization::IssueKind::NoPermittedCut)
            );
        }
        let mut forged = response;
        let AnalysisResult::IvtLinearization(report) = &mut forged.result else {
            unreachable!()
        };
        report.candidates[0].state = ivt_linearization::CandidateState::EligibleUnderDeclarations;
        assert!(Response::decode(&raw(b"CGR\0", &forged)).is_err());
    }
    let mut stale = request.clone();
    stale.revision = "0".repeat(64);
    assert!(stale.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &stale)).is_err());
    let mut stale_catalog = request.clone();
    let AnalysisJob::AssessIvtLinearization { request, .. } = &mut stale_catalog.job else {
        unreachable!()
    };
    request.selection.catalog_content_hash = "0".repeat(64);
    assert!(stale_catalog.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &stale_catalog)).is_err());
    let mut previous = encoded;
    previous[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert!(protocol::execution::process_request(&previous).is_err());
}
