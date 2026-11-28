//! Current binary boundary for source-bound documentary IVT template checks.
use collagene_core::{
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
    },
    rna_chemistry::{Provenance, ivt},
};
use collagene_model::{Molecule, Topology};

fn inputs() -> (Molecule, ivt::Request) {
    let source = Molecule::dna(
        "Literal IVT template",
        "TAATACGACTCACTATAAGTTAAAA",
        Topology::Linear,
    )
    .unwrap();
    let declaration = ivt::Request {
        source_revision: source.revision(),
        profile: ivt::profiles().unwrap().remove(0).reference,
        transcript: ivt::Transcript {
            start: 17,
            length: 8,
            coding_strand: ivt::CodingStrand::Forward,
        },
        cap: ivt::CapChoice::CleanCapAgN7113,
        initiator: ivt::InitiatorRequirement::NoAdditionalRequirement,
        linearization: Some(ivt::Linearization {
            coding_cut_offset: 8,
            template_cut_offset: 8,
            provenance: Provenance {
                statement: "Supplied geometry, no physical cutting inferred".into(),
                reference: None,
            },
        }),
        poly_a: Some(ivt::PolyARegion {
            start_offset: 4,
            end_offset: 8,
        }),
        provenance: Provenance {
            statement: "Explicit original coding-strand interval".into(),
            reference: None,
        },
    };
    (source, declaration)
}
fn raw<T: serde::Serialize>(magic: &[u8; 4], version: u16, value: &T) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&version.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    bytes
}

#[test]
fn ivt_discovery_and_assessment_roundtrip_share_the_exact_profile_and_source_contract() {
    let discovery = Request {
        id: 7,
        revision: ivt::catalogue_revision().unwrap(),
        job: AnalysisJob::ListIvtProfiles,
    };
    let mut session = protocol::session::AnalysisSession::default();
    let bytes = discovery.encode().unwrap();
    for response in [
        protocol::execution::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        let reply = Response::decode(&response).unwrap();
        assert_eq!(
            (reply.id, &reply.revision),
            (discovery.id, &discovery.revision)
        );
        let AnalysisResult::IvtProfiles(profiles) = reply.result else {
            panic!("IVT discovery result");
        };
        assert_eq!(profiles, ivt::profiles().unwrap());
    }
    let (source, declaration) = inputs();
    let before = source.clone();
    let request = Request {
        id: 8,
        revision: ivt::revision(&source, &declaration).unwrap(),
        job: AnalysisJob::AssessIvtTemplate {
            molecule: source.clone(),
            request: declaration.clone(),
        },
    };
    let bytes = request.encode().unwrap();
    assert_eq!(&bytes[4..6], &WORKER_SCHEMA_VERSION.to_le_bytes());
    for response in [
        protocol::execution::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        assert_eq!(&response[4..6], &WORKER_SCHEMA_VERSION.to_le_bytes());
        let reply = Response::decode(&response).unwrap();
        assert_eq!((reply.id, &reply.revision), (request.id, &request.revision));
        let AnalysisResult::IvtTemplate(report) = reply.result else {
            panic!("IVT assessment result");
        };
        report.validate_for(&source, &declaration).unwrap();
        assert_eq!(report.promoter.dna, "TAATACGACTCACTATA");
        assert_eq!(&report.initiator_context.dna[..2], "AG");
        assert_eq!(
            report.transcript_spans,
            vec![ivt::MappedSpan {
                start: 17,
                end: 25,
                reverse: false
            }]
        );
        assert_eq!(
            report.linearization.as_ref().unwrap().kind,
            ivt::EndKind::Blunt
        );
        let poly_a = report.poly_a.as_ref().unwrap();
        assert_eq!(
            (
                poly_a.adenosines,
                poly_a.other_concrete_bases,
                poly_a.ambiguous_bases
            ),
            (4, 0, 0)
        );
        assert!(poly_a.terminal && poly_a.segmentation_complete);
        assert!(
            report
                .findings
                .iter()
                .any(|finding| finding.rule == ivt::Rule::ConsensusPromoter
                    && finding.verdict == ivt::Verdict::Pass)
        );
        assert!(
            report
                .findings
                .iter()
                .any(|finding| finding.rule == ivt::Rule::CapInitiator
                    && finding.verdict == ivt::Verdict::Pass)
        );
    }
    assert_eq!(source, before);
}

#[test]
fn ivt_raw_worker_boundaries_refuse_stale_context_profile_report_and_unsupported_versions() {
    let (source, declaration) = inputs();
    let request = Request {
        id: 9,
        revision: ivt::revision(&source, &declaration).unwrap(),
        job: AnalysisJob::AssessIvtTemplate {
            molecule: source.clone(),
            request: declaration.clone(),
        },
    };
    let mut stale = request.clone();
    stale.revision = "0".repeat(64);
    assert!(stale.encode().is_err());
    assert!(
        protocol::execution::process_request(&raw(b"CGQ\0", WORKER_SCHEMA_VERSION, &stale))
            .is_err()
    );
    let mut forged = request.clone();
    let AnalysisJob::AssessIvtTemplate {
        request: declaration,
        ..
    } = &mut forged.job
    else {
        unreachable!()
    };
    declaration.profile.identity = "0".repeat(64);
    assert!(forged.encode().is_err());
    assert!(
        protocol::execution::process_request(&raw(b"CGQ\0", WORKER_SCHEMA_VERSION, &forged))
            .is_err()
    );
    let response = Response::decode(
        &protocol::execution::process_request(&request.encode().unwrap()).unwrap(),
    )
    .unwrap();
    let mut changed = response.clone();
    let AnalysisResult::IvtTemplate(report) = &mut changed.result else {
        unreachable!()
    };
    report.promoter.dna.replace_range(0..1, "C");
    assert!(Response::decode(&raw(b"CGR\0", WORKER_SCHEMA_VERSION, &changed)).is_err());
    assert!(
        protocol::execution::process_request(&raw(b"CGQ\0", WORKER_SCHEMA_VERSION - 1, &request))
            .is_err()
    );
    assert!(Response::decode(&raw(b"CGR\0", WORKER_SCHEMA_VERSION - 1, &response)).is_err());
    let discovery = Request {
        id: 10,
        revision: "0".repeat(64),
        job: AnalysisJob::ListIvtProfiles,
    };
    assert!(discovery.encode().is_err());
    assert!(
        protocol::execution::process_request(&raw(b"CGQ\0", WORKER_SCHEMA_VERSION, &discovery))
            .is_err()
    );
}

#[test]
fn ivt_conflicts_remain_typed_evidence_and_missing_geometry_is_explicit() {
    let (source, mut declaration) = inputs();
    declaration.initiator = ivt::InitiatorRequirement::Gg;
    declaration.linearization = None;
    declaration.poly_a = None;
    let request = Request {
        id: 11,
        revision: ivt::revision(&source, &declaration).unwrap(),
        job: AnalysisJob::AssessIvtTemplate {
            molecule: source.clone(),
            request: declaration.clone(),
        },
    };
    let result = Response::decode(
        &protocol::execution::process_request(&request.encode().unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::IvtTemplate(report) = result.result else {
        panic!("scientific conflict remains a report");
    };
    assert!(report.linearization.is_none() && report.poly_a.is_none());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == ivt::Rule::ConflictingRequirements
                && finding.verdict == ivt::Verdict::Fail)
    );
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == ivt::Rule::Linearization
                && finding.verdict == ivt::Verdict::Indeterminate)
    );
    report.validate_for(&source, &declaration).unwrap();
    let oversized = Molecule::dna(
        "Oversized for admitted IVT model",
        &"A".repeat(ivt::MAX_TEMPLATE_BASES + 1),
        Topology::Linear,
    )
    .unwrap();
    declaration.source_revision = oversized.revision();
    assert!(declaration.validate_for(&oversized).is_err());
}
