//! Binary screening preserves explicit mixtures, source pins and bounded uncertainty.
use collagene_core::{
    chemistry::{
        methylation::MethylationState,
        tracked_duplex::{GenerationAssumptions, checked_selection},
    },
    cloning::screening::{
        self, Component, DecoyKind, IdealSizeReadout, Limits, Method, Role, Subject, Test,
    },
    material::derivation::SourceOrigin,
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
    },
};
use collagene_model::{Molecule, Topology};

fn declaration() -> screening::Request {
    let components = [
        ("intended", "AAAAGAATTCTTTT"),
        ("parent", "CCCCCCCCCCCCCCCCCCCC"),
    ]
    .into_iter()
    .map(|(id, sequence)| {
        let molecule = Molecule::dna(id, sequence, Topology::Linear).unwrap();
        Component {
            id: id.into(),
            label: id.into(),
            revision: molecule.revision(),
            molecule,
            origin: SourceOrigin::OriginalReference,
            original_generation: Some(GenerationAssumptions {
                host: None,
                dam: MethylationState::Unmodified,
                dcm: MethylationState::Unmodified,
                cpg: MethylationState::Unmodified,
            }),
        }
    })
    .collect();
    screening::Request {
        components,
        subjects: vec![
            Subject {
                id: "intended".into(),
                label: "Intended".into(),
                role: Role::Intended,
                components: vec!["intended".into()],
            },
            Subject {
                id: "carryover".into(),
                label: "Parent carryover mixture".into(),
                role: Role::Decoy(DecoyKind::ParentCarryover),
                components: vec!["intended".into(), "parent".into()],
            },
        ],
        tests: vec![Test {
            id: "ecori".into(),
            label: "EcoRI".into(),
            method: Method::Digest {
                selection: checked_selection("EcoRI").unwrap(),
            },
        }],
        readout: IdealSizeReadout {
            counterpart_tolerance_twice_bp: 0,
        },
        limits: Limits::default(),
    }
}

fn request(declaration: screening::Request) -> Request {
    Request {
        id: 72,
        revision: screening::revision(&declaration).unwrap(),
        job: AnalysisJob::AssessScreeningStrategy {
            request: declaration,
        },
    }
}

fn raw<T: serde::Serialize>(magic: &[u8; 4], value: &T) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    bytes
}

#[test]
fn binary_and_session_screening_keep_literal_sizes_mixture_and_exact_context() {
    let declaration = declaration();
    let request = request(declaration.clone());
    let encoded = request.encode().unwrap();
    let mut session = protocol::session::AnalysisSession::default();
    for encoded_response in [
        protocol::execution::process_request(&encoded).unwrap(),
        session.process_request(&encoded).unwrap(),
    ] {
        let response = Response::decode(&encoded_response).unwrap();
        assert_eq!(response.id, request.id);
        assert_eq!(response.revision, request.revision);
        let AnalysisResult::ScreeningStrategy(report) = &response.result else {
            panic!("screening report")
        };
        report.validate_for(&declaration).unwrap();
        assert_eq!(report.tests[0].signatures[0].sizes_twice_bp, [14]);
        assert_eq!(report.tests[0].signatures[1].sizes_twice_bp, [14, 40]);
        assert_eq!(report.tests[0].separated_decoys, 1);
        assert_eq!(report.ranked_test_indices, [0]);
        assert_eq!(report.panel.test_indices, [0]);
        let mut forged = response;
        let AnalysisResult::ScreeningStrategy(report) = &mut forged.result else {
            unreachable!()
        };
        report.cells[0].sizes_twice_bp = vec![999];
        // A valid transport checksum cannot turn altered scientific evidence into a valid report.
        assert!(Response::decode(&raw(b"CGR\0", &forged)).is_err());
    }
}

#[test]
fn screening_transport_refuses_stale_sources_policy_and_obsolete_internal_schema() {
    let request = request(declaration());
    let mut stale = request.clone();
    stale.revision = "0".repeat(64);
    assert!(stale.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &stale)).is_err());
    let mut stale_source = request.clone();
    let AnalysisJob::AssessScreeningStrategy {
        request: declaration,
    } = &mut stale_source.job
    else {
        unreachable!()
    };
    declaration.components[0].revision = "0".repeat(64);
    assert!(stale_source.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &stale_source)).is_err());
    let mut old = request.encode().unwrap();
    old[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
}

#[test]
fn binary_limited_screening_retains_unresolved_comparisons_instead_of_inferred_absence() {
    let mut declaration = declaration();
    declaration.limits.max_cleavage_work = 1;
    let request = request(declaration.clone());
    let bytes = protocol::execution::process_request(&request.encode().unwrap()).unwrap();
    let response = Response::decode(&bytes).unwrap();
    let AnalysisResult::ScreeningStrategy(report) = response.result else {
        panic!("screening report")
    };
    report.validate_for(&declaration).unwrap();
    assert_eq!(report.cleavage_work, 1);
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        screening::Distinction::Unresolved
    );
    assert!(report.panel.separated_decoy_indices.is_empty());
}
