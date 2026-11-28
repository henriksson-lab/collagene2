//! Planner transport must retain source/target/options binding without promoting
//! suggested oligos into a verification report.
use collagene_core::{
    cloning::{
        overlap::primer_design::{
            self as overlap_primer_design, Fragment, Request as DesignRequest,
        },
        pcr::primers::DesignOptions,
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Vec<PlanInput>, Molecule, DesignRequest) {
    let inputs = [
        ("a", "ACGATTCGGTACCTAGGCTAACGTTAGCGATC"),
        ("b", "GCTTACCGATGACCTAGTACGGCATTCGAAGC"),
    ]
    .into_iter()
    .map(|(id, dna)| {
        let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
        PlanInput {
            id: id.into(),
            revision: molecule.revision(),
            molecule,
        }
    })
    .collect::<Vec<_>>();
    let target = Molecule::dna(
        "independent design",
        "ACGATTCGGTACCTAGGCTAACGTTAGCGATCGCTTACCGATGACCTAGTACGGCATTCGAAGC",
        Topology::Linear,
    )
    .unwrap();
    let request = DesignRequest {
        name: "source-derived homology".into(),
        fragments: inputs
            .iter()
            .map(|i| Fragment {
                role: i.id.clone(),
                input_id: i.id.clone(),
                start: 0,
                end: i.molecule.sequence.len(),
                reverse: false,
            })
            .collect(),
        overlap: 8,
        circular: false,
        primer_options: DesignOptions {
            min_length: 12,
            max_length: 12,
            max_tm_difference_c: 100.,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        profile_options: Default::default(),
    };
    (inputs, target, request)
}
fn message(
    inputs: Vec<PlanInput>,
    target: Molecule,
    request: DesignRequest,
    revision: String,
) -> Vec<u8> {
    protocol::request::Request::encode(&Request {
        id: 997,
        revision,
        job: AnalysisJob::DesignOverlapPrimers {
            inputs,
            target,
            request,
        },
    })
    .unwrap()
}
#[test]
fn new_appended_tags_are_roundtripped_and_rejected_under_pre_planner_headers() {
    let (inputs, target, request) = fixture();
    let revision = overlap_primer_design::design_revision(&inputs, &target, &request).unwrap();
    let bytes = message(inputs, target, request, revision.clone());
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        protocol::schema::WORKER_SCHEMA_VERSION
    );
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 95);
    let response = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 87);
    let decoded = protocol::response::Response::decode(&response).unwrap();
    assert_eq!(decoded.id, 997);
    assert_eq!(decoded.revision, revision);
    let AnalysisResult::OverlapPrimersDesigned(report) = decoded.result else {
        panic!("expected advisory planner report")
    };
    assert_eq!(report.context_revision, revision);
    assert_eq!(report.recipe.overlaps, vec![8]);
    assert_eq!(report.designs.len(), 2);
    assert_eq!(report.designs[0].pairs[0].reverse.primer.tail, "CGGTAAGC");
    assert!(report.designs[1].pairs[0].reverse.primer.tail.is_empty());
    for version in [1u16, 34] {
        let mut old = bytes.clone();
        old[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(protocol::execution::process_request(&old).is_err());
        let mut old = response.clone();
        old[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(protocol::response::Response::decode(&old).is_err());
    }
}
#[test]
fn changed_source_target_interval_or_options_cannot_reuse_an_old_context() {
    let (inputs, target, request) = fixture();
    let revision = overlap_primer_design::design_revision(&inputs, &target, &request).unwrap();
    let mut changed = request.clone();
    changed.overlap = 9;
    assert!(
        protocol::execution::process_request(&message(
            inputs.clone(),
            target.clone(),
            changed,
            revision.clone()
        ))
        .is_err()
    );
    let mut changed = request.clone();
    changed.fragments[0].end -= 1;
    assert!(
        protocol::execution::process_request(&message(
            inputs.clone(),
            target.clone(),
            changed,
            revision.clone()
        ))
        .is_err()
    );
    let wrong = target.replace(0, 1, "T").unwrap();
    assert!(
        protocol::execution::process_request(&message(
            inputs.clone(),
            wrong,
            request.clone(),
            revision.clone()
        ))
        .is_err()
    );
    let mut changed = inputs;
    changed[0].revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(&message(changed, target, request, revision)).is_err()
    );
}
