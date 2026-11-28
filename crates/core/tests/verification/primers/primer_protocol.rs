use collagene_core::{
    cloning::pcr::primers::{BindingOptions, DesignOptions, ThermoConditions},
    planning::plan::DesignPlan,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Error, Location, Molecule, Topology};
fn run(source: &Molecule, job: AnalysisJob) -> AnalysisResult {
    let request = Request {
        id: 71,
        revision: source.revision(),
        job,
    };
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(response.id, 71);
    assert_eq!(response.revision, source.revision());
    response.result
}
#[test]
fn primer_conditions_and_exact_binding_survive_binary_dispatch() {
    let source = Molecule::dna("primer worker", "CCCGTTGATTA", Topology::Linear).unwrap();
    let AnalysisResult::PrimerEvidence(report) = run(
        &source,
        AnalysisJob::AnalyzePrimer {
            molecule: source.clone(),
            primer: Primer::new("F", "CGTTGA", "AA").unwrap(),
            conditions: ThermoConditions {
                sodium_molar: 1.,
                total_strand_molar: 500e-9,
            },
            binding_options: BindingOptions::default(),
        },
    ) else {
        panic!("wrong response")
    };
    assert_eq!(report.bindings.source_revision, source.revision());
    assert_eq!(report.analysis.ordered_sequence, "AACGTTGA");
    assert_eq!(report.bindings.hits.len(), 1);
    assert_eq!(report.bindings.hits[0].start, 2);
    assert_eq!(report.bindings.hits[0].three_prime_exact, 6);
    let tm = report.analysis.tm.unwrap();
    assert!((tm.delta_h_kcal_per_mol + 41.2).abs() < 1e-9);
    assert!((tm.delta_s_cal_per_mol_kelvin_at_1m + 115.4).abs() < 1e-9);
}
#[test]
fn anchored_primer_pair_design_retains_tails_source_and_candidate_coordinates() {
    let source = Molecule::dna("design worker", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let options = DesignOptions {
        target_start: 0,
        target_end: 16,
        min_length: 8,
        max_length: 8,
        target_tm_c: 20.,
        max_tm_difference_c: 100.,
        min_gc_fraction: 0.,
        max_gc_fraction: 1.,
        max_homopolymer: 20,
        max_three_prime_complementarity: 20,
        forward_tail: "AA".into(),
        reverse_tail: "TT".into(),
        ..Default::default()
    };
    let AnalysisResult::PrimerDesign(report) = run(
        &source,
        AnalysisJob::DesignPrimers {
            molecule: source.clone(),
            options,
        },
    ) else {
        panic!("wrong response")
    };
    assert_eq!(report.source_revision, source.revision());
    assert_eq!(report.pairs.len(), 1);
    let pair = &report.pairs[0];
    assert_eq!(pair.forward.location, Location::range(0, 8));
    assert_eq!(
        pair.reverse.location,
        Location::Complement(Box::new(Location::range(8, 16)))
    );
    assert_eq!(pair.reverse.primer.ordered_sequence(), "TTGTTACCCC");
    assert_eq!(pair.product_length, 20);
}
#[test]
fn embedded_pdf_worker_reverifies_instead_of_accepting_a_claimed_status() {
    let target = Molecule::dna("PDF target", "ACGT", Topology::Linear).unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Empty invalid plan".into(),
        target_revision: target.revision(),
        inputs: vec![],
        steps: vec![],
    };
    let AnalysisResult::ProtocolPdf(report) = run(
        &target,
        AnalysisJob::ExportProtocolPdf {
            target: target.clone(),
            plan,
        },
    ) else {
        panic!("wrong response")
    };
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.target_revision, target.revision());
    assert!(report.bytes.starts_with(b"%PDF-"));
    assert!(report.pages > 0);
}
#[test]
fn new_variants_cannot_be_mislabeled_as_old_worker_schemas_and_stale_inputs_fail() {
    let source = Molecule::dna("binding", "ACGTACGT", Topology::Linear).unwrap();
    let request = Request {
        id: 3,
        revision: source.revision(),
        job: AnalysisJob::AnalyzePrimer {
            molecule: source.clone(),
            primer: Primer::new("F", "ACGT", "").unwrap(),
            conditions: Default::default(),
            binding_options: Default::default(),
        },
    };
    for old in [1u16, 2, 3] {
        let mut bytes = protocol::request::Request::encode(&request).unwrap();
        bytes[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(matches!(
            protocol::execution::process_request(&bytes),
            Err(Error::Unsupported(_))
        ));
    }
    let mut reply = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    reply[4..6].copy_from_slice(&3u16.to_le_bytes());
    assert!(matches!(
        protocol::response::Response::decode(&reply),
        Err(Error::Unsupported(_))
    ));
    let mut stale = request;
    stale.revision = "old".into();
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
}
