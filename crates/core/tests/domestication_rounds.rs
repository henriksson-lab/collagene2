//! Independent scattered-site genotypes; later PCRs use real preceding products.
use collagene_core::{
    Molecule, Topology, catalog_digest, codon_design, domestication as design,
    domestication_preparation as physical,
    protocol::{self, AnalysisJob, AnalysisResult},
    verification::{self, Operation, PlanInput, Source, Status},
};
const ORIGINAL: &str = "ACGTTGCACTGATCGAGTCCGGTCTCATGCTACGTAGGACGATCCAGTACGCTTATTCAGACCTGTAGATCGTCAGGTCTCTGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC";
const EXPECTED: &str = "ACGTTGCACTGATCGAGTCCGGCCTCATGCTACGTAGGACGATCCAGTACGCTTATTCAGACCTGTAGATCGTCAGGCCTCTGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC";
fn fixture() -> (PlanInput, physical::Request) {
    let source = Molecule::dna("Two scattered BsaI sites", ORIGINAL, Topology::Circular).unwrap();
    let request = design::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["bsai".into()],
        )
        .unwrap(),
        scope: vec![design::SourceRange { start: 0, end: 121 }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: Default::default(),
    };
    let patches = vec![
        design::Patch {
            position: 22,
            expected: "T".into(),
            replacement: "C".into(),
        },
        design::Patch {
            position: 77,
            expected: "T".into(),
            replacement: "C".into(),
        },
    ];
    let selection = design::Selection {
        context_revision: design::propose_revision(&source, &request).unwrap(),
        candidate_id: design::candidate_id(&request, &patches).unwrap(),
        patches,
    };
    let rounds = vec![
        physical::InverseRound {
            edit_start: 22,
            edit_end: 23,
            forward_length: 8,
            reverse_length: 8,
            reverse_tail_bases: 0,
        },
        physical::InverseRound {
            edit_start: 77,
            edit_end: 78,
            forward_length: 8,
            reverse_length: 8,
            reverse_tail_bases: 0,
        },
    ];
    (
        PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source,
        },
        physical::Request {
            design: request,
            selection,
            route: physical::Route::InversePcrRounds { rounds },
        },
    )
}
#[test]
fn scattered_changes_form_one_original_source_graph_with_actual_intermediate_consumption() {
    let (input, request) = fixture();
    let report = physical::materialize(&input, &request).unwrap();
    assert_eq!(report.target.sequence.as_str(), EXPECTED);
    report
        .validate_binding(&physical::revision(&input, &request).unwrap())
        .unwrap();
    assert_eq!(report.plan.inputs.len(), 1);
    assert_eq!(report.plan.inputs[0].molecule.sequence.as_str(), ORIGINAL);
    assert_eq!(report.plan.steps.len(), 6);
    assert!(matches!(
        &report.plan.steps[3],
        Operation::Pcr {
            template: Source::Step(2),
            ..
        }
    ));
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    assert!(verification::plan_target_matches(
        &report.plan,
        report.verification.products.last().unwrap(),
        &report.target
    ));
    let physical::MethodEvidence::InversePcrRounds { rounds } = &report.method else {
        panic!("round evidence")
    };
    assert_eq!(rounds[0].source_origin, 0);
    assert_eq!(rounds[1].source_origin, 22);
    assert_eq!(rounds[1].options.edit_start, 55);
    assert_eq!(rounds[0].forward.tail, "C");
    assert_eq!(rounds[1].forward.tail, "C");
    let html = verification::export_protocol_html(&report.plan, &report.target);
    assert!(html.contains("Domestication round1 forward"));
    assert!(html.contains("Domestication round2 forward"));
    assert!(
        collagene_core::pdf::export_protocol_pdf(&report.plan, &report.target)
            .unwrap()
            .bytes
            .starts_with(b"%PDF-")
    );
    let mut tampered = report.clone();
    let Operation::Pcr { template, .. } = &mut tampered.plan.steps[3] else {
        unreachable!()
    };
    *template = Source::Input(physical::INPUT_ID.into());
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    assert_eq!(
        verification::verify(&tampered.plan, &report.target).status,
        Status::Invalid
    );
}
#[test]
fn coverage_order_and_work_limits_are_explicit_not_automatic_clustering() {
    let (input, mut request) = fixture();
    let physical::Route::InversePcrRounds { rounds } = &mut request.route else {
        unreachable!()
    };
    rounds[1].edit_start = 22;
    rounds[1].edit_end = 23;
    assert!(physical::preflight(&input, &request).is_err());
    let (input, mut request) = fixture();
    let physical::Route::InversePcrRounds { rounds } = &mut request.route else {
        unreachable!()
    };
    rounds[1].edit_start = 50;
    rounds[1].edit_end = 51;
    assert!(physical::preflight(&input, &request).is_err());
    let (input, mut request) = fixture();
    let before = physical::revision(&input, &request).unwrap();
    let physical::Route::InversePcrRounds { rounds } = &mut request.route else {
        unreachable!()
    };
    rounds.swap(0, 1);
    assert_ne!(before, physical::revision(&input, &request).unwrap());
    assert_eq!(
        physical::materialize(&input, &request)
            .unwrap()
            .target
            .sequence
            .as_str(),
        EXPECTED
    );
}
#[test]
fn worker_roundtrip_retains_rounds_and_unchanged_independent_target() {
    let (input, request) = fixture();
    let revision = physical::revision(&input, &request).unwrap();
    let job = AnalysisJob::PrepareDomestication { input, request };
    let request = protocol::Request {
        id: 99,
        revision: revision.clone(),
        job,
    };
    let response = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::DomesticationPrepared(report) = response.result else {
        panic!("prepared domestication")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(report.target.sequence.as_str(), EXPECTED);
    assert_eq!(report.plan.inputs.len(), 1);
}

#[test]
fn partial_reverse_tail_and_original_origin_crossing_keep_actual_round_frames() {
    let (input, mut request) = fixture();
    let physical::Route::InversePcrRounds { rounds } = &mut request.route else {
        unreachable!()
    };
    // Whole first site is replaced through both primer tails, including the
    // changed base in the reverse tail. The next original interval wraps121→0.
    for round in rounds.iter_mut() {
        round.forward_length = 18;
        round.reverse_length = 18;
    }
    rounds[0].edit_start = 20;
    rounds[0].edit_end = 26;
    rounds[0].reverse_tail_bases = 3;
    rounds[1].edit_start = 75;
    rounds[1].edit_end = 131;
    rounds[1].reverse_tail_bases = 4;
    let report = physical::materialize(&input, &request).unwrap();
    assert_eq!(report.target.sequence.as_str(), EXPECTED);
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    let physical::MethodEvidence::InversePcrRounds { rounds } = &report.method else {
        unreachable!()
    };
    assert_eq!(rounds[0].forward.tail, "CTC");
    assert_eq!(rounds[0].reverse.tail, "GCC");
    assert_eq!(rounds[1].source_origin, 23);
    assert_eq!(
        (rounds[1].options.edit_start, rounds[1].options.edit_end),
        (52, 108)
    );
    assert!(
        report.verification.products[2]
            .sequence
            .as_str()
            .starts_with("CTCATGCTACGTAGGAC")
    );
    let intermediate = Molecule::dna("Independent first round", "ACGTTGCACTGATCGAGTCCGGCCTCATGCTACGTAGGACGATCCAGTACGCTTATTCAGACCTGTAGATCGTCAGGTCTCTGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC", Topology::Circular).unwrap();
    assert!(verification::plan_target_matches(
        &report.plan,
        &report.verification.products[2],
        &intermediate
    ));
    assert!(verification::plan_target_matches(
        &report.plan,
        report.verification.products.last().unwrap(),
        &report.target
    ));
    let mut tampered = report.clone();
    let physical::MethodEvidence::InversePcrRounds { rounds } = &mut tampered.method else {
        unreachable!()
    };
    rounds[1].source_origin = 22;
    assert!(tampered.validate_binding(&report.context_revision).is_err());
}
