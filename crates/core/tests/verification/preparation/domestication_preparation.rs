//! Actual original-template materialization, independent literal target oracles,
//! and explicit refusal of routes that cannot encode every selected change.
use collagene_core::{
    cloning::{
        golden_gate::domestication::{self as design, preparation as physical},
        pcr::mutagenesis,
        restriction::catalog as catalog_digest,
    },
    planning::plan::PlanInput,
    sequence::codon_design,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};

fn fixture(dna: &str, patches: Vec<(usize, &str, &str)>) -> (PlanInput, physical::Request) {
    let source = Molecule::dna("Original template", dna, Topology::Circular).unwrap();
    let n = source.sequence.len();
    let request = design::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["bsai".into()],
        )
        .unwrap(),
        scope: vec![design::SourceRange { start: 0, end: n }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: design::Limits::default(),
    };
    let patches = patches
        .into_iter()
        .map(|(position, expected, replacement)| design::Patch {
            position,
            expected: expected.into(),
            replacement: replacement.into(),
        })
        .collect::<Vec<_>>();
    let selection = design::Selection {
        context_revision: design::propose_revision(&source, &request).unwrap(),
        candidate_id: design::candidate_id(&request, &patches).unwrap(),
        patches,
    };
    (
        PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source,
        },
        physical::Request {
            origin: None,
            design: request,
            selection,
            route: physical::Route::InversePcr {
                edit_start: 10,
                edit_end: 11,
                forward_length: 8,
                reverse_length: 8,
                reverse_tail_bases: 0,
            },
        },
    )
}
fn simple() -> (PlanInput, physical::Request) {
    fixture("AGTCCGATGGTCTCGTAGGACTTCAGACCTGTA", vec![(10, "T", "C")])
}
fn repin(input: &PlanInput, request: &mut physical::Request) {
    request.selection.context_revision =
        design::propose_revision(&input.molecule, &request.design).unwrap();
    request.selection.candidate_id =
        design::candidate_id(&request.design, &request.selection.patches).unwrap();
}
#[test]
fn a_checked_literal_site_change_is_supplied_by_real_primers_from_original_dna() {
    let (input, request) = simple();
    let original = input.molecule.clone();
    let mut expected = original.sequence.as_str().as_bytes().to_vec();
    expected[10] = b'C';
    let expected = String::from_utf8(expected).unwrap();
    let bound = physical::revision(&input, &request, None).unwrap();
    let report = physical::materialize(&input, &request, None).unwrap();
    report.validate_binding(&bound).unwrap();
    assert_eq!(report.target.sequence.as_str(), expected);
    assert!(
        report
            .check
            .sites_before
            .iter()
            .any(|site| site.in_scope && !site.retained)
    );
    assert!(
        report
            .check
            .sites_after
            .iter()
            .all(|site| !site.in_scope || site.retained)
    );
    let physical::MethodEvidence::InversePcr {
        options,
        forward,
        reverse,
        ..
    } = &report.method
    else {
        panic!("inverse evidence");
    };
    assert_eq!(forward.tail, "C");
    assert!(reverse.tail.is_empty());
    assert_eq!(
        report.plan.inputs[0].molecule.sequence.as_str(),
        original.sequence.as_str()
    );
    assert_eq!(report.plan.inputs[0].revision, original.revision());
    assert_ne!(report.target.revision(), original.revision());
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(report.verification.reconstructed_count(), 3);
    assert!(verification::engine::plan_target_matches(
        &report.plan,
        report
            .verification
            .product(2)
            .expect("selected step must reconstruct"),
        &report.target
    ));
    assert_eq!(input.molecule, original);
    // A reference-only domesticated sequence is not a replacement material stock.
    assert!(
        mutagenesis::compile_inverse_pcr(&report.target, &report.target, None, options).is_err()
    );
}
#[test]
fn origin_spanning_edits_keep_original_coordinates_and_encode_both_tails() {
    let dna = "TCTCAGTCCGATGCTACGTAGGACTTCAGACCTGTAGG";
    let n = dna.len();
    let (input, mut request) = fixture(dna, vec![(0, "T", "C"), (n - 1, "G", "A")]);
    request.route = physical::Route::InversePcr {
        edit_start: n - 2,
        edit_end: n + 2,
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 2,
    };
    let report = physical::materialize(&input, &request, None).unwrap();
    let mut expected = dna.as_bytes().to_vec();
    expected[0] = b'C';
    expected[n - 1] = b'A';
    assert_eq!(
        report.target.sequence.as_str(),
        String::from_utf8(expected).unwrap()
    );
    let physical::MethodEvidence::InversePcr {
        options,
        forward,
        reverse,
        ..
    } = &report.method
    else {
        panic!("inverse evidence");
    };
    assert_eq!(options.inserted, "GACC");
    assert_eq!(forward.tail, "CC");
    assert_eq!(reverse.tail, "TC");
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(report.plan.inputs[0].molecule.sequence.as_str(), dna);
}
#[test]
fn stale_original_candidate_and_out_of_interval_changes_are_refused() {
    let (input, request) = simple();
    let mut wrong = input.clone();
    wrong.id = "recoded template".into();
    assert!(physical::preflight(&wrong, &request, None).is_err());
    wrong = input.clone();
    wrong.revision = "a".repeat(64);
    assert!(physical::preflight(&wrong, &request, None).is_err());
    let mut stale = request.clone();
    stale.selection.candidate_id = "a".repeat(64);
    assert!(physical::preflight(&input, &stale, None).is_err());
    let mut uncovered = request.clone();
    uncovered.route = physical::Route::InversePcr {
        edit_start: 1,
        edit_end: 2,
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    assert!(physical::preflight(&input, &uncovered, None).is_err());
    let mut forged = request.clone();
    forged.selection.patches[0].expected = "A".into();
    repin(&input, &mut forged);
    assert!(physical::materialize(&input, &forged, None).is_err());
    let mut outside = request.clone();
    outside.selection.patches[0].position = usize::MAX;
    repin(&input, &mut outside);
    assert!(physical::preflight(&input, &outside, None).is_err());
}
#[test]
fn distant_substitutions_and_a_201_base_ordered_oligo_need_another_explicit_route() {
    // Keep enough original DNA outside the distant replacement for both
    // annealing domains; this fixture isolates the ordered-oligo limit.
    let dna = format!("GGTCTC{}GGTCTC{}", "ACGT".repeat(120), "AGTC".repeat(8));
    let second = 6 + 4 * 120;
    let (input, mut distant) = fixture(&dna, vec![(2, "T", "C"), (second + 2, "T", "C")]);
    distant.route = physical::Route::InversePcr {
        edit_start: 2,
        edit_end: second + 3,
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    let error = physical::preflight(&input, &distant, None)
        .unwrap_err()
        .to_string();
    assert!(error.contains("200"));
    let (input, mut long) = fixture(&dna, vec![(2, "T", "C")]);
    long.route = physical::Route::InversePcr {
        edit_start: 0,
        edit_end: 193,
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    assert!(
        physical::preflight(&input, &long, None)
            .unwrap_err()
            .to_string()
            .contains("200")
    );
    long.route = physical::Route::InversePcr {
        edit_start: 0,
        edit_end: 193,
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 1,
    };
    physical::preflight(&input, &long, None).unwrap(); // exactly200, still only a declaration
}
#[test]
fn independent_target_and_method_evidence_tampering_fail_under_current_compiler() {
    let (input, request) = simple();
    let report = physical::materialize(&input, &request, None).unwrap();
    let bound = report.context_revision.clone();
    let physical::MethodEvidence::InversePcr { options, .. } = &report.method else {
        panic!("inverse evidence");
    };
    let mut changed = report.clone();
    if let physical::MethodEvidence::InversePcr { forward, .. } = &mut changed.method {
        forward.tail = "A".into();
    }
    assert!(changed.validate_binding(&bound).is_err());
    changed = report.clone();
    if let physical::MethodEvidence::InversePcr {
        forward_location, ..
    } = &mut changed.method
    {
        *forward_location = collagene_model::Location::range(0, 8);
    }
    assert!(changed.validate_binding(&bound).is_err());
    changed = report.clone();
    if let physical::MethodEvidence::InversePcr {
        retained_template_length,
        ..
    } = &mut changed.method
    {
        *retained_template_length += 1;
    }
    assert!(changed.validate_binding(&bound).is_err());
    changed = report.clone();
    if let physical::MethodEvidence::InversePcr { forward, .. } = &mut changed.method {
        forward.tail = "A".into();
    }
    if let collagene_core::planning::plan::Operation::Pcr { forward, .. } =
        &mut changed.plan.steps[0]
    {
        forward.tail = "A".into();
    }
    changed.verification.plan_revision =
        collagene_core::planning::plan::DesignPlan::revision(&changed.plan);
    assert!(changed.validate_binding(&bound).is_err());
    changed = report.clone();
    changed.target = Molecule::dna("Wrong target", "ACGTACGTACGTACGT", Topology::Circular).unwrap();
    assert!(changed.validate_binding(&bound).is_err());
    let original =
        mutagenesis::design_inverse_pcr(&input.molecule, &report.target, None, options).unwrap();
    assert_eq!(
        serde_json::to_value(&original.plan).unwrap(),
        serde_json::to_value(&report.plan).unwrap()
    );
    assert_eq!(original.verification.status, report.verification.status);
    let wrong =
        mutagenesis::compile_inverse_pcr(&input.molecule, &changed.target, None, options).unwrap();
    assert_eq!(
        verification::engine::verify(&wrong.plan, &changed.target).status,
        Status::Invalid
    );
}
