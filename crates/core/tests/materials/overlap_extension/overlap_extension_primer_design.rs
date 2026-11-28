//! Independent literal-edit expectations and source-geometry-selected oligos.
use collagene_core::{
    cloning::{
        assembly::edit_intent::Edit,
        overlap::{
            extension::primers::{self as design, Options, Report, Request},
            primer_design::Fragment,
        },
    },
    planning::plan::{Operation, Source},
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Sequence, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EDITED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const EDITED_REVERSE: &str = "TGACGATCTACAGGTCTGAAAATACGTAGCATCGGACT";
fn source() -> Molecule {
    Molecule::dna("Unchanged source40", ORIGINAL, Topology::Linear).unwrap()
}
fn options() -> Options {
    let mut options = Options {
        homology_flank: 8,
        ..Default::default()
    };
    options.primer_options.min_length = 8;
    options.primer_options.max_length = 8;
    options.primer_options.min_gc_fraction = 0.;
    options.primer_options.max_gc_fraction = 1.;
    options.primer_options.max_tm_difference_c = 100.;
    options.primer_options.max_homopolymer = 200;
    options.primer_options.max_three_prime_complementarity = 200;
    options
}
fn edit(start: usize, end: usize, expected: &str, replacement: &str) -> Edit {
    Edit {
        role: "insert".into(),
        start,
        end,
        expected: expected.into(),
        replacement: replacement.into(),
    }
}
fn request() -> Request {
    Request {
        template_strand: collagene_core::chemistry::Strand::Top,
        material_binding: None,
        fragment: Fragment {
            role: "insert".into(),
            input_id: "original".into(),
            start: 0,
            end: 40,
            reverse: false,
        },
        edits: vec![edit(16, 20, "GGAC", "TT")],
        options: options(),
    }
}
fn consistent(report: &Report, letters: &str) {
    report.validate_binding(&report.context_revision).unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(report.expected.sequence.as_str(), letters);
    assert_eq!(report.verification.products.len(), 4);
    assert_eq!(report.verification.products[3].sequence.as_str(), letters);
    assert!(report.search.iter().all(|s| s.unique_geometry));
    assert_eq!(report.designs.len(), 3);
}
#[test]
fn literal_replacement_selects_three_native_calls_and_verifies_actual_four_step_material() {
    let source = source();
    let original = serde_json::to_value(&source).unwrap();
    let req = request();
    let report = design::design(&source, &req, None).unwrap();
    consistent(&report, EDITED);
    assert_eq!(
        report.context_revision,
        design::design_revision(&source, &req, None).unwrap()
    );
    assert_eq!(report.source_revision, source.revision());
    assert_eq!(
        report
            .verification
            .products
            .iter()
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![26, 30, 38, 38]
    );
    assert_eq!(report.preparation.left.reverse.tail, "GGTCTGAAAA");
    assert_eq!(report.preparation.right.forward.tail, "GCTACGTATT");
    assert_eq!(report.preparation.extension.overlap, 18);
    assert_eq!(report.designs[0].options.target_start, 0);
    assert_eq!(report.designs[0].options.target_end, 16);
    assert_eq!(report.designs[1].options.target_start, 20);
    assert_eq!(report.designs[1].options.target_end, 40);
    assert_eq!(report.designs[2].options.target_start, 0);
    assert_eq!(report.designs[2].options.target_end, 40);
    assert!(
        report
            .designs
            .iter()
            .all(|d| d.source_revision == source.revision()
                && d.options.conditions == req.options.primer_options.conditions)
    );
    assert!(
        matches!(&report.plan.steps[0],Operation::Pcr{template:Source::Input(id),..}if id=="original")
    );
    assert!(
        matches!(&report.plan.steps[1],Operation::Pcr{template:Source::Input(id),..}if id=="original")
    );
    assert!(matches!(
        &report.plan.steps[3],
        Operation::Pcr {
            template: Source::Step(2),
            ..
        }
    ));
    assert_eq!(serde_json::to_value(&source).unwrap(), original);
    // A separate wrong target does not change the proposed source/oligo graph.
    let wrong = Molecule::dna("Independent unedited intent", ORIGINAL, Topology::Linear).unwrap();
    let mut wrong_plan = report.plan.clone();
    wrong_plan.target_revision = wrong.revision();
    let wrong_report = verification::engine::verify(&wrong_plan, &wrong);
    assert_eq!(wrong_report.status, Status::Invalid);
    assert_eq!(
        wrong_report.products.last().unwrap().sequence.as_str(),
        EDITED
    );
}
#[test]
fn reverse_contribution_edits_remain_source_oriented_and_declared_outer_tails_are_preserved() {
    let mut source = source();
    source.sequence = source.sequence.reverse_complement().unwrap();
    let mut req = request();
    req.fragment.reverse = true;
    req.options.outer_forward_tail = "GGA".into();
    req.options.outer_reverse_tail = "TCC".into();
    let report = design::design(&source, &req, None).unwrap();
    consistent(&report, "GGATGACGATCTACAGGTCTGAAAATACGTAGCATCGGACTGGA");
    assert_eq!(report.preparation.outer.forward.tail, "GGA");
    assert_eq!(report.preparation.outer.reverse.tail, "TCC");
    let core = &report.expected.sequence.as_str()[3..41];
    assert_eq!(core, EDITED_REVERSE);
    assert_eq!(
        Sequence::dna(core)
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str(),
        EDITED
    );
    assert!(req.fragment.reverse);
    let mut no_tails = req;
    no_tails.options.outer_forward_tail.clear();
    no_tails.options.outer_reverse_tail.clear();
    consistent(
        &design::design(&source, &no_tails, None).unwrap(),
        EDITED_REVERSE,
    );
}
#[test]
fn partial_and_origin_wrapping_intervals_copy_only_the_declared_original_contribution() {
    let text = format!("CCGGAATT{ORIGINAL}GGCAACCT");
    let source = Molecule::dna("Larger original", &text, Topology::Linear).unwrap();
    let mut req = request();
    req.fragment.start = 8;
    req.fragment.end = 48;
    let report = design::design(&source, &req, None).unwrap();
    consistent(&report, EDITED);
    assert_eq!(report.plan.inputs[0].molecule.sequence.len(), 56);
    let rotated = format!("{}{}", &text[30..], &text[..30]);
    let circular = Molecule::dna(
        "Original circle rotated for display",
        &rotated,
        Topology::Circular,
    )
    .unwrap();
    req.fragment.start = 34;
    req.fragment.end = 74;
    let report = design::design(&circular, &req, None).unwrap();
    consistent(&report, EDITED);
    assert!(report.designs[1].options.target_end > circular.sequence.len());
    assert_eq!(report.plan.inputs[0].molecule.topology, Topology::Circular);
    assert_eq!(report.expected.topology, Topology::Linear);
}
#[test]
fn multiple_disjoint_snp_indel_insert_delete_literals_define_one_small_central_bridge() {
    let source = source();
    for (edits, expected) in [
        (
            vec![edit(16, 20, "GGAC", "TT"), edit(23, 24, "A", "G")],
            "AGTCCGATGCTACGTATTTTCGGACCTGTAGATCGTCA",
        ),
        (
            vec![edit(20, 20, "", "AC")],
            "AGTCCGATGCTACGTAGGACACTTCAGACCTGTAGATCGTCA",
        ),
        (
            vec![edit(16, 20, "GGAC", "")],
            "AGTCCGATGCTACGTATTCAGACCTGTAGATCGTCA",
        ),
        (
            vec![edit(18, 19, "A", "T")],
            "AGTCCGATGCTACGTAGGTCTTCAGACCTGTAGATCGTCA",
        ),
    ] {
        let mut req = request();
        req.edits = edits;
        let report = design::design(&source, &req, None).unwrap();
        consistent(&report, expected);
    }
}
#[test]
fn native_ranked_search_skips_competing_short_domains_without_reading_expected_product() {
    let left = "AGTCCGATGCAGTCCGATTACGATTCGA";
    let text = format!("{left}GGACTTCAGACCTGTAGATCGTCA");
    let source = Molecule::dna(
        "Repeated short start, distinct longer endpoint",
        &text,
        Topology::Linear,
    )
    .unwrap();
    let mut req = request();
    req.fragment.end = text.len();
    req.edits = vec![edit(left.len(), left.len() + 4, "GGAC", "TT")];
    req.options.primer_options.max_length = 10;
    req.options.primer_options.target_tm_c = -50.;
    let report = design::design(&source, &req, None).unwrap();
    consistent(&report, &format!("{left}TTTTCAGACCTGTAGATCGTCA"));
    assert!(report.search[0].examined > 1, "{:?}", report.search);
    assert!(report.search[0].selected_index > 0);
    assert!(report.search[0].rejected > 0);
    assert!(report.search[2].unique_geometry);
}
#[test]
fn competing_or_unavailable_material_stays_diagnostic_and_terminal_scan_caps_are_not_success() {
    let repeated = Molecule::dna(
        "Two unchanged copies",
        &format!("{ORIGINAL}{ORIGINAL}"),
        Topology::Linear,
    )
    .unwrap();
    let report = design::design(&repeated, &request(), None).unwrap();
    assert_ne!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(report.verification.products.is_empty());
    assert!(!report.search[0].unique_geometry);
    assert!(!report.search[2].unique_geometry);
    assert_eq!(report.expected.sequence.as_str(), EDITED);
    let mut source = source();
    source.metadata.insert(
        "duplex_model".into(),
        "Unverified saved material summary".into(),
    );
    assert!(
        design::design(&source, &request(), None).is_err(),
        "Marked saved material cannot bypass its exact native binding as an original reference"
    );
    assert_eq!(
        source.metadata["duplex_model"],
        "Unverified saved material summary"
    );
    let mut req = request();
    req.options.max_comparisons = 18;
    let report = design::design(&self::source(), &req, None).unwrap();
    assert_eq!(report.verification.status, Status::Incomplete);
    assert_eq!(report.verification.products.len(), 4);
    assert_eq!(report.expected.sequence.as_str(), EDITED);
}
#[test]
fn wrong_literals_insertion_conflicts_boundary_edits_oligo_and_work_limits_fail_before_design() {
    let source = source();
    for kind in 0..8 {
        let mut req = request();
        match kind {
            0 => req.edits[0].expected = "AAAA".into(),
            1 => req.edits.push(edit(16, 16, "", "A")),
            2 => req.edits = vec![edit(0, 1, "A", "T")],
            3 => req.options.max_oligo_length = 16,
            4 => req.options.primer_options.forward_tail = "AC".into(),
            5 => req.options.primer_options.conditions.sodium_molar = f64::NAN,
            6 => req.options.outer_reverse_tail = "NNNN".into(),
            _ => req.options.max_candidates = 0,
        }
        assert!(
            design::preflight(&source, &req, None).is_err(),
            "invalid request {kind}"
        );
    }
    let mut req = request();
    req.options.primer_options.min_gc_fraction = 1.;
    req.options.primer_options.max_gc_fraction = 1.;
    assert!(design::design(&source, &req, None).is_err());
    let big = Molecule::dna(
        "Huge original search outside the requested contribution",
        &format!("{ORIGINAL}{}", "A".repeat(500_000)),
        Topology::Linear,
    )
    .unwrap();
    assert!(design::preflight(&big, &request(), None).is_err());
}
#[test]
fn documentary_context_shape_and_selected_pair_tampering_cannot_rebind_a_worker_receipt() {
    let report = design::design(&source(), &request(), None).unwrap();
    let decoded: Report = serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    decoded.validate_binding(&report.context_revision).unwrap();
    assert!(decoded.validate_binding(&"0".repeat(64)).is_err());
    for kind in 0..6 {
        let mut bad = report.clone();
        match kind {
            0 => bad.expected.sequence = Sequence::dna(ORIGINAL).unwrap(),
            1 => bad.preparation.right.forward.tail.push('A'),
            2 => bad.search[0].selected_index = usize::MAX,
            3 => bad.designs[1].options.target_start += 1,
            4 => bad.plan.inputs[0]
                .molecule
                .name
                .push_str("different source identity"),
            _ => bad.request.edits[0].replacement = "GG".into(),
        }
        assert!(
            bad.validate_binding(&report.context_revision).is_err(),
            "tamper{kind}"
        );
    }
}
