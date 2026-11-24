use collagene_core::{
    Molecule, Topology,
    assembly_recipe::Preparation,
    overlap_primer_design as d,
    primers::DesignOptions,
    verification::{self as v, PlanInput, Status},
};
const A: &str = "ACGATTCGGTACCTAGGCTAACGTTAGCGATC";
const B: &str = "GCTTACCGATGACCTAGTACGGCATTCGAAGC";
const AB: &str = "ACGATTCGGTACCTAGGCTAACGTTAGCGATCGCTTACCGATGACCTAGTACGGCATTCGAAGC";
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture(circular: bool) -> (Vec<PlanInput>, Molecule, d::Request) {
    let inputs = vec![input("a", A), input("b", B)];
    let target = Molecule::dna(
        "independently supplied target",
        AB,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap();
    let request = d::Request {
        name: "deterministic explicit sources".into(),
        fragments: vec![
            d::Fragment {
                role: "left".into(),
                input_id: "a".into(),
                start: 0,
                end: 32,
                reverse: false,
            },
            d::Fragment {
                role: "right".into(),
                input_id: "b".into(),
                start: 0,
                end: 32,
                reverse: false,
            },
        ],
        overlap: 8,
        circular,
        primer_options: DesignOptions {
            min_length: 12,
            max_length: 12,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_tm_difference_c: 100.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        profile_options: Default::default(),
    };
    (inputs, target, request)
}
#[test]
fn literal_linear_and_circular_targets_use_only_next_prefix_and_original_inputs() {
    for circular in [false, true] {
        let (inputs, target, request) = fixture(circular);
        let before = serde_json::to_value(&inputs).unwrap();
        let report = d::design(&inputs, &target, &request).unwrap();
        report
            .validate_binding(&d::design_revision(&inputs, &target, &request).unwrap())
            .unwrap();
        assert_eq!(report.designs.len(), 2);
        assert_eq!(report.designs[0].pairs[0].reverse.primer.tail, "CGGTAAGC");
        assert!(report.designs[0].pairs[0].forward.primer.tail.is_empty());
        assert_eq!(
            report.designs[0].pairs[0].forward.primer.name,
            "left-forward"
        );
        assert_eq!(
            report.designs[1].pairs[0].reverse.primer.tail,
            if circular { "CGAATCGT" } else { "" }
        );
        assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
        let compiled =
            collagene_core::assembly_recipe::compile(&inputs, &target, &report.recipe).unwrap();
        let verification = v::verify(&compiled.plan, &target);
        assert_eq!(
            verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            verification.findings
        );
        assert_eq!(
            verification.products.last().unwrap().sequence,
            target.sequence
        );
    }
}
#[test]
fn reversed_upstream_core_places_tail_on_original_forward_primer() {
    let (inputs, _, mut request) = fixture(false);
    request.fragments[0].reverse = true;
    let target = Molecule::dna(
        "independent reverse-oriented target",
        "GATCGCTAACGTTAGCCTAGGTACCGAATCGTGCTTACCGATGACCTAGTACGGCATTCGAAGC",
        Topology::Linear,
    )
    .unwrap();
    let report = d::design(&inputs, &target, &request).unwrap();
    let pair = &report.designs[0].pairs[0];
    assert_eq!(pair.forward.primer.tail, "CGGTAAGC");
    assert!(pair.reverse.primer.tail.is_empty());
    let compiled =
        collagene_core::assembly_recipe::compile(&inputs, &target, &report.recipe).unwrap();
    let verification = v::verify(&compiled.plan, &target);
    assert_eq!(
        verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        verification.findings
    );
}
#[test]
fn actual_circular_source_interval_crosses_origin_without_rotating_source_or_target() {
    let (mut inputs, _, mut request) = fixture(false);
    inputs[0].molecule.topology = Topology::Circular;
    inputs[0].revision = inputs[0].molecule.revision();
    request.fragments[0].start = 24;
    request.fragments[0].end = 56;
    let target = Molecule::dna(
        "independent origin-crossing target",
        "TAGCGATCACGATTCGGTACCTAGGCTAACGTGCTTACCGATGACCTAGTACGGCATTCGAAGC",
        Topology::Linear,
    )
    .unwrap();
    let original = serde_json::to_value(&inputs[0]).unwrap();
    let report = d::design(&inputs, &target, &request).unwrap();
    let compiled =
        collagene_core::assembly_recipe::compile(&inputs, &target, &report.recipe).unwrap();
    let verification = v::verify(&compiled.plan, &target);
    assert_eq!(
        verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        verification.findings
    );
    assert_eq!(
        verification.products[0].sequence.as_str(),
        "TAGCGATCACGATTCGGTACCTAGGCTAACGTGCTTACCG"
    );
    assert_eq!(serde_json::to_value(&inputs[0]).unwrap(), original);
}
#[test]
fn target_mismatch_is_advisory_but_bad_ranges_tails_bounds_and_receipts_are_rejected() {
    let (inputs, target, request) = fixture(false);
    let original = d::design(&inputs, &target, &request).unwrap();
    let wrong = Molecule::dna(
        "different independent target",
        &"A".repeat(64),
        Topology::Linear,
    )
    .unwrap();
    let report = d::design(&inputs, &wrong, &request).unwrap();
    assert_ne!(original.context_revision, report.context_revision);
    assert_eq!(
        serde_json::to_value(&original.recipe).unwrap(),
        serde_json::to_value(&report.recipe).unwrap()
    );
    let compiled =
        collagene_core::assembly_recipe::compile(&inputs, &wrong, &report.recipe).unwrap();
    assert_eq!(v::verify(&compiled.plan, &wrong).status, Status::Invalid);
    let mut bad = request.clone();
    bad.primer_options.forward_tail = "A".into();
    assert!(d::design(&inputs, &target, &bad).is_err());
    bad = request.clone();
    bad.fragments[0].end = 33;
    assert!(d::design(&inputs, &target, &bad).is_err());
    bad = request.clone();
    bad.overlap = 81;
    assert!(d::design(&inputs, &target, &bad).is_err());
    bad = request.clone();
    bad.fragments = vec![bad.fragments[0].clone(); 33];
    assert!(d::design(&inputs, &target, &bad).is_err());
    let mut stale = inputs.clone();
    stale[0].revision = "0".repeat(64);
    assert!(d::design(&stale, &target, &request).is_err());
    let mut forged = original.clone();
    if let Preparation::ExactPcr { reverse, .. } = &mut forged.recipe.fragments[0].preparation {
        reverse.tail = "AAAAAAAA".into();
    }
    assert!(forged.validate().is_err());
    bad = request;
    bad.primer_options.min_gc_fraction = 1.;
    bad.primer_options.max_gc_fraction = 1.;
    assert!(d::design(&inputs, &target, &bad).is_err());
}
#[test]
fn repeated_pcr_geometries_are_not_promoted_by_a_ranked_pair() {
    let source = "AAAACCCCAAAACCCCAAAACCCC";
    let inputs = vec![input("repeat", source)];
    let target = Molecule::dna(
        "independent repeatedtarget",
        "AAAACCCCAAAACCCC",
        Topology::Linear,
    )
    .unwrap();
    let request = d::Request {
        name: "repeated template".into(),
        fragments: vec![d::Fragment {
            role: "copy".into(),
            input_id: "repeat".into(),
            start: 0,
            end: 16,
            reverse: false,
        }],
        overlap: 1,
        circular: false,
        primer_options: DesignOptions {
            min_length: 8,
            max_length: 8,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_tm_difference_c: 100.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        profile_options: Default::default(),
    };
    let report = d::design(&inputs, &target, &request).unwrap();
    assert!(report.warnings.iter().any(|w| w.contains("unique PCR")));
    let compiled =
        collagene_core::assembly_recipe::compile(&inputs, &target, &report.recipe).unwrap();
    let verification = v::verify(&compiled.plan, &target);
    assert_ne!(
        verification.status,
        Status::ConsistentUnderAssumptions,
        "repeated geometries must stay independently checked"
    );
}
