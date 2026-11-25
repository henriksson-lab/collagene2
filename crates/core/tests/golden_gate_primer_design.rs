use collagene_core::{
    Molecule, Sequence, Topology, assembly_recipe::Preparation, golden_gate_primer_design as d,
    verification::PlanInput,
};
const SOURCE: &str = "GCTTACCGATGACCTAGTACGGCATTCGAAGC";
const REVERSED: &str = "GCTTCGAATGCCGTACTAGGTCATCGGTAAGC";
fn fixture(reverse: bool, enzyme: d::GoldenGateEnzyme) -> (Vec<PlanInput>, Molecule, d::Request) {
    let molecule = Molecule::dna("original", SOURCE, Topology::Linear).unwrap();
    let inputs = vec![PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    }];
    let target = Molecule::dna(
        "independent",
        if reverse { REVERSED } else { SOURCE },
        Topology::Linear,
    )
    .unwrap();
    let k = if enzyme == d::GoldenGateEnzyme::SapI {
        3
    } else {
        4
    };
    let request = d::Request {
        name: "literal source preparation".into(),
        fragments: vec![d::Fragment {
            role: "insert".into(),
            input_id: "source".into(),
            start: 0,
            end: 32,
            reverse,
            forward_length: 8,
            reverse_length: 8,
            right_overhang: if k == 3 { "CTG" } else { "CTGA" }.into(),
        }],
        enzyme,
        clamp: "ACGT".into(),
        spacer: if enzyme == d::GoldenGateEnzyme::BbsI {
            "AA"
        } else {
            "A"
        }
        .into(),
    };
    (inputs, target, request)
}
#[test]
fn literal_forward_bsa_i_and_sap_i_preserve_top_and_actual_bottom_fusions() {
    for enzyme in [
        d::GoldenGateEnzyme::BsaI,
        d::GoldenGateEnzyme::SapI,
        d::GoldenGateEnzyme::BsmBI,
        d::GoldenGateEnzyme::BbsI,
    ] {
        let (inputs, target, r) = fixture(false, enzyme);
        let before = serde_json::to_value(&inputs).unwrap();
        let report = d::design(&inputs, &target, &r).unwrap();
        report
            .validate_binding(&d::preview_revision(&inputs, &target, &r).unwrap())
            .unwrap();
        let row = &report.rows[0];
        assert_eq!(row.oriented_top, SOURCE);
        assert_eq!(row.forward.annealing.as_str(), "GCTTACCG");
        assert_eq!(row.reverse.annealing.as_str(), "GCTTCGAA");
        assert_eq!(row.forward_tail.fusion, "");
        assert_eq!(
            row.right.overhang,
            if enzyme == d::GoldenGateEnzyme::SapI {
                "CAG"
            } else {
                "TCAG"
            }
        );
        assert_eq!(
            row.left.overhang,
            if enzyme == d::GoldenGateEnzyme::SapI {
                "GCT"
            } else {
                "GCTT"
            }
        );
        assert_eq!(row.fragment_index, 1);
        assert_eq!(row.selection.enzymes.len(), 1);
        assert!(matches!(
            report.fragments[0].preparation,
            Preparation::PcrCatalogDigest { repair: None, .. }
        ));
        assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
        if enzyme == d::GoldenGateEnzyme::BsaI {
            assert_eq!(row.forward.tail, "ACGTGGTCTCA");
            assert_eq!(row.reverse.tail, "ACGTGGTCTCATCAG");
        }
    }
}
#[test]
fn reverse_rows_use_original_bottom_geometry_not_rc_of_original_full_top() {
    for enzyme in [d::GoldenGateEnzyme::BsaI, d::GoldenGateEnzyme::SapI] {
        let (inputs, target, r) = fixture(true, enzyme);
        let report = d::design(&inputs, &target, &r).unwrap();
        let row = &report.rows[0];
        assert_eq!(row.oriented_top, REVERSED);
        assert_eq!(row.forward.annealing.as_str(), "GCTTACCG");
        if enzyme == d::GoldenGateEnzyme::BsaI {
            assert_eq!(row.forward_tail.fusion, "TCAG");
            assert_eq!(row.reverse_tail.fusion, "GCTT");
            assert_eq!(row.reverse.annealing.as_str(), "CGAATGCC");
            assert_eq!(row.left.overhang, "GCTT");
            assert_eq!(row.right.overhang, "TCAG");
        } else {
            assert_eq!(row.forward_tail.fusion, "CAG");
            assert_eq!(row.reverse_tail.fusion, "GCT");
            assert_eq!(row.left.overhang, "GCT");
            assert_eq!(row.right.overhang, "CAG");
        }
    }
}
#[test]
fn circular_origin_interval_and_wrong_target_do_not_choose_primers() {
    let m = Molecule::dna(
        "circle",
        "ACGATTCGGTACCTAGGCTAACGTTAGCGATC",
        Topology::Circular,
    )
    .unwrap();
    let inputs = vec![PlanInput {
        id: "source".into(),
        revision: m.revision(),
        molecule: m,
    }];
    let (_, target, mut r) = fixture(false, d::GoldenGateEnzyme::BsaI);
    r.fragments[0].start = 24;
    r.fragments[0].end = 48;
    let report = d::design(&inputs, &target, &r).unwrap();
    assert_eq!(report.rows[0].oriented_top, "TAGCGATCACGATTCGGTACCTAG");
    let wrong = Molecule::dna("other target", "AAAA", Topology::Circular).unwrap();
    let second = d::design(&inputs, &wrong, &r).unwrap();
    assert_ne!(report.context_revision, second.context_revision);
    assert_eq!(
        serde_json::to_value(&report.fragments).unwrap(),
        serde_json::to_value(&second.fragments).unwrap()
    );
    assert!(report.validate_binding(&second.context_revision).is_err());
}
#[test]
fn internal_sites_binding_alternatives_and_invalid_bounds_withhold_proposals() {
    let (mut inputs, target, mut r) = fixture(false, d::GoldenGateEnzyme::BsaI);
    let m = Molecule::dna(
        "internal",
        "GCTTACCGGGTCTCAATACGGCATTCGAAGC",
        Topology::Linear,
    )
    .unwrap();
    r.fragments[0].end = m.sequence.len();
    inputs[0].revision = m.revision();
    inputs[0].molecule = m;
    assert!(d::design(&inputs, &target, &r).is_err());
    let (inputs, target, mut r) = fixture(false, d::GoldenGateEnzyme::BsaI);
    r.spacer = "AA".into();
    assert!(d::design(&inputs, &target, &r).is_err());
    r.spacer = "A".into();
    r.fragments[0].right_overhang = "NNNN".into();
    assert!(d::design(&inputs, &target, &r).is_err());
    r.fragments[0].right_overhang = "CTGA".into();
    r.fragments[0].end = usize::MAX;
    assert!(d::design(&inputs, &target, &r).is_err());
    r.fragments[0].end = 32;
    let mut repeated = inputs.clone();
    repeated[0].molecule.sequence = Sequence::dna(&format!("{SOURCE}{SOURCE}")).unwrap();
    repeated[0].revision = repeated[0].molecule.revision();
    assert!(d::design(&repeated, &target, &r).is_err());
}
#[test]
fn source_and_report_budget_context_and_geometry_guards() {
    let (mut inputs, target, r) = fixture(false, d::GoldenGateEnzyme::BsaI);
    let report = d::design(&inputs, &target, &r).unwrap();
    let mut forged = report.clone();
    forged.rows[0].reverse_tail.fusion = "AAAA".into();
    assert!(forged.validate().is_err());
    let mut forged = report.clone();
    forged.fragments[0].reverse = true;
    assert!(forged.validate().is_err());
    inputs[0].revision = "0".repeat(64);
    assert!(d::design(&inputs, &target, &r).is_err());
    let (inputs, target, mut r) = fixture(false, d::GoldenGateEnzyme::BsaI);
    r.fragments = vec![r.fragments[0].clone(); 33];
    assert!(d::source_preflight(&[(&inputs[0].id, &inputs[0].molecule)], &target, &r).is_err());
    let (mut inputs, target, r) = fixture(false, d::GoldenGateEnzyme::BsaI);
    inputs[0]
        .molecule
        .metadata
        .insert("duplex_model".into(), "Reference only".into());
    inputs[0].revision = inputs[0].molecule.revision();
    assert!(d::design(&inputs, &target, &r).is_err());
}
