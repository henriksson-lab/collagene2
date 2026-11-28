use collagene_core::{analysis::*, chemistry::*, verification::*};
use collagene_model::{Error, Molecule, Sequence, Topology};
fn dna(bases: &str, topology: Topology) -> Duplex {
    Duplex::from_molecule(&Molecule::dna("fixture", bases, topology).unwrap()).unwrap()
}
fn enzyme(name: &str) -> Enzyme {
    starter_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .unwrap()
}
fn cut(bases: &str, name: &str, topology: Topology) -> Vec<DigestFragment> {
    digest(&dna(bases, topology), &[enzyme(name)]).unwrap()
}

#[test]
fn independent_five_prime_and_three_prime_strand_fixtures() {
    let eco = cut("TTTGAATTCAAA", "EcoRI", Topology::Linear);
    assert_eq!(eco.len(), 2);
    assert_eq!(
        (
            eco[0].duplex.top.as_str(),
            eco[0].duplex.bottom.as_str(),
            eco[0].duplex.bottom_start
        ),
        ("TTTG", "AATTCAAA", 0)
    );
    assert_eq!(
        (
            eco[1].duplex.top.as_str(),
            eco[1].duplex.bottom.as_str(),
            eco[1].duplex.bottom_start
        ),
        ("AATTCAAA", "TTTG", 4)
    );
    let right = eco[0].duplex.ends().unwrap().1;
    let left = eco[1].duplex.ends().unwrap().0;
    assert_eq!(
        (right.polarity, right.overhang.as_str()),
        (Polarity::FivePrime, "AATT")
    );
    assert_eq!(left, right);
    assert_eq!(right.chemistry, Terminus::CLEAVED);
    assert_eq!(eco[0].duplex.left, Terminus::UNKNOWN);
    let psti = cut("AAACTGCAGTTT", "PstI", Topology::Linear);
    assert_eq!(
        (psti[0].duplex.top.as_str(), psti[0].duplex.bottom.as_str()),
        ("AAACTGCA", "GTTT")
    );
    assert_eq!(
        (
            psti[1].duplex.top.as_str(),
            psti[1].duplex.bottom.as_str(),
            psti[1].duplex.bottom_start
        ),
        ("GTTT", "AAACTGCA", -4)
    );
    assert_eq!(
        psti[0].duplex.ends().unwrap().1.polarity,
        Polarity::ThreePrime
    );
    assert_eq!(psti[1].duplex.ends().unwrap().0.overhang, "TGCA");
    let blunt = cut("AAAGATATCTTT", "EcoRV", Topology::Linear);
    assert_eq!(blunt[0].duplex.top.as_str(), "AAAGAT");
    assert_eq!(blunt[0].duplex.bottom.as_str(), "ATCTTT");
    assert_eq!(blunt[0].duplex.ends().unwrap().1.polarity, Polarity::Blunt);
}

#[test]
fn type_iis_nonpalindromic_overhangs_and_reverse_sites() {
    let bases = "AAAGGTCTCTAGTCCCCC";
    let fragments = cut(bases, "BsaI", Topology::Linear);
    assert_eq!(fragments[0].duplex.top.len(), 10);
    assert_eq!(fragments[0].duplex.bottom.len(), 14);
    assert_eq!(fragments[0].duplex.ends().unwrap().1.overhang, "GACT");
    assert_eq!(fragments[1].duplex.ends().unwrap().0.overhang, "AGTC");
    let source = dna(bases, Topology::Linear);
    let reverse = source.reverse().unwrap();
    let fragments = digest(&reverse, &[enzyme("BsaI")]).unwrap();
    assert_eq!(fragments[0].duplex.top.len(), 4);
    assert_eq!(fragments[0].duplex.bottom.len(), 8);
    assert_eq!(fragments[0].duplex.ends().unwrap().1.overhang, "AGTC");
    assert_eq!(fragments[1].duplex.ends().unwrap().0.overhang, "GACT");
    assert_eq!(
        ligate(
            &fragments
                .iter()
                .map(|f| f.duplex.clone())
                .collect::<Vec<_>>(),
            false,
            false
        )
        .unwrap(),
        reverse
    );
}

#[test]
fn circular_origin_cut_rotates_both_strands_without_losing_bases() {
    let source = dna("AATTCAAAAG", Topology::Circular);
    let fragments = digest(&source, &[enzyme("EcoRI")]).unwrap();
    assert_eq!(fragments.len(), 1);
    let fragment = &fragments[0];
    assert_eq!((fragment.top_start, fragment.bottom_start), (0, 4));
    assert_eq!(fragment.duplex.top.as_str(), "AATTCAAAAG");
    assert_eq!(fragment.duplex.bottom.as_str(), "AATTCTTTTG");
    let product = ligate(&[fragment.duplex.clone()], true, false).unwrap();
    assert_eq!(product.top, source.top);
    assert_eq!(product.bottom, source.bottom);
    assert_eq!(product.bottom_start, 0);
    assert!(product.nicks.is_empty());
}

#[test]
fn digest_ligation_round_trip_for_every_origin_and_both_strand_orientations() {
    for (name, bases) in [
        ("EcoRI", "AAAAGAATTCTTTTGAATTCCCC"),
        ("PstI", "AAACTGCAGTTTCTGCAGCCC"),
        ("EcoRV", "AAAGATATCTTTGATATCCCC"),
        ("BsaI", "AAAGGTCTCTAGTCCCCC"),
    ] {
        for offset in 0..bases.len() {
            let rotated = format!("{}{}", &bases[offset..], &bases[..offset]);
            let original = dna(&rotated, Topology::Circular);
            for source in [original.clone(), original.reverse().unwrap()] {
                let fragments = digest(&source, &[enzyme(name)]).unwrap();
                let result = ligate(
                    &fragments
                        .iter()
                        .map(|f| f.duplex.clone())
                        .collect::<Vec<_>>(),
                    true,
                    false,
                )
                .unwrap();
                assert!(
                    rotation_offset(source.top.as_str(), result.top.as_str()).is_some(),
                    "{name} origin {offset}"
                );
                assert_eq!(result.bottom, result.top.reverse_complement().unwrap());
                assert_eq!(source.top.len(), result.top.len());
                assert!(result.nicks.is_empty());
                assert_eq!(result.reverse().unwrap().reverse().unwrap(), result);
            }
        }
    }
}

#[test]
fn phosphorylation_and_nicks_are_not_silently_repaired() {
    let vector = cut("AAAAGAATTCTTTT", "EcoRI", Topology::Circular)
        .remove(0)
        .duplex
        .treat(EndTreatment::Dephosphorylate)
        .unwrap();
    let insert = cut("CCCCGAATTCGGGG", "EcoRI", Topology::Circular)
        .remove(0)
        .duplex;
    assert!(ligate(&[vector.clone()], true, true).is_err());
    assert!(ligate(&[vector.clone(), insert.clone()], true, false).is_err());
    let product = ligate(&[vector.clone(), insert.clone()], true, true).unwrap();
    assert_eq!(
        product.nicks,
        vec![
            Nick {
                strand: Strand::Bottom,
                position: 18
            },
            Nick {
                strand: Strand::Top,
                position: 0
            }
        ]
    );
    assert_eq!(product.summary().unwrap().covalently_closed_in_model, false);
    assert!(digest(&product, &[enzyme("EcoRI")]).is_err());
    assert!(product.treat(EndTreatment::Phosphorylate).is_err());
    // A saved molecule is not an authoritative replacement for the originating chemistry graph.
    assert!(Duplex::from_molecule(&product.molecule("nicked").unwrap()).is_err());
    let repaired_terminals = vector.treat(EndTreatment::Phosphorylate).unwrap();
    let sealed = ligate(&[repaired_terminals, insert], true, false).unwrap();
    assert!(sealed.nicks.is_empty());
    let unknown = dna("ACGT", Topology::Linear);
    assert!(matches!(
        ligate(&[unknown.clone(), unknown], false, true),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn incompatible_polarity_noncomplementary_overhangs_and_unknown_chemistry_fail() {
    let five = cut("TTTGAATTCAAA", "EcoRI", Topology::Linear)
        .remove(0)
        .duplex;
    let three = cut("AAACTGCAGTTT", "PstI", Topology::Linear)
        .remove(1)
        .duplex;
    assert!(ligate(&[five, three], false, true).is_err());
    let a = cut("AAAGGTCTCTAGTCCCCC", "BsaI", Topology::Linear)
        .remove(0)
        .duplex;
    let b = cut("AAAGGTCTCTAAAACCCC", "BsaI", Topology::Linear)
        .remove(1)
        .duplex;
    assert!(ligate(&[a, b], false, true).is_err());
    let mut invalid = dna("ACGT", Topology::Linear);
    invalid.bottom = Sequence::dna("AAAA").unwrap();
    assert!(invalid.validate().is_err());
    assert!(invalid.reverse().is_err());
}

#[test]
fn unsupported_cut_geometry_stays_explicit_and_uncut_circles_stay_circular() {
    let uncut = digest(&dna("AAA", Topology::Circular), &[enzyme("EcoRI")]).unwrap();
    assert_eq!(uncut.len(), 1);
    assert_eq!(uncut[0].duplex.topology, Topology::Circular);
    assert!(ligate(&[uncut[0].duplex.clone()], true, false).is_err());
    let source = dna("AAAAGAATTCTTTT", Topology::Linear);
    let a = Enzyme {
        name: "A".into(),
        recognition: "GAATTC".into(),
        top_cut: 1,
        bottom_cut: 5,
    };
    let b = Enzyme {
        name: "B".into(),
        recognition: "GAATTC".into(),
        top_cut: 1,
        bottom_cut: 3,
    };
    assert!(matches!(
        digest(&source, &[a, b]),
        Err(Error::Unsupported(_))
    ));
    let source = dna("GAATTCGAATTC", Topology::Circular);
    let crossing = Enzyme {
        name: "crossing".into(),
        recognition: "GAATTC".into(),
        top_cut: 1,
        bottom_cut: 8,
    };
    assert!(matches!(
        digest(&source, &[crossing]),
        Err(Error::Unsupported(_))
    ));
    let at_end = Enzyme {
        name: "edge".into(),
        recognition: "GAATTC".into(),
        top_cut: 0,
        bottom_cut: 4,
    };
    assert!(matches!(
        digest(&dna("GAATTCTTT", Topology::Linear), &[at_end]),
        Err(Error::Unsupported(_))
    ));
}

fn plan_fixture() -> (DesignPlan, Molecule) {
    let vector = Molecule::dna("vector", "AAAAGAATTCTTTT", Topology::Circular).unwrap();
    let insert = Molecule::dna("insert", "CCCCGAATTCGGGG", Topology::Circular).unwrap();
    let target =
        Molecule::dna("target", "AATTCTTTTAAAAGAATTCGGGGCCCCG", Topology::Circular).unwrap();
    let plan = DesignPlan {
        name: "Restriction ligation".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "vector".into(),
                revision: vector.revision(),
                molecule: vector,
            },
            PlanInput {
                id: "insert".into(),
                revision: insert.revision(),
                molecule: insert,
            },
        ],
        steps: vec![
            Operation::Digest {
                template: Source::Input("vector".into()),
                enzymes: vec!["EcoRI".into()],
                fragment_index: 0,
            },
            Operation::Digest {
                template: Source::Input("insert".into()),
                enzymes: vec!["EcoRI".into()],
                fragment_index: 0,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                ],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    (plan, target)
}

#[test]
fn typed_plan_reconstructs_duplex_graph_and_retains_intermediates_on_failure() {
    let (mut plan, target) = plan_fixture();
    let report = verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.len(), 3);
    assert_eq!(report.products[2].sequence, target.sequence);
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.contains("0 retained nicks"))
    );
    let html = export_protocol_html(&plan, &target);
    assert!(
        html.contains("Complete restriction digest") && html.contains("Intended ordered ligation")
    );
    plan.steps[2] = Operation::Ligate {
        fragments: vec![
            OrientedSource {
                source: Source::Step(0),
                reverse: false,
            },
            OrientedSource {
                source: Source::Input("insert".into()),
                reverse: false,
            },
        ],
        circular: true,
        allow_nicked: false,
        options: Default::default(),
    };
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid, "{:?}", report.findings);
    assert_eq!(report.products.len(), 2);
    plan.steps[0] = Operation::Digest {
        template: Source::Input("vector".into()),
        enzymes: vec!["Unmodeled".into()],
        fragment_index: 0,
    };
    assert_eq!(verify(&plan, &target).status, Status::Incomplete);
}
