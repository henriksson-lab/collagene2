//! Independent hand-specified strand diagrams. Expected products are not derived by digest().
use collagene_core::{analysis::Primer, chemistry::*, verification::*, *};

fn duplex(top: &str, bottom: &str, start: i64) -> Duplex {
    Duplex {
        top: Sequence::dna(top).unwrap(),
        bottom: Sequence::dna(bottom).unwrap(),
        bottom_start: start,
        topology: Topology::Linear,
        left: Terminus::CLEAVED,
        right: Terminus::CLEAVED,
        nicks: vec![],
    }
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}

#[test]
fn hand_drawn_nonpalindromic_five_and_three_prime_junctions_preserve_both_strands() {
    // 5′ case, top/reference coordinates:
    // A top:    GGGA        B top:        CTTT
    // A bottom: CCCTGA      B bottom:       AA       (bottom diagrams are 3′→5′)
    // A right overhang is 5′-AG; B left overhang is 5′-CT.
    let a = duplex("GGGA", "AGTCCC", 0);
    let b = duplex("CTTT", "AA", 2);
    let product = ligate(&[a.clone(), b.clone()], false, false).unwrap();
    assert_eq!(product.top.as_str(), "GGGACTTT");
    assert_eq!(product.bottom.as_str(), "AAAGTCCC");
    assert_eq!(product.bottom_start, 0);
    assert!(product.nicks.is_empty());
    assert_eq!(a.ends().unwrap().1.overhang, "AG");
    assert_eq!(b.ends().unwrap().0.overhang, "CT");
    // 3′ case: the same final duplex, but top A extends through CT and bottom B
    // extends left through its complementary 3′-AG. No reference bases are removed.
    let mut a = duplex("GGGACT", "TCCC", 0);
    let b = duplex("TT", "AAAG", -2);
    assert_eq!(a.ends().unwrap().1.polarity, Polarity::ThreePrime);
    assert_eq!(b.ends().unwrap().0.overhang, "AG");
    a.right.five_prime = FivePrime::Hydroxyl;
    let product = ligate(&[a, b], false, true).unwrap();
    assert_eq!(product.top.as_str(), "GGGACTTT");
    assert_eq!(product.bottom.as_str(), "AAAGTCCC");
    // Only the bottom junction lacks phosphate. Its boundary is 4, not top boundary 6.
    assert_eq!(
        product.nicks,
        vec![Nick {
            strand: Strand::Bottom,
            position: 4
        }]
    );
}

#[test]
fn hand_drawn_asymmetric_closure_has_strand_specific_origin_nicks() {
    // Top is CTGGGA (6 nt). Bottom in forward reference coordinates covers [2,8)
    // as GGGACT, making complementary CT/AG 5′ overhangs at the two ends.
    let mut fragment = duplex("CTGGGA", "AGTCCC", 2);
    fragment.left.five_prime = FivePrime::Hydroxyl;
    let circle = ligate(&[fragment.clone()], true, true).unwrap();
    assert_eq!(circle.top.as_str(), "CTGGGA");
    assert_eq!(circle.bottom.as_str(), "TCCCAG");
    assert_eq!(
        circle.nicks,
        vec![Nick {
            strand: Strand::Top,
            position: 0
        }]
    );
    fragment.left = Terminus::CLEAVED;
    fragment.right.five_prime = FivePrime::Hydroxyl;
    let circle = ligate(&[fragment], true, true).unwrap();
    assert_eq!(
        circle.nicks,
        vec![Nick {
            strand: Strand::Bottom,
            position: 2
        }]
    );
    // Equivalent 3′ overhang diagram: bottom forward coverage [-2,4) is CTGGGA.
    let mut three = duplex("GGGACT", "TCCCAG", -2);
    three.right.five_prime = FivePrime::Hydroxyl;
    let circle = ligate(&[three], true, true).unwrap();
    assert_eq!(circle.bottom.as_str(), "AGTCCC");
    assert_eq!(
        circle.nicks,
        vec![Nick {
            strand: Strand::Bottom,
            position: 4
        }]
    );
}

#[test]
fn reverse_maps_asymmetric_ends_and_nicks_using_actual_strand_extent() {
    let mut d = duplex("AACCGGTT", "GTAACCGG", 2);
    d.left = Terminus::PCR;
    d.right = Terminus::CLEAVED;
    d.nicks = vec![
        Nick {
            strand: Strand::Top,
            position: 3,
        },
        Nick {
            strand: Strand::Bottom,
            position: 5,
        },
    ];
    d.validate().unwrap();
    let r = d.reverse().unwrap();
    assert_eq!(r.top.as_str(), "GTAACCGG");
    assert_eq!(r.bottom.as_str(), "AACCGGTT");
    assert_eq!(r.bottom_start, 2);
    assert_eq!(r.left, Terminus::CLEAVED);
    assert_eq!(r.right, Terminus::PCR);
    assert_eq!(
        r.nicks,
        vec![
            Nick {
                strand: Strand::Bottom,
                position: 7
            },
            Nick {
                strand: Strand::Top,
                position: 5
            }
        ]
    );
    assert_eq!(r.ends().unwrap().0.overhang, "GT");
    assert_eq!(r.ends().unwrap().1.overhang, "AA");
}

#[test]
fn public_result_apis_reject_extreme_coordinates_without_panics() {
    for coordinate in [i64::MIN, i64::MAX, -10_000_001, 10_000_001] {
        let mut d = duplex("ACGT", "ACGT", 0);
        d.bottom_start = coordinate;
        assert!(std::panic::catch_unwind(|| d.validate()).unwrap().is_err());
        assert!(std::panic::catch_unwind(|| d.ends()).unwrap().is_err());
        assert!(std::panic::catch_unwind(|| d.reverse()).unwrap().is_err());
        assert!(
            std::panic::catch_unwind(|| d.treat(EndTreatment::Phosphorylate))
                .unwrap()
                .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| ligate(&[d.clone()], true, true))
                .unwrap()
                .is_err()
        );
    }
}

#[test]
fn staggered_physical_pcr_templates_need_a_two_strand_model_not_a_top_projection() {
    let source = Molecule::dna("BsaI source", "AAAGGTCTCTAGTCCCCC", Topology::Linear).unwrap();
    let target = Molecule::dna("PCR target", "AAAGGTCTCTAGTC", Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "PCR from staggered digest".into(),
        target_revision: target.revision(),
        inputs: vec![input("source", source)],
        steps: vec![
            Operation::Digest {
                template: Source::Input("source".into()),
                enzymes: vec!["BsaI".into()],
                fragment_index: 0,
            },
            Operation::Pcr {
                template: Source::Step(0),
                forward: Primer::new("F", "AAA", "").unwrap(),
                reverse: Primer::new("R", "GACT", "").unwrap(),
            },
        ],
    };
    let report = verify(&plan, &target);
    // The longer bottom strand can seed F extension to make a 14-nt top template.
    // This is not equivalent to PCR on the 10-nt top-only projection. Until modeled,
    // the verdict must be incomplete instead of incorrectly proving no product.
    assert_eq!(report.products[0].sequence.as_str(), "AAAGGTCTCT");
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
}

#[test]
fn saved_physical_projection_cannot_bypass_overlap_end_state_checks() {
    let physical = duplex("GGGA", "AGTCCC", 0)
        .molecule("saved staggered fragment")
        .unwrap();
    let other = Molecule::dna("other", "GAAT", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "GGGAAT", Topology::Linear).unwrap();
    let mut plan = DesignPlan {
        name: "No chemistry laundering through native projection".into(),
        target_revision: target.revision(),
        inputs: vec![input("physical", physical), input("other", other)],
        steps: vec![Operation::Overlap {
            fragments: vec![
                Source::Input("physical".into()),
                Source::Input("other".into()),
            ],
            overlaps: vec![2],
            circular: false,
        }],
    };
    assert_eq!(verify(&plan, &target).status, Status::Incomplete);
    // Native/imported known-but-unmapped chemistry gets the same treatment.
    plan.inputs[0].molecule.metadata.clear();
    plan.inputs[0]
        .molecule
        .metadata
        .insert("source:snapgene:chemistry_unmapped".into(), "true".into());
    plan.inputs[0].revision = plan.inputs[0].molecule.revision();
    assert_eq!(verify(&plan, &target).status, Status::Incomplete);
}

#[test]
fn invalid_or_duplicated_physical_states_cannot_publish_positive_evidence() {
    let mut invalid = duplex("ACGT", "AAAA", 0);
    invalid.topology = Topology::Circular;
    assert!(invalid.validate().is_err());
    assert!(invalid.summary().is_err());
    let mut duplicate = duplex("ACGT", "ACGT", 0);
    duplicate.nicks = vec![
        Nick {
            strand: Strand::Top,
            position: 1,
        },
        Nick {
            strand: Strand::Top,
            position: 1,
        },
    ];
    assert!(duplicate.validate().is_err());
    assert!(duplicate.summary().is_err());
}

#[test]
fn aggregate_nick_budget_is_explicit_instead_of_dropping_bonds_from_the_model() {
    let mut substrate = duplex(&"A".repeat(10_002), &"T".repeat(10_002), 0);
    substrate.nicks = (1..=10_000)
        .map(|position| Nick {
            strand: Strand::Top,
            position,
        })
        .collect();
    substrate.validate().unwrap();
    assert!(matches!(
        ligate(&[substrate.clone(), substrate], false, true),
        Err(Error::TooLarge)
    ));
}
