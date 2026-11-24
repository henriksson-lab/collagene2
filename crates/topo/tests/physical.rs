use collagene_model::Sequence;
use collagene_topo::*;
fn strand(id: &str, dna: &str, charged: bool) -> Strand {
    Strand {
        role: "literal-input".into(),
        id: id.into(),
        sequence: Sequence::dna(dna).unwrap(),
        five_prime: if charged {
            FivePrime::Phosphate
        } else {
            FivePrime::Hydroxyl
        },
        three_prime: if charged {
            ThreePrime::Phosphotyrosyl {
                enzyme: Enzyme::VacciniaTypeIBTyr274,
                activation_id: format!("enzyme-{id}"),
            }
        } else {
            ThreePrime::Hydroxyl
        },
    }
}
fn pair(a: usize, x: usize, b: usize, y: usize, n: usize) -> PairSpan {
    PairSpan {
        strand_a: a,
        start_a: x,
        strand_b: b,
        start_b: y,
        length: n,
    }
}
fn ta() -> Substrates {
    Substrates {
        strands: vec![
            strand("vector-top", "AGGGGCATCCCTT", true),
            strand("vector-bottom", "AGGGATGCCCCTT", true),
            strand("insert-top", "CGAATGA", false),
            strand("insert-bottom", "CATTCGA", false),
        ],
        pairing: vec![pair(0, 0, 1, 0, 12), pair(2, 0, 3, 0, 6)],
    }
}
fn directional() -> Substrates {
    Substrates {
        strands: vec![
            strand("vector-top", "AAGGGGCATCCCTT", true),
            strand("vector-bottom", "GGTGAAGGGATGCCCCTT", true),
            strand("insert-top", "CACCATGACG", false),
            strand("insert-bottom", "CGTCATGGTG", false),
        ],
        pairing: vec![pair(0, 0, 1, 4, 14), pair(2, 0, 3, 0, 10)],
    }
}
fn transfers(geometry: Geometry) -> [Transfer; 2] {
    [
        Transfer {
            donor: 0,
            acceptor: 2,
            geometry,
        },
        Transfer {
            donor: 1,
            acceptor: 3,
            geometry: if geometry == Geometry::DirectionalCacc {
                Geometry::Flush
            } else {
                geometry
            },
        },
    ]
}
fn bases(r: &Report) -> Vec<&str> {
    r.product
        .chains
        .iter()
        .map(|c| c.sequence.as_str())
        .collect()
}
#[test]
fn ta_both_strands_two_bonds_two_unsealed_nicks() {
    let s = ta();
    let t = transfers(Geometry::ThreePrimeTA);
    let r = simulate(&s, &t, &Limits::default()).unwrap();
    assert!(r.all_requested_transfers_applied);
    assert_eq!(
        bases(&r),
        vec!["AGGGGCATCCCTTCGAATGA", "AGGGATGCCCCTTCATTCGA"]
    );
    assert_eq!(r.source_bases, 40);
    assert_eq!(
        r.product.released_activations,
        vec!["enzyme-vector-top", "enzyme-vector-bottom"]
    );
    assert_eq!(r.product.nicks.len(), 2);
    assert!(r.product.unpaired.is_empty());
    assert_eq!(r.product.components.len(), 1);
    assert!(!r.product.components[0].covalently_closed);
    for c in &r.product.chains {
        assert!(matches!(
            c.ends,
            ChainEnds::Linear {
                five_prime: FivePrime::Phosphate,
                three_prime: ThreePrime::Hydroxyl
            }
        ));
    }
    audit_sources(&s, &r).unwrap();
    replay(&s, &t, &Limits::default(), &r).unwrap();
}
#[test]
fn ta_orientation_alternatives_keep_actual_sources() {
    let s = ta();
    let r = orientations(&s, [0, 1], [2, 3], false, &Limits::default()).unwrap();
    assert_eq!(r.candidates.len(), 2);
    assert!(r.orientation_search_complete);
    assert!(!r.biochemical_mixture_exhaustive);
    assert!(
        r.candidates
            .iter()
            .all(|r| r.all_requested_transfers_applied)
    );
    assert_eq!(
        bases(&r.candidates[1]),
        vec!["AGGGGCATCCCTTCATTCGA", "AGGGATGCCCCTTCGAATGA"]
    );
    for c in &r.candidates {
        audit_sources(&s, c).unwrap();
    }
}
#[test]
fn directional_four_actual_bottom_calls_are_retained_as_attached_flap() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let r = simulate(&s, &t, &Limits::default()).unwrap();
    assert!(r.all_requested_transfers_applied, "{:?}", r.steps);
    assert_eq!(
        bases(&r),
        vec!["AAGGGGCATCCCTTCACCATGACG", "GGTGAAGGGATGCCCCTTCGTCATGGTG"]
    );
    assert_eq!(r.source_bases, 52);
    assert_eq!(r.product.nicks.len(), 1);
    assert_eq!(r.product.unpaired.len(), 1);
    let flap = &r.product.unpaired[0];
    assert_eq!((flap.chain, flap.start, flap.end), (1, 24, 28));
    assert!(flap.terminal_three);
    assert_eq!(flap.displaced_by, vec![0]);
    assert_eq!(
        &r.product.chains[1].sequence.as_str()[flap.start..flap.end],
        "GGTG"
    );
    assert_eq!(r.steps[0].removed_pairs.len(), 4);
    assert_eq!(r.steps[0].added_pairs.len(), 4);
    assert!(r.product.components.iter().all(|c| !c.covalently_closed));
    audit_sources(&s, &r).unwrap();
}
#[test]
fn directional_reverse_is_unsupported_not_claimed_absent_and_partial_bonds_survive() {
    let r = orientations(&directional(), [0, 1], [2, 3], true, &Limits::default()).unwrap();
    assert!(r.candidates[0].all_requested_transfers_applied);
    let reverse = &r.candidates[1];
    assert_eq!(reverse.steps[0].outcome, Outcome::Unsupported);
    assert_eq!(reverse.steps[1].outcome, Outcome::Applied);
    assert!(reverse.unsupported_or_unknown);
    assert_eq!(
        reverse.product.released_activations,
        vec!["enzyme-vector-bottom"]
    );
    assert!(reverse.product.chains.iter().any(|c| matches!(
        c.ends,
        ChainEnds::Linear {
            three_prime: ThreePrime::Phosphotyrosyl { .. },
            ..
        }
    )));
    audit_sources(&directional(), reverse).unwrap();
}
#[test]
fn single_junction_does_not_seal_opposite_backbone_or_drop_other_overhang() {
    let s = Substrates {
        strands: vec![
            strand("donor", "GACGCCCTT", true),
            strand("donor-opposite", "AGGGCGTC", false),
            strand("acceptor", "CGTACGGA", false),
            strand("acceptor-opposite", "CCGTACGA", false),
        ],
        pairing: vec![pair(0, 0, 1, 0, 8), pair(2, 0, 3, 0, 7)],
    };
    let r = simulate(
        &s,
        &[Transfer {
            donor: 0,
            acceptor: 2,
            geometry: Geometry::ThreePrimeTA,
        }],
        &Limits::default(),
    )
    .unwrap();
    assert!(r.all_requested_transfers_applied);
    assert_eq!(bases(&r), vec!["GACGCCCTTCGTACGGA", "AGGGCGTC", "CCGTACGA"]);
    assert_eq!(r.product.nicks.len(), 1);
    assert_eq!(
        r.product
            .unpaired
            .iter()
            .map(|u| u.end - u.start)
            .sum::<usize>(),
        1
    );
    assert_eq!(r.source_bases, 33);
}
#[test]
fn five_prime_phosphate_blocks_transfer_and_unknown_chemistry_stays_unknown() {
    for five in [FivePrime::Phosphate, FivePrime::Blocked, FivePrime::Unknown] {
        let mut s = ta();
        s.strands[2].five_prime = five;
        let r = simulate(
            &s,
            &transfers(Geometry::ThreePrimeTA)[..1],
            &Limits::default(),
        )
        .unwrap();
        assert_eq!(
            r.steps[0].outcome,
            if five == FivePrime::Unknown {
                Outcome::Unknown
            } else {
                Outcome::BlockedUnderDeclaredModel
            }
        );
        assert_eq!(r.product.chains.len(), 4);
        assert!(r.product.released_activations.is_empty());
        audit_sources(&s, &r).unwrap();
    }
}
#[test]
fn motif_or_sequence_text_never_invents_activation() {
    let mut s = ta();
    s.strands[0].three_prime = ThreePrime::Hydroxyl;
    let r = simulate(
        &s,
        &transfers(Geometry::ThreePrimeTA)[..1],
        &Limits::default(),
    )
    .unwrap();
    assert_eq!(r.steps[0].outcome, Outcome::BlockedUnderDeclaredModel);
    assert!(r.product.released_activations.is_empty());
    let mut s = ta();
    s.strands[0].three_prime = ThreePrime::Phosphotyrosyl {
        enzyme: Enzyme::Uncharacterized,
        activation_id: "unknown".into(),
    };
    assert_eq!(
        simulate(
            &s,
            &transfers(Geometry::ThreePrimeTA)[..1],
            &Limits::default()
        )
        .unwrap()
        .steps[0]
            .outcome,
        Outcome::Unsupported
    );
}
#[test]
fn impossible_pairing_duplicates_and_reused_activation_are_invalid() {
    let mut s = ta();
    s.pairing.push(s.pairing[0].clone());
    assert!(simulate(&s, &[], &Limits::default()).is_err());
    let mut s = ta();
    s.strands[1].sequence = Sequence::dna("AGGGATGCCCCTC").unwrap();
    assert!(simulate(&s, &[], &Limits::default()).is_ok());
    s.strands[1].sequence = Sequence::dna("CGGGATGCCCCTT").unwrap();
    assert!(simulate(&s, &[], &Limits::default()).is_err());
    let mut s = ta();
    s.strands[1].three_prime = s.strands[0].three_prime.clone();
    assert!(preflight(&s, &[], &Limits::default()).is_err());
}
#[test]
fn source_target_or_claim_mutation_cannot_replay() {
    let s = ta();
    let t = transfers(Geometry::ThreePrimeTA);
    let r = simulate(&s, &t, &Limits::default()).unwrap();
    let mut forged = r.clone();
    forged.product.chains[0].sequence = Sequence::dna("AAGGGCATCCCTTCGAATGA").unwrap();
    assert!(audit_sources(&s, &forged).is_err());
    let mut forged = r.clone();
    forged.product.nicks.clear();
    assert!(replay(&s, &t, &Limits::default(), &forged).is_err());
    let mut changed = s.clone();
    changed.strands[0].role = "other source".into();
    assert!(replay(&changed, &t, &Limits::default(), &r).is_err());
}
#[test]
fn all_declared_limits_fail_explicitly_and_never_become_absence() {
    let s = ta();
    let t = transfers(Geometry::ThreePrimeTA);
    for mode in 0..8 {
        let mut l = Limits::default();
        match mode {
            0 => l.max_source_bases = 39,
            1 => l.max_strands = 3,
            2 => l.max_pair_spans = 1,
            3 => l.max_steps = 1,
            4 => l.max_work_cells = 1,
            5 => l.max_work_bytes = 1,
            6 => l.max_output_bases = 39,
            _ => l.max_report_bytes = 1,
        };
        assert!(
            matches!(simulate(&s, &t, &l), Err(Error::Limit { .. })),
            "limit mode{mode}"
        );
    }
    let l = Limits {
        max_output_bases: 79,
        ..Limits::default()
    };
    assert!(orientations(&s, [0, 1], [2, 3], false, &l).is_err());
}

#[test]
fn templated_nick_reseal_can_close_one_activated_strand_without_inventing_opposite_repair() {
    let mut top = strand("nicked-ring-top", "GACGCCCTT", true);
    top.five_prime = FivePrime::Hydroxyl;
    let s = Substrates {
        strands: vec![top, strand("opposite", "CGTCAAGGG", false)],
        pairing: vec![pair(0, 0, 1, 0, 4), pair(0, 4, 1, 4, 5)],
    };
    let before = simulate(&s, &[], &Limits::default()).unwrap();
    assert_eq!(before.product.nicks.len(), 2);
    let r = simulate(
        &s,
        &[Transfer {
            donor: 0,
            acceptor: 0,
            geometry: Geometry::NickReseal,
        }],
        &Limits::default(),
    )
    .unwrap();
    assert!(r.all_requested_transfers_applied);
    let circular = r
        .product
        .chains
        .iter()
        .find(|c| matches!(c.ends, ChainEnds::CovalentlyCircular))
        .unwrap();
    assert_eq!(circular.sequence.as_str(), "GACGCCCTT");
    assert!(
        r.product
            .chains
            .iter()
            .any(|c| c.sequence.as_str() == "CGTCAAGGG"
                && matches!(c.ends, ChainEnds::Linear { .. }))
    );
    assert_eq!(r.product.nicks.len(), 1);
    assert!(!r.product.components[0].covalently_closed);
    assert_eq!(r.source_bases, 18);
    audit_sources(&s, &r).unwrap();
}
#[test]
fn exact_templated_nick_ligation_retains_all_three_source_strands() {
    let s = Substrates {
        strands: vec![
            strand("donor", "GACGCCCTT", true),
            strand("acceptor", "CGTACGGA", false),
            strand("continuous-template", "TCCGTACGAAGGGCGTC", false),
        ],
        pairing: vec![pair(0, 0, 2, 8, 9), pair(1, 0, 2, 0, 8)],
    };
    let r = simulate(
        &s,
        &[Transfer {
            donor: 0,
            acceptor: 1,
            geometry: Geometry::NickReseal,
        }],
        &Limits::default(),
    )
    .unwrap();
    assert!(r.all_requested_transfers_applied);
    assert_eq!(bases(&r), vec!["GACGCCCTTCGTACGGA", "TCCGTACGAAGGGCGTC"]);
    assert!(r.product.nicks.is_empty());
    assert_eq!(r.source_bases, 34);
    assert!(
        r.product
            .chains
            .iter()
            .all(|c| matches!(c.ends, ChainEnds::Linear { .. }))
    );
}
#[test]
fn displaced_flap_survives_independent_transfer_order_and_unreacted_species() {
    let mut s = directional();
    s.strands
        .push(strand("unreacted-free-oligo", "TACGA", false));
    let mut t = transfers(Geometry::DirectionalCacc);
    t.swap(0, 1);
    let r = simulate(&s, &t, &Limits::default()).unwrap();
    assert!(r.all_requested_transfers_applied);
    assert_eq!(
        bases(&r),
        vec![
            "AAGGGGCATCCCTTCACCATGACG",
            "GGTGAAGGGATGCCCCTTCGTCATGGTG",
            "TACGA"
        ]
    );
    assert_eq!(r.source_bases, 57);
    assert_eq!(r.product.components.len(), 2);
    assert!(
        r.product
            .unpaired
            .iter()
            .any(|u| u.displaced_by == vec![1] && u.end - u.start == 4)
    );
    assert!(
        r.product
            .unpaired
            .iter()
            .any(|u| u.chain == 2 && u.start == 0 && u.end == 5)
    );
    audit_sources(&s, &r).unwrap();
}
#[test]
fn ambiguous_pairing_is_explicit_unknown_but_impossible_iupac_pairing_is_invalid() {
    let mut s = ta();
    s.strands[0].sequence = Sequence::dna("AGGGNCATCCCTT").unwrap();
    let r = simulate(&s, &transfers(Geometry::ThreePrimeTA), &Limits::default()).unwrap();
    assert!(!r.fully_specified_dna);
    assert!(r.steps.iter().all(|e| e.outcome == Outcome::Unknown));
    assert_eq!(r.product.chains.len(), 4);
    assert_eq!(r.product.chains[0].sequence.as_str(), "AGGGNCATCCCTT");
    audit_sources(&s, &r).unwrap();
    s.strands[0].sequence = Sequence::dna("RGGGGCATCCCTT").unwrap(); // R cannot pair with the declared opposite terminal T? R includes A, so possible.
    assert!(simulate(&s, &[], &Limits::default()).is_ok());
    s.strands[0].sequence = Sequence::dna("YGGGGCATCCCTT").unwrap();
    assert!(simulate(&s, &[], &Limits::default()).is_err());
}
#[test]
fn actual_strand_permutations_preserve_products_without_reverse_complement_synthesis() {
    let old = ta();
    let permutation = [3usize, 1, 0, 2];
    let mut inverse = [0; 4];
    for (n, &p) in permutation.iter().enumerate() {
        inverse[p] = n;
    }
    let s = Substrates {
        strands: permutation
            .iter()
            .map(|&p| old.strands[p].clone())
            .collect(),
        pairing: old
            .pairing
            .iter()
            .map(|p| PairSpan {
                strand_a: inverse[p.strand_a],
                strand_b: inverse[p.strand_b],
                ..p.clone()
            })
            .collect(),
    };
    let t = transfers(Geometry::ThreePrimeTA).map(|t| Transfer {
        donor: inverse[t.donor],
        acceptor: inverse[t.acceptor],
        ..t
    });
    let r = simulate(&s, &t, &Limits::default()).unwrap();
    assert!(r.all_requested_transfers_applied);
    let mut actual = bases(&r);
    actual.sort();
    let mut expected = vec!["AGGGGCATCCCTTCGAATGA", "AGGGATGCCCCTTCATTCGA"];
    expected.sort();
    assert_eq!(actual, expected);
    assert_eq!(r.product.nicks.len(), 2);
    audit_sources(&s, &r).unwrap();
}
#[test]
fn consumed_activation_cannot_form_a_second_bond() {
    let s = ta();
    let step = transfers(Geometry::ThreePrimeTA)[0].clone();
    let r = simulate(&s, &[step.clone(), step], &Limits::default()).unwrap();
    assert_eq!(r.steps[0].outcome, Outcome::Applied);
    assert_eq!(r.steps[1].outcome, Outcome::BlockedUnderDeclaredModel);
    assert_eq!(r.product.released_activations.len(), 1);
    audit_sources(&s, &r).unwrap();
}

#[test]
fn complete_four_base_displacement_is_a_separate_retained_strand_not_an_attached_flap() {
    let s = Substrates {
        strands: vec![
            strand("donor", "GACGCCCTT", true),
            strand("opposite", "GGTGAAGGGCGTC", false),
            strand("incoming", "CACC", false),
            strand("displaced", "GGTG", false),
        ],
        pairing: vec![pair(0, 0, 1, 4, 9), pair(2, 0, 3, 0, 4)],
    };
    let r = simulate(
        &s,
        &[Transfer {
            donor: 0,
            acceptor: 2,
            geometry: Geometry::DirectionalCacc,
        }],
        &Limits::default(),
    )
    .unwrap();
    assert!(r.all_requested_transfers_applied);
    assert_eq!(bases(&r), vec!["GACGCCCTTCACC", "GGTGAAGGGCGTC", "GGTG"]);
    assert_eq!(r.source_bases, 30);
    assert_eq!(r.product.components.len(), 2);
    assert!(r.product.nicks.is_empty());
    let d = r.steps[0].displacement.as_ref().unwrap();
    assert_eq!(d.kind, DisplacementKind::FreeDisplacedChain);
    assert_eq!(d.five_prime_continuation, None);
    assert_eq!(d.continuation_pair, None);
    assert_eq!(
        d.calls_five_to_three,
        (0..4)
            .map(|index| Call { strand: 3, index })
            .collect::<Vec<_>>()
    );
    let free = &r.product.chains[2];
    assert!(matches!(
        free.ends,
        ChainEnds::Linear {
            five_prime: FivePrime::Hydroxyl,
            three_prime: ThreePrime::Hydroxyl
        }
    ));
    assert_eq!(
        free.sources,
        vec![SourcePlacement {
            source_strand: 3,
            start: 0,
            end: 4,
            output_start: 0
        }]
    );
    let unpaired = r.product.unpaired.iter().find(|u| u.chain == 2).unwrap();
    assert!(unpaired.terminal_five && unpaired.terminal_three);
    assert_eq!((unpaired.start, unpaired.end), (0, 4));
    assert_eq!(unpaired.displaced_by, vec![0]);
    audit_sources(&s, &r).unwrap();
}
#[test]
fn unpaired_continuation_is_not_mislabeled_a_duplex_anchored_flap() {
    let s = Substrates {
        strands: vec![
            strand("donor", "GACGCCCTT", true),
            strand("opposite", "GGTGAAGGGCGTC", false),
            strand("incoming", "CACC", false),
            strand("tail", "AGGTG", false),
        ],
        pairing: vec![pair(0, 0, 1, 4, 9), pair(2, 0, 3, 1, 4)],
    };
    let r = simulate(
        &s,
        &[Transfer {
            donor: 0,
            acceptor: 2,
            geometry: Geometry::DirectionalCacc,
        }],
        &Limits::default(),
    )
    .unwrap();
    let d = r.steps[0].displacement.as_ref().unwrap();
    assert_eq!(d.kind, DisplacementKind::UnpairedTailExtension);
    assert_eq!(
        d.five_prime_continuation,
        Some(Call {
            strand: 3,
            index: 0
        })
    );
    assert_eq!(d.continuation_pair, None);
    assert_eq!(r.product.chains[2].sequence.as_str(), "AGGTG");
    assert_eq!(r.product.components.len(), 2);
    assert_eq!(r.source_bases, 31);
    audit_sources(&s, &r).unwrap();
    let anchored = simulate(
        &directional(),
        &transfers(Geometry::DirectionalCacc),
        &Limits::default(),
    )
    .unwrap();
    let a = anchored.steps[0].displacement.as_ref().unwrap();
    assert_eq!(a.kind, DisplacementKind::DuplexAnchoredThreePrimeFlap);
    assert_eq!(
        a.five_prime_continuation,
        Some(Call {
            strand: 3,
            index: 5
        })
    );
    assert_eq!(
        a.continuation_pair,
        Some(Call {
            strand: 2,
            index: 4
        })
    );
    assert_eq!(
        a.calls_five_to_three,
        (6..10)
            .map(|index| Call { strand: 3, index })
            .collect::<Vec<_>>()
    );
}
