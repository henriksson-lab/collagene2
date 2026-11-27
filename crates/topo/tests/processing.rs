use collagene_model::Sequence;
use collagene_topo::{nick_sealing as nick, processing as p, *};
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
fn transfers(geometry: Geometry) -> [Step; 2] {
    [
        Step::Transfer(Transfer {
            donor: Call {
                strand: 0,
                index: if geometry == Geometry::DirectionalCacc {
                    13
                } else {
                    12
                },
            },
            acceptor: Call {
                strand: 2,
                index: 0,
            },
            geometry,
        }),
        Step::Transfer(Transfer {
            donor: Call {
                strand: 1,
                index: if geometry == Geometry::DirectionalCacc {
                    17
                } else {
                    12
                },
            },
            acceptor: Call {
                strand: 3,
                index: 0,
            },
            geometry: if geometry == Geometry::DirectionalCacc {
                Geometry::Flush
            } else {
                geometry
            },
        }),
    ]
}
fn call(strand: usize, index: usize) -> Call {
    Call { strand, index }
}
fn cut() -> p::Action {
    p::Action::CleaveThreePrimeFlap{bond:p::Bond{upstream:call(3,5),downstream:call(3,6)},displacement_step:0,profile:p::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,assumption:"Explicit ideal hydrolysis of this anchored flap junction; no named nuclease or host repair claim".into()}
}
fn seal(three: Call, five: Call) -> p::Action {
    p::Action::SealNick {
        selection: nick::NickSelection { three, five },
        ligase_assumption:
            "Explicit active ligase/cofactor hypothesis, only continuous-template3OH5P nick sealing"
                .into(),
    }
}
fn request(actions: Vec<p::Action>) -> Vec<Step> {
    actions.into_iter().map(Step::Process).collect()
}
fn full() -> Vec<Step> {
    request(vec![
        cut(),
        seal(call(2, 9), call(0, 0)),
        seal(call(3, 5), call(1, 0)),
    ])
}
fn base_count(product: &ProductGraph) -> usize {
    product.chains.iter().map(|c| c.sequence.len()).sum()
}
#[test]
fn literal_cut_partitions_insert_bottom_into_retained_six_and_free_four_calls() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let r = request(vec![cut()]);
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_eq!(simulate(&s, &t, &l).unwrap().product.chains.len(), 2);
    assert_eq!(report.product.chains.len(), 3);
    assert_eq!(simulate(&s, &t, &l).unwrap().product.nicks.len(), 1);
    assert_eq!(report.product.nicks.len(), 2);
    assert_eq!(report.product.components.len(), 2);
    assert_eq!(base_count(&report.product), 52);
    let free = report
        .product
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert_eq!(
        free.sources,
        vec![SourcePlacement {
            source_strand: 3,
            start: 6,
            end: 10,
            output_start: 0
        }]
    );
    assert_eq!(
        free.ends,
        ChainEnds::Linear {
            five_prime: FivePrime::Phosphate,
            three_prime: ThreePrime::Hydroxyl
        }
    );
    let bottom = report
        .product
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTGAAGGGATGCCCCTTCGTCAT")
        .unwrap();
    assert_eq!(
        bottom.sources[1],
        SourcePlacement {
            source_strand: 3,
            start: 0,
            end: 6,
            output_start: 18
        }
    );
    assert_eq!(
        bottom.ends,
        ChainEnds::Linear {
            five_prime: FivePrime::Phosphate,
            three_prime: ThreePrime::Hydroxyl
        }
    );
    assert_eq!(
        report.product.source_pairing,
        simulate(&s, &t, &l).unwrap().product.source_pairing
    );
    assert_eq!(
        report.product.released_activations,
        vec!["enzyme-vector-top", "enzyme-vector-bottom"]
    );
    let StepEvidence::Process(p::ActionEvidence::FlapCleavage(e)) = &report.steps[2] else {
        panic!()
    };
    assert_eq!(e.outcome, Outcome::Applied);
    assert_eq!(
        e.flap_calls,
        vec![call(3, 6), call(3, 7), call(3, 8), call(3, 9)]
    );
    assert_eq!(e.created_three_prime, Some(ThreePrime::Hydroxyl));
    assert_eq!(e.created_five_prime, Some(FivePrime::Phosphate));
    audit_source_partition(&s, &report.product, &l).unwrap();
    replay(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
        &report,
    )
    .unwrap();
}
#[test]
fn explicit_cut_and_two_seals_close_main_48_calls_and_keep_free_4_component() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let r = full();
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert!(report.all_requested_steps_applied);
    assert_eq!(base_count(&report.product), 52);
    assert_eq!(report.product.chains.len(), 3);
    assert!(report.product.nicks.is_empty());
    let closed: Vec<_> = report
        .product
        .chains
        .iter()
        .filter(|c| matches!(c.ends, ChainEnds::CovalentlyCircular))
        .collect();
    assert_eq!(closed.len(), 2);
    assert!(closed.iter().all(|c| c.sequence.len() == 24));
    assert!(
        closed
            .iter()
            .any(|c| c.sequence.as_str() == "AAGGGGCATCCCTTCACCATGACG")
    );
    let component = report
        .product
        .components
        .iter()
        .find(|c| c.covalently_closed)
        .unwrap();
    assert_eq!(component.paired_nucleotides, 48);
    assert_eq!(component.unpaired_nucleotides, 0);
    let free = report
        .product
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert!(matches!(
        free.ends,
        ChainEnds::Linear {
            five_prime: FivePrime::Phosphate,
            three_prime: ThreePrime::Hydroxyl
        }
    ));
    assert!(
        report
            .product
            .components
            .iter()
            .any(|c| !c.covalently_closed && c.unpaired_nucleotides == 4)
    );
    let StepEvidence::Process(p::ActionEvidence::NickSealing(last)) = &report.steps[4] else {
        panic!()
    };
    assert_eq!(last.selection.three, call(3, 5));
    assert_eq!(
        last.three_prime,
        ThreePrime::Hydroxyl,
        "internal newly cut endpoint must use explicit hydrolysis chemistry"
    );
    assert_eq!(last.outcome, Outcome::Applied);
}
#[test]
fn declared_order_duplicate_cut_and_wrong_actual_bond_do_not_imply_repair() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let early = request(vec![
        seal(call(3, 5), call(1, 0)),
        cut(),
        seal(call(2, 9), call(0, 0)),
    ]);
    let result = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&early).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_eq!(
        result.steps[2].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(result.steps[3].outcome(), Outcome::Applied);
    assert_eq!(result.product.nicks.len(), 1);
    assert!(!result.all_requested_steps_applied);
    assert_ne!(
        batch_revision(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&early).iter().cloned())
                .collect::<Vec<Step>>(),
            &l
        )
        .unwrap(),
        batch_revision(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&full()).iter().cloned())
                .collect::<Vec<Step>>(),
            &l
        )
        .unwrap()
    );
    let duplicate = request(vec![cut(), cut()]);
    let result = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&duplicate).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_eq!(
        result.steps[3].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(
        result.product,
        simulate(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&request(vec![cut()])).iter().cloned())
                .collect::<Vec<Step>>(),
            &l
        )
        .unwrap()
        .product
    );
    let mut wrong = cut();
    let p::Action::CleaveThreePrimeFlap { bond, .. } = &mut wrong else {
        panic!()
    };
    bond.upstream = call(3, 4);
    bond.downstream = call(3, 5);
    let result = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&request(vec![wrong])).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_eq!(result.steps[2].outcome(), Outcome::Unsupported);
    assert_eq!(result.product, simulate(&s, &t, &l).unwrap().product);
    let mut wrong = cut();
    let p::Action::CleaveThreePrimeFlap {
        displacement_step, ..
    } = &mut wrong
    else {
        panic!()
    };
    *displacement_step = 1;
    let result = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&request(vec![wrong])).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_eq!(result.steps[2].outcome(), Outcome::Unsupported);
    assert_eq!(result.product, simulate(&s, &t, &l).unwrap().product);
}
#[test]
fn distal_blocked_chemistry_is_retained_and_missing_phosphate_is_never_created_for_sealing() {
    let mut s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    s.strands[3].three_prime = ThreePrime::Blocked;
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&full()).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    let free = report
        .product
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert_eq!(
        free.ends,
        ChainEnds::Linear {
            five_prime: FivePrime::Phosphate,
            three_prime: ThreePrime::Blocked
        }
    );
    assert!(
        report.all_requested_steps_applied,
        "created upstream3OH is distinct from preserved original distal blocked terminus"
    );
    let mut s = directional();
    s.strands[1].five_prime = FivePrime::Hydroxyl;
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&full()).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_eq!(
        report.steps[4].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(report.product.nicks.len(), 1);
    let StepEvidence::Process(p::ActionEvidence::NickSealing(e)) = &report.steps[4] else {
        panic!()
    };
    assert_eq!(e.three_prime, ThreePrime::Hydroxyl);
    assert_eq!(e.five_prime, FivePrime::Hydroxyl);
}
#[test]
fn every_source_call_and_created_endpoint_is_checked_without_hiding_excised_material() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let r = full();
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    let free = report
        .product
        .chains
        .iter()
        .position(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    let mut missing = report.product.clone();
    missing.chains.remove(free);
    assert!(audit_source_partition(&s, &missing, &l).is_err());
    let mut duplicate = report.product.clone();
    duplicate.chains.push(duplicate.chains[free].clone());
    assert!(audit_source_partition(&s, &duplicate, &l).is_err());
    let mut gap = report.product.clone();
    gap.chains[free].sources[0].start = 7;
    assert!(audit_source_partition(&s, &gap, &l).is_err());
    let mut changed = report.product.clone();
    changed.chains[free].sequence = Sequence::dna("GGTA").unwrap();
    assert!(audit_source_partition(&s, &changed, &l).is_err());
    let mut changed = report.clone();
    let StepEvidence::Process(p::ActionEvidence::FlapCleavage(e)) = &mut changed.steps[2] else {
        panic!()
    };
    e.created_five_prime = Some(FivePrime::Hydroxyl);
    assert!(changed.validate_binding(&report.batch_revision).is_err());
    let mut changed = report.clone();
    changed.product.chains[free].ends = ChainEnds::Linear {
        five_prime: FivePrime::Hydroxyl,
        three_prime: ThreePrime::Hydroxyl,
    };
    assert!(
        validate_report(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &l,
            &changed
        )
        .is_err()
    );
    let mut changed = r.clone();
    let Step::Process(p::Action::CleaveThreePrimeFlap { assumption, .. }) = &mut changed[0] else {
        panic!()
    };
    assumption.push_str(" altered");
    assert!(
        replay(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&changed).iter().cloned())
                .collect::<Vec<Step>>(),
            &l,
            &report
        )
        .is_err()
    );
}
#[test]
fn unanchored_displacements_and_unknown_activation_do_not_gain_an_anchored_cut() {
    let mut s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    // Removing the old complementary6-call duplex continuation is a different
    // declared physical geometry, not evidence for deleting the incoming tail.
    s.pairing[1].length = 4;
    s.pairing[1].start_a = 0;
    s.pairing[1].start_b = 6;
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&request(vec![cut()])).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_ne!(report.steps[2].outcome(), Outcome::Applied);
    assert_eq!(report.product, simulate(&s, &t, &l).unwrap().product);
    let mut s = directional();
    s.strands[0].three_prime = ThreePrime::Unknown;
    let report = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&request(vec![cut()])).iter().cloned())
            .collect::<Vec<Step>>(),
        &l,
    )
    .unwrap();
    assert_ne!(report.steps[2].outcome(), Outcome::Applied);
    assert_eq!(report.product, simulate(&s, &t, &l).unwrap().product);
}
#[test]
fn one_transcript_bounds_steps_source_output_and_work_before_execution() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let transcript = t.into_iter().chain(full()).collect::<Vec<_>>();
    let mut bound = l.clone();
    bound.max_steps = 4;
    assert!(preflight(&s, &transcript, &bound).is_err());
    bound.max_steps = 5;
    assert!(preflight(&s, &transcript, &bound).is_ok());
    bound.max_output_bases = 51;
    assert!(preflight(&s, &transcript, &bound).is_err());
    bound = l.clone();
    bound.max_work_bytes = 1;
    assert!(simulate(&s, &transcript, &bound).is_err());
    let mut bad = transcript.clone();
    let Step::Process(p::Action::CleaveThreePrimeFlap { bond, .. }) = &mut bad[2] else {
        panic!()
    };
    bond.upstream.strand = 99;
    assert!(preflight(&s, &bad, &l).is_err());
    let report = simulate(&s, &transcript, &l).unwrap();
    let mut other = s.clone();
    other.strands[0].id.push_str(" changed");
    assert!(replay(&other, &transcript, &l, &report).is_err());
    let ta = ta();
    let steps = transfers(Geometry::ThreePrimeTA);
    let ta_report = simulate(&ta, &steps, &l).unwrap();
    assert_eq!(ta_report.product.nicks.len(), 2);
    assert!(
        ta_report
            .steps
            .iter()
            .all(|step| matches!(step, StepEvidence::Transfer(_)))
    );
}
