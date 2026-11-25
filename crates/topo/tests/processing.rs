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
fn request(actions: Vec<p::Action>) -> p::Request {
    p::Request {
        actions,
        limits: p::ProcessingLimits::default(),
    }
}
fn full() -> p::Request {
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
    let report = p::simulate(&s, &t, &l, &r).unwrap();
    assert_eq!(report.parent.product.chains.len(), 2);
    assert_eq!(report.product.chains.len(), 3);
    assert_eq!(report.parent.product.nicks.len(), 1);
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
        report.parent.product.source_pairing
    );
    assert_eq!(
        report.product.released_activations,
        vec!["enzyme-vector-top", "enzyme-vector-bottom"]
    );
    let p::ActionEvidence::FlapCleavage(e) = &report.steps[0] else {
        panic!()
    };
    assert_eq!(e.outcome, Outcome::Applied);
    assert_eq!(
        e.flap_calls,
        vec![call(3, 6), call(3, 7), call(3, 8), call(3, 9)]
    );
    assert_eq!(e.created_three_prime, Some(ThreePrime::Hydroxyl));
    assert_eq!(e.created_five_prime, Some(FivePrime::Phosphate));
    p::audit_source_partition(&s, &report.product, &l, &r.limits).unwrap();
    p::replay(&s, &t, &l, &r, &report).unwrap();
}
#[test]
fn explicit_cut_and_two_seals_close_main_48_calls_and_keep_free_4_component() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let r = full();
    let report = p::simulate(&s, &t, &l, &r).unwrap();
    assert!(report.all_requested_actions_applied);
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
    let p::ActionEvidence::NickSealing(last) = &report.steps[2] else {
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
    let result = p::simulate(&s, &t, &l, &early).unwrap();
    assert_eq!(
        result.steps[0].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(result.steps[1].outcome(), Outcome::Applied);
    assert_eq!(result.product.nicks.len(), 1);
    assert!(!result.all_requested_actions_applied);
    assert_ne!(
        p::revision(&s, &t, &l, &early).unwrap(),
        p::revision(&s, &t, &l, &full()).unwrap()
    );
    let duplicate = request(vec![cut(), cut()]);
    let result = p::simulate(&s, &t, &l, &duplicate).unwrap();
    assert_eq!(
        result.steps[1].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(
        result.product,
        p::simulate(&s, &t, &l, &request(vec![cut()]))
            .unwrap()
            .product
    );
    let mut wrong = cut();
    let p::Action::CleaveThreePrimeFlap { bond, .. } = &mut wrong else {
        panic!()
    };
    bond.upstream = call(3, 4);
    bond.downstream = call(3, 5);
    let result = p::simulate(&s, &t, &l, &request(vec![wrong])).unwrap();
    assert_eq!(result.steps[0].outcome(), Outcome::Unsupported);
    assert_eq!(result.product, result.parent.product);
    let mut wrong = cut();
    let p::Action::CleaveThreePrimeFlap {
        displacement_step, ..
    } = &mut wrong
    else {
        panic!()
    };
    *displacement_step = 1;
    let result = p::simulate(&s, &t, &l, &request(vec![wrong])).unwrap();
    assert_eq!(result.steps[0].outcome(), Outcome::Unsupported);
    assert_eq!(result.product, result.parent.product);
}
#[test]
fn distal_blocked_chemistry_is_retained_and_missing_phosphate_is_never_created_for_sealing() {
    let mut s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    s.strands[3].three_prime = ThreePrime::Blocked;
    let report = p::simulate(&s, &t, &l, &full()).unwrap();
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
        report.all_requested_actions_applied,
        "created upstream3OH is distinct from preserved original distal blocked terminus"
    );
    let mut s = directional();
    s.strands[1].five_prime = FivePrime::Hydroxyl;
    let report = p::simulate(&s, &t, &l, &full()).unwrap();
    assert_eq!(
        report.steps[2].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(report.product.nicks.len(), 1);
    let p::ActionEvidence::NickSealing(e) = &report.steps[2] else {
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
    let report = p::simulate(&s, &t, &l, &r).unwrap();
    let free = report
        .product
        .chains
        .iter()
        .position(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    let mut missing = report.product.clone();
    missing.chains.remove(free);
    assert!(p::audit_source_partition(&s, &missing, &l, &r.limits).is_err());
    let mut duplicate = report.product.clone();
    duplicate.chains.push(duplicate.chains[free].clone());
    assert!(p::audit_source_partition(&s, &duplicate, &l, &r.limits).is_err());
    let mut gap = report.product.clone();
    gap.chains[free].sources[0].start = 7;
    assert!(p::audit_source_partition(&s, &gap, &l, &r.limits).is_err());
    let mut changed = report.product.clone();
    changed.chains[free].sequence = Sequence::dna("GGTA").unwrap();
    assert!(p::audit_source_partition(&s, &changed, &l, &r.limits).is_err());
    let mut changed = report.clone();
    let p::ActionEvidence::FlapCleavage(e) = &mut changed.steps[0] else {
        panic!()
    };
    e.created_five_prime = Some(FivePrime::Hydroxyl);
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut changed = report.clone();
    changed.product.chains[free].ends = ChainEnds::Linear {
        five_prime: FivePrime::Hydroxyl,
        three_prime: ThreePrime::Hydroxyl,
    };
    assert!(p::validate_report(&s, &t, &l, &r, &changed).is_err());
    let mut changed = r.clone();
    let p::Action::CleaveThreePrimeFlap { assumption, .. } = &mut changed.actions[0] else {
        panic!()
    };
    assumption.push_str(" altered");
    assert!(p::replay(&s, &t, &l, &changed, &report).is_err());
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
    let report = p::simulate(&s, &t, &l, &request(vec![cut()])).unwrap();
    assert_ne!(report.steps[0].outcome(), Outcome::Applied);
    assert_eq!(report.product, report.parent.product);
    let mut s = directional();
    s.strands[0].three_prime = ThreePrime::Unknown;
    let report = p::simulate(&s, &t, &l, &request(vec![cut()])).unwrap();
    assert_ne!(report.steps[0].outcome(), Outcome::Applied);
    assert_eq!(report.product, report.parent.product);
}
#[test]
fn source_and_processing_caps_are_explicit_and_empty_actions_preserve_legacy_graphs() {
    let s = directional();
    let t = transfers(Geometry::DirectionalCacc);
    let l = Limits::default();
    let r = request(vec![]);
    let legacy = simulate(&s, &t, &l).unwrap();
    let result = p::simulate(&s, &t, &l, &r).unwrap();
    assert_eq!(
        serde_json::to_vec(&legacy).unwrap(),
        serde_json::to_vec(&result.parent).unwrap()
    );
    assert_eq!(result.product, legacy.product);
    let mut limit = l.clone();
    limit.max_steps = 4;
    assert!(p::preflight(&s, &t, &limit, &full()).is_err());
    let mut bound = full();
    bound.limits.max_product_chains = 4;
    assert!(p::preflight(&s, &t, &l, &bound).is_err());
    let mut bound = full();
    bound.limits.max_actions = 2;
    assert!(p::request_preflight(&bound, &l).is_err());
    let mut bound = full();
    bound.limits.max_mapping_segments = 1;
    assert!(p::simulate(&s, &t, &l, &bound).is_err());
    let mut changed = s.clone();
    changed.strands[0].id.push_str(" changed identity");
    assert!(p::replay(&changed, &t, &l, &r, &result).is_err());
    // The old TA path remains available unchanged, without any inferred cuts.
    let ta = ta();
    let tt = transfers(Geometry::ThreePrimeTA);
    let old = simulate(&ta, &tt, &l).unwrap();
    let new = p::simulate(&ta, &tt, &l, &r).unwrap();
    assert_eq!(new.parent, old);
    assert_eq!(new.product, old.product);
}
