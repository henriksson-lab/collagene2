use collagene_model::Sequence;
use collagene_topo::{nick_sealing as n, *};
fn strand(id: &str, dna: &str, charged: bool) -> Strand {
    Strand {
        role: "source".into(),
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
                activation_id: id.into(),
            }
        } else {
            ThreePrime::Hydroxyl
        },
    }
}
fn pair(a: usize, x: usize, b: usize, y: usize, length: usize) -> PairSpan {
    PairSpan {
        strand_a: a,
        start_a: x,
        strand_b: b,
        start_b: y,
        length,
    }
}
fn call(strand: usize, index: usize) -> Call {
    Call { strand, index }
}
fn nick(a: usize, x: usize, b: usize, y: usize) -> n::NickSelection {
    n::NickSelection {
        three: call(a, x),
        five: call(b, y),
    }
}
fn request(nicks: Vec<n::NickSelection>) -> n::Request {
    n::Request{ligase_assumption:"Active T4 DNA ligase with suitable ATP-containing buffer; ideal exact-template sealing only".into(),nicks}
}
fn simple() -> Substrates {
    let mut right = strand("right", "GT", false);
    right.five_prime = FivePrime::Phosphate;
    Substrates {
        strands: vec![
            strand("left", "AC", false),
            right,
            strand("template", "ACGT", false),
        ],
        pairing: vec![pair(0, 0, 2, 2, 2), pair(1, 0, 2, 0, 2)],
    }
}
fn ta() -> (Substrates, Vec<Transfer>) {
    (
        Substrates {
            strands: vec![
                strand("vector-top", "AGGGGCATCCCTT", true),
                strand("vector-bottom", "AGGGATGCCCCTT", true),
                strand("insert-top", "CGAATGA", false),
                strand("insert-bottom", "CATTCGA", false),
            ],
            pairing: vec![pair(0, 0, 1, 0, 12), pair(2, 0, 3, 0, 6)],
        },
        vec![
            Transfer {
                donor: 0,
                acceptor: 2,
                geometry: Geometry::ThreePrimeTA,
            },
            Transfer {
                donor: 1,
                acceptor: 3,
                geometry: Geometry::ThreePrimeTA,
            },
        ],
    )
}
#[test]
fn literal_ta_one_then_two_nicks_close_only_backbones() {
    let (s, t) = ta();
    let l = Limits::default();
    let first = nick(2, 6, 0, 0);
    let second = nick(3, 6, 1, 0);
    let one = n::simulate(&s, &t, &l, &request(vec![first.clone()])).unwrap();
    assert_eq!(one.parent.product.nicks.len(), 2);
    assert_eq!(one.product.nicks.len(), 1);
    assert_eq!(one.steps[0].outcome, Outcome::Applied);
    assert_eq!(
        one.product
            .chains
            .iter()
            .filter(|c| matches!(c.ends, ChainEnds::CovalentlyCircular))
            .count(),
        1
    );
    let r = request(vec![first, second]);
    let both = n::simulate(&s, &t, &l, &r).unwrap();
    assert!(both.all_requested_nicks_sealed);
    assert!(both.product.nicks.is_empty());
    assert_eq!(both.source_bases, 40);
    assert_eq!(both.product.chains.len(), 2);
    assert!(
        both.product
            .chains
            .iter()
            .all(|c| matches!(c.ends, ChainEnds::CovalentlyCircular))
    );
    assert_eq!(both.product.components.len(), 1);
    assert!(both.product.components[0].covalently_closed);
    assert_eq!(
        both.product.source_pairing,
        both.parent.product.source_pairing
    );
    assert_eq!(
        both.product.released_activations,
        vec!["vector-top", "vector-bottom"]
    );
    let mut seq: Vec<_> = both
        .product
        .chains
        .iter()
        .map(|c| c.sequence.as_str())
        .collect();
    seq.sort();
    assert_eq!(seq, vec!["AGGGATGCCCCTTCATTCGA", "AGGGGCATCCCTTCGAATGA"]);
    n::replay(&s, &t, &l, &r, &both).unwrap();
}
#[test]
fn exact_linear_nick_and_consumed_endpoint_have_distinct_evidence() {
    let s = simple();
    let r = request(vec![nick(0, 1, 1, 0), nick(0, 1, 1, 0)]);
    let p = n::simulate(&s, &[], &Limits::default(), &r).unwrap();
    assert_eq!(p.steps[0].outcome, Outcome::Applied);
    assert_eq!(p.steps[1].outcome, Outcome::BlockedUnderDeclaredModel);
    assert!(!p.all_requested_nicks_sealed);
    assert_eq!(p.product.chains.len(), 2);
    assert!(p.product.nicks.is_empty());
    assert!(
        p.product
            .chains
            .iter()
            .all(|c| matches!(c.ends, ChainEnds::Linear { .. }))
    );
    assert_eq!(
        p.product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        8
    );
}
#[test]
fn five_hydroxyl_unknown_and_three_activation_are_never_repaired() {
    for (five, outcome) in [
        (FivePrime::Hydroxyl, Outcome::BlockedUnderDeclaredModel),
        (FivePrime::Unknown, Outcome::Unknown),
        (FivePrime::Blocked, Outcome::BlockedUnderDeclaredModel),
    ] {
        let mut s = simple();
        s.strands[1].five_prime = five;
        let p = n::simulate(
            &s,
            &[],
            &Limits::default(),
            &request(vec![nick(0, 1, 1, 0)]),
        )
        .unwrap();
        assert_eq!(p.steps[0].outcome, outcome);
        assert_eq!(p.product, p.parent.product);
    }
    let mut s = simple();
    s.strands[0].three_prime = ThreePrime::Phosphotyrosyl {
        enzyme: Enzyme::VacciniaTypeIBTyr274,
        activation_id: "unused-activation".into(),
    };
    let p = n::simulate(
        &s,
        &[],
        &Limits::default(),
        &request(vec![nick(0, 1, 1, 0)]),
    )
    .unwrap();
    assert_eq!(p.steps[0].outcome, Outcome::Unsupported);
    assert!(p.product.released_activations.is_empty());
    assert_eq!(p.product, p.parent.product);
}
#[test]
fn unpaired_discontinuous_wrong_and_ambiguous_templates_remain_unsealed() {
    let mut unpaired = simple();
    unpaired.pairing.clear();
    let mut split = simple();
    split.strands[2].sequence = Sequence::dna("AC").unwrap();
    split.strands.push(strand("template-left", "GT", false));
    split.pairing = vec![pair(0, 0, 3, 0, 2), pair(1, 0, 2, 0, 2)];
    for s in [unpaired, split] {
        let p = n::simulate(
            &s,
            &[],
            &Limits::default(),
            &request(vec![nick(0, 1, 1, 0)]),
        )
        .unwrap();
        assert_eq!(p.steps[0].outcome, Outcome::Unsupported);
        assert_eq!(p.product, p.parent.product);
    }
    let p = n::simulate(
        &simple(),
        &[],
        &Limits::default(),
        &request(vec![nick(1, 1, 0, 0)]),
    )
    .unwrap();
    assert_ne!(p.steps[0].outcome, Outcome::Applied);
    let mut amb = simple();
    amb.strands[0].sequence = Sequence::dna("AN").unwrap();
    let p = n::simulate(
        &amb,
        &[],
        &Limits::default(),
        &request(vec![nick(0, 1, 1, 0)]),
    )
    .unwrap();
    assert_eq!(p.steps[0].outcome, Outcome::Unknown);
    assert_eq!(p.product, p.parent.product);
}
#[test]
fn directional_cacc_flap_and_source_calls_are_retained() {
    let s = Substrates {
        strands: vec![
            strand("vtop", "AAGGGGCATCCCTT", true),
            strand("vbottom", "GGTGAAGGGATGCCCCTT", true),
            strand("itop", "CACCATGACG", false),
            strand("ibottom", "CGTCATGGTG", false),
        ],
        pairing: vec![pair(0, 0, 1, 4, 14), pair(2, 0, 3, 0, 10)],
    };
    let t = vec![
        Transfer {
            donor: 0,
            acceptor: 2,
            geometry: Geometry::DirectionalCacc,
        },
        Transfer {
            donor: 1,
            acceptor: 3,
            geometry: Geometry::Flush,
        },
    ];
    let parent = simulate(&s, &t, &Limits::default()).unwrap();
    assert!(parent.steps[0].displacement.is_some());
    let r = request(
        parent
            .product
            .nicks
            .iter()
            .map(|x| n::NickSelection {
                three: x.three,
                five: x.five,
            })
            .collect(),
    );
    let p = n::simulate(&s, &t, &Limits::default(), &r).unwrap();
    assert_eq!(p.source_bases, 52);
    assert_eq!(p.parent.steps, parent.steps);
    assert_eq!(p.product.source_pairing, parent.product.source_pairing);
    assert_eq!(
        p.product.released_activations,
        parent.product.released_activations
    );
    assert_eq!(
        p.parent.steps[0]
            .displacement
            .as_ref()
            .unwrap()
            .calls_five_to_three,
        vec![call(3, 6), call(3, 7), call(3, 8), call(3, 9)]
    );
    assert_eq!(
        p.product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
    assert!(
        p.product
            .unpaired
            .iter()
            .any(|r| r.displaced_by.contains(&0))
    );
}
#[test]
fn report_tampering_stale_tokens_and_bad_source_ownership_fail() {
    let s = simple();
    let r = request(vec![nick(0, 1, 1, 0)]);
    let l = Limits::default();
    let p = n::simulate(&s, &[], &l, &r).unwrap();
    let mut restored = p.clone();
    restored.product.nicks = p.parent.product.nicks.clone();
    assert!(restored.validate_binding(&p.context_revision).is_err());
    let mut altered = p.clone();
    altered.product.chains[0].sequence = Sequence::dna("AAAA").unwrap();
    assert!(n::validate_report(&s, &[], &l, &r, &altered).is_err());
    let mut other = r.clone();
    other.ligase_assumption.push_str(" changed");
    assert!(n::validate_report(&s, &[], &l, &other, &p).is_err());
    let mut foreign = s.clone();
    foreign.strands[0].id = "different-source".into();
    assert!(n::validate_report(&foreign, &[], &l, &r, &p).is_err());
    let json = serde_json::to_vec(&p).unwrap();
    let restored: n::Report = serde_json::from_slice(&json).unwrap();
    assert_eq!(restored, p);
    n::replay(&s, &[], &l, &r, &restored).unwrap();
}
#[test]
fn combined_steps_allocations_and_calls_are_bounded_before_execution() {
    let (s, t) = ta();
    let r = request(vec![nick(2, 6, 0, 0), nick(3, 6, 1, 0)]);
    let mut l = Limits::default();
    l.max_steps = 3;
    assert!(n::preflight(&s, &t, &l, &r).is_err());
    l.max_steps = 4;
    assert!(n::preflight(&s, &t, &l, &r).is_ok());
    l.max_output_bases = 79;
    assert!(n::preflight(&s, &t, &l, &r).is_err());
    l.max_output_bases = 80;
    l.max_work_bytes = 1;
    assert!(n::preflight(&s, &t, &l, &r).is_err());
    assert!(
        n::preflight(
            &s,
            &t,
            &Limits::default(),
            &request(vec![nick(99, 0, 0, 0)])
        )
        .is_err()
    );
    assert!(
        n::request_preflight(
            &n::Request {
                ligase_assumption: " ".into(),
                nicks: vec![]
            },
            &Limits::default()
        )
        .is_err()
    );
    let p = n::simulate(&simple(), &[], &Limits::default(), &request(vec![])).unwrap();
    assert!(!p.all_requested_nicks_sealed);
    assert_eq!(p.product, p.parent.product);
}
