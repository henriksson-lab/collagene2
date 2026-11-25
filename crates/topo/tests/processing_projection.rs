use collagene_model::Sequence;
use collagene_topo::{processing as actions, processing_projection as p, projection as contour, *};
fn call(strand: usize, index: usize) -> Call {
    Call { strand, index }
}
fn fixture() -> (Substrates, Vec<Transfer>) {
    let dna = [
        "AAGGGGCATCCCTT",
        "GGTGAAGGGATGCCCCTT",
        "CACCATGACG",
        "CGTCATGGTG",
    ];
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| Strand {
            role: "source".into(),
            id: format!("original{i}"),
            sequence: Sequence::dna(dna).unwrap(),
            five_prime: if i < 2 {
                FivePrime::Phosphate
            } else {
                FivePrime::Hydroxyl
            },
            three_prime: if i < 2 {
                ThreePrime::Phosphotyrosyl {
                    enzyme: Enzyme::VacciniaTypeIBTyr274,
                    activation_id: format!("enzyme{i}"),
                }
            } else {
                ThreePrime::Hydroxyl
            },
        })
        .collect();
    (
        Substrates {
            strands,
            pairing: vec![
                PairSpan {
                    strand_a: 0,
                    start_a: 0,
                    strand_b: 1,
                    start_b: 4,
                    length: 14,
                },
                PairSpan {
                    strand_a: 2,
                    start_a: 0,
                    strand_b: 3,
                    start_b: 0,
                    length: 10,
                },
            ],
        },
        vec![
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
        ],
    )
}
fn request(seal: bool) -> actions::Request {
    let mut list = vec![actions::Action::CleaveThreePrimeFlap {
        bond: actions::Bond {
            upstream: call(3, 5),
            downstream: call(3, 6),
        },
        displacement_step: 0,
        profile: actions::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,
        assumption: "Explicit ideal hydrolytic junction cut, no named enzyme claim".into(),
    }];
    if seal {
        for (three, five) in [(call(2, 9), call(0, 0)), (call(3, 5), call(1, 0))] {
            list.push(actions::Action::SealNick {
                selection: nick_sealing::NickSelection { three, five },
                ligase_assumption: "Ideal templated3OH5P ligase hypothesis".into(),
            });
        }
    }
    actions::Request {
        actions: list,
        limits: actions::ProcessingLimits::default(),
    }
}
fn select(s: &Substrates, t: &[Transfer], r: &actions::Request, seed: Call) -> p::Selection {
    let fresh = actions::simulate(s, t, &Limits::default(), r).unwrap();
    let component = fresh
        .product
        .components
        .iter()
        .position(|c| {
            c.chains.iter().any(|&i| {
                fresh.product.chains[i].sources.iter().any(|m| {
                    m.source_strand == seed.strand && m.start <= seed.index && seed.index < m.end
                })
            })
        })
        .unwrap();
    p::Selection {
        processing_context_revision: fresh.context_revision,
        component,
        reference_seed: seed,
    }
}
fn view(report: &p::ProjectionReport) -> &contour::DuplexView {
    match &report.outcome {
        contour::Outcome::Proven(v) => v,
        _ => panic!("withheld {:?}", report.outcome),
    }
}
#[test]
fn literal_cut_and_seals_prove_main48_calls_and_keep_excised_four() {
    let (s, t) = fixture();
    let r = request(true);
    let selection = select(&s, &t, &r, call(0, 0));
    let report = p::project(
        &s,
        &t,
        &Limits::default(),
        &r,
        &selection,
        &contour::Options::default(),
    )
    .unwrap();
    let v = view(&report);
    assert_eq!(v.top.as_str(), "AAGGGGCATCCCTTCACCATGACG");
    assert_eq!(v.topology, contour::ContourTopology::Circular);
    assert!(v.nicks.is_empty());
    assert_eq!(report.processing.source_bases, 52);
    assert_eq!(report.processing.product.components.len(), 2);
    assert!(
        report
            .processing
            .product
            .components
            .iter()
            .any(|c| c.paired_nucleotides == 48
                && c.unpaired_nucleotides == 0
                && c.covalently_closed)
    );
    let free = report
        .processing
        .product
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert_eq!(free.sources[0].source_strand, 3);
    assert_eq!((free.sources[0].start, free.sources[0].end), (6, 10));
    assert_eq!(
        free.ends,
        ChainEnds::Linear {
            five_prime: FivePrime::Phosphate,
            three_prime: ThreePrime::Hydroxyl
        }
    );
    assert_eq!(
        report
            .processing
            .parent
            .product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
    p::replay(
        &s,
        &t,
        &Limits::default(),
        &r,
        &selection,
        &contour::Options::default(),
        &report,
    )
    .unwrap();
}
#[test]
fn newly_exposed_internal_source_three_oh_nick_is_preserved_before_sealing() {
    let (s, t) = fixture();
    let r = request(false);
    let selection = select(&s, &t, &r, call(0, 0));
    let report = p::project(
        &s,
        &t,
        &Limits::default(),
        &r,
        &selection,
        &contour::Options::default(),
    )
    .unwrap();
    assert_eq!(view(&report).nicks.len(), 2);
    assert!(view(&report).nicks.iter().any(|n| n.three == call(3, 5)
        && n.five == call(1, 0)
        && n.three_prime == ThreePrime::Hydroxyl
        && n.five_prime == FivePrime::Phosphate));
    assert!(
        report
            .processing
            .product
            .chains
            .iter()
            .any(|c| c.sequence.as_str() == "GGTG")
    );
}
#[test]
fn reverse_seed_proves_actual_source_orientation_and_free_component_is_withheld() {
    let (s, t) = fixture();
    let r = request(true);
    let selection = select(&s, &t, &r, call(1, 0));
    let report = p::project(
        &s,
        &t,
        &Limits::default(),
        &r,
        &selection,
        &contour::Options::default(),
    )
    .unwrap();
    assert_eq!(view(&report).origin, call(1, 0));
    assert_eq!(view(&report).top.as_str(), "GGTGAAGGGATGCCCCTTCGTCAT");
    let free = select(&s, &t, &r, call(3, 6));
    let withheld = p::project(
        &s,
        &t,
        &Limits::default(),
        &r,
        &free,
        &contour::Options::default(),
    )
    .unwrap();
    assert_eq!(
        withheld.outcome,
        contour::Outcome::Withheld(contour::Withheld::UnpairedNucleotides { count: 4 })
    );
    assert_eq!(withheld.processing.source_bases, 52);
}
#[test]
fn absent_cut_and_strict_blocked_actions_never_invent_repaired_contours() {
    let (s, t) = fixture();
    let mut no_cut = request(false);
    no_cut.actions.clear();
    let selection = select(&s, &t, &no_cut, call(0, 0));
    let report = p::project(
        &s,
        &t,
        &Limits::default(),
        &no_cut,
        &selection,
        &contour::Options::default(),
    )
    .unwrap();
    assert_eq!(
        report.outcome,
        contour::Outcome::Withheld(contour::Withheld::UnpairedNucleotides { count: 4 })
    );
    let mut repeated = request(true);
    repeated.actions.push(repeated.actions[0].clone());
    let selection = select(&s, &t, &repeated, call(0, 0));
    let mut options = contour::Options::default();
    options.require_all_requested_transfers_applied = true;
    let blocked = p::project(&s, &t, &Limits::default(), &repeated, &selection, &options).unwrap();
    assert_eq!(
        blocked.outcome,
        contour::Outcome::Withheld(contour::Withheld::RequestedTransfersNotApplied)
    );
    assert_eq!(
        blocked.processing.steps[3].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(
        blocked
            .processing
            .product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
}
#[test]
fn changed_requests_foreign_versions_maps_and_dropped_flaps_do_not_replay() {
    let (s, t) = fixture();
    let r = request(true);
    let selection = select(&s, &t, &r, call(0, 0));
    let options = contour::Options::default();
    let report = p::project(&s, &t, &Limits::default(), &r, &selection, &options).unwrap();
    let mut wrong = report.clone();
    wrong.engine.push('x');
    assert!(wrong.validate_binding(&report.context_revision).is_err());
    let mut wrong = report.clone();
    wrong
        .processing
        .product
        .chains
        .retain(|c| c.sequence.as_str() != "GGTG");
    assert!(
        p::validate_report(&s, &t, &Limits::default(), &r, &selection, &options, &wrong).is_err()
    );
    let mut changed = r.clone();
    changed.actions.swap(0, 1);
    assert!(p::revision(&s, &t, &Limits::default(), &changed, &selection, &options).is_err());
    let mut wrong = report.clone();
    if let contour::Outcome::Proven(v) = &mut wrong.outcome {
        v.top = Sequence::dna("AAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
    }
    assert!(
        p::validate_report(&s, &t, &Limits::default(), &r, &selection, &options, &wrong).is_err()
    );
}
