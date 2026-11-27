use collagene_model::Sequence;
use collagene_topo::{nick_sealing as n, projection as p, projection as contour, *};
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
fn ta() -> (Substrates, Vec<Step>) {
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
            Step::Transfer(Transfer {
                donor: Call {
                    strand: 0,
                    index: 12,
                },
                acceptor: Call {
                    strand: 2,
                    index: 0,
                },
                geometry: Geometry::ThreePrimeTA,
            }),
            Step::Transfer(Transfer {
                donor: Call {
                    strand: 1,
                    index: 12,
                },
                acceptor: Call {
                    strand: 3,
                    index: 0,
                },
                geometry: Geometry::ThreePrimeTA,
            }),
        ],
    )
}
fn request(count: usize) -> Vec<Step> {
    [(call(2, 6), call(0, 0)), (call(3, 6), call(1, 0))]
        .into_iter()
        .take(count)
        .map(|(three, five)| {
            Step::Process(processing::Action::SealNick {
                selection: n::NickSelection { three, five },
                ligase_assumption: "Explicit active T4 DNA ligase in suitable ATP buffer".into(),
            })
        })
        .collect()
}
fn selection(s: &Substrates, steps: &[Step], seed: Call) -> p::Selection {
    p::Selection {
        source_batch_revision: batch_revision(s, steps, &Limits::default()).unwrap(),
        component: 0,
        reference_seed: seed,
    }
}
fn view(r: &p::ProjectionReport) -> &contour::DuplexView {
    match &r.outcome {
        contour::Outcome::Proven(v) => v,
        _ => panic!("expected paired contour: {:?}", r.outcome),
    }
}
#[test]
fn fresh_ta_processed_contour_retains_one_nick_then_proves_closed_backbones() {
    let (s, t) = ta();
    let l = Limits::default();
    let o = contour::Options::default();
    for (count, nicks) in [(0, 2), (1, 1), (2, 0)] {
        let r = request(count);
        let select = selection(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            call(0, 0),
        );
        let result = p::project(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &l,
            &select,
            &o,
        )
        .unwrap();
        let v = view(&result);
        assert_eq!(v.topology, contour::ContourTopology::Circular);
        assert_eq!(v.top.as_str(), "AGGGGCATCCCTTCGAATGA");
        assert_eq!(v.bottom.as_str(), "TCATTCGAAGGGATGCCCCT");
        assert_eq!(v.nicks.len(), nicks);
        assert_eq!(result.physical.product.nicks.len(), nicks);
        assert_eq!(result.physical.source_bases, 40);
        assert_eq!(
            simulate(&s, &t, &Limits::default())
                .unwrap()
                .product
                .nicks
                .len(),
            2
        );
        p::validate_report(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &l,
            &select,
            &o,
            &result,
        )
        .unwrap();
        p::replay(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &l,
            &select,
            &o,
            &result,
        )
        .unwrap();
    }
}
#[test]
fn origin_and_reverse_seed_use_literal_source_calls_without_target_guidance() {
    let (s, t) = ta();
    let r = request(2);
    let o = contour::Options::default();
    for (seed, letters) in [
        (call(2, 3), "ATGAAGGGGCATCCCTTCGA"),
        (call(1, 0), "AGGGATGCCCCTTCATTCGA"),
    ] {
        let select = selection(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            seed,
        );
        let report = p::project(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &Limits::default(),
            &select,
            &o,
        )
        .unwrap();
        let v = view(&report);
        assert_eq!(v.origin, seed);
        assert_eq!(v.seed_reference_position, 0);
        assert_eq!(v.top.as_str(), letters);
        assert_eq!(v.source_maps[0].top_strand, seed.strand);
        assert_eq!(v.source_maps[0].top_start, seed.index);
        p::validate_report(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &Limits::default(),
            &select,
            &o,
            &report,
        )
        .unwrap();
    }
}
#[test]
fn other_components_and_unreacted_source_calls_are_never_deleted() {
    let (mut s, t) = ta();
    s.strands.push(strand("free-byproduct", "GGAC", false));
    let r = request(2);
    let fresh = simulate(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
    )
    .unwrap();
    let component = fresh
        .product
        .components
        .iter()
        .position(|component| {
            component.chains.iter().any(|&chain| {
                fresh.product.chains[chain]
                    .sources
                    .iter()
                    .any(|source| source.source_strand == 0 && source.start == 0 && source.end > 0)
            })
        })
        .unwrap();
    assert_eq!(fresh.product.components[component].paired_nucleotides, 40);
    assert_eq!(fresh.product.components[component].unpaired_nucleotides, 0);
    let select = p::Selection {
        source_batch_revision: fresh.batch_revision,
        component,
        reference_seed: call(0, 0),
    };
    let result = p::project(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
        &select,
        &contour::Options::default(),
    )
    .unwrap();
    assert_eq!(view(&result).top.len(), 20);
    assert_eq!(result.physical.source_bases, 44);
    assert_eq!(result.physical.product.components.len(), 2);
    assert!(
        result
            .physical
            .product
            .components
            .iter()
            .any(|c| { c.paired_nucleotides == 0 && c.unpaired_nucleotides == 4 })
    );
    assert!(
        result
            .physical
            .product
            .chains
            .iter()
            .any(|c| c.sequence.as_str() == "GGAC")
    );
    assert_eq!(
        result
            .physical
            .product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        44
    );
}
#[test]
fn strict_requested_policy_includes_blocked_sealing_as_well_as_original_transfers() {
    let (mut s, t) = ta();
    s.strands[0].five_prime = FivePrime::Hydroxyl;
    let r = request(2);
    let select = selection(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        call(0, 0),
    );
    let mut o = contour::Options::default();
    let relaxed = p::project(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
        &select,
        &o,
    )
    .unwrap();
    assert_eq!(
        relaxed.physical.steps[2].outcome(),
        Outcome::BlockedUnderDeclaredModel
    );
    assert!(!view(&relaxed).all_requested_steps_applied);
    assert_eq!(view(&relaxed).nicks.len(), 1);
    o.require_all_requested_steps_applied = true;
    let strict = p::project(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
        &select,
        &o,
    )
    .unwrap();
    assert_eq!(
        strict.outcome,
        contour::Outcome::Withheld(contour::Withheld::RequestedStepsNotApplied)
    );
    assert_eq!(strict.physical.product, relaxed.physical.product);
}
#[test]
fn directional_flap_stays_withheld_after_explicit_sealing() {
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
        Step::Transfer(Transfer {
            donor: Call {
                strand: 0,
                index: s.strands[0].sequence.len() - 1,
            },
            acceptor: Call {
                strand: 2,
                index: 0,
            },
            geometry: Geometry::DirectionalCacc,
        }),
        Step::Transfer(Transfer {
            donor: Call {
                strand: 1,
                index: s.strands[1].sequence.len() - 1,
            },
            acceptor: Call {
                strand: 3,
                index: 0,
            },
            geometry: Geometry::Flush,
        }),
    ];
    let r = request(0);
    let select = selection(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        call(0, 0),
    );
    let result = p::project(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
        &select,
        &contour::Options::default(),
    )
    .unwrap();
    assert_eq!(
        result.outcome,
        contour::Outcome::Withheld(contour::Withheld::UnpairedNucleotides { count: 4 })
    );
    assert_eq!(
        result
            .physical
            .product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
    assert_eq!(
        (match &result.physical.steps[0] {
            StepEvidence::Transfer(e) => e,
            _ => panic!("transfer evidence"),
        })
        .displacement
        .as_ref()
        .unwrap()
        .calls_five_to_three,
        vec![call(3, 6), call(3, 7), call(3, 8), call(3, 9)]
    );
}
#[test]
fn imported_positive_reports_versions_and_old_selection_contexts_do_not_replay() {
    let (s, t) = ta();
    let r = request(2);
    let select = selection(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        call(0, 0),
    );
    let o = contour::Options::default();
    let result = p::project(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
        &select,
        &o,
    )
    .unwrap();
    let mut altered = result.clone();
    altered.engine.push_str("-future");
    assert!(altered.validate_binding(&result.batch_revision).is_err());
    let mut changed = result.clone();
    if let contour::Outcome::Proven(v) = &mut changed.outcome {
        v.top = Sequence::dna("AAAAAAAAAAAAAAAAAAAA").unwrap();
    }
    assert!(
        p::validate_report(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &Limits::default(),
            &select,
            &o,
            &changed
        )
        .is_err()
    );
    assert!(
        p::replay(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &Limits::default(),
            &select,
            &o,
            &changed
        )
        .is_err()
    );
    let mut stale = select.clone();
    stale.source_batch_revision = batch_revision(&s, &t, &Limits::default()).unwrap();
    assert!(
        p::batch_revision(
            &s,
            &(&t)
                .iter()
                .cloned()
                .chain((&r).iter().cloned())
                .collect::<Vec<Step>>(),
            &Limits::default(),
            &stale,
            &o
        )
        .is_err()
    );
    let encoded = serde_json::to_vec(&result).unwrap();
    let restored: p::ProjectionReport = serde_json::from_slice(&encoded).unwrap();
    assert_eq!(restored, result);
    p::replay(
        &s,
        &(&t)
            .iter()
            .cloned()
            .chain((&r).iter().cloned())
            .collect::<Vec<Step>>(),
        &Limits::default(),
        &select,
        &o,
        &restored,
    )
    .unwrap();
}
#[test]
fn projection_budget_withholds_only_the_view_and_retains_complete_reaction() {
    let (s, t) = ta();
    let select = selection(&s, &t, call(0, 0));
    let mut options = contour::Options::default();
    let mut too_small = Limits::default();
    too_small.max_output_bases = 39;
    assert!(matches!(
        p::project(&s, &t, &too_small, &select, &options),
        Err(contour::Error::Engine(Error::Limit {
            kind: Budget::OutputBases,
            ..
        }))
    ));
    let full = p::project(&s, &t, &Limits::default(), &select, &options).unwrap();
    assert_eq!(view(&full).nicks.len(), 2);
    options.limits.max_pairs = 19;
    let bounded = p::project(&s, &t, &Limits::default(), &select, &options).unwrap();
    assert!(matches!(
        bounded.outcome,
        contour::Outcome::Withheld(contour::Withheld::Budget {
            kind: contour::LimitKind::Pairs,
            ..
        })
    ));
    assert_eq!(bounded.physical, full.physical);
}
