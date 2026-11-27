use collagene_model::Sequence;
use collagene_topo::{projection as p, *};
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
fn steps(directional: bool) -> [Step; 2] {
    [
        Step::Transfer(Transfer {
            donor: Call {
                strand: 0,
                index: if directional { 13 } else { 12 },
            },
            acceptor: Call {
                strand: 2,
                index: 0,
            },
            geometry: if directional {
                Geometry::DirectionalCacc
            } else {
                Geometry::ThreePrimeTA
            },
        }),
        Step::Transfer(Transfer {
            donor: Call {
                strand: 1,
                index: if directional { 17 } else { 12 },
            },
            acceptor: Call {
                strand: 3,
                index: 0,
            },
            geometry: if directional {
                Geometry::Flush
            } else {
                Geometry::ThreePrimeTA
            },
        }),
    ]
}
fn selection(s: &Substrates, t: &[Step], component: usize, seed: Call) -> p::Selection {
    p::Selection {
        source_batch_revision: batch_revision(s, t, &Limits::default()).unwrap(),
        component,
        reference_seed: seed,
    }
}
fn view(r: &p::ProjectionReport) -> &p::DuplexView {
    let p::Outcome::Proven(v) = &r.outcome else {
        panic!("{:?}", r.outcome)
    };
    v
}
fn audit_maps(s: &Substrates, v: &p::DuplexView) {
    let mut covered = 0;
    for m in &v.source_maps {
        assert_eq!(m.reference_start, covered);
        for k in 0..m.length {
            assert_eq!(
                v.top.as_str().as_bytes()[covered + k],
                s.strands[m.top_strand].sequence.as_str().as_bytes()[m.top_start + k]
            );
            assert_eq!(
                v.bottom.as_str().as_bytes()[v.top.len() - 1 - covered - k],
                s.strands[m.bottom_strand].sequence.as_str().as_bytes()
                    [m.bottom_start + m.length - 1 - k]
            );
        }
        covered += m.length;
    }
    assert_eq!(covered, v.top.len());
}
#[test]
fn ta_reference_cycle_retains_both_actual_nicks_and_their_chemistry() {
    let s = ta();
    let t = steps(false);
    let select = selection(
        &s,
        &t,
        0,
        Call {
            strand: 0,
            index: 0,
        },
    );
    let r = p::project(&s, &t, &Limits::default(), &select, &p::Options::default()).unwrap();
    let v = view(&r);
    assert_eq!(v.topology, p::ContourTopology::Circular);
    assert_eq!(v.top.as_str(), "AGGGGCATCCCTTCGAATGA");
    assert_eq!(v.bottom.as_str(), "TCATTCGAAGGGATGCCCCT");
    assert!(v.ends.is_none());
    assert_eq!(
        v.origin,
        Call {
            strand: 0,
            index: 0
        }
    );
    assert_eq!(v.seed_reference_position, 0);
    assert_eq!(v.nicks.len(), 2);
    assert!(v.nicks.iter().any(|n| n.strand == p::ReferenceStrand::Top
        && n.boundary == 0
        && n.three
            == Call {
                strand: 2,
                index: 6
            }
        && n.five
            == Call {
                strand: 0,
                index: 0
            }));
    assert!(
        v.nicks
            .iter()
            .any(|n| n.strand == p::ReferenceStrand::Bottom
                && n.boundary == 12
                && n.three
                    == Call {
                        strand: 3,
                        index: 6
                    }
                && n.five
                    == Call {
                        strand: 1,
                        index: 0
                    })
    );
    assert!(
        v.nicks
            .iter()
            .all(|n| n.five_prime == FivePrime::Phosphate && n.three_prime == ThreePrime::Hydroxyl)
    );
    audit_maps(&s, v);
    assert_eq!(r.physical, simulate(&s, &t, &Limits::default()).unwrap());
    p::replay(
        &s,
        &t,
        &Limits::default(),
        &select,
        &p::Options::default(),
        &r,
    )
    .unwrap();
}
#[test]
fn circular_origin_and_reference_sense_follow_source_seed_not_a_target() {
    let s = ta();
    let t = steps(false);
    let literal = "AGGGGCATCCCTTCGAATGA";
    for index in [0, 3, 7, 12] {
        let r = p::project(
            &s,
            &t,
            &Limits::default(),
            &selection(&s, &t, 0, Call { strand: 0, index }),
            &p::Options::default(),
        )
        .unwrap();
        let v = view(&r);
        assert_eq!(
            v.top.as_str(),
            format!("{}{}", &literal[index..], &literal[..index])
        );
        assert!(
            v.nicks
                .iter()
                .any(|n| n.strand == p::ReferenceStrand::Top && n.boundary == (20 - index) % 20)
        );
        assert!(v.nicks.iter().any(
            |n| n.strand == p::ReferenceStrand::Bottom && n.boundary == (12 + 20 - index) % 20
        ));
        audit_maps(&s, v);
    }
    let r = p::project(
        &s,
        &t,
        &Limits::default(),
        &selection(
            &s,
            &t,
            0,
            Call {
                strand: 1,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert_eq!(view(&r).top.as_str(), "AGGGATGCCCCTTCATTCGA");
    audit_maps(&s, view(&r));
}
#[test]
fn linear_origin_is_physical_left_endpoint_even_for_middle_seed() {
    let s = Substrates {
        strands: vec![
            strand("top", "ACGTAC", false),
            strand("bottom", "GTACGT", false),
        ],
        pairing: vec![pair(0, 0, 1, 0, 6)],
    };
    let r = p::project(
        &s,
        &[],
        &Limits::default(),
        &selection(
            &s,
            &[],
            0,
            Call {
                strand: 0,
                index: 3,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    let v = view(&r);
    assert_eq!(v.topology, p::ContourTopology::Linear);
    assert_eq!(v.top.as_str(), "ACGTAC");
    assert_eq!(v.bottom.as_str(), "GTACGT");
    assert_eq!(
        v.origin,
        Call {
            strand: 0,
            index: 0
        }
    );
    assert_eq!(v.seed_reference_position, 3);
    assert!(v.nicks.is_empty());
    let e = v.ends.as_ref().unwrap();
    assert_eq!(
        (e.left.five_prime, &e.left.three_prime),
        (FivePrime::Hydroxyl, &ThreePrime::Hydroxyl)
    );
    assert!(v.fully_known_terminal_chemistry);
    audit_maps(&s, v);
}
#[test]
fn nick_resealed_circle_does_not_erase_opposite_nick() {
    let mut top = strand("top", "GACGCCCTT", true);
    top.five_prime = FivePrime::Hydroxyl;
    let s = Substrates {
        strands: vec![top, strand("bottom", "CGTCAAGGG", false)],
        pairing: vec![pair(0, 0, 1, 0, 4), pair(0, 4, 1, 4, 5)],
    };
    let t = [Step::Transfer(Transfer {
        donor: Call {
            strand: 0,
            index: s.strands[0].sequence.len() - 1,
        },
        acceptor: Call {
            strand: 0,
            index: 0,
        },
        geometry: Geometry::NickReseal,
    })];
    let r = p::project(
        &s,
        &t,
        &Limits::default(),
        &selection(
            &s,
            &t,
            0,
            Call {
                strand: 0,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    let v = view(&r);
    assert_eq!(v.topology, p::ContourTopology::Circular);
    assert_eq!(v.top.as_str(), "GACGCCCTT");
    assert_eq!(v.bottom.as_str(), "AAGGGCGTC");
    assert_eq!(v.nicks.len(), 1);
    assert_eq!(v.nicks[0].strand, p::ReferenceStrand::Bottom);
    assert_eq!(v.nicks[0].boundary, 4);
    assert!(!r.physical.product.components[0].covalently_closed);
    audit_maps(&s, v);
}
#[test]
fn directional_flap_is_withheld_without_losing_a_single_source_call() {
    let s = directional();
    let t = steps(true);
    let r = p::project(
        &s,
        &t,
        &Limits::default(),
        &selection(
            &s,
            &t,
            0,
            Call {
                strand: 0,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert_eq!(
        r.outcome,
        p::Outcome::Withheld(p::Withheld::UnpairedNucleotides { count: 4 })
    );
    assert_eq!(
        r.physical.product.chains[1].sequence.as_str(),
        "GGTGAAGGGATGCCCCTTCGTCATGGTG"
    );
    assert_eq!(
        r.physical.product.unpaired[0].end - r.physical.product.unpaired[0].start,
        4
    );
    audit_sources(&s, &r.physical).unwrap();
}
#[test]
fn separate_free_strand_is_retained_even_when_another_component_has_a_duplex_view() {
    let s = Substrates {
        strands: vec![
            strand("donor", "GACGCCCTT", true),
            strand("opposite", "GGTGAAGGGCGTC", false),
            strand("incoming", "CACC", false),
            strand("free", "GGTG", false),
        ],
        pairing: vec![pair(0, 0, 1, 4, 9), pair(2, 0, 3, 0, 4)],
    };
    let t = [Step::Transfer(Transfer {
        donor: Call {
            strand: 0,
            index: s.strands[0].sequence.len() - 1,
        },
        acceptor: Call {
            strand: 2,
            index: 0,
        },
        geometry: Geometry::DirectionalCacc,
    })];
    let physical = simulate(&s, &t, &Limits::default()).unwrap();
    let duplex = physical
        .product
        .components
        .iter()
        .position(|c| c.paired_nucleotides > 0)
        .unwrap();
    let free = physical
        .product
        .components
        .iter()
        .position(|c| c.paired_nucleotides == 0)
        .unwrap();
    let r = p::project(
        &s,
        &t,
        &Limits::default(),
        &selection(
            &s,
            &t,
            duplex,
            Call {
                strand: 0,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert_eq!(view(&r).top.as_str(), "GACGCCCTTCACC");
    assert_eq!(r.physical.product.components.len(), 2);
    assert_eq!(r.physical.product.chains[2].sequence.as_str(), "GGTG");
    let r = p::project(
        &s,
        &t,
        &Limits::default(),
        &selection(
            &s,
            &t,
            free,
            Call {
                strand: 3,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert_eq!(
        r.outcome,
        p::Outcome::Withheld(p::Withheld::UnpairedNucleotides { count: 4 })
    );
}
#[test]
fn four_way_branch_and_paired_track_backbone_link_cannot_be_flattened() {
    let s = Substrates {
        strands: vec![
            strand("a", "AAA", false),
            strand("b", "TTT", false),
            strand("c", "T", false),
            strand("d", "T", false),
            strand("e", "A", false),
            strand("f", "A", false),
        ],
        pairing: vec![
            pair(0, 1, 1, 1, 1),
            pair(0, 0, 2, 0, 1),
            pair(0, 2, 3, 0, 1),
            pair(1, 0, 4, 0, 1),
            pair(1, 2, 5, 0, 1),
        ],
    };
    let r = p::project(
        &s,
        &[],
        &Limits::default(),
        &selection(
            &s,
            &[],
            0,
            Call {
                strand: 0,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert!(matches!(
        r.outcome,
        p::Outcome::Withheld(p::Withheld::ConflictingContinuation { .. })
    ));
    assert_eq!(r.physical.source_bases, 10);
    let s = Substrates {
        strands: vec![strand("a", "CCCTT", true), strand("b", "AAGGG", false)],
        pairing: vec![pair(0, 0, 1, 0, 5)],
    };
    let t = [Step::Transfer(Transfer {
        donor: Call {
            strand: 0,
            index: s.strands[0].sequence.len() - 1,
        },
        acceptor: Call {
            strand: 1,
            index: 0,
        },
        geometry: Geometry::Flush,
    })];
    let r = p::project(
        &s,
        &t,
        &Limits::default(),
        &selection(
            &s,
            &t,
            0,
            Call {
                strand: 0,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert!(matches!(
        r.outcome,
        p::Outcome::Withheld(p::Withheld::BackboneCrossesPairedTracks { .. })
    ));
    assert_eq!(r.physical.source_bases, 10);
    audit_sources(&s, &r.physical).unwrap();
}
#[test]
fn unknown_chemistry_remains_explicit_and_strict_policy_is_separately_bound() {
    let mut s = Substrates {
        strands: vec![
            strand("top", "ACGTAC", false),
            strand("bottom", "GTACGT", false),
        ],
        pairing: vec![pair(0, 0, 1, 0, 6)],
    };
    s.strands[0].five_prime = FivePrime::Unknown;
    let select = selection(
        &s,
        &[],
        0,
        Call {
            strand: 0,
            index: 0,
        },
    );
    let a = p::project(&s, &[], &Limits::default(), &select, &p::Options::default()).unwrap();
    assert!(!view(&a).fully_known_terminal_chemistry);
    assert_eq!(
        view(&a).ends.as_ref().unwrap().left.five_prime,
        FivePrime::Unknown
    );
    let strict = p::Options {
        require_known_terminal_chemistry: true,
        ..p::Options::default()
    };
    let b = p::project(&s, &[], &Limits::default(), &select, &strict).unwrap();
    assert_eq!(
        b.outcome,
        p::Outcome::Withheld(p::Withheld::UnknownTerminalChemistry)
    );
    assert_ne!(a.batch_revision, b.batch_revision);
    let mut ambiguous = s.clone();
    ambiguous.strands[0].sequence = Sequence::dna("ACNTAC").unwrap();
    ambiguous.strands[1].sequence = Sequence::dna("GTANGT").unwrap();
    let r = p::project(
        &ambiguous,
        &[],
        &Limits::default(),
        &selection(
            &ambiguous,
            &[],
            0,
            Call {
                strand: 0,
                index: 0,
            },
        ),
        &p::Options::default(),
    )
    .unwrap();
    assert!(matches!(
        r.outcome,
        p::Outcome::Withheld(p::Withheld::AmbiguousNucleotide { .. })
    ));
}
#[test]
fn source_seed_limits_and_imported_projection_claims_never_restore_authority() {
    let s = ta();
    let t = steps(false);
    let select = selection(
        &s,
        &t,
        0,
        Call {
            strand: 0,
            index: 0,
        },
    );
    let original = p::project(&s, &t, &Limits::default(), &select, &p::Options::default()).unwrap();
    let mut changed = s.clone();
    changed.strands[0].role = "other material".into();
    assert!(
        p::project(
            &changed,
            &t,
            &Limits::default(),
            &select,
            &p::Options::default()
        )
        .is_err()
    );
    for mode in 0..4 {
        let mut options = p::Options::default();
        match mode {
            0 => options.limits.max_pairs = 1,
            1 => options.limits.max_work_cells = 1,
            2 => options.limits.max_work_bytes = 1,
            _ => options.limits.max_mapping_segments = 1,
        };
        let r = p::project(&s, &t, &Limits::default(), &select, &options).unwrap();
        assert!(matches!(
            r.outcome,
            p::Outcome::Withheld(p::Withheld::Budget { .. })
        ));
        assert_eq!(r.physical, original.physical);
        assert_ne!(r.batch_revision, original.batch_revision);
    }
    let mut forged = original.clone();
    let p::Outcome::Proven(v) = &mut forged.outcome else {
        panic!()
    };
    v.nicks.clear();
    assert!(
        p::replay(
            &s,
            &t,
            &Limits::default(),
            &select,
            &p::Options::default(),
            &forged
        )
        .is_err()
    );
}

#[test]
fn geometric_projection_does_not_claim_an_unapplied_transfer_succeeded() {
    let s = Substrates {
        strands: vec![
            strand("top", "ACGTAC", false),
            strand("bottom", "GTACGT", false),
        ],
        pairing: vec![pair(0, 0, 1, 0, 6)],
    };
    let t = [Step::Transfer(Transfer {
        donor: Call {
            strand: 0,
            index: s.strands[0].sequence.len() - 1,
        },
        acceptor: Call {
            strand: 1,
            index: 0,
        },
        geometry: Geometry::Flush,
    })];
    let select = selection(
        &s,
        &t,
        0,
        Call {
            strand: 0,
            index: 0,
        },
    );
    let ordinary = p::project(&s, &t, &Limits::default(), &select, &p::Options::default()).unwrap();
    assert!(!view(&ordinary).all_requested_steps_applied);
    let strict = p::Options {
        require_all_requested_steps_applied: true,
        ..p::Options::default()
    };
    let withheld = p::project(&s, &t, &Limits::default(), &select, &strict).unwrap();
    assert_eq!(
        withheld.outcome,
        p::Outcome::Withheld(p::Withheld::RequestedStepsNotApplied)
    );
    assert_eq!(ordinary.physical, withheld.physical);
    assert_ne!(ordinary.batch_revision, withheld.batch_revision);
}

#[test]
fn activated_and_unknown_nick_ends_are_retained_without_repair() {
    let mut s = Substrates {
        strands: vec![
            strand("top", "GACGCCCTT", true),
            strand("bottom", "CGTCAAGGG", false),
        ],
        pairing: vec![pair(0, 0, 1, 0, 4), pair(0, 4, 1, 4, 5)],
    };
    s.strands[0].five_prime = FivePrime::Hydroxyl;
    let select = selection(
        &s,
        &[],
        0,
        Call {
            strand: 0,
            index: 0,
        },
    );
    let r = p::project(&s, &[], &Limits::default(), &select, &p::Options::default()).unwrap();
    let v = view(&r);
    assert_eq!(v.topology, p::ContourTopology::Circular);
    assert_eq!(v.nicks.len(), 2);
    assert!(v.nicks.iter().any(|n| n.strand == p::ReferenceStrand::Top
        && n.boundary == 0
        && n.three_prime == s.strands[0].three_prime
        && n.five_prime == FivePrime::Hydroxyl));
    s.strands[0].three_prime = ThreePrime::Unknown;
    let select = selection(
        &s,
        &[],
        0,
        Call {
            strand: 0,
            index: 0,
        },
    );
    let r = p::project(&s, &[], &Limits::default(), &select, &p::Options::default()).unwrap();
    assert!(!view(&r).fully_known_terminal_chemistry);
    assert!(
        view(&r)
            .nicks
            .iter()
            .any(|n| n.three_prime == ThreePrime::Unknown)
    );
    let strict = p::Options {
        require_known_terminal_chemistry: true,
        ..p::Options::default()
    };
    let r = p::project(&s, &[], &Limits::default(), &select, &strict).unwrap();
    assert_eq!(
        r.outcome,
        p::Outcome::Withheld(p::Withheld::UnknownTerminalChemistry)
    );
}
