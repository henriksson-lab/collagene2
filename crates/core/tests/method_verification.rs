use collagene_core::{
    analysis::Primer,
    chemistry::{AnnealOptions, EndTreatment, Terminus, ThreePrime},
    overlap_profiles::ProfileOptions,
    tailing::TailingBase,
    verification::*,
    *,
};
fn dna(name: &str, s: &str, t: Topology) -> Molecule {
    Molecule::dna(name, s, t).unwrap()
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn source(s: &str) -> Source {
    Source::Input(s.into())
}
fn step(n: usize) -> Source {
    Source::Step(n)
}
fn anneal(top: &str, bottom: &str, chemistry: Terminus) -> Operation {
    Operation::AnnealOligos {
        top: source(top),
        bottom: source(bottom),
        options: AnnealOptions {
            min_paired_bases: 4,
            top_chemistry: chemistry,
            bottom_chemistry: chemistry,
            ..Default::default()
        },
        candidate_index: 0,
    }
}
fn ta_plan(base: TailingBase, insert: TailingBase, chemistry: Terminus) -> DesignPlan {
    let target = dna("target", "ACGCTGGTAA", Topology::Circular);
    DesignPlan {
        name: "literal TA".into(),
        target_revision: target.revision(),
        inputs: vec![
            input("vt", dna("vt", "ACGC", Topology::Linear)),
            input("vb", dna("vb", "GCGT", Topology::Linear)),
            input("it", dna("it", "GGTA", Topology::Linear)),
            input("ib", dna("ib", "TACC", Topology::Linear)),
        ],
        steps: vec![
            anneal("vt", "vb", chemistry),
            Operation::Tail3Prime {
                template: step(0),
                base,
            },
            anneal("it", "ib", chemistry),
            Operation::Tail3Prime {
                template: step(2),
                base: insert,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: step(1),
                        reverse: false,
                    },
                    OrientedSource {
                        source: step(3),
                        reverse: false,
                    },
                ],
                circular: true,
                allow_nicked: false,
            },
        ],
    }
}
#[test]
fn literal_ta_junctions_include_both_added_bases_and_are_independently_bound() {
    let target = dna("target", "ACGCTGGTAA", Topology::Circular);
    let p = ta_plan(TailingBase::T, TailingBase::A, Terminus::CLEAVED);
    let r = verify(&p, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products[1].sequence.as_str(), "ACGCT");
    assert_eq!(r.products[3].sequence.as_str(), "GGTAA");
    assert_eq!(r.products[4].sequence.len(), 10);
    assert!(
        r.findings
            .iter()
            .any(|s| s.contains("no source-template ancestry"))
    );
    let wrong = dna("wrong", "ACGCGGTA", Topology::Circular);
    let mut wrong_plan = p.clone();
    wrong_plan.target_revision = wrong.revision();
    assert_eq!(verify(&wrong_plan, &wrong).status, Status::Invalid);
    assert_eq!(verify(&p, &wrong).status, Status::Invalid);
    let mut mismatch = p.clone();
    mismatch.steps[3] = Operation::Tail3Prime {
        template: step(2),
        base: TailingBase::C,
    };
    assert_eq!(verify(&mismatch, &target).status, Status::Invalid);
}
#[test]
fn phosphorylation_and_orientation_are_explicit_in_ta_gc_candidates() {
    let target = dna("target", "ACGCTGGTAA", Topology::Circular);
    let mut p = ta_plan(TailingBase::T, TailingBase::A, Terminus::PCR);
    assert_eq!(verify(&p, &target).status, Status::Invalid);
    p.steps.pop();
    p.steps.extend([
        Operation::TreatEnds {
            template: step(1),
            treatment: EndTreatment::Phosphorylate,
        },
        Operation::TreatEnds {
            template: step(3),
            treatment: EndTreatment::Phosphorylate,
        },
        Operation::Ligate {
            fragments: vec![
                OrientedSource {
                    source: step(4),
                    reverse: false,
                },
                OrientedSource {
                    source: step(5),
                    reverse: false,
                },
            ],
            circular: true,
            allow_nicked: false,
        },
    ]);
    assert_eq!(
        verify(&p, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let reverse = dna("reverse", "ACGCTTACCA", Topology::Circular);
    p.target_revision = reverse.revision();
    if let Operation::Ligate { fragments, .. } = p.steps.last_mut().unwrap() {
        fragments[1].reverse = true;
    }
    let r = verify(&p, &reverse);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    let gc = dna("gc", "ACGCCGGTAG", Topology::Circular);
    let mut g = ta_plan(TailingBase::C, TailingBase::G, Terminus::CLEAVED);
    g.target_revision = gc.revision();
    assert_eq!(verify(&g, &gc).status, Status::ConsistentUnderAssumptions);
}
#[test]
fn preview_metadata_cannot_grant_chemistry_and_unknown_or_blocked_tails_remain_incomplete() {
    let source_m = dna("plain", "ACGC", Topology::Linear);
    let target = dna("tailed", "ACGCA", Topology::Linear);
    let mut plan = DesignPlan {
        name: "unknown".into(),
        target_revision: target.revision(),
        inputs: vec![input("s", source_m)],
        steps: vec![Operation::Tail3Prime {
            template: source("s"),
            base: TailingBase::A,
        }],
    };
    assert_eq!(verify(&plan, &target).status, Status::Incomplete);
    let mut p = ta_plan(
        TailingBase::T,
        TailingBase::A,
        Terminus {
            three_prime: ThreePrime::Blocked,
            ..Terminus::CLEAVED
        },
    );
    let target_circle = dna("target", "ACGCTGGTAA", Topology::Circular);
    assert_eq!(verify(&p, &target_circle).status, Status::Incomplete);
    p.steps[0] = anneal("vt", "vb", Terminus::CLEAVED);
    let r = verify(&p, &target_circle);
    let preview = r.products[1].clone();
    plan.inputs = vec![input("s", preview)];
    assert_eq!(verify(&plan, &target).status, Status::Incomplete);
}
#[test]
fn digest_and_unique_pcr_tailing_previews_replay_real_sources() {
    let circle = dna("circle", "GATATCAAAA", Topology::Circular);
    let recipe = TailSubstrate::Digest {
        enzymes: vec!["EcoRV".into()],
        fragment_index: 0,
        blunt_before: false,
    };
    let r = simulate_tailing(&circle, &recipe, TailingBase::T).unwrap();
    assert_eq!(r.tailing.source.top.as_str(), "ATCAAAAGAT");
    assert_eq!(r.tailing.duplex.top.as_str(), "ATCAAAAGATT");
    assert_eq!(r.tailing.duplex.bottom_start, -1);
    r.tailing.validate().unwrap();
    let pcr = dna("template", "ACGACCCCGGGGTAAC", Topology::Linear);
    let recipe = TailSubstrate::ExactPcr {
        forward: Primer::new("f", "ACGA", "").unwrap(),
        reverse: Primer::new("r", "GTTA", "").unwrap(),
    };
    let r = simulate_tailing(&pcr, &recipe, TailingBase::A).unwrap();
    assert_eq!(r.tailing.duplex.top.as_str(), "ACGACCCCGGGGTAACA");
    assert_eq!(
        r.tailing.duplex.left.five_prime,
        chemistry::FivePrime::Hydroxyl
    );
}
#[test]
fn oriented_profiles_reconstruct_and_map_source_features_without_enforcing_recommendations() {
    let first = dna("first", "AAAACC", Topology::Linear);
    let mut second = dna("reverse source", "AACCGG", Topology::Linear);
    second.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Range {
            start: 0,
            end: 2,
            before: false,
            after: false,
        },
        qualifiers: vec![("label".into(), Some("original AA".into()))],
        needs_review: false,
    });
    let target = dna("target", "AAAACCGGTT", Topology::Linear);
    let plan = DesignPlan {
        name: "oriented short homology".into(),
        target_revision: target.revision(),
        inputs: vec![input("a", first), input("b", second)],
        steps: vec![Operation::ProfiledOverlap {
            fragments: vec![
                OrientedSource {
                    source: source("a"),
                    reverse: false,
                },
                OrientedSource {
                    source: source("b"),
                    reverse: true,
                },
            ],
            overlaps: vec![2],
            circular: false,
            options: ProfileOptions {
                scan_repeated_overlaps: false,
                ..Default::default()
            },
        }],
    };
    let r = verify(&plan, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert!(r.findings.iter().any(|s| s.contains("advisory")));
    let f = r.products[0]
        .features
        .iter()
        .find(|f| f.label() == "original AA")
        .unwrap();
    assert_eq!(f.location.spans(), vec![(8, 10, true)]);
    assert_eq!(f.location.extract(&r.products[0].sequence).unwrap(), "AA");
    let wrong = dna("wrong", "AAAACCGGTA", Topology::Linear);
    let mut p = plan;
    p.target_revision = wrong.revision();
    assert_eq!(verify(&p, &wrong).status, Status::Invalid);
}
#[test]
fn json_cannot_silently_discard_requested_tailing_or_primer_chemistry() {
    let mut value = serde_json::to_value(Operation::Tail3Prime {
        template: source("s"),
        base: TailingBase::A,
    })
    .unwrap();
    value["Tail3Prime"]["strand"] = serde_json::json!("Top");
    assert!(serde_json::from_value::<Operation>(value).is_err());
    let mut primer = serde_json::to_value(Primer::new("f", "ACGA", "").unwrap()).unwrap();
    primer["five_prime"] = serde_json::json!("Phosphate");
    assert!(serde_json::from_value::<Primer>(primer).is_err());
    let mut seq = serde_json::to_value(Sequence::dna("ACGA").unwrap()).unwrap();
    seq["modification"] = serde_json::json!("blocked");
    assert!(serde_json::from_value::<Sequence>(seq).is_err());
}
