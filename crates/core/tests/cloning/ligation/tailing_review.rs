//! Independent verifier review: literal two-strand products and annotation origins.
use collagene_core::{
    chemistry::{AnnealOptions, DuplexSummary, FivePrime, Nick, Strand, Terminus},
    cloning::ligation::tailing::TailingBase,
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    verification::{self, evidence::Status},
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn input(id: &str, sequence: &str, label: Option<&str>) -> PlanInput {
    let mut molecule = Molecule::dna(id, sequence, Topology::Linear).unwrap();
    if let Some(label) = label {
        molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, sequence.len()),
            qualifiers: vec![("label".into(), Some(label.into()))],
            needs_review: false,
        });
    }
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn source(id: &str) -> Source {
    Source::Input(id.into())
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
fn plan(
    target: &str,
    vector_base: TailingBase,
    insert_base: TailingBase,
    vector_chemistry: Terminus,
    insert_chemistry: Terminus,
    reverse: bool,
    allow_nicked: bool,
) -> (DesignPlan, Molecule) {
    let target =
        Molecule::dna("independently authored target", target, Topology::Circular).unwrap();
    (
        DesignPlan {
            name: "literal source-tail-ligate review".into(),
            target_revision: target.revision(),
            inputs: vec![
                input("vector_top", "ACGC", Some("vector original")),
                input("vector_bottom", "GCGT", None),
                input("insert_top", "GGTA", Some("insert original")),
                input("insert_bottom", "TACC", None),
            ],
            steps: vec![
                anneal("vector_top", "vector_bottom", vector_chemistry),
                Operation::Tail3Prime {
                    template: Source::Step(0),
                    base: vector_base,
                },
                anneal("insert_top", "insert_bottom", insert_chemistry),
                Operation::Tail3Prime {
                    template: Source::Step(2),
                    base: insert_base,
                },
                Operation::Ligate {
                    fragments: vec![
                        OrientedSource {
                            source: Source::Step(1),
                            reverse: false,
                        },
                        OrientedSource {
                            source: Source::Step(3),
                            reverse,
                        },
                    ],
                    circular: true,
                    allow_nicked,
                    options: Default::default(),
                },
            ],
        },
        target,
    )
}
#[test]
fn forward_and_reversed_ta_targets_keep_original_features_off_the_two_added_bases() {
    for (reverse, sequence) in [(false, "ACGCTGGTAA"), (true, "ACGCTTACCA")] {
        let (p, t) = plan(
            sequence,
            TailingBase::T,
            TailingBase::A,
            Terminus::CLEAVED,
            Terminus::CLEAVED,
            reverse,
            false,
        );
        let r = verification::engine::verify(&p, &t);
        assert_eq!(
            r.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            r.findings
        );
        let product = r
            .final_product()
            .expect("final declared step must reconstruct");
        assert_eq!(product.sequence.as_str(), sequence);
        assert_eq!(product.features.len(), 2);
        let vector = product
            .features
            .iter()
            .find(|f| f.label() == "vector original")
            .unwrap();
        assert_eq!(vector.location.spans(), vec![(0, 4, false)]);
        assert_eq!(vector.location.extract(&product.sequence).unwrap(), "ACGC");
        let insert = product
            .features
            .iter()
            .find(|f| f.label() == "insert original")
            .unwrap();
        assert_eq!(insert.location.spans(), vec![(5, 9, reverse)]);
        assert_eq!(insert.location.extract(&product.sequence).unwrap(), "GGTA");
        for feature in &product.features {
            assert!(
                feature
                    .location
                    .spans()
                    .iter()
                    .all(|(a, b, _)| !(*a <= 4 && 4 < *b) && !(*a <= 9 && 9 < *b))
            );
        }
        let summary: DuplexSummary =
            serde_json::from_str(&product.metadata["duplex_model"]).unwrap();
        assert!(summary.covalently_closed_in_model);
        assert!(summary.nicks.is_empty());
        assert!(!r.exhaustive_alternative_search);
        let wrong = Molecule::dna("naive no-tail target", "ACGCGGTA", Topology::Circular).unwrap();
        let mut bad = p.clone();
        bad.target_revision = wrong.revision();
        assert_eq!(
            verification::engine::verify(&bad, &wrong).status,
            Status::Invalid
        );
    }
}
#[test]
fn gc_phosphates_and_nick_acceptance_are_not_inferred_from_tails_or_target_matches() {
    let (p, t) = plan(
        "ACGCCGGTAG",
        TailingBase::C,
        TailingBase::G,
        Terminus::PCR,
        Terminus::CLEAVED,
        false,
        true,
    );
    let r = verification::engine::verify(&p, &t);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    let product = r
        .final_product()
        .expect("final declared step must reconstruct");
    let summary: DuplexSummary = serde_json::from_str(&product.metadata["duplex_model"]).unwrap();
    assert_eq!(
        summary.nicks,
        vec![
            Nick {
                strand: Strand::Bottom,
                position: 4
            },
            Nick {
                strand: Strand::Top,
                position: 0
            }
        ]
    );
    assert!(!summary.covalently_closed_in_model);
    let mut sealed_required = p.clone();
    if let Some(Operation::Ligate { allow_nicked, .. }) = sealed_required.steps.last_mut() {
        *allow_nicked = false;
    }
    assert_eq!(
        verification::engine::verify(&sealed_required, &t).status,
        Status::Invalid
    );
    let (p, t) = plan(
        "ACGCCGGTAG",
        TailingBase::C,
        TailingBase::G,
        Terminus::PCR,
        Terminus::PCR,
        false,
        true,
    );
    assert_eq!(verification::engine::verify(&p, &t).status, Status::Invalid);
    let mut unknown = Terminus::CLEAVED;
    unknown.five_prime = FivePrime::Unknown;
    let (p, t) = plan(
        "ACGCCGGTAG",
        TailingBase::C,
        TailingBase::G,
        unknown,
        Terminus::CLEAVED,
        false,
        true,
    );
    assert_eq!(
        verification::engine::verify(&p, &t).status,
        Status::Incomplete
    );
}
#[test]
fn imported_physical_metadata_cannot_establish_accessible_tailing_ends() {
    let target = Molecule::dna("target", "ACGCA", Topology::Linear).unwrap();
    for forged in [false, true] {
        let mut source = input("source", "ACGC", None);
        if forged {
            source.molecule.metadata.insert(
                "duplex_model".into(),
                "{\"ends\":\"knownOH\",\"tailing_verified\":true}".into(),
            );
            source
                .molecule
                .metadata
                .insert("tailing:base".into(), "A".into());
            source.revision = source.molecule.revision();
        }
        let p = DesignPlan {
            name: "no imported authority".into(),
            target_revision: target.revision(),
            inputs: vec![source],
            steps: vec![Operation::Tail3Prime {
                template: Source::Input("source".into()),
                base: TailingBase::A,
            }],
        };
        let r = verification::engine::verify(&p, &target);
        assert_eq!(r.status, Status::Incomplete);
        assert!((r.reconstructed_count() == 0));
        assert!(
            !r.findings
                .iter()
                .any(|f| f.message.contains("matches the separately bound target"))
        );
    }
}
#[test]
fn prior_trusted_linear_ligation_nicks_survive_tailing_with_unchanged_five_prime_state() {
    let target = Molecule::dna("nicked linear target", "ACGCGGTAA", Topology::Linear).unwrap();
    let p = DesignPlan {
        name: "retain a known nick".into(),
        target_revision: target.revision(),
        inputs: vec![
            input("a", "ACGC", None),
            input("ar", "GCGT", None),
            input("b", "GGTA", None),
            input("br", "TACC", None),
        ],
        steps: vec![
            anneal("a", "ar", Terminus::CLEAVED),
            anneal("b", "br", Terminus::PCR),
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
                circular: false,
                allow_nicked: true,
                options: Default::default(),
            },
            Operation::Tail3Prime {
                template: Source::Step(2),
                base: TailingBase::A,
            },
        ],
    };
    let r = verification::engine::verify(&p, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    let summary: DuplexSummary = serde_json::from_str(
        &r.final_product()
            .expect("final declared step must reconstruct")
            .metadata["duplex_model"],
    )
    .unwrap();
    assert_eq!(
        summary.nicks,
        vec![Nick {
            strand: Strand::Top,
            position: 4
        }]
    );
    let (left, right) = summary.ends.unwrap();
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Hydroxyl);
}
#[test]
fn reversed_compound_source_feature_preserves_biological_order_after_either_ligation_orientation() {
    for (reverse, sequence) in [(false, "ACGCTGGTAA"), (true, "ACGCTTACCA")] {
        let (mut p, t) = plan(
            sequence,
            TailingBase::T,
            TailingBase::A,
            Terminus::CLEAVED,
            Terminus::CLEAVED,
            reverse,
            false,
        );
        p.inputs[2].molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::Complement(Box::new(Location::Join(vec![
                Location::range(0, 1),
                Location::range(2, 4),
            ]))),
            qualifiers: vec![("label".into(), Some("reverse compound original".into()))],
            needs_review: false,
        });
        p.inputs[2].revision = p.inputs[2].molecule.revision();
        let r = verification::engine::verify(&p, &t);
        assert_eq!(
            r.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            r.findings
        );
        let product = r
            .final_product()
            .expect("final declared step must reconstruct");
        let feature = product
            .features
            .iter()
            .find(|f| f.label() == "reverse compound original")
            .unwrap();
        assert_eq!(feature.location.extract(&product.sequence).unwrap(), "TAC");
        assert!(
            feature
                .location
                .spans()
                .iter()
                .all(|(a, b, _)| !(*a <= 4 && 4 < *b) && !(*a <= 9 && 9 < *b))
        );
    }
}
