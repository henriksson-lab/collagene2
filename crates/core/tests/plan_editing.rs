//! Graph edits retain producer identity; scientific recipe validity is separate.
use collagene_core::{
    analysis::Primer,
    chemistry::{EndRepairMethod, EndTreatment},
    plan_editing as edit,
    verification::{Operation, OrientedSource, Source, operation_sources, operation_sources_mut},
};
fn input(id: &str) -> Source {
    Source::Input(id.into())
}
fn step(i: usize) -> Source {
    Source::Step(i)
}
fn inputs() -> Vec<String> {
    vec!["original-a".into(), "original-b".into()]
}
fn graph() -> Vec<Operation> {
    vec![
        Operation::Pcr {
            template: input("original-a"),
            forward: Primer::new("A forward", "ACGT", "GG").unwrap(),
            reverse: Primer::new("A reverse", "TGCA", "CC").unwrap(),
        },
        Operation::Digest {
            template: input("original-b"),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 2,
        },
        Operation::TreatEnds {
            template: step(0),
            treatment: EndTreatment::Phosphorylate,
        },
        Operation::Ligate {
            fragments: vec![
                OrientedSource {
                    source: step(2),
                    reverse: false,
                },
                OrientedSource {
                    source: step(1),
                    reverse: true,
                },
            ],
            circular: true,
            allow_nicked: false,
            options: Default::default(),
        },
        Operation::RepairEnds {
            template: step(3),
            method: EndRepairMethod::Blunt,
        },
    ]
}
fn label(op: &Operation) -> &str {
    match op {
        Operation::Pcr { .. } => "A",
        Operation::Digest { .. } => "B",
        Operation::TreatEnds { .. } => "C",
        Operation::Ligate { .. } => "D",
        Operation::RepairEnds { .. } => "E",
        _ => "extra",
    }
}
fn edges(ops: &[Operation]) -> Vec<String> {
    let mut edges = ops
        .iter()
        .flat_map(|op| {
            operation_sources(op).into_iter().map(move |source| {
                format!(
                    "{}<-{}",
                    label(op),
                    match source {
                        Source::Input(id) => format!("input:{id}"),
                        Source::Step(index) => label(&ops[*index]).into(),
                    }
                )
            })
        })
        .collect::<Vec<_>>();
    edges.sort();
    edges
}
fn bytes(ops: &[Operation]) -> Vec<u8> {
    serde_json::to_vec(ops).unwrap()
}
#[test]
fn moving_independent_branch_preserves_literal_dag_edges_orientation_and_pcr_parameters() {
    let original = graph();
    let frozen = bytes(&original);
    let moved = edit::moved(&inputs(), &original, 1, 0).unwrap();
    assert_eq!(
        moved.iter().map(label).collect::<Vec<_>>(),
        vec!["B", "A", "C", "D", "E"]
    );
    let expected = vec![
        "A<-input:original-a",
        "B<-input:original-b",
        "C<-A",
        "D<-B",
        "D<-C",
        "E<-D",
    ];
    assert_eq!(edges(&moved), expected);
    assert_eq!(bytes(&original), frozen);
    let Operation::TreatEnds { template, .. } = &moved[2] else {
        panic!()
    };
    assert!(matches!(template, Source::Step(1)));
    let Operation::Ligate {
        fragments,
        circular,
        allow_nicked,
        ..
    } = &moved[3]
    else {
        panic!()
    };
    assert!(*circular);
    assert!(!allow_nicked);
    assert!(fragments[1].reverse);
    assert!(matches!(fragments[1].source, Source::Step(0)));
    let Operation::Pcr {
        forward, reverse, ..
    } = &moved[1]
    else {
        panic!()
    };
    assert_eq!(forward.name, "A forward");
    assert_eq!(forward.ordered_sequence(), "GGACGT");
    assert_eq!(reverse.ordered_sequence(), "CCTGCA");
    let restored = edit::moved(&inputs(), &moved, 0, 1).unwrap();
    assert_eq!(bytes(&restored), frozen);
}
#[test]
fn moving_consumer_before_producer_or_producer_after_consumer_fails_atomically() {
    let original = graph();
    let frozen = bytes(&original);
    assert!(edit::moved(&inputs(), &original, 3, 0).is_err());
    assert!(edit::moved(&inputs(), &original, 0, 4).is_err());
    assert!(edit::moved(&inputs(), &original, usize::MAX, 0).is_err());
    assert!(edit::moved(&inputs(), &original, 0, original.len()).is_err());
    assert_eq!(bytes(&original), frozen);
    assert_eq!(
        bytes(&edit::moved(&inputs(), &original, 2, 2).unwrap()),
        frozen
    );
}
#[test]
fn deletion_refuses_live_producer_and_shifts_independent_remaining_edges() {
    let original = graph();
    let frozen = bytes(&original);
    assert!(edit::removed(&inputs(), &original, 1).is_err());
    assert!(edit::removed(&inputs(), &original, 0).is_err());
    assert_eq!(bytes(&original), frozen);
    // Remove B from a graph where that branch is explicitly no longer consumed.
    let mut independent = graph();
    let Operation::Ligate { fragments, .. } = &mut independent[3] else {
        panic!()
    };
    fragments.pop();
    let removed = edit::removed(&inputs(), &independent, 1).unwrap();
    assert_eq!(
        removed.iter().map(label).collect::<Vec<_>>(),
        vec!["A", "C", "D", "E"]
    );
    assert_eq!(
        edges(&removed),
        vec!["A<-input:original-a", "C<-A", "D<-C", "E<-D"]
    );
    let sink_removed = edit::removed(&inputs(), &original, 4).unwrap();
    assert_eq!(sink_removed.len(), 4);
    assert_eq!(
        edges(&sink_removed),
        vec![
            "A<-input:original-a",
            "B<-input:original-b",
            "C<-A",
            "D<-B",
            "D<-C"
        ]
    );
    let empty = edit::removed(&inputs(), &original[..1], 0).unwrap();
    assert!(empty.is_empty());
    edit::validate_steps(&[], &empty).unwrap();
}
#[test]
fn insertion_uses_pre_insertion_indices_and_does_not_retarget_existing_consumers() {
    let original = graph();
    let frozen = bytes(&original);
    let new = Operation::TreatEnds {
        template: step(1),
        treatment: EndTreatment::Dephosphorylate,
    };
    let inserted = edit::inserted(&inputs(), &original, 2, new).unwrap();
    let Operation::TreatEnds {
        template,
        treatment,
    } = &inserted[2]
    else {
        panic!()
    };
    assert!(matches!(template, Source::Step(1)));
    assert_eq!(*treatment, EndTreatment::Dephosphorylate);
    let Operation::Ligate { fragments, .. } = &inserted[4] else {
        panic!()
    };
    assert!(matches!(fragments[0].source, Source::Step(3)));
    assert!(matches!(fragments[1].source, Source::Step(1)));
    let Operation::RepairEnds { template, .. } = &inserted[5] else {
        panic!()
    };
    assert!(matches!(template, Source::Step(4)));
    assert_eq!(bytes(&original), frozen);
    edit::validate_steps(&inputs(), &inserted).unwrap();
    assert!(
        edit::inserted(
            &inputs(),
            &original,
            2,
            Operation::TreatEnds {
                template: step(2),
                treatment: EndTreatment::Phosphorylate
            }
        )
        .is_err(),
        "the old step2 has not been produced before insertion position2"
    );
    let at_start = edit::inserted(
        &inputs(),
        &original,
        0,
        Operation::Digest {
            template: input("original-a"),
            enzymes: vec![],
            fragment_index: 0,
        },
    )
    .unwrap();
    let Operation::TreatEnds { template, .. } = &at_start[3] else {
        panic!()
    };
    assert!(matches!(template, Source::Step(1)));
    let appended = edit::inserted(
        &inputs(),
        &original,
        5,
        Operation::TreatEnds {
            template: step(4),
            treatment: EndTreatment::Phosphorylate,
        },
    )
    .unwrap();
    assert_eq!(appended.len(), 6);
    assert_eq!(bytes(&appended[..5]), frozen);
}
#[test]
fn replacement_retains_node_identity_and_separates_graph_edits_from_scientific_claims() {
    let original = graph();
    let frozen = bytes(&original);
    let changed = edit::replaced(
        &inputs(),
        &original,
        1,
        Operation::Digest {
            template: input("original-a"),
            enzymes: vec![],
            fragment_index: 999,
        },
    )
    .unwrap();
    // Empty enzyme settings are an unfinished editable declaration. The graph
    // helper must neither simulate it nor claim its chemistry is valid.
    let Operation::Digest {
        enzymes,
        fragment_index,
        ..
    } = &changed[1]
    else {
        panic!()
    };
    assert!(enzymes.is_empty());
    assert_eq!(*fragment_index, 999);
    let Operation::Ligate { fragments, .. } = &changed[3] else {
        panic!()
    };
    assert!(matches!(fragments[1].source, Source::Step(1)));
    assert!(
        edit::replaced(
            &inputs(),
            &original,
            1,
            Operation::TreatEnds {
                template: step(1),
                treatment: EndTreatment::Phosphorylate
            }
        )
        .is_err()
    );
    assert!(
        edit::replaced(
            &inputs(),
            &original,
            1,
            Operation::Digest {
                template: input("missing"),
                enzymes: vec!["EcoRI".into()],
                fragment_index: 0
            }
        )
        .is_err()
    );
    assert_eq!(bytes(&original), frozen);
}
#[test]
fn borrowed_ids_unknown_inputs_forward_edges_and_native_graph_caps_are_explicit() {
    edit::validate_steps(&[], &[]).unwrap();
    for ids in [
        vec!["".into()],
        vec![" ".into()],
        vec!["a".into(), "a".into()],
        vec!["x".repeat(129)],
    ] {
        assert!(edit::validate_steps(&ids, &[]).is_err());
    }
    assert!(edit::validate_steps(&["original-a".into()], &graph()).is_err());
    let invalid = vec![Operation::TreatEnds {
        template: step(0),
        treatment: EndTreatment::Phosphorylate,
    }];
    assert!(edit::validate_steps(&inputs(), &invalid).is_err());
    let large = vec![
        Operation::Digest {
            template: input("original-a"),
            enzymes: vec![],
            fragment_index: 0
        };
        1000
    ];
    edit::validate_steps(&inputs(), &large).unwrap();
    assert!(edit::inserted(&inputs(), &large, 1000, large[0].clone()).is_err());
    assert!(edit::validate_steps(&inputs(), &vec![large[0].clone(); 1001]).is_err());
    let references = vec![Operation::Overlap {
        fragments: vec![input("original-a"); 1001],
        overlaps: vec![],
        circular: false,
    }];
    assert!(edit::validate_steps(&inputs(), &references).is_err());
}
#[test]
fn shared_mutable_visitor_edits_actual_ordered_sources_without_touching_parameters() {
    let mut operation = Operation::Ligate {
        fragments: vec![
            OrientedSource {
                source: input("original-a"),
                reverse: true,
            },
            OrientedSource {
                source: input("original-b"),
                reverse: false,
            },
        ],
        circular: false,
        allow_nicked: true,
        options: collagene_core::ligation_outcomes::Options {
            max_products: 7,
            max_path_checks: 9,
            ..Default::default()
        },
    };
    for source in operation_sources_mut(&mut operation) {
        if let Source::Input(id) = source {
            *source = Source::Step(if id == "original-a" { 2 } else { 4 });
        }
    }
    let Operation::Ligate {
        fragments,
        circular,
        allow_nicked,
        options,
    } = &operation
    else {
        panic!()
    };
    assert!(fragments[0].reverse);
    assert!(!fragments[1].reverse);
    assert!(!circular);
    assert!(*allow_nicked);
    assert_eq!(options.max_products, 7);
    assert_eq!(options.max_path_checks, 9);
    assert!(matches!(
        operation_sources(&operation).as_slice(),
        [Source::Step(2), Source::Step(4)]
    ));
}

#[test]
fn replacement_can_repair_one_invalid_imported_dependency_without_changing_other_edges() {
    let ids = inputs();
    let mut broken = graph();
    *operation_sources_mut(&mut broken[0])[0] = Source::Step(4);
    assert!(edit::validate_steps(&ids, &broken).is_err());
    let fixed = edit::replaced(&ids, &broken, 0, graph()[0].clone()).unwrap();
    edit::validate_steps(&ids, &fixed).unwrap();
    assert_eq!(
        serde_json::to_value(&fixed).unwrap(),
        serde_json::to_value(graph()).unwrap()
    );
    assert!(edit::replaced(&ids, &broken, 1, graph()[1].clone()).is_err());
    assert!(matches!(operation_sources(&broken[0])[0], Source::Step(4)));
}
