//! Independent source-choice, exact-event and literal-target fixtures.
use collagene_core::{
    cloning::assembly::recipe::Preparation,
    documents::workspace::Workspace,
    material::lineage::EventBinding,
    planning::combinatorial::{
        self as design, Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct,
        SelectedMatrix, TargetDeclaration,
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, Design) {
    let w = Workspace::from_documents(
        "source alternatives",
        vec![
            Molecule::dna("backbone A", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("backbone B", "GGCC", Topology::Linear).unwrap(),
            Molecule::dna("insert C", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("insert D", "CCAA", Topology::Linear).unwrap(),
            Molecule::dna("independent wrong target", "GGTT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let choice = |id: &str, index| Choice {
        anneal_bottom: None,
        id: id.into(),
        label: id.into(),
        source: binding(&w, index),
        selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            start: 0,
            end: 4,
            reverse: false,
        },
        preparation_reverse: false,
        generation: None,
        preparation: Some(Preparation::WholeLinear),
        treatment: None,
    };
    let d = Design {
        material_bindings: vec![],
        id: "library-a".into(),
        name: "Two by two library".into(),
        circular: false,
        bins: vec![
            Bin {
                id: "backbone".into(),
                label: "Backbone".into(),
                choices: vec![choice("backbone-a", 0), choice("backbone-b", 1)],
            },
            Bin {
                id: "insert".into(),
                label: "Insert".into(),
                choices: vec![choice("insert-c", 2), choice("insert-d", 3)],
            },
        ],
        method: Method::Ligation {
            allow_nicked: false,
        },
    };
    (w, d)
}
fn construct() -> SelectedConstruct {
    SelectedConstruct {
        overrides: None,
        id: "construct-ac".into(),
        name: "AC construct".into(),
        choices: vec![
            ChoiceSelection {
                bin_id: "backbone".into(),
                choice_id: "backbone-a".into(),
            },
            ChoiceSelection {
                bin_id: "insert".into(),
                choice_id: "insert-c".into(),
            },
        ],
        target: TargetDeclaration::OriginalRegions {
            substitutions: vec![],
            edits: vec![],
            junction_edits: vec![],
        },
    }
}
#[test]
fn cartesian_rows_are_choices_not_reaction_products() {
    let (_, d) = fixture();
    let page = design::matrix_page(&d, 0, 4).unwrap();
    assert_eq!(page.total, 4);
    let pairs = page
        .combinations
        .iter()
        .map(|row| {
            row.choices
                .iter()
                .map(|c| c.choice_id.as_str())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        pairs,
        vec![
            vec!["backbone-a", "insert-c"],
            vec!["backbone-a", "insert-d"],
            vec!["backbone-b", "insert-c"],
            vec!["backbone-b", "insert-d"]
        ]
    );
    let partial = design::matrix_page(&d, 2, 1).unwrap();
    assert_eq!(partial.combinations[0].id, page.combinations[2].id);
    let mut renamed = d.clone();
    renamed.name = "Renamed library".into();
    renamed.bins[0].choices[0].label = "Display rename".into();
    assert_eq!(
        design::matrix_page(&renamed, 0, 1).unwrap().combinations[0].id,
        page.combinations[0].id
    );
}
#[test]
fn original_contributions_declare_targets_independently_of_preparations() {
    let (w, d) = fixture();
    let row = construct();
    let resolved = design::resolve(&d, &row, &w).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), "AACCTTAA");
    assert!(resolved.target_binding.is_none());
    assert!(resolved.intent.is_some());
    assert_eq!(
        resolved
            .bindings
            .iter()
            .map(|b| b.event.clone())
            .collect::<Vec<_>>(),
        vec![binding(&w, 0), binding(&w, 2)]
    );
    let mut different = row.clone();
    different.target = TargetDeclaration::Existing(binding(&w, 4));
    let wrong = design::resolve(&d, &different, &w).unwrap();
    assert_eq!(wrong.target.sequence.as_str(), "GGTT");
    assert_eq!(wrong.target_binding, Some(binding(&w, 4)));
    assert_ne!(resolved.declaration_revision, wrong.declaration_revision);
    assert_eq!(resolved.combination_id, wrong.combination_id);
    // Resolve compiles a declaration; it deliberately supplies no success badge.
}

#[test]
fn internal_edit_changes_only_declared_row_and_never_repairs_physical_preparation() {
    use collagene_core::{
        chemistry::Terminus, cloning::assembly::edit_intent::Edit, verification::evidence::Status,
    };
    let (w, mut d) = fixture();
    for bin in &mut d.bins {
        for choice in &mut bin.choices {
            choice.preparation = Some(Preparation::DeclaredDuplex {
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED,
            });
        }
    }
    let original = construct();
    let mut edited = original.clone();
    edited.id = "construct-edited".into();
    edited.target = TargetDeclaration::OriginalRegions {
        junction_edits: vec![],
        substitutions: vec![],
        edits: vec![Edit {
            role: "fragment-backbone".into(),
            start: 2,
            end: 4,
            expected: "CC".into(),
            replacement: "TAG".into(),
        }],
    };
    let source_bytes = serde_json::to_vec(&w).unwrap();
    let prepared = design::evaluate(&d, &original, &w).unwrap();
    assert_eq!(prepared.resolved.target.sequence.as_str(), "AACCTTAA");
    assert_eq!(
        prepared.verification().status,
        Status::ConsistentUnderAssumptions
    );
    let changed = design::evaluate(&d, &edited, &w).unwrap();
    assert_eq!(changed.resolved.target.sequence.as_str(), "AATAGTTAA");
    assert_eq!(changed.verification().status, Status::Invalid);
    assert_eq!(
        changed
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "AACCTTAA"
    );
    assert_ne!(
        prepared.resolved.declaration_revision,
        changed.resolved.declaration_revision
    );
    assert_eq!(
        design::resolve(&d, &original, &w)
            .unwrap()
            .target
            .sequence
            .as_str(),
        "AACCTTAA"
    );
    assert_eq!(serde_json::to_vec(&w).unwrap(), source_bytes);
    let TargetDeclaration::OriginalRegions { edits, .. } = &mut edited.target else {
        unreachable!()
    };
    edits[0].expected = "GG".into();
    assert!(design::resolve(&d, &edited, &w).is_err());
}

#[test]
fn current_library_declarations_require_explicit_internal_edit_vectors() {
    let current = serde_json::to_value(construct()).unwrap();
    assert!(serde_json::from_value::<SelectedConstruct>(current.clone()).is_ok());
    for field in ["substitutions", "edits"] {
        let mut malformed = current.clone();
        malformed["target"]["OriginalRegions"]
            .as_object_mut()
            .unwrap()
            .remove(field);
        assert!(serde_json::from_value::<SelectedConstruct>(malformed).is_err());
    }
}
#[test]
fn wrong_event_owner_and_duplicate_ids_are_rejected() {
    let (w, d) = fixture();
    let mut owner = d.clone();
    owner.bins[0].choices[0].source.document_id = w.documents[1].id.clone();
    assert!(design::resolve(&owner, &construct(), &w).is_err());
    let mut duplicate = d.clone();
    duplicate.bins[1].choices[0].id = "backbone-a".into();
    assert!(duplicate.validate().is_err());
    let row = construct();
    let mut duplicate_row = row.clone();
    duplicate_row.id = "different-row-id".into();
    assert!(
        design::validate_matrix(&SelectedMatrix {
            design: d,
            constructs: vec![row, duplicate_row]
        })
        .is_ok()
    );
}
#[test]
fn unselected_unfinished_choice_and_unrelated_rich_document_do_not_block_selected_construct() {
    let (w, mut d) = fixture();
    d.bins[0].choices[1].preparation = None;
    let mut docs = w
        .documents
        .iter()
        .map(|doc| w.revisions[&doc.head].clone())
        .collect::<Vec<_>>();
    let mut unrelated =
        Molecule::dna("unrelated rich reference", "ACGT", Topology::Linear).unwrap();
    unrelated.features = (0..10_000)
        .map(|i| Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, 1),
            qualifiers: vec![("label".into(), Some(format!("feature {i}")))],
            needs_review: false,
        })
        .collect();
    docs.push(unrelated);
    let rich = Workspace::from_documents("rich workspace", docs, 0).unwrap();
    let resolved = design::resolve(&d, &construct(), &rich).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), "AACCTTAA");
    let mut unfinished = construct();
    unfinished.choices[0].choice_id = "backbone-b".into();
    assert!(design::resolve(&d, &unfinished, &rich).is_err());
}

#[test]
fn checked_mixed_radix_pages_and_stable_operation_ids_keep_distinct_choices() {
    assert_eq!(design::cartesian_size(&[2, 3]).unwrap(), 6);
    assert_eq!(
        design::cartesian_page(&[2, 3], 2, 3).unwrap(),
        vec![vec![0, 2], vec![1, 0], vec![1, 1]]
    );
    assert!(design::cartesian_page(&[2, 3], 6, 1).unwrap().is_empty());
    for counts in [vec![], vec![0], vec![65], vec![64; 5], vec![1; 13]] {
        assert!(design::cartesian_size(&counts).is_err());
    }
    assert!(design::cartesian_page(&[2], u128::MAX, 1).is_err());
    assert!(design::cartesian_page(&[2], 0, 65).is_err());
    let a = design::source_role_ids("library", "construct-a", "bin", "choice-a").unwrap();
    let b = design::source_role_ids("library", "construct-a", "bin", "choice-b").unwrap();
    let c = design::source_role_ids("library", "construct-b", "bin", "choice-a").unwrap();
    assert_ne!(a.0, b.0);
    assert_ne!(a.0, c.0);
    assert_eq!(a.1, "fragment-bin");
    assert!(a.0.len() < 128);
}

#[test]
fn fresh_existing_verifier_preserves_a_wrong_target_and_changed_contribution() {
    use collagene_core::{chemistry::Terminus, verification::evidence::Status};
    let (w, mut d) = fixture();
    for bin in &mut d.bins {
        for choice in &mut bin.choices {
            choice.preparation = Some(Preparation::DeclaredDuplex {
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED,
            });
        }
    }
    let row = construct();
    let evaluated = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(evaluated.resolved.target.sequence.as_str(), "AACCTTAA");
    assert_eq!(
        evaluated.verification().status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        evaluated.verification().findings
    );
    let compiled = design::compile(&d, &row, &w).unwrap();
    assert_eq!(
        compiled.plan.target_revision,
        evaluated.resolved.target.revision()
    );
    assert_eq!(
        compiled.plan.inputs[0].id,
        design::source_role_ids(&d.id, &row.id, "backbone", "backbone-a")
            .unwrap()
            .0
    );
    let mut wrong = row.clone();
    wrong.target = TargetDeclaration::Existing(binding(&w, 4));
    let mismatch = design::evaluate(&d, &wrong, &w).unwrap();
    assert_eq!(mismatch.verification().status, Status::Invalid);
    assert_eq!(mismatch.resolved.target.sequence.as_str(), "GGTT");
    assert_eq!(
        mismatch
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "AACCTTAA"
    );
    let mut cropped = d.clone();
    let collagene_core::cloning::assembly::intent::Selector::ReferenceRange { start, .. } =
        &mut cropped.bins[0].choices[0].selector
    else {
        unreachable!()
    };
    *start = 1;
    let changed = design::evaluate(&cropped, &row, &w).unwrap();
    assert_eq!(changed.resolved.target.sequence.as_str(), "ACCTTAA");
    assert_eq!(changed.verification().status, Status::Invalid);
    assert_eq!(
        changed.resolved.combination_id,
        evaluated.resolved.combination_id
    );
    assert_ne!(
        changed.resolved.declaration_revision,
        evaluated.resolved.declaration_revision
    );
}

#[test]
fn historical_choice_keeps_its_event_without_rewinding_later_source_heads() {
    use collagene_core::documents::history::DocumentHistory;
    let (w, mut d) = fixture();
    let row = construct();
    let original = design::resolve(&d, &row, &w).unwrap();
    let mut docs = DocumentHistory::restore_all(&w).unwrap();
    docs[0]
        .commit(
            Molecule::dna("edited source", "TCCC", Topology::Linear).unwrap(),
            "Later source change",
        )
        .unwrap();
    let mut changed = w.clone();
    changed.documents.clear();
    changed.history.clear();
    for doc in docs {
        doc.contribute(&mut changed).unwrap();
    }
    changed.validate().unwrap();
    let head = changed.documents[0].head.clone();
    let historical = design::resolve(&d, &row, &changed).unwrap();
    assert_eq!(historical.target.sequence.as_str(), "AACCTTAA");
    assert_eq!(historical.bindings[0].event, binding(&w, 0));
    assert_eq!(changed.documents[0].head, head);
    assert_eq!(
        historical.declaration_revision,
        original.declaration_revision
    );
    d.bins[0].choices[0].source = binding(&changed, 0);
    let updated = design::resolve(&d, &row, &changed).unwrap();
    assert_eq!(updated.target.sequence.as_str(), "TCCCTTAA");
    assert_ne!(
        updated.declaration_revision,
        historical.declaration_revision
    );
    assert_eq!(updated.combination_id, historical.combination_id);
}

#[test]
fn circular_unrolled_contribution_and_explicit_reverse_keep_source_geometry() {
    use collagene_core::sequence::analysis::Primer;
    let m = Molecule::dna("circle", "ACGTAA", Topology::Circular).unwrap();
    let w = Workspace::from_documents("circle choices", vec![m], 0).unwrap();
    let d = Design {
        material_bindings: vec![],
        id: "circle-design".into(),
        name: "Circle design".into(),
        circular: false,
        bins: vec![Bin {
            id: "part".into(),
            label: "Part".into(),
            choices: vec![Choice {
                anneal_bottom: None,
                id: "crossing".into(),
                label: "Crossing".into(),
                source: binding(&w, 0),
                selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: 4,
                    end: 8,
                    reverse: true,
                },
                preparation_reverse: true,
                generation: None,
                preparation: Some(Preparation::ExactPcr {
                    forward: Primer::new("F", "AAAC", "").unwrap(),
                    reverse: Primer::new("R", "GTTT", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                }),
                treatment: None,
            }],
        }],
        method: Method::Ligation {
            allow_nicked: false,
        },
    };
    let row = SelectedConstruct {
        overrides: None,
        id: "reverse-crossing".into(),
        name: "Reverse crossing".into(),
        choices: vec![ChoiceSelection {
            bin_id: "part".into(),
            choice_id: "crossing".into(),
        }],
        target: TargetDeclaration::OriginalRegions {
            substitutions: vec![],
            edits: vec![],
            junction_edits: vec![],
        },
    };
    let resolved = design::resolve(&d, &row, &w).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), "GTTT");
    let mut overrun = d.clone();
    let collagene_core::cloning::assembly::intent::Selector::ReferenceRange { end, .. } =
        &mut overrun.bins[0].choices[0].selector
    else {
        unreachable!()
    };
    *end = 11;
    assert!(design::resolve(&overrun, &row, &w).is_err());
    let mut unsupported = d.clone();
    unsupported.method = Method::OrderedTypeIis {
        enzyme: collagene_core::cloning::golden_gate::GoldenGateEnzyme::BsaI,
        allow_nicked: false,
    };
    assert!(design::resolve(&unsupported, &row, &w).is_err());
}
