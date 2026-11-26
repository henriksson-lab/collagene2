//! Saved reference documents are not stock: material must be reconstructed from
//! the exact retained run and an uninterrupted name-only descendant path.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    chemistry::{Strand, Terminus},
    cloning_outputs::{self, Commit, Options},
    history::DocumentHistory,
    methylation::MethylationState,
    modified_duplex::{self, CallState, Chemical, GenerationAssumptions},
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::Workspace,
};

fn binding(workspace: &Workspace, index: usize) -> EventBinding {
    let document = &workspace.documents[index];
    EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event.clone().unwrap(),
        revision: document.head.clone(),
    }
}
fn pcr() -> Commit {
    let source = Molecule::dna(
        "Original PCR template",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent expected PCR",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Fresh saved PCR material".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("Forward", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("Reverse", "GTTACCCC", "").unwrap(),
        }],
    };
    let workspace =
        Workspace::from_documents("PCR source and independent target", vec![source, target], 1)
            .unwrap();
    cloning_outputs::finalize(
        &workspace,
        &plan,
        binding(&workspace, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&workspace, 0),
        }],
        &Options {
            name: "Named PCR fragment".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
}
fn selected(commit: &Commit, step_index: usize) -> Selection {
    Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index,
        event: commit.receipt.construct.clone(),
    }
}
fn changed_named_output(
    commit: &Commit,
    mutate: impl FnOnce(&mut DocumentHistory),
) -> (Workspace, Selection) {
    let mut workspace = commit.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&workspace).unwrap();
    let named = histories
        .iter_mut()
        .find(|history| history.id == commit.receipt.construct.document_id)
        .unwrap();
    mutate(named);
    let event = EventBinding {
        document_id: named.id.clone(),
        event_id: named.head_event().into(),
        revision: named.head_revision().into(),
    };
    workspace.documents.clear();
    workspace.revisions.clear();
    workspace.history.clear();
    for history in &histories {
        history.contribute(&mut workspace).unwrap();
    }
    workspace.validate().unwrap();
    let mut selection = selected(commit, 0);
    selection.event = event;
    (workspace, selection)
}

#[test]
fn named_pcr_reference_reconstructs_actual_hydroxyl_duplex_without_mutating_workspace() {
    let committed = pcr();
    let selection = selected(&committed, 0);
    let before = run_lineage::workspace_fingerprint(&committed.workspace).unwrap();
    assert_eq!(
        run_material::originating_products(&committed.workspace, &selection.event).unwrap(),
        vec![selection.clone()]
    );
    let fresh = run_material::reconstruct(&committed.workspace, &selection).unwrap();
    assert_eq!(fresh.reference().name, "Named PCR fragment");
    assert_eq!(fresh.reference().sequence.as_str(), "ACGACCCCGGGGTAAC");
    assert_eq!(fresh.reference().revision(), selection.event.revision);
    let tracked = fresh.tracked();
    tracked.validate().unwrap();
    assert_eq!(tracked.duplex.top.as_str(), "ACGACCCCGGGGTAAC");
    assert_eq!(tracked.duplex.bottom.as_str(), "GTTACCCCGGGGTCGT");
    assert_eq!(tracked.duplex.topology, Topology::Linear);
    assert_eq!(tracked.duplex.bottom_start, 0);
    assert_eq!(tracked.duplex.left, Terminus::PCR);
    assert_eq!(tracked.duplex.right, Terminus::PCR);
    assert!(tracked.duplex.nicks.is_empty() && tracked.generation.is_empty());
    let modifications = modified_duplex::summary(tracked).unwrap();
    assert_eq!(modifications.known_modified, 0);
    assert_eq!(modifications.unknown_eligible, 0);
    assert_eq!(fresh.report().status, Status::ConsistentUnderAssumptions);
    assert_eq!(fresh.selection().run_id, selection.run_id);
    assert_eq!(fresh.selection().step_index, 0);
    assert_eq!(fresh.selection().event, selection.event);
    assert_eq!(fresh.context_revision().len(), 64);
    assert_eq!(
        run_lineage::workspace_fingerprint(&committed.workspace).unwrap(),
        before
    );
}

#[test]
fn named_circular_pcr_ligation_retains_both_unsealed_origin_nicks() {
    use collagene_core::{chemistry::EndTreatment, verification::OrientedSource};
    let left = Molecule::dna("Left original", "TTGCAATCGTCGATGCTAGC", Topology::Linear).unwrap();
    let right = Molecule::dna("Right original", "ATCGTCGATGCTAGCGGTAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent nicked circle",
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC",
        Topology::Circular,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "One phosphorylated and one hydroxyl PCR fragment".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left.revision(),
                molecule: left.clone(),
            },
            PlanInput {
                id: "right".into(),
                revision: right.revision(),
                molecule: right.clone(),
            },
        ],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("left".into()),
                forward: Primer::new("Left forward", "TTGCAATC", "").unwrap(),
                reverse: Primer::new("Left reverse", "GCTAGCAT", "").unwrap(),
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Pcr {
                template: Source::Input("right".into()),
                forward: Primer::new("Right forward", "ATCGTCGA", "").unwrap(),
                reverse: Primer::new("Right reverse", "GTACCGCT", "").unwrap(),
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(2),
                        reverse: false,
                    },
                ],
                circular: true,
                allow_nicked: true,
            },
        ],
    };
    let workspace = Workspace::from_documents(
        "PCR and retained-nick ligation",
        vec![left, right, target],
        2,
    )
    .unwrap();
    let report = collagene_core::verification::verify(
        &plan,
        workspace.head(&workspace.active_document).unwrap(),
    );
    assert!(
        report.status != Status::Invalid && report.products.len() == 4,
        "{:?}",
        report
    );
    assert_eq!(
        report.products.last().unwrap().sequence.as_str(),
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC"
    );
    let committed = cloning_outputs::finalize(
        &workspace,
        &plan,
        binding(&workspace, 2),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: binding(&workspace, 0),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: binding(&workspace, 1),
            },
        ],
        &Options {
            name: "Named nicked circle".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let before = run_lineage::workspace_fingerprint(&committed.workspace).unwrap();
    let selection = selected(&committed, 3);
    let fresh = run_material::reconstruct(&committed.workspace, &selection).unwrap();
    let tracked = fresh.tracked();
    tracked.validate().unwrap();
    assert_eq!(fresh.reference().name, "Named nicked circle");
    assert_eq!(
        tracked.duplex.top.as_str(),
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC"
    );
    assert_eq!(
        tracked.duplex.bottom.as_str(),
        "GTACCGCTAGCATCGACGATGCTAGCATCGACGATTGCAA"
    );
    assert_eq!(tracked.duplex.topology, Topology::Circular);
    assert_eq!(tracked.duplex.bottom_start, 0);
    // Circles have no outer termini. The two actual failed bonds remain nicks.
    assert_eq!(tracked.duplex.left, Terminus::UNKNOWN);
    assert_eq!(tracked.duplex.right, Terminus::UNKNOWN);
    assert_eq!(tracked.duplex.nicks.len(), 2);
    assert!(
        tracked
            .duplex
            .nicks
            .iter()
            .any(|nick| nick.strand == Strand::Top && nick.position == 20)
    );
    assert!(
        tracked
            .duplex
            .nicks
            .iter()
            .any(|nick| nick.strand == Strand::Bottom && nick.position == 0)
    );
    assert!(!tracked.duplex.summary().unwrap().covalently_closed_in_model);
    assert!(
        run_material::originating_products(&committed.workspace, &selection.event)
            .unwrap()
            .contains(&selection)
    );
    assert_eq!(
        run_lineage::workspace_fingerprint(&committed.workspace).unwrap(),
        before
    );
}

#[test]
fn blocked_noti_reconstruction_preserves_four_literal_cpg_positions_and_generation() {
    let source = Molecule::dna("Original CpG circle", "AAGCGGCCGCAA", Topology::Circular).unwrap();
    let target = Molecule::dna(
        "Independent blocked circle",
        "AAGCGGCCGCAA",
        Topology::Circular,
    )
    .unwrap();
    let assumptions = GenerationAssumptions {
        dam: MethylationState::Unmodified,
        dcm: MethylationState::Unmodified,
        cpg: MethylationState::FullyModified,
    };
    let plan = DesignPlan {
        name: "Declared CpG generation and blocked NotI".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "circle".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("circle".into()),
                assumptions: assumptions.clone(),
            },
            Operation::CheckedDigest {
                template: Source::Step(0),
                enzymes: vec!["NotI".into()],
                profile_revision: modified_duplex::profile_revision(),
                fragment_index: 0,
            },
        ],
    };
    let workspace = Workspace::from_documents(
        "Modification-preserving run",
        vec![source.clone(), target],
        1,
    )
    .unwrap();
    let committed = cloning_outputs::finalize(
        &workspace,
        &plan,
        binding(&workspace, 1),
        vec![InputEventBinding {
            input_id: "circle".into(),
            event: binding(&workspace, 0),
        }],
        &Options {
            name: "Named blocked circle".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let before = run_lineage::workspace_fingerprint(&committed.workspace).unwrap();
    let fresh = run_material::reconstruct(&committed.workspace, &selected(&committed, 1)).unwrap();
    let tracked = fresh.tracked();
    tracked.validate().unwrap();
    assert_eq!(fresh.reference().name, "Named blocked circle");
    assert_eq!(tracked.duplex.top.as_str(), "AAGCGGCCGCAA");
    assert_eq!(tracked.duplex.bottom.as_str(), "TTGCGGCCGCTT");
    assert_eq!(tracked.duplex.topology, Topology::Circular);
    assert!(tracked.duplex.nicks.is_empty());
    for strand in [Strand::Top, Strand::Bottom] {
        let positions = (0..12)
            .filter(|&position| {
                tracked.state_at(strand, Chemical::C5C, position).unwrap()
                    == Some(CallState::Modified)
            })
            .collect::<Vec<_>>();
        assert_eq!(positions, vec![3, 7]);
    }
    assert_eq!(modified_duplex::summary(tracked).unwrap().known_modified, 4);
    assert_eq!(tracked.generation.len(), 1);
    assert_eq!(tracked.generation[0].source_revision, source.revision());
    assert_eq!(tracked.generation[0].assumptions, assumptions);
    assert_eq!(fresh.report().status, Status::ConsistentUnderAssumptions);
    assert_eq!(
        run_lineage::workspace_fingerprint(&committed.workspace).unwrap(),
        before
    );
}

#[test]
fn wrong_run_step_stale_event_and_unrelated_equal_reference_cannot_supply_material() {
    let committed = pcr();
    let before = run_lineage::workspace_fingerprint(&committed.workspace).unwrap();
    for invalid in [
        {
            let mut value = selected(&committed, 0);
            value.run_id = "0".repeat(64);
            value
        },
        selected(&committed, usize::MAX),
        {
            let mut value = selected(&committed, 0);
            value.event.event_id = "stale-event".into();
            value
        },
        {
            let mut value = selected(&committed, 0);
            value.event.revision = "0".repeat(64);
            value
        },
    ] {
        assert!(run_material::reconstruct(&committed.workspace, &invalid).is_err());
    }
    assert_eq!(
        run_lineage::workspace_fingerprint(&committed.workspace).unwrap(),
        before
    );
    let mut workspace = committed.workspace.clone();
    let identical = workspace
        .head(&committed.receipt.construct.document_id)
        .unwrap()
        .clone();
    let unrelated = DocumentHistory::new("unrelated-identical-reference", identical).unwrap();
    let mut selection = selected(&committed, 0);
    selection.event = EventBinding {
        document_id: unrelated.id.clone(),
        event_id: unrelated.head_event().into(),
        revision: unrelated.head_revision().into(),
    };
    unrelated.contribute(&mut workspace).unwrap();
    assert_eq!(
        selection.event.revision,
        committed.receipt.construct.revision
    );
    let before = run_lineage::workspace_fingerprint(&workspace).unwrap();
    assert!(run_material::reconstruct(&workspace, &selection).is_err());
    assert!(
        run_material::originating_products(&workspace, &selection.event)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        run_lineage::workspace_fingerprint(&workspace).unwrap(),
        before
    );
}

#[test]
fn name_only_descendants_are_allowed_but_reverted_real_edits_break_material_ancestry() {
    let committed = pcr();
    let (renamed, selection) = changed_named_output(&committed, |history| {
        let mut name = history.molecule.clone();
        name.name = "Renamed again".into();
        history.commit(name, "Changed display name only").unwrap();
    });
    let before = run_lineage::workspace_fingerprint(&renamed).unwrap();
    let fresh = run_material::reconstruct(&renamed, &selection).unwrap();
    assert_eq!(fresh.reference().name, "Renamed again");
    assert_eq!(fresh.tracked().duplex.left, Terminus::PCR);
    assert_eq!(
        run_lineage::workspace_fingerprint(&renamed).unwrap(),
        before
    );
    for metadata in [false, true] {
        let (changed, selection) = changed_named_output(&committed, |history| {
            let baseline = history.molecule.clone();
            let altered = if metadata {
                let mut altered = baseline.clone();
                altered
                    .metadata
                    .insert("operator-note".into(), "Temporary changed assertion".into());
                altered
            } else {
                baseline.replace(0, 1, "T").unwrap()
            };
            history
                .commit(altered, "Real temporary sequence/metadata edit")
                .unwrap();
            history
                .commit(baseline, "Reverted exactly to the old reference")
                .unwrap();
        });
        assert_eq!(
            selection.event.revision,
            committed.receipt.construct.revision
        );
        assert_ne!(
            selection.event.event_id,
            committed.receipt.construct.event_id
        );
        let before = run_lineage::workspace_fingerprint(&changed).unwrap();
        assert!(run_material::reconstruct(&changed, &selection).is_err());
        assert!(
            run_material::originating_products(&changed, &selection.event)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            run_lineage::workspace_fingerprint(&changed).unwrap(),
            before
        );
    }
}

#[test]
fn documentary_origin_and_forged_preview_summary_do_not_promote_reference_only_output() {
    let source = Molecule::dna("Original reference", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent reference-only target",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Sequence-only overlap reference".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "reference".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("reference".into())],
            overlaps: vec![],
            circular: false,
        }],
    };
    let report = collagene_core::verification::verify(&plan, &target);
    assert_eq!(report.products.len(), 1);
    let workspace =
        Workspace::from_documents("Reference-only saved run", vec![source, target], 1).unwrap();
    let committed = cloning_outputs::finalize(
        &workspace,
        &plan,
        binding(&workspace, 1),
        vec![InputEventBinding {
            input_id: "reference".into(),
            event: binding(&workspace, 0),
        }],
        &Options {
            name: "Named reference-only result".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let selection = selected(&committed, 0);
    let before = run_lineage::workspace_fingerprint(&committed.workspace).unwrap();
    assert!(
        run_material::originating_products(&committed.workspace, &selection.event)
            .unwrap()
            .contains(&selection)
    );
    assert!(run_material::reconstruct(&committed.workspace, &selection).is_err());
    assert_eq!(
        run_lineage::workspace_fingerprint(&committed.workspace).unwrap(),
        before
    );

    // A physically plausible summary copied from a real independent PCR is
    // still only edited metadata on this reference-only descendant.
    let actual = pcr();
    let physical = run_material::reconstruct(&actual.workspace, &selected(&actual, 0)).unwrap();
    let forged = serde_json::to_string(&physical.tracked().duplex.summary().unwrap()).unwrap();
    let (changed, selection) = changed_named_output(&committed, |history| {
        let mut reference = history.molecule.clone();
        reference.metadata.insert("duplex_model".into(), forged);
        history
            .commit(reference, "Added a plausible physical preview summary")
            .unwrap();
    });
    let before = run_lineage::workspace_fingerprint(&changed).unwrap();
    assert!(
        run_material::originating_products(&changed, &selection.event)
            .unwrap()
            .is_empty()
    );
    assert!(run_material::reconstruct(&changed, &selection).is_err());
    assert_eq!(
        run_lineage::workspace_fingerprint(&changed).unwrap(),
        before
    );
}
