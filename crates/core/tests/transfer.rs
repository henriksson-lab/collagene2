use collagene_core::{transfer::*, verification::*, *};

fn fixture() -> (PlanningWorkspace, Candidate) {
    let left = Molecule::dna("left", "AAAACCCC", Topology::Linear).unwrap();
    let right = Molecule::dna("right", "CCCCGGGG", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AAAACCCCGGGG", Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Independent target".into(),
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
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![4],
            circular: false,
        }],
    };
    (
        PlanningWorkspace {
            target_revision: target.revision(),
            documents: vec![left, right, target],
        },
        Candidate::new(plan),
    )
}

#[test]
fn binary_exchange_preserves_exact_document_revisions_and_origin_annotations() {
    let molecule =
        collagene_core::formats::import_text(include_str!("../../../fixtures/origin-spanning.gb"))
            .unwrap()
            .molecules
            .remove(0);
    let workspace = PlanningWorkspace {
        target_revision: molecule.revision(),
        documents: vec![molecule.clone()],
    };
    let bytes = save_workspace(&workspace).unwrap();
    assert!(bytes.starts_with(b"CGD\0"));
    let restored = load_workspace(&bytes).unwrap();
    assert_eq!(restored.target_revision, molecule.revision());
    assert_eq!(restored.documents[0], molecule);
    assert_eq!(restored.documents[0].topology, Topology::Circular);
    let mut corrupt = bytes;
    let last = corrupt.len() - 1;
    corrupt[last] ^= 1;
    assert!(load_workspace(&corrupt).is_err());
}

#[test]
fn candidate_reconstructs_only_against_matching_current_documents() {
    let (workspace, candidate) = fixture();
    let bytes = save_candidate(&candidate).unwrap();
    assert!(bytes.starts_with(b"CGP\0"));
    let restored = load_candidate(&bytes).unwrap();
    let target = &workspace.documents[2];
    let plan = restored.bind(target, &workspace.documents).unwrap();
    assert_eq!(
        verify(&plan, target).status,
        Status::ConsistentUnderAssumptions
    );
    let mut edited = workspace.documents.clone();
    edited[0] = edited[0].replace(0, 1, "T").unwrap();
    assert!(
        restored
            .bind(target, &edited)
            .unwrap_err()
            .to_string()
            .contains("missing or has changed")
    );
    let mut renamed = workspace.documents.clone();
    renamed[0].name = "renamed".into();
    assert!(restored.bind(target, &renamed).is_err());
    let changed_target = target.replace(0, 1, "T").unwrap();
    assert!(
        restored
            .bind(&changed_target, &workspace.documents)
            .unwrap_err()
            .to_string()
            .contains("active target differs")
    );
    assert!(restored.bind(target, &workspace.documents[1..]).is_err());
}

#[test]
fn packages_reject_rebound_snapshots_and_do_not_promote_bad_science() {
    let (mut workspace, mut candidate) = fixture();
    workspace.target_revision = "missing".into();
    assert!(save_workspace(&workspace).is_err());
    workspace.target_revision = workspace.documents[2].revision();
    workspace.documents.push(workspace.documents[0].clone());
    assert!(save_workspace(&workspace).is_err());
    candidate.plan.name = "tampered".into();
    assert!(save_candidate(&candidate).is_err());
    candidate = Candidate::new(candidate.plan);
    candidate.plan.inputs[0].revision = "stale".into();
    candidate.plan_revision = plan_revision(&candidate.plan);
    assert!(save_candidate(&candidate).is_err());
    let (workspace, mut candidate) = fixture();
    candidate.plan.steps[0] = Operation::Overlap {
        fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
        overlaps: vec![3],
        circular: false,
    };
    let bad = Candidate::new(candidate.plan);
    // A candidate package records a proposal, not a scientific verdict.
    let imported = load_candidate(&save_candidate(&bad).unwrap()).unwrap();
    let bound = imported
        .bind(&workspace.documents[2], &workspace.documents)
        .unwrap();
    assert_eq!(
        verify(&bound, &workspace.documents[2]).status,
        Status::Invalid
    );
}

#[test]
fn invalid_step_indices_render_diagnostics_without_integer_overflow() {
    let (workspace, mut candidate) = fixture();
    candidate.plan.steps = vec![Operation::Overlap {
        fragments: vec![Source::Step(usize::MAX)],
        overlaps: vec![],
        circular: false,
    }];
    let target = &workspace.documents[2];
    assert_eq!(verify(&candidate.plan, target).status, Status::Invalid);
    let html = export_protocol_html(&candidate.plan, target);
    assert!(html.contains("<strong>Invalid</strong>"));
    assert!(html.contains("invalid step index"));
}
