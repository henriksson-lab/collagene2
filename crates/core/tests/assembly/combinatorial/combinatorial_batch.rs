use collagene_core::{
    chemistry::Terminus,
    cloning::assembly::recipe::Preparation,
    documents::{
        history::DocumentHistory,
        workspace::{Workspace, WorkspaceAttachment},
    },
    material::lineage::EventBinding,
    planning::combinatorial::{
        Bin, Choice, ChoiceSelection, ConstructOverrides, Design, Method, PreparationOverride,
        SelectedConstruct, TargetDeclaration,
        batch::{self, Request, RowRequest, RowState, Session},
    },
    verification::evidence::Status,
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
fn fixture() -> (Workspace, Request) {
    let w = Workspace::from_documents(
        "literal batch",
        vec![
            Molecule::dna("left", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("right", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("independent wrong target", "GGTT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let preparation = Some(Preparation::DeclaredDuplex {
        left: Terminus::CLEAVED,
        right: Terminus::CLEAVED,
    });
    let design = Design {
        material_bindings: vec![],
        id: "library".into(),
        name: "Literal two-part batch".into(),
        circular: false,
        method: Method::Ligation {
            allow_nicked: false,
        },
        bins: vec![
            Bin {
                id: "left-bin".into(),
                label: "Left".into(),
                choices: vec![Choice {
                    generation: None,
                    id: "left-choice".into(),
                    label: "Left choice".into(),
                    source: binding(&w, 0),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: 4,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation: preparation.clone(),
                    treatment: None,
                }],
            },
            Bin {
                id: "right-bin".into(),
                label: "Right".into(),
                choices: vec![Choice {
                    generation: None,
                    id: "right-choice".into(),
                    label: "Right choice".into(),
                    source: binding(&w, 1),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: 4,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation,
                    treatment: None,
                }],
            },
        ],
    };
    let row = SelectedConstruct {
        id: "good".into(),
        name: "Literal AACCTTAA".into(),
        choices: vec![
            ChoiceSelection {
                bin_id: "left-bin".into(),
                choice_id: "left-choice".into(),
            },
            ChoiceSelection {
                bin_id: "right-bin".into(),
                choice_id: "right-choice".into(),
            },
        ],
        target: TargetDeclaration::OriginalRegions {
            substitutions: vec![],
            edits: vec![],
            junction_edits: vec![],
        },
        overrides: None,
    };
    let mut wrong = row.clone();
    wrong.id = "wrong-target".into();
    wrong.name = "Keep the wrong independent target".into();
    wrong.target = TargetDeclaration::Existing(binding(&w, 2));
    let mut unfinished = row.clone();
    unfinished.id = "missing-preparation".into();
    unfinished.name = "Unfinished per-row preparation".into();
    unfinished.overrides = Some(ConstructOverrides {
        preparations: vec![PreparationOverride {
            preparation_reverse: false,
            generation: None,
            bin_id: "left-bin".into(),
            preparation: None,
            treatment: None,
        }],
        method: None,
        circular: None,
    });
    (
        w,
        Request {
            design: Some(design),
            rows: vec![
                RowRequest::Ready(row),
                RowRequest::Ready(wrong),
                RowRequest::Ready(unfinished),
                RowRequest::DraftError {
                    id: "invalid-edit".into(),
                    name: "Incomplete editable primer".into(),
                    choices: vec![],
                    message: "Finish the primer domain before preparing this row.".into(),
                },
            ],
        },
    )
}
fn run(workspace: Workspace, request: Request) -> Session {
    let mut session = Session::begin(workspace, request).unwrap();
    let context = session.summary().context_revision;
    let total = session.summary().total_rows;
    for index in 0..total {
        session.advance(&context, index).unwrap();
    }
    session
}
#[test]
fn exact_rows_continue_after_wrong_target_and_unprepared_edits() {
    let (w, r) = fixture();
    let before = serde_json::to_vec(&w).unwrap();
    let token = batch::begin_revision(&w, &r).unwrap();
    let mut session = Session::begin(w.clone(), r).unwrap();
    assert_eq!(session.summary().context_revision, token);
    assert!(session.finish_record().is_err());
    assert!(session.advance(&"0".repeat(64), 0).is_err());
    assert!(session.advance(&token, 1).is_err());
    let first = session.advance(&token, 0).unwrap();
    first.validate().unwrap();
    let first = first.row.unwrap();
    assert_eq!(first.status, Some(Status::ConsistentUnderAssumptions));
    assert_eq!(first.completed_steps, first.planned_steps);
    assert!(first.final_product_revision.is_some());
    let retry = session.advance(&token, 0).unwrap();
    assert_eq!(retry.next_index, 1);
    assert_eq!(
        serde_json::to_vec(&retry.row.unwrap()).unwrap(),
        serde_json::to_vec(&first).unwrap()
    );
    let wrong = session.advance(&token, 1).unwrap().row.unwrap();
    assert_eq!(wrong.status, Some(Status::Invalid));
    assert!(wrong.final_product_revision.is_some());
    assert_eq!(wrong.target_revision, Some(binding(&w, 2).revision));
    assert_eq!(
        session.advance(&token, 2).unwrap().row.unwrap().state,
        RowState::Error
    );
    assert_eq!(
        session.advance(&token, 3).unwrap().row.unwrap().state,
        RowState::DraftError
    );
    let record = session.finish_record().unwrap();
    record.validate().unwrap();
    assert_eq!(record.rows.len(), 4);
    assert_eq!(serde_json::to_vec(&w).unwrap(), before);
    let commit = session.capture_to(&w, &token).unwrap();
    commit
        .validate_binding(&batch::capture_revision(&w, &token).unwrap())
        .unwrap();
    assert!(commit.validate_binding(&"0".repeat(64)).is_err());
    assert_eq!(batch::list(&commit.workspace).unwrap()[0].evaluated_rows, 2);
}
#[test]
fn capture_keeps_later_heads_and_unrelated_native_rich_documents() {
    let (w, r) = fixture();
    let session = run(w.clone(), r);
    let token = session.summary().context_revision;
    let mut docs = DocumentHistory::restore_all(&w).unwrap();
    docs[0]
        .commit(
            Molecule::dna("later edited left", "TCCC", Topology::Linear).unwrap(),
            "Later edit",
        )
        .unwrap();
    let mut current = w.clone();
    current.documents.clear();
    current.history.clear();
    for doc in docs {
        doc.contribute(&mut current).unwrap();
    }
    let mut rich = Molecule::dna("unrelated rich reference", "ACGT", Topology::Linear).unwrap();
    rich.features = (0..10_000)
        .map(|index| Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, 1),
            qualifiers: vec![("label".into(), Some(format!("feature {index}")))],
            needs_review: false,
        })
        .collect();
    let rich = DocumentHistory::new("unrelated-native", rich).unwrap();
    rich.contribute(&mut current).unwrap();
    current.active_document = "unrelated-native".into();
    current.validate().unwrap();
    let heads = serde_json::to_vec(&current.documents).unwrap();
    let history = serde_json::to_vec(&current.history).unwrap();
    let active = current.active_document.clone();
    let commit = session.capture_to(&current, &token).unwrap();
    assert_eq!(
        serde_json::to_vec(&commit.workspace.documents).unwrap(),
        heads
    );
    assert_eq!(
        serde_json::to_vec(&commit.workspace.history).unwrap(),
        history
    );
    assert_eq!(commit.workspace.active_document, active);
    let replay = batch::replay(&commit.workspace, &commit.attachment_hash, 0).unwrap();
    assert_eq!(
        replay.current.status,
        Some(Status::ConsistentUnderAssumptions)
    );
    assert_eq!(replay.current.input_bindings[0].event, binding(&w, 0));
    assert_eq!(
        commit.workspace.documents[0].head,
        current.documents[0].head
    );
}
#[test]
fn native_records_deduplicate_and_archived_assertions_need_fresh_replay() {
    let (w, r) = fixture();
    let session = run(w.clone(), r.clone());
    let token = session.summary().context_revision;
    let saved = session.into_capture(&token).unwrap();
    let retry = run(saved.workspace.clone(), r);
    let retry_token = retry.summary().context_revision;
    let repeated = retry.into_capture(&retry_token).unwrap();
    assert_eq!(repeated.attachment_hash, saved.attachment_hash);
    assert_eq!(
        repeated.workspace.attachments.len(),
        saved.workspace.attachments.len()
    );
    let mut record = batch::load(&saved.workspace, &saved.attachment_hash).unwrap();
    record.engine_version = "archived engine assertion".into();
    record.rows[1].status = Some(Status::ConsistentUnderAssumptions);
    record.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            "combinatorial-batch-record-2",
            &record.engine_version,
            &record.request,
            &record.rows,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let bytes = serde_json::to_vec(&record).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut archive = saved.workspace.clone();
    archive.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: batch::KIND.into(),
            name: "Archived assertion".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let opened = batch::load(&archive, &hash).unwrap();
    assert_eq!(opened.engine_version, "archived engine assertion");
    assert_eq!(
        opened.rows[1].status,
        Some(Status::ConsistentUnderAssumptions)
    );
    let fresh = batch::replay(&archive, &hash, 1).unwrap();
    assert_eq!(fresh.current.status, Some(Status::Invalid));
    assert_eq!(fresh.current_engine, batch::ENGINE_VERSION);
}
#[test]
fn identical_choices_keep_per_construct_topology_and_method_overrides() {
    let (w, mut request) = fixture();
    let RowRequest::Ready(base) = request.rows[0].clone() else {
        unreachable!()
    };
    let mut circle = base.clone();
    circle.id = "circle".into();
    circle.name = "Circular variant".into();
    circle.overrides = Some(ConstructOverrides {
        preparations: vec![],
        method: None,
        circular: Some(true),
    });
    let mut overlap = base.clone();
    overlap.id = "wrong-overlap".into();
    overlap.name = "Explicit incompatible overlap".into();
    overlap.overrides = Some(ConstructOverrides {
        preparations: vec![],
        method: Some(Method::Overlap {
            finishing: collagene_core::cloning::assembly::recipe::OverlapFinishing::ReferenceOnly,
            overlaps: vec![1],
            options: Default::default(),
        }),
        circular: None,
    });
    request.rows = vec![
        RowRequest::Ready(base),
        RowRequest::Ready(circle),
        RowRequest::Ready(overlap),
    ];
    let record = run(w, request).finish_record().unwrap();
    assert_eq!(
        record.rows[0].status,
        Some(Status::ConsistentUnderAssumptions)
    );
    assert_eq!(
        record.rows[1].status,
        Some(Status::ConsistentUnderAssumptions)
    );
    assert_eq!(record.rows[2].status, Some(Status::Invalid));
    assert_eq!(record.rows[0].combination_id, record.rows[1].combination_id);
    assert_ne!(
        record.rows[0].declaration_revision,
        record.rows[1].declaration_revision
    );
    assert_ne!(
        record.rows[0].input_bindings[0].input_id,
        record.rows[1].input_bindings[0].input_id
    );
}
#[test]
fn unavailable_event_owner_is_a_local_error_and_not_a_saved_source_attestation() {
    let (w, mut r) = fixture();
    let mut bad = r.design.as_ref().unwrap().bins[0].choices[0].clone();
    bad.id = "wrong-owner-choice".into();
    bad.source.document_id = w.documents[1].id.clone();
    r.design.as_mut().unwrap().bins[0].choices.push(bad);
    let RowRequest::Ready(mut row) = r.rows[0].clone() else {
        unreachable!()
    };
    row.id = "wrong-owner".into();
    row.choices[0].choice_id = "wrong-owner-choice".into();
    r.rows = vec![r.rows[0].clone(), RowRequest::Ready(row)];
    let session = run(w, r);
    let token = session.summary().context_revision;
    let commit = session.into_capture(&token).unwrap();
    let record = batch::load(&commit.workspace, &commit.attachment_hash).unwrap();
    assert_eq!(record.rows[0].state, RowState::Evaluated);
    assert_eq!(record.rows[1].state, RowState::Error);
    assert!(record.rows[1].input_bindings.is_empty());
}
#[test]
fn all_draft_errors_and_aggregate_budget_exhaustion_remain_explicit() {
    let (w, r) = fixture();
    let all = Request {
        design: None,
        rows: vec![r.rows[3].clone()],
    };
    let session = run(w.clone(), all);
    let token = session.summary().context_revision;
    let saved = session.into_capture(&token).unwrap();
    assert_eq!(
        batch::load(&saved.workspace, &saved.attachment_hash)
            .unwrap()
            .rows[0]
            .state,
        RowState::DraftError
    );
    let RowRequest::Ready(base) = r.rows[0].clone() else {
        unreachable!()
    };
    let mut rows = vec![];
    for index in 0..42 {
        let mut row = base.clone();
        row.id = format!("budget-row-{index}");
        row.overrides = Some(ConstructOverrides {
            preparations: vec![],
            method: Some(Method::Overlap {
                finishing:
                    collagene_core::cloning::assembly::recipe::OverlapFinishing::ReferenceOnly,
                overlaps: vec![1],
                options: collagene_core::cloning::overlap::profiles::ProfileOptions {
                    max_comparisons: 50_000_000,
                    ..Default::default()
                },
            }),
            circular: None,
        });
        rows.push(RowRequest::Ready(row));
    }
    let record = run(
        w,
        Request {
            design: r.design,
            rows,
        },
    )
    .finish_record()
    .unwrap();
    assert_eq!(
        record
            .rows
            .iter()
            .filter(|r| r.state == RowState::Evaluated)
            .count(),
        39
    );
    assert!(record.rows[39..].iter().all(|r| r.state == RowState::Error
        && r.status.is_none()
        && r.error.as_ref().unwrap().contains("Not computed")));
}
