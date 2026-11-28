//! Literal source/target/event fixtures for private atomic output release.
use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::{
    chemistry::{EndTreatment, Terminus},
    cloning::assembly::recipe::Preparation,
    documents::{
        history::DocumentHistory,
        transfer::Candidate,
        workspace::{Workspace, WorkspaceAttachment},
    },
    material::{
        lineage::{self as run_lineage, EventBinding},
        outputs::Options,
    },
    planning::{
        combinatorial::{
            Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
            batch,
            release::{self, ReleaseSession, Request, Selection},
        },
        plan::{DesignPlan, Operation, PlanInput, Source},
    },
    sequence::analysis::Primer,
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
fn fixture(wrong_second: bool, pcr: bool, rich: bool) -> (Workspace, Request) {
    let (left, right) = if pcr {
        ("TTGCAATCGTCGATGCTAGC", "ATCGTCGATGCTAGCGGTAC")
    } else {
        ("AACC", "TTAA")
    };
    let mut models = vec![
        Molecule::dna("left", left, Topology::Linear).unwrap(),
        Molecule::dna("right", right, Topology::Linear).unwrap(),
        Molecule::dna("wrong independent target", "GGTT", Topology::Linear).unwrap(),
    ];
    if rich {
        let mut m = Molecule::dna("unrelated rich reference", "ACGT", Topology::Linear).unwrap();
        m.features = (0..10_000)
            .map(|i| Feature {
                kind: "misc_feature".into(),
                location: Location::range(0, 1),
                qualifiers: vec![("label".into(), Some(format!("feature {i}")))],
                needs_review: false,
            })
            .collect();
        models.push(m);
    }
    let mut workspace = Workspace::from_documents("release fixture", models, 0).unwrap();
    let prepare = |i| {
        if pcr {
            let (f, r) = if i == 0 {
                ("TTGCAATC", "GCTAGCAT")
            } else {
                ("ATCGTCGA", "GTACCGCT")
            };
            Preparation::ExactPcr {
                forward: Primer::new("F", f, "").unwrap(),
                reverse: Primer::new("R", r, "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            }
        } else {
            Preparation::DeclaredDuplex {
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED,
            }
        }
    };
    let source_design = Design {
        material_bindings: vec![],
        id: "release-library".into(),
        name: "Release literal inputs".into(),
        circular: false,
        method: Method::Ligation {
            allow_nicked: false,
        },
        bins: (0..2)
            .map(|i| Bin {
                id: format!("bin-{i}"),
                label: format!("Bin {i}"),
                choices: vec![Choice {
                    anneal_bottom: None,
                    generation: None,
                    id: format!("choice-{i}"),
                    label: format!("Choice {i}"),
                    source: binding(&workspace, i),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: if pcr { 20 } else { 4 },
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation: Some(prepare(i)),
                    treatment: pcr.then_some(EndTreatment::Phosphorylate),
                }],
            })
            .collect(),
    };
    let selections = (0..2)
        .map(|i| ChoiceSelection {
            bin_id: format!("bin-{i}"),
            choice_id: format!("choice-{i}"),
        })
        .collect::<Vec<_>>();
    let rows = (0..2)
        .map(|i| {
            batch::RowRequest::Ready(SelectedConstruct {
                standard: None,
                id: format!("construct-{i}"),
                name: format!("Declared construct {i}"),
                choices: selections.clone(),
                target: if i == 1 && wrong_second {
                    TargetDeclaration::Existing(binding(&workspace, 2))
                } else {
                    TargetDeclaration::OriginalRegions {
                        substitutions: vec![],
                        edits: vec![],
                        junction_edits: vec![],
                    }
                },
                overrides: None,
            })
        })
        .collect::<Vec<_>>();
    let request = batch::Request {
        design: Some(source_design),
        rows,
    };
    // Historical uncomputed summaries are deliberately sufficient input. Release
    // must execute fresh native finalization regardless of any saved verdict.
    let evidence = request
        .rows
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let batch::RowRequest::Ready(row) = row else {
                unreachable!()
            };
            batch::RowEvidence {
                index,
                construct_id: row.id.clone(),
                name: row.name.clone(),
                state: batch::RowState::Error,
                status: None,
                error: Some("Historical row not computed.".into()),
                combination_id: None,
                declaration_revision: None,
                input_bindings: vec![],
                target_binding: None,
                plan_revision: None,
                target_revision: None,
                final_product_revision: None,
                completed_steps: 0,
                planned_steps: 0,
                preparation_complete: false,
                tailing_profile: None,
                findings: vec![],
                assumptions: vec![],
                exhaustive_alternative_search: None,
            }
        })
        .collect::<Vec<_>>();
    let engine = "archived uncomputed request".to_string();
    let context = blake3::hash(
        &serde_json::to_vec(&("combinatorial-batch-record-6", &engine, &request, &evidence))
            .unwrap(),
    )
    .to_hex()
    .to_string();
    let record = batch::Record {
        schema_version: 6,
        context_revision: context,
        engine_version: engine,
        request,
        rows: evidence,
    };
    record.validate().unwrap();
    let bytes = serde_json::to_vec(&record).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    workspace.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: batch::KIND.into(),
            name: "Historical batch request".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    if rich {
        let molecule = workspace.head(&workspace.documents[3].id).unwrap().clone();
        let plan = DesignPlan {
            standard: None,
            name: "Unrelated supported rich candidate".into(),
            target_revision: molecule.revision(),
            inputs: vec![PlanInput {
                id: "unrelated".into(),
                revision: molecule.revision(),
                molecule,
            }],
            steps: vec![Operation::Overlap {
                fragments: vec![Source::Input("unrelated".into())],
                overlaps: vec![],
                circular: false,
            }],
        };
        workspace.candidates.push(Candidate::new(plan));
    }
    workspace.validate().unwrap();
    let release = Request {
        batch_hash: hash,
        rows: (0..2)
            .map(|i| Selection {
                construct_id: format!("construct-{i}"),
                options: Options {
                    name: format!("Released construct {i}"),
                    show_prepared: false,
                    save_primers: pcr,
                    allow_incomplete: false,
                },
            })
            .collect(),
    };
    (workspace, release)
}
#[test]
fn all_selected_rows_release_atomically_with_source_run_and_primer_links() {
    let (workspace, request) = fixture(false, true, false);
    let source = serde_json::to_vec(&workspace).unwrap();
    let fp = run_lineage::workspace_fingerprint(&workspace).unwrap();
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    let inspected = session.inspect(0).unwrap();
    inspected.validate().unwrap();
    assert_eq!(
        inspected.product.sequence.as_str(),
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC"
    );
    assert_eq!(
        inspected.target.sequence.as_str(),
        inspected.product.sequence.as_str()
    );
    assert_ne!(inspected.product.name, inspected.target.name);
    let mut tampered = inspected.clone();
    tampered.product.name = "not the retained event".into();
    assert!(tampered.validate().is_err());
    assert!(session.inspect(1).is_err());
    assert!(
        session
            .summary()
            .row
            .unwrap()
            .output
            .primer_attachment
            .is_some()
    );
    assert_eq!(serde_json::to_vec(&workspace).unwrap(), source);
    session.advance(&context, 1).unwrap();
    let commit = session.finish_to(&workspace, &context).unwrap();
    commit.validate_binding(&context, &fp).unwrap();
    let loaded = release::load(&commit.workspace, &commit.manifest_hash).unwrap();
    loaded
        .validate_current_declarations(&commit.workspace)
        .unwrap();
    let mut documentary = loaded.clone();
    documentary.rows[0].declaration_revision = "0".repeat(64);
    documentary.validate(&commit.workspace).unwrap();
    assert!(
        documentary
            .validate_current_declarations(&commit.workspace)
            .is_err()
    );
    assert_eq!(release::list(&commit.workspace).unwrap().len(), 1);
    let asset = &commit.workspace.attachments[&commit.manifest_hash];
    let manifest: release::Manifest = serde_json::from_slice(&asset.bytes).unwrap();
    assert_eq!(manifest.rows.len(), 2);
    for row in manifest.rows {
        let construct = &commit.workspace.revisions[&row.output.construct.revision];
        assert_eq!(
            construct.sequence.as_str(),
            "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC"
        );
        let hash = row.output.primer_attachment.unwrap();
        let table: collagene_core::formats::oligos::OligoTable =
            serde_json::from_slice(&commit.workspace.attachments[&hash].bytes).unwrap();
        assert_eq!(table.records.len(), 4);
        let event = commit
            .workspace
            .history
            .iter()
            .find(|event| event.id == row.target.event_id)
            .unwrap();
        assert_eq!(
            event
                .parents
                .iter()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            row.inputs
                .iter()
                .map(|input| input.event.event_id.clone())
                .collect()
        );
    }
}

#[test]
fn paired_oligo_release_keeps_primary_then_bottom_roles_in_fresh_runs_and_current_manifests() {
    use collagene_core::{chemistry::AnnealOptions, cloning::assembly::recipe};
    let (mut workspace, mut request) = fixture(false, false, false);
    let mut record = batch::load(&workspace, &request.batch_hash).unwrap();
    let mut options = AnnealOptions::default();
    options.min_paired_bases = 4;
    options.top_chemistry = Terminus::CLEAVED;
    options.bottom_chemistry = Terminus::CLEAVED;
    let batch::RowRequest::Ready(first_row) = &record.request.rows[0] else {
        unreachable!()
    };
    let construct_id = first_row.id.clone();
    let design = record.request.design.as_mut().unwrap();
    for (index, bin) in design.bins.iter_mut().enumerate() {
        let choice = &mut bin.choices[0];
        let primary = collagene_core::planning::combinatorial::source_role_ids(
            &design.id,
            &construct_id,
            &bin.id,
            &choice.id,
        )
        .unwrap()
        .0;
        choice.anneal_bottom = Some(binding(&workspace, if index == 0 { 2 } else { 1 }));
        choice.preparation = Some(Preparation::AnnealedOligos {
            bottom_input_id: recipe::anneal_bottom_input_id(&primary).unwrap(),
            options: options.clone(),
            candidate_index: 0,
        });
    }
    workspace.attachments.remove(&request.batch_hash).unwrap();
    let mut evaluated = batch::Session::begin(workspace.clone(), record.request).unwrap();
    let progress = evaluated.summary();
    for index in 0..progress.total_rows {
        evaluated
            .advance(&progress.context_revision, index)
            .unwrap();
    }
    let captured = evaluated.into_capture(&progress.context_revision).unwrap();
    workspace = captured.workspace;
    request.batch_hash = captured.attachment_hash;
    let before = serde_json::to_value(&workspace).unwrap();
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    assert_eq!(
        session.inspect(0).unwrap().product.sequence.as_str(),
        "AACCTTAA"
    );
    session.advance(&context, 1).unwrap();
    let released = session.finish_to(&workspace, &context).unwrap();
    let manifest = release::load(&released.workspace, &released.manifest_hash).unwrap();
    manifest
        .validate_current_declarations(&released.workspace)
        .unwrap();
    for row in &manifest.rows {
        assert_eq!(row.inputs.len(), 4);
        assert_eq!(
            row.inputs
                .iter()
                .map(|i| i.event.clone())
                .collect::<Vec<_>>(),
            vec![
                binding(&workspace, 0),
                binding(&workspace, 1),
                binding(&workspace, 2),
                binding(&workspace, 1)
            ]
        );
        assert_ne!(row.inputs[1].input_id, row.inputs[3].input_id);
        let run = run_lineage::load_run(&released.workspace, &row.output.run.run_id).unwrap();
        assert_eq!(
            run.inputs
                .iter()
                .map(|i| (&i.input_id, &i.event))
                .collect::<Vec<_>>(),
            row.inputs
                .iter()
                .map(|i| (&i.input_id, &i.event))
                .collect::<Vec<_>>()
        );
        let target = released
            .workspace
            .history
            .iter()
            .find(|e| e.id == row.target.event_id)
            .unwrap();
        assert_eq!(
            target.parents,
            vec![
                binding(&workspace, 0).event_id,
                binding(&workspace, 1).event_id,
                binding(&workspace, 2).event_id
            ]
        );
    }
    let mut obsolete = manifest.clone();
    obsolete.schema_version = 1;
    assert!(obsolete.validate(&released.workspace).is_err());
    let mut unsupported = released.workspace.clone();
    unsupported
        .attachments
        .get_mut(&released.manifest_hash)
        .unwrap()
        .kind = "combinatorial-release-v1".into();
    assert!(release::list(&unsupported).is_err());
    assert!(release::load(&unsupported, &released.manifest_hash).is_err());
    // Wire capacity is 24 roles, independent of the 12 contribution limit.
    // These synthetic rows are decode checks, never fresh release authority.
    let mut encoded = serde_json::to_value(&manifest.rows[0]).unwrap();
    let inputs = (0..24)
        .map(|index| {
            let mut input = encoded["inputs"][index % 4].clone();
            input["input_id"] = serde_json::json!(format!("role-{index}"));
            input
        })
        .collect::<Vec<_>>();
    encoded["inputs"] = serde_json::json!(inputs);
    assert_eq!(
        serde_json::from_value::<release::ReleasedRow>(encoded.clone())
            .unwrap()
            .inputs
            .len(),
        24
    );
    let extra = encoded["inputs"][0].clone();
    encoded["inputs"].as_array_mut().unwrap().push(extra);
    assert!(serde_json::from_value::<release::ReleasedRow>(encoded).is_err());
    assert_eq!(serde_json::to_value(&workspace).unwrap(), before);
}
#[test]
fn a_late_wrong_target_failure_discards_the_entire_successful_prefix() {
    let (workspace, request) = fixture(true, false, false);
    let before = serde_json::to_vec(&workspace).unwrap();
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    assert!(session.advance(&context, 1).is_err());
    assert!(session.summary().failed.is_some());
    assert!(session.inspect(0).is_err());
    assert!(session.into_commit(&context).is_err());
    assert_eq!(serde_json::to_vec(&workspace).unwrap(), before);
    assert_eq!(
        workspace
            .attachments
            .values()
            .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
            .count(),
        0
    );
    // A failure after scientific reconstruction, while making requested outputs,
    // also discards the prefix: these duplex-declaration plans contain no primers.
    let (workspace, mut request) = fixture(false, false, false);
    request.rows[1].options.save_primers = true;
    let original = serde_json::to_vec(&workspace).unwrap();
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    assert!(session.advance(&context, 1).is_err());
    assert!(session.into_commit(&context).is_err());
    assert_eq!(serde_json::to_vec(&workspace).unwrap(), original);
}
#[test]
fn changed_live_workspace_prevents_installation_and_safe_archive_names_never_use_raw_ids() {
    let (workspace, request) = fixture(false, false, false);
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    session.advance(&context, 1).unwrap();
    let mut docs = DocumentHistory::restore_all(&workspace).unwrap();
    let edit_end = docs[0].molecule.sequence.len();
    docs[0]
        .commit(
            vec![
                Edit::ReplaceSequence {
                    start: 0,
                    end: edit_end,
                    replacement: "TCCC".into(),
                },
                Edit::Rename {
                    name: "later edit".into(),
                },
                Edit::SetTopology {
                    topology: Topology::Linear,
                },
            ],
            "Later source edit",
        )
        .unwrap();
    let mut current = workspace.clone();
    current.documents.clear();
    current.history.clear();
    for doc in docs {
        doc.contribute(&mut current).unwrap();
    }
    assert!(session.finish_to(&current, &context).is_err());
    let path = release::archive_entry(0, "..", "gb").unwrap();
    assert!(path.starts_with("constructs/0001-"));
    assert_eq!(path.split('/').count(), 2);
    assert!(!path.contains(".."));
    assert!(release::archive_entry(0, "a", "../../secret").is_err());
    let (_, request) = fixture(false, false, false);
    let mut too_many = serde_json::to_value(request).unwrap();
    let row = too_many["rows"][0].clone();
    too_many["rows"] = serde_json::Value::Array(vec![row; 257]);
    assert!(serde_json::from_value::<release::Request>(too_many).is_err());
}
#[test]
fn unrelated_native_rich_documents_and_candidates_do_not_block_small_cloning_release() {
    // Requires staged native-scope patch. Relevant plan/target preflight remains unchanged.
    let (workspace, request) = fixture(false, false, true);
    let rich = workspace.documents[3].clone();
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    session.advance(&context, 1).unwrap();
    let commit = session.into_commit(&context).unwrap();
    assert!(
        commit
            .workspace
            .documents
            .iter()
            .any(|d| d.id == rich.id && d.head == rich.head)
    );
    assert_eq!(
        commit.workspace.revisions[&rich.head].features.len(),
        10_000
    );
    assert_eq!(
        commit.workspace.candidates.len(),
        workspace.candidates.len()
    );
}

#[test]
fn actual_selected_rich_source_still_receives_scientific_bounds() {
    let (mut workspace, mut request) = fixture(false, false, true);
    let mut record = batch::load(&workspace, &request.batch_hash).unwrap();
    record.request.design.as_mut().unwrap().bins[0].choices[0].source = binding(&workspace, 3);
    record.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            "combinatorial-batch-record-6",
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
    workspace.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: batch::KIND.into(),
            name: "Relevant rich source".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    request.batch_hash = hash;
    let before = serde_json::to_vec(&workspace).unwrap();
    let mut session = ReleaseSession::begin(workspace.clone(), request).unwrap();
    let context = session.summary().context_revision;
    assert!(session.advance(&context, 0).is_err());
    assert!(session.into_commit(&context).is_err());
    assert_eq!(serde_json::to_vec(&workspace).unwrap(), before);
}
