use collagene_core::documents::{
    history::transition::Transition, reference_edit::Operation as Edit,
};
use collagene_core::{
    chemistry::Terminus,
    documents::workspace::{self, RevisionEvent, Workspace, WorkspaceAttachment},
    material::{
        derivation as material,
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Molecule, Topology};
const DNA: &str = "AGTCCGATGCTAGAATTCGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn pin(b: &EventBinding) -> Selection {
    Selection {
        document_id: b.document_id.clone(),
        revision: b.revision.clone(),
        head_event: Some(b.event_id.clone()),
    }
}
fn fixture() -> (Workspace, EventBinding, String, String, String) {
    let stock = Molecule::dna("Original PCR source", DNA, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent full PCR target", DNA, Topology::Linear).unwrap();
    let private = Molecule::dna(
        "PRIVATE UNRELATED DOCUMENT",
        "CCCCAAAAGGGGTTTT",
        Topology::Circular,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Local unpublished workspace",
        vec![stock.clone(), target.clone(), private.clone()],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Fresh unique PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let mut run = run_lineage::record_run(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        false,
    )
    .unwrap();
    let product = &run.receipt.products[0];
    let source = EventBinding {
        document_id: product.document_id.clone(),
        event_id: product.event_id.clone(),
        revision: product.revision.clone(),
    };
    let material = run_material::Selection {
        run_id: run.receipt.run_id.clone(),
        step_index: 0,
        event: source.clone(),
    };
    // Installation is an explicit operation on the active exact saved PCR
    // document, rather than on the independent target left active by capture.
    run.workspace.active_document = source.document_id.clone();
    let origin = collagene_core::material::capsule::capsule(&run.workspace, &material).unwrap();
    let c = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let eco = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection = collagene_core::cloning::restriction::catalog::select(&c, &[eco]).unwrap();
    let installed = material::prepare_install(
        &run.workspace,
        &material::Request {
            sources: vec![material::BoundSource {
                role: material::SourceRole::Template,
                event: source,
                origin: material::SourceOrigin::Saved(run_material::MaterialOrigin::Run(origin)),
            }],
            derivation: material::Derivation::CatalogDigest {
                conditions: None,
                selection,
                options:
                    collagene_core::cloning::restriction::inspection::InspectionOptions::default(),
                fragment_index: 0,
            },
            name: "Selected PCR EcoRI fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let product = installed.receipt.product.clone();
    let record = installed.receipt.record_hash.clone();
    let mut w = installed.workspace;
    // Later private source edits cannot be chosen instead of the captured source
    // event merely because the owning document now points somewhere else.
    let old = event(&w, 0);
    let previous_source = w.revisions.get(&old.revision).unwrap();
    let actions = vec![
        Edit::ReplaceSequence {
            start: 0,
            end: previous_source.sequence.len(),
            replacement: "TTTTGGGGCCCCAAAA".into(),
        },
        Edit::Rename {
            name: "PRIVATE NEW SOURCE HEAD".into(),
        },
        Edit::SetTopology {
            topology: Topology::Linear,
        },
    ];
    let edited = Edit::apply_all(previous_source, &actions).unwrap();
    let revision = edited.revision();
    let e = RevisionEvent::new(
        &old.document_id,
        &revision,
        vec![old.event_id.clone()],
        "Private edit after digest",
        Transition::ReferenceEdit {
            source: old.clone(),
            actions,
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    w.revisions.insert(revision.clone(), edited);
    w.history.push(e.clone());
    w.documents[0].undo.push(old.revision.clone());
    w.documents[0].undo_events.push(old.event_id);
    w.documents[0].head = revision.clone();
    w.documents[0].head_event = Some(e.id);
    let bytes = b"private presentation state".to_vec();
    w.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        WorkspaceAttachment {
            kind: "oligo-table-v1".into(),
            name: "Private UI".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    w.validate().unwrap();
    (w, product, record, private.revision(), revision)
}
#[test]
fn actual_pcr_digest_fragment_public_and_private_closure_preserve_fresh_material_without_private_current_heads()
 {
    let (w, product, record, private, new_head) = fixture();
    let before = serde_json::to_value(&w).unwrap();
    let fresh = material::reconstruct(&w, &product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), "AGTCCGATGCTAG");
    assert_eq!(fresh.tracked().duplex.top.as_str(), "AGTCCGATGCTAG");
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), "AATTCTAGCATCGGACT");
    assert_eq!(fresh.tracked().duplex.left, Terminus::PCR);
    assert_eq!(fresh.tracked().duplex.right, Terminus::CLEAVED);
    for projected in [
        sharing::project_document_closure(&w, &pin(&product)).unwrap(),
        sharing::project_private_document_closure(&w, &pin(&product)).unwrap(),
    ] {
        assert!(projected.inventory.attachments.contains(&record));
        assert!(!projected.workspace.revisions.contains_key(&private));
        assert!(!projected.workspace.revisions.contains_key(&new_head));
        assert!(projected.workspace.candidates.is_empty());
        assert!(
            projected
                .workspace
                .attachments
                .values()
                .all(|a| a.kind != "oligo-table-v1")
        );
        assert!(
            projected
                .workspace
                .documents
                .iter()
                .all(|d| d.undo.is_empty() && d.redo.is_empty())
        );
        let reopened =
            workspace::load_workspace(&workspace::save_workspace(&projected.workspace).unwrap())
                .unwrap();
        let again = material::reconstruct(&reopened, &product).unwrap();
        assert_eq!(again.reference(), fresh.reference());
        assert_eq!(again.tracked(), fresh.tracked());
        assert_eq!(again.incomplete(), fresh.incomplete());
        let capsule = material::capsule(&reopened, &product).unwrap();
        let replay = material::replay(&capsule, again.reference()).unwrap();
        assert_eq!(replay.tracked(), fresh.tracked());
    }
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
}

#[test]
fn named_digest_setup_closure_retains_exact_run_and_derived_origins_without_private_setups() {
    use collagene_core::cloning::restriction::{catalog, inspection::InspectionOptions, setups};
    use collagene_core::material::{capsule, derivation::SourceOrigin};
    let (workspace, derived, _, private_revision, _) = fixture();
    let run_record = workspace
        .attachments
        .values()
        .find(|asset| asset.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
        .map(|asset| run_lineage::load_run_record(&asset.bytes).unwrap())
        .unwrap();
    let run_product = &run_record.products[0];
    let run_event = EventBinding {
        document_id: run_product.document_id.clone(),
        event_id: run_product.event_id.clone(),
        revision: run_product.revision.clone(),
    };
    let selections = [
        run_material::MaterialSelection::Run(run_material::Selection {
            run_id: run_record.run_id.clone(),
            step_index: run_product.step_index,
            event: run_event,
        }),
        run_material::MaterialSelection::Derived(derived),
    ];
    let embedded = catalog::builtin_catalog().unwrap();
    let enzymes = catalog::select(&embedded, &["ecori".into()]).unwrap();
    for selection in selections {
        let mut origins =
            capsule::documentary_origins(&workspace, &[selection.clone()], &[]).unwrap();
        let recipe = setups::Recipe {
            source: selection.event().clone(),
            origin: SourceOrigin::Saved(origins.materials.remove(0)),
            selection: enzymes.clone(),
            options: InspectionOptions::default(),
            cleavage_mode: setups::CleavageMode::AssumedCompleteMethylationUnassessed,
            conditions: None,
        };
        let source = &workspace.revisions[&recipe.source.revision];
        let inspection = recipe.inspect(source).unwrap();
        let (mut staged, hash) = setups::capture(
            &embedded,
            &workspace,
            &recipe,
            "Reusable saved digest",
            &inspection,
        )
        .unwrap();
        let private_document = staged
            .documents
            .iter()
            .find(|doc| doc.head == private_revision)
            .unwrap();
        let private_recipe = setups::Recipe {
            source: EventBinding {
                document_id: private_document.id.clone(),
                event_id: private_document.head_event.clone().unwrap(),
                revision: private_revision.clone(),
            },
            origin: SourceOrigin::OriginalReference,
            selection: enzymes.clone(),
            options: InspectionOptions::default(),
            cleavage_mode: setups::CleavageMode::AssumedCompleteMethylationUnassessed,
            conditions: None,
        };
        let private_report = private_recipe
            .inspect(&staged.revisions[&private_revision])
            .unwrap();
        let (with_private, private_setup) = setups::capture(
            &embedded,
            &staged,
            &private_recipe,
            "Private unrelated digest",
            &private_report,
        )
        .unwrap();
        staged = with_private;
        let before = serde_json::to_value(&staged).unwrap();
        for projection in [
            sharing::project_document_closure(&staged, &pin(&recipe.source)).unwrap(),
            sharing::project_private_document_closure(&staged, &pin(&recipe.source)).unwrap(),
        ] {
            assert!(projection.inventory.attachments.contains(&hash));
            assert!(!projection.inventory.attachments.contains(&private_setup));
            assert!(
                !projection
                    .workspace
                    .revisions
                    .contains_key(&private_revision)
            );
            let restored = workspace::load_workspace(
                &workspace::save_workspace(&projection.workspace).unwrap(),
            )
            .unwrap();
            let record = setups::load(&restored, &hash).unwrap();
            record.recipe.validate_for_workspace(&restored).unwrap();
            let fresh = record
                .recipe
                .inspect(&restored.revisions[&record.recipe.source.revision])
                .unwrap();
            assert_eq!(fresh.context_revision, inspection.context_revision);
            assert_eq!(
                serde_json::to_value(&restored.attachments[&hash]).unwrap(),
                serde_json::to_value(&staged.attachments[&hash]).unwrap()
            );
        }
        assert_eq!(serde_json::to_value(&staged).unwrap(), before);
    }
}

#[test]
fn obsolete_named_digest_setup_kind_cannot_enter_document_publication() {
    let (mut workspace, product, _, _, _) = fixture();
    let bytes = b"{}".to_vec();
    workspace.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        WorkspaceAttachment {
            kind: "named-digest-setup-v1".into(),
            name: "Obsolete setup".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    assert!(sharing::inspect_document_closure(&workspace, &pin(&product)).is_err());
    assert!(sharing::inspect_private_document_closure(&workspace, &pin(&product)).is_err());
    assert!(workspace::save_workspace(&workspace).is_err());
}
#[test]
fn unrelated_selected_document_never_publishes_derived_record_or_its_private_origin_runs() {
    let (w, _, record, private, _) = fixture();
    let document = w.documents.iter().find(|d| d.head == private).unwrap();
    let selected = Selection {
        document_id: document.id.clone(),
        revision: document.head.clone(),
        head_event: document.head_event.clone(),
    };
    let out = sharing::project_document_closure(&w, &selected).unwrap();
    assert_eq!(out.workspace.documents.len(), 1);
    assert_eq!(out.workspace.revisions.len(), 1);
    assert_eq!(
        out.workspace.revisions[&private].name,
        "PRIVATE UNRELATED DOCUMENT"
    );
    assert!(!out.inventory.attachments.contains(&record));
    assert!(out.workspace.attachments.is_empty());
}
#[test]
fn orphan_reserved_producer_corrupt_record_and_wrong_product_binding_refuse_instead_of_opaque_fallback()
 {
    let (w, product, hash, _, _) = fixture();
    let mut orphan = w.clone();
    orphan.attachments.remove(&hash);
    assert!(sharing::project_document_closure(&orphan, &pin(&product)).is_err());
    let mut corrupt = w.clone();
    let mut a = corrupt.attachments.remove(&hash).unwrap();
    a.bytes = b"{}".to_vec();
    corrupt
        .attachments
        .insert(blake3::hash(&a.bytes).to_hex().to_string(), a);
    assert!(sharing::project_document_closure(&corrupt, &pin(&product)).is_err());
    let mut wrong = w.clone();
    let mut a = wrong.attachments.remove(&hash).unwrap();
    let mut record = material::load_record(&a.bytes).unwrap();
    record.product.document_id = wrong.documents[2].id.clone();
    a.bytes = material::save_record(&record).unwrap();
    wrong
        .attachments
        .insert(blake3::hash(&a.bytes).to_hex().to_string(), a);
    assert!(sharing::project_private_document_closure(&wrong, &pin(&product)).is_err());
}

#[test]
fn later_pcr_run_closes_nested_saved_digest_recipe_and_replays_after_publication() {
    let (mut w, product, record, private, new_head) = fixture();
    let origin = material::capsule(&w, &product).unwrap();
    let source = w.revisions[&product.revision].clone();
    let target = Molecule::dna(
        "Independent downstream PCR target",
        "AGTCCGATGCTAG",
        Topology::Linear,
    )
    .unwrap();
    let id = "independent-downstream-target".to_string();
    let revision = target.revision();
    let target_event = RevisionEvent::new(
        &id,
        &revision,
        vec![],
        "Independent target import",
        Transition::IndependentReference,
    )
    .unwrap();
    w.documents.push(workspace::WorkspaceDocument {
        id: id.clone(),
        head: revision.clone(),
        head_event: Some(target_event.id.clone()),
        undo: vec![],
        redo: vec![],
        undo_events: vec![],
        redo_events: vec![],
    });
    w.revisions.insert(revision.clone(), target.clone());
    w.history.push(target_event.clone());
    // The baseline PCR engine requires non-overlapping annealing sites.
    // Literal six-base sites 0..6 and 7..13 leave one intervening source base.
    let forward = Primer::new("F", "AGTCCG", "").unwrap();
    let reverse = Primer::new("R", "CTAGCA", "").unwrap();
    let geometries =
        collagene_core::sequence::analysis::simulate_pcr_geometry(&source, &forward, &reverse)
            .unwrap();
    assert_eq!(geometries.len(), 1);
    assert_eq!(geometries[0].product.sequence.as_str(), "AGTCCGATGCTAG");
    assert_eq!(
        (
            geometries[0].product.forward_start,
            geometries[0].product.reverse_start
        ),
        (0, 7)
    );
    let plan = DesignPlan {
        name: "PCR from saved digest material".into(),
        target_revision: revision.clone(),
        inputs: vec![PlanInput {
            id: "digest".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![
            Operation::ReconstructMaterial {
                template: Source::Input("digest".into()),
                origin,
            },
            Operation::Pcr {
                template: Source::Step(0),
                forward,
                reverse,
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        ],
    };
    let committed = run_lineage::record_run(
        &w,
        &plan,
        EventBinding {
            document_id: id,
            event_id: target_event.id,
            revision,
        },
        vec![InputEventBinding {
            input_id: "digest".into(),
            event: product,
        }],
        false,
    )
    .unwrap();
    let p = committed.receipt.products.last().unwrap();
    let derived = EventBinding {
        document_id: p.document_id.clone(),
        event_id: p.event_id.clone(),
        revision: p.revision.clone(),
    };
    let selected = run_material::Selection {
        run_id: committed.receipt.run_id,
        step_index: 1,
        event: derived.clone(),
    };
    let before = serde_json::to_vec(&committed.workspace).unwrap();
    for export in [
        sharing::project_document_closure(&committed.workspace, &pin(&derived)).unwrap(),
        sharing::project_private_document_closure(&committed.workspace, &pin(&derived)).unwrap(),
    ] {
        assert!(export.inventory.attachments.contains(&record));
        assert!(!export.workspace.revisions.contains_key(&private));
        assert!(!export.workspace.revisions.contains_key(&new_head));
        let archive = workspace::save_workspace(&export.workspace).unwrap();
        let reopened = workspace::load_workspace(&archive).unwrap();
        let fresh = run_material::reconstruct(&reopened, &selected).unwrap();
        assert_eq!(fresh.tracked().duplex.top.as_str(), "AGTCCGATGCTAG");
    }
    assert_eq!(serde_json::to_vec(&committed.workspace).unwrap(), before);
}

#[test]
fn newest_editable_drafts_are_presentation_only_not_publication_dependencies() {
    let (mut w, product, _, _, _) = fixture();
    let bytes = br#"{"version":14,"private_unfinished_text":"DO NOT PUBLISH"}"#.to_vec();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_KIND.into(),
            name: "private draft".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    for projected in [
        sharing::project_document_closure(&w, &pin(&product)).unwrap(),
        sharing::project_private_document_closure(&w, &pin(&product)).unwrap(),
    ] {
        assert!(!projected.inventory.attachments.contains(&hash));
        assert!(!projected.workspace.attachments.contains_key(&hash));
    }
}
