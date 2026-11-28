//! Retained digest recipes are historical evidence, never a freshness receipt.
use collagene_core::{
    cloning::restriction::{
        catalog::{self as catalog_digest, catalog},
        conditions::{self, ChosenFormulation, EvidenceStatus, ReactionConditions, Unassessed},
        history as restriction_history,
        inspection::InspectionOptions,
        setups::{self as digest_setups, CleavageMode, Record},
    },
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace, WorkspaceAttachment},
    },
    material::{
        derivation::{self, SourceOrigin, digest::DigestPreview},
        lineage::EventBinding,
    },
};
use collagene_model::{Molecule, Topology};

fn binding(frame: &Workspace, index: usize) -> EventBinding {
    let doc = &frame.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn recipe(source: &EventBinding, report: &DigestPreview) -> digest_setups::Recipe {
    digest_setups::Recipe {
        source: source.clone(),
        origin: SourceOrigin::OriginalReference,
        selection: report.selection.clone(),
        options: report.inspection.options.clone(),
        cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
        conditions: None,
    }
}
fn fixture() -> (catalog::Catalog, Workspace, EventBinding, DigestPreview) {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let molecule = Molecule::dna(
        "original template",
        "TTTGAATTCAAAGGATCCAAA",
        Topology::Circular,
    )
    .unwrap();
    let frame =
        Workspace::from_documents("saved digests", vec![molecule.clone(), molecule], 0).unwrap();
    let source = binding(&frame, 0);
    let selection = catalog_digest::select(&catalog, &["ecori".into(), "bamhi".into()]).unwrap();
    let report = derivation::digest::inspect_digest(
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &SourceOrigin::OriginalReference,
        &selection,
        &InspectionOptions::default(),
    )
    .unwrap();
    (catalog, frame, source, report)
}
fn reaction() -> ReactionConditions {
    ReactionConditions {
        schema_version: 1,
        formulations: vec![
            ChosenFormulation {
                enzyme_id: "ecori".into(),
                formulation_id: Some("neb-r3101".into()),
            },
            ChosenFormulation {
                enzyme_id: "bamhi".into(),
                formulation_id: None,
            },
        ],
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37_000,
        incubation_minutes: 60,
        minimum_activity_percent: 75,
        activity_context: Some("NEB unit-assay relative activity".into()),
        temperature_context: Some("vendor recommended incubation".into()),
    }
}
fn insert(frame: &Workspace, value: &serde_json::Value) -> (Workspace, String) {
    let bytes = serde_json::to_vec(value).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut staged = frame.clone();
    staged.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: digest_setups::KIND.into(),
            name: value["name"].as_str().unwrap_or("invalid").into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    (staged, hash)
}
fn advance(frame: &Workspace) -> Workspace {
    let mut histories = DocumentHistory::restore_all(frame).unwrap();
    histories[0]
        .commit(
            Molecule::dna(
                "later template",
                "CCCGGATCCAAAGAATTCTTT",
                Topology::Circular,
            )
            .unwrap(),
            "later edit",
        )
        .unwrap();
    let mut next = frame.clone();
    next.documents.clear();
    next.history.clear();
    for history in histories {
        history.contribute(&mut next).unwrap();
    }
    next.validate().unwrap();
    next
}

#[test]
fn capture_roundtrips_native_archive_and_changes_only_one_attachment() {
    let (catalog, frame, source, report) = fixture();
    let original = serde_json::to_value(&frame).unwrap();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Diagnostic double digest",
        &report,
    )
    .unwrap();
    assert_eq!(serde_json::to_value(&frame).unwrap(), original);
    let mut without = saved.clone();
    without.attachments.remove(&hash);
    assert_eq!(serde_json::to_value(without).unwrap(), original);
    let record = digest_setups::load(&saved, &hash).unwrap();
    assert_eq!(record.name, "Diagnostic double digest");
    assert_eq!(record.recipe.source, source);
    assert_eq!(record.recipe.selection, report.selection);
    assert_eq!(record.recipe.options, report.inspection.options);
    assert_eq!(
        record.recipe.cleavage_mode,
        CleavageMode::AssumedCompleteMethylationUnassessed
    );
    assert_eq!(
        serde_json::to_value(&record.report).unwrap(),
        serde_json::to_value(&report).unwrap()
    );
    let listed = digest_setups::list(&saved).unwrap();
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].attachment_hash, hash);
    assert!(listed[0].scan_complete && listed[0].products_determined);
    let archive = workspace::save_workspace(&saved).unwrap();
    let restored = workspace::load_workspace(&archive).unwrap();
    assert_eq!(
        serde_json::to_value(digest_setups::load(&restored, &hash).unwrap()).unwrap(),
        serde_json::to_value(record).unwrap()
    );
    let (deduplicated, again) = digest_setups::capture(
        &catalog,
        &saved,
        &recipe(&source, &report),
        "Diagnostic double digest",
        &report,
    )
    .unwrap();
    assert_eq!(again, hash);
    assert_eq!(deduplicated.attachments.len(), saved.attachments.len());
}

#[test]
fn historical_setup_survives_head_advance_and_does_not_alias_equal_dna_events() {
    let (catalog, frame, source, report) = fixture();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Original digest",
        &report,
    )
    .unwrap();
    let later = advance(&saved);
    assert_ne!(later.documents[0].head, source.revision);
    assert_eq!(
        digest_setups::load(&later, &hash).unwrap().recipe.source,
        source
    );
    assert_eq!(
        digest_setups::list(&later).unwrap()[0].source_name,
        "original template"
    );
    let original = serde_json::to_value(&later).unwrap();
    assert!(
        digest_setups::capture(
            &catalog,
            &later,
            &recipe(&binding(&later, 0), &report),
            "Wrong head",
            &report
        )
        .is_err()
    );
    let alias = binding(&later, 1);
    assert_eq!(alias.revision, source.revision);
    assert!(
        digest_setups::capture(
            &catalog,
            &later,
            &recipe(&alias, &report),
            "Wrong event",
            &report
        )
        .is_err()
    );
    assert_eq!(serde_json::to_value(&later).unwrap(), original);
    let alias_report = derivation::digest::inspect_digest(
        &restriction_history::source_role(&alias).unwrap(),
        &later.revisions[&alias.revision],
        &SourceOrigin::OriginalReference,
        &report.selection,
        &report.inspection.options,
    )
    .unwrap();
    let (both, alias_hash) = digest_setups::capture(
        &catalog,
        &later,
        &recipe(&alias, &alias_report),
        "Original digest",
        &alias_report,
    )
    .unwrap();
    assert_ne!(alias_hash, hash);
    assert_eq!(digest_setups::list(&both).unwrap().len(), 2);
    assert_eq!(both.documents[0].head, later.documents[0].head);
}

#[test]
fn gel_recipe_retains_actual_ordered_geometry_catalog_and_limits() {
    let (catalog, frame, source, report) = fixture();
    let mut options = report.inspection.options.clone();
    options.max_comparisons = 12345;
    let selected = catalog_digest::select(&catalog, &["bamhi".into(), "ecori".into()]).unwrap();
    let report = derivation::digest::inspect_digest(
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &SourceOrigin::OriginalReference,
        &selected,
        &options,
    )
    .unwrap();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Gel lane recipe",
        &report,
    )
    .unwrap();
    let recipe = digest_setups::gel_recipe(&saved, &hash).unwrap();
    assert_eq!(recipe.source, source);
    assert_eq!(recipe.selection, selected);
    assert_eq!(recipe.options, options);
    assert_eq!(recipe.selection.enzymes[0].id, "bamhi");
    assert_eq!(recipe.selection.enzymes[0].enzyme.top_cut, 1);
    assert_eq!(recipe.selection.enzymes[0].enzyme.bottom_cut, 5);
    let json = serde_json::to_value(digest_setups::load(&saved, &hash).unwrap()).unwrap();
    assert!(json.get("fresh").is_none());
    assert!(json.get("verified").is_none());
    assert!(json["recipe"].get("methylation_unmodified").is_none());
}

#[test]
fn partial_scan_can_be_saved_without_inventing_determined_products() {
    let (catalog, frame, source, prior) = fixture();
    let options = InspectionOptions {
        max_comparisons: 1,
        ..Default::default()
    };
    let report = derivation::digest::inspect_digest(
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &SourceOrigin::OriginalReference,
        &prior.selection,
        &options,
    )
    .unwrap();
    assert!(!report.inspection.scan_complete && !report.inspection.products_determined);
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Incomplete exploration",
        &report,
    )
    .unwrap();
    let loaded = digest_setups::load(&saved, &hash).unwrap();
    assert!(!loaded.report.inspection.scan_complete);
    assert!(!loaded.report.inspection.products_determined);
    assert!(loaded.report.inspection.fragments.is_empty());
    assert!(!digest_setups::list(&saved).unwrap()[0].products_determined);
}

#[test]
fn historical_catalog_pins_remain_readable_but_cannot_rebind_to_another_catalog() {
    let (catalog, frame, source, report) = fixture();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Old catalog recipe",
        &report,
    )
    .unwrap();
    let mut edition = catalog.edition().clone();
    edition.edition = "future catalog".into();
    let future = catalog::Catalog::from_edition(catalog::seal(edition).unwrap()).unwrap();
    assert!(digest_setups::load(&saved, &hash).is_ok());
    assert!(
        digest_setups::capture(
            &future,
            &saved,
            &recipe(&source, &report),
            "Changed catalog",
            &report
        )
        .is_err()
    );
    assert!(
        digest_setups::gel_recipe(&saved, &hash)
            .unwrap()
            .selection
            .validate_for(&future)
            .is_err()
    );
}

#[test]
fn external_catalog_documentary_setup_loads_without_installed_membership_but_cannot_compute() {
    let (catalog, frame, source, report) = fixture();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Documentary external pin",
        &report,
    )
    .unwrap();
    let mut record = digest_setups::load(&saved, &hash).unwrap();
    let mut edition = catalog.edition().clone();
    edition.edition = "independent unavailable catalog edition".into();
    let unavailable = catalog::Catalog::from_edition(catalog::seal(edition).unwrap()).unwrap();
    record.recipe.selection =
        catalog_digest::select(&unavailable, &["ecori".into(), "bamhi".into()]).unwrap();
    record.report.selection = record.recipe.selection.clone();
    record.recipe.conditions = Some(reaction());
    record.condition_report = Some(
        conditions::assess(
            &unavailable,
            &record.recipe.selection,
            record.recipe.conditions.as_ref().unwrap(),
        )
        .unwrap(),
    );
    // Documentary declaration identity is independent of installed-catalog lookup.
    record.report.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            derivation::ENGINE_VERSION,
            "digest-inspection",
            collagene_core::verification::ENGINE_VERSION,
            restriction_history::source_role(&source).unwrap().as_str(),
            &frame.revisions[&source.revision],
            &record.recipe.origin,
            &record.recipe.selection,
            &record.recipe.options,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let (documentary, external_hash) = insert(&frame, &serde_json::to_value(&record).unwrap());
    let reopened = digest_setups::load(&documentary, &external_hash).unwrap();
    assert_eq!(
        reopened.recipe.selection.catalog_edition,
        "independent unavailable catalog edition"
    );
    let condition_report = reopened.condition_report.as_ref().unwrap();
    condition_report
        .validate_binding(
            &conditions::documentary_revision(
                &reopened.recipe.selection,
                reopened.recipe.conditions.as_ref().unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    assert!(
        condition_report
            .validate_for(
                &catalog,
                &reopened.recipe.selection,
                reopened.recipe.conditions.as_ref().unwrap()
            )
            .is_err()
    );
    assert!(
        reopened
            .recipe
            .inspection_revision(&frame.revisions[&source.revision])
            .is_err()
    );
    assert!(
        reopened
            .recipe
            .inspect(&frame.revisions[&source.revision])
            .is_err()
    );
    let restored =
        workspace::load_workspace(&workspace::save_workspace(&documentary).unwrap()).unwrap();
    assert!(digest_setups::load(&restored, &external_hash).is_ok());
}

#[test]
fn reaction_setup_retains_typed_documentary_evidence_without_changing_ideal_digest() {
    let (catalog, frame, source, report) = fixture();
    let mut chosen = recipe(&source, &report);
    let ideal_revision = chosen
        .inspection_revision(&frame.revisions[&source.revision])
        .unwrap();
    chosen.conditions = Some(reaction());
    assert_eq!(
        chosen
            .inspection_revision(&frame.revisions[&source.revision])
            .unwrap(),
        ideal_revision
    );
    assert_eq!(
        serde_json::to_value(chosen.inspect(&frame.revisions[&source.revision]).unwrap()).unwrap(),
        serde_json::to_value(&report).unwrap()
    );
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &chosen, "Recorded reaction", &report).unwrap();
    let record = digest_setups::load(&saved, &hash).unwrap();
    assert_eq!(record.schema_version, 3);
    assert_eq!(record.recipe.source, source);
    assert_eq!(record.recipe.conditions, chosen.conditions);
    let evidence = record.condition_report.as_ref().unwrap();
    evidence
        .validate_for(
            &catalog,
            &chosen.selection,
            chosen.conditions.as_ref().unwrap(),
        )
        .unwrap();
    assert_eq!(evidence.documentary_status, EvidenceStatus::Unknown);
    assert_eq!(
        evidence.formulations[0].activity.status,
        EvidenceStatus::Known
    );
    assert_eq!(
        evidence.formulations[1].activity.status,
        EvidenceStatus::Unknown
    );
    assert!(evidence.formulations[1].formulation.is_none());
    assert_eq!(evidence.formulations[0].star_risk, Unassessed::Unknown);
    assert_eq!(
        evidence.formulations[0].incubation_duration,
        Unassessed::Unknown
    );
    let restored = workspace::load_workspace(&workspace::save_workspace(&saved).unwrap()).unwrap();
    let reopened = digest_setups::load(&restored, &hash).unwrap();
    assert_eq!(reopened.condition_report, record.condition_report);
    assert_eq!(reopened.recipe.source, source);
    let original_preview = serde_json::to_value(&record.report).unwrap();
    chosen.conditions.as_mut().unwrap().incubation_minutes = 30;
    let (second, other) =
        digest_setups::capture(&catalog, &saved, &chosen, "Recorded reaction", &report).unwrap();
    assert_ne!(hash, other);
    let changed = digest_setups::load(&second, &other).unwrap();
    assert_eq!(
        serde_json::to_value(changed.report).unwrap(),
        original_preview
    );
    assert_ne!(
        changed.condition_report.as_ref().unwrap().request_revision,
        evidence.request_revision
    );
    assert_eq!(second.documents.len(), frame.documents.len());
    assert_eq!(second.history.len(), frame.history.len());
    assert_eq!(second.revisions, frame.revisions);
}

#[test]
fn reaction_attachment_refuses_unpaired_options_stale_bindings_tampering_and_v2() {
    let (catalog, frame, source, report) = fixture();
    let mut chosen = recipe(&source, &report);
    chosen.conditions = Some(reaction());
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &chosen, "Reaction", &report).unwrap();
    let original = serde_json::to_value(digest_setups::load(&saved, &hash).unwrap()).unwrap();
    let mut cases = Vec::new();
    let mut v = original.clone();
    v["recipe"]["conditions"]["incubation_minutes"] = 30.into();
    cases.push(v);
    let mut v = original.clone();
    v["condition_report"] = serde_json::Value::Null;
    cases.push(v);
    let mut v = original.clone();
    v["recipe"]["conditions"] = serde_json::Value::Null;
    cases.push(v);
    let mut v = original.clone();
    v["condition_report"]["formulations"][0]["activity"]["status"] = "Unknown".into();
    cases.push(v);
    let mut v = original.clone();
    v["condition_report"]["sources"][0]["publisher"] = "forged".into();
    cases.push(v);
    let mut v = original.clone();
    v["schema_version"] = 2.into();
    cases.push(v);
    for value in cases {
        let (invalid, hash) = insert(&frame, &value);
        assert!(digest_setups::load(&invalid, &hash).is_err());
        assert!(workspace::save_workspace(&invalid).is_err());
    }
    let mut obsolete = saved.clone();
    obsolete.attachments.get_mut(&hash).unwrap().kind = "named-digest-setup-v2".into();
    assert!(digest_setups::load(&obsolete, &hash).is_err());
    assert!(workspace::save_workspace(&obsolete).is_err());
    let mut malformed = chosen;
    malformed.conditions.as_mut().unwrap().formulations[0].formulation_id =
        Some("neb-r0141".into());
    assert!(
        digest_setups::capture(&catalog, &frame, &malformed, "Wrong formulation", &report).is_err()
    );
}

#[test]
fn tampered_source_options_geometry_version_mode_and_evidence_flags_fail() {
    let (catalog, frame, source, report) = fixture();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Recorded digest",
        &report,
    )
    .unwrap();
    let original = serde_json::to_value(digest_setups::load(&saved, &hash).unwrap()).unwrap();
    let mut cases = Vec::new();
    let mut value = original.clone();
    value["schema_version"] = serde_json::json!(1);
    cases.push(value);
    let mut value = original.clone();
    value["fresh"] = serde_json::json!(true);
    cases.push(value);
    let mut value = original.clone();
    value["recipe"]["options"]["max_sites"] = serde_json::json!(1);
    cases.push(value);
    let mut value = original.clone();
    value["recipe"]["selection"]["enzymes"][0]["enzyme"]["top_cut"] = serde_json::json!(2);
    cases.push(value);
    let mut value = original.clone();
    value["recipe"]["source"] = serde_json::to_value(binding(&frame, 1)).unwrap();
    cases.push(value);
    let mut value = original.clone();
    value["recipe"]["cleavage_mode"] = serde_json::json!("MethylationConfirmedUnmodified");
    cases.push(value);
    let mut value = original.clone();
    value["report"]["inspection"]["sites"][0]["spans"][0]["end"] = serde_json::json!(1000);
    cases.push(value);
    for value in cases {
        let (corrupted, hash) = insert(&frame, &value);
        let before = serde_json::to_value(&corrupted).unwrap();
        assert!(digest_setups::load(&corrupted, &hash).is_err());
        assert!(digest_setups::list(&corrupted).is_err());
        assert!(digest_setups::delete(&corrupted, &hash).is_err());
        assert_eq!(serde_json::to_value(&corrupted).unwrap(), before);
    }
    let mut wrong_type = saved.clone();
    wrong_type.attachments.get_mut(&hash).unwrap().kind = "unrelated".into();
    assert!(digest_setups::load(&wrong_type, &hash).is_err());
    let mut wrong_name = saved.clone();
    wrong_name.attachments.get_mut(&hash).unwrap().name = "renamed outside record".into();
    assert!(digest_setups::load(&wrong_name, &hash).is_err());
}

#[test]
fn delete_only_removes_valid_setup_and_keeps_original_document_history() {
    let (catalog, frame, source, report) = fixture();
    let (one, hash) =
        digest_setups::capture(&catalog, &frame, &recipe(&source, &report), "One", &report)
            .unwrap();
    let (two, other) =
        digest_setups::capture(&catalog, &one, &recipe(&source, &report), "Two", &report).unwrap();
    let deleted = digest_setups::delete(&two, &hash).unwrap();
    assert!(digest_setups::load(&deleted, &hash).is_err());
    assert!(digest_setups::load(&deleted, &other).is_ok());
    assert_eq!(
        serde_json::to_value(&deleted.documents).unwrap(),
        serde_json::to_value(&frame.documents).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&deleted.history).unwrap(),
        serde_json::to_value(&frame.history).unwrap()
    );
    assert_eq!(deleted.revisions, frame.revisions);
    let empty = digest_setups::delete(&deleted, &other).unwrap();
    assert_eq!(
        serde_json::to_value(empty).unwrap(),
        serde_json::to_value(frame).unwrap()
    );
}

#[test]
fn workspace_attachment_capacity_is_atomic_and_idempotent_capture_still_works() {
    let (catalog, frame, source, report) = fixture();
    let (mut full, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Existing",
        &report,
    )
    .unwrap();
    for index in 1..workspace::MAX_ATTACHMENTS {
        let bytes = format!("unrelated {index}").into_bytes();
        let id = blake3::hash(&bytes).to_hex().to_string();
        full.attachments.insert(
            id,
            WorkspaceAttachment {
                kind: "unrelated".into(),
                name: format!("Unrelated {index}"),
                media_type: "text/plain".into(),
                bytes,
            },
        );
    }
    full.validate().unwrap();
    let before = serde_json::to_value(&full).unwrap();
    assert!(
        digest_setups::capture(
            &catalog,
            &full,
            &recipe(&source, &report),
            "Would exceed native capacity",
            &report
        )
        .is_err()
    );
    assert_eq!(serde_json::to_value(&full).unwrap(), before);
    let (retry, repeated) = digest_setups::capture(
        &catalog,
        &full,
        &recipe(&source, &report),
        "Existing",
        &report,
    )
    .unwrap();
    assert_eq!(repeated, hash);
    assert_eq!(retry.attachments.len(), workspace::MAX_ATTACHMENTS);
    for name in ["", " ", "ends with space ", "bad\nname"] {
        assert!(
            digest_setups::capture(&catalog, &frame, &recipe(&source, &report), name, &report)
                .is_err()
        );
    }
    assert!(
        digest_setups::capture(
            &catalog,
            &frame,
            &recipe(&source, &report),
            &"x".repeat(257),
            &report
        )
        .is_err()
    );
    let _: Record = digest_setups::load(&retry, &hash).unwrap();
}

#[test]
fn unrelated_large_annotation_document_does_not_restrict_saved_digest_operations() {
    let (catalog, prior, source, report) = fixture();
    let original = prior.revisions[&source.revision].clone();
    let mut unrelated =
        Molecule::dna("unrelated annotated reference", "ACGT", Topology::Linear).unwrap();
    unrelated.features = vec![
        collagene_model::Feature {
            kind: "misc_feature".into(),
            location: collagene_model::Location::range(0, 4),
            qualifiers: Vec::new(),
            needs_review: false,
        };
        10_000
    ];
    let frame =
        Workspace::from_documents("native annotation capacity", vec![original, unrelated], 0)
            .unwrap();
    frame.validate().unwrap();
    let source = binding(&frame, 0);
    let report = derivation::digest::inspect_digest(
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &SourceOrigin::OriginalReference,
        &report.selection,
        &report.inspection.options,
    )
    .unwrap();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &recipe(&source, &report),
        "Unrelated features preserved",
        &report,
    )
    .unwrap();
    assert_eq!(digest_setups::list(&saved).unwrap().len(), 1);
    assert_eq!(
        digest_setups::load(&saved, &hash).unwrap().recipe.source,
        source
    );
    let restored = workspace::load_workspace(&workspace::save_workspace(&saved).unwrap()).unwrap();
    assert_eq!(
        restored.revisions[&restored.documents[1].head]
            .features
            .len(),
        10_000
    );
    assert_eq!(digest_setups::list(&restored).unwrap().len(), 1);
    assert!(digest_setups::load(&restored, &hash).is_ok());
    assert!(
        digest_setups::delete(&restored, &hash)
            .unwrap()
            .attachments
            .is_empty()
    );
}

fn saved_material_recipe(derived: bool) -> (Workspace, digest_setups::Recipe) {
    use collagene_core::{
        chemistry::Strand,
        material::{
            capsule,
            lineage::InputEventBinding,
            outputs,
            run::{MaterialOrigin, Selection},
        },
        planning::plan::{DesignPlan, Operation, PlanInput, Source},
        sequence::analysis::Primer,
    };
    let source = Molecule::dna("Original six bases", "AAACCC", Topology::Linear).unwrap();
    let forward = Primer::new("Forward", "AAA", "CCGAATTC").unwrap();
    let reverse = Primer::new("Reverse", "GGG", "GGCTCGAG").unwrap();
    let (frame, event, origin) = if derived {
        let frame =
            Workspace::from_documents("Target-free material", vec![source.clone()], 0).unwrap();
        let previews = collagene_core::verification::pcr::simulate_pcr_documents(
            &source,
            &forward,
            &reverse,
            Strand::Top,
        )
        .unwrap();
        assert_eq!(previews.len(), 1);
        let request = derivation::Request {
            sources: vec![derivation::BoundSource {
                role: derivation::SourceRole::Template,
                event: binding(&frame, 0),
                origin: SourceOrigin::OriginalReference,
            }],
            derivation: derivation::Derivation::ExactPcr {
                forward,
                reverse,
                template_strand: Strand::Top,
                selection: derivation::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: previews[0].molecule.revision(),
                    geometry: previews[0].geometry.clone(),
                },
            },
            name: "Derived digest substrate".into(),
            save_primers: false,
        };
        let prepared = derivation::prepare_install(&frame, &request).unwrap();
        let event = prepared.receipt.product;
        let origin =
            MaterialOrigin::Derived(derivation::capsule(&prepared.workspace, &event).unwrap());
        (prepared.workspace, event, origin)
    } else {
        let target = Molecule::dna(
            "Independent target",
            "CCGAATTCAAACCCCTCGAGCC",
            Topology::Linear,
        )
        .unwrap();
        let plan = DesignPlan {
            name: "Source-backed PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "template".into(),
                revision: source.revision(),
                molecule: source.clone(),
            }],
            steps: vec![Operation::Pcr {
                template: Source::Input("template".into()),
                forward,
                reverse,
                template_strand: Strand::Top,
            }],
        };
        let frame = Workspace::from_documents("Run substrate", vec![source, target], 1).unwrap();
        let prepared = outputs::finalize(
            &frame,
            &plan,
            binding(&frame, 1),
            vec![InputEventBinding {
                input_id: "template".into(),
                event: binding(&frame, 0),
            }],
            &outputs::Options {
                name: "Saved digest substrate".into(),
                show_prepared: false,
                save_primers: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let event = prepared.receipt.construct;
        let origin = MaterialOrigin::Run(
            capsule::capsule(
                &prepared.workspace,
                &Selection {
                    run_id: prepared.receipt.run.run_id,
                    step_index: 0,
                    event: event.clone(),
                },
            )
            .unwrap(),
        );
        (prepared.workspace, event, origin)
    };
    assert_eq!(
        frame.revisions[&event.revision].sequence.as_str(),
        "CCGAATTCAAACCCCTCGAGCC"
    );
    let selected = catalog_digest::select(
        &catalog_digest::builtin_catalog().unwrap(),
        &["ecori".into()],
    )
    .unwrap();
    (
        frame,
        digest_setups::Recipe {
            source: event,
            origin: SourceOrigin::Saved(origin),
            selection: selected,
            options: InspectionOptions::default(),
            cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
            conditions: None,
        },
    )
}

#[test]
fn saved_run_and_target_free_material_setup_restart_gel_uses_literal_real_strand_lengths() {
    use collagene_core::output::gel::{
        self, Calibration, CalibrationPoint, LaneInput, Preparation,
    };
    for derived in [false, true] {
        let (frame, recipe) = saved_material_recipe(derived);
        let source = &frame.revisions[&recipe.source.revision];
        assert!(collagene_core::chemistry::Duplex::from_molecule(source).is_err());
        let report = recipe.inspect(source).unwrap();
        assert_eq!(report.inspection.fragments.len(), 2);
        assert_eq!(
            report
                .inspection
                .fragments
                .iter()
                .map(|f| (f.duplex.top_length, f.duplex.bottom_length))
                .collect::<Vec<_>>(),
            vec![(3, 7), (19, 15)]
        );
        let (saved, hash) = digest_setups::capture(
            &catalog_digest::builtin_catalog().unwrap(),
            &frame,
            &recipe,
            "Saved source diagnostic",
            &report,
        )
        .unwrap();
        let restarted =
            workspace::load_workspace(&workspace::save_workspace(&saved).unwrap()).unwrap();
        let retained = digest_setups::load(&restarted, &hash).unwrap();
        assert_eq!(
            serde_json::to_value(&retained.recipe.origin).unwrap(),
            serde_json::to_value(&recipe.origin).unwrap()
        );
        let lane = LaneInput {
            label: "Fresh saved material".into(),
            molecule: restarted.revisions[&recipe.source.revision].clone(),
            preparation: Preparation::CatalogDigest {
                recipe: digest_setups::gel_recipe(&restarted, &hash).unwrap(),
            },
            mass_ng: Some(22.),
        };
        let gel = gel::predict_gel(
            &[lane],
            &Calibration {
                label: "Literal ladder".into(),
                resolution: 0.01,
                points: vec![
                    CalibrationPoint {
                        size_bp: 1.,
                        migration: 0.9,
                    },
                    CalibrationPoint {
                        size_bp: 1000.,
                        migration: 0.1,
                    },
                ],
            },
        )
        .unwrap();
        assert!(gel.lanes[0].error.is_none(), "{:?}", gel.lanes[0].error);
        assert_eq!(
            gel.lanes[0]
                .fragments
                .iter()
                .map(|f| f.size_bp)
                .collect::<Vec<_>>(),
            vec![5., 17.]
        );
        let mut reference_only = retained.recipe.clone();
        reference_only.origin = SourceOrigin::OriginalReference;
        let recognition_only = reference_only
            .inspect(&restarted.revisions[&recipe.source.revision])
            .unwrap();
        assert!(recognition_only.source_material.is_none());
        assert!(recognition_only.incomplete);
        assert!(!recognition_only.inspection.products_determined);
        assert!(recognition_only.inspection.fragments.is_empty());
        assert!(recognition_only.fragment_modifications.is_empty());
        assert!(
            digest_setups::capture(
                &catalog_digest::builtin_catalog().unwrap(),
                &restarted,
                &reference_only,
                "No material shortcut",
                &retained.report
            )
            .is_err()
        );
    }
}

#[test]
fn saved_setup_requires_exact_origin_event_and_retained_ancestry() {
    for derived in [false, true] {
        let (frame, recipe) = saved_material_recipe(derived);
        let report = recipe
            .inspect(&frame.revisions[&recipe.source.revision])
            .unwrap();
        let mut no_material = report.clone();
        no_material.source_material = None;
        no_material.incomplete = true;
        no_material.inspection.products_determined = false;
        no_material.inspection.fragments.clear();
        no_material.fragment_modifications.clear();
        no_material
            .validate_binding(&no_material.context_revision)
            .unwrap();
        assert!(
            digest_setups::capture(
                &catalog_digest::builtin_catalog().unwrap(),
                &frame,
                &recipe,
                "Saved replay cannot be recognition-only",
                &no_material,
            )
            .is_err()
        );
        let before = workspace::save_workspace(&frame).unwrap();
        let mut wrong_event = recipe.clone();
        wrong_event.source.event_id = binding(&frame, 0).event_id;
        assert!(wrong_event.validate_for_workspace(&frame).is_err());
        assert!(
            digest_setups::capture(
                &catalog_digest::builtin_catalog().unwrap(),
                &frame,
                &wrong_event,
                "Wrong event",
                &report
            )
            .is_err()
        );
        let mut missing_parent = frame.clone();
        let parent_event = binding(&frame, 0).event_id;
        missing_parent
            .history
            .retain(|event| event.id != parent_event);
        assert!(recipe.validate_for_workspace(&missing_parent).is_err());
        assert_eq!(workspace::save_workspace(&frame).unwrap(), before);
    }
}
