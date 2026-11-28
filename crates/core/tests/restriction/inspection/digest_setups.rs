//! Retained digest recipes are historical evidence, never a freshness receipt.
use collagene_core::{
    cloning::restriction::{
        catalog::{self as catalog_digest, CatalogRestrictionInspection, catalog},
        history as restriction_history,
        inspection::InspectionOptions,
        setups::{self as digest_setups, CleavageMode, Record},
    },
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace, WorkspaceAttachment},
    },
    material::lineage::EventBinding,
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
fn fixture() -> (
    catalog::Catalog,
    Workspace,
    EventBinding,
    CatalogRestrictionInspection,
) {
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
    let report = catalog_digest::inspect(
        &catalog,
        &selection,
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &InspectionOptions::default(),
    )
    .unwrap();
    (catalog, frame, source, report)
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
        &source,
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
        &source,
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
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &source, "Original digest", &report).unwrap();
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
        digest_setups::capture(&catalog, &later, &binding(&later, 0), "Wrong head", &report)
            .is_err()
    );
    let alias = binding(&later, 1);
    assert_eq!(alias.revision, source.revision);
    assert!(digest_setups::capture(&catalog, &later, &alias, "Wrong event", &report).is_err());
    assert_eq!(serde_json::to_value(&later).unwrap(), original);
    let alias_report = catalog_digest::inspect(
        &catalog,
        &report.selection,
        &restriction_history::source_role(&alias).unwrap(),
        &later.revisions[&alias.revision],
        &report.inspection.options,
    )
    .unwrap();
    let (both, alias_hash) =
        digest_setups::capture(&catalog, &later, &alias, "Original digest", &alias_report).unwrap();
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
    let report = catalog_digest::inspect(
        &catalog,
        &selected,
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &options,
    )
    .unwrap();
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &source, "Gel lane recipe", &report).unwrap();
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
    let report = catalog_digest::inspect(
        &catalog,
        &prior.selection,
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &options,
    )
    .unwrap();
    assert!(!report.inspection.scan_complete && !report.inspection.products_determined);
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &source, "Incomplete exploration", &report)
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
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &source, "Old catalog recipe", &report).unwrap();
    let mut edition = catalog.edition().clone();
    edition.edition = "future catalog".into();
    let future = catalog::Catalog::from_edition(catalog::seal(edition).unwrap()).unwrap();
    assert!(digest_setups::load(&saved, &hash).is_ok());
    assert!(digest_setups::capture(&future, &saved, &source, "Changed catalog", &report).is_err());
    assert!(
        digest_setups::gel_recipe(&saved, &hash)
            .unwrap()
            .selection
            .validate_for(&future)
            .is_err()
    );
}

#[test]
fn tampered_source_options_geometry_version_mode_and_evidence_flags_fail() {
    let (catalog, frame, source, report) = fixture();
    let (saved, hash) =
        digest_setups::capture(&catalog, &frame, &source, "Recorded digest", &report).unwrap();
    let original = serde_json::to_value(digest_setups::load(&saved, &hash).unwrap()).unwrap();
    let mut cases = Vec::new();
    let mut value = original.clone();
    value["schema_version"] = serde_json::json!(2);
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
    let (one, hash) = digest_setups::capture(&catalog, &frame, &source, "One", &report).unwrap();
    let (two, other) = digest_setups::capture(&catalog, &one, &source, "Two", &report).unwrap();
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
    let (mut full, hash) =
        digest_setups::capture(&catalog, &frame, &source, "Existing", &report).unwrap();
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
            &source,
            "Would exceed native capacity",
            &report
        )
        .is_err()
    );
    assert_eq!(serde_json::to_value(&full).unwrap(), before);
    let (retry, repeated) =
        digest_setups::capture(&catalog, &full, &source, "Existing", &report).unwrap();
    assert_eq!(repeated, hash);
    assert_eq!(retry.attachments.len(), workspace::MAX_ATTACHMENTS);
    for name in ["", " ", "ends with space ", "bad\nname"] {
        assert!(digest_setups::capture(&catalog, &frame, &source, name, &report).is_err());
    }
    assert!(digest_setups::capture(&catalog, &frame, &source, &"x".repeat(257), &report).is_err());
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
    let report = catalog_digest::inspect(
        &catalog,
        &report.selection,
        &restriction_history::source_role(&source).unwrap(),
        &frame.revisions[&source.revision],
        &report.inspection.options,
    )
    .unwrap();
    let (saved, hash) = digest_setups::capture(
        &catalog,
        &frame,
        &source,
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
