use collagene_core::{
    Molecule, Topology, catalog_digest,
    enzyme_sets::{self, Set},
    restriction_inspection::EnzymeSpec,
    workspace::{Workspace, WorkspaceAttachment},
};
fn frame() -> Workspace {
    Workspace::from_documents(
        "enzyme preferences",
        vec![Molecule::dna("source", "GAATTCAAGCTT", Topology::Circular).unwrap()],
        0,
    )
    .unwrap()
}
fn choices() -> Vec<EnzymeSpec> {
    vec![
        EnzymeSpec {
            name: "HindIII".into(),
            recognition: "AAGCTT".into(),
            top_cut: 1,
            bottom_cut: 5,
        },
        EnzymeSpec {
            name: "EcoRI".into(),
            recognition: "GAATTC".into(),
            top_cut: 1,
            bottom_cut: 5,
        },
    ]
}
fn encoded(frame: &Workspace, set: &serde_json::Value) -> (Workspace, String) {
    let bytes = serde_json::to_vec(set).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut staged = frame.clone();
    staged.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: enzyme_sets::KIND.into(),
            name: set["name"].as_str().unwrap().into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    (staged, hash)
}
#[test]
fn names_replace_atomically_preserving_order_and_every_other_workspace_field() {
    let original = frame();
    let set = enzyme_sets::build("  Common cutters  ", &choices(), "catalog-a").unwrap();
    let (stored, hash) = enzyme_sets::store(&original, &set).unwrap();
    assert_eq!(enzyme_sets::load(&stored, &hash).unwrap(), set);
    assert_eq!(
        enzyme_sets::load(&stored, &hash).unwrap().choices,
        choices()
    );
    assert_eq!(
        enzyme_sets::list(&stored).unwrap(),
        vec![(hash.clone(), set.clone())]
    );
    let (same, same_hash) = enzyme_sets::store(&stored, &set).unwrap();
    assert_eq!(hash, same_hash);
    assert_eq!(
        serde_json::to_value(&same).unwrap(),
        serde_json::to_value(&stored).unwrap()
    );
    let empty = enzyme_sets::build("COMMON CUTTERS", &[], "catalog-a").unwrap();
    let (replaced, new_hash) = enzyme_sets::store(&stored, &empty).unwrap();
    assert_eq!(replaced.attachments.len(), 1);
    assert!(!replaced.attachments.contains_key(&hash));
    assert!(
        enzyme_sets::load(&replaced, &new_hash)
            .unwrap()
            .choices
            .is_empty()
    );
    let deleted = enzyme_sets::delete(&replaced, &new_hash).unwrap();
    assert_eq!(
        serde_json::to_value(&original).unwrap(),
        serde_json::to_value(&deleted).unwrap()
    );
    let mut without_preferences = stored.clone();
    without_preferences.attachments.clear();
    assert_eq!(
        serde_json::to_value(&original).unwrap(),
        serde_json::to_value(without_preferences).unwrap()
    );
    assert!(original.attachments.is_empty());
}
#[test]
fn missing_changed_and_unchanged_geometry_are_distinct_even_after_revision_change() {
    let set = enzyme_sets::build("two", &choices(), "old").unwrap();
    let mut changed = choices();
    changed[0].top_cut = 2;
    changed.pop();
    let r = enzyme_sets::resolve(&set, &changed, "new").unwrap();
    assert_eq!(r.missing, vec!["EcoRI"]);
    assert_eq!(r.changed, vec!["HindIII"]);
    assert!(r.available_names.is_empty());
    assert!(r.catalog_changed);
    let r = enzyme_sets::resolve(&set, &choices(), "new").unwrap();
    assert_eq!(r.available_names, vec!["HindIII", "EcoRI"]);
    assert!(r.changed.is_empty());
    assert!(r.missing.is_empty());
    assert!(r.catalog_changed);
}
#[test]
fn stable_ids_cannot_be_replaced_by_identical_names_or_geometry() {
    let c = catalog_digest::builtin_catalog().unwrap();
    let mut set = enzyme_sets::build_for_catalog("pinned", &choices(), &c).unwrap();
    assert_eq!(
        enzyme_sets::resolve_for_catalog(&set, &c)
            .unwrap()
            .available_names,
        vec!["HindIII", "EcoRI"]
    );
    assert!(enzyme_sets::resolve(&set, &choices(), &enzyme_sets::catalog_revision(&c)).is_err());
    set.catalog_pin.as_mut().unwrap().enzyme_ids[0] = "missing-stable-identity".into();
    let r = enzyme_sets::resolve_for_catalog(&set, &c).unwrap();
    assert_eq!(r.missing, vec!["HindIII"]);
    assert_eq!(r.available_names, vec!["EcoRI"]);
    let mut set = enzyme_sets::build_for_catalog("pinned", &choices(), &c).unwrap();
    set.choices[0].bottom_cut += 1;
    assert_eq!(
        enzyme_sets::resolve_for_catalog(&set, &c).unwrap().changed,
        vec!["HindIII"]
    );
}
#[test]
fn catalog_edition_change_keeps_unchanged_exact_id_matches_available() {
    let c = catalog_digest::builtin_catalog().unwrap();
    let set = enzyme_sets::build_for_catalog("pinned", &choices(), &c).unwrap();
    let mut edition = c.edition().clone();
    edition.edition.push_str("-new");
    let current = catalog_digest::catalog::Catalog::from_edition(edition).unwrap();
    let r = enzyme_sets::resolve_for_catalog(&set, &current).unwrap();
    assert!(r.catalog_changed);
    assert_eq!(r.available_names, vec!["HindIII", "EcoRI"]);
    assert!(r.changed.is_empty());
}
#[test]
fn malformed_or_tampered_records_fail_instead_of_being_hidden() {
    let original = frame();
    let set = enzyme_sets::build("valid", &choices(), "old").unwrap();
    let (mut stored, hash) = enzyme_sets::store(&original, &set).unwrap();
    stored.attachments.get_mut(&hash).unwrap().bytes.push(b' ');
    assert!(enzyme_sets::list(&stored).is_err());
    let mut value = serde_json::to_value(&set).unwrap();
    value["schema_version"] = 99.into();
    let (bad, hash) = encoded(&original, &value);
    assert!(enzyme_sets::load(&bad, &hash).is_err());
    assert!(enzyme_sets::list(&bad).is_err());
    value = serde_json::to_value(&set).unwrap();
    value["choices"][0]["recognition"] = "INVALID".into();
    let (bad, hash) = encoded(&original, &value);
    assert!(enzyme_sets::load(&bad, &hash).is_err());
    value = serde_json::to_value(&set).unwrap();
    value["unexpected"] = true.into();
    let (bad, hash) = encoded(&original, &value);
    assert!(enzyme_sets::load(&bad, &hash).is_err());
    let (mut stored, hash) = enzyme_sets::store(&original, &set).unwrap();
    stored.attachments.get_mut(&hash).unwrap().name = "other".into();
    assert!(enzyme_sets::load(&stored, &hash).is_err());
    assert!(enzyme_sets::delete(&original, &"0".repeat(64)).is_err());
}
#[test]
fn bounded_records_and_named_set_capacity_reject_before_mutation() {
    assert!(enzyme_sets::build(" ", &[], "edition").is_err());
    assert!(enzyme_sets::build(&"a".repeat(257), &[], "edition").is_err());
    let mut duplicate = choices();
    duplicate.push(duplicate[0].clone());
    assert!(enzyme_sets::build("duplicate", &duplicate, "edition").is_err());
    assert!(
        enzyme_sets::build(
            "overlarge",
            &vec![choices()[0].clone(); enzyme_sets::MAX_CHOICES + 1],
            "edition"
        )
        .is_err()
    );
    let mut stored = frame();
    // Construct valid independently hashed attachment records directly, avoiding
    // repeated workspace clones while filling the boundary fixture.
    for i in 0..enzyme_sets::MAX_SETS {
        let set = enzyme_sets::build(&format!("set {i}"), &[], "edition").unwrap();
        let bytes = serde_json::to_vec(&set).unwrap();
        let hash = blake3::hash(&bytes).to_hex().to_string();
        stored.attachments.insert(
            hash,
            WorkspaceAttachment {
                kind: enzyme_sets::KIND.into(),
                name: set.name,
                media_type: "application/json".into(),
                bytes,
            },
        );
    }
    assert_eq!(
        enzyme_sets::list(&stored).unwrap().len(),
        enzyme_sets::MAX_SETS
    );
    assert!(
        enzyme_sets::store(
            &stored,
            &enzyme_sets::build("extra", &[], "edition").unwrap()
        )
        .is_err()
    );
    assert!(
        enzyme_sets::store(
            &stored,
            &enzyme_sets::build("set 0", &choices(), "edition").unwrap()
        )
        .is_ok()
    );
    let invalid = Set {
        schema_version: 1,
        name: "bad".into(),
        choices: vec![],
        catalog_revision: "edition".into(),
        catalog_pin: Some(enzyme_sets::CatalogPin {
            edition: "edition".into(),
            content_hash: "0".repeat(64),
            enzyme_ids: vec![],
        }),
    };
    assert!(enzyme_sets::store(&frame(), &invalid).is_err());
}

#[test]
fn native_archive_roundtrip_and_deletion_leave_unrelated_assets_intact() {
    let mut original = frame();
    let bytes = b"unrelated retained evidence".to_vec();
    let asset_hash = blake3::hash(&bytes).to_hex().to_string();
    original.attachments.insert(
        asset_hash.clone(),
        WorkspaceAttachment {
            kind: "other-evidence".into(),
            name: "evidence".into(),
            media_type: "text/plain".into(),
            bytes,
        },
    );
    let c = catalog_digest::builtin_catalog().unwrap();
    let set = enzyme_sets::build_for_catalog("archive preferences", &choices(), &c).unwrap();
    let (stored, hash) = enzyme_sets::store(&original, &set).unwrap();
    let bytes = collagene_core::workspace::save_workspace(&stored).unwrap();
    let restored = collagene_core::workspace::load_workspace(&bytes).unwrap();
    assert_eq!(enzyme_sets::load(&restored, &hash).unwrap(), set);
    assert!(enzyme_sets::delete(&restored, &asset_hash).is_err());
    let deleted = enzyme_sets::delete(&restored, &hash).unwrap();
    assert_eq!(
        serde_json::to_value(deleted).unwrap(),
        serde_json::to_value(original).unwrap()
    );
}
