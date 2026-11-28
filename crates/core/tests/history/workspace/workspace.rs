//! Archive fixtures include an independently written Stored-ZIP encoder, not just
//! save/load against the same implementation. No filesystem extraction is used.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::{
        history::transition::Transition, reference_edit::Operation as ReferenceEdit,
        transfer::Candidate, workspace::*,
    },
    material::lineage::EventBinding,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
};
use collagene_model::{Error, Feature, Location, Molecule, Topology};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    io::{Cursor, Read, Write},
};
use zip::{CompressionMethod, ZipArchive, ZipWriter, write::SimpleFileOptions};
fn hash(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}
fn crc(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for b in bytes {
        crc ^= *b as u32;
        for _ in 0..8 {
            crc = (crc >> 1) ^ if crc & 1 != 0 { 0xedb88320 } else { 0 };
        }
    }
    !crc
}
fn u16v(out: &mut Vec<u8>, v: u16) {
    out.extend(v.to_le_bytes());
}
fn u32v(out: &mut Vec<u8>, v: u32) {
    out.extend(v.to_le_bytes());
}
fn stored(entries: &[(String, Vec<u8>)]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut central = Vec::new();
    for (name, data) in entries {
        let offset = out.len() as u32;
        let check = crc(data);
        out.extend(b"PK\x03\x04");
        for v in [20, 0, 0, 0, 33] {
            u16v(&mut out, v);
        }
        for v in [check, data.len() as u32, data.len() as u32] {
            u32v(&mut out, v);
        }
        u16v(&mut out, name.len() as u16);
        u16v(&mut out, 0);
        out.extend(name.as_bytes());
        out.extend(data);
        central.extend(b"PK\x01\x02");
        for v in [20, 20, 0, 0, 0, 33] {
            u16v(&mut central, v);
        }
        for v in [check, data.len() as u32, data.len() as u32] {
            u32v(&mut central, v);
        }
        for v in [name.len() as u16, 0, 0, 0, 0] {
            u16v(&mut central, v);
        }
        u32v(&mut central, 0);
        u32v(&mut central, offset);
        central.extend(name.as_bytes());
    }
    let start = out.len();
    out.extend(&central);
    out.extend(b"PK\x05\x06");
    for v in [0, 0, entries.len() as u16, entries.len() as u16] {
        u16v(&mut out, v);
    }
    u32v(&mut out, central.len() as u32);
    u32v(&mut out, start as u32);
    u16v(&mut out, 0);
    out
}
fn unpack(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut zip = ZipArchive::new(Cursor::new(bytes)).unwrap();
    (0..zip.len())
        .map(|i| {
            let mut f = zip.by_index(i).unwrap();
            let name = f.name().into();
            let mut v = Vec::new();
            f.read_to_end(&mut v).unwrap();
            (name, v)
        })
        .collect()
}
fn mutate_manifest(bytes: &[u8], edit: impl FnOnce(&mut Value)) -> Vec<u8> {
    let mut files = unpack(bytes);
    let entry = files
        .iter_mut()
        .find(|(n, _)| n == "manifest.json")
        .unwrap();
    let mut value: Value = serde_json::from_slice(&entry.1).unwrap();
    edit(&mut value);
    entry.1 = serde_json::to_vec(&value).unwrap();
    stored(&files)
}
fn simple() -> Workspace {
    Workspace::from_documents(
        "Simple",
        vec![Molecule::dna("sequence", "AAAACCCCGGGG", Topology::Circular).unwrap()],
        0,
    )
    .unwrap()
}
fn manual_fixture() -> Vec<u8> {
    let name = b"manual";
    let bases = b"ACGT";
    let features = b"[]";
    let transition = br#""IndependentReference""#;
    let molecule_json=br#"{"name":"manual","sequence":{"alphabet":"Dna","bases":"ACGT"},"topology":"Circular","features":[],"metadata":{},"original_genbank":null}"#;
    let revision = hash(molecule_json);
    let record=serde_json::to_vec(&json!({"name":hash(name),"sequence":hash(bases),"alphabet":"Dna","topology":"Circular","features":hash(features),"metadata":{},"original_genbank":null})).unwrap();
    let mut blobs = BTreeMap::new();
    for bytes in [
        name.as_slice(),
        bases.as_slice(),
        features.as_slice(),
        record.as_slice(),
        transition.as_slice(),
    ] {
        blobs.insert(hash(bytes), bytes.to_vec());
    }
    // Independently specify a current native history root and its content ID.
    let mut event_digest = blake3::Hasher::new();
    event_digest.update(b"collagene-history-event-v2\0");
    event_digest.update(
        &serde_json::to_vec(&(
            "manual-doc",
            &revision,
            Vec::<String>::new(),
            "Imported manual fixture",
            "IndependentReference",
        ))
        .unwrap(),
    );
    let event_id = format!("event-{}", event_digest.finalize().to_hex());
    let manifest=serde_json::to_vec(&json!({"format":"collagene-workspace","version":2,"molecule_schema":1,"plan_schema":PLAN_SCHEMA_VERSION,"name":"Manual fixture","active_document":"manual-doc","documents":[{"id":"manual-doc","head":revision,"undo":[],"redo":[],"head_event":event_id,"undo_events":[],"redo_events":[]}],"history":[{"id":event_id,"document_id":"manual-doc","revision":revision,"parents":[],"operation":"Imported manual fixture","transition_blob":hash(transition)}],"revisions":{revision:hash(&record)},"candidates":[],"unresolved_candidate_targets":[],"attachments":{},"blobs":blobs.iter().map(|(h,b)|(h.clone(),json!({"length":b.len()}))).collect::<BTreeMap<_,_>>()})).unwrap();
    let mut entries = vec![("manifest.json".into(), manifest)];
    entries.extend(blobs.into_iter().map(|(h, b)| (format!("blobs/{h}"), b)));
    stored(&entries)
}
#[test]
fn independent_stored_zip_and_manifest_fixture_recovers_exact_molecule() {
    let bytes = manual_fixture();
    assert!(is_workspace_archive(&bytes));
    let workspace = load_workspace(&bytes).unwrap();
    assert_eq!(workspace.name, "Manual fixture");
    assert_eq!(
        workspace.head("manual-doc").unwrap().sequence.as_str(),
        "ACGT"
    );
    assert_eq!(
        workspace.head("manual-doc").unwrap().topology,
        Topology::Circular
    );
    assert_eq!(workspace.history.len(), 1);
    assert_eq!(
        workspace.documents[0].head_event.as_deref(),
        Some(workspace.history[0].id.as_str())
    );
    assert_eq!(
        load_workspace(&save_workspace(&workspace).unwrap())
            .unwrap()
            .head("manual-doc")
            .unwrap()
            .revision(),
        workspace.head("manual-doc").unwrap().revision()
    );
}
#[test]
fn current_native_archive_requires_explicit_exact_navigation_events() {
    let bytes = manual_fixture();
    let missing = mutate_manifest(&bytes, |m| m["documents"][0]["head_event"] = Value::Null);
    assert!(load_workspace(&missing).is_err());
    for stack in ["undo", "redo"] {
        let missing = mutate_manifest(&bytes, |m| {
            let revision = m["documents"][0]["head"].clone();
            m["documents"][0][stack] = json!([revision]);
        });
        assert!(
            load_workspace(&missing).is_err(),
            "Revision-only {stack} stack was accepted"
        );
    }
}

#[test]
fn branch_events_repeated_content_exact_navigation_and_opaque_evidence_survive() {
    let mut root = Molecule::dna("Original", "AAAACCCCGGGG", Topology::Circular).unwrap();
    root.metadata.insert(
        "trace_original_base64".into(),
        STANDARD.encode([0, 128, 255, 7, 9]),
    );
    root.metadata
        .insert("unknown".into(), "line one\nline two\0end".into());
    root.original_genbank = Some("Original record evidence\n//\n".into());
    root.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(10, 12),
            Location::range(0, 2),
        ]))),
        qualifiers: vec![
            ("note".into(), Some("first\nsecond".into())),
            ("note".into(), None),
        ],
        needs_review: true,
    });
    let mut workspace =
        Workspace::from_documents("Branches", vec![root.clone(), root.clone()], 0).unwrap();
    assert_eq!(workspace.revisions.len(), 1);
    assert_ne!(workspace.documents[0].id, workspace.documents[1].id);
    let doc = workspace.documents[0].id.clone();
    let e0 = workspace.history[0].id.clone();
    let one = root.replace(4, 5, "T").unwrap();
    let two = root.replace(5, 6, "A").unwrap();
    let e1 = RevisionEvent::new(
        &doc,
        one.revision(),
        vec![e0.clone()],
        "Replace 4..5 with T",
        Transition::ReferenceEdit {
            source: EventBinding {
                document_id: doc.clone(),
                event_id: e0.clone(),
                revision: root.revision(),
            },
            actions: vec![ReferenceEdit::ReplaceSequence {
                start: 4,
                end: 5,
                replacement: "T".into(),
            }],
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    let e2 = RevisionEvent::new(
        &doc,
        two.revision(),
        vec![e0.clone()],
        "Alternative branch",
        Transition::ReferenceEdit {
            source: EventBinding {
                document_id: doc.clone(),
                event_id: e0.clone(),
                revision: root.revision(),
            },
            actions: vec![ReferenceEdit::ReplaceSequence {
                start: 5,
                end: 6,
                replacement: "A".into(),
            }],
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    let e3 = RevisionEvent::new(
        &doc,
        root.revision(),
        vec![e1.id.clone()],
        "Return to original bases",
        Transition::ReferenceEdit {
            source: EventBinding {
                document_id: doc.clone(),
                event_id: e1.id.clone(),
                revision: one.revision(),
            },
            actions: vec![ReferenceEdit::ReplaceSequence {
                start: 4,
                end: 5,
                replacement: "C".into(),
            }],
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    assert_ne!(e0, e3.id);
    workspace.revisions.insert(one.revision(), one.clone());
    workspace.revisions.insert(two.revision(), two.clone());
    workspace.documents[0].head = root.revision();
    workspace.documents[0].head_event = Some(e3.id.clone());
    workspace.documents[0].undo = vec![root.revision(), one.revision()];
    workspace.documents[0].undo_events = vec![e0, e1.id.clone()];
    workspace.documents[0].redo = vec![two.revision()];
    workspace.documents[0].redo_events = vec![e2.id.clone()];
    workspace.history.extend([e1, e2, e3]);
    let attachment = WorkspaceAttachment {
        kind: "oligo_table".into(),
        name: "oligos.csv".into(),
        media_type: "application/json".into(),
        bytes: br#"{"rows":["ACGT"]}"#.to_vec(),
    };
    workspace
        .attachments
        .insert(hash(&attachment.bytes), attachment);
    let bytes = save_workspace(&workspace).unwrap();
    assert_eq!(bytes, save_workspace(&workspace).unwrap());
    let restored = load_workspace(&bytes).unwrap();
    assert_eq!(
        serde_json::to_value(&restored).unwrap(),
        serde_json::to_value(&workspace).unwrap()
    );
    assert!(
        restored
            .head(&doc)
            .unwrap()
            .original_genbank
            .as_ref()
            .unwrap()
            .ends_with("//\n")
    );
}
#[test]
fn candidate_snapshots_deduplicate_sources_and_preserve_unresolved_targets_as_data() {
    let mut workspace = simple();
    let source = workspace.head("document-1").unwrap().clone();
    let candidate = Candidate::new(DesignPlan {
        name: "Unbound draft".into(),
        target_revision: "f".repeat(64),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("source".into())],
            overlaps: vec![1],
            circular: true,
        }],
    });
    workspace.candidates.push(candidate.clone());
    workspace.active_candidate = Some(candidate.plan_revision.clone());
    assert_eq!(
        workspace.unresolved_candidate_targets(),
        vec!["f".repeat(64)]
    );
    let bytes = save_workspace(&workspace).unwrap();
    let restored = load_workspace(&bytes).unwrap();
    assert_eq!(
        restored.candidates[0].plan_revision,
        candidate.plan_revision
    );
    assert_eq!(
        restored.active_candidate,
        Some(candidate.plan_revision.clone())
    );
    assert_eq!(restored.revisions.len(), 1);
    assert_eq!(
        restored.unresolved_candidate_targets(),
        vec!["f".repeat(64)]
    );
    let inconsistent = mutate_manifest(&bytes, |m| m["unresolved_candidate_targets"] = json!([]));
    assert!(load_workspace(&inconsistent).is_err());
}
#[test]
fn large_original_attachments_are_stored_once_across_revision_history() {
    let evidence = (0..384 * 1024)
        .map(|i| (i * 73 % 251) as u8)
        .collect::<Vec<_>>();
    let encoded = STANDARD.encode(&evidence);
    let mut source = Molecule::dna("r0", "ACGTACGT", Topology::Linear).unwrap();
    source
        .metadata
        .insert("trace_original_base64".into(), encoded.clone());
    source.original_genbank = Some("original\n".repeat(10_000));
    let mut workspace =
        Workspace::from_documents("Evidence dedup", vec![source.clone()], 0).unwrap();
    let mut parent = workspace.history[0].id.clone();
    for i in 1..8 {
        let old_event = parent.clone();
        let mut next = source.clone();
        next.name = format!("r{i}");
        let binding = EventBinding {
            document_id: "document-1".into(),
            event_id: parent.clone(),
            revision: workspace.documents[0].head.clone(),
        };
        let event = RevisionEvent::new(
            "document-1",
            next.revision(),
            vec![parent],
            "Rename",
            Transition::ReferenceEdit {
                source: binding,
                actions: vec![ReferenceEdit::Rename {
                    name: next.name.clone(),
                }],
                supplemental_sources: vec![],
            },
        )
        .unwrap();
        parent = event.id.clone();
        let old_head = workspace.documents[0].head.clone();
        workspace.documents[0].undo.push(old_head);
        workspace.documents[0].undo_events.push(old_event);
        workspace.documents[0].head = next.revision();
        workspace.revisions.insert(next.revision(), next);
        workspace.history.push(event);
    }
    workspace.documents[0].head_event = Some(parent);
    let bytes = save_workspace(&workspace).unwrap();
    assert!(bytes.len() < 600_000, "archive has {} bytes", bytes.len());
    let restored = load_workspace(&bytes).unwrap();
    assert_eq!(restored.revisions.len(), 8);
    assert_eq!(restored.documents[0].undo_events.len(), 7);
    let mut navigation =
        collagene_core::documents::history::DocumentHistory::restore_all(&restored)
            .unwrap()
            .remove(0);
    for _ in 0..7 {
        navigation.undo().unwrap();
    }
    assert_eq!(navigation.molecule.name, "r0");
    assert!(
        restored
            .revisions
            .values()
            .all(|m| m.metadata["trace_original_base64"] == encoded)
    );
    let entries = unpack(&bytes);
    assert_eq!(
        entries
            .iter()
            .filter(|(_, bytes)| bytes == &evidence)
            .count(),
        1
    );
}
#[test]
fn invalid_revision_graph_navigation_and_attachment_hashes_fail_before_save() {
    let original = simple();
    let mut workspace = original.clone();
    workspace.active_document = "missing".into();
    assert!(save_workspace(&workspace).is_err());
    let mut workspace = original.clone();
    workspace.documents[0].head = "0".repeat(64);
    assert!(save_workspace(&workspace).is_err());
    let mut workspace = original.clone();
    workspace.history[0]
        .parents
        .push(format!("event-{}", "0".repeat(64)));
    assert!(save_workspace(&workspace).is_err());
    let mut workspace = original.clone();
    workspace.history.push(workspace.history[0].clone());
    assert!(save_workspace(&workspace).is_err());
    let mut workspace = original.clone();
    workspace.documents[0].undo_events = vec![workspace.history[0].id.clone()];
    assert!(save_workspace(&workspace).is_err());
    let mut workspace = original.clone();
    workspace.attachments.insert(
        "0".repeat(64),
        WorkspaceAttachment {
            kind: "evidence".into(),
            name: "bytes".into(),
            media_type: "application/octet-stream".into(),
            bytes: vec![1],
        },
    );
    assert!(save_workspace(&workspace).is_err());
    let mut workspace = original;
    let revision = workspace.documents[0].head.clone();
    workspace.revisions.get_mut(&revision).unwrap().name = "Tampered".into();
    assert!(save_workspace(&workspace).is_err());
}
#[test]
fn manifest_versions_unknown_fields_duplicate_json_and_blob_tampering_are_rejected() {
    let bytes = save_workspace(&simple()).unwrap();
    let future = mutate_manifest(&bytes, |m| m["version"] = json!(3));
    assert!(is_workspace_archive(&future));
    assert!(matches!(
        load_workspace(&future),
        Err(Error::Unsupported(_))
    ));
    let unknown = mutate_manifest(&bytes, |m| m["authority"] = json!("verified"));
    assert!(load_workspace(&unknown).is_err());
    let mut files = unpack(&bytes);
    let manifest = &mut files
        .iter_mut()
        .find(|(n, _)| n == "manifest.json")
        .unwrap()
        .1;
    let mut duplicate = b"{\"version\":999,".to_vec();
    duplicate.extend(&manifest[1..]);
    *manifest = duplicate;
    assert!(load_workspace(&stored(&files)).is_err());
    let mut files = unpack(&bytes);
    files
        .iter_mut()
        .find(|(n, _)| n.starts_with("blobs/"))
        .unwrap()
        .1
        .push(b'X');
    assert!(load_workspace(&stored(&files)).is_err());
    let extra = mutate_manifest(&bytes, |m| m["blobs"]["0".repeat(64)] = json!({"length":0}));
    assert!(load_workspace(&extra).is_err());
}
#[test]
fn zip_duplicate_names_paths_hidden_local_data_encryption_and_extreme_offsets_reject() {
    let bytes = save_workspace(&simple()).unwrap();
    let entries = unpack(&bytes);
    let mut duplicate = entries.clone();
    duplicate.push(duplicate[0].clone());
    assert!(load_workspace(&stored(&duplicate)).is_err());
    let mut path = entries.clone();
    path.push(("../escape".into(), vec![1]));
    assert!(load_workspace(&stored(&path)).is_err());
    let mut encrypted = bytes.clone();
    let directory = encrypted
        .windows(4)
        .rposition(|s| s == b"PK\x05\x06")
        .unwrap();
    let central = u32::from_le_bytes(
        encrypted[directory + 16..directory + 20]
            .try_into()
            .unwrap(),
    ) as usize;
    encrypted[central + 8..central + 10].copy_from_slice(&1u16.to_le_bytes());
    encrypted[6..8].copy_from_slice(&1u16.to_le_bytes());
    assert!(load_workspace(&encrypted).is_err());
    let mut wrong = bytes.clone();
    wrong[central + 42..central + 46].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(load_workspace(&wrong).is_err());
    let mut too_large = bytes.clone();
    too_large[central + 24..central + 28]
        .copy_from_slice(&((MAX_MANIFEST_BYTES + 1) as u32).to_le_bytes());
    assert!(matches!(load_workspace(&too_large), Err(Error::TooLarge)));
    let mut leading = vec![0];
    leading.extend(bytes);
    assert!(load_workspace(&leading).is_err());
}
#[test]
fn bounded_deflate_input_is_supported_but_bad_crc_and_symlinks_are_not() {
    let source = save_workspace(&simple()).unwrap();
    let files = unpack(&source);
    let mut writer = ZipWriter::new(Cursor::new(Vec::new()));
    for (name, data) in &files {
        writer
            .start_file(
                name,
                SimpleFileOptions::default().compression_method(CompressionMethod::Deflated),
            )
            .unwrap();
        writer.write_all(data).unwrap();
    }
    let deflated = writer.finish().unwrap().into_inner();
    assert_eq!(
        load_workspace(&deflated)
            .unwrap()
            .head("document-1")
            .unwrap()
            .sequence
            .as_str(),
        "AAAACCCCGGGG"
    );
    let mut bad_crc = source.clone();
    let local = b"manifest.json";
    let at = bad_crc
        .windows(local.len())
        .position(|b| b == local)
        .unwrap()
        + local.len();
    bad_crc[at] ^= 1;
    assert!(load_workspace(&bad_crc).is_err());
    let mut link = source;
    let footer = link.windows(4).rposition(|s| s == b"PK\x05\x06").unwrap();
    let central = u32::from_le_bytes(link[footer + 16..footer + 20].try_into().unwrap()) as usize;
    link[central + 38..central + 42].copy_from_slice(&(0o120777u32 << 16).to_le_bytes());
    assert!(load_workspace(&link).is_err());
}
#[test]
fn truncation_and_workspace_dimensions_are_bounded() {
    let bytes = manual_fixture();
    for end in 0..bytes.len() {
        assert!(
            load_workspace(&bytes[..end]).is_err(),
            "accepted prefix {end}"
        );
    }
    let mut workspace = simple();
    workspace.documents[0].undo = vec![workspace.documents[0].head.clone(); 1001];
    assert!(save_workspace(&workspace).is_err());
    assert!(Workspace::from_documents("empty", vec![], 0).is_err());
    assert!(
        RevisionEvent::new(
            "../bad",
            "0".repeat(64),
            vec![],
            "x",
            Transition::IndependentReference
        )
        .is_err()
    );
    assert!(!is_workspace_archive(b">not a zip\nACGT\n"));
}

#[test]
fn listed_but_unreachable_blobs_and_silently_ignored_feature_fields_are_rejected() {
    let bytes = save_workspace(&simple()).unwrap();
    let mut files = unpack(&bytes);
    let extra = b"hidden but correctly hashed evidence".to_vec();
    let extra_hash = hash(&extra);
    let manifest_entry = files
        .iter_mut()
        .find(|(name, _)| name == "manifest.json")
        .unwrap();
    let mut manifest: Value = serde_json::from_slice(&manifest_entry.1).unwrap();
    manifest["blobs"][&extra_hash] = json!({"length": extra.len()});
    manifest_entry.1 = serde_json::to_vec(&manifest).unwrap();
    files.push((format!("blobs/{extra_hash}"), extra));
    assert!(load_workspace(&stored(&files)).is_err());

    let mut molecule = Molecule::dna("Feature", "ACGT", Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(0, 4),
        qualifiers: vec![],
        needs_review: false,
    });
    let workspace = Workspace::from_documents("Unknown typed field", vec![molecule], 0).unwrap();
    let mut files: BTreeMap<_, _> = unpack(&save_workspace(&workspace).unwrap())
        .into_iter()
        .collect();
    let mut manifest: Value =
        serde_json::from_slice(&files.remove("manifest.json").unwrap()).unwrap();
    let revision = workspace.documents[0].head.clone();
    let record_hash = manifest["revisions"][&revision]
        .as_str()
        .unwrap()
        .to_owned();
    let mut record: Value =
        serde_json::from_slice(&files.remove(&format!("blobs/{record_hash}")).unwrap()).unwrap();
    let feature_hash = record["features"].as_str().unwrap().to_owned();
    let mut features: Value =
        serde_json::from_slice(&files.remove(&format!("blobs/{feature_hash}")).unwrap()).unwrap();
    features[0]["future_semantics"] = json!("must not disappear");
    let new_features = serde_json::to_vec(&features).unwrap();
    let new_features_hash = hash(&new_features);
    record["features"] = json!(new_features_hash);
    let new_record = serde_json::to_vec(&record).unwrap();
    let new_record_hash = hash(&new_record);
    manifest["revisions"][&revision] = json!(new_record_hash);
    let blob_manifest = manifest["blobs"].as_object_mut().unwrap();
    blob_manifest.remove(&record_hash);
    blob_manifest.remove(&feature_hash);
    blob_manifest.insert(
        new_features_hash.clone(),
        json!({"length": new_features.len()}),
    );
    blob_manifest.insert(new_record_hash.clone(), json!({"length": new_record.len()}));
    files.insert(format!("blobs/{new_features_hash}"), new_features);
    files.insert(format!("blobs/{new_record_hash}"), new_record);
    files.insert(
        "manifest.json".into(),
        serde_json::to_vec(&manifest).unwrap(),
    );
    let error = load_workspace(&stored(&files.into_iter().collect::<Vec<_>>())).unwrap_err();
    assert!(
        error.to_string().contains("unsupported or noncanonical"),
        "{error}"
    );
}

#[test]
fn repair_candidate_roundtrips_with_only_the_current_workspace_plan_schema() {
    let mut workspace = simple();
    let source = workspace.head(&workspace.active_document).unwrap().clone();
    workspace.candidates.push(Candidate::new(DesignPlan {
        name: "Repair candidate".into(),
        target_revision: source.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::RepairEnds {
            template: Source::Input("source".into()),
            method: collagene_core::chemistry::EndRepairMethod::Blunt,
        }],
    }));
    let current = save_workspace(&workspace).unwrap();
    let restored = load_workspace(&current).unwrap();
    assert_eq!(
        serde_json::to_value(&restored.candidates).unwrap(),
        serde_json::to_value(&workspace.candidates).unwrap()
    );
    for unsupported in [PLAN_SCHEMA_VERSION - 1, PLAN_SCHEMA_VERSION + 1] {
        let mislabeled = mutate_manifest(&current, |manifest| {
            manifest["plan_schema"] = json!(unsupported)
        });
        assert!(matches!(
            load_workspace(&mislabeled),
            Err(Error::Unsupported(_))
        ));
    }
}
#[test]
fn tailing_candidate_preserves_current_schema_and_literal_operation() {
    let mut w = simple();
    let source = w.head(&w.active_document).unwrap().clone();
    w.candidates.push(Candidate::new(DesignPlan {
        name: "tail".into(),
        target_revision: source.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Tail3Prime {
            template: Source::Input("s".into()),
            base: collagene_core::cloning::ligation::tailing::TailingBase::A,
        }],
    }));
    let b = save_workspace(&w).unwrap();
    let restored = load_workspace(&b).unwrap();
    assert_eq!(
        serde_json::to_value(&restored.candidates).unwrap(),
        serde_json::to_value(&w.candidates).unwrap()
    );
    for version in [PLAN_SCHEMA_VERSION - 1, PLAN_SCHEMA_VERSION + 1] {
        let changed = mutate_manifest(&b, |m| m["plan_schema"] = json!(version));
        assert!(matches!(
            load_workspace(&changed),
            Err(Error::Unsupported(_))
        ));
    }
}

#[test]
fn native_zip_preserves_import_scale_features() {
    let mut molecule = Molecule::dna("Import-scale annotation", "ACGT", Topology::Linear).unwrap();
    molecule.features = vec![
        Feature {
            kind: "x".into(),
            location: Location::range(0, 1),
            qualifiers: vec![],
            needs_review: false,
        };
        100_000
    ];
    let workspace = Workspace::from_documents("Annotation archive", vec![molecule], 0).unwrap();
    let bytes = save_workspace(&workspace).unwrap();
    let loaded = load_workspace(&bytes).unwrap();
    assert_eq!(
        loaded.head(&loaded.active_document).unwrap().features.len(),
        100_000
    );
    assert_eq!(save_workspace(&loaded).unwrap(), bytes);
}

#[test]
fn native_zip_preserves_full_document_and_retained_revision_capacities() {
    let molecule = Molecule::dna("Source", "ACGT", Topology::Linear).unwrap();
    let mut workspace =
        Workspace::from_documents("Retained history", vec![molecule.clone(); MAX_DOCUMENTS], 0)
            .unwrap();
    for index in 1..MAX_REVISIONS {
        let mut retained = molecule.clone();
        retained.name = format!("Retained revision {index}");
        workspace.revisions.insert(retained.revision(), retained);
    }
    let bytes = save_workspace(&workspace).unwrap();
    let loaded = load_workspace(&bytes).unwrap();
    assert_eq!(loaded.documents.len(), MAX_DOCUMENTS);
    assert_eq!(loaded.revisions.len(), MAX_REVISIONS);
    assert_eq!(save_workspace(&loaded).unwrap(), bytes);
}

#[test]
fn native_v2_history_uses_deduplicated_closed_transition_blob_references() {
    let source = Molecule::dna("Two independent references", "ACGT", Topology::Linear).unwrap();
    let workspace =
        Workspace::from_documents("Transition blobs", vec![source.clone(), source], 0).unwrap();
    let bytes = save_workspace(&workspace).unwrap();
    let files = unpack(&bytes);
    let manifest: Value =
        serde_json::from_slice(&files.iter().find(|(n, _)| n == "manifest.json").unwrap().1)
            .unwrap();
    assert_eq!(manifest["version"], 2);
    let first = &manifest["history"][0]["transition_blob"];
    assert_eq!(first, &manifest["history"][1]["transition_blob"]);
    assert!(manifest["history"][0].get("transition").is_none());
    let transition = br#""IndependentReference""#;
    let expected_transition_hash = hash(transition);
    assert_eq!(first.as_str(), Some(expected_transition_hash.as_str()));
    assert_eq!(
        files
            .iter()
            .filter(|(_, b)| b.as_slice() == transition)
            .count(),
        1
    );
    let reopened = load_workspace(&bytes).unwrap();
    assert_eq!(
        serde_json::to_value(&reopened.history).unwrap(),
        serde_json::to_value(&workspace.history).unwrap()
    );
    let missing = mutate_manifest(&bytes, |manifest| {
        manifest["history"][0]["transition_blob"] = json!("0".repeat(64))
    });
    assert!(load_workspace(&missing).is_err());
    let inline = mutate_manifest(&bytes, |manifest| {
        manifest["history"][0]["transition"] = json!("IndependentReference");
    });
    assert!(load_workspace(&inline).is_err());
    let unsupported = mutate_manifest(&bytes, |manifest| manifest["version"] = json!(1));
    assert!(matches!(
        load_workspace(&unsupported),
        Err(Error::Unsupported(_))
    ));

    // Rehash an unsupported transition payload and close every blob reference.
    // Hash-valid ZIP contents still cannot introduce an unknown transition.
    let mut files = files.into_iter().collect::<BTreeMap<_, _>>();
    let mut manifest: Value =
        serde_json::from_slice(&files.remove("manifest.json").unwrap()).unwrap();
    let original = first.as_str().unwrap();
    files.remove(&format!("blobs/{original}"));
    manifest["blobs"].as_object_mut().unwrap().remove(original);
    let unknown=br#"{"Analysis":{"origin":{"producer":"CodonDesign","record":{"Context":{"revision":"0000000000000000000000000000000000000000000000000000000000000000"}},"authority":true}}}"#.to_vec();
    let altered = hash(&unknown);
    for entry in manifest["history"].as_array_mut().unwrap() {
        entry["transition_blob"] = json!(altered);
    }
    manifest["blobs"][&altered] = json!({"length":unknown.len()});
    files.insert(format!("blobs/{altered}"), unknown);
    files.insert(
        "manifest.json".into(),
        serde_json::to_vec(&manifest).unwrap(),
    );
    assert!(load_workspace(&stored(&files.into_iter().collect::<Vec<_>>())).is_err());
}

#[test]
fn typed_reference_action_larger_than_manifest_limit_roundtrips_in_its_own_blob() {
    use collagene_core::documents::history::DocumentHistory;
    let source = Molecule::dna("Large explicit reference edit", "ACGT", Topology::Linear).unwrap();
    let mut workspace =
        Workspace::from_documents("Large action archive", vec![source.clone()], 0).unwrap();
    let mut document = DocumentHistory::new("large-action", source).unwrap();
    // Preserve an actual user literal above16MiB while normalization produces
    // one valid base, staying within the native10M sequence bound.
    document
        .commit(
            vec![ReferenceEdit::ReplaceSequence {
                start: 0,
                end: 4,
                replacement: format!("{}T", " ".repeat(MAX_MANIFEST_BYTES + 1)),
            }],
            "Literal large replacement",
        )
        .unwrap();
    workspace.documents.clear();
    workspace.revisions.clear();
    workspace.history.clear();
    document.contribute(&mut workspace).unwrap();
    workspace.active_document = document.id.clone();
    let bytes = save_workspace(&workspace).unwrap();
    let files = unpack(&bytes);
    let manifest = &files
        .iter()
        .find(|(name, _)| name == "manifest.json")
        .unwrap()
        .1;
    assert!(
        manifest.len() < 4096,
        "typed action bytes must stay out of the manifest"
    );
    let manifest: Value = serde_json::from_slice(manifest).unwrap();
    let event = manifest["history"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["operation"] == "Literal large replacement")
        .unwrap();
    let path = format!("blobs/{}", event["transition_blob"].as_str().unwrap());
    let transition = &files.iter().find(|(name, _)| name == &path).unwrap().1;
    assert!(transition.len() > MAX_MANIFEST_BYTES && transition.len() <= MAX_ENTRY_BYTES);
    drop(files);
    let loaded = load_workspace(&bytes).unwrap();
    assert_eq!(loaded.head("large-action").unwrap(), &document.molecule);
    assert_eq!(
        serde_json::to_value(&loaded.history).unwrap(),
        serde_json::to_value(&workspace.history).unwrap()
    );
}
