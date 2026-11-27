use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    protocol,
    run_lineage::{EventBinding, InputEventBinding},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
use collagene_lab_library as store;
use collagene_library_host::{self as host, ArchiveInventory, PublicationScope};
use collagene_library_protocol as wire;
use std::collections::BTreeSet;
fn native(name: &str, text: &str) -> wire::Payload {
    let bytes =
        protocol::save_document(&Molecule::dna(name, text, Topology::Linear).unwrap()).unwrap();
    wire::Payload {
        kind: wire::PayloadKind::NativeDocument,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: None,
        sharing: None,
    }
}
fn archive(frame: &Workspace, index: usize, sharing: wire::WorkspaceSharing) -> wire::Payload {
    let bytes = workspace::save_workspace(frame).unwrap();
    let doc = &frame.documents[index];
    wire::Payload {
        kind: wire::PayloadKind::WorkspaceArchive,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: Some(wire::NativeSelection {
            document_id: wire::NativeId::new(&doc.id).unwrap(),
            revision: wire::Digest::new(&doc.head).unwrap(),
            head_event: doc
                .head_event
                .as_ref()
                .map(|id| wire::NativeId::new(id).unwrap()),
        }),
        sharing: Some(sharing),
    }
}
fn setup() -> (store::Library, store::LabAccess) {
    let mut db = store::Library::open_in_memory().unwrap();
    let lab = db.create_lab("Joint lab").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    (db, access)
}
fn call(
    db: &mut store::Library,
    access: &store::LabAccess,
    request: wire::Request,
    scope: Option<&PublicationScope>,
) -> wire::Outcome {
    let request = wire::RequestEnvelope {
        request_id: 7,
        request,
    };
    let bytes = wire::encode_request(&request).unwrap();
    let reply = host::dispatch_frame(db, access, &bytes, scope).unwrap();
    let reply = wire::decode_response_for(&reply, &request).unwrap();
    reply
        .validate_for_in_context(
            &request,
            &wire::LabId::new(access.lab_id().as_str()).unwrap(),
        )
        .unwrap();
    reply.outcome
}
fn success(outcome: wire::Outcome) -> wire::Success {
    match outcome {
        wire::Outcome::Success(s) => s,
        other => panic!("{other:?}"),
    }
}
fn failed(outcome: wire::Outcome, code: wire::ErrorCode) -> wire::Failure {
    match outcome {
        wire::Outcome::Failure(f) => {
            assert_eq!(f.code, code);
            f
        }
        other => panic!("{other:?}"),
    }
}
fn import(
    db: &mut store::Library,
    access: &store::LabAccess,
    payload: wire::Payload,
    scope: Option<&PublicationScope>,
) -> wire::Document {
    match success(call(
        db,
        access,
        wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("Plasmid").unwrap(),
            payload,
        },
        scope,
    )) {
        wire::Success::Imported(d) => d,
        other => panic!("{other:?}"),
    }
}
#[test]
fn binary_native_save_cas_and_pinned_open_preserve_unsaved_divergence() {
    let (mut db, access) = setup();
    let original = native("Original", "AGTCCGATGCTACGTA");
    let doc = import(&mut db, &access, original.clone(), None);
    let next = native("Edited", "AGTCCGATGCTACGTT");
    let changed = match success(call(
        &mut db,
        &access,
        wire::Request::SaveRevision {
            id: doc.id.clone(),
            base_revision: doc.head.id.clone(),
            payload: next,
        },
        None,
    )) {
        wire::Success::Saved { document, revision } => {
            assert_eq!(revision.parent, Some(doc.head.id.clone()));
            document
        }
        other => panic!("{other:?}"),
    };
    let conflict = failed(
        call(
            &mut db,
            &access,
            wire::Request::SaveRevision {
                id: doc.id.clone(),
                base_revision: doc.head.id.clone(),
                payload: native("Stale competing edit", "AGTCCGATGCTACGTC"),
            },
            None,
        ),
        wire::ErrorCode::Conflict,
    );
    assert!(
        matches!(conflict.conflict,Some(wire::Conflict::Document{current}) if current.head.id==changed.head.id)
    );
    match success(call(
        &mut db,
        &access,
        wire::Request::ListRevisions {
            id: doc.id.clone(),
            page: Default::default(),
        },
        None,
    )) {
        wire::Success::Revisions(page) => assert_eq!(page.revisions.len(), 2),
        other => panic!("{other:?}"),
    }
    match success(call(
        &mut db,
        &access,
        wire::Request::OpenDocument {
            id: doc.id,
            revision: Some(doc.head.id),
        },
        None,
    )) {
        wire::Success::Opened {
            payload, revision, ..
        } => {
            assert_eq!(payload, original);
            assert_eq!(revision.number, 1);
            assert_eq!(
                protocol::load_document(payload.bytes.as_slice())
                    .unwrap()
                    .sequence
                    .as_str(),
                "AGTCCGATGCTACGTA"
            );
        }
        other => panic!("{other:?}"),
    }
}
#[test]
fn folders_metadata_conflicts_search_and_lab_context_are_native_scoped() {
    let (mut db, access) = setup();
    let folder = match success(call(
        &mut db,
        &access,
        wire::Request::CreateFolder {
            parent: None,
            name: wire::Text::new("Shared plasmids").unwrap(),
        },
        None,
    )) {
        wire::Success::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let moved = match success(call(
        &mut db,
        &access,
        wire::Request::RenameFolder {
            id: folder.id.clone(),
            expected_generation: folder.generation,
            name: wire::Text::new("Validated constructs").unwrap(),
        },
        None,
    )) {
        wire::Success::RenamedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let f = failed(
        call(
            &mut db,
            &access,
            wire::Request::DeleteEmptyFolder {
                id: folder.id.clone(),
                expected_generation: folder.generation,
            },
            None,
        ),
        wire::ErrorCode::Conflict,
    );
    assert!(
        matches!(f.conflict,Some(wire::Conflict::Folder{current}) if current.generation==moved.generation)
    );
    let doc = import(&mut db, &access, native("DNA", "ACGTACGT"), None);
    let changed = match success(call(
        &mut db,
        &access,
        wire::Request::MoveDocument {
            id: doc.id.clone(),
            expected_generation: doc.generation,
            folder: Some(folder.id.clone()),
        },
        None,
    )) {
        wire::Success::MovedDocument(d) => d,
        other => panic!("{other:?}"),
    };
    assert_eq!(changed.folder, Some(folder.id.clone()));
    failed(
        call(
            &mut db,
            &access,
            wire::Request::DeleteEmptyFolder {
                id: folder.id.clone(),
                expected_generation: moved.generation,
            },
            None,
        ),
        wire::ErrorCode::NotEmpty,
    );
    match success(call(
        &mut db,
        &access,
        wire::Request::Search {
            query: wire::Text::new("Plasmid").unwrap(),
            scope: wire::SearchScope::Folder {
                folder: Some(folder.id.clone()),
                recursive: true,
            },
            page: Default::default(),
        },
        None,
    )) {
        wire::Success::Searched(page) => assert_eq!(page.entries.len(), 1),
        other => panic!("{other:?}"),
    }
    let other = db.create_lab("Other lab").unwrap();
    let other = db.trusted_host_access(&other.id).unwrap();
    failed(
        call(
            &mut db,
            &other,
            wire::Request::OpenDocument {
                id: doc.id,
                revision: None,
            },
            None,
        ),
        wire::ErrorCode::NotFound,
    );
    failed(
        call(
            &mut db,
            &other,
            wire::Request::Browse {
                parent: Some(folder.id),
                page: Default::default(),
            },
            None,
        ),
        wire::ErrorCode::NotFound,
    );
}
#[test]
fn correct_byte_hash_cannot_bypass_native_loader_or_exact_archive_selection() {
    let (mut db, access) = setup();
    let mut invalid = native("DNA", "ACGT");
    let bytes = b"not a native document".to_vec();
    invalid.archive_hash = wire::Digest::of(&bytes);
    invalid.bytes = wire::Blob::new(bytes).unwrap();
    let failure = failed(
        call(
            &mut db,
            &access,
            wire::Request::ImportDocument {
                folder: None,
                name: wire::Text::new("Bad").unwrap(),
                payload: invalid,
            },
            None,
        ),
        wire::ErrorCode::Invalid,
    );
    assert!(!failure.message.contains("SQL"));
    let frame = Workspace::from_documents(
        "Exact document",
        vec![Molecule::dna("DNA", "ACGTACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let good = archive(&frame, 0, wire::WorkspaceSharing::ExplicitWholeWorkspace);
    let mut wrong = good.clone();
    wrong.selection.as_mut().unwrap().document_id = wire::NativeId::new("absent-document").unwrap();
    assert!(host::review_payload(&wrong).is_err());
    let mut wrong = good.clone();
    wrong.selection.as_mut().unwrap().revision = wire::Digest::of(b"another revision");
    assert!(host::review_payload(&wrong).is_err());
    let mut wrong = good.clone();
    wrong.selection.as_mut().unwrap().head_event = None;
    assert!(host::review_payload(&wrong).is_err());
    let mut wrong = good.clone();
    wrong.selection.as_mut().unwrap().head_event =
        Some(wire::NativeId::new("event-forged").unwrap());
    assert!(host::review_payload(&wrong).is_err());
    failed(
        call(
            &mut db,
            &access,
            wire::Request::ImportDocument {
                folder: None,
                name: wire::Text::new("No native scope").unwrap(),
                payload: good,
            },
            None,
        ),
        wire::ErrorCode::Forbidden,
    );
    let page = match success(call(
        &mut db,
        &access,
        wire::Request::Browse {
            parent: None,
            page: Default::default(),
        },
        None,
    )) {
        wire::Success::Browsed(p) => p,
        other => panic!("{other:?}"),
    };
    assert!(page.entries.is_empty());
    assert_eq!(page.catalog_generation, 0);
}
#[test]
fn closure_scope_rejects_unpermitted_extra_data_and_binds_exact_lab_bytes_and_sharing() {
    let (mut db, access) = setup();
    let frame = Workspace::from_documents(
        "Only shared plasmid",
        vec![Molecule::dna("DNA", "ACGTACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let permitted = ArchiveInventory {
        documents: frame.documents.iter().map(|d| d.id.clone()).collect(),
        revisions: frame.revisions.keys().cloned().collect(),
        events: frame.history.iter().map(|e| e.id.clone()).collect(),
        candidates: BTreeSet::new(),
        attachments: BTreeSet::new(),
    };
    let payload = archive(&frame, 0, wire::WorkspaceSharing::DocumentClosure);
    let review = host::review_payload(&payload).unwrap();
    let scope = PublicationScope::document_closure(&access, &review, &permitted).unwrap();
    assert!(PublicationScope::explicit_whole_workspace(&access, &review).is_err());
    let mut broader = frame.clone();
    let bytes = b"unrelated private note".to_vec();
    broader.attachments.insert(
        wire::Digest::of(&bytes).as_str().into(),
        workspace::WorkspaceAttachment {
            kind: "private-note".into(),
            name: "Unrelated".into(),
            media_type: "text/plain".into(),
            bytes,
        },
    );
    let broader = archive(&broader, 0, wire::WorkspaceSharing::DocumentClosure);
    let broader_review = host::review_payload(&broader).unwrap();
    assert!(PublicationScope::document_closure(&access, &broader_review, &permitted).is_err());
    failed(
        call(
            &mut db,
            &access,
            wire::Request::ImportDocument {
                folder: None,
                name: wire::Text::new("Broader").unwrap(),
                payload: broader,
            },
            Some(&scope),
        ),
        wire::ErrorCode::Forbidden,
    );
    let other = db.create_lab("Other").unwrap();
    let other = db.trusted_host_access(&other.id).unwrap();
    failed(
        call(
            &mut db,
            &other,
            wire::Request::ImportDocument {
                folder: None,
                name: wire::Text::new("Wrong lab").unwrap(),
                payload: payload.clone(),
            },
            Some(&scope),
        ),
        wire::ErrorCode::Forbidden,
    );
    let doc = import(&mut db, &access, payload.clone(), Some(&scope));
    match success(call(
        &mut db,
        &access,
        wire::Request::OpenDocument {
            id: doc.id,
            revision: None,
        },
        None,
    )) {
        wire::Success::Opened {
            payload: opened, ..
        } => assert_eq!(opened, payload),
        other => panic!("{other:?}"),
    }
}
#[test]
fn ordinary_save_scope_is_derived_and_rejects_a_whole_archive_labeled_as_one_plasmid() {
    let (mut db, access) = setup();
    assert!(
        host::document_publication_scope(&access, &native("Bare DNA", "ACGT"))
            .unwrap()
            .is_none()
    );
    let frame = Workspace::from_documents(
        "Local workspace",
        vec![
            Molecule::dna("Shared plasmid", "AGTCCGAT", Topology::Linear).unwrap(),
            Molecule::dna(
                "Unpublished unrelated plasmid",
                "TCAGGCTA",
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let selected = &frame.documents[0];
    let selection = collagene_library_sharing::Selection {
        document_id: selected.id.clone(),
        revision: selected.head.clone(),
        head_event: selected.head_event.clone(),
    };
    let full = archive(&frame, 0, wire::WorkspaceSharing::DocumentClosure);
    assert!(host::document_publication_scope(&access, &full).is_err());
    let export = collagene_library_sharing::document_closure(&frame, &selection).unwrap();
    let at = export
        .workspace
        .documents
        .iter()
        .position(|d| d.id == selected.id)
        .unwrap();
    let payload = archive(
        &export.workspace,
        at,
        wire::WorkspaceSharing::DocumentClosure,
    );
    let scope = host::document_publication_scope(&access, &payload)
        .unwrap()
        .unwrap();
    import(&mut db, &access, payload, Some(&scope));
    let whole = archive(&frame, 0, wire::WorkspaceSharing::ExplicitWholeWorkspace);
    assert!(host::document_publication_scope(&access, &whole).is_err());
    let review = host::review_payload(&whole).unwrap();
    let explicit_scope = PublicationScope::explicit_whole_workspace(&access, &review).unwrap();
    let mut other_db = store::Library::open_in_memory().unwrap();
    let lab = other_db.create_lab("Explicit whole").unwrap();
    let other_access = other_db.trusted_host_access(&lab.id).unwrap();
    // Explicit sharing is still bound to its original lab capability.
    failed(
        call(
            &mut other_db,
            &other_access,
            wire::Request::ImportDocument {
                folder: None,
                name: wire::Text::new("Wrong scope").unwrap(),
                payload: whole.clone(),
            },
            Some(&explicit_scope),
        ),
        wire::ErrorCode::Forbidden,
    );
    let explicit_scope =
        PublicationScope::explicit_whole_workspace(&other_access, &review).unwrap();
    import(&mut other_db, &other_access, whole, Some(&explicit_scope));
}

fn event(frame: &Workspace, index: usize) -> EventBinding {
    let d = &frame.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn real_saved_run_archive_retains_ancestry_while_native_snapshot_never_restores_material() {
    let source = Molecule::dna(
        "Original stock",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Separate required product",
        source.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let frame =
        Workspace::from_documents("PCR origin", vec![source.clone(), target.clone()], 0).unwrap();
    let plan = DesignPlan {
        name: "Actual source PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let saved = collagene_core::cloning_outputs::finalize(
        &frame,
        &plan,
        event(&frame, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&frame, 0),
        }],
        &collagene_core::cloning_outputs::Options {
            name: "Saved PCR stock".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selection = collagene_core::run_material::Selection {
        run_id: saved.receipt.run.run_id.clone(),
        step_index: 0,
        event: saved.receipt.construct.clone(),
    };
    let closure_selection = collagene_library_sharing::Selection {
        document_id: selection.event.document_id.clone(),
        revision: selection.event.revision.clone(),
        head_event: Some(selection.event.event_id.clone()),
    };
    let export =
        collagene_library_sharing::document_closure(&saved.workspace, &closure_selection).unwrap();
    let index = export
        .workspace
        .documents
        .iter()
        .position(|d| d.id == selection.event.document_id)
        .unwrap();
    let payload = archive(
        &export.workspace,
        index,
        wire::WorkspaceSharing::DocumentClosure,
    );
    let (mut db, access) = setup();
    let scope = host::document_publication_scope(&access, &payload)
        .unwrap()
        .unwrap();
    let doc = import(&mut db, &access, payload, Some(&scope));
    let opened = match success(call(
        &mut db,
        &access,
        wire::Request::OpenDocument {
            id: doc.id,
            revision: None,
        },
        None,
    )) {
        wire::Success::Opened { payload, .. } => payload,
        other => panic!("{other:?}"),
    };
    let loaded = workspace::load_workspace(opened.bytes.as_slice()).unwrap();
    assert_eq!(loaded.history.len(), export.workspace.history.len());
    assert_eq!(loaded.attachments.len(), export.workspace.attachments.len());
    let fresh = collagene_core::run_material::reconstruct(&loaded, &selection).unwrap();
    assert_eq!(fresh.report().status, Status::ConsistentUnderAssumptions);
    let molecule = &loaded.revisions[&selection.event.revision];
    let bytes = protocol::save_document(molecule).unwrap();
    let isolated = protocol::load_document(&bytes).unwrap();
    assert_eq!(isolated.revision(), molecule.revision());
    let isolated = Workspace::from_documents("Reference snapshot only", vec![isolated], 0).unwrap();
    assert!(collagene_core::run_material::reconstruct(&isolated, &selection).is_err());
}
#[test]
fn malformed_frames_and_unknown_actor_claims_never_reach_storage() {
    let (mut db, access) = setup();
    let request = wire::RequestEnvelope {
        request_id: 1,
        request: wire::Request::Browse {
            parent: None,
            page: Default::default(),
        },
    };
    let mut frame = wire::encode_request(&request).unwrap();
    frame.push(0);
    assert!(host::dispatch_frame(&mut db, &access, &frame, None).is_err());
    let claimed = serde_json::json!({"request_id":1,"request":{"Browse":{"parent":null,"page":{"offset":0,"limit":20,"expected_catalog_generation":null},"actor":"administrator","lab":"other"}}});
    assert!(serde_json::from_value::<wire::RequestEnvelope>(claimed).is_err());
    let page = db.browse(&access, None, 0, 20, None).unwrap();
    assert_eq!(page.catalog_generation, 0);
    assert!(page.entries.is_empty());
}
