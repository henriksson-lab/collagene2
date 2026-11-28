use collagene_core::{
    documents::workspace::{self, Workspace},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_desktop_host::{DesktopHost, WindowAccess};
use collagene_desktop_protocol as wire;
use collagene_model::{Molecule, Topology};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
use wire::library::{self as library, Blob, Digest};
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "collagene-desktop-host-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        Self(root)
    }
}
impl Drop for Disk {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn operation(n: u64) -> wire::OperationId {
    wire::OperationId::new(format!("op-{n:032x}")).unwrap()
}
fn request(n: u64, command: wire::Command) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: n,
        operation_id: operation(n),
        command,
    }
}
fn call(
    host: &mut DesktopHost,
    window: &WindowAccess,
    n: u64,
    command: wire::Command,
) -> wire::Outcome {
    let request = request(n, command);
    let frame = wire::encode_request(&request).unwrap();
    let reply = host.dispatch_raw(window, &frame).unwrap();
    wire::decode_response_for(&reply, &request).unwrap().outcome
}
fn success(outcome: wire::Outcome) -> wire::Success {
    match outcome {
        wire::Outcome::Success(s) => s,
        other => panic!("{other:?}"),
    }
}
fn metadata(host: &mut DesktopHost, w: &WindowAccess, n: u64) -> wire::WorkingMetadata {
    match success(call(host, w, n, wire::Command::WorkingMetadata)) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    }
}
fn archive(name: &str) -> Vec<u8> {
    let molecule = Molecule::dna(
        name,
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let mut frame =
        Workspace::from_documents("Complete private working state", vec![molecule], 0).unwrap();
    let bytes = format!("{{\"private unfinished draft\":\"{name}\"}}").into_bytes();
    frame.attachments.insert(
        Digest::of(&bytes).as_str().into(),
        workspace::WorkspaceAttachment {
            kind: "assembly-draft-v10".into(),
            name: "Private draft".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    workspace::save_workspace(&frame).unwrap()
}
fn read(
    host: &mut DesktopHost,
    w: &WindowAccess,
    n: u64,
    slot: wire::WorkingSlot,
    version: &Option<Digest>,
) -> Vec<u8> {
    match success(call(
        host,
        w,
        n,
        wire::Command::ReadWorking {
            slot,
            expected_version: version.clone(),
        },
    )) {
        wire::Success::WorkingArchive { archive, .. } => archive.into_vec(),
        other => panic!("{other:?}"),
    }
}
#[test]
fn current_caller_backup_preferences_restart_and_stale_windows_use_independent_cas() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let first = host.authorize_window().unwrap();
    let second = host.authorize_window().unwrap();
    let bootstrap = match success(call(&mut host, &first, 1, wire::Command::Bootstrap)) {
        wire::Success::Bootstrap(b) => b,
        other => panic!("{other:?}"),
    };
    assert!(bootstrap.working.version.is_none());
    assert!(bootstrap.preferences.version.is_none());
    let a = archive("On-disk current A");
    let b = archive("Replacement B");
    let unsaved = archive("Unsaved in-memory caller backup");
    let c = archive("Current-only C");
    let saved = match success(call(
        &mut host,
        &first,
        2,
        wire::Command::SaveCurrent {
            expected_version: None,
            archive: Blob::new(a).unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    let replaced = match success(call(
        &mut host,
        &first,
        3,
        wire::Command::ReplaceWorking {
            expected_version: saved.version.clone(),
            current: Blob::new(b).unwrap(),
            previous: Blob::new(unsaved.clone()).unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    let current = match success(call(
        &mut host,
        &first,
        4,
        wire::Command::SaveCurrent {
            expected_version: replaced.version.clone(),
            archive: Blob::new(c.clone()).unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    assert_eq!(
        read(
            &mut host,
            &first,
            5,
            wire::WorkingSlot::Previous,
            &current.version
        ),
        unsaved
    );
    assert!(
        matches!(call(&mut host,&second,6,wire::Command::SaveCurrent{expected_version:saved.version,archive:Blob::new(archive("stale write")).unwrap()}),wire::Outcome::Failure(wire::Failure::Conflict{current:wire::Conflict::Working(m)}) if m.version==current.version)
    );
    let prefs = match success(call(
        &mut host,
        &first,
        7,
        wire::Command::SavePreferences {
            expected_version: None,
            bytes: Blob::new(br#"{"zoom":1.5}"#.to_vec()).unwrap(),
        },
    )) {
        wire::Success::Preferences(p) => p,
        other => panic!("{other:?}"),
    };
    assert!(
        matches!(call(&mut host,&second,8,wire::Command::SavePreferences{expected_version:None,bytes:Blob::new(vec![]).unwrap()}),wire::Outcome::Failure(wire::Failure::Conflict{current:wire::Conflict::Preferences(p)}) if p.version==prefs.version)
    );
    assert_eq!(metadata(&mut host, &first, 9).version, current.version);
    drop(host);
    let mut reopened = DesktopHost::open(&disk.0).unwrap();
    let w = reopened.authorize_window().unwrap();
    let state = metadata(&mut reopened, &w, 10);
    assert_eq!(state, current);
    assert_eq!(
        read(
            &mut reopened,
            &w,
            11,
            wire::WorkingSlot::Current,
            &state.version
        ),
        c
    );
    assert_eq!(
        read(
            &mut reopened,
            &w,
            12,
            wire::WorkingSlot::Previous,
            &state.version
        ),
        unsaved
    );
    match success(call(&mut reopened, &w, 13, wire::Command::ReadPreferences)) {
        wire::Success::Preferences(p) => assert_eq!(p, prefs),
        other => panic!("{other:?}"),
    }
    let same = match success(call(
        &mut reopened,
        &w,
        14,
        wire::Command::SaveCurrent {
            expected_version: state.version.clone(),
            archive: Blob::new(c).unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    assert_ne!(same.version, state.version);
}
#[test]
fn invalid_native_pair_is_atomic_and_corrupt_current_cannot_replace_good_recovery() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    let current = archive("Current");
    let previous = archive("Known-good previous");
    let state = match success(call(
        &mut host,
        &w,
        1,
        wire::Command::ReplaceWorking {
            expected_version: None,
            current: Blob::new(current.clone()).unwrap(),
            previous: Blob::new(previous.clone()).unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    assert!(matches!(
        call(
            &mut host,
            &w,
            2,
            wire::Command::ReplaceWorking {
                expected_version: state.version.clone(),
                current: Blob::new(archive("Must roll back")).unwrap(),
                previous: Blob::new(b"not an archive".to_vec()).unwrap()
            }
        ),
        wire::Outcome::Failure(wire::Failure::Invalid { .. })
    ));
    assert_eq!(metadata(&mut host, &w, 3), state);
    let db = rusqlite::Connection::open(disk.0.join("working.sqlite")).unwrap();
    let mut corrupt = current;
    *corrupt.last_mut().unwrap() ^= 1;
    db.execute(
        "UPDATE working SET current_archive=?1 WHERE singleton=1",
        [&corrupt],
    )
    .unwrap();
    drop(db);
    assert!(matches!(
        call(
            &mut host,
            &w,
            4,
            wire::Command::ReadWorking {
                slot: wire::WorkingSlot::Current,
                expected_version: state.version.clone()
            }
        ),
        wire::Outcome::Failure(wire::Failure::Invalid { .. })
    ));
    assert_eq!(
        read(
            &mut host,
            &w,
            5,
            wire::WorkingSlot::Previous,
            &state.version
        ),
        previous
    );
    let recovered = match success(call(
        &mut host,
        &w,
        6,
        wire::Command::RecoverPrevious {
            expected_version: state.version.unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    assert_eq!(
        read(
            &mut host,
            &w,
            7,
            wire::WorkingSlot::Current,
            &recovered.version
        ),
        previous
    );
    assert_eq!(
        read(
            &mut host,
            &w,
            8,
            wire::WorkingSlot::Previous,
            &recovered.version
        ),
        previous
    );
    drop(host);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    assert_eq!(metadata(&mut host, &w, 9), recovered);
}
#[test]
fn native_window_scope_cancel_before_dispatch_and_control_lane_do_not_mutate_storage() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let a = host.authorize_window().unwrap();
    let b = host.authorize_window().unwrap();
    let control = host.control();
    let target = operation(100);
    let cancel = request(
        1,
        wire::Command::Cancel {
            operation_id: target.clone(),
        },
    );
    let bytes = wire::encode_request(&cancel).unwrap();
    let reply = control.dispatch_raw(&a, &bytes).unwrap().unwrap();
    assert!(matches!(
        wire::decode_response_for(&reply, &cancel).unwrap().outcome,
        wire::Outcome::Success(wire::Success::Cancelled {
            status: wire::OperationState::Cancelled,
            ..
        })
    ));
    let command = wire::Command::SaveCurrent {
        expected_version: None,
        archive: Blob::new(archive("cancelled")).unwrap(),
    };
    assert!(matches!(
        call(&mut host, &a, 100, command.clone()),
        wire::Outcome::Failure(wire::Failure::Cancelled)
    ));
    assert!(metadata(&mut host, &a, 2).version.is_none());
    assert!(matches!(
        call(&mut host, &b, 100, command),
        wire::Outcome::Success(_)
    ));
    let big = vec![0u8; wire::MAX_CONTROL_BYTES + 1];
    assert!(control.dispatch_raw(&a, &big).unwrap().is_none());
    control.revoke_window(&a).unwrap();
    assert!(
        host.dispatch_raw(
            &a,
            &wire::encode_request(&request(3, wire::Command::Bootstrap)).unwrap()
        )
        .is_err()
    );
    let other = Disk::new();
    let mut other = DesktopHost::open(&other.0).unwrap();
    assert!(
        other
            .dispatch_raw(
                &b,
                &wire::encode_request(&request(3, wire::Command::Bootstrap)).unwrap()
            )
            .is_err()
    );
}
#[test]
fn local_library_mutations_use_structured_native_rejections_and_private_workspace_loader() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    let name = library::Text::new("Local plasmids").unwrap();
    let folder = match success(call(
        &mut host,
        &w,
        1,
        wire::Command::Library(library::Request::CreateFolder {
            parent: None,
            name: name.clone(),
        }),
    )) {
        wire::Success::Library(library::Success::CreatedFolder(f)) => f,
        other => panic!("{other:?}"),
    };
    assert!(matches!(
        call(
            &mut host,
            &w,
            2,
            wire::Command::Library(library::Request::CreateFolder { parent: None, name })
        ),
        wire::Outcome::Failure(wire::Failure::Library(library::Failure {
            code: library::ErrorCode::AlreadyExists,
            ..
        }))
    ));
    let renamed = match success(call(
        &mut host,
        &w,
        3,
        wire::Command::Library(library::Request::RenameFolder {
            id: folder.id.clone(),
            expected_generation: folder.generation,
            name: library::Text::new("Edited locally").unwrap(),
        }),
    )) {
        wire::Success::Library(library::Success::RenamedFolder(f)) => f,
        other => panic!("{other:?}"),
    };
    assert!(matches!(
        call(
            &mut host,
            &w,
            4,
            wire::Command::Library(library::Request::RenameFolder {
                id: folder.id,
                expected_generation: folder.generation,
                name: library::Text::new("Stale").unwrap()
            })
        ),
        wire::Outcome::Failure(wire::Failure::Library(library::Failure {
            code: library::ErrorCode::Conflict,
            ..
        }))
    ));
    assert!(renamed.generation > folder.generation);
    let body = archive("Full private local library archive");
    let loaded = workspace::load_workspace(&body).unwrap();
    let selected = &loaded.documents[0];
    let payload = library::Payload {
        kind: library::PayloadKind::WorkspaceArchive,
        archive_hash: Digest::of(&body),
        bytes: Blob::new(body.clone()).unwrap(),
        selection: Some(library::NativeSelection {
            document_id: library::NativeId::new(&selected.id).unwrap(),
            revision: Digest::new(&selected.head).unwrap(),
            head_event: selected
                .head_event
                .as_ref()
                .map(|id| library::NativeId::new(id).unwrap()),
        }),
        sharing: Some(library::WorkspaceSharing::ExplicitWholeWorkspace),
    };
    let doc = match success(call(
        &mut host,
        &w,
        5,
        wire::Command::Library(library::Request::ImportDocument {
            folder: None,
            name: library::Text::new("Private design").unwrap(),
            payload,
        }),
    )) {
        wire::Success::Library(library::Success::Imported(d)) => d,
        other => panic!("{other:?}"),
    };
    match success(call(
        &mut host,
        &w,
        6,
        wire::Command::Library(library::Request::OpenDocument {
            id: doc.id,
            revision: None,
        }),
    )) {
        wire::Success::Library(library::Success::Opened { payload, .. }) => {
            assert_eq!(payload.bytes.as_slice(), body)
        }
        other => panic!("{other:?}"),
    }
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
fn native_run_history_and_fromrun_material_survive_working_restart() {
    let disk = Disk::new();
    let source = Molecule::dna(
        "Source",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent target",
        source.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let frame = Workspace::from_documents(
        "Private source events",
        vec![source.clone(), target.clone()],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let recorded = run_lineage::record_run(
        &frame,
        &plan,
        event(&frame, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&frame, 0),
        }],
        false,
    )
    .unwrap();
    let product = &recorded.receipt.products[0];
    let origin = run_material::Selection {
        run_id: recorded.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let bytes = workspace::save_workspace(&recorded.workspace).unwrap();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    let state = match success(call(
        &mut host,
        &w,
        1,
        wire::Command::SaveCurrent {
            expected_version: None,
            archive: Blob::new(bytes).unwrap(),
        },
    )) {
        wire::Success::WorkingMetadata(m) => m,
        other => panic!("{other:?}"),
    };
    drop(host);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    let body = read(&mut host, &w, 2, wire::WorkingSlot::Current, &state.version);
    let loaded = workspace::load_workspace(&body).unwrap();
    assert_eq!(loaded.history.len(), recorded.workspace.history.len());
    assert_eq!(
        run_material::reconstruct(&loaded, &origin)
            .unwrap()
            .report()
            .status,
        Status::ConsistentUnderAssumptions
    );
}

#[cfg(unix)]
#[test]
fn unix_new_roots_and_sqlite_files_are_private_and_insecure_existing_storage_is_rejected() {
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    let disk = Disk::new();
    let host = DesktopHost::open(&disk.0).unwrap();
    let root = std::fs::metadata(&disk.0).unwrap();
    assert_eq!(root.permissions().mode() & 0o777, 0o700);
    assert_eq!(root.uid(), unsafe { libc::geteuid() });
    for name in [
        "repository.sqlite",
        "working.sqlite",
        "repository.sqlite-wal",
        "repository.sqlite-shm",
        "working.sqlite-wal",
        "working.sqlite-shm",
    ] {
        if let Ok(meta) = std::fs::metadata(disk.0.join(name)) {
            assert_eq!(meta.permissions().mode() & 0o777, 0o600, "{name}");
            assert_eq!(meta.uid(), unsafe { libc::geteuid() });
        }
    }
    drop(host);
    std::fs::set_permissions(&disk.0, std::fs::Permissions::from_mode(0o755)).unwrap();
    assert!(DesktopHost::open(&disk.0).is_err());
    std::fs::set_permissions(&disk.0, std::fs::Permissions::from_mode(0o700)).unwrap();
    let file = disk.0.join("repository.sqlite");
    std::fs::set_permissions(&file, std::fs::Permissions::from_mode(0o644)).unwrap();
    let paused = DesktopHost::open(&disk.0)
        .expect("unsafe catalog permissions must leave private working recovery accessible");
    assert!(
        !paused.catalog_available(),
        "unsafe catalog is never opened as authoritative material"
    );
    drop(paused);
    std::fs::set_permissions(&file, std::fs::Permissions::from_mode(0o600)).unwrap();
    assert!(DesktopHost::open(&disk.0).is_ok());
}
#[test]
fn native_mutation_replay_after_restart_retains_original_outcome_and_rejects_intent_reuse() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    let molecule = Molecule::dna("Native import", "ACGTACGT", Topology::Linear).unwrap();
    let bytes = collagene_core::protocol::save_document(&molecule).unwrap();
    let payload = library::Payload {
        kind: library::PayloadKind::NativeDocument,
        archive_hash: Digest::of(&bytes),
        bytes: Blob::new(bytes).unwrap(),
        selection: None,
        sharing: None,
    };
    let original = wire::Command::Library(library::Request::ImportDocument {
        folder: None,
        name: library::Text::new("Original free name").unwrap(),
        payload: payload.clone(),
    });
    let first = match success(call(&mut host, &w, 500, original.clone())) {
        wire::Success::Library(library::Success::Imported(d)) => d,
        other => panic!("{other:?}"),
    };
    success(call(
        &mut host,
        &w,
        501,
        wire::Command::Library(library::Request::RenameDocument {
            id: first.id.clone(),
            expected_generation: first.generation,
            name: library::Text::new("Renamed after original commit").unwrap(),
        }),
    ));
    drop(host);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    match success(call(&mut host, &w, 500, original)) {
        wire::Success::Library(library::Success::Imported(d)) => assert_eq!(d, first),
        other => panic!("{other:?}"),
    }
    match success(call(
        &mut host,
        &w,
        502,
        wire::Command::Library(library::Request::Browse {
            parent: None,
            page: Default::default(),
        }),
    )) {
        wire::Success::Library(library::Success::Browsed(page)) => {
            assert_eq!(page.entries.len(), 1);
            assert!(
                matches!(&page.entries[0],library::Entry::Document(d) if d.id==first.id&&d.name.as_str()=="Renamed after original commit")
            );
        }
        other => panic!("{other:?}"),
    }
    drop(host);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let w = host.authorize_window().unwrap();
    let changed = wire::Command::Library(library::Request::ImportDocument {
        folder: None,
        name: library::Text::new("Different intent").unwrap(),
        payload,
    });
    assert!(matches!(
        call(&mut host, &w, 500, changed),
        wire::Outcome::Failure(wire::Failure::Invalid { .. })
    ));
}
