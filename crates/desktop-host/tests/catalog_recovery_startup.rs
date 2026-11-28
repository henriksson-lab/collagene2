//! Real temporary native state; no shell, network, browser or credential service.
use collagene_core::documents::workspace::{self, Workspace};
use collagene_desktop_host::{DesktopHost, WindowAccess, catalog_backup::DesktopCatalogBackups};
use collagene_desktop_protocol::{self as w, library as b};
use collagene_model::{Molecule, Topology};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self(std::env::temp_dir().join(format!(
                "collagene-editor-recovery-{}-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            )))
    }
}
impl Drop for Disk {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn op(n: u64) -> w::OperationId {
    w::OperationId::new(format!("op-{n:032x}")).unwrap()
}
fn request(n: u64, c: w::Command) -> w::RequestEnvelope {
    w::RequestEnvelope {
        request_id: n,
        operation_id: op(n),
        command: c,
    }
}
fn host_call(host: &mut DesktopHost, access: &WindowAccess, n: u64, c: w::Command) -> w::Outcome {
    let q = request(n, c);
    let bytes = w::encode_request(&q).unwrap();
    let reply = host.dispatch_raw(access, &bytes).unwrap();
    w::decode_response_for(&reply, &q).unwrap().outcome
}
fn backup_call(
    backups: &DesktopCatalogBackups,
    access: &WindowAccess,
    n: u64,
    c: w::CatalogBackupCommand,
) -> w::Outcome {
    let q = request(n, w::Command::CatalogBackup(c));
    let bytes = w::encode_request(&q).unwrap();
    let reply = backups.dispatch_raw(access, &bytes).unwrap();
    w::decode_response_for(&reply, &q).unwrap().outcome
}
fn good(o: w::Outcome) -> w::Success {
    match o {
        w::Outcome::Success(s) => s,
        x => panic!("unexpected failure {x:?}"),
    }
}
fn bootstrap(host: &mut DesktopHost, access: &WindowAccess, n: u64) -> w::Bootstrap {
    match good(host_call(host, access, n, w::Command::Bootstrap)) {
        w::Success::Bootstrap(b) => b,
        x => panic!("bootstrap {x:?}"),
    }
}
fn archive(name: &str) -> Vec<u8> {
    let m = Molecule::dna(
        name,
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let mut frame = Workspace::from_documents("Private current editor", vec![m], 0).unwrap();
    let bytes = format!("{{\"unfinished\":\"{name}\"}}").into_bytes();
    frame.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        workspace::WorkspaceAttachment {
            kind: "private-recovery-fixture-v1".into(),
            name: "Unpublished editor draft".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    workspace::save_workspace(&frame).unwrap()
}
struct Seed {
    disk: Disk,
    current: Vec<u8>,
    previous: Vec<u8>,
    meta: w::WorkingMetadata,
    prefs: w::Preferences,
    accepted: w::CatalogActivationReceipt,
}
fn seed() -> Seed {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let current = archive("Latest private editor");
    let previous = archive("Unsaved in-memory backup");
    let meta = match good(host_call(
        &mut host,
        &access,
        1,
        w::Command::ReplaceWorking {
            expected_version: None,
            current: b::Blob::new(current.clone()).unwrap(),
            previous: b::Blob::new(previous.clone()).unwrap(),
        },
    )) {
        w::Success::WorkingMetadata(m) => m,
        _ => panic!("working"),
    };
    let prefs = match good(host_call(
        &mut host,
        &access,
        2,
        w::Command::SavePreferences {
            expected_version: None,
            bytes: b::Blob::new(b"{\"zoom\":3}".to_vec()).unwrap(),
        },
    )) {
        w::Success::Preferences(p) => p,
        _ => panic!("prefs"),
    };
    let backups = host.open_catalog_backup().unwrap();
    let backup = match good(backup_call(
        &backups,
        &access,
        10,
        w::CatalogBackupCommand::Create,
    )) {
        w::Success::CatalogBackup(w::CatalogBackupSuccess::Created(r)) => r,
        _ => panic!("backup"),
    };
    good(backup_call(
        &backups,
        &access,
        11,
        w::CatalogBackupCommand::RestoreCopy {
            id: backup.id.clone(),
            expected_hash: backup.archive_hash.clone(),
        },
    ));
    let review = match good(backup_call(
        &backups,
        &access,
        12,
        w::CatalogBackupCommand::ReviewRecovery {
            id: op(11),
            expected_hash: backup.archive_hash,
        },
    )) {
        w::Success::CatalogBackup(w::CatalogBackupSuccess::Reviewed(r)) => r,
        _ => panic!("review"),
    };
    let accepted = match good(backup_call(
        &backups,
        &access,
        13,
        w::CatalogBackupCommand::AdoptRecovery { review },
    )) {
        w::Success::CatalogBackup(w::CatalogBackupSuccess::AdoptionPending(r)) => r,
        _ => panic!("accepted"),
    };
    assert!(!host.catalog_available());
    drop(backups);
    drop(host);
    Seed {
        disk,
        current,
        previous,
        meta,
        prefs,
        accepted,
    }
}
fn assert_editor(
    seed: &Seed,
    host: &mut DesktopHost,
    access: &WindowAccess,
    n: u64,
) -> w::Bootstrap {
    let boot = bootstrap(host, access, n);
    assert_eq!(boot.working, seed.meta);
    assert_eq!(boot.preferences, seed.prefs);
    for (offset, slot, expected) in [
        (1, w::WorkingSlot::Current, &seed.current),
        (2, w::WorkingSlot::Previous, &seed.previous),
    ] {
        match good(host_call(
            host,
            access,
            n + offset,
            w::Command::ReadWorking {
                slot,
                expected_version: seed.meta.version.clone(),
            },
        )) {
            w::Success::WorkingArchive { archive, .. } => {
                assert_eq!(&*archive, expected.as_slice())
            }
            x => panic!("read {x:?}"),
        }
    }
    boot
}
fn cancel(seed: &Seed) -> w::CatalogBackupCommand {
    w::CatalogBackupCommand::CancelActivation {
        accepted_operation: seed.accepted.operation_id.clone(),
        review_revision: seed.accepted.review.review_revision.clone(),
        activation_generation: seed.accepted.review.activation_generation,
    }
}
#[test]
fn changed_source_opens_editor_only_and_exact_cancel_survives_retry_restart_without_rewriting_branches()
 {
    let seed = seed();
    let mut old =
        collagene_repository::Repository::open(seed.disk.0.join("repository.sqlite")).unwrap();
    old.mutate_local(
        collagene_lab_library::mutations::MutationOperation::CreateFolder {
            parent: None,
            name: "Late private original branch".into(),
        },
    )
    .unwrap();
    drop(old);
    let mut host = DesktopHost::open(&seed.disk.0)
        .expect("catalog failure must not block working editor startup");
    assert!(!host.catalog_available());
    let access = host.authorize_window().unwrap();
    let boot = assert_editor(&seed, &mut host, &access, 100);
    assert!(boot.library_context.is_none());
    assert!(
        matches!(boot.catalog,w::CatalogStatus::Paused{activation:Some(w::CatalogActivationStatus::PendingRestart(ref r)),..} if r==&seed.accepted)
    );
    assert!(host.open_repository_connection().is_err());
    let backups = host
        .open_catalog_backup()
        .expect("root-only recovery control is always available");
    let wrong = w::CatalogBackupCommand::CancelActivation {
        accepted_operation: seed.accepted.operation_id.clone(),
        review_revision: b::Digest::of(b"wrong exact review"),
        activation_generation: seed.accepted.review.activation_generation,
    };
    assert!(matches!(
        backup_call(&backups, &access, 110, wrong),
        w::Outcome::Failure(_)
    ));
    assert!(
        !seed
            .disk
            .0
            .join("catalog-activation/00000000000000000001.cancelled.json")
            .exists()
    );
    let cancelled = good(backup_call(&backups, &access, 111, cancel(&seed)));
    assert_eq!(
        good(backup_call(&backups, &access, 111, cancel(&seed))),
        cancelled
    );
    let w::Success::CatalogBackup(w::CatalogBackupSuccess::ActivationCancelled(receipt)) =
        &cancelled
    else {
        panic!("cancelled")
    };
    assert_eq!(receipt.accepted, seed.accepted);
    assert_eq!(receipt.generation, 1);
    let boot = assert_editor(&seed, &mut host, &access, 120);
    assert!(boot.library_context.is_none());
    let mut fabricated = boot.clone();
    fabricated.library_context = Some(b::LabContext {
        lab_id: seed.accepted.review.current_lab.clone(),
        name: b::Text::new("Fabricated available catalog").unwrap(),
    });
    assert!(
        fabricated.validate().is_err(),
        "paused bootstrap may not invent a catalog context"
    );
    assert!(matches!(
        boot.catalog,
        w::CatalogStatus::Paused {
            activation: Some(w::CatalogActivationStatus::Cancelled(_)),
            ..
        }
    ));
    assert!(!host.catalog_available());
    let mut next_review = seed.accepted.review.clone();
    next_review.activation_generation = 1;
    next_review.review_revision = next_review.revision().unwrap();
    assert!(
        matches!(
            backup_call(
                &backups,
                &access,
                112,
                w::CatalogBackupCommand::AdoptRecovery {
                    review: next_review
                }
            ),
            w::Outcome::Failure(_)
        ),
        "cancellation cannot admit another selection before restart"
    );
    let same_bootstrap = bootstrap(&mut host, &access, 100);
    assert!(
        matches!(
            same_bootstrap.catalog,
            w::CatalogStatus::Paused {
                activation: Some(w::CatalogActivationStatus::Cancelled(_)),
                ..
            }
        ),
        "cached pending bootstrap must refresh after cancellation"
    );
    assert!(matches!(
        host_call(
            &mut host,
            &access,
            130,
            w::Command::Library(b::Request::CreateFolder {
                parent: None,
                name: b::Text::new("Must remain blocked").unwrap()
            })
        ),
        w::Outcome::Failure(_)
    ));
    drop(backups);
    drop(host);
    let mut restarted = DesktopHost::open(&seed.disk.0).unwrap();
    assert!(restarted.catalog_available());
    let access = restarted.authorize_window().unwrap();
    let boot = assert_editor(&seed, &mut restarted, &access, 140);
    assert!(boot.library_context.is_some());
    assert!(matches!(
        boot.catalog,
        w::CatalogStatus::Ready {
            activation: w::CatalogActivationStatus::Cancelled(_)
        }
    ));
    let backups = restarted.open_catalog_backup().unwrap();
    assert_eq!(
        good(backup_call(&backups, &access, 111, cancel(&seed))),
        cancelled,
        "lost reply survives a process restart"
    );
    assert!(
        backups.catalog_session_ready(),
        "historical cancellation retry does not pause the current generation"
    );
    let c = rusqlite::Connection::open_with_flags(
        seed.disk.0.join("repository.sqlite"),
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )
    .unwrap();
    let count: u64 = c
        .query_row(
            "SELECT count(*) FROM folders WHERE name='Late private original branch'",
            [],
            |r| r.get(0),
        )
        .unwrap();
    assert_eq!(count, 1);
    assert!(
        seed.disk
            .0
            .join("recoveries")
            .join(op(11).as_str())
            .join("library.sqlite")
            .is_file()
    );
}
#[test]
fn corrupt_target_can_be_cancelled_without_opening_that_database_and_editor_stays_available() {
    let seed = seed();
    let target = seed
        .disk
        .0
        .join("recoveries")
        .join(op(11).as_str())
        .join("library.sqlite");
    let mut bytes = fs::read(&target).unwrap();
    bytes[0] ^= 1;
    fs::write(&target, &bytes).unwrap();
    let mut host = DesktopHost::open(&seed.disk.0).unwrap();
    assert!(!host.catalog_available());
    let access = host.authorize_window().unwrap();
    let boot = assert_editor(&seed, &mut host, &access, 200);
    assert!(boot.library_context.is_none());
    let backups = host.open_catalog_backup().unwrap();
    good(backup_call(&backups, &access, 210, cancel(&seed)));
    assert_eq!(
        fs::read(&target).unwrap(),
        bytes,
        "cancellation never repairs/rewrites a branch"
    );
    assert_editor(&seed, &mut host, &access, 220);
    drop(backups);
    drop(host);
    let restarted = DesktopHost::open(&seed.disk.0).unwrap();
    assert!(restarted.catalog_available());
    assert_eq!(fs::read(&target).unwrap(), bytes);
}
#[test]
fn corrupt_activation_journal_does_not_block_editor_but_refuses_cancellation_and_catalog() {
    let seed = seed();
    let accepted = seed
        .disk
        .0
        .join("catalog-activation/00000000000000000001.accepted.json");
    let mut bytes = fs::read(&accepted).unwrap();
    bytes[0] ^= 1;
    fs::write(&accepted, &bytes).unwrap();
    let mut host = DesktopHost::open(&seed.disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    assert!(!host.catalog_available());
    let boot = assert_editor(&seed, &mut host, &access, 300);
    assert!(boot.library_context.is_none());
    assert!(matches!(
        boot.catalog,
        w::CatalogStatus::Paused {
            activation: None,
            ..
        }
    ));
    let backups = host
        .open_catalog_backup()
        .expect("corrupt journal must not prevent root-only facade");
    assert!(matches!(
        backup_call(&backups, &access, 310, cancel(&seed)),
        w::Outcome::Failure(_)
    ));
    assert_eq!(fs::read(&accepted).unwrap(), bytes);
    assert!(
        !seed
            .disk
            .0
            .join("catalog-activation/00000000000000000001.cancelled.json")
            .exists()
    );
    let prefs = good(host_call(
        &mut host,
        &access,
        311,
        w::Command::SavePreferences {
            expected_version: seed.prefs.version.clone(),
            bytes: b::Blob::new(b"{\"editorOnly\":true}".to_vec()).unwrap(),
        },
    ));
    assert!(matches!(prefs, w::Success::Preferences(_)));
    assert!(host.open_exports().is_ok());
}
#[test]
fn bootstrap_and_exact_cancel_use_the_only_current_native_frame() {
    for command in [
        w::Command::Bootstrap,
        w::Command::CatalogBackup(w::CatalogBackupCommand::CancelActivation {
            accepted_operation: op(1),
            review_revision: b::Digest::of(b"review"),
            activation_generation: 0,
        }),
    ] {
        let q = request(400, command);
        let bytes = w::encode_request(&q).unwrap();
        assert_eq!(u16::from_le_bytes([bytes[4], bytes[5]]), w::SCHEMA_VERSION);
        w::decode_request(&bytes).unwrap();
        for unsupported in [0, 1, w::SCHEMA_VERSION - 1, w::SCHEMA_VERSION + 1] {
            let mut bad = bytes.clone();
            bad[4..6].copy_from_slice(&unsupported.to_le_bytes());
            assert!(w::decode_request(&bad).is_err());
            assert!(w::decode_control(&bad).is_err());
        }
    }
    let q = request(402, w::Command::Bootstrap);
    let bytes = w::encode_request(&q).unwrap();
    assert!(w::decode_control(&bytes).unwrap().is_none());
    let mut corrupt = bytes;
    let i = corrupt.len() - 1;
    corrupt[i] ^= 1;
    assert!(
        w::decode_request(&corrupt).is_err(),
        "ordinary fast-path hint never grants checksum authority"
    );
}

#[test]
fn cancellation_never_undoes_an_already_applied_catalog_selection() {
    let seed = seed();
    let mut host = DesktopHost::open(&seed.disk.0).unwrap();
    assert!(host.catalog_available());
    let access = host.authorize_window().unwrap();
    let boot = assert_editor(&seed, &mut host, &access, 500);
    assert!(matches!(
        boot.catalog,
        w::CatalogStatus::Ready {
            activation: w::CatalogActivationStatus::ActiveAfterRestart(_)
        }
    ));
    let backups = host.open_catalog_backup().unwrap();
    assert!(matches!(
        backup_call(&backups, &access, 510, cancel(&seed)),
        w::Outcome::Failure(_)
    ));
    assert!(host.catalog_available());
    assert!(
        backups.catalog_session_ready(),
        "refused cancellation does not pause a healthy catalog session"
    );
    assert!(
        !seed
            .disk
            .0
            .join("catalog-activation/00000000000000000001.cancelled.json")
            .exists()
    );
}
