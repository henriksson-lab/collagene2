use collagene_desktop_host::{DesktopExports, DesktopHost, WindowAccess};
use collagene_desktop_protocol as wire;
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
use wire::library::{Blob, Digest, Text};
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        Self(std::env::temp_dir().join(format!(
                "collagene-native-export-{}-{}-{}",
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
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn op(n: u64) -> wire::OperationId {
    wire::OperationId::new(format!("op-{n:032x}")).unwrap()
}
fn request(n: u64, name: &str, bytes: &[u8]) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: n,
        operation_id: op(n),
        command: wire::Command::ExportFile {
            suggested_name: Text::new(name).unwrap(),
            bytes: Blob::new(bytes.to_vec()).unwrap(),
        },
    }
}
fn call(
    service: &DesktopExports,
    access: &WindowAccess,
    request: &wire::RequestEnvelope,
) -> wire::Outcome {
    let frame = wire::encode_request(request).unwrap();
    let reply = service.dispatch_raw(access, &frame).unwrap();
    wire::decode_response_for(&reply, request).unwrap().outcome
}
fn receipt(outcome: wire::Outcome) -> (String, Digest, u64) {
    match outcome {
        wire::Outcome::Success(wire::Success::Exported { name, hash, bytes }) => {
            (name.to_string(), hash, bytes)
        }
        other => panic!("{other:?}"),
    }
}
#[test]
fn exact_biological_bytes_and_receipts_survive_restart_and_retry_without_duplicates() {
    let disk = Disk::new();
    let bytes = b"LOCUS exact\r\nORIGIN\n 1 aaccggtt\n//\n\0\xff";
    let req = request(1, "Plasmid sequence.gb", bytes);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let service = host.open_exports().unwrap();
    let first = receipt(call(&service, &access, &req));
    assert_eq!(first.1, Digest::of(bytes));
    assert_eq!(first.2, bytes.len() as u64);
    assert_eq!(
        std::fs::read(disk.0.join("exports").join(&first.0)).unwrap(),
        bytes
    );
    assert_eq!(receipt(call(&service, &access, &req)), first);
    drop(service);
    drop(host);
    let mut restarted = DesktopHost::open(&disk.0).unwrap();
    let access = restarted.authorize_window().unwrap();
    let service = restarted.open_exports().unwrap();
    assert_eq!(receipt(call(&service, &access, &req)), first);
    let mut changed = req.clone();
    changed.request_id = 999;
    changed.command = wire::Command::ExportFile {
        suggested_name: Text::new("changed.gb").unwrap(),
        bytes: Blob::new(bytes.to_vec()).unwrap(),
    };
    // Changed intent also fails after a new process/control registry.
    drop(service);
    drop(restarted);
    let mut restarted = DesktopHost::open(&disk.0).unwrap();
    let access = restarted.authorize_window().unwrap();
    let service = restarted.open_exports().unwrap();
    assert!(matches!(
        call(&service, &access, &changed),
        wire::Outcome::Failure(wire::Failure::Invalid { .. })
    ));
    let published = std::fs::read_dir(disk.0.join("exports"))
        .unwrap()
        .filter_map(Result::ok)
        .filter(|e| e.file_type().unwrap().is_file())
        .count();
    assert_eq!(published, 1);
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        assert_eq!(
            std::fs::metadata(disk.0.join("exports").join(&first.0))
                .unwrap()
                .permissions()
                .mode()
                & 0o777,
            0o600
        );
    }
}
#[test]
fn incomplete_and_foreign_files_are_never_accepted_or_overwritten_and_orphans_are_ignored() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let service = host.open_exports().unwrap();
    let req = request(2, "result.pdf", b"complete payload");
    let name = wire::export_file_name(&op(2), "result.pdf").unwrap();
    let path = disk.0.join("exports").join(name.as_str());
    let mut options = std::fs::OpenOptions::new();
    options.write(true).create_new(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.mode(0o600);
    }
    std::io::Write::write_all(&mut options.open(&path).unwrap(), b"partial").unwrap();
    assert!(matches!(
        call(&service, &access, &req),
        wire::Outcome::Failure(wire::Failure::Invalid { .. })
    ));
    assert_eq!(std::fs::read(&path).unwrap(), b"partial");
    let opdir = disk.0.join("exports/.operations").join(op(2).as_str());
    std::fs::write(opdir.join(".partial-abandoned"), b"old partial bytes").unwrap();
    // Simulate explicit operator removal of the mismatched foreign output, not
    // an automatic overwrite. A fresh host recovers only from exact request bytes.
    drop(service);
    drop(host);
    std::fs::remove_file(&path).unwrap();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let service = host.open_exports().unwrap();
    receipt(call(&service, &access, &req));
    assert_eq!(std::fs::read(&path).unwrap(), b"complete payload");
    assert_eq!(
        std::fs::read(opdir.join(".partial-abandoned")).unwrap(),
        b"old partial bytes"
    );
}
#[test]
fn scoped_controls_cancel_before_export_and_private_exports_do_not_need_editor_host_lock() {
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let service = host.open_exports().unwrap();
    let control = host.control();
    let req = request(3, "cancelled.csv", b"a,b\n");
    let cancel = wire::RequestEnvelope {
        request_id: 90,
        operation_id: op(90),
        command: wire::Command::Cancel {
            operation_id: op(3),
        },
    };
    control
        .dispatch_raw(&access, &wire::encode_request(&cancel).unwrap())
        .unwrap()
        .unwrap();
    assert!(matches!(
        call(&service, &access, &req),
        wire::Outcome::Failure(wire::Failure::Cancelled)
    ));
    assert!(
        !disk
            .0
            .join("exports/.operations")
            .join(op(3).as_str())
            .exists()
    );
    let host = std::sync::Mutex::new(host);
    let _editor_guard = host.lock().unwrap();
    // Separate native service can run while working/editor storage is locked.
    receipt(call(
        &service,
        &access,
        &request(4, "saved.fasta", b">x\nACGT\n"),
    ));
    let foreign_disk = Disk::new();
    let mut foreign = DesktopHost::open(&foreign_disk.0).unwrap();
    let foreign_window = foreign.authorize_window().unwrap();
    assert!(
        service
            .dispatch_raw(
                &foreign_window,
                &wire::encode_request(&request(5, "denied", b"x")).unwrap()
            )
            .is_err()
    );
}
#[cfg(unix)]
#[test]
fn export_symlinks_and_insecure_directories_fail_closed() {
    use std::os::unix::fs::{PermissionsExt, symlink};
    let disk = Disk::new();
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let service = host.open_exports().unwrap();
    let outside = disk.0.join("unrelated");
    std::fs::write(&outside, b"private unrelated file").unwrap();
    let name = wire::export_file_name(&op(6), "out.gb").unwrap();
    symlink(&outside, disk.0.join("exports").join(name.as_str())).unwrap();
    assert!(matches!(
        call(&service, &access, &request(6, "out.gb", b"new")),
        wire::Outcome::Failure(_)
    ));
    assert_eq!(std::fs::read(&outside).unwrap(), b"private unrelated file");
    drop(service);
    std::fs::set_permissions(
        disk.0.join("exports"),
        std::fs::Permissions::from_mode(0o755),
    )
    .unwrap();
    assert!(host.open_exports().is_err());
}

#[test]
fn two_native_connections_cannot_publish_different_payloads_under_one_operation() {
    let disk = Disk::new();
    let mut first = DesktopHost::open(&disk.0).unwrap();
    let mut second = DesktopHost::open(&disk.0).unwrap();
    let first_access = first.authorize_window().unwrap();
    let second_access = second.authorize_window().unwrap();
    let one = first.open_exports().unwrap();
    let two = second.open_exports().unwrap();
    let barrier = std::sync::Arc::new(std::sync::Barrier::new(2));
    let request_a = request(8, "shared.gb", b"first exact payload");
    let request_b = request(8, "shared.gb", b"second exact payload");
    let launch = |service: DesktopExports,
                  access: WindowAccess,
                  request: wire::RequestEnvelope,
                  barrier: std::sync::Arc<std::sync::Barrier>| {
        std::thread::spawn(move || {
            barrier.wait();
            call(&service, &access, &request)
        })
    };
    let a = launch(one, first_access, request_a, barrier.clone());
    let b = launch(two, second_access, request_b, barrier);
    let outcomes = [a.join().unwrap(), b.join().unwrap()];
    assert_eq!(
        outcomes
            .iter()
            .filter(|o| matches!(o, wire::Outcome::Success(_)))
            .count(),
        1
    );
    assert_eq!(
        outcomes
            .iter()
            .filter(|o| matches!(o, wire::Outcome::Failure(wire::Failure::Invalid { .. })))
            .count(),
        1
    );
    let name = wire::export_file_name(&op(8), "shared.gb").unwrap();
    let bytes = std::fs::read(disk.0.join("exports").join(name.as_str())).unwrap();
    assert!(bytes == b"first exact payload" || bytes == b"second exact payload");
    // A persisted descriptor binds the content even if the display name is unchanged.
    drop(first);
    drop(second);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let service = host.open_exports().unwrap();
    assert!(matches!(
        call(
            &service,
            &access,
            &request(8, "shared.gb", b"third changed payload")
        ),
        wire::Outcome::Failure(wire::Failure::Invalid { .. })
    ));
}
