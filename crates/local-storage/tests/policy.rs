use collagene_local_storage::StorageRoot;
use std::{
    fs,
    io::{Read, Seek, SeekFrom, Write},
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};
static SERIAL: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "collagene-storage-test-{}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
#[test]
fn root_must_already_be_a_directory_and_is_canonical() {
    let f = Fixture::new();
    assert!(StorageRoot::existing("").is_err());
    let missing = f.0.join("missing");
    assert!(StorageRoot::existing(&missing).is_err());
    assert!(!missing.exists());
    let file = f.0.join("ordinary-file");
    fs::write(&file, b"unchanged").unwrap();
    assert!(StorageRoot::existing(&file).is_err());
    let folder = f.0.join("safe area ü");
    fs::create_dir(&folder).unwrap();
    let root = StorageRoot::existing(folder.join(".")).unwrap();
    assert_eq!(root.canonical_path(), fs::canonicalize(&folder).unwrap());
    assert_eq!(fs::read_dir(&folder).unwrap().count(), 0);
}
#[test]
fn isolated_files_are_retained_until_last_handle_then_cleaned_without_deleting_root() {
    let f = Fixture::new();
    fs::write(f.0.join("user-file"), b"keep").unwrap();
    let root = StorageRoot::existing(&f.0).unwrap();
    let first = root.session().unwrap();
    let second = root.session().unwrap();
    assert_ne!(first.namespace(), second.namespace());
    assert_ne!(first.diagnostic_path(), second.diagnostic_path());
    assert_eq!(
        first.diagnostic_path().parent().unwrap(),
        root.canonical_path()
    );
    let path = first.diagnostic_path();
    let mut file = first.create_file().unwrap();
    file.file().write_all(b"bounded managed bytes").unwrap();
    file.file().seek(SeekFrom::Start(0)).unwrap();
    let mut bytes = vec![];
    file.file().read_to_end(&mut bytes).unwrap();
    assert_eq!(bytes, b"bounded managed bytes");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        assert_eq!(
            fs::metadata(&path).unwrap().permissions().mode() & 0o777,
            0o700
        );
        let name = fs::read_dir(&path).unwrap().next().unwrap().unwrap().path();
        assert_eq!(
            fs::metadata(name).unwrap().permissions().mode() & 0o777,
            0o600
        );
    }
    drop(first);
    assert!(path.exists());
    drop(file);
    assert!(!path.exists());
    drop(second);
    assert_eq!(fs::read_dir(&f.0).unwrap().count(), 1);
    assert_eq!(fs::read(f.0.join("user-file")).unwrap(), b"keep");
}
#[cfg(unix)]
#[test]
fn configured_symlinks_are_rejected_including_trailing_dot() {
    use std::os::unix::fs::symlink;
    let f = Fixture::new();
    let outside = Fixture::new();
    let link = f.0.join("link");
    symlink(&outside.0, &link).unwrap();
    assert!(StorageRoot::existing(&link).is_err());
    assert!(StorageRoot::existing(link.join(".")).is_err());
    assert_eq!(fs::read_dir(&outside.0).unwrap().count(), 0);
}
#[cfg(unix)]
#[test]
fn preexisting_file_symlink_cannot_be_opened_or_overwritten() {
    use std::os::unix::fs::symlink;
    let f = Fixture::new();
    let outside = Fixture::new();
    let sentinel = outside.0.join("sentinel");
    fs::write(&sentinel, b"private outside").unwrap();
    let root = StorageRoot::existing(&f.0).unwrap();
    let session = root.session().unwrap();
    symlink(&sentinel, session.diagnostic_path().join("managed-0.bin")).unwrap();
    assert!(session.create_file().is_err());
    assert_eq!(fs::read(&sentinel).unwrap(), b"private outside");
    drop(session);
    assert_eq!(fs::read(&sentinel).unwrap(), b"private outside");
}
#[cfg(unix)]
#[test]
fn session_path_replacement_never_redirects_new_io_or_cleanup() {
    use std::os::unix::fs::symlink;
    let f = Fixture::new();
    let outside = Fixture::new();
    let sentinel = outside.0.join("managed-0.bin");
    fs::write(&sentinel, b"outside sentinel").unwrap();
    let root = StorageRoot::existing(&f.0).unwrap();
    let session = root.session().unwrap();
    let path = session.diagnostic_path();
    let mut file = session.create_file().unwrap();
    file.file().write_all(b"owned").unwrap();
    let moved = f.0.join("moved-session");
    fs::rename(&path, &moved).unwrap();
    symlink(&outside.0, &path).unwrap();
    assert!(session.create_file().is_err());
    drop(file);
    assert_eq!(fs::read_dir(&moved).unwrap().count(), 0);
    drop(session);
    assert_eq!(fs::read(&sentinel).unwrap(), b"outside sentinel");
    assert_eq!(fs::read_dir(&outside.0).unwrap().count(), 1);
}
#[cfg(unix)]
#[test]
fn root_path_replacement_does_not_replace_the_opened_root_capability() {
    use std::os::unix::fs::symlink;
    let f = Fixture::new();
    let configured = f.0.join("configured");
    let moved = f.0.join("original-root");
    fs::create_dir(&configured).unwrap();
    let outside = Fixture::new();
    let root = StorageRoot::existing(&configured).unwrap();
    fs::rename(&configured, &moved).unwrap();
    symlink(&outside.0, &configured).unwrap();
    let session = root.session().unwrap();
    let name = PathBuf::from(session.diagnostic_path().file_name().unwrap());
    assert!(moved.join(name).exists());
    let mut file = session.create_file().unwrap();
    file.file().write_all(b"still original directory").unwrap();
    assert_eq!(fs::read_dir(&outside.0).unwrap().count(), 0);
    drop(file);
    drop(session);
    assert_eq!(fs::read_dir(&moved).unwrap().count(), 0);
}

#[cfg(unix)]
#[test]
fn observed_regular_file_replacement_is_never_deleted_on_cleanup() {
    let f = Fixture::new();
    let root = StorageRoot::existing(&f.0).unwrap();
    let session = root.session().unwrap();
    let mut file = session.create_file().unwrap();
    file.file().write_all(b"managed original").unwrap();
    let path = session.diagnostic_path();
    let owned = path.join("managed-0.bin");
    let moved = path.join("renamed-owned.bin");
    fs::rename(&owned, &moved).unwrap();
    fs::write(&owned, b"unrelated replacement").unwrap();
    drop(file);
    assert_eq!(fs::read(&owned).unwrap(), b"unrelated replacement");
    assert_eq!(fs::read(&moved).unwrap(), b"managed original");
    drop(session);
    assert_eq!(fs::read(&owned).unwrap(), b"unrelated replacement");
}
#[cfg(windows)]
#[test]
fn windows_owned_file_denies_delete_sharing_and_uses_handle_cleanup() {
    let f = Fixture::new();
    let root = StorageRoot::existing(&f.0).unwrap();
    let session = root.session().unwrap();
    let file = session.create_file().unwrap();
    let path = session.diagnostic_path();
    assert!(fs::rename(path.join("managed-0.bin"), path.join("replacement.bin")).is_err());
    drop(file);
    assert!(fs::read_dir(&path).unwrap().next().is_none());
    drop(session);
    assert!(!path.exists());
}
