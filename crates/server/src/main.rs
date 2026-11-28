use collagene_lab_library::{LabId, Library, StorageOptions, backup};
use collagene_server::{
    AccessKey, AssetSnapshot, Error, Limits, MemberRegistry, Result, ServerState,
};
use std::{
    collections::BTreeMap,
    net::SocketAddr,
    path::{Path, PathBuf},
};

#[tokio::main]
async fn main() {
    if let Err(error) = run().await {
        eprintln!("Server could not complete the requested action: {error}");
        std::process::exit(1);
    }
}
async fn run() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let command = args.next().ok_or(Error::Configuration(USAGE))?;
    if matches!(command.as_str(), "help" | "--help" | "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    let mut options = BTreeMap::new();
    while let Some(flag) = args.next() {
        if !flag.starts_with("--") {
            return Err(Error::Configuration(USAGE));
        }
        let value = args.next().ok_or(Error::Configuration(USAGE))?;
        if options.insert(flag, value).is_some() {
            return Err(Error::Configuration("duplicate option"));
        }
    }
    match command.as_str() {
        "keygen" => {
            allowed(&options, &["--key-file"])?;
            AccessKey::generate_file(Path::new(required(&options, "--key-file")?))?;
            println!("Access key file created; its contents were not printed.");
        }
        "init" => {
            allowed(
                &options,
                &["--storage-root", "--lab-name", "--max-lab-bytes"],
            )?;
            let root = storage_root(Path::new(required(&options, "--storage-root")?), true)?;
            let mut database = Library::open_with_options(
                database_path(&root, false)?,
                storage_options(&options)?,
            )?;
            let lab = database.create_lab(required(&options, "--lab-name")?)?;
            println!(
                "Created shared lab {}. Retain this lab ID for trusted server configuration.",
                lab.id
            );
        }
        "backup" | "verify-backup" | "restore-new" => {
            backup_command(&command, &options)?;
        }
        "serve" => {
            allowed(
                &options,
                &[
                    "--storage-root",
                    "--lab-id",
                    "--members-file",
                    "--assets",
                    "--bind",
                    "--origin",
                    "--max-concurrent",
                    "--max-lab-bytes",
                    "--mcp-storage-root",
                ],
            )?;
            let bind: SocketAddr = options
                .get("--bind")
                .map(String::as_str)
                .unwrap_or("127.0.0.1:8084")
                .parse()
                .map_err(|_| Error::Configuration("bind address must be IP:port"))?;
            if !bind.ip().is_loopback() && !options.contains_key("--origin") {
                return Err(Error::Configuration(
                    "non-loopback bind requires explicit --origin",
                ));
            }
            if bind.port() == 0 {
                return Err(Error::Configuration("serve requires a fixed port"));
            }
            let origin = options
                .get("--origin")
                .cloned()
                .unwrap_or_else(|| format!("http://{bind}"));
            let root = storage_root(Path::new(required(&options, "--storage-root")?), false)?;
            let assets_path = Path::new(required(&options, "--assets")?);
            if std::fs::symlink_metadata(assets_path)?
                .file_type()
                .is_symlink()
            {
                return Err(Error::Configuration("asset root symlink"));
            }
            let assets_root = assets_path.canonicalize()?;
            if root.starts_with(&assets_root) || assets_root.starts_with(&root) {
                return Err(Error::Configuration(
                    "storage and asset roots must be separate",
                ));
            }
            let members = MemberRegistry::read_file(
                Path::new(required(&options, "--members-file")?),
                &assets_root,
            )?;
            let assets = AssetSnapshot::load(&assets_root)?;
            let file_count = assets.file_count();
            let database = Library::open_with_options(
                database_path(&root, true)?,
                storage_options(&options)?,
            )?;
            let lab_id = LabId::parse(required(&options, "--lab-id")?)?;
            let access = database.trusted_host_access(&lab_id)?;
            let mut limits = Limits::default();
            if let Some(value) = options.get("--max-concurrent") {
                limits.max_concurrent_library_requests = value
                    .parse()
                    .map_err(|_| Error::Configuration("max concurrent"))?;
            }
            let mut state = ServerState::new(database, access, members, &origin, assets, limits)?;
            if let Some(path) = options.get("--mcp-storage-root") {
                let mcp_root = collagene_local_storage::StorageRoot::existing(path)?;
                state = state.with_hosted_mcp(mcp_root, Default::default())?;
            }
            let shutdown = state.clone();
            let listener = tokio::net::TcpListener::bind(bind).await?;
            println!(
                "Shared-library server on {bind}; {file_count} immutable static assets; lab {lab_id}. Individual member credentials are required for API requests."
            );
            axum::serve(listener, collagene_server::router(state))
                .with_graceful_shutdown(async move {
                    let _ = tokio::signal::ctrl_c().await;
                    shutdown.shutdown_hosted_mcp();
                })
                .await?;
        }
        _ => return Err(Error::Configuration(USAGE)),
    }
    Ok(())
}
fn backup_command(command: &str, options: &BTreeMap<String, String>) -> Result<()> {
    let permitted = match command {
        "backup" => vec![
            "--storage-root",
            "--lab-id",
            "--backup-id",
            "--max-backup-bytes",
            "--max-lab-bytes",
        ],
        "verify-backup" => vec![
            "--storage-root",
            "--lab-id",
            "--backup-id",
            "--max-backup-bytes",
        ],
        "restore-new" => vec![
            "--storage-root",
            "--lab-id",
            "--backup-id",
            "--max-backup-bytes",
            "--destination-root",
        ],
        _ => return Err(Error::Configuration(USAGE)),
    };
    allowed(options, &permitted)?;
    let root = storage_root(Path::new(required(options, "--storage-root")?), false)?;
    let path = backup_path(
        &root,
        required(options, "--backup-id")?,
        command == "backup",
    )?;
    let mut policy = backup::Policy::single_lab(LabId::parse(required(options, "--lab-id")?)?);
    if let Some(n) = options.get("--max-backup-bytes") {
        policy.max_database_bytes = n
            .parse()
            .map_err(|_| Error::Configuration("backup byte bound"))?;
    }
    let manifest = match command {
        "backup" => {
            let db =
                Library::open_with_options(database_path(&root, true)?, storage_options(options)?)?;
            backup::create(&db, &path, &policy)?
        }
        "verify-backup" => backup::verify(&path, &policy)?,
        "restore-new" => backup::restore_new(
            &path,
            Path::new(required(options, "--destination-root")?),
            &policy,
        )?,
        _ => return Err(Error::Configuration(USAGE)),
    };
    println!(
        "{} completed: {} checked database bytes; {} lab(s); immutable backup ID {}. No access-key contents were read or printed.",
        command,
        manifest.database_bytes,
        manifest.labs.len(),
        required(options, "--backup-id")?
    );
    println!(
        "Recorded archive quota: {} bytes per lab. Restore does not change runtime server settings; supply the intended --max-lab-bytes when serving the recovered database.",
        manifest.max_lab_archive_bytes
    );
    Ok(())
}
fn backup_path(root: &Path, id: &str, create: bool) -> Result<PathBuf> {
    if id.is_empty()
        || id.len() > 64
        || !id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_'))
    {
        return Err(Error::Configuration(
            "backup ID must be1–64 ASCII letters, digits, hyphens or underscores",
        ));
    }
    let directory = root.join("backups");
    let existed = directory.exists();
    storage_root(&directory, create)?;
    #[cfg(unix)]
    if create && !existed {
        std::fs::File::open(root)?.sync_all()?;
    }
    let _ = existed;
    Ok(directory.join(format!("backup-{id}.cglbak")))
}
fn allowed(options: &BTreeMap<String, String>, names: &[&str]) -> Result<()> {
    if options.keys().any(|key| !names.contains(&key.as_str())) {
        return Err(Error::Configuration("unknown option"));
    }
    Ok(())
}
fn required<'a>(options: &'a BTreeMap<String, String>, name: &str) -> Result<&'a str> {
    options
        .get(name)
        .map(String::as_str)
        .filter(|value| !value.is_empty())
        .ok_or(Error::Configuration(USAGE))
}
fn storage_options(options: &BTreeMap<String, String>) -> Result<StorageOptions> {
    let mut storage = StorageOptions::default();
    if let Some(value) = options.get("--max-lab-bytes") {
        storage.max_lab_archive_bytes = value
            .parse()
            .map_err(|_| Error::Configuration("lab byte quota"))?;
    }
    Ok(storage)
}
fn storage_root(path: &Path, create: bool) -> Result<PathBuf> {
    if !path.exists() && create {
        let mut builder = std::fs::DirBuilder::new();
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            builder.mode(0o700);
        }
        builder.create(path)?;
    }
    let metadata = std::fs::symlink_metadata(path)?;
    if !metadata.is_dir() || metadata.file_type().is_symlink() {
        return Err(Error::Configuration(
            "storage root must be a real directory",
        ));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if metadata.mode() & 0o077 != 0 {
            return Err(Error::Configuration("storage root must be owner-only"));
        }
    }
    Ok(path.canonicalize()?)
}
fn database_path(root: &Path, must_exist: bool) -> Result<PathBuf> {
    let path = root.join("library.sqlite");
    for suffix in ["library.sqlite", "library.sqlite-wal", "library.sqlite-shm"] {
        match std::fs::symlink_metadata(root.join(suffix)) {
            Ok(metadata) if !metadata.is_file() || metadata.file_type().is_symlink() => {
                return Err(Error::Configuration(
                    "database/sidecars must be regular files",
                ));
            }
            Ok(_) => (),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => (),
            Err(error) => return Err(error.into()),
        }
    }
    if must_exist && !path.is_file() {
        return Err(Error::Configuration(
            "initialize the lab database before serving",
        ));
    }
    Ok(path)
}
const USAGE: &str = "Usage:\n  collagene-server keygen --key-file PATH\n  collagene-server init --storage-root DIRECTORY --lab-name NAME [--max-lab-bytes N]\n  collagene-server backup --storage-root DIRECTORY --lab-id ID --backup-id NAME [--max-backup-bytes N] [--max-lab-bytes N]\n  collagene-server verify-backup --storage-root DIRECTORY --lab-id ID --backup-id NAME [--max-backup-bytes N]\n  collagene-server restore-new --storage-root DIRECTORY --lab-id ID --backup-id NAME --destination-root NEW_DIRECTORY [--max-backup-bytes N]\n  collagene-server serve --storage-root DIRECTORY --lab-id ID --members-file PATH --assets DIRECTORY [--bind 127.0.0.1:8084] [--origin https://your-host] [--max-concurrent N] [--max-lab-bytes N] [--mcp-storage-root DIRECTORY]";

#[cfg(test)]
mod backup_cli_tests {
    use super::*;
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
    struct Temp(PathBuf);
    impl Temp {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "collagene-backup-cli-{}-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos(),
                NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            ));
            storage_root(&path, true).unwrap();
            Self(path)
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    #[test]
    fn native_backup_verify_restore_commands_use_fixed_confined_paths_and_never_overwrite() {
        let temp = Temp::new();
        let mut db = Library::open(temp.0.join("library.sqlite")).unwrap();
        let lab = db.create_lab("Temporary operator fixture").unwrap();
        drop(db);
        let mut options = BTreeMap::from([
            (
                "--storage-root".into(),
                temp.0.to_string_lossy().into_owned(),
            ),
            ("--lab-id".into(), lab.id.to_string()),
            ("--backup-id".into(), "nightly_001".into()),
        ]);
        backup_command("backup", &options).unwrap();
        let backup = backup_path(&temp.0, "nightly_001", false).unwrap();
        let before = std::fs::read(&backup).unwrap();
        backup_command("verify-backup", &options).unwrap();
        assert!(backup_command("backup", &options).is_err());
        assert_eq!(std::fs::read(&backup).unwrap(), before);
        std::fs::rename(
            temp.0.join("library.sqlite"),
            temp.0.join("offline-original.sqlite"),
        )
        .unwrap();
        backup_command("verify-backup", &options).unwrap();
        assert!(
            !temp.0.join("library.sqlite").exists(),
            "verification never opens or recreates the source database"
        );
        let destination = temp.0.join("new-recovery");
        options.insert(
            "--destination-root".into(),
            destination.to_string_lossy().into_owned(),
        );
        backup_command("restore-new", &options).unwrap();
        assert!(backup_command("restore-new", &options).is_err());
        let restored = Library::open(destination.join("library.sqlite")).unwrap();
        let access = restored.trusted_host_access(&lab.id).unwrap();
        assert_eq!(
            restored.lab(&access).unwrap().name,
            "Temporary operator fixture"
        );
        assert!(
            backup_command("backup", &options).is_err(),
            "restore-only option is rejected for backup"
        );
        for id in ["", "../other", "a/b", "a\\b", ".", "..\0", "with space"] {
            assert!(backup_path(&temp.0, id, false).is_err());
        }
    }
    #[cfg(unix)]
    #[test]
    fn backup_directory_symlink_is_not_followed() {
        use std::os::unix::fs::symlink;
        let temp = Temp::new();
        let outside = Temp::new();
        symlink(&outside.0, temp.0.join("backups")).unwrap();
        assert!(backup_path(&temp.0, "snapshot", true).is_err());
        assert_eq!(std::fs::read_dir(&outside.0).unwrap().count(), 0);
    }
}
