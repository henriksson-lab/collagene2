use collagene_lab_library::{LabId, Library, StorageOptions};
use collagene_server::{AccessKey, AssetSnapshot, Error, Limits, Result, ServerState};
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
        "serve" => {
            allowed(
                &options,
                &[
                    "--storage-root",
                    "--lab-id",
                    "--key-file",
                    "--assets",
                    "--bind",
                    "--origin",
                    "--max-concurrent",
                    "--max-lab-bytes",
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
            let key_path = Path::new(required(&options, "--key-file")?);
            if key_path.canonicalize()?.starts_with(&assets_root) {
                return Err(Error::Configuration(
                    "access key must be outside static asset root",
                ));
            }
            let key = AccessKey::read(key_path)?;
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
            let state = ServerState::new(database, access, key, &origin, assets, limits)?;
            let listener = tokio::net::TcpListener::bind(bind).await?;
            println!(
                "Shared-library server on {bind}; {file_count} immutable static assets; lab {lab_id}. Access key is required for library requests."
            );
            axum::serve(listener, collagene_server::router(state))
                .with_graceful_shutdown(async {
                    let _ = tokio::signal::ctrl_c().await;
                })
                .await?;
        }
        _ => return Err(Error::Configuration(USAGE)),
    }
    Ok(())
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
const USAGE: &str = "Usage:\n  collagene-server keygen --key-file PATH\n  collagene-server init --storage-root DIRECTORY --lab-name NAME [--max-lab-bytes N]\n  collagene-server serve --storage-root DIRECTORY --lab-id ID --key-file PATH --assets DIRECTORY [--bind 127.0.0.1:8084] [--origin https://your-host] [--max-concurrent N] [--max-lab-bytes N]";
