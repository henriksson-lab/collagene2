//! Native operator configuration. No IPC command accepts or changes these paths.
use serde::{Deserialize, Serialize};
use std::{
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const MAX_SETTINGS_BYTES: u64 = 8192;
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Settings {
    version: u16,
    storage_root: PathBuf,
}

pub fn command_line() -> Result<Option<PathBuf>> {
    let args = std::env::args_os().skip(1).collect::<Vec<_>>();
    match args.as_slice() {
        [] => Ok(None),
        [option, path] if option == "--storage-root" && !path.is_empty() => {
            Ok(Some(PathBuf::from(path)))
        }
        _ => Err("Usage: collagene-desktop [--storage-root DIRECTORY]".into()),
    }
}
pub fn select(
    config_dir: &Path,
    default_root: &Path,
    override_root: Option<PathBuf>,
) -> Result<PathBuf> {
    if let Some(root) = override_root {
        return Ok(root);
    }
    let path = config_dir.join("desktop.json");
    let metadata = match fs::symlink_metadata(&path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            return Ok(default_root.into());
        }
        Err(error) => return Err(error.into()),
    };
    if metadata.file_type().is_symlink()
        || !metadata.is_file()
        || metadata.len() > MAX_SETTINGS_BYTES
    {
        return Err("Desktop settings must be a bounded regular file".into());
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if metadata.mode() & 0o077 != 0 {
            return Err("Desktop settings must be owner-only".into());
        }
    }
    let mut bytes = Vec::new();
    fs::File::open(path)?
        .take(MAX_SETTINGS_BYTES + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_SETTINGS_BYTES {
        return Err("Desktop settings exceed their size limit".into());
    }
    let settings: Settings = serde_json::from_slice(&bytes)?;
    if settings.version != 1 || !settings.storage_root.is_absolute() {
        return Err("Unsupported desktop storage-root configuration".into());
    }
    Ok(settings.storage_root)
}
pub fn persist(config_dir: &Path, root: &Path) -> Result<()> {
    let root = root.canonicalize()?;
    let settings = Settings {
        version: 1,
        storage_root: root,
    };
    let bytes = serde_json::to_vec(&settings)?;
    if bytes.len() as u64 > MAX_SETTINGS_BYTES {
        return Err("Desktop settings exceed their size limit".into());
    }
    let mut builder = fs::DirBuilder::new();
    builder.recursive(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        builder.mode(0o700);
    }
    builder.create(config_dir)?;
    if fs::symlink_metadata(config_dir)?.file_type().is_symlink() {
        return Err("Desktop configuration directory cannot be a symlink".into());
    }
    let path = config_dir.join("desktop.json");
    if fs::symlink_metadata(&path).is_ok_and(|m| m.file_type().is_symlink() || !m.is_file()) {
        return Err("Desktop settings must be a regular file".into());
    }
    let mut nonce = [0u8; 16];
    getrandom::fill(&mut nonce).map_err(|_| "Cannot initialize atomic desktop settings")?;
    let suffix = nonce.iter().map(|b| format!("{b:02x}")).collect::<String>();
    let temporary = config_dir.join(format!(".desktop-{suffix}.tmp"));
    let result = (|| -> Result<()> {
        let mut options = fs::OpenOptions::new();
        options.create_new(true).write(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        let mut file = options.open(&temporary)?;
        file.write_all(&bytes)?;
        file.sync_all()?;
        drop(file);
        fs::rename(&temporary, &path)?;
        #[cfg(unix)]
        {
            fs::File::open(config_dir)?.sync_all()?;
        }
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(temporary);
    }
    result
}
