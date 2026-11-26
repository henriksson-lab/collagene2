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

pub const HELP: &str = "Usage: collagene-desktop [--storage-root DIRECTORY] [--allow-loopback-http]\n       collagene-desktop --mcp [--storage-root DIRECTORY] [--mcp-document ID@REVISION ...]\n\nDesktop storage roots are operator configuration. MCP runs over stdio in this same\nexecutable, without opening a window. Saved MCP document selections require an\nexplicit root and exact library revision; MCP results remain session-only.";
#[derive(Debug)]
pub struct GuiOptions {
    pub storage_root: Option<PathBuf>,
    pub allow_loopback_http: bool,
}
#[derive(Debug)]
pub struct McpOptions {
    pub storage_root: Option<PathBuf>,
    pub documents: Vec<String>,
}
#[derive(Debug)]
pub enum Launch {
    Gui(GuiOptions),
    Mcp(McpOptions),
    Help,
    Version,
}
pub fn command_line() -> Result<Launch> {
    parse(std::env::args_os().skip(1))
}
fn parse(args: impl IntoIterator<Item = std::ffi::OsString>) -> Result<Launch> {
    let mut args = args.into_iter();
    let mut root = None;
    let mut mcp = false;
    let mut allow_loopback_http = false;
    let mut documents = Vec::new();
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            return Ok(Launch::Help);
        }
        if arg == "--version" || arg == "-V" {
            return Ok(Launch::Version);
        }
        if arg == "--mcp" {
            if mcp {
                return Err("--mcp may be supplied only once".into());
            }
            mcp = true;
        } else if arg == "--allow-loopback-http" {
            if allow_loopback_http {
                return Err("Loopback HTTP may be enabled only once".into());
            }
            allow_loopback_http = true;
        } else if arg == "--storage-root"
            || arg
                .to_str()
                .is_some_and(|s| s.starts_with("--storage-root="))
        {
            let path = match arg.to_str().and_then(|s| s.strip_prefix("--storage-root=")) {
                Some(s) if !s.is_empty() => std::ffi::OsString::from(s),
                Some(_) => return Err("--storage-root requires a nonempty directory".into()),
                None => args
                    .next()
                    .filter(|s| !s.is_empty() && !s.to_string_lossy().starts_with("--"))
                    .ok_or("--storage-root requires a directory")?,
            };
            if root.replace(PathBuf::from(path)).is_some() {
                return Err("--storage-root may be supplied only once".into());
            }
        } else if arg == "--mcp-document"
            || arg
                .to_str()
                .is_some_and(|s| s.starts_with("--mcp-document="))
        {
            let selection = match arg.to_str().and_then(|s| s.strip_prefix("--mcp-document=")) {
                Some(s) => s.to_owned(),
                None => args
                    .next()
                    .and_then(|s| s.into_string().ok())
                    .ok_or("--mcp-document requires ID@REVISION")?,
            };
            if selection.is_empty() || selection.len() > 256 || documents.len() >= 32 {
                return Err("Select at most32 saved documents as ID@REVISION".into());
            }
            documents.push(selection);
        } else {
            return Err("Unknown desktop argument; use --help".into());
        }
    }
    if mcp {
        if allow_loopback_http {
            return Err(
                "--allow-loopback-http applies only to desktop database synchronization".into(),
            );
        }
        if !documents.is_empty() && root.is_none() {
            return Err("Saved MCP selections require --storage-root".into());
        }
        Ok(Launch::Mcp(McpOptions {
            storage_root: root,
            documents,
        }))
    } else {
        if !documents.is_empty() {
            return Err("--mcp-document requires --mcp".into());
        }
        Ok(Launch::Gui(GuiOptions {
            storage_root: root,
            allow_loopback_http,
        }))
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
