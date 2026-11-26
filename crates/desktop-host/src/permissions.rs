//! The configured root is native shell policy. Unix permissions are checked here;
//! Windows requires a shell-selected user-private AppData ACL (std has no portable
//! DACL inspection API). This never treats a frontend path as authority.
use crate::{HostError, HostResult};
use std::{
    fs,
    path::{Path, PathBuf},
};
fn denied() -> HostError {
    HostError::Invalid(
        "Private desktop storage must be owned by this user and inaccessible to other users.",
    )
}
#[cfg(unix)]
fn unix_private(metadata: &fs::Metadata, directory: bool) -> HostResult<()> {
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    // OS effective identity is authoritative; environment variables are not.
    let uid = unsafe { libc::geteuid() };
    let wanted = if directory { 0o700 } else { 0o600 };
    if metadata.uid() != uid || metadata.permissions().mode() & 0o777 != wanted {
        return Err(denied());
    }
    Ok(())
}
fn check(path: &Path, directory: bool) -> HostResult<()> {
    let metadata = fs::symlink_metadata(path).map_err(|_| HostError::Unavailable)?;
    if metadata.file_type().is_symlink()
        || if directory {
            !metadata.is_dir()
        } else {
            !metadata.is_file()
        }
    {
        return Err(denied());
    }
    #[cfg(unix)]
    unix_private(&metadata, directory)?;
    Ok(())
}
pub(crate) fn root(configured: &Path) -> HostResult<PathBuf> {
    if !configured.exists() {
        let mut builder = fs::DirBuilder::new();
        builder.recursive(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            builder.mode(0o700);
        }
        builder
            .create(configured)
            .map_err(|_| HostError::Unavailable)?;
    }
    check(configured, true)?;
    let canonical = fs::canonicalize(configured).map_err(|_| HostError::Unavailable)?;
    check(&canonical, true)?;
    Ok(canonical)
}
pub(crate) fn database(root: &Path, name: &str) -> HostResult<PathBuf> {
    let path = root.join(name);
    // Primary files are created0600 before SQLite can apply a default022 umask.
    // create_new also refuses an already present symlink/racing substitution.
    let mut options = fs::OpenOptions::new();
    options.read(true).write(true).create_new(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.mode(0o600).custom_flags(libc::O_NOFOLLOW);
    }
    match options.open(&path) {
        Ok(file) => drop(file),
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
        Err(_) => return Err(HostError::Unavailable),
    }
    check(&path, false)?;
    for suffix in ["-wal", "-shm"] {
        let sidecar = root.join(format!("{name}{suffix}"));
        match fs::symlink_metadata(&sidecar) {
            Ok(_) => check(&sidecar, false)?,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(_) => return Err(HostError::Unavailable),
        }
    }
    Ok(path)
}
