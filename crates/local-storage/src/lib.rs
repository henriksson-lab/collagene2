//! Managed native storage beneath an explicitly selected directory capability.
//!
//! Only initialization uses an ambient host path. All later creation and cleanup
//! use retained directory handles and generated single-component names. This is an
//! application I/O boundary, not an OS sandbox against other same-user processes.
use cap_fs_ext::DirExt;
use cap_std::{
    ambient_authority,
    fs::{Dir, DirBuilder, OpenOptions},
};
use std::{
    fs::File,
    io,
    path::{Component, Path, PathBuf},
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
    },
    time::{SystemTime, UNIX_EPOCH},
};
static SESSION_SERIAL: AtomicU64 = AtomicU64::new(0);

/// An existing canonical directory selected by trusted native configuration.
/// No filesystem path is accepted from a browser or MCP tool through this type.
#[derive(Clone, Debug)]
pub struct StorageRoot {
    inner: Arc<Root>,
}
#[derive(Debug)]
struct Root {
    directory: Dir,
    canonical: PathBuf,
}
impl StorageRoot {
    /// Open an existing folder. The final component may not be a symlink/junction.
    /// Parent aliases are resolved once (e.g. macOS /tmp); later I/O never reopens
    /// an absolute path, so a symlink replacement cannot redirect it elsewhere.
    pub fn existing(path: impl AsRef<Path>) -> io::Result<Self> {
        let path = path.as_ref();
        if path.as_os_str().is_empty() {
            return Err(invalid("storage root path is empty"));
        }
        let normalized = path.components().collect::<PathBuf>();
        let metadata = std::fs::symlink_metadata(&normalized)?;
        if metadata.file_type().is_symlink() || !metadata.is_dir() {
            return Err(invalid(
                "storage root must be an existing directory, not a file or symlink",
            ));
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::MetadataExt;
            if metadata.file_attributes() & 0x400 != 0 {
                return Err(invalid(
                    "storage root may not be a reparse point or junction",
                ));
            }
        }
        let canonical = std::fs::canonicalize(&normalized)?;
        let directory =
            if let (Some(parent), Some(name)) = (canonical.parent(), canonical.file_name()) {
                Dir::open_ambient_dir(parent, ambient_authority())?.open_dir_nofollow(name)?
            } else {
                Dir::open_ambient_dir(&canonical, ambient_authority())?
            };
        #[cfg(unix)]
        {
            use cap_std::fs::MetadataExt as CapMetadataExt;
            use std::os::unix::fs::MetadataExt as StdMetadataExt;
            let opened = directory.dir_metadata()?;
            if opened.dev() != metadata.dev() || opened.ino() != metadata.ino() {
                return Err(io::Error::new(
                    io::ErrorKind::PermissionDenied,
                    "storage root changed while opening it",
                ));
            }
        }
        Ok(Self {
            inner: Arc::new(Root {
                directory,
                canonical,
            }),
        })
    }
    /// Preserve the default native behavior: an isolated child of the OS temp folder.
    pub fn temporary() -> io::Result<Self> {
        Self::existing(std::fs::canonicalize(std::env::temp_dir())?)
    }
    /// Diagnostic location only. Managed I/O never uses this path to open a file.
    pub fn canonical_path(&self) -> &Path {
        &self.inner.canonical
    }
    pub fn session(&self) -> io::Result<ManagedSession> {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|_| io::Error::other("invalid storage clock"))?
            .as_nanos();
        let serial = SESSION_SERIAL
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
            .map_err(|_| io::Error::other("managed session counter exhausted"))?;
        let namespace = format!("{}-{nanos}-{serial}", std::process::id());
        let name = PathBuf::from(format!("collagene2-session-{namespace}"));
        let builder = DirBuilder::new();
        #[cfg(unix)]
        let builder = {
            use cap_std::fs::DirBuilderExt;
            let mut builder = builder;
            builder.mode(0o700);
            builder
        };
        self.inner.directory.create_dir_with(&name, &builder)?;
        let directory = match self.inner.directory.open_dir_nofollow(&name) {
            Ok(d) => d,
            Err(error) => {
                let _ = self.inner.directory.remove_dir(&name);
                return Err(error);
            }
        };
        Ok(ManagedSession {
            inner: Arc::new(Session {
                root: self.inner.clone(),
                directory: Some(directory),
                name,
                namespace,
                serial: AtomicU64::new(0),
            }),
        })
    }
}
fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message)
}
#[derive(Clone, Debug)]
pub struct ManagedSession {
    inner: Arc<Session>,
}
#[derive(Debug)]
struct Session {
    root: Arc<Root>,
    directory: Option<Dir>,
    name: PathBuf,
    namespace: String,
    serial: AtomicU64,
}
impl Session {
    fn directory(&self) -> &Dir {
        self.directory.as_ref().expect("live session directory")
    }
    fn check_name(&self) -> io::Result<()> {
        let current = self.root.directory.open_dir_nofollow(&self.name)?;
        #[cfg(unix)]
        {
            use cap_std::fs::MetadataExt;
            let actual = current.dir_metadata()?;
            let expected = self.directory().dir_metadata()?;
            if actual.dev() != expected.dev() || actual.ino() != expected.ino() {
                return Err(io::Error::new(
                    io::ErrorKind::PermissionDenied,
                    "managed session directory was replaced",
                ));
            }
        }
        #[cfg(not(unix))]
        drop(current);
        Ok(())
    }
}
impl ManagedSession {
    pub fn namespace(&self) -> &str {
        &self.inner.namespace
    }
    /// Diagnostic location only, never an instruction to reopen it using std::fs.
    pub fn diagnostic_path(&self) -> PathBuf {
        self.inner.root.canonical.join(&self.inner.name)
    }
    /// Allocate an opaque new file; callers cannot choose a path or overwrite a file.
    pub fn create_file(&self) -> io::Result<ManagedFile> {
        self.inner.check_name()?;
        let serial = self
            .inner
            .serial
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
            .map_err(|_| io::Error::other("managed file counter exhausted"))?;
        let name = PathBuf::from(format!("managed-{serial}.bin"));
        debug_assert!(matches!(
            name.components().next(),
            Some(Component::Normal(_))
        ));
        let mut options = OpenOptions::new();
        options.create_new(true).read(true).write(true);
        #[cfg(unix)]
        {
            use cap_std::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        #[cfg(windows)]
        {
            use cap_std::fs::OpenOptionsExt;
            // Generic read/write + DELETE; share read/write, never share DELETE.
            // Windows deletes this exact owned file on final handle close, without
            // reopening a pathname that could name an unrelated replacement.
            options
                .access_mode(0x8000_0000 | 0x4000_0000 | 0x0001_0000)
                .share_mode(0x1 | 0x2)
                .custom_flags(0x0400_0000); // DELETE_ON_CLOSE
        }
        let file = self
            .inner
            .directory()
            .open_with(&name, &options)?
            .into_std();
        Ok(ManagedFile {
            file: Some(file),
            _session: self.clone(),
            _name: name,
        })
    }
}
/// An owned file handle, not a host path. Unix unlinks its matching entry while
/// still open; Windows uses delete-on-close against the pinned handle.
#[derive(Debug)]
pub struct ManagedFile {
    file: Option<File>,
    _session: ManagedSession,
    _name: PathBuf,
}
impl ManagedFile {
    pub fn file(&mut self) -> &mut File {
        self.file.as_mut().expect("live managed file")
    }
}
impl Drop for ManagedFile {
    fn drop(&mut self) {
        #[cfg(unix)]
        {
            use cap_std::fs::MetadataExt as CapMetadataExt;
            use std::os::unix::fs::MetadataExt as StdMetadataExt;
            let owned = self.file.as_ref().and_then(|f| f.metadata().ok());
            let current = self
                ._session
                .inner
                .directory()
                .symlink_metadata(&self._name)
                .ok();
            if let (Some(owned), Some(current)) = (owned, current)
                && owned.dev() == current.dev()
                && owned.ino() == current.ino()
            {
                // Keep the handle alive through unlink and never delete an
                // observed replacement. Same-user concurrent renames are not
                // an OS sandbox guarantee; the directory capability still
                // confines this operation to the managed session.
                let _ = self
                    ._session
                    .inner
                    .directory()
                    .remove_file_or_symlink(&self._name);
            }
        }
        // On Windows DELETE_ON_CLOSE handles cleanup, with SHARE_DELETE denied.
        drop(self.file.take());
    }
}
impl Drop for Session {
    fn drop(&mut self) {
        let owns_name = self.check_name().is_ok();
        #[cfg(unix)]
        {
            if owns_name {
                let _ = self.root.directory.remove_dir(&self.name);
            }
            drop(self.directory.take());
        }
        #[cfg(not(unix))]
        {
            // Windows directory handles deny SHARE_DELETE. Close before the
            // nonrecursive removal; cap-relative cleanup cannot traverse outside.
            drop(self.directory.take());
            if owns_name {
                let _ = self.root.directory.remove_dir(&self.name);
            }
        }
    }
}
