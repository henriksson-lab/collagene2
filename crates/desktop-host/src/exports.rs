//! Operator-rooted, no-overwrite exports. No frontend path is ever opened.
use crate::{HostControl, HostError, HostResult, JobStart, WindowAccess};
use cap_fs_ext::{DirExt, FollowSymlinks, OpenOptionsFollowExt};
use cap_std::{
    ambient_authority,
    fs::{Dir, DirBuilder, OpenOptions},
};
use collagene_desktop_protocol as wire;
use serde::{Deserialize, Serialize};
use std::{
    io::{Read, Write},
    path::Path,
    sync::atomic::{AtomicU64, Ordering},
};
use wire::library::Digest;

pub struct DesktopExports {
    directory: Dir,
    operations: Dir,
    control: HostControl,
}
#[derive(Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
struct Intent {
    version: u8,
    suggested_name: String,
    name: String,
    hash: Digest,
    bytes: u64,
}
fn unavailable<T>(_: T) -> HostError {
    HostError::Unavailable
}
fn private_file(file: &std::fs::File) -> HostResult<()> {
    let metadata = file.metadata().map_err(unavailable)?;
    if !metadata.is_file() {
        return Err(HostError::Invalid(
            "Export storage requires regular private files.",
        ));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::{MetadataExt, PermissionsExt};
        if metadata.uid() != unsafe { libc::geteuid() }
            || metadata.permissions().mode() & 0o777 != 0o600
        {
            return Err(HostError::Invalid(
                "Export storage requires regular private files.",
            ));
        }
    }
    Ok(())
}
fn sync_directory(dir: &Dir) -> HostResult<()> {
    // Windows has no portable directory FlushFileBuffers guarantee. File data
    // is flushed there; the documented namespace power-loss limit remains.
    #[cfg(unix)]
    {
        // A directory capability may use O_PATH on Linux. Open a real readable
        // descriptor relative to it; O_PATH descriptors cannot be fsynced.
        dir.open(".")
            .map_err(unavailable)?
            .sync_all()
            .map_err(unavailable)?;
    }
    #[cfg(not(unix))]
    {
        let _ = dir;
    }
    Ok(())
}
fn private_directory(parent: &Dir, name: &str) -> HostResult<Dir> {
    let mut builder = DirBuilder::new();
    #[cfg(unix)]
    {
        use cap_std::fs::DirBuilderExt;
        builder.mode(0o700);
    }
    match parent.create_dir_with(name, &builder) {
        Ok(()) => sync_directory(parent)?,
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
        Err(e) => return Err(unavailable(e)),
    }
    let dir = parent.open_dir_nofollow(name).map_err(unavailable)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::{MetadataExt, PermissionsExt};
        let metadata = dir
            .try_clone()
            .map_err(unavailable)?
            .into_std_file()
            .metadata()
            .map_err(unavailable)?;
        if metadata.uid() != unsafe { libc::geteuid() }
            || metadata.permissions().mode() & 0o777 != 0o700
        {
            return Err(HostError::Invalid(
                "Export directories must remain private.",
            ));
        }
    }
    sync_directory(parent)?;
    Ok(dir)
}
fn read_private(dir: &Dir, name: &str, limit: usize) -> HostResult<Vec<u8>> {
    let mut options = OpenOptions::new();
    options.read(true).follow(FollowSymlinks::No);
    let mut file = dir
        .open_with(name, &options)
        .map_err(unavailable)?
        .into_std();
    private_file(&file)?;
    if file.metadata().map_err(unavailable)?.len() > limit as u64 {
        return Err(HostError::Invalid(
            "Export receipt or payload exceeds its bound.",
        ));
    }
    read_file(&mut file, limit)
}
fn read_file(file: &mut std::fs::File, limit: usize) -> HostResult<Vec<u8>> {
    private_file(file)?;
    if file.metadata().map_err(unavailable)?.len() > limit as u64 {
        return Err(HostError::Invalid(
            "Export receipt or payload exceeds its bound.",
        ));
    }
    let mut bytes = Vec::new();
    Read::by_ref(file)
        .take(limit as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(unavailable)?;
    if bytes.len() > limit {
        return Err(HostError::Invalid(
            "Export receipt or payload exceeds its bound.",
        ));
    }
    Ok(bytes)
}
struct Temporary<'a> {
    dir: &'a Dir,
    name: String,
    file: std::fs::File,
}
impl Temporary<'_> {
    fn publish(&self, destination: &Dir, name: &str) -> std::io::Result<()> {
        self.dir.hard_link(&self.name, destination, name)
    }
}
impl Drop for Temporary<'_> {
    fn drop(&mut self) {
        #[cfg(unix)]
        {
            use cap_std::fs::MetadataExt as CapMetadataExt;
            use std::os::unix::fs::MetadataExt as StdMetadataExt;
            if let (Ok(owned), Ok(current)) =
                (self.file.metadata(), self.dir.symlink_metadata(&self.name))
                && owned.dev() == current.dev()
                && owned.ino() == current.ino()
            {
                // Retain the owned handle through identity check and unlink.
                // This protects observed replacement; it is not a same-user
                // concurrent rename sandbox or an atomic conditional unlink.
                let _ = self.dir.remove_file(&self.name);
            }
        }
        #[cfg(windows)]
        {
            use std::os::windows::io::AsRawHandle;
            use windows_sys::Win32::Storage::FileSystem::{
                FILE_DISPOSITION_INFO, FileDispositionInfo, SetFileInformationByHandle,
            };
            let info = FILE_DISPOSITION_INFO { DeleteFile: true };
            // Mark the retained creation handle, never reopen a pathname. DELETE
            // access was requested at creation; share DELETE allows hard-link
            // publication and does not redirect this owned cleanup handle. Do this
            // only now: delete-pending temps cannot be hard-link published.
            unsafe {
                let _ = SetFileInformationByHandle(
                    self.file.as_raw_handle(),
                    FileDispositionInfo,
                    (&info as *const FILE_DISPOSITION_INFO).cast(),
                    std::mem::size_of::<FILE_DISPOSITION_INFO>() as u32,
                );
            }
        }
    }
}
fn temporary<'a>(dir: &'a Dir, bytes: &[u8]) -> HostResult<Temporary<'a>> {
    static NEXT: AtomicU64 = AtomicU64::new(1);
    for _ in 0..128 {
        let name = format!(
            ".partial-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        );
        let mut options = OpenOptions::new();
        options
            .read(true)
            .write(true)
            .create_new(true)
            .follow(FollowSymlinks::No);
        #[cfg(unix)]
        {
            use cap_std::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        #[cfg(windows)]
        {
            use cap_std::fs::OpenOptionsExt;
            options
                .access_mode(0x8000_0000 | 0x4000_0000 | 0x0001_0000)
                .share_mode(0x1 | 0x2 | 0x4); // READ/WRITE/DELETE sharing for publication
        }
        match dir.open_with(&name, &options) {
            Ok(file) => {
                let mut temporary = Temporary {
                    dir,
                    name,
                    file: file.into_std(),
                };
                private_file(&temporary.file)?;
                temporary.file.write_all(bytes).map_err(unavailable)?;
                temporary.file.sync_all().map_err(unavailable)?;
                return Ok(temporary);
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(unavailable(e)),
        }
    }
    Err(HostError::Unavailable)
}
impl DesktopExports {
    pub(crate) fn open(root: &Path, control: HostControl) -> HostResult<Self> {
        let root = crate::permissions::root(root)?;
        let root = Dir::open_ambient_dir(root, ambient_authority()).map_err(unavailable)?;
        let directory = private_directory(&root, "exports")?;
        let operations = private_directory(&directory, ".operations")?;
        Ok(Self {
            directory,
            operations,
            control,
        })
    }
    fn write(
        &self,
        op: &wire::OperationId,
        suggested: &str,
        bytes: &[u8],
    ) -> HostResult<wire::Success> {
        let name = wire::export_file_name(op, suggested).map_err(|_| HostError::Wire)?;
        let intent = Intent {
            version: 1,
            suggested_name: suggested.into(),
            name: name.to_string(),
            hash: Digest::of(bytes),
            bytes: bytes.len() as u64,
        };
        let operation = private_directory(&self.operations, op.as_str())?;
        let encoded = serde_json::to_vec(&intent).map_err(unavailable)?;
        // An immutable intent is durable before publishing any output. Changing
        // bytes OR suggested name under the same op ID fails even after restart.
        let reservation = temporary(&operation, &encoded)?;
        match reservation.publish(&operation, "intent.json") {
            Ok(()) => sync_directory(&operation)?,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(e) => return Err(unavailable(e)),
        }
        let existing = read_private(&operation, "intent.json", 4096)?;
        let existing: Intent = serde_json::from_slice(&existing).map_err(|_| {
            HostError::Invalid("Export intent receipt is invalid; preserve existing files.")
        })?;
        if existing != intent {
            return Err(HostError::Invalid(
                "Export operation intent changed; preserve the original export.",
            ));
        }
        sync_directory(&operation)?;
        drop(reservation);
        let payload = temporary(&operation, bytes)?;
        // Atomic no-overwrite publication, even across two native processes.
        match payload.publish(&self.directory, name.as_str()) {
            Ok(()) => (),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(e) => return Err(unavailable(e)),
        }
        // Verify and flush the same opened no-follow file, retaining its handle
        // throughout; do not validate one pathname lookup and flush another.
        let mut options = OpenOptions::new();
        options.read(true).write(true).follow(FollowSymlinks::No);
        let mut file = self
            .directory
            .open_with(name.as_str(), &options)
            .map_err(unavailable)?
            .into_std();
        let published = read_file(&mut file, wire::MAX_EXPORT_BYTES)?;
        if published.len() != bytes.len() || Digest::of(&published) != intent.hash {
            return Err(HostError::Invalid(
                "An existing export does not match its receipt; no file was overwritten.",
            ));
        }
        file.sync_all().map_err(unavailable)?;
        drop(file);
        sync_directory(&self.directory)?;
        drop(payload);
        Ok(wire::Success::Exported {
            name,
            hash: intent.hash,
            bytes: intent.bytes,
        })
    }
    pub fn dispatch_raw(&self, access: &WindowAccess, bytes: &[u8]) -> HostResult<Vec<u8>> {
        self.control.check_window(access)?;
        let request = wire::decode_request(bytes).map_err(|_| HostError::Wire)?;
        let wire::Command::ExportFile {
            suggested_name,
            bytes: payload,
        } = &request.command
        else {
            return Err(HostError::Wire);
        };
        let job = match self
            .control
            .begin_job(access, &request.operation_id, Digest::of(bytes))
        {
            Err(error) => {
                let response = wire::ResponseEnvelope {
                    request_id: request.request_id,
                    operation_id: request.operation_id.clone(),
                    outcome: wire::Outcome::Failure(crate::failure(error)),
                };
                response
                    .validate_for(&request)
                    .map_err(|_| HostError::Wire)?;
                return wire::encode_response_for_request(
                    &response,
                    &request,
                    u16::from_le_bytes([bytes[4], bytes[5]]),
                )
                .map_err(|_| HostError::Wire);
            }

            Ok(JobStart::Cached(reply)) => {
                wire::decode_response_for(&reply, &request).map_err(|_| HostError::Wire)?;
                return Ok(reply);
            }
            Ok(JobStart::Run(job)) => job,
        };
        let result = job
            .accept()
            .and_then(|()| self.write(&request.operation_id, suggested_name, payload));
        let outcome = match result {
            Ok(value) => wire::Outcome::Success(value),
            Err(error) => wire::Outcome::Failure(crate::failure(error)),
        };
        let terminal = match &outcome {
            wire::Outcome::Success(_) => wire::OperationState::Completed,
            wire::Outcome::Failure(wire::Failure::Cancelled) => wire::OperationState::Cancelled,
            _ => wire::OperationState::Failed,
        };
        let response = wire::ResponseEnvelope {
            request_id: request.request_id,
            operation_id: request.operation_id.clone(),
            outcome,
        };
        response
            .validate_for(&request)
            .map_err(|_| HostError::Wire)?;
        let reply = wire::encode_response_for_request(
            &response,
            &request,
            u16::from_le_bytes([bytes[4], bytes[5]]),
        )
        .map_err(|_| HostError::Wire)?;
        // Durable bytes and their immutable intent are the outcome authority.
        // A cache bookkeeping failure after commit must not turn success into a
        // claimed cancellation/failure; a retry revalidates the same owned file.
        let _ = job.finish(terminal, &reply);
        Ok(reply)
    }
}

#[cfg(test)]
mod cleanup_tests {
    use super::*;
    struct Disk(std::path::PathBuf);
    impl Drop for Disk {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    fn directory() -> (Disk, Dir) {
        static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(1);
        let path = std::env::temp_dir().join(format!(
            "collagene-export-cleanup-{}-{}-{}",
            std::process::id(),
            NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let root = crate::permissions::root(&path).unwrap();
        let dir = Dir::open_ambient_dir(root, ambient_authority()).unwrap();
        (Disk(path), dir)
    }
    #[test]
    fn temporary_handle_cleanup_keeps_the_published_hardlink() {
        let (_disk, dir) = directory();
        let temp = temporary(&dir, b"exact final bytes").unwrap();
        let name = temp.name.clone();
        temp.publish(&dir, "published.bin").unwrap();
        drop(temp);
        assert!(dir.symlink_metadata(&name).is_err());
        assert_eq!(
            read_private(&dir, "published.bin", 64).unwrap(),
            b"exact final bytes"
        );
    }
    #[cfg(any(unix, windows))]
    #[test]
    fn observed_replacement_is_not_removed_by_temporary_cleanup() {
        let (_disk, dir) = directory();
        let temp = temporary(&dir, b"owned original").unwrap();
        let name = temp.name.clone();
        dir.rename(&name, &dir, "renamed-owned.bin").unwrap();
        let mut options = OpenOptions::new();
        options.write(true).create_new(true);
        #[cfg(unix)]
        {
            use cap_std::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        dir.open_with(&name, &options)
            .unwrap()
            .write_all(b"foreign replacement sentinel")
            .unwrap();
        drop(temp);
        assert_eq!(
            read_private(&dir, &name, 64).unwrap(),
            b"foreign replacement sentinel"
        );
        #[cfg(unix)]
        assert_eq!(
            read_private(&dir, "renamed-owned.bin", 64).unwrap(),
            b"owned original"
        );
        #[cfg(windows)]
        assert!(dir.symlink_metadata("renamed-owned.bin").is_err());
    }
}
