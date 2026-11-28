//! Framework-independent native host for one bounded raw desktop IPC command.
//! Shell/window trust and the private storage root are supplied natively, never
//! accepted from request data. Local operations need no network or login.
#![cfg(not(target_arch = "wasm32"))]
pub mod catalog_backup;
mod control;
mod exports;
pub use exports::DesktopExports;
mod library;
mod permissions;
mod working;
use collagene_desktop_protocol as wire;
pub use collagene_desktop_protocol::{MAX_REQUEST_BYTES, MAX_RESPONSE_BYTES, OperationState};
use collagene_lab_library as store;
use collagene_library_protocol as base;
pub use control::{HostControl, JobStart, ScopedJob, WindowAccess};
use serde::{Serialize, de::DeserializeOwned};
use std::{
    io::Write,
    path::{Path, PathBuf},
    sync::atomic::{AtomicBool, Ordering},
};
#[derive(Debug, thiserror::Error)]
pub enum HostError {
    #[error("desktop window access is not valid")]
    Window,
    #[error("invalid desktop request framing")]
    Wire,
    #[error("{0}")]
    Invalid(&'static str),
    #[error("native desktop persistence unavailable")]
    Unavailable,
    #[error("private working storage is corrupt; inspect previous recovery data")]
    Corrupt,
    #[error("working state changed in another window")]
    WorkingConflict(wire::WorkingMetadata),
    #[error("view preferences changed in another window")]
    PreferencesConflict(wire::Preferences),
    #[error("local library operation rejected")]
    Library(base::Failure),
    #[error("operation cancelled before committing")]
    Cancelled,
}
pub type HostResult<T> = std::result::Result<T, HostError>;

/// Trusted native operator entry point for a saved-library session. Unlike the
/// editor host this does not open or initialize working/recovery state.
/// No IPC/MCP request may supply the root passed here.
pub fn open_local_repository(
    configured_root: impl AsRef<Path>,
) -> HostResult<collagene_repository::Repository> {
    let root = permissions::root(configured_root.as_ref())?;
    if catalog_backup::catalog_activation_pending(&root)? {
        return Err(HostError::Invalid(
            "Restart the desktop before reading the accepted catalog recovery.",
        ));
    }
    let path = catalog_backup::current_catalog_path(&root)?;
    collagene_repository::Repository::open_read_only(path).map_err(|_| HostError::Unavailable)
}
impl From<rusqlite::Error> for HostError {
    fn from(_: rusqlite::Error) -> Self {
        Self::Unavailable
    }
}
struct LocalCatalog {
    path: PathBuf,
    generation: u64,
    repository: collagene_repository::Repository,
    reads: store::Library,
    local_access: store::LabAccess,
}
pub struct DesktopHost {
    catalog: Option<LocalCatalog>,
    working: working::Working,
    control: HostControl,
    _root: PathBuf,
    // Last field: release process lease only after every owned connection closes.
    _catalog_lease: std::sync::Arc<catalog_backup::CatalogLease>,
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > 4 << 20)
        {
            return Err(std::io::Error::other("desktop metadata bound"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn convert<T: Serialize, U: DeserializeOwned>(value: &T) -> HostResult<U> {
    let mut out = Bounded(vec![]);
    serde_json::to_writer(&mut out, value).map_err(|_| HostError::Unavailable)?;
    serde_json::from_slice(&out.0).map_err(|_| HostError::Unavailable)
}
fn failure(error: HostError) -> wire::Failure {
    match error{
    HostError::Invalid(message)=>wire::Failure::Invalid{message:base::Text::new(message).expect("static bounded desktop message")},
    HostError::WorkingConflict(meta)=>wire::Failure::Conflict{current:wire::Conflict::Working(meta)},
    HostError::PreferencesConflict(prefs)=>wire::Failure::Conflict{current:wire::Conflict::Preferences(prefs)},
    HostError::Cancelled=>wire::Failure::Cancelled,
    HostError::Library(error)=>wire::Failure::Library(error),
    HostError::Corrupt=>wire::Failure::Invalid{message:base::Text::new("Private working data failed native validation; previous recovery can be inspected separately.").expect("bounded static recovery message")},
    _=>wire::Failure::Unavailable,
}
}
fn commit_gate(
    control: &HostControl,
    access: &WindowAccess,
    op: &wire::OperationId,
    accepted: &AtomicBool,
) -> HostResult<()> {
    control.committing(access, op)?;
    accepted.store(true, Ordering::Release);
    Ok(())
}
fn execution_outcome(
    control: &HostControl,
    access: &WindowAccess,
    op: &wire::OperationId,
    accepted: bool,
    result: HostResult<wire::Success>,
) -> wire::Outcome {
    if !accepted && control.running(access, op).is_err() {
        wire::Outcome::Failure(wire::Failure::Cancelled)
    } else {
        match result {
            Ok(success) => wire::Outcome::Success(success),
            Err(error) => wire::Outcome::Failure(failure(error)),
        }
    }
}
impl DesktopHost {
    pub fn open(configured_root: impl AsRef<Path>) -> HostResult<Self> {
        // The shell supplies a configured private root; request data never does.
        let root = permissions::root(configured_root.as_ref())?;
        let catalog_lease = catalog_backup::CatalogLease::acquire(&root)?;
        // Working data is opened independently, before any optional catalog startup.
        let working_path = permissions::database(&root, "working.sqlite")?;
        let working = working::Working::open(&working_path)?;
        let catalog = (|| -> HostResult<LocalCatalog> {
            let path = catalog_backup::startup_catalog_path(&root, &catalog_lease)?;
            let mut repository = collagene_repository::Repository::open(&path)
                .map_err(|_| HostError::Unavailable)?;
            let mut reads = store::Library::open(&path).map_err(|_| HostError::Unavailable)?;
            repository
                .retain_native_owner(catalog_lease.clone())
                .map_err(|_| HostError::Unavailable)?;
            reads
                .retain_native_owner(catalog_lease.clone())
                .map_err(|_| HostError::Unavailable)?;
            let local_access = reads
                .local_repository_access("offline-plasmids", "Local plasmid repository")
                .map_err(|_| HostError::Unavailable)?;
            let generation = catalog_backup::catalog_generation(&root)?;
            Ok(LocalCatalog {
                path,
                generation,
                repository,
                reads,
                local_access,
            })
        })()
        .ok();
        Ok(Self {
            catalog,
            working,
            control: HostControl::new(),
            _root: root,
            _catalog_lease: catalog_lease,
        })
    }
    pub fn catalog_available(&self) -> bool {
        matches!(self.catalog_status(), wire::CatalogStatus::Ready { .. })
    }
    fn catalog_context(&self) -> HostResult<base::LabContext> {
        let c = self.catalog.as_ref().ok_or(HostError::Unavailable)?;
        catalog_backup::check_catalog_connection(&self._root, &c.path)?;
        let lab = c
            .reads
            .lab(&c.local_access)
            .map_err(|_| HostError::Unavailable)?;
        base::validate_name(&lab.name).map_err(|_| HostError::Unavailable)?;
        Ok(base::LabContext {
            lab_id: base::LabId::new(lab.id.as_str()).map_err(|_| HostError::Unavailable)?,
            name: base::Text::new(lab.name).map_err(|_| HostError::Unavailable)?,
        })
    }
    fn catalog_status(&self) -> wire::CatalogStatus {
        let status = catalog_backup::catalog_session_status(
            &self._root,
            self.catalog.as_ref().map(|c| c.generation),
        );
        if matches!(status, wire::CatalogStatus::Ready { .. }) && self.catalog_context().is_err() {
            return wire::CatalogStatus::Paused{message:base::Text::new("Saved catalog connection is unavailable. Editor recovery remains available; restart after resolving catalog integrity.").expect("bounded static catalog status"),activation:None};
        }
        status
    }

    pub fn authorize_window(&mut self) -> HostResult<WindowAccess> {
        self.control.authorize()
    }
    /// A second trusted native connection permits remote I/O to proceed without
    /// owning the working-file host mutex. IPC cannot choose this fixed path.
    pub fn open_catalog_backup(&self) -> HostResult<catalog_backup::DesktopCatalogBackups> {
        let lab = self
            .catalog
            .as_ref()
            .and_then(|c| c.reads.lab(&c.local_access).ok().map(|l| l.id));
        catalog_backup::DesktopCatalogBackups::open_recovery_with_lease(
            &self._root,
            lab.as_ref(),
            self.control.clone(),
            self._catalog_lease.clone(),
            self.catalog.as_ref().map(|c| c.generation),
        )
    }
    pub fn open_repository_connection(&self) -> HostResult<collagene_repository::Repository> {
        if !self.catalog_available() {
            return Err(HostError::Invalid(
                "Saved catalog is unavailable; editor recovery remains available.",
            ));
        }
        let catalog = self.catalog.as_ref().ok_or(HostError::Unavailable)?;
        catalog_backup::check_catalog_connection(&self._root, &catalog.path)?;
        let mut repository = collagene_repository::Repository::open(&catalog.path)
            .map_err(|_| HostError::Unavailable)?;
        repository
            .retain_native_owner(self._catalog_lease.clone())
            .map_err(|_| HostError::Unavailable)?;
        Ok(repository)
    }

    pub fn open_exports(&self) -> HostResult<DesktopExports> {
        DesktopExports::open(&self._root, self.control.clone())
    }
    pub fn control(&self) -> HostControl {
        self.control.clone()
    }
    fn execute(
        &mut self,
        access: &WindowAccess,
        request: &wire::RequestEnvelope,
        accepted: &AtomicBool,
    ) -> HostResult<wire::Success> {
        let op = &request.operation_id;
        let control = self.control.clone();
        use wire::{Command as C, Success as S};
        Ok(match &request.command {
            C::CatalogBackup(_) => {
                return Err(HostError::Invalid("Use the catalog backup service lane."));
            }
            C::Bootstrap => {
                let catalog = self.catalog_status();
                let library_context = if matches!(catalog, wire::CatalogStatus::Ready { .. }) {
                    self.catalog_context().ok()
                } else {
                    None
                };
                let catalog = if library_context.is_none()
                    && matches!(catalog, wire::CatalogStatus::Ready { .. })
                {
                    wire::CatalogStatus::Paused{message:base::Text::new("Saved catalog context is unavailable. Working editor recovery remains available.").expect("bounded static message"),activation:None}
                } else {
                    catalog
                };
                S::Bootstrap(wire::Bootstrap {
                    library_context,
                    working: self.working.metadata()?,
                    preferences: self.working.preferences()?,
                    catalog,
                })
            }

            C::Library(command) => {
                if !self.catalog_available() {
                    return Err(HostError::Invalid(
                        "Saved catalog is paused; editor recovery remains available.",
                    ));
                }
                self.local_library(access, op, command, accepted)
                    .map(S::Library)?
            }
            C::WorkingMetadata => S::WorkingMetadata(self.working.metadata()?),
            C::ReadWorking {
                slot,
                expected_version,
            } => {
                let (version, bytes) = self.working.read(*slot, expected_version)?;
                S::WorkingArchive {
                    slot: *slot,
                    version,
                    archive: base::Blob::new(bytes).map_err(|_| HostError::Unavailable)?,
                }
            }
            C::SaveCurrent {
                expected_version,
                archive,
            } => S::WorkingMetadata(self.working.save(expected_version, archive, None, || {
                commit_gate(&control, access, op, accepted)
            })?),
            C::ReplaceWorking {
                expected_version,
                current,
                previous,
            } => S::WorkingMetadata(self.working.save(
                expected_version,
                current,
                Some(previous.as_slice()),
                || commit_gate(&control, access, op, accepted),
            )?),
            C::RecoverPrevious { expected_version } => {
                S::WorkingMetadata(self.working.recover(expected_version, || {
                    commit_gate(&control, access, op, accepted)
                })?)
            }
            C::ReadPreferences => S::Preferences(self.working.preferences()?),
            C::SavePreferences {
                expected_version,
                bytes,
            } => S::Preferences(self.working.save_preferences(expected_version, bytes, || {
                commit_gate(&control, access, op, accepted)
            })?),
            C::Cancel { .. } | C::OperationStatus { .. } | C::ExportFile { .. } => {
                return Err(HostError::Wire);
            }
        })
    }
    pub fn dispatch_raw(&mut self, access: &WindowAccess, bytes: &[u8]) -> HostResult<Vec<u8>> {
        self.control.check(access)?;
        if let Some(reply) = self.control.dispatch_raw(access, bytes)? {
            return Ok(reply);
        }
        if bytes.len() > MAX_REQUEST_BYTES {
            return Err(HostError::Wire);
        }
        let request = wire::decode_request(bytes).map_err(|_| HostError::Wire)?;
        if matches!(request.command, wire::Command::ExportFile { .. }) {
            return Err(HostError::Wire);
        }
        // Frame checksum already validated; use its digest without hashing or
        // allocating the complete128MiB replacement body again.
        let mut hash = String::with_capacity(64);
        for byte in &bytes[10..42] {
            use std::fmt::Write;
            let _ = write!(&mut hash, "{byte:02x}");
        }
        let digest = base::Digest::new(hash).map_err(|_| HostError::Wire)?;
        let begin = self.control.begin(access, &request.operation_id, digest);
        let started = matches!(&begin, Ok(control::Begin::Run));
        let outcome = match begin {
            Ok(control::Begin::Cached(reply)) => {
                if matches!(&request.command, wire::Command::Bootstrap) {
                    // Bootstrap is an observational read. A cached Ready/Pending
                    // receipt must not hide a later accepted/cancelled generation.
                    let result = self.execute(access, &request, &AtomicBool::new(false));
                    let response = wire::ResponseEnvelope {
                        request_id: request.request_id,
                        operation_id: request.operation_id.clone(),
                        outcome: match result {
                            Ok(s) => wire::Outcome::Success(s),
                            Err(e) => wire::Outcome::Failure(failure(e)),
                        },
                    };
                    return wire::encode_response_for_frame(&response, bytes)
                        .map_err(|_| HostError::Wire);
                }
                return Ok(reply);
            }
            Ok(control::Begin::Run) => {
                let accepted = AtomicBool::new(false);
                let result = self.execute(access, &request, &accepted);
                execution_outcome(
                    &self.control,
                    access,
                    &request.operation_id,
                    accepted.load(Ordering::Acquire),
                    result,
                )
            }
            Err(error) => wire::Outcome::Failure(failure(error)),
        };
        let state = match &outcome {
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
        let reply =
            wire::encode_response_for_request(&response, &request).map_err(|_| HostError::Wire)?;
        if started {
            self.control
                .finish(access, &request.operation_id, state, &reply)?;
        }
        Ok(reply)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_lab_library::mutations as durable;
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(1);
    fn op(n: u64) -> wire::OperationId {
        wire::OperationId::new(format!("op-{n:032x}")).unwrap()
    }
    fn cancel(
        control: &HostControl,
        access: &WindowAccess,
        target: wire::OperationId,
    ) -> wire::OperationState {
        let request = wire::RequestEnvelope {
            request_id: 1,
            operation_id: op(999),
            command: wire::Command::Cancel {
                operation_id: target,
            },
        };
        let frame = wire::encode_request(&request).unwrap();
        let response = control.dispatch_raw(access, &frame).unwrap().unwrap();
        match wire::decode_response_for(&response, &request)
            .unwrap()
            .outcome
        {
            wire::Outcome::Success(wire::Success::Cancelled { status, .. }) => status,
            other => panic!("{other:?}"),
        }
    }
    #[test]
    fn commit_gate_rolls_back_real_repository_effect_and_too_late_cancel_does_not_undo_working_write()
     {
        let root = std::env::temp_dir().join(format!(
            "collagene-native-gate-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut host = DesktopHost::open(&root).unwrap();
        let window = host.authorize_window().unwrap();
        let control = host.control();
        let id = op(1);
        control
            .begin(&window, &id, base::Digest::of(b"real local create"))
            .unwrap();
        let result = host
            .catalog
            .as_mut()
            .unwrap()
            .repository
            .mutate_local_outcome_with_gate(
                durable::MutationOperation::CreateFolder {
                    parent: None,
                    name: "Must roll back".into(),
                },
                || {
                    assert_eq!(
                        cancel(&control, &window, id.clone()),
                        wire::OperationState::CancelRequested
                    );
                    control.committing(&window, &id).is_ok()
                },
            );
        assert!(matches!(
            result,
            Err(collagene_repository::Error::Cancelled)
        ));
        let catalog = host.catalog.as_mut().unwrap();
        assert!(
            catalog
                .reads
                .browse(&catalog.local_access, None, 0, 20, None)
                .unwrap()
                .entries
                .is_empty()
        );
        let molecule = collagene_model::Molecule::dna(
            "Working data",
            "ACGTACGT",
            collagene_model::Topology::Linear,
        )
        .unwrap();
        let frame = collagene_core::documents::workspace::Workspace::from_documents(
            "Private workspace",
            vec![molecule],
            0,
        )
        .unwrap();
        let archive = collagene_core::documents::workspace::save_workspace(&frame).unwrap();
        let id = op(2);
        control
            .begin(&window, &id, base::Digest::of(b"real working save"))
            .unwrap();
        let saved = host
            .working
            .save(&None, &archive, None, || {
                control.committing(&window, &id)?;
                assert_eq!(
                    cancel(&control, &window, id.clone()),
                    wire::OperationState::Committing
                );
                Ok(())
            })
            .unwrap();
        assert!(saved.version.is_some());
        assert_eq!(
            host.working
                .read(wire::WorkingSlot::Current, &saved.version)
                .unwrap()
                .1,
            archive
        );
        drop(host);
        let reopened = DesktopHost::open(&root).unwrap();
        assert_eq!(reopened.working.metadata().unwrap(), saved);
        drop(reopened);
        let _ = std::fs::remove_dir_all(root);
    }
    #[test]
    fn committed_result_survives_window_revocation_and_replay_never_reruns_cancel_gate() {
        let root = std::env::temp_dir().join(format!(
            "collagene-native-revocation-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut host = DesktopHost::open(&root).unwrap();
        let window = host.authorize_window().unwrap();
        let control = host.control();
        let id = op(50);
        control
            .begin(&window, &id, base::Digest::of(b"commit then revoke"))
            .unwrap();
        let accepted = AtomicBool::new(false);
        let frame = collagene_core::documents::workspace::Workspace::from_documents(
            "Committed private workspace",
            vec![
                collagene_model::Molecule::dna(
                    "DNA",
                    "ACGTACGT",
                    collagene_model::Topology::Linear,
                )
                .unwrap(),
            ],
            0,
        )
        .unwrap();
        let archive = collagene_core::documents::workspace::save_workspace(&frame).unwrap();
        let saved = host
            .working
            .save(&None, &archive, None, || {
                commit_gate(&control, &window, &id, &accepted)?;
                control.revoke_window(&window)?;
                Ok(())
            })
            .unwrap();
        assert!(
            matches!(execution_outcome(&control,&window,&id,accepted.load(Ordering::Acquire),Ok(wire::Success::WorkingMetadata(saved.clone()))),wire::Outcome::Success(wire::Success::WorkingMetadata(m)) if m==saved)
        );
        assert!(matches!(
            execution_outcome(&control, &window, &id, true, Err(HostError::Unavailable)),
            wire::Outcome::Failure(wire::Failure::Unavailable)
        ));
        let mutation = durable::Mutation {
            id: durable::MutationId::parse("mut-aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa").unwrap(),
            operation: durable::MutationOperation::CreateFolder {
                parent: None,
                name: "Stable operation".into(),
            },
        };
        assert!(matches!(
            host.catalog
                .as_mut()
                .unwrap()
                .repository
                .mutate_local_reliable_with_gate(mutation.clone(), || false),
            Err(collagene_repository::Error::Cancelled)
        ));
        let first = host
            .catalog
            .as_mut()
            .unwrap()
            .repository
            .mutate_local_reliable_with_gate(mutation.clone(), || true)
            .unwrap();
        let replay = host
            .catalog
            .as_mut()
            .unwrap()
            .repository
            .mutate_local_reliable_with_gate(mutation, || panic!("receipt replay reran callback"))
            .unwrap();
        assert_eq!(first, replay);
        drop(host);
        let reopened = DesktopHost::open(&root).unwrap();
        assert_eq!(reopened.working.metadata().unwrap(), saved);
        drop(reopened);
        let _ = std::fs::remove_dir_all(root);
    }
    #[test]
    fn shell_state_can_move_to_blocking_worker_but_window_control_is_shareable() {
        fn send<T: Send>() {}
        fn sync<T: Sync>() {}
        send::<DesktopHost>();
        send::<WindowAccess>();
        sync::<WindowAccess>();
        send::<HostControl>();
        sync::<HostControl>();
    }
}
