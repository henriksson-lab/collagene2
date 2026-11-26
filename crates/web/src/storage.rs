//! Workspace storage is owned by the native repository in Tauri and IndexedDB
//! in a browser. UI components do not depend on either database implementation.
#[path = "storage_browser.rs"]
mod browser;
use crate::desktop_rpc;
use collagene_desktop_protocol as wire;
use futures_util::lock::Mutex;
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};

#[derive(Clone)]
pub enum Database {
    Browser(web_sys::IdbDatabase),
    Native {
        state: Rc<Mutex<Native>>,
        context: Option<wire::library::LabContext>,
    },
}
impl Database {
    pub fn native_context(&self) -> Option<wire::library::LabContext> {
        match self {
            Self::Browser(_) => None,
            Self::Native { context, .. } => context.clone(),
        }
    }

    pub fn close(&self) {
        if let Self::Browser(db) = self {
            db.close();
        }
    }
}
pub struct Native {
    working: Option<wire::library::Digest>,
    preferences: Option<wire::library::Digest>,
}
#[derive(Clone)]
pub enum Transaction {
    Browser(web_sys::IdbTransaction),
    Native {
        operation: wire::OperationId,
        cancelled: Rc<Cell<bool>>,
    },
}
impl Transaction {
    pub fn abort(&self) -> Result<(), String> {
        match self {
            Self::Browser(tx) => tx
                .abort()
                .map_err(|_| "Cannot abort browser transaction".into()),
            Self::Native {
                operation,
                cancelled,
            } => {
                cancelled.set(true);
                let operation_id = operation.clone();
                wasm_bindgen_futures::spawn_local(async move {
                    let _ = desktop_rpc::call(wire::Command::Cancel { operation_id }).await;
                });
                Ok(())
            }
        }
    }
}
pub async fn open() -> Result<Database, String> {
    if desktop_rpc::is_desktop() {
        match desktop_rpc::call(wire::Command::BootstrapRecovery).await? {
            wire::Success::BootstrapRecovery(b) => Ok(Database::Native {
                state: Rc::new(Mutex::new(Native {
                    working: b.working.version,
                    preferences: b.preferences.version,
                })),
                context: b.library_context,
            }),
            _ => Err("Unexpected local repository bootstrap".into()),
        }
    } else {
        browser::open().await.map(Database::Browser)
    }
}
pub async fn load(db: &Database) -> Result<Option<Vec<u8>>, String> {
    match db {
        Database::Browser(db) => browser::load(db).await,
        Database::Native { state, .. } => {
            let mut state = state.lock().await;
            let meta = match desktop_rpc::call(wire::Command::WorkingMetadata).await? {
                wire::Success::WorkingMetadata(m) => m,
                _ => return Err("Unexpected workspace metadata".into()),
            };
            if meta.current.is_none() {
                state.working = meta.version;
                return Ok(None);
            }
            match desktop_rpc::call(wire::Command::ReadWorking {
                slot: wire::WorkingSlot::Current,
                expected_version: meta.version,
            })
            .await?
            {
                wire::Success::WorkingArchive {
                    archive, version, ..
                } => {
                    state.working = Some(version);
                    Ok(Some(archive.to_vec()))
                }
                _ => Err("Unexpected local workspace".into()),
            }
        }
    }
}
pub async fn load_previous(db: &Database) -> Result<Option<Vec<u8>>, String> {
    match db {
        Database::Browser(db) => browser::load_previous(db).await,
        Database::Native { state, .. } => {
            let state = state.lock().await;
            // Do not advance this window's write baseline when reading recovery.
            let meta = match desktop_rpc::call(wire::Command::WorkingMetadata).await? {
                wire::Success::WorkingMetadata(m) => m,
                _ => return Err("Unexpected workspace metadata".into()),
            };
            if meta.previous.is_none() {
                return Ok(None);
            }
            if meta.version != state.working {
                return Err("Workspace changed in another window; reload before recovery".into());
            }
            match desktop_rpc::call(wire::Command::ReadWorking {
                slot: wire::WorkingSlot::Previous,
                expected_version: state.working.clone(),
            })
            .await?
            {
                wire::Success::WorkingArchive { archive, .. } => Ok(Some(archive.to_vec())),
                _ => Err("Unexpected workspace recovery copy".into()),
            }
        }
    }
}
pub async fn save(
    db: &Database,
    bytes: &[u8],
    guard: impl Fn() -> Result<(), String>,
) -> Result<(), String> {
    match db {
        Database::Browser(db) => {
            guard()?;
            browser::save(db, bytes).await
        }
        Database::Native { state, .. } => {
            let mut state = state.lock().await;
            let archive = wire::library::Blob::new(bytes.to_vec()).map_err(|e| e.to_string())?;
            let request = desktop_rpc::request(wire::Command::SaveCurrent {
                expected_version: state.working.clone(),
                archive,
            })?;
            match desktop_rpc::dispatch_guarded(&request, guard).await? {
                wire::Outcome::Success(wire::Success::WorkingMetadata(meta)) => {
                    state.working = meta.version;
                    Ok(())
                }
                wire::Outcome::Failure(failure) => Err(desktop_rpc::failure(&failure)),
                _ => Err("Unexpected local save response".into()),
            }
        }
    }
}
pub async fn replace_workspace(
    db: &Database,
    current: &[u8],
    previous: &[u8],
) -> Result<(), String> {
    replace_workspace_cancellable(db, current, previous, Rc::new(RefCell::new(None))).await
}
pub async fn replace_workspace_cancellable(
    db: &Database,
    current: &[u8],
    previous: &[u8],
    pending: Rc<RefCell<Option<Transaction>>>,
) -> Result<(), String> {
    match db {
        Database::Browser(db) => {
            browser::replace_workspace_cancellable(db, current, previous, pending).await
        }
        Database::Native { state, .. } => {
            let mut request = desktop_rpc::request(wire::Command::ReplaceWorking {
                expected_version: None,
                current: wire::library::Blob::new(current.to_vec()).map_err(|e| e.to_string())?,
                previous: wire::library::Blob::new(previous.to_vec()).map_err(|e| e.to_string())?,
            })?;
            let cancelled = Rc::new(Cell::new(false));
            *pending.borrow_mut() = Some(Transaction::Native {
                operation: request.operation_id.clone(),
                cancelled: cancelled.clone(),
            });
            let mut state = state.lock().await;
            if cancelled.get() {
                pending.borrow_mut().take();
                return Err("Workspace replacement cancelled".into());
            }
            if let wire::Command::ReplaceWorking {
                expected_version, ..
            } = &mut request.command
            {
                *expected_version = state.working.clone();
            }
            let result = desktop_rpc::dispatch(&request).await;
            pending.borrow_mut().take();
            match result? {
                // A committed result remains authoritative even if cancellation
                // arrived after the host's commit point.
                wire::Outcome::Success(wire::Success::WorkingMetadata(meta)) => {
                    state.working = meta.version;
                    Ok(())
                }
                wire::Outcome::Failure(failure) => Err(desktop_rpc::failure(&failure)),
                _ => Err("Unexpected workspace replacement response".into()),
            }
        }
    }
}
pub async fn load_view_preferences(db: &Database) -> Result<Option<Vec<u8>>, String> {
    match db {
        Database::Browser(db) => browser::load_view_preferences(db).await,
        Database::Native { state, .. } => {
            let mut state = state.lock().await;
            match desktop_rpc::call(wire::Command::ReadPreferences).await? {
                wire::Success::Preferences(p) => {
                    state.preferences = p.version.clone();
                    Ok(p.version.map(|_| p.bytes.to_vec()))
                }
                _ => Err("Unexpected display preferences".into()),
            }
        }
    }
}
pub async fn save_view_preferences(db: &Database, bytes: &[u8]) -> Result<(), String> {
    match db {
        Database::Browser(db) => browser::save_view_preferences(db, bytes).await,
        Database::Native { state, .. } => {
            let mut state = state.lock().await;
            match desktop_rpc::call(wire::Command::SavePreferences {
                expected_version: state.preferences.clone(),
                bytes: wire::library::Blob::new(bytes.to_vec()).map_err(|e| e.to_string())?,
            })
            .await?
            {
                wire::Success::Preferences(p) => {
                    state.preferences = p.version;
                    Ok(())
                }
                _ => Err("Unexpected preference save response".into()),
            }
        }
    }
}
