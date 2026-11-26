use crate::{HostError, HostResult};
use collagene_desktop_protocol as wire;
use collagene_library_protocol::Digest;
use std::{
    collections::BTreeMap,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, AtomicU64, Ordering},
    },
};
const MAX_ENTRIES: usize = 2048;
const MAX_CACHED_BYTES: usize = 4 << 20;
/// Opaque shell capability. Native creation is not authentication of arbitrary
/// web content; the shell grants it only to its approved local application window.
#[derive(Clone)]
pub struct WindowAccess {
    brand: Arc<()>,
    id: u64,
    alive: Arc<AtomicBool>,
}
struct Entry {
    state: wire::OperationState,
    seen_job: bool,
    digest: Option<Digest>,
    reply: Option<Vec<u8>>,
    order: u64,
}
#[derive(Default)]
struct Registry {
    entries: BTreeMap<(u64, String), Entry>,
    order: u64,
}
struct Inner {
    brand: Arc<()>,
    next_window: AtomicU64,
    registry: Mutex<Registry>,
}
#[derive(Clone)]
pub struct HostControl {
    inner: Arc<Inner>,
}
pub(crate) enum Begin {
    Run,
    Cached(Vec<u8>),
}
impl HostControl {
    pub(crate) fn new() -> Self {
        Self {
            inner: Arc::new(Inner {
                brand: Arc::new(()),
                next_window: AtomicU64::new(1),
                registry: Mutex::new(Registry::default()),
            }),
        }
    }
    pub(crate) fn authorize(&self) -> HostResult<WindowAccess> {
        let id = self.inner.next_window.fetch_add(1, Ordering::Relaxed);
        if id == 0 || id == u64::MAX {
            return Err(HostError::Unavailable);
        }
        Ok(WindowAccess {
            brand: self.inner.brand.clone(),
            id,
            alive: Arc::new(AtomicBool::new(true)),
        })
    }
    pub(crate) fn check(&self, access: &WindowAccess) -> HostResult<()> {
        if !Arc::ptr_eq(&access.brand, &self.inner.brand) || !access.alive.load(Ordering::Acquire) {
            return Err(HostError::Window);
        }
        Ok(())
    }
    pub fn revoke_window(&self, access: &WindowAccess) -> HostResult<()> {
        self.check(access)?;
        access.alive.store(false, Ordering::Release);
        let mut registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        registry
            .entries
            .retain(|(window, _), _| *window != access.id);
        Ok(())
    }
    fn room(registry: &mut Registry) -> HostResult<()> {
        if registry.entries.len() < MAX_ENTRIES {
            return Ok(());
        }
        let oldest = registry
            .entries
            .iter()
            .filter(|(_, entry)| {
                entry.seen_job
                    && matches!(
                        entry.state,
                        wire::OperationState::Completed
                            | wire::OperationState::Cancelled
                            | wire::OperationState::Failed
                    )
            })
            .min_by_key(|(_, e)| e.order)
            .map(|(key, _)| key.clone());
        if let Some(key) = oldest {
            registry.entries.remove(&key);
            Ok(())
        } else {
            Err(HostError::Unavailable)
        }
    }
    pub(crate) fn begin(
        &self,
        access: &WindowAccess,
        id: &wire::OperationId,
        digest: Digest,
    ) -> HostResult<Begin> {
        self.check(access)?;
        let mut registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        let key = (access.id, id.as_str().into());
        if let Some(entry) = registry.entries.get_mut(&key) {
            if matches!(
                entry.state,
                wire::OperationState::Cancelled | wire::OperationState::CancelRequested
            ) {
                entry.seen_job = true;
                return Err(HostError::Cancelled);
            }
            if entry.digest.as_ref() != Some(&digest) {
                return Err(HostError::Invalid(
                    "Operation ID was already used for another request.",
                ));
            }
            if let Some(reply) = &entry.reply {
                return Ok(Begin::Cached(reply.clone()));
            }
            return Err(HostError::Invalid(
                "Operation already started; inspect current state before retrying.",
            ));
        }
        Self::room(&mut registry)?;
        registry.order = registry.order.wrapping_add(1);
        let order = registry.order;
        registry.entries.insert(
            key,
            Entry {
                state: wire::OperationState::Running,
                seen_job: true,
                digest: Some(digest),
                reply: None,
                order,
            },
        );
        Ok(Begin::Run)
    }
    pub(crate) fn running(&self, access: &WindowAccess, id: &wire::OperationId) -> HostResult<()> {
        self.check(access)?;
        let registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        match registry
            .entries
            .get(&(access.id, id.as_str().into()))
            .map(|e| e.state)
        {
            Some(wire::OperationState::Running | wire::OperationState::Committing) => Ok(()),
            Some(wire::OperationState::Cancelled | wire::OperationState::CancelRequested) => {
                Err(HostError::Cancelled)
            }
            _ => Err(HostError::Unavailable),
        }
    }
    /// The single cancellation linearization point. SQLite work before this
    /// gate remains rollback-able; after acceptance cancellation cannot undo it.
    pub(crate) fn committing(
        &self,
        access: &WindowAccess,
        id: &wire::OperationId,
    ) -> HostResult<()> {
        self.check(access)?;
        let mut registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        let entry = registry
            .entries
            .get_mut(&(access.id, id.as_str().into()))
            .ok_or(HostError::Unavailable)?;
        match entry.state {
            wire::OperationState::Running => {
                entry.state = wire::OperationState::Committing;
                Ok(())
            }
            wire::OperationState::Cancelled | wire::OperationState::CancelRequested => {
                Err(HostError::Cancelled)
            }
            _ => Err(HostError::Unavailable),
        }
    }
    pub(crate) fn finish(
        &self,
        access: &WindowAccess,
        id: &wire::OperationId,
        state: wire::OperationState,
        reply: &[u8],
    ) -> HostResult<()> {
        let mut registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        let cached = registry
            .entries
            .values()
            .filter_map(|e| e.reply.as_ref())
            .map(Vec::len)
            .sum::<usize>();
        if let Some(entry) = registry.entries.get_mut(&(access.id, id.as_str().into())) {
            entry.state = state;
            if reply.len() <= 32768 && cached + reply.len() <= MAX_CACHED_BYTES {
                entry.reply = Some(reply.to_vec());
            }
        }
        Ok(())
    }
    pub fn dispatch_raw(&self, access: &WindowAccess, bytes: &[u8]) -> HostResult<Option<Vec<u8>>> {
        self.check(access)?;
        if bytes.len() > wire::MAX_CONTROL_BYTES {
            return Ok(None);
        }
        let Some(request) = wire::decode_control(bytes).map_err(|_| HostError::Wire)? else {
            return Ok(None);
        };
        let mut registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        let (target, cancel) = match &request.control {
            wire::Control::Cancel { operation_id } => (operation_id, true),
            wire::Control::OperationStatus { operation_id } => (operation_id, false),
        };
        let key = (access.id, target.as_str().into());
        if cancel && !registry.entries.contains_key(&key) {
            Self::room(&mut registry)?;
            registry.order = registry.order.wrapping_add(1);
            let order = registry.order;
            registry.entries.insert(
                key.clone(),
                Entry {
                    state: wire::OperationState::Cancelled,
                    seen_job: false,
                    digest: None,
                    reply: None,
                    order,
                },
            );
        }
        let state = if let Some(entry) = registry.entries.get_mut(&key) {
            if cancel
                && matches!(
                    entry.state,
                    wire::OperationState::Queued | wire::OperationState::Running
                )
            {
                entry.state = wire::OperationState::CancelRequested;
            }
            entry.state
        } else {
            wire::OperationState::Unknown
        };
        let success = if cancel {
            wire::Success::Cancelled {
                operation_id: target.clone(),
                status: state,
            }
        } else {
            wire::Success::Status {
                operation_id: target.clone(),
                status: state,
            }
        };
        let response = wire::ResponseEnvelope {
            request_id: request.request_id,
            operation_id: request.operation_id,
            outcome: wire::Outcome::Success(success),
        };
        Ok(Some(
            wire::encode_response(&response).map_err(|_| HostError::Wire)?,
        ))
    }
}
