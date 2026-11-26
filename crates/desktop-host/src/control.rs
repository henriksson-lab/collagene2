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
/// Trusted native caller has already validated the external protocol frame.
/// Replies cached here must be checked, secret-free protocol responses.
pub enum JobStart {
    Cached(Vec<u8>),
    Run(ScopedJob),
}
/// A branded window/op/digest lease, not raw repository or filesystem authority.
/// It is deliberately neither Clone, Debug nor serializable.
pub struct ScopedJob {
    control: HostControl,
    access: WindowAccess,
    operation: wire::OperationId,
    digest: Digest,
    accepted: AtomicBool,
    finished: bool,
}
impl ScopedJob {
    /// Linearizes acceptance with cancellation/revocation under the registry lock.
    /// Repeated acceptance of this same already-accepted guard is harmless.
    pub fn accept(&self) -> HostResult<()> {
        if self.accepted() {
            return Ok(());
        }
        let mut registry = self
            .control
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        if self.accepted() {
            return Ok(());
        }
        self.control.check(&self.access)?;
        let entry = registry
            .entries
            .get_mut(&(self.access.id, self.operation.as_str().into()))
            .ok_or(HostError::Unavailable)?;
        if entry.digest.as_ref() != Some(&self.digest) {
            return Err(HostError::Invalid("Operation intent changed."));
        }
        match entry.state {
            wire::OperationState::Running => {
                entry.state = wire::OperationState::Committing;
                self.accepted.store(true, Ordering::Release);
                Ok(())
            }
            wire::OperationState::Cancelled | wire::OperationState::CancelRequested => {
                Err(HostError::Cancelled)
            }
            _ => Err(HostError::Unavailable),
        }
    }
    pub fn accepted(&self) -> bool {
        self.accepted.load(Ordering::Acquire)
    }
    pub fn cancelled_before_acceptance(&self) -> bool {
        if self.accepted() {
            return false;
        }
        self.control.running(&self.access, &self.operation).is_err() && !self.accepted()
    }
    /// The caller validates its own protocol and excludes credentials before
    /// handing a reply to this bounded shared cache. Large replies are not cached.
    pub fn finish(
        mut self,
        state: wire::OperationState,
        checked_secret_free_reply: &[u8],
    ) -> HostResult<()> {
        if !matches!(
            state,
            wire::OperationState::Completed
                | wire::OperationState::Failed
                | wire::OperationState::Cancelled
        ) {
            return Err(HostError::Invalid(
                "External job requires a terminal state.",
            ));
        }
        if self.accepted() && state == wire::OperationState::Cancelled {
            return Err(HostError::Invalid(
                "Accepted work cannot be reported as rolled back.",
            ));
        }
        let mut registry = self
            .control
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        let cached = registry
            .entries
            .values()
            .filter_map(|entry| entry.reply.as_ref())
            .map(Vec::len)
            .sum::<usize>();
        let entry = registry
            .entries
            .get_mut(&(self.access.id, self.operation.as_str().into()))
            .ok_or(HostError::Unavailable)?;
        if entry.digest.as_ref() != Some(&self.digest) {
            return Err(HostError::Invalid("Operation intent changed."));
        }
        if !self.accepted()
            && (!self.access.alive.load(Ordering::Acquire)
                || matches!(
                    entry.state,
                    wire::OperationState::Cancelled | wire::OperationState::CancelRequested
                ))
            && state != wire::OperationState::Cancelled
        {
            return Err(HostError::Cancelled);
        }
        if !matches!(
            entry.state,
            wire::OperationState::Running
                | wire::OperationState::Committing
                | wire::OperationState::CancelRequested
                | wire::OperationState::Cancelled
        ) {
            return Err(HostError::Invalid("Operation already finished."));
        }
        entry.state = state;
        if checked_secret_free_reply.len() <= 32768
            && cached + checked_secret_free_reply.len() <= MAX_CACHED_BYTES
        {
            entry.reply = Some(checked_secret_free_reply.to_vec());
        }
        self.finished = true;
        Ok(())
    }
}
impl Drop for ScopedJob {
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        let Ok(mut registry) = self.control.inner.registry.lock() else {
            return;
        };
        let Some(entry) = registry
            .entries
            .get_mut(&(self.access.id, self.operation.as_str().into()))
        else {
            return;
        };
        if entry.digest.as_ref() != Some(&self.digest)
            || !matches!(
                entry.state,
                wire::OperationState::Running
                    | wire::OperationState::Committing
                    | wire::OperationState::CancelRequested
            )
        {
            return;
        }
        entry.state = if !self.accepted()
            && (!self.access.alive.load(Ordering::Acquire)
                || entry.state == wire::OperationState::CancelRequested)
        {
            wire::OperationState::Cancelled
        } else {
            wire::OperationState::Failed
        };
        entry.reply = None;
    }
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
        let mut registry = self
            .inner
            .registry
            .lock()
            .map_err(|_| HostError::Unavailable)?;
        self.check(access)?;
        // Revocation and accept/commit now share one lock: either it prevents
        // acceptance, or an already-accepted native outcome remains authoritative.
        access.alive.store(false, Ordering::Release);
        for ((window, _), entry) in registry.entries.iter_mut() {
            if *window == access.id
                && matches!(
                    entry.state,
                    wire::OperationState::Queued | wire::OperationState::Running
                )
            {
                entry.state = wire::OperationState::CancelRequested;
            }
        }
        registry.entries.retain(|(window, _), entry| {
            *window != access.id
                || (entry.seen_job
                    && matches!(
                        entry.state,
                        wire::OperationState::Committing | wire::OperationState::CancelRequested
                    ))
        });
        // Retain in-flight entries until their owner terminalizes them. Revoked
        // capabilities cannot query them, and ordinary bounded eviction applies.
        Ok(())
    }
    pub fn begin_job(
        &self,
        access: &WindowAccess,
        operation: &wire::OperationId,
        validated_frame_digest: Digest,
    ) -> HostResult<JobStart> {
        match self.begin(access, operation, validated_frame_digest.clone())? {
            Begin::Cached(reply) => Ok(JobStart::Cached(reply)),
            Begin::Run => Ok(JobStart::Run(ScopedJob {
                control: self.clone(),
                access: access.clone(),
                operation: operation.clone(),
                digest: validated_frame_digest,
                accepted: AtomicBool::new(false),
                finished: false,
            })),
        }
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
        self.check(access)?;
        let key = (access.id, id.as_str().into());
        if let Some(entry) = registry.entries.get_mut(&key) {
            if entry
                .digest
                .as_ref()
                .is_some_and(|previous| previous != &digest)
            {
                return Err(HostError::Invalid(
                    "Operation ID was already used for another request.",
                ));
            }
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
        self.check(access)?;
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
        self.check(access)?;
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

#[cfg(test)]
mod external_job_tests {
    use super::*;
    fn op(n: u64) -> wire::OperationId {
        wire::OperationId::new(format!("op-{n:032x}")).unwrap()
    }
    fn start(control: &HostControl, window: &WindowAccess, n: u64) -> ScopedJob {
        match control
            .begin_job(
                window,
                &op(n),
                Digest::of(format!("validated-dsq-frame-{n}").as_bytes()),
            )
            .unwrap()
        {
            JobStart::Run(job) => job,
            JobStart::Cached(_) => panic!("unexpected replay"),
        }
    }
    fn cancel(
        control: &HostControl,
        window: &WindowAccess,
        target: wire::OperationId,
    ) -> wire::OperationState {
        let request = wire::RequestEnvelope {
            request_id: 44,
            operation_id: op(99999),
            command: wire::Command::Cancel {
                operation_id: target,
            },
        };
        let reply = control
            .dispatch_raw(window, &wire::encode_request(&request).unwrap())
            .unwrap()
            .unwrap();
        let wire::Outcome::Success(wire::Success::Cancelled { status, .. }) =
            wire::decode_response_for(&reply, &request).unwrap().outcome
        else {
            panic!("wrong control reply")
        };
        status
    }
    fn status(control: &HostControl, window: &WindowAccess, n: u64) -> wire::OperationState {
        control.inner.registry.lock().unwrap().entries[&(window.id, op(n).as_str().into())].state
    }
    #[test]
    fn dtq_control_cancels_external_job_before_acceptance_and_revocation_blocks_new_authority() {
        let control = HostControl::new();
        let window = control.authorize().unwrap();
        let job = start(&control, &window, 1);
        assert_eq!(
            cancel(&control, &window, op(1)),
            wire::OperationState::CancelRequested
        );
        assert!(job.cancelled_before_acceptance());
        assert!(matches!(job.accept(), Err(HostError::Cancelled)));
        job.finish(
            wire::OperationState::Cancelled,
            b"checked secret-free cancelled reply",
        )
        .unwrap();
        assert_eq!(
            status(&control, &window, 1),
            wire::OperationState::Cancelled
        );
        let revoked = control.authorize().unwrap();
        let job = start(&control, &revoked, 2);
        control.revoke_window(&revoked).unwrap();
        assert!(job.cancelled_before_acceptance());
        assert!(matches!(job.accept(), Err(HostError::Window)));
        drop(job);
        assert_eq!(
            status(&control, &revoked, 2),
            wire::OperationState::Cancelled
        );
        assert!(matches!(
            control.begin_job(&revoked, &op(3), Digest::of(b"validated")),
            Err(HostError::Window)
        ));
    }
    #[test]
    fn accepted_external_outcome_survives_late_cancel_and_window_revocation() {
        let control = HostControl::new();
        let window = control.authorize().unwrap();
        let job = start(&control, &window, 10);
        job.accept().unwrap();
        job.accept().unwrap();
        assert!(job.accepted());
        assert_eq!(
            cancel(&control, &window, op(10)),
            wire::OperationState::Committing
        );
        control.revoke_window(&window).unwrap();
        assert!(!job.cancelled_before_acceptance());
        job.finish(
            wire::OperationState::Completed,
            b"checked secret-free completed receipt",
        )
        .unwrap();
        let registry = control.inner.registry.lock().unwrap();
        let entry = &registry.entries[&(window.id, op(10).as_str().into())];
        assert_eq!(entry.state, wire::OperationState::Completed);
        assert_eq!(
            entry.reply.as_deref(),
            Some(b"checked secret-free completed receipt".as_slice())
        );
    }
    #[test]
    fn changed_intent_and_foreign_windows_cannot_overwrite_original_checked_reply() {
        let control = HostControl::new();
        let window = control.authorize().unwrap();
        let other = HostControl::new();
        let foreign = other.authorize().unwrap();
        assert!(matches!(
            control.begin_job(&foreign, &op(20), Digest::of(b"validated")),
            Err(HostError::Window)
        ));
        let digest = Digest::of(b"namespace-bound validated DSQ frame");
        let JobStart::Run(job) = control.begin_job(&window, &op(20), digest.clone()).unwrap()
        else {
            panic!("new job missing")
        };
        job.accept().unwrap();
        job.finish(wire::OperationState::Completed, b"checked original receipt")
            .unwrap();
        assert!(matches!(
            control.begin_job(&window, &op(20), Digest::of(b"changed frame")),
            Err(HostError::Invalid(_))
        ));
        let JobStart::Cached(reply) = control.begin_job(&window, &op(20), digest).unwrap() else {
            panic!("original receipt was overwritten")
        };
        assert_eq!(reply, b"checked original receipt");
    }
    #[test]
    fn abandoned_guard_terminalizes_without_claiming_an_accepted_effect_rolled_back() {
        let control = HostControl::new();
        let window = control.authorize().unwrap();
        drop(start(&control, &window, 30));
        assert_eq!(status(&control, &window, 30), wire::OperationState::Failed);
        let job = start(&control, &window, 31);
        job.accept().unwrap();
        drop(job);
        assert_eq!(status(&control, &window, 31), wire::OperationState::Failed);
        let job = start(&control, &window, 32);
        job.accept().unwrap();
        assert!(matches!(
            job.finish(wire::OperationState::Cancelled, b"incorrect rollback claim"),
            Err(HostError::Invalid(_))
        ));
        assert_eq!(status(&control, &window, 32), wire::OperationState::Failed);
        let job = start(&control, &window, 33);
        cancel(&control, &window, op(33));
        drop(job);
        assert_eq!(
            status(&control, &window, 33),
            wire::OperationState::Cancelled
        );
    }
    #[test]
    fn external_jobs_use_the_existing_bounded_reply_cache_and_drop_cleans_revoked_tombstones() {
        let control = HostControl::new();
        let window = control.authorize().unwrap();
        let job = start(&control, &window, 40);
        job.finish(wire::OperationState::Completed, &vec![0u8; 32769])
            .unwrap();
        assert!(
            control.inner.registry.lock().unwrap().entries[&(window.id, op(40).as_str().into())]
                .reply
                .is_none()
        );
        cancel(&control, &window, op(41));
        assert!(
            control
                .inner
                .registry
                .lock()
                .unwrap()
                .entries
                .contains_key(&(window.id, op(41).as_str().into()))
        );
        control.revoke_window(&window).unwrap();
        assert!(
            !control
                .inner
                .registry
                .lock()
                .unwrap()
                .entries
                .contains_key(&(window.id, op(41).as_str().into()))
        );
    }
}
