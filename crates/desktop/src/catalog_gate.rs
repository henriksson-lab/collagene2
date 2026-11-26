//! Process-local catalog admission. Working-editor files do not take this gate.
use std::sync::{
    RwLock, RwLockReadGuard, RwLockWriteGuard,
    atomic::{AtomicBool, Ordering},
};
#[derive(Default)]
pub(crate) struct CatalogGate {
    pending: AtomicBool,
    gate: RwLock<()>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Error {
    Busy,
    PendingRestart,
    Unavailable,
}
impl CatalogGate {
    pub fn new(paused: bool) -> Self {
        Self {
            pending: AtomicBool::new(paused),
            gate: RwLock::new(()),
        }
    }
    pub fn read(&self) -> Result<RwLockReadGuard<'_, ()>, Error> {
        let lease = self.gate.try_read().map_err(|e| match e {
            std::sync::TryLockError::WouldBlock => Error::Busy,
            _ => Error::Unavailable,
        })?;
        if self.pending.load(Ordering::Acquire) {
            return Err(Error::PendingRestart);
        }
        Ok(lease)
    }
    pub fn status(&self) -> Result<RwLockReadGuard<'_, ()>, Error> {
        self.gate.try_read().map_err(|e| match e {
            std::sync::TryLockError::WouldBlock => Error::Busy,
            _ => Error::Unavailable,
        })
    }
    pub fn adoption(&self) -> Result<Adoption<'_>, Error> {
        let lease = self.gate.try_write().map_err(|e| match e {
            std::sync::TryLockError::WouldBlock => Error::Busy,
            _ => Error::Unavailable,
        })?;
        Ok(Adoption {
            owner: self,
            _lease: lease,
        })
    }
    pub fn pending(&self) -> bool {
        self.pending.load(Ordering::Acquire)
    }
}
pub(crate) struct Adoption<'a> {
    owner: &'a CatalogGate,
    _lease: RwLockWriteGuard<'a, ()>,
}
impl Adoption<'_> {
    /// Call only after native durable acceptance or discovering a retained
    /// pending activation. Error paths do not claim adoption succeeded.
    pub fn mark_pending(&self) {
        self.owner.pending.store(true, Ordering::Release);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn adoption_refuses_running_catalog_jobs_then_rejects_old_sessions_after_acceptance() {
        let gate = CatalogGate::default();
        let old = gate.read().unwrap();
        assert!(matches!(gate.adoption(), Err(Error::Busy)));
        drop(old);
        {
            let adoption = gate.adoption().unwrap();
            assert!(matches!(gate.read(), Err(Error::Busy)));
            adoption.mark_pending();
        }
        assert!(gate.pending());
        assert!(matches!(gate.read(), Err(Error::PendingRestart)));
        // A known adoption/status retry can still acquire the exclusive gate.
        assert!(gate.adoption().is_ok());
    }
    #[test]
    fn editor_only_startup_blocks_catalog_but_allows_status_and_cancel_admission() {
        let gate = CatalogGate::new(true);
        assert!(matches!(gate.read(), Err(Error::PendingRestart)));
        assert!(gate.status().is_ok());
        let cancellation = gate.adoption().unwrap();
        cancellation.mark_pending();
        drop(cancellation);
        assert!(matches!(gate.read(), Err(Error::PendingRestart)));
    }
    #[test]
    fn failed_or_cancelled_adoption_does_not_activate_a_recovery() {
        let gate = CatalogGate::default();
        {
            let _review = gate.adoption().unwrap();
        }
        assert!(!gate.pending());
        assert!(gate.read().is_ok());
    }
}
