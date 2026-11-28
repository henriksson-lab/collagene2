//! A cancelled pre-transaction operation cannot unlock a subsequent commit.
use std::cell::Cell;

#[derive(Default)]
pub struct CommitLease(Cell<Option<u64>>);
impl CommitLease {
    pub fn acquire(&self, owner: u64) -> Result<(), String> {
        if self.0.get().is_some() {
            return Err("A workspace commit is already pending.".into());
        }
        self.0.set(Some(owner));
        Ok(())
    }
    pub fn release(&self, owner: u64) -> bool {
        if self.0.get() != Some(owner) {
            return false;
        }
        self.0.set(None);
        true
    }
    pub fn cancel_unstarted(&self) -> bool {
        self.0.take().is_some()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cancelled_encoder_completion_cannot_unlock_a_new_workspace_transaction() {
        let lease = CommitLease::default();
        lease.acquire(11).unwrap();
        assert!(lease.cancel_unstarted());
        lease.acquire(12).unwrap();
        // The old cancelled future finishes after a new operation has begun.
        assert!(!lease.release(11));
        assert!(lease.acquire(13).is_err());
        assert!(lease.release(12));
        lease.acquire(13).unwrap();
        assert!(lease.release(13));
    }
}
