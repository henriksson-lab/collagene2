//! Opaque native review lifetime and ordered vault/SQLite remote retirement.
use crate::*;
use repository::remote_removal::RemovalReview;
use std::collections::{BTreeMap, VecDeque};
pub(crate) struct Reviews {
    seed: [u8; 16],
    next: u64,
    items: BTreeMap<wire::ReviewId, RemovalReview>,
    order: VecDeque<wire::ReviewId>,
}
impl Reviews {
    pub fn new() -> std::result::Result<Self, HostError> {
        let mut seed = [0; 16];
        getrandom::fill(&mut seed).map_err(|_| HostError::Unavailable)?;
        Ok(Self {
            seed,
            next: 0,
            items: BTreeMap::new(),
            order: VecDeque::new(),
        })
    }
    pub fn clear(&mut self) {
        self.items.clear();
        self.order.clear();
    }
    fn remove(&mut self, id: &wire::ReviewId) {
        self.items.remove(id);
        self.order.retain(|held| held != id);
    }
    fn clear_remote(&mut self, remote: &str) {
        let ids: Vec<_> = self
            .items
            .iter()
            .filter(|(_, r)| r.description().remote.id == remote)
            .map(|(id, _)| id.clone())
            .collect();
        for id in ids {
            self.remove(&id);
        }
    }
    fn insert(&mut self, review: RemovalReview) -> Result<wire::ReviewId> {
        if self.items.len() >= 4 {
            if let Some(old) = self.order.pop_front() {
                self.items.remove(&old);
            }
        }
        self.next = self.next.checked_add(1).ok_or_else(Fault::unavailable)?;
        let hash = library::Digest::of(&[self.seed.as_slice(), &self.next.to_le_bytes()].concat());
        let id = wire::ReviewId::new(format!("review-{}", &hash.as_str()[..32]))
            .map_err(|_| Fault::unavailable())?;
        self.items.insert(id.clone(), review);
        self.order.push_back(id.clone());
        Ok(id)
    }
}
fn stale() -> Fault {
    Fault::error(
        wire::ErrorCode::ReviewChanged,
        "Remote removal review expired or changed. Review the exact current database configuration again; no key was forgotten.",
    )
}
impl<A: SyncAdapter> SyncHost<A> {
    pub(crate) fn execute_removal(
        &mut self,
        command: &wire::Command,
        gate: &mut impl FnMut() -> bool,
    ) -> Result<Option<wire::Success>> {
        use wire::{Command as C, Success as S};
        Ok(Some(match command {
            C::PrepareRemoteRemoval { remote } => {
                self.remote(remote)?;
                let review = self
                    .repository
                    .prepare_remote_removal(remote.as_str())
                    .map_err(repository_error)?;
                let description = convert(review.description())?;
                let digest = review.digest();
                accepted(gate)?;
                let id = self.removals.insert(review)?;
                S::RemoteRemovalPrepared {
                    review: wire::RemoteRemovalReview {
                        id,
                        description,
                        digest,
                    },
                }
            }
            C::ConfirmRemoteRemoval { remote, review } => {
                let held = self.removals.items.get(review).cloned().ok_or_else(stale)?;
                if held.description().remote.id != remote.as_str() {
                    return Err(stale());
                }
                let fresh = self
                    .repository
                    .prepare_remote_removal(remote.as_str())
                    .map_err(|_| stale())?;
                if fresh.digest() != held.digest() {
                    return Err(stale());
                }
                let approved = held.approve();
                let scope = held.description().remote.clone();
                accepted(gate)?;
                // Even a refused vault delete may already clear the active
                // session credential. Existing publication acknowledgements for
                // this exact remote must expire before that credential change.
                self.reviews.clear_remote(remote.as_str());
                // Credentials::delete clears matching session state then addresses the exact
                // old native vault scope. A missing entry is already the desired state.
                match self.adapter.delete_credential(&scope) {
                    Ok(()) | Err(credentials::Error::Missing) => {}
                    Err(_) => {
                        return Err(Fault::error(
                            wire::ErrorCode::CredentialForgetFailed,
                            "The native key could not be forgotten. Database configuration, bindings and local designs were retained; its session key may already be cleared. Unlock the vault and review before retrying.",
                        ));
                    }
                }
                let result = self.repository.remove_remote_reviewed(&approved);
                let receipt = match result {
                    Ok(r) => r,
                    Err(error) => {
                        // A completed metadata commit remains authoritative even if a subsequent
                        // reload failed. Recheck its exact durable receipt before reporting doubt.
                        if let Ok(Some(r)) = self.repository.remote_removal_receipt(remote.as_str())
                        {
                            if r.review_digest == held.digest() {
                                r
                            } else {
                                return Err(Fault::error(
                                    wire::ErrorCode::RemovalChangedAfterForget,
                                    "The key was forgotten, but local unlink could not be confirmed. Check its durable outcome before another action; local designs were not deleted.",
                                ));
                            }
                        } else {
                            let known_uncommitted = matches!(
                                error,
                                repository::Error::Stale
                                    | repository::Error::Pending
                                    | repository::Error::Transaction(
                                        store::mutations::MutationError::CheckpointConflict { .. }
                                    )
                            );
                            return Err(Fault::error(
                                wire::ErrorCode::RemovalChangedAfterForget,
                                if known_uncommitted {
                                    "The key was forgotten, but concurrent local changes prevented unlinking. Database configuration, bindings and local designs remain; prepare a fresh review before retrying."
                                } else {
                                    "The key was forgotten, but local unlink could not be confirmed. Check its durable outcome before another action; local designs were not deleted."
                                },
                            ));
                        }
                    }
                };
                self.removals.clear_remote(remote.as_str());
                self.reviews.clear_remote(remote.as_str());
                S::RemoteRemoved(convert(&receipt)?)
            }
            C::DiscardRemoteRemovalReview { review } => {
                self.removals.remove(review);
                S::RemoteRemovalReviewDiscarded {
                    review: review.clone(),
                }
            }
            C::CheckRemoteRemovalOutcome { remote } => {
                let receipt = self
                    .repository
                    .remote_removal_receipt(remote.as_str())
                    .map_err(repository_error)?
                    .as_ref()
                    .map(convert)
                    .transpose()?;
                S::RemoteRemovalOutcome {
                    remote: remote.clone(),
                    receipt,
                }
            }
            _ => return Ok(None),
        }))
    }
}
