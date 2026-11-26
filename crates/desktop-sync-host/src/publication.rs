//! Native-only review lifetime and source-bound queue admission.
use crate::*;
use repository::{folder_publication as metadata, publication_review::PushReview};
use std::collections::{BTreeMap, VecDeque};
pub(crate) struct Reviews {
    seed: [u8; 16],
    next: u64,
    items: BTreeMap<wire::ReviewId, PushReview>,
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
    pub(crate) fn clear_remote(&mut self, remote: &str) {
        let ids: Vec<_> = self
            .items
            .iter()
            .filter(|(_, review)| review.remote().id == remote)
            .map(|(id, _)| id.clone())
            .collect();
        for id in ids {
            self.remove(&id);
        }
    }
    fn remove(&mut self, id: &wire::ReviewId) {
        self.items.remove(id);
        self.order.retain(|held| held != id);
    }
    fn insert(&mut self, review: PushReview) -> Result<wire::ReviewId> {
        if self.items.len() >= 4 {
            if let Some(old) = self.order.pop_front() {
                self.items.remove(&old);
            }
        }
        self.next = self.next.checked_add(1).ok_or_else(Fault::unavailable)?;
        let digest =
            library::Digest::of(&[self.seed.as_slice(), &self.next.to_le_bytes()].concat());
        let id = wire::ReviewId::new(format!("review-{}", &digest.as_str()[..32]))
            .map_err(|_| Fault::unavailable())?;
        self.order.push_back(id.clone());
        self.items.insert(id.clone(), review);
        Ok(id)
    }
}
fn error(e: repository::Error) -> Fault {
    if matches!(e, repository::Error::Stale) {
        Fault::error(
            wire::ErrorCode::ReviewChanged,
            "Saved source or destination changed; prepare and review it again.",
        )
    } else {
        repository_error(e)
    }
}
fn choices(
    rows: &library::Rows<wire::PushChoice, 128>,
) -> sharing::sbol_subset::PublicationChoices {
    let mut out = sharing::sbol_subset::PublicationChoices::new();
    for row in rows.iter() {
        out.entry(row.original_hash.as_str().into())
            .or_default()
            .insert(
                row.component_iri.as_str().into(),
                row.sequences.iter().map(|s| s.as_str().into()).collect(),
            );
    }
    out
}
impl<A: SyncAdapter> SyncHost<A> {
    pub(crate) fn execute_publication(
        &mut self,
        c: &wire::Command,
        gate: &mut impl FnMut() -> bool,
    ) -> Result<Option<wire::Success>> {
        use wire::{Command as C, Success as S};
        Ok(Some(match c {
            C::PrepareMetadataPublication { target, remote } => {
                self.remote(remote)?;
                let target: metadata::PublicationTarget = convert(target)?;
                let review = self
                    .repository
                    .prepare_metadata_publication(&target, remote.as_str())
                    .map_err(error)?;
                S::MetadataPrepared(wire::MetadataReview {
                    description: convert(review.description())?,
                    digest: review.digest(),
                })
            }
            C::QueueMetadataPublication { review } => {
                let description: metadata::PublicationDescription = convert(&review.description)?;
                self.remote(&review.description.remote.id)?;
                accepted(gate)?;
                let id = self
                    .repository
                    .queue_metadata_publication_checked(&description, &review.digest)
                    .map_err(error)?;
                S::Queued(self.pending(&id)?)
            }
            C::InspectPushChoices {
                local_document,
                expected_local_revision,
                remote,
            } => {
                self.remote(remote)?;
                let id = document(local_document)?;
                let revision = revision(expected_local_revision)?;
                self.saved(&id, &revision)?;
                let rows = self
                    .repository
                    .push_review_choices(&id, &revision, remote.as_str())
                    .map_err(error)?;
                S::PushChoices(wire::PushChoices {
                    local_document: local_document.clone(),
                    local_revision: expected_local_revision.clone(),
                    remote: remote.clone(),
                    requirements: library::Rows::new(
                        rows.iter()
                            .map(|r| {
                                Ok(wire::PushChoiceRequirement {
                                    original_hash: library::Digest::new(&r.original_hash)
                                        .map_err(|_| Fault::too_large())?,
                                    component_iri: library::Text::new(&r.component_iri)
                                        .map_err(|_| Fault::too_large())?,
                                    available_sequences: convert(&r.available_sequences)?,
                                    required_sequences: convert(&r.required_sequences)?,
                                })
                            })
                            .collect::<Result<Vec<_>>>()?,
                    )
                    .map_err(|_| Fault::too_large())?,
                })
            }
            C::PrepareReviewedPush {
                local_document,
                expected_local_revision,
                remote,
                choices: rows,
            } => {
                self.remote(remote)?;
                let id = document(local_document)?;
                let rev = revision(expected_local_revision)?;
                self.saved(&id, &rev)?;
                let review = self
                    .repository
                    .prepare_push_review(&id, &rev, remote.as_str(), &choices(rows))
                    .map_err(error)?;
                let entries = review
                    .review()
                    .entries()
                    .iter()
                    .map(|r| {
                        Ok(wire::PushReviewEntry {
                            original_hash: library::Digest::new(&r.original_hash)
                                .map_err(|_| Fault::too_large())?,
                            subset_hash: library::Digest::new(&r.subset_hash)
                                .map_err(|_| Fault::too_large())?,
                            retained_objects: r
                                .retained_objects
                                .try_into()
                                .map_err(|_| Fault::too_large())?,
                            omitted_objects: r
                                .omissions
                                .unrelated_objects
                                .try_into()
                                .map_err(|_| Fault::too_large())?,
                            omitted_sequence_edges: r
                                .omissions
                                .alternative_sequence_edges
                                .try_into()
                                .map_err(|_| Fault::too_large())?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                let checkpoint = review.checkpoint().as_ref().map(convert).transpose()?;
                let entries = library::Rows::new(entries).map_err(|_| Fault::too_large())?;
                let configuration = convert(review.remote())?;
                accepted(gate)?;
                let review_id = self.reviews.insert(review)?;
                S::PushPrepared(wire::PushPublicationReview {
                    id: review_id,
                    local_document: local_document.clone(),
                    local_revision: expected_local_revision.clone(),
                    remote: remote.clone(),
                    configuration,
                    checkpoint,
                    entries,
                })
            }
            C::QueueReviewedPush {
                local_document,
                expected_local_revision,
                remote,
                review,
            } => {
                let held = self.reviews.items.get(review).ok_or_else(|| {
                    Fault::error(
                        wire::ErrorCode::ReviewChanged,
                        "Publication review expired; prepare it again.",
                    )
                })?;
                if held.local_document().as_str() != local_document.as_str()
                    || held.local_revision().as_str() != expected_local_revision.as_str()
                    || held.remote().id != remote.as_str()
                {
                    return Err(Fault::invalid("Publication review source does not match."));
                }
                let approved = held.approve();
                accepted(gate)?;
                let id = self
                    .repository
                    .queue_reviewed_push(&approved)
                    .map_err(error)?;
                self.reviews.remove(review);
                S::Queued(self.pending(&id)?)
            }
            C::DiscardPushReview { review } => {
                self.reviews.remove(review);
                S::PushReviewDiscarded {
                    review: review.clone(),
                }
            }
            C::ListDeletedFolderPublications {
                remote,
                page,
                expected_checkpoint,
            } => {
                self.remote(remote)?;
                let all = self
                    .repository
                    .deleted_folder_publications(remote.as_str())
                    .map_err(error)?;
                let checkpoint = self.repository.checkpoint_version();
                if page.offset > 0
                    && checkpoint != expected_checkpoint.as_ref().map(convert).transpose()?
                {
                    return Err(Fault::error(
                        wire::ErrorCode::ReviewChanged,
                        "Deleted folder list changed; refresh it.",
                    ));
                }
                let start = usize::from(page.offset);
                let end = start.saturating_add(usize::from(page.limit)).min(all.len());
                if start > all.len() {
                    return Err(Fault::invalid("Deleted folder page is outside the list."));
                }
                S::DeletedFolderPublications(wire::DeletedFolderPage {
                    remote: remote.clone(),
                    items: convert(&all[start..end].to_vec())?,
                    next_offset: (end < all.len()).then(|| end.try_into().ok()).flatten(),
                    checkpoint: checkpoint.as_ref().map(convert).transpose()?,
                })
            }
            _ => return Ok(None),
        }))
    }
}
