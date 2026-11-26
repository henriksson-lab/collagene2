//! Pure native remote-removal review/receipt scope checks.
use collagene_desktop_sync_protocol as wire;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PinnedRemoval {
    pub remote: wire::Remote,
    pub review: wire::RemoteRemovalReview,
}
fn configuration(remote: &wire::Remote, scope: &wire::PublicationRemote) -> bool {
    remote.id == scope.id
        && remote.name == scope.name
        && remote.endpoint == scope.endpoint
        && remote.expected_lab == scope.lab
}
pub fn pin(
    remote: &wire::Remote,
    review: wire::RemoteRemovalReview,
) -> Result<PinnedRemoval, String> {
    remote
        .validate()
        .map_err(|_| "Invalid configured database.")?;
    review
        .description
        .validate()
        .map_err(|_| "Invalid removal review.")?;
    if !configuration(remote, &review.description.remote) {
        return Err("Removal review belongs to another database configuration.".into());
    }
    Ok(PinnedRemoval {
        remote: remote.clone(),
        review,
    })
}
pub fn is_current(pinned: &PinnedRemoval, selected: Option<&wire::Remote>) -> bool {
    selected == Some(&pinned.remote)
}
pub fn validate_receipt(
    remote: &wire::Remote,
    review: Option<&PinnedRemoval>,
    receipt: &wire::RemovalReceipt,
) -> Result<(), String> {
    receipt.validate().map_err(|_| "Invalid removal receipt.")?;
    if !configuration(remote, &receipt.remote) {
        return Err("Removal outcome belongs to another database configuration.".into());
    }
    if let Some(p) = review {
        if p.remote != *remote
            || p.review.digest != receipt.review_digest
            || p.review.description.local_lab != receipt.local_lab
            || p.review.description.document_bindings != receipt.removed_document_bindings
            || p.review.description.folder_bindings != receipt.removed_folder_bindings
            || p.review.description.completed_tree_jobs != receipt.removed_completed_tree_jobs
        {
            return Err("Removal outcome does not match the reviewed scope.".into());
        }
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    use wire::library;
    fn remote() -> wire::Remote {
        wire::Remote {
            id: wire::RemoteId::new(format!("rem-{:032x}", 1)).unwrap(),
            name: library::Text::new("Lab").unwrap(),
            endpoint: library::Text::new("https://lab.example").unwrap(),
            expected_lab: library::LabId::new(format!("lab-{:032x}", 1)).unwrap(),
        }
    }
    fn review(r: &wire::Remote) -> wire::RemoteRemovalReview {
        wire::RemoteRemovalReview {
            id: wire::ReviewId::new(format!("review-{:032x}", 1)).unwrap(),
            digest: library::Digest::of(b"native reviewed scope"),
            description: wire::RemovalDescription {
                local_lab: library::LabId::new(format!("lab-{:032x}", 2)).unwrap(),
                remote: wire::PublicationRemote {
                    id: r.id.clone(),
                    name: r.name.clone(),
                    endpoint: r.endpoint.clone(),
                    lab: r.expected_lab.clone(),
                },
                expected_checkpoint: Some(library::Digest::of(b"checkpoint")),
                document_bindings: 2,
                folder_bindings: 3,
                completed_tree_jobs: 1,
            },
        }
    }
    #[test]
    fn selection_and_full_configuration_are_bound_before_confirmation() {
        let r = remote();
        let p = pin(&r, review(&r)).unwrap();
        assert!(is_current(&p, Some(&r)));
        let mut changed = r.clone();
        changed.endpoint = library::Text::new("https://other.example").unwrap();
        assert!(!is_current(&p, Some(&changed)));
        assert!(pin(&changed, p.review.clone()).is_err());
        assert!(!is_current(&p, None));
    }
    #[test]
    fn outcome_must_match_review_counts_local_context_digest_and_complete_remote() {
        let r = remote();
        let p = pin(&r, review(&r)).unwrap();
        let receipt = wire::RemovalReceipt {
            schema_version: 1,
            local_lab: p.review.description.local_lab.clone(),
            remote: p.review.description.remote.clone(),
            review_digest: p.review.digest.clone(),
            removed_document_bindings: 2,
            removed_folder_bindings: 3,
            removed_completed_tree_jobs: 1,
        };
        assert!(validate_receipt(&r, Some(&p), &receipt).is_ok());
        for mode in 0..4 {
            let mut wrong = receipt.clone();
            match mode {
                0 => wrong.review_digest = library::Digest::of(b"wrong"),
                1 => wrong.local_lab = r.expected_lab.clone(),
                2 => wrong.removed_folder_bindings += 1,
                _ => wrong.remote.endpoint = library::Text::new("https://wrong.example").unwrap(),
            };
            assert!(validate_receipt(&r, Some(&p), &wrong).is_err());
        }
        assert!(validate_receipt(&r, None, &receipt).is_ok());
    }
}
