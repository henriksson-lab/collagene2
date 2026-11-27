//! Pure backup-stage and receipt pinning. No file path or archive enters this UI.
use collagene_desktop_protocol::{self as wire, library};
#[derive(Clone)]
pub(crate) struct VerificationPin {
    context: library::LabContext,
    id: wire::OperationId,
    expected_hash: Option<library::Digest>,
    basename: library::Text<128>,
}
#[derive(Clone)]
pub(crate) struct RestoreReview {
    context: library::LabContext,
    receipt: wire::CatalogBackupReceipt,
}
fn correct_lab(
    receipt: &wire::CatalogBackupReceipt,
    context: &library::LabContext,
) -> Result<(), String> {
    receipt
        .validate()
        .map_err(|_| "Native backup receipt is malformed; refresh known backups.".to_string())?;
    if receipt.lab != context.lab_id || receipt.lab_count != 1 {
        return Err("Backup receipt belongs to a different local catalog scope; no restore review was enabled.".into());
    }
    Ok(())
}
impl VerificationPin {
    pub fn new(
        entry: &wire::CatalogBackupEntry,
        context: &library::LabContext,
    ) -> Result<Self, String> {
        if entry.basename.as_str() != format!("{}.cglbak", entry.id.as_str()) {
            return Err("Known backup name does not match its managed identity.".into());
        }
        if let Some(receipt) = &entry.receipt {
            if receipt.id != entry.id || receipt.basename != entry.basename {
                return Err(
                    "Known backup entry and receipt disagree. Refresh the backup list.".into(),
                );
            }
            correct_lab(receipt, context)?;
        }
        Ok(Self {
            context: context.clone(),
            id: entry.id.clone(),
            basename: entry.basename.clone(),
            expected_hash: entry.receipt.as_ref().map(|r| r.archive_hash.clone()),
        })
    }
    pub fn command(&self) -> wire::CatalogBackupCommand {
        wire::CatalogBackupCommand::Verify {
            id: self.id.clone(),
            expected_hash: self.expected_hash.clone(),
        }
    }
    pub fn accept(
        &self,
        receipt: wire::CatalogBackupReceipt,
        current: &library::LabContext,
    ) -> Result<RestoreReview, String> {
        if current != &self.context
            || receipt.id != self.id
            || receipt.basename != self.basename
            || self
                .expected_hash
                .as_ref()
                .is_some_and(|h| h != &receipt.archive_hash)
        {
            return Err("Verified backup or local context changed. Refresh and verify the selected backup again.".into());
        }
        correct_lab(&receipt, current)?;
        Ok(RestoreReview {
            context: current.clone(),
            receipt,
        })
    }
}
impl RestoreReview {
    pub fn receipt(&self) -> &wire::CatalogBackupReceipt {
        &self.receipt
    }
    pub fn command(
        &self,
        current: &library::LabContext,
    ) -> Result<wire::CatalogBackupCommand, String> {
        if current != &self.context {
            return Err("The local catalog context changed; verify the backup again.".into());
        }
        Ok(wire::CatalogBackupCommand::RestoreCopy {
            id: self.receipt.id.clone(),
            expected_hash: self.receipt.archive_hash.clone(),
        })
    }
    pub fn accept(
        &self,
        receipt: &wire::CatalogRecoveryReceipt,
        current: &library::LabContext,
    ) -> Result<(), String> {
        if current != &self.context
            || receipt.backup != self.receipt
            || receipt.directory_basename.as_str() != receipt.operation_id.as_str()
            || receipt.database_basename.as_str() != "library.sqlite"
        {
            return Err(
                "Restored-copy receipt does not match this verified backup/context.".into(),
            );
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (library::LabContext, wire::CatalogBackupEntry) {
        let context = library::LabContext {
            lab_id: library::LabId::new(format!("lab-{:032x}", 1)).unwrap(),
            name: library::Text::new("Local catalog").unwrap(),
        };
        let receipt = wire::CatalogBackupReceipt {
            id: wire::OperationId::new(format!("op-{:032x}", 2)).unwrap(),
            basename: library::Text::new(format!("op-{:032x}.cglbak", 2)).unwrap(),
            archive_hash: library::Digest::new("11".repeat(32)).unwrap(),
            archive_bytes: 4096,
            database_hash: library::Digest::new("22".repeat(32)).unwrap(),
            database_bytes: 2048,
            schema_hash: library::Digest::new("33".repeat(32)).unwrap(),
            schema_version: library::DATABASE_SCHEMA_VERSION,
            lab: context.lab_id.clone(),
            lab_count: 1,
            documents: 3,
            folders: 2,
            revisions: 5,
            mutation_receipts: 4,
            created_unix_seconds: 1795680000,
            request_checkpoint: None,
        };
        (
            context,
            wire::CatalogBackupEntry {
                id: receipt.id.clone(),
                basename: receipt.basename.clone(),
                receipt: Some(receipt),
            },
        )
    }
    #[test]
    fn known_completed_entry_pins_verified_hash_before_enabling_restore() {
        let (context, entry) = fixture();
        let pin = VerificationPin::new(&entry, &context).unwrap();
        assert!(matches!(
            pin.command(),
            wire::CatalogBackupCommand::Verify {
                expected_hash: Some(_),
                ..
            }
        ));
        let mut changed = entry.receipt.clone().unwrap();
        changed.archive_hash = library::Digest::new("44".repeat(32)).unwrap();
        assert!(pin.accept(changed, &context).is_err());
        let review = pin.accept(entry.receipt.unwrap(), &context).unwrap();
        assert!(
            matches!(review.command(&context).unwrap(),wire::CatalogBackupCommand::RestoreCopy{expected_hash,..} if expected_hash==review.receipt().archive_hash)
        );
    }
    #[test]
    fn unfinished_known_intent_can_be_verified_but_never_restored_directly() {
        let (context, mut entry) = fixture();
        let actual = entry.receipt.take().unwrap();
        let pin = VerificationPin::new(&entry, &context).unwrap();
        assert!(matches!(
            pin.command(),
            wire::CatalogBackupCommand::Verify {
                expected_hash: None,
                ..
            }
        ));
        let mut wrong = actual.clone();
        wrong.id = wire::OperationId::new(format!("op-{:032x}", 9)).unwrap();
        assert!(pin.accept(wrong, &context).is_err());
        let review = pin.accept(actual, &context).unwrap();
        assert!(review.command(&context).is_ok());
    }
    #[test]
    fn local_scope_changes_and_recovery_receipt_substitution_fail_closed() {
        let (context, entry) = fixture();
        let pin = VerificationPin::new(&entry, &context).unwrap();
        let actual = entry.receipt.unwrap();
        let mut foreign = context.clone();
        foreign.lab_id = library::LabId::new(format!("lab-{:032x}", 8)).unwrap();
        assert!(pin.accept(actual.clone(), &foreign).is_err());
        let review = pin.accept(actual.clone(), &context).unwrap();
        assert!(review.command(&foreign).is_err());
        let mut recovered = wire::CatalogRecoveryReceipt {
            operation_id: wire::OperationId::new(format!("op-{:032x}", 3)).unwrap(),
            backup: actual,
            directory_basename: library::Text::new(format!("op-{:032x}", 3)).unwrap(),
            database_basename: library::Text::new("library.sqlite").unwrap(),
        };
        assert!(review.accept(&recovered, &context).is_ok());
        recovered.backup.archive_hash = library::Digest::new("99".repeat(32)).unwrap();
        assert!(review.accept(&recovered, &context).is_err());
    }
}

/// A managed restored copy remains inactive until the native host accepts a
/// checked review. Neither its display name nor its ID alone grants adoption.
#[derive(Clone)]
pub(crate) struct RecoveryChoice {
    context: library::LabContext,
    recovery: wire::CatalogRecoveryReceipt,
}
#[derive(Clone)]
pub(crate) struct AdoptionReview {
    context: library::LabContext,
    review: wire::CatalogRecoveryReview,
}
impl RecoveryChoice {
    pub fn new(
        recovery: wire::CatalogRecoveryReceipt,
        context: &library::LabContext,
    ) -> Result<Self, String> {
        correct_lab(&recovery.backup, context)?;
        if recovery.directory_basename.as_str() != recovery.operation_id.as_str()
            || recovery.database_basename.as_str() != "library.sqlite"
        {
            return Err("Recovery copy does not name a known managed catalog.".into());
        }
        Ok(Self {
            context: context.clone(),
            recovery,
        })
    }
    pub fn command(
        &self,
        current: &library::LabContext,
    ) -> Result<wire::CatalogBackupCommand, String> {
        if current != &self.context {
            return Err("Local catalog changed; select and review the recovery copy again.".into());
        }
        Ok(wire::CatalogBackupCommand::ReviewRecovery {
            id: self.recovery.operation_id.clone(),
            expected_hash: self.recovery.backup.archive_hash.clone(),
        })
    }
    pub fn accept(
        &self,
        review: wire::CatalogRecoveryReview,
        current: &library::LabContext,
    ) -> Result<AdoptionReview, String> {
        review
            .validate()
            .map_err(|_| "Native recovery review has an invalid scope or digest.".to_string())?;
        if current != &self.context
            || review.current_lab != current.lab_id
            || review.target_recovery.as_ref() != Some(&self.recovery)
            || review.target
                != (wire::CatalogLocation::Recovery {
                    id: self.recovery.operation_id.clone(),
                })
        {
            return Err("Recovery review belongs to a different copy or current catalog.".into());
        }
        Ok(AdoptionReview {
            context: current.clone(),
            review,
        })
    }
}
impl AdoptionReview {
    pub fn original(
        review: wire::CatalogRecoveryReview,
        current: &library::LabContext,
    ) -> Result<Self, String> {
        review
            .validate()
            .map_err(|_| "Original-library review has an invalid scope or digest.".to_string())?;
        if review.current_lab != current.lab_id
            || review.target != wire::CatalogLocation::Original
            || review.target_recovery.is_some()
        {
            return Err(
                "Original-library review belongs to another target or catalog scope.".into(),
            );
        }
        Ok(Self {
            context: current.clone(),
            review,
        })
    }
    pub fn review(&self) -> &wire::CatalogRecoveryReview {
        &self.review
    }
    pub fn command(
        &self,
        current: &library::LabContext,
    ) -> Result<wire::CatalogBackupCommand, String> {
        if current != &self.context || current.lab_id != self.review.current_lab {
            return Err(
                "The current catalog scope changed; review recovery adoption again.".into(),
            );
        }
        Ok(wire::CatalogBackupCommand::AdoptRecovery {
            review: self.review.clone(),
        })
    }
    pub fn accept_pending(
        &self,
        receipt: &wire::CatalogActivationReceipt,
        current: &library::LabContext,
    ) -> Result<(), String> {
        if current != &self.context || receipt.review != self.review {
            return Err("Native adoption receipt does not match the exact reviewed copy and current catalog.".into());
        }
        Ok(())
    }
}
/// Cancellation can be reviewed without a healthy catalog or a fabricated lab.
/// It names exactly the accepted journal intent already shown by native status.
#[derive(Clone)]
pub(crate) struct CancellationPin {
    accepted: wire::CatalogActivationReceipt,
}
impl CancellationPin {
    pub fn new(accepted: wire::CatalogActivationReceipt) -> Result<Self, String> {
        accepted
            .review
            .validate()
            .map_err(|_| "Pending activation review is malformed.".to_string())?;
        accepted
            .review
            .activation_generation
            .checked_add(1)
            .ok_or("Pending activation generation overflow.")?;
        Ok(Self { accepted })
    }
    pub fn command(
        &self,
        current: &wire::CatalogActivationReceipt,
    ) -> Result<wire::CatalogBackupCommand, String> {
        if current != &self.accepted {
            return Err(
                "Pending library selection changed; inspect native status before cancelling."
                    .into(),
            );
        }
        Ok(wire::CatalogBackupCommand::CancelActivation {
            accepted_operation: self.accepted.operation_id.clone(),
            review_revision: self.accepted.review.review_revision.clone(),
            activation_generation: self.accepted.review.activation_generation,
        })
    }
    pub fn accept(&self, receipt: &wire::CatalogCancellationReceipt) -> Result<(), String> {
        receipt
            .validate()
            .map_err(|_| "Native cancellation receipt is malformed.".to_string())?;
        if receipt.accepted != self.accepted
            || Some(receipt.generation) != self.accepted.review.activation_generation.checked_add(1)
        {
            return Err(
                "Cancellation receipt does not match the exact pending library selection.".into(),
            );
        }
        Ok(())
    }
}
#[cfg(test)]
mod adoption_tests {
    use super::*;
    fn values() -> (library::LabContext, wire::CatalogRecoveryReceipt) {
        let context = library::LabContext {
            lab_id: library::LabId::new(format!("lab-{:032x}", 1)).unwrap(),
            name: library::Text::new("Current local catalog").unwrap(),
        };
        let backup = wire::CatalogBackupReceipt {
            id: wire::OperationId::new(format!("op-{:032x}", 2)).unwrap(),
            basename: library::Text::new(format!("op-{:032x}.cglbak", 2)).unwrap(),
            archive_hash: library::Digest::new("11".repeat(32)).unwrap(),
            archive_bytes: 4096,
            database_hash: library::Digest::new("22".repeat(32)).unwrap(),
            database_bytes: 2048,
            schema_hash: library::Digest::new("33".repeat(32)).unwrap(),
            schema_version: library::DATABASE_SCHEMA_VERSION,
            lab: context.lab_id.clone(),
            lab_count: 1,
            documents: 3,
            folders: 2,
            revisions: 4,
            mutation_receipts: 0,
            created_unix_seconds: 1795680000,
            request_checkpoint: None,
        };
        (
            context,
            wire::CatalogRecoveryReceipt {
                operation_id: wire::OperationId::new(format!("op-{:032x}", 3)).unwrap(),
                backup,
                directory_basename: library::Text::new(format!("op-{:032x}", 3)).unwrap(),
                database_basename: library::Text::new("library.sqlite").unwrap(),
            },
        )
    }
    #[test]
    fn explicit_adoption_is_pinned_to_checked_copy_scope_generation_and_receipt() {
        let (context, recovery) = values();
        let choice = RecoveryChoice::new(recovery.clone(), &context).unwrap();
        assert!(matches!(
            choice.command(&context).unwrap(),
            wire::CatalogBackupCommand::ReviewRecovery { .. }
        ));
        let mut review = wire::CatalogRecoveryReview {
            target: wire::CatalogLocation::Recovery {
                id: recovery.operation_id.clone(),
            },
            target_recovery: Some(recovery.clone()),
            current_catalog: wire::CatalogLocation::Original,
            activation_generation: 0,
            current_database_hash: library::Digest::new("77".repeat(32)).unwrap(),
            target_snapshot: wire::CatalogSnapshot {
                database_hash: recovery.backup.database_hash.clone(),
                database_bytes: recovery.backup.database_bytes,
                schema_hash: recovery.backup.schema_hash.clone(),
                schema_version: library::DATABASE_SCHEMA_VERSION,
                documents: 3,
                folders: 2,
                revisions: 4,
                mutation_receipts: 0,
                catalog_generation: 5,
                checkpoint: None,
            },
            current_lab: context.lab_id.clone(),
            current_generation: 7,
            current_checkpoint: Some(library::Digest::new("44".repeat(32)).unwrap()),
            review_revision: library::Digest::new("55".repeat(32)).unwrap(),
        };
        review.review_revision = review.revision().unwrap();
        let approved = choice.accept(review.clone(), &context).unwrap();
        assert!(matches!(
            approved.command(&context).unwrap(),
            wire::CatalogBackupCommand::AdoptRecovery { .. }
        ));
        let receipt = wire::CatalogActivationReceipt {
            operation_id: wire::OperationId::new(format!("op-{:032x}", 4)).unwrap(),
            review,
        };
        approved.accept_pending(&receipt, &context).unwrap();
        let cancellation = CancellationPin::new(receipt.clone()).unwrap();
        assert!(
            matches!(cancellation.command(&receipt).unwrap(), wire::CatalogBackupCommand::CancelActivation { accepted_operation, .. } if accepted_operation == receipt.operation_id)
        );
        let cancelled = wire::CatalogCancellationReceipt {
            operation_id: wire::OperationId::new(format!("op-{:032x}", 8)).unwrap(),
            accepted: receipt.clone(),
            generation: receipt.review.activation_generation + 1,
        };
        cancellation.accept(&cancelled).unwrap();
        let mut wrong = cancelled.clone();
        wrong.generation += 1;
        assert!(cancellation.accept(&wrong).is_err());
        let mut other = receipt.clone();
        other.operation_id = wire::OperationId::new(format!("op-{:032x}", 9)).unwrap();
        assert!(cancellation.command(&other).is_err());
        wrong = cancelled.clone();
        wrong.accepted = other;
        assert!(cancellation.accept(&wrong).is_err());
        let mut changed = receipt.clone();
        changed.review.current_generation += 1;
        assert!(approved.accept_pending(&changed, &context).is_err());
        let mut changed = receipt.clone();
        changed.review.review_revision = library::Digest::new("66".repeat(32)).unwrap();
        assert!(approved.accept_pending(&changed, &context).is_err());
        let mut old = receipt.clone();
        old.review.activation_generation += 1;
        assert!(approved.accept_pending(&old, &context).is_err());
        let mut rollback = approved.review().clone();
        rollback.current_catalog = rollback.target.clone();
        rollback.target = wire::CatalogLocation::Original;
        rollback.target_recovery = None;
        rollback.activation_generation = 1;
        rollback.review_revision = rollback.revision().unwrap();
        let original = AdoptionReview::original(rollback.clone(), &context).unwrap();
        assert!(
            matches!(original.command(&context).unwrap(),wire::CatalogBackupCommand::AdoptRecovery{review} if review.target==wire::CatalogLocation::Original)
        );
        assert!(choice.accept(rollback, &context).is_err());
        let mut foreign = context.clone();
        foreign.lab_id = library::LabId::new(format!("lab-{:032x}", 9)).unwrap();
        assert!(approved.command(&foreign).is_err());
    }
}
