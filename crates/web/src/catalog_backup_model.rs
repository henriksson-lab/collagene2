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
            schema_version: 1,
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
