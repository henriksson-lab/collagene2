//! Saved-catalog recovery metadata only; frontend requests contain no paths or archives.
use crate::{Error, OperationId, library};
use library::{Digest, LabId, Rows, Text};
use serde::{Deserialize, Serialize};
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
pub const MAX_CATALOG_BACKUP_REQUEST_BYTES: usize = 16 << 10;
pub const MAX_CATALOG_BACKUP_RESPONSE_BYTES: usize = 256 << 10;
dto! {pub enum CatalogBackupCommand {
 Create,List{offset:u32,limit:u16,expected_catalog:Option<Digest>},
 Verify{id:OperationId,expected_hash:Option<Digest>},
 RestoreCopy{id:OperationId,expected_hash:Digest},
 ListRecoveries{offset:u32,limit:u16,expected_catalog:Option<Digest>},
 ReviewRecovery{id:OperationId,expected_hash:Digest}, ReviewOriginal, AdoptRecovery{review:CatalogRecoveryReview}, ActivationStatus,
 CancelActivation{accepted_operation:OperationId,review_revision:Digest,activation_generation:u64},
}}
dto! {pub struct CatalogBackupReceipt {
 pub id:OperationId,pub basename:Text<128>,pub archive_hash:Digest,pub archive_bytes:u64,
 pub database_hash:Digest,pub database_bytes:u64,pub schema_hash:Digest,pub schema_version:u32,
 pub lab:LabId,pub lab_count:u32,pub documents:u64,pub folders:u64,pub revisions:u64,
 pub mutation_receipts:u64,pub created_unix_seconds:u64,pub request_checkpoint:Option<Digest>,
}}
dto! {pub struct CatalogBackupEntry {pub id:OperationId,pub basename:Text<128>,pub receipt:Option<CatalogBackupReceipt>}}
dto! {pub struct CatalogBackupPage {pub entries:Rows<CatalogBackupEntry,200>,pub next_offset:Option<u32>,pub catalog_revision:Digest}}
dto! {pub struct CatalogRecoveryReceipt {pub operation_id:OperationId,pub backup:CatalogBackupReceipt,pub directory_basename:Text<128>,pub database_basename:Text<128>}}
dto! {pub struct CatalogRecoveryPage {pub entries:Rows<CatalogRecoveryReceipt,200>,pub next_offset:Option<u32>,pub catalog_revision:Digest}}
dto! {pub enum CatalogLocation {Original,Recovery{id:OperationId}}}
dto! {pub struct CatalogSnapshot {pub database_hash:Digest,pub database_bytes:u64,pub schema_hash:Digest,pub schema_version:u32,pub documents:u64,pub folders:u64,pub revisions:u64,pub mutation_receipts:u64,pub catalog_generation:u64,pub checkpoint:Option<Digest>}}
dto! {pub struct CatalogRecoveryReview {pub target:CatalogLocation,pub target_recovery:Option<CatalogRecoveryReceipt>,pub current_catalog:CatalogLocation,pub activation_generation:u64,pub current_lab:LabId,pub current_generation:u64,pub current_checkpoint:Option<Digest>,pub current_database_hash:Digest,pub target_snapshot:CatalogSnapshot,pub review_revision:Digest}}
dto! {pub struct CatalogActivationReceipt {pub operation_id:OperationId,pub review:CatalogRecoveryReview}}
dto! {pub struct CatalogCancellationReceipt {pub operation_id:OperationId,pub accepted:CatalogActivationReceipt,pub generation:u64}}
dto! {pub enum CatalogActivationStatus {None,PendingRestart(CatalogActivationReceipt),ActiveAfterRestart(CatalogActivationReceipt),Cancelled(CatalogCancellationReceipt)}}
impl CatalogCancellationReceipt {
    pub fn validate(&self) -> Result<(), Error> {
        self.accepted.review.validate()?;
        if self.generation != self.accepted.review.activation_generation + 1
            || self.operation_id == self.accepted.operation_id
        {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}
impl CatalogActivationStatus {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::None => Ok(()),
            Self::PendingRestart(r) | Self::ActiveAfterRestart(r) => r.review.validate(),
            Self::Cancelled(r) => r.validate(),
        }
    }
}
impl CatalogRecoveryReview {
    pub fn revision(&self) -> Result<Digest, Error> {
        let mut value = self.clone();
        value.review_revision = Digest::of(b"catalog-activation-review-v3");
        let bytes = postcard::to_allocvec(&value).map_err(|_| Error::Malformed)?;
        Ok(Digest::of(&bytes))
    }
    pub fn validate(&self) -> Result<(), Error> {
        match (&self.target, &self.target_recovery) {
            (CatalogLocation::Original, None) => {}
            (CatalogLocation::Recovery { id }, Some(r)) => {
                CatalogBackupSuccess::Restored(r.clone()).validate()?;
                if id != &r.operation_id || self.current_lab != r.backup.lab {
                    return Err(Error::Mismatch);
                }
            }
            _ => return Err(Error::Mismatch),
        }
        if self.target == self.current_catalog
            || self.activation_generation >= 4096
            || self.target_snapshot.database_bytes == 0
            || self.target_snapshot.schema_version != library::DATABASE_SCHEMA_VERSION
            || self.revision()? != self.review_revision
        {
            return Err(Error::Mismatch);
        }
        Ok(())
    }
}
dto! {pub enum CatalogBackupSuccess {Created(CatalogBackupReceipt),Listed(CatalogBackupPage),Verified(CatalogBackupReceipt),Restored(CatalogRecoveryReceipt),Reviewed(CatalogRecoveryReview),AdoptionPending(CatalogActivationReceipt),ActivationStatus(CatalogActivationStatus),RecoveriesListed(CatalogRecoveryPage),ActivationCancelled(CatalogCancellationReceipt)}}
impl CatalogBackupCommand {
    pub fn validate(&self) -> Result<(), Error> {
        if let Self::List {
            offset,
            limit,
            expected_catalog,
        }
        | Self::ListRecoveries {
            offset,
            limit,
            expected_catalog,
        } = self
        {
            if *limit == 0
                || *limit > 200
                || *offset > 4096
                || (*offset > 0 && expected_catalog.is_none())
            {
                return Err(Error::Invalid("catalog backup page"));
            }
        }
        if let Self::CancelActivation {
            activation_generation,
            ..
        } = self
        {
            if *activation_generation >= 4096 {
                return Err(Error::Invalid("activation generation"));
            }
        }
        if let Self::AdoptRecovery { review } = self {
            review.validate()?;
        }
        Ok(())
    }
}
impl CatalogBackupReceipt {
    pub fn validate(&self) -> Result<(), Error> {
        crate::validate_export_name(&self.basename)?;
        if self.basename.as_str() != format!("{}.cglbak", self.id.as_str())
            || self.lab_count != 1
            || self.archive_bytes == 0
            || self.database_bytes == 0
            || self.archive_bytes <= self.database_bytes
            || self.schema_version != library::DATABASE_SCHEMA_VERSION
        {
            return Err(Error::Invalid("catalog backup receipt"));
        }
        Ok(())
    }
}
impl CatalogBackupSuccess {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::ActivationCancelled(r) => r.validate()?,
            Self::ActivationStatus(CatalogActivationStatus::Cancelled(r)) => r.validate()?,
            Self::RecoveriesListed(p) => {
                let mut seen = std::collections::BTreeSet::new();
                for r in p.entries.iter() {
                    CatalogBackupSuccess::Restored(r.clone()).validate()?;
                    if !seen.insert(&r.operation_id) {
                        return Err(Error::Mismatch);
                    }
                }
            }
            Self::Reviewed(r) => r.validate()?,
            Self::AdoptionPending(r) => r.review.validate()?,
            Self::ActivationStatus(CatalogActivationStatus::None) => {}
            Self::ActivationStatus(
                CatalogActivationStatus::PendingRestart(r)
                | CatalogActivationStatus::ActiveAfterRestart(r),
            ) => r.review.validate()?,
            Self::Created(r) | Self::Verified(r) => r.validate()?,
            Self::Restored(r) => {
                r.backup.validate()?;
                crate::validate_export_name(&r.directory_basename)?;
                if r.directory_basename.as_str() != r.operation_id.as_str()
                    || r.database_basename.as_str() != "library.sqlite"
                {
                    return Err(Error::Invalid("recovery receipt"));
                }
            }
            Self::Listed(p) => {
                let mut seen = std::collections::BTreeSet::new();
                for e in p.entries.iter() {
                    if !seen.insert(&e.id)
                        || e.basename.as_str() != format!("{}.cglbak", e.id.as_str())
                    {
                        return Err(Error::Invalid("backup listing"));
                    }
                    if let Some(r) = &e.receipt {
                        r.validate()?;
                        if r.id != e.id || r.basename != e.basename {
                            return Err(Error::Mismatch);
                        }
                    }
                }
            }
        }
        Ok(())
    }
    pub fn validate_for(
        &self,
        c: &CatalogBackupCommand,
        operation: &OperationId,
    ) -> Result<(), Error> {
        self.validate()?;
        let matches = match (c, self) {
            (
                CatalogBackupCommand::CancelActivation {
                    accepted_operation,
                    review_revision,
                    activation_generation,
                },
                Self::ActivationCancelled(r),
            ) => {
                &r.operation_id == operation
                    && &r.accepted.operation_id == accepted_operation
                    && &r.accepted.review.review_revision == review_revision
                    && r.accepted.review.activation_generation == *activation_generation
            }
            (CatalogBackupCommand::ReviewOriginal, Self::Reviewed(r)) => {
                r.target == CatalogLocation::Original
            }
            (
                CatalogBackupCommand::ListRecoveries {
                    offset,
                    limit,
                    expected_catalog,
                },
                Self::RecoveriesListed(p),
            ) => {
                p.entries.len() <= usize::from(*limit)
                    && p.next_offset.is_none_or(|n| n > *offset)
                    && (*offset == 0 || expected_catalog.as_ref() == Some(&p.catalog_revision))
            }
            (CatalogBackupCommand::ReviewRecovery { id, expected_hash }, Self::Reviewed(r)) => r
                .target_recovery
                .as_ref()
                .is_some_and(|v| id == &v.operation_id && expected_hash == &v.backup.archive_hash),
            (CatalogBackupCommand::AdoptRecovery { review }, Self::AdoptionPending(r)) => {
                review == &r.review && operation == &r.operation_id
            }
            (CatalogBackupCommand::ActivationStatus, Self::ActivationStatus(_)) => true,
            (CatalogBackupCommand::Create, Self::Created(r)) => &r.id == operation,
            (CatalogBackupCommand::Verify { id, expected_hash }, Self::Verified(r)) => {
                id == &r.id && expected_hash.as_ref().is_none_or(|h| h == &r.archive_hash)
            }
            (CatalogBackupCommand::RestoreCopy { id, expected_hash }, Self::Restored(r)) => {
                &r.operation_id == operation
                    && id == &r.backup.id
                    && expected_hash == &r.backup.archive_hash
            }
            (
                CatalogBackupCommand::List {
                    offset,
                    limit,
                    expected_catalog,
                },
                Self::Listed(p),
            ) => {
                p.entries.len() <= usize::from(*limit)
                    && p.next_offset.is_none_or(|n| n > *offset)
                    && (*offset == 0 || expected_catalog.as_ref() == Some(&p.catalog_revision))
            }
            _ => false,
        };
        if matches {
            Ok(())
        } else {
            Err(Error::Mismatch)
        }
    }
}
