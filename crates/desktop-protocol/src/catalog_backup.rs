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
dto! {pub enum CatalogBackupSuccess {Created(CatalogBackupReceipt),Listed(CatalogBackupPage),Verified(CatalogBackupReceipt),Restored(CatalogRecoveryReceipt)}}
impl CatalogBackupCommand {
    pub fn validate(&self) -> Result<(), Error> {
        if let Self::List {
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
            || self.schema_version == 0
        {
            return Err(Error::Invalid("catalog backup receipt"));
        }
        Ok(())
    }
}
impl CatalogBackupSuccess {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
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
