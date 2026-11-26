//! Metadata-only local remote retirement. No vault secret or biological record.
use crate::{Error, PublicationRemote, ReviewId, library};
use library::{Digest, LabId};
use serde::{Deserialize, Serialize};
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct RemovalDescription {pub local_lab:LabId,pub remote:PublicationRemote,pub expected_checkpoint:Option<Digest>,pub document_bindings:u32,pub folder_bindings:u32,pub completed_tree_jobs:u16}}
dto! {pub struct RemoteRemovalReview {pub id:ReviewId,pub description:RemovalDescription,pub digest:Digest}}
dto! {pub struct RemovalReceipt {pub schema_version:u16,pub local_lab:LabId,pub remote:PublicationRemote,pub review_digest:Digest,pub removed_document_bindings:u32,pub removed_folder_bindings:u32,pub removed_completed_tree_jobs:u16}}
impl RemovalDescription {
    pub fn validate(&self) -> Result<(), Error> {
        self.remote.validate()?;
        if self.document_bindings > 100000
            || self.folder_bindings > 100000
            || self.completed_tree_jobs > 128
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
impl RemovalReceipt {
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(Error::Invalid("removal receipt version"));
        }
        RemovalDescription {
            local_lab: self.local_lab.clone(),
            remote: self.remote.clone(),
            expected_checkpoint: None,
            document_bindings: self.removed_document_bindings,
            folder_bindings: self.removed_folder_bindings,
            completed_tree_jobs: self.removed_completed_tree_jobs,
        }
        .validate()
    }
}
