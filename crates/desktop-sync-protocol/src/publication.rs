//! Bounded local review data. No biological archive, credential or path crosses DSQ.
use crate::{Error, RemoteId, ReviewId, library};
use library::{Digest, DocumentId, FolderId, LabId, Rows, Text};
use serde::{Deserialize, Serialize};
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum PublicationTarget {Folder(FolderId),DocumentMetadata(DocumentId),DeletedFolder(FolderId)}}
dto! {pub struct PublicationRemote {pub id:RemoteId,pub name:Text<256>,pub endpoint:Text<2048>,pub lab:LabId}}
dto! {pub enum MetadataAction {
 CreateFolder{parent:Option<FolderId>,name:Text<1024>},
 RenameFolder{id:FolderId,expected_generation:u64,name:Text<1024>},
 MoveFolder{id:FolderId,expected_generation:u64,parent:Option<FolderId>},
 RenameDocument{id:DocumentId,expected_generation:u64,name:Text<1024>},
 MoveDocument{id:DocumentId,expected_generation:u64,folder:Option<FolderId>},
 DeleteEmptyFolder{id:FolderId,expected_generation:u64},
}}
dto! {pub struct PublicationDescription {
 pub local_lab:LabId,pub target:PublicationTarget,pub remote:PublicationRemote,
 pub expected_checkpoint:Option<Digest>,pub local_generation:Option<u64>,
 pub local_name:Option<Text<1024>>,pub local_parent:Option<FolderId>,
 pub remote_parent:Option<FolderId>,pub parent_binding:Option<library::Folder>,pub action:MetadataAction,
}}
dto! {pub struct MetadataReview {pub description:PublicationDescription,pub digest:Digest}}
dto! {pub struct PushChoice {pub original_hash:Digest,pub component_iri:Text<4096>,pub sequences:Rows<Text<4096>,128>}}
dto! {pub struct PushChoiceRequirement {pub original_hash:Digest,pub component_iri:Text<4096>,pub available_sequences:Rows<Text<4096>,128>,pub required_sequences:Rows<Text<4096>,128>}}
dto! {pub struct PushChoices {pub local_document:DocumentId,pub local_revision:library::RevisionId,pub remote:RemoteId,pub requirements:Rows<PushChoiceRequirement,128>}}
dto! {pub struct PushReviewEntry {pub original_hash:Digest,pub subset_hash:Digest,pub retained_objects:u32,pub omitted_objects:u32,pub omitted_sequence_edges:u32}}
dto! {pub struct PushPublicationReview {pub id:ReviewId,pub local_document:DocumentId,pub local_revision:library::RevisionId,pub remote:RemoteId,pub configuration:PublicationRemote,pub checkpoint:Option<Digest>,pub entries:Rows<PushReviewEntry,128>}}

impl PublicationDescription {
    pub fn validate(&self) -> Result<(), Error> {
        super::name(&self.remote.name)?;
        super::endpoint(&self.remote.endpoint)?;
        if let Some(name) = &self.local_name {
            library::validate_name(name)?;
        }
        if let Some(parent) = &self.parent_binding {
            parent.validate()?;
            if self.remote_parent.as_ref() != Some(&parent.id) {
                return Err(Error::Mismatch);
            }
        }
        match &self.action {
            MetadataAction::CreateFolder { name, .. }
            | MetadataAction::RenameFolder { name, .. }
            | MetadataAction::RenameDocument { name, .. } => library::validate_name(name)?,
            _ => (),
        }
        Ok(())
    }
}
pub fn validate_choices(choices: &Rows<PushChoice, 128>) -> Result<(), Error> {
    let mut seen = std::collections::BTreeSet::new();
    for c in choices.iter() {
        if c.component_iri.is_empty()
            || c.sequences.is_empty()
            || !seen.insert((&c.original_hash, &c.component_iri))
        {
            return Err(Error::Invalid("sequence choices"));
        }
        let mut sequences = std::collections::BTreeSet::new();
        for s in c.sequences.iter() {
            if s.is_empty() || !sequences.insert(s) {
                return Err(Error::Invalid("sequence alternatives"));
            }
        }
    }
    Ok(())
}

dto! {pub struct DeletedFolderPublication {pub local_folder:FolderId,pub remote_folder:library::Folder}}
dto! {pub struct DeletedFolderPage {pub remote:RemoteId,pub items:Rows<DeletedFolderPublication,200>,pub next_offset:Option<u16>,pub checkpoint:Option<Digest>}}

impl PublicationRemote {
    pub fn validate(&self) -> Result<(), Error> {
        super::name(&self.name)?;
        super::endpoint(&self.endpoint)
    }
}
