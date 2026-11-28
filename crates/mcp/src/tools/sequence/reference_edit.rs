//! Immutable reference edits use the canonical core operation owner.
use crate::session::references::DocumentRef;
use collagene_core::documents::reference_edit::Operation;
use schemars::JsonSchema;
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub document: DocumentRef,
    pub operation: Operation,
}
