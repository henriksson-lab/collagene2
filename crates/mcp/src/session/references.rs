use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DocumentRef {
    /// Immutable session document ID, returned by import_document.
    pub id: String,
    /// Exact expected revision, returned by import_document or list_documents.
    pub revision: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CollectionRef {
    pub id: String,
    pub revision: String,
}
/// A reference to immutable session evidence; exact revisions are always required.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EvidenceRef {
    pub id: String,
    pub revision: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceUploadRef {
    pub upload_id: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceRef {
    pub id: String,
    pub revision: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlanRef {
    /// Exact immutable plan revision returned by propose_plan.
    pub plan_revision: String,
}

/// MCP callers pass small immutable event selections, never closed capsule bytes.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum MaterialDerivationOrigin {
    OriginalReference,
    Saved(collagene_core::material::run::MaterialSelection),
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialOriginWorkspace {
    pub role: String,
    pub workspace: WorkspaceRef,
}
