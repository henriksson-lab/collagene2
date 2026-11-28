//! Explicit catalog/selection/options-bound documentary reaction assessment.
use super::catalog_tools;
use crate::server::{
    CollageneServer,
    response::{bad, bounded_json_limit, output},
};
use collagene_core::cloning::restriction::{
    catalog::CatalogSelection,
    conditions::{self, ReactionConditions},
};
use rmcp::{ErrorData, model::CallToolResult};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

const AUTHORITY: &str = "Documentary formulation/assay compatibility under the exact declared options only. No sequence-specific cleavage, methylation permission, reaction yield, star-cleavage probability or experimental success is established.";

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssessmentRequest {
    pub selection: CatalogSelection,
    pub options: ReactionConditions,
}

pub(super) async fn assess(
    server: &CollageneServer,
    args: AssessmentRequest,
) -> Result<CallToolResult, ErrorData> {
    let owner = server.clone();
    server
        .compute(move || {
            let catalog = catalog_tools::builtin()?;
            let report =
                conditions::assess(&catalog, &args.selection, &args.options).map_err(bad)?;
            let revision =
                conditions::revision(&catalog, &args.selection, &args.options).map_err(bad)?;
            let artifact_id = owner.session.text_artifact(
                bounded_json_limit(&report, conditions::MAX_REPORT_BYTES)?,
                "application/json",
                "restriction-conditions.json",
            )?;
            output(json!({
                "identity":catalog_tools::identity(&catalog),
                "revision":revision,
                "batch_revision":report.batch_revision,
                "documentary_status":report.documentary_status,
                "selected_enzyme_count":args.selection.enzymes.len(),
                "common_condition_count":report.common_conditions.len(),
                "suggestions_complete":report.suggestions_complete,
                "artifact_id":artifact_id,"encoding":"utf-8",
                "authority":AUTHORITY,"sequence_assessed":false,
                "experimental_success_established":false,"persistence":"session only"
            }))
        })
        .await
}
