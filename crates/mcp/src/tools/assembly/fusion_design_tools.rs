use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::tools::assembly::intents::BoundInput;
use collagene_core::planning::plan::PlanInput;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;
use std::collections::BTreeSet;

use collagene_core::cloning::{
    golden_gate::fusion_design as fusion, ligation::fidelity as ligation_fidelity,
};

#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SearchQuery {
    pub inputs: Vec<BoundInput>,
    pub request: fusion::Request,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterializeQuery {
    pub inputs: Vec<BoundInput>,
    pub request: fusion::Request,
    pub selection: fusion::Selection,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FidelityQuery {
    pub dataset: ligation_fidelity::Dataset,
    pub ends: Vec<String>,
    #[serde(default)]
    pub physical_competitors: Option<Vec<String>>,
}
impl CollageneServer {
    pub(crate) fn fusion_sources(
        &self,
        bindings: &[BoundInput],
        request: &fusion::Request,
    ) -> Result<Vec<PlanInput>, ErrorData> {
        if bindings.is_empty() || bindings.len() > 32 {
            return Err(bad("choose1..32 original source bindings"));
        }
        let store = self.session.store.lock().unwrap();
        let mut seen = BTreeSet::new();
        let mut sources = Vec::with_capacity(bindings.len());
        for binding in bindings {
            if binding.id.trim().is_empty()
                || binding.id.len() > 128
                || !seen.insert(binding.id.as_str())
            {
                return Err(bad("invalid or duplicate original source ID"));
            }
            sources.push((binding.id.as_str(), store.document(&binding.document)?));
        }
        fusion::preflight(&sources, request).map_err(bad)?;
        Ok(sources
            .into_iter()
            .map(|(id, molecule)| PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule: molecule.clone(),
            })
            .collect())
    }
    pub(crate) async fn fusion_search(
        &self,
        args: SearchQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let inputs = server.fusion_sources(&args.inputs, &args.request)?;
            let report = fusion::inspect(&inputs, &args.request).map_err(bad)?;
            report.validate_binding(&fusion::inspect_revision(&inputs, &args.request).map_err(bad)?).map_err(bad)?;
            let artifact = server.session.text_artifact(
                bounded_json_limit(&json!({"inputs":args.inputs,"report":report}),32*1024*1024)?,
                "application/json", "golden-gate-fusions.json",
            )?;
            output(json!({"state":"advisory-placement-search","verified":false,
                "context_revision":report.context_revision,"candidate_count":report.candidates.len(),
                "candidates":report.candidates.iter().take(20).collect::<Vec<_>>(),
                "states_examined":report.states_examined,"search_complete":report.search_complete,
                "target_revision":report.original_intent.target.revision(),"warnings":report.warnings,
                "artifact_id":artifact,
                "authority":"Mathematical fusion constraints only. Explicitly materialize one candidate, then freshly reconstruct its complete recipe against the unchanged intended target. No empirical ranking, PCR feasibility or experimental success is established by this search."}))
        }).await
    }
    pub(crate) async fn fusion_materialize(
        &self,
        args: MaterializeQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let inputs = server.fusion_sources(&args.inputs, &args.request)?;
            let report = fusion::materialize(&inputs, &args.request, &args.selection).map_err(bad)?;
            report.validate_binding(&fusion::materialize_revision(&inputs, &args.request, &args.selection).map_err(bad)?).map_err(bad)?;
            let artifact = server.session.text_artifact(
                bounded_json_limit(&json!({"inputs":args.inputs,"report":report}),32*1024*1024)?,
                "application/json", "golden-gate-fusion-preparations.json",
            )?;
            output(json!({"state":"advisory-preparation","verified":false,
                "context_revision":report.context_revision,"parent_context":report.parent_context,
                "target_revision":report.original_intent.target.revision(),"selection":report.selection,
                "recipe":report.recipe,"placements":report.placements,"warnings":report.warnings,
                "artifact_id":artifact,
                "authority":"Original-source primers and pinned digest preparations for the explicit fusion choice. The independent intended DNA is unchanged. Freshly check the full recipe with propose_ligation_assembly; this does not establish cycling-mixture uniqueness or experimental success."}))
        }).await
    }
    pub(crate) async fn fidelity_score(
        &self,
        args: FidelityQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let report = if let Some(extras) = &args.physical_competitors {
                let report = ligation_fidelity::score_competition(&args.dataset, &args.ends, extras).map_err(bad)?;
                report.validate_binding(&ligation_fidelity::competition_revision(&args.dataset, &args.ends, extras).map_err(bad)?).map_err(bad)?;
                serde_json::to_value(report).map_err(bad)?
            } else {
                let report = ligation_fidelity::score(&args.dataset, &args.ends).map_err(bad)?;
                report.validate_binding(&ligation_fidelity::revision(&args.dataset, &args.ends).map_err(bad)?).map_err(bad)?;
                serde_json::to_value(report).map_err(bad)?
            };
            let artifact = server.session.text_artifact(
                bounded_json_limit(&json!({"dataset":args.dataset,"report":report}),16*1024*1024)?,
                "application/json", "ligation-fidelity.json",
            )?;
            output(json!({"state":"conditional-empirical-estimate","verified":false,
                "report":report,"artifact_id":artifact,
                "authority":"Declared directed ligation-count model under the supplied dataset conditions. Missing data remains unknown. This estimate does not validate physical compatibility, measured reaction yield or experimental success."}))
        }).await
    }
}
