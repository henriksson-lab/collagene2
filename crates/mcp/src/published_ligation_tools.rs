//! Published data and mathematical ranking remain separate from reaction authority.
use super::*;
use collagene_core::{
    golden_gate_fusion_rank as ranking, ligation_profiles as profiles, topo_preparation,
};
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListQuery {}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScoreQuery {
    pub profile: profiles::ProfileRef,
    pub ends: Vec<String>,
    #[serde(default)]
    pub physical_competitors: Vec<String>,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RankQuery {
    pub inputs: Vec<BoundInput>,
    pub request: ranking::Request,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoQuery {
    pub sources: Vec<TopoSourceRef>,
    pub declaration: topo_preparation::PreparationRequest,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NickQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: collagene_core::topo_nick_sealing::Request,
}
impl CollageneServer {
    pub(super) async fn topo_seal_nicks(
        &self,
        args: NickQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            if args.request.sealing.nicks.len()>16 || args.request.sealing.ligase_assumption.len()>4096 {
                return Err(bad("Bounded explicit nick and ligase declarations required"));
            }
            let (inputs,_) = server.topo_sources(&args.sources,&args.request.preparation,None)?;
            let report=collagene_core::topo_nick_sealing::inspect(&inputs,&args.request).map_err(bad)?;
            report.validate_binding(&collagene_core::topo_nick_sealing::revision(&inputs,&args.request).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"report":report}),32*1024*1024)?,"application/json","topo-nick-sealing.json")?;
            output(json!({"state":"declared-downstream-graph","verified":false,"context_revision":report.context_revision,
                "sources":args.sources,"steps":report.sealing.steps,"physical":json!({"chains":report.sealing.product.chains.len(),"components":report.sealing.product.components.len(),"nicks":report.sealing.product.nicks.len(),"unpaired_runs":report.sealing.product.unpaired.len(),"released_activations":report.sealing.product.released_activations}),
                "artifact_id":artifact,"authority":"Fresh explicit ligase-step reconstruction. Every original strand/component/displaced call is retained. No kinase, flap deletion, host repair, graph-to-duplex authority or independent-target acceptance is inferred."}))
        }).await
    }
    pub(super) async fn published_profiles_list(
        &self,
        _args: ListQuery,
    ) -> Result<CallToolResult, ErrorData> {
        self.compute(move||output(json!({"catalog_revision":profiles::catalog_revision().map_err(bad)?,
            "profiles":profiles::list().map_err(bad)?,"verified":false,
            "authority":"Published experimental frequency data under its exact cited conditions, separately attributed CC-BY data. Profiles do not establish actual supplied enzyme formulation, reaction yield or physical compatibility."}))).await
    }
    pub(super) async fn published_profile_score(
        &self,
        args: ScoreQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let report=profiles::score_with_competitors(&args.profile,&args.ends,&args.physical_competitors).map_err(bad)?;
            report.validate_binding(&profiles::score_revision_with_competitors(&args.profile,&args.ends,&args.physical_competitors).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&report,16*1024*1024)?,"application/json","published-ligation-frequency-score.json")?;
            output(json!({"state":"conditional-published-frequency-estimate","verified":false,"report":report,"artifact_id":artifact}))
        }).await
    }
    pub(super) async fn published_fusion_rank(
        &self,
        args: RankQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let inputs=server.fusion_sources(&args.inputs,&args.request.design)?;
            let report=ranking::inspect(&inputs,&args.request).map_err(bad)?;
            report.validate_binding(&ranking::inspect_revision(&inputs,&args.request).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"inputs":args.inputs,"report":report}),32*1024*1024)?,"application/json","ranked-golden-gate-fusions.json")?;
            output(json!({"state":"advisory-fusion-ranking","verified":false,"context_revision":report.context_revision,
                "selection_context_revision":report.search.context_revision,"metadata":report.metadata,
                "candidate_count":report.ranked.len(),"ranked":report.ranked.iter().take(20).collect::<Vec<_>>(),
                "search_complete":report.search.search_complete,"limits":report.limits,"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Ranks only the retained bounded mathematical fusion sets under the explicit published experiment. No recipe, intended DNA or physical candidate is chosen or changed. Explicitly materialize one original candidate ID using selection_context_revision and freshly verify its full preparation."}))
        }).await
    }
    pub(super) async fn topo_compile_preparation(
        &self,
        args: TopoQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let initial=topo_preparation::compile(&args.declaration).map_err(bad)?;
            let (inputs,_)=server.topo_sources(&args.sources,&initial.request,None)?;
            let compiled=topo_preparation::compile_preparation(&inputs,&args.declaration).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"declaration":args.declaration,"compiled":compiled}),32*1024*1024)?,"application/json","topo-preparation-declaration.json")?;
            output(json!({"state":"declared-preparation","verified":false,"sources":args.sources,"compiled":compiled,"artifact_id":artifact,
                "authority":"Source-backed preparation declaration only. Use inspect_topo_orientations with the returned request, strand indexes and directional flag; inspect actual transfers/components before choosing a source-call seed. A displaced flap can withhold a paired contour. No activation, repair or closed clone is inferred."}))
        }).await
    }
}
