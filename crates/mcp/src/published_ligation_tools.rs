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
    pub declaration: topo_graph_tools::Declaration,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
    #[serde(deserialize_with = "Option::deserialize")]
    pub template: Option<topo_graph_tools::BoundTemplate>,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NickQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Sealing,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProcessedQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Processed,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
impl CollageneServer {
    pub(super) async fn topo_processed_product(
        &self,
        args: ProcessedQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            bounded_json_limit(&args.request,2*1024*1024)?;
            let (preparation,_)=server.hydrate_topo(&args.request.treatment.preparation,args.origin_workspace.as_ref(),&[])?;
            let native=collagene_core::topo_processed_product::Request{treatment:collagene_core::topo_nick_sealing::Request{preparation,sealing:args.request.treatment.sealing.clone()},selection:args.request.selection.clone(),options:args.request.options.clone()};

            if args.request.treatment.sealing.nicks.len()>16 || args.request.treatment.sealing.ligase_assumption.len()>4096 {
                return Err(bad("Bounded source-backed treatment required"));
            }
            let (inputs,_) = server.topo_sources(&args.sources,&native.treatment.preparation,None)?;
            let report=collagene_core::topo_processed_product::inspect(&inputs,&native).map_err(bad)?;
            report.validate_binding(&collagene_core::topo_processed_product::revision(&inputs,&native).map_err(bad)?).map_err(bad)?;
            let summary=report.duplex.as_ref().map(|d|d.summary()).transpose().map_err(bad)?;
            let (artifact,document)=server.candidate_document_artifact_bounded(&json!({"sources":args.sources,"origin_workspace":args.origin_workspace,"report":report}),report.product.clone(),"processed-topo-product.json",32*1024*1024)?;
            output(json!({"state":"source-reconstructed-processed-contour","verified":false,"context_revision":report.context_revision,
                "contour_proven":matches!(report.projection.outcome,collagene_core::topo::physical::projection::Outcome::Proven(_)),"withheld":match &report.projection.outcome { collagene_core::topo::physical::projection::Outcome::Withheld(reason)=>Some(reason),_=>None },"duplex":summary,"bridge_notes":report.bridge_notes,"document":document,
                "artifact_id":artifact,"authority":"Fresh original-source transfer, explicit ligase treatment and selected-component paired-contour proof. Full parent/treated graphs and other components remain in evidence. Internal modifications remain unknown; an exported sequence document cannot restore physical state. No independent-target acceptance, isolation or experimental mixture winner is established."}))
        }).await
    }
    pub(super) async fn topo_seal_nicks(
        &self,
        args: NickQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            bounded_json_limit(&args.request,2*1024*1024)?;
            let (preparation,_)=server.hydrate_topo(&args.request.preparation,args.origin_workspace.as_ref(),&[])?;
            let native=collagene_core::topo_nick_sealing::Request{preparation,sealing:args.request.sealing.clone()};

            if args.request.sealing.nicks.len()>16 || args.request.sealing.ligase_assumption.len()>4096 {
                return Err(bad("Bounded explicit nick and ligase declarations required"));
            }
            let (inputs,_) = server.topo_sources(&args.sources,&native.preparation,None)?;
            let report=collagene_core::topo_nick_sealing::inspect(&inputs,&native).map_err(bad)?;
            report.validate_binding(&collagene_core::topo_nick_sealing::revision(&inputs,&native).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"origin_workspace":args.origin_workspace,"report":report}),32*1024*1024)?,"application/json","topo-nick-sealing.json")?;
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
            let (declaration,template)=server.hydrate_topo_declaration(&args.declaration,args.origin_workspace.as_ref(),args.template.as_ref())?;
            let compiled=topo_preparation::compile(&declaration).map_err(bad)?;
            {let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut refs=vec![];let mut sources=vec![];
                for source in &args.sources {if source.role.trim().is_empty()||source.role.len()>256||!seen.insert(source.role.as_str()){return Err(bad("unique bounded preparation source roles required"));}refs.push(store.document(&source.document)?);sources.push(verification::TopoSource{role:source.role.clone(),source:verification::Source::Input(source.role.clone())});}
                collagene_core::topo_run_preparation::source_inputs(&sources,&refs,&compiled.request,template.as_ref()).map_err(bad)?;
            }
            let artifact=server.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"declaration":args.declaration,"origin_workspace":args.origin_workspace,"template":args.template,"compiled":compiled}),32*1024*1024)?,"application/json","topo-preparation-declaration.json")?;
            let compact=json!({"request":topo_graph_tools::Preparation::from_native(&compiled.request),"vector_strands":compiled.vector_strands,"insert_strands":compiled.insert_strands,"directional":compiled.directional,"assumptions":compiled.assumptions});
            output(json!({"state":"declared-preparation","verified":false,"sources":args.sources,"origin_workspace":args.origin_workspace,"template":args.template,"compiled":compact,"artifact_id":artifact,"authority":"Current source-bound compiler instructions only. Saved graph/template records are owned native declarations; full chemistry is freshly replayed by inspection. Resume reuses existing insert calls and does not synthesize another insert or activation."}))
        }).await
    }
}
