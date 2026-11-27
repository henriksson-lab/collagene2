//! Fresh explicit graph actions; declarations never accept imported material.
use super::*;
use collagene_core::{topo::physical, topo_actions as actions};
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Actions,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProductQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::ActionProduct,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
impl CollageneServer {
    pub(super) async fn inspect_actions(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            bounded_json_limit(&args.request,2*1024*1024)?;
            let (preparation,_)=server.hydrate_topo(&args.request.preparation,args.origin_workspace.as_ref(),&[])?;
            let native=actions::ActionsRequest{preparation,processing:args.request.processing.clone()};

            physical::processing::request_preflight(&args.request.processing,&native.preparation.engine_limits).map_err(bad)?;
            let (inputs,_)=server.topo_sources(&args.sources,&native.preparation,None)?;
            let report=actions::inspect_actions(&inputs,&native).map_err(bad)?;
            report.validate_binding(&actions::revision_actions(&inputs,&native).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"origin_workspace":args.origin_workspace,"report":report}),32*1024*1024)?,"application/json","topo-action-inspection.json")?;
            output(json!({"state":"declared-action-graph","verified":false,"context_revision":report.context_revision,
                "processing_context_revision":report.processing.context_revision,"source_calls":report.processing.source_bases,
                "steps":report.processing.steps,"component_count":report.processing.product.components.len(),
                "components":report.processing.product.components.iter().enumerate().take(20).map(|(i,c)|json!({"component":i,"chains":c.chains,"paired_calls":c.paired_nucleotides,"unpaired_calls":c.unpaired_nucleotides,"covalently_closed":c.covalently_closed})).collect::<Vec<_>>(),
                "remaining_nicks":report.processing.product.nicks,"artifact_id":artifact,
                "authority":"Explicit ideal hydrolysis/sealing under declared assumptions, not a named nuclease, kit, host-repair or efficacy claim. Original and excised calls/components remain retained; no target chooses a cut or a product."}))
        }).await
    }
    pub(super) async fn project_actions(
        &self,
        args: ProductQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            bounded_json_limit(&args.request,2*1024*1024)?;
            let (preparation,_)=server.hydrate_topo(&args.request.preparation,args.origin_workspace.as_ref(),&[])?;
            let native=actions::Request{preparation,processing:args.request.processing.clone(),selection:args.request.selection.clone(),options:args.request.options.clone()};

            actions::request_preflight(&native).map_err(bad)?;
            let (inputs,_)=server.topo_sources(&args.sources,&native.preparation,None)?;
            let report=actions::inspect(&inputs,&native).map_err(bad)?;
            report.validate_binding(&actions::revision(&inputs,&native).map_err(bad)?).map_err(bad)?;
            let summary=report.duplex.as_ref().map(|d|d.summary()).transpose().map_err(bad)?;
            let (artifact,document)=server.candidate_document_artifact_bounded(&json!({"sources":args.sources,"origin_workspace":args.origin_workspace,"report":report}),report.product.clone(),"topo-action-product.json",32*1024*1024)?;
            output(json!({"state":"source-reconstructed-action-contour","verified":false,"context_revision":report.context_revision,
                "contour_proven":matches!(report.projection.outcome,physical::projection::Outcome::Proven(_)),
                "withheld":match &report.projection.outcome{physical::projection::Outcome::Withheld(r)=>Some(r),_=>None},
                "duplex":summary,"component_count":report.projection.processing.product.components.len(),"source_calls":report.projection.processing.source_bases,
                "bridge_notes":report.bridge_notes,"document":document,"artifact_id":artifact,
                "authority":"Fresh source graph/action reconstruction and explicit component/source-seed proof. Excised fragments remain full evidence. No isolation, experimental success or imported-document material authority is established. Use literal TopoActions plan declarations and independent target verification for downstream operations."}))
        }).await
    }
}
