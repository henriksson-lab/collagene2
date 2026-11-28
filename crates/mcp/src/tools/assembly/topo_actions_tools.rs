use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::references::WorkspaceRef;
use crate::tools::assembly::TopoSourceRef;
use crate::tools::assembly::topo_graph_tools;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;
use serde_json::json;

use collagene_core::cloning::topo::{actions, physical};
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProductQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::ActionProduct,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
impl CollageneServer {
    pub(crate) async fn project_actions(
        &self,
        args: ProductQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            bounded_json_limit(&args.request,2*1024*1024)?;
            let (preparation,_)=server.hydrate_topo(&args.request.preparation,args.origin_workspace.as_ref(),&[])?;
            let native=actions::Request{preparation,selection:args.request.selection.clone(),options:args.request.options.clone()};

            actions::request_preflight(&native).map_err(bad)?;
            let (inputs,_)=server.topo_sources(&args.sources,&native.preparation,None)?;
            let report=actions::inspect(&inputs,&native).map_err(bad)?;
            report.validate_binding(&actions::revision(&inputs,&native).map_err(bad)?).map_err(bad)?;
            let summary=report.duplex.as_ref().map(|d|d.summary()).transpose().map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{sources:&'a [TopoSourceRef],origin_workspace:&'a Option<WorkspaceRef>,report:&'a actions::Report}
            let (artifact,document)=server.session.candidate_document_artifact_bounded(&Evidence{sources:&args.sources,origin_workspace:&args.origin_workspace,report:&report},report.product.clone(),"topo-action-product.json",32*1024*1024)?;
            output(json!({"state":"source-reconstructed-action-contour","verified":false,"context_revision":report.context_revision,
                "contour_proven":matches!(report.projection.outcome,physical::projection::Outcome::Proven(_)),
                "withheld":match &report.projection.outcome{physical::projection::Outcome::Withheld(r)=>Some(r),_=>None},
                "duplex":summary,"component_count":report.projection.physical.product.components.len(),"source_calls":report.projection.physical.source_bases,
                "bridge_notes":report.bridge_notes,"document":document,"artifact_id":artifact,
                "authority":"Fresh chronological source graph reconstruction and explicit component/source-seed proof. Excised fragments remain full evidence. No isolation, experimental success or imported-document material authority is established. Use literal TopoActions plan declarations and independent target verification for downstream operations."}))
        }).await
    }
}
