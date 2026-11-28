//! Native-history CDS evidence over an exact immutable session workspace.
use crate::{
    server::{
        CollageneServer,
        response::{bad, bounded_json_limit, output},
    },
    session::references::WorkspaceRef,
};
use collagene_core::construct_qc::causal;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;
fn page_size() -> usize {
    1
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CausalQcQuery {
    pub workspace: WorkspaceRef,
    pub request: causal::Request,
    #[serde(default)]
    pub pair_start: usize,
    #[serde(default = "page_size")]
    pub max_pairs: usize,
}
#[tool_router(router=causal_qc_router,vis="pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Compare explicitly paired source/product CDS annotations through freshly replayed retained native history in one exact immutable workspace. Supply exact event bindings, current annotation assumptions and declared genetic code. Returns per-pair typed native-history, annotation/coding correspondence, junction-register and definite internal-stop evidence with explicit completeness, native transition pins and bounded previews; complete bounded JSON is retained as an immutable artifact. Unretained analysis previews or missing/unsupported correspondence remain indeterminate, never a guessed sequence match. No overall biological pass, expression/function prediction, material authority, source edits, history writes or implicit database access."
    )]
    pub(crate) async fn assess_construct_causal_qc(
        &self,
        Parameters(args): Parameters<CausalQcQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=4).contains(&args.max_pairs) {
            return Err(bad("causal QC pair page size must be 1..4"));
        }
        bounded_json_limit(&args.request, causal::MAX_REQUEST_BYTES)?;
        let server = self.clone();
        self.compute(move ||{
            let workspace=server.session.restore_workspace(&args.workspace)?;
            args.request.validate(&workspace).map_err(bad)?;
            let context=causal::revision(&workspace,&args.request).map_err(bad)?;
            let report=causal::assess(&workspace,&args.request).map_err(bad)?;
            report.validate_binding(&context).map_err(bad)?;
            report.validate_request(&workspace,&args.request).map_err(bad)?;
            if args.pair_start>report.pairs.len(){return Err(bad("pair offset exceeds retained causal QC pairs"));}
            let end=args.pair_start.saturating_add(args.max_pairs).min(report.pairs.len());
            let pairs=report.pairs[args.pair_start..end].iter().map(|pair|json!({
                "id":pair.id,"source":pair.source,"product":pair.product,
                "source_cds":pair.source_cds,"product_cds":pair.product_cds,"code":pair.code,
                "annotation_correspondence":pair.annotation_correspondence,"mapping_complete":pair.mapping_complete,
                "finding_count":pair.findings.len(),"findings":pair.findings.iter().take(8).collect::<Vec<_>>(),
                "findings_preview_truncated":pair.findings.len()>8,
                "register_segment_count":pair.register_segments.len(),"register_segments":pair.register_segments.iter().take(8).collect::<Vec<_>>(),
                "register_preview_truncated":pair.register_segments.len()>8,
                "internal_stop_count":pair.new_stops.len(),"new_stops":pair.new_stops.iter().take(8).collect::<Vec<_>>(),
                "stop_preview_truncated":pair.new_stops.len()>8,
            })).collect::<Vec<_>>();
            let mut summary=json!({"workspace":args.workspace,"engine_version":report.engine_version,
                "context_revision":report.context_revision,"workspace_fingerprint":report.workspace_fingerprint,
                "request_revision":report.request_revision,"report_revision":report.report_revision,
                "pair_count":report.pairs.len(),"pair_start":args.pair_start,"next_start":end,"done":end==report.pairs.len(),
                "page_truncated":args.pair_start!=0 || end!=report.pairs.len(),"pairs":pairs,
                "work":report.work,"coding_calls":report.coding_calls,"history_events":report.history_events,
                "statement":report.statement,"persistence":"session evidence only"});
            // Reject oversized summaries before publishing an artifact. Paging
            // never implies that omitted native evidence was absent or complete.
            bounded_json_limit(&summary,1024*1024-1024)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&report,causal::MAX_REPORT_BYTES)?,"application/json","construct-causal-qc.json")?;
            summary["artifact_id"]=json!(artifact);
            output(summary)
        }).await
    }
}
