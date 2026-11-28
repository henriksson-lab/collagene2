use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::references::WorkspaceRef;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

use collagene_core::cloning::assembly::library_history as history;
fn page_size() -> usize {
    20
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListQuery {
    pub workspace: WorkspaceRef,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "page_size")]
    pub limit: usize,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OpenQuery {
    pub workspace: WorkspaceRef,
    pub attachment_hash: String,
}
impl CollageneServer {
    pub(crate) async fn assembly_library_history_list(
        &self,
        args: ListQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 64 {
            return Err(bad("library history pages contain1..64 records"));
        }
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;
            let rows = history::list(&frame).map_err(bad)?;
            if args.offset > rows.len() { return Err(bad("library history page starts beyond the retained list")); }
            let total = rows.len();
            let page = rows.into_iter().skip(args.offset).take(args.limit).collect::<Vec<_>>();
            output(json!({"workspace":args.workspace,"records":page,"total":total,"offset":args.offset,
                "next_offset":args.offset.saturating_add(args.limit).min(total),"state":"historical-edit-recovery","verified":false,
                "authority":"Immutable raw edit recovery only. Declared retained source pins do not establish complete host draft validation or scientific authority. No draft, source head or design is restored."}))
        }).await
    }
    pub(crate) async fn assembly_library_history_open(
        &self,
        args: OpenQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;
            let record = history::load(&frame, &args.attachment_hash).map_err(bad)?;
            // Do not serialize before_draft as integer JSON: publish its exact
            // UTF-8 bytes separately, together with the complete typed metadata.
            let metadata = json!({"attachment_hash":args.attachment_hash,"schema_version":record.schema_version,
                "before_active_hash":record.before_active_hash,"before_hash":record.before_hash,"before_bytes":record.before_draft.len(),
                "after_kind":record.after_kind,"after_hash":record.after_hash,"library_id":record.library_id,
                "before_next_id":record.before_next_id,"after_next_id":record.after_next_id,"allocation_highwater":record.allocation_highwater,
                "proposal_revision":record.proposal_revision,"source_workspace_fingerprint":record.source_workspace_fingerprint,
                "sources":record.sources,"impacts":record.impacts,"state":"historical-edit-recovery","verified":false,
                "authority":"Exact historical recovery bytes and retained declared pins only. Opening does not apply a draft, rewind allocation identities, validate all host draft semantics or restore scientific authority."});
            let summary = json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"library_id":record.library_id,
                "before_hash":record.before_hash,"after_hash":record.after_hash,"allocation_highwater":record.allocation_highwater,
                "source_count":record.sources.len(),"impact_count":record.impacts.len(),"before_bytes":record.before_draft.len(),
                "state":"historical-edit-recovery","verified":false});
            let raw = String::from_utf8(record.before_draft).map_err(bad)?;
            let (metadata_artifact_id,before_draft_artifact_id) = server.session.paired_text_artifacts([
                (bounded_json_limit(&metadata,4*1024*1024)?,"application/json","assembly-library-history.json"),
                (raw,"application/json","assembly-library-before-draft.json"),
            ])?;
            let mut response = summary;
            response["metadata_artifact_id"] = json!(metadata_artifact_id);
            response["before_draft_artifact_id"] = json!(before_draft_artifact_id);
            output(response)
        }).await
    }
}
