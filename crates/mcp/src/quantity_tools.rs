//! Immutable source-bound documentary arithmetic, separate from science status.
use super::*;
use collagene_core::protocol_quantities as quantities;
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Query {
    pub workspace: WorkspaceRef,
    pub task: Task,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum Task {
    Preview {
        request: quantities::Request,
    },
    Capture {
        request: quantities::Request,
        expected_context_revision: String,
    },
    List {
        offset: usize,
        limit: usize,
    },
    Open {
        quantity_hash: String,
    },
    Export {
        quantity_hash: String,
        format: quantities::ExportFormat,
    },
}
impl CollageneServer {
    pub(super) async fn reaction_quantities(
        &self,
        args: Query,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;
            match args.task {
                Task::Preview{request}=>{
                    let expected=quantities::preview_revision(&frame,&request).map_err(bad)?;
                    let preview=quantities::preview(&frame,&request).map_err(bad)?;
                    preview.validate_binding(&expected).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&preview,quantities::MAX_BYTES)?,"application/json","reaction-quantities.json")?;
                    output(json!({"state":"documentary quantity preview","verified":false,"context_revision":expected,"arithmetic_complete":preview.complete,"reaction_count":preview.reactions.len(),"finding_count":preview.findings.len(),"findings":preview.findings.iter().take(40).collect::<Vec<_>>(),"findings_truncated":preview.findings.len()>40,"warnings":preview.warnings,"artifact_id":artifact,"authority":"Explicit user quantities and exact arithmetic only. Completeness is not scientific verification, reaction readiness, stock availability, isolation or experimental success."}))
                }
                Task::Capture{request,expected_context_revision}=>{
                    let expected=quantities::preview_revision(&frame,&request).map_err(bad)?;
                    if expected!=expected_context_revision{return Err(bad("Quantity source, operation or declaration changed; preview again before capture."));}
                    let commit=quantities::capture(&frame,&request).map_err(bad)?;
                    commit.validate_binding(&expected).map_err(bad)?;
                    let text=bounded_json_limit(&commit.preview,quantities::MAX_BYTES)?;
                    let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,"source_workspace_fingerprint":commit.source_workspace_fingerprint,"arithmetic_complete":commit.preview.complete,"reaction_count":commit.preview.reactions.len(),"state":"documentary","verified":false});
                    let archive=workspace::save_workspace(&commit.workspace).map_err(bad)?;
                    server.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(archive),"reaction-quantities.json","Source-bound documentary arithmetic only; unchanged biological heads and historical assessments. No reagent stock, experimental execution or scientific verification authority.")
                }
                Task::List{offset,limit}=>{
                    if limit==0||limit>200||offset>quantities::MAX_RECORDS{return Err(bad("Invalid quantity-document page bounds."));}
                    let expected=quantities::list_revision(&frame).map_err(bad)?;
                    let listed=quantities::list(&frame).map_err(bad)?;
                    listed.validate_binding(&expected).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&listed,quantities::MAX_BYTES)?,"application/json","reaction-quantity-catalog.json")?;
                    let total=listed.entries.len();
                    let page=listed.entries.iter().skip(offset).take(limit).map(|s|json!({"quantity_hash":s.attachment_hash,"name":s.name,"plan_revision":s.plan_revision,"arithmetic_complete":s.complete,"source_available":s.available,"finding_count":s.findings.len(),"findings":s.findings.iter().take(10).collect::<Vec<_>>(),"findings_truncated":s.findings.len()>10,"reaction_count":s.reaction_count})).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"context_revision":expected,"offset":offset,"total":total,"next_offset":(offset+page.len()<total).then_some(offset+page.len()),"quantities":page,"artifact_id":artifact,"authority":"Retained documentary quantity plans, not scientific assessments or stocks. Exact source/target event pins are retained in the full bounded catalog artifact."}))
                }
                Task::Open{quantity_hash}=>{
                    let expected=quantities::open_revision(&frame,&quantity_hash).map_err(bad)?;
                    let opened=quantities::open(&frame,&quantity_hash).map_err(bad)?;
                    opened.validate_binding(&expected).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&opened.preview,quantities::MAX_BYTES)?,"application/json","reaction-quantities.json")?;
                    output(json!({"quantity_hash":quantity_hash,"context_revision":expected,"historical_context_revision":opened.preview.context_revision,"arithmetic_complete":opened.preview.complete,"reaction_count":opened.preview.reactions.len(),"finding_count":opened.preview.findings.len(),"findings":opened.preview.findings.iter().take(40).collect::<Vec<_>>(),"findings_truncated":opened.preview.findings.len()>40,"warnings":opened.preview.warnings,"artifact_id":artifact,"state":"documentary","verified":false,"authority":"Exact retained native source events and arithmetic were checked; historical or current reaction success is not asserted."}))
                }
                Task::Export{quantity_hash,format}=>{
                    let expected=quantities::export_revision(&frame,&quantity_hash,format).map_err(bad)?;
                    let exported=quantities::export(&frame,&quantity_hash,format).map_err(bad)?;
                    exported.validate_binding(&expected).map_err(bad)?;
                    let length=exported.bytes.len();
                    let (mime,filename)=match format{quantities::ExportFormat::Csv=>("text/csv","reaction-quantities.csv"),quantities::ExportFormat::Html=>("text/html","reaction-quantities.html"),quantities::ExportFormat::Pdf=>("application/pdf","reaction-quantities.pdf"),quantities::ExportFormat::ProtocolHtml=>("text/html","cloning-protocol-quantities.html"),quantities::ExportFormat::ProtocolPdf=>("application/pdf","cloning-protocol-quantities.pdf")};
                    if exported.media_type!=mime||exported.filename!=filename{return Err(bad("Quantity export representation differs from requested format."));}
                    let artifact=server.binary_artifact(exported.bytes,mime,filename)?;
                    let combined=matches!(format,quantities::ExportFormat::ProtocolHtml|quantities::ExportFormat::ProtocolPdf);
                    output(json!({"quantity_hash":quantity_hash,"context_revision":expected,"source_workspace_fingerprint":exported.source_workspace_fingerprint,"plan_revision":exported.plan_revision,"target":exported.target,"inputs":exported.inputs,"arithmetic_complete":exported.arithmetic_complete,"molecular":exported.molecular,"artifact_id":artifact,"filename":exported.filename,"mime_type":exported.media_type,"encoding":"base64","bytes":length,"content_hash":exported.content_hash,"warnings":exported.warnings,"state":if combined{"freshly checked cloning protocol with declared quantities"}else{"documentary quantity companion"},"verified":false,"authority":if combined{"Molecular status is a fresh native consistency result for this exact retained plan/target, separate from declared arithmetic completeness. Neither establishes stock availability, reaction readiness, isolation or experimental success."}else{"Quantity arithmetic companion only; no molecular verification, recommended reaction conditions or laboratory execution is asserted."}}))
                }
            }
        }).await
    }
}
