pub mod quantity_tools;
pub mod unit_tools;

use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::Artifact;
use crate::session::ArtifactData;
use crate::session::MAX_ITEMS;
use crate::session::references::PlanRef;
use crate::session::references::WorkspaceRef;

use collagene_core::documents::workspace;
use collagene_core::formats::oligos as oligo_formats;
use collagene_core::verification;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use rmcp::tool;
use rmcp::tool_router;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProtocolLayoutQuery {
    pub workspace: WorkspaceRef,
    pub task: ProtocolLayoutTask,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum ProtocolLayoutTask {
    Preview {
        request: collagene_core::output::layout::Request,
    },
    Capture {
        request: collagene_core::output::layout::Request,
        expected_context_revision: String,
    },
    List {
        offset: usize,
        limit: usize,
    },
    Open {
        layout_hash: String,
    },
    Export {
        layout_hash: String,
        format: collagene_core::output::layout::ExportFormat,
    },
    ReviewPublication {
        layout_hash: String,
    },
    ProjectPublication {
        layout_hash: String,
        expected_scope_revision: String,
    },
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum PlanOligoFormat {
    Csv,
    Tsv,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlanOligoExport {
    pub plan_revision: String,
    pub format: PlanOligoFormat,
    #[serde(default)]
    pub options: collagene_core::planning::oligos::Options,
}
#[tool_router(router = protocol_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = r#"Plan explicit source-bound cloning reaction quantities using exact retained candidate/saved-run plan, independent target and ordered input event pins. Preview returns exact arithmetic plus missing-field diagnostics, not chemistry, stock availability, recommended recipes or experimental readiness. Capture requires reviewed context and atomically returns a NEW immutable workspace, full quantity artifact and archive without changing biological heads/status. List/Open inspect saved plans. Export formats Csv/Html/Pdf are quantity companions; ProtocolHtml/ProtocolPdf combine complete cloning instructions, fresh native molecular consistency and declared quantity tables. Molecular status and arithmetic completeness remain separate; neither establishes experimental success. Ordinary plasmid sharing omits this sidecar.

Exact quantity JSON: Volume={"unit":"nL","value":{"numerator":"20000","denominator":"1"}} is20uL. Rational numerator/denominator are decimal unsigned64-bit STRINGS; denominator positive. Never use float, exponent or NaN. Canonical mass concentration unit is"ng/nL":100ng/uL=value1/10. Canonical molar concentration unit is"pmol/nL":10uM=value1/100. Never interchange dimensional tags. Blank direct volume is{"Direct":{"volume":null}}, not zero. MassDilution has {"recipe":{"schema_version":1,"stock":<mass concentration or null>,"target":<mass concentration or null>,"final_volume":<Volume or null>},"quantization":null}; MolarDilution uses molar concentrations. Explicit quantization={"increment":<Volume>,"mode":"ExactOnly"|"Floor"|"Ceiling"|"NearestTiesUp"}; inspect achieved concentration/deviation. Nullable request fields must be present.

MINIMAL_PREVIEW={"workspace":{"id":"<workspace id>","revision":"<workspace revision>"},"task":{"Preview":{"request":{"name":"Explicit PCR quantities","plan_source":{"SavedRun":{"run_id":"<retained run id>"}},"plan_revision":"<plan revision>","target":{"document_id":"<target document>","event_id":"<original target event>","revision":"<target revision>"},"inputs":[{"input_id":"template","event":{"document_id":"<source document>","event_id":"<original source event>","revision":"<source revision>"}}],"reactions":[{"id":"pcr","step_index":0,"final_volume":{"unit":"nL","value":{"numerator":"20000","denominator":"1"}},"fill":{"ToFinal":{"component_id":"water","label":"Operator-declared water"}},"components":[{"id":"template-volume","label":"User-declared template volume","source":{"Input":"template"},"volume":{"Direct":{"volume":{"unit":"nL","value":{"numerator":"2000","denominator":"1"}}}}}]}]}}}} END_PREVIEW

Replace ALL angle-bracket pins with exact native values, never names or equal-DNA substitutions. This arithmetic example leaves all other PCR reagents unspecified;18uL fill is an operator declaration, NOT a recommended reaction recipe. Arithmetic completeness is not scientific readiness. Use plan_source="Candidate" only for an exactly retained candidate and explicitly declared events. Capture adds expected_context_revision from Preview. Open/Export use quantity_hash from capture; Export format Csv,Html,Pdf,ProtocolHtml,ProtocolPdf. Combined outputs include a fresh current-engine molecular receipt, actual findings/assumptions in the instructions and exact original target/input events. They never reuse historical scientific status as fresh acceptance. List requires offset and limit1..200."#
    )]
    pub(crate) async fn protocol_quantities(
        &self,
        Parameters(args): Parameters<quantity_tools::Query>,
    ) -> Result<CallToolResult, ErrorData> {
        self.reaction_quantities(args).await
    }
    #[tool(
        description = "Plan explicit96/384-well documentary cloning assignments using exact retained batch declarations or released run/product events. Preview never records; Capture requires the reviewed context and creates a NEW immutable workspace with a bounded layout attachment, preserving all biological heads and historical status. List/Open inspect saved layouts. Export emits native CSV/HTML/PDF planning documents with exact identifiers and missing-quantity notices. No automatic placement, inventory, schedules, instrument commands, fabricated amounts or experimental success. Ordinary plasmid publication omits plate layouts; complete plate publication requires a separate explicit closure review."
    )]
    pub(crate) async fn protocol_layout(
        &self,
        Parameters(args): Parameters<ProtocolLayoutQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::output::layout;
            let frame=server.session.restore_workspace(&args.workspace)?;
            match args.task {
                ProtocolLayoutTask::Preview{request}=>{
                    let expected=layout::preview_revision(&frame,&request).map_err(bad)?;
                    let report=layout::preview(&frame,&request).map_err(bad)?;
                    report.validate_binding(&expected).map_err(bad)?;
                    let text=bounded_json_limit(&report,layout::MAX_BYTES)?;
                    let artifact=server.session.text_artifact(text,"application/json","protocol-layout.json")?;
                    output(json!({"context_revision":report.context_revision,"assignment_count":report.resolved.len(),"warnings":report.warnings,"artifact_id":artifact,"state":"documentary preview","verified":false,"authority":"Well assignments do not verify reactions or laboratory execution. Quantities remain missing unless explicitly provided by a separate typed planning workflow."}))
                }
                ProtocolLayoutTask::Capture{request,expected_context_revision}=>{
                    let expected=layout::preview_revision(&frame,&request).map_err(bad)?;
                    if expected!=expected_context_revision{return Err(bad("Layout source or assignments changed; preview again before capture."));}
                    let commit=layout::capture(&frame,&request).map_err(bad)?;
                    commit.validate_binding(&expected).map_err(bad)?;
                    let text=bounded_json_limit(&commit.layout,layout::MAX_BYTES)?;
                    let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,"source_workspace_fingerprint":commit.source_workspace_fingerprint,"assignment_count":commit.layout.resolved.len(),"state":"documentary","verified":false});
                    let bytes=workspace::save_workspace(&commit.workspace).map_err(bad)?;
                    server.session.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(bytes),"protocol-layout.json","Documentary well assignments only; no laboratory execution or scientific verification authority.")
                }
                ProtocolLayoutTask::List{offset,limit}=>{
                    if limit==0||limit>200||offset>layout::MAX_LAYOUTS{return Err(bad("Invalid layout page bounds."));}
                    let mut entries=Vec::new();
                    for (hash,asset) in frame.attachments.iter().filter(|(_,a)|a.kind==layout::KIND){
                        if entries.len()>=layout::MAX_LAYOUTS{return Err(bad("Too many layout attachments."));}
                        let report=layout::load(&frame,hash).map_err(bad)?;
                        entries.push(json!({"layout_hash":hash,"name":asset.name,"context_revision":report.context_revision,"plate":report.request.plate,"assignment_count":report.request.assignments.len()}));
                    }
                    let total=entries.len();let page=entries.into_iter().skip(offset).take(limit).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"offset":offset,"total":total,"next_offset":(offset+page.len()<total).then_some(offset+page.len()),"layouts":page,"authority":"Historical documentary layouts only."}))
                }
                ProtocolLayoutTask::Open{layout_hash}=>{
                    let report=layout::load(&frame,&layout_hash).map_err(bad)?;
                    let artifact=server.session.text_artifact(bounded_json_limit(&report,layout::MAX_BYTES)?,"application/json","protocol-layout.json")?;
                    output(json!({"layout_hash":layout_hash,"context_revision":report.context_revision,"assignment_count":report.resolved.len(),"artifact_id":artifact,"state":"documentary","verified":false}))
                }
                ProtocolLayoutTask::ReviewPublication{layout_hash}=>{
                    let review=collagene_library_sharing::inspect_protocol_layout_closure(&frame,&layout_hash).map_err(bad)?;
                    let text=bounded_json_limit(&review,32*1024*1024)?;
                    let scope_revision=blake3::hash(text.as_bytes()).to_hex().to_string();
                    let artifact=server.session.text_artifact(text,"application/json","protocol-layout-publication.json")?;
                    output(json!({"layout_hash":layout_hash,"scope_revision":scope_revision,"documents":review.inventory.documents.len(),"events":review.inventory.events.len(),"attachments":review.inventory.attachments.len(),"omitted":review.inventory.omitted,"notices":review.notices,"artifact_id":artifact,"authority":"Explicit layout scope review; no archive or publication has been made."}))
                }
                ProtocolLayoutTask::ProjectPublication{layout_hash,expected_scope_revision}=>{
                    let review=collagene_library_sharing::inspect_protocol_layout_closure(&frame,&layout_hash).map_err(bad)?;
                    let text=bounded_json_limit(&review,32*1024*1024)?;
                    if blake3::hash(text.as_bytes()).to_hex().as_str()!=expected_scope_revision{return Err(bad("Layout publication scope changed; review the full dependency counts again."));}
                    let projected=collagene_library_sharing::project_protocol_layout_closure(&frame,&layout_hash).map_err(bad)?;
                    let bytes=workspace::save_workspace(&projected.workspace).map_err(bad)?;
                    server.session.commit_workspace_outputs_named(projected.workspace,Some((text,json!({"layout_hash":layout_hash,"scope_revision":expected_scope_revision,"authority":"Explicit complete layout projection, including indivisible batch/release dependencies. No server upload or experimental authority."}))),Some(bytes),"protocol-layout-publication.json","Reviewed explicit layout dependency projection only; no server upload, experiment or scientific verification.")
                }
                ProtocolLayoutTask::Export{layout_hash,format}=>{
                    let expected=layout::export_revision(&frame,&layout_hash,format).map_err(bad)?;
                    let exported=layout::export(&frame,&layout_hash,format).map_err(bad)?;
                    exported.validate_binding(&expected).map_err(bad)?;
                    let (mime,filename)=match format{layout::ExportFormat::Csv=>("text/csv","protocol-layout.csv"),layout::ExportFormat::Html=>("text/html","protocol-layout.html"),layout::ExportFormat::Pdf=>("application/pdf","protocol-layout.pdf")};
                    let length=exported.bytes.len();
                    let artifact=server.session.binary_artifact(exported.bytes,mime,filename)?;
                    output(json!({"layout_hash":layout_hash,"context_revision":expected,"artifact_id":artifact,"filename":filename,"mime_type":mime,"encoding":"base64","bytes":length,"warnings":exported.warnings,"authority":"Documentary planning export, not experimental execution evidence."}))
                }
            }
        }).await
    }
    #[tool(
        description = "Calculate exact single-stock dilution from decimal-string quantities and explicit mass or molar units. Missing values and mixed dimensions are errors. Optional explicit stock-volume quantization reports absolute rounding deviation, refill diluent and achieved concentration; omitted quantization never rounds. Assumes C1*V1=C2*V2, additive volumes and zero-concentration diluent. No stock inventory lookup, solvent compatibility or reaction verification."
    )]
    pub(crate) fn calculate_dilution(
        &self,
        Parameters(args): Parameters<unit_tools::DilutionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(unit_tools::calculate_dilution(args)?)
    }
    #[tool(
        description = "Balance up to 128 uniquely named declared component volumes against an explicit final volume, returning exact specified and remaining fill-diluent volumes. Uses decimal strings and nL/uL/mL unit enums; missing component volumes, overfill and bounded arithmetic overflow remain errors. Additive volume arithmetic only, without concentration, inventory or chemical-compatibility inference."
    )]
    pub(crate) fn balance_mixture(
        &self,
        Parameters(args): Parameters<unit_tools::MixtureQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(unit_tools::balance_mixture(args)?)
    }
    #[tool(
        description = "Calculate duplex-DNA molecule-ratio masses, amounts and optional stock volumes for one weighed anchor and up to 127 explicitly named targets. Requires a caller-declared positive g/mol/base-pair mass coefficient, model identity and assumptions; no implicit 660 coefficient or recommended ratio. Exact decimal-string inputs and rational-string outputs, with explicit missing-stock/invalid/overflow errors. Amounts count molecules, not ends; no sequence chemistry, stock inventory or reaction success is verified."
    )]
    pub(crate) fn calculate_dna_molar_ratios(
        &self,
        Parameters(args): Parameters<unit_tools::DnaRatioQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(unit_tools::calculate_dna_molar_ratios(args)?)
    }
    #[tool(
        description = "Generate an HTML protocol artifact from an exact saved plan revision. Always re-verifies at export, visibly includes invalid/incomplete findings and assumptions. Read the artifact in bounded chunks; no files are written."
    )]
    pub(crate) async fn export_protocol(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.export_protocol_inner(args))
            .await
    }
    #[tool(
        description = "Export actual declared oligos from an exact saved plan and its stored independent target as spreadsheet-safe CSV or TSV plus complete use/domain/tail/mismatch evidence JSON. Exact ordered DNA/modifications may deduplicate; every original use and name remains separate. No PCR/target verification is run, no primer is designed, and downstream phosphorylation never becomes a supplier order modification. Missing/unrepresentable original order data is explicit complete:false/unavailable; do not order from a partial sheet. No quantities, stocks or supplier codes are invented."
    )]
    pub(crate) async fn export_plan_oligos(
        &self,
        Parameters(args): Parameters<PlanOligoExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let delimiter = match args.format {
            PlanOligoFormat::Csv => oligo_formats::Delimiter::Comma,
            PlanOligoFormat::Tsv => oligo_formats::Delimiter::Tab,
        };
        let filename = match args.format {
            PlanOligoFormat::Csv => "plan-oligos.csv",
            PlanOligoFormat::Tsv => "plan-oligos.tsv",
        };
        let server = self.clone();
        self.compute(move || {
            let exported={let store=server.session.store.lock().unwrap();let saved=store.plans.get(&args.plan_revision).ok_or_else(||bad("unknown plan revision"))?;
                let expected=collagene_core::planning::oligos::export_revision(&saved.plan,&saved.target,&args.options,delimiter).map_err(bad)?;
                let exported=collagene_core::planning::oligos::export(&saved.plan,&saved.target,&args.options,delimiter).map_err(bad)?;
                exported.validate_binding(&expected).map_err(bad)?;exported};
            if exported.report.plan_revision!=args.plan_revision{return Err(bad("stored plan revision changed"));}
            let result=json!({"plan_revision":exported.report.plan_revision,"target_revision":exported.report.target_revision,
                "context_revision":exported.context_revision,"options":exported.report.options,"verified":false,
                "complete":exported.report.complete,"unavailable":exported.report.unavailable,
                "oligo_count":exported.report.oligos.records.len(),"use_count":exported.report.uses.len(),
                "uses":exported.report.uses.iter().take(100).collect::<Vec<_>>(),"uses_truncated":exported.report.uses.len()>100,
                "warnings":exported.report.warnings,"export_warnings":exported.export.warnings,
                "authority":"Declared order/use evidence only. PCR products, independent target agreement and experimental success require separate verification. Incomplete sheets are review evidence, not all-required orders; declared starting termini are not supplier modification requests."});
            let report_text=bounded_json_limit(&exported.report,32*1024*1024)?;
            let (sheet,report)=server.session.paired_text_artifacts([(exported.export.text,"text/plain",filename),(report_text,"application/json","plan-oligos.json")])?;
            let mut result=result;result["table_artifact_id"]=json!(sheet);result["report_artifact_id"]=json!(report);result["filename"]=json!(filename);output(result)
        }).await
    }
    #[tool(
        description = "Generate an embedded vector PDF protocol from an exact saved plan revision, re-verifying all steps and the independent target at export. Invalid/incomplete drafts remain conspicuous. Returns a bounded base64 artifact; read chunks, concatenate and decode to a PDF. No files, browser or external renderer are invoked."
    )]
    pub(crate) async fn export_protocol_pdf(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let saved=server.session.store.lock().unwrap().plans.get(&args.plan_revision).cloned().ok_or_else(||bad("unknown plan revision"))?;
            let report=collagene_core::output::pdf::export_protocol_pdf(&saved.plan,&saved.target).map_err(bad)?;
            let id=server.session.binary_artifact(report.bytes.clone(),"application/pdf","cloning-protocol.pdf")?;
            output(json!({"artifact_id":id,"filename":"cloning-protocol.pdf","mime_type":"application/pdf","encoding":"base64","bytes":report.bytes.len(),"pages":report.pages,"status":report.status,"plan_revision":report.plan_revision,"target_revision":report.target_revision,"engine_version":report.engine_version,"warnings":report.warnings}))
        }).await
    }
    pub(crate) fn export_protocol_inner(&self, args: PlanRef) -> Result<CallToolResult, ErrorData> {
        let saved = self
            .session
            .store
            .lock()
            .unwrap()
            .plans
            .get(&args.plan_revision)
            .cloned()
            .ok_or_else(|| bad("unknown plan revision"))?;
        let html = verification::output::export_protocol_html(&saved.plan, &saved.target);
        let id = blake3::hash(html.as_bytes()).to_hex().to_string();
        let length = html.len();
        let mut store = self.session.store.lock().unwrap();
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(length)?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Text(html),
                    mime_type: "text/html; charset=utf-8",
                    encoding: "utf-8",
                    filename: "cloning-protocol.html",
                },
            );
        }
        output(
            json!({"artifact_id":id,"mime_type":"text/html; charset=utf-8","bytes":length,"plan_revision":args.plan_revision}),
        )
    }
}
