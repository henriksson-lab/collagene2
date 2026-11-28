use crate::session::references::MaterialDerivationOrigin;

use crate::server::CollageneServer;
use crate::server::EmptyRequest;
use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::Artifact;
use crate::session::ArtifactData;
use crate::session::MAX_ITEMS;
use crate::session::MAX_UPLOAD_CHUNK;
use crate::session::MAX_UPLOADS;
use crate::session::WorkspaceUpload;
use crate::session::references::DocumentRef;
use crate::session::references::WorkspaceRef;
use crate::session::references::WorkspaceUploadRef;

use base64::engine::general_purpose::STANDARD;
use collagene_core::sequence::analysis::Primer;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;

use rmcp::tool;
use rmcp::tool_router;
use serde_json::json;

use base64::Engine as _;
use collagene_core::documents::transfer;
use collagene_core::documents::workspace;
use collagene_core::verification;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MapExport {
    pub document: DocumentRef,
    pub format: collagene_core::output::raster::Format,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportWorkspace {
    /// Base64 of a .cg2design snapshot exported by the web workspace; at most 8 MiB decoded.
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BeginWorkspaceUpload {
    /// Exact decoded native archive byte length, at most 64 MiB.
    pub expected_bytes: usize,
    /// Lower-case BLAKE3 digest of the complete decoded archive.
    pub blake3: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendWorkspaceUpload {
    pub upload_id: String,
    /// Exact next decoded byte offset. No overlap, gaps or implicit replay.
    pub offset: usize,
    /// Independently base64-encoded chunk, at most 256 KiB decoded.
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunMaterialOriginsQuery {
    pub workspace: WorkspaceRef,
    pub event: collagene_core::material::lineage::EventBinding,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialOriginQuery {
    pub workspace: WorkspaceRef,
    pub event: collagene_core::material::lineage::EventBinding,
    /// None discovers choices; Some explicitly replays that exact retained origin.
    pub selection: Option<collagene_core::material::run::MaterialSelection>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialPcrChoice {
    pub selected_index: usize,
    pub expected_preview_revision: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum MaterialDerivationOperation {
    AnnealedOligos {
        options: collagene_core::chemistry::AnnealOptions,
        selection: collagene_core::material::derivation::AnnealSelection,
    },
    TerminalProcessing {
        action: collagene_core::material::derivation::TerminalAction,
    },
    OverlapExtensionPcr {
        preparation: collagene_core::cloning::overlap::extension::pcr::Request,
        selection: collagene_core::material::derivation::OverlapExtensionSelection,
    },
    ExactPcr {
        forward: Primer,
        reverse: Primer,
        selection: MaterialPcrChoice,
        template_strand: collagene_core::chemistry::Strand,
    },
    MismatchPcr {
        forward: Primer,
        reverse: Primer,
        options: collagene_core::cloning::pcr::mutagenesis::PcrOptions,
        selection: MaterialPcrChoice,
    },
    EndRepair {
        method: collagene_core::chemistry::EndRepairMethod,
    },
    CatalogDigest {
        selection: collagene_core::cloning::restriction::catalog::CatalogSelection,
        options: collagene_core::cloning::restriction::inspection::InspectionOptions,
        fragment_index: usize,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialDerivationSource {
    pub role: collagene_core::material::derivation::SourceRole,
    pub event: collagene_core::material::lineage::EventBinding,
    pub origin: MaterialDerivationOrigin,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialDerivationInstallRequest {
    pub sources: Vec<MaterialDerivationSource>,
    pub derivation: MaterialDerivationOperation,
    pub name: String,
    pub save_primers: bool,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum MaterialDerivationTask {
    PreviewAnnealedOligos {
        top: collagene_core::material::lineage::EventBinding,
        bottom: collagene_core::material::lineage::EventBinding,
        options: collagene_core::chemistry::AnnealOptions,
    },
    PreviewTerminalProcessing {
        source: collagene_core::material::lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        action: collagene_core::material::derivation::TerminalAction,
    },
    PreviewOverlapExtensionPcr {
        source: collagene_core::material::lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        preparation: collagene_core::cloning::overlap::extension::pcr::Request,
    },
    PreviewPcr {
        source: collagene_core::material::lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        forward: Primer,
        reverse: Primer,
        template_strand: collagene_core::chemistry::Strand,
    },
    PreviewMismatchPcr {
        source: collagene_core::material::lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        forward: Primer,
        reverse: Primer,
        options: collagene_core::cloning::pcr::mutagenesis::PcrOptions,
    },
    PreviewRepair {
        source: collagene_core::material::lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        method: collagene_core::chemistry::EndRepairMethod,
    },
    PreviewDigest {
        source: collagene_core::material::lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        selection: collagene_core::cloning::restriction::catalog::CatalogSelection,
        options: collagene_core::cloning::restriction::inspection::InspectionOptions,
    },
    Install {
        request: MaterialDerivationInstallRequest,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialDerivationQuery {
    pub workspace: WorkspaceRef,
    pub task: MaterialDerivationTask,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunMaterialPreparationQuery {
    pub workspace: WorkspaceRef,
    pub selection: collagene_core::material::run::Selection,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecordWorkspaceRun {
    pub workspace: WorkspaceRef,
    pub plan_revision: String,
    /// Explicit target event/document/revision binding; never inferred by sequence hash.
    pub target_binding: collagene_core::material::lineage::EventBinding,
    /// Explicit plan input IDs and their exact historical workspace events.
    pub input_bindings: Vec<collagene_core::material::lineage::InputEventBinding>,
    #[serde(default)]
    pub allow_incomplete: bool,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FinalizeCloningOutputs {
    pub workspace: WorkspaceRef,
    pub plan_revision: String,
    /// Exact independent target event, never selected by DNA equality.
    pub target_binding: collagene_core::material::lineage::EventBinding,
    /// Exact original event binding for each plan-local input ID.
    pub input_bindings: Vec<collagene_core::material::lineage::InputEventBinding>,
    #[serde(default)]
    pub options: collagene_core::material::outputs::Options,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReplayWorkspaceRun {
    pub workspace: WorkspaceRef,
    pub run_id: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportArchive {
    /// Base64 native .cg2workspace bytes; <=8 MiB decoded, else use upload tools.
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadReviewQuery {
    pub workspace: WorkspaceRef,
    pub request: collagene_core::review::reads::review::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecordReadCheck {
    pub workspace: WorkspaceRef,
    pub reference: collagene_core::material::lineage::EventBinding,
    pub reads: Vec<collagene_core::material::lineage::EventBinding>,
    pub options: collagene_core::review::reads::trace_validation::TraceValidationOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListReadChecks {
    pub workspace: WorkspaceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReplayReadCheck {
    pub workspace: WorkspaceRef,
    pub attachment_hash: String,
}
#[tool_router(router = workspace_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Import a full native .cg2workspace archive from base64 (<=8 MiB decoded convenience). Histories, candidates, source attachments and exact revisions are retained atomically. Candidates remain unverified data. For up to 64 MiB archives use bounded begin/append/finalize_workspace_upload."
    )]
    pub(crate) async fn import_workspace_archive(
        &self,
        Parameters(args): Parameters<ImportArchive>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > (8 * 1024 * 1024usize).div_ceil(3) * 4 {
            return Err(bad("use chunked upload for archives above 8 MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > 8 * 1024 * 1024 {
                return Err(bad(
                    "small native archive import exceeds 8 MiB; use chunked upload",
                ));
            }
            server.session.import_archive_bytes(&bytes)
        })
        .await
    }
    #[tool(
        description = "Export a complete immutable native workspace with exact history, undo/redo, candidate data and attachments. Supply exact imported workspace id/revision. Returns a base64 artifact read in bounded read_artifact chunks; no path/write or trusted verification is granted."
    )]
    pub(crate) async fn export_workspace_archive(
        &self,
        Parameters(args): Parameters<WorkspaceRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let value={let store=server.session.store.lock().unwrap();if args.id!=args.revision{return Err(bad("stale workspace revision"));}store.workspaces.get(&args.id).ok_or_else(||bad("unknown workspace"))?.restore(&store)?};
            let bytes=workspace::save_workspace(&value).map_err(bad)?;let length=bytes.len();let digest=blake3::hash(&bytes).to_hex().to_string();
            let artifact=server.session.binary_artifact(bytes,"application/octet-stream","workspace.cg2workspace")?;
            output(json!({"workspace":args.id,"revision":args.revision,"bytes":length,"blake3":digest,"encoding":"base64","artifact_id":artifact,"filename":"workspace.cg2workspace","candidates_verified":false}))
        }).await
    }
    #[tool(
        description = "Finalize a saved cloning plan against explicit exact historical input/independent-target events in an immutable uploaded native workspace. Freshly reconstructs and retains ALL native step history; options control final construct name, prepared-output presentation, faithful saved oligo attachment and explicit incomplete acceptance. Atomically returns a NEW session workspace, native archive artifact and full receipt JSON; original workspace and sources are unchanged. Hiding preparations never deletes evidence. No browser/hosted writes, guessed DNA event mapping, experimental-material claim or imported verification authority."
    )]
    pub(crate) async fn finalize_cloning_outputs(
        &self,
        Parameters(args): Parameters<FinalizeCloningOutputs>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.input_bindings.len() > 128 {
            return Err(bad("at most128 explicit original input event bindings"));
        }
        args.options.validate().map_err(bad)?;
        let binding_bounded = |binding: &collagene_core::material::lineage::EventBinding| {
            binding.document_id.len() <= 128
                && binding.event_id.len() <= 128
                && binding.revision.len() == 64
        };
        if !binding_bounded(&args.target_binding)
            || args
                .input_bindings
                .iter()
                .any(|binding| binding.input_id.len() > 128 || !binding_bounded(&binding.event))
        {
            return Err(bad("oversized explicit event binding"));
        }
        let server = self.clone();
        self.compute(move || {
            let value = server.session.restore_workspace(&args.workspace)?;
            let plan = {
                let store = server.session.store.lock().unwrap();
                let saved = store
                    .plans
                    .get(&args.plan_revision)
                    .ok_or_else(|| bad("unknown saved plan revision"))?;
                if saved.target.revision() != args.target_binding.revision {
                    return Err(bad(
                        "chosen target event differs from the saved independent target",
                    ));
                }
                saved.plan.clone()
            };
            let expected = collagene_core::material::outputs::finalize_revision(
                &value,
                &plan,
                args.target_binding.clone(),
                args.input_bindings.clone(),
                &args.options,
            )
            .map_err(bad)?;
            let commit = collagene_core::material::outputs::finalize(
                &value,
                &plan,
                args.target_binding,
                args.input_bindings,
                &args.options,
            )
            .map_err(bad)?;
            commit.validate_binding(&expected).map_err(bad)?;
            let text = bounded_json_limit(&commit, 32 * 1024 * 1024)?;
            let receipt = serde_json::to_value(&commit.receipt).map_err(bad)?;
            let native = workspace::save_workspace(&commit.workspace).map_err(bad)?;
            server.session.commit_workspace_outputs(
                commit.workspace,
                Some((text, receipt)),
                Some(native),
            )
        })
        .await
    }
    #[tool(
        description = "Freshly check sequencing reads against an exact construct reference event under explicit native-quality, trimming, placement and support rules, and retain the complete evidence in a NEW immutable native workspace. Reference and read event bindings are explicit; equal DNA never chooses events. Original sequences, heads and chromatograms are unchanged. Atomically returns workspace plus full JSON/native artifacts. Saved results are historical evidence, not restored current verification or experimental sample confirmation."
    )]
    pub(crate) async fn record_read_check(
        &self,
        Parameters(args): Parameters<RecordReadCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.reads.is_empty() || args.reads.len() > 32 {
            return Err(bad("choose1..32 exact read event bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let frame=server.session.restore_workspace(&args.workspace)?;
            let expected=collagene_core::review::reads::checks::capture_revision(&frame,&args.reference,&args.reads,&args.options).map_err(bad)?;
            let commit=collagene_core::review::reads::checks::capture(&frame,args.reference,args.reads,&args.options).map_err(bad)?;
            commit.validate_binding(&expected).map_err(bad)?;
            let saved=collagene_core::review::reads::checks::load(&commit.workspace,&commit.attachment_hash).map_err(bad)?;
            let text=bounded_json_limit(&saved,17*1024*1024)?;
            let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,
                "reference":saved.reference,"read_count":saved.reads.len(),"recorded_status":saved.report.status,
                "recorded_engine":saved.report.engine_version,"state":"historical","verified":false});
            let native=workspace::save_workspace(&commit.workspace).map_err(bad)?;
            server.session.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"read-check.json",
                "Historical source-event-bound read evidence. Fresh replay reconstructs current evidence; no experimental sample identity or whole-construct confirmation is asserted.")
        }).await
    }
    #[tool(
        description = "List saved construct-owned read checks in an immutable native workspace. Every result retains its exact reference event and read count. Recorded engine/status are historical fields; this does not rerun reads or restore verification authority."
    )]
    pub(crate) async fn list_read_checks(
        &self,
        Parameters(args): Parameters<ListReadChecks>,
    ) -> Result<CallToolResult, ErrorData> {
        let frame = self.session.restore_workspace(&args.workspace)?;
        let checks = collagene_core::review::reads::checks::list(&frame).map_err(bad)?;
        output(
            json!({"workspace":args.workspace,"checks":checks,"state":"historical","verified":false}),
        )
    }
    #[tool(
        description = "Replay one saved read check by resolving its exact retained reference/read events and freshly recomputing trimming, placement, coverage and native call evidence with the current engine. Does not change source heads, raw chromatograms, saved evidence or the native workspace. Current computational consistency remains distinct from empirical sample identity."
    )]
    pub(crate) async fn replay_read_check(
        &self,
        Parameters(args): Parameters<ReplayReadCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.attachment_hash.len() != 64 {
            return Err(bad("an exact saved read-check attachment hash is required"));
        }
        let server = self.clone();
        self.compute(move || {
            let frame=server.session.restore_workspace(&args.workspace)?;
            let expected=collagene_core::review::reads::checks::replay_revision(&frame,&args.attachment_hash).map_err(bad)?;
            let checked=collagene_core::review::reads::checks::replay(&frame,&args.attachment_hash).map_err(bad)?;
            checked.validate_binding(&expected).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&checked,32*1024*1024)?,"application/json","read-check-replay.json")?;
            output(json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"context_revision":checked.context_revision,
                "reference":checked.saved.reference,"reads":checked.saved.reads,"state":"fresh-computational-evidence",
                "engine_version":checked.report.engine_version,"status":checked.report.status,
                "eligible_coverage":checked.report.eligible_coverage,"reference_coverage_fraction":checked.report.reference_coverage_fraction,
                "artifact_id":artifact,"verified":false,"authority":"Fresh computational evidence under recorded rules; not experimental sample confirmation."}))
        }).await
    }
    #[tool(
        description = "Preview explicit decisions on a saved construct read check after fresh replay. Nominated bases, calls, deletion evidence and insertion read/boundary/sequence must match current eligible evidence exactly; index-only or unsupported edits are refused. Retain/unresolved decisions do not edit DNA. Original reads/reference remain unchanged; any revised reference is a proposal requiring a new read check, not measured sequence or experimental confirmation."
    )]
    pub(crate) async fn preview_read_review(
        &self,
        Parameters(args): Parameters<ReadReviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;
            let preview=collagene_core::review::reads::review::preview(&frame,&args.request).map_err(bad)?;
            let expected=collagene_core::review::reads::review::preview_revision(&frame,&args.request).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&preview,32*1024*1024)?,"application/json","read-review-preview.json")?;
            let target=preview.target.as_ref().map(|m|json!({"name":m.name,"length":m.sequence.len(),"topology":m.topology,"revision":m.revision()}));
            output(json!({"state":"proposal","verified":false,"preview":{"context_revision":preview.context_revision,"request":preview.request,
                "reference":preview.reference,"reads":preview.reads,"target":target,"edits":preview.edits,"warnings":preview.warnings},"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Freshly replay a saved read check, validate explicit review decisions and atomically retain the decisions and optional derived reference in a NEW immutable native workspace. Exact original reference/read events become computational parents. Source heads and raw calls remain unchanged. Human choices never establish a measurement or current positive read verdict; a derived reference requires a new independent read check. Retain/unresolved-only requests save decisions without creating a reference."
    )]
    pub(crate) async fn commit_read_review(
        &self,
        Parameters(args): Parameters<ReadReviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;
            let expected=collagene_core::review::reads::review::commit_revision(&frame,&args.request).map_err(bad)?;
            let commit=collagene_core::review::reads::review::commit(&frame,&args.request).map_err(bad)?;
            commit.validate_binding(&expected).map_err(bad)?;
            let saved=collagene_core::review::reads::review::load(&commit.workspace,&commit.attachment_hash).map_err(bad)?;
            let text=bounded_json_limit(&saved,32*1024*1024)?;
            let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,"target":commit.target,"state":"historical","verified":false});
            let native=workspace::save_workspace(&commit.workspace).map_err(bad)?;
            server.session.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"read-review.json",
                "Recorded caller decisions and optional proposed reference only. Unselected differences remain unreviewed; source calls are unchanged and the revised reference requires a new read check.")
        }).await
    }
    #[tool(
        description = "List immutable saved decisions on construct read checks and their optional derived-reference events. Historical review records do not restore current validation, experimental confirmation or an approved whole-construct sequence."
    )]
    pub(crate) async fn list_read_reviews(
        &self,
        Parameters(args): Parameters<ListReadChecks>,
    ) -> Result<CallToolResult, ErrorData> {
        let frame = self.session.restore_workspace(&args.workspace)?;
        let reviews = collagene_core::review::reads::review::list(&frame).map_err(bad)?;
        output(
            json!({"workspace":args.workspace,"reviews":reviews,"state":"historical","verified":false}),
        )
    }
    #[tool(
        description = "Reconstruct an exact session plan against explicit historical target/input event bindings and retain all products in a NEW immutable session workspace. Does not alter the original, browser workspace or any file. Incomplete evidence requires explicit allow_incomplete; failed reconstructions never record products. Returns a full JSON artifact and new workspace reference."
    )]
    pub(crate) async fn record_workspace_run(
        &self,
        Parameters(args): Parameters<RecordWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let value = server.session.restore_workspace(&args.workspace)?;
            let plan = server
                .session
                .store
                .lock()
                .unwrap()
                .plans
                .get(&args.plan_revision)
                .ok_or_else(|| bad("unknown plan revision"))?
                .plan
                .clone();
            let commit = collagene_core::material::lineage::record_run(
                &value,
                &plan,
                args.target_binding,
                args.input_bindings,
                args.allow_incomplete,
            )
            .map_err(bad)?;
            let text = bounded_json(&commit)?;
            let receipt = serde_json::to_value(&commit.receipt).map_err(bad)?;
            server
                .session
                .commit_workspace_artifact(commit.workspace, Some((text, receipt)))
        })
        .await
    }
    #[tool(
        description = "List historical cloning runs and exact source events in an immutable workspace. Previews cap100 runs/events and20 products per run; full JSON artifact retains every binding. Stored runs never establish current verification authority."
    )]
    pub(crate) async fn list_workspace_runs(
        &self,
        Parameters(args): Parameters<WorkspaceRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let value = server.session.restore_workspace(&args)?;
            let runs = collagene_core::material::lineage::list_runs(&value).map_err(bad)?;
            let full = json!({"runs":runs,"events":value.history});
            let artifact = server.session.analysis_artifact(&full)?;
            let previews = runs.iter().take(100).map(|run|json!({"run_id":run.run_id,"attachment_hash":run.attachment_hash,"plan_name":run.plan_name.chars().take(256).collect::<String>(),"target":run.target,"product_count":run.products.len(),"products":run.products.iter().take(20).collect::<Vec<_>>(),"products_truncated":run.products.len()>20,"recorded_status":run.recorded_status,"recorded_engine":run.recorded_engine,"input_count":run.input_count})).collect::<Vec<_>>();
            let events = value.history.iter().take(100).map(|event|json!({"document_id":event.document_id,"event_id":event.id,"revision":event.revision,"operation":event.operation.chars().take(256).collect::<String>()})).collect::<Vec<_>>();
            output(json!({"workspace":args,"run_count":runs.len(),"runs":previews,"runs_truncated":runs.len()>100,"event_count":value.history.len(),"events":events,"events_truncated":value.history.len()>100,"artifact_id":artifact,"authority":"historical data; replay recomputes current evidence"}))
        }).await
    }
    #[tool(
        description = "Prepare one exact saved native run for current review using its immutable workspace/run ID. Resolves original input and independent-target events even when current heads changed. Atomically returns full documentary JSON plus a binary candidate artifact and unverified session plan usable with verify_plan. Does not replay/adopt historical status or material, replace current targets, alter workspace heads, or write a native database."
    )]
    pub(crate) async fn prepare_saved_run_review(
        &self,
        Parameters(args): Parameters<ReplayWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.run_id.is_empty() || args.run_id.len() > 128 {
            return Err(bad("explicit bounded native run ID required"));
        }
        let server = self.clone();
        self.compute(move || {
            let workspace = server.session.restore_workspace(&args.workspace)?;
            let expected =
                collagene_core::material::review::request_revision(&workspace, &args.run_id)
                    .map_err(bad)?;
            let review =
                collagene_core::material::review::prepare(&workspace, &args.run_id).map_err(bad)?;
            review.validate_request_binding(&expected).map_err(bad)?;
            server
                .session
                .saved_review_artifacts(review, args.workspace)
        })
        .await
    }
    #[tool(
        description = "Export the canonical feature map of one exact immutable session document as SVG or PNG. Select Svg or Png with explicit pixel width (16–4096, at most 4 million pixels). PNG uses our bundled font and embedded source identity. Complete artifact, exact source revision and renderer warnings are returned. No paths, external renderer, uploaded SVG, restriction overlay or physical/material authority."
    )]
    pub(crate) async fn export_map(
        &self,
        Parameters(args): Parameters<MapExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let reference = args.document;
            let source = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&reference)?;
                collagene_model::bounds::model_preflight(&[source]).map_err(bad)?;
                source.clone()
            };
            let rendered = collagene_core::output::map_export::Request {
                source, source_revision: reference.revision.clone(), restriction_sites: None, format: args.format,
            }.render().map_err(bad)?;
            rendered.validate_binding(&reference.revision).map_err(bad)?;
            let size = rendered.image.bytes().len();
            let artifact = blake3::hash(rendered.image.bytes()).to_hex().to_string();
            let (data, mime_type, encoding, filename) = match rendered.image {
                collagene_core::output::raster::Image::Svg(text) => (ArtifactData::Text(text), "image/svg+xml", "utf-8", "map.svg"),
                collagene_core::output::raster::Image::Png(image) => (ArtifactData::Binary(image.bytes), "image/png", "base64", "map.png"),
            };
            let result = output(json!({"document":reference,"source_revision":rendered.source_revision,"format":args.format,"artifact_id":artifact,"filename":filename,"encoding":encoding,"mime_type":mime_type,"bytes":size,"warnings":rendered.warnings,"warning_count":rendered.warnings.len(),"warnings_truncated":false,"restriction_overlay":"none_requested","scope":"Canonical features of this exact native revision; no restriction-site overlay requested.","authority":"Documentary display only. Map export does not verify a cloning design, prove physical material, alter a source or write a native catalog."}))?;
            let mut store = server.session.store.lock().unwrap();
            if let Some(existing) = store.artifacts.get(&artifact) {
                if !matches!((&existing.data, &data), (ArtifactData::Text(_), ArtifactData::Text(_)) | (ArtifactData::Binary(_), ArtifactData::Binary(_)))
                    || existing.mime_type != mime_type || existing.filename != filename || existing.encoding != encoding {
                    return Err(bad("map artifact already has another representation"));
                }
            } else {
                if store.artifacts.len() >= MAX_ITEMS { return Err(bad("at most 128 artifacts; no map published")); }
                store.reserve(size)?;
                store.artifacts.insert(artifact, Artifact { data, mime_type, encoding, filename });
            }
            Ok(result)
        }).await
    }
    #[tool(
        description = "Replay a retained run from exact historical sources using the current embedded verifier. Returns bounded product previews and a complete JSON evidence artifact. Engine/status/product comparisons are explicit. Does not edit a workspace, trust a stored report or invoke external programs."
    )]
    pub(crate) async fn replay_workspace_run(
        &self,
        Parameters(args): Parameters<ReplayWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let value = server.session.restore_workspace(&args.workspace)?;
            let report = collagene_core::material::lineage::replay_run(&value, &args.run_id).map_err(bad)?;
            report.validate_for_workspace(&value).map_err(bad)?;
            let full = serde_json::to_value(&report).map_err(bad)?;
            let artifact = server.session.analysis_artifact(&full)?;
            let products = report.current_report.products.iter().take(100).map(|p| json!({"revision":p.revision(),"name":p.name.chars().take(256).collect::<String>(),"length":p.sequence.len()})).collect::<Vec<_>>();
            output(json!({"workspace":args.workspace,"run_id":report.run_id,"source_workspace_fingerprint":report.source_workspace_fingerprint,"recorded_status":report.recorded_status,"recorded_engine":report.recorded_engine,"current_status":report.current_report.status,"current_engine":report.current_report.engine_version,"engine_matches":report.engine_matches,"products_match":report.products_match,"status_matches":report.status_matches,"products":products,"product_count":report.current_report.products.len(),"products_truncated":report.current_report.products.len()>100,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Read immutable historical reference-design and processed physical graph attachments for a retained run. Checks exact run/plan/engine/source/target/product-event bindings but does not replay chemistry or restore physical authority. Older runs may have no historical graph; fresh replay is a separate tool."
    )]
    pub(crate) async fn inspect_workspace_run_evidence(
        &self,
        Parameters(args): Parameters<ReplayWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let workspace=server.session.restore_workspace(&args.workspace)?;
            let record=collagene_core::material::lineage::load_run(&workspace,&args.run_id).map_err(bad)?;
            let evidence=collagene_core::material::evidence::download(&workspace,&args.run_id).map_err(bad)?;
            let artifact=evidence.bytes.as_ref().map(|bytes|server.session.text_artifact(String::from_utf8(bytes.clone()).map_err(bad)?,"application/json","historical-cloning-graphs.json")).transpose()?;
            output(json!({"state":"historical","verified":false,"run_id":record.run_id,"recorded_engine":record.assessment.engine_version,"scientific_step_count":evidence.step_count,"artifact_id":artifact,"authority":"Stored historical graph evidence only. Fresh replay of exact original inputs is required for current material state or verification."}))
        }).await
    }
    #[tool(
        description = "Reserve a session upload for exact native workspace bytes. Length 1..64 MiB, lower-case BLAKE3 required. At most four uploads; full expected length is charged to the shared 64 MiB session budget immediately. No filesystem or external programs."
    )]
    pub(crate) fn begin_workspace_upload(
        &self,
        Parameters(args): Parameters<BeginWorkspaceUpload>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.expected_bytes == 0
            || args.expected_bytes > workspace::MAX_ARCHIVE_BYTES
            || args.blake3.len() != 64
            || !args
                .blake3
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(bad("invalid archive length or lower-case BLAKE3"));
        }
        let mut store = self.session.store.lock().unwrap();
        if store.uploads.len() >= MAX_UPLOADS {
            return Err(bad("at most four outstanding workspace uploads"));
        }
        let serial = store
            .upload_serial
            .checked_add(1)
            .ok_or_else(|| bad("upload counter exhausted"))?;
        store.reserve(args.expected_bytes)?;
        let mut bytes = Vec::new();
        if let Err(error) = bytes.try_reserve_exact(args.expected_bytes) {
            store.used -= args.expected_bytes;
            return Err(bad(format!("upload allocation failed: {error}")));
        }
        store.upload_serial = serial;
        let id = format!("workspace-upload-{serial}");
        store.uploads.insert(
            id.clone(),
            WorkspaceUpload {
                expected: args.expected_bytes,
                hash: args.blake3,
                bytes,
            },
        );
        output(
            json!({"upload_id":id,"next_offset":0,"expected_bytes":args.expected_bytes,"max_decoded_chunk_bytes":MAX_UPLOAD_CHUNK}),
        )
    }
    #[tool(
        description = "Append one independently base64-encoded native workspace chunk, <=256 KiB decoded. Offset must equal the exact next decoded offset; gaps/overlaps/replays/overflow are rejected without changing the upload. No raw archive is parsed before finalize."
    )]
    pub(crate) fn append_workspace_upload(
        &self,
        Parameters(args): Parameters<AppendWorkspaceUpload>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > MAX_UPLOAD_CHUNK.div_ceil(3) * 4 {
            return Err(bad("upload chunk exceeds 256 KiB decoded"));
        }
        let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
        if bytes.is_empty() || bytes.len() > MAX_UPLOAD_CHUNK {
            return Err(bad("upload chunk must contain 1..256 KiB decoded bytes"));
        }
        let mut store = self.session.store.lock().unwrap();
        let upload = store
            .uploads
            .get_mut(&args.upload_id)
            .ok_or_else(|| bad("unknown upload"))?;
        if args.offset != upload.bytes.len()
            || bytes.len() > upload.expected.saturating_sub(upload.bytes.len())
        {
            return Err(bad(
                "incorrect next offset or chunk exceeds expected archive length",
            ));
        }
        upload.bytes.extend_from_slice(&bytes);
        output(
            json!({"upload_id":args.upload_id,"next_offset":upload.bytes.len(),"expected_bytes":upload.expected,"complete":upload.bytes.len()==upload.expected}),
        )
    }
    #[tool(
        description = "Finalize an exactly complete native workspace upload. Verifies length and BLAKE3 before parsing, then imports all revisions/history/candidates/attachments atomically under core64MiB archive/128MiB expanded and session64MiB limits. Incomplete uploads remain appendable; completed attempts consume/release upload storage even on hash/parse/storage failure. Imported candidates are never authoritative verification."
    )]
    pub(crate) async fn finalize_workspace_upload(
        &self,
        Parameters(args): Parameters<WorkspaceUploadRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let upload = {
                let mut store = server.session.store.lock().unwrap();
                let upload = store
                    .uploads
                    .get(&args.upload_id)
                    .ok_or_else(|| bad("unknown upload"))?;
                if upload.bytes.len() != upload.expected {
                    return Err(bad("upload incomplete; append remaining bytes or abort"));
                }
                store.uploads.remove(&args.upload_id).unwrap()
            };
            let decoded = if blake3::hash(&upload.bytes).to_hex().as_str() != upload.hash {
                Err(bad("archive BLAKE3 mismatch; upload consumed"))
            } else {
                workspace::load_workspace(&upload.bytes).map_err(bad)
            };
            // The reservation remains charged throughout hashing/decompression. Imports
            // are committed only after all core validation, with no partial documents.
            server.session.store.lock().unwrap().used -= upload.expected;
            drop(upload);
            server.session.commit_workspace(decoded?)
        })
        .await
    }
    #[tool(
        description = "Abort a pending workspace upload and release its full reserved session bytes. No documents or candidates are imported by abort."
    )]
    pub(crate) fn abort_workspace_upload(
        &self,
        Parameters(args): Parameters<WorkspaceUploadRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let mut store = self.session.store.lock().unwrap();
        let upload = store
            .uploads
            .remove(&args.upload_id)
            .ok_or_else(|| bad("unknown upload"))?;
        store.used -= upload.expected;
        output(json!({"upload_id":args.upload_id,"aborted":true,"released_bytes":upload.expected}))
    }
    #[tool(
        description = "List immutable imported workspace summaries with unverified candidate metadata. Exact workspace id/revision is required to export; molecule revisions are shared and accessible through paginated list_documents."
    )]
    pub(crate) fn list_workspaces(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.session.store.lock().unwrap();
        output(
            json!({"workspaces":store.workspaces.iter().map(|(id,w)|w.summary(id)).collect::<Vec<_>>()}),
        )
    }
    #[tool(
        description = "Import an exact binary .cg2design snapshot from the browser, including annotations/topology/revisions. Supply base64 file bytes, at most 8 MiB decoded. Returns the independently selected target and immutable document references; no plans are applied."
    )]
    pub(crate) async fn import_workspace(
        &self,
        Parameters(args): Parameters<ImportWorkspace>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > 12 * 1024 * 1024 {
            return Err(bad("planning snapshot exceeds MCP exchange limit"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > 8 * 1024 * 1024 {
                return Err(bad(
                    "planning snapshots are currently limited to 8 MiB over MCP",
                ));
            }
            let workspace = transfer::load_workspace(&bytes).map_err(bad)?;
            let target = DocumentRef {
                id: workspace.target_revision.clone(),
                revision: workspace.target_revision,
            };
            let mut result = server
                .session
                .insert_documents(workspace.documents, vec![])?;
            if let Some(value) = result.structured_content.as_mut() {
                value["target"] = json!(target);
            }
            // Rebuild the standard text fallback as well as structuredContent.
            output(result.structured_content.unwrap())
        })
        .await
    }
    #[tool(
        description = "Preview or install target-free native annealed oligo adapters, exact PCR, intentional-mismatch PCR, overlap-extension PCR, catalog digestion, end treatment, single-base 3′ tailing or saved-material end repair from exact current source events in an immutable workspace. Template methods require one template role; annealed adapters require explicit top_oligo and bottom_oligo roles, which may refer to the same self-complementary event. Preview tasks return bounded metadata plus complete reference/geometry artifacts; they create no material documents. Install freshly reconstructs actual source strands, validates the literal selected PCR geometry or actual digest fragment and atomically returns a NEW session workspace, receipt and native archive, with optional newly performed PCR primers. No independent target, isolation, yield, experimental success, database write or filesystem authority is fabricated. Pass saved-origin MaterialSelection handles and explicit PCR preview index/revision or OE context/product revisions only; closed capsules and whole amplicon geometry are resolved natively, never copied into requests. Saved Run/Derived handles must match their exact retained source events; marked reference-only documents cannot bypass replay as plain originals."
    )]
    pub(crate) async fn material_derivation(
        &self,
        Parameters(args): Parameters<MaterialDerivationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::material::derivation as native;
            #[derive(Serialize)]
            struct MaterialPreviewEvidence<'a,T:Serialize,P:Serialize>{workspace:&'a WorkspaceRef,source:&'a collagene_core::material::lineage::EventBinding,parameters:P,preview:&'a T}
            let mut workspace=server.session.restore_workspace(&args.workspace)?;
            match args.task {
                MaterialDerivationTask::PreviewAnnealedOligos{top,bottom,options}=>{
                    for event in [&top,&bottom] {
                        server.session.material_derivation_source(&workspace,event,&native::SourceOrigin::OriginalReference)?;
                    }
                    let context=native::annealing::annealed_oligos_preview_revision(&workspace,&top,&bottom,&options).map_err(bad)?;
                    let preview=native::annealing::preview_annealed_oligos(&workspace,&top,&bottom,&options).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    #[derive(Serialize)]
                    struct AnnealedMaterialEvidence<'a>{workspace:&'a WorkspaceRef,preview:&'a native::annealing::AnnealedOligosPreview}
                    let artifact=server.session.text_artifact(bounded_json_limit(&AnnealedMaterialEvidence{workspace:&args.workspace,preview:&preview},32*1024*1024)?,"application/json","material-annealed-oligos-preview.json")?;
                    let candidates=preview.preview.annealing.candidates.iter().zip(&preview.products).take(100).map(|(candidate,product)|{
                        Ok(json!({"candidate_id":candidate.id,"paired_bases":candidate.paired_bases,"duplex":candidate.duplex.summary().map_err(bad)?,"product_revision":product.revision(),"selection":native::AnnealSelection{expected_preview_revision:context.clone(),candidate_id:candidate.id.clone(),expected_product_revision:product.revision()}}))
                    }).collect::<Result<Vec<_>,ErrorData>>()?;
                    output(json!({"workspace":args.workspace,"sources":preview.sources,"context_revision":context,"complete":preview.preview.annealing.complete,"unique":preview.preview.annealing.unique,"count":preview.products.len(),"candidates":candidates,"candidates_truncated":preview.products.len()>100,"incomplete":preview.incomplete,"findings":preview.findings,"artifact_id":artifact,"authority":"Explicit standard-unmodified single-strand oligo hypotheses bound to both exact source events. Native annealing retains the selected geometry, alternatives and search limitations. Preview creates no material, isolation or covalent sealing. Installation independently replays both sources and the selected native candidate; no intended assembly target is invented."}))
                },
                MaterialDerivationTask::PreviewTerminalProcessing{source,origin,action}=>{
                    let origin=server.session.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    server.session.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::terminal::terminal_processing_preview_revision(&workspace,&source,&origin,&action).map_err(bad)?;
                    let preview=native::terminal::preview_terminal_processing(&workspace,&source,&origin,&action).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.session.text_artifact(bounded_json_limit(&MaterialPreviewEvidence{workspace:&args.workspace,source:&source,parameters:&action,preview:&preview},32*1024*1024)?,"application/json","material-terminal-processing-preview.json")?;
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"action":action,"product":{"revision":preview.product.revision(),"length":preview.product.sequence.len(),"topology":preview.product.topology},"source_duplex":preview.source_duplex,"duplex":preview.duplex,"source_modifications":preview.source_modifications,"modifications":preview.modifications,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Target-free explicit terminal treatment or single-base 3′ tailing preview on freshly reconstructed source strands. Native end-geometry guards and positional modifications remain in force. Installation independently executes the selected action. No ligation, target agreement, nick sealing, new primers, isolation or experimental yield is established."}))
                },
                MaterialDerivationTask::PreviewOverlapExtensionPcr{source,origin,preparation}=>{
                    let origin=server.session.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    server.session.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::overlap_extension::overlap_extension_preview_revision(&workspace,&source,&origin,&preparation).map_err(bad)?;
                    let preview=native::overlap_extension::preview_overlap_extension(&workspace,&source,&origin,&preparation).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.session.text_artifact(bounded_json_limit(&MaterialPreviewEvidence{workspace:&args.workspace,source:&source,parameters:&preparation,preview:&preview},32*1024*1024)?,"application/json","material-overlap-extension-pcr-preview.json")?;
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"reusable":preview.reusable(),"selection":preview.selection().ok(),"product":preview.product.as_ref().map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Source-bound target-free OE PCR preview only. Both source PCR branches, reciprocal extension and outer PCR are freshly executed. Explicit installation repeats native execution and retains all six new primer uses; no independent target, isolation, yield or experimental success is established."}))
                },
                MaterialDerivationTask::PreviewPcr{source,origin,forward,reverse,template_strand}=>{
                    let origin=server.session.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.session.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::pcr::pcr_preview_revision(&molecule,&origin,&forward,&reverse,template_strand).map_err(bad)?;
                    let documents=native::pcr::preview_pcr(&molecule,&origin,&forward,&reverse,template_strand).map_err(bad)?;
                    #[derive(Serialize)]
                    struct PcrEvidence<'a>{workspace:&'a WorkspaceRef,source:&'a collagene_core::material::lineage::EventBinding,context_revision:&'a str,forward:&'a Primer,reverse:&'a Primer,template_strand:collagene_core::chemistry::Strand,documents:&'a[verification::pcr::PcrDocument]}
                    let text=bounded_json_limit(&PcrEvidence{workspace:&args.workspace,source:&source,context_revision:&context,forward:&forward,reverse:&reverse,template_strand,documents:&documents},32*1024*1024)?;
                    let artifact=server.session.text_artifact(text,"application/json","material-pcr-preview.json")?;
                    let rows=documents.iter().take(100).enumerate().map(|(index,d)|json!({"index":index,"preview_revision":d.molecule.revision(),"length":d.molecule.sequence.len(),"reusable":documents.len()==1 && d.molecule.sequence.as_str().bytes().all(|b|b"ACGT".contains(&b)) && forward.ordered_sequence().bytes().chain(reverse.ordered_sequence().bytes()).all(|b|b"ACGT".contains(&b)),"template_strand":d.geometry.template_strand,"forward_start":d.geometry.forward_start,"reverse_start":d.geometry.reverse_start,"crosses_origin":d.geometry.crosses_origin,"forward_primer":d.geometry.forward_primer,"reverse_primer":d.geometry.reverse_primer})).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"count":documents.len(),"previews":rows,"previews_truncated":documents.len()>100,"artifact_id":artifact,"authority":"Reference PCR previews only. Use the explicit preview index and revision for Install; the server recomputes full geometry internally. Full geometry remains optional artifact evidence. Multiple products do not establish isolation or reusable material."}))
                },
                MaterialDerivationTask::PreviewMismatchPcr{source,origin,forward,reverse,options}=>{
                    let origin=server.session.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.session.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::pcr::mismatch_pcr_preview_revision(&molecule,&origin,&forward,&reverse,&options).map_err(bad)?;
                    let preview=native::pcr::preview_mismatch_pcr(&molecule,&origin,&forward,&reverse,&options).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.session.text_artifact(bounded_json_limit(&MaterialPreviewEvidence{workspace:&args.workspace,source:&source,parameters:json!({"forward":forward,"reverse":reverse,"options":options}),preview:&preview},32*1024*1024)?,"application/json","material-mismatch-pcr-preview.json")?;
                    let report=&preview.preview.report;
                    let reusable=preview.reusable();
                    let rows=preview.preview.documents.iter().zip(&report.products).take(100).enumerate().map(|(index,(document,geometry))|json!({"index":index,"preview_revision":document.revision(),"length":document.sequence.len(),"certain":geometry.certain,"overlapping_domains":geometry.overlapping_domains,"forward_start":geometry.forward_start,"reverse_start":geometry.reverse_start,"crosses_origin":geometry.crosses_origin,"forward_primer":geometry.forward_primer_index,"reverse_primer":geometry.reverse_primer_index,"substitution_count":geometry.substitutions.len(),"reusable":reusable&&index==0})).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"count":preview.preview.documents.len(),"reusable":reusable,"previews":rows,"previews_truncated":preview.preview.documents.len()>100,"search_complete":report.search_complete,"enumeration_complete":report.enumeration_complete,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Source-bound mismatch-PCR reference preview with full literal ordered primers and explicit mismatch/3-prime/search policy. Install freshly recomputes the selected geometry and independently refuses unsupported or competing physical products; this preview grants no material or isolation."}))
                },
                MaterialDerivationTask::PreviewRepair{source,origin,method}=>{
                    let origin=server.session.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.session.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::end_repair::repair_preview_revision(&molecule,&origin,method).map_err(bad)?;
                    let preview=native::end_repair::preview_repair(&molecule,&origin,method).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.session.text_artifact(bounded_json_limit(&MaterialPreviewEvidence{workspace:&args.workspace,source:&source,parameters:json!({"method":method}),preview:&preview},32*1024*1024)?,"application/json","material-end-repair-preview.json")?;
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"method":method,"length":preview.molecule.sequence.len(),"preview_revision":preview.molecule.revision(),"source_duplex":preview.source_duplex,"source_modifications":preview.source_modifications,"repair":preview.repair,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Fresh saved-material strand repair under the explicit method; actual synthesized/removed bases and inherited findings remain evidence. No independent target, isolated stock, experimental success or saved product is asserted by this preview."}))
                },
                MaterialDerivationTask::PreviewDigest{source,origin,selection,options}=>{
                    let origin=server.session.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.session.material_derivation_source(&workspace,&source,&origin)?;
                    let role=collagene_core::cloning::restriction::history::source_role(&source).map_err(bad)?;
                    let context=native::digest::digest_inspection_revision(&role,&molecule,&origin,&selection,&options).map_err(bad)?;
                    let preview=native::digest::inspect_digest(&role,&molecule,&origin,&selection,&options).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.session.text_artifact(bounded_json_limit(&preview,32*1024*1024)?,"application/json","material-digest-preview.json")?;
                    output(json!({"workspace":args.workspace,"source":source,"source_role":role,"context_revision":context,"products_determined":preview.inspection.products_determined,"scan_complete":preview.inspection.scan_complete,"fragment_count":preview.inspection.fragments.len(),"fragments":preview.inspection.fragments.iter().take(100).collect::<Vec<_>>(),"fragments_truncated":preview.inspection.fragments.len()>100,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Fresh actual source-strand catalog inspection under declared complete cleavage. This creates no material document, selected physical stock or independent target assessment."}))
                },
                MaterialDerivationTask::Install{request}=>{
                    if request.name.trim().is_empty()||request.name.trim()!=request.name||request.name.len()>256||request.name.chars().any(char::is_control){return Err(bad("Declare a bounded nonempty material product name."));}
                    use collagene_core::material::derivation::SourceRole;
                    let expected_roles=if matches!(&request.derivation,MaterialDerivationOperation::AnnealedOligos{..}) {vec![SourceRole::TopOligo,SourceRole::BottomOligo]} else {vec![SourceRole::Template]};
                    if request.sources.len()!=expected_roles.len()||request.sources.iter().zip(&expected_roles).any(|(source,role)|source.role!=*role){return Err(bad("Declare the exact canonical source roles for the selected material operation."));}
                    if matches!(&request.derivation,MaterialDerivationOperation::AnnealedOligos{..})&&request.sources.iter().any(|source|!matches!(&source.origin,MaterialDerivationOrigin::OriginalReference)){return Err(bad("Annealing requires explicit original single-strand oligo declarations; saved duplex material is not implicitly denatured."));}
                    let mut sources=Vec::with_capacity(request.sources.len());
                    for source in request.sources {
                        let origin=server.session.resolve_material_derivation_origin(&workspace,&source.event,&source.origin)?;
                        server.session.material_derivation_source(&workspace,&source.event,&origin)?;
                        sources.push(native::BoundSource{role:source.role,event:source.event,origin});
                    }
                    let template_context=||{
                        let source=sources.first().filter(|source|source.role==SourceRole::Template).ok_or_else(||bad("This operation requires the exact template source role."))?;
                        let molecule=server.session.material_derivation_source(&workspace,&source.event,&source.origin)?;
                        Ok::<_,ErrorData>((molecule,&source.origin))
                    };
                    let derivation=match request.derivation {
                        MaterialDerivationOperation::ExactPcr{forward,reverse,selection,template_strand}=>{
                            let (molecule,origin)=template_context()?;
                            if selection.selected_index>=10000||selection.expected_preview_revision.len()!=64||!selection.expected_preview_revision.bytes().all(|b|b.is_ascii_digit()||(b'a'..=b'f').contains(&b)){return Err(bad("Invalid explicit PCR preview index/revision."));}
                            let previews=native::pcr::preview_pcr(&molecule,&origin,&forward,&reverse,template_strand).map_err(bad)?;
                            let chosen=previews.into_iter().nth(selection.selected_index).ok_or_else(||bad("The explicitly selected PCR preview no longer exists."))?;
                            if chosen.molecule.revision()!=selection.expected_preview_revision {return Err(bad("The explicitly selected PCR preview changed; preview again."));}
                            native::Derivation::ExactPcr{forward,reverse,template_strand,selection:native::PcrSelection{selected_index:selection.selected_index,expected_preview_revision:selection.expected_preview_revision,geometry:chosen.geometry}}
                        },
                        MaterialDerivationOperation::MismatchPcr{forward,reverse,options,selection}=>{
                            let (molecule,origin)=template_context()?;
                            if selection.selected_index>=10000||selection.expected_preview_revision.len()!=64||!selection.expected_preview_revision.bytes().all(|b|b.is_ascii_digit()||(b'a'..=b'f').contains(&b)){return Err(bad("Invalid explicit mismatch-PCR preview index/revision."));}
                            let preview=native::pcr::preview_mismatch_pcr(&molecule,&origin,&forward,&reverse,&options).map_err(bad)?;
                            let chosen=preview.preview.documents.get(selection.selected_index).ok_or_else(||bad("The explicitly selected mismatch-PCR preview no longer exists."))?;
                            if chosen.revision()!=selection.expected_preview_revision{return Err(bad("The explicitly selected mismatch-PCR preview changed; preview again."));}
                            let geometry=preview.preview.report.products.get(selection.selected_index).cloned().ok_or_else(||bad("The selected mismatch-PCR geometry is unavailable."))?;
                            native::Derivation::MismatchPcr{forward,reverse,options,selection:native::MismatchPcrSelection{selected_index:selection.selected_index,expected_preview_revision:selection.expected_preview_revision,geometry}}
                        },
                        MaterialDerivationOperation::OverlapExtensionPcr{preparation,selection}=>native::Derivation::OverlapExtensionPcr{preparation,selection},
                        MaterialDerivationOperation::TerminalProcessing{action}=>native::Derivation::TerminalProcessing{action},
                        MaterialDerivationOperation::AnnealedOligos{options,selection}=>native::Derivation::AnnealedOligos{options,selection},
                        MaterialDerivationOperation::EndRepair{method}=>native::Derivation::EndRepair{method},
                        MaterialDerivationOperation::CatalogDigest{selection,options,fragment_index}=>native::Derivation::CatalogDigest{selection,options,fragment_index},
                    };
                    let request=native::Request{sources,derivation,name:request.name,save_primers:request.save_primers};
                    // Selecting the requested retained source in this local copy
                    // changes presentation only; the stored input is immutable.
                    workspace.active_document=request.sources[0].event.document_id.clone();
                    let expected=native::revision(&workspace,&request).map_err(bad)?;
                    let installed=native::prepare_install(&workspace,&request).map_err(bad)?;
                    installed.receipt.validate_binding(&expected).map_err(bad)?;
                    installed.receipt.validate_archive(&installed.workspace).map_err(bad)?;
                    let receipt=serde_json::to_value(&installed.receipt).map_err(bad)?;
                    let text=bounded_json_limit(&installed.receipt,32*1024*1024)?;
                    let archive=workspace::save_workspace(&installed.workspace).map_err(bad)?;
                    let result=server.session.commit_workspace_outputs_named(installed.workspace,Some((text,receipt)),Some(archive),"material-derivation.json","Target-free material derivation in a NEW immutable session workspace only; every reuse reconstructs native source history again.")?;
                    let mut value=result.structured_content.ok_or_else(||bad("missing installed workspace receipt"))?;
                    value["source_workspace"]=json!(args.workspace);
                    value["source_selection"]=json!("Only the derived working copy selects the requested source tab before installation; all source biological heads and the input session workspace remain unchanged.");
                    value["authority"]=json!("Target-free declared-model material derivation saved to a NEW immutable session workspace/archive only. No independent target was invented or assessed, and no experimental isolation/yield, local catalog write, filesystem write or external publication occurred. Reuse freshly replays the complete native record.");
                    output(value)
                },
            }
        }).await
    }
    #[tool(
        description = "Discover or explicitly replay a saved run product or generic target-free material derivation origin for an exact event in an immutable native workspace. No origin is automatically selected. Returns bounded origin choices, findings and complete immutable inspection JSON including any selected closed documentary origin. Every later material preparation must freshly replay again; imported summaries and equal DNA are not material authority. Creates no documents, candidates, targets or stock."
    )]
    pub(crate) async fn inspect_material_origin(
        &self,
        Parameters(args): Parameters<MaterialOriginQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            if args.selection.as_ref().is_some_and(|selection|selection.event()!=&args.event) {
                return Err(bad("selected origin belongs to a different exact event"));
            }
            let workspace=server.session.restore_workspace(&args.workspace)?;
            let expected=collagene_core::material::run::material_inspection_revision(&workspace,&args.event,&args.selection).map_err(bad)?;
            let inspection=collagene_core::material::run::inspect_material(&workspace,&args.event,&args.selection).map_err(bad)?;
            inspection.validate_binding(&expected).map_err(bad)?;
            if let Some(selected)=&inspection.selected {
                selected.validate_for_workspace(&workspace,&args.event).map_err(bad)?;
                if args.selection.as_ref()!=Some(&selected.selection()) { return Err(bad("replayed origin differs from the explicit selection")); }
            }
            let text=bounded_json_limit(&inspection,32*1024*1024)?;
            let artifact=server.session.text_artifact(text,"application/json","material-origin-inspection.json")?;
            output(json!({"workspace":args.workspace,"event":inspection.event,"context_revision":inspection.context_revision,"engine_version":inspection.engine_version,"origins":inspection.origins,"selection":args.selection,"origin_available":inspection.selected.is_some(),"findings":inspection.findings,"artifact_id":artifact,"authority":"A closed documentary replay request is returned only after fresh native reconstruction. Its serialized bytes confer no private material authority; every later preparation independently replays against exact retained sources."}))
        }).await
    }
    #[tool(
        description = "Discover documentary originating run/step choices for an exact saved DNA event in an explicitly supplied immutable native workspace. Returns every bounded origin selection and immutable full inspection JSON. No origin is automatically selected, no chemistry is restored and no document, candidate or target is created. Equal DNA content is not an event identity."
    )]
    pub(crate) async fn list_run_material_origins(
        &self,
        Parameters(args): Parameters<RunMaterialOriginsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let workspace=server.session.restore_workspace(&args.workspace)?;
            let expected=collagene_core::material::run::inspection_revision(&workspace,&args.event,&None).map_err(bad)?;
            let inspection=collagene_core::material::run::inspect(&workspace,&args.event,&None).map_err(bad)?;
            inspection.validate_binding(&expected).map_err(bad)?;
            let text=bounded_json_limit(&inspection,32*1024*1024)?;
            let artifact=server.session.text_artifact(text,"application/json","run-material-origins.json")?;
            output(json!({"workspace":args.workspace,"event":inspection.event,"context_revision":inspection.context_revision,"engine_version":inspection.engine_version,"origins":inspection.origins,"selected":false,"findings":inspection.findings,"artifact_id":artifact,"authority":"Documentary origin discovery only. No private material, successful replay or independent target acceptance is imported."}))
        }).await
    }
    #[tool(
        description = "Explicitly select an exact originating run/step/event in an immutable native workspace and freshly build its bounded durable reconstruction capsule. Rechecks original source/target/product event ownership and current deterministic replay. Returns documentary capsule JSON, not transferable private tracked material or scientific authority. Creates no document/candidate/target and performs no automatic origin choice. Existing source-bound recipes can consume the capsule through explicit FromRun preparation and freshly verify again."
    )]
    pub(crate) async fn prepare_run_material(
        &self,
        Parameters(args): Parameters<RunMaterialPreparationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let workspace=server.session.restore_workspace(&args.workspace)?;
            let selected=Some(args.selection.clone());
            let expected=collagene_core::material::run::inspection_revision(&workspace,&args.selection.event,&selected).map_err(bad)?;
            let inspection=collagene_core::material::run::inspect(&workspace,&args.selection.event,&selected).map_err(bad)?;
            inspection.validate_binding(&expected).map_err(bad)?;
            if inspection.selected.as_ref().is_none_or(|capsule|capsule.selection!=args.selection){return Err(bad("explicit selected origin produced no matching reconstruction capsule"));}
            let text=bounded_json_limit(&inspection,32*1024*1024)?;
            let artifact=server.session.text_artifact(text,"application/json","run-material-preparation.json")?;
            output(json!({"workspace":args.workspace,"event":inspection.event,"context_revision":inspection.context_revision,"engine_version":inspection.engine_version,"selection":args.selection,"capsule_available":true,"findings":inspection.findings,"artifact_id":artifact,"authority":"Freshly built documentary reconstruction request only; the JSON capsule does not carry private tracked material. Every later recipe independently reconstructs from the exact retained original run inputs and compares its separate target."}))
        }).await
    }
}
