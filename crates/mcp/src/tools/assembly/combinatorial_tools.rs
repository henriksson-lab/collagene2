use crate::server::CollageneServer;

use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::references::WorkspaceRef;
use collagene_core::documents::workspace;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;
use serde_json::Value;
use serde_json::json;
use std::collections::BTreeSet;

use collagene_core::planning::combinatorial::{self as design, batch};
use collagene_core::planning::combinatorial::{export as archive, release};

/// Session-resolved source selectors: callers never supply closed material bindings.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LibraryMaterialBinding {
    pub event: collagene_core::material::lineage::EventBinding,
    pub selection: collagene_core::material::run::MaterialSelection,
}
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LibraryDesign {
    pub id: String,
    pub name: String,
    pub circular: bool,
    pub bins: Vec<design::Bin>,
    pub method: design::Method,
    pub material_bindings: Vec<LibraryMaterialBinding>,
}
impl LibraryDesign {
    fn into_native(self, frame: &workspace::Workspace) -> Result<design::Design, ErrorData> {
        bounded_json_limit(&self, 32 * 1024 * 1024)?;
        if self.bins.is_empty()
            || self.bins.len() > design::MAX_BINS
            || self.material_bindings.len() > design::MAX_CHOICES
            || self.bins.iter().any(|b| b.choices.len() > 64)
            || self.bins.iter().map(|b| b.choices.len()).sum::<usize>() > design::MAX_CHOICES
        {
            return Err(bad("bounded source library required"));
        }
        let mut pins = BTreeSet::new();
        for b in &self.material_bindings {
            b.selection.validate().map_err(bad)?;
            if b.selection.event() != &b.event
                || !pins.insert((&b.event.document_id, &b.event.event_id, &b.event.revision))
            {
                return Err(bad(
                    "duplicate or mismatched exact material source selector",
                ));
            }
        }
        let selections = self
            .material_bindings
            .iter()
            .map(|b| b.selection.clone())
            .collect::<Vec<_>>();
        let origins =
            collagene_core::material::capsule::documentary_origins(frame, &selections, &[])
                .map_err(bad)?
                .materials;
        if origins.len() != self.material_bindings.len() {
            return Err(bad("native material hydration count differs"));
        }
        let native = design::Design {
            id: self.id,
            name: self.name,
            circular: self.circular,
            bins: self.bins,
            method: self.method,
            material_bindings: self
                .material_bindings
                .into_iter()
                .zip(origins)
                .map(|(b, origin)| design::LibraryMaterialBinding {
                    event: b.event,
                    origin,
                })
                .collect(),
        };
        design::source_preflight(&native, frame).map_err(bad)?;
        Ok(native)
    }
}
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LibraryRequest {
    pub design: Option<LibraryDesign>,
    pub rows: Vec<batch::RowRequest>,
}
impl LibraryRequest {
    fn bounds(&self) -> Result<(), ErrorData> {
        if self.rows.is_empty() || self.rows.len() > design::MAX_SELECTED {
            return Err(bad("select1–256 current library rows"));
        }
        bounded_json_limit(self, batch::MAX_RECORD_BYTES)?;
        Ok(())
    }
    fn into_native(self, frame: &workspace::Workspace) -> Result<batch::Request, ErrorData> {
        self.bounds()?;
        let native = batch::Request {
            design: self.design.map(|d| d.into_native(frame)).transpose()?,
            rows: self.rows,
        };
        native.validate().map_err(bad)?;
        Ok(native)
    }
}
fn page_size() -> usize {
    20
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MatrixQuery {
    pub workspace: WorkspaceRef,
    pub design: LibraryDesign,
    #[serde(default)]
    pub first: u128,
    #[serde(default = "page_size")]
    pub limit: usize,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectQuery {
    pub workspace: WorkspaceRef,
    pub design: LibraryDesign,
    pub construct: design::SelectedConstruct,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BatchQuery {
    pub workspace: WorkspaceRef,
    pub matrix: LibraryRequest,
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
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReplayQuery {
    pub workspace: WorkspaceRef,
    pub attachment_hash: String,
    pub row_index: usize,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReleaseQuery {
    pub workspace: WorkspaceRef,
    pub request: release::Request,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportQuery {
    pub workspace: WorkspaceRef,
    pub manifest_hash: String,
    #[serde(default)]
    pub options: archive::Options,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReleaseRef {
    pub workspace: WorkspaceRef,
    pub manifest_hash: String,
}

fn hash(value: &str) -> Result<(), ErrorData> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(bad("an exact lower-case batch attachment hash is required"));
    }
    Ok(())
}
fn summaries(resolved: &design::Resolved) -> Value {
    json!({"construct_id":resolved.construct_id,"combination_id":resolved.combination_id,
        "declaration_revision":resolved.declaration_revision,"input_bindings":resolved.bindings,
        "inputs":resolved.inputs.iter().map(|input|json!({"input_id":input.id,"revision":input.revision,"name":input.molecule.name,"length":input.molecule.sequence.len(),"topology":input.molecule.topology})).collect::<Vec<_>>(),
        "target":{"name":resolved.target.name,"revision":resolved.target.revision(),"length":resolved.target.sequence.len(),"topology":resolved.target.topology,"binding":resolved.target_binding},
        "state":"independent-declaration-and-recipe","verified":false})
}
fn row_summary(row: &batch::RowEvidence) -> Value {
    let error = row
        .error
        .as_ref()
        .map(|error| error.chars().take(512).collect::<String>());
    json!({"index":row.index,"construct_id":row.construct_id,"name":row.name,"state":row.state,"status":row.status,"error":error,"error_abbreviated":row.error.as_ref().is_some_and(|error|error.chars().count()>512),"combination_id":row.combination_id,"declaration_revision":row.declaration_revision,"input_count":row.input_bindings.len(),"target_binding":row.target_binding,"plan_revision":row.plan_revision,"target_revision":row.target_revision,"final_product_revision":row.final_product_revision,"completed_steps":row.completed_steps,"planned_steps":row.planned_steps,"preparation_complete":row.preparation_complete})
}
impl CollageneServer {
    pub(crate) async fn combinatorial_release(
        &self,
        args: ReleaseQuery,
    ) -> Result<CallToolResult, ErrorData> {
        args.request.validate().map_err(bad)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;
            let fingerprint=collagene_core::material::lineage::workspace_fingerprint(&frame).map_err(bad)?;
            let context=release::revision(&frame,&args.request).map_err(bad)?;
            let mut session=release::ReleaseSession::begin(frame,args.request).map_err(bad)?;
            let mut progress=session.summary();
            if progress.context_revision!=context{return Err(bad("private release differs from its exact native workspace request"));}
            for index in 0..progress.total_rows{
                if progress.next_index!=index||progress.failed.is_some(){return Err(bad("private release cannot skip or retain a failed prefix"));}
                progress=session.advance(&context,index).map_err(bad)?;
            }
            if !progress.finished||progress.failed.is_some(){return Err(bad("no outputs were released because the private release did not complete"));}
            let commit=session.into_commit(&context).map_err(bad)?;
            commit.validate_binding(&context,&fingerprint).map_err(bad)?;
            let manifest=release::load(&commit.workspace,&commit.manifest_hash).map_err(bad)?;
            let text=bounded_json_limit(&manifest,release::MAX_MANIFEST_BYTES)?;
            let receipt=json!({"manifest_hash":commit.manifest_hash,"context_revision":commit.context_revision,"source_workspace_fingerprint":commit.source_workspace_fingerprint,"batch_hash":manifest.batch_hash,"recorded_engine":manifest.engine_version,"released_rows":manifest.rows,"state":"historical-output-lineage","verified":false});
            let native=workspace::save_workspace(&commit.workspace).map_err(bad)?;
            let committed=server.session.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"combinatorial-release.json","Freshly reconstructed selected rows released atomically into a new immutable native workspace. Manifest/run assertions are historical on reopening; experimental success is not asserted.")?;
            let mut value=committed.structured_content.ok_or_else(||bad("missing immutable release response"))?;
            value["authority"]=json!("Every selected row was reconstructed against its independent target before atomic native output release. Reopening the manifest retains historical lineage, not current scientific authority or experimental success.");
            output(value)
        }).await
    }
    pub(crate) async fn combinatorial_archive(
        &self,
        args: ExportQuery,
    ) -> Result<CallToolResult, ErrorData> {
        hash(&args.manifest_hash)?;
        args.options.validate().map_err(bad)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;
            let generated=archive::export(&frame,&args.manifest_hash,&args.options).map_err(bad)?;
            let mut response=json!({"workspace":args.workspace,"manifest_hash":generated.manifest_hash,"context_revision":generated.context_revision,"blake3":generated.blake3,"bytes":generated.bytes.len(),"files":generated.files,"warnings":generated.warnings,"mime_type":generated.mime,"filename":"combinatorial-outputs.zip","verified":false,"state":"portable-output-exchange"});
            let artifact=server.session.binary_artifact(generated.bytes,"application/zip","combinatorial-outputs.zip")?;
            response["artifact_id"]=json!(artifact);response["encoding"]=json!("base64");output(response)
        }).await
    }
    pub(crate) async fn combinatorial_release_list(
        &self,
        args: ListQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 64 {
            return Err(bad("release pages contain1..64 records"));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let all=release::list(&frame).map_err(bad)?;
            if args.offset>all.len(){return Err(bad("release page starts beyond the retained list"));}
            let total=all.len();let rows=all.into_iter().skip(args.offset).take(args.limit).collect::<Vec<_>>();
            output(json!({"workspace":args.workspace,"releases":rows,"total":total,"offset":args.offset,"next_offset":(args.offset+args.limit).min(total),"state":"historical","verified":false}))
        }).await
    }
    pub(crate) async fn combinatorial_release_open(
        &self,
        args: ReleaseRef,
    ) -> Result<CallToolResult, ErrorData> {
        hash(&args.manifest_hash)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let manifest=release::load(&frame,&args.manifest_hash).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&manifest,release::MAX_MANIFEST_BYTES)?,"application/json","combinatorial-release.json")?;
            output(json!({"workspace":args.workspace,"manifest_hash":args.manifest_hash,"context_revision":manifest.context_revision,"batch_hash":manifest.batch_hash,"source_workspace_fingerprint":manifest.source_workspace_fingerprint,"recorded_engine":manifest.engine_version,"rows":manifest.rows,"artifact_id":artifact,"state":"historical-output-lineage","verified":false,"authority":"Exact retained native output events and source lineage only. Opening does not restore fresh scientific authority."}))
        }).await
    }
    pub(crate) async fn combinatorial_evaluate(
        &self,
        args: BatchQuery,
    ) -> Result<CallToolResult, ErrorData> {
        args.matrix.bounds()?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let matrix=args.matrix.into_native(&frame)?;
            let context=batch::begin_revision(&frame,&matrix).map_err(bad)?;
            let install=batch::capture_revision(&frame,&context).map_err(bad)?;
            let mut session=batch::Session::begin(frame,matrix).map_err(bad)?;
            let mut progress=session.summary();
            if progress.context_revision!=context{return Err(bad("batch session differs from its exact native input context"));}
            let total=progress.total_rows;
            for index in 0..total{
                if progress.next_index!=index{return Err(bad("batch progression differs from the selected construct order"));}
                progress=session.advance(&context,index).map_err(bad)?;
            }
            if !progress.finished||progress.next_index!=total{return Err(bad("batch did not finish every explicitly selected row"));}
            let commit=session.into_capture(&context).map_err(bad)?;
            commit.validate_binding(&install).map_err(bad)?;
            let saved=batch::load(&commit.workspace,&commit.attachment_hash).map_err(bad)?;
            let text=bounded_json_limit(&saved,32*1024*1024)?;
            let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":saved.context_revision,"recorded_engine":saved.engine_version,"recorded_row_count":saved.rows.len(),"rows":saved.rows.iter().map(row_summary).collect::<Vec<_>>(),"state":"historical","verified":false});
            let native=workspace::save_workspace(&commit.workspace).map_err(bad)?;
            let committed=server.session.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"combinatorial-batch.json",
                "Historical source-event-bound batch evidence. Explicit replay reconstructs current row evidence; no experimental success or bulk product release is asserted.")?;
            // The shared output publisher also handles finalized construct
            // archives. This archive contains batch evidence only, so return
            // matching authority text in both structured and text payloads.
            let mut value=committed.structured_content.ok_or_else(||bad("missing committed batch response"))?;
            value["authority"]=json!("A freshly computed batch was retained as historical evidence in a new immutable workspace. No source head changed; no bulk product or experimental verification was released.");
            output(value)
        }).await
    }
    pub(crate) async fn combinatorial_list(
        &self,
        args: ListQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 64 {
            return Err(bad("saved batch pages contain1..64 records"));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let records=batch::list(&frame).map_err(bad)?;let count=records.len();
            if args.offset>count{return Err(bad("batch page starts beyond the saved record list"));}
            let page=records.into_iter().skip(args.offset).take(args.limit).collect::<Vec<_>>();
            output(json!({"workspace":args.workspace,"batches":page,"total":count,"offset":args.offset,"next_offset":(args.offset+args.limit).min(count),"state":"historical","verified":false}))
        }).await
    }
    pub(crate) async fn combinatorial_open(
        &self,
        args: OpenQuery,
    ) -> Result<CallToolResult, ErrorData> {
        hash(&args.attachment_hash)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let saved=batch::load(&frame,&args.attachment_hash).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&saved,32*1024*1024)?,"application/json","combinatorial-batch.json")?;
            output(json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"context_revision":saved.context_revision,"recorded_engine":saved.engine_version,"recorded_row_count":saved.rows.len(),"rows":saved.rows.iter().map(row_summary).collect::<Vec<_>>(),"artifact_id":artifact,"state":"historical","verified":false,"authority":"Saved outcomes are historical assertions. Opening does not restore fresh evaluation or authorize output release."}))
        }).await
    }
    pub(crate) async fn combinatorial_replay(
        &self,
        args: ReplayQuery,
    ) -> Result<CallToolResult, ErrorData> {
        hash(&args.attachment_hash)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let checked=batch::replay(&frame,&args.attachment_hash,args.row_index).map_err(bad)?;
            checked.validate().map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&checked,32*1024*1024)?,"application/json","combinatorial-batch-replay.json")?;
            output(json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"row_index":args.row_index,"recorded":row_summary(&checked.recorded),"current":row_summary(&checked.current),"current_engine":checked.current_engine,"artifact_id":artifact,"state":"fresh-computational-evidence","verified":false,"authority":"Recorded and freshly reconstructed row evidence are separate. Full findings remain in the artifact. No source, batch record or native workspace was changed; experimental success and bulk output release are not asserted."}))
        }).await
    }
    pub(crate) async fn combinatorial_matrix(
        &self,
        args: MatrixQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > design::MAX_PAGE {
            return Err(bad("matrix pages contain1..64 combinations"));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let declared=args.design.into_native(&frame)?;
            let page=design::matrix_page(&declared,args.first,args.limit).map_err(bad)?;
            output(json!({"workspace":args.workspace,"design_id":declared.id,"design_name":declared.name,"page":page,
                "state":"source-bound-declarations","verified":false,"authority":"Exact retained source events and declared native ownership were checked; actual material state/ranges are freshly checked when resolving a construct. Matrix enumeration does not check preparations, junctions, target agreement or experimental success."}))
        }).await
    }
    pub(crate) async fn combinatorial_pool_inspect(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;let declared=args.design.into_native(&frame)?;
            let resolved = design::resolve_pool(&declared, &args.construct, &frame).map_err(bad)?;
            let design::Recipe::GoldenGatePool { request, source_revision, .. } = &resolved.recipe else {
                return Err(bad("selected construct is not a Golden Gate reaction pool"));
            };
            let report = collagene_core::cloning::golden_gate::reaction::inspect(&resolved.inputs, &resolved.target, request).map_err(bad)?;
            let evidence = json!({"workspace":args.workspace,"design":declared,"construct":args.construct,
                "source_revision":source_revision,"construct_id":resolved.construct_id,
                "combination_id":resolved.combination_id,"declaration_revision":resolved.declaration_revision,
                "input_bindings":resolved.bindings,"target_binding":resolved.target_binding,
                "report":report,"verified":false,"state":"fresh-pool-discovery"});
            let artifact = server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","combinatorial-pool.json")?;
            output(json!({"workspace":args.workspace,"source_revision":source_revision,
                "construct_id":resolved.construct_id,"combination_id":resolved.combination_id,
                "declaration_revision":resolved.declaration_revision,"input_bindings":resolved.bindings,
                "target_binding":resolved.target_binding,"context_revision":report.context_revision,
                "assessments":report.assessments,"warnings":report.warnings,"artifact_id":artifact,
                "state":"fresh-pool-discovery","verified":false,
                "authority":"Finite native reaction discovery only. Choose an explicit index and both candidate keys with this source_revision. Target correspondence never chooses a candidate. Batch evaluation and output release reconstruct the selected path afresh; imported reports do not restore authority."}))
        }).await
    }
    pub(crate) async fn combinatorial_multisite_inspect(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;let declared=args.design.into_native(&frame)?;
            let report = design::inspect_multisite(&declared, &args.construct, &frame).map_err(bad)?;
            report.validate_binding(&design::inspect_multisite_revision(&declared, &args.construct, &frame).map_err(bad)?).map_err(bad)?;
            let candidates = report.batch.reference.candidates.iter().take(25).map(|candidate| json!({
                "candidate_id":candidate.id,"sites":candidate.sites,"pairs":candidate.pairs,
                "components":candidate.products.iter().enumerate().map(|(i,p)|json!({"component_id":i,"length":p.molecule.sequence.len(),"topology":p.molecule.topology,"reference_revision":p.molecule.revision()})).collect::<Vec<_>>()
            })).collect::<Vec<_>>();
            let artifact = server.session.text_artifact(bounded_json_limit(&json!({
                "workspace":args.workspace,"design":declared,"construct":args.construct,"inspection":report,
                "state":"fresh-reference-network-discovery","verified":false
            }),32*1024*1024)?,"application/json","combinatorial-multisite-reference.json")?;
            output(json!({"workspace":args.workspace,"context_revision":report.context_revision,"source_revision":report.source_revision,
                "batch_revision":report.batch.batch_revision,"complete":report.batch.reference.complete,
                "candidate_count":report.batch.reference.candidates.len(),"candidates":candidates,
                "candidates_truncated":report.batch.reference.candidates.len()>25,
                "warnings":report.batch.warnings,"unsupported":report.batch.reference.unsupported,"limits":report.batch.reference.limits,
                "artifact_id":artifact,"state":"fresh-reference-network-discovery","verified":false,
                "authority":"Independent target intent precedes reconstruction. Explicitly choose a candidate/component with this source_revision, then evaluate/release through the existing library workflow. All reciprocal components remain in typed evidence. Reference agreement does not establish physical chemistry."}))
        }).await
    }
    pub(crate) async fn combinatorial_inspect(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.session.restore_workspace(&args.workspace)?;let declared=args.design.into_native(&frame)?;
            let resolved=design::resolve(&declared,&args.construct,&frame).map_err(bad)?;
            let mut summary=summaries(&resolved);
            let intent=match &resolved.intent{Some(design::IndependentIntent::Regions(report))=>json!({"OriginalRegions":report}),Some(design::IndependentIntent::Junctions(report))=>json!({"OriginalJunctions":report}),Some(design::IndependentIntent::Substitutions(report))=>json!({"OriginalSubstitutions":report}),Some(design::IndependentIntent::Edits(report))=>json!({"OriginalEdits":report}),None=>Value::Null};
            let recipe=match &resolved.recipe{design::Recipe::Overlap{recipe,finishing}=>json!({"Overlap":{"recipe":recipe,"finishing":finishing}}),design::Recipe::Tailing{recipe}=>json!({"Tailing":{"recipe":recipe}}),design::Recipe::Ligation{recipe,cleavage_policy}=>json!({"Ligation":{"recipe":recipe,"cleavage_policy":cleavage_policy}}),design::Recipe::GoldenGatePool{request,selection,source_revision}=>json!({"GoldenGatePool":{"request":request,"selection":selection,"source_revision":source_revision}}),design::Recipe::MultiSiteReference{request,selection,source_revision,sources}=>json!({"MultiSiteReference":{"request":request,"selection":selection,"source_revision":source_revision,"sources":sources}})};
            let evidence=json!({"workspace":args.workspace,"design":declared,"construct":args.construct,"construct_id":resolved.construct_id,"combination_id":resolved.combination_id,"declaration_revision":resolved.declaration_revision,"input_bindings":resolved.bindings,"inputs":resolved.inputs,"independent_target":resolved.target,"target_binding":resolved.target_binding,"independent_intent":intent,"recipe":recipe,"verified":false});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","combinatorial-construct.json")?;
            summary["workspace"]=serde_json::to_value(args.workspace).map_err(bad)?;summary["artifact_id"]=json!(artifact);
            summary["authority"]=json!("Independent target and explicit preparation recipe only. Use fresh batch evaluation to reconstruct operations and check target agreement; inspection does not establish experimental success.");
            output(summary)
        }).await
    }
}
