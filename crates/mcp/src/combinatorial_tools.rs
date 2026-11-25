//! Exact native-workspace source alternatives and fresh selected batch checks.
//! All computation runs outside the immutable session store lock.
use super::*;
use collagene_core::{combinatorial_batch as batch, combinatorial_design as design};

fn page_size() -> usize {
    20
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MatrixQuery {
    pub workspace: WorkspaceRef,
    pub design: design::Design,
    #[serde(default)]
    pub first: u128,
    #[serde(default = "page_size")]
    pub limit: usize,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectQuery {
    pub workspace: WorkspaceRef,
    pub design: design::Design,
    pub construct: design::SelectedConstruct,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BatchQuery {
    pub workspace: WorkspaceRef,
    pub matrix: batch::Request,
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
// Matrix rows are declarations, but their claimed source ownership and literal
// original intervals can still be checked without running preparation science.
fn attest_sources(
    frame: &workspace::Workspace,
    declared: &design::Design,
) -> Result<(), ErrorData> {
    design::validate_design(declared).map_err(bad)?;
    let mut checked = BTreeSet::new();
    for bin in &declared.bins {
        for choice in &bin.choices {
            let b = &choice.source;
            if !frame.documents.iter().any(|doc| doc.id == b.document_id) {
                return Err(bad("an alternative source document is not retained"));
            }
            let event = frame
                .history
                .iter()
                .find(|event| event.id == b.event_id)
                .ok_or_else(|| bad("an alternative source event is not retained"))?;
            if event.document_id != b.document_id || event.revision != b.revision {
                return Err(bad("alternative source document/event/revision mismatch"));
            }
            let molecule = frame
                .revisions
                .get(&b.revision)
                .ok_or_else(|| bad("an alternative source revision is not retained"))?;
            if checked.insert(&b.revision) && molecule.revision() != b.revision {
                return Err(bad(
                    "alternative source content differs from the retained revision",
                ));
            }
            let n = molecule.sequence.len();
            if molecule.sequence.alphabet() != collagene_core::Alphabet::Dna
                || choice.start >= n
                || choice.end <= choice.start
                || choice.end - choice.start > n
                || (molecule.topology == collagene_core::Topology::Linear && choice.end > n)
            {
                return Err(bad(
                    "an alternative contribution is outside its original DNA source or exceeds one turn",
                ));
            }
        }
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
    pub(super) async fn combinatorial_evaluate(
        &self,
        args: BatchQuery,
    ) -> Result<CallToolResult, ErrorData> {
        args.matrix.validate().map_err(bad)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;
            let context=batch::begin_revision(&frame,&args.matrix).map_err(bad)?;
            let install=batch::capture_revision(&frame,&context).map_err(bad)?;
            let mut session=batch::Session::begin(frame,args.matrix).map_err(bad)?;
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
            let committed=server.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"combinatorial-batch.json",
                "Historical source-event-bound batch evidence. Explicit replay reconstructs current row evidence; no experimental success or bulk product release is asserted.")?;
            // The shared output publisher also handles finalized construct
            // archives. This archive contains batch evidence only, so return
            // matching authority text in both structured and text payloads.
            let mut value=committed.structured_content.ok_or_else(||bad("missing committed batch response"))?;
            value["authority"]=json!("A freshly computed batch was retained as historical evidence in a new immutable workspace. No source head changed; no bulk product or experimental verification was released.");
            output(value)
        }).await
    }
    pub(super) async fn combinatorial_list(
        &self,
        args: ListQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 64 {
            return Err(bad("saved batch pages contain1..64 records"));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;let records=batch::list(&frame).map_err(bad)?;let count=records.len();
            if args.offset>count{return Err(bad("batch page starts beyond the saved record list"));}
            let page=records.into_iter().skip(args.offset).take(args.limit).collect::<Vec<_>>();
            output(json!({"workspace":args.workspace,"batches":page,"total":count,"offset":args.offset,"next_offset":(args.offset+args.limit).min(count),"state":"historical","verified":false}))
        }).await
    }
    pub(super) async fn combinatorial_open(
        &self,
        args: OpenQuery,
    ) -> Result<CallToolResult, ErrorData> {
        hash(&args.attachment_hash)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;let saved=batch::load(&frame,&args.attachment_hash).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&saved,32*1024*1024)?,"application/json","combinatorial-batch.json")?;
            output(json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"context_revision":saved.context_revision,"recorded_engine":saved.engine_version,"recorded_row_count":saved.rows.len(),"rows":saved.rows.iter().map(row_summary).collect::<Vec<_>>(),"artifact_id":artifact,"state":"historical","verified":false,"authority":"Saved outcomes are historical assertions. Opening does not restore fresh evaluation or authorize output release."}))
        }).await
    }
    pub(super) async fn combinatorial_replay(
        &self,
        args: ReplayQuery,
    ) -> Result<CallToolResult, ErrorData> {
        hash(&args.attachment_hash)?;
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;let checked=batch::replay(&frame,&args.attachment_hash,args.row_index).map_err(bad)?;
            checked.validate().map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&checked,32*1024*1024)?,"application/json","combinatorial-batch-replay.json")?;
            output(json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"row_index":args.row_index,"recorded":row_summary(&checked.recorded),"current":row_summary(&checked.current),"current_engine":checked.current_engine,"artifact_id":artifact,"state":"fresh-computational-evidence","verified":false,"authority":"Recorded and freshly reconstructed row evidence are separate. Full findings remain in the artifact. No source, batch record or native workspace was changed; experimental success and bulk output release are not asserted."}))
        }).await
    }
    pub(super) async fn combinatorial_matrix(
        &self,
        args: MatrixQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > design::MAX_PAGE {
            return Err(bad("matrix pages contain1..64 combinations"));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;attest_sources(&frame,&args.design)?;
            let page=design::matrix_page(&args.design,args.first,args.limit).map_err(bad)?;
            output(json!({"workspace":args.workspace,"design_id":args.design.id,"design_name":args.design.name,"page":page,
                "state":"source-bound-declarations","verified":false,"authority":"Exact retained source events and original contribution intervals were resolved. Matrix enumeration does not check preparations, junctions, target agreement or experimental success."}))
        }).await
    }
    pub(super) async fn combinatorial_inspect(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;
            let resolved=design::resolve(&args.design,&args.construct,&frame).map_err(bad)?;
            let mut summary=summaries(&resolved);
            let intent=match &resolved.intent{Some(design::IndependentIntent::Regions(report))=>json!({"OriginalRegions":report}),Some(design::IndependentIntent::Junctions(report))=>json!({"OriginalJunctions":report}),None=>Value::Null};
            let recipe=match &resolved.recipe{design::Recipe::Overlap(recipe)=>json!({"Overlap":recipe}),design::Recipe::Ligation(recipe)=>json!({"Ligation":recipe})};
            let evidence=json!({"workspace":args.workspace,"design":args.design,"construct":args.construct,"construct_id":resolved.construct_id,"combination_id":resolved.combination_id,"declaration_revision":resolved.declaration_revision,"input_bindings":resolved.bindings,"inputs":resolved.inputs,"independent_target":resolved.target,"target_binding":resolved.target_binding,"independent_intent":intent,"recipe":recipe,"verified":false});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","combinatorial-construct.json")?;
            summary["workspace"]=serde_json::to_value(args.workspace).map_err(bad)?;summary["artifact_id"]=json!(artifact);
            summary["authority"]=json!("Independent target and explicit preparation recipe only. Use fresh batch evaluation to reconstruct operations and check target agreement; inspection does not establish experimental success.");
            output(summary)
        }).await
    }
}
