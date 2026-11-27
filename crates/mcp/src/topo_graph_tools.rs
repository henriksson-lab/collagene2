//! Compact native graph selections. Supplied reports/capsules are never session authority.
use super::*;
use collagene_core::{
    run_material,
    topo::{self, physical},
    topo_graph_history as graph,
};
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum Material {
    Prepared {
        declaration_id: String,
        assumption: String,
        strands: Vec<topo::PreparedStrand>,
        pairing: Vec<physical::PairSpan>,
    },
    ExactPcr {
        material_id: String,
        template_role: String,
        forward: collagene_core::analysis::Primer,
        reverse: collagene_core::analysis::Primer,
        tail: Option<collagene_core::tailing::TailingBase>,
    },
    SavedGraph {
        material_id: String,
        selection: graph::Selection,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Recipe {
    pub schema_version: u16,
    pub materials: Vec<Material>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Preparation {
    pub recipe: Recipe,
    pub transfers: Vec<physical::Transfer>,
    pub engine_limits: physical::Limits,
    pub projection_options: physical::projection::Options,
}
impl Preparation {
    pub(super) fn selections(&self) -> Result<Vec<graph::Selection>, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        if self.recipe.materials.is_empty()
            || self.recipe.materials.len() > 16
            || self.transfers.len() > 16
        {
            return Err(bad("bounded current TOPO materials/transfers required"));
        }
        Ok(self
            .recipe
            .materials
            .iter()
            .filter_map(|m| {
                if let Material::SavedGraph { selection, .. } = m {
                    Some(selection.clone())
                } else {
                    None
                }
            })
            .collect())
    }
    pub(super) fn from_native(request: &topo::Request) -> Self {
        let materials = request
            .recipe
            .materials
            .iter()
            .map(|m| match m {
                topo::MaterialRecipe::Prepared {
                    declaration_id,
                    assumption,
                    strands,
                    pairing,
                } => Material::Prepared {
                    declaration_id: declaration_id.clone(),
                    assumption: assumption.clone(),
                    strands: strands.clone(),
                    pairing: pairing.clone(),
                },
                topo::MaterialRecipe::ExactPcr {
                    material_id,
                    template_role,
                    forward,
                    reverse,
                    tail,
                } => Material::ExactPcr {
                    material_id: material_id.clone(),
                    template_role: template_role.clone(),
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    tail: *tail,
                },
                topo::MaterialRecipe::SavedGraph {
                    material_id,
                    origin,
                } => Material::SavedGraph {
                    material_id: material_id.clone(),
                    selection: origin.selection.clone(),
                },
            })
            .collect();
        Self {
            recipe: Recipe {
                schema_version: request.recipe.schema_version,
                materials,
            },
            transfers: request.transfers.clone(),
            engine_limits: request.engine_limits.clone(),
            projection_options: request.projection_options.clone(),
        }
    }
    fn native(
        &self,
        graphs: &mut impl Iterator<Item = graph::Capsule>,
    ) -> Result<topo::Request, ErrorData> {
        let mut materials = Vec::with_capacity(self.recipe.materials.len());
        for material in &self.recipe.materials {
            materials.push(match material {
                Material::Prepared {
                    declaration_id,
                    assumption,
                    strands,
                    pairing,
                } => topo::MaterialRecipe::Prepared {
                    declaration_id: declaration_id.clone(),
                    assumption: assumption.clone(),
                    strands: strands.clone(),
                    pairing: pairing.clone(),
                },
                Material::ExactPcr {
                    material_id,
                    template_role,
                    forward,
                    reverse,
                    tail,
                } => topo::MaterialRecipe::ExactPcr {
                    material_id: material_id.clone(),
                    template_role: template_role.clone(),
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    tail: *tail,
                },
                Material::SavedGraph {
                    material_id,
                    selection,
                } => {
                    let origin = graphs
                        .next()
                        .ok_or_else(|| bad("native graph hydration count differs"))?;
                    if &origin.selection != selection {
                        return Err(bad("native graph hydration selection differs"));
                    }
                    topo::MaterialRecipe::SavedGraph {
                        material_id: material_id.clone(),
                        origin,
                    }
                }
            });
        }
        Ok(topo::Request {
            recipe: topo::GenerationRecipe {
                schema_version: self.recipe.schema_version,
                materials,
            },
            transfers: self.transfers.clone(),
            engine_limits: self.engine_limits.clone(),
            projection_options: self.projection_options.clone(),
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Sealing {
    pub preparation: Preparation,
    pub sealing: physical::nick_sealing::Request,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Processed {
    pub treatment: Sealing,
    pub selection: physical::processed_projection::Selection,
    pub options: physical::projection::Options,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Actions {
    pub preparation: Preparation,
    pub processing: physical::processing::Request,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ActionProduct {
    pub preparation: Preparation,
    pub processing: physical::processing::Request,
    pub selection: physical::processing_projection::Selection,
    pub options: physical::projection::Options,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum PreparationSource {
    Original {
        vector: collagene_core::topo_preparation::PreparedVectorDeclaration,
        insert: collagene_core::topo_preparation::ExactPcrDeclaration,
    },
    SavedVector {
        material_id: String,
        selection: graph::Selection,
        vector_strands: [usize; 2],
        insert: collagene_core::topo_preparation::ExactPcrDeclaration,
    },
    ResumeGraph {
        material_id: String,
        selection: graph::Selection,
        vector_strands: [usize; 2],
        insert_strands: [usize; 2],
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Declaration {
    pub mode: collagene_core::topo_preparation::PreparationMode,
    pub source: PreparationSource,
    pub engine_limits: physical::Limits,
    pub projection_options: physical::projection::Options,
}
impl Declaration {
    fn selection(&self) -> Vec<graph::Selection> {
        match &self.source {
            PreparationSource::Original { .. } => vec![],
            PreparationSource::SavedVector { selection, .. }
            | PreparationSource::ResumeGraph { selection, .. } => vec![selection.clone()],
        }
    }
    fn native(
        &self,
        origin: Option<graph::Capsule>,
    ) -> Result<collagene_core::topo_preparation::PreparationRequest, ErrorData> {
        use collagene_core::topo_preparation::{PreparationRequest, PreparationSource as Native};
        let source = match &self.source {
            PreparationSource::Original { vector, insert } => Native::Original {
                vector: vector.clone(),
                insert: insert.clone(),
            },
            PreparationSource::SavedVector {
                material_id,
                selection,
                vector_strands,
                insert,
            } => {
                let origin =
                    origin.ok_or_else(|| bad("saved vector graph selection was not hydrated"))?;
                if origin.selection != *selection {
                    return Err(bad("saved vector graph selection changed"));
                }
                Native::SavedVector {
                    material_id: material_id.clone(),
                    origin,
                    vector_strands: *vector_strands,
                    insert: insert.clone(),
                }
            }
            PreparationSource::ResumeGraph {
                material_id,
                selection,
                vector_strands,
                insert_strands,
            } => {
                let origin =
                    origin.ok_or_else(|| bad("partial graph selection was not hydrated"))?;
                if origin.selection != *selection {
                    return Err(bad("partial graph selection changed"));
                }
                Native::ResumeGraph {
                    material_id: material_id.clone(),
                    origin,
                    vector_strands: *vector_strands,
                    insert_strands: *insert_strands,
                }
            }
        };
        Ok(PreparationRequest {
            mode: self.mode,
            source,
            engine_limits: self.engine_limits.clone(),
            projection_options: self.projection_options.clone(),
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoundTemplate {
    pub role: String,
    pub selection: run_material::MaterialSelection,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SourceBinding {
    pub role: String,
    pub event: collagene_core::run_lineage::EventBinding,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin: Option<run_material::MaterialSelection>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CaptureRequest {
    pub name: String,
    pub bindings: Vec<SourceBinding>,
    pub request: Preparation,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Query {
    pub workspace: WorkspaceRef,
    pub task: Task,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum Task {
    List { offset: usize, limit: usize },
    Inspect { selection: graph::Selection },
    Capture { request: CaptureRequest },
}
impl CollageneServer {
    /// One shared documentary construction phase for every current native root.
    /// Enclosing native task alone freshly executes the retained+new transcript.
    pub(super) fn hydrate_topo(
        &self,
        request: &Preparation,
        workspace_ref: Option<&WorkspaceRef>,
        materials: &[run_material::MaterialSelection],
    ) -> Result<(topo::Request, Vec<run_material::MaterialOrigin>), ErrorData> {
        let selected = request.selections()?;
        let roots = if selected.is_empty() && materials.is_empty() {
            run_material::DocumentaryOrigins {
                materials: vec![],
                graphs: vec![],
            }
        } else {
            let frame = self.restore_workspace(
                workspace_ref.ok_or_else(|| bad("saved TOPO sources require origin_workspace"))?,
            )?;
            run_material::documentary_origins(&frame, materials, &selected).map_err(bad)?
        };
        let mut graphs = roots.graphs.into_iter();
        let native = request.native(&mut graphs)?;
        if graphs.next().is_some() {
            return Err(bad("unused native graph hydration"));
        }
        Ok((native, roots.materials))
    }
    pub(super) fn hydrate_topo_declaration(
        &self,
        query: &Declaration,
        workspace_ref: Option<&WorkspaceRef>,
        template: Option<&BoundTemplate>,
    ) -> Result<
        (
            collagene_core::topo_preparation::PreparationRequest,
            Option<collagene_core::topo_run_preparation::BoundTemplate>,
        ),
        ErrorData,
    > {
        bounded_json_limit(query, 2 * 1024 * 1024)?;
        if template.is_some_and(|t| t.role.trim().is_empty() || t.role.len() > 256) {
            return Err(bad("bounded saved PCR-template role required"));
        }
        let graphs = query.selection();
        let materials = template
            .map(|t| vec![t.selection.clone()])
            .unwrap_or_default();
        let mut origins = if graphs.is_empty() && materials.is_empty() {
            run_material::DocumentaryOrigins {
                graphs: vec![],
                materials: vec![],
            }
        } else {
            let frame = self.restore_workspace(
                workspace_ref.ok_or_else(|| bad("saved preparation requires origin_workspace"))?,
            )?;
            run_material::documentary_origins(&frame, &materials, &graphs).map_err(bad)?
        };
        let declaration = query.native(origins.graphs.pop())?;
        let template = if let Some(template) = template {
            Some(collagene_core::topo_run_preparation::BoundTemplate {
                role: template.role.clone(),
                origin: origins
                    .materials
                    .pop()
                    .ok_or_else(|| bad("saved PCR template was not hydrated"))?,
            })
        } else {
            None
        };
        Ok((declaration, template))
    }
    pub(super) async fn topo_graph_dispatch(
        &self,
        args: Query,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;
            match args.task {
                Task::List{offset,limit}=>{
                    if offset>graph::MAX_RECORDS||limit==0||limit>200{return Err(bad("bounded native graph page required"));}
                    let expected=graph::list_revision(&frame).map_err(bad)?;let listed=graph::list(&frame).map_err(bad)?;listed.validate_binding(&expected).map_err(bad)?;
                    let total=listed.entries.len();let end=offset.saturating_add(limit).min(total);let rows=listed.entries.get(offset.min(total)..end).unwrap_or(&[]);
                    let artifact=server.text_artifact(bounded_json_limit(&listed,32*1024*1024)?,"application/json","native-topo-graph-catalog.json")?;
                    let entries=rows.iter().map(|row|json!({"selection":row.selection,"name":row.name,"source_event_count":row.source_events.len(),"source_bases":row.source_bases,"chains":row.chains,"components":row.components,"remaining_activations":row.remaining_activations,"released_activations":row.released_activations,"incomplete":row.incomplete,"finding_count":row.findings.len(),"findings":row.findings.iter().take(8).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"findings_truncated":row.findings.len()>8||row.findings.iter().any(|s|s.chars().count()>1024),"assumption_count":row.assumptions.len()})).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"context_revision":expected,"source_workspace_fingerprint":listed.source_workspace_fingerprint,"total":total,"entries":entries,"artifact_id":artifact,"next_offset":if end<total{Some(end)}else{None},"authority":"Recorded whole-graph summaries only; full source pins/findings/assumptions in artifact. Listing does not replay material, select a component, establish activation manufacture or experimental success."}))
                }
                Task::Inspect{selection}=>{
                    let expected=graph::inspect_revision(&frame,&selection).map_err(bad)?;let inspection=graph::inspect(&frame,&selection).map_err(bad)?;inspection.validate_binding(&expected).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&inspection,32*1024*1024)?,"application/json","native-topo-graph.json")?;
                    output(json!({"workspace":args.workspace,"selection":selection,"context_revision":expected,"source_workspace_fingerprint":inspection.source_workspace_fingerprint,"summary":inspection.summary,"artifact_id":artifact,"verified":false,"authority":"Fresh native execution of declared sources and exact elementary transcript. Entire graph/calls/byproducts and inherited assumptions retained; no ordinary-duplex restoration, physical isolation or experimental success."}))
                }
                Task::Capture{request}=>{
                    bounded_json_limit(&request,2*1024*1024)?;
                    if request.bindings.len()>32{return Err(bad("native graph source binding bound"));}
                    let selections=request.bindings.iter().filter_map(|b|b.origin.clone()).collect::<Vec<_>>();
                    for binding in &request.bindings {if binding.origin.as_ref().is_some_and(|o|o.event()!=&binding.event){return Err(bad("graph PCR-template selection names another exact event"));}}
                    let (prepared,origins)=server.hydrate_topo(&request.request,Some(&args.workspace),&selections)?;
                    let mut origins=origins.into_iter();let mut bindings=Vec::with_capacity(request.bindings.len());
                    for binding in &request.bindings {bindings.push(graph::SourceBinding{role:binding.role.clone(),event:binding.event.clone(),origin:if binding.origin.is_some(){Some(origins.next().ok_or_else(||bad("native template hydration differs"))?)}else{None}});}
                    let native=graph::CaptureRequest{name:request.name,bindings,request:prepared};
                    let expected=graph::revision(&frame,&native).map_err(bad)?;
                    let capture=graph::prepare_capture(&frame,&native).map_err(bad)?;capture.validate_binding(&expected).map_err(bad)?;
                    let text=bounded_json_limit(&capture.receipt,4*1024*1024)?;
                    let receipt=json!({"record_hash":capture.receipt.record_hash,"selection":capture.receipt.selection,"context_revision":expected,"summary":capture.receipt.summary,"state":"native-declared-elementary-graph","verified":false});
                    let archive=workspace::save_workspace(&capture.workspace).map_err(bad)?;
                    server.commit_workspace_outputs_named(capture.workspace,Some((text,receipt)),Some(archive),"topo-graph-capture.json","Target-free declared elementary graph. Biological heads unchanged; this is neither observed activated stock, ordinary duplex material nor an accepted cloning result.")
                }
            }
        }).await
    }
}
