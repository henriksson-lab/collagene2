pub mod assembly_library_history_tools;
pub mod combinatorial_tools;
pub mod domestication_tools;
pub mod fusion_design_tools;
pub mod intents;
pub mod multisite_recombination_tools;
pub mod published_ligation_tools;
pub mod recombination_tools;
pub mod topo_actions_tools;
pub mod topo_graph_tools;

use crate::server::CollageneServer;
use crate::session::Artifact;
use crate::session::ArtifactData;
use crate::session::MAX_ITEMS;
use crate::session::StoredPlan;
use crate::session::references::DocumentRef;
use crate::session::references::EvidenceRef;
use crate::session::references::MaterialOriginWorkspace;
use crate::session::references::WorkspaceRef;
use crate::tools::assembly::intents::BoundInput;

use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::bounded_json_limit;

use crate::server::response::output;

use collagene_core::documents::transfer;
use collagene_core::planning::plan::PlanInput;
use collagene_model::Molecule;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use rmcp::tool;
use rmcp::tool_router;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;
use serde_json::Value;
use serde_json::json;
use std::collections::BTreeSet;

use collagene_core::verification;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoSourceRef {
    pub role: String,
    pub document: DocumentRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoSimulateQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoOrientationQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    /// Exact prepared-strand indexes; not source-document indexes or a kit guess.
    pub vector: [usize; 2],
    pub insert: [usize; 2],
    pub directional: bool,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoProjectQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    pub selection: collagene_core::cloning::topo::Selection,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoTargetQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    pub selection: collagene_core::cloning::topo::Selection,
    pub target: DocumentRef,
    /// Every allowance is explicit; the target never chooses a component/seed.
    pub policy: collagene_core::cloning::topo::AcceptancePolicy,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
/// Saved template context is attested against an exact session-owned archive;
/// the task can never supply an alleged reconstructed duplex.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum SavedTopoTask {
    Orientations {
        request: topo_graph_tools::Preparation,
        vector: [usize; 2],
        insert: [usize; 2],
        directional: bool,
    },
    Project {
        request: topo_graph_tools::Preparation,
        selection: collagene_core::cloning::topo::Selection,
    },
    Check {
        request: topo_graph_tools::Preparation,
        selection: collagene_core::cloning::topo::Selection,
        target: DocumentRef,
        policy: collagene_core::cloning::topo::AcceptancePolicy,
    },
    Simulate {
        request: topo_graph_tools::Preparation,
    },
    ProcessActions {
        request: topo_graph_tools::ActionProduct,
    },
}
impl SavedTopoTask {
    pub(crate) fn preparation(&self) -> &topo_graph_tools::Preparation {
        match self {
            Self::Orientations { request, .. }
            | Self::Project { request, .. }
            | Self::Check { request, .. }
            | Self::Simulate { request } => request,
            Self::ProcessActions { request } => &request.preparation,
        }
    }
    pub(crate) fn native(
        &self,
        target: Option<Molecule>,
        preparation: collagene_core::cloning::topo::Request,
    ) -> Result<collagene_core::cloning::topo::run_preparation::Task, ErrorData> {
        use collagene_core::cloning::topo::run_preparation::Task;
        Ok(match self {
            Self::Orientations {
                request: _,
                vector,
                insert,
                directional,
            } => Task::Orientations {
                request: preparation,
                vector: *vector,
                insert: *insert,
                directional: *directional,
            },
            Self::Project {
                request: _,
                selection,
            } => Task::Project {
                request: preparation,
                selection: selection.clone(),
            },
            Self::Check {
                request: _,
                selection,
                policy,
                ..
            } => Task::Check {
                request: preparation,
                selection: selection.clone(),
                target: target.ok_or_else(|| bad("Independent target reference is required"))?,
                policy: policy.clone(),
            },
            Self::Simulate { request: _ } => Task::Simulate {
                request: preparation,
            },
            Self::ProcessActions { request } => Task::ProcessActions {
                request: collagene_core::cloning::topo::actions::Request {
                    preparation,
                    selection: request.selection.clone(),
                    options: request.options.clone(),
                },
            },
        })
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SavedTopoInspectionQuery {
    pub inputs: Vec<BoundInput>,
    #[serde(deserialize_with = "Option::deserialize")]
    pub template: Option<topo_graph_tools::BoundTemplate>,
    pub origin_workspace: WorkspaceRef,
    pub task: SavedTopoTask,
}
#[derive(Serialize)]
struct TopoEvidence<'a, T: Serialize> {
    origin_workspace: &'a Option<WorkspaceRef>,
    source_bindings: &'a [TopoSourceRef],
    request: &'a collagene_core::cloning::topo::Request,
    target_requirement: Option<&'a DocumentRef>,
    report: &'a T,
    binding_scope: &'static str,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentQuery {
    pub alignment: EvidenceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct HomologyRecipePreviewQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub recipe: collagene_core::cloning::assembly::recipe::Recipe,
    /// Explicit native strand resection/fill/sealing policy; not inferred from a vendor label.
    pub policy: collagene_core::cloning::overlap::assembly::Options,
    #[serde(default)]
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnealingRecipePreviewQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub recipe: collagene_core::cloning::assembly::recipe::Recipe,
    /// Explicit native 3′ resection/annealing policy, with no fill or sealing; not inferred from a vendor label.
    pub policy: collagene_core::cloning::overlap::annealing::Options,
    #[serde(default)]
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapAssemblyProposal {
    /// Exact immutable source documents, identified by recipe-local input IDs.
    pub inputs: Vec<BoundInput>,
    /// Independently supplied intended DNA; never inferred from prepared fragments.
    pub target: DocumentRef,
    /// Ordered preparations, orientation, explicit homology lengths and method profile.
    pub recipe: collagene_core::cloning::assembly::recipe::Recipe,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PreparedGoldenGatePoolQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub request: collagene_core::cloning::golden_gate::prepared::Request,
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PreparedGoldenGateSelectionQuery {
    pub artifact_id: String,
    pub selection: collagene_core::cloning::golden_gate::prepared::Selection,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGateReactionQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub request: collagene_core::cloning::golden_gate::reaction::Request,
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct GoldenGateReactionEvidence {
    kind: String,
    inputs: Vec<BoundInput>,
    target: DocumentRef,
    report: collagene_core::cloning::golden_gate::reaction::Report,
    origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PreparedGoldenGateEvidence {
    kind: String,
    inputs: Vec<BoundInput>,
    target: DocumentRef,
    report: collagene_core::cloning::golden_gate::prepared::Report,
    origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LigationAssemblyProposal {
    /// Original immutable DNA sources; preparations never replace these inputs.
    pub inputs: Vec<BoundInput>,
    /// Independent desired reference, not a sequence inferred from a ligated product.
    pub target: DocumentRef,
    /// Explicit ordered physical preparations, end treatments and ligation policy.
    pub recipe: collagene_core::cloning::ligation::recipe::Recipe,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TailingAssemblyProposal {
    /// Original immutable DNA sources; preparations never replace these inputs.
    pub inputs: Vec<BoundInput>,
    /// Independent desired reference, not a sequence inferred from a ligated product.
    pub target: DocumentRef,
    /// Explicit TA/GC profile, original-source preparations, per-row keep/add-tail actions, terminal treatments and ligation policy.
    pub recipe: collagene_core::cloning::ligation::tailing_recipe::Recipe,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssemblyJunctionQuery {
    /// Exact immutable saved candidate; its separately bound target is reused.
    pub plan_revision: String,
    /// Explicit bounded sequence flanks and genetic-code inspection policy.
    pub options: collagene_core::cloning::assembly::junctions::Options,
}
pub(crate) fn topo_summary(report: &collagene_core::cloning::topo::physical::Report) -> Value {
    json!({"engine":report.engine,"model":report.model,"physical_batch_revision":report.batch_revision,
        "work_complete":report.work_complete,"all_requested_steps_applied":report.all_requested_steps_applied,
        "unsupported_or_unknown":report.unsupported_or_unknown,"fully_specified_dna":report.fully_specified_dna,
        "step_count":report.steps.len(),"component_count":report.product.components.len(),"chain_count":report.product.chains.len(),
        "nick_count":report.product.nicks.len(),"unpaired_run_count":report.product.unpaired.len(),"released_activation_count":report.product.released_activations.len(),
        "components":report.product.components.iter().take(64).enumerate().map(|(i,c)|json!({"component":i,"chains":c.chains,"paired_nucleotides":c.paired_nucleotides,"unpaired_nucleotides":c.unpaired_nucleotides,"covalently_closed":c.covalently_closed})).collect::<Vec<_>>(),
        "components_truncated":report.product.components.len()>64,
        "chain_lengths":report.product.chains.iter().take(64).map(|c|c.sequence.len()).collect::<Vec<_>>(),
        "chains_truncated":report.product.chains.len()>64})
}
pub(crate) fn topo_annotation_summary(
    evidence: Option<&collagene_core::cloning::topo::annotations::AnnotationEvidence>,
) -> Value {
    match evidence {
        None => {
            json!({"state":"withheld_without_paired_reference","authority":"No annotation coordinates or chemical state inferred."})
        }
        Some(evidence) => {
            json!({"state":"source_correspondence_only","features":evidence.features.len(),"sources":evidence.sources.iter().map(|source|json!({"role":source.role,"revision":source.source_revision,"projection_supported":source.projection_supported,"projected":source.projected,"dropped_from_selected_reference":source.dropped,"unmapped_original_features":source.unmapped.len(),"warnings":source.warnings.iter().take(8).collect::<Vec<_>>() })).collect::<Vec<_>>(),"authority":"Review disrupted/PCR-derived annotations; retained source features are not chemical ancestry or restored physical authority."})
        }
    }
}
pub(crate) fn topo_projection_summary(product: &collagene_core::cloning::topo::Product) -> Value {
    use collagene_core::cloning::topo::physical::projection::Outcome;
    match &product.projection.outcome {
        Outcome::Proven(view) => {
            json!({"state":"paired_contour","topology":view.topology,"length":view.top.len(),
            "nick_count":view.nicks.len(),"fully_known_terminal_chemistry":view.fully_known_terminal_chemistry,
            "all_requested_steps_applied":view.all_requested_steps_applied,
            "reference_available":product.reference.is_some(),"source_map_count":view.source_maps.len()})
        }
        Outcome::Withheld(reason) => {
            json!({"state":"withheld","reason":reason,"reference_available":false})
        }
    }
}
impl CollageneServer {
    pub(crate) fn topo_sources(
        &self,
        sources: &[TopoSourceRef],
        request: &collagene_core::cloning::topo::Request,
        target: Option<&DocumentRef>,
    ) -> Result<(Vec<collagene_core::cloning::topo::Input>, Option<Molecule>), ErrorData> {
        if sources.len() > 32
            || request.recipe.materials.is_empty()
            || request.recipe.materials.len() > 16
            || request.steps.len() > 16
        {
            return Err(bad(
                "TOPO requires bounded explicit sources and current native recipe/materials/steps",
            ));
        }
        let mut roles = BTreeSet::new();
        for source in sources {
            if source.role.trim().is_empty()
                || source.role.len() > 256
                || source.role.chars().any(char::is_control)
                || !roles.insert(source.role.as_str())
            {
                return Err(bad(
                    "TOPO roles must be unique nonblank bounded identifiers",
                ));
            }
        }
        let store = self.session.store.lock().unwrap();
        let borrowed = sources
            .iter()
            .map(|s| store.document(&s.document))
            .collect::<Result<Vec<_>, _>>()?;
        let target = target.map(|t| store.document(t)).transpose()?;
        let mut aggregate = borrowed.clone();
        aggregate.extend(target);
        collagene_model::bounds::model_preflight(&aggregate).map_err(bad)?;
        let mut bases = 0usize;
        for molecule in &borrowed {
            if molecule.sequence.alphabet() != collagene_model::Alphabet::Dna {
                return Err(bad("TOPO requires explicit DNA sources"));
            }
            bases = bases
                .checked_add(molecule.sequence.len())
                .filter(|n| *n <= 100_000)
                .ok_or_else(|| bad("TOPO source aggregate exceeds100000 bases"))?;
        }
        if target.is_some_and(|m| {
            m.sequence.len() > 100_000 || m.sequence.alphabet() != collagene_model::Alphabet::Dna
        }) {
            return Err(bad("TOPO independent target requires bounded DNA"));
        }
        let inputs = sources
            .iter()
            .zip(borrowed)
            .map(|(s, m)| collagene_core::cloning::topo::Input {
                role: s.role.clone(),
                molecule: m.clone(),
            })
            .collect::<Vec<_>>();
        collagene_core::cloning::topo::source_preflight(&inputs, request).map_err(bad)?;
        Ok((inputs, target.cloned()))
    }
}

#[tool_router(router = assembly_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Freshly reconstruct every explicitly selected saved combinatorial row against its independent target and atomically release named native constructs, retained preparation/run events and optional oligo collections in a NEW immutable workspace. Per-row output options are explicit. A late failure discards the entire private prefix and publishes no workspace/artifact. Returns full historical release manifest plus native archive; imported reports never authorize release, and experimental success is not asserted."
    )]
    pub(crate) async fn release_combinatorial_outputs(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ReleaseQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_release(args).await
    }
    #[tool(
        description = "Export a saved native combinatorial release as a bounded ZIP of selected GenBank/FASTA, declared-primer CSV/TSV and HTML/PDF protocols. Optional construct IDs select a subset. Stable ordinal/hash entry paths never use names as paths. Protocol files require fresh independent-declaration and native-product checks; sequence/primer-only exchange remains clearly historical. All requested files or an explicit failure, with complete format losses and exact source/run/output events in manifest.json. No subprocess, arbitrary filesystem path or experimental-success claim."
    )]
    pub(crate) async fn export_combinatorial_outputs(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_archive(args).await
    }
    #[tool(
        description = "List retained native combinatorial release manifests with bounded paging. The records identify actual native output events and historical lineage; listing does not recompute science or restore current authority."
    )]
    pub(crate) async fn list_combinatorial_releases(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_release_list(args).await
    }
    #[tool(
        description = "Open one native combinatorial release manifest as historical JSON, including exact input/target/output events, per-row options and retained run identities. Leaves heads unchanged and never treats imported manifest status as fresh scientific authority."
    )]
    pub(crate) async fn open_combinatorial_release(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ReleaseRef>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_release_open(args).await
    }
    #[tool(
        description = "Page a declared alternative-source matrix in an imported native workspace, with exact retained source-event ownership and original interval checks. Uses stable bin/choice identities and bounded lazy enumeration (at most64 rows), not reaction-product enumeration. Does not test preparation, junctions, independent target agreement or experimental success."
    )]
    pub(crate) async fn page_combinatorial_design(
        &self,
        Parameters(args): Parameters<combinatorial_tools::MatrixQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_matrix(args).await
    }
    #[tool(
        description = "Resolve one explicitly selected construct's exact retained source events, original contributions, per-construct overrides, independent target declaration and typed native preparation recipe. Returns complete source-bound JSON evidence. This inspection does not reconstruct the causal assembly or establish verification; use fresh batch evaluation for that."
    )]
    pub(crate) async fn inspect_combinatorial_construct(
        &self,
        Parameters(args): Parameters<combinatorial_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_inspect(args).await
    }
    #[tool(
        description = "List immutable assembly-library edit recovery records in an imported native workspace, with bounded paging and stable historical hashes/allocation counters. Declared source pins are retained native event references, not complete host-draft attestation or scientific evidence. Does not apply recovery, rewind identities or change source heads."
    )]
    pub(crate) async fn list_assembly_library_history(
        &self,
        Parameters(args): Parameters<assembly_library_history_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.assembly_library_history_list(args).await
    }
    #[tool(
        description = "Open an exact immutable assembly-library recovery hash. Returns complete metadata JSON and the exact raw before-draft JSON in separate artifacts, preserving unsaved text without integer-array expansion. Historical recovery only: no apply/restore action, allocation rewind, complete host-draft validation or scientific authority is granted."
    )]
    pub(crate) async fn open_assembly_library_history(
        &self,
        Parameters(args): Parameters<assembly_library_history_tools::OpenQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.assembly_library_history_open(args).await
    }
    #[tool(
        description = "Discover the finite native Golden Gate reaction pool for one source-bound combinatorial construct, including whole sources, complete PCR amplicons, discarded fragments, carryover and explicit preparation substrates. Returns exact source/target event bindings and source_revision, candidate indices and both sequence keys; complete native report is retained in the JSON artifact. Missing selection is allowed and never defaults to candidate zero. Independent target correspondence is advisory, never a selection rule. Explicit selection followed by fresh batch evaluation and release is required; imported reports confer no authority."
    )]
    pub(crate) async fn inspect_combinatorial_pool(
        &self,
        Parameters(args): Parameters<combinatorial_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_pool_inspect(args).await
    }
    #[tool(
        description = "Inspect whole-parent multi-site BP/LR reference networks for one exact-source library construct. Resolve its independently declared target or freeze original contribution regions/junctions before reaction reconstruction. Returns a source_revision for explicit component selection, bounded candidate summaries and complete typed batch evidence. Missing selection never means candidate zero; matching the target never selects a pairing. Reference-only scope remains incomplete, with no physical recombination or experimental success claim."
    )]
    pub(crate) async fn inspect_combinatorial_multisite(
        &self,
        Parameters(args): Parameters<combinatorial_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_multisite_inspect(args).await
    }
    #[tool(
        description = "Freshly evaluate explicitly selected constructs with the existing deterministic preparation compilers/verifier and retain a named immutable batch record in a NEW native workspace. Stable design/construct/bin/choice identities and exact original events are preserved. Per-construct preparation/method/topology overrides are explicit; unresolved draft rows remain errors. Imported reports never restore fresh authority. Original source heads remain unchanged; native archive and complete historical JSON are returned atomically. This does not release bulk outputs or assert experimental success."
    )]
    pub(crate) async fn evaluate_combinatorial_batch(
        &self,
        Parameters(args): Parameters<combinatorial_tools::BatchQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_evaluate(args).await
    }
    #[tool(
        description = "List retained combinatorial batch records in an imported native workspace with bounded paging. Recorded statuses are historical evidence; listing never recomputes preparations or restores current verification authority."
    )]
    pub(crate) async fn list_combinatorial_batches(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_list(args).await
    }
    #[tool(
        description = "Open one immutable combinatorial batch record as complete historical JSON. Exact source/target bindings, original choices, per-construct settings and recorded outcomes remain inspectable without changing current heads. Opening does not rerun science or establish fresh authority."
    )]
    pub(crate) async fn open_combinatorial_batch(
        &self,
        Parameters(args): Parameters<combinatorial_tools::OpenQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_open(args).await
    }
    #[tool(
        description = "Freshly replay one explicitly selected row from a saved combinatorial batch using its exact original source events, independent target and preparation settings. Returns recorded and current row evidence separately. Leaves the workspace, batch record and source heads unchanged; current computational consistency does not establish experimental success or release bulk products."
    )]
    pub(crate) async fn replay_combinatorial_batch(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ReplayQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_replay(args).await
    }
    #[tool(
        description = "Freshly reconstruct TOPO strand preparation from explicit immutable source roles, source intervals or exact PCR/tailing, declared end chemistry/activation, pairing and requested transfers. No imported molecule metadata or foreign physical product grants activation/state. Returns bounded component/nick/unpaired/transfer summaries and a full <=32MiB source/request-bound JSON artifact, retaining every component. A declared graph computation is not kit efficacy, kinetics, manufacturing conservation or experimental success."
    )]
    pub(crate) async fn simulate_topo(
        &self,
        Parameters(args): Parameters<TopoSimulateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,_)=server.topo_sources(&args.sources,&native,None)?;
            let report=collagene_core::cloning::topo::simulate(&inputs,&native).map_err(bad)?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:None,report:&report,
                binding_scope:"Exact immutable MCP document roles/revisions; native event parents require the typed plan/run workflow."};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,collagene_core::cloning::topo::MAX_REPORT_BYTES)?,"application/json","topo-simulation.json")?;
            output(json!({"sources":args.sources,"engine_version":report.engine_version,"batch_revision":report.batch_revision,
                "source_identities":report.sources,"preparation_assumptions":report.preparation.assumptions,
                "unrepresented_reference_interval_count":report.preparation.not_represented.len(),"physical":topo_summary(&report.physical),
                "artifact_id":artifact,"authority":"Fresh computation under explicit preparation/activation declarations; no inferred kit state, winner or experimental validation."}))
        }).await
    }
    #[tool(
        description = "Inspect both declared TOPO vector/insert orientations from freshly reconstructed sources. Explicit prepared-strand indexes and directional flag; input transfers must be empty. Retains all physical components, failed transfers and alternative orientation graphs in full <=32MiB JSON. Returns exact orientation transfer requests and fresh source batch pins for explicit later projection, without choosing a winner or using target DNA; two orientations are not an exhaustive biochemical mixture search."
    )]
    pub(crate) async fn inspect_topo_orientations(
        &self,
        Parameters(args): Parameters<TopoOrientationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,_)=server.topo_sources(&args.sources,&native,None)?;
            let report=collagene_core::cloning::topo::inspect_orientations(&inputs,&native,args.vector,args.insert,args.directional).map_err(bad)?;
            let candidates=report.reference.candidates.iter().zip(&report.candidate_requests).enumerate().map(|(index,(physical,request))|{
                Ok::<_,ErrorData>(json!({"index":index,"batch_revision":collagene_core::cloning::topo::batch_revision(&inputs,request).map_err(bad)?,"steps":request.steps,"physical":topo_summary(physical)}))
            }).collect::<Result<Vec<_>,_>>()?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:None,report:&report,
                binding_scope:"Exact immutable MCP roles/revisions; complete candidate requests are retained, never target-selected."};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,collagene_core::cloning::topo::MAX_REPORT_BYTES)?,"application/json","topo-orientations.json")?;
            output(json!({"sources":args.sources,"engine_version":report.engine_version,"batch_revision":report.batch_revision,
                "orientation_search_complete":report.reference.orientation_search_complete,"biochemical_mixture_exhaustive":report.reference.biochemical_mixture_exhaustive,
                "candidate_count":candidates.len(),"candidates":candidates,"artifact_id":artifact,
                "authority":"All declared orientation alternatives retained; no preferred winner, target-derived selection or exhaustive mixture claim."}))
        }).await
    }
    #[tool(
        description = "Freshly reconstruct one explicitly selected TOPO component and source-call seed pinned to the source/recipe batch. Projection is withheld unless a supported paired contour can be proven; retains full physical graph, nicks, unpaired DNA, residual components and unapplied transfers. A projected immutable sequence document carries conservative source-coordinate annotations and remains a chemistry-unmapped view, atomically installed with <=32MiB full graph/annotation evidence; it cannot substitute for graph replay or independent target comparison. Imported reports/products are never accepted as physical state."
    )]
    pub(crate) async fn project_topo_product(
        &self,
        Parameters(args): Parameters<TopoProjectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,_)=server.topo_sources(&args.sources,&native,None)?;
            let product=collagene_core::cloning::topo::product(&inputs,&native,&args.selection).map_err(bad)?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:None,report:&product,
                binding_scope:"Source-backed physical reconstruction; sequence view is not standalone chemistry authority."};
            let (artifact,document)=server.session.candidate_document_artifact_bounded(&evidence,product.reference.clone(),"topo-projection.json",collagene_core::cloning::topo::MAX_REPORT_BYTES)?;
            output(json!({"sources":args.sources,"selection":args.selection,"engine_version":product.engine_version,
                "batch_revision":product.batch_revision,"projection":topo_projection_summary(&product),
                "physical":topo_summary(&product.projection.physical),"annotations":topo_annotation_summary(product.annotations.as_ref()),"document":document,"artifact_id":artifact,
                "authority":"Explicit component/seed only; no target-guided path selection, nick repair or standalone physical authority from the projected reference."}))
        }).await
    }
    #[tool(
        description = "Freshly reconstruct the explicitly selected TOPO component/seed and compare an independently supplied immutable target under a fully explicit AcceptancePolicy. Target DNA is a noncausal requirement and never chooses components, seeds, orientation or product sequence. Exact equality cannot silently seal nicks, specify unknown terminal chemistry, discard other components or apply failed transfers. Returns conditional accepted/rejected/withheld outcome and complete source/request/target-bound <=32MiB JSON; no experimental success guarantee."
    )]
    pub(crate) async fn check_topo_target(
        &self,
        Parameters(args): Parameters<TopoTargetQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,target)=server.topo_sources(&args.sources,&native,Some(&args.target))?;
            let target=target.ok_or_else(||bad("explicit independent TOPO target required"))?;
            let report=collagene_core::cloning::topo::check_target(&inputs,&native,&args.selection,&target,&args.policy).map_err(bad)?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:Some(&args.target),report:&report,
                binding_scope:"Independent target is a noncausal requirement, not an ancestry parent or source of generated DNA."};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,collagene_core::cloning::topo::MAX_REPORT_BYTES)?,"application/json","topo-target-check.json")?;
            output(json!({"sources":args.sources,"target":args.target,"selection":args.selection,"policy":args.policy,
                "engine_version":report.engine_version,"batch_revision":report.batch_revision,"target_revision":report.target_revision,
                "outcome":report.outcome,"projection":topo_projection_summary(&report.product),
                "physical":topo_summary(&report.product.projection.physical),"artifact_id":artifact,"annotations":topo_annotation_summary(report.product.annotations.as_ref()),
                "authority":"Conditional comparison under declared physical graph and explicit allowances, not a repaired molecule or experimental verification."}))
        }).await
    }
    #[tool(
        description = "Inspect explicitly profiled BP/LR reciprocal reference-sequence alternatives from exactly two ordered immutable source roles. Mandatory full profile snapshot and resource limits; all reciprocal components and uncertainty retained in immutable JSON. Bounded candidate page only, no automatic winner or product installation, no chemical activity claims."
    )]
    pub(crate) async fn inspect_recombination(
        &self,
        Parameters(args): Parameters<recombination_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.recombination_inspect(args).await
    }
    #[tool(
        description = "Freshly compare independently supplied target DNA against all explicitly profiled BP/LR reference components. Preserves possible sites, unsupported geometry, ambiguity and incomplete-search flags; no-match claims are scoped to the declared full-site profile. The target is noncausal and never supplies product DNA."
    )]
    pub(crate) async fn check_recombination_target(
        &self,
        Parameters(args): Parameters<recombination_tools::CheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.recombination_check(args).await
    }
    #[tool(
        description = "Freshly reconstruct one exact candidate/component pinned to its source/profile batch and require a fully specified independent target match. Atomically installs its immutable reference-only document and full checked JSON; all reciprocal component identities retained. No physical reaction, strand state, efficacy or experimental winner is inferred. Native history uses the typed plan/run workflow."
    )]
    pub(crate) async fn select_recombination_component(
        &self,
        Parameters(args): Parameters<recombination_tools::SelectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.recombination_select(args).await
    }
    #[tool(
        description = "Inspect simultaneous multi-parent reciprocal BP/LR reference alternatives from 2..5 exact immutable source roles using a complete caller-supplied finite profile and explicit specificity set. Optional native source pins are attested by owner/event/revision and literal molecule equality. Returns bounded candidate/component summaries and full typed batch/annotation JSON; no curated biological profile, automatic winner or physical activity claim."
    )]
    pub(crate) async fn inspect_multisite_recombination(
        &self,
        Parameters(args): Parameters<multisite_recombination_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.multisite_inspect(args).await
    }
    #[tool(
        description = "Freshly reconstruct one explicitly selected multisite reference component pinned to its exact source/profile/options batch. Atomically retains the projected reference document and full JSON with exact source bindings, reciprocal identities and conservative annotations. No independent target match, repaired state or physical reaction is inferred; this is reference-only material."
    )]
    pub(crate) async fn select_multisite_recombination_component(
        &self,
        Parameters(args): Parameters<multisite_recombination_tools::ProductQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.multisite_select(args).await
    }
    #[tool(
        description = "Freshly compare one explicitly selected multisite reference component against an independent immutable target. Target never determines site matching, component selection or product DNA. Retains full typed evidence, mismatch and incomplete/unsupported search diagnostics; does not install documents or establish chemical success."
    )]
    pub(crate) async fn check_multisite_recombination_target(
        &self,
        Parameters(args): Parameters<multisite_recombination_tools::CheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.multisite_check(args).await
    }
    #[tool(
        description = "Search bounded signed Golden Gate fusion-placement windows while retaining the original independently declared intended DNA and fixed manual preparations. Mathematical palindrome/end-distance constraints are not ligation fidelity or PCR feasibility. Returns complete candidate evidence as an artifact; explicitly select and materialize a candidate before fresh recipe verification."
    )]
    pub(crate) async fn search_golden_gate_fusions(
        &self,
        Parameters(args): Parameters<fusion_design_tools::SearchQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.fusion_search(args).await
    }
    #[tool(
        description = "Freshly materialize an explicitly selected source-bound Golden Gate fusion placement using original-template PCR and pinned Type IIS digests. Retains original intended DNA, records primer-added prefixes/suffixes, and preserves fixed preparations. Advisory only: verify the full returned ligation recipe against the original independent target before publication."
    )]
    pub(crate) async fn materialize_golden_gate_fusions(
        &self,
        Parameters(args): Parameters<fusion_design_tools::MaterializeQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.fusion_materialize(args).await
    }
    #[tool(
        description = "Estimate conditional ligation fidelity from an explicitly supplied, provenance-bound directed count matrix and selected fusion ends. Reports conditions and missing-data uncertainty; does not infer empirical data from exact complementarity or claim experimental yield. No external optimizer or runtime software is called."
    )]
    pub(crate) async fn score_ligation_fidelity(
        &self,
        Parameters(args): Parameters<fusion_design_tools::FidelityQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.fidelity_score(args).await
    }
    #[tool(
        description = "List exact versioned published ligation-frequency profiles, including formulation, conditions, source/conversion hashes, citation and separate CC-BY attribution. These data do not establish the actual preparation conditions or reaction efficiency."
    )]
    pub(crate) async fn list_ligation_profiles(
        &self,
        Parameters(args): Parameters<published_ligation_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.published_profiles_list(args).await
    }
    #[tool(
        description = "Score an explicit set of fusion ends using one exact published symmetric pair-frequency profile. Preserves literal published weights, zeros, palindromic/unknown cases and full conditions. No raw independent directional reads, measured yield or chemical verification is inferred."
    )]
    pub(crate) async fn score_published_ligation_profile(
        &self,
        Parameters(args): Parameters<published_ligation_tools::ScoreQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.published_profile_score(args).await
    }
    #[tool(
        description = "Freshly search source-bound Golden Gate fusion placements once and rank retained sets using one exact published profile prepared once. Requires an explicit conditional formulation assumption; never substitutes formulation from a generic enzyme name, selects a reaction winner, changes intended DNA or applies primers. Complete source/profile-bound evidence is returned as a bounded artifact."
    )]
    pub(crate) async fn rank_golden_gate_fusions(
        &self,
        Parameters(args): Parameters<published_ligation_tools::RankQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.published_fusion_rank(args).await
    }
    #[tool(
        description = "Freshly inspect TOPO insert PCR from one deliberately selected saved native material origin. Requires immutable source references, an explicit template role, a small explicit Run/Derived material selection and its exact retained native workspace. Replays the origin before actual new PCR; graph metadata or matching DNA cannot substitute for material authority. Select an explicit orientations, projection, independent-target Check, actions inspection or processed-actions task. Check resolves an exact immutable session target reference, never replaces it with a projected product. Returns bounded full typed evidence without saving any candidate or document; kit activation, transfers, processing and acceptance remain explicit assumptions."
    )]
    pub(crate) async fn inspect_saved_topo(
        &self,
        Parameters(args): Parameters<SavedTopoInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.len() > 32 {
            return Err(bad("bounded explicit TOPO input roles required"));
        }
        bounded_json_limit(&args.task, 2 * 1024 * 1024)?;
        if args
            .template
            .as_ref()
            .is_some_and(|t| t.role.trim().is_empty() || t.role.len() > 256)
        {
            return Err(bad("bounded saved PCR-template role required"));
        }
        if args.template.is_none() && args.task.preparation().selections()?.is_empty() {
            return Err(bad("choose an explicit saved graph or PCR-template origin"));
        }
        let server = self.clone();
        self.compute(move||{
            let selected=args.template.iter().map(|t|t.selection.clone()).collect::<Vec<_>>();
            let (preparation,mut materials)=server.hydrate_topo(args.task.preparation(),Some(&args.origin_workspace),&selected)?;
            let template=if let Some(t)=&args.template {
                if t.role.trim().is_empty()||t.role.len()>256||materials.len()!=1{return Err(bad("one exact saved PCR-template role required"));}
                Some(collagene_core::cloning::topo::run_preparation::BoundTemplate{role:t.role.clone(),origin:materials.remove(0)})
            }else{None};
            let (inputs,target)={let store=server.session.store.lock().unwrap();let mut ids=BTreeSet::new();let mut refs=vec![];let mut sources=vec![];
                for b in &args.inputs {if b.id.trim().is_empty()||b.id.len()>128||!ids.insert(b.id.as_str()){return Err(bad("invalid or duplicate saved TOPO role"));}refs.push(store.document(&b.document)?);sources.push(collagene_core::planning::plan::TopoSource{role:b.id.clone(),source:collagene_core::planning::plan::Source::Input(b.id.clone())});}
                let target=if let SavedTopoTask::Check{target,..}=&args.task {let m=store.document(target)?;collagene_core::cloning::topo::worker::target_preflight(m).map_err(bad)?;Some(m.clone())}else{None};
                let inputs=collagene_core::cloning::topo::run_preparation::source_inputs(&sources,&refs,&preparation,template.as_ref()).map_err(bad)?;
                (inputs,target)
            };
            let task=args.task.native(target,preparation)?;
            let context=collagene_core::cloning::topo::run_preparation::revision(&inputs,template.as_ref(),&task).map_err(bad)?;
            let report=collagene_core::cloning::topo::run_preparation::inspect(&inputs,template.as_ref(),&task).map_err(bad)?;
            report.validate_for(&inputs,template.as_ref(),&task,&context).map_err(bad)?;
            #[derive(Serialize)]struct Evidence<'a>{inputs:&'a[BoundInput],template:&'a Option<topo_graph_tools::BoundTemplate>,origin_workspace:&'a WorkspaceRef,task:&'a SavedTopoTask,inspection:&'a collagene_core::cloning::topo::run_preparation::Inspection}
            let evidence=Evidence{inputs:&args.inputs,template:&args.template,origin_workspace:&args.origin_workspace,task:&args.task,inspection:&report};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","saved-topo-inspection.json")?;
            output(json!({"context_revision":context,"engine_version":collagene_core::cloning::topo::run_preparation::ENGINE_VERSION,"template":args.template,"origin_workspace":args.origin_workspace,"origin_incomplete":report.origin_incomplete,"origin_findings":report.origin_findings,"origin_findings_total":report.origin_findings_total,"origin_findings_truncated":report.origin_findings_truncated,"artifact_id":artifact,"verified":false,"authority":"Exact retained native graph and optional PCR-template origins freshly reconstructed; whole graph/prefix/byproducts retained. No duplicate insert synthesis, imported report authority, isolation or experimental success."}))
        }).await
    }
    #[tool(
        description = "List, freshly inspect, or explicitly capture a whole native TOPO graph in a pinned workspace. Source/event/record selections are small and session-attested; no supplied graph report or metadata becomes material. Capture returns a NEW immutable archive/receipt without changing heads or live application state. Retains partial-transfer history, activations, nicks, byproducts and inherited limitations. Retains the exact ordered transfer, flap-cleavage and nick-sealing transcript; no inventory, activation manufacture, isolation or experiment claim."
    )]
    pub(crate) async fn topo_graph(
        &self,
        Parameters(args): Parameters<topo_graph_tools::Query>,
    ) -> Result<CallToolResult, ErrorData> {
        self.topo_graph_dispatch(args).await
    }
    #[tool(
        description = "Compile explicit source-backed TA or directional CACC TOPO preparation declarations into the shared native recipe. Preserves actual ordered primers, vector strand geometry, terminal states and activation assumptions. Does not simulate activation or repair, add a tail silently, select transfers or invent a closed clone. Reuse inspect_topo_orientations and existing projection/target checks with the returned fields."
    )]
    pub(crate) async fn compile_topo_preparation(
        &self,
        Parameters(args): Parameters<published_ligation_tools::TopoQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.topo_compile_preparation(args).await
    }

    #[tool(
        description = "Freshly reconstruct source-bound chemical actions and prove a selected paired contour. Uses exact processing context/component/source seed and preserves the excised fragment, all parent/child graphs, created end chemistry and Unknown modifications. An annotated sequence view cannot restore private material; use the TopoActions plan operation with an independent target for full verification."
    )]
    pub(crate) async fn project_topo_actions(
        &self,
        Parameters(args): Parameters<topo_actions_tools::ProductQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.project_actions(args).await
    }
    #[tool(
        description = "Inspect source-bound restriction sites, explicit retained/protected ranges and coding-feature constraints for internal-site domestication. Uses exact catalog/feature/source pins and reports unsupported geometry or translation semantics rather than assuming a silent change. Reference inspection only; no physical stock is made."
    )]
    pub(crate) async fn inspect_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::DesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_inspect(args).await
    }
    #[tool(
        description = "Propose bounded synonymous source-coordinate substitutions that remove scoped unwanted restriction sites while preserving declared sites, protected bases and affected supported coding features. Positional minimum-edit constraints are inside the shared native search. Pruning remains explicit; a candidate is reference intent, not prepared DNA."
    )]
    pub(crate) async fn propose_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::DesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_propose(args).await
    }
    #[tool(
        description = "Independently check explicit source-bound domestication substitutions, preserving original annotation coordinates and qualifiers and rescanning affected proteins and the entire DNA for scoped unwanted/new sites. Retains the independently declared target and evidence atomically. This check confers no physical-preparation or experimental authority."
    )]
    pub(crate) async fn check_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::CheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_check(args).await
    }
    #[tool(
        description = "Materialize a literal checked domestication target using an explicit exact-source route: inverse PCR on a circle, sequential inverse-PCR rounds, or whole-source overlap-extension PCR on linear/circular DNA. Request origin is required-nullable; saved Run/Derived sources require the exact retained origin_workspace before planning. OE performs two actual source-template PCR branches, physical extension and outer PCR; circles then require explicit phosphorylation and ligation. All substitutions, coding/protected constraints and scoped sites are checked before planning. Atomically retains the independently frozen target, one canonical native plan and full route/oligo/search evidence for fresh verify_plan and HTML/PDF/order export. Unresolved PCR geometry retains a nonready diagnostic candidate; unsupported sources, settings and failed coverage are refused. No invented recoded stock or imported report authority."
    )]
    pub(crate) async fn prepare_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::PrepareQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_prepare(args).await
    }
    #[tool(
        description = "Preview an explicit physical homology-assembly recipe from original immutable sources and a separately bound intended target. Reconstructs actual preparations, oriented source strands, declared5-prime resection, complementary fill and sealing under an explicit finite policy; a vendor Gibson/HiFi/InFusion label supplies advisory sequence guidance only. Recipe overlaps/topology must exactly match the physical policy, never be silently overwritten. Retains products, rejected/ambiguous chemistry and bounded-search findings. Every FromRun role requires its exact native workspace and fresh closed-origin replay. Returns complete bounded plan/physical/junction evidence and concise summaries; it saves no candidate, source, target or product document. To keep a design, propose the retained typed operations against these same inputs/target and use normal verification/review."
    )]
    pub(crate) async fn preview_homology_recipe(
        &self,
        Parameters(args): Parameters<HomologyRecipePreviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty()
            || args.inputs.len() > collagene_core::cloning::assembly::recipe::MAX_ROWS
            || args.origin_workspaces.len() > collagene_core::cloning::assembly::recipe::MAX_ROWS
        {
            return Err(bad(
                "homology recipe sources/workspace bindings exceed their bounded row scope",
            ));
        }
        let origins = args
            .recipe
            .fragments
            .iter()
            .filter(|fragment| {
                matches!(
                    fragment.preparation,
                    collagene_core::cloning::assembly::recipe::Preparation::FromRun { .. }
                )
            })
            .count();
        if origins != args.origin_workspaces.len() {
            return Err(bad(
                "supply exactly one native workspace for each saved-source role",
            ));
        }
        let mut seen = BTreeSet::new();
        for binding in &args.origin_workspaces {
            if binding.role.is_empty()
                || binding.role.len() > 128
                || !seen.insert(binding.role.as_str())
                || !args.recipe.fragments.iter().any(|fragment| {
                    fragment.role == binding.role
                        && matches!(
                            fragment.preparation,
                            collagene_core::cloning::assembly::recipe::Preparation::FromRun { .. }
                        )
                })
            {
                return Err(bad(
                    "homology origin workspace role is invalid, duplicated or unused",
                ));
            }
        }
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)={let store=server.session.store.lock().unwrap();let target=store.document(&args.target)?;let mut ids=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {if input.id.trim().is_empty()||input.id.len()>128||!ids.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::cloning::assembly::recipe::homology_source_preflight(&sources,target,&args.recipe,&args.policy).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())};
            for fragment in &args.recipe.fragments {
                if let collagene_core::cloning::assembly::recipe::Preparation::FromRun{origin,..}=&fragment.preparation {
                    let workspace=args.origin_workspaces.iter().find(|binding|binding.role==fragment.role).ok_or_else(||bad("saved homology source has no exact workspace binding"))?;
                    let source=inputs.iter().find(|input|input.id==fragment.input_id).ok_or_else(||bad("saved homology source input unavailable"))?;
                    server.session.attest_run_source(&workspace.workspace,origin,&source.molecule)?;
                }
            }
            let preview=collagene_core::cloning::assembly::recipe::preview_homology(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            let context=collagene_core::cloning::assembly::recipe::homology_preview_revision(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            preview.validate_binding(&context).map_err(bad)?;
            let diagnostics=collagene_core::cloning::assembly::diagnostics::overlap(&preview).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{inputs:&'a[BoundInput],target:&'a DocumentRef,recipe:&'a collagene_core::cloning::assembly::recipe::Recipe,policy:&'a collagene_core::cloning::overlap::assembly::Options,origin_workspaces:&'a[MaterialOriginWorkspace],preview:&'a collagene_core::cloning::assembly::recipe::Preview,diagnostics:&'a collagene_core::cloning::assembly::diagnostics::Report}
            let evidence=Evidence{inputs:&args.inputs,target:&args.target,recipe:&args.recipe,policy:&args.policy,origin_workspaces:&args.origin_workspaces,preview:&preview,diagnostics:&diagnostics};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","physical-homology-recipe.json")?;
            let product=preview.verification.final_product();
            let physical=preview.verification.physical_processing_reports.iter().find_map(|step|match &step.evidence{verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(report) if step.step_index==preview.compiled.assembly_step=>Some(report),_=>None});
            output(json!({"state":"preview","verified":false,"candidate_saved":false,"inputs":args.inputs,"target":args.target,
                "context_revision":context,"plan_revision":preview.verification.plan_revision,"target_revision":preview.verification.target_revision,
                "engine_version":preview.verification.engine_version,"status":preview.verification.status,"preparation_complete":preview.preparation_complete,
                "product":product.map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                "physical":physical.map(|report|json!({"engine_version":report.engine_version,"paired_product_available":report.product.is_some(),"reusable_material_available":report.product.as_ref().is_some_and(|product|product.nicks.is_empty()),"terminal_search_complete":report.terminal_search_complete,"incomplete":report.incomplete,"comparisons":report.comparisons,"source_geometry_revisions":report.source_geometry_revisions})),
                "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                "findings_truncated":preview.verification.findings.len()>100,
                "policy":args.policy,"diagnostics":diagnostics,"artifact_id":artifact,"plan_artifact_path":"preview.compiled.plan",
                "authority":"Fresh explicit physical-model preview and independent-target comparison only. Paired reference DNA can remain diagnostic when joining termini cannot seal; it does not grant reusable material. Eligible material requires fresh native reconstruction under the same policy and explicit saved-origin replay for reuse. This response saves no candidate or product and makes no experimental-success, isolation or mixture-uniqueness claim."
            }))
        }).await
    }
    #[tool(
        description = "Preview explicit 3-prime resection and homology annealing from original immutable sources and an independent target. Requires actual preparations, per-source left/right strand cuts, explicit junction lengths and bounded pairing search. Preserves remaining strands, source coordinates, unfilled gaps and unsealed nicks. The joined reference always remains incomplete and cannot confer reusable assembled material; no polymerase, ligase or host processing is inferred. FromRun roles require exact native workspace attestations and fresh closed-origin replay. Returns complete bounded typed plan, physical and junction evidence without saving any candidate or document."
    )]
    pub(crate) async fn preview_annealing_recipe(
        &self,
        Parameters(args): Parameters<AnnealingRecipePreviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty()
            || args.inputs.len() > collagene_core::cloning::assembly::recipe::MAX_ROWS
            || args.origin_workspaces.len() > collagene_core::cloning::assembly::recipe::MAX_ROWS
        {
            return Err(bad(
                "homology recipe sources/workspace bindings exceed their bounded row scope",
            ));
        }
        let origins = args
            .recipe
            .fragments
            .iter()
            .filter(|fragment| {
                matches!(
                    fragment.preparation,
                    collagene_core::cloning::assembly::recipe::Preparation::FromRun { .. }
                )
            })
            .count();
        if origins != args.origin_workspaces.len() {
            return Err(bad(
                "supply exactly one native workspace for each saved-source role",
            ));
        }
        let mut seen = BTreeSet::new();
        for binding in &args.origin_workspaces {
            if binding.role.is_empty()
                || binding.role.len() > 128
                || !seen.insert(binding.role.as_str())
                || !args.recipe.fragments.iter().any(|fragment| {
                    fragment.role == binding.role
                        && matches!(
                            fragment.preparation,
                            collagene_core::cloning::assembly::recipe::Preparation::FromRun { .. }
                        )
                })
            {
                return Err(bad(
                    "homology origin workspace role is invalid, duplicated or unused",
                ));
            }
        }
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)={let store=server.session.store.lock().unwrap();let target=store.document(&args.target)?;let mut ids=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {if input.id.trim().is_empty()||input.id.len()>128||!ids.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::cloning::assembly::recipe::annealing_source_preflight(&sources,target,&args.recipe,&args.policy).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())};
            for fragment in &args.recipe.fragments {
                if let collagene_core::cloning::assembly::recipe::Preparation::FromRun{origin,..}=&fragment.preparation {
                    let workspace=args.origin_workspaces.iter().find(|binding|binding.role==fragment.role).ok_or_else(||bad("saved homology source has no exact workspace binding"))?;
                    let source=inputs.iter().find(|input|input.id==fragment.input_id).ok_or_else(||bad("saved homology source input unavailable"))?;
                    server.session.attest_run_source(&workspace.workspace,origin,&source.molecule)?;
                }
            }
            let preview=collagene_core::cloning::assembly::recipe::preview_annealing(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            let context=collagene_core::cloning::assembly::recipe::annealing_preview_revision(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            preview.validate_binding(&context).map_err(bad)?;
            let diagnostics=collagene_core::cloning::assembly::diagnostics::overlap(&preview).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{inputs:&'a[BoundInput],target:&'a DocumentRef,recipe:&'a collagene_core::cloning::assembly::recipe::Recipe,policy:&'a collagene_core::cloning::overlap::annealing::Options,origin_workspaces:&'a[MaterialOriginWorkspace],preview:&'a collagene_core::cloning::assembly::recipe::Preview,diagnostics:&'a collagene_core::cloning::assembly::diagnostics::Report}
            let evidence=Evidence{inputs:&args.inputs,target:&args.target,recipe:&args.recipe,policy:&args.policy,origin_workspaces:&args.origin_workspaces,preview:&preview,diagnostics:&diagnostics};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","three-prime-annealing-recipe.json")?;
            let product=preview.verification.final_product();
            let physical=preview.verification.physical_processing_reports.iter().find_map(|step|match &step.evidence{verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(report) if step.step_index==preview.compiled.assembly_step=>Some(report),_=>None});
            output(json!({"state":"preview","verified":false,"candidate_saved":false,"inputs":args.inputs,"target":args.target,
                "context_revision":context,"plan_revision":preview.verification.plan_revision,"target_revision":preview.verification.target_revision,
                "engine_version":preview.verification.engine_version,"status":preview.verification.status,"preparation_complete":preview.preparation_complete,
                "product":product.map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                "physical":physical.map(|report|json!({"engine_version":report.engine_version,"paired_diagnostic_available":report.paired_diagnostic.is_some(),"reusable_material_available":false,"terminal_search_complete":report.terminal_search_complete,"incomplete":report.incomplete,"comparisons":report.comparisons,"source_geometry_revisions":report.source_geometry_revisions})),
                "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                "findings_truncated":preview.verification.findings.len()>100,
                "policy":args.policy,"diagnostics":diagnostics,"artifact_id":artifact,"plan_artifact_path":"preview.compiled.plan",
                "authority":"Fresh 3-prime resection and homology-annealing preview against an independent target. The joined reference remains incomplete and documentary; all gaps and nicks remain, no fill, sealing or host activity is inferred, and no reusable assembled physical material is created. This response saves no candidate or product and makes no experimental-success, isolation or mixture-uniqueness claim."
            }))
        }).await
    }
    #[tool(
        description = "Compile an explicit ordered overlap cloning recipe from immutable source document revisions and an independent target. Each fragment declares WholeLinear, an explicit DeclaredDuplex original-source terminal assumption, complete Digest or exact CatalogDigest selection with chosen fragment/optional repair, or ExactPcr with actual ordered oligos, or PcrCatalogDigest chaining actual original-source PCR into a pinned catalog digest; orientations, homology lengths, topology and published-method advisory profile are explicit. Freshly reconstructs preparations and verifies the resulting typed plan. Returns additive typed diagnostics locating failed/blocked preparation rows, junction evidence and final target mismatch, without parsing findings. Atomically saves an immutable candidate and full preparation/junction/verification JSON; no guessed tails, automatic winner, source mutation or laboratory guarantee. Catalog-digest preparations require an explicit nullable conditions field: null leaves setup unspecified; selected formulations and assay choices are retained and freshly reassessed without changing cleavage reconstruction. Unknown or below-threshold documentary constraints remain incomplete, not evidence of experimental failure or product absence. Reuse verify_plan and HTML/PDF exports with its plan_revision. Explicit MismatchPcr retains both ordered primers and its declared mismatch/3′-exact/search policy. OverlapExtensionPcr retains two actual PCR branches from the same source, selected physical overlap extension and separate outer-primer PCR; every ordered oligo and terminal-pair search policy is explicit. OverlapExtensionPcrCatalogDigest appends exact catalog-pinned cleavage of that actual outer-PCR product, explicit fragment index and optional repair; inspect_overlap_extension_digest exposes the same real product without requiring a target. FromRun preserves a deliberately selected saved origin before downstream preparation. Mutation intent and the target remain independently declared; no mutation is selected to fit a simulated product."
    )]
    pub(crate) async fn propose_overlap_assembly(
        &self,
        Parameters(args): Parameters<OverlapAssemblyProposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad("overlap recipe requires 1..128 explicit source inputs"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.session.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len()>128 || !seen.insert(bound.id.as_str()) {
                        return Err(bad("empty, oversized or duplicate recipe input ID"));
                    }
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::cloning::assembly::recipe::source_preflight(&sources,target,&args.recipe).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,molecule)|PlanInput {
                    id:id.into(),revision:molecule.revision(),molecule:molecule.clone()
                }).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let preview = collagene_core::cloning::assembly::recipe::preview(&inputs,&target,&args.recipe).map_err(bad)?;
            let expected = collagene_core::cloning::assembly::recipe::preview_revision(&inputs,&target,&args.recipe).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            // Typed locations come from the checked execution prefix, never finding strings.
            let diagnostics = collagene_core::cloning::assembly::diagnostics::overlap(&preview).map_err(bad)?;
            if diagnostics.context_revision != expected { return Err(bad("assembly diagnostics belong to another preview context")); }
            // Preserve all existing top-level preview fields, adding diagnostics only.
            #[derive(Serialize)]
            struct Evidence<'a> {
                #[serde(flatten)]
                preview: &'a collagene_core::cloning::assembly::recipe::Preview,
                diagnostics: &'a collagene_core::cloning::assembly::diagnostics::Report,
            }
            // Serialize all outputs and capacity-check before either is published.
            let text = bounded_json_limit(&Evidence { preview: &preview, diagnostics: &diagnostics },32*1024*1024)?;
            let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
            let revision = collagene_core::planning::plan::DesignPlan::revision(&preview.compiled.plan);
            if revision.is_empty() {return Err(bad("compiled plan has no bounded canonical revision"));}
            let plan_bytes = bounded_json(&preview.compiled.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("candidate storage accounting overflow"))?;
            let result = json!({"plan_revision":revision,"state":"candidate","verified":false,
                "target":args.target,"inputs":args.inputs,"context_revision":preview.context_revision,
                "preparation_complete":preview.preparation_complete,"diagnostics":diagnostics,"assembly_step":preview.compiled.assembly_step,
                "verification":{"status":preview.verification.status,"engine_version":preview.verification.engine_version,
                    "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                    "findings_truncated":preview.verification.findings.len()>100,
                    "assumptions":preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),
                    "assumptions_truncated":preview.verification.assumptions.len()>100,
                    "exhaustive_alternative_search":preview.verification.exhaustive_alternative_search},
                "artifact_id":artifact_id,
                "authority":"Fresh deterministic preparation and declared ordered-plan evidence; saving remains a candidate, method guidance is advisory, and sequence agreement is not laboratory success or exhaustive assembly uniqueness. Full rows, prepared products and junction evidence are in the artifact."});
            let response = output(result)?;
            let mut store=server.session.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id) {
                if !matches!(existing.data,ArtifactData::Text(_)) || existing.mime_type!="application/json" || existing.filename!="overlap-assembly-recipe.json" {
                    return Err(bad("assembly evidence already has another artifact representation"));
                }
            }
            if (new_plan && store.plans.len()>=MAX_ITEMS) || (new_artifact && store.artifacts.len()>=MAX_ITEMS) {
                return Err(bad("candidate/artifact item limit; no assembly outputs created"));
            }
            let bytes=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("session storage accounting overflow"))?;
            store.reserve(bytes)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:preview.compiled.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"overlap-assembly-recipe.json"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Compile an explicit ordered ligation cloning recipe from immutable DNA source revisions and an independent target. Each fragment declares WholeLinear, explicit DeclaredDuplex terminal assumptions on an original intact linear DNA input, complete Digest or catalog-pinned CatalogDigest with optional end repair, or actual ExactPcr preparation, explicit terminal treatment and orientation; circular and allow_nicked policies are explicit. Freshly reconstructs physical strands, oriented ends, junctions and retained nicks through the deterministic verifier. Returns additive typed row/step/junction diagnostics without treating available endpoints as proof of compatibility. Atomically stores a candidate and full preparation/physical/verification JSON, including failed drafts. No metadata-derived chemistry authority, automatic treatment/repair, experimental guarantee or source mutation. Catalog-digest preparations require an explicit nullable conditions field: null leaves setup unspecified; selected formulations and assay choices are retained and freshly reassessed without changing cleavage reconstruction. Unknown or below-threshold documentary constraints remain incomplete, not evidence of experimental failure or product absence. Reuse verify_plan and HTML/PDF exports with plan_revision. Explicit MismatchPcr retains both ordered primers and its declared mismatch/3′-exact/search policy. OverlapExtensionPcr retains two actual PCR branches from the same source, selected physical overlap extension and separate outer-primer PCR; every ordered oligo and terminal-pair search policy is explicit. OverlapExtensionPcrCatalogDigest appends exact catalog-pinned cleavage of that actual outer-PCR product, explicit fragment index and optional repair; inspect_overlap_extension_digest exposes the same real product without requiring a target. FromRun preserves a deliberately selected saved origin before downstream preparation. Mutation intent and the target remain independently declared; no mutation is selected to fit a simulated product."
    )]
    pub(crate) async fn propose_ligation_assembly(
        &self,
        Parameters(args): Parameters<LigationAssemblyProposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad(
                "ligation recipe requires 1..128 explicit source inputs",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.session.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len()>128 || !seen.insert(bound.id.as_str()) {
                        return Err(bad("empty, oversized or duplicate recipe input ID"));
                    }
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::cloning::ligation::recipe::source_preflight(&sources,target,&args.recipe).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,molecule)|PlanInput {
                    id:id.into(),revision:molecule.revision(),molecule:molecule.clone()
                }).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let preview = collagene_core::cloning::ligation::recipe::preview(&inputs,&target,&args.recipe, collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete).map_err(bad)?;
            let expected = collagene_core::cloning::ligation::recipe::preview_revision(&inputs,&target,&args.recipe, collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            // Typed locations come from the checked execution prefix, never finding strings.
            let diagnostics = collagene_core::cloning::assembly::diagnostics::ligation(&preview).map_err(bad)?;
            if diagnostics.context_revision != expected { return Err(bad("assembly diagnostics belong to another preview context")); }
            // Preserve all existing top-level preview fields, adding diagnostics only.
            #[derive(Serialize)]
            struct Evidence<'a> {
                #[serde(flatten)]
                preview: &'a collagene_core::cloning::ligation::recipe::Preview,
                diagnostics: &'a collagene_core::cloning::assembly::diagnostics::Report,
            }
            // Serialize all outputs and capacity-check before either is published.
            let text = bounded_json_limit(&Evidence { preview: &preview, diagnostics: &diagnostics },32*1024*1024)?;
            let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
            let revision = collagene_core::planning::plan::DesignPlan::revision(&preview.compiled.plan);
            if revision.is_empty() {return Err(bad("compiled plan has no bounded canonical revision"));}
            let plan_bytes = bounded_json(&preview.compiled.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("candidate storage accounting overflow"))?;
            let result = json!({"plan_revision":revision,"state":"candidate","verified":false,
                "target":args.target,"inputs":args.inputs,"context_revision":preview.context_revision,
                "preparation_complete":preview.preparation_complete,"diagnostics":diagnostics,"ligation_step":preview.compiled.ligation_step,
                "fragment_count":preview.compiled.rows.len(),"junction_count":preview.junctions.len(),
                "product_physical":preview.product_physical.as_ref().map(|physical|json!({"top_length":physical.top_length,"bottom_length":physical.bottom_length,"topology":physical.topology,"nick_count":physical.nicks.len(),"covalently_closed_in_model":physical.covalently_closed_in_model})),
                "verification":{"status":preview.verification.status,"engine_version":preview.verification.engine_version,
                    "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                    "findings_truncated":preview.verification.findings.len()>100,
                    "assumptions":preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),
                    "assumptions_truncated":preview.verification.assumptions.len()>100,
                    "exhaustive_alternative_search":preview.verification.exhaustive_alternative_search},
                "artifact_id":artifact_id,
                "authority":"Fresh deterministic strand/end/junction reconstruction under explicit preparation and terminal-treatment assumptions. Saving remains a candidate; sequence agreement is not laboratory success or exhaustive mixture uniqueness. Full physical components, junctions and nick evidence are retained in the artifact."});
            let response = output(result)?;
            let mut store=server.session.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id) {
                if !matches!(existing.data,ArtifactData::Text(_)) || existing.mime_type!="application/json" || existing.filename!="ligation-assembly-recipe.json" {
                    return Err(bad("assembly evidence already has another artifact representation"));
                }
            }
            if (new_plan && store.plans.len()>=MAX_ITEMS) || (new_artifact && store.artifacts.len()>=MAX_ITEMS) {
                return Err(bad("candidate/artifact item limit; no assembly outputs created"));
            }
            let bytes=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("session storage accounting overflow"))?;
            store.reserve(bytes)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:preview.compiled.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"ligation-assembly-recipe.json"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Compile explicit source-bound TA/GC cloning against an independently supplied target. Each row retains actual prepared ends or explicitly adds one A/T/G/C nucleotide at both 3′ ends, with separate optional end treatment and intended ligation order/orientation. The selected TA/GC profile is assessed separately from generic ligation verification: blunt ligation success does not establish a matching TA profile. Saved origins are freshly reconstructed, retaining original primers, modifications and incomplete evidence. No implicit polymerase, kit, phosphate, host repair, probabilities or experimental isolation is assumed. Retains an immutable candidate and full evidence; verify_plan and HTML/PDF export use its plan_revision."
    )]
    pub(crate) async fn propose_tailing_assembly(
        &self,
        Parameters(args): Parameters<TailingAssemblyProposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad("tailing recipe requires 1..128 explicit source inputs"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.session.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len()>128 || !seen.insert(bound.id.as_str()) {
                        return Err(bad("empty, oversized or duplicate recipe input ID"));
                    }
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::cloning::ligation::tailing_recipe::source_preflight(&sources,target,&args.recipe).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,molecule)|PlanInput {
                    id:id.into(),revision:molecule.revision(),molecule:molecule.clone()
                }).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let preview = collagene_core::cloning::ligation::tailing_recipe::preview(&inputs,&target,&args.recipe).map_err(bad)?;
            let expected = collagene_core::cloning::ligation::tailing_recipe::preview_revision(&inputs,&target,&args.recipe).map_err(bad)?;
            preview.validate_for(&inputs,&target,&args.recipe,&expected).map_err(bad)?;
            // Serialize bounded complete evidence before publishing either output.
            let text = bounded_json_limit(&preview,32*1024*1024)?;
            let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
            let revision = collagene_core::planning::plan::DesignPlan::revision(&preview.compiled.plan);
            if revision.is_empty() {return Err(bad("compiled plan has no bounded canonical revision"));}
            let plan_bytes = bounded_json(&preview.compiled.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("candidate storage accounting overflow"))?;
            let result = json!({"plan_revision":revision,"state":"candidate","verified":false,
                "target":args.target,"inputs":args.inputs,"context_revision":preview.context_revision,
                "profile_status":preview.profile_status,"profile_findings":preview.profile_findings,
                "preparation_complete":preview.preparation_complete,"ligation_step":preview.compiled.ligation_step,
                "fragment_count":preview.compiled.rows.len(),"junction_count":preview.junctions.len(),
                "product_physical":preview.product_physical.as_ref().map(|physical|json!({"top_length":physical.top_length,"bottom_length":physical.bottom_length,"topology":physical.topology,"nick_count":physical.nicks.len(),"covalently_closed_in_model":physical.covalently_closed_in_model})),
                "verification":{"status":preview.verification.status,"engine_version":preview.verification.engine_version,
                    "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                    "findings_truncated":preview.verification.findings.len()>100,
                    "assumptions":preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),
                    "assumptions_truncated":preview.verification.assumptions.len()>100,
                    "exhaustive_alternative_search":preview.verification.exhaustive_alternative_search},
                "artifact_id":artifact_id,
                "authority":"Fresh deterministic strand/end/junction reconstruction under explicit preparation and terminal-treatment assumptions. Saving remains a candidate; profile matching does not replace ligation verification, and sequence agreement is not laboratory success or exhaustive mixture uniqueness. Full physical components, junctions and nick evidence are retained in the artifact."});
            let response = output(result)?;
            let mut store=server.session.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id) {
                if !matches!(existing.data,ArtifactData::Text(_)) || existing.mime_type!="application/json" || existing.filename!="tailing-assembly-recipe.json" {
                    return Err(bad("assembly evidence already has another artifact representation"));
                }
            }
            if (new_plan && store.plans.len()>=MAX_ITEMS) || (new_artifact && store.artifacts.len()>=MAX_ITEMS) {
                return Err(bad("candidate/artifact item limit; no assembly outputs created"));
            }
            let bytes=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("session storage accounting overflow"))?;
            store.reserve(bytes)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:preview.compiled.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"tailing-assembly-recipe.json"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Inspect finite Golden Gate alternatives from freshly reconstructed explicit preparations and an independent target. Reuses actual private PCR/digest/treatment strand state; no serialized chemistry summary supplies a substrate. Redigests the prepared pool under the stated complete-digest, one-copy, orientation and search limits. Includes subsets/self-background and recutting candidates, preserves incomplete search and original ordered findings, and withholds enumeration when required physical preparation is missing. Returns all compact candidate summaries and an immutable full report, never an automatic winner or exhaustive real-mixture/yield claim."
    )]
    pub(crate) async fn inspect_prepared_golden_gate_pool(
        &self,
        Parameters(args): Parameters<PreparedGoldenGatePoolQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 original pool source bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.session.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {
                    if input.id.trim().is_empty() || input.id.len()>128 || !seen.insert(input.id.as_str()) { return Err(bad("invalid or duplicate pool source ID")); }
                    sources.push((input.id.as_str(),store.document(&input.document)?));
                }
                collagene_core::cloning::golden_gate::prepared::preflight(&sources,target,&args.request).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())
            };
            server.session.attest_material_origins(&inputs,&args.origin_workspaces,&collagene_core::cloning::golden_gate::prepared::material_origins(&args.request))?;
            let report = collagene_core::cloning::golden_gate::prepared::inspect(&inputs,&target,&args.request).map_err(bad)?;
            let expected = collagene_core::cloning::golden_gate::prepared::inspect_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let candidates = report.assessments.iter().map(|a| {
                let c=report.pool.as_ref().and_then(|pool|pool.candidates.get(a.candidate_index)).ok_or_else(||bad("pool assessment has no bound candidate"))?;
                Ok(json!({"candidate_index":a.candidate_index,"directed_sequence_key":a.directed_sequence_key,"physical_sequence_key":a.physical_sequence_key,
                    "target_matches":a.target_matches,"source_roles":a.source_roles,"length":c.duplex.top.len(),"topology":c.duplex.topology,
                    "restriction_site_free":c.restriction_site_free,"recut_sites":c.recut_sites.len(),"recut_complete":c.cleavage.complete,"recut_products_determined":c.cleavage.products_determined}))
            }).collect::<Result<Vec<_>,ErrorData>>()?;
            #[derive(Serialize)]
            struct Evidence<'a> {kind:&'static str,inputs:&'a [BoundInput],target:&'a DocumentRef,origin_workspaces:&'a[MaterialOriginWorkspace],report:&'a collagene_core::cloning::golden_gate::prepared::Report}
            let text = bounded_json_limit(&Evidence {kind:"prepared-golden-gate-pool-v1",inputs:&args.inputs,target:&args.target,origin_workspaces:&args.origin_workspaces,report:&report},32*1024*1024)?;
            let artifact = server.session.text_artifact(text,"application/json","prepared-golden-gate-pool.json")?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,
                "plan_revision":report.preview.verification.plan_revision,"target_revision":report.preview.verification.target_revision,
                "inputs":args.inputs,"target":args.target,"preparation_complete":report.preview.preparation_complete,"missing_physical":report.missing_physical,"products_determined":report.pool.as_ref().map(|p|p.products_determined),"profile_revision":report.pool.as_ref().map(|p|p.profile_revision.clone()),"cleavage_comparisons":report.pool.as_ref().map(|p|p.comparisons),"cleavage_budget_scope":"Preparative checked cuts plus Golden Gate source/candidate scans in this native plan; saved-ancestry replay has separate bounded scopes.",
                "ordered_verification":{"status":report.preview.verification.status,"engine_version":report.preview.verification.engine_version,
                    "findings":report.preview.verification.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":report.preview.verification.findings.len()>100,
                    "assumptions":report.preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),"assumptions_truncated":report.preview.verification.assumptions.len()>100},
                "enumeration_complete":report.pool.as_ref().is_some_and(|p|p.enumeration_complete),"exhaustive_mixture_search":false,
                "states_examined":report.pool.as_ref().map_or(0,|p|p.states_examined),"truncation_reasons":report.pool.as_ref().map(|p|&p.truncation_reasons),
                "stable_physical_product_count":report.pool.as_ref().map_or(0,|p|p.stable_physical_product_count),"candidate_count":candidates.len(),"display_truncated":false,
                "roles":report.preview.compiled.rows.iter().map(|r|json!({"index":r.preparation.row_index,"role":r.preparation.role,"input_id":r.preparation.input_id,"reverse_in_ordered_recipe":r.preparation.reverse})).collect::<Vec<_>>(),
                "candidates":candidates,"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Fresh preparation and bounded one-copy pool inspection only; no automatic candidate, laboratory outcome or exhaustive mixture claim. Select an exact index and both keys, then freshly verify its plan."}))
        }).await
    }
    #[tool(
        description = "Select one exact site-free candidate from an immutable prepared Golden Gate pool artifact by context, index and both sequence keys. Preserves original inputs and every preparation/treatment step, replacing only final ordered ligation with the existing GoldenGate operation. Atomically saves an unverified plan and binary cg2plan artifact. Does not trust serialized physical state, install product DNA or assert target agreement. verify_plan must freshly reconstruct the selected plan against the same independent target."
    )]
    pub(crate) async fn select_prepared_golden_gate_candidate(
        &self,
        Parameters(args): Parameters<PreparedGoldenGateSelectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (evidence,target) = {
                let store = server.session.store.lock().unwrap();
                let artifact=store.artifacts.get(&args.artifact_id).ok_or_else(||bad("unknown prepared pool artifact"))?;
                if artifact.filename!="prepared-golden-gate-pool.json" || artifact.mime_type!="application/json" { return Err(bad("artifact is not prepared pool evidence")); }
                let ArtifactData::Text(text)=&artifact.data else {return Err(bad("prepared pool evidence must be JSON"));};
                if text.len()>32*1024*1024 {return Err(bad("prepared pool evidence exceeds32MiB"));}
                let evidence: PreparedGoldenGateEvidence=serde_json::from_str(text).map_err(bad)?;
                if evidence.kind!="prepared-golden-gate-pool-v1" || evidence.inputs.len()!=evidence.report.preview.compiled.plan.inputs.len() { return Err(bad("prepared pool source context changed")); }
                evidence.report.validate_binding(&args.selection.context_revision).map_err(bad)?;
                let target=store.document(&evidence.target)?;
                let mut sources=Vec::with_capacity(evidence.inputs.len());
                for (bound,input) in evidence.inputs.iter().zip(&evidence.report.preview.compiled.plan.inputs) {
                    let source=store.document(&bound.document)?;
                    if input.id!=bound.id || input.revision!=bound.document.revision || source.revision()!=input.revision { return Err(bad("prepared pool original source differs from retained evidence")); }
                    sources.push((bound.id.as_str(),source));
                }
                collagene_core::cloning::golden_gate::prepared::preflight(&sources,target,&evidence.report.request).map_err(bad)?;
                if target.revision()!=evidence.report.preview.verification.target_revision { return Err(bad("prepared pool independent target changed")); }
                let target=target.clone();
                (evidence,target)
            };
            server.session.attest_material_origins(&evidence.report.preview.compiled.plan.inputs,&evidence.origin_workspaces,&collagene_core::cloning::golden_gate::prepared::material_origins(&evidence.report.request))?;
            let selected=collagene_core::cloning::golden_gate::prepared::select(&evidence.report,&args.selection).map_err(bad)?;
            let expected=collagene_core::cloning::golden_gate::prepared::selection_revision(&evidence.report,&args.selection).map_err(bad)?;
            selected.validate_binding(&expected).map_err(bad)?;
            let revision=collagene_core::planning::plan::DesignPlan::revision(&selected.plan);
            let bytes=transfer::save_candidate(&transfer::Candidate::new(selected.plan.clone())).map_err(bad)?;
            let artifact_id=blake3::hash(&bytes).to_hex().to_string();
            let plan_bytes=bounded_json(&selected.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("selected candidate storage overflow"))?;
            let response=output(json!({"state":"candidate","verified":false,"plan_revision":revision,"context_revision":selected.context_revision,"parent_context":selected.parent_context,
                "target":evidence.target,"inputs":evidence.inputs,"selection":selected.selection,"warnings":selected.warnings,
                "candidate_artifact_id":artifact_id,"artifact_id":artifact_id,"filename":"candidate.cg2plan","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len()}))?;
            let mut store=server.session.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(old)=store.artifacts.get(&artifact_id) {
                if !matches!(old.data,ArtifactData::Binary(_)) || old.filename!="candidate.cg2plan" || old.mime_type!="application/octet-stream" { return Err(bad("candidate bytes already have another artifact representation")); }
            }
            if new_plan&&store.plans.len()>=MAX_ITEMS || new_artifact&&store.artifacts.len()>=MAX_ITEMS {return Err(bad("candidate/artifact item limit; no pool selection created"));}
            let total=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact {bytes.len()}else{0}).ok_or_else(||bad("selected candidate storage overflow"))?;
            store.reserve(total)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:selected.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Binary(bytes),mime_type:"application/octet-stream",encoding:"base64",filename:"candidate.cg2plan"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Inspect an explicitly declared Golden Gate reaction from whole original circular/linear DNA, complete PCR amplicons or selected preparations. Checked positional cleavage includes all established fragments of every declared substrate and reports proven intact circular carryover separately. Unknown or impaired source cleavage withholds the entire requested pool; candidate recut evidence retains uncertainty. Actual preparation strands and enumeration come from one native pass; metadata cannot supply physical state. Retains recuttable products, backgrounds and search limits. Index0 probe diagnostics do not nominate a candidate or establish kinetics, repeated-copy stoichiometry, methylation clearance or experimental yield."
    )]
    pub(crate) async fn inspect_golden_gate_reaction(
        &self,
        Parameters(args): Parameters<GoldenGateReactionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 original pool source bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.session.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {
                    if input.id.trim().is_empty() || input.id.len()>128 || !seen.insert(input.id.as_str()) { return Err(bad("invalid or duplicate pool source ID")); }
                    sources.push((input.id.as_str(),store.document(&input.document)?));
                }
                collagene_core::cloning::golden_gate::reaction::preflight(&sources,target,&args.request).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())
            };
            server.session.attest_material_origins(&inputs,&args.origin_workspaces,&collagene_core::cloning::golden_gate::reaction::material_origins(&args.request))?;
            let report = collagene_core::cloning::golden_gate::reaction::inspect(&inputs,&target,&args.request).map_err(bad)?;
            let expected = collagene_core::cloning::golden_gate::reaction::inspect_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let candidates = report.assessments.iter().map(|a| {
                let c=report.pool.as_ref().and_then(|pool|pool.candidates.get(a.candidate_index)).ok_or_else(||bad("pool assessment has no bound candidate"))?;
                Ok(json!({"candidate_index":a.candidate_index,"directed_sequence_key":a.directed_sequence_key,"physical_sequence_key":a.physical_sequence_key,
                    "target_matches":a.target_matches,"source_roles":a.source_roles,"length":c.duplex.top.len(),"topology":c.duplex.topology,
                    "restriction_site_free":c.restriction_site_free,"recut_sites":c.recut_sites.len(),"recut_complete":c.cleavage.complete,"recut_products_determined":c.cleavage.products_determined}))
            }).collect::<Result<Vec<_>,ErrorData>>()?;
            #[derive(Serialize)]
            struct Evidence<'a> {kind:&'static str,inputs:&'a [BoundInput],target:&'a DocumentRef,origin_workspaces:&'a[MaterialOriginWorkspace],report:&'a collagene_core::cloning::golden_gate::reaction::Report}
            let text = bounded_json_limit(&Evidence {kind:"golden-gate-reaction-v1",inputs:&args.inputs,target:&args.target,origin_workspaces:&args.origin_workspaces,report:&report},32*1024*1024)?;
            let artifact = server.session.text_artifact(text,"application/json","golden-gate-reaction.json")?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,
                "probe_plan_revision":report.verification.plan_revision,"target_revision":report.verification.target_revision,
                "inputs":args.inputs,"target":args.target,"preparation_complete":report.verification.steps.iter().take(report.probe_plan.steps.len().saturating_sub(1)).all(|step|step.state==verification::evidence::StepState::Reconstructed&&step.product.is_some()),"pool_available":report.pool.is_some(),"products_determined":report.pool.as_ref().map(|p|p.products_determined),"profile_revision":report.pool.as_ref().map(|p|p.profile_revision.clone()),"cleavage_comparisons":report.pool.as_ref().map(|p|p.comparisons),"cleavage_budget_scope":"Preparative checked cuts plus Golden Gate source/candidate scans in this native plan; saved-ancestry replay has separate bounded scopes.",
                "probe_verification":{"status":report.verification.status,"engine_version":report.verification.engine_version,
                    "findings":report.verification.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":report.verification.findings.len()>100,
                    "assumptions":report.verification.assumptions.iter().take(100).collect::<Vec<_>>(),"assumptions_truncated":report.verification.assumptions.len()>100},
                "enumeration_complete":report.pool.as_ref().is_some_and(|p|p.enumeration_complete),"exhaustive_mixture_search":false,
                "states_examined":report.pool.as_ref().map_or(0,|p|p.states_examined),"truncation_reasons":report.pool.as_ref().map(|p|&p.truncation_reasons),
                "stable_physical_product_count":report.pool.as_ref().map_or(0,|p|p.stable_physical_product_count),"candidate_count":candidates.len(),"display_truncated":false,
                "roles":report.materials.iter().map(|r|json!({"index":r.index,"role":r.role,"input_id":r.input_id,"source":r.source})).collect::<Vec<_>>(),
                "carryover_sources":report.pool.as_ref().map(|p|&p.carryover_sources),
                "candidates":candidates,"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Fresh declared-reaction inspection only; the probe selects no saved winner and establishes no laboratory outcome or exhaustive mixture claim. Select an exact index and both keys, then freshly verify its plan."}))
        }).await
    }
    #[tool(
        description = "Select one exact site-free candidate from an immutable whole-reaction Golden Gate artifact by context, index and both keys. Preserves original inputs, all explicit preparations and the declared reaction substrates; changes only the final GoldenGate candidate index. Atomically saves an unverified plan and binary cg2plan artifact. Does not trust serialized physical state, install product DNA or assert target agreement. verify_plan must freshly reconstruct the selected plan against the same independent target."
    )]
    pub(crate) async fn select_golden_gate_reaction_candidate(
        &self,
        Parameters(args): Parameters<PreparedGoldenGateSelectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (evidence,target) = {
                let store = server.session.store.lock().unwrap();
                let artifact=store.artifacts.get(&args.artifact_id).ok_or_else(||bad("unknown prepared pool artifact"))?;
                if artifact.filename!="golden-gate-reaction.json" || artifact.mime_type!="application/json" { return Err(bad("artifact is not prepared pool evidence")); }
                let ArtifactData::Text(text)=&artifact.data else {return Err(bad("prepared pool evidence must be JSON"));};
                if text.len()>32*1024*1024 {return Err(bad("prepared pool evidence exceeds32MiB"));}
                let evidence: GoldenGateReactionEvidence=serde_json::from_str(text).map_err(bad)?;
                if evidence.kind!="golden-gate-reaction-v1" || evidence.inputs.len()!=evidence.report.probe_plan.inputs.len() { return Err(bad("prepared pool source context changed")); }
                evidence.report.validate_binding(&args.selection.context_revision).map_err(bad)?;
                let target=store.document(&evidence.target)?;
                let mut sources=Vec::with_capacity(evidence.inputs.len());
                for (bound,input) in evidence.inputs.iter().zip(&evidence.report.probe_plan.inputs) {
                    let source=store.document(&bound.document)?;
                    if input.id!=bound.id || input.revision!=bound.document.revision || source.revision()!=input.revision { return Err(bad("prepared pool original source differs from retained evidence")); }
                    sources.push((bound.id.as_str(),source));
                }
                collagene_core::cloning::golden_gate::reaction::preflight(&sources,target,&evidence.report.request).map_err(bad)?;
                if target.revision()!=evidence.report.verification.target_revision { return Err(bad("prepared pool independent target changed")); }
                let target=target.clone();
                (evidence,target)
            };
            server.session.attest_material_origins(&evidence.report.probe_plan.inputs,&evidence.origin_workspaces,&collagene_core::cloning::golden_gate::reaction::material_origins(&evidence.report.request))?;
            let selected=collagene_core::cloning::golden_gate::reaction::select(&evidence.report,&args.selection).map_err(bad)?;
            let expected=collagene_core::cloning::golden_gate::reaction::selection_revision(&evidence.report,&args.selection).map_err(bad)?;
            selected.validate_binding(&expected).map_err(bad)?;
            let revision=collagene_core::planning::plan::DesignPlan::revision(&selected.plan);
            let bytes=transfer::save_candidate(&transfer::Candidate::new(selected.plan.clone())).map_err(bad)?;
            let artifact_id=blake3::hash(&bytes).to_hex().to_string();
            let plan_bytes=bounded_json(&selected.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("selected candidate storage overflow"))?;
            let response=output(json!({"state":"candidate","verified":false,"plan_revision":revision,"context_revision":selected.context_revision,"parent_context":selected.parent_context,
                "target":evidence.target,"inputs":evidence.inputs,"selection":selected.selection,"warnings":selected.warnings,
                "candidate_artifact_id":artifact_id,"artifact_id":artifact_id,"filename":"candidate.cg2plan","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len()}))?;
            let mut store=server.session.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(old)=store.artifacts.get(&artifact_id) {
                if !matches!(old.data,ArtifactData::Binary(_)) || old.filename!="candidate.cg2plan" || old.mime_type!="application/octet-stream" { return Err(bad("candidate bytes already have another artifact representation")); }
            }
            if new_plan&&store.plans.len()>=MAX_ITEMS || new_artifact&&store.artifacts.len()>=MAX_ITEMS {return Err(bad("candidate/artifact item limit; no pool selection created"));}
            let total=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact {bytes.len()}else{0}).ok_or_else(||bad("selected candidate storage overflow"))?;
            store.reserve(total)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:selected.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Binary(bytes),mime_type:"application/octet-stream",encoding:"base64",filename:"candidate.cg2plan"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Freshly inspect DNA junctions of an exact saved overlap, ordered ligation or selected Golden Gate candidate and its stored independent target. Golden Gate boundaries and source mappings come from its actual selected digest path, captured in the same verifier pass; recipe order and target DNA never determine them. Explicit bounded flanks, genetic code, six-frame choices and CDS inspection options. Reports actual product boundary/homology/orientation/circular-origin windows, boundary-crossing codons, conditional source-CDS evidence and physical ends/nicks where supported. Full target verification and unsupported final-method diagnostics remain visible. Immutable JSON evidence only: no source edits, candidate promotion, inferred fusion function or experimental guarantee."
    )]
    pub(crate) async fn inspect_assembly_junctions(
        &self,
        Parameters(args): Parameters<AssemblyJunctionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            // Borrow actual saved source/target snapshots. No source or plan copy
            // precedes the core's bounded native-model preflight.
            let (plan,target,expected)={let store=server.session.store.lock().unwrap();
                let saved=store.plans.get(&args.plan_revision).ok_or_else(||bad("unknown plan revision"))?;
                collagene_core::cloning::assembly::junctions::source_preflight(&saved.plan,&saved.target).map_err(bad)?;
                let expected=collagene_core::cloning::assembly::junctions::inspection_revision(&saved.plan,&saved.target,&args.options).map_err(bad)?;
                (saved.plan.clone(),saved.target.clone(),expected)};
            let report=collagene_core::cloning::assembly::junctions::inspect_plan(&plan,&target,&args.options).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&report,32*1024*1024)?,"application/json","assembly-junctions.json")?;
            let junctions=report.junctions.iter().take(10).map(|junction|json!({
                "index":junction.index,"upstream":junction.upstream,"downstream":junction.downstream,
                "closing":junction.closing,"boundary":junction.boundary,"homology":junction.homology,
                "window":junction.window,
                "frames":junction.frames.iter().take(6).map(|frame|json!({"spec":frame.spec,
                    "crossing_codon":frame.crossing_codon,
                    "protein_preview":frame.translation.protein.chars().take(100).collect::<String>(),
                    "protein_truncated":frame.translation.protein.len()>100})).collect::<Vec<_>>(),
                "cds_count":junction.cds.len(),"cds_omitted":junction.cds_omitted,
                "cds":junction.cds.iter().take(10).map(|cds|json!({"feature_index":cds.feature_index,
                    "label":cds.label.chars().take(256).collect::<String>(),"label_truncated":cds.label.chars().count()>256,
                    "needs_review":cds.needs_review,"spans_boundary":cds.spans_boundary,
                    "warning":cds.warning.as_ref().map(|warning|warning.chars().take(1024).collect::<String>()),
                    "warning_truncated":cds.warning.as_ref().is_some_and(|warning|warning.chars().count()>1024)})).collect::<Vec<_>>(),
                "cds_preview_truncated":junction.cds.len()>10,
                "physical":junction.physical.as_ref().map(|physical|json!({"upstream_end":physical.upstream_end,"downstream_end":physical.downstream_end,"nick_count":physical.product_nicks.len()}))
            })).collect::<Vec<_>>();
            output(json!({"plan_revision":report.plan_revision,"target_revision":report.target_revision,
                "context_revision":report.context_revision,"engine_version":report.engine_version,"options":report.options,
                "verification_status":report.verification.status,"verification_findings":report.verification.findings.iter().take(100).collect::<Vec<_>>(),
                "verification_findings_truncated":report.verification.findings.len()>100,
                "product_revision":report.product_revision,"product_length":report.product_length,"product_topology":report.product_topology,
                "method":report.method,"complete":report.complete,"omitted_junctions":report.omitted_junctions,
                "junction_count":report.junctions.len(),"junctions":junctions,"junction_preview_truncated":report.junctions.len()>10,
                "warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100 || report.warnings.iter().take(100).any(|w|w.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Fresh declared-model junction DNA and conditional annotation/frame evidence only. Local codons do not prove a functional fusion or alter full target verification; saved evidence never promotes this candidate or repairs physical strands."}))
        }).await
    }
}
