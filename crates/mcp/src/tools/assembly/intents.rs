use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::MAX_ITEMS;
use crate::session::MAX_STORE_BYTES;
use crate::session::StoredPlan;
use crate::session::references::DocumentRef;
use crate::session::references::PlanRef;
use crate::session::references::WorkspaceRef;
use collagene_core::documents::transfer;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::planning::plan::Operation;
use collagene_core::planning::plan::PlanInput;
use collagene_core::verification;
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

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoundInput {
    /// Plan-local source ID used by steps, e.g. left or template.
    pub id: String,
    pub document: DocumentRef,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Proposal {
    pub name: String,
    /// Independent desired target. Never inferred from a candidate product.
    pub target: DocumentRef,
    pub inputs: Vec<BoundInput>,
    /// Typed operations refer to input IDs or earlier steps. OverlapExtension explicitly denatures/reanneals two intact linear tracked duplexes and extends their selected terminal pairing; use actual PCR oligos before and after it. Reference overlap assemblies cannot substitute for physical PCR templates.
    pub steps: Vec<Operation>,
}
/// MCP accepts small exact native selections, never a caller's closed capsule.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IntentMaterialBinding {
    pub input_id: String,
    pub event: collagene_core::material::lineage::EventBinding,
    pub selection: collagene_core::material::run::MaterialSelection,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssemblyIntentRequest {
    pub name: String,
    pub fragments: Vec<collagene_core::cloning::assembly::intent::Contribution>,
    pub circular: bool,
    pub material_bindings: Vec<IntentMaterialBinding>,
}
impl AssemblyIntentRequest {
    pub(crate) fn native(
        &self,
        workspace: Option<&collagene_core::documents::workspace::Workspace>,
    ) -> Result<collagene_core::cloning::assembly::intent::Request, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        for f in &self.fragments {
            if f.role.trim().is_empty()
                || f.role.len() > 128
                || f.input_id.trim().is_empty()
                || f.input_id.len() > 128
            {
                return Err(bad("bounded source contribution role/input required"));
            }
        }
        if self.fragments.is_empty()
            || self.fragments.len() > 32
            || self.material_bindings.len() > 32
            || self.name.len() > 1024
        {
            return Err(bad(
                "bounded explicit contributions and material bindings required",
            ));
        }
        let mut ids = BTreeSet::new();
        for b in &self.material_bindings {
            if b.input_id.trim().is_empty()
                || b.input_id.len() > 128
                || !ids.insert(b.input_id.as_str())
                || b.selection.event() != &b.event
            {
                return Err(bad(
                    "material selection must name its exact unique input/event",
                ));
            }
            b.selection.validate().map_err(bad)?;
        }
        let origins = if self.material_bindings.is_empty() {
            vec![]
        } else {
            let workspace = workspace
                .ok_or_else(|| bad("native material selection requires origin_workspace"))?;
            let selections = self
                .material_bindings
                .iter()
                .map(|b| b.selection.clone())
                .collect::<Vec<_>>();
            collagene_core::material::capsule::documentary_origins(workspace, &selections, &[])
                .map_err(bad)?
                .materials
        };
        Ok(collagene_core::cloning::assembly::intent::Request {
            name: self.name.clone(),
            fragments: self.fragments.clone(),
            circular: self.circular,
            material_bindings: self
                .material_bindings
                .iter()
                .zip(origins)
                .map(|(binding, origin)| {
                    collagene_core::cloning::assembly::intent::MaterialBinding {
                        input_id: binding.input_id.clone(),
                        event: binding.event.clone(),
                        origin,
                    }
                })
                .collect(),
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SubstitutedIntentRequest {
    pub base: AssemblyIntentRequest,
    pub junction_edits: Vec<collagene_core::cloning::assembly::junction_edit::Edit>,
    pub substitutions: Vec<collagene_core::cloning::assembly::substitution_intent::Substitution>,
}
impl SubstitutedIntentRequest {
    pub(crate) fn native(
        &self,
        w: Option<&collagene_core::documents::workspace::Workspace>,
    ) -> Result<collagene_core::cloning::assembly::substitution_intent::Request, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        if self.junction_edits.len() > 32 || self.substitutions.len() > 256 {
            return Err(bad("intent edit bound"));
        }
        Ok(
            collagene_core::cloning::assembly::substitution_intent::Request {
                base: self.base.native(w)?,
                junction_edits: self.junction_edits.clone(),
                substitutions: self.substitutions.clone(),
            },
        )
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContributionEditIntentRequest {
    pub base: SubstitutedIntentRequest,
    pub edits: Vec<collagene_core::cloning::assembly::edit_intent::Edit>,
}
impl ContributionEditIntentRequest {
    pub(crate) fn native(
        &self,
        w: Option<&collagene_core::documents::workspace::Workspace>,
    ) -> Result<collagene_core::cloning::assembly::edit_intent::Request, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        if self.edits.len() > 256 {
            return Err(bad("intent edit bound"));
        }
        Ok(collagene_core::cloning::assembly::edit_intent::Request {
            base: self.base.native(w)?,
            edits: self.edits.clone(),
        })
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssemblyIntentQuery {
    /// Exact original source snapshots; local IDs match request fragments.
    pub inputs: Vec<BoundInput>,
    /// Declared source regions, order and orientation only. No simulator-derived target.
    pub request: AssemblyIntentRequest,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedOverlapPrimerRequest {
    pub intent: ContributionEditIntentRequest,
    pub overlap: usize,
    pub options: collagene_core::cloning::overlap::extension::primers::Options,
    pub profile_options: collagene_core::cloning::overlap::profiles::ProfileOptions,
    pub finishing: collagene_core::cloning::assembly::recipe::OverlapFinishing,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SubstitutedAssemblyIntentQuery {
    /// Exact immutable MCP source snapshots, including retained historical revisions.
    pub inputs: Vec<BoundInput>,
    /// Independent regions, junction edits and expected-base substitutions; never primer-derived.
    pub request: SubstitutedIntentRequest,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContributionEditIntentQuery {
    /// Exact immutable original source snapshots, never products or guessed native parents.
    pub inputs: Vec<BoundInput>,
    /// Independent construct-oriented original contribution intervals and literal edits.
    pub request: ContributionEditIntentRequest,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedAssemblyIntentQuery {
    pub inputs: Vec<BoundInput>,
    pub base: AssemblyIntentRequest,
    pub edits: Vec<collagene_core::cloning::assembly::junction_edit::Edit>,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedAssemblyPrimersQuery {
    pub inputs: Vec<BoundInput>,
    pub request: collagene_core::cloning::assembly::junction_edit::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SelectiveAssemblyPrimersQuery {
    pub inputs: Vec<BoundInput>,
    pub request: collagene_core::cloning::assembly::selective_edit::Request,
}
#[tool_router(router = assembly_intents_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Export an exact saved candidate as a binary .cg2plan artifact for browser review. Contains no trusted verification status and never changes open browser documents. Read all base64 chunks with read_artifact, concatenate their text, then decode to file bytes."
    )]
    pub(crate) async fn export_candidate(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let saved=server.session.store.lock().unwrap().plans.get(&args.plan_revision).cloned().ok_or_else(||bad("unknown plan revision"))?;
            let bytes=transfer::save_candidate(&transfer::Candidate::new(saved.plan)).map_err(bad)?;
            let id=server.session.binary_artifact(bytes.clone(),"application/octet-stream","candidate.cg2plan")?;
            output(json!({"artifact_id":id,"filename":"candidate.cg2plan","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len(),"plan_revision":args.plan_revision,"state":"candidate"}))
        }).await
    }
    #[tool(
        description = "Build an intended reference from explicit regions of actual immutable source documents, their order, orientation and desired topology. No pre-existing target is required and no PCR, digest, ligation or assembly simulator supplies expected DNA. Atomically retains an unverified intended-reference target document and full source-placement JSON. Selection is design intent, not evidence a preparation can produce it or chemical/experimental material. Pass returned target DocumentRef with the same inputs to design_overlap_primers and propose_overlap_assembly for separate fresh preparation and verification."
    )]
    pub(crate) async fn build_assembly_intent(
        &self,
        Parameters(args): Parameters<AssemblyIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad("assembly intent requires 1..128 exact source bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.session.restore_workspace(r)).transpose()?;
            let native=args.request.native(origin_workspace.as_ref())?;

            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()) {return Err(bad("empty, oversized or duplicate intent source ID"));}
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::cloning::assembly::intent::preflight(&sources,&native,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,molecule)|PlanInput{id:id.into(),revision:molecule.revision(),molecule:molecule.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::cloning::assembly::intent::build(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::cloning::assembly::intent::intent_revision(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{origin_workspace:&'a Option<WorkspaceRef>,inputs:&'a [BoundInput],report:&'a collagene_core::cloning::assembly::intent::Report,interpretation:&'static str}
            let evidence=Evidence{origin_workspace:&args.origin_workspace,inputs:&args.inputs,report:&report,
                interpretation:"Declared source-region/order/orientation intent. Saved source strands are freshly reconstructed from exact native origins; the requested new preparation is never simulated to choose expected DNA. This target is not prepared material, an observed sequence or an independently reconstructed cloning product. Verify actual molecular steps separately against this frozen intent."};
            let (artifact,target_document)=server.session.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"assembly-intent.json",32*1024*1024)?;
            let target_document=target_document.ok_or_else(||bad("intended-reference document was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,
                "inputs":args.inputs,"context_revision":report.context_revision,"request":args.request,
                "target":{"id":target_document["id"],"revision":target_document["revision"]},
                "target_summary":target_document,"intended_length":report.target.sequence.len(),"topology":report.target.topology,
                "fragment_count":report.fragments.len(),"fragments":report.fragments,
                "warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|w|w.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Unverified intended-reference document only. Original source snapshots are unchanged. Region selection does not imply PCR amplification, restriction cuts, compatible ends, sealing, chemistry or experimental success; actual preparation must be reconstructed separately."}))
        }).await
    }
    #[tool(
        description = "Build an independently declared assembly reference from exact immutable source snapshots, explicit region/order/orientation, junction trims/inserts and expected-base role-relative substitutions. Mutation intent comes from the declaration, never primers or simulated products. Checks original expected bases before copying/building and atomically retains an unverified intended-reference session document plus complete mapping JSON. Repeated source contributions retain distinct roles. Original/historical source snapshots remain unchanged. Material selectors are resolved from explicit native event/state bindings in origin_workspace; reference-only rows remain independent source-document declarations. No new preparation, material or successful verification is created."
    )]
    pub(crate) async fn build_substituted_assembly_intent(
        &self,
        Parameters(args): Parameters<SubstitutedAssemblyIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad(
                "substituted intent requires1..32 exact original source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.session.restore_workspace(r)).transpose()?;
            let native=args.request.native(origin_workspace.as_ref())?;

            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs{if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()){return Err(bad("invalid or duplicate original substitution source ID"));}sources.push((bound.id.as_str(),store.document(&bound.document)?));}
                collagene_core::cloning::assembly::substitution_intent::source_preflight(&sources,&native,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::cloning::assembly::substitution_intent::build(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::cloning::assembly::substitution_intent::intent_revision(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]struct Evidence<'a>{origin_workspace:&'a Option<WorkspaceRef>,inputs:&'a[BoundInput],report:&'a collagene_core::cloning::assembly::substitution_intent::Report,binding_scope:&'static str}
            let scope="Exact immutable source roles/revisions; material rows additionally attest explicit native events/state through the pinned origin_workspace. Reference-only rows declare correspondence, never material ancestry.";
            let evidence=Evidence{origin_workspace:&args.origin_workspace,inputs:&args.inputs,report:&report,binding_scope:scope};
            let(artifact,target)=server.session.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"substituted-assembly-intent.json",32*1024*1024)?;
            let target=target.ok_or_else(||bad("substituted intended reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,"inputs":args.inputs,"context_revision":report.context_revision,"request":args.request,
                "target":{"id":target["id"],"revision":target["revision"]},"target_summary":target,"intended_length":report.target.sequence.len(),"topology":report.target.topology,
                "fragments":report.fragments,"applied":report.applied,"warnings":report.warnings,"artifact_id":artifact,"binding_scope":scope,
                "authority":"Independent source-bound design requirement only; original sources are unchanged. No primer-derived target, physical preparation, experimental success or native event ancestry is inferred."}))
        }).await
    }
    #[tool(
        description = "Build an independently declared assembly reference from exact immutable source snapshots, explicit region/order/orientation, junction trims/inserts, expected-base substitutions and literal contribution insertions/deletions/replacements. Mutation intent comes from the declaration, never primers or simulated products. Checks original expected literal bases before copying/building and atomically retains an unverified intended-reference session document plus complete mapping JSON. Repeated source contributions retain distinct roles. Original/historical source snapshots remain unchanged. Material selectors are resolved from explicit native event/state bindings in origin_workspace; reference-only rows remain independent source-document declarations. No new preparation, material or successful verification is created."
    )]
    pub(crate) async fn build_contribution_edit_intent(
        &self,
        Parameters(args): Parameters<ContributionEditIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad(
                "contribution edit intent requires1..32 exact original source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.session.restore_workspace(r)).transpose()?;
            let native=args.request.native(origin_workspace.as_ref())?;

            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs{if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()){return Err(bad("invalid or duplicate original contribution edit source ID"));}sources.push((bound.id.as_str(),store.document(&bound.document)?));}
                collagene_core::cloning::assembly::edit_intent::source_preflight(&sources,&native,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::cloning::assembly::edit_intent::build(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::cloning::assembly::edit_intent::intent_revision(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            report.validate_sources(&inputs,&expected,origin_workspace.as_ref()).map_err(bad)?;
            #[derive(Serialize)]struct Evidence<'a>{origin_workspace:&'a Option<WorkspaceRef>,inputs:&'a[BoundInput],report:&'a collagene_core::cloning::assembly::edit_intent::Report,binding_scope:&'static str}
            let scope="Exact immutable source roles/revisions; material rows additionally attest explicit native events/state through the pinned origin_workspace. Reference-only rows declare correspondence, never material ancestry.";
            let evidence=Evidence{origin_workspace:&args.origin_workspace,inputs:&args.inputs,report:&report,binding_scope:scope};
            let(artifact,target)=server.session.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"contribution-edit-intent.json",32*1024*1024)?;
            let target=target.ok_or_else(||bad("contribution-edited intended reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,"inputs":args.inputs,"context_revision":report.context_revision,"request":args.request,
                "target":{"id":target["id"],"revision":target["revision"]},"target_summary":target,"intended_length":report.target.sequence.len(),"topology":report.target.topology,
                "fragments":report.fragments,"applied":report.applied,"warnings":report.warnings,"artifact_id":artifact,"binding_scope":scope,
                "authority":"Independent source-bound design requirement only; original sources are unchanged. No primer-derived target, physical preparation, experimental success or native event ancestry is inferred."}))
        }).await
    }
    #[tool(
        description = "Build a declared edited assembly reference from exact original source regions and explicit junction deletions/insertions, independently of PCR or assembly simulation. Each edit trims bases from the upstream oriented core and downstream oriented core, then inserts literal DNA at that junction. Retains original order/orientation and an explicit circular origin. Atomically returns an unverified intended target document and complete JSON. Source documents stay unchanged; physical preparation must be reconstructed separately."
    )]
    pub(crate) async fn build_edited_assembly_intent(
        &self,
        Parameters(args): Parameters<EditedAssemblyIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.session.restore_workspace(r)).transpose()?;
            let native=args.base.native(origin_workspace.as_ref())?;

            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::cloning::assembly::junction_edit::intent_preflight(&sources,&native,&args.edits,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let intent=collagene_core::cloning::assembly::junction_edit::build_intent(&inputs,&native,&args.edits,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::cloning::assembly::junction_edit::intent_revision(&inputs,&native,&args.edits,origin_workspace.as_ref()).map_err(bad)?;intent.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"origin_workspace":args.origin_workspace,"intent":intent});
            let (artifact,target)=server.session.candidate_document_artifact_bounded(&evidence,Some(intent.target.clone()),"edited-assembly-intent.json",32*1024*1024)?;
            let target=target.ok_or_else(||bad("edited intended reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,"inputs":args.inputs,
                "context_revision":intent.context_revision,"base":args.base,"edits":intent.edits,
                "target":{"id":target["id"],"revision":target["revision"]},"target_summary":target,"fragments":intent.fragments,
                "warnings":intent.warnings,"artifact_id":artifact,"authority":"Declared source-region edits only. Original sources are unchanged; primer production, assembly outcome and experimental success remain separate evidence."}))
        }).await
    }
    #[tool(
        description = "Propose original-source PCR primers for an explicitly junction-edited overlap assembly. Uses the declared trimmed source cores and literal inserted DNA to form the independent intended reference before designing primers. Added primer tails include the requested junction DNA plus next-core homology; full oligo metrics and existing filters are recomputed. Returns an advisory explicit-PCR recipe, no automatic application or verified cloning result. Reconstruct that recipe against the independently built edited target using propose_overlap_assembly."
    )]
    pub(crate) async fn design_edited_assembly_primers(
        &self,
        Parameters(args): Parameters<EditedAssemblyPrimersQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::cloning::assembly::junction_edit::preflight(&sources,&args.request).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::cloning::assembly::junction_edit::design(&inputs,&args.request).map_err(bad)?;
            let expected=collagene_core::cloning::assembly::junction_edit::design_revision(&inputs,&args.request).map_err(bad)?;report.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"report":report});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","edited-assembly-primers.json")?;
            let ordered=report.recipe.fragments.iter().map(|f|match &f.preparation{collagene_core::cloning::assembly::recipe::Preparation::ExactPcr{forward,reverse, .. }=>Ok(json!({"role":f.role,"input_id":f.input_id,"reverse_fragment":f.reverse,"forward":forward,"reverse":reverse,"forward_ordered":forward.ordered_sequence(),"reverse_ordered":reverse.ordered_sequence()})),_=>Err(bad("edited primer proposal lost its exact original-source PCR"))}).collect::<Result<Vec<_>,_>>()?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,"intent_context_revision":report.intent.context_revision,
                "target_revision":report.intent.target.revision(),"recipe":report.recipe,"ordered_oligos":ordered,"design_count":report.designs.len(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Propose original-source PCR primers for explicitly selected rows of a junction-edited overlap assembly, preserving all retained preparation DTOs exactly, including manual primers and pinned digests. Retained preparations may be incompatible: fresh full-recipe verification is required. Uses the declared trimmed source cores and literal inserted DNA to form the independent intended reference before designing primers. Added primer tails include the requested junction DNA plus next-core homology; full oligo metrics and existing filters are recomputed. Returns an advisory explicit-PCR recipe, no automatic application or verified cloning result. Reconstruct that recipe against the independently built edited target using propose_overlap_assembly."
    )]
    pub(crate) async fn design_selective_assembly_primers(
        &self,
        Parameters(args): Parameters<SelectiveAssemblyPrimersQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::cloning::assembly::selective_edit::preflight(&sources,&args.request).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::cloning::assembly::selective_edit::design(&inputs,&args.request).map_err(bad)?;
            let expected=collagene_core::cloning::assembly::selective_edit::design_revision(&inputs,&args.request).map_err(bad)?;report.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"report":report});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","selective-assembly-primers.json")?;
            let ordered=report.recipe.fragments.iter().zip(&report.request.rows).filter_map(|(f,row)|match row {
                collagene_core::cloning::assembly::selective_edit::Row::Retain(_) => None,
                collagene_core::cloning::assembly::selective_edit::Row::DesignPcr{..} => Some(match &f.preparation {
                    collagene_core::cloning::assembly::recipe::Preparation::ExactPcr{forward,reverse, .. }=>Ok(json!({"role":f.role,"input_id":f.input_id,"reverse_fragment":f.reverse,"forward":forward,"reverse":reverse,"forward_ordered":forward.ordered_sequence(),"reverse_ordered":reverse.ordered_sequence()})),
                    _=>Err(bad("selected primer proposal lost its exact original-source PCR"))
                })
            }).collect::<Result<Vec<_>,_>>()?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,"intent_context_revision":report.intent.context_revision,
                "target_revision":report.intent.target.revision(),"recipe":report.recipe,"ordered_oligos":ordered,"design_count":ordered.len(),"retained_rows":report.retained_rows,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Save an immutable candidate plan binding an independent target and input revisions. DeclareSourceModifications accepts only an explicit Input source with Dam/Dcm/CpG generation assumptions; CheckedDigest pins an exact CatalogSelection and explicit inspection limits, consumes freshly reconstructed tracked strands and reports the fixed current positional profile. DeclareSourceTermini requires an original intact linear DNA Input and explicit left/right5-prime/3-prime assumptions; it cannot promote Step, ssDNA, circular, imported physical or reference-only states, infer chemistry from metadata or repair a graph. OverlapExtension requires two freshly reconstructed intact linear duplexes, explicit orientations and terminal pairing/search options; it reciprocally extends A.top/B.bottom and can feed actual outer-primer PCR. A sequence-only Overlap/ProfiledOverlap cannot supply that physical template. Saving never verifies or applies a plan. Duplicate submissions are idempotent."
    )]
    pub(crate) fn propose_plan(
        &self,
        Parameters(args): Parameters<Proposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.name.trim().is_empty()
            || args.name.len() > 1024
            || args.steps.is_empty()
            || args.steps.len() > 100
            || args.inputs.len() > 128
            || serde_json::to_vec(&args).map_err(bad)?.len() > 2 * 1024 * 1024
        {
            return Err(bad(
                "plan requires a name, 1–100 steps, at most 128 inputs and at most 2 MiB",
            ));
        }
        let mut store = self.session.store.lock().unwrap();
        let mut estimated = serde_json::to_vec(&args).map_err(bad)?.len();
        estimated = estimated.saturating_add(
            serde_json::to_vec(store.document(&args.target)?)
                .map_err(bad)?
                .len(),
        );
        for input in &args.inputs {
            estimated = estimated.saturating_add(
                serde_json::to_vec(store.document(&input.document)?)
                    .map_err(bad)?
                    .len(),
            );
            if estimated > MAX_STORE_BYTES {
                return Err(bad(
                    "candidate input materialization exceeds session storage budget",
                ));
            }
        }
        let target = store.document(&args.target)?.clone();
        let mut seen = std::collections::BTreeSet::new();
        let mut inputs = Vec::new();
        for input in args.inputs {
            if input.id.trim().is_empty() || input.id.len() > 128 || !seen.insert(input.id.clone())
            {
                return Err(bad("empty or duplicate input ID"));
            }
            inputs.push(PlanInput {
                id: input.id,
                revision: input.document.revision.clone(),
                molecule: store.document(&input.document)?.clone(),
            });
        }
        let plan = DesignPlan {
            name: args.name,
            target_revision: args.target.revision.clone(),
            inputs,
            steps: args.steps,
        };
        let revision = collagene_core::planning::plan::DesignPlan::revision(&plan);
        if !store.plans.contains_key(&revision) {
            if store.plans.len() >= MAX_ITEMS {
                return Err(bad("at most 128 candidate plans"));
            }
            let bytes = serde_json::to_vec(&plan).map_err(bad)?.len()
                + serde_json::to_vec(&target).map_err(bad)?.len();
            store.reserve(bytes)?;
            store
                .plans
                .insert(revision.clone(), StoredPlan { plan, target });
        }
        output(
            json!({"plan_revision":revision,"state":"candidate","verified":false,"target":args.target}),
        )
    }
    #[tool(
        description = "Reconstruct an exact saved plan revision against its separately bound target with the scientific core. Returns findings, assumptions and bounded product summaries. OverlapExtension freshly reconstructs the selected reciprocal extension; competing or unsearched terminal pairings remain explicit findings. Never trusts AI claimed products or status."
    )]
    pub(crate) async fn verify_plan(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.verify_plan_inner(args)).await
    }
    pub(crate) fn verify_plan_inner(&self, args: PlanRef) -> Result<CallToolResult, ErrorData> {
        let saved = self
            .session
            .store
            .lock()
            .unwrap()
            .plans
            .get(&args.plan_revision)
            .cloned()
            .ok_or_else(|| bad("unknown plan revision"))?;
        let report = verification::engine::verify(&saved.plan, &saved.target);
        let scientific_step_count = report.physical_processing_reports.len()
            + report.reference_step_reports.len()
            + report
                .step_reports
                .iter()
                .filter(|step| step.topo.is_some())
                .count();
        let reference_step_count = report.reference_step_reports.len();
        let artifact = if scientific_step_count > 0
            || reference_step_count > 0
            || !report.step_reports.is_empty()
        {
            Some(self.session.text_artifact(
                bounded_json_limit(&report, 32 * 1024 * 1024)?,
                "application/json",
                "verified-plan-scientific-evidence.json",
            )?)
        } else {
            None
        };
        let products:Vec<_>=report.products.iter().enumerate().map(|(i,m)|json!({"step":i,"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology,"duplex":m.metadata.get("duplex_model").and_then(|s|serde_json::from_str::<Value>(s).ok()),"tracked_modifications":m.metadata.get("tracked_modifications").and_then(|s|serde_json::from_str::<Value>(s).ok())})).collect();
        output(
            json!({"plan_revision":report.plan_revision,"target_revision":report.target_revision,"engine_version":report.engine_version,"status":report.status,"findings":report.findings,"assumptions":report.assumptions,"products":products,"exhaustive_alternative_search":report.exhaustive_alternative_search,"scientific_step_count":scientific_step_count,"reference_step_count":reference_step_count,"artifact_id":artifact}),
        )
    }
}
