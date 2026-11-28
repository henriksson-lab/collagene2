use crate::tools::default_inspection_page;

use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::Artifact;
use crate::session::ArtifactData;
use crate::session::MAX_ITEMS;
use crate::session::StoredPlan;
use crate::session::references::DocumentRef;
use crate::session::references::WorkspaceRef;
use crate::tools::assembly::intents::BoundInput;
use crate::tools::assembly::intents::EditedOverlapPrimerRequest;
use crate::tools::catalog::catalog_tools;

use collagene_core::chemistry;
use collagene_core::cloning::pcr::primers;
use collagene_core::documents::transfer;
use collagene_core::planning::plan::PlanInput;
use collagene_core::sequence::analysis;
use collagene_core::sequence::analysis::Primer;
use collagene_core::verification;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use rmcp::tool;
use rmcp::tool_router;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeSet;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SavedPrimerQuery {
    pub document: DocumentRef,
    /// Exact zero-based primer_bind feature index in this document revision.
    pub feature_index: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PrimerQuery {
    pub document: DocumentRef,
    pub primer: Primer,
    pub conditions: primers::ThermoConditions,
    pub binding_options: primers::BindingOptions,
    /// Optional second oligo for a contiguous complementarity screen.
    #[serde(default)]
    pub partner: Option<Primer>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PrimerDesignQuery {
    pub document: DocumentRef,
    pub options: primers::DesignOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TailingQuery {
    pub document: DocumentRef,
    pub recipe: verification::terminal::TailSubstrate,
    pub base: collagene_core::cloning::ligation::tailing::TailingBase,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapProfileQuery {
    pub documents: Vec<DocumentRef>,
    pub reverse: Vec<bool>,
    pub overlaps: Vec<usize>,
    pub circular: bool,
    pub options: collagene_core::cloning::overlap::profiles::ProfileOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EndRepairQuery {
    pub document: DocumentRef,
    pub enzymes: Vec<String>,
    pub fragment_index: usize,
    pub method: chemistry::EndRepairMethod,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnealQuery {
    pub top: DocumentRef,
    pub bottom: DocumentRef,
    pub options: chemistry::AnnealOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGateQuery {
    pub documents: Vec<DocumentRef>,
    pub options: collagene_core::cloning::golden_gate::GoldenGateOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MismatchPcrQuery {
    pub document: DocumentRef,
    pub forward: Primer,
    pub reverse: Primer,
    pub options: collagene_core::cloning::pcr::mutagenesis::PcrOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InversePcrQuery {
    pub document: DocumentRef,
    pub target: DocumentRef,
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin: Option<collagene_core::material::run::MaterialSelection>,
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin_workspace: Option<WorkspaceRef>,
    pub options: collagene_core::cloning::pcr::mutagenesis::InversePcrOptions,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoundGelLane {
    pub label: String,
    pub document: DocumentRef,
    pub preparation: collagene_core::output::gel::Preparation,
    pub mass_ng: Option<f64>,
    /// Required for CatalogDigest: native retained owner/event/revision proof.
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GelQuery {
    pub lanes: Vec<BoundGelLane>,
    pub calibration: collagene_core::output::gel::Calibration,
    pub image_format: collagene_core::output::raster::Format,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DiagnosticQuery {
    pub documents: Vec<DocumentRef>,
    pub options: collagene_core::output::gel::DiagnosticOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Restriction {
    pub document: DocumentRef,
    /// One name from the embedded starter catalog, e.g. EcoRI or BsaI.
    pub enzyme: String,
}
/// Inspect actual PCR-product digest geometry while keeping the original source immutable.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PcrDigestInspectionQuery {
    pub document: DocumentRef,
    pub request: collagene_core::cloning::pcr::digest::Request,
    /// Explicit mismatch-PCR policy; omitted keeps the exact-PCR inspector unchanged.
    #[serde(default)]
    pub pcr_options: Option<collagene_core::cloning::pcr::mutagenesis::PcrOptions>,
    /// Explicit documentary origin for a saved-run source; never inferred from metadata.
    #[serde(default)]
    pub origin: Option<collagene_core::material::run::MaterialOrigin>,
    /// Required exact retained native workspace when a saved-run origin is supplied.
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
    #[serde(default)]
    pub fragment_start: usize,
    #[serde(default = "default_inspection_page")]
    pub max_fragments: usize,
}
/// Inspect the actual final outer-PCR product of a complete physical OE preparation.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapExtensionDigestInspectionQuery {
    pub document: DocumentRef,
    pub request: collagene_core::cloning::overlap::extension::digest::Request,
    #[serde(default)]
    pub origin: Option<collagene_core::material::run::MaterialOrigin>,
    /// Required exact retained native workspace when an origin is supplied.
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
    #[serde(default)]
    pub fragment_start: usize,
    #[serde(default = "default_inspection_page")]
    pub max_fragments: usize,
}
/// Complete catalog-bound restriction evidence from one immutable source.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "source", deny_unknown_fields)]
pub enum RestrictionInspectionEvidence {
    Original {
        inspection: collagene_core::cloning::restriction::catalog::CatalogRestrictionInspection,
    },
    SavedRun {
        inspection: collagene_core::cloning::restriction::run_digest::Inspection,
    },
}
/// Bounded recognition evidence and assumed-complete digest summaries for one immutable source.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RestrictionInspectionQuery {
    pub document: DocumentRef,
    /// Stable source role; included in every report binding independently of document identity.
    pub source_role: String,
    /// Explicit saved material source; requires its exact native workspace.
    #[serde(default)]
    pub origin: Option<collagene_core::material::capsule::Capsule>,
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
    /// One to 64 exact catalog enzyme IDs or uniquely resolved names/aliases. Formulations are separate.
    pub enzymes: Vec<String>,
    /// Optional exact current catalog edition/hash/ordered-ID/geometry pin; stale or changed pins fail.
    pub catalog_pin: Option<collagene_core::cloning::restriction::catalog::CatalogSelection>,
    #[serde(default)]
    pub options: collagene_core::cloning::restriction::inspection::InspectionOptions,
    #[serde(default)]
    pub site_start: usize,
    #[serde(default = "default_inspection_page")]
    pub max_sites: usize,
    #[serde(default)]
    pub fragment_start: usize,
    #[serde(default = "default_inspection_page")]
    pub max_fragments: usize,
}
/// Recognition/cut-coordinate evidence only; never constructs a physical digest.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RestrictionScanQuery {
    pub document: DocumentRef,
    /// Stable explicit source role, independently bound to source revision.
    pub source_role: String,
    /// One to 64 exact catalog IDs or uniquely resolved names/aliases, in order.
    pub enzymes: Vec<String>,
    /// Optional exact edition/hash/ordered-ID/geometry pin; changed pins fail.
    pub catalog_pin: Option<collagene_core::cloning::restriction::catalog::CatalogSelection>,
    #[serde(default)]
    pub options: collagene_core::cloning::restriction::inspection::InspectionOptions,
    #[serde(default)]
    pub site_start: usize,
    #[serde(default = "default_inspection_page")]
    pub max_sites: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DigestQuery {
    pub document: DocumentRef,
    pub enzymes: Vec<String>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Pcr {
    pub document: DocumentRef,
    pub forward: Primer,
    pub reverse: Primer,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapPrimerDesignQuery {
    /// Actual original templates with exact immutable source revisions.
    pub inputs: Vec<BoundInput>,
    /// Independently supplied intended DNA; never manufactured from chosen intervals.
    pub target: DocumentRef,
    /// Explicit ordered template intervals, orientation and bounded heuristic controls.
    pub request: collagene_core::cloning::overlap::primer_design::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapExtensionPrimerDesignQuery {
    /// Immutable original or exact saved native source snapshot; saved sources require their owning workspace.
    pub document: DocumentRef,
    /// Construct-oriented contribution in the explicit stored-strand frame, independently declared edits and bounded primer policy.
    pub request: collagene_core::cloning::overlap::extension::primers::Request,
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedOverlapPrimerDesignQuery {
    pub inputs: Vec<BoundInput>,
    pub request: EditedOverlapPrimerRequest,
    /// One immutable native archive owning every selected saved-source role; null for originals.
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGateJunctionPrimersQuery {
    pub inputs: Vec<BoundInput>,
    pub request: collagene_core::cloning::golden_gate::junction_edit::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MutationIntentQuery {
    pub document: DocumentRef,
    pub request: collagene_core::cloning::pcr::mutation_intent::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGatePrimerDesignQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub request: collagene_core::cloning::golden_gate::primers::Request,
}
#[tool_router(router = cloning_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Reconstruct explicit single-nucleotide tailing of BOTH DNA strand3-prime ends from a complete digest (optionally blunted) or a unique exact-PCR product. Requires actual blunt linear duplex and accessibleOH; retains five-prime chemistry and nicks. This idealendpoint transform does not predict enzyme/kit activity and cannot establish physical authority when imported as metadata. Full immutable source/recipe/strand report artifact included; use explicit Tail3Prime in a targetbound plan."
    )]
    pub(crate) async fn simulate_tailing(
        &self,
        Parameters(args): Parameters<TailingQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let preview=verification::terminal::simulate_tailing(&source,&args.recipe,args.base).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"tailing":preview}))?;
            output(json!({"document":args.document,"source_revision":preview.source_revision,"recipe":preview.recipe,"base":preview.base,"engine_version":verification::ENGINE_VERSION,"top_length":preview.tailing.duplex.top.len(),"bottom_length":preview.tailing.duplex.bottom.len(),"bottom_start":preview.tailing.duplex.bottom_start,"additions":preview.tailing.additions,"warnings":preview.tailing.warnings,"artifact_id":artifact,"authority":"Explicit ideal tailing evidence only. Imported preview metadata cannot establish physical end state, template ancestry, enzyme activity or target verification."}))
        }).await
    }
    #[tool(
        description = "Check explicitly oriented ordered DNA overlap junctions against the declared Gibson/HiFi/In-Fusion-style profile and bounded repeat evidence. Source orientation, exact homology and length/GC findings are retained. Vendor recommendations are advisory: this neither simulates proprietary reagents nor establishes experimental compatibility, uniqueness or physical substrate chemistry. Full source/options-bound JSON artifact included; reconstruct ProfiledOverlap against an independent target to verify the candidate."
    )]
    pub(crate) async fn check_overlap_profile(
        &self,
        Parameters(args): Parameters<OverlapProfileQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let sources={let store=server.session.store.lock().unwrap();args.documents.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?};
            let report=collagene_core::cloning::overlap::profiles::check_oriented_fragments(&sources,&args.reverse,&args.overlaps,args.circular,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"documents":args.documents,"engine_version":verification::ENGINE_VERSION,"overlap_profile":report}))?;
            let mut preview=report.clone();preview.junctions.truncate(50);for junction in &mut preview.junctions{junction.occurrences.truncate(20);if junction.homology.len()>600{junction.homology.truncate(600);}}
            output(json!({"documents":args.documents,"reverse":args.reverse,"engine_version":verification::ENGINE_VERSION,"geometry_valid":report.geometry_valid,"recommendation_compliant":report.recommendation_compliant,"complete":report.complete,"exhaustive_alternative_search":false,"comparisons":report.comparisons,"junction_count":report.junctions.len(),"display_truncated":report.junctions.len()>50||report.junctions.iter().any(|j|j.occurrences.len()>20||j.homology.len()>600),"junctions":preview.junctions,"findings":report.findings.iter().take(100).collect::<Vec<_>>(),"warnings":report.warnings.iter().take(100).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Reconstruct a complete restriction digest and selected fragment, then apply explicit 3-prime fill-in/exonuclease end repair. Preserves 5-prime chemistry and internal nicks. Imported duplex metadata never grants physical authority. Returns bounded sequence/end summaries and full source-bound report artifact; use Digest then RepairEnds in a separately target-bound plan for verification."
    )]
    pub(crate) async fn simulate_end_repair(
        &self,
        Parameters(args): Parameters<EndRepairQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let preview=verification::terminal::simulate_end_repair(&source,&args.enzymes,args.fragment_index,args.method).map_err(bad)?;
            let summary=preview.repair.duplex.summary().map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"repair":preview}))?;
            output(json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"length":preview.molecule.sequence.len(),"sequence_preview":preview.molecule.sequence.as_str().chars().take(200).collect::<String>(),"sequence_preview_truncated":preview.molecule.sequence.len()>200,"summary":summary,"top_added_3prime":preview.repair.top_added_3prime.chars().take(200).collect::<String>(),"top_added_count":preview.repair.top_added_3prime.len(),"top_removed_count":preview.repair.top_removed_3prime.len(),"bottom_added_count":preview.repair.bottom_added_3prime.len(),"bottom_removed_count":preview.repair.bottom_removed_3prime.len(),"annotation_count":preview.molecule.features.len(),"findings":preview.findings.iter().take(100).collect::<Vec<_>>(),"artifact_id":artifact,"authority":"Read-only preview. Physical intermediates must be reconstructed in the bound plan; this artifact is not an imported end-chemistry authority."}))
        }).await
    }
    #[tool(
        description = "Enumerate bounded exact antiparallel placements of two explicitly declared linear 5-prime-to-3-prime ssDNA oligos with supplied terminal chemistry. Does not assume automatic longest/unique pairing, phosphorylation, ligation or experimental success. Returns alternatives, completeness and a full immutable-source-bound artifact. Use an explicit AnnealOligos candidate in a separately bound plan and reverify."
    )]
    pub(crate) async fn simulate_oligo_annealing(
        &self,
        Parameters(args): Parameters<AnnealQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(top,bottom)={let store=server.session.store.lock().unwrap();(store.document(&args.top)?.clone(),store.document(&args.bottom)?.clone())};
            let preview=verification::terminal::simulate_oligo_annealing(&top,&bottom,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"top":args.top,"bottom":args.bottom,"engine_version":verification::ENGINE_VERSION,"annealing":preview}))?;
            output(json!({"top":args.top,"bottom":args.bottom,"engine_version":verification::ENGINE_VERSION,"candidate_count":preview.annealing.candidates.len(),"complete":preview.annealing.complete,"unique":preview.annealing.unique,"comparisons":preview.annealing.comparisons,"options":preview.annealing.options,"warnings":preview.annealing.warnings,"display_truncated":preview.annealing.candidates.len()>50,"candidates":preview.annealing.candidates.iter().enumerate().take(50).map(|(index,c)|json!({"candidate_index":index,"id":c.id,"bottom_start":c.duplex.bottom_start,"top_length":c.duplex.top.len(),"bottom_length":c.duplex.bottom.len(),"paired_start":c.paired_start,"paired_end":c.paired_end,"paired_bases":c.paired_bases})).collect::<Vec<_>>(),"artifact_id":artifact,"authority":"Terminal chemistry is an explicit caller assumption. Base pairing creates no covalent bonds and imported metadata cannot establish a duplex state."}))
        }).await
    }
    #[tool(
        description = "Enumerate finite one-copy Golden Gate fragment assemblies using actual tracked-position cleavage checks. These bare originals have unknown modification state; unresolved cleavage withholds the requested pool. Use inspect_golden_gate_reaction for explicit original-generation or saved-material/PCR declarations. Distinguishes stable site-free candidates, recut paths, uncut carryover, source histories and budget completeness. This is not an exhaustive kinetic reaction-mixture search. Full JSON evidence artifact included; candidate index can be used in a separately target-bound GoldenGate plan step."
    )]
    pub(crate) async fn simulate_golden_gate(
        &self,
        Parameters(args): Parameters<GoldenGateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let sources={let store=server.session.store.lock().unwrap();args.documents.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?};
            let report=collagene_core::cloning::golden_gate::simulate(&sources,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"documents":args.documents,"engine_version":verification::ENGINE_VERSION,"golden_gate":report}))?;
            output(json!({"documents":args.documents,"candidate_count":report.candidates.len(),"profile_revision":report.profile_revision,"products_determined":report.products_determined,"cleavage_comparisons":report.comparisons,"stable_physical_product_count":report.stable_physical_product_count,"enumeration_complete":report.enumeration_complete,"exhaustive_mixture_search":report.exhaustive_mixture_search,"states_examined":report.states_examined,"carryover_sources":report.carryover_sources,"truncation_reasons":report.truncation_reasons,"warnings":report.warnings,"display_truncated":report.candidates.len()>50,"candidates":report.candidates.iter().enumerate().take(50).map(|(i,c)|json!({"candidate_index":i,"length":c.duplex.top.len(),"topology":c.duplex.topology,"path":c.path,"restriction_site_free":c.restriction_site_free,"recut_site_count":c.recut_sites.len(),"recut_complete":c.cleavage.complete,"recut_products_determined":c.cleavage.products_determined,"directed_sequence_key":c.directed_sequence_key,"physical_sequence_key":c.physical_sequence_key})).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Reconstruct intentional-mismatch PCR from immutable template and actual tailed oligos under an explicit mismatch/3-prime/work policy. Returns introduced changes, competing primer roles, overlap/ambiguity diagnostics and search completeness. Overlapping-domain/uncertain products remain previews. Full report artifact is computed, not an AI assertion."
    )]
    pub(crate) async fn simulate_mismatch_pcr(
        &self,
        Parameters(args): Parameters<MismatchPcrQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=collagene_core::cloning::pcr::mutagenesis::simulate_pcr(&source,&args.forward,&args.reverse,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"forward":args.forward,"reverse":args.reverse,"pcr":report}))?;
            output(json!({"document":args.document,"product_count":report.products.len(),"search_complete":report.search_complete,"enumeration_complete":report.enumeration_complete,"unsupported_geometries":report.unsupported_geometries,"conflicting_overlaps":report.conflicting_overlaps,"warnings":report.warnings,"display_truncated":report.products.len()>100,"products":report.products.iter().enumerate().take(100).map(|(i,p)|json!({"product_index":i,"length":p.sequence.len(),"forward_start":p.forward_start,"reverse_start":p.reverse_start,"forward_primer_index":p.forward_primer_index,"reverse_primer_index":p.reverse_primer_index,"introduced_changes":p.substitutions.iter().take(20).collect::<Vec<_>>(),"change_count":p.substitutions.len(),"certain":p.certain,"overlapping_domains":p.overlapping_domains})).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Design inverse-PCR replacement/insertion/deletion primers from an explicit circular source and independently supplied circular target. origin and origin_workspace are required nullable fields: null/null declares an original reference hypothesis; a saved Run/Derived selection requires its exact immutable native workspace and fresh source-event attestation. Reconstructs the origin when present, then actual PCR, phosphorylation and circular ligation. Nicked/unsupported templates refuse without repair or metadata authority. Atomically stores an unverified candidate and full source-bound artifact; inherited incomplete findings remain explicit. Verify the returned plan_revision before export/use; no persistent project writes or experimental success are implied."
    )]
    pub(crate) async fn design_inverse_pcr(
        &self,
        Parameters(args): Parameters<InversePcrQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            if args.origin.is_some()!=args.origin_workspace.is_some(){return Err(bad("Saved inverse-PCR origin and exact native workspace must be provided together."));}
            let frame=args.origin_workspace.as_ref().map(|w|server.session.restore_workspace(w)).transpose()?;
            let origin=match (&args.origin,&frame){
                (Some(selection),Some(frame))=>Some(server.session.resolve_material_origin(frame,selection.event(),selection)?),
                (None,None)=>None,
                _=>return Err(bad("Incomplete inverse-PCR origin binding.")),
            };
            let (source,target,expected)={
                let store=server.session.store.lock().unwrap();let source=store.document(&args.document)?;let target=store.document(&args.target)?;
                if let (Some(frame),Some(origin))=(&frame,&origin){if frame.revisions.get(&origin.event().revision)!=Some(source){return Err(bad("Inverse-PCR source differs from the exact retained material event."));}}
                let expected=collagene_core::cloning::pcr::mutagenesis::inverse_pcr_revision(source,target,origin.as_ref(),&args.options).map_err(bad)?;
                (source.clone(),target.clone(),expected)
            };
            let design=collagene_core::cloning::pcr::mutagenesis::design_inverse_pcr(&source,&target,origin.as_ref(),&args.options).map_err(bad)?;
            design.validate_binding(&expected).map_err(bad)?;
            transfer::Candidate::new(design.plan.clone()).validate().map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{document:&'a DocumentRef,target:&'a DocumentRef,origin:&'a Option<collagene_core::material::run::MaterialSelection>,origin_workspace:&'a Option<WorkspaceRef>,source_event:Option<&'a collagene_core::material::lineage::EventBinding>,engine_version:&'static str,inverse_pcr:&'a collagene_core::cloning::pcr::mutagenesis::InversePcrDesign}
            let text=bounded_json_limit(&Evidence{document:&args.document,target:&args.target,origin:&args.origin,origin_workspace:&args.origin_workspace,source_event:origin.as_ref().map(|o|o.event()),engine_version:verification::ENGINE_VERSION,inverse_pcr:&design},32*1024*1024)?;
            let artifact_id=blake3::hash(text.as_bytes()).to_hex().to_string();let plan_revision=collagene_core::planning::plan::DesignPlan::revision(&design.plan);
            let plan_bytes=bounded_json_limit(&design.plan,32*1024*1024)?.len().checked_add(bounded_json_limit(&target,32*1024*1024)?.len()).ok_or_else(||bad("Inverse-PCR storage accounting overflow."))?;
            let result=output(json!({"document":args.document,"target":args.target,"origin":args.origin,"origin_workspace":args.origin_workspace,"source_event":origin.as_ref().map(|o|o.event()),"context_revision":expected,"plan_revision":plan_revision,"verified":false,"forward":design.forward,"reverse":design.reverse,"source_revision":design.source_revision,"target_revision":design.target_revision,"computed_status":design.verification.status,"findings":design.verification.findings,"warnings":design.warnings,"artifact_id":artifact_id}))?;
            let mut store=server.session.store.lock().unwrap();let new_plan=!store.plans.contains_key(&plan_revision);let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id){if !matches!(existing.data,ArtifactData::Text(_))||existing.filename!="sequence-analysis.json"||existing.mime_type!="application/json"{return Err(bad("Inverse-PCR artifact has another representation; no candidate was published."));}}
            if (new_plan&&store.plans.len()>=MAX_ITEMS)||(new_artifact&&store.artifacts.len()>=MAX_ITEMS){return Err(bad("Inverse-PCR candidate/artifact capacity reached; nothing was published."));}
            let total=(if new_plan{plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("Inverse-PCR storage accounting overflow."))?;store.reserve(total)?;
            if new_plan{store.plans.insert(plan_revision,StoredPlan{plan:design.plan,target});}
            if new_artifact{store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"sequence-analysis.json"});}
            Ok(result)
        }).await
    }
    #[tool(
        description = "Predict a reference gel from revision-bound uncut, complete-digest or PCR lanes under an explicit ladder/calibration. CatalogDigest lanes require a native workspace reference proving the recipe's exact retained document/event/revision and a matching session document; historical sources need not be current heads. Native inspection is recomputed. Circular/unsupported profiles remain unpositioned. Select image_format Svg or Png with explicit pixel width. JSON retains declared lane recipes and references; rendered images remain predictions. No observed-gel/OpenGel integration."
    )]
    pub(crate) async fn predict_gel(
        &self,
        Parameters(args): Parameters<GelQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            args.image_format.validate().map_err(bad)?;
            if args.lanes.is_empty() || args.lanes.len() > 32 {
                return Err(bad("gel needs 1–32 source lanes"));
            }
            args.calibration.validate().map_err(bad)?;
            for lane in &args.lanes {
                if lane.label.len() > 256 || lane.label.chars().any(char::is_control)
                    || lane.mass_ng.is_some_and(|mass| !mass.is_finite() || mass <= 0. || mass > 1e9)
                {
                    return Err(bad("invalid gel lane label or loaded mass"));
                }
                if let collagene_core::output::gel::Preparation::CatalogDigest { recipe } = &lane.preparation {
                    recipe.validate().map_err(bad)?;
                    let reference = lane.workspace.as_ref().ok_or_else(|| bad(
                        "CatalogDigest requires a native workspace reference proving its exact retained source event"
                    ))?;
                    if reference.id.len() != 64 || reference.revision.len() != 64 {
                        return Err(bad("invalid pinned gel workspace reference"));
                    }
                }
            }
            // Resolve immutable session documents together, then release this
            // lock before restore_workspace acquires the same store mutex.
            let lanes = {
                let store = server.session.store.lock().unwrap();
                args.lanes.iter().map(|lane| Ok(collagene_core::output::gel::LaneInput {
                    label: lane.label.clone(),
                    molecule: store.document(&lane.document)?.clone(),
                    preparation: lane.preparation.clone(),
                    mass_ng: lane.mass_ng,
                })).collect::<Result<Vec<_>, ErrorData>>()?
            };
            let mut groups = std::collections::BTreeMap::<String, (WorkspaceRef, Vec<usize>)>::new();
            for (index, lane) in args.lanes.iter().enumerate() {
                if !matches!(&lane.preparation, collagene_core::output::gel::Preparation::CatalogDigest { .. }) { continue; }
                let reference = lane.workspace.as_ref().ok_or_else(|| bad("pinned gel workspace reference is unavailable"))?;
                if reference.id != reference.revision { return Err(bad("stale workspace revision")); }
                groups.entry(reference.id.clone()).or_insert_with(|| (reference.clone(), Vec::new())).1.push(index);
            }
            // Restore each distinct workspace once, then drop it before loading
            // another: shared large histories are never cloned once per lane.
            for (_, (reference, indices)) in groups {
                let frame = server.session.restore_workspace(&reference)?;
                for index in indices {
                    let lane = &lanes[index];
                    let collagene_core::output::gel::Preparation::CatalogDigest { recipe } = &lane.preparation else { continue; };
                    let source = &recipe.source;
                    if !frame.documents.iter().any(|document| document.id == source.document_id) {
                        return Err(bad("pinned gel source owner is not retained in this workspace"));
                    }
                    let event = frame.history.iter().find(|event| event.id == source.event_id)
                        .ok_or_else(|| bad("the exact pinned gel source event is not retained"))?;
                    if event.document_id != source.document_id || event.revision != source.revision {
                        return Err(bad("pinned gel source owner/event/revision mismatch"));
                    }
                    let original = frame.revisions.get(&source.revision)
                        .ok_or_else(|| bad("pinned gel historical molecule is not retained"))?;
                    if original != &lane.molecule || lane.molecule.revision() != source.revision {
                        return Err(bad("session lane document differs from the exact historical pinned gel source"));
                    }
                    recipe.validate_for_workspace(&frame).map_err(bad)?;
                }
            }
            let report=collagene_core::output::gel::predict_gel(&lanes,&args.calibration).map_err(bad)?;
            let image = collagene_core::output::gel_export::Request::new(report.clone(), args.image_format).map_err(bad)?.render().map_err(bad)?;
            let image_artifact = blake3::hash(image.image.bytes()).to_hex().to_string();
            let (image_data, image_mime, image_encoding, image_filename) = match image.image {
                collagene_core::output::raster::Image::Svg(svg) => (ArtifactData::Text(svg), "image/svg+xml", "utf-8", "predicted-gel.svg"),
                collagene_core::output::raster::Image::Png(png) => (ArtifactData::Binary(png.bytes), "image/png", "base64", "predicted-gel.png"),
            };
            let evidence = bounded_json_limit(&json!({"documents":args.lanes.iter().map(|l|&l.document).collect::<Vec<_>>(),"lanes":args.lanes,"engine_version":verification::ENGINE_VERSION,"gel":report,"image_format":args.image_format,"image_context_revision":image.context_revision,"image_warnings":image.warnings}),32*1024*1024)?;
            let artifact = blake3::hash(evidence.as_bytes()).to_hex().to_string();
            if artifact == image_artifact { return Err(bad("gel artifact identities collide")); }
            let result = output(json!({"model":report.model,"lane_count":report.lanes.len(),"lanes":report.lanes.iter().map(|l|json!({"label":l.label,"source_revision":l.source_revision,"fragment_count":l.fragments.len(),"resolved_band_count":l.bands.len(),"unpositioned_fragment_count":l.unpositioned_fragments.len(),"complete_under_model":l.complete_under_model,"error":l.error,"warnings":l.warnings})).collect::<Vec<_>>(),"artifact_id":artifact,"image_artifact_id":image_artifact,"image_format":args.image_format,"image_context_revision":image.context_revision,"image_warnings":image.warnings}))?;
            // Reserve and publish the complete image/report pair atomically. A large PNG
            // must not leave an orphan artifact when the report cannot fit the session.
            let mut store = server.session.store.lock().unwrap();
            if let Some(existing) = store.artifacts.get(&image_artifact) {
                if !matches!((&existing.data,&image_data),(ArtifactData::Text(_),ArtifactData::Text(_))|(ArtifactData::Binary(_),ArtifactData::Binary(_)))
                    || existing.mime_type != image_mime || existing.encoding != image_encoding || existing.filename != image_filename {
                    return Err(bad("gel image artifact has another representation; nothing published"));
                }
            }
            if let Some(existing) = store.artifacts.get(&artifact) {
                if !matches!(existing.data, ArtifactData::Text(_)) || existing.mime_type != "application/json" || existing.encoding != "utf-8" || existing.filename != "sequence-analysis.json" {
                    return Err(bad("gel evidence artifact has another representation; nothing published"));
                }
            }
            let new_image = !store.artifacts.contains_key(&image_artifact);
            let new_report = !store.artifacts.contains_key(&artifact);
            if store.artifacts.len() + usize::from(new_image) + usize::from(new_report) > MAX_ITEMS { return Err(bad("gel artifact capacity reached; nothing published")); }
            let image_len = match &image_data { ArtifactData::Text(text) => text.len(), ArtifactData::Binary(bytes) => bytes.len() };
            let total = (if new_image {image_len}else{0}).checked_add(if new_report {evidence.len()}else{0}).ok_or_else(||bad("gel artifact size overflow"))?;
            store.reserve(total)?;
            if new_image { store.artifacts.insert(image_artifact, Artifact { data:image_data, mime_type:image_mime, encoding:image_encoding, filename:image_filename }); }
            if new_report { store.artifacts.insert(artifact, Artifact { data:ArtifactData::Text(evidence), mime_type:"application/json", encoding:"utf-8", filename:"sequence-analysis.json" }); }
            Ok(result)
        }).await
    }

    #[tool(
        description = "Rank bounded single/double restriction digests that distinguish immutable candidate plasmids under an explicit gel calibration/resolution. Reports distinguished, co-migrating and unknown pairs separately, with evaluated/work bounds and exhaustive-search flag. Unsupported profiles do not become negative specificity evidence. Full JSON evidence artifact included."
    )]
    pub(crate) async fn design_diagnostic_digest(
        &self,
        Parameters(args): Parameters<DiagnosticQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let sources={let store=server.session.store.lock().unwrap();args.documents.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?};
            let report=collagene_core::output::gel::diagnostic_digests(&sources,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"documents":args.documents,"engine_version":verification::ENGINE_VERSION,"diagnostic_digests":report}))?;
            output(json!({"documents":args.documents,"source_revisions":report.source_revisions,"evaluated_digests":report.evaluated_digests,"work_bases":report.work_bases,"exhaustive_search":report.exhaustive_search,"warnings":report.warnings,"display_truncated":report.candidates.len()>20,"candidates":report.candidates.iter().take(20).map(|c|json!({"enzymes":c.enzymes,"distinguished_pairs":c.distinguished_pairs,"indistinguishable_pairs":c.indistinguishable_pairs,"unknown_pairs":c.unknown_pairs,"minimum_diagnostic_separation":c.minimum_diagnostic_separation})).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Resolve a saved primer_bind annotation at its exact immutable document revision and zero-based feature index. Uses shared strict oligo/annealing/5-prime-tail consistency checks; missing fields may be inferred only from exact supported strand-oriented geometry, never fuzzy/unresolved/reviewed geometry. Returns the actual ordered Primer, sequence origin and warnings, preserving deliberate mismatches/IUPAC and historical source qualifiers as evidence. No binding, PCR, thermodynamic or experimental proof; pass the returned primer explicitly to simulate_pcr/analyze_primer."
    )]
    pub(crate) async fn resolve_primer(
        &self,
        Parameters(args): Parameters<SavedPrimerQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let resolved={let store=server.session.store.lock().unwrap();let source=store.document(&args.document)?;
                collagene_core::cloning::pcr::reuse::resolve_at_revision(source,args.feature_index,&args.document.revision).map_err(bad)?};
            output(json!({"document":args.document,"resolved":resolved,
                "authority":"Saved oligo identity only. Current template binding and exact PCR products must be recomputed independently; source location and historical analysis are not amplification evidence."}))
        }).await
    }
    #[tool(
        description = "Analyze a tailed DNA primer under explicit concentration/sodium conditions, screen contiguous self/partner complementarity, and find bounded binding sites on an immutable template. Certain 3-prime matches, mismatches and ambiguous compatibility remain separate. This is not Primer3, empirical structure prediction or genome-wide specificity. Full JSON evidence artifact included."
    )]
    pub(crate) async fn analyze_primer(
        &self,
        Parameters(args): Parameters<PrimerQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let analysis=primers::analyze_primer(&args.primer,&args.conditions).map_err(bad)?;
            let bindings=primers::find_binding_sites(&source,&args.primer,&args.binding_options).map_err(bad)?;
            let heterodimer=args.partner.as_ref().map(|p|primers::heterodimer(&args.primer,p)).transpose().map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"analysis":analysis,"bindings":bindings,"partner":args.partner,"heterodimer":heterodimer}))?;
            output(json!({"document":args.document,"analysis":analysis,"heterodimer":heterodimer,"bindings":bindings.hits.iter().take(100).collect::<Vec<_>>(),"binding_count":bindings.hits.len(),"search_truncated":bindings.truncated,"display_truncated":bindings.hits.len()>100,"binding_options":bindings.options,"warnings":bindings.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Design bounded PCR primer pairs anchored to an explicit reference interval (circular end may be unwrapped). Uses recorded nearest-neighbor conditions, GC/length/homopolymer/complementarity filters and separated synthetic tails. Candidates require binding-alternative/PCR and laboratory review; no automatic verification status or document edits. Full JSON artifact retains source/options and ranked pairs."
    )]
    pub(crate) async fn design_primers(
        &self,
        Parameters(args): Parameters<PrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=primers::design_primers(&source,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"design":report}))?;
            output(json!({"document":args.document,"options":report.options,"pairs":report.pairs.iter().take(20).collect::<Vec<_>>(),"count":report.pairs.len(),"display_truncated":report.pairs.len()>20,"evaluated_pairs":report.evaluated_pairs,"rejected_pairs":report.rejected_pairs,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Scan an immutable DNA source for restriction recognition spans and reference top/bottom cut boundaries without constructing duplexes, fragments or products. Select ordered exact catalog IDs or uniquely resolved names/aliases, optionally pinning edition/hash/full geometry. Certain and possible recognition are separate; unavailable cuts retain raw boundary diagnostics. Returns bounded site pages, validated source/geometry/options evidence and full JSON artifact. Per-enzyme observed counts are lower bounds when scan_complete=false; exact counts and single-cutter classification are null on incomplete scans, and possible recognition withholds single-cutter certainty. Scan completion describes sequence search only, not enzyme activity, methylation, buffer compatibility or empirical cleavage. Use inspect_restriction separately only when assumed-complete physical digest summaries are wanted."
    )]
    pub(crate) async fn scan_restriction(
        &self,
        Parameters(args): Parameters<RestrictionScanQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_sites) {
            return Err(bad("restriction site page size must be 1..100"));
        }
        if args.enzymes.is_empty() || args.enzymes.len() > 64 {
            return Err(bad("select 1..64 embedded enzymes"));
        }
        let catalog = catalog_tools::builtin()?;
        let ids = args
            .enzymes
            .iter()
            .map(|name| match catalog.lookup_enzyme(name) {
                collagene_enzyme_catalog::Lookup::Unique(enzyme) => Ok(enzyme.id.clone()),
                collagene_enzyme_catalog::Lookup::Missing => {
                    Err(bad("unknown enzyme identity; call lookup_enzyme_catalog"))
                }
                collagene_enzyme_catalog::Lookup::Ambiguous(_) => {
                    Err(bad("ambiguous enzyme alias; select exact catalog IDs"))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let selection =
            collagene_core::cloning::restriction::catalog::select(&catalog, &ids).map_err(bad)?;
        if let Some(pin) = &args.catalog_pin {
            pin.validate_for(&catalog).map_err(bad)?;
            if pin != &selection {
                return Err(bad("catalog pin and ordered enzyme selection differ"));
            }
        }
        let enzymes = selection.resolved_enzymes(&catalog).map_err(bad)?;
        let server = self.clone();
        self.compute(move|| {
            let molecule={
                let store=server.session.store.lock().unwrap();
                let molecule=store.document(&args.document)?;
                collagene_core::cloning::restriction::inspection::source_preflight(&args.source_role,molecule,&enzymes,&args.options).map_err(bad)?;
                molecule.clone()
            };
            let report=collagene_core::cloning::restriction::inspection::scan(&args.source_role,&molecule,&enzymes,&args.options).map_err(bad)?;
            report.validate_bindings_for(&args.source_role,&molecule,&enzymes,&args.options).map_err(bad)?;
            if args.site_start>report.sites.len() {return Err(bad("restriction page offset exceeds retained evidence count"));}
            let mut counts=vec![(0usize,0usize);enzymes.len()];
            for site in &report.sites {
                let entry=&mut counts[site.enzyme_index];
                if site.certain {entry.0+=1;}else{entry.1+=1;}
            }
            let counts:Vec<_>=counts.into_iter().enumerate().map(|(i,(certain,possible))|json!({
                "enzyme_index":i,"catalog_id":selection.enzymes[i].id,"enzyme":enzymes[i].name,
                "counts_complete":report.scan_complete,"certain_observed":certain,"possible_observed":possible,
                "observed_site_count":certain+possible,
                "exact_site_count":(report.scan_complete && possible==0).then_some(certain),
                "single_cutter":if report.scan_complete && possible==0 {Some(certain==1)}else{None},
            })).collect();
            let authority="Recognition/cut-coordinate search evidence only; no duplex or fragment construction, no experimental cleavage or reaction-condition authority. Observed counts are lower bounds when scan is incomplete; possible IUPAC matches do not establish exact biological site counts.";
            let artifact_id=server.session.text_artifact(bounded_json(&json!({"catalog_pin":selection,"report":report,"enzyme_counts":counts,"authority":authority}))?,"application/json","restriction-scan.json")?;
            let end=args.site_start.saturating_add(args.max_sites).min(report.sites.len());
            output(json!({
                "document":args.document,"source_role":report.source_role,"source_revision":report.source_revision,
                "source_length":report.source_length,"topology":report.topology,"engine_version":report.engine_version,
                "batch_revision":report.batch_revision,"geometry_revision":report.geometry_revision,
                "options_revision":report.options_revision,"options":report.options,"enzymes":report.enzymes,
                "catalog_pin":selection,"scan_complete":report.scan_complete,"products_determined":false,
                "fragment_count":0,"ref_top_lengths":null,"comparisons":report.comparisons,
                "limits":report.limits,"diagnostics":report.diagnostics,"site_count":report.sites.len(),
                "enzyme_counts":counts,"site_start":args.site_start,"next_site_start":end,
                "sites_done":end==report.sites.len(),"sites":&report.sites[args.site_start..end],
                "artifact_id":artifact_id,"authority":authority,
            }))
        }).await
    }
    #[tool(
        description = "Inspect catalog-pinned cleavage of the actual outer-PCR product from an explicit four-step overlap-extension preparation. Freshly reconstructs both original-source PCR branches, selected physical strand extension and outer PCR before scanning/cutting its actual product; no intended target or claimed product is an input. Retains each stage's geometry counts, terminal-pair uncertainty, source/origin pins and full immutable evidence. Optional saved origin requires its exact native workspace and fresh replay; imported physical summaries never grant substrate authority. Returns a bounded fragment page and explicit selectable/incomplete status; it creates no stock, candidate, target or chosen cloning fragment."
    )]
    pub(crate) async fn inspect_overlap_extension_digest(
        &self,
        Parameters(args): Parameters<OverlapExtensionDigestInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_fragments) {
            return Err(bad(
                "overlap-extension digest fragment page size must be1..100",
            ));
        }
        if args.origin.is_some() && args.workspace.is_none() {
            return Err(bad(
                "saved overlap-extension source requires its exact native workspace reference",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let workspace = args.workspace.as_ref().map(|reference|server.session.restore_workspace(reference)).transpose()?;
            let source = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_core::cloning::overlap::extension::digest::source_preflight(source, args.origin.as_ref(), &args.request, workspace.as_ref()).map_err(bad)?;
                source.clone()
            };
            let expected = collagene_core::cloning::overlap::extension::digest::inspection_revision(&source, args.origin.as_ref(), &args.request, workspace.as_ref()).map_err(bad)?;
            let preview = collagene_core::cloning::overlap::extension::digest::inspect(&source, args.origin.as_ref(), &args.request, workspace.as_ref()).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            let count = preview.inspection.as_ref().map_or(0, |inspection|inspection.fragments.len());
            if args.fragment_start > count { return Err(bad("overlap-extension digest fragment offset exceeds retained evidence")); }
            let end = args.fragment_start.saturating_add(args.max_fragments).min(count);
            let artifact = server.session.text_artifact(bounded_json_limit(&preview, 32 * 1024 * 1024)?, "application/json", "overlap-extension-digest-inspection.json")?;
            output(json!({
                "document":args.document,"mode":"overlap-extension-pcr","engine_version":preview.engine_version,
                "context_revision":preview.context_revision,"source_revision":preview.source_revision,
                "origin_revision":preview.origin_revision,"origin_incomplete":preview.origin_incomplete,
                "selectable":preview.selectable(),"incomplete":preview.incomplete,
                "amplicon":preview.amplicon.as_ref().map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                "pcr_stages":preview.pcr_stages.iter().take(3).map(|stage|json!({
                    "stage":stage.stage,"template_revision":stage.template_revision,"geometry_count":stage.geometries.len(),
                    "product_revision":stage.product_revision,"product_length":stage.product_length
                })).collect::<Vec<_>>(),
                "extension":preview.extension.as_ref().map(|extension|json!({
                    "homology_length":extension.homology.len(),"terminal_pairings":extension.terminal_pairings.overlaps.iter().take(100).collect::<Vec<_>>(),
                    "pairings_truncated":extension.terminal_pairings.overlaps.len()>100,
                    "complete":extension.terminal_pairings.complete,"unique":extension.terminal_pairings.unique,
                    "comparisons":extension.terminal_pairings.comparisons
                })),
                "products_determined":preview.inspection.as_ref().is_some_and(|i|i.products_determined),
                "scan_complete":preview.inspection.as_ref().map(|i|i.scan_complete),
                "fragment_count":count,"fragment_start":args.fragment_start,"next_fragment_start":end,"fragments_done":end==count,
                "fragments":preview.inspection.as_ref().map(|i|&i.fragments[args.fragment_start..end]),
                "findings":preview.findings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "findings_truncated":preview.findings.len()>100||preview.findings.iter().take(100).any(|s|s.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Fresh source-bound four-step preparation and declared complete-cleavage inspection only. No target agreement, fragment isolation, stock, candidate or experimental success is inferred. Keep original-source bindings and incomplete evidence when explicitly selecting a fragment for a separate cloning recipe."
            }))
        }).await
    }
    #[tool(
        description = "Inspect actual restriction fragments of a freshly reconstructed unique exact PCR product from an immutable original source, explicit ordered primers and a full catalog pin. Primer 5′ extensions participate in the actual product and catalog scan. Ambiguous/no PCR geometry withholds physical fragment inspection; uncertain or partial cleavage retains diagnostics. Returns a paged actual end-chemistry summary and immutable full JSON evidence. Optional saved-run origin requires an explicit native workspace proving exact original run/source/target events before fresh reconstruction. Optional pcr_options enables the explicit mismatch-PCR engine with its own complete typed artifact and visible incomplete/conditional-singleton evidence; omitting it retains exact-PCR behavior. It never chooses a fragment, creates a plan, copies a result into an independent target, publishes a product document or infers amplification/cleavage success."
    )]
    pub(crate) async fn inspect_pcr_digest(
        &self,
        Parameters(args): Parameters<PcrDigestInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_fragments) {
            return Err(bad("PCR digest fragment page size must be1..100"));
        }
        if args.origin.is_some() && args.workspace.is_none() {
            return Err(bad(
                "saved-run PCR digest inspection requires its exact native workspace reference",
            ));
        }
        if let Some(pcr_options) = args.pcr_options {
            let request = collagene_core::cloning::pcr::mismatch_digest::Request {
                base: args.request,
                pcr_options,
            };
            let server = self.clone();
            return self.compute(move||{
                let workspace=args.workspace.as_ref().map(|reference|server.session.restore_workspace(reference)).transpose()?;
                let source={let store=server.session.store.lock().unwrap();let source=store.document(&args.document)?;
                    collagene_core::cloning::pcr::mismatch_digest::source_preflight(source,args.origin.as_ref(),&request,workspace.as_ref()).map_err(bad)?;source.clone()};
                let expected=collagene_core::cloning::pcr::mismatch_digest::inspection_revision(&source,args.origin.as_ref(),&request,workspace.as_ref()).map_err(bad)?;
                let preview=collagene_core::cloning::pcr::mismatch_digest::inspect(&source,args.origin.as_ref(),&request,workspace.as_ref()).map_err(bad)?;
                preview.validate_binding(&expected).map_err(bad)?;
                let count=preview.inspection.as_ref().map_or(0,|r|r.fragments.len());
                if args.fragment_start>count{return Err(bad("mismatch-PCR digest fragment offset exceeds retained evidence"));}
                let end=args.fragment_start.saturating_add(args.max_fragments).min(count);
                let result=json!({"document":args.document,"mode":"mismatch-pcr","engine_version":preview.engine_version,"context_revision":preview.context_revision,"source_revision":preview.source_revision,
                    "origin_revision":preview.origin_revision,"origin_incomplete":preview.origin_incomplete,"origin_findings_total":preview.origin_findings_total,"origin_findings_truncated":preview.origin_findings_truncated,"pcr_options":preview.request.pcr_options,"incomplete":preview.incomplete,
                    "pcr_candidate_count":preview.pcr.report.products.len(),"pcr_search_complete":preview.pcr.report.search_complete,"pcr_enumeration_complete":preview.pcr.report.enumeration_complete,
                    "pcr_candidates":preview.pcr.report.products.iter().take(100).enumerate().map(|(index,p)|json!({"index":index,"forward_start":p.forward_start,"reverse_start":p.reverse_start,"length":p.sequence.len(),"certain":p.certain,"overlapping_domains":p.overlapping_domains,"substitution_count":p.substitutions.len()})).collect::<Vec<_>>(),"pcr_candidates_truncated":preview.pcr.report.products.len()>100,
                    "amplicon":preview.amplicon.as_ref().map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                    "products_determined":preview.inspection.as_ref().is_some_and(|i|i.products_determined),"scan_complete":preview.inspection.as_ref().map(|i|i.scan_complete),
                    "fragment_count":count,"fragment_start":args.fragment_start,"next_fragment_start":end,"fragments_done":end==count,"fragments":preview.inspection.as_ref().map(|i|&i.fragments[args.fragment_start..end]),
                    "findings":preview.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":preview.findings.len()>100,
                    "authority":"Explicit mismatch-PCR policy and fresh catalog cleavage inspection only. Conditional singleton and incomplete-search evidence remain explicit. The report creates no independent target, selected cloning product, physical stock or experimental success."});
                let artifact=server.session.text_artifact(bounded_json_limit(&preview,32*1024*1024)?,"application/json","mismatch-pcr-digest-inspection.json")?;
                let mut result=result;result["artifact_id"]=json!(artifact);output(result)
            }).await;
        }
        let server = self.clone();
        self.compute(move||{
            let frame=args.workspace.as_ref().map(|reference|server.session.restore_workspace(reference)).transpose()?;
            let source={let store=server.session.store.lock().unwrap();let source=store.document(&args.document)?;
                if let Some(origin)=&args.origin{collagene_core::cloning::pcr::digest::preflight_from_material(source,origin,&args.request,frame.as_ref().ok_or_else(||bad("Saved PCR digest source workspace is unavailable"))?).map_err(bad)?;}
                else{collagene_core::cloning::pcr::digest::source_preflight(source,&args.request).map_err(bad)?;}source.clone()};
            let (ordinary,from_material)=if let Some(origin)=&args.origin {
                let frame=frame.as_ref().ok_or_else(||bad("Saved PCR digest source workspace is unavailable"))?;
                let expected=collagene_core::cloning::pcr::digest::inspection_revision_from_material(&source,origin,&args.request,frame).map_err(bad)?;
                let preview=collagene_core::cloning::pcr::digest::inspect_from_material(&source,origin,&args.request,frame).map_err(bad)?;
                preview.validate_binding(&expected).map_err(bad)?;(None,Some(preview))
            }else{
                let expected=collagene_core::cloning::pcr::digest::inspection_revision(&source,&args.request).map_err(bad)?;
                let preview=collagene_core::cloning::pcr::digest::inspect(&source,&args.request).map_err(bad)?;
                preview.validate_binding(&expected).map_err(bad)?;(Some(preview),None)
            };
            let preview=from_material.as_ref().map(|r|&r.preview).or(ordinary.as_ref()).ok_or_else(||bad("PCR digest inspection unavailable"))?;
            let context=from_material.as_ref().map(|r|r.context_revision.as_str()).unwrap_or(&preview.context_revision);
            let origin_revision=from_material.as_ref().map(|r|&r.origin_revision);
            let count=preview.inspection.as_ref().map_or(0,|i|i.fragments.len());
            if args.fragment_start>count{return Err(bad("PCR digest fragment offset exceeds retained evidence"));}
            let end=args.fragment_start.saturating_add(args.max_fragments).min(count);
            let result=json!({"document":args.document,"engine_version":preview.engine_version,"context_revision":context,"origin_revision":origin_revision,"source_revision":preview.source_revision,
                "origin_incomplete":from_material.as_ref().map(|r|r.origin_incomplete),"origin_findings_total":from_material.as_ref().map(|r|r.origin_findings_total),"origin_findings_truncated":from_material.as_ref().map(|r|r.origin_findings_truncated),
                "pcr_geometry_count":preview.pcr_geometries.len(),"pcr_geometries":preview.pcr_geometries.iter().take(100).collect::<Vec<_>>(),"pcr_geometries_truncated":preview.pcr_geometries.len()>100,
                "amplicon_length":preview.amplicon.as_ref().map(|m|m.sequence.len()),"products_determined":preview.inspection.as_ref().is_some_and(|i|i.products_determined),
                "scan_complete":preview.inspection.as_ref().map(|i|i.scan_complete),"limits":preview.inspection.as_ref().map(|i|&i.limits),
                "fragment_count":count,"fragment_start":args.fragment_start,"next_fragment_start":end,"fragments_done":end==count,
                "fragments":preview.inspection.as_ref().map(|i|&i.fragments[args.fragment_start..end]),
                "findings":preview.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":preview.findings.len()>100,
                "authority":"Fresh explicitly bound original or saved-material PCR and assumed-complete catalog digestion evidence only. Inspect actual outer PCR and cut-derived termini; no selected product, target acceptance, reaction yield or imported physical authority is created."});
            let (text,filename)=if let Some(material)=&from_material{(bounded_json_limit(material,32*1024*1024)?,"material-pcr-digest-inspection.json")}else{(bounded_json_limit(preview,32*1024*1024)?,"pcr-digest-inspection.json")};
            let artifact=server.session.text_artifact(text,"application/json",filename)?;
            let mut result=result;result["artifact_id"]=json!(artifact);output(result)
        }).await
    }
    #[tool(
        description = "Inspect an immutable DNA source with one or more exact catalog enzyme identities or uniquely resolved names. One artifact_id contains complete tagged Original or SavedRun evidence with full edition/hash/ordered-identity/geometry binding; inline pages remain bounded summaries. Preserves certain/possible recognition sites, exact origin-split spans, cuts outside substrate, bounded-scan limits and separate scan_complete/products_determined flags. Returns paged site and duplex-fragment summaries plus a complete immutable JSON artifact, bound to source role/revision, exact enzyme geometry and options. Products assume complete cleavage without blocking methylation or conditions; reports are evidence, not experimental activity or permission to release imported products."
    )]
    pub(crate) async fn inspect_restriction(
        &self,
        Parameters(args): Parameters<RestrictionInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_sites) || !(1..=100).contains(&args.max_fragments) {
            return Err(bad("restriction site/fragment page sizes must be 1..100"));
        }
        if args.enzymes.is_empty() || args.enzymes.len() > 64 {
            return Err(bad("select 1..64 embedded enzymes"));
        }
        if args.origin.is_some() && args.workspace.is_none() {
            return Err(bad(
                "saved-run catalog inspection requires its exact native workspace reference",
            ));
        }
        let catalog = catalog_tools::builtin()?;
        let ids = args
            .enzymes
            .iter()
            .map(|name| match catalog.lookup_enzyme(name) {
                collagene_enzyme_catalog::Lookup::Unique(enzyme) => Ok(enzyme.id.clone()),
                collagene_enzyme_catalog::Lookup::Missing => {
                    Err(bad("unknown enzyme identity; call lookup_enzyme_catalog"))
                }
                collagene_enzyme_catalog::Lookup::Ambiguous(_) => {
                    Err(bad("ambiguous enzyme alias; select exact catalog IDs"))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let selection =
            collagene_core::cloning::restriction::catalog::select(&catalog, &ids).map_err(bad)?;
        if let Some(pin) = &args.catalog_pin {
            pin.validate_for(&catalog).map_err(bad)?;
            if pin != &selection {
                return Err(bad("catalog pin and ordered enzyme selection differ"));
            }
        }
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                if let Some(origin)=&args.origin {collagene_core::cloning::restriction::run_digest::source_preflight(&catalog,&selection,&args.source_role,molecule,origin,&args.options).map_err(bad)?;}
                else{collagene_core::cloning::restriction::catalog::source_preflight(&catalog, &selection, &args.source_role, molecule, &args.options).map_err(bad)?;}
                molecule.clone()
            };
            let(ordinary,from_run)=if let Some(origin)=&args.origin {
                let workspace=args.workspace.as_ref().ok_or_else(||bad("saved-run catalog source workspace unavailable"))?;
                server.session.attest_run_source(workspace,origin,&molecule)?;
                let report=collagene_core::cloning::restriction::run_digest::inspect(&catalog,&selection,&args.source_role,&molecule,origin,&args.options).map_err(bad)?;
                report.validate_bindings_for(&catalog,&selection,&args.source_role,&molecule,origin,&args.options).map_err(bad)?;
                (None,Some(report))
            }else{
                let pinned=collagene_core::cloning::restriction::catalog::inspect(&catalog,&selection,&args.source_role,&molecule,&args.options).map_err(bad)?;
                pinned.validate_bindings_for(&catalog,&selection,&args.source_role,&molecule,&args.options).map_err(bad)?;(Some(pinned),None)
            };
            let pinned=from_run.as_ref().map(|r|&r.inspection).or(ordinary.as_ref()).ok_or_else(||bad("catalog inspection unavailable"))?;
            let report=&pinned.inspection;
            if args.site_start > report.sites.len() || args.fragment_start > report.fragments.len() {
                return Err(bad("restriction page offset exceeds retained evidence count"));
            }
            let site_end = args.site_start.saturating_add(args.max_sites).min(report.sites.len());
            let fragment_end = args.fragment_start.saturating_add(args.max_fragments).min(report.fragments.len());
            let mut response=json!({
                "document":args.document,"engine_version":report.engine_version,"source_role":report.source_role,
                "source_revision":report.source_revision,"source_length":report.source_length,"topology":report.topology,
                "batch_revision":report.batch_revision,"geometry_revision":report.geometry_revision,
                "options_revision":report.options_revision,"enzymes":report.enzymes,"options":report.options,
                "scan_complete":report.scan_complete,"products_determined":report.products_determined,
                "comparisons":report.comparisons,"limits":report.limits,"diagnostics":report.diagnostics,
                "ref_top_lengths":report.ref_top_lengths,"site_count":report.sites.len(),
                "site_start":args.site_start,"next_site_start":site_end,"sites_done":site_end==report.sites.len(),
                "sites":&report.sites[args.site_start..site_end],"fragment_count":report.fragments.len(),
                "fragment_start":args.fragment_start,"next_fragment_start":fragment_end,
                "fragments_done":fragment_end==report.fragments.len(),"fragments":&report.fragments[args.fragment_start..fragment_end],
                "catalog_pin":pinned.selection,"catalog_batch_revision":pinned.batch_revision,
                "catalog_engine_version":pinned.engine_version,
                "authority":"Recognition evidence and assumed-complete cleavage summaries only; possible/partial/edge evidence cannot authorize a digest product. Reconstruct selected products from current source inputs."
            });
            if let Some(run)=&from_run {
                response["context_revision"]=json!(run.context_revision);
                response["origin_revision"]=json!(run.origin_revision);
                response["origin_status"]=json!(run.origin_status);
                response["origin_findings"]=json!(run.findings);
                response["authority"]=json!("Fresh saved-run tracked material and catalog cleavage inspection only. Serialized evidence never grants material, isolation, target acceptance or permission to install a reference-only digest product.");
            }
            let complete_evidence = match (ordinary, from_run) {
                (Some(inspection), None) => RestrictionInspectionEvidence::Original {
                    inspection,
                },
                (None, Some(inspection)) => RestrictionInspectionEvidence::SavedRun {
                    inspection,
                },
                _ => return Err(bad("catalog inspection source evidence is inconsistent")),
            };
            let artifact_id = server.session.text_artifact(
                bounded_json_limit(&complete_evidence, 32 * 1024 * 1024)?,
                "application/json",
                "catalog-restriction-inspection.json",
            )?;
            response["artifact_id"] = json!(artifact_id);
            output(response)
        }).await
    }
    #[tool(
        description = "Compute restriction sites and reference-strand fragment lengths for one embedded enzyme and exact document revision."
    )]
    pub(crate) async fn analyze_restriction(
        &self,
        Parameters(args): Parameters<Restriction>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.analyze_restriction_inner(args))
            .await
    }
    pub(crate) fn analyze_restriction_inner(
        &self,
        args: Restriction,
    ) -> Result<CallToolResult, ErrorData> {
        let m = self
            .session
            .store
            .lock()
            .unwrap()
            .document(&args.document)?
            .clone();
        let enzyme = analysis::starter_enzymes()
            .into_iter()
            .find(|e| e.name.eq_ignore_ascii_case(&args.enzyme))
            .ok_or_else(|| bad("unknown enzyme; call list_enzymes"))?;
        let sites = analysis::restriction_sites(&m, &enzyme).map_err(bad)?;
        let lengths = analysis::digest_lengths(&m, &sites).map_err(bad)?;
        // Sites are bounded separately to avoid flooding a model context.
        output(
            json!({"document":args.document,"enzyme":enzyme,"site_count":sites.len(),"sites":sites.iter().take(100).collect::<Vec<_>>(),"sites_truncated":sites.len()>100,"fragment_count":lengths.len(),"fragment_lengths":lengths.iter().take(100).collect::<Vec<_>>(),"fragments_truncated":lengths.len()>100,"model":"reference-strand lengths; full duplex end chemistry and methylation are not modeled"}),
        )
    }
    #[tool(
        description = "Compute complete duplex digestion using embedded enzyme names. Reports separate top/bottom lengths, 5′/3′ overhangs and terminal chemistry in ascending top-cut order. Imported unresolved chemical evidence is not trusted; methylation/conditions are not modeled. Returns at most 100 summaries with an explicit truncation flag."
    )]
    pub(crate) async fn digest_fragments(
        &self,
        Parameters(args): Parameters<DigestQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let catalog=analysis::starter_enzymes();
            let enzymes=args.enzymes.iter().map(|name|catalog.iter().find(|e|e.name.eq_ignore_ascii_case(name)).cloned().ok_or_else(||bad("unknown enzyme; call list_enzymes"))).collect::<Result<Vec<_>,_>>()?;
            let fragments=chemistry::digest(&chemistry::Duplex::from_molecule(&molecule).map_err(bad)?,&enzymes).map_err(bad)?;
            let summaries=fragments.iter().take(100).enumerate().map(|(index,f)|f.duplex.summary().map(|summary|json!({"fragment_index":index,"top_start":f.top_start,"bottom_start":f.bottom_start,"duplex":summary}))).collect::<Result<Vec<_>,_>>().map_err(bad)?;
            output(json!({"document":args.document,"fragment_count":fragments.len(),"fragments":summaries,"truncated":fragments.len()>100,"assumption":"complete cleavage, no blocking methylation; unknown free-end chemistry remains unknown"}))
        }).await
    }
    #[tool(
        description = "Compute exact PCR products from an immutable template. Returns bounded summaries, considers competing oligo orientations and single-oligo products; mismatches/secondary structure are outside this model."
    )]
    pub(crate) async fn simulate_pcr(
        &self,
        Parameters(args): Parameters<Pcr>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.simulate_pcr_inner(args)).await
    }
    pub(crate) fn simulate_pcr_inner(&self, args: Pcr) -> Result<CallToolResult, ErrorData> {
        let m = self
            .session
            .store
            .lock()
            .unwrap()
            .document(&args.document)?
            .clone();
        let products = analysis::simulate_pcr(&m, &args.forward, &args.reverse).map_err(bad)?;
        let summaries: Vec<_> = products.iter().take(100).map(|p|json!({"length":p.sequence.len(),"sequence_hash":blake3::hash(p.sequence.as_str().as_bytes()).to_hex().to_string()})).collect();
        output(
            json!({"document":args.document,"product_count":products.len(),"products":summaries,"truncated":products.len()>100,"model":"exact non-overlapping binding sites within one template turn; no experimental guarantee"}),
        )
    }
    #[tool(
        description = "Propose explicitly selected original-source PCR and catalog-pinned Type IIS digest preparations for a junction-edited Golden Gate design. Builds independent source-core trim/insert intent before any PCR/digest. Actual tails distinguish fusion and inserted DNA, and fresh original-source PCR/cut/end geometry checks reject extra sites or nonunique binding. Retained rows, named oligos, catalog pins, treatments and repairs stay exact. This is an advisory ordered preparation/ligation recipe, not a cycling, domestication, fusion-optimization or experimental result. Build the edited intent independently and freshly check the returned complete recipe using propose_ligation_assembly."
    )]
    pub(crate) async fn design_golden_gate_junction_primers(
        &self,
        Parameters(args): Parameters<GoldenGateJunctionPrimersQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let inputs={let store=server.session.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::cloning::golden_gate::junction_edit::preflight(&sources,&args.request).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::cloning::golden_gate::junction_edit::design(&inputs,&args.request).map_err(bad)?;
            let expected=collagene_core::cloning::golden_gate::junction_edit::design_revision(&inputs,&args.request).map_err(bad)?;report.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"report":report});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","golden-gate-junction-primers.json")?;
            let ordered=report.recipe.fragments.iter().zip(&report.request.rows).filter_map(|(fragment,row)|{let f=&fragment.source;match row {
                collagene_core::cloning::golden_gate::junction_edit::Row::Retain(_) => None,
                collagene_core::cloning::golden_gate::junction_edit::Row::DesignPcr{..} => Some(match &f.preparation {
                    collagene_core::cloning::assembly::recipe::Preparation::PcrCatalogDigest{forward,reverse,..}=>Ok(json!({"role":f.role,"input_id":f.input_id,"reverse_fragment":f.reverse,"forward":forward,"reverse":reverse,"forward_ordered":forward.ordered_sequence(),"reverse_ordered":reverse.ordered_sequence()})),
                    _=>Err(bad("selected Golden Gate proposal lost its original-source PCR/digest chain"))
                })
            }}).collect::<Result<Vec<_>,_>>()?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,"intent_context_revision":report.intent.context_revision,
                "target_revision":report.intent.target.revision(),"recipe":report.recipe,"ordered_oligos":ordered,"design_count":ordered.len(),"retained_rows":report.retained_rows,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Build an unverified intended reference by explicitly replacing, inserting or deleting a single source interval. Uses only immutable source DNA and the requested edit, independently of any PCR or cloning simulator. Supports linear DNA and at most one circular turn with explicit origin shift. Atomically retains the intended target and full source-bound JSON. Does not mutate the source or claim physical material, successful cloning or observed DNA; design and verify inverse-PCR or another preparation separately against this frozen reference."
    )]
    pub(crate) async fn build_mutation_intent(
        &self,
        Parameters(args): Parameters<MutationIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let source = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_core::cloning::pcr::mutation_intent::source_preflight(source,&args.request).map_err(bad)?;
                source.clone()
            };
            let report = collagene_core::cloning::pcr::mutation_intent::build(&source,&args.request).map_err(bad)?;
            let expected = collagene_core::cloning::pcr::mutation_intent::intent_revision(&source,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a> { document:&'a DocumentRef, report:&'a collagene_core::cloning::pcr::mutation_intent::Report }
            let evidence = Evidence{document:&args.document,report:&report};
            let (artifact, target) = server.session.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"mutation-intent.json",32*1024*1024)?;
            let target = target.ok_or_else(||bad("intended mutation reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,
                "document":args.document,"context_revision":report.context_revision,"request":report.request,
                "source_revision":report.source_revision,"target":{"id":target["id"],"revision":target["revision"]},
                "target_summary":target,"deleted_length":report.deleted.len(),"source_origin_shift":report.origin_shift,
                "warnings":report.warnings,"artifact_id":artifact,
                "authority":"Declared reference edit only. The source is unchanged; physical preparation, target agreement and experimental success require separate evidence."}))
        }).await
    }
    #[tool(
        description = "Plan bounded anchored PCR primers for explicitly ordered source intervals and orientations, using exact immutable DNA documents and an independently supplied target. Generates one-sided homology tails from the next actual source interval, never from target DNA, plus a practical overlap recipe and ordered oligos. Existing primer heuristic/conditions and Gibson/HiFi/InFusion advisory profile are explicit. Planner output is advisory and unverified, not a saved candidate or proof of amplification, specificity or target agreement. Pass its recipe and exact inputs/target to propose_overlap_assembly for fresh preparation and verification."
    )]
    pub(crate) async fn design_overlap_primers(
        &self,
        Parameters(args): Parameters<OverlapPrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad(
                "overlap-primer design requires 1..128 actual source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target)={let store=server.session.store.lock().unwrap();let target=store.document(&args.target)?;
                let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()) {return Err(bad("empty, oversized or duplicate template input ID"));}
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::cloning::overlap::primer_design::preflight(&sources,target,&args.request).map_err(bad)?;
                let inputs=sources.into_iter().map(|(id,molecule)|PlanInput{id:id.into(),revision:molecule.revision(),molecule:molecule.clone()}).collect::<Vec<_>>();
                (inputs,target.clone())};
            let report=collagene_core::cloning::overlap::primer_design::design(&inputs,&target,&args.request).map_err(bad)?;
            let expected=collagene_core::cloning::overlap::primer_design::design_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let ordered_oligos=report.recipe.fragments.iter().map(|fragment| {
                match &fragment.preparation {
                    collagene_core::cloning::assembly::recipe::Preparation::ExactPcr{forward,reverse, .. }=>Ok(json!({
                        "role":fragment.role,"input_id":fragment.input_id,"reverse_fragment":fragment.reverse,
                        "forward":{"primer":forward,"ordered_sequence":forward.ordered_sequence()},
                        "reverse":{"primer":reverse,"ordered_sequence":reverse.ordered_sequence()}})),
                    _=>Err(bad("planner did not retain explicit exact-PCR preparation"))
                }
            }).collect::<Result<Vec<_>,_>>()?;
            #[derive(Serialize)]
            struct Evidence<'a>{inputs:&'a [BoundInput],target:&'a DocumentRef,report:&'a collagene_core::cloning::overlap::primer_design::Report,interpretation:&'static str}
            let evidence=Evidence{inputs:&args.inputs,target:&args.target,report:&report,
                interpretation:"Advisory anchored primer heuristic and source-derived one-sided homology only. Independent target identity is retained but agreement, amplification and uniqueness require fresh compiler/verifier reconstruction; no source or candidate is promoted."};
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","overlap-primer-design.json")?;
            output(json!({"state":"advisory","verified":false,"inputs":args.inputs,"target":args.target,
                "context_revision":report.context_revision,"source_revisions":report.source_revisions,
                "recipe":report.recipe,"ordered_oligos":ordered_oligos,
                "design_count":report.designs.len(),"designs":report.designs.iter().take(32).map(|design|json!({
                    "source_revision":design.source_revision,"candidate_pair_count":design.pairs.len(),
                    "evaluated_pairs":design.evaluated_pairs,"rejected_pairs":design.rejected_pairs,
                    "warnings":design.warnings.iter().take(10).map(|w|w.chars().take(512).collect::<String>()).collect::<Vec<_>>(),
                    "warnings_truncated":design.warnings.len()>10||design.warnings.iter().take(10).any(|w|w.chars().count()>512)})).collect::<Vec<_>>(),
                "design_preview_truncated":report.designs.len()>32,
                "warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|w|w.chars().count()>1024),
                "artifact_id":artifact,"authority":"Advisory design only. Use the returned recipe with propose_overlap_assembly against the same exact inputs and independent target; method guidance and primer heuristics do not verify a cloning result."}))
        }).await
    }
    #[tool(
        description = "Propose primers for an entire explicitly edited overlap assembly from immutable original source documents and a typed independent contribution/junction/substitution edit intent. Unedited rows use exact PCR in their explicitly bound original or native strand frame; internally edited rows use actual two-branch overlap-extension PCR. Outgoing tails come from literal junction inserts plus the next edited source-core prefix, never a supplied product sequence. Freshly verifies the whole composed recipe against the separately built edit intent and retains failures/uncertainty. The single intent material-binding roster contains compact exact input/event Run/Derived selections, including actual Top/Bottom strands; all bindings require one exact origin_workspace. Closed recipes are hydrated jointly without replay, then the native proposal freshly reconstructs them. Returns a read-only advisory recipe, ordered oligos and bounded status summaries with full evidence; no source, target document, candidate or material status is installed. Retain the intended target separately with build_contribution_edit_intent using the same intent, then use propose_overlap_assembly and normal review/finalization."
    )]
    pub(crate) async fn design_edited_overlap_primers(
        &self,
        Parameters(args): Parameters<EditedOverlapPrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 exact sources"));
        }
        bounded_json_limit(&args.request, 2 * 1024 * 1024)?;
        if !args.request.intent.base.base.material_bindings.is_empty()
            && args.origin_workspace.is_none()
        {
            return Err(bad(
                "material source bindings require one exact origin_workspace",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=args.origin_workspace.as_ref().map(|wr|server.session.restore_workspace(wr)).transpose()?;
            let request=collagene_core::cloning::assembly::edited_primers::Request {
                intent:args.request.intent.native(frame.as_ref())?,overlap:args.request.overlap,
                options:args.request.options.clone(),profile_options:args.request.profile_options.clone(),finishing:args.request.finishing.clone(),
            };
            let inputs = {
                let store = server.session.store.lock().unwrap();
                let mut ids = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {
                    if input.id.trim().is_empty() || input.id.len() > 128 || !ids.insert(input.id.as_str()) {
                        return Err(bad("invalid or duplicate original source ID"));
                    }
                    sources.push((input.id.as_str(), store.document(&input.document)?));
                }
                // Native source_preflight jointly attests the single input/event
                // roster against borrowed immutable session documents before cloning.
                collagene_core::cloning::assembly::edited_primers::source_preflight(&sources, &request, frame.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,molecule)|PlanInput { id:id.into(),revision:molecule.revision(),molecule:molecule.clone() }).collect::<Vec<_>>()
            };
            let expected = collagene_core::cloning::assembly::edited_primers::design_revision(&inputs, &request, frame.as_ref()).map_err(bad)?;
            let report = collagene_core::cloning::assembly::edited_primers::design(&inputs, &request, frame.as_ref()).map_err(bad)?;
            report.validate_sources(&inputs, &expected, frame.as_ref()).map_err(bad)?;
            let mut ordered = Vec::new();
            for fragment in &report.recipe.fragments {
                use collagene_core::cloning::assembly::recipe::{Preparation, RunPreparation};
                let pairs:Vec<(&str,&Primer,&Primer,collagene_core::chemistry::Strand)> = match &fragment.preparation {
                    Preparation::ExactPcr {forward,reverse,template_strand} => vec![("exact_pcr",forward,reverse,*template_strand)],
                    Preparation::OverlapExtensionPcr {request} => vec![("left_pcr",&request.left.forward,&request.left.reverse,request.left.template_strand),("right_pcr",&request.right.forward,&request.right.reverse,request.right.template_strand),("outer_pcr",&request.outer.forward,&request.outer.reverse,request.outer.template_strand)],
                    Preparation::FromRun {preparation:RunPreparation::ExactPcr {forward,reverse,template_strand},..} | Preparation::FromDerived {preparation:RunPreparation::ExactPcr {forward,reverse,template_strand},..} => vec![("exact_pcr",forward,reverse,*template_strand)],
                    Preparation::FromRun {preparation:RunPreparation::OverlapExtensionPcr {request},..} | Preparation::FromDerived {preparation:RunPreparation::OverlapExtensionPcr {request},..} => vec![("left_pcr",&request.left.forward,&request.left.reverse,request.left.template_strand),("right_pcr",&request.right.forward,&request.right.reverse,request.right.template_strand),("outer_pcr",&request.outer.forward,&request.outer.reverse,request.outer.template_strand)],
                    _ => return Err(bad("edited assembly primer proposal lost its explicit exact/OE PCR preparation")),
                };
                for (stage,forward,reverse,template_strand) in pairs {
                    for (direction,primer) in [("forward",forward),("reverse",reverse)] {
                        ordered.push(json!({"role":fragment.role,"input_id":fragment.input_id,"reverse_fragment":fragment.reverse,"stage":stage,"template_strand":template_strand,"direction":direction,"primer":primer,"ordered_sequence":primer.ordered_sequence()}));
                    }
                }
            }
            #[derive(Serialize)]
            struct Evidence<'a> { inputs:&'a [BoundInput],origin_workspace:&'a Option<WorkspaceRef>,request:&'a EditedOverlapPrimerRequest,report:&'a collagene_core::cloning::assembly::edited_primers::Report,interpretation:&'static str }
            let evidence = Evidence { inputs:&args.inputs,origin_workspace:&args.origin_workspace,request:&args.request,report:&report,
                interpretation:"Read-only whole-construct primer proposal. Independent intended DNA is built only from exact original source contributions and explicit edit declarations. Primer tails and selected native source geometry never select or replace that intended target; whole-recipe reconstruction remains diagnostic until separate review/finalization." };
            let artifact = server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","edited-overlap-primer-design.json")?;
            output(json!({
                "state":"advisory","verified":false,"inputs":args.inputs,"origin_workspace":args.origin_workspace,
                "engine_version":report.engine_version,"context_revision":report.context_revision,"source_revisions":report.source_revisions,
                "intent_context_revision":report.intent.context_revision,
                "target_summary":{"revision":report.intent.target.revision(),"length":report.intent.target.sequence.len(),"topology":report.intent.target.topology},
                "recipe":report.recipe,"ordered_oligos":ordered,"tails":report.tails,
                "verification":{"status":report.preview.verification.status,"plan_revision":report.preview.verification.plan_revision,
                    "target_revision":report.preview.verification.target_revision,"product_count":report.preview.verification.products.len(),
                    "findings":report.preview.verification.findings.iter().take(50).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                    "findings_truncated":report.preview.verification.findings.len()>50||report.preview.verification.findings.iter().take(50).any(|s|s.chars().count()>1024)},
                "rows":report.rows.iter().take(32).map(|row|json!({"role":row.role,"input_id":row.input_id,"input_revision":row.input_revision,
                    "edit_count":row.edits.len(),"design_count":row.designs.len(),"template_strand":row.template_strand,"search":row.search})).collect::<Vec<_>>(),
                "warnings":report.warnings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|s|s.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Read-only advisory recipe. Retain the independently declared target with build_contribution_edit_intent using the same inputs/request.intent, then freshly verify this recipe against that target through propose_overlap_assembly. No source, target document, candidate or physical material is installed."
            }))
        }).await
    }
    #[tool(
        description = "Propose six ordered oligos for a source-bound overlap-extension PCR preparation from an immutable original or exact saved Run/Derived source, an explicit stored-strand interval/orientation and independently declared literal edits. Saved sources require an exact material binding and owning origin_workspace; actual strand frames are reconstructed natively. Uses anchored primer heuristics on unchanged native source domains, then freshly reconstructs two PCR branches, selected reciprocal extension and outer-primer PCR against the separately constructed edit intent. Returns an advisory preparation, bounded oligo/status summaries and full source/request/intent/reconstruction evidence. Prepared DNA stays in stored-source orientation; apply request.fragment.reverse at subsequent assembly. No intended product is accepted as a source, no stock is changed, and no candidate, target document or material status is installed. Physical-state metadata and historical reports cannot bootstrap an original PCR template. Reuse the returned preparation in the ordinary assembly recipe and verify the full independent target."
    )]
    pub(crate) async fn design_overlap_extension_primers(
        &self,
        Parameters(args): Parameters<OverlapExtensionPrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let workspace = args.origin_workspace.as_ref().map(|reference| server.session.restore_workspace(reference)).transpose()?;
            let source = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_core::cloning::overlap::extension::primers::preflight(source, &args.request, workspace.as_ref()).map_err(bad)?;
                source.clone()
            };
            let report = collagene_core::cloning::overlap::extension::primers::design(&source, &args.request, workspace.as_ref()).map_err(bad)?;
            let expected = collagene_core::cloning::overlap::extension::primers::design_revision(&source, &args.request, workspace.as_ref()).map_err(bad)?;
            report.validate_sources(&source, &expected, workspace.as_ref()).map_err(bad)?;
            let mut ordered_oligos = Vec::with_capacity(6);
            for (stage, pair) in [("left_pcr", &report.preparation.left), ("right_pcr", &report.preparation.right), ("outer_pcr", &report.preparation.outer)] {
                for (direction, primer) in [("forward", &pair.forward), ("reverse", &pair.reverse)] {
                    ordered_oligos.push(json!({"stage":stage,"template_strand":pair.template_strand,"direction":direction,"primer":primer,"ordered_sequence":primer.ordered_sequence()}));
                }
            }
            #[derive(Serialize)]
            struct Evidence<'a> {
                document: &'a DocumentRef,
                report: &'a collagene_core::cloning::overlap::extension::primers::Report,
                interpretation: &'static str,
            }
            let evidence = Evidence { document: &args.document, report: &report,
                interpretation: "Read-only native primer proposal from explicit immutable source coordinates and independent literal edits. Expected DNA comes from edit intent, not simulated products. The retained local reconstruction is diagnostic; no stock, candidate, target document or material authority is installed. Reuse the preparation in a full independently targeted assembly before review/finalization." };
            let artifact = server.session.text_artifact(bounded_json_limit(&evidence, 32 * 1024 * 1024)?, "application/json", "overlap-extension-primer-design.json")?;
            output(json!({
                "state":"advisory", "verified":false, "document":args.document,
                "context_revision":report.context_revision,"source_revision":report.source_revision,"engine_version":report.engine_version,
                "fragment":report.request.fragment,"reverse_at_assembly":report.request.fragment.reverse,
                "preparation":report.preparation,"ordered_oligos":ordered_oligos,
                "expected":{"revision":report.expected.revision(),"length":report.expected.sequence.len(),"topology":report.expected.topology,"orientation":"stored_source"},
                "verification":{"status":report.verification.status,"plan_revision":report.verification.plan_revision,"target_revision":report.verification.target_revision,
                    "product_count":report.verification.products.len(),"products":report.verification.products.iter().take(4).map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})).collect::<Vec<_>>(),
                    "findings":report.verification.findings.iter().take(30).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                    "findings_truncated":report.verification.findings.len()>30||report.verification.findings.iter().take(30).any(|s|s.chars().count()>1024)},
                "design_count":report.designs.len(),
                "search":report.search.iter().take(3).map(|search|json!({
                    "call":search.call,"candidates":search.candidates,"examined":search.examined,"rejected":search.rejected,
                    "selected_index":search.selected_index,"capped":search.capped,"unique_geometry":search.unique_geometry,
                    "notes":search.notes.iter().take(10).map(|s|s.chars().take(512).collect::<String>()).collect::<Vec<_>>(),
                    "notes_truncated":search.notes.len()>10||search.notes.iter().take(10).any(|s|s.chars().count()>512)
                })).collect::<Vec<_>>(),
                "warnings":report.warnings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|s|s.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Advisory source-bound proposal only. The source and independent edits remain unchanged; use its explicit preparation in a full intended assembly for fresh native verification and review."
            }))
        }).await
    }
    #[tool(
        description = "Propose Type IIS site-bearing primers from explicit original source intervals, orientations, domain lengths and right fusion ends. Supports actual BsaI/BsmBI/BbsI/SapI geometry, including three-base SapI ends. Clamps and cleavage spacers are explicit DNA; target DNA does not select primers. Returns original-source PCR followed by pinned catalog digestion preparations and inspectable oligo segments. Advisory proposals do not install a candidate or establish experimental success or finite-mixture uniqueness. Pass fragments to propose_ligation_assembly for fresh reconstruction against the same independent target."
    )]
    pub(crate) async fn design_golden_gate_primers(
        &self,
        Parameters(args): Parameters<GoldenGatePrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad(
                "Type IIS primer design requires 1..128 original source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs, target) = {
                let store = server.session.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len() > 128 || !seen.insert(bound.id.as_str()) { return Err(bad("empty, oversized or duplicate source ID")); }
                    sources.push((bound.id.as_str(), store.document(&bound.document)?));
                }
                collagene_core::cloning::golden_gate::primers::source_preflight(&sources, target, &args.request).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let report = collagene_core::cloning::golden_gate::primers::design(&inputs,&target,&args.request).map_err(bad)?;
            let expected = collagene_core::cloning::golden_gate::primers::preview_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a> { inputs:&'a [BoundInput], target:&'a DocumentRef, report:&'a collagene_core::cloning::golden_gate::primers::Report }
            let evidence = Evidence{inputs:&args.inputs,target:&args.target,report:&report};
            let artifact = server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","golden-gate-primer-design.json")?;
            output(json!({"state":"advisory","verified":false,"inputs":args.inputs,"target":args.target,
                "context_revision":report.context_revision,"request":report.request,"fragments":report.fragments,
                "rows":report.rows.iter().map(|row|json!({
                    "role":row.role,"input_id":row.input_id,"input_revision":row.input_revision,
                    "forward":row.forward,"reverse":row.reverse,"forward_tail":row.forward_tail,"reverse_tail":row.reverse_tail,
                    "selection":row.selection,"fragment_index":row.fragment_index,
                    "oriented_length":row.oriented_top.len(),"oriented_top_preview":row.oriented_top.chars().take(160).collect::<String>(),
                    "oriented_top_truncated":row.oriented_top.len()>160,"left":row.left,"right":row.right,
                })).collect::<Vec<_>>(),"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Explicit source-derived oligo proposals and preparation requests. Fresh ordered ligation verification against the independent target is separate; no claim of reaction yield or exhaustive mixture uniqueness."}))
        }).await
    }
}
