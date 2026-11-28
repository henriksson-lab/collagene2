//! Request ownership and outbound validation.
use super::{
    envelope::{MAX_WORKSPACE_MESSAGE_BYTES, encode_with_limit},
    schema::AnalysisJob,
};
use crate::formats::MAX_IMPORT_BYTES;
use collagene_model::Error;
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Request {
    pub id: u32,
    pub revision: String,
    pub job: AnalysisJob,
}

impl Request {
    pub fn encode(&self) -> Result<Vec<u8>, Error> {
        if let AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            ..
        } = &self.job
        {
            crate::pairwise_export::preflight(record, reference_name, query_name)?;
        }
        if let AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } = &self.job
        {
            crate::pairwise_search::preflight(report, pattern, options)?;
        }
        if let AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } = &self.job
        {
            crate::pairwise_context::preflight(report, reference, query, *translate)?;
        }
        match &self.job {
            AnalysisJob::InversePcr {
                molecule,
                target,
                origin,
                options,
            } => {
                crate::mutagenesis::inverse_pcr_revision(
                    molecule,
                    target,
                    origin.as_ref(),
                    options,
                )?;
            }

            AnalysisJob::PreviewProtocolLayout { workspace, request }
            | AnalysisJob::CaptureProtocolLayout { workspace, request } => {
                crate::protocol_layout::preview_revision(workspace, request)?;
            }
            AnalysisJob::PreviewProtocolQuantities { workspace, request }
            | AnalysisJob::CaptureProtocolQuantities { workspace, request } => {
                crate::protocol_quantities::preview_revision(workspace, request)?;
            }
            AnalysisJob::ListTopoGraphs { workspace } => {
                crate::topo_graph_history::list_revision(workspace)?;
            }
            AnalysisJob::CaptureTopoGraph { workspace, request } => {
                crate::topo_graph_history::preflight(workspace, request)?;
            }
            AnalysisJob::InspectTopoGraph {
                workspace,
                selection,
            } => {
                crate::topo_graph_history::inspect_revision(workspace, selection)?;
            }
            AnalysisJob::ListProtocolQuantities { workspace } => {
                crate::protocol_quantities::list_revision(workspace)?;
            }
            AnalysisJob::OpenProtocolQuantities {
                workspace,
                quantity_hash,
            } => {
                crate::protocol_quantities::open_revision(workspace, quantity_hash)?;
            }
            AnalysisJob::ExportProtocolQuantities {
                workspace,
                quantity_hash,
                format,
            } => {
                crate::protocol_quantities::export_revision(workspace, quantity_hash, *format)?;
            }

            AnalysisJob::ExportProtocolLayout {
                workspace,
                layout_hash,
                format,
            } => {
                crate::protocol_layout::export_revision(workspace, layout_hash, *format)?;
            }

            AnalysisJob::InspectMaterialOrigin {
                workspace,
                event,
                selection,
            } => {
                crate::run_material::material_inspection_revision(workspace, event, selection)?;
            }
            AnalysisJob::MaterialPcrDocuments {
                molecule,
                origin,
                forward,
                reverse,
                template_strand,
            } => {
                crate::material_pcr::revision(
                    molecule,
                    origin,
                    forward,
                    reverse,
                    *template_strand,
                )?;
            }
            AnalysisJob::PrepareMaterialPcrOutputs { origin, request } => {
                crate::material_pcr::outputs_revision(origin, request)?;
            }
            AnalysisJob::ReviewSequenceExport { molecule, format }
            | AnalysisJob::PrepareSequenceExport {
                molecule, format, ..
            } => {
                crate::sequence_export::revision(molecule, format)?;
            }
            AnalysisJob::ExportScfTrace { trace, options } => {
                crate::scf_export::revision(trace, options)?;
            }
            AnalysisJob::InspectMaterialMismatchPcr {
                source,
                origin,
                forward,
                reverse,
                options,
            } => {
                crate::material_derivation::mismatch_pcr_preview_revision(
                    source, origin, forward, reverse, options,
                )?;
            }
            AnalysisJob::InspectMaterialRepair {
                source,
                origin,
                method,
            } => {
                crate::material_derivation::repair_preview_revision(source, origin, *method)?;
            }
            AnalysisJob::InspectMaterialOverlapExtension {
                workspace,
                source,
                origin,
                preparation,
            } => {
                crate::material_derivation::overlap_extension_preview_revision(
                    workspace,
                    source,
                    origin,
                    preparation,
                )?;
            }
            AnalysisJob::InspectMaterialTerminalProcessing {
                workspace,
                source,
                origin,
                action,
            } => {
                crate::material_derivation::terminal_processing_preview_revision(
                    workspace, source, origin, action,
                )?;
            }
            AnalysisJob::InspectMaterialAnnealedOligos {
                workspace,
                top,
                bottom,
                options,
            } => {
                crate::material_derivation::annealed_oligos_preview_revision(
                    workspace, top, bottom, options,
                )?;
            }
            AnalysisJob::AnalyzeAlignmentAntibody {
                alignment,
                source,
                options,
                limits,
                selection,
            } => {
                crate::alignment_antibody::request_revision(
                    alignment, source, options, limits, selection,
                )?;
            }
            AnalysisJob::AnalyzeAntibody {
                molecule,
                options,
                limits,
                selection,
            } => {
                crate::antibody_numbering::domain_tracks::request_revision(
                    molecule, options, limits, selection,
                )
                .map_err(|e| Error::Format(e.to_string()))?;
            }
            AnalysisJob::ScreenProteinLiabilities { molecule, options } => {
                crate::protein_liabilities::request_revision(molecule, options)
                    .map_err(|e| Error::Format(e.to_string()))?;
            }
            AnalysisJob::AnalyzeRnaChemistry {
                molecule,
                declaration,
            } => {
                crate::rna_chemistry::declaration_revision(molecule, declaration)
                    .map_err(|e| Error::Format(e.to_string()))?;
            }
            AnalysisJob::PrepareStructure {
                bytes,
                format,
                selection,
            } => {
                crate::structure_view::revision(bytes, *format, selection)?;
            }
            AnalysisJob::FindSimilarSequences {
                query,
                subjects,
                options,
            } => {
                crate::sequence_similarity::preflight(query, subjects, options)?;
            }
            AnalysisJob::InstallMaterialDerivation { workspace, request } => {
                crate::material_derivation::preflight(workspace, request)?;
            }
            AnalysisJob::PrepareSavedRunReview { workspace, run_id } => {
                crate::run_review::request_revision(workspace, run_id)?;
            }
            AnalysisJob::RenderMapSvg { request } => request.validate()?,
            AnalysisJob::InspectMaterialDigest {
                source_role,
                source,
                origin,
                selection,
                options,
            } => {
                crate::material_derivation::digest_inspection_revision(
                    source_role,
                    source,
                    origin,
                    selection,
                    options,
                )?;
            }

            AnalysisJob::InspectTopoRun {
                inputs,
                template,
                task,
            } => {
                crate::topo_run_preparation::preflight(inputs, template.as_ref(), task)?;
            }
            AnalysisJob::PreviewAnnealingRecipe {
                inputs,
                target,
                recipe,
                policy,
            } => {
                if inputs.len() > crate::assembly_recipe::MAX_ROWS {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|input| (input.id.as_str(), &input.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_recipe::annealing_source_preflight(
                    &sources, target, recipe, policy,
                )?;
            }
            AnalysisJob::PreviewHomologyRecipe {
                inputs,
                target,
                recipe,
                policy,
            } => {
                if inputs.len() > crate::assembly_recipe::MAX_ROWS {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|input| (input.id.as_str(), &input.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_recipe::homology_source_preflight(
                    &sources, target, recipe, policy,
                )?;
            }
            AnalysisJob::DesignEditedOverlapPrimers {
                inputs,
                request,
                workspace,
            } => {
                if inputs.len() > 32 {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|input| (input.id.as_str(), &input.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_edited_primer_design::source_preflight(
                    &sources,
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::DesignOverlapExtensionPrimers {
                source,
                request,
                workspace,
            } => {
                crate::overlap_extension_primer_design::preflight(
                    source,
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::InspectOverlapExtensionDigest {
                source,
                origin,
                request,
                workspace,
            } => {
                crate::overlap_extension_digest_inspection::source_preflight(
                    source,
                    origin.as_ref(),
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::InspectMismatchPcrDigest {
                molecule,
                origin,
                request,
                workspace,
            } => {
                crate::mismatch_pcr_digest_inspection::source_preflight(
                    molecule,
                    origin.as_ref(),
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::BuildAssemblyIntent {
                inputs,
                request,
                workspace,
            } => {
                if inputs.len() > 32 {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_intent::preflight(&sources, request, workspace.as_ref())?;
            }
            AnalysisJob::BuildEditedAssemblyIntent {
                inputs,
                base,
                edits,
                workspace,
            } => {
                if inputs.len() > 32 {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_junction_edit::intent_preflight(
                    &sources,
                    base,
                    edits,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::BuildContributionEditIntent {
                inputs,
                request,
                workspace,
            } => {
                if inputs.len() > 32 {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_edit_intent::source_preflight(
                    &sources,
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::BuildSubstitutedAssemblyIntent {
                inputs,
                request,
                workspace,
            } => {
                if inputs.len() > 32 {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>();
                crate::assembly_substitution_intent::source_preflight(
                    &sources,
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::InspectRunCatalogRestriction {
                source_role,
                molecule,
                origin,
                selection,
                options,
            } => {
                let catalog = crate::catalog_digest::builtin_catalog()?;
                crate::run_catalog_digest::source_preflight(
                    &catalog,
                    selection,
                    source_role,
                    molecule,
                    origin,
                    options,
                )?;
            }
            AnalysisJob::InspectMaterialPcrDigest {
                molecule,
                origin,
                request,
                workspace,
            } => {
                crate::pcr_digest_inspection::preflight_from_material(
                    molecule, origin, request, workspace,
                )?;
            }
            AnalysisJob::InspectRunMaterial {
                workspace,
                event,
                selection,
            } => {
                crate::run_material::inspection_revision(workspace, event, selection)?;
            }
            AnalysisJob::InspectPcrDigest { molecule, request } => {
                crate::pcr_digest_inspection::source_preflight(molecule, request)?;
            }
            AnalysisJob::PreviewTailingRecipe {
                inputs,
                target,
                recipe,
            } => {
                crate::tailing_recipe::preflight(inputs, target, recipe)?;
            }
            AnalysisJob::InspectCombinatorialMultiSite {
                workspace,
                design,
                construct,
            } => {
                crate::combinatorial_design::inspect_multisite_revision(
                    design, construct, workspace,
                )?;
            }
            AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
                crate::multisite_recombination::preflight(inputs, request)?;
            }
            AnalysisJob::ProjectMultiSiteRecombination {
                inputs,
                request,
                selection,
            } => {
                crate::multisite_recombination::product_revision(inputs, request, selection)?;
            }
            AnalysisJob::CheckMultiSiteRecombination {
                inputs,
                request,
                selection,
                target,
            } => {
                crate::multisite_recombination::check_revision(inputs, request, selection, target)?;
            }
            AnalysisJob::ProcessTopoActions { inputs, request } => {
                crate::topo_actions::preflight(inputs, request)?
            }

            AnalysisJob::ListLigationProfiles => {
                crate::ligation_profiles::catalog_revision()?;
            }
            AnalysisJob::ScorePublishedLigation {
                profile,
                ends,
                physical_competitors,
            } => {
                crate::ligation_profiles::score_revision_with_competitors(
                    profile,
                    ends,
                    physical_competitors,
                )?;
            }
            AnalysisJob::ScoreLigationCompetition {
                dataset,
                ends,
                physical_competitors,
            } => {
                crate::ligation_fidelity::competition_revision(
                    dataset,
                    ends,
                    physical_competitors,
                )?;
            }
            AnalysisJob::RankGoldenGateFusions { inputs, request } => {
                crate::golden_gate_fusion_rank::preflight(
                    &inputs
                        .iter()
                        .map(|i| (i.id.as_str(), &i.molecule))
                        .collect::<Vec<_>>(),
                    request,
                )?;
            }
            AnalysisJob::InspectDomestication { source, request }
            | AnalysisJob::ProposeDomestication { source, request }
            | AnalysisJob::CheckDomestication {
                source, request, ..
            } => {
                crate::domestication::preflight(source, request)?;
            }
            AnalysisJob::PrepareDomestication {
                input,
                request,
                workspace,
            } => {
                crate::domestication_preparation::preflight(input, request, workspace.as_ref())?;
            }
            AnalysisJob::SearchGoldenGateFusions { inputs, request }
            | AnalysisJob::MaterializeGoldenGateFusions {
                inputs, request, ..
            } => {
                crate::golden_gate_fusion_design::preflight(
                    &inputs
                        .iter()
                        .map(|input| (input.id.as_str(), &input.molecule))
                        .collect::<Vec<_>>(),
                    request,
                )?;
            }
            AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
                crate::ligation_fidelity::revision(dataset, ends)?;
            }
            AnalysisJob::BeginCombinatorialBatch { workspace, matrix } => {
                crate::combinatorial_batch::begin_revision(workspace, matrix)?;
            }
            AnalysisJob::CaptureCombinatorialBatch {
                workspace,
                current_context,
            } => {
                crate::combinatorial_batch::capture_revision(workspace, current_context)?;
            }
            AnalysisJob::BeginCombinatorialRelease { workspace, request } => {
                crate::combinatorial_release::preflight(workspace, request)?;
            }
            AnalysisJob::ExportCombinatorialOutputs {
                workspace,
                manifest_hash,
                options,
            } => {
                crate::combinatorial_export::preflight(workspace, manifest_hash, options)?;
            }
            _ => {}
        }
        let maximum = if matches!(
            self.job,
            AnalysisJob::SaveWorkspace { .. }
                | AnalysisJob::SimulateTopo { .. }
                | AnalysisJob::InspectTopoOrientations { .. }
                | AnalysisJob::ProjectTopoProduct { .. }
                | AnalysisJob::CheckTopoTarget { .. }
                | AnalysisJob::ProcessTopoActions { .. }
                | AnalysisJob::ListTopoGraphs { .. }
                | AnalysisJob::CaptureTopoGraph { .. }
                | AnalysisJob::InspectTopoGraph { .. }
                | AnalysisJob::BuildAssemblyIntent { .. }
                | AnalysisJob::BuildEditedAssemblyIntent { .. }
                | AnalysisJob::BuildSubstitutedAssemblyIntent { .. }
                | AnalysisJob::BuildContributionEditIntent { .. }
                | AnalysisJob::LoadWorkspace { .. }
                | AnalysisJob::RecordCloningRun { .. }
                | AnalysisJob::ReplayCloningRun { .. }
                | AnalysisJob::LoadRunEvidence { .. }
                | AnalysisJob::FinalizeCloningOutputs { .. }
                | AnalysisJob::CaptureReadCheck { .. }
                | AnalysisJob::ReplayReadCheck { .. }
                | AnalysisJob::PreviewReadReview { .. }
                | AnalysisJob::CommitReadReview { .. }
                | AnalysisJob::InstallDigestFragment { .. }
                | AnalysisJob::InstallCatalogDigestFragment { .. }
                | AnalysisJob::BeginCombinatorialBatch { .. }
                | AnalysisJob::CaptureCombinatorialBatch { .. }
                | AnalysisJob::BeginCombinatorialRelease { .. }
                | AnalysisJob::FinishCombinatorialRelease { .. }
                | AnalysisJob::ExportCombinatorialOutputs { .. }
                | AnalysisJob::InspectCombinatorialMultiSite { .. }
                | AnalysisJob::InspectRunMaterial { .. }
                | AnalysisJob::PreviewAnnealingRecipe { .. }
                | AnalysisJob::DesignEditedOverlapPrimers { .. }
                | AnalysisJob::DesignOverlapExtensionPrimers { .. }
                | AnalysisJob::InspectMaterialPcrDigest { .. }
                | AnalysisJob::InspectMismatchPcrDigest { .. }
                | AnalysisJob::InspectOverlapExtensionDigest { .. }
                | AnalysisJob::InspectTopoRun { .. }
                | AnalysisJob::PrepareDomestication { .. }
                | AnalysisJob::InspectMaterialOrigin { .. }
                | AnalysisJob::MaterialPcrDocuments { .. }
                | AnalysisJob::PrepareMaterialPcrOutputs { .. }
                | AnalysisJob::ReviewSequenceExport { .. }
                | AnalysisJob::PrepareSequenceExport { .. }
                | AnalysisJob::ExportScfTrace { .. }
                | AnalysisJob::InspectMaterialMismatchPcr { .. }
                | AnalysisJob::InspectMaterialRepair { .. }
                | AnalysisJob::InstallMaterialDerivation { .. }
                | AnalysisJob::InspectMaterialOverlapExtension { .. }
                | AnalysisJob::InspectMaterialTerminalProcessing { .. }
                | AnalysisJob::InspectMaterialAnnealedOligos { .. }
                | AnalysisJob::FindSimilarSequences { .. }
                | AnalysisJob::InspectMaterialDigest { .. }
                | AnalysisJob::PrepareSavedRunReview { .. }
                | AnalysisJob::PreviewProtocolLayout { .. }
                | AnalysisJob::CaptureProtocolLayout { .. }
                | AnalysisJob::ExportProtocolLayout { .. }
                | AnalysisJob::PreviewProtocolQuantities { .. }
                | AnalysisJob::CaptureProtocolQuantities { .. }
                | AnalysisJob::OpenProtocolQuantities { .. }
                | AnalysisJob::ExportProtocolQuantities { .. }
                | AnalysisJob::ListProtocolQuantities { .. }
        ) {
            MAX_WORKSPACE_MESSAGE_BYTES
        } else {
            MAX_IMPORT_BYTES
        };
        encode_with_limit(b"CGQ\0", self, maximum)
    }
}
