//! Request ownership and outbound validation.
use super::{
    envelope::{MAX_WORKSPACE_MESSAGE_BYTES, encode_with_limit},
    schema::AnalysisJob,
};
use crate::formats::sequence::MAX_IMPORT_BYTES;
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
        if self.revision.len() > 128 {
            return Err(Error::TooLarge);
        }
        if let AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            ..
        } = &self.job
        {
            crate::alignment::export::preflight(record, reference_name, query_name)?;
        }
        if let AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } = &self.job
        {
            crate::alignment::search::preflight(report, pattern, options)?;
        }
        if let AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } = &self.job
        {
            crate::alignment::context::preflight(report, reference, query, *translate)?;
        }
        match &self.job {
            AnalysisJob::PreviewAnnotationMerge { molecule, request }
            | AnalysisJob::ApplyAnnotationMerge { molecule, request } => {
                if crate::sequence::annotations::merge::revision(molecule, request)?
                    != self.revision
                {
                    return Err(Error::Format(
                        "annotation merge source/settings context changed".into(),
                    ));
                }
            }
            AnalysisJob::AssessDesignConsensus { workspace, request } => {
                if crate::alignment::design_consensus::revision(workspace, request)?
                    != self.revision
                {
                    return Err(Error::Format(
                        "Design/consensus source or settings context changed".into(),
                    ));
                }
            }
            AnalysisJob::Align {
                molecule,
                query,
                options,
            } => {
                let expected = crate::alignment::review::revision(
                    molecule,
                    query,
                    &crate::alignment::review::Options::Stored(*options),
                )?;
                if expected != self.revision {
                    return Err(Error::Format(
                        "Pairwise alignment sources/settings context changed".into(),
                    ));
                }
            }
            AnalysisJob::CompareSequences {
                reference,
                query,
                options,
            } => {
                if crate::alignment::review::revision(reference, query, options)? != self.revision {
                    return Err(Error::Format(
                        "Pairwise comparison sources/settings context changed".into(),
                    ));
                }
            }
            AnalysisJob::AssessScreeningStrategy { request } => {
                if crate::cloning::screening::revision(request)? != self.revision {
                    return Err(Error::Format(
                        "Screening source/settings context changed".into(),
                    ));
                }
            }
            AnalysisJob::AssessConstructCausalQc { workspace, request } => {
                if crate::construct_qc::causal::revision(workspace, request)? != self.revision {
                    return Err(Error::Format(
                        "Causal construct-QC workspace/request context changed".into(),
                    ));
                }
            }
            AnalysisJob::ListVectorClassProfiles => {
                if crate::construct_qc::vector_class::catalogue_revision()? != self.revision {
                    return Err(Error::Format(
                        "Vector-class catalogue context changed".into(),
                    ));
                }
            }
            AnalysisJob::AssessVectorClass { molecule, request } => {
                if crate::construct_qc::vector_class::revision(molecule, request)? != self.revision
                {
                    return Err(Error::Format(
                        "Vector-class source/declaration context changed".into(),
                    ));
                }
            }
            AnalysisJob::AssessIvtLinearization { molecule, request } => {
                if crate::cloning::restriction::ivt_linearization::revision(molecule, request)?
                    != self.revision
                {
                    return Err(Error::Format(
                        "IVT linearization source/declaration context changed".into(),
                    ));
                }
            }
            AnalysisJob::AnalyzeBaseEditing { molecule, request } => {
                if crate::sequence::crispr::base_editing::revision(molecule, request)?
                    != self.revision
                {
                    return Err(Error::Format(
                        "Base-editing source/declaration context changed".into(),
                    ));
                }
            }
            AnalysisJob::ListIvtProfiles => {
                if crate::rna_chemistry::ivt::catalogue_revision()
                    .map_err(|error| Error::Format(error.to_string()))?
                    != self.revision
                {
                    return Err(Error::Format(
                        "IVT profile catalogue context changed".into(),
                    ));
                }
            }
            AnalysisJob::AssessIvtTemplate { molecule, request } => {
                if crate::rna_chemistry::ivt::revision(molecule, request)
                    .map_err(|error| Error::Format(error.to_string()))?
                    != self.revision
                {
                    return Err(Error::Format(
                        "IVT source/declaration context changed".into(),
                    ));
                }
            }
            AnalysisJob::AnalyzeConstructQc { molecule, options } => {
                if crate::construct_qc::revision(molecule, options)? != self.revision {
                    return Err(Error::Format("Construct QC source/options differ".into()));
                }
            }
            AnalysisJob::InspectAssemblyStandard { request } => {
                if crate::cloning::assembly::standards::revision(request)? != self.revision {
                    return Err(Error::Format("Assembly-standard request differs".into()));
                }
            }
            AnalysisJob::SearchOligoCollection { request } => {
                if crate::cloning::pcr::collection_binding::revision(request)? != self.revision {
                    return Err(Error::Format("Oligo collection request differs".into()));
                }
            }
            AnalysisJob::AssessRestrictionConditions { selection, options } => {
                let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
                if crate::cloning::restriction::conditions::revision(&catalog, selection, options)?
                    != self.revision
                {
                    return Err(Error::Format("Restriction reaction context differs".into()));
                }
            }
            AnalysisJob::InversePcr {
                molecule,
                target,
                origin,
                options,
            } => {
                crate::cloning::pcr::mutagenesis::inverse_pcr_revision(
                    molecule,
                    target,
                    origin.as_ref(),
                    options,
                )?;
            }

            AnalysisJob::PreviewProtocolLayout { workspace, request }
            | AnalysisJob::CaptureProtocolLayout { workspace, request } => {
                crate::output::layout::preview_revision(workspace, request)?;
            }
            AnalysisJob::PreviewProtocolQuantities { workspace, request }
            | AnalysisJob::CaptureProtocolQuantities { workspace, request } => {
                crate::output::quantities::preview_revision(workspace, request)?;
            }
            AnalysisJob::ListTopoGraphs { workspace } => {
                crate::cloning::topo::graph_history::list_revision(workspace)?;
            }
            AnalysisJob::CaptureTopoGraph { workspace, request } => {
                crate::cloning::topo::graph_history::preflight(workspace, request)?;
            }
            AnalysisJob::InspectTopoGraph {
                workspace,
                selection,
            } => {
                crate::cloning::topo::graph_history::inspect_revision(workspace, selection)?;
            }
            AnalysisJob::ListProtocolQuantities { workspace } => {
                crate::output::quantities::list_revision(workspace)?;
            }
            AnalysisJob::OpenProtocolQuantities {
                workspace,
                quantity_hash,
            } => {
                crate::output::quantities::open_revision(workspace, quantity_hash)?;
            }
            AnalysisJob::ExportProtocolQuantities {
                workspace,
                quantity_hash,
                format,
            } => {
                crate::output::quantities::export_revision(workspace, quantity_hash, *format)?;
            }

            AnalysisJob::ExportProtocolLayout {
                workspace,
                layout_hash,
                format,
            } => {
                crate::output::layout::export_revision(workspace, layout_hash, *format)?;
            }

            AnalysisJob::InspectMaterialOrigin {
                workspace,
                event,
                selection,
            } => {
                crate::material::run::material_inspection_revision(workspace, event, selection)?;
            }
            AnalysisJob::MaterialPcrDocuments {
                molecule,
                origin,
                forward,
                reverse,
                template_strand,
            } => {
                crate::cloning::pcr::material::revision(
                    molecule,
                    origin,
                    forward,
                    reverse,
                    *template_strand,
                )?;
            }
            AnalysisJob::PrepareMaterialPcrOutputs { origin, request } => {
                crate::cloning::pcr::material::outputs_revision(origin, request)?;
            }
            AnalysisJob::ReviewSequenceExport { molecule, format }
            | AnalysisJob::PrepareSequenceExport {
                molecule, format, ..
            } => {
                crate::formats::sequence::export::revision(molecule, format)?;
            }
            AnalysisJob::ExportScfTrace { molecule, options } => {
                crate::formats::scf::document_revision(molecule, options)?;
            }
            AnalysisJob::InspectMaterialMismatchPcr {
                source,
                origin,
                forward,
                reverse,
                options,
            } => {
                crate::material::derivation::pcr::mismatch_pcr_preview_revision(
                    source, origin, forward, reverse, options,
                )?;
            }
            AnalysisJob::InspectMaterialRepair {
                source,
                origin,
                method,
            } => {
                crate::material::derivation::end_repair::repair_preview_revision(
                    source, origin, *method,
                )?;
            }
            AnalysisJob::InspectMaterialOverlapExtension {
                workspace,
                source,
                origin,
                preparation,
            } => {
                crate::material::derivation::overlap_extension::overlap_extension_preview_revision(
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
                crate::material::derivation::terminal::terminal_processing_preview_revision(
                    workspace, source, origin, action,
                )?;
            }
            AnalysisJob::InspectMaterialAnnealedOligos {
                workspace,
                top,
                bottom,
                options,
            } => {
                crate::material::derivation::annealing::annealed_oligos_preview_revision(
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
                crate::protein::antibody_alignment::request_revision(
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
                crate::protein::structure::revision(bytes, *format, selection)?;
            }
            AnalysisJob::FindSimilarSequences {
                query,
                subjects,
                options,
            } => {
                crate::sequence::similarity::preflight(query, subjects, options)?;
            }
            AnalysisJob::InstallMaterialDerivation { workspace, request } => {
                crate::material::derivation::preflight(workspace, request)?;
            }
            AnalysisJob::PrepareSavedRunReview { workspace, run_id } => {
                crate::material::review::request_revision(workspace, run_id)?;
            }
            AnalysisJob::PrepareImageEvidence {
                name,
                association,
                original,
            } => {
                if original.is_empty()
                    || original.len() > crate::documents::image_evidence::MAX_SOURCE_BYTES
                {
                    return Err(Error::TooLarge);
                }
                let hash = blake3::hash(original).to_hex().to_string();
                if association.preparation_revision(name, &hash)? != self.revision {
                    return Err(Error::Format("Image preparation context differs".into()));
                }
            }
            AnalysisJob::ReopenImageEvidence { prepared } => {
                if prepared.summary()?.reopen_revision()? != self.revision {
                    return Err(Error::Format("Image reopening context differs".into()));
                }
            }
            AnalysisJob::RenderMap { request } => request.validate()?,
            AnalysisJob::RenderGelImage { request } => {
                if request.revision()? != self.revision {
                    return Err(Error::Format("Gel image request context differs".into()));
                }
            }
            AnalysisJob::MaterializeMaterialDigestFragment { source, request } => {
                request.revision(source)?;
            }
            AnalysisJob::InspectMaterialDigest {
                source_role,
                source,
                origin,
                selection,
                options,
            } => {
                crate::material::derivation::digest::digest_inspection_revision(
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
                crate::cloning::topo::run_preparation::preflight(inputs, template.as_ref(), task)?;
            }
            AnalysisJob::PreviewAnnealingRecipe {
                inputs,
                target,
                recipe,
                policy,
            } => {
                if inputs.len() > crate::cloning::assembly::recipe::MAX_ROWS {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|input| (input.id.as_str(), &input.molecule))
                    .collect::<Vec<_>>();
                crate::cloning::assembly::recipe::annealing_source_preflight(
                    &sources, target, recipe, policy,
                )?;
            }
            AnalysisJob::PreviewHomologyRecipe {
                inputs,
                target,
                recipe,
                policy,
            } => {
                if inputs.len() > crate::cloning::assembly::recipe::MAX_ROWS {
                    return Err(Error::TooLarge);
                }
                let sources = inputs
                    .iter()
                    .map(|input| (input.id.as_str(), &input.molecule))
                    .collect::<Vec<_>>();
                crate::cloning::assembly::recipe::homology_source_preflight(
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
                crate::cloning::assembly::edited_primers::source_preflight(
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
                crate::cloning::overlap::extension::primers::preflight(
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
                crate::cloning::overlap::extension::digest::source_preflight(
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
                crate::cloning::pcr::mismatch_digest::source_preflight(
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
                crate::cloning::assembly::intent::preflight(&sources, request, workspace.as_ref())?;
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
                crate::cloning::assembly::junction_edit::intent_preflight(
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
                crate::cloning::assembly::edit_intent::source_preflight(
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
                crate::cloning::assembly::substitution_intent::source_preflight(
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
                let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
                crate::cloning::restriction::run_digest::source_preflight(
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
                crate::cloning::pcr::digest::preflight_from_material(
                    molecule, origin, request, workspace,
                )?;
            }
            AnalysisJob::InspectRunMaterial {
                workspace,
                event,
                selection,
            } => {
                crate::material::run::inspection_revision(workspace, event, selection)?;
            }
            AnalysisJob::InspectPcrDigest { molecule, request } => {
                crate::cloning::pcr::digest::source_preflight(molecule, request)?;
            }
            AnalysisJob::PreviewTailingRecipe {
                inputs,
                target,
                recipe,
            } => {
                crate::cloning::ligation::tailing_recipe::preflight(inputs, target, recipe)?;
            }
            AnalysisJob::InspectCombinatorialMultiSite {
                workspace,
                design,
                construct,
            } => {
                crate::planning::combinatorial::inspect_multisite_revision(
                    design, construct, workspace,
                )?;
            }
            AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
                crate::cloning::recombination::multisite::preflight(inputs, request)?;
            }
            AnalysisJob::ProjectMultiSiteRecombination {
                inputs,
                request,
                selection,
            } => {
                crate::cloning::recombination::multisite::product_revision(
                    inputs, request, selection,
                )?;
            }
            AnalysisJob::CheckMultiSiteRecombination {
                inputs,
                request,
                selection,
                target,
            } => {
                crate::cloning::recombination::multisite::check_revision(
                    inputs, request, selection, target,
                )?;
            }
            AnalysisJob::ProcessTopoActions { inputs, request } => {
                crate::cloning::topo::actions::preflight(inputs, request)?
            }

            AnalysisJob::ListLigationProfiles => {
                crate::cloning::ligation::profiles::catalog_revision()?;
            }
            AnalysisJob::ScorePublishedLigation {
                profile,
                ends,
                physical_competitors,
            } => {
                crate::cloning::ligation::profiles::score_revision_with_competitors(
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
                crate::cloning::ligation::fidelity::competition_revision(
                    dataset,
                    ends,
                    physical_competitors,
                )?;
            }
            AnalysisJob::RankGoldenGateFusions { inputs, request } => {
                crate::cloning::golden_gate::fusion_rank::preflight(
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
                crate::cloning::golden_gate::domestication::preflight(source, request)?;
            }
            AnalysisJob::PrepareDomestication {
                input,
                request,
                workspace,
            } => {
                crate::cloning::golden_gate::domestication::preparation::preflight(
                    input,
                    request,
                    workspace.as_ref(),
                )?;
            }
            AnalysisJob::SearchGoldenGateFusions { inputs, request }
            | AnalysisJob::MaterializeGoldenGateFusions {
                inputs, request, ..
            } => {
                crate::cloning::golden_gate::fusion_design::preflight(
                    &inputs
                        .iter()
                        .map(|input| (input.id.as_str(), &input.molecule))
                        .collect::<Vec<_>>(),
                    request,
                )?;
            }
            AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
                crate::cloning::ligation::fidelity::revision(dataset, ends)?;
            }
            AnalysisJob::BeginCombinatorialBatch { workspace, matrix } => {
                crate::planning::combinatorial::batch::begin_revision(workspace, matrix)?;
            }
            AnalysisJob::CaptureCombinatorialBatch {
                workspace,
                current_context,
            } => {
                crate::planning::combinatorial::batch::capture_revision(
                    workspace,
                    current_context,
                )?;
            }
            AnalysisJob::BeginCombinatorialRelease { workspace, request } => {
                crate::planning::combinatorial::release::preflight(workspace, request)?;
            }
            AnalysisJob::ExportCombinatorialOutputs {
                workspace,
                manifest_hash,
                options,
            } => {
                crate::planning::combinatorial::export::preflight(
                    workspace,
                    manifest_hash,
                    options,
                )?;
            }
            _ => {}
        }
        let maximum = if matches!(
            self.job,
            AnalysisJob::PreviewAnnotationMerge { .. }
                | AnalysisJob::ApplyAnnotationMerge { .. }
                | AnalysisJob::AssessConstructCausalQc { .. }
                | AnalysisJob::SaveWorkspace { .. }
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
                | AnalysisJob::InspectMaterialMismatchPcr { .. }
                | AnalysisJob::InspectMaterialRepair { .. }
                | AnalysisJob::InstallMaterialDerivation { .. }
                | AnalysisJob::InspectMaterialOverlapExtension { .. }
                | AnalysisJob::InspectMaterialTerminalProcessing { .. }
                | AnalysisJob::InspectMaterialAnnealedOligos { .. }
                | AnalysisJob::FindSimilarSequences { .. }
                | AnalysisJob::InspectMaterialDigest { .. }
                | AnalysisJob::PredictPinnedGel { .. }
                | AnalysisJob::MaterializeMaterialDigestFragment { .. }
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
        let encoded = encode_with_limit(b"CGQ\0", self, maximum);
        match encoded {
            Err(Error::TooLarge)
                if matches!(
                    self.job,
                    AnalysisJob::PrepareTraceDisplay { .. } | AnalysisJob::ExportScfTrace { .. }
                ) =>
            {
                Err(Error::Unsupported(
                    "trace document and worker framing exceed the ordinary 32 MiB request limit"
                        .into(),
                ))
            }
            other => other,
        }
    }
}
