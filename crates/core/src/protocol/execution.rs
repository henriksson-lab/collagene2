//! Stateless worker execution against the current typed schema.
use super::{
    bindings::{alignment_revision, overlap_profile_revision},
    envelope::{
        MAX_WORKSPACE_MESSAGE_BYTES, WORKSPACE_REQUEST_SELECTORS, decode_with_limit,
        encode_with_limit, worker_limit,
    },
    request::Request,
    response::Response,
    schema::{AnalysisJob, AnalysisResult, AnnotationExportReport, PrimerEvidence},
    validation::validate_verification_report,
};
use crate::formats::MAX_IMPORT_BYTES;
use collagene_model::Error;
pub fn process_request(bytes: &[u8]) -> Result<Vec<u8>, Error> {
    let maximum = worker_limit(b"CGQ\0", bytes, WORKSPACE_REQUEST_SELECTORS)?;

    let r: Request = decode_with_limit(b"CGQ\0", bytes, maximum)?;

    match &r.job {
        AnalysisJob::SelectGoldenGateReaction { report, .. } => {
            validate_verification_report(&report.verification)?
        }
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, .. } => {
            validate_verification_report(&report.preview.verification)?
        }
        AnalysisJob::LoadRunEvidence { run_id, .. } if run_id.len() != 64 => {
            return Err(Error::Format("Invalid run identity".into()));
        }
        _ => {}
    }
    if let AnalysisJob::PredictGel { lanes, .. } = &r.job {
        if lanes.iter().any(|lane| {
            matches!(
                lane.preparation,
                crate::gel::Preparation::CatalogDigest { .. }
            )
        }) {
            return Err(Error::Format(
                "Catalog-pinned lanes require the pinned-gel worker operation.".into(),
            ));
        }
    }
    let mol = match &r.job {
        AnalysisJob::InspectMaterialOrigin {
            workspace,
            event,
            selection,
        } => {
            if crate::run_material::material_inspection_revision(workspace, event, selection)?
                != r.revision
            {
                return Err(Error::Format("Material source/selection changed".into()));
            }
            None
        }
        AnalysisJob::MaterialPcrDocuments {
            molecule,
            origin,
            forward,
            reverse,
            template_strand,
        } => {
            if crate::material_pcr::revision(molecule, origin, forward, reverse, *template_strand)?
                != r.revision
            {
                return Err(Error::Format("Material PCR source/primers changed".into()));
            }
            None
        }
        AnalysisJob::PrepareMaterialPcrOutputs { origin, request } => {
            if crate::material_pcr::outputs_revision(origin, request)? != r.revision {
                return Err(Error::Format("Material PCR output choices changed".into()));
            }
            None
        }
        AnalysisJob::ReviewSequenceExport { molecule, format } => {
            if crate::sequence_export::revision(molecule, format)? != r.revision {
                return Err(Error::Format(
                    "Sequence export source/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareSequenceExport {
            molecule,
            format,
            expected_context_revision,
            ..
        } => {
            if expected_context_revision != &r.revision
                || crate::sequence_export::revision(molecule, format)? != r.revision
            {
                return Err(Error::Format(
                    "Reviewed sequence export context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ExportScfTrace { trace, options } => {
            if crate::scf_export::revision(trace, options)? != r.revision {
                return Err(Error::Format(
                    "Trace export source/view/options changed".into(),
                ));
            }
            None
        }

        AnalysisJob::PrepareSavedRunReview { workspace, run_id } => {
            if crate::run_review::request_revision(workspace, run_id)? != r.revision {
                return Err(Error::Format(
                    "Saved-run review workspace/run changed".into(),
                ));
            }
            None
        }
        AnalysisJob::RenderMapSvg { request } => {
            request.validate()?;
            if r.revision.len() != 64
                || !r
                    .revision
                    .bytes()
                    .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
            {
                return Err(Error::Format(
                    "Map export capture context is invalid".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialMismatchPcr {
            source,
            origin,
            forward,
            reverse,
            options,
        } => {
            if crate::material_derivation::mismatch_pcr_preview_revision(
                source, origin, forward, reverse, options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Material mismatch-PCR source/origin/primers/policy changed".into(),
                ));
            }
            None
        }

        AnalysisJob::InversePcr {
            molecule,
            target,
            origin,
            options,
        } => {
            if crate::mutagenesis::inverse_pcr_revision(molecule, target, origin.as_ref(), options)?
                != r.revision
            {
                return Err(Error::Format(
                    "Inverse PCR source, origin, independent target or edit settings changed"
                        .into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialRepair {
            source,
            origin,
            method,
        } => {
            if crate::material_derivation::repair_preview_revision(source, origin, *method)?
                != r.revision
            {
                return Err(Error::Format(
                    "Material repair source/origin/policy changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialOverlapExtension {
            workspace,
            source,
            origin,
            preparation,
        } => {
            if crate::material_derivation::overlap_extension_preview_revision(
                workspace,
                source,
                origin,
                preparation,
            )? != r.revision
            {
                return Err(Error::Format(
                    "OE material source/origin/preparation changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialTerminalProcessing {
            workspace,
            source,
            origin,
            action,
        } => {
            if crate::material_derivation::terminal_processing_preview_revision(
                workspace, source, origin, action,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Terminal-processing material source/origin/action changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialAnnealedOligos {
            workspace,
            top,
            bottom,
            options,
        } => {
            if crate::material_derivation::annealed_oligos_preview_revision(
                workspace, top, bottom, options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Annealed-material source events/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AnalyzeAlignmentAntibody {
            alignment,
            source,
            options,
            limits,
            selection,
        } => {
            if crate::alignment_antibody::request_revision(
                alignment, source, options, limits, selection,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Alignment antibody source/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AnalyzeAntibody {
            molecule,
            options,
            limits,
            selection,
        } => {
            if crate::antibody_numbering::domain_tracks::request_revision(
                molecule, options, limits, selection,
            )
            .map_err(|e| Error::Format(e.to_string()))?
                != r.revision
            {
                return Err(Error::Format(
                    "Antibody source/options/limits changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ScreenProteinLiabilities { molecule, options } => {
            if crate::protein_liabilities::request_revision(molecule, options)
                .map_err(|e| Error::Format(e.to_string()))?
                != r.revision
            {
                return Err(Error::Format(
                    "Protein liability source/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AnalyzeRnaChemistry {
            molecule,
            declaration,
        } => {
            if crate::rna_chemistry::request_revision(molecule, declaration)
                .map_err(|e| Error::Format(e.to_string()))?
                != r.revision
            {
                return Err(Error::Format("RNA source/declaration changed".into()));
            }
            None
        }
        AnalysisJob::PrepareStructure {
            bytes,
            format,
            selection,
        } => {
            if crate::structure_view::revision(bytes, *format, selection)? != r.revision {
                return Err(Error::Format("Structure source/options changed".into()));
            }
            None
        }
        AnalysisJob::FindSimilarSequences {
            query,
            subjects,
            options,
        } => {
            if crate::sequence_similarity::context_revision(query, subjects, options)? != r.revision
            {
                return Err(Error::Format(
                    "Sequence similarity sources/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InstallMaterialDerivation { workspace, request } => {
            if crate::material_derivation::revision(workspace, request)? != r.revision {
                return Err(Error::Format(
                    "Material derivation workspace/source/recipe changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialDigest {
            source_role,
            source,
            origin,
            selection,
            options,
        } => {
            if crate::material_derivation::digest_inspection_revision(
                source_role,
                source,
                origin,
                selection,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Material digest source/origin/options changed".into(),
                ));
            }
            None
        }

        AnalysisJob::FinalizeCloningOutputs {
            workspace,
            plan,
            target,
            inputs,
            options,
        } => {
            if crate::cloning_outputs::finalize_revision(
                workspace,
                plan,
                target.clone(),
                inputs.clone(),
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Cloning output source/options context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SelectDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            expected_batch_revision,
            choice,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::digest_cut_selection::selection_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
                expected_batch_revision,
                choice,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Digest cut selection source/boundaries changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ExportPlanOligos {
            plan,
            target,
            options,
            delimiter,
        } => {
            if crate::plan_oligos::export_revision(plan, target, options, *delimiter)? != r.revision
            {
                return Err(Error::Format("Plan oligo export context changed".into()));
            }
            None
        }
        AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin,
            request,
            workspace,
        } => {
            if crate::mismatch_pcr_digest_inspection::inspection_revision(
                molecule,
                origin.as_ref(),
                request,
                workspace.as_ref(),
            )? != r.revision
            {
                return Err(Error::Format(
                    "Mismatch PCR digest original source/origin/policy/catalog context changed"
                        .into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewAnnealingRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => {
            if crate::assembly_recipe::annealing_preview_revision(inputs, target, recipe, policy)?
                != r.revision
            {
                return Err(Error::Format("3-prime resection/annealing recipe source, target, preparation or explicit physical policy changed".into()));
            }
            None
        }
        AnalysisJob::PreviewHomologyRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => {
            if crate::assembly_recipe::homology_preview_revision(inputs, target, recipe, policy)?
                != r.revision
            {
                return Err(Error::Format("Homology recipe source, target, preparation or explicit physical policy changed".into()));
            }
            None
        }
        AnalysisJob::DesignEditedOverlapPrimers {
            inputs,
            request,
            workspace,
        } => {
            if crate::assembly_edited_primer_design::design_revision(
                inputs,
                request,
                workspace.as_ref(),
            )? != r.revision
            {
                return Err(Error::Format(
                    "Edited overlap primer original sources, intent, origins or policy changed"
                        .into(),
                ));
            }
            None
        }
        AnalysisJob::InspectOverlapExtensionDigest {
            source,
            origin,
            request,
            workspace,
        } => {
            if crate::overlap_extension_digest_inspection::inspection_revision(
                source,
                origin.as_ref(),
                request,
                workspace.as_ref(),
            )? != r.revision
            {
                return Err(Error::Format(
                    "Overlap-extension digest source, origin, oligos or catalog/options changed"
                        .into(),
                ));
            }
            None
        }
        AnalysisJob::DesignOverlapExtensionPrimers {
            source,
            request,
            workspace,
        } => {
            if crate::overlap_extension_primer_design::design_revision(
                source,
                request,
                workspace.as_ref(),
            )? != r.revision
            {
                return Err(Error::Format(
                    "Overlap-extension primer source, independent edits or options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildContributionEditIntent {
            inputs,
            request,
            workspace,
        } => {
            if crate::assembly_edit_intent::intent_revision(inputs, request, workspace.as_ref())?
                != r.revision
            {
                return Err(Error::Format(
                    "Contribution edit original source/declaration changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildSubstitutedAssemblyIntent {
            inputs,
            request,
            workspace,
        } => {
            if crate::assembly_substitution_intent::intent_revision(
                inputs,
                request,
                workspace.as_ref(),
            )? != r.revision
            {
                return Err(Error::Format(
                    "Substituted assembly original source/declaration changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildAssemblyIntent {
            inputs,
            request,
            workspace,
        } => {
            if crate::assembly_intent::intent_revision(inputs, request, workspace.as_ref())?
                != r.revision
            {
                return Err(Error::Format(
                    "assembly intent source/region/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildEditedAssemblyIntent {
            inputs,
            base,
            edits,
            workspace,
        } => {
            if crate::assembly_junction_edit::intent_revision(
                inputs,
                base,
                edits,
                workspace.as_ref(),
            )? != r.revision
            {
                return Err(Error::Format(
                    "Edited assembly source/intent changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignEditedAssemblyPrimers { inputs, request } => {
            if crate::assembly_junction_edit::design_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Edited assembly source/primers changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request } => {
            if crate::selective_junction_edit::design_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Selective assembly source/primers changed".into(),
                ));
            }
            None
        }

        AnalysisJob::PreviewProtocolQuantities { workspace, request }
        | AnalysisJob::CaptureProtocolQuantities { workspace, request } => {
            if crate::protocol_quantities::preview_revision(workspace, request)? != r.revision {
                return Err(Error::Format(
                    "Reaction quantity source/recipe context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ListTopoGraphs { workspace } => {
            if crate::topo_graph_history::list_revision(workspace)? != r.revision {
                return Err(Error::Format("TOPO graph listing workspace changed".into()));
            }
            None
        }
        AnalysisJob::CaptureTopoGraph { workspace, request } => {
            if crate::topo_graph_history::revision(workspace, request)? != r.revision {
                return Err(Error::Format(
                    "TOPO graph capture source/request changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectTopoGraph {
            workspace,
            selection,
        } => {
            if crate::topo_graph_history::inspect_revision(workspace, selection)? != r.revision {
                return Err(Error::Format(
                    "TOPO graph selection/workspace changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ListProtocolQuantities { workspace } => {
            if crate::protocol_quantities::list_revision(workspace)? != r.revision {
                return Err(Error::Format(
                    "Reaction quantity catalog/workspace changed".into(),
                ));
            }
            None
        }
        AnalysisJob::OpenProtocolQuantities {
            workspace,
            quantity_hash,
        } => {
            if crate::protocol_quantities::open_revision(workspace, quantity_hash)? != r.revision {
                return Err(Error::Format(
                    "Reaction quantity attachment/workspace changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ExportProtocolQuantities {
            workspace,
            quantity_hash,
            format,
        } => {
            if crate::protocol_quantities::export_revision(workspace, quantity_hash, *format)?
                != r.revision
            {
                return Err(Error::Format(
                    "Reaction quantity export source/format changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewProtocolLayout { workspace, request }
        | AnalysisJob::CaptureProtocolLayout { workspace, request } => {
            if crate::protocol_layout::preview_revision(workspace, request)? != r.revision {
                return Err(Error::Format(
                    "Protocol layout source or assignment context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ExportProtocolLayout {
            workspace,
            layout_hash,
            format,
        } => {
            if crate::protocol_layout::export_revision(workspace, layout_hash, *format)?
                != r.revision
            {
                return Err(Error::Format(
                    "Protocol layout export source or format changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectGoldenGateReaction {
            inputs,
            target,
            request,
        } => {
            if crate::golden_gate_reaction::inspect_revision(inputs, target, request)? != r.revision
            {
                return Err(Error::Format("Golden Gate reaction context changed".into()));
            }
            None
        }
        AnalysisJob::PreparePcrOutputs { request } => {
            if crate::pcr_outputs::revision(request)? != r.revision {
                return Err(Error::Format("PCR output selection changed".into()));
            }
            None
        }
        AnalysisJob::CopyReferenceFragment { request } => {
            if crate::reference_transfer::copy_revision(request)? != r.revision {
                return Err(Error::Format(
                    "Reference copy source or selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PasteReferenceFragment { request } => {
            if crate::reference_transfer::paste_revision(request)? != r.revision {
                return Err(Error::Format(
                    "Reference paste sources or selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } => {
            if crate::pairwise_context::revision(report, reference, query, *translate)?
                != r.revision
            {
                return Err(Error::Format("Pairwise annotation context changed.".into()));
            }
            None
        }
        AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } => {
            if crate::pairwise_search::revision(report, pattern, options)? != r.revision {
                return Err(Error::Format("Pairwise search context changed.".into()));
            }
            None
        }
        AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            fresh,
            format,
        } => {
            if crate::pairwise_export::revision(
                record,
                reference_name,
                query_name,
                *fresh,
                *format,
            )? != r.revision
            {
                return Err(Error::Format("Pairwise export context changed.".into()));
            }
            None
        }
        AnalysisJob::PredictPinnedGel {
            molecule,
            lanes,
            calibration,
        } => {
            if lanes.is_empty() || lanes.len() > 32 {
                return Err(Error::TooLarge);
            }
            calibration.validate()?;
            molecule.validate()?;
            let catalog = crate::catalog_digest::builtin_catalog()?;
            let mut bases = 0usize;
            for lane in lanes {
                lane.molecule.validate()?;
                bases = bases
                    .checked_add(lane.molecule.sequence.len())
                    .filter(|n| *n <= MAX_IMPORT_BYTES)
                    .ok_or(Error::TooLarge)?;
                if let crate::gel::Preparation::CatalogDigest { recipe } = &lane.preparation {
                    collagene_model::bounds::model_preflight(&[&lane.molecule])?;
                    recipe.validate()?;
                    recipe.selection.validate_for(&catalog)?;
                    if lane.molecule.revision() != recipe.source.revision {
                        return Err(Error::Format("Pinned gel source revision changed.".into()));
                    }
                }
            }
            Some(molecule)
        }
        AnalysisJob::TranslateReference { request } => {
            if crate::translation_reference::inspect_revision(request)? != r.revision {
                return Err(Error::Format(
                    "Translation source or selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareTranslatedReference { request, name } => {
            if crate::translation_reference::output_revision(request, name)? != r.revision {
                return Err(Error::Format(
                    "Translation source, selection or output name changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CompareSequences {
            reference,
            query,
            options,
        } => {
            if crate::pairwise_review::revision(reference, query, options)? != r.revision {
                return Err(Error::Format(
                    "Comparison sources or options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareMismatchPcrOutputs { request } => {
            if crate::mismatch_pcr_outputs::revision(request)? != r.revision {
                return Err(Error::Format(
                    "Mismatch-PCR output selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SelectGoldenGateReaction { report, selection } => {
            if crate::golden_gate_reaction::selection_revision(report, selection)? != r.revision {
                return Err(Error::Format(
                    "Golden Gate reaction selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectTopoRun {
            inputs,
            template,
            task,
        } => {
            if crate::topo_run_preparation::revision(inputs, template.as_ref(), task)? != r.revision
            {
                return Err(Error::Format(
                    "Saved TOPO source, origin or inspection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProcessTopoActions { inputs, request } => {
            if crate::topo_actions::revision(inputs, request)? != r.revision {
                return Err(Error::Format("TOPO actions product context changed".into()));
            }
            None
        }

        AnalysisJob::ListLigationProfiles => {
            if crate::ligation_profiles::catalog_revision()? != r.revision {
                return Err(Error::Format(
                    "Published ligation catalog context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ScorePublishedLigation {
            profile,
            ends,
            physical_competitors,
        } => {
            if crate::ligation_profiles::score_revision_with_competitors(
                profile,
                ends,
                physical_competitors,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Published ligation profile/end context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ScoreLigationCompetition {
            dataset,
            ends,
            physical_competitors,
        } => {
            if crate::ligation_fidelity::competition_revision(dataset, ends, physical_competitors)?
                != r.revision
            {
                return Err(Error::Format("Ligation competition context changed".into()));
            }
            None
        }
        AnalysisJob::RankGoldenGateFusions { inputs, request } => {
            if crate::golden_gate_fusion_rank::inspect_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Fusion ranking source/profile context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectDomestication { source, request } => {
            if crate::domestication::inspect_revision(source, request)? != r.revision {
                return Err(Error::Format(
                    "Domestication inspection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProposeDomestication { source, request } => {
            if crate::domestication::propose_revision(source, request)? != r.revision {
                return Err(Error::Format(
                    "Domestication proposal context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckDomestication {
            source,
            request,
            selection,
        } => {
            if crate::domestication::check_revision(source, request, selection)? != r.revision {
                return Err(Error::Format(
                    "Domestication substitution context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareDomestication {
            input,
            request,
            workspace,
        } => {
            if crate::domestication_preparation::revision(input, request, workspace.as_ref())?
                != r.revision
            {
                return Err(Error::Format(
                    "Domestication preparation context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SearchGoldenGateFusions { inputs, request } => {
            if crate::golden_gate_fusion_design::inspect_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Golden Gate fusion search context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::MaterializeGoldenGateFusions {
            inputs,
            request,
            selection,
        } => {
            if crate::golden_gate_fusion_design::materialize_revision(inputs, request, selection)?
                != r.revision
            {
                return Err(Error::Format("Golden Gate fusion selection changed".into()));
            }
            None
        }
        AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
            if crate::ligation_fidelity::revision(dataset, ends)? != r.revision {
                return Err(Error::Format(
                    "Ligation fidelity data or ends changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request } => {
            if crate::golden_gate_junction_edit::design_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Golden Gate junction source/primers changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request,
        } => {
            if crate::prepared_golden_gate::inspect_revision(inputs, target, request)? != r.revision
            {
                return Err(Error::Format(
                    "Prepared Golden Gate pool context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, selection } => {
            if crate::prepared_golden_gate::selection_revision(report, selection)? != r.revision {
                return Err(Error::Format(
                    "Prepared Golden Gate selection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewReadReview { workspace, request } => {
            if crate::read_review::preview_revision(workspace, request)? != r.revision {
                return Err(Error::Format("Read review preview context changed".into()));
            }
            None
        }
        AnalysisJob::CommitReadReview { workspace, request } => {
            if crate::read_review::commit_revision(workspace, request)? != r.revision {
                return Err(Error::Format("Read review commit context changed".into()));
            }
            None
        }
        AnalysisJob::CaptureReadCheck {
            workspace,
            reference,
            reads,
            options,
        } => {
            if crate::read_checks::capture_revision(workspace, reference, reads, options)?
                != r.revision
            {
                return Err(Error::Format("Read-check capture context changed".into()));
            }
            None
        }
        AnalysisJob::ReplayReadCheck {
            workspace,
            attachment_hash,
        } => {
            if crate::read_checks::replay_revision(workspace, attachment_hash)? != r.revision {
                return Err(Error::Format("Read-check replay context changed".into()));
            }
            None
        }
        AnalysisJob::BuildMutationIntent { molecule, request } => {
            if crate::mutation_intent::intent_revision(molecule, request)? != r.revision {
                return Err(Error::Format("Mutation intent source/edit changed".into()));
            }
            None
        }
        AnalysisJob::DesignGoldenGatePrimers {
            inputs,
            target,
            request,
        } => {
            if crate::golden_gate_primer_design::preview_revision(inputs, target, request)?
                != r.revision
            {
                return Err(Error::Format(
                    "Golden Gate primer sources/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignOverlapPrimers {
            inputs,
            target,
            request,
        } => {
            if crate::overlap_primer_design::design_revision(inputs, target, request)? != r.revision
            {
                return Err(Error::Format(
                    "overlap primer sources/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectAssemblyJunctions {
            plan,
            target,
            options,
        } => {
            if crate::assembly_junctions::inspection_revision(plan, target, options)? != r.revision
            {
                return Err(Error::Format(
                    "assembly junction plan/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRunMaterial {
            workspace,
            event,
            selection,
        } => {
            if crate::run_material::inspection_revision(workspace, event, selection)? != r.revision
            {
                return Err(Error::Format(
                    "Run material workspace/event/origin selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRunCatalogRestriction {
            source_role,
            molecule,
            origin,
            selection,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::run_catalog_digest::inspection_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                origin,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Saved-run catalog inspection origin/source/pin/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectMaterialPcrDigest {
            molecule,
            origin,
            request,
            workspace,
        } => {
            if crate::pcr_digest_inspection::inspection_revision_from_material(
                molecule, origin, request, workspace,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Saved-run PCR digest origin/source/primer/catalog context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectPcrDigest { molecule, request } => {
            if crate::pcr_digest_inspection::inspection_revision(molecule, request)? != r.revision {
                return Err(Error::Format(
                    "PCR digest original source, primers or catalog options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        } => {
            if crate::tailing_recipe::preview_revision(inputs, target, recipe)? != r.revision {
                return Err(Error::Format(
                    "tailing preparation source/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => {
            if crate::ligation_recipe::preview_revision(inputs, target, recipe, *policy)?
                != r.revision
            {
                return Err(Error::Format(
                    "ligation preparation source/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        } => {
            if crate::assembly_recipe::preview_revision(inputs, target, recipe)? != r.revision {
                return Err(Error::Format(
                    "overlap preparation source/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BeginCombinatorialBatch { .. }
        | AnalysisJob::AdvanceCombinatorialBatch { .. }
        | AnalysisJob::CaptureCombinatorialBatch { .. } => {
            return Err(Error::Unsupported(
                "combinatorial batches require a persistent private worker session".into(),
            ));
        }
        AnalysisJob::BeginCombinatorialRelease { .. }
        | AnalysisJob::AdvanceCombinatorialRelease { .. }
        | AnalysisJob::InspectCombinatorialRelease { .. }
        | AnalysisJob::FinishCombinatorialRelease { .. } => {
            return Err(Error::Unsupported(
                "combinatorial release requires a persistent private worker session".into(),
            ));
        }
        AnalysisJob::ExportCombinatorialOutputs {
            workspace,
            manifest_hash,
            options,
        } => {
            if crate::combinatorial_export::revision(workspace, manifest_hash, options)?
                != r.revision
            {
                return Err(Error::Format("Combinatorial export inputs changed.".into()));
            }
            None
        }
        AnalysisJob::PrepareVariants { .. }
        | AnalysisJob::VariantPage { .. }
        | AnalysisJob::CheckVariantReference { .. }
        | AnalysisJob::PreviewVariantAllele { .. }
        | AnalysisJob::ExportVariants { .. }
        | AnalysisJob::DisposeVariants { .. } => {
            return Err(Error::Unsupported(
                "variant operations require a persistent worker session".into(),
            ));
        }
        AnalysisJob::SimulateTopo { inputs, request } => {
            if crate::topo_worker::simulation_revision(inputs, request)? != r.revision {
                return Err(Error::Format("TOPO simulation context changed".into()));
            }
            None
        }
        AnalysisJob::InspectTopoOrientations {
            inputs,
            request,
            vector,
            insert,
            directional,
        } => {
            if crate::topo_worker::orientations_revision(
                inputs,
                request,
                *vector,
                *insert,
                *directional,
            )? != r.revision
            {
                return Err(Error::Format("TOPO orientation context changed".into()));
            }
            None
        }
        AnalysisJob::ProjectTopoProduct {
            inputs,
            request,
            selection,
        } => {
            if crate::topo_worker::projection_revision(inputs, request, selection)? != r.revision {
                return Err(Error::Format("TOPO projection context changed".into()));
            }
            None
        }
        AnalysisJob::CheckTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => {
            if crate::topo_worker::target_revision(inputs, request, selection, target, policy)?
                != r.revision
            {
                return Err(Error::Format("TOPO target/policy context changed".into()));
            }
            None
        }
        AnalysisJob::PrepareReadAlignments { .. }
        | AnalysisJob::ReadAlignmentPage { .. }
        | AnalysisJob::CheckReadAlignment { .. }
        | AnalysisJob::ExportReadAlignments { .. }
        | AnalysisJob::DisposeReadAlignments { .. } => {
            return Err(Error::Unsupported(
                "read alignment operations require a persistent worker session".into(),
            ));
        }
        AnalysisJob::InspectCombinatorialMultiSite {
            workspace,
            design,
            construct,
        } => {
            if crate::combinatorial_design::inspect_multisite_revision(
                design, construct, workspace,
            )? != r.revision
            {
                return Err(Error::Format(
                    "combinatorial multisite inspection source/declaration context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
            if crate::multisite_recombination::revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "multisite sources/profile/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection,
        } => {
            if crate::multisite_recombination::product_revision(inputs, request, selection)?
                != r.revision
            {
                return Err(Error::Format("multisite selected component changed".into()));
            }
            None
        }
        AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target,
        } => {
            if crate::multisite_recombination::check_revision(inputs, request, selection, target)?
                != r.revision
            {
                return Err(Error::Format(
                    "multisite independent target context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRecombination { inputs, request } => {
            if crate::recombination::batch_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "recombination source/profile/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckRecombinationTarget {
            inputs,
            request,
            target,
            policy,
        } => {
            if crate::recombination::check_revision(inputs, request, target, policy)? != r.revision
            {
                return Err(Error::Format("recombination target binding changed".into()));
            }
            None
        }
        AnalysisJob::SelectRecombinationComponent {
            inputs,
            request,
            selection,
            target,
            policy,
        } => {
            if crate::recombination::checked_product_revision(
                inputs, request, selection, target, policy,
            )? != r.revision
            {
                return Err(Error::Format(
                    "recombination selected component binding changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectCatalogRestriction {
            selection,
            source_role,
            molecule,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::catalog_digest::batch_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "catalog restriction source/pin/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::MaterializeCatalogDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::catalog_digest::product_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "catalog product source/pin/options/index changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InstallCatalogDigestFragment {
            selection,
            workspace,
            source,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::catalog_digest::installation_revision(
                &catalog,
                selection,
                workspace,
                source,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "catalog installation workspace/source/pin/options/index changed".into(),
                ));
            }
            None
        }

        AnalysisJob::InstallDigestFragment {
            workspace,
            source,
            enzymes,
            options,
            fragment_index,
        } => {
            if crate::restriction_history::installation_revision(
                workspace,
                source,
                enzymes,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "digest installation workspace/source/geometry/options/index changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        }
        | AnalysisJob::ScanRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => {
            if crate::restriction_inspection::batch_revision(
                source_role,
                molecule,
                enzymes,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "restriction inspection source/role/geometry/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::MaterializeDigestFragment {
            source_role,
            molecule,
            enzymes,
            options,
            fragment_index,
        } => {
            if crate::restriction_inspection::product_revision(
                source_role,
                molecule,
                enzymes,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "digest product source/role/geometry/options/index changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BeginGenomeFasta { .. }
        | AnalysisJob::PushGenomeFasta { .. }
        | AnalysisJob::FinishGenomeFasta { .. }
        | AnalysisJob::BeginGenomeNative { .. }
        | AnalysisJob::PushGenomeNativeBlock { .. }
        | AnalysisJob::FinishGenomeNative { .. }
        | AnalysisJob::AbortGenomeImport { .. }
        | AnalysisJob::RegisterGenome { .. }
        | AnalysisJob::GenomeIntervalGuides { .. }
        | AnalysisJob::StartGenomeSearch { .. }
        | AnalysisJob::AdvanceGenomeSearch { .. }
        | AnalysisJob::PageGenomeSearch { .. }
        | AnalysisJob::CancelGenomeJob { .. }
        | AnalysisJob::ReadGenomeInterval { .. } => None,
        AnalysisJob::DesignEditTemplate {
            source,
            target,
            options,
            guide,
        } => {
            if crate::edit_template::batch_revision(source, target, options, guide.as_ref())?
                != r.revision
            {
                return Err(Error::Format(
                    "edit-template source/target/options/guide changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckEditTemplate {
            source,
            target,
            options,
            guide,
            donor,
        } => {
            if crate::edit_template::check_batch_revision(
                source,
                target,
                options,
                guide.as_ref(),
                donor,
            )? != r.revision
            {
                return Err(Error::Format(
                    "edit-template candidate/source/target/options/guide changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignCrisprGuides { molecule, options } => {
            if crate::crispr::guide_revision(molecule, options)? != r.revision {
                return Err(Error::Format("CRISPR guide source/options changed".into()));
            }
            None
        }
        AnalysisJob::SearchCrisprOffTargets {
            query,
            references,
            options,
        } => {
            if crate::crispr::off_target_revision(query, references, options)? != r.revision {
                return Err(Error::Format(
                    "CRISPR off-target query/references/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AssessMethylation { molecule, options } => {
            if crate::methylation::batch_revision(molecule, options)? != r.revision {
                return Err(Error::Format(
                    "methylation source/options/profile changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AssembleContigs { reads, options } => {
            if crate::contig::batch_revision(reads, options)? != r.revision {
                return Err(Error::Format("contig batch source/options changed".into()));
            }
            None
        }
        AnalysisJob::RecordCloningRun { workspace, .. }
        | AnalysisJob::ReplayCloningRun { workspace, .. }
        | AnalysisJob::LoadRunEvidence { workspace, .. } => {
            if crate::run_lineage::workspace_fingerprint(workspace)? != r.revision {
                return Err(Error::Format("run request source workspace changed".into()));
            }
            None
        }
        AnalysisJob::CachedViewTracks { .. } => None,
        AnalysisJob::CheckOverlapProfile {
            molecules,
            reverse,
            overlaps,
            circular,
            options,
        } => {
            if overlap_profile_revision(molecules, reverse, overlaps, *circular, options)?
                != r.revision
            {
                return Err(Error::Format(
                    "overlap profile request source/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ImportAlignment { bytes, .. } | AnalysisJob::LoadAlignment { bytes } => {
            if blake3::hash(bytes).to_hex().as_str() != r.revision {
                return Err(Error::Format("alignment input digest mismatch".into()));
            }
            None
        }
        AnalysisJob::ExportAlignment { document, .. }
        | AnalysisJob::ExtractAlignmentRows { document, .. }
        | AnalysisJob::SaveAlignment { document } => {
            if alignment_revision(document)? != r.revision {
                return Err(Error::Format("alignment document revision mismatch".into()));
            }
            None
        }
        AnalysisJob::DotPlot {
            molecule,
            query,
            options,
        } => {
            if crate::dotplot::source_revisions(molecule, query, options)?.0 != r.revision {
                return Err(Error::Format(
                    "dot plot request source revision mismatch".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignCodons { molecule, options }
        | AnalysisJob::CheckCodonDesign {
            source: molecule,
            options,
            ..
        } => {
            // Bound source evidence and derive the required protein before hashing.
            // These jobs must not take the generic molecule-revision path first.
            if crate::codon_design::source_revision(molecule, options)? != r.revision {
                return Err(Error::Format(
                    "codon request source revision mismatch".into(),
                ));
            }
            None
        }
        AnalysisJob::TailPreview { molecule, .. }
        | AnalysisJob::ProteinProperties { molecule, .. }
        | AnalysisJob::Search { molecule, .. }
        | AnalysisJob::PrepareTracks { molecule }
        | AnalysisJob::ViewTracks { molecule, .. }
        | AnalysisJob::SearchDocument { molecule, .. }
        | AnalysisJob::Restriction { molecule, .. }
        | AnalysisJob::Pcr { molecule, .. }
        | AnalysisJob::PcrDocuments { molecule, .. }
        | AnalysisJob::AnalyzePrimer { molecule, .. }
        | AnalysisJob::DesignPrimers { molecule, .. }
        | AnalysisJob::GoldenGate { molecule, .. }
        | AnalysisJob::MismatchPcr { molecule, .. }
        | AnalysisJob::PredictGel { molecule, .. }
        | AnalysisJob::DiagnosticDigests { molecule, .. }
        | AnalysisJob::MismatchPcrPreview { molecule, .. }
        | AnalysisJob::DuplexDigest { molecule, .. }
        | AnalysisJob::Translate { molecule, .. }
        | AnalysisJob::TranslateCds { molecule, .. }
        | AnalysisJob::FindOrfs { molecule, .. }
        | AnalysisJob::Align { molecule, .. }
        | AnalysisJob::ValidateTraces { molecule, .. }
        | AnalysisJob::ScanFeatureLibrary { molecule, .. }
        | AnalysisJob::ApplyFeatureLibrary { molecule, .. }
        | AnalysisJob::ImportAnnotations { molecule, .. }
        | AnalysisJob::ApplyAnnotations { molecule, .. }
        | AnalysisJob::ExportAnnotations { molecule, .. }
        | AnalysisJob::RepairPreview { molecule, .. }
        | AnalysisJob::AnnealPreview { molecule, .. }
        | AnalysisJob::AnnealedDocument { molecule, .. } => Some(molecule),
        AnalysisJob::Verify { target, .. }
        | AnalysisJob::ExportProtocolPdf { target, .. }
        | AnalysisJob::ExportProtocolHtml { target, .. } => Some(target),
        AnalysisJob::SaveWorkspace { workspace } => {
            workspace.validate()?;
            let head = workspace.head(&workspace.active_document)?;
            if head.revision() != r.revision {
                return Err(Error::Format(
                    "workspace worker input revision mismatch".into(),
                ));
            }
            None
        }
        AnalysisJob::LoadWorkspace { bytes } => {
            if bytes.len() > crate::workspace::MAX_ARCHIVE_BYTES {
                return Err(Error::TooLarge);
            }
            if blake3::hash(bytes).to_hex().as_str() != r.revision {
                return Err(Error::Format(
                    "workspace archive worker digest mismatch".into(),
                ));
            }
            None
        }
    };
    if let Some(mol) = mol {
        mol.validate()?;
        if mol.revision() != r.revision {
            return Err(Error::Format("worker input revision mismatch".into()));
        }
    }
    let result = match r.job {
        AnalysisJob::InspectMaterialOrigin {
            workspace,
            event,
            selection,
        } => crate::run_material::inspect_material(&workspace, &event, &selection)
            .map(AnalysisResult::MaterialOriginInspection),
        AnalysisJob::MaterialPcrDocuments {
            molecule,
            origin,
            forward,
            reverse,
            template_strand,
        } => crate::material_pcr::simulate(&molecule, &origin, &forward, &reverse, template_strand)
            .map(AnalysisResult::MaterialPcrDocuments),
        AnalysisJob::PrepareMaterialPcrOutputs { origin, request } => {
            crate::material_pcr::prepare_outputs(&origin, &request)
                .map(AnalysisResult::MaterialPcrOutputsPrepared)
        }
        AnalysisJob::ReviewSequenceExport { molecule, format } => {
            crate::sequence_export::review(&molecule, &format)
                .map(AnalysisResult::SequenceExportAssessment)
        }
        AnalysisJob::PrepareSequenceExport {
            molecule,
            format,
            expected_context_revision,
            allow_omissions,
        } => crate::sequence_export::prepare(
            &molecule,
            &format,
            &expected_context_revision,
            allow_omissions,
        )
        .map(AnalysisResult::SequenceExportPrepared),
        AnalysisJob::ExportScfTrace { trace, options } => {
            crate::scf_export::export(&trace, &options).map(AnalysisResult::ScfTraceExported)
        }
        AnalysisJob::PrepareSavedRunReview { workspace, run_id } => {
            crate::run_review::prepare(&workspace, &run_id)
                .map(AnalysisResult::SavedRunReviewPrepared)
        }
        AnalysisJob::RenderMapSvg { request } => {
            crate::svg_export::render(&request).map(AnalysisResult::MapSvgExported)
        }
        AnalysisJob::InspectMaterialMismatchPcr {
            source,
            origin,
            forward,
            reverse,
            options,
        } => crate::material_derivation::preview_mismatch_pcr(
            &source, &origin, &forward, &reverse, &options,
        )
        .map(AnalysisResult::MaterialMismatchPcrPreview),
        AnalysisJob::InspectMaterialRepair {
            source,
            origin,
            method,
        } => crate::material_derivation::preview_repair(&source, &origin, method)
            .map(AnalysisResult::MaterialRepairPreview),
        AnalysisJob::InspectMaterialOverlapExtension {
            workspace,
            source,
            origin,
            preparation,
        } => crate::material_derivation::preview_overlap_extension(
            &workspace,
            &source,
            &origin,
            &preparation,
        )
        .map(AnalysisResult::MaterialOverlapExtensionPreview),
        AnalysisJob::InspectMaterialTerminalProcessing {
            workspace,
            source,
            origin,
            action,
        } => crate::material_derivation::preview_terminal_processing(
            &workspace, &source, &origin, &action,
        )
        .map(AnalysisResult::MaterialTerminalProcessingPreview),
        AnalysisJob::InspectMaterialAnnealedOligos {
            workspace,
            top,
            bottom,
            options,
        } => {
            crate::material_derivation::preview_annealed_oligos(&workspace, &top, &bottom, &options)
                .map(AnalysisResult::MaterialAnnealedOligosPreview)
        }
        AnalysisJob::AnalyzeAlignmentAntibody {
            alignment,
            source,
            options,
            limits,
            selection,
        } => crate::alignment_antibody::analyze(
            &alignment,
            &source,
            &options,
            &limits,
            &selection,
            &mut || false,
        )
        .map(AnalysisResult::AlignmentAntibody),
        AnalysisJob::AnalyzeAntibody {
            molecule,
            options,
            limits,
            selection,
        } => crate::antibody_numbering::domain_tracks::analyze(
            &molecule,
            &options,
            &limits,
            &selection,
            &mut || false,
        )
        .map(AnalysisResult::Antibody)
        .map_err(|e| Error::Format(e.to_string())),
        AnalysisJob::ScreenProteinLiabilities { molecule, options } => {
            crate::protein_liabilities::screen(&molecule, &options, &mut || false)
                .map(AnalysisResult::ProteinLiabilities)
                .map_err(|e| Error::Format(e.to_string()))
        }
        AnalysisJob::AnalyzeRnaChemistry {
            molecule,
            declaration,
        } => crate::rna_chemistry::analyze(&molecule, &declaration)
            .map(AnalysisResult::RnaChemistry)
            .map_err(|e| Error::Format(e.to_string())),
        AnalysisJob::PrepareStructure {
            bytes,
            format,
            selection,
        } => crate::structure_view::prepare(&bytes, format, selection)
            .map(AnalysisResult::StructurePrepared),
        AnalysisJob::FindSimilarSequences {
            query,
            subjects,
            options,
        } => crate::sequence_similarity::search(&query, &subjects, &options, &mut || false)
            .map(AnalysisResult::SequenceSimilarity),
        AnalysisJob::InstallMaterialDerivation { workspace, request } => {
            crate::material_derivation::install(&workspace, &request)
                .map(AnalysisResult::MaterialDerivationInstalled)
        }
        AnalysisJob::InspectMaterialDigest {
            source_role,
            source,
            origin,
            selection,
            options,
        } => crate::material_derivation::inspect_digest(
            &source_role,
            &source,
            &origin,
            &selection,
            &options,
        )
        .map(AnalysisResult::MaterialDigestInspection),

        AnalysisJob::InspectTopoRun {
            inputs,
            template,
            task,
        } => crate::topo_run_preparation::inspect(&inputs, template.as_ref(), &task)
            .map(AnalysisResult::TopoRunInspection),
        AnalysisJob::ProcessTopoActions { inputs, request } => {
            crate::topo_actions::inspect(&inputs, &request).map(AnalysisResult::TopoActionsProduct)
        }

        AnalysisJob::LoadRunEvidence { workspace, run_id } => {
            crate::run_evidence::download(&workspace, &run_id)
                .map(AnalysisResult::HistoricalRunEvidence)
        }

        AnalysisJob::ListLigationProfiles => {
            crate::ligation_profiles::list().map(AnalysisResult::LigationProfiles)
        }
        AnalysisJob::ScorePublishedLigation {
            profile,
            ends,
            physical_competitors,
        } => {
            crate::ligation_profiles::score_with_competitors(&profile, &ends, &physical_competitors)
                .map(AnalysisResult::PublishedLigationScore)
        }
        AnalysisJob::ScoreLigationCompetition {
            dataset,
            ends,
            physical_competitors,
        } => crate::ligation_fidelity::score_competition(&dataset, &ends, &physical_competitors)
            .map(AnalysisResult::LigationCompetition),
        AnalysisJob::RankGoldenGateFusions { inputs, request } => {
            crate::golden_gate_fusion_rank::inspect(&inputs, &request)
                .map(AnalysisResult::RankedGoldenGateFusions)
        }
        AnalysisJob::InspectDomestication { source, request } => {
            crate::domestication::inspect(&source, &request)
                .map(AnalysisResult::DomesticationInspected)
        }
        AnalysisJob::ProposeDomestication { source, request } => {
            crate::domestication::propose(&source, &request)
                .map(AnalysisResult::DomesticationProposed)
        }
        AnalysisJob::CheckDomestication {
            source,
            request,
            selection,
        } => crate::domestication::check(&source, &request, &selection).map(|checked| {
            AnalysisResult::DomesticationChecked {
                report: checked.report().clone(),
                target: checked.target().clone(),
            }
        }),
        AnalysisJob::PrepareDomestication {
            input,
            request,
            workspace,
        } => crate::domestication_preparation::materialize(&input, &request, workspace.as_ref())
            .map(AnalysisResult::DomesticationPrepared),
        AnalysisJob::SearchGoldenGateFusions { inputs, request } => {
            crate::golden_gate_fusion_design::inspect(&inputs, &request)
                .map(AnalysisResult::GoldenGateFusions)
        }
        AnalysisJob::MaterializeGoldenGateFusions {
            inputs,
            request,
            selection,
        } => crate::golden_gate_fusion_design::materialize(&inputs, &request, &selection)
            .map(AnalysisResult::GoldenGateFusionMaterialized),
        AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
            crate::ligation_fidelity::score(&dataset, &ends).map(AnalysisResult::LigationFidelity)
        }
        AnalysisJob::FinalizeCloningOutputs {
            workspace,
            plan,
            target,
            inputs,
            options,
        } => crate::cloning_outputs::finalize(&workspace, &plan, target, inputs, &options)
            .map(AnalysisResult::CloningOutputsFinalized),
        AnalysisJob::SelectDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            expected_batch_revision,
            choice,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::digest_cut_selection::preview(
                &catalog,
                &selection,
                &source_role,
                &molecule,
                &options,
                &expected_batch_revision,
                &choice,
            )
            .map(AnalysisResult::DigestFragmentSelected)
        }
        AnalysisJob::ExportPlanOligos {
            plan,
            target,
            options,
            delimiter,
        } => crate::plan_oligos::export(&plan, &target, &options, delimiter)
            .map(AnalysisResult::PlanOligosExported),
        AnalysisJob::ExportProtocolHtml { plan, target } => Ok(AnalysisResult::ProtocolHtml(
            crate::verification::export_protocol_html(&plan, &target),
        )),
        AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin,
            request,
            workspace,
        } => crate::mismatch_pcr_digest_inspection::inspect(
            &molecule,
            origin.as_ref(),
            &request,
            workspace.as_ref(),
        )
        .map(AnalysisResult::MismatchPcrDigestInspection),
        AnalysisJob::PreviewAnnealingRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => crate::assembly_recipe::preview_annealing(&inputs, &target, &recipe, &policy)
            .map(AnalysisResult::OverlapRecipePreview),
        AnalysisJob::PreviewHomologyRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => crate::assembly_recipe::preview_homology(&inputs, &target, &recipe, &policy)
            .map(AnalysisResult::OverlapRecipePreview),
        AnalysisJob::DesignEditedOverlapPrimers {
            inputs,
            request,
            workspace,
        } => crate::assembly_edited_primer_design::design(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::EditedOverlapPrimersDesigned),
        AnalysisJob::InspectOverlapExtensionDigest {
            source,
            origin,
            request,
            workspace,
        } => crate::overlap_extension_digest_inspection::inspect(
            &source,
            origin.as_ref(),
            &request,
            workspace.as_ref(),
        )
        .map(AnalysisResult::OverlapExtensionDigestInspection),
        AnalysisJob::DesignOverlapExtensionPrimers {
            source,
            request,
            workspace,
        } => crate::overlap_extension_primer_design::design(&source, &request, workspace.as_ref())
            .map(AnalysisResult::OverlapExtensionPrimersDesigned),
        AnalysisJob::BuildContributionEditIntent {
            inputs,
            request,
            workspace,
        } => crate::assembly_edit_intent::build(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::ContributionEditIntent),
        AnalysisJob::BuildSubstitutedAssemblyIntent {
            inputs,
            request,
            workspace,
        } => crate::assembly_substitution_intent::build(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::SubstitutedAssemblyIntent),
        AnalysisJob::BuildAssemblyIntent {
            inputs,
            request,
            workspace,
        } => crate::assembly_intent::build(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::AssemblyIntent),
        AnalysisJob::BuildEditedAssemblyIntent {
            inputs,
            base,
            edits,
            workspace,
        } => {
            crate::assembly_junction_edit::build_intent(&inputs, &base, &edits, workspace.as_ref())
                .map(AnalysisResult::EditedAssemblyIntent)
        }
        AnalysisJob::DesignEditedAssemblyPrimers { inputs, request } => {
            crate::assembly_junction_edit::design(&inputs, &request)
                .map(AnalysisResult::EditedAssemblyPrimers)
        }
        AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request } => {
            crate::selective_junction_edit::design(&inputs, &request)
                .map(AnalysisResult::SelectiveAssemblyPrimers)
        }

        AnalysisJob::PreviewProtocolQuantities { workspace, request } => {
            crate::protocol_quantities::preview(&workspace, &request)
                .map(AnalysisResult::ProtocolQuantitiesPreview)
        }
        AnalysisJob::CaptureProtocolQuantities { workspace, request } => {
            crate::protocol_quantities::capture(&workspace, &request)
                .map(AnalysisResult::ProtocolQuantitiesCaptured)
        }
        AnalysisJob::ListTopoGraphs { workspace } => {
            crate::topo_graph_history::list(&workspace).map(AnalysisResult::TopoGraphsListed)
        }
        AnalysisJob::CaptureTopoGraph { workspace, request } => {
            crate::topo_graph_history::prepare_capture(&workspace, &request)
                .map(AnalysisResult::TopoGraphCaptured)
        }
        AnalysisJob::InspectTopoGraph {
            workspace,
            selection,
        } => crate::topo_graph_history::inspect(&workspace, &selection)
            .map(AnalysisResult::TopoGraphInspected),
        AnalysisJob::ListProtocolQuantities { workspace } => {
            crate::protocol_quantities::list(&workspace)
                .map(AnalysisResult::ProtocolQuantitiesListed)
        }
        AnalysisJob::OpenProtocolQuantities {
            workspace,
            quantity_hash,
        } => crate::protocol_quantities::open(&workspace, &quantity_hash)
            .map(AnalysisResult::ProtocolQuantitiesOpened),
        AnalysisJob::ExportProtocolQuantities {
            workspace,
            quantity_hash,
            format,
        } => crate::protocol_quantities::export(&workspace, &quantity_hash, format)
            .map(AnalysisResult::ProtocolQuantitiesExported),
        AnalysisJob::PreviewProtocolLayout { workspace, request } => {
            crate::protocol_layout::preview(&workspace, &request)
                .map(AnalysisResult::ProtocolLayoutPreview)
        }
        AnalysisJob::CaptureProtocolLayout { workspace, request } => {
            crate::protocol_layout::capture(&workspace, &request)
                .map(AnalysisResult::ProtocolLayoutCaptured)
        }
        AnalysisJob::ExportProtocolLayout {
            workspace,
            layout_hash,
            format,
        } => crate::protocol_layout::export(&workspace, &layout_hash, format)
            .map(AnalysisResult::ProtocolLayoutExported),
        AnalysisJob::InspectGoldenGateReaction {
            inputs,
            target,
            request,
        } => crate::golden_gate_reaction::inspect(&inputs, &target, &request)
            .map(AnalysisResult::GoldenGateReaction),
        AnalysisJob::PreparePcrOutputs { request } => {
            crate::pcr_outputs::prepare(&request).map(AnalysisResult::PcrOutputs)
        }
        AnalysisJob::PrepareMismatchPcrOutputs { request } => {
            crate::mismatch_pcr_outputs::prepare(&request).map(AnalysisResult::MismatchPcrOutputs)
        }
        AnalysisJob::CompareSequences {
            reference,
            query,
            options,
        } => crate::pairwise_review::compare(&reference, &query, &options)
            .map(AnalysisResult::ReviewedComparison),
        AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            fresh,
            format,
        } => crate::pairwise_export::export(&record, &reference_name, &query_name, fresh, format)
            .map(AnalysisResult::ExportedPairwiseComparison),
        AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } => crate::pairwise_search::search(&report, &pattern, options)
            .map(AnalysisResult::SearchedPairwiseComparison),
        AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } => crate::pairwise_context::prepare(&report, &reference, &query, translate)
            .map(AnalysisResult::PreparedPairwiseContext),
        AnalysisJob::CopyReferenceFragment { request } => {
            crate::reference_transfer::copy(&request).map(AnalysisResult::ReferenceFragment)
        }
        AnalysisJob::PasteReferenceFragment { request } => {
            crate::reference_transfer::paste(&request).map(AnalysisResult::ReferencePaste)
        }
        AnalysisJob::TranslateReference { request } => {
            crate::translation_reference::inspect(&request).map(AnalysisResult::TranslatedReference)
        }
        AnalysisJob::PrepareTranslatedReference { request, name } => {
            crate::translation_reference::prepare(&request, &name)
                .map(AnalysisResult::TranslatedOutputs)
        }
        AnalysisJob::SelectGoldenGateReaction { report, selection } => {
            crate::golden_gate_reaction::select(&report, &selection)
                .map(AnalysisResult::GoldenGateReactionCandidate)
        }
        AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request } => {
            crate::golden_gate_junction_edit::design(&inputs, &request)
                .map(AnalysisResult::GoldenGateJunctionPrimers)
        }
        AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request,
        } => crate::prepared_golden_gate::inspect(&inputs, &target, &request)
            .map(AnalysisResult::PreparedGoldenGatePool),
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, selection } => {
            crate::prepared_golden_gate::select(&report, &selection)
                .map(AnalysisResult::PreparedGoldenGateCandidate)
        }
        AnalysisJob::PreviewReadReview { workspace, request } => {
            crate::read_review::preview(&workspace, &request).map(AnalysisResult::ReadReviewPreview)
        }
        AnalysisJob::CommitReadReview { workspace, request } => {
            crate::read_review::commit(&workspace, &request)
                .map(AnalysisResult::ReadReviewCommitted)
        }
        AnalysisJob::CaptureReadCheck {
            workspace,
            reference,
            reads,
            options,
        } => crate::read_checks::capture(&workspace, reference, reads, &options)
            .map(AnalysisResult::ReadCheckCaptured),
        AnalysisJob::ReplayReadCheck {
            workspace,
            attachment_hash,
        } => crate::read_checks::replay(&workspace, &attachment_hash)
            .map(AnalysisResult::ReadCheckReplayed),
        AnalysisJob::BuildMutationIntent { molecule, request } => {
            crate::mutation_intent::build(&molecule, &request).map(AnalysisResult::MutationIntent)
        }
        AnalysisJob::DesignGoldenGatePrimers {
            inputs,
            target,
            request,
        } => crate::golden_gate_primer_design::design(&inputs, &target, &request)
            .map(AnalysisResult::GoldenGatePrimersDesigned),
        AnalysisJob::DesignOverlapPrimers {
            inputs,
            target,
            request,
        } => crate::overlap_primer_design::design(&inputs, &target, &request)
            .map(AnalysisResult::OverlapPrimersDesigned),
        AnalysisJob::InspectAssemblyJunctions {
            plan,
            target,
            options,
        } => crate::assembly_junctions::inspect_plan(&plan, &target, &options)
            .map(AnalysisResult::AssemblyJunctionInspection),
        AnalysisJob::InspectRunMaterial {
            workspace,
            event,
            selection,
        } => crate::run_material::inspect(&workspace, &event, &selection)
            .map(AnalysisResult::RunMaterialInspection),
        AnalysisJob::InspectRunCatalogRestriction {
            source_role,
            molecule,
            origin,
            selection,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::run_catalog_digest::inspect(
                &catalog,
                &selection,
                &source_role,
                &molecule,
                &origin,
                &options,
            )
            .map(AnalysisResult::RunCatalogRestrictionInspection)
        }
        AnalysisJob::InspectMaterialPcrDigest {
            molecule,
            origin,
            request,
            workspace,
        } => crate::pcr_digest_inspection::inspect_from_material(
            &molecule, &origin, &request, &workspace,
        )
        .map(AnalysisResult::MaterialPcrDigestInspection),
        AnalysisJob::InspectPcrDigest { molecule, request } => {
            crate::pcr_digest_inspection::inspect(&molecule, &request)
                .map(AnalysisResult::PcrDigestInspection)
        }
        AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        } => crate::tailing_recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::TailingRecipePreview),
        AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => crate::ligation_recipe::preview(&inputs, &target, &recipe, policy)
            .map(AnalysisResult::LigationRecipePreview),
        AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        } => crate::assembly_recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::OverlapRecipePreview),
        AnalysisJob::PrepareVariants { .. }
        | AnalysisJob::VariantPage { .. }
        | AnalysisJob::CheckVariantReference { .. }
        | AnalysisJob::PreviewVariantAllele { .. }
        | AnalysisJob::ExportVariants { .. }
        | AnalysisJob::DisposeVariants { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::SimulateTopo { inputs, request } => {
            crate::topo::simulate(&inputs, &request).map(AnalysisResult::TopoSimulation)
        }
        AnalysisJob::InspectTopoOrientations {
            inputs,
            request,
            vector,
            insert,
            directional,
        } => crate::topo::inspect_orientations(&inputs, &request, vector, insert, directional)
            .map(AnalysisResult::TopoOrientations),
        AnalysisJob::ProjectTopoProduct {
            inputs,
            request,
            selection,
        } => {
            crate::topo::product(&inputs, &request, &selection).map(AnalysisResult::TopoProjection)
        }
        AnalysisJob::CheckTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => crate::topo::check_target(&inputs, &request, &selection, &target, &policy)
            .map(AnalysisResult::TopoTarget),
        AnalysisJob::PrepareReadAlignments { .. }
        | AnalysisJob::ReadAlignmentPage { .. }
        | AnalysisJob::CheckReadAlignment { .. }
        | AnalysisJob::ExportReadAlignments { .. }
        | AnalysisJob::DisposeReadAlignments { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::InspectCombinatorialMultiSite {
            workspace,
            design,
            construct,
        } => crate::combinatorial_design::inspect_multisite(&design, &construct, &workspace)
            .map(AnalysisResult::CombinatorialMultiSiteInspection),
        AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
            crate::multisite_recombination::simulate(&inputs, &request)
                .map(AnalysisResult::MultiSiteRecombination)
        }
        AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection,
        } => crate::multisite_recombination::product(&inputs, &request, &selection)
            .map(AnalysisResult::MultiSiteRecombinationProduct),
        AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target,
        } => crate::multisite_recombination::check_target(&inputs, &request, &selection, &target)
            .map(AnalysisResult::MultiSiteRecombinationChecked),
        AnalysisJob::InspectRecombination { inputs, request } => {
            crate::recombination::simulate(&inputs, &request)
                .map(AnalysisResult::RecombinationInspection)
        }
        AnalysisJob::CheckRecombinationTarget {
            inputs,
            request,
            target,
            policy,
        } => crate::recombination::check_target(&inputs, &request, &target, &policy)
            .map(AnalysisResult::RecombinationTargetCheck),
        AnalysisJob::SelectRecombinationComponent {
            inputs,
            request,
            selection,
            target,
            policy,
        } => crate::recombination::checked_product(&inputs, &request, &selection, &target, &policy)
            .map(AnalysisResult::RecombinationProduct),
        AnalysisJob::InspectCatalogRestriction {
            selection,
            source_role,
            molecule,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::catalog_digest::inspect(&catalog, &selection, &source_role, &molecule, &options)
                .map(AnalysisResult::CatalogRestrictionInspection)
        }
        AnalysisJob::MaterializeCatalogDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::catalog_digest::product(
                &catalog,
                &selection,
                &source_role,
                &molecule,
                &options,
                fragment_index,
            )
            .map(AnalysisResult::CatalogDigestProduct)
        }
        AnalysisJob::InstallCatalogDigestFragment {
            selection,
            workspace,
            source,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::catalog_digest::install(
                &catalog,
                &selection,
                &workspace,
                &source,
                &options,
                fragment_index,
            )
            .map(AnalysisResult::CatalogDigestInstalled)
        }

        AnalysisJob::InstallDigestFragment {
            workspace,
            source,
            enzymes,
            options,
            fragment_index,
        } => crate::restriction_history::install(
            &workspace,
            &source,
            &enzymes,
            &options,
            fragment_index,
        )
        .map(AnalysisResult::DigestInstalled),
        AnalysisJob::ScanRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => crate::restriction_inspection::scan(&source_role, &molecule, &enzymes, &options)
            .map(AnalysisResult::RestrictionInspection),
        AnalysisJob::InspectRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => crate::restriction_inspection::inspect(&source_role, &molecule, &enzymes, &options)
            .map(AnalysisResult::RestrictionInspection),
        AnalysisJob::MaterializeDigestFragment {
            source_role,
            molecule,
            enzymes,
            options,
            fragment_index,
        } => crate::restriction_inspection::product(
            &source_role,
            &molecule,
            &enzymes,
            &options,
            fragment_index,
        )
        .map(AnalysisResult::DigestProduct),
        AnalysisJob::BeginGenomeFasta { .. }
        | AnalysisJob::PushGenomeFasta { .. }
        | AnalysisJob::FinishGenomeFasta { .. }
        | AnalysisJob::BeginGenomeNative { .. }
        | AnalysisJob::PushGenomeNativeBlock { .. }
        | AnalysisJob::FinishGenomeNative { .. }
        | AnalysisJob::AbortGenomeImport { .. }
        | AnalysisJob::RegisterGenome { .. }
        | AnalysisJob::GenomeIntervalGuides { .. }
        | AnalysisJob::StartGenomeSearch { .. }
        | AnalysisJob::AdvanceGenomeSearch { .. }
        | AnalysisJob::PageGenomeSearch { .. }
        | AnalysisJob::CancelGenomeJob { .. }
        | AnalysisJob::ReadGenomeInterval { .. } => Err(Error::Unsupported(
            "genome operations require a persistent private worker session".into(),
        )),
        AnalysisJob::DesignEditTemplate {
            source,
            target,
            options,
            guide,
        } => crate::edit_template::design(&source, &target, &options, guide.as_ref())
            .map(AnalysisResult::EditTemplate),
        AnalysisJob::CheckEditTemplate {
            source,
            target,
            options,
            guide,
            donor,
        } => crate::edit_template::check(&source, &target, &options, guide.as_ref(), &donor)
            .map(AnalysisResult::EditTemplateCheck),
        AnalysisJob::DesignCrisprGuides { molecule, options } => {
            crate::crispr::design_guides(&molecule, &options).map(AnalysisResult::CrisprGuides)
        }
        AnalysisJob::SearchCrisprOffTargets {
            query,
            references,
            options,
        } => crate::crispr::search_off_targets(&query, &references, &options)
            .map(AnalysisResult::CrisprOffTargets),
        AnalysisJob::AssessMethylation { molecule, options } => {
            crate::methylation::assess(&molecule, &options).map(AnalysisResult::Methylation)
        }
        AnalysisJob::AssembleContigs { reads, options } => {
            crate::contig::assemble(&reads, &options).map(AnalysisResult::Contigs)
        }
        AnalysisJob::DotPlot {
            molecule,
            query,
            options,
        } => crate::dotplot::analyze(&molecule, &query, &options).map(AnalysisResult::DotPlot),
        AnalysisJob::DesignCodons { molecule, options } => {
            crate::codon_design::design(&molecule, &options).map(AnalysisResult::CodonDesign)
        }
        AnalysisJob::CheckCodonDesign {
            source,
            candidate,
            options,
        } => crate::codon_design::check(&source, &candidate, &options)
            .map(AnalysisResult::CodonCheck),
        AnalysisJob::RecordCloningRun {
            workspace,
            plan,
            target_binding,
            input_bindings,
            allow_incomplete,
        } => crate::run_lineage::record_run(
            &workspace,
            &plan,
            target_binding,
            input_bindings,
            allow_incomplete,
        )
        .map(AnalysisResult::RunCommit),
        AnalysisJob::ReplayCloningRun { workspace, run_id } => {
            crate::run_lineage::replay_run(&workspace, &run_id).map(AnalysisResult::RunReplay)
        }
        AnalysisJob::TailPreview {
            molecule,
            recipe,
            base,
        } => crate::verification::simulate_tailing(&molecule, &recipe, base)
            .map(AnalysisResult::TailPreview),
        AnalysisJob::CheckOverlapProfile {
            molecules,
            reverse,
            overlaps,
            circular,
            options,
        } => crate::overlap_profiles::check_oriented_fragments(
            &molecules, &reverse, &overlaps, circular, &options,
        )
        .map(AnalysisResult::OverlapProfile),
        AnalysisJob::ProteinProperties { molecule, options } => {
            crate::protein_properties::analyze_protein(&molecule, &options)
                .map(AnalysisResult::ProteinProperties)
        }
        AnalysisJob::ImportAlignment {
            bytes,
            format,
            alphabet,
            name,
        } => crate::alignment_formats::import_alignment(&bytes, format, alphabet, &name)
            .map(AnalysisResult::AlignmentDocument),
        AnalysisJob::LoadAlignment { bytes } => {
            crate::alignment_formats::load_alignment(&bytes).map(AnalysisResult::AlignmentDocument)
        }
        AnalysisJob::ExportAlignment { document, format } => {
            crate::alignment_formats::export_alignment(&document, format)
                .map(AnalysisResult::AlignmentExport)
        }
        AnalysisJob::ExtractAlignmentRows { document, row_ids } => document
            .extract_rows(&row_ids)
            .map(AnalysisResult::AlignmentRows),
        AnalysisJob::SaveAlignment { document } => {
            crate::alignment_formats::save_alignment(&document).map(AnalysisResult::AlignmentSaved)
        }
        AnalysisJob::BeginCombinatorialBatch { .. }
        | AnalysisJob::AdvanceCombinatorialBatch { .. }
        | AnalysisJob::CaptureCombinatorialBatch { .. } => Err(Error::Unsupported(
            "combinatorial batches require a persistent private worker session".into(),
        )),
        AnalysisJob::BeginCombinatorialRelease { .. }
        | AnalysisJob::AdvanceCombinatorialRelease { .. }
        | AnalysisJob::InspectCombinatorialRelease { .. }
        | AnalysisJob::FinishCombinatorialRelease { .. } => Err(Error::Unsupported(
            "combinatorial release requires a persistent private worker session".into(),
        )),
        AnalysisJob::ExportCombinatorialOutputs {
            workspace,
            manifest_hash,
            options,
        } => crate::combinatorial_export::export(&workspace, &manifest_hash, &options)
            .map(AnalysisResult::CombinatorialOutputsArchive),
        AnalysisJob::PrepareTracks { .. } | AnalysisJob::CachedViewTracks { .. } => Err(
            Error::Unsupported("cached track jobs require an AnalysisSession".into()),
        ),
        AnalysisJob::ViewTracks {
            molecule,
            window,
            options,
        } => crate::view_tracks::view_tracks(&molecule, window, &options)
            .map(AnalysisResult::ViewTracks),
        AnalysisJob::SearchDocument { molecule, options } => {
            crate::search::search_document(&molecule, &options).map(AnalysisResult::DocumentSearch)
        }
        AnalysisJob::Search { molecule, query } => {
            crate::analysis::find_matches(&molecule, &query, true).map(AnalysisResult::Search)
        }
        AnalysisJob::Restriction { molecule, enzyme } => {
            crate::analysis::restriction_sites(&molecule, &enzyme).map(AnalysisResult::Restriction)
        }
        AnalysisJob::Pcr {
            molecule,
            forward,
            reverse,
        } => crate::analysis::simulate_pcr(&molecule, &forward, &reverse).map(AnalysisResult::Pcr),
        AnalysisJob::RepairPreview {
            molecule,
            enzymes,
            fragment_index,
            method,
        } => crate::verification::simulate_end_repair(&molecule, &enzymes, fragment_index, method)
            .map(AnalysisResult::RepairPreview),
        AnalysisJob::AnnealPreview {
            molecule,
            bottom,
            options,
        } => crate::verification::simulate_oligo_annealing(&molecule, &bottom, &options)
            .map(AnalysisResult::AnnealPreview),
        AnalysisJob::AnnealedDocument {
            molecule,
            bottom,
            preview,
            candidate_index,
        } => crate::verification::annealed_document(&molecule, &bottom, &preview, candidate_index)
            .map(AnalysisResult::PhysicalPreviewDocument),
        AnalysisJob::PcrDocuments {
            molecule,
            forward,
            reverse,
            template_strand,
        } => crate::verification::simulate_pcr_documents(
            &molecule,
            &forward,
            &reverse,
            template_strand,
        )
        .map(AnalysisResult::PcrDocuments),
        AnalysisJob::AnalyzePrimer {
            molecule,
            primer,
            conditions,
            binding_options,
        } => crate::primers::analyze_primer(&primer, &conditions).and_then(|analysis| {
            crate::primers::find_binding_sites(&molecule, &primer, &binding_options).map(
                |bindings| AnalysisResult::PrimerEvidence(PrimerEvidence { analysis, bindings }),
            )
        }),
        AnalysisJob::DesignPrimers { molecule, options } => {
            crate::primers::design_primers(&molecule, &options).map(AnalysisResult::PrimerDesign)
        }
        AnalysisJob::ExportProtocolPdf { target, plan } => {
            crate::pdf::export_protocol_pdf(&plan, &target).map(AnalysisResult::ProtocolPdf)
        }
        AnalysisJob::GoldenGate {
            templates, options, ..
        } => crate::golden_gate::simulate(&templates, &options).map(AnalysisResult::GoldenGate),
        AnalysisJob::MismatchPcr {
            molecule,
            forward,
            reverse,
            options,
        } => crate::mutagenesis::simulate_pcr(&molecule, &forward, &reverse, &options)
            .map(AnalysisResult::MismatchPcr),
        AnalysisJob::InversePcr {
            molecule,
            target,
            origin,
            options,
        } => crate::mutagenesis::design_inverse_pcr(&molecule, &target, origin.as_ref(), &options)
            .map(AnalysisResult::InversePcr),
        AnalysisJob::PredictGel {
            lanes, calibration, ..
        } => crate::gel::predict_gel(&lanes, &calibration).map(AnalysisResult::Gel),
        AnalysisJob::PredictPinnedGel {
            lanes, calibration, ..
        } => crate::gel::predict_gel(&lanes, &calibration).map(AnalysisResult::PinnedGel),
        AnalysisJob::DiagnosticDigests {
            candidates,
            options,
            ..
        } => crate::gel::diagnostic_digests(&candidates, &options)
            .map(AnalysisResult::DiagnosticDigest),
        AnalysisJob::MismatchPcrPreview {
            molecule,
            forward,
            reverse,
            options,
        } => crate::verification::simulate_mutagenesis_preview(
            &molecule, &forward, &reverse, &options,
        )
        .map(AnalysisResult::MismatchPcrPreview),
        AnalysisJob::DuplexDigest { molecule, enzymes } => {
            crate::chemistry::Duplex::from_molecule(&molecule)
                .and_then(|d| crate::chemistry::digest(&d, &enzymes))
                .map(AnalysisResult::DuplexDigest)
        }
        AnalysisJob::Translate { molecule, options } => {
            crate::translation::translate(&molecule.sequence, &options)
                .map(AnalysisResult::Translation)
        }
        AnalysisJob::TranslateCds {
            molecule,
            feature_index,
        } => molecule
            .features
            .get(feature_index)
            .ok_or_else(|| Error::Coordinates("CDS annotation no longer exists".into()))
            .and_then(|f| crate::translation::translate_cds(&molecule, f))
            .map(AnalysisResult::CdsTranslation),
        AnalysisJob::FindOrfs { molecule, options } => {
            crate::translation::find_orfs(&molecule, &options).map(AnalysisResult::Orfs)
        }
        AnalysisJob::Align {
            molecule,
            query,
            options,
        } => query
            .validate()
            .and_then(|_| crate::alignment::align(&molecule.sequence, &query.sequence, &options))
            .map(AnalysisResult::Alignment),
        AnalysisJob::Verify { target, plan } => Ok(AnalysisResult::Verification(
            crate::verification::verify(&plan, &target),
        )),
        AnalysisJob::SaveWorkspace { workspace } => {
            crate::workspace::save_workspace(&workspace).map(AnalysisResult::WorkspaceSaved)
        }
        AnalysisJob::LoadWorkspace { bytes } => {
            crate::workspace::load_workspace(&bytes).map(AnalysisResult::WorkspaceLoaded)
        }
        AnalysisJob::ValidateTraces {
            molecule,
            reads,
            options,
        } => crate::trace_validation::validate_reads(&molecule, &reads, &options)
            .map(AnalysisResult::TraceValidation),
        AnalysisJob::ScanFeatureLibrary {
            molecule,
            library,
            options,
        } => crate::feature_library::scan(&molecule, &library, &options)
            .map(AnalysisResult::FeatureLibraryScan),
        AnalysisJob::ApplyFeatureLibrary {
            molecule,
            library,
            report,
            selected,
        } => crate::feature_library::apply_selected(&molecule, &library, &report, &selected)
            .map(AnalysisResult::AnnotatedDocument),
        AnalysisJob::ImportAnnotations {
            molecule,
            text,
            format,
            seqid,
        } => crate::annotation_formats::import_annotations(&molecule, &text, format, &seqid)
            .map(AnalysisResult::AnnotationImport),
        AnalysisJob::ApplyAnnotations {
            molecule,
            report,
            selected,
        } => crate::annotation_formats::apply_annotations(&molecule, &report, &selected)
            .map(AnalysisResult::AnnotatedDocument),
        AnalysisJob::ExportAnnotations {
            molecule,
            format,
            seqid,
        } => crate::annotation_formats::export_annotations(&molecule, format, &seqid).map(|r| {
            AnalysisResult::AnnotationExport(AnnotationExportReport {
                text: r.text,
                warnings: r.warnings,
            })
        }),
    }
    .unwrap_or_else(|e| AnalysisResult::Error(e.to_string()));
    let maximum = if matches!(
        result,
        AnalysisResult::WorkspaceSaved(_)
            | AnalysisResult::TopoRunInspection(_)
            | AnalysisResult::WorkspaceLoaded(_)
            | AnalysisResult::RunCommit(_)
            | AnalysisResult::RunReplay(_)
            | AnalysisResult::CloningOutputsFinalized(_)
            | AnalysisResult::ReadCheckCaptured(_)
            | AnalysisResult::ReadReviewCommitted(_)
            | AnalysisResult::DigestInstalled(_)
            | AnalysisResult::CatalogDigestInstalled(_)
            | AnalysisResult::CombinatorialOutputsArchive(_)
            | AnalysisResult::MaterialPcrOutputsPrepared(_)
            | AnalysisResult::SequenceExportAssessment(_)
            | AnalysisResult::SequenceExportPrepared(_)
            | AnalysisResult::ScfTraceExported(_)
            | AnalysisResult::MaterialMismatchPcrPreview(_)
            | AnalysisResult::MaterialRepairPreview(_)
            | AnalysisResult::MaterialOverlapExtensionPreview(_)
            | AnalysisResult::MaterialTerminalProcessingPreview(_)
            | AnalysisResult::MaterialAnnealedOligosPreview(_)
            | AnalysisResult::MaterialDerivationInstalled(_)
            | AnalysisResult::MaterialDigestInspection(_)
            | AnalysisResult::SavedRunReviewPrepared(_)
            | AnalysisResult::ProtocolLayoutCaptured(_)
            | AnalysisResult::ProtocolQuantitiesCaptured(_)
            | AnalysisResult::TopoGraphCaptured(_)
            | AnalysisResult::StructurePrepared(_)
    ) {
        MAX_WORKSPACE_MESSAGE_BYTES
    } else {
        MAX_IMPORT_BYTES
    };
    encode_with_limit(
        b"CGR\0",
        &Response {
            id: r.id,
            revision: r.revision,
            result,
        },
        maximum,
    )
}
