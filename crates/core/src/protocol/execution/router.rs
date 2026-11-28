//! Exhaustive typed request routing. Admission and actual operations stay together in each domain.
use super::{
    assembly, context::RequestContext, documents, library, material, output, preparation,
    restriction, review, sequence, topology, tracks,
};
use crate::protocol::schema::{AnalysisJob, AnalysisResult};
use collagene_model::Error;

pub(super) fn dispatch(
    job: AnalysisJob,
    context: &RequestContext,
) -> Result<AnalysisResult, Error> {
    match job {
        AnalysisJob::PrepareTraceDisplay { molecule } => {
            sequence::prepare_trace_display(context, molecule)
        }
        AnalysisJob::Search { molecule, query } => sequence::search(context, molecule, query),
        AnalysisJob::Restriction { molecule, enzyme } => {
            preparation::restriction(context, molecule, enzyme)
        }
        AnalysisJob::Pcr {
            molecule,
            forward,
            reverse,
        } => preparation::pcr(context, molecule, forward, reverse),
        AnalysisJob::Verify { target, plan } => output::verify(context, target, plan),
        AnalysisJob::DuplexDigest { molecule, enzymes } => {
            preparation::duplex_digest(context, molecule, enzymes)
        }
        AnalysisJob::Translate { molecule, options } => {
            sequence::translate(context, molecule, options)
        }
        AnalysisJob::TranslateCds {
            molecule,
            feature_index,
        } => sequence::translate_cds(context, molecule, feature_index),
        AnalysisJob::FindOrfs { molecule, options } => {
            sequence::find_orfs(context, molecule, options)
        }
        AnalysisJob::Align {
            molecule,
            query,
            options,
        } => review::align(context, molecule, query, options),
        AnalysisJob::PcrDocuments {
            molecule,
            forward,
            reverse,
            template_strand,
        } => preparation::pcr_documents(context, molecule, forward, reverse, template_strand),
        AnalysisJob::AnalyzePrimer {
            molecule,
            primer,
            conditions,
            binding_options,
        } => preparation::analyze_primer(context, molecule, primer, conditions, binding_options),
        AnalysisJob::DesignPrimers { molecule, options } => {
            preparation::design_primers(context, molecule, options)
        }
        AnalysisJob::ExportProtocolPdf { target, plan } => {
            output::export_protocol_pdf(context, target, plan)
        }
        AnalysisJob::GoldenGate {
            molecule,
            templates,
            options,
        } => preparation::golden_gate(context, molecule, templates, options),
        AnalysisJob::MismatchPcr {
            molecule,
            forward,
            reverse,
            options,
        } => preparation::mismatch_pcr(context, molecule, forward, reverse, options),
        AnalysisJob::InversePcr {
            molecule,
            target,
            origin,
            options,
        } => preparation::inverse_pcr(context, molecule, target, origin, options),
        AnalysisJob::PredictGel {
            molecule,
            lanes,
            calibration,
        } => output::predict_gel(context, molecule, lanes, calibration),
        AnalysisJob::DiagnosticDigests {
            molecule,
            candidates,
            options,
        } => output::diagnostic_digests(context, molecule, candidates, options),
        AnalysisJob::MismatchPcrPreview {
            molecule,
            forward,
            reverse,
            options,
        } => preparation::mismatch_pcr_preview(context, molecule, forward, reverse, options),
        AnalysisJob::SaveWorkspace { workspace } => documents::save_workspace(context, workspace),
        AnalysisJob::LoadWorkspace { bytes } => documents::load_workspace(context, bytes),
        AnalysisJob::ValidateTraces {
            molecule,
            reads,
            options,
        } => review::validate_traces(context, molecule, reads, options),
        AnalysisJob::ScanFeatureLibrary {
            molecule,
            library,
            options,
        } => sequence::scan_feature_library(context, molecule, library, options),
        AnalysisJob::ApplyFeatureLibrary {
            molecule,
            library,
            report,
            selected,
        } => sequence::apply_feature_library(context, molecule, library, report, selected),
        AnalysisJob::ImportAnnotations {
            molecule,
            text,
            format,
            seqid,
        } => sequence::import_annotations(context, molecule, text, format, seqid),
        AnalysisJob::ApplyAnnotations {
            molecule,
            report,
            selected,
        } => sequence::apply_annotations(context, molecule, report, selected),
        AnalysisJob::ExportAnnotations {
            molecule,
            format,
            seqid,
        } => sequence::export_annotations(context, molecule, format, seqid),
        AnalysisJob::RepairPreview {
            molecule,
            enzymes,
            fragment_index,
            method,
        } => preparation::repair_preview(context, molecule, enzymes, fragment_index, method),
        AnalysisJob::AnnealPreview {
            molecule,
            bottom,
            options,
        } => preparation::anneal_preview(context, molecule, bottom, options),
        AnalysisJob::AnnealedDocument {
            molecule,
            bottom,
            preview,
            candidate_index,
        } => preparation::annealed_document(context, molecule, bottom, preview, candidate_index),
        AnalysisJob::ViewTracks {
            molecule,
            window,
            options,
        } => tracks::view_tracks(context, molecule, window, options),
        AnalysisJob::SearchDocument { molecule, options } => {
            sequence::search_document(context, molecule, options)
        }
        AnalysisJob::PrepareTracks { molecule } => tracks::prepare_tracks(context, molecule),
        AnalysisJob::CachedViewTracks { .. } => tracks::cached_view_tracks(context),
        AnalysisJob::ImportAlignment {
            bytes,
            format,
            alphabet,
            name,
        } => review::import_alignment(context, bytes, format, alphabet, name),
        AnalysisJob::ExportAlignment { document, format } => {
            review::export_alignment(context, document, format)
        }
        AnalysisJob::ExtractAlignmentRows { document, row_ids } => {
            review::extract_alignment_rows(context, document, row_ids)
        }
        AnalysisJob::SaveAlignment { document } => review::save_alignment(context, document),
        AnalysisJob::LoadAlignment { bytes } => review::load_alignment(context, bytes),
        AnalysisJob::ProteinProperties { molecule, options } => {
            sequence::protein_properties(context, molecule, options)
        }
        AnalysisJob::TailPreview {
            molecule,
            recipe,
            base,
        } => preparation::tail_preview(context, molecule, recipe, base),
        AnalysisJob::CheckOverlapProfile {
            molecules,
            reverse,
            overlaps,
            circular,
            options,
        } => assembly::check_overlap_profile(
            context, molecules, reverse, overlaps, circular, options,
        ),
        AnalysisJob::RecordCloningRun {
            workspace,
            plan,
            target_binding,
            input_bindings,
            allow_incomplete,
        } => material::record_cloning_run(
            context,
            workspace,
            plan,
            target_binding,
            input_bindings,
            allow_incomplete,
        ),
        AnalysisJob::ReplayCloningRun { workspace, run_id } => {
            material::replay_cloning_run(context, workspace, run_id)
        }
        AnalysisJob::DesignCodons { molecule, options } => {
            sequence::design_codons(context, molecule, options)
        }
        AnalysisJob::CheckCodonDesign {
            source,
            candidate,
            options,
        } => sequence::check_codon_design(context, source, candidate, options),
        AnalysisJob::DotPlot {
            molecule,
            query,
            options,
        } => review::dot_plot(context, molecule, query, options),
        AnalysisJob::AssembleContigs { reads, options } => {
            review::assemble_contigs(context, reads, options)
        }
        AnalysisJob::AssessMethylation { molecule, options } => {
            preparation::assess_methylation(context, molecule, options)
        }
        AnalysisJob::DesignCrisprGuides { molecule, options } => {
            sequence::design_crispr_guides(context, molecule, options)
        }
        AnalysisJob::SearchCrisprOffTargets {
            query,
            references,
            options,
        } => sequence::search_crispr_off_targets(context, query, references, options),
        AnalysisJob::DesignEditTemplate {
            source,
            target,
            options,
            guide,
        } => sequence::design_edit_template(context, source, target, options, guide),
        AnalysisJob::CheckEditTemplate {
            source,
            target,
            options,
            guide,
            donor,
        } => sequence::check_edit_template(context, source, target, options, guide, donor),
        AnalysisJob::BeginGenomeFasta { .. } => tracks::begin_genome_fasta(context),
        AnalysisJob::PushGenomeFasta { .. } => tracks::push_genome_fasta(context),
        AnalysisJob::FinishGenomeFasta { .. } => tracks::finish_genome_fasta(context),
        AnalysisJob::BeginGenomeNative { .. } => tracks::begin_genome_native(context),
        AnalysisJob::PushGenomeNativeBlock { .. } => tracks::push_genome_native_block(context),
        AnalysisJob::FinishGenomeNative { .. } => tracks::finish_genome_native(context),
        AnalysisJob::AbortGenomeImport { .. } => tracks::abort_genome_import(context),
        AnalysisJob::RegisterGenome { .. } => tracks::register_genome(context),
        AnalysisJob::GenomeIntervalGuides { .. } => tracks::genome_interval_guides(context),
        AnalysisJob::StartGenomeSearch { .. } => tracks::start_genome_search(context),
        AnalysisJob::AdvanceGenomeSearch { .. } => tracks::advance_genome_search(context),
        AnalysisJob::PageGenomeSearch { .. } => tracks::page_genome_search(context),
        AnalysisJob::CancelGenomeJob { .. } => tracks::cancel_genome_job(context),
        AnalysisJob::ReadGenomeInterval { .. } => tracks::read_genome_interval(context),
        AnalysisJob::InspectRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => restriction::inspect_restriction(context, source_role, molecule, enzymes, options),
        AnalysisJob::MaterializeDigestFragment {
            source_role,
            molecule,
            enzymes,
            options,
            fragment_index,
        } => restriction::materialize_digest_fragment(
            context,
            source_role,
            molecule,
            enzymes,
            options,
            fragment_index,
        ),
        AnalysisJob::InstallDigestFragment {
            workspace,
            source,
            enzymes,
            options,
            fragment_index,
        } => restriction::install_digest_fragment(
            context,
            workspace,
            source,
            enzymes,
            options,
            fragment_index,
        ),
        AnalysisJob::InspectCatalogRestriction {
            selection,
            source_role,
            molecule,
            options,
        } => restriction::inspect_catalog_restriction(
            context,
            selection,
            source_role,
            molecule,
            options,
        ),
        AnalysisJob::MaterializeCatalogDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            fragment_index,
        } => restriction::materialize_catalog_digest_fragment(
            context,
            selection,
            source_role,
            molecule,
            options,
            fragment_index,
        ),
        AnalysisJob::InstallCatalogDigestFragment {
            selection,
            workspace,
            source,
            options,
            fragment_index,
        } => restriction::install_catalog_digest_fragment(
            context,
            selection,
            workspace,
            source,
            options,
            fragment_index,
        ),
        AnalysisJob::InspectRecombination { inputs, request } => {
            topology::inspect_recombination(context, inputs, request)
        }
        AnalysisJob::CheckRecombinationTarget {
            inputs,
            request,
            target,
            policy,
        } => topology::check_recombination_target(context, inputs, request, target, policy),
        AnalysisJob::SelectRecombinationComponent {
            inputs,
            request,
            selection,
            target,
            policy,
        } => topology::select_recombination_component(
            context, inputs, request, selection, target, policy,
        ),
        AnalysisJob::ScanRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => restriction::scan_restriction(context, source_role, molecule, enzymes, options),
        AnalysisJob::PrepareReadAlignments { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::ReadAlignmentPage { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::CheckReadAlignment { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::ExportReadAlignments { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::DisposeReadAlignments { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::SimulateTopo { inputs, request } => {
            topology::simulate_topo(context, inputs, request)
        }
        AnalysisJob::InspectTopoOrientations {
            inputs,
            request,
            vector,
            insert,
            directional,
        } => topology::inspect_topo_orientations(
            context,
            inputs,
            request,
            vector,
            insert,
            directional,
        ),
        AnalysisJob::ProjectTopoProduct {
            inputs,
            request,
            selection,
        } => topology::project_topo_product(context, inputs, request, selection),
        AnalysisJob::CheckTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => topology::check_topo_target(context, inputs, request, selection, target, policy),
        AnalysisJob::PrepareVariants { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::VariantPage { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::CheckVariantReference { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::PreviewVariantAllele { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::ExportVariants { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::DisposeVariants { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        } => assembly::preview_overlap_recipe(context, inputs, target, recipe),
        AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => assembly::preview_ligation_recipe(context, inputs, target, recipe, policy),
        AnalysisJob::InspectAssemblyJunctions {
            plan,
            target,
            options,
        } => assembly::inspect_assembly_junctions(context, plan, target, options),
        AnalysisJob::DesignOverlapPrimers {
            inputs,
            target,
            request,
        } => assembly::design_overlap_primers(context, inputs, target, request),
        AnalysisJob::BuildAssemblyIntent {
            inputs,
            request,
            workspace,
        } => assembly::build_assembly_intent(context, inputs, request, workspace),
        AnalysisJob::SelectDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            expected_batch_revision,
            choice,
        } => restriction::select_digest_fragment(
            context,
            selection,
            source_role,
            molecule,
            options,
            expected_batch_revision,
            choice,
        ),
        AnalysisJob::ExportPlanOligos {
            plan,
            target,
            options,
            delimiter,
        } => output::export_plan_oligos(context, plan, target, options, delimiter),
        AnalysisJob::ExportProtocolHtml { plan, target } => {
            output::export_protocol_html(context, plan, target)
        }
        AnalysisJob::FinalizeCloningOutputs {
            workspace,
            plan,
            target,
            inputs,
            options,
        } => material::finalize_cloning_outputs(context, workspace, plan, target, inputs, options),
        AnalysisJob::DesignGoldenGatePrimers {
            inputs,
            target,
            request,
        } => assembly::design_golden_gate_primers(context, inputs, target, request),
        AnalysisJob::BuildMutationIntent { molecule, request } => {
            assembly::build_mutation_intent(context, molecule, request)
        }
        AnalysisJob::CaptureReadCheck {
            workspace,
            reference,
            reads,
            options,
        } => review::capture_read_check(context, workspace, reference, reads, options),
        AnalysisJob::ReplayReadCheck {
            workspace,
            attachment_hash,
        } => review::replay_read_check(context, workspace, attachment_hash),
        AnalysisJob::PreviewReadReview { workspace, request } => {
            review::preview_read_review(context, workspace, request)
        }
        AnalysisJob::CommitReadReview { workspace, request } => {
            review::commit_read_review(context, workspace, request)
        }
        AnalysisJob::BuildEditedAssemblyIntent {
            inputs,
            base,
            edits,
            workspace,
        } => assembly::build_edited_assembly_intent(context, inputs, base, edits, workspace),
        AnalysisJob::DesignEditedAssemblyPrimers { inputs, request } => {
            assembly::design_edited_assembly_primers(context, inputs, request)
        }
        AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request } => {
            assembly::design_selective_assembly_primers(context, inputs, request)
        }
        AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request } => {
            assembly::design_golden_gate_junction_primers(context, inputs, request)
        }
        AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request,
        } => assembly::inspect_prepared_golden_gate_pool(context, inputs, target, request),
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, selection } => {
            assembly::select_prepared_golden_gate_candidate(context, report, selection)
        }
        AnalysisJob::InspectGoldenGateReaction {
            inputs,
            target,
            request,
        } => assembly::inspect_golden_gate_reaction(context, inputs, target, request),
        AnalysisJob::SelectGoldenGateReaction { report, selection } => {
            assembly::select_golden_gate_reaction(context, report, selection)
        }
        AnalysisJob::PreparePcrOutputs { request } => {
            preparation::prepare_pcr_outputs(context, request)
        }
        AnalysisJob::PrepareMismatchPcrOutputs { request } => {
            preparation::prepare_mismatch_pcr_outputs(context, request)
        }
        AnalysisJob::CompareSequences {
            reference,
            query,
            options,
        } => review::compare_sequences(context, reference, query, options),
        AnalysisJob::CopyReferenceFragment { request } => {
            documents::copy_reference_fragment(context, request)
        }
        AnalysisJob::PasteReferenceFragment { request } => {
            documents::paste_reference_fragment(context, request)
        }
        AnalysisJob::TranslateReference { request } => {
            documents::translate_reference(context, request)
        }
        AnalysisJob::PrepareTranslatedReference { request, name } => {
            documents::prepare_translated_reference(context, request, name)
        }
        AnalysisJob::PredictPinnedGel {
            molecule,
            lanes,
            calibration,
        } => output::predict_pinned_gel(context, molecule, lanes, calibration),
        AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            fresh,
            format,
        } => review::export_pairwise_comparison(
            context,
            record,
            reference_name,
            query_name,
            fresh,
            format,
        ),
        AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } => review::search_pairwise_comparison(context, report, pattern, options),
        AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } => review::prepare_pairwise_context(context, report, reference, query, translate),
        AnalysisJob::BeginCombinatorialBatch { .. } => Err(Error::Unsupported(
            "combinatorial batches require a persistent private worker session".into(),
        )),
        AnalysisJob::AdvanceCombinatorialBatch { .. } => Err(Error::Unsupported(
            "combinatorial batches require a persistent private worker session".into(),
        )),
        AnalysisJob::CaptureCombinatorialBatch { .. } => Err(Error::Unsupported(
            "combinatorial batches require a persistent private worker session".into(),
        )),
        AnalysisJob::BeginCombinatorialRelease { .. } => Err(Error::Unsupported(
            "combinatorial release requires a persistent private worker session".into(),
        )),
        AnalysisJob::AdvanceCombinatorialRelease { .. } => Err(Error::Unsupported(
            "combinatorial release requires a persistent private worker session".into(),
        )),
        AnalysisJob::InspectCombinatorialRelease { .. } => Err(Error::Unsupported(
            "combinatorial release requires a persistent private worker session".into(),
        )),
        AnalysisJob::FinishCombinatorialRelease { .. } => Err(Error::Unsupported(
            "combinatorial release requires a persistent private worker session".into(),
        )),
        AnalysisJob::ExportCombinatorialOutputs {
            workspace,
            manifest_hash,
            options,
        } => library::export_combinatorial_outputs(context, workspace, manifest_hash, options),
        AnalysisJob::SearchGoldenGateFusions { inputs, request } => {
            assembly::search_golden_gate_fusions(context, inputs, request)
        }
        AnalysisJob::MaterializeGoldenGateFusions {
            inputs,
            request,
            selection,
        } => assembly::materialize_golden_gate_fusions(context, inputs, request, selection),
        AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
            preparation::score_ligation_fidelity(context, dataset, ends)
        }
        AnalysisJob::InspectDomestication { source, request } => {
            assembly::inspect_domestication(context, source, request)
        }
        AnalysisJob::ProposeDomestication { source, request } => {
            assembly::propose_domestication(context, source, request)
        }
        AnalysisJob::CheckDomestication {
            source,
            request,
            selection,
        } => assembly::check_domestication(context, source, request, selection),
        AnalysisJob::PrepareDomestication {
            input,
            request,
            workspace,
        } => assembly::prepare_domestication(context, input, request, workspace),
        AnalysisJob::ListLigationProfiles => preparation::list_ligation_profiles(context),
        AnalysisJob::ScorePublishedLigation {
            profile,
            ends,
            physical_competitors,
        } => preparation::score_published_ligation(context, profile, ends, physical_competitors),
        AnalysisJob::RankGoldenGateFusions { inputs, request } => {
            assembly::rank_golden_gate_fusions(context, inputs, request)
        }
        AnalysisJob::ScoreLigationCompetition {
            dataset,
            ends,
            physical_competitors,
        } => preparation::score_ligation_competition(context, dataset, ends, physical_competitors),
        AnalysisJob::LoadRunEvidence { workspace, run_id } => {
            material::load_run_evidence(context, workspace, run_id)
        }
        AnalysisJob::ProcessTopoActions { inputs, request } => {
            topology::process_topo_actions(context, inputs, request)
        }
        AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
            topology::simulate_multi_site_recombination(context, inputs, request)
        }
        AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection,
        } => topology::project_multi_site_recombination(context, inputs, request, selection),
        AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target,
        } => topology::check_multi_site_recombination(context, inputs, request, selection, target),
        AnalysisJob::InspectCombinatorialMultiSite {
            workspace,
            design,
            construct,
        } => library::inspect_combinatorial_multi_site(context, workspace, design, construct),
        AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        } => assembly::preview_tailing_recipe(context, inputs, target, recipe),
        AnalysisJob::InspectPcrDigest { molecule, request } => {
            preparation::inspect_pcr_digest(context, molecule, request)
        }
        AnalysisJob::InspectRunMaterial {
            workspace,
            event,
            selection,
        } => material::inspect_run_material(context, workspace, event, selection),
        AnalysisJob::InspectMaterialPcrDigest {
            molecule,
            origin,
            request,
            workspace,
        } => {
            preparation::inspect_material_pcr_digest(context, molecule, origin, request, workspace)
        }
        AnalysisJob::InspectRunCatalogRestriction {
            source_role,
            molecule,
            origin,
            selection,
            options,
        } => restriction::inspect_run_catalog_restriction(
            context,
            source_role,
            molecule,
            origin,
            selection,
            options,
        ),
        AnalysisJob::BuildSubstitutedAssemblyIntent {
            inputs,
            request,
            workspace,
        } => assembly::build_substituted_assembly_intent(context, inputs, request, workspace),
        AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin,
            request,
            workspace,
        } => {
            preparation::inspect_mismatch_pcr_digest(context, molecule, origin, request, workspace)
        }
        AnalysisJob::BuildContributionEditIntent {
            inputs,
            request,
            workspace,
        } => assembly::build_contribution_edit_intent(context, inputs, request, workspace),
        AnalysisJob::DesignOverlapExtensionPrimers {
            source,
            request,
            workspace,
        } => preparation::design_overlap_extension_primers(context, source, request, workspace),
        AnalysisJob::InspectOverlapExtensionDigest {
            source,
            origin,
            request,
            workspace,
        } => preparation::inspect_overlap_extension_digest(
            context, source, origin, request, workspace,
        ),
        AnalysisJob::DesignEditedOverlapPrimers {
            inputs,
            request,
            workspace,
        } => assembly::design_edited_overlap_primers(context, inputs, request, workspace),
        AnalysisJob::PreviewHomologyRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => assembly::preview_homology_recipe(context, inputs, target, recipe, policy),
        AnalysisJob::PreviewAnnealingRecipe {
            inputs,
            target,
            recipe,
            policy,
        } => assembly::preview_annealing_recipe(context, inputs, target, recipe, policy),
        AnalysisJob::InspectTopoRun {
            inputs,
            template,
            task,
        } => topology::inspect_topo_run(context, inputs, template, task),
        AnalysisJob::InspectMaterialOrigin {
            workspace,
            event,
            selection,
        } => material::inspect_material_origin(context, workspace, event, selection),
        AnalysisJob::MaterialPcrDocuments {
            molecule,
            origin,
            forward,
            reverse,
            template_strand,
        } => material::material_pcr_documents(
            context,
            molecule,
            origin,
            forward,
            reverse,
            template_strand,
        ),
        AnalysisJob::PrepareMaterialPcrOutputs { origin, request } => {
            material::prepare_material_pcr_outputs(context, origin, request)
        }
        AnalysisJob::ReviewSequenceExport { molecule, format } => {
            documents::review_sequence_export(context, molecule, format)
        }
        AnalysisJob::PrepareSequenceExport {
            molecule,
            format,
            expected_context_revision,
            allow_omissions,
        } => documents::prepare_sequence_export(
            context,
            molecule,
            format,
            expected_context_revision,
            allow_omissions,
        ),
        AnalysisJob::ExportScfTrace { molecule, options } => {
            documents::export_scf_trace(context, molecule, options)
        }
        AnalysisJob::InstallMaterialDerivation { workspace, request } => {
            material::install_material_derivation(context, workspace, request)
        }
        AnalysisJob::InspectMaterialDigest {
            source_role,
            source,
            origin,
            selection,
            options,
        } => material::inspect_material_digest(
            context,
            source_role,
            source,
            origin,
            selection,
            options,
        ),
        AnalysisJob::PrepareSavedRunReview { workspace, run_id } => {
            material::prepare_saved_run_review(context, workspace, run_id)
        }
        AnalysisJob::PrepareImageEvidence {
            name,
            association,
            original,
        } => output::prepare_image_evidence(context, name, association, original),
        AnalysisJob::ReopenImageEvidence { prepared } => {
            output::reopen_image_evidence(context, prepared)
        }
        AnalysisJob::RenderMap { request } => output::render_map(context, request),
        AnalysisJob::RenderGelImage { request } => output::render_gel_image(context, request),
        AnalysisJob::InspectMaterialMismatchPcr {
            source,
            origin,
            forward,
            reverse,
            options,
        } => material::inspect_material_mismatch_pcr(
            context, source, origin, forward, reverse, options,
        ),
        AnalysisJob::InspectMaterialRepair {
            source,
            origin,
            method,
        } => material::inspect_material_repair(context, source, origin, method),
        AnalysisJob::PreviewProtocolLayout { workspace, request } => {
            output::preview_protocol_layout(context, workspace, request)
        }
        AnalysisJob::CaptureProtocolLayout { workspace, request } => {
            output::capture_protocol_layout(context, workspace, request)
        }
        AnalysisJob::ExportProtocolLayout {
            workspace,
            layout_hash,
            format,
        } => output::export_protocol_layout(context, workspace, layout_hash, format),
        AnalysisJob::PreviewProtocolQuantities { workspace, request } => {
            output::preview_protocol_quantities(context, workspace, request)
        }
        AnalysisJob::CaptureProtocolQuantities { workspace, request } => {
            output::capture_protocol_quantities(context, workspace, request)
        }
        AnalysisJob::OpenProtocolQuantities {
            workspace,
            quantity_hash,
        } => output::open_protocol_quantities(context, workspace, quantity_hash),
        AnalysisJob::ExportProtocolQuantities {
            workspace,
            quantity_hash,
            format,
        } => output::export_protocol_quantities(context, workspace, quantity_hash, format),
        AnalysisJob::ListProtocolQuantities { workspace } => {
            output::list_protocol_quantities(context, workspace)
        }
        AnalysisJob::ListTopoGraphs { workspace } => topology::list_topo_graphs(context, workspace),
        AnalysisJob::CaptureTopoGraph { workspace, request } => {
            topology::capture_topo_graph(context, workspace, request)
        }
        AnalysisJob::InspectTopoGraph {
            workspace,
            selection,
        } => topology::inspect_topo_graph(context, workspace, selection),
        AnalysisJob::InspectMaterialOverlapExtension {
            workspace,
            source,
            origin,
            preparation,
        } => material::inspect_material_overlap_extension(
            context,
            workspace,
            source,
            origin,
            preparation,
        ),
        AnalysisJob::InspectMaterialTerminalProcessing {
            workspace,
            source,
            origin,
            action,
        } => material::inspect_material_terminal_processing(
            context, workspace, source, origin, action,
        ),
        AnalysisJob::InspectMaterialAnnealedOligos {
            workspace,
            top,
            bottom,
            options,
        } => material::inspect_material_annealed_oligos(context, workspace, top, bottom, options),
        AnalysisJob::FindSimilarSequences {
            query,
            subjects,
            options,
        } => review::find_similar_sequences(context, query, subjects, options),
        AnalysisJob::PrepareStructure {
            bytes,
            format,
            selection,
        } => output::prepare_structure(context, bytes, format, selection),
        AnalysisJob::AnalyzeRnaChemistry {
            molecule,
            declaration,
        } => sequence::analyze_rna_chemistry(context, molecule, declaration),
        AnalysisJob::ScreenProteinLiabilities { molecule, options } => {
            sequence::screen_protein_liabilities(context, molecule, options)
        }
        AnalysisJob::AnalyzeAntibody {
            molecule,
            options,
            limits,
            selection,
        } => sequence::analyze_antibody(context, molecule, options, limits, selection),
        AnalysisJob::AnalyzeAlignmentAntibody {
            alignment,
            source,
            options,
            limits,
            selection,
        } => sequence::analyze_alignment_antibody(
            context, alignment, source, options, limits, selection,
        ),
        AnalysisJob::MaterializeMaterialDigestFragment { source, request } => {
            material::materialize_material_digest_fragment(context, source, request)
        }
    }
}
