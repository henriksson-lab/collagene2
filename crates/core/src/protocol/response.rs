//! Response ownership and validation after bounded decoding.
use super::{
    bindings::profile_binding,
    envelope::{decode_with_limit, worker_limit},
    schema::AnalysisResult,
    validation::validate_verification_report,
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Response {
    pub id: u32,
    pub revision: String,
    pub result: AnalysisResult,
}

impl Response {
    pub fn decode(bytes: &[u8]) -> Result<Response, Error> {
        // Read the selector only after checking its bounds; decode still validates
        // magic/version/length/checksum and the complete selected payload layout.
        let maximum = worker_limit(
            b"CGR\0",
            bytes,
            &[
                20, 21, 40, 41, 62, 65, 92, 95, 98, 120, 124, 125, 157, 158, 159, 160, 161, 162,
                163, 164, 165, 166, 168, 169, 171, 174, 176, 178, 179, 180, 182,
            ],
        )?;

        let response: Response = decode_with_limit(b"CGR\0", bytes, maximum)?;
        match &response.result {
            AnalysisResult::TopoGraphsListed(r) => r.validate_binding(&response.revision)?,
            AnalysisResult::TopoGraphCaptured(r) => r.validate_binding(&response.revision)?,
            AnalysisResult::TopoGraphInspected(r) => r.validate_binding(&response.revision)?,
            AnalysisResult::ProtocolQuantitiesListed(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ProtocolQuantitiesPreview(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ProtocolQuantitiesCaptured(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ProtocolQuantitiesOpened(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ProtocolQuantitiesExported(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ProtocolLayoutPreview(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::ProtocolLayoutCaptured(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::ProtocolLayoutExported(report) => {
                report.validate_binding(&response.revision)?
            }

            AnalysisResult::SavedRunReviewPrepared(report) => {
                report.validate_request_binding(&response.revision)?;
            }
            AnalysisResult::MapSvgExported(report) => {
                report.validate_binding(&report.source_revision)?;
            }
            AnalysisResult::MaterialOriginInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialMismatchPcrPreview(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialRepairPreview(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialOverlapExtensionPreview(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialTerminalProcessingPreview(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialAnnealedOligosPreview(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::AlignmentAntibody(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::Antibody(report) => report
                .validate_binding(&response.revision)
                .map_err(|e| Error::Format(e.to_string()))?,
            AnalysisResult::ProteinLiabilities(report) => report
                .validate_binding(&response.revision)
                .map_err(|e| Error::Format(e.to_string()))?,
            AnalysisResult::RnaChemistry(report) => report
                .validate_binding(&response.revision)
                .map_err(|e| Error::Format(e.to_string()))?,
            AnalysisResult::StructurePrepared(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::SequenceSimilarity(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialDerivationInstalled(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialDigestInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::SequenceExportAssessment(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::SequenceExportPrepared(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ScfTraceExported(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialPcrOutputsPrepared(report) => {
                report.validate_binding(&response.revision)?;
            }
            _ => {}
        }

        match &response.result {
            AnalysisResult::TopoRunInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::TopoActionsProduct(report) => {
                report.validate_binding(&response.revision)?
            }

            AnalysisResult::HistoricalRunEvidence(value) => {
                if value.run_id.len() != 64
                    || value.recorded_engine.len() > 128
                    || value.step_count > 128
                    || value
                        .bytes
                        .as_ref()
                        .is_some_and(|b| b.len() > crate::run_evidence::MAX_BYTES)
                    || value.bytes.is_none() != (value.step_count == 0)
                {
                    return Err(Error::Format(
                        "Invalid historical scientific-evidence download".into(),
                    ));
                }
            }

            AnalysisResult::LigationProfiles(items) => {
                let expected = crate::ligation_profiles::list()?;
                if items.len() != expected.len()
                    || response.revision != crate::ligation_profiles::catalog_revision()?
                    || crate::assembly_recipe::hash(items)?
                        != crate::assembly_recipe::hash(&expected)?
                {
                    return Err(Error::Format("Published ligation catalog changed".into()));
                }
            }
            AnalysisResult::PublishedLigationScore(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::RankedGoldenGateFusions(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::LigationCompetition(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::DomesticationInspected(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::DomesticationProposed(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::DomesticationChecked { report, target } => {
                report.validate_binding(&response.revision)?;
                target.validate()?;
                if report.target_revision != target.revision() {
                    return Err(Error::Format("Domestication reference changed".into()));
                }
            }
            AnalysisResult::DomesticationPrepared(report) => {
                validate_verification_report(&report.verification)?;

                report.validate_binding(&response.revision)?
            }
            AnalysisResult::GoldenGateFusions(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::GoldenGateFusionMaterialized(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::LigationFidelity(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::EditedAssemblyIntent(intent) => {
                intent.validate_binding(&response.revision)?
            }
            AnalysisResult::EditedAssemblyPrimers(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::SelectiveAssemblyPrimers(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::GoldenGateJunctionPrimers(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::PreparedGoldenGatePool(report) => {
                validate_verification_report(&report.preview.verification)?;

                report.validate_binding(&response.revision)?
            }
            AnalysisResult::PreparedGoldenGateCandidate(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::PcrOutputs(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MismatchPcrOutputs(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ReviewedComparison(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ExportedPairwiseComparison(export) => export.validate()?,
            AnalysisResult::SearchedPairwiseComparison(matches) => matches.validate()?,
            AnalysisResult::PreparedPairwiseContext(context) => {
                context.validate()?;
                if context.context_revision != response.revision {
                    return Err(Error::Format("Pairwise context response changed.".into()));
                }
            }
            AnalysisResult::CombinatorialBatchStarted(progress)
            | AnalysisResult::CombinatorialBatchAdvanced(progress) => {
                progress.validate()?;
                if progress.context_revision != response.revision {
                    return Err(Error::Format(
                        "Combinatorial batch response context changed.".into(),
                    ));
                }
            }
            AnalysisResult::CombinatorialBatchCaptured(commit) => {
                commit.validate_binding(&response.revision)?;
            }
            AnalysisResult::CombinatorialReleaseStarted(progress)
            | AnalysisResult::CombinatorialReleaseAdvanced(progress) => {
                progress.validate()?;
                if progress.context_revision != response.revision {
                    return Err(Error::Format(
                        "Combinatorial release response context changed.".into(),
                    ));
                }
            }
            AnalysisResult::CombinatorialReleaseInspected(inspection) => {
                inspection.validate()?;
                if inspection.context_revision != response.revision {
                    return Err(Error::Format(
                        "Combinatorial inspection response context changed.".into(),
                    ));
                }
            }
            AnalysisResult::CombinatorialReleaseFinished(commit) => {
                commit
                    .validate_binding(&response.revision, &commit.source_workspace_fingerprint)?;
            }
            AnalysisResult::CombinatorialOutputsArchive(archive) => {
                archive.validate()?;
                if archive.context_revision != response.revision {
                    return Err(Error::Format(
                        "Combinatorial archive response context changed.".into(),
                    ));
                }
            }
            AnalysisResult::ReferenceFragment(fragment) => {
                fragment.validate_binding(&response.revision)?;
            }
            AnalysisResult::ReferencePaste(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::TranslatedReference(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::TranslatedOutputs(outputs) => {
                outputs.validate_binding(&response.revision)?;
            }
            AnalysisResult::GoldenGateReaction(report) => {
                validate_verification_report(&report.verification)?;
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::GoldenGateReactionCandidate(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::ReadReviewPreview(preview) => {
                preview.validate_binding(&response.revision)?
            }
            AnalysisResult::ReadReviewCommitted(commit) => {
                commit.validate_binding(&response.revision)?;
            }
            AnalysisResult::ReadCheckCaptured(commit) => {
                commit.validate_binding(&response.revision)?;
            }
            AnalysisResult::ReadCheckReplayed(checked) => {
                checked.validate_binding(&response.revision)?
            }
            AnalysisResult::MutationIntent(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::GoldenGatePrimersDesigned(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::CloningOutputsFinalized(commit) => {
                commit.validate_binding(&response.revision)?;
            }
            AnalysisResult::DigestFragmentSelected(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::PlanOligosExported(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::AssemblyIntent(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::OverlapPrimersDesigned(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::OverlapExtensionPrimersDesigned(report) => {
                validate_verification_report(&report.verification)?;
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::OverlapExtensionDigestInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::EditedOverlapPrimersDesigned(report) => {
                validate_verification_report(&report.preview.verification)?;
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::AssemblyJunctionInspection(report) => {
                validate_verification_report(&report.verification)?;
                if matches!(
                    report.method,
                    Some(crate::assembly_junctions::Method::GoldenGate { .. })
                ) {}
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::MismatchPcrDigestInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::ContributionEditIntent(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::SubstitutedAssemblyIntent(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::RunCatalogRestrictionInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MaterialPcrDigestInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::RunMaterialInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::PcrDigestInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::TailingRecipePreview(report) => {
                validate_verification_report(&report.verification)?;
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::LigationRecipePreview(report) => {
                validate_verification_report(&report.verification)?;
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::OverlapRecipePreview(report) => {
                validate_verification_report(&report.verification)?;
                report.validate()?;
                if report.context_revision != response.revision {
                    return Err(Error::Format(
                        "overlap preparation response context changed".into(),
                    ));
                }
            }
            AnalysisResult::VariantsPrepared(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format("variant response context changed".into()));
                }
            }
            AnalysisResult::VariantPage(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format("variant response context changed".into()));
                }
            }
            AnalysisResult::VariantReferenceChecked(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format("variant response context changed".into()));
                }
            }
            AnalysisResult::VariantAllelePreview(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format("variant response context changed".into()));
                }
            }
            AnalysisResult::VariantsExported(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format("variant response context changed".into()));
                }
            }
            AnalysisResult::VariantsDisposed(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format("variant response context changed".into()));
                }
            }
            AnalysisResult::TopoSimulation(report) => {
                crate::topo_worker::validate_simulation(report, &response.revision)?
            }
            AnalysisResult::TopoOrientations(report) => {
                crate::topo_worker::validate_orientations(report, &response.revision)?
            }
            AnalysisResult::TopoProjection(report) => {
                crate::topo_worker::validate_projection(report, &response.revision)?
            }
            AnalysisResult::TopoTarget(report) => {
                crate::topo_worker::validate_target(report, &response.revision)?
            }
            AnalysisResult::ReadAlignmentsPrepared(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format(
                        "read alignment preparation response mismatch".into(),
                    ));
                }
            }
            AnalysisResult::ReadAlignmentPage(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format(
                        "read alignment page response mismatch".into(),
                    ));
                }
            }
            AnalysisResult::ReadAlignmentChecked(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format(
                        "read alignment check response mismatch".into(),
                    ));
                }
            }
            AnalysisResult::ReadAlignmentsExported(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format(
                        "read alignment export response mismatch".into(),
                    ));
                }
            }
            AnalysisResult::ReadAlignmentsDisposed(report) => {
                report.validate()?;
                if report.revision != response.revision {
                    return Err(Error::Format(
                        "read alignment disposal response mismatch".into(),
                    ));
                }
            }
            AnalysisResult::CombinatorialMultiSiteInspection(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MultiSiteRecombination(report) => {
                report.validate_binding(&response.revision)?;
            }
            AnalysisResult::MultiSiteRecombinationProduct(report) => {
                report.validate_binding(&report.selection.batch_revision)?;
                if report.context_revision()? != response.revision {
                    return Err(Error::Format(
                        "multisite selected-component reply context mismatch".into(),
                    ));
                }
            }
            AnalysisResult::MultiSiteRecombinationChecked(report) => {
                report.validate_binding(&report.product.selection.batch_revision)?;
                if report.context_revision()? != response.revision {
                    return Err(Error::Format(
                        "multisite independent-target reply context mismatch".into(),
                    ));
                }
            }
            AnalysisResult::RecombinationInspection(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "recombination response binding mismatch".into(),
                    ));
                }
            }
            AnalysisResult::RecombinationTargetCheck(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "recombination response binding mismatch".into(),
                    ));
                }
            }
            AnalysisResult::RecombinationProduct(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "recombination response binding mismatch".into(),
                    ));
                }
            }
            AnalysisResult::CatalogRestrictionInspection(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "catalog restriction batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::CatalogDigestProduct(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "catalog product batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::CatalogDigestInstalled(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "catalog installation batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::DigestInstalled(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "installed digest batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::RestrictionInspection(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "restriction inspection batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::DigestProduct(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "digest product batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::GenomeImportStarted(_)
            | AnalysisResult::GenomeImportProgress(_)
            | AnalysisResult::GenomeFastaFinished(_)
            | AnalysisResult::GenomeRegistered(_)
            | AnalysisResult::GenomeIntervalGuides(_)
            | AnalysisResult::GenomeSearchState(_)
            | AnalysisResult::GenomeSearchPage(_)
            | AnalysisResult::GenomeDisposed(_)
            | AnalysisResult::GenomeIntervalRead(_) => {
                crate::genome_session::validate_result(&response.result, &response.revision)?
            }
            AnalysisResult::InversePcr(report) => {
                report.validate_binding(&response.revision)?;
                validate_verification_report(&report.verification)?;
            }
            AnalysisResult::Verification(report) => {
                validate_verification_report(report)?;
            }
            AnalysisResult::EditTemplate(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "edit-template batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::EditTemplateCheck(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "edit-template check batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::CrisprGuides(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format("CRISPR guide batch revision mismatch".into()));
                }
            }
            AnalysisResult::CrisprOffTargets(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format(
                        "CRISPR off-target batch revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::Methylation(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format("methylation batch revision mismatch".into()));
                }
            }
            AnalysisResult::Contigs(report) => {
                report.validate()?;
                if report.batch_revision != response.revision {
                    return Err(Error::Format("contig batch revision mismatch".into()));
                }
            }
            AnalysisResult::DotPlot(report) => {
                report.validate()?;
                if report.reference_revision != response.revision {
                    return Err(Error::Format("dot plot source revision mismatch".into()));
                }
            }
            AnalysisResult::CodonDesign(report) => {
                report.validate()?;
                if report.source_revision != response.revision {
                    return Err(Error::Format(
                        "codon design source revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::CodonCheck(report) => {
                report.validate()?;
                if report.source_revision != response.revision {
                    return Err(Error::Format("codon check source revision mismatch".into()));
                }
            }
            AnalysisResult::RunCommit(commit) => {
                commit.validate()?;
                if commit.receipt.source_workspace_fingerprint != response.revision {
                    return Err(Error::Format("run commit source workspace mismatch".into()));
                }
            }
            AnalysisResult::RunReplay(report) => {
                validate_verification_report(&report.current_report)?;
                report.validate()?;
                if report.source_workspace_fingerprint != response.revision {
                    return Err(Error::Format("run replay source workspace mismatch".into()));
                }
            }
            AnalysisResult::TailPreview(preview) => {
                preview.molecule.validate()?;
                preview.tailing.validate()?;
                if preview.source_revision != response.revision
                    || preview.base != preview.tailing.base
                    || preview.molecule.sequence != preview.tailing.duplex.top
                {
                    return Err(Error::Format(
                        "tailing preview response binding mismatch".into(),
                    ));
                }
            }
            AnalysisResult::OverlapProfile(report) => {
                report.validate()?;
                if profile_binding(
                    &report.source_revisions,
                    &report.orientations,
                    &report.overlaps,
                    report.circular,
                    &report.options,
                )? != response.revision
                {
                    return Err(Error::Format(
                        "overlap profile response binding mismatch".into(),
                    ));
                }
            }
            AnalysisResult::ProteinProperties(report) => {
                report.validate()?;
                if report.source_revision != response.revision {
                    return Err(Error::Format(
                        "protein-property response revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::AlignmentDocument(document) => {
                document.validate()?;
            }
            AnalysisResult::AlignmentRows(rows) => {
                if rows.len() > crate::alignment_formats::MAX_ALIGNMENT_ROWS {
                    return Err(Error::TooLarge);
                }
                for row in rows {
                    row.validate()?;
                }
            }
            AnalysisResult::AlignmentSaved(archive) => {
                crate::alignment_formats::load_alignment(archive)?;

                if blake3::hash(archive).to_hex().as_str() != response.revision {
                    return Err(Error::Format(
                        "saved alignment response revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::AlignmentExport(report) => {
                if report.text.len() > crate::alignment_formats::MAX_ALIGNMENT_BYTES
                    || report.losses.len() > crate::alignment_formats::MAX_ALIGNMENT_ANNOTATIONS
                    || report.renamings.len() > crate::alignment_formats::MAX_ALIGNMENT_ROWS
                {
                    return Err(Error::TooLarge);
                }
            }
            AnalysisResult::TracksPrepared(revision) => {
                if revision != &response.revision
                    || revision.len() != 64
                    || !revision.bytes().all(|b| b.is_ascii_hexdigit())
                {
                    return Err(Error::Format(
                        "prepared-track response revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::DocumentSearch(report) => {
                report.validate()?;
                if report.source_revision != response.revision {
                    return Err(Error::Format("search response revision mismatch".into()));
                }
            }
            AnalysisResult::ViewTracks(report) => {
                report.validate()?;
                if report.source_revision != response.revision {
                    return Err(Error::Format("track response revision mismatch".into()));
                }
            }
            AnalysisResult::RepairPreview(preview) => {
                preview.molecule.validate()?;
                preview.repair.source.validate()?;
                preview.repair.duplex.validate()?;
                if preview.molecule.sequence != preview.repair.duplex.top {
                    return Err(Error::Format(
                        "repair preview sequence/duplex mismatch".into(),
                    ));
                }
            }
            AnalysisResult::AnnealPreview(preview) => {
                for candidate in &preview.annealing.candidates {
                    candidate.duplex.validate()?;
                }
            }
            AnalysisResult::AnnotatedDocument(molecule)
            | AnalysisResult::PhysicalPreviewDocument(molecule) => molecule.validate()?,
            AnalysisResult::WorkspaceLoaded(workspace) => workspace.validate()?,
            AnalysisResult::WorkspaceSaved(data) => {
                if data.len() > crate::workspace::MAX_ARCHIVE_BYTES {
                    return Err(Error::TooLarge);
                }
            }
            _ => {}
        }
        Ok(response)
    }
}
