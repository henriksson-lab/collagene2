//! Response ownership and validation after bounded decoding.
use super::{
    bindings::profile_binding,
    envelope::{decode_with_limit, encode_with_limit, worker_limit},
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

const LARGE_REPLY_SELECTORS: &[u32] = &[
    20, 21, 40, 41, 62, 65, 92, 95, 98, 120, 124, 125, 157, 158, 159, 160, 161, 162, 163, 165, 166,
    168, 169, 171, 174, 176, 178, 179, 180, 182, 187, 188,
];
/// Bounded serialization of only the id/revision/selector prefix. Bulk fields
/// are never serialized into a temporary full-size buffer to choose a budget.
struct ReplyPrefix<'a>(&'a mut Vec<u8>);
impl postcard::ser_flavors::Flavor for ReplyPrefix<'_> {
    type Output = ();
    fn try_push(&mut self, byte: u8) -> Result<(), postcard::Error> {
        if self.0.len() == 256 {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.0.push(byte);
        Ok(())
    }
    fn try_extend(&mut self, bytes: &[u8]) -> Result<(), postcard::Error> {
        let count = bytes.len().min(256 - self.0.len());
        self.0.extend_from_slice(&bytes[..count]);
        if count < bytes.len() {
            Err(postcard::Error::SerializeBufferFull)
        } else {
            Ok(())
        }
    }
    fn finalize(self) -> Result<(), postcard::Error> {
        Ok(())
    }
}
impl Response {
    /// Resource-only admission from at most 256 borrowed prefix bytes and the
    /// separately observed whole frame length. Full checksum/typed validation
    /// still runs in decode after a client admits the complete buffer.
    pub fn admitted_frame_limit(prefix: &[u8], frame_length: usize) -> Result<usize, Error> {
        worker_limit(b"CGR\0", prefix, frame_length, LARGE_REPLY_SELECTORS)
    }
    pub fn encode(&self) -> Result<Vec<u8>, Error> {
        if self.revision.len() > 128 {
            return Err(Error::TooLarge);
        }
        let mut prefix = Vec::with_capacity(256);
        match postcard::serialize_with_flavor(self, ReplyPrefix(&mut prefix)) {
            Ok(()) | Err(postcard::Error::SerializeBufferFull) => {}
            Err(e) => return Err(Error::Format(e.to_string())),
        }
        let ((_, _, selector), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&prefix)
            .map_err(|e| Error::Format(e.to_string()))?;
        let maximum = if LARGE_REPLY_SELECTORS.contains(&selector) {
            super::envelope::MAX_WORKSPACE_MESSAGE_BYTES
        } else {
            crate::formats::sequence::MAX_IMPORT_BYTES
        };
        encode_with_limit(b"CGR\0", self, maximum)
    }
    pub fn decode(bytes: &[u8]) -> Result<Response, Error> {
        let maximum = Self::admitted_frame_limit(bytes, bytes.len())?;
        let response: Response = decode_with_limit(b"CGR\0", bytes, maximum)?;
        match &response.result {
            AnalysisResult::ImageEvidencePrepared(prepared) => {
                prepared.validate_preparation_binding(&response.revision)?;
            }
            AnalysisResult::ImageEvidencePreview(preview) => {
                preview.validate_transport_binding(&response.revision)?;
            }
            AnalysisResult::TraceDisplayPrepared(data) => {
                data.validate_binding(&response.revision)?
            }
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
            AnalysisResult::GelImageExported(report) => {
                report.validate_binding(&response.revision)?
            }
            AnalysisResult::MapExported(report) => {
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
            AnalysisResult::MaterialDigestProduct(report) => {
                report.validate()?;
                if report.context_revision != response.revision {
                    return Err(Error::Format(
                        "Digest export response context changed".into(),
                    ));
                }
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
                        .is_some_and(|b| b.len() > crate::material::evidence::MAX_BYTES)
                    || value.bytes.is_none() != (value.step_count == 0)
                {
                    return Err(Error::Format(
                        "Invalid historical scientific-evidence download".into(),
                    ));
                }
            }

            AnalysisResult::LigationProfiles(items) => {
                let expected = crate::cloning::ligation::profiles::list()?;
                if items.len() != expected.len()
                    || response.revision != crate::cloning::ligation::profiles::catalog_revision()?
                    || crate::encoding::bounded_json::hash(items)?
                        != crate::encoding::bounded_json::hash(&expected)?
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
                    Some(crate::cloning::assembly::junctions::Method::GoldenGate { .. })
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
                crate::cloning::topo::worker::validate_simulation(report, &response.revision)?
            }
            AnalysisResult::TopoOrientations(report) => {
                crate::cloning::topo::worker::validate_orientations(report, &response.revision)?
            }
            AnalysisResult::TopoProjection(report) => {
                crate::cloning::topo::worker::validate_projection(report, &response.revision)?
            }
            AnalysisResult::TopoTarget(report) => {
                crate::cloning::topo::worker::validate_target(report, &response.revision)?
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
                crate::tracks::genome_session::validate_result(
                    &response.result,
                    &response.revision,
                )?
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
                if rows.len() > crate::formats::alignment::MAX_ALIGNMENT_ROWS {
                    return Err(Error::TooLarge);
                }
                for row in rows {
                    row.validate()?;
                }
            }
            AnalysisResult::AlignmentSaved(archive) => {
                crate::formats::alignment::load_alignment(archive)?;

                if blake3::hash(archive).to_hex().as_str() != response.revision {
                    return Err(Error::Format(
                        "saved alignment response revision mismatch".into(),
                    ));
                }
            }
            AnalysisResult::AlignmentExport(report) => {
                if report.text.len() > crate::formats::alignment::MAX_ALIGNMENT_BYTES
                    || report.losses.len() > crate::formats::alignment::MAX_ALIGNMENT_ANNOTATIONS
                    || report.renamings.len() > crate::formats::alignment::MAX_ALIGNMENT_ROWS
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
                if data.len() > crate::documents::workspace::MAX_ARCHIVE_BYTES {
                    return Err(Error::TooLarge);
                }
            }
            _ => {}
        }
        Ok(response)
    }
}

#[cfg(test)]
mod prefix_tests {
    use super::*;
    use crate::formats::sequence::MAX_IMPORT_BYTES;
    use crate::protocol::{envelope::MAX_WORKSPACE_MESSAGE_BYTES, schema::WORKER_SCHEMA_VERSION};
    fn prefix(selector: u32, revision: &str, payload_len: usize) -> (Vec<u8>, usize) {
        let body = postcard::to_allocvec(&(u32::MAX, revision, selector)).unwrap();
        let length = payload_len.max(body.len());
        let mut prefix = b"CGR\0".to_vec();
        prefix.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
        prefix.extend_from_slice(&(length as u32).to_le_bytes());
        prefix.extend_from_slice(&[0; 32]);
        prefix.extend(body);
        (prefix, 42 + length)
    }
    #[test]
    fn resource_prefix_admission_retains_ordinary_and_explicit_existing_large_selectors() {
        let (ordinary, length) = prefix(0, &"a".repeat(64), 100);
        assert_eq!(
            Response::admitted_frame_limit(&ordinary, length).unwrap(),
            MAX_IMPORT_BYTES
        );
        // Workspace, combinatorial and trace display replies already declare
        // their larger profiles. An arbitrary selector cannot acquire that cap.
        for selector in [20, 124, 188] {
            let (prefix, length) = prefix(selector, &"a".repeat(64), MAX_IMPORT_BYTES + 1);
            assert_eq!(
                Response::admitted_frame_limit(&prefix, length).unwrap(),
                MAX_WORKSPACE_MESSAGE_BYTES
            );
        }
        for selector in [0, 999] {
            let (prefix, length) = prefix(selector, &"a".repeat(64), MAX_IMPORT_BYTES + 1);
            assert!(matches!(
                Response::admitted_frame_limit(&prefix, length),
                Err(Error::TooLarge)
            ));
        }
    }
    #[test]
    fn prefix_rejects_wrong_namespace_version_length_and_oversized_claims_before_bulk_copy() {
        let (good, length) = prefix(188, &"a".repeat(64), MAX_IMPORT_BYTES + 1);
        let mut wrong = good.clone();
        wrong[0] = b'X';
        assert!(Response::admitted_frame_limit(&wrong, length).is_err());
        let mut wrong = good.clone();
        wrong[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
        assert!(Response::admitted_frame_limit(&wrong, length).is_err());
        assert!(Response::admitted_frame_limit(&good, length + 1).is_err());
        assert!(Response::admitted_frame_limit(&good[..42], length).is_err());
        let (too_large, length) = prefix(188, "a", MAX_WORKSPACE_MESSAGE_BYTES + 1);
        assert!(matches!(
            Response::admitted_frame_limit(&too_large, length),
            Err(Error::TooLarge)
        ));
    }
    #[test]
    fn current_revision_boundary_is_explicit_and_prefix_admission_never_skips_body_validation() {
        let response = Response {
            id: u32::MAX,
            revision: "a".repeat(128),
            result: AnalysisResult::Error("bounded diagnostic".into()),
        };
        let bytes = response.encode().unwrap();
        assert_eq!(
            Response::admitted_frame_limit(&bytes[..bytes.len().min(256)], bytes.len()).unwrap(),
            MAX_IMPORT_BYTES
        );
        assert_eq!(
            Response::decode(&bytes).unwrap().revision,
            response.revision
        );
        let mut corrupt = bytes;
        let end = corrupt.len() - 1;
        corrupt[end] ^= 1;
        assert!(
            Response::admitted_frame_limit(&corrupt[..corrupt.len().min(256)], corrupt.len())
                .is_ok()
        );
        assert!(Response::decode(&corrupt).is_err());
        let oversized = Response {
            revision: "a".repeat(129),
            ..response
        };
        assert!(matches!(oversized.encode(), Err(Error::TooLarge)));
        let (prefix, length) = prefix(0, &"a".repeat(129), 1000);
        assert!(matches!(
            Response::admitted_frame_limit(&prefix, length),
            Err(Error::TooLarge)
        ));
    }
}
