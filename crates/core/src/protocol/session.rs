//! Persistent worker session and the state it exclusively owns.
use super::{
    envelope::{HEADER, WORKSPACE_REQUEST_SELECTORS, decode_with_limit, worker_limit},
    execution::process_request,
    request::Request,
    response::Response,
    schema::{AnalysisJob, AnalysisResult},
};
use collagene_model::Error;
/// A persistent worker retains a validated immutable track source/index. Replacing
/// it is atomic; viewport jobs carry only options/window and the expected revision.
#[derive(Default)]
pub struct AnalysisSession {
    tracks: Option<crate::tracks::view::PreparedTracks>,
    genomes: crate::tracks::genome_session::GenomeSession,
    read_alignments: crate::review::reads::session::Session,
    variants: crate::review::variants::session::Session,
    combinatorial: Option<crate::planning::combinatorial::batch::Session>,
    release: Option<crate::planning::combinatorial::release::ReleaseSession>,
}
impl AnalysisSession {
    pub fn process_request(&mut self, bytes: &[u8]) -> Result<Vec<u8>, Error> {
        let maximum = worker_limit(b"CGQ\0", bytes, bytes.len(), WORKSPACE_REQUEST_SELECTORS)?;
        // Route by a borrowed prefix so ordinary jobs are not fully decoded twice.
        // Stable append-only track, genome and read-alignment variants retain
        // private worker state; all other jobs use ordinary stateless dispatch.
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
            .map_err(|e| Error::Format(e.to_string()))?;
        if !matches!(tag, 32 | 33 | 53..=66 | 77..=81 | 86..=91 | 126..=132) {
            return process_request(bytes);
        }

        let r: Request = decode_with_limit(b"CGQ\0", bytes, maximum)?;
        let result = match r.job {
            AnalysisJob::BeginCombinatorialRelease { workspace, request } => {
                let pending = crate::planning::combinatorial::release::ReleaseSession::begin(
                    workspace, request,
                )?;
                let progress = pending.summary();
                if progress.context_revision != r.revision {
                    return Err(Error::Format(
                        "Combinatorial release source or options changed.".into(),
                    ));
                }
                self.release = Some(pending);
                AnalysisResult::CombinatorialReleaseStarted(progress)
            }
            AnalysisJob::AdvanceCombinatorialRelease { context, index } => {
                if context != r.revision {
                    return Err(Error::Format(
                        "Combinatorial release advance context changed.".into(),
                    ));
                }
                let session = self.release.as_mut().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial release has been started.".into())
                })?;
                let progress = match session.advance(&context, index) {
                    Ok(progress) => progress,
                    Err(error) => {
                        // A failed reconstruction poisons the complete stage.
                        // A rejected out-of-order request has not executed a row.
                        if session.summary().failed.is_some() {
                            self.release = None;
                        }
                        return Err(error);
                    }
                };
                if progress.failed.is_some() {
                    self.release = None;
                }
                AnalysisResult::CombinatorialReleaseAdvanced(progress)
            }
            AnalysisJob::InspectCombinatorialRelease { context, index } => {
                if context != r.revision {
                    return Err(Error::Format(
                        "Combinatorial release inspection context changed.".into(),
                    ));
                }
                let session = self.release.as_ref().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial release has been started.".into())
                })?;
                if session.summary().context_revision != context {
                    return Err(Error::Format(
                        "Combinatorial release inspection session changed.".into(),
                    ));
                }
                AnalysisResult::CombinatorialReleaseInspected(session.inspect(index)?)
            }
            AnalysisJob::FinishCombinatorialRelease {
                workspace,
                current_context,
            } => {
                if current_context != r.revision
                    || self
                        .release
                        .as_ref()
                        .is_none_or(|session| session.summary().context_revision != current_context)
                {
                    return Err(Error::Format(
                        "Combinatorial release finish context changed.".into(),
                    ));
                }
                // Consuming the private stage makes every failed finish atomic:
                // no decoded prefix or subsequent retry can publish its products.
                let session = self.release.take().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial release has been started.".into())
                })?;
                AnalysisResult::CombinatorialReleaseFinished(
                    session.finish_to(&workspace, &current_context)?,
                )
            }
            AnalysisJob::BeginCombinatorialBatch { workspace, matrix } => {
                let pending =
                    crate::planning::combinatorial::batch::Session::begin(workspace, matrix)?;
                let progress = pending.summary();
                if progress.context_revision != r.revision {
                    return Err(Error::Format(
                        "Combinatorial batch source or matrix changed.".into(),
                    ));
                }
                self.combinatorial = Some(pending);
                AnalysisResult::CombinatorialBatchStarted(progress)
            }
            AnalysisJob::AdvanceCombinatorialBatch { context, index } => {
                if context != r.revision {
                    return Err(Error::Format(
                        "Combinatorial batch advance context changed.".into(),
                    ));
                }
                let session = self.combinatorial.as_mut().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial batch has been started.".into())
                })?;
                AnalysisResult::CombinatorialBatchAdvanced(session.advance(&context, index)?)
            }
            AnalysisJob::CaptureCombinatorialBatch {
                workspace,
                current_context,
            } => {
                if crate::planning::combinatorial::batch::capture_revision(
                    &workspace,
                    &current_context,
                )? != r.revision
                {
                    return Err(Error::Format(
                        "Combinatorial batch capture workspace changed.".into(),
                    ));
                }
                let session = self.combinatorial.as_ref().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial batch has been started.".into())
                })?;
                AnalysisResult::CombinatorialBatchCaptured(
                    session.capture_to(&workspace, &current_context)?,
                )
            }
            job @ (AnalysisJob::PrepareVariants { .. }
            | AnalysisJob::VariantPage { .. }
            | AnalysisJob::CheckVariantReference { .. }
            | AnalysisJob::PreviewVariantAllele { .. }
            | AnalysisJob::ExportVariants { .. }
            | AnalysisJob::DisposeVariants { .. }) => self.variants.process(job, &r.revision)?,
            job @ (AnalysisJob::PrepareReadAlignments { .. }
            | AnalysisJob::ReadAlignmentPage { .. }
            | AnalysisJob::CheckReadAlignment { .. }
            | AnalysisJob::ExportReadAlignments { .. }
            | AnalysisJob::DisposeReadAlignments { .. }) => {
                self.read_alignments.process(job, &r.revision)?
            }
            job @ (AnalysisJob::BeginGenomeFasta { .. }
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
            | AnalysisJob::ReadGenomeInterval { .. }) => self.genomes.process(job, &r.revision)?,
            AnalysisJob::PrepareTracks { molecule } => {
                let prepared =
                    crate::tracks::view::PreparedTracks::new(std::sync::Arc::new(molecule))?;
                if prepared.source_revision() != r.revision {
                    return Err(Error::Format("track preparation revision mismatch".into()));
                }
                self.tracks = Some(prepared);
                AnalysisResult::TracksPrepared(r.revision.clone())
            }
            AnalysisJob::CachedViewTracks { window, options } => {
                let prepared = self.tracks.as_ref().ok_or_else(|| {
                    Error::Unsupported("cached track source has not been prepared".into())
                })?;
                AnalysisResult::ViewTracks(prepared.report_for_revision(
                    &r.revision,
                    window,
                    &options,
                )?)
            }
            _ => return process_request(bytes),
        };
        Response {
            id: r.id,
            revision: r.revision,
            result,
        }
        .encode()
    }
}
