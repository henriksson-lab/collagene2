//! Bounded stateless worker entrypoint, typed domain routing and current response serialization.
mod assembly;
mod context;
mod documents;
mod library;
mod material;
mod output;
mod preparation;
mod restriction;
mod review;
mod router;
mod sequence;
mod topology;
mod tracks;

use super::{
    envelope::{
        MAX_WORKSPACE_MESSAGE_BYTES, WORKSPACE_REQUEST_SELECTORS, decode_with_limit,
        encode_with_limit, worker_limit,
    },
    request::Request,
    response::Response,
    schema::AnalysisResult,
};
use crate::formats::sequence::MAX_IMPORT_BYTES;
use collagene_model::Error;

pub fn process_request(bytes: &[u8]) -> Result<Vec<u8>, Error> {
    let maximum = worker_limit(b"CGQ\0", bytes, WORKSPACE_REQUEST_SELECTORS)?;
    let Request { id, revision, job } = decode_with_limit(b"CGQ\0", bytes, maximum)?;
    let context = context::RequestContext {
        revision: &revision,
    };
    let result = router::dispatch(job, &context)?;
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
            | AnalysisResult::MaterialDigestProduct(_)
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
            id,
            revision,
            result,
        },
        maximum,
    )
}
