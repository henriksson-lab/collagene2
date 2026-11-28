//! Source-bound domain admission; no second scientific implementation.
use super::context::RequestContext;
use crate::{
    cloning::{assembly::standards, pcr::collection_binding},
    construct_qc,
    protocol::schema::AnalysisResult,
};
use collagene_model::{Error, Molecule};
fn bound(context: &RequestContext, revision: String) -> Result<(), Error> {
    if revision != *context.revision {
        return Err(Error::Format("Design-check request context changed".into()));
    }
    context.admit(None)
}
pub(super) fn construct_qc(
    context: &RequestContext,
    molecule: Molecule,
    options: construct_qc::Options,
) -> Result<AnalysisResult, Error> {
    bound(context, construct_qc::revision(&molecule, &options)?)?;
    Ok(construct_qc::assess(&molecule, &options)
        .map(AnalysisResult::ConstructQc)
        .unwrap_or_else(|e| AnalysisResult::Error(e.to_string())))
}
pub(super) fn assembly_standard(
    context: &RequestContext,
    request: standards::Request,
) -> Result<AnalysisResult, Error> {
    bound(context, standards::revision(&request)?)?;
    Ok(standards::assess(&request)
        .map(AnalysisResult::AssemblyStandard)
        .unwrap_or_else(|e| AnalysisResult::Error(e.to_string())))
}
pub(super) fn oligo_collection(
    context: &RequestContext,
    request: collection_binding::Request,
) -> Result<AnalysisResult, Error> {
    bound(context, collection_binding::revision(&request)?)?;
    Ok(collection_binding::assess(&request)
        .map(AnalysisResult::OligoCollectionBinding)
        .unwrap_or_else(|e| AnalysisResult::Error(e.to_string())))
}
