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

pub(super) fn vector_profiles(context: &RequestContext) -> Result<AnalysisResult, Error> {
    bound(context, construct_qc::vector_class::catalogue_revision()?)?;
    Ok(AnalysisResult::VectorClassProfiles(
        construct_qc::vector_class::catalogue()?,
    ))
}
pub(super) fn vector_class(
    context: &RequestContext,
    molecule: Molecule,
    request: construct_qc::vector_class::Request,
) -> Result<AnalysisResult, Error> {
    bound(
        context,
        construct_qc::vector_class::revision(&molecule, &request)?,
    )?;
    Ok(construct_qc::vector_class::assess(&molecule, &request)
        .map(AnalysisResult::VectorClass)
        .unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

/// The scientific owner resolves native events and replays the selected history.
pub(super) fn causal_qc(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: construct_qc::causal::Request,
) -> Result<AnalysisResult, Error> {
    bound(
        context,
        construct_qc::causal::revision(&workspace, &request)?,
    )?;
    let result = construct_qc::causal::assess(&workspace, &request).and_then(|report| {
        report.validate_request(&workspace, &request)?;
        Ok(report)
    });
    Ok(result
        .map(AnalysisResult::ConstructCausalQc)
        .unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

/// Model only explicitly named products/mixtures with the existing scientific owners.
pub(super) fn screening(
    context: &RequestContext,
    request: crate::cloning::screening::Request,
) -> Result<AnalysisResult, Error> {
    bound(context, crate::cloning::screening::revision(&request)?)?;
    Ok(crate::cloning::screening::assess(&request)
        .map(AnalysisResult::ScreeningStrategy)
        .unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
