//! Source/context admission for the single native IVT check owner.
use super::context::RequestContext;
use crate::{protocol::schema::AnalysisResult, rna_chemistry::ivt};
use collagene_model::{Error, Molecule};

pub(super) fn profiles(context: &RequestContext) -> Result<AnalysisResult, Error> {
    if ivt::catalogue_revision().map_err(|error| Error::Format(error.to_string()))?
        != *context.revision
    {
        return Err(Error::Format(
            "IVT profile catalogue context changed".into(),
        ));
    }
    context.admit(None)?;
    ivt::profiles()
        .map(AnalysisResult::IvtProfiles)
        .map_err(|error| Error::Format(error.to_string()))
}

pub(super) fn assess(
    context: &RequestContext,
    molecule: Molecule,
    request: ivt::Request,
) -> Result<AnalysisResult, Error> {
    if ivt::revision(&molecule, &request).map_err(|error| Error::Format(error.to_string()))?
        != *context.revision
    {
        return Err(Error::Format(
            "IVT source/declaration context changed".into(),
        ));
    }
    context.admit(None)?;
    Ok(ivt::assess(&molecule, &request)
        .map(AnalysisResult::IvtTemplate)
        .unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
