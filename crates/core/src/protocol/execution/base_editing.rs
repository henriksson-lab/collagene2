//! Exact source and declared finite-editor-context admission.
use super::context::RequestContext;
use crate::{protocol::schema::AnalysisResult, sequence::crispr::base_editing};
use collagene_model::{Error, Molecule};

pub(super) fn assess(
    context: &RequestContext,
    molecule: Molecule,
    request: base_editing::Request,
) -> Result<AnalysisResult, Error> {
    if base_editing::revision(&molecule, &request)? != *context.revision {
        return Err(Error::Format(
            "Base-editing source/declaration context changed".into(),
        ));
    }
    context.admit(None)?;
    Ok(base_editing::assess(&molecule, &request)
        .map(AnalysisResult::BaseEditing)
        .unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
