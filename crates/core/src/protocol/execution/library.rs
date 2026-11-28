//! Combinatorial design, release and required-session admission.
use super::context::RequestContext;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Error, Molecule};

pub(super) fn export_combinatorial_outputs(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    manifest_hash: String,
    options: crate::planning::combinatorial::export::Options,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let manifest_hash = &manifest_hash;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::planning::combinatorial::export::revision(workspace, manifest_hash, options)?
                != (*context.revision)
            {
                return Err(Error::Format("Combinatorial export inputs changed.".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::planning::combinatorial::export::export(&workspace, &manifest_hash, &options)
            .map(AnalysisResult::CombinatorialOutputsArchive);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_combinatorial_multi_site(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    design: crate::planning::combinatorial::Design,
    construct: crate::planning::combinatorial::SelectedConstruct,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let design = &design;
        let construct = &construct;

        let molecule: Option<&Molecule> = {
            if crate::planning::combinatorial::inspect_multisite_revision(
                design, construct, workspace,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "combinatorial multisite inspection source/declaration context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::planning::combinatorial::inspect_multisite(&design, &construct, &workspace)
            .map(AnalysisResult::CombinatorialMultiSiteInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
