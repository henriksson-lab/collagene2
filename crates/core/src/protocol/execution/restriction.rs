//! Restriction site inspection and fragment materialization handlers.
use super::context::RequestContext;
use crate::protocol::schema::AnalysisResult;
use crate::sequence::analysis::Enzyme;
use collagene_model::{Error, Molecule};

pub(super) fn assess_conditions(
    context: &RequestContext,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    options: crate::cloning::restriction::conditions::ReactionConditions,
) -> Result<AnalysisResult, Error> {
    use crate::cloning::restriction::{catalog, conditions};
    let installed = catalog::builtin_catalog()?;
    if conditions::revision(&installed, &selection, &options)? != *context.revision {
        return Err(Error::Format("Restriction reaction context changed".into()));
    }
    context.admit(None)?;
    Ok(conditions::assess(&installed, &selection, &options)
        .map(AnalysisResult::RestrictionConditions)
        .unwrap_or_else(|e| AnalysisResult::Error(e.to_string())))
}

pub(super) fn inspect_restriction(
    context: &RequestContext,
    source_role: String,
    molecule: Molecule,
    enzymes: Vec<Enzyme>,
    options: crate::cloning::restriction::inspection::InspectionOptions,
) -> Result<AnalysisResult, Error> {
    {
        let source_role = &source_role;
        let molecule = &molecule;
        let enzymes = &enzymes;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::cloning::restriction::inspection::batch_revision(
                source_role,
                molecule,
                enzymes,
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "restriction inspection source/role/geometry/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::restriction::inspection::inspect(
        &source_role,
        &molecule,
        &enzymes,
        &options,
    )
    .map(AnalysisResult::RestrictionInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn materialize_digest_fragment(
    context: &RequestContext,
    source_role: String,
    molecule: Molecule,
    enzymes: Vec<Enzyme>,
    options: crate::cloning::restriction::inspection::InspectionOptions,
    fragment_index: usize,
) -> Result<AnalysisResult, Error> {
    {
        let source_role = &source_role;
        let molecule = &molecule;
        let enzymes = &enzymes;
        let options = &options;
        let fragment_index = &fragment_index;

        let molecule: Option<&Molecule> = {
            if crate::cloning::restriction::inspection::product_revision(
                source_role,
                molecule,
                enzymes,
                options,
                *fragment_index,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "digest product source/role/geometry/options/index changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::restriction::inspection::product(
        &source_role,
        &molecule,
        &enzymes,
        &options,
        fragment_index,
    )
    .map(AnalysisResult::DigestProduct);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn install_digest_fragment(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    source: crate::material::lineage::EventBinding,
    enzymes: Vec<Enzyme>,
    options: crate::cloning::restriction::inspection::InspectionOptions,
    fragment_index: usize,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let source = &source;
        let enzymes = &enzymes;
        let options = &options;
        let fragment_index = &fragment_index;

        let molecule: Option<&Molecule> = {
            if crate::cloning::restriction::history::installation_revision(
                workspace,
                source,
                enzymes,
                options,
                *fragment_index,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "digest installation workspace/source/geometry/options/index changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::restriction::history::install(
        &workspace,
        &source,
        &enzymes,
        &options,
        fragment_index,
    )
    .map(AnalysisResult::DigestInstalled);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_catalog_restriction(
    context: &RequestContext,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    source_role: String,
    molecule: Molecule,
    options: crate::cloning::restriction::inspection::InspectionOptions,
) -> Result<AnalysisResult, Error> {
    {
        let selection = &selection;
        let source_role = &source_role;
        let molecule = &molecule;
        let options = &options;

        let molecule: Option<&Molecule> = {
            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
            if crate::cloning::restriction::catalog::batch_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "catalog restriction source/pin/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
        crate::cloning::restriction::catalog::inspect(
            &catalog,
            &selection,
            &source_role,
            &molecule,
            &options,
        )
        .map(AnalysisResult::CatalogRestrictionInspection)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn materialize_catalog_digest_fragment(
    context: &RequestContext,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    source_role: String,
    molecule: Molecule,
    options: crate::cloning::restriction::inspection::InspectionOptions,
    fragment_index: usize,
) -> Result<AnalysisResult, Error> {
    {
        let selection = &selection;
        let source_role = &source_role;
        let molecule = &molecule;
        let options = &options;
        let fragment_index = &fragment_index;

        let molecule: Option<&Molecule> = {
            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
            if crate::cloning::restriction::catalog::product_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
                *fragment_index,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "catalog product source/pin/options/index changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
        crate::cloning::restriction::catalog::product(
            &catalog,
            &selection,
            &source_role,
            &molecule,
            &options,
            fragment_index,
        )
        .map(AnalysisResult::CatalogDigestProduct)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn install_catalog_digest_fragment(
    context: &RequestContext,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    workspace: crate::documents::workspace::Workspace,
    source: crate::material::lineage::EventBinding,
    options: crate::cloning::restriction::inspection::InspectionOptions,
    fragment_index: usize,
) -> Result<AnalysisResult, Error> {
    {
        let selection = &selection;
        let workspace = &workspace;
        let source = &source;
        let options = &options;
        let fragment_index = &fragment_index;

        let molecule: Option<&Molecule> = {
            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
            if crate::cloning::restriction::catalog::installation_revision(
                &catalog,
                selection,
                workspace,
                source,
                options,
                *fragment_index,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "catalog installation workspace/source/pin/options/index changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
        crate::cloning::restriction::catalog::install(
            &catalog,
            &selection,
            &workspace,
            &source,
            &options,
            fragment_index,
        )
        .map(AnalysisResult::CatalogDigestInstalled)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn scan_restriction(
    context: &RequestContext,
    source_role: String,
    molecule: Molecule,
    enzymes: Vec<Enzyme>,
    options: crate::cloning::restriction::inspection::InspectionOptions,
) -> Result<AnalysisResult, Error> {
    {
        let source_role = &source_role;
        let molecule = &molecule;
        let enzymes = &enzymes;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::cloning::restriction::inspection::batch_revision(
                source_role,
                molecule,
                enzymes,
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "restriction inspection source/role/geometry/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::restriction::inspection::scan(&source_role, &molecule, &enzymes, &options)
            .map(AnalysisResult::RestrictionInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn select_digest_fragment(
    context: &RequestContext,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    source_role: String,
    molecule: Molecule,
    options: crate::cloning::restriction::inspection::InspectionOptions,
    expected_batch_revision: String,
    choice: crate::cloning::restriction::cut_selection::Choice,
) -> Result<AnalysisResult, Error> {
    {
        let selection = &selection;
        let source_role = &source_role;
        let molecule = &molecule;
        let options = &options;
        let expected_batch_revision = &expected_batch_revision;
        let choice = &choice;

        let molecule: Option<&Molecule> = {
            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
            if crate::cloning::restriction::cut_selection::selection_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
                expected_batch_revision,
                choice,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Digest cut selection source/boundaries changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
        crate::cloning::restriction::cut_selection::preview(
            &catalog,
            &selection,
            &source_role,
            &molecule,
            &options,
            &expected_batch_revision,
            &choice,
        )
        .map(AnalysisResult::DigestFragmentSelected)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_run_catalog_restriction(
    context: &RequestContext,
    source_role: String,
    molecule: Molecule,
    origin: crate::material::capsule::Capsule,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    options: crate::cloning::restriction::inspection::InspectionOptions,
) -> Result<AnalysisResult, Error> {
    {
        let source_role = &source_role;
        let molecule = &molecule;
        let origin = &origin;
        let selection = &selection;
        let options = &options;

        let molecule: Option<&Molecule> = {
            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
            if crate::cloning::restriction::run_digest::inspection_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                origin,
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Saved-run catalog inspection origin/source/pin/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
        crate::cloning::restriction::run_digest::inspect(
            &catalog,
            &selection,
            &source_role,
            &molecule,
            &origin,
            &options,
        )
        .map(AnalysisResult::RunCatalogRestrictionInspection)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
