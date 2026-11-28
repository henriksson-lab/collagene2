//! Verified protocol output, predicted gels and presentation artifacts.
use super::context::RequestContext;
use crate::formats::sequence::MAX_IMPORT_BYTES;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Error, Molecule};

pub(super) fn verify(
    context: &RequestContext,
    target: Molecule,
    plan: crate::planning::plan::DesignPlan,
) -> Result<AnalysisResult, Error> {
    {
        let target = &target;

        let molecule: Option<&Molecule> = Some(target);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Ok(AnalysisResult::Verification(
        crate::verification::engine::verify(&plan, &target),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_protocol_pdf(
    context: &RequestContext,
    target: Molecule,
    plan: crate::planning::plan::DesignPlan,
) -> Result<AnalysisResult, Error> {
    {
        let target = &target;

        let molecule: Option<&Molecule> = Some(target);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::output::pdf::export_protocol_pdf(&plan, &target).map(AnalysisResult::ProtocolPdf)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn predict_gel(
    context: &RequestContext,
    molecule: Molecule,
    lanes: Vec<crate::output::gel::LaneInput>,
    calibration: crate::output::gel::Calibration,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let lanes = &lanes;
        if lanes.iter().any(|lane| {
            matches!(
                lane.preparation,
                crate::output::gel::Preparation::CatalogDigest { .. }
            )
        }) {
            return Err(Error::Format(
                "Catalog-pinned lanes require the pinned-gel worker operation.".into(),
            ));
        }
        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::output::gel::predict_gel(&lanes, &calibration).map(AnalysisResult::Gel);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn diagnostic_digests(
    context: &RequestContext,
    molecule: Molecule,
    candidates: Vec<Molecule>,
    options: crate::output::gel::DiagnosticOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::output::gel::diagnostic_digests(&candidates, &options)
            .map(AnalysisResult::DiagnosticDigest);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_plan_oligos(
    context: &RequestContext,
    plan: crate::planning::plan::DesignPlan,
    target: Molecule,
    options: crate::planning::oligos::Options,
    delimiter: crate::formats::oligos::Delimiter,
) -> Result<AnalysisResult, Error> {
    {
        let plan = &plan;
        let target = &target;
        let options = &options;
        let delimiter = &delimiter;

        let molecule: Option<&Molecule> = {
            if crate::planning::oligos::export_revision(plan, target, options, *delimiter)?
                != (*context.revision)
            {
                return Err(Error::Format("Plan oligo export context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::planning::oligos::export(&plan, &target, &options, delimiter)
            .map(AnalysisResult::PlanOligosExported);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_protocol_html(
    context: &RequestContext,
    plan: crate::planning::plan::DesignPlan,
    target: Molecule,
) -> Result<AnalysisResult, Error> {
    {
        let target = &target;

        let molecule: Option<&Molecule> = Some(target);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Ok(AnalysisResult::ProtocolHtml(
        crate::verification::output::export_protocol_html(&plan, &target),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn predict_pinned_gel(
    context: &RequestContext,
    molecule: Molecule,
    lanes: Vec<crate::output::gel::LaneInput>,
    calibration: crate::output::gel::Calibration,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let lanes = &lanes;
        let calibration = &calibration;

        let molecule: Option<&Molecule> = {
            if lanes.is_empty() || lanes.len() > 32 {
                return Err(Error::TooLarge);
            }
            calibration.validate()?;
            molecule.validate()?;
            let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
            let mut bases = 0usize;
            for lane in lanes {
                lane.molecule.validate()?;
                bases = bases
                    .checked_add(lane.molecule.sequence.len())
                    .filter(|n| *n <= MAX_IMPORT_BYTES)
                    .ok_or(Error::TooLarge)?;
                if let crate::output::gel::Preparation::CatalogDigest { recipe } = &lane.preparation
                {
                    collagene_model::bounds::model_preflight(&[&lane.molecule])?;
                    recipe.validate()?;
                    recipe.selection.validate_for(&catalog)?;
                    if lane.molecule.revision() != recipe.source.revision {
                        return Err(Error::Format("Pinned gel source revision changed.".into()));
                    }
                    recipe.inspection_revision(&lane.molecule)?;
                }
            }
            Some(molecule)
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::output::gel::predict_gel(&lanes, &calibration).map(AnalysisResult::PinnedGel);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn render_map(
    context: &RequestContext,
    request: crate::output::map_export::Request,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;

        let molecule: Option<&Molecule> = {
            request.validate()?;
            if (*context.revision).len() != 64
                || !context
                    .revision
                    .bytes()
                    .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
            {
                return Err(Error::Format(
                    "Map export capture context is invalid".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { request.render().map(AnalysisResult::MapExported) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_protocol_layout(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::output::layout::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::output::layout::preview_revision(workspace, request)? != (*context.revision) {
                return Err(Error::Format(
                    "Protocol layout source or assignment context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::output::layout::preview(&workspace, &request)
            .map(AnalysisResult::ProtocolLayoutPreview)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn capture_protocol_layout(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::output::layout::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::output::layout::preview_revision(workspace, request)? != (*context.revision) {
                return Err(Error::Format(
                    "Protocol layout source or assignment context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::output::layout::capture(&workspace, &request)
            .map(AnalysisResult::ProtocolLayoutCaptured)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_protocol_layout(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    layout_hash: String,
    format: crate::output::layout::ExportFormat,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let layout_hash = &layout_hash;
        let format = &format;

        let molecule: Option<&Molecule> = {
            if crate::output::layout::export_revision(workspace, layout_hash, *format)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Protocol layout export source or format changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::output::layout::export(&workspace, &layout_hash, format)
            .map(AnalysisResult::ProtocolLayoutExported);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_protocol_quantities(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::output::quantities::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::output::quantities::preview_revision(workspace, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Reaction quantity source/recipe context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::output::quantities::preview(&workspace, &request)
            .map(AnalysisResult::ProtocolQuantitiesPreview)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn capture_protocol_quantities(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::output::quantities::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::output::quantities::preview_revision(workspace, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Reaction quantity source/recipe context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::output::quantities::capture(&workspace, &request)
            .map(AnalysisResult::ProtocolQuantitiesCaptured)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn open_protocol_quantities(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    quantity_hash: String,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let quantity_hash = &quantity_hash;

        let molecule: Option<&Molecule> = {
            if crate::output::quantities::open_revision(workspace, quantity_hash)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Reaction quantity attachment/workspace changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::output::quantities::open(&workspace, &quantity_hash)
            .map(AnalysisResult::ProtocolQuantitiesOpened);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_protocol_quantities(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    quantity_hash: String,
    format: crate::output::quantities::ExportFormat,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let quantity_hash = &quantity_hash;
        let format = &format;

        let molecule: Option<&Molecule> = {
            if crate::output::quantities::export_revision(workspace, quantity_hash, *format)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Reaction quantity export source/format changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::output::quantities::export(&workspace, &quantity_hash, format)
            .map(AnalysisResult::ProtocolQuantitiesExported);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn list_protocol_quantities(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::output::quantities::list_revision(workspace)? != (*context.revision) {
                return Err(Error::Format(
                    "Reaction quantity catalog/workspace changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::output::quantities::list(&workspace).map(AnalysisResult::ProtocolQuantitiesListed)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_structure(
    context: &RequestContext,
    bytes: Vec<u8>,
    format: crate::protein::structure::coordinates::Format,
    selection: Option<crate::protein::structure::coordinates::scene::Selection>,
) -> Result<AnalysisResult, Error> {
    {
        let bytes = &bytes;
        let format = &format;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::protein::structure::revision(bytes, *format, selection)?
                != (*context.revision)
            {
                return Err(Error::Format("Structure source/options changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::protein::structure::prepare(&bytes, format, selection)
            .map(AnalysisResult::StructurePrepared);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

/// Documentary rendering is bound to the full original report and explicit image settings.
pub(super) fn render_gel_image(
    context: &RequestContext,
    request: crate::output::gel_export::Request,
) -> Result<AnalysisResult, Error> {
    if request.revision()? != *context.revision {
        return Err(Error::Format("Gel export context changed".into()));
    }
    context.admit(None)?;
    Ok(request
        .render()
        .map(AnalysisResult::GelImageExported)
        .unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_image_evidence(
    context: &RequestContext,
    name: String,
    association: crate::documents::image_evidence::Binding,
    original: Vec<u8>,
) -> Result<AnalysisResult, Error> {
    use crate::documents::image_evidence::{MAX_SOURCE_BYTES, Prepared};
    if original.is_empty() || original.len() > MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    let hash = blake3::hash(&original).to_hex().to_string();
    if association.preparation_revision(&name, &hash)? != *context.revision {
        return Err(Error::Format("Image preparation context changed".into()));
    }
    context.admit(None)?;
    Prepared::new(&name, association, original).map(AnalysisResult::ImageEvidencePrepared)
}
pub(super) fn reopen_image_evidence(
    context: &RequestContext,
    prepared: crate::documents::image_evidence::Prepared,
) -> Result<AnalysisResult, Error> {
    if prepared.summary()?.reopen_revision()? != *context.revision {
        return Err(Error::Format("Image reopening context changed".into()));
    }
    context.admit(None)?;
    prepared
        .reopen()
        .map(|opened| AnalysisResult::ImageEvidencePreview(opened.into_preview()))
}
