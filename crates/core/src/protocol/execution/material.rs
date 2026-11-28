//! Exact native material ancestry, reconstruction and saved-run handlers.
use super::context::RequestContext;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Error, Molecule};

pub(super) fn record_cloning_run(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    plan: crate::planning::plan::DesignPlan,
    target_binding: crate::material::lineage::EventBinding,
    input_bindings: Vec<crate::material::lineage::InputEventBinding>,
    allow_incomplete: bool,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::material::lineage::workspace_fingerprint(workspace)? != (*context.revision) {
                return Err(Error::Format("run request source workspace changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::material::lineage::record_run(
        &workspace,
        &plan,
        target_binding,
        input_bindings,
        allow_incomplete,
    )
    .map(AnalysisResult::RunCommit);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn replay_cloning_run(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    run_id: String,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::material::lineage::workspace_fingerprint(workspace)? != (*context.revision) {
                return Err(Error::Format("run request source workspace changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::material::lineage::replay_run(&workspace, &run_id).map(AnalysisResult::RunReplay)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn finalize_cloning_outputs(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    plan: crate::planning::plan::DesignPlan,
    target: crate::material::lineage::EventBinding,
    inputs: Vec<crate::material::lineage::InputEventBinding>,
    options: crate::material::outputs::Options,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let plan = &plan;
        let target = &target;
        let inputs = &inputs;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::material::outputs::finalize_revision(
                workspace,
                plan,
                target.clone(),
                inputs.clone(),
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Cloning output source/options context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::outputs::finalize(&workspace, &plan, target, inputs, &options)
            .map(AnalysisResult::CloningOutputsFinalized);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn load_run_evidence(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    run_id: String,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let run_id = &run_id;
        if run_id.len() != 64 {
            return Err(Error::Format("Invalid run identity".into()));
        }
        let molecule: Option<&Molecule> = {
            if crate::material::lineage::workspace_fingerprint(workspace)? != (*context.revision) {
                return Err(Error::Format("run request source workspace changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::material::evidence::download(&workspace, &run_id)
            .map(AnalysisResult::HistoricalRunEvidence)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_run_material(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    event: crate::material::lineage::EventBinding,
    selection: Option<crate::material::run::Selection>,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let event = &event;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::material::run::inspection_revision(workspace, event, selection)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Run material workspace/event/origin selection changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::run::inspect(&workspace, &event, &selection)
            .map(AnalysisResult::RunMaterialInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_origin(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    event: crate::material::lineage::EventBinding,
    selection: Option<crate::material::run::MaterialSelection>,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let event = &event;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::material::run::material_inspection_revision(workspace, event, selection)?
                != (*context.revision)
            {
                return Err(Error::Format("Material source/selection changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::run::inspect_material(&workspace, &event, &selection)
            .map(AnalysisResult::MaterialOriginInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn material_pcr_documents(
    context: &RequestContext,
    molecule: Molecule,
    origin: crate::material::run::MaterialOrigin,
    forward: crate::sequence::analysis::Primer,
    reverse: crate::sequence::analysis::Primer,
    template_strand: crate::chemistry::Strand,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let origin = &origin;
        let forward = &forward;
        let reverse = &reverse;
        let template_strand = &template_strand;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::material::revision(
                molecule,
                origin,
                forward,
                reverse,
                *template_strand,
            )? != (*context.revision)
            {
                return Err(Error::Format("Material PCR source/primers changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::pcr::material::simulate(
        &molecule,
        &origin,
        &forward,
        &reverse,
        template_strand,
    )
    .map(AnalysisResult::MaterialPcrDocuments);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_material_pcr_outputs(
    context: &RequestContext,
    origin: crate::material::run::MaterialOrigin,
    request: crate::cloning::pcr::outputs::Request,
) -> Result<AnalysisResult, Error> {
    {
        let origin = &origin;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::material::outputs_revision(origin, request)?
                != (*context.revision)
            {
                return Err(Error::Format("Material PCR output choices changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::pcr::material::prepare_outputs(&origin, &request)
            .map(AnalysisResult::MaterialPcrOutputsPrepared)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn install_material_derivation(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::material::derivation::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::revision(workspace, request)? != (*context.revision) {
                return Err(Error::Format(
                    "Material derivation workspace/source/recipe changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::material::derivation::install(&workspace, &request)
            .map(AnalysisResult::MaterialDerivationInstalled)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_digest(
    context: &RequestContext,
    source_role: String,
    source: Molecule,
    origin: crate::material::derivation::SourceOrigin,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    options: crate::cloning::restriction::inspection::InspectionOptions,
) -> Result<AnalysisResult, Error> {
    {
        let source_role = &source_role;
        let source = &source;
        let origin = &origin;
        let selection = &selection;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::digest::digest_inspection_revision(
                source_role,
                source,
                origin,
                selection,
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Material digest source/origin/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::derivation::digest::inspect_digest(
            &source_role,
            &source,
            &origin,
            &selection,
            &options,
        )
        .map(AnalysisResult::MaterialDigestInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_saved_run_review(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    run_id: String,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let run_id = &run_id;

        let molecule: Option<&Molecule> = {
            if crate::material::review::request_revision(workspace, run_id)? != (*context.revision)
            {
                return Err(Error::Format(
                    "Saved-run review workspace/run changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::material::review::prepare(&workspace, &run_id)
            .map(AnalysisResult::SavedRunReviewPrepared)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_mismatch_pcr(
    context: &RequestContext,
    source: Molecule,
    origin: crate::material::derivation::SourceOrigin,
    forward: crate::sequence::analysis::Primer,
    reverse: crate::sequence::analysis::Primer,
    options: crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let origin = &origin;
        let forward = &forward;
        let reverse = &reverse;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::pcr::mismatch_pcr_preview_revision(
                source, origin, forward, reverse, options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Material mismatch-PCR source/origin/primers/policy changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::derivation::pcr::preview_mismatch_pcr(
            &source, &origin, &forward, &reverse, &options,
        )
        .map(AnalysisResult::MaterialMismatchPcrPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_repair(
    context: &RequestContext,
    source: Molecule,
    origin: crate::material::derivation::SourceOrigin,
    method: crate::chemistry::EndRepairMethod,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let origin = &origin;
        let method = &method;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::end_repair::repair_preview_revision(
                source, origin, *method,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Material repair source/origin/policy changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::derivation::end_repair::preview_repair(&source, &origin, method)
            .map(AnalysisResult::MaterialRepairPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_overlap_extension(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    source: crate::material::lineage::EventBinding,
    origin: crate::material::derivation::SourceOrigin,
    preparation: crate::cloning::overlap::extension::pcr::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let source = &source;
        let origin = &origin;
        let preparation = &preparation;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::overlap_extension::overlap_extension_preview_revision(
                workspace,
                source,
                origin,
                preparation,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "OE material source/origin/preparation changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::derivation::overlap_extension::preview_overlap_extension(
            &workspace,
            &source,
            &origin,
            &preparation,
        )
        .map(AnalysisResult::MaterialOverlapExtensionPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_terminal_processing(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    source: crate::material::lineage::EventBinding,
    origin: crate::material::derivation::SourceOrigin,
    action: crate::material::derivation::TerminalAction,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let source = &source;
        let origin = &origin;
        let action = &action;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::terminal::terminal_processing_preview_revision(
                workspace, source, origin, action,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Terminal-processing material source/origin/action changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::derivation::terminal::preview_terminal_processing(
            &workspace, &source, &origin, &action,
        )
        .map(AnalysisResult::MaterialTerminalProcessingPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_annealed_oligos(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    top: crate::material::lineage::EventBinding,
    bottom: crate::material::lineage::EventBinding,
    options: crate::chemistry::AnnealOptions,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let top = &top;
        let bottom = &bottom;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::material::derivation::annealing::annealed_oligos_preview_revision(
                workspace, top, bottom, options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Annealed-material source events/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::material::derivation::annealing::preview_annealed_oligos(
            &workspace, &top, &bottom, &options,
        )
        .map(AnalysisResult::MaterialAnnealedOligosPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn materialize_material_digest_fragment(
    context: &RequestContext,
    source: Molecule,
    request: crate::material::derivation::digest_export::Request,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if request.revision(source)? != (*context.revision) {
                return Err(Error::Format(
                    "Material digest export context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = request
        .materialize(&source)
        .map(AnalysisResult::MaterialDigestProduct);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
