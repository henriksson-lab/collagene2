//! Independent assembly intent, preparation and candidate-selection handlers.
use super::context::RequestContext;
use crate::protocol::bindings::overlap_profile_revision;
use crate::protocol::schema::AnalysisResult;
use crate::protocol::validation::validate_verification_report;
use collagene_model::{Error, Molecule};

pub(super) fn check_overlap_profile(
    context: &RequestContext,
    molecules: Vec<Molecule>,
    reverse: Vec<bool>,
    overlaps: Vec<usize>,
    circular: bool,
    options: crate::cloning::overlap::profiles::ProfileOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecules = &molecules;
        let reverse = &reverse;
        let overlaps = &overlaps;
        let circular = &circular;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if overlap_profile_revision(molecules, reverse, overlaps, *circular, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "overlap profile request source/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::overlap::profiles::check_oriented_fragments(
            &molecules, &reverse, &overlaps, circular, &options,
        )
        .map(AnalysisResult::OverlapProfile);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_overlap_recipe(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    recipe: crate::cloning::assembly::recipe::Recipe,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let recipe = &recipe;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::recipe::preview_revision(inputs, target, recipe)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "overlap preparation source/target/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::OverlapRecipePreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_ligation_recipe(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    recipe: crate::cloning::ligation::recipe::Recipe,
    policy: crate::cloning::assembly::recipe::CleavagePolicy,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let recipe = &recipe;
        let policy = &policy;

        let molecule: Option<&Molecule> = {
            if crate::cloning::ligation::recipe::preview_revision(inputs, target, recipe, *policy)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "ligation preparation source/target/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::ligation::recipe::preview(&inputs, &target, &recipe, policy)
            .map(AnalysisResult::LigationRecipePreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_assembly_junctions(
    context: &RequestContext,
    plan: crate::planning::plan::DesignPlan,
    target: Molecule,
    options: crate::cloning::assembly::junctions::Options,
) -> Result<AnalysisResult, Error> {
    {
        let plan = &plan;
        let target = &target;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::junctions::inspection_revision(plan, target, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "assembly junction plan/target/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::junctions::inspect_plan(&plan, &target, &options)
            .map(AnalysisResult::AssemblyJunctionInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_overlap_primers(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    request: crate::cloning::overlap::primer_design::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::overlap::primer_design::design_revision(inputs, target, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "overlap primer sources/target/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::overlap::primer_design::design(&inputs, &target, &request)
            .map(AnalysisResult::OverlapPrimersDesigned);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn build_assembly_intent(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::assembly::intent::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::intent::intent_revision(
                inputs,
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "assembly intent source/region/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::intent::build(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::AssemblyIntent);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_golden_gate_primers(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    request: crate::cloning::golden_gate::primers::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::primers::preview_revision(inputs, target, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Golden Gate primer sources/target/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::golden_gate::primers::design(&inputs, &target, &request)
            .map(AnalysisResult::GoldenGatePrimersDesigned);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn build_mutation_intent(
    context: &RequestContext,
    molecule: Molecule,
    request: crate::cloning::pcr::mutation_intent::Request,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::mutation_intent::intent_revision(molecule, request)?
                != (*context.revision)
            {
                return Err(Error::Format("Mutation intent source/edit changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::pcr::mutation_intent::build(&molecule, &request)
            .map(AnalysisResult::MutationIntent)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn build_edited_assembly_intent(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    base: crate::cloning::assembly::intent::Request,
    edits: Vec<crate::cloning::assembly::junction_edit::Edit>,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let base = &base;
        let edits = &edits;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::junction_edit::intent_revision(
                inputs,
                base,
                edits,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Edited assembly source/intent changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::junction_edit::build_intent(
            &inputs,
            &base,
            &edits,
            workspace.as_ref(),
        )
        .map(AnalysisResult::EditedAssemblyIntent);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_edited_assembly_primers(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::assembly::junction_edit::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::junction_edit::design_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Edited assembly source/primers changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::assembly::junction_edit::design(&inputs, &request)
            .map(AnalysisResult::EditedAssemblyPrimers)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_selective_assembly_primers(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::assembly::selective_edit::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::selective_edit::design_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Selective assembly source/primers changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::assembly::selective_edit::design(&inputs, &request)
            .map(AnalysisResult::SelectiveAssemblyPrimers)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_golden_gate_junction_primers(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::golden_gate::junction_edit::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::junction_edit::design_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Golden Gate junction source/primers changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::junction_edit::design(&inputs, &request)
            .map(AnalysisResult::GoldenGateJunctionPrimers)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_prepared_golden_gate_pool(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    request: crate::cloning::golden_gate::prepared::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::prepared::inspect_revision(inputs, target, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Prepared Golden Gate pool context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::golden_gate::prepared::inspect(&inputs, &target, &request)
            .map(AnalysisResult::PreparedGoldenGatePool);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn select_prepared_golden_gate_candidate(
    context: &RequestContext,
    report: Box<crate::cloning::golden_gate::prepared::Report>,
    selection: crate::cloning::golden_gate::prepared::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let report = &report;
        let selection = &selection;
        validate_verification_report(&report.preview.verification)?;
        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::prepared::selection_revision(report, selection)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Prepared Golden Gate selection context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::prepared::select(&report, &selection)
            .map(AnalysisResult::PreparedGoldenGateCandidate)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_golden_gate_reaction(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    request: crate::cloning::golden_gate::reaction::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::reaction::inspect_revision(inputs, target, request)?
                != (*context.revision)
            {
                return Err(Error::Format("Golden Gate reaction context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::golden_gate::reaction::inspect(&inputs, &target, &request)
            .map(AnalysisResult::GoldenGateReaction);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn select_golden_gate_reaction(
    context: &RequestContext,
    report: Box<crate::cloning::golden_gate::reaction::Report>,
    selection: crate::cloning::golden_gate::prepared::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let report = &report;
        let selection = &selection;
        validate_verification_report(&report.verification)?;
        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::reaction::selection_revision(report, selection)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Golden Gate reaction selection changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::reaction::select(&report, &selection)
            .map(AnalysisResult::GoldenGateReactionCandidate)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn search_golden_gate_fusions(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::golden_gate::fusion_design::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::fusion_design::inspect_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Golden Gate fusion search context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::fusion_design::inspect(&inputs, &request)
            .map(AnalysisResult::GoldenGateFusions)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn materialize_golden_gate_fusions(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::golden_gate::fusion_design::Request,
    selection: crate::cloning::golden_gate::fusion_design::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::fusion_design::materialize_revision(
                inputs, request, selection,
            )? != (*context.revision)
            {
                return Err(Error::Format("Golden Gate fusion selection changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::golden_gate::fusion_design::materialize(&inputs, &request, &selection)
            .map(AnalysisResult::GoldenGateFusionMaterialized);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_domestication(
    context: &RequestContext,
    source: Molecule,
    request: crate::cloning::golden_gate::domestication::Request,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::domestication::inspect_revision(source, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Domestication inspection context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::domestication::inspect(&source, &request)
            .map(AnalysisResult::DomesticationInspected)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn propose_domestication(
    context: &RequestContext,
    source: Molecule,
    request: crate::cloning::golden_gate::domestication::Request,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::domestication::propose_revision(source, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Domestication proposal context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::domestication::propose(&source, &request)
            .map(AnalysisResult::DomesticationProposed)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn check_domestication(
    context: &RequestContext,
    source: Molecule,
    request: crate::cloning::golden_gate::domestication::Request,
    selection: crate::cloning::golden_gate::domestication::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let request = &request;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::domestication::check_revision(
                source, request, selection,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Domestication substitution context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::golden_gate::domestication::check(
        &source, &request, &selection,
    )
    .map(|checked| AnalysisResult::DomesticationChecked {
        report: checked.report().clone(),
        target: checked.target().clone(),
    });
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_domestication(
    context: &RequestContext,
    input: crate::planning::plan::PlanInput,
    request: crate::cloning::golden_gate::domestication::preparation::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let input = &input;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::domestication::preparation::revision(
                input,
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Domestication preparation context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::golden_gate::domestication::preparation::materialize(
            &input,
            &request,
            workspace.as_ref(),
        )
        .map(AnalysisResult::DomesticationPrepared);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn rank_golden_gate_fusions(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::golden_gate::fusion_rank::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::golden_gate::fusion_rank::inspect_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Fusion ranking source/profile context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::golden_gate::fusion_rank::inspect(&inputs, &request)
            .map(AnalysisResult::RankedGoldenGateFusions)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_tailing_recipe(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    recipe: crate::cloning::ligation::tailing_recipe::Recipe,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let recipe = &recipe;

        let molecule: Option<&Molecule> = {
            if crate::cloning::ligation::tailing_recipe::preview_revision(inputs, target, recipe)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "tailing preparation source/target/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::ligation::tailing_recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::TailingRecipePreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn build_substituted_assembly_intent(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::assembly::substitution_intent::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::substitution_intent::intent_revision(
                inputs,
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Substituted assembly original source/declaration changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::substitution_intent::build(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::SubstitutedAssemblyIntent);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn build_contribution_edit_intent(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::assembly::edit_intent::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::edit_intent::intent_revision(
                inputs,
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Contribution edit original source/declaration changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::edit_intent::build(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::ContributionEditIntent);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_edited_overlap_primers(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    request: crate::cloning::assembly::edited_primers::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::edited_primers::design_revision(
                inputs,
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Edited overlap primer original sources, intent, origins or policy changed"
                        .into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::assembly::edited_primers::design(&inputs, &request, workspace.as_ref())
            .map(AnalysisResult::EditedOverlapPrimersDesigned)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_homology_recipe(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    recipe: crate::cloning::assembly::recipe::Recipe,
    policy: crate::cloning::overlap::assembly::Options,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let recipe = &recipe;
        let policy = &policy;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::recipe::homology_preview_revision(
                inputs, target, recipe, policy,
            )? != (*context.revision)
            {
                return Err(Error::Format("Homology recipe source, target, preparation or explicit physical policy changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::assembly::recipe::preview_homology(&inputs, &target, &recipe, &policy)
            .map(AnalysisResult::OverlapRecipePreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_annealing_recipe(
    context: &RequestContext,
    inputs: Vec<crate::planning::plan::PlanInput>,
    target: Molecule,
    recipe: crate::cloning::assembly::recipe::Recipe,
    policy: crate::cloning::overlap::annealing::Options,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let target = &target;
        let recipe = &recipe;
        let policy = &policy;

        let molecule: Option<&Molecule> = {
            if crate::cloning::assembly::recipe::annealing_preview_revision(
                inputs, target, recipe, policy,
            )? != (*context.revision)
            {
                return Err(Error::Format("3-prime resection/annealing recipe source, target, preparation or explicit physical policy changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::assembly::recipe::preview_annealing(&inputs, &target, &recipe, &policy)
            .map(AnalysisResult::OverlapRecipePreview)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
