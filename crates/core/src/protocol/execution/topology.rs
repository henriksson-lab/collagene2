//! Physical TOPO graphs and bounded reference-recombination handlers.
use super::context::RequestContext;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Error, Molecule};

pub(super) fn inspect_recombination(
    context: &RequestContext,
    inputs: Vec<crate::cloning::recombination::reference::Input>,
    request: crate::cloning::recombination::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::recombination::batch_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "recombination source/profile/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::recombination::simulate(&inputs, &request)
            .map(AnalysisResult::RecombinationInspection)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn check_recombination_target(
    context: &RequestContext,
    inputs: Vec<crate::cloning::recombination::reference::Input>,
    request: crate::cloning::recombination::Request,
    target: Molecule,
    policy: crate::cloning::recombination::reference::TargetPolicy,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let target = &target;
        let policy = &policy;

        let molecule: Option<&Molecule> = {
            if crate::cloning::recombination::check_revision(inputs, request, target, policy)?
                != (*context.revision)
            {
                return Err(Error::Format("recombination target binding changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::recombination::check_target(&inputs, &request, &target, &policy)
            .map(AnalysisResult::RecombinationTargetCheck);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn select_recombination_component(
    context: &RequestContext,
    inputs: Vec<crate::cloning::recombination::reference::Input>,
    request: crate::cloning::recombination::Request,
    selection: crate::cloning::recombination::Selection,
    target: Molecule,
    policy: crate::cloning::recombination::reference::TargetPolicy,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let selection = &selection;
        let target = &target;
        let policy = &policy;

        let molecule: Option<&Molecule> = {
            if crate::cloning::recombination::checked_product_revision(
                inputs, request, selection, target, policy,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "recombination selected component binding changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::recombination::checked_product(
        &inputs, &request, &selection, &target, &policy,
    )
    .map(AnalysisResult::RecombinationProduct);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn simulate_topo(
    context: &RequestContext,
    inputs: Vec<crate::cloning::topo::Input>,
    request: crate::cloning::topo::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::worker::simulation_revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format("TOPO simulation context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { crate::cloning::topo::simulate(&inputs, &request).map(AnalysisResult::TopoSimulation) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_topo_orientations(
    context: &RequestContext,
    inputs: Vec<crate::cloning::topo::Input>,
    request: crate::cloning::topo::Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let vector = &vector;
        let insert = &insert;
        let directional = &directional;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::worker::orientations_revision(
                inputs,
                request,
                *vector,
                *insert,
                *directional,
            )? != (*context.revision)
            {
                return Err(Error::Format("TOPO orientation context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::topo::inspect_orientations(&inputs, &request, vector, insert, directional)
            .map(AnalysisResult::TopoOrientations);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn project_topo_product(
    context: &RequestContext,
    inputs: Vec<crate::cloning::topo::Input>,
    request: crate::cloning::topo::Request,
    selection: crate::cloning::topo::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::worker::projection_revision(inputs, request, selection)?
                != (*context.revision)
            {
                return Err(Error::Format("TOPO projection context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::topo::product(&inputs, &request, &selection)
            .map(AnalysisResult::TopoProjection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn check_topo_target(
    context: &RequestContext,
    inputs: Vec<crate::cloning::topo::Input>,
    request: crate::cloning::topo::Request,
    selection: crate::cloning::topo::Selection,
    target: Molecule,
    policy: crate::cloning::topo::AcceptancePolicy,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let selection = &selection;
        let target = &target;
        let policy = &policy;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::worker::target_revision(
                inputs, request, selection, target, policy,
            )? != (*context.revision)
            {
                return Err(Error::Format("TOPO target/policy context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::topo::check_target(&inputs, &request, &selection, &target, &policy)
            .map(AnalysisResult::TopoTarget);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn process_topo_actions(
    context: &RequestContext,
    inputs: Vec<crate::cloning::topo::Input>,
    request: crate::cloning::topo::actions::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::actions::revision(inputs, request)? != (*context.revision) {
                return Err(Error::Format("TOPO actions product context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::topo::actions::inspect(&inputs, &request)
            .map(AnalysisResult::TopoActionsProduct)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn simulate_multi_site_recombination(
    context: &RequestContext,
    inputs: Vec<crate::cloning::recombination::reference::Input>,
    request: crate::cloning::recombination::multisite::Request,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::recombination::multisite::revision(inputs, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "multisite sources/profile/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::recombination::multisite::simulate(&inputs, &request)
            .map(AnalysisResult::MultiSiteRecombination)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn project_multi_site_recombination(
    context: &RequestContext,
    inputs: Vec<crate::cloning::recombination::reference::Input>,
    request: crate::cloning::recombination::multisite::Request,
    selection: crate::cloning::recombination::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::cloning::recombination::multisite::product_revision(
                inputs, request, selection,
            )? != (*context.revision)
            {
                return Err(Error::Format("multisite selected component changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::recombination::multisite::product(&inputs, &request, &selection)
            .map(AnalysisResult::MultiSiteRecombinationProduct);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn check_multi_site_recombination(
    context: &RequestContext,
    inputs: Vec<crate::cloning::recombination::reference::Input>,
    request: crate::cloning::recombination::multisite::Request,
    selection: crate::cloning::recombination::Selection,
    target: Molecule,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let request = &request;
        let selection = &selection;
        let target = &target;

        let molecule: Option<&Molecule> = {
            if crate::cloning::recombination::multisite::check_revision(
                inputs, request, selection, target,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "multisite independent target context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::recombination::multisite::check_target(
            &inputs, &request, &selection, &target,
        )
        .map(AnalysisResult::MultiSiteRecombinationChecked);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_topo_run(
    context: &RequestContext,
    inputs: Vec<crate::cloning::topo::Input>,
    template: Option<crate::cloning::topo::run_preparation::BoundTemplate>,
    task: crate::cloning::topo::run_preparation::Task,
) -> Result<AnalysisResult, Error> {
    {
        let inputs = &inputs;
        let template = &template;
        let task = &task;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::run_preparation::revision(inputs, template.as_ref(), task)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Saved TOPO source, origin or inspection context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::topo::run_preparation::inspect(&inputs, template.as_ref(), &task)
            .map(AnalysisResult::TopoRunInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn list_topo_graphs(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::graph_history::list_revision(workspace)? != (*context.revision)
            {
                return Err(Error::Format("TOPO graph listing workspace changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::topo::graph_history::list(&workspace).map(AnalysisResult::TopoGraphsListed)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn capture_topo_graph(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::cloning::topo::graph_history::CaptureRequest,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::graph_history::revision(workspace, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "TOPO graph capture source/request changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::topo::graph_history::prepare_capture(&workspace, &request)
            .map(AnalysisResult::TopoGraphCaptured)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_topo_graph(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    selection: crate::cloning::topo::graph_history::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::cloning::topo::graph_history::inspect_revision(workspace, selection)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "TOPO graph selection/workspace changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::topo::graph_history::inspect(&workspace, &selection)
            .map(AnalysisResult::TopoGraphInspected);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
