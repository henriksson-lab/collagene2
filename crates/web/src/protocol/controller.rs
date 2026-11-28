use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::protocol;
use collagene_core::verification;
use collagene_core::{
    documents::workspace, material::lineage as run_lineage,
    output::quantities as protocol_quantities,
};
use collagene_model::Molecule;
use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct ProtocolController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl ProtocolController {
    pub(crate) fn quantity_binding(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_events: RwSignal<
            Option<(
                String,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<(), Result<crate::protocol::protocol_quantities_model::Binding, String>> {
        let session = self.session;
        let WorkspaceSession {
            workspace_epoch, ..
        } = self.session;
        Callback::new(
            move |()| -> Result<crate::protocol::protocol_quantities_model::Binding, String> {
                last_plan.with_untracked(|saved|->Result<_,String>{
            let (plan,target)=saved.as_ref().ok_or("Reconstruct the exact plan before declaring quantities")?;
            let (revision, epoch, event) = last_plan_target
                .get_untracked()
                .ok_or("Reconstruct against the exact intended target first")?;
            if epoch != workspace_epoch.get_untracked() || target.revision() != event.revision {
                return Err("The quantity target context changed.".into());
            }
            let inputs=last_plan_events.get_untracked().filter(|(id,_)|id==&revision).map(|(_,pins)|pins).ok_or("Choose exact input events in the plan editor, then verify before declaring quantities.")?;
            let retained = candidate
                .with_untracked(|c| c.as_ref().is_some_and(|c| c.plan_revision == revision))
                || stored_candidates
                    .with_untracked(|all| all.iter().any(|c| c.plan_revision == revision));
            if !retained {
                return Err(
                    "Open or keep this exact candidate before declaring source-bound quantities."
                        .into(),
                );
            }
            Ok(crate::protocol::protocol_quantities_model::Binding {
                plan_source: collagene_core::output::quantities::PlanSource::Candidate,
                plan_revision: revision,
                target: event,
                inputs,
                step_count: plan.steps.len(),
                operation_kinds: plan.steps.iter().map(std::mem::discriminant).collect(),
                component_sources: plan
                    .steps
                    .iter()
                    .map(|op| {
                        collagene_core::planning::plan::Operation::sources(op)
                            .into_iter()
                            .cloned()
                            .collect()
                    })
                    .collect(),
            })
            })
            },
        )
    }
}
