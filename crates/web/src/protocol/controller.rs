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

/// Retained run summaries and their current workspace presentation.
#[derive(Clone, Copy)]
pub(crate) struct RunHistoryState {
    pub(crate) recorded_runs: RwSignal<Vec<collagene_core::material::lineage::RunSummary>>,
    pub(crate) run_page: RwSignal<usize>,
    pub(crate) pending_run_event:
        RwSignal<Option<std::sync::Arc<crate::history::run_event_model::Snapshot>>>,
    pub(crate) hidden_products: RwSignal<std::collections::BTreeSet<String>>,
    pub(crate) show_prepared_history: RwSignal<bool>,
}
impl RunHistoryState {
    pub(crate) fn new() -> Self {
        Self {
            recorded_runs: RwSignal::new(
                Vec::<collagene_core::material::lineage::RunSummary>::new(),
            ),
            run_page: RwSignal::new(0usize),
            pending_run_event: RwSignal::new(
                None::<std::sync::Arc<crate::history::run_event_model::Snapshot>>,
            ),
            hidden_products: RwSignal::new(std::collections::BTreeSet::<String>::new()),
            show_prepared_history: RwSignal::new(false),
        }
    }
    pub(crate) fn restore(
        &self,
        runs: Vec<collagene_core::material::lineage::RunSummary>,
        hidden_products: std::collections::BTreeSet<String>,
    ) {
        self.recorded_runs.set(runs);
        self.hidden_products.set(hidden_products);
        self.run_page.set(0);
        self.pending_run_event.set(None);
        self.show_prepared_history.set(false);
    }
}

/// Candidate declarations and transient source-bound verification evidence.
#[derive(Clone, Copy)]
pub(crate) struct PlanState {
    pub(crate) candidate: RwSignal<Option<collagene_core::documents::transfer::Candidate>>,
    pub(crate) stored_candidates: RwSignal<Vec<collagene_core::documents::transfer::Candidate>>,
    pub(crate) candidate_events: RwSignal<
        Option<(
            String,
            u64,
            collagene_core::material::lineage::EventBinding,
            Vec<collagene_core::material::lineage::InputEventBinding>,
        )>,
    >,
    pub(crate) plan_report:
        RwSignal<Option<collagene_core::verification::evidence::VerificationReport>>,
    pub(crate) last_plan: RwSignal<
        Option<(
            collagene_core::planning::plan::DesignPlan,
            collagene_model::Molecule,
        )>,
    >,
    pub(crate) last_plan_target:
        RwSignal<Option<(String, u64, collagene_core::material::lineage::EventBinding)>>,
    pub(crate) plan_generation: RwSignal<u32>,
}
impl PlanState {
    pub(crate) fn new() -> Self {
        Self {
            candidate: RwSignal::new(None::<transfer::Candidate>),
            stored_candidates: RwSignal::new(Vec::<transfer::Candidate>::new()),
            candidate_events: RwSignal::new(
                None::<(
                    String,
                    u64,
                    collagene_core::material::lineage::EventBinding,
                    Vec<collagene_core::material::lineage::InputEventBinding>,
                )>,
            ),
            plan_report: RwSignal::new(
                None::<collagene_core::verification::evidence::VerificationReport>,
            ),
            last_plan: RwSignal::new(None::<(DesignPlan, Molecule)>),
            last_plan_target: RwSignal::new(
                None::<(String, u64, collagene_core::material::lineage::EventBinding)>,
            ),
            plan_generation: RwSignal::new(0u32),
        }
    }
    pub(crate) fn invalidate_evidence(&self) {
        self.plan_generation
            .update(|generation| *generation = generation.wrapping_add(1));
        self.plan_report.set(None);
        self.last_plan.set(None);
        self.last_plan_target.set(None);
    }
    pub(crate) fn restore(
        &self,
        candidate: Option<collagene_core::documents::transfer::Candidate>,
        candidates: Vec<collagene_core::documents::transfer::Candidate>,
    ) {
        self.invalidate_evidence();
        self.stored_candidates.set(candidates);
        self.candidate_events.set(None);
        self.candidate.set(candidate);
    }
}
