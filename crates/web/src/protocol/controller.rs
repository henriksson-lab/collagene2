use crate::platform::worker;
use crate::workspace::session::current_event_binding;
use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::protocol;
use collagene_model::Molecule;
use leptos::prelude::*;
use wasm_bindgen_futures::spawn_local;

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
    pub(crate) fn run_plan(
        self,
        state: PlanState,
    ) -> impl Fn(
        DesignPlan,
        Molecule,
        Option<Vec<collagene_core::material::lineage::InputEventBinding>>,
        Option<EventBinding>,
    ) + Copy
    + 'static {
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        let PlanState {
            plan_report,
            last_plan,
            last_plan_target,
            last_plan_events,
            plan_generation,
            ..
        } = state;
        move |plan: DesignPlan,
              target: Molecule,
              exact_sources: Option<Vec<collagene_core::material::lineage::InputEventBinding>>,
              exact_target: Option<collagene_core::material::lineage::EventBinding>| {
            let expected_epoch = workspace_epoch.get_untracked();
            let expected_target = exact_target.or_else(|| {
                documents.with_untracked(|docs| {
                    docs.get(active.get_untracked()).map(current_event_binding)
                })
            });
            let expected_sources = documents.with_untracked(|docs| {
                docs.iter()
                    .map(|d| {
                        (
                            d.id.clone(),
                            d.head_event().to_string(),
                            d.head_revision().to_string(),
                        )
                    })
                    .collect::<Vec<_>>()
            });
            state.invalidate_evidence();
            let revision = target.revision();
            let generation = plan_generation.get_untracked();
            let request = protocol::request::Request {
                id: 0,
                revision: revision.clone(),
                job: protocol::schema::AnalysisJob::Verify {
                    target: target.clone(),
                    plan: plan.clone(),
                },
            };
            message.set("Reconstructing the assembly in background…".into());
            spawn_local(async move {
                let result = worker::analyze(request).await;
                if plan_generation.try_get_untracked() != Some(generation)
                    || workspace_epoch.get_untracked() != expected_epoch
                    || !expected_target.as_ref().is_some_and(|binding| {
                        documents.with_untracked(|docs| {
                            crate::protocol::plan_editing_sources::resolve_target(
                                docs, binding, &revision,
                            )
                            .is_ok()
                        })
                    })
                    || documents.with_untracked(|docs| {
                        docs.iter()
                            .map(|d| {
                                (
                                    d.id.clone(),
                                    d.head_event().to_string(),
                                    d.head_revision().to_string(),
                                )
                            })
                            .collect::<Vec<_>>()
                    }) != expected_sources
                {
                    return;
                }
                match result {
                    Ok(protocol::response::Response {
                        result: protocol::schema::AnalysisResult::Verification(report),
                        ..
                    }) => {
                        plan_report.set(Some(report));
                        last_plan_events.set(exact_sources.map(|bindings| {
                            (
                                collagene_core::planning::plan::DesignPlan::revision(&plan),
                                bindings,
                            )
                        }));
                        last_plan_target.set(expected_target.map(|event| {
                            (
                                collagene_core::planning::plan::DesignPlan::revision(&plan),
                                expected_epoch,
                                event,
                            )
                        }));
                        last_plan.set(Some((plan, target)));
                        message.set("Assembly reconstruction finished.".into());
                    }
                    Ok(protocol::response::Response {
                        result: protocol::schema::AnalysisResult::Error(error),
                        ..
                    })
                    | Err(error) => message.set(error),
                    _ => message.set("Unexpected worker response".into()),
                }
            });
        }
    }
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
    pub(crate) last_plan_events: RwSignal<
        Option<(
            String,
            Vec<collagene_core::material::lineage::InputEventBinding>,
        )>,
    >,
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
            last_plan_events: RwSignal::new(
                None::<(
                    String,
                    Vec<collagene_core::material::lineage::InputEventBinding>,
                )>,
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
        self.last_plan_events.set(None);
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
