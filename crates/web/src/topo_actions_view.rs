//! Explicit ordered physical processing after a freshly inspected TOPO path.
//! No target-driven deletion, implicit repair, or advisory-plan handoff.
use crate::worker::CachedWorker;
use collagene_core::{
    protocol::AnalysisResult,
    topo::{
        self,
        physical::{self, processing},
    },
    topo_actions,
};
use leptos::prelude::*;
use std::{
    io::{self, Write},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, PartialEq, Eq)]
pub struct Frame {
    pub epoch: u64,
    pub sources: Vec<(String, collagene_core::run_lineage::EventBinding)>,
    pub target: collagene_core::run_lineage::EventBinding,
    pub inspection_revision: String,
    pub inspection_generation: u32,
    pub orientation: usize,
    pub graph: Option<collagene_core::topo_graph_history::Selection>,
}
/// An advisory declaration, never a positive material receipt.
pub struct OpenActions {
    pub frame: Frame,
    pub request: topo_actions::Request,
    pub expected_revision: String,
}
/// Exact reviewed request plus a cheap source/form guard. No imported report
/// or selected contour becomes graph authority; capture reruns the native recipe.
#[derive(Clone)]
pub struct ReviewedReaction {
    pub frame: Frame,
    pub request: Arc<topo::Request>,
    pub guard: Callback<(), Result<(), String>>,
}
pub struct Capture {
    pub inputs: Vec<topo::Input>,
    pub template: Option<collagene_core::topo_run_preparation::BoundTemplate>,
    pub request: topo::Request,
    pub physical: physical::Report,
}
#[derive(Clone, PartialEq, Eq, Default)]
struct Form {
    actions: Vec<processing::Action>,
    hydrolysis: String,
    ligase: String,
    upstream_strand: String,
    upstream_call: String,
    downstream_strand: String,
    downstream_call: String,
    displacement: String,
}
#[derive(Clone, PartialEq, Eq)]
struct ContourForm {
    component: String,
    strand: String,
    call: String,
    require_all: bool,
    require_known: bool,
}
impl Default for ContourForm {
    fn default() -> Self {
        Self {
            component: String::new(),
            strand: String::new(),
            call: String::new(),
            require_all: true,
            require_known: true,
        }
    }
}
struct Inspected {
    frame: Frame,
    form: Form,
    request: Arc<topo::Request>,
    report: Arc<topo::TopoReport>,
}
struct Projected {
    inspection: Arc<Inspected>,
    choice: ContourForm,
    report: Arc<topo_actions::Report>,
}
struct Proposal {
    frame: Frame,
    form: Form,
    action: processing::Action,
}
fn call_text(call: physical::Call) -> String {
    format!("source strand {}, call {}", call.strand, call.index + 1)
}
fn action_text(action: &processing::Action) -> String {
    match action {
        processing::Action::CleaveThreePrimeFlap {
            bond,
            displacement_step,
            assumption,
            ..
        } => format!(
            "Hydrolytic 3′ flap cleavage: bond {} → {}; recorded displacement transfer {}. Assumption: {}",
            call_text(bond.upstream),
            call_text(bond.downstream),
            displacement_step + 1,
            assumption
        ),
        processing::Action::SealNick {
            selection,
            ligase_assumption,
        } => format!(
            "Seal nick {} (3′) → {} (5′). Ligase assumption: {}",
            call_text(selection.three),
            call_text(selection.five),
            ligase_assumption
        ),
    }
}
fn outcome_text(step: &processing::ActionEvidence) -> String {
    match step {
        processing::ActionEvidence::FlapCleavage(e) => format!(
            "Action {} · {:?}: {}. {} detached source calls; created 3′ {:?} and 5′ {:?}; anchor template {}.",
            e.step + 1,
            e.outcome,
            e.reason,
            e.flap_calls.len(),
            e.created_three_prime,
            e.created_five_prime,
            e.anchor_template
                .map(call_text)
                .unwrap_or_else(|| "not established".into())
        ),
        processing::ActionEvidence::NickSealing(e) => format!(
            "Action {} · {:?}: {}. Actual 3′ {:?}, 5′ {:?}; paired templates {} / {}.",
            e.step + 1,
            e.outcome,
            e.reason,
            e.three_prime,
            e.five_prime,
            e.template_three
                .map(call_text)
                .unwrap_or_else(|| "not established".into()),
            e.template_five
                .map(call_text)
                .unwrap_or_else(|| "not established".into())
        ),
    }
}
pub(crate) fn step_text(step: &physical::StepEvidence) -> String {
    match step {
        physical::StepEvidence::Transfer(e) => format!(
            "Step {} transfer {} → {} · {:?} · {} · displacement {:?}",
            e.step + 1,
            call_text(e.donor_call),
            call_text(e.acceptor_call),
            e.outcome,
            e.reason,
            e.displacement.as_ref().map(|d| d.kind)
        ),
        physical::StepEvidence::Process(e) => outcome_text(e),
    }
}
struct Json(Vec<u8>);
impl Write for Json {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > 32 * 1024 * 1024)
        {
            return Err(io::Error::other("Actions evidence exceeds 32 MiB"));
        }
        self.0
            .try_reserve(b.len())
            .map_err(|_| io::Error::other("Evidence allocation failed"))?;
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
#[component]
fn Rows(label: &'static str, total: usize, row: Callback<usize, String>) -> impl IntoView {
    let page = RwSignal::new(0usize);
    view! {<details><summary>{format!("{label} ({total})")}</summary>
        <p>{move||format!("Page {} of {}",page.get()+1,total.div_ceil(40).max(1))}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous"</button>
        <button class="text-tool" disabled={move||(page.get()+1)*40>=total} on:click=move|_|page.update(|p|*p+=1)>"Next"</button>
        {move||(page.get()*40..((page.get()+1)*40).min(total)).map(|i|view!{<p>{row.run(i)}</p>}).collect_view()}
    </details>}
}
fn graph_view(graph: Arc<physical::ProductGraph>) -> impl IntoView {
    let components = graph.clone();
    let chains = graph.clone();
    let nicks = graph.clone();
    let runs = graph.clone();
    view! {<section>
        <p>{format!("{} components · {} remaining nicks · {} unpaired runs. Every free chain remains in the graph.",graph.components.len(),graph.nicks.len(),graph.unpaired.len())}</p>
        <Rows label="All retained components" total=graph.components.len() row=Callback::new(move|i:usize|{let c=&components.components[i];format!("Component {i}: {} paired calls / {} unpaired calls; chains {:?}; covalently closed {}",c.paired_nucleotides,c.unpaired_nucleotides,c.chains,c.covalently_closed)})/>
        <Rows label="Actual retained chains" total=graph.chains.len() row=Callback::new(move|i:usize|{let c=&chains.chains[i];let end=c.sequence.len().min(200);format!("Chain {i}: {} calls, {:?}. Stored 5′→3′ sequence{}: {}",c.sequence.len(),c.ends,if end<c.sequence.len(){" (first 200 calls; full sequence in evidence JSON)"}else{""},&c.sequence.as_str()[..end])})/>
        <Rows label="Remaining nicks and chemistry" total=graph.nicks.len() row=Callback::new(move|i:usize|{let n=&nicks.nicks[i];format!("{} ({:?}, 3′) → {} ({:?}, 5′)",call_text(n.three),n.three_prime,call_text(n.five),n.five_prime)})/>
        <Rows label="Unpaired runs and displaced material" total=graph.unpaired.len() row=Callback::new(move|i:usize|{let r=&runs.unpaired[i];format!("Chain {} calls {}–{}; 5′ terminal {}; 3′ terminal {}; displaced by transfers {:?}",r.chain,r.start,r.end,r.terminal_five,r.terminal_three,r.displaced_by)})/>
    </section>}
}
#[component]
pub fn TopoActionsView(
    frame: Memo<Option<Frame>>,
    capture: Callback<Frame, Result<Capture, String>>,
    on_open: Callback<OpenActions, Result<(), String>>,
    on_reviewed: Callback<Option<Arc<ReviewedReaction>>>,
    locked: Memo<bool>,
) -> impl IntoView {
    let captured = RwSignal::new(None::<(Frame, Arc<Capture>)>);
    let form = RwSignal::new(Form::default());
    let proposal = RwSignal::new(None::<Arc<Proposal>>);
    let inspected = RwSignal::new(None::<Arc<Inspected>>);
    let choice = RwSignal::new(ContourForm::default());
    let projected = RwSignal::new(None::<Arc<Projected>>);
    let origin_review = RwSignal::new(None::<crate::topo_preparation::OriginReview>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cancel = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = frame.get();
        origin_review.set(None);
        on_reviewed.run(None);
        cancel.run(());
        captured.set(None);
        form.set(Form::default());
        proposal.set(None);
        inspected.set(None);
        projected.set(None);
        choice.set(ContourForm::default());
        status.set(String::new());
    });
    Effect::new(move |_| {
        if locked.get() {
            cancel.run(());
            proposal.set(None);
            projected.set(None);
        }
    });
    let action_key = Memo::new(move |_| form.with(|f| f.actions.clone()));
    Effect::new(move |_| {
        let _ = action_key.get();
        on_reviewed.run(None);
        cancel.run(());
        proposal.set(None);
        inspected.set(None);
        projected.set(None);
        choice.set(ContourForm::default());
    });
    Effect::new(move |_| {
        let _ = choice.get();
        cancel.run(());
        projected.set(None);
    });
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
    });
    let current_capture = move || {
        captured
            .get()
            .filter(|(key, _)| !locked.get() && frame.get().as_ref() == Some(key))
    };
    let current_inspection = move || {
        inspected.get().filter(|value| {
            !locked.get()
                && frame.get().as_ref() == Some(&value.frame)
                && form.with(|f| f.actions == value.form.actions)
        })
    };
    let current_projection = move || {
        projected.get().filter(|value| {
            choice.get() == value.choice
                && current_inspection().is_some_and(|i| Arc::ptr_eq(&i, &value.inspection))
        })
    };
    let load = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some(key) = frame.get_untracked() else {
            return;
        };
        match capture.run(key.clone()) {
            Ok(value) => {
                cancel.run(());
                form.set(Form::default());
                on_reviewed.run(None);
                inspected.set(None);
                projected.set(None);
                captured.set(Some((key, Arc::new(value))));
                status.set(
                    "Original path loaded. No processing actions have been applied or selected."
                        .into(),
                );
            }
            Err(e) => status.set(e),
        }
    });
    let add = Callback::new(move |action: processing::Action| {
        if locked.get_untracked() || busy.get_untracked() || current_capture().is_none() {
            return;
        }
        let Some((_, base)) = current_capture() else {
            return;
        };
        if form.with_untracked(|f| {
            f.actions.len().saturating_add(base.physical.steps.len())
                >= base.physical.limits.max_steps
        }) {
            status.set("The complete reaction has reached its explicit step limit.".into());
            return;
        }
        form.update(|f| f.actions.push(action));
        status.set(
            "Action added to the ordered declaration. Inspect to reconstruct actual chemistry."
                .into(),
        );
    });
    let propose = Callback::new(move |step: usize| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some((key, value)) = current_capture() else {
            return;
        };
        let raw = form.get_untracked();
        if raw.hydrolysis.trim().is_empty() || raw.hydrolysis.len() > 4096 {
            status.set("Declare the ideal hydrolytic cleavage assumption (1–4096 bytes).".into());
            return;
        }
        match processing::flap_proposals(&value.physical){
            Ok(values)=>match values.into_iter().find(|p|p.displacement_step==step){
                Some(p)=>proposal.set(Some(Arc::new(Proposal{frame:key,form:raw.clone(),action:processing::Action::CleaveThreePrimeFlap{bond:p.bond,displacement_step:p.displacement_step,profile:processing::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,assumption:raw.hydrolysis}}))),
                None=>status.set("No actual recorded anchored 3′ flap proposal exists for that transfer.".into()),
            },Err(e)=>status.set(e.to_string()),
        }
    });
    let apply = Callback::new(move |()| {
        let Some(p) = proposal.get_untracked() else {
            return;
        };
        if frame.get_untracked().as_ref() != Some(&p.frame)
            || form.get_untracked() != p.form
            || locked.get_untracked()
            || busy.get_untracked()
        {
            return;
        }
        add.run(p.action.clone());
        proposal.set(None);
    });
    let advanced = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() || current_capture().is_none() {
            return;
        }
        let f = form.get_untracked();
        let parsed = (|| -> Result<processing::Action, String> {
            if f.hydrolysis.trim().is_empty() || f.hydrolysis.len() > 4096 {
                return Err("Declare an explicit ideal hydrolytic cleavage assumption.".into());
            }
            let number = |value: &str| {
                value
                    .parse::<usize>()
                    .map_err(|_| "Enter valid source-call coordinates.".to_string())
            };
            let one = |value: &str| -> Result<usize, String> {
                number(value)?
                    .checked_sub(1)
                    .ok_or("Call and transfer numbers start at one.".into())
            };
            Ok(processing::Action::CleaveThreePrimeFlap {
                bond: processing::Bond {
                    upstream: physical::Call {
                        strand: number(&f.upstream_strand)?,
                        index: one(&f.upstream_call)?,
                    },
                    downstream: physical::Call {
                        strand: number(&f.downstream_strand)?,
                        index: one(&f.downstream_call)?,
                    },
                },
                displacement_step: one(&f.displacement)?,
                profile: processing::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,
                assumption: f.hydrolysis,
            })
        })();
        match parsed {
            Ok(action) => add.run(action),
            Err(e) => status.set(e),
        }
    });
    let add_nick = Callback::new(move |(expected, nick): (Form, physical::Nick)| {
        if form.get_untracked() != expected || locked.get_untracked() || busy.get_untracked() {
            return;
        }
        if expected.ligase.trim().is_empty() || expected.ligase.len() > 4096 {
            status.set("Declare a ligase assumption (1–4096 bytes) before adding a nick.".into());
            return;
        }
        add.run(processing::Action::SealNick {
            selection: physical::nick_sealing::NickSelection {
                three: nick.three,
                five: nick.five,
            },
            ligase_assumption: expected.ligase,
        });
    });
    let inspect_alive = alive.clone();
    let inspect = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some((key, _)) = current_capture() else {
            return;
        };
        let raw = form.get_untracked();
        let fresh = match capture.run(key.clone()) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let mut request = fresh.request;
        request
            .steps
            .extend(raw.actions.iter().cloned().map(physical::Step::Process));
        let revision = match crate::topo_preparation::task_revision(
            &fresh.inputs,
            fresh.template.as_ref(),
            &collagene_core::topo_run_preparation::Task::Simulate {
                request: request.clone(),
            },
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        cancel.run(());
        inspected.set(None);
        on_reviewed.run(None);
        projected.set(None);
        let executor = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        worker.set_value(Some(executor.clone()));
        busy.set(true);
        let ticket = generation.get_untracked();
        let alive = inspect_alive.clone();
        status
            .set("Reconstructing original preparation and every declared action in order…".into());
        spawn_local(async move {
            let reply = crate::topo_preparation::execute(
                executor.clone(),
                ticket,
                fresh.inputs,
                fresh.template,
                collagene_core::topo_run_preparation::Task::Simulate {
                    request: request.clone(),
                },
            )
            .await;
            executor.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            if locked.get_untracked()
                || frame.get_untracked().as_ref() != Some(&key)
                || form.with_untracked(|f| f.actions != raw.actions)
            {
                return;
            }
            busy.set(false);
            worker.set_value(None);
            let reply = reply.map(|(value, review)| {
                origin_review.set(review);
                value
            });
            match reply {
                Ok(result) => match result {
                    AnalysisResult::TopoSimulation(report) => {
                        match collagene_core::topo_worker::validate_simulation(&report, &revision) {
                            Ok(()) => {
                                let value = Arc::new(Inspected {
                                    frame: key,
                                    form: raw,
                                    request: Arc::new(request),
                                    report: Arc::new(report),
                                });
                                inspected.set(Some(value.clone()));
                                let held = value.clone();
                                let guard_alive = alive.clone();
                                let guard = Callback::new(move |()| -> Result<(), String> {
                                    if !guard_alive.load(Ordering::Relaxed)
                                        || frame.get_untracked().as_ref() != Some(&held.frame)
                                        || form.with_untracked(|f| f.actions != held.form.actions)
                                        || inspected
                                            .get_untracked()
                                            .is_none_or(|current| !Arc::ptr_eq(&current, &held))
                                    {
                                        return Err("The reviewed reaction source, steps or inspection changed".into());
                                    }
                                    Ok(())
                                });
                                on_reviewed.run(Some(Arc::new(ReviewedReaction {
                                    frame: value.frame.clone(),
                                    request: value.request.clone(),
                                    guard,
                                })));
                                status.set("Ordered processing inspected. Review every retained component and choose a contour explicitly.".into());
                            }
                            Err(e) => status.set(e.to_string()),
                        }
                    }
                    AnalysisResult::Error(e) => status.set(e),
                    _ => status.set("Unexpected actions inspection response.".into()),
                },
                Err(e) => status.set(e),
            }
        });
    });
    let project_alive = alive.clone();
    let project = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some(i) = current_inspection() else {
            return;
        };
        let raw = choice.get_untracked();
        let parsed = (|| -> Result<(usize, physical::Call), String> {
            let component = raw
                .component
                .parse::<usize>()
                .map_err(|_| "Choose an actual processed component.")?;
            let strand = raw
                .strand
                .parse::<usize>()
                .map_err(|_| "Choose a prepared source strand.")?;
            let index = raw
                .call
                .parse::<usize>()
                .ok()
                .and_then(|n| n.checked_sub(1))
                .ok_or("Choose a one-based source call.")?;
            if component >= i.report.physical.product.components.len() {
                return Err("The component changed.".into());
            }
            if i.report
                .preparation
                .substrates
                .strands
                .get(strand)
                .is_none_or(|s| index >= s.sequence.len())
            {
                return Err("The source call is outside the preparation.".into());
            }
            Ok((component, physical::Call { strand, index }))
        })();
        let (component, seed) = match parsed {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let fresh = match capture.run(i.frame.clone()) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let request = topo_actions::Request {
            preparation: i.request.as_ref().clone(),
            selection: physical::projection::Selection {
                source_batch_revision: i.report.physical.batch_revision.clone(),
                component,
                reference_seed: seed,
            },
            options: physical::projection::Options {
                require_all_requested_steps_applied: raw.require_all,
                require_known_terminal_chemistry: raw.require_known,
                ..Default::default()
            },
        };
        let revision = match crate::topo_preparation::task_revision(
            &fresh.inputs,
            fresh.template.as_ref(),
            &collagene_core::topo_run_preparation::Task::ProcessActions {
                request: request.clone(),
            },
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        cancel.run(());
        projected.set(None);
        let executor = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        worker.set_value(Some(executor.clone()));
        busy.set(true);
        let ticket = generation.get_untracked();
        let alive = project_alive.clone();
        status.set("Reconstructing actions and proving the explicitly selected contour from actual source calls…".into());
        spawn_local(async move {
            let reply = crate::topo_preparation::execute(
                executor.clone(),
                ticket,
                fresh.inputs,
                fresh.template,
                collagene_core::topo_run_preparation::Task::ProcessActions { request },
            )
            .await;
            executor.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            if locked.get_untracked()
                || choice.get_untracked() != raw
                || !current_inspection().is_some_and(|current| Arc::ptr_eq(&current, &i))
            {
                return;
            }
            busy.set(false);
            worker.set_value(None);
            let reply = reply.map(|(value, review)| {
                origin_review.set(review);
                value
            });
            match reply {
                Ok(result) => match result {
                    AnalysisResult::TopoActionsProduct(report) => {
                        match report.validate_binding(&revision) {
                            Ok(()) => {
                                projected.set(Some(Arc::new(Projected {
                                    inspection: i,
                                    choice: raw,
                                    report: Arc::new(report),
                                })));
                                status.set("Actual action graph and selected contour reviewed. No target acceptance or cloning candidate was created.".into());
                            }
                            Err(e) => status.set(e.to_string()),
                        }
                    }
                    AnalysisResult::Error(e) => status.set(e),
                    _ => status.set("Unexpected processed-actions response.".into()),
                },
                Err(e) => status.set(e),
            }
        });
    });
    let open_alive = alive.clone();
    let open = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some(value) = current_projection() else {
            return;
        };
        if !can_open(&value.report) {
            return;
        }
        let declaration = OpenActions {
            frame: value.inspection.frame.clone(),
            request: value.report.request.clone(),
            expected_revision: value.report.context_revision.clone(),
        };
        // Opening can dispose this view. Only an error may update local state.
        if let Err(error) = on_open.run(declaration) {
            if open_alive.load(Ordering::Relaxed) {
                status.set(error);
            }
        }
    });
    let export_alive = alive.clone();
    let export = Callback::new(move |()| {
        #[derive(serde::Serialize)]
        struct Evidence<'a, T: serde::Serialize> {
            interpretation: &'static str,
            sources: &'a [(String, collagene_core::run_lineage::EventBinding)],
            target: &'a collagene_core::run_lineage::EventBinding,
            epoch: u64,
            inspection_revision: &'a str,
            inspection_generation: u32,
            orientation: usize,
            material_origin: &'a crate::topo_view::OriginEvidence,
            /// Simulation lacks its recipe; projection already retains it in report.request.
            request: Option<&'a topo::Request>,
            report: &'a T,
        }
        let mut output = Json(vec![]);
        let captured_origin_review = origin_review.get_untracked();
        let prepared = (|| -> Result<_, String> {
            if let Some(p) = current_projection() {
                let key = &p.inspection.frame;
                let fresh = capture.run(key.clone())?;
                let origin = crate::topo_view::origin_evidence(
                    fresh.template.as_ref().map(|o| o.origin.selection()),
                    key.graph.clone(),
                    captured_origin_review.as_ref(),
                )?;
                let evidence = Evidence {
                    interpretation: "Source-bound ordered physical actions. Saved origin limitations remain inherited; no isolation, host repair or target acceptance is inferred.",
                    sources: &key.sources,
                    target: &key.target,
                    epoch: key.epoch,
                    inspection_revision: &key.inspection_revision,
                    inspection_generation: key.inspection_generation,
                    orientation: key.orientation,
                    material_origin: &origin,
                    request: None,
                    report: p.report.as_ref(),
                };
                let result = serde_json::to_writer(&mut output, &evidence);
                Ok((
                    result,
                    key.clone(),
                    p.report.context_revision.clone(),
                    true,
                    origin,
                ))
            } else if let Some(i) = current_inspection() {
                let key = &i.frame;
                let fresh = capture.run(key.clone())?;
                let origin = crate::topo_view::origin_evidence(
                    fresh.template.as_ref().map(|o| o.origin.selection()),
                    key.graph.clone(),
                    captured_origin_review.as_ref(),
                )?;
                let evidence = Evidence {
                    interpretation: "Source-bound ordered physical actions. Saved origin limitations remain inherited; no isolation, host repair or target acceptance is inferred.",
                    sources: &key.sources,
                    target: &key.target,
                    epoch: key.epoch,
                    inspection_revision: &key.inspection_revision,
                    inspection_generation: key.inspection_generation,
                    orientation: key.orientation,
                    material_origin: &origin,
                    request: Some(i.request.as_ref()),
                    report: i.report.as_ref(),
                };
                let result = serde_json::to_writer(&mut output, &evidence);
                Ok((
                    result,
                    key.clone(),
                    i.report.batch_revision.clone(),
                    false,
                    origin,
                ))
            } else {
                Err("Inspect the current ordered actions before export.".into())
            }
        })();
        let (result, source_frame, expected_revision, projected, saved_origin) = match prepared {
            Ok(v) => v,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let result = result.map_err(|e| e.to_string()).and_then(|()| {
            crate::app::prepare_download(
                "topo-ordered-actions-evidence.json",
                "application/json",
                &output.0,
            )
        });
        match result {
            Ok(artifact) => {
                let alive = export_alive.clone();
                let check = move || {
                    if !alive.load(Ordering::Relaxed) || status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    let revision = if projected {
                        current_projection().map(|p| p.report.context_revision.clone())
                    } else {
                        current_inspection().map(|i| i.report.batch_revision.clone())
                    };
                    if revision.as_deref() != Some(expected_revision.as_str()) {
                        return Err("Ordered-action evidence changed.".into());
                    }
                    let fresh = capture.run(source_frame.clone())?;
                    if crate::topo_view::origin_evidence(
                        fresh.template.as_ref().map(|o| o.origin.selection()),
                        source_frame.graph.clone(),
                        origin_review.get_untracked().as_ref(),
                    )? != saved_origin
                    {
                        return Err(
                            "The saved material origin or inherited limitations changed.".into(),
                        );
                    }
                    Ok(())
                };
                let completion = check.clone();
                let completion_alive = export_alive.clone();
                artifact.complete_guarded(check, move |result| {
                    if !completion_alive.load(Ordering::Relaxed) || status.try_get_untracked().is_none(){return;}
                    status.set(match result {
                        Ok(()) if completion().is_ok()=>"Captured TOPO reaction evidence saved.".into(),
                        Ok(())=>"Captured TOPO reaction evidence saved; the current source selection has changed.".into(),
                        Err(error) => error,
                    });
                });
            }
            Err(error) => status.set(error),
        }
    });
    view! {<details class="topo-actions"><summary>"Directional flap processing and ordered physical actions"</summary>
        {move||origin_review.get().filter(|r|r.incomplete).map(|r|view!{<aside class="helper"><p>"Retained source limits remain incomplete. Current transfers and processing do not close the source’s search limits or establish experimental isolation."</p><details><summary>"Retained origin findings"</summary><p>{format!("Showing {} of {} findings{}. The retained origin recipe remains available for full replay.",r.findings.len(),r.total,if r.truncated{"; displayed findings are truncated"}else{""})}</p><ul>{r.findings.into_iter().map(|s|view!{<li>{s}</li>}).collect_view()}</ul></details></aside>})}
        <p class="helper">"Declare ideal hydrolytic flap cleavage and nick sealing explicitly. Detached material remains as a free component; target DNA never determines what to remove. These actions do not model a repair kit or establish laboratory success."</p>
        <button class="text-tool" disabled=move||locked.get()||busy.get()||frame.get().is_none() on:click=move|_|load.run(())>"Load original selected path for action review"</button>
        <label>"Ideal hydrolytic cleavage assumption"<textarea rows="2" maxlength="4096" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.hydrolysis.clone()) on:input=move|e|form.update(|f|f.hydrolysis=event_target_value(&e))/></label>
        {move||current_capture().map(|(_,value)|match processing::flap_proposals(&value.physical){
            Ok(proposals)=>view!{<section><h4>"Actual recorded anchored 3′ flaps"</h4>{if proposals.is_empty(){"No eligible recorded anchored 3′ flap proposal."}else{"Each proposal is advisory; review and explicitly add it before native reconstruction."}}{proposals.into_iter().map(|p|{let step=p.displacement_step;view!{<div><p>{format!("Transfer {} displaced {} actual source calls; proposed bond {} → {}. Cleavage retains the displaced fragment in the physical graph.",step+1,p.flap_calls.len(),call_text(p.bond.upstream),call_text(p.bond.downstream))}</p><button class="text-tool" disabled=move||locked.get()||busy.get() on:click=move|_|propose.run(step)>"Review this flap-cut proposal"</button></div>}}).collect_view()}</section>}.into_any(),
            Err(e)=>view!{<p role="alert">{e.to_string()}</p>}.into_any(),
        })}
        {move||proposal.get().filter(|p|frame.get().as_ref()==Some(&p.frame)&&form.get()==p.form).map(|p|view!{<section><h4>"Proposed action — not applied"</h4><p>{action_text(&p.action)}</p><button class="button secondary" disabled=move||locked.get()||busy.get() on:click=move|_|apply.run(())>"Add this exact cut to ordered actions"</button></section>})}
        <details><summary>"Advanced explicit source-call cut declaration"</summary><p class="helper">"Source strands are zero-based; calls and displacement transfers are one-based. Native processing must prove the declared original anchored-flap geometry."</p>
            <label>"Upstream source strand"<input type="number" min="0" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.upstream_strand.clone()) on:input=move|e|form.update(|f|f.upstream_strand=event_target_value(&e))/></label>
            <label>"Upstream source call"<input type="number" min="1" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.upstream_call.clone()) on:input=move|e|form.update(|f|f.upstream_call=event_target_value(&e))/></label>
            <label>"Downstream source strand"<input type="number" min="0" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.downstream_strand.clone()) on:input=move|e|form.update(|f|f.downstream_strand=event_target_value(&e))/></label>
            <label>"Downstream source call"<input type="number" min="1" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.downstream_call.clone()) on:input=move|e|form.update(|f|f.downstream_call=event_target_value(&e))/></label>
            <label>"Recorded displacement transfer"<input type="number" min="1" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.displacement.clone()) on:input=move|e|form.update(|f|f.displacement=event_target_value(&e))/></label>
            <button class="text-tool" disabled=move||locked.get()||busy.get() on:click=move|_|advanced.run(())>"Add explicit cut declaration"</button>
        </details>
        <p class="helper">"Assumption and advanced coordinate fields prepare a new declaration; existing actions change only when explicitly added, removed or reordered."</p>
        <label>"Ligase assumption for a new nick action"<textarea rows="2" maxlength="4096" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.ligase.clone()) on:input=move|e|form.update(|f|f.ligase=event_target_value(&e))/></label>
        {move||current_capture().map(|(_,value)|{
            let expected=form.get();let nicks=current_inspection().map(|i|i.report.physical.product.nicks.clone()).unwrap_or_else(||value.physical.product.nicks.clone());
            view!{<details><summary>{format!("Explicit available nick choices ({})",nicks.len())}</summary><p>"None are added automatically. After a cut, inspect actions again to obtain its actual remaining nicks."</p>{nicks.into_iter().map(|nick|{let raw=expected.clone();let text=format!("{} ({:?}, 3′) → {} ({:?}, 5′)",call_text(nick.three),nick.three_prime,call_text(nick.five),nick.five_prime);view!{<div><p>{text}</p><button class="text-tool" disabled=move||locked.get()||busy.get() on:click=move|_|add_nick.run((raw.clone(),nick.clone()))>"Add this nick to ordered actions"</button></div>}}).collect_view()}</details>}
        })}
        <h4>"Ordered declarations — require native inspection"</h4>
        {move||form.get().actions.into_iter().enumerate().map(|(i,action)|view!{<div><p>{format!("{}. {}",i+1,action_text(&action))}</p><button class="text-tool" disabled=move||locked.get()||busy.get()||i==0 on:click=move|_|form.update(|f|{if i>0&&i<f.actions.len(){f.actions.swap(i,i-1)}})>"Move earlier"</button><button class="text-tool" disabled=move||locked.get()||busy.get()||form.with(|f|i+1>=f.actions.len()) on:click=move|_|form.update(|f|{if i+1<f.actions.len(){f.actions.swap(i,i+1)}})>"Move later"</button><button class="text-tool" disabled=move||locked.get()||busy.get() on:click=move|_|form.update(|f|{if i<f.actions.len(){f.actions.remove(i);}})>"Remove action"</button></div>}).collect_view()}
        <button class="button secondary" disabled=move||locked.get()||busy.get()||current_capture().is_none() on:click=move|_|inspect.run(())>"Inspect these ordered actions"</button>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{cancel.run(());status.set("Actions calculation cancelled.".into());}>"Cancel actions calculation"</button></Show>
        <p role="status">{move||status.get()}</p>
        {move||current_inspection().map(|i|{
            let steps=i.clone();let assumptions=i.clone();
            view!{<section><h4>"Actual ordered reaction graph"</h4><p>{format!("Original source calls: {} · all requested steps applied: {} · unsupported or unknown: {} · work complete: {}",i.report.physical.source_bases,i.report.physical.all_requested_steps_applied,i.report.physical.unsupported_or_unknown,i.report.physical.work_complete)}</p><p>{format!("Reaction bounds: {} steps; {} work cells; {} bytes of work memory; {} bytes of report.",i.report.physical.limits.max_steps,i.report.physical.limits.max_work_cells,i.report.physical.limits.max_work_bytes,i.report.physical.limits.max_report_bytes)}</p><Rows label="Step outcomes and actual chemistry" total=i.report.physical.steps.len() row=Callback::new(move|index:usize|step_text(&steps.report.physical.steps[index]))/><Rows label="Reaction assumptions and limitations" total=i.report.physical.assumptions.len() row=Callback::new(move|index:usize|assumptions.report.physical.assumptions[index].clone())/>{graph_view(Arc::new(i.report.physical.product.clone()))}</section>}

        })}
        {move||current_inspection().map(|i|{
            let strands=i.clone();view!{<section><h4>"Explicit source-derived contour"</h4>
                <label>"Processed component"<select disabled=move||locked.get()||busy.get() prop:value=move||choice.with(|c|c.component.clone()) on:change=move|e|choice.update(|c|c.component=event_target_value(&e))><option value="">"Choose actual processed component…"</option>{i.report.physical.product.components.iter().enumerate().map(|(index,c)|view!{<option value=index.to_string()>{format!("Component {index}: {} paired / {} unpaired calls",c.paired_nucleotides,c.unpaired_nucleotides)}</option>}).collect_view()}</select></label>
                <label>"Prepared source strand"<select disabled=move||locked.get()||busy.get() prop:value=move||choice.with(|c|c.strand.clone()) on:change=move|e|choice.update(|c|{c.strand=event_target_value(&e);c.call.clear();})><option value="">"Choose prepared source strand…"</option>{strands.report.preparation.substrates.strands.iter().enumerate().map(|(index,s)|view!{<option value=index.to_string()>{format!("{index}: {} / {} · {} calls",s.role,s.id,s.sequence.len())}</option>}).collect_view()}</select></label>
                <label>"Source call (one based)"<input type="number" min="1" disabled=move||locked.get()||busy.get() prop:value=move||choice.with(|c|c.call.clone()) on:input=move|e|choice.update(|c|c.call=event_target_value(&e))/></label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get()||busy.get() prop:checked=move||choice.with(|c|c.require_all) on:change=move|e|choice.update(|c|c.require_all=event_target_checked(&e))/>"Require every requested transfer and action to apply"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get()||busy.get() prop:checked=move||choice.with(|c|c.require_known) on:change=move|e|choice.update(|c|c.require_known=event_target_checked(&e))/>"Require known terminal chemistry"</label>
                <button class="button secondary" disabled=move||locked.get()||busy.get() on:click=move|_|project.run(())>"Reconstruct selected action contour"</button>
            </section>}
        })}
        {move||current_projection().map(|p|projection_view(p.report.clone()))}
        <p class="helper">"Opening an advisory plan preserves the original sources, explicit action order and independent target. Verify the plan before creating outputs."</p>
        <button class="button secondary" disabled=move||locked.get()||busy.get()||!current_projection().is_some_and(|p|can_open(&p.report)) on:click=move|_|open.run(())>"Open action candidate for verification"</button>
        <button class="text-tool" disabled=move||current_inspection().is_none()||busy.get() on:click=move|_|export.run(())>"Download full action graph evidence JSON"</button>
    </details>}
}
fn can_open(report: &topo_actions::Report) -> bool {
    matches!(
        &report.projection.outcome,
        physical::projection::Outcome::Proven(_)
    ) && report.duplex.is_some()
        && report.tracked.is_some()
        && report.product.is_some()
}
fn projection_view(report: Arc<topo_actions::Report>) -> AnyView {
    use physical::projection::Outcome;
    let notes = report.clone();
    let outcome=match &report.projection.outcome{
        Outcome::Proven(v)=>view!{<section><h4>"Source-derived paired contour proven"</h4><p>{format!("{:?}: actual top {} / bottom {} bases; {} retained nicks; known terminal chemistry {}. Faithful duplex {}; tracked bridge {}; source-derived reference {}.",v.topology,v.top.len(),v.bottom.len(),v.nicks.len(),v.fully_known_terminal_chemistry,report.duplex.is_some(),report.tracked.is_some(),report.product.is_some())}</p></section>}.into_any(),
        Outcome::Withheld(reason)=>view!{<section><h4>"Action contour withheld"</h4><p>{format!("{reason:?}")}</p></section>}.into_any(),
    };
    view!{<section>{outcome}<p>"All detached strands and other components remain present. This review does not verify a target or establish experimental success. An explicitly opened candidate still requires fresh plan verification."</p><Rows label="Bridge limitations" total=report.bridge_notes.len() row=Callback::new(move|i:usize|notes.bridge_notes[i].clone())/></section>}.into_any()
}
