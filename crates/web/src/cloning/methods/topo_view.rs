//! Explicit TA/directional preparation and fresh physical TOPO inspection. Source
//! activation is declared, never inferred from a supplier label or DNA record.
use crate::cloning::methods::topo_preparation;
use crate::cloning::methods::topo_preparation::CaptureRequest;
use crate::cloning::methods::topo_preparation::Captured;
use crate::cloning::methods::topo_preparation::Draft;
use crate::io::download::prepare_download;
use crate::platform::worker::CachedWorker;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    cloning::topo::{
        self,
        physical::{self, projection},
        preparation::PreparationMode,
        worker as topo_worker,
    },
    documents::transfer::Candidate,
    material::lineage::{EventBinding, InputEventBinding},
    protocol::schema::AnalysisResult,
};
use collagene_model::Alphabet;
use leptos::prelude::*;
use serde::Serialize;
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
struct Context {
    draft: Draft,
    template_origin: Option<collagene_core::material::run::MaterialSelection>,
    sources: Vec<(String, EventBinding)>,
    target: SourceInfo,
    epoch: u64,
}
#[derive(Clone, Copy, Default, PartialEq, Eq)]
struct Policies {
    nicks: bool,
    unknown: bool,
    unapplied: bool,
    other: bool,
    rotation: bool,
    reverse: bool,
}
impl Policies {
    fn policy(self) -> topo::AcceptancePolicy {
        topo::AcceptancePolicy {
            allow_nicks: self.nicks,
            allow_unknown_terminal_chemistry: self.unknown,
            allow_unapplied_steps: self.unapplied,
            allow_other_components: self.other,
            allow_circular_rotation: self.rotation,
            allow_reverse_complement: self.reverse,
        }
    }
}
#[derive(Clone, PartialEq, Eq)]
struct Choice {
    orientation: String,
    component: String,
    strand: String,
    call: String,
    policies: Policies,
}
struct Bound {
    ticket: u32,
    context: Context,
    captured: Arc<Captured>,
    request: topo::Request,
    report: Arc<topo::OrientationInspection>,
}
struct Projected {
    bound: Arc<Bound>,
    choice: Choice,
    request: topo::Request,
    selection: topo::Selection,
    product: Arc<topo::Product>,
}
struct Checked {
    projected: Arc<Projected>,
    report: Arc<topo::TargetAssessment>,
}
/// Documentary origin caveats accompany every exported inner TOPO report.
#[derive(Clone, PartialEq, Eq, Serialize)]
pub(crate) struct OriginEvidence {
    template: Option<collagene_core::material::run::MaterialSelection>,
    graph: Option<collagene_core::cloning::topo::graph_history::Selection>,
    review: Option<OriginReviewEvidence>,
}
#[derive(Clone, PartialEq, Eq, Serialize)]
struct OriginReviewEvidence {
    incomplete: bool,
    findings: Vec<String>,
    findings_total: usize,
    findings_truncated: bool,
}
pub(crate) fn origin_evidence(
    selection: Option<collagene_core::material::run::MaterialSelection>,
    graph: Option<collagene_core::cloning::topo::graph_history::Selection>,
    review: Option<&topo_preparation::OriginReview>,
) -> Result<OriginEvidence, String> {
    if (selection.is_some() || graph.is_some()) != review.is_some() {
        return Err("Reinspect the saved material before exporting its origin limitations.".into());
    }
    Ok(OriginEvidence {
        template: selection,
        graph,
        review: review.map(|review| OriginReviewEvidence {
            incomplete: review.incomplete,
            findings: review.findings.clone(),
            findings_total: review.total,
            findings_truncated: review.truncated,
        }),
    })
}
struct BoundedJson(Vec<u8>);
impl Write for BoundedJson {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > 16 * 1024 * 1024)
        {
            return Err(io::Error::other("TOPO evidence exceeds16MiB"));
        }
        self.0
            .try_reserve(bytes.len())
            .map_err(|_| io::Error::other("TOPO evidence allocation"))?;
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn json_download(value: &impl Serialize) -> Result<crate::io::download::ExportArtifact, String> {
    let mut output = BoundedJson(vec![]);
    serde_json::to_writer(&mut output, value).map_err(|e| e.to_string())?;
    prepare_download("topo-source-evidence.json", "application/json", &output.0)
}
fn number(text: &str, label: &str) -> Result<usize, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("Choose a valid {label}."))
}
fn physical_text(report: &physical::Report) -> String {
    let mut lines = vec![format!(
        "Work complete: {}; all requested reaction steps applied: {}; unsupported/unknown: {}; concrete DNA: {}",
        report.work_complete,
        report.all_requested_steps_applied,
        report.unsupported_or_unknown,
        report.fully_specified_dna
    )];
    for step in &report.steps {
        lines.push(crate::cloning::methods::topo_actions_view::step_text(step));
    }
    for (i, component) in report.product.components.iter().enumerate() {
        lines.push(format!(
            "Component {i}: chains {:?}; paired calls {}; unpaired calls {}; covalently closed {}",
            component.chains,
            component.paired_nucleotides,
            component.unpaired_nucleotides,
            component.covalently_closed
        ));
    }
    lines.push(format!("Open nicks: {}; unpaired runs/flaps/free strands: {}; released activation identities: {:?}",report.product.nicks.len(),report.product.unpaired.len(),report.product.released_activations));
    for nick in report.product.nicks.iter().take(16) {
        lines.push(format!(
            "Nick {:?} ({:?}) → {:?} ({:?})",
            nick.three, nick.three_prime, nick.five, nick.five_prime
        ));
    }
    for run in report.product.unpaired.iter().take(16) {
        lines.push(format!(
            "Unpaired chain {}: {}..{}; 5′terminal {}; 3′terminal {}; displaced by {:?}",
            run.chain, run.start, run.end, run.terminal_five, run.terminal_three, run.displaced_by
        ));
    }
    let text = lines.join("\n");
    if text.chars().count() > 4000 {
        format!(
            "{}\nPreview clipped; download bounded full evidence.",
            text.chars().take(4000).collect::<String>()
        )
    } else {
        text
    }
}
fn projection_text(product: &topo::Product) -> String {
    match &product.projection.outcome {
        projection::Outcome::Proven(view) => format!(
            "Paired contour: {:?}; top {} bases / actual bottom {} bases; source seed {:?}; seed reference position {}; {} retained nicks; ends {:?}; known terminal chemistry {}. A circular contour is not covalent closure or host repair.",
            view.topology,
            view.top.len(),
            view.bottom.len(),
            view.origin,
            view.seed_reference_position,
            view.nicks.len(),
            view.ends,
            view.fully_known_terminal_chemistry
        ),
        projection::Outcome::Withheld(reason) => format!(
            "Contour withheld: {reason:?}. Full graph and all reciprocal components remain available; no target-derived DNA or repair is substituted."
        ),
    }
}
#[component]
fn AnnotationReview(product: Arc<topo::Product>) -> impl IntoView {
    let count = product
        .reference
        .as_ref()
        .map_or(0, |reference| reference.features.len());
    let reviewed = product.reference.as_ref().map_or(0, |reference| {
        reference
            .features
            .iter()
            .filter(|feature| feature.needs_review)
            .count()
    });
    view! {<section class="topo-annotation-review"><h3>"Source annotation correspondence"</h3>
        <p class="helper">"Annotations follow the selected top-reference source calls only. Bottom pairing is not a second source copy. Projected coordinates do not establish inherited chemistry or biological validity; unmapped originals remain in the evidence."</p>
        <p>{format!("{count} projected reference annotations; {reviewed} marked for review.")}</p>
        {product.annotations.as_ref().map(|evidence|view!{<div>
            {evidence.warnings.iter().take(16).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
            {evidence.sources.iter().map(|source|view!{<details><summary>{format!("{} · supported {} · projected {} · dropped {} · retained unmapped {}",source.role,source.projection_supported,source.projected,source.dropped,source.unmapped.len())}</summary><code>{source.source_revision.clone()}</code>
                <p class="helper">"Dropped counts refer to coordinate projection; the same original can remain in retained unmapped evidence. Candidate correspondence maps must not be applied when projection is unsupported."</p>
                {source.warnings.iter().take(16).map(|warning|view!{<p>{warning.clone()}</p>}).collect_view()}
                {source.unmapped.iter().take(8).map(|item|view!{<p>{format!("Original feature {} · {} · {}",item.source_index,item.feature.label().chars().take(80).collect::<String>(),item.reason.chars().take(400).collect::<String>())}</p>}).collect_view()}
                <p class="helper">{format!("Showing at most8 retained originals and16 warnings here; download full bounded evidence for all{} originals and exact locations/qualifiers.",source.unmapped.len())}</p>
            </details>}).collect_view()}
        </div>})}
        {product.annotations.is_none().then(||view!{<p>"No annotation projection is released because no paired reference contour is available."</p>})}
    </section>}
}
#[component]
fn SourcePicker(
    label: &'static str,
    role: &'static str,
    draft: RwSignal<Draft>,
    sources: Memo<Vec<SourceInfo>>,
    locked: Memo<bool>,
) -> impl IntoView {
    view! {<label>{label}<select aria-label=label disabled=move||locked.get() prop:value=move||draft.with(|d|match role{"vector_top"=>d.vector_top.clone(),"vector_bottom"=>d.vector_bottom.clone(),_=>d.insert_template.clone()}) on:change=move|event|{let value=event_target_value(&event);draft.update(|d|match role{"vector_top"=>d.vector_top=value,"vector_bottom"=>d.vector_bottom=value,_=>d.insert_template=value});}><option value="">"Choose an actual 5′→3′ DNA source…"</option>{move||sources.get().into_iter().filter(|s|s.alphabet==Alphabet::Dna).map(|s|view!{<option value=s.binding.document_id>{format!("{} · {} bases",s.name,s.length)}</option>}).collect_view()}</select></label>}
}

#[component]
pub fn TopoPanel(
    target: Memo<SourceInfo>,
    sources: Memo<Vec<SourceInfo>>,
    capture: Callback<CaptureRequest, Result<Captured, String>>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    on_candidate: Callback<(Candidate, Vec<InputEventBinding>)>,
    get_workspace: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    validate: Callback<collagene_core::documents::workspace::Workspace, Result<(), String>>,
    on_restore: Callback<collagene_core::documents::workspace::Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
    let origin_review = RwSignal::new(None::<topo_preparation::OriginReview>);
    let template_origin = RwSignal::new(None::<crate::history::run_origin_model::PinnedOrigin>);
    let template_source = Memo::new(move |_| {
        let id = draft.with(|d| d.insert_template.clone());
        sources.with(|s| s.iter().find(|s| s.binding.document_id == id).cloned())
    });
    let selected_origin = Memo::new(move |_| template_origin.get());
    let on_origin = Callback::new(
        move |pin: Option<crate::history::run_origin_model::PinnedOrigin>| {
            if let Some(pin) = &pin {
                pin.as_ref()
                    .validate_for_workspace(&get_workspace.run(())?, pin.event())
                    .map_err(|e| e.to_string())?;
                if template_source
                    .get_untracked()
                    .is_none_or(|s| s.binding != *pin.event())
                {
                    return Err("The insert template event changed.".into());
                }
            }
            template_origin.set(pin);
            Ok(())
        },
    );
    let cacc_proposal = RwSignal::new(None::<(Context, String)>);
    let orientation = RwSignal::new(String::new());
    let component = RwSignal::new(String::new());
    let strand = RwSignal::new(String::new());
    let call = RwSignal::new(String::new());
    let policies = RwSignal::new(Policies::default());
    let reviewed = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let bound = RwSignal::new(None::<Arc<Bound>>);
    let projected = RwSignal::new(None::<Arc<Projected>>);
    let checked = RwSignal::new(None::<Arc<Checked>>);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let context = Memo::new(move |_| {
        let mut d = draft.get();
        let target = target.get();
        let all = sources.get();
        if target.alphabet != Alphabet::Dna {
            return Err("Choose an independent DNA target document.".to_string());
        }
        let find = |id: &str| {
            all.iter()
                .find(|s| s.binding.document_id == id && s.alphabet == Alphabet::Dna)
                .ok_or_else(|| "Choose every requested current DNA source document.".to_string())
        };
        if d.source_mode == topo_preparation::SourceMode::Original {
            d.vector_top_length = find(&d.vector_top)?.length;
            d.vector_bottom_length = find(&d.vector_bottom)?.length;
        } else {
            if d.graph.is_none() {
                return Err("Choose a saved graph source".into());
            }
            if d.vector_strands.iter().any(String::is_empty)
                || (d.source_mode == topo_preparation::SourceMode::ResumeGraph
                    && d.insert_strands.iter().any(String::is_empty))
            {
                return Err("Inspect the saved graph and choose existing strand handles".into());
            }
        }
        let roles = topo_preparation::external_sources(&d);
        let mut bytes = 0usize;
        let mut bindings = vec![];
        for (role, id) in roles {
            let source = find(&id)?;
            bytes = bytes
                .checked_add(source.length)
                .filter(|&n| n <= 100_000)
                .ok_or("Selected sources exceed the100,000-base preparation limit.")?;
            bindings.push((role, source.binding.clone()));
        }
        if target.length > 100_000 {
            return Err("The independent target exceeds the100,000-base check limit.".into());
        }
        if let Some(pin) = template_origin
            .get()
            .filter(|_| d.source_mode != topo_preparation::SourceMode::ResumeGraph)
        {
            if bindings
                .iter()
                .find(|(r, _)| r == "insert_template")
                .is_none_or(|(_, e)| *e != *pin.event())
            {
                return Err("Select the saved origin again for this exact template event.".into());
            }
        }
        Ok(Context {
            template_origin: if d.source_mode == topo_preparation::SourceMode::ResumeGraph {
                None
            } else {
                template_origin.get().map(|p| p.selection())
            },
            sources: bindings,
            target,
            epoch: epoch.get(),
            draft: d,
        })
    });
    let choice = Memo::new(move |_| Choice {
        orientation: orientation.get(),
        component: component.get(),
        strand: strand.get(),
        call: call.get(),
        policies: policies.get(),
    });
    let cancel = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = context.get();
        origin_review.set(None);
        cancel.run(());
        bound.set(None);
        projected.set(None);
        checked.set(None);
        reviewed.set(false);
        cacc_proposal.set(None);
        orientation.set(String::new());
        component.set(String::new());
        strand.set(String::new());
        call.set(String::new());
    });
    Effect::new(move || {
        if locked.get() {
            cancel.run(());
            projected.set(None);
            checked.set(None);
            reviewed.set(false);
        }
    });
    Effect::new(move |_| {
        let _ = choice.get();
        cancel.run(());
        projected.set(None);
        checked.set(None);
        reviewed.set(false);
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
    });
    let start = move || -> Result<(Rc<CachedWorker>, u32), String> {
        cancel.run(());
        let worker = Rc::new(CachedWorker::new()?);
        active.set_value(Some(worker.clone()));
        busy.set(true);
        Ok((worker, generation.get_untracked()))
    };
    let current_bound = move || {
        bound
            .get()
            .filter(|b| context.with(|c| c.as_ref().ok() == Some(&b.context)) && !locked.get())
    };
    let reviewed_reaction =
        RwSignal::new(None::<Arc<crate::cloning::methods::topo_actions_view::ReviewedReaction>>);
    let nick_locked = Memo::new(move |_| locked.get() || busy.get());
    let nick_frame =
        Memo::new(move |_| {
            let b = bound
                .get()
                .filter(|b| context.with(|c| c.as_ref().ok() == Some(&b.context)))?;
            let index = orientation.get().parse::<usize>().ok()?;
            b.report.candidate_requests.get(index)?;
            b.report.reference.candidates.get(index)?;
            Some(crate::cloning::methods::topo_actions_view::Frame {
                epoch: b.context.epoch,
                sources: b.context.sources.clone(),
                target: b.context.target.binding.clone(),
                inspection_revision: b.report.batch_revision.clone(),
                inspection_generation: b.ticket,
                orientation: index,
                graph: b.context.draft.graph.clone().filter(|_| {
                    b.context.draft.source_mode != topo_preparation::SourceMode::Original
                }),
            })
        });
    let actions_capture=Callback::new(move|expected:crate::cloning::methods::topo_actions_view::Frame|->Result<crate::cloning::methods::topo_actions_view::Capture,String>{
        if nick_locked.get_untracked()||nick_frame.get_untracked().as_ref()!=Some(&expected){
            return Err("The exact original TOPO path or target context changed.".into());
        }
        let b=current_bound().ok_or("The original inspected source context changed.")?;
        let request=b.report.candidate_requests.get(expected.orientation).ok_or("The inspected orientation is unavailable.")?;
        let physical=b.report.reference.candidates.get(expected.orientation).ok_or("The inspected graph is unavailable.")?;
        let fresh=capture.run(CaptureRequest{draft:b.context.draft.clone(),sources:b.context.sources.clone(),template_origin:b.context.template_origin.clone(),target:b.context.target.binding.clone(),request:Some(request.clone())})?;
        if nick_frame.get_untracked().as_ref()!=Some(&expected){return Err("The source context changed during capture.".into());}
        Ok(crate::cloning::methods::topo_actions_view::Capture{template:fresh.template,inputs:fresh.inputs,request:request.clone(),physical:physical.clone()})
    });
    let open_actions = Callback::new(
        move |open: crate::cloning::methods::topo_actions_view::OpenActions| -> Result<(), String> {
            if nick_locked.get_untracked()
                || nick_frame.get_untracked().as_ref() != Some(&open.frame)
            {
                return Err("The exact action source path or independent target changed.".into());
            }
            let b = current_bound().ok_or("The original TOPO context is unavailable.")?;
            if b.ticket != open.frame.inspection_generation
                || b.context.epoch != open.frame.epoch
                || b.context.target.binding != open.frame.target
                || b.context.sources != open.frame.sources
                || orientation.get_untracked().parse::<usize>().ok() != Some(open.frame.orientation)
            {
                return Err("The inspected source events, target or orientation changed.".into());
            }
            let inspected = b
                .report
                .candidate_requests
                .get(open.frame.orientation)
                .ok_or("The exact inspected orientation is unavailable.")?;
            let fresh = capture.run(CaptureRequest {
                draft: b.context.draft.clone(),
                sources: b.context.sources.clone(),
                template_origin: b.context.template_origin.clone(),
                target: b.context.target.binding.clone(),
                request: Some(inspected.clone()),
            })?;
            let revision = topo_preparation::task_revision(
                &fresh.inputs,
                fresh.template.as_ref(),
                &collagene_core::cloning::topo::run_preparation::Task::ProcessActions {
                    request: open.request.clone(),
                },
            )
            .map_err(|e| e.to_string())?;
            if revision != open.expected_revision
                || topo_preparation::validate_action_append(inspected, &open.request.preparation)
                    .is_err()
            {
                return Err("The ordered action declaration differs from the exact inspected source orientation.".into());
            }
            if context.get_untracked().as_ref().ok() != Some(&b.context)
                || nick_frame.get_untracked().as_ref() != Some(&open.frame)
            {
                return Err(
                    "The source or independent target changed during action candidate capture."
                        .into(),
                );
            }
            let preferred = b
                .context
                .sources
                .iter()
                .map(|(role, event)| InputEventBinding {
                    input_id: role.clone(),
                    event: event.clone(),
                })
                .collect();
            let proposed_plan = collagene_core::cloning::topo::run_preparation::compile(
                &fresh.inputs,
                &fresh.target,
                fresh.template.as_ref(),
                &collagene_core::cloning::topo::run_preparation::Task::ProcessActions {
                    request: open.request,
                },
            )
            .map_err(|e| e.to_string())?;
            let proposed = Candidate::new(proposed_plan);
            proposed.validate().map_err(|e| e.to_string())?;
            on_candidate.run((proposed, preferred));
            // The root callback can dispose this panel.
            Ok(())
        },
    );
    let current_projection = move || {
        projected.get().filter(|p| {
            context.with(|c| c.as_ref().ok() == Some(&p.bound.context))
                && choice.with(|c| *c == p.choice)
                && !locked.get()
        })
    };
    let current_check = move || {
        checked.get().filter(|c| {
            context.with(|ctx| ctx.as_ref().ok() == Some(&c.projected.bound.context))
                && choice.with(|choice| *choice == c.projected.choice)
                && !locked.get()
        })
    };
    let inspect_alive = alive.clone();
    let inspect = Callback::new(move |()| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let ctx = match context.get_untracked() {
            Ok(ctx) => ctx,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let captured = match capture.run(CaptureRequest {
            draft: ctx.draft.clone(),
            sources: ctx.sources.clone(),
            template_origin: ctx.template_origin.clone(),
            target: ctx.target.binding.clone(),
            request: None,
        }) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let directional = captured.preparation.directional;
        let vector = captured.preparation.vector_strands;
        let insert = captured.preparation.insert_strands;
        let request = captured.preparation.request.clone();
        if let Err(error) = topo_worker::target_preflight(&captured.target) {
            status.set(error.to_string());
            return;
        }
        let revision = match topo_preparation::task_revision(
            &captured.inputs,
            captured.template.as_ref(),
            &collagene_core::cloning::topo::run_preparation::Task::Orientations {
                request: request.clone(),
                vector,
                insert,
                directional,
            },
        ) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let captured = Arc::new(captured);
        let (worker, ticket) = match start() {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        bound.set(None);
        projected.set(None);
        checked.set(None);
        reviewed.set(false);
        status.set(if directional {"Preparing actual blunt PCR and both directional CACC orientations; flap or contour withholding remains explicit."}else{"Preparing explicit strands and inspecting both declared TA orientations…"}.into());
        let alive = inspect_alive.clone();
        spawn_local(async move {
            let result = topo_preparation::execute(
                worker,
                ticket,
                captured.inputs.clone(),
                captured.template.clone(),
                collagene_core::cloning::topo::run_preparation::Task::Orientations {
                    request: request.clone(),
                    vector,
                    insert,
                    directional,
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            active.set_value(None);
            if context.get_untracked().as_ref().ok() != Some(&ctx) {
                return;
            }
            let result = result.map(|(value, review)| {
                origin_review.set(review);
                value
            });
            match result {
                Ok(AnalysisResult::TopoOrientations(report)) => {
                    match topo_worker::validate_orientations(&report, &revision) {
                        Ok(()) => {
                            bound.set(Some(Arc::new(Bound {
                                ticket,
                                context: ctx,
                                captured,
                                request,
                                report: Arc::new(report),
                            })));
                            status.set("Both declared paths inspected. Choose the request, physical component and source seed explicitly.".into());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Unexpected TOPO orientation response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let project_alive = alive.clone();
    let project = Callback::new(move |()| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let Some(b) = current_bound() else {
            return;
        };
        let chosen = choice.get_untracked();
        let selected = (|| -> Result<(topo::Request, topo::Selection), String> {
            let which = number(&chosen.orientation, "orientation request")?;
            let request = b
                .report
                .candidate_requests
                .get(which)
                .cloned()
                .ok_or("Choose an inspected orientation request.")?;
            let index = number(&chosen.component, "component")?;
            let physical = b
                .report
                .reference
                .candidates
                .get(which)
                .ok_or("Orientation report changed.")?;
            if index >= physical.product.components.len() {
                return Err("Choose a physical component.".into());
            }
            let seed_strand = number(&chosen.strand, "source strand")?;
            let seed_call = number(&chosen.call, "one-based source call")?
                .checked_sub(1)
                .ok_or("Source call numbering starts at1.")?;
            if b.report
                .preparation
                .substrates
                .strands
                .get(seed_strand)
                .is_none_or(|s| seed_call >= s.sequence.len())
            {
                return Err("The source seed is outside its original prepared strand.".into());
            }
            let batch_revision = topo_preparation::batch_revision(
                &b.captured.inputs,
                b.captured.template.as_ref(),
                &request,
            )
            .map_err(|e| e.to_string())?;
            Ok((
                request,
                topo::Selection {
                    batch_revision,
                    component: index,
                    reference_seed: physical::Call {
                        strand: seed_strand,
                        index: seed_call,
                    },
                },
            ))
        })();
        let (request, selection) = match selected {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let revision = match topo_preparation::task_revision(
            &b.captured.inputs,
            b.captured.template.as_ref(),
            &collagene_core::cloning::topo::run_preparation::Task::Project {
                request: request.clone(),
                selection: selection.clone(),
            },
        ) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let (worker, ticket) = match start() {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        projected.set(None);
        checked.set(None);
        reviewed.set(false);
        status.set(
            "Freshly reconstructing the selected physical component and paired contour…".into(),
        );
        let alive = project_alive.clone();
        spawn_local(async move {
            let result = topo_preparation::execute(
                worker,
                ticket,
                b.captured.inputs.clone(),
                b.captured.template.clone(),
                collagene_core::cloning::topo::run_preparation::Task::Project {
                    request: request.clone(),
                    selection: selection.clone(),
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            active.set_value(None);
            if context.get_untracked().as_ref().ok() != Some(&b.context)
                || choice.get_untracked() != chosen
            {
                return;
            }
            let result = result.map(|(value, review)| {
                origin_review.set(review);
                value
            });
            match result {
                Ok(AnalysisResult::TopoProjection(value)) => {
                    match topo_worker::validate_projection(&value, &revision) {
                        Ok(()) => {
                            projected.set(Some(Arc::new(Projected {
                                bound: b,
                                choice: chosen,
                                request,
                                selection,
                                product: Arc::new(value),
                            })));
                            status.set(String::new());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Unexpected TOPO projection response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let target_alive = alive.clone();
    let check = Callback::new(move |()| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let Some(p) = current_projection() else {
            return;
        };
        let policy = p.choice.policies.policy();
        let revision = match topo_preparation::task_revision(
            &p.bound.captured.inputs,
            p.bound.captured.template.as_ref(),
            &collagene_core::cloning::topo::run_preparation::Task::Check {
                request: p.request.clone(),
                selection: p.selection.clone(),
                target: p.bound.captured.target.clone(),
                policy: policy.clone(),
            },
        ) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let (worker, ticket) = match start() {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        checked.set(None);
        reviewed.set(false);
        status.set("Reconstructing again and checking the independent target under only the selected acceptance policies…".into());
        let alive = target_alive.clone();
        spawn_local(async move {
            let result = topo_preparation::execute(
                worker,
                ticket,
                p.bound.captured.inputs.clone(),
                p.bound.captured.template.clone(),
                collagene_core::cloning::topo::run_preparation::Task::Check {
                    request: p.request.clone(),
                    selection: p.selection.clone(),
                    target: p.bound.captured.target.clone(),
                    policy,
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            active.set_value(None);
            if context.get_untracked().as_ref().ok() != Some(&p.bound.context)
                || choice.get_untracked() != p.choice
            {
                return;
            }
            let result = result.map(|(value, review)| {
                origin_review.set(review);
                value
            });
            match result {
                Ok(AnalysisResult::TopoTarget(value)) => {
                    match topo_worker::validate_target(&value, &revision) {
                        Ok(()) => {
                            checked.set(Some(Arc::new(Checked {
                                projected: p,
                                report: Arc::new(value),
                            })));
                            status.set(String::new());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Unexpected TOPO target response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let candidate = Callback::new(move |()| {
        if busy.get_untracked() || locked.get_untracked() || !reviewed.get_untracked() {
            return;
        }
        let Some(c) = current_check() else {
            return;
        };
        if !matches!(
            &c.report.outcome,
            topo::TargetOutcome::AcceptedUnderDeclaredGraph
        ) {
            return;
        }
        if !c.report.product.projection.physical.fully_specified_dna {
            status.set("The full physical graph contains unresolved DNA; a final checked candidate requires fully specified source calls.".into());
            return;
        }
        let b = &c.projected.bound;
        let result = (|| -> Result<_, String> {
            let fresh = capture.run(CaptureRequest {
                draft: b.context.draft.clone(),
                sources: b.context.sources.clone(),
                target: b.context.target.binding.clone(),
                template_origin: b.context.template_origin.clone(),
                request: Some(c.projected.request.clone()),
            })?;
            let task = collagene_core::cloning::topo::run_preparation::Task::Check {
                request: c.projected.request.clone(),
                selection: c.projected.selection.clone(),
                target: fresh.target.clone(),
                policy: c.report.policy.clone(),
            };
            let plan = collagene_core::cloning::topo::run_preparation::compile(
                &fresh.inputs,
                &fresh.target,
                fresh.template.as_ref(),
                &task,
            )
            .map_err(|e| e.to_string())?;
            if context.get_untracked().as_ref().ok() != Some(&b.context)
                || choice.get_untracked() != c.projected.choice
            {
                return Err(
                    "The graph, source, target or checked policy changed before candidate capture"
                        .into(),
                );
            }
            let candidate = Candidate::new(plan);
            candidate.validate().map_err(|e| e.to_string())?;
            let preferred = b
                .context
                .sources
                .iter()
                .map(|(role, event)| InputEventBinding {
                    input_id: role.clone(),
                    event: event.clone(),
                })
                .collect();
            Ok((candidate, preferred))
        })();
        match result {
            Ok(value) => on_candidate.run(value),
            Err(error) => status.set(error),
        }
    });
    let export_alive = alive.clone();
    let export = Callback::new(move |()| {
        let Some(b) = current_bound() else {
            return;
        };
        let p = current_projection();
        let c = current_check();
        let expected_choice = choice.get_untracked();
        let saved_origin = match origin_evidence(
            b.context.template_origin.clone(),
            topo_preparation::active_graph(&b.context.draft),
            origin_review.get_untracked().as_ref(),
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        #[derive(Serialize)]
        struct Evidence<'a> {
            interpretation: &'static str,
            sources: &'a [(String, EventBinding)],
            material_origin: &'a OriginEvidence,
            target: &'a EventBinding,
            request: &'a topo::Request,
            orientations: &'a topo::OrientationInspection,
            selection: Option<&'a topo::Selection>,
            projection: Option<&'a topo::Product>,
            assessment: Option<&'a topo::TargetAssessment>,
        }
        let evidence = Evidence {
            interpretation: "Saved source-bound declared-graph evidence. Reconstruct with original source events before reuse; no kit efficacy, physical-target authority, host repair or imported-chemistry restoration.",
            sources: &b.context.sources,
            material_origin: &saved_origin,
            target: &b.context.target.binding,
            request: &b.request,
            orientations: &b.report,
            selection: p.as_ref().map(|p| &p.selection),
            projection: p.as_ref().map(|p| p.product.as_ref()),
            assessment: c.as_ref().map(|c| c.report.as_ref()),
        };
        match json_download(&evidence) {
            Ok(artifact) => {
                let alive = export_alive.clone();
                let check = move || -> Result<(), String> {
                    if !alive.load(Ordering::Relaxed) || status.try_get_untracked().is_none() {
                        return Err("TOPO export review closed.".into());
                    }
                    if !current_bound().is_some_and(|current| Arc::ptr_eq(&current, &b))
                        || choice.get_untracked() != expected_choice
                    {
                        return Err(
                            "The TOPO sources, material origin, target or selected policy changed."
                                .into(),
                        );
                    }
                    let same_projection = match (current_projection(), &p) {
                        (None, None) => true,
                        (Some(current), Some(expected)) => Arc::ptr_eq(&current, expected),
                        _ => false,
                    };
                    let same_check = match (current_check(), &c) {
                        (None, None) => true,
                        (Some(current), Some(expected)) => Arc::ptr_eq(&current, expected),
                        _ => false,
                    };
                    if !same_projection
                        || !same_check
                        || origin_evidence(
                            b.context.template_origin.clone(),
                            topo_preparation::active_graph(&b.context.draft),
                            origin_review.get_untracked().as_ref(),
                        )? != saved_origin
                    {
                        return Err(
                            "The selected TOPO evidence or inherited limitations changed.".into(),
                        );
                    }
                    capture
                        .run(CaptureRequest {
                            draft: b.context.draft.clone(),
                            sources: b.context.sources.clone(),
                            target: b.context.target.binding.clone(),
                            template_origin: b.context.template_origin.clone(),
                            request: Some(b.request.clone()),
                        })
                        .map(|_| ())
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() {
                        return;
                    }
                    status.set(match result {
                        Ok(()) => "TOPO evidence saved.".into(),
                        Err(error) => error,
                    });
                });
            }
            Err(error) => status.set(error),
        }
    });
    let propose_cacc = move |_| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let result = (|| -> Result<(Context, String), String> {
            let current = context.get_untracked()?;
            if current.draft.mode != PreparationMode::DirectionalCacc {
                return Err(
                    "Choose directional preparation before proposing a CACC prefix.".into(),
                );
            }
            let proposed = topo_preparation::cacc_prefix_proposal(&current.draft)?;
            Ok((current, proposed.forward_sequence))
        })();
        match result {
            Ok(proposal) => {
                cacc_proposal.set(Some(proposal));
                status.set("CACC prefix proposed only. Review the exact ordered forward oligo and explicitly Apply; PCR and transfer eligibility remain unchecked.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let apply_cacc = move |_| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let result = (|| -> Result<String, String> {
            let (before, ordered) = cacc_proposal
                .get_untracked()
                .ok_or("Propose a CACC prefix first.")?;
            if context.get_untracked()? != before {
                return Err(
                    "Source events or primer declarations changed; propose the prefix again."
                        .into(),
                );
            }
            Ok(ordered)
        })();
        match result {
            Ok(ordered) => {
                draft.update(|value| value.forward_sequence = ordered);
                cacc_proposal.set(None);
                status.set("The explicit CACC prefix was applied to the forward oligo. Reinspect actual PCR and both orientations; the independent target was not changed.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let graph_initial =
        Memo::new(move |_| draft.with(|d| d.source_mode == topo_preparation::SourceMode::Original));
    let graph_suffix = Memo::new(move |_| {
        let current = bound.get()?;
        if context.with(|c| c.as_ref().ok() != Some(&current.context)) {
            return None;
        }
        let selected = orientation.get().parse::<usize>().ok()?;
        Some(current.report.candidate_requests.get(selected)?.steps.len())
    });
    let graph_reviewed = Memo::new(move |_| {
        reviewed_reaction
            .get()
            .is_some_and(|r| nick_frame.get().as_ref() == Some(&r.frame) && r.guard.run(()).is_ok())
    });
    let graph_locked = Memo::new(move |_| locked.get() || busy.get());
    let graph_capture=Callback::new(move|(kind,name):(crate::cloning::methods::topo_graph_model::CaptureKind,String)|->Result<crate::cloning::methods::topo_graph_view::CaptureCandidate,String>{
        if locked.get_untracked()||busy.get_untracked(){return Err("Wait for current TOPO work to finish".into());}
        match kind {
            crate::cloning::methods::topo_graph_model::CaptureKind::ReviewedReaction=>{
                let reviewed=reviewed_reaction.get_untracked().ok_or("Inspect the current full reaction before saving")?;
                reviewed.guard.run(())?;
                let b=bound.get_untracked().ok_or("The exact source inspection is unavailable")?;
                if context.get_untracked().as_ref().ok()!=Some(&b.context)||nick_frame.get_untracked().as_ref()!=Some(&reviewed.frame){return Err("The reviewed reaction source context changed".into());}
                let original=b.report.candidate_requests.get(reviewed.frame.orientation).ok_or("The inspected orientation changed")?;
                topo_preparation::validate_action_append(original,reviewed.request.as_ref())?;
                let request=reviewed.request.as_ref().clone();
                let fresh=capture.run(CaptureRequest{draft:b.context.draft.clone(),sources:b.context.sources.clone(),target:b.context.target.binding.clone(),request:Some(request.clone()),template_origin:b.context.template_origin.clone()})?;
                let bindings=b.context.sources.iter().map(|(role,event)|collagene_core::cloning::topo::graph_history::SourceBinding{role:role.clone(),event:event.clone(),origin:fresh.template.as_ref().filter(|template|template.role==*role).map(|template|template.origin.clone())}).collect();
                let expected=b.context.clone();let ticket=b.ticket;let chosen=reviewed.frame.orientation;
                let guard=Callback::new(move|()|->Result<(),String>{
                    reviewed.guard.run(())?;
                    if context.get_untracked().as_ref().ok()!=Some(&expected)||orientation.get_untracked().parse::<usize>().ok()!=Some(chosen)||bound.get_untracked().is_none_or(|b|b.ticket!=ticket){return Err("The source events, reviewed transcript or independent target changed before saving".into());}
                    Ok(())
                });
                Ok(crate::cloning::methods::topo_graph_view::CaptureCandidate{request:collagene_core::cloning::topo::graph_history::CaptureRequest{name,bindings,request},guard})
            }
            crate::cloning::methods::topo_graph_model::CaptureKind::InitialVector=>{
                let before=draft.get_untracked();let captured_epoch=epoch.get_untracked();let frame=get_workspace.run(())?;
                let request=crate::cloning::methods::topo_graph_model::initial_request(&frame,&before,name)?;
                let expected=request.bindings.iter().map(|binding|(binding.role.clone(),binding.event.clone())).collect::<Vec<_>>();
                let guard=Callback::new(move|()|->Result<(),String>{if epoch.get_untracked()!=captured_epoch||draft.get_untracked()!=before{return Err("The vector declaration changed before saving".into());}if !sources.with_untracked(|all|expected.iter().all(|(_,event)|all.iter().any(|source|source.binding==*event))){return Err("A vector source event changed before saving".into());}Ok(())});
                Ok(crate::cloning::methods::topo_graph_view::CaptureCandidate{request,guard})
            }
            crate::cloning::methods::topo_graph_model::CaptureKind::Reaction{pending_steps}=>{
                let b=bound.get_untracked().ok_or("Inspect and choose an elementary orientation first")?;
                if context.get_untracked().as_ref().ok()!=Some(&b.context){return Err("The exact graph/source context changed".into());}
                let selected=orientation.get_untracked().parse::<usize>().map_err(|_|"Choose an inspected orientation")?;
                let mut request=b.report.candidate_requests.get(selected).ok_or("The inspected orientation is unavailable")?.clone();
                if pending_steps>request.steps.len(){return Err("The selected reaction prefix exceeds the declared transcript".into());}
                request.steps.truncate(pending_steps);
                let fresh=capture.run(CaptureRequest{draft:b.context.draft.clone(),sources:b.context.sources.clone(),target:b.context.target.binding.clone(),request:Some(request.clone()),template_origin:b.context.template_origin.clone()})?;
                let bindings=b.context.sources.iter().map(|(role,event)|collagene_core::cloning::topo::graph_history::SourceBinding{role:role.clone(),event:event.clone(),origin:fresh.template.as_ref().filter(|template|template.role==*role).map(|template|template.origin.clone())}).collect();
                let expected=b.context.clone();let ticket=b.ticket;
                let guard=Callback::new(move|()|->Result<(),String>{if context.get_untracked().as_ref().ok()!=Some(&expected)||orientation.get_untracked().parse::<usize>().ok()!=Some(selected)||bound.get_untracked().is_none_or(|b|b.ticket!=ticket){return Err("The graph, source events or elementary orientation changed before saving".into());}Ok(())});
                Ok(crate::cloning::methods::topo_graph_view::CaptureCandidate{request:collagene_core::cloning::topo::graph_history::CaptureRequest{name,bindings,request},guard})
            }
        }
    });
    view! {
        <div class="inspector-content topo-panel"><h2>"TOPO cloning"</h2>
            <p class="helper">"Explicit prepared-strand declarations and elementary phosphotyrosyl transfers. DNA records do not prove vector activation. Both insert orientations, reciprocal components, nicks and displaced strands remain visible. Directional flaps can withhold a paired contour; no repair or closed reference DNA is invented."</p>
            <label>"Preparation geometry"<select disabled=move||locked.get()||busy.get() prop:value=move||draft.with(|value|if value.mode==PreparationMode::Ta{"ta"}else{"directional"}) on:change=move|event|draft.update(|value|value.mode=if event_target_value(&event)=="directional"{PreparationMode::DirectionalCacc}else{PreparationMode::Ta})><option value="ta">"TA — explicitly A-tailed PCR"</option><option value="directional">"Directional CACC — actual blunt PCR"</option></select></label>
            <p><strong>"Independent target: "</strong>{move||target.with(|target|format!("{} · {} bases · document {} / event {}",target.name,target.length,target.binding.document_id,target.binding.event_id))}</p>
            <label>"Source mode"<select disabled=move||locked.get() prop:value=move||draft.with(|d|match d.source_mode{topo_preparation::SourceMode::Original=>"original",topo_preparation::SourceMode::SavedVector=>"saved-vector",topo_preparation::SourceMode::ResumeGraph=>"resume"}) on:change=move|e|{let mode=match event_target_value(&e).as_str(){"original"=>topo_preparation::SourceMode::Original,"saved-vector"=>topo_preparation::SourceMode::SavedVector,"resume"=>topo_preparation::SourceMode::ResumeGraph,_=>return};draft.update(|d|d.source_mode=mode);}>
            <option value="original">"Original vector declaration + new PCR insert"</option><option value="saved-vector">"Saved initial vector graph + new PCR insert"</option><option value="resume">"Continue saved reaction graph · no new PCR"</option></select></label>
            <Show when=move||draft.with(|d|d.source_mode!=topo_preparation::SourceMode::Original)><crate::cloning::methods::topo_graph_view::GraphSourcePicker draft epoch locked snapshot=get_workspace/></Show>
            <Show when=move||draft.with(|d|d.source_mode==topo_preparation::SourceMode::ResumeGraph)><p class="helper">"Continue the retained reaction using its existing insert strands. No new PCR, insert copy or activation is added."</p></Show>
            <Show when=move||draft.with(|d|d.source_mode==topo_preparation::SourceMode::Original)><details open><summary>"Original activated-vector declaration"</summary>
            <SourcePicker label="TOPO vector top strand source" role="vector_top" draft=draft sources=sources locked=locked/>
            <SourcePicker label="TOPO vector bottom strand source" role="vector_bottom" draft=draft sources=sources locked=locked/>
            <p class="helper">{move||draft.with(|value|if value.mode==PreparationMode::Ta{"Choose both actual vector strands written5′→3′: equal lengths, paired over all but one base, with one explicit3′ overhang. No circle is excised or activated from a label."}else{"Choose both actual vector strands written5′→3′. The bottom strand is four bases longer and begins with stored GGTG (not the reversed drawing GTGG); its paired core begins at offset4. Actual bases and pairing are checked natively."})}</p>
            <label>"Vector preparation assumption"<textarea aria-label="TOPO preparation assumption" rows="3" maxlength="4096" disabled=move||locked.get() prop:value=move||draft.with(|d|d.preparation_assumption.clone()) on:input=move|e|draft.update(|d|d.preparation_assumption=event_target_value(&e))/></label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Declare activated TOPO vector chemistry" disabled=move||locked.get() prop:checked=move||draft.with(|d|d.activation_confirmed) on:change=move|e|draft.update(|d|d.activation_confirmed=event_target_checked(&e))/>"I explicitly declare vector 5′ phosphate and 3′ vaccinia Tyr274-phosphotyrosyl termini; sequence alone does not establish activation."</label>
            </details></Show>
            <Show when=move||draft.with(|d|d.source_mode!=topo_preparation::SourceMode::ResumeGraph)><details open><summary>"New PCR insert"</summary>
            <SourcePicker label="TOPO insert PCR template" role="insert_template" draft=draft sources=sources locked=locked/>
            {move||origin_review.get().filter(|r|r.incomplete).map(|r|view!{<aside class="helper"><p>"The saved source's originating check is incomplete. This inspection reconstructs actual material under those retained assumptions; it does not close origin search limits or establish a uniquely isolated experimental product. Final plan verification remains incomplete."</p><details><summary>"Retained origin findings"</summary><p>{format!("Showing {} of {} findings{}. The retained origin recipe remains available for full replay.",r.findings.len(),r.total,if r.truncated{"; displayed findings are truncated"}else{""})}</p><ul>{r.findings.into_iter().map(|s|view!{<li>{s}</li>}).collect_view()}</ul></details></aside>})}
            <crate::history::run_origin_view::SourceOriginPicker source=template_source epoch get_workspace selected=selected_origin on_select=on_origin/>
            <Show when=move||template_origin.get().is_some()><p class="helper">"The saved insert material is reconstructed before new PCR. Vector graph history remains separate. Unsupported staggered or nicked templates are reported without inserting repair steps."</p></Show>
            <label>"Forward ordered oligo 5′→3′"<input aria-label="TOPO forward ordered primer" maxlength="1024" disabled=move||locked.get() prop:value=move||draft.with(|d|d.forward_sequence.clone()) on:input=move|e|draft.update(|d|d.forward_sequence=event_target_value(&e))/></label>
            <label>"Forward 3′ annealing length"<input type="number" min="1" max="1024" aria-label="TOPO forward annealing length" disabled=move||locked.get() prop:value=move||draft.with(|d|if d.forward_annealing==0{String::new()}else{d.forward_annealing.to_string()}) on:input=move|e|draft.update(|d|d.forward_annealing=event_target_value(&e).parse().unwrap_or(0))/></label>
            <label>"Reverse ordered oligo 5′→3′"<input aria-label="TOPO reverse ordered primer" maxlength="1024" disabled=move||locked.get() prop:value=move||draft.with(|d|d.reverse_sequence.clone()) on:input=move|e|draft.update(|d|d.reverse_sequence=event_target_value(&e))/></label>
            <label>"Reverse 3′ annealing length"<input type="number" min="1" max="1024" aria-label="TOPO reverse annealing length" disabled=move||locked.get() prop:value=move||draft.with(|d|if d.reverse_annealing==0{String::new()}else{d.reverse_annealing.to_string()}) on:input=move|e|draft.update(|d|d.reverse_annealing=event_target_value(&e).parse().unwrap_or(0))/></label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Declare single A tail on both PCR strands" disabled=move||locked.get() prop:checked=move||draft.with(|d|d.add_a_tail) on:change=move|e|draft.update(|d|d.add_a_tail=event_target_checked(&e))/>"Explicitly synthesize one A on each PCR strand's 3′ end after exact PCR."</label>
            <p class="helper">"The declared annealing suffix is split from each whole ordered oligo; the remaining5′ prefix stays literal. TA requires an explicitly selected A addition. Directional preparation requires no terminal addition and actual PCR top DNA beginning CACC; no CACC or phosphorylation is inferred. Exact binding and product alternatives are checked natively."</p>
            <Show when=move||draft.with(|value|value.mode==PreparationMode::DirectionalCacc)>
                <button class="text-tool" disabled=move||locked.get()||busy.get()||context.get().is_err() on:click=propose_cacc>"Propose literal CACC forward-primer prefix"</button>
                {move||cacc_proposal.get().map(|(_,ordered)|view!{<section><h4>"Forward oligo prefix proposal — not applied"</h4><pre>{ordered}</pre><p class="helper">"Four new literal bases will be supplied by this ordered primer. The existing annealing length and independent target stay unchanged; native PCR and transfer checks may still reject the design."</p><button class="button secondary" disabled=move||locked.get()||busy.get() on:click=apply_cacc>"Apply this exact forward oligo"</button></section>})}
            </Show>
            </details></Show>
            <crate::cloning::methods::topo_graph_view::GraphCaptureControls capture=graph_capture initial=graph_initial reviewed=graph_reviewed pending_steps=graph_suffix locked=graph_locked snapshot=get_workspace validate on_restore epoch db workspace_busy save_generation message/>
            <button class="button primary full" disabled=move||busy.get()||locked.get()||context.get().is_err() on:click=move|_|inspect.run(())>{move||draft.with(|value|if value.mode==PreparationMode::Ta{"Inspect both TA orientations"}else{"Inspect both directional orientations"})}</button>
            <Show when=move||busy.get()><button class="button secondary full" on:click=move|_|{cancel.run(());status.set("TOPO calculation cancelled; no late result will be used.".into());}>"Cancel TOPO"</button></Show>
            <p class="helper" role="status">{move||if status.get().is_empty(){context.get().err().unwrap_or_default()}else{status.get()}}</p>
            {move||current_bound().map(|b|view!{
                <section class="topo-orientations"><h3>"Both declared paths"</h3><p>{format!("Orientation enumeration complete: {}. Exhaustive biochemical mixture: {}.",b.report.reference.orientation_search_complete,b.report.reference.biochemical_mixture_exhaustive)}</p>
                    {b.report.reference.candidates.iter().enumerate().map(|(index,report)|view!{<details><summary>{format!("Path {} · {} components · {} nicks · {} unpaired runs",index+1,report.product.components.len(),report.product.nicks.len(),report.product.unpaired.len())}</summary><pre>{physical_text(report)}</pre></details>}).collect_view()}
                    <label>"Explicit orientation request"<select aria-label="TOPO orientation request" prop:value=move||orientation.get() on:change=move|e|{orientation.set(event_target_value(&e));component.set(String::new());strand.set(String::new());call.set(String::new());}><option value="">"Choose an inspected path…"</option>{b.report.candidate_requests.iter().enumerate().map(|(index,_)|view!{<option value=index.to_string()>{format!("Path {} · insert strand {} joins top vector",index+1,if index==0{2}else{3})}</option>}).collect_view()}</select></label>
                    <label>"Physical component (zero based)"<select aria-label="TOPO component" prop:value=move||component.get() on:change=move|e|component.set(event_target_value(&e))><option value="">"Choose a component…"</option>{move||current_bound().and_then(|b|orientation.get().parse::<usize>().ok().and_then(|index|b.report.reference.candidates.get(index).map(|report|report.product.components.iter().enumerate().map(|(i,c)|view!{<option value=i.to_string()>{format!("Component {i}: {} paired / {} unpaired calls",c.paired_nucleotides,c.unpaired_nucleotides)}</option>}).collect_view()))).into_view()}</select></label>
                    <label>"Reference seed source strand (zero based)"<select aria-label="TOPO source seed strand" prop:value=move||strand.get() on:change=move|e|strand.set(event_target_value(&e))><option value="">"Choose a prepared source strand…"</option>{b.report.preparation.substrates.strands.iter().enumerate().map(|(i,s)|view!{<option value=i.to_string()>{format!("{i}: {} / {} · {} calls",s.role,s.id,s.sequence.len())}</option>}).collect_view()}</select></label>
                    <label>"Reference seed call (one based)"<input type="number" min="1" aria-label="TOPO source seed call" prop:value=move||call.get() on:input=move|e|call.set(event_target_value(&e))/></label>
                    <p class="helper">"The source seed chooses actual strand direction and, for a circle, the reference origin. It must belong to the selected paired component. No component is presumed to be the experimental winner."</p>
                    <button class="button secondary full" disabled=move||busy.get()||locked.get()||orientation.with(String::is_empty)||component.with(String::is_empty)||strand.with(String::is_empty)||call.with(String::is_empty) on:click=move|_|project.run(())>"Project selected physical component"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|export.run(())>"Download TOPO evidence JSON"</button>
                </section>
            })}

            <crate::cloning::methods::topo_actions_view::TopoActionsView frame=nick_frame capture=actions_capture on_open=open_actions on_reviewed=Callback::new(move|value|reviewed_reaction.set(value)) locked=nick_locked/>
            <fieldset class="topo-acceptance"><legend>"Explicit independent-target acceptance policy"</legend><p class="helper">"Every permission starts off. Accepting a retained limitation does not repair it or prove laboratory success; changing policy requires fresh projection and target checking."</p>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.nicks) on:change=move|e|policies.update(|p|p.nicks=event_target_checked(&e))/>"Accept explicitly retained nicks"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.unknown) on:change=move|e|policies.update(|p|p.unknown=event_target_checked(&e))/>"Accept unknown terminal chemistry"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.unapplied) on:change=move|e|policies.update(|p|p.unapplied=event_target_checked(&e))/>"Accept requested reaction steps that did not apply"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.other) on:change=move|e|policies.update(|p|p.other=event_target_checked(&e))/>"Accept other physical components remaining"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.rotation) on:change=move|e|policies.update(|p|p.rotation=event_target_checked(&e))/>"Allow circular-origin rotation when comparing target"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.reverse) on:change=move|e|policies.update(|p|p.reverse=event_target_checked(&e))/>"Allow reverse-complement target comparison"</label>
            </fieldset>
            {move||current_projection().map(|p|view!{<section class="topo-projection"><h3>"Fresh selected component"</h3><p>{projection_text(&p.product)}</p><AnnotationReview product=p.product.clone()/><details><summary>"Retained full-graph summary"</summary><pre>{physical_text(&p.product.projection.physical)}</pre></details><button class="button primary full" disabled=move||busy.get()||locked.get() on:click=move|_|check.run(())>"Check independent target"</button></section>})}
            {move||current_check().map(|c|{let accepted=matches!(&c.report.outcome,topo::TargetOutcome::AcceptedUnderDeclaredGraph)&&c.report.product.projection.physical.fully_specified_dna;view!{<section class="topo-target"><h3>"Fresh independent-target assessment"</h3><p>{format!("{:?}",c.report.outcome)}</p><p class="helper">"Agreement applies only to this source/event/recipe/component and explicit policy. A final checked TOPO plan can verify as consistent under declared graph assumptions when the full DNA is specified. This does not predict experimental success, repair retained nicks or flaps, validate annotation biology, or promote the graph into downstream generic chemistry."</p><label class="checkbox-label"><input type="checkbox" prop:checked=move||reviewed.get() on:change=move|e|reviewed.set(event_target_checked(&e))/>"I reviewed the preparation assumptions, both paths, all components and nicks, and accept this explicit policy for final target checking."</label><button class="button secondary full" disabled=move||busy.get()||locked.get()||!accepted||!reviewed.get() on:click=move|_|candidate.run(())>"Review checked TOPO candidate in Plan"</button></section>}})}
        </div>
    }
}
