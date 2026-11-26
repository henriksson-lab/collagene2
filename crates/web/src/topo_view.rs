//! Explicit TA/directional preparation and fresh physical TOPO inspection. Source
//! activation is declared, never inferred from a supplier label or DNA record.
use crate::{
    app::prepare_download as download,
    digest::SourceInfo,
    topo_preparation::{self, Draft},
    worker::CachedWorker,
};
use collagene_core::{
    Alphabet, Molecule,
    protocol::{AnalysisJob, AnalysisResult, Request as WireRequest},
    run_lineage::{EventBinding, InputEventBinding},
    topo::{
        self,
        physical::{self, projection},
    },
    topo_preparation::PreparationMode,
    topo_worker,
    transfer::Candidate,
    verification::{DesignPlan, Operation, PlanInput, Source, TopoSource},
};
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

/// The app resolves these exact document events and guards borrowed source/target
/// models before cloning. Equal names or equal DNA are not event identities.
pub struct CaptureRequest {
    pub sources: Vec<(String, EventBinding)>,
    pub target: EventBinding,
    pub request: topo::Request,
}
pub struct Captured {
    pub inputs: Vec<topo::Input>,
    pub target: Molecule,
}
#[derive(Clone, PartialEq, Eq)]
struct Context {
    draft: Draft,
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
            allow_unapplied_transfers: self.unapplied,
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
    product: Arc<topo::AnnotatedProduct>,
}
struct Checked {
    projected: Arc<Projected>,
    report: Arc<topo::AnnotatedAssessment>,
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
fn json_download(value: &impl Serialize) -> Result<crate::app::ExportArtifact, String> {
    let mut output = BoundedJson(vec![]);
    serde_json::to_writer(&mut output, value).map_err(|e| e.to_string())?;
    download("topo-source-evidence.json", "application/json", &output.0)
}
fn number(text: &str, label: &str) -> Result<usize, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("Choose a valid {label}."))
}
fn physical_text(report: &physical::Report) -> String {
    let mut lines = vec![format!(
        "Work complete: {}; all requested transfers applied: {}; unsupported/unknown: {}; concrete DNA: {}",
        report.work_complete,
        report.all_requested_transfers_applied,
        report.unsupported_or_unknown,
        report.fully_specified_dna
    )];
    for step in &report.steps {
        lines.push(format!(
            "Transfer {}: {:?} → {:?} · {:?} · {} · displacement {:?}",
            step.step + 1,
            step.donor_call,
            step.acceptor_call,
            step.outcome,
            step.reason,
            step.displacement.as_ref().map(|d| d.kind)
        ));
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
fn projection_text(product: &topo::AnnotatedProduct) -> String {
    match &product.physical.projection.outcome {
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
fn AnnotationReview(product: Arc<topo::AnnotatedProduct>) -> impl IntoView {
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
async fn analyze(
    worker: Rc<CachedWorker>,
    id: u32,
    revision: String,
    job: AnalysisJob,
) -> Result<AnalysisResult, String> {
    let result = worker
        .ask(WireRequest {
            id,
            revision: revision.clone(),
            job,
        })
        .await;
    worker.close();
    let response = result?;
    if response.id != id || response.revision != revision {
        return Err("TOPO response has a different input binding.".into());
    }
    match response.result {
        AnalysisResult::Error(error) => Err(error),
        value => Ok(value),
    }
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
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
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
                .ok_or_else(|| "Choose all three current DNA source documents.".to_string())
        };
        d.vector_top_length = find(&d.vector_top)?.length;
        d.vector_bottom_length = find(&d.vector_bottom)?.length;
        let (roles, _) = topo_preparation::compile_request(&d)?;
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
        Ok(Context {
            draft: d,
            sources: bindings,
            target,
            epoch: epoch.get(),
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
        let _ = locked.get();
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
    let nick_locked = Memo::new(move |_| locked.get() || busy.get());
    let nick_frame = Memo::new(move |_| {
        let b = current_bound()?;
        if nick_locked.get() {
            return None;
        }
        let index = orientation.get().parse::<usize>().ok()?;
        b.report.candidate_requests.get(index)?;
        b.report.reference.candidates.get(index)?;
        Some(crate::topo_nick_view::Frame {
            epoch: b.context.epoch,
            sources: b.context.sources.clone(),
            target: b.context.target.binding.clone(),
            inspection_revision: b.report.batch_revision.clone(),
            inspection_generation: b.ticket,
            orientation: index,
        })
    });
    let nick_capture = Callback::new(move |expected: crate::topo_nick_view::Frame| {
        if nick_locked.get_untracked() || nick_frame.get_untracked().as_ref() != Some(&expected) {
            return Err("The inspected source path changed; reload nicks.".into());
        }
        let b = current_bound().ok_or("The exact TOPO source context changed.")?;
        if b.ticket != expected.inspection_generation
            || orientation.get_untracked().parse::<usize>().ok() != Some(expected.orientation)
        {
            return Err("The TOPO orientation changed.".into());
        }
        let request = b
            .report
            .candidate_requests
            .get(expected.orientation)
            .ok_or("No inspected orientation request.")?;
        let physical = b
            .report
            .reference
            .candidates
            .get(expected.orientation)
            .ok_or("No inspected physical graph.")?;
        Ok(crate::topo_nick_view::Capture {
            inputs: b.captured.inputs.clone(),
            request: request.clone(),
            nicks: physical.product.nicks.clone(),
        })
    });
    let actions_capture=Callback::new(move|expected:crate::topo_nick_view::Frame|->Result<crate::topo_actions_view::Capture,String>{
        if nick_locked.get_untracked()||nick_frame.get_untracked().as_ref()!=Some(&expected){
            return Err("The exact original TOPO path or target context changed.".into());
        }
        let b=current_bound().ok_or("The original inspected source context changed.")?;
        let request=b.report.candidate_requests.get(expected.orientation).ok_or("The inspected orientation is unavailable.")?;
        let physical=b.report.reference.candidates.get(expected.orientation).ok_or("The inspected graph is unavailable.")?;
        let fresh=capture.run(CaptureRequest{sources:b.context.sources.clone(),target:b.context.target.binding.clone(),request:request.clone()})?;
        if nick_frame.get_untracked().as_ref()!=Some(&expected){return Err("The source context changed during capture.".into());}
        Ok(crate::topo_actions_view::Capture{inputs:fresh.inputs,request:request.clone(),physical:physical.clone()})
    });
    let open_processed = Callback::new(
        move |open: crate::topo_nick_view::OpenProcessed| -> Result<(), String> {
            if nick_locked.get_untracked()
                || nick_frame.get_untracked().as_ref() != Some(&open.frame)
            {
                return Err(
                    "The exact processed TOPO source path or independent target changed.".into(),
                );
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
                .ok_or("The inspected orientation no longer exists.")?;
            // The root callback independently resolves exact current source events
            // and the active independent target before cloning borrowed models.
            let fresh = capture.run(CaptureRequest {
                sources: b.context.sources.clone(),
                target: b.context.target.binding.clone(),
                request: inspected.clone(),
            })?;
            let revision =
                collagene_core::topo_processed_product::revision(&fresh.inputs, &open.request)
                    .map_err(|e| e.to_string())?;
            if revision != open.expected_revision
                || serde_json::to_vec(inspected).map_err(|e| e.to_string())?
                    != serde_json::to_vec(&open.request.treatment.preparation)
                        .map_err(|e| e.to_string())?
            {
                return Err("The processed declaration differs from the exact inspected source orientation.".into());
            }
            if context.get_untracked().as_ref().ok() != Some(&b.context)
                || nick_frame.get_untracked().as_ref() != Some(&open.frame)
            {
                return Err(
                    "The source or independent target changed during candidate capture.".into(),
                );
            }
            let inputs = fresh
                .inputs
                .into_iter()
                .map(|input| PlanInput {
                    id: input.role,
                    revision: input.molecule.revision(),
                    molecule: input.molecule,
                })
                .collect();
            let sources = b
                .context
                .sources
                .iter()
                .map(|(role, _)| TopoSource {
                    role: role.clone(),
                    source: Source::Input(role.clone()),
                })
                .collect();
            let preferred = b
                .context
                .sources
                .iter()
                .map(|(role, event)| InputEventBinding {
                    input_id: role.clone(),
                    event: event.clone(),
                })
                .collect();
            let proposed = Candidate::new(DesignPlan {
                name: "Processed TOPO source-bound candidate".into(),
                target_revision: fresh.target.revision(),
                inputs,
                steps: vec![Operation::ProcessedTopo {
                    sources,
                    request: open.request,
                }],
            });
            proposed.validate().map_err(|e| e.to_string())?;
            on_candidate.run((proposed, preferred));
            // The callback may dispose this panel. Do not access local signals.
            Ok(())
        },
    );
    let open_actions = Callback::new(
        move |open: crate::topo_actions_view::OpenActions| -> Result<(), String> {
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
                sources: b.context.sources.clone(),
                target: b.context.target.binding.clone(),
                request: inspected.clone(),
            })?;
            let revision = collagene_core::topo_actions::revision(&fresh.inputs, &open.request)
                .map_err(|e| e.to_string())?;
            if revision != open.expected_revision
                || serde_json::to_vec(inspected).map_err(|e| e.to_string())?
                    != serde_json::to_vec(&open.request.preparation).map_err(|e| e.to_string())?
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
            let inputs = fresh
                .inputs
                .into_iter()
                .map(|input| PlanInput {
                    id: input.role,
                    revision: input.molecule.revision(),
                    molecule: input.molecule,
                })
                .collect();
            let sources = b
                .context
                .sources
                .iter()
                .map(|(role, _)| TopoSource {
                    role: role.clone(),
                    source: Source::Input(role.clone()),
                })
                .collect();
            let preferred = b
                .context
                .sources
                .iter()
                .map(|(role, event)| InputEventBinding {
                    input_id: role.clone(),
                    event: event.clone(),
                })
                .collect();
            let proposed = Candidate::new(DesignPlan {
                name: "TOPO ordered physical actions candidate".into(),
                target_revision: fresh.target.revision(),
                inputs,
                steps: vec![Operation::TopoActions {
                    sources,
                    request: open.request,
                }],
            });
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
        let compiled = match topo_preparation::compile_request(&ctx.draft) {
            Ok((_, compiled)) => compiled,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let directional = compiled.directional;
        let vector = compiled.vector_strands;
        let insert = compiled.insert_strands;
        let request = compiled.request;
        let captured = match capture.run(CaptureRequest {
            sources: ctx.sources.clone(),
            target: ctx.target.binding.clone(),
            request: request.clone(),
        }) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Err(error) = topo_worker::source_preflight(&captured.inputs, &request)
            .and_then(|_| topo_worker::target_preflight(&captured.target))
        {
            status.set(error.to_string());
            return;
        }
        let revision = match topo_worker::orientations_revision(
            &captured.inputs,
            &request,
            vector,
            insert,
            directional,
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
            let result = analyze(
                worker,
                ticket,
                revision.clone(),
                AnalysisJob::InspectTopoOrientations {
                    inputs: captured.inputs.clone(),
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
            let batch_revision = topo_worker::simulation_revision(&b.captured.inputs, &request)
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
        let revision = match topo_worker::annotated_projection_revision(
            &b.captured.inputs,
            &request,
            &selection,
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
            let result = analyze(
                worker,
                ticket,
                revision.clone(),
                AnalysisJob::ProjectAnnotatedTopo {
                    inputs: b.captured.inputs.clone(),
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
            match result {
                Ok(AnalysisResult::AnnotatedTopoProjection(value)) => {
                    match topo_worker::validate_annotated_projection(&value, &revision) {
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
        let revision = match topo_worker::annotated_target_revision(
            &p.bound.captured.inputs,
            &p.request,
            &p.selection,
            &p.bound.captured.target,
            &policy,
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
            let result = analyze(
                worker,
                ticket,
                revision.clone(),
                AnalysisJob::CheckAnnotatedTopoTarget {
                    inputs: p.bound.captured.inputs.clone(),
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
            match result {
                Ok(AnalysisResult::AnnotatedTopoTarget(value)) => {
                    match topo_worker::validate_annotated_target(&value, &revision) {
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
        if !c
            .report
            .product
            .physical
            .projection
            .physical
            .fully_specified_dna
        {
            status.set("The full physical graph contains unresolved DNA; a final checked candidate requires fully specified source calls.".into());
            return;
        }
        let b = &c.projected.bound;
        let inputs = b
            .captured
            .inputs
            .iter()
            .map(|input| PlanInput {
                id: input.role.clone(),
                revision: input.molecule.revision(),
                molecule: input.molecule.clone(),
            })
            .collect();
        let sources = b
            .context
            .sources
            .iter()
            .map(|(role, _)| TopoSource {
                role: role.clone(),
                source: Source::Input(role.clone()),
            })
            .collect();
        let preferred = b
            .context
            .sources
            .iter()
            .map(|(role, event)| InputEventBinding {
                input_id: role.clone(),
                event: event.clone(),
            })
            .collect();
        let candidate = Candidate::new(DesignPlan {
            name: "Explicit TA TOPO candidate".into(),
            target_revision: b.captured.target.revision(),
            inputs,
            steps: vec![Operation::CheckedTopo {
                sources,
                request: c.projected.request.clone(),
                selection: c.projected.selection.clone(),
                policy: c.report.policy.clone(),
            }],
        });
        if let Err(error) = candidate.validate() {
            status.set(error.to_string());
            return;
        }
        on_candidate.run((candidate, preferred));
    });
    let export = Callback::new(move |()| {
        let Some(b) = current_bound() else {
            return;
        };
        let p = current_projection();
        let c = current_check();
        #[derive(Serialize)]
        struct Evidence<'a> {
            interpretation: &'static str,
            sources: &'a [(String, EventBinding)],
            target: &'a EventBinding,
            request: &'a topo::Request,
            orientations: &'a topo::OrientationInspection,
            selection: Option<&'a topo::Selection>,
            projection: Option<&'a topo::AnnotatedProduct>,
            assessment: Option<&'a topo::AnnotatedAssessment>,
        }
        let evidence = Evidence {
            interpretation: "Saved source-bound declared-graph evidence. Reconstruct with original source events before reuse; no kit efficacy, physical-target authority, host repair or imported-chemistry restoration.",
            sources: &b.context.sources,
            target: &b.context.target.binding,
            request: &b.request,
            orientations: &b.report,
            selection: p.as_ref().map(|p| &p.selection),
            projection: p.as_ref().map(|p| p.product.as_ref()),
            assessment: c.as_ref().map(|c| c.report.as_ref()),
        };
        crate::app::complete_download(json_download(&evidence), move |result| {
            if status.try_get_untracked().is_none() {
                return;
            }
            if let Err(error) = result {
                status.set(error);
            }
        });
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
    view! {
        <div class="inspector-content topo-panel"><h2>"TOPO cloning"</h2>
            <p class="helper">"Explicit prepared-strand declarations and elementary phosphotyrosyl transfers. DNA records do not prove vector activation. Both insert orientations, reciprocal components, nicks and displaced strands remain visible. Directional flaps can withhold a paired contour; no repair or closed reference DNA is invented."</p>
            <label>"Preparation geometry"<select disabled=move||locked.get()||busy.get() prop:value=move||draft.with(|value|if value.mode==PreparationMode::Ta{"ta"}else{"directional"}) on:change=move|event|draft.update(|value|value.mode=if event_target_value(&event)=="directional"{PreparationMode::DirectionalCacc}else{PreparationMode::Ta})><option value="ta">"TA — explicitly A-tailed PCR"</option><option value="directional">"Directional CACC — actual blunt PCR"</option></select></label>
            <p><strong>"Independent target: "</strong>{move||target.with(|target|format!("{} · {} bases · document {} / event {}",target.name,target.length,target.binding.document_id,target.binding.event_id))}</p>
            <SourcePicker label="TOPO vector top strand source" role="vector_top" draft=draft sources=sources locked=locked/>
            <SourcePicker label="TOPO vector bottom strand source" role="vector_bottom" draft=draft sources=sources locked=locked/>
            <p class="helper">{move||draft.with(|value|if value.mode==PreparationMode::Ta{"Choose both actual vector strands written5′→3′: equal lengths, paired over all but one base, with one explicit3′ overhang. No circle is excised or activated from a label."}else{"Choose both actual vector strands written5′→3′. The bottom strand is four bases longer and begins with stored GGTG (not the reversed drawing GTGG); its paired core begins at offset4. Actual bases and pairing are checked natively."})}</p>
            <label>"Vector preparation assumption"<textarea aria-label="TOPO preparation assumption" rows="3" maxlength="4096" disabled=move||locked.get() prop:value=move||draft.with(|d|d.preparation_assumption.clone()) on:input=move|e|draft.update(|d|d.preparation_assumption=event_target_value(&e))/></label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Declare activated TOPO vector chemistry" disabled=move||locked.get() prop:checked=move||draft.with(|d|d.activation_confirmed) on:change=move|e|draft.update(|d|d.activation_confirmed=event_target_checked(&e))/>"I explicitly declare vector 5′ phosphate and 3′ vaccinia Tyr274-phosphotyrosyl termini; sequence alone does not establish activation."</label>
            <SourcePicker label="TOPO insert PCR template" role="insert_template" draft=draft sources=sources locked=locked/>
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
            <crate::topo_nick_view::TopoNickView frame=nick_frame capture=nick_capture on_open=open_processed locked=nick_locked/>
            <crate::topo_actions_view::TopoActionsView frame=nick_frame capture=actions_capture on_open=open_actions locked=nick_locked/>
            <fieldset class="topo-acceptance"><legend>"Explicit independent-target acceptance policy"</legend><p class="helper">"Every permission starts off. Accepting a retained limitation does not repair it or prove laboratory success; changing policy requires fresh projection and target checking."</p>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.nicks) on:change=move|e|policies.update(|p|p.nicks=event_target_checked(&e))/>"Accept explicitly retained nicks"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.unknown) on:change=move|e|policies.update(|p|p.unknown=event_target_checked(&e))/>"Accept unknown terminal chemistry"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.unapplied) on:change=move|e|policies.update(|p|p.unapplied=event_target_checked(&e))/>"Accept requested transfers that did not apply"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.other) on:change=move|e|policies.update(|p|p.other=event_target_checked(&e))/>"Accept other physical components remaining"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.rotation) on:change=move|e|policies.update(|p|p.rotation=event_target_checked(&e))/>"Allow circular-origin rotation when comparing target"</label>
                <label class="checkbox-label"><input type="checkbox" disabled=move||locked.get() prop:checked=move||policies.with(|p|p.reverse) on:change=move|e|policies.update(|p|p.reverse=event_target_checked(&e))/>"Allow reverse-complement target comparison"</label>
            </fieldset>
            {move||current_projection().map(|p|view!{<section class="topo-projection"><h3>"Fresh selected component"</h3><p>{projection_text(&p.product)}</p><AnnotationReview product=p.product.clone()/><details><summary>"Retained full-graph summary"</summary><pre>{physical_text(&p.product.physical.projection.physical)}</pre></details><button class="button primary full" disabled=move||busy.get()||locked.get() on:click=move|_|check.run(())>"Check independent target"</button></section>})}
            {move||current_check().map(|c|{let accepted=matches!(&c.report.outcome,topo::TargetOutcome::AcceptedUnderDeclaredGraph)&&c.report.product.physical.projection.physical.fully_specified_dna;view!{<section class="topo-target"><h3>"Fresh independent-target assessment"</h3><p>{format!("{:?}",c.report.outcome)}</p><p class="helper">"Agreement applies only to this source/event/recipe/component and explicit policy. A final checked TOPO plan can verify as consistent under declared graph assumptions when the full DNA is specified. This does not predict experimental success, repair retained nicks or flaps, validate annotation biology, or promote the graph into downstream generic chemistry."</p><label class="checkbox-label"><input type="checkbox" prop:checked=move||reviewed.get() on:change=move|e|reviewed.set(event_target_checked(&e))/>"I reviewed the preparation assumptions, both paths, all components and nicks, and accept this explicit policy for final target checking."</label><button class="button secondary full" disabled=move||busy.get()||locked.get()||!accepted||!reviewed.get() on:click=move|_|candidate.run(())>"Review checked TOPO candidate in Plan"</button></section>}})}
        </div>
    }
}
