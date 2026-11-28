//! Collapsed vector-graph controls beneath the shared construct-QC declarations.
use super::model::{BindingDraft, Capture, Form, SpacingDraft};
use crate::{
    biology::construct_qc_model::Draft as QcDraft, io::download::ExportArtifact,
    platform::worker::CachedWorker, sources::identity::SourceInfo,
};
use collagene_core::{
    construct_qc::vector_class as vector,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

fn binding_text(
    form: RwSignal<Form>,
    key: u64,
    label: &'static str,
    get: fn(&BindingDraft) -> &str,
    set: fn(&mut BindingDraft, String),
    locked: RwSignal<bool>,
) -> impl IntoView {
    view! {
        <label>{label}<input maxlength="128" disabled=move ||locked.get()
            prop:value=move ||form.with(|form|form.bindings.iter().find(|row|row.key==key).map(|row|get(row).to_owned()).unwrap_or_default())
            on:input=move |event|form.update(|form|if let Some(row)=form.bindings.iter_mut().find(|row|row.key==key){set(row,event_target_value(&event));})/></label>
    }
}
fn binding_row(
    form: RwSignal<Form>,
    parent: RwSignal<QcDraft>,
    rows: Memo<Vec<(usize, String, bool)>>,
    binding: BindingDraft,
    locked: RwSignal<bool>,
) -> impl IntoView {
    let key = binding.key;
    view! {
        <article class="qc-vector-binding">
            <strong>{binding.role}</strong>
            {binding_text(form,key,"Binding ID",|row|&row.id,|row,value|row.id=value,locked)}
            <label>"Source annotation"
                <select disabled=move ||locked.get()
                    prop:value=move ||form.with(|form|form.bindings.iter().find(|row|row.key==key).map(|row|row.feature.clone()).unwrap_or_default())
                    on:change=move |event|form.update(|form|if let Some(row)=form.bindings.iter_mut().find(|row|row.key==key){row.feature=event_target_value(&event);})>
                    <option value="">"No binding"</option>
                    {move ||rows.with(|rows|rows.iter().map(|(index,label,review)|view!{
                        <option value=index.to_string() disabled=*review>{label.clone()}</option>
                    }).collect_view())}
                </select>
            </label>
            <label>"Existing pinned reference"
                <select disabled=move ||locked.get()
                    prop:value=move ||form.with(|form|form.bindings.iter().find(|row|row.key==key).map(|row|row.reference.clone()).unwrap_or_default())
                    on:change=move |event|form.update(|form|if let Some(row)=form.bindings.iter_mut().find(|row|row.key==key){row.reference=event_target_value(&event);})>
                    <option value="">"Choose pinned reference"</option>
                    {move ||parent.with(|parent|parent.options.references.iter().take(64).map(|reference|view!{
                        <option value=reference.id.clone()>{format!("{} · {} bases",reference.id,reference.sequence.len())}</option>
                    }).collect_view())}
                </select>
            </label>
            <p class="helper">{move ||{
                let selected = form.with(|form|form.bindings.iter().find(|row|row.key==key).map(|row|(row.feature.clone(),row.reference.clone())));
                match selected {
                    Some((feature,reference)) if !feature.is_empty() => {
                        let assumed = feature.parse::<usize>().ok().is_some_and(|index|parent.with(|parent|parent.assumed.contains(&index)));
                        let available = parent.with(|parent|parent.options.references.iter().any(|candidate|candidate.id==reference));
                        if !assumed {"This annotation has no current role assumption. Use the shared annotation controls above; the role remains indeterminate."}
                        else if !available {"The selected reference is unavailable. Choose a current pinned reference above; no identity is inferred."}
                        else {"Uses your current annotation assumption and exact pinned reference. Biological role is not certified."}
                    }
                    _ => "Leave both annotation and reference unset for an undeclared role. Missing annotation is not biological absence.",
                }
            }}</p>
            <button class="text-tool" disabled=move ||locked.get()
                on:click=move |_|form.update(|form|form.bindings.retain(|row|row.key!=key))>"Remove binding row"</button>
        </article>
    }
}
fn binding_choice(
    form: RwSignal<Form>,
    label: &'static str,
    get: fn(&Form) -> &str,
    set: fn(&mut Form, String),
    locked: RwSignal<bool>,
) -> impl IntoView {
    view! {
        <label>{label}<select disabled=move ||locked.get()
            prop:value=move ||form.with(|form|get(form).to_owned())
            on:change=move |event|form.update(|form|set(form,event_target_value(&event)))>
            <option value="">"Undeclared"</option>
            {move ||form.with(|form|form.active_bindings().map(|row|view!{
                <option value=row.id.clone()>{format!("{} · {}",row.id,row.role)}</option>
            }).collect_view())}
        </select></label>
    }
}
fn spacing_text(
    form: RwSignal<Form>,
    key: u64,
    label: &'static str,
    get: fn(&SpacingDraft) -> &str,
    set: fn(&mut SpacingDraft, String),
    numeric: bool,
    locked: RwSignal<bool>,
) -> impl IntoView {
    view! {
        <label>{label}<input maxlength=if numeric {"20"} else {"1024"}
            inputmode=if numeric {"numeric"} else {"text"} disabled=move ||locked.get()
            prop:value=move ||form.with(|form|form.spacing.iter().find(|row|row.key==key).map(|row|get(row).to_owned()).unwrap_or_default())
            on:input=move |event|form.update(|form|if let Some(row)=form.spacing.iter_mut().find(|row|row.key==key){set(row,event_target_value(&event));})/></label>
    }
}
fn spacing_row(form: RwSignal<Form>, row: SpacingDraft, locked: RwSignal<bool>) -> impl IntoView {
    let key = row.key;
    let choice = move |label: &'static str, upstream: bool| {
        view! {
            <label>{label}<select disabled=move ||locked.get()
                prop:value=move ||form.with(|form|form.spacing.iter().find(|row|row.key==key).map(|row|if upstream {row.upstream.clone()} else {row.downstream.clone()}).unwrap_or_default())
                on:change=move |event|form.update(|form|if let Some(row)=form.spacing.iter_mut().find(|row|row.key==key){if upstream {row.upstream=event_target_value(&event);}else{row.downstream=event_target_value(&event);}})>
                <option value="">"Choose binding"</option>
                {move ||form.with(|form|form.active_bindings().map(|row|view!{
                    <option value=row.id.clone()>{row.id.clone()}</option>
                }).collect_view())}
            </select></label>
        }
    };
    view! {
        <article class="qc-vector-spacing">
            <strong>{row.id}</strong>
            {choice("Upstream binding",true)}{choice("Downstream binding",false)}
            {spacing_text(form,key,"Minimum gap (bases)",|row|&row.minimum,|row,value|row.minimum=value,true,locked)}
            {spacing_text(form,key,"Maximum gap (optional)",|row|&row.maximum,|row,value|row.maximum=value,true,locked)}
            {spacing_text(form,key,"Reason for this constraint",|row|&row.rationale,|row,value|row.rationale=value,false,locked)}
            <button class="text-tool" disabled=move ||locked.get()
                on:click=move |_|form.update(|form|form.spacing.retain(|row|row.key!=key))>"Remove spacing constraint"</button>
        </article>
    }
}

#[component]
pub(crate) fn VectorClassPanel(
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    draft: RwSignal<QcDraft>,
    rows: Memo<Vec<(usize, String, bool)>>,
    on_select: Callback<Vec<(usize, usize)>>,
    workspace_busy: RwSignal<bool>,
) -> impl IntoView {
    let profiles = StoredValue::new(vector::catalogue().map_err(|error| error.to_string()));
    let form = RwSignal::new(Form::default());
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<vector::Report>)>);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        source.get();
        epoch.get();
        form.set(Form::default());
    });
    Effect::new(move |_| {
        let _ = (
            source.get(),
            epoch.get(),
            draft.get(),
            form.get(),
            workspace_busy.get(),
        );
        generation.update(|token| *token = token.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        page.set(0);
        status.set("Declare a profile and exact graph bindings, then run the checks.".into());
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = worker.get_value() {
                worker.close();
            }
        });
    }
    let guard = Callback::new(move |capture: Capture| -> Result<(), String> {
        if workspace_busy.try_get_untracked() != Some(false) {
            return Err("The workspace is changing; try again afterward.".into());
        }
        capture.guard(
            &source
                .try_get_untracked()
                .ok_or("The vector-check panel was closed.")?,
            epoch
                .try_get_untracked()
                .ok_or("The vector-check panel was closed.")?,
            &draft
                .try_get_untracked()
                .ok_or("The vector-check panel was closed.")?,
            &form
                .try_get_untracked()
                .ok_or("The vector-check panel was closed.")?,
        )?;
        Ok(())
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            let molecule = get_molecule.run(());
            let prepared = Capture::prepare(
                source.get_untracked(),
                epoch.get_untracked(),
                draft.get_untracked(),
                form.get_untracked(),
                &molecule,
            );
            let capture = match prepared {
                Ok(capture) => capture,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            if let Err(error) = guard.run(capture.clone()) {
                status.set(error);
                return;
            }
            let current_worker = match CachedWorker::new() {
                Ok(worker) => Rc::new(worker),
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            worker.set_value(Some(current_worker.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            report.set(None);
            status.set("Checking the declared vector graph…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let response = current_worker
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::AssessVectorClass {
                            molecule,
                            request: capture.request.clone(),
                        },
                    })
                    .await;
                current_worker.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let checked = (|| {
                    guard.run(capture.clone())?;
                    match response?.result {
                        AnalysisResult::VectorClass(value) => {
                            capture.validate_report(&value)?;
                            Ok(value)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected vector-class response.".into()),
                    }
                })();
                match checked {
                    Ok(value) => {
                        report.set(Some((capture, Arc::new(value))));
                        status.set("Declared graph findings ready; review assumptions and unresolved premises.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let export_alive = alive.clone();
    let export = move |_| {
        let Some((capture, value)) = report.get_untracked() else {
            return;
        };
        let prepared = (|| {
            guard.run(capture.clone())?;
            let bytes =
                serde_json::to_vec_pretty(value.as_ref()).map_err(|error| error.to_string())?;
            ExportArtifact::prepare("vector-class-qc.json", "application/json", &bytes)
        })();
        match prepared {
            Ok(artifact) => {
                let alive = export_alive.clone();
                artifact.complete_guarded(
                    move || {
                        if !alive.load(Ordering::Relaxed) {
                            return Err("The vector-check panel was closed.".into());
                        }
                        guard.run(capture.clone())
                    },
                    move |result| {
                        if let Err(error) = result {
                            let _ = status.try_set(error);
                        }
                    },
                );
            }
            Err(error) => status.set(error),
        }
    };
    let navigate =
        Callback::new(move |(capture, spans): (Capture, Vec<(usize, usize)>)| {
            match guard.run(capture) {
                Ok(()) => on_select.run(spans),
                Err(error) => status.set(error),
            }
        });
    view! {
        <details class="qc-vector-class">
            <summary>"Declared vector class"</summary>
            <p class="helper">"Uses the annotation assumptions, pinned references and literal copy requirements above. A matching element sequence does not certify its biological role."</p>
            <label>"Declared profile"
                <select disabled=move ||workspace_busy.get() prop:value=move ||form.with(|form|form.profile.clone())
                    on:change=move |event|{
                        let id=event_target_value(&event);
                        if id.is_empty(){form.set(Form::default());return;}
                        let selected=profiles.with_value(|profiles|profiles.as_ref().map_err(Clone::clone).and_then(|profiles|profiles.iter().find(|profile|profile.reference.id==id).ok_or_else(||"Unknown installed profile.".into()).and_then(Form::for_profile)));
                        match selected{Ok(selected)=>form.set(selected),Err(error)=>status.set(error)}
                    }>
                    <option value="">"Choose a sourced profile"</option>
                    {profiles.with_value(|profiles|profiles.as_ref().ok().map(|profiles|profiles.iter().map(|profile|view!{
                        <option value=profile.reference.id.clone()>{profile.definition.name.clone()}</option>
                    }).collect_view()))}
                </select>
            </label>
            {move ||profiles.with_value(|profiles|profiles.as_ref().ok().and_then(|profiles|profiles.iter().find(|profile|profile.reference.id==form.with(|form|form.profile.clone()))).map(|profile|{
                let role_buttons=profile.definition.roles.iter().map(|role|{
                    let id=role.id.clone();
                    view!{<button class="text-tool" title=role.description.clone() disabled=move ||workspace_busy.get()||form.with(|form|form.bindings.len()>=64)
                        on:click=move |_|{let result=form.try_update(|form|form.add_binding(&id));if let Some(Err(error))=result{status.set(error);}}>{format!("Add {} binding",role.id)}</button>}
                }).collect_view();
                view!{<div><p class="helper">{profile.definition.scope.clone()}</p>
                    <p class="helper">{format!("Boundary roles: {} → {}",profile.definition.start_role,profile.definition.end_role)}</p>
                    <details><summary>"Source and scope"</summary>
                        {profile.definition.sources.iter().map(|source|view!{<p class="helper"><a href=source.url.clone() target="_blank" rel="noopener noreferrer">{source.id.clone()}</a>{format!(" · {} · {}",source.locator,source.rights)}</p>}).collect_view()}
                        {profile.definition.limitations.iter().map(|text|view!{<p class="helper">{text.clone()}</p>}).collect_view()}
                    </details>
                    <details><summary>"Role bindings"</summary>
                        {move ||form.get().bindings.into_iter().map(|binding|binding_row(form,draft,rows,binding,workspace_busy)).collect_view()}
                        {role_buttons}
                    </details>
                    {binding_choice(form,"First boundary binding",|form|&form.first,|form,value|form.first=value,workspace_busy)}
                    {binding_choice(form,"Last boundary binding",|form|&form.last,|form,value|form.last=value,workspace_busy)}
                    <label>"Declared graph axis"<select disabled=move ||workspace_busy.get() prop:value=move ||if form.with(|form|form.reverse){"reverse"}else{"forward"}
                        on:change=move |event|form.update(|form|form.reverse=event_target_value(&event)=="reverse")>
                        <option value="forward">"Forward"</option><option value="reverse">"Reverse"</option>
                    </select></label>
                    <p class="helper">"Sizes include the selected boundary elements, with interior length reported separately. This does not substitute total plasmid length or predict packaging."</p>
                    <details><summary>"Additional explicit spacing constraints"</summary>
                        {move ||form.get().spacing.into_iter().map(|row|spacing_row(form,row,workspace_busy)).collect_view()}
                        <button class="text-tool" disabled=move ||workspace_busy.get()||form.with(|form|form.spacing.len()>=32)
                            on:click=move |_|{let result=form.try_update(Form::add_spacing);if let Some(Err(error))=result{status.set(error);}}>"Add spacing constraint"</button>
                    </details>
                </div>}
            }))}
            <p class="helper">{move ||draft.with(|draft|format!("Shared premises: {} assumed annotations · {} pinned references · {} literal-copy requirements",draft.assumed.len(),draft.options.references.len(),draft.options.required_parts.len()))}</p>
            <button class="text-tool" disabled=move ||busy.get()||workspace_busy.get()||source.with(|source|source.alphabet!=Alphabet::Dna)||form.with(|form|form.profile.is_empty()) on:click=run>"Check declared graph"</button>
            <p class="helper" role="status">{move ||status.get()}</p>
            {move ||report.get().map(|(capture,value)|{
                let total=value.findings.len();
                let first=page.get().saturating_mul(25);
                let findings=value.findings.iter().skip(first).take(25).map(|finding|{
                    let spans=finding.coordinates.iter().map(|span|(span.start,span.end)).collect::<Vec<_>>();
                    let selectable=!spans.is_empty();let saved=capture.clone();
                    let sources=finding.source_ids.iter().filter_map(|id|value.profile.definition.sources.iter().find(|source|source.id==*id)).map(|source|view!{<a href=source.url.clone() target="_blank" rel="noopener noreferrer">{source.id.clone()}</a>}).collect_view();
                    view!{<article class="qc-finding"><strong>{format!("{:?} · {:?}",finding.rule,finding.verdict)}</strong><p>{finding.message.clone()}</p>
                        <p class="helper">{format!("Bindings {:?} · annotation indices (1-based) {:?}",finding.binding_ids,finding.feature_indices.iter().map(|index|index+1).collect::<Vec<_>>())}</p>{sources}
                        <button class="text-tool" disabled=move ||!selectable||workspace_busy.get()
                            on:click=move |_|navigate.run((saved.clone(),spans.clone()))>"Select resolved coordinates"</button>
                    </article>}
                }).collect_view();
                let boundaries=value.boundaries.as_ref().map(|boundary|{
                    let spans=boundary.source_spans.iter().map(|span|(span.start,span.end)).collect::<Vec<_>>();let saved=capture.clone();
                    view!{<p>{format!("Declared boundary span {} bases · interior {} bases · source {} bases",boundary.outer_bases,boundary.interior_bases,value.source_length)}</p>
                        <button class="text-tool" disabled=move ||workspace_busy.get() on:click=move |_|navigate.run((saved.clone(),spans.clone()))>"Select boundary span"</button>}
                });
                view!{<div><p class="helper">{value.statement.clone()}</p>{boundaries}
                    <p class="helper">{format!("{} typed graph findings · {} exact reference premises",total,value.bindings.len())}</p>{findings}
                    <div class="pagination"><button class="text-tool" disabled=move ||page.get()==0 on:click=move |_|page.update(|page|*page=page.saturating_sub(1))>"Previous findings"</button>
                        <button class="text-tool" disabled=move ||page.get().saturating_add(1).saturating_mul(25)>=total on:click=move |_|page.update(|page|*page+=1)>"Next findings"</button></div>
                    <details><summary>"Literal-copy checks and descriptive size context"</summary>
                        {value.literal_parts.iter().map(|part|view!{<p>{format!("{:?} · {}: {}",part.verdict,part.id,part.reason)}</p>}).collect_view()}
                        {value.literal_searches.iter().map(|search|view!{<p class="helper">{format!("Reference {} · search complete {} · {} retained intervals · work {}",search.reference_id,search.complete,search.matches.len(),search.work)}</p>}).collect_view()}
                        {value.profile.definition.numeric_context.iter().map(|context|view!{<p class="helper">{format!("{}: {}–{} {}. {} {}",context.metric,context.minimum,context.maximum,context.unit,context.population,context.interpretation)}</p>}).collect_view()}
                    </details>
                    <button class="text-tool" disabled=move ||workspace_busy.get() on:click=export.clone()>"Export full vector-check evidence (JSON)"</button>
                </div>}
            })}
        </details>
    }
}
