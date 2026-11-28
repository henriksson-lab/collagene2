//! Opt-in documentary formulation/assay review inside Virtual digest.
use super::reaction_conditions_model::{Capture, Draft, celsius};
use crate::{
    history::run_origin_model::PinnedOrigin, platform::worker::CachedWorker,
    sources::identity::SourceInfo,
};
use collagene_core::{
    cloning::restriction::{
        catalog::{CatalogSelection, catalog},
        conditions::{self, CommonCondition, ConditionsReport, EvidenceStatus, ReactionConditions},
    },
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[component]
pub(super) fn ReactionConditionsPanel(
    catalog: StoredValue<Result<Arc<catalog::Catalog>, String>>,
    selection: Memo<Result<CatalogSelection, String>>,
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    origin: RwSignal<Option<PinnedOrigin>>,
    conditions: RwSignal<Option<ReactionConditions>>,
    error: RwSignal<Option<String>>,
    loaded: RwSignal<u64>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::from_saved(conditions.get_untracked()));
    let report = RwSignal::new(None::<Arc<ConditionsReport>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u64);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let stop = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
        });
    }
    Effect::new(move |_| {
        let _ = loaded.get();
        draft.set(Draft::from_saved(conditions.get_untracked()));
    });
    Effect::new(move |_| {
        if let Ok(s) = selection.get() {
            let old = draft.get_untracked();
            let mut next = old.clone();
            next.keep_selected(&s);
            if next != old {
                draft.set(next);
            }
        }
    });
    let options = Memo::new(move |_| -> Result<Option<ReactionConditions>, String> {
        let d = draft.get();
        if !d.enabled {
            return Ok(None);
        }
        d.options(catalog.get_value()?.as_ref(), &selection.get()?)
    });
    Effect::new(move |_| match options.get() {
        Ok(value) => {
            error.set(None);
            conditions.set(value);
        }
        Err(message) => {
            error.set(Some(message));
            conditions.set(None);
        }
    });
    Effect::new(move |_| {
        let _ = (
            source.get(),
            epoch.get(),
            origin.get(),
            selection.get(),
            draft.get(),
        );
        stop.run(());
        report.set(None);
        status.set(String::new());
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let c = catalog.get_value()?;
                let selected = selection.get_untracked()?;
                let options = options
                    .get_untracked()?
                    .ok_or("Enable reaction conditions first.")?;
                let revision =
                    conditions::revision(&c, &selected, &options).map_err(|e| e.to_string())?;
                let capture = Capture {
                    source: source.get_untracked().binding,
                    epoch: epoch.get_untracked(),
                    selection: selected,
                    options,
                };
                let source_snapshot = source.get_untracked();
                let origin_snapshot = origin.get_untracked();
                stop.run(());
                report.set(None);
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                busy.set(true);
                Ok((
                    capture,
                    source_snapshot,
                    origin_snapshot,
                    revision,
                    w,
                    generation.get_untracked(),
                ))
            })();
            let (capture, source_snapshot, origin_snapshot, revision, w, expected_generation) =
                match prepared {
                    Ok(p) => p,
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                };
            status.set("Reviewing exact formulation evidence…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: 1,
                        revision: revision.clone(),
                        job: AnalysisJob::AssessRestrictionConditions {
                            selection: capture.selection.clone(),
                            options: capture.options.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.try_get_untracked() != Some(expected_generation)
                    || source.get_untracked() != source_snapshot
                    || origin.get_untracked() != origin_snapshot
                    || !selection.get_untracked().is_ok_and(|selected| {
                        options.get_untracked().is_ok_and(|configured| {
                            configured.as_ref().is_some_and(|configured| {
                                capture.current(
                                    &source.get_untracked().binding,
                                    epoch.get_untracked(),
                                    &selected,
                                    configured,
                                )
                            })
                        })
                    })
                {
                    return;
                }
                w.close();
                worker.set_value(None);
                busy.set(false);
                let result = (|| -> Result<ConditionsReport, String> {
                    let response = response?;
                    if response.id != 1 || response.revision != revision {
                        return Err("Condition response does not match this request.".into());
                    }
                    let AnalysisResult::RestrictionConditions(result) = response.result else {
                        return Err("Unexpected reaction-condition response.".into());
                    };
                    result
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if result.selection != capture.selection || result.options != capture.options {
                        return Err("Condition response choices differ.".into());
                    }
                    Ok(result)
                })();
                match result {
                    Ok(value) => {
                        status.set("Current documentary assessment; cleavage products remain a separate calculation.".into());
                        report.set(Some(Arc::new(value)));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let choose = Callback::new(move |value: CommonCondition| {
        draft.update(|d| {
            d.buffer = value.buffer;
            d.temperature = celsius(value.temperature_milli_c);
            d.activity_context = value.activity_context;
            d.temperature_context = value.temperature_context;
        })
    });
    view! {<details class="digest-reaction-conditions"><summary>"Reaction conditions"</summary>
        <label class="checkbox-label"><input type="checkbox" prop:checked=move||draft.get().enabled on:change=move|e|draft.update(|d|d.enabled=event_target_checked(&e))/><span>"Retain an explicit reaction setup"</span></label>
        <Show when=move||draft.get().enabled>
            <p class="helper">"Choose the exact vendor formulation for each enzyme. Unselected formulations remain unknown; ordinary enzymes are never replaced by HF variants."</p>
            {move||selection.get().ok().map(|selection|selection.enzymes.into_iter().map(|enzyme| {
                let enzyme_id=enzyme.id.clone();let change=enzyme.id.clone();
                let forms=catalog.with_value(|c|c.as_ref().ok().map(|c|c.edition().formulations.iter().filter(|f|f.enzyme==enzyme.id).cloned().collect::<Vec<_>>()).unwrap_or_default());
                view!{<label>{format!("{} formulation",enzyme.enzyme.name)}<select prop:value=move||draft.with(|d|d.formulations.get(&enzyme_id).cloned().flatten().unwrap_or_default()) on:change=move|e|{let value=event_target_value(&e);draft.update(|d|{d.formulations.insert(change.clone(),if value.is_empty(){None}else{Some(value)});});}>
                    <option value="">"Not selected · unknown"</option>{forms.into_iter().map(|f|view!{<option value=f.id>{format!("{} · {} {}",f.name,f.vendor,f.product_code)}</option>}).collect_view()}
                </select></label>}
            }).collect_view())}
            <label>"Buffer"<input maxlength="128" list="digest-condition-buffer-list" prop:value=move||draft.get().buffer on:input=move|e|draft.update(|d|d.buffer=event_target_value(&e)) placeholder="Choose a recorded buffer or enter its name"/></label>
            <datalist id="digest-condition-buffer-list">{move||catalog.with_value(|c|c.as_ref().ok().and_then(|c|selection.get().ok().map(|s|draft.get().buffers(c,&s))).unwrap_or_default().into_iter().map(|value|view!{<option value=value/>}).collect_view())}</datalist>
            <label>"Temperature (°C)"<input type="number" min="0" max="100" step="0.001" prop:value=move||draft.get().temperature on:input=move|e|draft.update(|d|d.temperature=event_target_value(&e))/></label>
            <label>"Incubation (minutes)"<input type="number" min="1" max="100800" step="1" prop:value=move||draft.get().minutes on:input=move|e|draft.update(|d|d.minutes=event_target_value(&e))/></label>
            <label>"Minimum recorded activity (%)"<input type="number" min="1" max="100" step="1" prop:value=move||draft.get().minimum on:input=move|e|draft.update(|d|d.minimum=event_target_value(&e))/></label>
            <details><summary>"Assay contexts"</summary><p class="helper">"Leave a context unspecified only when the catalog evidence resolves uniquely. Multiple contexts require an explicit choice."</p>
                <label>"Buffer-activity context"<select prop:value=move||draft.get().activity_context on:change=move|e|draft.update(|d|d.activity_context=event_target_value(&e))><option value="">"Unspecified"</option>{move||context_options(catalog,selection,draft,false)}</select></label>
                <label>"Incubation-temperature context"<select prop:value=move||draft.get().temperature_context on:change=move|e|draft.update(|d|d.temperature_context=event_target_value(&e))><option value="">"Unspecified"</option>{move||context_options(catalog,selection,draft,true)}</select></label>
            </details>
            <p class="digest-warning">{move||error.get().unwrap_or_default()}</p>
            <div class="digest-actions"><button class="button secondary" disabled=move||busy.get()||!matches!(options.get(),Ok(Some(_))) on:click=run.clone()>"Assess reaction conditions"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|{stop.run(());status.set("Reaction assessment cancelled.".into());}>"Cancel assessment"</button></div>
            <p class="helper" role="status">{move||status.get()}</p>
            {move||report.get().map(|report|view!{<ConditionsEvidence report on_choose=Some(choose)/>})}
        </Show>
    </details>}
}
fn context_options(
    c: StoredValue<Result<Arc<catalog::Catalog>, String>>,
    s: Memo<Result<CatalogSelection, String>>,
    d: RwSignal<Draft>,
    temperature: bool,
) -> impl IntoView {
    c.with_value(|c| {
        let draft = d.get();
        let current = if temperature {
            draft.temperature_context.clone()
        } else {
            draft.activity_context.clone()
        };
        let mut values = c
            .as_ref()
            .ok()
            .and_then(|c| s.get().ok().map(|s| draft.contexts(c, &s, temperature)))
            .unwrap_or_default();
        let unavailable = !current.is_empty() && !values.contains(&current);
        if unavailable {
            values.push(current.clone());
        }
        values
            .into_iter()
            .map(|value| {
                let label = if unavailable && value == current {
                    format!("{value} · not recorded for these choices")
                } else {
                    value.clone()
                };
                view! {<option value=value>{label}</option>}
            })
            .collect_view()
    })
}
fn status_label(value: EvidenceStatus) -> &'static str {
    match value {
        EvidenceStatus::Known => "Recorded evidence meets this request",
        EvidenceStatus::BelowThreshold => "Below the requested activity threshold",
        EvidenceStatus::Unknown => "Unknown",
        EvidenceStatus::Conflict => "Conflicting records",
        EvidenceStatus::ContextRequired => "Choose an assay context",
        EvidenceStatus::TemperatureMismatch => "Recorded temperature differs",
    }
}
#[component]
pub(super) fn ConditionsEvidence(
    report: Arc<ConditionsReport>,
    on_choose: Option<Callback<CommonCondition>>,
) -> impl IntoView {
    let suggestion = RwSignal::new(0usize);
    let saved = StoredValue::new(report.clone());
    view! {<section class="digest-condition-results"><h4>"Documentary condition evidence"</h4><p><strong>{status_label(report.documentary_status)}</strong></p>
        <p>{format!("{} · {} °C · {} minutes declared · minimum recorded activity {}%",report.options.buffer,celsius(report.options.temperature_milli_c),report.options.incubation_minutes,report.options.minimum_activity_percent)}</p>
        <p class="helper">"Recorded unit-assay activity is not percent cleavage or reaction success. Star activity and incubation-duration suitability remain unknown; methylation is assessed separately."</p>
        {report.formulations.iter().map(|row|{
            let name=report.selection.enzymes.iter().find(|e|e.id==row.enzyme_id).map(|e|e.enzyme.name.clone()).unwrap_or_else(||row.enzyme_id.clone());
            let formulation=row.formulation.as_ref().map(|f|format!("{} · {} {}",f.name,f.vendor,f.product_code)).unwrap_or_else(||"No formulation selected · unknown".into());
            view!{<article><strong>{format!("{name} · {formulation}")}</strong><p>{format!("Buffer activity: {}. Incubation temperature: {}.",status_label(row.activity.status),status_label(row.temperature.status))}</p>
                <details><summary>"Formulation qualifiers and citations"</summary>
                    {row.formulation.as_ref().map(|form|{
                        let facts=form.evidence.iter().map(|fact|view!{<section><p>{format!("{} · context: {} · {}",condition_label(&fact.condition),fact.context,value_label(&fact.value))}</p>{fact.qualifiers.iter().map(|q|view!{<p class="helper">{q.clone()}</p>}).collect_view()}{citations(&report.sources,&fact.citations)}</section>}).collect_view();
                        view!{<div><p>{form.formulation_scope.clone()}</p>{facts}{citations(&report.sources,&form.citations)}</div>}
                    })}
                    <p class="helper">"No selected record can establish star-activity risk or that the declared incubation duration is sufficient."</p>
                </details></article>}
        }).collect_view()}
        <details><summary>"Common recorded buffer and temperature combinations"</summary>
            <p class="helper">{format!("{} combinations · catalogue search {}. Suggestions apply only to the exact chosen formulations and unit-assay contexts; they do not establish a successful joint reaction.",report.common_conditions.len(),if report.suggestions_complete{"complete"}else{"limited"})}</p>
            {(!report.common_conditions.is_empty()).then(||view!{<select aria-label="Recorded common reaction conditions" on:change=move|e|suggestion.set(event_target_value(&e).parse().unwrap_or(usize::MAX))>{report.common_conditions.iter().enumerate().map(|(i,c)|view!{<option value=i.to_string()>{format!("{} · {} °C · {} / {}",c.buffer,celsius(c.temperature_milli_c),c.activity_context,c.temperature_context)}</option>}).collect_view()}</select>})}
            {on_choose.map(|choose|view!{<button class="text-tool" disabled=report.common_conditions.is_empty() on:click=move|_|{if let Some(value)=saved.with_value(|r|r.common_conditions.get(suggestion.get_untracked()).cloned()){choose.run(value);}}>"Use selected recorded conditions"</button>})}
        </details>
    </section>}
}
fn condition_label(condition: &catalog::Condition) -> String {
    match condition {
        catalog::Condition::IncubationTemperature => "Incubation temperature".into(),
        catalog::Condition::BufferActivity {
            buffer,
            temperature_milli_c,
        } => format!(
            "Activity in {buffer} at {} °C",
            celsius(*temperature_milli_c)
        ),
        catalog::Condition::Methylation { system } => format!("{system:?} methylation profile"),
    }
}
fn value_label(value: &catalog::FactValue) -> String {
    match value {
        catalog::FactValue::TemperatureMilliC(value) => format!("{} °C", celsius(*value)),
        catalog::FactValue::Percent(value) => format!("{value}% unit-assay activity"),
        catalog::FactValue::LessThanPercent(value) => {
            format!("less than {value}% unit-assay activity")
        }
        catalog::FactValue::Unknown(reason) => format!("Unknown: {reason}"),
        catalog::FactValue::Methylation(value) => format!("{value:?}"),
    }
}
fn citations(
    sources: &[catalog::Source],
    citations: &[catalog::Citation],
) -> impl IntoView + use<> {
    citations.iter().filter_map(|citation|sources.iter().find(|source|source.id==citation.source).filter(|source|source.url.starts_with("https://")).map(|source|view!{<p><a href=source.url.clone() target="_blank" rel="noopener noreferrer">{format!("{} · {}",source.publisher,citation.locator)}</a><span>{format!(" · retrieved {}",source.retrieved_on)}</span></p>})).collect_view()
}
