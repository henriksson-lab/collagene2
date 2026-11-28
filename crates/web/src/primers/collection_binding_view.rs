//! Worker-backed oligo collection binding; all source actions retain exact context.
use super::collection_binding_model::{Capture, Draft, Field, MAX_TEXT_BYTES, spans, status_label};
use crate::{platform::worker::CachedWorker, sources::identity::SourceInfo};
use collagene_core::{
    cloning::pcr::collection_binding::{self, Report, TmValue},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::Molecule;
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
pub fn OligoCollectionPanel(
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    on_select: Callback<Vec<(usize, usize)>>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
    let result = RwSignal::new(None::<(Capture, Arc<Report>)>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
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
        let _ = (source.get(), epoch.get(), draft.get());
        stop.run(());
        result.set(None);
        status.set(String::new());
    });
    let guard = {
        let alive = alive.clone();
        Callback::new(move |capture: Capture| -> Result<(), String> {
            if !alive.load(Ordering::Relaxed) {
                return Err("Collection panel closed.".into());
            }
            if !capture.current(
                &source.get_untracked(),
                epoch.get_untracked(),
                generation.get_untracked(),
                &draft.get_untracked(),
                &get_molecule.run(()).revision(),
            ) {
                return Err(
                    "Source, collection, options or calculation changed; search again.".into(),
                );
            }
            if !result.with_untracked(|r| {
                r.as_ref().is_some_and(|(current, report)| {
                    current == &capture && report.context_revision == capture.context_revision
                })
            }) {
                return Err("Binding results are no longer current.".into());
            }
            Ok(())
        })
    };
    let run = {
        let alive = alive.clone();
        move |_| {
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let (collection, options, budget, conditions) = draft.get_untracked().prepare()?;
                let request = collection_binding::Request {
                    source: get_molecule.run(()),
                    collection,
                    options,
                    budget,
                    conditions,
                };
                let revision = collection_binding::revision(&request).map_err(|e| e.to_string())?;
                let snapshot = source.get_untracked();
                if snapshot.binding.revision != request.source.revision() {
                    return Err("Selected source event is not the current molecule.".into());
                }
                stop.run(());
                result.set(None);
                let capture = Capture {
                    source: snapshot,
                    epoch: epoch.get_untracked(),
                    generation: generation.get_untracked(),
                    draft: draft.get_untracked(),
                    context_revision: revision.clone(),
                    source_revision: request.source.revision(),
                };
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                busy.set(true);
                Ok((request, capture, w, revision))
            })();
            let (request, capture, w, revision) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            status.set("Searching the collection in a background worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: capture.generation,
                        revision: revision.clone(),
                        job: AnalysisJob::SearchOligoCollection { request },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.try_get_untracked() != Some(capture.generation)
                {
                    return;
                }
                if !capture.current(
                    &source.get_untracked(),
                    epoch.get_untracked(),
                    generation.get_untracked(),
                    &draft.get_untracked(),
                    &get_molecule.run(()).revision(),
                ) {
                    w.close();
                    worker.set_value(None);
                    busy.set(false);
                    return;
                }
                w.close();
                worker.set_value(None);
                busy.set(false);
                let checked = (|| -> Result<Report, String> {
                    let response = response?;
                    if response.id != capture.generation || response.revision != revision {
                        return Err(
                            "Worker response does not match this collection request.".into()
                        );
                    }
                    match response.result {
                        AnalysisResult::OligoCollectionBinding(report) => {
                            report
                                .validate_binding(&revision)
                                .map_err(|e| e.to_string())?;
                            if report.source_revision != capture.source_revision {
                                return Err("Worker source revision differs.".into());
                            }
                            Ok(report)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected collection-binding response.".into()),
                    }
                })();
                match checked {
                    Ok(report) => {
                        status.set(if report.search.complete {
                            "Collection search complete within the selected ungapped model/policy."
                        } else {
                            "Search incomplete. Inspect each oligo's completion state; unsearched results do not establish absence."
                        }.into());
                        result.set(Some((capture, Arc::new(report))));
                    }
                    Err(e) => {
                        let _ = status.try_set(e);
                    }
                }
            });
        }
    };
    let select = Callback::new(move |index: usize| {
        let checked = (|| -> Result<(), String> {
            let (capture, report) = result
                .try_get_untracked()
                .flatten()
                .ok_or("Search the collection first.")?;
            guard.run(capture.clone())?;
            report
                .validate_binding(&capture.context_revision)
                .map_err(|e| e.to_string())?;
            let hit = report
                .search
                .hits
                .get(index)
                .ok_or("Binding is unavailable.")?;
            let selection = spans(hit);
            if selection.is_empty() {
                return Err("Binding has no source columns.".into());
            }
            guard.run(capture)?;
            on_select.run(selection);
            Ok(())
        })();
        if let Err(e) = checked {
            let _ = status.try_set(e);
        }
    });
    let export_lifetime = RwSignal::new(());
    let export = Callback::new(move |()| {
        let prepared = (|| -> Result<_, String> {
            let (capture, report) = result
                .try_get_untracked()
                .flatten()
                .ok_or("Search the collection first.")?;
            guard.run(capture.clone())?;
            report
                .validate_binding(&capture.context_revision)
                .map_err(|e| e.to_string())?;
            let bytes = super::collection_binding_model::export_json(&capture, report.as_ref())?;
            let artifact = crate::io::download::ExportArtifact::prepare(
                "oligo-collection-binding.json",
                "application/json",
                &bytes,
            )?;
            Ok((artifact, capture))
        })();
        match prepared {
            Ok((artifact, capture)) => {
                let completion_capture = capture.clone();
                artifact.complete_guarded(
                    move || {
                        if export_lifetime.try_get_untracked().is_none() {
                            return Err("Collection panel closed.".into());
                        }
                        guard.run(capture.clone())
                    },
                    move |outcome| {
                        if export_lifetime.try_get_untracked().is_none()
                            || guard.run(completion_capture.clone()).is_err()
                        {
                            return;
                        }
                        if let Err(e) = outcome {
                            status.set(e);
                        }
                    },
                );
            }
            Err(e) => {
                let _ = status.try_set(e);
            }
        }
    });
    view! {
        <section class="primer-collection-panel">
            <h3>"Oligo collection binding"</h3>
            <p class="helper">
                "Search both strands for ungapped 3′-anchored footprints. Alignment scores rank sequence matches; they do not predict PCR extension, specificity or experimental success."
            </p>
            <label>
                "Oligos (FASTA or tab-separated rows)"
                <textarea
                    rows="8"
                    maxlength=MAX_TEXT_BYTES
                    prop:value=move || draft.get().text
                    on:input=move |e| draft.update(|d| d.text=event_target_value(&e))
                    placeholder=">forward tail=2\nGGACGTACGTACGT\n>reverse\nACGTACGTACGT"
                />
            </label>
            <p class="helper">
                "FASTA headers may end in tail=N. Rows use name, ordered DNA, optional declared 5′ tail length, separated by tabs. Declared tails and additional unpaired 5′ bases are reported separately. Maximum 1,024 oligos, 200 bases each, 200,000 total bases."
            </p>
            <div class="form-grid">
                {Field::MAIN.into_iter()
                    .map(|field| view! { <NumericInput draft field/> })
                    .collect_view()}
            </div>
            <p class="helper">
                "Zero required 3′ pairing keeps terminal mismatches visible. IUPAC results report possible/certain pairing and conservative mismatch ranges; the mismatch limit applies to each concrete realization."
            </p>
            <details>
                <summary>"Alignment scoring and resource limits"</summary>
                <div class="form-grid">
                    {Field::SCORES.into_iter()
                        .chain(Field::BUDGET)
                        .map(|field| view! { <NumericInput draft field/> })
                        .collect_view()}
                </div>
                <p class="helper">
                    "Every best-score suffix tie is retained at each endpoint, or that endpoint is omitted with an explicit limit. Ambiguous envelopes can include alternatives that never win in any concrete realization. Resource limits are shared across this collection."
                </p>
            </details>
            <details>
                <summary>"Thermodynamic conditions"</summary>
                <div class="form-grid">
                    {Field::THERMO.into_iter()
                        .map(|field| view! { <NumericInput draft field/> })
                        .collect_view()}
                </div>
                <p class="helper">
                    "Sodium-only perfect-duplex or admitted single internal mismatch model (pH7). Whole-oligo Tm assumes a hypothetical perfectly complementary partner, including tails. Annealed-stretch Tm is unavailable for ambiguous, multiple-mismatch, near-end or other unsupported duplexes."
                </p>
            </details>
            <div class="digest-actions">
                <button
                    class="button secondary"
                    disabled=move || busy.get()
                    on:click=run
                >
                    "Search collection"
                </button>
                <button
                    class="text-tool"
                    disabled=move || !busy.get()
                    on:click=move |_| {
                        stop.run(());
                        result.set(None);
                        status.set("Collection search cancelled.".into());
                    }
                >
                    "Cancel"
                </button>
            </div>
            <p class="helper" role="status">{move || status.get()}</p>
            {move || result.get().map(|(capture, report)| {
                let names=capture.draft.prepare()
                    .map(|(c, _, _, _)| c.oligos.into_iter()
                        .map(|o| o.name)
                        .collect::<Vec<_>>())
                    .unwrap_or_default();
                view! { <CollectionResults report names select export/> }
            })}
        </section>
    }
}
#[component]
fn NumericInput(draft: RwSignal<Draft>, field: Field) -> impl IntoView {
    view! {
        <label>
            {field.label()}
            <input
                type="number"
                step="any"
                prop:value=move || draft.get().value(field)
                on:input=move |e| draft.update(|d| d.set(field, event_target_value(&e)))
            />
        </label>
    }
}
fn tm(value: &TmValue) -> String {
    match value {
        TmValue::Available(report) => format!("{:.1} °C ({})", report.tm_c, report.model),
        TmValue::MismatchAvailable(report) => format!(
            "{:.1} °C ({}; profile {})",
            report.estimate.tm_c, report.estimate.model, report.parameters.profile_revision
        ),
        TmValue::Unavailable(reason) => format!("Unavailable: {reason}"),
    }
}
#[component]
fn CollectionResults(
    report: Arc<Report>,
    names: Vec<String>,
    select: Callback<usize>,
    export: Callback<()>,
) -> impl IntoView {
    let progress=report.search.oligos.iter().enumerate()
        .map(|(i, p)| view! {
            <tr>
                <td>{format!("{} · {}", p.oligo_id, names.get(i).cloned().unwrap_or_default())}</td>
                <td>{status_label(p.status, p.alignments_retained)}</td>
                <td>{format!("{}/{}", p.endpoints_examined, p.endpoints_total)}</td>
            </tr>
        })
        .collect_view();
    let total = report.search.hits.len();
    let pages = total.div_ceil(50).max(1);
    let page = RwSignal::new(0usize);
    let usage = format!(
        "Work charged: {}; hits retained: {}; output allocation charge: {} bytes.",
        report
            .search
            .budget_after
            .work
            .saturating_sub(report.search.budget_before.work),
        total,
        report
            .search
            .budget_after
            .output_charge
            .saturating_sub(report.search.budget_before.output_charge)
    );
    let warnings = report
        .warnings
        .iter()
        .map(|w| view! { <p class="digest-warning">{w.clone()}</p> })
        .collect_view();
    let report = StoredValue::new(report);
    let names = StoredValue::new(names);
    let hits = move || {
        report.with_value(|report| report.search.hits.iter().enumerate()
            .skip(page.get()*50)
            .take(50)
            .map(|(index, h)| {
                let thermo=report.thermodynamics.iter().find(|t| t.hit_index==index);
                let annealing=thermo.map(|t| tm(&t.annealing))
                    .unwrap_or_else(|| "Unavailable: no model result".into());
                let whole=thermo.map(|t| tm(&t.whole_oligo))
                    .unwrap_or_else(|| "Unavailable: no model result".into());
                let name=names.with_value(|names| names.get(h.oligo_index).cloned().unwrap_or_default());
                let description=format!(
                    "{} · {} · {} · {} bp{} · mismatches {}–{} · score {}–{}",
                    h.oligo_id,
                    name,
                    if h.orientation==collagene_oligo_binding::Orientation::Forward {
                        "forward"
                    } else {
                        "reverse"
                    },
                    h.length,
                    if h.crosses_origin { " · crosses origin" } else { "" },
                    h.minimum_mismatches,
                    h.maximum_mismatches,
                    h.minimum_score,
                    h.maximum_score
                );
                let coordinates=spans(h).into_iter()
                    .map(|(a, b)| format!("{}–{}", a+1, b))
                    .collect::<Vec<_>>()
                    .join(", ");
                let columns=h.columns.iter()
                    .map(|c| format!(
                        "oligo {} {} ↔ template {} · source {} {} · {:?}",
                        c.oligo_index+1,
                        c.oligo_base,
                        c.template_base,
                        c.source_position+1,
                        c.source_base,
                        c.pairing
                    ))
                    .collect::<Vec<_>>()
                    .join("\n");
                view! {
                    <details>
                        <summary>{description}</summary>
                        <p>{format!(
                            "Source {coordinates} (1-based inclusive); ordered oligo annealing {}–{}. Declared tail {} bases; additional unpaired prefix {} bases. 3′ pairing: {} certain / {} possible contiguous bases. Admissible in every represented realization: {}. Ranking: {:?}.",
                            h.annealing_start+1,
                            h.annealing_start+h.length,
                            h.declared_tail_length,
                            h.extra_unpaired_length,
                            h.three_prime_certain,
                            h.three_prime_possible,
                            h.certainly_admissible,
                            h.score_scope
                        )}</p>
                        <p>{format!("Actual annealed-stretch Tm: {annealing}")}</p>
                        <p>{format!("Hypothetical whole-oligo perfect-duplex Tm: {whole}")}</p>
                        <button class="text-tool" on:click=move |_| select.run(index)>
                            "Select source footprint"
                        </button>
                        <details>
                            <summary>"Aligned columns (oligo 5′ → 3′)"</summary>
                            <pre>{columns}</pre>
                        </details>
                    </details>
                }
            })
            .collect_view())
    };
    view! {
        <div class="primer-collection-results">
            <table>
                <thead>
                    <tr>
                        <th>"Oligo"</th>
                        <th>"Completion"</th>
                        <th>"Endpoints examined"</th>
                    </tr>
                </thead>
                <tbody>{progress}</tbody>
            </table>
            <p class="helper">{usage}</p>
            {warnings}
            <div class="digest-actions">
                <button
                    class="text-tool"
                    disabled=move || page.get()==0
                    on:click=move |_| page.update(|p| *p=p.saturating_sub(1))
                >
                    "Previous hits"
                </button>
                <span>{move || format!(
                    "Hit page {} of {} · up to 50 alignments per page",
                    page.get()+1,
                    pages
                )}</span>
                <button
                    class="text-tool"
                    disabled=move || page.get()+1>=pages
                    on:click=move |_| page.update(|p| *p=(*p+1).min(pages-1))
                >
                    "Next hits"
                </button>
            </div>
            <div>{hits}</div>
            <button class="button secondary" on:click=move |_| export.run(())>
                "Download complete binding evidence (JSON)"
            </button>
        </div>
    }
}
