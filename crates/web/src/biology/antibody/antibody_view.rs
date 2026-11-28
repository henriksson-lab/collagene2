//! Conditional antibody reference analysis, kept off the UI thread.
use crate::biology::antibody::antibody_model;
use crate::biology::antibody::antibody_model::Capture;
use crate::biology::antibody::antibody_overlay::Overlay;
use crate::biology::antibody::antibody_overlay::origin;
use crate::biology::antibody::antibody_overlay::region;
use crate::platform::worker::CachedWorker;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    antibody_numbering::domain_tracks as n,
    material::lineage::EventBinding,
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
#[component]
pub fn AntibodyPanel(
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    prepare: Callback<EventBinding, Result<Molecule, String>>,
    on_select: Callback<(Capture, usize, usize), Result<(), String>>,
    on_overlay: Callback<Option<Arc<Overlay>>, Result<(), String>>,
) -> impl IntoView {
    let references = RwSignal::new(
        collagene_core::antibody_numbering::v_mapping::ReferenceSelection::human_all(),
    );
    let mouse = Memo::new(move |_| {
        references.with(|r| {
            r.species == collagene_core::antibody_numbering::references::Species::MouseC57Bl6
        })
    });
    let scheme = RwSignal::new(n::NumberingScheme::Imgt);
    let definition = RwSignal::new(n::RegionDefinition::Imgt);
    let choice = Memo::new(move |_| n::Selection {
        numbering_scheme: scheme.get(),
        region_definition: definition.get(),
    });
    let minimum = RwSignal::new("300".to_string());
    let substitutions = RwSignal::new("1".to_string());
    let loop_max = RwSignal::new("60".to_string());
    let retained = RwSignal::new("32".to_string());
    let settings = Memo::new(move |_| -> Result<_, String> {
        let reference_scope = references.get();
        antibody_model::validate_display_scope(&reference_scope, &choice.get())?;
        Ok(antibody_model::preset(
            reference_scope,
            minimum
                .get()
                .parse()
                .map_err(|_| "Enter an integer V score.")?,
            substitutions
                .get()
                .parse()
                .map_err(|_| "Enter a whole FR4 substitution count.")?,
            loop_max
                .get()
                .parse()
                .map_err(|_| "Enter a whole maximum CDR3 length.")?,
            retained
                .get()
                .parse()
                .map_err(|_| "Enter a whole candidate limit.")?,
        ))
    });
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<n::Report>)>);
    let candidate = RwSignal::new(0usize);
    let page = RwSignal::new(0usize);
    Effect::new(move |_| {
        candidate.get();
        page.set(0);
    });
    let show_overlay = RwSignal::new(true);
    Effect::new(move |_| {
        let value = report.get();
        let index = candidate.get();
        let visible = show_overlay.get();
        let _ = choice.get();
        let result = (|| {
            if !visible {
                return on_overlay.run(None);
            }
            let Some((capture, r)) = value else {
                return on_overlay.run(None);
            };
            let (o, l) = settings.get()?;
            capture.guard(
                epoch.get(),
                &source.get().binding,
                &o,
                &l,
                &choice.get_untracked(),
            )?;
            if r.tracks.is_empty() {
                return on_overlay.run(None);
            }
            let overlay = Overlay::new(&capture, &r, index)?;
            on_overlay.run(Some(Arc::new(overlay)))
        })();
        if let Err(e) = result {
            let _ = on_overlay.run(None);
            status.set(e);
        }
    });

    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (
            source.get(),
            epoch.get(),
            settings.get(),
            choice.get(),
            workspace_busy.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        candidate.set(0);
        page.set(0);
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
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            let prepared = (|| {
                let (o, l) = settings.get_untracked()?;
                let s = source.get_untracked();
                let m = prepare.run(s.binding.clone())?;
                let c = Capture::prepare(
                    epoch.get_untracked(),
                    s.binding,
                    &m,
                    o,
                    l,
                    choice.get_untracked(),
                )?;
                Ok::<_, String>((m, c))
            })();
            let (m, c) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Some(old) = worker.get_value() {
                old.close();
            }
            worker.set_value(Some(w.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            report.set(None);
            busy.set(true);
            status.set("Analyzing selected V/J reference candidates…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: c.context_revision.clone(),
                        job: AnalysisJob::AnalyzeAntibody {
                            molecule: m,
                            options: c.options.clone(),
                            limits: c.limits.clone(),
                            selection: c.selection.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                let result = (|| {
                    let (o, l) = settings.get_untracked()?;
                    c.guard(
                        epoch.get_untracked(),
                        &source.get_untracked().binding,
                        &o,
                        &l,
                        &choice.get_untracked(),
                    )?;
                    match result?.result {
                        AnalysisResult::Antibody(r) => {
                            c.validate_report(&r)?;
                            Ok(r)
                        }
                        AnalysisResult::Error(e) => Err(e),
                        _ => Err("Unexpected antibody analysis response.".into()),
                    }
                })();
                match result {
                    Ok(r) => {
                        status.set(format!("{} retained candidates. Conditional reference evidence; review unresolved assignments.",r.tracks.len()));
                        report.set(Some((c, Arc::new(r))));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    view! {
        <section class="inspector-content" aria-label="Antibody reference analysis" style="overflow-wrap:anywhere">
            <details>
                <summary>
                    "Antibody numbering and CDRs"
                </summary>
                <p class="helper">
                    "Explicit V/J reference scope · up to 2,000 residues · conditional numbering and independently chosen regions. Unresolved residues remain unassigned. Published profile assignments remain conditional; this is not a calibrated classifier or independently certified numbering."
                </p>
                <details>
                    <summary>
                        "Reference scope and numbering"
                    </summary>
                    <super::reference_controls::ReferenceControls references status/>
                    <label>
                        "Numbering scheme"
                        <select
                            prop:value=move || format!("{:?}", scheme.get())
                            on:change=move |e| match antibody_model::numbering_scheme(&event_target_value(&e)) {
                        Ok(value) => scheme.set(value),
                        Err(error) => status.set(error),
                    }>
                            <option value="Imgt">
                                "IMGT"
                            </option>
                            <option value="Kabat1991" disabled=move || mouse.get()>
                                "Kabat (1991 profiles)"
                            </option>
                            <option value="Chothia1997" disabled=move || mouse.get()>
                                "Chothia (1997 profiles)"
                            </option>
                        </select>
                    </label>
                    <label>
                        "Region definition"
                        <select
                            prop:value=move || format!("{:?}", definition.get())
                            on:change=move |e| match antibody_model::region_definition(&event_target_value(&e)) {
                        Ok(value) => definition.set(value),
                        Err(error) => status.set(error),
                    }>
                            <option value="Imgt">
                                "IMGT"
                            </option>
                            <option value="Kabat1991" disabled=move || mouse.get()>
                                "Kabat (1991)"
                            </option>
                            <option value="ChothiaCanonical" disabled=move || mouse.get()>
                                "Chothia (canonical charts)"
                            </option>
                            <option value="North2011" disabled=move || mouse.get()>
                                "North (2011)"
                            </option>
                        </select>
                    </label>
                    <p class="helper">
                        "Number labels and CDR definitions are independent. Unsupported species/scheme combinations refuse analysis; unresolved windows are not guessed."
                    </p>
                    <p class="helper" role="alert">
                        {move || settings.with(|r|r.as_ref().err().cloned().unwrap_or_default())}
                    </p>
                </details>
                <details>
                    <summary>
                        "Analysis settings"
                    </summary>
                    <p class="helper">
                        "V weights: observed 5, frequency 8, unseen −6; affine gaps −12/−2. Minimum CDR3 length 2. Score thresholds are heuristic."
                    </p>
                    <label>
                        "Minimum V score"
                        <input
                            type="number"
                            prop:value=move ||minimum.get()
                            on:input=move |e|minimum.set(event_target_value(&e))/>
                    </label>
                    <label>
                        "Allowed FR4 substitutions (0–3)"
                        <input
                            type="number"
                            min="0"
                            max="3"
                            prop:value=move ||substitutions.get()
                            on:input=move |e|substitutions.set(event_target_value(&e))/>
                    </label>
                    <label>
                        "Maximum CDR3 residues"
                        <input
                            type="number"
                            min="2"
                            max="2000"
                            prop:value=move ||loop_max.get()
                            on:input=move |e|loop_max.set(event_target_value(&e))/>
                    </label>
                    <label>
                        "Retained candidates (0–64)"
                        <input
                            type="number"
                            min="0"
                            max="64"
                            prop:value=move ||retained.get()
                            on:input=move |e|retained.set(event_target_value(&e))/>
                    </label>
                    <p class="helper">
                        "Zero retained candidates requests counts only. At most 128 optimal V endpoints per locus are retained; omitted endpoints and suboptimal paths are not searched for products."
                    </p>
                </details>
                <button class="button secondary full" disabled=move ||busy.get()||workspace_busy.get() on:click=run>
                    "Analyze antibody candidates"
                </button>
                <Show when=move ||busy.get()>
                    <button
                        class="text-tool"
                        on:click=move |_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Analysis cancelled.".into());}>
                        "Cancel"
                    </button>
                </Show>
                <p class="helper" role="status">
                    {move || status.get()}
                </p>
                <Show when=move ||report.with(|r|r.is_some())>
                    {
                        move || report.with(|r|r.as_ref().map(|(capture, report)|view! {
                            <p class="helper">
                                {antibody_model::reference_label(&capture.options.references)}
                            </p>
                            <super::candidate_tracks::CandidateTracks report = report.clone()selected = candidate source_length = capture
                            .source_length/>
                        }))
                    }
                    <p class="helper">
                        {
                            move || report.with(|r|r.as_ref().map(|(_, r)|format!("{} qualifying pairs within retained V ends. V endpoints truncated: {} · pairs truncated: {}.",
                            r.evidence.qualified_pairs_in_retained_v, r.evidence.v_retention_truncated,
                            r.evidence.pair_retention_truncated)).unwrap_or_default())
                        }
                    </p>
                    <label>
                        <input
                            type="checkbox"
                            prop:checked=move ||show_overlay.get()
                            on:change=move |e|show_overlay.set(event_target_checked(&e))/>
                        <span>
                            "Show selected candidate in sequence"
                        </span>
                    </label>
                    <label>
                        "Candidate"
                        <select
                            prop:value=move ||candidate.get().to_string()
                            on:change=move |e|{candidate.set(event_target_value(&e).parse().unwrap_or(0));page.set(0);}>
                            {
                                move || report.with(|r|r.as_ref().map(|(_, r)|r.tracks.iter().enumerate()
                                .map(|(i, t)|view! {
                                    <option value=i.to_string()>
                                        {
                                            format!("{} · {:?} · supported {}–{}", i + 1, t.locus, t.first_supported_source_index + 1,
                                            t.source_end)
                                        }
                                    </option>
                                }).collect_view()))
                            }
                        </select>
                    </label>
                    {
                        move || report.with(|r|r.as_ref().map(|(_, r)|crate::biology::antibody::antibody_overlay::cdr_summaries(r,
                        candidate.get()).into_iter().map(|s|view! {
                            <p class="helper">
                                {s}
                            </p>
                        }).collect_view()))
                    }
                    <details>
                        <summary>
                            "Mapping diagnostics"
                        </summary>
                        {
                            move || report.with(|r|r.as_ref().map(|(_, r)|crate::biology::antibody::antibody_overlay::diagnostics(r,
                            candidate.get()).into_iter().map(|s|view! {
                                <p class="helper">
                                    {s}
                                </p>
                            }).collect_view()))
                        }
                    </details>
                    <p class="helper">
                        {
                            move || report.with(|r|r.as_ref().and_then(|(_, r)|r.tracks.get(candidate
                            .get())).map(|t|format!("{} unnumbered residues · {} unresolved regions · position 1: {}. Coordinates below are one-based source residues.", t.unassigned_residues, t.unassigned_regions,
                            t.invariant_position_one.map(|i|(i + 1).to_string()).unwrap_or_else(|| "unresolved N terminus"
                            .into()))).unwrap_or_default())
                        }
                    </p>
                    <table class="data-table">
                        <thead>
                            <tr>
                                <th>
                                    "Source"
                                </th>
                                <th>
                                    "AA"
                                </th>
                                <th>
                                    {
                                        move || crate::biology::antibody::antibody_overlay::scheme_label(&scheme
                                        .get())
                                    }
                                </th>
                                <th>
                                    "Region"
                                </th>
                            </tr>
                        </thead>
                        <tbody>
                            {
                                move || report.with(|r|r.as_ref().and_then(|(c, r)|r.tracks.get(candidate
                                .get()).map(|t|t.residues.iter().skip(page.get() * 60).take(60).map(|row| {
                                    let index = row.source_index;
                                    let capture = c.clone();
                                    let assignment = row.assignment.as_ref();
                                    view! {
                                        <tr>
                                            <td>
                                                <button
                                                    class="text-tool"
                                                    disabled=move ||workspace_busy.get()
                                                    on:click=move |_|{let result=(||{let (o,l)=settings.get_untracked()?;capture.guard(epoch.get_untracked(),&source.get_untracked().binding,&o,&l, &choice.get_untracked())?;on_select.run((capture.clone(),index,index+1))})();if let Err(e)=result{status.set(e);}}>
                                                    {index + 1}
                                                </button>
                                            </td>
                                            <td>
                                                {row.residue.to_string()}
                                            </td>
                                            <td title=assignment.map(|a|origin(&a.evidence)).unwrap_or("No supported assignment")>
                                                {
                                                    assignment.map(|a|a.position.to_string()).unwrap_or_else(|| "—"
                                                    .into())
                                                }
                                            </td>
                                            <td>
                                                {row.region.as_ref().map(region).unwrap_or("unresolved")}
                                            </td>
                                        </tr>
                                    }
                                }).collect_view())))
                            }
                        </tbody>
                    </table>
                    <button
                        class="text-tool"
                        disabled=move ||page.get()==0
                        on:click=move |_|page.update(|p|*p=p.saturating_sub(1))>
                        "Previous residues"
                    </button>
                    <button
                        class="text-tool"
                        disabled=move ||report.with(|r|r.as_ref().and_then(|(_,r)|r.tracks.get(candidate.get())).is_none_or(|t|(page.get()+1)*60>=t.residues.len()))
                        on:click=move |_|page.update(|p|*p+=1)>
                        "Next residues"
                    </button>
                    <button
                        class="button secondary full"
                        disabled=move ||workspace_busy.get()
                        on:click=move |_|{let result=(||{let (c,r)=report.get_untracked().ok_or("Analyze the source first.")?;let (o,l)=settings.get_untracked()?;c.guard(epoch.get_untracked(),&source.get_untracked().binding,&o,&l, &choice.get_untracked())?;serde_json::to_vec_pretty(r.as_ref()).map_err(|e|e.to_string())})();match result {Ok(bytes)=>crate::io::download::download_with_completion("antibody-analysis.json","application/json",&bytes,move |r|if status.try_get_untracked().is_some(){if let Err(e)=r{status.set(e);}}),Err(e)=>status.set(e)}}>
                        "Download full analysis"
                    </button>
                </Show>
            </details>
        </section>
    }
}
