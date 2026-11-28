//! Collapsed restriction-alternative review beneath the existing IVT declaration.
use super::{
    linearization_model::{Capture, Draft, Snapshot},
    model::Draft as IvtDraft,
};
use crate::{
    cloning::methods::checked_cloning::{AssumptionSignals, SourceAssumptions},
    history::run_origin_model::PinnedOrigin,
    platform::worker::CachedWorker,
    restriction::{reaction_conditions::ReactionConditionsPanel, reaction_conditions_model},
    sources::identity::SourceInfo,
};
use collagene_core::{
    cloning::restriction::{catalog, conditions::ReactionConditions, ivt_linearization as native},
    documents::workspace::Workspace,
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

fn toggle(
    draft: RwSignal<Draft>,
    label: &'static str,
    get: fn(&Draft) -> bool,
    set: fn(&mut Draft, bool),
) -> impl IntoView {
    view! {
        <label class="checkbox-label">
            <input type="checkbox" prop:checked=move || draft.with(get)
                on:change=move |event| draft.update(|d| set(d, event_target_checked(&event)))/>
            {label}
        </label>
    }
}
fn limit(
    draft: RwSignal<Draft>,
    label: &'static str,
    get: fn(&Draft) -> &str,
    set: fn(&mut Draft, String),
) -> impl IntoView {
    view! {
        <label>{label}
            <input inputmode="numeric" maxlength="20"
                prop:value=move || draft.with(|d| get(d).to_string())
                on:input=move |event| draft.update(|d| set(d, event_target_value(&event)))/>
        </label>
    }
}

#[component]
pub(crate) fn IvtLinearizationPanel(
    molecule: Memo<Molecule>,
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    ivt_draft: RwSignal<IvtDraft>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
    let assumptions = AssumptionSignals::new();
    let origin = RwSignal::new(None::<PinnedOrigin>);
    let condition_draft = RwSignal::new(reaction_conditions_model::Draft::default());
    let conditions = RwSignal::new(None::<ReactionConditions>);
    let conditions_error = RwSignal::new(None::<String>);
    let catalog = StoredValue::new(
        catalog::builtin_catalog()
            .map(Arc::new)
            .map_err(|e| e.to_string()),
    );
    let selection = Memo::new(move |_| {
        let catalog = catalog.get_value()?;
        catalog::select(catalog.as_ref(), &draft.get().enzymes).map_err(|e| e.to_string())
    });
    let condition_source = Memo::new(move |_| Some(source.get()));
    let condition_origin = Memo::new(move |_| origin.get());
    let snapshot = Memo::new(move |_| {
        assumptions.track();
        Snapshot {
            binding: source.get().binding,
            epoch: epoch.get(),
            ivt: ivt_draft.get(),
            picker: draft.get(),
            generation: assumptions.current(),
            origin_revision: origin.get().as_ref().map(|pin| pin.revision().to_owned()),
            conditions: condition_draft.get(),
        }
    });
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<native::Report>)>);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get());
        origin.set(None);
        assumptions.host.set(None);
        assumptions
            .dam
            .set(collagene_core::chemistry::methylation::MethylationState::Unknown);
        assumptions
            .dcm
            .set(collagene_core::chemistry::methylation::MethylationState::Unknown);
        assumptions
            .cpg
            .set(collagene_core::chemistry::methylation::MethylationState::Unknown);
        condition_draft.set(reaction_conditions_model::Draft::default());
    });
    Effect::new(move |_| {
        let _ = snapshot.get();
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        status.set("Choose alternative enzymes and source assumptions, then inspect.".into());
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
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            if let Some(error) = conditions_error.get_untracked() {
                status.set(error);
                return;
            }
            let prepared = molecule.with_untracked(|m| {
                let selected = origin.get_untracked();
                let capture = Capture::prepare(
                    snapshot.get_untracked(),
                    m,
                    selected.as_ref().map(PinnedOrigin::as_ref),
                )?;
                Ok::<_, String>((m.clone(), capture))
            });
            let (molecule, capture) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let active = match CachedWorker::new() {
                Ok(value) => Rc::new(value),
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            worker.set_value(Some(active.clone()));
            busy.set(true);
            report.set(None);
            status.set("Inspecting alternative single-enzyme cut hypotheses…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let response = active
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::AssessIvtLinearization {
                            molecule,
                            request: capture.request.clone(),
                        },
                    })
                    .await;
                active.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                let received = (|| {
                    capture.guard(&snapshot.get_untracked())?;
                    if workspace_busy.get_untracked() {
                        return Err(
                            "The workspace changed while restriction alternatives were inspected."
                                .into(),
                        );
                    }
                    match response?.result {
                        AnalysisResult::IvtLinearization(value) => {
                            capture.validate_report(&value)?;
                            Ok(value)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected IVT linearization response.".into()),
                    }
                })();
                match received {
                    Ok(value) => {
                        report.set(Some((capture, Arc::new(value))));
                        status.set("Alternative evidence ready. Eligibility is conditional; no enzyme is selected automatically.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let export = move |_| {
        if workspace_busy.get_untracked() {
            return;
        }
        let Some((capture, value)) = report.get_untracked() else {
            return;
        };
        if let Err(error) = capture.guard(&snapshot.get_untracked()) {
            status.set(error);
            return;
        }
        let bytes = match serde_json::to_vec(value.as_ref()) {
            Ok(bytes) if bytes.len() <= native::MAX_REPORT_BYTES => bytes,
            Ok(_) => {
                status.set("The report exceeds its export bound.".into());
                return;
            }
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        wasm_bindgen_futures::spawn_local(async move {
            let result = crate::io::download::download_guarded(
                "ivt-linearization.json",
                "application/json",
                &bytes,
                move || {
                    if workspace_busy.try_get_untracked() != Some(false) {
                        return Err("The workspace is changing; export afterward.".into());
                    }
                    let current = snapshot
                        .try_get_untracked()
                        .ok_or("The IVT panel was closed.")?;
                    capture.guard(&current)
                },
            )
            .await;
            if let Err(error) = result {
                let _ = status.try_set(error);
            }
        });
    };
    view! {
        <details class="biology-section ivt-linearization">
            <summary>"Restriction linearization alternatives"</summary>
            <p class="helper">"Uses the transcript declaration above. Leave its supplied cut offsets empty. Each enzyme is a separate digestion hypothesis; promoter retention and exact runoff are checked without choosing a winner."</p>
            <crate::history::run_origin_view::SourceOriginPicker
                source=condition_source epoch get_workspace selected=condition_origin
                on_select=Callback::new(move |pin:Option<PinnedOrigin>| {
                    if pin.as_ref().is_some_and(|p|*p.event()!=source.get_untracked().binding) {
                        return Err("Saved material belongs to another exact source event.".into());
                    }
                    origin.set(pin);
                    Ok(())
                })/>
            <Show when=move || origin.with(Option::is_none)>
                <SourceAssumptions side="IVT template" signals=assumptions/>
            </Show>
            <Show when=move || origin.with(Option::is_some)>
                <p class="helper">"Saved material is freshly reconstructed from its native origin; host assumptions above do not replace its positional history."</p>
            </Show>
            <details>
                <summary>"Alternative enzymes"</summary>
                <button class="button secondary" on:click=move |_| {
                    match catalog.get_value() {
                        Ok(catalog) => draft.update(|d|d.enzymes=catalog.edition().enzymes.iter().map(|e|e.id.clone()).collect()),
                        Err(error) => status.set(error),
                    }
                }>"Select all installed enzymes"</button>
                <button class="button secondary" on:click=move |_| draft.update(|d|d.enzymes.clear())>"Clear enzymes"</button>
                <div class="digest-enzyme-picker">
                    {catalog.with_value(|value| value.as_ref().ok().map(|catalog| {
                        catalog.edition().enzymes.iter().map(|enzyme| {
                            let id = enzyme.id.clone();
                            let check = id.clone();
                            let name = format!("{} · {}",enzyme.name,enzyme.geometry.recognition);
                            view! {
                                <label class="checkbox-label">
                                    <input type="checkbox"
                                        prop:checked=move || draft.with(|d|d.enzymes.contains(&check))
                                        on:change=move |event| {
                                            let checked=event_target_checked(&event);
                                            draft.update(|d| {
                                                d.enzymes.retain(|old|old!=&id);
                                                if checked { d.enzymes.push(id.clone()); }
                                            });
                                        }/>
                                    {name}
                                </label>
                            }
                        }).collect_view()
                    }))}
                </div>
                <p class="helper">{move || format!("{} alternatives selected",draft.with(|d|d.enzymes.len()))}</p>
            </details>
            <ReactionConditionsPanel catalog selection source=condition_source epoch origin=condition_origin
                draft=Signal::from(condition_draft) on_change=Callback::new(move |value|condition_draft.set(value))
                conditions error=conditions_error/>
            <details>
                <summary>"Eligibility policy and limits"</summary>
                {toggle(draft,"Allow blunt runoff ends",|d|d.allow_blunt,|d,v|d.allow_blunt=v)}
                {toggle(draft,"Allow 5′ overhangs",|d|d.allow_five_prime,|d,v|d.allow_five_prime=v)}
                {toggle(draft,"Allow 3′ overhangs",|d|d.allow_three_prime,|d,v|d.allow_three_prime=v)}
                {toggle(draft,"Require known documentary reaction conditions",|d|d.require_conditions,|d,v|d.require_conditions=v)}
                {toggle(draft,"Require all IVT checks resolved",|d|d.require_resolved_ivt,|d,v|d.require_resolved_ivt=v)}
                {toggle(draft,"Require complete source history",|d|d.require_complete_origin,|d,v|d.require_complete_origin=v)}
                {limit(draft,"Shared recognition comparisons (maximum 100000000)",|d|&d.comparisons,|d,v|d.comparisons=v)}
                {limit(draft,"Retained sites per enzyme (maximum 64)",|d|&d.sites,|d,v|d.sites=v)}
                <p class="helper">"A limited search cannot establish absent cuts. Unknowns and multiple cuts remain visible. Star activity and incubation-duration compatibility are unassessed. Intact paired source duplexes only; existing nicks or staggered source ends are unsupported."</p>
            </details>
            <button class="button primary" disabled=move || busy.get()||workspace_busy.get()
                on:click=run>"Inspect restriction alternatives"</button>
            <p class="helper" role="status">{move || status.get()}</p>
            {move || report.get().map(|(_,report)| {
                view! {
                    <p class="helper">{format!("{} alternatives · {} shared comparisons · source history incomplete: {}",
                        report.candidates.len(),report.comparisons,report.source_history_incomplete)}</p>
                    {report.candidates.iter().map(|candidate| {
                        let mapping = candidate.mapping.as_ref().map(|m|format!(
                            "Source top/bottom bonds {}/{} · coding/template offsets {}/{} from +1 · signed stagger {}",
                            m.top_bond,m.bottom_bond,m.coding_cut_offset,m.template_cut_offset,m.signed_template_minus_coding,
                        ));
                        view! {
                            <details>
                                <summary>{format!("{} · {:?}",candidate.enzyme_id,candidate.state)}</summary>
                                <p class="helper">{format!("Scan complete: {} · products determined: {} · sites: {} · fragments: {}",
                                    candidate.cleavage.complete,candidate.cleavage.products_determined,
                                    candidate.cleavage.sites.len(),candidate.fragments.len())}</p>
                                <p class="helper">{mapping}</p>
                                {candidate.issues.iter().map(|issue|view! {
                                    <p class="helper">{format!("{:?}: {:?}",issue.severity,issue.kind)}</p>
                                }).collect_view()}
                                {candidate.cleavage.sites.iter().map(|site|view! {
                                    <p class="helper">{format!("{} recognition base {} · reverse {} · top/bottom bonds {:?}/{:?} · {:?}",
                                        site.enzyme,site.recognition_start+1,site.reverse,site.top_cut,site.bottom_cut,site.outcome)}</p>
                                }).collect_view()}
                                <details>
                                    <summary>"All predicted fragment geometries"</summary>
                                    {candidate.fragments.iter().map(|fragment|view! {
                                        <p class="helper">{format!(
                                            "Fragment {} · source top/bottom boundaries {}/{} · retained IVT template: {}. {}",
                                            fragment.index+1,fragment.source_top_start,fragment.source_bottom_start,
                                            candidate.retained_fragment==Some(fragment.index),fragment.geometry.describe(),
                                        )}</p>
                                    }).collect_view()}
                                </details>
                                {candidate.conditions.as_ref().map(|conditions|view! {
                                    <p class="helper">{format!("Documentary conditions: {:?} · source record {} · star/duration unassessed",
                                        conditions.report.documentary_status,conditions.context_revision)}</p>
                                })}
                                {candidate.ivt.as_ref().map(|ivt|view! {
                                    <div>{ivt.findings.iter().map(|finding|view! {
                                        <p class="helper">{format!("{:?} · {:?}: {}",finding.rule,finding.verdict,finding.message)}</p>
                                    }).collect_view()}
                                        <details><summary>"IVT evidence sources"</summary>
                                            {ivt.profile.sources.iter().map(|source|view! {
                                                <p class="helper"><a href=source.url.clone()
                                                    target="_blank" rel="noopener noreferrer">{source.id.clone()}</a>
                                                    {format!(" · {}",source.locator)}</p>
                                            }).collect_view()}
                                        </details>
                                    </div>
                                })}
                                {candidate.cleavage.limits.iter().map(|text|view! {<p class="helper">{text.clone()}</p>}).collect_view()}
                            </details>
                        }
                    }).collect_view()}
                    <p class="helper">{report.statement.clone()}</p>
                }
            })}
            <button class="button secondary" disabled=move || report.with(Option::is_none)||workspace_busy.get()
                on:click=export>"Download complete alternative evidence"</button>
        </details>
    }
}
