//! Revision-bound physical end repair and explicit oligo-annealing candidates.
use crate::worker;
use collagene_core::{
    Molecule,
    chemistry::{self, AnnealOptions, Duplex, EndRepairMethod, FivePrime, Terminus, ThreePrime},
    protocol::{AnalysisJob, AnalysisResult, Request},
    transfer::Candidate,
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source},
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;
use web_sys::HtmlInputElement;
fn short(sequence: &str) -> String {
    let preview = sequence.chars().take(600).collect::<String>();
    if sequence.len() > 600 {
        format!("{preview}… (download retains the full sequence)")
    } else {
        preview
    }
}
fn describe(duplex: &Duplex) -> String {
    format!(
        "Top 5′→3′: {}\nBottom 5′→3′: {}\n{}",
        short(duplex.top.as_str()),
        short(duplex.bottom.as_str()),
        duplex
            .summary()
            .map(|summary| summary.describe())
            .unwrap_or_else(|error| error.to_string())
    )
}
#[component]
pub fn EndRepairPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
) -> impl IntoView {
    let mode = RwSignal::new("repair");
    let source_index = RwSignal::new(0usize);
    let bottom_index = RwSignal::new(0usize);
    let source = Memo::new(move |_| get_query.run(source_index.get()));
    let bottom = Memo::new(move |_| get_query.run(bottom_index.get()));
    let enzymes = RwSignal::new("EcoRI".to_string());
    let fragment = RwSignal::new("1".to_string());
    let method = RwSignal::new(EndRepairMethod::Blunt);
    let minimum = RwSignal::new("12".to_string());
    let budget = RwSignal::new("20000000".to_string());
    let maximum = RwSignal::new("100".to_string());
    let top_phosphate = RwSignal::new(false);
    let bottom_phosphate = RwSignal::new(false);
    let close = RwSignal::new(false);
    let phosphorylate = RwSignal::new(false);
    let repair = RwSignal::new(None::<Arc<verification::RepairPreview>>);
    let anneal = RwSignal::new(None::<Arc<verification::AnnealPreview>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let page = RwSignal::new(0usize);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    Effect::new(move |_| {
        let _ = (
            molecule.with(Molecule::revision),
            source.get().map(|m| m.revision()),
            bottom.get().map(|m| m.revision()),
            mode.get(),
            source_index.get(),
            bottom_index.get(),
            enzymes.get(),
            fragment.get(),
            method.get(),
            minimum.get(),
            budget.get(),
            maximum.get(),
            top_phosphate.get(),
            bottom_phosphate.get(),
        );
        generation.update(|id| *id = id.wrapping_add(1));
        repair.set(None);
        anneal.set(None);
        status.set(String::new());
        busy.set(false);
        page.set(0);
    });
    let run_alive = alive.clone();
    let run = move |_| {
        let Some(top) = source.get_untracked() else {
            status.set("Choose an available source document.".into());
            return;
        };
        let revision = top.revision();
        let target_revision = molecule.with_untracked(Molecule::revision);
        let job = if mode.get_untracked() == "repair" {
            let enzymes = enzymes
                .get_untracked()
                .split(',')
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>();
            let Some(index) = fragment
                .get_untracked()
                .parse::<usize>()
                .ok()
                .and_then(|n| n.checked_sub(1))
            else {
                status.set("Digest fragment number must start at one.".into());
                return;
            };
            if enzymes.is_empty() {
                status.set("Supply at least one restriction enzyme so the physical ends are reconstructed from a digest.".into());
                return;
            }
            AnalysisJob::RepairPreview {
                molecule: top,
                enzymes,
                fragment_index: index,
                method: method.get_untracked(),
            }
        } else {
            let Some(bottom) = bottom.get_untracked() else {
                status.set("Choose the bottom oligo document.".into());
                return;
            };
            let options = (|| {
                Ok::<_, String>(AnnealOptions {
                    min_paired_bases: minimum
                        .get_untracked()
                        .parse()
                        .map_err(|_| "Invalid paired-base minimum")?,
                    max_comparisons: budget
                        .get_untracked()
                        .parse()
                        .map_err(|_| "Invalid comparison budget")?,
                    max_candidates: maximum
                        .get_untracked()
                        .parse()
                        .map_err(|_| "Invalid candidate limit")?,
                    top_chemistry: Terminus {
                        five_prime: if top_phosphate.get_untracked() {
                            FivePrime::Phosphate
                        } else {
                            FivePrime::Hydroxyl
                        },
                        three_prime: ThreePrime::Hydroxyl,
                    },
                    bottom_chemistry: Terminus {
                        five_prime: if bottom_phosphate.get_untracked() {
                            FivePrime::Phosphate
                        } else {
                            FivePrime::Hydroxyl
                        },
                        three_prime: ThreePrime::Hydroxyl,
                    },
                })
            })();
            let options = match options {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            AnalysisJob::AnnealPreview {
                molecule: top,
                bottom,
                options,
            }
        };
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        status.set("Reconstructing both strands in a Rust worker…".into());
        repair.set(None);
        anneal.set(None);
        let alive = run_alive.clone();
        spawn_local(async move {
            let result = worker::analyze(Request {
                id: token,
                revision: revision.clone(),
                job,
            })
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            if source
                .get_untracked()
                .is_none_or(|source| source.revision() != revision)
                || molecule.with_untracked(Molecule::revision) != target_revision
            {
                return;
            }
            busy.set(false);
            match result {
                Ok(response) if response.id == token && response.revision == revision => {
                    match response.result {
                        AnalysisResult::RepairPreview(preview)
                            if preview.source_revision == revision =>
                        {
                            status.set("End geometry reconstructed. Inspect added/removed bases and terminal chemistry; this preview does not verify the active target.".into());
                            repair.set(Some(Arc::new(preview)));
                        }
                        AnalysisResult::AnnealPreview(preview)
                            if preview.top_revision == revision
                                && bottom.get_untracked().is_some_and(|bottom| {
                                    bottom.revision() == preview.bottom_revision
                                }) =>
                        {
                            status.set(format!("{} annealing geometries · search complete {} · unique within this model {}. Inspect alternatives; no stability or experimental-success claim.",preview.annealing.candidates.len(),preview.annealing.complete,preview.annealing.unique));
                            anneal.set(Some(Arc::new(preview)));
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Unexpected physical-preview response.".into()),
                    }
                }
                Ok(_) => status.set("Physical-preview source binding changed.".into()),
                Err(error) => status.set(error),
            }
        });
    };
    let repair_candidate = move |_| {
        let (Some(preview), Some(input)) = (repair.get_untracked(), source.get_untracked()) else {
            return;
        };
        if input.revision() != preview.source_revision {
            status.set("Source changed; reconstruct the repair first.".into());
            return;
        }
        let target = molecule.get_untracked();
        let mut steps = vec![
            Operation::Digest {
                template: Source::Input("source".into()),
                enzymes: preview.enzymes.clone(),
                fragment_index: preview.fragment_index,
            },
            Operation::RepairEnds {
                template: Source::Step(0),
                method: preview.method,
            },
        ];
        let mut product = 1;
        if phosphorylate.get_untracked() {
            steps.push(Operation::TreatEnds {
                template: Source::Step(product),
                treatment: chemistry::EndTreatment::Phosphorylate,
            });
            product = steps.len() - 1;
        }
        if close.get_untracked() {
            steps.push(Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(product),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
            });
        }
        on_candidate.run(Candidate::new(DesignPlan {
            name: format!("Digest / end repair to {}", target.name),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "source".into(),
                revision: input.revision(),
                molecule: input,
            }],
            steps,
        }));
    };
    let anneal_candidate = Callback::new(move |index: usize| {
        let (Some(preview), Some(top), Some(bottom)) = (
            anneal.get_untracked(),
            source.get_untracked(),
            bottom.get_untracked(),
        ) else {
            return;
        };
        if top.revision() != preview.top_revision
            || bottom.revision() != preview.bottom_revision
            || index >= preview.annealing.candidates.len()
        {
            status.set("Oligo source or candidate changed; reconstruct annealing first.".into());
            return;
        }
        let target = molecule.get_untracked();
        let mut steps = vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: preview.annealing.options.clone(),
            candidate_index: index,
        }];
        let mut product = 0;
        if phosphorylate.get_untracked() {
            steps.push(Operation::TreatEnds {
                template: Source::Step(product),
                treatment: chemistry::EndTreatment::Phosphorylate,
            });
            product = steps.len() - 1;
        }
        if close.get_untracked() {
            steps.push(Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(product),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
            });
        }
        on_candidate.run(Candidate::new(DesignPlan {
            name: format!("Annealed oligos to {}", target.name),
            target_revision: target.revision(),
            inputs: vec![
                PlanInput {
                    id: "top".into(),
                    revision: top.revision(),
                    molecule: top,
                },
                PlanInput {
                    id: "bottom".into(),
                    revision: bottom.revision(),
                    molecule: bottom,
                },
            ],
            steps,
        }));
    });
    let open_annealed = Callback::new(move |index: usize| {
        let (Some(preview), Some(top), Some(bottom)) = (
            anneal.get_untracked(),
            source.get_untracked(),
            bottom.get_untracked(),
        ) else {
            return;
        };
        let revision = top.revision();
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        status.set("Rechecking annealing evidence in a worker…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let response = worker::analyze(Request {
                id: token,
                revision: revision.clone(),
                job: AnalysisJob::AnnealedDocument {
                    molecule: top,
                    bottom,
                    preview: preview.as_ref().clone(),
                    candidate_index: index,
                },
            })
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            busy.set(false);
            match response {
                Ok(response) if response.id==token&&response.revision==revision=>match response.result {
                    AnalysisResult::PhysicalPreviewDocument(document)=>match on_products.run(vec![document]){Ok(())=>status.set("Annealed reference preview opened. Reconstruct its originating plan for later physical operations.".into()),Err(error)=>status.set(error)},
                    AnalysisResult::Error(error)=>status.set(error),_=>status.set("Unexpected annealed-document response.".into()),
                },
                Ok(_)=>status.set("Annealed-document source binding changed.".into()),Err(error)=>status.set(error),
            }
        });
    });
    view! {<div class="inspector-content end-repair-panel"><h2>"End repair and annealing"</h2><p class="helper">"Reconstruct actual strands and end chemistry. The active document remains the independently chosen target; previews never choose or replace that target."</p>
        <label>"Physical operation"<select aria-label="End repair mode" prop:value=move||mode.get() on:change=move|event|mode.set(if event_target_value(&event)=="anneal"{"anneal"}else{"repair"})><option value="repair">"Restriction digest → end repair"</option><option value="anneal">"Anneal two oligos"</option></select></label>
        <label>{move||if mode.get()=="anneal"{"Top oligo (ordered 5′→3′)"}else{"Source before restriction digest"}}<select aria-label="End repair source" prop:value=move||source_index.get().to_string() on:change=move|event|source_index.set(event_target_value(&event).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(index,name)|view!{<option value=index.to_string()>{name}</option>}).collect_view()}</select></label>
        <Show when=move||mode.get()=="repair"><label>"Restriction enzymes (comma separated)"<input aria-label="End repair enzymes" prop:value=move||enzymes.get() on:input=move|event|enzymes.set(event_target_value(&event))/></label><label>"Digest fragment number"<input aria-label="End repair fragment" type="number" min="1" prop:value=move||fragment.get() on:input=move|event|fragment.set(event_target_value(&event))/></label><p class="helper">"Fragment numbers follow ascending top-reference cuts, not fragment size. Imported physical-intermediate metadata cannot establish the starting duplex."</p><label>"Repair method"<select aria-label="End repair method" prop:value=move||match method.get(){EndRepairMethod::Blunt=>"blunt",EndRepairMethod::FillFivePrimeOverhangs=>"fill",EndRepairMethod::RemoveThreePrimeOverhangs=>"chew"} on:change=move|event|method.set(match event_target_value(&event).as_str(){"fill"=>EndRepairMethod::FillFivePrimeOverhangs,"chew"=>EndRepairMethod::RemoveThreePrimeOverhangs,_=>EndRepairMethod::Blunt})><option value="blunt">"Blunt: fill 5′ overhangs / remove 3′ overhangs"</option><option value="fill">"Fill recessed 3′ ends from 5′ overhangs"</option><option value="chew">"Remove protruding 3′ ends"</option></select></label></Show>
        <Show when=move||mode.get()=="anneal"><label>"Bottom oligo (ordered 5′→3′)"<select aria-label="Annealing bottom oligo" prop:value=move||bottom_index.get().to_string() on:change=move|event|bottom_index.set(event_target_value(&event).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(index,name)|view!{<option value=index.to_string()>{name}</option>}).collect_view()}</select></label><label>"Minimum paired bases"<input aria-label="Annealing minimum paired bases" type="number" min="1" prop:value=move||minimum.get() on:input=move|event|minimum.set(event_target_value(&event))/></label><label>"Comparison budget"<input aria-label="Annealing comparison budget" type="number" prop:value=move||budget.get() on:input=move|event|budget.set(event_target_value(&event))/></label><label>"Maximum alternatives"<input aria-label="Annealing maximum alternatives" type="number" prop:value=move||maximum.get() on:input=move|event|maximum.set(event_target_value(&event))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Top oligo ordered phosphorylated" prop:checked=move||top_phosphate.get() on:change=move|event|top_phosphate.set(event_target::<HtmlInputElement>(&event).checked())/>"Top oligo ordered with 5′ phosphate"</label><label class="checkbox-label"><input type="checkbox" aria-label="Bottom oligo ordered phosphorylated" prop:checked=move||bottom_phosphate.get() on:change=move|event|bottom_phosphate.set(event_target::<HtmlInputElement>(&event).checked())/>"Bottom oligo ordered with 5′ phosphate"</label><p class="helper">"Both oligos are supplied as linear single-strand material with free 3′ OH. Phosphates are explicit input declarations. Only concrete complementary overlap is modeled; no mismatch repair, hairpin, self-association or thermodynamic ranking."</p></Show>
        <button class="button primary full" disabled=move||busy.get() on:click=run>"Compute physical intermediates"</button><p class="end-repair-status" role="status">{move||status.get()}</p>
        <h3>"Candidate plan options"</h3><label class="checkbox-label"><input type="checkbox" aria-label="Phosphorylate physical product" prop:checked=move||phosphorylate.get() on:change=move|event|phosphorylate.set(event_target::<HtmlInputElement>(&event).checked())/>"Add a terminal-phosphorylation step"</label><label class="checkbox-label"><input type="checkbox" aria-label="Circular ligation after physical operation" prop:checked=move||close.get() on:change=move|event|close.set(event_target::<HtmlInputElement>(&event).checked())/>"Add circular ligation with no retained nicks"</label>
        {move||repair.get().map(|preview|{let full=preview.clone();let opened=preview.clone();view!{<section class="end-repair-report"><h3>"Repaired duplex"</h3><pre>{describe(&preview.repair.duplex)}</pre><p>{format!("Top 3′ added '{}' / removed '{}'; bottom 3′ added '{}' / removed '{}'",preview.repair.top_added_3prime,preview.repair.top_removed_3prime,preview.repair.bottom_added_3prime,preview.repair.bottom_removed_3prime)}</p><details><summary>"Source, ancestry and assumptions"</summary><pre>{describe(&preview.repair.source)}</pre><ul>{preview.findings.iter().chain(&preview.repair.warnings).cloned().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details><button class="button primary full" on:click=repair_candidate>"Review repaired-product candidate"</button><button class="button secondary full" on:click=move|_|match on_products.run(vec![opened.molecule.clone()]){Ok(())=>status.set("Repaired reference preview opened; reconstruct its source graph before later physical operations.".into()),Err(error)=>status.set(error)}>"Open repaired reference preview"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(full.as_ref()){Ok(bytes)=>if let Err(error)=crate::app::download("end-repair.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}>"Download end-repair evidence"</button></section>}})}
        {move||anneal.get().map(|preview|{let full=preview.clone();let count=preview.annealing.candidates.len();let start=page.get().min(count.saturating_sub(1)/50)*50;view!{<section class="annealing-report"><h3>{format!("{} annealing candidates",count)}</h3><p>{format!("Search complete {} · unique {} · {} comparisons",preview.annealing.complete,preview.annealing.unique,preview.annealing.comparisons)}</p><details><summary>"Annealing assumptions and limits"</summary><ul>{preview.annealing.warnings.iter().cloned().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details>{preview.annealing.candidates.iter().enumerate().skip(start).take(50).map(|(index,candidate)|view!{<article class="annealing-candidate"><h4>{format!("Candidate {} · {} paired bases",index+1,candidate.paired_bases)}</h4><pre>{describe(&candidate.duplex)}</pre><button class="button primary full" on:click=move|_|anneal_candidate.run(index)>"Review selected annealing candidate"</button><button class="button secondary full" on:click=move|_|open_annealed.run(index)>"Open annealed reference preview"</button></article>}).collect_view()}<div class="annotation-pages"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|page|*page=page.saturating_sub(1))>"Previous annealing candidates"</button><button class="button secondary" disabled=move||{(page.get()+1)*50>=count} on:click=move|_|page.update(|page|*page+=1)>"Next annealing candidates"</button></div><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(full.as_ref()){Ok(bytes)=>if let Err(error)=crate::app::download("oligo-annealing.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}>"Download annealing evidence"</button></section>}})}
    </div>}
}
