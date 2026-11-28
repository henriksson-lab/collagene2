//! Revision-bound physical end repair and explicit oligo-annealing candidates.
use crate::platform::worker;
use crate::restriction::material_derivation_model::AnnealedReview;
use crate::restriction::material_derivation_model::SourcePin;
use collagene_core::{
    chemistry::{self, AnnealOptions, Duplex, EndRepairMethod, FivePrime, Terminus, ThreePrime},
    documents::transfer::Candidate,
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification,
};
use collagene_model::Molecule;
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
    source_info: Callback<usize, Option<crate::sources::identity::SourceInfo>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    on_material: Callback<
        crate::restriction::digest_product::MaterialOpenRequest,
        Result<(), String>,
    >,
) -> impl IntoView {
    let mode = RwSignal::new("repair");
    let source_index = RwSignal::new(0usize);
    let terminal_action = RwSignal::new(String::new());
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
    let repair = RwSignal::new(None::<Arc<verification::terminal::RepairPreview>>);
    let anneal = RwSignal::new(None::<Arc<AnnealedReview>>);
    let supplied_ssdna = RwSignal::new(false);
    let material_name = RwSignal::new("Annealed oligo material".to_string());
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
            supplied_ssdna.get(),
            epoch.get(),
            source_info.run(source_index.get()),
            source_info.run(bottom_index.get()),
        );
        generation.update(|id| *id = id.wrapping_add(1));
        repair.set(None);
        anneal.set(None);
        status.set(String::new());
        busy.set(false);
        page.set(0);
    });
    let anneal_options = Callback::new(move |()| {
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
    });
    let guard_alive = alive.clone();
    let guard_anneal = Callback::new(move |review: Arc<AnnealedReview>| -> Result<(), String> {
        if !guard_alive.load(Ordering::Relaxed) {
            return Err("The oligo material review was closed.".into());
        }
        if mode.get_untracked() != "anneal" {
            return Err("The selected physical operation changed.".into());
        }
        let top = source
            .get_untracked()
            .ok_or("The top oligo source is unavailable.")?;
        let bottom = bottom
            .get_untracked()
            .ok_or("The bottom oligo source is unavailable.")?;
        let top_info = source_info
            .run(source_index.get_untracked())
            .ok_or("Choose the exact top oligo event.")?;
        let bottom_info = source_info
            .run(bottom_index.get_untracked())
            .ok_or("Choose the exact bottom oligo event.")?;
        review.validate_current(
            &get_workspace.run(())?,
            &top_info.binding,
            &bottom_info.binding,
            epoch.get_untracked(),
            &top,
            &bottom,
            &anneal_options.run(())?,
            supplied_ssdna.get_untracked(),
        )
    });
    let run_alive = alive.clone();
    let run = move |_| {
        let Some(top) = source.get_untracked() else {
            status.set("Choose an available source document.".into());
            return;
        };
        let source_revision = top.revision();
        let mut anneal_context = None;
        let target_revision = molecule.with_untracked(Molecule::revision);
        let (revision, job) = if mode.get_untracked() == "repair" {
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
            (
                source_revision.clone(),
                AnalysisJob::RepairPreview {
                    molecule: top,
                    enzymes,
                    fragment_index: index,
                    method: method.get_untracked(),
                },
            )
        } else {
            let Some(bottom) = bottom.get_untracked() else {
                status.set("Choose the bottom oligo document.".into());
                return;
            };
            let options = anneal_options.run(());
            let options = match options {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let prepared = (|| -> Result<_, String> {
                if !supplied_ssdna.get_untracked() {
                    return Err("Declare both supplied strands as standard unmodified single-strand oligos before previewing annealing.".into());
                }
                let frame = get_workspace.run(())?;
                let top_info = source_info
                    .run(source_index.get_untracked())
                    .ok_or("Choose the exact top oligo event.")?;
                let bottom_info = source_info
                    .run(bottom_index.get_untracked())
                    .ok_or("Choose the exact bottom oligo event.")?;
                let top_pin = SourcePin::capture(
                    top_info.binding.clone(),
                    epoch.get_untracked(),
                    &top,
                    None,
                )?;
                let bottom_pin = SourcePin::capture(
                    bottom_info.binding.clone(),
                    epoch.get_untracked(),
                    &bottom,
                    None,
                )?;
                let revision =
                    collagene_core::material::derivation::annealed_oligos_preview_revision(
                        &frame,
                        &top_info.binding,
                        &bottom_info.binding,
                        &options,
                    )
                    .map_err(|e| e.to_string())?;
                Ok((
                    revision,
                    frame,
                    top_pin,
                    bottom_pin,
                    top_info.binding,
                    bottom_info.binding,
                ))
            })();
            let (revision, frame, top_pin, bottom_pin, top_event, bottom_event) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            anneal_context = Some((frame.clone(), top_pin, bottom_pin, options.clone()));
            (
                revision,
                AnalysisJob::InspectMaterialAnnealedOligos {
                    workspace: frame,
                    top: top_event,
                    bottom: bottom_event,
                    options,
                },
            )
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
                .is_none_or(|source| source.revision() != source_revision)
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
                        AnalysisResult::MaterialAnnealedOligosPreview(preview) => {
                            let reviewed = (|| -> Result<_, String> {
                                let (frame, top, bottom, options) = anneal_context
                                    .ok_or("The two-source oligo context is missing.")?;
                                let review = Arc::new(AnnealedReview::capture(
                                    &frame,
                                    top,
                                    bottom,
                                    &options,
                                    supplied_ssdna.get_untracked(),
                                    preview,
                                )?);
                                guard_anneal.run(review.clone())?;
                                Ok(review)
                            })();
                            match reviewed {
                                Ok(review) => {
                                    let preview = &review.preview().preview;
                                    status.set(format!("{} annealing geometries · search complete {} · unique within this model {}. Review the selected duplex and alternatives before saving material.",preview.annealing.candidates.len(),preview.annealing.complete,preview.annealing.unique));
                                    anneal.set(Some(review));
                                }
                                Err(e) => status.set(e),
                            }
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
                options: Default::default(),
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
        if top.revision() != preview.preview().preview.top_revision
            || bottom.revision() != preview.preview().preview.bottom_revision
            || index >= preview.preview().preview.annealing.candidates.len()
        {
            status.set("Oligo source or candidate changed; reconstruct annealing first.".into());
            return;
        }
        if let Err(e) = guard_anneal.run(preview.clone()) {
            status.set(e);
            return;
        }
        let preview = &preview.preview().preview;
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
                options: Default::default(),
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
        let Some(review) = anneal.get_untracked() else {
            return;
        };
        let result = (|| -> Result<(), String> {
            guard_anneal.run(review.clone())?;
            let product = review
                .preview()
                .products
                .get(index)
                .ok_or("Choose a reported annealing candidate.")?
                .clone();
            on_products.run(vec![product])
        })();
        match result {Ok(())=>status.set("Opened the selected annealed reference preview. Its retained native material recipe is created only by the separate reusable-material action.".into()),Err(e)=>status.set(e)}
    });
    let material_alive = alive.clone();
    let open_material = Callback::new(move |index: usize| {
        let result = (|| -> Result<(), String> {
            let review = anneal
                .get_untracked()
                .ok_or("Preview both supplied oligos first.")?;
            guard_anneal.run(review.clone())?;
            let candidate = review
                .preview()
                .preview
                .annealing
                .candidates
                .get(index)
                .ok_or("Choose a reported annealing candidate.")?;
            let selected_name = material_name.get_untracked();
            let request = review.request(&get_workspace.run(())?, &candidate.id, &selected_name)?;
            let held = review.clone();
            let alive = material_alive.clone();
            let guard = Callback::new(move |()| -> Result<(), String> {
                if !alive.load(Ordering::Relaxed) {
                    return Err("The oligo material review was closed.".into());
                }
                guard_anneal.run(held.clone())?;
                if material_name.get_untracked() != selected_name
                    || !anneal.with_untracked(|current| {
                        current
                            .as_ref()
                            .is_some_and(|current| Arc::ptr_eq(current, &held))
                    })
                {
                    return Err(
                        "The selected oligo material or its name changed; review again.".into(),
                    );
                }
                Ok(())
            });
            on_material.run(crate::restriction::digest_product::MaterialOpenRequest {
                request,
                epoch: review.epoch(),
                guard,
            })
        })();
        if let Err(e) = result {
            status.set(e);
        }
    });
    view! {<div class="inspector-content end-repair-panel"><h2>"DNA ends and annealing"</h2><p class="helper">{move||if matches!(mode.get(),"saved"|"terminal"){"Reconstruct actual strands and end chemistry from the selected source, then review the reusable material product."}else{"Reconstruct actual strands and end chemistry. Assembly candidates retain the active document as their independently chosen target."}}</p>
        <label>"Physical operation"<select aria-label="End repair mode" prop:value=move||mode.get() on:change=move|event|mode.set(match event_target_value(&event).as_str(){"anneal"=>"anneal","saved"=>"saved","terminal"=>"terminal",_=>"repair"})><option value="saved">"Repair selected saved material directly"</option><option value="terminal">"Tail or treat source material directly"</option><option value="repair">"Restriction digest → end repair"</option><option value="anneal">"Anneal two oligos"</option></select></label>
        <label>{move||if mode.get()=="anneal"{"Top oligo (ordered 5′→3′)"}else if mode.get()=="saved"{"Saved source material"}else if mode.get()=="terminal"{"Source material"}else{"Source before restriction digest"}}<select aria-label="End repair source" prop:value=move||source_index.get().to_string() on:change=move|event|source_index.set(event_target_value(&event).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(index,name)|view!{<option value=index.to_string()>{name}</option>}).collect_view()}</select></label>
        <Show when=move||mode.get()=="saved"><label>"Explicit end-repair method"<select aria-label="Saved material end repair method" prop:value=move||match method.get(){EndRepairMethod::Blunt=>"blunt",EndRepairMethod::FillFivePrimeOverhangs=>"fill",EndRepairMethod::RemoveThreePrimeOverhangs=>"chew"} on:change=move|e|method.set(match event_target_value(&e).as_str(){"fill"=>EndRepairMethod::FillFivePrimeOverhangs,"chew"=>EndRepairMethod::RemoveThreePrimeOverhangs,_=>EndRepairMethod::Blunt})><option value="blunt">"Blunt: fill 5′ overhangs / remove 3′ overhangs"</option><option value="fill">"Fill recessed 3′ ends from 5′ overhangs"</option><option value="chew">"Remove protruding 3′ ends"</option></select></label>
            <crate::restriction::reusable_material_view::ReusableMaterialView source source_info=Memo::new(move|_|source_info.run(source_index.get())) settings=Callback::new(move|()|Ok(crate::restriction::reusable_material_view::Settings::Repair{method:method.get_untracked()})) settings_revision=Memo::new(move|_|format!("{:?}",method.get())) kind=crate::restriction::reusable_material_view::Kind::EndRepair origin=RwSignal::new(None) epoch get_workspace on_open=on_material/>
        </Show>
        <Show when=move||mode.get()=="terminal">
            <label>"Explicit terminal action"<select aria-label="Reusable material terminal action" prop:value=move||terminal_action.get() on:change=move|e|terminal_action.set(event_target_value(&e))>
                <option value="">"Choose an action…"</option>
                <option value="phosphorylate">"Phosphorylate accessible 5′ ends"</option>
                <option value="dephosphorylate">"Dephosphorylate accessible 5′ ends"</option>
                <option value="tail-a">"Add one A to both 3′ ends (dA tailing)"</option>
                <option value="tail-t">"Add one T to both 3′ ends"</option>
                <option value="tail-g">"Add one G to both 3′ ends"</option>
                <option value="tail-c">"Add one C to both 3′ ends"</option>
            </select></label>
            <p class="helper">"Only the selected terminal action is performed. Review the actual strands and ends; tailing requires a blunt duplex with known free 3′ OH termini. This records source-derived material without selecting a complementary partner or an assembly target."</p>
            <crate::restriction::reusable_material_view::ReusableMaterialView source source_info=Memo::new(move|_|source_info.run(source_index.get()))
                settings=Callback::new(move|()|{
                    use collagene_core::{material::derivation::TerminalAction, chemistry::EndTreatment, cloning::ligation::tailing::TailingBase};
                    let action=match terminal_action.get_untracked().as_str() {
                        "phosphorylate"=>TerminalAction::TreatEnds{treatment:EndTreatment::Phosphorylate},
                        "dephosphorylate"=>TerminalAction::TreatEnds{treatment:EndTreatment::Dephosphorylate},
                        "tail-a"=>TerminalAction::Tail3Prime{base:TailingBase::A},
                        "tail-t"=>TerminalAction::Tail3Prime{base:TailingBase::T},
                        "tail-g"=>TerminalAction::Tail3Prime{base:TailingBase::G},
                        "tail-c"=>TerminalAction::Tail3Prime{base:TailingBase::C},
                        _=>return Err("Choose the explicit terminal action before previewing material.".into()),
                    };
                    Ok(crate::restriction::reusable_material_view::Settings::TerminalProcessing{action})
                }) settings_revision=Memo::new(move|_|terminal_action.get()) kind=crate::restriction::reusable_material_view::Kind::TerminalProcessing origin=RwSignal::new(None) epoch get_workspace on_open=on_material/>
        </Show>
        <Show when=move||!matches!(mode.get(),"saved"|"terminal")>
        <Show when=move||mode.get()=="repair"><label>"Restriction enzymes (comma separated)"<input aria-label="End repair enzymes" prop:value=move||enzymes.get() on:input=move|event|enzymes.set(event_target_value(&event))/></label><label>"Digest fragment number"<input aria-label="End repair fragment" type="number" min="1" prop:value=move||fragment.get() on:input=move|event|fragment.set(event_target_value(&event))/></label><p class="helper">"Fragment numbers follow ascending top-reference cuts, not fragment size. Imported physical-intermediate metadata cannot establish the starting duplex."</p><label>"Repair method"<select aria-label="End repair method" prop:value=move||match method.get(){EndRepairMethod::Blunt=>"blunt",EndRepairMethod::FillFivePrimeOverhangs=>"fill",EndRepairMethod::RemoveThreePrimeOverhangs=>"chew"} on:change=move|event|method.set(match event_target_value(&event).as_str(){"fill"=>EndRepairMethod::FillFivePrimeOverhangs,"chew"=>EndRepairMethod::RemoveThreePrimeOverhangs,_=>EndRepairMethod::Blunt})><option value="blunt">"Blunt: fill 5′ overhangs / remove 3′ overhangs"</option><option value="fill">"Fill recessed 3′ ends from 5′ overhangs"</option><option value="chew">"Remove protruding 3′ ends"</option></select></label></Show>
        <Show when=move||mode.get()=="anneal"><label>"Bottom oligo (ordered 5′→3′)"<select aria-label="Annealing bottom oligo" prop:value=move||bottom_index.get().to_string() on:change=move|event|bottom_index.set(event_target_value(&event).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(index,name)|view!{<option value=index.to_string()>{name}</option>}).collect_view()}</select></label><label>"Minimum paired bases"<input aria-label="Annealing minimum paired bases" type="number" min="1" prop:value=move||minimum.get() on:input=move|event|minimum.set(event_target_value(&event))/></label><label>"Comparison budget"<input aria-label="Annealing comparison budget" type="number" prop:value=move||budget.get() on:input=move|event|budget.set(event_target_value(&event))/></label><label>"Maximum alternatives"<input aria-label="Annealing maximum alternatives" type="number" prop:value=move||maximum.get() on:input=move|event|maximum.set(event_target_value(&event))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Top oligo ordered phosphorylated" prop:checked=move||top_phosphate.get() on:change=move|event|top_phosphate.set(event_target::<HtmlInputElement>(&event).checked())/>"Top oligo ordered with 5′ phosphate"</label><label class="checkbox-label"><input type="checkbox" aria-label="Bottom oligo ordered phosphorylated" prop:checked=move||bottom_phosphate.get() on:change=move|event|bottom_phosphate.set(event_target::<HtmlInputElement>(&event).checked())/>"Bottom oligo ordered with 5′ phosphate"</label><label class="checkbox-label"><input aria-label="Declare supplied single-strand oligos" type="checkbox" prop:checked=move||supplied_ssdna.get() on:change=move|e|supplied_ssdna.set(event_target_checked(&e))/>"Both selected references describe supplied standard unmodified single-strand oligos"</label><p class="helper">"Declare the actual oligos and 5′ groups explicitly. Both 3′ groups are free OH under this selected input hypothesis. Saved duplex material cannot supply a strand by implication."</p><p class="helper">"Both oligos are supplied as linear single-strand material with free 3′ OH. Phosphates are explicit input declarations. Only concrete complementary overlap is modeled; no mismatch repair, hairpin, self-association or thermodynamic ranking."</p></Show>
        <button class="button primary full" disabled=move||busy.get() on:click=run.clone()>"Compute physical intermediates"</button><p class="end-repair-status" role="status">{move||status.get()}</p>
        <h3>"Candidate plan options"</h3><label class="checkbox-label"><input type="checkbox" aria-label="Phosphorylate physical product" prop:checked=move||phosphorylate.get() on:change=move|event|phosphorylate.set(event_target::<HtmlInputElement>(&event).checked())/>"Add a terminal-phosphorylation step"</label><label class="checkbox-label"><input type="checkbox" aria-label="Circular ligation after physical operation" prop:checked=move||close.get() on:change=move|event|close.set(event_target::<HtmlInputElement>(&event).checked())/>"Add circular ligation with no retained nicks"</label>
        {move||repair.get().map(|preview|{let full=preview.clone();let opened=preview.clone();view!{<section class="end-repair-report"><h3>"Repaired duplex"</h3><pre>{describe(&preview.repair.duplex)}</pre><p>{format!("Top 3′ added '{}' / removed '{}'; bottom 3′ added '{}' / removed '{}'",preview.repair.top_added_3prime,preview.repair.top_removed_3prime,preview.repair.bottom_added_3prime,preview.repair.bottom_removed_3prime)}</p><details><summary>"Source, ancestry and assumptions"</summary><pre>{describe(&preview.repair.source)}</pre><ul>{preview.findings.iter().chain(&preview.repair.warnings).cloned().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details><button class="button primary full" on:click=repair_candidate>"Review repaired-product candidate"</button><button class="button secondary full" on:click=move|_|match on_products.run(vec![opened.molecule.clone()]){Ok(())=>status.set("Repaired reference preview opened; reconstruct its source graph before later physical operations.".into()),Err(error)=>status.set(error)}>"Open repaired reference preview"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(full.as_ref()){Ok(bytes)=>crate::io::download::download_with_completion("end-repair.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}>"Download end-repair evidence"</button></section>}})}
        {move||anneal.get().map(|review|{let preview=review.preview().preview.clone();let full=Arc::new(review.preview().clone());let count=preview.annealing.candidates.len();let start=page.get().min(count.saturating_sub(1)/50)*50;view!{<section class="annealing-report"><h3>{format!("{} annealing candidates",count)}</h3><label>"Reusable product name"<input maxlength="256" aria-label="Annealed material name" prop:value=move||material_name.get() on:input=move|e|material_name.set(event_target_value(&e))/></label><p>{format!("Search complete {} · unique {} · {} comparisons",preview.annealing.complete,preview.annealing.unique,preview.annealing.comparisons)}</p><details><summary>"Annealing assumptions and limits"</summary><ul>{preview.annealing.warnings.iter().cloned().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details>{preview.annealing.candidates.iter().enumerate().skip(start).take(50).map(|(index,candidate)|view!{<article class="annealing-candidate"><h4>{format!("Candidate {} · {} paired bases",index+1,candidate.paired_bases)}</h4><pre>{describe(&candidate.duplex)}</pre><button class="button primary full" on:click=move|_|anneal_candidate.run(index)>"Review selected annealing candidate"</button><button class="button secondary full" on:click=move|_|open_annealed.run(index)>"Open annealed reference preview"</button><button class="button primary full" on:click=move|_|open_material.run(index)>"Review reusable annealed material"</button></article>}).collect_view()}<div class="annotation-pages"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|page|*page=page.saturating_sub(1))>"Previous annealing candidates"</button><button class="button secondary" disabled=move||{(page.get()+1)*50>=count} on:click=move|_|page.update(|page|*page+=1)>"Next annealing candidates"</button></div><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(full.as_ref()){Ok(bytes)=>crate::io::download::download_with_completion("oligo-annealing.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}>"Download annealing evidence"</button></section>}})}
        </Show>
    </div>}
}
