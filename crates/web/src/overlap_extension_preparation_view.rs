//! Compact manual oligo entry for the explicit two-branch preparation graph.
use crate::cloning_methods_model::{
    ExtensionPrimerInput, ExtensionPrimerPairInput, MethodsState, OverlapExtensionPcrInput,
};
use leptos::prelude::*;

fn pair(raw: &OverlapExtensionPcrInput, index: usize) -> &ExtensionPrimerPairInput {
    match index {
        1 => &raw.right,
        2 => &raw.outer,
        _ => &raw.left,
    }
}
fn pair_mut(raw: &mut OverlapExtensionPcrInput, index: usize) -> &mut ExtensionPrimerPairInput {
    match index {
        1 => &mut raw.right,
        2 => &mut raw.outer,
        _ => &mut raw.left,
    }
}
#[derive(Clone, Copy)]
enum PrimerField {
    Name,
    Annealing,
    Tail,
}
#[component]
fn PrimerFields(
    config: RwSignal<MethodsState>,
    index: usize,
    pair_index: usize,
    reverse: bool,
) -> impl IntoView {
    let direction = if reverse { "Reverse" } else { "Forward" };
    let current = move || {
        config.with(|c| {
            let raw = c
                .rows
                .get(index)
                .and_then(|r| r.overlap_extension_pcr.as_ref());
            let default = OverlapExtensionPcrInput::default();
            let p = pair(raw.unwrap_or(&default), pair_index);
            if reverse {
                p.reverse.clone()
            } else {
                p.forward.clone()
            }
        })
    };
    let update = Callback::new(move |(field, value): (PrimerField, String)| {
        config.update(|c| {
            if let Some(row) = c.rows.get_mut(index) {
                let raw = row
                    .overlap_extension_pcr
                    .get_or_insert_with(Default::default);
                let p = pair_mut(raw, pair_index);
                let primer: &mut ExtensionPrimerInput = if reverse {
                    &mut p.reverse
                } else {
                    &mut p.forward
                };
                match field {
                    PrimerField::Name => primer.name = value,
                    PrimerField::Annealing => primer.annealing = value,
                    PrimerField::Tail => primer.tail = value,
                }
            }
        });
    });
    let label = move |field: &str| {
        format!(
            "Fragment {} {} PCR {} {}",
            index + 1,
            ["left", "right", "outer"][pair_index.min(2)],
            direction.to_lowercase(),
            field
        )
    };
    view! {<fieldset class="extension-primer-fields">
        <legend>{format!("{direction} oligo · 5′→3′")}</legend>
        <label>"Name"<input maxlength="256" aria-label=move||label("name") prop:value=move||current().name on:input=move|e|update.run((PrimerField::Name,event_target_value(&e)))/></label>
        <label>"Annealing domain"<input maxlength="1000" spellcheck="false" autocomplete="off" aria-label=move||label("annealing domain") prop:value=move||current().annealing on:input=move|e|update.run((PrimerField::Annealing,event_target_value(&e)))/></label>
        <label>"Added 5′ tail"<input maxlength="1000" spellcheck="false" autocomplete="off" aria-label=move||label("five prime tail") prop:value=move||current().tail on:input=move|e|update.run((PrimerField::Tail,event_target_value(&e)))/></label>
    </fieldset>}
}
#[derive(Clone, Copy)]
enum Setting {
    Overlap,
    Minimum,
    Comparisons,
    Candidates,
}
#[component]
pub fn OverlapExtensionPreparationFields(
    config: RwSignal<MethodsState>,
    index: usize,
) -> impl IntoView {
    // This changes only which pair is visible; every raw oligo stays in the row.
    let active_pair = RwSignal::new(0usize);
    let raw = move || {
        config.with(|c| {
            c.rows
                .get(index)
                .and_then(|r| r.overlap_extension_pcr.clone())
                .unwrap_or_default()
        })
    };
    let update = Callback::new(move |(field, value): (Setting, String)| {
        config.update(|c| {
            if let Some(row) = c.rows.get_mut(index) {
                let raw = row
                    .overlap_extension_pcr
                    .get_or_insert_with(Default::default);
                match field {
                    Setting::Overlap => raw.overlap = value,
                    Setting::Minimum => raw.min_paired_bases = value,
                    Setting::Comparisons => raw.max_comparisons = value,
                    Setting::Candidates => raw.max_candidates = value,
                }
            }
        });
    });
    view! {<section class="overlap-extension-preparation" aria-label=format!("Fragment {} overlap-extension PCR",index+1)>
        <p class="helper">"Two PCRs use the same original source. Their products undergo reciprocal extension, followed by a separate outer-primer PCR. Enter every annealing domain and 5′ tail explicitly; your intended target stays separate."</p>
        <label>"Edit primer pair"<select aria-label=format!("Fragment {} extension primer pair",index+1) prop:value=move||active_pair.get().to_string() on:change=move|e|{if let Ok(pair)=event_target_value(&e).parse::<usize>(){if pair<3{active_pair.set(pair);}}}>
            <option value="0">"1 · Left source PCR"</option><option value="1">"2 · Right source PCR"</option><option value="2">"3 · Outer PCR after extension"</option>
        </select></label>
        <For each=move||vec![active_pair.get()] key=|pair|*pair children=move|pair_index|view!{<div class="primer-grid"><PrimerFields config index pair_index reverse=false/><PrimerFields config index pair_index reverse=true/></div>}/>
        <div class="primer-settings-grid">
            <label>"Selected terminal homology (bp)"<input type="number" min="1" aria-label=format!("Fragment {} selected extension homology",index+1) prop:value=move||raw().overlap on:input=move|e|update.run((Setting::Overlap,event_target_value(&e)))/></label>
            <label>"Minimum terminal pairing (bp)"<input type="number" min="1" aria-label=format!("Fragment {} minimum extension pairing",index+1) prop:value=move||raw().min_paired_bases on:input=move|e|update.run((Setting::Minimum,event_target_value(&e)))/></label>
        </div>
        <details><summary>"Terminal-pair search limits"</summary><div class="primer-settings-grid">
            <label>"Maximum comparisons"<input type="number" min="1" max="50000000" aria-label=format!("Fragment {} extension comparison limit",index+1) prop:value=move||raw().max_comparisons on:input=move|e|update.run((Setting::Comparisons,event_target_value(&e)))/></label>
            <label>"Maximum retained pairings"<input type="number" min="1" max="1000" aria-label=format!("Fragment {} extension candidate limit",index+1) prop:value=move||raw().max_candidates on:input=move|e|update.run((Setting::Candidates,event_target_value(&e)))/></label>
        </div><p class="helper">"Competing terminal pairings and capped searches remain incomplete evidence, including after outer PCR."</p></details>
        <details><summary>"Review all six ordered oligos"</summary><pre>{move||{
            let raw=raw();
            [&raw.left,&raw.right,&raw.outer].into_iter().zip(["Left source PCR","Right source PCR","Outer PCR"]).map(|(p,label)|format!("{label}\n{}: {}{}\n{}: {}{}",p.forward.name,p.forward.tail,p.forward.annealing,p.reverse.name,p.reverse.tail,p.reverse.annealing)).collect::<Vec<_>>().join("\n\n")
        }}</pre></details>
    </section>}
}
