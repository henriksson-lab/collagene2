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
    raw: Memo<OverlapExtensionPcrInput>,
    on_change: Callback<OverlapExtensionPcrInput>,
    context: String,
    pair_index: usize,
    reverse: bool,
) -> impl IntoView {
    let direction = if reverse { "Reverse" } else { "Forward" };
    let current = move || {
        let raw = raw.get();
        let p = pair(&raw, pair_index);
        if reverse {
            p.reverse.clone()
        } else {
            p.forward.clone()
        }
    };
    let update = Callback::new(move |(field, value): (PrimerField, String)| {
        let mut next = raw.get_untracked();
        let p = pair_mut(&mut next, pair_index);
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
        on_change.run(next);
    });
    let context = StoredValue::new(context);
    let label = move |field: &str| {
        format!(
            "{} {} PCR {} {}",
            context.get_value(),
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
    let raw = Memo::new(move |_| {
        config.with(|c| {
            c.rows
                .get(index)
                .and_then(|r| r.overlap_extension_pcr.clone())
                .unwrap_or_default()
        })
    });
    let on_change = Callback::new(move |next: OverlapExtensionPcrInput| {
        config.update(|c| {
            if let Some(row) = c.rows.get_mut(index) {
                row.overlap_extension_pcr = Some(next);
            }
        })
    });
    let allow_bottom = Memo::new(move |_| {
        config.with(|c| c.rows.get(index).is_some_and(|r| r.run_origin.is_some()))
    });
    view! {<OverlapExtensionOligoFields raw on_change allow_bottom context=format!("Fragment {}",index+1)/>}
}

/// The same six literal oligos serve assembly preparation and target-free material.
#[component]
pub fn OverlapExtensionOligoFields(
    raw: Memo<OverlapExtensionPcrInput>,
    on_change: Callback<OverlapExtensionPcrInput>,
    allow_bottom: Memo<bool>,
    context: String,
) -> impl IntoView {
    let active_pair = RwSignal::new(0usize);
    let context = StoredValue::new(context);
    let update = Callback::new(move |(field, value): (Setting, String)| {
        let mut next = raw.get_untracked();
        match field {
            Setting::Overlap => next.overlap = value,
            Setting::Minimum => next.min_paired_bases = value,
            Setting::Comparisons => next.max_comparisons = value,
            Setting::Candidates => next.max_candidates = value,
        }
        on_change.run(next);
    });
    view! {<section class="overlap-extension-preparation" aria-label=format!("{} overlap-extension PCR",context.get_value())>
        <p class="helper">"Two PCRs use the same exact source event and selected native strand frame. Saved Run/Derived material is freshly reconstructed before PCR; bottom-only calls are never converted into a new stock. Their products undergo reciprocal extension, followed by outer PCR on the synthesized product. Enter every annealing domain and 5′ tail explicitly; your intended target stays separate."</p>
        <label>"Edit primer pair"<select aria-label=format!("{} extension primer pair",context.get_value()) prop:value=move||active_pair.get().to_string() on:change=move|e|{if let Ok(pair)=event_target_value(&e).parse::<usize>(){if pair<3{active_pair.set(pair);}}}>
            <option value="0">"1 · Left source PCR"</option><option value="1">"2 · Right source PCR"</option><option value="2">"3 · Outer PCR after extension"</option>
        </select></label>
        <Show when={move||active_pair.get()<2}><label>"Source PCR template strand"<select aria-label=format!("{} extension source template strand",context.get_value()) prop:value=move||match pair(&raw.get(),active_pair.get()).template_strand{collagene_core::chemistry::Strand::Top=>"top",collagene_core::chemistry::Strand::Bottom=>"bottom"} on:change=move|event|{let strand=match event_target_value(&event).as_str(){"bottom"=>collagene_core::chemistry::Strand::Bottom,_=>collagene_core::chemistry::Strand::Top};let mut next=raw.get_untracked();pair_mut(&mut next,active_pair.get_untracked()).template_strand=strand;on_change.run(next);}><option value="top">"Top source strand"</option><option value="bottom" disabled=move||!allow_bottom.get()>"Saved bottom source strand"</option></select></label></Show><Show when=move||active_pair.get()==2><p class="helper">"Outer PCR uses the newly synthesized extension product's top frame."</p></Show>
        <For each=move||vec![active_pair.get()] key=|pair|*pair children=move|pair_index|view!{<div class="primer-grid"><PrimerFields raw on_change context=context.get_value() pair_index reverse=false/><PrimerFields raw on_change context=context.get_value() pair_index reverse=true/></div>}/>
        <div class="primer-settings-grid">
            <label>"Selected terminal homology (bp)"<input type="number" min="1" aria-label=format!("{} selected extension homology",context.get_value()) prop:value=move||raw.get().overlap on:input=move|e|update.run((Setting::Overlap,event_target_value(&e)))/></label>
            <label>"Minimum terminal pairing (bp)"<input type="number" min="1" aria-label=format!("{} minimum extension pairing",context.get_value()) prop:value=move||raw.get().min_paired_bases on:input=move|e|update.run((Setting::Minimum,event_target_value(&e)))/></label>
        </div>
        <details><summary>"Terminal-pair search limits"</summary><div class="primer-settings-grid">
            <label>"Maximum comparisons"<input type="number" min="1" max="50000000" aria-label=format!("{} extension comparison limit",context.get_value()) prop:value=move||raw.get().max_comparisons on:input=move|e|update.run((Setting::Comparisons,event_target_value(&e)))/></label>
            <label>"Maximum retained pairings"<input type="number" min="1" max="1000" aria-label=format!("{} extension candidate limit",context.get_value()) prop:value=move||raw.get().max_candidates on:input=move|e|update.run((Setting::Candidates,event_target_value(&e)))/></label>
        </div><p class="helper">"Competing terminal pairings and capped searches remain incomplete evidence, including after outer PCR."</p></details>
        <details><summary>"Review all six ordered oligos"</summary><pre>{move||{
            let raw=raw.get();
            [&raw.left,&raw.right,&raw.outer].into_iter().zip(["Left source PCR","Right source PCR","Outer PCR"]).map(|(p,label)|format!("{label}\n{}: {}{}\n{}: {}{}",p.forward.name,p.forward.tail,p.forward.annealing,p.reverse.name,p.reverse.tail,p.reverse.annealing)).collect::<Vec<_>>().join("\n\n")
        }}</pre></details>
    </section>}
}

/// Target-free material uses the same literal six-oligo form as assembly preparation.
#[component]
pub fn ReusableOverlapExtensionMaterialView(
    source: Memo<Option<collagene_core::Molecule>>,
    source_info: Memo<Option<crate::source_identity::SourceInfo>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<collagene_core::workspace::Workspace, String>>,
    on_open: Callback<crate::digest_product::MaterialOpenRequest, Result<(), String>>,
) -> impl IntoView {
    let raw = RwSignal::new(OverlapExtensionPcrInput::default());
    let origin = RwSignal::new(None::<crate::run_origin_model::PinnedOrigin>);
    let settings_revision =
        Memo::new(move |_| serde_json::to_string(&raw.get()).unwrap_or_default());
    view! {<details class="physical-evidence"><summary>"Overlap-extension PCR material"</summary>
        <p class="helper">"Declare the two source PCRs and the outer PCR. Reconstruct their product from the selected source without choosing an assembly target. Review the source origin, actual product and assumptions before opening reusable material."</p>
        <OverlapExtensionOligoFields raw=Memo::new(move|_|raw.get()) on_change=Callback::new(move|next|raw.set(next)) allow_bottom=Memo::new(move|_|origin.get().is_some()) context="Reusable material".to_string()/>
        <crate::reusable_material_view::ReusableMaterialView source source_info
            settings=Callback::new(move|()|Ok(crate::reusable_material_view::Settings::OverlapExtension{preparation:raw.get_untracked().request()?}))
            settings_revision kind=crate::reusable_material_view::Kind::OverlapExtensionPcr origin epoch get_workspace on_open/>
    </details>}
}
