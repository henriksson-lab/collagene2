//! Independent ligation/tailing intent edits; preparation settings are retained.
use crate::{assembly_draft_model::Draft, cloning_methods_model::Mode};
use collagene_core::assembly_junction_edit::Edit;
use leptos::prelude::*;

#[component]
pub fn LiteralIntentEditor(draft: RwSignal<Draft>) -> impl IntoView {
    let selected = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let count = move || {
        draft.with(|d| {
            d.slots.len().saturating_sub(usize::from(
                !d.circular && !matches!(d.methods.mode, Mode::Ligation | Mode::Tailing),
            ))
        })
    };
    let index = move || selected.get().min(count().saturating_sub(1));
    let terminal = move || draft.with(|d| !d.circular && index() + 1 == d.slots.len());
    let current = move || {
        draft
            .with(|d| {
                d.junction_edits
                    .iter()
                    .find(|e| e.junction == index())
                    .cloned()
            })
            .unwrap_or(Edit {
                junction: index(),
                remove_left: 0,
                remove_right: 0,
                inserted: String::new(),
            })
    };
    let change = Callback::new(move |(field, value): (u8, String)| {
        let mut next = draft.get_untracked();
        let i = selected.get_untracked().min(count().saturating_sub(1));
        if next.slots.is_empty()
            || i >= count()
            || (!matches!(next.methods.mode, Mode::Ligation | Mode::Tailing)
                && next.substitutions.is_empty()
                && next.edits.is_empty())
        {
            status.set("Choose a declared junction in the current construct.".into());
            return;
        }
        let mut edit = next
            .junction_edits
            .iter()
            .find(|e| e.junction == i)
            .cloned()
            .unwrap_or(Edit {
                junction: i,
                remove_left: 0,
                remove_right: 0,
                inserted: String::new(),
            });
        if field == 2 {
            if value.len() > 200 {
                status.set("Each literal insertion is limited to200 DNA characters.".into());
                return;
            }
            edit.inserted = value.to_ascii_uppercase();
        } else {
            let amount = if value.trim().is_empty() {
                0
            } else {
                match value.parse::<usize>() {
                    Ok(n) => n,
                    Err(_) => {
                        status.set("Enter a nonnegative whole number of bases.".into());
                        return;
                    }
                }
            };
            if field == 0 {
                edit.remove_left = amount;
            } else {
                edit.remove_right = amount;
            }
        }
        if !next.circular && i + 1 == next.slots.len() {
            edit.remove_right = 0;
        }
        next.junction_edits.retain(|e| e.junction != i);
        if edit.remove_left != 0 || edit.remove_right != 0 || !edit.inserted.is_empty() {
            next.junction_edits.push(edit);
            next.junction_edits.sort_by_key(|e| e.junction);
        }
        match next.validate() {
            Ok(()) => {
                draft.set(next);
                status.set("Intended DNA declaration updated. Review it again, then reconstruct the preparation and verify the physical product.".into());
            }
            Err(error) => status.set(error),
        }
    });
    view!{<details class="assembly-junction-options"><summary>"Declare intended junction DNA"</summary>
        <p class="helper">"Declare the desired reference DNA, then prepare fragments that produce it. These declarations retain your primers, digest choices, strand chemistry and selected annealing geometry. Verification compares the reconstructed physical product with this independently reviewed reference."</p>
        <Show when={move||count()>0} fallback=||view!{<p class="helper">"This single linear contribution has no internal junction."</p>}>
        <label>"Junction or final suffix"<select aria-label="Intended assembly junction" prop:value=move||index().to_string() on:change=move|e|if let Ok(i)=event_target_value(&e).parse(){selected.set(i);status.set(String::new());}>{move||draft.with(|d|d.slots.iter().enumerate().take(count()).map(|(i,slot)|{
            let label=if !d.circular&&i+1==d.slots.len(){format!("{} · final outer 3′ suffix",slot.label)}else{format!("{} → {}{}",slot.label,d.slots[(i+1)%d.slots.len()].label,if i+1==d.slots.len(){" · circle closure"}else{""})};
            view!{<option value=i.to_string()>{label}</option>}
        }).collect_view())}</select></label>
        <label>"Trim upstream intended contribution (bases)"<input type="number" min="0" aria-label="Intended upstream trim" prop:value=move||current().remove_left on:input=move|e|change.run((0,event_target_value(&e)))/></label>
        <Show when=move||!terminal()><label>"Trim downstream intended contribution (bases)"<input type="number" min="0" aria-label="Intended downstream trim" prop:value=move||current().remove_right on:input=move|e|change.run((1,event_target_value(&e)))/></label></Show>
        <label>{move||if terminal(){"Intended final top-strand suffix 5′→3′"}else{"Intended inserted DNA at this junction 5′→3′"}}<textarea maxlength="200" aria-label="Intended literal junction DNA" prop:value=move||current().inserted on:input=move|e|change.run((2,event_target_value(&e)))></textarea></label>
        </Show>
        <p class="helper">"Trims apply to the selected source contributions in construct orientation. Review intended construct checks that each region retains at least one source base. A final linear suffix has no downstream contribution to trim."</p>
        <Show when=move||!draft.with(|d|d.junction_edits.is_empty())><details><summary>"Review all declared edits"</summary>{move||draft.with(|d|d.junction_edits.iter().map(|e|view!{<p>{format!("Boundary {} · trim {} upstream / {} downstream · insert {}",e.junction+1,e.remove_left,e.remove_right,if e.inserted.is_empty(){"(none)"}else{&e.inserted})}</p>}).collect_view())}</details></Show>
        <p class="helper" role="status">{move||status.get()}</p>
    </details>}.into_any()
}
