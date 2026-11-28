use leptos::prelude::*;

use web_sys::HtmlInputElement;
#[component]
pub(crate) fn RestrictionFields(
    side: &'static str,
    enzymes: RwSignal<String>,
    fragment: RwSignal<String>,
    treatment: RwSignal<String>,
    repair: RwSignal<String>,
    reverse: RwSignal<bool>,
) -> impl IntoView {
    view! {
        <h4>{format!("{side} substrate")}</h4>
        <label>{format!("{side} enzymes (comma separated)")}<input aria-label=format!("{side} enzymes") prop:value=move||enzymes.get() on:input=move|e|enzymes.set(event_target_value(&e))/></label>
        <label>{format!("{side} digest fragment")}<input aria-label=format!("{side} fragment number") type="number" min="1" prop:value=move||fragment.get() on:input=move|e|fragment.set(event_target_value(&e))/></label>
        <label>{format!("{side} end repair")}<select aria-label=format!("{side} end repair") prop:value=move||repair.get() on:change=move|e|repair.set(event_target_value(&e))><option value="none">"No end repair"</option><option value="blunt">"Blunt: fill 5′ overhangs and remove 3′ overhangs"</option><option value="fill">"Fill 5′ overhangs only"</option><option value="chew">"Remove 3′ overhangs only"</option></select></label>
        <label>{format!("{side} end treatment")}<select aria-label=format!("{side} end treatment") prop:value=move||treatment.get() on:change=move|e|treatment.set(event_target_value(&e))><option value="none">"No treatment"</option><option value="phosphorylate">"Phosphorylate 5′ termini"</option><option value="dephosphorylate">"Dephosphorylate 5′ termini"</option></select></label>
        <label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse {side} fragment") prop:checked=move||reverse.get() on:change=move|e|reverse.set(event_target::<HtmlInputElement>(&e).checked())/>{format!("Reverse {side} fragment")}</label>
    }
}
