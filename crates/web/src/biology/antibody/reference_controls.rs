//! Compact, explicit reference scope; species is never inferred from a score.
use super::antibody_model;
use collagene_core::antibody_numbering::{
    references::Species,
    v_mapping::{Locus, ReferenceSelection},
};
use leptos::prelude::*;

#[component]
pub(super) fn ReferenceControls(
    references: RwSignal<ReferenceSelection>,
    status: RwSignal<String>,
) -> impl IntoView {
    let mouse = Memo::new(move |_| references.with(|r| r.species == Species::MouseC57Bl6));
    view! {
        <label>
            "Reference species / strain"
            <select
                prop:value=move || references.with(|r| format!("{:?}", r.species))
                on:change=move |event| {
                    match event_target_value(&event).as_str() {
                        "Human" => references.update(|r| r.species = Species::Human),
                        "MouseC57Bl6" => references.set(ReferenceSelection {
                            species: Species::MouseC57Bl6,
                            loci: vec![Locus::Heavy],
                        }),
                        _ => status.set("Unsupported antibody reference species or strain.".into()),
                    }
                }
            >
                <option value="Human">"Human"</option>
                <option value="MouseC57Bl6">"Mouse C57BL/6 · IGH only"</option>
            </select>
        </label>
        <Show when=move || !mouse.get()>
            <fieldset>
                <legend>"Selected human reference loci"</legend>
                {[("IGH · heavy", Locus::Heavy), ("IGK · kappa", Locus::Kappa), ("IGL · lambda", Locus::Lambda)]
                    .into_iter()
                    .map(|(label, locus)| {
                        let checked_locus = locus.clone();
                        view! {
                            <label class="check-row">
                                <input
                                    type="checkbox"
                                    prop:checked=move || references.with(|r| r.loci.contains(&checked_locus))
                                    on:change=move |event| {
                                        let checked = event_target_checked(&event);
                                        references.update(|r| {
                                            r.loci.retain(|value| value != &locus);
                                            if checked { r.loci.push(locus.clone()); }
                                        });
                                    }
                                />
                                {label}
                            </label>
                        }
                    })
                    .collect_view()}
            </fieldset>
        </Show>
        <p class="helper">{move || antibody_model::reference_label(&references.get())}</p>
        <Show when=move || mouse.get()>
            <p class="helper">
                "Mouse coverage is the pinned C57BL/6 IGH reference set only. Light chains, other strains and species, and non-IMGT numbering/regions are unsupported. Existing choices are not silently converted."
            </p>
        </Show>
        <p class="helper" role="alert">
            {move || references.get().validate().err().map(|e| e.to_string()).unwrap_or_default()}
        </p>
    }
}
