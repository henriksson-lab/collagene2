use super::model::{Candidate, Draft};
use crate::{
    cloning::methods::checked_cloning::{AssumptionSignals, SourceAssumptions},
    history::{run_origin_model::PinnedOrigin, run_origin_view::SourceOriginPicker},
};
use collagene_core::{cloning::screening as native, documents::workspace::Workspace};
use leptos::prelude::*;
use std::sync::Arc;

#[component]
pub(super) fn ComponentEditor(
    id: String,
    draft: RwSignal<Draft>,
    workspace: RwSignal<Option<Arc<Workspace>>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<Workspace, String>>,
) -> impl IntoView {
    let id = StoredValue::new(id);
    let source = Memo::new(move |_| {
        let current = draft.with(|d| {
            d.components
                .iter()
                .find(|c| c.id == id.get_value())
                .cloned()
        })?;
        workspace.with(|w| {
            w.as_ref()
                .and_then(|w| super::model::info(w, &current.binding.event_id))
        })
    });
    let origin = Memo::new(move |_| {
        draft.with(|d| {
            d.components
                .iter()
                .find(|c| c.id == id.get_value())
                .and_then(|c| c.origin.clone())
        })
    });
    let assumptions = AssumptionSignals::new();
    if let Some(current) = draft.with_untracked(|d| {
        d.components
            .iter()
            .find(|c| c.id == id.get_value())
            .cloned()
    }) {
        assumptions.host.set(current.generation.host);
        assumptions.dam.set(current.generation.dam);
        assumptions.dcm.set(current.generation.dcm);
        assumptions.cpg.set(current.generation.cpg);
    }
    Effect::new(move |_| {
        assumptions.track();
        draft.update(|d| {
            if let Some(component) = d.components.iter_mut().find(|c| c.id == id.get_value()) {
                component.generation = assumptions.current();
            }
        });
    });
    let select_origin = Callback::new(move |pin: Option<PinnedOrigin>| {
        let source = source
            .get_untracked()
            .ok_or("The exact retained source is unavailable.")?;
        if pin
            .as_ref()
            .is_some_and(|pin| pin.event() != &source.binding)
        {
            return Err("Saved material belongs to another exact source event.".into());
        }
        draft.update(|d| {
            if let Some(component) = d.components.iter_mut().find(|c| c.id == id.get_value()) {
                component.origin = pin;
                component.declare_generation = false;
            }
        });
        Ok(())
    });
    view! {
        <details class="biology-section">
            <summary>{move || draft.with(|d| d.components.iter().find(|c| c.id == id.get_value()).map(|c| c.label.clone()).unwrap_or_default())}</summary>
            <label>"Component name"
                <input maxlength="256"
                    prop:value=move || draft.with(|d| d.components.iter().find(|c| c.id == id.get_value()).map(|c| c.label.clone()).unwrap_or_default())
                    on:input=move |event| draft.update(|d| {
                        if let Some(component) = d.components.iter_mut().find(|c| c.id == id.get_value()) {
                            component.label = event_target_value(&event);
                        }
                    })/>
            </label>
            <p class="helper">{move || source.get().map(|s| format!("{} · {} bp · event {}", s.name, s.length, s.binding.event_id)).unwrap_or_else(|| "Exact source unavailable".into())}</p>
            <SourceOriginPicker source epoch get_workspace selected=origin on_select=select_origin/>
            <Show when=move || origin.get().is_none()>
                <label class="checkbox-label">
                    <input type="checkbox"
                        prop:checked=move || draft.with(|d| d.components.iter().find(|c| c.id == id.get_value()).is_some_and(|c| c.declare_generation))
                        on:change=move |event| draft.update(|d| {
                            if let Some(component) = d.components.iter_mut().find(|c| c.id == id.get_value()) {
                                component.declare_generation = event_target_checked(&event);
                            }
                        })/>
                    "Declare original source generation (otherwise modifications stay unknown)"
                </label>
                <Show when=move || draft.with(|d| {
                    d.components.iter().find(|c| c.id == id.get_value())
                        .is_some_and(|c| c.declare_generation)
                })>
                    <SourceAssumptions side="Screening component" signals=assumptions/>
                </Show>
            </Show>
            <button on:click=move |_| draft.update(|d| d.remove_component(&id.get_value()))>"Remove component"</button>
        </details>
    }
}

pub(super) fn role_key(role: &native::Role) -> &'static str {
    match role {
        native::Role::Intended => "intended",
        native::Role::Decoy(native::DecoyKind::EmptyVector) => "empty",
        native::Role::Decoy(native::DecoyKind::WrongOrientation) => "orientation",
        native::Role::Decoy(native::DecoyKind::TandemInsert) => "tandem",
        native::Role::Decoy(native::DecoyKind::ParentCarryover) => "parent",
        native::Role::Decoy(native::DecoyKind::Other) => "other",
    }
}
pub(super) fn parse_role(value: &str) -> native::Role {
    native::Role::Decoy(match value {
        "empty" => native::DecoyKind::EmptyVector,
        "orientation" => native::DecoyKind::WrongOrientation,
        "tandem" => native::DecoyKind::TandemInsert,
        "parent" => native::DecoyKind::ParentCarryover,
        _ => native::DecoyKind::Other,
    })
}
#[component]
pub(super) fn SubjectEditor(id: String, draft: RwSignal<Draft>) -> impl IntoView {
    let id = StoredValue::new(id);
    let intended = Memo::new(move |_| {
        draft.with(|d| {
            d.subjects
                .iter()
                .find(|s| s.id == id.get_value())
                .is_some_and(|s| s.role == native::Role::Intended)
        })
    });
    view! {
        <details class="biology-section">
            <summary>{move || draft.with(|d| d.subjects.iter().find(|s| s.id == id.get_value()).map(|s| s.label.clone()).unwrap_or_default())}</summary>
            <label>"Outcome name"
                <input maxlength="256"
                    prop:value=move || draft.with(|d| d.subjects.iter().find(|s| s.id == id.get_value()).map(|s| s.label.clone()).unwrap_or_default())
                    on:input=move |event| draft.update(|d| {
                        if let Some(subject) = d.subjects.iter_mut().find(|s| s.id == id.get_value()) { subject.label = event_target_value(&event); }
                    })/>
            </label>
            <Show when=move || !intended.get() fallback=|| view! { <p>"Intended outcome"</p> }>
                <label>"Alternative role"
                    <select prop:value=move || draft.with(|d| d.subjects.iter().find(|s| s.id == id.get_value()).map(|s| role_key(&s.role)).unwrap_or("other"))
                        on:change=move |event| draft.update(|d| {
                            if let Some(subject) = d.subjects.iter_mut().find(|s| s.id == id.get_value()) { subject.role = parse_role(&event_target_value(&event)); }
                        })>
                        <option value="empty">"Empty vector"</option>
                        <option value="orientation">"Wrong orientation"</option>
                        <option value="tandem">"Tandem insert"</option>
                        <option value="parent">"Parent carryover"</option>
                        <option value="other">"Other named alternative"</option>
                    </select>
                </label>
            </Show>
            <p class="helper">"Select every actual component in this outcome. Multiple selections represent a mixture; no concentrations or new DNA are inferred."</p>
            <For each=move || draft.with(|d| d.components.iter().map(|c| (c.id.clone(), c.label.clone())).collect::<Vec<_>>())
                key=|(id, _)| id.clone()
                children=move |(component_id, _label)| {
                    let component_id = StoredValue::new(component_id);
                    view! {
                        <label class="checkbox-label">
                            <input type="checkbox"
                                prop:checked=move || draft.with(|d| d.subjects.iter().find(|s| s.id == id.get_value()).is_some_and(|s| s.components.contains(&component_id.get_value())))
                                on:change=move |event| draft.update(|d| {
                                    if let Some(subject) = d.subjects.iter_mut().find(|s| s.id == id.get_value()) {
                                        let key = component_id.get_value();
                                        subject.components.retain(|c| c != &key);
                                        if event_target_checked(&event) { subject.components.push(key); }
                                    }
                                })/>
                            {move || draft.with(|d| {
                                d.components.iter().find(|c| c.id == component_id.get_value())
                                    .map(|c| c.label.clone()).unwrap_or_default()
                            })}
                        </label>
                    }
                }/>
            <Show when=move || !intended.get()>
                <button on:click=move |_| draft.update(|d| d.subjects.retain(|s| s.id != id.get_value()))>"Remove outcome"</button>
            </Show>
        </details>
    }
}

#[component]
pub(super) fn CandidateEditor(index: usize, draft: RwSignal<Draft>) -> impl IntoView {
    // Browsing the catalogue is not a scientific selection. Only the user's
    // chosen 1–4 IDs become a checked selection when preparing the request.
    let enzymes = StoredValue::new(super::model::enzyme_choices());
    view! {
        <details class="biology-section">
            <summary>{move || draft.with(|d| {
                d.candidates.get(index).map(|candidate| match candidate {
                    Candidate::Digest { label, .. } => label.clone(),
                    Candidate::Pcr(pcr) => pcr.name.clone(),
                }).unwrap_or_default()
            })}</summary>
            <label>"Test name"
                <input maxlength="256"
                    prop:value=move || draft.with(|d| {
                        d.candidates.get(index).map(|candidate| match candidate {
                            Candidate::Digest { label, .. } => label.clone(),
                            Candidate::Pcr(pcr) => pcr.name.clone(),
                        }).unwrap_or_default()
                    })
                    on:input=move |event| draft.update(|d| match d.candidates.get_mut(index) {
                        Some(Candidate::Digest { label, .. }) => *label = event_target_value(&event),
                        Some(Candidate::Pcr(p)) => p.name = event_target_value(&event), _ => (),
                    })/>
            </label>
            {move || draft.with(|d| match d.candidates.get(index) {
                Some(Candidate::Digest { .. }) => {
                    view! {
                        <p class="helper">"Choose 1–4 exact enzymes for this digestion hypothesis. Add another candidate to compare separate digests."</p>
                        {enzymes.with_value(|value| value.as_ref().err().map(|error| view! {
                            <p role="alert">{format!("Restriction catalogue unavailable: {error}")}</p>
                        }))}
                        {enzymes.get_value().unwrap_or_default().into_iter().map(move |(id, name)| {
                            let enzyme_id = StoredValue::new(id);
                            view! {
                                <label class="checkbox-label">
                                    <input type="checkbox"
                                        prop:checked=move || draft.with(|d| {
                                            matches!(d.candidates.get(index), Some(Candidate::Digest { enzymes, .. }) if enzymes.contains(&enzyme_id.get_value()))
                                        })
                                        disabled=move || draft.with(|d| {
                                            matches!(d.candidates.get(index), Some(Candidate::Digest { enzymes, .. }) if enzymes.len() >= 4 && !enzymes.contains(&enzyme_id.get_value()))
                                        })
                                        on:change=move |event| draft.update(|d| {
                                            if let Some(Candidate::Digest { enzymes, .. }) = d.candidates.get_mut(index) {
                                                let key = enzyme_id.get_value();
                                                enzymes.retain(|e| e != &key);
                                                if event_target_checked(&event) && enzymes.len() < 4 { enzymes.push(key); }
                                            }
                                        })/>
                                    {name}
                                </label>
                            }
                        }).collect_view()}
                    }.into_any()
                }
                Some(Candidate::Pcr(_)) => view! {
                    <p class="helper">"Enter ordered 5′→3′ oligos as annealing sequence plus optional declared 5′ tail. Both primer roles and every modeled exact binding geometry are considered."</p>
                    {[
                        ("Forward annealing", 0), ("Forward 5′ tail", 1),
                        ("Reverse annealing", 2), ("Reverse 5′ tail", 3),
                    ].into_iter().map(move |(label, field)| view! {
                        <label>{label}
                            <input maxlength="200"
                                prop:value=move || draft.with(|d| match d.candidates.get(index) {
                                    Some(Candidate::Pcr(p)) => match field {
                                        0 => p.forward.clone(),
                                        1 => p.forward_tail.clone(),
                                        2 => p.reverse.clone(),
                                        _ => p.reverse_tail.clone(),
                                    },
                                    _ => String::new(),
                                })
                                on:input=move |event| draft.update(|d| {
                                    if let Some(Candidate::Pcr(p)) = d.candidates.get_mut(index) {
                                        let value = event_target_value(&event);
                                        match field {
                                            0 => p.forward = value,
                                            1 => p.forward_tail = value,
                                            2 => p.reverse = value,
                                            _ => p.reverse_tail = value,
                                        }
                                    }
                                })/>
                        </label>
                    }).collect_view()}
                    <label class="checkbox-label">
                        <input type="checkbox"
                            prop:checked=move || draft.with(|d| matches!(d.candidates.get(index), Some(Candidate::Pcr(p)) if p.bottom))
                            on:change=move |event| draft.update(|d| {
                                if let Some(Candidate::Pcr(p)) = d.candidates.get_mut(index) { p.bottom = event_target_checked(&event); }
                            })/>
                        "Use Bottom template strand"
                    </label>
                }.into_any(),
                None => ().into_any(),
            })}
            <button on:click=move |_| draft.update(|d| { if index < d.candidates.len() { d.candidates.remove(index); } })>"Remove test"</button>
        </details>
    }
}
