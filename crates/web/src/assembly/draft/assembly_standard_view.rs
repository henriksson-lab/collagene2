use super::{assembly_standard_model::StandardDraft, model::Draft};
use collagene_core::cloning::assembly::standards::{self, declaration::ProfileRef};
use leptos::prelude::*;

#[component]
pub(crate) fn AssemblyStandardControls(
    draft: RwSignal<Draft>,
    locked: Memo<bool>,
) -> impl IntoView {
    let catalogue = StoredValue::new(standards::catalogue().map_err(|error| error.to_string()));
    let choice = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let bind = Callback::new(move |()| {
        if locked.get_untracked() {
            return;
        }
        let result = draft.with_untracked(|current| {
            let profiles = catalogue.get_value()?;
            let selected = profiles
                .iter()
                .find(|p| p.profile_id == choice.get_untracked())
                .ok_or("Choose an installed scoped part profile.")?;
            let roles = current.methods.role_ids(current.slots.len())?;
            let standard = StandardDraft::capture(
                current,
                ProfileRef {
                    profile_id: selected.profile_id.clone(),
                    profile_identity: selected.profile_identity.clone(),
                },
            )?;
            Ok::<_, String>((standard, roles))
        });
        match result {
            Ok((standard, roles)) => {
                draft.update(|d| {
                    d.methods.source_roles = Some(roles);
                    d.methods.standard = Some(standard);
                });
                status.set("Exact source roles bound. Declare a slot and both oriented body boundaries for every part before method review.".into());
            }
            Err(error) => status.set(error),
        }
    });
    view! {
        <details class="assembly-source-details">
            <summary>"Declared part standard (advisory)"</summary>
            <p class="helper">
                "Review the whole ordered source-part set under a pinned format profile. "
                "This checks scoped grammar and sequence facts, not vector compliance, actual ends, coding function or experimental success."
            </p>
            <label>"Installed source-part profile"
                <select disabled=move || locked.get() prop:value=move || choice.get()
                    on:change=move |event| choice.set(event_target_value(&event))>
                    <option value="">"Choose a scoped profile…"</option>
                    {catalogue.with_value(|value| value.as_ref().ok().map(|profiles| profiles.iter().map(|p| view! {
                        <option value=p.profile_id.clone()>{p.name.clone()}</option>
                    }).collect_view()))}
                </select>
            </label>
            {catalogue.with_value(|value| value.as_ref().err().map(|error| view! { <p role="alert">{error.clone()}</p> }))}
            {move || catalogue.with_value(|value| value.as_ref().ok().and_then(|profiles| {
                profiles.iter().find(|p| p.profile_id == choice.get()).map(|profile| view! {
                    <p class="helper">{profile.scope.clone()}</p>
                    <details><summary>"Profile source and data identity"</summary>
                        <p>{format!("{} · profile {} · catalogue {} / {}", profile.profile_id, profile.profile_identity, profile.catalog_pin.edition, profile.catalog_pin.content_hash)}</p>
                        {profile.definition.limitations.iter().map(|text| view! { <p class="helper">{text.clone()}</p> }).collect_view()}
                    </details>
                })
            }))}
            <button disabled=move || locked.get() on:click=move |_| bind.run(())>
                "Bind profile to current source roles (clears prior slot/body declarations)"
            </button>
            <button disabled=move || locked.get() || draft.with(|d| d.methods.standard.is_none())
                on:click=move |_| draft.update(|d| d.methods.standard = None)>"Clear declared standard"</button>
            <p role="status">{move || status.get()}</p>
            {move || draft.with(|d| d.methods.standard.as_ref().map(|standard| {
                let state = standard.validate_sources(d).and_then(|_| standard.declaration(&d.methods)).map(|_| "Declaration ready for fresh ordered source-format review.".to_string()).unwrap_or_else(|error| error);
                view! { <p class="helper">{state}</p> }
            }))}
            <For each=move || draft.with(|d| {
                if d.methods.standard.is_none() { return vec![]; }
                d.methods.role_ids(d.slots.len()).map(|roles| roles.into_iter().map(|r| r.role).collect()).unwrap_or_default()
            }) key=|role| role.clone() children=move |role| {
                let role = StoredValue::new(role);
                view! {
                    <details class="assembly-source-details">
                        <summary>{move || draft.with(|d| {
                            let roles = d.methods.role_ids(d.slots.len()).unwrap_or_default();
                            roles.iter().position(|r| r.role == role.get_value()).and_then(|i| d.slots.get(i)).map(|slot| slot.label.clone()).unwrap_or_else(|| role.get_value())
                        })}</summary>
                        <p class="helper">{move || draft.with(|d| {
                            d.methods.standard.as_ref()
                                .and_then(|s| s.rows.iter().find(|r| r.role == role.get_value()))
                                .map(|r| format!("Exact event {} · input {} · declared interval {:?}", r.event.event_id, r.input_id, r.selector))
                                .unwrap_or_else(|| "This role is unbound; rebind or clear the declaration.".into())
                        })}</p>
                        <label>"Declared grammar slot"
                            <select disabled=move || locked.get()
                                prop:value=move || draft.with(|d| d.methods.standard.as_ref().and_then(|s| s.rows.iter().find(|r| r.role == role.get_value())).map(|r| r.slot.clone()).unwrap_or_default())
                                on:change=move |event| draft.update(|d| {
                                    if let Some(row) = d.methods.standard.as_mut().and_then(|s| s.rows.iter_mut().find(|r| r.role == role.get_value())) { row.slot = event_target_value(&event); }
                                })>
                                <option value="">"Declare this part's role…"</option>
                                {move || draft.with(|d| {
                                    let id = d.methods.standard.as_ref().map(|s| s.profile.profile_id.clone());
                                    catalogue.with_value(|value| value.as_ref().ok().and_then(|profiles| {
                                        profiles.iter().find(|p| Some(&p.profile_id) == id.as_ref()).map(|p| p.definition.grammar.slots.iter().map(|slot| view! {
                                            <option value=slot.id.clone()>{slot.name.clone()}</option>
                                        }).collect_view())
                                    }))
                                })}
                            </select>
                        </label>
                        {[("Body start (zero-based oriented-part offset)", true), ("Body end (exclusive oriented-part offset)", false)].into_iter().map(move |(label, start)| view! {
                            <label>{label}
                                <input inputmode="numeric" maxlength="20" disabled=move || locked.get()
                                    prop:value=move || draft.with(|d| {
                                        d.methods.standard.as_ref()
                                            .and_then(|s| s.rows.iter().find(|r| r.role == role.get_value()))
                                            .map(|r| if start { r.body_start.clone() } else { r.body_end.clone() })
                                            .unwrap_or_default()
                                    })
                                    on:input=move |event| draft.update(|d| {
                                        if let Some(row) = d.methods.standard.as_mut().and_then(|s| s.rows.iter_mut().find(|r| r.role == role.get_value())) {
                                            if start { row.body_start = event_target_value(&event); } else { row.body_end = event_target_value(&event); }
                                        }
                                    })/>
                            </label>
                        }).collect_view()}
                    </details>
                }
            }/>
            <p class="helper">
                "No slot or body boundaries are guessed from names, gene arrows or primers. "
                "PCR tails/clamps are neither inserted nor stripped. Reordering uses exact stable roles; changed sources, orientation or preparations must be explicitly rebound. "
                "Fresh format findings are shown with the method preview and retained in its protocol."
            </p>
        </details>
    }
}
