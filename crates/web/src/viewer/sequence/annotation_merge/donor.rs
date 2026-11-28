//! Explicit retained donor events and a declared full-source coordinate transform.
use super::view::MergeReview;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    documents::{reference_edit::Operation, workspace::Workspace},
    material::lineage::EventBinding,
    sequence::annotations::merge::{Source, Transform},
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{collections::BTreeSet, sync::Arc};

#[derive(Clone)]
struct Choice {
    event: EventBinding,
    label: String,
    short_event: String,
}

#[component]
pub(crate) fn DonorAnnotations(
    molecule: Memo<Molecule>,
    info: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    on_commit: Callback<(EventBinding, Vec<Operation>), Result<(), String>>,
) -> impl IntoView {
    let loaded = RwSignal::new(None::<Arc<Workspace>>);
    let events = RwSignal::new(Vec::<Choice>::new());
    let event = RwSignal::new(String::new());
    let selected = RwSignal::new(BTreeSet::<usize>::new());
    let origin = RwSignal::new("0".to_string());
    let reverse = RwSignal::new(false);
    let event_page = RwSignal::new(0usize);
    let feature_page = RwSignal::new(0usize);
    let generation = RwSignal::new(0u32);
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let _ = (info.get(), epoch.get());
        loaded.set(None);
        events.set(Vec::new());
        event.set(String::new());
        selected.set(BTreeSet::new());
        event_page.set(0);
        feature_page.set(0);
    });
    Effect::new(move |_| {
        let _ = (event.get(), selected.get(), origin.get(), reverse.get());
        generation.update(|generation| *generation = generation.wrapping_add(1));
    });
    let load = Callback::new(move |()| {
        if locked.get_untracked() {
            return;
        }
        match get_workspace.run(()) {
            Ok(workspace) => {
                let mut choices: Vec<Choice> = workspace
                    .documents
                    .iter()
                    .filter_map(|document| {
                        Some(Choice {
                            event: EventBinding {
                                document_id: document.id.clone(),
                                event_id: document.head_event.clone()?,
                                revision: document.head.clone(),
                            },
                            short_event: String::new(),
                            label: workspace
                                .revisions
                                .get(&document.head)?
                                .name
                                .chars()
                                .take(96)
                                .collect(),
                        })
                    })
                    .collect();
                for index in 0..choices.len() {
                    let id = &choices[index].event.event_id;
                    let mut width = 14.min(id.len());
                    while width < id.len()
                        && choices.iter().any(|other| {
                            other.event.event_id != *id
                                && other.event.event_id.starts_with(&id[..width])
                        })
                    {
                        width += 1;
                    }
                    let short = id[..width].to_owned();
                    choices[index].short_event = short;
                }
                events.set(choices);
                loaded.set(Some(Arc::new(workspace)));
                event.set(String::new());
                selected.set(BTreeSet::new());
                generation.update(|generation| *generation = generation.wrapping_add(1));
                status.set("Choose an exact retained donor event and explicit annotations. Names do not establish correspondence.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let prepare = Callback::new(move |()| -> Result<(Source, Vec<String>), String> {
        if locked.get_untracked() {
            return Err("Wait for the workspace update.".into());
        }
        let selected_event = events
            .with_untracked(|events| {
                events
                    .iter()
                    .find(|choice| choice.event.event_id == event.get_untracked())
                    .map(|choice| choice.event.clone())
            })
            .ok_or("Choose a retained donor event.")?;
        let indices = selected.get_untracked().into_iter().collect::<Vec<_>>();
        let transform = Transform {
            origin: origin
                .get_untracked()
                .parse()
                .map_err(|_| "Enter a zero-based donor origin.")?,
            reverse: reverse.get_untracked(),
        };
        // Resolve again from the actual current workspace on this action. The
        // displayed retained snapshot never authorizes a missing or rebound event.
        let workspace = get_workspace.run(())?;
        let source = Source::from_donor(&workspace, &selected_event, transform, &indices)
            .map_err(|error| error.to_string())?;
        let ids = source.proposal_ids().map_err(|error| error.to_string())?;
        Ok((source, ids))
    });
    let input = Memo::new(move |_| {
        let mut settings = selected
            .get()
            .into_iter()
            .map(|index| index.to_string())
            .collect::<Vec<_>>();
        settings.extend([event.get(), origin.get(), reverse.get().to_string()]);
        (generation.get(), settings)
    });
    view! {
        <details>
            <summary>"Merge annotations from a retained reference"</summary>
            <p class="helper">
                "The complete donor DNA must match after the selected orientation and origin.
                Only annotations are merged; no DNA is replaced."
            </p>
            <button class="text-tool" disabled=move || locked.get()
                on:click=move |_| load.run(())>"Load retained reference choices"</button>
            <label>
                "Donor event"
                <select prop:value=move || event.get() on:change=move |change| {
                    event.set(event_target_value(&change));
                    selected.set(BTreeSet::new());
                    feature_page.set(0);
                }>
                    <option value="">"Choose a source"</option>
                    {move || events.get().into_iter().skip(event_page.get()*100).take(100).map(|choice| {
                        view! { <option value=choice.event.event_id.clone() title=choice.event.event_id.clone()>{format!("{} · {}", choice.label, choice.short_event)}</option> }
                    }).collect_view()}
                </select>
            </label>
            <button class="text-tool" disabled=move || event_page.get() == 0
                on:click=move |_| event_page.update(|page| *page=page.saturating_sub(1))>"Previous donors"</button>
            <button class="text-tool" disabled=move || event_page.get()+1 >= events.get().len().div_ceil(100)
                on:click=move |_| event_page.update(|page| *page+=1)>"Next donors"</button>
            <label>
                <input type="checkbox" prop:checked=move || reverse.get()
                    on:change=move |change| reverse.set(event_target_checked(&change))/>
                "Reverse complement donor before origin rotation"
            </label>
            <label>
                "Origin after optional reverse (zero-based)"
                <input type="number" min="0" prop:value=move || origin.get()
                    on:input=move |change| origin.set(event_target_value(&change))/>
            </label>
            {move || loaded.with(|workspace| workspace.as_ref().map(|workspace| {
                let selected_event = events.with(|events| events.iter()
                    .find(|choice| choice.event.event_id == event.get())
                    .map(|choice| choice.event.clone()));
                let donor = selected_event.as_ref().and_then(|event| workspace.revisions.get(&event.revision));
                let Some(donor) = donor else { return view! { <p>"Choose a donor to inspect annotations."</p> }.into_any(); };
                let pages = donor.features.len().div_ceil(100).max(1);
                let current = feature_page.get().min(pages-1);
                let rows = donor.features.iter().enumerate().skip(current*100).take(100).map(|(index, feature)| {
                    let label = format!("{} · {} · {}", index+1, feature.kind, feature.label().chars().take(96).collect::<String>());
                    view! {
                        <label class="checkbox-label">
                            <input type="checkbox" prop:checked=move || selected.with(|selected| selected.contains(&index))
                                on:change=move |change| selected.update(|selected| {
                                    if event_target_checked(&change) { selected.insert(index); } else { selected.remove(&index); }
                                })/>
                            {label}
                        </label>
                    }
                }).collect_view();
                view! {
                    {rows}
                    <button class="text-tool" disabled=current==0
                        on:click=move |_| feature_page.set(current.saturating_sub(1))>"Previous donor annotations"</button>
                    <span>{format!("Page {} of {pages}", current+1)}</span>
                    <button class="text-tool" disabled=current+1>=pages
                        on:click=move |_| feature_page.set((current+1).min(pages-1))>"Next donor annotations"</button>
                }.into_any()
            }))}
            <MergeReview molecule info epoch locked input prepare on_commit/>
            <p class="helper" role="status">{move || status.get()}</p>
        </details>
    }
}
