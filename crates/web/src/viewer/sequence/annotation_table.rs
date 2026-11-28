//! Source-bound batch annotation editing and separate local visibility controls.
use collagene_core::sequence::{
    annotations::batch::{ColorEdit, Request},
    editing,
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::collections::BTreeSet;

const PAGE_SIZE: usize = 100;

#[derive(Clone, PartialEq, Eq)]
struct SourceKey {
    context: String,
    revision: String,
}

#[derive(Clone, PartialEq, Eq)]
struct Row {
    index: usize,
    label: String,
    kind: String,
    location: String,
    review: bool,
    qualifiers: usize,
}

fn short_label(label: &str) -> String {
    let mut chars = label.chars();
    let preview = chars.by_ref().take(100).collect::<String>();
    if chars.next().is_some() {
        format!("{preview}…")
    } else {
        preview
    }
}

#[component]
pub fn AnnotationTable(
    molecule: Memo<Molecule>,
    context: Memo<String>,
    locked: Memo<bool>,
    selected_feature: RwSignal<Option<usize>>,
    hidden: Memo<BTreeSet<usize>>,
    on_visibility: Callback<(Vec<usize>, bool), Result<(), String>>,
    on_inspect: Callback<usize>,
    on_edit: Callback<Option<usize>>,
    on_apply: Callback<Request, Result<(), String>>,
) -> impl IntoView {
    let identity = Memo::new(move |_| SourceKey {
        context: context.get(),
        revision: molecule.with(Molecule::revision),
    });
    let checked = RwSignal::new((identity.get_untracked(), BTreeSet::<usize>::new()));
    let query = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let kind = RwSignal::new(String::new());
    let flip = RwSignal::new(false);
    let apply_color = RwSignal::new(false);
    let color = RwSignal::new("#4f86c6".to_string());
    let status = RwSignal::new(String::new());
    let selection = Memo::new(move |_| {
        let current = identity.get();
        checked.with(|(key, indices)| {
            if *key == current {
                indices.clone()
            } else {
                BTreeSet::new()
            }
        })
    });
    Effect::new(move |_| {
        let key = identity.get();
        if checked.with_untracked(|(old, _)| *old != key) {
            checked.set((key, BTreeSet::new()));
            page.set(0);
            status.set("Source changed; batch selection cleared.".into());
        }
    });
    let rows = Memo::new(move |_| {
        molecule.with(|m| {
            m.features
                .iter()
                .enumerate()
                .map(|(index, feature)| Row {
                    index,
                    label: feature.label().to_string(),
                    kind: feature.kind.clone(),
                    location: editing::location_text(&feature.location)
                        .unwrap_or_else(|e| format!("Unresolved location: {e}")),
                    review: feature.needs_review,
                    qualifiers: feature.qualifiers.len(),
                })
                .collect::<Vec<_>>()
        })
    });
    let matching = Memo::new(move |_| {
        let needle = query.get().trim().to_lowercase();
        rows.with(|rows| {
            rows.iter()
                .filter(|row| {
                    needle.is_empty()
                        || row.label.to_lowercase().contains(&needle)
                        || row.kind.to_lowercase().contains(&needle)
                })
                .cloned()
                .collect::<Vec<_>>()
        })
    });
    Effect::new(move |_| {
        let _ = query.get();
        page.set(0);
    });
    let page_rows = Memo::new(move |_| {
        matching.with(|rows| {
            rows.iter()
                .skip(page.get().saturating_mul(PAGE_SIZE))
                .take(PAGE_SIZE)
                .cloned()
                .collect::<Vec<_>>()
        })
    });
    let mark = Callback::new(move |(indices, value): (Vec<usize>, bool)| {
        let key = identity.get_untracked();
        checked.update(|(saved, current)| {
            if *saved != key {
                *saved = key;
                current.clear();
            }
            for index in indices {
                if value {
                    current.insert(index);
                } else {
                    current.remove(&index);
                }
            }
        });
    });
    let visibility = Callback::new(move |(indices, hide): (Vec<usize>, bool)| {
        if indices.is_empty() {
            status.set("Check one or more annotations first.".into());
            return;
        }
        status.set(match on_visibility.run((indices, hide)) {
            Ok(()) => {
                if hide {
                    "Checked annotations hidden in views.".into()
                } else {
                    "Checked annotations shown in views.".into()
                }
            }
            Err(error) => error,
        });
    });
    let edit_batch = Callback::new(move |reset_color: bool| {
        if locked.get_untracked() {
            status.set("This source is read-only.".into());
            return;
        }
        let indices = selection.get_untracked().into_iter().collect::<Vec<_>>();
        if indices.is_empty() {
            status.set("Check one or more annotations first.".into());
            return;
        }
        let kind = if reset_color {
            None
        } else {
            let value = kind.get_untracked().trim().to_string();
            (!value.is_empty()).then_some(value)
        };
        let flip_strand = !reset_color && flip.get_untracked();
        let color = if reset_color {
            ColorEdit::Clear
        } else if apply_color.get_untracked() {
            ColorEdit::Set(color.get_untracked())
        } else {
            ColorEdit::Keep
        };
        if kind.is_none() && !flip_strand && matches!(color, ColorEdit::Keep) {
            status.set("Choose a feature type, strand flip or color change.".into());
            return;
        }
        let request = Request {
            source_revision: identity.get_untracked().revision,
            indices,
            kind,
            flip_strand,
            color,
        };
        status.set(match on_apply.run(request) {
            Ok(()) => "Annotation changes applied together.".into(),
            Err(error) => error,
        });
    });
    view! {
        <details class="annotation-table" open>
            <summary>"Annotation table"</summary>
            <div class="annotation-table-controls">
                <label>"Find annotations"<input aria-label="Find annotations by label or type" placeholder="Label or feature type" prop:value=move||query.get() on:input=move|e|{query.set(event_target_value(&e));}/></label>
                <p class="helper">{move||format!("{} matching of {} annotations · {} checked", matching.with(Vec::len), rows.with(Vec::len), selection.with(BTreeSet::len))}</p>
                <button class="button secondary" disabled=move||{locked.get()} on:click=move|_|{on_edit.run(None);}>"New annotation"</button>
                <div class="method-row-actions">
                    <label><input type="checkbox" aria-label="Check all annotations on this page" prop:checked=move||{page_rows.with(|rows|!rows.is_empty() && selection.with(|set|rows.iter().all(|row|set.contains(&row.index))))} on:change=move|e|{mark.run((page_rows.with_untracked(|rows|rows.iter().map(|row|row.index).collect()),event_target_checked(&e)));}/>"Check page"</label>
                    <button class="text-tool" on:click=move|_|{mark.run((matching.with_untracked(|rows|rows.iter().map(|row|row.index).collect()),true));}>"Check all matching"</button>
                    <button class="text-tool" on:click=move|_|{checked.set((identity.get_untracked(),BTreeSet::new()));}>"Clear checked"</button>
                </div>
                <div class="method-row-actions">
                    <button class="text-tool" disabled=move||{selection.with(BTreeSet::is_empty)} on:click=move|_|{visibility.run((selection.get_untracked().into_iter().collect(),false));}>"Show checked"</button>
                    <button class="text-tool" disabled=move||{selection.with(BTreeSet::is_empty)} on:click=move|_|{visibility.run((selection.get_untracked().into_iter().collect(),true));}>"Hide checked"</button>
                </div>
                <p class="helper">"Checks stay with the original annotations across search and pages. Visibility only changes the views."</p>
            </div>
            <div class="annotation-table-scroll" style="overflow:auto;max-height:32rem">
                <table><thead><tr><th>"Check"</th><th>"Annotation"</th><th>"Location"</th><th>"View"</th><th>"Actions"</th></tr></thead><tbody>
                    {move||page_rows.get().into_iter().map(|row|{
                        let index=row.index;
                        let label=short_label(&row.label);
                        view!{<tr class:active=move||{selected_feature.get()==Some(index)}>
                            <td><input type="checkbox" aria-label=format!("Check annotation {}",index+1) prop:checked=move||{selection.with(|set|set.contains(&index))} on:change=move|e|{mark.run((vec![index],event_target_checked(&e)));}/></td>
                            <td><strong title=row.label>{label}</strong><small>{format!("{} · {} qualifiers{}",row.kind,row.qualifiers,if row.review{" · review required"}else{""})}</small></td>
                            <td><span class="annotation-table-location">{row.location}</span></td>
                            <td><label><input type="checkbox" aria-label=format!("Show annotation {}",index+1) prop:checked=move||{hidden.with(|set|!set.contains(&index))} on:change=move|e|{visibility.run((vec![index],!event_target_checked(&e)));event_target::<web_sys::HtmlInputElement>(&e).set_checked(hidden.with_untracked(|set|!set.contains(&index)));}/>"Show"</label></td>
                            <td><button class="text-tool" on:click=move|_|{on_inspect.run(index);}>"Inspect"</button><button class="text-tool" disabled=move||{locked.get()} on:click=move|_|{on_edit.run(Some(index));}>"Edit"</button></td>
                        </tr>}
                    }).collect_view()}
                </tbody></table>
                <Show when=move||{matching.with(Vec::is_empty)}><p class="helper">"No annotations match this search."</p></Show>
            </div>
            <div class="method-row-actions">
                <button class="text-tool" disabled=move||{page.get()==0} on:click=move|_|{page.update(|page|*page=page.saturating_sub(1));}>"Previous page"</button>
                <span>{move||matching.with(|rows|{let count=rows.len();if count==0{"0 annotations".into()}else{format!("{}–{} of {}",page.get()*PAGE_SIZE+1,((page.get()+1)*PAGE_SIZE).min(count),count)}})}</span>
                <button class="text-tool" disabled=move||{matching.with(|rows|(page.get()+1)*PAGE_SIZE>=rows.len())} on:click=move|_|{page.update(|page|*page+=1);}>"Next page"</button>
            </div>
            <fieldset disabled=move||{locked.get()}><legend>"Edit checked annotations"</legend>
                <label>"Feature type (optional)"<input aria-label="Batch annotation type" placeholder="Keep current types" maxlength="256" prop:value=move||kind.get() on:input=move|e|{kind.set(event_target_value(&e));}/></label>
                <label><input type="checkbox" aria-label="Flip checked annotation strands" prop:checked=move||flip.get() on:change=move|e|{flip.set(event_target_checked(&e));}/>"Flip strand"</label>
                <p class="helper">"Flip the whole biological feature, including compound part order. Each checked feature is flipped from its own current orientation."</p>
                <label><input type="checkbox" aria-label="Apply annotation color" prop:checked=move||apply_color.get() on:change=move|e|{apply_color.set(event_target_checked(&e));}/>"Apply color"</label>
                <label>"Color"<input type="color" aria-label="Batch annotation color" prop:value=move||color.get() on:input=move|e|{color.set(event_target_value(&e));}/></label>
                <button class="button primary" disabled=move||{selection.with(BTreeSet::is_empty)} on:click=move|_|{edit_batch.run(false);}>"Apply to checked"</button>
                <button class="text-tool" disabled=move||{selection.with(BTreeSet::is_empty)} on:click=move|_|{edit_batch.run(true);}>"Reset checked colors"</button>
            </fieldset>
            <p class="helper" role="status">{move||status.get()}</p>
        </details>
    }.into_any()
}
