//! Explicit general edits in the active contribution, kept separate from preparation.
use crate::{
    assembly_draft_model::Draft, assembly_edit_model as model,
    assembly_workspace::SelectionNavigation, source_identity::SourceInfo,
};
use leptos::prelude::*;

#[component]
pub fn ContributionEditEditor(
    draft: RwSignal<Draft>,
    active_source: Memo<SourceInfo>,
    selected_region: Callback<(), Option<(usize, usize)>>,
    on_navigate: Callback<SelectionNavigation>,
) -> impl IntoView {
    let start = RwSignal::new("1".to_owned());
    let end = RwSignal::new("1".to_owned());
    let replacement = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let context = Memo::new(move |_| {
        draft.with(|d| {
            d.slots.get(d.active_slot).map(|s| {
                (
                    d.active_slot,
                    s.captured.as_ref().map(|c| c.source.binding.clone()),
                    s.start,
                    s.end,
                    s.reverse,
                )
            })
        })
    });
    Effect::new(move |_| {
        let _ = context.get();
        start.set("1".into());
        end.set("1".into());
        replacement.set(String::new());
        status.set(String::new());
        page.set(0);
    });
    let add = move |_| {
        let current = draft.get_untracked();
        match model::declare(
            &current,
            current.active_slot,
            &start.get_untracked(),
            &end.get_untracked(),
            &replacement.get_untracked(),
        ) {
            Ok(next) => {
                draft.set(next);
                replacement.set(String::new());
                status.set("Edit declared from the original reference. Review intended DNA, then verify the physical preparation.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let from_selection = move |_| {
        let result = (|| {
            let d = draft.get_untracked();
            let slot = d
                .slots
                .get(d.active_slot)
                .ok_or("Choose a source contribution.")?;
            let capture = slot.captured.as_ref().ok_or("Pin a source region first.")?;
            if capture.source.binding != active_source.get_untracked().binding {
                return Err("Open the exact pinned source before copying a selection.".to_owned());
            }
            let (a, b) = selected_region
                .run(())
                .ok_or("Select exact bases or place a caret in Sequence.")?;
            model::selection_interval(&d, d.active_slot, a, b)
        })();
        match result {
            Ok((a, b)) => {
                start.set((a + 1).to_string());
                end.set((b + 1).to_string());
                status.set("Exact source boundaries copied in construct orientation. Enter inserted DNA, or leave it empty to delete selected bases.".into());
            }
            Err(error) => status.set(error),
        }
    };
    view! {<details class="assembly-internal-edits"><summary>{move||draft.with(|d|format!("Internal insertions, deletions and replacements ({})",d.edits.len()))}</summary>
        <p class="helper">"Declare literal changes before preparing DNA. Boundaries are one-based within the original contribution in construct orientation: start is included, end is excluded. Equal boundaries insert DNA; empty inserted DNA deletes the interval. For example, start 3 and end 5 replace bases 3–4."</p>
        <Show when=move||draft.with(|d|d.alternatives.is_some())><p>"Internal edits currently require an ordinary construct."</p></Show>
        <Show when=move||draft.with(|d|d.alternatives.is_none())>
        <div class="primer-grid">
            <label>"Start boundary"<input type="number" min="1" maxlength="20" aria-label="Internal edit start boundary" prop:value=move||start.get() on:input=move|e|start.set(event_target_value(&e))/></label>
            <label>"End boundary (excluded)"<input type="number" min="1" maxlength="20" aria-label="Internal edit end boundary" prop:value=move||end.get() on:input=move|e|end.set(event_target_value(&e))/></label>
        </div>
        <label>"Inserted DNA (empty for deletion)"<textarea maxlength="100000" rows="2" spellcheck="false" aria-label="Internal edit inserted DNA" prop:value=move||replacement.get() on:input=move|e|replacement.set(event_target_value(&e))/></label>
        <p class="helper">{move||{
            let d=draft.get();
            let original=(|| {
                let a=start.get().parse::<usize>().ok()?.checked_sub(1)?;
                let b=end.get().parse::<usize>().ok()?.checked_sub(1)?;
                model::source_text(&d,d.active_slot,a,b).ok()
            })();
            match original {
                Some(text) if text.is_empty()=>"Insertion between original bases; no source letters removed.".to_owned(),
                Some(text)=>format!("Original reference: {}{} · {} bases. Exact expected letters are retained with the declaration.",text.chars().take(80).collect::<String>(),if text.len()>80{"…"}else{""},text.len()),
                None=>"Choose valid boundaries within the pinned contribution.".to_owned(),
            }
        }}</p>
        <div class="method-row-actions"><button class="button secondary" on:click=from_selection>"Use exact sequence selection"</button><button class="button primary" on:click=add>"Declare edit"</button></div>
        </Show>
        <p role="status">{move||status.get()}</p>
        <ul>{move||draft.with(|d|d.edits.iter().skip(page.get().saturating_mul(20)).take(20).cloned().collect::<Vec<_>>()).into_iter().map(|edit|{
            let remove=edit.clone();let navigate=edit.clone();
            let label=draft.with_untracked(|d|d.methods.role_ids(d.slots.len()).ok().and_then(|roles|roles.iter().position(|r|r.role==edit.role)).and_then(|i|d.slots.get(i)).map(|s|s.label.clone()).unwrap_or_else(||edit.role.clone()));
            let description=if edit.start==edit.end {format!("insert {} bases at boundary {}",edit.replacement.len(),edit.start+1)}else if edit.replacement.is_empty(){format!("delete bases {}–{}",edit.start+1,edit.end)}else{format!("replace bases {}–{} with {} bases",edit.start+1,edit.end,edit.replacement.len())};
            view!{<li><span>{format!("{label} · {description}")}</span>
                <button class="text-tool" on:click=move|_|{
                    let d=draft.get_untracked();
                    let result=(||{
                        let roles=d.methods.role_ids(d.slots.len())?;
                        let index=roles.iter().position(|r|r.role==navigate.role).ok_or("The source contribution is unavailable.")?;
                        let (start,end)=model::source_interval(&d,index,navigate.start,navigate.end)?;
                        let slot=&d.slots[index];
                        let capture=slot.captured.as_ref().ok_or("The pinned source is unavailable.")?;
                        Ok::<_,String>(SelectionNavigation{source:capture.source.binding.clone(),start,end,reverse:slot.reverse})
                    })();
                    match result {Ok(nav)=>on_navigate.run(nav),Err(error)=>status.set(error)}
                }>"View original interval"</button>
                <button class="text-tool" on:click=move|_|{
                    draft.update(|d|if let Some(i)=d.edits.iter().position(|e|e==&remove){d.edits.remove(i);});
                    page.set(0);
                    status.set("Edit removed. Review the intended construct again.".into());
                }>"Remove"</button>
            </li>}
        }).collect_view()}</ul>
        <Show when=move||draft.with(|d|d.edits.len()>20)><div class="method-row-actions"><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous edits"</button><button class="text-tool" disabled=move||draft.with(|d|page.get().saturating_add(1).saturating_mul(20)>=d.edits.len()) on:click=move|_|page.update(|p|*p=p.saturating_add(1))>"Next edits"</button></div></Show>
        <p class="helper">"Changing a source region or orientation clears its edits. Reordering retains their source roles. Primer proposals do not supply intended DNA; enter the physical preparation separately."</p>
    </details>}
}
