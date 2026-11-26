//! Compact, explicit internal intent editing beside the ordinary source slots.
use crate::{
    assembly_draft_model::Draft, assembly_substitution_model as model,
    assembly_workspace::SelectionNavigation, source_identity::SourceInfo,
};
use leptos::prelude::*;

#[component]
pub fn SubstitutionEditor(
    draft: RwSignal<Draft>,
    active_source: Memo<SourceInfo>,
    selected_region: Callback<(), Option<(usize, usize)>>,
    on_navigate: Callback<SelectionNavigation>,
) -> impl IntoView {
    let position = RwSignal::new("1".to_owned());
    let replacement = RwSignal::new("A".to_owned());
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
        position.set("1".into());
        status.set(String::new());
        page.set(0);
    });
    let add = move |_| match model::declare(
        &draft.get_untracked(),
        draft.with_untracked(|d| d.active_slot),
        &position.get_untracked(),
        &replacement.get_untracked(),
    ) {
        Ok(next) => {
            draft.set(next);
            status.set("Substitution declared. Review the intended construct, then verify its preparation.".into());
        }
        Err(error) => status.set(error),
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
                return Err(
                    "Open this contribution’s pinned source before copying a position.".to_owned(),
                );
            }
            let (start, end) = selected_region
                .run(())
                .ok_or("Select exactly one base in the sequence pane.")?;
            if start.checked_add(1) != Some(end) {
                return Err("Select exactly one base in the sequence pane.".to_owned());
            }
            model::selection_offset(&d, d.active_slot, start).map(|offset| offset + 1)
        })();
        match result {
            Ok(value) => {
                position.set(value.to_string());
                status.set("Position copied from the exact source selection. Choose its intended replacement.".into());
            }
            Err(error) => status.set(error),
        }
    };
    view! {<details class="assembly-internal-edits"><summary>{move||draft.with(|d|format!("Internal substitutions ({})",d.substitutions.len()))}</summary>
        <p class="helper">"Declare changes in the active fragment before choosing primers. Positions are one-based within its original contribution, in construct orientation. Junction trimming cannot silently remove an edit."</p>
        <Show when=move||draft.with(|d|d.alternatives.is_some())><p>"Internal substitutions currently require an ordinary construct. Library rows retain their existing junction-edit workflow."</p></Show>
        <Show when=move||draft.with(|d|d.alternatives.is_none())>
        <div class="primer-grid"><label>"Contribution position"<input type="number" min="1" maxlength="20" aria-label="Internal substitution contribution position" prop:value=move||position.get() on:input=move|e|position.set(event_target_value(&e))/></label>
        <label>"Intended replacement"<select aria-label="Internal substitution replacement" prop:value=move||replacement.get() on:change=move|e|replacement.set(event_target_value(&e))><option value="A">"A"</option><option value="C">"C"</option><option value="G">"G"</option><option value="T">"T"</option></select></label></div>
        <p class="helper">{move||{
            let d=draft.get();
            position.get().parse::<usize>().ok().and_then(|p|p.checked_sub(1)).and_then(|offset|model::source_base(&d,d.active_slot,offset).ok())
                .map(|(original,expected)|format!("Original base: {expected} · source position {}. Expected letters are stored with the declaration.",original+1))
                .unwrap_or_else(||"Choose a position in the pinned source contribution.".into())
        }}</p>
        <div class="method-row-actions"><button class="button secondary" on:click=from_selection>"Use single-base sequence selection"</button><button class="button primary" on:click=add>"Declare substitution"</button></div>
        </Show>
        <p role="status">{move||status.get()}</p>
        <ul>{move||draft.with(|d|d.substitutions.iter().skip(page.get().saturating_mul(20)).take(20).cloned().collect::<Vec<_>>()).into_iter().map(|edit|{
            let remove=edit.clone();let navigate=edit.clone();
            let label=draft.with_untracked(|d|d.methods.role_ids(d.slots.len()).ok().and_then(|roles|roles.iter().position(|r|r.role==edit.role)).and_then(|i|d.slots.get(i)).map(|s|s.label.clone()).unwrap_or_else(||edit.role.clone()));
            view!{<li><span>{format!("{} · position {} · {} → {}",label,edit.offset+1,edit.expected,edit.replacement)}</span>
                <button class="text-tool" on:click=move|_|{
                    let d=draft.get_untracked();
                    let result=(|| {
                        let roles=d.methods.role_ids(d.slots.len())?;
                        let index=roles.iter().position(|r|r.role==navigate.role).ok_or("The source contribution is unavailable.")?;
                        let (source_position,_)=model::source_base(&d,index,navigate.offset)?;
                        let slot=&d.slots[index];let capture=slot.captured.as_ref().ok_or("The original source is unavailable.")?;
                        Ok::<_,String>(SelectionNavigation{source:capture.source.binding.clone(),start:source_position,end:source_position+1,reverse:slot.reverse})
                    })();
                    match result{Ok(nav)=>on_navigate.run(nav),Err(error)=>status.set(error)}
                }>"View source base"</button>
                <button class="text-tool" on:click=move|_|{draft.update(|d|{if let Some(index)=d.substitutions.iter().position(|e|e==&remove){d.substitutions.remove(index);}});page.set(0);status.set("Substitution removed. Review the intended construct again.".into());}>"Remove"</button>
            </li>}
        }).collect_view()}</ul>
        <Show when=move||draft.with(|d|d.substitutions.len()>20)><div class="method-row-actions"><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous edits"</button><button class="text-tool" disabled={move||draft.with(|d|page.get().saturating_add(1).saturating_mul(20)>=d.substitutions.len())} on:click=move|_|page.update(|p|*p=p.saturating_add(1))>"Next edits"</button></div></Show>
        <p class="helper">"Replacing a source region or reversing a fragment clears its internal substitutions. Reordering retains edits with their source. Primer settings never populate this declaration."</p>
    </details>}
}
