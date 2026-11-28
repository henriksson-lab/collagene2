//! Compact, explicit internal intent editing beside the ordinary source slots.
use crate::assembly::draft::assembly_workspace::SelectionNavigation;
use crate::assembly::draft::model::Draft;
use crate::assembly::edit::assembly_substitution_model as model;
use crate::sources::identity::SourceInfo;
use leptos::prelude::*;

#[component]
pub fn SubstitutionEditor(
    draft: RwSignal<Draft>,
    active_source: Memo<SourceInfo>,
    selected_region: Callback<(), Option<(usize, usize)>>,
    on_navigate: Callback<SelectionNavigation>,
    source_review: Memo<
        Option<std::sync::Arc<crate::assembly::edit::assembly_edit_model::SourceReview>>,
    >,
) -> impl IntoView {
    let position = RwSignal::new("1".to_owned());
    let replacement = RwSignal::new("A".to_owned());
    let status = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let context = Memo::new(move |_| {
        draft.with(|d| {
            d.slots.get(d.active_slot).map(|s| {
                (
                    crate::assembly::draft::assembly_alternatives::identity(d),
                    d.active_slot,
                    s.captured.as_ref().map(|c| c.source.binding.clone()),
                    s.selector.clone(),
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
        source_review.get_untracked().as_deref(),
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
            model::selection_offset(
                &d,
                d.active_slot,
                start,
                source_review.get_untracked().as_deref(),
            )
            .map(|offset| offset + 1)
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
        <Show when=move||draft.with(|d|d.alternatives.is_some())><p class="helper">"These edits apply only to the selected construct."</p></Show>
        <div>
        <div class="primer-grid"><label>"Contribution position"<input type="number" min="1" maxlength="20" aria-label="Internal substitution contribution position" prop:value=move||position.get() on:input=move|e|position.set(event_target_value(&e))/></label>
        <label>"Intended replacement"<select aria-label="Internal substitution replacement" prop:value=move||replacement.get() on:change=move|e|replacement.set(event_target_value(&e))><option value="A">"A"</option><option value="C">"C"</option><option value="G">"G"</option><option value="T">"T"</option></select></label></div>
        <p class="helper">{move||{
            let d=draft.get();
            position.get().parse::<usize>().ok().and_then(|p|p.checked_sub(1)).and_then(|offset|model::source_base(&d,d.active_slot,offset, source_review.get().as_deref()).ok())
                .map(|(original,expected)|format!("Original base: {expected} · {}. Expected letters are stored with the declaration.",match original{collagene_core::cloning::assembly::intent::ContributionPosition::Reference{position,..}=>format!("reference position {}",position+1),collagene_core::cloning::assembly::intent::ContributionPosition::Material{strand,offset}=>format!("actual {} strand position {}",match strand{collagene_core::chemistry::Strand::Top=>"top",collagene_core::chemistry::Strand::Bottom=>"bottom"},offset+1)}))
                .unwrap_or_else(||"Choose a position in the pinned source contribution.".into())
        }}</p>
        <div class="method-row-actions"><button class="button secondary" on:click=from_selection>"Use single-base sequence selection"</button><button class="button primary" on:click=add>"Declare substitution"</button></div>
        </div>
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
                        let regions=crate::assembly::edit::assembly_edit_model::top_regions(&d,index,navigate.offset,navigate.offset+1,source_review.get_untracked().as_deref())?;
                        let slot=&d.slots[index];let capture=slot.captured.as_ref().ok_or("The original source is unavailable.")?;
                        Ok::<_,String>(SelectionNavigation{source:capture.source.binding.clone(),regions,description:crate::assembly::draft::model::selector_label(&slot.selector)})
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
