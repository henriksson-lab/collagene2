//! Progressive alternative-source and selected-construct controls in the existing
//! dock. Every action remains a guarded edit of the same assembly draft.
use crate::{
    assembly_alternatives::{self as alternatives, Choice},
    assembly_draft_model::Draft,
};
use leptos::prelude::*;
use std::sync::Arc;

const COMBINATIONS_PER_PAGE: u128 = 20;
const ROWS_PER_PAGE: usize = 10;

enum Action {
    Enable,
    SelectBin(usize),
    CaptureChoice,
    Save,
    Add(Vec<String>),
    Open(String),
    RemoveChoice(String),
    RemoveRow(String),
}
type Apply = Callback<(Arc<Draft>, u64, Action)>;

fn clipped(text: &str, limit: usize) -> String {
    let mut chars = text.chars();
    let mut out: String = chars.by_ref().take(limit).collect();
    if chars.next().is_some() {
        out.push_str("… (abbreviated; inspect the source/preparation controls)");
    }
    out
}
fn choice_label(choice: &Choice) -> String {
    let Some(capture) = choice.slot.captured.as_ref() else {
        return "Unassigned source · use the existing source controls before preparation".into();
    };
    format!(
        "{} · {}",
        clipped(&capture.source.name, 120),
        crate::assembly_draft_model::selector_label(&choice.slot.selector)
    )
}
fn dispatch(saved: Arc<Draft>, captured_epoch: u64, apply: Apply) -> Callback<Action> {
    Callback::new(move |action| apply.run((saved.clone(), captured_epoch, action)))
}

#[component]
pub fn AssemblyAlternativesView(draft: RwSignal<Draft>, epoch: RwSignal<u64>) -> impl IntoView {
    // Captured molecules are shared Arc snapshots. Clone the editable draft once
    // per change, rather than copying source DNA for every displayed matrix row.
    let snapshot = Memo::new(move |_| Arc::new(draft.get()));
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let _ = epoch.get();
        status.set(String::new());
    });
    let apply = Callback::new(
        move |(saved, captured_epoch, action): (Arc<Draft>, u64, Action)| {
            let result = (|| -> Result<(Draft, &'static str), String> {
                if epoch.get_untracked() != captured_epoch
                    || !draft.with_untracked(|current| current == saved.as_ref())
                {
                    return Err("The assembly or workspace changed; use its current alternative controls again.".into());
                }
                let (next, message) = match action {
                    Action::Enable => {
                        let seed = format!(
                            "{}:{}:{}",
                            js_sys::Date::now(),
                            js_sys::Math::random(),
                            captured_epoch
                        );
                        let hash = blake3::hash(seed.as_bytes()).to_hex().to_string();
                        let id = format!("lib-{}", &hash[..12]);
                        (
                            alternatives::enable(&saved, &id)?,
                            "Source alternatives enabled. Existing slots became ordered bins and the current assembly became the first selected construct.",
                        )
                    }
                    Action::SelectBin(index) => {
                        if index >= saved.slots.len() {
                            return Err("The selected source bin is unavailable.".into());
                        }
                        let mut next = saved.as_ref().clone();
                        next.active_slot = index;
                        (
                            next,
                            "Active source bin selected. Use the existing source controls, then explicitly capture the edited source as a new alternative.",
                        )
                    }
                    Action::CaptureChoice => (
                        alternatives::capture_active_choice(&saved)?,
                        "Captured the active slot as a new immutable source alternative and selected it for this construct. Other source edits still need their own explicit capture.",
                    ),
                    Action::Save => (
                        alternatives::save_active(&saved)?,
                        "Current construct settings saved. Binding, preparation and verification results must still be computed explicitly.",
                    ),
                    Action::Add(choices) => (
                        alternatives::add_row(&saved, choices)?,
                        "Added and opened a distinct construct with the selected source combination. Review its name, primers, junctions and intended reference, then compute its preparation explicitly.",
                    ),
                    Action::Open(id) => (
                        alternatives::open_row(&saved, &id)?,
                        "Saved the current construct's editable settings and opened the selected construct. Previous scientific results are not restored.",
                    ),
                    Action::RemoveChoice(id) => (
                        alternatives::remove_choice(&saved, &id)?,
                        "Removed the unused source alternative. Retained construct identities and source choices were preserved.",
                    ),
                    Action::RemoveRow(id) => (
                        alternatives::remove_row(&saved, &id)?,
                        "Removed the selected construct row. Other constructs and source alternatives were preserved.",
                    ),
                };
                // SelectBin is the only direct local projection; model actions
                // also validate, and all installations retain the same boundary.
                next.validate()?;
                if epoch.get_untracked() != captured_epoch
                    || !draft.with_untracked(|current| current == saved.as_ref())
                {
                    return Err(
                        "The assembly changed during this edit; nothing was installed.".into(),
                    );
                }
                Ok((next, message))
            })();
            match result {
                Ok((next, message)) => {
                    draft.set(next);
                    status.set(message.into());
                }
                Err(error) => status.set(error),
            }
        },
    );
    view! {
        <details class="assembly-alternatives">
            <summary>"Alternative sources and selected constructs"</summary>
            <p class="helper">"Use ordered source bins for interchangeable backbones, promoters, inserts or other parts. Choose combinations explicitly, then edit and review each construct in this same assembly dock. These are design alternatives, separate from reaction-product enumeration."</p>
            {move||{
                let saved=snapshot.get();
                let captured_epoch=epoch.get();
                let action=dispatch(saved.clone(),captured_epoch,apply);
                if let Some(library)=&saved.alternatives {
                    let active=library.rows.iter().find(|row|row.id==library.active_row);
                    let pending=alternatives::require_active_choice_match(&saved).err();
                    view!{<section>
                        <p>{format!("Design {} · {} ordered bins · {} selected constructs · current {}",library.id,library.bins.len(),library.rows.len(),active.map(|row|row.id.as_str()).unwrap_or("unavailable"))}</p>
                        <p>{format!("Current working construct name: {}",clipped(&saved.name,256))}</p>
                        {pending.map(|message|view!{<p role="alert" class="helper">{message}</p>})}
                        <div class="method-row-actions">
                            <button class="button secondary" on:click=move|_|action.run(Action::Save)>"Save current construct settings"</button>
                            <button class="button secondary" on:click=move|_|action.run(Action::CaptureChoice)>"Capture active slot as new alternative"</button>
                        </div>
                        <p class="helper">"To add an alternative, select its bin below, use the existing source selection/interval/orientation and preparation controls, then capture it explicitly. Changes to primers and junctions remain per construct; selecting another construct saves those edits first. A changed source must be captured before saving, switching or running preparation."</p>
                    </section>}.into_any()
                } else {
                    view!{<section><p class="helper">"Set the ordinary backbone/insert slots and preparation first. Enabling alternatives keeps that current design as the first construct and turns each slot into an ordered source bin."</p>
                        <button class="button secondary" on:click=move|_|action.run(Action::Enable)>"Enable source alternatives for this assembly"</button>
                    </section>}.into_any()
                }
            }}
            <Show when=move||snapshot.with(|saved|saved.alternatives.is_some())>
                <AlternativeBins snapshot epoch apply/>
                <CombinationMatrix snapshot epoch apply/>
                <SelectedConstructs snapshot epoch apply/>
            </Show>
            <p class="helper" role="status">{move||status.get()}</p>
        </details>
    }.into_any()
}

#[component]
fn AlternativeBins(
    snapshot: Memo<Arc<Draft>>,
    epoch: RwSignal<u64>,
    apply: Apply,
) -> impl IntoView {
    view! {<details open><summary>"Ordered source bins"</summary>
        {move||{
            let saved=snapshot.get();let library=saved.alternatives.as_ref()?;
            let action=dispatch(saved.clone(),epoch.get(),apply);
            Some(library.bins.iter().enumerate().map(|(index,bin)|{
                let current=library.rows.iter().find(|row|row.id==library.active_row).and_then(|row|row.choices.get(index));
                view!{<article class="assembly-alternative-bin">
                    <h4>{format!("{} · {} · {} choices",index+1,bin.label,bin.choices.len())}</h4>
                    <details><summary>"Bin identity"</summary><code>{bin.id.clone()}</code></details>
                    <button class="text-tool" aria-pressed=(saved.active_slot==index).to_string() on:click=move|_|action.run(Action::SelectBin(index))>"Edit this bin in source controls"</button>
                    {bin.choices.iter().map(|choice|{
                        let id=choice.id.clone();
                        let used=library.rows.iter().filter(|row|row.choices.iter().any(|chosen|chosen==&choice.id)).count();
                        let selected=current==Some(&choice.id);
                        let reference=choice.slot.captured.as_ref().map(|capture|format!("Document {} · event {} · revision {}",capture.source.binding.document_id,capture.source.binding.event_id,capture.source.binding.revision));
                        view!{<details class="assembly-alternative-choice"><summary>{format!("{}{} · {} selected constructs",choice_label(choice),if selected{" · current choice"}else{""},used)}</summary>
                            <p><code>{choice.id.clone()}</code></p><p>{clipped(&choice.row.summary(),600)}</p>
                            {reference.map(|reference|view!{<p class="helper">{reference}</p>})}
                            <button class="text-tool" disabled={used>0||bin.choices.len()<=1} on:click=move|_|action.run(Action::RemoveChoice(id.clone()))>"Remove unused alternative"</button>
                        </details>}
                    }).collect_view()}
                </article>}
            }).collect_view())
        }}
    </details>}.into_any()
}

#[component]
fn CombinationMatrix(
    snapshot: Memo<Arc<Draft>>,
    epoch: RwSignal<u64>,
    apply: Apply,
) -> impl IntoView {
    let page = RwSignal::new(0u128);
    let jump = RwSignal::new("1".to_string());
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let _ = epoch.get();
        page.set(0);
        jump.set("1".into());
        status.set(String::new());
    });
    let size = Memo::new(move |_| {
        snapshot.with(|saved| {
            saved
                .alternatives
                .as_ref()
                .map(alternatives::matrix_size)
                .transpose()
        })
    });
    let move_page = Callback::new(move |next: u128| {
        let Some(total) = size.get_untracked().ok().flatten() else {
            return;
        };
        let last = total.saturating_sub(1) / COMBINATIONS_PER_PAGE;
        page.set(next.min(last));
        jump.set((page.get_untracked() + 1).to_string());
    });
    view! {<details><summary>"Choose construct combinations"</summary>
        <p class="helper">"Only twenty combinations are generated at a time. Add the combinations you want; the software does not build the full Cartesian library. Adding the same sources again creates a distinct construct for different manual design settings."</p>
        {move||{
            let saved=snapshot.get();let library=saved.alternatives.as_ref()?;
            let action=dispatch(saved.clone(),epoch.get(),apply);
            let result=(||->Result<_,String>{
                let total=size.get()?.ok_or("Source alternatives are unavailable.")?;
                let current=page.get().min(total.saturating_sub(1)/COMBINATIONS_PER_PAGE);
                let combinations=alternatives::combinations_page(library,current*COMBINATIONS_PER_PAGE,COMBINATIONS_PER_PAGE as usize)?;
                Ok((total,current,combinations))
            })();
            Some(match result {
                Err(error)=>view!{<p role="alert">{error}</p>}.into_any(),
                Ok((total,current,combinations))=>{
                    let pages=total.div_ceil(COMBINATIONS_PER_PAGE).max(1);
                    let rows=combinations.into_iter().enumerate().map(|(offset,choices)|{
                        let label=library.bins.iter().zip(&choices).map(|(bin,id)|{
                            let choice=bin.choices.iter().find(|choice|choice.id==*id);
                            format!("{}: {}",bin.label,choice.map(choice_label).unwrap_or_else(||"unavailable choice".into()))
                        }).collect::<Vec<_>>().join(" · ");
                        let selected=library.rows.iter().filter(|row|row.choices==choices).count();
                        view!{<article class="assembly-combination"><p>{format!("Combination {} · {label}",current*COMBINATIONS_PER_PAGE+offset as u128+1)}</p>
                            <p class="helper">{format!("{selected} selected constructs currently use these sources.")}</p>
                            <button class="button secondary" disabled={library.rows.len()>=alternatives::MAX_ROWS} on:click=move|_|action.run(Action::Add(choices.clone()))>"Add as new construct"</button>
                        </article>}
                    }).collect_view();
                    view!{<section><p>{format!("{total} possible source combinations · page {} of {pages}",current+1)}</p>{rows}
                        <div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|move_page.run(current.saturating_sub(1))>"Previous combinations"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|move_page.run(current+1)>"Next combinations"</button></div>
                    </section>}.into_any()
                }
            })
        }}
        <label>"Combination page (decimal)"<input inputmode="numeric" maxlength="39" prop:value=move||jump.get() on:input=move|event|jump.set(event_target_value(&event))/></label>
        <button class="text-tool" on:click=move|_|{let result=(||->Result<u128,String>{let value=jump.get_untracked().parse::<u128>().map_err(|_|"Enter a positive decimal page number.")?;let total=size.get_untracked()?.ok_or("Source alternatives are unavailable.")?;let pages=total.div_ceil(COMBINATIONS_PER_PAGE).max(1);if value==0||value>pages{return Err(format!("Enter a combination page from 1 to {pages}."));}Ok(value-1)})();match result{Ok(value)=>{move_page.run(value);status.set(String::new());},Err(error)=>status.set(error)}}>"Go to combination page"</button>
        <p class="helper" role="status">{move||status.get()}</p>
    </details>}.into_any()
}

#[component]
fn SelectedConstructs(
    snapshot: Memo<Arc<Draft>>,
    epoch: RwSignal<u64>,
    apply: Apply,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    Effect::new(move |_| {
        let _ = epoch.get();
        page.set(0);
    });
    view!{<details open><summary>"Selected constructs"</summary>
        <p class="helper">"Open a row to use the existing preparation, primer/junction editing, independent intended reference, product review and output controls. Only editable settings are retained here; scientific success is not restored."</p>
        {move||{
            let saved=snapshot.get();let library=saved.alternatives.as_ref()?;
            let action=dispatch(saved.clone(),epoch.get(),apply);
            let pages=library.rows.len().div_ceil(ROWS_PER_PAGE).max(1);let current=page.get().min(pages-1);
            let rows=library.rows.iter().skip(current*ROWS_PER_PAGE).take(ROWS_PER_PAGE).map(|row|{
                let active=row.id==library.active_row;let open_id=row.id.clone();let remove_id=row.id.clone();
                let choices=library.bins.iter().zip(&row.choices).map(|(bin,id)|{
                    let choice=bin.choices.iter().find(|choice|choice.id==*id);
                    format!("{}: {} ({id})",bin.label,choice.map(choice_label).unwrap_or_else(||"unavailable choice".into()))
                }).collect::<Vec<_>>().join(" · ");
                view!{<article class="assembly-selected-construct"><h4>{format!("{}{}",if active{&saved.name}else{&row.edits.name},if active{" · current working construct"}else{""})}</h4>
                    <p><code>{row.id.clone()}</code></p><p class="helper">{choices}</p>
                    <div class="method-row-actions"><button class="button secondary" aria-pressed=active.to_string() on:click=move|_|action.run(Action::Open(open_id.clone()))>"Save current and open this construct"</button><button class="text-tool" disabled={library.rows.len()<=1} on:click=move|_|action.run(Action::RemoveRow(remove_id.clone()))>"Remove construct row"</button></div>
                </article>}
            }).collect_view();
            Some(view!{<section>{rows}<p>{format!("{} selected constructs · page {} of {pages}",library.rows.len(),current+1)}</p><div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous constructs"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next constructs"</button></div></section>})
        }}
    </details>}.into_any()
}
