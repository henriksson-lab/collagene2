//! Historical decisions belong to the original checked construct. Loading a
//! decision never promotes it to fresh read evidence.
use collagene_core::{read_review, run_lineage::EventBinding, workspace::Workspace};
use leptos::prelude::*;
#[component]
pub fn ReadReviewHistory(
    reference: Memo<EventBinding>,
    snapshot: Callback<(), Result<Workspace, String>>,
    epoch: RwSignal<u64>,
    on_open: Callback<String, Result<(), String>>,
    on_reference: Callback<crate::trace_validation::ReferenceNavigation, Result<(), String>>,
) -> impl IntoView {
    let records = RwSignal::new(Vec::<read_review::Summary>::new());
    let page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let owner = reference.get();
        let _ = epoch.get();
        page.set(0);
        match snapshot
            .run(())
            .and_then(|frame| read_review::list(&frame).map_err(|e| e.to_string()))
        {
            Ok(all) => {
                records.set(
                    all.into_iter()
                        .filter(|r| r.reference.document_id == owner.document_id)
                        .collect(),
                );
                status.set(String::new());
            }
            Err(error) => {
                records.set(Vec::new());
                status.set(error);
            }
        }
    });
    view! {<details class="read-review-history"><summary>{move||format!("Recorded review decisions ({})",records.with(Vec::len))}</summary>
        <p class="helper">"These are historical caller choices. Unselected differences remain unreviewed, and a derived reference requires a new read check. Opening decisions does not restore fresh evidence."</p>
        {move||records.get().into_iter().skip(page.get()*20).take(20).map(|record|{
            let hash=record.attachment_hash.clone();
            view!{<article><p>{format!("{} decisions · {} proposed changes · reference event {}",record.decision_count,record.change_count,record.reference.event_id)}</p>
                <button class="button secondary" on:click=move|_|if let Err(error)=on_open.run(hash.clone()){status.set(error)}>"Open recorded decisions"</button>
                {record.target.map(|target|view!{<button class="text-tool" on:click=move|_|if let Err(error)=on_reference.run(crate::trace_validation::ReferenceNavigation{reference:target.clone(),selection:None,epoch:epoch.get_untracked()}){status.set(error)}>"Open derived reference event"</button>})}
            </article>}
        }).collect_view()}
        <p>{move||format!("Page {} of {}",page.get()+1,records.with(|r|r.len().max(1).div_ceil(20)))}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous recorded reviews"</button>
        <button class="text-tool" disabled=move||records.with(|r|(page.get()+1)*20>=r.len()) on:click=move|_|page.update(|p|*p+=1)>"Next recorded reviews"</button>
        <p role="status">{move||status.get()}</p>
    </details>}
}
