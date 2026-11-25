//! Historical editable recovery data. Viewing it never restores scientific evidence.
use crate::app::download;
use collagene_core::{assembly_library_history as history, workspace::Workspace};
use leptos::prelude::*;
use std::sync::Arc;

struct HumanSummary {
    name: String,
    bin_labels: Vec<String>,
    constructs: usize,
}
struct Loaded {
    hash: String,
    epoch: u64,
    attachment_bytes: usize,
    record: history::Record,
    summary: Option<HumanSummary>,
}
const HISTORY_PAGE: usize = 10;
const IMPACT_PAGE: usize = 20;
const SOURCE_PAGE: usize = 40;
fn page_count(total: usize, limit: usize) -> usize {
    total.max(1).div_ceil(limit)
}

#[component]
fn Pager(page: RwSignal<usize>, total: usize, limit: usize, label: &'static str) -> impl IntoView {
    let count = page_count(total, limit);
    view! {<div class="evidence-pager"><p>{move||format!("{} {} · page {} of {}",total,label,page.get().min(count-1)+1,count)}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|value|*value=value.saturating_sub(1))>{format!("Previous {label}")}</button>
        <button class="text-tool" disabled=move||page.get()+1>=count on:click=move|_|page.update(|value|*value=(*value+1).min(count-1))>{format!("Next {label}")}</button>
    </div>}
}

#[component]
pub fn LibraryHistoryView(
    snapshot: Callback<(), Result<Workspace, String>>,
    epoch: RwSignal<u64>,
    busy: RwSignal<bool>,
    on_restore_proposal: Callback<String>,
) -> impl IntoView {
    let entries = RwSignal::new(Vec::<history::Summary>::new());
    let listed_epoch = RwSignal::new(None::<u64>);
    let selected = RwSignal::new(None::<Arc<Loaded>>);
    let list_page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let _ = epoch.get();
        entries.set(Vec::new());
        listed_epoch.set(None);
        selected.set(None);
        list_page.set(0);
    });
    let refresh = move |_| {
        if busy.get_untracked() {
            return;
        }
        let saved_epoch = epoch.get_untracked();
        let outcome = snapshot
            .run(())
            .and_then(|frame| history::list(&frame).map_err(|error| error.to_string()));
        if epoch.get_untracked() != saved_epoch {
            return;
        }
        match outcome {
            Ok(values) => {
                let total = values.len();
                entries.set(values);
                listed_epoch.set(Some(saved_epoch));
                selected.set(None);
                list_page.set(0);
                status.set(format!("{total} saved library edits. These records contain historical editable data, not fresh scientific results."));
            }
            Err(error) => {
                entries.set(Vec::new());
                listed_epoch.set(None);
                selected.set(None);
                status.set(error);
            }
        }
    };
    let open = Callback::new(move |hash: String| {
        let outcome = (|| -> Result<Arc<Loaded>, String> {
            if busy.get_untracked() || listed_epoch.get_untracked() != Some(epoch.get_untracked()) {
                return Err(
                    "Refresh library history after the workspace operation completes.".into(),
                );
            }
            if !entries
                .with_untracked(|values| values.iter().any(|entry| entry.attachment_hash == hash))
            {
                return Err("Choose a record from the current workspace history list.".into());
            }
            let saved_epoch = epoch.get_untracked();
            let frame = snapshot.run(())?;
            let record = history::load(&frame, &hash).map_err(|error| error.to_string())?;
            let attachment_bytes = frame
                .attachments
                .get(&hash)
                .ok_or("Recovery record is no longer attached.")?
                .bytes
                .len();
            if epoch.get_untracked() != saved_epoch {
                return Err("Workspace changed while opening recovery history.".into());
            }
            let summary = crate::assembly_draft_model::decode_version(&record.before_draft, 3)
                .ok()
                .and_then(|draft| {
                    draft.alternatives.as_ref().map(|library| HumanSummary {
                        name: draft.name.clone(),
                        bin_labels: library.bins.iter().map(|bin| bin.label.clone()).collect(),
                        constructs: library.rows.len(),
                    })
                });
            Ok(Arc::new(Loaded {
                hash,
                epoch: saved_epoch,
                attachment_bytes,
                record,
                summary,
            }))
        })();
        match outcome {
            Ok(value) => {
                selected.set(Some(value));
                status.set("Historical recovery selected. Preview restoration creates a separate proposal for explicit review.".into());
            }
            Err(error) => {
                selected.set(None);
                status.set(error);
            }
        }
    });
    let guard = Callback::new(move |saved: Arc<Loaded>| -> Result<(), String> {
        if busy.get_untracked()
            || epoch.get_untracked() != saved.epoch
            || !selected.with_untracked(|value| {
                value
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &saved))
            })
        {
            return Err(
                "The workspace or selected recovery changed; refresh and select it again.".into(),
            );
        }
        let frame = snapshot.run(())?;
        let asset = frame
            .attachments
            .get(&saved.hash)
            .ok_or("The recovery record is no longer attached.")?;
        if asset.kind != history::KIND
            || asset.media_type != history::MEDIA_TYPE
            || asset.bytes.len() != saved.attachment_bytes
            || blake3::hash(&asset.bytes).to_hex().as_str() != saved.hash
        {
            return Err("The selected recovery attachment changed.".into());
        }
        Ok(())
    });
    let before_download = Callback::new(move |saved: Arc<Loaded>| {
        let result = guard.run(saved.clone()).and_then(|()| {
            download(
                &format!(
                    "assembly-library-{}-{}.json",
                    saved.record.library_id, saved.record.before_hash
                ),
                "application/json",
                &saved.record.before_draft,
            )
        });
        status.set(match result {Ok(())=>"Exact historical working-draft JSON downloaded; this does not install or verify it.".into(),Err(error)=>error});
    });
    let preview = Callback::new(move |saved: Arc<Loaded>| {
        let result = guard.run(saved.clone()).and_then(|()| {
            if saved.summary.is_none() {
                Err("Raw historical data; full draft validation failed. Download the exact raw JSON to inspect it.".into())
            } else {
                Ok(())
            }
        });
        match result {
            Ok(()) => on_restore_proposal.run(saved.hash.clone()),
            Err(error) => status.set(error),
        }
    });
    view!{<details class="assembly-library-history"><summary>"Saved library edits and recovery"</summary>
        <p class="helper">"Historical editable data preserves removed sources and manual settings. It restores no verification result, candidate authority or experimental evidence. Preview a restoration before applying any change."</p>
        <button class="button secondary" disabled=move||busy.get() on:click=refresh>"Refresh saved library edits"</button>
        {move||{
            let values=entries.get();let total=values.len();let start=list_page.get().min(page_count(total,HISTORY_PAGE)-1)*HISTORY_PAGE;
            view!{<section>{values.into_iter().skip(start).take(HISTORY_PAGE).map(|entry|{
                let hash=entry.attachment_hash.clone();
                view!{<article><p>{format!("Saved edit · {} affected constructs",entry.impact_count)}</p>
                    <small>{format!("Record {}",entry.attachment_hash.chars().take(12).collect::<String>())}</small>
                    <button class="text-tool" disabled=move||busy.get()||listed_epoch.get()!=Some(epoch.get()) on:click=move|_|open.run(hash.clone())>"Inspect saved edit"</button>
                </article>}
            }).collect_view()}<Pager page=list_page total limit=HISTORY_PAGE label="saved edits"/></section>}
        }}
        {move||selected.get().map(|saved|view!{<HistoryDetails saved busy on_download=before_download on_preview=preview/>})}
        <p role="status">{move||status.get()}</p>
    </details>}.into_any()
}

#[component]
fn HistoryDetails(
    saved: Arc<Loaded>,
    busy: RwSignal<bool>,
    on_download: Callback<Arc<Loaded>>,
    on_preview: Callback<Arc<Loaded>>,
) -> impl IntoView {
    let impact_page = RwSignal::new(0usize);
    let source_page = RwSignal::new(0usize);
    let download = saved.clone();
    let preview = saved.clone();
    let impact_total = saved.record.impacts.len();
    let source_total = saved.record.sources.len();
    let valid_summary = saved.summary.is_some();
    let human = saved.summary.as_ref().map(|summary|view!{
        <h3>{summary.name.clone()}</h3>
        <p>{format!("Saved before this edit · {} constructs · ordered source bins: {}",summary.constructs,summary.bin_labels.join(" → "))}</p>
    });
    let metadata = view! {<section>
        {human}
        <Show when=move||!valid_summary><p>"Raw historical data; full draft validation failed"</p></Show>
        <p>{format!("{} bytes of exact working-draft JSON · {} affected constructs · {} retained source bindings",saved.record.before_draft.len(),impact_total,source_total)}</p>
        <button class="text-tool" disabled=move||busy.get() on:click=move|_|on_download.run(download.clone())>"Download exact before-draft JSON"</button>
        <button class="button secondary" disabled=move||busy.get()||!valid_summary on:click=move|_|on_preview.run(preview.clone())>"Preview restoration"</button>
        <details><summary>"Exact recovery metadata"</summary>
            <p>{format!("Library {} · record schema {} · record bytes {} · exact raw working-draft bytes {}",saved.record.library_id,saved.record.schema_version,saved.attachment_bytes,saved.record.before_draft.len())}</p>
            <p>"Immutable record hash"</p><code>{saved.hash.clone()}</code>
            <p>"Before persisted active-draft hash"</p><code>{saved.record.before_active_hash.clone()}</code>
            <p>"Before exact raw working-draft hash"</p><code>{saved.record.before_hash.clone()}</code>
            <p>"After active-draft hash"</p><code>{saved.record.after_hash.clone()}</code>
            <p>{format!("After attachment kind {} · allocation counters {} → {} · retained high-water mark {}",saved.record.after_kind,saved.record.before_next_id,saved.record.after_next_id,saved.record.allocation_highwater)}</p>
            <p>"Reviewed proposal revision"</p><code>{saved.record.proposal_revision.clone()}</code>
            <p>"Recorded source workspace fingerprint"</p><code>{saved.record.source_workspace_fingerprint.clone()}</code>
        </details>
    </section>};
    let saved = StoredValue::new(saved);
    view!{<section class="library-history-record">{metadata}
        <details><summary>"Recorded construct impacts"</summary>
            {move||saved.with_value(|saved|saved.record.impacts.iter().skip(impact_page.get()*IMPACT_PAGE).take(IMPACT_PAGE).map(|impact|view!{
                <article><h4>{impact.name.clone()}</h4><p>{format!("Construct {} · removed {} · added {}",impact.construct_id,impact.construct_removed,impact.construct_added)}</p>
                    <p>{format!("Primer review bins: {}",if impact.primer_review_bins.is_empty(){"none".into()}else{impact.primer_review_bins.join(", ")})}</p>
                    <p>{format!("Removed preparation bins: {}",if impact.removed_preparation_bins.is_empty(){"none".into()}else{impact.removed_preparation_bins.join(", ")})}</p>
                    <p>{format!("Junction declarations retained {} · prior declarations retained in recovery {} · new junctions (one-based): {}",impact.retained_junctions,impact.recovered_junctions,if impact.new_junctions.is_empty(){"none".into()}else{impact.new_junctions.iter().map(|index|(index+1).to_string()).collect::<Vec<_>>().join(", ")})}</p>
                    <p>{format!("Raw overlap profile changed {} · role order changed {} · prior candidate intent removed/replaced {}",impact.overlap_text_changed,impact.role_order_changed,impact.candidate_cleared)}</p>
                </article>
            }).collect_view())}
            <Pager page=impact_page total=impact_total limit=IMPACT_PAGE label="construct impacts"/>
        </details>
        <details><summary>"Recorded retained source bindings"</summary>
            <p class="helper">"These are historical document/event/revision declarations. Restoration separately checks every captured source in the raw draft."</p>
            {move||saved.with_value(|saved|saved.record.sources.iter().skip(source_page.get()*SOURCE_PAGE).take(SOURCE_PAGE).map(|source|view!{
                <article><p>"Document"</p><code>{source.document_id.clone()}</code><p>"Event"</p><code>{source.event_id.clone()}</code><p>"Exact revision"</p><code>{source.revision.clone()}</code></article>
            }).collect_view())}
            <Pager page=source_page total=source_total limit=SOURCE_PAGE label="source bindings"/>
        </details>
    </section>}.into_any()
}
