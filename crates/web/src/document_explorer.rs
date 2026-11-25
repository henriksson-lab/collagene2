//! Presentation-only document navigation. IDs, histories and scientific state stay with the app.
use collagene_core::{Alphabet, Topology};
use leptos::prelude::*;
use wasm_bindgen::{JsCast, closure::Closure};

#[derive(Clone, PartialEq, Eq)]
pub struct DocumentItem {
    pub id: String,
    pub name: String,
    pub length: usize,
    pub alphabet: Alphabet,
    pub topology: Topology,
    /// Actual recorded provider/accession or known imported format; None when unknown.
    pub source: Option<String>,
}
const PAGE: usize = 100;
fn index(
    documents: Memo<Vec<DocumentItem>>,
) -> Memo<std::collections::BTreeMap<String, DocumentItem>> {
    Memo::new(move |_| {
        documents.with(|items| {
            items
                .iter()
                .map(|item| (item.id.clone(), item.clone()))
                .collect()
        })
    })
}
fn item_title(item: &DocumentItem) -> String {
    format!(
        "{}\n{}\nDocument ID: {}",
        item.name,
        metadata(item),
        item.id
    )
}

fn kind(item: &DocumentItem) -> &'static str {
    match (item.alphabet, item.topology) {
        (Alphabet::Protein, _) => "Protein",
        (Alphabet::Rna, Topology::Circular) => "Circular RNA",
        (Alphabet::Rna, _) => "Linear RNA",
        (_, Topology::Circular) => "Circular DNA",
        _ => "Linear DNA",
    }
}
fn metadata(item: &DocumentItem) -> String {
    let base = format!(
        "{} {} · {}",
        item.length,
        item.alphabet.short_unit(),
        kind(item)
    );
    match item.source.as_deref().filter(|s| !s.trim().is_empty()) {
        Some(source) => format!("{base} · {source}"),
        None => base,
    }
}
fn icon(alphabet: Alphabet, topology: Topology) -> impl IntoView {
    let path = match (alphabet, topology) {
        (Alphabet::Protein, _) => "M2 13L5 3L8 13M3 9H7M9 13L12 3L15 13M10 9H14",
        (_, Topology::Circular) => "M8 2A6 6 0 1 1 2 8M2 8V3M2 8H7",
        _ => "M2 8H14M2 4V12M14 4V12",
    };
    view! {<svg class="explorer-document-icon" viewBox="0 0 16 16" aria-hidden="true"><path d=path/></svg>}
}
fn activate(id: String, current: Memo<String>, locked: Memo<bool>, callback: Callback<String>) {
    if !locked.get_untracked() && current.with_untracked(|active| active != &id) {
        callback.run(id);
    }
}
fn tab_id(id: &str) -> String {
    format!("document-tab-{}", blake3::hash(id.as_bytes()).to_hex())
}
fn show_tab(id: &str, focus: bool) {
    let Some(element) = web_sys::window()
        .and_then(|w| w.document())
        .and_then(|d| d.get_element_by_id(&tab_id(id)))
        .and_then(|e| e.dyn_into::<web_sys::HtmlElement>().ok())
    else {
        return;
    };
    if focus {
        let _ = element.focus();
    }
    // The button now lives inside a presentation shell. Scroll the actual
    // tablist, using viewport-relative rectangles rather than offset parents.
    if let Some(scroller) = element.closest(".document-tabs").ok().flatten() {
        let rect = element
            .parent_element()
            .map(|shell| shell.get_bounding_client_rect())
            .unwrap_or_else(|| element.get_bounding_client_rect());
        let viewport = scroller.get_bounding_client_rect();
        let at = f64::from(scroller.scroll_left());
        let width = f64::from(scroller.client_width());
        let left = rect.left() - viewport.left() - f64::from(scroller.client_left()) + at;
        let right = left + rect.width();
        if left < at {
            scroller.set_scroll_left(left.floor().max(0.) as i32);
        } else if right > at + width {
            scroller.set_scroll_left((right - width).ceil().max(0.) as i32);
        }
    }
}

#[component]
pub fn DocumentExplorer(
    documents: Memo<Vec<DocumentItem>>,
    active_document: Memo<String>,
    locked: Memo<bool>,
    on_activate: Callback<String>,
    on_open: Callback<()>,
) -> impl IntoView {
    let indexed = index(documents);
    let query = RwSignal::new(String::new());
    let sort = RwSignal::new("workspace".to_string());
    let page = RwSignal::new(0usize);
    let filtered = Memo::new(move |_| {
        let needle = query.get().trim().to_lowercase();
        let mut rows = documents.with(|items| {
            items
                .iter()
                .filter(|item| {
                    needle.is_empty()
                        || item.name.to_lowercase().contains(&needle)
                        || item.id.to_lowercase().contains(&needle)
                        || item
                            .source
                            .as_deref()
                            .is_some_and(|s| s.to_lowercase().contains(&needle))
                        || kind(item).to_lowercase().contains(&needle)
                })
                .cloned()
                .collect::<Vec<_>>()
        });
        match sort.get().as_str() {
            "name" => rows.sort_by_cached_key(|d| (d.name.to_lowercase(), d.id.clone())),
            "length" => rows.sort_by_cached_key(|d| {
                (
                    std::cmp::Reverse(d.length),
                    d.name.to_lowercase(),
                    d.id.clone(),
                )
            }),
            "source" => rows.sort_by_cached_key(|d| {
                (
                    d.source.is_none(),
                    d.source.as_deref().unwrap_or("").to_lowercase(),
                    d.name.to_lowercase(),
                    d.id.clone(),
                )
            }),
            _ => (),
        }
        rows
    });
    Effect::new(move |_| {
        let active = active_document.get();
        let rows = filtered.get();
        let target = rows
            .iter()
            .position(|d| d.id == active)
            .map(|i| i / PAGE)
            .unwrap_or(0);
        page.set(target);
    });
    let visible = Memo::new(move |_| {
        filtered.with(|rows| {
            rows.iter()
                .skip(page.get() * PAGE)
                .take(PAGE)
                .cloned()
                .collect::<Vec<_>>()
        })
    });
    view! {
        <aside class="document-explorer" aria-label="Document explorer">
            <div class="explorer-heading"><h2>"Documents"<span class="explorer-count">{move||documents.with(Vec::len)}</span></h2><button class="explorer-icon-button" aria-label="Open sequence file" title="Open sequence file" disabled=move||locked.get() on:click=move|_|{if !locked.get_untracked(){on_open.run(());}}>"+"</button></div>
            <div class="explorer-controls">
                <div class="explorer-search"><svg viewBox="0 0 16 16" aria-hidden="true"><circle cx="7" cy="7" r="4.5"/><path d="M10.5 10.5L14 14"/></svg><input type="search" aria-label="Filter documents" placeholder="Search documents" prop:value=move||query.get() on:input=move|e|query.set(event_target_value(&e))/></div>
                <label class="explorer-sort"><span>"Sort"</span><select aria-label="Sort documents" prop:value=move||sort.get() on:change=move|e|sort.set(event_target_value(&e))><option value="workspace">"Workspace order"</option><option value="name">"Name · A–Z"</option><option value="length">"Length · largest first"</option><option value="source">"Source"</option></select></label>
            </div>
            <div class="explorer-document-list" aria-label="Workspace documents">
                <For each=move||visible.get() key=|item|item.id.clone() children=move|item:DocumentItem|{
                    let selected=item.id.clone();let aria=item.id.clone();let target=item.id.clone();let row_id=item.id.clone();
                    let row=Memo::new(move|_|indexed.with(|items|items.get(&row_id).cloned()));
                    view!{<button class="document-entry" data-document-id=item.id.clone() data-document-name=move||row.with(|r|r.as_ref().map(|d|d.name.clone()).unwrap_or_default()) data-document-length=move||row.with(|r|r.as_ref().map(|d|d.length.to_string()).unwrap_or_default()) class:active=move||active_document.with(|id|id==&selected) aria-current=move||active_document.with(|id|if id==&aria{Some("page")}else{None}) title=move||row.with(|r|r.as_ref().map(item_title).unwrap_or_default()) disabled=move||locked.get() on:click=move|_|activate(target.clone(),active_document,locked,on_activate)>{move||row.get().map(|d|icon(d.alphabet,d.topology))}<span class="explorer-document-text"><strong>{move||row.with(|r|r.as_ref().map(|d|d.name.clone()).unwrap_or_default())}</strong><small>{move||row.with(|r|r.as_ref().map(metadata).unwrap_or_default())}</small></span></button>}
                }/>
                <Show when=move||filtered.with(Vec::is_empty)><p class="explorer-empty">{move||if documents.with(Vec::is_empty){"Open a sequence file to start."}else{"No matching documents."}}</p></Show>
            </div>
            <div class="explorer-list-footer"><span role="status">{move||filtered.with(|rows|if rows.is_empty(){"0 documents".into()}else{format!("{}–{} of {}",page.get()*PAGE+1,((page.get()+1)*PAGE).min(rows.len()),rows.len())})}</span><Show when=move||filtered.with(|rows|rows.len()>PAGE)><div class="explorer-pages"><button aria-label="Previous documents" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"‹"</button><button aria-label="Next documents" disabled=move||{(page.get()+1)*PAGE>=filtered.with(Vec::len)} on:click=move|_|page.update(|p|*p+=1)>"›"</button></div></Show></div>
            <button class="explorer-open" disabled=move||locked.get() on:click=move|_|{if !locked.get_untracked(){on_open.run(());}}>"+ Open another sequence"</button>
        </aside>
    }
}

#[component]
pub fn DocumentTabs(
    documents: Memo<Vec<DocumentItem>>,
    active_document: Memo<String>,
    locked: Memo<bool>,
    on_activate: Callback<String>,
    on_close: Callback<String>,
) -> impl IntoView {
    let indexed = index(documents);
    let alive = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, std::sync::atomic::Ordering::Relaxed));
    let closing_alive = alive.clone();
    Effect::new(move |_| {
        let id = active_document.get();
        let _ = documents.with(Vec::len);
        let alive = alive.clone();
        let frame = Closure::once_into_js(move || {
            if alive.load(std::sync::atomic::Ordering::Relaxed)
                && active_document.try_get_untracked().as_deref() == Some(id.as_str())
            {
                show_tab(&id, false);
            }
        });
        if let Some(window) = web_sys::window() {
            let _ = window.request_animation_frame(frame.unchecked_ref());
        }
    });
    let close_generation = RwSignal::new(0u32);
    let close_tab = Callback::new(move |id: String| {
        if locked.get_untracked()
            || !documents.with_untracked(|items| items.iter().any(|item| item.id == id))
        {
            return;
        }
        let was_active = active_document.with_untracked(|active| active == &id);
        on_close.run(id.clone());
        // Closing a background tab must not move focus out of the current task.
        if !was_active {
            return;
        }
        let expected = active_document.get_untracked();
        let expected_tabs = documents
            .with_untracked(|items| items.iter().map(|item| item.id.clone()).collect::<Vec<_>>());
        let ticket = close_generation.get_untracked().wrapping_add(1);
        close_generation.set(ticket);
        let alive = closing_alive.clone();
        let frame = Closure::once_into_js(move || {
            if !alive.load(std::sync::atomic::Ordering::Relaxed)
                || locked.try_get_untracked() != Some(false)
                || close_generation.try_get_untracked() != Some(ticket)
                || active_document.try_get_untracked().as_deref() != Some(expected.as_str())
                || documents.try_with_untracked(|items| {
                    items
                        .iter()
                        .map(|item| item.id.as_str())
                        .eq(expected_tabs.iter().map(String::as_str))
                        && !items.iter().any(|item| item.id == id)
                        && (expected.is_empty() || items.iter().any(|item| item.id == expected))
                }) != Some(true)
            {
                return;
            }
            if expected.is_empty() {
                if let Some(empty) = web_sys::window()
                    .and_then(|window| window.document())
                    .and_then(|document| document.get_element_by_id("no-open-documents"))
                    .and_then(|element| element.dyn_into::<web_sys::HtmlElement>().ok())
                {
                    let _ = empty.focus();
                }
            } else {
                show_tab(&expected, true);
            }
        });
        if let Some(window) = web_sys::window() {
            let _ = window.request_animation_frame(frame.unchecked_ref());
        }
    });
    view! {
        <div class="document-tabs" role="tablist" aria-label="Open sequence documents">
            <For each=move||documents.get() key=|item|item.id.clone() children=move|item:DocumentItem|{
                let selected=item.id.clone();let aria=item.id.clone();let focus=item.id.clone();let target=item.id.clone();let key_id=item.id.clone();
                let middle_id=item.id.clone();let close_id=item.id.clone();
                let row_id=item.id.clone();let row=Memo::new(move|_|indexed.with(|items|items.get(&row_id).cloned()));
                view!{<div class="document-tab-shell" role="presentation" class:active=move||active_document.with(|id|id==&selected) on:pointerdown=move|event:web_sys::PointerEvent|{if event.button()==1&&!locked.get_untracked(){event.prevent_default();}} on:auxclick=move|event:web_sys::MouseEvent|{
                    if event.button()==1&&!locked.get_untracked(){event.prevent_default();close_tab.run(middle_id.clone());}
                }>
                    <button class="document-tab" role="tab" aria-controls="active-sequence-document" id=tab_id(&item.id) data-document-id=item.id.clone() aria-selected=move||active_document.with(|id|id==&aria).to_string() tabindex=move||active_document.with(|id|if id==&focus{0}else{-1}) title=move||row.with(|r|r.as_ref().map(item_title).unwrap_or_default()) disabled=move||locked.get() on:click=move|_|activate(target.clone(),active_document,locked,on_activate) on:keydown=move|e:web_sys::KeyboardEvent|{
                        if locked.get_untracked()||e.alt_key()||e.ctrl_key()||e.meta_key(){return;}
                        let key=e.key();
                        if key=="Delete"{e.prevent_default();close_tab.run(key_id.clone());return;}
                        if !matches!(key.as_str(),"ArrowLeft"|"ArrowRight"|"Home"|"End"){return;}e.prevent_default();
                        let next=documents.with_untracked(|items|{if items.is_empty(){return None;}let at=items.iter().position(|d|d.id==key_id)?;let index=match key.as_str(){"Home"=>0,"End"=>items.len()-1,"ArrowLeft"=>(at+items.len()-1)%items.len(),_=>(at+1)%items.len()};Some(items[index].id.clone())});
                        if let Some(id)=next{activate(id.clone(),active_document,locked,on_activate);show_tab(&id,true);}
                    }>{move||row.get().map(|d|icon(d.alphabet,d.topology))}<span>{move||row.with(|r|r.as_ref().map(|d|d.name.clone()).unwrap_or_default())}</span></button>
                    <button class="document-tab-close" tabindex="-1" disabled=move||locked.get() aria-label=move||row.with(|r|format!("Close {} tab",r.as_ref().map(|item|item.name.as_str()).unwrap_or("document"))) title="Close tab · document remains in workspace" on:pointerdown=move|event:web_sys::PointerEvent|event.prevent_default() on:click=move|event|{event.stop_propagation();close_tab.run(close_id.clone());}>"×"</button>
                </div>}
            }/>
        </div>
    }
}
