//! Bounded event-DAG navigation and retained workspace collections.
#[path = "history_graph.rs"]
mod graph;
use collagene_core::history::HistorySummary;
use collagene_model::Molecule;
use graph::{GraphIndex, PAGE_SIZE, ROW_HEIGHT};
use leptos::prelude::*;
use std::collections::BTreeMap;
use wasm_bindgen::{JsCast, closure::Closure};

fn short(id: &str) -> String {
    let label: String = id.chars().take(12).collect();
    if id.chars().count() > 12 {
        format!("{label}…")
    } else {
        label
    }
}
fn label(text: &str, limit: usize) -> String {
    let shown = text.chars().take(limit).collect::<String>();
    if text.chars().count() > limit {
        format!("{shown}…")
    } else {
        shown
    }
}
fn focus_event(id: String) {
    let Some(window) = web_sys::window() else {
        return;
    };
    let callback = Closure::once_into_js(move || {
        if let Some(node) = web_sys::window()
            .and_then(|w| w.document())
            .and_then(|d| d.get_element_by_id(&format!("history-node-{id}")))
            .and_then(|e| e.dyn_into::<web_sys::HtmlElement>().ok())
        {
            let _ = node.focus();
        }
    });
    let _ = window.request_animation_frame(callback.unchecked_ref());
}
fn node_description(node: &graph::GraphNode, operation: &str) -> String {
    let mut text = format!("{} · event {}", operation, short(&node.event_id));
    if node.current {
        text.push_str(" · current event");
    } else if node.ancestor {
        text.push_str(" · ancestor of current event");
    } else {
        text.push_str(" · other retained branch or root");
    }
    if let Some(distance) = node.undo_distance {
        text.push_str(&format!(" · undo distance {distance}"));
    }
    if let Some(distance) = node.redo_distance {
        text.push_str(&format!(" · redo distance {distance}"));
    }
    if node.branch_point {
        text.push_str(" · branch point");
    }
    if node.same_content_events > 1 {
        text.push_str(&format!(
            " · same molecule revision in {} events",
            node.same_content_events
        ));
    }
    text
}

#[component]
pub fn HistoryPanel(
    molecule: Memo<Molecule>,
    events: Memo<Vec<HistorySummary>>,
    current_event: Memo<String>,
    on_restore_revision: Callback<String, Result<(), String>>,
    on_fork_revision: Callback<String, Result<(), String>>,
    candidates: Memo<Vec<(String, String, String, bool)>>,
    on_open_candidate: Callback<String>,
    on_remove_candidate: Callback<String>,
    collections: Memo<Vec<(usize, String, usize)>>,
    on_open_oligos: Callback<usize>,
    #[prop(optional, into)] navigation: Option<Memo<(Vec<String>, Vec<String>)>>,
    #[prop(optional, into)] external_events: Option<Memo<Vec<(String, String, String)>>>,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let selected = RwSignal::new(None::<String>);
    let event_query = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let effective_page = Memo::new(move |_| {
        page.get()
            .min(events.with(|items| items.len().saturating_sub(1) / PAGE_SIZE))
    });
    let index = Memo::new_owning(move |_| {
        (
            events.with(|items| GraphIndex::new(items, &current_event.get())),
            true,
        )
    });
    let chosen = Memo::new(move |_| {
        selected
            .get()
            .filter(|id| events.with(|items| items.iter().any(|e| e.event_id == *id)))
            .unwrap_or_else(|| current_event.get())
    });
    let graph_page = Memo::new_owning(move |_| {
        let (undo, redo) = navigation.map(|n| n.get()).unwrap_or_default();
        (
            index.with(|i| {
                i.as_ref()
                    .ok()
                    .map(|i| i.page(effective_page.get(), &undo, &redo))
            }),
            true,
        )
    });
    let tab_stop = Memo::new(move |_| {
        let selected = chosen.get();
        graph_page.with(|page| {
            page.as_ref().and_then(|page| {
                page.nodes
                    .iter()
                    .find(|n| n.event_id == selected)
                    .or_else(|| page.nodes.first())
                    .map(|n| n.event_id.clone())
            })
        })
    });
    let external = Memo::new_owning(move |_| {
        (
            external_events
                .map(|m| {
                    m.get()
                        .into_iter()
                        .take(graph::MAX_GRAPH_EVENTS)
                        .map(|(id, name, revision)| (id, (label(&name, 160), revision)))
                        .collect::<BTreeMap<_, _>>()
                })
                .unwrap_or_default(),
            true,
        )
    });
    let jump = Callback::new(move |id: String| {
        let position = index
            .with_untracked(|i| i.as_ref().ok().and_then(|i| i.display_position(&id)))
            .or_else(|| {
                events.with_untracked(|items| {
                    items
                        .iter()
                        .position(|e| e.event_id == id)
                        .map(|i| items.len() - 1 - i)
                })
            });
        if let Some(position) = position {
            page.set(position / PAGE_SIZE);
            selected.set(Some(id.clone()));
            focus_event(id);
        }
    });
    let find_event = move || {
        let id = event_query.get_untracked().trim().to_owned();
        if events.with_untracked(|items| items.iter().any(|e| e.event_id == id)) {
            jump.run(id);
            status.set(
                "Retained event selected for inspection; the active molecule has not changed."
                    .into(),
            );
        } else {
            status.set("No retained event with that exact identity in this document. Open its source document first.".into());
        }
    };
    let key_select = Callback::new(move |(id, key): (String, String)| {
        let target = index.with_untracked(|i| {
            let i = i.as_ref().ok()?;
            let at = i.display_position(&id)?;
            let count = events.with_untracked(Vec::len);
            let next = match key.as_str() {
                "ArrowUp" => at.saturating_sub(1),
                "ArrowDown" => (at + 1).min(count.saturating_sub(1)),
                "Home" => 0,
                "End" => count.saturating_sub(1),
                _ => return None,
            };
            i.event_at(next).map(str::to_owned)
        });
        if let Some(id) = target {
            jump.run(id);
        }
    });
    let selected_summary = Memo::new_owning(move |_| {
        (
            events.with(|items| items.iter().find(|e| e.event_id == chosen.get()).cloned()),
            true,
        )
    });
    view! {
        <style>".history-panel .history-graph{margin:16px 0;border-top:1px solid #dce5df;padding-top:12px}.history-panel .history-graph h4{font-size:13px;margin:0}.history-event-lookup{display:flex;gap:6px;align-items:flex-end;margin:10px 0}.history-event-lookup label{flex:1;min-width:0;font-size:11px;margin:0;gap:5px}.history-event-lookup input{min-width:0;font-size:11px;padding:7px}.history-event-lookup button{font-size:11px;padding:7px}.history-graph-heading{display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:6px}.history-panel .history-graph .helper{font-size:12px;line-height:1.6;margin:8px 0}.history-graph-scroll{max-height:340px;overflow:auto;border:1px solid #dce5df;border-radius:6px;background:#fcfdfb;overscroll-behavior:contain}.history-graph-grid{display:grid;align-items:start}.history-graph-grid>svg{display:block}.history-graph-buttons{padding-top:12px;padding-bottom:12px}.history-graph-event{display:flex;flex-direction:column;justify-content:center;gap:3px;height:46px;width:100%;border:0;background:transparent;text-align:left;padding:4px 8px;color:#456351;line-height:1.2}.history-graph-event>span{font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:100%}.history-graph-event>small{font-size:10px;color:#7b8b7d;overflow:hidden;white-space:nowrap;text-overflow:ellipsis;max-width:100%}.history-graph-event:hover,.history-graph-event.is-selected{background:#e9f3ec}.history-graph-event:focus-visible{outline:2px solid #298a73;outline-offset:-2px}.history-graph-dot{fill:#fff;stroke:#879c90;stroke-width:2;cursor:pointer}.history-graph-dot.is-ancestor{fill:#99bca8;stroke:#5a8f73}.history-graph-dot.is-current{fill:#267956;stroke:#267956}.history-graph-dot.is-selected{stroke:#123e2b;stroke-width:3}.history-graph-edge.is-continuation{stroke-dasharray:3 3}.history-panel .history-graph-inspection{border:1px solid #dce5df;border-radius:6px;padding:10px;margin-top:10px;background:#f8fbf7}.history-panel .history-graph-inspection p,.history-graph-parent-list{font-size:12px;line-height:1.6}.history-graph-parent-list{padding-left:16px}.history-graph-parent-list li{margin:8px 0}.history-graph-parent-list code{display:block;font-size:10px;overflow-wrap:anywhere}.history-graph-inspection .history-actions button{white-space:normal}.history-graph-inspection .history-actions{gap:6px}.history-panel .history-graph-inspection h4{overflow-wrap:anywhere}.history-panel{overflow-wrap:anywhere}.history-panel .history-card{border:1px solid #dce5df;border-radius:8px;padding:12px;margin:10px 0;background:#fff}.history-panel .history-card.is-current{border-color:#378568;background:#f3faf6}.history-panel .history-card h4{margin:0 0 8px;font-size:14px}.history-panel .history-card p{margin:6px 0;font-size:12px;line-height:1.6}.history-panel .history-actions{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}.history-panel .history-actions button{font-size:12px}.history-panel .history-identifiers{font-size:11px;line-height:1.6}.history-panel .history-identifiers code{white-space:normal;overflow-wrap:anywhere}.history-panel .history-current{font-size:11px;font-weight:700;color:#256e51}.history-panel .history-pagination{display:flex;align-items:center;flex-wrap:wrap;gap:8px;margin:12px 0}.history-panel .history-pagination span{font-size:12px}.history-panel .history-status{font-size:13px;line-height:1.6}.history-panel h3{margin-top:24px}"</style>
        <section class="history-panel" aria-label="Document history and retained collections">
            <h2>"History"</h2>
            <p>{move || molecule.with(|m| format!("{} · {} sequence symbols", m.name, m.sequence.len()))}</p>
            <details class="history-identifiers"><summary>"Current document identifiers"</summary><p>"Revision: "<code>{move || molecule.with(Molecule::revision)}</code></p><p>"Event: "<code>{move || current_event.get()}</code></p></details>
            <p class="helper">"History records document revisions and their parent events. Scientific verification is assessed separately for the exact revision and cloning plan; restoring an event does not establish a verification result."</p>
            <p class="history-status" role="status">{move || status.get()}</p>
            <h3>"Retained revision branches"</h3>
            <Show when=move || events.with(Vec::is_empty)><p class="helper">"No retained events for this document."</p></Show>
            <section class="history-graph" aria-label="Retained event branch graph">
                <div class="history-graph-heading"><h4>"Event branch graph"</h4><button class="text-tool" on:click=move |_|jump.run(current_event.get_untracked())>"Show current event"</button></div>
                <div class="history-event-lookup"><label>"Event ID"<input aria-label="History event ID" maxlength="80" placeholder="event-…" prop:value=move||event_query.get() on:input=move|e|event_query.set(event_target_value(&e)) on:keydown=move|e|if e.key()=="Enter"{e.prevent_default();find_event();}/></label><button class="button secondary" on:click=move |_|find_event()>"Find event"</button></div>
                <Show when=move||navigation.is_none()><p class="history-status">"Exact navigation stack metadata is unavailable; no Undo/Redo positions are inferred from ancestry."</p></Show>
                <p class="helper">"Arrows point from an event to its recorded parents. Green nodes belong to the current ancestry; other nodes are retained branches or roots. Double rings mark branch points. Dashed links continue outside this page or document."</p>
                <p class="helper">"Undo and Redo labels describe the actual navigation stacks, separately from ancestry. Event identity stays separate when two events contain the same molecule revision. Use Up/Down, Home or End on a graph button to inspect events."</p>
                {move || index.with(|i|i.as_ref().err().cloned()).map(|message|view!{<p class="history-status" role="status">{message}</p>})}
                {move || graph_page.get().map(|plot|{
                    let width=plot.width();let height=plot.height();let omitted=plot.omitted_edges;
                    let circles=plot.nodes.iter().map(|node|{
                        let id=node.event_id.clone();let click=id.clone();let chosen_id=id.clone();let x=12+node.lane*20;let y=12+node.row*ROW_HEIGHT+ROW_HEIGHT/2;
                        let class=if node.current{"history-graph-dot is-current"}else if node.ancestor{"history-graph-dot is-ancestor"}else{"history-graph-dot"};
                        let title=events.with_untracked(|items|node_description(node,&items[node.event_index].operation));
                        let branch=node.branch_point;
                        view!{<g><title>{title}</title>{branch.then(||view!{<circle cx=x cy=y r="10" fill="none" stroke="#9eafa5" stroke-width="1"/>})}<circle class=class class:is-selected=move||chosen.get()==chosen_id cx=x cy=y r="6" data-event-id=id on:click=move |_|jump.run(click.clone())/></g>}
                    }).collect_view();
                    let arrows=plot.edges.into_iter().map(|edge|{
                        let x=12+edge.from_lane*20;let y=12+edge.from_row*ROW_HEIGHT+ROW_HEIGHT/2;
                        let (path,continuation,description)=if let Some((row,lane))=edge.to{let px=12+lane*20;let py=12+row*ROW_HEIGHT+ROW_HEIGHT/2;let direction=if py>=y{1f64}else{-1.};let from=y as f64+direction*7.;let to=py as f64-direction*8.;let middle=(from+to)/2.;(format!("M{x} {from} C{x} {middle}, {px} {middle}, {px} {to}"),"visible",format!("Parent event {}",edge.parent))}else{let description=if edge.parent_in_document{format!("Parent event {} continues on another graph page",edge.parent)}else{external.with_untracked(|m|m.get(&edge.parent).map(|(name,revision)|format!("Parent event {} belongs to document {name}, revision {revision}",edge.parent))).unwrap_or_else(||format!("Parent event {} is outside this document; its display metadata is unavailable",edge.parent))};(format!("M{x} {} Q{x} {}, 4 {}",y+7,y+16,y+16),if edge.parent_in_document{"page"}else{"document"},description)};
                        view!{<path class="history-graph-edge" class:is-continuation=continuation!="visible" d=path fill="none" stroke="#7c9d8b" stroke-width="1.5" marker-end="url(#history-parent-arrow)" data-child=edge.child data-parent=edge.parent data-continuation=continuation><title>{description}</title></path>}
                    }).collect_view();
                    let buttons=plot.nodes.into_iter().map(|node|{
                        let summary=events.with_untracked(|items|items[node.event_index].clone());let description=node_description(&node,&summary.operation);let id=node.event_id.clone();let clicked=id.clone();let key_id=id.clone();let pressed=id.clone();let tab=id.clone();let selected_class=id.clone();
                        let mut badges=Vec::new();if node.current{badges.push("Current event".to_string());}else if node.ancestor{badges.push("Current ancestry".into());}else{badges.push("Other branch/root".into());}if let Some(distance)=node.undo_distance{badges.push(if distance==1{"Undo next".into()}else{format!("Undo {distance}")});}if let Some(distance)=node.redo_distance{badges.push(if distance==1{"Redo next".into()}else{format!("Redo {distance}")});}if node.branch_point{badges.push("Branch point".into());}if node.same_content_events>1{badges.push(format!("Same revision · {} events",node.same_content_events));}
                        view!{<button id=format!("history-node-{id}") class="history-graph-event" class:is-selected=move||chosen.get()==selected_class data-event-id=id data-revision=summary.revision data-current=node.current.to_string() data-undo-distance=node.undo_distance.map(|n|n.to_string()) data-redo-distance=node.redo_distance.map(|n|n.to_string()) aria-label=description.clone() title=description aria-pressed=move||(chosen.get()==pressed).to_string() tabindex=move||if tab_stop.get().as_deref()==Some(tab.as_str()){"0"}else{"-1"} on:click=move |_|selected.set(Some(clicked.clone())) on:keydown=move|e|{let key=e.key();if matches!(key.as_str(),"ArrowUp"|"ArrowDown"|"Home"|"End"){e.prevent_default();key_select.run((key_id.clone(),key));}}><span>{label(&summary.operation,34)}</span><small>{badges.join(" · ")}</small></button>}
                    }).collect_view();
                    view!{<div class="history-graph-scroll"><div class="history-graph-grid" style:grid-template-columns=format!("{width}px minmax(170px, 1fr)") style:min-width=format!("{}px",width+170)><svg width=width height=height viewBox=format!("0 0 {width} {height}") aria-hidden="true"><defs><marker id="history-parent-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#7c9d8b"/></marker></defs>{arrows}{circles}</svg><div class="history-graph-buttons">{buttons}</div></div></div><Show when=move||{omitted>0}><p class="history-status">{format!("{omitted} parent links exceed this page's 2,000-edge drawing limit. Every parent identifier remains in the event details.")}</p></Show>}
                })}
                {move || selected_summary.get().map(|summary|{
                    let id=summary.event_id.clone();let restore=id.clone();let fork=id.clone();let current=summary.event_id==current_event.get();let parents=summary.parents;let parent_count=parents.len();let revision=summary.revision.clone();let matching=candidates.with(|items|items.iter().filter(|(_,_,target,_)|*target==revision).map(|(_,name,_,_)|label(name,160)).collect::<Vec<_>>());let has_matching=!matching.is_empty();
                    view!{<article class="history-graph-inspection" aria-label="Selected history event" data-event-id=id><h4>{summary.operation}</h4><p>{format!("{} sequence symbols · {} recorded parent{}",summary.sequence_length,parent_count,if parent_count==1{""}else{"s"})}</p><details class="history-identifiers"><summary>"Selected event identifiers"</summary><p>"Event: "<code>{summary.event_id}</code></p><p>"Molecule revision: "<code>{summary.revision}</code></p></details><ul class="history-graph-parent-list">{parents.into_iter().map(|parent|{
                        let local=events.with_untracked(|items|items.iter().any(|e|e.event_id==parent));let target=parent.clone();let metadata=external.with_untracked(|m|m.get(&parent).cloned());let text=if local{"Parent in this document".to_string()}else if let Some((name,revision))=metadata{format!("Parent in document {name} · revision {}",short(&revision))}else{"Parent outside this document; display metadata unavailable".into()};
                        view!{<li><span>{text}</span><code>{parent}</code>{local.then(||view!{<button class="text-tool" on:click=move |_|jump.run(target.clone())>"Show parent event"</button>})}</li>}
                    }).collect_view()}</ul><Show when=move||has_matching><p class="helper">{format!("This molecule revision is a saved candidate target: {}. This is a content match, not evidence that this event produced or verified the candidate.",matching.join(", "))}</p></Show><div class="history-actions"><button class="button secondary" disabled=current on:click=move |_|{match on_restore_revision.run(restore.clone()){Ok(())=>status.set("Selected event restored. Other branches remain retained.".into()),Err(e)=>status.set(e)}}>"Restore selected event"</button><button class="button secondary" on:click=move |_|{match on_fork_revision.run(fork.clone()){Ok(())=>status.set("Selected event opened as a new document with its recorded source parent.".into()),Err(e)=>status.set(e)}}>"Open selected event as new document"</button></div></article>}
                })}
            </section>
            <div class="history-pagination"><button class="button secondary" disabled=move || effective_page.get()==0 on:click=move |_| page.set(effective_page.get_untracked().saturating_sub(1))>"Previous 100 events"</button><span>{move || events.with(|items| if items.is_empty() {"0 events".into()} else {format!("Events {}–{} of {} · parents below their descendants",effective_page.get()*100+1,((effective_page.get()+1)*100).min(items.len()),items.len())})}</span><button class="button secondary" disabled=move || {(effective_page.get()+1)*100>=events.with(Vec::len)} on:click=move |_| page.set(effective_page.get_untracked()+1)>"Next 100 events"</button></div>
            {move || {
                let head = current_event.get();
                let order=index.with(|i|i.as_ref().ok().map(|i|i.display_indices(effective_page.get()))).unwrap_or_else(||events.with(|items|(0..items.len()).rev().skip(effective_page.get()*PAGE_SIZE).take(PAGE_SIZE).collect()));
                let summaries=events.with(|items|order.into_iter().map(|i|items[i].clone()).collect::<Vec<_>>());
                summaries.into_iter().map(|event| {
                    let is_current = event.event_id == head;
                    let restore = event.event_id.clone();
                    let fork = event.event_id.clone();
                    let operation = event.operation;
                    view! {
                        <article class=if is_current {"history-card is-current"} else {"history-card"}>
                            <h4>{operation}</h4>
                            <Show when=move || is_current><span class="history-current">"Current event"</span></Show>
                            <p>{format!("{} sequence symbols · revision {}", event.sequence_length,short(&event.revision))}</p>
                            <p>{if event.parents.is_empty() {"Root event".into()} else {format!("Parent events: {}",event.parents.iter().map(|id|short(id)).collect::<Vec<_>>().join(", "))}}</p>
                            <details class="history-identifiers"><summary>"Full revision and branch identifiers"</summary><p>"Event: "<code>{event.event_id}</code></p><p>"Revision: "<code>{event.revision}</code></p><ul>{event.parents.into_iter().map(|parent|view!{<li>"Parent event: "<code>{parent}</code></li>}).collect_view()}</ul></details>
                            <div class="history-actions"><button class="button secondary" disabled=is_current on:click=move |_| {match on_restore_revision.run(restore.clone()){Ok(())=>status.set("Revision restored. The previous branch remains retained.".into()),Err(e)=>status.set(e)}}>"Restore revision"</button><button class="button secondary" on:click=move |_| {match on_fork_revision.run(fork.clone()){Ok(())=>status.set("Revision opened as a new document with its recorded source parent.".into()),Err(e)=>status.set(e)}}>"Open revision as new document"</button></div>
                        </article>
                    }
                }).collect_view()
            }}
            <h3>"Saved cloning candidates"</h3><p class="helper">"Saved candidates are proposals. Open one for review and verify against the exact active inputs and independently chosen target; history navigation does not establish scientific verification."</p>
            <Show when=move || candidates.with(Vec::is_empty)><p class="helper">"No saved cloning candidates."</p></Show>
            {move || candidates.get().into_iter().map(|(plan_revision,name,target_revision,active)| {
                let open = plan_revision.clone(); let remove = plan_revision.clone();
                view!{<article class="history-card"><h4>{name}</h4><Show when=move || active><span class="history-current">"Open for review"</span></Show><p>{format!("Plan {} · target {}",short(&plan_revision),short(&target_revision))}</p><details class="history-identifiers"><summary>"Full candidate identifiers"</summary><p>"Plan revision: "<code>{plan_revision}</code></p><p>"Target revision: "<code>{target_revision}</code></p></details><div class="history-actions"><button class="button secondary" on:click=move |_| on_open_candidate.run(open.clone())>"Open candidate for review"</button><button class="button secondary" on:click=move |_| on_remove_candidate.run(remove.clone())>"Remove saved candidate"</button></div></article>}
            }).collect_view()}
            <h3>"Retained oligo collections"</h3>
            <Show when=move || collections.with(Vec::is_empty)><p class="helper">"No retained oligo collections."</p></Show>
            {move || collections.get().into_iter().map(|(index,label,count)| view!{<article class="history-card"><h4>{label}</h4><p>{format!("{} source rows",count)}</p><button class="button secondary" on:click=move |_| on_open_oligos.run(index)>"Open oligo collection"</button></article>}).collect_view()}
        </section>
    }
}
