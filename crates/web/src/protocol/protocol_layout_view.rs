//! Documentary plate assignment; no reaction approval or physical transfer.
use crate::platform::worker::CachedWorker;
use crate::protocol::protocol_layout_model as model;
use collagene_core::{
    documents::workspace::Workspace,
    output::layout::{self as native, Assignment, Plate, Well},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
#[component]
pub fn ProtocolLayoutView(
    batch_hash: String,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    on_close: Callback<()>,
) -> impl IntoView {
    let initial = (|| {
        let frame = snapshot.run(())?;
        let scope = model::Scope::capture(&frame, epoch.get_untracked(), batch_hash.clone())?;
        let choices = model::choices(&frame, &batch_hash)?;
        let layouts = model::saved_layouts(&frame, &batch_hash)?;
        Ok::<_, String>((scope, choices, layouts))
    })();
    let saved_layouts = RwSignal::new(
        initial
            .as_ref()
            .map(|(_, _, layouts)| layouts.clone())
            .unwrap_or_default(),
    );
    let assignment_page = RwSignal::new(0usize);
    let choices = StoredValue::new(
        initial
            .as_ref()
            .map(|(_, c, _)| c.clone())
            .unwrap_or_default(),
    );
    let scope = StoredValue::new(initial.as_ref().ok().map(|(s, _, _)| s.clone()));
    let status = RwSignal::new(initial.err().unwrap_or_default());
    let name = RwSignal::new("Construct plate layout".to_string());
    let plate = RwSignal::new(None::<Plate>);
    let fill_kind = RwSignal::new("");
    let replicates = RwSignal::new("1".to_string());
    let assignments = RwSignal::new(Vec::<Assignment>::new());
    let proposed = RwSignal::new(None::<Vec<Assignment>>);
    let selected = RwSignal::new(0usize);
    let source = RwSignal::new(0usize);
    let label = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let report = RwSignal::new(None::<Arc<native::Layout>>);
    let saved = RwSignal::new(None::<String>);
    let executor = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let nonce = StoredValue::new_local(Rc::new(Cell::new(0u32)));
    let transaction = StoredValue::new_local(Rc::new(RefCell::new(
        None::<crate::platform::storage::Transaction>,
    )));
    let lease = StoredValue::new_local(Rc::new(
        crate::platform::commit_lease::CommitLease::default(),
    ));
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = executor.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|t| t.borrow().clone()) {
            let _ = tx.abort();
        } else if lease.with_value(|l| l.cancel_unstarted()) {
            workspace_busy.set(false);
        }
    });
    let current_request = Callback::new(move |()| {
        model::request(
            name.get_untracked(),
            plate.get_untracked(),
            assignments.get_untracked(),
        )
    });
    let form_key = move || serde_json::to_vec(&current_request.run(())?).map_err(|e| e.to_string());
    Effect::new(move || {
        let _ = (name.get(), plate.get(), assignments.get());
        report.set(None);
        saved.set(None);
    });
    let run_alive = alive.clone();
    let run = Callback::new(move |capture: bool| {
        if pending.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let ready = (|| -> Result<_, String> {
            let frame = snapshot.run(())?;
            let scope = scope.get_value().ok_or("No valid saved batch context")?;
            scope.permit(&frame, epoch.get_untracked())?;
            let request = current_request.run(())?;
            let context = native::preview_revision(&frame, &request).map_err(|e| e.to_string())?;
            if capture
                && report
                    .get_untracked()
                    .as_ref()
                    .is_none_or(|r| r.context_revision != context)
            {
                return Err("Review these exact assignments before saving.".into());
            }
            let database = if capture {
                Some(
                    db.get_value()
                        .ok_or("Native recovery storage is unavailable")?,
                )
            } else {
                None
            };
            let key = form_key()?;
            let worker = Rc::new(CachedWorker::new()?);
            Ok((frame, scope, request, context, database, key, worker))
        })();
        let (frame, pinned, request, context, database, key, w) = match ready {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let token = nonce.with_value(|n| {
            let next = n.get().wrapping_add(1);
            n.set(next);
            next
        });
        if capture {
            if let Err(e) = lease.with_value(|l| l.acquire(token as u64)) {
                status.set(e);
                return;
            }
            workspace_busy.set(true);
            save_generation.update(|g| *g = g.wrapping_add(1));
        }
        pending.set(true);
        executor.set_value(Some(w.clone()));
        let alive = run_alive.clone();
        let pending_tx = transaction.get_value();
        let lock = lease.get_value();
        let counter = nonce.get_value();
        spawn_local(async move {
            let guard = || -> Result<(), String> {
                if !alive.load(Ordering::Relaxed) || counter.get() != token {
                    return Err("Layout operation cancelled.".into());
                }
                pinned.permit(&snapshot.run(())?, epoch.get_untracked())?;
                if form_key()? != key {
                    return Err("Layout assignments changed; review again.".into());
                }
                Ok(())
            };
            let result = async {
                let job = if capture {
                    AnalysisJob::CaptureProtocolLayout {
                        workspace: frame.clone(),
                        request,
                    }
                } else {
                    AnalysisJob::PreviewProtocolLayout {
                        workspace: frame.clone(),
                        request,
                    }
                };
                let reply = w
                    .ask(Request {
                        id: token,
                        revision: context.clone(),
                        job,
                    })
                    .await?;
                guard()?;
                if reply.id != token || reply.revision != context {
                    return Err("Layout response differs from its request.".into());
                }
                match reply.result {
                    AnalysisResult::ProtocolLayoutPreview(layout) if !capture => {
                        layout
                            .validate_binding(&context)
                            .map_err(|e| e.to_string())?;
                        report.set(Some(Arc::new(layout)));
                        Ok(())
                    }
                    AnalysisResult::ProtocolLayoutCaptured(commit) if capture => {
                        commit
                            .validate_binding(&context)
                            .map_err(|e| e.to_string())?;
                        if commit.source_workspace_fingerprint != pinned.fingerprint {
                            return Err("Layout capture source differs.".into());
                        }
                        let hash = commit.attachment_hash.clone();
                        let layout = Arc::new(commit.layout.clone());
                        let listed = model::saved_layouts(&commit.workspace, &pinned.batch_hash)?;
                        let mut rebound =
                            model::Scope::capture(&commit.workspace, 0, pinned.batch_hash.clone())?;
                        let next = crate::workspace::workspace_commit::persist(
                            database.as_ref().ok_or("Missing recovery storage")?,
                            commit.workspace,
                            frame,
                            validate,
                            w.clone(),
                            guard,
                            pending_tx,
                        )
                        .await?;
                        on_restore.run(next);
                        if alive.load(Ordering::Relaxed) {
                            rebound.epoch = epoch.get_untracked();
                            saved_layouts.set(listed);
                            scope.set_value(Some(rebound));
                            report.set(Some(layout));
                            saved.set(Some(hash));
                        }
                        Ok(())
                    }
                    AnalysisResult::Error(e) => Err(e),
                    _ => Err("Unexpected layout response.".into()),
                }
            }
            .await;
            w.close();
            if capture && lock.release(token as u64) {
                workspace_busy.set(false);
            }
            if alive.load(Ordering::Relaxed) && counter.get() == token {
                pending.set(false);
                executor.set_value(None);
                status.set(match result {
                    Ok(()) => {
                        if capture {
                            "Layout saved. Source documents are unchanged.".into()
                        } else {
                            "Assignments checked. Reaction amounts remain unset.".into()
                        }
                    }
                    Err(e) => e,
                });
            }
        });
    });
    let export_alive = alive.clone();
    let export = Callback::new(move |format: native::ExportFormat| {
        if pending.get_untracked() {
            return;
        }
        let result = (|| -> Result<_, String> {
            let hash = saved
                .get_untracked()
                .ok_or("Save the current layout before export")?;
            let frame = snapshot.run(())?;
            let pin = scope.get_value().ok_or("Missing layout scope")?;
            pin.permit(&frame, epoch.get_untracked())?;
            let context = native::export_revision(&frame, &hash, format.clone())
                .map_err(|e| e.to_string())?;
            let w = Rc::new(CachedWorker::new()?);
            Ok((hash, frame, pin, context, w, form_key()?))
        })();
        let (hash, frame, pin, context, w, key) = match result {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        pending.set(true);
        executor.set_value(Some(w.clone()));
        let alive = export_alive.clone();
        spawn_local(async move {
            let result = async {
                let response = w
                    .ask(Request {
                        id: 1,
                        revision: context.clone(),
                        job: AnalysisJob::ExportProtocolLayout {
                            workspace: frame,
                            layout_hash: hash.clone(),
                            format,
                        },
                    })
                    .await?;
                if !alive.load(Ordering::Relaxed) {
                    return Err("Layout closed.".into());
                }
                pin.permit(&snapshot.run(())?, epoch.get_untracked())?;
                if form_key()? != key {
                    return Err("Layout changed before export.".into());
                }
                if response.id != 1 || response.revision != context {
                    return Err("Layout export response mismatch.".into());
                }
                let output = match response.result {
                    AnalysisResult::ProtocolLayoutExported(v) => v,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected layout export result".into()),
                };
                output
                    .validate_binding(&context)
                    .map_err(|e| e.to_string())?;
                let artifact = crate::io::download::prepare_download(
                    &output.filename,
                    &output.media_type,
                    &output.bytes,
                )?;
                if output.layout_hash != hash || output.format != format {
                    return Err("Export refers to another layout.".into());
                }
                let warnings = output.warnings;
                artifact
                    .save(move || {
                        if !alive.load(Ordering::Relaxed) {
                            return Err("Layout closed before export admission.".into());
                        }
                        pin.permit(&snapshot.run(())?, epoch.get_untracked())?;
                        if form_key()? != key {
                            return Err("Layout changed before export admission.".into());
                        }
                        Ok(())
                    })
                    .await?;
                Ok::<_, String>(warnings)
            }
            .await;
            w.close();
            if pending.try_get_untracked().is_some() {
                pending.set(false);
                executor.set_value(None);
                status.set(match result {
                    Ok(warnings) => format!("Layout export acknowledged. {}", warnings.join(" · ")),
                    Err(e) => e,
                });
            }
        });
    });
    view! {<section class="protocol-layout-editor"><h3>"Plate layout"</h3><p class="helper">"Saved verification is historical evidence, separate from plate assignment validity. Resolving a well does not freshly verify or approve its construct."</p><p class="helper">"Assign constructs to wells. Reaction amounts are not set in this layout. Empty wells remain empty."</p>
    <button class="text-tool" disabled=move||pending.get() on:click=move|_|on_close.run(())>"Close layout editor"</button>
    <fieldset disabled=move||pending.get()||workspace_busy.get()><details><summary>"Saved layouts for this exact batch"</summary><p class="helper">"Load assignments explicitly, then review before saving an edited layout. No existing layout is overwritten."</p>{move||saved_layouts.get().into_iter().map(|(hash,label)|view!{<button class="text-tool" disabled=move||pending.get() on:click=move|_|{let result=snapshot.run(()).and_then(|w|native::load(&w,&hash).map_err(|e|e.to_string()));match result{Ok(layout)=>{name.set(layout.request.name);plate.set(Some(layout.request.plate));assignments.set(layout.request.assignments);proposed.set(None);status.set("Saved assignments loaded; review before saving or exporting this edited declaration.".into());},Err(e)=>status.set(e)}}>{label}</button>}).collect_view()}</details>
    <label>"Layout name"<input maxlength="64" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
    <label>"Plate format"<select prop:value=move||match plate.get(){None=>"",Some(Plate::Wells96)=>"96",Some(Plate::Wells384)=>"384"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"96"=>Some(Plate::Wells96),"384"=>Some(Plate::Wells384),_=>None};if !assignments.get_untracked().is_empty(){status.set("Clear assignments before changing plate format.".into());return;}plate.set(value);proposed.set(None);selected.set(0);} ><option value="">"Choose explicitly"</option><option value="96">"96 wells (8 × 12)"</option><option value="384">"384 wells (16 × 24)"</option></select></label>
    <details><summary>"Deliberate row-major fill preview"</summary><p class="helper">"Choose either planned declarations or released products before previewing. Distinct released histories remain distinct identities. Names are display only. No placement occurs until Apply preview."</p><label>"Fill identity kind"<select prop:value=move||fill_kind.get() on:change=move|e|{fill_kind.set(match event_target_value(&e).as_str(){"planned"=>"planned","released"=>"released",_=>""});proposed.set(None);}><option value="">"Choose explicitly"</option><option value="planned">"Planned construct declarations"</option><option value="released">"Actual released product events"</option></select></label><label>"Replicates per construct"<input type="number" min="1" max="384" prop:value=move||replicates.get() on:input=move|e|{replicates.set(event_target_value(&e));proposed.set(None);}/></label><button class="button secondary" on:click=move|_|{let result=plate.get_untracked().ok_or("Choose a plate first".to_string()).and_then(|p|choices.with_value(|c|{let filtered=model::fill_choices(c,fill_kind.get_untracked())?;model::preview_fill(p,&filtered,&replicates.get_untracked())}));match result{Ok(v)=>proposed.set(Some(v)),Err(e)=>status.set(e)}}>"Preview fill"</button>{move||proposed.get().map(|items|{let count=items.len();view!{<p>{format!("Proposed {count} explicit assignments; Apply replaces current assignments.")}</p><button class="button secondary" on:click=move|_|{assignments.set(items.clone());proposed.set(None);}>"Apply preview"</button>}})}</details>
    {move||plate.get().map(|p|{let(rows,columns)=model::dimensions(p);view!{<div class="protocol-well-grid" style=format!("display:grid;grid-template-columns:repeat({columns},minmax(24px,1fr));gap:2px;overflow:auto")>{(0..rows*columns).map(|i|{view!{<button class="text-tool" class:active=move||selected.get()==i title=move||assignments.with(|all|all.iter().find(|a|usize::from(a.well.row)*columns+usize::from(a.well.column)==i).map(|a|format!("{} · {}",a.id,a.label)).unwrap_or_else(||"Empty well".into())) on:click=move|_|selected.set(i)>{model::well_name(&Well{row:(i/columns)as u8,column:(i%columns)as u8})}{move||if assignments.with(|all|all.iter().any(|a|usize::from(a.well.row)*columns+usize::from(a.well.column)==i)){" •"}else{""}}</button>}}).collect_view()}</div>}})}
    <label>"Exact saved construct"<select prop:value=move||source.get().to_string() on:change=move|e|source.set(event_target_value(&e).parse().unwrap_or(0))>{choices.with_value(|all|all.iter().enumerate().map(|(i,c)|view!{<option value=i.to_string()>{format!("{} · {}",c.label.chars().take(40).collect::<String>(),model::compact_identity(&c.source))}</option>}).collect_view())}</select></label><label>"Optional presentation label"<input maxlength="64" prop:value=move||label.get() on:input=move|e|label.set(event_target_value(&e))/></label>
    <button class="button secondary" on:click=move|_|{let Some(p)=plate.get_untracked()else{status.set("Choose a plate first.".into());return;};let(_,columns)=model::dimensions(p);let i=selected.get_untracked();let Some(choice)=choices.with_value(|all|all.get(source.get_untracked()).cloned())else{status.set("No exact declared construct is selected.".into());return;};assignments.update(|all|{all.retain(|a|usize::from(a.well.row)*columns+usize::from(a.well.column)!=i);all.push(Assignment{id:format!("well-assignment-{}",i+1),well:Well{row:(i/columns)as u8,column:(i%columns)as u8},label:label.get_untracked(),source:choice.source});});}>"Assign selected well"</button>
    <button class="text-tool" on:click=move|_|if let Some(p)=plate.get_untracked(){let(_,columns)=model::dimensions(p);assignments.update(|all|all.retain(|a|usize::from(a.well.row)*columns+usize::from(a.well.column)!=selected.get_untracked()));}>"Clear selected well"</button><button class="text-tool" on:click=move|_|{assignments.set(vec![]);proposed.set(None);}>"Clear all assignments"</button>
    <details><summary>"Exact well assignments"</summary>{move||{let all=assignments.get();let pages=all.len().max(1).div_ceil(24);let page=assignment_page.get().min(pages-1);view!{<div>{all.iter().skip(page*24).take(24).map(|a|view!{<p>{format!("{} · {} · label {:?} · {}",model::well_name(&a.well),a.id,a.label,model::source_identity(&a.source))}</p>}).collect_view()}<button class="text-tool" disabled={page==0} on:click=move|_|assignment_page.set(page.saturating_sub(1))>"Previous well assignments"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|assignment_page.set((page+1).min(pages-1))>"Next well assignments"</button></div>}}}</details>
    <p>{move||format!("{} assigned wells. Quantities: missing (identity-only layout).",assignments.with(Vec::len))}</p><button class="button secondary" disabled=move||pending.get()||workspace_busy.get() on:click=move|_|run.run(false)>"Review assignments"</button><button class="button primary" disabled=move||pending.get()||workspace_busy.get()||report.get().is_none() on:click=move|_|run.run(true)>"Save layout"</button>
    {move||report.get().map(|r|view!{<p>{format!("{} assignments checked",r.resolved.len())}</p><ul>{r.warnings.iter().map(|w|view!{<li>{w.clone()}</li>}).collect_view()}</ul>})}
    <Show when=move||saved.get().is_some()><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::Csv)>"Export CSV"</button><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::Html)>"Export HTML"</button><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::Pdf)>"Export PDF"</button></Show></fieldset><p role="status">{move||status.get()}</p></section>}
}
