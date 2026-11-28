//! Documentary quantities: exact arithmetic, immutable event pins and native ACK.
use crate::platform::worker::CachedWorker;
use crate::protocol::protocol_quantities_model as model;
use collagene_core::{
    documents::workspace::Workspace,
    output::quantities as native,
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
pub fn ReactionQuantitiesView(
    capture: Callback<(), Result<model::Binding, String>>,
    identity: Memo<(String, String)>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    // Capture only small immutable identities for display; workspace validation
    // and fingerprints are deferred to an explicit action.
    let initial = capture.run(());
    let binding = StoredValue::new(initial.as_ref().ok().map(|b| Arc::new(b.clone())));
    let scope = StoredValue::new(None::<model::Scope>);
    let records = RwSignal::new(Vec::<(String, String, bool, Vec<native::Finding>)>::new());
    let status = RwSignal::new(initial.err().unwrap_or_default());
    let form = RwSignal::new(model::Form::default());
    let reviewed = RwSignal::new(None::<(model::Form, Arc<native::Preview>)>);
    let historical = RwSignal::new(None::<(String, Arc<native::Preview>)>);
    let pending = RwSignal::new(false);
    let revision_source = RwSignal::new(None::<String>);
    let copy_source = RwSignal::new(String::new());
    let copy_targets = RwSignal::new(Vec::<String>::new());
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let counter = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let tx = StoredValue::new_local(Rc::new(RefCell::new(
        None::<crate::platform::storage::Transaction>,
    )));
    let lease = StoredValue::new_local(Rc::new(
        crate::platform::commit_lease::CommitLease::default(),
    ));
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(t) = tx.with_value(|t| t.borrow().clone()) {
            let _ = t.abort();
        } else if lease.with_value(|l| l.cancel_unstarted()) {
            workspace_busy.set(false);
        }
    });
    let check = Callback::new(move |()| -> Result<_, String> {
        let current = capture.run(())?;
        let initial = binding
            .get_value()
            .ok_or("No exact quantity plan context is available")?;
        if !model::same_binding(&initial, &current)? {
            return Err(
                "The displayed plan or source/target events changed; reopen quantities.".into(),
            );
        }
        let frame = snapshot.run(())?;
        let held = if let Some(held) = scope.get_value() {
            held
        } else {
            let held = model::Scope::capture(&frame, &current, epoch.get_untracked())?;
            scope.set_value(Some(held.clone()));
            held
        };
        held.permit(&frame, &current, epoch.get_untracked())?;
        Ok((current, frame))
    });
    Effect::new(move || {
        let _ = (identity.get(), epoch.get());
        counter.with_value(|n| n.set(n.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        reviewed.set(None);
        historical.set(None);
        revision_source.set(None);
        copy_source.set(String::new());
        copy_targets.set(Vec::new());
        records.set(Vec::new());
        scope.set_value(None);
        // Parent identity/epoch changes invalidate proof, but this mounted form
        // can capture the new exact context without retaining its old receipt.
        match capture.run(()) {
            Ok(current) => binding.set_value(Some(Arc::new(current))),
            Err(error) => {
                binding.set_value(None);
                status.set(error);
            }
        }
        pending.set(false);
    });
    let run_alive = alive.clone();
    let run = Callback::new(move |save: bool| {
        if pending.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let preparation = (|| -> Result<_, String> {
            let (b, frame) = check.run(())?;
            let draft = form.get_untracked();
            let request = draft.request(&b)?;
            let context = native::preview_revision(&frame, &request).map_err(|e| e.to_string())?;
            if save
                && !reviewed.with_untracked(|r| {
                    r.as_ref()
                        .is_some_and(|(f, p)| f == &draft && p.context_revision == context)
                })
            {
                return Err("Review these quantities before saving.".into());
            }
            let database = if save {
                Some(
                    db.get_value()
                        .ok_or("Native recovery storage is unavailable")?,
                )
            } else {
                None
            };
            let w = Rc::new(CachedWorker::new()?);
            Ok((frame, draft, request, context, database, w))
        })();
        let (frame, draft, request, context, database, w) = match preparation {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let token = counter.with_value(|n| {
            let next = n.get().wrapping_add(1);
            n.set(next);
            next
        });
        if save {
            if let Err(e) = lease.with_value(|l| l.acquire(token)) {
                status.set(e);
                return;
            }
            workspace_busy.set(true);
            save_generation.update(|v| *v = v.wrapping_add(1));
        }
        pending.set(true);
        worker.set_value(Some(w.clone()));
        let alive = run_alive.clone();
        let nonce = counter.get_value();
        let lock = lease.get_value();
        let pending_tx = tx.get_value();
        spawn_local(async move {
            let guard = || -> Result<(), String> {
                if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                    return Err("Quantity operation cancelled.".into());
                }
                check.run(())?;
                if form.get_untracked() != draft {
                    return Err("Quantity declarations changed; review again.".into());
                }
                Ok(())
            };
            let result = async {
                let job = if save {
                    AnalysisJob::CaptureProtocolQuantities {
                        workspace: frame.clone(),
                        request,
                    }
                } else {
                    AnalysisJob::PreviewProtocolQuantities {
                        workspace: frame.clone(),
                        request,
                    }
                };
                let response = w
                    .ask(Request {
                        id: token as u32,
                        revision: context.clone(),
                        job,
                    })
                    .await?;
                guard()?;
                if response.id != token as u32 || response.revision != context {
                    return Err("Quantity response context mismatch.".into());
                }
                match response.result {
                    AnalysisResult::ProtocolQuantitiesPreview(p) if !save => {
                        p.validate_binding(&context).map_err(|e| e.to_string())?;
                        reviewed.set(Some((draft.clone(), Arc::new(p))));
                        Ok(())
                    }
                    AnalysisResult::ProtocolQuantitiesCaptured(commit) if save => {
                        commit
                            .validate_binding(&context)
                            .map_err(|e| e.to_string())?;
                        if commit.source_workspace_fingerprint
                            != scope.get_value().ok_or("Lost quantity scope")?.fingerprint
                        {
                            return Err("Quantity capture source changed.".into());
                        }
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
                        // Restoration deliberately clears current science and
                        // may dispose this view. Feedback belongs to the parent.
                        message.set("Quantities saved. Reopen the saved run or review the retained candidate to open or export them.".into());
                        Ok(())
                    }
                    AnalysisResult::Error(e) => Err(e),
                    _ => Err("Unexpected quantity worker result.".into()),
                }
            }
            .await;
            w.close();
            if save && lock.release(token) {
                workspace_busy.set(false);
            }
            if alive.load(Ordering::Relaxed) && nonce.get() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok(()) => {
                        if save {
                            "Quantities saved.".into()
                        } else {
                            "Quantities reviewed. Missing amounts remain unspecified.".into()
                        }
                    }
                    Err(e) => e,
                });
            }
        });
    });
    let list_alive = alive.clone();
    let list = move |_| {
        if pending.get_untracked() {
            return;
        }
        let ready = (|| -> Result<_, String> {
            let (b, frame) = check.run(())?;
            let context = native::list_revision(&frame).map_err(|e| e.to_string())?;
            Ok((b, frame, context, Rc::new(CachedWorker::new()?)))
        })();
        let (b, frame, context, w) = match ready {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let token = counter.with_value(|n| {
            let t = n.get().wrapping_add(1);
            n.set(t);
            t
        });
        let nonce = counter.get_value();
        let alive = list_alive.clone();
        pending.set(true);
        worker.set_value(Some(w.clone()));
        spawn_local(async move {
            let result=async{let response=w.ask(Request{id:token as u32,revision:context.clone(),job:AnalysisJob::ListProtocolQuantities{workspace:frame}}).await?;
            if !alive.load(Ordering::Relaxed)||nonce.get()!=token{return Err("Quantity listing cancelled.".into());}check.run(())?;if response.id!=token as u32||response.revision!=context{return Err("Quantity list context changed.".into());}
            let listed=match response.result{AnalysisResult::ProtocolQuantitiesListed(v)=>v,AnalysisResult::Error(e)=>return Err(e),_=>return Err("Unexpected quantity list result.".into())};listed.validate_binding(&context).map_err(|e|e.to_string())?;
            let total=listed.entries.len();let mut matching=vec![];for row in listed.entries{if model::summary_matches(&b,&row)?{matching.push((row.attachment_hash,row.name,row.available,row.findings));}}
            let shown=matching.len();records.set(matching);Ok(format!("{shown} of {total} saved quantity documents match this exact plan source and target/source events."))
        }.await;
            w.close();
            if alive.load(Ordering::Relaxed) && nonce.get() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok(v) => v,
                    Err(e) => e,
                });
            }
        });
    };
    let open_alive = alive.clone();
    let open = Callback::new(move |hash: String| {
        if pending.get_untracked() {
            return;
        }
        let start = (|| -> Result<_, String> {
            let (_, frame) = check.run(())?;
            let context = native::open_revision(&frame, &hash).map_err(|e| e.to_string())?;
            Ok((frame, context, Rc::new(CachedWorker::new()?)))
        })();
        let (frame, context, w) = match start {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        pending.set(true);
        worker.set_value(Some(w.clone()));
        let token = counter.with_value(|n| {
            let next = n.get().wrapping_add(1);
            n.set(next);
            next
        });
        let nonce = counter.get_value();
        let alive = open_alive.clone();
        spawn_local(async move {
            let result = async {
                let response = w
                    .ask(Request {
                        id: token as u32,
                        revision: context.clone(),
                        job: AnalysisJob::OpenProtocolQuantities {
                            workspace: frame,
                            quantity_hash: hash.clone(),
                        },
                    })
                    .await?;
                if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                    return Err("Quantity open cancelled.".into());
                }
                check.run(())?;
                if response.id != token as u32 || response.revision != context {
                    return Err("Quantity open context mismatch.".into());
                }
                let opened = match response.result {
                    AnalysisResult::ProtocolQuantitiesOpened(v) => v,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected quantity open result.".into()),
                };
                opened
                    .validate_binding(&context)
                    .map_err(|e| e.to_string())?;
                if !model::request_matches(&capture.run(())?, &opened.preview.request)? {
                    return Err(
                        "Saved quantities belong to another exact plan/source/target binding."
                            .into(),
                    );
                }
                if opened.quantity_hash != hash {
                    return Err("Another quantity document was returned.".into());
                }
                historical.set(Some((hash, Arc::new(opened.preview))));
                Ok(())
            }
            .await;
            w.close();
            if alive.load(Ordering::Relaxed) && nonce.get() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok(()) => "Saved quantities opened.".into(),
                    Err(e) => e,
                });
            }
        });
    });
    let revise = move |_| {
        if pending.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let result = (|| -> Result<_, String> {
            let (binding, _) = check.run(())?;
            let (hash, preview) = historical
                .get_untracked()
                .ok_or("Open a saved quantity document first.")?;
            let restored = model::Form::from_request(&preview.request, &binding)?;
            Ok((hash, restored))
        })();
        match result {
            Ok((hash, restored)) => {
                form.set(restored);
                reviewed.set(None);
                revision_source.set(Some(hash));
                copy_targets.set(vec![]);
                status.set("Saved declarations loaded for editing. Review and save a new report; the original remains retained.".into());
            }
            Err(e) => status.set(e),
        }
    };
    let apply_quantities = move |_| {
        if pending.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let result = (|| -> Result<_, String> {
            let (binding, _) = check.run(())?;
            let mut draft = form.get_untracked();
            draft.apply_to_reactions(
                &binding,
                &copy_source.get_untracked(),
                &copy_targets.get_untracked(),
            )?;
            Ok(draft)
        })();
        match result {
            Ok(draft) => {
                form.set(draft);
                reviewed.set(None);
                status.set("Declarations copied to the selected reactions. Review quantities before saving; molecular verification is unchanged.".into());
            }
            Err(e) => status.set(e),
        }
    };
    let export_alive = alive.clone();
    let export = Callback::new(move |format: native::ExportFormat| {
        if pending.get_untracked() {
            return;
        }
        let start = (|| -> Result<_, String> {
            let (hash, preview) = historical
                .get_untracked()
                .ok_or("Open or save a quantity document first")?;
            let (_, frame) = check.run(())?;
            let context =
                native::export_revision(&frame, &hash, format).map_err(|e| e.to_string())?;
            let combined = matches!(
                format,
                native::ExportFormat::ProtocolHtml | native::ExportFormat::ProtocolPdf
            );
            let fingerprint = combined
                .then(|| collagene_core::material::lineage::workspace_fingerprint(&frame))
                .transpose()
                .map_err(|e| e.to_string())?;
            Ok((
                hash,
                preview,
                frame,
                context,
                fingerprint,
                Rc::new(CachedWorker::new()?),
            ))
        })();
        let (hash, preview, frame, context, fingerprint, w) = match start {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let nonce = counter.get_value();
        let token = nonce.get().wrapping_add(1);
        nonce.set(token);
        pending.set(true);
        worker.set_value(Some(w.clone()));
        let alive = export_alive.clone();
        spawn_local(async move {
            let result = async {
                let response = w
                    .ask(Request {
                        id: 1,
                        revision: context.clone(),
                        job: AnalysisJob::ExportProtocolQuantities {
                            workspace: frame,
                            quantity_hash: hash.clone(),
                            format,
                        },
                    })
                    .await?;
                if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                    return Err("Quantity export context changed.".into());
                }
                check.run(())?;
                if response.id != 1 || response.revision != context {
                    return Err("Quantity export context mismatch.".into());
                }
                let output = match response.result {
                    AnalysisResult::ProtocolQuantitiesExported(v) => v,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected quantity export result.".into()),
                };
                output
                    .validate_binding(&context)
                    .map_err(|e| e.to_string())?;
                output.validate_for_preview(&preview).map_err(|e| e.to_string())?;
                if output.source_workspace_fingerprint != fingerprint {
                    return Err("Combined protocol workspace context changed.".into());
                }
                if output.quantity_hash != hash || output.format != format {
                    return Err("Exported quantity identity/format changed.".into());
                }
                let summary = match &output.molecular {
                    Some(molecular) => {
                        let status = match molecular.status {
                            collagene_core::verification::evidence::Status::Invalid => "inconsistent",
                            collagene_core::verification::evidence::Status::Incomplete => "incomplete",
                            collagene_core::verification::evidence::Status::ConsistentUnderAssumptions => "consistent under stated assumptions",
                        };
                        format!("Protocol exported. Sequence check: {status}. Quantity arithmetic complete: {}.", output.arithmetic_complete)
                    }
                    None => format!("Quantities exported. Arithmetic complete: {}.", output.arithmetic_complete),
                };
                let artifact = crate::io::download::prepare_download(
                    &output.filename,
                    &output.media_type,
                    &output.bytes,
                )?;
                let admission_nonce = nonce.clone();
                artifact
                    .save(move || {
                        if !alive.load(Ordering::Relaxed) || admission_nonce.get() != token {
                            return Err("Quantity export context changed before admission.".into());
                        }
                        let (_, current) = check.run(())?;
                        if let Some(expected) = &fingerprint {
                            if collagene_core::material::lineage::workspace_fingerprint(&current).map_err(|e| e.to_string())? != *expected {
                                return Err("Combined protocol workspace changed before export admission.".into());
                            }
                        }
                        Ok(())
                    })
                    .await?;
                Ok::<_, String>((summary, output.warnings))
            }
            .await;
            w.close();
            if pending.try_get_untracked().is_some() && nonce.get() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok((summary, warnings)) => format!("{summary} {}", warnings.join(" · ")),
                    Err(e) => e,
                });
            }
        });
    });
    view! {<details class="assembly-fragment-options"><summary>"Reaction quantities"</summary><p class="helper">"Amounts stay tied to this plan. Missing amounts remain unspecified; these calculations do not change sequence verification."</p>
          <details><summary>"Exact plan and event identities"</summary>{binding.with_value(|b|b.as_ref().map(|b|view!{<p>{format!("Plan {} · independent target document {} · event {} · revision {}",b.plan_revision,b.target.document_id,b.target.event_id,b.target.revision)}</p><ul>{b.inputs.iter().map(|i|view!{<li>{format!("{}: document {} · event {} · revision {}",i.input_id,i.event.document_id,i.event.event_id,i.event.revision)}</li>}).collect_view()}</ul>}))}</details>
    <fieldset disabled=move||pending.get()||workspace_busy.get()><label>"Quantity document name"<input maxlength="256" prop:value=move||form.with(|f|f.name.clone()) on:input=move|e|form.update(|f|f.name=event_target_value(&e))/></label><button class="button secondary" on:click=move|_|{let mut result=Ok(());form.update(|f|result=f.add_reaction());if let Err(e)=result{status.set(e);}}>"Add explicit reaction"</button>
       <For each=move||form.with(|f|f.reactions.iter().map(|r|r.id.clone()).collect::<Vec<_>>()) key=|id|id.clone() children=move|reaction_id|view!{<ReactionEditor form reaction_id binding=binding.get_value()/>}/>
       <details><summary>"Apply quantities to selected reactions"</summary>
       <p class="helper">"Copy the source reaction's amount, concentration, rounding and fill declarations. Destination rows are replaced only when every selected operation has the same type and each material operand maps unambiguously by position. No compatibility, concentration measurement or experimental success is inferred."</p>
       <label>"Copy from reaction"<select aria-label="Quantity copy source" prop:value=move||copy_source.get() on:change=move|e|{copy_source.set(event_target_value(&e));copy_targets.set(vec![]);}><option value="">"Choose explicitly"</option>{move||form.with(|f|f.reactions.iter().map(|r|view!{<option value=r.id.clone()>{format!("{} · operation {}",r.id,r.step.parse::<usize>().map(|v|(v+1).to_string()).unwrap_or_else(|_|"unselected".into()))}</option>}).collect_view())}</select></label>
       <For each=move||form.with(|f|f.reactions.iter().map(|r|r.id.clone()).collect::<Vec<_>>()) key=|id|id.clone() children=move|id|{
            let held=StoredValue::new(id);
            view!{<label class="checkbox-label"><input type="checkbox" disabled=move||held.get_value()==copy_source.get() prop:checked=move||copy_targets.with(|v|v.contains(&held.get_value())) on:change=move|e|{let selected=event_target_checked(&e);copy_targets.update(|values|{values.retain(|id|id!=&held.get_value());if selected{values.push(held.get_value());}});}/>{format!("Replace {} declarations",held.get_value())}</label>}
       }/>
       <button class="button secondary" disabled=move||copy_source.get().is_empty()||copy_targets.get().is_empty() on:click=apply_quantities>"Apply to selected reactions"</button></details>
       <p class="helper">"Amounts accept decimal values or exact fractions such as 1/3. Blank values stay unspecified."</p>
       <button class="button secondary" on:click=move|_|run.run(false)>"Review quantities"</button><button class="button primary" disabled=move||reviewed.with(|r|r.as_ref().is_none_or(|(f,_)|f!=&form.get())) on:click=move|_|run.run(true)>{move||if revision_source.get().is_some(){"Save as new quantity report"}else{"Save quantities"}}</button></fieldset>
       {move||reviewed.get().map(|(_,preview)|view!{<QuantityReport preview/>})}
       <details><summary>"Saved quantities"</summary><button class="button secondary" disabled=move||pending.get() on:click=list>"Load saved quantities"</button>{move||records.get().into_iter().map(|(hash,name,available,findings)|view!{<article><button class="text-tool" disabled={move||pending.get()||!available} on:click=move|_|open.run(hash.clone())>{name}</button><p>{if available{"Plan and source events available"}else{"Unavailable: retained plan or source events are missing. Saved quantities remain retained."}}</p>{findings.into_iter().map(|f|view!{<p>{format!("{}: {}",f.path,f.message)}</p>}).collect_view()}</article>}).collect_view()}</details>
       {move||historical.get().map(|(_,preview)|view!{<h4>"Saved quantities"</h4><QuantityReport preview/>})}
       <Show when=move||historical.get().is_some()><button class="button secondary" disabled=move||pending.get()||workspace_busy.get() on:click=revise>"Revise saved quantities as new report"</button><p class="helper">"Export quantities alone, or include the actual protocol instructions with a fresh sequence check. Quantity completeness and sequence verification remain separate."</p><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::ProtocolHtml)>"Export protocol + quantities HTML"</button><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::ProtocolPdf)>"Export protocol + quantities PDF"</button><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::Csv)>"Export quantities CSV"</button><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::Html)>"Export quantities HTML"</button><button class="button secondary" disabled=move||pending.get() on:click=move|_|export.run(native::ExportFormat::Pdf)>"Export quantities PDF"</button></Show><p role="status">{move||status.get()}</p></details>}
}
fn edit_reaction(
    form: RwSignal<model::Form>,
    id: &str,
    edit: impl FnOnce(&mut model::ReactionForm),
) {
    form.update(|f| {
        if let Some(r) = f.reactions.iter_mut().find(|r| r.id == id) {
            edit(r)
        }
    });
}
fn edit_component(
    form: RwSignal<model::Form>,
    reaction: &str,
    id: &str,
    edit: impl FnOnce(&mut model::ComponentForm),
) {
    edit_reaction(form, reaction, |r| {
        if let Some(c) = r.components.iter_mut().find(|c| c.id == id) {
            edit(c)
        }
    });
}
#[component]
fn ReactionEditor(
    form: RwSignal<model::Form>,
    reaction_id: String,
    binding: Option<Arc<model::Binding>>,
) -> impl IntoView {
    let id = StoredValue::new(reaction_id);
    let note = RwSignal::new(String::new());
    let binding = StoredValue::new(binding);
    let row = move || form.with(|f| f.reactions.iter().find(|r| r.id == id.get_value()).cloned());
    view! {<details open><summary>{id.get_value()}</summary><label>"Canonical operation"<select prop:value=move||row().map(|r|r.step).unwrap_or_default() on:change=move|e|edit_reaction(form,&id.get_value(),|r|r.step=event_target_value(&e))><option value="">"Choose explicitly"</option>{binding.with_value(|b|b.as_ref().map(|b|(0..b.step_count).map(|i|view!{<option value=i.to_string()>{format!("Operation {}",i+1)}</option>}).collect_view()))}</select></label><label>"Final reaction volume"<input maxlength="64" prop:value=move||row().map(|r|r.total).unwrap_or_default() on:input=move|e|edit_reaction(form,&id.get_value(),|r|r.total=event_target_value(&e))/></label><label>"Volume unit"<select prop:value=move||row().map(|r|r.unit).unwrap_or_default() on:change=move|e|edit_reaction(form,&id.get_value(),|r|r.unit=event_target_value(&e))><option value="nL">"nL"</option><option value="uL">"µL"</option><option value="mL">"mL"</option></select></label>
    <label class="checkbox-label"><input type="checkbox" prop:checked=move||row().is_some_and(|r|r.fill) on:change=move|e|edit_reaction(form,&id.get_value(),|r|r.fill=event_target_checked(&e))/>"Explicitly allocate remaining volume to a named fill component"</label><Show when=move||row().is_some_and(|r|r.fill)><label>"Fill component ID"<input maxlength="128" prop:value=move||row().map(|r|r.fill_id).unwrap_or_default() on:input=move|e|edit_reaction(form,&id.get_value(),|r|r.fill_id=event_target_value(&e))/></label><label>"Fill label (no default reagent)"<input maxlength="256" prop:value=move||row().map(|r|r.fill_label).unwrap_or_default() on:input=move|e|edit_reaction(form,&id.get_value(),|r|r.fill_label=event_target_value(&e))/></label></Show>
    <For each=move||row().map(|r|r.components.into_iter().map(|c|c.id).collect::<Vec<_>>()).unwrap_or_default() key=|id|id.clone() children=move|component_id|view!{<ComponentEditor form reaction_id=id.get_value() component_id binding=binding.get_value()/>}/>
    <button class="button secondary" on:click=move|_|{let mut result:Result<(),String>=Err("Reaction no longer exists".into());form.update(|f|if let Some(i)=f.reactions.iter().position(|r|r.id==id.get_value()){result=f.add_component(i)});if let Err(e)=result{note.set(e);}}>"Add component with blank amount"</button><button class="text-tool" on:click=move|_|form.update(|f|f.reactions.retain(|r|r.id!=id.get_value()))>"Remove reaction"</button><p role="status">{move||note.get()}</p></details>}
}
#[component]
fn ComponentEditor(
    form: RwSignal<model::Form>,
    reaction_id: String,
    component_id: String,
    binding: Option<Arc<model::Binding>>,
) -> impl IntoView {
    let reaction = StoredValue::new(reaction_id);
    let id = StoredValue::new(component_id);
    let binding = StoredValue::new(binding);
    let component = move || {
        form.with(|f| {
            f.reactions
                .iter()
                .find(|r| r.id == reaction.get_value())
                .and_then(|r| r.components.iter().find(|c| c.id == id.get_value()))
                .cloned()
        })
    };
    let row_step = move || {
        form.with(|f| {
            f.reactions
                .iter()
                .find(|r| r.id == reaction.get_value())
                .and_then(|r| r.step.parse::<usize>().ok())
        })
    };
    view! {<details><summary>{id.get_value()}</summary><label>"Component label"<input maxlength="256" prop:value=move||component().map(|c|c.label).unwrap_or_default() on:input=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.label=event_target_value(&e))/></label><label>"Optional planned material reference"<select prop:value=move||component().map(|c|c.source).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.source=event_target_value(&e))><option value="">"No biological reference"</option>{move||binding.with_value(|b|b.as_ref().and_then(|b|row_step().and_then(|step|b.component_sources.get(step))).map(|sources|sources.iter().map(|source|match source{collagene_core::planning::plan::Source::Input(input)=>view!{<option value=format!("input:{input}")>{format!("Input {input}")}</option>},collagene_core::planning::plan::Source::Step(step)=>view!{<option value=format!("step:{step}")>{format!("Planned operation {} output",step+1)}</option>}}).collect_view()))}</select></label><p class="helper">"A planned source reference does not assert that material was isolated or that a concentration was measured."</p>
        <label>"Amount declaration"<select prop:value=move||component().map(|c|c.mode).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|{c.mode=event_target_value(&e);c.concentration_unit=if c.mode=="mass"{"ng/uL".into()}else{"nM".into()};})><option value="direct">"Direct volume"</option><option value="mass">"Declared mass-concentration dilution"</option><option value="molar">"Declared molar-concentration dilution"</option></select></label>
        <Show when=move||component().is_some_and(|c|c.mode=="direct")><label>"Component volume (blank = missing)"<input maxlength="64" prop:value=move||component().map(|c|c.volume).unwrap_or_default() on:input=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.volume=event_target_value(&e))/></label></Show>
        <label>"Volume / increment unit"<select prop:value=move||component().map(|c|c.unit).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.unit=event_target_value(&e))><option value="nL">"nL"</option><option value="uL">"µL"</option><option value="mL">"mL"</option></select></label>
        <Show when=move||component().is_some_and(|c|c.mode!="direct")><p class="helper">"Stock and target concentrations are declarations. Declare the dilution final volume explicitly; it must match the reaction total for a complete calculation. Only the computed stock addition counts toward the mixture. No automatic diluent is inserted."</p><label>"Dilution final volume (blank = missing)"<input maxlength="64" aria-label="Dilution final volume" prop:value=move||component().map(|c|c.dilution_volume).unwrap_or_default() on:input=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.dilution_volume=event_target_value(&e))/></label>
        <label>"Dilution volume unit"<select prop:value=move||component().map(|c|c.dilution_unit).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.dilution_unit=event_target_value(&e))><option value="nL">"nL"</option><option value="uL">"µL"</option><option value="mL">"mL"</option></select></label>
        <button class="button secondary" on:click=move|_|{let values=form.with_untracked(|f|f.reactions.iter().find(|r|r.id==reaction.get_value()).map(|r|(r.total.clone(),r.unit.clone())));if let Some((volume,unit))=values{edit_component(form,&reaction.get_value(),&id.get_value(),|c|{c.dilution_volume=volume;c.dilution_unit=unit;});}}>"Copy declared reaction total"</button>
    <label>"Stock concentration"<input maxlength="64" prop:value=move||component().map(|c|c.stock).unwrap_or_default() on:input=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.stock=event_target_value(&e))/></label><label>"Target concentration"<input maxlength="64" prop:value=move||component().map(|c|c.target).unwrap_or_default() on:input=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.target=event_target_value(&e))/></label><Show when=move||component().is_some_and(|c|c.mode=="mass")><label>"Mass concentration unit"<select prop:value=move||component().map(|c|c.concentration_unit).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.concentration_unit=event_target_value(&e))><option value="ng/uL">"ng/µL"</option><option value="ng/nL">"ng/nL"</option></select></label></Show><Show when=move||component().is_some_and(|c|c.mode=="molar")><label>"Molar concentration unit"<select prop:value=move||component().map(|c|c.concentration_unit).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.concentration_unit=event_target_value(&e))><option value="nM">"nM"</option><option value="uM">"µM"</option><option value="pmol/nL">"pmol/nL"</option></select></label></Show>
        <label>"Rounding policy"<select prop:value=move||component().map(|c|c.rounding).unwrap_or_default() on:change=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.rounding=event_target_value(&e))><option value="">"Exact arithmetic; do not quantize"</option><option value="exact">"Require exact increment multiple"</option><option value="floor">"Round down"</option><option value="ceiling">"Round up"</option><option value="nearest">"Nearest; ties up"</option></select></label><Show when=move||component().is_some_and(|c|!c.rounding.is_empty())><label>"Explicit volume increment"<input maxlength="64" prop:value=move||component().map(|c|c.increment).unwrap_or_default() on:input=move|e|edit_component(form,&reaction.get_value(),&id.get_value(),|c|c.increment=event_target_value(&e))/></label></Show></Show>
        <button class="text-tool" on:click=move|_|edit_reaction(form,&reaction.get_value(),|r|r.components.retain(|c|c.id!=id.get_value()))>"Remove component"</button></details>}
}
#[component]
fn QuantityReport(preview: Arc<native::Preview>) -> impl IntoView {
    let saved = StoredValue::new(preview);
    let page = RwSignal::new(0usize);
    let total = saved.with_value(|p| p.reactions.len());
    let pages = total.max(1).div_ceil(8);
    let finding_page = RwSignal::new(0usize);
    let finding_count = saved.with_value(|p| p.findings.len());
    let finding_pages = finding_count.max(1).div_ceil(40);
    view! {<section><p>{saved.with_value(|p|format!("{} · arithmetic complete {} · context {}",p.plan_name,p.complete,p.context_revision))}</p>{saved.with_value(|p|p.warnings.iter().map(|w|view!{<p class="helper">{w.clone()}</p>}).collect_view())}<details><summary>{format!("Calculation findings ({finding_count})")}</summary>
    {move||saved.with_value(|p|p.findings.iter().skip(finding_page.get().min(finding_pages-1)*40).take(40).map(|f|view!{<p>{format!("{}: {}",f.path,f.message)}</p>}).collect_view())}
    <Show when={move||finding_pages>1}><p>{move||format!("Page {} of {}",finding_page.get().min(finding_pages-1)+1,finding_pages)}</p><button class="text-tool" disabled=move||finding_page.get()==0 on:click=move|_|finding_page.update(|p|*p=p.saturating_sub(1))>"Previous findings"</button><button class="text-tool" disabled=move||finding_page.get()+1>=finding_pages on:click=move|_|finding_page.update(|p|*p=(*p+1).min(finding_pages-1))>"Next findings"</button></Show></details>
    {move||saved.with_value(|p|p.reactions.iter().enumerate().skip(page.get().min(pages-1)*8).take(8).map(|(index,_)|view!{<ReactionReport preview=saved.get_value() index/>}).collect_view())}<p>{move||format!("{} reactions · page {} of {}",total,page.get().min(pages-1)+1,pages)}</p><button class="text-tool" disabled={move||page.get()==0} on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous reaction calculations"</button><button class="text-tool" disabled={move||page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next reaction calculations"</button></section>}
}
#[component]
fn ReactionReport(preview: Arc<native::Preview>, index: usize) -> impl IntoView {
    let saved = StoredValue::new(preview);
    let page = RwSignal::new(0usize);
    let total = saved.with_value(|p| p.reactions[index].components.len());
    let pages = total.max(1).div_ceil(16);
    view! {<details><summary>{saved.with_value(|p|format!("{} · operation {}",p.reactions[index].id,p.reactions[index].step_index+1))}</summary>{saved.with_value(|p|{let r=&p.reactions[index];view!{<p>{format!("Arithmetic complete {}",r.complete)}</p>{r.findings.iter().map(|f|view!{<p>{format!("{}: {}",f.path,f.message)}</p>}).collect_view()}{r.fill.as_ref().map(|f|view!{<p>{format!("Explicit fill {} ({}): {}",f.label,f.component_id,model::amount(f.volume))}</p>})}{r.mixture.as_ref().map(|m|view!{<p>{format!("Specified volume {} · remaining volume {}",model::amount(m.specified_volume),model::amount(m.remaining_diluent))}</p>})}}})}
    {move||saved.with_value(|p|p.reactions[index].components.iter().skip(page.get().min(pages-1)*16).take(16).map(|c|view!{<article><p>{format!("{} · {}",c.component_id,c.volume.map(model::amount).unwrap_or_else(||"Volume missing".into()))}</p>{c.findings.iter().map(|f|view!{<p>{format!("{}: {}",f.path,f.message)}</p>}).collect_view()}{c.calculation.as_ref().map(|value|view!{<details><summary>"Exact dilution / rounding / achieved concentration"</summary><pre>{dilution_summary(value)}</pre></details>})}</article>}).collect_view())}
    <p>{move||format!("{} components · page {} of {}",total,page.get()+1,pages)}</p><button class="text-tool" disabled={move||page.get()==0} on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous component calculations"</button><button class="text-tool" disabled={move||page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next component calculations"</button></details>}
}

fn dilution_summary(result: &native::DilutionResult) -> String {
    use collagene_core::output::quantities::units;
    fn mass(value: Option<units::MassConcentration>) -> String {
        value
            .map(|v| {
                v.in_units(units::MassUnit::Nanogram, units::VolumeUnit::Microliter)
                    .map(|n| format!("{n} ng/µL"))
                    .unwrap_or_else(|_| format!("{} ng/nL", v.canonical()))
            })
            .unwrap_or_else(|| "unspecified".into())
    }
    fn molar(value: Option<units::MolarConcentration>) -> String {
        value
            .map(|v| {
                v.in_unit(units::MolarUnit::Nanomolar)
                    .map(|n| format!("{n} nM"))
                    .unwrap_or_else(|_| format!("{} pmol/nL", v.canonical()))
            })
            .unwrap_or_else(|| "unspecified".into())
    }
    fn exact<C: Copy>(
        plan: &units::DilutionPlan<C>,
        concentration: fn(Option<C>) -> String,
    ) -> String {
        format!(
            "Stock concentration: {}\nTarget concentration: {}\nFinal volume: {}\nExact stock addition: {}\nRemaining volume before other components: {}",
            concentration(plan.recipe.stock),
            concentration(plan.recipe.target),
            plan.recipe
                .final_volume
                .map(model::amount)
                .unwrap_or_else(|| "unspecified".into()),
            model::amount(plan.stock_volume),
            model::amount(plan.diluent_volume)
        )
    }
    fn rounded<C: Copy>(
        plan: &units::QuantizedDilution<C>,
        concentration: fn(Option<C>) -> String,
    ) -> String {
        format!(
            "{}\nRounding: {:?} to {}\nRounded stock addition: {}\nAbsolute volume deviation: {}\nAchieved concentration: {}\nRequested concentration reached: {}\nRemaining volume after rounding, before other components: {}",
            exact(&plan.exact, concentration),
            plan.stock.mode,
            model::amount(plan.stock.increment),
            model::amount(plan.stock.rounded),
            model::amount(plan.stock.difference),
            concentration(Some(plan.achieved)),
            plan.target_reached,
            model::amount(plan.diluent_volume)
        )
    }
    match result {
        native::DilutionResult::Mass(p) => exact(p, mass),
        native::DilutionResult::Molar(p) => exact(p, molar),
        native::DilutionResult::QuantizedMass(p) => rounded(p, mass),
        native::DilutionResult::QuantizedMolar(p) => rounded(p, molar),
    }
}
