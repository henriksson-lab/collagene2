//! Explicit event bindings for retained, replayable cloning products.
use crate::{
    app::{encode_workspace, prepare_download as download},
    run_event_model::{Request as EventRequest, Role},
    storage, worker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::{self, EventBinding, InputEventBinding},
    verification::{DesignPlan, Operation, Source},
    workspace::Workspace,
};
use leptos::prelude::*;
use std::{
    collections::BTreeMap,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

fn event_label(workspace: &Workspace, event_id: &str) -> String {
    workspace
        .history
        .iter()
        .find(|event| event.id == event_id)
        .map(|event| {
            let name = workspace
                .revisions
                .get(&event.revision)
                .map(|m| m.name.as_str())
                .unwrap_or("Document");
            format!(
                "{name} · document {} · event {}",
                event.document_id, event.id
            )
        })
        .unwrap_or_else(|| event_id.to_string())
}
fn binding(workspace: &Workspace, event_id: &str, revision: &str) -> Result<EventBinding, String> {
    let event = workspace
        .history
        .iter()
        .find(|e| e.id == event_id && e.revision == revision)
        .ok_or("Choose an exact retained event matching this input revision.")?;
    Ok(EventBinding {
        document_id: event.document_id.clone(),
        event_id: event.id.clone(),
        revision: event.revision.clone(),
    })
}

#[component]
pub fn RunCapturePanel(
    plan: DesignPlan,
    #[prop(optional)] preferred_bindings: Vec<InputEventBinding>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let initial = snapshot.run(());
    let status = RwSignal::new(initial.as_ref().err().cloned().unwrap_or_default());
    let selected = RwSignal::new(BTreeMap::<String, String>::new());
    if let Ok(frame) = &initial {
        selected.update(|values| {
            for input in &plan.inputs {
                let matches = frame
                    .documents
                    .iter()
                    .filter(|doc| doc.head == input.revision)
                    .filter_map(|doc| doc.head_event.as_ref())
                    .collect::<Vec<_>>();
                if matches.len() == 1 {
                    values.insert(input.id.clone(), matches[0].clone());
                }
            }
        });
    }
    // A manual design can preserve its exact selected source events, including
    // distinct documents with identical contents. Validate each preference.
    if let Ok(frame) = &initial {
        selected.update(|values| {
            for preferred in &preferred_bindings {
                if plan.inputs.iter().any(|input| {
                    input.id == preferred.input_id && input.revision == preferred.event.revision
                }) && binding(frame, &preferred.event.event_id, &preferred.event.revision)
                    .is_ok_and(|event| event == preferred.event)
                {
                    values.insert(preferred.input_id.clone(), preferred.event.event_id.clone());
                }
            }
        });
    }
    let initial = StoredValue::new(initial.ok());
    let plan = StoredValue::new(plan);
    let page = RwSignal::new(0usize);
    let allow_incomplete = RwSignal::new(false);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, Ordering::Relaxed));
    let record = move |_| {
        if busy.get_untracked() {
            return;
        }
        let Some(database) = db.get_value() else {
            status.set("Browser storage is unavailable; the run has not been recorded.".into());
            return;
        };
        busy.set(true);
        save_generation.update(|g| *g = g.wrapping_add(1));
        let result = (|| {
            let frame = snapshot.run(())?;
            let plan = plan.get_value();
            let target = frame
                .documents
                .iter()
                .find(|d| d.id == frame.active_document)
                .ok_or("No independent active target document.")?;
            if target.head != plan.target_revision {
                return Err("The active target has changed. Reopen the candidate for its intended independent target before recording.".into());
            }
            let target_binding = binding(
                &frame,
                target
                    .head_event
                    .as_deref()
                    .ok_or("Target has no exact history event.")?,
                &plan.target_revision,
            )?;
            let choices = selected.get_untracked();
            let input_bindings=plan.inputs.iter().map(|input|{
                let event=choices.get(&input.id).ok_or_else(||format!("Choose the source event for {}. Identical sequence content does not identify a document.",input.id))?;
                Ok(InputEventBinding { input_id:input.id.clone(), event:binding(&frame,event,&input.revision)? })
            }).collect::<Result<Vec<_>,String>>()?;
            let fingerprint =
                run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
            Ok((frame, plan, target_binding, input_bindings, fingerprint))
        })();
        let (frame, plan, target_binding, input_bindings, fingerprint) = match result {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                busy.set(false);
                return;
            }
        };
        let original_epoch = epoch.get_untracked();
        let alive = alive.clone();
        let allow_incomplete = allow_incomplete.get_untracked();
        status.set("Reconstructing every product and checking source events…".into());
        spawn_local(async move {
            let result = async {
                let response = worker::analyze_cancellable(
                    Request {
                        id: u32::MAX - 20,
                        revision: fingerprint.clone(),
                        job: AnalysisJob::RecordCloningRun {
                            workspace: frame.clone(),
                            plan,
                            target_binding,
                            input_bindings,
                            allow_incomplete,
                        },
                    },
                    alive.clone(),
                )
                .await?;
                if response.id != u32::MAX - 20 || response.revision != fingerprint {
                    return Err("Run response does not match the submitted workspace.".into());
                }
                let commit = match response.result {
                    AnalysisResult::RunCommit(commit) => commit,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected run recording response.".into()),
                };
                if commit.receipt.source_workspace_fingerprint != fingerprint {
                    return Err("Run receipt does not bind the submitted source workspace.".into());
                }
                if commit.receipt.recorded_engine != collagene_core::verification::ENGINE_VERSION {
                    return Err("Run recording came from a different scientific engine. Reconstruct with the current worker before installing products.".into());
                }
                let guard = || -> Result<(), String> {
                    if !alive.load(Ordering::Relaxed) || epoch.get_untracked() != original_epoch {
                        return Err(
                            "Workspace changed or run dialog closed; nothing was recorded.".into(),
                        );
                    }
                    let fresh = snapshot.run(())?;
                    if run_lineage::workspace_fingerprint(&fresh).map_err(|e| e.to_string())?
                        != fingerprint
                    {
                        return Err(
                            "Workspace changed during reconstruction; nothing was recorded.".into(),
                        );
                    }
                    Ok(())
                };
                guard()?;
                validate.run(commit.workspace.clone())?;
                let current = encode_workspace(commit.workspace.clone(), u32::MAX - 21).await?;
                let previous = encode_workspace(frame, u32::MAX - 22).await?;
                guard()?;
                storage::replace_workspace(&database, &current, &previous).await?;
                Ok::<_, String>(commit.workspace)
            }
            .await;
            match result {
                Ok(frame) => {
                    loaded.set(true);
                    message.set("Cloning run recorded with reconstructed products and exact source events. The independent target is unchanged. Replay is required for current verification.".into());
                    on_restore.run(frame);
                }
                Err(error) => {
                    if alive.load(Ordering::Relaxed) {
                        status.set(format!("Run was not recorded. {error}"));
                    }
                }
            }
            busy.set(false);
        });
    };
    view! {
        <section class="run-capture-panel">
            <p>"Retain every reconstructed step as a document with its exact source events. Recording checks the plan again. It does not establish experimental provenance."</p>
            <p class="helper">"The active target is an independent comparison and remains active. It is not automatically a causal parent. Identical-content documents require an explicit source choice."</p>
            {move||initial.get_value().map(|frame|{
                let plan=plan.get_value();
                let target=frame.documents.iter().find(|d|d.id==frame.active_document).map(|d|format!("Independent target: document {} · event {}",d.id,d.head_event.as_deref().unwrap_or("missing"))).unwrap_or_else(||"No active target".into());
                view!{<p class="run-target-binding">{target}</p>
                    <div class="run-source-bindings">{plan.inputs.iter().skip(page.get()*25).take(25).map(|input|{
                        let id=input.id.clone();let change_id=id.clone();let text_id=id.clone();let selected_id=id.clone();let revision=input.revision.clone();
                        let heads=frame.documents.iter().filter(|doc|doc.head==revision).filter_map(|doc|doc.head_event.clone()).collect::<Vec<_>>();
                        let options=heads.iter().chain(frame.history.iter().filter(|e|e.revision==revision&&!heads.contains(&e.id)).map(|e|&e.id)).take(100).map(|event|(event.clone(),event_label(&frame,event))).collect::<Vec<_>>();
                        let label=format!("Source event for {id}");let exact_label=format!("Exact source event for {id}");
                        view!{<div class="run-source-binding"><strong>{id}</strong><p class="helper">{format!("Revision {revision}")}</p>
                            <select aria-label=label prop:value=move||selected.get().get(&selected_id).cloned().unwrap_or_default() on:change=move|event|{let value=event_target_value(&event);selected.update(|v|{v.insert(change_id.clone(),value);});}>
                                <option value="">"Choose a document and event explicitly"</option>
                                {options.into_iter().map(|(value,label)|view!{<option value=value>{label}</option>}).collect_view()}
                            </select>
                            <input aria-label=exact_label placeholder="Or paste an exact retained event ID" on:input=move|event|{let value=event_target_value(&event);selected.update(|v|{v.insert(text_id.clone(),value);});}/>
                            <p class="helper">"At most 100 matching events are listed; an exact event ID can select any matching retained event."</p>
                        </div>}
                    }).collect_view()}</div>
                }
            })}
            <button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous source inputs"</button>
            <button class="button secondary" disabled=move||{(page.get()+1)*25>=plan.with_value(|p|p.inputs.len())} on:click=move|_|page.update(|p|*p+=1)>"Next source inputs"</button>
            <label><input type="checkbox" aria-label="Allow incomplete run recording" on:change=move|event|allow_incomplete.set(event_target_checked(&event))/>"Allow a fully reconstructed, target-matching plan with incomplete checks; retain its incomplete status"</label>
            <button class="button primary" disabled=move||busy.get() on:click=record>"Record cloning run"</button>
            <p role="status">{move||status.get()}</p>
        </section>
    }
}

/// Export limits apply while serializing, before an unbounded JSON buffer exists.
struct BoundedJson(Vec<u8>);
impl std::io::Write for BoundedJson {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > (8 * 1024 * 1024usize).saturating_sub(self.0.len()) {
            return Err(std::io::Error::other(
                "Replay JSON exceeds the 8 MiB download limit",
            ));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn operation_title(operation: &Operation) -> &'static str {
    match operation {
        Operation::ReconstructRunProduct { .. } => "Reconstruct saved material",
        Operation::CatalogDigest { .. } => "Catalog-pinned digest",
        Operation::DeclareSourceTermini { .. } => "Declare original source termini",
        Operation::Pcr { .. } => "PCR",
        Operation::OverlapExtension { .. } => "Reciprocal overlap extension",
        Operation::PcrWithMismatches { .. } => "PCR with mismatches",
        Operation::Digest { .. } => "Digest under cleavage assumptions",
        Operation::DeclareSourceModifications { .. } => "Declare original source modifications",
        Operation::CheckedDigest { .. } => "Checked restriction digest",
        Operation::CheckedTopo { .. } => "TOPO graph with explicit target acceptance",
        Operation::Topo { .. } => "TOPO physical strand reconstruction",
        Operation::ProcessedTopo { .. } => {
            "Processed TOPO source reconstruction and faithful material bridge"
        }
        Operation::TopoActions { .. } => "TOPO ordered source-backed physical actions",
        Operation::MultiSiteRecombineReference { .. } => "Multi-parent reference exchange",
        Operation::RecombineReference { .. } => "Reference-only BP/LR exchange",
        Operation::TreatEnds { .. } => "Treat ends",
        Operation::RepairEnds { .. } => "Repair ends",
        Operation::Ligate { .. } => "Ligate",
        Operation::Tail3Prime { .. } => "3′ tailing",
        Operation::AnnealOligos { .. } => "Anneal oligos",
        Operation::GoldenGate { .. } => "Golden Gate",
        Operation::SelectedGoldenGate { .. } => "Golden Gate with exact candidate keys",
        Operation::HomologyAssembly { .. } => "Physical homology assembly",
        Operation::ResectThreePrime { .. } => "Explicit 3′ source resection",
        Operation::HomologyAnnealing { .. } => {
            "3′ resection and unsealed 5′-tail annealing · incomplete reference"
        }
        Operation::Overlap { .. } | Operation::ProfiledOverlap { .. } => "Overlap assembly",
    }
}
fn reference_only_steps(plan: &DesignPlan) -> Vec<bool> {
    let mut flags = Vec::with_capacity(plan.steps.len());
    for operation in &plan.steps {
        let direct = matches!(
            operation,
            Operation::Topo { .. }
                | Operation::CheckedTopo { .. }
                | Operation::RecombineReference { .. }
                | Operation::MultiSiteRecombineReference { .. }
        );
        let inherited = collagene_core::verification::operation_sources(operation)
            .into_iter()
            .any(|source| matches!(source,Source::Step(i) if flags.get(*i)==Some(&true)));
        flags.push(direct || inherited);
    }
    flags
}

#[component]
fn RunProvenanceOverview(
    record: StoredValue<Result<run_lineage::RunRecord, String>>,
    epoch: RwSignal<u64>,
    opened_epoch: u64,
    open_event: Callback<(Role, EventBinding)>,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let reference_only = StoredValue::new(record.with_value(|r| {
        r.as_ref()
            .map(|r| reference_only_steps(&r.plan))
            .unwrap_or_default()
    }));
    let count = record.with_value(|r| r.as_ref().map_or(0, |r| r.products.len()));
    let pages = count.max(1).div_ceil(10);
    view! {<Show when=move||{count>0}><section aria-label="Recorded assembly provenance">
        <h4>"Assembly provenance"</h4><p class="helper">"Historical source versions → declared operation → retained product version. Select a version to inspect it read-only. The independent target is listed separately below."</p>
        <For each=move||record.with_value(|r|r.as_ref().map(|r|r.products.iter().skip(page.get().min(pages-1)*10).take(10).map(|p|p.step_index).collect::<Vec<_>>()).unwrap_or_default()) key=|index|*index children=move|step_index|view!{<RunProvenanceStep record epoch opened_epoch open_event step_index reference_only/>}/>
        <p class="helper">{move||format!("Provenance steps {}–{} of {count}",page.get().min(pages-1)*10+1,((page.get().min(pages-1)+1)*10).min(count))}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous provenance steps"</button>
        <button class="text-tool" disabled=move||page.get()+1>=pages on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next provenance steps"</button>
    </section></Show>}
}
#[component]
fn RunProvenanceStep(
    record: StoredValue<Result<run_lineage::RunRecord, String>>,
    epoch: RwSignal<u64>,
    opened_epoch: u64,
    open_event: Callback<(Role, EventBinding)>,
    step_index: usize,
    reference_only: StoredValue<Vec<bool>>,
) -> impl IntoView {
    let source_page = RwSignal::new(0usize);
    view! {<article style="border:1px solid var(--border, #d9dfdf);border-radius:8px;padding:10px;margin:8px 0;">
        <div style="display:flex;flex-wrap:wrap;align-items:center;gap:10px;">
            <div style="flex:2 1 220px;min-width:0;"><strong>"Declared sources"</strong>
                {move||record.with_value(|value|{
                    let Some(r)=value.as_ref().ok() else{return ().into_any()};
                    let Some(operation)=r.plan.steps.get(step_index) else{return view!{<p>"Operation unavailable"</p>}.into_any()};
                    let dependencies=collagene_core::verification::operation_sources(operation);let count=dependencies.len();let pages=count.max(1).div_ceil(6);let current=source_page.get().min(pages-1);
                    let links=dependencies.into_iter().skip(current*6).take(6).map(|source|{
                        let resolved=match source {
                            Source::Input(id)=>r.inputs.iter().find(|input|input.input_id==*id).map(|input|(
                                r.plan.inputs.iter().find(|source|source.id==*id).map(|source|format!("{} · input",source.molecule.name)).unwrap_or_else(||format!("Input {id}")),Role::Input{input_id:id.clone()},input.event.clone()
                            )),
                            Source::Step(index)=>r.products.iter().find(|p|p.step_index==*index).map(|p|(
                                format!("Step {} product",index+1),Role::Product{step_index:*index},EventBinding{document_id:p.document_id.clone(),event_id:p.event_id.clone(),revision:p.revision.clone()}
                            )),
                        };
                        match resolved {
                            Some((label,role,binding))=>{
                                let title=format!("{role:?} · exact document {} · event {} · revision {}",binding.document_id,binding.event_id,binding.revision);
                                view!{<button class="text-tool" style="display:block;max-width:100%;overflow-wrap:anywhere;text-align:left;" title=title disabled=move||epoch.get()!=opened_epoch on:click=move|_|open_event.run((role.clone(),binding.clone()))>{label}</button>}.into_any()
                            },
                            None=>view!{<p class="helper">"Exact retained dependency unavailable"</p>}.into_any(),
                        }
                    }).collect_view();
                    view!{<div>{links}<Show when=move||{count>6}><p class="helper">{format!("Source references {}–{} of {count}",current*6+1,((current+1)*6).min(count))}</p><button class="text-tool" disabled=current==0 on:click=move|_|source_page.set(current.saturating_sub(1))>"Previous sources"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|source_page.set((current+1).min(pages-1))>"Next sources"</button></Show></div>}.into_any()
                })}
            </div>
            <span aria-hidden="true">"→"</span>
            <div style="flex:2 1 200px;min-width:0;"><strong>{record.with_value(|r|r.as_ref().ok().and_then(|r|r.plan.steps.get(step_index)).map(|op|format!("Step {} · {}",step_index+1,operation_title(op))).unwrap_or_else(||"Operation unavailable".into()))}</strong>
                <p class="helper">{if reference_only.with_value(|flags|flags.get(step_index).copied().unwrap_or(false)){"Reference-only ancestry · historical"}else{"Declared preparation / assembly · historical"}}</p>
            </div>
            <span aria-hidden="true">"→"</span>
            <div style="flex:2 1 200px;min-width:0;"><strong>"Retained product"</strong>
                {record.with_value(|value|value.as_ref().ok().and_then(|r|r.products.iter().find(|p|p.step_index==step_index)).map(|product|{
                    let binding=EventBinding{document_id:product.document_id.clone(),event_id:product.event_id.clone(),revision:product.revision.clone()};
                    let label=format!("Open step {} product",step_index+1);let identity=format!("{} · {}",product.document_id,product.event_id);
                    view!{<p class="helper" style="overflow-wrap:anywhere;">{identity}</p><button class="text-tool" disabled=move||epoch.get()!=opened_epoch on:click=move|_|open_event.run((Role::Product{step_index},binding.clone()))>{label}</button>}
                }))}
            </div>
        </div>
    </article>}
}

#[component]
pub fn RunDetailsPanel(
    snapshot: Callback<(), Result<Workspace, String>>,
    run_id: String,
    epoch: RwSignal<u64>,
    on_open_event: Callback<EventRequest, Result<(), String>>,
) -> impl IntoView {
    let opened_epoch = epoch.get_untracked();
    let frame = snapshot.run(());
    let record =
        frame.and_then(|frame| run_lineage::load_run(&frame, &run_id).map_err(|e| e.to_string()));
    let record = StoredValue::new(record);
    let run_id = StoredValue::new(run_id);
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let replay = RwSignal::new(None::<run_lineage::ReplayReport>);
    let generation = RwSignal::new(0u64);
    let evidence_worker = StoredValue::new_local(None::<Rc<worker::CachedWorker>>);
    let evidence_pending = RwSignal::new(false);
    let product_page = RwSignal::new(0usize);
    let input_page = RwSignal::new(0usize);
    let open_event = Callback::new(move |(role, binding): (Role, EventBinding)| {
        if epoch.get_untracked() != opened_epoch {
            status.set("The workspace changed; reopen this cloning run before viewing its retained versions.".into());
            return;
        }
        // Bind navigation to the stored role/event, never a matching current head.
        // The host resolves membership and retained source content again before
        // opening its read-only view; navigation does not replay this assessment.
        if let Err(error) = on_open_event.run(EventRequest {
            run_id: run_id.get_value(),
            binding,
            epoch: opened_epoch,
            role,
        }) {
            status.set(error);
        }
    });
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = evidence_worker.get_value() {
            worker.close();
        }
        evidence_worker.set_value(None);
    });
    Effect::new(move || {
        let _ = epoch.get();
        generation.update(|g| *g = g.wrapping_add(1));
        replay.set(None);
        pending.set(false);
        evidence_pending.set(false);
        if let Some(worker) = evidence_worker.get_value() {
            worker.close();
        }
        evidence_worker.set_value(None);
    });
    let replay_alive = alive.clone();
    let replay_run = move |_| {
        if pending.get_untracked() {
            return;
        }
        let frame = match snapshot.run(()) {
            Ok(frame) => frame,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let fingerprint = match run_lineage::workspace_fingerprint(&frame) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let run_id = run_id.get_value();
        let original_epoch = epoch.get_untracked();
        generation.update(|g| *g = g.wrapping_add(1));
        let job_generation = generation.get_untracked();
        replay.set(None);
        pending.set(true);
        status.set("Replaying the retained plan against its exact source events…".into());
        let alive = replay_alive.clone();
        spawn_local(async move {
            let result = worker::analyze_cancellable(
                Request {
                    id: u32::MAX - 23,
                    revision: fingerprint.clone(),
                    job: AnalysisJob::ReplayCloningRun {
                        workspace: frame,
                        run_id,
                    },
                },
                alive.clone(),
            )
            .await;
            if !alive.load(Ordering::Relaxed)
                || epoch.get_untracked() != original_epoch
                || generation.get_untracked() != job_generation
            {
                return;
            }
            let current_workspace = snapshot.run(());
            let unchanged = current_workspace
                .as_ref()
                .ok()
                .and_then(|frame| run_lineage::workspace_fingerprint(frame).ok())
                .is_some_and(|v| v == fingerprint);
            if !unchanged {
                replay.set(None);
                pending.set(false);
                status.set(
                    "Workspace changed during replay. Replay again for current evidence.".into(),
                );
                return;
            }
            match result {
                Ok(response)
                    if response.id == u32::MAX - 23 && response.revision == fingerprint =>
                {
                    match response.result {
                        AnalysisResult::RunReplay(report)
                            if report.source_workspace_fingerprint == fingerprint
                                && current_workspace.as_ref().is_ok_and(|frame| {
                                    report.validate_for_workspace(frame).is_ok()
                                }) =>
                        {
                            replay.set(Some(report));
                            status.set("Fresh replay completed in this session. This checks source-event integrity and deterministic reconstruction; it is not experimental provenance.".into());
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Unexpected or stale replay response.".into()),
                    }
                }
                Ok(_) => status.set("Replay response does not match this workspace.".into()),
                Err(error) => status.set(error),
            }
            pending.set(false);
        });
    };
    let export = move |_| {
        let Some(report) = replay.get_untracked() else {
            return;
        };
        let fresh = snapshot.run(()).and_then(|frame| {
            run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())
        });
        if !fresh.is_ok_and(|v| v == report.source_workspace_fingerprint) {
            replay.set(None);
            status
                .set("Workspace changed. Replay again before downloading current evidence.".into());
            return;
        }
        let mut bytes = BoundedJson(Vec::new());
        let result = serde_json::to_writer_pretty(&mut bytes, &report)
            .map_err(|e| e.to_string())
            .and_then(|()| download("cloning-run-replay.json", "application/json", &bytes.0));
        crate::app::complete_download(result, move |result| {
            if status.try_get_untracked().is_some() {
                if let Err(error) = result {
                    status.set(error);
                }
            }
        });
    };
    let historical_alive = alive.clone();
    let download_history = Callback::new(move |()| {
        if pending.get_untracked() || epoch.get_untracked() != opened_epoch {
            return;
        }
        let frame = match snapshot.run(()) {
            Ok(frame) => frame,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let fingerprint = match run_lineage::workspace_fingerprint(&frame) {
            Ok(value) => value,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let requested_run = run_id.get_value();
        let executor = match worker::CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        generation.update(|g| *g = g.wrapping_add(1));
        let ticket = generation.get_untracked();
        evidence_worker.set_value(Some(executor.clone()));
        evidence_pending.set(true);
        pending.set(true);
        status.set("Loading recorded reference and physical graph evidence. This does not replay or validate the cloning result.".into());
        let alive = historical_alive.clone();
        spawn_local(async move {
            let reply = executor
                .ask(Request {
                    id: u32::MAX - 24,
                    revision: fingerprint.clone(),
                    job: AnalysisJob::LoadProcessedRunEvidence {
                        workspace: frame,
                        run_id: requested_run.clone(),
                    },
                })
                .await;
            executor.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != ticket
                || epoch.get_untracked() != opened_epoch
            {
                return;
            }
            pending.set(false);
            evidence_pending.set(false);
            evidence_worker.set_value(None);
            let current = snapshot.run(()).and_then(|frame| {
                run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())
            });
            if !current.is_ok_and(|value| value == fingerprint)
                || run_id.get_value() != requested_run
            {
                status.set(
                    "Workspace changed; reopen this run before downloading historical evidence."
                        .into(),
                );
                return;
            }
            match reply {
                Ok(reply) if reply.id == u32::MAX - 24 && reply.revision == fingerprint => {
                    match reply.result {
                        AnalysisResult::HistoricalProcessedRunEvidence(value) => {
                            let identity = record.with_value(|record| {
                                record.as_ref().is_ok_and(|record| {
                                    value.run_id == record.run_id
                                        && value.recorded_engine == record.assessment.engine_version
                                })
                            });
                            if !identity || value.run_id != requested_run || value.step_count > 64 {
                                status.set(
                                    "Historical evidence has a different run identity.".into(),
                                );
                                return;
                            }
                            match value.bytes {
                            Some(bytes)if bytes.len()<=32*1024*1024&&value.step_count>0=>{
                                crate::app::download_with_completion("recorded-cloning-evidence.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}match result{
                                    Ok(())=>status.set("Recorded graph evidence downloaded. This historical assessment does not restore physical authority or establish laboratory success.".into()),
                                    Err(e)=>status.set(e),
                                }})
                            },
                            None if value.step_count==0=>status.set("This run has no retained reference or physical graph attachment. Older records are not reconstructed by this download action; use explicit replay for fresh session evidence.".into()),
                            _=>status.set("Historical graph evidence exceeds its bound or has an inconsistent step count.".into()),
                        }
                        }
                        AnalysisResult::Error(e) => status.set(e),
                        _ => status.set("Unexpected historical graph response.".into()),
                    }
                }
                Ok(_) => status
                    .set("Historical graph response has a different workspace binding.".into()),
                Err(e) => status.set(e),
            }
        });
    });
    let cancel_history = Callback::new(move |()| {
        if !evidence_pending.get_untracked() {
            return;
        }
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = evidence_worker.get_value() {
            worker.close();
        }
        evidence_worker.set_value(None);
        evidence_pending.set(false);
        pending.set(false);
        status.set("Historical evidence download cancelled.".into());
    });
    view! {
        <section class="run-details-panel">
            <p class="run-authority">"Recorded history · unverified on opening"</p>
            <p class="helper">"Saved status and engine describe a historical assessment. Fresh replay is session evidence only. Neither an imported record nor a replay proves which physical reagents were used."</p>
            <RunProvenanceOverview record epoch opened_epoch open_event/>
            {move||match record.get_value(){
                Err(error)=>view!{<p role="alert">{error}</p>}.into_any(),
                Ok(record)=>{
                    let count=record.products.len();let input_count=record.inputs.len();
                    let target=record.target.clone();
                    view!{<h3>{record.plan.name.clone()}</h3>
                        <p>{format!("Run {} · recorded {:?} · engine {}",record.run_id,record.assessment.status,record.assessment.engine_version)}</p>
                        <p class="run-target-binding">{format!("Independent target: document {} · event {} · revision {}",record.target.document_id,record.target.event_id,record.target.revision)}</p>
                        <button class="button secondary" disabled=move||epoch.get()!=opened_epoch on:click=move|_|open_event.run((Role::Target,target.clone()))>"Open target version"</button>
                        <h4>"Source bindings"</h4>
                        <div class="run-retained-inputs">{record.inputs.into_iter().skip(input_page.get()*25).take(25).map(|input|{
                            let binding=input.event.clone();
                            let input_id=input.input_id.clone();
                            view!{<div class="run-retained-input"><p>{format!("{}: document {} · event {} · revision {}",input.input_id,input.event.document_id,input.event.event_id,input.event.revision)}</p>
                                <button class="button secondary" disabled=move||epoch.get()!=opened_epoch on:click=move|_|open_event.run((Role::Input{input_id:input_id.clone()},binding.clone()))>"Open source version"</button>
                            </div>}
                        }).collect_view()}</div>
                        <button class="button secondary" disabled=move||input_page.get()==0 on:click=move|_|input_page.update(|p|*p=p.saturating_sub(1))>"Previous retained inputs"</button>
                        <button class="button secondary" disabled=move||{(input_page.get()+1)*25>=input_count} on:click=move|_|input_page.update(|p|*p+=1)>"Next retained inputs"</button>
                        <h4>"Reconstructed product steps"</h4>
                        <div class="run-products">{record.products.into_iter().skip(product_page.get()*25).take(25).map(|product|{
                            let binding=EventBinding{document_id:product.document_id.clone(),event_id:product.event_id.clone(),revision:product.revision.clone()};
                            let step_index=product.step_index;
                            let operation_details=record.plan.steps.get(product.step_index).map(crate::app::describe_step).unwrap_or_default();
                            let operation=record.plan.steps.get(product.step_index).map(operation_title).unwrap_or("Missing operation");
                            view!{<div class="run-product" data-step-index=product.step_index data-document-id=product.document_id.clone() data-event-id=product.event_id.clone()>
                                <strong>{format!("Step {} · {operation}",product.step_index+1)}</strong>
                                <details><summary>"Operation parameters and source references"</summary><pre style="white-space:pre-wrap;overflow-wrap:anywhere">{operation_details}</pre></details>
                                <p>{format!("Document {} · event {} · revision {}",product.document_id,product.event_id,product.revision)}</p>
                                <button class="button secondary" disabled=move||epoch.get()!=opened_epoch on:click=move|_|open_event.run((Role::Product{step_index},binding.clone()))>"Open product version"</button>
                            </div>}
                        }).collect_view()}</div>
                        <button class="button secondary" disabled=move||product_page.get()==0 on:click=move|_|product_page.update(|p|*p=p.saturating_sub(1))>"Previous product steps"</button>
                        <button class="button secondary" disabled=move||{(product_page.get()+1)*25>=count} on:click=move|_|product_page.update(|p|*p+=1)>"Next product steps"</button>
                    }.into_any()
                }
            }}
            <p class="helper">"Recorded reference and physical graphs are historical evidence, separate from fresh replay. Downloading them restores no scientific or physical-material authority."</p>
            <button class="button secondary" disabled=move||pending.get()||epoch.get()!=opened_epoch on:click=move|_|download_history.run(())>"Download recorded cloning evidence JSON"</button>
            <Show when=move||evidence_pending.get()><button class="text-tool" on:click=move|_|cancel_history.run(())>"Cancel historical evidence download"</button></Show>
            <button class="button primary" disabled=move||pending.get() on:click=replay_run>"Replay cloning run"</button>
            <button class="button secondary" disabled=move||replay.get().is_none()||pending.get() on:click=export>"Download run replay JSON"</button>
            {move||replay.get().map(|report|view!{<div class="run-replay-result"><strong>{format!("Fresh replay: {:?}",report.current_report.status)}</strong><p>{format!("Engine matches recorded: {} · products match recorded: {} · status matches recorded: {}",report.engine_matches,report.products_match,report.status_matches)}</p>{report.current_report.findings.into_iter().map(|finding|view!{<p>{finding}</p>}).collect_view()}{report.current_report.products.iter().enumerate().take(25).filter_map(|(i,p)|crate::checked_cloning::modification_summary(p).map(|summary|view!{<details class="checked-modification-summary"><summary>{format!("Fresh step {} · positional modification evidence",i+1)}</summary><p>{summary}</p></details>})).collect_view()}</div>})}
            <p role="status">{move||status.get()}</p>
        </section>
    }
}
