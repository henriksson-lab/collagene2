//! Source-bin reference exchange inside the ordinary assembly intent/review flow.
use crate::{
    assembly_batch_projection,
    assembly_draft_model::Draft,
    cloning_methods_model::{Mode, OverlapPreparation, Treatment},
    multisite_form::{ComponentIntent, Form},
    worker::CachedWorker,
};
use collagene_core::{
    combinatorial_design as design, multisite_recombination as science,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    recombination::reference,
    run_lineage,
    transfer::Candidate,
    verification::{DesignPlan, Operation, PlanInput, Source, TopoSource},
    workspace::Workspace,
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::{JsFuture, spawn_local};

#[component]
pub fn MultisiteSettings(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    busy: RwSignal<bool>,
) -> impl IntoView {
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, Ordering::Relaxed));
    let text = Memo::new(move |_| {
        draft.with(|d| {
            d.methods
                .multisite
                .as_ref()
                .map(|f| f.profile_json.clone())
                .unwrap_or_default()
        })
    });
    let parsed = Memo::new(move |_| -> Result<reference::multisite::Profile, String> {
        let p: reference::multisite::Profile =
            serde_json::from_str(&text.get()).map_err(|e| e.to_string())?;
        p.validate().map_err(|e| e.to_string())?;
        Ok(p)
    });
    let import = move |event| {
        let input = event_target::<web_sys::HtmlInputElement>(&event);
        let Some(file) = input.files().and_then(|f| f.get(0)) else {
            return;
        };
        if file.size() > 512.0 * 1024.0 {
            status.set("Reference profiles must be at most512KiB.".into());
            return;
        }
        let before = Arc::new(draft.get_untracked());
        let saved_epoch = epoch.get_untracked();
        generation.update(|g| *g = g.wrapping_add(1));
        let token = generation.get_untracked();
        let alive = alive.clone();
        spawn_local(async move {
            let value = JsFuture::from(file.text())
                .await
                .map_err(|_| "Could not read profile file.".to_string())
                .and_then(|v| v.as_string().ok_or("Profile text unavailable.".into()))
                .and_then(|s| {
                    if s.len() > 512 * 1024 {
                        return Err("Profile exceeds512KiB.".into());
                    }
                    let p: reference::multisite::Profile =
                        serde_json::from_str(&s).map_err(|e| e.to_string())?;
                    p.validate().map_err(|e| e.to_string())?;
                    Ok(s)
                });
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            if busy.get_untracked()
                || epoch.get_untracked() != saved_epoch
                || !draft.with_untracked(|d| d == before.as_ref())
            {
                status.set("The draft changed while loading; choose the profile again.".into());
                return;
            }
            match value {
                Ok(s) => {
                    draft.update(|d| {
                        let f = d.methods.multisite.get_or_insert_with(Form::default);
                        f.profile_json = s;
                        f.reviewed = false;
                        f.selection = None;
                    });
                    status.set("Profile imported. Review its literal contexts, sources and limitations before accepting.".into());
                }
                Err(e) => status.set(e),
            }
        });
    };
    view! {<section class="assembly-multisite-settings"><h3>"Reference recombination setup"</h3><p class="helper">"Declare BP/LR site contexts for the whole original parents. Your selected source regions define a separate intended reference. This workflow does not digest, amplify, modify termini or establish physical recombinase activity."</p>
        <button class="button secondary" disabled=move||busy.get() on:click=move|_|draft.update(|d|{
            if d.methods.multisite.is_none(){d.methods.multisite=Some(Form::default());}
            for row in &mut d.methods.rows{row.preparation=OverlapPreparation::WholeLinear;row.treatment=Treatment::Keep;row.repair=None;}
            d.methods.allow_nicked=false;
            if let Some(f)=&mut d.methods.multisite{f.selection=None;}
        })>"Use whole original parents without preparation or treatment"</button>
        <p class="helper">"This explicit choice preserves source regions, orientation and saved primer text. It clears the selected physical preparation, end repair and treatment for this reference design."</p>
        <div class="method-row-actions"><button class="button secondary" disabled=move||busy.get() on:click=move|_|{draft.update(|d|d.methods.multisite=Some(Form::default()));status.set("Loaded registered1/2 contexts. Review their specific source scope.".into());}>"Registered1/2 source contexts"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|match Form::three_entry(){Ok(form)=>{draft.update(|d|d.methods.multisite=Some(form));status.set("Loaded exact228/383/302 + VectorII LR contexts. This finite profile is not a universal kit model; p5E-MCS is an MCS cassette, not a promoter.".into());},Err(e)=>status.set(e)}>"Three-entry source-edition LR contexts"</button></div>
        <label>"Reference reaction"<select disabled=move||busy.get() prop:value=move||draft.with(|d|if d.methods.multisite.as_ref().is_some_and(|f|f.reaction==reference::Reaction::BP){"BP"}else{"LR"}) on:change=move|e|draft.update(|d|{let f=d.methods.multisite.get_or_insert_with(Form::default);f.reaction=if event_target_value(&e)=="BP"{reference::Reaction::BP}else{reference::Reaction::LR};f.selection=None;})><option value="LR">"LR · L/R reciprocal reference exchange"</option><option value="BP">"BP · B/P reciprocal reference exchange"</option></select></label>
        <label>"Declared specificity IDs"<input maxlength="128" disabled=move||busy.get() prop:value=move||draft.with(|d|d.methods.multisite.as_ref().map(|f|f.specificities.clone()).unwrap_or_default()) on:input=move|e|draft.update(|d|{let f=d.methods.multisite.get_or_insert_with(Form::default);f.specificities=event_target_value(&e);f.selection=None;})/></label>
        <p class="helper">"The default profile covers1 and2. The optional three-entry preset covers only its exact deposited parent contexts; other sources require reviewed compatible contexts."</p>
        <details><summary>"Review profile sources, exact contexts and limitations"</summary>{move||match parsed.get(){Ok(p)=>view!{<p>{format!("{} · edition {} · {} full contexts",p.id,p.edition,p.variants.len())}</p>{p.limitations.into_iter().map(|s|view!{<p class="helper">{s}</p>}).collect_view()}{p.variants.into_iter().map(|v|view!{<details><summary>{format!("{} · {:?}{} · {} bases · core offset {}",v.id,v.family,v.subtype,v.sequence.len(),v.overlap_start)}</summary><pre>{v.sequence}</pre><p>{format!("Source record {}",v.source_id)}</p></details>}).collect_view()}{p.sources.into_iter().map(|s|view!{<p><a href=s.url target="_blank" rel="noopener noreferrer">{s.vector_name}</a>{format!(" · {}..{} · {} · SHA256 {}",s.start_1based,s.end_1based,if s.reverse_extraction{"reverse extraction"}else{"forward extraction"},s.raw_sha256)}</p>}).collect_view()}}.into_any(),Err(e)=>view!{<p role="alert">{format!("Profile not ready: {e}")}</p>}.into_any()}}</details>
        <label>"Import reviewed full-context profile JSON"<input type="file" accept=".json,application/json" disabled=move||busy.get() on:change=import/></label>
        <label class="checkbox-label"><input type="checkbox" disabled=move||busy.get() prop:checked=move||draft.with(|d|d.methods.multisite.as_ref().is_some_and(|f|f.reviewed)) on:change=move|e|draft.update(|d|{let f=d.methods.multisite.get_or_insert_with(Form::default);f.reviewed=event_target_checked(&e);f.selection=None;})/>"I reviewed this declared profile and its finite context scope"</label>
        {move||draft.with(|d|d.methods.multisite.as_ref().and_then(|f|f.selection.as_ref()).map(|s|view!{<p class="helper">{format!("Saved intent: candidate {}, component {} · source context {}… . Inspect again for fresh evidence; this saved choice is not a restored verification.",s.candidate_id+1,s.component_id+1,s.source_revision.chars().take(12).collect::<String>())}</p>}))}
        <p class="helper">"Limits: at most5 parents,16 specificity pairs,256 candidates,32 components per candidate and32MiB evidence. Search limits and possible sites remain explicit; target matching cannot close incomplete evidence."</p><p role="status">{move||status.get()}</p>
    </section>}
}

struct Evidence {
    before: Arc<Draft>,
    epoch: u64,
    fingerprint: String,
    inputs: Vec<PlanInput>,
    report: design::MultiSiteInspection,
}
#[component]
pub fn MultisiteReview(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    busy: RwSignal<bool>,
    snapshot: Callback<(), Result<Workspace, String>>,
    capture: Callback<(), Result<Vec<PlanInput>, String>>,
    guard: Callback<(), Result<(), String>>,
    target_revision: String,
    on_review: Callback<Candidate>,
    on_finalize: Callback<Candidate, Result<(), String>>,
) -> impl IntoView {
    let target_revision = StoredValue::new(target_revision);
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let cancel = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        pending.set(false);
    });
    Effect::new(move |_| {
        let _ = (draft.get(), epoch.get());
        cancel.run(());
        evidence.set(None);
        page.set(0);
    });
    let current = Callback::new(move |saved: Arc<Evidence>| -> Result<(), String> {
        guard.run(())?;
        if busy.get_untracked()
            || epoch.get_untracked() != saved.epoch
            || !draft.with_untracked(|d| d == saved.before.as_ref())
        {
            return Err(
                "Source, profile, preparation or intended reference changed; inspect again.".into(),
            );
        }
        if run_lineage::workspace_fingerprint(&snapshot.run(())?).map_err(|e| e.to_string())?
            != saved.fingerprint
        {
            return Err("Workspace changed; inspect the reference exchange again.".into());
        }
        Ok(())
    });
    let inspect_alive = alive.clone();
    let inspect = move |_| {
        if pending.get_untracked() || busy.get_untracked() {
            return;
        }
        let prepared = (|| -> Result<_, String> {
            guard.run(())?;
            let before = Arc::new(draft.get_untracked());
            if before.methods.mode != Mode::MultiSiteReference {
                return Err("Choose multi-parent reference exchange first.".into());
            }
            let inputs = capture.run(())?;
            let (catalog, construct) = assembly_batch_projection::active(&before)?;
            let frame = snapshot.run(())?;
            let fingerprint =
                run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
            let revision = design::inspect_multisite_revision(&catalog, &construct, &frame)
                .map_err(|e| e.to_string())?;
            let w = Rc::new(CachedWorker::new()?);
            Ok((
                before,
                inputs,
                fingerprint,
                revision,
                w,
                AnalysisJob::InspectCombinatorialMultiSite {
                    workspace: frame,
                    design: catalog,
                    construct,
                },
            ))
        })();
        let (before, inputs, fingerprint, revision, w, job) = match prepared {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        cancel.run(());
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        worker.set_value(Some(w.clone()));
        pending.set(true);
        evidence.set(None);
        let saved_epoch = epoch.get_untracked();
        let alive = inspect_alive.clone();
        status.set("Reconstructing declared reciprocal reference components against the already frozen independent intent…".into());
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let result = (|| -> Result<Arc<Evidence>, String> {
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Reference response binding changed.".into());
                }
                let report = match response.result {
                    AnalysisResult::CombinatorialMultiSiteInspection(report) => report,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected reference inspection response.".into()),
                };
                report
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                let saved = Arc::new(Evidence {
                    before,
                    epoch: saved_epoch,
                    fingerprint,
                    inputs,
                    report,
                });
                current.run(saved.clone())?;
                Ok(saved)
            })();
            match result {
                Ok(saved) => {
                    status.set("All retained reciprocal components are available. Explicitly choose the component; the independent target never chooses for you.".into());
                    evidence.set(Some(saved));
                }
                Err(e) => status.set(e),
            }
        });
    };
    let apply = Callback::new(
        move |(saved, selection, finalize): (Arc<Evidence>, science::Selection, bool)| {
            let result = (|| -> Result<(), String> {
                current.run(saved.clone())?;
                let document = saved
                    .report
                    .batch
                    .documents
                    .iter()
                    .find(|d| d.selection == selection)
                    .ok_or("Selected reciprocal component is unavailable.")?;
                let mut next = saved.before.as_ref().clone();
                let form = next
                    .methods
                    .multisite
                    .as_mut()
                    .ok_or("Reference profile unavailable.")?;
                let request = form.request()?;
                form.selection = Some(ComponentIntent {
                    batch_revision: selection.batch_revision.clone(),
                    candidate_id: selection.candidate_id,
                    component_id: selection.component_id,
                    source_revision: saved.report.source_revision.clone(),
                });
                let roles = next.methods.role_ids(saved.inputs.len())?;
                if document.selection.batch_revision != saved.report.batch.batch_revision {
                    return Err("Component source binding differs.".into());
                }
                let plan = DesignPlan {
                    name: format!(
                        "{} · {:?} reference exchange",
                        next.name, request.reaction.reaction
                    ),
                    target_revision: target_revision.get_value(),
                    inputs: saved.inputs.clone(),
                    steps: vec![Operation::MultiSiteRecombineReference {
                        sources: roles
                            .iter()
                            .map(|r| TopoSource {
                                role: r.role.clone(),
                                source: Source::Input(r.input_id.clone()),
                            })
                            .collect(),
                        request,
                        selection,
                    }],
                };
                next = if next.alternatives.is_some() {
                    crate::assembly_alternatives::save_active(&next)?
                } else {
                    next.validate()?;
                    next
                };
                draft.set(next);
                let candidate = Candidate::new(plan);
                if finalize {
                    on_finalize.run(candidate)?;
                } else {
                    on_review.run(candidate);
                }
                Ok(())
            })();
            match result{Ok(())=>status.set("Explicit component intent saved and passed to the existing independent-target plan review. Reference-only incompleteness remains visible.".into()),Err(e)=>status.set(e)}
        },
    );
    view! {<section class="assembly-multisite-review"><h3>"Inspect whole-parent reference exchange"</h3><p class="helper">"The intended reference above was frozen from your original source-region choices before this computation. Both desired and reciprocal components remain inspectable. A saved component choice is intent, not verified chemistry."</p><button class="button primary" disabled=move||busy.get()||pending.get() on:click=inspect>"Inspect reciprocal reference components"</button><button class="text-tool" disabled=move||!pending.get() on:click=move|_|cancel.run(())>"Cancel reference inspection"</button><p role="status">{move||status.get()}</p>
        {move||evidence.get().map(|saved|{
            let batch=&saved.report.batch;let count=batch.documents.len();let pages=count.max(1).div_ceil(8);let start=page.get().min(pages-1)*8;
            let notes=batch.reference.warnings.iter().chain(&batch.reference.limits).chain(&batch.reference.unsupported).chain(&batch.warnings).map(|s|view!{<p class="helper">{s.clone()}</p>}).collect_view();
            let docs=batch.documents.iter().skip(start).take(8).map(|d|{let selected=d.selection.clone();let saved=saved.clone();let for_finalize=saved.clone();let finalized=selected.clone();let preview=d.molecule.sequence.as_str().chars().take(240).collect::<String>();view!{<article class="candidate-review"><h4>{format!("Candidate {} · component {}",selected.candidate_id+1,selected.component_id+1)}</h4><p>{format!("{} bases · {:?} · {} reciprocal components retained · {} mapped features",d.molecule.sequence.len(),d.molecule.topology,d.reciprocal_components.len(),d.molecule.features.len())}</p><pre>{format!("{}{}",preview,if d.molecule.sequence.len()>240{"… (preview only)"}else{""})}</pre><details><summary>"Exact identity and reciprocal component inventory"</summary><p>{format!("Reference revision {}",d.reference_product_revision)}</p>{d.reciprocal_components.iter().map(|c|view!{<p>{format!("Component {} · {} bases · {:?} · {}",c.component_id+1,c.length,c.topology,c.reference_revision)}</p>}).collect_view()}{d.warnings.iter().map(|s|view!{<p>{s.clone()}</p>}).collect_view()}</details><button class="button primary" disabled=move||busy.get() on:click=move|_|apply.run((saved.clone(),selected.clone(),false))>"Use this component and review plan"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|apply.run((for_finalize.clone(),finalized.clone(),true))>"Use this component and prepare outputs"</button></article>}}).collect_view();
            view!{<p>{format!("{} candidates · {} reciprocal documents · declared search {}",batch.reference.candidates.len(),count,if batch.reference.complete{"complete within its scope"}else{"incomplete"})}</p>{notes}{docs}<div class="method-row-actions"><button class="text-tool" disabled=page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous components"</button><span>{format!("Page {} of {}",page.get()+1,pages)}</span><button class="text-tool" disabled=page.get()+1>=pages on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next components"</button></div>}
        })}
    </section>}
}
