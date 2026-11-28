//! A construct's finite reaction choice is editable intent, never restored proof.
use crate::assembly::draft::assembly_alternatives;
use crate::assembly::draft::model::Draft;
use crate::assembly::library::assembly_batch_projection;
use crate::cloning::preparation::golden_gate_pool_form::CandidateIntent;
use crate::cloning::preparation::golden_gate_pool_form::PoolForm;
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::golden_gate::reaction,
    documents::workspace::Workspace,
    material::lineage as run_lineage,
    planning::combinatorial::{self as design, Recipe, batch::RowRequest},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

struct Evidence {
    draft: Draft,
    epoch: u64,
    fingerprint: String,
    source_revision: String,
    report: reaction::Report,
    pool: Option<Arc<collagene_core::cloning::golden_gate::GoldenGateReport>>,
}

#[component]
pub fn AssemblyPoolView(
    draft: RwSignal<Draft>,
    snapshot: Callback<(), Result<Workspace, String>>,
    epoch: RwSignal<u64>,
    busy: RwSignal<bool>,
    on_review: Callback<
        crate::assembly::draft::assembly_workspace::ReviewRequest,
        Result<(), String>,
    >,
    on_finalize: Callback<
        crate::assembly::draft::assembly_workspace::ReviewRequest,
        Result<(), String>,
    >,
) -> impl IntoView {
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let note_page = RwSignal::new(0usize);
    let probe_page = RwSignal::new(0usize);
    let inspected = RwSignal::new(None::<usize>);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let cancel = Callback::new(move |()| {
        generation.update(|v| *v = v.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
    });
    Effect::new(move |_| {
        let _ = (draft.get(), epoch.get());
        cancel.run(());
        evidence.set(None);
        inspected.set(None);
        page.set(0);
        note_page.set(0);
        probe_page.set(0);
    });
    let guard = Callback::new(move |saved: Arc<Evidence>| -> Result<(), String> {
        if busy.get_untracked()
            || epoch.get_untracked() != saved.epoch
            || !draft.with_untracked(|current| current == &saved.draft)
        {
            return Err(
                "The construct, source or options changed; inspect the reaction again.".into(),
            );
        }
        if run_lineage::workspace_fingerprint(&snapshot.run(())?).map_err(|e| e.to_string())?
            != saved.fingerprint
        {
            return Err("The workspace changed; inspect the construct's reaction again.".into());
        }
        Ok(())
    });
    let enable = move |_| {
        let result = (|| -> Result<Draft, String> {
            if busy.get_untracked() {
                return Err("Finish the current workspace operation first.".into());
            }
            let mut next = draft.get_untracked();
            assembly_alternatives::require_active_choice_match(&next)?;
            if next.alternatives.is_none()
                || next.methods.mode
                    != crate::cloning::preparation::cloning_methods_model::Mode::GoldenGate
            {
                return Err("Open a Golden Gate construct from the source matrix first.".into());
            }
            if next.methods.golden_gate_pool.is_none() {
                next.methods.golden_gate_pool = Some(PoolForm {
                    choices: vec![0; next.methods.rows.len()],
                    circular: next.circular,
                    ..Default::default()
                });
            }
            assembly_alternatives::save_active(&next)
        })();
        match result {
            Ok(next) => {
                draft.set(next);
                status.set("Finite pool settings enabled for this construct. Choose substrates and search limits in its assembly controls, then inspect here.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let inspect_alive = alive.clone();
    let inspect = move |_| {
        if pending.get_untracked() || busy.get_untracked() {
            return;
        }
        let prepared = (|| -> Result<_, String> {
            let original = draft.get_untracked();
            let active = original
                .alternatives
                .as_ref()
                .ok_or("Create a source matrix and open a construct first.")?
                .active_row
                .clone();
            let projection = assembly_batch_projection::project(&original)?;
            let catalog = projection
                .request
                .design
                .as_ref()
                .ok_or("Capture source alternatives first.")?;
            let row = projection
                .request
                .rows
                .iter()
                .find_map(|row| match row {
                    RowRequest::Ready(row) if row.id == active => Some(Ok(row)),
                    RowRequest::DraftError { id, message, .. } if id == &active => {
                        Some(Err(message.clone()))
                    }
                    _ => None,
                })
                .ok_or("The active construct is not in the matrix.")??;
            let frame = snapshot.run(())?;
            let fingerprint =
                run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
            let resolved = design::resolve_pool(catalog, row, &frame).map_err(|e| e.to_string())?;
            let Recipe::GoldenGatePool {
                request,
                source_revision,
                ..
            } = resolved.recipe
            else {
                return Err("Enable finite pool settings for this construct first.".into());
            };
            let revision = reaction::inspect_revision(&resolved.inputs, &resolved.target, &request)
                .map_err(|e| e.to_string())?;
            let worker = Rc::new(CachedWorker::new()?);
            Ok((
                original,
                epoch.get_untracked(),
                fingerprint,
                source_revision,
                revision,
                worker,
                AnalysisJob::InspectGoldenGateReaction {
                    inputs: resolved.inputs,
                    target: resolved.target,
                    request,
                },
            ))
        })();
        let (original, saved_epoch, fingerprint, source_revision, revision, w, job) = match prepared
        {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        cancel.run(());
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        worker.set_value(Some(w.clone()));
        pending.set(true);
        evidence.set(None);
        inspected.set(None);
        status.set(
            "Reconstructing this construct's declared substrates and finite pool in Rust…".into(),
        );
        let alive = inspect_alive.clone();
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
                    return Err("Reaction response differs from its request.".into());
                }
                let report = match response.result {
                    AnalysisResult::GoldenGateReaction(report) => report,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected construct reaction response.".into()),
                };
                report
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                let pool = report.pool.as_ref().map(|pool| Arc::new(pool.clone()));
                let bound = Arc::new(Evidence {
                    draft: original,
                    epoch: saved_epoch,
                    fingerprint,
                    source_revision,
                    report,
                    pool,
                });
                guard.run(bound.clone())?;
                Ok(bound)
            })();
            match result {
                Ok(bound) => {
                    let withheld = bound.pool.is_none();
                    evidence.set(Some(bound));
                    status.set(if withheld{"Pool withheld: inspect preparation and probe findings. No candidate is available to choose."}else{"Review alternatives explicitly. Candidate zero and a target-matching candidate are never selected automatically."}.into());
                }
                Err(error) => status.set(error),
            }
        });
    };
    let choose = Callback::new(move |index: usize| {
        let result = (|| -> Result<Draft, String> {
            let saved = evidence
                .get_untracked()
                .ok_or("Inspect the current construct first.")?;
            guard.run(saved.clone())?;
            let candidate = saved
                .report
                .pool
                .as_ref()
                .and_then(|pool| pool.candidates.get(index))
                .ok_or("The candidate is unavailable.")?;
            if !candidate.restriction_site_free {
                return Err(
                    "A candidate retaining or recreating recognition sites is diagnostic only."
                        .into(),
                );
            }
            let intent = CandidateIntent {
                context_revision: saved.report.context_revision.clone(),
                source_revision: saved.source_revision.clone(),
                candidate_index: index,
                directed_sequence_key: candidate.directed_sequence_key.clone(),
                physical_sequence_key: candidate.physical_sequence_key.clone(),
            };
            assembly_alternatives::set_active_pool_candidate(&draft.get_untracked(), Some(intent))
        })();
        match result {
            Ok(next) => {
                draft.set(next);
                status.set("Candidate intent saved for this construct. Batch evaluation and output creation must reconstruct and check its exact sequence and physical keys afresh.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let clear = move |_| {
        if busy.get_untracked() {
            return;
        }
        match assembly_alternatives::set_active_pool_candidate(&draft.get_untracked(), None) {
            Ok(next) => {
                draft.set(next);
                status.set(
                    "Candidate intent cleared. Inspect and choose before batch evaluation.".into(),
                );
            }
            Err(error) => status.set(error),
        }
    };
    let ordered = move |_| {
        if busy.get_untracked() {
            return;
        }
        let result = (|| -> Result<Draft, String> {
            let mut next = draft.get_untracked();
            assembly_alternatives::require_active_choice_match(&next)?;
            next.methods.golden_gate_pool = None;
            assembly_alternatives::save_active(&next)
        })();
        match result {
            Ok(next) => {
                draft.set(next);
                status.set("Returned to ordered prepared assembly. Previous pool options and candidate intent were cleared explicitly.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let download = move |_| {
        let Some(saved) = evidence.get_untracked() else {
            status.set("Inspect the construct first.".into());
            return;
        };
        let result = (|| {
            guard.run(saved.clone())?;
            #[derive(serde::Serialize)]
            struct Export<'a> {
                source_revision: &'a str,
                report: &'a reaction::Report,
            }
            struct Bounded(Vec<u8>);
            impl std::io::Write for Bounded {
                fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                    if self
                        .0
                        .len()
                        .checked_add(bytes.len())
                        .is_none_or(|n| n > 32 * 1024 * 1024)
                    {
                        return Err(std::io::Error::other(
                            "Reaction evidence exceeds32MiB JSON export limit.",
                        ));
                    }
                    self.0.extend_from_slice(bytes);
                    Ok(bytes.len())
                }
                fn flush(&mut self) -> std::io::Result<()> {
                    Ok(())
                }
            }
            let mut bytes = Bounded(Vec::new());
            serde_json::to_writer(
                &mut bytes,
                &Export {
                    source_revision: &saved.source_revision,
                    report: &saved.report,
                },
            )
            .map_err(|e| e.to_string())?;
            guard.run(saved.clone())?;
            crate::io::download::ExportArtifact::prepare(
                "construct-reaction.json",
                "application/json",
                &bytes.0,
            )
        })();
        match result {
            Ok(artifact) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    guard.run(saved.clone())
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() { return; }
                    status.set(match result { Ok(()) => "Complete pinned reaction discovery saved. Imported evidence does not restore scientific authority.".into(), Err(error) => error });
                });
            }
            Err(error) => status.set(error),
        }
    };
    let current = Memo::new(move |_| {
        evidence.with(|saved| {
            saved.as_ref().is_some_and(|saved| {
                saved.epoch == epoch.get() && draft.with(|d| d == &saved.draft)
            })
        })
    });
    let handoff = Callback::new(move |finalize: bool| {
        let result = (|| -> Result<(), String> {
            if busy.get_untracked() {
                return Err("Finish the current workspace operation first.".into());
            }
            let original = draft.get_untracked();
            let active = original
                .alternatives
                .as_ref()
                .ok_or("Open a source-matrix construct first.")?
                .active_row
                .clone();
            let projected = assembly_batch_projection::project(&original)?;
            let catalog = projected
                .request
                .design
                .as_ref()
                .ok_or("Capture source alternatives first.")?;
            let construct = projected
                .request
                .rows
                .iter()
                .find_map(|row| match row {
                    RowRequest::Ready(row) if row.id == active => Some(Ok(row)),
                    RowRequest::DraftError { id, message, .. } if id == &active => {
                        Some(Err(message.clone()))
                    }
                    _ => None,
                })
                .ok_or("The active construct is unavailable.")??;
            let resolved = design::resolve_pool(catalog, construct, &snapshot.run(())?)
                .map_err(|e| e.to_string())?;
            let plan = resolved.compile().map_err(|e| e.to_string())?;
            let intent = match resolved
                .intent
                .ok_or("An independent original-source target is required.")?
            {
                design::IndependentIntent::Regions(intent) => {
                    crate::assembly::review::assembly_review::IntendedReference::Regions(intent)
                }
                design::IndependentIntent::Junctions(intent) => {
                    crate::assembly::review::assembly_review::IntendedReference::Junctions(intent)
                }
                design::IndependentIntent::Substitutions(intent) => {
                    crate::assembly::review::assembly_review::IntendedReference::Substitutions(
                        intent,
                    )
                }
                design::IndependentIntent::Edits(intent) => {
                    crate::assembly::review::assembly_review::IntendedReference::Edits(intent)
                }
            };
            let request = crate::assembly::draft::assembly_workspace::ReviewRequest {
                candidate: collagene_core::documents::transfer::Candidate::new(plan),
                sources: resolved
                    .bindings
                    .into_iter()
                    .map(|binding| binding.event)
                    .collect(),
                intent,
            };
            if finalize {
                on_finalize.run(request)
            } else {
                on_review.run(request)
            }
        })();
        status.set(match result{Ok(())=>"The exact selected reaction was submitted for fresh native reconstruction and independent-target review.".into(),Err(error)=>error});
    });
    view!{<details class="assembly-pool" aria-label="Construct reaction candidate"><summary>"Construct reaction candidate"</summary>
        <p class="helper">"Use this for a source-matrix construct with prepared fragments, whole originals or complete PCR amplicons. The finite one-copy pool includes competing digest pieces and circular carryover. A saved choice is intent; only fresh native reconstruction can check it."</p>
        <button class="text-tool" disabled=move||busy.get()||pending.get() on:click=enable>"Enable finite pool for this construct"</button>
        <button class="button secondary" disabled=move||busy.get()||pending.get() on:click=inspect>"Inspect this construct's reaction"</button>
        <button class="text-tool" disabled=move||busy.get()||pending.get() on:click=ordered>"Use ordered prepared assembly instead"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{cancel.run(());status.set("Reaction inspection cancelled.".into());}>"Cancel inspection"</button>
        <button class="button secondary" disabled=move||busy.get()||pending.get() on:click=move|_|handoff.run(false)>"Review selected reaction and junctions"</button>
        <button class="button primary" disabled=move||busy.get()||pending.get() on:click=move|_|handoff.run(true)>"Create outputs from selected reaction"</button>
        {move||draft.with(|d|d.methods.golden_gate_pool.as_ref().and_then(|pool|pool.candidate.clone())).map(|pin|view!{<p>{format!("Declared candidate {} · sequence {} · physical {}",pin.candidate_index+1,pin.directed_sequence_key,pin.physical_sequence_key)}</p><button class="text-tool" disabled=move||busy.get() on:click=clear>"Clear declared candidate"</button>})}
        {move||evidence.get().map(|saved|{
            let notes=&saved.report.warnings;let pages=notes.len().div_ceil(40).max(1);let current_page=note_page.get().min(pages-1);
            let warnings=notes.iter().skip(current_page*40).take(40).map(|note|view!{<li>{note.clone()}</li>}).collect_view();
            let probe_count=saved.report.verification.assumptions.len();
            let probe_pages=probe_count.div_ceil(40).max(1);let probe_current=probe_page.get().min(probe_pages-1);
            let probe_findings=crate::protocol::verification_view::report(&saved.report.verification);
            let probe_notes=saved.report.verification.assumptions.iter().map(|assumption|format!("Assumption: {assumption}"))
                .skip(probe_current*40).take(40).map(|note|view!{<li>{note}</li>}).collect_view();
            let pool=saved.pool.as_ref().map(|pool|{
                let pages=pool.candidates.len().div_ceil(20).max(1);let current_page=page.get().min(pages-1);
                let candidates=pool.candidates.iter().enumerate().skip(current_page*20).take(20).map(|(index,candidate)|{
                    let target=saved.report.assessments.iter().find(|a|a.candidate_index==index).map(|a|a.target_matches);
                    let allowed=candidate.restriction_site_free;
                    view!{<article><h4>{format!("Candidate {} · {} bases · {} fragments",index+1,candidate.duplex.top.len(),candidate.path.len())}</h4><p>{format!("Independent target match {target:?} · restriction-site-free {allowed}")}</p>
                        <button class="text-tool" disabled=move||!current.get() on:click=move|_|inspected.set(Some(index))>"Inspect top-strand reference"</button>
                        <button class="button secondary" disabled=move||busy.get()||!current.get()||!allowed on:click=move|_|choose.run(index)>"Use this candidate for this construct"</button>
                    </article>}
                }).collect_view();
                view!{<section><p>{format!("{} candidates · {} states · finite search completed {} · circular carryover roles {:?}",pool.candidates.len(),pool.states_examined,pool.enumeration_complete,pool.carryover_sources)}</p>
                    <crate::cloning::preparation::golden_gate_candidate_view::PoolSearchNotes pool=pool.clone()/>
                    {candidates}<p>{format!("Candidate page {} of {pages}",current_page+1)}</p><button class="text-tool" disabled=current_page==0 on:click=move|_|page.set(current_page.saturating_sub(1))>"Previous candidates"</button><button class="text-tool" disabled={current_page+1>=pages} on:click=move|_|page.set(current_page+1)>"Next candidates"</button>
                </section>}
            });
            view!{<section><ul>{warnings}</ul><p>{format!("Reaction note page {} of {pages}",current_page+1)}</p><button class="text-tool" disabled=current_page==0 on:click=move|_|note_page.set(current_page.saturating_sub(1))>"Previous reaction notes"</button><button class="text-tool" disabled={current_page+1>=pages} on:click=move|_|note_page.set(current_page+1)>"Next reaction notes"</button>
                <p>{format!("Diagnostic probe assessment: {:?}. This is not the selected construct's verification result.",saved.report.verification.status)}</p>
                {saved.pool.is_none().then(||view!{<p role="alert">"Finite pool withheld. Review failed or incomplete preparation steps; no winner or exhaustive absence is inferred."</p>})}
                <details><summary>"Fresh preparation and probe checks / assumptions"</summary>{probe_findings}<ul>{probe_notes}</ul><p>{format!("Assumptions page {} of {probe_pages}",probe_current+1)}</p><button class="text-tool" disabled=probe_current==0 on:click=move|_|probe_page.set(probe_current.saturating_sub(1))>"Previous probe notes"</button><button class="text-tool" disabled={probe_current+1>=probe_pages} on:click=move|_|probe_page.set(probe_current+1)>"Next probe notes"</button></details>{pool}
            </section>}
        })}
        {move||inspected.get().and_then(|index|evidence.get().and_then(|saved|saved.report.pool.as_ref().and_then(|pool|pool.candidates.get(index)).map(|candidate|(saved.clone(),index,candidate.duplex.top.as_str().to_string(),candidate.duplex.topology)))).map(|(saved,index,dna,topology)|{
            match collagene_model::Molecule::dna("Candidate top-strand reference",&dna,topology){
                Ok(product)=>view!{<p class="helper">"Read-only top-strand reference from this pinned discovery. Physical-state identity is retained separately in the candidate keys; this view does not assert laboratory success."</p><crate::cloning::preparation::golden_gate_candidate_view::GoldenGateCandidateView pool=saved.pool.as_ref().expect("selected candidate pool retained").clone() index/><crate::assembly::review::assembly_product_view::ReferenceCanvas product identity=format!("{}:{index}",saved.report.context_revision) current on_message=Callback::new(move|message|status.set(message))/>}.into_any(),
                Err(error)=>view!{<p>{error.to_string()}</p>}.into_any(),
            }
        })}
        <button class="text-tool" disabled=move||busy.get()||evidence.get().is_none() on:click=download>"Download complete reaction evidence"</button>
        <p role="status">{move||status.get()}</p>
    </details>}.into_any()
}
