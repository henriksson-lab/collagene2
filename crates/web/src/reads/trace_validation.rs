//! Explicit revision-bound read/reference evidence; no imported or AI-asserted verdict.
use crate::platform::worker::CachedWorker;
use collagene_core::{
    material::lineage::EventBinding,
    protocol::{
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    review::reads::trace_validation::{
        MissingQualityPolicy, SiteState, TraceValidationOptions, TraceValidationReport, TrimMode,
        TrimOptions,
    },
};
use collagene_model::{Alphabet, Molecule};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

fn number<T: std::str::FromStr>(text: String, label: &str) -> Result<T, String> {
    text.parse()
        .map_err(|_| format!("{label} is not a valid number"))
}
#[derive(Clone)]
pub struct ReferenceNavigation {
    pub reference: EventBinding,
    pub selection: Option<(usize, usize)>,
    pub epoch: u64,
}
#[derive(Clone)]
pub struct RecordRequest {
    pub reference: EventBinding,
    pub reads: Vec<EventBinding>,
    pub options: TraceValidationOptions,
    pub epoch: u64,
}
struct BoundEvidence {
    report: TraceValidationReport,
    reference: EventBinding,
    reads: Vec<EventBinding>,
    options: TraceValidationOptions,
    epoch: u64,
    indices: Vec<usize>,
    generation: u32,
}
fn consensus_document(
    source: &Molecule,
    report: &TraceValidationReport,
) -> Result<Molecule, String> {
    if source.revision() != report.reference_revision {
        return Err("Reference changed; run validation again before opening a consensus.".into());
    }
    let mut document = Molecule::dna(
        format!("{}_read_consensus", source.name),
        &report.consensus,
        source.topology,
    )
    .map_err(|e| e.to_string())?;
    document.metadata.insert(
        "read_validation:reference_revision".into(),
        report.reference_revision.clone(),
    );
    document.metadata.insert(
        "read_validation:source_revisions".into(),
        serde_json::to_string(&report.read_revisions).map_err(|e| e.to_string())?,
    );
    document.metadata.insert(
        "read_validation:report".into(),
        serde_json::to_string(report).map_err(|e| e.to_string())?,
    );
    document.metadata.insert("duplex_model".into(),"Support-only reference DNA; physical strand/end chemistry is unavailable. Reconstruct the original source graph for physical operations.".into());
    document.metadata.insert("read_validation:note".into(),"Support-only reference-coordinate consensus. N denotes uncovered, low/unknown confidence, unresolved placement, deletion or insufficient support; IUPAC can indicate conflicting eligible calls. Reference annotations are not copied. Retain original read documents/attachments and the independent reference; this is not experimental sample confirmation or physical chemistry.".into());
    collagene_core::protocol::save_document(&document).map_err(|e| e.to_string())?;
    Ok(document)
}

#[component]
pub fn TraceValidationPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    on_consensus: Callback<Molecule, Result<(), String>>,
    source_info: Callback<usize, Option<crate::sources::identity::SourceInfo>>,
    reference_info: Memo<crate::sources::identity::SourceInfo>,
    epoch: RwSignal<u64>,
    on_read: Callback<crate::reads::contigs::ReadNavigation>,
    on_record: Callback<RecordRequest, Result<(), String>>,
    on_reference: Callback<ReferenceNavigation, Result<(), String>>,
) -> impl IntoView {
    let _ = (selection, selected_feature); // Host applies guarded navigation after switching documents.
    let selected = RwSignal::new(Vec::<usize>::new());
    let mode = RwSignal::new("quality".to_string());
    let threshold = RwSignal::new("20".to_string());
    let window = RwSignal::new("5".to_string());
    let minimum_passing = RwSignal::new("4".to_string());
    let minimum_retained = RwSignal::new("20".to_string());
    let missing = RwSignal::new(false);
    let both = RwSignal::new(true);
    let support = RwSignal::new("1".to_string());
    let both_support = RwSignal::new(false);
    let matrix = RwSignal::new("64000000".to_string());
    let total = RwSignal::new("128000000".to_string());
    let options = Memo::new(move |_| {
        Ok::<_, String>(TraceValidationOptions {
            trim: TrimOptions {
                mode: if mode.get() == "all" {
                    TrimMode::AllCalls
                } else {
                    TrimMode::QualityWindows
                },
                raw_score_threshold: number(threshold.get(), "Confidence threshold")?,
                window_size: number(window.get(), "Window length")?,
                minimum_passing: number(minimum_passing.get(), "Passing calls per window")?,
                min_retained_bases: number(minimum_retained.get(), "Minimum retained length")?,
                missing_quality: if missing.get() {
                    MissingQualityPolicy::KeepUnscored
                } else {
                    MissingQualityPolicy::Reject
                },
            },
            both_strands: both.get(),
            min_support: number(support.get(), "Minimum source files")?,
            require_both_strands: both_support.get(),
            max_cells: number(matrix.get(), "Per-alignment cell budget")?,
            max_total_cells: number(total.get(), "Total alignment/proof cell budget")?,
            ..Default::default()
        })
    });
    let report = RwSignal::new(None::<Arc<BoundEvidence>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let read_context = Memo::new(move |_| {
        let _ = epoch.get();
        selected
            .get()
            .into_iter()
            .map(|index| source_info.run(index).map(|info| info.binding))
            .collect::<Option<Vec<_>>>()
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = active.get_value() {
                worker.close();
            }
        });
    }
    Effect::new(move |_| {
        let _ = (
            selected.get(),
            read_context.get(),
            reference_info.get(),
            options.get(),
            epoch.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        report.set(None);
        busy.set(false);
        status.set(String::new());
    });
    let guard = Callback::new(move |bound: Arc<BoundEvidence>| -> Result<(), String> {
        if epoch.get_untracked() != bound.epoch
            || generation.get_untracked() != bound.generation
            || selected.get_untracked() != bound.indices
            || read_context.get_untracked() != Some(bound.reads.clone())
            || reference_info.get_untracked().binding != bound.reference
            || options.get_untracked()? != bound.options
        {
            return Err("The reference, read events or settings changed. Run validation again before using this evidence.".into());
        }
        Ok(())
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let input = selected.get_untracked();
                if input.is_empty() || input.len() > 32 {
                    return Err("Choose 1–32 read documents.".into());
                }
                let reference = reference_info.get_untracked().binding;
                let bindings = read_context
                    .get_untracked()
                    .ok_or("A selected source event is no longer available.")?;
                let source = molecule.get_untracked();
                let reads = input
                    .iter()
                    .map(|&index| {
                        get_query
                            .run(index)
                            .ok_or("A selected read document is no longer open.")
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let mut borrowed = Vec::with_capacity(reads.len() + 1);
                borrowed.push(&source);
                borrowed.extend(reads.iter());
                collagene_model::bounds::model_preflight(&borrowed).map_err(|e| e.to_string())?;
                let revision = source.revision();
                let revisions = reads.iter().map(Molecule::revision).collect::<Vec<_>>();
                if reference.revision != revision
                    || bindings
                        .iter()
                        .map(|b| b.revision.clone())
                        .collect::<Vec<_>>()
                        != revisions
                {
                    return Err("Selected document/event identities do not match the original read/reference models.".into());
                }
                let opts = options.get_untracked()?;
                let worker = Rc::new(CachedWorker::new()?);
                Ok((
                    input, reference, bindings, source, reads, revision, revisions, opts, worker,
                ))
            })();
            let (input, reference, bindings, source, reads, revision, revisions, opts, w) =
                match prepared {
                    Ok(value) => value,
                    Err(error) => {
                        status.set(error);
                        return;
                    }
                };
            let captured_epoch = epoch.get_untracked();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            if let Some(old) = active.get_value() {
                old.close();
            }
            active.set_value(Some(w.clone()));
            busy.set(true);
            report.set(None);
            status.set(
                "Aligning event-bound reads and reconstructing quality/coverage evidence…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::ValidateTraces {
                            molecule: source,
                            reads,
                            options: opts.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || epoch.get_untracked() != captured_epoch
                    || selected.get_untracked() != input
                    || read_context.get_untracked() != Some(bindings.clone())
                    || reference_info.get_untracked().binding != reference
                    || options.get_untracked().as_ref() != Ok(&opts)
                {
                    return;
                }
                active.set_value(None);
                busy.set(false);
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result: AnalysisResult::TraceValidation(computed),
                    }) if id == token
                        && returned == revision
                        && computed.reference_revision == revision
                        && computed.read_revisions == revisions
                        && computed.options == opts
                        && computed.reads.len() == bindings.len()
                        && computed.engine_version
                            == collagene_core::verification::ENGINE_VERSION
                        && computed.model
                            == collagene_core::review::reads::trace_validation::MODEL
                        && computed
                            .reads
                            .iter()
                            .zip(&bindings)
                            .all(|(read, binding)| read.source_revision == binding.revision) =>
                    {
                        status.set("Read evidence computed for these exact document events under the recorded threshold and alignment rules; this is not experimental sample confirmation.".into());
                        report.set(Some(Arc::new(BoundEvidence {
                            report: computed,
                            reference,
                            reads: bindings,
                            options: opts,
                            epoch: captured_epoch,
                            indices: input,
                            generation: token,
                        })));
                    }
                    Ok(Response {
                        result: AnalysisResult::Error(error),
                        ..
                    })
                    | Err(error) => status.set(error),
                    _ => status.set(
                        "Unexpected read-validation worker response or source/options binding."
                            .into(),
                    ),
                }
            });
        }
    };
    view! {
        <style>".read-validation-panel .read-source{display:flex;align-items:flex-start;gap:8px;font-size:12px;line-height:1.5}.read-validation-panel .read-source input{width:auto;flex-shrink:0;margin-top:4px}.read-validation-panel .read-grid{display:grid;grid-template-columns:1fr 1fr;gap:8px}.read-validation-panel .read-grid label{margin:0}.read-validation-panel pre{font:11px/1.6 monospace;white-space:pre-wrap;overflow-wrap:anywhere;max-height:220px;overflow:auto;background:#f2f6f2;padding:8px}.read-validation-panel .read-source-list{max-height:240px;overflow:auto}.read-validation-panel .read-observation{padding:9px 0;border-bottom:1px solid #e1e7e0;font-size:12px;line-height:1.5}.read-validation-panel .read-table{overflow:auto}.read-validation-panel table{font-size:11px;white-space:nowrap;border-collapse:collapse}.read-validation-panel td,.read-validation-panel th{padding:5px;border-bottom:1px solid #e1e7e0}.read-validation-panel ul{padding-left:18px}.read-validation-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="inspector-content read-validation-panel"><h2>"Read validation"</h2><p class="helper">"The active document is the independent DNA reference. Select immutable trace or explicitly encoded FASTQ reads. Low-quality, ambiguous and uncovered positions remain distinct."</p>
        <details open><summary>"Read sources (up to 32)"</summary><div class="read-source-list">{move||query_options.get().into_iter().take(1000).map(|(index,name)|view!{<label class="read-source"><input type="checkbox" aria-label=format!("Validation read {index}") prop:checked=move||selected.with(|s|s.contains(&index)) disabled=move||selected.with(|s|s.len()>=32&&!s.contains(&index)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();selected.update(|s|{if checked&&!s.contains(&index){s.push(index);s.sort_unstable();}else if !checked{s.retain(|i|*i!=index);}});}/>{name.chars().take(160).collect::<String>()}</label>}).collect_view()}</div></details>
        <h3>"Confidence and trimming"</h3><p class="helper">"SCF/ABIF thresholds address native raw caller scales. FASTQ uses its declared Phred encoding. Neither is independently recalibrated; missing scores cannot support consensus."</p>
        <label>"Trim policy"<select aria-label="Read trim policy" prop:value=move||mode.get() on:change=move|e|mode.set(event_target_value(&e))><option value="quality">"First/last passing quality windows"</option><option value="all">"Retain all calls; still filter confidence"</option></select></label>
        <div class="read-grid">{[("Source-score threshold","Read score threshold",threshold),("Window calls","Read quality window",window),("Passing calls/window","Read minimum passing",minimum_passing),("Minimum retained calls","Read minimum retained",minimum_retained),("Minimum source files","Read minimum support",support)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria type="number" prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</div>
        <label class="read-source"><input type="checkbox" aria-label="Read retain unscored" prop:checked=move||missing.get() on:change=move|e|missing.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Align unscored traces as evidence only"</label><label class="read-source"><input type="checkbox" aria-label="Read both orientations" prop:checked=move||both.get() on:change=move|e|both.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Search both read orientations"</label><label class="read-source"><input type="checkbox" aria-label="Read require both strands" prop:checked=move||both_support.get() on:change=move|e|both_support.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Require source support on both reference strands"</label>
        <details><summary>"Bounded alignment work"</summary><label>"Cells per alignment"<input aria-label="Read alignment cell budget" type="number" prop:value=move||matrix.get() on:input=move|e|matrix.set(event_target_value(&e))/></label><label>"Total alignment/proof cells"<input aria-label="Read total cell budget" type="number" prop:value=move||total.get() on:input=move|e|total.set(event_target_value(&e))/></label><p class="helper">"Exact matches use an upper-bound shortcut with an exhaustive placement scan. Imperfect reads use affine alignment plus an optimal-path uniqueness proof. Tied/unproven coordinates never become eligible consensus coverage."</p></details>
        <button class="button primary full" disabled=move||busy.get()||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna)||selected.with(|s|s.is_empty()) on:click=run>"Validate selected reads"</button><p class="read-validation-status" role="status">{move||status.get()}</p>
        <button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(worker)=active.get_value(){worker.close();}active.set_value(None);busy.set(false);report.set(None);status.set("Read validation cancelled.".into());}>"Cancel read validation"</button>
        {move||report.get().map(|bound|view!{<EvidenceView bound guard molecule on_consensus on_read on_record on_reference on_message=Callback::new(move|message|status.set(message))/>})}
        </div>
    }
}

#[component]
fn Pager(page: RwSignal<usize>, total: usize, label: String) -> impl IntoView {
    let pages = total.max(1).div_ceil(50);
    let label = StoredValue::new(label);
    view! {<div class="method-row-actions" role="group" aria-label=label.get_value()>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous page"</button>
        <span>{move||format!("{} · page {} of {}",total,page.get().min(pages-1)+1,pages)}</span>
        <button class="text-tool" disabled=move||{page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next page"</button>
    </div>}
}
#[component]
fn EvidenceView(
    bound: Arc<BoundEvidence>,
    guard: Callback<Arc<BoundEvidence>, Result<(), String>>,
    molecule: Memo<Molecule>,
    on_consensus: Callback<Molecule, Result<(), String>>,
    on_read: Callback<crate::reads::contigs::ReadNavigation>,
    on_record: Callback<RecordRequest, Result<(), String>>,
    on_reference: Callback<ReferenceNavigation, Result<(), String>>,
    on_message: Callback<String>,
) -> impl IntoView {
    let bound = StoredValue::new(bound);
    let uncovered_page = RwSignal::new(0usize);
    let insertion_page = RwSignal::new(0usize);
    let site_page = RwSignal::new(0usize);
    let all_sites = RwSignal::new(false);
    Effect::new(move |_| {
        let _ = all_sites.get();
        site_page.set(0);
    });
    let reference = Callback::new(move |(start, end): (usize, usize)| {
        let outcome = bound.with_value(|b| {
            guard.run(b.clone())?;
            if start > end || end > b.report.consensus.len() {
                return Err("Reference coordinate is unavailable.".into());
            }
            on_reference.run(ReferenceNavigation {
                reference: b.reference.clone(),
                selection: Some((start, end)),
                epoch: b.epoch,
            })?;
            Ok::<_, String>(())
        });
        if let Err(error) = outcome {
            on_message.run(error);
        }
    });
    let read_call = Callback::new(move |(index, call_index): (usize, usize)| {
        let outcome = bound.with_value(|b| {
            guard.run(b.clone())?;
            let binding = b
                .reads
                .get(index)
                .ok_or("The original read event is unavailable.")?;
            on_read.run(crate::reads::contigs::ReadNavigation {
                binding: crate::reads::contigs::ReadBinding {
                    document_id: binding.document_id.clone(),
                    event_id: binding.event_id.clone(),
                    revision: binding.revision.clone(),
                },
                call_index,
                epoch: b.epoch,
            });
            Ok::<_, String>(())
        });
        if let Err(error) = outcome {
            on_message.run(error);
        }
    });
    let consensus = Callback::new(move |()| {
        bound.with_value(|b| {
            guard.run(b.clone())?;
            let mut document = consensus_document(&molecule.get_untracked(), &b.report)?;
            document.metadata.insert(
                "read_validation:reference_event".into(),
                serde_json::to_string(&b.reference).map_err(|e| e.to_string())?,
            );
            document.metadata.insert(
                "read_validation:read_events".into(),
                serde_json::to_string(&b.reads).map_err(|e| e.to_string())?,
            );
            Ok::<_, String>(document)
        })
    });
    let save = move |_| {
        let outcome = bound.with_value(|b| {
            guard.run(b.clone())?;
            on_record.run(RecordRequest {
                reference: b.reference.clone(),
                reads: b.reads.clone(),
                options: b.options.clone(),
                epoch: b.epoch,
            })
        });
        match outcome {Ok(())=>on_message.run("Fresh source-bound check recording requested. Saved evidence does not restore current verification after reopening.".into()),Err(error)=>on_message.run(error)}
    };
    let download = move |_| {
        let result = bound.with_value(|b| {
            guard.run(b.clone())?;
            // Preserve the existing top-level scientific report fields while
            // adding exact document identities for this current display.
            let mut value = serde_json::to_value(&b.report).map_err(|e| e.to_string())?;
            if let Some(object) = value.as_object_mut() {
                object.insert(
                    "source_events".into(),
                    serde_json::json!({"reference":b.reference,"reads":b.reads}),
                );
            }
            let bytes = serde_json::to_vec_pretty(&value).map_err(|e| e.to_string())?;
            crate::io::download::ExportArtifact::prepare(
                "read-validation.json",
                "application/json",
                &bytes,
            )
        });
        crate::io::download::complete_download(result, move |result| {
            if let Err(error) = result {
                on_message.run(error);
            }
        });
    };
    view! {<section class="read-validation-result">
        <h3>{bound.with_value(|b|format!("{:?}",b.report.status))}</h3>
        <p>{bound.with_value(|b|format!("Eligible coverage {} / {} bases ({:.2}%)",b.report.eligible_coverage,b.report.consensus.len(),b.report.reference_coverage_fraction*100.))}</p>
        <p class="helper">"Evidence is bound to the exact reference and ordered read document events. Distinct files do not establish distinct biological samples. Saving recomputes the check from retained originals."</p>
        <button class="button primary full" on:click=save>"Save check on this construct"</button>
        <details><summary>"Bound source events"</summary>{bound.with_value(|b|view!{<p>{format!("Reference document {} · event {}",b.reference.document_id,b.reference.event_id)}</p><ol>{b.reads.iter().map(|r|view!{<li>{format!("Document {} · event {} · revision {}",r.document_id,r.event_id,r.revision)}</li>}).collect_view()}</ol>})}</details>
        <ul>{bound.with_value(|b|b.report.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view())}</ul>
        <h4>"Uncovered / insufficiently supported ranges"</h4>
        {move||bound.with_value(|b|view!{<ul>{b.report.uncovered_ranges.iter().skip(uncovered_page.get()*50).take(50).map(|&(start,end)|view!{<li><button class="text-tool" on:click=move|_|reference.run((start,end))>{format!("{}–{} (1-based inclusive)",start+1,end)}</button></li>}).collect_view()}</ul>})}
        <Pager page=uncovered_page total=bound.with_value(|b|b.report.uncovered_ranges.len()) label="Uncovered reference ranges".to_string()/>
        <h4>"Support-only consensus"</h4><pre class="read-consensus-preview">{bound.with_value(|b|format!("{}{}",b.report.consensus.chars().take(600).collect::<String>(),if b.report.consensus.len()>600{"… (preview; export retains all bases)"}else{""}))}</pre>
        <p class="helper">"Uncovered and unresolved positions remain explicit. This reference does not assert experimental sample identity or physical end/strand chemistry."</p>
        <button class="button secondary full" on:click=move|_|match consensus.run(()).and_then(|m|on_consensus.run(m)){Ok(())=>on_message.run("Support-only consensus reference opened; original read evidence remains separate.".into()),Err(error)=>on_message.run(error)}>"Open support-only consensus"</button>
        <button class="button secondary full" on:click=move|_|{let outcome=consensus.run(()).and_then(|m|crate::io::download::ExportArtifact::prepare("read-consensus.fasta","text/plain",collagene_core::formats::sequence::export_fasta(&m).as_bytes()));crate::io::download::complete_download(outcome,move|result|if let Err(error)=result{on_message.run(error);});}>"Download consensus FASTA"</button>
        <button class="button secondary full" on:click=download>"Download read validation JSON"</button>
        <h4>"Read mapping and native evidence"</h4>
        {bound.with_value(|b|(0..b.report.reads.len()).map(|index|view!{<ReadEvidence bound=b.clone() index reference read_call/>}).collect_view())}
        <h4>"Insertion observations"</h4>
        {move||bound.with_value(|b|b.report.insertions.iter().enumerate().skip(insertion_page.get()*50).take(50).map(|(index,_)|view!{<InsertionEvidence bound=b.clone() index reference read_call/>}).collect_view())}
        <Pager page=insertion_page total=bound.with_value(|b|b.report.insertions.len()) label="Insertion observations".to_string()/>
        <p class="helper">"Reference insertion boundaries are zero-based carets, not existing bases. Inserted calls are not added to the reference-coordinate consensus. Individual call eligibility is separate from the minimum-file and strand-support rules."</p>
        <h4>"Differences and unresolved sites"</h4>
        <label class="read-source"><input type="checkbox" aria-label="Show all reported reference sites" prop:checked=move||all_sites.get() on:change=move|e|all_sites.set(event_target_checked(&e))/>"Include supported matches"</label>
        {move||bound.with_value(|b|{let all=all_sites.get();let count=b.report.sites.iter().filter(|s|all||s.state!=SiteState::SupportedMatch).count();view!{<div>{b.report.sites.iter().filter(|s|all||s.state!=SiteState::SupportedMatch).skip(site_page.get()*50).take(50).map(|site|{let position=site.position;view!{<article class="read-observation"><button class="text-tool" on:click=move|_|reference.run((position,position+1))>{format!("Base {} · {:?} · {} → {}",position+1,site.state,site.reference_base,site.consensus_base)}</button><p>{format!("Eligible files {:?} · F/R {}/{} · low {} · ambiguous {} · deletion observations {} (eligible {})",site.eligible_read_indices,site.forward_reads,site.reverse_reads,site.low_quality_calls,site.ambiguous_calls,site.deletion_observations,site.eligible_deletions)}</p></article>}}).collect_view()}<Pager page=site_page total=count label="Reported reference sites".to_string()/></div>}})}
    </section>}
}

fn optional<T: std::fmt::Display>(value: Option<T>) -> String {
    value.map(|v| v.to_string()).unwrap_or_else(|| "—".into())
}
#[component]
fn ReadEvidence(
    bound: Arc<BoundEvidence>,
    index: usize,
    reference: Callback<(usize, usize)>,
    read_call: Callback<(usize, usize)>,
) -> impl IntoView {
    let bound = StoredValue::new(bound);
    let page = RwSignal::new(0usize);
    view! {<details class="read-observation"><summary>{bound.with_value(|b|{let read=&b.report.reads[index];format!("{}: {} · unique placement {}",index+1,read.name,read.placement_unique_proven)})}</summary>
        {bound.with_value(|b|{let read=&b.report.reads[index];view!{<div><p>{format!("Quality source {:?} · original source reverse to reference {} · format {:?}",read.trim.as_ref().map(|t|t.quality_source),read.original_read_reverse_to_reference,read.read_format)}</p>{read.error.clone().map(|error|view!{<p role="alert">{error}</p>})}{read.trim.as_ref().map(|trim|view!{<p>{format!("Retained current calls [{}, {}) of {}; passing native-score calls {}; low/unknown {}",trim.start,trim.end,trim.original_calls,trim.passing_calls,trim.low_or_unknown_calls.len())}</p>})}<p>{format!("Optimal score proven {} · search complete {} · optimal paths (saturated) {} · paired reference coverage {:.2}%",read.optimality_proven,read.search_complete,read.optimal_paths_saturated,read.paired_reference_fraction*100.)}</p>{read.alignment.as_ref().map(|a|view!{<pre>{format!("Alignment preview (first 300 characters per line); all columns paged below.\nScore {} · CIGAR {}\nR {}\nQ {}",a.score,a.cigar.chars().take(300).collect::<String>(),a.reference_aligned.chars().take(300).collect::<String>(),a.query_aligned.chars().take(300).collect::<String>())}</pre>})}<ul>{read.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul></div>}})}
        <div class="read-table"><table><thead><tr><th>"Reference / boundary"</th><th>"Call / original (1-based)"</th><th>"Peak / original (zero-based)"</th><th>"Native score"</th><th>"Column / quality"</th><th>"Eligible"</th></tr></thead><tbody>
            {move||bound.with_value(|b|b.report.reads[index].columns.iter().skip(page.get()*50).take(50).map(|column|{let region=column.reference_position.map(|p|(p,p+1)).unwrap_or((column.reference_boundary,column.reference_boundary));let call=column.call_index;view!{<tr><td><button class="text-tool" on:click=move|_|reference.run(region)>{column.reference_position.map(|p|(p+1).to_string()).unwrap_or_else(||format!("Caret {}",column.reference_boundary))}</button></td><td><button class="text-tool" disabled=call.is_none() on:click=move|_|if let Some(call)=call{read_call.run((index,call))}>{format!("{} / {} · {}",optional(column.call_index.map(|n|n+1)),optional(column.stored_call_index.map(|n|n+1)),column.observed_base.map(|c|c.to_string()).unwrap_or_else(||"gap".into()))}</button></td><td>{format!("{} / {}",optional(column.sample_index),optional(column.stored_sample_index))}</td><td>{optional(column.raw_confidence)}</td><td>{format!("{} · score passes {} · flanks pass {}",match column.kind{collagene_core::review::reads::trace_validation::ColumnKind::Pair=>"paired",collagene_core::review::reads::trace_validation::ColumnKind::Insertion=>"insertion",collagene_core::review::reads::trace_validation::ColumnKind::Deletion=>"deletion"},column.quality_passes,column.flank_quality_passes)}</td><td>{column.eligible.to_string()}</td></tr>}}).collect_view())}
        </tbody></table></div>
        <Pager page total=bound.with_value(|b|b.report.reads[index].columns.len()) label=format!("Read {} alignment columns",index+1)/>
        <p class="helper">"Call buttons navigate the original retained read document at its current molecule call coordinate. Original attachment call/peak indices are separate evidence; deletions have no invented read call."</p>
    </details>}
}
#[component]
fn InsertionEvidence(
    bound: Arc<BoundEvidence>,
    index: usize,
    reference: Callback<(usize, usize)>,
    read_call: Callback<(usize, usize)>,
) -> impl IntoView {
    let bound = StoredValue::new(bound);
    let page = RwSignal::new(0usize);
    view! {<details class="read-observation"><summary>{bound.with_value(|b|{let insertion=&b.report.insertions[index];format!("Boundary {} · {} inserted calls · read {}",insertion.reference_boundary,insertion.call_indices.len(),insertion.read_index+1)})}</summary>
        {bound.with_value(|b|{let insertion=&b.report.insertions[index];let boundary=insertion.reference_boundary;view!{<button class="text-tool" on:click=move|_|reference.run((boundary,boundary))>"Show reference insertion caret"</button><p>{format!("Eligible call evidence {} · original source reverse {}",insertion.eligible,insertion.original_read_reverse_to_reference)}</p>}})}
        <div class="read-table"><table><thead><tr><th>"Insertion offset (1-based)"</th><th>"Base"</th><th>"Current read call (1-based)"</th></tr></thead><tbody>{move||bound.with_value(|b|{let insertion=&b.report.insertions[index];insertion.call_indices.iter().enumerate().skip(page.get()*50).take(50).map(|(offset,&call)|{let read=insertion.read_index;view!{<tr><td>{offset+1}</td><td>{insertion.sequence.as_bytes().get(offset).map(|b|char::from(*b).to_string()).unwrap_or_else(||"unavailable".into())}</td><td><button class="text-tool" on:click=move|_|read_call.run((read,call))>{call+1}</button></td></tr>}}).collect_view()})}</tbody></table></div>
        <Pager page total=bound.with_value(|b|b.report.insertions[index].call_indices.len()) label=format!("Insertion {} calls",index+1)/>
    </details>}
}
