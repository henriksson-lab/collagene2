//! Selected, source-bound mismatch-PCR reference outputs. Scientific preparation
//! stays in the worker; naming and ordering never imply target verification.
use crate::cloning::methods::mutagenesis::OutputOpen;
use crate::platform::worker::CachedWorker;
use crate::primers::pcr::pcr_overlay;
use crate::primers::pcr::pcr_overlay::Overlay;
use crate::primers::pcr::pcr_overlay::Publication;
use crate::primers::pcr::pcr_overlay::Role;
use collagene_core::{
    cloning::pcr::{
        mismatch_outputs as mismatch_pcr_outputs, mutagenesis::PcrReport, outputs::Options,
    },
    formats::oligos::{self as oligo_formats, Delimiter, TableExportMode},
    material::lineage::EventBinding,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, AtomicU32, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

static NEXT_OVERLAY_OWNER: AtomicU32 = AtomicU32::new(1);

pub struct Review {
    pub source: EventBinding,
    pub target: EventBinding,
    pub epoch: u64,
    pub source_molecule: Arc<Molecule>,
    pub forward: Primer,
    pub reverse: Primer,
    pub report: PcrReport,
    pub documents: Vec<Molecule>,
    pub findings: Vec<String>,
}
fn valid_name(name: &str) -> bool {
    !name.is_empty()
        && name.trim() == name
        && name.len() <= 256
        && !name.chars().any(char::is_control)
}
#[component]
pub fn MismatchPcrOutputsPanel(
    review: Arc<Review>,
    guard: Callback<(), Result<(), String>>,
    on_outputs: Callback<OutputOpen, Result<(), String>>,
    on_overlay: Callback<Publication>,
    on_navigate: Callback<(EventBinding, u64, Vec<(usize, usize)>), Result<(), String>>,
) -> impl IntoView {
    let count = review.report.products.len();
    let initial_name = if count == 1 {
        review
            .documents
            .first()
            .map(|d| d.name.clone())
            .unwrap_or_default()
    } else {
        String::new()
    };
    let selected = RwSignal::new((count == 1).then_some(0usize));
    let name = RwSignal::new(initial_name);
    let save_product = RwSignal::new(true);
    let save_primers = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let sequence_page = RwSignal::new(0usize);
    let substitutions_page = RwSignal::new(0usize);
    let findings_page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let bound = StoredValue::new(review);
    let owner = StoredValue::new(format!(
        "mismatch-{}",
        NEXT_OVERLAY_OWNER.fetch_add(1, Ordering::Relaxed)
    ));
    let projected = RwSignal::new(None::<Arc<Overlay>>);
    let project = Callback::new(move |()| -> Result<Arc<Overlay>, String> {
        guard.run(())?;
        let index = selected
            .get_untracked()
            .ok_or("Choose a reported mismatch-PCR geometry first.")?;
        bound
            .with_value(|review| {
                pcr_overlay::mismatch(
                    &review.source_molecule,
                    &review.source,
                    review.epoch,
                    &review.forward,
                    &review.reverse,
                    &review.report,
                    index,
                )
            })
            .map(Arc::new)
    });
    Effect::new(move |_| {
        let overlay = if selected.get().is_some() {
            match project.run(()) {
                Ok(overlay) => Some(overlay),
                Err(error) => {
                    status.set(error);
                    None
                }
            }
        } else {
            None
        };
        projected.set(overlay.clone());
        on_overlay.run(Publication {
            owner: owner.get_value(),
            overlay,
        });
    });
    on_cleanup(move || {
        on_overlay.run(Publication {
            owner: owner.get_value(),
            overlay: None,
        })
    });
    let navigate = Callback::new(move |role: Role| {
        let result = project.run(()).and_then(|overlay| {
            let spans = overlay
                .spans
                .iter()
                .filter(|span| span.role == role)
                .map(|span| (span.start, span.end))
                .collect::<Vec<_>>();
            if spans.is_empty() {
                return Err("The selected geometry has no mapped bases for this role.".into());
            }
            on_overlay.run(Publication {
                owner: owner.get_value(),
                overlay: Some(overlay.clone()),
            });
            on_navigate.run((overlay.source.clone(), overlay.epoch, spans))
        });
        if let Err(error) = result {
            status.set(error);
        }
    });
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
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
        busy.set(false);
        status.set("Output preparation cancelled. No outputs were created.".into());
    });
    let export = Callback::new(move |delimiter: Delimiter| {
        if busy.get_untracked() {
            return;
        }
        if let Err(error) = guard.run(()) {
            status.set(error);
            return;
        }
        let result = bound
            .with_value(|r| {
                collagene_core::planning::orders::from_primers(&[
                    r.forward.clone(),
                    r.reverse.clone(),
                ])
            })
            .and_then(|table| {
                oligo_formats::export_oligos(&table, delimiter, TableExportMode::SpreadsheetSafe)
            });
        if let Err(error) = guard.run(()) {
            status.set(error);
            return;
        }
        match result {
            Ok(report) => {
                let (name, mime) = match delimiter {
                    Delimiter::Comma => ("mismatch-pcr-primers.csv", "text/csv;charset=utf-8"),
                    Delimiter::Tab => (
                        "mismatch-pcr-primers.tsv",
                        "text/tab-separated-values;charset=utf-8",
                    ),
                };
                crate::io::download::download_with_completion(
                    name,
                    mime,
                    report.text.as_bytes(),
                    move |result| {
                        if status.try_get_untracked().is_none() {
                            return;
                        }
                        match result {
                            Ok(()) => status.set(report.warnings.join(" ")),
                            Err(error) => status.set(error),
                        }
                    },
                )
            }
            Err(error) => status.set(error.to_string()),
        }
    });
    let create = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                guard.run(())?;
                let index = selected
                    .get_untracked()
                    .ok_or("Choose a reported mismatch-PCR geometry.")?;
                let options = Options {
                    name: name.get_untracked(),
                    save_product: save_product.get_untracked(),
                    save_primers: save_primers.get_untracked(),
                };
                options.validate().map_err(|e| e.to_string())?;
                let request = bound.with_value(|r| -> Result<_, String> {
                    let product = r
                        .documents
                        .get(index)
                        .ok_or("The selected preview is unavailable; simulate again.")?;
                    Ok(mismatch_pcr_outputs::Request {
                        source: r.source_molecule.as_ref().clone(),
                        forward: r.forward.clone(),
                        reverse: r.reverse.clone(),
                        options: r.report.options.clone(),
                        selected_index: index,
                        expected_product_revision: product.revision(),
                        outputs: options.clone(),
                    })
                })?;
                let revision =
                    mismatch_pcr_outputs::revision(&request).map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((index, options, request, revision, w))
            })();
            let (index, options, request, revision, w) = match prepared {
                Ok(v) => v,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            worker.set_value(Some(w.clone()));
            busy.set(true);
            status.set(
                "Reconstructing the selected mismatch-PCR preview before creating outputs…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::PrepareMismatchPcrOutputs { request },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let outcome = (|| -> Result<(), String> {
                    guard.run(())?;
                    let response = result?;
                    if response.id != token || response.revision != revision {
                        return Err(
                            "The output response differs from the selected PCR request.".into()
                        );
                    }
                    let prepared = match response.result {
                        AnalysisResult::MismatchPcrOutputs(result) => result,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected mismatch-PCR output response.".into()),
                    };
                    prepared
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if prepared.selected_index != index
                        || prepared.product.is_some() != options.save_product
                        || prepared.primers.is_some() != options.save_primers
                    {
                        return Err("The prepared outputs differ from your selection.".into());
                    }
                    let r = bound.get_value();
                    guard.run(())?;
                    on_outputs.run(OutputOpen {
                        source: r.source.clone(),
                        target: r.target.clone(),
                        epoch: r.epoch,
                        revision,
                        options,
                        prepared,
                    })
                })();
                // A successful host commit may change the active document and
                // dispose this review before its callback returns.
                if !alive.load(Ordering::Relaxed) {
                    return;
                }
                match outcome {Ok(())=>status.set("Selected mismatch-PCR reference outputs created. No target verification or laboratory success is asserted.".into()),Err(error)=>status.set(error)}
            });
        }
    };
    let download_report = move |_| {
        if let Err(error) = guard.run(()) {
            status.set(error);
            return;
        }
        let result = bound.with_value(|r| serde_json::to_vec_pretty(&r.report));
        match result {
            Ok(bytes) => crate::io::download::download_with_completion(
                "mismatch-pcr.json",
                "application/json",
                &bytes,
                move |result| {
                    if status.try_get_untracked().is_none() {
                        return;
                    }
                    match result {
                        Ok(()) => {
                            status.set("Complete returned mismatch-PCR evidence downloaded.".into())
                        }
                        Err(error) => status.set(error),
                    }
                },
            ),
            Err(error) => status.set(error.to_string()),
        }
    };
    view! {
        <section class="pcr-output-review mutation-pcr-result">
            <h3>{format!("{count} conditional PCR geometries")}</h3>
            <p class="helper">{bound.with_value(|r|format!("Search complete: {} · complete under declared model: {} · pair checks {} · unsupported geometries {} · conflicting overlaps {} · ambiguous bindings {} · ambiguous products {} · rejected geometries {}",r.report.search_complete,r.report.enumeration_complete,r.report.pair_checks,r.report.unsupported_geometries,r.report.conflicting_overlaps,r.report.ambiguous_bindings,r.report.ambiguous_products,r.report.rejected_geometries))}</p>
            <p class="helper">"Choose the reference preview to retain. Binding rules, primer-derived substitutions and reference letters do not establish amplification, physical end chemistry or a match to the independent target. Incomplete searches remain incomplete."</p>
            <details><summary>"Original template and declared search"</summary><p>{bound.with_value(|r|format!("Document {} · event {} · revision {}",r.source.document_id,r.source.event_id,r.source.revision))}</p><pre>{bound.with_value(|r|serde_json::to_string_pretty(&r.report.options).unwrap_or_default())}</pre></details>
            <details><summary>"Model warnings"</summary>{bound.with_value(|r|r.report.warnings.clone()).into_iter().map(|warning|view!{<p class="helper">{warning}</p>}).collect_view()}</details>
            <details open><summary>"Named oligos · ordered 5′→3′"</summary>
                {bound.with_value(|r|[r.forward.clone(),r.reverse.clone()]).into_iter().enumerate().map(|(index,primer)|view!{
                    <article><strong>{format!("Oligo {}: {}",index+1,primer.name)}</strong><p class="sequence-preview">{primer.ordered_sequence()}</p><p class="helper">{format!("5′ tail: {} · annealing domain: {}",if primer.tail.is_empty(){"(none)"}else{&primer.tail},primer.annealing.as_str())}</p></article>
                }).collect_view()}
                <p class="helper">"The search uses both oligos in both roles, including same-oligo amplification. Order sheets retain both declared oligos, without inferred quantities or modifications."</p>
                <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get() on:click=move|_|export.run(Delimiter::Comma)>"Primer CSV"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|export.run(Delimiter::Tab)>"Primer TSV"</button></div>
            </details>
            {move||{
                let current=page.get().min(count.saturating_sub(1)/10);
                let pages=count.div_ceil(10).max(1);
                let rows=bound.with_value(|r|r.report.products.iter().enumerate().skip(current*10).take(10).map(|(index,geometry)|{
                    let default_name=r.documents.get(index).map(|d|d.name.clone()).unwrap_or_default();
                    let names=[r.forward.name.as_str(),r.reverse.name.as_str()];
                    let forward=names.get(geometry.forward_primer_index).copied().unwrap_or("unresolved oligo");
                    let reverse=names.get(geometry.reverse_primer_index).copied().unwrap_or("unresolved oligo");
                    view!{<article><label><input type="radio" name="mismatch-pcr-output-geometry" prop:checked=move||selected.get()==Some(index) disabled=move||busy.get() on:change=move|_|{selected.set(Some(index));name.set(default_name.clone());sequence_page.set(0);substitutions_page.set(0);status.set(String::new());}/>{format!("Product {} · {} bp · {} introduced substitutions",index+1,geometry.sequence.len(),geometry.substitutions.len())}</label>
                        <p class="helper">{format!("Forward-role {forward} at base {} · reverse-role {reverse} at base {} · {} template bases{} · overlapping domains: {} · concrete: {}",geometry.forward_start+1,geometry.reverse_start+1,geometry.template_length,if geometry.crosses_origin{" · crosses template origin"}else{""},geometry.overlapping_domains,geometry.certain)}</p>
                    </article>}
                }).collect_view());
                view!{<div>{rows}</div><div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous products"</button><span>{format!("Page {} of {pages}",current+1)}</span><button class="text-tool" disabled=move||{busy.get()||current+1>=pages} on:click=move|_|page.set((current+1).min(pages-1))>"Next products"</button></div>}
            }}
            {move||selected.get().and_then(|index|bound.with_value(|r|r.report.products.get(index).map(|geometry|{
                let length=geometry.sequence.len();
                let pages=length.div_ceil(500).max(1);let current=sequence_page.get().min(pages-1);
                let start=(current*500).min(length);let end=(start+500).min(length);
                let sequence=geometry.sequence.as_str()[start..end].to_string();
                let substitutions=geometry.substitutions.len();let change_pages=substitutions.div_ceil(30).max(1);let changes_page=substitutions_page.get().min(change_pages-1);
                let changed=geometry.substitutions.iter().skip(changes_page*30).take(30).map(|change|view!{<tr><td>{change.source_position+1}</td><td>{change.product_position+1}</td><td>{format!("{} → {}",change.template_base,change.product_base)}</td><td>{if change.template_ambiguous{"ambiguous source"}else{"concrete source"}}</td></tr>}).collect_view();
                view!{<details open><summary>"Selected reference and primer-introduced changes"</summary>
                    <p class="helper">{format!("Product bases {}–{} of {length}; page {} of {pages}. {} projected annotations.",start+1,end,current+1,r.documents.get(index).map(|d|d.features.len()).unwrap_or_default())}</p><pre class="sequence-preview">{sequence}</pre>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||current==0 on:click=move|_|sequence_page.set(current.saturating_sub(1))>"Previous bases"</button><button class="text-tool" disabled=move||{busy.get()||current+1>=pages} on:click=move|_|sequence_page.set((current+1).min(pages-1))>"Next bases"</button></div>
                    <p class="helper">{format!("Domain-relative mismatch offsets (one-based): forward {:?}; reverse {:?}. Ambiguous offsets: forward {:?}; reverse {:?}. Each domain is shown in its oligo's own 5′→3′ orientation.",geometry.forward_binding.mismatches.iter().map(|n|n+1).collect::<Vec<_>>(),geometry.reverse_binding.mismatches.iter().map(|n|n+1).collect::<Vec<_>>(),geometry.forward_binding.ambiguous_positions.iter().map(|n|n+1).collect::<Vec<_>>(),geometry.reverse_binding.ambiguous_positions.iter().map(|n|n+1).collect::<Vec<_>>())}</p>
                    <p>{format!("{substitutions} primer-introduced substitutions · page {} of {change_pages}",changes_page+1)}</p><table><thead><tr><th>"Template base"</th><th>"Product base"</th><th>"Change"</th><th>"Evidence"</th></tr></thead><tbody>{changed}</tbody></table>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||changes_page==0 on:click=move|_|substitutions_page.set(changes_page.saturating_sub(1))>"Previous changes"</button><button class="text-tool" disabled=move||{busy.get()||changes_page+1>=change_pages} on:click=move|_|substitutions_page.set((changes_page+1).min(change_pages-1))>"Next changes"</button></div>
                </details>}
            })))}
            <details><summary>"Annotation projection findings"</summary>{move||{
                let count=bound.with_value(|r|r.findings.len());let pages=count.div_ceil(20).max(1);let current=findings_page.get().min(pages-1);
                let findings=bound.with_value(|r|r.findings.iter().skip(current*20).take(20).cloned().map(|finding|view!{<p class="helper">{finding}</p>}).collect_view());
                view!{<div>{findings}</div><p>{format!("Findings across all returned products · page {} of {pages}",current+1)}</p><div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||current==0 on:click=move|_|findings_page.set(current.saturating_sub(1))>"Previous findings"</button><button class="text-tool" disabled=move||{busy.get()||current+1>=pages} on:click=move|_|findings_page.set((current+1).min(pages-1))>"Next findings"</button></div>}
            }}</details>
            {move||projected.get().map(|overlay|view!{
                <section class="pcr-placement-review"><h4>"Selected template placement"</h4>
                    <p class="helper">{overlay.summary.clone()}</p>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||{busy.get()} on:click=move|_|{navigate.run(Role::Amplified);}>"View template interval"</button><button class="text-tool" disabled=move||{busy.get()} on:click=move|_|{navigate.run(Role::Forward);}>"View forward-role primer"</button><button class="text-tool" disabled=move||{busy.get()} on:click=move|_|{navigate.run(Role::Reverse);}>"View reverse-role primer"</button></div>
                    <p class="helper">"Only this selected conditional geometry is overlaid. Mismatch and ambiguous bases retain their reported evidence; synthetic tails have no template coordinates. Other returned geometries remain separate."</p>
                    {overlay.tails.iter().cloned().map(|tail|view!{<p class="sequence-preview">{tail}</p>}).collect_view()}
                    {overlay.warnings.iter().cloned().map(|warning|view!{<p class="helper">{warning}</p>}).collect_view()}
                </section>
            })}
            <label>"Output name"<input maxlength="256" disabled=move||busy.get() prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
            <label><input type="checkbox" disabled=move||busy.get() prop:checked=move||save_product.get() on:change=move|e|save_product.set(event_target_checked(&e))/>"Create selected mismatch-PCR reference preview"</label>
            <label><input type="checkbox" disabled=move||busy.get() prop:checked=move||save_primers.get() on:change=move|e|save_primers.set(event_target_checked(&e))/>"Save declared pair as an oligo collection"</label>
            <button class="button primary full" disabled=move||busy.get()||selected.get().is_none()||!valid_name(&name.get())||(!save_product.get()&&!save_primers.get()) on:click=create>"Create selected outputs"</button>
            <button class="text-tool" disabled=move||!busy.get() on:click=move|_|cancel.run(())>"Cancel output preparation"</button>
            <button class="button secondary full" disabled=move||busy.get() on:click=download_report>"Download complete mismatch-PCR evidence"</button>
            <p class="helper" role="status">{move||status.get()}</p>
        </section>
    }.into_any()
}
