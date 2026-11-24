//! Explicit revision-bound read/reference evidence; no imported or AI-asserted verdict.
use crate::worker;
use collagene_core::{
    Alphabet, Molecule,
    protocol::{AnalysisJob, AnalysisResult, Request, Response},
    trace_validation::{
        MissingQualityPolicy, SiteState, TraceValidationOptions, TraceValidationReport, TrimMode,
        TrimOptions,
    },
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

fn number<T: std::str::FromStr>(text: String, label: &str) -> Result<T, String> {
    text.parse()
        .map_err(|_| format!("{label} is not a valid number"))
}
fn bound_reads(
    get_query: Callback<usize, Option<Molecule>>,
    indices: &[usize],
) -> Option<Vec<String>> {
    indices
        .iter()
        .map(|&i| get_query.run(i).map(|m| m.revision()))
        .collect()
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
) -> impl IntoView {
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
    let report = RwSignal::new(None::<Arc<TraceValidationReport>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    Effect::new(move |_| {
        let _ = (
            molecule.with(Molecule::revision),
            bound_reads(get_query, &selected.get()),
            options.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        busy.set(false);
        status.set(String::new());
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            let input = selected.get_untracked();
            let source = molecule.get_untracked();
            let reads = match input
                .iter()
                .map(|&i| {
                    get_query
                        .run(i)
                        .ok_or_else(|| "A read document is no longer open.".to_string())
                })
                .collect::<Result<Vec<_>, _>>()
            {
                Ok(reads) => reads,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            if reads.is_empty() {
                status.set("Choose at least one trace or quality-bearing FASTQ read.".into());
                return;
            }
            let opts = match options.get_untracked() {
                Ok(options) => options,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let revision = source.revision();
            let revisions = reads.iter().map(Molecule::revision).collect::<Vec<_>>();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Aligning reads and reconstructing quality/coverage evidence…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::ValidateTraces {
                        molecule: source,
                        reads,
                        options: opts,
                    },
                })
                .await;
                if !alive.load(Ordering::Relaxed) {
                    return;
                }
                if generation.get_untracked() != token
                    || molecule.with_untracked(Molecule::revision) != revision
                    || untrack(|| bound_reads(get_query, &input)) != Some(revisions.clone())
                {
                    return;
                }
                busy.set(false);
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result: AnalysisResult::TraceValidation(computed),
                    }) if id == token
                        && returned == revision
                        && computed.reference_revision == revision
                        && computed.read_revisions == revisions =>
                    {
                        status.set("Read evidence computed under the recorded threshold and alignment rules; this is not experimental sample confirmation.".into());
                        report.set(Some(Arc::new(computed)));
                    }
                    Ok(Response {
                        result: AnalysisResult::Error(error),
                        ..
                    })
                    | Err(error) => status.set(error),
                    _ => status.set("Unexpected read-validation worker response.".into()),
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
        {move||report.get().map(|computed|{
            let json=computed.clone();let consensus=computed.clone();let fasta=computed.clone();let status_text=format!("{:?}",computed.status);let uncovered=computed.uncovered_ranges.iter().take(100).map(|(a,b)|format!("{}–{}",a+1,b)).collect::<Vec<_>>().join(", ");
            view!{<section class="read-validation-result"><h3>{status_text}</h3><p>{format!("Eligible coverage {} / {} bases ({:.2}%)",computed.eligible_coverage,computed.consensus.len(),computed.reference_coverage_fraction*100.)}</p><p class="helper">{format!("Uncovered/insufficiently supported intervals (1-based): {}{}",if uncovered.is_empty(){"none"}else{&uncovered},if computed.uncovered_ranges.len()>100{" · first 100 shown"}else{""})}</p><ul>{computed.warnings.iter().cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>
            <h4>"Support-only consensus"</h4><pre class="read-consensus-preview">{format!("{}{}",computed.consensus.chars().take(600).collect::<String>(),if computed.consensus.len()>600{"… (preview; export retains all bases)"}else{""})}</pre>
            <button class="button secondary full" on:click=move|_|match consensus_document(&molecule.get_untracked(),&consensus).and_then(|m|on_consensus.run(m)){Ok(())=>status.set("Consensus reference opened. Source documents and original evidence remain separate.".into()),Err(error)=>status.set(error)}>"Open support-only consensus"</button>
            <button class="button secondary full" on:click=move|_|match consensus_document(&molecule.get_untracked(),&fasta){Ok(m)=>{let text=collagene_core::formats::export_fasta(&m);if let Err(error)=crate::app::download("read-consensus.fasta","text/plain",text.as_bytes()){status.set(error)}},Err(error)=>status.set(error)}>"Download consensus FASTA"</button>
            <button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(json.as_ref()){Ok(bytes)=>if let Err(error)=crate::app::download("read-validation.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}>"Download read validation JSON"</button>
            <h4>"Read mapping and native evidence"</h4>{computed.reads.iter().enumerate().map(|(i,read)|{let display=read;let alignment=read.alignment.as_ref().map(|a|format!("Score {} · CIGAR {}\nR {}\nQ {}",a.score,a.cigar.chars().take(300).collect::<String>(),a.reference_aligned.chars().take(300).collect::<String>(),a.query_aligned.chars().take(300).collect::<String>()));view!{<details class="read-observation"><summary>{format!("{}: {} · unique placement {}",i+1,display.name,display.placement_unique_proven)}</summary><p>{format!("{:?} · source {}… · orientation reverse {} · native format {:?}",display.trim.as_ref().map(|t|t.quality_source),display.source_revision.chars().take(12).collect::<String>(),display.original_read_reverse_to_reference,display.read_format)}</p>{display.error.clone().map(|e|view!{<p role="alert">{e}</p>})}{display.trim.as_ref().map(|t|view!{<p>{format!("Retained calls [{}, {}) of {}; passing native-score calls {}; low/unknown {}",t.start,t.end,t.original_calls,t.passing_calls,t.low_or_unknown_calls.len())}</p>})}<p>{format!("Optimal score proven {} · search complete {} · optimal paths (saturated) {} · paired reference coverage {:.2}%",display.optimality_proven,display.search_complete,display.optimal_paths_saturated,display.paired_reference_fraction*100.)}</p>{alignment.map(|a|view!{<pre>{a}</pre>})}<ul>{display.warnings.iter().cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><div class="read-table"><table><thead><tr><th>"Reference"</th><th>"Call / original"</th><th>"Peak / original"</th><th>"Score"</th><th>"Eligible"</th></tr></thead><tbody>{display.columns.iter().take(100).cloned().map(|c|view!{<tr><td>{c.reference_position.map(|p|(p+1).to_string()).unwrap_or_else(||format!("between at {}",c.reference_boundary))}</td><td>{format!("{:?} / {:?}",c.call_index,c.stored_call_index)}</td><td>{format!("{:?} / {:?}",c.sample_index,c.stored_sample_index)}</td><td>{format!("{:?}",c.raw_confidence)}</td><td>{c.eligible.to_string()}</td></tr>}).collect_view()}</tbody></table></div><p class="helper">"First 100 alignment columns; full calls, stored-call/sample coordinates and confidence evidence are in JSON."</p></details>}}).collect_view()}
            <h4>"Insertion observations"</h4>{computed.insertions.iter().take(100).map(|insertion|{let boundary=insertion.reference_boundary;let n=computed.consensus.len();view!{<article class="read-observation"><button class="text-tool" on:click=move|_|{let p=boundary.min(n.saturating_sub(1));selected_feature.set(None);selection.set(Some((p,p+1)));}>{format!("Boundary {} · insert {} · read {}",boundary,insertion.sequence.chars().take(100).collect::<String>(),insertion.read_index+1)}</button><p>{format!("Eligible call evidence {} · original source reverse {} · call indices {:?}",insertion.eligible,insertion.original_read_reverse_to_reference,insertion.call_indices.iter().take(100).collect::<Vec<_>>())}</p></article>}}).collect_view()}<p class="helper">"First 100 insertion observations. Eligibility applies to individual calls; conclusive support also requires the recorded minimum file and strand rules. Boundaries are zero-based; inserted bases are not added to the reference-coordinate consensus."</p>
            <h4>"Differences and unresolved sites"</h4>{computed.sites.iter().filter(|s|s.state!=SiteState::SupportedMatch).take(200).map(|site|{let position=site.position;let label=format!("Base {} · {:?} · {} → {}",position+1,site.state,site.reference_base,site.consensus_base);view!{<article class="read-observation"><button class="text-tool" on:click=move|_|{selected_feature.set(None);selection.set(Some((position,position+1)));}>{label}</button><p>{format!("Eligible files {:?} · F/R {}/{} · low {} · ambiguous {} · deletion observations {} (eligible {})",site.eligible_read_indices,site.forward_reads,site.reverse_reads,site.low_quality_calls,site.ambiguous_calls,site.deletion_observations,site.eligible_deletions)}</p></article>}}).collect_view()}<p class="helper">"First 200 nonmatching/unresolved sites. Uncovered intervals are listed above. Insertion/deletion observations and all per-call evidence remain in the full JSON report."</p>
            </section>}
        })}
        </div>
    }
}
