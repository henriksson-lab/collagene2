//! Worker-backed revision-bound pairwise comparison; no externally asserted verification state.
use collagene_core::{
    Molecule,
    alignment::{Alignment, AlignmentMode, AlignmentOptions, Scoring},
    analysis, protocol,
};
use leptos::prelude::*;
use serde::Serialize;
use std::{cell::Cell, rc::Rc};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::spawn_local;

const STYLE: &str = ".compare-panel .compare-scores{display:grid;grid-template-columns:1fr 1fr;gap:8px}.compare-panel .compare-scores label{margin:0}.compare-panel .compare-result{font-size:13px;line-height:1.6;overflow-wrap:anywhere}.compare-panel pre{font:11px/1.6 monospace;overflow:auto;max-height:380px;white-space:pre}.compare-panel .compare-warning{font-size:12px;color:#796a3e}.compare-panel .compare-difference{padding:8px 0;border-bottom:1px solid #edf1ed;font-size:12px}.compare-panel .compare-difference button{border:0;border-radius:4px;background:#edf4ef;color:#365147;padding:4px 7px}.compare-panel .compare-difference code{display:block;overflow-wrap:anywhere;white-space:pre-wrap}.compare-panel .compare-status{font-size:13px;line-height:1.5;overflow-wrap:anywhere}.compare-panel .compare-results-meta{font-size:12px;color:#718274}.compare-panel details{margin:12px 0}.compare-panel summary{cursor:pointer}.compare-panel .compare-stat-list{padding-left:18px}";
#[derive(Clone, Serialize)]
struct BoundAlignment {
    reference_revision: String,
    query_revision: String,
    reference_name: String,
    query_name: String,
    query_index: usize,
    alignment: Alignment,
}
fn preview(text: &str, maximum: usize) -> String {
    if text.len() <= maximum {
        text.into()
    } else {
        format!("{}…", &text[..maximum])
    }
}
fn aligned_preview(alignment: &Alignment) -> String {
    use std::fmt::Write;
    let mut output = String::new();
    let limit = alignment.reference_aligned.len().min(4000);
    let mut r = alignment.reference_start;
    let mut q = alignment.query_start;
    for offset in (0..limit).step_by(60) {
        let end = (offset + 60).min(limit);
        let reference = &alignment.reference_aligned[offset..end];
        let query = &alignment.query_aligned[offset..end];
        let marks: String = reference
            .bytes()
            .zip(query.bytes())
            .map(|(a, b)| {
                if a == b'-' || b == b'-' {
                    ' '
                } else {
                    let am = analysis::base_mask(a);
                    let bm = analysis::base_mask(b);
                    if am & bm == 0 {
                        '×'
                    } else if am.count_ones() == 1 && bm.count_ones() == 1 {
                        '|'
                    } else {
                        '?'
                    }
                }
            })
            .collect();
        let _ = writeln!(
            output,
            "R {r:>7} {reference}\n          {marks}\nQ {q:>7} {query}\n"
        );
        r += reference.bytes().filter(|b| *b != b'-').count();
        q += query.bytes().filter(|b| *b != b'-').count();
    }
    output
}
fn download(report: &BoundAlignment) -> Result<(), String> {
    let bytes = serde_json::to_vec_pretty(report).map_err(|e| e.to_string())?;
    let parts = js_sys::Array::new();
    parts.push(&js_sys::Uint8Array::from(bytes.as_slice()));
    let opts = web_sys::BlobPropertyBag::new();
    opts.set_type("application/json");
    let blob = web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &opts)
        .map_err(|_| "Cannot prepare comparison download")?;
    let url = web_sys::Url::create_object_url_with_blob(&blob)
        .map_err(|_| "Cannot create comparison download")?;
    let window = web_sys::window().ok_or("No browser window")?;
    let doc = window.document().ok_or("No browser document")?;
    let anchor = doc
        .create_element("a")
        .map_err(|_| "Cannot create comparison download link")?;
    anchor
        .set_attribute("href", &url)
        .map_err(|_| "Cannot set comparison download URL")?;
    anchor
        .set_attribute("download", "pairwise-comparison.json")
        .map_err(|_| "Cannot set comparison filename")?;
    let anchor = anchor.unchecked_into::<web_sys::HtmlElement>();
    if let Some(body) = doc.body() {
        let _ = body.append_child(&anchor);
    }
    anchor.click();
    anchor.remove();
    let cleanup = Closure::once_into_js(move || {
        let _ = web_sys::Url::revoke_object_url(&url);
    });
    let _ =
        window.set_timeout_with_callback_and_timeout_and_arguments_0(cleanup.unchecked_ref(), 1000);
    Ok(())
}

#[component]
pub fn ComparePanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    selection: RwSignal<Option<(usize, usize)>>,
) -> impl IntoView {
    let selected = RwSignal::new(0usize);
    let mode = RwSignal::new("global".to_owned());
    let match_score = RwSignal::new("2".to_owned());
    let ambiguity = RwSignal::new("0".to_owned());
    let mismatch = RwSignal::new("-3".to_owned());
    let gap_open = RwSignal::new("-4".to_owned());
    let gap_extend = RwSignal::new("-1".to_owned());
    let max_cells = RwSignal::new("16000000".to_owned());
    let options = Memo::new(move |_| {
        let number = |value: String, label: &str| {
            value
                .parse::<i32>()
                .map_err(|_| format!("{label} must be an integer"))
        };
        Ok::<_, String>(AlignmentOptions {
            mode: match mode.get().as_str() {
                "local" => AlignmentMode::Local,
                "query" => AlignmentMode::QueryInReference,
                _ => AlignmentMode::Global,
            },
            scoring: Scoring {
                match_score: number(match_score.get(), "Match score")?,
                ambiguous_score: number(ambiguity.get(), "Ambiguity score")?,
                mismatch_score: number(mismatch.get(), "Mismatch score")?,
                gap_open: number(gap_open.get(), "Gap opening")?,
                gap_extend: number(gap_extend.get(), "Gap extension")?,
            },
            max_cells: max_cells
                .get()
                .parse()
                .map_err(|_| "Matrix budget must be a positive integer".to_owned())?,
        })
    });
    let report = RwSignal::new(None::<BoundAlignment>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Rc::new(Cell::new(true));
    let lifetime = StoredValue::new_local(alive.clone());
    on_cleanup(move || lifetime.with_value(|alive| alive.set(false)));
    Effect::new(move |_| {
        let choices = query_options.get();
        let index = selected.get();
        if !choices.iter().any(|(i, _)| *i == index) {
            if let Some((first, _)) = choices.first() {
                selected.set(*first);
            }
        }
        let _ = molecule.get().revision();
        let _ = get_query.run(selected.get()).map(|m| m.revision());
        let _ = options.get();
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        status.set(String::new());
        busy.set(false);
    });
    let compute = move |_| {
        let reference = molecule.get_untracked();
        let query_index = selected.get_untracked();
        let Some(query) = untrack(|| get_query.run(query_index)) else {
            status.set("Choose an open query document.".into());
            return;
        };
        let parameters = match options.get_untracked() {
            Ok(o) => o,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let reference_revision = reference.revision();
        let query_revision = query.revision();
        let reference_name = reference.name.clone();
        let query_name = query.name.clone();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        report.set(None);
        busy.set(true);
        status.set("Aligning in the Rust worker…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let response = crate::worker::analyze(protocol::Request {
                id,
                revision: reference_revision.clone(),
                job: protocol::AnalysisJob::Align {
                    molecule: reference,
                    query,
                    options: parameters,
                },
            })
            .await;
            if !alive.get() {
                return;
            }
            if generation.get_untracked() != id {
                return;
            }
            if molecule.get_untracked().revision() != reference_revision
                || untrack(|| get_query.run(query_index))
                    .is_none_or(|q| q.revision() != query_revision)
                || options.get_untracked().ok() != Some(parameters)
            {
                busy.set(false);
                report.set(None);
                status.set("Sources/options changed; compare again.".into());
                return;
            }
            busy.set(false);
            match response {
                Ok(protocol::Response {
                    id: response_id,
                    revision,
                    result: protocol::AnalysisResult::Alignment(alignment),
                }) if response_id == id && revision == reference_revision => {
                    report.set(Some(BoundAlignment {
                        reference_revision,
                        query_revision,
                        reference_name,
                        query_name,
                        query_index,
                        alignment,
                    }));
                    status.set("Pairwise comparison completed. Review coverage, uncertainty and unaligned flanks.".into());
                }
                Ok(protocol::Response {
                    result: protocol::AnalysisResult::Error(error),
                    ..
                })
                | Err(error) => status.set(error),
                _ => status.set("Unexpected or unbound alignment response".into()),
            }
        });
    };
    let export = move |_| {
        let Some(result) = report.get_untracked() else {
            return;
        };
        if molecule.get_untracked().revision() != result.reference_revision
            || untrack(|| get_query.run(result.query_index))
                .is_none_or(|q| q.revision() != result.query_revision)
            || options.get_untracked().ok() != Some(result.alignment.options)
        {
            report.set(None);
            status.set("Sources/options changed; compare again before exporting.".into());
            return;
        }
        if let Err(error) = download(&result) {
            status.set(error);
        }
    };
    view! {
        <style>{STYLE}</style>
        <div class="inspector-content compare-panel">
            <h2>"Pairwise comparison"</h2>
            <p class="helper">"The active document is the reference. Align another open DNA/RNA document with explicit scoring; uncertain bases and uncovered flanks stay distinct."</p>
            <label>"Query document"<select aria-label="Query document" prop:value=move||selected.get().to_string() on:change=move|e|selected.set(event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <label>"Alignment mode"<select aria-label="Alignment mode" prop:value=move||mode.get() on:change=move|e|mode.set(event_target_value(&e))><option value="global">"Global · both end-to-end"</option><option value="local">"Local · matching regions"</option><option value="query">"Query end-to-end in reference"</option></select></label>
            <details><summary>"Scoring and compute budget"</summary><p class="helper">"Gap length k costs opening + (k−1) × extension. Compatible ambiguity is never counted as a certain match."</p><div class="compare-scores">{[("Certain match",match_score),("Compatible ambiguity",ambiguity),("Mismatch",mismatch),("Gap opening",gap_open),("Gap extension",gap_extend)].into_iter().map(|(label,field)|view!{<label>{label}<input type="number" aria-label=label prop:value=move||field.get() on:input=move|e|field.set(event_target_value(&e))/></label>}).collect_view()}</div><label>"Matrix cell budget"<select aria-label="Matrix cell budget" prop:value=move||max_cells.get() on:change=move|e|max_cells.set(event_target_value(&e))><option value="4000000">"4 million"</option><option value="16000000">"16 million"</option><option value="64000000">"64 million"</option></select></label></details>
            <button class="button primary full" disabled=move||busy.get() on:click=compute>"Compare sequences"</button>
            <p class="compare-status" role="status">{move||status.get()}</p>
            {move||report.get().map(|bound|{
                let result=bound.alignment;let stats=result.stats.clone();let preview=aligned_preview(&result);
                view!{<div class="compare-result"><h3>{format!("Score {}",result.score)}</h3><p class="compare-results-meta">{format!("{} → {}",bound.reference_name,bound.query_name)}</p><p>{format!("Reference [{}, {}); query [{}, {}).",result.reference_start,result.reference_end,result.query_start,result.query_end)}</p><p class="helper">"Coordinates are zero-based and half-open. Empty spans mark between-base boundaries."</p><ul class="compare-stat-list"><li>{format!("Paired coverage: reference {:.1}% · query {:.1}%",stats.reference_coverage*100.,stats.query_coverage*100.)}</li><li>{format!("Certain matches {} · compatible uncertain {} · mismatches {}",stats.certain_matches,stats.ambiguous_compatible,stats.mismatches)}</li><li>{format!("Query insertions {} nt · reference deletions {} nt",stats.inserted_query_bases,stats.deleted_reference_bases)}</li></ul><details><summary>"Aligned sequence preview"</summary><p class="helper">"R: reference; Q: query. | certain, ? compatible uncertainty, × mismatch. First 4,000 columns; full result is downloadable."</p><pre>{preview}</pre></details><details><summary>"CIGAR"</summary><code>{preview_cigar(&result.cigar)}</code><p class="helper">"= certain, M compatible ambiguity, X disjoint, I/D gaps."</p></details><h4>"Differences and uncertain columns"</h4>{result.differences.iter().take(200).map(|difference|{let (start,end)=(difference.reference_start,difference.reference_end);let label=format!("{:?}: R [{}, {}), Q [{}, {})",difference.kind,start,end,difference.query_start,difference.query_end);let before=preview_text(&difference.reference_bases);let after=preview_text(&difference.query_bases);view!{<div class="compare-difference"><button disabled=start==end on:click=move|_|selection.set(Some((start,end)))>{label}</button><code>{format!("R: {before}\nQ: {after}")}</code></div>}}).collect_view()}<Show when=move||{result.differences.len()>200}><p class="compare-warning">"Difference list limited to 200 entries; download the complete report."</p></Show>{result.warnings.into_iter().map(|w|view!{<p class="compare-warning">{w}</p>}).collect_view()}<p class="helper">"This view uses the supplied linear origins/orientation. It does not search circular rotations or the reverse complement automatically."</p></div>}
            })}
            <button class="button secondary full" disabled=move||report.get().is_none() on:click=export>"Download comparison JSON"</button>
        </div>
    }
}
fn preview_cigar(value: &str) -> String {
    preview(value, 4000)
}
fn preview_text(value: &str) -> String {
    if value.is_empty() {
        "∅".into()
    } else {
        preview(value, 80)
    }
}
