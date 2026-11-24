//! Revision-bound reference gel predictions and diagnostic digest ranking.
use crate::worker;
use collagene_core::{
    Molecule,
    analysis::Primer,
    gel::{
        self, Calibration, CalibrationPoint, DiagnosticOptions, DiagnosticReport, GelReport,
        LaneInput, Preparation,
    },
    protocol::{AnalysisJob, AnalysisResult, Request, Response},
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone)]
struct LaneChoice {
    key: u32,
    index: usize,
    label: String,
    preparation: Preparation,
    mass_ng: Option<f64>,
}
fn num<T: std::str::FromStr>(text: &str, label: &str) -> Result<T, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("{label} must be a valid number"))
}
fn names(text: &str) -> Vec<String> {
    text.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect()
}
fn calibration(text: &str, resolution: &str) -> Result<Calibration, String> {
    if text.len() > 8192 {
        return Err("Ladder exceeds its text limit".into());
    }
    let mut points = vec![];
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let words = line
            .split(|c: char| c == ',' || c.is_ascii_whitespace())
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>();
        if words.len() != 2 {
            return Err("Each ladder line needs size in bp, then normalized migration".into());
        }
        points.push(CalibrationPoint {
            size_bp: num(words[0], "Ladder size")?,
            migration: num(words[1], "Ladder migration")?,
        });
    }
    let result = Calibration {
        label: "User-supplied illustrative/calibration points; no run conditions inferred".into(),
        points,
        resolution: num(resolution, "Migration resolution")?,
    };
    result.validate().map_err(|e| e.to_string())?;
    Ok(result)
}
fn short(s: &str) -> String {
    s.chars().take(12).collect()
}
#[component]
pub fn GelPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
) -> impl IntoView {
    let source = RwSignal::new(0usize);
    let prep = RwSignal::new("digest".to_string());
    let enzymes = RwSignal::new("EcoRI".to_string());
    let mass = RwSignal::new(String::new());
    let f = RwSignal::new(String::new());
    let r = RwSignal::new(String::new());
    let ft = RwSignal::new(String::new());
    let rt = RwSignal::new(String::new());
    let ladder = RwSignal::new("100, 0.9\n1000, 0.5\n10000, 0.1".to_string());
    let resolution = RwSignal::new("0.02".to_string());
    let lanes = RwSignal::new(Vec::<LaneChoice>::new());
    let next_key = RwSignal::new(0u32);
    let diagnostic_indices = RwSignal::new(Vec::<usize>::new());
    let diagnostic_enzymes =
        RwSignal::new("EcoRI,BamHI,HindIII,EcoRV,SmaI,NotI,XbaI,BsaI,PstI".to_string());
    let double = RwSignal::new(true);
    let report = RwSignal::new(None::<GelReport>);
    let diagnostics = RwSignal::new(None::<DiagnosticReport>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    let binding = Memo::new(move |_| {
        let active = molecule.get().revision();
        let choices = lanes.get();
        let candidates = diagnostic_indices.get();
        let revisions = choices
            .iter()
            .map(|l| (l.key, l.index, get_query.run(l.index).map(|m| m.revision())))
            .collect::<Vec<_>>();
        let candidate_revisions = candidates
            .iter()
            .map(|i| (*i, get_query.run(*i).map(|m| m.revision())))
            .collect::<Vec<_>>();
        (active, revisions, candidate_revisions)
    });
    Effect::new(move |_| {
        let _ = (
            binding.get(),
            lanes.get(),
            ladder.get(),
            resolution.get(),
            diagnostic_enzymes.get(),
            double.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        diagnostics.set(None);
        status.set(String::new());
        busy.set(false);
    });
    let run = {
        let alive = alive.clone();
        move |job: AnalysisJob, expected: Vec<String>| {
            let revision = molecule.get_untracked().revision();
            let bound = binding.get_untracked();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Computing reference predictions in background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
                if !alive.load(Ordering::Relaxed)
                    || token != generation.get_untracked()
                    || binding.get_untracked() != bound
                {
                    return;
                }
                busy.set(false);
                match result {
                    Ok(Response {
                        id,
                        revision: returned,
                        result,
                    }) if id == token && returned == revision => match result {
                        AnalysisResult::Gel(result)
                            if result
                                .lanes
                                .iter()
                                .map(|l| l.source_revision.clone())
                                .collect::<Vec<_>>()
                                == expected =>
                        {
                            report.set(Some(result));
                            status.set("Gel prediction complete under the recorded calibration and cleavage assumptions.".into());
                        }
                        AnalysisResult::DiagnosticDigest(result)
                            if result.source_revisions == expected =>
                        {
                            diagnostics.set(Some(result));
                            status.set("Diagnostic digest ranking complete; inspect unknown pairs and search completeness.".into());
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Worker source bindings did not match the request".into()),
                    },
                    Err(error) => status.set(error),
                    _ => status.set("Unexpected gel worker response".into()),
                }
            });
        }
    };
    let add = move |_| {
        let Some(m) = get_query.run(source.get_untracked()) else {
            status.set("Choose an open source document".into());
            return;
        };
        if lanes.with_untracked(|l| l.len() >= 32) {
            status.set("At most 32 gel lanes".into());
            return;
        }
        let mass_ng = if mass.get_untracked().trim().is_empty() {
            None
        } else {
            match num::<f64>(&mass.get_untracked(), "Loaded mass") {
                Ok(v) if v.is_finite() && v > 0. && v <= 1e9 => Some(v),
                _ => {
                    status.set(
                        "Loaded mass must be positive and finite, or blank for unknown".into(),
                    );
                    return;
                }
            }
        };
        let preparation = match prep.get_untracked().as_str() {
            "uncut" => Preparation::Uncut,
            "pcr" => match Primer::new("F", &f.get_untracked(), &ft.get_untracked()).and_then(|f| {
                Primer::new("R", &r.get_untracked(), &rt.get_untracked()).map(|r| (f, r))
            }) {
                Ok((forward, reverse)) => Preparation::Pcr { forward, reverse },
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            },
            _ => Preparation::Digest {
                enzymes: names(&enzymes.get_untracked()),
            },
        };
        let key = next_key.get_untracked();
        next_key.set(key.wrapping_add(1));
        let label = m.name.chars().take(48).collect::<String>();
        lanes.update(|lanes| {
            lanes.push(LaneChoice {
                key,
                index: source.get_untracked(),
                label,
                preparation,
                mass_ng,
            })
        });
    };
    let predict = {
        let run = run.clone();
        move |_| {
            let calibration =
                match calibration(&ladder.get_untracked(), &resolution.get_untracked()) {
                    Ok(c) => c,
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                };
            let inputs = lanes
                .get_untracked()
                .into_iter()
                .map(|l| {
                    get_query
                        .run(l.index)
                        .ok_or_else(|| "A lane source is no longer open".to_string())
                        .map(|m| LaneInput {
                            label: l.label,
                            molecule: m,
                            preparation: l.preparation,
                            mass_ng: l.mass_ng,
                        })
                })
                .collect::<Result<Vec<_>, _>>();
            match inputs {
                Ok(lanes) => {
                    let revisions = lanes.iter().map(|l| l.molecule.revision()).collect();
                    run(
                        AnalysisJob::PredictGel {
                            molecule: molecule.get_untracked(),
                            lanes,
                            calibration,
                        },
                        revisions,
                    );
                }
                Err(error) => status.set(error),
            }
        }
    };
    let diagnose = move |_| {
        let calibration = match calibration(&ladder.get_untracked(), &resolution.get_untracked()) {
            Ok(c) => c,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let candidates = diagnostic_indices
            .get_untracked()
            .into_iter()
            .map(|i| {
                get_query
                    .run(i)
                    .ok_or_else(|| "A diagnostic source is no longer open".to_string())
            })
            .collect::<Result<Vec<_>, _>>();
        match candidates {
            Ok(candidates) => {
                let revisions = candidates.iter().map(Molecule::revision).collect();
                run(
                    AnalysisJob::DiagnosticDigests {
                        molecule: molecule.get_untracked(),
                        candidates,
                        options: DiagnosticOptions {
                            enzymes: names(&diagnostic_enzymes.get_untracked()),
                            include_double_digests: double.get_untracked(),
                            calibration,
                            ..Default::default()
                        },
                    },
                    revisions,
                );
            }
            Err(error) => status.set(error),
        }
    };
    view! {
        <style>".gel-panel textarea{font:11px/1.6 monospace}.gel-panel .gel-lane{padding:8px 0;border-bottom:1px solid #e1e7e0;font-size:12px;overflow-wrap:anywhere}.gel-panel .gel-check{display:flex;align-items:center;gap:8px}.gel-panel .gel-check input{width:auto}.gel-panel .gel-source-picker{max-height:180px;overflow:auto}.gel-panel .gel-svg{overflow:auto;max-height:440px;border:1px solid #dce5de}.gel-panel .gel-table-wrap{overflow:auto}.gel-panel table{border-collapse:collapse;font-size:11px;white-space:nowrap}.gel-panel th,.gel-panel td{padding:5px;border-bottom:1px solid #e1e7e0;text-align:left}.gel-panel .gel-candidate{padding:10px 0;border-bottom:1px solid #e1e7e0;font-size:12px;line-height:1.6}.gel-panel ul{padding-left:18px}.gel-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="inspector-content gel-panel"><h2>"Predicted gel"</h2><p class="helper">"Reference predictions only. Circular conformation and out-of-ladder fragments remain unplaced. No observed-gel comparison."</p>
            <details><summary>"Ladder and resolution"</summary><label>"Size bp, migration fraction (one point per line)"<textarea aria-label="Gel ladder calibration" rows="4" prop:value=move||ladder.get() on:input=move|e|ladder.set(event_target_value(&e))/></label><label>"Migration resolution (0–1 scale)"<input aria-label="Gel migration resolution" type="number" step="0.001" prop:value=move||resolution.get() on:input=move|e|resolution.set(event_target_value(&e))/></label><p class="helper">"Default ladder is illustrative, not an empirical gel calibration. Positions interpolate versus log size only between supplied points."</p></details>
            <h3>"Add a source lane"</h3><label>"Source document"<select aria-label="Gel source document" prop:value=move||source.get().to_string() on:change=move|e|source.set(event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <label>"Preparation"<select aria-label="Gel preparation" prop:value=move||prep.get() on:change=move|e|prep.set(event_target_value(&e))><option value="digest">"Complete restriction digest"</option><option value="uncut">"Uncut reference"</option><option value="pcr">"Possible exact-PCR products"</option></select></label>
            <Show when=move||prep.get()=="digest"><label>"Enzymes (comma separated)"<input aria-label="Gel digest enzymes" prop:value=move||enzymes.get() on:input=move|e|enzymes.set(event_target_value(&e))/></label></Show>
            <Show when=move||prep.get()=="pcr">{[("Forward annealing","Gel forward primer",f),("Forward tail","Gel forward tail",ft),("Reverse annealing","Gel reverse primer",r),("Reverse tail","Gel reverse tail",rt)].into_iter().map(|(name,aria,value)|view!{<label>{name}<input aria-label=aria prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</Show>
            <label>"Loaded DNA mass (ng; blank means unknown)"<input aria-label="Gel loaded mass" type="number" prop:value=move||mass.get() on:input=move|e|mass.set(event_target_value(&e))/></label><button class="button secondary full" on:click=add>"Add gel lane"</button>
            <div class="gel-lanes"><For each=move||lanes.get() key=|l|l.key children=move|lane|{let key=lane.key;view!{<div class="gel-lane"><strong>{lane.label}</strong><p>{format!("{:?} · mass {}",lane.preparation,lane.mass_ng.map(|m|format!("{m} ng")).unwrap_or_else(||"unknown".into()))}</p><button class="text-tool" on:click=move|_|lanes.update(|l|l.retain(|l|l.key!=key))>"Remove lane"</button></div>}}/></div>
            <button class="button primary full" disabled=move||busy.get()||lanes.with(|l|l.is_empty()) on:click=predict>"Predict gel lanes"</button><p class="gel-status" role="status">{move||status.get()}</p>
            {move||report.get().map(|result|{let svg=gel::render_svg(&result);let downloadable=result.clone();let for_svg=result.clone();view!{<section class="gel-result"><h3>"Predicted lanes"</h3>{match svg{Ok(svg)=>Some(view!{<div class="gel-svg" inner_html=svg/>}),Err(_)=>None}}
                <button class="button secondary full" on:click=move|_|match gel::render_svg(&for_svg){Ok(svg)=>if let Err(e)=crate::app::download("predicted-gel.svg","image/svg+xml",svg.as_bytes()){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download gel SVG"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&downloadable){Ok(bytes)=>if let Err(e)=crate::app::download("predicted-gel.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download gel report"</button>
                {result.lanes.into_iter().map(|lane|view!{<details class="gel-lane-report"><summary>{format!("{} · {} bands · {} unplaced",lane.label,lane.bands.len(),lane.unpositioned_fragments.len())}</summary><p>{format!("{} · source {}…",lane.preparation,short(&lane.source_revision))}</p>{lane.error.map(|e|view!{<p role="alert">{e}</p>})}<ul>{lane.warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><div class="gel-table-wrap"><table><thead><tr><th>"Fragment"</th><th>"Mean bp"</th><th>"Strands nt"</th><th>"Migration"</th><th>"Mass ng"</th></tr></thead><tbody>{lane.fragments.iter().take(200).enumerate().map(|(i,f)|view!{<tr><td>{i+1}</td><td>{format!("{:.2}",f.size_bp)}</td><td>{format!("{} / {}",f.top_nt,f.bottom_nt)}</td><td>{f.migration.map(|p|format!("{p:.4}")).unwrap_or_else(||format!("Unplaced {:?}",f.state))}</td><td>{f.amount_ng.map(|m|format!("{m:.3}")).unwrap_or_else(||"unknown".into())}</td></tr>}).collect_view()}</tbody></table></div><p class="helper">{format!("Showing {} of {} fragments; JSON retains the full returned lane.",lane.fragments.len().min(200),lane.fragments.len())}</p></details>}).collect_view()}
                <details><summary>"Prediction assumptions"</summary><ul>{result.assumptions.into_iter().map(|a|view!{<li>{a}</li>}).collect_view()}</ul></details></section>}})}
            <h3>"Find diagnostic digests"</h3><p class="helper">"Choose 2–16 candidate documents. Unknown circular/unsupported/out-of-ladder profiles never count as identification evidence."</p>
            <div class="gel-source-picker">{move||query_options.get().into_iter().map(|(i,name)|view!{<label class="gel-check"><input type="checkbox" aria-label=format!("Diagnostic candidate {i}") prop:checked=move||diagnostic_indices.with(|v|v.contains(&i)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();diagnostic_indices.update(|v|{if checked&&!v.contains(&i){v.push(i);v.sort();}else if !checked{v.retain(|n|*n!=i);}});}/>{name}</label>}).collect_view()}</div>
            <label>"Candidate enzymes"<input aria-label="Diagnostic digest enzymes" prop:value=move||diagnostic_enzymes.get() on:input=move|e|diagnostic_enzymes.set(event_target_value(&e))/></label><label class="gel-check"><input aria-label="Diagnostic double digests" type="checkbox" prop:checked=move||double.get() on:change=move|e|double.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Include double digests"</label>
            <button class="button primary full" disabled=move||busy.get()||diagnostic_indices.with(|v|v.len()<2||v.len()>16) on:click=diagnose>"Rank diagnostic digests"</button>
            {move||diagnostics.get().map(|result|{let downloaded=result.clone();let source_indices=diagnostic_indices.get_untracked();let source_labels=result.source_revisions.iter().enumerate().map(|(i,revision)|format!("{}: {} · {}…",i,source_indices.get(i).and_then(|i|get_query.run(*i)).map(|m|m.name).unwrap_or_else(||"missing source".into()),short(revision))).collect::<Vec<_>>();view!{<section class="gel-diagnostics"><h3>{format!("{} retained digest candidates",result.candidates.len())}</h3><ol class="gel-diagnostic-bindings" start="0">{source_labels.into_iter().map(|label|view!{<li>{label}</li>}).collect_view()}</ol><p>{format!("{} enzyme sets evaluated · exhaustive requested search: {}",result.evaluated_digests,result.exhaustive_search)}</p><ul>{result.warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>{result.candidates.into_iter().map(|candidate|{let names=candidate.enzymes.clone();let indices=source_indices.clone();view!{<article class="gel-candidate"><strong>{candidate.enzymes.join(" + ")}</strong><p>{format!("Distinguished {} pairs · indistinguishable {} · unknown {}",candidate.distinguished_pairs.len(),candidate.indistinguishable_pairs.len(),candidate.unknown_pairs.len())}</p><p>{format!("Candidate-index pairs (zero-based): {:?}",candidate.distinguished_pairs)}</p><button class="button secondary" on:click=move|_|{let mut choices=vec![];for index in &indices{if let Some(m)=get_query.run(*index){let key=next_key.get_untracked();next_key.set(key.wrapping_add(1));choices.push(LaneChoice{key,index:*index,label:m.name.chars().take(48).collect(),preparation:Preparation::Digest{enzymes:names.clone()},mass_ng:None});}}lanes.set(choices);}>"Use digest as gel lanes"</button></article>}}).collect_view()}<button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&downloaded){Ok(bytes)=>if let Err(e)=crate::app::download("diagnostic-digests.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download digest ranking"</button></section>}})}
        </div>
    }
}
