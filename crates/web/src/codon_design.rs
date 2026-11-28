//! Source-bound codon design with a separate independent candidate-check job.
use crate::worker::CachedWorker;
use collagene_core::{
    codon_design::{
        CheckReport, CodonWeight, DesignOptions, DesignReport, DesignStatus, GcBounds, WindowGc,
    },
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule, Topology};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

fn whole<T: std::str::FromStr>(text: &str, label: &str) -> Result<T, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {label}."))
}
fn percent(text: &str) -> Result<u16, String> {
    let number = text
        .trim()
        .parse::<f64>()
        .map_err(|_| "Enter GC percentages between 0 and 100.".to_string())?;
    if !number.is_finite()
        || !(0.0..=100.0).contains(&number)
        || ((number * 100.0).round() - number * 100.0).abs() > 0.000001
    {
        return Err("GC percentages accept at most two decimal places between 0 and 100.".into());
    }
    Ok((number * 100.0).round() as u16)
}
fn bases(text: &str) -> String {
    text.chars()
        .filter(|c| !c.is_ascii_whitespace())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}
fn preview(text: &str) -> String {
    let mut shown = text.chars().take(600).collect::<String>();
    if text.len() > 600 {
        shown.push_str("…\nPreview limited to 600 symbols; exports retain the complete sequence.");
    }
    shown
}
struct JsonOutput(Vec<u8>);
impl std::io::Write for JsonOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > 8 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "Codon report exceeds this browser's 8 MiB export limit",
            ));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
async fn compute(
    worker: &CachedWorker,
    source: &Molecule,
    options: &DesignOptions,
    supplied: Option<Molecule>,
    id: u32,
) -> Result<(Option<DesignReport>, Option<Molecule>, Option<CheckReport>), String> {
    let revision = source.revision();
    let (design, candidate) = if let Some(mut candidate) = supplied {
        // Bind portable evidence before checking, so the checker signs the exact
        // molecule revision later opened or exported. This metadata is not authority.
        if source.sequence.len() > collagene_core::codon_design::MAX_RESIDUES * 3 {
            return Err("Codon source exceeds this profile's coding-sequence bound.".into());
        }
        let protein = match source.sequence.alphabet() {
            Alphabet::Protein => source.sequence.as_str().to_string(),
            Alphabet::Dna => {
                collagene_core::translation::translate(
                    &source.sequence,
                    &collagene_core::translation::TranslationOptions {
                        genetic_code: options.genetic_code,
                        ..Default::default()
                    },
                )
                .map_err(|e| e.to_string())?
                .protein
            }
            Alphabet::Rna => return Err("Choose a protein or concrete DNA coding source.".into()),
        };
        candidate.metadata.insert(
            "codon_design:engine".into(),
            collagene_core::codon_design::CODON_ENGINE_VERSION.into(),
        );
        candidate
            .metadata
            .insert("codon_design:source_revision".into(), revision.clone());
        candidate.metadata.insert(
            "codon_design:options_revision".into(),
            collagene_core::codon_design::options_revision(options).map_err(|e| e.to_string())?,
        );
        candidate.metadata.insert(
            "codon_design:options".into(),
            serde_json::to_string(options).map_err(|e| e.to_string())?,
        );
        candidate
            .metadata
            .insert("codon_design:required_protein".into(), protein);
        candidate.metadata.insert(
            "codon_design:origin".into(),
            "User-supplied DNA checked against an independent source; not optimizer-produced"
                .into(),
        );
        (None, Some(candidate))
    } else {
        let response = worker
            .ask(Request {
                id,
                revision: revision.clone(),
                job: AnalysisJob::DesignCodons {
                    molecule: source.clone(),
                    options: options.clone(),
                },
            })
            .await?;
        let report = match response.result {
            AnalysisResult::CodonDesign(report) => report,
            AnalysisResult::Error(error) => return Err(error),
            _ => return Err("Unexpected codon design response.".into()),
        };
        report
            .validate_bindings_for_source(source, options)
            .map_err(|e| e.to_string())?;
        let candidate = report.candidate.clone();
        (Some(report), candidate)
    };
    let Some(candidate) = candidate else {
        return Ok((design, None, None));
    };
    // This is a separate worker request, not the producer's embedded candidate_check claim.
    let response = worker
        .ask(Request {
            id: id.wrapping_add(1),
            revision,
            job: AnalysisJob::CheckCodonDesign {
                source: source.clone(),
                candidate: candidate.clone(),
                options: options.clone(),
            },
        })
        .await?;
    let check = match response.result {
        AnalysisResult::CodonCheck(check) => check,
        AnalysisResult::Error(error) => return Err(error),
        _ => return Err("Unexpected codon checking response.".into()),
    };
    check
        .validate_bindings_for(source, &candidate, options)
        .map_err(|e| e.to_string())?;
    Ok((design, Some(candidate), Some(check)))
}

#[component]
pub fn CodonDesignPanel(
    molecule: Memo<Molecule>,
    source_identity: Memo<String>,
    epoch: RwSignal<u64>,
    on_open: Callback<Molecule, Result<(), String>>,
) -> impl IntoView {
    let code = RwSignal::new(1u8);
    let weights = RwSignal::new(String::new());
    let motifs = RwSignal::new(String::new());
    let left = RwSignal::new(String::new());
    let right = RwSignal::new(String::new());
    let min_gc = RwSignal::new("0".to_string());
    let max_gc = RwSignal::new("100".to_string());
    let homopolymer = RwSignal::new("12".to_string());
    let window_enabled = RwSignal::new(false);
    let window = RwSignal::new("30".to_string());
    let window_min = RwSignal::new("0".to_string());
    let window_max = RwSignal::new("100".to_string());
    let beam = RwSignal::new("128".to_string());
    let work = RwSignal::new(DesignOptions::default().max_work.to_string());
    let supplied = RwSignal::new(String::new());
    let options = Memo::new(move |_| {
        if weights.with(|s| s.len() > 4096)
            || motifs.with(|s| s.len() > 16384)
            || left.with(|s| s.len() > 8192)
            || right.with(|s| s.len() > 8192)
        {
            return Err("Codon options exceed the browser text input limit.".into());
        }
        let mut codon_weights = vec![];
        for line in weights.get().lines().filter(|line| !line.trim().is_empty()) {
            let fields = line
                .split(|c: char| c.is_ascii_whitespace() || c == ',' || c == '=')
                .filter(|v| !v.is_empty())
                .collect::<Vec<_>>();
            if fields.len() != 2 {
                return Err("Use one codon and integer weight per line, for example AAG 100. Unspecified codons have weight 1; zero excludes a codon.".into());
            }
            codon_weights.push(CodonWeight {
                codon: fields[0].to_ascii_uppercase(),
                weight: whole(fields[1], "codon weight")?,
            });
        }
        let o = DesignOptions {
            genetic_code: code.get(),
            codon_weights,
            forbidden_motifs: motifs
                .get()
                .split(|c: char| c.is_ascii_whitespace() || c == ',')
                .filter(|v| !v.is_empty())
                .map(|v| v.to_ascii_uppercase())
                .collect(),
            left_context: bases(&left.get()),
            right_context: bases(&right.get()),
            gc: GcBounds {
                minimum_bp: percent(&min_gc.get())?,
                maximum_bp: percent(&max_gc.get())?,
            },
            max_homopolymer: whole(&homopolymer.get(), "maximum homopolymer length")?,
            window_gc: if window_enabled.get() {
                Some(WindowGc {
                    window: whole(&window.get(), "GC window length")?,
                    bounds: GcBounds {
                        minimum_bp: percent(&window_min.get())?,
                        maximum_bp: percent(&window_max.get())?,
                    },
                })
            } else {
                None
            },
            beam_width: whole(&beam.get(), "retained search states")?,
            max_work: whole(&work.get(), "search work limit")?,
        };
        o.validate().map_err(|e| e.to_string())?;
        Ok::<DesignOptions, String>(o)
    });
    let revision = Memo::new(move |_| molecule.with(Molecule::revision));
    let design = RwSignal::new(None::<DesignReport>);
    let candidate = RwSignal::new(None::<Molecule>);
    let checked = RwSignal::new(None::<CheckReport>);
    let bound = RwSignal::new(None::<(String, String, u64, DesignOptions, String)>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (
            revision.get(),
            source_identity.get(),
            epoch.get(),
            options.get(),
            supplied.get(),
        );
        generation.update(|g| *g = g.wrapping_add(2));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        busy.set(false);
        design.set(None);
        candidate.set(None);
        checked.set(None);
        bound.set(None);
        status.set(String::new());
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = active.get_value() {
                worker.close();
            }
            active.set_value(None);
        });
    }
    let run = {
        let alive = alive.clone();
        Callback::new(move |check_supplied: bool| {
            let opts = match options.get_untracked() {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let supplied_text = supplied.get_untracked();
            let supplied = if check_supplied {
                if supplied_text.len() > 120_000 {
                    status.set(
                        "Supplied DNA exceeds this browser's 120,000-character input bound.".into(),
                    );
                    return;
                }
                match Molecule::dna("Checked codon candidate", &supplied_text, Topology::Linear) {
                    Ok(v) => Some(v),
                    Err(e) => {
                        status.set(e.to_string());
                        return;
                    }
                }
            } else {
                None
            };
            if let Some(worker) = active.get_value() {
                worker.close();
            }
            let worker = match CachedWorker::new() {
                Ok(v) => Rc::new(v),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            active.set_value(Some(worker.clone()));
            let source = molecule.get_untracked();
            let source_revision = revision.get_untracked();
            let identity = source_identity.get_untracked();
            let source_epoch = epoch.get_untracked();
            generation.update(|g| *g = g.wrapping_add(2));
            let token = generation.get_untracked();
            design.set(None);
            candidate.set(None);
            checked.set(None);
            bound.set(None);
            busy.set(true);
            status.set(if check_supplied {
                "Independently checking the supplied DNA against this source and its constraints…"
                    .into()
            } else {
                "Searching synonymous DNA candidates, then checking the result independently…"
                    .into()
            });
            let alive = alive.clone();
            spawn_local(async move {
                let result = compute(&worker, &source, &opts, supplied, token).await;
                worker.close();
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || revision.get_untracked() != source_revision
                    || source_identity.get_untracked() != identity
                    || epoch.get_untracked() != source_epoch
                    || options.get_untracked() != Ok(opts.clone())
                {
                    return;
                }
                busy.set(false);
                active.set_value(None);
                match result {
                    Ok((d, c, r)) => {
                        bound.set(Some((
                            identity,
                            source_revision,
                            source_epoch,
                            opts,
                            supplied_text,
                        )));
                        status.set(match r.as_ref(){Some(r)if r.complete&&r.satisfies==Some(true)=>"Candidate passed a separate source-bound translation and constraint check. This is sequence design evidence, not experimental validation.".into(),Some(_)=>"Candidate did not pass all required checks; inspect the violations and limits.".into(),None=>match d.as_ref().map(|v|v.status){Some(DesignStatus::Unsatisfiable)=>"No valid synonymous sequence exists under the declared constraints within this exhaustive search.".into(),_=>"No candidate was returned before the declared search limits. This does not establish impossibility.".into()}});
                        design.set(d);
                        candidate.set(c);
                        checked.set(r);
                    }
                    Err(e) => status.set(e),
                }
            });
        })
    };
    let approved = Memo::new(move |_| {
        checked.with(|r| {
            r.as_ref()
                .is_some_and(|r| r.complete && r.satisfies == Some(true))
        })
    });
    let can_run = Memo::new(move |_| {
        matches!(
            molecule.with(|m| m.sequence.alphabet()),
            Alphabet::Dna | Alphabet::Protein
        ) && options.get().is_ok()
    });
    let guard = Callback::new(move |()| -> Result<(), String> {
        let Some((identity, source_revision, source_epoch, opts, supplied_text)) =
            bound.get_untracked()
        else {
            return Err("Run codon analysis for this source before using the result.".into());
        };
        if source_identity.get_untracked() != identity
            || revision.get_untracked() != source_revision
            || epoch.get_untracked() != source_epoch
            || options.get_untracked() != Ok(opts)
            || supplied.get_untracked() != supplied_text
        {
            return Err(
                "Source or constraints changed. Run codon analysis again before using this result."
                    .into(),
            );
        }
        Ok(())
    });
    let open = Callback::new(move |()| {
        if let Err(e) = guard.run(()) {
            status.set(e);
            return;
        }
        if !approved.get_untracked() {
            return;
        }
        if let Some(value) = candidate.get_untracked() {
            if let Err(e) = on_open.run(value) {
                status.set(e);
            }
        }
    });
    let export = Callback::new(move |()| {
        if let Err(e) = guard.run(()) {
            status.set(e);
            return;
        }
        if !approved.get_untracked() {
            return;
        }
        if let Some(value) = candidate.get_untracked() {
            match collagene_core::formats::export_genbank(&value) {
                Ok(text) => crate::app::download_with_completion(
                    "codon-designed-dna.gb",
                    "text/plain",
                    text.as_bytes(),
                    move |result| {
                        if status.try_get_untracked().is_none() {
                            return;
                        }
                        if let Err(e) = result {
                            status.set(e)
                        }
                    },
                ),
                Err(e) => status.set(e.to_string()),
            }
        }
    });
    let report_download = Callback::new(move |()| {
        if let Err(e) = guard.run(()) {
            status.set(e);
            return;
        }
        let mut bytes = JsonOutput(vec![]);
        match serde_json::to_writer_pretty(
            &mut bytes,
            &serde_json::json!({"design":design.get_untracked(),"independent_check":checked.get_untracked()}),
        ) {
            Ok(()) => crate::app::download_with_completion(
                "codon-design-report.json",
                "application/json",
                &bytes.0,
                move |result| {
                    if status.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        status.set(e)
                    }
                },
            ),
            Err(e) => status.set(e.to_string()),
        }
    });
    view! {<section class="inspector-content codon-design-panel" aria-label="Codon design">
        <h2>"Codon design"</h2>
        <p class="helper">"Reverse translate a protein or redesign a complete, concrete DNA coding sequence. The required protein comes from the active source. A separate check must pass before the DNA can be opened or exported."</p>
        <p class="helper">"Translation uses internal-codon meanings, including a single terminal stop when present. Initiation, expression and synthesis success are not inferred. Context bases are fixed constraints and are not included in the returned coding sequence."</p>
        <label>"Genetic code"<select aria-label="Codon genetic code" prop:value=move||code.get().to_string() on:change=move|e|code.set(event_target_value(&e).parse().unwrap_or(1))>{collagene_core::translation::genetic_codes().into_iter().map(|c|view!{<option value=c.id.to_string()>{format!("{} · {}",c.id,c.name)}</option>}).collect_view()}</select></label>
        <details><summary>"Codon preferences"</summary><p class="helper">"Optional integer weights, one codon per line. Unspecified codons have weight 1; zero excludes a codon. Weights are a linear preference score, not CAI or an organism-specific expression model."</p><textarea aria-label="Codon weights" rows="4" placeholder="AAG 100\nAAA 1" prop:value=move||weights.get() on:input=move|e|weights.set(event_target_value(&e))></textarea></details>
        <label>"Forbidden DNA motifs"<textarea aria-label="Codon forbidden motifs" rows="3" placeholder="GAATTC\nGGTCTC" prop:value=move||motifs.get() on:input=move|e|motifs.set(event_target_value(&e))></textarea></label><p class="helper">"Separate motifs with spaces, commas or newlines. IUPAC motifs are checked on both strands, including context junctions."</p>
        <div class="codon-grid"><label>"Minimum coding GC (%)"<input aria-label="Codon minimum GC" type="number" min="0" max="100" step="0.01" prop:value=move||min_gc.get() on:input=move|e|min_gc.set(event_target_value(&e))/></label><label>"Maximum coding GC (%)"<input aria-label="Codon maximum GC" type="number" min="0" max="100" step="0.01" prop:value=move||max_gc.get() on:input=move|e|max_gc.set(event_target_value(&e))/></label></div>
        <label>"Maximum homopolymer length"<input aria-label="Codon maximum homopolymer" type="number" min="1" step="1" prop:value=move||homopolymer.get() on:input=move|e|homopolymer.set(event_target_value(&e))/></label>
        <details><summary>"Fixed contexts and local GC"</summary><label>"Left DNA context"<textarea aria-label="Codon left context" rows="2" prop:value=move||left.get() on:input=move|e|left.set(event_target_value(&e))></textarea></label><label>"Right DNA context"<textarea aria-label="Codon right context" rows="2" prop:value=move||right.get() on:input=move|e|right.set(event_target_value(&e))></textarea></label><label class="checkbox-label"><input aria-label="Codon window GC enabled" type="checkbox" prop:checked=move||window_enabled.get() on:change=move|e|window_enabled.set(event_target_checked(&e))/>"Check local GC across coding sequence and context"</label><Show when=move||window_enabled.get()><label>"Window length (bases)"<input aria-label="Codon GC window" type="number" min="1" step="1" prop:value=move||window.get() on:input=move|e|window.set(event_target_value(&e))/></label><div class="codon-grid"><label>"Minimum window GC (%)"<input aria-label="Codon window minimum GC" type="number" min="0" max="100" step="0.01" prop:value=move||window_min.get() on:input=move|e|window_min.set(event_target_value(&e))/></label><label>"Maximum window GC (%)"<input aria-label="Codon window maximum GC" type="number" min="0" max="100" step="0.01" prop:value=move||window_max.get() on:input=move|e|window_max.set(event_target_value(&e))/></label></div></Show></details>
        <p class="helper">"Only complete GC windows are evaluated; a construct shorter than the window contains no complete window."</p><details><summary>"Search limits"</summary><label>"Retained search states"<input aria-label="Codon beam width" type="number" min="1" step="1" prop:value=move||beam.get() on:input=move|e|beam.set(event_target_value(&e))/></label><label>"Maximum work"<input aria-label="Codon work limit" type="number" min="1" step="1" prop:value=move||work.get() on:input=move|e|work.set(event_target_value(&e))/></label><p class="helper">"Pruning or work exhaustion prevents a proof of global optimality or impossibility. Every returned candidate still requires its own complete check."</p></details>
        {move||options.get().err().map(|error|view!{<p role="alert" class="helper">{error}</p>})}
        <button class="button primary full" disabled=move||busy.get()||!can_run.get() on:click=move|_|run.run(false)>"Design coding DNA"</button>
        <details><summary>"Check an existing DNA candidate"</summary><textarea aria-label="Codon supplied DNA" rows="4" maxlength="120000" prop:value=move||supplied.get() on:input=move|e|supplied.set(event_target_value(&e))></textarea><button class="button secondary full" disabled=move||busy.get()||!can_run.get()||supplied.with(|s|s.trim().is_empty()) on:click=move|_|run.run(true)>"Check coding DNA"</button></details>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(2));if let Some(worker)=active.get_value(){worker.close();}active.set_value(None);busy.set(false);status.set("Codon design cancelled.".into());}>"Cancel codon design"</button></Show>
        <p class="helper codon-status" role="status">{move||status.get()}</p>
        {move||design.get().map(|report|view!{<div class="codon-design-result" data-status=format!("{:?}",report.status) data-search-complete=report.search_complete.to_string() data-optimality=report.optimality_proven.to_string()><h3>{format!("{:?}",report.status)}</h3><p>{format!("Search complete: {} · optimality proven: {} · work {} · pruned states {}",report.search_complete,report.optimality_proven,report.work_used,report.pruned_states)}</p><p>{format!("Preference score: {}",report.score.map(|s|s.to_string()).unwrap_or_else(||"No candidate".into()))}</p><h4>"Required protein"</h4><pre class="codon-target-protein">{preview(&report.target_protein)}</pre><details><summary>"Source and model"</summary><p>{format!("{} · source {} · options {}",report.engine_version,report.source_revision,report.options_revision)}</p></details>{report.limits.into_iter().chain(report.warnings).take(100).map(|text|view!{<p class="helper">{text}</p>}).collect_view()}</div>})}
        {move||checked.get().map(|report|view!{<div class="codon-check-result" data-complete=report.complete.to_string() data-satisfies=report.satisfies.map(|v|v.to_string()).unwrap_or_else(||"unknown".into())><h3>"Independent candidate check"</h3><p>{format!("Complete: {} · satisfies constraints: {}",report.complete,report.satisfies.map(|v|v.to_string()).unwrap_or_else(||"unknown".into()))}</p><p>{format!("{} coding bases · {} GC bases · score {}",report.coding_length,report.gc_count,report.score)}</p><h4>"Required protein"</h4><pre class="codon-check-target">{preview(&report.target_protein)}</pre><h4>"Candidate translation"</h4><pre class="codon-candidate-protein">{preview(&report.translated_protein)}</pre>{report.violations.into_iter().chain(report.warnings).take(100).map(|text|view!{<p class="helper">{text}</p>}).collect_view()}</div>})}
        {move||candidate.get().map(|value|view!{<details><summary>"DNA candidate preview"</summary><pre class="codon-candidate-dna">{preview(value.sequence.as_str())}</pre></details>})}
        <Show when=move||approved.get()><button class="button secondary full" on:click=move|_|open.run(())>"Open designed DNA"</button><button class="button secondary full" on:click=move|_|export.run(())>"Download designed GenBank"</button></Show>
        <Show when=move||design.get().is_some()||checked.get().is_some()><button class="button secondary full" on:click=move|_|report_download.run(())>"Download codon report"</button></Show>
    </section>}
}
