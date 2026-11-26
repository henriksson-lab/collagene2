//! Explicit background translation/ORF tools with source-event-bound outputs.
use crate::{digest::SourceInfo, worker};
use collagene_core::{
    Alphabet, Feature, Location, Molecule,
    protocol::{AnalysisJob, AnalysisResult, Request, Response},
    run_lineage::EventBinding,
    translation::{
        self, IncompletePolicy, InitiationPolicy, OrfOptions, OrfReport, StartPolicy,
        TranslationOptions,
    },
    translation_reference::{
        self, Outputs, Report as ReferenceReport, Request as TranslationRequest,
        Selection as TranslationSelection,
    },
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;
use web_sys::HtmlInputElement;

#[derive(Clone)]
pub struct ProteinOutput {
    pub source: EventBinding,
    pub epoch: u64,
    pub name: String,
    pub request: Arc<TranslationRequest>,
    pub outputs: Arc<Outputs>,
}
#[derive(Clone)]
pub struct AnnotationOutput {
    pub source: EventBinding,
    pub epoch: u64,
    pub name: String,
    pub request: Arc<TranslationRequest>,
    pub outputs: Arc<Outputs>,
}
#[derive(Clone)]
struct BoundProtein {
    source: EventBinding,
    epoch: u64,
    request: Arc<TranslationRequest>,
    report: Arc<ReferenceReport>,
}
#[derive(Clone)]
struct BoundOrfs {
    source: EventBinding,
    epoch: u64,
    report: Arc<OrfReport>,
}

fn preview(protein: &str) -> String {
    let shown: String = protein.chars().take(600).collect();
    if protein.len() > 600 {
        format!("{shown}…\nPreview limited to 600 symbols; download keeps the full protein.")
    } else {
        shown
    }
}
fn description(report: &ReferenceReport) -> String {
    if let Some(cds) = &report.cds {
        format!(
            "{} · code {} · codon_start {} · {} amino-acid symbols",
            cds.feature_label.chars().take(160).collect::<String>(),
            cds.translation.options.genetic_code,
            cds.codon_start,
            report.protein.len()
        )
    } else {
        let options = &report.translation.options;
        format!(
            "Code {} · frame {} · {} · {} amino-acid symbols",
            options.genetic_code,
            u16::from(options.frame) + 1,
            if options.reverse {
                "reverse"
            } else {
                "forward"
            },
            report.protein.len()
        )
    }
}
fn fasta_record(name: &str, protein: &str, description: &str) -> String {
    let name: String = name
        .chars()
        .map(|c| {
            if c.is_whitespace() || c.is_control() {
                '_'
            } else {
                c
            }
        })
        .collect();
    let description: String = description
        .chars()
        .map(|c| if c.is_control() { ' ' } else { c })
        .collect();
    let mut text = format!(">{name} [molecule=protein] {description}\n");
    for line in protein.as_bytes().chunks(80) {
        text.push_str(std::str::from_utf8(line).unwrap_or(""));
        text.push('\n');
    }
    text
}

#[component]
pub fn BiologyPanel(
    molecule: Memo<Molecule>,
    selected_feature: RwSignal<Option<usize>>,
    selection: RwSignal<Option<(usize, usize)>>,
    on_annotate: Callback<(String, Feature), Result<(), String>>,
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    selection_location: Callback<(), Result<Location, String>>,
    on_open_protein: Callback<ProteinOutput, Result<(), String>>,
    on_retain_cds: Callback<AnnotationOutput, Result<(), String>>,
    on_select_location: Callback<(String, Location), Result<(), String>>,
) -> impl IntoView {
    let code = RwSignal::new(1u8);
    let frame = RwSignal::new(0u8);
    let reverse = RwSignal::new(false);
    let scope_selection = RwSignal::new(false);
    let initiation = RwSignal::new(InitiationPolicy::Ordinary);
    let start_policy = RwSignal::new(StartPolicy::AtgOnly);
    let incomplete = RwSignal::new(IncompletePolicy::CompleteOnly);
    let minimum = RwSignal::new("20".to_string());
    let maximum = RwSignal::new("1000".to_string());
    let both_strands = RwSignal::new(true);
    let output_name = RwSignal::new("Translated protein".to_string());
    let protein_page = RwSignal::new(0usize);
    let orf_page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let result = RwSignal::new(None::<Arc<BoundProtein>>);
    let orfs = RwSignal::new(None::<Arc<BoundOrfs>>);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    let current_request = Callback::new(move |cds: bool| -> Result<TranslationRequest, String> {
        let m = molecule.get_untracked();
        if m.sequence.alphabet() == Alphabet::Protein {
            return Err("Choose a nucleotide reference.".into());
        }
        let selected = if cds {
            let feature_index = selected_feature
                .get_untracked()
                .filter(|i| {
                    m.features
                        .get(*i)
                        .is_some_and(|f| f.kind.eq_ignore_ascii_case("CDS"))
                })
                .ok_or("Select a CDS annotation first.")?;
            TranslationSelection::Cds { feature_index }
        } else {
            let location = if scope_selection.get_untracked() {
                selection_location.run(())?
            } else {
                Location::range(0, m.sequence.len())
            };
            TranslationSelection::Region {
                location,
                options: TranslationOptions {
                    genetic_code: code.get_untracked(),
                    frame: frame.get_untracked(),
                    reverse: reverse.get_untracked(),
                    initiation: initiation.get_untracked(),
                },
            }
        };
        Ok(TranslationRequest {
            source: m,
            selection: selected,
        })
    });
    let guard = Callback::new(move |saved: Arc<BoundProtein>| -> Result<(), String> {
        if source.get_untracked().binding != saved.source || epoch.get_untracked() != saved.epoch {
            return Err("The source event changed; translate again.".into());
        }
        let cds = matches!(&saved.request.selection, TranslationSelection::Cds { .. });
        if current_request.run(cds)? != *saved.request {
            return Err(
                "Translation scope, location, options or selected CDS changed; translate again."
                    .into(),
            );
        }
        Ok(())
    });
    Effect::new(move |_| {
        let _ = (
            source.get().binding,
            epoch.get(),
            molecule.with(Molecule::revision),
            code.get(),
            frame.get(),
            reverse.get(),
            initiation.get(),
            start_policy.get(),
            incomplete.get(),
            minimum.get(),
            maximum.get(),
            both_strands.get(),
            scope_selection.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        result.set(None);
        orfs.set(None);
        protein_page.set(0);
        orf_page.set(0);
        status.set(String::new());
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = (selected_feature.get(), selection.get());
        if let Some(saved) = result.get_untracked() {
            if guard.run(saved).is_err() {
                generation.update(|g| *g = g.wrapping_add(1));
                result.set(None);
                busy.set(false);
                status.set("Selection changed; translate again before creating outputs.".into());
            }
        }
    });
    let run_translation = {
        let alive = alive.clone();
        Callback::new(move |cds: bool| {
            let request = match current_request.run(cds) {
                Ok(request) => Arc::new(request),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let info = source.get_untracked();
            if info.binding.revision != request.source.revision() {
                status.set("Source identity changed; translate again.".into());
                return;
            }
            let revision = match translation_reference::inspect_revision(&request) {
                Ok(value) => value,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let captured_source = info.binding;
            let captured_epoch = epoch.get_untracked();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            result.set(None);
            protein_page.set(0);
            status.set("Translating the exact source location in background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::TranslateReference {
                        request: (*request).clone(),
                    },
                })
                .await;
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                let computed = (|| -> Result<Arc<BoundProtein>, String> {
                    let Response {
                        id,
                        revision: returned,
                        result: computed,
                    } = response?;
                    if id != token || returned != revision {
                        return Err("Translation response identity changed.".into());
                    }
                    let AnalysisResult::TranslatedReference(report) = computed else {
                        return Err(match computed {
                            AnalysisResult::Error(e) => e,
                            _ => "Unexpected translation worker response.".into(),
                        });
                    };
                    report
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    let saved = Arc::new(BoundProtein {
                        source: captured_source,
                        epoch: captured_epoch,
                        request,
                        report: Arc::new(report),
                    });
                    guard.run(saved.clone())?;
                    Ok(saved)
                })();
                match computed {
                    Ok(saved) => {
                        result.set(Some(saved));
                        status.set("Translation ready. Review every protein page and the exact location before choosing outputs.".into());
                    }
                    Err(e) => status.set(e),
                }
            });
        })
    };
    let do_translate = move |_| run_translation.run(false);
    let do_cds = move |_| run_translation.run(true);
    let prepare = {
        let alive = alive.clone();
        Callback::new(move |retain: bool| {
            let Some(saved) = result.get_untracked() else {
                status.set("Translate first.".into());
                return;
            };
            if let Err(e) = guard.run(saved.clone()) {
                status.set(e);
                return;
            }
            if retain
                && !matches!(
                    &saved.request.selection,
                    TranslationSelection::Region { .. }
                )
            {
                status.set("Selected CDS translation already has its source annotation.".into());
                return;
            }
            let name = output_name.get_untracked().trim().to_string();
            let revision = match translation_reference::output_revision(&saved.request, &name) {
                Ok(value) => value,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Freshly reconstructing the named translation output…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::PrepareTranslatedReference {
                        request: (*saved.request).clone(),
                        name: name.clone(),
                    },
                })
                .await;
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                let outcome = (|| -> Result<(), String> {
                    guard.run(saved.clone())?;
                    if output_name.get_untracked().trim() != name {
                        return Err("Output name changed; request the output again.".into());
                    }
                    if !result.with_untracked(|current| {
                        current
                            .as_ref()
                            .is_some_and(|current| Arc::ptr_eq(current, &saved))
                    }) {
                        return Err("Translation review changed; translate again.".into());
                    }
                    let Response {
                        id,
                        revision: returned,
                        result: computed,
                    } = response?;
                    if id != token || returned != revision {
                        return Err("Translated output response identity changed.".into());
                    }
                    let AnalysisResult::TranslatedOutputs(outputs) = computed else {
                        return Err(match computed {
                            AnalysisResult::Error(e) => e,
                            _ => "Unexpected translated output response.".into(),
                        });
                    };
                    outputs
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    let outputs = Arc::new(outputs);
                    if retain {
                        on_retain_cds.run(AnnotationOutput {
                            source: saved.source.clone(),
                            epoch: saved.epoch,
                            name,
                            request: saved.request.clone(),
                            outputs,
                        })
                    } else {
                        on_open_protein.run(ProteinOutput {
                            source: saved.source.clone(),
                            epoch: saved.epoch,
                            name,
                            request: saved.request.clone(),
                            outputs,
                        })
                    }
                })();
                status.set(match outcome {
                    Ok(()) => {
                        if retain {
                            "Retained the translated region as a CDS requiring review.".into()
                        } else {
                            "Opened the linked protein reference.".into()
                        }
                    }
                    Err(e) => e,
                });
            });
        })
    };
    let do_orfs = {
        let alive = alive.clone();
        move |_| {
            let (Ok(min), Ok(max)) = (
                minimum.get_untracked().parse::<usize>(),
                maximum.get_untracked().parse::<usize>(),
            ) else {
                status.set("ORF lengths and result limits must be non-negative integers.".into());
                return;
            };
            let m = molecule.get_untracked();
            let revision = m.revision();
            let captured_source = source.get_untracked().binding;
            let captured_epoch = epoch.get_untracked();
            if captured_source.revision != revision {
                status.set("Source changed; search again.".into());
                return;
            }
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            orfs.set(None);
            orf_page.set(0);
            status.set("Finding ORF candidates in background…".into());
            let alive = alive.clone();
            let options = OrfOptions {
                genetic_code: code.get_untracked(),
                start_policy: start_policy.get_untracked(),
                incomplete: incomplete.get_untracked(),
                min_amino_acids: min,
                max_orfs: max,
                both_strands: both_strands.get_untracked(),
            };
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::FindOrfs {
                        molecule: m,
                        options,
                    },
                })
                .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || source.get_untracked().binding != captured_source
                    || epoch.get_untracked() != captured_epoch
                {
                    return;
                }
                busy.set(false);
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result: AnalysisResult::Orfs(report),
                    }) if id == token && returned == revision => {
                        status.set(format!(
                            "{} ORF candidate(s){}; candidates do not establish expression.",
                            report.orfs.len(),
                            if report.truncated {
                                " · result limit reached"
                            } else {
                                ""
                            }
                        ));
                        orfs.set(Some(Arc::new(BoundOrfs {
                            source: captured_source,
                            epoch: captured_epoch,
                            report: Arc::new(report),
                        })));
                    }
                    Ok(Response {
                        result: AnalysisResult::Error(e),
                        ..
                    }) => status.set(e),
                    Err(e) => status.set(e),
                    _ => status.set("ORF response does not match the request.".into()),
                }
            });
        }
    };
    view! {
        <style>".biology-panel pre{white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.6 monospace;max-height:240px;overflow:auto;background:#f4f7f4;padding:10px;border-radius:5px}.biology-panel .biology-check{display:flex;align-items:flex-start;gap:8px}.biology-panel .biology-check input{width:auto;flex-shrink:0;margin-top:3px}.biology-panel .orf-entry{padding:10px 0;border-bottom:1px solid #e1e7e0;font-size:13px;line-height:1.5}.biology-panel .orf-actions{display:flex;gap:6px;flex-wrap:wrap}.biology-panel .orf-actions button{font-size:12px;padding:6px 8px}.biology-panel details{margin:10px 0}.biology-panel ul{padding-left:20px;overflow-wrap:anywhere}.biology-status{font-size:13px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="inspector-content biology-panel">
            <h2>"Translation & ORFs"</h2>
            <Show when=move||molecule.get().sequence.alphabet()==Alphabet::Protein><p class="helper">"Choose a DNA or RNA document for nucleotide translation and ORF search."</p></Show>
            <label>"Genetic code"<select aria-label="Genetic code" prop:value=move||code.get().to_string() on:change=move|e|code.set(event_target_value(&e).parse().unwrap_or(1))>{translation::genetic_codes().into_iter().map(|c|view!{<option value=c.id.to_string()>{format!("{} · {}",c.id,c.name)}</option>}).collect_view()}</select></label>
            <h3>"Read a frame"</h3>
            <label>"Translation scope"<select aria-label="Translation scope" prop:value=move||{if scope_selection.get(){"selection"}else{"whole"}} on:change=move|e|{scope_selection.set(event_target_value(&e)=="selection");}><option value="whole">"Whole reference"</option><option value="selection">"Current selection / selected feature"</option></select></label>
            <p class="helper">"Current selection captures the exact selected location, including compound joins and circular origin crossings. Use the selected-CDS action to honor its qualifiers."</p>
            <label>"Frame offset"<select aria-label="Translation frame" prop:value=move||frame.get().to_string() on:change=move|e|frame.set(event_target_value(&e).parse().unwrap_or(0))><option value="0">"1 · skip 0 bases"</option><option value="1">"2 · skip 1 base"</option><option value="2">"3 · skip 2 bases"</option></select></label>
            <label class="biology-check"><input type="checkbox" aria-label="Reverse translation" prop:checked=move||reverse.get() on:change=move|e|reverse.set(event_target::<HtmlInputElement>(&e).checked())/>"Reverse complement"</label>
            <label>"First codon"<select aria-label="Initiation policy" on:change=move|e|{initiation.set(match event_target_value(&e).as_str(){"recognized"=>InitiationPolicy::RecognizedStartAsMethionine,"required"=>InitiationPolicy::RequireRecognizedStart,_=>InitiationPolicy::Ordinary});}><option value="ordinary">"Ordinary internal-codon translation"</option><option value="recognized">"Recognized initiator → M"</option><option value="required">"Require a recognized initiator"</option></select></label>
            <button class="button primary full" disabled=move||busy.get()||molecule.get().sequence.alphabet()==Alphabet::Protein on:click=do_translate>"Translate frame"</button>
            <h3>"Selected coding feature"</h3>
            <p class="helper">{move||selected_feature.get().and_then(|i|molecule.get().features.get(i).cloned()).filter(|f|f.kind.eq_ignore_ascii_case("CDS")).map(|f|format!("Selected CDS: {}",f.label().chars().take(160).collect::<String>())).unwrap_or_else(||"Select a CDS in the feature list. CDS phase and genetic code come from its qualifiers.".into())}</p>
            <button class="button secondary full" disabled=move||busy.get()||!selected_feature.get().is_some_and(|i|molecule.get().features.get(i).is_some_and(|f|f.kind.eq_ignore_ascii_case("CDS"))) on:click=do_cds>"Translate selected CDS"</button>
            <p class="biology-status" role="status">{move||status.get()}</p>
            {move||result.get().map(|saved|{
                let report=StoredValue::new(saved.report.clone());
                view!{<section class="protein-result"><h3>"Computed protein"</h3>
                    <p class="helper">{report.with_value(|r|description(r))}</p>
                    <p class="helper">{report.with_value(|r|format!("Translated location: {}",collagene_core::editing::location_text(&r.effective_location).unwrap_or_else(|e|e.to_string())))}</p>
                    <p class="helper">"The effective location includes the selected orientation; the frame offset is applied along its biological sequence. Stops remain visible for region translations. This is a computed reference, not evidence of expression."</p>
                    <button class="text-tool" on:click=move|_|{
                        if let Some(saved)=result.get_untracked(){match guard.run(saved.clone()).and_then(|()|on_select_location.run((saved.source.revision.clone(),saved.report.effective_location.clone()))){Ok(())=>{},Err(e)=>status.set(e)}}
                    }>"Select full translated location"</button>
                    {report.with_value(|r|r.cds.clone()).map(|cds|view!{<p class="helper">{format!("Complete under this codon model: {} · internal stops: {} · stored translation: {}",cds.complete_under_model,cds.internal_stops,match cds.stored_translation_matches{Some(true)=>"matches",Some(false)=>"differs",None=>"not supplied"})}</p>})}
                    <pre class="protein-preview">{move||report.with_value(|r|{let start=protein_page.get().saturating_mul(600).min(r.protein.len());let end=start.saturating_add(600).min(r.protein.len());r.protein[start..end].to_string()})}</pre>
                    <div class="orf-actions"><button class="text-tool" disabled=move||{protein_page.get()==0} on:click=move|_|{protein_page.update(|p|*p=p.saturating_sub(1));}>"Previous protein page"</button>
                        <span>{move||report.with_value(|r|{let n=r.protein.len();if n==0{"0 amino-acid symbols".into()}else{format!("{}–{} of {}",protein_page.get()*600+1,((protein_page.get()+1)*600).min(n),n)}})}</span>
                        <button class="text-tool" disabled=move||{report.with_value(|r|(protein_page.get()+1)*600>=r.protein.len())} on:click=move|_|{protein_page.update(|p|*p+=1);}>"Next protein page"</button>
                    </div>
                    <ul>{report.with_value(|r|r.warnings.clone()).into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>
                    <label>"Output name"<input aria-label="Translated protein and CDS name" maxlength="1024" prop:value=move||output_name.get() on:input=move|e|{output_name.set(event_target_value(&e));}/></label>
                    <button class="button secondary full" disabled=move||{busy.get()||report.with_value(|r|r.protein.is_empty())} on:click=move|_|{
                        let Some(saved)=result.get_untracked()else{return;};
                        if let Err(e)=guard.run(saved.clone()){status.set(e);return;}
                        let text=fasta_record(&output_name.get_untracked(),&saved.report.protein,&description(&saved.report));
                        crate::app::download_with_completion("translation.fasta","text/plain",text.as_bytes(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e);}})
                    }>"Download full protein FASTA"</button>
                    <button class="button primary full" disabled=move||{busy.get()||report.with_value(|r|r.protein.is_empty())} on:click=move|_|{prepare.run(false);}>"Open linked protein reference"</button>
                    <Show when=move||{report.with_value(|r|matches!(&r.selection,TranslationSelection::Region{..}))}>
                        <p class="helper">"Retaining this region proposes a CDS with explicit frame, genetic code and initiation policy. It requires biological review; translating a frame does not establish that it is coding."</p>
                        <button class="button secondary full" disabled=move||{busy.get()} on:click=move|_|{prepare.run(true);}>"Retain translated region as CDS · review required"</button>
                    </Show>
                </section>}
            })}
            <h3>"Find open reading frames"</h3>
            <label>"Start codons"<select aria-label="ORF start policy" on:change=move|e|start_policy.set(if event_target_value(&e)=="code" {StartPolicy::CodeStarts}else{StartPolicy::AtgOnly})><option value="atg">"ATG/AUG only"</option><option value="code">"Recognized starts for this code"</option></select></label>
            <label>"Incomplete candidates"<select aria-label="ORF partial policy" on:change=move|e|{incomplete.set(match event_target_value(&e).as_str(){"linear"=>IncompletePolicy::IncludeLinearPartials,"all"=>IncompletePolicy::IncludeAllPartials,_=>IncompletePolicy::CompleteOnly});}><option value="complete">"Complete start-to-stop only"</option><option value="linear">"Include linear-edge partials"</option><option value="all">"Also include circular one-turn partials"</option></select></label>
            <label>"Minimum protein length"<input aria-label="ORF minimum amino acids" type="number" min="0" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Maximum candidates"<input aria-label="ORF maximum results" type="number" min="1" max="10000" prop:value=move||maximum.get() on:input=move|e|maximum.set(event_target_value(&e))/></label>
            <label class="biology-check"><input aria-label="ORF both strands" type="checkbox" prop:checked=move||both_strands.get() on:change=move|e|both_strands.set(event_target::<HtmlInputElement>(&e).checked())/>"Search both strands (six offset classes)"</label>
            <button class="button primary full" disabled=move||busy.get()||molecule.get().sequence.alphabet()==Alphabet::Protein on:click=do_orfs>"Find ORFs"</button>
            {move||orfs.get().map(|bound|{
                let report=bound.report.clone();
                let source_revision=report.source_revision.clone();
                view!{<section class="orf-result"><h3>{format!("{} ORF candidates{}",report.orfs.len(),if report.truncated {" · truncated"}else{""})}</h3><ul>{report.warnings.iter().cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>
                <p class="helper">"Every returned candidate is available through pages of 100. Selection and annotation retain the full compound location; predictions require review."</p>
                {report.orfs.iter().cloned().enumerate().skip(orf_page.get()*100).take(100).map(|(i,orf)|{
                    let selected=orf.clone();let annotation=orf.clone();let revision=source_revision.clone();let select_revision=source_revision.clone();let bound_source=bound.source.clone();let annotate_source=bound.source.clone();let bound_epoch=bound.epoch;let code=report.options.genetic_code;let start_policy=report.options.start_policy;
                    view!{<article class="orf-entry"><strong>{format!("ORF {} · {} aa · {}",i+1,orf.protein.len(),if orf.reverse {"reverse"}else{"forward"})}</strong><p>{format!("{} nt · {}{} · {:?}",orf.nucleotide_length,if orf.partial_5prime {"partial 5′ "}else{""},if orf.partial_3prime {"partial 3′"}else{"start/stop bounded"},orf.termination)}</p><details><summary>"Protein and location"</summary><pre>{format!("{}\nLocation (1-based inclusive): {}",preview(&orf.protein),collagene_core::editing::location_text(&orf.location).unwrap_or_else(|_|"unresolved".into()))}</pre></details><div class="orf-actions">
                        <button class="button secondary" on:click=move |_|{if source.get_untracked().binding!=bound_source||epoch.get_untracked()!=bound_epoch{status.set("ORF source event changed; search again.".into());return;}if let Err(e)=on_select_location.run((select_revision.clone(),selected.location.clone())){status.set(e);}}>{format!("Select ORF {}",i+1)}</button>
                        <button class="button secondary" on:click=move |_|{
                            if source.get_untracked().binding!=annotate_source||epoch.get_untracked()!=bound_epoch{status.set("ORF source event changed; search again.".into());return;}
                            let mut qualifiers=vec![("label".into(),Some(format!("Predicted ORF {}",i+1))),("transl_table".into(),Some(code.to_string())),("codon_start".into(),Some("1".into())),("note".into(),Some(format!("Predicted ORF; start policy {start_policy:?}; partial 5′: {}; partial 3′: {}; expression is not established.",annotation.partial_5prime,annotation.partial_3prime)))];
                            if annotation.partial_5prime||annotation.partial_3prime{qualifiers.push(("partial".into(),None));}
                            let feature=Feature{kind:"CDS".into(),location:annotation.location.clone(),qualifiers,needs_review:true};
                            match on_annotate.run((revision.clone(),feature)){Ok(())=>{},Err(e)=>status.set(e)}
                        }>{format!("Annotate ORF {}",i+1)}</button>
                    </div></article>}
                }).collect_view()}
                <div class="orf-actions"><button class="text-tool" disabled=move||{orf_page.get()==0} on:click=move|_|{orf_page.update(|p|*p=p.saturating_sub(1));}>"Previous ORFs"</button><span>{format!("Page {} of {}",orf_page.get()+1,report.orfs.len().div_ceil(100).max(1))}</span><button class="text-tool" disabled={(orf_page.get()+1)*100>=report.orfs.len()} on:click=move|_|{orf_page.update(|p|*p+=1);}>"Next ORFs"</button></div>
                <button class="button secondary full" disabled=move||{orfs.with(|r|r.as_ref().is_none_or(|r|r.report.orfs.is_empty()))} on:click=move |_|{
                    let Some(bound)=orfs.get_untracked()else{return;};if source.get_untracked().binding!=bound.source||epoch.get_untracked()!=bound.epoch{status.set("ORF source event changed; search again.".into());return;}let report=&bound.report;let name=molecule.get_untracked().name;
                    let mut text=String::new();for (i,orf) in report.orfs.iter().enumerate(){text.push_str(&fasta_record(&format!("{name}_ORF_{}",i+1),&orf.protein,&format!("code={} reverse={} partial5={} partial3={}",report.options.genetic_code,orf.reverse,orf.partial_5prime,orf.partial_3prime)));}
                    crate::app::download_with_completion("orf-proteins.fasta","text/plain",text.as_bytes(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e);}})
                }>"Download ORF proteins"</button></section>}
            })}
        </div>
    }.into_any()
}
