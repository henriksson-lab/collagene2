//! Explicit background translation/ORF tools with revision-bound results.
use crate::worker;
use collagene_core::{
    Alphabet, Feature, Molecule,
    protocol::{AnalysisJob, AnalysisResult, Request, Response},
    translation::{
        self, CdsTranslation, IncompletePolicy, InitiationPolicy, OrfOptions, OrfReport,
        StartPolicy, TranslationOptions, TranslationReport,
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
enum ProteinResult {
    Frame(TranslationReport),
    Cds(CdsTranslation),
}
impl ProteinResult {
    fn protein(&self) -> &str {
        match self {
            Self::Frame(report) => &report.protein,
            Self::Cds(report) => &report.protein,
        }
    }
    fn warnings(&self) -> &[String] {
        match self {
            Self::Frame(report) => &report.warnings,
            Self::Cds(report) => &report.warnings,
        }
    }
    fn description(&self) -> String {
        match self {
            Self::Frame(r) => format!(
                "Code {} · frame {} · {} · {} amino-acid symbols",
                r.options.genetic_code,
                u16::from(r.options.frame) + 1,
                if r.options.reverse {
                    "reverse"
                } else {
                    "forward"
                },
                r.protein.len()
            ),
            Self::Cds(r) => format!(
                "{} · code {} · codon_start {} · {} amino-acid symbols",
                r.feature_label.chars().take(160).collect::<String>(),
                r.translation.options.genetic_code,
                r.codon_start,
                r.protein.len()
            ),
        }
    }
}
fn preview(protein: &str) -> String {
    let shown: String = protein.chars().take(600).collect();
    if protein.chars().count() > 600 {
        format!("{shown}…\nPreview limited to 600 symbols; download keeps the full protein.")
    } else {
        shown
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
) -> impl IntoView {
    let code = RwSignal::new(1u8);
    let frame = RwSignal::new(0u8);
    let reverse = RwSignal::new(false);
    let initiation = RwSignal::new(InitiationPolicy::Ordinary);
    let start_policy = RwSignal::new(StartPolicy::AtgOnly);
    let incomplete = RwSignal::new(IncompletePolicy::CompleteOnly);
    let minimum = RwSignal::new("20".to_string());
    let maximum = RwSignal::new("1000".to_string());
    let both_strands = RwSignal::new(true);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let result = RwSignal::new(None::<ProteinResult>);
    let orfs = RwSignal::new(None::<OrfReport>);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    Effect::new(move |_| {
        let _ = (
            molecule.get().revision(),
            code.get(),
            frame.get(),
            reverse.get(),
            initiation.get(),
            start_policy.get(),
            incomplete.get(),
            minimum.get(),
            maximum.get(),
            both_strands.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        result.set(None);
        orfs.set(None);
        status.set(String::new());
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = selected_feature.get();
        if result.with_untracked(|r| matches!(r, Some(ProteinResult::Cds(_)))) {
            result.set(None);
        }
    });
    let run = {
        let alive = alive.clone();
        move |job: AnalysisJob, revision: String| {
            let requested_cds = if let AnalysisJob::TranslateCds { feature_index, .. } = &job {
                Some(*feature_index)
            } else {
                None
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Computing in background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
                if !alive.load(Ordering::Relaxed) {
                    return;
                }
                if generation.get_untracked() != token
                    || molecule.get_untracked().revision() != revision
                {
                    return;
                }
                busy.set(false);
                if requested_cds.is_some_and(|i| selected_feature.get_untracked() != Some(i)) {
                    status.set("Selected CDS changed; request its translation again.".into());
                    return;
                }
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result: computed,
                    }) if id == token && returned == revision => {
                        match computed {
                            AnalysisResult::Translation(report) => {
                                result.set(Some(ProteinResult::Frame(report)));
                                status.set("Frame translation complete; stops remain visible as * and unresolved codons as X.".into());
                            }
                            AnalysisResult::CdsTranslation(report) => {
                                result.set(Some(ProteinResult::Cds(report)));
                                status.set("CDS translated from its location and qualifiers; stored protein text was not used as the answer.".into());
                            }
                            AnalysisResult::Orfs(report) => {
                                status.set(format!("{} ORF candidate(s){}; candidates do not establish expression.",report.orfs.len(),if report.truncated {" · result limit reached"}else{""}));
                                orfs.set(Some(report));
                            }
                            AnalysisResult::Error(error) => status.set(error),
                            _ => status.set("Unexpected biology worker response".into()),
                        }
                    }
                    Ok(_) => {
                        status.set("Worker response does not match the requested revision".into())
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let do_translate = {
        let run = run.clone();
        move |_| {
            let m = molecule.get_untracked();
            let revision = m.revision();
            result.set(None);
            run(
                AnalysisJob::Translate {
                    molecule: m,
                    options: TranslationOptions {
                        genetic_code: code.get_untracked(),
                        frame: frame.get_untracked(),
                        reverse: reverse.get_untracked(),
                        initiation: initiation.get_untracked(),
                    },
                },
                revision,
            );
        }
    };
    let do_cds = {
        let run = run.clone();
        move |_| {
            let m = molecule.get_untracked();
            let Some(index) = selected_feature.get_untracked().filter(|i| {
                m.features
                    .get(*i)
                    .is_some_and(|f| f.kind.eq_ignore_ascii_case("CDS"))
            }) else {
                status.set("Select a CDS annotation first.".into());
                return;
            };
            let revision = m.revision();
            result.set(None);
            run(
                AnalysisJob::TranslateCds {
                    molecule: m,
                    feature_index: index,
                },
                revision,
            );
        }
    };
    let do_orfs = move |_| {
        let (Ok(min), Ok(max)) = (
            minimum.get_untracked().parse::<usize>(),
            maximum.get_untracked().parse::<usize>(),
        ) else {
            status.set("ORF lengths and result limits must be non-negative integers.".into());
            return;
        };
        let m = molecule.get_untracked();
        let revision = m.revision();
        orfs.set(None);
        run(
            AnalysisJob::FindOrfs {
                molecule: m,
                options: OrfOptions {
                    genetic_code: code.get_untracked(),
                    start_policy: start_policy.get_untracked(),
                    incomplete: incomplete.get_untracked(),
                    min_amino_acids: min,
                    max_orfs: max,
                    both_strands: both_strands.get_untracked(),
                },
            },
            revision,
        );
    };
    view! {
        <style>".biology-panel pre{white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.6 monospace;max-height:240px;overflow:auto;background:#f4f7f4;padding:10px;border-radius:5px}.biology-panel .biology-check{display:flex;align-items:flex-start;gap:8px}.biology-panel .biology-check input{width:auto;flex-shrink:0;margin-top:3px}.biology-panel .orf-entry{padding:10px 0;border-bottom:1px solid #e1e7e0;font-size:13px;line-height:1.5}.biology-panel .orf-actions{display:flex;gap:6px;flex-wrap:wrap}.biology-panel .orf-actions button{font-size:12px;padding:6px 8px}.biology-panel details{margin:10px 0}.biology-panel ul{padding-left:20px;overflow-wrap:anywhere}.biology-status{font-size:13px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="inspector-content biology-panel">
            <h2>"Translation & ORFs"</h2>
            <Show when=move||molecule.get().sequence.alphabet()==Alphabet::Protein><p class="helper">"Choose a DNA or RNA document for nucleotide translation and ORF search."</p></Show>
            <label>"Genetic code"<select aria-label="Genetic code" prop:value=move||code.get().to_string() on:change=move|e|code.set(event_target_value(&e).parse().unwrap_or(1))>{translation::genetic_codes().into_iter().map(|c|view!{<option value=c.id.to_string()>{format!("{} · {}",c.id,c.name)}</option>}).collect_view()}</select></label>
            <h3>"Read a frame"</h3>
            <label>"Frame offset"<select aria-label="Translation frame" prop:value=move||frame.get().to_string() on:change=move|e|frame.set(event_target_value(&e).parse().unwrap_or(0))><option value="0">"1 · skip 0 bases"</option><option value="1">"2 · skip 1 base"</option><option value="2">"3 · skip 2 bases"</option></select></label>
            <label class="biology-check"><input type="checkbox" aria-label="Reverse translation" prop:checked=move||reverse.get() on:change=move|e|reverse.set(event_target::<HtmlInputElement>(&e).checked())/>"Reverse complement"</label>
            <label>"First codon"<select aria-label="Initiation policy" on:change=move|e|initiation.set(match event_target_value(&e).as_str(){"recognized"=>InitiationPolicy::RecognizedStartAsMethionine,"required"=>InitiationPolicy::RequireRecognizedStart,_=>InitiationPolicy::Ordinary})><option value="ordinary">"Ordinary internal-codon translation"</option><option value="recognized">"Recognized initiator → M"</option><option value="required">"Require a recognized initiator"</option></select></label>
            <button class="button primary full" disabled=move||busy.get()||molecule.get().sequence.alphabet()==Alphabet::Protein on:click=do_translate>"Translate frame"</button>
            <h3>"Selected coding feature"</h3>
            <p class="helper">{move||selected_feature.get().and_then(|i|molecule.get().features.get(i).cloned()).filter(|f|f.kind.eq_ignore_ascii_case("CDS")).map(|f|format!("Selected CDS: {}",f.label().chars().take(160).collect::<String>())).unwrap_or_else(||"Select a CDS in the feature list. CDS phase and genetic code come from its qualifiers.".into())}</p>
            <button class="button secondary full" disabled=move||busy.get()||!selected_feature.get().is_some_and(|i|molecule.get().features.get(i).is_some_and(|f|f.kind.eq_ignore_ascii_case("CDS"))) on:click=do_cds>"Translate selected CDS"</button>
            <p class="biology-status" role="status">{move||status.get()}</p>
            {move||result.get().map(|report|view!{
                <section class="protein-result"><h3>"Computed protein"</h3><p class="helper">{report.description()}</p><pre class="protein-preview">{preview(report.protein())}</pre>
                {match &report {ProteinResult::Cds(cds)=>Some(view!{<p class="helper">{format!("Complete under this codon model: {} · internal stops: {} · stored translation: {}",cds.complete_under_model,cds.internal_stops,match cds.stored_translation_matches{Some(true)=>"matches",Some(false)=>"differs",None=>"not supplied"})}</p>}),_=>None}}
                <ul>{report.warnings().iter().cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>
                <button class="button secondary full" disabled=report.protein().is_empty() on:click=move |_|{
                    let Some(report)=result.get_untracked() else{return;};
                    let m=molecule.get_untracked();let text=fasta_record(&m.name,report.protein(),&report.description());
                    if let Err(e)=crate::app::download("translation.fasta","text/plain",text.as_bytes()){status.set(e);}
                }>"Download protein FASTA"</button></section>
            })}
            <h3>"Find open reading frames"</h3>
            <label>"Start codons"<select aria-label="ORF start policy" on:change=move|e|start_policy.set(if event_target_value(&e)=="code" {StartPolicy::CodeStarts}else{StartPolicy::AtgOnly})><option value="atg">"ATG/AUG only"</option><option value="code">"Recognized starts for this code"</option></select></label>
            <label>"Incomplete candidates"<select aria-label="ORF partial policy" on:change=move|e|incomplete.set(match event_target_value(&e).as_str(){"linear"=>IncompletePolicy::IncludeLinearPartials,"all"=>IncompletePolicy::IncludeAllPartials,_=>IncompletePolicy::CompleteOnly})><option value="complete">"Complete start-to-stop only"</option><option value="linear">"Include linear-edge partials"</option><option value="all">"Also include circular one-turn partials"</option></select></label>
            <label>"Minimum protein length"<input aria-label="ORF minimum amino acids" type="number" min="0" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Maximum candidates"<input aria-label="ORF maximum results" type="number" min="1" max="10000" prop:value=move||maximum.get() on:input=move|e|maximum.set(event_target_value(&e))/></label>
            <label class="biology-check"><input aria-label="ORF both strands" type="checkbox" prop:checked=move||both_strands.get() on:change=move|e|both_strands.set(event_target::<HtmlInputElement>(&e).checked())/>"Search both strands (six offset classes)"</label>
            <button class="button primary full" disabled=move||busy.get()||molecule.get().sequence.alphabet()==Alphabet::Protein on:click=do_orfs>"Find ORFs"</button>
            {move||orfs.get().map(|report|{
                let source_revision=report.source_revision.clone();
                view!{<section class="orf-result"><h3>{format!("{} ORF candidates{}",report.orfs.len(),if report.truncated {" · truncated"}else{""})}</h3><ul>{report.warnings.iter().cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>
                <p class="helper">"Displaying at most 100 candidates. Annotation stores the full compound location; predictions require review."</p>
                {report.orfs.into_iter().take(100).enumerate().map(|(i,orf)|{
                    let selected=orf.clone();let annotation=orf.clone();let revision=source_revision.clone();let code=report.options.genetic_code;let start_policy=report.options.start_policy;
                    view!{<article class="orf-entry"><strong>{format!("ORF {} · {} aa · {}",i+1,orf.protein.len(),if orf.reverse {"reverse"}else{"forward"})}</strong><p>{format!("{} nt · {}{} · {:?}",orf.nucleotide_length,if orf.partial_5prime {"partial 5′ "}else{""},if orf.partial_3prime {"partial 3′"}else{"start/stop bounded"},orf.termination)}</p><details><summary>"Protein and location"</summary><pre>{format!("{}\nLocation (1-based inclusive): {}",preview(&orf.protein),collagene_core::editing::location_text(&orf.location).unwrap_or_else(|_|"unresolved".into()))}</pre></details><div class="orf-actions">
                        <button class="button secondary" on:click=move |_|{let spans=selected.location.spans();if let Some((a,b,_))=spans.first(){selected_feature.set(None);selection.set(Some((*a,*b)));if spans.len()>1{status.set(format!("Compound ORF has {} spans; selection displays its first reference span. Annotation retains the full location.",spans.len()));}}}>{format!("Select ORF {}",i+1)}</button>
                        <button class="button secondary" on:click=move |_|{
                            let mut qualifiers=vec![("label".into(),Some(format!("Predicted ORF {}",i+1))),("transl_table".into(),Some(code.to_string())),("codon_start".into(),Some("1".into())),("note".into(),Some(format!("Predicted ORF; start policy {start_policy:?}; partial 5′: {}; partial 3′: {}; expression is not established.",annotation.partial_5prime,annotation.partial_3prime)))];
                            if annotation.partial_5prime||annotation.partial_3prime{qualifiers.push(("partial".into(),None));}
                            let feature=Feature{kind:"CDS".into(),location:annotation.location.clone(),qualifiers,needs_review:true};
                            match on_annotate.run((revision.clone(),feature)){Ok(())=>{},Err(e)=>status.set(e)}
                        }>{format!("Annotate ORF {}",i+1)}</button>
                    </div></article>}
                }).collect_view()}
                <button class="button secondary full" disabled=move||orfs.with(|r|r.as_ref().is_none_or(|r|r.orfs.is_empty())) on:click=move |_|{
                    let Some(report)=orfs.get_untracked()else{return;};let name=molecule.get_untracked().name;
                    let mut text=String::new();for (i,orf) in report.orfs.iter().enumerate(){text.push_str(&fasta_record(&format!("{name}_ORF_{}",i+1),&orf.protein,&format!("code={} reverse={} partial5={} partial3={}",report.options.genetic_code,orf.reverse,orf.partial_5prime,orf.partial_3prime)));}
                    if let Err(e)=crate::app::download("orf-proteins.fasta","text/plain",text.as_bytes()){status.set(e);}
                }>"Download ORF proteins"</button></section>}
            })}
        </div>
    }
}
