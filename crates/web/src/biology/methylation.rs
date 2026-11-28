//! Explicit-assumption restriction planning review, bound to one source event.
use crate::io::download::prepare_download;
use crate::platform::worker::CachedWorker;
use crate::reads::contigs::ReadBinding;
use crate::reads::contigs::ReadNavigation;
use collagene_core::{
    chemistry::methylation::{
        self, MethylationOptions, MethylationReport, MethylationState, Outcome, Sensitivity,
    },
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule};
use leptos::prelude::*;
use serde::Serialize;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, PartialEq, Eq)]
pub struct SourceInfo {
    pub binding: ReadBinding,
    pub name: String,
    pub length: usize,
    pub alphabet: Alphabet,
}
#[derive(Clone)]
struct BoundReview {
    report: Arc<MethylationReport>,
    source: SourceInfo,
    epoch: u64,
}
fn state_key(state: MethylationState) -> &'static str {
    match state {
        MethylationState::Unknown => "unknown",
        MethylationState::Unmodified => "unmodified",
        MethylationState::FullyModified => "fully",
        MethylationState::HemiOrMixed => "mixed",
    }
}
fn parse_state(value: &str) -> MethylationState {
    match value {
        "unmodified" => MethylationState::Unmodified,
        "fully" => MethylationState::FullyModified,
        "mixed" => MethylationState::HemiOrMixed,
        _ => MethylationState::Unknown,
    }
}
fn state_label(state: MethylationState) -> &'static str {
    match state {
        MethylationState::Unknown => "Unknown",
        MethylationState::Unmodified => "Assume unmodified",
        MethylationState::FullyModified => "Assume fully modified",
        MethylationState::HemiOrMixed => "Hemi / mixed (unsupported model)",
    }
}
fn outcome_label(outcome: Outcome) -> &'static str {
    match outcome {
        Outcome::PermittedUnderProfile => "Permitted under profile",
        Outcome::Blocked => "Blocked under assumptions",
        Outcome::Impaired => "Impaired under assumptions",
        Outcome::Unknown => "Unknown / unresolved",
    }
}
fn sensitivity_label(sensitivity: Sensitivity) -> &'static str {
    match sensitivity {
        Sensitivity::NotSensitive => "No reported sensitivity",
        Sensitivity::Blocked => "Blocked",
        Sensitivity::Impaired => "Impaired",
        Sensitivity::BlockedByOverlapping => "Blocked by an overlapping motif",
        Sensitivity::BlockedBySomeOverlaps => "Blocked in some overlap combinations",
        Sensitivity::ImpairedBySomeOverlaps => "Impaired in some overlap combinations",
        Sensitivity::Unsupported => "Unsupported / conflicting profile",
    }
}
fn boundary_label(boundary: Option<usize>) -> String {
    boundary
        .map(|n| format!("{n} bases from source start"))
        .unwrap_or_else(|| "unavailable in this source".into())
}
fn alphabet_label(alphabet: Alphabet) -> &'static str {
    match alphabet {
        Alphabet::Dna => "DNA",
        Alphabet::Rna => "RNA",
        Alphabet::Protein => "protein",
    }
}
fn whole(value: &str, label: &str) -> Result<usize, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {label}."))
}
struct Json(Vec<u8>);
impl std::io::Write for Json {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > 32 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "Methylation JSON exceeds the 32 MiB browser output bound",
            ));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
#[derive(Serialize)]
struct ExportReview<'a> {
    kind: &'static str,
    source_document_id: &'a str,
    source_event_id: &'a str,
    source_revision: &'a str,
    report: &'a MethylationReport,
}

#[component]
fn StateField(category: &'static str, value: RwSignal<MethylationState>) -> impl IntoView {
    view! {<label>{format!("{category} source assumption")}<select aria-label=format!("Methylation {category} assumption") prop:value=move||state_key(value.get()) on:change=move|e|value.set(parse_state(&event_target_value(&e)))>{[MethylationState::Unknown,MethylationState::Unmodified,MethylationState::FullyModified,MethylationState::HemiOrMixed].into_iter().map(|state|view!{<option value=state_key(state)>{state_label(state)}</option>}).collect_view()}</select></label>}
}

#[component]
pub fn MethylationPanel(
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<ReadNavigation>,
    on_digest: Callback<()>,
) -> impl IntoView {
    let catalog = StoredValue::new(methylation::profiles());
    let selected = RwSignal::new(
        catalog.with_value(|v| v.iter().map(|p| p.enzyme.clone()).collect::<Vec<_>>()),
    );
    let dam = RwSignal::new(MethylationState::Unknown);
    let dcm = RwSignal::new(MethylationState::Unknown);
    let cpg = RwSignal::new(MethylationState::Unknown);
    let cells = RwSignal::new("20000000".to_string());
    let site_limit = RwSignal::new("4096".to_string());
    let options = Memo::new(move |_| -> Result<MethylationOptions, String> {
        Ok(MethylationOptions {
            enzymes: selected.get(),
            dam: dam.get(),
            dcm: dcm.get(),
            cpg: cpg.get(),
            max_comparisons: whole(&cells.get(), "comparison budget")?,
            max_sites: whole(&site_limit.get(), "site limit")?,
        })
    });
    let data = RwSignal::new(None::<Arc<BoundReview>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
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
    let page = RwSignal::new(0usize);
    let focus = RwSignal::new(None::<usize>);
    let enzyme_filter = RwSignal::new(String::new());
    let outcome_filter = RwSignal::new("all".to_string());
    let expanded = RwSignal::new(false);
    Effect::new(move |_| {
        let _ = (source.get(), options.get(), epoch.get());
        generation.update(|v| *v = v.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        data.set(None);
        busy.set(false);
        page.set(0);
        focus.set(None);
        status.set(String::new());
    });
    let guard = Callback::new(move |bound: Arc<BoundReview>| -> Result<(), String> {
        if source.get_untracked() != bound.source
            || epoch.get_untracked() != bound.epoch
            || options.get_untracked() != Ok(bound.report.options.clone())
        {
            return Err("The source document, event or assumptions changed. Assess the current source again.".into());
        }
        Ok(())
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            let opts = match options.get_untracked() {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let info = source.get_untracked();
            if info.alphabet != Alphabet::Dna {
                status.set("Methylation assessment requires a DNA source.".into());
                return;
            }
            let molecule = get_molecule.run(());
            if let Err(error) = methylation::source_preflight(&molecule, &opts) {
                status.set(error.to_string());
                return;
            }
            if molecule.revision() != info.binding.revision {
                status.set("Source identity changed. Reopen this review.".into());
                return;
            }
            let revision = match methylation::batch_revision(&molecule, &opts) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            if let Some(w) = worker.get_value() {
                w.close();
            }
            let w = match CachedWorker::new() {
                Ok(v) => Rc::new(v),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            generation.update(|v| *v = v.wrapping_add(1));
            let token = generation.get_untracked();
            let expected_epoch = epoch.get_untracked();
            busy.set(true);
            data.set(None);
            status.set(
                "Inspecting recognition and methylation contexts in a background worker…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::AssessMethylation {
                            molecule: molecule.clone(),
                            options: opts.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || epoch.get_untracked() != expected_epoch
                    || source.get_untracked() != info
                    || options.get_untracked() != Ok(opts.clone())
                {
                    return;
                }
                busy.set(false);
                let result = (|| -> Result<MethylationReport, String> {
                    let response = result?;
                    if response.id != token || response.revision != revision {
                        return Err("Methylation worker returned a stale source binding.".into());
                    }
                    let AnalysisResult::Methylation(report) = response.result else {
                        return Err(match response.result {
                            AnalysisResult::Error(e) => e,
                            _ => "Unexpected methylation worker result.".into(),
                        });
                    };
                    report
                        .validate_bindings_for(&molecule, &opts)
                        .map_err(|e| e.to_string())?;
                    Ok(report)
                })();
                match result {
                    Ok(report) => {
                        status.set(if report.complete{format!("Context assessment completed: {} recognition sites. Outcomes remain conditional on the declared source assumptions.",report.sites.len())}else{format!("Incomplete context assessment: {} retained sites. Review resource limits before interpreting absent sites.",report.sites.len())});
                        data.set(Some(Arc::new(BoundReview {
                            report: Arc::new(report),
                            source: info,
                            epoch: expected_epoch,
                        })));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let cancel = move |_| {
        generation.update(|v| *v = v.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        status.set("Methylation review cancelled.".into());
    };
    let export = move |_| {
        let Some(bound) = data.get_untracked() else {
            return;
        };
        if let Err(e) = guard.run(bound.clone()) {
            status.set(e);
            return;
        }
        let mut bytes = Json(vec![]);
        let export = ExportReview {
            kind: "methylation-assessment-review-v1",
            source_document_id: &bound.source.binding.document_id,
            source_event_id: &bound.source.binding.event_id,
            source_revision: &bound.source.binding.revision,
            report: &bound.report,
        };
        let result = serde_json::to_writer_pretty(&mut bytes, &export)
            .map_err(|e| e.to_string())
            .and_then(|_| {
                prepare_download("methylation-assessment.json", "application/json", &bytes.0)
            });
        crate::io::download::complete_download(result, move |result| {
            if status.try_get_untracked().is_some() {
                if let Err(e) = result {
                    status.set(e);
                }
            }
        });
    };
    let filtered = Memo::new(move |_| {
        data.with(|bound| {
            bound
                .as_ref()
                .map(|b| {
                    b.report
                        .sites
                        .iter()
                        .enumerate()
                        .filter(|(_, s)| {
                            (enzyme_filter.get().is_empty() || s.enzyme == enzyme_filter.get())
                                && (outcome_filter.get() == "all"
                                    || format!("{:?}", s.outcome) == outcome_filter.get())
                        })
                        .map(|(i, _)| i)
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        })
    });
    Effect::new(move |_| {
        let _ = (enzyme_filter.get(), outcome_filter.get());
        page.set(0);
    });
    let navigate = Callback::new(move |at: usize| {
        let Some(bound) = data.get_untracked() else {
            return;
        };
        if let Err(e) = guard.run(bound.clone()) {
            status.set(e);
            return;
        }
        if at < bound.source.length {
            on_navigate.run(ReadNavigation {
                binding: bound.source.binding.clone(),
                call_index: at,
                epoch: bound.epoch,
            });
        }
    });
    view! {<section class="methylation-panel" class:methylation-expanded=move||expanded.get()>
        <div class="methylation-heading"><h2>"Methylation review"</h2><button class="button secondary" on:click=move|_|expanded.update(|v|*v=!*v)>{move||if expanded.get(){"Close expanded methylation"}else{"Expand methylation"}}</button></div>
        <p class="helper">"Review possible restriction-enzyme sensitivity under explicit Dam, Dcm and CpG source assumptions. Sequence alone does not identify modified bases. These assumptions do not change the source molecule or verify a digest."</p>
        <p class="methylation-source">{move||source.with(|s|format!("{} · {} bases · {}",s.name,s.length,alphabet_label(s.alphabet)))}</p>
        <StateField category="Dam" value=dam/><StateField category="Dcm" value=dcm/><StateField category="CpG" value=cpg/>
        <p class="helper">"Fully modified means the category's modeled target bases on both strands. Hemi-methylated or mixed populations require an unsupported model and remain unresolved where relevant. Unknown is the default; organism names and imported annotations do not establish these states."</p>
        <details class="methylation-enzyme-picker"><summary>"Enzymes and versioned profiles"</summary><button class="button secondary" on:click=move|_|selected.set(catalog.with_value(|v|v.iter().map(|p|p.enzyme.clone()).collect()))>"Select all methylation profiles"</button><button class="button secondary" on:click=move|_|selected.set(vec![])>"Clear methylation profiles"</button>{catalog.with_value(|v|v.iter().map(|p|{let name=p.enzyme.clone();let checked=name.clone();let changed=name.clone();view!{<label class="methylation-enzyme-choice"><input type="checkbox" aria-label=format!("Methylation enzyme {name}") prop:checked=move||selected.with(|v|v.contains(&checked)) on:change=move|e|{let yes=event_target_checked(&e);selected.update(|v|{if yes{if !v.contains(&changed){v.push(changed.clone());}}else{v.retain(|n|n!=&changed);}});}/><span>{format!("{} · {}",p.enzyme,p.recognition)}</span></label>}}).collect_view())}<p class="helper">"Small supplier-specific starter profiles. Enzyme variants and substitute suppliers are not interchangeable."</p></details>
        <details><summary>"Resource limits"</summary><label>"Comparison budget"<input aria-label="Methylation comparison budget" prop:value=move||cells.get() on:input=move|e|cells.set(event_target_value(&e))/></label><label>"Retained site limit"<input aria-label="Methylation site limit" prop:value=move||site_limit.get() on:input=move|e|site_limit.set(event_target_value(&e))/></label></details>
        <div class="methylation-actions"><button class="button primary" disabled=move||busy.get()||source.with(|s|s.alphabet!=Alphabet::Dna) on:click=run>"Assess methylation"</button><button class="button secondary" disabled=move||!busy.get() on:click=cancel>"Cancel methylation assessment"</button></div>
        <p class="methylation-status" role="status">{move||status.get()}</p>
        {move||data.get().map(|bound|{let complete=bound.report.complete;view!{<div class="methylation-result" data-complete=complete.to_string()>
            <h3>"Restriction planning review"</h3><p class="helper">"Permitted under a profile means no modeled methylation restriction was established under these assumptions. It does not predict complete cleavage, reaction yield, reagent conditions, star activity or experimental success. Existing cleavage-based plans and gels retain their own assumptions."</p>
            <p>{format!("{} sites · {} comparisons · {} context scan",bound.report.sites.len(),bound.report.comparisons,if complete{"complete"}else{"incomplete"})}</p>
            <div class="methylation-summary">{[Outcome::PermittedUnderProfile,Outcome::Blocked,Outcome::Impaired,Outcome::Unknown].into_iter().map(|outcome|{let count=bound.report.sites.iter().filter(|s|s.outcome==outcome).count();view!{<p>{format!("{}: {}",outcome_label(outcome),count)}</p>}}).collect_view()}</div>
            <button class="button secondary" on:click=export>"Download methylation report"</button><button class="button secondary" on:click=move|_|on_digest.run(())>"Review cleavage assumptions"</button>
            <details><summary>"Source binding and assumptions"</summary><p>{format!("Document {} · event {} · source revision {}",bound.source.binding.document_id,bound.source.binding.event_id,bound.source.binding.revision)}</p><p>{format!("Dam: {} · Dcm: {} · CpG: {}",state_label(bound.report.options.dam),state_label(bound.report.options.dcm),state_label(bound.report.options.cpg))}</p><p>{format!("Engine {} · profiles {} · profile revision {} · options {}",bound.report.engine_version,bound.report.profile_version,bound.report.profile_revision,bound.report.options_revision)}</p></details>
            <details><summary>"Profile identities and references"</summary>{bound.report.profiles.iter().map(|p|view!{<article class="methylation-profile"><strong>{p.enzyme.clone()}</strong><p>{format!("{} · {} · Dam / Dcm / CpG sensitivity: {}",p.supplier_identity,p.recognition,p.sensitivities.iter().map(|s|sensitivity_label(*s)).collect::<Vec<_>>().join(" / "))}</p>{p.notes.iter().cloned().map(|n|view!{<p>{n}</p>}).collect_view()}{p.sources.iter().cloned().map(|url|view!{<a href=url.clone() target="_blank" rel="noopener noreferrer">{url.clone()}</a>}).collect_view()}</article>}).collect_view()}</details>
            <h3>"Recognition sites and context evidence"</h3>
            <label>"Site enzyme filter"<select aria-label="Methylation site enzyme" prop:value=move||enzyme_filter.get() on:change=move|e|enzyme_filter.set(event_target_value(&e))><option value="">"All assessed enzymes"</option>{bound.report.profiles.iter().map(|p|view!{<option value=p.enzyme.clone()>{p.enzyme.clone()}</option>}).collect_view()}</select></label>
            <label>"Outcome filter"<select aria-label="Methylation site outcome" prop:value=move||outcome_filter.get() on:change=move|e|outcome_filter.set(event_target_value(&e))><option value="all">"All outcomes"</option>{[Outcome::PermittedUnderProfile,Outcome::Blocked,Outcome::Impaired,Outcome::Unknown].into_iter().map(|o|view!{<option value=format!("{o:?}")>{outcome_label(o)}</option>}).collect_view()}</select></label>
            <p>{move||format!("{} matching sites · page {}",filtered.with(Vec::len),page.get()+1)}</p>
            {move||data.get().map(|b|filtered.get().into_iter().skip(page.get()*25).take(25).filter_map(|index|b.report.sites.get(index).map(|s|{let index=index;view!{<button class="button secondary methylation-site-entry" data-site-index=s.index on:click=move|_|focus.set(Some(index))>{format!("{} · base {}{} · {} · {}",s.enzyme,s.recognition_start+1,if s.reverse{" reverse"}else{" forward"},if s.certain{"certain recognition"}else{"possible recognition"},outcome_label(s.outcome))}</button>}})).collect_view())}
            <div class="methylation-actions"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|v|*v=v.saturating_sub(1))>"Previous methylation sites"</button><button class="button secondary" disabled=move||{(page.get()+1)*25>=filtered.with(Vec::len)} on:click=move|_|page.update(|v|*v+=1)>"Next methylation sites"</button></div>
            {move||data.get().and_then(|b|focus.get().and_then(|i|b.report.sites.get(i).cloned())).map(|s|{let start=s.recognition_start;view!{<article class="methylation-site-detail"><h4>{format!("{} · {}",s.enzyme,outcome_label(s.outcome))}</h4><p>{format!("Recognition starts at source base {} (1-based) · length {} · {} · {}. Top cut boundary: {}; bottom cut boundary: {}. Cut coordinates count bases before each boundary, not base numbers.",s.recognition_start+1,s.recognition_length,if s.reverse{"reverse strand"}else{"forward strand"},if s.certain{"certain recognition"}else{"possible recognition"},boundary_label(s.top_cut),boundary_label(s.bottom_cut))}</p><button class="button secondary" on:click=move|_|navigate.run(start)>"Open recognition start"</button>{s.effects.into_iter().map(|effect|view!{<section class="methylation-effect"><h4>{format!("{:?}: {}",effect.category,outcome_label(effect.outcome))}</h4><p>{format!("{} · profile sensitivity: {}",state_label(effect.state),sensitivity_label(effect.sensitivity))}</p><p>{effect.reason}</p>{effect.overlaps.into_iter().map(|o|view!{<div class="methylation-context"><p>{format!("Motif {} · starts at source base {} · {} · original source bases {} (1-based)",o.motif,o.start+1,if o.certain{"certain motif"}else{"possible motif"},o.positions.iter().map(|p|(p+1).to_string()).collect::<Vec<_>>().join(", "))}</p>{o.modified_bases.into_iter().map(|base|{let position=base.position;view!{<p>{format!("Modeled target {} · original source base {} · {} strand",base.base,base.position+1,if base.reverse{"reverse"}else{"forward"})}<button class="text-tool" on:click=move|_|navigate.run(position)>"Open target base"</button></p>}}).collect_view()}</div>}).collect_view()}</section>}).collect_view()}</article>}})}
            {bound.report.limits.iter().map(|v|view!{<p class="methylation-limit" role="alert">{v.clone()}</p>}).collect_view()}{bound.report.warnings.iter().map(|v|view!{<p class="helper">{v.replace("PermittedUnderProfile","Permitted under profile")}</p>}).collect_view()}
        </div>}})}
    </section>}
}
