//! Explicit preparation graphs for TA/GC ligation and oriented overlap profiles.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    Molecule,
    analysis::Primer,
    chemistry::{EndRepairMethod, EndTreatment},
    overlap_profiles::{OverlapMethod, ProfileOptions, ProfileReport},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    tailing::TailingBase,
    transfer::Candidate,
    verification::{
        DesignPlan, Operation, OrientedSource, PlanInput, Source, TailPreview, TailSubstrate,
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

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Tailing,
    Overlap,
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Treatment {
    Keep,
    Phosphorylate,
    Dephosphorylate,
}
#[derive(Clone, PartialEq, Eq)]
pub struct TailInput {
    pub document: usize,
    pub pcr: bool,
    pub enzymes: String,
    pub fragment: String,
    pub blunt: bool,
    pub forward: String,
    pub reverse: String,
    pub forward_tail: String,
    pub reverse_tail: String,
    pub base: TailingBase,
    pub treatment: Treatment,
    pub reversed: bool,
}
impl TailInput {
    fn new(pcr: bool, base: TailingBase) -> Self {
        Self {
            document: 0,
            pcr,
            enzymes: "EcoRV".into(),
            fragment: "1".into(),
            blunt: false,
            forward: String::new(),
            reverse: String::new(),
            forward_tail: String::new(),
            reverse_tail: String::new(),
            base,
            treatment: Treatment::Keep,
            reversed: false,
        }
    }
    fn recipe(&self) -> Result<TailSubstrate, String> {
        if self.pcr {
            Ok(TailSubstrate::ExactPcr {
                forward: Primer::new("Forward", &self.forward, &self.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new("Reverse", &self.reverse, &self.reverse_tail)
                    .map_err(|e| e.to_string())?,
            })
        } else {
            let enzymes = self
                .enzymes
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect::<Vec<_>>();
            if enzymes.is_empty() {
                return Err("Choose explicit restriction enzymes".into());
            }
            let fragment_index = self
                .fragment
                .parse::<usize>()
                .ok()
                .and_then(|n| n.checked_sub(1))
                .ok_or("Fragment numbers start at 1")?;
            Ok(TailSubstrate::Digest {
                enzymes,
                fragment_index,
                blunt_before: self.blunt,
            })
        }
    }
}
#[derive(Clone, PartialEq, Eq)]
pub struct OverlapInput {
    pub document: usize,
    pub reverse: bool,
}
#[derive(Clone, PartialEq, Eq)]
pub struct MethodsState {
    pub mode: Mode,
    pub vector: TailInput,
    pub insert: TailInput,
    pub circular: bool,
    pub allow_nicked: bool,
    pub rows: Vec<OverlapInput>,
    pub overlaps: String,
    pub method: OverlapMethod,
    pub scan_repeats: bool,
    pub comparisons: String,
    pub occurrences: String,
}
impl Default for MethodsState {
    fn default() -> Self {
        Self {
            mode: Mode::Tailing,
            vector: TailInput::new(false, TailingBase::T),
            insert: TailInput::new(true, TailingBase::A),
            circular: true,
            allow_nicked: false,
            rows: vec![
                OverlapInput {
                    document: 0,
                    reverse: false,
                },
                OverlapInput {
                    document: 0,
                    reverse: false,
                },
            ],
            overlaps: "20".into(),
            method: OverlapMethod::HiFi,
            scan_repeats: true,
            comparisons: "5000000".into(),
            occurrences: "1000".into(),
        }
    }
}
impl MethodsState {
    fn input(&self, vector: bool) -> &TailInput {
        if vector { &self.vector } else { &self.insert }
    }
    fn input_mut(&mut self, vector: bool) -> &mut TailInput {
        if vector {
            &mut self.vector
        } else {
            &mut self.insert
        }
    }
    fn profile(&self) -> Result<(Vec<usize>, ProfileOptions), String> {
        let overlaps = self
            .overlaps
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<usize>()
                    .map_err(|_| "Supply comma-separated overlap lengths".to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let options = ProfileOptions {
            method: self.method,
            scan_repeated_overlaps: self.scan_repeats,
            max_comparisons: self
                .comparisons
                .parse()
                .map_err(|_| "Invalid comparison budget")?,
            max_occurrences: self
                .occurrences
                .parse()
                .map_err(|_| "Invalid occurrence limit")?,
        };
        options.validate().map_err(|e| e.to_string())?;
        Ok((overlaps, options))
    }
}
#[derive(Clone)]
struct Evidence {
    config: MethodsState,
    revisions: Vec<String>,
    target: String,
    epoch: u64,
    tails: Vec<TailPreview>,
    profile: Option<ProfileReport>,
}
fn guidance(minimum: Option<usize>, maximum: Option<usize>) -> String {
    match (minimum, maximum) {
        (Some(a), Some(b)) => format!("{a}–{b} bp"),
        (Some(a), None) => format!("at least {a} bp; no upper bound specified"),
        _ => "not established".into(),
    }
}
fn short(sequence: &str) -> String {
    let s = sequence.chars().take(300).collect::<String>();
    if s.len() < sequence.len() {
        format!("{s}…")
    } else {
        s
    }
}
fn append_preparation(
    steps: &mut Vec<Operation>,
    input_id: &str,
    input: &TailInput,
) -> Result<Source, String> {
    let source = Source::Input(input_id.into());
    match input.recipe()? {
        TailSubstrate::Digest {
            enzymes,
            fragment_index,
            blunt_before,
        } => {
            steps.push(Operation::Digest {
                template: source,
                enzymes,
                fragment_index,
            });
            if blunt_before {
                let template = Source::Step(steps.len() - 1);
                steps.push(Operation::RepairEnds {
                    template,
                    method: EndRepairMethod::Blunt,
                });
            }
        }
        TailSubstrate::ExactPcr { forward, reverse } => steps.push(Operation::Pcr {
            template: source,
            forward,
            reverse,
        }),
    };
    let template = Source::Step(steps.len() - 1);
    steps.push(Operation::Tail3Prime {
        template,
        base: input.base,
    });
    let mut source = Source::Step(steps.len() - 1);
    let treatment = match input.treatment {
        Treatment::Keep => None,
        Treatment::Phosphorylate => Some(EndTreatment::Phosphorylate),
        Treatment::Dephosphorylate => Some(EndTreatment::Dephosphorylate),
    };
    if let Some(treatment) = treatment {
        steps.push(Operation::TreatEnds {
            template: source,
            treatment,
        });
        source = Source::Step(steps.len() - 1);
    }
    Ok(source)
}
async fn ask(
    worker: &CachedWorker,
    job: AnalysisJob,
    revision: String,
    id: u32,
) -> Result<AnalysisResult, String> {
    let response = worker.ask(Request { id, revision, job }).await?;
    match response.result {
        AnalysisResult::Error(e) => Err(e),
        result => Ok(result),
    }
}
#[component]
pub fn CloningMethodsPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    config: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
) -> impl IntoView {
    let sources = Memo::new(move |_| {
        let c = config.get();
        let indexes = if c.mode == Mode::Tailing {
            vec![c.vector.document, c.insert.document]
        } else {
            c.rows.iter().map(|r| r.document).collect()
        };
        indexes
            .into_iter()
            .map(|i| get_query.run(i))
            .collect::<Option<Vec<_>>>()
    });
    let source_revisions = Memo::new(move |_| {
        sources.with(|s| {
            s.as_ref()
                .map(|m| m.iter().map(Molecule::revision).collect::<Vec<_>>())
        })
    });
    let target = Memo::new(move |_| molecule.with(Molecule::revision));
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (
            config.get(),
            source_revisions.get(),
            target.get(),
            epoch.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1));
        evidence.set(None);
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
        busy.set(false);
        status.set(String::new());
    });
    let compute = Callback::new(move |()| {
        let state = config.get_untracked();
        let Some(inputs) = sources.get_untracked() else {
            status.set("Choose available source documents".into());
            return;
        };
        let recipes = if state.mode == Mode::Tailing {
            match state
                .vector
                .recipe()
                .and_then(|v| state.insert.recipe().map(|i| vec![v, i]))
            {
                Ok(r) => r,
                Err(e) => {
                    status.set(e);
                    return;
                }
            }
        } else {
            vec![]
        };
        let profile = if state.mode == Mode::Overlap {
            match state.profile() {
                Ok(p) => Some(p),
                Err(e) => {
                    status.set(e);
                    return;
                }
            }
        } else {
            None
        };
        let worker = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        status.set("Reconstructing method evidence…".into());
        let bound_target = target.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let alive = alive.clone();
        spawn_local(async move {
            let revs = inputs.iter().map(Molecule::revision).collect::<Vec<_>>();
            let result = async {
                if let Some((overlaps, options)) = profile {
                    let reverse = state.rows.iter().map(|r| r.reverse).collect::<Vec<_>>();
                    let revision = protocol::overlap_profile_revision(
                        &inputs,
                        &reverse,
                        &overlaps,
                        state.circular,
                        &options,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::CheckOverlapProfile {
                            molecules: inputs,
                            reverse,
                            overlaps,
                            circular: state.circular,
                            options,
                        },
                        revision,
                        token,
                    )
                    .await?;
                    if let AnalysisResult::OverlapProfile(report) = response {
                        Ok((vec![], Some(report)))
                    } else {
                        Err("Unexpected overlap profile response".into())
                    }
                } else {
                    let mut previews = vec![];
                    for (i, (input, recipe)) in inputs.into_iter().zip(recipes).enumerate() {
                        let base = if i == 0 {
                            state.vector.base
                        } else {
                            state.insert.base
                        };
                        let result = ask(
                            &worker,
                            AnalysisJob::TailPreview {
                                molecule: input,
                                recipe,
                                base,
                            },
                            revs[i].clone(),
                            token.wrapping_add(i as u32),
                        )
                        .await
                        .map_err(|e| {
                            format!("{}: {e}", if i == 0 { "Vector" } else { "Insert" })
                        })?;
                        if let AnalysisResult::TailPreview(preview) = result {
                            previews.push(preview);
                        } else {
                            return Err("Unexpected tailing preview response".into());
                        }
                    }
                    Ok((previews, None))
                }
            }
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || source_revisions.get_untracked() != Some(revs.clone())
                || target.get_untracked() != bound_target
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match result {
                Ok((tails, profile)) => {
                    evidence.set(Some(Arc::new(Evidence {
                        config: state,
                        revisions: revs,
                        target: bound_target,
                        epoch: bound_epoch,
                        tails,
                        profile,
                    })));
                    status.set("Evidence ready. Review assumptions and build a candidate against the independently selected target.".into());
                }
                Err(e) => status.set(e),
            }
        });
    });
    let review = Callback::new(move |()| {
        let (Some(proof), Some(inputs)) = (evidence.get_untracked(), sources.get_untracked())
        else {
            return;
        };
        let state = config.get_untracked();
        let current_target = molecule.get_untracked();
        if proof.config != state
            || proof.revisions != inputs.iter().map(Molecule::revision).collect::<Vec<_>>()
            || proof.target != current_target.revision()
            || proof.epoch != epoch.get_untracked()
        {
            status.set(
                "Inputs changed; recompute method evidence before reviewing a candidate".into(),
            );
            return;
        }
        let mut steps = vec![];
        let mut bound_inputs = vec![];
        if state.mode == Mode::Tailing {
            for (id, mol) in ["vector", "insert"].into_iter().zip(inputs) {
                bound_inputs.push(PlanInput {
                    id: id.into(),
                    revision: mol.revision(),
                    molecule: mol,
                });
            }
            let vector = match append_preparation(&mut steps, "vector", &state.vector) {
                Ok(s) => s,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let insert = match append_preparation(&mut steps, "insert", &state.insert) {
                Ok(s) => s,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            steps.push(Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: vector,
                        reverse: state.vector.reversed,
                    },
                    OrientedSource {
                        source: insert,
                        reverse: state.insert.reversed,
                    },
                ],
                circular: state.circular,
                allow_nicked: state.allow_nicked,
            });
        } else {
            let (overlaps, options) = match state.profile() {
                Ok(p) => p,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let mut fragments = vec![];
            for (i, (mol, row)) in inputs.into_iter().zip(&state.rows).enumerate() {
                let id = format!("fragment-{}", i + 1);
                fragments.push(OrientedSource {
                    source: Source::Input(id.clone()),
                    reverse: row.reverse,
                });
                bound_inputs.push(PlanInput {
                    id,
                    revision: mol.revision(),
                    molecule: mol,
                });
            }
            steps.push(Operation::ProfiledOverlap {
                fragments,
                overlaps,
                circular: state.circular,
                options,
            });
        }
        on_candidate.run(Candidate::new(DesignPlan {
            name: format!(
                "{} to {}",
                if state.mode == Mode::Tailing {
                    "Declared tailing and ligation"
                } else {
                    "Profiled overlap assembly"
                },
                current_target.name
            ),
            target_revision: current_target.revision(),
            inputs: bound_inputs,
            steps,
        }));
    });
    view! {<div class="inspector-content cloning-methods-panel"><h2>"Cloning methods"</h2><p class="helper">"The active document is the independently chosen target. Method evidence never selects a target or establishes experimental success."</p>
        <label>"Method workflow"<select aria-label="Cloning workflow" prop:value=move||if config.get().mode==Mode::Tailing{"tailing"}else{"overlap"} on:change=move|e|config.update(|c|c.mode=if event_target_value(&e)=="overlap"{Mode::Overlap}else{Mode::Tailing})><option value="tailing">"TA / GC tailing and ligation"</option><option value="overlap">"Gibson / HiFi / In-Fusion overlap"</option></select></label>
        <Show when=move||config.get().mode==Mode::Tailing><p class="helper">"Append one declared nucleotide to both actual 3′ ends of an established blunt duplex. This is an ideal endpoint operation; polymerase activity, mixtures and proprietary kit chemistry are not predicted. TA/GC products include two added junction bases."</p><div class="method-shortcuts"><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::T;c.insert.base=TailingBase::A;})>"Use T/A pair"</button><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::C;c.insert.base=TailingBase::G;})>"Use C/G pair"</button></div><TailFields vector=true config=config query_options=query_options/><TailFields vector=false config=config query_options=query_options/><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained tailing ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Exact PCR has 5′ OH termini. A sealed ligation may need an explicit phosphorylation step. Factory-vector chemistry and host nick repair are not inferred."</p></Show>
        <Show when=move||config.get().mode==Mode::Overlap><label>"Overlap profile"<select aria-label="Overlap profile method" prop:value=move||format!("{:?}",config.get().method) on:change=move|e|config.update(|c|c.method=match event_target_value(&e).as_str(){"Gibson"=>OverlapMethod::Gibson,"InFusion"=>OverlapMethod::InFusion,_=>OverlapMethod::HiFi})><option value="Gibson">"Gibson"</option><option value="HiFi">"HiFi"</option><option value="InFusion">"In-Fusion"</option></select></label><p class="helper">"Fragments are listed in product order with explicit orientation. Homology geometry is checked exactly; published length guidance is advisory. No reaction completion or unique mixture product is established."</p>
            {move||config.get().rows.into_iter().enumerate().map(|(i,row)|view!{<fieldset class="overlap-source-row"><legend>{format!("Fragment {}",i+1)}</legend><label>"Source"<select aria-label=format!("Overlap source {}",i+1) prop:value=move||config.with(|c|c.rows.get(i).map(|r|r.document).unwrap_or(0).to_string()) on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(i){r.document=event_target_value(&e).parse().unwrap_or(0)})>{move||query_options.get().into_iter().map(|(n,label)|view!{<option value=n.to_string()>{label}</option>}).collect_view()}</select></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse overlap source {}",i+1) prop:checked=row.reverse on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(i){r.reverse=event_target_checked(&e)})/>"Reverse complement"</label><div class="method-row-actions"><button class="text-tool" disabled=move||i==0 on:click=move|_|config.update(|c|{if i>0&&i<c.rows.len(){c.rows.swap(i,i-1);}})>"Move up"</button><button class="text-tool" disabled=move||{config.with(|c|i+1>=c.rows.len())} on:click=move|_|config.update(|c|{if i+1<c.rows.len(){c.rows.swap(i,i+1);}})>"Move down"</button><button class="text-tool" disabled=move||config.with(|c|c.rows.len()<=1) on:click=move|_|config.update(|c|{if c.rows.len()>1&&i<c.rows.len(){c.rows.remove(i);}})>"Remove fragment"</button></div></fieldset>}).collect_view()}
            <button class="button secondary full" disabled=move||{config.with(|c|c.rows.len()>=12)} on:click=move|_|config.update(|c|c.rows.push(OverlapInput{document:0,reverse:false}))>"Add overlap fragment"</button><p class="helper">"This panel supports up to 12 ordered fragments; the typed core has its separate documented work bounds."</p><label>"Overlap lengths in junction order (bp)"<input aria-label="Profile overlap lengths" prop:value=move||config.get().overlaps on:input=move|e|config.update(|c|c.overlaps=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Scan repeated overlap words" prop:checked=move||config.get().scan_repeats on:change=move|e|config.update(|c|c.scan_repeats=event_target_checked(&e))/>"Search supplied fragments for repeated homology words"</label><label>"Repeat comparison budget"<input type="number" min="1" aria-label="Overlap repeat budget" prop:value=move||config.get().comparisons on:input=move|e|config.update(|c|c.comparisons=event_target_value(&e))/></label><label>"Maximum reported homology occurrences"<input type="number" min="1" aria-label="Overlap occurrence limit" prop:value=move||config.get().occurrences on:input=move|e|config.update(|c|c.occurrences=event_target_value(&e))/></label>
        </Show>
        <label class="checkbox-label"><input type="checkbox" aria-label="Method product circular" prop:checked=move||config.get().circular on:change=move|e|config.update(|c|c.circular=event_target_checked(&e))/>"Circular candidate product"</label><button class="button primary full" disabled=move||busy.get() on:click=move|_|compute.run(())>"Compute method evidence"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);status.set("Method analysis cancelled".into());}>"Cancel method analysis"</button><p class="method-status" role="status">{move||status.get()}</p>
        {move||evidence.get().map(|proof|{let report=proof.clone();view!{<div class="method-evidence">
            {proof.tails.iter().enumerate().map(|(i,preview)|{let opened=preview.molecule.clone();view!{<article class="tail-evidence"><h3>{format!("{} · {}-tailed duplex",if i==0{"Vector"}else{"Insert"},preview.base.as_str())}</h3><pre>{format!("Top 5′→3′: {}\nBottom 5′→3′: {}\n{}",short(preview.tailing.duplex.top.as_str()),short(preview.tailing.duplex.bottom.as_str()),preview.tailing.duplex.summary().map(|s|s.describe()).unwrap_or_else(|e|e.to_string()))}</pre><details><summary>"Source, additions and assumptions"</summary><pre>{format!("Source top 5′→3′: {}\nSource bottom 5′→3′: {}",short(preview.tailing.source.top.as_str()),short(preview.tailing.source.bottom.as_str()))}</pre>{preview.tailing.additions.iter().map(|a|view!{<p>{format!("New {:?} residue {}–{}: {} (no template ancestry)",a.strand,a.product_start+1,a.product_end,a.sequence)}</p>}).collect_view()}{preview.findings.iter().chain(&preview.tailing.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details><button class="button secondary full" on:click=move|_|match on_products.run(vec![opened.clone()]){Ok(())=>status.set("Reference preview opened; reconstruct its original graph for later physical operations.".into()),Err(e)=>status.set(e)}>{format!("Open {} reference preview",if i==0{"vector"}else{"insert"})}</button></article>}}).collect_view()}
            {proof.profile.as_ref().map(|p|view!{<article class="overlap-profile-evidence"><h3>"Overlap evidence"</h3><p>{format!("Geometry consistent: {} · analysis complete: {} · published length guidance: {}",p.geometry_valid,p.complete,match p.recommendation_compliant{Some(true)=>"met",Some(false)=>"outside published guidance",None=>"not established"})}</p><p class="helper">{format!("{} comparisons; no exhaustive alternative-product search",p.comparisons)}</p>{p.junctions.iter().map(|j|view!{<details class="profile-junction"><summary>{format!("Junction {} · {}bp · {}",j.index+1,j.overlap_bases,if j.exact_match{"exact end homology"}else{"incompatible homology"})}</summary><pre>{format!("Upstream {}\nDownstream {}",j.homology,j.downstream_homology)}</pre><p>{format!("GC: {} · advisory length: {} · repeat search complete: {}",j.gc_fraction.map(|f|format!("{:.1}%",f*100.)).unwrap_or_else(||"unavailable".into()),guidance(j.recommended_min,j.recommended_max),j.repeat_search_complete)}</p>{j.occurrences.iter().take(40).map(|o|view!{<p>{format!("Source {} · {}–{}{}{}",o.fragment+1,o.start+1,o.end,if o.reverse{" reverse"}else{" forward"},if o.assigned_junction{" · assigned end"}else{" · additional occurrence"})}</p>}).collect_view()}<p class="helper">"Up to 40 occurrence positions shown; full evidence download retains every reported position."</p></details>}).collect_view()}<details><summary>"Profile findings and assumptions"</summary>{p.findings.iter().chain(&p.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details></article>})}
            <button class="button primary full" on:click=move|_|review.run(())>"Review method candidate"</button><button class="button secondary full" on:click=move|_|{let result=if let Some(profile)=&report.profile{serde_json::to_vec_pretty(profile)}else{serde_json::to_vec_pretty(&report.tails)};match result{Ok(bytes)=>if let Err(e)=download("cloning-method-evidence.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}}>"Download method evidence"</button>
        </div>}})}
    </div>}
}

#[component]
fn TailFields(
    vector: bool,
    config: RwSignal<MethodsState>,
    query_options: Memo<Vec<(usize, String)>>,
) -> impl IntoView {
    let label = if vector { "Vector" } else { "Insert" };
    view! {<fieldset class="tail-source-fields"><legend>{label}</legend><label>"Source document"<select aria-label=format!("{label} tailing source") prop:value=move||config.with(|c|c.input(vector).document.to_string()) on:change=move|e|config.update(|c|c.input_mut(vector).document=event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label><label>"Preparation"<select aria-label=format!("{label} tailing preparation") prop:value=move||if config.with(|c|c.input(vector).pcr){"pcr"}else{"digest"} on:change=move|e|config.update(|c|c.input_mut(vector).pcr=event_target_value(&e)=="pcr")><option value="digest">"Restriction fragment"</option><option value="pcr">"Unique exact PCR product"</option></select></label>
    <Show when=move||!config.with(|c|c.input(vector).pcr)><label>"Restriction enzymes"<input aria-label=format!("{label} tailing enzymes") prop:value=move||config.with(|c|c.input(vector).enzymes.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).enzymes=event_target_value(&e))/></label><label>"Digest fragment (cut order, 1-based)"<input type="number" min="1" aria-label=format!("{label} tailing fragment") prop:value=move||config.with(|c|c.input(vector).fragment.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).fragment=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("{label} blunt before tailing") prop:checked=move||config.with(|c|c.input(vector).blunt) on:change=move|e|config.update(|c|c.input_mut(vector).blunt=event_target_checked(&e))/>"Blunt staggered ends before tailing"</label></Show>
    <Show when=move||config.with(|c|c.input(vector).pcr)><label>"Forward annealing sequence 5′→3′"<input aria-label=format!("{label} tailing forward primer") prop:value=move||config.with(|c|c.input(vector).forward.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward=event_target_value(&e))/></label><label>"Reverse annealing sequence 5′→3′"<input aria-label=format!("{label} tailing reverse primer") prop:value=move||config.with(|c|c.input(vector).reverse.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse=event_target_value(&e))/></label><label>"Forward primer 5′ extension"<input aria-label=format!("{label} tailing forward extension") prop:value=move||config.with(|c|c.input(vector).forward_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_tail=event_target_value(&e))/></label><label>"Reverse primer 5′ extension"<input aria-label=format!("{label} tailing reverse extension") prop:value=move||config.with(|c|c.input(vector).reverse_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_tail=event_target_value(&e))/></label></Show>
    <label>"One nucleotide added to BOTH 3′ ends"<select aria-label=format!("{label} tailing base") prop:value=move||config.with(|c|c.input(vector).base.as_str()) on:change=move|e|config.update(|c|c.input_mut(vector).base=match event_target_value(&e).as_str(){"C"=>TailingBase::C,"G"=>TailingBase::G,"T"=>TailingBase::T,_=>TailingBase::A})><option value="A">"A"</option><option value="C">"C"</option><option value="G">"G"</option><option value="T">"T"</option></select></label><label>"Explicit terminal treatment after tailing"<select aria-label=format!("{label} tailing terminal treatment") prop:value=move||config.with(|c|match c.input(vector).treatment{Treatment::Keep=>"keep",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"}) on:change=move|e|config.update(|c|c.input_mut(vector).treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep})><option value="keep">"Keep reconstructed 5′ chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse {label} tailing product") prop:checked=move||config.with(|c|c.input(vector).reversed) on:change=move|e|config.update(|c|c.input_mut(vector).reversed=event_target_checked(&e))/>"Reverse fragment for ligation"</label></fieldset>}
}
