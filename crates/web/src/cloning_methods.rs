//! Explicit preparation graphs for TA/GC ligation and oriented overlap profiles.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    Molecule,
    analysis::Primer,
    chemistry::{EndRepairMethod, EndTreatment},
    overlap_profiles::{OverlapMethod, ProfileOptions, ProfileReport},
    protocol::{AnalysisJob, AnalysisResult, Request},
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
    pub forward_name: String,
    pub reverse_name: String,
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
            forward_name: "Forward".into(),
            reverse_name: "Reverse".into(),
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
                forward: Primer::new(&self.forward_name, &self.forward, &self.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&self.reverse_name, &self.reverse, &self.reverse_tail)
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
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OverlapPreparation {
    WholeLinear,
    Digest,
    ExactPcr,
}
#[derive(Clone, PartialEq, Eq)]
pub struct OverlapInput {
    pub document: usize,
    pub reverse: bool,
    pub preparation: OverlapPreparation,
    pub enzymes: String,
    pub fragment: String,
    pub repair: Option<EndRepairMethod>,
    pub forward_name: String,
    pub reverse_name: String,
    pub forward: String,
    pub reverse_primer: String,
    pub forward_tail: String,
    pub reverse_tail: String,
}
impl Default for OverlapInput {
    fn default() -> Self {
        Self {
            document: 0,
            reverse: false,
            preparation: OverlapPreparation::WholeLinear,
            enzymes: "EcoRV".into(),
            fragment: "1".into(),
            repair: None,
            forward_name: "Forward".into(),
            reverse_name: "Reverse".into(),
            forward: String::new(),
            reverse_primer: String::new(),
            forward_tail: String::new(),
            reverse_tail: String::new(),
        }
    }
}
impl OverlapInput {
    fn summary(&self) -> String {
        let preparation = match self.preparation {
            OverlapPreparation::WholeLinear => {
                "Whole linear source (circular DNA requires preparation)".into()
            }
            OverlapPreparation::Digest => format!(
                "Restriction digest: {}; fragment {} (1-based); {}",
                self.enzymes,
                self.fragment,
                match self.repair {
                    None => "cohesive ends retained",
                    Some(EndRepairMethod::FillFivePrimeOverhangs) => "explicit 5′ overhang fill-in",
                    Some(EndRepairMethod::RemoveThreePrimeOverhangs) =>
                        "explicit 3′ overhang removal",
                    Some(EndRepairMethod::Blunt) => "explicit blunting",
                }
            ),
            OverlapPreparation::ExactPcr => format!(
                "Exact PCR: annealing {} / {} bp; 5′ extensions {} / {} bp",
                self.forward.len(),
                self.reverse_primer.len(),
                self.forward_tail.len(),
                self.reverse_tail.len()
            ),
        };
        format!(
            "{} · {}",
            preparation,
            if self.reverse {
                "reverse complement after preparation"
            } else {
                "forward after preparation"
            }
        )
    }
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
            rows: vec![OverlapInput::default(), OverlapInput::default()],
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
    overlap: Option<collagene_core::assembly_recipe::Preview>,
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
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
) -> impl IntoView {
    let source_indices = Memo::new(move |_| {
        config.with(|c| {
            if c.mode == Mode::Tailing {
                vec![c.vector.document, c.insert.document]
            } else {
                c.rows.iter().map(|r| r.document).collect()
            }
        })
    });
    let sources = Memo::new(move |_| {
        source_indices
            .get()
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
        let bound_target_molecule = molecule.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let alive = alive.clone();
        spawn_local(async move {
            let revs = inputs.iter().map(Molecule::revision).collect::<Vec<_>>();
            let result = async {
                if let Some((overlaps, options)) = profile {
                    let _ = (overlaps, options); // shared recipe reuses the exact parsed profile.
                    let (bound_inputs, recipe) =
                        overlap_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::assembly_recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewOverlapRecipe {
                            inputs: bound_inputs,
                            target: bound_target_molecule,
                            recipe,
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::OverlapRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((vec![], preview.profile.clone(), Some(preview)))
                    } else {
                        Err("Unexpected overlap preparation response".into())
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
                    Ok((previews, None, None))
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
                Ok((tails, profile, overlap)) => {
                    evidence.set(Some(Arc::new(Evidence {
                        config: state,
                        revisions: revs,
                        target: bound_target,
                        epoch: bound_epoch,
                        tails,
                        profile,
                        overlap,
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
            let (bound, recipe) = match overlap_request(&state, inputs, &current_target) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            match collagene_core::assembly_recipe::compile(&bound, &current_target, &recipe) {
                Ok(compiled) => {
                    on_candidate.run(Candidate::new(compiled.plan));
                    return;
                }
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            }
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
        <Show when=move||config.get().mode==Mode::Tailing><p class="helper">"Append one declared nucleotide to both actual 3′ ends of an established blunt duplex. This is an ideal endpoint operation; polymerase activity, mixtures and proprietary kit chemistry are not predicted. TA/GC products include two added junction bases."</p><div class="method-shortcuts"><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::T;c.insert.base=TailingBase::A;})>"Use T/A pair"</button><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::C;c.insert.base=TailingBase::G;})>"Use C/G pair"</button></div><TailFields vector=true config=config query_options=query_options primer_sources list_primers resolve_primer epoch/><TailFields vector=false config=config query_options=query_options primer_sources list_primers resolve_primer epoch/><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained tailing ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Exact PCR has 5′ OH termini. A sealed ligation may need an explicit phosphorylation step. Factory-vector chemistry and host nick repair are not inferred."</p></Show>
        <Show when=move||config.get().mode==Mode::Overlap><label>"Overlap profile"<select aria-label="Overlap profile method" prop:value=move||format!("{:?}",config.get().method) on:change=move|e|config.update(|c|c.method=match event_target_value(&e).as_str(){"Gibson"=>OverlapMethod::Gibson,"InFusion"=>OverlapMethod::InFusion,_=>OverlapMethod::HiFi})><option value="Gibson">"Gibson"</option><option value="HiFi">"HiFi"</option><option value="InFusion">"In-Fusion"</option></select></label><p class="helper">"Fragments are listed in product order with explicit orientation. Homology geometry is checked exactly; published length guidance is advisory. No reaction completion or unique mixture product is established."</p>
            <For each=move||config.with(|c|(0..c.rows.len()).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{<OverlapFields config index=i query_options primer_sources list_primers resolve_primer epoch/>}/>
            <button class="button secondary full" disabled=move||{config.with(|c|c.rows.len()>=12)} on:click=move|_|config.update(|c|c.rows.push(OverlapInput::default()))>"Add overlap fragment"</button><p class="helper">"This panel supports up to 12 ordered fragments; the typed core has its separate documented work bounds."</p><label>"Overlap lengths in junction order (bp)"<input aria-label="Profile overlap lengths" prop:value=move||config.get().overlaps on:input=move|e|config.update(|c|c.overlaps=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Scan repeated overlap words" prop:checked=move||config.get().scan_repeats on:change=move|e|config.update(|c|c.scan_repeats=event_target_checked(&e))/>"Search supplied fragments for repeated homology words"</label><label>"Repeat comparison budget"<input type="number" min="1" aria-label="Overlap repeat budget" prop:value=move||config.get().comparisons on:input=move|e|config.update(|c|c.comparisons=event_target_value(&e))/></label><label>"Maximum reported homology occurrences"<input type="number" min="1" aria-label="Overlap occurrence limit" prop:value=move||config.get().occurrences on:input=move|e|config.update(|c|c.occurrences=event_target_value(&e))/></label>
        </Show>
        <label class="checkbox-label"><input type="checkbox" aria-label="Method product circular" prop:checked=move||config.get().circular on:change=move|e|config.update(|c|c.circular=event_target_checked(&e))/>"Circular candidate product"</label><button class="button primary full" disabled=move||busy.get() on:click=move|_|compute.run(())>"Compute method evidence"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);status.set("Method analysis cancelled".into());}>"Cancel method analysis"</button><p class="method-status" role="status">{move||status.get()}</p>
        {move||evidence.get().map(|proof|{let report=proof.clone();view!{<div class="method-evidence">
            {proof.tails.iter().enumerate().map(|(i,preview)|{let opened=preview.molecule.clone();view!{<article class="tail-evidence"><h3>{format!("{} · {}-tailed duplex",if i==0{"Vector"}else{"Insert"},preview.base.as_str())}</h3><pre>{format!("Top 5′→3′: {}\nBottom 5′→3′: {}\n{}",short(preview.tailing.duplex.top.as_str()),short(preview.tailing.duplex.bottom.as_str()),preview.tailing.duplex.summary().map(|s|s.describe()).unwrap_or_else(|e|e.to_string()))}</pre><details><summary>"Source, additions and assumptions"</summary><pre>{format!("Source top 5′→3′: {}\nSource bottom 5′→3′: {}",short(preview.tailing.source.top.as_str()),short(preview.tailing.source.bottom.as_str()))}</pre>{preview.tailing.additions.iter().map(|a|view!{<p>{format!("New {:?} residue {}–{}: {} (no template ancestry)",a.strand,a.product_start+1,a.product_end,a.sequence)}</p>}).collect_view()}{preview.findings.iter().chain(&preview.tailing.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details><button class="button secondary full" on:click=move|_|match on_products.run(vec![opened.clone()]){Ok(())=>status.set("Reference preview opened; reconstruct its original graph for later physical operations.".into()),Err(e)=>status.set(e)}>{format!("Open {} reference preview",if i==0{"vector"}else{"insert"})}</button></article>}}).collect_view()}
            {proof.overlap.as_ref().map(|preview|view!{<article class="overlap-preparation-evidence"><h3>"Prepared fragments from original sources"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {:?}",preview.preparation_complete,preview.verification.status)}</p><p class="helper">"Prepared reference sequences below are reconstructed from the declared originals. They are not substituted as plan inputs, confirmed amplicons or physical reaction outcomes."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source 5′→3′ (before orientation): {}",short(fragment.sequence.as_str()))}</pre></details>}).collect_view()}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.profile.as_ref().map(|p|view!{<article class="overlap-profile-evidence"><h3>"Overlap evidence"</h3><p>{format!("Geometry consistent: {} · analysis complete: {} · published length guidance: {}",p.geometry_valid,p.complete,match p.recommendation_compliant{Some(true)=>"met",Some(false)=>"outside published guidance",None=>"not established"})}</p><p class="helper">{format!("{} comparisons; no exhaustive alternative-product search",p.comparisons)}</p>{p.junctions.iter().map(|j|view!{<details class="profile-junction"><summary>{format!("Junction {} · {}bp · {}",j.index+1,j.overlap_bases,if j.exact_match{"exact end homology"}else{"incompatible homology"})}</summary><pre>{format!("Upstream {}\nDownstream {}",j.homology,j.downstream_homology)}</pre><p>{format!("GC: {} · advisory length: {} · repeat search complete: {}",j.gc_fraction.map(|f|format!("{:.1}%",f*100.)).unwrap_or_else(||"unavailable".into()),guidance(j.recommended_min,j.recommended_max),j.repeat_search_complete)}</p>{j.occurrences.iter().take(40).map(|o|view!{<p>{format!("Source {} · {}–{}{}{}",o.fragment+1,o.start+1,o.end,if o.reverse{" reverse"}else{" forward"},if o.assigned_junction{" · assigned end"}else{" · additional occurrence"})}</p>}).collect_view()}<p class="helper">"Up to 40 occurrence positions shown; full evidence download retains every reported position."</p></details>}).collect_view()}<details><summary>"Profile findings and assumptions"</summary>{p.findings.iter().chain(&p.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details></article>})}
            <button class="button primary full" on:click=move|_|review.run(())>"Review method candidate"</button><button class="button secondary full" on:click=move|_|{let result=if let Some(preview)=&report.overlap{serde_json::to_vec_pretty(preview)}else if let Some(profile)=&report.profile{serde_json::to_vec_pretty(profile)}else{serde_json::to_vec_pretty(&report.tails)};match result{Ok(bytes)=>if let Err(e)=download("cloning-method-evidence.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}}>"Download method evidence"</button>
        </div>}})}
    </div>}
}

#[component]
fn TailFields(
    vector: bool,
    config: RwSignal<MethodsState>,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let label = if vector { "Vector" } else { "Insert" };
    view! {<fieldset class="tail-source-fields"><legend>{label}</legend><label>"Source document"<select aria-label=format!("{label} tailing source") prop:value=move||config.with(|c|c.input(vector).document.to_string()) on:change=move|e|config.update(|c|c.input_mut(vector).document=event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label><label>"Preparation"<select aria-label=format!("{label} tailing preparation") prop:value=move||if config.with(|c|c.input(vector).pcr){"pcr"}else{"digest"} on:change=move|e|config.update(|c|c.input_mut(vector).pcr=event_target_value(&e)=="pcr")><option value="digest">"Restriction fragment"</option><option value="pcr">"Unique exact PCR product"</option></select></label>
        <Show when=move||!config.with(|c|c.input(vector).pcr)><label>"Restriction enzymes"<input aria-label=format!("{label} tailing enzymes") prop:value=move||config.with(|c|c.input(vector).enzymes.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).enzymes=event_target_value(&e))/></label><label>"Digest fragment (cut order, 1-based)"<input type="number" min="1" aria-label=format!("{label} tailing fragment") prop:value=move||config.with(|c|c.input(vector).fragment.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).fragment=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("{label} blunt before tailing") prop:checked=move||config.with(|c|c.input(vector).blunt) on:change=move|e|config.update(|c|c.input_mut(vector).blunt=event_target_checked(&e))/>"Blunt staggered ends before tailing"</label></Show>
        <Show when=move||config.with(|c|c.input(vector).pcr)>
        <crate::primer_picker::PrimerPicker label=format!("{label} tailing forward") sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|{let r=c.input_mut(vector);r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
        <label>"Forward primer name"<input maxlength="256" aria-label=format!("{label} tailing forward primer name") prop:value=move||config.with(|c|c.input(vector).forward_name.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_name=event_target_value(&e))/></label>
        <crate::primer_picker::PrimerPicker label=format!("{label} tailing reverse") sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|{let r=c.input_mut(vector);r.reverse_name=primer.name;r.reverse=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
        <label>"Reverse primer name"<input maxlength="256" aria-label=format!("{label} tailing reverse primer name") prop:value=move||config.with(|c|c.input(vector).reverse_name.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_name=event_target_value(&e))/></label>
    <label>"Forward annealing sequence 5′→3′"<input aria-label=format!("{label} tailing forward primer") prop:value=move||config.with(|c|c.input(vector).forward.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward=event_target_value(&e))/></label><label>"Reverse annealing sequence 5′→3′"<input aria-label=format!("{label} tailing reverse primer") prop:value=move||config.with(|c|c.input(vector).reverse.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse=event_target_value(&e))/></label><label>"Forward primer 5′ extension"<input aria-label=format!("{label} tailing forward extension") prop:value=move||config.with(|c|c.input(vector).forward_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_tail=event_target_value(&e))/></label><label>"Reverse primer 5′ extension"<input aria-label=format!("{label} tailing reverse extension") prop:value=move||config.with(|c|c.input(vector).reverse_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_tail=event_target_value(&e))/></label><details><summary>"Review current primer domains"</summary><pre>{move||config.with(|c|{let r=c.input(vector);format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse,r.reverse,r.reverse_tail)})}</pre></details></Show>
        <label>"One nucleotide added to BOTH 3′ ends"<select aria-label=format!("{label} tailing base") prop:value=move||config.with(|c|c.input(vector).base.as_str()) on:change=move|e|config.update(|c|c.input_mut(vector).base=match event_target_value(&e).as_str(){"C"=>TailingBase::C,"G"=>TailingBase::G,"T"=>TailingBase::T,_=>TailingBase::A})><option value="A">"A"</option><option value="C">"C"</option><option value="G">"G"</option><option value="T">"T"</option></select></label><label>"Explicit terminal treatment after tailing"<select aria-label=format!("{label} tailing terminal treatment") prop:value=move||config.with(|c|match c.input(vector).treatment{Treatment::Keep=>"keep",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"}) on:change=move|e|config.update(|c|c.input_mut(vector).treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep})><option value="keep">"Keep reconstructed 5′ chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse {label} tailing product") prop:checked=move||config.with(|c|c.input(vector).reversed) on:change=move|e|config.update(|c|c.input_mut(vector).reversed=event_target_checked(&e))/>"Reverse fragment for ligation"</label></fieldset>}
}

#[component]
fn OverlapFields(
    config: RwSignal<MethodsState>,
    index: usize,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    view! {<fieldset class="overlap-source-row"><legend>{format!("Fragment {}",index+1)}</legend>
        <label>"Original source"<select aria-label=format!("Overlap source {}",index+1) prop:value=move||row().document.to_string() on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.document=event_target_value(&e).parse().unwrap_or(0)})>{move||query_options.get().into_iter().map(|(n,label)|view!{<option value=n.to_string()>{label}</option>}).collect_view()}</select></label>
        <label>"Prepare fragment"<select aria-label=format!("Overlap preparation {}",index+1) prop:value=move||match row().preparation{OverlapPreparation::WholeLinear=>"whole",OverlapPreparation::Digest=>"digest",OverlapPreparation::ExactPcr=>"pcr"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.preparation=match event_target_value(&e).as_str(){"digest"=>OverlapPreparation::Digest,"pcr"=>OverlapPreparation::ExactPcr,_=>OverlapPreparation::WholeLinear}})><option value="whole">"Whole linear DNA"</option><option value="digest">"Restriction digest fragment"</option><option value="pcr">"Exact PCR product"</option></select></label>
        <Show when=move||row().preparation==OverlapPreparation::WholeLinear><p class="helper">"Use an already linear source. Circular DNA must first be prepared by restriction digest or exact PCR; its origin is not an assumed cut."</p></Show>
        <Show when=move||row().preparation==OverlapPreparation::Digest>
            <label>"Restriction enzymes (comma-separated)"<input aria-label=format!("Overlap digest enzymes {}",index+1) prop:value=move||row().enzymes on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.enzymes=event_target_value(&e)})/></label>
            <label>"Digest fragment number (1-based)"<input type="number" min="1" aria-label=format!("Overlap digest fragment {}",index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.fragment=event_target_value(&e)})/></label>
            <label>"Optional explicit end repair"<select aria-label=format!("Overlap end repair {}",index+1) prop:value=move||match row().repair{None=>"none",Some(EndRepairMethod::FillFivePrimeOverhangs)=>"fill",Some(EndRepairMethod::RemoveThreePrimeOverhangs)=>"remove",Some(EndRepairMethod::Blunt)=>"blunt"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.repair=match event_target_value(&e).as_str(){"fill"=>Some(EndRepairMethod::FillFivePrimeOverhangs),"remove"=>Some(EndRepairMethod::RemoveThreePrimeOverhangs),"blunt"=>Some(EndRepairMethod::Blunt),_=>None}})><option value="none">"None — retain cohesive ends"</option><option value="fill">"Fill 5′ overhangs"</option><option value="remove">"Remove 3′ overhangs"</option><option value="blunt">"Blunt ends"</option></select></label>
            <p class="helper">"Fragments are reconstructed in core digest order using the stated complete-cleavage assumption. This is not evidence of enzyme activity, methylation permission or reagent compatibility."</p>
        </Show>
        <Show when=move||row().preparation==OverlapPreparation::ExactPcr>
            <crate::primer_picker::PrimerPicker label=format!("Overlap forward {}",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
            <label>"Forward primer name"<input maxlength="256" aria-label=format!("Overlap forward primer name {}",index+1) prop:value=move||row().forward_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=event_target_value(&e)})/></label>
            <crate::primer_picker::PrimerPicker label=format!("Overlap reverse {}",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=primer.name;r.reverse_primer=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
            <label>"Reverse primer name"<input maxlength="256" aria-label=format!("Overlap reverse primer name {}",index+1) prop:value=move||row().reverse_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=event_target_value(&e)})/></label>

            <label>"Forward annealing sequence 5′→3′"<input aria-label=format!("Overlap forward primer {}",index+1) prop:value=move||row().forward on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward=event_target_value(&e)})/></label>
            <label>"Forward primer 5′ extension"<input aria-label=format!("Overlap forward extension {}",index+1) prop:value=move||row().forward_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_tail=event_target_value(&e)})/></label>
            <label>"Reverse annealing sequence 5′→3′"<input aria-label=format!("Overlap reverse primer {}",index+1) prop:value=move||row().reverse_primer on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_primer=event_target_value(&e)})/></label>
            <label>"Reverse primer 5′ extension"<input aria-label=format!("Overlap reverse extension {}",index+1) prop:value=move||row().reverse_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_tail=event_target_value(&e)})/></label>
            <details><summary>"Review current primer domains"</summary><pre>{move||{let r=row();format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse_primer,r.reverse_primer,r.reverse_tail)}}</pre></details>
            <p class="helper">"Enter annealing domains and added tails separately. Both oligos are written 5′→3′; no primer design or unique amplicon is assumed."</p>
        </Show>
        <label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse overlap source {}",index+1) prop:checked=move||row().reverse on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse=event_target_checked(&e)})/>"Reverse complement prepared fragment"</label>
        <p class="helper">{move||row().summary()}</p>
        <div class="method-row-actions"><button class="text-tool" disabled=move||index==0 on:click=move|_|config.update(|c|{if index>0&&index<c.rows.len(){c.rows.swap(index,index-1);}})>"Move up"</button><button class="text-tool" disabled=move||config.with(|c|index+1>=c.rows.len()) on:click=move|_|config.update(|c|{if index+1<c.rows.len(){c.rows.swap(index,index+1);}})>"Move down"</button><button class="text-tool" disabled=move||config.with(|c|c.rows.len()<=1) on:click=move|_|config.update(|c|{if c.rows.len()>1&&index<c.rows.len(){c.rows.remove(index);}})>"Remove fragment"</button></div>
    </fieldset>}
}

fn overlap_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<(Vec<PlanInput>, collagene_core::assembly_recipe::Recipe), String> {
    use collagene_core::assembly_recipe::{FragmentSpec, Preparation, Recipe};
    if inputs.len() != state.rows.len() || inputs.is_empty() || inputs.len() > 12 {
        return Err("Choose 1–12 ordered original sources".into());
    }
    let (overlaps, options) = state.profile()?;
    let ids = (0..inputs.len())
        .map(|i| format!("source-{}", i + 1))
        .collect::<Vec<_>>();
    let mut fragments = Vec::with_capacity(inputs.len());
    for (i, row) in state.rows.iter().enumerate() {
        if [
            &row.enzymes,
            &row.forward,
            &row.reverse_primer,
            &row.forward_tail,
            &row.reverse_tail,
        ]
        .iter()
        .any(|s| s.len() > 65536)
        {
            return Err("Preparation fields exceed the 64 KiB field limit".into());
        }
        let preparation = match row.preparation {
            OverlapPreparation::WholeLinear => Preparation::WholeLinear,
            OverlapPreparation::Digest => {
                let enzymes = row
                    .enzymes
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect::<Vec<_>>();
                if enzymes.is_empty() || enzymes.len() > 64 || enzymes.iter().any(|e| e.len() > 256)
                {
                    return Err("Choose 1–64 explicit restriction enzyme names".into());
                }
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                Preparation::Digest {
                    enzymes,
                    fragment_index,
                    repair: row.repair,
                }
            }
            OverlapPreparation::ExactPcr => Preparation::ExactPcr {
                forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&row.reverse_name, &row.reverse_primer, &row.reverse_tail)
                    .map_err(|e| e.to_string())?,
            },
        };
        fragments.push(FragmentSpec {
            role: format!("fragment-{}", i + 1),
            input_id: ids[i].clone(),
            preparation,
            reverse: row.reverse,
        });
    }
    let recipe = Recipe {
        name: format!("Profiled overlap assembly to {}", target.name),
        fragments,
        overlaps,
        circular: state.circular,
        options,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::assembly_recipe::source_preflight(&borrowed, target, &recipe)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, recipe))
}
