//! Explicit downstream graph treatment. A sealed graph is not a checked target
//! or a duplex projection; all untreated components and limitations remain.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request as WireRequest},
    run_lineage::EventBinding,
    topo::{
        self,
        physical::{self, nick_sealing},
    },
    topo_nick_sealing,
};
use leptos::prelude::*;
use std::{
    io::{self, Write},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, PartialEq, Eq)]
pub struct Frame {
    pub epoch: u64,
    pub sources: Vec<(String, EventBinding)>,
    pub target: EventBinding,
    pub inspection_revision: String,
    pub inspection_generation: u32,
    pub orientation: usize,
}
pub struct Capture {
    pub inputs: Vec<topo::Input>,
    pub request: topo::Request,
    pub nicks: Vec<physical::Nick>,
}
#[derive(Clone, PartialEq, Eq)]
struct TreatmentForm {
    selected: Vec<bool>,
    assumption: String,
}
struct Evidence {
    frame: Frame,
    form: TreatmentForm,
    report: Arc<topo_nick_sealing::Report>,
}
struct Json(Vec<u8>);
impl Write for Json {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > 32 * 1024 * 1024)
        {
            return Err(io::Error::other("Nick treatment evidence exceeds 32 MiB"));
        }
        self.0
            .try_reserve(bytes.len())
            .map_err(|_| io::Error::other("Evidence allocation failed"))?;
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
#[component]
fn EvidencePage(label: &'static str, total: usize, row: Callback<usize, String>) -> impl IntoView {
    let page = RwSignal::new(0usize);
    view! { <details><summary>{format!("{label} ({total})")}</summary>
        <p>{move || format!("Page {} of {}", page.get()+1,total.div_ceil(40).max(1))}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous"</button>
        <button class="text-tool" disabled=move||(page.get()+1)*40>=total on:click=move|_|page.update(|p|*p+=1)>"Next"</button>
        {move||(page.get()*40..((page.get()+1)*40).min(total)).map(|index|view!{<pre>{row.run(index)}</pre>}).collect_view()}
    </details> }
}
fn call_text(call: physical::Call) -> String {
    format!("source strand {}, call {}", call.strand, call.index)
}
fn report_view(report: Arc<topo_nick_sealing::Report>) -> impl IntoView {
    // Each list renders only its current bounded page. Full evidence remains
    // available explicitly; no prefix is presented as the whole graph.
    let steps = report.clone();
    let components = report.clone();
    let nicks = report.clone();
    let unpaired = report.clone();
    let chains = report.clone();
    let activations = report.clone();
    let assumptions = report.clone();
    view! { <section>
        <h4>"Downstream graph — read-only evidence"</h4>
        <p>"This treatment does not project a duplex, compare a target or create an accepted cloning candidate. Flaps, remaining nicks, terminal activation and other components remain part of the graph."</p>
        <p>{format!("{} components · {} remaining nicks · {} unpaired runs",report.sealing.product.components.len(),report.sealing.product.nicks.len(),report.sealing.product.unpaired.len())}</p>
        <p>{format!("All selected nicks sealed: {} · unsupported or unknown: {} · original transfer work complete: {}",report.sealing.all_requested_nicks_sealed,report.sealing.unsupported_or_unknown,report.sealing.parent.work_complete)}</p>
        <EvidencePage label="Selected nick outcomes" total=report.sealing.steps.len() row=Callback::new(move|i:usize|{ let s=&steps.sealing.steps[i]; format!("Nick {} · {:?} · {}\n3′ {:?}; 5′ {:?}\nPaired template calls: {} / {}",i+1,s.outcome,s.reason,s.three_prime,s.five_prime,s.template_three.map(call_text).unwrap_or_else(||"not established".into()),s.template_five.map(call_text).unwrap_or_else(||"not established".into())) })/>
        <EvidencePage label="All components" total=report.sealing.product.components.len() row=Callback::new(move|i:usize|{let c=&components.sealing.product.components[i];format!("Component {i}: {} paired calls · {} unpaired calls · {} chains · covalently closed: {}",c.paired_nucleotides,c.unpaired_nucleotides,c.chains.len(),c.covalently_closed)})/>
        <EvidencePage label="Remaining nicks and chemistry" total=report.sealing.product.nicks.len() row=Callback::new(move|i:usize|{let n=&nicks.sealing.product.nicks[i];format!("Nick {i}: {} ({:?}) → {} ({:?})",call_text(n.three),n.three_prime,call_text(n.five),n.five_prime)})/>
        <EvidencePage label="Unpaired runs and flaps" total=report.sealing.product.unpaired.len() row=Callback::new(move|i:usize|{let r=&unpaired.sealing.product.unpaired[i];format!("Run {i}: chain {} · calls {}–{} · 5′ terminal: {} · 3′ terminal: {} · displaced by transfers {:?}",r.chain,r.start,r.end,r.terminal_five,r.terminal_three,r.displaced_by)})/>
        <EvidencePage label="All chain ends — residual activation" total=report.sealing.product.chains.len() row=Callback::new(move|i:usize|format!("Chain {i}: {:?}",chains.sealing.product.chains[i].ends))/>
        <EvidencePage label="Treatment assumptions" total=report.sealing.assumptions.len() row=Callback::new(move|i:usize|assumptions.sealing.assumptions[i].clone())/>
        <EvidencePage label="Released activation identities" total=report.sealing.product.released_activations.len() row=Callback::new(move|i:usize|activations.sealing.product.released_activations[i].clone())/>
    </section> }
}
#[component]
pub fn TopoNickView(
    frame: Memo<Option<Frame>>,
    capture: Callback<Frame, Result<Capture, String>>,
    locked: Memo<bool>,
) -> impl IntoView {
    let captured = RwSignal::new(None::<(Frame, Arc<Capture>)>);
    let form = RwSignal::new(TreatmentForm {
        selected: vec![],
        assumption: String::new(),
    });
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cancel = Callback::new(move |()| {
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(value) = worker.get_value() {
            value.close();
        }
        worker.set_value(None);
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = frame.get();
        let _ = locked.get();
        cancel.run(());
        captured.set(None);
        evidence.set(None);
        form.set(TreatmentForm {
            selected: vec![],
            assumption: String::new(),
        });
        status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = form.get();
        cancel.run(());
        evidence.set(None);
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(value) = worker.get_value() {
            value.close();
        }
        worker.set_value(None);
    });
    let load = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some(key) = frame.get_untracked() else {
            status.set("Choose a freshly inspected orientation first.".into());
            return;
        };
        match capture.run(key.clone()) {
            Ok(value) => {
                cancel.run(());
                evidence.set(None);
                form.set(TreatmentForm {
                    selected: vec![false; value.nicks.len()],
                    assumption: String::new(),
                });
                captured.set(Some((key, Arc::new(value))));
                status.set("No nicks selected. Select explicit source-call pairs and declare the ligase assumption.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let review_alive = alive.clone();
    let review = Callback::new(move |()| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let Some((key, value)) = captured.get_untracked() else {
            return;
        };
        if frame.get_untracked().as_ref() != Some(&key) {
            return;
        }
        let chosen = form.get_untracked();
        if chosen.assumption.trim().is_empty() || chosen.assumption.len() > 4096 {
            status.set("Declare a ligase assumption of 1–4096 bytes.".into());
            return;
        }
        if chosen.selected.len() != value.nicks.len() || !chosen.selected.iter().any(|v| *v) {
            status.set("Select at least one actual nick explicitly.".into());
            return;
        }
        // Re-attest the parent bound immediately before submitting any source
        // clones; a captured report is never accepted as treatment authority.
        let fresh = match capture.run(key.clone()) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if fresh.nicks != value.nicks {
            status.set("The actual nick list changed. Reload this path.".into());
            return;
        }
        let value = Arc::new(fresh);
        let request = topo_nick_sealing::Request {
            preparation: value.request.clone(),
            sealing: nick_sealing::Request {
                nicks: value
                    .nicks
                    .iter()
                    .zip(&chosen.selected)
                    .filter(|(_, yes)| **yes)
                    .map(|(nick, _)| nick_sealing::NickSelection {
                        three: nick.three,
                        five: nick.five,
                    })
                    .collect(),
                ligase_assumption: chosen.assumption.clone(),
            },
        };
        let revision = match topo_nick_sealing::revision(&value.inputs, &request) {
            Ok(r) => r,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        cancel.run(());
        let executor = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        worker.set_value(Some(executor.clone()));
        busy.set(true);
        evidence.set(None);
        let ticket = generation.get_untracked();
        let alive = review_alive.clone();
        status.set("Reviewing explicit nick treatment against the original sources…".into());
        spawn_local(async move {
            let response = executor
                .ask(WireRequest {
                    id: ticket,
                    revision: revision.clone(),
                    job: AnalysisJob::SealTopoNicks {
                        inputs: value.inputs.clone(),
                        request,
                    },
                })
                .await;
            executor.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            if locked.get_untracked()
                || frame.get_untracked().as_ref() != Some(&key)
                || form.get_untracked() != chosen
            {
                return;
            }
            busy.set(false);
            worker.set_value(None);
            match response {
                Ok(response) if response.id == ticket && response.revision == revision => {
                    match response.result {
                        AnalysisResult::TopoNicksSealed(report) => {
                            match report.validate_binding(&revision) {
                                Ok(()) => {
                                    evidence.set(Some(Arc::new(Evidence {
                                        frame: key,
                                        form: chosen,
                                        report: Arc::new(report),
                                    })));
                                    status.set("Nick treatment graph reviewed. No target acceptance or duplex handoff was performed.".into());
                                }
                                Err(error) => status.set(error.to_string()),
                            }
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Unexpected nick treatment response.".into()),
                    }
                }
                Ok(_) => status.set("Nick treatment response has a different binding.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let current = move || {
        evidence.get().filter(|e| {
            frame.get().as_ref() == Some(&e.frame) && form.get() == e.form && !locked.get()
        })
    };
    let export = Callback::new(move |()| {
        let Some(e) = current() else {
            return;
        };
        if let Err(error) = capture.run(e.frame.clone()) {
            status.set(error);
            return;
        }
        let mut bytes = Json(vec![]);
        let result = serde_json::to_writer(&mut bytes, e.report.as_ref())
            .map_err(|e| e.to_string())
            .and_then(|()| download("topo-nick-treatment.json", "application/json", &bytes.0));
        if let Err(error) = result {
            status.set(error);
        }
    });
    let nick_page = RwSignal::new(0usize);
    Effect::new(move |_| {
        let _ = frame.get();
        nick_page.set(0);
    });
    view! { <details class="topo-nick-treatment"><summary>"Downstream nick-sealing graph review"</summary>
        <p class="helper">"Explicitly seal supported paired nicks only. This does not trim flaps, phosphorylate ends, replace activation, choose a target or establish laboratory success."</p>
        <button class="text-tool" disabled=move||locked.get()||busy.get()||frame.get().is_none() on:click=move|_|load.run(())>"Load nicks from selected path"</button>
        {move||captured.get().filter(|(key,_)|frame.get().as_ref()==Some(key)).map(|(_,value)| {
            let total=value.nicks.len();
            view!{<section><p>{format!("{total} actual nicks; none selected by default")}</p>
                <button class="text-tool" disabled=move||nick_page.get()==0 on:click=move|_|nick_page.update(|p|*p=p.saturating_sub(1))>"Previous nicks"</button>
                <button class="text-tool" disabled=move||(nick_page.get()+1)*40>=total on:click=move|_|nick_page.update(|p|*p+=1)>"Next nicks"</button>
                {move||(nick_page.get()*40..((nick_page.get()+1)*40).min(total)).map(|i| {let text=format!("Nick {i}: {} ({:?}) → {} ({:?})",call_text(value.nicks[i].three),value.nicks[i].three_prime,call_text(value.nicks[i].five),value.nicks[i].five_prime);view!{<label class="checkbox-label"><input type="checkbox" disabled=move||locked.get()||busy.get() prop:checked=move||form.with(|f|f.selected.get(i).copied().unwrap_or(false)) on:change=move|e|form.update(|f|{if let Some(v)=f.selected.get_mut(i){*v=event_target_checked(&e)}})/>{text}</label>}}).collect_view()}
            </section>}
        })}
        <label>"Explicit ligase assumption"<textarea rows="2" maxlength="4096" disabled=move||locked.get()||busy.get() prop:value=move||form.with(|f|f.assumption.clone()) on:input=move|e|form.update(|f|f.assumption=event_target_value(&e))/></label>
        <button class="button secondary" disabled=move||locked.get()||busy.get()||captured.get().is_none() on:click=move|_|review.run(())>"Review selected nick treatment"</button>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{cancel.run(());status.set("Nick treatment cancelled.".into());}>"Cancel nick treatment"</button></Show>
        <p role="status">{move||status.get()}</p>
        {move||current().map(|e|report_view(e.report.clone()))}
        <button class="text-tool" disabled=move||current().is_none() on:click=move|_|export.run(())>"Download full nick treatment evidence JSON"</button>
    </details> }
}
