//! Snapshot diagnostics over the existing native pairwise engine, without history events.
use crate::{
    assembly_comparison_model as model,
    assembly_product_model::Snapshot,
    assembly_product_view::ReferenceCanvas,
    pairwise_pages::{PageIndex, render_page},
    worker::CachedWorker,
};
use collagene_core::{
    alignment::ComparisonOptions,
    pairwise_review,
    protocol::{AnalysisJob, AnalysisResult, Request},
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
struct Bound {
    capture: model::Capture,
    report: pairwise_review::Report,
    index: PageIndex,
}
fn limit(value: &str, name: &str) -> Result<usize, String> {
    if value.len() > 20 || !value.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Enter a positive whole number for {name}."));
    }
    value
        .parse::<usize>()
        .ok()
        .filter(|n| *n > 0)
        .ok_or_else(|| format!("Enter a positive whole number for {name}."))
}
fn clip(value: &str) -> String {
    if value.len() > 100 {
        format!("{}… ({} bases)", &value[..100], value.len())
    } else if value.is_empty() {
        "∅".into()
    } else {
        value.into()
    }
}
#[component]
pub(crate) fn AssemblyComparison(
    snapshot: Arc<Snapshot>,
    guard: Callback<(), Result<(), String>>,
    current: Memo<bool>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let snapshot = StoredValue::new(snapshot);
    let rotate = RwSignal::new(true);
    let reverse = RwSignal::new(false);
    let defaults = ComparisonOptions::default();
    let cells = RwSignal::new(defaults.alignment.max_cells.to_string());
    let total = RwSignal::new(defaults.max_total_cells.to_string());
    let origins = RwSignal::new(defaults.max_origins.to_string());
    let status = RwSignal::new(String::new());
    let proof = RwSignal::new(None::<Arc<Bound>>);
    let busy = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let differences_page = RwSignal::new(0usize);
    let selected = RwSignal::new(None::<usize>);
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
    let settings = Memo::new(move |_| {
        (
            rotate.get(),
            reverse.get(),
            cells.get(),
            total.get(),
            origins.get(),
            current.get(),
            epoch.get(),
        )
    });
    Effect::new(move |_| {
        let _ = settings.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        proof.set(None);
        busy.set(false);
        page.set(0);
        differences_page.set(0);
        selected.set(None);
    });
    let compute = move |_| {
        let captured = (|| {
            guard.run(())?;
            if !current.get_untracked() {
                return Err("Workspace changed; reconstruct the product again.".into());
            }
            let mut options = ComparisonOptions::default();
            options.both_query_strands = reverse.get_untracked();
            options.reference_rotations = rotate.get_untracked();
            options.alignment.max_cells = limit(&cells.get_untracked(), "alignment cells")?;
            options.max_total_cells = limit(&total.get_untracked(), "total cells")?;
            options.max_origins = limit(&origins.get_untracked(), "reference origins")?;
            snapshot.with_value(|s| model::capture(s, pairwise_review::Options::Whole(options)))
        })();
        let capture = match captured {
            Ok(c) => c,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let active = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(Some(active.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        let bound_settings = settings.get_untracked();
        let (reference, query) =
            snapshot.with_value(|s| (s.target.clone(), s.report.products.last().unwrap().clone()));
        busy.set(true);
        proof.set(None);
        selected.set(None);
        status.set("Comparing the captured actual product with the intended reference…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let response = active
                .ask(Request {
                    id: token,
                    revision: capture.pairwise_revision.clone(),
                    job: AnalysisJob::CompareSequences {
                        reference,
                        query,
                        options: capture.options,
                    },
                })
                .await;
            active.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || settings.get_untracked() != bound_settings
            {
                return;
            }
            worker.set_value(None);
            busy.set(false);
            let result = (|| {
                guard.run(())?;
                let response = response?;
                if response.id != token || response.revision != capture.pairwise_revision {
                    return Err("Comparison response identity changed.".into());
                }
                let report = match response.result {
                    AnalysisResult::ReviewedComparison(r) => r,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected comparison response.".into()),
                };
                snapshot.with_value(|s| model::validate(s, &capture, &report))?;
                Ok::<_, String>(report)
            })();
            match result {
                Ok(report) => {
                    let index = PageIndex::new(&report.comparison.alignment);
                    status.set("Sequence differences are diagnostic; the original cloning verdict and chemical assumptions are unchanged.".into());
                    page.set(0);
                    differences_page.set(0);
                    proof.set(Some(Arc::new(Bound {
                        capture,
                        report,
                        index,
                    })));
                }
                Err(error) => status.set(error),
            }
        });
    };
    let focus = Callback::new(move |index: usize| {
        let Some(bound) = proof.get_untracked() else {
            return;
        };
        let result = guard.run(()).and_then(|_| {
            snapshot
                .with_value(|s| model::selected_difference(s, &bound.capture, &bound.report, index))
        });
        match result {
            Ok(_) => {
                selected.set(Some(index));
                if let Some(column) = bound.index.difference_column(index) {
                    page.set(column / 600);
                }
            }
            Err(error) => status.set(error),
        }
    });
    let export = move |_| {
        let Some(bound) = proof.get_untracked() else {
            return;
        };
        let result = (|| {
            guard.run(())?;
            snapshot.with_value(|s| model::validate(s, &bound.capture, &bound.report))?;
            let data = serde_json::json!({"scope":"Read-only comparison of captured intended and reconstructed DNA; no history events or experimental verification","context_revision":bound.capture.context_revision,"plan_revision":bound.capture.plan_revision,"target_revision":bound.capture.target_revision,"product_revision":bound.capture.product_revision,"comparison":bound.report});
            let bytes = serde_json::to_vec(&data).map_err(|e| e.to_string())?;
            if bytes.len() > collagene_core::pairwise_history::MAX_SAVED_BYTES {
                return Err("Comparison export exceeds the native report limit.".into());
            }
            crate::app::prepare_download(
                "cloning-target-comparison.json",
                "application/json",
                &bytes,
            )
        })();
        match result {
            Ok(artifact) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    guard.run(())?;
                    snapshot.with_value(|s| model::validate(s, &bound.capture, &bound.report))
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() {
                        return;
                    }
                    status.set(match result {
                        Ok(()) => "Diagnostic comparison saved.".into(),
                        Err(error) => error,
                    });
                });
            }
            Err(error) => status.set(error),
        }
    };
    view! {<details class="assembly-product-details"><summary>"Compare actual DNA with intended target"</summary>
        <p>{snapshot.with_value(|s|{let p=s.report.products.last().unwrap();format!("Intended: {} nt / {:?} · actual: {} nt / {:?} · cloning verdict {:?}",s.target.sequence.len(),s.target.topology,p.sequence.len(),p.topology,s.report.status)})}</p>
        <p class="helper">"R is the intended reference; Q is the actual reconstructed product. Global sequence comparison explains differences. It does not change the cloning verdict, establish chemical compatibility or verify an experiment."</p>
        <details><summary>"Comparison policy and limits"</summary>
            <label><input type="checkbox" prop:checked=move||rotate.get() on:change=move|e|rotate.set(event_target_checked(&e))/>"Allow circular reference-origin rotation"</label>
            <label><input type="checkbox" prop:checked=move||reverse.get() on:change=move|e|reverse.set(event_target_checked(&e))/>"Also inspect reverse-complement query (diagnostic only)"</label>
            <div class="primer-grid"><label>"Cells per alignment"<input type="number" min="1" prop:value=move||cells.get() on:input=move|e|cells.set(event_target_value(&e))/></label><label>"Total cells"<input type="number" min="1" prop:value=move||total.get() on:input=move|e|total.set(event_target_value(&e))/></label><label>"Reference origin candidates"<input type="number" min="1" max="64" prop:value=move||origins.get() on:input=move|e|origins.set(event_target_value(&e))/></label></div>
        </details>
        <button class="button secondary" disabled=move||busy.get()||!current.get() on:click=compute>"Compare captured sequences"</button>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Comparison cancelled.".into());}>"Cancel comparison"</button></Show>
        <p class="helper" role="status">{move||if current.get(){status.get()}else{"Workspace changed; reconstruct the product again.".into()}}</p>
        {move||proof.get().map(|bound|{
            let a=&bound.report.comparison.alignment;let c=&bound.report.comparison;let stats=&a.stats;
            let count=a.differences.len();let pages=bound.index.page_count();let dp=count.div_ceil(20).max(1);let p=page.get().min(pages-1);let d=differences_page.get().min(dp-1);
            view!{<section>
                <p>{format!("{} difference runs · certain matches {} · compatible ambiguity {} · mismatches {} · inserted/deleted query bases {}/{}",count,stats.certain_matches,stats.ambiguous_compatible,stats.mismatches,stats.inserted_query_bases,stats.deleted_reference_bases)}</p>
                <p>{format!("Reference origin {} · query reverse complement {} · candidates {} · origin search exhaustive {} · orientation search exhaustive {} · best score proven {}",c.reference_origin,c.query_reversed,c.evaluated_candidates,c.exhaustive_reference_origins,c.exhaustive_query_orientations,c.optimality_proven)}</p>
                {c.warnings.iter().take(40).map(|text|view!{<p class="helper">{text.clone()}</p>}).collect_view()}
                <p class="helper">"Alignment R/Q offsets below are zero-based in the explicitly rotated/oriented sequences. Difference navigation maps them back to the original snapshots."</p>
                <pre>{render_page(a,&bound.index,p)}</pre>
                <div class="method-row-actions"><button class="text-tool" disabled=p==0 on:click=move|_|page.set(p.saturating_sub(1))>"Previous alignment page"</button><button class="text-tool" disabled={p+1>=pages} on:click=move|_|page.set((p+1).min(pages-1))>"Next alignment page"</button></div>
                {a.differences.iter().enumerate().skip(d*20).take(20).map(|(i,delta)|view!{<article><button class="text-tool" on:click=move|_|focus.run(i)>{format!("Difference {} · {:?} · R [{}, {}), Q [{}, {})",i+1,delta.kind,delta.reference_start,delta.reference_end,delta.query_start,delta.query_end)}</button><code>{format!("R: {}\nQ: {}",clip(&delta.reference_bases),clip(&delta.query_bases))}</code></article>}).collect_view()}
                <div class="method-row-actions"><button class="text-tool" disabled=d==0 on:click=move|_|differences_page.set(d.saturating_sub(1))>"Previous differences"</button><button class="text-tool" disabled={d+1>=dp} on:click=move|_|differences_page.set((d+1).min(dp-1))>"Next differences"</button></div>
                <button class="text-tool" disabled=move||!current.get() on:click=export>"Download diagnostic comparison JSON"</button>
            </section>}
        })}
        <For each={move||proof.get().and_then(|b|selected.get().map(|i|(b,i))).into_iter().collect::<Vec<_>>()} key={|(b,i)|(b.capture.context_revision.clone(),*i)} children={move|(bound,index)|{
            snapshot.with_value(|s|match model::selected_difference(s,&bound.capture,&bound.report,index){
                Ok((reference,query))=>view!{<section><p>{format!("Original intended intervals {:?}; actual interval {:?} · zero-based half-open, empty intervals are carets",reference,query)}</p><details open><summary>"Intended reference at difference"</summary><ReferenceCanvas product=s.target.clone() identity=format!("{}:intended:{index}",bound.capture.context_revision) current on_message=Callback::new(move|m|status.set(m)) initial_spans=reference/></details><details open><summary>"Actual product at difference"</summary><ReferenceCanvas product=s.report.products.last().unwrap().clone() identity=format!("{}:actual:{index}",bound.capture.context_revision) current on_message=Callback::new(move|m|status.set(m)) initial_spans=vec![query]/></details></section>}.into_any(),
                Err(error)=>view!{<p role="alert">{error}</p>}.into_any(),
            })
        }}/>
    </details>}
}
