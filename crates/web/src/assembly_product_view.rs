//! Read-only final-step reconstruction review. The parent supplies a fresh private
//! verifier result; a serialized report is never accepted as installation proof.
use crate::{
    assembly_junctions_view::AssemblyJunctionsPanel,
    assembly_product_model::{Snapshot, operation_label},
    renderer::MapRenderer,
    sequence_view::SequenceView,
};
use collagene_core::{
    Location, Molecule,
    verification::{self, DesignPlan, Status, VerificationReport},
};
use leptos::{html, prelude::*};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::JsCast;
use web_sys::HtmlCanvasElement;

fn clip(value: &str, limit: usize) -> String {
    let mut chars = value.chars();
    let text = chars.by_ref().take(limit).collect::<String>();
    if chars.next().is_some() {
        format!("{text}…")
    } else {
        text
    }
}
/// Private display request: the parent has checked the complete native receipt.
/// A new serial can revisit the same exact interval without recreating the map.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ReferenceNavigation {
    product_revision: String,
    parts: Vec<(usize, usize)>,
    serial: u64,
}

#[component]
pub fn AssemblyProductView(
    plan: DesignPlan,
    target: Molecule,
    report: VerificationReport,
    epoch: RwSignal<u64>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    #[prop(optional_no_strip)] guard: Option<Callback<(), Result<(), String>>>,
    #[prop(optional_no_strip)] on_finalize: Option<
        Callback<collagene_core::transfer::Candidate, Result<(), String>>,
    >,
) -> impl IntoView {
    let snapshot = StoredValue::new(Arc::new(Snapshot::new(plan, target, report)));
    let captured_epoch = epoch.get_untracked();
    let current = Memo::new(move |_| epoch.get() == captured_epoch);
    let status = RwSignal::new(String::new());
    let junctions_open = RwSignal::new(false);
    let navigation = RwSignal::new(None::<ReferenceNavigation>);
    let navigation_serial = RwSignal::new(0u64);
    let navigation_request = Memo::new(move |_| navigation.get());

    let intermediates_open = RwSignal::new(false);
    let inspected_step =
        RwSignal::new(snapshot.with_value(|s| s.prefix_available.saturating_sub(1)));
    let available = snapshot.with_value(|s| s.product_revision.is_some());
    let report_available = snapshot.with_value(|s| s.prefix_error.is_none());
    let output_available =
        snapshot.with_value(|s| s.product_revision.is_some() && s.report.status != Status::Invalid);
    let identity = Memo::new(move |_| {
        if current.get() {
            snapshot.with_value(|s| (s.plan_revision.clone(), s.target_revision.clone()))
        } else {
            (String::new(), String::new())
        }
    });
    let capture = Callback::new(move |()| {
        if let Some(guard) = guard {
            guard.run(())?;
        }
        snapshot.with_value(|s| {
            if !current.get_untracked() {
                return Err("The source context changed; reconstruct the plan again.".into());
            }
            collagene_core::assembly_junctions::source_preflight(&s.plan, &s.target)
                .map_err(|e| e.to_string())?;
            if s.report.engine_version != verification::ENGINE_VERSION
                || s.report.plan_revision != verification::plan_revision(&s.plan)
                || s.report.target_revision != s.target.revision()
                || s.plan.target_revision != s.report.target_revision
            {
                return Err("The reconstruction receipt differs from this plan or target.".into());
            }
            // Valid failed-final receipts can still export diagnostic protocols and
            // inspect the failed step. Opening/finalization requires an actual final
            // product separately; a successful preparation never stands in for it.
            Ok((s.plan.clone(), s.target.clone()))
        })
    });
    let navigate_junction = Callback::new(
        move |(report, index, kind): (
            Arc<collagene_core::assembly_junctions::Report>,
            usize,
            crate::assembly_junction_navigation_model::Kind,
        )|
              -> Result<(), String> {
            if let Some(guard) = guard {
                guard.run(())?;
            }
            if !current.get_untracked() {
                return Err("Reconstruct the current product before navigating junctions.".into());
            }
            let selected = snapshot.with_value(|s| {
                use crate::assembly_junction_navigation_model as model;
                let bound = model::capture(s, report.options.clone())?;
                let selected = model::select(s, &bound, &report, index, kind)?;
                model::validate_selection(s, &bound, &report, &selected)?;
                Ok::<_, String>(selected)
            })?;
            navigation_serial.update(|serial| *serial = serial.wrapping_add(1));
            navigation.set(Some(ReferenceNavigation {
                product_revision: selected.product_revision,
                parts: selected.parts,
                serial: navigation_serial.get_untracked(),
            }));
            status.set("Exact reconstructed-product bases selected in Sequence. Plasmid rotation and zoom are preserved.".into());
            Ok(())
        },
    );
    let export_identity = Memo::new(move |_| {
        if current.get() {
            snapshot.with_value(|s| {
                s.product_revision
                    .as_ref()
                    .map(|p| format!("{}:{}:{}", s.plan_revision, s.target_revision, p))
                    .unwrap_or_default()
            })
        } else {
            String::new()
        }
    });
    let export_actual = Callback::new(move |()| {
        capture.run(())?;
        snapshot.with_value(|s| {
            if s.product_revision.is_none() {
                return Err("No actual final product is available for export.".into());
            }
            s.report
                .products
                .last()
                .cloned()
                .ok_or_else(|| "The actual final product is unavailable.".into())
        })
    });
    let export_report = move |_| {
        let result = capture.run(()).and_then(|_| {
            snapshot.with_value(|s| {
                if s.prefix_error.is_some() {
                    return Err("The report context is not bound to this plan and target.".into());
                }
                let bytes = serde_json::to_vec_pretty(&s.report).map_err(|e| e.to_string())?;
                if bytes.len() > 64 * 1024 * 1024 {
                    return Err("The diagnostic JSON exceeds the 64MiB export profile.".into());
                }
                crate::app::prepare_download(
                    "cloning-reconstruction.json",
                    "application/json",
                    &bytes,
                )
            })
        });
        match result {
            Ok(artifact) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    capture.run(()).map(|_| ())
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() { return; }
                    status.set(match result { Ok(()) => "Complete reconstruction report saved, including retained intermediate and graph-only evidence. This is diagnostic evidence, not material authority.".into(), Err(error) => error });
                });
            }
            Err(error) => status.set(error),
        }
    };
    let open = move |_| {
        if !current.get_untracked() {
            status.set(
                "The workspace changed; reconstruct again before opening this product reference."
                    .into(),
            );
            return;
        }
        let result = capture.run(()).and_then(|_| {
            snapshot.with_value(|s| {
                if s.product_revision.is_none() || s.error.is_some() {
                    return Err("No bound final product is available.".into());
                }
                let final_product = s.report.products.last().ok_or("No final product")?;
                on_products.run(vec![final_product.clone()])
            })
        });
        match result{Ok(())=>status.set("Reconstructed reference opened. The original sources and intended target remain separate.".into()),Err(error)=>status.set(error)}
    };
    view! {<section class="assembly-product-view" aria-label="Reconstructed assembly product">
        <header class="assembly-product-header product-review-header"><div class="product-review-title"><h2>"Reconstructed product"</h2>
            <p class="product-review-result"><strong>{snapshot.with_value(|s|format!("{:?}",s.report.status))}</strong><span>{move||if current.get(){"Current reconstruction snapshot"}else{"Historical snapshot · reconstruct again"}}</span></p>
        </div><div class="product-review-actions">
            {on_finalize.map(|finalize|view!{<button class="button primary" disabled=move||!output_available||!current.get() on:click=move|_|{let result=capture.run(()).and_then(|(plan,_)|finalize.run(collagene_core::transfer::Candidate::new(plan)));if let Err(error)=result{status.set(error);}}>"Create assembly outputs"</button>})}
            <details class="product-export-menu" data-transient-menu=""><summary>"Export primers / steps"</summary><crate::assembly_outputs::AssemblyOutputs capture identity epoch/><Show when=move||available><crate::assembly_product_exports::AssemblyProductExports capture=export_actual identity=export_identity epoch/></Show><button class="text-tool" disabled=move||!report_available||!current.get() on:click=export_report>"Export full reconstruction JSON"</button></details>
        </div></header>
        {snapshot.with_value(|s|(s.product_revision.is_some()&&s.report.status==Status::Invalid).then(||view!{<p class="assembly-product-status" role="status">"Output creation is unavailable because the independent target check is invalid. Inspect the actual DNA and findings, then correct the design."</p>}))}
        <p class="helper product-review-scope">"Actual reconstructed reference DNA, separate from intended DNA. This view does not establish experimental success, exhaustive product uniqueness or newly prepared material."</p>
        <p class="assembly-product-status" role="status">{move||if !current.get(){"Workspace context changed. This snapshot is historical; reconstruct again before opening or inspecting it.".to_string()}else{status.get()}}</p>
        {snapshot.with_value(|s|s.error.clone()).map(|error|view!{<p class="assembly-product-status" role="alert">{error}</p>})}
        {snapshot.with_value(|s|if s.product_revision.is_some(){s.report.products.last().map(|product|{
            let product=product.clone();let product_identity=format!("{}:{}",s.plan_revision,s.product_revision.as_deref().unwrap_or_default());
            view!{<div>
                <p>{format!("Final step {} · {} bases · {:?} · {} annotations",s.plan.steps.len(),product.sequence.len(),product.topology,product.features.len())}</p>
                {(s.report.status==Status::Invalid).then(||view!{<p class="assembly-product-status" role="alert">"The plan/independent target was not accepted. The actual final reconstructed DNA is shown for diagnosis; it must not be mistaken for the intended target."</p>})}
                {(s.report.status==Status::Incomplete).then(||view!{<p class="helper">"Reconstruction is incomplete under the declared model. This reference view does not resolve its missing evidence or assumptions."</p>})}
                <ReferenceCanvas product=product identity=product_identity current=current navigation=Some(navigation_request) on_message=Callback::new(move|text|status.set(text))/>
                <button class="button secondary" disabled=move||!current.get() on:click=open>"Open reconstructed product reference"</button>
            </div>}
        })}else{None})}
        <Show when=move||available><crate::assembly_comparison_view::AssemblyComparison snapshot=snapshot.get_value() guard=Callback::new(move|()|capture.run(()).map(|_|())) current epoch/></Show>
        <details class="assembly-product-details" on:toggle=move|event:web_sys::Event|{
            let open=event.current_target().and_then(|target|target.dyn_into::<web_sys::HtmlElement>().ok()).is_some_and(|element|element.has_attribute("open"));
            if open {
                if let Err(error)=capture.run(()){status.set(error);intermediates_open.set(false);return;}
            }
            intermediates_open.set(open);
        }><summary>"Inspect preparation steps and intermediate DNA"</summary>
            {snapshot.with_value(|s|s.prefix_error.clone()).map(|error|view!{<p role="alert">{error}</p>})}
            <p class="helper">{snapshot.with_value(|s|format!("{} of {} step products available in this reconstruction. Intermediate views do not select the final output or restore physical material.",s.prefix_available,s.plan.steps.len()))}</p>
            {snapshot.with_value(|s|{
                let index=s.prefix_available;
                if s.prefix_error.is_some() || index>=s.plan.steps.len(){return None;}
                let graphs=s.report.step_reports.iter().filter(|r|r.step_index==index).count();
                let processed=s.report.processed_step_reports.iter().filter(|r|r.step_index==index).count();
                let physical=s.report.physical_processing_reports.iter().filter(|r|r.step_index==index).count();
                let reference=s.report.reference_step_reports.iter().filter(|r|r.step_index==index).count();
                let annotation=s.report.annotation_reports.iter().filter(|r|r.step_index==index).count();
                Some(view!{<p>{format!("Next step not reconstructed: {}. {}",index+1,operation_label(&s.plan.steps[index]))}</p><p class="helper">{format!("Retained typed evidence at this step: {graphs} graph, {processed} processed graph, {physical} physical processing, {reference} reference and {annotation} annotation reports. Export the full reconstruction JSON for complete evidence, including graph-only results with no DNA product.")}</p>})
            })}
            <Show when=move||snapshot.with_value(|s|s.prefix_available>0)>
                <label>"Reconstructed step"<select aria-label="Reconstructed preparation step" disabled=move||!current.get() prop:value=move||inspected_step.get().to_string() on:change=move|event|{
                    if let Ok(index)=event_target_value(&event).parse::<usize>() {
                        if let Err(error)=capture.run(()){status.set(error);intermediates_open.set(false);return;}
                    if snapshot.with_value(|s|s.product_at(index).is_ok()) {inspected_step.set(index);}
                    }
                }>{snapshot.with_value(|s|(0..s.prefix_available).map(|index|view!{<option value=index.to_string()>{format!("{}. {}{}",index+1,operation_label(&s.plan.steps[index]),if index+1==s.plan.steps.len(){" · final step"}else{" · intermediate"})}</option>}).collect_view())}</select></label>
                <div class="method-row-actions"><button class="text-tool" disabled=move||!current.get()||inspected_step.get()==0 on:click=move|_|{match capture.run(()){Ok(_)=>inspected_step.update(|index|*index=index.saturating_sub(1)),Err(error)=>{status.set(error);intermediates_open.set(false);}}}>"Previous step"</button><button class="text-tool" disabled={move||!current.get()||snapshot.with_value(|s|inspected_step.get().saturating_add(1)>=s.prefix_available)} on:click=move|_|{match capture.run(()){Ok(_)=>inspected_step.update(|index|*index=index.saturating_add(1)),Err(error)=>{status.set(error);intermediates_open.set(false);}}}>"Next step"</button></div>
                <For each=move||if intermediates_open.get()&&current.get(){vec![inspected_step.get()]}else{Vec::new()} key=|index|*index children=move|index|{
                    snapshot.with_value(|s|match s.product_at(index) {
                        Ok(product)=>{
                            let product=product.clone();
                            let identity=format!("{}:step{}:{}",s.plan_revision,index,product.revision());
                            let physical=s.descriptive_physical_at(index);
                            view!{<section aria-label="Read-only intermediate reconstruction">
                                <p><strong>{format!("Step {} · {}",index+1,operation_label(&s.plan.steps[index]))}</strong>{format!(" · {} bases · {:?} · {} annotations",product.sequence.len(),product.topology,product.features.len())}</p>
                                <ReferenceCanvas product identity current on_message=Callback::new(move|text|status.set(text))/>
                                {physical.map(|summary|view!{<details><summary>"Descriptive strand/end receipt"</summary>{crate::end_view::summary(&summary)}<p class="helper">"Read-only evidence from this reconstructed step. Detached metadata cannot establish material for another operation."</p></details>})}
                            </section>}.into_any()
                        },
                        Err(error)=>view!{<p role="alert">{error}</p>}.into_any(),
                    })
                }/>
            </Show>
        </details>
        {snapshot.with_value(|s|s.report.physical_processing_reports.iter().filter_map(|step|match &step.evidence {
            collagene_core::verification::PhysicalProcessingEvidence::HomologyAssembly(report)=>Some(homology_processing_view(step.step_index,report)),_=>None,
        }).collect_view())}
        <details class="assembly-product-details"><summary>"Sources, preparation history and verification findings"</summary>
            <p>{snapshot.with_value(|s|format!("Scientific engine: {}",s.report.engine_version))}</p>
            <p>{snapshot.with_value(|s|format!("Plan {} · independent target {} ({})",s.plan_revision,clip(&s.target.name,120),s.target_revision))}</p>
            <ol>{snapshot.with_value(|s|s.plan.inputs.iter().take(128).map(|input|view!{<li>{format!("{} · {} · {} bases · revision {}",clip(&input.id,128),clip(&input.molecule.name,120),input.molecule.sequence.len(),input.revision)}</li>}).collect_view())}</ol>
            <ol>{snapshot.with_value(|s|s.plan.steps.iter().take(100).map(|step|view!{<li>{operation_label(step)}</li>}).collect_view())}</ol>
            {snapshot.with_value(|s|(s.plan.steps.len()>100).then(||view!{<p class="helper">{format!("First100 of {} typed steps shown; the saved candidate/protocol retains the full plan.",s.plan.steps.len())}</p>}))}
            <ul>{snapshot.with_value(|s|s.report.findings.iter().take(100).map(|finding|view!{<li>{clip(finding,2048)}</li>}).collect_view())}</ul>
            <p class="helper">"Finding previews are bounded. Keep the full reconstruction JSON and HTML/PDF protocol for all findings and assumptions. Source annotation claims are not functional validation."</p>
        </details>
        {snapshot.with_value(|s|s.physical.as_ref().map(|physical|view!{<details class="assembly-product-details"><summary>{format!("Fresh report physical receipt · {} retained nicks",physical.nicks.len())}</summary>
            <p>{format!("Top {} nt · bottom {} nt · covalently closed in this model: {}",physical.top_length,physical.bottom_length,physical.covalently_closed_in_model)}</p>
            {physical.ends.as_ref().map(|(left,right)|view!{<p>{format!("Left {:?} / 5′ {:?} / 3′ {:?}; right {:?} / 5′ {:?} / 3′ {:?}",left.polarity,left.chemistry.five_prime,left.chemistry.three_prime,right.polarity,right.chemistry.five_prime,right.chemistry.three_prime)}</p>})}
            <ul>{physical.nicks.iter().take(100).map(|nick|view!{<li>{format!("{:?} strand at top-reference boundary {} (zero based)",nick.strand,nick.position)}</li>}).collect_view()}</ul>
            {(physical.nicks.len()>100).then(||view!{<p>"First100 nicks shown; complete state remains in the fresh report."</p>})}
            <p class="helper">"Descriptive output from this fresh reconstruction only. This receipt cannot restore physical authority when detached or imported."</p>
        </details>}))}
        <button class="text-tool" disabled=move||!available||!current.get() on:click=move|_|junctions_open.update(|open|*open=!*open)>{move||if junctions_open.get(){"Close fresh junction inspection"}else{"Inspect reconstructed junctions"}}</button>
        <Show when=move||junctions_open.get()&&available&&current.get()><AssemblyJunctionsPanel capture=capture identity=identity epoch=epoch on_navigate=Some(navigate_junction)/></Show>
    </section>}
}

fn feature_parts(location: &Location, length: usize) -> Result<Vec<(usize, usize)>, String> {
    let mut stack = vec![(location, 0usize)];
    let mut parts = Vec::new();
    let mut work = 0usize;
    while let Some((location, depth)) = stack.pop() {
        work += 1;
        if work > 4096 || depth >= 32 {
            return Err("Annotation navigation exceeds the bounded part limit; use the full annotation evidence.".into());
        }
        match location {
            Location::Range { start, end, .. } => {
                if *start < *end && *end <= length {
                    parts.push((*start, *end));
                }
            }
            Location::Complement(inner) => stack.push((inner, depth + 1)),
            Location::Join(children) | Location::Order(children) => {
                if children.len() + stack.len() + work > 4096 {
                    return Err("Annotation has too many parts for interactive navigation.".into());
                }
                stack.extend(children.iter().rev().map(|child| (child, depth + 1)));
            }
            Location::Unresolved(_) => {
                return Err(
                    "Unresolved annotation coordinates remain evidence; no position was invented."
                        .into(),
                );
            }
        }
    }
    Ok(parts)
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Sequence,
    Map,
    Split,
}
#[component]
pub(crate) fn ReferenceCanvas(
    product: Molecule,
    identity: String,
    current: Memo<bool>,
    on_message: Callback<String>,
    #[prop(optional)] initial_spans: Vec<(usize, usize)>,
    #[prop(optional_no_strip)] navigation: Option<Memo<Option<ReferenceNavigation>>>,
    #[prop(optional_no_strip)] pcr_overlay: Option<Memo<Option<Arc<crate::pcr_overlay::Overlay>>>>,
) -> impl IntoView {
    let visible_pcr_overlay = Memo::new(move |_| {
        if current.get() {
            pcr_overlay.and_then(|overlay| overlay.get())
        } else {
            None
        }
    });
    let source = StoredValue::new(Arc::new(product));
    // One memo-owned model for the existing virtual viewer. Source access in
    // this component stays borrowed; selection changes never clone full DNA.
    let molecule = Memo::new(move |_| source.with_value(|source| source.as_ref().clone()));
    let revision = StoredValue::new(source.with_value(|source| source.revision()));
    let identity = StoredValue::new(identity);
    let selection = RwSignal::new(if initial_spans.len() == 1 {
        initial_spans.first().copied()
    } else {
        None
    });
    let selected_feature = RwSignal::new(None::<usize>);
    let regions = RwSignal::new(if initial_spans.len() > 1 {
        Some((revision.get_value(), initial_spans))
    } else {
        None
    });
    let features_open = RwSignal::new(false);
    let feature_page = RwSignal::new(0usize);
    let mode = RwSignal::new(Mode::Split);
    let map_stats = RwSignal::new(crate::map_labels::LabelStats::default());
    let sequence_pane = NodeRef::<html::Div>::new();
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || cleanup_alive.store(false, Ordering::Relaxed));
    Effect::new(move |_| {
        let Some(request) = navigation.and_then(|value| value.get()) else {
            return;
        };
        if !current.get_untracked() {
            return;
        }
        let valid =
            source.with_value(|source| {
                let n = source.sequence.len();
                request.product_revision == revision.get_value()
                    && !request.parts.is_empty()
                    && request.parts.len() <= 1024
                    && request.parts.iter().all(|&(a, b)| a < b && b <= n)
                    && request.parts.iter().enumerate().all(|(i, &(a, b))| {
                        request.parts[..i].iter().all(|&(c, d)| !(a < d && c < b))
                    })
            });
        if !valid {
            on_message
                .run("The navigation interval differs from this exact product reference.".into());
            return;
        }
        batch(|| {
            selected_feature.set(None);
            selection.set(request.parts.first().copied());
            regions.set(Some((revision.get_value(), request.parts.clone())));
            if mode.get_untracked() == Mode::Map {
                mode.set(Mode::Split);
            }
        });
        // Wait for Map→Split layout before scrolling the existing sequence pane.
        // No molecule, renderer, camera or scientific computation is replaced.
        let alive = alive.clone();
        let callback = wasm_bindgen::closure::Closure::once_into_js(move || {
            if !alive.load(Ordering::Relaxed) || !current.get_untracked() {
                return;
            }
            if navigation.and_then(|n| n.get_untracked()).as_ref() != Some(&request) {
                return;
            }
            // A manual feature, Go to, Clear, base selection or Map action is
            // newer than this queued frame even when the request nonce stayed fixed.
            if mode.get_untracked() == Mode::Map
                || selected_feature.get_untracked().is_some()
                || selection.get_untracked() != request.parts.first().copied()
                || !regions.with_untracked(|value| {
                    value.as_ref().is_some_and(|(bound, parts)| {
                        bound == &request.product_revision && parts == &request.parts
                    })
                })
            {
                return;
            }
            if let Some(pane) = sequence_pane.get() {
                pane.scroll_into_view();
            }
            // Explicit notification also revisits unchanged ranges after manual scrolling.
            selection.set(request.parts.first().copied());
        });
        if let Some(window) = web_sys::window() {
            let _ = window.request_animation_frame(callback.unchecked_ref());
        }
    });
    let start = RwSignal::new(String::new());
    let end = RwSignal::new(String::new());
    let highlighted = Memo::new(move |_| {
        if let Some(index) = selected_feature.get() {
            return source.with_value(|source| {
                source
                    .features
                    .get(index)
                    .and_then(|f| feature_parts(&f.location, source.sequence.len()).ok())
                    .unwrap_or_default()
            });
        }
        if let Some((bound, spans)) = regions.get() {
            if bound == revision.get_value() {
                return spans;
            }
        }
        selection
            .get()
            .filter(|(a, b)| a < b)
            .into_iter()
            .collect::<Vec<_>>()
    });
    let show_feature = Callback::new(move |index: usize| {
        let parts = source.with_value(|source| {
            source
                .features
                .get(index)
                .ok_or("Annotation no longer exists".to_string())
                .and_then(|f| feature_parts(&f.location, source.sequence.len()))
        });
        match parts {
            Ok(parts) => {
                selected_feature.set(Some(index));
                selection.set(parts.first().copied());
                regions.set(Some((revision.get_value(), parts)));
                features_open.set(true);
            }
            Err(error) => {
                selected_feature.set(Some(index));
                selection.set(None);
                regions.set(None);
                features_open.set(true);
                on_message.run(error);
            }
        }
    });
    let canvas = NodeRef::<html::Canvas>::new();
    let labels = NodeRef::<html::Canvas>::new();
    let renderer = StoredValue::new_local(None::<MapRenderer>);
    let map_error = RwSignal::new(String::new());
    Effect::new(move |_| {
        let viewing = mode.get() != Mode::Sequence;
        let context_current = current.get();
        let parts = highlighted.get();
        let feature = selected_feature.get();
        if !context_current {
            if let Some(map) = renderer.get_value() {
                map.set_pcr_overlay(None);
                map.pointer_cancel();
            }
            return;
        }
        if !viewing {
            if let Some(map) = renderer.get_value() {
                map.pointer_cancel();
            }
            return;
        }
        let (Some(canvas), Some(labels)) = (canvas.get(), labels.get()) else {
            return;
        };
        if renderer.get_value().is_none() {
            let created = source.with_value(|source| {
                MapRenderer::new(
                    canvas.unchecked_ref::<HtmlCanvasElement>().clone(),
                    labels.unchecked_ref::<HtmlCanvasElement>().clone(),
                    source.as_ref().clone(),
                    Rc::new(move |stats| map_stats.set(stats)),
                )
            });
            match created {
                Ok(map) => {
                    map.set_interaction_context(&identity.get_value());
                    map.set_selection_mode(false);
                    renderer.set_value(Some(map));
                }
                Err(error) => {
                    map_error.set(error);
                    return;
                }
            }
        }
        if let Some(map) = renderer.get_value() {
            source.with_value(|source| {
                map.set_document_spans(source, &revision.get_value(), &parts, feature)
            });
            map.set_pcr_overlay(visible_pcr_overlay.get());
            map.queue_draw();
        }
    });
    on_cleanup(move || {
        if let Some(map) = renderer.get_value() {
            map.pointer_cancel();
        }
        renderer.set_value(None);
    });
    let navigate = move |_| {
        let parsed = (|| -> Result<(usize, usize), String> {
            let a = start
                .get_untracked()
                .parse::<usize>()
                .map_err(|_| "Start must be a one-based coordinate")?
                .checked_sub(1)
                .ok_or("Start begins at1")?;
            let b = end
                .get_untracked()
                .parse::<usize>()
                .map_err(|_| "End must be a one-based inclusive coordinate")?;
            if a >= b || b > source.with_value(|source| source.sequence.len()) {
                return Err("Choose a nonempty interval within the displayed reference.".into());
            }
            Ok((a, b))
        })();
        match parsed {
            Ok((a, b)) => {
                selected_feature.set(None);
                regions.set(None);
                selection.set(Some((a, b)));
                if mode.get_untracked() == Mode::Map {
                    mode.set(Mode::Split);
                }
            }
            Err(error) => on_message.run(error),
        }
    };
    let view_identity = Memo::new(move |_| identity.get_value());
    let locked = Memo::new(move |_| !current.get());
    let no_restrictions =
        Memo::new(
            move |_| None::<Arc<Vec<collagene_core::restriction_inspection::InspectionSite>>>,
        );
    let restriction_visible = RwSignal::new(false);
    view! {<div class="assembly-product-stage product-reference-shell">
        <div class="assembly-product-controls product-reference-toolbar" aria-label="Read-only reference view">
            <div class="product-view-tabs" role="group" aria-label="Product view">
                <button class="text-tool" aria-pressed=move||(mode.get()==Mode::Sequence).to_string() on:click=move|_|mode.set(Mode::Sequence)>"Sequence"</button>
                <button class="text-tool" aria-pressed=move||(mode.get()==Mode::Map).to_string() on:click=move|_|mode.set(Mode::Map)>"Map"</button>
                <button class="text-tool" aria-pressed=move||(mode.get()==Mode::Split).to_string() on:click=move|_|mode.set(Mode::Split)>"Split"</button>
            </div>
            <Show when=move||mode.get()!=Mode::Sequence><div class="product-map-zoom" role="group" aria-label="Product map zoom">
                <button class="text-tool" aria-label="Zoom product map out" title="Zoom out" disabled=move||!current.get() on:click=move|_|if let Some(map)=renderer.get_value(){map.zoom(1./1.15);}>"−"</button>
                <button class="text-tool" title="Fit product map" disabled=move||!current.get() on:click=move|_|if let Some(map)=renderer.get_value(){map.reset();}>"Fit"</button>
                <button class="text-tool" aria-label="Zoom product map in" title="Zoom in" disabled=move||!current.get() on:click=move|_|if let Some(map)=renderer.get_value(){map.zoom(1.15);}>"+"</button>
            </div></Show>
            <details class="product-interval-menu" data-transient-menu=""><summary>"Go to interval"</summary><div class="product-interval-fields">
                <label>"Start (1-based)"<input type="number" min="1" aria-label="Reference start" disabled=move||!current.get() prop:value=move||start.get() on:input=move|event|start.set(event_target_value(&event))/></label>
                <label>"End (inclusive)"<input type="number" min="1" aria-label="Reference end" disabled=move||!current.get() prop:value=move||end.get() on:input=move|event|end.set(event_target_value(&event))/></label>
                <button class="button secondary" disabled=move||!current.get() on:click=navigate>"Show exact interval"</button>
            </div></details>
            <button class="text-tool" aria-expanded=move||features_open.get().to_string() on:click=move|_|features_open.update(|value|*value=!*value)>"Annotations"</button>
            <button class="text-tool" disabled=move||highlighted.with(Vec::is_empty) on:click=move|_|{selected_feature.set(None);regions.set(None);selection.set(None);}>"Clear highlight"</button>
            <span class="product-selection-summary">{move||{let spans=highlighted.get();if let Some(&(a,b))=spans.first().filter(|&&(a,b)|spans.len()==1&&a<b){format!("{}–{} · {} bases",a+1,b,b-a)}else if !spans.is_empty(){format!("{} annotation parts highlighted",spans.len())}else{"Read-only reference".into()}}}</span>
        </div>
        <p class="helper product-view-help">"Select exact bases in Sequence or Go to interval. Drag a circular map to rotate; click a feature to inspect its annotated coordinates. Editing requires a separate document copy."</p>
        <div class="product-reference-panes" data-view=move||match mode.get(){Mode::Sequence=>"sequence",Mode::Map=>"map",Mode::Split=>"split"}>
        <div class="assembly-product-map" hidden=move||mode.get()==Mode::Sequence>
            <p role="status">{move||map_error.get()}</p>
            <div class="map-stage" style="min-height:260px;position:relative">
                <canvas node_ref=canvas aria-label="Read-only reference map; drag a circular map to rotate, click annotations to inspect" style="width:100%;height:100%;display:block"
                    on:pointerdown=move|event:web_sys::PointerEvent|{if !current.get_untracked()||event.button()!=0{return;}if let Some(canvas)=canvas.get(){let _=canvas.set_pointer_capture(event.pointer_id());}if let Some(map)=renderer.get_value(){map.pointer_down(event.client_x()as f64,event.client_y()as f64);}}
                    on:pointermove=move|event:web_sys::PointerEvent|if current.get_untracked(){if let Some(map)=renderer.get_value(){map.pointer_move(event.client_x()as f64,event.client_y()as f64);}}
                    on:pointerup=move|event:web_sys::PointerEvent|{if let Some(map)=renderer.get_value(){if current.get_untracked(){if let Some(hit)=map.pointer_up(event.client_x()as f64,event.client_y()as f64){if let Some(feature)=hit.feature_index{show_feature.run(feature);}}}else{map.pointer_cancel();}}if let Some(canvas)=canvas.get(){let _=canvas.release_pointer_capture(event.pointer_id());}}
                    on:pointercancel=move|_|if let Some(map)=renderer.get_value(){map.pointer_cancel();}
                    on:lostpointercapture=move|_|if let Some(map)=renderer.get_value(){map.pointer_cancel();}
                    on:wheel=move|event:web_sys::WheelEvent|{event.prevent_default();if !current.get_untracked(){return;}if let Some(map)=renderer.get_value(){map.zoom(if event.delta_y()>0.{0.95}else{1.05});}}/>
                <canvas node_ref=labels class="map-labels" aria-hidden="true" style="position:absolute;inset:0;width:100%;height:100%;pointer-events:none"/>
                <p class="map-geometry-note" role="status">"Overview simplified; browse annotations for omitted geometry."</p>
            </div>
        </div>
        <div node_ref=sequence_pane class="assembly-product-sequence" hidden=move||mode.get()==Mode::Map>
            <SequenceView pcr_overlay=Some(visible_pcr_overlay) molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted search_regions=regions
                on_show_features=Callback::new(move|()|features_open.set(true)) identity=view_identity locked=locked read_only=true shared_restrictions=false
                on_edit=Callback::new(|_|Err("This reference snapshot is read-only; fork or open a separate document copy to edit.".into()))
                on_history=Callback::new(|_|Err("Source history is unchanged by this read-only review.".into())) on_message=on_message
                restriction_sites=no_restrictions restriction_visible=restriction_visible
                on_show_restrictions=Callback::new(move|()|on_message.run("Restriction choices belong to the source/assembly preparation; this snapshot does not run an independent restriction scan.".into()))/>
        </div>
        </div>
        <Show when=move||map_stats.with(|s|s.omitted>0||s.abbreviated>0||s.lane_overflow>0)><p class="helper product-map-status" role="status">{move||map_stats.with(|s|{
            let mut notes=vec![format!("{} of {} map labels shown",s.displayed,s.total)];
            if s.abbreviated>0{notes.push(format!("{} abbreviated",s.abbreviated));}
            if s.lane_overflow>0{notes.push(format!("{} overlapping features hidden from the 3 map lanes",s.lane_overflow));}
            if s.unplaced>0{notes.push(format!("{} without local coordinates",s.unplaced));}
            notes.push("Browse Annotations for the complete retained features".into());notes.join(" · ")
        })}</p></Show>
        <Show when=move||features_open.get()><section class="assembly-product-features" aria-label="Reference annotations">
            <h3>"Reference annotations"</h3>
            <p>{source.with_value(|source|format!("{} annotations retained; projected or uncertain features may require review.",source.features.len()))}</p>
            <button class="text-tool" disabled=move||feature_page.get()==0 on:click=move|_|feature_page.update(|page|*page=page.saturating_sub(1))>"Previous annotations"</button>
            <button class="text-tool" disabled=move||source.with_value(|source|(feature_page.get()+1)*50>=source.features.len()) on:click=move|_|feature_page.update(|page|*page+=1)>"Next annotations"</button>
            <ul>{move||source.with_value(|source|source.features.iter().enumerate().skip(feature_page.get()*50).take(50).map(|(index,feature)|view!{<li><button class="text-tool" aria-pressed=move||(selected_feature.get()==Some(index)).to_string() on:click=move|_|show_feature.run(index)>{format!("{} · {}{}",clip(feature.label(),160),clip(&feature.kind,48),if feature.needs_review{" · review required"}else{""})}</button></li>}).collect_view())}</ul>
            {move||selected_feature.get().and_then(|index|source.with_value(|source|source.features.get(index).map(|feature|view!{<article class="annotation-details"><h4>{clip(feature.label(),256)}</h4><p>{format!("Reference feature index {} · {} · review {}",index,clip(&feature.kind,80),feature.needs_review)}</p><dl>{feature.qualifiers.iter().take(30).map(|(key,value)|view!{<dt>{clip(key,80)}</dt><dd>{value.as_deref().map(|value|clip(value,512)).unwrap_or_else(||"present".into())}</dd>}).collect_view()}</dl><p class="helper">"Bounded qualifier preview; the opened native reference retains complete annotations and source evidence. Partial/fuzzy/ordered drawing bounds do not establish a complete feature sequence; unresolved locations are not placed."</p></article>})))}
        </section></Show>
    </div>}.into_any()
}

/// Read-only processing evidence from a freshly bound complete verifier report.
/// Includes withheld outcomes; it never promotes a reference into material.
pub(crate) fn homology_processing_view(
    step_index: usize,
    report: &collagene_core::homology_assembly::Report,
) -> AnyView {
    use collagene_core::chemistry::{FivePrime, Strand, ThreePrime};
    let five = |v: FivePrime| match v {
        FivePrime::Phosphate => "phosphate",
        FivePrime::Hydroxyl => "hydroxyl",
        FivePrime::Unknown => "unknown",
    };
    let three = |v: ThreePrime| match v {
        ThreePrime::Hydroxyl => "hydroxyl",
        ThreePrime::Blocked => "blocked",
        ThreePrime::Unknown => "unknown",
    };
    let physical = report.product.as_ref().and_then(|p| p.summary().ok());
    let reusable = report.product.as_ref().is_some_and(|p| p.nicks.is_empty())
        && report
            .junctions
            .iter()
            .all(|j| j.top_sealed && j.bottom_sealed);
    let unsealed = report
        .junctions
        .iter()
        .any(|j| !j.top_sealed || !j.bottom_sealed);
    let complete = report.terminal_search_complete;
    let incomplete = report.incomplete;
    view!{<details class="assembly-product-details"><summary>{format!("Step {} · homology-assembly strand processing",step_index+1)}</summary>
        <p>{if reusable{"A duplex was reconstructed with every joining nick sealed under the declared processing assumptions."}else if physical.is_some(){"A paired diagnostic duplex is available; reusable physical material is withheld."}else{"A reusable duplex was withheld. An intended reference does not replace missing material."}}</p>
        {unsealed.then(||view!{<p role="alert">"An unsealed joining nick retains blocked or unknown terminal chemistry that the reusable material model cannot faithfully represent. Its coordinates and groups remain diagnostic evidence; this reference cannot be reused as saved physical material."</p>})}
        <p>{format!("Terminal search: {} · unresolved assumptions or alternatives: {} · {} comparisons",if complete{"complete within this adjacent-pair model"}else{"incomplete"},if incomplete{"retained"}else{"none reported within this model"},report.comparisons)}</p>
        <p>{format!("Product top-reference interval in the assembled contour: {}..{} of {} nt (zero-based, half-open). Joining coordinates below use the actual product top-reference origin.",report.reference_top_start,report.reference_top_end,report.reference_length)}</p>
        <p class="helper">"This generic model assumes declared 5′ resection, annealing, template-supported gap fill and chemistry-dependent nick sealing. It does not establish enzyme activity, vendor chemistry, isolation or a unique reaction outcome."</p>
        <details><summary>{format!("Removed 5′ source segments · {}",report.resected.len())}</summary>
            {report.resected.is_empty().then(||view!{<p>"No 5′ source segments removed."</p>})}
            <table><thead><tr><th>"Fragment"</th><th>"Strand"</th><th>"Original strand interval · zero-based half-open"</th><th>"Removed bases"</th></tr></thead><tbody>{report.resected.iter().map(|cut|view!{<tr><td>{cut.source+1}</td><td>{match cut.strand{Strand::Top=>"top",Strand::Bottom=>"bottom"}}</td><td>{format!("{}..{}",cut.start,cut.end)}</td><td>{cut.end.saturating_sub(cut.start)}</td></tr>}).collect_view()}</tbody></table>
        </details>
        {report.junctions.iter().enumerate().map(|(index,j)|view!{<details><summary>{format!("Junction {} · fragment {} → {} · {} bp homology",index+1,j.left_source+1,j.right_source+1,j.overlap)}</summary>
            <pre>{j.homology.clone()}</pre>
            <p>{format!("Homology at product top-reference interval {}..{} (zero-based, unrolled across circle closure)",j.homology_start,j.homology_end)}</p>
            <p>{format!("Joining gaps before processing: top {} nt · bottom {} nt",j.top_gap,j.bottom_gap)}</p>
            <p>{format!("Top join at product boundary {} · 3′ {} / 5′ {} · {}",j.top_nick,three(j.top_three_prime),five(j.top_five_prime),if j.top_sealed{"sealed under declared chemistry"}else{"nick remains or processing withheld"})}</p>
            <p>{format!("Bottom join at product boundary {} · 3′ {} / 5′ {} · {}",j.bottom_nick,three(j.bottom_three_prime),five(j.bottom_five_prime),if j.bottom_sealed{"sealed under declared chemistry"}else{"nick remains or processing withheld"})}</p>
            <p>{format!("Adjacent terminal search {} · other retained overlap lengths: {}",if j.terminal_search_complete{"complete within declared bounds"}else{"incomplete"},if j.alternatives.is_empty(){"none".into()}else{j.alternatives.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")})}</p>
        </details>}).collect_view()}
        {physical.map(|summary|view!{<details><summary>"Remaining product nicks and outer termini"</summary>{crate::end_view::summary(&summary)}</details>})}
        <details><summary>"Processing findings and limits"</summary>{report.warnings.iter().chain(&report.limits).map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</details>
        <p class="helper">{if reusable{"Save the actual reconstructed final product with its run, then select that exact saved material for later cloning. Imported or detached strand summaries cannot establish physical material."}else{"Diagnostic DNA and processing evidence can be inspected and retained. This outcome does not establish reusable material for a later cloning step."}}</p>
    </details>}.into_any()
}
