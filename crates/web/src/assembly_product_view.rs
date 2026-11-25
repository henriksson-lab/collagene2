//! Read-only final-step reconstruction review. The parent supplies a fresh private
//! verifier result; a serialized report is never accepted as installation proof.
use crate::{
    assembly_junctions_view::AssemblyJunctionsPanel, renderer::MapRenderer,
    sequence_view::SequenceView,
};
use collagene_core::{
    Location, Molecule,
    chemistry::DuplexSummary,
    verification::{self, DesignPlan, Operation, Source, Status, VerificationReport},
};
use leptos::{html, prelude::*};
use std::{rc::Rc, sync::Arc};
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
fn source_label(source: &Source) -> String {
    match source {
        Source::Input(id) => format!("original input {}", clip(id, 80)),
        Source::Step(index) => format!("step {}", index + 1),
    }
}
fn operation_label(operation: &Operation) -> String {
    match operation {
        Operation::Pcr {
            template,
            forward,
            reverse,
        } => format!(
            "Exact PCR of {} · {} / {}",
            source_label(template),
            clip(&forward.name, 60),
            clip(&reverse.name, 60)
        ),
        Operation::Digest {
            template,
            enzymes,
            fragment_index,
        } => format!(
            "Digest {} · {} · fragment {} (zero based)",
            source_label(template),
            enzymes
                .iter()
                .take(8)
                .map(|name| clip(name, 32))
                .collect::<Vec<_>>()
                .join(", "),
            fragment_index
        ),
        Operation::CatalogDigest {
            template,
            fragment_index,
            ..
        } => format!(
            "Catalog-pinned digest of {} · fragment {} (zero based)",
            source_label(template),
            fragment_index
        ),
        Operation::RepairEnds { template, method } => format!(
            "Explicit end repair of {} · {:?}",
            source_label(template),
            method
        ),
        Operation::TreatEnds {
            template,
            treatment,
        } => format!(
            "Explicit terminal treatment of {} · {:?}",
            source_label(template),
            treatment
        ),
        Operation::DeclareSourceTermini {
            template,
            left,
            right,
        } => format!(
            "Original duplex terminal assumption for {} · left {:?} · right {:?}",
            source_label(template),
            left,
            right
        ),
        Operation::ProfiledOverlap {
            fragments,
            circular,
            options,
            ..
        } => format!(
            "Ordered {:?} overlap · {} fragments · {}",
            options.method,
            fragments.len(),
            if *circular { "circular" } else { "linear" }
        ),
        Operation::Overlap {
            fragments,
            circular,
            ..
        } => format!(
            "Ordered exact overlap · {} fragments · {}",
            fragments.len(),
            if *circular { "circular" } else { "linear" }
        ),
        Operation::Ligate {
            fragments,
            circular,
            allow_nicked,
        } => format!(
            "Ordered ligation · {} fragments · {} · retained nicks {}",
            fragments.len(),
            if *circular { "circular" } else { "linear" },
            if *allow_nicked {
                "permitted"
            } else {
                "not permitted"
            }
        ),
        _ => "Other typed preparation; inspect the complete saved plan and protocol.".into(),
    }
}
struct Snapshot {
    plan: DesignPlan,
    target: Molecule,
    report: VerificationReport,
    plan_revision: String,
    target_revision: String,
    product_revision: Option<String>,
    error: Option<String>,
    physical: Option<DuplexSummary>,
}
impl Snapshot {
    fn new(plan: DesignPlan, target: Molecule, report: VerificationReport) -> Self {
        let validation = (|| -> Result<(String, String, String), String> {
            collagene_core::assembly_junctions::source_preflight(&plan, &target)
                .map_err(|e| e.to_string())?;
            let plan_revision = verification::plan_revision(&plan);
            let target_revision = target.revision();
            if plan_revision.is_empty()
                || report.plan_revision != plan_revision
                || report.target_revision != target_revision
                || plan.target_revision != target_revision
            {
                return Err("Plan or independent-target revision does not match this reconstruction receipt.".into());
            }
            if report.engine_version != verification::ENGINE_VERSION {
                return Err("This reconstruction uses another scientific engine. Reconstruct the current plan again before reviewing its final product.".into());
            }
            if plan.steps.is_empty() || report.products.len() != plan.steps.len() {
                return Err(format!(
                    "Final product unavailable: {} of {} step products were reconstructed. Leading preparation products are not the final assembly.",
                    report.products.len(),
                    plan.steps.len()
                ));
            }
            let product = report.products.last().ok_or("Final product unavailable")?;
            collagene_core::model::bounds::model_preflight(&[product])
                .map_err(|e| e.to_string())?;
            product.validate().map_err(|e| e.to_string())?;
            Ok((plan_revision, target_revision, product.revision()))
        })();
        let (plan_revision, target_revision, product_revision, error) = match validation {
            Ok((p, t, m)) => (p, t, Some(m), None),
            Err(error) => (
                report.plan_revision.clone(),
                report.target_revision.clone(),
                None,
                Some(error),
            ),
        };
        // Display a descriptive receipt from this fresh report only. It is never
        // used to reconstruct chemistry, authorize installation or choose DNA.
        let physical = if product_revision.is_some() {
            report
                .products
                .last()
                .and_then(|m| m.metadata.get("duplex_model"))
                .filter(|text| text.len() <= 1024 * 1024)
                .and_then(|text| serde_json::from_str::<DuplexSummary>(text).ok())
                .filter(|summary| {
                    summary.nicks.len() <= 10_000
                        && summary.ends.as_ref().is_none_or(|(a, b)| {
                            a.overhang.len() <= 10_000 && b.overhang.len() <= 10_000
                        })
                })
        } else {
            None
        };
        Self {
            plan,
            target,
            report,
            plan_revision,
            target_revision,
            product_revision,
            error,
            physical,
        }
    }
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
    let available = snapshot.with_value(|s| s.product_revision.is_some());
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
        <header class="assembly-product-header"><h2>"Reconstructed product"</h2>
            <p>{snapshot.with_value(|s|format!("Fresh final-step snapshot · {:?} · {}",s.report.status,s.report.engine_version))}</p>
            <p class="helper">"This view shows actual reconstructed reference DNA, separately from intended DNA. It does not establish experimental success, exhaustive product uniqueness or newly prepared material."</p>
        </header>
        <crate::assembly_outputs::AssemblyOutputs capture identity epoch/>
        {on_finalize.map(|finalize|view!{<button class="button primary" disabled=move||!available||!current.get() on:click=move|_|{let result=capture.run(()).and_then(|(plan,_)|finalize.run(collagene_core::transfer::Candidate::new(plan)));if let Err(error)=result{status.set(error);}}>"Create assembly outputs"</button>})}
        <p class="assembly-product-status" role="status">{move||if !current.get(){"Workspace context changed. This snapshot is historical; reconstruct again before opening or inspecting it.".to_string()}else{status.get()}}</p>
        {snapshot.with_value(|s|s.error.clone()).map(|error|view!{<p class="assembly-product-status" role="alert">{error}</p>})}
        {snapshot.with_value(|s|if s.product_revision.is_some(){s.report.products.last().map(|product|{
            let product=product.clone();let product_identity=format!("{}:{}",s.plan_revision,s.product_revision.as_deref().unwrap_or_default());
            view!{<div>
                <p>{format!("Final step {} · {} bases · {:?} · {} annotations",s.plan.steps.len(),product.sequence.len(),product.topology,product.features.len())}</p>
                {(s.report.status==Status::Invalid).then(||view!{<p class="assembly-product-status" role="alert">"The plan/independent target was not accepted. The actual final reconstructed DNA is shown for diagnosis; it must not be mistaken for the intended target."</p>})}
                {(s.report.status==Status::Incomplete).then(||view!{<p class="helper">"Reconstruction is incomplete under the declared model. This reference view does not resolve its missing evidence or assumptions."</p>})}
                <ProductCanvas product=product identity=product_identity current=current on_message=Callback::new(move|text|status.set(text))/>
                <button class="button secondary" disabled=move||!current.get() on:click=open>"Open reconstructed product reference"</button>
            </div>}
        })}else{None})}
        <details class="assembly-product-details"><summary>"Sources, preparation history and verification findings"</summary>
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
        <Show when=move||junctions_open.get()&&available&&current.get()><AssemblyJunctionsPanel capture=capture identity=identity epoch=epoch/></Show>
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
fn ProductCanvas(
    product: Molecule,
    identity: String,
    current: Memo<bool>,
    on_message: Callback<String>,
) -> impl IntoView {
    let source = StoredValue::new(Arc::new(product));
    // One memo-owned model for the existing virtual viewer. Source access in
    // this component stays borrowed; selection changes never clone full DNA.
    let molecule = Memo::new(move |_| source.with_value(|source| source.as_ref().clone()));
    let revision = StoredValue::new(source.with_value(|source| source.revision()));
    let identity = StoredValue::new(identity);
    let selection = RwSignal::new(None::<(usize, usize)>);
    let selected_feature = RwSignal::new(None::<usize>);
    let regions = RwSignal::new(None::<(String, Vec<(usize, usize)>)>);
    let features_open = RwSignal::new(false);
    let feature_page = RwSignal::new(0usize);
    let mode = RwSignal::new(Mode::Split);
    let map_select = RwSignal::new(false);
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
                    Rc::new(|_| {}),
                )
            });
            match created {
                Ok(map) => {
                    map.set_interaction_context(&identity.get_value());
                    map.set_selection_mode(map_select.get_untracked());
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
                return Err(
                    "Choose a nonempty interval within the actual reconstructed product.".into(),
                );
            }
            Ok((a, b))
        })();
        match parsed {
            Ok((a, b)) => {
                selected_feature.set(None);
                regions.set(None);
                selection.set(Some((a, b)));
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
    view! {<div class="assembly-product-stage">
        <div class="assembly-product-controls" role="group" aria-label="Reconstructed product view">
            <button class="text-tool" aria-pressed=move||(mode.get()==Mode::Sequence).to_string() on:click=move|_|mode.set(Mode::Sequence)>"Sequence"</button>
            <button class="text-tool" aria-pressed=move||(mode.get()==Mode::Map).to_string() on:click=move|_|mode.set(Mode::Map)>"Map"</button>
            <button class="text-tool" aria-pressed=move||(mode.get()==Mode::Split).to_string() on:click=move|_|mode.set(Mode::Split)>"Split"</button>
            <label>"Actual product start"<input type="number" min="1" aria-label="Reconstructed product start" prop:value=move||start.get() on:input=move|event|start.set(event_target_value(&event))/></label>
            <label>"End (inclusive)"<input type="number" min="1" aria-label="Reconstructed product end" prop:value=move||end.get() on:input=move|event|end.set(event_target_value(&event))/></label>
            <button class="text-tool" on:click=navigate>"Show actual interval"</button>
            <button class="text-tool" on:click=move|_|{selected_feature.set(None);regions.set(None);selection.set(None);}>"Clear highlight"</button>
            <button class="text-tool" on:click=move|_|features_open.update(|value|*value=!*value)>"Browse product annotations"</button>
        </div>
        <p class="helper">"Read-only reconstructed snapshot. Select bases or use map, coordinates and annotations to navigate; open a separate reference document to edit. Coordinates refer to actual final DNA, even when the intended target differs."</p>
        <div class="assembly-product-map" hidden=move||mode.get()==Mode::Sequence>
            <div class="map-controls"><button class="text-tool" on:click=move|_|if let Some(map)=renderer.get_value(){map.zoom(1./1.15);}>"Zoom out"</button><button class="text-tool" on:click=move|_|if let Some(map)=renderer.get_value(){map.reset();}>"Fit map"</button><button class="text-tool" on:click=move|_|if let Some(map)=renderer.get_value(){map.zoom(1.15);}>"Zoom in"</button>
                <button class="text-tool" aria-pressed=move||map_select.get().to_string() on:click=move|_|{map_select.set(true);if let Some(map)=renderer.get_value(){map.set_selection_mode(true);}}>"Select map range"</button>
                <button class="text-tool" aria-pressed=move||(!map_select.get()).to_string() on:click=move|_|{map_select.set(false);if let Some(map)=renderer.get_value(){map.set_selection_mode(false);}}>"Rotate / inspect"</button>
            </div>
            <p role="status">{move||map_error.get()}</p>
            <div class="map-stage" style="height:400px;min-height:260px;position:relative">
                <canvas node_ref=canvas aria-label="Actual reconstructed product map; drag to rotate or select, click annotations to inspect" style="width:100%;height:100%;display:block"
                    on:pointerdown=move|event:web_sys::PointerEvent|{if !current.get_untracked()||event.button()!=0{return;}if let Some(canvas)=canvas.get(){let _=canvas.set_pointer_capture(event.pointer_id());}if let Some(map)=renderer.get_value(){map.pointer_down(event.client_x()as f64,event.client_y()as f64);}}
                    on:pointermove=move|event:web_sys::PointerEvent|if current.get_untracked(){if let Some(map)=renderer.get_value(){map.pointer_move(event.client_x()as f64,event.client_y()as f64);}}
                    on:pointerup=move|event:web_sys::PointerEvent|{if let Some(map)=renderer.get_value(){if current.get_untracked(){if let Some(hit)=map.pointer_up(event.client_x()as f64,event.client_y()as f64){if let Some(feature)=hit.feature_index{show_feature.run(feature);}else{selected_feature.set(None);selection.set(if hit.spans.len()==1{Some(hit.span)}else{None});regions.set(Some((revision.get_value(),hit.spans)));}}}else{map.pointer_cancel();}}if let Some(canvas)=canvas.get(){let _=canvas.release_pointer_capture(event.pointer_id());}}
                    on:pointercancel=move|_|if let Some(map)=renderer.get_value(){map.pointer_cancel();}
                    on:lostpointercapture=move|_|if let Some(map)=renderer.get_value(){map.pointer_cancel();}
                    on:wheel=move|event:web_sys::WheelEvent|{event.prevent_default();if let Some(map)=renderer.get_value(){map.zoom(if event.delta_y()>0.{0.95}else{1.05});}}/>
                <canvas node_ref=labels class="map-labels" aria-hidden="true" style="position:absolute;inset:0;width:100%;height:100%;pointer-events:none"/>
                <p class="map-geometry-note" role="status">"Overview simplified; browse annotations for omitted geometry."</p>
            </div>
        </div>
        <div class="assembly-product-sequence" hidden=move||mode.get()==Mode::Map>
            <SequenceView molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted search_regions=regions
                on_show_features=Callback::new(move|()|features_open.set(true)) identity=view_identity locked=locked read_only=true shared_restrictions=false
                on_edit=Callback::new(|_|Err("This reconstructed snapshot is read-only; open a separate reference document to edit.".into()))
                on_history=Callback::new(|_|Err("Source history is unchanged by this read-only review.".into())) on_message=on_message
                restriction_sites=no_restrictions restriction_visible=restriction_visible
                on_show_restrictions=Callback::new(move|()|on_message.run("Restriction choices belong to the source/assembly preparation; this snapshot does not run an independent restriction scan.".into()))/>
        </div>
        <Show when=move||features_open.get()><section class="assembly-product-features" aria-label="Reconstructed product annotations">
            <h3>"Actual product annotations"</h3>
            <p>{source.with_value(|source|format!("{} annotations retained; projected or uncertain features may require review.",source.features.len()))}</p>
            <button class="text-tool" disabled=move||feature_page.get()==0 on:click=move|_|feature_page.update(|page|*page=page.saturating_sub(1))>"Previous annotations"</button>
            <button class="text-tool" disabled=move||source.with_value(|source|(feature_page.get()+1)*50>=source.features.len()) on:click=move|_|feature_page.update(|page|*page+=1)>"Next annotations"</button>
            <ul>{move||source.with_value(|source|source.features.iter().enumerate().skip(feature_page.get()*50).take(50).map(|(index,feature)|view!{<li><button class="text-tool" aria-pressed=move||(selected_feature.get()==Some(index)).to_string() on:click=move|_|show_feature.run(index)>{format!("{} · {}{}",clip(feature.label(),160),clip(&feature.kind,48),if feature.needs_review{" · review required"}else{""})}</button></li>}).collect_view())}</ul>
            {move||selected_feature.get().and_then(|index|source.with_value(|source|source.features.get(index).map(|feature|view!{<article class="annotation-details"><h4>{clip(feature.label(),256)}</h4><p>{format!("Original product feature index {} · {} · review {}",index,clip(&feature.kind,80),feature.needs_review)}</p><dl>{feature.qualifiers.iter().take(30).map(|(key,value)|view!{<dt>{clip(key,80)}</dt><dd>{value.as_deref().map(|value|clip(value,512)).unwrap_or_else(||"present".into())}</dd>}).collect_view()}</dl><p class="helper">"Bounded qualifier preview; the opened native reference retains complete annotations and source evidence. Partial/fuzzy/ordered drawing bounds do not establish a complete feature sequence; unresolved locations are not placed."</p></article>})))}
        </section></Show>
    </div>}
}
