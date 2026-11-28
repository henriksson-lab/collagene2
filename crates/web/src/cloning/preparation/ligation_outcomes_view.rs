//! Documentary view of the fresh verifier's bounded ligation outcomes.
use collagene_core::{
    cloning::ligation::outcomes, planning::plan::Source,
    verification::evidence::LigationOutcomeEvidence,
};
use leptos::prelude::*;
use std::sync::Arc;
fn path(value: &[outcomes::Placement]) -> String {
    value
        .iter()
        .map(|p| {
            format!(
                "{} ({})",
                p.key,
                if p.reverse { "reverse" } else { "forward" }
            )
        })
        .collect::<Vec<_>>()
        .join(" → ")
}
fn limit(value: &outcomes::Limit) -> &'static str {
    match value {
        outcomes::Limit::PathChecks => "path checks",
        outcomes::Limit::Products => "retained products",
        outcomes::Limit::ProductBases => "retained product bases",
        outcomes::Limit::ReportBytes => "report size",
        outcomes::Limit::Work => "work budget",
        outcomes::Limit::NativeMaterial => "native material support",
        outcomes::Limit::EnumerationFragments => "fragment enumeration scope",
    }
}
fn product_view(p: &outcomes::Product, title: String) -> AnyView {
    let actual = p
        .duplex
        .summary()
        .map(|s| s.describe())
        .unwrap_or_else(|e| format!("Geometry unavailable: {e}"));
    view!{<section class="ligation-outcome-product"><h5>{title}</h5><p>{path(&p.path)}</p><p>{format!("{} bases · {}",p.sequence.len(),match p.topology{collagene_model::Topology::Circular=>"circular",collagene_model::Topology::Linear=>"linear"})}</p><pre class="physical-geometry">{actual}</pre><p>{format!("Modified bases: {} known modified · {} unknown eligible · {} known unmodified",p.modifications.known_modified,p.modifications.unknown_eligible,p.modifications.known_unmodified)}</p><details><summary>"Actual strand sequences and material state"</summary><p class="helper">"Both strands are shown in their own 5′ to 3′ orientation. These documentary sequences do not grant material authority."</p><p>"Top strand"</p><pre style="white-space:pre-wrap;overflow-wrap:anywhere">{p.duplex.top.as_str().to_owned()}</pre><p>"Bottom strand"</p><pre style="white-space:pre-wrap;overflow-wrap:anywhere">{p.duplex.bottom.as_str().to_owned()}</pre><p>{format!("Material state {}",p.modifications.state_revision)}</p><p>{format!("{} retained generation records; full records remain in the reconstruction JSON and PDF attachment.",p.modifications.generation.len())}</p></details></section>}.into_any()
}
/// Parent supplies fresh verifier evidence, never an imported claim of isolation.
pub(crate) fn ligation_outcomes_view(
    step_index: usize,
    evidence: &LigationOutcomeEvidence,
) -> AnyView {
    if let Err(error) = evidence.validate_documentary() {
        return view!{<p class="helper">{format!("Ligation outcome evidence cannot be displayed: {error}")}</p>}.into_any();
    }
    let evidence = Arc::new(evidence.clone());
    let selected = RwSignal::new(0usize);
    let count = evidence.report.alternatives.len() + 1;
    let r = &evidence.report;
    let summary = format!(
        "Step {} · ligation alternatives ({})",
        step_index + 1,
        r.alternatives.len()
    );
    let scope = format!(
        "Exactly one copy of each of {} declared fragments; {} topology only. Enumeration fragment limit {}. {} paths checked · {} incompatible · {} unresolved · {} products withheld.",
        r.sources.len(),
        if r.request.circular {
            "circular"
        } else {
            "linear"
        },
        r.request.options.max_enumeration_fragments,
        r.paths_checked,
        r.incompatible_paths,
        r.unresolved_paths,
        r.withheld_products
    );
    let completeness = if r.exhaustive_within_scope {
        "Complete within this declared scope"
    } else if r.enumeration_complete {
        "Enumeration completed with unresolved paths"
    } else {
        "Enumeration incomplete; additional alternatives may exist"
    };
    let limits = if r.limits.is_empty() {
        "No enumeration limit reached.".to_owned()
    } else {
        format!(
            "Limits reached: {}.",
            r.limits.iter().map(limit).collect::<Vec<_>>().join(", ")
        )
    };
    let requested = format!(
        "Requested bounds: {} paths, {} products, {} retained bases, {} report bytes and {} work units. Actual retained bases {} and work {}.",
        r.request.options.max_path_checks,
        r.request.options.max_products,
        r.request.options.max_total_product_bases,
        r.request.options.max_report_bytes,
        r.request.options.max_work,
        r.retained_product_bases,
        r.work
    );
    let choices=(0..count).map(|i|view!{<option value=i.to_string()>{if i==0{"Intended path".to_owned()}else{format!("Alternative {i}")}}</option>}).collect_view();
    let warnings = r
        .warnings
        .iter()
        .cloned()
        .map(|text| view! {<li>{text}</li>})
        .collect_view();
    let unresolved = r
        .unresolved_examples
        .iter()
        .map(|v| view! {<li>{format!("{}: {}",path(&v.path),v.reason)}</li>})
        .collect_view();
    let unresolved_count = format!(
        "Unresolved examples shown {} of {}; excerpt truncated {}. Unresolved or withheld products are not absent alternatives.",
        r.unresolved_examples.len(),
        r.unresolved_paths,
        r.unresolved_examples_truncated
    );
    let bindings=evidence.bindings.iter().map(|b|{let source=match &b.source{Source::Input(id)=>format!("input {id}"),Source::Step(i)=>format!("step {}",i.saturating_add(1))};view!{<li>{format!("{} · {} · source {} · material {}",b.key,source,b.source_revision,b.material_revision)}</li>}}).collect_view();
    let context = r.context_revision.clone();
    view!{<details class="physical-evidence ligation-outcomes"><summary>{summary}</summary><p>{scope}</p><p>{completeness}</p><p>{limits}</p><p class="helper">"This is an inspection of intended and alternative paths, not product selection or a prediction of isolation, probability, yield or target success. Other copy numbers, partial assemblies, vector-only products and multimers are outside this scope."</p><label>"Inspect retained product"<select prop:value=move||selected.get().to_string() on:change=move|event|{if let Ok(i)=event_target_value(&event).parse::<usize>(){if i<count{selected.set(i)}}}>{choices}</select></label>{move||{let i=selected.get();let p=if i==0{Some(&evidence.report.intended)}else{evidence.report.alternatives.get(i-1)};p.map(|p|product_view(p,if i==0{"Intended path".into()}else{format!("Alternative {i}")}))}}<details><summary>"Enumeration scope, limits and sources"</summary><p>{requested}</p><ul>{warnings}</ul><p>{unresolved_count}</p><ul>{unresolved}</ul><ul>{bindings}</ul><p>{format!("Report context {context}")}</p></details><p class="helper">"The reconstruction JSON retains every reported product. Protocol PDFs print every retained product and embed the complete per-step JSON and strand FASTA; an oversized PDF is refused explicitly."</p></details>}.into_any()
}
