//! Bounded presentation of native findings. Navigation remains caller guarded.
use super::model::Capture;
use collagene_core::construct_qc::causal::{Finding, PairReport, RegisterSegment, StopEvidence};
use leptos::prelude::*;
type Select = Callback<(Capture, Vec<(usize, usize)>)>;

fn finding(finding: &Finding, capture: &Capture, select: Select) -> AnyView {
    let ranges = finding
        .product_coordinates
        .iter()
        .flat_map(|location| {
            location
                .spans()
                .into_iter()
                .map(|(start, end, _)| (start, end))
        })
        .collect::<Vec<_>>();
    let selectable = !ranges.is_empty();
    let capture = capture.clone();
    view! {
        <details class="qc-finding">
            <summary>{format!("{:?} · {:?}",finding.verdict,finding.rule)}</summary>
            <p>{finding.reason.clone()}</p>
            <p class="helper">
                {if finding.under_declared_annotation_assumptions {
                    "Conditional on your declared CDS roles."
                } else {
                    "No established annotation-role premise."
                }}
            </p>
            <button class="text-tool" disabled=!selectable
                on:click=move |_|select.run((capture.clone(),ranges.clone()))>
                "Select product coordinates"
            </button>
            {finding.citations.iter().map(|citation| view! {
                <p class="helper">{format!("{} · {}",citation.source,citation.locator)}</p>
            }).collect_view()}
        </details>
    }
    .into_any()
}
fn segment(segment: &RegisterSegment) -> AnyView {
    let summary = format!(
        "Product coding bases {}–{} · source anchors {}–{} · frame offset {} bases · recorded join {} · {} recorded causes",
        segment.product_coding_start + 1,
        segment.product_coding_end,
        segment.first_parent_coding_ordinal + 1,
        segment.last_parent_coding_ordinal + 1,
        segment.offset_mod3,
        if segment.native_junction_proven {
            "supported"
        } else {
            "unresolved"
        },
        segment.causes.len(),
    );
    view! { <p>{summary}</p> }.into_any()
}
fn stop(stop: &StopEvidence, capture: &Capture, select: Select) -> AnyView {
    let ranges = stop
        .product_coordinates
        .spans()
        .into_iter()
        .map(|(start, end, _)| (start, end))
        .collect::<Vec<_>>();
    let capture = capture.clone();
    view! {
        <p>{format!("{:?} · product codon {} ({}) · {}",
            stop.verdict,stop.product_codon_index+1,stop.product_codon,stop.reason)}</p>
        <button class="text-tool" on:click=move |_|select.run((capture.clone(),ranges.clone()))>
            "Select codon"
        </button>
    }
    .into_any()
}
pub(super) fn pair(pair: &PairReport, capture: &Capture, page: usize, select: Select) -> AnyView {
    view! {
        <article>
            <h4>"CDS history comparison"</h4>
            <p>{format!("Annotation correspondence: {:?} · {}",
                pair.annotation_correspondence.verdict,pair.annotation_correspondence.reason)}</p>
            <p class="helper">{format!(
                "Mapping complete: {} · genetic code {} · {} register segments · {} stop comparisons",
                pair.mapping_complete,pair.code.descriptor.id,
                pair.register_segments.len(),pair.new_stops.len())}</p>
            {pair.findings.iter().map(|value|finding(value,capture,select)).collect_view()}
            <details>
                <summary>"Reading frame evidence"</summary>
                <p class="helper">"Positions count coding bases from 1 after the declared CDS start offset. Source anchors may bracket changes; frame offsets are modulo 3. An offset alone does not establish how a DNA junction was made."</p>
                {pair.register_segments.iter().skip(page.saturating_mul(20)).take(20).map(segment).collect_view()}
            </details>
            <details>
                <summary>"Internal stop comparisons"</summary>
                {pair.new_stops.iter().skip(page.saturating_mul(20)).take(20)
                    .map(|value|stop(value,capture,select)).collect_view()}
            </details>
        </article>
    }.into_any()
}
