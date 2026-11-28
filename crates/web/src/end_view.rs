//! Small diagrams over fresh typed physical ends; no new reconstruction or pairing.
use crate::end_geometry::{self, Side};
use collagene_core::chemistry::{DuplexSummary, End, Polarity};
use collagene_model::Topology;
use leptos::prelude::*;

fn endpoint(end: Option<&End>, side: Side) -> AnyView {
    let side_name = match side {
        Side::Left => "Left end",
        Side::Right => "Right end",
    };
    let Some(end) = end else {
        return view!{<p class="physical-end-unavailable">{format!("{side_name}: physical geometry unavailable")}</p>}.into_any();
    };
    let layout = match end_geometry::layout(end, side) {
        Ok(layout) => layout,
        Err(error) => return view!{<p class="physical-end-unavailable" role="status">{format!("{side_name}: {error}")}</p>}.into_any(),
    };
    let x = |column: f64| 44. + 12. * column;
    let width = 88. + 12. * layout.width_columns as f64;
    let top_tip = x(layout.top_tip);
    let bottom_tip = x(layout.bottom_tip);
    let body_edge = x(layout.body_edge);
    let body_center = match side {
        Side::Left => (top_tip.max(bottom_tip) + body_edge) / 2.,
        Side::Right => (top_tip.min(bottom_tip) + body_edge) / 2.,
    };
    let shape = match end.polarity {
        Polarity::Blunt => "Blunt".to_owned(),
        Polarity::FivePrime => format!("5′ overhang · {} nt", layout.overhang_length),
        Polarity::ThreePrime => format!("3′ overhang · {} nt", layout.overhang_length),
    };
    let caption = format!(
        "{side_name} · {shape}{}",
        if layout.abbreviated {
            " · abbreviated"
        } else {
            ""
        }
    );
    view! {
        <figure class="physical-end-diagram">
            <svg viewBox=format!("0 0 {width} 108") style:min-width=format!("{width}px") role="img" aria-label=caption.clone()>
                <title>{caption.clone()}</title>
                <line x1=top_tip y1="40" x2=body_edge y2="40" class="physical-end-top"/>
                <line x1=bottom_tip y1="72" x2=body_edge y2="72" class="physical-end-bottom"/>
                <text x=top_tip y="15" text-anchor="middle" class="physical-end-terminus">{layout.top_terminus}</text>
                <text x=bottom_tip y="99" text-anchor="middle" class="physical-end-terminus">{layout.bottom_terminus}</text>
                <text x=body_center y="59" text-anchor="middle" class="physical-end-interior">"duplex"</text>
                {layout.glyphs.into_iter().map(|glyph| view!{
                    <text x=x(glyph.column) y=if glyph.top { "34" } else { "66" } text-anchor="middle" class="physical-end-base">{glyph.symbol.to_string()}</text>
                }).collect_view()}
            </svg>
            <figcaption>{caption}</figcaption>
        </figure>
    }.into_any()
}

pub fn pair(left: Option<&End>, right: Option<&End>) -> AnyView {
    view!{
        <div class="physical-end-pair" role="group" aria-label="Physical fragment ends in assembly orientation">
            {endpoint(left, Side::Left)}{endpoint(right, Side::Right)}
        </div>
        <p class="physical-end-key">"End geometry only; interior and nicks not drawn. Bottom letters run 3′→5′ from left to right. P: phosphate · OH: hydroxyl · ?: unknown."</p>
    }.into_any()
}

pub fn summary(summary: &DuplexSummary) -> AnyView {
    let drawing = match &summary.ends {
        Some((left, right)) => pair(Some(left), Some(right)),
        None if summary.topology == Topology::Circular => view!{
            <p class="physical-end-unavailable">"Circular duplex: no external end pair."</p>
        }.into_any(),
        None => view!{<p class="physical-end-unavailable">"Physical free-end geometry unavailable."</p>}.into_any(),
    };
    view!{
        {drawing}
        <p class="physical-end-nicks">{format!("Retained nicks in reconstructed model: {}", summary.nicks.len())}</p>
        <details class="physical-end-metadata"><summary>"Full physical state and nick positions"</summary><pre>{format!("Bottom reference offset: {}\n{}",summary.bottom_start,summary.describe())}</pre></details>
    }.into_any()
}

pub fn left(end: Option<&End>) -> AnyView {
    endpoint(end, Side::Left)
}
pub fn right(end: Option<&End>) -> AnyView {
    endpoint(end, Side::Right)
}
