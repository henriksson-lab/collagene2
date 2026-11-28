//! Transient source-coordinate PCR display. No analysis or molecule mutation.
use crate::primers::pcr::pcr_overlay::Kind;
use crate::primers::pcr::pcr_overlay::Overlay;
use crate::primers::pcr::pcr_overlay::Role;
use crate::viewer::map::map_labels;
use crate::viewer::map::map_labels::View;
use collagene_model::Topology;
use std::f64::consts::{PI, TAU};
use web_sys::CanvasRenderingContext2d as Ctx;

#[derive(Clone)]
pub struct Hit {
    pub x: f64,
    pub y: f64,
    pub width: f64,
    pub height: f64,
    pub start: usize,
    pub end: usize,
    pub label: String,
}
impl Hit {
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.x && x < self.x + self.width && y >= self.y && y < self.y + self.height
    }
}
pub fn context(overlay: &Overlay) -> String {
    format!(
        "{}:{}:{}",
        overlay.epoch, overlay.source.document_id, overlay.source.event_id
    )
}
pub fn kind_label(kind: Kind) -> &'static str {
    match kind {
        Kind::ExactProduct => "Selected exact PCR product",
        Kind::MismatchProduct => "Selected conditional mismatch-PCR product",
        Kind::CapturedSelection => "Captured selection comparison · binding not evaluated",
    }
}
pub fn color(role: Role) -> &'static str {
    match role {
        Role::Forward => "#2563a6",
        Role::Reverse => "#c36330",
        Role::Amplified => "#21866a",
        Role::Mismatch => "#c43e48",
        Role::Ambiguous => "#9a6c13",
    }
}
fn offset(role: Role) -> f64 {
    match role {
        Role::Forward => 0.49,
        Role::Reverse => 0.43,
        Role::Amplified => 0.37,
        Role::Mismatch => 0.49,
        Role::Ambiguous => 0.43,
    }
}
fn arrow(ctx: &Ctx, x: f64, y: f64, angle: f64) {
    ctx.begin_path();
    ctx.move_to(x, y);
    ctx.line_to(x - 7. * (angle - 0.5).cos(), y - 7. * (angle - 0.5).sin());
    ctx.line_to(x - 7. * (angle + 0.5).cos(), y - 7. * (angle + 0.5).sin());
    ctx.close_path();
    ctx.fill();
}
/// Paint only source-bound intervals. Added tails appear in the companion
/// legend, never stretched onto the template's coordinate system.
pub fn map(ctx: &Ctx, overlay: &Overlay, view: View, length: usize) -> Vec<Hit> {
    let mut hits = Vec::new();
    if length == 0 || view.width <= 0. || view.height <= 0. {
        return hits;
    }
    let radius = map_labels::radius(view);
    let n = length as f64;
    ctx.save();
    ctx.set_font("11px system-ui");
    let captured = overlay.kind == Kind::CapturedSelection;
    if captured {
        let _ = ctx.set_line_dash(&js_sys::Array::of2(&4.into(), &3.into()));
    }
    for span in &overlay.spans {
        ctx.set_stroke_style_str(color(span.role));
        ctx.set_fill_style_str(color(span.role));
        ctx.set_line_width(if matches!(span.role, Role::Mismatch | Role::Ambiguous) {
            7.
        } else {
            4.
        });
        let middle = (span.start + span.end) as f64 / 2.;
        let (x, y) = if view.topology == Topology::Circular {
            let r = radius * offset(span.role);
            let start = span.start as f64 / n * TAU - PI / 2. + view.rotation;
            let end = span.end as f64 / n * TAU - PI / 2. + view.rotation;
            ctx.begin_path();
            let _ = ctx.arc(view.width / 2., view.height / 2., r, start, end);
            ctx.stroke();
            if matches!(span.role, Role::Forward | Role::Reverse) {
                let a = if span.reverse { start } else { end };
                arrow(
                    ctx,
                    view.width / 2. + r * a.cos(),
                    view.height / 2. + r * a.sin(),
                    a + if span.reverse { -PI / 2. } else { PI / 2. },
                );
            }
            let a = middle / n * TAU - PI / 2. + view.rotation;
            (
                view.width / 2. + r * a.cos(),
                view.height / 2. + r * a.sin(),
            )
        } else {
            let y = view.height / 2.
                - match span.role {
                    Role::Forward | Role::Mismatch => 52.,
                    Role::Reverse | Role::Ambiguous => 34.,
                    Role::Amplified => 16.,
                };
            let x = view.width / 2. + (span.start as f64 / n * 1.6 - 0.8) * radius;
            let end = view.width / 2. + (span.end as f64 / n * 1.6 - 0.8) * radius;
            ctx.begin_path();
            ctx.move_to(x, y);
            ctx.line_to(end, y);
            ctx.stroke();
            if matches!(span.role, Role::Forward | Role::Reverse) {
                arrow(
                    ctx,
                    if span.reverse { x } else { end },
                    y,
                    if span.reverse { PI } else { 0. },
                );
            }
            ((x + end) / 2., y)
        };
        if matches!(span.role, Role::Forward | Role::Reverse) {
            let label = format!(
                "{} {}",
                if span.reverse { "←" } else { "→" },
                span.label.chars().take(48).collect::<String>()
            );
            let width = ctx
                .measure_text(&label)
                .map(|m| m.width())
                .unwrap_or(80.)
                .min(320.);
            ctx.set_text_align("center");
            let _ = ctx.fill_text_with_max_width(&label, x, y - 9., 320.);
            hits.push(Hit {
                x: x - width / 2.,
                y: y - 23.,
                width,
                height: 19.,
                start: span.start,
                end: span.end,
                label: span.label.clone(),
            });
        }
    }
    ctx.set_text_align("left");
    ctx.set_font("11px system-ui");
    ctx.set_fill_style_str("#52645f");
    let _ = ctx.fill_text_with_max_width(
        kind_label(overlay.kind),
        10.,
        view.height - 12.,
        (view.width - 20.).max(1.),
    );
    ctx.restore();
    hits
}
/// Arc/line hits are source spans, not feature identities. Label hits are
/// separately considered after native annotation/restriction labels.
pub fn map_span_hit(
    overlay: &Overlay,
    view: View,
    length: usize,
    x: f64,
    y: f64,
) -> Option<(usize, usize)> {
    if length == 0 {
        return None;
    }
    let radius = map_labels::radius(view);
    if radius <= 0. {
        return None;
    }
    let p = if view.topology == Topology::Circular {
        ((y - view.height / 2.).atan2(x - view.width / 2.) - view.rotation + PI / 2.)
            .rem_euclid(TAU)
            / TAU
    } else {
        ((x - view.width / 2.) / radius + 0.8) / 1.6
    };
    if !(0.0..1.0).contains(&p) {
        return None;
    }
    let pos = (p * length as f64).floor() as usize;
    overlay
        .spans
        .iter()
        .rev()
        .find(|s| {
            if pos < s.start || pos >= s.end {
                return false;
            }
            if view.topology == Topology::Circular {
                ((x - view.width / 2.).hypot(y - view.height / 2.) - radius * offset(s.role)).abs()
                    <= 8.
            } else {
                let lane = view.height / 2.
                    - match s.role {
                        Role::Forward | Role::Mismatch => 52.,
                        Role::Reverse | Role::Ambiguous => 34.,
                        Role::Amplified => 16.,
                    };
                (y - lane).abs() <= 8.
            }
        })
        .map(|s| (s.start, s.end))
}
/// Overlay beneath reference glyphs: amplicon extent and individual diagnostic
/// bases stay in template coordinates; primers receive two dedicated lanes.
pub fn sequence(
    ctx: &Ctx,
    overlay: &Overlay,
    start: usize,
    end: usize,
    x: f64,
    y: f64,
    cell: f64,
) -> Vec<Hit> {
    let mut hits = Vec::new();
    ctx.save();
    ctx.set_font("11px system-ui");
    ctx.set_text_align("left");
    if overlay.kind == Kind::CapturedSelection {
        let _ = ctx.set_line_dash(&js_sys::Array::of2(&4.into(), &3.into()));
    }
    for span in &overlay.spans {
        let left = span.start.max(start);
        let right = span.end.min(end);
        if left >= right {
            continue;
        }
        let sx = x + (left - start) as f64 * cell;
        let width = (right - left) as f64 * cell;
        ctx.set_fill_style_str(color(span.role));
        ctx.set_stroke_style_str(color(span.role));
        match span.role {
            Role::Amplified => {
                ctx.set_global_alpha(0.10);
                ctx.fill_rect(sx, y - 15., width, 39.);
                ctx.set_global_alpha(1.);
            }
            Role::Mismatch | Role::Ambiguous => {
                ctx.set_global_alpha(0.23);
                ctx.fill_rect(sx, if span.reverse { y + 5. } else { y - 15. }, width, 19.);
                ctx.set_global_alpha(1.);
                ctx.set_line_width(2.);
                ctx.stroke_rect(sx, if span.reverse { y + 5. } else { y - 15. }, width, 19.);
            }
            Role::Forward | Role::Reverse => {
                let lane_y = y + if span.role == Role::Forward { 35. } else { 53. };
                ctx.set_line_width(2.);
                ctx.begin_path();
                ctx.move_to(sx, lane_y);
                ctx.line_to(sx + width, lane_y);
                ctx.stroke();
                if (span.reverse && left == span.start) || (!span.reverse && right == span.end) {
                    arrow(
                        ctx,
                        if span.reverse { sx } else { sx + width },
                        lane_y,
                        if span.reverse { PI } else { 0. },
                    );
                }
                if width > 48. {
                    let text = format!(
                        "{} {}",
                        if span.reverse { "←" } else { "→" },
                        span.label.chars().take(48).collect::<String>()
                    );
                    let _ = ctx.fill_text_with_max_width(&text, sx, lane_y - 4., width);
                }
                hits.push(Hit {
                    x: sx,
                    y: lane_y - 15.,
                    width,
                    height: 20.,
                    start: span.start,
                    end: span.end,
                    label: span.label.clone(),
                });
            }
        }
    }
    ctx.restore();
    hits
}
