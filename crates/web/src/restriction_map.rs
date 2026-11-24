//! Cached restriction recognition/cut display. No science is recomputed here.
use crate::map_labels::{self, View};
use collagene_core::{Topology, restriction_inspection::InspectionSite};
use collagene_map_label_layout::{
    self as layout, Column, Label, Outcome, OverflowReason, Point, Rect, Side, Size, Viewport,
};
use std::{
    collections::BTreeMap,
    f64::consts::{PI, TAU},
    sync::Arc,
};
use web_sys::CanvasRenderingContext2d as Ctx;

const MAX_SITES: usize = 10_000;
const LABELS: usize = 64;
const FONT: &str = "12px system-ui";

pub struct Hit {
    pub index: usize,
    pub span: (usize, usize),
    pub rect: Rect,
    pub label: bool,
    segment: Option<(Point, Point)>,
}
impl Hit {
    pub fn contains(&self, x: f64, y: f64) -> bool {
        if x < self.rect.x || x >= self.rect.right() || y < self.rect.y || y >= self.rect.bottom() {
            return false;
        }
        let Some((a, b)) = self.segment else {
            return true;
        };
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let length = dx * dx + dy * dy;
        let t = if length > 0. {
            ((x - a.x) * dx + (y - a.y) * dy) / length
        } else {
            0.
        }
        .clamp(0., 1.);
        (x - a.x - t * dx).hypot(y - a.y - t * dy) <= 4.
    }
}
pub struct Frame {
    pub hits: Vec<Hit>,
    pub total: usize,
    pub markers: usize,
    pub labels: usize,
    pub limited: bool,
}
struct Site {
    index: usize,
    span: (usize, usize),
    ranges: Vec<(f64, f64)>,
    top: Option<f64>,
    bottom: Option<f64>,
    anchor: f64,
    rank: u64,
    certain: bool,
    text: String,
    width: f64,
    left: f64,
    ascent: f64,
    tooltip: String,
}
pub struct Cache {
    rows: Vec<Site>,
    by_id: BTreeMap<usize, usize>,
    total: usize,
    failed: bool,
    candidates: Vec<usize>,
    selected: Option<Option<usize>>,
}

impl Cache {
    pub fn prepare(sites: &Arc<Vec<InspectionSite>>, n: usize, circular: bool, ctx: &Ctx) -> Self {
        let mut cache = Self {
            rows: Vec::new(),
            by_id: BTreeMap::new(),
            total: sites.len(),
            failed: false,
            candidates: Vec::new(),
            selected: None,
        };
        let mut bytes = 0usize;
        if sites.len() > MAX_SITES || (n == 0 && !sites.is_empty()) {
            cache.failed = true;
            return cache;
        }
        for site in sites.iter() {
            bytes = bytes.saturating_add(site.enzyme.len());
            if bytes > 1024 * 1024
                || site.enzyme.len() > 1024
                || site.spans.is_empty()
                || site.spans.len() > 2
                || site.recognition_start >= n
                || site.recognition_length == 0
                || site.recognition_length > n
                || site.spans.iter().any(|s| s.start >= s.end || s.end > n)
                || site
                    .spans
                    .iter()
                    .try_fold(0usize, |total, s| total.checked_add(s.end - s.start))
                    != Some(site.recognition_length)
                || site.top_cut.is_some_and(|p| p > n || (circular && p == n))
                || site
                    .bottom_cut
                    .is_some_and(|p| p > n || (circular && p == n))
            {
                cache.failed = true;
                return cache;
            }
        }
        ctx.save();
        ctx.set_font(FONT);
        ctx.set_text_align("left");
        ctx.set_text_baseline("alphabetic");
        for site in sites.iter() {
            if cache.by_id.contains_key(&site.index) {
                cache.failed = true;
                break;
            }
            let text = format!(
                "{} {}{}",
                site.enzyme,
                site.recognition_start + 1,
                if site.certain { "" } else { " ?" }
            );
            let Ok(metrics) = ctx.measure_text(&text) else {
                cache.failed = true;
                break;
            };
            let left = metrics.actual_bounding_box_left().max(0.);
            let width = left
                + metrics
                    .width()
                    .max(metrics.actual_bounding_box_right())
                    .max(0.);
            let ascent = metrics.actual_bounding_box_ascent().max(0.);
            if !width.is_finite()
                || !ascent.is_finite()
                || width > layout::MAX_VIEWPORT
                || ascent + metrics.actual_bounding_box_descent().max(0.) > 16.
            {
                cache.failed = true;
                break;
            }
            let bounds = |cut: Option<usize>, raw: i64| {
                cut.map(|p| format!("{p} (boundary between bases)"))
                    .unwrap_or_else(|| format!("unavailable (raw reference boundary {raw})"))
            };
            let tooltip = format!(
                "{} · recognition {} · {} · top {} · bottom {} · {}",
                site.enzyme,
                site.spans
                    .iter()
                    .map(|p| format!("{}–{}", p.start + 1, p.end))
                    .collect::<Vec<_>>()
                    .join(", "),
                if site.reverse {
                    "reverse orientation"
                } else {
                    "forward orientation"
                },
                bounds(site.top_cut, site.top_boundary),
                bounds(site.bottom_cut, site.bottom_boundary),
                if site.certain {
                    "certain recognition; cleavage remains an assay assumption"
                } else {
                    "possible recognition; ambiguous source calls"
                }
            );
            let rank = site
                .top_cut
                .or(site.bottom_cut)
                .unwrap_or(site.recognition_start);
            cache.by_id.insert(site.index, cache.rows.len());
            cache.rows.push(Site {
                index: site.index,
                span: (site.spans[0].start, site.spans[0].end),
                ranges: site
                    .spans
                    .iter()
                    .map(|p| (p.start as f64 / n as f64, p.end as f64 / n as f64))
                    .collect(),
                top: site.top_cut.map(|p| p as f64 / n as f64),
                bottom: site.bottom_cut.map(|p| p as f64 / n as f64),
                anchor: rank as f64 / n as f64,
                rank: rank as u64,
                certain: site.certain,
                text,
                width,
                left,
                ascent,
                tooltip,
            });
        }
        ctx.restore();
        if cache.failed {
            cache.rows.clear();
            cache.by_id.clear();
        }
        cache
    }

    pub fn tooltip(&self, index: usize) -> Option<&str> {
        self.by_id
            .get(&index)
            .and_then(|&i| self.rows.get(i))
            .map(|s| s.tooltip.as_str())
    }

    fn candidates(&mut self, selected: Option<usize>) {
        if self.selected == Some(selected) {
            return;
        }
        let input: Vec<_> = self
            .rows
            .iter()
            .map(|s| Label {
                id: s.index as u64,
                stable_rank: s.rank,
                anchor: Point { x: 0., y: 0. },
                preferred_center: Point { x: 0., y: 0. },
                side: Side::Above,
                measured: Size {
                    width: 1.,
                    height: 1.,
                },
                priority: u32::from(selected == Some(s.index)),
            })
            .collect();
        self.candidates = layout::layout(
            &input,
            Viewport {
                width: 1000.,
                height: 1000.,
                padding: 0.,
            },
            &[Column {
                side: Side::Above,
                rect: Rect {
                    x: 0.,
                    y: 0.,
                    width: 1000.,
                    height: 1000.,
                },
            }],
            &layout::Options {
                max_inputs: MAX_SITES,
                max_displayed: LABELS,
                gap: 0.,
            },
        )
        .map(|r| {
            r.labels
                .iter()
                .enumerate()
                .filter_map(|(i, r)| {
                    (!matches!(
                        r.outcome,
                        Outcome::Overflow {
                            reason: OverflowReason::DisplayBudget,
                            ..
                        }
                    ))
                    .then_some(i)
                })
                .collect()
        })
        .unwrap_or_default();
        self.selected = Some(selected);
    }

    pub fn draw(&mut self, ctx: &Ctx, view: View, selected: Option<usize>) -> Frame {
        let mut frame = Frame {
            hits: Vec::with_capacity(self.rows.len().saturating_mul(2).saturating_add(LABELS)),
            total: self.total,
            markers: 0,
            labels: 0,
            limited: self.failed,
        };
        let radius = map_labels::radius(view);
        if self.failed || radius <= 0. || view.width < 80. || view.height < 80. {
            frame.limited |= self.total > 0;
            return frame;
        }
        self.candidates(selected);
        let (cx, cy) = (view.width / 2., view.height / 2.);
        let project = |fraction: f64, r: f64| {
            if view.topology == Topology::Circular {
                let a = fraction * TAU - PI / 2. + view.rotation;
                Point {
                    x: cx + a.cos() * r * radius,
                    y: cy + a.sin() * r * radius,
                }
            } else {
                Point {
                    x: cx + (fraction * 1.6 - 0.8) * radius,
                    y: cy,
                }
            }
        };
        let solid = js_sys::Array::new();
        let dashed = js_sys::Array::of2(&3.into(), &3.into());
        ctx.save();
        for site in &self.rows {
            let _ = ctx.set_line_dash(if site.certain {
                solid.as_ref()
            } else {
                dashed.as_ref()
            });
            ctx.set_line_width(if selected == Some(site.index) {
                3.
            } else {
                1.5
            });
            ctx.set_stroke_style_str("#8b5ca4");
            for &(a, b) in &site.ranges {
                ctx.begin_path();
                if view.topology == Topology::Circular {
                    let _ = ctx.arc(
                        cx,
                        cy,
                        radius * 0.64,
                        a * TAU - PI / 2. + view.rotation,
                        b * TAU - PI / 2. + view.rotation,
                    );
                } else {
                    let a = project(a, 0.);
                    let b = project(b, 0.);
                    ctx.move_to(a.x, cy - 19.);
                    ctx.line_to(b.x, cy - 19.);
                }
                ctx.stroke();
            }
            for (cut, top) in [(site.top, true), (site.bottom, false)] {
                let Some(cut) = cut else {
                    continue;
                };
                let (a, b) = if view.topology == Topology::Circular {
                    (
                        project(cut, if top { 0.65 } else { 0.59 }),
                        project(cut, if top { 0.72 } else { 0.65 }),
                    )
                } else {
                    let p = project(cut, 0.);
                    (
                        Point {
                            x: p.x,
                            y: cy + if top { -13. } else { 2. },
                        },
                        Point {
                            x: p.x,
                            y: cy + if top { 1. } else { 15. },
                        },
                    )
                };
                ctx.set_stroke_style_str(if top { "#2563a6" } else { "#b45c19" });
                ctx.begin_path();
                ctx.move_to(a.x, a.y);
                ctx.line_to(b.x, b.y);
                ctx.stroke();
                let rect = Rect {
                    x: a.x.min(b.x) - 4.,
                    y: a.y.min(b.y) - 4.,
                    width: (a.x - b.x).abs() + 8.,
                    height: (a.y - b.y).abs() + 8.,
                };
                if rect.right() >= 0.
                    && rect.bottom() >= 0.
                    && rect.x < view.width
                    && rect.y < view.height
                {
                    frame.hits.push(Hit {
                        index: site.index,
                        span: site.span,
                        rect,
                        label: false,
                        segment: Some((a, b)),
                    });
                }
            }
            frame.markers += 1;
        }
        let _ = ctx.set_line_dash(solid.as_ref());
        // RE text occupies central top/bottom bands. Feature labels occupy the
        // side columns (circle), or their own above/below bands (linear).
        let circular = view.topology == Topology::Circular;
        let inset = if circular {
            (view.width * 0.18).clamp(40., 100.) + 16.
        } else {
            12.
        };
        let above_bottom = (cy - 0.38 * radius - 8.).max(32.);
        let top_ok = circular || (12. + above_bottom) / 2. - 18. >= 44.;
        let bottom_ok = circular || view.height >= 360.;
        let mut columns = Vec::new();
        if view.width > inset * 2. + 20. && view.height >= 160. {
            if top_ok {
                columns.push(Column {
                    side: Side::Above,
                    rect: Rect {
                        x: inset,
                        y: 8.,
                        width: view.width - inset * 2.,
                        height: 32.,
                    },
                });
            }
            if bottom_ok {
                columns.push(Column {
                    side: Side::Below,
                    rect: Rect {
                        x: inset,
                        y: view.height - 44.,
                        width: view.width - inset * 2.,
                        height: 32.,
                    },
                });
            }
        }
        let mut input = Vec::with_capacity(self.candidates.len());
        for &i in &self.candidates {
            let site = &self.rows[i];
            let anchor = project(site.anchor, 0.66);
            let side = if top_ok && (!bottom_ok || site.index % 2 == 0) {
                Side::Above
            } else {
                Side::Below
            };
            let y = if side == Side::Above {
                24.
            } else {
                view.height - 28.
            };
            input.push(Label {
                id: site.index as u64,
                stable_rank: site.rank,
                anchor,
                preferred_center: Point { x: anchor.x, y },
                side,
                measured: Size {
                    width: site.width + 4.,
                    height: 20.,
                },
                priority: u32::from(selected == Some(site.index)),
            });
        }
        if let Ok(labels) = layout::layout(
            &input,
            Viewport {
                width: view.width,
                height: view.height,
                padding: 8.,
            },
            &columns,
            &layout::Options {
                max_inputs: LABELS,
                max_displayed: LABELS,
                gap: 4.,
            },
        ) {
            ctx.set_stroke_style_str("#b9b1c3");
            ctx.set_line_width(0.8);
            for result in &labels.labels {
                if let Outcome::Placed { leader, .. } = result.outcome {
                    ctx.begin_path();
                    ctx.move_to(leader.anchor.x, leader.anchor.y);
                    ctx.line_to(leader.endpoint.x, leader.endpoint.y);
                    ctx.stroke();
                }
            }
            ctx.set_font(FONT);
            ctx.set_text_align("left");
            ctx.set_text_baseline("alphabetic");
            for (ordinal, result) in labels.labels.iter().enumerate() {
                if let Outcome::Placed { rect, .. } = result.outcome {
                    let site = &self.rows[self.candidates[ordinal]];
                    ctx.set_fill_style_str("#fbfcfb");
                    ctx.fill_rect(rect.x, rect.y, rect.width, rect.height);
                    ctx.set_fill_style_str(if site.certain { "#674478" } else { "#8b5b25" });
                    if ctx
                        .fill_text(
                            &site.text,
                            rect.x + 2. + site.left,
                            rect.y + 2. + site.ascent,
                        )
                        .is_ok()
                    {
                        frame.hits.push(Hit {
                            index: site.index,
                            span: site.span,
                            rect,
                            label: true,
                            segment: None,
                        });
                        frame.labels += 1;
                    }
                }
            }
        }
        ctx.restore();
        frame.limited |= frame.markers < frame.total || frame.labels < frame.total;
        frame
    }
}
