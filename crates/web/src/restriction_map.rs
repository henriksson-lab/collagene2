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
const BOLD_FONT: &str = "600 12px system-ui";

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
    tooltip: String,
}
struct Group {
    members: Vec<usize>,
    text: String,
    width: f64,
    left: f64,
    ascent: f64,
}
pub struct Cache {
    rows: Vec<Site>,
    groups: Vec<Group>,
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
            groups: Vec::new(),
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
            let text = format!("{}{}", site.enzyme, if site.certain { "" } else { " ?" });
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
                tooltip,
            });
        }
        if !cache.failed {
            let mut by_anchor = BTreeMap::<u64, Vec<usize>>::new();
            for (i, site) in cache.rows.iter().enumerate() {
                by_anchor.entry(site.rank).or_default().push(i);
            }
            for members in by_anchor.into_values() {
                let mut names = members
                    .iter()
                    .map(|i| cache.rows[*i].text.as_str())
                    .collect::<Vec<_>>();
                names.sort_unstable();
                names.dedup();
                let mut text = names.iter().take(3).copied().collect::<Vec<_>>().join(", ");
                if names.len() > 3 {
                    text.push_str(&format!(", +{}", names.len() - 3));
                }
                let Ok(m) = ctx.measure_text(&text) else {
                    cache.failed = true;
                    break;
                };
                ctx.set_font(BOLD_FONT);
                let bold = ctx.measure_text(&text);
                ctx.set_font(FONT);
                let Ok(bold) = bold else {
                    cache.failed = true;
                    break;
                };
                let left = m
                    .actual_bounding_box_left()
                    .max(bold.actual_bounding_box_left())
                    .max(0.);
                let width = left
                    + m.width()
                        .max(m.actual_bounding_box_right())
                        .max(bold.width())
                        .max(bold.actual_bounding_box_right())
                        .max(0.);
                let ascent = m
                    .actual_bounding_box_ascent()
                    .max(bold.actual_bounding_box_ascent())
                    .max(0.);
                if !width.is_finite()
                    || !ascent.is_finite()
                    || width > layout::MAX_VIEWPORT
                    || ascent > 16.
                {
                    cache.failed = true;
                    break;
                }
                cache.groups.push(Group {
                    members,
                    text,
                    width,
                    left,
                    ascent,
                });
            }
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
            .groups
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let s = &self.rows[g.members[0]];
                Label {
                    id: i as u64,
                    stable_rank: s.rank,
                    anchor: Point { x: 0., y: 0. },
                    preferred_center: Point { x: 0., y: 0. },
                    side: Side::Above,
                    measured: Size {
                        width: 1.,
                        height: 1.,
                    },
                    priority: u32::from(
                        g.members
                            .iter()
                            .any(|i| selected == Some(self.rows[*i].index)),
                    ),
                }
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
        let circular = view.topology == Topology::Circular;
        let columns = vec![Column {
            side: Side::Above,
            rect: Rect {
                x: 12.,
                y: 12.,
                width: (view.width - 24.).max(1.),
                height: (cy - 44.).max(20.),
            },
        }];
        let mut input = Vec::with_capacity(self.candidates.len());
        for &i in &self.candidates {
            let group = &self.groups[i];
            let site = &self.rows[group.members[0]];
            let anchor = project(site.anchor, 0.66);
            input.push(Label {
                id: i as u64,
                stable_rank: site.rank,
                anchor,
                preferred_center: Point {
                    x: anchor.x,
                    y: cy - 45.,
                },
                side: Side::Above,
                measured: Size {
                    width: group.width + 4.,
                    height: 20.,
                },
                priority: u32::from(
                    group
                        .members
                        .iter()
                        .any(|i| selected == Some(self.rows[*i].index)),
                ),
            });
        }
        let viewport = Viewport {
            width: view.width,
            height: view.height,
            padding: 8.,
        };
        let options = layout::Options {
            max_inputs: LABELS,
            max_displayed: LABELS,
            gap: 4.,
        };
        let placed = if circular {
            layout::radial::outside_circle(
                &input,
                viewport,
                Point { x: cx, y: cy },
                radius * 0.66,
                &options,
            )
        } else {
            layout::layout(&input, viewport, &columns, &options)
        };
        if let Ok(labels) = placed {
            ctx.set_stroke_style_str("#b4b4b4");
            ctx.set_line_width(0.8);
            for result in &labels.labels {
                if let Outcome::Placed { leader, .. } = result.outcome {
                    ctx.begin_path();
                    ctx.move_to(leader.anchor.x, leader.anchor.y);
                    if circular {
                        ctx.line_to(leader.endpoint.x, leader.anchor.y);
                    }
                    ctx.line_to(leader.endpoint.x, leader.endpoint.y);
                    ctx.stroke();
                }
            }
            ctx.set_font(FONT);
            ctx.set_text_align("left");
            ctx.set_text_baseline("alphabetic");
            for (ordinal, result) in labels.labels.iter().enumerate() {
                if let Outcome::Placed { rect, .. } = result.outcome {
                    let group = &self.groups[self.candidates[ordinal]];
                    let site = &self.rows[*group
                        .members
                        .iter()
                        .find(|i| selected == Some(self.rows[**i].index))
                        .unwrap_or(&group.members[0])];
                    ctx.set_font(
                        if group
                            .members
                            .iter()
                            .any(|i| selected == Some(self.rows[*i].index))
                        {
                            BOLD_FONT
                        } else {
                            FONT
                        },
                    );
                    ctx.set_fill_style_str("#fbfcfb");
                    ctx.fill_rect(rect.x, rect.y, rect.width, rect.height);
                    ctx.set_fill_style_str(if selected == Some(site.index) {
                        "#151515"
                    } else if site.certain {
                        "#333333"
                    } else {
                        "#756047"
                    });
                    if ctx
                        .fill_text(
                            &group.text,
                            rect.x + 2. + group.left,
                            rect.y + 2. + group.ascent,
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
                        frame.labels += group.members.len();
                    }
                }
            }
        }
        ctx.restore();
        frame.limited |= frame.markers < frame.total || frame.labels < frame.total;
        frame
    }
}
