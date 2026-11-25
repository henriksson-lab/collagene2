//! Source-bound measured map labels. Full text and one actual-part anchor are
//! cached once; camera frames copy only small numeric layout inputs.
use crate::annotation_geometry::PreparedGeometry;
use collagene_core::{Molecule, Topology};
use collagene_map_label_layout as label_layout;
use label_layout::{
    Column, Label, Options, Outcome, OverflowReason, Point, Rect, Side, Size, Viewport,
};
use std::{
    collections::BTreeMap,
    f64::consts::{PI, TAU},
};
use web_sys::CanvasRenderingContext2d as Ctx;

const FONT: &str = "13px system-ui";
const LINE: f64 = 16.;
const TEXT_PAD: f64 = 2.;
const DISPLAY_BUDGET: usize = 128;
const MAX_NAME_BYTES: usize = 4 * 1024 * 1024;
const MAX_SINGLE_NAME_BYTES: usize = 4 * 1024;
const MAX_VARIANT_CHARS: usize = 512;

#[derive(Clone, Copy, Debug)]
pub struct View {
    pub width: f64,
    pub height: f64,
    pub rotation: f64,
    pub zoom: f64,
    pub topology: Topology,
}

#[derive(Clone, Default, Debug, PartialEq, Eq)]
pub struct LabelStats {
    pub total: usize,
    pub displayed: usize,
    pub omitted: usize,
    pub abbreviated: usize,
    pub unplaced: usize,
    pub lane_overflow: usize,
    pub limited: bool,
}

pub struct Hit {
    pub feature_index: usize,
    pub rect: Rect,
}

pub struct Frame {
    pub hits: Vec<Hit>,
    pub stats: LabelStats,
}

struct Name {
    feature_index: usize,
    lane: usize,
    midpoint: f64,
    span_fraction: f64,
    inline_color: &'static str,
    stable_rank: u64,
    text: String,
    width: f64,
    left: f64,
    ascent: f64,
    color: String,
}

struct Variant {
    width_key: u64,
    height_key: u64,
    lines: Vec<TextLine>,
    measured: Size,
    abbreviated: bool,
}

struct TextLine {
    text: String,
    left: f64,
    ascent: f64,
}
struct Metrics {
    width: f64,
    left: f64,
    ascent: f64,
}

pub struct Cache {
    names: Vec<Name>,
    total: usize,
    unplaced: usize,
    lane_overflow: usize,
    failed: bool,
    variants: BTreeMap<usize, Variant>,
    candidates: Vec<usize>,
    candidate_key: Option<(Option<usize>, usize)>,
}

fn valid_view(view: View) -> bool {
    view.width.is_finite()
        && view.height.is_finite()
        && view.rotation.is_finite()
        && view.zoom.is_finite()
        && view.width > 0.
        && view.height > 0.
        && view.width <= label_layout::MAX_VIEWPORT
        && view.height <= label_layout::MAX_VIEWPORT
        && (0.1..=10.).contains(&view.zoom)
}

/// Shared geometry scale: normalized feature arcs remain below .615r;
/// linear bars span +/- .8r. The host must use this same radius for GPU scale,
/// coordinate ticks, source-pointer projection and label anchors.
pub fn radius(view: View) -> f64 {
    if !valid_view(view) {
        return 0.;
    }
    let available = match view.topology {
        Topology::Circular => ((view.width - 220.) / 1.32).min((view.height - 76.) / 1.32),
        Topology::Linear => ((view.width - 32.) / 1.6).min(view.height * 0.9),
    };
    available.max(1.) * view.zoom
}

fn columns(view: View) -> Option<Vec<Column>> {
    if !valid_view(view) || view.width < 160. || view.height < 160. {
        return None;
    }
    let (w, h) = (view.width, view.height);
    Some(match view.topology {
        Topology::Circular => {
            let r = radius(view);
            let width = (r * 0.38).clamp(80., 170.);
            let cy = h / 2.;
            let cx = w / 2.;
            vec![
                Column {
                    side: Side::Left,
                    rect: Rect {
                        x: cx - width - 24.,
                        y: cy - r * 0.4,
                        width,
                        height: r * 0.8,
                    },
                },
                Column {
                    side: Side::Right,
                    rect: Rect {
                        x: cx + 24.,
                        y: cy - r * 0.4,
                        width,
                        height: r * 0.8,
                    },
                },
            ]
        }
        Topology::Linear => {
            let cy = h / 2.;
            let top = cy + 0.28 * radius(view) + 12.;
            let bottom = h - 72.;
            if bottom - top < 20. {
                return None;
            }
            vec![Column {
                side: Side::Below,
                rect: Rect {
                    x: 12.,
                    y: top,
                    width: w - 24.,
                    height: bottom - top,
                },
            }]
        }
    })
}

fn metrics(ctx: &Ctx, text: &str) -> Option<Metrics> {
    let measured = ctx.measure_text(text).ok()?;
    let (advance, left, right, ascent, descent) = (
        measured.width(),
        measured.actual_bounding_box_left(),
        measured.actual_bounding_box_right(),
        measured.actual_bounding_box_ascent(),
        measured.actual_bounding_box_descent(),
    );
    if [advance, left, right, ascent, descent]
        .iter()
        .any(|v| !v.is_finite())
        || advance < 0.
        || ascent.max(0.) + descent.max(0.) > LINE
    {
        return None;
    }
    // Left-aligned alphabetic metrics: reserve actual overhang, not just text
    // advance. Taller-than-line glyph runs are explicitly unsupported rather
    // than being painted through neighboring boxes.
    let left = left.max(0.);
    Some(Metrics {
        width: left + advance.max(right).max(0.),
        left,
        ascent: ascent.max(0.),
    })
}

fn measure(ctx: &Ctx, text: &str) -> Option<f64> {
    Some(metrics(ctx, text)?.width)
}

impl Cache {
    /// Complete cache or explicit whole-cache limitation, never a first-feature
    /// prefix. Bounds on all names/counts/ranges precede names/colors/measurement.
    /// The full source remains untouched and available to the inspector.
    pub fn prepare(
        molecule: &Molecule,
        geometry: &PreparedGeometry,
        lanes: &[usize],
        lane_overflow: usize,
        ctx: &Ctx,
    ) -> Self {
        let total = geometry
            .located
            .len()
            .saturating_add(geometry.unplaced.len())
            .saturating_add(lane_overflow);
        let mut cache = Self {
            names: Vec::new(),
            total,
            unplaced: geometry.unplaced.len(),
            lane_overflow,
            failed: false,
            variants: BTreeMap::new(),
            candidates: Vec::new(),
            candidate_key: None,
        };
        if lanes.len() != geometry.located.len() || lanes.iter().any(|n| *n >= 3) {
            cache.failed = true;
            return cache;
        }
        let mut bytes = 0usize;
        let mut spans = 0usize;
        let mut previous = None;
        if total > label_layout::MAX_INPUTS || molecule.sequence.is_empty() {
            cache.failed = true;
            return cache;
        }
        for located in &geometry.located {
            let Some(feature) = molecule.features.get(located.index) else {
                cache.failed = true;
                return cache;
            };
            let name = feature.label();
            bytes = bytes.saturating_add(name.len());
            spans = spans.saturating_add(located.spans.len());
            if bytes > MAX_NAME_BYTES
                || spans > crate::annotation_geometry::MAX_SPANS
                || name.len() > MAX_SINGLE_NAME_BYTES
                || previous.is_some_and(|value| located.index <= value)
                || located.spans.is_empty()
                || located
                    .spans
                    .iter()
                    .any(|&(a, b, _)| a >= b || b > molecule.sequence.len())
            {
                cache.failed = true;
                return cache;
            }
            previous = Some(located.index);
        }
        ctx.save();
        ctx.set_font(FONT);
        ctx.set_text_align("left");
        ctx.set_text_baseline("alphabetic");
        let mut names = Vec::with_capacity(geometry.located.len());
        for (ordinal, located) in geometry.located.iter().enumerate() {
            let feature = &molecule.features[located.index];
            // Longest actual local part avoids putting joined annotations' one
            // label anchor inside an unannotated gap. Ties retain tree order.
            let mut part = located.spans[0];
            for &candidate in &located.spans[1..] {
                if candidate.1 - candidate.0 > part.1 - part.0 {
                    part = candidate;
                }
            }
            let midpoint = part.0 as f64 + (part.1 - part.0) as f64 / 2.;
            let Some(measured) = metrics(ctx, feature.label()) else {
                cache.failed = true;
                break;
            };
            if measured.width + 4. > label_layout::MAX_VIEWPORT {
                cache.failed = true;
                break;
            }
            names.push(Name {
                feature_index: located.index,
                lane: lanes[ordinal],
                midpoint: midpoint / molecule.sequence.len() as f64,
                span_fraction: (part.1 - part.0) as f64 / molecule.sequence.len() as f64,
                inline_color: {
                    let c = crate::renderer::color_for_strand(feature, part.2);
                    if c[0] * 0.2126 + c[1] * 0.7152 + c[2] * 0.0722 < 0.60 {
                        "#ffffff"
                    } else {
                        "#173d35"
                    }
                },
                stable_rank: (part.0 + (part.1 - part.0) / 2) as u64,
                text: feature.label().to_owned(),
                width: measured.width,
                left: measured.left,
                ascent: measured.ascent,
                color: crate::renderer::css_color_for_strand(feature, part.2),
            });
        }
        ctx.restore();
        if !cache.failed {
            cache.names = names;
        }
        cache
    }

    fn initial_frame(&self) -> Frame {
        Frame {
            hits: Vec::new(),
            stats: LabelStats {
                total: self.total,
                omitted: self.total,
                unplaced: self.unplaced,
                lane_overflow: self.lane_overflow,
                limited: self.failed,
                ..LabelStats::default()
            },
        }
    }

    fn select_candidates(&mut self, selected: Option<usize>, rendered: usize) -> bool {
        let count = rendered.min(self.names.len());
        let key = (selected, count);
        if self.candidate_key == Some(key) {
            return true;
        }
        // The generic crate is the sole budget-selection implementation. Dummy
        // projection is internal presentation bookkeeping, not source evidence.
        // Stable ranks and explicit selection priority are the actual source
        // data; camera coordinates are deliberately absent from this pass.
        let labels: Vec<_> = self
            .names
            .iter()
            .take(count)
            .map(|name| Label {
                id: name.feature_index as u64,
                stable_rank: name.stable_rank,
                anchor: Point { x: 0., y: 0. },
                preferred_center: Point { x: 0., y: 0. },
                side: Side::Left,
                measured: Size {
                    width: 1.,
                    height: 1.,
                },
                priority: u32::from(selected == Some(name.feature_index)),
            })
            .collect();
        let options = Options {
            max_inputs: label_layout::MAX_INPUTS,
            max_displayed: DISPLAY_BUDGET,
            gap: 0.,
        };
        let Ok(budget) = label_layout::layout(
            &labels,
            Viewport {
                width: 1000.,
                height: 1000.,
                padding: 0.,
            },
            &[Column {
                side: Side::Left,
                rect: Rect {
                    x: 0.,
                    y: 0.,
                    width: 1000.,
                    height: 1000.,
                },
            }],
            &options,
        ) else {
            return false;
        };
        self.candidates = budget
            .labels
            .iter()
            .enumerate()
            .filter_map(|(i, row)| {
                (!matches!(
                    row.outcome,
                    Outcome::Overflow {
                        reason: OverflowReason::DisplayBudget,
                        ..
                    }
                ))
                .then_some(i)
            })
            .collect();
        self.candidate_key = Some(key);
        true
    }

    pub fn draw(
        &mut self,
        ctx: &Ctx,
        view: &View,
        selected_feature: Option<usize>,
        rendered_geometry_count: usize,
    ) -> Frame {
        let mut frame = self.initial_frame();
        if self.failed {
            return frame;
        }
        let Some(columns) = columns(*view) else {
            frame.stats.limited = true;
            return frame;
        };
        if !self.select_candidates(selected_feature, rendered_geometry_count) {
            frame.stats.limited = true;
            return frame;
        }
        let r = radius(*view);
        let (cx, cy) = (view.width / 2., view.height / 2.);
        let mut labels = Vec::with_capacity(self.candidates.len());
        let mut label_sources = Vec::with_capacity(self.candidates.len());
        for &candidate in &self.candidates {
            let name = &self.names[candidate];
            let layer = name.lane as f64;
            let (inner, outer) = label_layout::radial::feature_band(name.lane);
            let distance = (inner + outer) / 2. * r;
            let inline_width = if view.topology == Topology::Circular {
                name.span_fraction * TAU * distance * 0.82
            } else {
                name.span_fraction * 1.6 * r - 12.
            };
            if name.width + 8. < inline_width
                && (view.topology != Topology::Circular || name.width < distance * 0.65)
            {
                let angle = name.midpoint * TAU - PI / 2. + view.rotation;
                let (anchor, mut rotate) = if view.topology == Topology::Circular {
                    (
                        Point {
                            x: cx + angle.cos() * distance,
                            y: cy + angle.sin() * distance,
                        },
                        angle + PI / 2.,
                    )
                } else {
                    (
                        Point {
                            x: cx + (name.midpoint * 1.6 - 0.8) * r,
                            y: cy + (0.0725 + layer * 0.07) * r,
                        },
                        0.,
                    )
                };
                if rotate.cos() < 0. {
                    rotate += PI;
                }
                let width = name.width * rotate.cos().abs() + LINE * rotate.sin().abs();
                let height = name.width * rotate.sin().abs() + LINE * rotate.cos().abs();
                let rect = Rect {
                    x: anchor.x - width / 2.,
                    y: anchor.y - height / 2.,
                    width,
                    height,
                };
                let collision = frame.hits.iter().any(|h| {
                    h.rect.x < rect.right()
                        && h.rect.right() > rect.x
                        && h.rect.y < rect.bottom()
                        && h.rect.bottom() > rect.y
                });
                if !collision
                    && rect.x >= 0.
                    && rect.y >= 0.
                    && rect.right() <= view.width
                    && rect.bottom() <= view.height
                {
                    ctx.save();
                    ctx.set_font(FONT);
                    let _ = ctx.translate(anchor.x, anchor.y);
                    let _ = ctx.rotate(rotate);
                    ctx.set_text_align("center");
                    ctx.set_text_baseline("middle");
                    ctx.set_fill_style_str(name.inline_color);
                    let success = ctx.fill_text(&name.text, 0., 0.).is_ok();
                    ctx.restore();
                    if success {
                        frame.hits.push(Hit {
                            feature_index: name.feature_index,
                            rect,
                        });
                        frame.stats.displayed += 1;
                        continue;
                    }
                }
            }
            let (anchor, side, preferred_center) = match view.topology {
                Topology::Circular => {
                    let angle = name.midpoint * TAU - PI / 2. + view.rotation;
                    let distance = (inner + outer) / 2. * r;
                    let anchor = Point {
                        x: cx + angle.cos() * distance,
                        y: cy + angle.sin() * distance,
                    };
                    let side = if angle.cos() < 0. {
                        Side::Left
                    } else {
                        Side::Right
                    };
                    let column = columns
                        .iter()
                        .find(|c| c.side == side)
                        .expect("circular sides are constructed together");
                    (
                        anchor,
                        side,
                        Point {
                            x: column.rect.x + column.rect.width / 2.,
                            y: if (anchor.y - cy).abs() < 38. {
                                cy + if anchor.y < cy { -44. } else { 44. }
                            } else {
                                anchor.y
                            },
                        },
                    )
                }
                Topology::Linear => {
                    let anchor = Point {
                        x: cx + (name.midpoint * 1.6 - 0.8) * r,
                        y: cy + (0.0725 + layer * 0.07) * r,
                    };
                    let side = Side::Below;
                    let column = columns
                        .iter()
                        .find(|c| c.side == side)
                        .expect("linear sides are constructed together");
                    (
                        anchor,
                        side,
                        Point {
                            x: anchor.x,
                            y: column.rect.y + column.rect.height / 2.,
                        },
                    )
                }
            };
            label_sources.push(candidate);
            labels.push(Label {
                id: name.feature_index as u64,
                stable_rank: name.stable_rank,
                anchor,
                preferred_center,
                side,
                measured: Size {
                    width: name.width + TEXT_PAD * 2.,
                    height: LINE + TEXT_PAD * 2.,
                },
                priority: u32::from(selected_feature == Some(name.feature_index)),
            });
        }
        let options = Options {
            max_inputs: DISPLAY_BUDGET,
            max_displayed: DISPLAY_BUDGET,
            gap: 4.,
        };
        let viewport = Viewport {
            width: view.width,
            height: view.height,
            padding: 12.,
        };
        let Ok(first) = label_layout::layout(&labels, viewport, &columns, &options) else {
            frame.stats.limited = true;
            return frame;
        };
        // At most 128 current variants plus 128 prior variants. Rotation alone
        // neither changes budget IDs nor repeatedly measures the same text.
        if self.variants.len() > DISPLAY_BUDGET {
            self.variants.clear();
        }
        ctx.save();
        ctx.set_font(FONT);
        ctx.set_text_align("left");
        ctx.set_text_baseline("alphabetic");
        let mut applied = Vec::with_capacity(DISPLAY_BUDGET);
        for (index, result) in first.labels.iter().enumerate() {
            let Outcome::Overflow {
                reason: OverflowReason::TooWide,
                available: Some(available),
            } = result.outcome
            else {
                continue;
            };
            let source_index = label_sources[index];
            let valid = self.variants.get(&source_index).is_some_and(|v| {
                v.width_key == available.width.to_bits()
                    && v.height_key == line_capacity(available) as u64
            });
            if !valid {
                if let Some(variant) = variant(ctx, &self.names[source_index].text, available) {
                    self.variants.insert(source_index, variant);
                } else {
                    frame.stats.limited = true;
                    continue;
                }
            }
            if let Some(variant) = self.variants.get(&source_index) {
                labels[index].measured = variant.measured;
                applied.push(index);
            }
        }
        let final_layout = if applied.is_empty() {
            Some(first)
        } else {
            label_layout::layout(&labels, viewport, &columns, &options).ok()
        };
        let Some(final_layout) = final_layout else {
            ctx.restore();
            frame.stats.limited = true;
            return frame;
        };
        frame.hits.reserve(final_layout.displayed);
        // Leaders may cross labels in a valid geometric layout. Draw every
        // leader first, then opaque map-background text boxes, protecting text
        // from later leaders AND zoomed underlying GPU features. This is a
        // deliberate local occlusion policy, not a geometry exclusion claim.
        ctx.set_stroke_style_str("#bfcfca");
        ctx.set_line_width(0.8);
        for result in &final_layout.labels {
            if let Outcome::Placed { leader, .. } = result.outcome {
                ctx.begin_path();
                ctx.move_to(leader.anchor.x, leader.anchor.y);
                ctx.line_to(leader.endpoint.x, leader.endpoint.y);
                ctx.stroke();
            }
        }
        for (index, result) in final_layout.labels.iter().enumerate() {
            let Outcome::Placed { rect, .. } = result.outcome else {
                continue;
            };
            let source_index = label_sources[index];
            let name = &self.names[source_index];
            let variant = applied
                .binary_search(&index)
                .ok()
                .and_then(|_| self.variants.get(&source_index));
            ctx.set_fill_style_str("#fbfcfb");
            ctx.fill_rect(rect.x, rect.y, rect.width, rect.height);
            ctx.set_fill_style_str(&name.color);
            let success = if let Some(variant) = variant {
                variant.lines.iter().enumerate().all(|(line, text)| {
                    ctx.fill_text(
                        &text.text,
                        rect.x + TEXT_PAD + text.left,
                        rect.y + TEXT_PAD + line as f64 * LINE + text.ascent,
                    )
                    .is_ok()
                })
            } else {
                ctx.fill_text(
                    &name.text,
                    rect.x + TEXT_PAD + name.left,
                    rect.y + TEXT_PAD + name.ascent,
                )
                .is_ok()
            };
            if !success {
                frame.stats.limited = true;
                continue;
            }
            frame.stats.displayed += 1;
            frame.stats.abbreviated += usize::from(variant.is_some_and(|v| v.abbreviated));
            frame.hits.push(Hit {
                feature_index: name.feature_index,
                rect,
            });
        }
        ctx.restore();
        frame.stats.omitted = frame.stats.total.saturating_sub(frame.stats.displayed);
        frame.stats.limited |= frame.stats.omitted > 0;
        frame
    }
}

fn line_capacity(available: Size) -> usize {
    (((available.height - TEXT_PAD * 2.) / LINE).floor() as usize).min(2)
}

fn variant(ctx: &Ctx, text: &str, available: Size) -> Option<Variant> {
    let max_width = available.width - TEXT_PAD * 2.;
    let max_lines = line_capacity(available);
    if max_width <= 0. || max_lines == 0 {
        return None;
    }
    let chars: Vec<char> = text.chars().take(MAX_VARIANT_CHARS).collect();
    let complete = chars.iter().map(|c| c.len_utf8()).sum::<usize>() == text.len();
    let mut at = 0;
    let mut lines = Vec::with_capacity(max_lines);
    let mut abbreviated = !complete;
    let mut width = 0f64;
    for line in 0..max_lines {
        if at == chars.len() {
            break;
        }
        let rest: String = chars[at..].iter().collect();
        // If the bounded prefix fits but there is unseen source text, this is
        // already the final display line: append a measured ellipsis rather
        // than silently exhausting the prefix before reaching max_lines.
        let last = line + 1 == max_lines || (!complete && measure(ctx, &rest)? <= max_width);
        let suffix = if last && (!complete || measure(ctx, &rest)? > max_width) {
            "…"
        } else {
            ""
        };
        if measure(ctx, suffix)? > max_width {
            return None;
        }
        let (mut lo, mut hi) = (0, chars.len() - at);
        while lo < hi {
            let mid = (lo + hi).div_ceil(2);
            let mut candidate: String = chars[at..at + mid].iter().collect();
            candidate.push_str(suffix);
            if measure(ctx, &candidate)? <= max_width {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        if lo == 0 && suffix.is_empty() {
            return None;
        }
        if !last && at + lo < chars.len() {
            if let Some(split) = chars[at..at + lo]
                .iter()
                .rposition(|c| c.is_whitespace())
                .filter(|n| *n > 0)
            {
                lo = split;
            } else {
                // An indivisible name gets one explicit abbreviated line.
                let mut fitted: String = chars[at..at + lo].iter().collect();
                while !fitted.is_empty() && measure(ctx, &format!("{fitted}…"))? > max_width {
                    fitted.pop();
                }
                fitted.push('…');
                let measured = metrics(ctx, &fitted)?;
                width = width.max(measured.width);
                lines.push(TextLine {
                    text: fitted,
                    left: measured.left,
                    ascent: measured.ascent,
                });
                abbreviated = true;
                break;
            }
        }
        let mut fitted: String = chars[at..at + lo].iter().collect();
        fitted.push_str(suffix);
        let measured = metrics(ctx, &fitted)?;
        width = width.max(measured.width);
        lines.push(TextLine {
            text: fitted,
            left: measured.left,
            ascent: measured.ascent,
        });
        at += lo;
        while at < chars.len() && chars[at].is_whitespace() {
            at += 1;
        }
        if !suffix.is_empty() {
            abbreviated = true;
            break;
        }
    }
    if at < chars.len() {
        abbreviated = true;
    }
    if lines.is_empty() {
        return None;
    }
    Some(Variant {
        width_key: available.width.to_bits(),
        height_key: max_lines as u64,
        measured: Size {
            width: width + TEXT_PAD * 2.,
            height: lines.len() as f64 * LINE + TEXT_PAD * 2.,
        },
        lines,
        abbreviated,
    })
}
