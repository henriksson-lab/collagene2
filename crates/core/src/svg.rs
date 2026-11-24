//! Portable vector sequence maps built from a bounded, typed coordinate scene.
//! SVG contains no scripts, external resources, raw markup, or executable qualifiers.
use crate::{Alphabet, Error, Location, Molecule, Topology};
use std::{
    collections::BTreeMap,
    f64::consts::{FRAC_PI_2, TAU},
    fmt::Write,
};

pub const MAX_MAP_FEATURES: usize = 512;
pub const MAX_MAP_SPANS: usize = 4096;
pub const MAX_MAP_LABELS: usize = 32;
pub const MAX_SVG_BYTES: usize = 2 * 1024 * 1024;
const LANES: usize = 4;
const LABEL_CHARS: usize = 28;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureSpan {
    /// Zero-based inclusive reference start.
    pub start: usize,
    /// Zero-based exclusive reference end; a circular origin-crossing feature has separate spans.
    pub end: usize,
    pub reverse: bool,
    pub fuzzy_start: bool,
    pub fuzzy_end: bool,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MapFeature {
    pub source_index: usize,
    pub kind: String,
    pub label: String,
    pub color: &'static str,
    /// Biological concatenation order, reversing compound-part order under complement.
    pub spans: Vec<FeatureSpan>,
    pub lane: usize,
    /// Index in spans bearing the biological 3′ arrow. Uncertain or protein features have none.
    pub terminal_span: Option<usize>,
    pub fuzzy: bool,
    pub ordered: bool,
    pub unresolved: bool,
    pub partial: bool,
    pub needs_review: bool,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegendEntry {
    pub name: &'static str,
    pub color: &'static str,
    pub count: usize,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MapLabel {
    pub feature_index: usize,
    pub text: String,
    pub x: f64,
    pub y: f64,
    pub anchor_x: f64,
    pub anchor_y: f64,
    pub right_aligned: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MapScene {
    pub name: String,
    pub revision: String,
    pub length: usize,
    pub alphabet: Alphabet,
    pub topology: Topology,
    pub features: Vec<MapFeature>,
    pub labels: Vec<MapLabel>,
    pub legend: Vec<LegendEntry>,
    pub warnings: Vec<String>,
    pub omitted_features: usize,
    pub spans_omitted: bool,
    pub omitted_labels: usize,
    pub track_collisions: usize,
    /// Main feature-track y for linear scenes; lane offsets are 28 drawing units.
    pub linear_feature_y: f64,
}
fn category(kind: &str) -> (&'static str, &'static str) {
    match kind {
        "CDS" => ("Coding sequence", "#2b8775"),
        "promoter" => ("Promoter", "#d6913a"),
        "rep_origin" => ("Replication origin", "#666bb2"),
        "primer_bind" => ("Primer binding", "#4f8cbf"),
        "terminator" => ("Terminator", "#c1636d"),
        "source" => ("Source", "#b3c4bc"),
        _ => ("Other feature", "#7aa396"),
    }
}
fn abbreviated(text: &str, max: usize) -> String {
    let mut iter = text.chars();
    let mut out = String::new();
    for _ in 0..max {
        if let Some(c) = iter.next() {
            out.push(if c.is_whitespace() { ' ' } else { c });
        } else {
            return out;
        }
    }
    if iter.next().is_some() {
        out.push('…');
    }
    out
}
/// XML 1.0-safe text/attribute escaping. Invalid control characters become U+FFFD.
fn escape(text: &str) -> String {
    let mut out = String::new();
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            '\t' | '\n' | '\r' => out.push(c),
            c if matches!(c as u32,0x20..=0xd7ff|0xe000..=0xfffd|0x10000..=0x10ffff) => out.push(c),
            _ => out.push('�'),
        }
    }
    out
}
struct Flags {
    fuzzy: bool,
    ordered: bool,
    unresolved: bool,
    partial: bool,
}
fn collect(
    location: &Location,
    reverse: bool,
    budget: &mut usize,
    spans: &mut Vec<FeatureSpan>,
    flags: &mut Flags,
) {
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => {
            flags.fuzzy |= *before || *after;
            if *budget == 0 {
                flags.partial = true;
                return;
            }
            *budget -= 1;
            spans.push(FeatureSpan {
                start: *start,
                end: *end,
                reverse,
                fuzzy_start: *before,
                fuzzy_end: *after,
            });
        }
        Location::Complement(inner) => collect(inner, !reverse, budget, spans, flags),
        Location::Join(parts) | Location::Order(parts) => {
            flags.ordered |= matches!(location, Location::Order(_));
            if reverse {
                for part in parts.iter().rev() {
                    collect(part, reverse, budget, spans, flags);
                }
            } else {
                for part in parts {
                    collect(part, reverse, budget, spans, flags);
                }
            }
        }
        Location::Unresolved(_) => flags.unresolved = true,
    }
}
fn overlap(a: &FeatureSpan, b: &(usize, usize)) -> bool {
    a.start < b.1 && b.0 < a.end
}
fn coordinate(position: f64, length: usize, radius: f64) -> (f64, f64) {
    let angle = position / length as f64 * TAU - FRAC_PI_2;
    (500.0 + radius * angle.cos(), 350.0 + radius * angle.sin())
}
fn representative(feature: &MapFeature) -> Option<&FeatureSpan> {
    feature.spans.iter().max_by_key(|s| s.end - s.start)
}
fn label_text(feature: &MapFeature) -> String {
    format!(
        "{} {}{}",
        feature.source_index + 1,
        abbreviated(&feature.label, LABEL_CHARS),
        if feature.needs_review { " *" } else { "" }
    )
}
fn circular_labels(features: &[MapFeature], length: usize) -> Vec<MapLabel> {
    let mut sides = [Vec::new(), Vec::new()];
    for (i, feature) in features.iter().enumerate() {
        let Some(span) = representative(feature) else {
            continue;
        };
        let midpoint = (span.start + span.end) as f64 / 2.0;
        let radius = 205.0 + feature.lane as f64 * 14.0;
        let (x, y) = coordinate(midpoint, length, radius);
        let side = usize::from(x >= 500.0);
        if sides[side].len() < MAX_MAP_LABELS / 2 {
            sides[side].push((i, x, y));
        }
    }
    let mut labels = Vec::new();
    for (side, items) in sides.iter_mut().enumerate() {
        items.sort_by(|a, b| a.2.total_cmp(&b.2));
        let mut previous = 100.0;
        let mut positions = Vec::new();
        for item in items.iter() {
            let y = item.2.clamp(115.0, 595.0).max(previous + 24.0);
            positions.push(y);
            previous = y;
        }
        if positions.last().is_some_and(|y| *y > 595.0) {
            let mut next = 619.0;
            for y in positions.iter_mut().rev() {
                *y = (*y).min(next - 24.0);
                next = *y;
            }
        }
        for ((i, x, y), label_y) in items.iter().zip(positions) {
            labels.push(MapLabel {
                feature_index: *i,
                text: label_text(&features[*i]),
                x: if side == 0 { 230.0 } else { 770.0 },
                y: label_y,
                anchor_x: *x,
                anchor_y: *y,
                right_aligned: side == 0,
            });
        }
    }
    labels
}
fn linear_labels(features: &[MapFeature], length: usize) -> Vec<MapLabel> {
    let mut rows: Vec<Vec<(f64, f64)>> = vec![vec![]; 8];
    let mut labels = Vec::new();
    for (i, feature) in features.iter().enumerate() {
        if labels.len() >= MAX_MAP_LABELS {
            break;
        }
        let Some(span) = representative(feature) else {
            continue;
        };
        let x = 80.0 + (span.start + span.end) as f64 / 2.0 / length as f64 * 840.0;
        let text = label_text(feature);
        let width = (text.chars().count() as f64 * 7.0).clamp(40.0, 230.0);
        let center = x.clamp(20.0 + width / 2.0, 980.0 - width / 2.0);
        let bounds = (center - width / 2.0 - 8.0, center + width / 2.0 + 8.0);
        let Some(row) = rows
            .iter_mut()
            .position(|r| r.iter().all(|a| bounds.0 >= a.1 || a.0 >= bounds.1))
        else {
            continue;
        };
        rows[row].push(bounds);
        labels.push(MapLabel {
            feature_index: i,
            text,
            x: center,
            y: 110.0 + row as f64 * 20.0,
            anchor_x: x,
            anchor_y: 360.0 - feature.lane as f64 * 28.0,
            right_aligned: false,
        });
    }
    labels
}

/// Build bounded map data without reinterpreting feature coordinates or extracting sequence.
/// Unknown locations are flagged; they are never guessed onto the backbone.
pub fn layout_map(molecule: &Molecule) -> Result<MapScene, Error> {
    molecule.validate()?;
    let mut scene = MapScene {
        name: abbreviated(&molecule.name, 128),
        revision: molecule.revision(),
        length: molecule.sequence.len(),
        alphabet: molecule.sequence.alphabet(),
        topology: molecule.topology,
        features: vec![],
        labels: vec![],
        legend: vec![],
        warnings: vec![],
        omitted_features: molecule.features.len().saturating_sub(MAX_MAP_FEATURES),
        spans_omitted: false,
        omitted_labels: 0,
        track_collisions: 0,
        linear_feature_y: 360.0,
    };
    let mut budget = MAX_MAP_SPANS;
    let mut occupied: Vec<Vec<(usize, usize)>> = vec![vec![]; LANES];
    let mut legend = BTreeMap::new();
    for (index, feature) in molecule.features.iter().take(MAX_MAP_FEATURES).enumerate() {
        let mut flags = Flags {
            fuzzy: false,
            ordered: false,
            unresolved: false,
            partial: false,
        };
        let mut spans = Vec::new();
        collect(
            &feature.location,
            false,
            &mut budget,
            &mut spans,
            &mut flags,
        );
        scene.spans_omitted |= flags.partial;
        let (lane, collisions) = (0..LANES)
            .map(|lane| {
                (
                    lane,
                    spans
                        .iter()
                        .map(|s| occupied[lane].iter().filter(|b| overlap(s, b)).count())
                        .sum::<usize>(),
                )
            })
            .min_by_key(|p| p.1)
            .unwrap();
        if collisions > 0 {
            scene.track_collisions += 1;
        }
        occupied[lane].extend(spans.iter().map(|s| (s.start, s.end)));
        let (category, color) = category(&feature.kind);
        *legend.entry((category, color)).or_insert(0usize) += 1;
        let terminal_span = if scene.alphabet != Alphabet::Protein && !flags.partial {
            let terminal = crate::tracks::terminal_span(&feature.location);
            spans
                .iter()
                .rposition(|span| terminal == Some((span.start, span.end, span.reverse)))
        } else {
            None
        };
        scene.features.push(MapFeature {
            source_index: index,
            kind: abbreviated(&feature.kind, 64),
            label: abbreviated(feature.label(), 128),
            color,
            spans,
            lane,
            terminal_span,
            fuzzy: flags.fuzzy,
            ordered: flags.ordered,
            unresolved: flags.unresolved,
            partial: flags.partial,
            needs_review: feature.needs_review,
        });
    }
    scene.legend = legend
        .into_iter()
        .map(|((name, color), count)| LegendEntry { name, color, count })
        .collect();
    scene.labels = if scene.topology == Topology::Circular {
        circular_labels(&scene.features, scene.length)
    } else {
        linear_labels(&scene.features, scene.length)
    };
    if scene.topology == Topology::Linear {
        let last_label = scene.labels.iter().map(|l| l.y).fold(90.0, f64::max);
        let highest_lane = scene.features.iter().map(|f| f.lane).max().unwrap_or(0);
        scene.linear_feature_y = last_label + 50.0 + highest_lane as f64 * 28.0;
        for label in &mut scene.labels {
            label.anchor_y =
                scene.linear_feature_y - scene.features[label.feature_index].lane as f64 * 28.0;
        }
    }
    scene.omitted_labels = scene.features.len().saturating_sub(scene.labels.len());
    if scene.omitted_features > 0 {
        scene.warnings.push(format!(
            "{} additional features are not drawn ({}-feature export limit).",
            scene.omitted_features, MAX_MAP_FEATURES
        ));
    }
    if scene.spans_omitted {
        scene.warnings.push(format!("Additional compound spans are omitted ({}-span export limit); affected features have no direction arrow.",MAX_MAP_SPANS));
    }
    let count =
        |predicate: fn(&MapFeature) -> bool| scene.features.iter().filter(|f| predicate(f)).count();
    let unresolved = count(|f| f.unresolved);
    let fuzzy = count(|f| f.fuzzy);
    let ordered = count(|f| f.ordered);
    let review = count(|f| f.needs_review);
    if unresolved > 0 {
        scene.warnings.push(format!("{unresolved} features contain unresolved locations; unknown portions are not positioned."));
    }
    if fuzzy > 0 || ordered > 0 {
        scene.warnings.push(format!("Approximate boundaries: {fuzzy} features; order(...) relationships: {ordered}. Dashed geometry is indicative; no arrows are drawn."));
    }
    if review > 0 {
        scene.warnings.push(format!("{review} annotations require review after edits; labels carry * and geometry has a red outline."));
    }
    if scene.track_collisions > 0 {
        scene.warnings.push(format!(
            "{} features share occupied tracks; some bands may obscure one another.",
            scene.track_collisions
        ));
    }
    if scene.omitted_labels > 0 {
        scene.warnings.push(format!("{} features have no map label; labels/notes are abbreviated. Inspect the sequence document for complete annotation text.",scene.omitted_labels));
    }
    if scene.alphabet == Alphabet::Protein
        && scene
            .features
            .iter()
            .flat_map(|f| &f.spans)
            .any(|s| s.reverse)
    {
        scene.warnings.push("Protein complement declarations are retained as data, not interpreted as nucleotide strand orientation.".into());
    }
    Ok(scene)
}
fn path_arc(start: f64, end: f64, length: usize, radius: f64) -> String {
    let mut path = String::new();
    let (x, y) = coordinate(start, length, radius);
    let _ = write!(path, "M {x:.3} {y:.3}");
    // Split at half turns so a whole-molecule range remains a complete circle, not a zero arc.
    let count = (((end - start) / length as f64 * 2.0).ceil() as usize).max(1);
    for step in 1..=count {
        let position = start + (end - start) * step as f64 / count as f64;
        let (x, y) = coordinate(position, length, radius);
        let _ = write!(path, " A {radius:.3} {radius:.3} 0 0 1 {x:.3} {y:.3}");
    }
    path
}
fn description(feature: &MapFeature) -> String {
    let mut text = format!(
        "Feature {}: {} — {}",
        feature.source_index + 1,
        feature.kind,
        feature.label
    );
    for (i, span) in feature.spans.iter().take(8).enumerate() {
        let _ = write!(
            text,
            "{} {}{}..{}{} ({})",
            if i == 0 {
                "; one-based displayed ranges:"
            } else {
                ","
            },
            if span.fuzzy_start { "<" } else { "" },
            span.start + 1,
            if span.fuzzy_end { ">" } else { "" },
            span.end,
            if span.reverse { "reverse" } else { "forward" }
        );
    }
    if feature.spans.len() > 8 {
        text.push_str(", additional ranges not listed");
    }
    if feature.ordered {
        text.push_str("; order relationship is not an exact join");
    }
    if feature.unresolved {
        text.push_str("; unresolved portions are not positioned");
    }
    if feature.partial {
        text.push_str("; geometry omitted by export limit");
    }
    if feature.needs_review {
        text.push_str("; annotation review required");
    }
    text
}
fn draw_feature(output: &mut String, scene: &MapScene, feature: &MapFeature) {
    let uncertain = feature.fuzzy || feature.ordered || feature.unresolved || feature.partial;
    let dash = if uncertain {
        " stroke-dasharray=\"6 4\""
    } else {
        ""
    };
    let _ = write!(
        output,
        "<g data-feature-index=\"{}\"><title>{}</title>",
        feature.source_index,
        escape(&description(feature))
    );
    for (i, span) in feature.spans.iter().enumerate() {
        let path = if scene.topology == Topology::Circular {
            path_arc(
                span.start as f64,
                span.end as f64,
                scene.length,
                205.0 + feature.lane as f64 * 14.0,
            )
        } else {
            let y = scene.linear_feature_y - feature.lane as f64 * 28.0;
            format!(
                "M {:.3} {y:.3} L {:.3} {y:.3}",
                80.0 + span.start as f64 / scene.length as f64 * 840.0,
                80.0 + span.end as f64 / scene.length as f64 * 840.0
            )
        };
        if feature.needs_review {
            let _ = write!(
                output,
                "<path d=\"{path}\" fill=\"none\" stroke=\"#b14848\" stroke-width=\"15\"{dash}/>"
            );
        }
        let _ = write!(
            output,
            "<path data-start=\"{}\" data-end=\"{}\" data-reverse=\"{}\" data-fuzzy-start=\"{}\" data-fuzzy-end=\"{}\" d=\"{path}\" fill=\"none\" stroke=\"{}\" stroke-width=\"10\"{dash}/>",
            span.start, span.end, span.reverse, span.fuzzy_start, span.fuzzy_end, feature.color
        );
        if feature.terminal_span == Some(i) {
            let coordinate = if span.reverse { span.start } else { span.end };
            let (tip_x, tip_y, base_x, base_y, normal_x, normal_y) = if scene.topology
                == Topology::Circular
            {
                let radius = 205.0 + feature.lane as f64 * 14.0;
                let angle = coordinate as f64 / scene.length as f64 * TAU - FRAC_PI_2;
                let direction = if span.reverse { -1.0 } else { 1.0 };
                let (x, y) = self::coordinate(coordinate as f64, scene.length, radius);
                let length =
                    ((span.end - span.start) as f64 / scene.length as f64 * TAU * radius * 0.7)
                        .min(12.0);
                (
                    x,
                    y,
                    x + direction * angle.sin() * length,
                    y - direction * angle.cos() * length,
                    angle.cos() * 7.0,
                    angle.sin() * 7.0,
                )
            } else {
                let x = 80.0 + coordinate as f64 / scene.length as f64 * 840.0;
                let y = scene.linear_feature_y - feature.lane as f64 * 28.0;
                let direction = if span.reverse { -1.0 } else { 1.0 };
                let length =
                    ((span.end - span.start) as f64 / scene.length as f64 * 840.0 * 0.7).min(12.0);
                (x, y, x - direction * length, y, 0.0, 7.0)
            };
            let _ = write!(
                output,
                "<path data-terminal=\"{}\" data-direction=\"{}\" d=\"M {tip_x:.3} {tip_y:.3} L {:.3} {:.3} L {:.3} {:.3} Z\" fill=\"{}\"/>",
                coordinate,
                if span.reverse { "reverse" } else { "forward" },
                base_x + normal_x,
                base_y + normal_y,
                base_x - normal_x,
                base_y - normal_y,
                feature.color
            );
        }
    }
    output.push_str("</g>");
}
fn units(alphabet: Alphabet) -> &'static str {
    alphabet.map_unit()
}

/// Render a standalone safe SVG. The same typed scene can be reused by another renderer.
pub fn render_map(molecule: &Molecule) -> Result<String, Error> {
    let scene = layout_map(molecule)?;
    let axis_y = scene.linear_feature_y + 40.0;
    let legend_y = if scene.topology == Topology::Circular {
        654.0
    } else {
        axis_y + 68.0
    };
    let footer_y = legend_y + 66.0;
    let height = footer_y + 108.0 + scene.warnings.len() as f64 * 20.0;
    let name = abbreviated(&scene.name, 40);
    let label = match scene.alphabet {
        Alphabet::Dna => "DNA",
        Alphabet::Rna => "RNA",
        Alphabet::Protein => "protein",
    };
    let topology = if scene.topology == Topology::Circular {
        "circular"
    } else {
        "linear"
    };
    let mut output = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"{height:.0}\" viewBox=\"0 0 1000 {height:.0}\" role=\"img\" data-revision=\"{}\" aria-label=\"{} sequence map\"><title>{}</title><desc>{} {} {topology} {label}; {} features, {} labels. Coordinates in data attributes are zero-based and half-open. {}</desc><rect width=\"1000\" height=\"{height:.0}\" fill=\"white\"/><g font-family=\"system-ui, sans-serif\" fill=\"#233c35\"><text x=\"500\" y=\"42\" font-size=\"22\" font-weight=\"600\" text-anchor=\"middle\">{}</text>",
        scene.revision,
        escape(&name),
        escape(&scene.name),
        scene.length,
        units(scene.alphabet),
        scene.features.len(),
        scene.labels.len(),
        escape(&scene.warnings.join(" ")),
        escape(&name)
    );
    if scene.topology == Topology::Circular {
        output.push_str("<circle cx=\"500\" cy=\"350\" r=\"190\" fill=\"none\" stroke=\"#c7d5ce\" stroke-width=\"2\"/>");
        let _ = write!(
            output,
            "<text x=\"500\" y=\"344\" font-size=\"25\" text-anchor=\"middle\">{} {}</text><text x=\"500\" y=\"373\" font-size=\"15\" text-anchor=\"middle\">Circular {label}</text>",
            scene.length,
            units(scene.alphabet)
        );
        let tick_count = scene.length.min(8);
        for i in 0..tick_count {
            let position = (i * scene.length / tick_count) as f64;
            let (x1, y1) = coordinate(position, scene.length, 181.0);
            let (x2, y2) = coordinate(position, scene.length, 190.0);
            let (x, y) = coordinate(position, scene.length, 164.0);
            let _ = write!(
                output,
                "<path d=\"M {x1:.3} {y1:.3} L {x2:.3} {y2:.3}\" stroke=\"#8ca298\"/><text x=\"{x:.3}\" y=\"{:.3}\" text-anchor=\"middle\" font-size=\"11\" fill=\"#62796e\">{}</text>",
                y + 4.0,
                (position as usize + 1).min(scene.length)
            );
        }
    } else {
        let _ = write!(
            output,
            "<text x=\"500\" y=\"68\" text-anchor=\"middle\" font-size=\"14\">{} {} · linear {label}</text><path d=\"M 80 {axis_y:.3} L 920 {axis_y:.3}\" stroke=\"#c7d5ce\" stroke-width=\"3\"/>",
            scene.length,
            units(scene.alphabet)
        );
        let tick_count = scene.length.min(9);
        for i in 0..tick_count {
            let base = if tick_count == 1 {
                0
            } else {
                i * (scene.length - 1) / (tick_count - 1)
            };
            let x = 80.0 + (base as f64 + 0.5) / scene.length as f64 * 840.0;
            let number = base + 1;
            let _ = write!(
                output,
                "<path d=\"M {x:.3} {axis_y:.3} L {x:.3} {:.3}\" stroke=\"#8ca298\"/><text x=\"{x:.3}\" y=\"{:.3}\" text-anchor=\"middle\" font-size=\"11\" fill=\"#62796e\">{number}</text>",
                axis_y + 9.0,
                axis_y + 28.0
            );
        }
    }
    for feature in &scene.features {
        draw_feature(&mut output, &scene, feature);
    }
    for map_label in &scene.labels {
        let feature = &scene.features[map_label.feature_index];
        let _ = write!(
            output,
            "<path d=\"M {:.3} {:.3} L {:.3} {:.3}\" fill=\"none\" stroke=\"#a4b4ab\" stroke-width=\"0.8\"/>",
            map_label.anchor_x,
            map_label.anchor_y,
            map_label.x,
            map_label.y - 4.0
        );
        if scene.topology == Topology::Circular {
            let x = if map_label.right_aligned { 8.0 } else { 770.0 };
            let text_x = if map_label.right_aligned { 222.0 } else { 0.0 };
            let anchor = if map_label.right_aligned {
                "end"
            } else {
                "start"
            };
            let _ = write!(
                output,
                "<svg x=\"{x:.0}\" y=\"{:.3}\" width=\"222\" height=\"20\" overflow=\"hidden\"><text x=\"{text_x:.0}\" y=\"14\" text-anchor=\"{anchor}\" font-size=\"11.5\" fill=\"{}\"><title>{}</title>{}</text></svg>",
                map_label.y - 14.0,
                feature.color,
                escape(&description(feature)),
                escape(&map_label.text)
            );
        } else {
            let viewport_x = (map_label.x - 115.0).clamp(8.0, 762.0);
            let _ = write!(
                output,
                "<svg x=\"{viewport_x:.3}\" y=\"{:.3}\" width=\"230\" height=\"20\" overflow=\"hidden\"><text x=\"{:.3}\" y=\"14\" font-size=\"11.5\" text-anchor=\"middle\" fill=\"{}\"><title>{}</title>{}</text></svg>",
                map_label.y - 14.0,
                map_label.x - viewport_x,
                feature.color,
                escape(&description(feature)),
                escape(&map_label.text)
            );
        }
    }
    for (i, entry) in scene.legend.iter().enumerate() {
        let x = 24.0 + (i % 4) as f64 * 242.0;
        let y = legend_y + (i / 4) as f64 * 25.0;
        let _ = write!(
            output,
            "<rect x=\"{x:.0}\" y=\"{:.0}\" width=\"10\" height=\"10\" rx=\"2\" fill=\"{}\"/><text x=\"{:.0}\" y=\"{y:.0}\" font-size=\"12\">{} ({})</text>",
            y - 9.0,
            entry.color,
            x + 17.0,
            entry.name,
            entry.count
        );
    }
    let coordinates = if scene.topology == Topology::Circular {
        "Origin (before base 1) is at 12 o'clock; coordinates increase clockwise. Origin-spanning joins retain separate ranges."
    } else {
        "Sequence coordinates increase left to right; visible tick labels are one-based base/residue positions."
    };
    let _ = write!(
        output,
        "<text x=\"24\" y=\"{footer_y:.0}\" font-size=\"11\" fill=\"#62796e\">{coordinates}</text><text x=\"24\" y=\"{:.0}\" font-size=\"11\" fill=\"#62796e\">Arrows denote exact nucleotide feature orientation; protein, fuzzy, order, unresolved or partial features have no arrows.</text><text x=\"24\" y=\"{:.0}\" font-size=\"11\" fill=\"#62796e\">{} of {} features drawn; {} mapped spans; {} labels. Native annotation/evidence data is unchanged.</text>",
        footer_y + 20.0,
        footer_y + 40.0,
        scene.features.len(),
        molecule.features.len(),
        scene.features.iter().map(|f| f.spans.len()).sum::<usize>(),
        scene.labels.len()
    );
    let _ = write!(
        output,
        "<text x=\"24\" y=\"{:.0}\" font-size=\"11\" fill=\"#62796e\">Reference map only: duplex ends, nicks, methylation, secondary structure and experimental feasibility are not drawn.</text>",
        footer_y + 60.0
    );
    let _ = write!(
        output,
        "<text x=\"24\" y=\"{:.0}\" font-size=\"10\" fill=\"#62796e\">Molecule revision: {}</text>",
        footer_y + 82.0,
        scene.revision
    );
    for (i, warning) in scene.warnings.iter().enumerate() {
        let _ = write!(
            output,
            "<text x=\"24\" y=\"{:.0}\" font-size=\"11\" fill=\"#8c5542\">{}</text>",
            footer_y + 102.0 + i as f64 * 20.0,
            escape(warning)
        );
    }
    output.push_str("</g></svg>");
    if output.len() > MAX_SVG_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(output)
}
