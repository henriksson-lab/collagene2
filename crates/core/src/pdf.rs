//! Embedded, deterministic, printable protocol PDF. No filesystem or subprocess use.
//! PDF syntax/font objects use pdf-writer; scientific content is reverified here.
use crate::{
    Error, Molecule, Topology,
    chemistry::DuplexSummary,
    svg::{self, MapScene},
    verification::{self, DesignPlan, Operation, Source, Status},
};
use pdf_writer::types::{CidFontType, FontFlags, SystemInfo, UnicodeCmap};
use pdf_writer::{Content, Finish, Name, Pdf, Rect, Ref, Str, TextStr};
use std::{
    collections::{BTreeMap, BTreeSet},
    f64::consts::{FRAC_PI_2, TAU},
    io,
};
use ttf_parser::Face;

const FONT: &[u8] = include_bytes!("../assets/fonts/LiberationSans-Regular.ttf");
const MARGIN: f32 = 42.0;
const WIDTH: f32 = 595.0;
const HEIGHT: f32 = 842.0;
const TEXT_WIDTH: f32 = WIDTH - MARGIN * 2.0;
const BOTTOM: f32 = 56.0;
const TOP: f32 = 762.0;
pub const MAX_PDF_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PDF_PAGES: usize = 256;
pub const MAX_PDF_TEXT_BYTES: usize = 512 * 1024;
const MAX_INPUT_BYTES: usize = 32 * 1024 * 1024;
const MAX_GLYPHS: usize = 4096;
const GREEN: [f32; 3] = [0.08, 0.31, 0.27];
const INK: [f32; 3] = [0.10, 0.19, 0.17];
const RED: [f32; 3] = [0.65, 0.13, 0.12];

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PdfReport {
    pub bytes: Vec<u8>,
    pub pages: usize,
    pub warnings: Vec<String>,
    pub status: Status,
    pub plan_revision: String,
    pub target_revision: String,
    pub engine_version: String,
}

struct CountBytes(usize);
impl io::Write for CountBytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_INPUT_BYTES)
            .ok_or_else(|| io::Error::other("PDF source limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn preflight(plan: &DesignPlan, target: &Molecule) -> Result<(), Error> {
    if plan.inputs.len() > 1000 || plan.steps.len() > 1000 {
        return Err(Error::TooLarge);
    }
    target.validate()?;
    for input in &plan.inputs {
        input.molecule.validate()?;
    }
    let mut count = CountBytes(0);
    serde_json::to_writer(&mut count, plan).map_err(|_| Error::TooLarge)?;
    serde_json::to_writer(&mut count, target).map_err(|_| Error::TooLarge)?;
    Ok(())
}
fn source_name(source: &Source) -> String {
    match source {
        Source::Input(id) => format!("input {id}"),
        Source::Step(i) => i
            .checked_add(1)
            .map(|n| format!("product of step {n}"))
            .unwrap_or_else(|| format!("invalid step index {i}")),
    }
}
// TOPO request/evidence text shares the PDF text ceiling before collecting a
// String. Large valid graphs fail explicitly rather than truncate a recipe or
// allocate an unbounded JSON value merely to discover the PDF is too large.
struct TopoText(Vec<u8>);
impl io::Write for TopoText {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > MAX_PDF_TEXT_BYTES)
        {
            return Err(io::Error::other("TOPO PDF text limit"));
        }
        self.0
            .try_reserve(bytes.len())
            .map_err(|_| io::Error::other("TOPO PDF text allocation"))?;
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl TopoText {
    fn json(&mut self, value: &impl serde::Serialize) -> Result<(), Error> {
        serde_json::to_writer(self, value).map_err(|_| Error::TooLarge)
    }
    fn finish(self) -> Result<String, Error> {
        String::from_utf8(self.0).map_err(|error| Error::Format(error.to_string()))
    }
}
fn topo_step_text(
    sources: &[verification::TopoSource],
    request: &crate::topo::Request,
    selection: &crate::topo::Selection,
    policy: Option<&crate::topo::AcceptancePolicy>,
    plan: &DesignPlan,
    report: &verification::VerificationReport,
) -> Result<String, Error> {
    use io::Write as _;
    let mut out = TopoText(Vec::new());
    writeln!(out,"TOPO elementary strand-transfer reconstruction under explicitly declared preparation hypotheses. No supplier label or imported reference metadata proves activation, terminal chemistry, kit formulation or experimental success.").map_err(|_|Error::TooLarge)?;
    for source in sources {
        let revision = match &source.source {
            Source::Input(id) => plan
                .inputs
                .iter()
                .find(|input| input.id == *id)
                .map(|input| input.revision.clone()),
            Source::Step(index) => report.products.get(*index).map(Molecule::revision),
        };
        writeln!(
            out,
            "Source role {:?}: {}; revision {}",
            source.role,
            source_name(&source.source),
            revision
                .as_deref()
                .unwrap_or("unavailable; upstream reconstruction is not established")
        )
        .map_err(|_| Error::TooLarge)?;
    }
    // Preserve the complete typed request: strand source intervals/orientations,
    // preparation statements, pairing, transfers and both sets of work limits.
    writeln!(out, "Full source-role recipe and bounded request:").map_err(|_| Error::TooLarge)?;
    out.json(request)?;
    writeln!(out, "\nExplicit component/source-seed selection:").map_err(|_| Error::TooLarge)?;
    out.json(selection)?;
    for (index, material) in request.recipe.materials.iter().enumerate() {
        match material {
            crate::topo::MaterialRecipe::Prepared {
                declaration_id,
                assumption,
                ..
            } => {
                writeln!(
                    out,
                    "\nPrepared material {} / {}: {}",
                    index, declaration_id, assumption
                )
                .map_err(|_| Error::TooLarge)?;
            }
            crate::topo::MaterialRecipe::ExactPcr {
                material_id,
                template_role,
                forward,
                reverse,
                tail,
            } => {
                writeln!(out,"\nExact PCR material {} / {} from role {:?}. Ordered oligos 5′→3′: forward {} = {} (annealing {}; explicit 5′ tail {}); reverse {} = {} (annealing {}; explicit 5′ tail {}). Explicit subsequent 3′ tailing {:?}. These are ideal copied/synthesized strands, not inherited modifications or empirical amplification.",index,material_id,template_role,forward.name,forward.ordered_sequence(),forward.annealing.as_str(),forward.tail,reverse.name,reverse.ordered_sequence(),reverse.annealing.as_str(),reverse.tail,tail).map_err(|_|Error::TooLarge)?;
            }
        }
    }
    writeln!(out,"\nAll reciprocal components, open nicks, retained flaps/free strands and unreleased activations remain in the physical graph. A paired contour is not covalent closure or host repair. Kit conditions, quantities, enzyme recycling, mixture kinetics and laboratory efficacy are unspecified.").map_err(|_|Error::TooLarge)?;
    if let Some(policy) = policy {
        writeln!(out,"Final-only checked TOPO target policy: accept nicks {}; accept unknown terminal chemistry {}; accept unapplied transfers {}; accept other components {}; allow circular-origin rotation {}; allow reverse-complement comparison {}.",policy.allow_nicks,policy.allow_unknown_terminal_chemistry,policy.allow_unapplied_transfers,policy.allow_other_components,policy.allow_circular_rotation,policy.allow_reverse_complement).map_err(|_|Error::TooLarge)?;
        out.json(policy)?;
        writeln!(out,"\nThe final step reconstructs the graph and checks the independent target afresh. Acceptance establishes consistency under these declared assumptions only; it does not repair termini/nicks/flaps, establish experimental confidence or authorize downstream generic chemical transformations. Fully specified graph DNA is required. Annotation correspondence and review evidence do not validate biological function.").map_err(|_|Error::TooLarge)?;
    } else {
        writeln!(out,"Preview-only TOPO operation: no physical-target acceptance policy is supplied; plan verification remains incomplete.").map_err(|_|Error::TooLarge)?;
    }
    out.finish()
}
fn topo_evidence_text(product: &crate::topo::Product) -> Result<String, Error> {
    use io::Write as _;
    let mut out = TopoText(Vec::new());
    let physical = &product.projection.physical;
    writeln!(out,"Fresh TOPO evidence: adapter {} / {}; physical engine {} / model {}. Work complete {}; all requested transfers applied {}; unsupported/unknown {}; fully specified DNA {}.",product.engine_version,product.batch_revision,physical.engine,physical.model,physical.work_complete,physical.all_requested_transfers_applied,physical.unsupported_or_unknown,physical.fully_specified_dna).map_err(|_|Error::TooLarge)?;
    for assumption in product
        .preparation
        .assumptions
        .iter()
        .chain(physical.assumptions.iter())
    {
        writeln!(out, "Preparation/model assumption: {assumption}").map_err(|_| Error::TooLarge)?;
    }
    writeln!(
        out,
        "Actual transfer outcomes, source-call pair changes and displacement:"
    )
    .map_err(|_| Error::TooLarge)?;
    out.json(&physical.steps)?;
    writeln!(out, "\nEvery reconstructed component:").map_err(|_| Error::TooLarge)?;
    out.json(&physical.product.components)?;
    for (index, chain) in physical.product.chains.iter().enumerate() {
        writeln!(out,"\nChemical chain {index}: {} stored 5′→3′ bases; actual ends and immutable source placements:",chain.sequence.len()).map_err(|_|Error::TooLarge)?;
        out.json(&chain.ends)?;
        out.json(&chain.sources)?;
    }
    writeln!(
        out,
        "\nOpen nicks with actual source calls and terminal chemistry:"
    )
    .map_err(|_| Error::TooLarge)?;
    out.json(&physical.product.nicks)?;
    writeln!(
        out,
        "\nRetained unpaired runs, including flaps and free strands:"
    )
    .map_err(|_| Error::TooLarge)?;
    out.json(&physical.product.unpaired)?;
    writeln!(out, "\nExplicitly released activation identities:").map_err(|_| Error::TooLarge)?;
    out.json(&physical.product.released_activations)?;
    match &product.projection.outcome {
        crate::topo::physical::projection::Outcome::Proven(view) => {
            writeln!(out,"\nSelected paired contour {:?}: top {} bases, actual bottom {} bases; known terminal chemistry {}; all requested transfers applied {}. External ends and every retained nick follow:",view.topology,view.top.len(),view.bottom.len(),view.fully_known_terminal_chemistry,view.all_requested_transfers_applied).map_err(|_|Error::TooLarge)?;
            out.json(&view.ends)?;
            out.json(&view.nicks)?;
        }
        crate::topo::physical::projection::Outcome::Withheld(reason) => {
            writeln!(
                out,
                "\nPaired contour withheld; no replacement reference DNA or repair is fabricated:"
            )
            .map_err(|_| Error::TooLarge)?;
            out.json(reason)?;
        }
    }
    writeln!(out,"\nReference preview present: {}. This summary cannot restore physical authority from saved/imported producer metadata or establish experimental target acceptance.",product.reference.is_some()).map_err(|_|Error::TooLarge)?;
    out.finish()
}
fn status_text(status: Status) -> &'static str {
    match status {
        Status::Invalid => "DRAFT - INVALID: do not use as a validated protocol",
        Status::Incomplete => "DRAFT - INCOMPLETE: scientific checks are unresolved",
        Status::ConsistentUnderAssumptions => {
            "Consistent under assumptions - laboratory review required"
        }
    }
}

struct Layout {
    face: Face<'static>,
    pages: Vec<Content>,
    current: Content,
    y: f32,
    glyphs: BTreeMap<char, u16>,
    escaped: BTreeSet<char>,
    text_bytes: usize,
    status: Status,
}
impl Layout {
    fn new(status: Status) -> Result<Self, Error> {
        let mut layout = Self {
            face: Face::parse(FONT, 0)
                .map_err(|_| Error::Format("bundled PDF font is invalid".into()))?,
            pages: vec![],
            current: Content::new(),
            y: TOP,
            glyphs: BTreeMap::new(),
            escaped: BTreeSet::new(),
            text_bytes: 0,
            status,
        };
        layout.header()?;
        Ok(layout)
    }
    fn header(&mut self) -> Result<(), Error> {
        let draft = self.status != Status::ConsistentUnderAssumptions;
        let color = if draft { RED } else { GREEN };
        self.current
            .set_fill_rgb(0.95, 0.97, 0.95)
            .rect(MARGIN, 779.0, TEXT_WIDTH, 30.0)
            .fill_nonzero();
        self.text_at(MARGIN + 8.0, 790.0, status_text(self.status), 10.0, color)?;
        self.text_at(MARGIN, 822.0, "COLLAGENE 2 / CLONING PROTOCOL", 9.0, GREEN)?;
        Ok(())
    }
    fn next_page(&mut self) -> Result<(), Error> {
        if self.pages.len() + 1 >= MAX_PDF_PAGES {
            return Err(Error::TooLarge);
        }
        self.pages
            .push(std::mem::replace(&mut self.current, Content::new()));
        self.y = TOP;
        self.header()
    }
    fn ensure(&mut self, height: f32) -> Result<(), Error> {
        if self.y - height < BOTTOM {
            self.next_page()?;
        }
        Ok(())
    }
    fn width(&self, text: &str, size: f32) -> f32 {
        text.chars()
            .map(|c| {
                self.face
                    .glyph_index(c)
                    .and_then(|g| self.face.glyph_hor_advance(g))
                    .unwrap_or(0) as f32
            })
            .sum::<f32>()
            * size
            / self.face.units_per_em() as f32
    }
    fn normalize(&mut self, text: &str) -> Result<String, Error> {
        self.text_bytes = self
            .text_bytes
            .checked_add(text.len())
            .filter(|n| *n <= MAX_PDF_TEXT_BYTES)
            .ok_or(Error::TooLarge)?;
        let mut normalized = String::new();
        for c in text.chars() {
            // No bidi/complex-script shaping is claimed. Render unsupported scalars visibly,
            // rather than reorder names incorrectly or silently substitute an empty glyph.
            let simple = matches!(c as u32, 0x20..=0x2ff | 0x370..=0x52f | 0x1e00..=0x1eff | 0x2000..=0x206f | 0x2100..=0x22ff);
            let bidi = matches!(c as u32, 0x200b..=0x200f | 0x202a..=0x202e | 0x2060..=0x206f);
            if c == '\n' {
                normalized.push(c);
            } else if c == '\t' {
                normalized.push(' ');
            } else if simple && !bidi && self.face.glyph_index(c).is_some() {
                normalized.push(c);
            } else {
                self.escaped.insert(c);
                normalized.push_str(&format!("[U+{:04X}]", c as u32));
            }
            if normalized.len() > MAX_PDF_TEXT_BYTES {
                return Err(Error::TooLarge);
            }
        }
        Ok(normalized)
    }
    fn text_at(
        &mut self,
        x: f32,
        y: f32,
        text: &str,
        size: f32,
        color: [f32; 3],
    ) -> Result<(), Error> {
        let text = self.normalize(text)?;
        if text.contains('\n') {
            return Err(Error::Format("internal PDF line contains newline".into()));
        }
        let mut encoded = Vec::with_capacity(text.len() * 2);
        for c in text.chars() {
            let next = self.glyphs.len() + 1;
            if next > MAX_GLYPHS && !self.glyphs.contains_key(&c) {
                return Err(Error::TooLarge);
            }
            let cid = *self.glyphs.entry(c).or_insert(next as u16);
            encoded.extend_from_slice(&cid.to_be_bytes());
        }
        self.current.set_fill_rgb(color[0], color[1], color[2]);
        self.current
            .begin_text()
            .set_font(Name(b"F1"), size)
            .next_line(x, y)
            .show(Str(&encoded))
            .end_text();
        Ok(())
    }
    fn paragraph(&mut self, text: &str) -> Result<(), Error> {
        self.paragraph_size(text, 10.5, INK)
    }
    fn paragraph_size(&mut self, text: &str, size: f32, color: [f32; 3]) -> Result<(), Error> {
        let text = self.normalize(text)?;
        let line_height = size * 1.4;
        for paragraph in text.split('\n') {
            let mut line = String::new();
            for word in paragraph.split_whitespace() {
                let separator = if line.is_empty() { "" } else { " " };
                if self.width(&format!("{line}{separator}{word}"), size) <= TEXT_WIDTH {
                    line.push_str(separator);
                    line.push_str(word);
                    continue;
                }
                if !line.is_empty() {
                    self.ensure(line_height)?;
                    let y = self.y;
                    self.text_at(MARGIN, y, &line, size, color)?;
                    self.y -= line_height;
                    line.clear();
                }
                // Hashes, long oligos and source URLs may have no whitespace.
                for c in word.chars() {
                    let mut candidate = line.clone();
                    candidate.push(c);
                    if self.width(&candidate, size) > TEXT_WIDTH && !line.is_empty() {
                        self.ensure(line_height)?;
                        let y = self.y;
                        self.text_at(MARGIN, y, &line, size, color)?;
                        self.y -= line_height;
                        line.clear();
                    }
                    line.push(c);
                }
            }
            self.ensure(line_height)?;
            let y = self.y;
            self.text_at(MARGIN, y, &line, size, color)?;
            self.y -= line_height;
        }
        self.y -= 6.0;
        Ok(())
    }
    fn heading(&mut self, text: &str) -> Result<(), Error> {
        self.ensure(48.0)?;
        self.y -= 8.0;
        self.paragraph_size(text, 14.0, GREEN)
    }
    fn circle_arc(&mut self, center: (f32, f32), radius: f32, start: f64, end: f64) {
        let count = ((end - start).abs() / std::f64::consts::FRAC_PI_2)
            .ceil()
            .max(1.0) as usize;
        let delta = (end - start) / count as f64;
        let point = |a: f64| {
            (
                center.0 + radius * a.cos() as f32,
                center.1 + radius * a.sin() as f32,
            )
        };
        let (x, y) = point(start);
        self.current.move_to(x, y);
        for i in 0..count {
            let a = start + i as f64 * delta;
            let b = a + delta;
            let k = (4.0 / 3.0 * (delta / 4.0).tan()) as f32 * radius;
            let p = point(a);
            let q = point(b);
            self.current.cubic_to(
                p.0 - k * a.sin() as f32,
                p.1 + k * a.cos() as f32,
                q.0 + k * b.sin() as f32,
                q.1 - k * b.cos() as f32,
                q.0,
                q.1,
            );
        }
        self.current.stroke();
    }
    fn map(&mut self, scene: &MapScene) -> Result<(), Error> {
        self.next_page()?;
        self.heading("Intended target map (reference, not a proof of construction)")?;
        self.paragraph(&format!(
            "{} | {} symbols | {:?} {:?}",
            scene.name, scene.length, scene.topology, scene.alphabet
        ))?;
        // The scene has the same positions/lanes/strand terminal span as the SVG exporter.
        // PDF coordinates are bottom-up; the source scene is top-down.
        let scale = TEXT_WIDTH / 1000.0;
        let top = self.y - 4.0;
        let px = |x: f64| MARGIN + x as f32 * scale;
        let py = |y: f64| top - y as f32 * scale;
        let coordinate = |position: usize, radius: f64| {
            let a = position as f64 / scene.length as f64 * TAU - FRAC_PI_2;
            (px(500.0 + radius * a.cos()), py(350.0 + radius * a.sin()))
        };
        self.current
            .set_stroke_rgb(0.73, 0.80, 0.76)
            .set_line_width(1.2);
        if scene.topology == Topology::Circular {
            self.circle_arc((px(500.0), py(350.0)), 190.0 * scale, 0.0, TAU);
            let center = format!("{} symbols", scene.length);
            let w = self.width(&center, 12.0);
            self.text_at(px(500.0) - w / 2.0, py(350.0), &center, 12.0, INK)?;
        } else {
            let y = py(scene.linear_feature_y + 40.0);
            self.current
                .move_to(px(80.0), y)
                .line_to(px(920.0), y)
                .stroke();
            self.text_at(px(80.0), y - 15.0, "1", 8.0, INK)?;
            self.text_at(px(900.0), y - 15.0, &scene.length.to_string(), 8.0, INK)?;
        }
        for feature in &scene.features {
            let color = color(feature.color);
            let uncertain =
                feature.fuzzy || feature.ordered || feature.unresolved || feature.partial;
            for (i, span) in feature.spans.iter().enumerate() {
                if uncertain {
                    self.current.set_dash_pattern([3.0, 2.0], 0.0);
                } else {
                    self.current.set_dash_pattern([], 0.0);
                }
                let passes: Vec<([f32; 3], f32)> = if feature.needs_review {
                    vec![(RED, 7.5), (color, 5.0)]
                } else {
                    vec![(color, 5.0)]
                };
                for (stroke, width) in passes {
                    self.current
                        .set_stroke_rgb(stroke[0], stroke[1], stroke[2])
                        .set_line_width(width);
                    if scene.topology == Topology::Circular {
                        self.circle_arc(
                            (px(500.0), py(350.0)),
                            (205.0 + feature.lane as f32 * 14.0) * scale,
                            FRAC_PI_2 - span.start as f64 / scene.length as f64 * TAU,
                            FRAC_PI_2 - span.end as f64 / scene.length as f64 * TAU,
                        );
                    } else {
                        let y = py(scene.linear_feature_y - feature.lane as f64 * 28.0);
                        self.current
                            .move_to(
                                px(80.0 + span.start as f64 / scene.length as f64 * 840.0),
                                y,
                            )
                            .line_to(px(80.0 + span.end as f64 / scene.length as f64 * 840.0), y)
                            .stroke();
                    }
                }
                if feature.terminal_span == Some(i) {
                    let position = if span.reverse { span.start } else { span.end };
                    let (tip, base, normal) = if scene.topology == Topology::Circular {
                        let radius = 205.0 + feature.lane as f64 * 14.0;
                        let tip = coordinate(position, radius);
                        let angle = FRAC_PI_2 - position as f64 / scene.length as f64 * TAU;
                        let direction = if span.reverse { -1.0 } else { 1.0 };
                        let len = (((span.end - span.start) as f64 / scene.length as f64
                            * TAU
                            * radius
                            * 0.7)
                            .min(12.0)) as f32
                            * scale;
                        (
                            tip,
                            (
                                tip.0 - direction * angle.sin() as f32 * len,
                                tip.1 + direction * angle.cos() as f32 * len,
                            ),
                            (angle.cos() as f32 * 3.5, angle.sin() as f32 * 3.5),
                        )
                    } else {
                        let tip = (
                            px(80.0 + position as f64 / scene.length as f64 * 840.0),
                            py(scene.linear_feature_y - feature.lane as f64 * 28.0),
                        );
                        let direction = if span.reverse { -1.0 } else { 1.0 };
                        let len =
                            (((span.end - span.start) as f64 / scene.length as f64 * 840.0 * 0.7)
                                .min(12.0)) as f32
                                * scale;
                        (tip, (tip.0 - direction * len, tip.1), (0.0, 3.5))
                    };
                    self.current
                        .set_fill_rgb(color[0], color[1], color[2])
                        .move_to(tip.0, tip.1)
                        .line_to(base.0 + normal.0, base.1 + normal.1)
                        .line_to(base.0 - normal.0, base.1 - normal.1)
                        .close_path()
                        .fill_nonzero();
                }
            }
        }
        self.current.set_dash_pattern([], 0.0).set_line_width(0.5);
        for label in &scene.labels {
            let feature = &scene.features[label.feature_index];
            self.current
                .set_stroke_rgb(0.62, 0.70, 0.66)
                .move_to(px(label.anchor_x), py(label.anchor_y))
                .line_to(px(label.x), py(label.y - 4.0))
                .stroke();
            let label_text = self.normalize(&label.text)?;
            let size = 6.3;
            let width = self.width(&label_text, size);
            // Text is fitted in the same bounded label slots; full names remain in documents.
            let max_width = if scene.topology == Topology::Circular {
                222.0 * scale
            } else {
                230.0 * scale
            };
            let size = if width > max_width {
                size * max_width / width
            } else {
                size
            };
            let width = self.width(&label_text, size);
            let x = if scene.topology == Topology::Circular {
                if label.right_aligned {
                    px(230.0) - width
                } else {
                    px(770.0)
                }
            } else {
                px(label.x) - width / 2.0
            };
            self.text_at(x, py(label.y), &label_text, size, color(feature.color))?;
        }
        let map_bottom = if scene.topology == Topology::Circular {
            640.0
        } else {
            scene.linear_feature_y + 100.0
        };
        self.y = py(map_bottom);
        let legend = scene
            .legend
            .iter()
            .map(|e| format!("{} ({})", e.name, e.count))
            .collect::<Vec<_>>()
            .join(" | ");
        self.paragraph_size(&legend, 9.0, INK)?;
        self.paragraph_size("Dashed: uncertain/order/partial geometry. Red outline: review required. Arrow: biological terminal strand direction where exact. Labels may be abbreviated; full annotations remain in the sequence document.",9.0,INK)?;
        for warning in &scene.warnings {
            self.paragraph_size(warning, 9.0, RED)?;
        }
        Ok(())
    }
    fn finish(
        mut self,
        title: &str,
        plan_revision: &str,
        target_revision: &str,
    ) -> Result<(Vec<u8>, usize, Vec<String>), Error> {
        self.pages
            .push(std::mem::replace(&mut self.current, Content::new()));
        let count = self.pages.len();
        // Add final page count without disturbing the scientific content layout.
        for index in 0..count {
            let page = std::mem::replace(&mut self.pages[index], Content::new());
            self.current = page;
            self.text_at(
                MARGIN,
                30.0,
                &format!(
                    "Page {} of {} | Plan {} | Target {}",
                    index + 1,
                    count,
                    &plan_revision[..12],
                    &target_revision[..12]
                ),
                8.0,
                INK,
            )?;
            self.pages[index] = std::mem::replace(&mut self.current, Content::new());
        }
        let mut pdf = Pdf::new();
        let catalog = Ref::new(1);
        let pages = Ref::new(2);
        let font = Ref::new(3);
        let cid = Ref::new(4);
        let descriptor = Ref::new(5);
        let font_file = Ref::new(6);
        let cmap_id = Ref::new(7);
        let gid_map_id = Ref::new(8);
        let info = Ref::new(9);
        pdf.catalog(catalog).pages(pages).lang(TextStr("en"));
        pdf.document_info(info)
            .title(TextStr(title))
            .subject(TextStr(status_text(self.status)))
            .creator(TextStr("Collagene 2"))
            .producer(TextStr("Collagene embedded Rust protocol exporter"));
        let page_ids: Vec<_> = (0..count).map(|i| Ref::new(10 + i as i32 * 2)).collect();
        pdf.pages(pages)
            .kids(page_ids.iter().copied())
            .count(count as i32);
        let mut page_bytes = 0usize;
        for (i, content) in self.pages.into_iter().enumerate() {
            let content_id = Ref::new(11 + i as i32 * 2);
            let mut page = pdf.page(page_ids[i]);
            page.parent(pages)
                .media_box(Rect::new(0.0, 0.0, WIDTH, HEIGHT))
                .contents(content_id);
            page.resources().fonts().pair(Name(b"F1"), font);
            page.finish();
            let bytes = content.finish();
            page_bytes = page_bytes
                .checked_add(bytes.len())
                .filter(|n| *n <= MAX_PDF_BYTES - FONT.len())
                .ok_or(Error::TooLarge)?;
            pdf.stream(content_id, &bytes);
        }
        pdf.type0_font(font)
            .base_font(Name(b"LiberationSans"))
            .encoding_predefined(Name(b"Identity-H"))
            .descendant_font(cid)
            .to_unicode(cmap_id);
        let system = SystemInfo {
            registry: Str(b"Adobe"),
            ordering: Str(b"Identity"),
            supplement: 0,
        };
        let mut cid_font = pdf.cid_font(cid);
        cid_font
            .subtype(CidFontType::Type2)
            .base_font(Name(b"LiberationSans"))
            .system_info(system)
            .font_descriptor(descriptor)
            .cid_to_gid_map_stream(gid_map_id);
        let mut widths = vec![0.0; self.glyphs.len() + 1];
        let mut gid_map = vec![0u8; (self.glyphs.len() + 1) * 2];
        let scale = 1000.0 / self.face.units_per_em() as f32;
        let mut cmap = UnicodeCmap::new(Name(b"CollageneUnicode"), system);
        for (&character, &code) in &self.glyphs {
            let glyph = self
                .face
                .glyph_index(character)
                .ok_or_else(|| Error::Format("internal unsupported PDF glyph".into()))?;
            widths[code as usize] = self.face.glyph_hor_advance(glyph).unwrap_or(0) as f32 * scale;
            gid_map[code as usize * 2..code as usize * 2 + 2]
                .copy_from_slice(&glyph.0.to_be_bytes());
            cmap.pair(code, character);
        }
        cid_font.widths().consecutive(0, widths);
        cid_font.finish();
        let bbox = self.face.global_bounding_box();
        pdf.font_descriptor(descriptor)
            .name(Name(b"LiberationSans"))
            .flags(FontFlags::SYMBOLIC)
            .bbox(Rect::new(
                bbox.x_min as f32 * scale,
                bbox.y_min as f32 * scale,
                bbox.x_max as f32 * scale,
                bbox.y_max as f32 * scale,
            ))
            .italic_angle(0.0)
            .ascent(self.face.ascender() as f32 * scale)
            .descent(self.face.descender() as f32 * scale)
            .cap_height(self.face.capital_height().unwrap_or(self.face.ascender()) as f32 * scale)
            .stem_v(80.0)
            .font_file2(font_file);
        pdf.stream(font_file, FONT)
            .pair(Name(b"Length1"), FONT.len() as i32);
        pdf.stream(gid_map_id, &gid_map);
        let cmap_bytes = cmap.finish();
        let cmap_text = std::str::from_utf8(&cmap_bytes)
            .map_err(|_| Error::Format("internal PDF CMap encoding".into()))?
            .replace("/CMapType 0 def", "/CMapType 2 def");
        pdf.cmap(cmap_id, cmap_text.as_bytes());
        let bytes = pdf.finish();
        if bytes.len() > MAX_PDF_BYTES {
            return Err(Error::TooLarge);
        }
        let warnings = if self.escaped.is_empty() {
            vec![]
        } else {
            vec![format!(
                "PDF font/layout cannot represent {} distinct scalars without unsupported shaping or missing glyphs. They are displayed as [U+XXXX] escapes; original Unicode remains in source documents and title metadata. Examples: {}",
                self.escaped.len(),
                self.escaped
                    .iter()
                    .take(16)
                    .map(|c| format!("U+{:04X}", *c as u32))
                    .collect::<Vec<_>>()
                    .join(", ")
            )]
        };
        Ok((bytes, count, warnings))
    }
}
fn color(hex: &str) -> [f32; 3] {
    let component =
        |start| u8::from_str_radix(&hex[start..start + 2], 16).unwrap_or(0) as f32 / 255.0;
    [component(1), component(3), component(5)]
}

/// Reconstructs the plan against this exact target at export time. An externally
/// asserted report/status cannot be supplied. Invalid/incomplete plans export as
/// conspicuously marked drafts, including successfully reconstructed intermediates.
pub fn export_protocol_pdf(plan: &DesignPlan, target: &Molecule) -> Result<PdfReport, Error> {
    preflight(plan, target)?;
    let report = verification::verify(plan, target);
    let mut layout = Layout::new(report.status)?;
    layout.paragraph_size(&plan.name, 21.0, GREEN)?;
    layout.paragraph("Computed sequence/chemistry checks. Reaction conditions, quantities and experimental outcome require laboratory review.")?;
    layout.paragraph(&format!(
        "Target: {} | {} symbols | {:?} {:?}",
        target.name,
        target.sequence.len(),
        target.topology,
        target.sequence.alphabet()
    ))?;
    layout.heading("Starting materials")?;
    for input in &plan.inputs {
        layout.paragraph(&format!(
            "Input {}: {} | {} symbols | {:?}\nRevision {}",
            input.id,
            input.molecule.name,
            input.molecule.sequence.len(),
            input.molecule.topology,
            input.revision
        ))?;
        for key in [
            "source:provider",
            "source:canonical_url",
            "source:identifier",
            "source:accession_version",
            "source:retrieved_at",
            "source:content_blake3",
        ] {
            if let Some(value) = input.molecule.metadata.get(key) {
                layout.paragraph_size(&format!("{key}: {value}"), 9.0, INK)?;
            }
        }
    }
    layout.heading("Cloning steps")?;
    for (index, step) in plan.steps.iter().enumerate() {
        let body = match step {
            Operation::DeclareSourceTermini {
                template,
                left,
                right,
            } => format!(
                "Declare original intact linear duplex termini for {}: left 5prime {:?} / 3prime {:?}; right 5prime {:?} / 3prime {:?}. These are explicit source-material assumptions, not measured chemistry, imported proof, synthesis or repair of a physical intermediate.",
                source_name(template),
                left.five_prime,
                left.three_prime,
                right.five_prime,
                right.three_prime
            ),
            Operation::Topo {
                sources,
                request,
                selection,
            } => topo_step_text(sources, request, selection, None, plan, &report)?,
            Operation::CheckedTopo {
                sources,
                request,
                selection,
                policy,
            } => topo_step_text(sources, request, selection, Some(policy), plan, &report)?,
            Operation::RecombineReference {
                sources,
                request,
                selection,
            } => format!(
                "Reference-only {:?} exchange of {} and {}. Pinned full profile, data, work and annotation limits: {}. Explicit selected candidate/component: {}. All reciprocal components are reconstructed; reference sequence matching does not establish physical chemistry, biological repair, efficacy or experimental selection.",
                request.reaction,
                source_name(&sources[0]),
                source_name(&sources[1]),
                serde_json::to_string(request).map_err(|e| Error::Format(e.to_string()))?,
                serde_json::to_string(selection).map_err(|e| Error::Format(e.to_string()))?
            ),
            Operation::DeclareSourceModifications {
                template,
                assumptions,
            } => format!(
                "Declare source-generation modification assumptions for {}: Dam {:?}, Dcm {:?}, CpG {:?}. This is an explicit three-mechanism generation model; it is not measured methylation or an imported physical-state claim.",
                source_name(template),
                assumptions.dam,
                assumptions.dcm,
                assumptions.cpg
            ),
            Operation::CheckedDigest {
                template,
                enzymes,
                profile_revision,
                fragment_index,
            } => format!(
                "Checked digest of {} with {}; select fragment index {} (zero based) under positional sensitivity profile {}. Reconstructed tracked strand modifications must permit cleavage; unknown, mixed or unsupported effects cannot silently become an unchecked digest.",
                source_name(template),
                enzymes.join(", "),
                fragment_index,
                profile_revision
            ),
            Operation::PcrWithMismatches {
                template,
                forward,
                reverse,
                options,
            } => format!(
                "Intentional-mismatch PCR of {}. Oligos 5′→3′: {} (annealing {}; tail {}), {} (annealing {}; tail {}). Recorded binding/search policy: {}. Primer-derived changes are reconstructed; empirical amplification is not established.",
                source_name(template),
                forward.ordered_sequence(),
                forward.annealing.as_str(),
                forward.tail,
                reverse.ordered_sequence(),
                reverse.annealing.as_str(),
                reverse.tail,
                serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?
            ),
            Operation::GoldenGate {
                templates,
                options,
                candidate_index,
            } => format!(
                "Golden Gate with {:?}; inputs {}. Explicit zero-based candidate index {} under search policy {}. Complete cleavage, orientation, cohesive ends and recognition-site removal are checked. One-copy fragment search does not establish exhaustive mixture uniqueness; reaction conditions and quantities are unspecified.",
                options.enzyme,
                templates
                    .iter()
                    .map(source_name)
                    .collect::<Vec<_>>()
                    .join(", "),
                candidate_index,
                serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?
            ),
            Operation::Pcr {
                template,
                forward,
                reverse,
            } => format!(
                "PCR of {}. Primers 5′→3′:\n{}: {} (annealing {}; 5′ tail {})\n{}: {} (annealing {}; 5′ tail {}).\nReaction conditions, enzyme choice and quantities are unspecified.",
                source_name(template),
                forward.name,
                forward.ordered_sequence(),
                forward.annealing.as_str(),
                forward.tail,
                reverse.name,
                reverse.ordered_sequence(),
                reverse.annealing.as_str(),
                reverse.tail
            ),
            Operation::CatalogDigest {
                template,
                selection,
                fragment_index,
            } => format!(
                "Catalog-pinned complete digest of {}. Exact ordered identity/name/geometry selection: {}. Select zero-based fragment {} in top-reference cut ordering. The embedded catalog pin is checked; complete cleavage/conditions remain assumptions. Supplier documentary observations do not prove methylation compatibility, reaction kinetics or experimental completion.",
                source_name(template),
                serde_json::to_string(selection).unwrap_or_default(),
                fragment_index
            ),
            Operation::Digest {
                template,
                enzymes,
                fragment_index,
            } => format!(
                "Complete restriction digest of {} with {}. Select zero-based fragment index {} in ascending top-reference cut order. Cleavage assumes compatible conditions and no blocking methylation.",
                source_name(template),
                enzymes.join(", "),
                fragment_index
            ),
            Operation::TreatEnds {
                template,
                treatment,
            } => format!(
                "Terminal {treatment:?} treatment of {}. Internal nick repair is not assumed; reaction conditions and quantities are unspecified.",
                source_name(template)
            ),
            Operation::Tail3Prime { template, base } => format!(
                "Append one {} nucleotide to BOTH strand3′ ends of {}. Requires reconstructed blunt linear DNA and accessible3′-OH. Retain5′ chemistry and internal nicks. No template ancestry, enzyme kinetics or commercial kit activity is inferred.",
                base.as_str(),
                source_name(template)
            ),
            Operation::ProfiledOverlap {
                fragments,
                overlaps,
                circular,
                options,
            } => format!(
                "{:?} profiled ordered overlap assembly: {}. Homology lengths {:?}bp; {} product. Explicit orientations are checked. Vendor recommendations are advisory; reagent compatibility, completion and uniqueness remain unestablished. Options {}.",
                options.method,
                fragments
                    .iter()
                    .map(|f| format!(
                        "{} ({})",
                        source_name(&f.source),
                        if f.reverse {
                            "reverse complement"
                        } else {
                            "forward"
                        }
                    ))
                    .collect::<Vec<_>>()
                    .join(" -> "),
                overlaps,
                if *circular { "circular" } else { "linear" },
                serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?
            ),
            Operation::RepairEnds { template, method } => format!(
                "End repair {method:?} of {}. Reconstruct explicit 3′ fill-in/exonuclease changes. Five-prime chemistry and internal nicks remain unchanged; conditions, reagents and empirical success are not established.",
                source_name(template)
            ),
            Operation::AnnealOligos {
                top,
                bottom,
                options,
                candidate_index,
            } => format!(
                "Anneal declared 5′→3′ ssDNA inputs {} and {}.\nTop ordered oligo 5′→3′: {} (5′ {:?}; 3′ {:?}).\nBottom ordered oligo 5′→3′: {} (5′ {:?}; 3′ {:?}).\nSelect zero-based candidate {} under policy {}. Exact antiparallel pairing and placement alternatives are reconstructed. Terminal chemistry is an explicit assumption; annealing performs no ligation or implicit phosphorylation.",
                source_name(top),
                source_name(bottom),
                verification::protocol_oligo_sequence(plan, &report, top, index)
                    .unwrap_or("Unresolved or unsupported bound oligo sequence"),
                options.top_chemistry.five_prime,
                options.top_chemistry.three_prime,
                verification::protocol_oligo_sequence(plan, &report, bottom, index)
                    .unwrap_or("Unresolved or unsupported bound oligo sequence"),
                options.bottom_chemistry.five_prime,
                options.bottom_chemistry.three_prime,
                candidate_index,
                serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?
            ),
            Operation::Ligate {
                fragments,
                circular,
                allow_nicked,
            } => format!(
                "Ordered ligation: {}. {} product. Retained nicks {}. Exact end complementarity and terminal chemistry are checked; competing products and reaction conditions are not modeled.",
                fragments
                    .iter()
                    .map(|f| format!(
                        "{} ({})",
                        source_name(&f.source),
                        if f.reverse { "reversed" } else { "forward" }
                    ))
                    .collect::<Vec<_>>()
                    .join(", "),
                if *circular { "Circular" } else { "Linear" },
                if *allow_nicked {
                    "explicitly allowed"
                } else {
                    "not allowed"
                }
            ),
            Operation::Overlap {
                fragments,
                overlaps,
                circular,
            } => format!(
                "Ordered overlap assembly: {}. Junction overlaps {overlaps:?} bp; {} product. Reagent conditions and quantities are unspecified.",
                fragments
                    .iter()
                    .map(source_name)
                    .collect::<Vec<_>>()
                    .join(", "),
                if *circular { "circular" } else { "linear" }
            ),
        };
        layout.paragraph(&format!("Step {}. {body}", index + 1))?;
        if let Some(topo) = report
            .step_reports
            .iter()
            .find(|entry| entry.step_index == index)
            .and_then(|entry| entry.topo.as_ref())
        {
            layout.paragraph_size(&topo_evidence_text(topo)?, 9.0, INK)?;
        }
        if let Some(annotations) = report
            .annotation_reports
            .iter()
            .find(|entry| entry.step_index == index)
        {
            layout.paragraph_size("TOPO annotation correspondence follows actual top-reference source calls. It is not inherited chemical state or experimental validation.",9.0,INK)?;
            for source in &annotations.evidence.sources {
                layout.paragraph_size(&format!("Source role {} (revision {}): {} projected, {} outside this selected contour, {} retained only as original unmapped evidence. Projection supported: {}. Review all disrupted/PCR-derived annotations.",source.role,source.source_revision,source.projected,source.dropped,source.unmapped.len(),source.projection_supported),9.0,INK)?;
                for warning in source.warnings.iter().take(8) {
                    layout.paragraph_size(warning, 9.0, INK)?;
                }
            }
        }
        if let Some(product) = report.products.get(index) {
            layout.paragraph_size(
                &format!(
                    "Reconstructed: {} | {} symbols | {:?}\nProduct revision {}",
                    product.name,
                    product.sequence.len(),
                    product.topology,
                    product.revision()
                ),
                9.0,
                INK,
            )?;
            if let Some(summary) = product
                .metadata
                .get("duplex_model")
                .and_then(|s| serde_json::from_str::<DuplexSummary>(s).ok())
            {
                layout.paragraph_size(&summary.describe(), 9.0, INK)?;
            }
        } else {
            layout.paragraph_size("No reconstructed product for this step. See verification findings; later steps are not established.",9.0,RED)?;
        }
    }
    layout.heading("Verification findings")?;
    for finding in &report.findings {
        layout.paragraph(&format!("- {finding}"))?;
    }
    layout.heading("Assumptions and limits")?;
    for assumption in &report.assumptions {
        layout.paragraph(&format!("- {assumption}"))?;
    }
    layout.paragraph(&format!(
        "Exhaustive alternative-product search: {}",
        report.exhaustive_alternative_search
    ))?;
    layout.heading("Laboratory details to complete before use")?;
    layout.paragraph("Record reaction volumes, material quantities/concentrations, enzymes/reagents, buffer, temperature/time program, purification/selection and experimental confirmation. These are not supplied or validated by the sequence model.")?;
    layout.heading("Provenance")?;
    layout.paragraph(&format!(
        "Engine {}\nPlan {}\nTarget {}",
        report.engine_version, report.plan_revision, report.target_revision
    ))?;
    let mut map_warnings = Vec::new();
    match svg::layout_map(target) {
        Ok(scene) => {
            map_warnings.extend(scene.warnings.iter().cloned());
            layout.map(&scene)?;
        }
        Err(e) => {
            let warning = format!("Target map unavailable: {e}");
            layout.paragraph(&warning)?;
            map_warnings.push(warning);
        }
    }
    if !layout.escaped.is_empty() {
        layout.heading("PDF text rendering limits")?;
        layout.paragraph("Some Unicode scalars need shaping or a different font. They are visibly written as [U+XXXX] escapes, never silently removed. Full Unicode names remain in the original sequence documents and PDF title metadata. This export is not tagged PDF, PDF/A, or a claim of complete script shaping.")?;
    }
    let (bytes, pages, mut warnings) =
        layout.finish(&plan.name, &report.plan_revision, &report.target_revision)?;
    warnings.extend(map_warnings);
    Ok(PdfReport {
        bytes,
        pages,
        warnings,
        status: report.status,
        plan_revision: report.plan_revision,
        target_revision: report.target_revision,
        engine_version: report.engine_version,
    })
}
