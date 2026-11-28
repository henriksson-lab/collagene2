//! Declared ideal 5′-resection / exact annealing / gap-fill / nick-seal model.
//! See docs/HOMOLOGY_MATERIAL_SOURCES.md. This is neither a kinetic simulation
//! nor a reconstruction of a proprietary mix. There is no intended-target input.
use crate::{
    chemistry::{Duplex, FivePrime, Nick, Strand, StrandMapping, Terminus, ThreePrime},
    cloning::overlap::extension::SourcedMapping,
};
use collagene_model::{Error, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const ENGINE_VERSION: &str = "homology-assembly-1";
pub const MAX_SOURCE_BASES: usize = 2_000_000;
pub const MAX_COMPARISONS: usize = 100_000_000;
pub const MAX_REPORT_BYTES: usize = 32 << 20;
pub const MAX_SOURCES: usize = 32;
pub const MAX_CANDIDATES: usize = 1000;
macro_rules! dto {($($t:tt)*) => {#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)] #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))] #[serde(deny_unknown_fields)] $($t)*};}
dto! {pub struct EndCuts {pub left:usize,pub right:usize}}
dto! {pub enum Resection {ExactHomology,Explicit{cuts:Vec<EndCuts>}}}
dto! {pub enum OuterEndPolicy {JoiningEndsOnly,AllAccessibleEnds{left:usize,right:usize}}}
dto! {pub struct Options {
    pub overlaps:Vec<usize>,pub topology:Topology,pub resection:Resection,
    pub outer_end_policy:OuterEndPolicy,pub min_paired_bases:usize,
    pub max_comparisons:usize,pub max_candidates:usize,
}}
impl Default for Options {
    fn default() -> Self {
        Self {
            overlaps: vec![20],
            topology: Topology::Circular,
            resection: Resection::ExactHomology,
            outer_end_policy: OuterEndPolicy::JoiningEndsOnly,
            min_paired_bases: 12,
            max_comparisons: 20_000_000,
            max_candidates: 100,
        }
    }
}
dto! {pub struct ResectedSpan {pub source:usize,pub strand:Strand,pub start:usize,pub end:usize}}
dto! {pub struct Junction {
    pub left_source:usize,pub right_source:usize,pub overlap:usize,pub homology:String,
    pub top_gap:usize,pub bottom_gap:usize,
    /// Product top-reference boundaries, including any linear top-origin shift.
    pub top_nick:i64,pub bottom_nick:i64,
    pub homology_start:i64,pub homology_end:i64,
    pub top_five_prime:FivePrime,pub top_three_prime:ThreePrime,
    pub bottom_five_prime:FivePrime,pub bottom_three_prime:ThreePrime,
    pub top_sealed:bool,pub bottom_sealed:bool,
    /// Other terminal overlap lengths within the declared adjacent-pair search.
    pub alternatives:Vec<usize>,pub terminal_search_complete:bool,
}}
dto! {pub struct Report {
    pub engine_version:String,pub source_geometry_revisions:Vec<String>,pub source_lengths:Vec<usize>,pub options:Options,
    pub reference_length:usize,pub reference_top_start:usize,pub reference_top_end:usize,
    pub product:Option<Duplex>,pub mappings:Vec<SourcedMapping>,pub resected:Vec<ResectedSpan>,pub junctions:Vec<Junction>,
    pub comparisons:usize,pub terminal_search_complete:bool,pub incomplete:bool,pub warnings:Vec<String>,pub limits:Vec<String>,
}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("physical homology assembly: {s}"))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_REPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("homology evidence limit"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(v: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Count(0), v).map_err(|_| Error::TooLarge)
}
fn revision(v: &impl Serialize) -> Result<String, Error> {
    bounded(v)?;
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}
impl Options {
    pub fn validate_shape(&self, source_count: usize) -> Result<(), Error> {
        if !(1..=MAX_SOURCES).contains(&source_count)
            || self.overlaps.len()
                != source_count - usize::from(self.topology != Topology::Circular)
        {
            return Err(bad("source/junction count differs"));
        }
        if let Resection::Explicit { cuts } = &self.resection {
            if cuts.len() != source_count {
                return Err(bad("explicit cuts must name every source"));
            }
        }
        if self.overlaps.len() > MAX_SOURCES
            || self.min_paired_bases == 0
            || self
                .overlaps
                .iter()
                .any(|h| *h < self.min_paired_bases || *h > MAX_SOURCE_BASES / 2)
            || self.max_comparisons == 0
            || self.max_comparisons > MAX_COMPARISONS
            || self.max_candidates == 0
            || self.max_candidates > MAX_CANDIDATES
        {
            return Err(bad("invalid overlap, work, or candidate bounds"));
        }
        if let Resection::Explicit { cuts } = &self.resection {
            if cuts.is_empty()
                || cuts.len() > MAX_SOURCES
                || cuts
                    .iter()
                    .any(|c| c.left > MAX_SOURCE_BASES / 2 || c.right > MAX_SOURCE_BASES / 2)
            {
                return Err(Error::TooLarge);
            }
        }
        if let OuterEndPolicy::AllAccessibleEnds { left, right } = &self.outer_end_policy {
            if *left > MAX_SOURCE_BASES / 2 || *right > MAX_SOURCE_BASES / 2 {
                return Err(Error::TooLarge);
            }
            if self.topology == Topology::Circular {
                return Err(bad(
                    "circular products have no outer ends; select JoiningEndsOnly",
                ));
            }
        }
        Ok(())
    }
}
fn cuts(lengths: &[usize], o: &Options) -> Result<Vec<EndCuts>, Error> {
    let n = lengths.len();
    let circular = o.topology == Topology::Circular;
    if n < 1 || n > MAX_SOURCES || o.overlaps.len() != n - usize::from(!circular) {
        return Err(bad(
            "ordered sources and declared junction count differ (1..32 sources required)",
        ));
    }
    let outer = match o.outer_end_policy {
        OuterEndPolicy::JoiningEndsOnly => (0, 0),
        OuterEndPolicy::AllAccessibleEnds { left, right } => (left, right),
    };
    let mut c = match &o.resection {
        Resection::ExactHomology => (0..n)
            .map(|i| EndCuts {
                left: if i > 0 {
                    o.overlaps[i - 1]
                } else if circular {
                    o.overlaps[n - 1]
                } else {
                    outer.0
                },
                right: if i < n - 1 || circular {
                    o.overlaps[i]
                } else {
                    outer.1
                },
            })
            .collect(),
        Resection::Explicit { cuts } => {
            if cuts.len() != n {
                return Err(bad("explicit cuts must name every oriented source"));
            }
            cuts.clone()
        }
    };
    if !circular && (c[0].left != outer.0 || c[n - 1].right != outer.1) {
        return Err(bad(
            "explicit outer cuts differ from the declared outer-end policy",
        ));
    }
    for i in 0..n {
        if c[i]
            .left
            .checked_add(c[i].right)
            .is_none_or(|k| k >= lengths[i])
        {
            return Err(bad(
                "resection domains cross, exhaust a strand, or leave no original paired core",
            ));
        }
    }
    for (i, h) in o.overlaps.iter().enumerate() {
        let j = (i + 1) % n;
        if *h >= lengths[i].min(lengths[j]) || c[i].right < *h || c[j].left < *h {
            return Err(bad(
                "joining resection must expose the entire declared homology without exhausting an operand",
            ));
        }
    }
    // Keep explicit state construction separate from the ordinary end policy.
    c.shrink_to_fit();
    Ok(c)
}
/// Borrowed size/shape/work check, before any original strands are cloned.
pub fn preflight(sources: &[Duplex], o: &Options) -> Result<(), Error> {
    o.validate_shape(sources.len())?;
    preflight_borrowed(&sources.iter().collect::<Vec<_>>(), o)
}
/// The tracked bridge can check geometry and reserved work before DNA clones.
pub fn preflight_borrowed(sources: &[&Duplex], o: &Options) -> Result<(), Error> {
    o.validate_shape(sources.len())?;
    if sources.is_empty() || sources.len() > MAX_SOURCES {
        return Err(bad("1..32 ordered source duplexes required"));
    }
    let mut total = 0usize;
    let mut lengths = Vec::with_capacity(sources.len());
    for s in sources {
        total = total
            .checked_add(s.top.len())
            .and_then(|n| n.checked_add(s.bottom.len()))
            .filter(|n| *n <= MAX_SOURCE_BASES)
            .ok_or(Error::TooLarge)?;
        if s.topology != Topology::Linear
            || s.bottom_start != 0
            || s.top.len() != s.bottom.len()
            || !s.nicks.is_empty()
        {
            return Err(bad(
                "only intact paired linear source duplexes are supported",
            ));
        }
        lengths.push(s.top.len());
    }
    cuts(&lengths, o)?;
    let product = lengths
        .iter()
        .sum::<usize>()
        .checked_sub(o.overlaps.iter().sum())
        .filter(|n| *n > 0)
        .ok_or_else(|| bad("invalid overlapping product extent"))?;
    // Conservative shared reservation: geometry/source validation, tracked
    // validation and copied-piece checks, and product validation in both layers.
    let required = 3 * (total / 2) + 2 * product + o.overlaps.iter().sum::<usize>();
    if required > o.max_comparisons {
        return Err(bad(
            "comparison budget cannot establish source pairing, declared homologies and product pairing",
        ));
    }
    bounded(o)
}
#[derive(Clone, Copy)]
struct Cell {
    source: u8,
    index: u32,
    synthesized: bool,
}
const EMPTY: Cell = Cell {
    source: u8::MAX,
    index: 0,
    synthesized: false,
};
fn source_bytes(s: &Duplex, strand: Strand) -> &[u8] {
    match strand {
        Strand::Top => s.top.as_str().as_bytes(),
        Strand::Bottom => s.bottom.as_str().as_bytes(),
    }
}
fn cell_base(c: Cell, product_strand: Strand, sources: &[Duplex]) -> u8 {
    let source_strand = if c.synthesized {
        other(product_strand)
    } else {
        product_strand
    };
    let b = source_bytes(&sources[c.source as usize], source_strand)[c.index as usize];
    if c.synthesized { complement(b) } else { b }
}
fn other(s: Strand) -> Strand {
    match s {
        Strand::Top => Strand::Bottom,
        Strand::Bottom => Strand::Top,
    }
}
fn mapping_runs(cells: &[Cell], strand: Strand) -> Vec<SourcedMapping> {
    let mut out = vec![];
    let mut i = 0;
    while i < cells.len() {
        let c = cells[i];
        let mut end = i + 1;
        while end < cells.len() {
            let d = cells[end];
            let expected = if c.synthesized {
                (c.index as usize).checked_sub(end - i)
            } else {
                (c.index as usize).checked_add(end - i)
            };
            if d.source != c.source
                || d.synthesized != c.synthesized
                || expected != Some(d.index as usize)
            {
                break;
            }
            end += 1;
        }
        let first = c.index as usize;
        let last = cells[end - 1].index as usize;
        out.push(SourcedMapping {
            source: c.source as usize,
            mapping: StrandMapping {
                source_strand: if c.synthesized { other(strand) } else { strand },
                source_start: first.min(last),
                source_end: first.max(last) + 1,
                product_strand: strand,
                product_start: i,
                product_end: end,
                reverse: c.synthesized,
                synthesized: c.synthesized,
            },
        });
        i = end;
    }
    out
}
/// Fresh strand reconstruction. Report data alone never establishes material authority.
pub fn geometry(sources: &[Duplex], o: &Options) -> Result<Report, Error> {
    preflight(sources, o)?;
    for s in sources {
        s.validate()?;
    }
    let n = sources.len();
    let lengths: Vec<_> = sources.iter().map(|s| s.top.len()).collect();
    let c = cuts(&lengths, o)?;
    let circular = o.topology == Topology::Circular;
    let length = lengths.iter().sum::<usize>() - o.overlaps.iter().sum::<usize>();
    let mut offsets = vec![0usize; n];
    for i in 1..n {
        offsets[i] = offsets[i - 1] + lengths[i - 1] - o.overlaps[i - 1];
    }
    let left = if circular { 0 } else { c[0].left };
    let right = if circular { 0 } else { c[n - 1].right };
    if left.checked_add(right).is_none_or(|k| k >= length) {
        return Err(bad("outer resection leaves no paired product domain"));
    }
    let mut report=Report{engine_version:ENGINE_VERSION.into(),source_geometry_revisions:sources.iter().map(revision).collect::<Result<_,_>>()?,source_lengths:lengths.clone(),options:o.clone(),reference_length:length,reference_top_start:left,reference_top_end:length,product:None,mappings:vec![],resected:vec![],junctions:vec![],comparisons:3*lengths.iter().sum::<usize>()+2*length,terminal_search_complete:true,incomplete:false,warnings:vec!["Ideal declared 5′ resection, exact annealing, standard-dNTP gap filling and eligible nick sealing; not reaction kinetics, yield, host rescue, In-Fusion or SLIC chemistry.".into()],limits:vec![]};
    if !circular {
        report.warnings.push(match o.outer_end_policy {OuterEndPolicy::JoiningEndsOnly=>"JoiningEndsOnly explicitly exempts linear outer ends from exonuclease processing; this is an abstraction, not a claimed reagent selectivity.".into(),OuterEndPolicy::AllAccessibleEnds{..}=>"AllAccessibleEnds retains the declared outer 3′ overhangs; no template is invented to restore blunt ends.".into()});
    }
    let mut blocked_gap = false;
    for (i, h) in o.overlaps.iter().copied().enumerate() {
        let j = (i + 1) % n;
        for k in 0..h {
            report.comparisons += 1;
            if sources[i].top.as_str().as_bytes()[lengths[i] - h + k]
                != sources[j].top.as_str().as_bytes()[k]
            {
                return Err(bad("declared terminal homology is not exact"));
            }
        }
        let tg = c[j].left - h;
        let bg = c[i].right - h;
        let t3 = sources[i].right.three_prime;
        let b3 = sources[j].left.three_prime;
        let t5 = if c[j].left > 0 {
            FivePrime::Phosphate
        } else {
            sources[j].left.five_prime
        };
        let b5 = if c[i].right > 0 {
            FivePrime::Phosphate
        } else {
            sources[i].right.five_prime
        };
        blocked_gap |=
            (tg > 0 && t3 != ThreePrime::Hydroxyl) || (bg > 0 && b3 != ThreePrime::Hydroxyl);
        let boundary = |x: usize| -> i64 {
            if circular {
                (x % length) as i64
            } else {
                x as i64 - left as i64
            }
        };
        report.junctions.push(Junction {
            left_source: i,
            right_source: j,
            overlap: h,
            homology: sources[j].top.as_str()[..h].into(),
            top_gap: tg,
            bottom_gap: bg,
            top_nick: boundary(offsets[j] + c[j].left),
            bottom_nick: boundary(offsets[i] + lengths[i] - c[i].right),
            homology_start: boundary(offsets[j]),
            homology_end: boundary(offsets[j]) + h as i64,
            top_five_prime: t5,
            top_three_prime: t3,
            bottom_five_prime: b5,
            bottom_three_prime: b3,
            top_sealed: t5 == FivePrime::Phosphate && t3 == ThreePrime::Hydroxyl,
            bottom_sealed: b5 == FivePrime::Phosphate && b3 == ThreePrime::Hydroxyl,
            alternatives: vec![],
            terminal_search_complete: true,
        });
    }
    // Adjacent-pair terminal alternatives, bounded globally. Explicit declaration
    // remains the chosen geometry; repeats/caps are retained as uncertainty.
    'search: for i in 0..report.junctions.len() {
        let j = (i + 1) % n;
        for h in o.min_paired_bases..lengths[i].min(lengths[j]) {
            if h == o.overlaps[i] {
                continue;
            }
            let mut equal = true;
            for k in 0..h {
                if report.comparisons >= o.max_comparisons {
                    for pending in &mut report.junctions[i..] {
                        pending.terminal_search_complete = false;
                    }
                    report.terminal_search_complete = false;
                    report
                        .limits
                        .push("Global terminal-alternative comparison budget reached.".into());
                    break 'search;
                }
                report.comparisons += 1;
                if sources[i].top.as_str().as_bytes()[lengths[i] - h + k]
                    != sources[j].top.as_str().as_bytes()[k]
                {
                    equal = false;
                    break;
                }
            }
            if equal {
                if report
                    .junctions
                    .iter()
                    .map(|j| j.alternatives.len())
                    .sum::<usize>()
                    >= o.max_candidates
                {
                    for pending in &mut report.junctions[i..] {
                        pending.terminal_search_complete = false;
                    }
                    report.terminal_search_complete = false;
                    report
                        .limits
                        .push("Global alternative-retention budget reached.".into());
                    break 'search;
                }
                report.junctions[i].alternatives.push(h);
            }
        }
    }
    for i in 0..n {
        if c[i].left > 0 {
            report.resected.push(ResectedSpan {
                source: i,
                strand: Strand::Top,
                start: 0,
                end: c[i].left,
            });
        }
        if c[i].right > 0 {
            report.resected.push(ResectedSpan {
                source: i,
                strand: Strand::Bottom,
                start: 0,
                end: c[i].right,
            });
        }
    }
    report.incomplete = !report.terminal_search_complete
        || report
            .junctions
            .iter()
            .any(|j| !j.alternatives.is_empty() || !j.top_sealed || !j.bottom_sealed);
    if report.junctions.iter().any(|j| !j.alternatives.is_empty()) {
        report.warnings.push("Alternative adjacent terminal homologies exist; this product follows the explicit declaration, without a reaction-mixture uniqueness claim.".into());
    }
    if report
        .junctions
        .iter()
        .any(|j| !j.top_sealed || !j.bottom_sealed)
    {
        report.warnings.push("Unsealable or unknown junction chemistry remains as explicit nicks; no host rescue or phosphorylation of original ends is assumed.".into());
    }
    if blocked_gap {
        report.incomplete = true;
        report.limits.push(
            "A gap lacks an available, known 3′OH extension end; no product was fabricated.".into(),
        );
        bounded(&report)?;
        return Ok(report);
    }
    let mut top = vec![EMPTY; length];
    let mut bottom = vec![EMPTY; length];
    for i in 0..n {
        for (strand, cut, slots) in [
            (Strand::Top, c[i].left, &mut top),
            (Strand::Bottom, c[i].right, &mut bottom),
        ] {
            for k in cut..lengths[i] {
                let global = offsets[i]
                    + if strand == Strand::Top {
                        k
                    } else {
                        lengths[i] - 1 - k
                    };
                let at = if circular { global % length } else { global };
                if at >= length || slots[at].source != u8::MAX {
                    return Err(bad(
                        "retained strand pieces overlap or lie outside declared product",
                    ));
                }
                slots[at] = Cell {
                    source: i as u8,
                    index: k as u32,
                    synthesized: false,
                };
            }
        }
    }
    for at in 0..length {
        let need_top = at >= left;
        let need_bottom = at < length - right;
        if need_top && top[at].source == u8::MAX && (!need_bottom || bottom[at].source == u8::MAX)
            || need_bottom
                && bottom[at].source == u8::MAX
                && (!need_top || top[at].source == u8::MAX)
        {
            report.incomplete = true;
            report
                .limits
                .push("Resection left an untemplated gap; no product was fabricated.".into());
            bounded(&report)?;
            return Ok(report);
        }
    }
    for at in left..length {
        if top[at].source == u8::MAX {
            top[at] = Cell {
                synthesized: true,
                ..bottom[at]
            };
        }
    }
    for at in 0..length - right {
        if bottom[at].source == u8::MAX {
            bottom[at] = Cell {
                synthesized: true,
                ..top[at]
            };
        }
    }
    let top_cells = &top[left..];
    let bottom_cells: Vec<_> = bottom[..length - right].iter().rev().copied().collect();
    // Each strand is independently tiled from real surviving strands and actual
    // complementary templates. Never reverse-complement an assembled reference.
    let top_bytes: Vec<_> = top_cells
        .iter()
        .map(|c| cell_base(*c, Strand::Top, sources))
        .collect();
    let bottom_bytes: Vec<_> = bottom_cells
        .iter()
        .map(|c| cell_base(*c, Strand::Bottom, sources))
        .collect();
    let mut nicks = vec![];
    for j in &report.junctions {
        if !j.top_sealed {
            nicks.push(Nick {
                strand: Strand::Top,
                position: j.top_nick,
            });
        }
        if !j.bottom_sealed {
            nicks.push(Nick {
                strand: Strand::Bottom,
                position: j.bottom_nick,
            });
        }
    }
    let product = Duplex {
        top: Sequence::dna(std::str::from_utf8(&top_bytes).map_err(|_| bad("invalid top bases"))?)?,
        bottom: Sequence::dna(
            std::str::from_utf8(&bottom_bytes).map_err(|_| bad("invalid bottom bases"))?,
        )?,
        bottom_start: -(left as i64),
        topology: o.topology,
        left: Terminus {
            five_prime: if c[0].left > 0 {
                FivePrime::Phosphate
            } else {
                sources[0].left.five_prime
            },
            three_prime: sources[0].left.three_prime,
        },
        right: Terminus {
            five_prime: if c[n - 1].right > 0 {
                FivePrime::Phosphate
            } else {
                sources[n - 1].right.five_prime
            },
            three_prime: sources[n - 1].right.three_prime,
        },
        nicks,
    };
    report.mappings = mapping_runs(top_cells, Strand::Top);
    report
        .mappings
        .extend(mapping_runs(&bottom_cells, Strand::Bottom));
    report.product = Some(product);
    report.validate_documentary()?;
    Ok(report)
}
impl Report {
    /// Bounded documentary validation only. Does not replay source chemistry,
    /// establish source identity, or restore private tracked modification state.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        self.options.validate_shape(self.source_lengths.len())?;
        if self.engine_version != ENGINE_VERSION
            || self.source_lengths.is_empty()
            || self.source_lengths.len() > MAX_SOURCES
            || self.source_geometry_revisions.len() != self.source_lengths.len()
            || self
                .source_lengths
                .iter()
                .try_fold(0usize, |sum, n| sum.checked_add(*n))
                .is_none_or(|n| n > MAX_SOURCE_BASES / 2)
            || self.source_geometry_revisions.iter().any(|h| {
                h.len() != 64
                    || !h
                        .bytes()
                        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            })
            || self.junctions.len() != self.options.overlaps.len()
            || self.mappings.len() > MAX_SOURCES * 16
            || self.resected.len() > MAX_SOURCES * 2
            || self.warnings.len() > 100
            || self.limits.len() > 100
            || self
                .warnings
                .iter()
                .chain(&self.limits)
                .any(|s| s.len() > 4096)
            || self.comparisons > self.options.max_comparisons
        {
            return Err(bad("invalid documentary bounds/version"));
        }
        bounded(self)?;
        let cuts = cuts(&self.source_lengths, &self.options)?;
        let length =
            self.source_lengths.iter().sum::<usize>() - self.options.overlaps.iter().sum::<usize>();
        let circular = self.options.topology == Topology::Circular;
        let left = if circular { 0 } else { cuts[0].left };
        let right = if circular {
            0
        } else {
            cuts[cuts.len() - 1].right
        };
        if left.checked_add(right).is_none_or(|k| k >= length)
            || self.reference_length != length
            || self.reference_top_start != left
            || self.reference_top_end != length
        {
            return Err(bad("invalid documentary reference interval"));
        }
        let mut expected_resected = vec![];
        for (source, c) in cuts.iter().enumerate() {
            for (strand, end) in [(Strand::Top, c.left), (Strand::Bottom, c.right)] {
                if end > 0 {
                    expected_resected.push(ResectedSpan {
                        source,
                        strand,
                        start: 0,
                        end,
                    });
                }
            }
        }
        if self.resected != expected_resected {
            return Err(bad(
                "documentary resection evidence differs from declared cuts",
            ));
        }
        let mut offsets = vec![0usize; self.source_lengths.len()];
        for i in 1..offsets.len() {
            offsets[i] = offsets[i - 1] + self.source_lengths[i - 1] - self.options.overlaps[i - 1];
        }
        let boundary = |x: usize| -> i64 {
            if circular {
                (x % length) as i64
            } else {
                x as i64 - left as i64
            }
        };
        let mut alternatives = 0usize;
        for (i, j) in self.junctions.iter().enumerate() {
            alternatives = alternatives
                .checked_add(j.alternatives.len())
                .ok_or(Error::TooLarge)?;
            if j.left_source != i
                || j.right_source != (i + 1) % self.source_lengths.len()
                || j.overlap != self.options.overlaps[i]
                || j.homology.len() != j.overlap
                || !j
                    .homology
                    .bytes()
                    .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                || j.top_gap != cuts[j.right_source].left - j.overlap
                || j.bottom_gap != cuts[i].right - j.overlap
                || j.top_five_prime != FivePrime::Phosphate
                || j.bottom_five_prime != FivePrime::Phosphate
                || j.top_nick != boundary(offsets[j.right_source] + cuts[j.right_source].left)
                || j.bottom_nick != boundary(offsets[i] + self.source_lengths[i] - cuts[i].right)
                || j.homology_start != boundary(offsets[j.right_source])
                || j.homology_end != boundary(offsets[j.right_source]) + j.overlap as i64
                || j.alternatives.iter().any(|h| {
                    *h < self.options.min_paired_bases
                        || *h >= self.source_lengths[i].min(self.source_lengths[j.right_source])
                        || *h == j.overlap
                })
                || j.alternatives.windows(2).any(|w| w[0] >= w[1])
                || j.top_sealed
                    != (j.top_five_prime == FivePrime::Phosphate
                        && j.top_three_prime == ThreePrime::Hydroxyl)
                || j.bottom_sealed
                    != (j.bottom_five_prime == FivePrime::Phosphate
                        && j.bottom_three_prime == ThreePrime::Hydroxyl)
            {
                return Err(bad("invalid documentary junction evidence"));
            }
        }
        if alternatives > self.options.max_candidates
            || self.terminal_search_complete
                != self.junctions.iter().all(|j| j.terminal_search_complete)
            || ((!self.terminal_search_complete
                || alternatives > 0
                || self.product.is_none()
                || self
                    .junctions
                    .iter()
                    .any(|j| !j.top_sealed || !j.bottom_sealed))
                && !self.incomplete)
        {
            return Err(bad("documentary uncertainty was erased"));
        }
        for r in &self.resected {
            if r.source >= self.source_lengths.len()
                || r.start != 0
                || r.start >= r.end
                || r.end > self.source_lengths[r.source]
            {
                return Err(bad("invalid resected interval"));
            }
        }
        if let Some(p) = &self.product {
            if p.top.len() > MAX_SOURCE_BASES / 2
                || p.bottom.len() > MAX_SOURCE_BASES / 2
                || p.nicks.len() > MAX_SOURCES * 2
            {
                return Err(Error::TooLarge);
            }
            p.validate()?;
            if p.topology != self.options.topology
                || p.top.len() != length - left
                || p.bottom.len() != length - right
                || p.bottom_start != -(left as i64)
            {
                return Err(bad("documentary product topology/extents differ"));
            }
            let mut expected_nicks = vec![];
            for j in &self.junctions {
                if (j.top_gap > 0 && j.top_three_prime != ThreePrime::Hydroxyl)
                    || (j.bottom_gap > 0 && j.bottom_three_prime != ThreePrime::Hydroxyl)
                {
                    return Err(bad("documentary product invents unsupported gap fill"));
                }
                if !j.top_sealed {
                    expected_nicks.push(Nick {
                        strand: Strand::Top,
                        position: j.top_nick,
                    });
                }
                if !j.bottom_sealed {
                    expected_nicks.push(Nick {
                        strand: Strand::Bottom,
                        position: j.bottom_nick,
                    });
                }
            }
            if p.nicks != expected_nicks {
                return Err(bad("documentary product erased or added junction nicks"));
            }
            for m in &self.mappings {
                let x = &m.mapping;
                let size = if x.product_strand == Strand::Top {
                    p.top.len()
                } else {
                    p.bottom.len()
                };
                if m.source >= self.source_lengths.len()
                    || x.source_start >= x.source_end
                    || x.source_end > self.source_lengths[m.source]
                    || x.product_start >= x.product_end
                    || x.product_end > size
                    || x.source_end - x.source_start != x.product_end - x.product_start
                    || x.reverse != x.synthesized
                    || (x.source_strand != x.product_strand) != x.synthesized
                {
                    return Err(bad("invalid documentary strand mapping"));
                }
                let cut = if x.source_strand == Strand::Top {
                    cuts[m.source].left
                } else {
                    cuts[m.source].right
                };
                if x.source_start < cut {
                    return Err(bad("mapping reuses resected source material"));
                }
            }
            for strand in [Strand::Top, Strand::Bottom] {
                let mut spans: Vec<_> = self
                    .mappings
                    .iter()
                    .filter(|m| m.mapping.product_strand == strand)
                    .map(|m| (m.mapping.product_start, m.mapping.product_end))
                    .collect();
                spans.sort_unstable();
                let mut end = 0;
                for (a, b) in spans {
                    if a != end {
                        return Err(bad("product mapping coverage overlaps or is incomplete"));
                    }
                    end = b;
                }
                if end
                    != if strand == Strand::Top {
                        p.top.len()
                    } else {
                        p.bottom.len()
                    }
                {
                    return Err(bad("product mapping coverage is incomplete"));
                }
            }
        } else if !self.mappings.is_empty() {
            return Err(bad("withheld product has product mappings"));
        }
        Ok(())
    }
}
