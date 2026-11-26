//! Ordered annealing of explicitly3′-resected intact source duplexes.
//! Every surviving covalent strand remains separate. The contour is a reference
//! coordinate system, never a synthesized/sealed construct or host-rescue model.
use crate::{
    Error, Sequence, Topology,
    chemistry::{Duplex, FivePrime, Nick, Strand, ThreePrime},
    three_prime_resection::{self as resection, Cuts},
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "homology-annealing-1";
pub const MAX_SOURCES: usize = 32;
pub const MAX_SOURCE_BASES: usize = 2_000_000;
pub const MAX_COMPARISONS: usize = 100_000_000;
pub const MAX_CANDIDATES: usize = 1000;
macro_rules! dto {($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct Options {
    pub overlaps:Vec<usize>,pub topology:Topology,pub cuts:Vec<Cuts>,
    pub min_paired_bases:usize,pub max_comparisons:usize,pub max_candidates:usize,
}}
// A whole surviving covalent strand. Reference coordinates increase on top,
// while actual stored bottom5′→3′ runs from reference_end to reference_start.
// The final bottom segment may cross the circular origin (end>reference.len).
dto! {pub struct Segment {
    pub source:usize,pub strand:Strand,pub source_start:usize,pub source_end:usize,
    pub sequence:Sequence,pub reference_start:usize,pub reference_end:usize,
    pub five_prime:FivePrime,pub three_prime:ThreePrime,
}}
// A gap or nick joins adjacent *separate* covalent segments on one strand.
// start==end denotes a nick. Circle coordinates are unwrapped; use moduloL
// for display only, preserving the real interval length.
dto! {pub struct Discontinuity {
    pub strand:Strand,pub start:usize,pub end:usize,
    pub five_prime:FivePrime,pub three_prime:ThreePrime,
}}
dto! {pub struct Junction {
    pub left_source:usize,pub right_source:usize,pub overlap:usize,pub homology:String,
    pub homology_start:usize,pub homology_end:usize,
    pub top:Discontinuity,pub bottom:Discontinuity,
    pub alternatives:Vec<usize>,pub terminal_search_complete:bool,
}}
dto! {pub struct Report {
    pub engine_version:String,pub source_geometry_revisions:Vec<String>,pub options:Options,
    /// Original-source-derived coordinate reference only. Gaps are not filled.
    pub reference:Sequence,pub preparations:Vec<resection::Report>,
    pub segments:Vec<Segment>,pub junctions:Vec<Junction>,
    /// Only gapless strand coverage fits Duplex; retained nicks are diagnostic.
    /// This field never conveys a private reusable or covalently closed product.
    pub paired_diagnostic:Option<Duplex>,
    pub comparisons:usize,pub terminal_search_complete:bool,pub incomplete:bool,
    pub warnings:Vec<String>,pub limits:Vec<String>,
}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("homology annealing: {s}"))
}
impl Options {
    pub fn validate_shape(&self, n: usize) -> Result<(), Error> {
        if !(1..=MAX_SOURCES).contains(&n)
            || self.cuts.len() != n
            || self.overlaps.len() != n - usize::from(self.topology != Topology::Circular)
        {
            return Err(bad("ordered source, cut and junction counts differ"));
        }
        if self.min_paired_bases == 0
            || self.min_paired_bases > MAX_SOURCE_BASES / 2
            || self
                .overlaps
                .iter()
                .any(|h| *h < self.min_paired_bases || *h > MAX_SOURCE_BASES / 2)
            || self.max_comparisons == 0
            || self.max_comparisons > MAX_COMPARISONS
            || self.max_candidates == 0
            || self.max_candidates > MAX_CANDIDATES
        {
            return Err(bad("invalid pairing/work/candidate bounds"));
        }
        if self
            .cuts
            .iter()
            .any(|c| c.left > MAX_SOURCE_BASES / 2 || c.right > MAX_SOURCE_BASES / 2)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
// Ordered assembly currently accepts intact blunt operands. General free-end
// clipping remains available independently via three_prime_resection.
pub fn preflight_borrowed(sources: &[&Duplex], options: &Options) -> Result<(), Error> {
    options.validate_shape(sources.len())?;
    let mut bases = 0usize;
    for source in sources {
        bases = bases
            .checked_add(source.top.len())
            .and_then(|n| n.checked_add(source.bottom.len()))
            .filter(|n| *n <= MAX_SOURCE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    for (source, cuts) in sources.iter().zip(&options.cuts) {
        if source.bottom_start != 0 || source.top.len() != source.bottom.len() {
            return Err(bad(
                "ordered annealing requires initially blunt intact duplexes; unsupported end geometry must not be flattened",
            ));
        }
        resection::preflight(source, cuts)?;
    }
    for (i, &h) in options.overlaps.iter().enumerate() {
        let j = (i + 1) % sources.len();
        if h >= sources[i].top.len().min(sources[j].top.len())
            || options.cuts[i].right < h
            || options.cuts[j].left < h
        {
            return Err(bad(
                "each chosen homology must be proper and completely exposed in both5′tails",
            ));
        }
    }
    Ok(())
}
struct Budget {
    used: usize,
    max: usize,
}
impl Budget {
    fn equal(&mut self, a: &[u8], b: &[u8]) -> Option<bool> {
        for (&x, &y) in a.iter().zip(b) {
            if self.used == self.max {
                return None;
            }
            self.used += 1;
            if x != y {
                return Some(false);
            }
        }
        Some(true)
    }
}
fn projection_point(point: usize, len: usize, circular: bool) -> usize {
    if circular { point % len } else { point }
}
pub fn geometry(sources: &[Duplex], options: &Options) -> Result<Report, Error> {
    options.validate_shape(sources.len())?;
    geometry_borrowed(&sources.iter().collect::<Vec<_>>(), options)
}
// Borrowed sources avoid cloning full DNA before native bounds are checked.
pub fn geometry_borrowed(sources: &[&Duplex], options: &Options) -> Result<Report, Error> {
    preflight_borrowed(sources, options)?;
    let n = sources.len();
    let circular = options.topology == Topology::Circular;
    let mut budget = Budget {
        used: 0,
        max: options.max_comparisons,
    };
    // Validate every explicitly chosen junction before exploring alternatives.
    // Target letters never select a pairing or supply missing biological data.
    for (i, &h) in options.overlaps.iter().enumerate() {
        let j = (i + 1) % n;
        let a = sources[i].top.as_str().as_bytes();
        let b = sources[j].top.as_str().as_bytes();
        match budget.equal(&a[a.len() - h..], &b[..h]) {
            Some(true) => {}
            Some(false) => {
                return Err(bad(
                    "declared terminal homology differs between actual source strands",
                ));
            }
            None => {
                return Err(bad(
                    "work limit prevents validation of the selected homology",
                ));
            }
        }
    }
    let mut offsets = vec![0usize; n];
    for i in 1..n {
        offsets[i] = offsets[i - 1] + sources[i - 1].top.len() - options.overlaps[i - 1];
    }
    let last_end = offsets[n - 1] + sources[n - 1].top.len();
    let length = last_end - if circular { options.overlaps[n - 1] } else { 0 };
    if length == 0 || length > MAX_SOURCE_BASES / 2 {
        return Err(Error::TooLarge);
    }
    let mut text = String::with_capacity(length);
    for (i, source) in sources.iter().enumerate() {
        let trim = options.overlaps.get(i).copied().unwrap_or(0);
        text.push_str(&source.top.as_str()[..source.top.len() - trim]);
    }
    if text.len() != length {
        return Err(bad("source contour length accounting differs"));
    }
    let reference = Sequence::dna(&text)?;
    let mut preparations = Vec::with_capacity(n);
    let mut segments = Vec::with_capacity(2 * n);
    for (i, source) in sources.iter().enumerate() {
        let report = resection::geometry(source, &options.cuts[i])?;
        let p = &report.product;
        segments.push(Segment {
            source: i,
            strand: Strand::Top,
            source_start: 0,
            source_end: p.top.len(),
            sequence: p.top.clone(),
            reference_start: offsets[i],
            reference_end: offsets[i] + p.top.len(),
            five_prime: p.left.five_prime,
            three_prime: p.right.three_prime,
        });
        segments.push(Segment {
            source: i,
            strand: Strand::Bottom,
            source_start: 0,
            source_end: p.bottom.len(),
            sequence: p.bottom.clone(),
            reference_start: offsets[i] + options.cuts[i].left,
            reference_end: offsets[i] + source.top.len(),
            five_prime: p.right.five_prime,
            three_prime: p.left.three_prime,
        });
        preparations.push(report);
    }
    let mut junctions = Vec::with_capacity(options.overlaps.len());
    let mut candidate_count = 0usize;
    let mut search_complete = true;
    let mut limits = vec![];
    for (i, &h) in options.overlaps.iter().enumerate() {
        let j = (i + 1) % n;
        let a = sources[i].top.as_str().as_bytes();
        let b = sources[j].top.as_str().as_bytes();
        let mut alternatives = vec![];
        let mut complete = true;
        let maximum = options.cuts[i]
            .right
            .min(options.cuts[j].left)
            .min(a.len().min(b.len()) - 1);
        for size in options.min_paired_bases..=maximum {
            if size == h {
                continue;
            }
            match budget.equal(&a[a.len() - size..], &b[..size]) {
                None => {
                    complete = false;
                    break;
                }
                Some(false) => {}
                Some(true) => {
                    if candidate_count == options.max_candidates {
                        complete = false;
                        break;
                    }
                    candidate_count += 1;
                    alternatives.push(size);
                }
            }
        }
        search_complete &= complete;
        // Unwrap the closing right operand by one complete contour. Keeping
        // these actual boundaries avoids losing origin-spanning gap lengths.
        let right_offset = if j == 0 { length } else { offsets[j] };
        let left_end = offsets[i] + a.len();
        let top = Discontinuity {
            strand: Strand::Top,
            start: left_end - options.cuts[i].right,
            end: right_offset,
            five_prime: preparations[j].product.left.five_prime,
            three_prime: preparations[i].product.right.three_prime,
        };
        let bottom = Discontinuity {
            strand: Strand::Bottom,
            start: left_end,
            end: right_offset + options.cuts[j].left,
            five_prime: preparations[i].product.right.five_prime,
            three_prime: preparations[j].product.left.three_prime,
        };
        if top.start > top.end || bottom.start > bottom.end {
            return Err(bad(
                "resection leaves overlapping covalent strands instead of declared nicks/gaps",
            ));
        }
        junctions.push(Junction {
            left_source: i,
            right_source: j,
            overlap: h,
            homology: sources[j].top.as_str()[..h].into(),
            homology_start: right_offset,
            homology_end: left_end,
            top,
            bottom,
            alternatives,
            terminal_search_complete: complete,
        });
    }
    if !search_complete {
        limits.push("Bounded adjacent exposed-terminal search is incomplete; no mixture or pairing uniqueness is established.".into());
    }
    let gapless = junctions
        .iter()
        .all(|j| j.top.start == j.top.end && j.bottom.start == j.bottom.end);
    let paired_diagnostic = if gapless {
        let right_cut = if circular {
            0
        } else {
            options.cuts[n - 1].right
        };
        let left_cut = if circular { 0 } else { options.cuts[0].left };
        let mut nicks = vec![];
        for j in &junctions {
            nicks.push(Nick {
                strand: Strand::Top,
                position: projection_point(j.top.start, length, circular) as i64,
            });
            nicks.push(Nick {
                strand: Strand::Bottom,
                position: projection_point(j.bottom.start, length, circular) as i64,
            });
        }
        let product = Duplex {
            top: Sequence::dna(&reference.as_str()[..length - right_cut])?,
            bottom: Sequence::dna(&reference.as_str()[left_cut..])?.reverse_complement()?,
            bottom_start: left_cut as i64,
            topology: options.topology,
            left: preparations[0].product.left,
            right: preparations[n - 1].product.right,
            nicks,
        };
        product.validate()?;
        Some(product)
    } else {
        limits.push("Actual strand gaps remain unfilled. No Duplex product was manufactured from the coordinate reference.".into());
        None
    };
    if !junctions.is_empty() {
        limits.push("Joining strands remain separate at every nick/gap; no ligation, polymerase filling or cellular rescue is modeled. No reusable assembled material is released.".into());
    }
    let report=Report{engine_version:ENGINE_VERSION.into(),source_geometry_revisions:preparations.iter().map(|p|p.source_geometry_revision.clone()).collect(),options:options.clone(),reference,preparations,segments,junctions,paired_diagnostic,comparisons:budget.used,terminal_search_complete:search_complete,incomplete:!options.overlaps.is_empty()||!search_complete,warnings:vec!["Explicit generic3′ hydrolysis and exact5′overhang pairing only; not an In-Fusion formulation, resection-length prediction, experimental efficiency or E.coli rescue. Original5′ end groups remain unchanged.".into()],limits};
    // Streamed bounds include literal strand pieces; no report-sized temporary.
    resection::revision(&report)?;
    Ok(report)
}
// Exact fresh source-bound replay. Serialized diagnostics cannot restore a
// native tracked intermediate or authorize missing repair/host rescue.
pub fn replay(sources: &[Duplex], options: &Options, report: &Report) -> Result<(), Error> {
    if geometry(sources, options)? != *report {
        return Err(bad("report differs from fresh source-bound annealing"));
    }
    Ok(())
}
