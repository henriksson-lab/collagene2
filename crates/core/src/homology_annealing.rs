//! Ordered annealing of explicitly3′-resected intact source duplexes.
//! Every surviving covalent strand remains separate. The contour is a reference
//! coordinate system, never a synthesized/sealed construct or host-rescue model.
use crate::{
    Error, Sequence, Topology,
    chemistry::{Duplex, FivePrime, Nick, Strand, ThreePrime},
    three_prime_resection::{self as resection, Cuts},
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "homology-annealing-2";
pub const MAX_SOURCES: usize = 32;
pub const MAX_SOURCE_BASES: usize = 2_000_000;
pub const MAX_COMPARISONS: usize = 100_000_000;
pub const MAX_CANDIDATES: usize = 1000;
macro_rules! dto {($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum PlannedRescue {
    EcoliTransformation,
}
impl PlannedRescue {
    pub fn label(self) -> &'static str {
        "Planned E. coli transformation/rescue (not simulated)"
    }
}
// Original top-reference coverage only; bottom-only calls have no invented
// top annotation mapping. A circular contour may cross the origin.
dto! {pub struct SourcePlacement {
    pub source:usize,pub reference_start:usize,pub contour_length:usize,pub source_top_end:usize,
}}
dto! {pub struct Options {
    pub overlaps:Vec<usize>,pub topology:Topology,pub cuts:Vec<Cuts>,
    pub min_paired_bases:usize,pub max_comparisons:usize,pub max_candidates:usize,
    #[serde(deserialize_with="crate::assembly_recipe::required_option")]
    pub planned_rescue:Option<PlannedRescue>,
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
    pub fn validate_profile(
        &self,
        method: crate::overlap_profiles::OverlapMethod,
    ) -> Result<(), Error> {
        if method == crate::overlap_profiles::OverlapMethod::InFusion
            && self.planned_rescue != Some(PlannedRescue::EcoliTransformation)
        {
            return Err(bad(
                "In-Fusion requires an explicit planned E. coli transformation/rescue declaration; rescue is not simulated",
            ));
        }
        Ok(())
    }
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
// Coordinates run between the actual 5′ termini. Any original 3′ extension
// outside that contour must be explicitly removed, never silently clipped.
fn domain(p: &Duplex) -> Result<(usize, usize, usize), Error> {
    let start = usize::try_from(p.bottom_start)
        .map_err(|_| bad("remaining left 3′ extension requires explicit resection"))?;
    let end = start.checked_add(p.bottom.len()).ok_or(Error::TooLarge)?;
    if p.top.len() > end {
        return Err(bad(
            "remaining right 3′ extension requires explicit resection",
        ));
    }
    Ok((end, start, end - p.top.len()))
}
fn source_domain(s: &Duplex, c: &Cuts) -> Result<(usize, usize, usize), Error> {
    resection::preflight(s, c)?;
    let start = s
        .bottom_start
        .checked_add(c.left as i64)
        .ok_or(Error::TooLarge)?;
    let end = s
        .bottom_start
        .checked_add(s.bottom.len() as i64)
        .ok_or(Error::TooLarge)?;
    if start < 0 || end < (s.top.len() - c.right) as i64 {
        return Err(bad("remaining 3′ extension requires explicit resection"));
    }
    Ok((
        end as usize,
        start as usize,
        end as usize - (s.top.len() - c.right),
    ))
}
pub fn preflight_borrowed(sources: &[&Duplex], options: &Options) -> Result<(), Error> {
    options.validate_shape(sources.len())?;
    let mut bases = 0usize;
    for s in sources {
        bases = bases
            .checked_add(s.top.len())
            .and_then(|n| n.checked_add(s.bottom.len()))
            .filter(|n| *n <= MAX_SOURCE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    let domains = sources
        .iter()
        .zip(&options.cuts)
        .map(|(s, c)| source_domain(s, c))
        .collect::<Result<Vec<_>, _>>()?;
    for (i, &h) in options.overlaps.iter().enumerate() {
        let j = (i + 1) % sources.len();
        if h >= domains[i].0.min(domains[j].0) || h > domains[i].2 || h > domains[j].1 {
            return Err(bad(
                "each chosen homology must be proper and completely exposed in both5′tails",
            ));
        }
    }
    Ok(())
}
fn contour(p: &Duplex) -> Result<String, Error> {
    let (end, _, _) = domain(p)?;
    let mut text = p.top.as_str().to_owned();
    for x in p.top.len()..end {
        text.push(crate::model::complement(
            p.bottom.as_str().as_bytes()[end - 1 - x],
            crate::Alphabet::Dna,
        ) as char);
    }
    Ok(text)
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
    let prepared = sources
        .iter()
        .zip(&options.cuts)
        .map(|(s, c)| resection::geometry(s, c))
        .collect::<Result<Vec<_>, _>>()?;
    let contours = prepared
        .iter()
        .map(|p| contour(&p.product))
        .collect::<Result<Vec<_>, _>>()?;
    let domains = prepared
        .iter()
        .map(|p| domain(&p.product))
        .collect::<Result<Vec<_>, _>>()?;
    let circular = options.topology == Topology::Circular;
    let mut budget = Budget {
        used: 0,
        max: options.max_comparisons,
    };
    // Validate every explicitly chosen junction before exploring alternatives.
    // Target letters never select a pairing or supply missing biological data.
    for (i, &h) in options.overlaps.iter().enumerate() {
        let j = (i + 1) % n;
        let a = contours[i].as_bytes();
        let b = contours[j].as_bytes();
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
        offsets[i] = offsets[i - 1] + domains[i - 1].0 - options.overlaps[i - 1];
    }
    let last_end = offsets[n - 1] + domains[n - 1].0;
    let length = last_end - if circular { options.overlaps[n - 1] } else { 0 };
    if length == 0 || length > MAX_SOURCE_BASES / 2 {
        return Err(Error::TooLarge);
    }
    let mut text = String::with_capacity(length);
    for (i, source) in contours.iter().enumerate() {
        let trim = options.overlaps.get(i).copied().unwrap_or(0);
        text.push_str(&source[..source.len() - trim]);
    }
    if text.len() != length {
        return Err(bad("source contour length accounting differs"));
    }
    if domains.iter().any(|d| d.0 > length.saturating_mul(2)) {
        return Err(bad("source contour would require multiple circular wraps"));
    }
    let reference = Sequence::dna(&text)?;
    let mut preparations = Vec::with_capacity(n);
    let mut segments = Vec::with_capacity(2 * n);
    for (i, report) in prepared.into_iter().enumerate() {
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
            reference_start: offsets[i] + domains[i].1,
            reference_end: offsets[i] + domains[i].0,
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
        let a = contours[i].as_bytes();
        let b = contours[j].as_bytes();
        let mut alternatives = vec![];
        let mut complete = true;
        let maximum = domains[i].2.min(domains[j].1).min(a.len().min(b.len()) - 1);
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
            start: offsets[i] + preparations[i].product.top.len(),
            end: right_offset,
            five_prime: preparations[j].product.left.five_prime,
            three_prime: preparations[i].product.right.three_prime,
        };
        let bottom = Discontinuity {
            strand: Strand::Bottom,
            start: left_end,
            end: right_offset + domains[j].1,
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
            homology: contours[j][..h].into(),
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
            length - (offsets[n - 1] + preparations[n - 1].product.top.len())
        };
        let left_cut = if circular { 0 } else { domains[0].1 };
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

impl Report {
    pub fn source_placements(&self) -> Result<Vec<SourcePlacement>, Error> {
        self.validate_documentary()?;
        let mut at = 0;
        self.preparations
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let (length, _, _) = domain(&p.product)?;
                let placement = SourcePlacement {
                    source: i,
                    reference_start: at,
                    contour_length: length,
                    source_top_end: (p.product.top.len() + p.cuts.right).min(length),
                };
                at += length - self.options.overlaps.get(i).copied().unwrap_or(0);
                Ok(placement)
            })
            .collect()
    }

    /// Linear bounded literal/shape checks, without source recovery, terminal
    /// search, polymerase work or any conversion into private tracked material.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        let n = self.preparations.len();
        self.options.validate_shape(n)?;
        if self.engine_version != ENGINE_VERSION
            || self.source_geometry_revisions.len() != n
            || self.segments.len() != 2 * n
            || self.junctions.len() != self.options.overlaps.len()
            || self.reference.alphabet() != crate::Alphabet::Dna
            || self.reference.len() > MAX_SOURCE_BASES / 2
            || self.comparisons > self.options.max_comparisons
            || self.warnings.len() > 1000
            || self.limits.len() > 1000
            || self
                .warnings
                .iter()
                .chain(&self.limits)
                .any(|s| s.len() > 4096)
        {
            return Err(bad("report engine/shape/bounds mismatch"));
        }
        let mut lengths = Vec::with_capacity(n);
        let mut domains = Vec::with_capacity(n);
        let mut total = 0usize;
        for (i, p) in self.preparations.iter().enumerate() {
            let length = p
                .product
                .top
                .len()
                .checked_add(p.cuts.right)
                .ok_or(Error::TooLarge)?;
            if p.cuts != self.options.cuts[i]
                || p.source_geometry_revision != self.source_geometry_revisions[i]
            {
                return Err(bad(
                    "preparation geometry differs from declared source cuts",
                ));
            }
            total = total
                .checked_add(
                    length
                        .checked_add(p.product.bottom.len())
                        .and_then(|n| n.checked_add(p.cuts.left))
                        .ok_or(Error::TooLarge)?,
                )
                .filter(|n| *n <= MAX_SOURCE_BASES)
                .ok_or(Error::TooLarge)?;
            p.validate_documentary()?;
            let d = domain(&p.product)?;
            lengths.push(d.0);
            domains.push(d);
        }
        let circular = self.options.topology == Topology::Circular;
        let mut offsets = vec![0usize; n];
        for (i, &h) in self.options.overlaps.iter().enumerate() {
            let j = (i + 1) % n;
            if h >= lengths[i].min(lengths[j]) || domains[i].2 < h || domains[j].1 < h {
                return Err(bad("report homology exposure differs"));
            }
        }
        for i in 1..n {
            offsets[i] = offsets[i - 1] + lengths[i - 1] - self.options.overlaps[i - 1];
        }
        let length = offsets[n - 1] + lengths[n - 1]
            - if circular {
                self.options.overlaps[n - 1]
            } else {
                0
            };
        if length == 0
            || self.reference.len() != length
            || lengths.iter().any(|n| *n > length.saturating_mul(2))
        {
            return Err(bad("reference contour length differs from copied sources"));
        }
        let reference = self.reference.as_str().as_bytes();
        for i in 0..n {
            let p = &self.preparations[i].product;
            for (part, strand, sequence, start, end, five, three) in [
                (
                    0,
                    Strand::Top,
                    &p.top,
                    offsets[i],
                    offsets[i] + p.top.len(),
                    p.left.five_prime,
                    p.right.three_prime,
                ),
                (
                    1,
                    Strand::Bottom,
                    &p.bottom,
                    offsets[i] + domains[i].1,
                    offsets[i] + lengths[i],
                    p.right.five_prime,
                    p.left.three_prime,
                ),
            ] {
                let s = &self.segments[2 * i + part];
                if s.source != i
                    || s.strand != strand
                    || s.source_start != 0
                    || s.source_end != sequence.len()
                    || s.sequence != *sequence
                    || s.reference_start != start
                    || s.reference_end != end
                    || s.five_prime != five
                    || s.three_prime != three
                {
                    return Err(bad("surviving strand segment changed"));
                }
                for (k, &b) in sequence.as_str().as_bytes().iter().enumerate() {
                    let at = if strand == Strand::Top {
                        start + k
                    } else {
                        end - 1 - k
                    };
                    if (!circular && at >= length)
                        || if strand == Strand::Top {
                            reference[at % length] != b
                        } else {
                            crate::model::complement(reference[at % length], crate::Alphabet::Dna)
                                != b
                        }
                    {
                        return Err(bad(
                            "coordinate reference differs from an actual surviving source strand",
                        ));
                    }
                }
            }
        }
        let mut retained = 0usize;
        let mut all_complete = true;
        let mut gapless = true;
        for (i, junction) in self.junctions.iter().enumerate() {
            let j = (i + 1) % n;
            let h = self.options.overlaps[i];
            let right = if j == 0 { length } else { offsets[j] };
            let left = offsets[i] + lengths[i];
            let top = Discontinuity {
                strand: Strand::Top,
                start: offsets[i] + self.preparations[i].product.top.len(),
                end: right,
                five_prime: self.preparations[j].product.left.five_prime,
                three_prime: self.preparations[i].product.right.three_prime,
            };
            let bottom = Discontinuity {
                strand: Strand::Bottom,
                start: left,
                end: right + domains[j].1,
                five_prime: self.preparations[i].product.right.five_prime,
                three_prime: self.preparations[j].product.left.three_prime,
            };
            if junction.left_source != i
                || junction.right_source != j
                || junction.overlap != h
                || junction.homology_start != right
                || junction.homology_end != left
                || junction.top != top
                || junction.bottom != bottom
                || junction.homology != self.preparations[j].product.top.as_str()[..h]
            {
                return Err(bad("junction/nick/gap declaration changed"));
            }
            let a = &self.preparations[i].product.bottom.as_str().as_bytes()[..h];
            if a.iter()
                .zip(junction.homology.as_bytes().iter().rev())
                .any(|(&a, &b)| a != crate::model::complement(b, crate::Alphabet::Dna))
            {
                return Err(bad("declared homology differs from actual exposed tails"));
            }
            let maximum = domains[i]
                .2
                .min(domains[j].1)
                .min(lengths[i].min(lengths[j]) - 1);
            if junction
                .alternatives
                .iter()
                .any(|v| *v < self.options.min_paired_bases || *v > maximum || *v == h)
                || junction.alternatives.windows(2).any(|p| p[0] >= p[1])
            {
                return Err(bad("invalid terminal alternative lengths"));
            }
            retained = retained
                .checked_add(junction.alternatives.len())
                .filter(|n| *n <= self.options.max_candidates)
                .ok_or(Error::TooLarge)?;
            all_complete &= junction.terminal_search_complete;
            gapless &= top.start == top.end && bottom.start == bottom.end;
        }
        if self.terminal_search_complete != all_complete
            || self.incomplete != (!self.junctions.is_empty() || !all_complete)
            || self.paired_diagnostic.is_some() != gapless
        {
            return Err(bad("report completeness/diagnostic scope changed"));
        }
        if let Some(d) = &self.paired_diagnostic {
            let left = if circular { 0 } else { domains[0].1 };
            let right = if circular {
                0
            } else {
                length - (offsets[n - 1] + self.preparations[n - 1].product.top.len())
            };
            if d.topology != self.options.topology
                || d.bottom_start != left as i64
                || d.top.as_str() != &self.reference.as_str()[..length - right]
                || d.bottom.len() != length - left
                || d.left != self.preparations[0].product.left
                || d.right != self.preparations[n - 1].product.right
                || d.nicks.len() != 2 * self.junctions.len()
            {
                return Err(bad(
                    "diagnostic duplex domain differs from retained strands",
                ));
            }
            if d.bottom
                .as_str()
                .bytes()
                .zip(reference[left..].iter().rev())
                .any(|(a, &b)| a != crate::model::complement(b, crate::Alphabet::Dna))
            {
                return Err(bad("diagnostic bottom sequence changed"));
            }
            for (i, j) in self.junctions.iter().enumerate() {
                if d.nicks[2 * i]
                    != (Nick {
                        strand: Strand::Top,
                        position: projection_point(j.top.start, length, circular) as i64,
                    })
                    || d.nicks[2 * i + 1]
                        != (Nick {
                            strand: Strand::Bottom,
                            position: projection_point(j.bottom.start, length, circular) as i64,
                        })
                {
                    return Err(bad("diagnostic joining nicks changed or were sealed"));
                }
            }
            d.validate()?;
        }
        resection::revision(self)?;
        Ok(())
    }
}
