//! Geometry-only restriction-site classification relative to selected reference bases.
//!
//! This neither scans sequence nor establishes physical cleavage. In particular,
//! `Inside` never upgrades possible recognition evidence or an incomplete scan.

use crate::restriction_inspection::{InspectionSite, SiteDiagnostic};
use collagene_model::{Error, Location, Topology};

/// Maximum input chunks accepted when constructing a union. Chunks may overlap.
pub const MAX_REGION_SPANS: usize = 10_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    length: usize,
    topology: Topology,
    spans: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Basis {
    RecognitionBases,
    StrandCuts,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Relation {
    Inside,
    Outside,
    Crossing,
    Undetermined,
}

fn bad(message: &str) -> Error {
    Error::Unsupported(format!("restriction selection: {message}"))
}

/// Extract exact selected base spans in declaration order, ignoring strand.
/// Only nonempty `Range`, `Join` and `Complement` geometry is accepted;
/// fuzzy endpoints, `Order` and unresolved locations cannot define this filter.
/// Repeated/overlapping spans are allowed and are merged by `Region` constructors.
pub fn footprint(location: &Location, length: usize) -> Result<Vec<(usize, usize)>, Error> {
    if length == 0 || length > collagene_model::MAX_SEQUENCE_LEN {
        return Err(bad("invalid reference length"));
    }
    fn visit(
        location: &Location,
        length: usize,
        depth: usize,
        nodes: &mut usize,
        spans: &mut Vec<(usize, usize)>,
    ) -> Result<(), Error> {
        if depth > 32 || *nodes >= 100_000 {
            return Err(bad("selection location exceeds traversal bounds"));
        }
        *nodes += 1;
        match location {
            Location::Range {
                start,
                end,
                before,
                after,
            } => {
                if *before || *after || start >= end || *end > length {
                    return Err(bad("selection requires exact nonempty reference ranges"));
                }
                if spans.len() == MAX_REGION_SPANS {
                    return Err(bad("too many selection chunks"));
                }
                spans.push((*start, *end));
            }
            Location::Complement(inner) => visit(inner, length, depth + 1, nodes, spans)?,
            Location::Join(parts) => {
                if parts.is_empty() || parts.len() > 100_000 - *nodes {
                    return Err(bad("selection join is empty or exceeds traversal bounds"));
                }
                for part in parts {
                    visit(part, length, depth + 1, nodes, spans)?;
                }
            }
            Location::Order(_) | Location::Unresolved(_) => {
                return Err(bad(
                    "selection requires exact range/join/complement geometry",
                ));
            }
        }
        Ok(())
    }
    let mut spans = Vec::new();
    visit(location, length, 0, &mut 0, &mut spans)?;
    Ok(spans)
}

impl Region {
    /// Selected zero-based, half-open base intervals on a linear reference.
    /// Empty selections (including a caret) are errors. Overlapping and adjacent
    /// chunks are merged. A wrapping selection is supplied as two ordinary spans.
    pub fn new(
        length: usize,
        spans: impl IntoIterator<Item = (usize, usize)>,
    ) -> Result<Self, Error> {
        Self::with_topology(length, Topology::Linear, spans)
    }

    /// As `new`, but reference boundary 0 and boundary `length` are equivalent.
    pub fn circular(
        length: usize,
        spans: impl IntoIterator<Item = (usize, usize)>,
    ) -> Result<Self, Error> {
        Self::with_topology(length, Topology::Circular, spans)
    }

    pub fn with_topology(
        length: usize,
        topology: Topology,
        spans: impl IntoIterator<Item = (usize, usize)>,
    ) -> Result<Self, Error> {
        if length == 0 || length > collagene_model::MAX_SEQUENCE_LEN {
            return Err(bad("invalid reference length"));
        }
        let mut chunks = Vec::new();
        for (start, end) in spans {
            if chunks.len() == MAX_REGION_SPANS {
                return Err(bad("too many selection chunks"));
            }
            if start >= end || end > length {
                return Err(bad("select nonempty base ranges within the reference"));
            }
            chunks.push((start, end));
        }
        if chunks.is_empty() {
            return Err(bad("select a base range"));
        }
        chunks.sort_unstable();
        let mut union: Vec<(usize, usize)> = Vec::with_capacity(chunks.len());
        for (start, end) in chunks {
            if let Some(last) = union.last_mut().filter(|last| start <= last.1) {
                last.1 = last.1.max(end);
            } else {
                union.push((start, end));
            }
        }
        Ok(Self {
            length,
            topology,
            spans: union,
        })
    }

    pub fn length(&self) -> usize {
        self.length
    }
    pub fn topology(&self) -> Topology {
        self.topology
    }
    /// Sorted, disjoint, nonadjacent half-open base ranges.
    pub fn spans(&self) -> &[(usize, usize)] {
        &self.spans
    }

    fn contains_boundary(&self, boundary: usize) -> bool {
        // Cuts are reference boundaries, not bases: BOTH selection endpoints
        // are included. A circular origin belongs to either adjoining interval.
        self.spans.iter().any(|&(start, end)| {
            start <= boundary && boundary <= end
                || self.topology == Topology::Circular
                    && (boundary == 0 || boundary == self.length)
                    && (start == 0 || end == self.length)
        })
    }
}

/// Classify one borrowed inspection site without changing its evidence status.
///
/// Recognition uses base overlap. Strand cuts include both endpoints of every
/// selected range; circular 0 and `length` denote the same boundary. Both strand
/// cuts must be known and unambiguous. Missing cuts or outside-substrate /
/// ambiguous-strand diagnostics yield `Undetermined` for the cut basis.
///
/// Coordinates and the exact recognition footprint are validated before work.
/// The caller must retain the scan's completeness and each site's `certain`
/// flag when presenting counts. This function has no scan-completeness input.
pub fn classify(site: &InspectionSite, region: &Region, basis: &Basis) -> Result<Relation, Error> {
    validate_site(site, region)?;
    match basis {
        Basis::RecognitionBases => {
            let mut overlap = 0;
            for footprint in &site.spans {
                for &(start, end) in &region.spans {
                    overlap += footprint
                        .end
                        .min(end)
                        .saturating_sub(footprint.start.max(start));
                }
            }
            Ok(if overlap == 0 {
                Relation::Outside
            } else if overlap == site.recognition_length {
                Relation::Inside
            } else {
                Relation::Crossing
            })
        }
        Basis::StrandCuts => {
            if site.diagnostics.iter().any(|d| {
                matches!(
                    d,
                    SiteDiagnostic::CutsOutsideSubstrate | SiteDiagnostic::AmbiguousStrandGeometry
                )
            }) {
                return Ok(Relation::Undetermined);
            }
            let (Some(top), Some(bottom)) = (site.top_cut, site.bottom_cut) else {
                return Ok(Relation::Undetermined);
            };
            Ok(
                match (
                    region.contains_boundary(top),
                    region.contains_boundary(bottom),
                ) {
                    (true, true) => Relation::Inside,
                    (false, false) => Relation::Outside,
                    _ => Relation::Crossing,
                },
            )
        }
    }
}

fn validate_site(site: &InspectionSite, region: &Region) -> Result<(), Error> {
    let n = region.length;
    if site.spans.is_empty()
        || site.spans.len() > 2
        || site.diagnostics.len() > 3
        || site.recognition_start >= n
        || site.recognition_length == 0
        || site.recognition_length > 64
        || site.recognition_length > n
        || site.top_cut.is_some_and(|cut| cut > n)
        || site.bottom_cut.is_some_and(|cut| cut > n)
    {
        return Err(bad("invalid site geometry or bounds"));
    }
    let first = site.recognition_length.min(n - site.recognition_start);
    let wraps = first < site.recognition_length;
    if wraps && region.topology != Topology::Circular
        || site.spans.len() != if wraps { 2 } else { 1 }
        || site.spans[0].start != site.recognition_start
        || site.spans[0].end != site.recognition_start + first
        || wraps
            && (site.spans[1].start != 0 || site.spans[1].end != site.recognition_length - first)
    {
        return Err(bad("recognition spans disagree with reference coordinates"));
    }
    for (cut, raw) in [
        (site.top_cut, site.top_boundary),
        (site.bottom_cut, site.bottom_boundary),
    ] {
        if let Some(cut) = cut {
            let agrees = if region.topology == Topology::Circular {
                cut % n == raw.rem_euclid(n as i64) as usize
            } else {
                raw >= 0 && raw as u64 == cut as u64
            };
            if !agrees {
                return Err(bad("strand boundary disagrees with cut coordinate"));
            }
        }
    }
    Ok(())
}
