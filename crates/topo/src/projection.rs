//! Fresh, source-bound proof of an ordinary antiparallel paired contour.
//! This module never repairs a nick, deletes an unpaired call, or takes a target
//! sequence as an origin/orientation guide. The full physical graph is retained.
use crate::{
    Call, ChainEnds, FivePrime, Limits, PairSpan, ProductGraph, Report, Step, Substrates,
    ThreePrime,
};
use collagene_model::Sequence;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

pub const ENGINE: &str = "topo-paired-contour-proof-1";
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub source_batch_revision: String,
    pub component: usize,
    pub reference_seed: Call,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProjectionLimits {
    pub max_pairs: usize,
    pub max_work_cells: usize,
    pub max_work_bytes: usize,
    pub max_mapping_segments: usize,
    pub max_report_bytes: usize,
}
impl Default for ProjectionLimits {
    fn default() -> Self {
        Self {
            max_pairs: 50_000,
            max_work_cells: 16_000_000,
            max_work_bytes: 32 * 1024 * 1024,
            max_mapping_segments: 8192,
            max_report_bytes: 32 * 1024 * 1024,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    /// Geometry inspection can describe a retained-state branch even when an
    /// intended reaction did not apply. Strict cloning callers can require all.
    pub require_all_requested_steps_applied: bool,
    /// Unknown end chemistry is retained by default, never changed to OH/P.
    pub require_known_terminal_chemistry: bool,
    pub limits: ProjectionLimits,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ContourTopology {
    Linear,
    Circular,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ReferenceStrand {
    Top,
    Bottom,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct EndChemistry {
    pub five_prime: FivePrime,
    pub three_prime: ThreePrime,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DuplexEnds {
    pub left: EndChemistry,
    pub right: EndChemistry,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProjectedNick {
    pub strand: ReferenceStrand,
    pub boundary: usize,
    pub three: Call,
    pub five: Call,
    pub three_prime: ThreePrime,
    pub five_prime: FivePrime,
}
/// For reference coordinate r=start+k, top is source top_start+k and bottom is
/// source bottom_start+length-1-k. Both source strands use actual5′→3′ indices.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PairedSourceMap {
    pub reference_start: usize,
    pub length: usize,
    pub top_strand: usize,
    pub top_start: usize,
    pub bottom_strand: usize,
    pub bottom_start: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DuplexView {
    pub topology: ContourTopology,
    pub top: Sequence,
    /// Actual source calls aligned antiparallel to top, stored5′→3′. A circular
    /// view may rotate a chemical chain; nicks and full source graph stay intact.
    pub bottom: Sequence,
    pub origin: Call,
    pub seed_reference_position: usize,
    pub source_maps: Vec<PairedSourceMap>,
    pub nicks: Vec<ProjectedNick>,
    pub ends: Option<DuplexEnds>,
    pub fully_known_terminal_chemistry: bool,
    pub all_requested_steps_applied: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum LimitKind {
    Pairs,
    WorkCells,
    WorkBytes,
    MappingSegments,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Withheld {
    UnpairedNucleotides {
        count: usize,
    },
    AmbiguousNucleotide {
        call: Call,
    },
    RequestedStepsNotApplied,
    SeedOutsideComponent,
    ConflictingContinuation {
        call: Call,
        forward: bool,
        first: Call,
        second: Call,
    },
    BackboneCrossesPairedTracks {
        call: Call,
        other: Call,
    },
    NonInverseContinuation {
        call: Call,
    },
    NonAntiparallelCycle {
        call: Call,
    },
    IncompleteContour {
        covered_pairs: usize,
        component_pairs: usize,
    },
    UnknownTerminalChemistry,
    Budget {
        kind: LimitKind,
        required: usize,
        limit: usize,
    },
    InconsistentGraph,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Outcome {
    Proven(DuplexView),
    Withheld(Withheld),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProjectionReport {
    pub engine: String,
    /// Domain includes exact source/step/engine-limit binding AND selection,
    /// projection policy and projection limits. No independent target is read.
    pub batch_revision: String,
    pub selection: Selection,
    pub options: Options,
    pub physical: Report,
    pub outcome: Outcome,
}
#[derive(Debug)]
pub enum Error {
    Engine(crate::Error),
    InvalidSelection(&'static str),
    InvalidLimits,
    ReportBytes { limit: usize },
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Engine(e) => e.fmt(f),
            Self::InvalidSelection(s) => write!(f, "invalid paired-contour selection: {s}"),
            Self::InvalidLimits => write!(f, "invalid paired-contour limits"),
            Self::ReportBytes { limit } => write!(
                f,
                "paired-contour evidence exceeds {limit} bytes; no view/absence claim returned"
            ),
        }
    }
}
impl std::error::Error for Error {}
impl From<crate::Error> for Error {
    fn from(e: crate::Error) -> Self {
        Self::Engine(e)
    }
}
pub(crate) fn limits(l: &ProjectionLimits) -> Result<(), Error> {
    if l.max_pairs == 0
        || l.max_pairs > 500_000
        || l.max_work_cells == 0
        || l.max_work_cells > 128_000_000
        || l.max_work_bytes == 0
        || l.max_work_bytes > 128 * 1024 * 1024
        || l.max_mapping_segments == 0
        || l.max_mapping_segments > 16384
        || l.max_report_bytes == 0
        || l.max_report_bytes > 32 * 1024 * 1024
    {
        Err(Error::InvalidLimits)
    } else {
        Ok(())
    }
}
fn selection(
    s: &Substrates,
    t: &[Step],
    l: &Limits,
    select: &Selection,
    o: &Options,
) -> Result<String, Error> {
    limits(&o.limits)?;
    if select.source_batch_revision.len() != 64
        || !select
            .source_batch_revision
            .bytes()
            .all(|b| b.is_ascii_hexdigit())
        || select.component >= 64
        || select.reference_seed.strand >= s.strands.len()
    {
        return Err(Error::InvalidSelection(
            "bounded source/component/seed identity required",
        ));
    }
    crate::preflight(s, t, l)?;
    if select.reference_seed.index >= s.strands[select.reference_seed.strand].sequence.len() {
        return Err(Error::InvalidSelection(
            "seed outside immutable source strand",
        ));
    }
    let actual = crate::batch_revision(s, t, l)?;
    if select.source_batch_revision != actual {
        return Err(Error::InvalidSelection(
            "source/transfer/engine options changed",
        ));
    }
    Ok(actual)
}
pub fn batch_revision(
    s: &Substrates,
    t: &[Step],
    l: &Limits,
    select: &Selection,
    o: &Options,
) -> Result<String, Error> {
    let source = selection(s, t, l, select, o)?;
    Ok(crate::hash(&(ENGINE, source, select, o))?)
}

struct Graph {
    offsets: Vec<usize>,
    bases: Vec<u8>,
    pairs: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
    prev: Vec<Option<usize>>,
    member: Vec<bool>,
    five: BTreeMap<usize, FivePrime>,
    three: BTreeMap<usize, ThreePrime>,
}
impl Graph {
    fn call(&self, node: usize) -> Call {
        let strand = self.offsets.partition_point(|&x| x <= node) - 1;
        Call {
            strand,
            index: node - self.offsets[strand],
        }
    }
    fn flat(&self, c: Call) -> usize {
        self.offsets[c.strand] + c.index
    }
    fn new(
        s: &Substrates,
        product: &ProductGraph,
        n: usize,
        component: usize,
    ) -> Result<Self, Withheld> {
        let mut g = Self {
            offsets: vec![],
            bases: Vec::with_capacity(n),
            pairs: vec![None; n],
            next: vec![None; n],
            prev: vec![None; n],
            member: vec![false; n],
            five: BTreeMap::new(),
            three: BTreeMap::new(),
        };
        for strand in &s.strands {
            g.offsets.push(g.bases.len());
            g.bases
                .extend_from_slice(strand.sequence.as_str().as_bytes())
        }
        for PairSpan {
            strand_a,
            start_a,
            strand_b,
            start_b,
            length,
        } in &product.source_pairing
        {
            for k in 0..*length {
                let a = g.offsets[*strand_a] + start_a + k;
                let b = g.offsets[*strand_b] + start_b + length - 1 - k;
                if g.pairs[a].is_some() || g.pairs[b].is_some() {
                    return Err(Withheld::InconsistentGraph);
                }
                g.pairs[a] = Some(b);
                g.pairs[b] = Some(a)
            }
        }
        for (id, chain) in product.chains.iter().enumerate() {
            let selected = product.components[component].chains.contains(&id);
            let mut first = None;
            let mut last: Option<usize> = None;
            for p in &chain.sources {
                for k in p.start..p.end {
                    let current = g.offsets[p.source_strand] + k;
                    if selected {
                        g.member[current] = true
                    }
                    if let Some(previous) = last {
                        if g.next[previous].is_some() || g.prev[current].is_some() {
                            return Err(Withheld::InconsistentGraph);
                        }
                        g.next[previous] = Some(current);
                        g.prev[current] = Some(previous)
                    } else {
                        first = Some(current)
                    }
                    last = Some(current)
                }
            }
            let (Some(first), Some(last)) = (first, last) else {
                return Err(Withheld::InconsistentGraph);
            };
            match &chain.ends {
                ChainEnds::CovalentlyCircular => {
                    g.next[last] = Some(first);
                    g.prev[first] = Some(last)
                }
                ChainEnds::Linear {
                    five_prime,
                    three_prime,
                } => {
                    g.five.insert(first, *five_prime);
                    g.three.insert(last, three_prime.clone());
                }
            }
        }
        Ok(g)
    }
    fn consensus(&self, node: usize, forward: bool) -> Result<Option<usize>, Withheld> {
        let partner = self.pairs[node].ok_or(Withheld::InconsistentGraph)?;
        let a = if forward {
            self.next[node]
        } else {
            self.prev[node]
        };
        let opposite = if forward {
            self.prev[partner]
        } else {
            self.next[partner]
        };
        let b = opposite
            .map(|x| self.pairs[x].ok_or(Withheld::InconsistentGraph))
            .transpose()?;
        if let (Some(first), Some(second)) = (a, b)
            && first != second
        {
            return Err(Withheld::ConflictingContinuation {
                call: self.call(node),
                forward,
                first: self.call(first),
                second: self.call(second),
            });
        }
        let result = a.or(b);
        if let Some(other) = result {
            if other == node || Some(other) == self.pairs[node] {
                return Err(Withheld::BackboneCrossesPairedTracks {
                    call: self.call(node),
                    other: self.call(other),
                });
            }
            if !self.member[other] {
                return Err(Withheld::InconsistentGraph);
            }
        }
        Ok(result)
    }
}
fn budget(required: usize, limit: usize, kind: LimitKind) -> Result<(), Withheld> {
    if required > limit {
        Err(Withheld::Budget {
            kind,
            required,
            limit,
        })
    } else {
        Ok(())
    }
}
fn known(f: FivePrime, t: &ThreePrime) -> bool {
    f != FivePrime::Unknown && !matches!(t, ThreePrime::Unknown)
}
pub(crate) fn prove(
    s: &Substrates,
    product: &ProductGraph,
    source_bases: usize,
    all_requested_applied: bool,
    component: usize,
    reference_seed: Call,
    o: &Options,
) -> Result<DuplexView, Withheld> {
    let c = &product.components[component];
    if c.unpaired_nucleotides > 0 {
        return Err(Withheld::UnpairedNucleotides {
            count: c.unpaired_nucleotides,
        });
    }
    if o.require_all_requested_steps_applied && !all_requested_applied {
        return Err(Withheld::RequestedStepsNotApplied);
    }
    let pairs = c.paired_nucleotides / 2;
    if pairs == 0 || !c.paired_nucleotides.is_multiple_of(2) {
        return Err(Withheld::InconsistentGraph);
    }
    budget(pairs, o.limits.max_pairs, LimitKind::Pairs)?;
    // Additional projection workspace, not the engine's already-bounded report.
    let work = source_bases * 64 + c.paired_nucleotides * 8;
    let bytes = source_bases * 128
        + o.limits.max_mapping_segments.min(pairs) * 128
        + s.strands.len() * 4096
        + 8192;
    budget(work, o.limits.max_work_cells, LimitKind::WorkCells)?;
    budget(bytes, o.limits.max_work_bytes, LimitKind::WorkBytes)?;
    let g = Graph::new(s, product, source_bases, component)?;
    let seed = g.flat(reference_seed);
    if !g.member[seed] {
        return Err(Withheld::SeedOutsideComponent);
    }
    for (node, &member) in g.member.iter().enumerate() {
        if member {
            if !matches!(g.bases[node], b'A' | b'C' | b'G' | b'T') {
                return Err(Withheld::AmbiguousNucleotide { call: g.call(node) });
            }
            let p = g.pairs[node].ok_or(Withheld::InconsistentGraph)?;
            if !g.member[p] || crate::comp(g.bases[node]) != Some(g.bases[p]) {
                return Err(Withheld::InconsistentGraph);
            }
        }
    }
    let mut next = vec![None; source_bases];
    let mut prev = vec![None; source_bases];
    for (node, &member) in g.member.iter().enumerate() {
        if member {
            next[node] = g.consensus(node, true)?;
            prev[node] = g.consensus(node, false)?
        }
    }
    for (node, &member) in g.member.iter().enumerate() {
        if member
            && (next[node].is_some_and(|n| prev[n] != Some(node))
                || prev[node].is_some_and(|p| next[p] != Some(node)))
        {
            return Err(Withheld::NonInverseContinuation { call: g.call(node) });
        }
    }
    // Find a physical path endpoint, or prove that predecessor traversal returns
    // to this explicit seed. No target sequence/origin is consulted.
    let mut left = seed;
    let mut topology = ContourTopology::Linear;
    let mut found = false;
    for _ in 0..=pairs {
        match prev[left] {
            None => {
                found = true;
                break;
            }
            Some(p) if p == seed => {
                topology = ContourTopology::Circular;
                left = seed;
                found = true;
                break;
            }
            Some(p) => {
                if Some(p) == g.pairs[seed] {
                    return Err(Withheld::NonAntiparallelCycle { call: g.call(p) });
                }
                left = p
            }
        }
    }
    if !found {
        return Err(Withheld::IncompleteContour {
            covered_pairs: 0,
            component_pairs: pairs,
        });
    }
    let mut visited = vec![false; source_bases];
    let mut path = Vec::with_capacity(pairs);
    let mut current = left;
    loop {
        let partner = g.pairs[current].ok_or(Withheld::InconsistentGraph)?;
        let canonical = current.min(partner);
        if visited[canonical] {
            return Err(Withheld::NonAntiparallelCycle {
                call: g.call(current),
            });
        }
        visited[canonical] = true;
        path.push(current);
        if path.len() > pairs {
            return Err(Withheld::IncompleteContour {
                covered_pairs: path.len(),
                component_pairs: pairs,
            });
        }
        match next[current] {
            None => {
                if topology == ContourTopology::Circular {
                    return Err(Withheld::InconsistentGraph);
                }
                break;
            }
            Some(n) if n == left => {
                if topology != ContourTopology::Circular {
                    return Err(Withheld::InconsistentGraph);
                }
                break;
            }
            Some(n) => current = n,
        }
    }
    if path.len() != pairs {
        return Err(Withheld::IncompleteContour {
            covered_pairs: path.len(),
            component_pairs: pairs,
        });
    }
    let seed_reference_position =
        path.iter()
            .position(|&p| p == seed)
            .ok_or(Withheld::NonAntiparallelCycle {
                call: reference_seed,
            })?;
    let mut maps: Vec<PairedSourceMap> = vec![];
    let mut top = String::with_capacity(pairs);
    let mut bottom = String::with_capacity(pairs);
    for (i, &node) in path.iter().enumerate() {
        let opposite = g.pairs[node].ok_or(Withheld::InconsistentGraph)?;
        let a = g.call(node);
        let b = g.call(opposite);
        top.push(g.bases[node] as char);
        if let Some(last) = maps.last_mut()
            && last.top_strand == a.strand
            && last.top_start + last.length == a.index
            && last.bottom_strand == b.strand
            && b.index + 1 == last.bottom_start
        {
            last.length += 1;
            last.bottom_start = b.index;
        } else {
            budget(
                maps.len() + 1,
                o.limits.max_mapping_segments,
                LimitKind::MappingSegments,
            )?;
            maps.push(PairedSourceMap {
                reference_start: i,
                length: 1,
                top_strand: a.strand,
                top_start: a.index,
                bottom_strand: b.strand,
                bottom_start: b.index,
            });
        }
    }
    for &node in path.iter().rev() {
        bottom.push(g.bases[g.pairs[node].ok_or(Withheld::InconsistentGraph)?] as char)
    }
    let mut nicks = vec![];
    let boundaries = if topology == ContourTopology::Circular {
        pairs
    } else {
        pairs - 1
    };
    for i in 0..boundaries {
        let j = (i + 1) % pairs;
        let a = path[i];
        let b = path[j];
        let pa = g.pairs[a].ok_or(Withheld::InconsistentGraph)?;
        let pb = g.pairs[b].ok_or(Withheld::InconsistentGraph)?;
        for (side, three, five, bond) in [
            (ReferenceStrand::Top, a, b, g.next[a] == Some(b)),
            (ReferenceStrand::Bottom, pb, pa, g.next[pb] == Some(pa)),
        ] {
            if !bond {
                let t = g
                    .three
                    .get(&three)
                    .ok_or(Withheld::InconsistentGraph)?
                    .clone();
                let f = *g.five.get(&five).ok_or(Withheld::InconsistentGraph)?;
                nicks.push(ProjectedNick {
                    strand: side,
                    boundary: j,
                    three: g.call(three),
                    five: g.call(five),
                    three_prime: t,
                    five_prime: f,
                });
            }
        }
    }
    let ends = if topology == ContourTopology::Linear {
        let a = path[0];
        let b = *path.last().ok_or(Withheld::InconsistentGraph)?;
        let pa = g.pairs[a].ok_or(Withheld::InconsistentGraph)?;
        let pb = g.pairs[b].ok_or(Withheld::InconsistentGraph)?;
        Some(DuplexEnds {
            left: EndChemistry {
                five_prime: *g.five.get(&a).ok_or(Withheld::InconsistentGraph)?,
                three_prime: g.three.get(&pa).ok_or(Withheld::InconsistentGraph)?.clone(),
            },
            right: EndChemistry {
                five_prime: *g.five.get(&pb).ok_or(Withheld::InconsistentGraph)?,
                three_prime: g.three.get(&b).ok_or(Withheld::InconsistentGraph)?.clone(),
            },
        })
    } else {
        None
    };
    let fully_known = ends.as_ref().is_none_or(|e| {
        known(e.left.five_prime, &e.left.three_prime)
            && known(e.right.five_prime, &e.right.three_prime)
    }) && nicks.iter().all(|n| known(n.five_prime, &n.three_prime));
    if o.require_known_terminal_chemistry && !fully_known {
        return Err(Withheld::UnknownTerminalChemistry);
    }
    Ok(DuplexView {
        topology,
        top: Sequence::dna(&top).map_err(|_| Withheld::InconsistentGraph)?,
        bottom: Sequence::dna(&bottom).map_err(|_| Withheld::InconsistentGraph)?,
        origin: g.call(left),
        seed_reference_position,
        source_maps: maps,
        nicks,
        ends,
        fully_known_terminal_chemistry: fully_known,
        all_requested_steps_applied: all_requested_applied,
    })
}
/// Fresh reconstruction only. The physical report is moved intact into every
/// returned outcome; a withheld duplex view never erases a flap or byproduct.
pub fn project(
    s: &Substrates,
    t: &[Step],
    engine_limits: &Limits,
    select: &Selection,
    o: &Options,
) -> Result<ProjectionReport, Error> {
    let fresh = crate::simulate_with_checkpoints(s, t, engine_limits, &[])?;
    project_fresh(s, t, engine_limits, select, o, fresh)
}
/// The token has a private native constructor; deserialized graphs/reports
/// cannot be supplied as incoming projection authority.
pub fn project_fresh(
    s: &Substrates,
    t: &[Step],
    engine_limits: &Limits,
    select: &Selection,
    o: &Options,
    mut fresh: crate::Executed,
) -> Result<ProjectionReport, Error> {
    let batch_revision = batch_revision(s, t, engine_limits, select, o)?;
    if fresh.report().batch_revision != crate::batch_revision(s, t, engine_limits)? {
        return Err(Error::InvalidSelection(
            "fresh native execution source differs",
        ));
    }
    fresh.reserve_work(
        crate::mul(fresh.report().source_bases, 128)?,
        crate::mul(fresh.report().source_bases, 512)?,
    )?;
    let physical = fresh.into_report();
    if select.component >= physical.product.components.len() {
        return Err(Error::InvalidSelection(
            "component absent from fresh physical report",
        ));
    }
    let seed_in_component = physical.product.components[select.component]
        .chains
        .iter()
        .any(|&chain| {
            physical.product.chains[chain].sources.iter().any(|p| {
                p.source_strand == select.reference_seed.strand
                    && p.start <= select.reference_seed.index
                    && select.reference_seed.index < p.end
            })
        });
    let outcome = if !seed_in_component {
        Outcome::Withheld(Withheld::SeedOutsideComponent)
    } else {
        match prove(
            s,
            &physical.product,
            physical.source_bases,
            physical.all_requested_steps_applied,
            select.component,
            select.reference_seed,
            o,
        ) {
            Ok(v) => Outcome::Proven(v),
            Err(reason) => Outcome::Withheld(reason),
        }
    };
    let result = ProjectionReport {
        engine: ENGINE.into(),
        batch_revision,
        selection: select.clone(),
        options: o.clone(),
        physical,
        outcome,
    };
    if crate::report_size(&result, o.limits.max_report_bytes).is_err() {
        return Err(Error::ReportBytes {
            limit: o.limits.max_report_bytes,
        });
    }
    result.validate_binding(&result.batch_revision)?;
    Ok(result)
}
/// Imported projection reports cannot establish current geometry/chemistry.
/// Replay reconstructs both the original physical reaction and this contour proof.
pub fn replay(
    s: &Substrates,
    t: &[Step],
    engine_limits: &Limits,
    select: &Selection,
    o: &Options,
    r: &ProjectionReport,
) -> Result<(), Error> {
    let fresh = project(s, t, engine_limits, select, o)?;
    if &fresh != r {
        return Err(Error::InvalidSelection(
            "projection differs from fresh source-based reconstruction",
        ));
    }
    Ok(())
}

impl ProjectionReport {
    /// Readonly envelope/contour correspondence; does not reconstruct chemistry.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        limits(&self.options.limits)?;
        let physical = &self.physical;
        physical.validate_binding(&self.selection.source_batch_revision)?;
        if self.engine != ENGINE
            || self.batch_revision != expected
            || crate::hash(&(
                ENGINE,
                &physical.batch_revision,
                &self.selection,
                &self.options,
            ))? != expected
            || self.selection.component >= physical.product.components.len()
            || self.selection.reference_seed.strand >= physical.source_revisions.len()
            || self.selection.reference_seed.index >= physical.source_bases
        {
            return Err(Error::InvalidSelection(
                "canonical projection source/context changed",
            ));
        }
        if let Outcome::Proven(view) = &self.outcome {
            let n = view.top.len();
            if n == 0
                || n > self.options.limits.max_pairs
                || n.checked_mul(2).is_none_or(|n| n > physical.source_bases)
                || view.bottom.len() != n
                || view.seed_reference_position >= n
                || view.source_maps.len() > self.options.limits.max_mapping_segments
                || view.nicks.len() > physical.product.nicks.len()
                || view.all_requested_steps_applied != physical.all_requested_steps_applied
                || (self.options.require_all_requested_steps_applied
                    && !view.all_requested_steps_applied)
                || (self.options.require_known_terminal_chemistry
                    && !view.fully_known_terminal_chemistry)
            {
                return Err(Error::InvalidSelection(
                    "canonical contour bounds/policy changed",
                ));
            }
            let mut cursor = 0usize;
            for map in &view.source_maps {
                if map.reference_start != cursor
                    || map.length == 0
                    || map.top_strand >= physical.source_revisions.len()
                    || map.bottom_strand >= physical.source_revisions.len()
                    || map
                        .top_start
                        .checked_add(map.length)
                        .is_none_or(|n| n > physical.source_bases)
                    || map
                        .bottom_start
                        .checked_add(map.length)
                        .is_none_or(|n| n > physical.source_bases)
                {
                    return Err(Error::InvalidSelection("canonical source map changed"));
                }
                cursor = cursor
                    .checked_add(map.length)
                    .ok_or(Error::InvalidSelection("source map overflow"))?;
                if cursor > n {
                    return Err(Error::InvalidSelection("source map exceeds contour"));
                }
            }
            if cursor != n {
                return Err(Error::InvalidSelection("source map coverage changed"));
            }
            for nick in &view.nicks {
                if nick.boundary > n
                    || !physical.product.nicks.iter().any(|x| {
                        x.three == nick.three
                            && x.five == nick.five
                            && x.three_prime == nick.three_prime
                            && x.five_prime == nick.five_prime
                    })
                {
                    return Err(Error::InvalidSelection(
                        "projected nick absent from actual graph",
                    ));
                }
            }
        }
        crate::report_size(self, self.options.limits.max_report_bytes).map_err(|_| {
            Error::ReportBytes {
                limit: self.options.limits.max_report_bytes,
            }
        })?;
        Ok(())
    }
}
/// Borrowed source/map validation, not a scientific replay or material constructor.
pub fn validate_report(
    s: &Substrates,
    steps: &[Step],
    l: &Limits,
    select: &Selection,
    o: &Options,
    report: &ProjectionReport,
) -> Result<(), Error> {
    report.validate_binding(&batch_revision(s, steps, l, select, o)?)?;
    if report.selection != *select || report.options != *o {
        return Err(Error::InvalidSelection("selected contour/options changed"));
    }
    crate::validate_report(s, steps, l, &report.physical)?;
    if let Outcome::Proven(view) = &report.outcome {
        for map in &view.source_maps {
            let top = s
                .strands
                .get(map.top_strand)
                .and_then(|st| {
                    st.sequence
                        .as_str()
                        .as_bytes()
                        .get(map.top_start..map.top_start.checked_add(map.length)?)
                })
                .ok_or(Error::InvalidSelection(
                    "top source map exceeds actual source",
                ))?;
            let bottom = s
                .strands
                .get(map.bottom_strand)
                .and_then(|st| {
                    st.sequence
                        .as_str()
                        .as_bytes()
                        .get(map.bottom_start..map.bottom_start.checked_add(map.length)?)
                })
                .ok_or(Error::InvalidSelection(
                    "bottom source map exceeds actual source",
                ))?;
            if view
                .top
                .as_str()
                .as_bytes()
                .get(map.reference_start..map.reference_start + map.length)
                != Some(top)
            {
                return Err(Error::InvalidSelection(
                    "top contour differs from immutable calls",
                ));
            }
            for (k, base) in bottom.iter().rev().enumerate() {
                if crate::comp(*base) != Some(top[k]) {
                    return Err(Error::InvalidSelection("paired source calls differ"));
                }
            }
            let start = view.bottom.len() - map.reference_start - map.length;
            if view
                .bottom
                .as_str()
                .as_bytes()
                .get(start..start + map.length)
                != Some(bottom)
            {
                return Err(Error::InvalidSelection(
                    "bottom contour differs from immutable calls",
                ));
            }
        }
    }
    Ok(())
}
