//! Defined two-site reciprocal reference-sequence exchange, not experimental Gateway prediction.
//! Products retain every source reference base. Seven-base overlap provenance is co-sourced,
//! and does not identify physical-strand ancestry. Caller profiles are hypotheses, not authority.
use collagene_model::bounds::ModelBudget;
pub use collagene_model::bounds::{feature_preflight, model_preflight};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE: &str = "two-site-reciprocal-reference-exchange-v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Family {
    B,
    P,
    L,
    R,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Reaction {
    BP,
    LR,
}
impl Reaction {
    fn families(self) -> [Family; 2] {
        match self {
            Self::BP => [Family::B, Family::P],
            Self::LR => [Family::L, Family::R],
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SourceRecord {
    pub id: String,
    pub url: String,
    pub raw_sha256: String,
    pub vector_name: String,
    pub vector_length: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub reverse_extraction: bool,
    pub map_url: String,
    pub map_sha256: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SiteVariant {
    pub id: String,
    pub family: Family,
    pub subtype: u8,
    pub sequence: String,
    pub overlap_start: usize,
    pub source_id: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Profile {
    pub id: String,
    pub edition: String,
    pub sources: Vec<SourceRecord>,
    pub variants: Vec<SiteVariant>,
    pub limitations: Vec<String>,
}
// Shared borrowed capacities for two-parent and multi-parent profile owners.
// Provenance relationships, bases, specificity and snapshot pins remain semantic
// checks; this guard performs no hashing, cloning or sequence expansion.
fn profile_resource_preflight(
    id: &str,
    edition: &str,
    sources: &[SourceRecord],
    variants: &[SiteVariant],
    limitations: &[String],
) -> Result<(), Error> {
    if id.len() > 256
        || edition.len() > 256
        || variants.len() > 64
        || sources.len() > 64
        || limitations.len() > 64
        || limitations.iter().any(|text| text.len() > 4096)
        || sources.iter().any(|source| {
            source.id.len() > 256
                || source.vector_name.len() > 1024
                || source.url.len() > 4096
                || source.map_url.len() > 4096
                || source.map_sha256.len() > 64
                || source.raw_sha256.len() > 64
        })
        || variants.iter().any(|variant| {
            variant.id.len() > 256 || variant.source_id.len() > 256 || variant.sequence.len() > 4096
        })
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
impl Profile {
    pub fn resource_preflight(&self) -> Result<(), Error> {
        profile_resource_preflight(
            &self.id,
            &self.edition,
            &self.sources,
            &self.variants,
            &self.limitations,
        )
    }
    pub fn validate(&self) -> Result<(), Error> {
        self.resource_preflight()?;
        if self.id.trim().is_empty() || self.edition.trim().is_empty() || self.variants.is_empty() {
            return Err(Error::TooLarge);
        }
        let mut ids = BTreeSet::new();
        let mut source_ids = BTreeSet::new();
        for s in &self.sources {
            if s.id.trim().is_empty()
                || !s.map_url.starts_with("https://")
                || s.map_sha256.len() != 64
                || !s.map_sha256.bytes().all(|b| b.is_ascii_hexdigit())
                || !source_ids.insert(s.id.as_str())
                || !s.url.starts_with("https://")
                || s.raw_sha256.len() != 64
                || !s.raw_sha256.bytes().all(|b| b.is_ascii_hexdigit())
                || s.start_1based == 0
                || s.end_1based < s.start_1based
                || s.end_1based > s.vector_length
            {
                return Err(Error::Format("invalid profile source record".into()));
            }
        }
        for v in &self.variants {
            if v.id.trim().is_empty()
                || !ids.insert(v.id.as_str())
                || !source_ids.contains(v.source_id.as_str())
                || !(1..=2).contains(&v.subtype)
                || v.sequence.len() < 7
                || !v.sequence.bytes().all(|b| b"ACGT".contains(&b))
                || v.overlap_start
                    .checked_add(7)
                    .is_none_or(|n| n > v.sequence.len())
            {
                return Err(Error::Format(
                    "invalid exact attachment-site variant".into(),
                ));
            }
            let source = self
                .sources
                .iter()
                .find(|s| s.id == v.source_id)
                .expect("source checked");
            if source.end_1based - source.start_1based + 1 != v.sequence.len() {
                return Err(Error::Format(
                    "profile source span differs from exact variant length".into(),
                ));
            }
            let overlap = &v.sequence[v.overlap_start..v.overlap_start + 7];
            if overlap != if v.subtype == 1 { "GTACAAA" } else { "GTACAAG" } {
                return Err(Error::Unsupported(
                    "this engine requires sourced att1/att2 seven-base overlaps".into(),
                ));
            }
        }
        encoded_len(self, 131072)?;
        Ok(())
    }
    pub fn revision(&self) -> Result<String, Error> {
        self.validate()?;
        hash(self)
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Input {
    pub role: String,
    pub molecule: Molecule,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Limits {
    pub max_input_bases: usize,
    pub max_comparisons: u64,
    pub max_sites: usize,
    pub max_combinations: usize,
    pub max_candidates: usize,
    pub max_components_per_candidate: usize,
    pub max_product_bases: usize,
    pub max_report_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_input_bases: 1_000_000,
            max_comparisons: 20_000_000,
            max_sites: 256,
            max_combinations: 4096,
            max_candidates: 64,
            max_components_per_candidate: 4,
            max_product_bases: 4_000_000,
            max_report_bytes: 16 * 1024 * 1024,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SiteHit {
    pub input: usize,
    pub variant_id: String,
    pub family: Family,
    pub subtype: u8,
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub certain: bool,
    /// Boundary before the overlap in canonical orientation; for reverse sites this is
    /// physically AFTER the seven-base reference overlap, not its left coordinate.
    pub switch_boundary: usize,
    pub overlap_positions: Vec<usize>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SegmentMap {
    pub output_start: usize,
    pub output_end: usize,
    pub input: usize,
    pub source_start: usize,
    pub source_end: usize,
    pub reverse: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Junction {
    pub output_boundary: usize,
    pub product_overlap_positions: Vec<usize>,
    pub subtype: u8,
    pub family: Family,
    pub reverse: bool,
    pub parents: [SiteHit; 2],
    pub matched_output_variants: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Product {
    pub molecule: Molecule,
    pub segments: Vec<SegmentMap>,
    pub junctions: Vec<Junction>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Candidate {
    pub id: usize,
    pub sites: [SiteHit; 4],
    pub products: Vec<Product>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Report {
    pub engine: String,
    pub batch_revision: String,
    pub profile_revision: String,
    pub reaction: Reaction,
    pub limits: Limits,
    pub sites: Vec<SiteHit>,
    pub candidates: Vec<Candidate>,
    pub comparisons: u64,
    pub combinations: usize,
    pub graph_plans: usize,
    pub components_withheld: usize,
    pub materializations: usize,
    pub materialized_bases: usize,
    pub unsupported_combinations: usize,
    pub complete: bool,
    pub diagnostics: Vec<String>,
    pub assumptions: Vec<String>,
}
// Counting serialization stops at the cap without allocating an untrusted JSON copy.
fn encoded_len<T: Serialize + ?Sized>(v: &T, cap: usize) -> Result<usize, Error> {
    struct Budget {
        used: usize,
        cap: usize,
    }
    impl std::io::Write for Budget {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            let end = self
                .used
                .checked_add(buf.len())
                .ok_or_else(|| std::io::Error::other("JSON limit"))?;
            if end > self.cap {
                return Err(std::io::Error::other("JSON limit"));
            }
            self.used = end;
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut writer = Budget { used: 0, cap };
    serde_json::to_writer(&mut writer, v).map_err(|_| Error::TooLarge)?;
    Ok(writer.used)
}
fn fmt(e: impl std::fmt::Display) -> Error {
    Error::Format(e.to_string())
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    Ok(blake3::hash(&serde_json::to_vec(v).map_err(fmt)?)
        .to_hex()
        .to_string())
}
fn rc(s: &str) -> String {
    Sequence::dna(s)
        .expect("validated DNA")
        .reverse_complement()
        .expect("DNA")
        .as_str()
        .to_owned()
}
fn mask(b: u8) -> u8 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
pub fn limits_preflight(l: &Limits) -> Result<(), Error> {
    if l.max_input_bases == 0
        || l.max_input_bases > 2_000_000
        || l.max_comparisons == 0
        || l.max_comparisons > 100_000_000
        || l.max_sites == 0
        || l.max_sites > 1024
        || l.max_combinations == 0
        || l.max_combinations > 65536
        || l.max_candidates == 0
        || l.max_candidates > 256
        || l.max_components_per_candidate == 0
        || l.max_components_per_candidate > 4
        || l.max_product_bases == 0
        || l.max_product_bases > 8_000_000
        || l.max_report_bytes < 1024
        || l.max_report_bytes > 32 * 1024 * 1024
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub fn preflight(inputs: &[Input], profile: &Profile, l: &Limits) -> Result<(), Error> {
    preflight_models(inputs, profile, l, None)
}
/// Validate borrowed sources and an optional independent target before copying them
/// into owned inputs. Source role identity is checked separately by `preflight`.
pub fn source_preflight(
    molecules: [&Molecule; 2],
    profile: &Profile,
    l: &Limits,
    target: Option<&Molecule>,
) -> Result<(), Error> {
    limits_preflight(l)?;
    profile.validate()?;
    let mut total = 0usize;
    for molecule in molecules {
        total = total
            .checked_add(molecule.sequence.len())
            .filter(|n| *n <= l.max_input_bases)
            .ok_or(Error::TooLarge)?;
        if molecule.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Unsupported("DNA required".into()));
        }
    }
    if let Some(t) = target {
        if t.sequence.len() > l.max_input_bases {
            return Err(Error::TooLarge);
        }
        if t.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Unsupported("DNA target required".into()));
        }
        model_preflight(&[molecules[0], molecules[1], t])?;
    } else {
        model_preflight(&molecules)?;
    }
    // All borrowed trees/strings, including the target, are bounded before recursion.
    for molecule in molecules {
        molecule.validate()?;
    }
    if let Some(t) = target {
        t.validate()?;
        encoded_len(t, 16 * 1024 * 1024)?;
    }
    encoded_len(&molecules, 16 * 1024 * 1024)?;
    Ok(())
}

fn preflight_models(
    inputs: &[Input],
    profile: &Profile,
    l: &Limits,
    target: Option<&Molecule>,
) -> Result<(), Error> {
    limits_preflight(l)?;
    profile.validate()?;
    if inputs.len() != 2 {
        return Err(Error::Unsupported(
            "exactly two source molecules are required".into(),
        ));
    }
    let mut roles = BTreeSet::new();
    for i in inputs {
        if i.role.len() > 256 || i.role.trim().is_empty() || !roles.insert(i.role.as_str()) {
            return Err(Error::Format(
                "distinct bounded source roles required".into(),
            ));
        }
    }
    source_preflight(
        [&inputs[0].molecule, &inputs[1].molecule],
        profile,
        l,
        target,
    )?;
    encoded_len(inputs, 16 * 1024 * 1024)?;
    Ok(())
}
pub fn batch_revision(
    inputs: &[Input],
    profile: &Profile,
    r: Reaction,
    l: &Limits,
) -> Result<String, Error> {
    preflight(inputs, profile, l)?;
    hash(&(
        ENGINE,
        inputs
            .iter()
            .map(|i| (&i.role, i.molecule.revision()))
            .collect::<Vec<_>>(),
        profile,
        r,
        l,
    ))
}
fn positions(start: usize, len: usize, n: usize, reverse: bool) -> Vec<usize> {
    let mut p = (0..len).map(|j| (start + j) % n).collect::<Vec<_>>();
    if reverse {
        p.reverse()
    }
    p
}
fn scan(inputs: &[Input], profile: &Profile, r: Reaction, l: &Limits) -> (Vec<SiteHit>, u64, bool) {
    let mut hits = vec![];
    let mut comparisons = 0;
    for (input, i) in inputs.iter().enumerate() {
        let dna = i.molecule.sequence.as_str().as_bytes();
        let n = dna.len();
        for v in profile
            .variants
            .iter()
            .filter(|v| r.families().contains(&v.family))
        {
            if v.sequence.len() > n {
                continue;
            }
            for reverse in [false, true] {
                let pattern = if reverse {
                    rc(&v.sequence)
                } else {
                    v.sequence.clone()
                };
                let starts = if i.molecule.topology == Topology::Circular {
                    n
                } else {
                    n - pattern.len() + 1
                };
                for start in 0..starts {
                    let mut possible = true;
                    let mut certain = true;
                    for (j, b) in pattern.bytes().enumerate() {
                        if comparisons == l.max_comparisons {
                            return (hits, comparisons, false);
                        }
                        comparisons += 1;
                        let actual = dna[(start + j) % n];
                        if mask(actual) & mask(b) == 0 {
                            possible = false;
                            break;
                        }
                        certain &= actual == b;
                    }
                    if possible {
                        if hits.len() == l.max_sites {
                            return (hits, comparisons, false);
                        }
                        let offset = if reverse {
                            v.sequence.len() - v.overlap_start - 7
                        } else {
                            v.overlap_start
                        };
                        let boundary = if reverse {
                            start + v.sequence.len() - v.overlap_start
                        } else {
                            start + v.overlap_start
                        };
                        hits.push(SiteHit {
                            input,
                            variant_id: v.id.clone(),
                            family: v.family,
                            subtype: v.subtype,
                            start,
                            length: v.sequence.len(),
                            reverse,
                            certain,
                            switch_boundary: if i.molecule.topology == Topology::Circular {
                                boundary % n
                            } else {
                                boundary
                            },
                            overlap_positions: positions(start + offset, 7, n, reverse),
                        });
                    }
                }
            }
        }
    }
    (hits, comparisons, true)
}
fn disjoint(a: &SiteHit, b: &SiteHit, n: usize) -> bool {
    // Half-open cyclic arcs are disjoint iff neither reaches the other's start.
    // This is constant-time even at the full combination cap.
    let distance = (b.start + n - a.start) % n;
    distance >= a.length && n - distance >= b.length
}
#[derive(Clone)]
struct Edge {
    input: usize,
    start: usize,
    end: usize,
    wrap: bool,
    a: usize,
    b: usize,
}
fn append_edge(
    edge: &Edge,
    reverse: bool,
    inputs: &[Input],
    sequence: &mut String,
    maps: &mut Vec<SegmentMap>,
) {
    let source = inputs[edge.input].molecule.sequence.as_str();
    let n = source.len();
    let mut parts = if edge.wrap {
        vec![(edge.start, n), (0, edge.end)]
    } else {
        vec![(edge.start, edge.end)]
    };
    if reverse {
        parts.reverse()
    }
    for (a, b) in parts {
        if a == b {
            continue;
        }
        let start = sequence.len();
        let s = &source[a..b];
        if reverse {
            sequence.push_str(&rc(s))
        } else {
            sequence.push_str(s)
        }
        maps.push(SegmentMap {
            output_start: start,
            output_end: sequence.len(),
            input: edge.input,
            source_start: a,
            source_end: b,
            reverse,
        });
    }
}
fn family(left: Family, right: Family, r: Reaction) -> Result<Family, Error> {
    match (r, left, right) {
        (Reaction::BP, Family::P, Family::B) => Ok(Family::L),
        (Reaction::BP, Family::B, Family::P) => Ok(Family::R),
        (Reaction::LR, Family::L, Family::R) => Ok(Family::P),
        (Reaction::LR, Family::R, Family::L) => Ok(Family::B),
        _ => Err(Error::Format("incompatible junction families".into())),
    }
}
fn exact_at(s: &str, start: usize, p: &str, circular: bool) -> bool {
    if p.len() > s.len() || (!circular && start + p.len() > s.len()) {
        return false;
    }
    p.bytes()
        .enumerate()
        .all(|(j, b)| s.as_bytes()[(start + j) % s.len()] == b)
}
struct Graph {
    edges: Vec<Edge>,
    joins: BTreeMap<usize, usize>,
    port_hit: BTreeMap<usize, (SiteHit, bool)>,
    node_edge: BTreeMap<usize, usize>,
    components: usize,
}
fn plan_graph(inputs: &[Input], sites: &[SiteHit; 4]) -> Result<Option<Graph>, Error> {
    // Each source contributes its two cut boundaries. Directed reference segments become
    // undirected duplex edges; recombination rewires only their ports, not their DNA.
    let mut edges = vec![];
    let mut hit_ports = BTreeMap::new();
    for input in 0..2 {
        let mut selected = sites
            .iter()
            .filter(|s| s.input == input)
            .collect::<Vec<_>>();
        selected.sort_by_key(|s| s.switch_boundary);
        if selected.len() != 2 || selected[0].switch_boundary == selected[1].switch_boundary {
            return Ok(None);
        }
        let n = inputs[input].molecule.sequence.len();
        let c0 = selected[0].switch_boundary;
        let c1 = selected[1].switch_boundary;
        for (j, h) in selected.iter().enumerate() {
            hit_ports.insert(
                (input, h.subtype),
                (input * 4 + j * 2, input * 4 + j * 2 + 1),
            );
        }
        if inputs[input].molecule.topology == Topology::Linear {
            edges.push(Edge {
                input,
                start: 0,
                end: c0,
                wrap: false,
                a: 8 + input * 2,
                b: input * 4,
            });
            edges.push(Edge {
                input,
                start: c0,
                end: c1,
                wrap: false,
                a: input * 4 + 1,
                b: input * 4 + 2,
            });
            edges.push(Edge {
                input,
                start: c1,
                end: n,
                wrap: false,
                a: input * 4 + 3,
                b: 9 + input * 2,
            });
        } else {
            edges.push(Edge {
                input,
                start: c0,
                end: c1,
                wrap: false,
                a: input * 4 + 1,
                b: input * 4 + 2,
            });
            edges.push(Edge {
                input,
                start: c1,
                end: c0,
                wrap: true,
                a: input * 4 + 3,
                b: input * 4,
            });
        }
    }
    let mut joins = BTreeMap::new();
    let mut port_hit = BTreeMap::new();
    for h in sites {
        let (left, right) = hit_ports[&(h.input, h.subtype)];
        port_hit.insert(left, (h.clone(), false));
        port_hit.insert(right, (h.clone(), true));
    }
    for subtype in 1..=2 {
        let a = sites
            .iter()
            .find(|h| h.input == 0 && h.subtype == subtype)
            .unwrap();
        let b = sites
            .iter()
            .find(|h| h.input == 1 && h.subtype == subtype)
            .unwrap();
        let (al, ar) = hit_ports[&(0, subtype)];
        let (bl, br) = hit_ports[&(1, subtype)];
        let pairs = if a.reverse == b.reverse {
            [(al, br), (ar, bl)]
        } else {
            [(al, bl), (ar, br)]
        };
        for (x, y) in pairs {
            joins.insert(x, y);
            joins.insert(y, x);
        }
    }
    let mut node_edge = BTreeMap::new();
    for (idx, e) in edges.iter().enumerate() {
        node_edge.insert(e.a, idx);
        node_edge.insert(e.b, idx);
    }
    // At most six source segments. Walk only their port indices; no source DNA is
    // copied until component and cumulative materialization budgets have passed.
    let mut visited = BTreeSet::new();
    let mut components = 0;
    let starts = node_edge
        .keys()
        .filter(|n| !joins.contains_key(n))
        .copied()
        .chain(node_edge.keys().copied())
        .collect::<Vec<_>>();
    for first in starts {
        if visited.contains(&node_edge[&first]) {
            continue;
        }
        components += 1;
        let mut current = first;
        loop {
            let index = node_edge[&current];
            if !visited.insert(index) {
                if current != first {
                    return Err(Error::Format("non-cycle graph plan".into()));
                }
                break;
            }
            let edge = &edges[index];
            let exit = if current == edge.a { edge.b } else { edge.a };
            let Some(next) = joins.get(&exit) else { break };
            current = *next;
        }
    }
    Ok(Some(Graph {
        edges,
        joins,
        port_hit,
        node_edge,
        components,
    }))
}
fn materialize(
    inputs: &[Input],
    graph: &Graph,
    r: Reaction,
    profile: &Profile,
) -> Result<Option<Vec<Product>>, Error> {
    let Graph {
        edges,
        joins,
        port_hit,
        node_edge,
        ..
    } = graph;
    let mut used = BTreeSet::new();
    let mut products = vec![];
    let starts = node_edge
        .keys()
        .filter(|n| !joins.contains_key(n))
        .copied()
        .chain(node_edge.keys().copied())
        .collect::<Vec<_>>();
    for first in starts {
        if used.contains(&node_edge[&first]) {
            continue;
        }
        let topology = if joins.contains_key(&first) {
            Topology::Circular
        } else {
            Topology::Linear
        };
        let mut current = first;
        let mut sequence = String::new();
        let mut maps = vec![];
        let mut junctions = vec![];
        loop {
            let idx = node_edge[&current];
            if !used.insert(idx) {
                if current != first {
                    return Err(Error::Format("non-cycle product graph".into()));
                }
                break;
            }
            let e = &edges[idx];
            let reverse = current == e.b;
            append_edge(e, reverse, inputs, &mut sequence, &mut maps);
            let exit = if reverse { e.a } else { e.b };
            let Some(&next) = joins.get(&exit) else { break };
            let (from, from_right) = &port_hit[&exit];
            let (to, _) = &port_hit[&next];
            let reverse_junction = *from_right != from.reverse;
            let (left, right) = if reverse_junction {
                (to.family, from.family)
            } else {
                (from.family, to.family)
            };
            junctions.push(Junction {
                output_boundary: sequence.len(),
                product_overlap_positions: vec![],
                subtype: from.subtype,
                family: family(left, right, r)?,
                reverse: reverse_junction,
                parents: [from.clone(), to.clone()],
                matched_output_variants: vec![],
            });
            current = next;
        }
        if sequence.is_empty() {
            return Ok(None);
        }
        let n = sequence.len();
        for j in &mut junctions {
            j.output_boundary %= if topology == Topology::Circular {
                n
            } else {
                n + 1
            };
            let start = if j.reverse {
                if j.output_boundary >= 7 {
                    j.output_boundary - 7
                } else if topology == Topology::Circular {
                    n + j.output_boundary - 7
                } else {
                    return Err(Error::Format(
                        "junction overlap before linear product".into(),
                    ));
                }
            } else {
                j.output_boundary
            };
            if topology == Topology::Linear && start + 7 > n {
                return Err(Error::Format(
                    "junction overlap after linear product".into(),
                ));
            }
            j.product_overlap_positions = positions(start, 7, n, j.reverse);
            for v in profile
                .variants
                .iter()
                .filter(|v| v.family == j.family && v.subtype == j.subtype)
            {
                let p = if j.reverse {
                    rc(&v.sequence)
                } else {
                    v.sequence.clone()
                };
                let off = if j.reverse {
                    v.sequence.len() - v.overlap_start - 7
                } else {
                    v.overlap_start
                };
                let candidate_start = if start >= off {
                    Some(start - off)
                } else if topology == Topology::Circular {
                    Some((n + start - off % n) % n)
                } else {
                    None
                };
                if let Some(s) = candidate_start {
                    if exact_at(&sequence, s, &p, topology == Topology::Circular) {
                        j.matched_output_variants.push(v.id.clone())
                    }
                }
            }
        }
        let molecule = Molecule::dna(
            format!("reciprocal product {}", products.len() + 1),
            &sequence,
            topology,
        )?;
        products.push(Product {
            molecule,
            segments: maps,
            junctions,
        });
    }
    Ok(Some(products))
}
/// Fresh reconstruction from exact inputs. `complete` means the bounded search of the
/// declared full-site profile completed; it is never a whole-Gateway applicability claim.
pub fn recombine(
    inputs: &[Input],
    profile: &Profile,
    reaction: Reaction,
    limits: &Limits,
) -> Result<Report, Error> {
    let revision = batch_revision(inputs, profile, reaction, limits)?;
    let (sites, comparisons, scan_complete) = scan(inputs, profile, reaction, limits);
    let mut report=Report{engine:ENGINE.into(),batch_revision:revision,profile_revision:profile.revision()?,reaction,limits:limits.clone(),sites,candidates:vec![],comparisons,combinations:0,graph_plans:0,components_withheld:0,materializations:0,materialized_bases:0,unsupported_combinations:0,complete:scan_complete,diagnostics:vec![],assumptions:vec![
        "Defined full-site-profile two-site reciprocal reference exchange; not universal Gateway applicability".into(),
        "All reciprocal components retained; no selection, host, expression, cellular repair or efficiency prediction".into(),
        "Seven-base overlap source maps are co-sourced reference provenance, not physical-strand ancestry".into(),
        "Exact full-site recognition; unprofiled, truncated or mutated arms are unsupported, not proved inactive".into(),
        "Features and producer metadata are not projected; inspect explicit source-segment maps".into()]};
    if !scan_complete {
        report
            .diagnostics
            .push("Attachment-site scan incomplete because a work or site limit was reached".into())
    }
    if report.sites.iter().any(|h| !h.certain) {
        report.diagnostics.push("IUPAC-compatible possible sites retained as evidence; no product released from those alternatives".into())
    }
    let groups = [(0, 1), (0, 2), (1, 1), (1, 2)].map(|(input, subtype)| {
        report
            .sites
            .iter()
            .filter(|h| h.input == input && h.subtype == subtype && h.certain)
            .cloned()
            .collect::<Vec<_>>()
    });
    'search: for a in &groups[0] {
        for b in &groups[1] {
            for c in &groups[2] {
                for d in &groups[3] {
                    if report.combinations == limits.max_combinations {
                        report.complete = false;
                        report
                            .diagnostics
                            .push("Site-combination work limit reached".into());
                        break 'search;
                    }
                    report.combinations += 1;
                    let families = reaction.families();
                    if a.family == c.family
                        || b.family == d.family
                        || ![a.family, b.family, c.family, d.family]
                            .iter()
                            .all(|x| families.contains(x))
                    {
                        continue;
                    }
                    // The declared two-site cassette exchange aligns both canonical
                    // site directions coherently. Odd orientation parity is a different
                    // reaction graph outside this profile, even if two paths can be drawn.
                    if (a.reverse != c.reverse) != (b.reverse != d.reverse) {
                        report.unsupported_combinations += 1;
                        if !report
                            .diagnostics
                            .iter()
                            .any(|x| x.starts_with("Unsupported orientation"))
                        {
                            report.diagnostics.push("Unsupported orientation/topology outside the defined paired-orientation exchange model".into());
                        }
                        continue;
                    }
                    if !disjoint(a, b, inputs[0].molecule.sequence.len())
                        || !disjoint(c, d, inputs[1].molecule.sequence.len())
                    {
                        report.unsupported_combinations += 1;
                        if !report
                            .diagnostics
                            .iter()
                            .any(|d| d.starts_with("Unsupported overlapping"))
                        {
                            report.diagnostics.push(
                                "Unsupported overlapping full attachment-site combination".into(),
                            );
                        }
                        continue;
                    }
                    if report.candidates.len() == limits.max_candidates {
                        report.complete = false;
                        report.diagnostics.push(
                            "Candidate limit reached; further alternatives not reconstructed"
                                .into(),
                        );
                        break 'search;
                    }
                    let total = inputs
                        .iter()
                        .map(|i| i.molecule.sequence.len())
                        .sum::<usize>();
                    let selected = [a.clone(), b.clone(), c.clone(), d.clone()];
                    report.graph_plans += 1;
                    let Some(graph) = plan_graph(inputs, &selected)? else {
                        report.unsupported_combinations += 1;
                        continue;
                    };
                    if graph.components > limits.max_components_per_candidate {
                        report.components_withheld += 1;
                        report.complete = false;
                        if !report
                            .diagnostics
                            .iter()
                            .any(|d| d.starts_with("Component limit"))
                        {
                            report.diagnostics.push("Component limit reached; graph counted before any DNA materialization; candidate withheld atomically".into());
                        }
                        continue;
                    }
                    let Some(charged) = report
                        .materialized_bases
                        .checked_add(total)
                        .filter(|n| *n <= limits.max_product_bases)
                    else {
                        report.complete = false;
                        report
                            .diagnostics
                            .push("Cumulative attempted materialization-base limit reached".into());
                        break 'search;
                    };
                    // Charge BEFORE copying, even if later report-size/shape checks withhold this attempt.
                    report.materialized_bases = charged;
                    report.materializations += 1;
                    let Some(products) = materialize(inputs, &graph, reaction, profile)? else {
                        report.unsupported_combinations += 1;
                        report.diagnostics.push("Unsupported orientation/topology: empty DNA component in declared graph".into());
                        continue;
                    };
                    for product in &products {
                        audit_product(inputs, product)?
                    }
                    if products
                        .iter()
                        .map(|p| p.molecule.sequence.len())
                        .sum::<usize>()
                        != total
                    {
                        return Err(Error::Format(
                            "global reference-length conservation failed".into(),
                        ));
                    }
                    let candidate = Candidate {
                        id: report.candidates.len(),
                        sites: selected,
                        products,
                    };
                    // Serialize the prospective candidate before installing it. Hard input/profile/site
                    // bounds constrain this temporary allocation; persistent reports never exceed the cap.
                    let added = encoded_len(&candidate, 32 * 1024 * 1024)?;
                    let current = encoded_len(&report, 32 * 1024 * 1024)?;
                    if current
                        .checked_add(added + 2)
                        .is_none_or(|x| x > limits.max_report_bytes)
                    {
                        report.complete = false;
                        report
                            .diagnostics
                            .push("Report-byte limit reached; candidate withheld".into());
                        break 'search;
                    }
                    report.candidates.push(candidate);
                }
            }
        }
    }
    if report.candidates.is_empty() {
        report.diagnostics.push("No product within this explicitly defined full-site profile and exchange model; this does not establish absence of other Gateway-compatible contexts".into())
    }
    if report
        .candidates
        .iter()
        .flat_map(|c| &c.products)
        .flat_map(|p| &p.junctions)
        .any(|j| j.matched_output_variants.is_empty())
    {
        report.diagnostics.push("At least one source-derived hybrid does not match any registered full output-site variant; missing distal DNA has not been restored".into())
    }
    // Maps and all source DNA are accounted for independently of the graph traversal.
    for candidate in &report.candidates {
        audit_conservation(inputs, candidate)?
    }
    report.validate()?;
    Ok(report)
}
/// Independent base-by-base source-map audit. This checks descriptive maps, not the
/// biological applicability of the supplied profile or a persisted positive report.
pub fn audit_product(inputs: &[Input], product: &Product) -> Result<(), Error> {
    if inputs.len() != 2 || product.molecule.sequence.len() > 2_000_000 {
        return Err(Error::TooLarge);
    }
    let mut budget = ModelBudget::default();
    budget.molecule(&product.molecule)?;
    for i in inputs {
        budget.molecule(&i.molecule)?;
    }
    product.molecule.validate()?;
    if product.molecule.sequence.alphabet() != Alphabet::Dna
        || inputs
            .iter()
            .any(|i| i.molecule.sequence.alphabet() != Alphabet::Dna)
    {
        return Err(Error::Unsupported("DNA maps required".into()));
    }
    let mut cursor = 0usize;
    let bytes = product.molecule.sequence.as_str().as_bytes();
    for m in &product.segments {
        if m.input >= inputs.len()
            || m.output_start != cursor
            || m.source_start >= m.source_end
            || m.source_end > inputs[m.input].molecule.sequence.len()
            || m.output_end.checked_sub(m.output_start) != Some(m.source_end - m.source_start)
            || m.output_end > bytes.len()
        {
            return Err(Error::Format("invalid product segment map".into()));
        }
        let source = inputs[m.input].molecule.sequence.as_str();
        let expected = if m.reverse {
            rc(&source[m.source_start..m.source_end])
        } else {
            source[m.source_start..m.source_end].into()
        };
        if &bytes[m.output_start..m.output_end] != expected.as_bytes() {
            return Err(Error::Format(
                "product does not reconstruct from its source map".into(),
            ));
        }
        cursor = m.output_end;
    }
    if cursor != bytes.len() {
        return Err(Error::Format("product map is incomplete".into()));
    }
    for j in &product.junctions {
        if j.product_overlap_positions.len() != 7
            || j.parents
                .iter()
                .any(|h| h.overlap_positions.len() != 7 || h.input >= inputs.len())
        {
            return Err(Error::Format("invalid co-sourced overlap evidence".into()));
        }
        let overlap = if j.subtype == 1 {
            b"GTACAAA"
        } else if j.subtype == 2 {
            b"GTACAAG"
        } else {
            return Err(Error::Format("unknown subtype".into()));
        };
        for (k, &position) in j.product_overlap_positions.iter().enumerate() {
            let Some(&b) = bytes.get(position) else {
                return Err(Error::Format("overlap outside product".into()));
            };
            let actual = if j.reverse {
                rc(&(b as char).to_string()).as_bytes()[0]
            } else {
                b
            };
            if actual != overlap[k] {
                return Err(Error::Format(
                    "product junction does not preserve its canonical overlap".into(),
                ));
            }
        }
        for parent in &j.parents {
            let dna = inputs[parent.input].molecule.sequence.as_str().as_bytes();
            for (k, &position) in parent.overlap_positions.iter().enumerate() {
                let Some(&b) = dna.get(position) else {
                    return Err(Error::Format("parent overlap outside input".into()));
                };
                let actual = if parent.reverse {
                    rc(&(b as char).to_string()).as_bytes()[0]
                } else {
                    b
                };
                if actual != overlap[k] {
                    return Err(Error::Format("co-source overlap does not agree".into()));
                }
            }
        }
    }
    Ok(())
}
/// Every original reference position must occur once across all products. Overlap
/// co-source evidence is not a second copy of output DNA and is excluded from this sum.
pub fn audit_conservation(inputs: &[Input], candidate: &Candidate) -> Result<(), Error> {
    if candidate.products.is_empty() || candidate.products.len() > 4 {
        return Err(Error::Format(
            "all bounded reciprocal components required".into(),
        ));
    }
    for (input, i) in inputs.iter().enumerate() {
        let mut intervals = candidate
            .products
            .iter()
            .flat_map(|p| &p.segments)
            .filter(|m| m.input == input)
            .map(|m| (m.source_start, m.source_end))
            .collect::<Vec<_>>();
        intervals.sort_unstable();
        let mut end = 0;
        for (a, b) in intervals {
            if a != end || b < a {
                return Err(Error::Format(
                    "source reference not conserved exactly once".into(),
                ));
            }
            end = b;
        }
        if end != i.molecule.sequence.len() {
            return Err(Error::Format(
                "source reference missing from reciprocal products".into(),
            ));
        }
    }
    Ok(())
}
impl Report {
    pub fn unique_under_declared_profile(&self) -> bool {
        self.complete
            && self.unsupported_combinations == 0
            && self.candidates.len() == 1
            && !self.sites.iter().any(|h| !h.certain)
            && !self
                .diagnostics
                .iter()
                .any(|s| s.starts_with("Unsupported orientation"))
    }
    /// Bounded structural check only. Imported reports do not grant current authority.
    pub fn validate(&self) -> Result<(), Error> {
        limits_preflight(&self.limits)?;
        if self.engine != ENGINE
            || self.batch_revision.len() != 64
            || self.profile_revision.len() != 64
            || self.sites.len() > self.limits.max_sites
            || self.candidates.len() > self.limits.max_candidates
            || self.comparisons > self.limits.max_comparisons
            || self.combinations > self.limits.max_combinations
            || self.graph_plans > self.combinations
            || self.components_withheld > self.graph_plans
            || self.materializations > self.graph_plans
            || self.materialized_bases > self.limits.max_product_bases
            || self.unsupported_combinations > self.combinations
        {
            return Err(Error::TooLarge);
        }
        if ![&self.batch_revision, &self.profile_revision]
            .iter()
            .all(|h| h.bytes().all(|b| b.is_ascii_hexdigit()))
        {
            return Err(Error::Format("invalid report revision".into()));
        }
        let valid_hit = |h: &SiteHit| {
            h.input < 2
                && (1..=2).contains(&h.subtype)
                && h.variant_id.len() <= 256
                && !h.variant_id.is_empty()
                && (7..=4096).contains(&h.length)
                && h.start < self.limits.max_input_bases
                && h.switch_boundary <= self.limits.max_input_bases
                && h.overlap_positions.len() == 7
                && h.overlap_positions
                    .iter()
                    .all(|&p| p < self.limits.max_input_bases)
        };
        if self.sites.iter().any(|h| !valid_hit(h)) {
            return Err(Error::Format("invalid report site shape".into()));
        }
        if self.assumptions.len() > 16 || self.diagnostics.len() > self.limits.max_combinations + 32
        {
            return Err(Error::TooLarge);
        }
        let mut text_bytes = 0usize;
        for text in self.assumptions.iter().chain(&self.diagnostics) {
            if text.len() > 4096 {
                return Err(Error::TooLarge);
            }
            text_bytes = text_bytes
                .checked_add(text.len())
                .filter(|n| *n <= 4 * 1024 * 1024)
                .ok_or(Error::TooLarge)?;
        }
        let mut shape_items = self.sites.len();
        let mut alias_items = 0usize;
        let mut aggregate = 0usize;
        let mut model_budget = ModelBudget::default();
        for (id, c) in self.candidates.iter().enumerate() {
            if c.id != id
                || c.products.is_empty()
                || c.products.len() > self.limits.max_components_per_candidate
                || c.sites
                    .iter()
                    .any(|h| !h.certain || !valid_hit(h) || !self.sites.contains(h))
            {
                return Err(Error::Format("invalid candidate shape".into()));
            }
            for p in &c.products {
                // Exact producer shape: annotation projection is deliberately absent.
                // Reject arbitrary native trees/claims before traversing or serializing them.
                if !p.molecule.features.is_empty()
                    || !p.molecule.metadata.is_empty()
                    || p.molecule.original_genbank.is_some()
                {
                    return Err(Error::TooLarge);
                }
                model_budget.molecule(&p.molecule)?;
                p.molecule.validate()?;
                let n = p.molecule.sequence.len();
                aggregate = aggregate.checked_add(n).ok_or(Error::TooLarge)?;
                if p.molecule.sequence.alphabet() != Alphabet::Dna
                    || p.segments.len() > 8
                    || p.segments.is_empty()
                    || p.junctions.is_empty()
                    || p.junctions.len() > 4
                {
                    return Err(Error::Format("invalid reciprocal product shape".into()));
                }
                let mut cursor = 0;
                for m in &p.segments {
                    if m.input >= 2
                        || m.output_start != cursor
                        || m.output_end > n
                        || m.source_end > self.limits.max_input_bases
                        || m.source_start >= m.source_end
                        || m.output_end.checked_sub(m.output_start)
                            != Some(m.source_end - m.source_start)
                    {
                        return Err(Error::Format("invalid report map shape".into()));
                    }
                    cursor = m.output_end;
                }
                if cursor != n {
                    return Err(Error::Format("report map does not cover product".into()));
                }
                shape_items = shape_items
                    .checked_add(p.segments.len() + p.junctions.len())
                    .filter(|n| *n <= 100_000)
                    .ok_or(Error::TooLarge)?;
                for j in &p.junctions {
                    alias_items = alias_items
                        .checked_add(j.matched_output_variants.len())
                        .filter(|n| *n <= 65536)
                        .ok_or(Error::TooLarge)?;
                    shape_items = shape_items
                        .checked_add(j.matched_output_variants.len())
                        .filter(|n| *n <= 100_000)
                        .ok_or(Error::TooLarge)?;
                    if j.output_boundary > n
                        || !(1..=2).contains(&j.subtype)
                        || j.product_overlap_positions.len() != 7
                        || j.product_overlap_positions.iter().any(|&p| p >= n)
                        || j.parents
                            .iter()
                            .any(|h| !valid_hit(h) || !c.sites.contains(h))
                        || j.parents[0].input == j.parents[1].input
                        || j.parents.iter().any(|h| h.subtype != j.subtype)
                        || j.matched_output_variants.len() > 64
                        || j.matched_output_variants.iter().any(|id| id.len() > 256)
                    {
                        return Err(Error::Format("invalid report junction shape".into()));
                    }
                }
            }
        }
        if aggregate > self.limits.max_product_bases || aggregate > self.materialized_bases {
            return Err(Error::TooLarge);
        }
        encoded_len(self, self.limits.max_report_bytes)?;
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        inputs: &[Input],
        profile: &Profile,
        reaction: Reaction,
        limits: &Limits,
    ) -> Result<(), Error> {
        self.validate()?;
        if self.batch_revision != batch_revision(inputs, profile, reaction, limits)?
            || self.profile_revision != profile.revision()?
            || self.reaction != reaction
            || &self.limits != limits
        {
            return Err(Error::Format("report binding mismatch".into()));
        }
        Ok(())
    }
}
pub fn replay(
    inputs: &[Input],
    profile: &Profile,
    reaction: Reaction,
    limits: &Limits,
    report: &Report,
) -> Result<(), Error> {
    report.validate_bindings_for(inputs, profile, reaction, limits)?;
    if &recombine(inputs, profile, reaction, limits)? != report {
        return Err(Error::Format(
            "report differs from fresh reconstruction".into(),
        ));
    }
    Ok(())
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TargetPolicy {
    pub allow_reverse_complement: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum TargetApplicability {
    DeclaredFullSiteProfileOnly,
    Unresolved,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TargetMatch {
    pub candidate: usize,
    pub product: usize,
    pub reverse: bool,
    pub target_origin: usize,
    /// False for a literal IUPAC representation match; it is not specified DNA.
    pub fully_specified_dna: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TargetCheck {
    pub batch_revision: String,
    pub target_revision: String,
    pub matches: Vec<TargetMatch>,
    /// Completion of search ONLY within the declared exact full-site profile/model.
    pub search_complete: bool,
    pub applicability: TargetApplicability,
    pub possible_sites: bool,
    pub unsupported_geometry: bool,
    pub ambiguous_source_dna: bool,
    pub ambiguous_target_dna: bool,
    pub no_match_proven_under_model: bool,
    pub diagnostics: Vec<String>,
    pub unique_under_declared_profile: bool,
    pub policy: TargetPolicy,
    pub warnings: Vec<String>,
}
/// Compute the exact target-check request binding without searching for products.
/// This identifies a request; it does not validate an imported answer.
pub fn check_batch_revision(
    inputs: &[Input],
    profile: &Profile,
    reaction: Reaction,
    limits: &Limits,
    target: &Molecule,
    policy: TargetPolicy,
) -> Result<String, Error> {
    preflight_models(inputs, profile, limits, Some(target))?;
    hash(&(
        "target-check-v2",
        batch_revision(inputs, profile, reaction, limits)?,
        target.revision(),
        policy,
    ))
}

/// Target is an independently supplied requirement. No target DNA is copied into products.
pub fn check_target(
    inputs: &[Input],
    profile: &Profile,
    reaction: Reaction,
    limits: &Limits,
    target: &Molecule,
    policy: TargetPolicy,
) -> Result<TargetCheck, Error> {
    preflight_models(inputs, profile, limits, Some(target))?;
    let report = recombine(inputs, profile, reaction, limits)?;
    let mut matches = vec![];
    for candidate in &report.candidates {
        for (product, p) in candidate.products.iter().enumerate() {
            if p.molecule.topology != target.topology
                || p.molecule.sequence.len() != target.sequence.len()
            {
                continue;
            }
            for reverse in [false, true] {
                if reverse && !policy.allow_reverse_complement {
                    continue;
                }
                let target_sequence = if reverse {
                    rc(target.sequence.as_str())
                } else {
                    target.sequence.as_str().to_owned()
                };
                let source = p.molecule.sequence.as_str();
                let origin = if target.topology == Topology::Circular {
                    let mut doubled = String::with_capacity(source.len() * 2);
                    doubled.push_str(source);
                    doubled.push_str(source);
                    doubled[..doubled.len() - 1].find(&target_sequence)
                } else {
                    (source == target_sequence).then_some(0)
                };
                if let Some(target_origin) = origin {
                    matches.push(TargetMatch {
                        candidate: candidate.id,
                        product,
                        reverse,
                        target_origin,
                        fully_specified_dna: p
                            .molecule
                            .sequence
                            .as_str()
                            .bytes()
                            .all(|b| b"ACGT".contains(&b))
                            && target
                                .sequence
                                .as_str()
                                .bytes()
                                .all(|b| b"ACGT".contains(&b)),
                    })
                }
            }
        }
    }
    let possible_sites = report.sites.iter().any(|h| !h.certain);
    let unsupported_geometry = report.unsupported_combinations > 0;
    let ambiguous_source_dna = inputs.iter().any(|i| {
        i.molecule
            .sequence
            .as_str()
            .bytes()
            .any(|b| !b"ACGT".contains(&b))
    });
    let ambiguous_target_dna = target
        .sequence
        .as_str()
        .bytes()
        .any(|b| !b"ACGT".contains(&b));
    let unresolved = !report.complete
        || possible_sites
        || unsupported_geometry
        || ambiguous_source_dna
        || ambiguous_target_dna;
    let no_match_proven_under_model = matches.is_empty() && !unresolved;
    let mut warnings=vec!["Sequence/topology consistency under the declared exact full-site profile is not experimental efficacy, universal Gateway applicability or biological repair verification".into()];
    if possible_sites {
        warnings.push("Possible IUPAC-compatible attachment sites were not materialized; absence of a target is unresolved".into());
    }
    if unsupported_geometry {
        warnings.push("At least one compatible site combination has unmodeled geometry; absence of a target is unresolved".into());
    }
    if ambiguous_source_dna || ambiguous_target_dna {
        warnings.push("Ambiguous DNA symbols represent sets of bases; a literal representation match does not establish fully specified DNA".into());
    }
    Ok(TargetCheck {
        batch_revision: hash(&(
            "target-check-v2",
            &report.batch_revision,
            target.revision(),
            policy,
        ))?,
        target_revision: target.revision(),
        matches,
        search_complete: report.complete,
        applicability: if unresolved {
            TargetApplicability::Unresolved
        } else {
            TargetApplicability::DeclaredFullSiteProfileOnly
        },
        possible_sites,
        unsupported_geometry,
        ambiguous_source_dna,
        ambiguous_target_dna,
        no_match_proven_under_model,
        diagnostics: report.diagnostics.clone(),
        unique_under_declared_profile: report.unique_under_declared_profile(),
        policy,
        warnings,
    })
}
/// Explicit opt-in manufacturer vector-context profile, never a universal default.
pub fn manufacturer_context_profile() -> Result<Profile, Error> {
    let p: Profile =
        serde_json::from_str(include_str!("../data/manufacturer-context-v1.json")).map_err(fmt)?;
    p.validate()?;
    Ok(p)
}

/// Bounded whole-network multi-parent reference exchange.
pub mod multisite;
