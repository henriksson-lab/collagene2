//! Bounded simultaneous multi-parent reference port exchange. All components and
//! every original reference base remain; no physical-strand or kit outcome is inferred.
use super::{
    Edge, Graph, Input, Junction, Limits, Product, Reaction, SiteHit, SiteVariant, SourceRecord,
    append_edge, disjoint, encoded_len, exact_at, family, hash, mask, positions, rc,
};
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE: &str = "simultaneous-multisite-reference-exchange-v1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Profile{pub id:String,pub edition:String,pub sources:Vec<SourceRecord>,pub variants:Vec<SiteVariant>,pub limitations:Vec<String>}}
dto! {pub struct Request{pub profile:Profile,pub reaction:Reaction,pub specificities:Vec<u8>,pub limits:Limits}}
dto! {pub struct Candidate{pub id:usize,pub sites:Vec<SiteHit>,pub pairs:Vec<[usize;2]>,pub products:Vec<Product>}}
dto! {pub struct Report{pub engine:String,pub batch_revision:String,pub profile_revision:String,pub reaction:Reaction,pub input_revisions:Vec<String>,pub complete:bool,pub comparisons:u64,pub combinations:usize,pub sites:Vec<SiteHit>,pub candidates:Vec<Candidate>,pub warnings:Vec<String>,pub limits:Vec<String>,pub unsupported:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Format(format!("Multi-site reference exchange: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit())
}
impl Profile {
    pub fn validate(&self) -> Result<(), Error> {
        if self.id.trim().is_empty()
            || self.id.len() > 256
            || self.edition.trim().is_empty()
            || self.edition.len() > 256
            || self.sources.is_empty()
            || self.sources.len() > 64
            || self.variants.is_empty()
            || self.variants.len() > 64
            || self.limitations.len() > 64
            || self.limitations.iter().any(|s| s.len() > 4096)
        {
            return Err(Error::TooLarge);
        }
        let mut sources = BTreeSet::new();
        for s in &self.sources {
            if s.id.trim().is_empty()
                || s.id.len() > 256
                || !sources.insert(s.id.as_str())
                || s.vector_name.len() > 1024
                || s.url.len() > 4096
                || s.map_url.len() > 4096
                || !s.url.starts_with("https://")
                || !s.map_url.starts_with("https://")
                || !digest(&s.raw_sha256)
                || !digest(&s.map_sha256)
                || s.start_1based == 0
                || s.end_1based < s.start_1based
                || s.end_1based > s.vector_length
            {
                return Err(bad("invalid exact full-context source record"));
            }
        }
        let mut ids = BTreeSet::new();
        let mut cores = BTreeMap::<u8, &str>::new();
        for v in &self.variants {
            if v.id.trim().is_empty()
                || v.id.len() > 256
                || !ids.insert(v.id.as_str())
                || v.source_id.len() > 256
                || !sources.contains(v.source_id.as_str())
                || !(1..=16).contains(&v.subtype)
                || v.sequence.len() < 9
                || v.sequence.len() > 4096
                || !v.sequence.bytes().all(|b| b"ACGT".contains(&b))
                || v.overlap_start == 0
                || v.overlap_start
                    .checked_add(7)
                    .is_none_or(|end| end >= v.sequence.len())
            {
                return Err(bad(
                    "full concrete site contexts with sourced flanking arms and an internal seven-base core are required",
                ));
            }
            let s = self
                .sources
                .iter()
                .find(|s| s.id == v.source_id)
                .expect("source checked");
            if s.end_1based - s.start_1based + 1 != v.sequence.len() {
                return Err(bad("source span differs from the exact full-site context"));
            }
            let core = &v.sequence[v.overlap_start..v.overlap_start + 7];
            if cores
                .insert(v.subtype, core)
                .is_some_and(|previous| previous != core)
            {
                return Err(bad(
                    "registered families/variants of one specificity disagree on their canonical seven-base core",
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
/// Finite, independently curated attachment-context facts from the cited source
/// editions. Includes eight LR contexts and the known nested full-P3 competitor;
/// not a universal kit, biochemical activity or manufacturer-endorsement model.
pub fn three_entry_source_contexts() -> Result<Profile, Error> {
    let profile: Profile =
        serde_json::from_str(include_str!("../data/three-entry-source-contexts.json"))
            .map_err(super::fmt)?;
    profile.validate()?;
    Ok(profile)
}

fn limits(l: &Limits) -> Result<(), Error> {
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
        || l.max_components_per_candidate > 32
        || l.max_product_bases == 0
        || l.max_product_bases > 8_000_000
        || l.max_report_bytes < 1024
        || l.max_report_bytes > 32 * 1024 * 1024
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub fn preflight(inputs: &[Input], r: &Request) -> Result<(), Error> {
    limits(&r.limits)?;
    r.profile.validate()?;
    if !(2..=5).contains(&inputs.len()) || r.specificities.is_empty() || r.specificities.len() > 16
    {
        return Err(bad(
            "two to five parents and one to sixteen explicit specificities are required",
        ));
    }
    let mut specificities = BTreeSet::new();
    for s in &r.specificities {
        if !(1..=16).contains(s) || !specificities.insert(*s) {
            return Err(bad("requested specificities must be distinct IDs1..16"));
        }
        for f in r.reaction.families() {
            if !r
                .profile
                .variants
                .iter()
                .any(|v| v.subtype == *s && v.family == f)
            {
                return Err(bad(
                    "requested specificity lacks a registered cognate family context",
                ));
            }
        }
    }
    let mut roles = BTreeSet::new();
    let mut bases = 0usize;
    for i in inputs {
        if i.role.trim().is_empty()
            || i.role.len() > 256
            || !roles.insert(i.role.as_str())
            || i.molecule.sequence.alphabet() != Alphabet::Dna
            || i.molecule.sequence.is_empty()
        {
            return Err(bad("distinct bounded nonempty DNA parents required"));
        }
        bases = bases
            .checked_add(i.molecule.sequence.len())
            .filter(|n| *n <= r.limits.max_input_bases)
            .ok_or(Error::TooLarge)?;
    }
    super::model_preflight(&inputs.iter().map(|i| &i.molecule).collect::<Vec<_>>())?;
    for i in inputs {
        i.molecule.validate()?;
    }
    encoded_len(inputs, 16 * 1024 * 1024)?;
    Ok(())
}
pub fn revision(inputs: &[Input], r: &Request) -> Result<String, Error> {
    preflight(inputs, r)?;
    hash(&(
        ENGINE,
        inputs
            .iter()
            .map(|i| (&i.role, i.molecule.revision()))
            .collect::<Vec<_>>(),
        r,
    ))
}
fn scan(inputs: &[Input], r: &Request) -> (Vec<SiteHit>, u64, bool) {
    let mut hits = Vec::new();
    let mut comparisons = 0;
    // Classify every registered family before reaction filtering. A short derived
    // L context can be literally embedded in an unrecombined P context.
    for (input, i) in inputs.iter().enumerate() {
        let dna = i.molecule.sequence.as_str().as_bytes();
        let n = dna.len();
        for v in &r.profile.variants {
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
                        if comparisons == r.limits.max_comparisons {
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
                        if hits.len() == r.limits.max_sites {
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
fn graph(inputs: &[Input], sites: &[SiteHit], pairs: &[[usize; 2]]) -> Result<Graph, Error> {
    let mut edges = Vec::new();
    let mut joins = BTreeMap::new();
    let mut port_hit = BTreeMap::new();
    let terminals = sites.len() * 2;
    for (idx, h) in sites.iter().enumerate() {
        port_hit.insert(idx * 2, (h.clone(), false));
        port_hit.insert(idx * 2 + 1, (h.clone(), true));
    }
    for (input, parent) in inputs.iter().enumerate() {
        let n = parent.molecule.sequence.len();
        let mut selected = sites
            .iter()
            .enumerate()
            .filter(|(_, h)| h.input == input)
            .collect::<Vec<_>>();
        selected.sort_by_key(|(_, h)| h.switch_boundary);
        let end_left = terminals + input * 2;
        let end_right = end_left + 1;
        if selected.is_empty() {
            edges.push(Edge {
                input,
                start: 0,
                end: n,
                wrap: false,
                a: end_left,
                b: end_right,
            });
            if parent.molecule.topology == Topology::Circular {
                joins.insert(end_left, end_right);
                joins.insert(end_right, end_left);
            }
            continue;
        }
        if parent.molecule.topology == Topology::Linear {
            edges.push(Edge {
                input,
                start: 0,
                end: selected[0].1.switch_boundary,
                wrap: false,
                a: end_left,
                b: selected[0].0 * 2,
            });
        }
        for adjacent in selected.windows(2) {
            edges.push(Edge {
                input,
                start: adjacent[0].1.switch_boundary,
                end: adjacent[1].1.switch_boundary,
                wrap: false,
                a: adjacent[0].0 * 2 + 1,
                b: adjacent[1].0 * 2,
            });
        }
        let first = selected[0];
        let last = selected[selected.len() - 1];
        if parent.molecule.topology == Topology::Circular {
            edges.push(Edge {
                input,
                start: last.1.switch_boundary,
                end: first.1.switch_boundary,
                wrap: true,
                a: last.0 * 2 + 1,
                b: first.0 * 2,
            });
        } else {
            edges.push(Edge {
                input,
                start: last.1.switch_boundary,
                end: n,
                wrap: false,
                a: last.0 * 2 + 1,
                b: end_right,
            });
        }
    }
    for [a, b] in pairs {
        let x = &sites[*a];
        let y = &sites[*b];
        let (al, ar) = (*a * 2, *a * 2 + 1);
        let (bl, br) = (*b * 2, *b * 2 + 1);
        let reconnect = if x.reverse == y.reverse {
            [(al, br), (ar, bl)]
        } else {
            [(al, bl), (ar, br)]
        };
        for (a, b) in reconnect {
            if joins.insert(a, b).is_some() || joins.insert(b, a).is_some() {
                return Err(bad("site port reused by more than one exchange"));
            }
        }
    }
    let mut node_edge = BTreeMap::new();
    for (i, e) in edges.iter().enumerate() {
        if node_edge.insert(e.a, i).is_some() || node_edge.insert(e.b, i).is_some() {
            return Err(bad("source port assigned twice"));
        }
    }
    let mut visited = BTreeSet::new();
    let mut components = 0;
    let starts = node_edge
        .keys()
        .filter(|p| !joins.contains_key(p))
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
            let idx = node_edge[&current];
            if !visited.insert(idx) {
                if current != first {
                    return Err(bad("non-cycle multi-site graph"));
                }
                break;
            }
            let e = &edges[idx];
            let exit = if current == e.a { e.b } else { e.a };
            let Some(next) = joins.get(&exit) else {
                break;
            };
            current = *next;
        }
    }
    Ok(Graph {
        edges,
        joins,
        port_hit,
        node_edge,
        components,
    })
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
            if !port_hit.contains_key(&exit) {
                current = next;
                continue;
            }
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
// Exact declared-context classification, not a universal family-recognition model.
fn classified(inputs: &[Input], hits: &[SiteHit]) -> (Vec<usize>, bool, bool) {
    // Group once by physical overlap. Comparing only group members avoids a
    // cubic all-hit scan and repeated seven-position heap allocations for aliases.
    let mut groups = BTreeMap::<(usize, u8, [usize; 7]), Vec<usize>>::new();
    for (index, hit) in hits.iter().enumerate() {
        let mut core = [0; 7];
        core.copy_from_slice(&hit.overlap_positions);
        core.sort_unstable();
        groups
            .entry((hit.input, hit.subtype, core))
            .or_default()
            .push(index);
    }
    let contains = |a: &SiteHit, b: &SiteHit| {
        let n = inputs[a.input].molecule.sequence.len();
        a.length > b.length && (b.start + n - a.start) % n + b.length <= a.length
    };
    let mut eligible = Vec::new();
    let mut ambiguous = false;
    let mut aliases = false;
    for group in groups.values() {
        let maximal = group
            .iter()
            .copied()
            .filter(|index| {
                !group
                    .iter()
                    .any(|outer| hits[*outer].certain && contains(&hits[*outer], &hits[*index]))
            })
            .collect::<Vec<_>>();
        let Some(first) = maximal.iter().copied().find(|i| hits[*i].certain) else {
            continue;
        };
        let first_hit = &hits[first];
        if maximal.iter().any(|i| {
            let h = &hits[*i];
            !h.certain || h.family != first_hit.family || h.reverse != first_hit.reverse
        }) {
            ambiguous = true;
            continue;
        }
        let representative = maximal
            .iter()
            .copied()
            .min_by(|a, b| {
                hits[*a]
                    .variant_id
                    .cmp(&hits[*b].variant_id)
                    .then(hits[*a].start.cmp(&hits[*b].start))
            })
            .expect("certain maximal context present");
        aliases |= maximal.len() > 1;
        eligible.push(representative);
    }
    // Preserve scan ordering and candidate IDs independently of map key order.
    eligible.sort_unstable();
    (eligible, ambiguous, aliases)
}
/// Independently reconstruct every mapped reference base and its co-sourced overlap.
/// This audits reference bookkeeping, not biochemical applicability or physical atoms.
pub fn audit_candidate(inputs: &[Input], profile: &Profile, c: &Candidate) -> Result<(), Error> {
    if !(2..=5).contains(&inputs.len())
        || c.products.is_empty()
        || c.products.len() > 32
        || c.sites.len() > 32
    {
        return Err(Error::TooLarge);
    }
    // Source cuts create at most one edge per selected site plus two per
    // source; splitting wrapping edges adds at most one map per source. The
    // loose 74-map bound covers all 32-site/five-parent graphs before allocation.
    if c.products.iter().any(|p| {
        p.segments.is_empty()
            || p.segments.len() > 74
            || p.junctions.len() > 32
            || !p.molecule.features.is_empty()
            || !p.molecule.metadata.is_empty()
            || p.molecule.original_genbank.is_some()
            || p.molecule.name.len() > 256
            || p.molecule.sequence.is_empty()
            || p.molecule.sequence.alphabet() != Alphabet::Dna
    }) {
        return Err(Error::TooLarge);
    }
    let mut input_bases = 0usize;
    for input in inputs {
        input_bases = input_bases
            .checked_add(input.molecule.sequence.len())
            .filter(|n| *n <= 2_000_000)
            .ok_or(Error::TooLarge)?;
        if input.molecule.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("DNA source maps required"));
        }
    }
    super::model_preflight(&inputs.iter().map(|i| &i.molecule).collect::<Vec<_>>())?;
    for input in inputs {
        input.molecule.validate()?;
    }
    profile.validate()?;
    let mut total = 0usize;
    for p in &c.products {
        let bytes = p.molecule.sequence.as_str().as_bytes();
        total = total
            .checked_add(bytes.len())
            .filter(|n| *n <= 2_000_000)
            .ok_or(Error::TooLarge)?;
        p.molecule.validate()?;
        let mut cursor = 0;
        for m in &p.segments {
            if m.input >= inputs.len()
                || m.output_start != cursor
                || m.source_start >= m.source_end
                || m.source_end > inputs[m.input].molecule.sequence.len()
                || m.output_end > bytes.len()
                || m.output_end.checked_sub(m.output_start) != Some(m.source_end - m.source_start)
            {
                return Err(bad("invalid source segment"));
            }
            let s = &inputs[m.input].molecule.sequence.as_str()[m.source_start..m.source_end];
            let expected = if m.reverse { rc(s) } else { s.to_string() };
            if bytes[m.output_start..m.output_end] != *expected.as_bytes() {
                return Err(bad("source segment does not reconstruct product"));
            }
            cursor = m.output_end;
        }
        if cursor != bytes.len() {
            return Err(bad("source map leaves unmapped output"));
        }
        for j in &p.junctions {
            let v = profile
                .variants
                .iter()
                .find(|v| v.subtype == j.subtype)
                .ok_or_else(|| bad("unknown specificity"))?;
            let core = &v.sequence.as_bytes()[v.overlap_start..v.overlap_start + 7];
            if j.product_overlap_positions.len() != 7
                || j.parents
                    .iter()
                    .any(|h| h.input >= inputs.len() || h.overlap_positions.len() != 7)
            {
                return Err(bad("invalid overlap evidence"));
            }
            for (k, position) in j.product_overlap_positions.iter().enumerate() {
                let b = *bytes
                    .get(*position)
                    .ok_or_else(|| bad("output overlap outside product"))?;
                let b = if j.reverse { complement(b) } else { b };
                if b != core[k] {
                    return Err(bad("output core differs from declared specificity"));
                }
            }
            for h in &j.parents {
                for (k, position) in h.overlap_positions.iter().enumerate() {
                    let b = *inputs[h.input]
                        .molecule
                        .sequence
                        .as_str()
                        .as_bytes()
                        .get(*position)
                        .ok_or_else(|| bad("source overlap outside parent"))?;
                    let b = if h.reverse { complement(b) } else { b };
                    if b != core[k] {
                        return Err(bad("co-sourced overlap disagrees"));
                    }
                }
            }
        }
    }
    for (input, parent) in inputs.iter().enumerate() {
        let mut intervals = c
            .products
            .iter()
            .flat_map(|p| &p.segments)
            .filter(|s| s.input == input)
            .map(|s| (s.source_start, s.source_end))
            .collect::<Vec<_>>();
        intervals.sort_unstable();
        let mut cursor = 0;
        for (a, b) in intervals {
            if a != cursor || b <= a {
                return Err(bad("source reference repeated or lost"));
            }
            cursor = b;
        }
        if cursor != parent.molecule.sequence.len() {
            return Err(bad("source reference not conserved exactly once"));
        }
    }
    Ok(())
}
struct Search<'a> {
    inputs: &'a [Input],
    request: &'a Request,
    groups: Vec<(Vec<usize>, Vec<usize>)>,
    report: Report,
    spent_bases: usize,
    stop: bool,
}
impl Search<'_> {
    fn halt(&mut self, note: &str) {
        self.report.complete = false;
        if !self.report.limits.iter().any(|s| s == note) {
            self.report.limits.push(note.into());
        }
        self.stop = true;
    }
    fn visit(&mut self, depth: usize, chosen: &mut Vec<usize>) -> Result<(), Error> {
        if self.stop {
            return Ok(());
        }
        if depth == self.groups.len() {
            if self.report.candidates.len() == self.request.limits.max_candidates {
                self.halt("Candidate limit reached; further combinations were not retained");
                return Ok(());
            }
            let sites = chosen
                .iter()
                .map(|i| self.report.sites[*i].clone())
                .collect::<Vec<_>>();
            let pairs = (0..sites.len() / 2)
                .map(|i| [i * 2, i * 2 + 1])
                .collect::<Vec<_>>();
            let graph = graph(self.inputs, &sites, &pairs)?;
            if graph.components > self.request.limits.max_components_per_candidate {
                self.report.complete = false;
                if self.report.limits.is_empty() {
                    self.report.limits.push("An entire candidate was withheld because its component count exceeds the declared limit".into());
                }
                return Ok(());
            }
            let bases = self
                .inputs
                .iter()
                .map(|i| i.molecule.sequence.len())
                .sum::<usize>();
            if self
                .spent_bases
                .checked_add(bases)
                .is_none_or(|n| n > self.request.limits.max_product_bases)
            {
                self.halt("Aggregate materialized-reference-base limit reached before copying candidate DNA");
                return Ok(());
            }
            self.spent_bases += bases;
            let Some(products) = materialize(
                self.inputs,
                &graph,
                self.request.reaction,
                &self.request.profile,
            )?
            else {
                self.report.complete = false;
                if !self
                    .report
                    .unsupported
                    .iter()
                    .any(|s| s.starts_with("A selected exchange creates"))
                {
                    self.report.unsupported.push("A selected exchange creates an empty reference component; the whole candidate is withheld".into());
                }
                return Ok(());
            };
            let candidate = Candidate {
                id: self.report.candidates.len(),
                sites,
                pairs,
                products,
            };
            audit_candidate(self.inputs, &self.request.profile, &candidate)?;
            // Reserve room for the finite limit/warning diagnostics before installation.
            let Ok(candidate_bytes) = encoded_len(&candidate, self.request.limits.max_report_bytes)
            else {
                self.halt("Report-byte limit reached; whole candidate withheld");
                return Ok(());
            };
            if encoded_len(&self.report, self.request.limits.max_report_bytes)?
                .checked_add(candidate_bytes + 1024)
                .is_none_or(|n| n > self.request.limits.max_report_bytes)
            {
                self.halt("Report-byte limit reached; whole candidate withheld");
                return Ok(());
            }
            self.report.candidates.push(candidate);
            return Ok(());
        }
        let (left, right) = &self.groups[depth];
        let left_len = left.len();
        let right_len = right.len();
        for a in 0..left_len {
            for b in 0..right_len {
                if self.stop {
                    return Ok(());
                }
                if self.report.combinations == self.request.limits.max_combinations {
                    self.halt("Pair-binding search-state limit reached (partial and complete assignments counted)");
                    return Ok(());
                }
                self.report.combinations += 1;
                let x = self.groups[depth].0[a];
                let y = self.groups[depth].1[b];
                let hx = &self.report.sites[x];
                let hy = &self.report.sites[y];
                if hx.input == hy.input
                    && !disjoint(hx, hy, self.inputs[hx.input].molecule.sequence.len())
                {
                    continue;
                }
                if chosen.iter().any(|i| {
                    let h = &self.report.sites[*i];
                    [hx, hy].iter().any(|other| {
                        h.input == other.input
                            && !disjoint(h, other, self.inputs[h.input].molecule.sequence.len())
                    })
                }) {
                    continue;
                }
                chosen.extend([x, y]);
                self.visit(depth + 1, chosen)?;
                chosen.truncate(chosen.len() - 2);
            }
        }
        Ok(())
    }
}
/// Reconstruct all explicitly selected cognate exchanges in one whole-reference graph.
/// No target, selection marker, host repair, expression or material-yield model is used.
pub fn recombine(inputs: &[Input], request: &Request) -> Result<Report, Error> {
    let batch_revision = revision(inputs, request)?;
    let (sites, comparisons, complete) = scan(inputs, request);
    let mut report=Report{engine:ENGINE.into(),batch_revision,profile_revision:request.profile.revision()?,reaction:request.reaction,input_revisions:inputs.iter().map(|i|i.molecule.revision()).collect(),complete,comparisons,combinations:0,sites,candidates:vec![],warnings:vec![
        "Declared exact-context simultaneous reference exchange only; not a physical-atom, commercial-kit, host-repair, selection or efficiency prediction".into(),
        "All reciprocal components, including unreacted parents, are retained; seven-base overlap co-sourcing is reference provenance, not extra output DNA".into(),
        "Family classification uses maximal certain full physical contexts at one canonical overlap; strict contained contexts are suppressed, equal/incomparable conflicting maxima remain unsupported".into(),
        "combinations counts attempted pair-binding search states, including partial assignments; complete refers only to the declared bounded reference model".into(),
        "Feature projection and biological applicability require separate source-bound review".into()],limits:vec![],unsupported:vec![]};
    if !complete {
        report.limits.push("Global all-family scan incomplete; no candidate enumeration is safe because unseen larger contexts may change classification".into());
        report.validate_binding(&report.batch_revision)?;
        encoded_len(&report, request.limits.max_report_bytes)?;
        return Ok(report);
    }
    let (eligible, ambiguous, aliases) = classified(inputs, &report.sites);
    if ambiguous {
        report.complete = false;
        report.unsupported.push("Conflicting or uncertain maximal full contexts at a canonical overlap were excluded from reaction selection".into());
    }
    if aliases {
        report.warnings.push("Equivalent same-family reactive loci are enumerated once using the lexically first maximal variant ID; all scanned context aliases remain visible".into());
    }
    if report.sites.iter().any(|h| !h.certain) {
        report.complete = false;
        report.unsupported.push("Possible attachment sites containing ambiguous source DNA are not certain reactive bindings".into());
    }
    let families = request.reaction.families();
    let groups = request
        .specificities
        .iter()
        .map(|s| {
            let group = |family| {
                eligible
                    .iter()
                    .copied()
                    .filter(|i| {
                        let h = &report.sites[*i];
                        h.subtype == *s && h.family == family
                    })
                    .collect::<Vec<_>>()
            };
            (group(families[0]), group(families[1]))
        })
        .collect::<Vec<_>>();
    if groups.iter().any(|(a, b)| a.is_empty() || b.is_empty()) {
        report.unsupported.push("At least one requested specificity has no certain unambiguous cognate family pair in these parents; no complete exchange combination was invented".into());
        report.validate_binding(&report.batch_revision)?;
        encoded_len(&report, request.limits.max_report_bytes)?;
        return Ok(report);
    }
    let mut search = Search {
        inputs,
        request,
        groups,
        report,
        spent_bases: 0,
        stop: false,
    };
    search.visit(0, &mut vec![])?;
    if search.report.candidates.is_empty() {
        search.report.warnings.push("No candidate under this declared profile and selected specificities does not prove absence of another recombination context or outcome".into());
    }
    if search
        .report
        .candidates
        .iter()
        .flat_map(|c| &c.products)
        .flat_map(|p| &p.junctions)
        .any(|j| j.matched_output_variants.is_empty())
    {
        search.report.warnings.push("At least one source-derived output hybrid does not match a registered full output context; missing distal sequence was not synthesized".into());
    }
    search
        .report
        .validate_binding(&search.report.batch_revision)?;
    encoded_len(&search.report, request.limits.max_report_bytes)?;
    Ok(search.report)
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b,
    }
}
impl Report {
    /// Bounded historical shape/binding check, never a substitute for fresh recombination.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine != ENGINE
            || self.batch_revision != expected
            || !digest(expected)
            || !digest(&self.profile_revision)
            || !(2..=5).contains(&self.input_revisions.len())
            || self.input_revisions.iter().any(|s| !digest(s))
            || self.sites.len() > 1024
            || self.candidates.len() > 256
            || self.comparisons > 100_000_000
            || self.combinations > 65536
        {
            return Err(bad("invalid report header or stale expected context"));
        }
        for notes in [&self.warnings, &self.limits, &self.unsupported] {
            if notes.len() > 256 || notes.iter().any(|s| s.len() > 4096) {
                return Err(Error::TooLarge);
            }
        }
        let valid_hit = |h: &SiteHit| {
            h.input < self.input_revisions.len()
                && (1..=16).contains(&h.subtype)
                && !h.variant_id.is_empty()
                && h.variant_id.len() <= 256
                && (9..=4096).contains(&h.length)
                && h.start < 2_000_000
                && h.switch_boundary <= 2_000_000
                && h.overlap_positions.len() == 7
                && h.overlap_positions.iter().all(|p| *p < 2_000_000)
        };
        if self.sites.iter().any(|h| !valid_hit(h)) {
            return Err(bad("invalid scanned site shape"));
        }
        let mut bases = 0usize;
        let mut items = self.sites.len();
        for (id, c) in self.candidates.iter().enumerate() {
            if c.id != id
                || c.sites.is_empty()
                || c.sites.len() > 32
                || c.sites.len() != c.pairs.len() * 2
                || c.pairs.len() > 16
                || c.products.is_empty()
                || c.products.len() > 32
                || c.sites
                    .iter()
                    .any(|h| !h.certain || !valid_hit(h) || !self.sites.contains(h))
            {
                return Err(bad("invalid candidate shape"));
            }
            let mut indices = BTreeSet::new();
            let mut specificities = BTreeSet::new();
            for pair in &c.pairs {
                if pair
                    .iter()
                    .any(|i| *i >= c.sites.len() || !indices.insert(*i))
                {
                    return Err(bad("candidate reuses a site"));
                }
                let (a, b) = (&c.sites[pair[0]], &c.sites[pair[1]]);
                let families = self.reaction.families();
                if a.subtype != b.subtype
                    || a.family != families[0]
                    || b.family != families[1]
                    || !specificities.insert(a.subtype)
                {
                    return Err(bad("noncognate or repeated specificity pair"));
                }
            }
            let mut junctions = 0usize;
            for p in &c.products {
                if !p.molecule.features.is_empty()
                    || !p.molecule.metadata.is_empty()
                    || p.molecule.original_genbank.is_some()
                    || p.molecule.sequence.is_empty()
                    || p.molecule.sequence.alphabet() != Alphabet::Dna
                    || p.molecule.name.len() > 256
                {
                    return Err(bad("unexpected product model shape"));
                }
                let n = p.molecule.sequence.len();
                bases = bases
                    .checked_add(n)
                    .filter(|n| *n <= 8_000_000)
                    .ok_or(Error::TooLarge)?;
                if p.segments.is_empty() || p.segments.len() > 74 || p.junctions.len() > 32 {
                    return Err(Error::TooLarge);
                }
                p.molecule.validate()?;
                items = items
                    .checked_add(p.segments.len() + p.junctions.len())
                    .filter(|n| *n <= 100_000)
                    .ok_or(Error::TooLarge)?;
                let mut cursor = 0;
                for m in &p.segments {
                    if m.input >= self.input_revisions.len()
                        || m.output_start != cursor
                        || m.output_end > n
                        || m.source_end > 2_000_000
                        || m.source_start >= m.source_end
                        || m.output_end.checked_sub(m.output_start)
                            != Some(m.source_end - m.source_start)
                    {
                        return Err(bad("invalid reference map shape"));
                    }
                    cursor = m.output_end;
                }
                if cursor != n {
                    return Err(bad("reference maps do not cover product"));
                }
                junctions += p.junctions.len();
                for j in &p.junctions {
                    if j.output_boundary > n
                        || !(1..=16).contains(&j.subtype)
                        || j.product_overlap_positions.len() != 7
                        || j.product_overlap_positions.iter().any(|i| *i >= n)
                        || j.parents.iter().any(|h| {
                            !valid_hit(h) || !c.sites.contains(h) || h.subtype != j.subtype
                        })
                        || j.matched_output_variants.len() > 64
                        || j.matched_output_variants
                            .iter()
                            .any(|s| s.is_empty() || s.len() > 256)
                    {
                        return Err(bad("invalid hybrid junction shape"));
                    }
                    let pair = c
                        .pairs
                        .iter()
                        .find(|pair| {
                            let a = &c.sites[pair[0]];
                            let b = &c.sites[pair[1]];
                            j.parents == [a.clone(), b.clone()]
                                || j.parents == [b.clone(), a.clone()]
                        })
                        .ok_or_else(|| bad("junction is not a selected cognate pair"))?;
                    let _ = pair;
                    items = items
                        .checked_add(j.matched_output_variants.len())
                        .filter(|n| *n <= 100_000)
                        .ok_or(Error::TooLarge)?;
                }
            }
            if junctions != c.pairs.len() * 2 {
                return Err(bad("not every reciprocal hybrid junction was retained"));
            }
        }
        encoded_len(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}
/// A historical report is not authority; replay reconstructs from source and request.
pub fn replay(inputs: &[Input], request: &Request, report: &Report) -> Result<(), Error> {
    report.validate_binding(&revision(inputs, request)?)?;
    if &recombine(inputs, request)? != report {
        return Err(bad("fresh replay differs from historical evidence"));
    }
    Ok(())
}
