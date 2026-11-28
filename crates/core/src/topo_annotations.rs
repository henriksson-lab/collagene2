//! Annotation correspondence through fresh TOPO preparation and paired-contour
//! top-call provenance. This is not chemical ancestry or biological validation.
use crate::{
    annotations::{self, MappingSegment},
    topo::{
        self, Origin, PreparationEvidence,
        physical::projection::{ContourTopology, DuplexView},
    },
};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

const MAX_MAPS: usize = 4096;
const MAX_VIEW_MAPS: usize = 8192;
const MAX_WORK: usize = 1_000_000;
const MAX_COPY_BYTES: usize = 32 * 1024 * 1024;
const MAX_FEATURE_BYTES: usize = 4 * 1024 * 1024;
const MAX_OUTPUT_BYTES: usize = 16 * 1024 * 1024;
const MAX_WARNINGS: usize = 64;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct UnmappedAnnotation {pub source_index:usize,pub feature:Feature,pub reason:String}}
dto! {pub struct RoleEvidence {
    pub role:String,pub source_revision:String,
    /// Composed candidates. Do not apply when projection_supported is false.
    pub maps:Vec<MappingSegment>,pub projection_supported:bool,pub projected:usize,
    pub unmapped:Vec<UnmappedAnnotation>,pub dropped:usize,pub warnings:Vec<String>,
}}
dto! {pub struct AnnotationEvidence {pub sources:Vec<RoleEvidence>,pub features:Vec<Feature>,pub warnings:Vec<String>}}

fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO annotation projection: {s}"))
}
fn end(start: usize, length: usize, cap: usize) -> Result<usize, Error> {
    start
        .checked_add(length)
        .filter(|e| *e <= cap)
        .ok_or_else(|| bad("range outside source/call coordinates"))
}
fn charge(work: &mut usize, amount: usize) -> Result<(), Error> {
    *work = work
        .checked_add(amount)
        .filter(|v| *v <= MAX_WORK)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
fn warning(out: &mut Vec<String>, text: &str) {
    if out.len() < MAX_WARNINGS {
        out.push(text.chars().take(1024).collect());
    } else if out.len() == MAX_WARNINGS {
        out.push("Further annotation warnings omitted; inspect source annotations and retained evidence.".into());
    }
}
struct Count {
    bytes: usize,
    cap: usize,
}
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= self.cap)
            .ok_or_else(|| std::io::Error::other("bounded TOPO annotations"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(value: &impl Serialize, cap: usize) -> Result<usize, Error> {
    let mut c = Count { bytes: 0, cap };
    serde_json::to_writer(&mut c, value).map_err(|_| Error::TooLarge)?;
    Ok(c.bytes)
}
fn nodes(location: &Location) -> usize {
    match location {
        Location::Complement(l) => 1 + nodes(l),
        Location::Join(v) | Location::Order(v) => 1 + v.iter().map(nodes).sum::<usize>(),
        _ => 1,
    }
}
fn unresolved(location: &Location) -> bool {
    match location {
        Location::Unresolved(_) => true,
        Location::Complement(l) => unresolved(l),
        Location::Join(v) | Location::Order(v) => v.iter().any(unresolved),
        _ => false,
    }
}

fn preflight(
    inputs: &[topo::Input],
    p: &PreparationEvidence,
    v: &DuplexView,
    d: &Molecule,
) -> Result<(), Error> {
    if inputs.is_empty()
        || inputs.len() > 32
        || p.substrates.strands.len() > 32
        || p.maps.len() > MAX_MAPS
        || p.substrates.pairing.len() > 1024
        || p.not_represented.len() > 4096
        || p.assumptions.len() > 256
        || v.source_maps.len() > MAX_VIEW_MAPS
        || v.top.len() > 100_000
        || v.nicks.len() > 100_000
    {
        return Err(Error::TooLarge);
    }
    let mut borrowed: Vec<_> = inputs.iter().map(|i| &i.molecule).collect();
    borrowed.push(d);
    collagene_model::bounds::model_preflight(&borrowed)?;
    let mut roles = BTreeSet::new();
    let mut bases = 0usize;
    for input in inputs {
        if input.role.trim().is_empty()
            || input.role.len() > 256
            || !roles.insert(input.role.as_str())
            || input.molecule.sequence.alphabet() != Alphabet::Dna
        {
            return Err(bad("unique bounded DNA source roles required"));
        }
        bases = bases
            .checked_add(input.molecule.sequence.len())
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
    }
    let mut prepared = 0usize;
    for strand in &p.substrates.strands {
        if strand.sequence.alphabet() != Alphabet::Dna || strand.sequence.is_empty() {
            return Err(bad("prepared strand must be nonempty DNA"));
        }
        prepared = prepared
            .checked_add(strand.sequence.len())
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
    }
    if d.sequence.alphabet() != Alphabet::Dna
        || d.sequence != v.top
        || d.sequence.is_empty()
        || d.topology
            != match v.topology {
                ContourTopology::Linear => Topology::Linear,
                ContourTopology::Circular => Topology::Circular,
            }
        || v.bottom.len() != v.top.len()
    {
        return Err(bad("destination differs from paired-contour top reference"));
    }
    if !v
        .top
        .as_str()
        .bytes()
        .zip(v.bottom.as_str().bytes().rev())
        .all(|(a, b)| a == collagene_model::complement(b, Alphabet::Dna))
    {
        return Err(bad(
            "view bottom is not the top-reference reverse complement",
        ));
    }
    let count = p
        .maps
        .len()
        .checked_mul(v.source_maps.len())
        .and_then(|n| n.checked_add(prepared))
        .and_then(|n| n.checked_add(v.top.len()))
        .ok_or(Error::TooLarge)?;
    if count > MAX_WORK {
        return Err(Error::TooLarge);
    }
    Ok(())
}

/// Emit up to two mappings for a subinterval of one source turn. Source bounds
/// are forward coordinates; reverse changes destination reading direction only.
fn compose(
    source: &Molecule,
    start: usize,
    length: usize,
    reverse: bool,
    offset: usize,
    count: usize,
    destination: usize,
    out: &mut Vec<MappingSegment>,
) -> Result<(), Error> {
    let n = source.sequence.len();
    if start >= n
        || length == 0
        || length > n
        || offset.checked_add(count).is_none_or(|e| e > length)
        || (source.topology == Topology::Linear && start.checked_add(length).is_none_or(|e| e > n))
    {
        return Err(bad("invalid source correspondence interval"));
    }
    let first = if reverse {
        (start + length - 1 - offset) % n
    } else {
        (start + offset) % n
    };
    let take = count.min(if reverse { first + 1 } else { n - first });
    if out.len() + 2 > MAX_MAPS {
        return Err(Error::TooLarge);
    }
    if take > 0 {
        out.push(MappingSegment {
            source_start: if reverse { first + 1 - take } else { first },
            source_end: if reverse { first + 1 } else { first + take },
            destination_start: destination,
            reverse,
        });
    }
    let remaining = count - take;
    if remaining > 0 {
        out.push(MappingSegment {
            source_start: if reverse { n - remaining } else { 0 },
            source_end: if reverse { n } else { remaining },
            destination_start: destination + take,
            reverse,
        });
    }
    Ok(())
}

/// Recompute annotation correspondence from supplied fresh preparation/view.
/// Serialized evidence is historical data only and cannot authorize projections.
/// No destination/source data is mutated, and bottom pairing contributes no
/// duplicate reference ancestry. Repeated source-copy roles are withheld intact.
pub fn project(
    inputs: &[topo::Input],
    preparation: &PreparationEvidence,
    view: &DuplexView,
    destination: &Molecule,
) -> Result<AnnotationEvidence, Error> {
    preflight(inputs, preparation, view, destination)?;
    let roles: BTreeMap<_, _> = inputs
        .iter()
        .enumerate()
        .map(|(i, s)| (s.role.as_str(), i))
        .collect();
    let mut per_role: Vec<Vec<MappingSegment>> = vec![Vec::new(); inputs.len()];
    let mut pcr = vec![false; inputs.len()];
    let mut work = 0;
    // Validate complete nonoverlapping prepared-map coverage before composing.
    let mut prepared_ranges: Vec<Vec<(usize, usize)>> =
        vec![Vec::new(); preparation.substrates.strands.len()];
    for map in &preparation.maps {
        let strand = preparation
            .substrates
            .strands
            .get(map.strand)
            .ok_or_else(|| bad("prepared strand index"))?;
        if map.length == 0 {
            return Err(bad("empty prepared map"));
        }
        let e = end(map.strand_start, map.length, strand.sequence.len())?;
        prepared_ranges[map.strand].push((map.strand_start, e));
        match &map.origin {
            Origin::ReferenceCorrespondence(i) => {
                let source = inputs
                    .get(
                        *roles
                            .get(i.source_role.as_str())
                            .ok_or_else(|| bad("unknown origin source role"))?,
                    )
                    .unwrap();
                if map.length != i.length {
                    return Err(bad("prepared/reference interval lengths differ"));
                }
                let mut checks = Vec::new();
                compose(
                    &source.molecule,
                    i.start,
                    i.length,
                    i.reverse,
                    0,
                    i.length,
                    0,
                    &mut checks,
                )?;
                for check in checks {
                    for k in 0..check.source_end - check.source_start {
                        let call = source.molecule.sequence.as_str().as_bytes()[if check.reverse {
                            check.source_end - 1 - k
                        } else {
                            check.source_start + k
                        }];
                        let call = if check.reverse {
                            collagene_model::complement(call, Alphabet::Dna)
                        } else {
                            call
                        };
                        if strand.sequence.as_str().as_bytes()
                            [map.strand_start + check.destination_start + k]
                            != call
                        {
                            return Err(bad(
                                "prepared DNA differs from declared source correspondence",
                            ));
                        }
                    }
                }
                charge(&mut work, map.length)?;
            }
            Origin::ExactPcrCopy {
                template_role,
                template_start,
                template_length,
                forward_tail_length,
                reverse_tail_length,
                reverse,
                ..
            } => {
                let source = &inputs[*roles
                    .get(template_role.as_str())
                    .ok_or_else(|| bad("unknown PCR template role"))?]
                .molecule;
                let length = forward_tail_length
                    .checked_add(*template_length)
                    .and_then(|n| n.checked_add(*reverse_tail_length))
                    .ok_or(Error::TooLarge)?;
                if map.length != length {
                    return Err(bad("PCR copy length/tails differ"));
                }
                let tail = if *reverse {
                    *reverse_tail_length
                } else {
                    *forward_tail_length
                };
                let mut checks = Vec::new();
                compose(
                    source,
                    *template_start,
                    *template_length,
                    *reverse,
                    0,
                    *template_length,
                    tail,
                    &mut checks,
                )?;
                for check in checks {
                    for k in 0..check.source_end - check.source_start {
                        let call = source.sequence.as_str().as_bytes()[if check.reverse {
                            check.source_end - 1 - k
                        } else {
                            check.source_start + k
                        }];
                        let call = if check.reverse {
                            collagene_model::complement(call, Alphabet::Dna)
                        } else {
                            call
                        };
                        if strand.sequence.as_str().as_bytes()
                            [map.strand_start + check.destination_start + k]
                            != call
                        {
                            return Err(bad(
                                "PCR template correspondence differs from prepared DNA",
                            ));
                        }
                    }
                }
                charge(&mut work, *template_length)?;
            }
            Origin::SynthesizedTail { base, .. } => {
                if map.length != 1
                    || &strand.sequence.as_str()[map.strand_start..e] != base.as_str()
                {
                    return Err(bad("synthesized tail differs from prepared base"));
                }
            }
        }
    }
    for (i, ranges) in prepared_ranges.iter_mut().enumerate() {
        ranges.sort_unstable();
        let mut cursor = 0;
        for &(a, b) in ranges.iter() {
            if a != cursor {
                return Err(bad("prepared map coverage overlaps or has gaps"));
            }
            cursor = b;
        }
        if cursor != preparation.substrates.strands[i].sequence.len() {
            return Err(bad("incomplete prepared strand provenance"));
        }
    }
    let mut top_calls: Vec<Vec<(usize, usize)>> =
        vec![Vec::new(); preparation.substrates.strands.len()];
    let mut cursor = 0;
    for run in &view.source_maps {
        if run.length == 0 || run.reference_start != cursor {
            return Err(bad("view maps must cover the reference exactly in order"));
        }
        cursor = end(run.reference_start, run.length, destination.sequence.len())?;
        let top = preparation
            .substrates
            .strands
            .get(run.top_strand)
            .ok_or_else(|| bad("top strand index"))?;
        let bottom = preparation
            .substrates
            .strands
            .get(run.bottom_strand)
            .ok_or_else(|| bad("bottom strand index"))?;
        let te = end(run.top_start, run.length, top.sequence.len())?;
        let be = end(run.bottom_start, run.length, bottom.sequence.len())?;
        top_calls[run.top_strand].push((run.top_start, te));
        if top.sequence.as_str()[run.top_start..te]
            != destination.sequence.as_str()[run.reference_start..cursor]
        {
            return Err(bad("actual top calls disagree with destination"));
        }
        for k in 0..run.length {
            if destination.sequence.as_str().as_bytes()[run.reference_start + k]
                != collagene_model::complement(
                    bottom.sequence.as_str().as_bytes()[be - 1 - k],
                    Alphabet::Dna,
                )
            {
                return Err(bad("paired bottom calls disagree with destination"));
            }
        }
        charge(&mut work, run.length)?;
        for map in &preparation.maps {
            charge(&mut work, 1)?;
            if map.strand != run.top_strand {
                continue;
            }
            let (role, start, length, reverse, offset) = match &map.origin {
                Origin::ReferenceCorrespondence(i) => (
                    i.source_role.as_str(),
                    i.start,
                    i.length,
                    i.reverse,
                    map.strand_start,
                ),
                Origin::ExactPcrCopy {
                    template_role,
                    template_start,
                    template_length,
                    forward_tail_length,
                    reverse_tail_length,
                    reverse,
                    ..
                } => (
                    template_role.as_str(),
                    *template_start,
                    *template_length,
                    *reverse,
                    map.strand_start
                        + if *reverse {
                            *reverse_tail_length
                        } else {
                            *forward_tail_length
                        },
                ),
                Origin::SynthesizedTail { .. } => continue,
            };
            let lo = run.top_start.max(offset);
            let hi = te.min(offset + length);
            if lo >= hi {
                continue;
            }
            let role = *roles
                .get(role)
                .ok_or_else(|| bad("unknown composed source role"))?;
            compose(
                &inputs[role].molecule,
                start,
                length,
                reverse,
                lo - offset,
                hi - lo,
                run.reference_start + lo - run.top_start,
                &mut per_role[role],
            )?;
            pcr[role] |= matches!(&map.origin, Origin::ExactPcrCopy { .. });
        }
    }
    if cursor != destination.sequence.len() {
        return Err(bad("view maps omit destination calls"));
    }
    for ranges in &mut top_calls {
        ranges.sort_unstable();
        if ranges.windows(2).any(|p| p[0].1 > p[1].0) {
            return Err(bad(
                "one prepared top call maps to multiple destination positions",
            ));
        }
    }
    // Before ANY source/feature clone, reserve conservative copy/expansion work.
    let mut feature_bytes = 0usize;
    let mut expansion = 0usize;
    let mut total_maps = 0usize;
    for (input, maps) in inputs.iter().zip(&per_role) {
        feature_bytes = feature_bytes
            .checked_add(size(&input.molecule.features, MAX_FEATURE_BYTES)?)
            .filter(|n| *n <= MAX_FEATURE_BYTES)
            .ok_or(Error::TooLarge)?;
        let count = input
            .molecule
            .features
            .iter()
            .map(|f| nodes(&f.location))
            .sum::<usize>();
        expansion = expansion
            .checked_add(count.checked_mul(maps.len() + 1).ok_or(Error::TooLarge)?)
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
        total_maps = total_maps
            .checked_add(maps.len())
            .filter(|n| *n <= MAX_MAPS)
            .ok_or(Error::TooLarge)?;
        charge(
            &mut work,
            count.checked_mul(maps.len() + 1).ok_or(Error::TooLarge)?,
        )?;
    }
    let copies = feature_bytes
        .checked_mul(3)
        .and_then(|n| n.checked_add(expansion * 128))
        .and_then(|n| n.checked_add(total_maps * 256))
        .ok_or(Error::TooLarge)?;
    if copies > MAX_COPY_BYTES {
        return Err(Error::TooLarge);
    }
    let mut evidence=AnnotationEvidence{sources:Vec::new(),features:Vec::new(),warnings:vec!["Annotation correspondence follows selected top-reference source calls only. It establishes neither inherited atoms/chemical state nor biological feature validity; bottom pairing is not duplicate annotation ancestry.".into()]};
    for (role, input) in inputs.iter().enumerate() {
        let mut maps = std::mem::take(&mut per_role[role]);
        maps.sort_by_key(|m| (m.source_start, m.source_end, m.destination_start, m.reverse));
        let repeated = maps.windows(2).any(|p| p[0].source_end > p[1].source_start);
        let mut result = RoleEvidence {
            role: input.role.clone(),
            source_revision: input.molecule.revision(),
            maps,
            projection_supported: !repeated,
            projected: 0,
            unmapped: Vec::new(),
            dropped: 0,
            warnings: Vec::new(),
        };
        if repeated {
            warning(
                &mut result.warnings,
                "Repeated source-to-multiple-destination correspondence is unsupported; all annotations for this role remain read-only source evidence.",
            );
        }
        if result.maps.is_empty() {
            warning(
                &mut result.warnings,
                "No selected top-reference correspondence exists for this role; pairing on the bottom strand is not duplicate annotation ancestry.",
            );
        }
        if pcr[role] {
            warning(
                &mut result.warnings,
                "PCR template correspondence only: primer tails and synthesized tails have no template annotation ancestry. Primer annealing lengths are absent from this provenance DTO, so every projected annotation from this role requires review.",
            );
        }
        let mut shadow = Molecule {
            name: "TOPO annotation source".into(),
            sequence: input.molecule.sequence.clone(),
            topology: input.molecule.topology,
            features: Vec::with_capacity(1),
            metadata: Default::default(),
            original_genbank: None,
        };
        for (source_index, original) in input.molecule.features.iter().enumerate() {
            let withheld = if repeated {
                Some("Unsupported repeated source-copy correspondence")
            } else if result.maps.is_empty() {
                Some(
                    "No selected top-reference correspondence; bottom-only pairing does not confer reference ancestry",
                )
            } else if unresolved(&original.location) {
                Some("Unresolved original location cannot be asserted as a destination annotation")
            } else {
                None
            };
            if let Some(reason) = withheld {
                let mut feature = original.clone();
                feature.needs_review = true;
                result.unmapped.push(UnmappedAnnotation {
                    source_index,
                    feature,
                    reason: reason.into(),
                });
                if result.maps.is_empty() {
                    result.dropped += 1;
                }
                continue;
            }
            shadow.features.clear();
            shadow.features.push(original.clone());
            match annotations::project_features_with_topology(
                &shadow,
                &result.maps,
                destination.sequence.len(),
                destination.topology,
            ) {
                Ok(projected) => {
                    for message in projected.warnings {
                        warning(&mut result.warnings, &message);
                    }
                    if projected.features.is_empty() {
                        let mut feature = shadow.features.pop().unwrap();
                        feature.needs_review = true;
                        result.unmapped.push(UnmappedAnnotation{source_index,feature,reason:if projected.dropped_features>0{"Original feature is outside selected top-reference source intervals"}else{"Coordinate projection withheld this original feature"}.into()});
                        result.dropped += projected.dropped_features;
                    } else {
                        for mut feature in projected.features {
                            if pcr[role] {
                                feature.needs_review = true;
                                feature.qualifiers.retain(|(key, _)| {
                                    !key.eq_ignore_ascii_case("translation")
                                        && !key.eq_ignore_ascii_case("codon_start")
                                });
                                feature.qualifiers.push(("note".into(),Some("TOPO PCR-derived template correspondence; primer/tail changes and biological validity require review. Recalculate translation/frame assertions.".into())));
                            }
                            result.projected += 1;
                            evidence.features.push(feature);
                        }
                    }
                }
                Err(error) => {
                    let mut feature = shadow.features.pop().unwrap();
                    feature.needs_review = true;
                    result.unmapped.push(UnmappedAnnotation {
                        source_index,
                        feature,
                        reason: error.to_string().chars().take(1024).collect(),
                    });
                    warning(
                        &mut result.warnings,
                        "An unprojectable annotation remains read-only source evidence; supported annotations are retained separately.",
                    );
                }
            }
        }
        evidence.sources.push(result);
    }
    size(&evidence, MAX_OUTPUT_BYTES)?;
    Ok(evidence)
}
