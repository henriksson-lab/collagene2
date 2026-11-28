//! Independently implemented coordinate projections. A segment maps source reference bases
//! into a derived document; biological feature order is retained, not sorted by display position.
use crate::{Alphabet, Error, Feature, Location, MAX_SEQUENCE_LEN, Molecule, Topology};
use serde::{Deserialize, Serialize};
pub(crate) mod editing;
pub mod insdc;
mod qualifiers;

const MAX_SEGMENTS: usize = 4096;
const MAX_WORK: usize = 1_000_000;
const MAX_WARNINGS: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct MappingSegment {
    pub source_start: usize,
    pub source_end: usize,
    pub destination_start: usize,
    /// Source bases are reverse-complemented in this destination interval.
    pub reverse: bool,
}
impl MappingSegment {
    fn len(&self) -> usize {
        self.source_end - self.source_start
    }
    fn destination_end(&self) -> usize {
        self.destination_start + self.len()
    }
}
#[derive(Debug, Clone)]
pub struct Projection {
    pub features: Vec<Feature>,
    pub warnings: Vec<String>,
    /// Original, unprojectable source features, not annotations asserted to exist in the product.
    /// Callers constructing a product must retain these as evidence or retain their source graph.
    pub unmapped_features: Vec<Feature>,
    pub dropped_features: usize,
}
fn bad(message: impl Into<String>) -> Error {
    Error::Coordinates(message.into())
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}
fn warn(warnings: &mut Vec<String>, value: impl Into<String>) {
    if warnings.len() < MAX_WARNINGS {
        warnings.push(value.into());
    } else if warnings.len() == MAX_WARNINGS {
        warnings.push("Further annotation projection warnings omitted; inspect the source annotations and review flags.".into());
    }
}
fn budget(work: &mut usize) -> Result<(), Error> {
    *work = work.saturating_add(1);
    if *work > MAX_WORK {
        Err(unsupported(
            "annotation projection exceeds its one-million-node/intersection budget",
        ))
    } else {
        Ok(())
    }
}
fn mappings(
    source: &Molecule,
    segments: &[MappingSegment],
    destination_length: usize,
    allow_destination_overlap: bool,
) -> Result<Vec<MappingSegment>, Error> {
    source.validate()?;
    if destination_length == 0
        || destination_length > MAX_SEQUENCE_LEN
        || segments.len() > MAX_SEGMENTS
    {
        return Err(Error::TooLarge);
    }
    for segment in segments {
        if segment.source_start >= segment.source_end || segment.source_end > source.sequence.len()
        {
            return Err(bad("source mapping interval is empty or out of bounds"));
        }
        if segment
            .destination_start
            .checked_add(segment.source_end - segment.source_start)
            .is_none_or(|end| end > destination_length)
        {
            return Err(bad("destination mapping interval is out of bounds"));
        }
        if segment.reverse && source.sequence.alphabet() == Alphabet::Protein {
            return Err(unsupported(
                "protein annotation projection has no reverse-complement operation",
            ));
        }
    }
    let mut source_order = segments.to_vec();
    source_order.sort_by_key(|s| (s.source_start, s.source_end));
    if source_order
        .windows(2)
        .any(|p| p[0].source_end > p[1].source_start)
    {
        return Err(bad(
            "overlapping source mappings require separate per-copy projection calls",
        ));
    }
    let mut destination_order = segments.to_vec();
    destination_order.sort_by_key(|s| s.destination_start);
    if !allow_destination_overlap
        && destination_order
            .windows(2)
            .any(|p| p[0].destination_end() > p[1].destination_start)
    {
        return Err(bad("destination mapping intervals overlap"));
    }
    Ok(source_order)
}
fn identity(segments: &[MappingSegment], source_length: usize, destination_length: usize) -> bool {
    source_length == destination_length
        && segments
            .iter()
            .all(|s| !s.reverse && s.destination_start == s.source_start)
        && segments.iter().map(MappingSegment::len).sum::<usize>() == source_length
}
fn has_unresolved(location: &Location) -> bool {
    match location {
        Location::Unresolved(_) => true,
        Location::Complement(l) => has_unresolved(l),
        Location::Join(parts) | Location::Order(parts) => parts.iter().any(has_unresolved),
        _ => false,
    }
}
fn uncertain(location: &Location) -> bool {
    match location {
        Location::Range { before, after, .. } => *before || *after,
        Location::Complement(l) => uncertain(l),
        Location::Join(parts) => parts.iter().any(uncertain),
        Location::Order(_) | Location::Unresolved(_) => true,
    }
}
fn complement(location: Location) -> Location {
    match location {
        Location::Complement(inner) => *inner,
        other => Location::Complement(Box::new(other)),
    }
}
fn compound(parts: Vec<Location>, order: bool) -> Option<Location> {
    if parts.is_empty() {
        None
    } else if parts.len() == 1 && !order {
        parts.into_iter().next()
    } else if order {
        Some(Location::Order(parts))
    } else {
        Some(Location::Join(parts))
    }
}
#[derive(Clone, Copy)]
struct Atom {
    source_start: usize,
    source_end: usize,
    destination_start: usize,
    destination_end: usize,
    source_reverse: bool,
    destination_reverse: bool,
}
fn intersections(
    start: usize,
    end: usize,
    segments: &[MappingSegment],
) -> impl Iterator<Item = (usize, usize, &MappingSegment)> {
    let first = segments.partition_point(|s| s.source_end <= start);
    segments[first..]
        .iter()
        .take_while(move |s| s.source_start < end)
        .map(move |s| (start.max(s.source_start), end.min(s.source_end), s))
}
fn project_location(
    location: &Location,
    segments: &[MappingSegment],
    work: &mut usize,
) -> Result<(Option<Location>, bool), Error> {
    budget(work)?;
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => {
            let pieces: Vec<_> = intersections(*start, *end, segments).collect();
            let covered = pieces.iter().map(|(a, b, _)| b - a).sum::<usize>();
            let mut out = Vec::with_capacity(pieces.len());
            for (i, &(a, b, segment)) in pieces.iter().enumerate() {
                budget(work)?;
                let lo_partial = if a == *start {
                    *before
                } else {
                    i == 0 || pieces[i - 1].1 != a
                };
                let hi_partial = if b == *end {
                    *after
                } else {
                    i + 1 == pieces.len() || pieces[i + 1].0 != b
                };
                let location = if segment.reverse {
                    complement(Location::Range {
                        start: segment.destination_start + segment.source_end - b,
                        end: segment.destination_start + segment.source_end - a,
                        before: hi_partial,
                        after: lo_partial,
                    })
                } else {
                    Location::Range {
                        start: segment.destination_start + a - segment.source_start,
                        end: segment.destination_start + b - segment.source_start,
                        before: lo_partial,
                        after: hi_partial,
                    }
                };
                out.push(location);
            }
            Ok((compound(out, false), covered < end - start))
        }
        Location::Complement(inner) => {
            let (mapped, partial) = project_location(inner, segments, work)?;
            Ok((mapped.map(complement), partial))
        }
        Location::Join(parts) | Location::Order(parts) => {
            let mut out = Vec::new();
            let mut partial = false;
            for part in parts {
                let (mapped, cut) = project_location(part, segments, work)?;
                partial |= cut;
                if let Some(mapped) = mapped {
                    out.push(mapped);
                }
            }
            Ok((
                compound(out, matches!(location, Location::Order(_))),
                partial,
            ))
        }
        Location::Unresolved(_) => Err(unsupported(
            "unresolved source location cannot be projected",
        )),
    }
}
fn atoms(
    location: &Location,
    segments: &[MappingSegment],
    work: &mut usize,
) -> Result<Vec<Atom>, Error> {
    budget(work)?;
    match location {
        Location::Range { start, end, .. } => {
            let mut out = Vec::new();
            for (a, b, segment) in intersections(*start, *end, segments) {
                budget(work)?;
                let (destination_start, destination_end) = if segment.reverse {
                    (
                        segment.destination_start + segment.source_end - b,
                        segment.destination_start + segment.source_end - a,
                    )
                } else {
                    (
                        segment.destination_start + a - segment.source_start,
                        segment.destination_start + b - segment.source_start,
                    )
                };
                out.push(Atom {
                    source_start: a,
                    source_end: b,
                    destination_start,
                    destination_end,
                    source_reverse: false,
                    destination_reverse: segment.reverse,
                });
            }
            Ok(out)
        }
        Location::Complement(inner) => {
            let mut out = atoms(inner, segments, work)?;
            out.reverse();
            for a in &mut out {
                a.source_reverse = !a.source_reverse;
                a.destination_reverse = !a.destination_reverse;
            }
            Ok(out)
        }
        Location::Join(parts) | Location::Order(parts) => {
            let mut out = Vec::new();
            for p in parts {
                out.extend(atoms(p, segments, work)?);
            }
            Ok(out)
        }
        Location::Unresolved(_) => Ok(vec![]),
    }
}
fn adjacent(
    a_start: usize,
    a_end: usize,
    b_start: usize,
    b_end: usize,
    a_reverse: bool,
    b_reverse: bool,
    length: usize,
    topology: Topology,
) -> bool {
    if a_reverse != b_reverse {
        return false;
    }
    if a_reverse {
        a_start == b_end || topology == Topology::Circular && a_start == 0 && b_end == length
    } else {
        a_end == b_start || topology == Topology::Circular && a_end == length && b_start == 0
    }
}
fn disrupted(
    location: &Location,
    segments: &[MappingSegment],
    source: &Molecule,
    destination_length: usize,
    destination_topology: Topology,
    work: &mut usize,
) -> Result<bool, Error> {
    let parts = atoms(location, segments, work)?;
    if let Some(first) = parts.first() {
        let direction = first.source_reverse ^ first.destination_reverse;
        if parts
            .iter()
            .any(|p| (p.source_reverse ^ p.destination_reverse) != direction)
        {
            return Ok(true);
        }
        // A previously coherent spliced feature must not acquire a new exon reversal/order.
        // Circular reference traversal permits one origin crossing, never arbitrary permutations.
        if parts
            .iter()
            .all(|p| p.source_reverse == first.source_reverse)
        {
            let source_wraps = parts
                .windows(2)
                .filter(|p| {
                    if first.source_reverse {
                        p[1].source_end > p[0].source_start
                    } else {
                        p[0].source_end > p[1].source_start
                    }
                })
                .count();
            let source_limit = usize::from(source.topology == Topology::Circular);
            if source_wraps <= source_limit {
                let destination_wraps = parts
                    .windows(2)
                    .filter(|p| {
                        if first.destination_reverse {
                            p[1].destination_end > p[0].destination_start
                        } else {
                            p[0].destination_end > p[1].destination_start
                        }
                    })
                    .count();
                if destination_wraps > usize::from(destination_topology == Topology::Circular) {
                    return Ok(true);
                }
            }
        }
    }
    Ok(parts.windows(2).any(|p| {
        adjacent(
            p[0].source_start,
            p[0].source_end,
            p[1].source_start,
            p[1].source_end,
            p[0].source_reverse,
            p[1].source_reverse,
            source.sequence.len(),
            source.topology,
        ) && !adjacent(
            p[0].destination_start,
            p[0].destination_end,
            p[1].destination_start,
            p[1].destination_end,
            p[0].destination_reverse,
            p[1].destination_reverse,
            destination_length,
            destination_topology,
        )
    }))
}
fn mark_partial(location: &mut Location) {
    match location {
        Location::Range { before, after, .. } => {
            *before = true;
            *after = true;
        }
        Location::Complement(l) => mark_partial(l),
        Location::Join(parts) | Location::Order(parts) => {
            for p in parts {
                mark_partial(p);
            }
        }
        Location::Unresolved(_) => {}
    }
}
fn short_label(feature: &Feature) -> String {
    feature.label().chars().take(128).collect()
}
pub fn coordinate_qualifier(key: &str) -> bool {
    ["transl_except", "anticodon", "rpt_unit_range", "coded_by"]
        .iter()
        .any(|name| key.eq_ignore_ascii_case(name))
}

/// Deactivate biological assertions when an annotation is manually reinterpreted.
/// The source values remain bounded evidence, never active claims on the new feature.
/// All admission checks complete before replacing the caller's feature.
pub fn deactivate_assertions(feature: &mut Feature) -> Result<(), Error> {
    crate::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
    let mut work = 0;
    qualifiers::preflight_feature(feature, &mut work)?;
    let values = qualifiers::deactivated(feature, &mut work)?;
    let mut result = feature.clone();
    result.qualifiers = values;
    result.needs_review = true;
    crate::bounds::feature_preflight(&[std::slice::from_ref(&result)])?;
    *feature = result;
    Ok(())
}

/// Project into a linear destination. Use the topology-aware variant for a circular product.
pub fn project_features(
    source: &Molecule,
    segments: &[MappingSegment],
    destination_length: usize,
) -> Result<Projection, Error> {
    project_features_with_topology(source, segments, destination_length, Topology::Linear)
}
/// Maps are one-to-one within a call. Project repeated source copies separately, then combine.
/// This maps annotations only; the caller must independently construct/validate destination bases.
pub fn project_features_with_topology(
    source: &Molecule,
    segments: &[MappingSegment],
    destination_length: usize,
    destination_topology: Topology,
) -> Result<Projection, Error> {
    project_with_maps(
        source,
        segments,
        destination_length,
        destination_topology,
        false,
    )
}

/// Internal homology-collapse mapping. Source intervals must remain disjoint. Each
/// feature is projected only when its retained bases map injectively; otherwise
/// its original annotation remains reviewed evidence, never an exact product gene.
pub fn project_collapsing_features(
    source: &Molecule,
    segments: &[MappingSegment],
    destination_length: usize,
    destination_topology: Topology,
) -> Result<Projection, Error> {
    project_with_maps(
        source,
        segments,
        destination_length,
        destination_topology,
        true,
    )
}

fn project_with_maps(
    source: &Molecule,
    segments: &[MappingSegment],
    destination_length: usize,
    destination_topology: Topology,
    allow_collapse: bool,
) -> Result<Projection, Error> {
    let mut work = 0;
    qualifiers::preflight(source, &mut work)?;
    let segments = mappings(source, segments, destination_length, allow_collapse)?;
    if destination_topology == Topology::Circular && source.sequence.alphabet() == Alphabet::Protein
    {
        return Err(unsupported("circular protein products are not supported"));
    }
    let is_identity = identity(&segments, source.sequence.len(), destination_length);
    let mut projection = Projection {
        features: vec![],
        warnings: vec![],
        unmapped_features: vec![],
        dropped_features: 0,
    };
    for original in &source.features {
        budget(&mut work)?;
        if has_unresolved(&original.location) && !is_identity {
            let mut evidence = original.clone();
            evidence.needs_review = true;
            projection.unmapped_features.push(evidence);
            warn(
                &mut projection.warnings,
                format!(
                    "Source feature '{}' has an unresolved location; retained as unmapped evidence, not a destination annotation.",
                    short_label(original)
                ),
            );
            continue;
        }
        // A circular self-overlap identifies two source homology copies. Trim
        // the map to the union of this feature's intervals before checking that
        // distinct source bases do not claim the same product reference base.
        let feature_segments;
        let segments = if allow_collapse {
            let mut spans: Vec<_> = original
                .location
                .spans()
                .into_iter()
                .map(|(a, b, _)| (a, b))
                .collect();
            spans.sort_unstable();
            let mut merged: Vec<(usize, usize)> = Vec::new();
            for (a, b) in spans {
                budget(&mut work)?;
                if let Some(last) = merged.last_mut().filter(|last| a <= last.1) {
                    last.1 = last.1.max(b);
                } else {
                    merged.push((a, b));
                }
            }
            let mut trimmed = Vec::new();
            for segment in &segments {
                for &(a, b) in &merged {
                    budget(&mut work)?;
                    let lo = a.max(segment.source_start);
                    let hi = b.min(segment.source_end);
                    if lo < hi {
                        trimmed.push(MappingSegment {
                            source_start: lo,
                            source_end: hi,
                            destination_start: segment.destination_start
                                + if segment.reverse {
                                    segment.source_end - hi
                                } else {
                                    lo - segment.source_start
                                },
                            reverse: segment.reverse,
                        });
                    }
                }
            }
            let mut destinations = trimmed.clone();
            destinations.sort_by_key(|s| s.destination_start);
            if destinations
                .windows(2)
                .any(|p| p[0].destination_end() > p[1].destination_start)
            {
                let mut evidence = original.clone();
                evidence.needs_review = true;
                projection.unmapped_features.push(evidence);
                warn(
                    &mut projection.warnings,
                    format!(
                        "Source feature '{}' spans collapsed homology copies; retained as reviewed unmapped evidence because distinct source bases share product coordinates.",
                        short_label(original)
                    ),
                );
                continue;
            }
            feature_segments = trimmed;
            feature_segments.as_slice()
        } else {
            segments.as_slice()
        };
        let (mapped, partial) = if has_unresolved(&original.location) {
            (Some(original.location.clone()), false)
        } else {
            project_location(&original.location, segments, &mut work)?
        };
        let Some(mut location) = mapped else {
            projection.dropped_features += 1;
            continue;
        };
        let broken = if has_unresolved(&original.location) {
            false
        } else {
            disrupted(
                &original.location,
                segments,
                source,
                destination_length,
                destination_topology,
                &mut work,
            )?
        };
        if partial || broken {
            mark_partial(&mut location);
        }
        let review = partial || broken || uncertain(&original.location);
        let projected_qualifiers = qualifiers::project(
            original,
            source,
            segments,
            destination_length,
            is_identity && !partial && !broken,
            review,
            &mut work,
        )?;
        let mut feature = Feature {
            kind: original.kind.clone(),
            location,
            qualifiers: projected_qualifiers,
            needs_review: original.needs_review || review,
        };
        if review {
            feature.qualifiers.push(("note".into(),Some(if broken {"Collagene: source annotation is interrupted by new reference boundaries/junctions; biological continuity requires review."}else if partial {"Collagene: source annotation is only partly retained in this product; endpoints and CDS phase require review."}else{"Collagene: source location remains fuzzy, ordered or unresolved; no exact biological extraction was established."}.into())));
            warn(
                &mut projection.warnings,
                format!(
                    "Feature '{}' requires review after projection{}.",
                    short_label(original),
                    if broken {
                        " because biological continuity is interrupted"
                    } else if partial {
                        " because part of it was removed"
                    } else {
                        " because the source location is uncertain"
                    }
                ),
            );
        }
        feature
            .location
            .validate(destination_length)
            .map_err(|error| {
                unsupported(format!(
                    "projected annotation cannot be represented within location limits: {error}"
                ))
            })?;
        projection.features.push(feature);
    }
    if projection.dropped_features > 0 {
        warn(
            &mut projection.warnings,
            format!(
                "{} source feature(s) lie wholly outside the retained source intervals and were omitted.",
                projection.dropped_features
            ),
        );
    }
    Ok(projection)
}
