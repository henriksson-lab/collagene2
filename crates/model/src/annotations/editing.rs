//! Native edits retain endpoint-envelope semantics. Inserting inside a range
//! expands that range; source-subset projection would instead exclude new bases.
use super::{MappingSegment, budget, qualifiers, uncertain};
use crate::{Error, Feature, Location, Molecule};

/// Non-deserializable receipt emitted by the native edit feature loop.
/// Ordinals address the original source and align with retained product features.
#[derive(Debug)]
pub struct ReplacementProjection {
    pub map: Vec<MappingSegment>,
    pub source_feature_indices: Vec<usize>,
    pub dropped_source_feature_indices: Vec<usize>,
    pub start: usize,
    pub end: usize,
    pub inserted_length: usize,
    pub unchanged: bool,
}
pub(crate) struct Edit {
    start: usize,
    end: usize,
    inserted: usize,
}
impl Edit {
    pub(crate) fn new(start: usize, end: usize, inserted: usize) -> Self {
        Self {
            start,
            end,
            inserted,
        }
    }
    fn map_location(
        &self,
        location: &Location,
        work: &mut usize,
    ) -> Result<(Option<Location>, bool), Error> {
        budget(work)?;
        let shift = |position: usize| position - (self.end - self.start) + self.inserted;
        match location {
            Location::Range {
                start,
                end,
                before,
                after,
            } => {
                if *end <= self.start {
                    return Ok((Some(location.clone()), false));
                }
                if *start >= self.end {
                    return Ok((
                        Some(Location::Range {
                            start: shift(*start),
                            end: shift(*end),
                            before: *before,
                            after: *after,
                        }),
                        false,
                    ));
                }
                let start = if *start < self.start {
                    *start
                } else {
                    self.start
                };
                let end = if *end > self.end {
                    shift(*end)
                } else {
                    self.start + self.inserted
                };
                Ok((
                    (start < end).then_some(Location::Range {
                        start,
                        end,
                        before: *before,
                        after: *after,
                    }),
                    true,
                ))
            }
            Location::Complement(inner) => {
                let (mapped, touched) = self.map_location(inner, work)?;
                Ok((
                    mapped.map(|location| Location::Complement(Box::new(location))),
                    touched,
                ))
            }
            Location::Join(parts) | Location::Order(parts) => {
                let mut mapped = Vec::with_capacity(parts.len());
                let mut touched = false;
                for part in parts {
                    let (location, changed) = self.map_location(part, work)?;
                    touched |= changed;
                    if let Some(location) = location {
                        mapped.push(location);
                    }
                }
                Ok((
                    if mapped.is_empty() {
                        None
                    } else if matches!(location, Location::Join(_)) {
                        Some(Location::Join(mapped))
                    } else {
                        Some(Location::Order(mapped))
                    },
                    touched,
                ))
            }
            Location::Unresolved(_) => Err(Error::Unsupported(
                "editing documents with unresolved feature locations".into(),
            )),
        }
    }
    fn retained_segments(&self, length: usize) -> Vec<MappingSegment> {
        let mut segments = Vec::with_capacity(2);
        if self.start > 0 {
            segments.push(MappingSegment {
                source_start: 0,
                source_end: self.start,
                destination_start: 0,
                reverse: false,
            });
        }
        if self.end < length {
            segments.push(MappingSegment {
                source_start: self.end,
                source_end: length,
                destination_start: self.start + self.inserted,
                reverse: false,
            });
        }
        segments
    }
    pub(crate) fn features(
        &self,
        source: &Molecule,
        destination_length: usize,
        capture: bool,
    ) -> Result<(Vec<Feature>, Option<ReplacementProjection>), Error> {
        let length = source.sequence.len();
        if self.start > self.end
            || self.end > length
            || length - (self.end - self.start) + self.inserted != destination_length
        {
            return Err(Error::Coordinates(
                "invalid native annotation edit geometry".into(),
            ));
        }
        let mut work = 0;
        qualifiers::preflight(source, &mut work)?;
        let segments = self.retained_segments(length);
        let mut result = Vec::with_capacity(source.features.len());
        let mut ordinals = Vec::new();
        let mut dropped = Vec::new();
        for (index, original) in source.features.iter().enumerate() {
            let (location, touched) = self.map_location(&original.location, &mut work)?;
            let Some(location) = location else {
                if capture {
                    dropped.push(index);
                }
                continue;
            };
            let review = touched
                || uncertain(&original.location)
                || !original.location.is_single_range()
                    && super::disrupted(
                        &original.location,
                        &segments,
                        source,
                        destination_length,
                        source.topology,
                        &mut work,
                    )?;
            let qualifiers = qualifiers::project(
                original,
                source,
                &segments,
                destination_length,
                false,
                review,
                &mut work,
            )?;
            let feature = Feature {
                kind: original.kind.clone(),
                location,
                qualifiers,
                needs_review: original.needs_review || review,
            };
            feature.location.validate(destination_length)?;
            result.push(feature);
            if capture {
                ordinals.push(index);
            }
        }
        let projection = capture.then(|| ReplacementProjection {
            map: segments,
            source_feature_indices: ordinals,
            dropped_source_feature_indices: dropped,
            start: self.start,
            end: self.end,
            inserted_length: self.inserted,
            unchanged: false,
        });
        Ok((result, projection))
    }
}
