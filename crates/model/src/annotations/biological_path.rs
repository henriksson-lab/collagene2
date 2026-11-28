//! Exact biological call order over the canonical annotation atom traversal.
//! This is coordinate evidence, not verification that an annotation is biological truth.
use super::{atoms, uncertain, unsupported};
use crate::{Error, Location, MAX_SEQUENCE_LEN, annotations::MappingSegment};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BiologicalCall {
    pub position: usize,
    pub reverse: bool,
}
impl BiologicalCall {
    /// Adjacent reference call in this biological orientation; one origin wrap
    /// is allowed only for an explicitly circular source.
    pub fn next(self, source_length: usize, topology: crate::Topology) -> Option<Self> {
        if self.position >= source_length {
            return None;
        }
        let position = if self.reverse {
            self.position
                .checked_sub(1)
                .or_else(|| (topology == crate::Topology::Circular).then_some(source_length - 1))?
        } else if self.position + 1 < source_length {
            self.position + 1
        } else if topology == crate::Topology::Circular {
            0
        } else {
            return None;
        };
        Some(Self {
            position,
            reverse: self.reverse,
        })
    }
}
#[derive(Clone, Copy, Debug)]
pub struct PathLimits {
    pub max_calls: usize,
    pub max_spans: usize,
}
impl Default for PathLimits {
    fn default() -> Self {
        Self {
            max_calls: 1_000_000,
            max_spans: 4096,
        }
    }
}
#[derive(Clone, Debug)]
struct Span {
    start: usize,
    end: usize,
    reverse: bool,
    prefix: usize,
}
/// Both indexes retain biological ordinal rather than sorting the biological path.
/// Repeated identical source calls are refused instead of choosing one occurrence.
#[derive(Clone, Debug)]
pub struct BiologicalPath {
    spans: Vec<Span>,
    inverse: Vec<usize>,
    length: usize,
}
impl BiologicalPath {
    pub fn from_location(
        location: &Location,
        source_length: usize,
        limits: PathLimits,
    ) -> Result<Self, Error> {
        if source_length == 0
            || source_length > MAX_SEQUENCE_LEN
            || !(1..=MAX_SEQUENCE_LEN).contains(&limits.max_calls)
            || !(1..=4096).contains(&limits.max_spans)
        {
            return Err(Error::TooLarge);
        }
        crate::bounds::ModelBudget::default().location(location)?;
        location.validate(source_length)?;
        if uncertain(location) {
            return Err(unsupported(
                "an exact biological path requires resolved, nonfuzzy, ordered coordinates",
            ));
        }
        let map = [MappingSegment {
            source_start: 0,
            source_end: source_length,
            destination_start: 0,
            reverse: false,
        }];
        let mut work = 0;
        let native = atoms(location, &map, &mut work)?;
        if native.is_empty() || native.len() > limits.max_spans {
            return Err(Error::TooLarge);
        }
        let mut length = 0usize;
        let mut spans = Vec::with_capacity(native.len());
        for atom in native {
            let prefix = length;
            length = length
                .checked_add(atom.source_end - atom.source_start)
                .filter(|n| *n <= limits.max_calls)
                .ok_or(Error::TooLarge)?;
            spans.push(Span {
                start: atom.source_start,
                end: atom.source_end,
                reverse: atom.source_reverse,
                prefix,
            });
        }
        let mut inverse = (0..spans.len()).collect::<Vec<_>>();
        inverse.sort_unstable_by_key(|index| {
            let span = &spans[*index];
            (span.reverse, span.start, span.end)
        });
        for pair in inverse.windows(2) {
            let a = &spans[pair[0]];
            let b = &spans[pair[1]];
            if a.reverse == b.reverse && a.end > b.start {
                return Err(unsupported(
                    "a source call occurs repeatedly in the biological path",
                ));
            }
        }
        Ok(Self {
            spans,
            inverse,
            length,
        })
    }
    pub fn len(&self) -> usize {
        self.length
    }
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }
    pub fn call_at(&self, ordinal: usize) -> Option<BiologicalCall> {
        if ordinal >= self.length {
            return None;
        }
        let index = self
            .spans
            .partition_point(|span| span.prefix <= ordinal)
            .checked_sub(1)?;
        let span = &self.spans[index];
        let offset = ordinal - span.prefix;
        Some(BiologicalCall {
            position: if span.reverse {
                span.end - 1 - offset
            } else {
                span.start + offset
            },
            reverse: span.reverse,
        })
    }
    pub fn ordinal_of(&self, call: BiologicalCall) -> Option<usize> {
        let key = (call.reverse, call.position);
        let index = self
            .inverse
            .partition_point(|index| {
                let span = &self.spans[*index];
                (span.reverse, span.start) <= key
            })
            .checked_sub(1)?;
        let span = &self.spans[self.inverse[index]];
        if span.reverse != call.reverse || call.position >= span.end {
            return None;
        }
        Some(
            span.prefix
                + if span.reverse {
                    span.end - 1 - call.position
                } else {
                    call.position - span.start
                },
        )
    }
    pub fn calls(&self) -> impl ExactSizeIterator<Item = BiologicalCall> + '_ {
        (0..self.length).map(|ordinal| self.call_at(ordinal).expect("admitted biological ordinal"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complemented_join_keeps_biological_order_and_exact_inverse_index() {
        let location = Location::Complement(Box::new(Location::Join(vec![
            Location::range(1, 4),
            Location::range(7, 9),
        ])));
        let path = BiologicalPath::from_location(&location, 10, PathLimits::default()).unwrap();
        let calls = path.calls().collect::<Vec<_>>();
        assert_eq!(
            calls.iter().map(|call| call.position).collect::<Vec<_>>(),
            vec![8, 7, 3, 2, 1]
        );
        for (ordinal, call) in calls.into_iter().enumerate() {
            assert!(call.reverse);
            assert_eq!(path.ordinal_of(call), Some(ordinal));
            assert_eq!(
                path.ordinal_of(BiologicalCall {
                    reverse: false,
                    ..call
                }),
                None
            );
        }
        assert_eq!(path.call_at(5), None);
    }
    #[test]
    fn repeated_call_uncertain_order_and_limits_do_not_supply_a_path() {
        for location in [
            Location::Join(vec![Location::range(0, 3), Location::range(2, 4)]),
            Location::Order(vec![Location::range(0, 3)]),
            Location::Range {
                start: 0,
                end: 3,
                before: true,
                after: false,
            },
        ] {
            assert!(BiologicalPath::from_location(&location, 4, PathLimits::default()).is_err());
        }
        assert!(
            BiologicalPath::from_location(
                &Location::range(0, 4),
                4,
                PathLimits {
                    max_calls: 3,
                    max_spans: 1
                }
            )
            .is_err()
        );
    }
    #[test]
    fn projection_dispositions_keep_original_ordinals_for_all_three_outcomes() {
        use crate::{Feature, Molecule, Topology};
        let mut source = Molecule::dna("literal", "AACCGGTT", Topology::Linear).unwrap();
        for location in [
            Location::range(2, 4),
            Location::range(6, 8),
            Location::Unresolved("remote:1..2".into()),
        ] {
            source.features.push(Feature {
                kind: "misc_feature".into(),
                location,
                qualifiers: vec![],
                needs_review: true,
            });
        }
        let map = [MappingSegment {
            source_start: 2,
            source_end: 6,
            destination_start: 0,
            reverse: false,
        }];
        let projection =
            super::super::project_features_with_topology(&source, &map, 4, Topology::Linear)
                .unwrap();
        assert_eq!(projection.source_feature_indices, vec![0]);
        assert_eq!(projection.unmapped_source_feature_indices, vec![2]);
        assert_eq!(projection.dropped_source_feature_indices, vec![1]);
        assert_eq!(projection.features.len(), 1);
        assert_eq!(projection.unmapped_features.len(), 1);
        assert_eq!(projection.dropped_features, 1);
    }
}
