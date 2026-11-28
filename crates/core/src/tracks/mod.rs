//! Bounded interval index for visible reference annotation tracks.
use collagene_model::{Error, Molecule};
use std::collections::BTreeSet;
pub const MAX_INDEX_SPANS: usize = 1_000_000;
#[derive(Debug, Clone, PartialEq, Eq)]
struct Span {
    start: usize,
    end: usize,
    feature: usize,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureIndex {
    spans: Vec<Span>,
    prefix_end: Vec<usize>,
}
impl FeatureIndex {
    pub fn new(molecule: &Molecule) -> Result<Self, Error> {
        molecule.validate()?;
        let mut spans = Vec::new();
        for (feature, f) in molecule
            .features
            .iter()
            .enumerate()
            .filter(|(_, f)| f.kind != "source")
        {
            for (start, end, _) in f.location.spans() {
                if spans.len() >= MAX_INDEX_SPANS {
                    return Err(Error::TooLarge);
                }
                spans.push(Span {
                    start,
                    end,
                    feature,
                });
            }
        }
        spans.sort_by_key(|s| (s.start, s.end, s.feature));
        let mut maximum = 0;
        let prefix_end = spans
            .iter()
            .map(|s| {
                maximum = maximum.max(s.end);
                maximum
            })
            .collect();
        Ok(Self { spans, prefix_end })
    }
    /// Original feature order, with each compound feature returned once. The limit
    /// bounds displayed tracks, not an assertion that other annotations are absent.
    pub fn visible(&self, start: usize, end: usize, limit: usize) -> Vec<usize> {
        if start >= end || limit == 0 {
            return Vec::new();
        }
        let first = self.prefix_end.partition_point(|&maximum| maximum <= start);
        let mut selected = BTreeSet::new();
        for span in self.spans[first..]
            .iter()
            .take_while(|s| s.start < end)
            .filter(|s| s.end > start)
        {
            selected.insert(span.feature);
            if selected.len() > limit {
                selected.pop_last();
            }
        }
        selected.into_iter().collect()
    }
}

/// The final biological span, with its direction, only when every part has an
/// exact ordered location. Unknown-strand/order/fuzzy evidence must not acquire
/// a directional arrow merely because its reference intervals are drawable.
pub fn terminal_span(location: &collagene_model::Location) -> Option<(usize, usize, bool)> {
    fn exact(location: &collagene_model::Location) -> bool {
        match location {
            collagene_model::Location::Range {
                before: false,
                after: false,
                ..
            } => true,
            collagene_model::Location::Complement(inner) => exact(inner),
            collagene_model::Location::Join(parts) => !parts.is_empty() && parts.iter().all(exact),
            _ => false,
        }
    }
    fn last(location: &collagene_model::Location, reverse: bool) -> Option<(usize, usize, bool)> {
        match location {
            collagene_model::Location::Range { start, end, .. } => Some((*start, *end, reverse)),
            collagene_model::Location::Complement(inner) => last(inner, !reverse),
            collagene_model::Location::Join(parts) => last(
                if reverse {
                    parts.first()?
                } else {
                    parts.last()?
                },
                reverse,
            ),
            _ => None,
        }
    }
    exact(location).then(|| last(location, false)).flatten()
}

pub mod crispr;
pub mod genome;
pub mod genome_session;
pub mod view;
