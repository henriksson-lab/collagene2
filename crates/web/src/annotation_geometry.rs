//! Bounded local annotation geometry shared by static display caches.
//! Coordinates are source ranges, not biological extraction/validity claims.
use collagene_core::{Location, Molecule};

pub(crate) const MAX_FEATURES: usize = collagene_view_layout::MAX_ITEMS;
pub(crate) const MAX_SPANS: usize = collagene_view_layout::MAX_INTERVALS;
pub(crate) const MAX_NODES: usize = 1_000_000;
const MAX_DEPTH: usize = 32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocatedFeature {
    pub index: usize,
    pub primer: bool,
    /// Local ranges in source-tree order, with nested complement parity.
    pub spans: Vec<(usize, usize, bool)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct PreparedGeometry {
    pub located: Vec<LocatedFeature>,
    /// Non-source features having no known local Range part.
    pub unplaced: Vec<usize>,
}

#[derive(Default)]
struct Counts {
    nodes: usize,
    spans: usize,
}

fn preflight_location(
    location: &Location,
    length: usize,
    depth: usize,
    counts: &mut Counts,
) -> Result<(), String> {
    if depth >= MAX_DEPTH {
        return Err("Annotation layout location depth exceeds 32 levels".into());
    }
    counts.nodes = counts
        .nodes
        .checked_add(1)
        .ok_or("Annotation layout node count overflow")?;
    if counts.nodes > MAX_NODES {
        return Err("Annotation layout node limit exceeded".into());
    }
    match location {
        Location::Range { start, end, .. } => {
            counts.spans = counts
                .spans
                .checked_add(1)
                .ok_or("Annotation layout span count overflow")?;
            if counts.spans > MAX_SPANS {
                return Err("Annotation layout span limit exceeded".into());
            }
            if start >= end || *end > length {
                return Err(
                    "Annotation layout range is empty, reversed or outside the source".into(),
                );
            }
        }
        Location::Complement(inner) => preflight_location(inner, length, depth + 1, counts)?,
        Location::Join(parts) | Location::Order(parts) => {
            // Reject a huge flat child list before iterating it. No intermediate
            // .spans(), cloned location tree, or Items vector is constructed.
            if parts.len() > MAX_NODES.saturating_sub(counts.nodes) {
                return Err("Annotation layout node limit exceeded".into());
            }
            for part in parts {
                preflight_location(part, length, depth + 1, counts)?;
            }
        }
        Location::Unresolved(_) => {}
    }
    Ok(())
}

fn collect_spans(location: &Location, reverse: bool, spans: &mut Vec<(usize, usize, bool)>) {
    // The entire borrowed tree has passed all resource/coordinate bounds.
    match location {
        Location::Range { start, end, .. } => spans.push((*start, *end, reverse)),
        Location::Complement(inner) => collect_spans(inner, !reverse, spans),
        Location::Join(parts) | Location::Order(parts) => {
            for part in parts {
                collect_spans(part, reverse, spans);
            }
        }
        Location::Unresolved(_) => {}
    }
}

/// Copy only bounded known local annotation ranges. No DNA hashing, extraction,
/// label/qualifier cloning or molecule validation is performed. Fuzzy endpoints,
/// joined/order constructs and partial unresolved locations are display evidence,
/// never a claim of biological validity. All source feature indices are retained;
/// case-insensitive source annotations are omitted from display.
///
/// Complete borrowed feature/node/span/depth/range preflight runs BEFORE any
/// output vector or span allocation. Resource counts include ignored source
/// annotations, preserving the former annotation-layout error contract.
pub fn prepare(molecule: &Molecule) -> Result<PreparedGeometry, String> {
    if molecule.features.len() > MAX_FEATURES {
        return Err("Annotation layout feature limit exceeded".into());
    }
    let mut counts = Counts::default();
    for feature in &molecule.features {
        preflight_location(&feature.location, molecule.sequence.len(), 0, &mut counts)?;
    }
    let mut result = PreparedGeometry::default();
    for (index, feature) in molecule.features.iter().enumerate() {
        if feature.kind.eq_ignore_ascii_case("source") {
            continue;
        }
        let mut spans = Vec::new();
        collect_spans(&feature.location, false, &mut spans);
        if spans.is_empty() {
            result.unplaced.push(index);
        } else {
            result.located.push(LocatedFeature {
                index,
                primer: feature.kind.eq_ignore_ascii_case("primer_bind"),
                spans,
            });
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Topology, model::Feature};

    fn range(a: usize, b: usize) -> Location {
        Location::range(a, b)
    }
    fn feature(kind: &str, location: Location) -> Feature {
        Feature {
            kind: kind.into(),
            location,
            qualifiers: vec![("label".into(), Some("retained full label".into()))],
            needs_review: true,
        }
    }
    fn molecule(features: Vec<Feature>) -> Molecule {
        let mut m =
            Molecule::dna("literal geometry", &"ACGT".repeat(25), Topology::Circular).unwrap();
        m.features = features;
        m.original_genbank = Some("immutable original evidence".into());
        m.metadata
            .insert("fixture".into(), "literal nonmutation evidence".into());
        m
    }
    fn checked(m: &Molecule) -> PreparedGeometry {
        let before = m.clone();
        let result = prepare(m).unwrap();
        assert_eq!(*m, before);
        result
    }

    #[test]
    fn original_indices_classification_order_and_cardinality() {
        let m = molecule(vec![
            feature("SoUrCe", range(0, 100)),
            feature("CDS", range(80, 90)),
            feature("PrImEr_BiNd", range(1, 4)),
            feature("misc_feature", Location::Unresolved("remote:1..4".into())),
            feature("source", Location::Unresolved("ignored".into())),
            feature("misc_feature", Location::Join(vec![])),
        ]);
        let result = checked(&m);
        assert_eq!(
            result.located,
            vec![
                LocatedFeature {
                    index: 1,
                    primer: false,
                    spans: vec![(80, 90, false)]
                },
                LocatedFeature {
                    index: 2,
                    primer: true,
                    spans: vec![(1, 4, false)]
                }
            ]
        );
        assert_eq!(result.unplaced, vec![3, 5]);
        assert_eq!(result.located.len() + result.unplaced.len(), 4);
        assert_eq!(checked(&molecule(vec![])), PreparedGeometry::default());
    }

    #[test]
    fn complement_parity_join_order_partial_unresolved_and_fuzzy_spans() {
        let tree = Location::Complement(Box::new(Location::Join(vec![
            range(90, 100),
            Location::Unresolved("remote:1..4".into()),
            Location::Complement(Box::new(Location::Order(vec![
                range(20, 30),
                Location::Complement(Box::new(range(40, 50))),
            ]))),
            Location::Range {
                start: 0,
                end: 10,
                before: true,
                after: true,
            },
        ])));
        let result = checked(&molecule(vec![feature("CDS", tree)]));
        assert_eq!(
            result.located[0].spans,
            vec![
                (90, 100, true),
                (20, 30, false),
                (40, 50, true),
                (0, 10, true)
            ]
        );
        assert!(result.unplaced.is_empty());
    }

    #[test]
    fn invalid_range_and_depth_preserve_previous_error_contract() {
        for location in [range(1, 1), range(2, 1), range(90, 101)] {
            let m = molecule(vec![feature("CDS", location)]);
            let before = m.clone();
            assert!(prepare(&m).unwrap_err().contains("range"));
            assert_eq!(m, before);
        }
        let mut tree = range(0, 1);
        for _ in 0..31 {
            tree = Location::Complement(Box::new(tree));
        }
        assert!(prepare(&molecule(vec![feature("CDS", tree.clone())])).is_ok());
        tree = Location::Complement(Box::new(tree));
        assert!(
            prepare(&molecule(vec![feature("CDS", tree)]))
                .unwrap_err()
                .contains("depth")
        );
    }

    #[test]
    fn borrowed_aggregate_count_bounds_include_ignored_source_nodes() {
        let mut m = molecule(vec![]);
        m.features = vec![feature("source", Location::Unresolved(String::new())); MAX_FEATURES + 1];
        assert!(prepare(&m).unwrap_err().contains("feature limit"));
        m.features = vec![feature(
            "SoUrCe",
            Location::Join(vec![range(0, 1); MAX_SPANS + 1]),
        )];
        assert!(prepare(&m).unwrap_err().contains("span limit"));
        m.features = vec![feature(
            "source",
            Location::Order(vec![Location::Unresolved(String::new()); MAX_NODES]),
        )];
        assert!(prepare(&m).unwrap_err().contains("node limit"));
    }

    #[test]
    fn ten_thousand_features_keep_exact_ranges_and_source_evidence() {
        let mut m = Molecule::dna("1Mb", &"ACGT".repeat(250_000), Topology::Linear).unwrap();
        m.features = (0..10_000)
            .map(|i| {
                feature(
                    if i % 2 == 0 { "CDS" } else { "primer_bind" },
                    range(i * 100 + 10, i * 100 + 30),
                )
            })
            .collect();
        let result = checked(&m);
        assert_eq!(result.located.len(), 10_000);
        assert!(result.unplaced.is_empty());
        for (i, located) in result.located.iter().enumerate() {
            assert_eq!(located.index, i);
            assert_eq!(located.primer, i % 2 == 1);
            assert_eq!(located.spans, vec![(i * 100 + 10, i * 100 + 30, false)]);
        }
    }
}
