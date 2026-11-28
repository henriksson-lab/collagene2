//! Whole-source display lanes. This does not validate annotation biology, hash
//! source DNA, or derive layout from asynchronously returned viewport tracks.
#[cfg(test)]
use crate::annotation_geometry::MAX_NODES;
use crate::annotation_geometry::{MAX_FEATURES, MAX_SPANS};
#[cfg(test)]
use collagene_model::Location;
use collagene_model::Molecule;
use collagene_view_layout::{Interval, Item, Limits, Placement};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Group {
    pub lane_count: usize,
    pub placements: BTreeMap<usize, usize>,
    pub overflow: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct AnnotationLayout {
    pub features: Group,
    pub primers: Group,
    /// Non-source annotations with no known local Range part.
    pub unplaced: Vec<usize>,
}

fn group(items: &[Item], source_length: u64, cap: usize) -> Result<Group, String> {
    let limits = Limits {
        max_items: MAX_FEATURES,
        max_intervals: MAX_SPANS,
        max_lanes: cap,
        max_work: collagene_view_layout::MAX_WORK,
    };
    let packed =
        collagene_view_layout::pack(items, source_length, &limits).map_err(|e| e.to_string())?;
    let mut result = Group {
        lane_count: packed.lane_count,
        ..Group::default()
    };
    for assignment in packed.assignments {
        let id = usize::try_from(assignment.id)
            .map_err(|_| "Annotation layout feature ID exceeds usize")?;
        match assignment.placement {
            Placement::Lane(lane) => {
                result.placements.insert(id, lane);
            }
            Placement::Overflow => result.overflow.push(id),
        }
    }
    Ok(result)
}

/// Prepare one fixed whole-source layout, using original native feature indices.
/// Case-insensitive `source` annotations are ignored for display; `primer_bind`
/// annotations have a separate lane budget. Partial known local locations are
/// displayed without claiming unresolved parts or fuzzy endpoints are validated.
/// All feature/tree counts are checked on borrowed data before constructing any
/// span vectors / Items. Qualifiers, labels, sequence and metadata are untouched.
pub fn prepare(
    molecule: &Molecule,
    feature_cap: usize,
    primer_cap: usize,
) -> Result<AnnotationLayout, String> {
    prepare_visible(molecule, feature_cap, primer_cap, &BTreeSet::new())
}

/// Display preferences retain native feature IDs and do not consume lanes or
/// overflow slots. Hidden source trees still receive borrowed resource checks.
pub fn prepare_visible(
    molecule: &Molecule,
    feature_cap: usize,
    primer_cap: usize,
    hidden: &BTreeSet<usize>,
) -> Result<AnnotationLayout, String> {
    if feature_cap > collagene_view_layout::MAX_LANES
        || primer_cap > collagene_view_layout::MAX_LANES
    {
        return Err("Annotation layout lane cap exceeds the supported limit".into());
    }
    let geometry = crate::annotation_geometry::prepare(molecule)?;
    let mut features = Vec::new();
    let mut primers = Vec::new();
    for located in geometry
        .located
        .into_iter()
        .filter(|f| !hidden.contains(&f.index))
    {
        let item = Item {
            id: located.index as u64,
            spans: located
                .spans
                .into_iter()
                .map(|(start, end, _)| Interval {
                    start: start as u64,
                    end: end as u64,
                })
                .collect(),
        };
        if located.primer {
            primers.push(item);
        } else {
            features.push(item);
        }
    }
    Ok(AnnotationLayout {
        features: group(&features, molecule.sequence.len() as u64, feature_cap)?,
        primers: group(&primers, molecule.sequence.len() as u64, primer_cap)?,
        unplaced: geometry
            .unplaced
            .into_iter()
            .filter(|i| !hidden.contains(i))
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Feature, Topology};

    fn range(start: usize, end: usize) -> Location {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        }
    }
    fn feature(kind: &str, location: Location) -> Feature {
        Feature {
            kind: kind.into(),
            location,
            qualifiers: vec![("label".into(), Some("literal annotation".into()))],
            needs_review: false,
        }
    }
    fn molecule(features: Vec<Feature>) -> Molecule {
        let mut m = Molecule::dna("literal", &"ACGT".repeat(25), Topology::Circular).unwrap();
        m.features = features;
        m.original_genbank = Some("original evidence remains unchanged".into());
        m
    }
    fn checked(m: &Molecule, features: usize, primers: usize) -> AnnotationLayout {
        let before = m.clone();
        let layout = prepare(m, features, primers).unwrap();
        assert_eq!(*m, before);
        layout
    }

    #[test]
    fn none_single_and_distant_annotations_use_actual_lanes() {
        assert_eq!(
            checked(&molecule(vec![]), 3, 2),
            AnnotationLayout::default()
        );
        let m = molecule(vec![
            feature("misc_feature", range(1, 5)),
            feature("CDS", range(50, 60)),
        ]);
        let l = checked(&m, 3, 2);
        assert_eq!(l.features.lane_count, 1);
        assert_eq!(l.features.placements, BTreeMap::from([(0, 0), (1, 0)]));
        assert_eq!(l.primers, Group::default());
        let l = checked(&molecule(vec![feature("misc_feature", range(1, 5))]), 3, 2);
        assert_eq!(l.features.lane_count, 1);
    }

    #[test]
    fn overlap_overflow_and_zero_caps_retain_native_indices() {
        let m = molecule(vec![
            feature("source", range(0, 100)),
            feature("CDS", range(0, 10)),
            feature("misc_feature", range(1, 9)),
            feature("CDS", range(2, 8)),
        ]);
        let l = checked(&m, 2, 1);
        assert_eq!(l.features.placements, BTreeMap::from([(1, 0), (2, 1)]));
        assert_eq!(l.features.overflow, vec![3]);
        let l = checked(&m, 0, 0);
        assert_eq!(l.features.lane_count, 0);
        assert_eq!(l.features.overflow, vec![1, 2, 3]);
    }

    #[test]
    fn joined_circular_split_complement_order_and_partial_locations() {
        let m = molecule(vec![
            feature(
                "CDS",
                Location::Complement(Box::new(Location::Join(vec![range(90, 100), range(0, 10)]))),
            ),
            feature(
                "misc_feature",
                Location::Order(vec![range(10, 40), range(40, 90)]),
            ),
            feature(
                "misc_feature",
                Location::Join(vec![
                    Location::Unresolved("remote:1..4".into()),
                    range(20, 30),
                ]),
            ),
        ]);
        let l = checked(&m, 2, 0);
        assert_eq!(
            l.features.placements,
            BTreeMap::from([(0, 0), (1, 0), (2, 1)])
        );
        assert!(l.unplaced.is_empty());
    }

    #[test]
    fn primers_sources_and_wholly_unplaced_are_separate() {
        let m = molecule(vec![
            feature("SoUrCe", range(0, 100)),
            feature("PrImEr_BiNd", range(0, 10)),
            feature("CDS", range(0, 10)),
            feature("primer_bind", Location::Unresolved("remote:1..8".into())),
            feature("misc_feature", Location::Join(vec![])),
        ]);
        let l = checked(&m, 1, 1);
        assert_eq!(l.features.placements, BTreeMap::from([(2, 0)]));
        assert_eq!(l.primers.placements, BTreeMap::from([(1, 0)]));
        assert_eq!(l.unplaced, vec![3, 4]);
    }

    #[test]
    fn hidden_features_free_lanes_preserve_ids_and_do_not_count_as_overflow() {
        let m = molecule(vec![
            feature("CDS", range(0, 10)),
            feature("misc_feature", range(1, 9)),
            feature("primer_bind", range(0, 10)),
            feature("primer_bind", range(1, 9)),
            feature("misc_feature", Location::Unresolved("remote:1..3".into())),
        ]);
        let before = m.clone();
        let hidden = BTreeSet::from([0, 2, 4, usize::MAX]);
        let layout = prepare_visible(&m, 1, 1, &hidden).unwrap();
        assert_eq!(layout.features.placements, BTreeMap::from([(1, 0)]));
        assert_eq!(layout.primers.placements, BTreeMap::from([(3, 0)]));
        assert!(layout.features.overflow.is_empty());
        assert!(layout.primers.overflow.is_empty());
        assert!(layout.unplaced.is_empty());
        assert_eq!(m, before);
        let all = prepare_visible(&m, 1, 1, &BTreeSet::new()).unwrap();
        assert_eq!(all, prepare(&m, 1, 1).unwrap());
        assert_eq!(all.features.overflow, vec![1]);
        assert_eq!(all.primers.overflow, vec![3]);
        assert_eq!(all.unplaced, vec![4]);
        let none = prepare_visible(&m, 1, 1, &(0..m.features.len()).collect()).unwrap();
        assert_eq!(none, AnnotationLayout::default());
    }

    #[test]
    fn invalid_ranges_lane_caps_and_deep_trees_fail_without_mutation() {
        for location in [range(4, 4), range(5, 4), range(99, 101)] {
            let m = molecule(vec![feature("CDS", location)]);
            let before = m.clone();
            assert!(prepare(&m, 3, 2).is_err());
            assert_eq!(m, before);
        }
        let mut location = range(0, 1);
        for _ in 0..31 {
            location = Location::Complement(Box::new(location));
        }
        assert!(prepare(&molecule(vec![feature("CDS", location.clone())]), 3, 2).is_ok());
        location = Location::Complement(Box::new(location));
        assert!(
            prepare(&molecule(vec![feature("CDS", location)]), 3, 2)
                .unwrap_err()
                .contains("depth")
        );
        assert!(prepare(&molecule(vec![]), collagene_view_layout::MAX_LANES + 1, 0).is_err());
    }

    #[test]
    fn borrowed_count_caps_reject_before_items_or_span_cloning() {
        let mut m = molecule(vec![]);
        m.features = vec![feature("source", Location::Unresolved(String::new())); MAX_FEATURES + 1];
        assert!(prepare(&m, 3, 2).unwrap_err().contains("feature limit"));
        m.features = vec![feature(
            "CDS",
            Location::Join(vec![range(0, 1); MAX_SPANS + 1]),
        )];
        assert!(prepare(&m, 3, 2).unwrap_err().contains("span limit"));
        m.features = vec![feature(
            "CDS",
            Location::Join(vec![Location::Unresolved(String::new()); MAX_NODES]),
        )];
        assert!(prepare(&m, 3, 2).unwrap_err().contains("node limit"));
    }

    #[test]
    fn ten_thousand_annotations_on_one_megabase_keep_source_unchanged() {
        let mut m =
            Molecule::dna("one_megabase", &"ACGT".repeat(250_000), Topology::Linear).unwrap();
        m.features = (0..10_000)
            .map(|i| feature("misc_feature", range(i * 100 + 10, i * 100 + 30)))
            .collect();
        let l = checked(&m, 3, 2);
        assert_eq!(l.features.lane_count, 1);
        assert_eq!(l.features.placements.len(), 10_000);
        assert!(l.features.placements.values().all(|&lane| lane == 0));
        assert!(l.features.overflow.is_empty() && l.unplaced.is_empty());
    }
}
