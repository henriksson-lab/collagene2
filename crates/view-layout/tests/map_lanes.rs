//! Map tracks use this same cached whole-feature packing, never annotation ordinal modulo.
use collagene_view_layout::*;
use std::collections::BTreeMap;
fn item(id: u64, spans: &[(u64, u64)]) -> Item {
    Item {
        id,
        spans: spans
            .iter()
            .map(|&(start, end)| Interval { start, end })
            .collect(),
    }
}
fn map(items: &[Item]) -> Layout {
    pack(
        items,
        100,
        &Limits {
            max_lanes: 3,
            ..Limits::default()
        },
    )
    .unwrap()
}
fn lanes(layout: Layout) -> BTreeMap<u64, Placement> {
    layout
        .assignments
        .into_iter()
        .map(|a| (a.id, a.placement))
        .collect()
}
#[test]
fn separated_genes_share_the_backbone_track_regardless_of_source_annotation_order() {
    let items = vec![
        item(90, &[(70, 90)]),
        item(2, &[(0, 10)]),
        item(35, &[(40, 55)]),
        item(7, &[(10, 35)]),
    ];
    let a = lanes(map(&items));
    assert!(a.values().all(|p| *p == Placement::Lane(0)));
    let mut reverse = items;
    reverse.reverse();
    assert_eq!(a, lanes(map(&reverse)));
}
#[test]
fn actual_overlap_separates_tracks_and_excess_overlap_withholds_the_whole_feature() {
    let items = vec![
        item(40, &[(20, 60)]),
        item(20, &[(20, 60)]),
        item(30, &[(20, 60)]),
        item(10, &[(20, 60)]),
    ];
    let result = map(&items);
    assert_eq!(result.lane_count, 3);
    let a = lanes(result);
    assert_eq!(a[&10], Placement::Lane(0));
    assert_eq!(a[&20], Placement::Lane(1));
    assert_eq!(a[&30], Placement::Lane(2));
    assert_eq!(a[&40], Placement::Overflow);
    // Overflow does not occupy any lane or interfere with a later disjoint gene.
    let mut with_gap = items;
    with_gap.push(item(50, &[(65, 80)]));
    assert_eq!(lanes(map(&with_gap))[&50], Placement::Lane(0));
}
#[test]
fn origin_split_and_compound_parts_keep_one_lane_and_do_not_reserve_their_gaps() {
    let items = vec![
        item(1, &[(90, 100), (0, 10)]),
        item(2, &[(5, 15)]),
        item(3, &[(35, 45), (55, 65)]),
        item(4, &[(46, 54)]),
    ];
    let result = map(&items);
    let a = lanes(result.clone());
    assert_eq!(a[&1], Placement::Lane(0));
    assert_eq!(a[&2], Placement::Lane(1));
    assert_eq!(a[&3], Placement::Lane(0));
    assert_eq!(a[&4], Placement::Lane(0));
    let origin = result.assignments.iter().find(|r| r.id == 1).unwrap();
    assert_eq!(
        origin.spans,
        vec![
            Interval { start: 0, end: 10 },
            Interval {
                start: 90,
                end: 100
            }
        ]
    );
}
