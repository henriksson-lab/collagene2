use collagene_core::tracks::FeatureIndex;
use collagene_model::{Feature, Location, Molecule, Topology};
fn f(kind: &str, location: Location) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    }
}
#[test]
fn interval_queries_match_independent_original_order_scan_for_compound_and_long_features() {
    let mut source = Molecule::dna("index", &"A".repeat(300), Topology::Circular).unwrap();
    source.features.push(f("source", Location::range(0, 300)));
    source
        .features
        .push(f("misc_feature", Location::range(1, 299)));
    for i in 0..120 {
        let start = i * 71 % 280;
        source.features.push(f(
            "misc_feature",
            Location::Complement(Box::new(Location::Join(vec![
                Location::range(start, start + 10),
                Location::range((start + 20) % 280, (start + 20) % 280 + 10),
            ]))),
        ));
    }
    source.features.push(f(
        "misc_feature",
        Location::Unresolved("remote:1..4".into()),
    ));
    let index = FeatureIndex::new(&source).unwrap();
    for start in 0..300 {
        for width in [1, 5, 20, 300] {
            for limit in [0, 1, 2, 10] {
                let end = (start + width).min(300);
                let expected = source
                    .features
                    .iter()
                    .enumerate()
                    .filter(|(_, f)| {
                        f.kind != "source"
                            && f.location
                                .spans()
                                .iter()
                                .any(|(a, b, _)| *a < end && *b > start)
                    })
                    .map(|(i, _)| i)
                    .take(limit)
                    .collect::<Vec<_>>();
                assert_eq!(index.visible(start, end, limit), expected);
            }
        }
    }
    assert!(index.visible(20, 20, 2).is_empty());
}
#[test]
fn viewport_at_far_end_does_not_lose_annotations_beyond_a_fixed_prefix() {
    let mut source = Molecule::dna("dense", &"A".repeat(1_000_000), Topology::Linear).unwrap();
    source.features = (0..10_000)
        .map(|i| f("misc_feature", Location::range(i * 100, i * 100 + 50)))
        .collect();
    let index = FeatureIndex::new(&source).unwrap();
    assert_eq!(index.visible(999_900, 999_950, 2), vec![9999]);
    assert_eq!(index.visible(999_950, 1_000_000, 2), Vec::<usize>::new());
}

#[test]
fn directional_heads_follow_biological_order_and_unknown_strands_have_no_arrow_in_either_scene() {
    use collagene_core::tracks::terminal_span;
    let ordered = Location::Join(vec![Location::range(8, 12), Location::range(0, 3)]);
    let reversed = Location::Complement(Box::new(ordered.clone()));
    let fuzzy = Location::Range {
        start: 2,
        end: 5,
        before: true,
        after: false,
    };
    let cases = vec![
        (ordered, Some((0, 3, false))),
        (reversed, Some((8, 12, true))),
        (Location::Order(vec![Location::range(2, 5)]), None),
        (
            Location::Complement(Box::new(Location::Order(vec![Location::range(2, 5)]))),
            None,
        ),
        (fuzzy, None),
        (
            Location::Join(vec![
                Location::range(0, 3),
                Location::Unresolved("remote:4..6".into()),
            ]),
            None,
        ),
    ];
    for (location, expected) in cases {
        assert_eq!(terminal_span(&location), expected);
        let mut molecule = Molecule::dna("direction", "ACGTACGTACGT", Topology::Circular).unwrap();
        molecule.features.push(f("misc_feature", location));
        let scene = collagene_core::svg::layout_map(&molecule).unwrap();
        let displayed = scene.features[0].terminal_span.map(|index| {
            let span = &scene.features[0].spans[index];
            (span.start, span.end, span.reverse)
        });
        assert_eq!(displayed, expected);
    }
}
