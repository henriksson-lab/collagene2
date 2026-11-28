use collagene_core::svg::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![("label".into(), Some("Synthetic feature".into()))],
        needs_review: false,
    }
}
fn molecule(topology: Topology) -> Molecule {
    Molecule::dna("Synthetic map", "AAAACCCCGGGGTTTT", topology).unwrap()
}

#[test]
fn origin_spanning_complements_preserve_biological_part_order_and_terminal_direction() {
    let mut m = molecule(Topology::Circular);
    let join = Location::Join(vec![Location::range(12, 16), Location::range(0, 4)]);
    m.features = vec![
        feature(join.clone()),
        feature(Location::Complement(Box::new(join))),
    ];
    let scene = layout_map(&m).unwrap();
    let forward = &scene.features[0];
    let reverse = &scene.features[1];
    assert_eq!(
        forward
            .spans
            .iter()
            .map(|s| (s.start, s.end, s.reverse))
            .collect::<Vec<_>>(),
        vec![(12, 16, false), (0, 4, false)]
    );
    assert_eq!(
        reverse
            .spans
            .iter()
            .map(|s| (s.start, s.end, s.reverse))
            .collect::<Vec<_>>(),
        vec![(0, 4, true), (12, 16, true)]
    );
    assert_eq!(forward.terminal_span, Some(1));
    assert_eq!(reverse.terminal_span, Some(1));
    let svg = render_map(&m).unwrap();
    let document = roxmltree::Document::parse(&svg).unwrap();
    let arrows: Vec<_> = document
        .descendants()
        .filter(|n| n.attribute("data-terminal").is_some())
        .collect();
    assert_eq!(arrows.len(), 2);
    assert_eq!(
        (
            arrows[0].attribute("data-terminal"),
            arrows[0].attribute("data-direction")
        ),
        (Some("4"), Some("forward"))
    );
    assert_eq!(
        (
            arrows[1].attribute("data-terminal"),
            arrows[1].attribute("data-direction")
        ),
        (Some("12"), Some("reverse"))
    );
    assert_eq!(
        document
            .descendants()
            .filter(|n| n.attribute("data-start").is_some())
            .count(),
        4
    );
    assert_eq!(
        document.root_element().attribute("data-revision"),
        Some(m.revision().as_str())
    );
}

#[test]
fn uncertain_ordered_and_unresolved_features_are_never_given_exact_arrows() {
    let mut m = molecule(Topology::Linear);
    let mut review = feature(Location::range(9, 12));
    review.needs_review = true;
    m.features = vec![
        feature(Location::Range {
            start: 0,
            end: 4,
            before: true,
            after: true,
        }),
        feature(Location::Order(vec![
            Location::range(5, 6),
            Location::range(7, 8),
        ])),
        feature(Location::Join(vec![
            Location::range(12, 14),
            Location::Unresolved("REMOTE:1..2".into()),
        ])),
        review,
    ];
    let scene = layout_map(&m).unwrap();
    assert!(scene.features[0].fuzzy);
    assert!(scene.features[1].ordered);
    assert!(scene.features[2].unresolved);
    assert!(
        scene.features[..3]
            .iter()
            .all(|f| f.terminal_span.is_none())
    );
    assert_eq!(scene.features[2].spans.len(), 1);
    assert!(scene.labels.iter().any(|label| label.text.ends_with('*')));
    assert!(scene.warnings.iter().any(|s| s.contains("unresolved")));
    let svg = render_map(&m).unwrap();
    let document = roxmltree::Document::parse(&svg).unwrap();
    assert_eq!(
        document
            .descendants()
            .filter(|n| n.attribute("data-terminal").is_some())
            .count(),
        1
    );
    assert!(
        document
            .descendants()
            .any(|n| n.attribute("stroke-dasharray") == Some("6 4"))
    );
    assert!(
        document
            .descendants()
            .any(|n| n.attribute("stroke") == Some("#b14848"))
    );
    assert!(!svg.contains("REMOTE:1..2"));
}

#[test]
fn svg_is_xml_safe_and_has_no_executable_or_external_markup() {
    let mut m = molecule(Topology::Circular);
    m.name = "\"><script>alert('x')</script>&\u{b}β".into();
    let mut f = feature(Location::range(0, 16));
    f.kind = "\" onload=\"evil()".into();
    f.qualifiers=vec![("label".into(),Some("<foreignObject><img src='https://example.invalid/x' onerror='evil()'/></foreignObject>&\u{0}東京".into()))];
    m.features.push(f);
    let svg = render_map(&m).unwrap();
    let document = roxmltree::Document::parse(&svg).unwrap();
    assert_eq!(
        document.root_element().tag_name().namespace(),
        Some("http://www.w3.org/2000/svg")
    );
    for node in document.descendants().filter(|n| n.is_element()) {
        assert!(!matches!(
            node.tag_name().name(),
            "script" | "foreignObject" | "image" | "img" | "a" | "style"
        ));
        for attribute in node.attributes() {
            assert!(!attribute.name().to_ascii_lowercase().starts_with("on"));
            assert!(!matches!(attribute.name(), "href" | "src"));
        }
    }
    assert!(!svg.contains('\u{0}') && !svg.contains('\u{b}'));
    assert!(svg.contains("&lt;script&gt;") && svg.contains("&amp;"));
    assert!(!svg.contains("<!DOCTYPE"));
}

#[test]
fn adjacent_half_open_ranges_can_share_a_track_but_overlaps_cannot() {
    let mut m = molecule(Topology::Linear);
    m.features = vec![
        feature(Location::range(0, 4)),
        feature(Location::range(4, 8)),
        feature(Location::range(2, 6)),
    ];
    let scene = layout_map(&m).unwrap();
    assert_eq!(scene.features[0].lane, scene.features[1].lane);
    assert_ne!(scene.features[0].lane, scene.features[2].lane);
    assert_eq!(scene.track_collisions, 0);
    let svg = render_map(&m).unwrap();
    let document = roxmltree::Document::parse(&svg).unwrap();
    let intervals: Vec<_> = document
        .descendants()
        .filter(|n| n.attribute("data-start").is_some())
        .map(|n| {
            (
                n.attribute("data-start").unwrap(),
                n.attribute("data-end").unwrap(),
            )
        })
        .collect();
    assert_eq!(intervals, vec![("0", "4"), ("4", "8"), ("2", "6")]);
}

#[test]
fn complete_circle_ranges_are_drawn_with_two_arcs_and_empty_maps_are_valid() {
    let mut m = molecule(Topology::Circular);
    assert!(roxmltree::Document::parse(&render_map(&m).unwrap()).is_ok());
    m.features.push(feature(Location::range(0, 16)));
    let svg = render_map(&m).unwrap();
    let document = roxmltree::Document::parse(&svg).unwrap();
    let full = document
        .descendants()
        .find(|n| n.attribute("data-start") == Some("0"))
        .unwrap();
    assert_eq!(full.attribute("data-end"), Some("16"));
    assert_eq!(full.attribute("d").unwrap().matches(" A ").count(), 2);
    let tiny = Molecule::dna("one base", "A", Topology::Linear).unwrap();
    let tiny_svg = render_map(&tiny).unwrap();
    let document = roxmltree::Document::parse(&tiny_svg).unwrap();
    assert_eq!(
        document
            .descendants()
            .filter(|n| n.tag_name().name() == "text" && n.text() == Some("1"))
            .count(),
        1
    );
}

#[test]
fn rna_and_protein_maps_have_correct_units_and_protein_has_no_strand_arrows() {
    let mut m = molecule(Topology::Linear);
    m.sequence = Sequence::new(Alphabet::Rna, "ACGUACGUACGUACGU").unwrap();
    assert!(render_map(&m).unwrap().contains("16 nt"));
    m.sequence = Sequence::new(Alphabet::Protein, "ACDEFGHIKLMNPQRS").unwrap();
    m.features
        .push(feature(Location::Complement(Box::new(Location::range(
            0, 4,
        )))));
    let scene = layout_map(&m).unwrap();
    assert!(scene.features[0].terminal_span.is_none());
    assert!(
        scene
            .warnings
            .iter()
            .any(|s| s.contains("Protein complement"))
    );
    let svg = render_map(&m).unwrap();
    assert!(svg.contains("16 aa") && !svg.contains("data-terminal="));
    assert!(svg.contains("duplex ends, nicks, methylation"));
    m.topology = Topology::Circular;
    assert!(render_map(&m).is_err());
}

#[test]
fn feature_span_and_label_limits_are_visible_and_deterministic() {
    let mut m = molecule(Topology::Circular);
    m.features = (0..MAX_MAP_FEATURES + 1)
        .map(|_| feature(Location::range(0, 16)))
        .collect();
    let scene = layout_map(&m).unwrap();
    assert_eq!(scene.features.len(), MAX_MAP_FEATURES);
    assert_eq!(scene.omitted_features, 1);
    assert!(scene.labels.len() <= MAX_MAP_LABELS);
    assert_eq!(
        scene.omitted_labels,
        scene.features.len() - scene.labels.len()
    );
    assert!(
        scene
            .warnings
            .iter()
            .any(|s| s.contains("additional features"))
    );
    assert_eq!(render_map(&m).unwrap(), render_map(&m).unwrap());
    assert!(render_map(&m).unwrap().len() < MAX_SVG_BYTES);
    let mut full = molecule(Topology::Linear);
    full.features.push(feature(Location::Join(vec![
        Location::range(0, 1);
        MAX_MAP_SPANS
    ])));
    let exact = layout_map(&full).unwrap();
    assert!(!exact.spans_omitted);
    assert_eq!(exact.features[0].spans.len(), MAX_MAP_SPANS);
    full.features[0].location = Location::Join(vec![Location::range(0, 1); MAX_MAP_SPANS + 1]);
    let partial = layout_map(&full).unwrap();
    assert!(partial.spans_omitted);
    assert!(partial.features[0].partial && partial.features[0].terminal_span.is_none());
    assert_eq!(partial.features[0].spans.len(), MAX_MAP_SPANS);
    assert!(
        partial
            .warnings
            .iter()
            .any(|s| s.contains("compound spans"))
    );
}

#[test]
fn malformed_coordinates_fail_before_rendering_and_labels_stay_inside_plot() {
    let mut m = molecule(Topology::Circular);
    m.features = (0..100)
        .map(|i| feature(Location::range(i % 16, (i % 16) + 1)))
        .collect();
    let scene = layout_map(&m).unwrap();
    for label in &scene.labels {
        assert!(
            label.x.is_finite()
                && label.y.is_finite()
                && label.anchor_x.is_finite()
                && label.anchor_y.is_finite()
        );
        assert!((100.0..=600.0).contains(&label.y));
    }
    for right in [false, true] {
        let mut ys: Vec<_> = scene
            .labels
            .iter()
            .filter(|l| l.right_aligned == right)
            .map(|l| l.y)
            .collect();
        ys.sort_by(f64::total_cmp);
        assert!(ys.windows(2).all(|p| p[1] - p[0] >= 23.99));
    }
    m.features[0].location = Location::range(0, 17);
    assert!(render_map(&m).is_err());
}
