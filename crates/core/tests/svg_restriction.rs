//! Independent literal coordinate fixtures; no scan/digest algorithm is used to
//! construct the expected geometry or evidence passed to the vector renderer.
use collagene_core::{
    Error, Molecule, Topology,
    restriction_inspection::{InspectionSite, ReferenceSpan},
    svg::*,
};

fn molecule(topology: Topology) -> Molecule {
    Molecule::dna("Literal restriction map", "AAAAGAATTCAAAAAA", topology).unwrap()
}
fn site() -> InspectionSite {
    InspectionSite {
        index: 7,
        enzyme_index: 0,
        enzyme: "EcoRI".into(),
        recognition_start: 4,
        recognition_length: 6,
        reverse: false,
        certain: true,
        spans: vec![ReferenceSpan { start: 4, end: 10 }],
        top_boundary: 5,
        bottom_boundary: 9,
        top_cut: Some(5),
        bottom_cut: Some(9),
        diagnostics: vec![],
    }
}

#[test]
fn linear_recognition_and_both_strand_boundaries_use_existing_axis_coordinates() {
    let m = molecule(Topology::Linear);
    let before = m.clone();
    let text = render_map_with_restriction_sites(&m, &m.revision(), &[site()]).unwrap();
    let document = roxmltree::Document::parse(&text).unwrap();
    let group = document
        .descendants()
        .find(|n| n.attribute("data-restriction-index") == Some("7"))
        .unwrap();
    let cut = |strand| {
        group
            .children()
            .find(|n| n.attribute("data-restriction-cut") == Some(strand))
            .unwrap()
    };
    assert_eq!(cut("top").attribute("data-position"), Some("5"));
    assert_eq!(
        cut("top").attribute("d"),
        Some("M 342.500 172.000 L 342.500 183.000")
    );
    assert_eq!(cut("bottom").attribute("data-position"), Some("9"));
    assert_eq!(
        cut("bottom").attribute("d"),
        Some("M 552.500 184.000 L 552.500 196.000")
    );
    let recognition = group
        .children()
        .find(|n| n.attribute("data-restriction-recognition") == Some("true"))
        .unwrap();
    assert_eq!(recognition.attribute("data-start"), Some("4"));
    assert_eq!(recognition.attribute("data-end"), Some("10"));
    assert_eq!(
        recognition.attribute("d"),
        Some("M 290.000 168.000 L 605.000 168.000")
    );
    assert_eq!(m, before);
    assert!(text.len() <= MAX_SVG_BYTES);
    assert!(!render_map(&m).unwrap().contains("data-restriction-index"));
}

#[test]
fn circular_origin_split_and_possible_reverse_sites_preserve_optional_cuts() {
    let m = molecule(Topology::Circular);
    let mut s = site();
    s.recognition_start = 14;
    s.recognition_length = 4;
    s.spans = vec![
        ReferenceSpan { start: 14, end: 16 },
        ReferenceSpan { start: 0, end: 2 },
    ];
    s.reverse = true;
    s.certain = false;
    s.top_boundary = 16;
    s.top_cut = Some(0);
    s.bottom_boundary = 18;
    s.bottom_cut = None;
    let text = render_map_with_restriction_sites(&m, &m.revision(), &[s]).unwrap();
    let document = roxmltree::Document::parse(&text).unwrap();
    let group = document
        .descendants()
        .find(|n| n.attribute("data-restriction-index") == Some("7"))
        .unwrap();
    assert_eq!(group.attribute("data-reverse"), Some("true"));
    assert_eq!(group.attribute("data-certain"), Some("false"));
    let ranges: Vec<_> = group
        .children()
        .filter(|n| n.has_attribute("data-restriction-recognition"))
        .map(|n| {
            (
                n.attribute("data-start").unwrap(),
                n.attribute("data-end").unwrap(),
            )
        })
        .collect();
    assert_eq!(ranges, vec![("14", "16"), ("0", "2")]);
    let cuts: Vec<_> = group
        .children()
        .filter(|n| n.has_attribute("data-restriction-cut"))
        .collect();
    assert_eq!(cuts.len(), 1);
    assert_eq!(cuts[0].attribute("data-restriction-cut"), Some("top"));
    assert_eq!(
        cuts[0].attribute("d"),
        Some("M 500.000 166.000 L 500.000 152.000")
    );
    assert_eq!(cuts[0].attribute("stroke-dasharray"), Some("4 3"));
    let title = group
        .children()
        .find(|n| n.has_tag_name("title"))
        .unwrap()
        .text()
        .unwrap();
    assert!(
        title.contains("15–16, 1–2")
            && title.contains("possible recognition")
            && title.contains("unavailable (raw boundary 18)")
    );
}

#[test]
fn stale_invalid_and_duplicate_site_evidence_is_rejected() {
    let m = molecule(Topology::Linear);
    assert!(render_map_with_restriction_sites(&m, &"0".repeat(64), &[site()]).is_err());
    assert!(render_map_with_restriction_sites(&m, &m.revision(), &[site(), site()]).is_err());
    for change in 0..4 {
        let mut s = site();
        match change {
            0 => s.spans[0].end = 11,
            1 => s.top_cut = Some(17),
            2 => s.top_boundary = 6,
            _ => s.recognition_start = 16,
        }
        assert!(render_map_with_restriction_sites(&m, &m.revision(), &[s]).is_err());
    }
    assert!(matches!(
        render_map_with_restriction_sites(
            &m,
            &m.revision(),
            &vec![site(); MAX_RESTRICTION_SITES + 1]
        ),
        Err(Error::TooLarge)
    ));
}

#[test]
fn hostile_names_are_only_xml_text_and_full_names_remain_accessible() {
    let m = molecule(Topology::Linear);
    let mut s = site();
    s.enzyme = "E<&\"'><script>alert(1)</script>\u{1}".into();
    let text = render_map_with_restriction_sites(&m, &m.revision(), &[s]).unwrap();
    let document = roxmltree::Document::parse(&text).unwrap();
    assert!(!document.descendants().any(|n| n.has_tag_name("script")));
    assert!(!text.contains("\u{1}"));
    assert!(
        document
            .descendants()
            .filter(|n| n.has_tag_name("title"))
            .any(|n| n
                .text()
                .is_some_and(|t| t.contains("E<&\"'><script>alert(1)</script>�")))
    );
    assert!(!text.contains("href=") && !text.contains("onload=") && !text.contains("<script>"));
}

#[test]
fn labels_are_evenly_sampled_without_dropping_markers_and_output_cap_fails_closed() {
    let m = molecule(Topology::Circular);
    let sites: Vec<_> = (0..100)
        .map(|index| InspectionSite { index, ..site() })
        .collect();
    let text = render_map_with_restriction_sites(&m, &m.revision(), &sites).unwrap();
    let document = roxmltree::Document::parse(&text).unwrap();
    assert_eq!(
        document
            .descendants()
            .filter(|n| n.has_attribute("data-restriction-index"))
            .count(),
        100
    );
    let labels: Vec<_> = document
        .descendants()
        .filter_map(|n| n.attribute("data-restriction-label-index"))
        .map(|s| s.parse::<usize>().unwrap())
        .collect();
    assert_eq!(labels.len(), 64);
    assert_eq!(labels[0], 0);
    assert_eq!(*labels.last().unwrap(), 99);
    assert!(labels.contains(&50));
    assert!(
        document
            .descendants()
            .any(|n| n.attribute("data-restriction-labels-omitted") == Some("36"))
    );
    let many: Vec<_> = (0..MAX_RESTRICTION_SITES)
        .map(|index| InspectionSite { index, ..site() })
        .collect();
    assert!(matches!(
        render_map_with_restriction_sites(&m, &m.revision(), &many),
        Err(Error::TooLarge)
    ));
    assert!(render_map(&m).is_ok());
}
