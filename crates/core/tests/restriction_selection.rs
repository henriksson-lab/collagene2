use collagene_core::{
    Location, Molecule, Topology,
    analysis::{self, Enzyme},
    restriction_inspection::{
        self as inspection, InspectionOptions, InspectionSite, ReferenceSpan, SiteDiagnostic,
    },
    restriction_selection::{self as selection, Basis::*, MAX_REGION_SPANS, Region, Relation::*},
};

fn enzyme(name: &str) -> Enzyme {
    analysis::starter_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .unwrap()
}
fn sites(sequence: &str, topology: Topology, enzyme: Enzyme) -> Vec<InspectionSite> {
    let m = Molecule::dna("selection fixture", sequence, topology).unwrap();
    inspection::scan("source", &m, &[enzyme], &InspectionOptions::default())
        .unwrap()
        .sites
}
fn linear_site() -> InspectionSite {
    sites("TTTGAATTCAAA", Topology::Linear, enzyme("EcoRI")).remove(0)
}

#[test]
fn recognition_and_type_iis_cuts_are_distinct_geometry() {
    let site = sites("AAAGGTCTCTAGTCCCCCAA", Topology::Linear, enzyme("BsaI")).remove(0);
    let footprint = Region::new(20, [(3, 9)]).unwrap();
    assert_eq!(
        selection::classify(&site, &footprint, &RecognitionBases).unwrap(),
        Inside
    );
    assert_eq!(
        selection::classify(&site, &footprint, &StrandCuts).unwrap(),
        Outside
    );
    let cuts = Region::new(20, [(10, 14)]).unwrap();
    assert_eq!(
        selection::classify(&site, &cuts, &RecognitionBases).unwrap(),
        Outside
    );
    assert_eq!(
        selection::classify(&site, &cuts, &StrandCuts).unwrap(),
        Inside
    );
}

#[test]
fn union_overlap_counts_each_base_once_and_preserves_crossing() {
    let site = linear_site();
    let all = Region::new(12, [(5, 9), (3, 6), (3, 4)]).unwrap();
    assert_eq!(all.spans(), &[(3, 9)]);
    assert_eq!(
        selection::classify(&site, &all, &RecognitionBases).unwrap(),
        Inside
    );
    let partial = Region::new(12, [(0, 5), (8, 12)]).unwrap();
    assert_eq!(
        selection::classify(&site, &partial, &RecognitionBases).unwrap(),
        Crossing
    );
    let outside = Region::new(12, [(0, 3), (9, 12)]).unwrap();
    assert_eq!(
        selection::classify(&site, &outside, &RecognitionBases).unwrap(),
        Outside
    );
    let adjacent = Region::new(12, [(6, 9), (3, 6)]).unwrap();
    assert_eq!(adjacent.spans(), &[(3, 9)]);
}

#[test]
fn cuts_include_both_selection_endpoints() {
    let site = linear_site();
    assert_eq!((site.top_cut, site.bottom_cut), (Some(4), Some(8)));
    assert_eq!(
        selection::classify(&site, &Region::new(12, [(4, 8)]).unwrap(), &StrandCuts).unwrap(),
        Inside
    );
    assert_eq!(
        selection::classify(&site, &Region::new(12, [(4, 7)]).unwrap(), &StrandCuts).unwrap(),
        Crossing
    );
    assert_eq!(
        selection::classify(&site, &Region::new(12, [(5, 7)]).unwrap(), &StrandCuts).unwrap(),
        Outside
    );
}

#[test]
fn circular_origin_wrapping_footprint_and_selection() {
    let site = sites("AATTCAAG", Topology::Circular, enzyme("EcoRI")).remove(0);
    assert_eq!(site.recognition_start, 7);
    assert_eq!(site.spans.len(), 2);
    let region = Region::circular(8, [(7, 8), (0, 5)]).unwrap();
    assert_eq!(region.spans(), &[(0, 5), (7, 8)]);
    assert_eq!(
        selection::classify(&site, &region, &RecognitionBases).unwrap(),
        Inside
    );
    // EcoRI top cut is exactly the origin, included by the final base interval.
    assert_eq!(site.top_cut, Some(0));
    let flank = Region::circular(8, [(7, 8)]).unwrap();
    assert_eq!(
        selection::classify(&site, &flank, &StrandCuts).unwrap(),
        Crossing
    );
    assert_eq!(
        selection::classify(&site, &flank, &RecognitionBases).unwrap(),
        Crossing
    );
    assert!(
        selection::classify(&site, &Region::new(8, [(7, 8)]).unwrap(), &RecognitionBases).is_err()
    );
}

#[test]
fn circular_full_reference_covers_every_cut_and_site() {
    let site = sites("AATTCAAG", Topology::Circular, enzyme("EcoRI")).remove(0);
    let region = Region::circular(8, [(0, 8)]).unwrap();
    for basis in [RecognitionBases, StrandCuts] {
        assert_eq!(selection::classify(&site, &region, &basis).unwrap(), Inside);
    }
}

#[test]
fn outside_linear_and_missing_cuts_are_undetermined() {
    let edge = sites("AAAGGTCTC", Topology::Linear, enzyme("BsaI")).remove(0);
    let region = Region::new(9, [(0, 9)]).unwrap();
    assert_eq!(
        selection::classify(&edge, &region, &RecognitionBases).unwrap(),
        Inside
    );
    assert_eq!(
        selection::classify(&edge, &region, &StrandCuts).unwrap(),
        Undetermined
    );
    let mut missing = linear_site();
    missing.top_cut = None;
    assert_eq!(
        selection::classify(&missing, &Region::new(12, [(0, 12)]).unwrap(), &StrandCuts).unwrap(),
        Undetermined
    );
}

#[test]
fn possible_recognition_never_becomes_certain() {
    let site = sites("TTTGAANTCAAA", Topology::Linear, enzyme("EcoRI")).remove(0);
    assert!(!site.certain);
    assert!(
        site.diagnostics
            .contains(&SiteDiagnostic::PossibleRecognition)
    );
    let before = serde_json::to_string(&site).unwrap();
    for basis in [RecognitionBases, StrandCuts] {
        assert_eq!(
            selection::classify(&site, &Region::new(12, [(0, 12)]).unwrap(), &basis).unwrap(),
            Inside
        );
    }
    assert_eq!(serde_json::to_string(&site).unwrap(), before);
}

#[test]
fn ambiguous_declared_strand_geometry_is_not_a_cut_claim() {
    let mut unusual = enzyme("EcoRI");
    unusual.bottom_cut = 6;
    let site = sites("TTTGAATTCAAA", Topology::Linear, unusual).remove(0);
    assert!(
        site.diagnostics
            .contains(&SiteDiagnostic::AmbiguousStrandGeometry)
    );
    let region = Region::new(12, [(0, 12)]).unwrap();
    assert_eq!(
        selection::classify(&site, &region, &RecognitionBases).unwrap(),
        Inside
    );
    assert_eq!(
        selection::classify(&site, &region, &StrandCuts).unwrap(),
        Undetermined
    );
}

#[test]
fn no_sites_and_partial_scan_do_not_create_regional_completeness() {
    assert!(sites("AAAACCCC", Topology::Linear, enzyme("EcoRI")).is_empty());
    let m = Molecule::dna("partial", "TTTGAATTCAAATTTGAATTCAAA", Topology::Linear).unwrap();
    let options = InspectionOptions {
        max_sites: 1,
        ..Default::default()
    };
    let report = inspection::scan("source", &m, &[enzyme("EcoRI")], &options).unwrap();
    assert!(!report.scan_complete);
    let region = Region::new(m.sequence.len(), [(0, m.sequence.len())]).unwrap();
    assert_eq!(
        selection::classify(&report.sites[0], &region, &RecognitionBases).unwrap(),
        Inside
    );
    assert!(!report.scan_complete);
}

#[test]
fn invalid_regions_and_chunk_bounds_are_rejected() {
    assert!(Region::new(12, []).is_err());
    for spans in [vec![(3, 3)], vec![(6, 2)], vec![(0, 13)]] {
        assert!(Region::new(12, spans).is_err());
    }
    assert!(Region::new(0, [(0, 1)]).is_err());
    assert!(Region::new(usize::MAX, [(0, 1)]).is_err());
    assert!(Region::new(12, std::iter::repeat_n((0, 1), MAX_REGION_SPANS + 1)).is_err());
    assert_eq!(
        Region::new(12, std::iter::repeat_n((0, 1), MAX_REGION_SPANS))
            .unwrap()
            .spans(),
        &[(0, 1)]
    );
}

#[test]
fn malformed_borrowed_site_is_rejected_before_classification() {
    let region = Region::new(12, [(0, 12)]).unwrap();
    let original = linear_site();
    let mut cases = Vec::new();
    let mut s = original.clone();
    s.spans.clear();
    cases.push(s);
    let mut s = original.clone();
    s.spans = vec![ReferenceSpan { start: 3, end: 9 }; 3];
    cases.push(s);
    let mut s = original.clone();
    s.spans[0].end = usize::MAX;
    cases.push(s);
    let mut s = original.clone();
    s.recognition_start = usize::MAX;
    cases.push(s);
    let mut s = original.clone();
    s.recognition_length = 0;
    cases.push(s);
    let mut s = original.clone();
    s.recognition_length = usize::MAX;
    cases.push(s);
    let mut s = original.clone();
    s.top_cut = Some(13);
    cases.push(s);
    let mut s = original.clone();
    s.top_cut = Some(3);
    cases.push(s);
    let mut s = original;
    s.diagnostics = vec![SiteDiagnostic::PossibleRecognition; 4];
    cases.push(s);
    for site in cases {
        for basis in [RecognitionBases, StrandCuts] {
            assert!(selection::classify(&site, &region, &basis).is_err());
        }
    }
}

#[test]
fn exact_location_footprint_preserves_parts_then_normalizes_union() {
    let location = Location::Complement(Box::new(Location::Join(vec![
        Location::range(8, 12),
        Location::range(0, 4),
        Location::range(2, 5),
    ])));
    let spans = selection::footprint(&location, 12).unwrap();
    assert_eq!(spans, vec![(8, 12), (0, 4), (2, 5)]);
    assert_eq!(
        Region::circular(12, spans).unwrap().spans(),
        &[(0, 5), (8, 12)]
    );
}

#[test]
fn footprint_rejects_fuzzy_unknown_empty_and_bounded_invalid_geometry() {
    let fuzzy = Location::Range {
        start: 1,
        end: 4,
        before: true,
        after: false,
    };
    for location in [
        fuzzy,
        Location::Order(vec![Location::range(0, 3)]),
        Location::Unresolved("1^2".into()),
        Location::Join(vec![]),
        Location::range(1, 1),
        Location::range(0, 13),
    ] {
        assert!(selection::footprint(&location, 12).is_err());
    }
    let mut deep = Location::range(0, 1);
    for _ in 0..33 {
        deep = Location::Complement(Box::new(deep));
    }
    assert!(selection::footprint(&deep, 12).is_err());
    let wide = Location::Join(vec![Location::range(0, 1); MAX_REGION_SPANS + 1]);
    assert!(selection::footprint(&wide, 12).is_err());
    assert!(selection::footprint(&Location::range(0, 1), 0).is_err());
}
