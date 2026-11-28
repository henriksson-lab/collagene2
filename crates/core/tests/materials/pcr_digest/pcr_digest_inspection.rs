//! Literal primer-added restriction sites and the chemistry of their actual cuts.
use collagene_core::{
    chemistry::{FivePrime, ThreePrime},
    cloning::{pcr::digest as p, restriction::catalog as catalog_digest},
    sequence::analysis::Primer,
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Molecule, p::Request) {
    let source = Molecule::dna("Original6bp template", "AAACCC", Topology::Linear).unwrap();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let request = p::Request {
        template_strand: collagene_core::chemistry::Strand::Top,
        source_role: "original-template".into(),
        forward: Primer::new("F", "AAA", "CCGAATTC").unwrap(),
        reverse: Primer::new("R", "GGG", "GGGAATTC").unwrap(),
        selection: catalog_digest::select(&catalog, &[id]).unwrap(),
        options: Default::default(),
    };
    (source, request)
}
#[test]
fn unique_primer_added_sites_preserve_outer_oh_and_internal_cut_phosphate() {
    let (source, request) = fixture();
    let before = serde_json::to_vec(&source).unwrap();
    let preview = p::inspect(&source, &request).unwrap();
    preview
        .validate_binding(&p::inspection_revision(&source, &request).unwrap())
        .unwrap();
    assert_eq!(preview.pcr_geometries.len(), 1);
    let amplicon = preview.amplicon.as_ref().unwrap();
    // F tail8 + original6 + reverse-complement of R tail8 =22 bases.
    assert_eq!(amplicon.sequence.as_str(), "CCGAATTCAAACCCGAATTCCC");
    assert_eq!(preview.pcr_geometries[0].product_length, 22);
    let report = preview.inspection.as_ref().unwrap();
    assert!(report.scan_complete && report.products_determined);
    assert_eq!(
        report
            .sites
            .iter()
            .map(|s| s.recognition_start)
            .collect::<Vec<_>>(),
        vec![2, 14]
    );
    assert_eq!(report.fragments.len(), 3);
    let ends = report
        .fragments
        .iter()
        .map(|f| f.duplex.ends.as_ref().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(ends[0].0.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(ends[0].1.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[1].0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[1].1.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[2].0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[2].1.chemistry.five_prime, FivePrime::Hydroxyl);
    for (left, right) in ends {
        assert_eq!(left.chemistry.three_prime, ThreePrime::Hydroxyl);
        assert_eq!(right.chemistry.three_prime, ThreePrime::Hydroxyl);
    }
    assert_eq!(
        report.fragments[1].duplex.ends.as_ref().unwrap().0.overhang,
        "AATT"
    );
    assert_eq!(serde_json::to_vec(&source).unwrap(), before);
    let decoded: p::Preview =
        serde_json::from_slice(&serde_json::to_vec(&preview).unwrap()).unwrap();
    decoded.validate_binding(&preview.context_revision).unwrap();
}
#[test]
fn zero_and_competing_pcr_geometries_never_offer_digest_stock_or_choose_one() {
    let (_, request) = fixture();
    for (letters, minimum) in [("TTTTTT", 0), ("AAACCCAAACCC", 2)] {
        let source =
            Molecule::dna("Different literal template", letters, Topology::Linear).unwrap();
        let preview = p::inspect(&source, &request).unwrap();
        if minimum == 0 {
            assert!(preview.pcr_geometries.is_empty());
        } else {
            assert!(preview.pcr_geometries.len() >= minimum);
        }
        assert!(preview.amplicon.is_none() && preview.inspection.is_none());
        preview.validate().unwrap();
    }
}
#[test]
fn partial_scan_and_fragment_limits_withhold_the_entire_product_set() {
    let (source, mut request) = fixture();
    request.options.max_comparisons = 1;
    let preview = p::inspect(&source, &request).unwrap();
    assert!(preview.amplicon.is_some());
    let report = preview.inspection.unwrap();
    assert!(!report.scan_complete && !report.products_determined);
    assert!(report.fragments.is_empty());
    request.options.max_comparisons = 50_000_000;
    request.options.max_fragments = 2;
    let report = p::inspect(&source, &request).unwrap().inspection.unwrap();
    assert!(report.scan_complete && !report.products_determined);
    assert!(report.fragments.is_empty());
}
#[test]
fn stale_pins_primers_sources_conditions_and_imported_physical_state_are_not_authority() {
    let (mut source, mut request) = fixture();
    let original = p::inspection_revision(&source, &request).unwrap();
    let mut preview = p::inspect(&source, &request).unwrap();
    request.options.max_comparisons -= 1;
    assert_ne!(original, p::inspection_revision(&source, &request).unwrap());
    assert!(
        preview
            .validate_binding(&p::inspection_revision(&source, &request).unwrap())
            .is_err()
    );
    request.selection.catalog_content_hash = "0".repeat(64);
    assert!(p::source_preflight(&source, &request).is_err());
    let (_, request) = fixture();
    source.sequence = collagene_model::Sequence::dna("AAACCA").unwrap();
    assert_ne!(original, p::inspection_revision(&source, &request).unwrap());
    source.metadata.insert(
        "duplex_model".into(),
        "untrusted imported physical state".into(),
    );
    assert!(p::source_preflight(&source, &request).is_err());
    preview.request.forward.tail.push('A');
    assert!(preview.validate_binding(&original).is_err());
    let (source, mut request) = fixture();
    request.forward.tail = "A".repeat(1001);
    assert!(p::source_preflight(&source, &request).is_err());
}

#[test]
fn two_different_primer_added_enzymes_keep_the_literal_insert_and_distinct_sticky_ends() {
    let (source, mut request) = fixture();
    request.reverse = Primer::new("R-XhoI", "GGG", "GGCTCGAG").unwrap();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let ids = ["EcoRI", "XhoI"]
        .iter()
        .map(|name| {
            catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.name == *name)
                .unwrap()
                .id
                .clone()
        })
        .collect::<Vec<_>>();
    request.selection = catalog_digest::select(&catalog, &ids).unwrap();
    let preview = p::inspect(&source, &request).unwrap();
    preview
        .validate_binding(&p::inspection_revision(&source, &request).unwrap())
        .unwrap();
    // Literal F tail CCGAATTC + original AAACCC + RC(GGCTCGAG)=CTCGAGCC.
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCAAACCCCTCGAGCC"
    );
    let report = preview.inspection.as_ref().unwrap();
    assert!(report.scan_complete && report.products_determined);
    assert_eq!(report.sites.len(), 2);
    // EcoRI G^AATTC begins at2; XhoI C^TCGAG begins at14.
    assert_eq!(report.ref_top_lengths.as_deref(), Some(&[12, 7, 3][..]));
    assert_eq!(report.fragments.len(), 3);
    assert_eq!(
        report
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![3, 12, 7]
    );
    let (left, right) = report.fragments[1].duplex.ends.as_ref().unwrap();
    assert_eq!(left.overhang, "AATT");
    assert_eq!(right.overhang, "TCGA");
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(left.chemistry.three_prime, ThreePrime::Hydroxyl);
    assert_eq!(right.chemistry.three_prime, ThreePrime::Hydroxyl);
    let (outer_left, _) = report.fragments[0].duplex.ends.as_ref().unwrap();
    let (_, outer_right) = report.fragments[2].duplex.ends.as_ref().unwrap();
    assert_eq!(outer_left.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(outer_right.chemistry.five_prime, FivePrime::Hydroxyl);
}
