use collagene_core::sequence::annotations::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn molecule(sequence: &str, topology: Topology) -> Molecule {
    Molecule::dna("fixture", sequence, topology).unwrap()
}
fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("gene".into())),
            ("translation".into(), Some("MK".into())),
            ("codon_start".into(), Some("1".into())),
        ],
        needs_review: false,
    }
}
fn segment(a: usize, b: usize, d: usize, reverse: bool) -> MappingSegment {
    MappingSegment {
        source_start: a,
        source_end: b,
        destination_start: d,
        reverse,
    }
}

#[test]
fn reverse_complement_preserves_biological_feature_sequence_and_frame_assertions() {
    let mut source = molecule("CCCATGAAATAACCCGGG", Topology::Linear);
    source.features.push(feature(Location::range(3, 12)));
    let result = reverse_complement(&source).unwrap();
    assert_eq!(result.sequence.as_str(), "CCCGGGTTATTTCATGGG");
    assert_eq!(
        result.features[0].location,
        Location::Complement(Box::new(Location::range(6, 15)))
    );
    assert_eq!(
        result.features[0]
            .location
            .extract(&result.sequence)
            .unwrap(),
        "ATGAAATAA"
    );
    assert!(!result.features[0].needs_review);
    assert!(
        result.features[0]
            .qualifiers
            .contains(&("translation".into(), Some("MK".into())))
    );
    let twice = reverse_complement(&result).unwrap();
    assert_eq!(twice.sequence, source.sequence);
    assert_eq!(twice.features, source.features);
}

#[test]
fn reverse_compound_origin_order_is_biological_not_sorted_by_destination() {
    let mut source = molecule("ACGTCCTAGGATTCGA", Topology::Circular);
    source
        .features
        .push(feature(Location::Complement(Box::new(Location::Join(
            vec![Location::range(12, 16), Location::range(0, 4)],
        )))));
    let result = reverse_complement(&source).unwrap();
    assert_eq!(result.sequence.as_str(), "TCGAATCCTAGGACGT");
    assert_eq!(
        source.features[0]
            .location
            .extract(&source.sequence)
            .unwrap(),
        "ACGTTCGA"
    );
    assert_eq!(
        result.features[0]
            .location
            .extract(&result.sequence)
            .unwrap(),
        "ACGTTCGA"
    );
    assert!(result.features[0].location.spans().iter().all(|s| !s.2));
    assert!(!result.features[0].needs_review);
}

#[test]
fn origin_rotation_preserves_intact_features_but_linear_cut_marks_interruption() {
    let mut source = molecule("AAAACCCCGGGGTTTT", Topology::Circular);
    source.features.push(feature(Location::range(4, 10)));
    let rotated = rotate_origin(&source, 6).unwrap();
    assert_eq!(rotated.sequence.as_str(), "CCGGGGTTTTAAAACC");
    assert_eq!(
        rotated.features[0].location,
        Location::Join(vec![Location::range(14, 16), Location::range(0, 4)])
    );
    assert_eq!(
        rotated.features[0]
            .location
            .extract(&rotated.sequence)
            .unwrap(),
        "CCCCGG"
    );
    assert!(!rotated.features[0].needs_review);
    let cut = extract_region(&source, 6, 22).unwrap();
    assert_eq!(cut.sequence, rotated.sequence);
    assert_eq!(cut.topology, Topology::Linear);
    assert!(cut.features[0].needs_review);
    assert!(cut.features[0].location.extract(&cut.sequence).is_err());
    assert!(
        !cut.features[0]
            .qualifiers
            .iter()
            .any(|(k, _)| k == "translation" || k == "codon_start")
    );
    assert!(
        cut.features[0]
            .qualifiers
            .iter()
            .any(|(k, v)| k == "note" && v.as_ref().is_some_and(|v| v.contains("/translation=MK")))
    );
}

#[test]
fn wrapped_extraction_and_clipped_compounds_do_not_invent_complete_cds() {
    let mut source = molecule("AAAACCCCGGGGTTTT", Topology::Circular);
    source.features.push(feature(Location::Join(vec![
        Location::range(12, 16),
        Location::range(0, 4),
    ])));
    let product = extract_region(&source, 12, 20).unwrap();
    assert_eq!(product.sequence.as_str(), "TTTTAAAA");
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "TTTTAAAA"
    );
    assert!(!product.features[0].needs_review);
    let mut linear = molecule("AAACCCGGGTTT", Topology::Linear);
    linear.features.push(feature(Location::Join(vec![
        Location::range(0, 3),
        Location::range(6, 9),
    ])));
    let clipped = extract_region(&linear, 6, 12).unwrap();
    assert_eq!(clipped.sequence.as_str(), "GGGTTT");
    assert!(clipped.features[0].needs_review);
    assert!(
        clipped.features[0]
            .location
            .extract(&clipped.sequence)
            .is_err()
    );
    assert!(
        !clipped.features[0]
            .qualifiers
            .iter()
            .any(|(k, _)| k == "codon_start")
    );
}

#[test]
fn uncertain_and_unmapped_locations_are_never_promoted_to_exact_product_annotations() {
    let mut source = molecule("AAAACCCCGGGG", Topology::Linear);
    source.features.push(feature(Location::Order(vec![
        Location::range(0, 2),
        Location::range(4, 6),
    ])));
    source
        .features
        .push(feature(Location::Unresolved("REMOTE.1:1..3".into())));
    source.features.push(feature(Location::range(9, 12)));
    let p = project_features(&source, &[segment(0, 6, 0, false)], 6).unwrap();
    assert_eq!(p.features.len(), 1);
    assert!(p.features[0].needs_review);
    assert!(matches!(p.features[0].location, Location::Order(_)));
    assert_eq!(p.unmapped_features.len(), 1);
    assert_eq!(p.dropped_features, 1);
    assert!(p.unmapped_features[0].needs_review);
    let product = extract_region(&source, 0, 6).unwrap();
    let evidence: Vec<Feature> =
        serde_json::from_str(&product.metadata["annotation_mapping:unmapped_features"]).unwrap();
    assert_eq!(evidence[0].location, source.features[1].location);
    assert!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .is_err()
    );
}

#[test]
fn mapping_bounds_ambiguity_and_unsupported_coordinate_qualifiers_fail_explicitly() {
    let mut source = molecule("AAAACCCC", Topology::Linear);
    source.features.push(feature(Location::range(0, 4)));
    assert!(
        project_features(
            &source,
            &[segment(0, 4, 0, false), segment(3, 5, 4, false)],
            8
        )
        .is_err()
    );
    assert!(
        project_features(
            &source,
            &[segment(0, 4, 0, false), segment(4, 8, 3, false)],
            8
        )
        .is_err()
    );
    assert!(project_features(&source, &[segment(0, 4, usize::MAX, false)], 8).is_err());
    source.features[0]
        .qualifiers
        .push(("rpt_unit_range".into(), Some("1..3".into())));
    assert!(reverse_complement(&source).is_err());
    assert!(project_features(&source, &[segment(0, 8, 0, false)], 8).is_ok());
    let omitted = project_features(&source, &[segment(4, 8, 0, false)], 4).unwrap();
    assert!(omitted.features.is_empty());
    assert_eq!(omitted.dropped_features, 1);
    assert!(rotate_origin(&source, 1).is_err());
    assert!(extract_region(&source, 0, 0).is_err());
    assert!(extract_region(&source, 7, 10).is_err());
}

#[test]
fn derived_documents_preserve_original_headers_without_reusing_stale_reference_coordinates() {
    let mut source = molecule("AAAACCCCGGGGTTTT", Topology::Circular); // palindromic: sequence bytes alone cannot detect coordinate transform
    let original = "LOCUS       old 16 bp DNA circular\nREFERENCE   1  (bases 2 to 6)\nORIGIN\n        1 aaaaccccggggtttt\n//\n";
    source.original_genbank = Some(original.into());
    source.metadata.insert(
        "embl.original_record".into(),
        "old EMBL header with RP 2-6".into(),
    );
    source.metadata.insert(
        "source:canonical_url".into(),
        "https://www.ncbi.nlm.nih.gov/nuccore/SYN000001.1".into(),
    );
    source.features.push(feature(Location::range(2, 6)));
    let product = reverse_complement(&source).unwrap();
    assert_eq!(product.sequence, source.sequence);
    assert!(product.original_genbank.is_none());
    assert!(!product.metadata.contains_key("embl.original_record"));
    assert_eq!(
        product.metadata["annotation_mapping:original_genbank"],
        original
    );
    assert_eq!(
        product.metadata["annotation_mapping:original_embl"],
        "old EMBL header with RP 2-6"
    );
    assert_eq!(
        product.metadata["source:canonical_url"],
        source.metadata["source:canonical_url"]
    );
    let gb = collagene_core::formats::sequence::export_genbank(&product).unwrap();
    assert!(!gb.lines().any(|line| line.starts_with("REFERENCE")));
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "AACC"
    );
    assert_eq!(
        collagene_core::protocol::load_document(
            &collagene_core::protocol::save_document(&product).unwrap()
        )
        .unwrap(),
        product
    );
}

#[test]
fn changing_reference_topology_marks_origin_crossing_feature_and_does_not_claim_ligation() {
    let mut source = molecule("AAAACCCCGGGGTTTT", Topology::Circular);
    source.features.push(feature(Location::Join(vec![
        Location::range(12, 16),
        Location::range(0, 4),
    ])));
    let linear = change_reference_topology(&source, Topology::Linear).unwrap();
    assert_eq!(linear.sequence, source.sequence);
    assert!(linear.features[0].needs_review);
    assert!(
        linear.features[0]
            .location
            .extract(&linear.sequence)
            .is_err()
    );
    assert!(
        linear.metadata["annotation_mapping:reference_topology_note"]
            .contains("no physical junction")
    );
    assert_eq!(
        change_reference_topology(&source, Topology::Circular).unwrap(),
        source
    );
    let protein = Molecule {
        sequence: Sequence::new(Alphabet::Protein, "MKW").unwrap(),
        ..molecule("AAA", Topology::Linear)
    };
    assert!(change_reference_topology(&protein, Topology::Circular).is_err());
    let mut physical = molecule("AAACCC", Topology::Linear);
    physical
        .metadata
        .insert("duplex_model".into(), "source physical evidence".into());
    assert!(reverse_complement(&physical).is_err());
    assert!(change_reference_topology(&physical, Topology::Circular).is_err());
}

#[test]
fn reverse_partial_end_flags_and_new_internal_insertions_stay_reviewed() {
    let mut source = molecule("AACCGGTT", Topology::Linear);
    source.features.push(feature(Location::Range {
        start: 1,
        end: 5,
        before: true,
        after: false,
    }));
    let projection = project_features(&source, &[segment(0, 8, 0, true)], 8).unwrap();
    assert_eq!(
        projection.features[0].location,
        Location::Complement(Box::new(Location::Range {
            start: 3,
            end: 7,
            before: false,
            after: true
        }))
    );
    assert!(projection.features[0].needs_review);
    let mut exact = molecule("AAACCC", Topology::Linear);
    exact.features.push(feature(Location::range(0, 6)));
    // Product contains an inserted two-base region between previously adjacent source bases.
    let projected = project_features(
        &exact,
        &[segment(0, 3, 0, false), segment(3, 6, 5, false)],
        8,
    )
    .unwrap();
    assert!(projected.features[0].needs_review);
    assert!(
        !projected.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation")
    );
    assert!(
        projected.features[0]
            .location
            .extract(&Sequence::dna("AAATTCCC").unwrap())
            .is_err()
    );
    // A source intron gap was already an explicitly spliced join; merely retaining it is not a new interruption.
    let mut spliced = molecule("AAATTTCCC", Topology::Linear);
    spliced.features.push(feature(Location::Join(vec![
        Location::range(0, 3),
        Location::range(6, 9),
    ])));
    let projection = project_features(&spliced, &[segment(0, 9, 2, false)], 11).unwrap();
    assert!(!projection.features[0].needs_review);
    assert_eq!(
        projection.features[0]
            .location
            .extract(&Sequence::dna("GGAAATTTCCC").unwrap())
            .unwrap(),
        "AAACCC"
    );
}

fn simple_scf() -> Vec<u8> {
    let mut bytes = vec![0u8; 160];
    bytes[..4].copy_from_slice(b".scf");
    for (offset, value) in [(4, 2u32), (8, 128), (12, 2), (20, 2), (24, 136), (40, 1)] {
        bytes[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
    }
    bytes[36..40].copy_from_slice(b"2.00");
    bytes[128..136].copy_from_slice(&[10, 0, 0, 0, 0, 20, 0, 0]);
    bytes[136..148].copy_from_slice(&[0, 0, 0, 0, 40, 0, 0, 0, b'A', 0, 0, 0]);
    bytes[148..160].copy_from_slice(&[0, 0, 0, 1, 0, 40, 0, 0, b'C', 0, 0, 0]);
    bytes
}
#[test]
fn full_trace_reverse_updates_evidence_view_and_region_extract_does_not_fake_remapped_signals() {
    let original = simple_scf();
    let trace = collagene_core::review::reads::traces::import_trace(&original).unwrap();
    let mut source = trace.attached_molecule().unwrap();
    source.features.push(feature(Location::range(0, 1)));
    let reverse = reverse_complement(&source).unwrap();
    assert_eq!(reverse.sequence.as_str(), "GT");
    let view = collagene_core::review::reads::traces::trace_from_molecule(&reverse).unwrap();
    assert_eq!(
        view.orientation,
        collagene_core::review::reads::traces::TraceOrientation::ReverseComplement
    );
    assert_eq!(view.channels.g, vec![20, 0]);
    assert_eq!(view.channels.t, vec![0, 10]);
    assert_eq!(view.original_bytes, original);
    assert_eq!(
        reverse.features[0]
            .location
            .extract(&reverse.sequence)
            .unwrap(),
        "A"
    );
    let back = reverse_complement(&reverse).unwrap();
    assert_eq!(back.sequence, source.sequence);
    assert_eq!(back.features, source.features);
    assert_eq!(
        collagene_core::review::reads::traces::trace_from_molecule(&back)
            .unwrap()
            .orientation,
        collagene_core::review::reads::traces::TraceOrientation::AsStored
    );
    let cropped = extract_region(&source, 0, 1).unwrap();
    assert!(collagene_core::review::reads::traces::trace_from_molecule(&cropped).is_err());
    assert_eq!(
        collagene_core::review::reads::traces::original_trace_bytes(&cropped).unwrap(),
        original
    );
    assert!(change_reference_topology(&source, Topology::Circular).is_err());
}

#[test]
fn remapping_does_not_preserve_translation_assertions_across_new_exon_permutations() {
    let mut source = molecule("AAATTTCCC", Topology::Linear);
    source.features.push(feature(Location::Join(vec![
        Location::range(0, 3),
        Location::range(6, 9),
    ])));
    let swapped = project_features(
        &source,
        &[segment(0, 3, 6, false), segment(6, 9, 0, false)],
        9,
    )
    .unwrap();
    assert!(swapped.features[0].needs_review);
    assert!(
        swapped.features[0]
            .location
            .extract(&Sequence::dna("CCCTTTAAA").unwrap())
            .is_err()
    );
    assert!(
        !swapped.features[0]
            .qualifiers
            .iter()
            .any(|(k, _)| k == "translation")
    );
    let mixed = project_features(
        &source,
        &[segment(0, 3, 0, false), segment(6, 9, 6, true)],
        9,
    )
    .unwrap();
    assert!(mixed.features[0].needs_review);
    assert!(
        mixed.features[0]
            .location
            .extract(&Sequence::dna("AAATTTGGG").unwrap())
            .is_err()
    );
}

#[test]
fn moving_legacy_header_evidence_does_not_erase_explicit_single_stranded_state() {
    let mut legacy = molecule("AAGC", Topology::Linear);
    legacy.original_genbank = Some(
        "LOCUS       legacy 4 bp ss-DNA linear SYN 01-JAN-2026\nORIGIN\n        1 aagc\n//\n"
            .into(),
    );
    assert!(!legacy.metadata.contains_key("genbank.molecule_type"));
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &legacy
    ));
    let reverse = reverse_complement(&legacy).unwrap();
    assert_eq!(reverse.sequence.as_str(), "GCTT");
    assert_eq!(reverse.metadata["genbank.molecule_type"], "ss-DNA");
    assert!(reverse.original_genbank.is_none());
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &reverse
    ));
    assert!(collagene_core::chemistry::Duplex::from_molecule(&reverse).is_err());
}
