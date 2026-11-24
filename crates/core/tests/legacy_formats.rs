use collagene_core::{
    chemistry::Duplex,
    legacy_formats::{self as l},
    *,
};
fn xml(sequence: &str, feature: &str) -> String {
    format!(
        "<f><sequence name=\"test\" su=\"{sequence}\" sl=\"{}\" circular=\"true\" notes=\"literal &amp; notes\"/>{feature}</f>",
        sequence
            .bytes()
            .map(|base| model::complement(base, Alphabet::Dna) as char)
            .collect::<String>()
    )
}
fn annotation(from: i64, to: i64, direction: &str) -> String {
    format!(
        "<annotation name=\"feature\" from=\"{from}\" to=\"{to}\" colr=\"1\" colg=\"2\" colb=\"3\" orientation=\"{direction}\"/>"
    )
}
fn xdna(sequence: &[u8], kind: u8, circle: bool, comment: &[u8]) -> Vec<u8> {
    let mut bytes = vec![0; 112];
    bytes[1] = kind;
    bytes[2] = u8::from(circle);
    bytes[28..32].copy_from_slice(&(sequence.len() as u32).to_be_bytes());
    bytes[96..100].copy_from_slice(&(comment.len() as u32).to_be_bytes());
    bytes[111] = 255;
    bytes.extend_from_slice(sequence);
    bytes.extend_from_slice(comment);
    bytes
}
fn p(bytes: &mut Vec<u8>, value: &str) {
    bytes.push(value.len() as u8);
    bytes.extend_from_slice(value.as_bytes());
}
#[test]
fn original_serializer_schema_preserves_notes_names_features_and_primers() {
    let bytes = include_bytes!("../../../fixtures/legacy/serializer-shaped-project.collagene");
    let report = l::import_collagene_xml(bytes).unwrap();
    assert!(l::is_collagene_xml(bytes));
    assert_eq!(report.molecules.len(), 2);
    let molecule = &report.molecules[0];
    assert_eq!(molecule.name, "legacy_circle");
    assert_eq!(
        molecule.metadata["description"],
        "Retained & original notes"
    );
    assert_eq!(molecule.features.len(), 4);
    assert_eq!(
        molecule.features[0]
            .location
            .extract(&molecule.sequence)
            .unwrap(),
        "GTGC"
    );
    assert!(matches!(molecule.features[1].location, Location::Order(_)));
    assert!(
        molecule.features[1]
            .location
            .extract(&molecule.sequence)
            .is_err()
    );
    assert_eq!(molecule.features[2].location.spans(), vec![(0, 4, false)]);
    assert_eq!(molecule.features[3].location.spans(), vec![(4, 8, true)]);
    assert_eq!(
        molecule.features[3]
            .qualifiers
            .iter()
            .find(|(key, _)| key == "oligo_sequence")
            .unwrap()
            .1
            .as_deref(),
        Some("GCTA")
    );
    assert!(molecule.features.iter().all(|f| f.needs_review));
    assert_eq!(l::original_legacy_bytes(molecule).unwrap(), bytes);
    assert_eq!(report.molecules[0].name, report.molecules[1].name);
    assert_ne!(
        report.molecules[0].metadata["legacy:record_index"],
        report.molecules[1].metadata["legacy:record_index"]
    );
}
#[test]
fn aligned_gap_columns_are_projected_without_fabricating_duplex_chemistry() {
    let bytes = include_bytes!("../../../fixtures/legacy/serializer-shaped-project.collagene");
    let report = l::import_collagene_xml(bytes).unwrap();
    let molecule = &report.molecules[1];
    assert_eq!(molecule.sequence.as_str(), "ACGTTAC");
    assert_eq!(
        molecule.features[0]
            .location
            .extract(&molecule.sequence)
            .unwrap(),
        "GTTAC"
    );
    assert!(matches!(
        molecule.features[1].location,
        Location::Unresolved(_)
    ));
    assert!(Duplex::from_molecule(molecule).is_err());
    assert_eq!(molecule.metadata["legacy:upper_aligned"], "AC GT TAC");
    assert_eq!(molecule.metadata["legacy:lower_aligned"], "TGCCAAATG");
    let reopened = protocol::load_document(&protocol::save_document(molecule).unwrap()).unwrap();
    assert_eq!(reopened.revision(), molecule.revision());
    assert_eq!(l::original_legacy_bytes(&reopened).unwrap(), bytes);
}
#[test]
fn legacy_pair_mismatch_and_iupac_source_bug_never_grant_intact_duplex_authority() {
    for (upper, lower) in [("ACGT", "AAAA"), ("SWNN", "WSNN")] {
        let bytes = format!(
            "<f><sequence name=\"evidence\" su=\"{upper}\" sl=\"{lower}\" circular=\"false\"/></f>"
        );
        let molecule = l::import_collagene_xml(bytes.as_bytes())
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(molecule.metadata["genbank.molecule_type"], "ms-DNA");
        assert!(Duplex::from_molecule(&molecule).is_err());
    }
    let correct = xml("SWNN", "");
    let molecule = l::import_collagene_xml(correct.as_bytes())
        .unwrap()
        .molecules
        .remove(0);
    assert!(!molecule.metadata.contains_key("genbank.molecule_type")); // still no invented terminal phosphorylation
}
#[test]
fn coordinates_match_independent_circular_forward_slices_and_literal_reverse_complements() {
    let reference = "ACGTTAGC";
    let literal_complements = |text: &str| {
        text.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => unreachable!(),
            })
            .collect::<String>()
    };
    let mut cases = 0;
    for start in 0..8 {
        for length in 1..=8 {
            for reverse in [false, true] {
                let end = (start + length) % 8;
                let end = if end == 0 { 8 } else { end };
                if start == end {
                    continue;
                }
                let bytes = xml(
                    reference,
                    &annotation(
                        start as i64,
                        end as i64,
                        if reverse { "rev" } else { "fwd" },
                    ),
                );
                let molecule = l::import_collagene_xml(bytes.as_bytes())
                    .unwrap()
                    .molecules
                    .remove(0);
                let expected = (0..length)
                    .map(|offset| reference.as_bytes()[(start + offset) % 8] as char)
                    .collect::<String>();
                let expected = if reverse {
                    literal_complements(&expected)
                } else {
                    expected
                };
                assert_eq!(
                    molecule.features[0]
                        .location
                        .extract(&molecule.sequence)
                        .unwrap(),
                    expected,
                    "{start} {length} {reverse}"
                );
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 114);
}
#[test]
fn empty_unoriented_overflow_and_unknown_schema_are_not_silently_repaired() {
    let empty = xml("ACGT", &annotation(2, 2, "fwd"));
    assert!(matches!(
        l::import_collagene_xml(empty.as_bytes()).unwrap().molecules[0].features[0].location,
        Location::Unresolved(_)
    ));
    let excessive = xml("ACGT", &annotation(i64::MIN, 2, "fwd"));
    assert!(matches!(
        l::import_collagene_xml(excessive.as_bytes())
            .unwrap()
            .molecules[0]
            .features[0]
            .location,
        Location::Unresolved(_)
    ));
    for input in [
        "<madgene/>",
        "<madgene><future/></madgene>",
        "<f><sequence name=\"x\" su=\"AC\" sl=\"T\" circular=\"false\"/></f>",
        "<!DOCTYPE f [<!ENTITY x 'AAAA'>]><f>&x;</f>",
        "<f xmlns=\"future\"><sequence/></f>",
    ] {
        assert!(
            l::import_collagene_xml(input.as_bytes()).is_err(),
            "{input}"
        );
    }
    let unknown = xml("ACGT", "<future><operation/></future>");
    let report = l::import_collagene_xml(unknown.as_bytes()).unwrap();
    assert!(report.warnings.iter().any(|w| w.contains("future")));
    assert_eq!(
        l::original_legacy_bytes(&report.molecules[0]).unwrap(),
        unknown.as_bytes()
    );
}
#[test]
fn xml_single_writer_f_and_old_detector_sequence_wrapper_are_explicit() {
    let input = xml("ACGT", "");
    let molecule = l::import_collagene_xml(input.as_bytes())
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(molecule.metadata["description"], "literal & notes");
    let wrapped = input
        .replacen("<f>", "<sequence>", 1)
        .replacen("</f>", "</sequence>", 1);
    assert!(l::is_collagene_xml(wrapped.as_bytes()));
    assert_eq!(
        l::import_collagene_xml(wrapped.as_bytes())
            .unwrap()
            .molecules[0]
            .sequence
            .as_str(),
        "ACGT"
    );
    assert!(!l::is_collagene_xml(b"<sequence su=\"ACGT\"/>"));
}
#[test]
fn independent_xdna_base_headers_support_dna_rna_protein_and_origin_evidence() {
    for (kind, alphabet, bases) in [
        (1, Alphabet::Dna, "ACGT"),
        (2, Alphabet::Dna, "ARYN"),
        (3, Alphabet::Rna, "ACGU"),
        (4, Alphabet::Protein, "MKW*"),
    ] {
        let mut bytes = xdna(bases.as_bytes(), kind, false, b"source comment");
        bytes[32..36].copy_from_slice(&2u32.to_be_bytes());
        let report = l::import_xdna(&bytes, "file name").unwrap();
        let molecule = &report.molecules[0];
        assert_eq!(molecule.sequence.alphabet(), alphabet);
        assert_eq!(molecule.sequence.as_str(), bases);
        assert_eq!(molecule.metadata["description"], "source comment");
        assert_eq!(molecule.metadata["legacy:xdna_origin_offset"], "2");
        assert_eq!(l::original_legacy_bytes(molecule).unwrap(), bytes);
    }
}
#[test]
fn xdna_extended_features_use_inclusive_reversed_and_wrapped_coordinates() {
    let mut bytes = xdna(b"ACGTTAGC", 1, true, b"title");
    bytes.push(0);
    p(&mut bytes, "0");
    p(&mut bytes, "0");
    bytes.push(3);
    for (name, start, end, strand, expected) in [
        ("forward", "2", "5", 1, "CGTT"),
        ("reverse", "8", "5", 0, "GCTA"),
        ("origin", "7", "2", 1, "GCAC"),
    ] {
        p(&mut bytes, name);
        p(&mut bytes, "literal note\rgene=\"geneA\"");
        p(&mut bytes, "CDS");
        p(&mut bytes, start);
        p(&mut bytes, end);
        bytes.extend_from_slice(&[strand, 1, 0, 1]);
        p(&mut bytes, "20,40,220");
        let _ = expected;
    }
    let molecule = l::import_xdna(&bytes, "fixture")
        .unwrap()
        .molecules
        .remove(0);
    for (feature, expected) in molecule.features.iter().zip(["CGTT", "GCTA", "GCAC"]) {
        assert_eq!(
            feature.location.extract(&molecule.sequence).unwrap(),
            expected
        );
        assert!(feature.needs_review);
        assert!(
            feature
                .qualifiers
                .iter()
                .any(|(key, value)| key == "color" && value.as_deref() == Some("#1428dc"))
        );
        assert!(!feature.qualifiers.iter().any(|(key, _)| key == "gene"));
    }
}
#[test]
fn xdna_overhangs_are_retained_without_silent_blunting_or_phosphorylation() {
    let mut bytes = xdna(b"ACGT", 1, false, b"");
    bytes.push(0);
    p(&mut bytes, "2");
    bytes.extend_from_slice(b"AA");
    p(&mut bytes, "-1");
    bytes.extend_from_slice(b"C");
    bytes.push(0);
    let molecule = l::import_xdna(&bytes, "overhangs")
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(molecule.metadata["legacy:xdna_right_length"], "2");
    assert_eq!(molecule.metadata["legacy:xdna_left_length"], "-1");
    assert!(Duplex::from_molecule(&molecule).is_err());
    assert_eq!(molecule.sequence.as_str(), "ACGT");
}
#[test]
fn xdna_truncation_unknown_flags_lengths_encoding_and_trailing_data_fail_atomically() {
    let valid = xdna(b"ACGT", 1, false, b"note");
    for length in 0..valid.len() {
        assert!(l::import_xdna(&valid[..length], "bad").is_err(), "{length}");
    }
    for (offset, value) in [(0, 1), (1, 9), (2, 3), (111, 0)] {
        let mut bad = valid.clone();
        bad[offset] = value;
        assert!(l::import_xdna(&bad, "bad").is_err());
    }
    let mut bad = valid.clone();
    bad[28..32].copy_from_slice(&u32::MAX.to_be_bytes());
    assert!(l::import_xdna(&bad, "bad").is_err());
    let mut bad = valid.clone();
    bad[112] = 255;
    assert!(l::import_xdna(&bad, "bad").is_err());
    let mut bad = valid.clone();
    bad.extend_from_slice(&[0, 1, b'0', 1, b'0', 0, 99]);
    assert!(l::import_xdna(&bad, "bad").is_err());
}
#[test]
fn original_legacy_evidence_rejects_hash_tampering() {
    let bytes = xdna(b"ACGT", 1, false, b"");
    let mut molecule = l::import_xdna(&bytes, "file").unwrap().molecules.remove(0);
    molecule
        .metadata
        .insert(l::ORIGINAL_BASE64_KEY.into(), "QUJDRA==".into());
    assert!(l::original_legacy_bytes(&molecule).is_err());
}
#[test]
fn malicious_primer_target_arithmetic_and_nonzero_xdna_origin_do_not_guess_coordinates() {
    for target in [i64::MAX, i64::MIN] {
        let body = format!(
            "<primer name=\"overflow\" sequence=\"ACGT\" target=\"{target}\" orientation=\"{}\"/>",
            if target > 0 { "rev" } else { "fwd" }
        );
        assert!(l::import_collagene_xml(xml("ACGT", &body).as_bytes()).is_err());
    }
    let mut bytes = xdna(b"ACGT", 1, false, b"");
    bytes[32..36].copy_from_slice(&1u32.to_be_bytes());
    bytes.push(0);
    p(&mut bytes, "0");
    p(&mut bytes, "0");
    bytes.push(1);
    for value in ["feature", "note", "misc_feature", "1", "2"] {
        p(&mut bytes, value);
    }
    bytes.extend_from_slice(&[1, 1, 0, 1]);
    p(&mut bytes, "0,0,0");
    let molecule = l::import_xdna(&bytes, "offset")
        .unwrap()
        .molecules
        .remove(0);
    assert!(matches!(
        molecule.features[0].location,
        Location::Unresolved(_)
    ));
    assert_eq!(molecule.sequence.as_str(), "ACGT");
    assert!(molecule.features[0].qualifiers.iter().any(|(key,value)|key=="legacy:xdna_coordinates"&&value.as_deref()==Some("1,2")));
}

#[test]
fn lower_row_is_aligned_complement_not_ordered_reverse_complement() {
    for (lower, unmapped) in [("TGCAATCG", false), ("GCTAACGT", true)] {
        let source = format!(
            "<f><sequence name=\"orientation\" su=\"ACGTTAGC\" sl=\"{lower}\" circular=\"false\"/></f>"
        );
        let molecule = l::import_collagene_xml(source.as_bytes())
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(molecule.metadata["legacy:lower_aligned"], lower);
        assert_eq!(
            molecule.metadata.get("genbank.molecule_type").is_some(),
            unmapped
        );
        if !unmapped {
            let duplex = Duplex::from_molecule(&molecule).unwrap();
            assert_eq!(duplex.bottom.as_str(), "GCTAACGT");
        } else {
            assert!(Duplex::from_molecule(&molecule).is_err());
        }
    }
    let deeply_nested = format!(
        "<f><sequence name=\"x\" su=\"AC\" sl=\"TG\" circular=\"false\"/>{}{} </f>",
        "<extra>".repeat(65),
        "</extra>".repeat(65)
    );
    assert!(l::import_collagene_xml(deeply_nested.as_bytes()).is_err());
}

#[test]
fn uninterpreted_material_fields_cannot_erase_unknown_physical_states() {
    for document in [
        "<f><sequence name=\"unknown\" su=\"ACGT\" sl=\"TGCA\" circular=\"false\" nicked=\"true\"/></f>",
        "<f><sequence name=\"unknown\" su=\"ACGT\" sl=\"TGCA\" circular=\"false\"/><nick strand=\"top\" at=\"2\"/></f>",
    ] {
        let report = l::import_collagene_xml(document.as_bytes()).unwrap();
        let source = &report.molecules[0];
        assert_eq!(source.metadata["source:legacy:chemistry_unmapped"], "true");
        assert!(Duplex::from_molecule(source).is_err());
        assert_eq!(
            l::original_legacy_bytes(source).unwrap(),
            document.as_bytes()
        );
    }
}
