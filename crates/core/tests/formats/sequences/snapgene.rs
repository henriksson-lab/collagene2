//! Independently authored packets and XML, not outputs of another cloning program.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::snapgene::{ORIGINAL_BASE64_KEY, import_snapgene, original_snapgene_bytes};
use collagene_model::{Location, Topology};

fn packet(kind: u8, data: &[u8]) -> Vec<u8> {
    let mut result = vec![kind];
    result.extend_from_slice(&(data.len() as u32).to_be_bytes());
    result.extend_from_slice(data);
    result
}
fn cookie() -> Vec<u8> {
    let mut data = b"SnapGene".to_vec();
    data.extend_from_slice(&1u16.to_be_bytes());
    data.extend_from_slice(&15u16.to_be_bytes());
    data.extend_from_slice(&15u16.to_be_bytes());
    packet(9, &data)
}
fn file(flags: u8, bases: &str, packets: &[(u8, &[u8])]) -> Vec<u8> {
    let mut bytes = cookie();
    let mut data = vec![flags];
    data.extend_from_slice(bases.as_bytes());
    bytes.extend(packet(0, &data));
    for (kind, data) in packets {
        bytes.extend(packet(*kind, data));
    }
    bytes
}

#[test]
fn dna_topology_and_source_evidence_are_exact_and_recoverable_after_edits() {
    for (flags, topology) in [(2, Topology::Linear), (3, Topology::Circular)] {
        let bytes = file(flags, "aCgTRYN", &[]);
        let report = import_snapgene(&bytes).unwrap();
        assert_eq!(report.molecules.len(), 1);
        let molecule = &report.molecules[0];
        assert_eq!(molecule.sequence.as_str(), "ACGTRYN");
        assert_eq!(molecule.topology, topology);
        assert!(
            report
                .warnings
                .iter()
                .any(|w| w.contains("case was normalized"))
        );
        assert_eq!(original_snapgene_bytes(molecule).unwrap(), bytes);
        let saved = collagene_core::protocol::save_document(molecule).unwrap();
        let restored = collagene_core::protocol::load_document(&saved).unwrap();
        assert_eq!(original_snapgene_bytes(&restored).unwrap(), bytes);
        let changed = molecule.replace(0, 1, "T").unwrap();
        assert_eq!(original_snapgene_bytes(&changed).unwrap(), bytes);
        let mut corrupt = molecule.clone();
        corrupt
            .metadata
            .insert(ORIGINAL_BASE64_KEY.into(), STANDARD.encode(b"wrong bytes"));
        assert!(original_snapgene_bytes(&corrupt).is_err());
    }
}

#[test]
fn origin_spanning_reverse_compound_features_and_duplicate_qualifiers() {
    let xml=br##"<?xml version="1.0"?><Features><Feature name="reverse &amp; joined" type="CDS" directionality="2"><Segment range="13-16" color="#00aa77"/><Segment range="1-4" color="#00aa77"/><Q name="note"><V text="first"/><V text="second"/></Q><Q name="codon_start"><V int="1"/></Q><Q name="pseudo"/><Q name="product"><V predef="illustrative protein"/></Q></Feature><Feature name="wrap" type="misc_feature" directionality="1"><Segment range="14-3"/></Feature><Feature name="one base" directionality="1"><Segment range="5-5"/></Feature></Features>"##;
    let report = import_snapgene(&file(3, "AAAACCCCGGGGTTTT", &[(10, xml)])).unwrap();
    let features = &report.molecules[0].features;
    assert_eq!(features.len(), 3);
    assert!(matches!(features[0].location, Location::Complement(_)));
    assert_eq!(
        features[0]
            .location
            .extract(&report.molecules[0].sequence)
            .unwrap(),
        "TTTTAAAA"
    );
    assert_eq!(
        features[1]
            .location
            .extract(&report.molecules[0].sequence)
            .unwrap(),
        "TTTAAA"
    );
    assert_eq!(features[2].location, Location::range(4, 5));
    assert_eq!(features[0].label(), "reverse & joined");
    let notes: Vec<_> = features[0]
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "note")
        .collect();
    assert_eq!(notes.len(), 2);
    assert_eq!(notes[0].1.as_deref(), Some("first"));
    assert_eq!(notes[1].1.as_deref(), Some("second"));
    assert!(features[0].qualifiers.contains(&("pseudo".into(), None)));
    assert!(
        features[0]
            .qualifiers
            .contains(&("ApEinfo_fwdcolor".into(), Some("#00aa77".into())))
    );
    assert!(!features[0].needs_review);
}

#[test]
fn notes_are_untrusted_source_assertions_and_do_not_grant_verified_status() {
    let notes=br#"<Notes><AccessionNumber>ACC123</AccessionNumber><UseCustomMapLabel>1</UseCustomMapLabel><CustomMapLabel>my plasmid</CustomMapLabel><Description>&lt;script&gt;data only&lt;/script&gt;</Description><ConfirmedExperimentally>1</ConfirmedExperimentally><Custom>unmapped</Custom></Notes>"#;
    let report = import_snapgene(&file(
        2,
        "ACGT",
        &[(6, notes), (0xfe, b"unmapped attachment")],
    ))
    .unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.name, "my plasmid");
    assert_eq!(
        molecule.metadata["source:snapgene:note:Description"],
        "<script>data only</script>"
    );
    assert_eq!(
        molecule.metadata["source:snapgene:note:ConfirmedExperimentally"],
        "1"
    );
    assert!(!molecule.metadata.contains_key("verified"));
    assert!(report.warnings.iter().any(|w| w.contains("packet 254")));
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("notes are retained"))
    );
}

#[test]
fn gaps_and_bidirectional_features_are_not_claimed_as_exact_biological_locations() {
    let xml=br#"<Features><Feature name="gap" directionality="1"><Segment range="1-3"/><Segment range="4-6" type="gap"/></Feature><Feature name="both" directionality="3"><Segment range="2-5"/></Feature></Features>"#;
    let report = import_snapgene(&file(2, "ACGTACGT", &[(10, xml)])).unwrap();
    let molecule = &report.molecules[0];
    assert!(matches!(
        molecule.features[0].location,
        Location::Unresolved(_)
    ));
    assert!(
        molecule.features[0]
            .location
            .extract(&molecule.sequence)
            .is_err()
    );
    assert!(molecule.features.iter().all(|f| f.needs_review));
    assert!(report.warnings.iter().any(|w| w.contains("unresolved")));
    assert!(report.warnings.iter().any(|w| w.contains("bidirectional")));
}

#[test]
fn primers_history_and_chemical_flags_are_retained_but_not_promoted_to_chemistry() {
    let bytes = file(
        3 | 4,
        "ACGTACGT",
        &[
            (5, b"<Primers/>"),
            (7, b"opaque history"),
            (8, b"<AdditionalSequenceProperties/>"),
            (16, b"trace container"),
        ],
    );
    let report = import_snapgene(&bytes).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(
        molecule.metadata["source:snapgene:chemistry_unmapped"],
        "true"
    );
    for message in ["primers", "history", "end/overhang", "traces"] {
        assert!(report.warnings.iter().any(|w| w.contains(message)));
    }
    assert_eq!(original_snapgene_bytes(molecule).unwrap(), bytes);
}

#[test]
fn truncated_packets_duplicates_unknown_sequence_encodings_and_invalid_bases_fail() {
    let good = file(2, "ACGT", &[]);
    for stop in 0..good.len() {
        assert!(
            import_snapgene(&good[..stop]).is_err(),
            "accepted prefix {stop}"
        );
    }
    let mut bad = good.clone();
    bad[20..24].copy_from_slice(&u32::MAX.to_be_bytes());
    assert!(import_snapgene(&bad).is_err());
    let mut duplicate = good.clone();
    duplicate.extend(packet(0, b"\x02ACGT"));
    assert!(import_snapgene(&duplicate).is_err());
    assert!(import_snapgene(&file(0, "ACGT", &[])).is_err());
    assert!(import_snapgene(&file(2 | 0x80, "ACGT", &[])).is_err());
    assert!(import_snapgene(&file(2, "AC GT", &[])).is_err());
    assert!(import_snapgene(&file(2, "AC!T", &[])).is_err());
    assert!(import_snapgene(&file(2, "ACGT", &[(1, b"compressed")])).is_err());
    let mut protein = good.clone();
    protein[13..15].copy_from_slice(&2u16.to_be_bytes());
    assert!(import_snapgene(&protein).is_err());
    let mut future = good;
    future[15..17].copy_from_slice(&16u16.to_be_bytes());
    assert!(import_snapgene(&future).is_err());
}

#[test]
fn xml_entities_nesting_and_invalid_coordinate_semantics_are_rejected() {
    let dtd =
        br#"<!DOCTYPE Features [<!ENTITY x SYSTEM "file:///etc/passwd">]><Features>&x;</Features>"#;
    assert!(import_snapgene(&file(2, "ACGT", &[(10, dtd)])).is_err());
    let deep = format!(
        "<Features>{}{}<Feature><Segment range=\"1-2\"/></Feature>{}</Features>",
        "<A>".repeat(40),
        "",
        "</A>".repeat(40)
    );
    assert!(import_snapgene(&file(2, "ACGT", &[(10, deep.as_bytes())])).is_err());
    for xml in [
        "<Features><Feature><Segment range=\"0-3\"/></Feature></Features>",
        "<Features><Feature><Segment range=\"1-5\"/></Feature></Features>",
        "<Features><Feature><Segment range=\"4-2\"/></Feature></Features>",
        "<Features><Feature><Segment range=\"1-2\"/><Q name=\"x\"><V text=\"a\" int=\"3\"/></Q></Feature></Features>",
    ] {
        assert!(import_snapgene(&file(2, "ACGT", &[(10, xml.as_bytes())])).is_err());
    }
    assert!(
        import_snapgene(&file(
            2,
            "ACGT",
            &[(6, b"<WrongRoot/>"), (10, b"<Features/>")]
        ))
        .is_err()
    );
}

#[test]
fn packets_need_not_follow_dna_and_fidelity_warning_truncation_is_explicit() {
    let mut bytes = cookie();
    bytes.extend(packet(10,b"<Features><Feature name=\"first\" directionality=\"1\"><Segment range=\"1-2\"/></Feature></Features>"));
    bytes.extend(packet(
        6,
        b"<Notes><AccessionNumber>ordered independently</AccessionNumber></Notes>",
    ));
    bytes.extend(packet(0, b"\x02ACGT"));
    for kind in 40..168 {
        bytes.extend(packet(kind, b"opaque"));
    }
    let report = import_snapgene(&bytes).unwrap();
    assert_eq!(report.molecules[0].name, "ordered independently");
    assert_eq!(
        report.molecules[0].features[0].location,
        Location::range(0, 2)
    );
    assert!(report.warnings.len() <= 100);
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("warnings were truncated"))
    );
    let counts: serde_json::Value =
        serde_json::from_str(&report.molecules[0].metadata["source:snapgene:packet_counts"])
            .unwrap();
    assert_eq!(counts["167"], 1);
    assert_eq!(
        original_snapgene_bytes(&report.molecules[0]).unwrap(),
        bytes
    );
}
