//! Independently specified bytes, XML and biological extraction oracles.
//! No vendor application implementation or exported fixture is used as oracle.
use collagene_core::{
    snapgene,
    snapgene_export::{self as writer, LossKind, Policy},
};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
fn packets(bytes: &[u8]) -> Vec<(u8, &[u8])> {
    assert_eq!(bytes[0], 9);
    assert_eq!(&bytes[1..5], &[0, 0, 0, 14]);
    let mut at = 19;
    let mut out = vec![];
    while at < bytes.len() {
        let kind = bytes[at];
        let n = ((bytes[at + 1] as usize) << 24)
            | ((bytes[at + 2] as usize) << 16)
            | ((bytes[at + 3] as usize) << 8)
            | (bytes[at + 4] as usize);
        at += 5;
        out.push((kind, &bytes[at..at + n]));
        at += n;
    }
    assert_eq!(at, bytes.len());
    out
}
fn feature(kind: &str, location: Location, qualifiers: Vec<(&str, Option<&str>)>) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        qualifiers: qualifiers
            .into_iter()
            .map(|(k, v)| (k.into(), v.map(str::to_string)))
            .collect(),
        needs_review: false,
    }
}
#[test]
fn literal_cookie_uncompressed_current_dna_and_notes_have_independent_exact_bytes() {
    let m = Molecule::dna("Literal circle", "ACGTRY", Topology::Circular).unwrap();
    let out = writer::export_snapgene(&m, Policy::Strict).unwrap();
    assert_eq!(out.extension, ".dna");
    assert_eq!(out.source_revision, m.revision());
    assert!(out.losses.is_empty());
    assert_eq!(
        &out.bytes[..19],
        &[
            9, 0, 0, 0, 14, b'S', b'n', b'a', b'p', b'G', b'e', b'n', b'e', 0, 1, 0, 15, 0, 15
        ]
    );
    let parsed = packets(&out.bytes);
    assert_eq!(parsed.iter().map(|p| p.0).collect::<Vec<_>>(), vec![0, 6]);
    assert_eq!(parsed[0].1, b"\x03ACGTRY");
    assert_eq!(parsed[1].1,b"<Notes><CustomMapLabel>Literal circle</CustomMapLabel><UseCustomMapLabel>1</UseCustomMapLabel><ConfirmedExperimentally>0</ConfirmedExperimentally></Notes>");
    let imported = snapgene::import_snapgene(&out.bytes).unwrap();
    let read = &imported.molecules[0];
    assert_eq!(read.name, m.name);
    assert_eq!(read.sequence, m.sequence);
    assert_eq!(read.topology, m.topology);
    let mut linear = m;
    linear.topology = Topology::Linear;
    let out = writer::export_snapgene(&linear, Policy::Strict).unwrap();
    assert_eq!(packets(&out.bytes)[0].1, b"\x02ACGTRY");
}
#[test]
fn reverse_origin_join_and_individually_reversed_parts_preserve_literal_biological_order_and_qualifiers()
 {
    let mut m = Molecule::dna("Reverse source", "AGTCCGATGCTA", Topology::Circular).unwrap();
    m.features = vec![
        feature(
            "CDS",
            Location::Complement(Box::new(Location::Join(vec![
                Location::range(8, 12),
                Location::range(0, 3),
            ]))),
            vec![
                ("label", Some("Reporter")),
                ("note", Some("A&B\t\"literal\"\nsecond")),
                ("note", Some("duplicate value")),
                ("pseudo", None),
                ("ApEinfo_fwdcolor", Some("#Ab12EF")),
                ("ApEinfo_revcolor", Some("#Ab12EF")),
            ],
        ),
        feature(
            "misc_feature",
            Location::Join(vec![
                Location::Complement(Box::new(Location::range(8, 12))),
                Location::Complement(Box::new(Location::range(0, 3))),
            ]),
            vec![("label", Some("Independent reversed parts"))],
        ),
    ];
    let out = writer::export_snapgene(&m, Policy::Strict).unwrap();
    let packet = packets(&out.bytes).into_iter().find(|p| p.0 == 10).unwrap();
    let text = std::str::from_utf8(packet.1).unwrap();
    assert!(text.contains("&#9;"));
    assert!(text.contains("&#10;"));
    assert!(text.contains("&amp;"));
    let xml = roxmltree::Document::parse(text).unwrap();
    let features = xml
        .root_element()
        .children()
        .filter(|n| n.has_tag_name("Feature"))
        .collect::<Vec<_>>();
    assert_eq!(features[0].attribute("directionality"), Some("2"));
    let ranges = |n: roxmltree::Node<'_, '_>| {
        n.children()
            .filter(|s| s.has_tag_name("Segment"))
            .map(|s| s.attribute("range").unwrap().to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(ranges(features[0]), vec!["9-12", "1-3"]);
    assert_eq!(ranges(features[1]), vec!["1-3", "9-12"]);
    assert!(
        features[0]
            .children()
            .filter(|s| s.has_tag_name("Segment"))
            .all(|s| s.attribute("color") == Some("#Ab12EF"))
    );
    let imported = snapgene::import_snapgene(&out.bytes).unwrap();
    let read = &imported.molecules[0];
    assert_eq!(
        read.features[0].location.extract(&read.sequence).unwrap(),
        "ACTTAGC"
    );
    assert_eq!(
        read.features[1].location.extract(&read.sequence).unwrap(),
        "TAGCACT"
    );
    let notes = read.features[0]
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "note")
        .map(|(_, v)| v.as_deref())
        .collect::<Vec<_>>();
    assert_eq!(
        notes,
        vec![Some("A&B\t\"literal\"\nsecond"), Some("duplicate value")]
    );
    assert!(
        read.features[0]
            .qualifiers
            .iter()
            .any(|(k, v)| k == "pseudo" && v.is_none())
    );
}
#[test]
fn linear_protein_uses_type_two_packet_twenty_one_and_refuses_circular_protein_or_rna() {
    let mut m = Molecule::dna("Annotated protein", "ACGT", Topology::Linear).unwrap();
    m.sequence = Sequence::new(Alphabet::Protein, "MKWVTFISLL*").unwrap();
    m.features.push(feature(
        "domain",
        Location::range(1, 5),
        vec![("label", Some("Domain"))],
    ));
    let out = writer::export_snapgene(&m, Policy::Strict).unwrap();
    assert_eq!(out.extension, ".prot");
    assert_eq!(&out.bytes[13..15], &[0, 2]);
    let p = packets(&out.bytes);
    assert_eq!(p[0].0, 21);
    assert_eq!(p[0].1, b"\x00MKWVTFISLL*");
    assert_eq!(p.iter().map(|p| p.0).collect::<Vec<_>>(), vec![21, 6, 10]);
    let imported = snapgene::import_snapgene(&out.bytes).unwrap();
    assert_eq!(imported.molecules[0].sequence.alphabet(), Alphabet::Protein);
    assert_eq!(imported.molecules[0].sequence, m.sequence);
    assert_eq!(
        imported.molecules[0].features[0].location,
        m.features[0].location
    );
    m.topology = Topology::Circular;
    for policy in [Policy::Strict, Policy::AllowOmissions] {
        assert!(writer::export_snapgene(&m, policy).is_err())
    }
    m.topology = Topology::Linear;
    m.sequence = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
    m.features.clear();
    assert!(writer::export_snapgene(&m, Policy::AllowOmissions).is_err());
}
#[test]
fn explicit_losses_do_not_silently_normalize_colors_arrows_fuzzy_features_or_metadata() {
    let mut m = Molecule::dna("Explicit losses", "ACGTTGCA", Topology::Linear).unwrap();
    m.features.push(feature(
        "CDS",
        Location::range(0, 4),
        vec![
            ("label", Some("Readable label")),
            ("ApEinfo_fwdcolor", Some("#FF0000")),
            ("ApEinfo_revcolor", Some("#0000FF")),
            ("snapgene_directionality", Some("2")),
            ("snapgene_directionality", Some("1")),
            ("note", Some("<b>unmapped rich text</b>")),
        ],
    ));
    m.features[0].needs_review = true;
    m.features.push(feature(
        "misc_feature",
        Location::Range {
            start: 4,
            end: 8,
            before: true,
            after: false,
        },
        vec![],
    ));
    m.metadata
        .insert("definition".into(), "First literal definition".into());
    m.metadata.insert(
        "description".into(),
        "Not merged into first definition".into(),
    );
    m.metadata.insert(
        "private:display".into(),
        "no implicit native extension".into(),
    );
    m.original_genbank = Some("Original source text retained in native only".into());
    let before = m.clone();
    let review = writer::review(&m).unwrap();
    for kind in [
        LossKind::FeatureColor,
        LossKind::FeatureDisplay,
        LossKind::FeatureReview,
        LossKind::QualifierEncoding,
        LossKind::FeatureLocation,
        LossKind::Metadata,
        LossKind::OriginalEvidence,
    ] {
        assert!(
            review.losses.iter().any(|l| l.kind == kind),
            "missing typed loss: {kind:?}"
        )
    }
    assert!(writer::export_snapgene(&m, Policy::Strict).is_err());
    let out = writer::export_snapgene(&m, Policy::AllowOmissions).unwrap();
    assert_eq!(m, before);
    let xml = packets(&out.bytes).into_iter().find(|p| p.0 == 10).unwrap();
    let doc = roxmltree::Document::parse(std::str::from_utf8(xml.1).unwrap()).unwrap();
    let features = doc
        .root_element()
        .children()
        .filter(|n| n.has_tag_name("Feature"))
        .collect::<Vec<_>>();
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].attribute("directionality"), Some("1"));
    assert!(
        features[0]
            .children()
            .filter(|n| n.has_tag_name("Segment"))
            .all(|n| n.attribute("color").is_none())
    );
    assert!(std::str::from_utf8(xml.1).unwrap().contains("#FF0000"));
    assert!(std::str::from_utf8(xml.1).unwrap().contains("#0000FF"));
    assert!(
        !std::str::from_utf8(xml.1)
            .unwrap()
            .contains("unmapped rich text")
    );
    let encoded = serde_json::to_vec(&out.losses).unwrap();
    let decoded: Vec<writer::Loss> = serde_json::from_slice(&encoded).unwrap();
    assert_eq!(decoded, out.losses);
}
#[test]
fn physical_reference_only_and_modified_protein_state_cannot_be_erased_by_allowing_omissions() {
    for key in [
        "duplex_model",
        "tracked_modifications",
        "source:declared_termini",
        "source:snapgene:chemistry_unmapped",
        "source:recombination:reference_only",
        "topo:actions_context",
    ] {
        let mut m = Molecule::dna("Marked DNA", "ACGT", Topology::Linear).unwrap();
        m.metadata.insert(key.into(), "true".into());
        for p in [Policy::Strict, Policy::AllowOmissions] {
            assert!(writer::export_snapgene(&m, p).is_err(), "{key}")
        }
    }
    let mut m = Molecule::dna("Producer marker", "ACGT", Topology::Linear).unwrap();
    m.original_genbank = Some(
        "COMMENT     Collagene physical product v1: unverified sequence-only evidence\n".into(),
    );
    assert!(writer::export_snapgene(&m, Policy::AllowOmissions).is_err());
    m.original_genbank = None;
    m.metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    assert!(writer::export_snapgene(&m, Policy::AllowOmissions).is_err());
    m.metadata.clear();
    m.sequence = Sequence::new(Alphabet::Protein, "MCAK").unwrap();
    m.metadata.insert(
        "protein:modified".into(),
        "native unsupported processing state".into(),
    );
    assert!(writer::export_snapgene(&m, Policy::AllowOmissions).is_err());
    m.metadata.clear();
    m.features.push(feature(
        "DISULFID",
        Location::Order(vec![Location::range(1, 2), Location::range(3, 4)]),
        vec![],
    ));
    assert!(
        writer::export_snapgene(&m, Policy::AllowOmissions).is_err(),
        "an unsupported chemical annotation cannot disappear into a free chain"
    );
}
#[test]
fn edited_import_never_reemits_old_sequence_history_or_opaque_packets_and_retains_exact_original_source()
 {
    fn add(b: &mut Vec<u8>, tag: u8, p: &[u8]) {
        b.push(tag);
        b.extend_from_slice(&(p.len() as u32).to_be_bytes());
        b.extend_from_slice(p)
    }
    let mut original = vec![
        9, 0, 0, 0, 14, b'S', b'n', b'a', b'p', b'G', b'e', b'n', b'e', 0, 1, 0, 15, 0, 15,
    ];
    add(&mut original, 0, b"\x03ACGT");
    add(
        &mut original,
        7,
        b"opaque historical state must not be copied",
    );
    add(&mut original, 99, b"unknown cached source packet");
    let report = snapgene::import_snapgene(&original).unwrap();
    let mut m = report.molecules[0].clone();
    let old_revision = m.revision();
    m.name = "Edited source".into();
    m.sequence = Sequence::dna("TTACGTRYN").unwrap();
    let input = m.clone();
    assert!(writer::export_snapgene(&m, Policy::Strict).is_err());
    let out = writer::export_snapgene(&m, Policy::AllowOmissions).unwrap();
    assert_ne!(old_revision, out.source_revision);
    assert_eq!(m, input);
    assert_eq!(snapgene::original_snapgene_bytes(&m).unwrap(), original);
    let p = packets(&out.bytes);
    assert_eq!(p.iter().map(|p| p.0).collect::<Vec<_>>(), vec![0, 6]);
    assert_eq!(p[0].1, b"\x03TTACGTRYN");
    assert!(!out.bytes.windows(16).any(|s| s == b"opaque historica"));
    assert!(
        out.losses
            .iter()
            .any(|l| l.kind == LossKind::OriginalEvidence)
    );
}
#[test]
fn invalid_xml_oversized_native_models_and_unresolved_locations_never_produce_silent_partial_files()
{
    let mut m = Molecule::dna("Control\0name", "ACGT", Topology::Linear).unwrap();
    assert!(writer::export_snapgene(&m, Policy::AllowOmissions).is_err());
    m.name = "Bounds".into();
    m.features = vec![feature("misc_feature", Location::range(0, 1), vec![]); 4097];
    assert!(matches!(writer::review(&m), Err(Error::TooLarge)));
    m.features = vec![feature(
        "misc_feature",
        Location::Unresolved("remote:1..4".into()),
        vec![],
    )];
    assert!(writer::export_snapgene(&m, Policy::Strict).is_err());
    let out = writer::export_snapgene(&m, Policy::AllowOmissions).unwrap();
    assert_eq!(out.losses.len(), 1);
    assert!(
        String::from_utf8(
            packets(&out.bytes)
                .into_iter()
                .find(|p| p.0 == 10)
                .unwrap()
                .1
                .to_vec()
        )
        .unwrap()
        .contains("<Features></Features>")
    );
    m.features = vec![feature(
        "misc_feature",
        Location::range(0, 4),
        vec![("note", Some("bad\u{1}control"))],
    )];
    let out = writer::export_snapgene(&m, Policy::AllowOmissions).unwrap();
    assert!(
        out.losses
            .iter()
            .any(|l| l.kind == LossKind::QualifierEncoding)
    );
    assert!(
        !packets(&out.bytes)
            .into_iter()
            .find(|p| p.0 == 10)
            .unwrap()
            .1
            .contains(&1)
    );
}
