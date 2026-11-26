//! Literal protein bytes independent of the writer are the primary oracle.
use collagene_core::{
    snapgene,
    snapgene_export::{self, Policy},
    *,
};
fn packet(bytes: &mut Vec<u8>, tag: u8, data: &[u8]) {
    bytes.push(tag);
    bytes.extend_from_slice(&(data.len() as u32).to_be_bytes());
    bytes.extend_from_slice(data)
}
fn literal(flags: u8, residues: &[u8], extra: &[(u8, &[u8])]) -> Vec<u8> {
    let mut b = vec![
        9, 0, 0, 0, 14, b'S', b'n', b'a', b'p', b'G', b'e', b'n', b'e', 0, 2, 0, 15, 0, 15,
    ];
    let mut seq = vec![flags];
    seq.extend_from_slice(residues);
    packet(&mut b, 21, &seq);
    for (tag, payload) in extra {
        packet(&mut b, *tag, payload)
    }
    b
}
#[test]
fn independent_literal_protein_cookie_packet_and_aa_features_retain_exact_source_and_no_chemical_authority()
 {
    let notes=br#"<Notes><CustomMapLabel>Literal protein</CustomMapLabel><UseCustomMapLabel>1</UseCustomMapLabel><ConfirmedExperimentally>1</ConfirmedExperimentally></Notes>"#;
    let features=br#"<Features><Feature name="Domain" type="domain" directionality="1"><Segment range="2-5"/><Q name="note"><V text="first"/><V text="second"/></Q><Q name="experimental_claim"/></Feature></Features>"#;
    let bytes = literal(0, b"mKwVtfISLL*", &[(6, notes), (10, features)]);
    let report = snapgene::import_snapgene(&bytes).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.name, "Literal protein");
    assert_eq!(m.sequence.alphabet(), Alphabet::Protein);
    assert_eq!(m.sequence.as_str(), "MKWVTFISLL*");
    assert_eq!(m.topology, Topology::Linear);
    assert_eq!(m.features[0].location, Location::range(1, 5));
    assert_eq!(m.features[0].location.extract(&m.sequence).unwrap(), "KWVT");
    assert!(
        m.features[0]
            .qualifiers
            .iter()
            .any(|(k, v)| k == "experimental_claim" && v.is_none())
    );
    assert_eq!(
        m.features[0]
            .qualifiers
            .iter()
            .filter(|(k, _)| k == "note")
            .count(),
        2
    );
    assert_eq!(m.metadata["source:snapgene:protein_flags"], "0");
    assert!(!m.metadata.contains_key("source:snapgene:dna_flags"));
    assert_eq!(
        m.metadata["source:snapgene:note:ConfirmedExperimentally"],
        "1"
    );
    assert!(m.metadata["protein:authority"].contains("documentary"));
    assert!(!m.metadata.contains_key("verified"));
    assert!(report.warnings.iter().any(|s| s.contains("No PTM")));
    assert_eq!(snapgene::original_snapgene_bytes(m).unwrap(), bytes);
    let loaded = collagene_core::protocol::load_document(
        &collagene_core::protocol::save_document(m).unwrap(),
    )
    .unwrap();
    assert_eq!(snapgene::original_snapgene_bytes(&loaded).unwrap(), bytes);
    // Alphabet follows the cookie, not a heuristic that mistakes protein A/C/G/T.
    let acgt = snapgene::import_snapgene(&literal(0, b"ACGT", &[])).unwrap();
    assert_eq!(acgt.molecules[0].sequence.alphabet(), Alphabet::Protein);
}
#[test]
fn conflicting_compressed_duplicate_circular_or_invalid_protein_packets_are_refused() {
    for flags in [1, 2, 4, 8, 16, 128] {
        assert!(
            snapgene::import_snapgene(&literal(flags, b"MKWV", &[])).is_err(),
            "flags{flags}"
        )
    }
    for extra in [
        vec![(0, b"\x02ACGT".as_slice())],
        vec![(1, b"unsupported compressed bytes".as_slice())],
        vec![(21, b"\x00ACGT".as_slice())],
        vec![(32, b"\x00ACGU".as_slice())],
    ] {
        assert!(snapgene::import_snapgene(&literal(0, b"MKWV", &extra)).is_err())
    }
    let mut wrong = literal(0, b"MKWV", &[]);
    wrong[13..15].copy_from_slice(&1u16.to_be_bytes());
    assert!(snapgene::import_snapgene(&wrong).is_err());
    let mut wrong = literal(0, b"MKWV", &[]);
    wrong[13..15].copy_from_slice(&7u16.to_be_bytes());
    assert!(snapgene::import_snapgene(&wrong).is_err());
    for residues in [b"".as_slice(), b"MK WV", b"MK+V", b"MK\xffV"] {
        assert!(snapgene::import_snapgene(&literal(0, residues, &[])).is_err())
    }
    let mut no_sequence = literal(0, b"MKWV", &[]);
    no_sequence.truncate(19);
    assert!(snapgene::import_snapgene(&no_sequence).is_err());
    let mut truncated = literal(0, b"MKWV", &[]);
    truncated.pop();
    assert!(snapgene::import_snapgene(&truncated).is_err());
}
#[test]
fn protein_feature_coordinates_are_residues_and_reverse_wrap_or_chemical_interval_guesses_are_not_allowed()
 {
    for xml in [
        br#"<Features><Feature directionality="2"><Segment range="1-3"/></Feature></Features>"#
            .as_slice(),
        br#"<Features><Feature directionality="1"><Segment range="9-3"/></Feature></Features>"#,
        br#"<Features><Feature directionality="1"><Segment range="1-12"/></Feature></Features>"#,
    ] {
        assert!(snapgene::import_snapgene(&literal(0, b"MKWVTFISLL*", &[(10, xml)])).is_err())
    }
    let features=br#"<Features><Feature type="domain" directionality="1"><Segment range="2-3"/><Segment range="7-9"/></Feature><Feature type="DISULFID" directionality="1"><Segment range="2-9"/></Feature></Features>"#;
    let r = snapgene::import_snapgene(&literal(0, b"MCKVTFISCL*", &[(10, features)])).unwrap();
    let m = &r.molecules[0];
    assert_eq!(
        m.features[0].location,
        Location::Join(vec![Location::range(1, 3), Location::range(6, 9)])
    );
    assert!(matches!(m.features[1].location, Location::Unresolved(_)));
    assert!(m.features[1].needs_review);
    assert!(
        r.warnings
            .iter()
            .any(|s| s.contains("no contiguous peptide"))
    );
    assert!(
        collagene_core::protein_properties::analyze_protein(m, &Default::default())
            .unwrap()
            .mass
            .is_none()
    );
}
#[test]
fn extra_properties_and_opaque_protein_history_remain_source_only_and_block_free_chain_chemistry() {
    let bytes = literal(
        0,
        b"MCKVTFISCL",
        &[
            (8, b"opaque additional source properties"),
            (7, b"opaque protein source history"),
            (99, b"unmapped source bytes"),
        ],
    );
    let r = snapgene::import_snapgene(&bytes).unwrap();
    let m = &r.molecules[0];
    assert_eq!(m.metadata["source:snapgene:chemistry_unmapped"], "true");
    assert!(m.metadata.contains_key("protein:modifications"));
    assert!(
        collagene_core::protein_properties::analyze_protein(m, &Default::default())
            .unwrap()
            .mass
            .is_none()
    );
    assert!(snapgene_export::export_snapgene(m, Policy::AllowOmissions).is_err());
    assert_eq!(snapgene::original_snapgene_bytes(m).unwrap(), bytes);
    assert!(r.warnings.iter().any(|w| w.contains("cloning history")));
    assert!(r.warnings.iter().any(|w| w.contains("packet 99")));
}
#[test]
fn current_protein_writer_roundtrip_is_secondary_to_literals_and_originals_survive_later_edits() {
    let m = Molecule {
        name: "Current protein".into(),
        sequence: Sequence::new(Alphabet::Protein, "MKWVTFISLL*").unwrap(),
        topology: Topology::Linear,
        features: vec![Feature {
            kind: "domain".into(),
            location: Location::range(1, 5),
            qualifiers: vec![("label".into(), Some("Domain".into()))],
            needs_review: false,
        }],
        metadata: Default::default(),
        original_genbank: None,
    };
    let out = snapgene_export::export_snapgene(&m, Policy::Strict).unwrap();
    assert_eq!(out.extension, ".prot");
    let r = snapgene::import_snapgene(&out.bytes).unwrap();
    let mut imported = r.molecules[0].clone();
    assert_eq!(imported.sequence, m.sequence);
    assert_eq!(imported.features[0].location, m.features[0].location);
    assert_eq!(imported.name, m.name);
    imported.sequence = Sequence::new(Alphabet::Protein, "MKAA").unwrap();
    imported.features.clear();
    assert_eq!(
        snapgene::original_snapgene_bytes(&imported).unwrap(),
        out.bytes
    );
}
