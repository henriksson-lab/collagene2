//! Literal byte/profile assertions precede reader roundtrips. Our reader alone
//! cannot establish compatibility with proprietary producer applications.
use collagene_core::formats::original::{
    self as legacy_formats,
    xdna_export::{self as x, LossKind, Policy},
};
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn dna() -> Molecule {
    Molecule::dna("Untruncated native name", "ACGTTAGC", Topology::Circular).unwrap()
}
fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("literal label".into())),
            ("note".into(), Some("literal note".into())),
            ("color".into(), Some("#1428dc".into())),
        ],
        needs_review: true,
    }
}
fn p(bytes: &mut Vec<u8>, s: &str) {
    bytes.push(u8::try_from(s.len()).unwrap());
    bytes.extend_from_slice(s.as_bytes())
}
fn header(kind: u8, circular: bool, sequence: &str, comment: &str) -> Vec<u8> {
    let mut bytes = vec![0; 112];
    bytes[1] = kind;
    bytes[2] = u8::from(circular);
    bytes[28..32].copy_from_slice(&u32::try_from(sequence.len()).unwrap().to_be_bytes());
    bytes[96..100].copy_from_slice(&u32::try_from(comment.len()).unwrap().to_be_bytes());
    bytes[111] = 255;
    bytes.extend_from_slice(sequence.as_bytes());
    bytes.extend_from_slice(comment.as_bytes());
    bytes
}
#[test]
fn plain_literal_header_and_zero_overhang_extension_match_the_selected_binary_contract() {
    let mut m = Molecule::dna("Literal current name", "ACGT", Topology::Linear).unwrap();
    m.metadata
        .insert("description".into(), "ASCII comment".into());
    let report = x::export_xdna(&m, Policy::AllowOmissions).unwrap();
    let mut literal = header(1, false, "ACGT", "ASCII comment");
    literal.extend_from_slice(&[0, 1, b'0', 1, b'0', 0]);
    assert_eq!(report.bytes, literal);
    assert_eq!(report.extension, ".xdna");
    assert_eq!(report.suggested_name, m.name);
    assert_eq!(report.source_revision, m.revision());
    assert_eq!(report.losses.len(), 1);
    assert_eq!(report.losses[0].kind, LossKind::Name);
    assert!(
        x::export_xdna(&m, Policy::Strict).is_err(),
        "Names are filename-context only, never silently claimed embedded"
    );
    let imported = legacy_formats::import_xdna(&literal, &m.name)
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(imported.sequence, m.sequence);
    assert_eq!(imported.topology, Topology::Linear);
    assert_eq!(imported.metadata["description"], "ASCII comment");
}
#[test]
fn exact_forward_reverse_and_origin_wrap_have_literal_inclusive_v0_coordinates() {
    let mut m = dna();
    m.features = vec![
        feature(Location::range(1, 5)),
        feature(Location::Complement(Box::new(Location::range(4, 8)))),
        feature(Location::Join(vec![
            Location::range(6, 8),
            Location::range(0, 2),
        ])),
    ];
    let report = x::export_xdna(&m, Policy::AllowOmissions).unwrap();
    let mut literal = header(1, true, "ACGTTAGC", "");
    literal.extend_from_slice(&[0, 1, b'0', 1, b'0', 3]);
    for (a, b, strand) in [("2", "5", 1), ("8", "5", 0), ("7", "2", 1)] {
        for field in ["literal label", "literal note", "CDS", a, b] {
            p(&mut literal, field)
        }
        literal.extend_from_slice(&[strand, 1, 0, 1]);
        p(&mut literal, "20,40,220");
    }
    assert_eq!(report.bytes, literal);
    let read = legacy_formats::import_xdna(&report.bytes, &m.name)
        .unwrap()
        .molecules
        .remove(0);
    for (f, expected) in read.features.iter().zip(["CGTT", "GCTA", "GCAC"]) {
        assert_eq!(f.location.extract(&read.sequence).unwrap(), expected);
        assert!(f.needs_review)
    }
    assert_eq!(
        read.features[0]
            .qualifiers
            .iter()
            .find(|q| q.0 == "color")
            .unwrap()
            .1
            .as_deref(),
        Some("#1428dc")
    );
    // A complement of the same wrap is biologically distinct and must traverse
    // both parts in reverse order, preserving exactly the literal reverse read.
    m.features = vec![feature(Location::Complement(Box::new(Location::Join(
        vec![Location::range(6, 8), Location::range(0, 2)],
    ))))];
    let bytes = x::export_xdna(&m, Policy::AllowOmissions).unwrap().bytes;
    let read = legacy_formats::import_xdna(&bytes, &m.name)
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(
        read.features[0].location.extract(&read.sequence).unwrap(),
        "GTGC"
    );
}
#[test]
fn protein_type4_is_linear_and_rna_or_circular_protein_are_explicitly_refused() {
    let mut m = dna();
    m.topology = Topology::Linear;
    m.sequence = Sequence::new(Alphabet::Protein, "MKWVTF").unwrap();
    m.features = vec![feature(Location::range(0, 3))];
    let report = x::export_xdna(&m, Policy::AllowOmissions).unwrap();
    assert_eq!(report.bytes[1], 4);
    assert_eq!(report.bytes[2], 0);
    let read = legacy_formats::import_xdna(&report.bytes, "Protein reference")
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(read.sequence.alphabet(), Alphabet::Protein);
    assert_eq!(read.sequence.as_str(), "MKWVTF");
    assert_eq!(read.features[0].location.spans(), vec![(0, 3, false)]);
    m.topology = Topology::Circular;
    assert!(x::review(&m).is_err());
    m.topology = Topology::Linear;
    m.sequence = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
    m.features.clear();
    assert!(x::review(&m).is_err());
}
#[test]
fn losses_are_explicit_without_truncation_join_flattening_or_color_spelling_changes() {
    let mut m = dna();
    let mut f = feature(Location::range(0, 3));
    f.needs_review = false;
    f.qualifiers = vec![
        ("label".into(), Some("x".repeat(256))),
        ("note".into(), Some("Unicode α".into())),
        ("color".into(), Some("#ABCDEF".into())),
        ("gene".into(), Some("typed gene".into())),
    ];
    m.features = vec![
        f,
        feature(Location::Join(vec![
            Location::range(0, 1),
            Location::range(3, 5),
        ])),
    ];
    m.metadata
        .insert("unmapped".into(), "preserve natively".into());
    let review = x::review(&m).unwrap();
    for kind in [
        LossKind::Name,
        LossKind::FeatureReview,
        LossKind::FeatureText,
        LossKind::FeatureColor,
        LossKind::Qualifier,
        LossKind::FeatureLocation,
        LossKind::Metadata,
    ] {
        assert!(
            review.losses.iter().any(|loss| loss.kind == kind),
            "{kind:?}"
        )
    }
    let report = x::export_xdna(&m, Policy::AllowOmissions).unwrap();
    let read = legacy_formats::import_xdna(&report.bytes, "Visible omissions")
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(read.features.len(), 1);
    assert_eq!(
        read.features[0]
            .qualifiers
            .iter()
            .find(|q| q.0 == "label")
            .unwrap()
            .1
            .as_deref(),
        Some("")
    );
    assert!(
        !String::from_utf8_lossy(&report.bytes).contains("x".repeat(255).as_str()),
        "No prefix truncation"
    );
    assert!(x::export_xdna(&m, Policy::Strict).is_err());
}
#[test]
fn physical_reference_only_and_modified_protein_evidence_cannot_be_omitted_into_clean_references() {
    for key in [
        "duplex_model",
        "tracked_modifications",
        "source:declared_termini",
        "source:recombination:reference_only",
        "topo:physical_graph_required",
        "legacy:xdna_right_length",
    ] {
        let mut m = dna();
        m.metadata.insert(key.into(), "true".into());
        assert!(x::export_xdna(&m, Policy::AllowOmissions).is_err(), "{key}");
    }
    let mut m = dna();
    m.original_genbank = Some(
        "COMMENT     Collagene physical product v1: unverified sequence-only evidence\n".into(),
    );
    assert!(x::review(&m).is_err());
    let mut m = dna();
    m.topology = Topology::Linear;
    m.sequence = Sequence::new(Alphabet::Protein, "MKWVTF").unwrap();
    let mut f = feature(Location::range(0, 1));
    f.kind = "MOD_RES".into();
    m.features.push(f);
    assert!(x::export_xdna(&m, Policy::AllowOmissions).is_err());
}
#[test]
fn imported_protein_crosslinks_and_nonstandard_residues_cannot_be_omitted() {
    // CRC literal independently calculated with forward ISO polynomial0x1b
    // over each byte's low-to-high bits and a final bit reflection.
    for (kind, location) in [("CROSSLNK", "1..4"), ("NON_STD", "2")] {
        let text = format!(
            "ID   TEST_PROT               Unreviewed;         6 AA.\nAC   P12345;\nFT   {kind:<15} {location}\nFT                   /note=\"retained chemical evidence\"\nSQ   SEQUENCE   6 AA;  0 MW;  69D1B2D37336F000 CRC64;\n     MKWVTF\n//\n"
        );
        let imported = collagene_core::formats::protein::import_bytes(
            text.as_bytes(),
            collagene_core::formats::protein::Profile::UniProt,
        )
        .unwrap();
        let molecule = &imported.molecules[0];
        assert_eq!(molecule.features[0].kind, kind);
        assert!(
            x::export_xdna(molecule, Policy::AllowOmissions).is_err(),
            "{kind}"
        );
        if kind == "CROSSLNK" {
            let mut uncertain = molecule.clone();
            uncertain.features[0].location = Location::Range {
                start: 0,
                end: 4,
                before: true,
                after: true,
            };
            assert!(
                collagene_core::formats::snapgene::export::export_snapgene(
                    &uncertain,
                    collagene_core::formats::snapgene::export::Policy::AllowOmissions
                )
                .is_err()
            );
        }
    }
}
#[test]
fn regenerated_current_edits_preserve_original_source_bytes_without_replaying_old_padding_or_flags()
{
    let mut original = header(1, false, "ACGT", "Original note");
    original[4] = 91;
    original.extend_from_slice(&[77, 1, b'0', 1, b'0', 0]);
    let imported = legacy_formats::import_xdna(&original, "Original filename")
        .unwrap()
        .molecules
        .remove(0);
    let current = imported.replace(1, 2, "TT").unwrap();
    assert_eq!(current.sequence.as_str(), "ATTGT");
    let unchanged = current.clone();
    let report = x::export_xdna(&current, Policy::AllowOmissions).unwrap();
    assert_eq!(current, unchanged);
    assert_eq!(
        legacy_formats::original_legacy_bytes(&current).unwrap(),
        original
    );
    assert_eq!(report.bytes[4], 0);
    assert_eq!(&report.bytes[112..117], b"ATTGT");
    assert_eq!(
        report.bytes[117 + "Original note".len()],
        0,
        "Regenerate extension marker rather than replay77"
    );
    assert!(
        report
            .losses
            .iter()
            .any(|l| l.kind == LossKind::OriginalEvidence)
    );
    let read = legacy_formats::import_xdna(&report.bytes, "New context")
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(read.sequence.as_str(), "ATTGT");
}
#[test]
fn byte_field_and_feature_count_limits_never_wrap_or_silently_discard_a_tail() {
    let mut m = dna();
    let mut f = feature(Location::range(0, 1));
    f.qualifiers = vec![("label".into(), Some("a".repeat(255)))];
    m.features.push(f);
    let report = x::export_xdna(&m, Policy::AllowOmissions).unwrap();
    let read = legacy_formats::import_xdna(&report.bytes, "Bounded label")
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(
        read.features[0].qualifiers[0].1.as_deref(),
        Some("a".repeat(255).as_str())
    );
    m.features = vec![feature(Location::range(0, 1)); 256];
    assert!(x::review(&m).is_err());
    assert!(x::export_xdna(&m, Policy::AllowOmissions).is_err());
}
