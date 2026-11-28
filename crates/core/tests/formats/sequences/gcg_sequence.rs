use collagene_core::gcg_sequence::{self as gcg, ExportOptions, ImportOptions};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
fn options(alphabet: Alphabet) -> ImportOptions {
    ImportOptions {
        alphabet,
        fallback_name: "Explicit fallback".into(),
    }
}
// Wisconsin GetSeq manual's independent eleven-call/4920 oracle. Documentary
// identifier is our own fixture label; checksum does not cover header text.
const LITERAL:&[u8]=b"!!NA_SEQUENCE 1.0\n\nExample  Length: 11  September 29, 1998 18:05  Type: N  Check: 4920  ..\n\n       1  GATTCCGATT G\n";
#[test]
fn literal_wisconsin_oracle_imports_exact_calls_header_and_original_bytes() {
    assert_eq!(gcg::checksum("GATTCCGATTG").unwrap(), 4920);
    assert_eq!(gcg::checksum("ACGTACGTACGT").unwrap(), 5688);
    assert_eq!(gcg::checksum("acgtACGTacgt").unwrap(), 5688);
    assert!(gcg::is_gcg(LITERAL));
    let report = gcg::import_gcg(LITERAL, &options(Alphabet::Dna)).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.name, "Example");
    assert_eq!(m.sequence.as_str(), "GATTCCGATTG");
    assert_eq!(m.metadata["gcg:date_label"], "September 29, 1998 18:05");
    assert_eq!(gcg::original_bytes(m).unwrap(), LITERAL);
    assert!(m.features.is_empty());
    let out = gcg::export_gcg(m, &ExportOptions::default()).unwrap();
    let again = gcg::import_gcg(out.text.as_bytes(), &options(Alphabet::Dna)).unwrap();
    assert_eq!(again.molecules[0].sequence, m.sequence);
    assert_eq!(gcg::original_bytes(m).unwrap(), LITERAL);
}
#[test]
fn independent_checksum_wrap_rna_and_real_protein_stop_are_preserved() {
    assert_eq!(gcg::checksum(&"A".repeat(58)).unwrap(), 7510);
    assert_eq!(gcg::checksum("AUGCU").unwrap(), 1141);
    assert_eq!(gcg::checksum("aaa*").unwrap(), 558);
    let rna=b"!!NA_SEQUENCE 1.0\nRNA Length: 5  November 26, 2026 00:00  Type: N  Check: 1141 ..\n1 AUGCU\n";
    assert!(gcg::import_gcg(rna, &options(Alphabet::Dna)).is_err());
    let m = gcg::import_gcg(rna, &options(Alphabet::Rna))
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(m.sequence.alphabet(), Alphabet::Rna);
    let protein=b"!!AA_SEQUENCE 1.0\nProtein Length: 4  November 26, 2026 00:00  Type: P  Check: 558 ..\n1 aaa*\n";
    assert!(gcg::import_gcg(protein, &options(Alphabet::Dna)).is_err());
    let m = gcg::import_gcg(protein, &options(Alphabet::Protein))
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(m.sequence.as_str(), "AAA*");
    assert_eq!(m.sequence.len(), 4);
    assert!(
        gcg::export_gcg(&m, &ExportOptions::default())
            .unwrap()
            .text
            .contains("AAA*")
    );
}
#[test]
fn corrupt_counts_types_coordinates_alignment_and_annotated_sequences_are_not_silently_stripped() {
    let text = std::str::from_utf8(LITERAL).unwrap();
    for bad in [
        text.replace("4920", "4921"),
        text.replace("Length: 11", "Length: 10"),
        text.replace("Type: N", "Type: P"),
        text.replace("       1", "       2"),
        text.replace("GATTCCGATT G", "GATT-CGATT G"),
        text.replace("GATTCCGATT G", "GATTCCGATT G >annotation<"),
        format!("{text}12\n"),
        format!("{text}{text}"),
    ] {
        assert!(
            gcg::import_gcg(bad.as_bytes(), &options(Alphabet::Dna)).is_err(),
            "{bad}"
        );
    }
    let msf = b"!!NA_MULTIPLE_ALIGNMENT 1.0\n MSF: 11 Type: N Check: 4920 ..\n";
    assert!(!gcg::is_gcg(msf));
    assert!(gcg::import_gcg(msf, &options(Alphabet::Dna)).is_err());
    let absent = text.replace("Type: N", "");
    assert!(gcg::import_gcg(absent.as_bytes(), &options(Alphabet::Dna)).is_err());
}
#[test]
fn synthetic_export_requires_explicit_identity_and_documentary_date_and_reports_loss() {
    let mut m = Molecule::dna("Circular native source", "ACGTACGT", Topology::Circular).unwrap();
    assert!(gcg::export_gcg(&m, &ExportOptions::default()).is_err());
    let opts = ExportOptions {
        record_id: Some("Declared_ID".into()),
        date_label: Some("November 26, 2026 12:00".into()),
    };
    let out = gcg::export_gcg(&m, &opts).unwrap();
    assert!(out.warnings.iter().any(|w| w.contains("Circular topology")));
    assert!(out.warnings.iter().any(|w| w.contains("cloning history")));
    let r = gcg::import_gcg(out.text.as_bytes(), &options(Alphabet::Dna)).unwrap();
    assert_eq!(r.molecules[0].name, "Declared_ID");
    assert_eq!(r.molecules[0].sequence, m.sequence);
    m.metadata
        .insert("source:legacy:chemistry_unmapped".into(), "true".into());
    assert!(gcg::export_gcg(&m, &opts).is_err());
    let mut protein = Molecule::dna("Protein", "A", Topology::Linear).unwrap();
    protein.sequence = Sequence::new(Alphabet::Protein, "AAA*").unwrap();
    assert!(
        gcg::export_gcg(&protein, &opts)
            .unwrap()
            .text
            .contains("Type: P")
    );
}
#[test]
fn legacy_typed_checksum_alias_empty_source_id_and_immutable_original_hash_have_explicit_behavior()
{
    let text = std::str::from_utf8(LITERAL)
        .unwrap()
        .replace("!!NA_SEQUENCE 1.0\n", "")
        .replace("Example  Length:", "Length:")
        .replace("Check:", "Checksum:");
    assert!(gcg::is_gcg(text.as_bytes()));
    let report = gcg::import_gcg(text.as_bytes(), &options(Alphabet::Dna)).unwrap();
    let mut m = report.molecules[0].clone();
    assert_eq!(m.name, "Explicit fallback");
    assert!(report.warnings.iter().any(|w| w.contains("legacy")));
    assert!(gcg::export_gcg(&m, &ExportOptions::default()).is_err());
    m.metadata
        .insert("gcg:original:hash".into(), "00".repeat(32));
    assert!(gcg::original_bytes(&m).is_err());
    let bad = ExportOptions {
        record_id: Some("bad id".into()),
        date_label: Some("November 26, 2026".into()),
    };
    assert!(gcg::export_gcg(&m, &bad).is_err());
}
#[test]
fn source_header_and_sequence_limits_fail_without_partial_document() {
    assert!(gcg::import_gcg(&vec![b'A'; gcg::MAX_BYTES + 1], &options(Alphabet::Dna)).is_err());
    let mut header = vec![b'x'; 65537];
    header.extend_from_slice(LITERAL);
    assert!(gcg::import_gcg(&header, &options(Alphabet::Dna)).is_err());
    assert!(gcg::checksum("AC-GT").is_err());
    assert!(gcg::checksum("é").is_err());
}

#[test]
fn actual_declared_termini_product_and_all_restrictive_documentary_shapes_cannot_be_laundered() {
    use collagene_core::{
        chemistry::Terminus,
        verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    };
    let stock = Molecule::dna("Original intact reference", "ACGTACGT", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent same-sequence target",
        "ACGTACGT",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actually declare native end chemistry".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock.clone(),
        }],
        steps: vec![Operation::DeclareSourceTermini {
            template: Source::Input("stock".into()),
            left: Terminus::CLEAVED,
            right: Terminus::CLEAVED,
        }],
    };
    let report = verification::verify(&plan, &target);
    assert_ne!(report.status, Status::Invalid, "{:?}", report.findings);
    assert_eq!(report.products.len(), 1);
    let product = &report.products[0];
    assert!(product.metadata.contains_key("source:declared_termini"));
    let options = ExportOptions {
        record_id: Some("Declared".into()),
        date_label: Some("November 26, 2026 12:00".into()),
    };
    assert!(gcg::export_gcg(product, &options).is_err());
    assert!(gcg::export_gcg(&stock, &options).is_ok());
    for key in [
        "source:declared_termini",
        "source:restriction_inspection:recipe",
        "source:restriction_inspection:duplex_summary",
        "restriction_catalog:selection",
        "source:restriction_catalog:selection",
        "tracked_modifications",
        "protein:modifications",
    ] {
        let mut evidence = stock.clone();
        if key == "protein:modifications" {
            evidence.sequence = Sequence::new(Alphabet::Protein, "ACGT").unwrap();
        }
        evidence.metadata.insert(
            key.into(),
            "documentary evidence is not a fresh native state".into(),
        );
        assert!(gcg::export_gcg(&evidence, &options).is_err(), "{key}");
    }
    for marker in [
        "Collagene physical product v1: unverified sequence-only evidence",
        "Collagene reference-only v1: unverified sequence-only evidence",
    ] {
        let mut evidence = stock.clone();
        evidence.original_genbank = Some(format!("COMMENT     {marker}\n"));
        assert!(gcg::export_gcg(&evidence, &options).is_err());
    }
    // The shared gate must not turn into linear-unambiguous-DNA declaration
    // preflight: all ordinary GCG alphabet/topology profiles remain available.
    let mut rna = stock.clone();
    rna.sequence = Sequence::new(Alphabet::Rna, "ACGUN").unwrap();
    assert!(gcg::export_gcg(&rna, &options).is_ok());
    let mut circular = stock.clone();
    circular.topology = Topology::Circular;
    circular.sequence = Sequence::new(Alphabet::Dna, "ACGTN").unwrap();
    assert!(gcg::export_gcg(&circular, &options).is_ok());
}

#[test]
fn protein_chemical_features_refuse_sequence_only_laundering_while_ordinary_domains_report_loss() {
    use collagene_model::{Feature, Location};
    let mut protein = Molecule::dna("Annotated protein", "A", Topology::Linear).unwrap();
    protein.sequence = Sequence::new(Alphabet::Protein, "CCAAA").unwrap();
    let options = ExportOptions {
        record_id: Some("Protein".into()),
        date_label: Some("November 26, 2026 12:00".into()),
    };
    for kind in [
        "MOD_RES",
        "DISULFID",
        "modified_residue",
        "cross-link",
        "cross_link",
        "bond",
        "LIPID",
        "CARBOHYD",
    ] {
        let mut modified = protein.clone();
        modified.features.push(Feature {
            kind: kind.into(),
            location: if kind == "DISULFID" {
                Location::Join(vec![Location::range(0, 1), Location::range(1, 2)])
            } else {
                Location::range(2, 3)
            },
            qualifiers: vec![(
                "note".into(),
                Some("Retained chemical annotation; not an unmodified free chain".into()),
            )],
            needs_review: false,
        });
        assert!(gcg::export_gcg(&modified, &options).is_err(), "{kind}");
    }
    protein.features.push(Feature {
        kind: "DOMAIN".into(),
        location: Location::range(0, 5),
        qualifiers: vec![("label".into(), Some("ordinary sequence domain".into()))],
        needs_review: false,
    });
    let exported = gcg::export_gcg(&protein, &options).unwrap();
    assert!(exported.warnings.iter().any(|w| w.contains("features")));
    let imported = gcg::import_gcg(exported.text.as_bytes(), &options_for_protein()).unwrap();
    assert_eq!(imported.molecules[0].sequence, protein.sequence);
    assert!(imported.molecules[0].features.is_empty());
    fn options_for_protein() -> ImportOptions {
        ImportOptions {
            alphabet: Alphabet::Protein,
            fallback_name: "Protein".into(),
        }
    }
}
