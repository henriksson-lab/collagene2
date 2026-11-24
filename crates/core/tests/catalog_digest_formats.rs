//! Catalog identities survive annotated text as bounded producer claims, never chemical authority.
use collagene_core::{
    Molecule, Topology,
    catalog_digest::*,
    chemistry,
    formats::{self, TextFormat},
    restriction_inspection::InspectionOptions,
};
fn product_molecule() -> (Molecule, CatalogSelection) {
    let c = builtin_catalog().unwrap();
    let s = select(&c, &["haeiii".into()]).unwrap();
    let m = Molecule::dna("catalog source", "AAGGCCTT", Topology::Linear).unwrap();
    let p = product(&c, &s, "source-role", &m, &InspectionOptions::default(), 0).unwrap();
    (p.product.molecule, s)
}
#[test]
fn genbank_and_embl_preserve_full_catalog_pin_without_restoring_physical_authority() {
    let (m, s) = product_molecule();
    let pin = m.metadata[RECIPE_KEY].clone();
    for format in [TextFormat::GenBank, TextFormat::Embl] {
        let output = formats::export_with_report(&m, format).unwrap();
        assert!(output.text.contains("CG-CATALOG-DIGEST-v1"));
        assert!(output.text.contains("unverified producer"));
        assert!(output.warnings.iter().any(|w| w.contains("Catalog-pinned")));
        let imported = formats::import_text(&output.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(imported.sequence.as_str(), "AAGG");
        assert_eq!(imported.metadata[IMPORTED_RECIPE_KEY], pin);
        assert_eq!(
            imported.metadata["source:restriction_catalog:engine"],
            m.metadata[ENGINE_KEY]
        );
        assert_eq!(
            imported.metadata["source:restriction_catalog:batch_revision"],
            m.metadata[PRODUCT_REVISION_KEY]
        );
        assert_eq!(unverified_selection(&imported).unwrap(), Some(s.clone()));
        assert!(!imported.metadata.contains_key(RECIPE_KEY));
        assert_eq!(
            imported.metadata["source:legacy:chemistry_unmapped"],
            "true"
        );
        assert!(chemistry::Duplex::from_molecule(&imported).is_err());
        let repeated = formats::export_with_report(&imported, format).unwrap();
        let second = formats::import_text(&repeated.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(second.metadata[IMPORTED_RECIPE_KEY], pin);
        assert!(chemistry::Duplex::from_molecule(&second).is_err());
    }
    let fasta = formats::export_with_report(&m, TextFormat::Fasta).unwrap();
    assert!(!fasta.text.contains("CG-CATALOG"));
    assert!(fasta.warnings.iter().any(|w| w.contains("Catalog-pinned")));
}
#[test]
fn malformed_and_oversized_catalog_comments_remain_explicitly_unverified_bounded_claims() {
    let (m, _) = product_molecule();
    let gb = formats::export_with_report(&m, TextFormat::GenBank)
        .unwrap()
        .text;
    let corrupt = gb
        .lines()
        .map(|line| {
            if line.starts_with("COMMENT     CG-CATALOG-DIGEST-v1 ") {
                "COMMENT     CG-CATALOG-DIGEST-v1 0001/0002 eA=="
            } else {
                line
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let imported = formats::import_text(&corrupt).unwrap().molecules.remove(0);
    assert!(!imported.metadata.contains_key(IMPORTED_RECIPE_KEY));
    assert_eq!(
        imported.metadata["source:restriction_inspection:invalid_producer_claims"],
        "true"
    );
    assert!(chemistry::Duplex::from_molecule(&imported).is_err());
    let mut oversized =
        Molecule::dna("opaque oversized producer", "ACGT", Topology::Linear).unwrap();
    oversized
        .metadata
        .insert(RECIPE_KEY.into(), "é".repeat(20_000));
    for format in [TextFormat::GenBank, TextFormat::Embl] {
        let output = formats::export_with_report(&oversized, format).unwrap();
        assert!(output.text.contains("producer data truncated"));
        assert!(
            output.text.len() < 128 * 1024,
            "32KiB claim plus base64/COMMENT framing: {} bytes",
            output.text.len()
        );
        let imported = formats::import_text(&output.text)
            .unwrap()
            .molecules
            .remove(0);
        assert!(imported.metadata[IMPORTED_RECIPE_KEY].len() <= 32 * 1024);
        assert!(unverified_selection(&imported).is_err());
        assert!(chemistry::Duplex::from_molecule(&imported).is_err());
    }
}
