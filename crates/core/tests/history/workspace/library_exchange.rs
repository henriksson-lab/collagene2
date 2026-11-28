use collagene_core::{
    feature_library::{FeatureLibrary, ScanOptions, apply_selected, scan},
    formats::{self, TextFormat},
    *,
};
use collagene_model::{Feature, Location, Molecule, Topology};

#[test]
fn library_machine_identifiers_survive_genbank_and_embl_physical_wrapping() {
    let mut source = Molecule::dna("source", "ACGTT", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "promoter".into(),
        location: Location::range(0, 5),
        qualifiers: vec![("label".into(), Some("Template".into()))],
        needs_review: false,
    });
    let mut library = FeatureLibrary::new("Lab library", "1").unwrap();
    library.add_feature(&source, 0).unwrap();
    let target = Molecule::dna("target", "GGGACGTTCCC", Topology::Linear).unwrap();
    let report = scan(&target, &library, &ScanOptions::default()).unwrap();
    assert_eq!(report.hits.len(), 1);
    let annotated =
        apply_selected(&target, &library, &report, &[report.hits[0].id.clone()]).unwrap();
    for format in [TextFormat::GenBank, TextFormat::Embl] {
        let output = formats::export_with_report(&annotated, format).unwrap();
        assert!(
            output
                .warnings
                .iter()
                .any(|warning| warning.contains("Complete library templates"))
        );
        let reopened = formats::import_text(&output.text)
            .unwrap()
            .molecules
            .remove(0);
        assert!(
            reopened.features[0].needs_review,
            "inferred annotation lost its review flag"
        );
        for (key, value) in annotated.features[0].qualifiers.iter().filter(|(key, _)| {
            matches!(
                key.as_str(),
                "library_entry_id"
                    | "library_entry_revision"
                    | "library_revision"
                    | "library_source_revision"
                    | "library_proposal_id"
            )
        }) {
            assert_eq!(
                reopened.features[0]
                    .qualifiers
                    .iter()
                    .find(|(found, _)| found == key)
                    .unwrap()
                    .1,
                *value,
                "{key}"
            );
        }
    }
}
#[test]
fn text_exports_report_omitted_quality_and_annotation_sources() {
    let read = collagene_core::oligo_formats::import_fastq(
        b"@read\nACGT\n+\nIIII\n",
        collagene_core::oligo_formats::QualityEncoding::Phred33,
    )
    .unwrap()
    .molecules
    .remove(0);
    let target = Molecule::dna("ref", "ACGT", Topology::Linear).unwrap();
    let annotations = collagene_core::annotation_formats::import_annotations(
        &target,
        "ref\t0\t4\tname\t0\t+\n",
        collagene_core::annotation_formats::AnnotationFormat::Bed,
        "ref",
    )
    .unwrap();
    let annotated =
        collagene_core::annotation_formats::apply_annotations(&target, &annotations, &[0]).unwrap();
    for format in [TextFormat::GenBank, TextFormat::Embl, TextFormat::Fasta] {
        assert!(
            formats::export_with_report(&read, format)
                .unwrap()
                .warnings
                .iter()
                .any(|warning| warning.contains("per-base qualities"))
        );
        assert!(
            formats::export_with_report(&annotated, format)
                .unwrap()
                .warnings
                .iter()
                .any(|warning| warning.contains("Original GFF3/BED"))
        );
    }
}

#[test]
fn original_annotation_download_survives_edits_and_rejects_tampered_metadata() {
    use collagene_core::annotation_formats::{self as a, AnnotationFormat};
    let source = "ref\t0\t3\tshort\t0\t+\n";
    let reference = Molecule::dna("ref", "ACGTACGT", Topology::Linear).unwrap();
    let report = a::import_annotations(&reference, source, AnnotationFormat::Bed, "ref").unwrap();
    let annotated = a::apply_annotations(&reference, &report, &[0]).unwrap();
    let edited = annotated.replace(1, 2, "TT").unwrap();
    assert_eq!(
        a::original_annotation_source(&edited, &report.source_blake3).unwrap(),
        source.as_bytes()
    );
    let reopened = protocol::load_document(&protocol::save_document(&edited).unwrap()).unwrap();
    assert_eq!(
        a::original_annotation_source(&reopened, &report.source_blake3).unwrap(),
        source.as_bytes()
    );
    let mut tampered = reopened;
    tampered.metadata.insert(
        format!("annotation_import:{}:original_base64", report.source_blake3),
        "QUJDRA==".into(),
    );
    assert!(a::original_annotation_source(&tampered, &report.source_blake3).is_err());
}
