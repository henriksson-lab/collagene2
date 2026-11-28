use collagene_core::sequence::annotations::merge;
use collagene_core::{
    formats::sequence::{self as formats, TextFormat},
    protocol,
    sequence::annotations::library::{FeatureLibrary, ScanOptions, scan},
};
use collagene_model::{Feature, Location, Molecule, Topology};

fn merge_request(
    target: &Molecule,
    source: merge::Source,
    selections: Vec<merge::Selection>,
) -> merge::Request {
    merge::Request {
        schema_version: 1,
        target_revision: target.revision(),
        source,
        selections,
        limits: Default::default(),
    }
}
#[test]
fn library_machine_identifiers_survive_genbank_and_embl_physical_wrapping() {
    let mut source = Molecule::dna("source", "ACGTT", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "promoter".into(),
        location: Location::range(0, 5),
        qualifiers: vec![("label".into(), Some("Template".into()))],
        needs_review: false,
    });
    let mut library = FeatureLibrary::new("Lab library α with spaces", "1").unwrap();
    library.add_feature(&source, 0).unwrap();
    let target = Molecule::dna("target", "GGGACGTTCCC", Topology::Linear).unwrap();
    let report = scan(&target, &library, &ScanOptions::default()).unwrap();
    assert_eq!(report.hits.len(), 1);
    let annotated = merge::apply(
        &target,
        &merge_request(
            &target,
            merge::Source::Library {
                library: (&library).clone(),
                report: (&report).clone(),
            },
            (&[report.hits[0].id.clone()])
                .iter()
                .map(|id| merge::Selection {
                    proposal_id: id.clone(),
                    decision: merge::Decision::Append,
                })
                .collect(),
        ),
    )
    .map(|applied| applied.molecule)
    .unwrap();
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
        let original = merge::Ledger::from_feature(&annotated.features[0]).unwrap();
        let restored = merge::Ledger::from_feature(&reopened.features[0]).unwrap();
        assert_eq!(original.contributions.len(), 1);
        assert_eq!(restored, original);
    }
}
#[test]
fn text_exports_report_omitted_quality_and_annotation_sources() {
    let read = collagene_core::formats::oligos::import_fastq(
        b"@read\nACGT\n+\nIIII\n",
        collagene_core::formats::oligos::QualityEncoding::Phred33,
    )
    .unwrap()
    .molecules
    .remove(0);
    let target = Molecule::dna("ref", "ACGT", Topology::Linear).unwrap();
    let annotations = collagene_core::formats::annotations::import_annotations(
        &target,
        "ref\t0\t4\tname\t0\t+\n",
        collagene_core::formats::annotations::AnnotationFormat::Bed,
        "ref",
    )
    .unwrap();
    let annotated = merge::apply(
        &target,
        &merge_request(
            &target,
            merge::Source::Dataset {
                report: (&annotations).clone(),
            },
            (&[0])
                .iter()
                .map(|index| merge::Selection {
                    proposal_id: (&annotations)
                        .proposals
                        .get(*index)
                        .map(|proposal| proposal.id.clone())
                        .unwrap_or_else(|| format!("unavailable-dataset-proposal-{index}")),
                    decision: merge::Decision::Append,
                })
                .collect(),
        ),
    )
    .map(|applied| applied.molecule)
    .unwrap();
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
    use collagene_core::formats::annotations::{self as a, AnnotationFormat};
    let source = "ref\t0\t3\tshort\t0\t+\n";
    let reference = Molecule::dna("ref", "ACGTACGT", Topology::Linear).unwrap();
    let report = a::import_annotations(&reference, source, AnnotationFormat::Bed, "ref").unwrap();
    let annotated = merge::apply(
        &reference,
        &merge_request(
            &reference,
            merge::Source::Dataset {
                report: (&report).clone(),
            },
            (&[0])
                .iter()
                .map(|index| merge::Selection {
                    proposal_id: (&report)
                        .proposals
                        .get(*index)
                        .map(|proposal| proposal.id.clone())
                        .unwrap_or_else(|| format!("unavailable-dataset-proposal-{index}")),
                    decision: merge::Decision::Append,
                })
                .collect(),
        ),
    )
    .map(|applied| applied.molecule)
    .unwrap();
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
