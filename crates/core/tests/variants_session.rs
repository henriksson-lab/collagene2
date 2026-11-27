use collagene_core::{
    Molecule, Topology,
    protocol::{self, AnalysisJob as Job, AnalysisResult as Value, AnalysisSession, Request},
    variant_reference::{Options, Outcome},
    variants_session::*,
};
const VCF:&[u8]=b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=7>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\t.\tCC\tCTC\t.\tPASS\t.\n";
fn run(s: &mut AnalysisSession, revision: String, job: Job) -> Value {
    let bytes = protocol::encode_request(&Request {
        id: 19,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    let response = s.process_request(&bytes).unwrap();
    let r = protocol::decode_response(&response).unwrap();
    assert_eq!(r.revision, revision);
    r.result
}
#[test]
fn private_source_pages_fresh_reference_and_preview_then_exact_original_export() {
    let mut s = AnalysisSession::default();
    let Value::VariantsPrepared(p) = run(
        &mut s,
        prepare_revision(VCF, FileFormat::Vcf).unwrap(),
        Job::PrepareVariants {
            bytes: VCF.to_vec(),
            format: FileFormat::Vcf,
        },
    ) else {
        panic!("prepared")
    };
    assert_eq!(p.summary.records, 1);
    assert_eq!(p.summary.contigs[0].length, Some(7));
    let source = p.summary.source_hash;
    let Value::VariantPage(page) = run(
        &mut s,
        page_revision(&source, 0, 100).unwrap(),
        Job::VariantPage {
            source_hash: source.clone(),
            start: 0,
            count: 100,
        },
    ) else {
        panic!("page")
    };
    assert_eq!(page.rows[0].position_1based, Some(3));
    let reference = Molecule::dna("manual mapping", "AACCGTT", Topology::Linear).unwrap();
    let options = Options {
        source_blake3: source.clone(),
        record_index: 0,
        reference_name: "chr1".into(),
        expected_reference_revision: reference.revision(),
        alternate_index: 1,
    };
    let Value::VariantReferenceChecked(checked) = run(
        &mut s,
        check_revision(&source, &reference, &options).unwrap(),
        Job::CheckVariantReference {
            source_hash: source.clone(),
            reference: reference.clone(),
            options: options.clone(),
        },
    ) else {
        panic!("check")
    };
    checked.validate_bindings(&source, &options).unwrap();
    assert_eq!(checked.report.outcome, Outcome::ReferenceAlleleMatches);
    let Value::VariantAllelePreview(preview) = run(
        &mut s,
        preview_revision(&source, &reference, &options).unwrap(),
        Job::PreviewVariantAllele {
            source_hash: source.clone(),
            reference,
            options: options.clone(),
        },
    ) else {
        panic!("preview")
    };
    preview.validate_bindings(&source, &options).unwrap();
    assert_eq!(preview.preview.molecule.sequence.as_str(), "AACTCGTT");
    let Value::VariantsExported(export) = run(
        &mut s,
        export_revision(&source, ExportFormat::Original, VCF.len()).unwrap(),
        Job::ExportVariants {
            source_hash: source.clone(),
            format: ExportFormat::Original,
            max_output_bytes: VCF.len(),
        },
    ) else {
        panic!("export")
    };
    assert_eq!(export.bytes, VCF);
    run(
        &mut s,
        dispose_revision(&source).unwrap(),
        Job::DisposeVariants {
            source_hash: source.clone(),
        },
    );
    let job = Job::VariantPage {
        source_hash: source.clone(),
        start: 0,
        count: 1,
    };
    let frame = protocol::encode_request(&Request {
        id: 1,
        revision: page_revision(&source, 0, 1).unwrap(),
        job,
    })
    .unwrap();
    assert!(s.process_request(&frame).is_err());
}
#[test]
fn append_only_tags_and_legacy_version_cannot_prepare_or_restore_private_state() {
    let request = Request {
        id: 1,
        revision: prepare_revision(VCF, FileFormat::Vcf).unwrap(),
        job: Job::PrepareVariants {
            bytes: VCF.to_vec(),
            format: FileFormat::Vcf,
        },
    };
    let mut frame = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&frame[42..]).unwrap();
    assert_eq!(tag, 86);
    frame[4..6].copy_from_slice(&28u16.to_le_bytes());
    assert!(AnalysisSession::default().process_request(&frame).is_err());
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
    let reference = Molecule::dna("ref", "AACCGTT", Topology::Linear).unwrap();
    let options = Options {
        source_blake3: "0".repeat(64),
        record_index: 0,
        reference_name: "chr1".into(),
        expected_reference_revision: reference.revision(),
        alternate_index: 1,
    };
    assert!(check_revision(&"1".repeat(64), &reference, &options).is_err());
}
