//! GTF adds a nested format selector, not a new outer job or report layout.
use collagene_core::{
    Molecule, Topology,
    annotation_formats::{self, AnnotationFormat},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};

const GTF: &str = "chr_demo\tfixture\texon\t2\t5\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";\n";
fn source() -> Molecule {
    Molecule::dna("actual reference", "AACCGGTT", Topology::Linear).unwrap()
}
fn request(job: AnalysisJob) -> Vec<u8> {
    protocol::encode_request(&Request {
        id: 91,
        revision: source().revision(),
        job,
    })
    .unwrap()
}
fn historical(mut bytes: Vec<u8>, version: u16) -> Vec<u8> {
    bytes[4..6].copy_from_slice(&version.to_le_bytes());
    bytes
}
fn response(result: AnalysisResult, version: u16) -> Vec<u8> {
    let value = Response {
        id: 91,
        revision: source().revision(),
        result,
    };
    let payload = postcard::to_allocvec(&value).unwrap();
    let mut bytes = b"CGR\0".to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    bytes
}

#[test]
fn gtf_import_apply_export_requires_thirty_without_changing_old_tags() {
    let molecule = source();
    let imported =
        annotation_formats::import_annotations(&molecule, GTF, AnnotationFormat::Gtf, "chr_demo")
            .unwrap();
    let jobs = vec![
        AnalysisJob::ImportAnnotations {
            molecule: molecule.clone(),
            text: GTF.into(),
            format: AnnotationFormat::Gtf,
            seqid: "chr_demo".into(),
        },
        AnalysisJob::ApplyAnnotations {
            molecule: molecule.clone(),
            report: imported.clone(),
            selected: vec![0],
        },
        AnalysisJob::ExportAnnotations {
            molecule: molecule.clone(),
            format: AnnotationFormat::Gtf,
            seqid: "chr_demo".into(),
        },
    ];
    for (offset, job) in jobs.into_iter().enumerate() {
        let bytes = request(job);
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 24 + offset as u32);
        let output = protocol::process_request(&bytes).unwrap();
        assert!(!matches!(
            protocol::decode_response(&output).unwrap().result,
            AnalysisResult::Error(_)
        ));
        for old in [7, 28, 29] {
            assert!(protocol::process_request(&historical(bytes.clone(), old)).is_err());
        }
    }
    assert_eq!(
        imported.proposals[0].feature.location.spans(),
        vec![(1, 5, true)]
    );
    assert_eq!(imported.original_bytes().unwrap(), GTF.as_bytes());
    for old in [7, 29] {
        assert!(
            protocol::decode_response(&response(
                AnalysisResult::AnnotationImport(imported.clone()),
                old
            ))
            .is_err()
        );
    }
    assert!(
        protocol::decode_response(&response(AnalysisResult::AnnotationImport(imported), 30))
            .is_ok()
    );
}

#[test]
fn historical_gff_and_bed_keep_their_original_binary_selectors() {
    assert_eq!(
        postcard::to_allocvec(&AnnotationFormat::Gff3).unwrap(),
        vec![0]
    );
    assert_eq!(
        postcard::to_allocvec(&AnnotationFormat::Bed).unwrap(),
        vec![1]
    );
    assert_eq!(
        postcard::to_allocvec(&AnnotationFormat::Gtf).unwrap(),
        vec![2]
    );
    for (format, text) in [
        (
            AnnotationFormat::Gff3,
            "##gff-version 3\nchr_demo\tfixture\texon\t2\t5\t.\t-\t.\tID=e1\n",
        ),
        (AnnotationFormat::Bed, "chr_demo\t1\t5\texon\t0\t-\n"),
    ] {
        let bytes = historical(
            request(AnalysisJob::ImportAnnotations {
                molecule: source(),
                text: text.into(),
                format,
                seqid: "chr_demo".into(),
            }),
            7,
        );
        let output = protocol::process_request(&bytes).unwrap();
        let AnalysisResult::AnnotationImport(report) =
            protocol::decode_response(&output).unwrap().result
        else {
            panic!("expected annotation report")
        };
        assert_eq!(report.format, format);
        assert_eq!(
            report.proposals[0].feature.location.spans(),
            vec![(1, 5, true)]
        );
    }
}
