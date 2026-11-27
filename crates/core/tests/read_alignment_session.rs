use collagene_core::{
    Molecule, Topology,
    protocol::{self, AnalysisJob as Job, AnalysisResult as ResultValue, AnalysisSession, Request},
    read_alignment_reference::Options,
    read_alignment_session::*,
};
const SAM: &str = "@SQ\tSN:ref\tLN:8\nr1\t0\tref\t1\t30\t4M\t*\t0\t0\tACGT\t!III\tNM:i:0\nr2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
fn encoded(revision: String, job: Job) -> Vec<u8> {
    protocol::encode_request(&Request {
        id: 17,
        revision,
        job,
    })
    .unwrap()
}
fn run(session: &mut AnalysisSession, revision: String, job: Job) -> (ResultValue, u32) {
    let request = encoded(revision.clone(), job);
    assert_eq!(
        &request[4..6],
        &protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = session.process_request(&request).unwrap();
    let ((id, binding, tag), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!((id, binding), (17, revision.as_str()));
    (protocol::decode_response(&response).unwrap().result, tag)
}
fn prepare(session: &mut AnalysisSession, bytes: &[u8]) -> Prepared {
    let revision = prepare_revision(bytes, FileFormat::Sam).unwrap();
    let (value, tag) = run(
        session,
        revision,
        Job::PrepareReadAlignments {
            bytes: bytes.to_vec(),
            format: FileFormat::Sam,
        },
    );
    assert_eq!(tag, 69);
    let ResultValue::ReadAlignmentsPrepared(value) = value else {
        panic!("expected prepared")
    };
    value
}

#[test]
fn complete_private_source_pages_fresh_checks_exports_and_disposal_share_exact_binding() {
    let mut session = AnalysisSession::default();
    let prepared = prepare(&mut session, SAM.as_bytes());
    assert_eq!(prepared.summary.records, 2);
    assert_eq!(prepared.summary.flag_unmapped, 1);
    assert_eq!(
        prepared.summary.source_blake3,
        blake3::hash(SAM.as_bytes()).to_hex().as_str()
    );
    let source = prepared.summary.source_blake3;
    let (page, tag) = run(
        &mut session,
        page_revision(&source, 0, 2).unwrap(),
        Job::ReadAlignmentPage {
            source_hash: source.clone(),
            start: 0,
            count: 2,
        },
    );
    assert_eq!(tag, 70);
    let ResultValue::ReadAlignmentPage(page) = page else {
        panic!("expected page")
    };
    assert_eq!(page.records[0].sequence.as_deref(), Some("ACGT"));
    assert_eq!(page.records[0].qualities, Some(vec![0, 40, 40, 40]));
    assert!(page.records[1].sequence.is_none() && page.records[1].qualities.is_none());
    let reference = Molecule::dna("manual source", "ACGTACGT", Topology::Linear).unwrap();
    let options = Options {
        reference_name: "ref".into(),
        expected_reference_revision: reference.revision(),
        ..Options::default()
    };
    let revision = check_revision(&source, 0, &reference, &options).unwrap();
    let (checked, tag) = run(
        &mut session,
        revision.clone(),
        Job::CheckReadAlignment {
            source_hash: source.clone(),
            index: 0,
            reference: reference.clone(),
            options: options.clone(),
        },
    );
    assert_eq!(tag, 71);
    let ResultValue::ReadAlignmentChecked(checked) = checked else {
        panic!("expected checked")
    };
    checked
        .validate_bindings(
            &revision,
            &source,
            0,
            &page.records[0].original_blake3,
            &reference.revision(),
            &options,
        )
        .unwrap();
    assert!(checked.check.report.complete);
    assert_eq!(checked.check.report.summary.certain_agreements, 4);
    let (export, tag) = run(
        &mut session,
        export_revision(&source, FileFormat::Sam, 4096).unwrap(),
        Job::ExportReadAlignments {
            source_hash: source.clone(),
            format: FileFormat::Sam,
            max_output_bytes: 4096,
        },
    );
    assert_eq!(tag, 72);
    let ResultValue::ReadAlignmentsExported(export) = export else {
        panic!("expected export")
    };
    assert_eq!(export.bytes, SAM.as_bytes());
    let (_, tag) = run(
        &mut session,
        dispose_revision(&source).unwrap(),
        Job::DisposeReadAlignments {
            source_hash: source.clone(),
        },
    );
    assert_eq!(tag, 73);
    assert!(
        session
            .process_request(&encoded(
                page_revision(&source, 0, 1).unwrap(),
                Job::ReadAlignmentPage {
                    source_hash: source,
                    start: 0,
                    count: 1
                }
            ))
            .is_err()
    );
}

#[test]
fn previews_are_explicitly_truncated_without_filling_missing_calls_and_export_never_publishes_partial_bytes()
 {
    let tags = (b'a'..=b'b')
        .flat_map(|a| (b'A'..=b'T').map(move |b| format!("\t{}{}:i:1", a as char, b as char)))
        .collect::<String>();
    let text = format!(
        "@SQ\tSN:ref\tLN:130\nr\t0\tref\t1\t30\t{}\t*\t0\t0\t{}\t{}{}\n",
        "1M1I".repeat(130),
        "A".repeat(260),
        "!".repeat(260),
        tags
    );
    let mut session = AnalysisSession::default();
    let prepared = prepare(&mut session, text.as_bytes());
    let source = prepared.summary.source_blake3;
    let (page, _) = run(
        &mut session,
        page_revision(&source, 0, 100).unwrap(),
        Job::ReadAlignmentPage {
            source_hash: source.clone(),
            start: 0,
            count: 100,
        },
    );
    let ResultValue::ReadAlignmentPage(page) = page else {
        panic!("expected page")
    };
    let record = &page.records[0];
    assert_eq!(
        (
            record.sequence.as_ref().unwrap().len(),
            record.sequence_length
        ),
        (200, Some(260))
    );
    assert_eq!(
        (
            record.qualities.as_ref().unwrap().len(),
            record.qualities_length
        ),
        (200, Some(260))
    );
    assert_eq!((record.cigar.len(), record.cigar_count), (64, 260));
    assert_eq!((record.tags.len(), record.tag_count), (32, 40));
    assert!(
        record.truncated.sequence
            && record.truncated.qualities
            && record.truncated.cigar
            && record.truncated.tags
    );
    assert!(
        session
            .process_request(&encoded(
                export_revision(&source, FileFormat::Sam, 16).unwrap(),
                Job::ExportReadAlignments {
                    source_hash: source.clone(),
                    format: FileFormat::Sam,
                    max_output_bytes: 16
                }
            ))
            .is_err()
    );
    assert!(
        session
            .process_request(&encoded(
                page_revision(&source, 0, 1).unwrap(),
                Job::ReadAlignmentPage {
                    source_hash: source,
                    start: 0,
                    count: 1
                }
            ))
            .is_ok()
    );
}

#[test]
fn failed_imports_and_decoded_summaries_cannot_replace_or_recreate_private_source_authority() {
    let mut session = AnalysisSession::default();
    let prepared = prepare(&mut session, SAM.as_bytes());
    let source = prepared.summary.source_blake3.clone();
    let bad = b"@SQ\tSN:ref\tLN:8\nbroken\n";
    assert!(
        session
            .process_request(&encoded(
                prepare_revision(bad, FileFormat::Sam).unwrap(),
                Job::PrepareReadAlignments {
                    bytes: bad.to_vec(),
                    format: FileFormat::Sam
                }
            ))
            .is_err()
    );
    let page = encoded(
        page_revision(&source, 0, 1).unwrap(),
        Job::ReadAlignmentPage {
            source_hash: source.clone(),
            start: 0,
            count: 1,
        },
    );
    assert!(session.process_request(&page).is_ok());
    let historical: Prepared =
        serde_json::from_slice(&serde_json::to_vec(&prepared).unwrap()).unwrap();
    historical.validate().unwrap();
    assert!(AnalysisSession::default().process_request(&page).is_err());
    assert!(
        protocol::process_request(&page)
            .unwrap_err()
            .to_string()
            .contains("persistent worker")
    );
    assert!(page_revision(&source, 0, 101).is_err());
}

#[test]
fn current_job_tags_and_unsupported_header_admission_remain_explicit() {
    for tag in 77u32..=81 {
        let payload = postcard::to_allocvec(&(17u32, "context", tag)).unwrap();
        let mut bytes = b"CGQ\0".to_vec();
        bytes.extend(23u16.to_le_bytes());
        bytes.extend((payload.len() as u32).to_le_bytes());
        bytes.extend(blake3::hash(&payload).as_bytes());
        bytes.extend(payload);
        for result in [
            protocol::process_request(&bytes),
            AnalysisSession::default().process_request(&bytes),
        ] {
            assert!(matches!(result, Err(collagene_core::Error::Unsupported(_))));
        }
    }
    let source = blake3::hash(SAM.as_bytes()).to_hex().to_string();
    let jobs = [
        Job::PrepareReadAlignments {
            bytes: SAM.as_bytes().to_vec(),
            format: FileFormat::Sam,
        },
        Job::ReadAlignmentPage {
            source_hash: source.clone(),
            start: 0,
            count: 1,
        },
        Job::CheckReadAlignment {
            source_hash: source.clone(),
            index: 0,
            reference: Molecule::dna("ref", "ACGTACGT", Topology::Linear).unwrap(),
            options: Options::default(),
        },
        Job::ExportReadAlignments {
            source_hash: source.clone(),
            format: FileFormat::Sam,
            max_output_bytes: 4096,
        },
        Job::DisposeReadAlignments {
            source_hash: source,
        },
    ];
    for (offset, job) in jobs.into_iter().enumerate() {
        let bytes = encoded("context".into(), job);
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 77 + offset as u32);
    }
}
