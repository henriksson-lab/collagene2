//! Literal export selectors preserve old binary prefixes and bind source/options.
use collagene_core::{
    Molecule, Topology, gcg_sequence,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    scf_export,
    sequence_export::{self, Format},
    traces,
};

fn tag(value: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(value).unwrap())
        .unwrap()
        .0
}
fn run(
    job: AnalysisJob,
    revision: String,
    request_tag: u32,
    result_tag: u32,
) -> protocol::Response {
    assert_eq!(tag(&job), request_tag);
    let request = protocol::encode_request(&Request {
        id: 17,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    let mut old = request.clone();
    old[4..6].copy_from_slice(&89u16.to_le_bytes());
    assert_eq!(&old[6..], &request[6..]);
    assert!(protocol::process_request(&old).is_err());
    let bytes = protocol::process_request(&request).unwrap();
    let response = protocol::decode_response(&bytes).unwrap();
    assert_eq!(response.revision, revision);
    assert_eq!(tag(&response.result), result_tag, "{:?}", response.result);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&89u16.to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(protocol::decode_response(&old).is_err());
    response
}
#[test]
fn sequence_export_review_and_preparation_are_source_bound_append_only_profiles() {
    let source = Molecule::dna("Actual reference", "ACGTCGAT", Topology::Circular).unwrap();
    let before = source.clone();
    for format in [
        Format::SnapGene,
        Format::Xdna,
        Format::Gcg {
            options: gcg_sequence::ExportOptions {
                record_id: Some("REFERENCE".into()),
                date_label: Some("26-NOV-2026".into()),
            },
        },
    ] {
        let context = sequence_export::revision(&source, &format).unwrap();
        let reviewed = run(
            AnalysisJob::ReviewSequenceExport {
                molecule: source.clone(),
                format: format.clone(),
            },
            context.clone(),
            174,
            164,
        );
        let AnalysisResult::SequenceExportAssessment(review) = reviewed.result else {
            unreachable!()
        };
        review.validate_binding(&context).unwrap();
        assert_eq!(review.source_revision, source.revision());
        if !review.losses.is_empty() {
            assert!(sequence_export::prepare(&source, &format, &context, false).is_err());
        }
        let result = run(
            AnalysisJob::PrepareSequenceExport {
                molecule: source.clone(),
                format: format.clone(),
                expected_context_revision: context.clone(),
                allow_omissions: true,
            },
            context.clone(),
            175,
            165,
        );
        let AnalysisResult::SequenceExportPrepared(prepared) = result.result else {
            unreachable!()
        };
        prepared.validate_binding(&context).unwrap();
        assert_eq!(prepared.assessment.source_revision, source.revision());
        match &format {
            Format::SnapGene => assert_eq!(&prepared.bytes[5..13], b"SnapGene"),
            Format::Xdna => {
                assert_eq!(&prepared.bytes[..3], &[0, 1, 1]);
                assert_eq!(&prepared.bytes[112..120], b"ACGTCGAT");
            }
            Format::Gcg { .. } => {
                let text = std::str::from_utf8(&prepared.bytes).unwrap();
                assert!(text.contains("Length: 8"));
                assert!(text.contains("Check:"));
            }
        }
        assert!(sequence_export::prepare(&source, &format, &"0".repeat(64), true).is_err());
    }
    assert_eq!(source, before);
}
fn source_trace() -> traces::Trace {
    // SCF2.02: four interleaved sample rows, two base records, actual per-channel accuracies.
    let mut b = vec![0u8; 181];
    let word = |b: &mut [u8], at: usize, n: u32| b[at..at + 4].copy_from_slice(&n.to_be_bytes());
    b[..4].copy_from_slice(b".scf");
    for (at, n) in [
        (4, 4),
        (8, 128),
        (12, 2),
        (24, 144),
        (28, 13),
        (32, 168),
        (40, 1),
        (44, 2),
    ] {
        word(&mut b, at, n);
    }
    b[36..40].copy_from_slice(b"2.02");
    b[128..144].copy_from_slice(&[0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23, 30, 31, 32, 33]);
    word(&mut b, 144, 1);
    b[148..153].copy_from_slice(&[1, 2, 3, 4, b'A']);
    word(&mut b, 156, 3);
    b[160..165].copy_from_slice(&[5, 6, 7, 8, b'C']);
    b[168..].copy_from_slice(b"NAME=actual\n\0");
    traces::import_trace(&b).unwrap()
}
#[test]
fn trace_worker_pins_view_precision_and_actual_signal_and_refuses_fabrication() {
    let trace = source_trace();
    let before = trace.clone();
    let options = scf_export::Options {
        precision: scf_export::SamplePrecision::SixteenBit,
        loss_policy: scf_export::LossPolicy::AllowDocumentaryOmissions,
    };
    let context = scf_export::revision(&trace, &options).unwrap();
    let response = run(
        AnalysisJob::ExportScfTrace {
            trace: trace.clone(),
            options: options.clone(),
        },
        context.clone(),
        176,
        166,
    );
    let AnalysisResult::ScfTraceExported(report) = response.result else {
        unreachable!()
    };
    report.validate_binding(&context).unwrap();
    assert_eq!(report.source_hash, trace.source_hash);
    assert_eq!(report.sample_size, 2);
    assert_eq!(&report.bytes[36..40], b"3.00");
    assert_eq!(
        traces::import_trace(&report.bytes).unwrap().channels,
        trace.channels
    );
    assert!(report.validate_binding(&"0".repeat(64)).is_err());
    let mut changed = trace.clone();
    changed.name.push_str(" changed");
    assert_ne!(scf_export::revision(&changed, &options).unwrap(), context);
    changed = trace.clone();
    changed.channels.a[0] = 99;
    assert!(scf_export::export(&changed, &options).is_err());
    let mut other = options.clone();
    other.precision = scf_export::SamplePrecision::EightBit;
    assert_ne!(scf_export::revision(&trace, &other).unwrap(), context);
    assert_eq!(trace, before);
}
