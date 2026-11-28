use collagene_core::{
    formats::scf,
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    review::reads::traces::{self, display::PlotData},
};
use collagene_model::Sequence;
fn source() -> traces::Trace {
    let comments = b"NAME=source\n\0";
    let mut bytes = vec![0u8; 216 + comments.len()];
    bytes[..4].copy_from_slice(b".scf");
    for (at, value) in [
        (4, 8u32),
        (8, 128),
        (12, 2),
        (24, 192),
        (28, comments.len() as u32),
        (32, 216),
        (40, 2),
        (44, 2),
    ] {
        bytes[at..at + 4].copy_from_slice(&value.to_be_bytes());
    }
    bytes[36..40].copy_from_slice(b"2.02");
    for sample in 0..8 {
        for channel in 0..4 {
            let at = 128 + (sample * 4 + channel) * 2;
            bytes[at..at + 2].copy_from_slice(&((sample * 8 + channel) as u16).to_be_bytes());
        }
    }
    for (i, peak, call) in [(0, 1u32, b'A'), (1, 5, b'C')] {
        let at = 192 + i * 12;
        bytes[at..at + 4].copy_from_slice(&peak.to_be_bytes());
        bytes[at + 4..at + 8].copy_from_slice(&[
            1 + i as u8 * 4,
            2 + i as u8 * 4,
            3 + i as u8 * 4,
            4 + i as u8 * 4,
        ]);
        bytes[at + 8] = call;
    }
    bytes[216..].copy_from_slice(comments);
    traces::import_trace(&bytes).unwrap()
}
#[test]
fn stateless_and_session_display_preserve_original_evidence_and_oriented_samples() {
    let trace = source();
    for input in [trace.clone(), trace.reverse_complement().unwrap()] {
        let molecule = input.attached_molecule().unwrap();
        let revision = molecule.revision();
        let request = Request {
            id: 29,
            revision: revision.clone(),
            job: AnalysisJob::PrepareTraceDisplay {
                molecule: molecule.clone(),
            },
        };
        let bytes = request.encode().unwrap();
        let mut session = protocol::session::AnalysisSession::default();
        for reply in [
            protocol::execution::process_request(&bytes).unwrap(),
            session.process_request(&bytes).unwrap(),
        ] {
            let response = Response::decode(&reply).unwrap();
            assert_eq!(response.id, 29);
            let AnalysisResult::TraceDisplayPrepared(data) = response.result else {
                panic!("display result")
            };
            data.validate_binding(&revision).unwrap();
            assert_eq!(data.trace.channels, input.channels);
            assert_eq!(data.trace.sequence, input.sequence);
            assert_eq!(data.trace.qualities, input.qualities);
            assert_eq!(data.trace.source_hash, trace.source_hash);
            assert_eq!(data.minimum, 0);
            assert_eq!(data.maximum, 59);
            let json = serde_json::to_value(&data).unwrap();
            assert!(json["trace"].get("original_bytes").is_none());
            assert!(json["trace"].get("metadata").is_none());
            match input.orientation {
                traces::TraceOrientation::AsStored => {
                    assert_eq!(data.trace.peak_positions, [1, 5]);
                    assert_eq!(data.trace.channels.a, [0, 8, 16, 24, 32, 40, 48, 56]);
                }
                traces::TraceOrientation::ReverseComplement => {
                    assert_eq!(data.trace.peak_positions, [2, 6]);
                    assert_eq!(data.trace.channels.a, [59, 51, 43, 35, 27, 19, 11, 3]);
                }
            }
        }
        assert_eq!(
            traces::original_trace_bytes(&molecule).unwrap(),
            trace.original_bytes
        );
    }
}
#[test]
fn stale_sequence_wrong_source_and_malformed_display_do_not_reuse_old_calls() {
    let molecule = source().attached_molecule().unwrap();
    let revision = molecule.revision();
    let data = PlotData::prepare(&molecule, &mut || false).unwrap();
    let mut changed = molecule.clone();
    changed.sequence = Sequence::dna("AT").unwrap();
    let stale = Request {
        id: 1,
        revision: revision.clone(),
        job: AnalysisJob::PrepareTraceDisplay {
            molecule: changed.clone(),
        },
    }
    .encode()
    .unwrap();
    assert!(protocol::execution::process_request(&stale).is_err());
    let current = Request {
        id: 1,
        revision: changed.revision(),
        job: AnalysisJob::PrepareTraceDisplay {
            molecule: changed.clone(),
        },
    }
    .encode()
    .unwrap();
    assert!(matches!(
        Response::decode(&protocol::execution::process_request(&current).unwrap())
            .unwrap()
            .result,
        AnalysisResult::Error(_)
    ));
    assert_eq!(
        traces::original_trace_bytes(&changed).unwrap(),
        source().original_bytes
    );
    let mut bad = data.clone();
    bad.trace.peak_positions = [7, 2].into();
    assert!(bad.validate_binding(&revision).is_err());
    let mut bad = data.clone();
    bad.trace.channels.c.pop();
    assert!(bad.validate_binding(&revision).is_err());
    let mut bad = data.clone();
    bad.minimum = -1;
    assert!(bad.validate_binding(&revision).is_err());
    let mut bad = data.clone();
    bad.trace.called_symbols[0] = b'T';
    assert!(bad.validate_binding(&revision).is_err());
    assert!(data.validate_binding(&"0".repeat(64)).is_err());
    assert!(PlotData::prepare(&molecule, &mut || true).is_err());
}
#[test]
fn regenerated_export_reparses_document_instead_of_accepting_display_arrays() {
    let trace = source();
    let molecule = trace.attached_molecule().unwrap();
    let mut display = PlotData::prepare(&molecule, &mut || false).unwrap();
    display.trace.channels.a[0] = 999;
    let options = scf::Options {
        precision: scf::SamplePrecision::Smallest,
        loss_policy: scf::LossPolicy::AllowDocumentaryOmissions,
    };
    let context = scf::document_revision(&molecule, &options).unwrap();
    let request = Request {
        id: 3,
        revision: context.clone(),
        job: AnalysisJob::ExportScfTrace {
            molecule: molecule.clone(),
            options,
        },
    }
    .encode()
    .unwrap();
    let AnalysisResult::ScfTraceExported(export) =
        Response::decode(&protocol::execution::process_request(&request).unwrap())
            .unwrap()
            .result
    else {
        panic!("SCF export")
    };
    export.validate_binding(&context).unwrap();
    assert_eq!(export.source_revision, molecule.revision());
    assert_eq!(
        traces::import_trace(&export.report.bytes).unwrap().channels,
        trace.channels
    );
    assert_ne!(trace.channels, display.trace.channels);
}

#[test]
fn trace_document_requests_do_not_gain_the_large_reply_budget() {
    let mut source = source().attached_molecule().unwrap();
    source.metadata.insert(
        "oversized_untrusted_field".into(),
        "x".repeat(collagene_core::formats::sequence::MAX_IMPORT_BYTES),
    );
    let request = Request {
        id: 4,
        revision: "a".repeat(64),
        job: AnalysisJob::PrepareTraceDisplay { molecule: source },
    };
    assert!(request.encode().unwrap_err().to_string().contains("32 MiB"));
    // A malicious current-version frame cannot grant its ordinary selector the
    // larger response budget by constructing a valid length/checksum itself.
    let payload = postcard::to_allocvec(&request).unwrap();
    let mut frame = b"CGQ\0".to_vec();
    frame.extend_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    frame.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    frame.extend_from_slice(blake3::hash(&payload).as_bytes());
    frame.extend(payload);
    assert!(matches!(
        protocol::execution::process_request(&frame),
        Err(collagene_model::Error::TooLarge)
    ));
    assert!(matches!(
        protocol::session::AnalysisSession::default().process_request(&frame),
        Err(collagene_model::Error::TooLarge)
    ));
}
