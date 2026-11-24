//! Native/WASM worker boundary contracts for durable workspace I/O and trace validation.
use collagene_core::{
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    trace_validation::{TraceValidationOptions, TrimMode},
    traces,
    workspace::*,
    *,
};
use std::collections::BTreeMap;
fn source() -> Molecule {
    Molecule::dna("reference", "ACGT", Topology::Linear).unwrap()
}
fn workspace() -> Workspace {
    Workspace::from_documents("worker workspace", vec![source()], 0).unwrap()
}
fn execute(request: &Request) -> Response {
    let bytes = protocol::encode_request(request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    response
}
fn envelope(version: u16, magic: &[u8; 4], payload: &[u8]) -> Vec<u8> {
    let mut bytes = magic.to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(payload).as_bytes());
    bytes.extend(payload);
    bytes
}
#[test]
fn workspace_save_load_jobs_preserve_document_identity_history_and_request_binding() {
    let mut workspace = workspace();
    let old = workspace.documents[0].head.clone();
    let edited = workspace
        .head("document-1")
        .unwrap()
        .replace(1, 2, "T")
        .unwrap();
    let event = RevisionEvent::new(
        "document-1",
        edited.revision(),
        vec![workspace.history[0].id.clone()],
        "Substitution",
    )
    .unwrap();
    workspace.documents[0].undo.push(old);
    workspace.documents[0]
        .undo_events
        .push(workspace.history[0].id.clone());
    workspace.documents[0].head = edited.revision();
    workspace.documents[0].head_event = Some(event.id.clone());
    workspace.history.push(event);
    workspace.revisions.insert(edited.revision(), edited);
    let request = Request {
        id: 0xfff123,
        revision: workspace
            .head(&workspace.active_document)
            .unwrap()
            .revision(),
        job: AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
    };
    let AnalysisResult::WorkspaceSaved(bytes) = execute(&request).result else {
        panic!("wrong save result")
    };
    assert!(is_workspace_archive(&bytes));
    let request = Request {
        id: 91,
        revision: blake3::hash(&bytes).to_hex().to_string(),
        job: AnalysisJob::LoadWorkspace { bytes },
    };
    let AnalysisResult::WorkspaceLoaded(restored) = execute(&request).result else {
        panic!("wrong load result")
    };
    assert_eq!(
        serde_json::to_value(&restored).unwrap(),
        serde_json::to_value(&workspace).unwrap()
    );
}
#[test]
fn workspace_stale_bindings_corruption_and_malformed_data_do_not_become_loaded_state() {
    let workspace = workspace();
    let request = Request {
        id: 1,
        revision: "stale".into(),
        job: AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
    };
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
    let bytes = save_workspace(&workspace).unwrap();
    let request = Request {
        id: 2,
        revision: "stale".into(),
        job: AnalysisJob::LoadWorkspace { bytes },
    };
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
    let bytes = b"not a workspace".to_vec();
    let request = Request {
        id: 3,
        revision: blake3::hash(&bytes).to_hex().to_string(),
        job: AnalysisJob::LoadWorkspace { bytes },
    };
    assert!(matches!(execute(&request).result, AnalysisResult::Error(_)));
    let request = Request {
        id: 4,
        revision: source().revision(),
        job: AnalysisJob::SaveWorkspace { workspace },
    };
    let mut bytes = protocol::encode_request(&request).unwrap();
    *bytes.last_mut().unwrap() ^= 1;
    assert!(protocol::process_request(&bytes).is_err());
}
#[test]
fn new_workspace_jobs_results_and_trace_jobs_cannot_claim_old_schemas() {
    let workspace = workspace();
    let request = Request {
        id: 5,
        revision: source().revision(),
        job: AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
    };
    for version in 1u16..=5 {
        let mut bytes = protocol::encode_request(&request).unwrap();
        bytes[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(matches!(
            protocol::process_request(&bytes),
            Err(Error::Unsupported(_))
        ));
    }
    let result = Response {
        id: 5,
        revision: source().revision(),
        result: AnalysisResult::WorkspaceLoaded(workspace.clone()),
    };
    for version in 1u16..=5 {
        assert!(
            protocol::decode_response(&envelope(
                version,
                b"CGR\0",
                &postcard::to_allocvec(&result).unwrap()
            ))
            .is_err()
        );
    }
    let mut invalid = workspace;
    invalid.active_document = "missing".into();
    let response = Response {
        id: 6,
        revision: source().revision(),
        result: AnalysisResult::WorkspaceLoaded(invalid),
    };
    assert!(
        protocol::decode_response(&envelope(
            6,
            b"CGR\0",
            &postcard::to_allocvec(&response).unwrap()
        ))
        .is_err()
    );
    // An original variant remains readable when honestly tagged as the prior schema.
    let request = Request {
        id: 7,
        revision: source().revision(),
        job: AnalysisJob::Search {
            molecule: source(),
            query: "CG".into(),
        },
    };
    let mut bytes = protocol::encode_request(&request).unwrap();
    bytes[4..6].copy_from_slice(&5u16.to_le_bytes());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert!(matches!(response.result, AnalysisResult::Search(_)));
}
#[test]
fn large_workspace_roundtrip_has_its_own_cap_without_enlarging_document_or_scientific_jobs() {
    let mut workspace = workspace();
    for value in [7u8, 9] {
        let bytes = vec![value; 16 * 1024 * 1024 + 1];
        let digest = blake3::hash(&bytes).to_hex().to_string();
        workspace.attachments.insert(
            digest,
            WorkspaceAttachment {
                kind: "evidence".into(),
                name: format!("binary-{value}"),
                media_type: "application/octet-stream".into(),
                bytes,
            },
        );
    }
    let request = Request {
        id: 0xffff,
        revision: source().revision(),
        job: AnalysisJob::SaveWorkspace { workspace },
    };
    let encoded = protocol::encode_request(&request).unwrap();
    assert!(encoded.len() > 32 * 1024 * 1024);
    let raw = protocol::process_request(&encoded).unwrap();
    assert!(raw.len() > 32 * 1024 * 1024);
    let AnalysisResult::WorkspaceSaved(bytes) = protocol::decode_response(&raw).unwrap().result
    else {
        panic!("wrong result")
    };
    let request = Request {
        id: 12,
        revision: blake3::hash(&bytes).to_hex().to_string(),
        job: AnalysisJob::LoadWorkspace { bytes },
    };
    let AnalysisResult::WorkspaceLoaded(restored) = execute(&request).result else {
        panic!("wrong result")
    };
    assert_eq!(restored.attachments.len(), 2);
    assert_eq!(
        restored
            .attachments
            .values()
            .map(|a| a.bytes.len())
            .sum::<usize>(),
        32 * 1024 * 1024 + 2
    );
    let mut ordinary = source();
    ordinary.metadata = BTreeMap::from([("large".into(), "x".repeat(32 * 1024 * 1024))]);
    assert!(matches!(
        protocol::save_document(&ordinary),
        Err(Error::TooLarge)
    ));
    let request = Request {
        id: 13,
        revision: ordinary.revision(),
        job: AnalysisJob::Search {
            molecule: ordinary,
            query: "A".into(),
        },
    };
    assert!(matches!(
        protocol::encode_request(&request),
        Err(Error::TooLarge)
    ));
}
#[test]
fn oversized_nonworkspace_and_old_large_headers_fail_before_payload_deserialization() {
    // Independently constructed prefix: varint id=1, revision UTF8 length=1/'x', enum tag.
    // Remaining zeros are intentionally not a valid Search/Workspace payload.
    for (magic, tag) in [(b"CGQ\0", 0u8), (b"CGR\0", 2u8)] {
        let mut payload = vec![1, 1, b'x', tag];
        payload.resize(32 * 1024 * 1024 + 1, 0);
        assert!(matches!(
            if magic == b"CGQ\0" {
                protocol::process_request(&envelope(6, magic, &payload)).map(|_| ())
            } else {
                protocol::decode_response(&envelope(6, magic, &payload)).map(|_| ())
            },
            Err(Error::TooLarge)
        ));
    }
    let mut payload = vec![1, 1, b'x', 19];
    payload.resize(32 * 1024 * 1024 + 1, 0);
    assert!(matches!(
        protocol::process_request(&envelope(5, b"CGQ\0", &payload)),
        Err(Error::TooLarge)
    ));
    let mut short = b"CGQ\0".to_vec();
    short.extend(6u16.to_le_bytes());
    short.extend(((128 * 1024 * 1024) + 1u32).to_le_bytes());
    short.extend([0u8; 32]);
    assert!(matches!(
        protocol::process_request(&short),
        Err(Error::TooLarge)
    ));
}
fn trace_fixture() -> Molecule {
    let mut bytes = vec![0; 128];
    bytes[..4].copy_from_slice(b".scf");
    for (offset, value) in [(4, 4u32), (8, 128), (12, 4), (24, 160), (40, 2), (44, 2)] {
        bytes[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
    }
    bytes[36..40].copy_from_slice(b"2.02");
    for sample in 0..4 {
        for channel in 0..4 {
            bytes.extend((if sample == channel { 80u16 } else { 0 }).to_be_bytes());
        }
    }
    for (i, base) in b"ACGT".iter().enumerate() {
        bytes.extend((i as u32).to_be_bytes());
        for channel in 0..4 {
            bytes.push(if channel == i { 40 } else { 0 });
        }
        bytes.push(*base);
        bytes.extend([0; 3]);
    }
    traces::import_trace(&bytes)
        .unwrap()
        .attached_molecule()
        .unwrap()
}
#[test]
fn trace_validation_worker_binds_reference_and_preserves_raw_quality_options() {
    let options = TraceValidationOptions {
        both_strands: false,
        trim: collagene_core::trace_validation::TrimOptions {
            mode: TrimMode::AllCalls,
            min_retained_bases: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let read = trace_fixture();
    let request = Request {
        id: 0x30,
        revision: source().revision(),
        job: AnalysisJob::ValidateTraces {
            molecule: source(),
            reads: vec![read.clone()],
            options: options.clone(),
        },
    };
    let AnalysisResult::TraceValidation(report) = execute(&request).result else {
        panic!("wrong trace result")
    };
    assert_eq!(report.reference_revision, source().revision());
    assert_eq!(report.read_revisions, vec![read.revision()]);
    assert_eq!(report.options, options);
    assert_eq!(report.consensus, "ACGT");
    assert_eq!(report.eligible_coverage, 4);
    let mut bytes = protocol::encode_request(&request).unwrap();
    bytes[4..6].copy_from_slice(&5u16.to_le_bytes());
    assert!(matches!(
        protocol::process_request(&bytes),
        Err(Error::Unsupported(_))
    ));
}
