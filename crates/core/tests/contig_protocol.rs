//! Source-batch binding and append-only schema gates for de novo assembly.
use collagene_core::{
    Molecule, Topology,
    contig::*,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
fn fixture() -> (Vec<ReadInput>, ContigOptions) {
    let reads = [
        ("left", "AGTCGACTA"),
        ("middle", "ACTAGGCTT"),
        ("right", "GCTTACG"),
    ]
    .into_iter()
    .map(|(id, bases)| ReadInput {
        id: id.into(),
        molecule: Molecule::dna(id, bases, Topology::Linear).unwrap(),
    })
    .collect();
    let mut options = ContigOptions::default();
    options.trim.min_retained_bases = 1;
    options.min_overlap = 4;
    options.max_error_bp = 0;
    options.both_strands = false;
    (reads, options)
}
fn raw(magic: &[u8; 4], version: u16, payload: &[u8]) -> Vec<u8> {
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(payload).as_bytes());
    out.extend(payload);
    out
}
#[test]
fn assemble_contigs_is_bound_to_ordered_read_roles_revisions_and_all_options() {
    let (reads, options) = fixture();
    let revision = batch_revision(&reads, &options).unwrap();
    let request = Request {
        id: 301,
        revision: revision.clone(),
        job: AnalysisJob::AssembleContigs {
            reads: reads.clone(),
            options: options.clone(),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(response.id, 301);
    assert_eq!(response.revision, revision);
    let AnalysisResult::Contigs(report) = response.result else {
        panic!("wrong response")
    };
    report.validate_bindings_for(&reads, &options).unwrap();
    assert!(
        report
            .contigs
            .iter()
            .any(|c| c.sequence == "AGTCGACTAGGCTTACG"),
        "{report:?}"
    );
    let mut altered = reads.clone();
    altered[0].id = "different-role".into();
    assert_ne!(batch_revision(&altered, &options).unwrap(), revision);
    assert!(report.validate_bindings_for(&altered, &options).is_err());
    let mut reordered = reads.clone();
    reordered.swap(0, 1);
    assert_ne!(batch_revision(&reordered, &options).unwrap(), revision);
    assert!(report.validate_bindings_for(&reordered, &options).is_err());
    let mut wrong_source = reads.clone();
    wrong_source[0].molecule = Molecule::dna("left", "AGTCGACTT", Topology::Linear).unwrap();
    assert_ne!(batch_revision(&wrong_source, &options).unwrap(), revision);
    assert!(
        report
            .validate_bindings_for(&wrong_source, &options)
            .is_err()
    );
    let mut changed = options.clone();
    changed.min_support += 1;
    assert!(report.validate_bindings_for(&reads, &changed).is_err());
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let wrong = Response {
        id: 301,
        revision: "0".repeat(64),
        result: AnalysisResult::Contigs(report),
    };
    let payload = postcard::to_allocvec(&wrong).unwrap();
    assert!(
        protocol::decode_response(&raw(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &payload))
            .is_err()
    );
}
#[test]
fn contig_job_and_result_require_schema_fifteen_before_body_deserialization() {
    for (magic, tag) in [(b"CGQ\0", 47u32), (b"CGR\0", 45u32)] {
        let prefix = postcard::to_allocvec(&(302u32, "source-batch", tag)).unwrap();
        for version in 1..15 {
            let bytes = raw(magic, version, &prefix);
            let error = if magic == b"CGQ\0" {
                protocol::process_request(&bytes).unwrap_err()
            } else {
                protocol::decode_response(&bytes).unwrap_err()
            };
            assert!(
                matches!(error, collagene_core::Error::Unsupported(_)),
                "{error}"
            );
        }
        let mut large = prefix;
        large.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
        let bytes = raw(magic, protocol::WORKER_SCHEMA_VERSION, &large);
        let error = if magic == b"CGQ\0" {
            protocol::process_request(&bytes).unwrap_err()
        } else {
            protocol::decode_response(&bytes).unwrap_err()
        };
        assert!(matches!(error, collagene_core::Error::TooLarge));
    }
}
#[test]
fn assembly_source_bounds_precede_batch_hashing_and_stale_revision_comparison() {
    let (mut reads, options) = fixture();
    reads[0].molecule.metadata.insert(
        "oversized source evidence".into(),
        "x".repeat(20 * 1024 * 1024),
    );
    let request = Request {
        id: 303,
        revision: "0".repeat(64),
        job: AnalysisJob::AssembleContigs { reads, options },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert!(matches!(
        protocol::process_request(&bytes),
        Err(collagene_core::Error::TooLarge)
    ));
}

#[test]
fn recognized_assembly_evidence_requires_schema_fifteen_in_inline_workspaces() {
    use collagene_core::workspace::{Workspace, WorkspaceAttachment};
    let molecule = Molecule::dna("source", "ACGT", Topology::Linear).unwrap();
    let mut workspace =
        Workspace::from_documents("assembly evidence", vec![molecule.clone()], 0).unwrap();
    let evidence = b"{}".to_vec();
    let key = blake3::hash(&evidence).to_hex().to_string();
    workspace.attachments.insert(
        key,
        WorkspaceAttachment {
            kind: "read-assembly-v1".into(),
            name: "untrusted producer record".into(),
            media_type: "application/json".into(),
            bytes: evidence,
        },
    );
    // The minimum version must be checked before semantic interpretation of a
    // new attachment, even when the producer's record would later be invalid.
    let request = Request {
        id: 304,
        revision: molecule.revision(),
        job: AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
    };
    let bytes = raw(b"CGQ\0", 14, &postcard::to_allocvec(&request).unwrap());
    let error = protocol::process_request(&bytes).unwrap_err();
    assert!(
        matches!(error, collagene_core::Error::Unsupported(_)),
        "{error}"
    );
    let response = Response {
        id: 304,
        revision: molecule.revision(),
        result: AnalysisResult::WorkspaceLoaded(workspace),
    };
    let bytes = raw(b"CGR\0", 14, &postcard::to_allocvec(&response).unwrap());
    let error = protocol::decode_response(&bytes).unwrap_err();
    assert!(
        matches!(error, collagene_core::Error::Unsupported(_)),
        "{error}"
    );
}
