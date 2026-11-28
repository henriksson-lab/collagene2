use collagene_core::{
    documents::workspace::Workspace,
    formats::oligos::{self as oligo_formats, QualityEncoding},
    material::lineage::EventBinding,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    review::reads::{
        checks as read_checks,
        trace_validation::{TraceValidationOptions, TrimMode, ValidationStatus},
    },
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> Request {
    let reference = Molecule::dna("construct", "TTTCGATACGGG", Topology::Linear).unwrap();
    let read = oligo_formats::import_fastq(b"@read\nCGATAC\n+\nIIIIII\n", QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    let workspace = Workspace::from_documents("sequencing", vec![reference, read], 0).unwrap();
    let reference = binding(&workspace, 0);
    let reads = vec![binding(&workspace, 1)];
    let mut options = TraceValidationOptions::default();
    options.trim.mode = TrimMode::AllCalls;
    options.trim.min_retained_bases = 1;
    let revision = read_checks::capture_revision(&workspace, &reference, &reads, &options).unwrap();
    Request {
        id: 420,
        revision,
        job: AnalysisJob::CaptureReadCheck {
            workspace,
            reference,
            reads,
            options,
        },
    }
}
#[test]
fn captured_check_and_fresh_historical_replay_have_exact_stable_tags_and_reject_old_headers() {
    let request = fixture();
    let bytes = protocol::request::Request::encode(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 103);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&41u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 95);
    let AnalysisResult::ReadCheckCaptured(commit) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected retained evidence")
    };
    commit.validate_binding(&request.revision).unwrap();
    let mut old = bytes;
    old[4..6].copy_from_slice(&41u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let revision =
        read_checks::replay_revision(&commit.workspace, &commit.attachment_hash).unwrap();
    let replay = Request {
        id: 421,
        revision: revision.clone(),
        job: AnalysisJob::ReplayReadCheck {
            workspace: commit.workspace,
            attachment_hash: commit.attachment_hash,
        },
    };
    let bytes = protocol::request::Request::encode(&replay).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 104);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&41u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 96);
    let AnalysisResult::ReadCheckReplayed(checked) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected fresh reconstruction")
    };
    checked.validate_binding(&revision).unwrap();
    assert_eq!(checked.report.status, ValidationStatus::Incomplete);
    assert_eq!(checked.report.eligible_coverage, 6);
    assert_eq!(checked.report.consensus, "NNNCGATACNNN");
}
#[test]
fn a_different_reference_event_cannot_reuse_a_pending_capture_token() {
    let mut request = fixture();
    if let AnalysisJob::CaptureReadCheck { reference, .. } = &mut request.job {
        reference.event_id = format!("event-{}", "0".repeat(64));
    }
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap()
        )
        .is_err()
    );
}
