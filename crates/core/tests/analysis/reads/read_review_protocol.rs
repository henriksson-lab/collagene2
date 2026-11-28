use collagene_core::{
    oligo_formats::{self, QualityEncoding},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    read_checks,
    read_review::{self, Choice, Decision, EvidenceKey},
    run_lineage::EventBinding,
    trace_validation::{TraceValidationOptions, TrimMode},
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};
fn bind(w: &Workspace, n: usize) -> EventBinding {
    let d = &w.documents[n];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, read_review::Request) {
    let reference = Molecule::dna("reference", "TTTCGATACGGG", Topology::Linear).unwrap();
    let read = oligo_formats::import_fastq(
        b"@observed\nTTTCGTTACGGG\n+\nIIIIIIIIIIII\n",
        QualityEncoding::Phred33,
    )
    .unwrap()
    .molecules
    .remove(0);
    let w = Workspace::from_documents("review", vec![reference, read], 0).unwrap();
    let mut options = TraceValidationOptions::default();
    options.trim.mode = TrimMode::AllCalls;
    options.trim.min_retained_bases = 1;
    let capture = read_checks::capture(&w, bind(&w, 0), vec![bind(&w, 1)], &options).unwrap();
    (
        capture.workspace,
        read_review::Request {
            attachment_hash: capture.attachment_hash,
            name: "proposed reference".into(),
            decisions: vec![Decision {
                key: EvidenceKey::Site { position: 5 },
                choice: Choice::UseBase {
                    read_index: 0,
                    call_index: 5,
                    base: 'T',
                },
                reason: "Explicit observed T; recheck derived reference independently.".into(),
            }],
        },
    )
}
#[test]
fn review_preview_and_commit_use_exact_tags_and_reject_pre_43_headers() {
    let (workspace, request) = fixture();
    let revision = read_review::preview_revision(&workspace, &request).unwrap();
    let bytes = protocol::request::Request::encode(&Request {
        id: 430,
        revision: revision.clone(),
        job: AnalysisJob::PreviewReadReview {
            workspace: workspace.clone(),
            request: request.clone(),
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 105);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&42u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 97);
    let AnalysisResult::ReadReviewPreview(preview) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected explicit proposal")
    };
    preview.validate_binding(&revision).unwrap();
    assert_eq!(preview.target.unwrap().sequence.as_str(), "TTTCGTTACGGG");
    let mut old = bytes;
    old[4..6].copy_from_slice(&42u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let revision = read_review::commit_revision(&workspace, &request).unwrap();
    let bytes = protocol::request::Request::encode(&Request {
        id: 431,
        revision: revision.clone(),
        job: AnalysisJob::CommitReadReview {
            workspace: workspace.clone(),
            request,
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 106);
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 98);
    let AnalysisResult::ReadReviewCommitted(commit) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected immutable decision record")
    };
    commit.validate_binding(&revision).unwrap();
    let target = commit.target.unwrap();
    assert_eq!(
        commit
            .workspace
            .head(&target.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "TTTCGTTACGGG"
    );
    assert_eq!(
        commit.workspace.documents[0].head,
        workspace.documents[0].head
    );
}
#[test]
fn pending_review_cannot_be_reused_for_a_different_expected_allele() {
    let (workspace, mut request) = fixture();
    let revision = read_review::preview_revision(&workspace, &request).unwrap();
    if let Choice::UseBase { base, .. } = &mut request.decisions[0].choice {
        *base = 'A';
    }
    let bytes = protocol::request::Request::encode(&Request {
        id: 432,
        revision,
        job: AnalysisJob::PreviewReadReview { workspace, request },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&bytes).is_err());
}
