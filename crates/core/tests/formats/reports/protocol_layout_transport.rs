#[path = "../../common/protocol_layout.rs"]
mod fixture;
use collagene_core::{
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    protocol_layout as layout,
};
fn execute(job: AnalysisJob, revision: String) -> AnalysisResult {
    let q = Request {
        id: 17,
        revision: revision.clone(),
        job,
    };
    let bytes = protocol::request::Request::encode(&q).unwrap();
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        protocol::schema::WORKER_SCHEMA_VERSION
    );
    let raw = protocol::execution::process_request(&bytes).unwrap();
    let first = protocol::response::Response::decode(&raw).unwrap();
    let raw = protocol::session::AnalysisSession::default()
        .process_request(&bytes)
        .unwrap();
    let second = protocol::response::Response::decode(&raw).unwrap();
    assert_eq!(first.id, 17);
    assert_eq!(second.id, 17);
    assert_eq!(first.revision, revision);
    assert_eq!(second.revision, revision);
    assert_eq!(
        serde_json::to_value(&first.result).unwrap(),
        serde_json::to_value(&second.result).unwrap()
    );
    second.result
}
#[test]
fn actual_layout_preview_capture_and_three_documentary_exports_bind_exact_sources() {
    let (workspace, request) = fixture::fixture();
    let before = serde_json::to_vec(&workspace).unwrap();
    let context = layout::preview_revision(&workspace, &request).unwrap();
    assert!(matches!(
        execute(
            AnalysisJob::PreviewProtocolLayout {
                workspace: workspace.clone(),
                request: request.clone()
            },
            context.clone()
        ),
        AnalysisResult::ProtocolLayoutPreview(_)
    ));
    let AnalysisResult::ProtocolLayoutCaptured(c) = execute(
        AnalysisJob::CaptureProtocolLayout {
            workspace: workspace.clone(),
            request: request.clone(),
        },
        context.clone(),
    ) else {
        panic!("expected capture")
    };
    c.validate_binding(&context).unwrap();
    assert_eq!(
        serde_json::to_value(&c.workspace.documents).unwrap(),
        serde_json::to_value(&workspace.documents).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&c.workspace.history).unwrap(),
        serde_json::to_value(&workspace.history).unwrap()
    );
    for format in [
        layout::ExportFormat::Csv,
        layout::ExportFormat::Html,
        layout::ExportFormat::Pdf,
    ] {
        let revision = layout::export_revision(&c.workspace, &c.attachment_hash, format).unwrap();
        let AnalysisResult::ProtocolLayoutExported(e) = execute(
            AnalysisJob::ExportProtocolLayout {
                workspace: c.workspace.clone(),
                layout_hash: c.attachment_hash.clone(),
                format,
            },
            revision.clone(),
        ) else {
            panic!("expected export")
        };
        e.validate_binding(&revision).unwrap();
        assert!(!e.bytes.is_empty());
    }
    assert_eq!(serde_json::to_vec(&workspace).unwrap(), before);
}
#[test]
fn changed_assignment_context_and_missing_dependencies_do_not_capture() {
    let (w, r) = fixture::fixture();
    let context = layout::preview_revision(&w, &r).unwrap();
    let mut changed = r;
    changed.assignments[0].well.column = 1;
    let q = Request {
        id: 2,
        revision: context,
        job: AnalysisJob::CaptureProtocolLayout {
            workspace: w.clone(),
            request: changed,
        },
    };
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&q).unwrap())
            .is_err()
    );
    let mut missing = w;
    missing.attachments.clear();
    assert!(layout::preview(&missing, &fixture::fixture().1).is_err());
}
