//! Fresh private bulk finalization and output export across the binary boundary.
use collagene_core::{
    chemistry::Terminus,
    cloning::assembly::recipe::Preparation,
    documents::workspace::Workspace,
    material::{
        lineage::{self as run_lineage, EventBinding},
        outputs::Options as OutputOptions,
    },
    planning::combinatorial::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
        batch::{self, RowRequest},
        export,
        release::{self, Selection},
    },
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
        session::AnalysisSession,
    },
};
use collagene_model::{Error, Molecule, Topology};
fn binding(workspace: &Workspace, index: usize) -> EventBinding {
    let doc = &workspace.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn fixture() -> (Workspace, String) {
    let frame = Workspace::from_documents(
        "bulk output sources",
        vec![
            Molecule::dna("backbone", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("insert A", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("insert B", "CCAA", Topology::Linear).unwrap(),
            Molecule::dna("independent wrong target", "GGGG", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let choice = |id: &str, index| Choice {
        generation: None,
        id: id.into(),
        label: id.into(),
        source: binding(&frame, index),
        selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            start: 0,
            end: 4,
            reverse: false,
        },
        preparation_reverse: false,
        preparation: Some(Preparation::DeclaredDuplex {
            left: Terminus::CLEAVED,
            right: Terminus::CLEAVED,
        }),
        treatment: None,
    };
    let design = Design {
        material_bindings: vec![],
        id: "wire-release".into(),
        name: "Release fixture".into(),
        circular: false,
        bins: vec![
            Bin {
                id: "backbone".into(),
                label: "Backbone".into(),
                choices: vec![choice("backbone-a", 0)],
            },
            Bin {
                id: "insert".into(),
                label: "Insert".into(),
                choices: vec![choice("insert-a", 1), choice("insert-b", 2)],
            },
        ],
        method: Method::Ligation {
            allow_nicked: false,
        },
    };
    let row = |id: &str, insert: &str, target| {
        RowRequest::Ready(SelectedConstruct {
            id: id.into(),
            name: id.into(),
            choices: vec![
                ChoiceSelection {
                    bin_id: "backbone".into(),
                    choice_id: "backbone-a".into(),
                },
                ChoiceSelection {
                    bin_id: "insert".into(),
                    choice_id: insert.into(),
                },
            ],
            target,
            overrides: None,
        })
    };
    let matrix = batch::Request {
        design: Some(design),
        rows: vec![
            row(
                "construct-a",
                "insert-a",
                TargetDeclaration::OriginalRegions {
                    substitutions: vec![],
                    edits: vec![],
                    junction_edits: vec![],
                },
            ),
            row(
                "construct-b",
                "insert-b",
                TargetDeclaration::OriginalRegions {
                    substitutions: vec![],
                    edits: vec![],
                    junction_edits: vec![],
                },
            ),
            row(
                "construct-wrong",
                "insert-a",
                TargetDeclaration::Existing(binding(&frame, 3)),
            ),
        ],
    };
    let mut session = batch::Session::begin(frame, matrix).unwrap();
    let context = session.summary().context_revision;
    for index in 0..3 {
        session.advance(&context, index).unwrap();
    }
    let commit = session.into_capture(&context).unwrap();
    (commit.workspace, commit.attachment_hash)
}
fn selection(id: &str, name: &str) -> Selection {
    Selection {
        construct_id: id.into(),
        options: OutputOptions {
            name: name.into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    }
}
fn begin(frame: &Workspace, batch_hash: &str, fail: bool) -> Request {
    let request = release::Request {
        batch_hash: batch_hash.into(),
        rows: vec![
            selection("construct-a", "First published construct"),
            selection(
                if fail {
                    "construct-wrong"
                } else {
                    "construct-b"
                },
                "Second published construct",
            ),
        ],
    };
    Request {
        id: 131,
        revision: release::revision(frame, &request).unwrap(),
        job: AnalysisJob::BeginCombinatorialRelease {
            workspace: frame.clone(),
            request,
        },
    }
}
fn advance(context: &str, index: usize) -> Request {
    Request {
        id: 132 + index as u32,
        revision: context.into(),
        job: AnalysisJob::AdvanceCombinatorialRelease {
            context: context.into(),
            index,
        },
    }
}
fn inspect(context: &str, index: usize) -> Request {
    Request {
        id: 133,
        revision: context.into(),
        job: AnalysisJob::InspectCombinatorialRelease {
            context: context.into(),
            index,
        },
    }
}
fn finish(frame: &Workspace, context: &str) -> Request {
    Request {
        id: 134,
        revision: context.into(),
        job: AnalysisJob::FinishCombinatorialRelease {
            workspace: frame.clone(),
            current_context: context.into(),
        },
    }
}
fn tag(wire: &[u8]) -> u32 {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    tag
}
fn ask(session: &mut AnalysisSession, request: &Request) -> (Vec<u8>, Response) {
    let reply = session
        .process_request(&protocol::request::Request::encode(request).unwrap())
        .unwrap();
    let response = protocol::response::Response::decode(&reply).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    (reply, response)
}
fn response_wire(response: &Response) -> Vec<u8> {
    let body = postcard::to_allocvec(response).unwrap();
    let mut wire = b"CGR\0".to_vec();
    wire.extend_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire.extend_from_slice(&(body.len() as u32).to_le_bytes());
    wire.extend_from_slice(blake3::hash(&body).as_bytes());
    wire.extend(body);
    wire
}
fn completed() -> (Workspace, Request, Vec<Vec<u8>>, release::Commit) {
    let (frame, hash) = fixture();
    let request = begin(&frame, &hash, false);
    let mut session = AnalysisSession::default();
    let (mut replies, _) = {
        let (wire, response) = ask(&mut session, &request);
        (vec![wire], response)
    };
    for index in 0..2 {
        replies.push(ask(&mut session, &advance(&request.revision, index)).0);
    }
    replies.push(ask(&mut session, &inspect(&request.revision, 0)).0);
    let (wire, response) = ask(&mut session, &finish(&frame, &request.revision));
    replies.push(wire);
    let AnalysisResult::CombinatorialReleaseFinished(commit) = response.result else {
        panic!("missing release commit")
    };
    (frame, request, replies, commit)
}

#[test]
fn fresh_rows_are_inspectable_privately_and_only_a_complete_session_can_publish() {
    let (frame, hash) = fixture();
    let request = begin(&frame, &hash, false);
    let fingerprint = run_lineage::workspace_fingerprint(&frame).unwrap();
    let mut session = AnalysisSession::default();
    assert_eq!(
        tag(&protocol::request::Request::encode(&request).unwrap()),
        129
    );
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap()
        )
        .is_err()
    );
    let (wire, response) = ask(&mut session, &request);
    assert_eq!(tag(&wire), 121);
    let AnalysisResult::CombinatorialReleaseStarted(progress) = response.result else {
        panic!("missing release start")
    };
    assert_eq!(
        (progress.next_index, progress.total_rows, progress.finished),
        (0, 2, false)
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&inspect(&request.revision, 0)).unwrap()
            )
            .is_err()
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&advance(&request.revision, 1)).unwrap()
            )
            .is_err()
    );
    for (index, dna, name) in [
        (0, "AACCTTAA", "First published construct"),
        (1, "AACCCCAA", "Second published construct"),
    ] {
        let (wire, response) = ask(&mut session, &advance(&request.revision, index));
        assert_eq!(tag(&wire), 122);
        let AnalysisResult::CombinatorialReleaseAdvanced(progress) = response.result else {
            panic!("missing release row")
        };
        assert!(progress.failed.is_none());
        assert_eq!(progress.next_index, index + 1);
        assert_eq!(progress.finished, index == 1);
        let (wire, response) = ask(&mut session, &inspect(&request.revision, index));
        assert_eq!(tag(&wire), 123);
        let AnalysisResult::CombinatorialReleaseInspected(inspection) = response.result else {
            panic!("missing private models")
        };
        assert_eq!(inspection.index, index);
        assert_eq!(inspection.product.sequence.as_str(), dna);
        assert_eq!(inspection.target.sequence.as_str(), dna);
        assert_eq!(inspection.product.name, name);
        assert_eq!(
            inspection.product.revision(),
            inspection.row.output.construct.revision
        );
    }
    let (wire, response) = ask(&mut session, &finish(&frame, &request.revision));
    assert_eq!(tag(&wire), 124);
    let AnalysisResult::CombinatorialReleaseFinished(commit) = response.result else {
        panic!("missing commit")
    };
    commit
        .validate_binding(&request.revision, &fingerprint)
        .unwrap();
    assert!(commit.workspace.documents.len() > frame.documents.len());
    for original in &frame.documents {
        let preserved = commit
            .workspace
            .documents
            .iter()
            .find(|doc| doc.id == original.id)
            .unwrap();
        assert_eq!(preserved.head, original.head);
        assert_eq!(preserved.head_event, original.head_event);
    }
    let manifest = release::load(&commit.workspace, &commit.manifest_hash).unwrap();
    assert_eq!(manifest.rows.len(), 2);
    assert_eq!(manifest.rows[0].construct_id, "construct-a");
    assert_eq!(manifest.rows[1].construct_id, "construct-b");
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&finish(&frame, &request.revision)).unwrap()
            )
            .is_err()
    );
}

#[test]
fn a_failed_row_or_stale_finish_destroys_the_private_prefix() {
    let (frame, hash) = fixture();
    let request = begin(&frame, &hash, true);
    let mut session = AnalysisSession::default();
    ask(&mut session, &request);
    ask(&mut session, &advance(&request.revision, 0));
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&advance(&request.revision, 1)).unwrap()
            )
            .is_err()
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&inspect(&request.revision, 0)).unwrap()
            )
            .is_err()
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&finish(&frame, &request.revision)).unwrap()
            )
            .is_err()
    );
    assert!(release::list(&frame).unwrap().is_empty());
    let request = begin(&frame, &hash, false);
    ask(&mut session, &request);
    ask(&mut session, &advance(&request.revision, 0));
    ask(&mut session, &advance(&request.revision, 1));
    let mut changed = frame.clone();
    changed.name = "A different live workspace".into();
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&finish(&changed, &request.revision)).unwrap()
            )
            .is_err()
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&finish(&frame, &request.revision)).unwrap()
            )
            .is_err()
    );
}

#[test]
fn request_context_and_native_inspection_receipts_reject_stale_options_and_models() {
    let (frame, hash) = fixture();
    let request = begin(&frame, &hash, false);
    let mut session = AnalysisSession::default();
    let mut wrong = request.clone();
    wrong.revision = "0".repeat(64);
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
    let mut wrong = request.clone();
    if let AnalysisJob::BeginCombinatorialRelease { request, .. } = &mut wrong.job {
        request.rows[0].options.name = "Changed requested name".into();
    }
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
    ask(&mut session, &request);
    ask(&mut session, &advance(&request.revision, 0));
    let (_, response) = ask(&mut session, &inspect(&request.revision, 0));
    let mut wrong = response.clone();
    wrong.revision = "0".repeat(64);
    assert!(protocol::response::Response::decode(&response_wire(&wrong)).is_err());
    let mut wrong = response;
    if let AnalysisResult::CombinatorialReleaseInspected(inspection) = &mut wrong.result {
        inspection.product.sequence = collagene_model::Sequence::dna("AAAA").unwrap();
    }
    assert!(protocol::response::Response::decode(&response_wire(&wrong)).is_err());
    let mut fresh = AnalysisSession::default();
    assert!(
        fresh
            .process_request(
                &protocol::request::Request::encode(&inspect(&request.revision, 0)).unwrap()
            )
            .is_err()
    );
}

#[test]
fn schema57_full_and_selector_only_prefixes_cannot_decode_new_release_payloads() {
    let (frame, request, replies, commit) = completed();
    let exports = export::Options {
        genbank: true,
        fasta: true,
        primer_csv: false,
        primer_tsv: false,
        protocol_html: true,
        protocol_pdf: false,
        construct_ids: None,
    };
    let archive_request = Request {
        id: 135,
        revision: export::revision(&commit.workspace, &commit.manifest_hash, &exports).unwrap(),
        job: AnalysisJob::ExportCombinatorialOutputs {
            workspace: commit.workspace.clone(),
            manifest_hash: commit.manifest_hash.clone(),
            options: exports,
        },
    };
    let mut session = AnalysisSession::default();
    let (archive_reply, _) = ask(&mut session, &archive_request);
    let requests = [
        request.clone(),
        advance(&request.revision, 0),
        inspect(&request.revision, 0),
        finish(&frame, &request.revision),
        archive_request,
    ];
    let responses = [
        replies[0].clone(),
        replies[1].clone(),
        replies[3].clone(),
        replies[4].clone(),
        archive_reply,
    ];
    for (request, response, job, result) in requests
        .into_iter()
        .zip(responses)
        .zip([
            (129u32, 121u32),
            (130, 122),
            (131, 123),
            (132, 124),
            (133, 125),
        ])
        .map(|((request, response), (job, result))| (request, response, job, result))
    {
        for (mut wire, is_response, selector) in [
            (
                protocol::request::Request::encode(&request).unwrap(),
                false,
                job,
            ),
            (response, true, result),
        ] {
            assert_eq!(tag(&wire), selector);
            wire[4..6].copy_from_slice(&57u16.to_le_bytes());
            let rejected = if is_response {
                protocol::response::Response::decode(&wire).map(|_| ())
            } else {
                session.process_request(&wire).map(|_| ())
            };
            assert!(matches!(rejected, Err(Error::Unsupported(_))));
            let body = postcard::to_allocvec(&(1u32, "context", selector)).unwrap();
            wire.truncate(42);
            wire[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
            wire.extend(body);
            let rejected = if is_response {
                protocol::response::Response::decode(&wire).map(|_| ())
            } else {
                session.process_request(&wire).map(|_| ())
            };
            assert!(matches!(rejected, Err(Error::Unsupported(_))));
        }
    }
}

#[test]
fn released_manifest_exports_through_the_worker_and_checks_archive_integrity() {
    let (_, _, _, commit) = completed();
    let options = export::Options {
        genbank: true,
        fasta: true,
        primer_csv: false,
        primer_tsv: false,
        protocol_html: true,
        protocol_pdf: false,
        construct_ids: None,
    };
    let request = Request {
        id: 135,
        revision: export::revision(&commit.workspace, &commit.manifest_hash, &options).unwrap(),
        job: AnalysisJob::ExportCombinatorialOutputs {
            workspace: commit.workspace.clone(),
            manifest_hash: commit.manifest_hash.clone(),
            options,
        },
    };
    let (wire, response) = ask(&mut AnalysisSession::default(), &request);
    assert_eq!(tag(&wire), 125);
    let AnalysisResult::CombinatorialOutputsArchive(archive) = &response.result else {
        panic!("missing portable bundle")
    };
    assert_eq!(archive.extension, "zip");
    assert_eq!(archive.mime, "application/zip");
    assert!(archive.bytes.starts_with(b"PK"));
    assert_eq!(archive.manifest_hash, commit.manifest_hash);
    assert!(archive.files.len() >= 4);
    assert!(
        archive
            .files
            .iter()
            .any(|file| file.construct_id.as_deref() == Some("construct-a"))
    );
    assert!(
        archive
            .files
            .iter()
            .any(|file| file.construct_id.as_deref() == Some("construct-b"))
    );
    let mut wrong = response;
    if let AnalysisResult::CombinatorialOutputsArchive(archive) = &mut wrong.result {
        archive.bytes[0] ^= 1;
    }
    assert!(protocol::response::Response::decode(&response_wire(&wrong)).is_err());
}

#[test]
fn progress_and_inspection_keep_ordinary_limits_and_only_explicit_workspace_variants_expand() {
    let maximum = collagene_core::formats::sequence::MAX_IMPORT_BYTES;
    let mut wire = vec![0; 42 + maximum + 1];
    wire[4..6].copy_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire[6..10].copy_from_slice(&((maximum + 1) as u32).to_le_bytes());
    for (magic, selector, is_response) in [
        (b"CGQ\0", 130u32, false),
        (b"CGQ\0", 131, false),
        (b"CGR\0", 121, true),
        (b"CGR\0", 122, true),
        (b"CGR\0", 123, true),
    ] {
        wire[..4].copy_from_slice(magic);
        let body = postcard::to_allocvec(&(1u32, "context", selector)).unwrap();
        wire[42..42 + body.len()].copy_from_slice(&body);
        let rejected = if is_response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            AnalysisSession::default()
                .process_request(&wire)
                .map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::TooLarge)));
    }
    for (magic, selector, is_response) in [
        (b"CGQ\0", 129u32, false),
        (b"CGQ\0", 132, false),
        (b"CGQ\0", 133, false),
        (b"CGR\0", 124, true),
        (b"CGR\0", 125, true),
    ] {
        wire[..4].copy_from_slice(magic);
        let body = postcard::to_allocvec(&(1u32, "context", selector)).unwrap();
        wire[42..42 + body.len()].copy_from_slice(&body);
        let rejected = if is_response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            AnalysisSession::default()
                .process_request(&wire)
                .map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Format(_))));
    }
}
