//! Current-format review/export boundaries preserve exact scientific intent.
use collagene_core::{
    Molecule, Topology,
    analysis::{self, Primer},
    protocol::{self, AnalysisJob, AnalysisResult},
    restriction_inspection,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_review, svg_export,
    verification::{self, DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
};
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved() -> (Workspace, String) {
    let source = Molecule::dna("Stock", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent expected PCR",
        "AGTCCGATGCTACGTA",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Saved exact intent",
        vec![source.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Retained PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCC", "").unwrap(),
            reverse: Primer::new("R", "TACGT", "").unwrap(),
        }],
    };
    let result = run_lineage::record_run(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        false,
    )
    .unwrap();
    (result.workspace, result.receipt.run_id)
}
fn tag(value: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(value).unwrap())
        .unwrap()
        .0
}
fn call(job: AnalysisJob, revision: String, job_tag: u32, reply_tag: u32) -> protocol::Response {
    assert_eq!(tag(&job), job_tag);
    let encoded = protocol::encode_request(&protocol::Request {
        id: 93,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &encoded[6..]);
    assert!(protocol::process_request(&old).is_err());
    let encoded = protocol::process_request(&encoded).unwrap();
    let response = protocol::decode_response(&encoded).unwrap();
    assert_eq!(response.id, 93);
    assert_eq!(response.revision, revision);
    assert_eq!(tag(&response.result), reply_tag, "{:?}", response.result);
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &encoded[6..]);
    assert!(protocol::decode_response(&old).is_err());
    response
}
#[test]
fn review178_reply168_preserves_exact_independent_intent_without_adopting_status() {
    let (w, run_id) = saved();
    let before = serde_json::to_value(&w).unwrap();
    let revision = run_review::request_revision(&w, &run_id).unwrap();
    let response = call(
        AnalysisJob::PrepareSavedRunReview {
            workspace: w.clone(),
            run_id: run_id.clone(),
        },
        revision.clone(),
        178,
        168,
    );
    let AnalysisResult::SavedRunReviewPrepared(review) = response.result else {
        unreachable!()
    };
    review.validate_request_binding(&revision).unwrap();
    review.validate_for_workspace(&w).unwrap();
    assert_eq!(review.target.sequence.as_str(), "AGTCCGATGCTACGTA");
    assert_eq!(review.target_event, event(&w, 1));
    assert_eq!(review.inputs[0].event, event(&w, 0));
    let mut forged = review.clone();
    forged.candidate.plan.name.push_str(" altered");
    assert!(forged.validate_request_binding(&revision).is_err());
    let wrong = protocol::Request {
        id: 179,
        revision: "0".repeat(64),
        job: AnalysisJob::PrepareSavedRunReview {
            workspace: w.clone(),
            run_id,
        },
    };
    assert!(protocol::process_request(&protocol::encode_request(&wrong).unwrap()).is_err());
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
}
#[test]
fn map179_reply169_keeps_safe_xml_and_exact_supplied_cut_coordinates() {
    let source = Molecule::dna("A<&>", "CCGAATTCTT", Topology::Circular).unwrap();
    let enzyme = analysis::starter_enzymes()
        .into_iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    let sites = restriction_inspection::inspect("map", &source, &[enzyme], &Default::default())
        .unwrap()
        .sites;
    assert_eq!(sites.len(), 1);
    let request = svg_export::Request {
        source_revision: source.revision(),
        source: source.clone(),
        restriction_sites: Some(sites),
    };
    let direct = svg_export::render(&request).unwrap();
    let response = call(
        AnalysisJob::RenderMapSvg {
            request: request.clone(),
        },
        "a".repeat(64),
        179,
        169,
    );
    let AnalysisResult::MapSvgExported(map) = response.result else {
        unreachable!()
    };
    map.validate_binding(&source.revision()).unwrap();
    assert_eq!(map.svg, direct.svg);
    assert!(map.svg.contains("A&lt;&amp;&gt;"));
    assert!(map.svg.contains("data-restriction-total=\"1\""));
    assert!(map.svg.contains("data-recognition-start=\"2\""));
    assert!(map.svg.contains("data-top-boundary=\"3\""));
    assert!(map.svg.contains("data-bottom-boundary=\"7\""));
    assert!(!map.svg.contains("<script"));
    assert!(!map.warnings.is_empty());
    let mut stale = request;
    stale.source_revision = "b".repeat(64);
    assert!(svg_export::render(&stale).is_err());
}
#[test]
fn current_plan_data_replays_only_through_current_transport() {
    let (w, id) = saved();
    let request = protocol::Request {
        id: 1,
        revision: run_lineage::workspace_fingerprint(&w).unwrap(),
        job: AnalysisJob::ReplayCloningRun {
            workspace: w,
            run_id: id,
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(&bytes[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    let reply = protocol::process_request(&bytes).unwrap();
    assert_eq!(&reply[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    assert!(matches!(
        protocol::decode_response(&reply).unwrap().result,
        AnalysisResult::RunReplay(_)
    ));
    let mut unsupported = bytes.clone();
    unsupported[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert!(protocol::process_request(&unsupported).is_err());
}

#[test]
fn session_and_stateless_admission_keep_review_and_map_size_policies_distinct() {
    let (workspace, run_id) = saved();
    let source = workspace.head(&workspace.active_document).unwrap().clone();
    let review_revision = run_review::request_revision(&workspace, &run_id).unwrap();
    let jobs = [
        (
            AnalysisJob::PrepareSavedRunReview { workspace, run_id },
            review_revision,
            true,
        ),
        (
            AnalysisJob::RenderMapSvg {
                request: svg_export::Request {
                    source_revision: source.revision(),
                    source: source.clone(),
                    restriction_sites: None,
                },
            },
            source.revision(),
            false,
        ),
    ];
    let mut session = protocol::AnalysisSession::default();
    for (job, revision, is_review) in jobs {
        let request = protocol::Request {
            id: 800,
            revision,
            job,
        };
        let mut frame = protocol::encode_request(&request).unwrap();
        // Valid current payload plus trailing bytes: review's larger snapshot
        // allowance must reach complete-payload validation. Map export must
        // refuse the ordinary-message size before decoding its body.
        frame.resize(42 + collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
        let length = (frame.len() - 42) as u32;
        frame[6..10].copy_from_slice(&length.to_le_bytes());
        let checksum = blake3::hash(&frame[42..]);
        frame[10..42].copy_from_slice(checksum.as_bytes());
        for outcome in [
            protocol::process_request(&frame),
            session.process_request(&frame),
        ] {
            if is_review {
                assert!(
                    matches!(outcome, Err(collagene_core::Error::Format(_))),
                    "{outcome:?}"
                );
            } else {
                assert!(
                    matches!(outcome, Err(collagene_core::Error::TooLarge)),
                    "{outcome:?}"
                );
            }
        }
    }
}
