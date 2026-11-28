//! Current-format review/export boundaries preserve exact scientific intent.
use collagene_core::{
    cloning::restriction::inspection as restriction_inspection,
    documents::workspace::Workspace,
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        review as run_review,
    },
    output::map_export,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        self,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::{self, Primer},
};
use collagene_model::{Molecule, Topology};
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
        standard: None,
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
            template_strand: collagene_core::chemistry::Strand::Top,
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
fn call(
    job: AnalysisJob,
    revision: String,
    job_tag: u32,
    reply_tag: u32,
) -> protocol::response::Response {
    assert_eq!(tag(&job), job_tag);
    let encoded = protocol::request::Request::encode(&protocol::request::Request {
        id: 93,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&(protocol::schema::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &encoded[6..]);
    assert!(protocol::execution::process_request(&old).is_err());
    let encoded = protocol::execution::process_request(&encoded).unwrap();
    let response = protocol::response::Response::decode(&encoded).unwrap();
    assert_eq!(response.id, 93);
    assert_eq!(response.revision, revision);
    assert_eq!(tag(&response.result), reply_tag, "{:?}", response.result);
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&(protocol::schema::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(&old[6..], &encoded[6..]);
    assert!(protocol::response::Response::decode(&old).is_err());
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
        173,
        163,
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
    let wrong = protocol::request::Request {
        id: 179,
        revision: "0".repeat(64),
        job: AnalysisJob::PrepareSavedRunReview {
            workspace: w.clone(),
            run_id,
        },
    };
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
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
    let request = map_export::Request {
        source_revision: source.revision(),
        source: source.clone(),
        restriction_sites: Some(sites),
        format: collagene_core::output::raster::Format::Svg,
    };
    let direct = request.render().unwrap();
    let response = call(
        AnalysisJob::RenderMap {
            request: request.clone(),
        },
        "a".repeat(64),
        174,
        164,
    );
    let AnalysisResult::MapExported(map) = response.result else {
        unreachable!()
    };
    map.validate_binding(&source.revision()).unwrap();
    assert_eq!(map.image.bytes(), direct.image.bytes());
    let collagene_core::output::raster::Image::Svg(svg) = &map.image else {
        panic!("expected SVG")
    };
    assert!(svg.contains("A&lt;&amp;&gt;"));
    assert!(svg.contains("data-restriction-total=\"1\""));
    assert!(svg.contains("data-recognition-start=\"2\""));
    assert!(svg.contains("data-top-boundary=\"3\""));
    assert!(svg.contains("data-bottom-boundary=\"7\""));
    assert!(!svg.contains("<script"));
    assert!(!map.warnings.is_empty());
    let mut stale = request;
    stale.source_revision = "b".repeat(64);
    assert!(stale.render().is_err());
}
#[test]
fn current_plan_data_replays_only_through_current_transport() {
    let (w, id) = saved();
    let request = protocol::request::Request {
        id: 1,
        revision: run_lineage::workspace_fingerprint(&w).unwrap(),
        job: AnalysisJob::ReplayCloningRun {
            workspace: w,
            run_id: id,
        },
    };
    let bytes = protocol::request::Request::encode(&request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let reply = protocol::execution::process_request(&bytes).unwrap();
    assert_eq!(
        &reply[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(matches!(
        protocol::response::Response::decode(&reply).unwrap().result,
        AnalysisResult::RunReplay(_)
    ));
    let mut unsupported = bytes.clone();
    unsupported[4..6].copy_from_slice(&(protocol::schema::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert!(protocol::execution::process_request(&unsupported).is_err());
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
            AnalysisJob::RenderMap {
                request: map_export::Request {
                    source_revision: source.revision(),
                    source: source.clone(),
                    restriction_sites: None,
                    format: collagene_core::output::raster::Format::Svg,
                },
            },
            source.revision(),
            false,
        ),
    ];
    let mut session = protocol::session::AnalysisSession::default();
    for (job, revision, is_review) in jobs {
        let request = protocol::request::Request {
            id: 800,
            revision,
            job,
        };
        let mut frame = protocol::request::Request::encode(&request).unwrap();
        // Valid current payload plus trailing bytes: review's larger snapshot
        // allowance must reach complete-payload validation. Map export must
        // refuse the ordinary-message size before decoding its body.
        frame.resize(
            42 + collagene_core::formats::sequence::MAX_IMPORT_BYTES + 1,
            0,
        );
        let length = (frame.len() - 42) as u32;
        frame[6..10].copy_from_slice(&length.to_le_bytes());
        let checksum = blake3::hash(&frame[42..]);
        frame[10..42].copy_from_slice(checksum.as_bytes());
        for outcome in [
            protocol::execution::process_request(&frame),
            session.process_request(&frame),
        ] {
            if is_review {
                assert!(
                    matches!(outcome, Err(collagene_model::Error::Format(_))),
                    "{outcome:?}"
                );
            } else {
                assert!(
                    matches!(outcome, Err(collagene_model::Error::TooLarge)),
                    "{outcome:?}"
                );
            }
        }
    }
}
