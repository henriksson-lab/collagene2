//! Current worker frames bind actual native strand windows, not RC(displayed top).
use collagene_core::{
    Molecule, Topology,
    assembly_intent::{self, Selector},
    cloning_outputs,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    run_lineage::InputEventBinding,
    run_material::{self, MaterialOrigin, MaterialSelection},
    verification::{DesignPlan, Operation, Source},
    workspace::{self, Workspace},
};
#[path = "common/actual_strand.rs"]
mod native_fixture;
use native_fixture::{event, fixture, input};
#[test]
fn native_bottom_window_current_binary_and_session_reply_keep_independent_literal() {
    let (w, _, source, intent) = fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let revision = assembly_intent::intent_revision(&[source.clone()], &intent, Some(&w)).unwrap();
    let request = Request {
        id: 10401,
        revision: revision.clone(),
        job: AnalysisJob::BuildAssemblyIntent {
            inputs: vec![source],
            request: intent.clone(),
            workspace: Some(w.clone()),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(
        u16::from_le_bytes(bytes[4..6].try_into().unwrap()),
        protocol::WORKER_SCHEMA_VERSION
    );
    let response = protocol::process_request(&bytes).unwrap();
    let AnalysisResult::AssemblyIntent(report) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("bound intended reference");
    };
    assert_eq!(report.target.sequence.as_str(), "AATTCGTAAG");
    assert_ne!(report.target.sequence.as_str(), "CGTAAGAATT");
    assert_eq!(
        report.fragments[0].source_event,
        intent.material_bindings.first().map(|b| b.event.clone())
    );
    report.validate_binding(&revision).unwrap();
    let mut session = protocol::AnalysisSession::default();
    let result = session.process_request(&bytes).unwrap();
    assert!(matches!(
        protocol::decode_response(&result).unwrap().result,
        AnalysisResult::AssemblyIntent(_)
    ));
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let mut missing = serde_json::to_value(&request.job).unwrap();
    missing["BuildAssemblyIntent"]
        .as_object_mut()
        .unwrap()
        .remove("workspace");
    assert!(serde_json::from_value::<AnalysisJob>(missing).is_err());
    let mut no_workspace = request.clone();
    if let AnalysisJob::BuildAssemblyIntent { workspace, .. } = &mut no_workspace.job {
        *workspace = None;
    }
    assert!(protocol::encode_request(&no_workspace).is_err());
    let mut stale = request;
    if let AnalysisJob::BuildAssemblyIntent { request, .. } = &mut stale.job {
        let Selector::MaterialStrandRange { end, .. } = &mut request.fragments[0].selector else {
            unreachable!()
        };
        *end = 9;
    }
    let stale = protocol::encode_request(&stale).unwrap();
    assert!(protocol::process_request(&stale).is_err());
}
#[test]
fn documentary_construction_does_not_replay_or_promote_reference_only_history() {
    let source = Molecule::dna("Reference source", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent reference",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Reference-only native origin",
        vec![source.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Reference overlap".into(),
        target_revision: target.revision(),
        inputs: vec![input("source", source)],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("source".into())],
            overlaps: vec![],
            circular: false,
        }],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Documentary output".into(),
            allow_incomplete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 0,
        event: saved.receipt.construct,
    };
    let origins = run_material::documentary_origins(
        &saved.workspace,
        &[MaterialSelection::Run(selection.clone())],
        &[],
    )
    .unwrap()
    .materials;
    assert_eq!(origins.len(), 1);
    let reference = &saved.workspace.revisions[&selection.event.revision];
    assert!(origins[0].replay(reference).is_err());
    assert!(run_material::capsule(&saved.workspace, &selection).is_err());
    let roots = vec![MaterialSelection::Run(selection); 33];
    assert!(matches!(
        run_material::documentary_origins(&saved.workspace, &roots, &[]),
        Err(collagene_core::Error::TooLarge)
    ));
    assert!(matches!(origins[0], MaterialOrigin::Run(_)));
}
#[test]
fn all_four_workspace_intent_selectors_share_stateless_and_session_large_admission() {
    let mut payload = vec![0; collagene_core::formats::MAX_IMPORT_BYTES + 1];
    for tag in [96u32, 107, 159, 161, 0] {
        let prefix = postcard::to_allocvec(&(10402u32, "bounded-workspace", tag)).unwrap();
        payload[..64].fill(0);
        payload[..prefix.len()].copy_from_slice(&prefix);
        let mut frame = b"CGQ\0".to_vec();
        frame.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
        frame.extend((payload.len() as u32).to_le_bytes());
        frame.extend(blake3::hash(&payload).as_bytes());
        frame.extend_from_slice(&payload);
        let stateless = protocol::process_request(&frame).unwrap_err();
        let session = protocol::AnalysisSession::default()
            .process_request(&frame)
            .unwrap_err();
        // These are deliberately invalid bodies. Only genuine workspace selectors
        // reach native current-layout decoding; ordinary jobs retain their32MiB cap.
        assert_eq!(
            matches!(stateless, collagene_core::Error::TooLarge),
            tag == 0
        );
        assert_eq!(matches!(session, collagene_core::Error::TooLarge), tag == 0);
    }
}
