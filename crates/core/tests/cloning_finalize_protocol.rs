use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    cloning_outputs::{self, Options},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    run_lineage::{self, EventBinding, InputEventBinding},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
};
fn request() -> Request {
    let source = Molecule::dna("original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna("independent target", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "literal exact PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("forward", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("reverse", "GTTACCCC", "").unwrap(),
        }],
    };
    let workspace = Workspace::from_documents("finalization", vec![source, target], 1).unwrap();
    let event = |index: usize| {
        let d = &workspace.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    };
    let target = event(1);
    let inputs = vec![InputEventBinding {
        input_id: "source".into(),
        event: event(0),
    }];
    let options = Options {
        name: "Saved construct".into(),
        save_primers: true,
        ..Default::default()
    };
    let revision = cloning_outputs::finalize_revision(
        &workspace,
        &plan,
        target.clone(),
        inputs.clone(),
        &options,
    )
    .unwrap();
    Request {
        id: 190,
        revision,
        job: AnalysisJob::FinalizeCloningOutputs {
            workspace,
            plan,
            target,
            inputs,
            options,
        },
    }
}
#[test]
fn output_commit_roundtrips_with_replayable_history_and_rejects_old_headers() {
    let request = request();
    let bytes = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 100);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&38u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 92);
    let decoded = protocol::decode_response(&bytes).unwrap();
    assert_eq!(decoded.id, request.id);
    let AnalysisResult::CloningOutputsFinalized(commit) = decoded.result else {
        panic!("expected atomic output commit")
    };
    commit.validate_binding(&request.revision).unwrap();
    assert_eq!(
        commit
            .workspace
            .head(&commit.receipt.construct.document_id)
            .unwrap()
            .name,
        "Saved construct"
    );
    assert!(commit.receipt.primer_attachment.is_some());
    assert!(
        run_lineage::replay_run(&commit.workspace, &commit.receipt.run.run_id)
            .unwrap()
            .products_match
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&38u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
}
#[test]
fn changed_output_options_cannot_reuse_an_old_workspace_context() {
    let mut request = request();
    if let AnalysisJob::FinalizeCloningOutputs { options, .. } = &mut request.job {
        options.show_prepared = true;
    }
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
}
