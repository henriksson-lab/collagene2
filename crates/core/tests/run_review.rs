use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_review,
    verification::{self, DesignPlan, Operation, PlanInput, Source},
    workspace::{RevisionEvent, Workspace, WorkspaceAttachment},
};
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, String) {
    let original = Molecule::dna(
        "Original physical stock",
        "AGTCCGATGCTACGTA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Separately declared intended PCR",
        "AGTCCGATGCTACGTA",
        Topology::Linear,
    )
    .unwrap();
    let alias = original.clone();
    let w = Workspace::from_documents(
        "Original bindings and unrelated same-DNA alias",
        vec![original.clone(), target.clone(), alias],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Retained exact PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: original.revision(),
            molecule: original,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCC", "").unwrap(),
            reverse: Primer::new("R", "TACGT", "").unwrap(),
        }],
    };
    let saved = run_lineage::record_run(
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
    (saved.workspace, saved.receipt.run_id)
}
fn advance(w: &mut Workspace, index: usize, name: &str, dna: &str) {
    let molecule = Molecule::dna(name, dna, Topology::Linear).unwrap();
    let d = &mut w.documents[index];
    let e = RevisionEvent::new(
        d.id.clone(),
        molecule.revision(),
        vec![d.head_event.clone().unwrap()],
        "Independent subsequent edit",
    )
    .unwrap();
    d.head = molecule.revision();
    d.head_event = Some(e.id.clone());
    w.revisions.insert(molecule.revision(), molecule);
    w.history.push(e);
}
#[test]
fn review_uses_original_events_and_independent_target_without_moving_edited_heads() {
    let (mut w, id) = fixture();
    let original = run_lineage::load_run(&w, &id).unwrap();
    advance(&mut w, 0, "Later edited stock", "TTTTCCCC");
    advance(&mut w, 1, "Later target draft", "AAAAGGGG");
    w.validate().unwrap();
    let before = serde_json::to_value(&w).unwrap();
    let prepared = run_review::prepare(&w, &id).unwrap();
    assert_eq!(prepared.target.sequence.as_str(), "AGTCCGATGCTACGTA");
    assert_eq!(prepared.target_event, original.target);
    assert_eq!(prepared.inputs, original.inputs);
    assert_eq!(
        serde_json::to_value(&prepared.candidate.plan).unwrap(),
        serde_json::to_value(&original.plan).unwrap()
    );
    assert_eq!(prepared.candidate.plan_revision, original.plan_revision);
    assert!(prepared.recorded_engine_matches);
    prepared.validate_for_workspace(&w).unwrap();
    assert_eq!(before, serde_json::to_value(&w).unwrap());
    advance(&mut w, 2, "Private branch edit", "ACACACAC");
    assert!(prepared.validate_for_workspace(&w).is_err());
}
#[test]
fn equal_snapshot_alias_cannot_replace_original_input_or_target_event() {
    let (w, id) = fixture();
    let prepared = run_review::prepare(&w, &id).unwrap();
    assert_eq!(
        w.head(&w.documents[0].id).unwrap().revision(),
        w.head(&w.documents[2].id).unwrap().revision()
    );
    let mut changed = prepared.clone();
    changed.inputs[0].event = event(&w, 2);
    assert!(changed.validate_for_workspace(&w).is_err());
    let mut changed = prepared.clone();
    changed.target_event = event(&w, 2);
    assert!(changed.validate_for_workspace(&w).is_err());
    let mut missing = w.clone();
    missing
        .history
        .retain(|e| e.id != prepared.inputs[0].event.event_id);
    assert!(run_review::prepare(&missing, &id).is_err());
    let mut missing = w.clone();
    missing
        .history
        .retain(|e| e.id != prepared.target_event.event_id);
    assert!(run_review::prepare(&missing, &id).is_err());
    let mut altered = prepared.clone();
    altered.candidate.plan.steps.clear();
    assert!(altered.validate_for_workspace(&w).is_err());
}
#[test]
fn historical_engine_is_documentary_and_never_adopted_as_current_status() {
    let (mut w, id) = fixture();
    let mut r = run_lineage::load_run(&w, &id).unwrap();
    r.assessment.engine_version = "historical-test-engine/science-1".into();
    // Construct a structurally valid historical archive using the documented native identity.
    let identity = serde_json::to_vec(&(
        "cloning-run-v1",
        r.schema_version,
        r.plan_schema,
        &r.plan_revision,
        &r.plan,
        &r.target,
        &r.inputs,
        r.products
            .iter()
            .map(|p| (&p.step_index, &p.revision))
            .collect::<Vec<_>>(),
        &r.assessment,
    ))
    .unwrap();
    r.run_id = blake3::hash(&identity).to_hex().to_string();
    let p = &mut r.products[0];
    let old_event = p.event_id.clone();
    let e = RevisionEvent::new(
        p.document_id.clone(),
        p.revision.clone(),
        vec![r.inputs[0].event.event_id.clone()],
        format!("Cloning run {} step 1", r.run_id),
    )
    .unwrap();
    p.event_id = e.id.clone();
    w.history.retain(|e| e.id != old_event);
    w.history.push(e.clone());
    let d = w
        .documents
        .iter_mut()
        .find(|d| d.id == p.document_id)
        .unwrap();
    d.head_event = Some(e.id);
    w.attachments
        .retain(|_, a| a.kind != run_lineage::WORKSPACE_ATTACHMENT_KIND);
    let bytes = run_lineage::save_run_record(&r).unwrap();
    w.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        WorkspaceAttachment {
            kind: run_lineage::WORKSPACE_ATTACHMENT_KIND.into(),
            name: "Historical native run".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    w.validate().unwrap();
    let prepared = run_review::prepare(&w, &r.run_id).unwrap();
    assert!(!prepared.recorded_engine_matches);
    assert_eq!(
        prepared.recorded_assessment.engine_version,
        "historical-test-engine/science-1"
    );
    assert_eq!(
        prepared.current_verification_engine,
        verification::ENGINE_VERSION
    );
    assert!(prepared.warnings.iter().any(|w| w.contains("differs")));
    prepared.candidate.validate().unwrap();
    prepared.validate_for_workspace(&w).unwrap();
}
