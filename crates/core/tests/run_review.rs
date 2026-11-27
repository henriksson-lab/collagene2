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

#[test]
fn edited_saved_pcr_reverifies_original_target_and_records_new_run_without_rewinding_heads() {
    use collagene_core::{
        chemistry::{EndTreatment, FivePrime},
        plan_editing, run_material, transfer,
    };
    let (mut w, id) = fixture();
    let original = run_lineage::load_run(&w, &id).unwrap();
    advance(&mut w, 0, "Later edited stock", "TTTTCCCC");
    advance(&mut w, 1, "Later target draft", "AAAAGGGG");
    w.validate().unwrap();
    let before = serde_json::to_value(&w).unwrap();
    let prepared = run_review::prepare(&w, &id).unwrap();
    prepared.validate_for_workspace(&w).unwrap();
    assert_eq!(prepared.inputs, original.inputs);
    assert_eq!(prepared.target_event, original.target);
    assert_ne!(event(&w, 0), prepared.inputs[0].event);
    assert_ne!(event(&w, 1), prepared.target_event);
    let mut edited = prepared.candidate.plan.clone();
    let inputs = edited
        .inputs
        .iter()
        .map(|i| i.id.clone())
        .collect::<Vec<_>>();
    edited.steps = plan_editing::inserted(
        &inputs,
        &edited.steps,
        1,
        Operation::TreatEnds {
            template: Source::Step(0),
            treatment: EndTreatment::Phosphorylate,
        },
    )
    .unwrap();
    edited.name = "Explicit terminal phosphorylation after saved PCR review".into();
    let candidate = transfer::Candidate::new(edited);
    candidate.validate().unwrap();
    assert_ne!(candidate.plan_revision, prepared.candidate.plan_revision);
    assert_eq!(candidate.plan.target_revision, original.target.revision);
    assert_eq!(
        candidate.plan.inputs[0].revision,
        original.inputs[0].event.revision
    );
    let report = verification::verify(&candidate.plan, &prepared.target);
    assert_eq!(
        report.status,
        verification::Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.len(), 2);
    assert_eq!(report.products[1].sequence.as_str(), "AGTCCGATGCTACGTA");
    let html = verification::export_protocol_html(&candidate.plan, &prepared.target);
    assert!(html.contains("Phosphorylate") && html.contains("Original physical stock"));
    assert!(!html.contains("Later edited stock") && !html.contains("Later target draft"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&candidate.plan, &prepared.target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.plan_revision, candidate.plan_revision);
    assert_eq!(pdf.target_revision, prepared.target_event.revision);
    assert_eq!(pdf.status, report.status);
    // Recheck the unchanged native receipt at admission. The edited candidate
    // stays separate from the immutable saved review and its recorded verdict.
    prepared.validate_for_workspace(&w).unwrap();
    assert_eq!(before, serde_json::to_value(&w).unwrap());
    let committed = run_lineage::record_run(
        &w,
        &candidate.plan,
        prepared.target_event.clone(),
        prepared.inputs.clone(),
        false,
    )
    .unwrap();
    assert_ne!(committed.receipt.run_id, id);
    assert_eq!(
        serde_json::to_value(&committed.workspace.documents[..w.documents.len()]).unwrap(),
        before["documents"]
    );
    assert_eq!(
        serde_json::to_value(&committed.workspace.history[..w.history.len()]).unwrap(),
        before["history"]
    );
    assert_eq!(committed.workspace.active_document, w.active_document);
    for (hash, attachment) in &w.attachments {
        assert_eq!(
            serde_json::to_value(&committed.workspace.attachments[hash]).unwrap(),
            serde_json::to_value(attachment).unwrap()
        );
    }
    assert_eq!(
        serde_json::to_value(run_lineage::load_run(&committed.workspace, &id).unwrap()).unwrap(),
        serde_json::to_value(&original).unwrap()
    );
    let restarted = collagene_core::workspace::load_workspace(
        &collagene_core::workspace::save_workspace(&committed.workspace).unwrap(),
    )
    .unwrap();
    let new = run_lineage::load_run(&restarted, &committed.receipt.run_id).unwrap();
    assert_eq!(new.target, original.target);
    assert_eq!(new.inputs, original.inputs);
    let replay = run_lineage::replay_run(&restarted, &new.run_id).unwrap();
    assert!(replay.engine_matches && replay.products_match && replay.status_matches);
    assert_eq!(
        replay.current_report.target_revision,
        original.target.revision
    );
    let final_product = &new.products[1];
    let fresh = run_material::reconstruct(
        &restarted,
        &run_material::Selection {
            run_id: new.run_id,
            step_index: 1,
            event: EventBinding {
                document_id: final_product.document_id.clone(),
                event_id: final_product.event_id.clone(),
                revision: final_product.revision.clone(),
            },
        },
    )
    .unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), "AGTCCGATGCTACGTA");
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        FivePrime::Phosphate
    );
}

#[test]
fn saved_edit_receipt_refuses_equal_snapshot_event_aliases_and_stale_native_context() {
    use collagene_core::history::DocumentHistory;
    let (mut w, id) = fixture();
    let original = run_lineage::load_run(&w, &id).unwrap();
    let target_snapshot = w.revisions[&original.target.revision].clone();
    let alias = DocumentHistory::new("equal-target-alias", target_snapshot).unwrap();
    let alias_id = alias.id.clone();
    alias.contribute(&mut w).unwrap();
    let target_alias = w.documents.iter().position(|d| d.id == alias_id).unwrap();
    advance(&mut w, 0, "Later stock", "TTTTCCCC");
    advance(&mut w, 1, "Later intended draft", "AAAAGGGG");
    let prepared = run_review::prepare(&w, &id).unwrap();
    prepared.validate_for_workspace(&w).unwrap();
    assert_eq!(event(&w, 2).revision, prepared.inputs[0].event.revision);
    assert_eq!(
        event(&w, target_alias).revision,
        prepared.target_event.revision
    );
    let mut wrong = prepared.clone();
    wrong.inputs[0].event = event(&w, 2);
    assert!(wrong.validate_for_workspace(&w).is_err());
    let mut wrong = prepared.clone();
    wrong.target_event = event(&w, target_alias);
    assert!(wrong.validate_for_workspace(&w).is_err());
    let mut advanced = w.clone();
    advance(&mut advanced, 2, "Later private alias", "GCGCGCGC");
    assert!(prepared.validate_for_workspace(&advanced).is_err());
    // A valid re-encoding preserves the literal record/assessment and all native
    // event ownership, while changing the exact retained attachment byte pin.
    let mut rebound = w.clone();
    let mut attachment = rebound
        .attachments
        .remove(&prepared.run_attachment_hash)
        .unwrap();
    attachment.bytes.push(b'\n');
    let replacement_hash = blake3::hash(&attachment.bytes).to_hex().to_string();
    assert_ne!(replacement_hash, prepared.run_attachment_hash);
    rebound
        .attachments
        .insert(replacement_hash.clone(), attachment);
    rebound.validate().unwrap();
    let latest = run_review::prepare(&rebound, &id).unwrap();
    assert_eq!(latest.run_attachment_hash, replacement_hash);
    assert_eq!(
        serde_json::to_value(&latest.recorded_assessment).unwrap(),
        serde_json::to_value(&prepared.recorded_assessment).unwrap()
    );
    assert_eq!(latest.target_event, prepared.target_event);
    assert_eq!(latest.inputs, prepared.inputs);
    assert!(prepared.validate_for_workspace(&rebound).is_err());
}
