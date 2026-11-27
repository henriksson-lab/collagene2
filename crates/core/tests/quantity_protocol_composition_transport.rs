#[path = "common/protocol_quantities.rs"]
mod fixture;
#[path = "common/pdf.rs"]
mod pdf_readback;
use collagene_core::{
    Molecule, Topology,
    history::DocumentHistory,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    protocol_quantities as q, run_lineage,
    transfer::Candidate,
    verification::{self, Status},
    workspace,
};
fn export(w: &workspace::Workspace, hash: &str, format: q::ExportFormat) -> q::Exported {
    let context = q::export_revision(w, hash, format).unwrap();
    let request = Request {
        id: 718,
        revision: context.clone(),
        job: AnalysisJob::ExportProtocolQuantities {
            workspace: w.clone(),
            quantity_hash: hash.into(),
            format,
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        protocol::WORKER_SCHEMA_VERSION
    );
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 186);
    let a = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    let b = protocol::decode_response(
        &protocol::AnalysisSession::default()
            .process_request(&bytes)
            .unwrap(),
    )
    .unwrap();
    assert_eq!(a.id, 718);
    assert_eq!(a.revision, context);
    assert_eq!(
        serde_json::to_value(&a).unwrap(),
        serde_json::to_value(&b).unwrap()
    );
    let AnalysisResult::ProtocolQuantitiesExported(exported) = a.result else {
        panic!("quantity protocol export");
    };
    exported.validate_binding(&context).unwrap();
    assert_eq!(
        exported.content_hash,
        blake3::hash(&exported.bytes).to_hex().to_string()
    );
    exported
}
#[test]
fn full_instructions_and_quantities_bind_original_events_after_both_heads_advance() {
    let (w, r) = fixture::fixture();
    let captured = q::capture(&w, &r).unwrap();
    let mut w = workspace::load_workspace(&workspace::save_workspace(&captured.workspace).unwrap())
        .unwrap();
    let run_id = match &r.plan_source {
        q::PlanSource::SavedRun { run_id } => run_id,
        _ => unreachable!(),
    };
    let original = run_lineage::load_run(&w, run_id).unwrap();
    fixture::advance(&mut w, 0, "Later original stock head");
    fixture::advance(&mut w, 1, "Later independent target head");
    let before = workspace::save_workspace(&w).unwrap();
    for format in [q::ExportFormat::ProtocolHtml, q::ExportFormat::ProtocolPdf] {
        let exported = export(&w, &captured.attachment_hash, format);
        assert!(exported.arithmetic_complete);
        assert_eq!(
            exported.source_workspace_fingerprint.as_deref(),
            Some(run_lineage::workspace_fingerprint(&w).unwrap().as_str())
        );
        assert_eq!(exported.target, r.target);
        assert_eq!(exported.inputs, r.inputs);
        assert_eq!(exported.plan_revision, r.plan_revision);
        let molecular = exported.molecular.as_ref().unwrap();
        assert_eq!(molecular.status, Status::ConsistentUnderAssumptions);
        assert_eq!(molecular.engine_version, verification::ENGINE_VERSION);
        assert_eq!(molecular.plan_revision, r.plan_revision);
        assert_eq!(molecular.target_revision, r.target.revision);
        assert_eq!(molecular.report_revision.len(), 64);
        let text = if format == q::ExportFormat::ProtocolPdf {
            pdf_readback::check_xref(&exported.bytes);
            pdf_readback::extracted(&exported.bytes)
        } else {
            String::from_utf8(exported.bytes.clone()).unwrap()
        };
        assert!(text.contains("PCR"));
        assert!(text.contains("template-stock"));
        assert!(text.contains("2000/1 nL"));
        for id in [
            &captured.attachment_hash,
            &r.plan_revision,
            &r.target.event_id,
            &r.inputs[0].event.event_id,
        ] {
            let unwrapped = text.split_whitespace().collect::<String>();
            assert!(unwrapped.contains(id), "missing {id}");
        }
        let mut wrong = exported.clone();
        wrong.target = fixture::pin(&w, 3);
        assert!(wrong.validate_binding(&exported.context_revision).is_err());
        let mut json = serde_json::to_value(&exported).unwrap();
        json.as_object_mut().unwrap().remove("molecular");
        assert!(serde_json::from_value::<q::Exported>(json).is_err());
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    assert_eq!(
        serde_json::to_value(run_lineage::load_run(&w, run_id).unwrap()).unwrap(),
        serde_json::to_value(original).unwrap()
    );
    let companion = export(&w, &captured.attachment_hash, q::ExportFormat::Html);
    assert!(companion.molecular.is_none());
    assert!(companion.source_workspace_fingerprint.is_none());
    assert!(companion.arithmetic_complete);
}
#[test]
fn complete_arithmetic_does_not_rescue_invalid_molecules_and_missing_amounts_do_not_rewrite_science()
 {
    let (mut w, mut r) = fixture::fixture();
    let run = match &r.plan_source {
        q::PlanSource::SavedRun { run_id } => run_lineage::load_run(&w, run_id).unwrap(),
        _ => unreachable!(),
    };
    let wrong = Molecule::dna(
        "Independently supplied wrong target",
        "AGTCCGATGCTACGTT",
        Topology::Linear,
    )
    .unwrap();
    DocumentHistory::new("wrong-target", wrong.clone())
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let mut plan = run.plan;
    plan.target_revision = wrong.revision();
    let candidate = Candidate::new(plan.clone());
    w.candidates.push(candidate);
    r.plan_source = q::PlanSource::Candidate;
    r.plan_revision = verification::plan_revision(&plan);
    r.target = fixture::pin(&w, w.documents.len() - 1);
    let captured = q::capture(&w, &r).unwrap();
    let exported = export(
        &captured.workspace,
        &captured.attachment_hash,
        q::ExportFormat::ProtocolHtml,
    );
    assert!(exported.arithmetic_complete);
    assert_eq!(exported.molecular.as_ref().unwrap().status, Status::Invalid);
    let direct = verification::verify(&plan, &wrong);
    assert_eq!(direct.status, Status::Invalid);
    assert_eq!(
        exported.molecular.as_ref().unwrap().findings_total,
        direct.findings.len()
    );
    let (w, mut r) = fixture::fixture();
    r.reactions[0].final_volume = None;
    let captured = q::capture(&w, &r).unwrap();
    let exported = export(
        &captured.workspace,
        &captured.attachment_hash,
        q::ExportFormat::ProtocolPdf,
    );
    assert!(!exported.arithmetic_complete);
    assert_eq!(
        exported.molecular.as_ref().unwrap().status,
        Status::ConsistentUnderAssumptions
    );
    let text = pdf_readback::extracted(&exported.bytes);
    assert!(text.contains("missing") || text.contains("Missing"));
}
