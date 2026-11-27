#[path = "common/protocol_quantities.rs"]
mod fixture;
use collagene_core::{
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    protocol_quantities as q, run_lineage, workspace,
};
fn execute(job: AnalysisJob, revision: String, tag: u32) -> AnalysisResult {
    let bytes = protocol::encode_request(&Request {
        id: 291,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        protocol::WORKER_SCHEMA_VERSION
    );
    let ((id, revision_pin, selector), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(id, 291);
    assert_eq!(revision_pin, revision);
    assert_eq!(selector, tag);
    let a = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    let b = protocol::decode_response(
        &protocol::AnalysisSession::default()
            .process_request(&bytes)
            .unwrap(),
    )
    .unwrap();
    assert_eq!(a.id, 291);
    assert_eq!(a.revision, revision);
    assert_eq!(
        serde_json::to_value(&a).unwrap(),
        serde_json::to_value(&b).unwrap()
    );
    a.result
}
#[test]
fn native_saved_run_quantity_preview_capture_open_exports_keep_exact_events_and_status() {
    let (w, r) = fixture::fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let run_id = match &r.plan_source {
        q::PlanSource::SavedRun { run_id } => run_id,
        _ => unreachable!(),
    };
    let original = run_lineage::load_run(&w, run_id).unwrap();
    let ctx = q::preview_revision(&w, &r).unwrap();
    let AnalysisResult::ProtocolQuantitiesPreview(p) = execute(
        AnalysisJob::PreviewProtocolQuantities {
            workspace: w.clone(),
            request: r.clone(),
        },
        ctx.clone(),
        185,
    ) else {
        panic!("preview");
    };
    assert!(p.complete);
    assert_eq!(
        p.reactions[0].fill.as_ref().unwrap().volume,
        fixture::volume(15)
    );
    assert_eq!(p.request.inputs, r.inputs);
    assert_eq!(p.request.target, r.target);
    let AnalysisResult::ProtocolQuantitiesCaptured(c) = execute(
        AnalysisJob::CaptureProtocolQuantities {
            workspace: w.clone(),
            request: r.clone(),
        },
        ctx,
        186,
    ) else {
        panic!("capture");
    };
    c.validate_binding(&c.context_revision).unwrap();
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    assert_eq!(
        serde_json::to_value(&c.workspace.documents).unwrap(),
        serde_json::to_value(&w.documents).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&c.workspace.history).unwrap(),
        serde_json::to_value(&w.history).unwrap()
    );
    let listed_context = q::list_revision(&c.workspace).unwrap();
    let AnalysisResult::ProtocolQuantitiesListed(listed) = execute(
        AnalysisJob::ListProtocolQuantities {
            workspace: c.workspace.clone(),
        },
        listed_context.clone(),
        189,
    ) else {
        panic!("list");
    };
    listed.validate_binding(&listed_context).unwrap();
    assert_eq!(listed.entries.len(), 1);
    assert_eq!(listed.entries[0].target, r.target);
    assert_eq!(listed.entries[0].inputs, r.inputs);
    let mut restored =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    fixture::advance(&mut restored, 0, "New stock head");
    fixture::advance(&mut restored, 1, "New target head");
    let ctx = q::open_revision(&restored, &c.attachment_hash).unwrap();
    let AnalysisResult::ProtocolQuantitiesOpened(o) = execute(
        AnalysisJob::OpenProtocolQuantities {
            workspace: restored.clone(),
            quantity_hash: c.attachment_hash.clone(),
        },
        ctx.clone(),
        187,
    ) else {
        panic!("open");
    };
    o.validate_binding(&ctx).unwrap();
    assert_ne!(o.context_revision, o.preview.context_revision);
    assert_eq!(o.preview.request.inputs, r.inputs);
    assert_eq!(o.preview.request.target, r.target);
    assert!(o.preview.complete);
    assert_eq!(
        serde_json::to_value(&original.assessment).unwrap(),
        serde_json::to_value(&run_lineage::load_run(&restored, run_id).unwrap().assessment)
            .unwrap()
    );
    for format in [
        q::ExportFormat::Csv,
        q::ExportFormat::Html,
        q::ExportFormat::Pdf,
    ] {
        let ctx = q::export_revision(&restored, &c.attachment_hash, format).unwrap();
        let AnalysisResult::ProtocolQuantitiesExported(e) = execute(
            AnalysisJob::ExportProtocolQuantities {
                workspace: restored.clone(),
                quantity_hash: c.attachment_hash.clone(),
                format,
            },
            ctx.clone(),
            188,
        ) else {
            panic!("export");
        };
        e.validate_binding(&ctx).unwrap();
        assert_eq!(e.quantity_hash, c.attachment_hash);
        assert_eq!(e.content_hash, blake3::hash(&e.bytes).to_hex().to_string());
        if format == q::ExportFormat::Pdf {
            assert!(e.bytes.starts_with(b"%PDF-"));
        } else {
            assert!(
                String::from_utf8(e.bytes)
                    .unwrap()
                    .contains("PRIVATE EXPLICIT QUANTITY PLAN")
            );
        }
    }
}
#[test]
fn missing_quantities_are_drafts_and_equal_dna_events_or_stale_context_do_not_substitute() {
    let (w, r) = fixture::fixture();
    for target in [false, true] {
        let mut wrong = r.clone();
        if target {
            wrong.target = fixture::pin(&w, 3);
        } else {
            wrong.inputs[0].event = fixture::pin(&w, 2);
        }
        assert!(q::preview_revision(&w, &wrong).is_err());
    }
    let mut draft = r.clone();
    draft.reactions[0].final_volume = None;
    let ctx = q::preview_revision(&w, &draft).unwrap();
    let AnalysisResult::ProtocolQuantitiesPreview(p) = execute(
        AnalysisJob::PreviewProtocolQuantities {
            workspace: w.clone(),
            request: draft.clone(),
        },
        ctx.clone(),
        185,
    ) else {
        panic!("draft");
    };
    assert!(!p.complete);
    assert!(!p.findings.is_empty());
    assert!(p.reactions[0].fill.is_none());
    draft.reactions[0].components[0].volume = q::VolumeSpec::Direct {
        volume: Some(fixture::volume(4)),
    };
    let raw = protocol::encode_request(&Request {
        id: 1,
        revision: ctx,
        job: AnalysisJob::CaptureProtocolQuantities {
            workspace: w.clone(),
            request: draft,
        },
    })
    .unwrap();
    assert!(protocol::process_request(&raw).is_err());
    let mut over = r.clone();
    over.reactions[0].components[0].volume = q::VolumeSpec::Direct {
        volume: Some(fixture::volume(30)),
    };
    let expected = q::preview_revision(&w, &over).unwrap();
    let AnalysisResult::ProtocolQuantitiesPreview(p) = execute(
        AnalysisJob::PreviewProtocolQuantities {
            workspace: w,
            request: over,
        },
        expected,
        185,
    ) else {
        panic!("overfilled");
    };
    assert!(!p.complete);
    assert!(!p.findings.is_empty());
}
#[test]
fn binary_mass_and_molar_dilutions_preserve_explicit_rounding_and_achieved_dimension() {
    use q::units::{
        self, MassConcentration, MassUnit, MolarConcentration, MolarUnit, Rational, VolumeUnit,
    };
    let (w, mut r) = fixture::fixture();
    r.reactions[0].components[0].volume = q::VolumeSpec::MassDilution {
        recipe: units::DilutionRecipe {
            schema_version: units::SCHEMA_VERSION,
            stock: Some(
                MassConcentration::from_units(
                    Rational::integer(100),
                    MassUnit::Nanogram,
                    VolumeUnit::Microliter,
                )
                .unwrap(),
            ),
            target: Some(
                MassConcentration::from_units(
                    Rational::integer(10),
                    MassUnit::Nanogram,
                    VolumeUnit::Microliter,
                )
                .unwrap(),
            ),
            final_volume: Some(fixture::volume(20)),
        },
        quantization: Some(q::Quantization {
            increment: fixture::volume(3),
            mode: units::RoundingMode::Ceiling,
        }),
    };
    r.reactions[0].components[1].volume = q::VolumeSpec::MolarDilution {
        recipe: units::DilutionRecipe {
            schema_version: units::SCHEMA_VERSION,
            stock: Some(
                MolarConcentration::from_unit(Rational::integer(10), MolarUnit::Micromolar)
                    .unwrap(),
            ),
            target: Some(
                MolarConcentration::from_unit(Rational::integer(1), MolarUnit::Micromolar).unwrap(),
            ),
            final_volume: Some(fixture::volume(20)),
        },
        quantization: None,
    };
    let expected = q::preview_revision(&w, &r).unwrap();
    let AnalysisResult::ProtocolQuantitiesPreview(p) = execute(
        AnalysisJob::PreviewProtocolQuantities {
            workspace: w,
            request: r,
        },
        expected,
        185,
    ) else {
        panic!("dilution");
    };
    let q::DilutionResult::QuantizedMass(mass) =
        p.reactions[0].components[0].calculation.as_ref().unwrap()
    else {
        panic!("massdimension");
    };
    assert_eq!(mass.exact.stock_volume, fixture::volume(2));
    assert_eq!(mass.stock.rounded, fixture::volume(3));
    assert!(!mass.target_reached);
    assert_eq!(
        mass.achieved
            .in_units(MassUnit::Nanogram, VolumeUnit::Microliter)
            .unwrap(),
        Rational::integer(15)
    );
    let q::DilutionResult::Molar(molar) =
        p.reactions[0].components[1].calculation.as_ref().unwrap()
    else {
        panic!("molardimension");
    };
    assert_eq!(molar.stock_volume, fixture::volume(2));
    assert_eq!(
        p.reactions[0].fill.as_ref().unwrap().volume,
        fixture::volume(15)
    );
    assert!(
        p.complete,
        "Arithmetic definedness does not assert exact target concentration after explicit rounding"
    );
}
#[test]
fn candidate_deletion_keeps_discoverable_unavailable_receipt_without_reconstructing_a_fallback() {
    let (mut w, mut r) = fixture::fixture();
    let run = match &r.plan_source {
        q::PlanSource::SavedRun { run_id } => run_lineage::load_run(&w, run_id).unwrap(),
        _ => unreachable!(),
    };
    let candidate = collagene_core::transfer::Candidate::new(run.plan);
    w.active_candidate = Some(candidate.plan_revision.clone());
    w.candidates.push(candidate);
    r.plan_source = q::PlanSource::Candidate;
    let captured = q::capture(&w, &r).unwrap();
    let mut w = captured.workspace;
    w.candidates.clear();
    w.active_candidate = None;
    w.validate().unwrap();
    let ctx = q::list_revision(&w).unwrap();
    let AnalysisResult::ProtocolQuantitiesListed(listed) = execute(
        AnalysisJob::ListProtocolQuantities {
            workspace: w.clone(),
        },
        ctx.clone(),
        189,
    ) else {
        panic!("catalog");
    };
    listed.validate_binding(&ctx).unwrap();
    assert_eq!(listed.entries.len(), 1);
    assert!(!listed.entries[0].available);
    assert!(!listed.entries[0].findings.is_empty());
    assert_eq!(listed.entries[0].attachment_hash, captured.attachment_hash);
    assert!(q::load(&w, &captured.attachment_hash).is_err());
    assert!(q::open_revision(&w, &captured.attachment_hash).is_err());
    assert!(q::export_revision(&w, &captured.attachment_hash, q::ExportFormat::Html).is_err());
}
