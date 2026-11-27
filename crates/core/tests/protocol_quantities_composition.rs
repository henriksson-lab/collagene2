use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    history::DocumentHistory,
    protocol_quantities::{self as q, units as u, *},
    run_lineage::{self, EventBinding, InputEventBinding},
    transfer::Candidate,
    verification::{self, DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
fn v(n: u64, d: u64) -> u::Volume {
    u::Volume::from_unit(u::Rational::new(n, d).unwrap(), u::VolumeUnit::Microliter).unwrap()
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn base_fixture() -> (Workspace, DesignPlan, Request) {
    let stock = Molecule::dna("Literal stock16", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
    let target =
        Molecule::dna("Independent target16", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
    let mut w = Workspace::from_documents(
        "Declared quantities",
        vec![stock.clone(), target.clone(), stock.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Original-template PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCC", "").unwrap(),
            reverse: Primer::new("R", "TACGT", "").unwrap(),
        }],
    };
    w.candidates.push(Candidate::new(plan.clone()));
    w.validate().unwrap();
    let request = Request {
        name: "Explicit20uL planning".into(),
        plan_source: PlanSource::Candidate,
        plan_revision: verification::plan_revision(&plan),
        target: event(&w, 1),
        inputs: vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        reactions: vec![Reaction {
            id: "PCR".into(),
            step_index: 0,
            final_volume: Some(v(20, 1)),
            fill: FillPolicy::ToFinal {
                component_id: "fill".into(),
                label: "Explicit water".into(),
            },
            components: vec![
                Component {
                    id: "dna".into(),
                    label: "Declared DNA aliquot".into(),
                    source: Some(Source::Input("stock".into())),
                    volume: VolumeSpec::Direct {
                        volume: Some(v(2, 1)),
                    },
                },
                Component {
                    id: "buffer".into(),
                    label: "Explicit buffer quantity".into(),
                    source: None,
                    volume: VolumeSpec::Direct {
                        volume: Some(v(3, 1)),
                    },
                },
            ],
        }],
    };
    (w, plan, request)
}

#[path = "common/pdf.rs"]
mod pdf_readback;
use pdf_readback::{check_xref, extracted};
fn fixture() -> (Workspace, DesignPlan, Request) {
    let (mut w, mut plan, mut r) = base_fixture();
    plan.steps.push(Operation::TreatEnds {
        template: Source::Step(0),
        treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
    });
    plan.name = "Retained PCR followed by explicit phosphorylation".into();
    w.candidates = vec![Candidate::new(plan.clone())];
    r.plan_revision = verification::plan_revision(&plan);
    r.reactions[0].components[0].label = "DNA <aliquot> & exact".into();
    let mut second = r.reactions[0].clone();
    second.id = "Terminal treatment".into();
    second.step_index = 1;
    second.components[0].source = Some(Source::Step(0));
    r.reactions.push(second);
    let run =
        run_lineage::record_run(&w, &plan, r.target.clone(), r.inputs.clone(), false).unwrap();
    r.plan_source = PlanSource::SavedRun {
        run_id: run.receipt.run_id,
    };
    (run.workspace, plan, r)
}
fn normalized_pdf(bytes: &[u8]) -> String {
    extracted(bytes)
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect()
}
#[test]
fn saved_pcr_treatment_restart_exports_one_checked_protocol_with_exact_quantity_tables() {
    let (w, plan, r) = fixture();
    let original = workspace::save_workspace(&w).unwrap();
    let c = q::capture(&w, &r).unwrap();
    let restored =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    let target = restored.revisions.get(&r.target.revision).unwrap();
    let fresh = verification::verify(&plan, target);
    assert_eq!(
        fresh.status,
        verification::Status::ConsistentUnderAssumptions
    );
    for format in [ExportFormat::ProtocolHtml, ExportFormat::ProtocolPdf] {
        let e = q::export(&restored, &c.attachment_hash, format).unwrap();
        e.validate_binding(&q::export_revision(&restored, &c.attachment_hash, format).unwrap())
            .unwrap();
        e.validate_for_preview(&c.preview).unwrap();
        assert!(e.arithmetic_complete);
        let m = e.molecular.as_ref().unwrap();
        assert_eq!(m.status, fresh.status);
        assert_eq!(m.plan_revision, verification::plan_revision(&plan));
        assert_eq!(e.target, r.target);
        assert_eq!(
            e.source_workspace_fingerprint.as_deref(),
            Some(
                run_lineage::workspace_fingerprint(&restored)
                    .unwrap()
                    .as_str()
            )
        );
        if format == ExportFormat::ProtocolHtml {
            let text = String::from_utf8(e.bytes).unwrap();
            assert_eq!(text.matches("<!doctype html>").count(), 1);
            assert_eq!(text.matches("</html>").count(), 1);
            assert!(text.contains("<table>") && text.contains("Phosphorylate"));
            assert!(text.contains("AGTCC") && text.contains("TACGT"));
            assert!(
                text.contains("DNA &lt;aliquot&gt; &amp; exact") && !text.contains("DNA <aliquot>")
            );
            assert!(text.contains("15000/1 nL") && text.contains(&c.attachment_hash));
            assert!(
                text.contains(&r.target.event_id) && text.contains(&r.inputs[0].event.event_id)
            );
            assert!(text.contains("Unlisted operations and omitted reagents remain unspecified"));
        } else {
            check_xref(&e.bytes);
            let text = normalized_pdf(&e.bytes);
            for expected in [
                "Declaredreactionquantities",
                "Phosphorylate",
                "AGTCC",
                "TACGT",
                "15000/1nL",
                &c.attachment_hash,
                &r.target.event_id,
                &r.inputs[0].event.event_id,
            ] {
                assert!(text.contains(expected), "missing {expected}");
            }
        }
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), original);
    // Standalone molecular exports remain actual useful renderer paths.
    assert!(
        !verification::export_protocol_html(&plan, target)
            .contains("<h2>Declared reaction quantities</h2>")
    );
    assert!(
        collagene_core::pdf::export_protocol_pdf(&plan, target)
            .unwrap()
            .bytes
            .starts_with(b"%PDF-")
    );
}
#[test]
fn arithmetic_and_molecular_statuses_are_independent_and_neither_is_restored() {
    let (w, _, mut r) = fixture();
    r.reactions[0].components[0].volume = VolumeSpec::Direct { volume: None };
    let c = q::capture(&w, &r).unwrap();
    let e = q::export(&c.workspace, &c.attachment_hash, ExportFormat::ProtocolHtml).unwrap();
    assert!(!e.arithmetic_complete);
    assert_eq!(
        e.molecular.as_ref().unwrap().status,
        verification::Status::ConsistentUnderAssumptions
    );
    let (mut invalid, mut plan, mut r) = base_fixture();
    // Independent target is a frozen literal different from the real PCR output.
    let wrong = Molecule::dna(
        "Different independent target16",
        "TTTTGGGGAAAACCCC",
        Topology::Linear,
    )
    .unwrap();
    let mut histories = DocumentHistory::restore_all(&invalid).unwrap();
    histories[1]
        .commit(wrong, "Explicit independent target change")
        .unwrap();
    invalid.documents.clear();
    invalid.history.clear();
    for h in histories {
        h.contribute(&mut invalid).unwrap();
    }
    r.target = event(&invalid, 1);
    plan.target_revision = r.target.revision.clone();
    invalid.candidates = vec![Candidate::new(plan.clone())];
    r.plan_revision = verification::plan_revision(&plan);
    let c = q::capture(&invalid, &r).unwrap();
    let e = q::export(&c.workspace, &c.attachment_hash, ExportFormat::ProtocolPdf).unwrap();
    assert!(e.arithmetic_complete);
    assert_eq!(e.molecular.unwrap().status, verification::Status::Invalid);
    assert!(normalized_pdf(&e.bytes).contains("INVALID"));
}
#[test]
fn historic_exact_pins_survive_heads_while_changed_inputs_plan_or_hash_refuse() {
    let (w, _, r) = fixture();
    let c = q::capture(&w, &r).unwrap();
    let original = q::export(&c.workspace, &c.attachment_hash, ExportFormat::ProtocolHtml).unwrap();
    let companion =
        q::export_revision(&c.workspace, &c.attachment_hash, ExportFormat::Html).unwrap();
    let mut later = c.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&later).unwrap();
    for index in [0, 1] {
        histories[index]
            .commit(
                Molecule::dna("Unrelated later head", "TTTTGGGG", Topology::Linear).unwrap(),
                "User edit after recording",
            )
            .unwrap();
    }
    later.documents.clear();
    later.history.clear();
    for h in histories {
        h.contribute(&mut later).unwrap();
    }
    let e = q::export(&later, &c.attachment_hash, ExportFormat::ProtocolHtml).unwrap();
    assert_eq!(e.target, r.target);
    assert_eq!(e.inputs[0].event, r.inputs[0].event);
    assert_ne!(e.context_revision, original.context_revision);
    assert_eq!(
        q::export_revision(&later, &c.attachment_hash, ExportFormat::Html).unwrap(),
        companion
    );
    assert!(original.validate_binding(&e.context_revision).is_err());
    let mut alias = e.clone();
    alias.target = event(&later, 2);
    assert!(alias.validate_binding(&e.context_revision).is_err());
    let mut other = c.preview.clone();
    other.request.inputs[0].event = event(&later, 2);
    assert!(e.validate_for_preview(&other).is_err());
    let mut wrong = e.clone();
    wrong.quantity_hash = "0".repeat(64);
    assert!(wrong.validate_binding(&e.context_revision).is_err());
    let mut changed = verification::plan_revision(&later.candidates[0].plan);
    assert_eq!(changed, r.plan_revision);
    later.candidates[0].plan.steps.reverse();
    changed = verification::plan_revision(&later.candidates[0].plan);
    assert_ne!(changed, r.plan_revision);
    assert!(q::export(&later, &c.attachment_hash, ExportFormat::ProtocolHtml).is_err());
    // A changed retained Candidate graph cannot silently replace its original record.
    let (w, _, r) = base_fixture();
    let c = q::capture(&w, &r).unwrap();
    let mut missing = c.workspace;
    missing.candidates.clear();
    missing.active_candidate = None;
    assert!(q::export(&missing, &c.attachment_hash, ExportFormat::ProtocolHtml).is_err());
}
#[test]
fn combined_exact_and_quantized_dilutions_keep_achieved_values_and_quantity_identity() {
    let (w, _, mut r) = fixture();
    r.reactions.truncate(1);
    r.reactions[0].components.truncate(1);
    r.reactions[0].components[0].volume = VolumeSpec::MassDilution {
        recipe: u::DilutionRecipe {
            schema_version: u::SCHEMA_VERSION,
            stock: Some(
                u::MassConcentration::from_units(
                    u::Rational::integer(10),
                    u::MassUnit::Nanogram,
                    u::VolumeUnit::Microliter,
                )
                .unwrap(),
            ),
            target: Some(
                u::MassConcentration::from_units(
                    u::Rational::new(1, 10).unwrap(),
                    u::MassUnit::Nanogram,
                    u::VolumeUnit::Microliter,
                )
                .unwrap(),
            ),
            final_volume: Some(v(20, 1)),
        },
        quantization: None,
    };
    let c = q::capture(&w, &r).unwrap();
    assert_eq!(c.preview.reactions[0].components[0].volume, Some(v(1, 5)));
    let text = String::from_utf8(
        q::export(&c.workspace, &c.attachment_hash, ExportFormat::ProtocolHtml)
            .unwrap()
            .bytes,
    )
    .unwrap();
    assert!(text.contains("200/1 nL") && text.contains("19800/1 nL"));
    r.reactions[0].final_volume = Some(v(10, 1));
    r.reactions[0].components[0].volume = VolumeSpec::MolarDilution {
        recipe: u::DilutionRecipe {
            schema_version: u::SCHEMA_VERSION,
            stock: Some(
                u::MolarConcentration::from_unit(
                    u::Rational::integer(1000),
                    u::MolarUnit::Nanomolar,
                )
                .unwrap(),
            ),
            target: Some(
                u::MolarConcentration::from_unit(u::Rational::integer(15), u::MolarUnit::Nanomolar)
                    .unwrap(),
            ),
            final_volume: Some(v(10, 1)),
        },
        quantization: Some(Quantization {
            increment: v(1, 10),
            mode: u::RoundingMode::NearestTiesUp,
        }),
    };
    let c = q::capture(&w, &r).unwrap();
    let Some(DilutionResult::QuantizedMolar(calculation)) =
        &c.preview.reactions[0].components[0].calculation
    else {
        panic!("quantized molar calculation")
    };
    assert_eq!(calculation.stock.difference, v(1, 20));
    assert!(!calculation.target_reached);
    assert_eq!(
        calculation
            .achieved
            .in_unit(u::MolarUnit::Nanomolar)
            .unwrap(),
        u::Rational::integer(20)
    );
    let e = q::export(&c.workspace, &c.attachment_hash, ExportFormat::ProtocolPdf).unwrap();
    let text = normalized_pdf(&e.bytes);
    assert!(
        text.contains("achieved") && text.contains("difference") && text.contains("target_reached")
    );
    let mut forged = e.clone();
    forged.arithmetic_complete = false;
    assert!(forged.validate_for_preview(&c.preview).is_err());
    let mut forged = e;
    forged.molecular = None;
    assert!(forged.validate_binding(&forged.context_revision).is_err());
}
