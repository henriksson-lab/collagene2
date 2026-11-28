use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::{
    documents::{
        history::DocumentHistory,
        transfer::Candidate,
        workspace::{self, Workspace},
    },
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding},
    output::quantities::{self as q, units as u, *},
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
};
use collagene_model::{Molecule, Topology};
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
fn fixture() -> (Workspace, DesignPlan, Request) {
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
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    w.candidates.push(Candidate::new(plan.clone()));
    w.validate().unwrap();
    let request = Request {
        name: "Explicit20uL planning".into(),
        plan_source: PlanSource::Candidate,
        plan_revision: collagene_core::planning::plan::DesignPlan::revision(&plan),
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
#[test]
fn exact_balanced_quantities_capture_restart_and_three_exports_preserve_biology() {
    let (w, _, r) = fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let expected = q::preview_revision(&w, &r).unwrap();
    let p = q::preview(&w, &r).unwrap();
    assert!(p.complete);
    assert_eq!(
        p.reactions[0].mixture.as_ref().unwrap().specified_volume,
        v(5, 1)
    );
    assert_eq!(p.reactions[0].fill.as_ref().unwrap().volume, v(15, 1));
    let c = q::capture(&w, &r).unwrap();
    c.validate_binding(&expected).unwrap();
    let opened =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    let loaded = q::load(&opened, &c.attachment_hash).unwrap();
    assert_eq!(
        serde_json::to_value(loaded.request).unwrap(),
        serde_json::to_value(r).unwrap()
    );
    for format in [ExportFormat::Csv, ExportFormat::Html, ExportFormat::Pdf] {
        let e = q::export(&opened, &c.attachment_hash, format).unwrap();
        e.validate_binding(&q::export_revision(&opened, &c.attachment_hash, format).unwrap())
            .unwrap();
        assert!(!e.bytes.is_empty());
        if format == ExportFormat::Pdf {
            assert!(e.bytes.starts_with(b"%PDF-"));
        } else {
            let text = String::from_utf8(e.bytes).unwrap();
            assert!(text.contains("15000/1 nL") || text.contains("15000"));
            assert!(text.contains(&c.preview.request.target.event_id));
        }
    }
    let mut original = c.workspace.clone();
    original.attachments.remove(&c.attachment_hash);
    assert_eq!(workspace::save_workspace(&original).unwrap(), before);
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn missing_is_not_zero_unallocated_and_overfilled_recipes_remain_documentary_drafts() {
    let (w, _, mut r) = fixture();
    r.reactions[0].components[0].volume = VolumeSpec::Direct { volume: None };
    r.reactions[0].final_volume = None;
    let missing = q::preview(&w, &r).unwrap();
    assert!(!missing.complete);
    assert!(missing.reactions[0].fill.is_none());
    assert!(missing.reactions[0].mixture.is_none());
    assert!(!q::capture(&w, &r).unwrap().preview.complete);
    r.reactions[0].final_volume = Some(v(20, 1));
    r.reactions[0].components.clear();
    let empty = q::preview(&w, &r).unwrap();
    assert!(empty.complete);
    assert_eq!(empty.reactions[0].fill.as_ref().unwrap().volume, v(20, 1));
    assert!(
        empty
            .warnings
            .iter()
            .any(|s| s.contains("omitted reagents"))
    );
    r.reactions[0].fill = FillPolicy::Unallocated;
    assert!(!q::preview(&w, &r).unwrap().complete);
    r.reactions[0].components.push(Component {
        id: "too-much".into(),
        label: "Overfill".into(),
        source: None,
        volume: VolumeSpec::Direct {
            volume: Some(v(21, 1)),
        },
    });
    let over = q::preview(&w, &r).unwrap();
    assert!(!over.complete);
    assert!(over.reactions[0].mixture.is_none());
    assert!(
        over.findings
            .iter()
            .any(|f| f.message == u::Error::OverfilledMixture.to_string())
    );
}
#[test]
fn declared_mass_and_quantized_molar_dilutions_keep_stock_volume_achieved_and_deviation() {
    let (w, _, mut r) = fixture();
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
    let mass = q::preview(&w, &r).unwrap();
    assert!(mass.complete);
    assert_eq!(mass.reactions[0].components[0].volume, Some(v(1, 5)));
    assert_eq!(mass.reactions[0].fill.as_ref().unwrap().volume, v(99, 5));
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
    let rounded = q::preview(&w, &r).unwrap();
    assert!(rounded.complete);
    let Some(DilutionResult::QuantizedMolar(p)) = &rounded.reactions[0].components[0].calculation
    else {
        panic!("typed molar result required")
    };
    assert_eq!(p.stock.rounded, v(1, 5));
    assert_eq!(p.stock.difference, v(1, 20));
    assert!(!p.target_reached);
    assert_eq!(
        p.achieved.in_unit(u::MolarUnit::Nanomolar).unwrap(),
        u::Rational::integer(20)
    );
    r.reactions[0].final_volume = Some(v(11, 1));
    let mismatch = q::preview(&w, &r).unwrap();
    assert!(!mismatch.complete);
    assert!(mismatch.reactions[0].components[0].volume.is_none());
}
#[test]
fn candidate_fresh_alias_is_explicit_but_cannot_substitute_an_existing_receipt_or_native_run() {
    let (w, plan, mut r) = fixture();
    let p = q::preview(&w, &r).unwrap();
    let mut replaced = p.clone();
    replaced.request.inputs[0].event = event(&w, 2);
    assert!(replaced.validate_binding(&p.context_revision).is_err());
    r.inputs[0].event = event(&w, 2);
    let fresh = q::preview(&w, &r).unwrap();
    assert_ne!(fresh.context_revision, p.context_revision);
    assert!(
        fresh
            .warnings
            .iter()
            .any(|s| s.contains("caller declarations"))
    );
    let captured = run_lineage::record_run(
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
    r.plan_source = PlanSource::SavedRun {
        run_id: captured.receipt.run_id,
    };
    assert!(q::preview(&captured.workspace, &r).is_err());
}
#[test]
fn historic_source_and_target_heads_may_advance_but_new_preview_context_changes() {
    let (w, _, r) = fixture();
    let captured = q::capture(&w, &r).unwrap();
    let mut later = captured.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&later).unwrap();
    for index in [0, 1] {
        let edit_end = histories[index].molecule.sequence.len();
        histories[index]
            .commit(
                vec![
                    Edit::ReplaceSequence {
                        start: 0,
                        end: edit_end,
                        replacement: "TTTTGGGG".into(),
                    },
                    Edit::Rename {
                        name: "Later edited head".into(),
                    },
                    Edit::SetTopology {
                        topology: Topology::Linear,
                    },
                ],
                "Later user edit",
            )
            .unwrap();
    }
    later.documents.clear();
    later.history.clear();
    for h in histories {
        h.contribute(&mut later).unwrap();
    }
    later.validate().unwrap();
    let loaded = q::load(&later, &captured.attachment_hash).unwrap();
    assert_eq!(loaded.request.target, r.target);
    assert_eq!(loaded.request.inputs[0].event, r.inputs[0].event);
    assert_ne!(
        q::preview_revision(&later, &r).unwrap(),
        loaded.context_revision
    );
    let opened = q::open(&later, &captured.attachment_hash).unwrap();
    opened
        .validate_binding(&q::open_revision(&later, &captured.attachment_hash).unwrap())
        .unwrap();
    assert_ne!(opened.context_revision, opened.preview.context_revision);
}
#[test]
fn deleted_candidate_does_not_hide_an_independent_saved_run_quantity_record() {
    let (w, plan, r) = fixture();
    let run = run_lineage::record_run(&w, &plan, event(&w, 1), r.inputs.clone(), false).unwrap();
    let candidate = q::capture(&run.workspace, &r).unwrap();
    let mut saved_request = r.clone();
    saved_request.name = "Saved run quantities".into();
    saved_request.plan_source = PlanSource::SavedRun {
        run_id: run.receipt.run_id,
    };
    let saved = q::capture(&candidate.workspace, &saved_request).unwrap();
    let mut later = saved.workspace;
    later.candidates.clear();
    later.active_candidate = None;
    later.validate().unwrap();
    let listed = q::list(&later).unwrap();
    listed
        .validate_binding(&q::list_revision(&later).unwrap())
        .unwrap();
    assert_eq!(listed.entries.len(), 2);
    assert!(
        listed
            .entries
            .iter()
            .any(|e| e.attachment_hash == candidate.attachment_hash
                && !e.available
                && !e.findings.is_empty())
    );
    assert!(
        listed
            .entries
            .iter()
            .any(|e| e.attachment_hash == saved.attachment_hash && e.available)
    );
    assert!(q::open(&later, &candidate.attachment_hash).is_err());
    q::open(&later, &saved.attachment_hash).unwrap();
}
#[test]
fn current_canonical_record_and_actual_operation_material_references_are_required() {
    let (w, _, mut r) = fixture();
    r.reactions[0].components[0].source = Some(Source::Step(0));
    assert!(q::preview(&w, &r).is_err());
    r.reactions[0].components[0].source = Some(Source::Input("not-this-operation".into()));
    assert!(q::preview(&w, &r).is_err());
    r.reactions[0].components[0].source = None;
    let c = q::capture(&w, &r).unwrap();
    let asset = &c.workspace.attachments[&c.attachment_hash];
    q::load_record(&asset.bytes).unwrap();
    let mut whitespace = asset.bytes.clone();
    whitespace.push(b'\n');
    assert!(q::load_record(&whitespace).is_err());
    let mut forged = c.preview;
    forged.reactions[0].fill.as_mut().unwrap().volume = v(1, 1);
    assert!(forged.validate_binding(&forged.context_revision).is_err());
}
