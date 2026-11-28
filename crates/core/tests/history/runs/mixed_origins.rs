//! Literal saved PCR → target-free sticky digest → new PCR → saved run. Nested
//! declarations remain instructions; only fresh replay supplies tracked strands.
use collagene_core::{
    cloning::restriction::catalog as catalog_digest,
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace},
    },
    material::{
        capsule::MaterialOriginRef,
        derivation as digest,
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, MaterialOrigin, MaterialSelection},
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const CUT: &str = "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const AMPLICON: &str = "CCGTAGGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn digest_fixture() -> (Workspace, digest::Capsule, EventBinding) {
    let source = Molecule::dna("Original stock", STOCK, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent original PCR target", STOCK, Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Mixed origin literal",
        vec![source.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Original actual PCR".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source)],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("Original F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Original R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved full PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selected = run_material::Selection {
        run_id: saved.receipt.run.run_id.clone(),
        step_index: 0,
        event: saved.receipt.construct.clone(),
    };
    let origin = collagene_core::material::capsule::capsule(&saved.workspace, &selected).unwrap();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let prepared = digest::prepare_install(
        &saved.workspace,
        &digest::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: selected.event,
                origin: digest::SourceOrigin::Saved(
                    collagene_core::material::run::MaterialOrigin::Run(origin),
                ),
            }],
            derivation: digest::Derivation::CatalogDigest {
                conditions: None,
                selection: catalog_digest::select(&catalog, &[enzyme]).unwrap(),
                options: Default::default(),
                fragment_index: 1,
            },
            name: "Native selected digest fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let pin = prepared.receipt.product;
    let cap = digest::capsule(&prepared.workspace, &pin).unwrap();
    (prepared.workspace, cap, pin)
}
fn prepared_plan(
    w: &mut Workspace,
    cap: &digest::Capsule,
    pin: &EventBinding,
) -> (DesignPlan, EventBinding) {
    let target =
        Molecule::dna("Independent second PCR target", AMPLICON, Topology::Linear).unwrap();
    DocumentHistory::new("independent-second-target", target.clone())
        .unwrap()
        .contribute(w)
        .unwrap();
    let target_pin = event(w, w.documents.len() - 1);
    let source = w.revisions[&pin.revision].clone();
    assert_eq!(source.sequence.as_str(), CUT);
    let plan = DesignPlan {
        standard: None,
        name: "PCR from exact saved digest fragment".into(),
        target_revision: target.revision(),
        inputs: vec![input("fragment", source)],
        steps: vec![
            Operation::ReconstructMaterial {
                template: Source::Input("fragment".into()),
                origin: cap.clone(),
            },
            Operation::Pcr {
                template: Source::Step(0),
                forward: Primer::new("New F", "CCGTAGGA", "").unwrap(),
                reverse: Primer::new("New R", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        ],
    };
    (plan, target_pin)
}
#[test]
fn actual_sticky_fragment_then_new_pcr_survives_native_archive_and_mixed_run_replay() {
    let (mut w, cap, pin) = digest_fixture();
    let inspected =
        run_material::inspect_material(&w, &pin, &Some(MaterialSelection::Derived(pin.clone())))
            .unwrap();
    assert!(matches!(
        inspected.selected,
        Some(MaterialOrigin::Derived(_))
    ));
    inspected
        .validate_binding(&inspected.context_revision)
        .unwrap();
    let (plan, target) = prepared_plan(&mut w, &cap, &pin);
    let report = verification::engine::verify(&plan, &w.revisions[&target.revision]);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report
            .product(0)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        CUT
    );
    assert_eq!(
        report
            .product(1)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        AMPLICON
    );
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        target,
        vec![InputEventBinding {
            input_id: "fragment".into(),
            event: pin,
        }],
        &cloning_outputs::Options {
            name: "Saved second PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 1,
        event: saved.receipt.construct,
    };
    let fresh = run_material::reconstruct(&reopened, &selection).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), AMPLICON);
    assert!(fresh.tracked().duplex.nicks.is_empty());
    let native = collagene_core::material::capsule::capsule(&reopened, &selection).unwrap();
    MaterialOrigin::Run(native.clone())
        .validate_for_workspace(&reopened, &selection.event)
        .unwrap();
    let mut plans = Vec::new();
    collagene_core::material::capsule::visit_candidate_plans(&native.candidate, |candidate| {
        plans.push(candidate.plan.name.clone());
        Ok(())
    })
    .unwrap();
    assert_eq!(
        plans,
        vec![
            "PCR from exact saved digest fragment",
            "Original actual PCR"
        ]
    );
    let (reference, tracked, incomplete, findings, total, truncated) =
        digest::replay(&cap, &reopened.revisions[&cap.product.revision])
            .unwrap()
            .into_parts();
    assert_eq!(reference.sequence.as_str(), CUT);
    assert_eq!(tracked.duplex.bottom_start, 4);
    assert!(!incomplete);
    assert!(total >= findings.len());
    assert!(!truncated);
}
#[test]
fn equal_content_alias_and_detached_record_cannot_attest_the_selected_digest_event() {
    let (mut w, cap, pin) = digest_fixture();
    let alias =
        DocumentHistory::new("equal-content-alias", w.revisions[&pin.revision].clone()).unwrap();
    alias.contribute(&mut w).unwrap();
    let alias_pin = event(&w, w.documents.len() - 1);
    assert_eq!(alias_pin.revision, pin.revision);
    assert_ne!(alias_pin.event_id, pin.event_id);
    let (plan, _) = prepared_plan(&mut w, &cap, &pin);
    assert!(
        collagene_core::material::capsule::validate_input_origins(
            &plan,
            &[InputEventBinding {
                input_id: "fragment".into(),
                event: alias_pin.clone()
            }]
        )
        .is_err()
    );
    assert!(
        MaterialOrigin::Derived(cap.clone())
            .validate_for_workspace(&w, &alias_pin)
            .is_err()
    );
    let mut changed = cap.clone();
    changed.source_snapshots[0].molecule.name = "Different source snapshot".into();
    assert!(changed.validate_documentary().is_err());
    let key = w
        .attachments
        .iter()
        .find(|(_, a)| a.kind == digest::WORKSPACE_ATTACHMENT_KIND)
        .unwrap()
        .0
        .clone();
    w.attachments.remove(&key);
    assert!(
        MaterialOrigin::Derived(cap)
            .validate_for_workspace(&w, &pin)
            .is_err()
    );
}
#[test]
fn one_aggregate_ledger_charges_mixed_descendants_instead_of_resetting_each_child() {
    let (mut w, cap, pin) = digest_fixture();
    let record = digest::load_record(&cap.record_bytes).unwrap();
    let digest::SourceOrigin::Saved(MaterialOrigin::Run(parent)) =
        &record.request.sources[0].origin
    else {
        panic!("explicit run parent")
    };
    collagene_core::material::capsule::material_capsules_preflight([
        MaterialOriginRef::Derived(&cap),
        MaterialOriginRef::Run(parent),
    ])
    .unwrap();
    let (plan, target) = prepared_plan(&mut w, &cap, &pin);
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        target,
        vec![InputEventBinding {
            input_id: "fragment".into(),
            event: pin,
        }],
        &cloning_outputs::Options {
            name: "Budget actual saved PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 1,
        event: saved.receipt.construct,
    };
    let native = collagene_core::material::capsule::capsule(&saved.workspace, &selection).unwrap();
    native.validate_documentary().unwrap();
    // Each root is valid (two outer steps plus one native ancestor step).
    // Four hundred uses must not each obtain their own 1000-step allowance.
    assert!(
        collagene_core::material::capsule::material_capsules_preflight(
            (0..400).map(|_| MaterialOriginRef::Run(&native))
        )
        .is_err()
    );
}
