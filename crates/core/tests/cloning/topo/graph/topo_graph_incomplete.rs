//! Chemical correspondence cannot upgrade incomplete source history.
#[path = "../../../common/topo_graph.rs"]
mod literal;
use collagene_core::cloning::topo::physical as p;
use collagene_core::{
    cloning::topo::{self as t, graph_history as g, run_preparation as native},
    documents::{
        history::DocumentHistory,
        workspace::{Workspace, WorkspaceAttachment},
    },
    material::{
        lineage::InputEventBinding,
        outputs as cloning_outputs,
        run::{self as run_material, MaterialOrigin},
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
fn check(
    w: &Workspace,
    selection: &g::Selection,
) -> (t::Request, native::Task, DesignPlan, Molecule) {
    let request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::SavedGraph {
                material_id: "retained conditional graph".into(),
                origin: g::documentary_capsule(w, selection).unwrap(),
            }],
        },
        steps: vec![],
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    let target = w.revisions[&literal::event(w, 4).revision].clone();
    let task = native::Task::Check {
        request: request.clone(),
        selection: t::Selection {
            batch_revision: t::batch_revision(&[], &request).unwrap(),
            component: 0,
            reference_seed: p::Call {
                strand: 0,
                index: 0,
            },
        },
        target: target.clone(),
        policy: t::AcceptancePolicy {
            allow_nicks: true,
            allow_unapplied_steps: true,
            ..Default::default()
        },
    };
    let plan = native::compile(&[], &target, None, &task).unwrap();
    (request, task, plan, target)
}
fn forged_false_record(w: &Workspace, selection: &g::Selection) -> (Workspace, g::Selection) {
    let mut record = g::load(w, &selection.record_hash).unwrap();
    assert!(record.incomplete);
    record.incomplete = false;
    let bytes = g::save_record(&record).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut w = w.clone();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: g::WORKSPACE_ATTACHMENT_KIND.into(),
            name: "False completeness claim (diagnostic fixture)".into(),
            media_type: g::MEDIA_TYPE.into(),
            bytes,
        },
    );
    w.validate().unwrap();
    (
        w,
        g::Selection {
            record_hash: hash,
            state_revision: selection.state_revision.clone(),
        },
    )
}
#[test]
fn completed_ta_graph_with_blocked_third_transfer_stays_incomplete_despite_target_acceptance() {
    let first = literal::captured();
    let mut r = literal::partial(&first.workspace, &first.receipt.selection);
    r.request.steps.extend([
        p::Step::Transfer(p::Transfer {
            donor: p::Call {
                strand: 1,
                index: 12,
            },
            acceptor: p::Call {
                strand: 3,
                index: 0,
            },
            geometry: p::Geometry::ThreePrimeTA,
        }),
        p::Step::Transfer(p::Transfer {
            donor: p::Call {
                strand: 0,
                index: 12,
            },
            acceptor: p::Call {
                strand: 2,
                index: 0,
            },
            geometry: p::Geometry::ThreePrimeTA,
        }),
    ]);
    let saved = g::prepare_capture(&first.workspace, &r).unwrap();
    assert!(saved.receipt.summary.incomplete);
    let actual =
        g::replay(&g::documentary_capsule(&saved.workspace, &saved.receipt.selection).unwrap())
            .unwrap();
    assert!(actual.incomplete());
    assert_eq!(actual.report().physical.product.source_pairing.len(), 4);
    assert_eq!(
        actual.report().physical.product.chains[0].sequence.as_str(),
        literal::FINAL
    );
    assert_eq!(
        actual.report().physical.product.released_activations.len(),
        2
    );
    assert_eq!(
        actual.report().physical.steps[2].outcome(),
        p::Outcome::BlockedUnderDeclaredModel
    );
    let (_, task, plan, target) = check(&saved.workspace, &saved.receipt.selection);
    let inspected = native::inspect(&[], None, &task).unwrap();
    assert!(inspected.origin_incomplete);
    let native::Reply::Checked(accepted) = &inspected.reply else {
        panic!("checked task")
    };
    assert!(matches!(
        accepted.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    let verified = verification::engine::verify(&plan, &target);
    assert_eq!(verified.reconstructed_count(), 1);
    assert_eq!(
        verified
            .product(0)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        literal::FINAL
    );
    assert_eq!(verified.status, Status::Incomplete);
    assert!(
        verified
            .findings
            .iter()
            .any(|s| s.message.contains("incomplete saved native graph"))
    );
    let (forged, selection) = forged_false_record(&saved.workspace, &saved.receipt.selection);
    let (_, task, bad_plan, target) = check(&forged, &selection);
    assert!(native::inspect(&[], None, &task).is_err());
    assert_eq!(
        verification::engine::verify(&bad_plan, &target).status,
        Status::Invalid,
        "a framed false flag is not chemical source authority"
    );
}
#[test]
fn actual_capped_ordinary_origin_uncertainty_survives_new_pcr_and_graph_source_reuse() {
    use collagene_core::{
        chemistry::{AnnealOptions, Terminus},
        cloning::ligation::tailing::TailingBase,
        sequence::analysis::Primer,
    };
    let (mut w, vector) = literal::fixture();
    for (id, m) in [
        (
            "left-template",
            Molecule::dna("left template", "CGAATG", Topology::Linear).unwrap(),
        ),
        (
            "right-template",
            Molecule::dna("right ordered oligo", "CATTCG", Topology::Linear).unwrap(),
        ),
        (
            "parent-target",
            Molecule::dna(
                "Independent aligned annealed template",
                "CGAATG",
                Topology::Linear,
            )
            .unwrap(),
        ),
    ] {
        DocumentHistory::new(id, m)
            .unwrap()
            .contribute(&mut w)
            .unwrap();
    }
    let left = w.revisions[&literal::event(&w, 6).revision].clone();
    let right = w.revisions[&literal::event(&w, 7).revision].clone();
    let target = w.revisions[&literal::event(&w, 8).revision].clone();
    // CGAATG has no proper matching prefix/suffix, so the first retained
    // complementary placement is the actual full0-offset duplex. Candidate1
    // stops the positive-offset search; no uniqueness/complete-search claim.
    let options = AnnealOptions {
        min_paired_bases: 1,
        max_candidates: 1,
        top_chemistry: Terminus::PCR,
        bottom_chemistry: Terminus::PCR,
        ..Default::default()
    };
    let preview =
        verification::terminal::simulate_oligo_annealing(&left, &right, &options).unwrap();
    assert!(!preview.annealing.complete);
    assert_eq!(preview.annealing.candidates.len(), 1);
    let retained = &preview.annealing.candidates[0];
    assert_eq!(retained.paired_bases, 6);
    assert_eq!(retained.duplex.bottom_start, 0);
    assert_eq!(retained.duplex.top.as_str(), "CGAATG");
    assert_eq!(retained.duplex.bottom.as_str(), "CATTCG");
    assert!(retained.duplex.nicks.is_empty());
    let plan = DesignPlan {
        name: "Real capped annealing placement search".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left.revision(),
                molecule: left,
            },
            PlanInput {
                id: "right".into(),
                revision: right.revision(),
                molecule: right,
            },
        ],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("left".into()),
            bottom: Source::Input("right".into()),
            options,
            candidate_index: 0,
        }],
    };
    let report = verification::engine::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::Incomplete,
        "the real annealing placement search cap must actually be reached"
    );
    assert_eq!(report.reconstructed_count(), 1);
    assert_eq!(
        report
            .product(0)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        "CGAATG"
    );
    let captured = cloning_outputs::finalize(
        &w,
        &plan,
        literal::event(&w, 8),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: literal::event(&w, 6),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: literal::event(&w, 7),
            },
        ],
        &cloning_outputs::Options {
            name: "Conditional intact annealed PCR template".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: captured.receipt.run.run_id,
        step_index: 0,
        event: captured.receipt.construct,
    };
    let origin = MaterialOrigin::Run(
        collagene_core::material::capsule::capsule(&captured.workspace, &selection).unwrap(),
    );
    let source = &captured.workspace.revisions[&selection.event.revision];
    assert!(origin.replay(source).unwrap().incomplete());
    let mut request = vector;
    request.name = "Source-conditional new insert PCR".into();
    request.bindings.push(g::SourceBinding {
        role: "saved-template".into(),
        event: selection.event,
        origin: Some(origin),
    });
    request
        .request
        .recipe
        .materials
        .push(t::MaterialRecipe::ExactPcr {
            material_id: "new insert".into(),
            template_role: "saved-template".into(),
            forward: Primer::new("insert F", "CGA", "").unwrap(),
            reverse: Primer::new("insert R", "CAT", "").unwrap(),
            tail: Some(TailingBase::A),
        });
    request.request.steps = vec![
        p::Step::Transfer(p::Transfer {
            donor: p::Call {
                strand: 0,
                index: 12,
            },
            acceptor: p::Call {
                strand: 2,
                index: 0,
            },
            geometry: p::Geometry::ThreePrimeTA,
        }),
        p::Step::Transfer(p::Transfer {
            donor: p::Call {
                strand: 1,
                index: 12,
            },
            acceptor: p::Call {
                strand: 3,
                index: 0,
            },
            geometry: p::Geometry::ThreePrimeTA,
        }),
    ];
    let graph = g::prepare_capture(&captured.workspace, &request).unwrap();
    assert!(graph.receipt.summary.incomplete);
    let (_, task, plan, target) = check(&graph.workspace, &graph.receipt.selection);
    let inspected = native::inspect(&[], None, &task).unwrap();
    assert!(inspected.origin_incomplete);
    let native::Reply::Checked(accepted) = &inspected.reply else {
        panic!("checked task")
    };
    assert!(
        accepted
            .product
            .projection
            .physical
            .all_requested_steps_applied
    );
    assert!(matches!(
        accepted.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    assert_eq!(
        verification::engine::verify(&plan, &target).status,
        Status::Incomplete
    );
    let (forged, selection) = forged_false_record(&graph.workspace, &graph.receipt.selection);
    let (_, task, bad_plan, target) = check(&forged, &selection);
    assert!(native::inspect(&[], None, &task).is_err());
    assert_eq!(
        verification::engine::verify(&bad_plan, &target).status,
        Status::Invalid
    );
}
