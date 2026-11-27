//! Independent literal continuation through actual cleavage/sealing, including
//! all excised calls. No supplied graph or report is a reusable source.
use collagene_core::{
    Molecule, Topology,
    run_lineage::EventBinding,
    topo as t, topo_actions as material, topo_graph_history as history,
    topo_run_preparation as native,
    verification::{self, Status},
    workspace::{self, Workspace},
};
use t::physical as p;
const TARGET: &str = "AAGGGGCATCCCTTCACCATGACG";
fn call(strand: usize, index: usize) -> p::Call {
    p::Call { strand, index }
}
fn transfer(donor: p::Call, acceptor: p::Call, geometry: p::Geometry) -> p::Step {
    p::Step::Transfer(p::Transfer {
        donor,
        acceptor,
        geometry,
    })
}
fn cut() -> p::Step {
    p::Step::Process(p::processing::Action::CleaveThreePrimeFlap{
    bond:p::processing::Bond{upstream:call(3,5),downstream:call(3,6)},displacement_step:0,
    profile:p::processing::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,
    assumption:"Declared ideal hydrolysis at this actual anchored flap bond; no nuclease or host repair claim".into(),
})
}
fn seal(three: p::Call, five: p::Call) -> p::Step {
    p::Step::Process(p::processing::Action::SealNick{
    selection:p::nick_sealing::NickSelection{three,five},
    ligase_assumption:"Explicit active ligase/cofactor hypothesis for this exact continuous-template 3OH/5P nick".into(),
})
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let doc = &w.documents[i];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn strand(role: &str, length: usize, activated: bool) -> t::PreparedStrand {
    t::PreparedStrand {
        id: role.into(),
        segments: vec![t::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: if activated {
            p::FivePrime::Phosphate
        } else {
            p::FivePrime::Hydroxyl
        },
        three_prime: if activated {
            p::ThreePrime::Phosphotyrosyl {
                enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                activation_id: format!("declared-{role}"),
            }
        } else {
            p::ThreePrime::Hydroxyl
        },
    }
}
fn fixture() -> (Workspace, history::CaptureRequest, Vec<t::Input>) {
    let literals = [
        ("vt", "AAGGGGCATCCCTT"),
        ("vb", "GGTGAAGGGATGCCCCTT"),
        ("it", "CACCATGACG"),
        ("ib", "CGTCATGGTG"),
    ];
    let molecules = literals
        .iter()
        .map(|(name, dna)| Molecule::dna(*name, *dna, Topology::Linear).unwrap())
        .collect::<Vec<_>>();
    let w =
        Workspace::from_documents("Whole processed source graph", molecules.clone(), 0).unwrap();
    let inputs = literals
        .iter()
        .zip(molecules)
        .map(|((role, _), molecule)| t::Input {
            role: (*role).into(),
            molecule,
        })
        .collect::<Vec<_>>();
    let request=t::Request{recipe:t::GenerationRecipe{schema_version:1,materials:vec![t::MaterialRecipe::Prepared{
        declaration_id:"Literal two-ended activated vector and OH insert".into(),
        assumption:"Both vector phosphotyrosyl termini explicitly declared; DNA documents alone do not establish activation".into(),
        strands:vec![strand("vt",14,true),strand("vb",18,true),strand("it",10,false),strand("ib",10,false)],
        pairing:vec![p::PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:4,length:14},p::PairSpan{strand_a:2,start_a:0,strand_b:3,start_b:0,length:10}],
    }]},steps:vec![transfer(call(0,13),call(2,0),p::Geometry::DirectionalCacc),cut(),seal(call(3,5),call(1,0))],engine_limits:p::Limits::default(),projection_options:p::projection::Options::default()};
    let capture = history::CaptureRequest {
        name: "One transfer, actual flap hydrolysis and one nick sealed".into(),
        bindings: literals
            .iter()
            .enumerate()
            .map(|(i, (role, _))| history::SourceBinding {
                role: (*role).into(),
                event: event(&w, i),
                origin: None,
            })
            .collect(),
        request,
    };
    (w, capture, inputs)
}
fn resumed(w: &Workspace, selection: &history::Selection) -> t::Request {
    t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::SavedGraph {
                material_id: "Entire52-call processed reaction".into(),
                origin: history::documentary_capsule(w, selection).unwrap(),
            }],
        },
        steps: vec![
            transfer(call(1, 17), call(3, 0), p::Geometry::Flush),
            seal(call(2, 9), call(0, 0)),
        ],
        engine_limits: p::Limits::default(),
        projection_options: p::projection::Options::default(),
    }
}
fn assert_whole_closed_graph(product: &p::ProductGraph) {
    assert_eq!(
        product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
    assert_eq!(product.chains.len(), 3);
    assert_eq!(product.components.len(), 2);
    assert!(product.nicks.is_empty());
    let main = product
        .components
        .iter()
        .find(|c| c.covalently_closed)
        .unwrap();
    assert_eq!(main.paired_nucleotides, 48);
    assert_eq!(main.unpaired_nucleotides, 0);
    let closed = product
        .chains
        .iter()
        .filter(|c| matches!(c.ends, p::ChainEnds::CovalentlyCircular))
        .collect::<Vec<_>>();
    assert_eq!(closed.len(), 2);
    assert!(closed.iter().all(|c| c.sequence.len() == 24));
    assert!(closed.iter().any(|c| c.sequence.as_str() == TARGET));
    let free = product
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert_eq!(
        free.sources,
        vec![p::SourcePlacement {
            source_strand: 3,
            start: 6,
            end: 10,
            output_start: 0
        }]
    );
    assert_eq!(
        free.ends,
        p::ChainEnds::Linear {
            five_prime: p::FivePrime::Phosphate,
            three_prime: p::ThreePrime::Hydroxyl
        }
    );
    assert_eq!(product.released_activations.len(), 2);
}
#[test]
fn archived_post_flap_and_sealed_graph_continues_remaining_activation_without_losing_four_calls() {
    let (w, request, _) = fixture();
    let before = serde_json::to_vec(&w.documents).unwrap();
    let captured = history::prepare_capture(&w, &request).unwrap();
    assert_eq!(captured.receipt.summary.source_bases, 52);
    assert_eq!(captured.receipt.summary.remaining_activations.len(), 1);
    assert_eq!(captured.receipt.summary.released_activations.len(), 1);
    assert!(!captured.receipt.summary.incomplete);
    assert_eq!(
        serde_json::to_vec(&captured.workspace.documents).unwrap(),
        before
    );
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&captured.workspace).unwrap())
            .unwrap();
    let origin = history::documentary_capsule(&reopened, &captured.receipt.selection).unwrap();
    let fresh = history::replay(&origin).unwrap();
    assert_eq!(fresh.report().physical.steps.len(), 3);
    assert_eq!(
        fresh
            .report()
            .physical
            .product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
    assert!(
        matches!(&fresh.report().physical.steps[1],p::StepEvidence::Process(p::processing::ActionEvidence::FlapCleavage(e)) if e.outcome==p::Outcome::Applied)
    );
    assert!(
        matches!(&fresh.report().physical.steps[2],p::StepEvidence::Process(p::processing::ActionEvidence::NickSealing(e)) if e.outcome==p::Outcome::Applied)
    );
    let continuation = resumed(&reopened, &captured.receipt.selection);
    let report = t::simulate(&[], &continuation).unwrap();
    assert!(report.physical.all_requested_steps_applied);
    assert!(!report.incomplete);
    assert_eq!(report.physical.request.len(), 5);
    assert_eq!(report.preparation.retained_steps.len(), 3);
    assert_whole_closed_graph(&report.physical.product);
    let completed = history::prepare_capture(
        &reopened,
        &history::CaptureRequest {
            name: "Completed processed reaction, free flap retained".into(),
            bindings: vec![],
            request: continuation.clone(),
        },
    )
    .unwrap();
    let completed_reopened =
        workspace::load_workspace(&workspace::save_workspace(&completed.workspace).unwrap())
            .unwrap();
    let completed_origin =
        history::documentary_capsule(&completed_reopened, &completed.receipt.selection).unwrap();
    assert_whole_closed_graph(
        &history::replay(&completed_origin)
            .unwrap()
            .report()
            .physical
            .product,
    );
    let component = report
        .physical
        .product
        .components
        .iter()
        .position(|c| c.covalently_closed)
        .unwrap();
    let bridge = material::Request {
        preparation: continuation,
        selection: p::projection::Selection {
            source_batch_revision: report.physical.batch_revision.clone(),
            component,
            reference_seed: call(0, 0),
        },
        options: p::projection::Options {
            require_all_requested_steps_applied: true,
            ..Default::default()
        },
    };
    let material = material::inspect(&[], &bridge).unwrap();
    assert_eq!(material.duplex.as_ref().unwrap().top.as_str(), TARGET);
    assert!(material.tracked.is_some());
    assert_whole_closed_graph(&material.projection.physical.product);
    let target = Molecule::dna("Independent requirement", TARGET, Topology::Circular).unwrap();
    let plan = native::compile(
        &[],
        &target,
        None,
        &native::Task::ProcessActions { request: bridge },
    )
    .unwrap();
    let verified = verification::verify(&plan, &target);
    assert_eq!(
        verified.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        verified.findings
    );
    let wrong = Molecule::dna(
        "Independent wrong requirement",
        "AAGGGGCATCCCTTCACCATGACA",
        Topology::Circular,
    )
    .unwrap();
    let mut wrong_plan = plan.clone();
    wrong_plan.target_revision = wrong.revision();
    assert_eq!(
        verification::verify(&wrong_plan, &wrong).status,
        Status::Invalid
    );
    assert_eq!(
        serde_json::to_vec(&completed_reopened.documents).unwrap(),
        before
    );
}
#[test]
fn effective_new_end_chemistry_and_consumed_backbones_are_not_redeclared_as_activation() {
    let (_, mut request, inputs) = fixture();
    request.request.steps.truncate(2);
    request.request.steps.extend([
        transfer(call(3, 5), call(1, 0), p::Geometry::NickReseal),
        transfer(call(1, 17), call(3, 6), p::Geometry::Flush),
        seal(call(3, 5), call(1, 0)),
        transfer(call(1, 17), call(3, 0), p::Geometry::Flush),
        seal(call(2, 9), call(0, 0)),
        transfer(call(1, 17), call(3, 0), p::Geometry::Flush),
        cut(),
    ]);
    let report = t::simulate(&inputs, &request.request).unwrap();
    assert_eq!(
        report.physical.steps[2].outcome(),
        p::Outcome::BlockedUnderDeclaredModel,
        "newly hydrolyzed3′OH is not a TOPO activation"
    );
    assert_eq!(
        report.physical.steps[3].outcome(),
        p::Outcome::BlockedUnderDeclaredModel,
        "newly hydrolyzed5′P is not5′OH"
    );
    assert_eq!(report.physical.steps[4].outcome(), p::Outcome::Applied);
    assert_eq!(
        report.physical.steps[5].outcome(),
        p::Outcome::Applied,
        "failed5′P attempt must not consume the remaining activation"
    );
    assert_eq!(report.physical.steps[6].outcome(), p::Outcome::Applied);
    assert_eq!(
        report.physical.steps[7].outcome(),
        p::Outcome::BlockedUnderDeclaredModel
    );
    assert_eq!(
        report.physical.steps[8].outcome(),
        p::Outcome::BlockedUnderDeclaredModel
    );
    assert!(!report.physical.all_requested_steps_applied);
    assert!(report.incomplete);
    assert_whole_closed_graph(&report.physical.product);
    p::replay(
        &report.preparation.substrates,
        &report.physical.request,
        &request.request.engine_limits,
        &report.physical,
    )
    .unwrap();
}
#[test]
fn stale_projection_pins_record_state_and_duplicate_native_scope_refuse() {
    let (w, r, _) = fixture();
    let capture = history::prepare_capture(&w, &r).unwrap();
    let request = resumed(&capture.workspace, &capture.receipt.selection);
    let report = t::simulate(&[], &request).unwrap();
    let mut changed = request.clone();
    changed.steps.pop();
    assert!(
        material::inspect(
            &[],
            &material::Request {
                preparation: changed,
                selection: p::projection::Selection {
                    source_batch_revision: report.physical.batch_revision,
                    component: 0,
                    reference_seed: call(0, 0)
                },
                options: Default::default()
            }
        )
        .is_err()
    );
    let mut forged =
        history::documentary_capsule(&capture.workspace, &capture.receipt.selection).unwrap();
    forged.selection.state_revision = "0".repeat(64);
    assert!(history::replay(&forged).is_err());
    let mut duplicate = request.clone();
    duplicate
        .recipe
        .materials
        .push(duplicate.recipe.materials[0].clone());
    assert!(t::simulate(&[], &duplicate).is_err());
}
