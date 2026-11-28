//! Intended target verification remains separate from retained alternative paths.
use collagene_core::{
    chemistry::{EndTreatment, Terminus},
    cloning::ligation::outcomes::{self as ligation_outcomes, Limit},
    documents::workspace::{self, Workspace},
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs, physical_evidence as physical_run_evidence,
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput},
    sequence::analysis::Primer,
    verification::{
        self,
        evidence::{PhysicalProcessingEvidence, Status},
    },
};
use collagene_model::{Molecule, Sequence, Topology};
const INTENDED: &str = "AGTCCGGACTGA";
const ALTERNATIVE: &str = "AGTCCGTCAGTC";
fn fixture(target: &str) -> (Vec<Molecule>, DesignPlan) {
    let a = Molecule::dna("Original vector stock", "AGTCCG", Topology::Linear).unwrap();
    let b = Molecule::dna("Original insert stock", "GACTGA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent explicitly declared circular target",
        target,
        Topology::Circular,
    )
    .unwrap();
    let inputs = [("a", a.clone()), ("b", b.clone())]
        .into_iter()
        .map(|(id, molecule)| PlanInput {
            id: id.into(),
            revision: molecule.revision(),
            molecule,
        })
        .collect();
    let steps = vec![
        Operation::Pcr {
            template: collagene_core::planning::plan::Source::Input("a".into()),
            forward: Primer::new("Vector F", "AGT", "").unwrap(),
            reverse: Primer::new("Vector R", "CGG", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        Operation::TreatEnds {
            template: collagene_core::planning::plan::Source::Step(0),
            treatment: EndTreatment::Phosphorylate,
        },
        Operation::Pcr {
            template: collagene_core::planning::plan::Source::Input("b".into()),
            forward: Primer::new("Insert F", "GAC", "").unwrap(),
            reverse: Primer::new("Insert R", "TCA", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        Operation::TreatEnds {
            template: collagene_core::planning::plan::Source::Step(2),
            treatment: EndTreatment::Phosphorylate,
        },
        ligation_outcomes::operation(
            vec![
                OrientedSource {
                    source: collagene_core::planning::plan::Source::Step(1),
                    reverse: false,
                },
                OrientedSource {
                    source: collagene_core::planning::plan::Source::Step(3),
                    reverse: false,
                },
            ],
            true,
            false,
        ),
    ];
    let plan = DesignPlan {
        standard: None,
        name: "Intended blunt ligation with retained orientation alternative".into(),
        target_revision: target.revision(),
        inputs,
        steps,
    };
    (vec![a, b, target], plan)
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn intended_and_alternative_stay_extractable_after_real_run_capture_reload_and_native_replay() {
    let (docs, plan) = fixture(INTENDED);
    let w = Workspace::from_documents("Original sources and independent intent", docs, 2).unwrap();
    let report = verification::engine::verify(&plan, w.head(&w.documents[2].id).unwrap());
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report
            .product(4)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        INTENDED
    );
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: event(&w, 0),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: event(&w, 1),
            },
        ],
        &Default::default(),
    )
    .unwrap();
    let restored =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let historical = physical_run_evidence::load(&restored, &saved.receipt.run.run_id)
        .unwrap()
        .unwrap();
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) = &historical.steps[0].evidence
    else {
        panic!("missing retained outcomes")
    };
    assert_eq!(evidence.report.intended.sequence.as_str(), INTENDED);
    assert_eq!(
        evidence.report.alternatives[0].sequence.as_str(),
        ALTERNATIVE
    );
    assert!(!evidence.report.intended.path[1].reverse);
    assert!(evidence.report.alternatives[0].path[1].reverse);
    evidence.validate_documentary().unwrap();
    assert!(evidence.report.fasta().unwrap().contains(ALTERNATIVE));
    // A PCR job without the current ligation capability cannot hide the operation inside its opaque saved-origin
    // candidate, even when the enclosing selector itself predates worker92.
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id.clone(),
        step_index: 4,
        event: saved.receipt.construct.clone(),
    };
    let capsule = collagene_core::material::capsule::capsule(&restored, &selection).unwrap();
    let reference = restored.revisions[&selection.event.revision].clone();
    let forward = Primer::new("Current F", "AGT", "").unwrap();
    let reverse = Primer::new("Current R", "TCA", "").unwrap();
    let origin = run_material::MaterialOrigin::Run(capsule);
    let revision = collagene_core::cloning::pcr::material::revision(
        &reference,
        &origin,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    let request = collagene_core::protocol::request::Request {
        id: 172,
        revision,
        job: collagene_core::protocol::schema::AnalysisJob::MaterialPcrDocuments {
            molecule: reference,
            origin,
            forward,
            reverse,
            template_strand: collagene_core::chemistry::Strand::Top,
        },
    };
    let bytes = collagene_core::protocol::request::Request::encode(&request).unwrap();
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&91u16.to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(collagene_core::protocol::execution::process_request(&old).is_err());
    let response = collagene_core::protocol::execution::process_request(&bytes).unwrap();
    let response = collagene_core::protocol::response::Response::decode(&response).unwrap();
    assert!(matches!(
        response.result,
        collagene_core::protocol::schema::AnalysisResult::MaterialPcrDocuments(_)
    ));
    let fresh = run_material::reconstruct(
        &restored,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id,
            step_index: 4,
            event: saved.receipt.construct,
        },
    )
    .unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), INTENDED);
    assert!(fresh.tracked().duplex.nicks.is_empty());
    let html =
        verification::output::export_protocol_html(&plan, w.head(&w.documents[2].id).unwrap());
    assert!(
        html.contains(INTENDED)
            && html.contains(ALTERNATIVE)
            && html.contains("alternative_1")
            && html.contains("Full typed outcome")
    );
    // A new attachment hash alone cannot replace alternative evidence under an
    // unchanged product/run identity, even if the forged duplex is paired.
    let mut altered = restored;
    let key = altered
        .attachments
        .iter()
        .find(|(_, a)| a.kind == physical_run_evidence::KIND)
        .unwrap()
        .0
        .clone();
    let mut asset = altered.attachments.remove(&key).unwrap();
    let mut data = historical;
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) = &mut data.steps[0].evidence else {
        panic!("outcomes")
    };
    evidence.report.alternatives[0].sequence = Sequence::dna("CCCCCCCCCCCC").unwrap();
    evidence.report.alternatives[0].duplex.top = Sequence::dna("CCCCCCCCCCCC").unwrap();
    evidence.report.alternatives[0].duplex.bottom = Sequence::dna("GGGGGGGGGGGG").unwrap();
    asset.bytes = serde_json::to_vec(&data).unwrap();
    altered
        .attachments
        .insert(blake3::hash(&asset.bytes).to_hex().to_string(), asset);
    assert!(altered.validate().is_err());
}
#[test]
fn matching_alternative_never_rescues_wrong_intended_target_and_capped_search_does_not_block_capture()
 {
    let (docs, plan) = fixture(ALTERNATIVE);
    let report = verification::engine::verify(&plan, &docs[2]);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(
        report
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        INTENDED
    );
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("outcomes")
    };
    assert_eq!(
        evidence.report.alternatives[0].sequence.as_str(),
        ALTERNATIVE
    );
    let (docs, mut plan) = fixture(INTENDED);
    let Operation::Ligate { options, .. } = plan.steps.last_mut().unwrap() else {
        panic!("new default")
    };
    options.max_products = 1;
    let w = Workspace::from_documents("Capped alternatives, complete intended material", docs, 2)
        .unwrap();
    let report = verification::engine::verify(&plan, w.head(&w.documents[2].id).unwrap());
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("outcomes")
    };
    assert!(evidence.report.limits.contains(&Limit::Products));
    assert!(!evidence.report.enumeration_complete);
    assert_eq!(evidence.report.intended.sequence.as_str(), INTENDED);
    assert!(
        cloning_outputs::finalize(
            &w,
            &plan,
            event(&w, 2),
            vec![
                InputEventBinding {
                    input_id: "a".into(),
                    event: event(&w, 0)
                },
                InputEventBinding {
                    input_id: "b".into(),
                    event: event(&w, 1)
                }
            ],
            &Default::default()
        )
        .is_ok()
    );
}
#[test]
fn large_recipe_retains_intended_outcomes_and_native_history_without_enumeration() {
    use collagene_core::cloning::{
        assembly::recipe::{FragmentSpec, Preparation},
        ligation::recipe as ligation_recipe,
    };
    let original = Molecule::dna("Original stock", "AG", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent33-copy declared intent",
        &"AG".repeat(33),
        Topology::Circular,
    )
    .unwrap();
    let inputs = vec![PlanInput {
        id: "stock".into(),
        revision: original.revision(),
        molecule: original,
    }];
    let recipe = ligation_recipe::Recipe {
        standard: None,
        name: "Explicit large ordered path".into(),
        fragments: (0..33)
            .map(|i| ligation_recipe::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: format!("row-{i}"),
                    input_id: "stock".into(),
                    preparation: Preparation::DeclaredDuplex {
                        left: Terminus::CLEAVED,
                        right: Terminus::CLEAVED,
                    },
                    reverse: false,
                },
                treatment: None,
            })
            .collect(),
        circular: true,
        allow_nicked: false,
    };
    let preview = ligation_recipe::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert!(matches!(
        preview.compiled.plan.steps.last(),
        Some(Operation::Ligate { .. })
    ));
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence,
        target.sequence
    );
    assert_eq!(preview.verification.physical_processing_reports.len(), 1);
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) =
        &preview.verification.physical_processing_reports[0].evidence
    else {
        panic!("every ligation retains outcomes");
    };
    assert_eq!(evidence.report.intended.sequence, target.sequence);
    assert!(evidence.report.alternatives.is_empty());
    assert!(
        evidence
            .report
            .limits
            .contains(&Limit::EnumerationFragments)
    );
    assert!(!evidence.report.enumeration_complete);
    let workspace = Workspace::from_documents(
        "Large intended ligation",
        vec![inputs[0].molecule.clone(), target.clone()],
        1,
    )
    .unwrap();
    let captured = collagene_core::material::lineage::record_run(
        &workspace,
        &preview.compiled.plan,
        event(&workspace, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&workspace, 0),
        }],
        false,
    )
    .unwrap();
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&captured.workspace).unwrap())
            .unwrap();
    let stored = physical_run_evidence::load(&reopened, &captured.receipt.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(stored.steps.len(), 1);
    let PhysicalProcessingEvidence::LigationOutcomes(stored) = &stored.steps[0].evidence else {
        panic!("missing persisted outcomes");
    };
    assert_eq!(stored.report.intended.sequence, target.sequence);
    assert!(stored.report.limits.contains(&Limit::EnumerationFragments));
    preview.validate_binding(&preview.context_revision).unwrap();
}
