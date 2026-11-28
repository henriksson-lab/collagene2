//! Fresh typed processed-material science, downstream chemistry and independent intent.
const LINEAR_ECORI: &str = "AATTCAAGGGGCATCCCTTCG";
const ECORI_CIRCLE: &str = "AGGGGCATCCCTTCGAATTCA";
const CIRCLE: &str = "AGGGGCATCCCTTCGAATGA";
use collagene_core::{
    chemistry, modified_duplex, topo, topo_actions as adapter,
    verification::{self as v, DesignPlan, Operation, PlanInput, Source, Status, TopoSource},
};
use collagene_model::{Molecule, Topology};
use topo::physical::{self as p, processing as process, projection};
fn fixture(directional: bool) -> (Vec<topo::Input>, topo::Request) {
    let dna = if directional {
        [
            "AAGGGGCATCCCTT",
            "GGTGAAGGGATGCCCCTT",
            "CACCATGACG",
            "CGTCATGGTG",
        ]
    } else {
        ["AGGGGCATCCCTT", "AGGGATGCCCCTT", "CGAATGA", "CATTCGA"]
    };
    let inputs = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::Input {
            role: format!("original-{i}"),
            molecule: Molecule::dna(format!("source{i}"), dna, Topology::Linear).unwrap(),
        })
        .collect();
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![topo::SourceInterval {
                source_role: format!("original-{i}"),
                start: 0,
                length: dna.len(),
                reverse: false,
            }],
            five_prime: if i < 2 {
                p::FivePrime::Phosphate
            } else {
                p::FivePrime::Hydroxyl
            },
            three_prime: if i < 2 {
                p::ThreePrime::Phosphotyrosyl {
                    enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                    activation_id: format!("enzyme{i}"),
                }
            } else {
                p::ThreePrime::Hydroxyl
            },
        })
        .collect();
    let pairing = if directional {
        vec![
            p::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 4,
                length: 14,
            },
            p::PairSpan {
                strand_a: 2,
                start_a: 0,
                strand_b: 3,
                start_b: 0,
                length: 10,
            },
        ]
    } else {
        vec![
            p::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 0,
                length: 12,
            },
            p::PairSpan {
                strand_a: 2,
                start_a: 0,
                strand_b: 3,
                start_b: 0,
                length: 6,
            },
        ]
    };
    let preparation=topo::Request{
        recipe:topo::GenerationRecipe{schema_version:1,materials:vec![topo::MaterialRecipe::Prepared{declaration_id:"declared original material".into(),assumption:"Activation and supplied strand geometry are explicit hypotheses, never inferred from reference DNA".into(),strands,pairing}]},
        steps:vec![p::Step::Transfer(p::Transfer{donor: p::Call {strand:0,index:if directional { 13 } else { 12 }},acceptor:p::Call{strand:2,index:0},geometry:if directional{p::Geometry::DirectionalCacc}else{p::Geometry::ThreePrimeTA}}),p::Step::Transfer(p::Transfer{donor: p::Call {strand:1,index:if directional { 17 } else { 12 }},acceptor:p::Call{strand:3,index:0},geometry:if directional{p::Geometry::Flush}else{p::Geometry::ThreePrimeTA}})],
        engine_limits:p::Limits::default(),projection_options:projection::Options::default()
    };
    let mut preparation = preparation;
    for strand in [2, 3] {
        preparation
            .steps
            .push(p::Step::Process(process::Action::SealNick {
                selection: p::nick_sealing::NickSelection {
                    three: p::Call { strand, index: 6 },
                    five: p::Call {
                        strand: strand - 2,
                        index: 0,
                    },
                },
                ligase_assumption:
                    "Explicit ligase hypothesis: only templated3OH/5P nick sealing, no repair"
                        .into(),
            }));
    }
    (inputs, preparation)
}

fn selected(inputs: &[topo::Input], preparation: topo::Request) -> adapter::Request {
    let fresh = topo::simulate(inputs, &preparation).unwrap();
    adapter::Request {
        preparation,
        selection: projection::Selection {
            source_batch_revision: fresh.physical.batch_revision,
            component: 0,
            reference_seed: p::Call {
                strand: 0,
                index: 0,
            },
        },
        options: projection::Options {
            require_all_requested_steps_applied: true,
            require_known_terminal_chemistry: true,
            ..Default::default()
        },
    }
}

fn plan(inputs: &[topo::Input], request: adapter::Request, target: &Molecule) -> DesignPlan {
    DesignPlan {
        name: "Actual original TA material and independently specified product".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .map(|i| PlanInput {
                id: i.role.clone(),
                revision: i.molecule.revision(),
                molecule: i.molecule.clone(),
            })
            .collect(),
        steps: vec![Operation::TopoActions {
            sources: inputs
                .iter()
                .map(|i| TopoSource {
                    role: i.role.clone(),
                    source: Source::Input(i.role.clone()),
                })
                .collect(),
            request,
        }],
    }
}
fn digest() -> Operation {
    Operation::Digest {
        template: Source::Step(0),
        enzymes: vec!["EcoRI".into()],
        fragment_index: 0,
    }
}
fn ecorifixture() -> (Vec<topo::Input>, topo::Request) {
    let (mut inputs, mut treatment) = fixture(false);
    inputs[2].molecule =
        Molecule::dna("Original EcoRI insert top", "CGAATTCA", Topology::Linear).unwrap();
    inputs[3].molecule =
        Molecule::dna("Original EcoRI insert bottom", "GAATTCGA", Topology::Linear).unwrap();
    let topo::MaterialRecipe::Prepared {
        strands, pairing, ..
    } = &mut treatment.recipe.materials[0]
    else {
        panic!()
    };
    strands[2].segments[0].length = 8;
    strands[3].segments[0].length = 8;
    pairing[1].length = 7;
    for step in treatment.steps.iter_mut().skip(2) {
        let p::Step::Process(process::Action::SealNick { selection, .. }) = step else {
            panic!()
        };
        selection.three.index = 7;
    }
    (inputs, treatment)
}

#[test]
fn actual_sealed_ta_flows_to_ecori_digest_and_frozen_independent_target() {
    let (inputs, treatment) = ecorifixture();
    let request = selected(&inputs, treatment);
    let target = Molecule::dna(
        "Independently specified EcoRI linear output",
        LINEAR_ECORI,
        Topology::Linear,
    )
    .unwrap();
    let mut plan = plan(&inputs, request, &target);
    plan.steps.push(digest());
    let report = v::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.engine_version, v::ENGINE_VERSION);
    assert_eq!(report.products.len(), 2);
    assert_eq!(report.products[1].sequence.as_str(), LINEAR_ECORI);
    assert_eq!(report.physical_processing_reports.len(), 1);
    assert_eq!(report.physical_processing_reports[0].step_index, 0);

    let v::PhysicalProcessingEvidence::TopoActions(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("canonical action evidence")
    };
    assert_eq!(evidence.projection.physical.steps.len(), 4);
    assert_eq!(evidence.projection.physical.product.nicks.len(), 0);
    let d = evidence.duplex.as_ref().unwrap();
    assert_eq!(d.top.as_str(), ECORI_CIRCLE);
    let summary: chemistry::DuplexSummary =
        serde_json::from_str(&report.products[1].metadata["duplex_model"]).unwrap();
    assert_eq!(summary.bottom_start, 4);
    let (left, right) = summary.ends.unwrap();
    assert_eq!(left.overhang, "AATT");
    assert_eq!(right.overhang, "AATT");
    let modifications: modified_duplex::ModificationSummary =
        serde_json::from_str(&report.products[1].metadata["tracked_modifications"]).unwrap();
    assert!(modifications.unknown_eligible > 0);
    assert_eq!(modifications.known_unmodified, 0);
    assert_eq!(modifications.known_modified, 0);
    let html = v::export_protocol_html(&plan, &target);
    assert!(html.contains("TOPO"));
    assert!(html.contains("Explicit ligase hypothesis"));
    assert_eq!(evidence.projection.physical.product.nicks.len(), 0);
    assert!(html.contains("not isolation"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, report.status);
    assert_eq!(pdf.plan_revision, report.plan_revision);
    assert_eq!(pdf.target_revision, target.revision());
    let encoded = postcard::to_allocvec(&plan.steps[0]).unwrap();
    assert_eq!(encoded[0], 19);
    let candidate = collagene_core::transfer::Candidate::new(plan.clone());
    let bytes = collagene_core::transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        collagene_core::protocol::schema::CANDIDATE_SCHEMA_VERSION
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&12u16.to_le_bytes());
    assert!(collagene_core::transfer::load_candidate(&old).is_err());
    let request = collagene_core::protocol::request::Request::encode(
        &collagene_core::protocol::request::Request {
            id: 65,
            revision: target.revision(),
            job: collagene_core::protocol::schema::AnalysisJob::Verify {
                plan: plan.clone(),
                target: target.clone(),
            },
        },
    )
    .unwrap();
    let mut old = request.clone();
    old[4..6].copy_from_slice(&64u16.to_le_bytes());
    assert!(collagene_core::protocol::execution::process_request(&old).is_err());
    let response = collagene_core::protocol::execution::process_request(&request).unwrap();
    let mut old = response.clone();
    old[4..6].copy_from_slice(&64u16.to_le_bytes());
    assert!(collagene_core::protocol::response::Response::decode(&old).is_err());
    let collagene_core::protocol::schema::AnalysisResult::Verification(rebuilt) =
        collagene_core::protocol::response::Response::decode(&response)
            .unwrap()
            .result
    else {
        panic!("fresh typed verifier")
    };
    assert_eq!(rebuilt.physical_processing_reports.len(), 1);
}
#[test]
fn unsealed_nick_is_retained_and_downstream_digest_cannot_erase_it() {
    let (inputs, mut treatment) = fixture(false);
    treatment.steps.truncate(3);
    let request = selected(&inputs, treatment);
    let target = Molecule::dna("Independent nicked contour", CIRCLE, Topology::Circular).unwrap();
    let mut plan = plan(&inputs, request, &target);
    let report = v::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let v::PhysicalProcessingEvidence::TopoActions(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("canonical action evidence")
    };
    assert_eq!(evidence.duplex.as_ref().unwrap().nicks.len(), 1);
    assert_eq!(evidence.projection.physical.product.nicks.len(), 1);
    plan.steps.push(digest());
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.products.len(), 1);
    assert_eq!(
        {
            let v::PhysicalProcessingEvidence::TopoActions(e) =
                &report.physical_processing_reports[0].evidence
            else {
                panic!()
            };
            e
        }
        .duplex
        .as_ref()
        .unwrap()
        .nicks
        .len(),
        1
    );
    assert!(
        report
            .findings
            .iter()
            .any(|s| s.contains("digest propagation of existing nicks"))
    );
}
#[test]
fn directional_flap_withholds_product_but_preserves_full_new_typed_evidence() {
    let (inputs, mut treatment) = fixture(true);
    treatment.steps.truncate(2);
    let fresh = topo::simulate(&inputs, &treatment).unwrap();
    for nick in &fresh.physical.product.nicks {
        treatment
            .steps
            .push(p::Step::Process(process::Action::SealNick {
                selection: p::nick_sealing::NickSelection {
                    three: nick.three,
                    five: nick.five,
                },
                ligase_assumption: "Explicit ideal ligase, no flap hydrolysis".into(),
            }));
    }
    let request = selected(&inputs, treatment);
    let target = Molecule::dna(
        "Independently specified desired contour",
        CIRCLE,
        Topology::Circular,
    )
    .unwrap();
    let plan = plan(&inputs, request, &target);
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.products.is_empty());
    assert!(report.step_reports.iter().all(|s| s.topo.is_none()));
    let v::PhysicalProcessingEvidence::TopoActions(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("canonical action evidence")
    };
    assert!(evidence.tracked.is_none());
    assert!(evidence.product.is_none());
    let graph = &evidence.projection.physical.product;
    assert_eq!(
        graph.chains.iter().map(|c| c.sequence.len()).sum::<usize>(),
        52
    );
    assert!(
        graph
            .unpaired
            .iter()
            .any(|u| &graph.chains[u.chain].sequence.as_str()[u.start..u.end] == "GGTG")
    );
    let html = v::export_protocol_html(&plan, &target);
    assert!(html.contains("TOPO"));
    assert!(html.contains("Faithful selected duplex available: false"));
}
#[test]
fn actual_material_does_not_select_or_rewrite_an_independently_wrong_target() {
    let (inputs, treatment) = fixture(false);
    let request = selected(&inputs, treatment);
    let target = Molecule::dna(
        "Wrong independently specified circle",
        &format!("{CIRCLE}A"),
        Topology::Circular,
    )
    .unwrap();
    let plan = plan(&inputs, request, &target);
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.products[0].sequence.as_str(), CIRCLE);
    assert_eq!(report.physical_processing_reports.len(), 1);
    assert!(
        matches!(&report.physical_processing_reports[0].evidence,v::PhysicalProcessingEvidence::TopoActions(e) if e.tracked.is_some())
    );
    assert!(
        report
            .findings
            .iter()
            .any(|s| s.contains("does not match the target"))
    );
}
#[test]
fn portable_product_summaries_and_tracked_step_ancestry_cannot_restore_material() {
    let (inputs, treatment) = fixture(false);
    let request = selected(&inputs, treatment);
    let target = Molecule::dna("Independent circle", CIRCLE, Topology::Circular).unwrap();
    let first = plan(&inputs, request.clone(), &target);
    let report = v::verify(&first, &target);
    let imported = report.products[0].clone();
    assert!(chemistry::has_unmapped_physical_evidence(&imported));
    let imported_plan = DesignPlan {
        name: "Portable summary is no material".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "portable".into(),
            revision: imported.revision(),
            molecule: imported,
        }],
        steps: vec![Operation::Digest {
            template: Source::Input("portable".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        }],
    };
    let report = v::verify(&imported_plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.products.is_empty());
    let mut nested = first.clone();
    let Operation::TopoActions { sources, .. } = &mut nested.steps[0] else {
        panic!()
    };
    sources[0].source = Source::Step(0);
    nested.steps.insert(
        0,
        Operation::DeclareSourceModifications {
            template: Source::Input(inputs[0].role.clone()),
            assumptions: modified_duplex::GenerationAssumptions::default(),
        },
    );
    let report = v::verify(&nested, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.physical_processing_reports.is_empty());
    assert!(
        report
            .findings
            .iter()
            .any(|s| { s.contains("Prepared activation hypotheses require original references") })
    );
    // Even a source-reference shape that contains imported physical provenance
    // cannot be promoted by naming it as an original role in the new operation.
    let mut tainted = first;
    tainted.inputs[0]
        .molecule
        .metadata
        .insert("duplex_model".into(), "imported graph summary".into());
    tainted.inputs[0].revision = tainted.inputs[0].molecule.revision();
    let report = v::verify(&tainted, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.products.is_empty());
    assert!(
        report
            .findings
            .iter()
            .any(|s| s.contains("imported physical/reference-only summaries"))
    );
}
#[test]
fn treatment_and_nested_search_caps_reject_before_plan_revision_or_source_clones() {
    let (inputs, treatment) = fixture(false);
    let request = selected(&inputs, treatment);
    let target = Molecule::dna("Independent", CIRCLE, Topology::Circular).unwrap();
    let mut plan = plan(&inputs, request, &target);
    let Operation::TopoActions { request, .. } = &mut plan.steps[0] else {
        panic!()
    };
    request.options.limits.max_pairs = 0;
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert!(report.plan_revision.is_empty());
    assert!(report.products.is_empty());
    let Operation::TopoActions { request, .. } = &mut plan.steps[0] else {
        panic!()
    };
    request.options.limits.max_pairs = 50_000;
    request.preparation.engine_limits.max_steps = 3;
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert!(report.plan_revision.is_empty());
}
