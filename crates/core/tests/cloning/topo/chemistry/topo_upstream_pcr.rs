//! Literal PCR-copy composition; upstream chemistry remains distinct from new DNA.
use collagene_core::cloning::topo::physical::{self as p, processing as a, projection};
use collagene_core::{
    chemistry::{
        methylation::MethylationState,
        tracked_duplex::{GenerationAssumptions, ModificationSummary},
    },
    cloning::topo::{self, actions},
    documents::{
        transfer,
        workspace::{self, Workspace},
    },
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        physical_evidence as physical_run_evidence,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source, TopoSource},
    protocol,
    sequence::analysis::Primer,
    verification::{
        self as v,
        evidence::{PhysicalProcessingEvidence, Status},
    },
};
use collagene_model::{Molecule, Topology};
const TARGET: &str = "AAGGGGCATCCCTTCACCATGACG";
fn call(strand: usize, index: usize) -> p::Call {
    p::Call { strand, index }
}
fn fixture() -> (DesignPlan, Molecule, Vec<topo::Input>) {
    let inputs = vec![
        ("vt", "AAGGGGCATCCCTT"),
        ("vb", "GGTGAAGGGATGCCCCTT"),
        ("template", "CACCATGACG"),
    ]
    .into_iter()
    .map(|(role, dna)| topo::Input {
        role: role.into(),
        molecule: Molecule::dna(role, dna, Topology::Linear).unwrap(),
    })
    .collect::<Vec<_>>();
    let mut preparation = topo::Request {
        recipe: topo::GenerationRecipe {
            schema_version: 1,
            materials: vec![
                topo::MaterialRecipe::Prepared {
                    declaration_id: "Original activated vector".into(),
                    assumption: "Declared activated synthetic vector, not measured supplier state"
                        .into(),
                    strands: (0..2)
                        .map(|i| topo::PreparedStrand {
                            id: format!("vector{i}"),
                            segments: vec![topo::SourceInterval {
                                source_role: inputs[i].role.clone(),
                                start: 0,
                                length: inputs[i].molecule.sequence.len(),
                                reverse: false,
                            }],
                            five_prime: p::FivePrime::Phosphate,
                            three_prime: p::ThreePrime::Phosphotyrosyl {
                                enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                                activation_id: format!("activation{i}"),
                            },
                        })
                        .collect(),
                    pairing: vec![p::PairSpan {
                        strand_a: 0,
                        start_a: 0,
                        strand_b: 1,
                        start_b: 4,
                        length: 14,
                    }],
                },
                topo::MaterialRecipe::ExactPcr {
                    material_id: "New insert PCR strands".into(),
                    template_role: "template".into(),
                    forward: Primer::new("F", "CACCA", "").unwrap(),
                    reverse: Primer::new("R", "CGTCA", "").unwrap(),
                    tail: None,
                },
            ],
        },
        steps: vec![
            p::Step::Transfer(p::Transfer {
                donor: p::Call {
                    strand: 0,
                    index: 13,
                },
                acceptor: p::Call {
                    strand: 2,
                    index: 0,
                },
                geometry: p::Geometry::DirectionalCacc,
            }),
            p::Step::Transfer(p::Transfer {
                donor: p::Call {
                    strand: 1,
                    index: 17,
                },
                acceptor: p::Call {
                    strand: 3,
                    index: 0,
                },
                geometry: p::Geometry::Flush,
            }),
        ],
        engine_limits: p::Limits::default(),
        projection_options: projection::Options::default(),
    };
    let processing = vec![
        a::Action::CleaveThreePrimeFlap {
            bond: a::Bond {
                upstream: call(3, 5),
                downstream: call(3, 6),
            },
            displacement_step: 0,
            profile: a::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,
            assumption: "Declared ideal hydrolysis; not host repair or a named nuclease".into(),
        },
        a::Action::SealNick {
            selection: p::nick_sealing::NickSelection {
                three: call(2, 9),
                five: call(0, 0),
            },
            ligase_assumption: "Declared active ligase and cofactors".into(),
        },
        a::Action::SealNick {
            selection: p::nick_sealing::NickSelection {
                three: call(3, 5),
                five: call(1, 0),
            },
            ligase_assumption: "Declared active ligase and cofactors".into(),
        },
    ];
    // Select the actual source component before any target is supplied. Metadata
    // on the later template changes its identity, not these physical source calls.
    preparation
        .steps
        .extend(processing.into_iter().map(p::Step::Process));
    let inspected = topo::simulate(&inputs, &preparation).unwrap();
    let graph = &inspected.physical.product;
    let chain = graph
        .chains
        .iter()
        .position(|chain| {
            chain
                .sources
                .iter()
                .any(|placement| placement.source_strand == 0 && placement.start == 0)
        })
        .unwrap();
    let component = graph
        .components
        .iter()
        .position(|component| component.chains.contains(&chain))
        .unwrap();
    let request = actions::Request {
        preparation,
        selection: projection::Selection {
            source_batch_revision: inspected.physical.batch_revision,
            component,
            reference_seed: call(0, 0),
        },
        options: projection::Options {
            require_all_requested_steps_applied: true,
            require_known_terminal_chemistry: true,
            ..Default::default()
        },
    };
    let target =
        Molecule::dna("Independent literal24bp target", TARGET, Topology::Circular).unwrap();
    let plan = DesignPlan {
        name: "Upstream modified template, new PCR, accounted TOPO actions".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .map(|i| PlanInput {
                id: i.role.clone(),
                revision: i.molecule.revision(),
                molecule: i.molecule.clone(),
            })
            .collect(),
        steps: vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("template".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::Unknown,
                    dcm: MethylationState::Unknown,
                    cpg: MethylationState::FullyModified,
                },
            },
            Operation::TopoActions {
                sources: inputs
                    .iter()
                    .enumerate()
                    .map(|(i, input)| TopoSource {
                        role: input.role.clone(),
                        source: if i == 2 {
                            Source::Step(0)
                        } else {
                            Source::Input(input.role.clone())
                        },
                    })
                    .collect(),
                request,
            },
        ],
    };
    (plan, target, inputs)
}
fn evidence(r: &v::evidence::VerificationReport) -> &actions::Report {
    assert_eq!(r.physical_processing_reports.len(), 1);
    assert_eq!(r.physical_processing_reports[0].step_index, 1);
    let PhysicalProcessingEvidence::TopoActions(e) = &r.physical_processing_reports[0].evidence
    else {
        panic!("expected retained TopoActions evidence")
    };
    e
}
#[test]
fn modified_upstream_template_is_copied_into_new_material_and_portable_summary_is_not_authority() {
    let (plan, target, mut inputs) = fixture();
    let r = v::engine::verify(&plan, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products.len(), 2);
    assert_eq!(r.products[1].sequence.as_str(), TARGET);
    let before: ModificationSummary =
        serde_json::from_str(&r.products[0].metadata["tracked_modifications"]).unwrap();
    assert!(before.known_modified > 0);
    let e = evidence(&r);
    assert_eq!(e.projection.physical.source_bases, 52);
    assert_eq!(
        e.projection
            .physical
            .product
            .chains
            .iter()
            .map(|s| s.sequence.len())
            .sum::<usize>(),
        52
    );
    assert!(
        e.projection
            .physical
            .product
            .chains
            .iter()
            .any(|s| s.sequence.as_str() == "GGTG")
    );
    let actual = e.tracked.as_ref().unwrap();
    let after = collagene_core::chemistry::tracked_duplex::summary(actual).unwrap();
    assert_eq!(after.known_modified, 0);
    assert!(after.unknown_eligible > 0);
    assert!(actual.generation.is_empty());
    assert_eq!(
        e.sources
            .iter()
            .find(|s| s.role == "template")
            .unwrap()
            .revision,
        r.products[0].revision()
    );
    assert!(
        e.preparation
            .maps
            .iter()
            .any(|m| matches!(m.origin, topo::Origin::ExactPcrCopy { .. }))
    );
    // Same marked product supplied as an ordinary input remains rejected by the
    // public adapter and by the plan verifier; no metadata is stripped to pass it.
    inputs[2].molecule = r.products[0].clone();
    assert!(actions::inspect(&inputs, &e.request).is_err());
    let mut imported = plan.clone();
    imported.inputs[2].molecule = r.products[0].clone();
    imported.inputs[2].revision = r.products[0].revision();
    if let Operation::TopoActions { sources, .. } = &mut imported.steps[1] {
        sources[2].source = Source::Input("template".into());
    }
    assert_ne!(
        v::engine::verify(&imported, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let wrong = Molecule::dna(
        "Wrong frozen target",
        "AAGGGGCATCCCTTCACCATGACC",
        Topology::Circular,
    )
    .unwrap();
    let mut p = plan.clone();
    p.target_revision = wrong.revision();
    assert_eq!(v::engine::verify(&p, &wrong).status, Status::Invalid);
    assert!(v::output::export_protocol_html(&plan, &target).contains("PCR"));
    assert!(
        collagene_core::output::pdf::export_protocol_pdf(&plan, &target)
            .unwrap()
            .bytes
            .starts_with(b"%PDF-")
    );
}
#[test]
fn prepared_correspondence_unused_roles_and_nicked_or_staggered_templates_cannot_enter_the_copy_route()
 {
    let (plan, target, _) = fixture();
    let mut mixed = plan.clone();
    if let Operation::TopoActions { request, .. } = &mut mixed.steps[1] {
        if let topo::MaterialRecipe::Prepared { strands, .. } =
            &mut request.preparation.recipe.materials[0]
        {
            strands[0].segments[0].source_role = "template".into();
            strands[0].segments[0].length = 10;
        }
    }
    let r = v::engine::verify(&mixed, &target);
    assert_ne!(r.status, Status::ConsistentUnderAssumptions);
    assert!(r.physical_processing_reports.is_empty());
    let mut unused = plan.clone();
    if let Operation::TopoActions { sources, .. } = &mut unused.steps[1] {
        sources.push(TopoSource {
            role: "unused-upstream".into(),
            source: Source::Step(0),
        });
    }
    assert_ne!(
        v::engine::verify(&unused, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let pcr = Operation::Pcr {
        template: Source::Input("template".into()),
        forward: Primer::new("Fresh F", "CACCA", "").unwrap(),
        reverse: Primer::new("Fresh R", "CGTCA", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let mut staggered = plan.clone();
    staggered.steps[0] = pcr.clone();
    staggered.steps.insert(
        1,
        Operation::Tail3Prime {
            template: Source::Step(0),
            base: collagene_core::cloning::ligation::tailing::TailingBase::A,
        },
    );
    if let Operation::TopoActions { sources, .. } = &mut staggered.steps[2] {
        sources[2].source = Source::Step(1);
    }
    let r = v::engine::verify(&staggered, &target);
    assert_eq!(r.products.len(), 2, "{:?}", r.findings);
    assert_ne!(r.status, Status::ConsistentUnderAssumptions);
    assert!(r.physical_processing_reports.is_empty());
    assert!(r.findings.iter().any(|f| f.contains("nicked or staggered")));
    let mut nicked = plan.clone();
    nicked.steps[0] = Operation::DeclareSourceTermini {
        template: Source::Input("template".into()),
        left: collagene_core::chemistry::Terminus::CLEAVED,
        right: collagene_core::chemistry::Terminus::PCR,
    };
    nicked.steps.insert(
        1,
        Operation::Ligate {
            fragments: vec![collagene_core::planning::plan::OrientedSource {
                source: Source::Step(0),
                reverse: false,
            }],
            circular: true,
            allow_nicked: true,
            options: Default::default(),
        },
    );
    if let Operation::TopoActions { sources, .. } = &mut nicked.steps[2] {
        sources[2].source = Source::Step(1);
    }
    let r = v::engine::verify(&nicked, &target);
    assert_eq!(r.products.len(), 2, "{:?}", r.findings);
    assert_ne!(r.status, Status::ConsistentUnderAssumptions);
    // The successful upstream ligation retains its actual nick/alternative
    // evidence even though TOPO refuses that template. No TOPO PCR product or
    // physical transfer graph is created for the rejected third step.
    assert_eq!(r.physical_processing_reports.len(), 1);
    let evidence = &r.physical_processing_reports[0];
    assert_eq!(evidence.step_index, 1);
    let v::evidence::PhysicalProcessingEvidence::LigationOutcomes(ligation) = &evidence.evidence
    else {
        panic!("upstream ligation evidence")
    };
    assert_eq!(ligation.report.intended.duplex.nicks.len(), 1);
    assert!(!r.physical_processing_reports.iter().any(|step| matches!(
        &step.evidence,
        v::evidence::PhysicalProcessingEvidence::TopoActions(_)
    )));
    assert!(r.findings.iter().any(|f| f.contains("nicked or staggered")));
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn binary_candidate_and_retained_history_pin_the_exact_upstream_product_event() {
    let (plan, target, inputs) = fixture();
    let candidate = transfer::Candidate::new(plan.clone());
    let mut bytes = transfer::save_candidate(&candidate).unwrap();
    transfer::load_candidate(&bytes).unwrap();
    bytes[4..6].copy_from_slice(&14u16.to_le_bytes());
    assert!(transfer::load_candidate(&bytes).is_err());
    let mut bytes = protocol::request::Request::encode(&protocol::request::Request {
        id: 101,
        revision: target.revision(),
        job: protocol::schema::AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
    })
    .unwrap();
    let reply = protocol::response::Response::decode(
        &protocol::execution::process_request(&bytes).unwrap(),
    )
    .unwrap();
    let protocol::schema::AnalysisResult::Verification(r) = reply.result else {
        panic!("verification")
    };
    assert_eq!(r.status, Status::ConsistentUnderAssumptions);
    bytes[4..6].copy_from_slice(&66u16.to_le_bytes());
    assert!(protocol::execution::process_request(&bytes).is_err());
    let mut docs = inputs
        .iter()
        .map(|i| i.molecule.clone())
        .collect::<Vec<_>>();
    docs.push(target.clone());
    let w = Workspace::from_documents("Upstream source event fixture", docs, 3).unwrap();
    let bindings = inputs
        .iter()
        .enumerate()
        .map(|(i, input)| InputEventBinding {
            input_id: input.role.clone(),
            event: event(&w, i),
        })
        .collect();
    let committed = run_lineage::record_run(&w, &plan, event(&w, 3), bindings, false).unwrap();
    let loaded =
        workspace::load_workspace(&workspace::save_workspace(&committed.workspace).unwrap())
            .unwrap();
    let record = run_lineage::load_run(&loaded, &committed.receipt.run_id).unwrap();
    assert_eq!(record.plan_schema, workspace::PLAN_SCHEMA_VERSION);
    let mut older_record = record.clone();
    older_record.plan_schema = 12;
    assert!(run_lineage::save_run_record(&older_record).is_err());
    let mut historical_request = protocol::request::Request::encode(&protocol::request::Request {
        id: 102,
        revision: run_lineage::workspace_fingerprint(&loaded).unwrap(),
        job: protocol::schema::AnalysisJob::LoadRunEvidence {
            workspace: loaded.clone(),
            run_id: record.run_id.clone(),
        },
    })
    .unwrap();
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(&historical_request).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        response.result,
        protocol::schema::AnalysisResult::HistoricalRunEvidence(_)
    ));
    historical_request[4..6].copy_from_slice(&66u16.to_le_bytes());
    assert!(protocol::execution::process_request(&historical_request).is_err());
    let saved = physical_run_evidence::load(&loaded, &record.run_id)
        .unwrap()
        .unwrap();
    let PhysicalProcessingEvidence::TopoActions(e) = &saved.steps[0].evidence else {
        panic!("expected retained TopoActions evidence")
    };
    assert_eq!(
        e.sources
            .iter()
            .find(|s| s.role == "template")
            .unwrap()
            .revision,
        record.products[0].revision
    );
    assert_eq!(record.products[0].step_index, 0);
    assert_eq!(
        serde_json::to_value(&saved.products).unwrap(),
        serde_json::to_value(&record.products).unwrap()
    );
    let replay = run_lineage::replay_run(&loaded, &record.run_id).unwrap();
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(replay.products_match);
    let bytes = collagene_core::material::evidence::download(&loaded, &record.run_id)
        .unwrap()
        .bytes
        .unwrap();
    assert!(
        String::from_utf8(bytes)
            .unwrap()
            .contains("scientific-run-bundle-v1")
    );
    // Re-keying a historical sidecar after changing its upstream identity does
    // not transplant it onto an equal-looking source or another product event.
    let mut forged = loaded.clone();
    let (hash, mut asset) = forged
        .attachments
        .iter()
        .find(|(_, a)| a.kind == physical_run_evidence::KIND)
        .map(|(h, a)| (h.clone(), a.clone()))
        .unwrap();
    let mut changed = saved;
    let PhysicalProcessingEvidence::TopoActions(e) = &mut changed.steps[0].evidence else {
        panic!("expected retained TopoActions evidence")
    };
    e.sources
        .iter_mut()
        .find(|s| s.role == "template")
        .unwrap()
        .revision = target.revision();
    asset.bytes = serde_json::to_vec(&changed).unwrap();
    forged.attachments.remove(&hash);
    forged
        .attachments
        .insert(blake3::hash(&asset.bytes).to_hex().to_string(), asset);
    assert!(forged.validate().is_err());
}
