//! Durable continuation reexecutes original declarations; saved molecule metadata
//! and documentary replay evidence are never deserialized into physical stock.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    chemistry::{EndTreatment, Strand, Terminus},
    cloning_outputs::{self, Commit, Options},
    history::DocumentHistory,
    methylation::MethylationState,
    modified_duplex::{self, CallState, Chemical, GenerationAssumptions},
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, Capsule, Selection},
    transfer::{self, Candidate},
    verification::{
        self, DesignPlan, Operation, OrientedSource, PhysicalProcessingEvidence, PlanInput, Source,
        Status,
    },
    workspace::Workspace,
};

fn binding(w: &Workspace, index: usize) -> EventBinding {
    let doc = &w.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn finalize_original(
    originals: Vec<Molecule>,
    target: Molecule,
    steps: Vec<Operation>,
    name: &str,
) -> Commit {
    let n = originals.len();
    let plan = DesignPlan {
        name: name.into(),
        target_revision: target.revision(),
        inputs: originals
            .iter()
            .enumerate()
            .map(|(i, m)| PlanInput {
                id: format!("src{i}"),
                revision: m.revision(),
                molecule: m.clone(),
            })
            .collect(),
        steps,
    };
    let mut models = originals;
    models.push(target);
    let w = Workspace::from_documents("Original source-backed run", models, n).unwrap();
    cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, n),
        (0..n)
            .map(|i| InputEventBinding {
                input_id: format!("src{i}"),
                event: binding(&w, i),
            })
            .collect(),
        &Options {
            name: name.into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap()
}
fn pcr_step(template: Source) -> Operation {
    Operation::Pcr {
        template,
        forward: Primer::new("Forward", "ACGACCCC", "").unwrap(),
        reverse: Primer::new("Reverse", "GTTACCCC", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn pcr() -> Commit {
    finalize_original(
        vec![Molecule::dna("PCR original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap()],
        Molecule::dna(
            "Independent PCR target",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap(),
        vec![pcr_step(Source::Input("src0".into()))],
        "Named saved PCR",
    )
}
fn modified_circle() -> Commit {
    finalize_original(
        vec![Molecule::dna("CpG original", "AAGCGGCCGCAA", Topology::Circular).unwrap()],
        Molecule::dna("Independent CpG target", "AAGCGGCCGCAA", Topology::Circular).unwrap(),
        vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("src0".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::Unmodified,
                    dcm: MethylationState::Unmodified,
                    cpg: MethylationState::FullyModified,
                },
            },
            Operation::CheckedDigest {
                selection: checked_selection(vec!["NotI".into()]),
                options: collagene_core::modified_duplex::checked_inspection_options(),
                template: Source::Step(0),
                fragment_index: 0,
            },
        ],
        "Named saved CpG circle",
    )
}
const NICKED: &str = "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC";
fn nicked_circle() -> Commit {
    finalize_original(
        vec![
            Molecule::dna("Left original", "TTGCAATCGTCGATGCTAGC", Topology::Linear).unwrap(),
            Molecule::dna("Right original", "ATCGTCGATGCTAGCGGTAC", Topology::Linear).unwrap(),
        ],
        Molecule::dna("Independent nicked target", NICKED, Topology::Circular).unwrap(),
        vec![
            Operation::Pcr {
                template: Source::Input("src0".into()),
                forward: Primer::new("Left forward", "TTGCAATC", "").unwrap(),
                reverse: Primer::new("Left reverse", "GCTAGCAT", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Pcr {
                template: Source::Input("src1".into()),
                forward: Primer::new("Right forward", "ATCGTCGA", "").unwrap(),
                reverse: Primer::new("Right reverse", "GTACCGCT", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(2),
                        reverse: false,
                    },
                ],
                circular: true,
                allow_nicked: true,
                options: Default::default(),
            },
        ],
        "Named saved nicked circle",
    )
}
fn selection(commit: &Commit, step: usize) -> Selection {
    Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: step,
        event: commit.receipt.construct.clone(),
    }
}
fn continuation(
    commit: &Commit,
    step: usize,
    target: Molecule,
    next: Operation,
) -> (Workspace, DesignPlan, EventBinding) {
    let origin = run_material::capsule(&commit.workspace, &selection(commit, step)).unwrap();
    origin.validate_documentary().unwrap();
    let saved = commit
        .workspace
        .head(&commit.receipt.construct.document_id)
        .unwrap()
        .clone();
    let plan = DesignPlan {
        name: "Fresh continuation of saved material".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "saved".into(),
            revision: saved.revision(),
            molecule: saved,
        }],
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin,
            },
            next,
        ],
    };
    let mut w = commit.workspace.clone();
    let document = DocumentHistory::new(
        format!("independent-continuation-target-{}", w.documents.len()),
        target,
    )
    .unwrap();
    let target = EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event().into(),
        revision: document.head_revision().into(),
    };
    document.contribute(&mut w).unwrap();
    w.active_document = target.document_id.clone();
    (w, plan, target)
}
fn origin_mut(plan: &mut DesignPlan) -> &mut Capsule {
    match &mut plan.steps[0] {
        Operation::ReconstructRunProduct { origin, .. } => origin,
        _ => panic!("missing replay step"),
    }
}
fn evidence(report: &verification::VerificationReport) -> &run_material::ReplayEvidence {
    report
        .physical_processing_reports
        .iter()
        .find_map(|step| match &step.evidence {
            PhysicalProcessingEvidence::RunProduct(evidence) => Some(evidence),
            _ => None,
        })
        .expect("complete typed replay evidence")
}

#[test]
fn saved_pcr_can_seed_new_pcr_after_binary_roundtrip_and_outer_run_finalization() {
    let original = pcr();
    let (w, plan, target) = continuation(
        &original,
        0,
        Molecule::dna(
            "Independent second PCR",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap(),
        pcr_step(Source::Step(0)),
    );
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let candidate = Candidate::new(plan.clone());
    let encoded = transfer::save_candidate(&candidate).unwrap();
    let decoded = transfer::load_candidate(&encoded).unwrap();
    assert_eq!(decoded.plan_revision, candidate.plan_revision);
    let report = verification::verify(&decoded.plan, w.head(&target.document_id).unwrap());
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.len(), 2);
    assert_eq!(report.products[1].sequence.as_str(), "ACGACCCCGGGGTAAC");
    let html =
        verification::export_protocol_html(&decoded.plan, w.head(&target.document_id).unwrap());
    assert!(html.contains("Use saved material"));
    assert!(html.contains(&original.receipt.run.run_id));
    let pdf = collagene_core::pdf::export_protocol_pdf(
        &decoded.plan,
        w.head(&target.document_id).unwrap(),
    )
    .unwrap();
    let printed = pdf_text(&pdf.bytes)
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    assert!(printed.contains("Reconstruct already-prepared saved material"));
    let oligos = collagene_core::plan_oligos::build(
        &decoded.plan,
        w.head(&target.document_id).unwrap(),
        &Default::default(),
    )
    .unwrap();
    assert_eq!(oligos.uses.len(), 2);
    assert!(oligos.uses.iter().all(|use_| use_.step_index == 1));
    let replay = evidence(&report);
    replay.validate_documentary().unwrap();
    assert_eq!(replay.origin.selection, selection(&original, 0));
    assert_eq!(replay.decode_report().unwrap().products.len(), 1);
    assert_eq!(replay.material_revision.len(), 64);
    let committed = cloning_outputs::finalize(
        &w,
        &decoded.plan,
        target,
        vec![InputEventBinding {
            input_id: "saved".into(),
            event: original.receipt.construct.clone(),
        }],
        &Options {
            name: "Second named PCR".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let run = run_lineage::load_run(&committed.workspace, &committed.receipt.run.run_id).unwrap();
    assert_eq!(run.inputs[0].event, original.receipt.construct);
    let historical = collagene_core::physical_run_evidence::load(&committed.workspace, &run.run_id)
        .unwrap()
        .unwrap();
    let saved_evidence = historical
        .steps
        .iter()
        .find_map(|step| match &step.evidence {
            PhysicalProcessingEvidence::RunProduct(saved) => Some(saved),
            _ => None,
        })
        .unwrap();
    saved_evidence.validate_documentary().unwrap();
    assert_eq!(
        saved_evidence.decode_report().unwrap().products[0]
            .sequence
            .as_str(),
        "ACGACCCCGGGGTAAC"
    );
    let fresh = run_material::reconstruct(&committed.workspace, &selection(&committed, 1)).unwrap();
    assert_eq!(fresh.tracked().duplex.left, Terminus::PCR);
    assert!(fresh.tracked().duplex.nicks.is_empty() && fresh.tracked().generation.is_empty());
    assert!(serde_json::to_value(&saved_evidence.origin).unwrap()["candidate"].is_string());
    let archive = collagene_core::workspace::save_workspace(&committed.workspace).unwrap();
    let loaded = collagene_core::workspace::load_workspace(&archive).unwrap();
    let restored = run_material::reconstruct(&loaded, &selection(&committed, 1)).unwrap();
    assert_eq!(restored.reference().sequence.as_str(), "ACGACCCCGGGGTAAC");
    assert_eq!(restored.tracked().duplex.left, Terminus::PCR);
    assert_eq!(restored.tracked().duplex.right, Terminus::PCR);
    assert!(restored.tracked().duplex.nicks.is_empty());
    assert!(restored.tracked().generation.is_empty());
    let restored_evidence = evidence(restored.report());
    restored_evidence.validate_documentary().unwrap();
    assert_eq!(restored_evidence.origin.selection, selection(&original, 0));
    assert_eq!(restored_evidence.decode_report().unwrap().products.len(), 1);
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}

#[test]
fn saved_modified_circle_remains_noti_blocked_in_durable_checked_digest() {
    let original = modified_circle();
    let (w, plan, target) = continuation(
        &original,
        1,
        Molecule::dna(
            "Independent still blocked target",
            "AAGCGGCCGCAA",
            Topology::Circular,
        )
        .unwrap(),
        Operation::CheckedDigest {
            selection: checked_selection(vec!["NotI".into()]),
            options: collagene_core::modified_duplex::checked_inspection_options(),
            template: Source::Step(0),
            fragment_index: 0,
        },
    );
    let report = verification::verify(&plan, w.head(&target.document_id).unwrap());
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products[1].sequence.as_str(), "AAGCGGCCGCAA");
    assert_eq!(report.products[1].topology, Topology::Circular);
    assert_eq!(evidence(&report).decode_report().unwrap().products.len(), 2);
    let committed = cloning_outputs::finalize(
        &w,
        &plan,
        target,
        vec![InputEventBinding {
            input_id: "saved".into(),
            event: original.receipt.construct.clone(),
        }],
        &Options {
            name: "Still CpG blocked".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let fresh = run_material::reconstruct(&committed.workspace, &selection(&committed, 1)).unwrap();
    for strand in [Strand::Top, Strand::Bottom] {
        assert_eq!(
            (0..12)
                .filter(
                    |&p| fresh.tracked().state_at(strand, Chemical::C5C, p).unwrap()
                        == Some(CallState::Modified)
                )
                .collect::<Vec<_>>(),
            vec![3, 7]
        );
    }
    assert_eq!(fresh.tracked().generation.len(), 1);
    let original_material =
        run_material::reconstruct(&original.workspace, &selection(&original, 1)).unwrap();
    assert_eq!(
        fresh.tracked().generation,
        original_material.tracked().generation
    );
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .known_modified,
        4
    );
}

#[test]
fn saved_nicked_circle_cannot_become_an_intact_pcr_template_or_lose_nicks() {
    let original = nicked_circle();
    let (w, mut plan, target) = continuation(
        &original,
        3,
        Molecule::dna(
            "Independent retained-nick circle",
            NICKED,
            Topology::Circular,
        )
        .unwrap(),
        Operation::Digest {
            template: Source::Step(0),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        },
    );
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let digest = plan.steps.pop().unwrap();
    let report = verification::verify(&plan, w.head(&target.document_id).unwrap());
    assert_eq!(report.products.len(), 1);
    assert_ne!(report.status, Status::Invalid, "{:?}", report.findings);
    let committed = cloning_outputs::finalize(
        &w,
        &plan,
        target.clone(),
        vec![InputEventBinding {
            input_id: "saved".into(),
            event: original.receipt.construct.clone(),
        }],
        &Options {
            name: "Reconstructed nicked circle".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let fresh = run_material::reconstruct(&committed.workspace, &selection(&committed, 0)).unwrap();
    assert_eq!(fresh.tracked().duplex.nicks.len(), 2);
    assert!(
        fresh
            .tracked()
            .duplex
            .nicks
            .iter()
            .any(|nick| nick.strand == Strand::Top && nick.position == 20)
    );
    assert!(
        fresh
            .tracked()
            .duplex
            .nicks
            .iter()
            .any(|nick| nick.strand == Strand::Bottom && nick.position == 0)
    );
    plan.steps.push(digest);
    let withheld_digest = verification::verify(&plan, w.head(&target.document_id).unwrap());
    assert_eq!(withheld_digest.status, Status::Incomplete);
    assert_eq!(withheld_digest.products.len(), 1);
    assert!(
        withheld_digest
            .findings
            .iter()
            .any(|finding| finding.contains("digest propagation of existing nicks"))
    );
    plan.steps[1] = Operation::Pcr {
        template: Source::Step(0),
        forward: Primer::new("Forward", "TTGCAATC", "").unwrap(),
        reverse: Primer::new("Reverse", "GTACCGCT", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let rejected = verification::verify(&plan, w.head(&target.document_id).unwrap());
    assert!(matches!(
        rejected.status,
        Status::Invalid | Status::Incomplete
    ));
    assert_eq!(rejected.products.len(), 1);
    assert!(
        rejected
            .findings
            .iter()
            .any(|finding| finding.to_ascii_lowercase().contains("nick"))
    );
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}

#[test]
fn capsule_target_source_product_checksum_and_schema_tampering_are_rejected() {
    let original = pcr();
    let (w, plan, target) = continuation(
        &original,
        0,
        Molecule::dna(
            "Independent continuation target",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap(),
        pcr_step(Source::Step(0)),
    );
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let mutations: Vec<Box<dyn Fn(&mut Capsule)>> = vec![
        Box::new(|origin| origin.target.sequence = collagene_core::Sequence::dna("TTTT").unwrap()),
        Box::new(|origin| origin.inputs[0].event.revision = "0".repeat(64)),
        Box::new(|origin| origin.products[0].revision = "0".repeat(64)),
        Box::new(|origin| {
            let last = origin.candidate.len() - 1;
            origin.candidate[last] ^= 0x40;
        }),
        Box::new(|origin| origin.schema_version = 0),
    ];
    for mutate in mutations {
        let mut changed = plan.clone();
        mutate(origin_mut(&mut changed));
        let result = verification::verify(&changed, w.head(&target.document_id).unwrap());
        assert!(
            matches!(result.status, Status::Invalid | Status::Incomplete),
            "{:?}",
            result.findings
        );
        assert!(
            result.products.is_empty(),
            "A rejected origin cannot fabricate any reconstructed product"
        );
    }
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}

#[test]
fn equal_reference_alias_is_not_the_capsules_exact_saved_source_event() {
    let original = pcr();
    let (mut w, plan, target) = continuation(
        &original,
        0,
        Molecule::dna(
            "Independent alias test",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap(),
        pcr_step(Source::Step(0)),
    );
    let alias = DocumentHistory::new(
        "unrelated-equal-saved-reference",
        plan.inputs[0].molecule.clone(),
    )
    .unwrap();
    let event = EventBinding {
        document_id: alias.id.clone(),
        event_id: alias.head_event().into(),
        revision: alias.head_revision().into(),
    };
    assert_eq!(event.revision, original.receipt.construct.revision);
    assert_ne!(event.document_id, original.receipt.construct.document_id);
    alias.contribute(&mut w).unwrap();
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    // Reconstruction alone checks chemistry; durable adoption also attests the
    // exact source event chosen by the host, without equal-DNA substitution.
    assert_ne!(
        verification::verify(&plan, w.head(&target.document_id).unwrap()).status,
        Status::Invalid
    );
    assert!(
        cloning_outputs::finalize(
            &w,
            &plan,
            target,
            vec![InputEventBinding {
                input_id: "saved".into(),
                event
            }],
            &Options {
                name: "Cannot adopt alias".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: true
            }
        )
        .is_err()
    );
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}

#[test]
fn saved_material_operation_and_report_require_candidate17_and_worker73() {
    use collagene_core::protocol::{self, AnalysisJob, AnalysisResult, Request};
    let original = pcr();
    let (w, plan, target) = continuation(
        &original,
        0,
        Molecule::dna(
            "Independent version target",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap(),
        pcr_step(Source::Step(0)),
    );
    let candidate = Candidate::new(plan.clone());
    let mut encoded = transfer::save_candidate(&candidate).unwrap();
    assert!(u16::from_le_bytes([encoded[4], encoded[5]]) >= 17);
    encoded[4..6].copy_from_slice(&16u16.to_le_bytes());
    assert!(transfer::load_candidate(&encoded).is_err());
    let request = Request {
        id: 902,
        revision: target.revision.clone(),
        job: AnalysisJob::Verify {
            target: w.head(&target.document_id).unwrap().clone(),
            plan,
        },
    };
    let encoded = protocol::encode_request(&request).unwrap();
    assert!(u16::from_le_bytes([encoded[4], encoded[5]]) >= 73);
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&72u16.to_le_bytes());
    if let Ok(response) = protocol::process_request(&old) {
        assert!(matches!(
            protocol::decode_response(&response).unwrap().result,
            AnalysisResult::Error(_)
        ));
    }
    let mut response = protocol::process_request(&encoded).unwrap();
    assert!(matches!(
        protocol::decode_response(&response).unwrap().result,
        AnalysisResult::Verification(_)
    ));
    response[4..6].copy_from_slice(&72u16.to_le_bytes());
    assert!(protocol::decode_response(&response).is_err());
}

#[test]
fn capsule_nesting_and_advertised_oversize_fail_before_unbounded_decoding() {
    let original = pcr();
    let mut current = run_material::capsule(&original.workspace, &selection(&original, 0)).unwrap();
    let saved = original
        .workspace
        .head(&original.receipt.construct.document_id)
        .unwrap()
        .clone();
    let mut accepted = 1;
    let mut rejected = false;
    for _ in 0..10 {
        let plan = DesignPlan {
            name: "One nested source-backed reconstruction".into(),
            target_revision: current.target.revision(),
            inputs: vec![PlanInput {
                id: "saved".into(),
                revision: saved.revision(),
                molecule: saved.clone(),
            }],
            steps: vec![Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin: current.clone(),
            }],
        };
        let Ok(candidate) = transfer::save_candidate(&Candidate::new(plan)) else {
            rejected = true;
            break;
        };
        let mut parent = current.clone();
        parent.candidate = candidate;
        parent.inputs = vec![InputEventBinding {
            input_id: "saved".into(),
            event: original.receipt.construct.clone(),
        }];
        if parent.validate_documentary().is_err() {
            rejected = true;
            break;
        }
        accepted += 1;
        current = parent;
    }
    assert!(rejected);
    assert_eq!(
        accepted, 8,
        "The documented eight-generation bound must allow the valid prefix and reject its next generation."
    );
    // An attacker can advertise an enormous nested CGP payload in a tiny
    // envelope; rejecting it must not require allocating that advertised body.
    let mut malformed =
        run_material::capsule(&original.workspace, &selection(&original, 0)).unwrap();
    malformed.candidate[6..10].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(malformed.candidate.len() < 64 * 1024);
    assert!(malformed.validate_documentary().is_err());
}

#[test]
fn context_pinned_processed_topo_origin_replays_before_pcr_instead_of_trusting_its_projection() {
    use collagene_core::{topo, topo_nick_sealing as sealing, topo_processed_product as processed};
    use topo::physical::{self as physical, nick_sealing, processed_projection, projection};
    let dna = ["AGGGGCATCCCTT", "AGGGATGCCCCTT", "CGAATGA", "CATTCGA"];
    let originals = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| {
            Molecule::dna(format!("Original TOPO strand{i}"), dna, Topology::Linear).unwrap()
        })
        .collect::<Vec<_>>();
    let inputs = originals
        .iter()
        .enumerate()
        .map(|(i, molecule)| topo::Input {
            role: format!("src{i}"),
            molecule: molecule.clone(),
        })
        .collect::<Vec<_>>();
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![topo::SourceInterval {
                source_role: format!("src{i}"),
                start: 0,
                length: dna.len(),
                reverse: false,
            }],
            five_prime: if i < 2 {
                physical::FivePrime::Phosphate
            } else {
                physical::FivePrime::Hydroxyl
            },
            three_prime: if i < 2 {
                physical::ThreePrime::Phosphotyrosyl {
                    enzyme: physical::Enzyme::VacciniaTypeIBTyr274,
                    activation_id: format!("activation{i}"),
                }
            } else {
                physical::ThreePrime::Hydroxyl
            },
        })
        .collect();
    let preparation=topo::Request {
        recipe:topo::GenerationRecipe{schema_version:1,materials:vec![topo::MaterialRecipe::Prepared {
            declaration_id:"Explicit original TA material".into(),assumption:"Source correspondence and vector activation are explicit operator hypotheses, not imported material authority".into(),strands,
            pairing:vec![physical::PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:0,length:12},physical::PairSpan{strand_a:2,start_a:0,strand_b:3,start_b:0,length:6}],
        }]},
        transfers:vec![physical::Transfer{donor:0,acceptor:2,geometry:physical::Geometry::ThreePrimeTA},physical::Transfer{donor:1,acceptor:3,geometry:physical::Geometry::ThreePrimeTA}],
        engine_limits:physical::Limits::default(),projection_options:projection::Options::default(),
    };
    let treatment = sealing::Request {
        preparation,
        sealing: nick_sealing::Request {
            ligase_assumption: "Explicit continuous-template3OH5P ligase hypothesis".into(),
            nicks: vec![
                nick_sealing::NickSelection {
                    three: physical::Call {
                        strand: 2,
                        index: 6,
                    },
                    five: physical::Call {
                        strand: 0,
                        index: 0,
                    },
                },
                nick_sealing::NickSelection {
                    three: physical::Call {
                        strand: 3,
                        index: 6,
                    },
                    five: physical::Call {
                        strand: 1,
                        index: 0,
                    },
                },
            ],
        },
    };
    let native = sealing::inspect(&inputs, &treatment).unwrap();
    let request = processed::Request {
        treatment,
        selection: processed_projection::Selection {
            sealing_context_revision: native.sealing.context_revision,
            component: 0,
            reference_seed: physical::Call {
                strand: 0,
                index: 0,
            },
        },
        options: projection::Options {
            require_all_requested_transfers_applied: true,
            require_known_terminal_chemistry: true,
            ..Default::default()
        },
    };
    let target = Molecule::dna(
        "Independent processed TA circle",
        "AGGGGCATCCCTTCGAATGA",
        Topology::Circular,
    )
    .unwrap();
    let original = finalize_original(
        originals,
        target,
        vec![Operation::ProcessedTopo {
            sources: (0..4)
                .map(|i| verification::TopoSource {
                    role: format!("src{i}"),
                    source: Source::Input(format!("src{i}")),
                })
                .collect(),
            request,
        }],
        "Named sealed TOPO circle",
    );
    let (w, plan, target) = continuation(
        &original,
        0,
        Molecule::dna(
            "Independent PCR from saved TOPO",
            "AGGGGCATCCCTTCGAATGA",
            Topology::Linear,
        )
        .unwrap(),
        Operation::Pcr {
            template: Source::Step(0),
            forward: Primer::new("TOPO forward", "AGGGGCAT", "").unwrap(),
            reverse: Primer::new("TOPO reverse", "TCATTCGA", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
    );
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let report = verification::verify(&plan, w.head(&target.document_id).unwrap());
    assert_ne!(report.status, Status::Invalid, "{:?}", report.findings);
    assert_eq!(report.products.len(), 2);
    assert_eq!(report.products[1].sequence.as_str(), "AGGGGCATCCCTTCGAATGA");
    let original_report = evidence(&report).decode_report().unwrap();
    assert_eq!(original_report.processed_step_reports.len(), 1);
    assert!(
        original_report.processed_step_reports[0]
            .evidence
            .duplex
            .as_ref()
            .unwrap()
            .nicks
            .is_empty()
    );
    // The literal same saved reference without explicit run reconstruction is
    // still a restricted projection and cannot silently become a PCR template.
    let mut detached = plan.clone();
    detached.steps = vec![Operation::Pcr {
        template: Source::Input("saved".into()),
        forward: Primer::new("TOPO forward", "AGGGGCAT", "").unwrap(),
        reverse: Primer::new("TOPO reverse", "TCATTCGA", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }];
    let detached_report = verification::verify(&detached, w.head(&target.document_id).unwrap());
    assert!(matches!(
        detached_report.status,
        Status::Invalid | Status::Incomplete
    ));
    assert!(detached_report.products.is_empty());
    let mut changed = plan.clone();
    let mut originating = transfer::load_candidate(&origin_mut(&mut changed).candidate).unwrap();
    let Operation::ProcessedTopo { request, .. } = &mut originating.plan.steps[0] else {
        panic!("expected processed origin");
    };
    request.selection.sealing_context_revision = "0".repeat(64);
    if let Ok(bytes) = transfer::save_candidate(&Candidate::new(originating.plan)) {
        origin_mut(&mut changed).candidate = bytes;
        let changed_report = verification::verify(&changed, w.head(&target.document_id).unwrap());
        assert!(matches!(
            changed_report.status,
            Status::Invalid | Status::Incomplete
        ));
        assert!(changed_report.products.is_empty());
    }
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}

// Independent PDF string/ToUnicode reader shared in approach with tests/pdf.rs.
fn literals(bytes: &[u8]) -> Vec<(Vec<u8>, usize)> {
    let mut result = vec![];
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'(' {
            i += 1;
            continue;
        }
        let mut out = vec![];
        let mut depth = 1;
        i += 1;
        while i < bytes.len() && depth > 0 {
            let b = bytes[i];
            i += 1;
            match b {
                b'\\' => {
                    if i >= bytes.len() {
                        break;
                    }
                    let escaped = bytes[i];
                    i += 1;
                    if escaped.is_ascii_digit() && escaped <= b'7' {
                        let mut n = escaped - b'0';
                        for _ in 0..2 {
                            if i < bytes.len() && (b'0'..=b'7').contains(&bytes[i]) {
                                n = n.wrapping_mul(8).wrapping_add(bytes[i] - b'0');
                                i += 1;
                            } else {
                                break;
                            }
                        }
                        out.push(n);
                    } else {
                        match escaped {
                            b'n' => out.push(b'\n'),
                            b'r' => out.push(b'\r'),
                            b't' => out.push(b'\t'),
                            b'b' => out.push(8),
                            b'f' => out.push(12),
                            b'\n' => {}
                            b'\r' => {
                                if bytes.get(i) == Some(&b'\n') {
                                    i += 1;
                                }
                            }
                            x => out.push(x),
                        }
                    }
                }
                b'(' => {
                    depth += 1;
                    out.push(b);
                }
                b')' => {
                    depth -= 1;
                    if depth > 0 {
                        out.push(b);
                    }
                }
                _ => out.push(b),
            }
        }
        if depth == 0 {
            result.push((out, i));
        }
    }
    result
}
fn pdf_text(bytes: &[u8]) -> String {
    let text = String::from_utf8_lossy(bytes);
    let mut map = std::collections::BTreeMap::<u16, String>::new();
    for line in text.lines() {
        let tokens: Vec<_> = line.split_whitespace().collect();
        if tokens.len() != 2
            || !tokens[0].starts_with('<')
            || tokens[0].len() != 6
            || !tokens[1].starts_with('<')
        {
            continue;
        }
        let Ok(cid) = u16::from_str_radix(&tokens[0][1..5], 16) else {
            continue;
        };
        let unicode = &tokens[1][1..tokens[1].len() - 1];
        if unicode.len() % 4 != 0 {
            continue;
        }
        let units: Option<Vec<_>> = (0..unicode.len())
            .step_by(4)
            .map(|i| u16::from_str_radix(&unicode[i..i + 4], 16).ok())
            .collect();
        if let Some(units) = units {
            if let Ok(s) = String::from_utf16(&units) {
                map.insert(cid, s);
            }
        }
    }
    let mut out = String::new();
    // Delimit real streams by their declared /Length, independent of embedded bytes.
    let mut start = 0;
    while let Some(relative) = bytes[start..].windows(8).position(|s| s == b"\nstream\n") {
        let content = start + relative + 8;
        let dictionary = String::from_utf8_lossy(&bytes[start..content]);
        let Some(length) = dictionary
            .rsplit("/Length ")
            .next()
            .and_then(|s| s.split_whitespace().next())
            .and_then(|s| s.parse::<usize>().ok())
        else {
            break;
        };
        if content + length > bytes.len() {
            break;
        }
        let stream = &bytes[content..content + length];
        if stream.windows(3).any(|s| s == b"BT\n") {
            for (s, end) in literals(stream) {
                if !stream[end..]
                    .iter()
                    .skip_while(|b| b.is_ascii_whitespace())
                    .copied()
                    .take(2)
                    .eq(b"Tj".iter().copied())
                {
                    continue;
                }
                for pair in s.chunks_exact(2) {
                    if let Some(t) = map.get(&u16::from_be_bytes([pair[0], pair[1]])) {
                        out.push_str(t);
                    }
                }
                out.push('\n');
            }
        }
        start = content + length;
    }
    out
}

fn checked_selection(names: Vec<String>) -> collagene_core::catalog_digest::CatalogSelection {
    let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
    let ids = names
        .iter()
        .map(|name| {
            catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.name.eq_ignore_ascii_case(name))
                .unwrap()
                .id
                .clone()
        })
        .collect::<Vec<_>>();
    collagene_core::catalog_digest::select(&catalog, &ids).unwrap()
}
