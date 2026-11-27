//! Independent directional source-call, whole-mixture and final-target oracles.
use collagene_core::{
    Molecule, Topology, chemistry, modified_duplex, protocol, topo, topo_actions as actions,
    transfer,
    verification::{
        self as v, DesignPlan, Operation, PhysicalProcessingEvidence, PlanInput, Source, Status,
        TopoSource,
    },
};
use topo::physical::{self as p, processing as process, processing_projection as pp, projection};
const ORIGINAL_MAIN: &str = "AAGGGGCATCCCTTCACCATGACG";
const ECO_MAIN: &str = "AAGGGGCATCCCTTCACCGAATTCA";
const ECO_TARGET: &str = "AATTCAAAGGGGCATCCCTTCACCG";
fn call(strand: usize, index: usize) -> p::Call {
    p::Call { strand, index }
}
fn cut(anchor: usize) -> process::Action {
    process::Action::CleaveThreePrimeFlap{bond:process::Bond{upstream:call(3,anchor),downstream:call(3,anchor+1)},displacement_step:0,profile:process::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,assumption:"Explicit ideal hydrolytic cut of this anchored four-call flap; no named nuclease, host repair or kit claim".into()}
}
fn seal(three: p::Call, five: p::Call) -> process::Action {
    process::Action::SealNick {
        selection: p::nick_sealing::NickSelection { three, five },
        ligase_assumption:
            "Explicit active ligase/cofactor hypothesis for continuous-template3OH5P nicks".into(),
    }
}
fn fixture(site: bool) -> (Vec<topo::Input>, actions::ActionsRequest) {
    let dna = if site {
        [
            "AAGGGGCATCCCTT",
            "GGTGAAGGGATGCCCCTT",
            "CACCGAATTCA",
            "TGAATTCGGTG",
        ]
    } else {
        [
            "AAGGGGCATCCCTT",
            "GGTGAAGGGATGCCCCTT",
            "CACCATGACG",
            "CGTCATGGTG",
        ]
    };
    let inputs = dna
        .iter()
        .enumerate()
        .map(|(i, s)| topo::Input {
            role: format!("original-{i}"),
            molecule: Molecule::dna(
                format!("Original directional strand{i}"),
                s,
                Topology::Linear,
            )
            .unwrap(),
        })
        .collect();
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, s)| topo::PreparedStrand {
            id: format!("source{i}"),
            segments: vec![topo::SourceInterval {
                source_role: format!("original-{i}"),
                start: 0,
                length: s.len(),
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
    let preparation=topo::Request{recipe:topo::GenerationRecipe{schema_version:1,materials:vec![topo::MaterialRecipe::Prepared{declaration_id:"Original directional material declaration".into(),assumption:"Exact source correspondence, activation and terminal state are operator hypotheses, never inferred supplier evidence".into(),strands,pairing:vec![p::PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:4,length:14},p::PairSpan{strand_a:2,start_a:0,strand_b:3,start_b:0,length:dna[2].len()}]}]},transfers:vec![p::Transfer{donor:0,acceptor:2,geometry:p::Geometry::DirectionalCacc},p::Transfer{donor:1,acceptor:3,geometry:p::Geometry::Flush}],engine_limits:p::Limits::default(),projection_options:projection::Options::default()};
    let (anchor, last) = if site { (6, 10) } else { (5, 9) };
    let processing = process::Request {
        actions: vec![
            cut(anchor),
            seal(call(2, last), call(0, 0)),
            seal(call(3, anchor), call(1, 0)),
        ],
        limits: process::ProcessingLimits::default(),
    };
    (
        inputs,
        actions::ActionsRequest {
            preparation,
            processing,
        },
    )
}
fn select(inputs: &[topo::Input], request: actions::ActionsRequest) -> actions::Request {
    // Discover actual source component independently of any desired target.
    let fresh = actions::inspect_actions(inputs, &request).unwrap();
    let graph = &fresh.processing.product;
    let chain = graph
        .chains
        .iter()
        .position(|c| {
            c.sources
                .iter()
                .any(|m| m.source_strand == 0 && m.start == 0)
        })
        .unwrap();
    let component = graph
        .components
        .iter()
        .position(|c| c.chains.contains(&chain))
        .unwrap();
    actions::Request {
        preparation: request.preparation,
        processing: request.processing,
        selection: pp::Selection {
            processing_context_revision: fresh.processing.context_revision,
            component,
            reference_seed: call(0, 0),
        },
        options: projection::Options {
            require_all_requested_transfers_applied: true,
            require_known_terminal_chemistry: true,
            ..Default::default()
        },
    }
}
fn plan(inputs: &[topo::Input], request: actions::Request, target: &Molecule) -> DesignPlan {
    DesignPlan {
        name: "Directional accounted flap actions with independent intent".into(),
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
fn evidence(report: &v::VerificationReport) -> &actions::Report {
    assert_eq!(report.physical_processing_reports.len(), 1);
    assert_eq!(report.physical_processing_reports[0].step_index, 0);
    let PhysicalProcessingEvidence::TopoActions(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("expected retained TopoActions evidence")
    };
    evidence
}
#[test]
fn actual54calls_to50paired_plus4free_then_ecori_digest_matches_frozen25bp_target() {
    let (inputs, r) = fixture(true);
    let r = select(&inputs, r);
    let target = Molecule::dna(
        "Independently specified EcoRI25bp product",
        ECO_TARGET,
        Topology::Linear,
    )
    .unwrap();
    let mut p = plan(&inputs, r, &target);
    p.steps.push(digest());
    let report = v::verify(&p, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.engine_version, v::ENGINE_VERSION);
    assert!(report.processed_step_reports.is_empty());
    assert!(report.step_reports.is_empty());
    assert_eq!(report.products.len(), 2);
    let e = evidence(&report);
    let graph = &e.projection.processing.product;
    assert_eq!(e.projection.processing.source_bases, 54);
    assert_eq!(
        graph.chains.iter().map(|c| c.sequence.len()).sum::<usize>(),
        54
    );
    assert_eq!(e.duplex.as_ref().unwrap().top.as_str(), ECO_MAIN);
    assert!(
        graph
            .components
            .iter()
            .any(|c| c.covalently_closed && c.paired_nucleotides == 50)
    );
    assert!(
        graph
            .components
            .iter()
            .any(|c| !c.covalently_closed && c.unpaired_nucleotides == 4)
    );
    let free = graph
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert_eq!(
        free.sources,
        vec![p::SourcePlacement {
            source_strand: 3,
            start: 7,
            end: 11,
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
    let cuts = chemistry::digest(
        e.duplex.as_ref().unwrap(),
        &[collagene_core::analysis::Enzyme {
            name: "Independent EcoRI cut geometry".into(),
            recognition: "GAATTC".into(),
            top_cut: 1,
            bottom_cut: 5,
        }],
    )
    .unwrap();
    assert_eq!(cuts.len(), 1);
    assert_eq!(cuts[0].top_start, 19);
    assert_eq!(cuts[0].bottom_start, 23);
    assert_eq!(report.products[1].sequence.as_str(), ECO_TARGET);
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
    let html = v::export_protocol_html(&p, &target);
    assert!(html.contains("hydrolytic"));
    assert!(html.contains("GGTG") || html.contains("excised"));
    assert!(html.contains("not physical isolation") || html.contains("not isolation"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&p, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, report.status);
    assert_eq!(pdf.plan_revision, report.plan_revision);
    assert_eq!(pdf.target_revision, target.revision());
}
#[test]
fn original52call_oracle_retains_free_fragment_and_wrong_target_never_selects_material() {
    let (inputs, r) = fixture(false);
    let r = select(&inputs, r);
    let target = Molecule::dna(
        "Original independent24bp contour",
        ORIGINAL_MAIN,
        Topology::Circular,
    )
    .unwrap();
    let mut p = plan(&inputs, r, &target);
    let report = v::verify(&p, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(evidence(&report).projection.processing.source_bases, 52);
    assert_eq!(report.products[0].sequence.as_str(), ORIGINAL_MAIN);
    assert!(
        evidence(&report)
            .projection
            .processing
            .product
            .components
            .iter()
            .any(|c| c.paired_nucleotides == 48)
    );
    let wrong = Molecule::dna(
        "Independent wrong specification",
        &format!("{ORIGINAL_MAIN}A"),
        Topology::Circular,
    )
    .unwrap();
    p.target_revision = wrong.revision();
    let report = v::verify(&p, &wrong);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.products[0].sequence.as_str(), ORIGINAL_MAIN);
    assert!(evidence(&report).tracked.is_some());
}
#[test]
fn no_cut_retains_anchored_flap_and_full_generic_typed_evidence_without_product_substitution() {
    let (inputs, mut r) = fixture(false);
    r.processing.actions.remove(0);
    let mut r = select(&inputs, r);
    r.options.require_all_requested_transfers_applied = false;
    let target = Molecule::dna(
        "Independent intended contour",
        ORIGINAL_MAIN,
        Topology::Circular,
    )
    .unwrap();
    let report = v::verify(&plan(&inputs, r, &target), &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.products.is_empty());
    let e = evidence(&report);
    assert!(e.tracked.is_none());
    assert!(e.product.is_none());
    let graph = &e.projection.processing.product;
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
    assert!(report.processed_step_reports.is_empty());
}
#[test]
fn ordered_and_blocked_actions_preserve_actual_branch_and_never_infer_phosphate() {
    let (inputs, mut r) = fixture(false);
    r.processing.actions.swap(0, 2);
    let mut r = select(&inputs, r);
    r.options.require_all_requested_transfers_applied = false;
    let target = Molecule::dna(
        "Independent retainednick contour",
        ORIGINAL_MAIN,
        Topology::Circular,
    )
    .unwrap();
    let report = v::verify(&plan(&inputs, r, &target), &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(evidence(&report).duplex.as_ref().unwrap().nicks.len(), 1);
    assert_eq!(
        evidence(&report).projection.processing.steps[0].outcome(),
        p::Outcome::BlockedUnderDeclaredModel
    );
    let (inputs, mut r) = fixture(false);
    let topo::MaterialRecipe::Prepared { strands, .. } = &mut r.preparation.recipe.materials[0]
    else {
        panic!()
    };
    strands[1].five_prime = p::FivePrime::Hydroxyl;
    let mut r = select(&inputs, r);
    r.options.require_all_requested_transfers_applied = false;
    r.options.require_known_terminal_chemistry = false;
    let report = v::verify(&plan(&inputs, r, &target), &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(evidence(&report).duplex.is_none());
    assert_eq!(
        evidence(&report).projection.processing.steps[2].outcome(),
        p::Outcome::BlockedUnderDeclaredModel
    );
}
#[test]
fn op20_candidate14_and_worker66_gates_preserve_original_source_context() {
    let (inputs, r) = fixture(false);
    let r = select(&inputs, r);
    let target = Molecule::dna("Independent", ORIGINAL_MAIN, Topology::Circular).unwrap();
    let p = plan(&inputs, r, &target);
    assert_eq!(postcard::to_allocvec(&p.steps[0]).unwrap()[0], 20);
    let candidate = transfer::Candidate::new(p.clone());
    let bytes = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        u16::from_le_bytes([bytes[4], bytes[5]]),
        protocol::CANDIDATE_SCHEMA_VERSION
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&13u16.to_le_bytes());
    assert!(transfer::load_candidate(&old).is_err());
    let encoded = protocol::encode_request(&protocol::Request {
        id: 2001,
        revision: target.revision(),
        job: protocol::AnalysisJob::Verify {
            plan: p.clone(),
            target: target.clone(),
        },
    })
    .unwrap();
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&65u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::process_request(&encoded).unwrap();
    let mut old = response.clone();
    old[4..6].copy_from_slice(&65u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let protocol::AnalysisResult::Verification(current) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("wrong result")
    };
    assert_eq!(current.status, Status::ConsistentUnderAssumptions);
    assert_eq!(evidence(&current).projection.processing.source_bases, 52);
    let mut stale = p.clone();
    stale.inputs[0]
        .molecule
        .metadata
        .insert("different original capture".into(), "true".into());
    let report = v::verify(&stale, &target);
    assert_eq!(report.status, Status::Invalid);
    assert!(report.products.is_empty());
    let Operation::TopoActions { sources, .. } = &p.steps[0] else {
        panic!()
    };
    assert_eq!(sources[0].role, "original-0");
    assert!(matches!(&sources[0].source,Source::Input(id) if id=="original-0"));
}
#[test]
fn nick_only_operation_cannot_cut_directional_flap_in_current_encoding() {
    let (inputs, r) = fixture(false);
    let sealing = p::nick_sealing::Request {
        ligase_assumption: "Only the old original exposed nick is selected; no flap action".into(),
        nicks: vec![p::nick_sealing::NickSelection {
            three: call(2, 9),
            five: call(0, 0),
        }],
    };
    let treatment = collagene_core::topo_nick_sealing::Request {
        preparation: r.preparation,
        sealing,
    };
    let fresh = collagene_core::topo_nick_sealing::inspect(&inputs, &treatment).unwrap();
    let old = collagene_core::topo_processed_product::Request {
        treatment,
        selection: p::processed_projection::Selection {
            sealing_context_revision: fresh.sealing.context_revision,
            component: 0,
            reference_seed: call(0, 0),
        },
        options: projection::Options::default(),
    };
    let target = Molecule::dna(
        "Independently specified24bp reference",
        ORIGINAL_MAIN,
        Topology::Circular,
    )
    .unwrap();
    let old_plan = DesignPlan {
        name: "Unchanged nick-only operation19".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .map(|i| PlanInput {
                id: i.role.clone(),
                revision: i.molecule.revision(),
                molecule: i.molecule.clone(),
            })
            .collect(),
        steps: vec![Operation::ProcessedTopo {
            sources: inputs
                .iter()
                .map(|i| TopoSource {
                    role: i.role.clone(),
                    source: Source::Input(i.role.clone()),
                })
                .collect(),
            request: old,
        }],
    };
    assert_eq!(postcard::to_allocvec(&old_plan.steps[0]).unwrap()[0], 19);
    let report = v::verify(&old_plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.processed_step_reports.len(), 1);
    assert!(report.physical_processing_reports.is_empty());
    assert!(report.products.is_empty());
    let candidate = transfer::save_candidate(&transfer::Candidate::new(old_plan.clone())).unwrap();
    assert_eq!(
        &candidate[4..6],
        &protocol::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    transfer::load_candidate(&candidate).unwrap();
    let encoded = protocol::encode_request(&protocol::Request {
        id: 1919,
        revision: target.revision(),
        job: protocol::AnalysisJob::Verify {
            plan: old_plan,
            target,
        },
    })
    .unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::process_request(&encoded).unwrap();
    let protocol::AnalysisResult::Verification(report) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("current nick-only processing result")
    };
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.processed_step_reports.len(), 1);
}
#[test]
fn native_plan12_and_generic_run_sidecar_keep_main_and_excised_source_events_historical() {
    use collagene_core::{
        physical_run_evidence, processed_run_evidence,
        run_lineage::{self, EventBinding, InputEventBinding},
        workspace::{self, Workspace},
    };
    fn event(w: &Workspace, index: usize) -> EventBinding {
        let d = &w.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    let (inputs, r) = fixture(true);
    let r = select(&inputs, r);
    let target = Molecule::dna(
        "Independent finalEcoRI product",
        ECO_TARGET,
        Topology::Linear,
    )
    .unwrap();
    let mut p = plan(&inputs, r, &target);
    p.steps.push(digest());
    let mut documents = inputs
        .iter()
        .map(|i| i.molecule.clone())
        .collect::<Vec<_>>();
    documents.push(target.clone());
    let frame =
        Workspace::from_documents("Original directional source events", documents, 4).unwrap();
    let bindings = (0..4)
        .map(|i| InputEventBinding {
            input_id: inputs[i].role.clone(),
            event: event(&frame, i),
        })
        .collect();
    let commit = run_lineage::record_run(&frame, &p, event(&frame, 4), bindings, false).unwrap();
    let archive = workspace::save_workspace(&commit.workspace).unwrap();
    let retained = workspace::load_workspace(&archive).unwrap();
    let record = run_lineage::load_run(&retained, &commit.receipt.run_id).unwrap();
    assert_eq!(record.plan_schema, workspace::PLAN_SCHEMA_VERSION);
    let saved = physical_run_evidence::load(&retained, &commit.receipt.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(saved.steps.len(), 1);
    assert_eq!(saved.target, event(&frame, 4));
    for (i, binding) in saved.inputs.iter().enumerate() {
        assert_eq!(binding.event, event(&frame, i));
    }
    let PhysicalProcessingEvidence::TopoActions(e) = &saved.steps[0].evidence else {
        panic!("expected retained TopoActions evidence")
    };
    assert_eq!(e.projection.processing.source_bases, 54);
    assert!(
        e.projection
            .processing
            .product
            .chains
            .iter()
            .any(|c| c.sequence.as_str() == "GGTG")
    );
    let replay = run_lineage::replay_run(&retained, &commit.receipt.run_id).unwrap();
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        evidence(&replay.current_report)
            .projection
            .processing
            .source_bases,
        54
    );
    let sidecar = retained
        .attachments
        .iter()
        .find(|(_, a)| a.kind == physical_run_evidence::KIND)
        .unwrap()
        .0
        .clone();
    assert!(
        processed_run_evidence::load(&retained, &commit.receipt.run_id)
            .unwrap()
            .is_none(),
        "new generic evidence never fabricates the old nick-only record"
    );
    let download = protocol::encode_request(&protocol::Request {
        id: 149,
        revision: run_lineage::workspace_fingerprint(&retained).unwrap(),
        job: protocol::AnalysisJob::LoadProcessedRunEvidence {
            workspace: retained.clone(),
            run_id: commit.receipt.run_id.clone(),
        },
    })
    .unwrap();
    let response = protocol::process_request(&download).unwrap();
    let protocol::AnalysisResult::HistoricalProcessedRunEvidence(download) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("historical bundle")
    };
    assert_eq!(download.step_count, 1);
    let bundle: serde_json::Value = serde_json::from_slice(&download.bytes.unwrap()).unwrap();
    assert_eq!(bundle["format"], "physical-processing-run-bundle-v1");
    assert_eq!(bundle["attachments"].as_array().unwrap().len(), 1);
    assert_eq!(
        bundle["attachments"][0]["kind"],
        physical_run_evidence::KIND
    );
    assert_eq!(bundle["attachments"][0]["blake3"], sidecar);
    let mut missing = retained.clone();
    missing.attachments.remove(&sidecar);
    assert!(missing.validate().is_err());
    let mut historic = saved.clone();
    let PhysicalProcessingEvidence::TopoActions(e) = &mut historic.steps[0].evidence else {
        panic!("expected retained TopoActions evidence")
    };
    e.engine_version = "historical-actions-inner-engine-v0".into();
    let mut historic_frame = retained.clone();
    let mut attachment = historic_frame.attachments.remove(&sidecar).unwrap();
    attachment.bytes = serde_json::to_vec(&historic).unwrap();
    historic_frame.attachments.insert(
        blake3::hash(&attachment.bytes).to_hex().to_string(),
        attachment,
    );
    historic_frame.validate().unwrap();
    let historical = physical_run_evidence::load(&historic_frame, &commit.receipt.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(historical.engine_version, saved.engine_version);
    let fresh = run_lineage::replay_run(&historic_frame, &commit.receipt.run_id).unwrap();
    assert_eq!(
        evidence(&fresh.current_report).engine_version,
        actions::ENGINE_VERSION
    );
    let portable = fresh.current_report.products[0].clone();
    assert!(chemistry::has_unmapped_physical_evidence(&portable));
    let imported = DesignPlan {
        name: "Historical reference is not isolated material".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "portable".into(),
            revision: portable.revision(),
            molecule: portable,
        }],
        steps: vec![Operation::Digest {
            template: Source::Input("portable".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        }],
    };
    assert_eq!(v::verify(&imported, &target).status, Status::Incomplete);
    // Relabelling the native manifest must not hide a new operation under plan11.
    use std::io::{Cursor, Read, Write};
    let mut input = zip::ZipArchive::new(Cursor::new(&archive)).unwrap();
    let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for i in 0..input.len() {
        let mut entry = input.by_index(i).unwrap();
        let name = entry.name().to_string();
        let mut bytes = vec![];
        entry.read_to_end(&mut bytes).unwrap();
        if name == "manifest.json" {
            let mut manifest: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(manifest["plan_schema"], workspace::PLAN_SCHEMA_VERSION);
            manifest["plan_schema"] = 11.into();
            bytes = serde_json::to_vec(&manifest).unwrap();
        }
        output
            .start_file(name, zip::write::SimpleFileOptions::default())
            .unwrap();
        output.write_all(&bytes).unwrap();
    }
    assert!(workspace::load_workspace(&output.finish().unwrap().into_inner()).is_err());
}
