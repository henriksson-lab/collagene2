//! Fresh typed processed-material science, downstream chemistry and independent intent.
use collagene_core::{
    Molecule, Topology, chemistry, modified_duplex, topo, topo_nick_sealing as treatment,
    topo_processed_product as adapter,
    verification::{self as v, DesignPlan, Operation, PlanInput, Source, Status, TopoSource},
};
use topo::physical::{self as p, nick_sealing as n, processed_projection as pp, projection};
fn fixture(directional: bool) -> (Vec<topo::Input>, treatment::Request) {
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
        transfers:vec![p::Transfer{donor:0,acceptor:2,geometry:if directional{p::Geometry::DirectionalCacc}else{p::Geometry::ThreePrimeTA}},p::Transfer{donor:1,acceptor:3,geometry:if directional{p::Geometry::Flush}else{p::Geometry::ThreePrimeTA}}],
        engine_limits:p::Limits::default(),projection_options:projection::Options::default()
    };
    let sealing = n::Request {
        ligase_assumption:
            "Explicit ligase hypothesis: only templated 3′OH/5′P nick sealing, no repair enzymes"
                .into(),
        nicks: vec![
            n::NickSelection {
                three: p::Call {
                    strand: 2,
                    index: 6,
                },
                five: p::Call {
                    strand: 0,
                    index: 0,
                },
            },
            n::NickSelection {
                three: p::Call {
                    strand: 3,
                    index: 6,
                },
                five: p::Call {
                    strand: 1,
                    index: 0,
                },
            },
        ],
    };
    (
        inputs,
        treatment::Request {
            preparation,
            sealing,
        },
    )
}
fn selected(inputs: &[topo::Input], treatment: treatment::Request) -> adapter::Request {
    let native = treatment::inspect(inputs, &treatment).unwrap();
    adapter::Request {
        treatment,
        selection: pp::Selection {
            sealing_context_revision: native.sealing.context_revision,
            component: 0,
            reference_seed: p::Call {
                strand: 0,
                index: 0,
            },
        },
        options: projection::Options {
            require_all_requested_transfers_applied: true,
            require_known_terminal_chemistry: true,
            ..Default::default()
        },
    }
}
const CIRCLE: &str = "AGGGGCATCCCTTCGAATGA";
const ECORI_CIRCLE: &str = "AGGGGCATCCCTTCGAATTCA";
const LINEAR_ECORI: &str = "AATTCAAGGGGCATCCCTTCG";
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
        steps: vec![Operation::ProcessedTopo {
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
fn ecorifixture() -> (Vec<topo::Input>, treatment::Request) {
    let (mut inputs, mut treatment) = fixture(false);
    inputs[2].molecule =
        Molecule::dna("Original EcoRI insert top", "CGAATTCA", Topology::Linear).unwrap();
    inputs[3].molecule =
        Molecule::dna("Original EcoRI insert bottom", "GAATTCGA", Topology::Linear).unwrap();
    let topo::MaterialRecipe::Prepared {
        strands, pairing, ..
    } = &mut treatment.preparation.recipe.materials[0]
    else {
        panic!()
    };
    strands[2].segments[0].length = 8;
    strands[3].segments[0].length = 8;
    pairing[1].length = 7;
    for nick in &mut treatment.sealing.nicks {
        nick.three.index = 7;
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
    assert!(report.engine_version.ends_with("/science-16"));
    assert_eq!(report.products.len(), 2);
    assert_eq!(report.products[1].sequence.as_str(), LINEAR_ECORI);
    assert_eq!(report.processed_step_reports.len(), 1);
    assert_eq!(report.processed_step_reports[0].step_index, 0);
    assert!(
        report.step_reports.is_empty(),
        "no fabricated legacy topo report"
    );
    let evidence = &report.processed_step_reports[0].evidence;
    assert_eq!(evidence.projection.sealing.parent.product.nicks.len(), 2);
    assert_eq!(evidence.projection.sealing.product.nicks.len(), 0);
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
    assert!(html.contains("Processed TOPO physical step evidence"));
    assert!(html.contains("Explicit ligase hypothesis"));
    assert!(html.contains("original graph 2 nicks; processed graph 0 nicks"));
    assert!(html.contains("not physical isolation"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, report.status);
    assert_eq!(pdf.plan_revision, report.plan_revision);
    assert_eq!(pdf.target_revision, target.revision());
    let encoded = postcard::to_allocvec(&plan.steps[0]).unwrap();
    assert_eq!(encoded[0], 19);
    let candidate = collagene_core::transfer::Candidate::new(plan.clone());
    let bytes = collagene_core::transfer::save_candidate(&candidate).unwrap();
    assert_eq!(u16::from_le_bytes([bytes[4], bytes[5]]), 13);
    let mut old = bytes;
    old[4..6].copy_from_slice(&12u16.to_le_bytes());
    assert!(collagene_core::transfer::load_candidate(&old).is_err());
    let request = collagene_core::protocol::encode_request(&collagene_core::protocol::Request {
        id: 65,
        revision: target.revision(),
        job: collagene_core::protocol::AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
    })
    .unwrap();
    let mut old = request.clone();
    old[4..6].copy_from_slice(&64u16.to_le_bytes());
    assert!(collagene_core::protocol::process_request(&old).is_err());
    let response = collagene_core::protocol::process_request(&request).unwrap();
    let mut old = response.clone();
    old[4..6].copy_from_slice(&64u16.to_le_bytes());
    assert!(collagene_core::protocol::decode_response(&old).is_err());
    let collagene_core::protocol::AnalysisResult::Verification(rebuilt) =
        collagene_core::protocol::decode_response(&response)
            .unwrap()
            .result
    else {
        panic!("fresh typed verifier")
    };
    assert_eq!(rebuilt.processed_step_reports.len(), 1);
}
#[test]
fn unsealed_nick_is_retained_and_downstream_digest_cannot_erase_it() {
    let (inputs, mut treatment) = fixture(false);
    treatment.sealing.nicks.truncate(1);
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
    let evidence = &report.processed_step_reports[0].evidence;
    assert_eq!(evidence.duplex.as_ref().unwrap().nicks.len(), 1);
    assert_eq!(evidence.projection.sealing.product.nicks.len(), 1);
    plan.steps.push(digest());
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.products.len(), 1);
    assert_eq!(
        report.processed_step_reports[0]
            .evidence
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
    let fresh = topo::simulate(&inputs, &treatment.preparation).unwrap();
    treatment.sealing.nicks = fresh
        .physical
        .product
        .nicks
        .iter()
        .map(|k| n::NickSelection {
            three: k.three,
            five: k.five,
        })
        .collect();
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
    assert!(report.step_reports.is_empty());
    let evidence = &report.processed_step_reports[0].evidence;
    assert!(evidence.tracked.is_none());
    assert!(evidence.product.is_none());
    let graph = &evidence.projection.sealing.product;
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
    assert!(html.contains("Processed TOPO physical step evidence"));
    assert!(html.contains("Faithful duplex bridge available: false"));
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
    assert_eq!(report.processed_step_reports.len(), 1);
    assert!(report.processed_step_reports[0].evidence.tracked.is_some());
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
    let Operation::ProcessedTopo { sources, .. } = &mut nested.steps[0] else {
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
    assert!(report.processed_step_reports.is_empty());
    assert!(
        report
            .findings
            .iter()
            .any(|s| s.contains("Tracked Step ancestry is not mapped"))
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
    let Operation::ProcessedTopo { request, .. } = &mut plan.steps[0] else {
        panic!()
    };
    request.options.limits.max_pairs = 0;
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert!(report.plan_revision.is_empty());
    assert!(report.products.is_empty());
    let Operation::ProcessedTopo { request, .. } = &mut plan.steps[0] else {
        panic!()
    };
    request.options.limits.max_pairs = 50_000;
    request.treatment.preparation.engine_limits.max_steps = 3;
    let report = v::verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert!(report.plan_revision.is_empty());
}
#[test]
fn immutable_run_graph_sidecar_survives_native_history_without_restoring_imported_material() {
    use collagene_core::{
        processed_run_evidence,
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
    let (inputs, treatment) = ecorifixture();
    let request = selected(&inputs, treatment);
    let target = Molecule::dna("Independent EcoRI target", LINEAR_ECORI, Topology::Linear).unwrap();
    let mut p = plan(&inputs, request, &target);
    p.steps.push(digest());
    let mut documents = inputs
        .iter()
        .map(|i| i.molecule.clone())
        .collect::<Vec<_>>();
    documents.push(target.clone());
    let w = Workspace::from_documents("Retained actual processed run", documents, 4).unwrap();
    let sources = (0..4)
        .map(|i| InputEventBinding {
            input_id: inputs[i].role.clone(),
            event: event(&w, i),
        })
        .collect();
    let commit = run_lineage::record_run(&w, &p, event(&w, 4), sources, false).unwrap();
    let bytes = workspace::save_workspace(&commit.workspace).unwrap();
    let retained = workspace::load_workspace(&bytes).unwrap();
    {
        use std::io::{Cursor, Read, Write};
        let mut input = zip::ZipArchive::new(Cursor::new(&bytes)).unwrap();
        let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for index in 0..input.len() {
            let mut entry = input.by_index(index).unwrap();
            let name = entry.name().to_owned();
            let mut data = vec![];
            entry.read_to_end(&mut data).unwrap();
            if name == "manifest.json" {
                let mut manifest: serde_json::Value = serde_json::from_slice(&data).unwrap();
                assert_eq!(manifest["plan_schema"], 11);
                manifest["plan_schema"] = 10.into();
                data = serde_json::to_vec(&manifest).unwrap();
            }
            output
                .start_file(name, zip::write::SimpleFileOptions::default())
                .unwrap();
            output.write_all(&data).unwrap();
        }
        assert!(workspace::load_workspace(&output.finish().unwrap().into_inner()).is_err());
    }

    let stored = processed_run_evidence::load(&retained, &commit.receipt.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(stored.steps.len(), 1);
    assert_eq!(stored.steps[0].step_index, 0);
    assert_eq!(stored.target, event(&w, 4));
    assert_eq!(stored.inputs.len(), 4);
    for (i, binding) in stored.inputs.iter().enumerate() {
        assert_eq!(binding.event, event(&w, i));
    }
    assert_eq!(stored.products.len(), 2);
    assert_eq!(
        stored.steps[0]
            .evidence
            .projection
            .sealing
            .parent
            .product
            .nicks
            .len(),
        2
    );
    assert_eq!(
        stored.steps[0]
            .evidence
            .projection
            .sealing
            .product
            .nicks
            .len(),
        0
    );
    let replayed = run_lineage::replay_run(&retained, &commit.receipt.run_id).unwrap();
    assert_eq!(
        replayed.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(replayed.current_report.processed_step_reports.len(), 1);
    let sidecar_hash = retained
        .attachments
        .iter()
        .find(|(_, a)| a.kind == processed_run_evidence::KIND)
        .unwrap()
        .0
        .clone();
    let download_bytes =
        collagene_core::protocol::encode_request(&collagene_core::protocol::Request {
            id: 149,
            revision: run_lineage::workspace_fingerprint(&retained).unwrap(),
            job: collagene_core::protocol::AnalysisJob::LoadProcessedRunEvidence {
                workspace: retained.clone(),
                run_id: commit.receipt.run_id.clone(),
            },
        })
        .unwrap();
    let mut old = download_bytes.clone();
    old[4..6].copy_from_slice(&64u16.to_le_bytes());
    assert!(collagene_core::protocol::process_request(&old).is_err());
    let response = collagene_core::protocol::process_request(&download_bytes).unwrap();
    let mut old = response.clone();
    old[4..6].copy_from_slice(&64u16.to_le_bytes());
    assert!(collagene_core::protocol::decode_response(&old).is_err());
    let collagene_core::protocol::AnalysisResult::HistoricalProcessedRunEvidence(download) =
        collagene_core::protocol::decode_response(&response)
            .unwrap()
            .result
    else {
        panic!("historical exact-byte download")
    };
    assert_eq!(download.run_id, commit.receipt.run_id);
    assert_eq!(download.step_count, 1);
    assert_eq!(
        download.bytes.unwrap(),
        retained.attachments[&sidecar_hash].bytes
    );
    // Historical inner algorithm labels and opaque contexts remain documentary:
    // future engine changes must not make a retained native archive unopenable.
    fn replace_evidence(
        w: &Workspace,
        hash: &str,
        e: &processed_run_evidence::Evidence,
    ) -> Workspace {
        let mut next = w.clone();
        let mut attachment = next.attachments.remove(hash).unwrap();
        attachment.bytes = serde_json::to_vec(e).unwrap();
        let key = blake3::hash(&attachment.bytes).to_hex().to_string();
        next.attachments.insert(key, attachment);
        next
    }
    let mut historical = stored.clone();
    historical.steps[0].evidence.engine_version = "historic-processed-inner-algorithm-v0".into();
    let historical_workspace = replace_evidence(&retained, &sidecar_hash, &historical);
    historical_workspace.validate().unwrap();
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&historical_workspace).unwrap())
            .unwrap();
    let historical = processed_run_evidence::load(&reopened, &commit.receipt.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(
        historical.engine_version, stored.engine_version,
        "outer recorded scientific engine remains exact"
    );
    assert_eq!(
        historical.steps[0].evidence.engine_version,
        "historic-processed-inner-algorithm-v0"
    );
    let fresh = run_lineage::replay_run(&reopened, &commit.receipt.run_id).unwrap();
    assert_eq!(
        fresh.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        fresh.current_report.processed_step_reports[0]
            .evidence
            .engine_version,
        adapter::ENGINE_VERSION
    );
    assert_ne!(
        fresh.current_report.processed_step_reports[0]
            .evidence
            .engine_version,
        historical.steps[0].evidence.engine_version
    );
    let mut changed_outer = stored.clone();
    changed_outer.engine_version = "forged-outer-engine".into();
    assert!(
        replace_evidence(&retained, &sidecar_hash, &changed_outer)
            .validate()
            .is_err()
    );
    let mut changed_request = stored.clone();
    changed_request.steps[0]
        .evidence
        .request
        .treatment
        .sealing
        .ligase_assumption
        .push_str(" changed declaration");
    assert!(
        replace_evidence(&retained, &sidecar_hash, &changed_request)
            .validate()
            .is_err()
    );
    let mut changed_event = stored.clone();
    changed_event.inputs[0].event = event(&w, 4);
    assert!(
        replace_evidence(&retained, &sidecar_hash, &changed_event)
            .validate()
            .is_err()
    );
    let mut missing = retained.clone();
    missing.attachments.remove(&sidecar_hash);
    assert!(missing.validate().is_err());
    let mut orphan = retained.clone();
    orphan
        .attachments
        .retain(|_, a| a.kind != run_lineage::WORKSPACE_ATTACHMENT_KIND);
    assert!(orphan.validate().is_err());
    let mut tampered = retained.clone();
    tampered.attachments.get_mut(&sidecar_hash).unwrap().bytes[0] ^= 1;
    assert!(tampered.validate().is_err());
    // Historical load provides reference documents and evidence only. An imported
    // run product retains restrictions; only replaying its plan restores private state.
    let product = &replayed.current_report.products[0];
    assert!(chemistry::has_unmapped_physical_evidence(product));
    let imported = DesignPlan {
        name: "Historical document is not material".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "historic".into(),
            revision: product.revision(),
            molecule: product.clone(),
        }],
        steps: vec![Operation::Digest {
            template: Source::Input("historic".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        }],
    };
    let report = v::verify(&imported, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.products.is_empty());
}
