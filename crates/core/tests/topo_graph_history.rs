//! Native saved graphs retain actual transfers and source calls; no report is stock.
#[path = "common/topo_graph.rs"]
mod literal;
use collagene_core::{
    Feature, Location, Molecule, Topology, cloning_outputs,
    history::DocumentHistory,
    plan_oligos, run_lineage, topo as t, topo_graph_history as g, topo_run_preparation as native,
    verification::{self, DesignPlan, Status},
    workspace::{self, Workspace, WorkspaceAttachment},
};
use t::physical as p;
fn reopened(w: &Workspace) -> Workspace {
    workspace::load_workspace(&workspace::save_workspace(w).unwrap()).unwrap()
}
fn request_product(request: &t::Request) -> t::Product {
    let selection = t::Selection {
        batch_revision: t::batch_revision(&[], request).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    };
    t::product(&[], request, &selection).unwrap()
}
fn checked(w: &Workspace, request: &t::Request) -> (DesignPlan, Molecule) {
    let target = w.revisions[&literal::event(w, 4).revision].clone();
    let task = native::Task::Check {
        request: request.clone(),
        selection: t::Selection {
            batch_revision: t::batch_revision(&[], request).unwrap(),
            component: 0,
            reference_seed: p::Call {
                strand: 0,
                index: 0,
            },
        },
        target: target.clone(),
        policy: t::AcceptancePolicy {
            allow_nicks: true,
            ..Default::default()
        },
    };
    (native::compile(&[], &target, None, &task).unwrap(), target)
}
#[test]
fn saved_source_events_survive_both_head_advances_and_reopen_without_ordinary_duplex_authority() {
    let (original, r) = literal::fixture();
    let before_docs = serde_json::to_vec(&original.documents).unwrap();
    let before_history = serde_json::to_vec(&original.history).unwrap();
    let first = g::prepare_capture(&original, &r).unwrap();
    let old = first.receipt.selection.clone();
    assert_eq!(
        serde_json::to_vec(&first.workspace.documents).unwrap(),
        before_docs
    );
    assert_eq!(
        serde_json::to_vec(&first.workspace.history).unwrap(),
        before_history
    );
    let mut w = first.workspace;
    let mut histories = DocumentHistory::restore_all(&w).unwrap();
    for history in histories.iter_mut().take(2) {
        let source = history.revision_molecule(history.head_event()).unwrap();
        let changed = source.replace(0, 1, "T").unwrap();
        history
            .commit(changed, "Explicitly advance source head")
            .unwrap();
    }
    w.documents.clear();
    w.history.clear();
    for history in &histories {
        history.contribute(&mut w).unwrap();
    }
    assert_ne!(serde_json::to_vec(&w.history).unwrap(), before_history);
    let w = reopened(&w);
    let inspected = g::inspect(&w, &old).unwrap();
    assert_eq!(inspected.report.physical.source_bases, 26);
    assert_eq!(inspected.summary.remaining_activations.len(), 2);
    assert!(inspected.summary.released_activations.is_empty());
    assert_eq!(
        inspected.report.preparation.substrates.strands[0]
            .sequence
            .as_str(),
        "AGGGGCATCCCTT"
    );
    assert_eq!(
        inspected.report.preparation.substrates.strands[1]
            .sequence
            .as_str(),
        "AGGGATGCCCCTT"
    );
    assert!(matches!(
        inspected.report.physical.product.chains[0].ends,
        p::ChainEnds::Linear {
            three_prime: p::ThreePrime::Phosphotyrosyl { .. },
            ..
        }
    ));
    let partial = g::prepare_capture(&w, &literal::partial(&w, &old)).unwrap();
    assert_eq!(partial.receipt.summary.source_bases, 40);
    assert_eq!(partial.receipt.summary.chains, 3);
    let request = literal::continuation(&partial.workspace, &partial.receipt.selection);
    let product = request_product(&request);
    assert_eq!(
        product.reference.as_ref().unwrap().sequence.as_str(),
        literal::FINAL
    );
    assert!(
        collagene_core::chemistry::Duplex::from_molecule(product.reference.as_ref().unwrap())
            .is_err()
    );
    assert_eq!(
        product.projection.physical.product.chains[1]
            .sequence
            .as_str(),
        "AGGGATGCCCCTTCATTCGA"
    );
    assert_eq!(product.projection.physical.product.nicks.len(), 2);
}
#[test]
fn orientation_resume_satisfies_exact_prefix_once_and_keeps_wrong_orientation_diagnostic() {
    let initial = literal::captured();
    let partial = g::prepare_capture(
        &initial.workspace,
        &literal::partial(&initial.workspace, &initial.receipt.selection),
    )
    .unwrap();
    let mut request = literal::continuation(&partial.workspace, &partial.receipt.selection);
    request.transfers.clear();
    let choices = t::inspect_orientations(&[], &request, [0, 1], [2, 3], false).unwrap();
    assert_eq!(choices.reference.candidates.len(), 2);
    assert_eq!(choices.reference.satisfied_prefix_steps[0], vec![0]);
    assert_eq!(
        choices.reference.pending_transfers[0],
        vec![p::Transfer {
            donor: 1,
            acceptor: 3,
            geometry: p::Geometry::ThreePrimeTA
        }]
    );
    let first = &choices.reference.candidates[0];
    assert_eq!(first.steps.len(), 2);
    assert_eq!(first.product.released_activations.len(), 2);
    assert_eq!(first.product.chains[0].sequence.as_str(), literal::FINAL);
    assert!(!choices.reference.candidates[1].all_requested_transfers_applied);
    assert!(
        choices.reference.candidates[1]
            .steps
            .iter()
            .any(|s| s.outcome != p::Outcome::Applied)
    );
    assert_eq!(
        choices.candidate_requests[0].transfers.len(),
        1,
        "prefix is owned by saved source, not appended twice"
    );
    let actual = request_product(&choices.candidate_requests[0]);
    assert_eq!(actual.projection.physical.product, first.product);
}
#[test]
fn renamed_declarations_labels_and_role_aliases_cannot_create_an_extra_activated_instance() {
    let (w, original) = literal::fixture();
    let first = g::prepare_capture(&w, &original).unwrap();
    let mut renamed = original.clone();
    renamed.name = "New label, same physical source calls".into();
    for (i, binding) in renamed.bindings.iter_mut().enumerate() {
        binding.role = format!("renamed-source-{i}");
    }
    if let t::MaterialRecipe::Prepared {
        declaration_id,
        assumption,
        strands,
        ..
    } = &mut renamed.request.recipe.materials[0]
    {
        *declaration_id = "cosmetic declaration rename".into();
        *assumption = "Same original hypotheses, reworded".into();
        for (i, strand) in strands.iter_mut().enumerate() {
            strand.id = format!("renamed-strand-{i}");
            strand.segments[0].source_role = renamed.bindings[i].role.clone();
            if let p::ThreePrime::Phosphotyrosyl { activation_id, .. } = &mut strand.three_prime {
                *activation_id = format!("cosmetic-adduct-label-{i}");
            }
        }
    }
    let second = g::prepare_capture(&first.workspace, &renamed).unwrap();
    assert_ne!(first.receipt.record_hash, second.receipt.record_hash);
    let a = g::documentary_capsule(&second.workspace, &first.receipt.selection).unwrap();
    let b = g::documentary_capsule(&second.workspace, &second.receipt.selection).unwrap();
    let doubled = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![
                t::MaterialRecipe::SavedGraph {
                    material_id: "first".into(),
                    origin: a.clone(),
                },
                t::MaterialRecipe::SavedGraph {
                    material_id: "renamed".into(),
                    origin: b,
                },
            ],
        },
        transfers: vec![],
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    assert!(t::simulate(&[], &doubled).is_err());
    // An unbound declaration label cannot manufacture a second saved instance by
    // recapturing the same original source events beside the already saved graph.
    let mut mixed = renamed;
    mixed.request.recipe.materials.insert(
        0,
        t::MaterialRecipe::SavedGraph {
            material_id: "old native source".into(),
            origin: a,
        },
    );
    assert!(g::prepare_capture(&second.workspace, &mixed).is_err());
}
#[test]
fn generated_role_allocator_never_collides_with_a_current_external_role() {
    let initial = literal::captured();
    let mut partial = literal::partial(&initial.workspace, &initial.receipt.selection);
    partial.bindings[0].role = "graph1-source0".into();
    if let t::MaterialRecipe::Prepared { strands, .. } = &mut partial.request.recipe.materials[1] {
        strands[0].segments[0].source_role = "graph1-source0".into();
    }
    let saved = g::prepare_capture(&initial.workspace, &partial).unwrap();
    let inspected = g::inspect(&saved.workspace, &saved.receipt.selection).unwrap();
    let roles = inspected
        .report
        .sources
        .iter()
        .map(|s| s.role.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(roles.len(), 4);
    assert!(roles.contains("graph1-source0"));
    assert!(roles.contains("graph1-source0-1"));
    assert_eq!(inspected.report.physical.source_bases, 40);
}
#[test]
fn nested_capture_final_run_and_protocols_keep_all_calls_without_historical_oligo_orders() {
    let initial = literal::captured();
    let partial = g::prepare_capture(
        &initial.workspace,
        &literal::partial(&initial.workspace, &initial.receipt.selection),
    )
    .unwrap();
    let full_request = g::CaptureRequest {
        name: "Both transfers retained".into(),
        bindings: vec![],
        request: literal::continuation(&partial.workspace, &partial.receipt.selection),
    };
    let full = g::prepare_capture(&partial.workspace, &full_request).unwrap();
    let w = reopened(&full.workspace);
    let resume = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::SavedGraph {
                material_id: "completed native graph".into(),
                origin: g::documentary_capsule(&w, &full.receipt.selection).unwrap(),
            }],
        },
        transfers: vec![],
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    let product = request_product(&resume);
    assert_eq!(product.projection.physical.source_bases, 40);
    assert_eq!(product.projection.physical.product.chains.len(), 2);
    assert_eq!(
        product
            .projection
            .physical
            .product
            .released_activations
            .len(),
        2
    );
    let (plan, target) = checked(&w, &resume);
    let verified = verification::verify(&plan, &target);
    assert_eq!(verified.status, Status::ConsistentUnderAssumptions);
    let orders = plan_oligos::build(&plan, &target, &Default::default()).unwrap();
    assert!(orders.oligos.records.is_empty());
    let html = verification::export_protocol_html(&plan, &target);
    assert!(html.contains("retained source transfers"));
    assert!(html.contains(literal::FINAL));
    assert!(!html.contains("record_bytes"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF"));
    let before = serde_json::to_vec(&w.documents).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        literal::event(&w, 4),
        vec![],
        &cloning_outputs::Options {
            name: "Explicit checked graph reference".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    assert_eq!(
        serde_json::to_vec(&commit.workspace.documents[..w.documents.len()]).unwrap(),
        before
    );
    let restart = reopened(&commit.workspace);
    let replay = run_lineage::replay_run(&restart, &commit.receipt.run.run_id).unwrap();
    assert!(replay.engine_matches);
    assert!(replay.products_match);
    // A graph plus its descendant cannot silently duplicate the inherited vector.
    let mut duplicated = resume;
    duplicated
        .recipe
        .materials
        .push(t::MaterialRecipe::SavedGraph {
            material_id: "same ancestor again".into(),
            origin: g::documentary_capsule(&w, &initial.receipt.selection).unwrap(),
        });
    assert!(t::simulate(&[], &duplicated).is_err());
}
#[test]
fn large_original_documentary_bytes_never_become_protocol_capsule_json() {
    let (original, r) = literal::fixture();
    let mut molecules = original
        .documents
        .iter()
        .map(|d| original.revisions[&d.head].clone())
        .collect::<Vec<_>>();
    molecules[0].original_genbank = Some("PRIVATE_SOURCE_DOCUMENTARY_PAYLOAD ".repeat(20_000));
    molecules[0].features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("retained vector segment".into()))],
        needs_review: false,
    });
    let w =
        Workspace::from_documents("Large raw source, small actual graph", molecules, 0).unwrap();
    let mut r = r;
    for (i, b) in r.bindings.iter_mut().enumerate() {
        b.event = literal::event(&w, i)
    }
    let first = g::prepare_capture(&w, &r).unwrap();
    let partial = g::prepare_capture(
        &first.workspace,
        &literal::partial(&first.workspace, &first.receipt.selection),
    )
    .unwrap();
    let request = literal::continuation(&partial.workspace, &partial.receipt.selection);
    let (plan, target) = checked(&partial.workspace, &request);
    let report = verification::verify(&plan, &target);
    assert_eq!(report.status, Status::ConsistentUnderAssumptions);
    let product = report.step_reports[0].topo.as_ref().unwrap();
    assert!(
        product
            .annotations
            .as_ref()
            .unwrap()
            .features
            .iter()
            .any(|f| f
                .qualifiers
                .iter()
                .any(|(_, v)| v.as_deref() == Some("retained vector segment")))
    );
    let html = verification::export_protocol_html(&plan, &target);
    assert!(!html.contains("PRIVATE_SOURCE_DOCUMENTARY_PAYLOAD"));
    assert!(!html.contains("record_bytes"));
    assert!(html.len() < 100_000);
    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(pdf.bytes.len() < 1_000_000);
}
#[test]
fn historical_documentary_record_remains_readable_but_fresh_replay_and_wrong_source_refuse() {
    let first = literal::captured();
    let mut record = g::load(&first.workspace, &first.receipt.record_hash).unwrap();
    record.scientific_engine = "historical-scientific-model".into();
    record.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            &record.engine_version,
            &record.scientific_engine,
            &record.source_workspace_fingerprint,
            &record.request,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let bytes = g::save_record(&record).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut w = first.workspace;
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: g::WORKSPACE_ATTACHMENT_KIND.into(),
            name: "Historical documentary graph".into(),
            media_type: g::MEDIA_TYPE.into(),
            bytes,
        },
    );
    let w = reopened(&w);
    let selection = g::Selection {
        record_hash: hash,
        state_revision: record.state_revision,
    };
    let capsule = g::documentary_capsule(&w, &selection).unwrap();
    assert!(g::replay(&capsule).is_err());
    assert!(g::load(&w, &selection.record_hash).is_ok());
    let mut wrong = capsule;
    wrong.source_snapshots[0].molecule.name = "Same DNA different full snapshot".into();
    assert!(wrong.validate_shallow().is_err());
    let current = g::documentary_capsule(&w, &first.receipt.selection).unwrap();
    let mut forged = current.clone();
    forged.selection.state_revision = "0".repeat(64);
    assert!(forged.validate_shallow().is_err());
    let mut malformed = g::load(&w, &first.receipt.record_hash).unwrap();
    malformed.context_revision = "0".repeat(64);
    assert!(g::save_record(&malformed).is_err());
}
#[test]
fn saved_vector_new_pcr_orders_only_current_pair_then_resume_orders_none() {
    use collagene_core::{
        analysis::Primer,
        tailing::TailingBase,
        topo_preparation::{
            self, ExactPcrDeclaration, PreparationMode, PreparationRequest, PreparationSource,
        },
    };
    let first = literal::captured();
    let mut w = first.workspace;
    let template = Molecule::dna(
        "Independent new insert PCR template",
        "CGAATG",
        Topology::Linear,
    )
    .unwrap();
    DocumentHistory::new("actual-pcr-template", template.clone())
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let pin = literal::event(&w, w.documents.len() - 1);
    let declaration = PreparationRequest {
        mode: PreparationMode::Ta,
        source: PreparationSource::SavedVector {
            material_id: "saved activated vector".into(),
            origin: g::documentary_capsule(&w, &first.receipt.selection).unwrap(),
            vector_strands: [0, 1],
            insert: ExactPcrDeclaration {
                material_id: "current insert PCR".into(),
                template_role: "template".into(),
                forward: Primer::new("current-forward", "CGA", "").unwrap(),
                reverse: Primer::new("current-reverse", "CAT", "").unwrap(),
                tail: Some(TailingBase::A),
            },
        },
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    let inputs = vec![t::Input {
        role: "template".into(),
        molecule: template.clone(),
    }];
    let mut construction = topo_preparation::compile_preparation(&inputs, &declaration)
        .unwrap()
        .request;
    construction.transfers = vec![
        p::Transfer {
            donor: 0,
            acceptor: 2,
            geometry: p::Geometry::ThreePrimeTA,
        },
        p::Transfer {
            donor: 1,
            acceptor: 3,
            geometry: p::Geometry::ThreePrimeTA,
        },
    ];
    let target = w.revisions[&literal::event(&w, 4).revision].clone();
    let task = native::Task::Check {
        request: construction.clone(),
        selection: t::Selection {
            batch_revision: t::batch_revision(&inputs, &construction).unwrap(),
            component: 0,
            reference_seed: p::Call {
                strand: 0,
                index: 0,
            },
        },
        target: target.clone(),
        policy: t::AcceptancePolicy {
            allow_nicks: true,
            ..Default::default()
        },
    };
    let plan = native::compile(&inputs, &target, None, &task).unwrap();
    assert_eq!(plan.inputs[0].revision, template.revision());
    assert_eq!(
        verification::verify(&plan, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let orders = plan_oligos::build(&plan, &target, &Default::default()).unwrap();
    assert_eq!(orders.oligos.records.len(), 2);
    assert!(
        orders
            .oligos
            .records
            .iter()
            .all(|r| matches!(r.raw_ordered_sequence.as_deref(), Some("CGA" | "CAT")))
    );
    construction.transfers.truncate(1);
    let partial = g::prepare_capture(
        &w,
        &g::CaptureRequest {
            name: "Current PCR and one transfer".into(),
            bindings: vec![g::SourceBinding {
                role: "template".into(),
                event: pin,
                origin: None,
            }],
            request: construction,
        },
    )
    .unwrap();
    assert_eq!(partial.receipt.summary.source_bases, 40);
    let request = literal::continuation(&partial.workspace, &partial.receipt.selection);
    let (resume, target) = checked(&partial.workspace, &request);
    assert!(resume.inputs.is_empty());
    assert_eq!(
        verification::verify(&resume, &target).status,
        Status::ConsistentUnderAssumptions
    );
    assert!(
        plan_oligos::build(&resume, &target, &Default::default())
            .unwrap()
            .oligos
            .records
            .is_empty()
    );
}
