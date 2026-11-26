//! Literal PCR-copy composition; upstream chemistry remains distinct from new DNA.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    methylation::MethylationState,
    modified_duplex::GenerationAssumptions,
    protocol,
    run_lineage::{self, EventBinding, InputEventBinding},
    topo, topo_actions as actions, transfer,
    verification::{self as v, DesignPlan, Operation, PlanInput, Source, Status, TopoSource},
    workspace::{self, Workspace},
};
use topo::physical::{self as p, processing as a, processing_projection as pp, projection};
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
    let preparation = topo::Request {
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
        transfers: vec![
            p::Transfer {
                donor: 0,
                acceptor: 2,
                geometry: p::Geometry::DirectionalCacc,
            },
            p::Transfer {
                donor: 1,
                acceptor: 3,
                geometry: p::Geometry::Flush,
            },
        ],
        engine_limits: p::Limits::default(),
        projection_options: projection::Options::default(),
    };
    let processing = a::Request {
        actions: vec![
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
        ],
        limits: a::ProcessingLimits::default(),
    };
    // Select the actual source component before any target is supplied. Metadata
    // on the later template changes its identity, not these physical source calls.
    let inspected = actions::inspect_actions(
        &inputs,
        &actions::ActionsRequest {
            preparation: preparation.clone(),
            processing: processing.clone(),
        },
    )
    .unwrap();
    let graph = &inspected.processing.product;
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
        processing,
        selection: pp::Selection {
            processing_context_revision: inspected.processing.context_revision,
            component,
            reference_seed: call(0, 0),
        },
        options: projection::Options {
            require_all_requested_transfers_applied: true,
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
use collagene_core::{cloning_outputs, run_material, topo_run_preparation as adapter};
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved(
    staggered: bool,
) -> (
    Workspace,
    Vec<topo::Input>,
    Molecule,
    run_material::Capsule,
    actions::Request,
) {
    let (old, target, mut inputs) = fixture();
    let Operation::TopoActions {
        request: mut req, ..
    } = old.steps[1].clone()
    else {
        panic!("actions")
    };
    let native = inputs[2].molecule.clone();
    let mut originplan = DesignPlan {
        name: "Original new PCR material".into(),
        target_revision: native.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: native.revision(),
            molecule: native.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("Saved F", "CACCA", "").unwrap(),
            reverse: Primer::new("Saved R", "CGTCA", "").unwrap(),
        }],
    };
    let expected = if staggered {
        originplan.steps.push(Operation::Tail3Prime {
            template: Source::Step(0),
            base: collagene_core::tailing::TailingBase::A,
        });
        Molecule::dna("Independent staggered top", "CACCATGACGA", Topology::Linear).unwrap()
    } else {
        native.clone()
    };
    originplan.target_revision = expected.revision();
    let w = Workspace::from_documents(
        "Native template, independent expected and unchanged vectors",
        vec![
            native,
            expected,
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target.clone(),
        ],
        4,
    )
    .unwrap();
    let committed = cloning_outputs::finalize(
        &w,
        &originplan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved insert template".into(),
            show_prepared: false,
            save_primers: true,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: committed.receipt.run.run_id.clone(),
        step_index: if staggered { 1 } else { 0 },
        event: committed.receipt.construct.clone(),
    };
    let capsule = run_material::capsule(&committed.workspace, &selection).unwrap();
    inputs[2].molecule = committed.workspace.revisions[&selection.event.revision].clone();
    if !staggered {
        let task = adapter::Task::InspectActions {
            request: actions::ActionsRequest {
                preparation: req.preparation.clone(),
                processing: req.processing.clone(),
            },
        };
        let inspected = adapter::inspect(&inputs, "template", &capsule, &task).unwrap();
        let adapter::Reply::ActionsInspected(report) = inspected.reply else {
            panic!("inspection")
        };
        req.selection.processing_context_revision = report.processing.context_revision;
    }
    (committed.workspace, inputs, target, capsule, req)
}
#[test]
fn saved_template_new_pcr_processed_product_finalization_archive_and_second_reuse() {
    let (w, inputs, target, origin, request) = saved(false);
    assert!(
        topo::source_preflight(&inputs, &request.preparation).is_err(),
        "saved sequence is not standalone material"
    );
    let plan = adapter::compile(&inputs, &target, "template", &origin, &request).unwrap();
    assert_eq!(plan.steps.len(), 2);
    let Operation::TopoActions { sources, .. } = &plan.steps[1] else {
        panic!("actions")
    };
    assert!(matches!(sources[2].source, Source::Step(0)));
    assert!(
        sources[..2]
            .iter()
            .all(|s| matches!(s.source, Source::Input(_)))
    );
    let checked = v::verify(&plan, &target);
    assert_eq!(
        checked.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        checked.findings
    );
    assert_eq!(checked.products.last().unwrap().sequence.as_str(), TARGET);
    let bindings = vec![
        InputEventBinding {
            input_id: "vt".into(),
            event: event(&w, 2),
        },
        InputEventBinding {
            input_id: "vb".into(),
            event: event(&w, 3),
        },
        InputEventBinding {
            input_id: "template".into(),
            event: origin.selection.event.clone(),
        },
    ];
    let finalrun = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 4),
        bindings,
        &cloning_outputs::Options {
            name: "Reusable literal TOPO circle".into(),
            show_prepared: true,
            save_primers: true,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let loaded =
        workspace::load_workspace(&workspace::save_workspace(&finalrun.workspace).unwrap())
            .unwrap();
    let next = run_material::Selection {
        run_id: finalrun.receipt.run.run_id.clone(),
        step_index: 1,
        event: finalrun.receipt.construct.clone(),
    };
    let fresh = run_material::reconstruct(&loaded, &next).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), TARGET);
    assert!(v::export_protocol_html(&plan, &target).contains("saved"));
    assert!(
        collagene_core::pdf::export_protocol_pdf(&plan, &target)
            .unwrap()
            .bytes
            .starts_with(b"%PDF-")
    );
}
#[test]
fn exact_capsule_and_role_pins_vector_prepared_rules_and_staggered_material_reject() {
    let (_, inputs, target, origin, request) = saved(false);
    let mut wrong = origin.clone();
    wrong.selection.event.event_id = "ev-unrelated".into();
    assert!(adapter::compile(&inputs, &target, "template", &wrong, &request).is_err());
    assert!(adapter::compile(&inputs, &target, "vt", &origin, &request).is_err());
    let mut mixed = request.clone();
    if let topo::MaterialRecipe::Prepared { strands, .. } =
        &mut mixed.preparation.recipe.materials[0]
    {
        strands[0].segments[0].source_role = "template".into();
        strands[0].segments[0].length = 10;
    }
    assert!(adapter::compile(&inputs, &target, "template", &origin, &mixed).is_err());
    let mut changed = inputs.clone();
    changed[2].molecule.name.push_str(" changed");
    assert!(adapter::compile(&changed, &target, "template", &origin, &request).is_err());
    let (_, staggered, _, cap, req) = saved(true);
    let task = adapter::Task::InspectActions {
        request: actions::ActionsRequest {
            preparation: req.preparation,
            processing: req.processing,
        },
    };
    assert!(adapter::inspect(&staggered, "template", &cap, &task).is_err());
    let wrong_target = Molecule::dna(
        "Independent wrong target",
        "AAGGGGCATCCCTTCACCATGACC",
        Topology::Circular,
    )
    .unwrap();
    let plan = adapter::compile(&inputs, &wrong_target, "template", &origin, &request).unwrap();
    assert_eq!(v::verify(&plan, &wrong_target).status, Status::Invalid);
}
fn wire<T: serde::Serialize>(magic: &[u8; 4], version: u16, payload: &T) -> Vec<u8> {
    let payload = postcard::to_allocvec(payload).unwrap();
    let mut b = magic.to_vec();
    b.extend_from_slice(&version.to_le_bytes());
    b.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    b.extend_from_slice(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    b
}
#[test]
fn worker89_wrapper_and_unchanged_inner_bindings_reject_old_worker_and_changed_declarations() {
    let (_, inputs, _, origin, request) = saved(false);
    let task = adapter::Task::InspectActions {
        request: actions::ActionsRequest {
            preparation: request.preparation.clone(),
            processing: request.processing.clone(),
        },
    };
    let revision = adapter::revision(&inputs, "template", &origin, &task).unwrap();
    let request = protocol::Request {
        id: 17,
        revision: revision.clone(),
        job: protocol::AnalysisJob::InspectTopoRun {
            inputs: inputs.clone(),
            template_role: "template".into(),
            origin: origin.clone(),
            task: task.clone(),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 169);
    assert!(protocol::process_request(&wire(b"CGQ\0", 88, &request)).is_err());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert!(protocol::decode_response(&wire(b"CGR\0", 88, &response)).is_err());
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&response.result).unwrap())
            .unwrap()
            .0,
        159
    );
    let protocol::AnalysisResult::TopoRunInspection(inspection) = response.result else {
        panic!("saved template inspection")
    };
    inspection
        .validate_for(&inputs, "template", &origin, &task, &revision)
        .unwrap();
    let adapter::Reply::ActionsInspected(inner) = &inspection.reply else {
        panic!("actions")
    };
    inner
        .validate_binding(
            &adapter::documentary_revision(&inputs, "template", &origin, &task).unwrap(),
        )
        .unwrap();
    assert_ne!(inner.context_revision, revision);
    assert!(
        inspection
            .validate_for(&inputs, "vt", &origin, &task, &revision)
            .is_err()
    );
    let mut changed = inputs.clone();
    changed[0].molecule.name.push('!');
    assert!(
        inspection
            .validate_for(&changed, "template", &origin, &task, &revision)
            .is_err()
    );
    let mut cap = origin.clone();
    cap.selection.run_id.push('x');
    assert!(
        inspection
            .validate_for(&inputs, "template", &cap, &task, &revision)
            .is_err()
    );
    let mut different = task.clone();
    if let adapter::Task::InspectActions { request } = &mut different {
        request.processing.actions.clear();
    }
    assert!(
        inspection
            .validate_for(&inputs, "template", &origin, &different, &revision)
            .is_err()
    );
}
#[test]
fn saved_template_orientation_projection_and_independent_check_share_exact_native_inspection() {
    let (_, inputs, target, origin, request) = saved(false);
    let mut prep = request.preparation.clone();
    prep.transfers.clear();
    let task = adapter::Task::Orientations {
        request: prep,
        vector: [0, 1],
        insert: [2, 3],
        directional: true,
    };
    let out = adapter::inspect(&inputs, "template", &origin, &task).unwrap();
    let adapter::Reply::Orientations(orientations) = out.reply else {
        panic!("orientations")
    };
    assert_eq!(orientations.reference.candidates.len(), 2);
    let task = adapter::Task::ProcessActions {
        request: request.clone(),
    };
    let out = adapter::inspect(&inputs, "template", &origin, &task).unwrap();
    let adapter::Reply::ActionsProcessed(product) = out.reply else {
        panic!("processed")
    };
    assert_eq!(product.product.as_ref().unwrap().sequence.as_str(), TARGET);
    let selection = topo::Selection {
        batch_revision: adapter::preparation_revision(
            &inputs,
            "template",
            &origin,
            &request.preparation,
        )
        .unwrap(),
        component: 0,
        reference_seed: call(0, 0),
    };
    let task = adapter::Task::Project {
        request: request.preparation.clone(),
        selection: selection.clone(),
    };
    adapter::inspect(&inputs, "template", &origin, &task).unwrap();
    let task = adapter::Task::Check {
        request: request.preparation,
        selection,
        target,
        policy: topo::AcceptancePolicy::default(),
    };
    adapter::inspect(&inputs, "template", &origin, &task).unwrap();
}
#[test]
fn embedded_topo_oligos_keep_literal_domains_tails_and_exact_source_without_ordering_prepared_vectors()
 {
    use collagene_core::plan_oligos::{self, Options};
    let (old, target, inputs) = fixture();
    let Operation::TopoActions { mut request, .. } = old.steps[1].clone() else {
        panic!("actions")
    };
    if let topo::MaterialRecipe::ExactPcr {
        forward,
        reverse,
        tail,
        ..
    } = &mut request.preparation.recipe.materials[1]
    {
        *forward = Primer::new("New literal forward", "CACCA", "GG").unwrap();
        *reverse = Primer::new("New literal reverse", "CGTCA", "TT").unwrap();
        *tail = Some(collagene_core::tailing::TailingBase::A);
    }
    let sources = inputs
        .iter()
        .map(|i| TopoSource {
            role: i.role.clone(),
            source: Source::Input(i.role.clone()),
        })
        .collect::<Vec<_>>();
    let selection = topo::Selection {
        batch_revision: "a".repeat(64),
        component: 0,
        reference_seed: call(0, 0),
    };
    let operations = vec![
        Operation::Topo {
            sources: sources.clone(),
            request: request.preparation.clone(),
            selection: selection.clone(),
        },
        Operation::CheckedTopo {
            sources: sources.clone(),
            request: request.preparation.clone(),
            selection,
            policy: topo::AcceptancePolicy::default(),
        },
        Operation::ProcessedTopo {
            sources: sources.clone(),
            request: collagene_core::topo_processed_product::Request {
                treatment: collagene_core::topo_nick_sealing::Request {
                    preparation: request.preparation.clone(),
                    sealing: p::nick_sealing::Request {
                        ligase_assumption: "Explicit documentary ligase declaration".into(),
                        nicks: vec![],
                    },
                },
                selection: p::processed_projection::Selection {
                    sealing_context_revision: "b".repeat(64),
                    component: 0,
                    reference_seed: call(0, 0),
                },
                options: projection::Options::default(),
            },
        },
        Operation::TopoActions { sources, request },
    ];
    for operation in operations {
        let plan = DesignPlan {
            name: "Declared embedded oligos; no physical success asserted".into(),
            target_revision: target.revision(),
            inputs: inputs
                .iter()
                .map(|i| PlanInput {
                    id: i.role.clone(),
                    revision: i.molecule.revision(),
                    molecule: i.molecule.clone(),
                })
                .collect(),
            steps: vec![operation],
        };
        let out = plan_oligos::build(&plan, &target, &Options::default()).unwrap();
        assert!(out.complete, "{:?}", out.unavailable);
        assert_eq!(out.uses.len(), 2);
        assert_eq!(out.oligos.records.len(), 2);
        assert!(out.uses.iter().all(
            |u| matches!(&u.source, Some(Source::Input(id)) if id == "template")
                && u.step_index == 0
        ));
        assert_eq!(
            out.uses[0].primer.as_ref().unwrap().ordered_sequence(),
            "GGCACCA"
        );
        assert_eq!(
            out.uses[1].primer.as_ref().unwrap().ordered_sequence(),
            "TTCGTCA"
        );
        assert_eq!(
            out.uses[0].primer.as_ref().unwrap().annealing.as_str(),
            "CACCA"
        );
        assert_eq!(out.uses[0].primer.as_ref().unwrap().tail, "GG");
        assert!(out.warnings.iter().any(|s| s.contains("post-PCR")));
        assert!(
            out.warnings
                .iter()
                .any(|s| s.contains("not an inferred synthesis order"))
        );
    }
}
#[test]
fn saved_template_orders_only_new_pcr_and_unsupported_prepared_or_missing_roles_cannot_be_flattened()
 {
    use collagene_core::plan_oligos::{self, Options};
    let (_, inputs, target, origin, request) = saved(false);
    let plan = adapter::compile(&inputs, &target, "template", &origin, &request).unwrap();
    let out = plan_oligos::build(&plan, &target, &Options::default()).unwrap();
    assert!(out.complete);
    assert_eq!(out.uses.len(), 2);
    assert!(
        out.uses
            .iter()
            .all(|u| u.step_index == 1 && matches!(&u.source, Some(Source::Step(0))))
    );
    assert_eq!(out.uses[0].primer_name, "F");
    assert_eq!(out.uses[1].primer_name, "R");
    assert!(!out.uses.iter().any(|u| u.primer_name.starts_with("Saved ")));
    assert!(
        out.warnings
            .iter()
            .any(|s| s.contains("Historical origin PCR primers"))
    );
    let mut missing = plan.clone();
    if let Operation::TopoActions { request, .. } = &mut missing.steps[1] {
        if let topo::MaterialRecipe::ExactPcr { template_role, .. } =
            &mut request.preparation.recipe.materials[1]
        {
            *template_role = "missing-role".into();
        }
    }
    assert!(plan_oligos::build(&missing, &target, &Options::default()).is_err());
    let mut derived = plan.clone();
    if let Operation::TopoActions { sources, .. } = &mut derived.steps[1] {
        sources[0].source = Source::Step(0);
    }
    let out = plan_oligos::build(&derived, &target, &Options::default()).unwrap();
    assert!(!out.complete);
    assert!(
        out.unavailable
            .iter()
            .any(|s| s.contains("derived material"))
    );
    assert_eq!(out.uses.len(), 2);
    let mut original_order = plan.clone();
    original_order.inputs[0]
        .molecule
        .metadata
        .insert("oligo:opaque".into(), "unmapped original order".into());
    original_order.inputs[0].revision = original_order.inputs[0].molecule.revision();
    let out = plan_oligos::build(&original_order, &target, &Options::default()).unwrap();
    assert!(!out.complete);
    assert_eq!(out.uses.len(), 2);
}
#[test]
fn actual_reusable_but_search_incomplete_origin_stays_conditional_through_new_topo_pcr() {
    use collagene_core::{
        chemistry::Terminus, homology_assembly, overlap_profiles, verification::OrientedSource,
    };
    let (_, target, mut inputs) = fixture();
    let left = Molecule::dna(
        "Left actual declared stock",
        "CACCATGACGAAAA",
        Topology::Linear,
    )
    .unwrap();
    let right = Molecule::dna("Right actual declared stock", "AAAATTGC", Topology::Linear).unwrap();
    let assembled = Molecule::dna(
        "Independent expected homology18",
        "CACCATGACGAAAATTGC",
        Topology::Linear,
    )
    .unwrap();
    let policy = homology_assembly::Options {
        overlaps: vec![4],
        topology: Topology::Linear,
        min_paired_bases: 4,
        max_comparisons: 106,
        ..Default::default()
    };
    let plan = DesignPlan {
        name: "Literal physical homology with finite terminal search".into(),
        target_revision: assembled.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: left.revision(),
                molecule: left.clone(),
            },
            PlanInput {
                id: "b".into(),
                revision: right.revision(),
                molecule: right.clone(),
            },
        ],
        steps: vec![
            Operation::DeclareSourceTermini {
                template: Source::Input("a".into()),
                left: Terminus::PCR,
                right: Terminus::PCR,
            },
            Operation::DeclareSourceTermini {
                template: Source::Input("b".into()),
                left: Terminus::PCR,
                right: Terminus::PCR,
            },
            Operation::HomologyAssembly {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                ],
                options: overlap_profiles::ProfileOptions::default(),
                policy,
            },
        ],
    };
    let report = v::verify(&plan, &assembled);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(
        report.products.last().unwrap().sequence.as_str(),
        "CACCATGACGAAAATTGC"
    );
    let w = Workspace::from_documents(
        "Actual conditional origin plus unchanged TOPO sources",
        vec![
            left,
            right,
            assembled,
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target.clone(),
        ],
        5,
    )
    .unwrap();
    let prior = cloning_outputs::finalize(
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
        &cloning_outputs::Options {
            name: "Conditional but reusable physical template".into(),
            show_prepared: true,
            save_primers: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let origin = run_material::capsule(
        &prior.workspace,
        &run_material::Selection {
            run_id: prior.receipt.run.run_id.clone(),
            step_index: 2,
            event: prior.receipt.construct.clone(),
        },
    )
    .unwrap();
    inputs[2].molecule = prior.workspace.revisions[&origin.selection.event.revision].clone();
    let (old, _, _) = fixture();
    let Operation::TopoActions { mut request, .. } = old.steps[1].clone() else {
        panic!("actions")
    };
    let task = adapter::Task::InspectActions {
        request: actions::ActionsRequest {
            preparation: request.preparation.clone(),
            processing: request.processing.clone(),
        },
    };
    let inspected = adapter::inspect(&inputs, "template", &origin, &task).unwrap();
    assert!(inspected.origin_incomplete);
    assert!(!inspected.origin_findings.is_empty());
    assert_eq!(inspected.origin_findings_total, report.findings.len());
    assert!(!inspected.origin_findings_truncated);
    let adapter::Reply::ActionsInspected(actual) = inspected.reply else {
        panic!("actions")
    };
    request.selection.processing_context_revision = actual.processing.context_revision;
    let newplan = adapter::compile(&inputs, &target, "template", &origin, &request).unwrap();
    let fresh = v::verify(&newplan, &target);
    assert_eq!(fresh.status, Status::Incomplete, "{:?}", fresh.findings);
    assert_eq!(fresh.products.last().unwrap().sequence.as_str(), TARGET);
    let bindings = vec![
        InputEventBinding {
            input_id: "vt".into(),
            event: event(&prior.workspace, 3),
        },
        InputEventBinding {
            input_id: "vb".into(),
            event: event(&prior.workspace, 4),
        },
        InputEventBinding {
            input_id: "template".into(),
            event: origin.selection.event.clone(),
        },
    ];
    let mut options = cloning_outputs::Options {
        name: "Explicit conditional TOPO design".into(),
        show_prepared: true,
        save_primers: true,
        allow_incomplete: false,
    };
    assert!(
        cloning_outputs::finalize(
            &prior.workspace,
            &newplan,
            event(&prior.workspace, 5),
            bindings.clone(),
            &options
        )
        .is_err()
    );
    options.allow_incomplete = true;
    let saved = cloning_outputs::finalize(
        &prior.workspace,
        &newplan,
        event(&prior.workspace, 5),
        bindings,
        &options,
    )
    .unwrap();
    assert!(saved.receipt.primer_attachment.is_some());
}
