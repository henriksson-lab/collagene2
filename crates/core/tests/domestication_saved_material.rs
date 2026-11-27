//! Real saved Run/Derived templates; independent genotypes precede native preparation.
use collagene_core::{
    Molecule, Sequence, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::{EndTreatment, Strand, Terminus},
    cloning_outputs, codon_design, domestication as design, domestication_preparation as prep,
    material_derivation as derivation, overlap_extension_primer_design as oe, plan_oligos,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGAATTCCAGACCTGTAGATCGTCA";
const TARGET: &str = "AGTCCGATGCTACGTAGACTTCCAGACCTGTAGATCGTCA";
const SCATTERED: &str = "ACGTTGCACTGATCGAGTCCGGTCTCATGCTACGTAGGACGATCCAGTACGCTTATTCAGACCTGTAGATCGTCAGGTCTCTGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC";
const SCATTERED_TARGET: &str = "ACGTTGCACTGATCGAGTCCGGCCTCATGCTACGTAGGACGATCCAGTACGCTTATTCAGACCTGTAGATCGTCAGGCCTCTGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC";
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn input(molecule: Molecule) -> PlanInput {
    PlanInput {
        id: prep::INPUT_ID.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn stock(dna: &str, topology: Topology) -> (Workspace, MaterialOrigin) {
    let mut source = Molecule::dna("Historical stock source", dna, Topology::Linear).unwrap();
    source.features.push(collagene_core::Feature {
        kind: "misc_feature".into(),
        location: collagene_core::Location::range(24, 30),
        qualifiers: vec![("label".into(), Some("Retained native source region".into()))],
        needs_review: false,
    });
    let target = Molecule::dna("Independent historical stock requirement", dna, topology).unwrap();
    let mut w = Workspace::from_documents(
        "Native source ancestry",
        vec![source.clone(), target.clone()],
        0,
    )
    .unwrap();
    // The long stock repeats its final8-base word;18-base domains disambiguate the actual whole-source amplification.
    let domain = if dna == SCATTERED { 18 } else { 8 };
    let mut steps = vec![Operation::Pcr {
        template: Source::Input("stock".into()),
        forward: Primer::new("Historical F", &dna[..domain], "").unwrap(),
        reverse: Primer::new(
            "Historical R",
            Sequence::dna(&dna[dna.len() - domain..])
                .unwrap()
                .reverse_complement()
                .unwrap()
                .as_str(),
            "",
        )
        .unwrap(),
        template_strand: Strand::Top,
    }];
    if topology == Topology::Circular {
        steps.push(Operation::TreatEnds {
            template: Source::Step(0),
            treatment: EndTreatment::Phosphorylate,
        });
        steps.push(Operation::Ligate {
            fragments: vec![OrientedSource {
                source: Source::Step(1),
                reverse: false,
            }],
            circular: true,
            allow_nicked: false,
            options: Default::default(),
        });
    }
    let plan = DesignPlan {
        name: "Actual upstream PCR stock".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps,
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved native stock".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    w = workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let capsule = run_material::capsule(
        &w,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id,
            step_index: plan.steps.len() - 1,
            event: saved.receipt.construct,
        },
    )
    .unwrap();
    (w, MaterialOrigin::Run(capsule))
}
fn derived(w: &Workspace, origin: &MaterialOrigin) -> (Workspace, MaterialOrigin) {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRV")
        .unwrap()
        .id
        .clone();
    let p = derivation::prepare_install(
        w,
        &derivation::Request {
            source: origin.event().clone(),
            origin: derivation::SourceOrigin::Saved(origin.clone()),
            derivation: derivation::Derivation::CatalogDigest {
                selection: catalog_digest::select(&catalog, &[id]).unwrap(),
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Native no-site derivative".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let w = workspace::load_workspace(&workspace::save_workspace(&p.workspace).unwrap()).unwrap();
    let origin = MaterialOrigin::Derived(derivation::capsule(&w, &p.receipt.product).unwrap());
    (w, origin)
}
fn options() -> oe::Options {
    let mut o = oe::Options {
        homology_flank: 8,
        min_paired_bases: 8,
        ..Default::default()
    };
    o.primer_options.min_length = 8;
    o.primer_options.max_length = 8;
    o.primer_options.min_gc_fraction = 0.;
    o.primer_options.max_gc_fraction = 1.;
    o.primer_options.max_tm_difference_c = 100.;
    o.primer_options.max_homopolymer = 200;
    o.primer_options.max_three_prime_complementarity = 200;
    o
}
fn request(
    source: &Molecule,
    origin: MaterialOrigin,
    enzyme: &str,
    positions: &[usize],
    route: prep::Route,
) -> prep::Request {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == enzyme)
        .unwrap()
        .id
        .clone();
    let design = design::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(&catalog, &[id]).unwrap(),
        scope: vec![design::SourceRange {
            start: 0,
            end: source.sequence.len(),
        }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: Default::default(),
    };
    let patches = positions
        .iter()
        .map(|&p| design::Patch {
            position: p,
            expected: source.sequence.as_str()[p..p + 1].into(),
            replacement: "C".into(),
        })
        .collect::<Vec<_>>();
    let selection = design::Selection {
        context_revision: design::propose_revision(source, &design).unwrap(),
        candidate_id: design::candidate_id(&design, &patches).unwrap(),
        patches,
    };
    prep::Request {
        design,
        selection,
        route,
        origin: Some(origin),
    }
}
fn inverse(position: usize) -> prep::Route {
    prep::Route::InversePcr {
        edit_start: position,
        edit_end: position + 1,
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 1,
    }
}
#[test]
fn run_and_derived_inverse_outputs_keep_independent_target_and_only_new_primers() {
    let (w, run) = stock(ORIGINAL, Topology::Circular);
    let (d, derived) = derived(&w, &run);
    for (w, origin) in [(w, run), (d, derived)] {
        let before = workspace::save_workspace(&w).unwrap();
        let source = w.revisions[&origin.event().revision].clone();
        let r = request(&source, origin.clone(), "EcoRI", &[18], inverse(18));
        let p = prep::materialize(&input(source.clone()), &r, Some(&w)).unwrap();
        assert_eq!(p.target.sequence.as_str(), TARGET);
        assert!(p.proposal_ready);
        assert!(
            p.target
                .features
                .iter()
                .any(|f| f.label() == "Retained native source region")
        );
        assert!(
            p.verification
                .products
                .last()
                .unwrap()
                .features
                .iter()
                .any(|f| f.label() == "Retained native source region")
        );
        assert_eq!(p.plan.steps.len(), 4);
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert!(matches!(
            &p.plan.steps[1],
            Operation::Pcr {
                template: Source::Step(0),
                ..
            }
        ));
        p.validate_binding(&p.context_revision).unwrap();
        let orders = plan_oligos::build(&p.plan, &p.target, &Default::default()).unwrap();
        assert!(orders.complete);
        assert_eq!(orders.uses.len(), 2);
        assert!(
            orders
                .uses
                .iter()
                .all(|u| u.step_index == 1 && !u.primer_name.contains("Historical"))
        );
        let mut frame = w.clone();
        let declared = collagene_core::history::DocumentHistory::new(
            "independent-domestication-target",
            p.target.clone(),
        )
        .unwrap();
        let target = EventBinding {
            document_id: declared.id.clone(),
            event_id: declared.head_event().into(),
            revision: p.target.revision(),
        };
        declared.contribute(&mut frame).unwrap();
        let saved = cloning_outputs::finalize(
            &frame,
            &p.plan,
            target,
            vec![InputEventBinding {
                input_id: prep::INPUT_ID.into(),
                event: origin.event().clone(),
            }],
            &cloning_outputs::Options {
                name: "New domesticated circle".into(),
                save_primers: true,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let reopened =
            workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap())
                .unwrap();
        for original in &w.documents {
            let retained = reopened
                .documents
                .iter()
                .find(|d| d.id == original.id)
                .unwrap();
            assert_eq!(retained.head, original.head);
            assert_eq!(retained.head_event, original.head_event);
        }
        let replay = run_lineage::replay_run(&reopened, &saved.receipt.run.run_id).unwrap();
        assert!(replay.products_match);
        assert_eq!(
            replay.current_report.status,
            Status::ConsistentUnderAssumptions
        );
        assert!(saved.receipt.primer_attachment.is_some());
        assert!(verification::export_protocol_html(&p.plan, &p.target).contains("Reconstruct"));
        assert!(
            collagene_core::pdf::export_protocol_pdf(&p.plan, &p.target)
                .unwrap()
                .bytes
                .starts_with(b"%PDF-")
        );
        assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    }
}
#[test]
fn saved_linear_and_circular_oe_use_actual_reconstruction_prefix_and_all_six_new_oligos() {
    for topology in [Topology::Linear, Topology::Circular] {
        let (w, run) = stock(ORIGINAL, topology);
        let (w, origin) = derived(&w, &run);
        let source = w.revisions[&origin.event().revision].clone();
        let r = request(
            &source,
            origin,
            "EcoRI",
            &[18],
            prep::Route::OverlapExtension { options: options() },
        );
        let p = prep::materialize(&input(source), &r, Some(&w)).unwrap();
        assert_eq!(p.target.sequence.as_str(), TARGET);
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert!(p.proposal_ready);
        assert_eq!(
            p.plan.steps.len(),
            if topology == Topology::Circular { 7 } else { 5 }
        );
        for index in [1, 2] {
            assert!(matches!(
                &p.plan.steps[index],
                Operation::Pcr {
                    template: Source::Step(0),
                    ..
                }
            ));
        }
        assert!(matches!(
            &p.plan.steps[4],
            Operation::Pcr {
                template: Source::Step(3),
                ..
            }
        ));
        p.validate_binding(&p.context_revision).unwrap();
        let orders = plan_oligos::build(&p.plan, &p.target, &Default::default()).unwrap();
        assert!(orders.complete);
        assert_eq!(orders.uses.len(), 6);
        assert!(
            orders
                .uses
                .iter()
                .all(|u| !u.primer_name.contains("Historical"))
        );
        let mut changed = p.clone();
        let prep::MethodEvidence::OverlapExtension { source_frame, .. } = &mut changed.method
        else {
            panic!()
        };
        source_frame.sequence = Sequence::dna("AAAAAAAA").unwrap();
        assert!(changed.validate_binding(&p.context_revision).is_err());
    }
}
#[test]
fn saved_rounds_keep_original_coordinates_and_consume_previous_actual_closed_products() {
    let (w, origin) = stock(SCATTERED, Topology::Circular);
    let source = w.revisions[&origin.event().revision].clone();
    let rounds = vec![
        prep::InverseRound {
            edit_start: 22,
            edit_end: 23,
            forward_length: 8,
            reverse_length: 8,
            reverse_tail_bases: 1,
        },
        prep::InverseRound {
            edit_start: 77,
            edit_end: 78,
            forward_length: 8,
            reverse_length: 8,
            reverse_tail_bases: 0,
        },
    ];
    let r = request(
        &source,
        origin,
        "BsaI",
        &[22, 77],
        prep::Route::InversePcrRounds { rounds },
    );
    let p = prep::materialize(&input(source), &r, Some(&w)).unwrap();
    assert_eq!(p.target.sequence.as_str(), SCATTERED_TARGET);
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(p.plan.steps.len(), 7);
    assert!(matches!(
        &p.plan.steps[4],
        Operation::Pcr {
            template: Source::Step(3),
            ..
        }
    ));
    let prep::MethodEvidence::InversePcrRounds { rounds } = &p.method else {
        panic!()
    };
    assert_eq!(rounds[1].source_origin, 23);
    assert_eq!(rounds[1].options.edit_start, 54);
    p.validate_binding(&p.context_revision).unwrap();
    let mut changed = p.clone();
    let Operation::Pcr { template, .. } = &mut changed.plan.steps[4] else {
        panic!()
    };
    *template = Source::Step(0);
    changed.verification.plan_revision = verification::plan_revision(&changed.plan);
    assert!(changed.validate_binding(&p.context_revision).is_err());
}
#[test]
fn owning_workspace_source_alias_and_response_tampering_never_become_material() {
    let (w, origin) = stock(ORIGINAL, Topology::Circular);
    let source = w.revisions[&origin.event().revision].clone();
    let r = request(&source, origin, "EcoRI", &[18], inverse(18));
    let i = input(source.clone());
    assert!(prep::revision(&i, &r, None).is_err());
    let foreign =
        Workspace::from_documents("Same DNA lacks native ownership", vec![source], 0).unwrap();
    assert!(prep::preflight(&i, &r, Some(&foreign)).is_err());
    let mut missing = r.clone();
    missing.origin = None;
    assert!(prep::materialize(&i, &missing, Some(&w)).is_err());
    let mut alias = r.clone();
    let Some(MaterialOrigin::Run(o)) = &mut alias.origin else {
        panic!()
    };
    o.selection.event.event_id.push_str("-alias");
    assert!(prep::revision(&i, &alias, Some(&w)).is_err());
    let p = prep::materialize(&i, &r, Some(&w)).unwrap();
    let mut changed = p.clone();
    changed.target = Molecule::dna("Wrong intended target", ORIGINAL, Topology::Circular).unwrap();
    assert!(changed.validate_binding(&p.context_revision).is_err());
    let mut json = serde_json::to_value(&r).unwrap();
    json.as_object_mut().unwrap().remove("origin");
    assert!(serde_json::from_value::<prep::Request>(json).is_err());
}
#[test]
fn actual_nicked_saved_circle_does_not_gain_a_fabricated_pcr_template() {
    let a = Molecule::dna("Actual left half", &ORIGINAL[..20], Topology::Linear).unwrap();
    let b = Molecule::dna("Actual right half", &ORIGINAL[20..], Topology::Linear).unwrap();
    let target = Molecule::dna("Independent nicked40", ORIGINAL, Topology::Circular).unwrap();
    let w = Workspace::from_documents(
        "Actual nicked source",
        vec![a.clone(), b.clone(), target.clone()],
        0,
    )
    .unwrap();
    let pair = |id: &str, s: &str| Operation::Pcr {
        template: Source::Input(id.into()),
        forward: Primer::new("Historical F", &s[..8], "").unwrap(),
        reverse: Primer::new(
            "Historical R",
            Sequence::dna(&s[s.len() - 8..])
                .unwrap()
                .reverse_complement()
                .unwrap()
                .as_str(),
            "",
        )
        .unwrap(),
        template_strand: Strand::Top,
    };
    let plan = DesignPlan {
        name: "One phosphorylated half only".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: a.revision(),
                molecule: a,
            },
            PlanInput {
                id: "b".into(),
                revision: b.revision(),
                molecule: b,
            },
        ],
        steps: vec![
            pair("a", &ORIGINAL[..20]),
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            pair("b", &ORIGINAL[20..]),
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
    };
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
        &cloning_outputs::Options {
            name: "Actual nicked source circle".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let w = saved.workspace;
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &w,
            &run_material::Selection {
                run_id: saved.receipt.run.run_id,
                step_index: 3,
                event: saved.receipt.construct,
            },
        )
        .unwrap(),
    );
    let source = w.revisions[&origin.event().revision].clone();
    assert!(
        !origin
            .replay(&source)
            .unwrap()
            .tracked()
            .duplex
            .nicks
            .is_empty()
    );
    let before = workspace::save_workspace(&w).unwrap();
    let r = request(&source, origin.clone(), "EcoRI", &[18], inverse(18));
    let p = prep::materialize(&input(source.clone()), &r, Some(&w)).unwrap();
    assert_eq!(p.verification.status, Status::Incomplete);
    assert_eq!(p.verification.products.len(), 1);
    assert!(
        p.verification
            .findings
            .iter()
            .any(|f| f.to_lowercase().contains("nick"))
    );
    let r = request(
        &source,
        origin,
        "EcoRI",
        &[18],
        prep::Route::OverlapExtension { options: options() },
    );
    assert!(prep::materialize(&input(source), &r, Some(&w)).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn genuine_incomplete_origin_stays_conditional_after_fresh_inverse_pcr() {
    const SOURCE: &str = "AGTCCGAATTCCACCATGACGAAAATTGCTCGGATCGTCA";
    const INTENDED: &str = "AGTCCGACTTCCACCATGACGAAAATTGCTCGGATCGTCA";
    let a = Molecule::dna("Declared first stock", "CACCATGACGAAAA", Topology::Linear).unwrap();
    let b = Molecule::dna("Declared second stock", "AAAATTGC", Topology::Linear).unwrap();
    let target = Molecule::dna("Independent final stock40", SOURCE, Topology::Circular).unwrap();
    let w = Workspace::from_documents(
        "Finite actual homology search",
        vec![a.clone(), b.clone(), target.clone()],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual bounded homology then new stock PCR".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: a.revision(),
                molecule: a,
            },
            PlanInput {
                id: "b".into(),
                revision: b.revision(),
                molecule: b,
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
                options: Default::default(),
                policy: collagene_core::homology_assembly::Options {
                    overlaps: vec![4],
                    topology: Topology::Linear,
                    min_paired_bases: 4,
                    max_comparisons: 106,
                    ..Default::default()
                },
            },
            Operation::Pcr {
                template: Source::Step(2),
                forward: Primer::new("Historical template F", "CACCATGA", "AGTCCGAATTC").unwrap(),
                reverse: Primer::new("Historical template R", "GCAATTTT", "TGACGATCCGA").unwrap(),
                template_strand: Strand::Top,
            },
            Operation::TreatEnds {
                template: Source::Step(3),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(4),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.products.len(), 6);
    assert!(verification::plan_target_matches(
        &plan,
        report.products.last().unwrap(),
        &target
    ));
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
        &cloning_outputs::Options {
            name: "Conditional but actual circle".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let w =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &w,
            &run_material::Selection {
                run_id: saved.receipt.run.run_id,
                step_index: 5,
                event: saved.receipt.construct,
            },
        )
        .unwrap(),
    );
    let source = w.revisions[&origin.event().revision].clone();
    let fresh = origin.replay(&source).unwrap();
    assert!(fresh.incomplete());
    assert!(fresh.tracked().duplex.nicks.is_empty());
    let r = request(&source, origin, "EcoRI", &[7], inverse(7));
    let p = prep::materialize(&input(source), &r, Some(&w)).unwrap();
    assert_eq!(p.target.sequence.as_str(), INTENDED);
    assert!(p.proposal_ready);
    assert_eq!(
        p.verification.status,
        Status::Incomplete,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(p.verification.products.len(), 4);
    assert!(verification::plan_target_matches(
        &p.plan,
        p.verification.products.last().unwrap(),
        &p.target
    ));
    p.validate_binding(&p.context_revision).unwrap();
    let orders = plan_oligos::build(&p.plan, &p.target, &Default::default()).unwrap();
    assert!(orders.complete);
    assert_eq!(orders.uses.len(), 2);
    assert!(
        orders
            .uses
            .iter()
            .all(|u| !u.primer_name.contains("Historical"))
    );
}
