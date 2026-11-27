//! Literal independent targets and native saved-template ancestry, not imported summaries.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::{EndTreatment, FivePrime},
    history::DocumentHistory,
    material_derivation as derived,
    modified_duplex::{self, ModificationSummary},
    mutagenesis::{self, InversePcrOptions},
    pdf, plan_oligos,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
const REPLACED: &str = "GGAGATGCTACGTAGGACTTCAGACCTGTAAGT";
fn dna(name: &str, bases: &str, topology: Topology) -> Molecule {
    Molecule::dna(name, bases, topology).unwrap()
}
fn pin(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn append(w: &mut Workspace, id: &str, molecule: Molecule) -> EventBinding {
    DocumentHistory::new(id, molecule)
        .unwrap()
        .contribute(w)
        .unwrap();
    pin(w, w.documents.len() - 1)
}
fn restart(w: &Workspace) -> Workspace {
    workspace::load_workspace(&workspace::save_workspace(w).unwrap()).unwrap()
}
fn options() -> InversePcrOptions {
    InversePcrOptions {
        edit_start: 3,
        edit_end: 5,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 1,
    }
}
fn fixture(nicked: bool) -> (Workspace, MaterialOrigin) {
    let mut left = dna("Original left", "AGTCCGATGCTACGTA", Topology::Linear);
    left.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(12, 16),
        qualifiers: vec![("label".into(), Some("Retained literal CGTA".into()))],
        needs_review: false,
    });
    let right = dna("Original right", "GGACTTCAGACCTGTA", Topology::Linear);
    let target = dna("Independent original circle", ORIGINAL, Topology::Circular);
    let w = Workspace::from_documents(
        "Native circle, no invented stock",
        vec![left.clone(), right.clone(), target.clone()],
        2,
    )
    .unwrap();
    let mut steps = vec![
        Operation::Pcr {
            template: Source::Input("left".into()),
            forward: Primer::new("Original LF", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Original LR", "TACGTAGC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        Operation::TreatEnds {
            template: Source::Step(0),
            treatment: EndTreatment::Phosphorylate,
        },
        Operation::Pcr {
            template: Source::Input("right".into()),
            forward: Primer::new("Original RF", "GGACTTCA", "").unwrap(),
            reverse: Primer::new("Original RR", "TACAGGTC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
    ];
    let right_step = if nicked {
        2
    } else {
        steps.push(Operation::TreatEnds {
            template: Source::Step(2),
            treatment: EndTreatment::Phosphorylate,
        });
        3
    };
    let ligation = collagene_core::ligation_outcomes::Options::default();
    steps.push(Operation::Ligate {
        fragments: vec![
            OrientedSource {
                source: Source::Step(1),
                reverse: false,
            },
            OrientedSource {
                source: Source::Step(right_step),
                reverse: false,
            },
        ],
        circular: true,
        allow_nicked: nicked,
        options: ligation,
    });
    let plan = DesignPlan {
        name: "Actual prepared circular template".into(),
        target_revision: target.revision(),
        inputs: vec![input("left", left), input("right", right)],
        steps,
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.products.len(),
        plan.steps.len(),
        "{:?}",
        report.findings
    );
    assert!(verification::plan_target_matches(
        &plan,
        report.products.last().unwrap(),
        &target
    ));
    assert_ne!(report.status, Status::Invalid, "{:?}", report.findings);
    let commit = run_lineage::record_run(
        &w,
        &plan,
        pin(&w, 2),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: pin(&w, 0),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: pin(&w, 1),
            },
        ],
        nicked,
    )
    .unwrap();
    let product = commit.receipt.products.last().unwrap();
    let selection = run_material::Selection {
        run_id: commit.receipt.run_id,
        step_index: product.step_index,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let saved = restart(&commit.workspace);
    let origin = MaterialOrigin::Run(run_material::capsule(&saved, &selection).unwrap());
    (saved, origin)
}
fn source(w: &Workspace, origin: &MaterialOrigin) -> Molecule {
    w.revisions[&origin.event().revision].clone()
}
fn prepare(
    w: &Workspace,
    origin: &MaterialOrigin,
    target: &Molecule,
    o: &InversePcrOptions,
) -> mutagenesis::InversePcrDesign {
    origin.validate_for_workspace(w, origin.event()).unwrap();
    let source = source(w, origin);
    let context = mutagenesis::inverse_pcr_revision(&source, target, Some(origin), o).unwrap();
    let design = mutagenesis::design_inverse_pcr(&source, target, Some(origin), o).unwrap();
    design.validate_binding(&context).unwrap();
    assert_eq!(design.context_revision, context);
    design
}
#[test]
fn actual_saved_circle_supports_literal_insertion_deletion_and_replacement() {
    let (w, origin) = fixture(false);
    let before = workspace::save_workspace(&w).unwrap();
    for (start, end, inserted, split, expected) in [
        (3, 5, "GGA", 1, REPLACED),
        (3, 5, "", 0, "GATGCTACGTAGGACTTCAGACCTGTAAGT"),
        (10, 10, "AACG", 2, "AGTCCGATGCAACGTACGTAGGACTTCAGACCTGTA"),
    ] {
        let o = InversePcrOptions {
            edit_start: start,
            edit_end: end,
            inserted: inserted.into(),
            reverse_tail_bases: split,
            ..options()
        };
        let target = dna("Independent mutation", expected, Topology::Circular);
        let design = prepare(&w, &origin, &target, &o);
        assert_eq!(
            design.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            design.verification.findings
        );
        assert_eq!(design.plan.steps.len(), 4);
        assert!(
            matches!(&design.plan.steps[0],Operation::ReconstructRunProduct{template:Source::Input(id),..} if id=="template")
        );
        assert!(matches!(
            &design.plan.steps[1],
            Operation::Pcr {
                template: Source::Step(0),
                ..
            }
        ));
        assert!(verification::plan_target_matches(
            &design.plan,
            design.verification.products.last().unwrap(),
            &target
        ));
        let feature = design
            .verification
            .products
            .last()
            .unwrap()
            .features
            .iter()
            .find(|f| f.label() == "Retained literal CGTA")
            .unwrap();
        assert_eq!(
            feature
                .location
                .extract(&design.verification.products.last().unwrap().sequence)
                .unwrap(),
            "CGTA"
        );
        let mods: ModificationSummary = serde_json::from_str(
            &design.verification.products[1].metadata["tracked_modifications"],
        )
        .unwrap();
        assert_eq!(mods.known_modified, 0);
        assert_eq!(mods.unknown_eligible, 0);
    }
    assert_eq!(before, workspace::save_workspace(&w).unwrap());
}
#[test]
fn cut_free_derived_circle_survives_archive_and_is_replayed_before_new_pcr() {
    let (w, run) = fixture(false);
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    let selection = catalog_digest::select(&catalog, &[enzyme.id.clone()]).unwrap();
    let original = source(&w, &run);
    assert!(!ORIGINAL.contains("GAATTC"));
    let origin = derived::SourceOrigin::Saved(run.clone());
    let preview = derived::inspect_digest(
        "circle",
        &original,
        &origin,
        &selection,
        &Default::default(),
    )
    .unwrap();
    assert!(preview.inspection.products_determined);
    assert_eq!(preview.inspection.fragments.len(), 1);
    let request = derived::Request {
        source: run.event().clone(),
        origin,
        derivation: derived::Derivation::CatalogDigest {
            selection,
            options: Default::default(),
            fragment_index: 0,
        },
        name: "Actually uncut native circle".into(),
        save_primers: false,
    };
    let installed = derived::prepare_install(&w, &request).unwrap();
    let w = restart(&installed.workspace);
    let origin = MaterialOrigin::Derived(derived::capsule(&w, &installed.receipt.product).unwrap());
    assert_eq!(source(&w, &origin).topology, Topology::Circular);
    let target = dna(
        "Independent derived-template replacement",
        REPLACED,
        Topology::Circular,
    );
    let design = prepare(&w, &origin, &target, &options());
    assert_eq!(
        design.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        design.verification.findings
    );
    assert!(matches!(
        &design.plan.steps[0],
        Operation::ReconstructMaterial { .. }
    ));
    assert_eq!(
        design
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "GAGATGCTACGTAGGACTTCAGACCTGTAAGTG"
    );
}
#[test]
fn new_run_exports_new_oligos_protocol_and_replays_without_changing_original_heads() {
    let (mut w, origin) = fixture(false);
    let old_run = run_lineage::list_runs(&w).unwrap();
    let old_heads = w
        .documents
        .iter()
        .map(|d| (d.id.clone(), d.head.clone(), d.head_event.clone()))
        .collect::<Vec<_>>();
    let target = dna(
        "Reviewed independently supplied replacement",
        REPLACED,
        Topology::Circular,
    );
    let target_pin = append(&mut w, "independent-mutant-target", target.clone());
    let design = prepare(&w, &origin, &target, &options());
    assert_eq!(design.forward.tail, "GA");
    assert_eq!(design.reverse.tail, "C");
    let orders = plan_oligos::build(&design.plan, &target, &Default::default()).unwrap();
    assert_eq!(orders.oligos.records.len(), 2);
    assert_eq!(orders.uses.len(), 2);
    assert!(orders.uses.iter().all(|u| u.step_index == 1));
    let html = verification::export_protocol_html(&design.plan, &target);
    assert!(html.contains("saved") && html.contains("Mutagenesis forward"));
    let pdf = pdf::export_protocol_pdf(&design.plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, design.verification.status);
    let saved = run_lineage::record_run(
        &w,
        &design.plan,
        target_pin,
        vec![InputEventBinding {
            input_id: "template".into(),
            event: origin.event().clone(),
        }],
        false,
    )
    .unwrap();
    let w = restart(&saved.workspace);
    for (id, head, event) in old_heads {
        let d = w.documents.iter().find(|d| d.id == id).unwrap();
        assert_eq!(d.head, head);
        assert_eq!(d.head_event, event);
    }
    for old in old_run {
        let retained = run_lineage::load_run(&w, &old.run_id).unwrap();
        assert_eq!(retained.assessment.status, old.recorded_status);
    }
    let replay = run_lineage::replay_run(&w, &saved.receipt.run_id).unwrap();
    assert!(replay.products_match && replay.engine_matches && replay.status_matches);
    let product = saved.receipt.products.last().unwrap();
    let selected = run_material::Selection {
        run_id: saved.receipt.run_id,
        step_index: product.step_index,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let fresh = MaterialOrigin::Run(run_material::capsule(&w, &selected).unwrap())
        .replay(&w.revisions[&product.revision])
        .unwrap();
    assert!(fresh.tracked().duplex.nicks.is_empty());
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .known_modified,
        0
    );
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .unknown_eligible,
        0
    );
    // The explicit phosphorylated linear intermediate is real, not inferred from final DNA.
    let linear: collagene_core::chemistry::DuplexSummary =
        serde_json::from_str(&replay.current_report.products[2].metadata["duplex_model"]).unwrap();
    let (left, right) = linear.ends.unwrap();
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
}
#[test]
fn wrong_target_and_tampered_documentary_receipts_cannot_inherit_success() {
    let (mut w, origin) = fixture(false);
    let source = source(&w, &origin);
    let alias = append(
        &mut w,
        "same-native-snapshot-different-event",
        source.clone(),
    );
    assert_eq!(alias.revision, origin.event().revision);
    assert_ne!(alias.event_id, origin.event().event_id);
    assert!(origin.validate_for_workspace(&w, &alias).is_err());
    let wrong = dna(
        "Independent incorrect target",
        "GGAGATGCTACGTAGGACTTCAGACCTGTCAGT",
        Topology::Circular,
    );
    let report = prepare(&w, &origin, &wrong, &options());
    assert_eq!(report.verification.status, Status::Invalid);
    assert_eq!(report.verification.products.len(), 4);
    let target = dna("Correct independent target", REPLACED, Topology::Circular);
    let design = prepare(&w, &origin, &target, &options());
    let context = design.context_revision.clone();
    let mut changed = design.clone();
    changed.forward.tail = "TT".into();
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = design.clone();
    changed.plan.steps[1] = Operation::Pcr {
        template: Source::Input("template".into()),
        forward: changed.forward.clone(),
        reverse: changed.reverse.clone(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = design.clone();
    changed.options.edit_end = 6;
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = design.clone();
    changed.verification.target_revision = wrong.revision();
    assert!(changed.validate_binding(&context).is_err());
    assert_ne!(
        mutagenesis::inverse_pcr_revision(&source, &wrong, Some(&origin), &options()).unwrap(),
        context
    );
    let mut damaged = origin.clone();
    if let MaterialOrigin::Run(c) = &mut damaged {
        c.selection.event.event_id = "substituted".into();
    }
    assert!(damaged.validate_for_workspace(&w, origin.event()).is_err());
}
#[test]
fn actual_nicked_circle_is_not_silently_repaired_or_amplified() {
    let (w, origin) = fixture(true);
    let original = source(&w, &origin);
    let before = original.clone();
    assert!(
        !origin
            .replay(&original)
            .unwrap()
            .tracked()
            .duplex
            .nicks
            .is_empty()
    );
    let target = dna("Independent mutation", REPLACED, Topology::Circular);
    let design = prepare(&w, &origin, &target, &options());
    assert_eq!(design.verification.status, Status::Incomplete);
    assert_eq!(
        design.verification.products.len(),
        1,
        "only actual reconstructed source, no PCR product"
    );
    assert!(
        design
            .verification
            .findings
            .iter()
            .any(|s| s.contains("nick propagation"))
    );
    assert_eq!(source(&w, &origin), before);
}
#[test]
fn reusable_capped_origin_keeps_incomplete_status_through_new_material() {
    // Alternative ligation enumeration does not make the intended reaction
    // incomplete. Instead cap a real physical-homology terminal search only
    // after the declared sealed circle is reconstructed from both strands.
    let left = dna("Homology left", "AGTCCGATGCTACGTAGGAC", Topology::Linear);
    let right = dna("Homology right", "GGACTTCAGACCTGTAAGTC", Topology::Linear);
    let original_target = dna("Independent original circle", ORIGINAL, Topology::Circular);
    let w = Workspace::from_documents(
        "Conditional native circular origin",
        vec![left.clone(), right.clone(), original_target.clone()],
        2,
    )
    .unwrap();
    let mut profile = collagene_core::overlap_profiles::ProfileOptions::default();
    profile.scan_repeated_overlaps = false;
    let plan = DesignPlan {
        name: "Physical homology with bounded terminal alternatives".into(),
        target_revision: original_target.revision(),
        inputs: vec![input("left", left), input("right", right)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("left".into()),
                forward: Primer::new("LF", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("LR", "GTCCTACG", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("right".into()),
                forward: Primer::new("RF", "GGACTTCA", "").unwrap(),
                reverse: Primer::new("RR", "GACTTACA", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
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
                options: profile,
                policy: collagene_core::homology_assembly::Options {
                    overlaps: vec![4, 4],
                    topology: Topology::Circular,
                    min_paired_bases: 4,
                    // Three strand checks per input base, paired product, and
                    // both declared 4bp joins; alternatives exceed this budget.
                    max_comparisons: 3 * (20 + 20) + 2 * 32 + 4 + 4,
                    max_candidates: 100,
                    ..Default::default()
                },
            },
        ],
    };
    let checked = verification::verify(&plan, &original_target);
    assert_eq!(checked.status, Status::Incomplete, "{:?}", checked.findings);
    assert_eq!(checked.products.len(), 3);
    assert!(verification::plan_target_matches(
        &plan,
        checked.products.last().unwrap(),
        &original_target
    ));
    let evidence = checked
        .physical_processing_reports
        .iter()
        .find_map(|step| match &step.evidence {
            verification::PhysicalProcessingEvidence::HomologyAssembly(report) => Some(report),
            _ => None,
        })
        .unwrap();
    assert!(evidence.incomplete);
    assert!(!evidence.terminal_search_complete);
    assert!(evidence.product.as_ref().unwrap().nicks.is_empty());
    let saved = run_lineage::record_run(
        &w,
        &plan,
        pin(&w, 2),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: pin(&w, 0),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: pin(&w, 1),
            },
        ],
        true,
    )
    .unwrap();
    let product = saved.receipt.products.last().unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run_id.clone(),
        step_index: product.step_index,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let w = restart(&saved.workspace);
    let origin = MaterialOrigin::Run(run_material::capsule(&w, &selection).unwrap());
    let source = source(&w, &origin);
    assert!(origin.replay(&source).unwrap().incomplete());
    let target = dna(
        "Independent capped-origin replacement",
        REPLACED,
        Topology::Circular,
    );
    let design = prepare(&w, &origin, &target, &options());
    assert_eq!(
        design.verification.products.len(),
        4,
        "{:?}",
        design.verification.findings
    );
    assert_eq!(design.verification.status, Status::Incomplete);
    assert!(verification::plan_target_matches(
        &design.plan,
        design.verification.products.last().unwrap(),
        &target
    ));
    assert!(
        design
            .verification
            .findings
            .iter()
            .any(|s| s.contains("incomplete") || s.contains("Incomplete"))
    );
}
#[test]
fn original_reference_guards_bounds_and_context_changes_stay_explicit() {
    let (w, origin) = fixture(false);
    let source = source(&w, &origin);
    let target = dna("Independent replacement", REPLACED, Topology::Circular);
    assert!(mutagenesis::compile_inverse_pcr(&source, &target, None, &options()).is_err());
    let plain = dna("Declared original circle", ORIGINAL, Topology::Circular);
    let compiled = mutagenesis::compile_inverse_pcr(&plain, &target, None, &options()).unwrap();
    assert_eq!(compiled.plan.steps.len(), 3);
    compiled
        .validate_binding(&compiled.context_revision)
        .unwrap();
    let lowercase = InversePcrOptions {
        inserted: "gga".into(),
        ..options()
    };
    let normalized = mutagenesis::compile_inverse_pcr(&plain, &target, None, &lowercase).unwrap();
    normalized
        .validate_binding(&normalized.context_revision)
        .unwrap();
    assert_eq!(normalized.options.inserted, "gga");
    assert_eq!(
        normalized.forward.ordered_sequence(),
        compiled.forward.ordered_sequence()
    );
    assert_ne!(normalized.context_revision, compiled.context_revision);
    let huge = InversePcrOptions {
        inserted: "A".repeat(1025),
        ..options()
    };
    assert!(matches!(
        mutagenesis::inverse_pcr_revision(&plain, &target, None, &huge),
        Err(collagene_core::Error::TooLarge)
    ));
    let changed = InversePcrOptions {
        reverse_tail_bases: 2,
        ..options()
    };
    assert_ne!(
        mutagenesis::inverse_pcr_revision(&plain, &target, None, &changed).unwrap(),
        compiled.context_revision
    );
    let linear = dna("Actual linear reference", ORIGINAL, Topology::Linear);
    assert!(mutagenesis::compile_inverse_pcr(&linear, &target, None, &options()).is_err());
}
