//! Actual standalone PCR/digest descendants retain source chemistry and intended ligation paths.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    chemistry::{EndTreatment, FivePrime, Strand, ThreePrime},
    material_derivation as native,
    methylation::MethylationState,
    modified_duplex::{self, CallState, Chemical, GenerationAssumptions},
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    tailing::TailingBase,
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
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
fn derived(w: &Workspace, event: &EventBinding) -> native::SourceOrigin {
    native::SourceOrigin::Saved(MaterialOrigin::Derived(native::capsule(w, event).unwrap()))
}
fn request(
    source: EventBinding,
    origin: native::SourceOrigin,
    action: native::TerminalAction,
) -> native::Request {
    native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: source,
            origin: origin,
        }],
        derivation: native::Derivation::TerminalProcessing { action },
        name: "Explicit terminal product".into(),
        save_primers: false,
    }
}
fn pcr(w: &Workspace, index: usize, f: &str, r: &str) -> native::PreparedInstallation {
    let source = pin(w, index);
    let molecule = &w.revisions[&source.revision];
    let forward = Primer::new("Historical PCR F", f, "").unwrap();
    let reverse = Primer::new("Historical PCR R", r, "").unwrap();
    let preview = native::preview_pcr(
        molecule,
        &native::SourceOrigin::OriginalReference,
        &forward,
        &reverse,
        Strand::Top,
    )
    .unwrap();
    assert_eq!(preview.len(), 1);
    native::prepare_install(
        w,
        &native::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: source,
                origin: native::SourceOrigin::OriginalReference,
            }],
            derivation: native::Derivation::ExactPcr {
                forward,
                reverse,
                template_strand: Strand::Top,
                selection: native::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: preview[0].molecule.revision(),
                    geometry: preview[0].geometry.clone(),
                },
            },
            name: "Known actual PCR material".into(),
            save_primers: true,
        },
    )
    .unwrap()
}
fn stock() -> Workspace {
    let mut source = Molecule::dna("Source six", "AAACCC", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(1, 4),
        qualifiers: vec![("label".into(), Some("Retained source feature".into()))],
        needs_review: false,
    });
    Workspace::from_documents("No target for processing", vec![source], 0).unwrap()
}
#[test]
fn all_four_explicit_single_base_tails_and_five_prime_treatments_keep_literal_actual_strands_after_restart()
 {
    let seed = pcr(&stock(), 0, "AAA", "GGG");
    let source = seed.receipt.product;
    let w = seed.workspace;
    for (base, top, bottom) in [
        (TailingBase::A, "AAACCCA", "GGGTTTA"),
        (TailingBase::T, "AAACCCT", "GGGTTTT"),
        (TailingBase::G, "AAACCCG", "GGGTTTG"),
        (TailingBase::C, "AAACCCC", "GGGTTTC"),
    ] {
        let action = native::TerminalAction::Tail3Prime { base };
        let origin = derived(&w, &source);
        let before = workspace::save_workspace(&w).unwrap();
        let preview = native::preview_terminal_processing(&w, &source, &origin, &action).unwrap();
        preview
            .validate_binding(
                &native::terminal_processing_preview_revision(&w, &source, &origin, &action)
                    .unwrap(),
            )
            .unwrap();
        assert_eq!(preview.product.sequence.as_str(), top);
        assert_eq!(preview.duplex.bottom_start, -1);
        assert_eq!(preview.duplex.top_length, 7);
        assert_eq!(preview.duplex.bottom_length, 7);
        assert!(preview.duplex.nicks.is_empty());
        assert!(!preview.incomplete);
        let processed = native::install(&w, &request(source.clone(), origin, action)).unwrap();
        assert!(processed.receipt.primer_attachment.is_none());
        assert_eq!(workspace::save_workspace(&w).unwrap(), before);
        let reopened = workspace::load_workspace(&processed.archive).unwrap();
        let fresh = native::reconstruct(&reopened, &processed.receipt.product).unwrap();
        assert_eq!(fresh.tracked().duplex.top.as_str(), top);
        assert_eq!(fresh.tracked().duplex.bottom.as_str(), bottom);
        assert_eq!(
            fresh.reference().features[0].location,
            Location::range(1, 4)
        );
        let action = native::TerminalAction::TreatEnds {
            treatment: EndTreatment::Phosphorylate,
        };
        let next = native::prepare_install(
            &reopened,
            &request(
                processed.receipt.product.clone(),
                derived(&reopened, &processed.receipt.product),
                action,
            ),
        )
        .unwrap();
        let fresh = native::reconstruct(&next.workspace, &next.receipt.product).unwrap();
        assert_eq!(fresh.tracked().duplex.top.as_str(), top);
        assert_eq!(fresh.tracked().duplex.bottom.as_str(), bottom);
        assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
        assert_eq!(
            fresh.tracked().duplex.right.five_prime,
            FivePrime::Phosphate
        );
        assert_eq!(
            fresh.tracked().duplex.left.three_prime,
            ThreePrime::Hydroxyl
        );
        let action = native::TerminalAction::TreatEnds {
            treatment: EndTreatment::Dephosphorylate,
        };
        let last = native::install(
            &next.workspace,
            &request(
                next.receipt.product.clone(),
                derived(&next.workspace, &next.receipt.product),
                action,
            ),
        )
        .unwrap();
        let reopened = workspace::load_workspace(&last.archive).unwrap();
        let fresh = native::reconstruct(&reopened, &last.receipt.product).unwrap();
        assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Hydroxyl);
        assert_eq!(fresh.tracked().duplex.right.five_prime, FivePrime::Hydroxyl);
        assert_eq!(fresh.tracked().duplex.top.as_str(), top);
    }
}
#[test]
fn independently_declared_ta_ligation_keeps_intended_path_and_alternative_products() {
    let w = Workspace::from_documents(
        "Standalone material before assembly",
        vec![
            Molecule::dna("Left", "AAACCC", Topology::Linear).unwrap(),
            Molecule::dna("Right", "GGGTTT", Topology::Linear).unwrap(),
            Molecule::dna(
                "Independent intended circle",
                "AAACCCTGGGTTTA",
                Topology::Circular,
            )
            .unwrap(),
        ],
        2,
    )
    .unwrap();
    let left = pcr(&w, 0, "AAA", "GGG");
    let right = pcr(&left.workspace, 1, "GGG", "AAA");
    let mut w = right.workspace;
    let mut events = vec![];
    for (source, base) in [
        (left.receipt.product, TailingBase::T),
        (right.receipt.product, TailingBase::A),
    ] {
        let tailed = native::prepare_install(
            &w,
            &request(
                source.clone(),
                derived(&w, &source),
                native::TerminalAction::Tail3Prime { base },
            ),
        )
        .unwrap();
        w = tailed.workspace;
        let treated = native::prepare_install(
            &w,
            &request(
                tailed.receipt.product.clone(),
                derived(&w, &tailed.receipt.product),
                native::TerminalAction::TreatEnds {
                    treatment: EndTreatment::Phosphorylate,
                },
            ),
        )
        .unwrap();
        w = treated.workspace;
        events.push(treated.receipt.product);
    }
    let w = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    let target = w.revisions[&pin(&w, 2).revision].clone();
    let mut steps = events
        .iter()
        .enumerate()
        .map(|(index, event)| Operation::ReconstructMaterial {
            template: Source::Input(format!("source{index}")),
            origin: native::capsule(&w, event).unwrap(),
        })
        .collect::<Vec<_>>();
    steps.push(collagene_core::ligation_outcomes::operation(
        vec![
            OrientedSource {
                source: Source::Step(0),
                reverse: false,
            },
            OrientedSource {
                source: Source::Step(1),
                reverse: false,
            },
        ],
        true,
        false,
    ));
    let plan = DesignPlan {
        name: "Declared TA orientation after standalone tails".into(),
        target_revision: target.revision(),
        inputs: events
            .iter()
            .enumerate()
            .map(|(i, e)| input(&format!("source{i}"), w.revisions[&e.revision].clone()))
            .collect(),
        steps,
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let outcomes = report
        .physical_processing_reports
        .iter()
        .find_map(|step| {
            if let verification::PhysicalProcessingEvidence::LigationOutcomes(e) = &step.evidence {
                Some(&e.report)
            } else {
                None
            }
        })
        .unwrap();
    assert_eq!(outcomes.intended.sequence.as_str(), "AAACCCTGGGTTTA");
    assert!(outcomes.enumeration_complete);
    assert!(!outcomes.alternatives.is_empty());
    let saved = run_lineage::record_run(
        &w,
        &plan,
        pin(&w, 2),
        events
            .iter()
            .enumerate()
            .map(|(i, e)| InputEventBinding {
                input_id: format!("source{i}"),
                event: e.clone(),
            })
            .collect(),
        false,
    )
    .unwrap();
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    assert!(reopened.attachments.len() > w.attachments.len());
    let html = verification::export_protocol_html(&plan, &target);
    assert!(html.contains("ligation") || html.contains("Ligation"));
}
#[test]
fn unknown_original_chemistry_wrong_origins_and_detached_preview_cannot_supply_tailing_stock() {
    let w = stock();
    let source = pin(&w, 0);
    let treatment = native::TerminalAction::TreatEnds {
        treatment: EndTreatment::Phosphorylate,
    };
    let tail = native::TerminalAction::Tail3Prime {
        base: TailingBase::A,
    };
    let report = native::preview_terminal_processing(
        &w,
        &source,
        &native::SourceOrigin::OriginalReference,
        &treatment,
    )
    .unwrap();
    assert_eq!(
        report.duplex.ends.as_ref().unwrap().0.chemistry.five_prime,
        FivePrime::Phosphate
    );
    assert_eq!(
        report.duplex.ends.as_ref().unwrap().0.chemistry.three_prime,
        ThreePrime::Unknown
    );
    assert!(
        native::preview_terminal_processing(
            &w,
            &source,
            &native::SourceOrigin::OriginalReference,
            &tail
        )
        .is_err()
    );
    let processed = native::prepare_install(
        &w,
        &request(
            source.clone(),
            native::SourceOrigin::OriginalReference,
            treatment,
        ),
    )
    .unwrap();
    assert!(
        native::preview_terminal_processing(
            &processed.workspace,
            &processed.receipt.product,
            &derived(&processed.workspace, &processed.receipt.product),
            &tail
        )
        .is_err()
    );
    let mut detached = report.product.clone();
    detached.metadata.remove("duplex_model");
    let detached = Workspace::from_documents("Detached preview", vec![detached], 0).unwrap();
    assert!(
        native::preview_terminal_processing(
            &detached,
            &pin(&detached, 0),
            &native::SourceOrigin::OriginalReference,
            &treatment
        )
        .is_err()
    );
    let mut altered = report.clone();
    altered.source_sequence_blake3 = "f".repeat(64);
    assert!(altered.validate_binding(&report.context_revision).is_err());
    let mut altered = report.clone();
    altered.action = tail;
    assert!(altered.validate_binding(&report.context_revision).is_err());
    let mut altered = report.clone();
    altered.modifications.known_modified += 1;
    assert!(altered.validate_binding(&report.context_revision).is_err());
    let mut req = request(source, native::SourceOrigin::OriginalReference, treatment);
    req.save_primers = true;
    assert!(native::install(&w, &req).is_err());
    let seed = pcr(&w, 0, "AAA", "GGG");
    let source = seed.receipt.product;
    let origin = derived(&seed.workspace, &source);
    assert!(native::preview_terminal_processing(&w, &source, &origin, &treatment).is_err());
    let mut alias = source.clone();
    alias.event_id = pin(&seed.workspace, 0).event_id;
    assert!(
        native::preview_terminal_processing(&seed.workspace, &alias, &origin, &treatment).is_err()
    );
    assert!(
        native::preview_terminal_processing(
            &seed.workspace,
            &source,
            &native::SourceOrigin::OriginalReference,
            &tail
        )
        .is_err()
    );
    let circle = Workspace::from_documents(
        "Closed source",
        vec![Molecule::dna("Circle", "AAACCC", Topology::Circular).unwrap()],
        0,
    )
    .unwrap();
    assert!(
        native::preview_terminal_processing(
            &circle,
            &pin(&circle, 0),
            &native::SourceOrigin::OriginalReference,
            &treatment
        )
        .is_err()
    );
    assert!(serde_json::from_str::<native::TerminalAction>("\"None\"").is_err());
}

#[test]
fn cut_created_oh_and_original_dam_calls_survive_tailing_and_phosphorylation_without_regeneration()
{
    let source = Molecule::dna(
        "Declared original Dam stock",
        "TTTGATATCAGATCCGATATCAAA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent real middle fragment",
        "ATCAGATCCGAT",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Native cleavage ancestry",
        vec![source.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Declared generation and actual two-sided EcoRV cleavage".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source)],
        steps: vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("stock".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::FullyModified,
                    dcm: MethylationState::Unmodified,
                    cpg: MethylationState::Unmodified,
                },
            },
            Operation::Digest {
                template: Source::Step(0),
                enzymes: vec!["EcoRV".into()],
                fragment_index: 1,
            },
        ],
    };
    let saved = run_lineage::record_run(
        &w,
        &plan,
        pin(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: pin(&w, 0),
        }],
        true,
    )
    .unwrap();
    let product = &saved.receipt.products[1];
    let source = EventBinding {
        document_id: product.document_id.clone(),
        event_id: product.event_id.clone(),
        revision: product.revision.clone(),
    };
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &saved.workspace,
            &run_material::Selection {
                run_id: saved.receipt.run_id,
                step_index: 1,
                event: source.clone(),
            },
        )
        .unwrap(),
    );
    let old = origin
        .replay(&saved.workspace.revisions[&source.revision])
        .unwrap();
    assert_eq!(old.tracked().duplex.top.as_str(), "ATCAGATCCGAT");
    assert_eq!(
        old.tracked()
            .state_at(Strand::Top, Chemical::N6A, 5)
            .unwrap(),
        Some(CallState::Modified)
    );
    let preview = native::preview_terminal_processing(
        &saved.workspace,
        &source,
        &native::SourceOrigin::Saved(origin.clone()),
        &native::TerminalAction::Tail3Prime {
            base: TailingBase::A,
        },
    )
    .unwrap();
    assert_eq!(preview.product.sequence.as_str(), "ATCAGATCCGATA");
    assert_eq!(
        preview.source_modifications.known_modified,
        preview.modifications.known_modified
    );
    assert_eq!(
        preview.modifications.generation,
        preview.source_modifications.generation
    );
    let tail = native::prepare_install(
        &saved.workspace,
        &request(
            source,
            native::SourceOrigin::Saved(origin),
            native::TerminalAction::Tail3Prime {
                base: TailingBase::A,
            },
        ),
    )
    .unwrap();
    let treated = native::install(
        &tail.workspace,
        &request(
            tail.receipt.product.clone(),
            derived(&tail.workspace, &tail.receipt.product),
            native::TerminalAction::TreatEnds {
                treatment: EndTreatment::Dephosphorylate,
            },
        ),
    )
    .unwrap();
    let reopened = workspace::load_workspace(&treated.archive).unwrap();
    let fresh = native::reconstruct(&reopened, &treated.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), "ATCAGATCCGATA");
    assert_eq!(fresh.tracked().duplex.bottom_start, -1);
    assert_eq!(
        fresh
            .tracked()
            .state_at(Strand::Top, Chemical::N6A, 5)
            .unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        fresh
            .tracked()
            .state_at(Strand::Top, Chemical::N6A, 12)
            .unwrap(),
        Some(CallState::Unmodified)
    );
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .known_modified,
        modified_duplex::summary(old.tracked())
            .unwrap()
            .known_modified
    );
}
#[test]
fn tailing_keeps_a_real_internal_nick_and_treatment_refuses_to_seal_it() {
    let a = Molecule::dna("A", "AGTCCG", Topology::Linear).unwrap();
    let b = Molecule::dna("B", "GACTGA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent nicked intent",
        "AGTCCGGACTGA",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Actual one-strand join",
        vec![a.clone(), b.clone(), target.clone()],
        2,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Asymmetric phosphorylation leaves one real nick".into(),
        target_revision: target.revision(),
        inputs: vec![input("a", a), input("b", b)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("a".into()),
                forward: Primer::new("AF", "AGT", "").unwrap(),
                reverse: Primer::new("AR", "CGG", "").unwrap(),
                template_strand: Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("b".into()),
                forward: Primer::new("BF", "GAC", "").unwrap(),
                reverse: Primer::new("BR", "TCA", "").unwrap(),
                template_strand: Strand::Top,
            },
            Operation::TreatEnds {
                template: Source::Step(1),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(2),
                        reverse: false,
                    },
                ],
                circular: false,
                allow_nicked: true,
                options: Default::default(),
            },
        ],
    };
    let saved = run_lineage::record_run(
        &w,
        &plan,
        pin(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: pin(&w, 0),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: pin(&w, 1),
            },
        ],
        true,
    )
    .unwrap();
    let product = &saved.receipt.products[3];
    let source = EventBinding {
        document_id: product.document_id.clone(),
        event_id: product.event_id.clone(),
        revision: product.revision.clone(),
    };
    let origin = native::SourceOrigin::Saved(MaterialOrigin::Run(
        run_material::capsule(
            &saved.workspace,
            &run_material::Selection {
                run_id: saved.receipt.run_id,
                step_index: 3,
                event: source.clone(),
            },
        )
        .unwrap(),
    ));
    let action = native::TerminalAction::Tail3Prime {
        base: TailingBase::A,
    };
    let preview =
        native::preview_terminal_processing(&saved.workspace, &source, &origin, &action).unwrap();
    assert_eq!(preview.source_duplex.nicks.len(), 1);
    assert_eq!(preview.duplex.nicks, preview.source_duplex.nicks);
    assert_eq!(preview.product.sequence.as_str(), "AGTCCGGACTGAA");
    let installed = native::install(&saved.workspace, &request(source, origin, action)).unwrap();
    let reopened = workspace::load_workspace(&installed.archive).unwrap();
    let fresh = native::reconstruct(&reopened, &installed.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.nicks, preview.duplex.nicks);
    assert!(
        native::preview_terminal_processing(
            &reopened,
            &installed.receipt.product,
            &derived(&reopened, &installed.receipt.product),
            &native::TerminalAction::TreatEnds {
                treatment: EndTreatment::Phosphorylate
            }
        )
        .is_err()
    );
}

#[test]
fn conditionally_capped_upstream_pcr_does_not_become_complete_after_terminal_processing() {
    let source = Molecule::dna(
        "Repeated original",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAAGTCCGATGCTACGTAGGACTTCAGACCTGTA",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents("Conditional source", vec![source.clone()], 0).unwrap();
    let forward = Primer::new("Changed F", "AGACCGAT", "").unwrap();
    let reverse = Primer::new("Changed R", "TCCAGGTC", "").unwrap();
    let mut options = collagene_core::mutagenesis::PcrOptions::default();
    options.binding.max_mismatches = 1;
    options.binding.min_three_prime_exact = 2;
    options.max_products = 1;
    let preview = native::preview_mismatch_pcr(
        &source,
        &native::SourceOrigin::OriginalReference,
        &forward,
        &reverse,
        &options,
    )
    .unwrap();
    assert!(preview.reusable());
    assert!(preview.incomplete);
    let seed = native::prepare_install(
        &w,
        &native::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: pin(&w, 0),
                origin: native::SourceOrigin::OriginalReference,
            }],
            derivation: native::Derivation::MismatchPcr {
                forward,
                reverse,
                options,
                selection: native::MismatchPcrSelection {
                    selected_index: 0,
                    expected_preview_revision: preview.preview.documents[0].revision(),
                    geometry: preview.preview.report.products[0].clone(),
                },
            },
            name: "Conditional singleton PCR".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let action = native::TerminalAction::Tail3Prime {
        base: TailingBase::A,
    };
    let source = seed.receipt.product;
    let origin = derived(&seed.workspace, &source);
    let preview =
        native::preview_terminal_processing(&seed.workspace, &source, &origin, &action).unwrap();
    assert!(preview.origin_incomplete);
    assert!(preview.incomplete);
    let mut forged = preview.clone();
    forged.incomplete = false;
    assert!(forged.validate_binding(&preview.context_revision).is_err());
    let installed = native::install(&seed.workspace, &request(source, origin, action)).unwrap();
    let reopened = workspace::load_workspace(&installed.archive).unwrap();
    assert!(installed.receipt.incomplete);
    assert!(
        native::reconstruct(&reopened, &installed.receipt.product)
            .unwrap()
            .incomplete()
    );
}
