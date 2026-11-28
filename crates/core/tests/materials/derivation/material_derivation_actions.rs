//! Literal intentional substitutions and real saved-strand repair, without targets
//! for the new derivations. Saved runs only establish their independently declared originals.
use collagene_core::{
    chemistry::{
        EndRepairMethod, EndTreatment, FivePrime, Strand, Terminus, ThreePrime,
        methylation::MethylationState,
        tracked_duplex::{self as modified_duplex, CallState, Chemical, GenerationAssumptions},
    },
    cloning::pcr::{material as material_pcr, mutagenesis::PcrOptions},
    documents::workspace::{self, Workspace},
    material::{
        derivation as d,
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run::{self as run_material, MaterialOrigin},
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
};
use collagene_model::{Feature, Location, Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
const MUTANT: &str = "AGACCGATGCTACGTAGGACTTCAGACCTGGA";
const TAILED: &str = "CCAGACCGATGCTACGTAGGACTTCAGACCTGGATC";
fn pin(w: &Workspace, index: usize) -> EventBinding {
    let doc = &w.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn options() -> PcrOptions {
    let mut o = PcrOptions::default();
    o.binding.max_mismatches = 1;
    o.binding.min_three_prime_exact = 2;
    o
}
fn oligos() -> (Primer, Primer) {
    (
        Primer::new("Actual mutant F", "AGACCGAT", "CC").unwrap(),
        Primer::new("Actual mutant R", "TCCAGGTC", "GA").unwrap(),
    )
}
fn mismatch_request(
    w: &Workspace,
    source: EventBinding,
    origin: d::SourceOrigin,
    f: Primer,
    r: Primer,
    o: PcrOptions,
    name: &str,
) -> d::Request {
    let m = &w.revisions[&source.revision];
    let preview = d::pcr::preview_mismatch_pcr(m, &origin, &f, &r, &o).unwrap();
    assert!(preview.reusable(), "{:?}", preview.preview.report.warnings);
    d::Request {
        sources: vec![collagene_core::material::derivation::BoundSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: source,
            origin: origin,
        }],
        derivation: d::Derivation::MismatchPcr {
            forward: f,
            reverse: r,
            options: o,
            selection: d::MismatchPcrSelection {
                selected_index: 0,
                expected_preview_revision: preview.preview.documents[0].revision(),
                geometry: preview.preview.report.products[0].clone(),
            },
        },
        name: name.into(),
        save_primers: true,
    }
}
#[test]
fn both_encoded_substitutions_and_full_tails_survive_archive_replay_and_downstream_pcr() {
    let mut source = Molecule::dna("Original annotated stock", ORIGINAL, Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 8),
        qualifiers: vec![
            ("translation".into(), Some("OLD".into())),
            ("codon_start".into(), Some("1".into())),
        ],
        needs_review: false,
    });
    let w = Workspace::from_documents("No intended target for derivation", vec![source.clone()], 0)
        .unwrap();
    let (f, r) = oligos();
    let preview = d::pcr::preview_mismatch_pcr(
        &source,
        &d::SourceOrigin::OriginalReference,
        &f,
        &r,
        &options(),
    )
    .unwrap();
    assert_eq!(preview.preview.report.products.len(), 1);
    assert!(preview.preview.report.enumeration_complete);
    assert!(preview.reusable());
    assert_eq!(preview.preview.documents[0].sequence.as_str(), TAILED);
    let mut detached_preview = preview.preview.documents[0].clone();
    detached_preview.metadata.remove("duplex_model");
    assert!(
        d::pcr::preview_pcr(
            &detached_preview,
            &d::SourceOrigin::OriginalReference,
            &f,
            &r,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    assert_eq!(
        preview.preview.report.products[0]
            .substitutions
            .iter()
            .map(|c| c.source_position)
            .collect::<Vec<_>>(),
        vec![2, 30]
    );
    let request = mismatch_request(
        &w,
        pin(&w, 0),
        d::SourceOrigin::OriginalReference,
        f,
        r,
        options(),
        "Explicit primer-encoded product",
    );
    let before = workspace::save_workspace(&w).unwrap();
    let installed = d::prepare_install(&w, &request).unwrap();
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    assert_eq!(installed.workspace.documents.len(), 2);
    assert!(
        run_lineage::list_runs(&installed.workspace)
            .unwrap()
            .is_empty()
    );
    installed
        .receipt
        .validate_archive(&installed.workspace)
        .unwrap();
    let saved =
        workspace::load_workspace(&workspace::save_workspace(&installed.workspace).unwrap())
            .unwrap();
    let fresh = d::reconstruct(&saved, &installed.receipt.product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), TAILED);
    assert!(!fresh.incomplete());
    assert_eq!(fresh.tracked().duplex.left, Terminus::PCR);
    assert_eq!(fresh.tracked().duplex.right, Terminus::PCR);
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .known_modified,
        0
    );
    let feature = &fresh.reference().features[0];
    assert!(feature.needs_review);
    assert_eq!(feature.location.spans(), vec![(2, 10, false)]);
    assert!(
        !feature
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation" || key == "codon_start")
    );
    let record = d::load_record(&saved.attachments[&installed.receipt.record_hash].bytes).unwrap();
    let table = d::primer_table(&record).unwrap().unwrap();
    assert_eq!(table.records.len(), 2);
    assert_eq!(
        table.records[0].primer.as_ref().unwrap().ordered_sequence(),
        "CCAGACCGAT"
    );
    assert_eq!(
        table.records[1].primer.as_ref().unwrap().ordered_sequence(),
        "GATCCAGGTC"
    );
    let cap = d::capsule(&saved, &installed.receipt.product).unwrap();
    let origin = MaterialOrigin::Derived(cap);
    let f = Primer::new("New inner F", "AGACCGAT", "").unwrap();
    let r = Primer::new("New inner R", "TCCAGGTC", "").unwrap();
    let products = material_pcr::simulate(
        fresh.reference(),
        &origin,
        &f,
        &r,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(products.len(), 1);
    assert_eq!(products[0].molecule.sequence.as_str(), MUTANT);
    let child = d::prepare_install(
        &saved,
        &d::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: installed.receipt.product,
                origin: d::SourceOrigin::Saved(origin),
            }],
            derivation: d::Derivation::ExactPcr {
                forward: f,
                reverse: r,
                selection: d::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: products[0].molecule.revision(),
                    geometry: products[0].geometry.clone(),
                },
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            name: "New PCR from saved mutant".into(),
            save_primers: true,
        },
    )
    .unwrap();
    let replay = d::reconstruct(&child.workspace, &child.receipt.product).unwrap();
    assert_eq!(replay.reference().sequence.as_str(), MUTANT);
    let record =
        d::load_record(&child.workspace.attachments[&child.receipt.record_hash].bytes).unwrap();
    let table = d::primer_table(&record).unwrap().unwrap();
    assert_eq!(table.records.len(), 2);
    assert_eq!(table.records[0].primer.as_ref().unwrap().tail, "");
}
fn saved_staggered() -> (Workspace, EventBinding, MaterialOrigin) {
    let mut stock = Molecule::dna(
        "Declared Dam original",
        "AGTCCGATCGTACAGATATC",
        Topology::Linear,
    )
    .unwrap();
    stock.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(4, 10),
        qualifiers: vec![],
        needs_review: false,
    });
    let target = Molecule::dna(
        "Independent shortened top reference",
        "AGTCCGATCGTAC",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Original untrimmed strands",
        vec![stock.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Declared Dam generation, actual EcoRV cleavage and top shortening".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", stock)],
        steps: vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("stock".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::FullyModified,
                    dcm: MethylationState::Unmodified,
                    cpg: MethylationState::Unmodified,
                },
            },
            // EcoRV's real GAT|ATC cleavage establishes the accessible3′OH.
            // The original outer left5′ group stays unknown; the cut right5′ is P.
            Operation::Digest {
                template: Source::Step(0),
                enzymes: vec!["EcoRV".into()],
                fragment_index: 0,
            },
            Operation::ResectThreePrime {
                source: Source::Step(1),
                cuts: collagene_core::cloning::overlap::resection::Cuts { left: 0, right: 4 },
            },
        ],
    };
    let reconstruction = collagene_core::verification::engine::verify(&plan, &target);
    assert_eq!(
        reconstruction.products.len(),
        plan.steps.len(),
        "actual EcoRV cleavage and resection must reconstruct every step: {:?}",
        reconstruction.findings
    );
    assert_eq!(
        reconstruction.products.last().unwrap().sequence.as_str(),
        "AGTCCGATCGTAC"
    );
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
    let p = &saved.receipt.products[2];
    let event = EventBinding {
        document_id: p.document_id.clone(),
        event_id: p.event_id.clone(),
        revision: p.revision.clone(),
    };
    let origin = collagene_core::material::capsule::capsule(
        &saved.workspace,
        &run_material::Selection {
            run_id: saved.receipt.run_id,
            step_index: 2,
            event: event.clone(),
        },
    )
    .unwrap();
    (saved.workspace, event, MaterialOrigin::Run(origin))
}
#[test]
fn direct_saved_strand_fill_preserves_old_modifications_and_groups_then_enables_declared_mismatch()
{
    let (w, event, origin) = saved_staggered();
    let source = &w.revisions[&event.revision];
    let old = origin.replay(source).unwrap();
    assert_eq!(
        (
            old.tracked().duplex.top.len(),
            old.tracked().duplex.bottom.len()
        ),
        (13, 17)
    );
    let f = Primer::new("New F", "AATCC", "").unwrap();
    let r = Primer::new("New R", "ATCTG", "").unwrap();
    assert!(
        d::pcr::preview_mismatch_pcr(
            source,
            &d::SourceOrigin::Saved(origin.clone()),
            &f,
            &r,
            &options()
        )
        .is_err()
    );
    let preview = d::end_repair::preview_repair(
        source,
        &d::SourceOrigin::Saved(origin.clone()),
        EndRepairMethod::FillFivePrimeOverhangs,
    )
    .unwrap();
    assert_eq!(preview.repair.top_added_3prime, "AGAT");
    let mut detached_preview = preview.molecule.clone();
    detached_preview.metadata.remove("duplex_model");
    assert!(
        d::pcr::preview_pcr(
            &detached_preview,
            &d::SourceOrigin::OriginalReference,
            &f,
            &r,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    assert_eq!(preview.molecule.sequence.as_str(), "AGTCCGATCGTACAGAT");
    let installed = d::prepare_install(
        &w,
        &d::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: event,
                origin: d::SourceOrigin::Saved(origin),
            }],
            derivation: d::Derivation::EndRepair {
                method: EndRepairMethod::FillFivePrimeOverhangs,
            },
            name: "Explicit fill from actual bottom template".into(),
            save_primers: false,
        },
    )
    .unwrap();
    assert!(installed.receipt.primer_attachment.is_none());
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&installed.workspace).unwrap())
            .unwrap();
    let fresh = d::reconstruct(&reopened, &installed.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.bottom_start, 0);
    assert_eq!(fresh.tracked().duplex.top.len(), 17);
    assert_eq!(fresh.tracked().duplex.bottom.len(), 17);
    assert_eq!(
        fresh.tracked().duplex.left.five_prime,
        old.tracked().duplex.left.five_prime
    );
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        old.tracked().duplex.right.five_prime
    );
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Unknown);
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        FivePrime::Phosphate
    );
    assert_eq!(
        fresh
            .tracked()
            .state_at(Strand::Top, Chemical::N6A, 6)
            .unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        fresh
            .tracked()
            .state_at(Strand::Top, Chemical::N6A, 15)
            .unwrap(),
        Some(CallState::Unmodified)
    );
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .known_modified,
        modified_duplex::summary(old.tracked())
            .unwrap()
            .known_modified
    );
    assert_eq!(
        fresh.reference().features[0].location.spans(),
        vec![(4, 10, false)]
    );
    let cap = d::capsule(&reopened, &installed.receipt.product).unwrap();
    let origin = d::SourceOrigin::Saved(MaterialOrigin::Derived(cap));
    let request = mismatch_request(
        &reopened,
        installed.receipt.product,
        origin,
        f,
        r,
        options(),
        "Explicit mutation after requested repair",
    );
    let result = d::prepare_install(&reopened, &request).unwrap();
    let fresh = d::reconstruct(&result.workspace, &result.receipt.product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), "AATCCGATCGTACAGAT");
    assert_eq!(
        modified_duplex::summary(fresh.tracked())
            .unwrap()
            .known_modified,
        0
    );
}
#[test]
fn capped_singleton_search_remains_conditional_and_competing_products_cannot_be_selected() {
    let text = format!("{ORIGINAL}CCCCGACCTGTA");
    let source = Molecule::dna("Two actual reverse sites", &text, Topology::Linear).unwrap();
    let w = Workspace::from_documents("Competing mismatch placements", vec![source.clone()], 0)
        .unwrap();
    let f = Primer::new("F", "AGACCGAT", "").unwrap();
    let r = Primer::new("R", "TCCAGGTC", "").unwrap();
    let full = d::pcr::preview_mismatch_pcr(
        &source,
        &d::SourceOrigin::OriginalReference,
        &f,
        &r,
        &options(),
    )
    .unwrap();
    assert!(full.preview.report.products.len() > 1);
    assert!(!full.reusable());
    let forged = d::Request {
        sources: vec![collagene_core::material::derivation::BoundSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: pin(&w, 0),
            origin: d::SourceOrigin::OriginalReference,
        }],
        derivation: d::Derivation::MismatchPcr {
            forward: f.clone(),
            reverse: r.clone(),
            options: options(),
            selection: d::MismatchPcrSelection {
                selected_index: 0,
                expected_preview_revision: full.preview.documents[0].revision(),
                geometry: full.preview.report.products[0].clone(),
            },
        },
        name: "No arbitrary isolation".into(),
        save_primers: false,
    };
    assert!(d::prepare_install(&w, &forged).is_err());
    let mut capped = options();
    capped.max_products = 1;
    let preview = d::pcr::preview_mismatch_pcr(
        &source,
        &d::SourceOrigin::OriginalReference,
        &f,
        &r,
        &capped,
    )
    .unwrap();
    assert!(preview.reusable());
    assert!(preview.incomplete);
    assert!(!preview.preview.report.enumeration_complete);
    assert_eq!(preview.preview.documents[0].sequence.as_str(), MUTANT);
    let request = mismatch_request(
        &w,
        pin(&w, 0),
        d::SourceOrigin::OriginalReference,
        f,
        r,
        capped,
        "Conditional declared path",
    );
    let saved = d::prepare_install(&w, &request).unwrap();
    assert!(saved.receipt.incomplete);
    let fresh = d::reconstruct(&saved.workspace, &saved.receipt.product).unwrap();
    assert!(fresh.incomplete());
    assert_eq!(fresh.reference().sequence.as_str(), MUTANT);
    let origin = d::SourceOrigin::Saved(MaterialOrigin::Derived(
        d::capsule(&saved.workspace, &saved.receipt.product).unwrap(),
    ));
    let repair = d::prepare_install(
        &saved.workspace,
        &d::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: saved.receipt.product,
                origin: origin,
            }],
            derivation: d::Derivation::EndRepair {
                method: EndRepairMethod::Blunt,
            },
            name: "Keep conditional limitations on blunt ends".into(),
            save_primers: false,
        },
    )
    .unwrap();
    assert!(repair.receipt.incomplete);
}
#[test]
fn mismatch_selection_policy_source_alias_and_new_primer_order_choices_are_pinned() {
    let stock = Molecule::dna("Original stock", ORIGINAL, Topology::Linear).unwrap();
    let w = Workspace::from_documents("Exact immutable event", vec![stock.clone()], 0).unwrap();
    let (f, r) = oligos();
    let request = mismatch_request(
        &w,
        pin(&w, 0),
        d::SourceOrigin::OriginalReference,
        f,
        r,
        options(),
        "Pinned mismatch",
    );
    let before = workspace::save_workspace(&w).unwrap();
    let mut altered = request.clone();
    let d::Derivation::MismatchPcr { selection, .. } = &mut altered.derivation else {
        unreachable!()
    };
    selection.geometry.sequence = collagene_model::Sequence::dna("AAAA").unwrap();
    assert!(d::prepare_install(&w, &altered).is_err());
    let mut altered = request.clone();
    let d::Derivation::MismatchPcr { options, .. } = &mut altered.derivation else {
        unreachable!()
    };
    options.binding.min_three_prime_exact = 9;
    assert!(d::preflight(&w, &altered).is_err());
    let mut changed = w.clone();
    let alias =
        collagene_core::documents::history::DocumentHistory::new("Equal-reference-alias", stock)
            .unwrap();
    alias.contribute(&mut changed).unwrap();
    let saved = d::prepare_install(&w, &request).unwrap();
    let cap = d::capsule(&saved.workspace, &saved.receipt.product).unwrap();
    let mut wrong = cap.clone();
    wrong.source_snapshots[0].molecule.name = "Substituted original".into();
    assert!(wrong.validate_documentary().is_err());
    assert!(
        MaterialOrigin::Derived(cap)
            .validate_for_workspace(&changed, &pin(&changed, 1))
            .is_err()
    );
    let repair = d::Request {
        sources: vec![collagene_core::material::derivation::BoundSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: pin(&w, 0),
            origin: d::SourceOrigin::OriginalReference,
        }],
        derivation: d::Derivation::EndRepair {
            method: EndRepairMethod::Blunt,
        },
        name: "No historical primer ordering".into(),
        save_primers: true,
    };
    assert!(d::preflight(&w, &repair).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn end_repair_preserves_real_internal_nick_and_does_not_enable_unsupported_nicked_mismatch() {
    let a = Molecule::dna("Original A", "AGTCCG", Topology::Linear).unwrap();
    let b = Molecule::dna("Original B", "GACTGA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent linear nicked intent",
        "AGTCCGGACTGA",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Explicit nicked source",
        vec![a.clone(), b.clone(), target.clone()],
        2,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "One eligible strand joins".into(),
        target_revision: target.revision(),
        inputs: vec![input("a", a), input("b", b)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("a".into()),
                forward: Primer::new("AF", "AGT", "").unwrap(),
                reverse: Primer::new("AR", "CGG", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("b".into()),
                forward: Primer::new("BF", "GAC", "").unwrap(),
                reverse: Primer::new("BR", "TCA", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
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
    let p = &saved.receipt.products[3];
    let event = EventBinding {
        document_id: p.document_id.clone(),
        event_id: p.event_id.clone(),
        revision: p.revision.clone(),
    };
    let origin = collagene_core::material::capsule::capsule(
        &saved.workspace,
        &run_material::Selection {
            run_id: saved.receipt.run_id,
            step_index: 3,
            event: event.clone(),
        },
    )
    .unwrap();
    let source = &saved.workspace.revisions[&event.revision];
    let preview = d::end_repair::preview_repair(
        source,
        &d::SourceOrigin::Saved(MaterialOrigin::Run(origin.clone())),
        EndRepairMethod::Blunt,
    )
    .unwrap();
    assert_eq!(preview.repair.source.nicks.len(), 1);
    assert_eq!(preview.repair.duplex.nicks, preview.repair.source.nicks);
    let repaired = d::prepare_install(
        &saved.workspace,
        &d::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: event,
                origin: d::SourceOrigin::Saved(MaterialOrigin::Run(origin)),
            }],
            derivation: d::Derivation::EndRepair {
                method: EndRepairMethod::Blunt,
            },
            name: "Explicitly retained nick".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let fresh = d::reconstruct(&repaired.workspace, &repaired.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.nicks.len(), 1);
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        FivePrime::Phosphate
    );
    assert_eq!(
        fresh.tracked().duplex.left.three_prime,
        ThreePrime::Hydroxyl
    );
    let origin = d::SourceOrigin::Saved(MaterialOrigin::Derived(
        d::capsule(&repaired.workspace, &repaired.receipt.product).unwrap(),
    ));
    assert!(
        d::pcr::preview_mismatch_pcr(
            fresh.reference(),
            &origin,
            &Primer::new("F", "AGTCCG", "").unwrap(),
            &Primer::new("R", "TCAGTC", "").unwrap(),
            &options()
        )
        .is_err()
    );
}
