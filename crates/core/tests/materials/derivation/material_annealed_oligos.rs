//! Independently literal supplied-oligo products, complete source ancestry and restart.
use collagene_core::{
    chemistry::{AnnealOptions, EndTreatment, FivePrime, Polarity, Terminus},
    documents::workspace::{self, Workspace},
    material::{
        derivation::{self as native, BoundSource, Derivation, SourceOrigin, SourceRole},
        lineage::EventBinding,
        run::MaterialOrigin,
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> Workspace {
    let mut top = Molecule::dna("Supplied top", "AATTCGTACT", Topology::Linear).unwrap();
    top.metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    top.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("Top overhang".into()))],
        needs_review: false,
    });
    Workspace::from_documents(
        "Two independently supplied oligos",
        vec![
            top,
            Molecule::dna("Supplied bottom", "AGTACG", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap()
}
fn options() -> AnnealOptions {
    AnnealOptions {
        min_paired_bases: 6,
        top_chemistry: Terminus::PCR,
        bottom_chemistry: Terminus::PCR,
        ..Default::default()
    }
}
fn request(
    w: &Workspace,
    top: EventBinding,
    bottom: EventBinding,
    options: AnnealOptions,
) -> native::Request {
    let preview = native::preview_annealed_oligos(w, &top, &bottom, &options).unwrap();
    let selection = preview
        .selection(&preview.preview.annealing.candidates[0].id)
        .unwrap();
    native::Request {
        sources: vec![
            BoundSource {
                role: SourceRole::TopOligo,
                event: top,
                origin: SourceOrigin::OriginalReference,
            },
            BoundSource {
                role: SourceRole::BottomOligo,
                event: bottom,
                origin: SourceOrigin::OriginalReference,
            },
        ],
        derivation: Derivation::AnnealedOligos { options, selection },
        name: "Reusable sticky adapter".into(),
        save_primers: false,
    }
}
#[test]
fn literal_sticky_adapter_retains_both_source_roles_and_replays_before_terminal_processing() {
    let w = fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let top = event(&w, 0);
    let bottom = event(&w, 1);
    let request = request(&w, top.clone(), bottom.clone(), options());
    let installed = native::install(&w, &request).unwrap();
    let reopened = workspace::load_workspace(&installed.archive).unwrap();
    installed.receipt.validate_archive(&reopened).unwrap();
    assert_eq!(installed.receipt.sources, request.source_bindings());
    assert!(installed.receipt.primer_attachment.is_none());
    let product = installed.receipt.product;
    let fresh = native::reconstruct(&reopened, &product).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), "AATTCGTACT");
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), "AGTACG");
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert!(!fresh.incomplete());
    let summary = fresh.tracked().duplex.summary().unwrap();
    let ends = summary.ends.unwrap();
    assert_eq!(ends.0.polarity, Polarity::FivePrime);
    assert_eq!(ends.0.overhang, "AATT");
    assert!(
        fresh
            .reference()
            .features
            .iter()
            .any(|f| f.label() == "Top overhang" && f.location == Location::range(0, 4))
    );
    let history = reopened
        .history
        .iter()
        .find(|h| h.id == product.event_id)
        .unwrap();
    assert_eq!(history.parents, vec![top.event_id, bottom.event_id]);
    let capsule = native::capsule(&reopened, &product).unwrap();
    assert_eq!(capsule.schema_version, 2);
    assert_eq!(capsule.source_snapshots.len(), 2);
    let next = native::prepare_install(
        &reopened,
        &native::Request {
            sources: vec![BoundSource {
                role: SourceRole::Template,
                event: product,
                origin: SourceOrigin::Saved(MaterialOrigin::Derived(capsule)),
            }],
            derivation: Derivation::TerminalProcessing {
                action: native::TerminalAction::TreatEnds {
                    treatment: EndTreatment::Phosphorylate,
                },
            },
            name: "Phosphorylated adapter".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let fresh = native::reconstruct(&next.workspace, &next.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert!(!fresh.incomplete());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn both_roles_can_use_one_self_complementary_event_without_duplicate_history_parents() {
    let w = Workspace::from_documents(
        "Self complementary oligo",
        vec![Molecule::dna("One supplied oligo", "AGCT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let e = event(&w, 0);
    let o = AnnealOptions {
        min_paired_bases: 4,
        ..options()
    };
    let request = request(&w, e.clone(), e.clone(), o);
    let p = native::prepare_install(&w, &request).unwrap();
    assert_eq!(p.receipt.sources.len(), 2);
    assert_eq!(
        p.workspace
            .history
            .iter()
            .find(|h| h.id == p.receipt.product.event_id)
            .unwrap()
            .parents,
        vec![e.event_id]
    );
    let fresh = native::reconstruct(&p.workspace, &p.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), "AGCT");
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), "AGCT");
    assert!(!fresh.incomplete());
}
#[test]
fn swapped_roles_changed_selection_and_detached_physical_previews_cannot_authorize_material() {
    let w = fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let top = event(&w, 0);
    let bottom = event(&w, 1);
    let request = request(&w, top.clone(), bottom.clone(), options());
    let mut changed = request.clone();
    changed.sources.swap(0, 1);
    assert!(native::preflight(&w, &changed).is_err());
    let mut changed = request.clone();
    changed.sources[1].event = top.clone();
    assert!(native::prepare_install(&w, &changed).is_err());
    let mut changed = request.clone();
    let Derivation::AnnealedOligos { selection, .. } = &mut changed.derivation else {
        unreachable!()
    };
    selection.expected_product_revision = "0".repeat(64);
    assert!(native::prepare_install(&w, &changed).is_err());
    let mut changed = request;
    changed.save_primers = true;
    assert!(native::preflight(&w, &changed).is_err());
    let preview = native::preview_annealed_oligos(&w, &top, &bottom, &options()).unwrap();
    let mut detached = preview.products[0].clone();
    detached.metadata.remove("duplex_model");
    let detached = Workspace::from_documents(
        "Detached documentary product",
        vec![detached, w.revisions[&bottom.revision].clone()],
        0,
    )
    .unwrap();
    assert!(
        native::preview_annealed_oligos(
            &detached,
            &event(&detached, 0),
            &event(&detached, 1),
            &options()
        )
        .is_err()
    );
    let p = native::prepare_install(&w, &request_for(&w)).unwrap();
    let mut capsule = native::capsule(&p.workspace, &p.receipt.product).unwrap();
    capsule.source_snapshots[1].molecule.name = "Equal letters with altered source".into();
    assert!(capsule.validate_documentary().is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
fn request_for(w: &Workspace) -> native::Request {
    request(w, event(w, 0), event(w, 1), options())
}
#[test]
fn selecting_a_repeated_or_capped_geometry_retains_alternatives_and_incompleteness_after_restart() {
    let w = Workspace::from_documents(
        "Repeated oligos",
        vec![
            Molecule::dna("Top", "AAAA", Topology::Linear).unwrap(),
            Molecule::dna("Bottom", "TTTT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    for max_candidates in [1, 20] {
        let o = AnnealOptions {
            min_paired_bases: 2,
            max_candidates,
            ..options()
        };
        let preview =
            native::preview_annealed_oligos(&w, &event(&w, 0), &event(&w, 1), &o).unwrap();
        assert!(preview.incomplete);
        assert!(!preview.preview.annealing.unique);
        if max_candidates > 1 {
            assert!(preview.products.len() > 1)
        } else {
            assert!(!preview.preview.annealing.complete)
        }
        let installed = native::install(&w, &request(&w, event(&w, 0), event(&w, 1), o)).unwrap();
        let reopened = workspace::load_workspace(&installed.archive).unwrap();
        assert!(installed.receipt.incomplete);
        assert!(
            native::reconstruct(&reopened, &installed.receipt.product)
                .unwrap()
                .incomplete()
        );
    }
}

#[test]
fn saved_annealed_adapters_supply_literal_sticky_and_blunt_junctions_to_later_ligation() {
    use collagene_core::{
        chemistry::ThreePrime,
        planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
        verification::{self, evidence::Status},
    };
    let mut w = Workspace::from_documents(
        "Four supplied strands before cloning",
        vec![
            Molecule::dna("Left top", "AATTCGTACT", Topology::Linear).unwrap(),
            Molecule::dna("Left bottom", "AGTACG", Topology::Linear).unwrap(),
            Molecule::dna("Right top", "AGTC", Topology::Linear).unwrap(),
            Molecule::dna("Right bottom", "AATTGACT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let mut o = options();
    o.top_chemistry = Terminus {
        five_prime: FivePrime::Phosphate,
        three_prime: ThreePrime::Hydroxyl,
    };
    o.bottom_chemistry = o.top_chemistry;
    let first =
        native::prepare_install(&w, &request(&w, event(&w, 0), event(&w, 1), o.clone())).unwrap();
    w = first.workspace;
    o.min_paired_bases = 4;
    let second = native::prepare_install(&w, &request(&w, event(&w, 2), event(&w, 3), o)).unwrap();
    w = second.workspace;
    let w = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    let products = [first.receipt.product, second.receipt.product];
    let inputs = products
        .iter()
        .enumerate()
        .map(|(i, e)| {
            let m = w.revisions[&e.revision].clone();
            PlanInput {
                id: format!("adapter{i}"),
                revision: m.revision(),
                molecule: m,
            }
        })
        .collect();
    let mut steps = products
        .iter()
        .enumerate()
        .map(|(i, e)| Operation::ReconstructMaterial {
            template: Source::Input(format!("adapter{i}")),
            origin: native::capsule(&w, e).unwrap(),
        })
        .collect::<Vec<_>>();
    steps.push(collagene_core::cloning::ligation::outcomes::operation(
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
    let target = Molecule::dna(
        "Independent circular target",
        "AATTCGTACTAGTC",
        Topology::Circular,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Source-backed annealed adapters into circle".into(),
        target_revision: target.revision(),
        inputs,
        steps,
    };
    let report = verification::engine::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let html = verification::output::export_protocol_html(&plan, &target);
    assert!(html.to_lowercase().contains("ligation"));
    let wrong = Molecule::dna(
        "Independently wrong target",
        "AATTCGTACTAGTA",
        Topology::Circular,
    )
    .unwrap();
    let mut wrong_plan = plan;
    wrong_plan.target_revision = wrong.revision();
    assert_eq!(
        verification::engine::verify(&wrong_plan, &wrong).status,
        Status::Invalid
    );
}
