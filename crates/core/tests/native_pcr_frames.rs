//! Independent actual13/17 strand oracles. No bottom Molecule replaces the saved stock.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    chemistry::{AnnealOptions, Strand, Terminus},
    cloning_outputs, material_derivation as d, material_pcr, modified_duplex,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
const TOP: &str = "AGTCCGATGCTAC";
const BOTTOM: &str = "GTAGGTAGCATCGGACT";
fn pin(w: &Workspace, i: usize) -> EventBinding {
    let doc = &w.documents[i];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn saved() -> (Workspace, EventBinding, MaterialOrigin) {
    saved_pair(TOP, BOTTOM)
}
fn saved_pair(top_dna: &str, bottom_dna: &str) -> (Workspace, EventBinding, MaterialOrigin) {
    let mut top = Molecule::dna("Independently ordered top", top_dna, Topology::Linear).unwrap();
    top.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, top_dna.len()),
        qualifiers: vec![("label".into(), Some("Original top calls only".into()))],
        needs_review: false,
    });
    let bottom =
        Molecule::dna("Independently ordered bottom", bottom_dna, Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent annealed top reference",
        top_dna,
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Actual unequal native source",
        vec![top.clone(), bottom.clone(), target.clone()],
        2,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual literal13/17 oligo annealing".into(),
        target_revision: target.revision(),
        inputs: vec![input("top", top), input("bottom", bottom)],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: AnnealOptions {
                min_paired_bases: 12,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        pin(&w, 2),
        vec![
            InputEventBinding {
                input_id: "top".into(),
                event: pin(&w, 0),
            },
            InputEventBinding {
                input_id: "bottom".into(),
                event: pin(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Saved actual13/17".into(),
            allow_incomplete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let event = saved.receipt.construct;
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &saved.workspace,
            &run_material::Selection {
                run_id: saved.receipt.run.run_id,
                step_index: 0,
                event: event.clone(),
            },
        )
        .unwrap(),
    );
    (saved.workspace, event, origin)
}
#[test]
fn bottom_seventeen_seeds_native_pcr_and_maps_only_real_original_top_annotations() {
    let (w, event, origin) = saved();
    let reopened = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    let source = &reopened.revisions[&event.revision];
    let before = workspace::save_workspace(&reopened).unwrap();
    let fresh = origin.replay(source).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), TOP);
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), BOTTOM);
    assert_eq!(fresh.tracked().duplex.bottom_start, 0);
    let ranked = collagene_core::primers::design_primers_view(
        collagene_core::analysis::TemplateView {
            sequence: &fresh.tracked().duplex.bottom,
            topology: Topology::Linear,
            template_strand: Strand::Bottom,
        },
        &collagene_core::primers::DesignOptions {
            target_start: 0,
            target_end: 17,
            min_length: 8,
            max_length: 8,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_tm_difference_c: 100.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        &source.revision(),
    )
    .unwrap();
    assert_eq!(ranked.pairs.len(), 1);
    assert_eq!(
        ranked.pairs[0].forward.primer.annealing.as_str(),
        "GTAGGTAG"
    );
    assert_eq!(
        ranked.pairs[0].reverse.primer.annealing.as_str(),
        "AGTCCGAT"
    );
    let f = Primer::new("Bottom F", "GTAGGTAG", "").unwrap();
    let r = Primer::new("Bottom R", "AGTCCGAT", "").unwrap();
    let previews = material_pcr::simulate(source, &origin, &f, &r, Strand::Bottom).unwrap();
    assert_eq!(previews.len(), 1);
    assert_eq!(previews[0].geometry.template_strand, Strand::Bottom);
    assert_eq!(previews[0].molecule.sequence.as_str(), BOTTOM);
    let feature = previews[0]
        .molecule
        .features
        .iter()
        .find(|f| f.label() == "Original top calls only")
        .unwrap();
    assert_eq!(feature.location.spans(), vec![(4, 17, true)]);
    // All13 original calls survive. Reversal changes coordinates/strand, not
    // feature completeness, so it preserves the original review state.
    let original = source
        .features
        .iter()
        .find(|f| f.label() == "Original top calls only")
        .unwrap();
    assert!(!original.needs_review);
    assert_eq!(feature.needs_review, original.needs_review);
    assert_eq!(
        feature
            .location
            .extract(&previews[0].molecule.sequence)
            .unwrap(),
        TOP
    );
    let prepared = d::prepare_install(
        &reopened,
        &d::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: event.clone(),
                origin: d::SourceOrigin::Saved(origin.clone()),
            }],
            derivation: d::Derivation::ExactPcr {
                forward: f,
                reverse: r,
                template_strand: Strand::Bottom,
                selection: d::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: previews[0].molecule.revision(),
                    geometry: previews[0].geometry.clone(),
                },
            },
            name: "Reusable actual bottom PCR".into(),
            save_primers: true,
        },
    )
    .unwrap();
    let result = d::reconstruct(&prepared.workspace, &prepared.receipt.product).unwrap();
    assert_eq!(result.reference().sequence.as_str(), BOTTOM);
    assert_eq!(result.tracked().duplex.top.len(), 17);
    assert_eq!(result.tracked().duplex.bottom.len(), 17);
    assert_eq!(result.tracked().duplex.left, Terminus::PCR);
    assert_eq!(result.tracked().duplex.right, Terminus::PCR);
    assert_eq!(
        modified_duplex::summary(result.tracked())
            .unwrap()
            .known_modified,
        0
    );
    assert!(prepared.receipt.primer_attachment.is_some());
    assert_eq!(prepared.workspace.revisions[&event.revision], *source);
    assert_eq!(workspace::save_workspace(&reopened).unwrap(), before);
}
#[test]
fn full_native_source_enumerates_both_thirteen_and_seventeen_without_target_selection() {
    let (w, event, origin) = saved();
    let source = &w.revisions[&event.revision];
    let f = Primer::new("F", "AGTC", "").unwrap();
    let r = Primer::new("R", "GTAG", "").unwrap();
    let previews = material_pcr::simulate(source, &origin, &f, &r, Strand::Bottom).unwrap();
    let mut sequences = previews
        .iter()
        .map(|p| p.molecule.sequence.as_str())
        .collect::<Vec<_>>();
    sequences.sort();
    let mut expected = vec!["GTAGCATCGGACT", BOTTOM];
    expected.sort();
    assert_eq!(sequences, expected);
    assert!(
        previews
            .iter()
            .all(|p| p.geometry.template_strand == Strand::Bottom)
    );
    let selected = &previews[0];
    assert!(
        d::prepare_install(
            &w,
            &d::Request {
                sources: vec![collagene_core::material_derivation::BoundSource {
                    role: collagene_core::material_derivation::SourceRole::Template,
                    event: event,
                    origin: d::SourceOrigin::Saved(origin.clone())
                }],
                derivation: d::Derivation::ExactPcr {
                    forward: f.clone(),
                    reverse: r.clone(),
                    template_strand: Strand::Bottom,
                    selection: d::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: selected.molecule.revision(),
                        geometry: selected.geometry.clone()
                    }
                },
                name: "Target cannot select one competing product".into(),
                save_primers: false
            }
        )
        .is_err()
    );
    assert!(matches!(
        material_pcr::simulate(source, &origin, &f, &r, Strand::Top),
        Err(collagene_core::Error::Unsupported(_))
    ));
}
#[test]
fn tails_are_synthesized_in_selected_frame_without_changing_the_source() {
    let (w, event, origin) = saved();
    let source = &w.revisions[&event.revision];
    let f = Primer::new("Tailed F", "GTAGGTAG", "GG").unwrap();
    let r = Primer::new("Tailed R", "AGTCCGAT", "TT").unwrap();
    let previews = material_pcr::simulate(source, &origin, &f, &r, Strand::Bottom).unwrap();
    assert_eq!(previews.len(), 1);
    assert_eq!(
        previews[0].molecule.sequence.as_str(),
        "GGGTAGGTAGCATCGGACTAA"
    );
    let feature = previews[0]
        .molecule
        .features
        .iter()
        .find(|f| f.label() == "Original top calls only")
        .unwrap();
    assert_eq!(feature.location.spans(), vec![(6, 19, true)]);
    let still = origin.replay(source).unwrap();
    assert_eq!(still.tracked().duplex.top.as_str(), TOP);
    assert_eq!(still.tracked().duplex.bottom.as_str(), BOTTOM);
}

#[test]
fn coverage_union_is_not_a_covalent_template_when_neither_strand_spans_the_product() {
    let (w, event, origin) = saved_pair("CCCCAGTCCGATGCTAC", BOTTOM);
    let source = &w.revisions[&event.revision];
    let fresh = origin.replay(source).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), "CCCCAGTCCGATGCTAC");
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), BOTTOM);
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    let forward = Primer::new("Left boundary", "CCCCAGTC", "").unwrap();
    let reverse = Primer::new("Right boundary", "GTAGGTAG", "").unwrap();
    for strand in [Strand::Top, Strand::Bottom] {
        assert!(matches!(
            material_pcr::simulate(source, &origin, &forward, &reverse, strand),
            Err(collagene_core::Error::Unsupported(_))
        ));
    }
}

#[test]
fn framed_step_pcr_orders_and_saved_collection_retain_bottom_declarations_after_restart() {
    let (mut w, event, origin) = saved();
    let source = w.revisions[&event.revision].clone();
    let target = Molecule::dna(
        "Independent full stored-bottom requirement",
        BOTTOM,
        Topology::Linear,
    )
    .unwrap();
    let target_id = "independent-bottom-target".to_owned();
    let target_event = workspace::RevisionEvent::new(
        &target_id,
        target.revision(),
        vec![],
        "Independent literal reference requirement",
    )
    .unwrap();
    let binding = EventBinding {
        document_id: target_id.clone(),
        event_id: target_event.id.clone(),
        revision: target.revision(),
    };
    w.revisions.insert(target.revision(), target.clone());
    w.history.push(target_event.clone());
    w.documents.push(workspace::WorkspaceDocument {
        id: target_id,
        head: target.revision(),
        head_event: Some(target_event.id),
        undo: vec![],
        redo: vec![],
        undo_events: vec![],
        redo_events: vec![],
    });
    let MaterialOrigin::Run(capsule) = origin else {
        unreachable!()
    };
    let plan = DesignPlan {
        name: "Actual Bottom-frame Step PCR with explicit order".into(),
        target_revision: target.revision(),
        inputs: vec![input("saved", source.clone())],
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin: capsule,
            },
            Operation::Pcr {
                template: Source::Step(0),
                forward: Primer::new("Actual Bottom F", "GTAGGTAG", "").unwrap(),
                reverse: Primer::new("Actual Bottom R", "AGTCCGAT", "").unwrap(),
                template_strand: Strand::Bottom,
            },
        ],
    };
    let report = collagene_core::plan_oligos::build(&plan, &target, &Default::default()).unwrap();
    assert!(report.complete);
    assert_eq!(report.uses.len(), 2);
    for usage in &report.uses {
        assert_eq!(usage.template_strand, Some(Strand::Bottom));
        assert!(matches!(usage.source, Some(Source::Step(0))));
    }
    let mut missing = serde_json::to_value(&report).unwrap();
    missing["uses"][0]
        .as_object_mut()
        .unwrap()
        .remove("template_strand");
    assert!(serde_json::from_value::<collagene_core::plan_oligos::Report>(missing).is_err());
    let mut bad = report.clone();
    bad.uses[0].template_strand = None;
    assert!(bad.validate().is_err());
    let committed = cloning_outputs::finalize(
        &w,
        &plan,
        binding,
        vec![InputEventBinding {
            input_id: "saved".into(),
            event: event.clone(),
        }],
        &cloning_outputs::Options {
            name: "Recorded actual Bottom PCR".into(),
            save_primers: true,
            allow_incomplete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&committed.workspace).unwrap())
            .unwrap();
    let hash = committed.receipt.primer_attachment.as_ref().unwrap();
    let table: collagene_core::oligo_formats::OligoTable =
        serde_json::from_slice(&reopened.attachments[hash].bytes).unwrap();
    assert_eq!(table.records.len(), 2);
    for row in &table.records {
        assert_eq!(row.original_fields[7], "Bottom");
    }
    let actual =
        collagene_core::run_lineage::replay_run(&reopened, &committed.receipt.run.run_id).unwrap();
    assert!(actual.products_match);
    assert!(actual.status_matches);
    assert_eq!(
        actual
            .current_report
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        BOTTOM
    );
    assert_eq!(reopened.revisions[&event.revision], source);
}
