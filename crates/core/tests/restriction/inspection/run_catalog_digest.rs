//! Literal saved materials are freshly replayed; displayed chemistry is never recovered from metadata.
use collagene_core::{
    analysis::Primer,
    catalog_digest,
    chemistry::{EndTreatment, FivePrime, Terminus},
    cloning_outputs::{self, Options},
    digest_cut_selection::{self, Boundary, Choice},
    methylation::MethylationState,
    modified_duplex::GenerationAssumptions,
    restriction_inspection::InspectionOptions,
    run_catalog_digest as d,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Capsule, Selection},
    verification::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let doc = &w.documents[i];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn saved(sources: Vec<Molecule>, target: Molecule, steps: Vec<Operation>) -> (Molecule, Capsule) {
    let inputs = sources
        .iter()
        .enumerate()
        .map(|(i, m)| PlanInput {
            id: format!("input-{i}"),
            revision: m.revision(),
            molecule: m.clone(),
        })
        .collect();
    let plan = DesignPlan {
        name: "Independent source preparation".into(),
        target_revision: target.revision(),
        inputs,
        steps,
    };
    let n = sources.len();
    let mut docs = sources;
    docs.push(target);
    let w = Workspace::from_documents("Source and intent", docs, n).unwrap();
    let expected = &w.revisions[&w.documents[n].head];
    let reconstruction = collagene_core::verification::verify(&plan, expected);
    assert_eq!(
        reconstruction.products.len(),
        plan.steps.len(),
        "Origin preparation failed before capture: {:?}",
        reconstruction.findings
    );
    assert!(
        reconstruction.products.last().is_some_and(|product| {
            collagene_core::verification::plan_target_matches(&plan, product, expected)
        }),
        "Origin differs from its independent literal target: {:?}",
        reconstruction.findings
    );
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, n),
        (0..n)
            .map(|i| InputEventBinding {
                input_id: format!("input-{i}"),
                event: binding(&w, i),
            })
            .collect(),
        &Options {
            name: "Named saved source".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let event = commit.receipt.construct;
    let source = commit.workspace.revisions[&event.revision].clone();
    let origin = run_material::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: plan.steps.len() - 1,
            event,
        },
    )
    .unwrap();
    (source, origin)
}
fn select(
    names: &[&str],
) -> (
    catalog_digest::catalog::Catalog,
    catalog_digest::CatalogSelection,
) {
    let c = catalog_digest::builtin_catalog().unwrap();
    let ids = names
        .iter()
        .map(|name| {
            c.edition()
                .enzymes
                .iter()
                .find(|e| e.name == *name)
                .unwrap()
                .id
                .clone()
        })
        .collect::<Vec<_>>();
    let s = catalog_digest::select(&c, &ids).unwrap();
    (c, s)
}
fn pcr() -> (Molecule, Capsule) {
    saved(
        vec![Molecule::dna("Original6bp", "AAACCC", Topology::Linear).unwrap()],
        Molecule::dna(
            "Independent22bp",
            "CCGAATTCAAACCCCTCGAGCC",
            Topology::Linear,
        )
        .unwrap(),
        vec![Operation::Pcr {
            template: Source::Input("input-0".into()),
            forward: Primer::new("F", "AAA", "CCGAATTC").unwrap(),
            reverse: Primer::new("R", "GGG", "GGCTCGAG").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    )
}
#[test]
fn named_saved_pcr_digest_retains_literal_cut_order_and_real_outer_hydroxyls() {
    let (source, origin) = pcr();
    let (c, s) = select(&["EcoRI", "XhoI"]);
    let options = InspectionOptions::default();
    let before = serde_json::to_vec(&source).unwrap();
    let report = d::inspect(&c, &s, "saved-source", &source, &origin, &options).unwrap();
    report
        .validate_bindings_for(&c, &s, "saved-source", &source, &origin, &options)
        .unwrap();
    let inner = &report.inspection.inspection;
    assert!(inner.scan_complete && inner.products_determined);
    assert_eq!(inner.sites.len(), 2);
    assert_eq!(
        inner
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![3, 12, 7]
    );
    let outer_left = &inner.fragments[0].duplex.ends.as_ref().unwrap().0;
    let outer_right = &inner.fragments[2].duplex.ends.as_ref().unwrap().1;
    assert_eq!(outer_left.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(outer_right.chemistry.five_prime, FivePrime::Hydroxyl);
    let (left, right) = inner.fragments[1].duplex.ends.as_ref().unwrap();
    assert_eq!((&*left.overhang, &*right.overhang), ("AATT", "TCGA"));
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(report.fragment_modifications.len(), 3);
    assert_eq!(report.source_modifications.known_modified, 0);
    let choice = Choice {
        first: Boundary::Cut { site_index: 0 },
        second: Boundary::Cut { site_index: 1 },
        complement: false,
    };
    assert_eq!(
        digest_cut_selection::selected_fragment(&report.inspection, &choice)
            .unwrap()
            .index,
        1
    );
    assert_eq!(before, serde_json::to_vec(&source).unwrap());
    let decoded: d::Inspection =
        serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    decoded.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn saved_circle_keeps_the_cut_overhangs_instead_of_unknown_reference_ends() {
    let (source, origin) = saved(
        vec![Molecule::dna("Original12bp", "GGAATTCTTACG", Topology::Linear).unwrap()],
        Molecule::dna("Independent circle", "GGAATTCTTACG", Topology::Circular).unwrap(),
        vec![
            Operation::DeclareSourceTermini {
                template: Source::Input("input-0".into()),
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    );
    let (c, s) = select(&["EcoRI"]);
    let report = d::inspect(
        &c,
        &s,
        "saved-circle",
        &source,
        &origin,
        &InspectionOptions::default(),
    )
    .unwrap();
    let inner = &report.inspection.inspection;
    assert_eq!(inner.topology, Topology::Circular);
    assert_eq!(inner.fragments.len(), 1);
    assert_eq!(inner.fragments[0].top_start, 2);
    let (left, right) = inner.fragments[0].duplex.ends.as_ref().unwrap();
    assert_eq!((&*left.overhang, &*right.overhang), ("AATT", "AATT"));
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
}
#[test]
fn inherited_cpg_positions_survive_conditional_not_i_cleavage_without_claiming_permission() {
    let source = Molecule::dna("Original CpG circle", "AAGCGGCCGCAA", Topology::Circular).unwrap();
    let target = Molecule::dna(
        "Independent blocked circle",
        "AAGCGGCCGCAA",
        Topology::Circular,
    )
    .unwrap();
    let (source, origin) = saved(
        vec![source],
        target,
        vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("input-0".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::Unmodified,
                    dcm: MethylationState::Unmodified,
                    cpg: MethylationState::FullyModified,
                },
            },
            Operation::CheckedDigest {
                selection: checked_selection(vec!["NotI".into()]),
                options: collagene_core::modified_duplex::checked_inspection_options(),
                template: Source::Step(0),
                fragment_index: 0,
            },
        ],
    );
    let (c, s) = select(&["NotI"]);
    let report = d::inspect(
        &c,
        &s,
        "modified-circle",
        &source,
        &origin,
        &InspectionOptions::default(),
    )
    .unwrap();
    assert_eq!(report.source_modifications.known_modified, 4);
    assert_eq!(report.source_modifications.generation.len(), 1);
    assert!(report.inspection.inspection.products_determined);
    assert_eq!(
        report
            .fragment_modifications
            .iter()
            .map(|m| m.known_modified)
            .sum::<usize>(),
        4
    );
    assert!(
        report
            .fragment_modifications
            .iter()
            .all(|m| m.generation == report.source_modifications.generation)
    );
    assert!(
        report
            .findings
            .iter()
            .any(|s| s.contains("does not assess methylation susceptibility"))
    );
}
#[test]
fn actual_retained_nicks_withhold_products_but_keep_recognition_scan() {
    let (source, origin) = saved(
        vec![
            Molecule::dna("Left6bp", "AAACCC", Topology::Linear).unwrap(),
            Molecule::dna("Right6bp", "GGGTTT", Topology::Linear).unwrap(),
        ],
        Molecule::dna(
            "Independent nicked circle",
            "AAACCCGGGTTT",
            Topology::Circular,
        )
        .unwrap(),
        vec![
            Operation::Pcr {
                template: Source::Input("input-0".into()),
                forward: Primer::new("LF", "AAA", "").unwrap(),
                reverse: Primer::new("LR", "GGG", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            // Each junction needs at least one real seal. Only the left PCR
            // fragment is phosphorylated, leaving one genuine nick per junction.
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Pcr {
                template: Source::Input("input-1".into()),
                forward: Primer::new("RF", "GGG", "").unwrap(),
                reverse: Primer::new("RR", "AAA", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
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
    );
    let (c, s) = select(&["SmaI"]);
    let report = d::inspect(
        &c,
        &s,
        "nicked-circle",
        &source,
        &origin,
        &InspectionOptions::default(),
    )
    .unwrap();
    let inner = &report.inspection.inspection;
    assert_eq!(report.source_duplex.nicks.len(), 2);
    assert!(!report.source_duplex.covalently_closed_in_model);
    assert!(inner.scan_complete);
    assert_eq!(inner.sites.len(), 1);
    assert!(!inner.products_determined);
    assert!(inner.fragments.is_empty() && report.fragment_modifications.is_empty());
    assert!(inner.diagnostics.iter().any(|s| s.contains("nicks")));
}
#[test]
fn partial_scans_stale_origins_and_detached_modification_summaries_are_rejected_or_withheld() {
    let (mut source, mut origin) = pcr();
    let (c, s) = select(&["EcoRI", "XhoI"]);
    let mut options = InspectionOptions::default();
    let expected =
        d::inspection_revision(&c, &s, "saved-source", &source, &origin, &options).unwrap();
    let report = d::inspect(&c, &s, "saved-source", &source, &origin, &options).unwrap();
    let mut tampered = report.clone();
    tampered.fragment_modifications.pop();
    assert!(tampered.validate_binding(&expected).is_err());
    options.max_comparisons = 1;
    let partial = d::inspect(&c, &s, "saved-source", &source, &origin, &options).unwrap();
    assert!(!partial.inspection.inspection.scan_complete);
    assert!(partial.fragment_modifications.is_empty());
    assert!(report.validate_binding(&partial.context_revision).is_err());
    origin.selection.event.revision = "0".repeat(64);
    assert!(d::source_preflight(&c, &s, "saved-source", &source, &origin, &options).is_err());
    let (_, origin) = pcr();
    source.name.push_str(" edited");
    assert!(d::source_preflight(&c, &s, "saved-source", &source, &origin, &options).is_err());
}

#[test]
fn saved_staggered_source_uses_actual_bottom_free_ends_for_both_cut_paths() {
    use collagene_core::chemistry::AnnealOptions;
    // Two ordered21mers pair over17bases. The top has a left AATT extension;
    // the bottom has its own right AATT extension, both written5′→3′ here.
    let top = Molecule::dna("Top21mer", "AATTGCTACGGAATTCTTACG", Topology::Linear).unwrap();
    let bottom = Molecule::dna("Bottom21mer", "AATTCGTAAGAATTCCGTAGC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent top21mer",
        "AATTGCTACGGAATTCTTACG",
        Topology::Linear,
    )
    .unwrap();
    let (source, origin) = saved(
        vec![top, bottom],
        target,
        vec![Operation::AnnealOligos {
            top: Source::Input("input-0".into()),
            bottom: Source::Input("input-1".into()),
            options: AnnealOptions {
                min_paired_bases: 17,
                top_chemistry: Terminus::CLEAVED,
                bottom_chemistry: Terminus::CLEAVED,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    );
    let (c, s) = select(&["EcoRI"]);
    let report = d::inspect(
        &c,
        &s,
        "staggered-source",
        &source,
        &origin,
        &InspectionOptions::default(),
    )
    .unwrap();
    assert_eq!(
        (
            report.source_duplex.top_length,
            report.source_duplex.bottom_start,
            report.source_duplex.bottom_length
        ),
        (21, 4, 21)
    );
    assert_eq!(
        report.source_duplex.bottom_start + report.source_duplex.bottom_length as i64,
        25
    );
    assert_eq!(report.inspection.inspection.sites.len(), 1);
    let left = Choice {
        first: Boundary::SourceStart,
        second: Boundary::Cut { site_index: 0 },
        complement: false,
    };
    let right = Choice {
        first: Boundary::Cut { site_index: 0 },
        second: Boundary::SourceEnd,
        complement: false,
    };
    let first = report.selected_fragment(&left).unwrap();
    assert_eq!(
        (
            first.index,
            first.top_start,
            first.bottom_start,
            first.duplex.top_length,
            first.duplex.bottom_length
        ),
        (0, 0, 4, 11, 11)
    );
    let last = report.selected_fragment(&right).unwrap();
    assert_eq!(
        (
            last.index,
            last.top_start,
            last.bottom_start,
            last.duplex.top_length,
            last.duplex.bottom_length
        ),
        (1, 11, 15, 10, 10)
    );
    assert!(digest_cut_selection::selected_fragment(&report.inspection, &left).is_err());
    assert!(digest_cut_selection::selected_fragment(&report.inspection, &right).is_err());
    let mut tampered = report.clone();
    tampered.source_duplex.bottom_start = i64::MAX;
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    let mut tampered = report.clone();
    tampered.source_duplex.bottom_length -= 1;
    assert!(tampered.validate_binding(&report.context_revision).is_err());
}

fn checked_selection(names: Vec<String>) -> collagene_core::catalog_digest::CatalogSelection {
    let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
    let ids = names
        .iter()
        .map(|name| {
            catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.name.eq_ignore_ascii_case(name))
                .unwrap()
                .id
                .clone()
        })
        .collect::<Vec<_>>();
    collagene_core::catalog_digest::select(&catalog, &ids).unwrap()
}
