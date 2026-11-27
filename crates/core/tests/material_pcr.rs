//! Actual native saved-material replay, not imported duplex-summary bootstrap.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest, cloning_outputs, homology_assembly, material_derivation as digest,
    material_pcr as m, pcr_outputs,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source},
    workspace::Workspace,
};
fn dna(name: &str, text: &str, topology: Topology) -> Molecule {
    Molecule::dna(name, text, topology).unwrap()
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn stock(incomplete: bool) -> (Workspace, Molecule, run_material::Capsule) {
    let (docs, plan, target) = if incomplete {
        let left = dna("Left actual stock", "AGTCCGATGCTA", Topology::Linear);
        let right = dna("Right actual stock", "GCTATTCGAATTCCAGTC", Topology::Linear);
        let target = dna(
            "Independent closed literal",
            "AGTCCGATGCTATTCGAATTCC",
            Topology::Circular,
        );
        let mut profile = collagene_core::overlap_profiles::ProfileOptions::default();
        profile.scan_repeated_overlaps = false;
        let plan = DesignPlan {
            name: "Finite terminal enumeration".into(),
            target_revision: target.revision(),
            inputs: vec![input("left", left.clone()), input("right", right.clone())],
            steps: vec![
                Operation::Pcr {
                    template: Source::Input("left".into()),
                    forward: Primer::new("LF", "AGTC", "").unwrap(),
                    reverse: Primer::new("LR", "TAGC", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                Operation::Pcr {
                    template: Source::Input("right".into()),
                    forward: Primer::new("RF", "GCTA", "").unwrap(),
                    reverse: Primer::new("RR", "GACT", "").unwrap(),
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
                    policy: homology_assembly::Options {
                        overlaps: vec![4, 4],
                        topology: Topology::Circular,
                        min_paired_bases: 4,
                        max_comparisons: 142,
                        ..Default::default()
                    },
                },
            ],
        };
        (vec![left, right, target], plan, 2)
    } else {
        let source = dna(
            "Actual original stock",
            "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
            Topology::Linear,
        );
        let target = dna(
            "Independent original full PCR",
            "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
            Topology::Linear,
        );
        let plan = DesignPlan {
            name: "Real original PCR".into(),
            target_revision: target.revision(),
            inputs: vec![input("stock", source.clone())],
            steps: vec![Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("Original F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Original R", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            }],
        };
        (vec![source, target], plan, 1)
    };
    let w = Workspace::from_documents("Original native source", docs, target).unwrap();
    let bindings = plan
        .inputs
        .iter()
        .enumerate()
        .map(|(i, p)| InputEventBinding {
            input_id: p.id.clone(),
            event: event(&w, i),
        })
        .collect();
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, target),
        bindings,
        &cloning_outputs::Options {
            allow_incomplete: incomplete,
            ..Default::default()
        },
    )
    .unwrap();
    let source = saved.workspace.revisions[&saved.receipt.construct.revision].clone();
    let origin = run_material::capsule(
        &saved.workspace,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id,
            step_index: plan.steps.len() - 1,
            event: saved.receipt.construct,
        },
    )
    .unwrap();
    (saved.workspace, source, origin)
}
fn cut(incomplete: bool) -> (Molecule, MaterialOrigin) {
    let (w, _, origin) = stock(incomplete);
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let installed = digest::prepare_install(
        &w,
        &digest::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: origin.selection.event.clone(),
                origin: digest::SourceOrigin::Saved(
                    collagene_core::run_material::MaterialOrigin::Run(origin),
                ),
            }],
            derivation: digest::Derivation::CatalogDigest {
                selection: catalog_digest::select(&c, &[id]).unwrap(),
                options: Default::default(),
                fragment_index: if incomplete { 0 } else { 1 },
            },
            name: "Native selected digest fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let event = installed.receipt.product;
    let source = installed.workspace.revisions[&event.revision].clone();
    let origin = MaterialOrigin::Derived(digest::capsule(&installed.workspace, &event).unwrap());
    (source, origin)
}
fn oligos(incomplete: bool) -> (Primer, Primer) {
    (
        Primer::new(
            "New F",
            if incomplete { "AATTCCAG" } else { "AATTCCGT" },
            "",
        )
        .unwrap(),
        Primer::new(
            "New R",
            if incomplete { "CGAATAGC" } else { "TGACGATC" },
            "",
        )
        .unwrap(),
    )
}
#[test]
fn actual_run_origin_previews_and_selected_outputs_recompute_with_full_origin_binding() {
    let (_, source, origin) = stock(false);
    let origin = MaterialOrigin::Run(origin);
    let forward = Primer::new("New F", "AGTCCGAT", "").unwrap();
    let reverse = Primer::new("New R", "TGACGATC", "").unwrap();
    let previews = m::simulate(
        &source,
        &origin,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(previews.len(), 1);
    assert_eq!(
        previews[0].molecule.sequence.as_str(),
        "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    let evidence: m::OriginEvidence =
        serde_json::from_str(&previews[0].molecule.metadata[m::ORIGIN_METADATA]).unwrap();
    assert!(!evidence.origin_incomplete);
    assert_eq!(evidence.selection, origin.selection());
    let request = pcr_outputs::Request {
        source,
        forward,
        reverse,
        selected_index: 0,
        expected_product_revision: previews[0].molecule.revision(),
        options: pcr_outputs::Options {
            name: "New selected PCR".into(),
            save_product: true,
            save_primers: true,
        },
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let context = m::outputs_revision(&origin, &request).unwrap();
    let selected = m::prepare_outputs(&origin, &request).unwrap();
    selected.validate_binding(&context).unwrap();
    assert_eq!(selected.product.as_ref().unwrap().name, "New selected PCR");
    assert_eq!(selected.primers.unwrap().records.len(), 2);
    assert!(
        selected
            .product
            .unwrap()
            .metadata
            .contains_key(m::ORIGIN_METADATA)
    );
    let mut stale = request.clone();
    stale.expected_product_revision = "00".repeat(32);
    assert!(m::prepare_outputs(&origin, &stale).is_err());
    let mut renamed = request;
    renamed.options.name = "Different save choice".into();
    assert_ne!(context, m::outputs_revision(&origin, &renamed).unwrap());
}
#[test]
fn actual_staggered_digest_is_used_without_stripping_source_markers_or_trusting_reference_only_input()
 {
    let (source, origin) = cut(false);
    let before = source.clone();
    let (forward, reverse) = oligos(false);
    assert!(
        verification::simulate_pcr_documents(
            &source,
            &forward,
            &reverse,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    let previews = m::simulate(
        &source,
        &origin,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(source, before);
    assert_eq!(previews.len(), 1);
    assert_eq!(
        previews[0].molecule.sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &previews[0].molecule
    ));
    let request = pcr_outputs::Request {
        source,
        forward,
        reverse,
        selected_index: 0,
        expected_product_revision: previews[0].molecule.revision(),
        options: pcr_outputs::Options {
            name: "Selected staggered-template PCR".into(),
            save_product: true,
            save_primers: false,
        },
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    assert!(pcr_outputs::prepare(&request).is_err());
    assert!(
        m::prepare_outputs(&origin, &request)
            .unwrap()
            .product
            .is_some()
    );
}
#[test]
fn reusable_but_incomplete_actual_origin_remains_explicitly_conditional_in_preview_and_output() {
    let (source, origin) = cut(true);
    let (forward, reverse) = oligos(true);
    let previews = m::simulate(
        &source,
        &origin,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(previews.len(), 1);
    assert_eq!(
        previews[0].molecule.sequence.as_str(),
        "AATTCCAGTCCGATGCTATTCG"
    );
    let evidence: m::OriginEvidence =
        serde_json::from_str(&previews[0].molecule.metadata[m::ORIGIN_METADATA]).unwrap();
    assert!(evidence.origin_incomplete);
    assert!(evidence.origin_findings_total >= evidence.origin_findings.len());
    assert!(
        previews[0]
            .findings
            .iter()
            .any(|s| s.contains("remains incomplete"))
    );
    let request = pcr_outputs::Request {
        source,
        forward,
        reverse,
        selected_index: 0,
        expected_product_revision: previews[0].molecule.revision(),
        options: pcr_outputs::Options {
            name: "Conditional PCR".into(),
            save_product: true,
            save_primers: false,
        },
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let selected = m::prepare_outputs(&origin, &request).unwrap();
    let saved: m::OriginEvidence =
        serde_json::from_str(&selected.product.unwrap().metadata[m::ORIGIN_METADATA]).unwrap();
    assert!(saved.origin_incomplete);
}
#[test]
fn changed_source_origin_and_unbounded_oligo_are_rejected_before_preview_selection() {
    let (mut source, origin) = cut(false);
    let (forward, reverse) = oligos(false);
    source.name = "Another revision".into();
    assert!(
        m::revision(
            &source,
            &origin,
            &forward,
            &reverse,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    let (source, mut origin) = cut(false);
    if let MaterialOrigin::Derived(d) = &mut origin {
        d.record_bytes[44] ^= 1; // Actual opaque recipe checksum must fail before any PCR selection.
    }
    assert!(
        m::simulate(
            &source,
            &origin,
            &forward,
            &reverse,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    let (_, _, run) = stock(false);
    assert!(
        m::simulate(
            &source,
            &MaterialOrigin::Run(run),
            &forward,
            &reverse,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
    let (source, origin) = cut(false);
    let long = Primer::new("Bounded F", &"A".repeat(1001), "").unwrap();
    assert!(
        m::revision(
            &source,
            &origin,
            &long,
            &reverse,
            collagene_core::chemistry::Strand::Top
        )
        .is_err()
    );
}
