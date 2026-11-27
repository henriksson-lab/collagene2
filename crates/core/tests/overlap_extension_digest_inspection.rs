//! Target-free scientific inspection retains actual OE material and uncertainty.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    assembly_recipe, chemistry, cloning_outputs, overlap_extension,
    overlap_extension_digest_inspection::{self as inspect, Request},
    overlap_extension_pcr::{PrimerPair, Request as Preparation},
    restriction_inspection::InspectionOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EDITED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
fn source() -> Molecule {
    Molecule::dna("Immutable source40", ORIGINAL, Topology::Linear).unwrap()
}
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("Forward", f, ft).unwrap(),
        reverse: Primer::new("Reverse", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn request() -> Request {
    let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
    Request {
        source_role: "insert".into(),
        preparation: Preparation {
            left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
            right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
            outer: pair("AGTCCGAT", "ATGCGAATTC", "TGACGATC", "ATGCGAATTC"),
            extension: overlap_extension::Options {
                overlap: 18,
                min_paired_bases: 8,
                ..Default::default()
            },
        },
        selection: collagene_core::catalog_digest::select(&catalog, &["ecori".into()]).unwrap(),
        options: InspectionOptions::default(),
    }
}
#[test]
fn three_native_singleton_pcrs_and_extension_supply_real_digest_cuts_and_annotation_mapping() {
    let mut source = source();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("left mark".into()))],
        needs_review: false,
    });
    let immutable = serde_json::to_value(&source).unwrap();
    let request = request();
    let report = inspect::inspect(&source, None, &request, None).unwrap();
    report
        .validate_binding(&inspect::inspection_revision(&source, None, &request, None).unwrap())
        .unwrap();
    assert!(!report.incomplete);
    assert!(report.selectable());
    assert_eq!(
        report
            .pcr_stages
            .iter()
            .map(|p| p.product_length.unwrap())
            .collect::<Vec<_>>(),
        vec![26, 30, 58]
    );
    let extension = report.extension.as_ref().unwrap();
    assert_eq!(extension.fill.duplex.top.as_str(), EDITED);
    let amplicon = report.amplicon.as_ref().unwrap();
    assert_eq!(
        amplicon.sequence.as_str(),
        "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT"
    );
    assert!(
        amplicon
            .features
            .iter()
            .any(|f| f.label() == "left mark" && f.location == Location::range(10, 14))
    );
    let digest = report.inspection.as_ref().unwrap();
    assert!(digest.products_determined);
    assert_eq!(
        digest
            .sites
            .iter()
            .map(|s| (s.top_cut, s.bottom_cut))
            .collect::<Vec<_>>(),
        vec![(Some(5), Some(9)), (Some(49), Some(53))]
    );
    assert_eq!(
        digest
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 44, 9]
    );
    assert!(report.origin_revision.is_none());
    assert_eq!(serde_json::to_value(&source).unwrap(), immutable);
}
#[test]
fn competing_source_or_outer_pcr_geometries_never_select_a_desired_product() {
    let repeated = Molecule::dna(
        "Two source copies",
        &format!("{ORIGINAL}{ORIGINAL}"),
        Topology::Linear,
    )
    .unwrap();
    let report = inspect::inspect(&repeated, None, &request(), None).unwrap();
    assert!(report.incomplete);
    assert!(!report.selectable());
    assert!(report.amplicon.is_none());
    assert!(report.extension.is_none());
    assert!(report.inspection.is_none());
    assert!(report.pcr_stages[0].geometries.len() > 1);
    let mut request = request();
    request.preparation.outer = pair("T", "", "A", "");
    let report = inspect::inspect(&source(), None, &request, None).unwrap();
    assert!(report.incomplete);
    assert!(report.extension.is_some());
    assert_eq!(report.pcr_stages.len(), 3);
    assert!(report.pcr_stages[2].geometries.len() > 1);
    assert!(report.amplicon.is_none());
    assert!(report.inspection.is_none());
}
#[test]
fn valid_selected_terminal_geometry_remains_conditional_but_actual_digest_is_selectable() {
    let mut request = request();
    request.preparation.extension.max_comparisons = 18;
    let report = inspect::inspect(&source(), None, &request, None).unwrap();
    assert!(report.incomplete);
    assert!(
        !report
            .extension
            .as_ref()
            .unwrap()
            .terminal_pairings
            .complete
    );
    assert!(report.amplicon.is_some());
    assert!(report.selectable());
    assert!(report.inspection.as_ref().unwrap().products_determined);
    let mut tampered = report.clone();
    tampered.incomplete = false;
    assert!(tampered.validate().is_err());
    request.preparation.extension.max_comparisons = 20_000_000;
    request.options.max_comparisons = 1;
    let capped = inspect::inspect(&source(), None, &request, None).unwrap();
    assert!(capped.incomplete);
    assert!(!capped.selectable());
    assert!(capped.amplicon.is_some());
}
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved(incomplete: bool) -> (Molecule, MaterialOrigin, Workspace) {
    let source = source();
    let mut steps = vec![];
    let target = if incomplete {
        let mut p = request().preparation;
        p.outer.forward.tail.clear();
        p.outer.reverse.tail.clear();
        p.extension.max_comparisons = 18;
        steps.extend([
            Operation::Pcr {
                template: Source::Input("original".into()),
                forward: p.left.forward,
                reverse: p.left.reverse,
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("original".into()),
                forward: p.right.forward,
                reverse: p.right.reverse,
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::OverlapExtension {
                left: collagene_core::verification::OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                },
                right: collagene_core::verification::OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                },
                options: p.extension,
            },
            Operation::Pcr {
                template: Source::Step(2),
                forward: p.outer.forward,
                reverse: p.outer.reverse,
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        ]);
        Molecule::dna("Independent originating edited38", EDITED, Topology::Linear).unwrap()
    } else {
        steps.push(Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        });
        Molecule::dna(
            "Independent originating stock40",
            ORIGINAL,
            Topology::Linear,
        )
        .unwrap()
    };
    let plan = DesignPlan {
        name: "Closed originating PCR graph".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps,
    };
    let w = Workspace::from_documents("Origin", vec![source, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved exact material".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: incomplete,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
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
    (source, MaterialOrigin::Run(origin), commit.workspace)
}
#[test]
fn exact_saved_origin_is_replayed_and_its_original_incompleteness_and_target_are_retained() {
    for incomplete in [false, true] {
        let (source, origin, workspace) = saved(incomplete);
        let pin = serde_json::to_value(&origin).unwrap();
        let request = request();
        assert!(chemistry::has_unmapped_physical_evidence(&source));
        assert!(inspect::source_preflight(&source, None, &request, None).is_err());
        let report = inspect::inspect(&source, Some(&origin), &request, Some(&workspace)).unwrap();
        assert_eq!(report.incomplete, incomplete);
        assert!(report.selectable());
        assert_eq!(report.origin_incomplete, Some(incomplete));
        let evidence = report.origin_evidence.as_ref().unwrap();
        assert_eq!(evidence.selection, origin.selection());
        assert_eq!(&evidence.source, origin.event());
        assert_eq!(evidence.incomplete, incomplete);
        assert!(evidence.findings_total >= evidence.findings.len());
        assert_eq!(
            report.origin_revision,
            Some(
                blake3::hash(&serde_json::to_vec(&origin).unwrap())
                    .to_hex()
                    .to_string()
            )
        );
        assert_eq!(serde_json::to_value(&origin).unwrap(), pin);
        let mut altered = source.clone();
        altered.name.push('!');
        assert!(
            inspect::source_preflight(&altered, Some(&origin), &request, Some(&workspace)).is_err()
        );
    }
}
#[test]
fn documentary_stage_catalog_context_and_source_binding_tampering_is_rejected() {
    let report = inspect::inspect(&source(), None, &request(), None).unwrap();
    let decoded: inspect::Preview =
        serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    decoded.validate_binding(&report.context_revision).unwrap();
    assert!(decoded.validate_binding(&"0".repeat(64)).is_err());
    for kind in 0..5 {
        let mut bad = report.clone();
        match kind {
            0 => bad.pcr_stages[0].product_length = Some(999),
            1 => bad.request.selection.catalog_content_hash = "0".repeat(64),
            2 => bad.pcr_stages[2].template_revision = "0".repeat(64),
            3 => bad.extension = None,
            _ => bad.source_revision = "0".repeat(64),
        }
        assert!(
            bad.validate_binding(&report.context_revision).is_err(),
            "tamper{kind}"
        );
    }
    // The exact source/pinned request also compiles as a normal native graph.
    let prep = assembly_recipe::Preparation::OverlapExtensionPcrCatalogDigest {
        request: request().preparation,
        selection: request().selection,
        fragment_index: 1,
        repair: None,
    };
    assert!(matches!(
        prep,
        assembly_recipe::Preparation::OverlapExtensionPcrCatalogDigest {
            fragment_index: 1,
            ..
        }
    ));
}

#[test]
fn derived_origin_owns_its_exact_event_and_inherits_limits_without_a_target_verdict() {
    use collagene_core::{material_derivation as d, workspace};
    let (_source, origin, w) = saved(true);
    let next = d::prepare_install(
        &w,
        &d::Request {
            source: origin.event().clone(),
            origin: d::SourceOrigin::Saved(origin),
            derivation: d::Derivation::EndRepair {
                method: chemistry::EndRepairMethod::Blunt,
            },
            name: "Exact native repaired descendant".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let event = next.receipt.product;
    let w =
        workspace::load_workspace(&workspace::save_workspace(&next.workspace).unwrap()).unwrap();
    let origin = MaterialOrigin::Derived(d::capsule(&w, &event).unwrap());
    let source = &w.revisions[&event.revision];
    let r = inspect::inspect(source, Some(&origin), &request(), Some(&w)).unwrap();
    r.validate_binding(&r.context_revision).unwrap();
    assert_eq!(r.origin_incomplete, Some(true));
    assert!(r.incomplete);
    assert!(r.selectable());
    assert_eq!(
        r.origin_evidence.as_ref().unwrap().selection,
        origin.selection()
    );
    assert!(!r.origin_evidence.as_ref().unwrap().findings.is_empty());
    assert!(inspect::inspect(source, Some(&origin), &request(), None).is_err());
    let mut foreign = w.clone();
    foreign
        .attachments
        .retain(|_, a| a.kind != d::WORKSPACE_ATTACHMENT_KIND);
    assert!(inspect::source_preflight(source, Some(&origin), &request(), Some(&foreign)).is_err());
    // The old source was never renamed, rewritten, or adopted as the derived event.
    assert_ne!(source.revision(), self::source().revision());
}

#[test]
fn bottom_only_source_calls_reach_all_six_oligos_and_real_ecori_ends() {
    use collagene_core::{
        assembly_edit_intent::Edit, chemistry::Strand, overlap_extension_primer_design as planner,
        overlap_primer_design::Fragment,
    };
    const BOTTOM: &str = "GGGGTGACGATCTACAGGTCTGAAGTCCTACGTAGCATCGGACT";
    let top = source();
    let bottom = Molecule::dna("Independently ordered bottom44", BOTTOM, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent paired top40", ORIGINAL, Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Native bottom44 source",
        vec![top.clone(), bottom.clone(), target.clone()],
        2,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual40/44 annealing".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "top".into(),
                revision: top.revision(),
                molecule: top,
            },
            PlanInput {
                id: "bottom".into(),
                revision: bottom.revision(),
                molecule: bottom,
            },
        ],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: chemistry::AnnealOptions {
                min_paired_bases: 40,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    };
    let installed = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 2),
        vec![
            InputEventBinding {
                input_id: "top".into(),
                event: binding(&w, 0),
            },
            InputEventBinding {
                input_id: "bottom".into(),
                event: binding(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Saved native40/44".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let event = installed.receipt.construct;
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &installed.workspace,
            &Selection {
                run_id: installed.receipt.run.run_id,
                step_index: 0,
                event: event.clone(),
            },
        )
        .unwrap(),
    );
    let source = &installed.workspace.revisions[&event.revision];
    let mut o = planner::Options::default();
    o.homology_flank = 4;
    o.min_paired_bases = 8;
    o.outer_forward_tail = "ATGCGAATTC".into();
    o.outer_reverse_tail = "ATGCGAATTC".into();
    o.primer_options.min_length = 8;
    o.primer_options.max_length = 8;
    o.primer_options.max_pairs = 1;
    o.primer_options.min_gc_fraction = 0.;
    o.primer_options.max_gc_fraction = 1.;
    o.primer_options.max_tm_difference_c = 100.;
    o.primer_options.max_homopolymer = 200;
    o.primer_options.max_three_prime_complementarity = 200;
    let proposal = planner::design(
        source,
        &planner::Request {
            fragment: Fragment {
                role: "insert".into(),
                input_id: "source".into(),
                start: 0,
                end: 44,
                reverse: false,
            },
            edits: vec![Edit {
                role: "insert".into(),
                start: 16,
                end: 20,
                expected: "GGTC".into(),
                replacement: "AA".into(),
            }],
            options: o,
            template_strand: Strand::Bottom,
            material_binding: Some(planner::BoundOrigin {
                event,
                origin: origin.clone(),
            }),
        },
        Some(&installed.workspace),
    )
    .unwrap();
    let r = Request {
        preparation: proposal.preparation,
        ..request()
    };
    let before = collagene_core::workspace::save_workspace(&installed.workspace).unwrap();
    let preview = inspect::inspect(source, Some(&origin), &r, Some(&installed.workspace)).unwrap();
    preview.validate_binding(&preview.context_revision).unwrap();
    assert!(!preview.incomplete, "{:?}", preview.findings);
    assert!(preview.selectable());
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "ATGCGAATTCGGGGTGACGATCTACAAATGAAGTCCTACGTAGCATCGGACTGAATTCGCAT"
    );
    assert!(
        preview.pcr_stages[0]
            .geometries
            .iter()
            .all(|g| g.template_strand == Strand::Bottom)
    );
    assert!(
        preview.pcr_stages[1]
            .geometries
            .iter()
            .all(|g| g.template_strand == Strand::Bottom)
    );
    assert!(
        preview.pcr_stages[2]
            .geometries
            .iter()
            .all(|g| g.template_strand == Strand::Top)
    );
    let digest = preview.inspection.as_ref().unwrap();
    assert_eq!(
        digest
            .sites
            .iter()
            .map(|s| (s.top_cut, s.bottom_cut))
            .collect::<Vec<_>>(),
        vec![(Some(5), Some(9)), (Some(53), Some(57))]
    );
    assert_eq!(
        digest
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 48, 9]
    );
    let ends = digest.fragments[1].duplex.ends.as_ref().unwrap();
    assert_eq!(ends.0.overhang, "AATT");
    assert_eq!(ends.1.overhang, "AATT");
    assert_eq!(ends.0.chemistry.five_prime, chemistry::FivePrime::Phosphate);
    assert_eq!(
        ends.1.chemistry.three_prime,
        chemistry::ThreePrime::Hydroxyl
    );
    assert_eq!(
        collagene_core::workspace::save_workspace(&installed.workspace).unwrap(),
        before
    );
}
