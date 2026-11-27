//! Target-free scientific inspection retains actual OE material and uncertainty.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    assembly_recipe, chemistry, cloning_outputs, overlap_extension,
    overlap_extension_digest_inspection::{self as inspect, Request},
    overlap_extension_pcr::{PrimerPair, Request as Preparation},
    restriction_inspection::InspectionOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Capsule, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
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
    let report = inspect::inspect(&source, None, &request).unwrap();
    report
        .validate_binding(&inspect::inspection_revision(&source, None, &request).unwrap())
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
    let report = inspect::inspect(&repeated, None, &request()).unwrap();
    assert!(report.incomplete);
    assert!(!report.selectable());
    assert!(report.amplicon.is_none());
    assert!(report.extension.is_none());
    assert!(report.inspection.is_none());
    assert!(report.pcr_stages[0].geometries.len() > 1);
    let mut request = request();
    request.preparation.outer = pair("T", "", "A", "");
    let report = inspect::inspect(&source(), None, &request).unwrap();
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
    let report = inspect::inspect(&source(), None, &request).unwrap();
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
    let capped = inspect::inspect(&source(), None, &request).unwrap();
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
fn saved(incomplete: bool) -> (Molecule, Capsule) {
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
    (source, origin)
}
#[test]
fn exact_saved_origin_is_replayed_and_its_original_incompleteness_and_target_are_retained() {
    for incomplete in [false, true] {
        let (source, origin) = saved(incomplete);
        let pin = serde_json::to_value(&origin).unwrap();
        let request = request();
        assert!(chemistry::has_unmapped_physical_evidence(&source));
        assert!(inspect::source_preflight(&source, None, &request).is_err());
        let report = inspect::inspect(&source, Some(&origin), &request).unwrap();
        assert_eq!(report.incomplete, incomplete);
        assert!(report.selectable());
        assert_eq!(
            report.origin_status,
            Some(if incomplete {
                Status::Incomplete
            } else {
                Status::ConsistentUnderAssumptions
            })
        );
        let evidence = report.origin_evidence.as_ref().unwrap();
        assert_eq!(evidence.selection, origin.selection);
        assert_eq!(evidence.target, origin.target_event);
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
        assert!(inspect::source_preflight(&altered, Some(&origin), &request).is_err());
    }
}
#[test]
fn documentary_stage_catalog_context_and_source_binding_tampering_is_rejected() {
    let report = inspect::inspect(&source(), None, &request()).unwrap();
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
