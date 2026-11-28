//! Library finishing is an explicit native route, independent of target intent.
use collagene_core::{
    chemistry::Strand,
    cloning::{
        assembly::{
            intent::Selector,
            recipe::{OverlapFinishing, Preparation, RunPreparation},
        },
        overlap::{
            annealing as homology_annealing, assembly as homology_assembly,
            profiles::{OverlapMethod, ProfileOptions},
            resection::Cuts,
        },
    },
    documents::workspace::{self, Workspace},
    material::{
        derivation as material,
        lineage::{self as run_lineage, EventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, MaterialOrigin},
    },
    planning::{
        combinatorial::{self as design, batch, export, release},
        plan::Operation,
    },
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const LEFT: &str = "AGTCCGATGCTACGTAGGAC";
const RIGHT: &str = "TTCAGACCTGTAGATCGTCA";
const CIRCLE: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn five_prime() -> OverlapFinishing {
    OverlapFinishing::HomologyAssembly {
        policy: homology_assembly::Options {
            overlaps: vec![8, 8],
            topology: Topology::Circular,
            min_paired_bases: 4,
            max_comparisons: 100_000,
            max_candidates: 20,
            ..Default::default()
        },
    }
}
fn three_prime() -> OverlapFinishing {
    OverlapFinishing::HomologyAnnealing {
        policy: homology_annealing::Options {
            overlaps: vec![8, 8],
            topology: Topology::Circular,
            cuts: vec![Cuts { left: 8, right: 8 }, Cuts { left: 8, right: 8 }],
            min_paired_bases: 4,
            max_comparisons: 100_000,
            max_candidates: 20,
            planned_rescue: None,
        },
    }
}
fn fixture(saved: bool) -> (Workspace, design::Design, design::SelectedConstruct) {
    let mut w = Workspace::from_documents(
        "Explicit native finishing",
        vec![
            Molecule::dna("Original left", LEFT, Topology::Linear).unwrap(),
            Molecule::dna("Original right", RIGHT, Topology::Linear).unwrap(),
            Molecule::dna(
                "Independent wrong target",
                &format!("C{}", &CIRCLE[1..]),
                Topology::Circular,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let left_f = Primer::new("Left F", "AGTCCGAT", "").unwrap();
    let left_r = Primer::new("Left R with downstream homology", "GTCCTACG", "GGTCTGAA").unwrap();
    let mut left_source = event(&w, 0);
    let mut bindings = vec![];
    let left_preparation = if saved {
        let reverse = Primer::new("Source R", "GTCCTACG", "").unwrap();
        let source = &w.revisions[&left_source.revision];
        let previews =
            verification::pcr::simulate_pcr_documents(source, &left_f, &reverse, Strand::Top)
                .unwrap();
        assert_eq!(previews.len(), 1);
        let installed = material::prepare_install(
            &w,
            &material::Request {
                sources: vec![collagene_core::material::derivation::BoundSource {
                    role: collagene_core::material::derivation::SourceRole::Template,
                    event: left_source,
                    origin: material::SourceOrigin::OriginalReference,
                }],
                derivation: material::Derivation::ExactPcr {
                    forward: left_f.clone(),
                    reverse,
                    template_strand: Strand::Top,
                    selection: material::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: previews[0].molecule.revision(),
                        geometry: previews[0].geometry.clone(),
                    },
                },
                name: "Retained actual left PCR".into(),
                save_primers: false,
            },
        )
        .unwrap();
        left_source = installed.receipt.product;
        w = installed.workspace;
        let origin = material::capsule(&w, &left_source).unwrap();
        bindings.push(design::LibraryMaterialBinding {
            event: left_source.clone(),
            origin: MaterialOrigin::Derived(origin.clone()),
        });
        Preparation::FromDerived {
            origin,
            preparation: RunPreparation::ExactPcr {
                forward: left_f,
                reverse: left_r,
                template_strand: Strand::Top,
            },
        }
    } else {
        Preparation::ExactPcr {
            forward: left_f,
            reverse: left_r,
            template_strand: Strand::Top,
        }
    };
    let d = design::Design {
        id: "finishing-library".into(),
        name: "Independent finishing choices".into(),
        circular: true,
        material_bindings: bindings,
        method: design::Method::Overlap {
            overlaps: vec![8, 8],
            options: ProfileOptions {
                scan_repeated_overlaps: false,
                ..Default::default()
            },
            finishing: five_prime(),
        },
        bins: vec![
            design::Bin {
                id: "left".into(),
                label: "Left".into(),
                choices: vec![design::Choice {
                    anneal_bottom: None,
                    id: "left-choice".into(),
                    label: "Left native source".into(),
                    source: left_source,
                    selector: Selector::ReferenceRange {
                        start: 0,
                        end: 20,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    generation: None,
                    treatment: None,
                    preparation: Some(left_preparation),
                }],
            },
            design::Bin {
                id: "right".into(),
                label: "Right".into(),
                choices: vec![design::Choice {
                    anneal_bottom: None,
                    id: "right-choice".into(),
                    label: "Right original".into(),
                    source: event(&w, 1),
                    selector: Selector::ReferenceRange {
                        start: 0,
                        end: 20,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    generation: None,
                    treatment: None,
                    preparation: Some(Preparation::ExactPcr {
                        forward: Primer::new("Right F", "TTCAGACC", "").unwrap(),
                        reverse: Primer::new(
                            "Right R with closure homology",
                            "TGACGATC",
                            "ATCGGACT",
                        )
                        .unwrap(),
                        template_strand: Strand::Top,
                    }),
                }],
            },
        ],
    };
    let row = design::SelectedConstruct {
        id: "native-circle".into(),
        name: "Independent literal forty".into(),
        choices: vec![
            design::ChoiceSelection {
                bin_id: "left".into(),
                choice_id: "left-choice".into(),
            },
            design::ChoiceSelection {
                bin_id: "right".into(),
                choice_id: "right-choice".into(),
            },
        ],
        target: design::TargetDeclaration::OriginalRegions {
            junction_edits: vec![],
            substitutions: vec![],
            edits: vec![],
        },
        overrides: None,
    };
    (w, d, row)
}
fn finishing(d: &mut design::Design) -> &mut OverlapFinishing {
    let design::Method::Overlap { finishing, .. } = &mut d.method else {
        unreachable!()
    };
    finishing
}
fn release_row(
    w: Workspace,
    d: design::Design,
    row: design::SelectedConstruct,
    incomplete: bool,
) -> release::Commit {
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows: vec![batch::RowRequest::Ready(row)],
        },
    )
    .unwrap();
    let ctx = session.summary().context_revision;
    session.advance(&ctx, 0).unwrap();
    let captured = session.into_capture(&ctx).unwrap();
    let restarted =
        workspace::load_workspace(&workspace::save_workspace(&captured.workspace).unwrap())
            .unwrap();
    let replay = batch::replay(&restarted, &captured.attachment_hash, 0).unwrap();
    assert_eq!(
        replay.current.status,
        Some(if incomplete {
            Status::Incomplete
        } else {
            Status::ConsistentUnderAssumptions
        })
    );
    let mut released = release::ReleaseSession::begin(
        restarted,
        release::Request {
            batch_hash: captured.attachment_hash,
            rows: vec![release::Selection {
                construct_id: "native-circle".into(),
                options: cloning_outputs::Options {
                    name: "Released explicit native route".into(),
                    show_prepared: false,
                    save_primers: true,
                    allow_incomplete: incomplete,
                },
            }],
        },
    )
    .unwrap();
    let ctx = released.summary().context_revision;
    released.advance(&ctx, 0).unwrap();
    released.into_commit(&ctx).unwrap()
}
#[test]
fn five_prime_saved_source_library_releases_native_material_and_exports_actual_protocol_after_restart()
 {
    let (w, d, row) = fixture(true);
    let before = workspace::save_workspace(&w).unwrap();
    let physical = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(physical.resolved.target.sequence.as_str(), CIRCLE);
    assert_eq!(
        physical.verification().status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        physical.verification().findings
    );
    assert!(matches!(
        physical.plan().steps.last(),
        Some(Operation::HomologyAssembly { .. })
    ));
    assert_eq!(
        physical
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        CIRCLE
    );
    let report = physical
        .verification()
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(r) => Some(r),
            _ => None,
        })
        .unwrap();
    assert!(report.product.as_ref().unwrap().nicks.is_empty());
    assert!(
        report
            .junctions
            .iter()
            .all(|j| j.top_sealed && j.bottom_sealed)
    );
    let mut reference = d.clone();
    *finishing(&mut reference) = OverlapFinishing::ReferenceOnly;
    let reference = design::evaluate(&reference, &row, &w).unwrap();
    assert!(matches!(
        reference.plan().steps.last(),
        Some(Operation::ProfiledOverlap { .. })
    ));
    assert_eq!(
        reference.resolved.target.revision(),
        physical.resolved.target.revision()
    );
    assert_ne!(
        reference.resolved.declaration_revision,
        physical.resolved.declaration_revision
    );
    assert_eq!(
        physical.resolved.declared_work() - reference.resolved.declared_work(),
        100_000
    );
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let heads = w.documents.clone();
    let committed = release_row(w, d, row, false);
    for head in heads {
        assert_eq!(
            committed
                .workspace
                .documents
                .iter()
                .find(|d| d.id == head.id)
                .unwrap()
                .head,
            head.head
        );
    }
    let manifest = release::load(&committed.workspace, &committed.manifest_hash).unwrap();
    let output = &manifest.rows[0].output;
    let replay = run_lineage::replay_run(&committed.workspace, &output.run.run_id).unwrap();
    assert!(replay.products_match);
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    let product = &output.construct;
    let capsule = collagene_core::material::capsule::capsule(
        &committed.workspace,
        &run_material::Selection {
            run_id: output.run.run_id.clone(),
            step_index: replay.current_report.products.len() - 1,
            event: product.clone(),
        },
    )
    .unwrap();
    let fresh = run_material::MaterialOrigin::Run(capsule)
        .replay(&committed.workspace.revisions[&product.revision])
        .unwrap();
    assert!(fresh.tracked().duplex.nicks.is_empty());
    assert_eq!(fresh.tracked().duplex.top.as_str(), CIRCLE);
    let exported = export::export(
        &committed.workspace,
        &committed.manifest_hash,
        &export::Options {
            genbank: false,
            protocol_html: true,
            protocol_pdf: true,
            primer_csv: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(exported.files.iter().any(|f| f.path.ends_with(".html")));
    assert!(exported.files.iter().any(|f| f.path.ends_with(".pdf")));
    assert!(exported.files.iter().any(|f| f.path.ends_with(".csv")));
}
#[test]
fn three_prime_library_keeps_real_segments_and_nicks_and_planned_rescue_never_creates_material() {
    let (w, mut d, row) = fixture(false);
    *finishing(&mut d) = three_prime();
    let preview = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(preview.resolved.target.sequence.as_str(), CIRCLE);
    assert_eq!(
        preview.verification().status,
        Status::Incomplete,
        "{:?}",
        preview.verification().findings
    );
    assert!(matches!(
        preview.plan().steps.last(),
        Some(Operation::HomologyAnnealing { .. })
    ));
    let report = preview
        .verification()
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(r) => Some(r),
            _ => None,
        })
        .unwrap();
    assert_eq!(report.reference.as_str(), CIRCLE);
    assert_eq!(report.segments.len(), 4);
    assert_eq!(report.paired_diagnostic.as_ref().unwrap().nicks.len(), 4);
    let committed = release_row(w.clone(), d.clone(), row.clone(), true);
    let manifest = release::load(&committed.workspace, &committed.manifest_hash).unwrap();
    let output = &manifest.rows[0].output;
    let replay = run_lineage::replay_run(&committed.workspace, &output.run.run_id).unwrap();
    assert!(replay.products_match);
    assert_eq!(replay.current_report.status, Status::Incomplete);
    assert!(
        collagene_core::material::capsule::capsule(
            &committed.workspace,
            &run_material::Selection {
                run_id: output.run.run_id.clone(),
                step_index: replay.current_report.products.len() - 1,
                event: output.construct.clone()
            }
        )
        .is_err(),
        "retained unsealed contour cannot become reusable duplex stock"
    );
    let design::Method::Overlap { options, .. } = &mut d.method else {
        unreachable!()
    };
    options.method = OverlapMethod::InFusion;
    assert!(
        d.validate_construct(&row).is_err(),
        "real InFusion-associated annealing requires explicit planned rescue"
    );
    let design::Method::Overlap {
        finishing: OverlapFinishing::HomologyAnnealing { policy },
        ..
    } = &mut d.method
    else {
        unreachable!()
    };
    policy.planned_rescue = Some(homology_annealing::PlannedRescue::EcoliTransformation);
    let rescued = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(rescued.verification().status, Status::Incomplete);
    let report = rescued
        .verification()
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(r) => Some(r),
            _ => None,
        })
        .unwrap();
    assert_eq!(report.paired_diagnostic.as_ref().unwrap().nicks.len(), 4);
}
#[test]
fn mismatched_policy_is_row_local_and_wrong_independent_target_is_never_replaced_by_product() {
    let (w, d, row) = fixture(false);
    let mut wrong = row.clone();
    wrong.id = "wrong-target".into();
    wrong.target = design::TargetDeclaration::Existing(event(&w, 2));
    let evaluated = design::evaluate(&d, &wrong, &w).unwrap();
    assert_eq!(evaluated.verification().status, Status::Invalid);
    assert_eq!(
        evaluated
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        CIRCLE
    );
    let mut malformed = row.clone();
    malformed.id = "wrong-policy".into();
    let mut policy = five_prime();
    let OverlapFinishing::HomologyAssembly { policy } = &mut policy else {
        unreachable!()
    };
    policy.topology = Topology::Linear;
    malformed.overrides = Some(design::ConstructOverrides {
        preparations: vec![],
        method: Some(design::Method::Overlap {
            overlaps: vec![8, 8],
            options: ProfileOptions::default(),
            finishing: OverlapFinishing::HomologyAssembly {
                policy: policy.clone(),
            },
        }),
        circular: None,
    });
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows: vec![
                batch::RowRequest::Ready(malformed),
                batch::RowRequest::Ready(row),
            ],
        },
    )
    .unwrap();
    let ctx = session.summary().context_revision;
    assert_eq!(
        session.advance(&ctx, 0).unwrap().row.unwrap().state,
        batch::RowState::Error
    );
    assert_eq!(
        session.advance(&ctx, 1).unwrap().row.unwrap().status,
        Some(Status::ConsistentUnderAssumptions)
    );
}
#[test]
fn finishing_is_required_and_physical_search_reservations_bound_the_existing_batch() {
    let (w, mut d, row) = fixture(false);
    let mut json = serde_json::to_value(&d).unwrap();
    json["method"]["Overlap"]
        .as_object_mut()
        .unwrap()
        .remove("finishing");
    assert!(serde_json::from_value::<design::Design>(json).is_err());
    let OverlapFinishing::HomologyAssembly { policy } = finishing(&mut d) else {
        unreachable!()
    };
    policy.max_comparisons = 100_000_000;
    let rows = (0..21)
        .map(|i| {
            let mut r = row.clone();
            r.id = format!("bounded-{i}");
            batch::RowRequest::Ready(r)
        })
        .collect();
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows,
        },
    )
    .unwrap();
    let ctx = session.summary().context_revision;
    for i in 0..21 {
        let r = session.advance(&ctx, i).unwrap().row.unwrap();
        if i < 19 {
            assert_eq!(r.status, Some(Status::ConsistentUnderAssumptions));
        } else {
            assert_eq!(r.state, batch::RowState::Error);
            assert!(
                r.error
                    .unwrap()
                    .contains("Not computed: aggregate batch work budget")
            );
        }
    }
}
