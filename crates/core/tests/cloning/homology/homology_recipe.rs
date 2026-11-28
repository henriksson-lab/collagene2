//! Independent literal products through prepare → assemble → save → reuse.
use collagene_core::{
    chemistry::{EndRepairMethod, FivePrime, ThreePrime},
    cloning::{
        assembly::recipe::{self, FragmentSpec, Preparation, RunPreparation},
        ligation::recipe as ligation_recipe,
        overlap::assembly::{Options, OuterEndPolicy},
        restriction::catalog as catalog_digest,
    },
    documents::workspace::Workspace,
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::Selection,
    },
    planning::plan::{Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Feature, Location, Molecule, Topology};
const LEFT: &str = "AGTCCGATGCTACGTAGGAC";
const RIGHT: &str = "TTCAGACCTGTAGATCGTCA";
const CIRCLE: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const REPAIRED: &str = "AATTCAGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCAGAATT";
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture(circular: bool) -> (Vec<PlanInput>, recipe::Recipe, Options) {
    let mut left = Molecule::dna("Original left template", LEFT, Topology::Linear).unwrap();
    left.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 8),
        qualifiers: vec![("label".into(), Some("Original left calls".into()))],
        needs_review: false,
    });
    let right = Molecule::dna("Original right template", RIGHT, Topology::Linear).unwrap();
    let inputs = vec![input("left", left), input("right", right)];
    let fragments = vec![
        FragmentSpec {
            generation: None,
            role: "left-role".into(),
            input_id: "left".into(),
            reverse: false,
            preparation: Preparation::ExactPcr {
                forward: Primer::new("Left F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Left R + downstream homology", "GTCCTACG", "GGTCTGAA")
                    .unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        },
        FragmentSpec {
            generation: None,
            role: "right-role".into(),
            input_id: "right".into(),
            reverse: false,
            preparation: Preparation::ExactPcr {
                forward: Primer::new("Right F", "TTCAGACC", "").unwrap(),
                reverse: Primer::new(
                    "Right R + circle closure",
                    "TGACGATC",
                    if circular { "ATCGGACT" } else { "" },
                )
                .unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        },
    ];
    let overlaps = if circular { vec![8, 8] } else { vec![8] };
    let mut profile = collagene_core::cloning::overlap::profiles::ProfileOptions::default();
    profile.scan_repeated_overlaps = false;
    let request = recipe::Recipe {
        name: "Declared physical homology route".into(),
        fragments,
        overlaps: overlaps.clone(),
        circular,
        options: profile,
    };
    let policy = Options {
        overlaps,
        topology: if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
        min_paired_bases: 4,
        max_comparisons: 100000,
        max_candidates: 20,
        ..Default::default()
    };
    (inputs, request, policy)
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn pcr_five_prime_oh_fragments_form_a_native_circle_then_closed_origin_supplies_pcr_and_digest() {
    let (inputs, request, policy) = fixture(true);
    let target = Molecule::dna("Independent literal circle40", CIRCLE, Topology::Circular).unwrap();
    let preview = recipe::preview_homology(&inputs, &target, &request, &policy).unwrap();
    preview
        .validate_binding(
            &recipe::homology_preview_revision(&inputs, &target, &request, &policy).unwrap(),
        )
        .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.assembly_step, 2);
    assert!(matches!(
        preview.compiled.plan.steps[2],
        Operation::HomologyAssembly { .. }
    ));
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        CIRCLE
    );
    let processing = preview
        .verification
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(r) => Some(r),
            _ => None,
        })
        .unwrap();
    assert_eq!(processing.junctions.len(), 2);
    assert!(processing.product.as_ref().unwrap().nicks.is_empty());
    for j in &processing.junctions {
        assert_eq!(j.top_five_prime, FivePrime::Phosphate);
        assert_eq!(j.bottom_five_prime, FivePrime::Phosphate);
        assert_eq!(j.top_three_prime, ThreePrime::Hydroxyl);
        assert!(j.top_sealed && j.bottom_sealed);
    }
    let w = Workspace::from_documents(
        "Original events and independent target",
        vec![
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target,
        ],
        2,
    )
    .unwrap();
    let target_pin = binding(&w, 2);
    let commit = cloning_outputs::finalize(
        &w,
        &preview.compiled.plan,
        target_pin.clone(),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: binding(&w, 0),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: binding(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Actual assembled circle".into(),
            show_prepared: false,
            save_primers: true,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
    let source = commit.workspace.revisions[&event.revision].clone();
    let origin = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 2,
            event: event.clone(),
        },
    )
    .unwrap();
    assert_eq!(origin.target_event, target_pin);
    let before = serde_json::to_value(&origin).unwrap();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection = catalog_digest::select(&catalog, &[id]).unwrap();
    let downstream = ligation_recipe::Recipe {
        name: "Reuse actual circle, PCR and cut real amplicon".into(),
        circular: false,
        allow_nicked: false,
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "actual-circle".into(),
                input_id: "circle".into(),
                reverse: false,
                preparation: Preparation::FromRun {
                    origin: origin.clone(),
                    preparation: RunPreparation::PcrCatalogDigest {
                        forward: Primer::new("Stock F + EcoRI", "AGTCCGAT", "CCGAATTC").unwrap(),
                        reverse: Primer::new("Stock R + EcoRI", "TGACGATC", "GGGAATTC").unwrap(),
                        selection,
                        fragment_index: 1,
                        repair: Some(EndRepairMethod::FillFivePrimeOverhangs),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    },
                },
            },
            treatment: None,
        }],
    };
    let independent = Molecule::dna("Independent repaired50", REPAIRED, Topology::Linear).unwrap();
    let p = ligation_recipe::preview(
        &[input("circle", source.clone())],
        &independent,
        &downstream,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(
        p.verification.products.last().unwrap().sequence.as_str(),
        REPAIRED
    );
    assert!(matches!(
        p.compiled.plan.steps[0],
        Operation::ReconstructRunProduct { .. }
    ));
    assert_eq!(serde_json::to_value(&origin).unwrap(), before);
    // Native export/import and retained replay include the physical operation.
    let archive = collagene_core::documents::workspace::save_workspace(&commit.workspace).unwrap();
    let restored = collagene_core::documents::workspace::load_workspace(&archive).unwrap();
    assert!(
        collagene_core::material::capsule::capsule(
            &restored,
            &Selection {
                run_id: commit.receipt.run.run_id,
                step_index: 2,
                event
            }
        )
        .is_ok()
    );
    let html = verification::output::export_protocol_html(
        &preview.compiled.plan,
        &w.revisions[&target_pin.revision],
    );
    assert!(html.contains("Physical homology assembly"));
    assert!(html.contains("5′"));
    let pdf = collagene_core::output::pdf::export_protocol_pdf(
        &preview.compiled.plan,
        &w.revisions[&target_pin.revision],
    )
    .unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, Status::ConsistentUnderAssumptions);
}
#[test]
fn reference_only_and_imported_product_summaries_cannot_acquire_homology_material_authority() {
    let (inputs, request, policy) = fixture(true);
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let reference = recipe::preview(&inputs, &target, &request).unwrap();
    assert_eq!(
        reference.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let mut plan = reference.compiled.plan.clone();
    plan.steps.push(Operation::Pcr {
        template: Source::Step(2),
        forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    });
    let linear = Molecule::dna("Independent PCR40", CIRCLE, Topology::Linear).unwrap();
    plan.target_revision = linear.revision();
    let unavailable = verification::engine::verify(&plan, &linear);
    assert_eq!(unavailable.status, Status::Incomplete);
    assert_eq!(
        unavailable.products.len(),
        3,
        "no downstream PCR product may be invented"
    );
    assert!(
        unavailable
            .findings
            .iter()
            .any(|f| f.contains("actual physical source state"))
    );
    let physical = recipe::preview_homology(&inputs, &target, &request, &policy).unwrap();
    let imported = physical.verification.products.last().unwrap().clone();
    let forged = collagene_core::planning::plan::DesignPlan {
        name: "Metadata is not originating material".into(),
        target_revision: linear.revision(),
        inputs: vec![input("copied", imported)],
        steps: vec![Operation::Pcr {
            template: Source::Input("copied".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let imported = verification::engine::verify(&forged, &linear);
    assert_eq!(imported.status, Status::Incomplete);
    assert!(
        imported.products.is_empty(),
        "producer metadata must not supply a PCR template"
    );
    assert!(
        imported
            .findings
            .iter()
            .any(|f| f.contains("fresh native origin"))
    );
    let mut mismatched = policy;
    mismatched.overlaps = vec![7, 7];
    assert!(recipe::homology_preflight(&inputs, &target, &request, &mismatched).is_err());
    let mut incompatible = request.clone();
    incompatible.options.method =
        collagene_core::cloning::overlap::profiles::OverlapMethod::InFusion;
    let (_, _, physical_policy) = fixture(true);
    assert!(recipe::homology_preflight(&inputs, &target, &incompatible, &physical_policy).is_err());
    let mut direct = physical.compiled.plan;
    if let Operation::HomologyAssembly { options, .. } = direct.steps.last_mut().unwrap() {
        options.method = collagene_core::cloning::overlap::profiles::OverlapMethod::InFusion;
    } else {
        panic!("physical assembly");
    }
    assert_eq!(
        verification::engine::verify(&direct, &target).status,
        Status::Invalid
    );
}
#[test]
fn linear_outer_resection_projects_only_actual_top_calls_and_preserves_real_overhangs() {
    let (inputs, request, mut policy) = fixture(false);
    policy.outer_end_policy = OuterEndPolicy::AllAccessibleEnds { left: 2, right: 3 };
    let target = Molecule::dna(
        "Independent actual top38",
        "TCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let p = recipe::preview_homology(&inputs, &target, &request, &policy).unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    let product = p.verification.products.last().unwrap();
    assert_eq!(product.sequence.as_str(), target.sequence.as_str());
    let clipped = product
        .features
        .iter()
        .find(|f| {
            f.qualifiers
                .iter()
                .any(|(k, v)| k == "label" && v.as_deref() == Some("Original left calls"))
        })
        .unwrap();
    assert!(clipped.needs_review);
    assert!(
        matches!(
            clipped.location,
            Location::Range {
                start: 0,
                end: 6,
                // The shared projector conservatively marks the entire
                // clipped feature partial (annotations::mark_partial).
                before: true,
                after: true
            }
        ),
        "{:?}",
        clipped.location
    );
    assert!(
        clipped.location.extract(&product.sequence).is_err(),
        "a fuzzy clipped feature is not a complete extracted feature"
    );
    assert_eq!(&product.sequence.as_str()[..6], "TCCGAT");
    let g = p
        .verification
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(r) => Some(r),
            _ => None,
        })
        .unwrap();
    let d = g.product.as_ref().unwrap();
    assert_eq!(d.top.len(), 38);
    assert_eq!(d.bottom.len(), 37);
    assert_eq!(d.bottom_start, -2);
    let ends = d.ends().unwrap();
    assert_eq!(ends.0.overhang.len(), 2);
    assert_eq!(ends.1.overhang.len(), 3);
}

#[test]
fn unsealed_diagnostic_strands_are_retained_but_cannot_become_saved_material() {
    use collagene_core::chemistry::Terminus;
    for three_prime in [ThreePrime::Unknown, ThreePrime::Blocked] {
        let (_, mut request, policy) = fixture(true);
        let inputs = vec![
            input(
                "left",
                Molecule::dna(
                    "Declared left stock",
                    "AGTCCGATGCTACGTAGGACTTCAGACC",
                    Topology::Linear,
                )
                .unwrap(),
            ),
            input(
                "right",
                Molecule::dna(
                    "Declared right stock",
                    "TTCAGACCTGTAGATCGTCAAGTCCGAT",
                    Topology::Linear,
                )
                .unwrap(),
            ),
        ];
        for row in &mut request.fragments {
            let end = Terminus {
                five_prime: FivePrime::Hydroxyl,
                three_prime,
            };
            row.preparation = Preparation::DeclaredDuplex {
                left: end,
                right: end,
            };
        }
        let target =
            Molecule::dna("Independent diagnostic circle", CIRCLE, Topology::Circular).unwrap();
        let p = recipe::preview_homology(&inputs, &target, &request, &policy).unwrap();
        assert_eq!(
            p.verification.status,
            Status::Incomplete,
            "{:?}",
            p.verification.findings
        );
        assert_eq!(
            p.verification.products.last().unwrap().sequence.as_str(),
            CIRCLE
        );
        assert!(
            p.verification
                .products
                .last()
                .unwrap()
                .metadata
                .contains_key(collagene_core::cloning::recombination::REFERENCE_ONLY_KEY)
        );
        let g = p
            .verification
            .physical_processing_reports
            .iter()
            .find_map(|s| match &s.evidence {
                verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(g) => Some(g),
                _ => None,
            })
            .unwrap();
        assert_eq!(g.product.as_ref().unwrap().nicks.len(), 4);
        assert!(
            g.junctions
                .iter()
                .all(|j| !j.top_sealed && !j.bottom_sealed && j.top_three_prime == three_prime)
        );
        let w = Workspace::from_documents(
            "Diagnostic source declarations",
            vec![
                inputs[0].molecule.clone(),
                inputs[1].molecule.clone(),
                target,
            ],
            2,
        )
        .unwrap();
        let commit = cloning_outputs::finalize(
            &w,
            &p.compiled.plan,
            binding(&w, 2),
            vec![
                InputEventBinding {
                    input_id: "left".into(),
                    event: binding(&w, 0),
                },
                InputEventBinding {
                    input_id: "right".into(),
                    event: binding(&w, 1),
                },
            ],
            &cloning_outputs::Options {
                name: "Retained incomplete diagnostic".into(),
                show_prepared: false,
                save_primers: false,
                allow_incomplete: true,
            },
        )
        .unwrap();
        assert!(
            collagene_core::material::capsule::capsule(
                &commit.workspace,
                &Selection {
                    run_id: commit.receipt.run.run_id,
                    step_index: 2,
                    event: commit.receipt.construct
                }
            )
            .is_err()
        );
        let mut downstream = p.compiled.plan;
        let expected =
            Molecule::dna("Independent downstream PCR", CIRCLE, Topology::Linear).unwrap();
        downstream.target_revision = expected.revision();
        downstream.steps.push(Operation::Pcr {
            template: Source::Step(2),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        });
        let unavailable = verification::engine::verify(&downstream, &expected);
        assert_eq!(unavailable.status, Status::Incomplete);
        assert_eq!(
            unavailable.products.len(),
            3,
            "diagnostic strands may not create downstream PCR material"
        );
        assert!(
            unavailable
                .findings
                .iter()
                .any(|f| f.contains("Reference-only ancestry"))
        );
    }
}
