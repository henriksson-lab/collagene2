//! Saved material is reconstructed before the ordinary exact-PCR/digest engine.
use collagene_core::{
    chemistry::{EndTreatment, FivePrime, ThreePrime},
    cloning::{pcr::digest as p, restriction::catalog as catalog_digest},
    documents::workspace::Workspace,
    material::{
        capsule::Capsule,
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs, run as run_material,
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
    verification,
};
use collagene_model::{Molecule, Topology};
fn saved(
    originals: Vec<Molecule>,
    target: Molecule,
    steps: Vec<Operation>,
) -> (Molecule, Capsule, Workspace) {
    let inputs = originals
        .iter()
        .enumerate()
        .map(|(i, m)| PlanInput {
            id: format!("src{i}"),
            revision: m.revision(),
            molecule: m.clone(),
        })
        .collect();
    let count = originals.len();
    let mut docs = originals;
    docs.push(target.clone());
    let workspace = Workspace::from_documents("Origin", docs, 0).unwrap();
    let event = |i: usize| {
        let d = &workspace.documents[i];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    };
    let step_index = steps.len() - 1;
    let plan = DesignPlan {
        standard: None,
        name: "Original reactions".into(),
        target_revision: target.revision(),
        inputs,
        steps,
    };
    let commit = cloning_outputs::finalize(
        &workspace,
        &plan,
        event(count),
        (0..count)
            .map(|i| InputEventBinding {
                input_id: format!("src{i}"),
                event: event(i),
            })
            .collect(),
        &cloning_outputs::Options {
            name: "Named saved product".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let source = commit
        .workspace
        .head(&commit.receipt.construct.document_id)
        .unwrap()
        .clone();
    let capsule = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &run_material::Selection {
            run_id: commit.receipt.run.run_id,
            step_index,
            event: commit.receipt.construct,
        },
    )
    .unwrap();
    (source, capsule, commit.workspace)
}
fn pcr(template: Source, forward: &str, reverse: &str) -> Operation {
    Operation::Pcr {
        template,
        forward: Primer::new("F", forward, "").unwrap(),
        reverse: Primer::new("R", reverse, "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn fixture() -> (Molecule, Capsule, p::Request, Workspace) {
    let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent original target",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let (source, origin, workspace) = saved(
        vec![source],
        target,
        vec![pcr(Source::Input("src0".into()), "ACGACCCC", "GTTACCCC")],
    );
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let request = p::Request {
        template_strand: collagene_core::chemistry::Strand::Top,
        source_role: "saved-template".into(),
        forward: Primer::new("Downstream F", "ACGACCCC", "CCGAATTC").unwrap(),
        reverse: Primer::new("Downstream R", "GTTACCCC", "GGGAATTC").unwrap(),
        selection: catalog_digest::select(&catalog, &[enzyme]).unwrap(),
        options: Default::default(),
    };
    (source, origin, request, workspace)
}
#[test]
fn fresh_saved_pcr_tail_digest_retains_outer_hydroxyl_and_internal_phosphate() {
    let (source, origin, request, workspace) = fixture();
    let before = serde_json::to_vec(&source).unwrap();
    let result = p::inspect_from_material(
        &source,
        &run_material::MaterialOrigin::Run(origin.clone()),
        &request,
        &workspace,
    )
    .unwrap();
    result
        .validate_binding(
            &p::inspection_revision_from_material(
                &source,
                &run_material::MaterialOrigin::Run(origin.clone()),
                &request,
                &workspace,
            )
            .unwrap(),
        )
        .unwrap();
    assert_eq!(result.origin_revision.len(), 64);
    assert_eq!(
        result.preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCACGACCCCGGGGTAACGAATTCCC"
    );
    let inspection = result.preview.inspection.as_ref().unwrap();
    assert!(inspection.scan_complete && inspection.products_determined);
    assert_eq!(inspection.fragments.len(), 3);
    let ends = inspection
        .fragments
        .iter()
        .map(|f| f.duplex.ends.as_ref().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(ends[0].0.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(ends[0].1.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[1].0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[1].1.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[2].0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[2].1.chemistry.five_prime, FivePrime::Hydroxyl);
    for (left, right) in ends {
        assert_eq!(left.chemistry.three_prime, ThreePrime::Hydroxyl);
        assert_eq!(right.chemistry.three_prime, ThreePrime::Hydroxyl);
    }
    assert_eq!(serde_json::to_vec(&source).unwrap(), before);
}
#[test]
fn competing_pcr_geometries_withhold_saved_material_digest_choices() {
    let (source, origin, mut request, workspace) = fixture();
    request.forward = Primer::new("Ambiguous F", "C", "").unwrap();
    request.reverse = Primer::new("Ambiguous R", "C", "").unwrap();
    let result = p::inspect_from_material(
        &source,
        &run_material::MaterialOrigin::Run(origin.clone()),
        &request,
        &workspace,
    )
    .unwrap();
    assert!(result.preview.pcr_geometries.len() > 1);
    assert!(result.preview.amplicon.is_none() && result.preview.inspection.is_none());
}
#[test]
fn changed_origin_request_and_source_cannot_reuse_saved_inspection_context() {
    let (source, origin, request, workspace) = fixture();
    let result = p::inspect_from_material(
        &source,
        &run_material::MaterialOrigin::Run(origin.clone()),
        &request,
        &workspace,
    )
    .unwrap();
    let mut other = origin.clone();
    other.products[0].revision = "0".repeat(64);
    assert!(
        p::inspection_revision_from_material(
            &source,
            &run_material::MaterialOrigin::Run(other.clone()),
            &request,
            &workspace
        )
        .is_err()
    );
    assert!(
        p::inspect_from_material(
            &source,
            &run_material::MaterialOrigin::Run(other.clone()),
            &request,
            &workspace
        )
        .is_err()
    );
    let mut changed = request.clone();
    changed.forward.tail = "AAGAATTC".into();
    assert!(
        result
            .validate_binding(
                &p::inspection_revision_from_material(
                    &source,
                    &run_material::MaterialOrigin::Run(origin.clone()),
                    &changed,
                    &workspace
                )
                .unwrap()
            )
            .is_err()
    );
    let mut changed_source = source.clone();
    changed_source.name.push_str(" changed");
    assert!(
        p::preflight_from_material(
            &changed_source,
            &run_material::MaterialOrigin::Run(origin.clone()),
            &request,
            &workspace
        )
        .is_err()
    );
}
#[test]
fn actual_retained_nicks_cannot_become_a_saved_pcr_template() {
    let (left, right) = ("TTGCAATCGTCGATGCTAGC", "ATCGTCGATGCTAGCGGTAC");
    let originals = [left, right]
        .into_iter()
        .enumerate()
        .map(|(i, dna)| Molecule::dna(format!("Parent{i}"), dna, Topology::Linear).unwrap())
        .collect();
    let target = Molecule::dna(
        "Independent nicked circle",
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC",
        Topology::Circular,
    )
    .unwrap();
    let (source, origin, workspace) = saved(
        originals,
        target,
        vec![
            pcr(Source::Input("src0".into()), "TTGCAATC", "GCTAGCAT"),
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            pcr(Source::Input("src1".into()), "ATCGTCGA", "GTACCGCT"),
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
    let (_, _, mut request, _) = fixture();
    request.forward = Primer::new("F", "TTGCAATC", "").unwrap();
    request.reverse = Primer::new("R", "GTACCGCT", "").unwrap();
    let err = p::inspect_from_material(
        &source,
        &run_material::MaterialOrigin::Run(origin.clone()),
        &request,
        &workspace,
    )
    .unwrap_err();
    assert!(matches!(&err, collagene_model::Error::Unsupported(_)));
    assert!(err.to_string().contains("nick"), "{err}");
    // The original native state, independently replayed, still carries both nicks.
    let plan = DesignPlan {
        standard: None,
        name: "Replay only".into(),
        target_revision: source.revision(),
        inputs: vec![PlanInput {
            id: "saved".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::ReconstructRunProduct {
            template: Source::Input("saved".into()),
            origin,
        }],
    };
    let report = verification::engine::verify(&plan, &source);
    assert_eq!(report.reconstructed_count(), 1);
    let verification::evidence::PhysicalProcessingEvidence::RunProduct(evidence) =
        &report.physical_processing_reports[0].evidence
    else {
        panic!("replay evidence")
    };
    let original = evidence.decode_report().unwrap();
    let summary: collagene_core::chemistry::DuplexSummary = serde_json::from_str(
        original
            .final_product()
            .unwrap()
            .metadata
            .get("duplex_model")
            .unwrap(),
    )
    .unwrap();
    assert_eq!(summary.nicks.len(), 2);
    assert!(
        original
            .findings
            .iter()
            .any(|f| f.message.to_ascii_lowercase().contains("nick"))
    );
}
