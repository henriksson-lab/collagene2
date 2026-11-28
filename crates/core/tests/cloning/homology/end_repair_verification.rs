use collagene_core::{
    chemistry::{self, AnnealOptions, EndRepairMethod, EndTreatment},
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    verification::{self, evidence::Status},
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn one_source(sequence: &str, enzyme: &str, target: &str) -> (DesignPlan, Molecule) {
    let source = Molecule::dna("cleaved source", sequence, Topology::Circular).unwrap();
    let target = Molecule::dna("independent repaired target", target, Topology::Circular).unwrap();
    (
        DesignPlan {
            name: "digest repair close".into(),
            target_revision: target.revision(),
            inputs: vec![input("source", source)],
            steps: vec![
                Operation::Digest {
                    template: Source::Input("source".into()),
                    enzymes: vec![enzyme.into()],
                    fragment_index: 0,
                },
                Operation::RepairEnds {
                    template: Source::Step(0),
                    method: EndRepairMethod::Blunt,
                },
                Operation::Ligate {
                    fragments: vec![OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    }],
                    circular: true,
                    allow_nicked: false,
                    options: Default::default(),
                },
            ],
        },
        target,
    )
}
#[test]
fn independent_five_prime_fill_and_three_prime_trim_targets_reconstruct_both_strands() {
    for (sequence, enzyme, expected, length) in [
        ("GAATTCAAAA", "EcoRI", "AATTCAAAAGAATT", 14),
        ("CTGCAGAAAA", "PstI", "GAAAAC", 6),
    ] {
        let (plan, target) = one_source(sequence, enzyme, expected);
        let report = verification::engine::verify(&plan, &target);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        let product = report.products.last().unwrap();
        assert_eq!(product.sequence.as_str(), expected);
        assert_eq!(product.sequence.len(), length);
        let summary: chemistry::DuplexSummary =
            serde_json::from_str(&product.metadata["duplex_model"]).unwrap();
        assert!(summary.covalently_closed_in_model);
        assert!(summary.nicks.is_empty());
        let wrong = Molecule::dna("wrong", "AAAAAA", Topology::Circular).unwrap();
        let mut wrong_plan = plan.clone();
        wrong_plan.target_revision = wrong.revision();
        assert_eq!(
            verification::engine::verify(&wrong_plan, &wrong).status,
            Status::Invalid
        );
        let html = verification::output::export_protocol_html(&plan, &target);
        assert!(html.contains("End repair"));
        assert!(
            collagene_core::output::pdf::export_protocol_pdf(&plan, &target)
                .unwrap()
                .bytes
                .starts_with(b"%PDF-")
        );
    }
}
#[test]
fn repair_preview_reconstructs_digest_and_keeps_clipped_features_as_review_evidence() {
    let mut source = Molecule::dna("annotated", "CTGCAGAAAA", Topology::Circular).unwrap();
    source.features = vec![
        Feature {
            kind: "misc_feature".into(),
            location: Location::range(6, 9),
            qualifiers: vec![("label".into(), Some("retained".into()))],
            needs_review: false,
        },
        Feature {
            kind: "CDS".into(),
            location: Location::range(1, 4),
            qualifiers: vec![("label".into(), Some("removed top-source annotation".into()))],
            needs_review: false,
        },
    ];
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Unresolved("remote:1..3".into()),
        qualifiers: vec![("label".into(), Some("unmapped ancestral annotation".into()))],
        needs_review: true,
    });
    let preview = verification::terminal::simulate_end_repair(
        &source,
        &["PstI".into()],
        0,
        EndRepairMethod::Blunt,
    )
    .unwrap();
    assert_eq!(preview.molecule.sequence.as_str(), "GAAAAC");
    assert_eq!(preview.source_revision, source.revision());
    assert_eq!(
        preview
            .molecule
            .features
            .iter()
            .find(|f| f.label() == "retained")
            .unwrap()
            .location,
        Location::range(1, 4)
    );
    assert!(
        !preview
            .molecule
            .features
            .iter()
            .any(|f| f.label() == "removed top-source annotation")
    );
    assert!(
        preview
            .molecule
            .metadata
            .keys()
            .any(|k| k.starts_with("annotation_mapping:unmapped:"))
    );
    assert!(
        verification::terminal::simulate_end_repair(&source, &[], 0, EndRepairMethod::Blunt)
            .is_err()
    );
    let mut forged = source.clone();
    forged
        .metadata
        .insert("duplex_model".into(), "{\"phosphorylated\":true}".into());
    assert!(
        verification::terminal::simulate_end_repair(
            &forged,
            &["PstI".into()],
            0,
            EndRepairMethod::Blunt
        )
        .is_err()
    );
}
fn anneal_plan(options: AnnealOptions, phosphorylate: bool) -> (DesignPlan, Molecule) {
    let top = Molecule::dna("top oligo", "ACGTCG", Topology::Linear).unwrap();
    let bottom = Molecule::dna("bottom oligo", "CGACGT", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "ACGTCG", Topology::Circular).unwrap();
    let mut steps = vec![Operation::AnnealOligos {
        top: Source::Input("top".into()),
        bottom: Source::Input("bottom".into()),
        options,
        candidate_index: 0,
    }];
    if phosphorylate {
        steps.push(Operation::TreatEnds {
            template: Source::Step(0),
            treatment: EndTreatment::Phosphorylate,
        });
    }
    steps.push(Operation::Ligate {
        fragments: vec![OrientedSource {
            source: Source::Step(steps.len() - 1),
            reverse: false,
        }],
        circular: true,
        allow_nicked: false,
        options: Default::default(),
    });
    (
        DesignPlan {
            name: "anneal phosphorylate close".into(),
            target_revision: target.revision(),
            inputs: vec![input("top", top), input("bottom", bottom)],
            steps,
        },
        target,
    )
}
#[test]
fn declared_oligos_need_explicit_phosphorylation_before_they_can_be_covalently_closed() {
    let options = AnnealOptions {
        min_paired_bases: 4,
        ..Default::default()
    };
    let (plan, target) = anneal_plan(options.clone(), true);
    let report = verification::engine::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.len(), 3);
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.contains("no covalent") || f.contains("no covalent bonds"))
    );
    let (no_phosphate, target) = anneal_plan(options, false);
    assert_eq!(
        verification::engine::verify(&no_phosphate, &target).status,
        Status::Invalid
    );
    let mut imported = plan.clone();
    imported.inputs[0]
        .molecule
        .metadata
        .insert("duplex_model".into(), "{}".into());
    imported.inputs[0].revision = imported.inputs[0].molecule.revision();
    assert_eq!(
        verification::engine::verify(&imported, &target).status,
        Status::Incomplete
    );
    let html = verification::output::export_protocol_html(&plan, &target);
    assert!(html.contains("Anneal declared"));
    assert!(html.contains("<code>ACGTCG</code>"));
    assert!(html.contains("<code>CGACGT</code>"));
    assert!(html.contains("5′ Hydroxyl; 3′ Hydroxyl"));
}
#[test]
fn anneal_annotations_use_bottom_strand_coordinates_and_report_tampering_is_rejected() {
    let top = Molecule::dna("top", "ACGTCG", Topology::Linear).unwrap();
    let mut bottom = Molecule::dna("bottom", "CGACGT", Topology::Linear).unwrap();
    bottom.features.push(Feature {
        kind: "promoter".into(),
        location: Location::range(1, 4),
        qualifiers: vec![("label".into(), Some("bottom feature".into()))],
        needs_review: false,
    });
    let options = AnnealOptions {
        min_paired_bases: 4,
        ..Default::default()
    };
    let preview =
        verification::terminal::simulate_oligo_annealing(&top, &bottom, &options).unwrap();
    assert!(preview.annealing.unique);
    let product = verification::terminal::annealed_document(&top, &bottom, &preview, 0).unwrap();
    assert_eq!(
        product.features[0].location,
        Location::Complement(Box::new(Location::range(2, 5)))
    );
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "GAC"
    );
    let mut forged = preview.clone();
    forged.annealing.candidates[0].duplex.bottom_start = 1;
    assert!(verification::terminal::annealed_document(&top, &bottom, &forged, 0).is_err());
    let mut ds = bottom.clone();
    ds.metadata
        .insert("genbank.molecule_type".into(), "ds-DNA".into());
    assert!(verification::terminal::simulate_oligo_annealing(&top, &ds, &options).is_err());
    let mut raw = bottom;
    raw.original_genbank = Some("LOCUS       explicit 6 bp ds-DNA linear\n".into());
    assert!(verification::terminal::simulate_oligo_annealing(&top, &raw, &options).is_err());
}
#[test]
fn repeated_or_capped_annealing_cannot_claim_unique_physical_verification() {
    let options = AnnealOptions {
        min_paired_bases: 2,
        ..Default::default()
    };
    let (mut plan, target) = anneal_plan(options, false);
    plan.steps.truncate(1);
    let linear = Molecule::dna("linear target", "ACGTCG", Topology::Linear).unwrap();
    plan.target_revision = linear.revision();
    // The full-match fixture may remain unique at this threshold; a repeated template has multiple offsets.
    plan.inputs[0] = input(
        "top",
        Molecule::dna("top", "AAAAAA", Topology::Linear).unwrap(),
    );
    plan.inputs[1] = input(
        "bottom",
        Molecule::dna("bottom", "TTTTTT", Topology::Linear).unwrap(),
    );
    let repeated = Molecule::dna("repeated target", "AAAAAA", Topology::Linear).unwrap();
    plan.target_revision = repeated.revision();
    let report = verification::engine::verify(&plan, &repeated);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    if let Operation::AnnealOligos { options, .. } = &mut plan.steps[0] {
        options.max_candidates = 1;
    }
    let report = verification::engine::verify(&plan, &repeated);
    assert_eq!(report.status, Status::Incomplete);
    let wrong = Molecule::dna("known wrong", "CCCCCC", Topology::Linear).unwrap();
    plan.target_revision = wrong.revision();
    assert_eq!(
        verification::engine::verify(&plan, &wrong).status,
        Status::Invalid
    );
    let _ = target;
}
