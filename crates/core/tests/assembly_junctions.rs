use collagene_core::{
    Feature, Location, Molecule, Topology, assembly_junctions as j,
    verification::{DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
};
fn molecule(name: &str, dna: &str) -> Molecule {
    Molecule::dna(name, dna, Topology::Linear).unwrap()
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn overlap(
    left: &str,
    right: &str,
    k: usize,
    target: &str,
    circular: bool,
) -> (DesignPlan, Molecule) {
    let mut target = molecule("separately supplied target", target);
    if circular {
        target.topology = Topology::Circular;
    }
    let p = DesignPlan {
        name: "literal ordered overlap".into(),
        target_revision: target.revision(),
        inputs: vec![
            input("left", molecule("left actual", left)),
            input("right", molecule("right actual", right)),
        ],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![k],
            circular,
        }],
    };
    (p, target)
}
#[test]
fn actual_shared_homology_boundary_window_relative_codon_and_genetic_code_are_explicit() {
    let (plan, target) = overlap("AAAT", "TGATTT", 1, "AAATGATTT", false);
    let mut opts = j::Options {
        flank_bases: 6,
        frames: vec![j::Frame {
            reverse: false,
            frame: 0,
        }],
        ..Default::default()
    };
    let report = j::inspect_plan(&plan, &target, &opts).unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        report.junctions[0].homology,
        vec![j::Range { start: 3, end: 4 }]
    );
    assert_eq!(report.junctions[0].boundary, 4);
    assert_eq!(report.junctions[0].window.dna, "AAATGATTT");
    assert!(report.junctions[0].window.clipped);
    assert!(
        report.complete,
        "Known linear boundaries limit the displayed flank without withholding any requested junction evidence"
    );
    let f = &report.junctions[0].frames[0];
    assert_eq!(f.translation.protein, "K*F");
    let c = f.crossing_codon.as_ref().unwrap();
    assert_eq!(
        (c.window_start, c.dna.as_str(), c.amino_acid),
        (3, "TGA", '*')
    );
    opts.genetic_code = 2;
    let mito = j::inspect_plan(&plan, &target, &opts).unwrap();
    assert_eq!(mito.junctions[0].frames[0].translation.protein, "KWF");
    let mut altered = report.clone();
    altered.junctions[0].window.dna = "CCCCCCCCC".into();
    assert!(altered.validate().is_err());
}
#[test]
fn circular_closing_homology_uses_origin_and_never_repeats_short_product() {
    let target = Molecule::dna("independent circle", "AAACCC", Topology::Circular).unwrap();
    let plan = DesignPlan {
        name: "self overlap".into(),
        target_revision: target.revision(),
        inputs: vec![input(
            "linear",
            molecule("actual closure arms", "AAACCCAAA"),
        )],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("linear".into())],
            overlaps: vec![3],
            circular: true,
        }],
    };
    let report = j::inspect_plan(&plan, &target, &j::Options::default()).unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let row = &report.junctions[0];
    assert_eq!(row.boundary, 0);
    assert_eq!(row.homology, vec![j::Range { start: 0, end: 3 }]);
    assert_eq!(
        row.window.parts,
        vec![j::Range { start: 3, end: 6 }, j::Range { start: 0, end: 3 }]
    );
    assert_eq!(row.window.dna, "CCCAAA");
    assert!(row.window.clipped);
    assert!(
        report.complete,
        "One full known circle is a complete inspection, even when the requested display flanks are longer"
    );
}
#[test]
fn ordinary_nine_hundred_base_cds_yields_only_local_coding_window() {
    let gene = format!("ATG{}TAA", "GCC".repeat(300));
    let mut right = molecule("actual complete coding source", &gene);
    right.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, gene.len()),
        qualifiers: vec![
            ("label".into(), Some("ordinary coding source".into())),
            ("codon_start".into(), Some("1".into())),
            ("transl_table".into(), Some("11".into())),
        ],
        needs_review: false,
    });
    let target = molecule("independent full expected target", &format!("TTTTTT{gene}"));
    let plan = DesignPlan {
        name: "coding overlap".into(),
        target_revision: target.revision(),
        inputs: vec![
            input("left", molecule("upstream actual", "TTTTTTATGGCCGCCGCCGCC")),
            input("right", right),
        ],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![15],
            circular: false,
        }],
    };
    let report = j::inspect_plan(
        &plan,
        &target,
        &j::Options {
            flank_bases: 12,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    let cds = &report.junctions[0].cds[0];
    assert!(cds.spans_boundary);
    let translated = cds
        .translation
        .as_ref()
        .expect("ordinary >600nt CDS should be inspected");
    assert_eq!(translated.genetic_code, 11);
    assert_eq!(translated.junction_coding_offset, 15);
    assert!(translated.dna.len() <= 27);
    assert_eq!(translated.dna.len(), translated.protein.len() * 3);
    assert!(translated.dna.len() < gene.len());
    assert!(translated.protein.bytes().all(|b| b == b'A'));
}
#[test]
fn failed_final_operation_never_exposes_preparation_as_final_but_wrong_target_does() {
    let mut target = molecule("independent target", "AAACCCGGG");
    let mut plan = DesignPlan {
        name: "failed unknown chemistry".into(),
        target_revision: target.revision(),
        inputs: vec![
            input("a", molecule("first", "AAACCC")),
            input("b", molecule("second", "GGGTTT")),
        ],
        steps: vec![
            Operation::TreatEnds {
                template: Source::Input("a".into()),
                treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Input("b".into()),
                        reverse: false,
                    },
                ],
                circular: false,
                allow_nicked: true,
                options: Default::default(),
            },
        ],
    };
    let report = j::inspect_plan(&plan, &target, &j::Options::default()).unwrap();
    assert_eq!(report.verification.products.len(), 1);
    assert!(report.product_revision.is_none());
    assert!(report.junctions.is_empty());
    assert!(!report.complete);
    let (p, _) = overlap("AAAT", "TGATTT", 1, "AAATGATTT", false);
    plan = p;
    target = molecule("wrong target", "CCCCCCCCC");
    plan.target_revision = target.revision();
    let report = j::inspect_plan(&plan, &target, &j::Options::default()).unwrap();
    assert_eq!(report.verification.status, Status::Invalid);
    assert_eq!(report.junctions[0].window.dna, "AAATGATTT");
    assert!(report.product_revision.is_some());
}
