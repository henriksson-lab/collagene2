//! Annotation expectations are hand specified from molecular strand diagrams.
use collagene_core::{
    analysis::{self, Primer},
    chemistry::EndTreatment,
    verification::*,
    *,
};

fn dna(sequence: &str, topology: Topology) -> Molecule {
    Molecule::dna("fixture", sequence, topology).unwrap()
}
fn feature(label: &str, location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some(label.into())),
            ("translation".into(), Some("SOURCE_ASSERTION".into())),
            ("codon_start".into(), Some("1".into())),
        ],
        needs_review: false,
    }
}
fn named<'a>(m: &'a Molecule, label: &str) -> &'a Feature {
    m.features
        .iter()
        .find(|f| f.qualifiers.contains(&("label".into(), Some(label.into()))))
        .unwrap_or_else(|| panic!("missing feature {label}: {:?}", m.features))
}
fn report(
    inputs: Vec<(&str, Molecule)>,
    steps: Vec<Operation>,
    sequence: &str,
    topology: Topology,
) -> VerificationReport {
    let target = dna(sequence, topology);
    verify(
        &DesignPlan {
            name: "annotation fixture".into(),
            target_revision: target.revision(),
            inputs: inputs
                .into_iter()
                .map(|(id, molecule)| PlanInput {
                    id: id.into(),
                    revision: molecule.revision(),
                    molecule,
                })
                .collect(),
            steps,
        },
        &target,
    )
}
fn consistent(r: &VerificationReport) -> &Molecule {
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    r.products.last().unwrap()
}
fn pcr(template: Source, f: (&str, &str), r: (&str, &str)) -> Operation {
    Operation::Pcr {
        template,
        forward: Primer::new("F", f.0, f.1).unwrap(),
        reverse: Primer::new("R", r.0, r.1).unwrap(),
    }
}
fn input(id: &str) -> Source {
    Source::Input(id.into())
}
fn assert_review(f: &Feature) {
    assert!(f.needs_review);
    assert!(
        !f.qualifiers
            .iter()
            .any(|(k, _)| k == "translation" || k == "codon_start")
    );
}

#[test]
fn pcr_tails_shift_both_annotation_strands_without_annotating_synthetic_bases() {
    let mut m = dna("ACGACCCCGGGGTAAC", Topology::Linear);
    m.features.push(feature("forward", Location::range(0, 4)));
    m.features.push(feature(
        "reverse",
        Location::Complement(Box::new(Location::range(12, 16))),
    ));
    let r = report(
        vec![("t", m)],
        vec![pcr(input("t"), ("ACGA", "CGT"), ("GTTA", "GAC"))],
        "CGTACGACCCCGGGGTAACGTC",
        Topology::Linear,
    );
    let p = consistent(&r);
    assert_eq!(named(p, "forward").location, Location::range(3, 7));
    assert_eq!(
        named(p, "reverse").location,
        Location::Complement(Box::new(Location::range(15, 19)))
    );
    assert_eq!(
        named(p, "reverse").location.extract(&p.sequence).unwrap(),
        "GTTA"
    );
    assert!(p.features.iter().all(|f| !f.needs_review));
}

#[test]
fn circular_pcr_keeps_origin_spanning_gene_and_marks_clipped_gene() {
    let mut m = dna("CCCCTAACGGGGACGA", Topology::Circular);
    m.features.push(feature(
        "origin",
        Location::Join(vec![Location::range(12, 16), Location::range(0, 4)]),
    ));
    m.features.push(feature("clipped", Location::range(6, 10)));
    let r = report(
        vec![("t", m)],
        vec![pcr(input("t"), ("ACGA", "TT"), ("GTTA", "CG"))],
        "TTACGACCCCTAACCG",
        Topology::Linear,
    );
    let p = consistent(&r);
    assert_eq!(
        named(p, "origin").location.extract(&p.sequence).unwrap(),
        "ACGACCCC"
    );
    assert!(!named(p, "origin").needs_review);
    assert_review(named(p, "clipped"));
}

#[test]
fn equal_pcr_sequences_with_distinct_template_origins_remain_ambiguous() {
    // The terminal G can come from template position 12 or the reverse primer tail.
    // Both oligos have ordered sequence CAAAA but distinct declared annealing domains.
    let mut m = dna("AAAACCCCTTTTG", Topology::Linear);
    m.features
        .push(feature("terminal", Location::range(12, 13)));
    let a = Primer::new("A", "AAAA", "C").unwrap();
    let b = Primer::new("B", "CAAAA", "").unwrap();
    let products = analysis::simulate_pcr_geometry(&m, &a, &b).unwrap();
    assert_eq!(products.len(), 2);
    assert!(
        products
            .iter()
            .all(|p| p.product.sequence.as_str() == "CAAAACCCCTTTTG")
    );
    let r = report(
        vec![("t", m)],
        vec![Operation::Pcr {
            template: input("t"),
            forward: a,
            reverse: b,
        }],
        "CAAAACCCCTTTTG",
        Topology::Linear,
    );
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    assert!(r.products.is_empty());
}

#[test]
fn circular_digest_marks_cut_gene_but_preserves_intact_gene() {
    let mut m = dna("AAAAGAATTCTTTT", Topology::Circular);
    m.features.push(feature("cut", Location::range(2, 8)));
    m.features.push(feature("intact", Location::range(8, 12)));
    let r = report(
        vec![("t", m)],
        vec![Operation::Digest {
            template: input("t"),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        }],
        "AATTCTTTTAAAAG",
        Topology::Linear,
    );
    let p = consistent(&r);
    assert_review(named(p, "cut"));
    assert_eq!(named(p, "intact").location, Location::range(3, 7));
    assert_eq!(
        named(p, "intact").location.extract(&p.sequence).unwrap(),
        "TCTT"
    );
    assert!(!named(p, "intact").needs_review);
}

#[test]
fn reversed_staggered_self_ligation_maps_top_overhang_after_closure() {
    // BsaI top cut 10, bottom cut 14. Original top AGTC is absent from the
    // reversed fragment's top reference but paired again after circular closure.
    let mut m = dna("AAAGGTCTCTAGTCCCCC", Topology::Circular);
    m.features
        .push(feature("overhang", Location::range(10, 14)));
    m.features.push(feature("paired", Location::range(14, 18)));
    let r = report(
        vec![("t", m)],
        vec![
            Operation::Digest {
                template: input("t"),
                enzymes: vec!["BsaI".into()],
                fragment_index: 0,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(0),
                    reverse: true,
                }],
                circular: true,
                allow_nicked: false,
            },
        ],
        "GACTAGAGACCTTTGGGG",
        Topology::Circular,
    );
    let p = consistent(&r);
    assert_eq!(
        named(p, "overhang").location,
        Location::Complement(Box::new(Location::range(0, 4)))
    );
    assert_eq!(
        named(p, "overhang").location.extract(&p.sequence).unwrap(),
        "AGTC"
    );
    assert_eq!(
        named(p, "paired").location,
        Location::Complement(Box::new(Location::range(14, 18)))
    );
    assert!(p.features.iter().all(|f| !f.needs_review));
}

#[test]
fn linear_overlap_deduplicates_identical_shared_arm_annotation() {
    let mut a = dna("AAAACCCC", Topology::Linear);
    let mut b = dna("CCCCGGGG", Topology::Linear);
    a.features.push(feature("shared", Location::range(4, 8)));
    b.features.push(feature("shared", Location::range(0, 4)));
    b.features.push(feature("right", Location::range(4, 8)));
    let r = report(
        vec![("a", a), ("b", b)],
        vec![Operation::Overlap {
            fragments: vec![input("a"), input("b")],
            overlaps: vec![4],
            circular: false,
        }],
        "AAAACCCCGGGG",
        Topology::Linear,
    );
    let p = consistent(&r);
    assert_eq!(p.features.len(), 2);
    assert_eq!(named(p, "shared").location, Location::range(4, 8));
    assert_eq!(named(p, "right").location, Location::range(8, 12));
}

#[test]
fn circular_overlap_wraps_last_fragment_without_breaking_intact_gene() {
    let a = dna("AAAACCCCGGGG", Topology::Linear);
    let mut b = dna("GGGGTTTTAAAA", Topology::Linear);
    b.features.push(feature("wrap", Location::range(6, 12)));
    let r = report(
        vec![("a", a), ("b", b)],
        vec![Operation::Overlap {
            fragments: vec![input("a"), input("b")],
            overlaps: vec![4, 4],
            circular: true,
        }],
        "AAAACCCCGGGGTTTT",
        Topology::Circular,
    );
    let p = consistent(&r);
    assert_eq!(
        named(p, "wrap").location,
        Location::Join(vec![Location::range(14, 16), Location::range(0, 4)])
    );
    assert_eq!(
        named(p, "wrap").location.extract(&p.sequence).unwrap(),
        "TTAAAA"
    );
    assert!(!named(p, "wrap").needs_review);
}

#[test]
fn self_circular_overlap_preserves_unique_features_and_reviews_collapsed_gene() {
    let mut m = dna("AAAACCCCGGGGAAAA", Topology::Linear);
    m.features.push(feature("prefix", Location::range(0, 6)));
    m.features.push(feature("suffix", Location::range(10, 16)));
    m.features
        .push(feature("collapsed", Location::range(0, 16)));
    let r = report(
        vec![("t", m)],
        vec![Operation::Overlap {
            fragments: vec![input("t")],
            overlaps: vec![4],
            circular: true,
        }],
        "AAAACCCCGGGG",
        Topology::Circular,
    );
    let p = consistent(&r);
    assert_eq!(
        named(p, "prefix").location.extract(&p.sequence).unwrap(),
        "AAAACC"
    );
    assert_eq!(
        named(p, "suffix").location.extract(&p.sequence).unwrap(),
        "GGAAAA"
    );
    assert!(!named(p, "prefix").needs_review);
    assert!(!named(p, "suffix").needs_review);
    // A 16-base source gene cannot become a fictitious exact 12-base gene.
    if let Some(f) = p.features.iter().find(|f| {
        f.qualifiers
            .contains(&("label".into(), Some("collapsed".into())))
    }) {
        assert_review(f);
    } else {
        assert!(
            p.metadata
                .iter()
                .any(|(k, v)| k.starts_with("annotation_mapping:") && v.contains("collapsed"))
        );
    }
}

#[test]
fn unresolved_features_stay_evidence_and_large_source_blobs_are_not_copied_per_step() {
    let mut m = dna("AAAAGAATTCTTTT", Topology::Circular);
    let unknown = feature("remote", Location::Unresolved("REMOTE.1:5..8".into()));
    m.features.push(unknown.clone());
    m.metadata
        .insert("immutable_source_blob".into(), "X".repeat(100_000));
    let source_revision = m.revision();
    let r = report(
        vec![("t", m)],
        vec![
            Operation::Digest {
                template: input("t"),
                enzymes: vec!["EcoRI".into()],
                fragment_index: 0,
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Dephosphorylate,
            },
            Operation::TreatEnds {
                template: Source::Step(1),
                treatment: EndTreatment::Phosphorylate,
            },
        ],
        "AATTCTTTTAAAAG",
        Topology::Linear,
    );
    consistent(&r);
    assert!(
        r.products
            .iter()
            .all(|p| p.features.is_empty() && !p.metadata.contains_key("immutable_source_blob"))
    );
    let evidence =
        &r.products[0].metadata[&format!("annotation_mapping:unmapped:{source_revision}")];
    let mut expected = unknown;
    expected.needs_review = true;
    assert_eq!(
        serde_json::from_str::<Vec<Feature>>(evidence).unwrap(),
        vec![expected]
    );
    assert!(serde_json::to_vec(&r.products).unwrap().len() < 20_000);
}
