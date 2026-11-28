use collagene_core::{
    assembly_recipe as a, catalog_digest as c,
    chemistry::{self, Duplex, EndRepairMethod, FivePrime, Polarity, ThreePrime},
    ligation_recipe as l,
    modified_duplex::{self, GenerationAssumptions, TrackedDuplex},
    verification::{self as v, DesignPlan, Operation, PlanInput, Source, Status},
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn selection() -> c::CatalogSelection {
    c::select(&c::builtin_catalog().unwrap(), &["nhei".into()]).unwrap()
}
fn circle(name: &str, dna: &str) -> Molecule {
    Molecule::dna(name, dna, Topology::Circular).unwrap()
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn single(m: Molecule) -> (DesignPlan, Molecule) {
    let target = Molecule::dna(
        "independent cut reference",
        "CTAGCCCCAAAAG",
        Topology::Linear,
    )
    .unwrap();
    (
        DesignPlan {
            name: "pinned nonstarter digest".into(),
            target_revision: target.revision(),
            inputs: vec![input("original", m)],
            steps: vec![Operation::CatalogDigest {
                template: Source::Input("original".into()),
                selection: selection(),
                fragment_index: 0,
            }],
        },
        target,
    )
}
#[test]
fn nonstarter_circle_digest_has_literal_both_strands_and_cross_origin_annotation() {
    assert!(
        !collagene_core::analysis::starter_enzymes()
            .iter()
            .any(|e| e.name == "NheI")
    );
    let mut source = circle("original circle", "AAAAGCTAGCCCC");
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Join(vec![Location::range(11, 13), Location::range(0, 2)]),
        qualifiers: vec![("label".into(), Some("source origin".into()))],
        needs_review: false,
    });
    let original = serde_json::to_value(&source).unwrap();
    let sel = selection();
    let enzymes = sel
        .resolved_enzymes(&c::builtin_catalog().unwrap())
        .unwrap();
    let tracked =
        TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(&source).unwrap()).unwrap();
    let fragments = modified_duplex::digest(&tracked, &enzymes).unwrap();
    assert_eq!(fragments.len(), 1);
    let d = &fragments[0].duplex;
    assert_eq!(d.top.as_str(), "CTAGCCCCAAAAG");
    assert_eq!(d.bottom.as_str(), "CTAGCTTTTGGGG");
    assert_eq!(d.bottom_start, 4);
    let (left, right) = d.ends().unwrap();
    assert_eq!(
        (left.polarity, right.polarity),
        (Polarity::FivePrime, Polarity::FivePrime)
    );
    assert_eq!(
        (left.overhang.as_str(), right.overhang.as_str()),
        ("CTAG", "CTAG")
    );
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.three_prime, ThreePrime::Hydroxyl);
    let (plan, target) = single(source.clone());
    let report = v::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products[0].sequence.as_str(), "CTAGCCCCAAAAG");
    assert_eq!(
        report.products[0].features[0]
            .location
            .extract(&report.products[0].sequence)
            .unwrap(),
        "CCAA"
    );
    assert_eq!(serde_json::to_value(source).unwrap(), original);
    assert_eq!(postcard::to_allocvec(&plan.steps[0]).unwrap()[0], 17);
}
#[test]
fn shared_pinned_preparation_feeds_actual_cohesive_ligation_against_independent_target() {
    let inputs = vec![
        input("a", circle("first source", "AAAAGCTAGCCCC")),
        input("b", circle("second source", "TTTTGCTAGCGGG")),
    ];
    let original = serde_json::to_value(&inputs).unwrap();
    let target = circle("independent expected circle", "CTAGCCCCAAAAGCTAGCGGGTTTTG");
    let recipe = l::Recipe {
        name: "nonstarter catalog ligation".into(),
        fragments: ["a", "b"]
            .into_iter()
            .map(|id| l::Fragment {
                source: a::FragmentSpec {
                    generation: None,
                    role: id.into(),
                    input_id: id.into(),
                    preparation: a::Preparation::CatalogDigest {
                        selection: selection(),
                        fragment_index: 0,
                        repair: None,
                    },
                    reverse: false,
                },
                treatment: None,
            })
            .collect(),
        circular: true,
        allow_nicked: false,
    };
    let preview = l::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert!(
        preview
            .product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    assert_eq!(preview.junctions.len(), 2);
    assert_eq!(
        preview.verification.products.last().unwrap().sequence,
        target.sequence
    );
    assert_eq!(serde_json::to_value(inputs).unwrap(), original);
    assert_eq!(
        postcard::to_allocvec(&recipe.fragments[0].source.preparation).unwrap()[0],
        4
    );
}
#[test]
fn optional_repair_is_explicit_and_uses_existing_tracked_end_transform() {
    let inputs = vec![input("a", circle("source", "AAAAGCTAGCCCC"))];
    let target = Molecule::dna(
        "independent filled reference",
        "CTAGCCCCAAAAGCTAG",
        Topology::Linear,
    )
    .unwrap();
    let recipe = a::Recipe {
        name: "explicit postdigest repair".into(),
        fragments: vec![a::FragmentSpec {
            generation: None,
            role: "arm".into(),
            input_id: "a".into(),
            preparation: a::Preparation::CatalogDigest {
                selection: selection(),
                fragment_index: 0,
                repair: Some(EndRepairMethod::Blunt),
            },
            reverse: false,
        }],
        overlaps: vec![],
        circular: false,
        options: Default::default(),
    };
    let preview = a::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(preview.compiled.rows[0].preparation_step, Some(0));
    assert_eq!(preview.compiled.rows[0].repair_step, Some(1));
    assert!(matches!(
        preview.compiled.plan.steps[1],
        Operation::RepairEnds {
            template: Source::Step(0),
            method: EndRepairMethod::Blunt
        }
    ));
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.prepared[0].sequence, target.sequence);
}
#[test]
fn tampered_catalog_pin_stale_original_and_detached_physical_evidence_are_rejected() {
    let (mut plan, target) = single(circle("source", "AAAAGCTAGCCCC"));
    if let Operation::CatalogDigest { selection, .. } = &mut plan.steps[0] {
        selection.enzymes[0].enzyme.top_cut = 2;
    }
    let bad = v::verify(&plan, &target);
    assert_eq!(bad.status, Status::Invalid);
    assert!(bad.products.is_empty());
    let (mut plan, target) = single(circle("source", "AAAAGCTAGCCCC"));
    plan.inputs[0].revision = "0".repeat(64);
    assert_eq!(v::verify(&plan, &target).status, Status::Invalid);
    let mut source = circle("detached physical source", "AAAAGCTAGCCCC");
    source
        .metadata
        .insert("duplex_model".into(), "forged known ends".into());
    let (plan, target) = single(source);
    let r = v::verify(&plan, &target);
    assert_eq!(r.status, Status::Incomplete);
    assert!(r.products.is_empty());
}
#[test]
fn catalog_pin_does_not_upgrade_unknown_tracked_methylation_to_complete_cleavage() {
    // BamHI GGATCC includes a real Dam GATC context. The NheI fixture above
    // has no eligible Dam/Dcm/CpG sites and correctly creates no uncertainty.
    let source = circle("unknown Dam source", "AAGGATCCTT");
    let target = Molecule::dna("independent BamHI cut", "GATCCTTAAG", Topology::Linear).unwrap();
    let selection = c::select(&c::builtin_catalog().unwrap(), &["bamhi".into()]).unwrap();
    let mut plan = DesignPlan {
        name: "catalog digest of unknown Dam material".into(),
        target_revision: target.revision(),
        inputs: vec![input("original", source)],
        steps: vec![Operation::CatalogDigest {
            template: Source::Input("original".into()),
            selection,
            fragment_index: 0,
        }],
    };
    plan.steps.insert(
        0,
        Operation::DeclareSourceModifications {
            template: Source::Input("original".into()),
            assumptions: GenerationAssumptions::default(),
        },
    );
    if let Operation::CatalogDigest { template, .. } = &mut plan.steps[1] {
        *template = Source::Step(0);
    }
    let r = v::verify(&plan, &target);
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    assert_eq!(r.products.len(), 2);
    assert!(
        r.findings
            .iter()
            .any(|s| s.contains("assumed-complete cleavage"))
    );
    let p = &r.products[1];
    let summary: chemistry::DuplexSummary =
        serde_json::from_str(&p.metadata["duplex_model"]).unwrap();
    assert_eq!(summary.ends.unwrap().0.overhang, "GATC");
}

#[test]
fn reversed_digest_intent_uses_actual_global_bottom_interval_not_reversed_top_region() {
    let inputs = vec![input("a", circle("source", "AAAAGCTAGCCCC"))];
    let target = Molecule::dna(
        "independent reversed actual top",
        "CTAGCTTTTGGGG",
        Topology::Linear,
    )
    .unwrap();
    let recipe = l::Recipe {
        name: "explicit reversed cohesive arm".into(),
        fragments: vec![l::Fragment {
            source: a::FragmentSpec {
                generation: None,
                role: "reversed".into(),
                input_id: "a".into(),
                preparation: a::Preparation::CatalogDigest {
                    selection: selection(),
                    fragment_index: 0,
                    repair: None,
                },
                reverse: true,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    };
    let mut plan = l::compile(
        &inputs,
        &target,
        &recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap()
    .plan;
    let report = v::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products[0].sequence.as_str(), "CTAGCCCCAAAAG");
    assert_eq!(report.products[1].sequence.as_str(), "CTAGCTTTTGGGG");
    let inspected = c::inspect(
        &c::builtin_catalog().unwrap(),
        &selection(),
        "a",
        &inputs[0].molecule,
        &c::InspectionOptions::default(),
    )
    .unwrap();
    let summary = &inspected.inspection.fragments[0];
    assert_eq!(summary.top_start, 5);
    assert_eq!(summary.bottom_start, 9);
    assert_eq!(summary.duplex.bottom_start, 4);
    assert_eq!(summary.duplex.bottom_length, 13);
    let declared = Location::Complement(Box::new(Location::Join(vec![
        Location::range(9, 13),
        Location::range(0, 9),
    ])))
    .extract(&inputs[0].molecule.sequence)
    .unwrap();
    assert_eq!(declared, "CTAGCTTTTGGGG");
    let wrong = Molecule::dna(
        "wrong naive reversed reference",
        "CTTTTGGGGCTAG",
        Topology::Linear,
    )
    .unwrap();
    plan.target_revision = wrong.revision();
    assert_eq!(v::verify(&plan, &wrong).status, Status::Invalid);
}
