use collagene_core::{
    chemistry::{EndRepairMethod, FivePrime, Polarity},
    cloning::{
        assembly::recipe::{self as a, FragmentSpec, Preparation},
        ligation::recipe as l,
        restriction::catalog as c,
    },
    planning::plan::{Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
fn input(source: &str) -> PlanInput {
    let molecule = Molecule::dna("original source", source, Topology::Linear).unwrap();
    PlanInput {
        id: "s".into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn recipe() -> l::Recipe {
    let selection = c::select(&c::builtin_catalog().unwrap(), &["bsai".into()]).unwrap();
    l::Recipe {
        name: "PCR-added sites then digest".into(),
        fragments: vec![l::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "s".into(),
                preparation: Preparation::PcrCatalogDigest {
                    forward: Primer::new("F", "ACGACCCC", "GGTCTCA").unwrap(),
                    reverse: Primer::new("R", "GTTACCCC", "GGTCTCA").unwrap(),
                    selection,
                    fragment_index: 1,
                    repair: None,
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                reverse: false,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    }
}
fn target(dna: &str) -> Molecule {
    Molecule::dna("independent literal target", dna, Topology::Linear).unwrap()
}
#[test]
fn pcr_added_bsai_sites_have_independent_top_bottom_geometry_and_original_source() {
    let inputs = vec![input("ACGACCCCGGGGTAAC")];
    let before = serde_json::to_value(&inputs).unwrap();
    let report = l::preview(
        &inputs,
        &target("ACGACCCCGGGG"),
        &recipe(),
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(report.verification.products.len(), 3);
    assert_eq!(
        report.verification.products[0].sequence.as_str(),
        "GGTCTCAACGACCCCGGGGTAACTGAGACC"
    );
    assert_eq!(report.prepared[0].sequence.as_str(), "ACGACCCCGGGG");
    let summary = report.physical[0].as_ref().unwrap();
    assert_eq!(
        (
            summary.top_length,
            summary.bottom_length,
            summary.bottom_start
        ),
        (12, 12, 4)
    );
    let (left, right) = summary.ends.as_ref().unwrap();
    assert_eq!(
        (left.polarity, right.polarity),
        (Polarity::FivePrime, Polarity::FivePrime)
    );
    assert_eq!(
        (left.overhang.as_str(), right.overhang.as_str()),
        ("ACGA", "GTTA")
    );
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    let row = &report.compiled.rows[0].preparation;
    assert_eq!(row.preparation_step, Some(0));
    assert_eq!(row.repair_step, None);
    assert!(matches!(row.prepared_source, Source::Step(1)));
    assert!(
        matches!(&report.compiled.plan.steps[0],Operation::Pcr{template:Source::Input(id),..} if id=="s")
    );
    assert!(matches!(
        &report.compiled.plan.steps[1],
        Operation::CatalogDigest {
            template: Source::Step(0),
            fragment_index: 1,
            ..
        }
    ));
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    report.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn optional_repair_is_a_separate_chain_step_and_overlap_receipt_accepts_chain_mapping() {
    let inputs = vec![input("ACGACCCCGGGGTAAC")];
    let mut recipe = recipe();
    if let Preparation::PcrCatalogDigest { repair, .. } =
        &mut recipe.fragments[0].source.preparation
    {
        *repair = Some(EndRepairMethod::Blunt);
    }
    let prepared = a::compile_preparations(
        &inputs,
        &target("ACGACCCCGGGG"),
        &recipe.name,
        &[recipe.fragments[0].source.clone()],
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(prepared.steps.len(), 3);
    assert_eq!(prepared.rows[0].preparation_step, Some(0));
    assert_eq!(prepared.rows[0].repair_step, Some(2));
    assert!(matches!(
        &prepared.steps[2],
        Operation::RepairEnds {
            template: Source::Step(1),
            ..
        }
    ));
    // Even when homology assembly cannot use a staggered product, decoded evidence
    // must retain the honest completed preparation rather than reject its row map.
    let mut source = recipe.fragments[0].source.clone();
    if let Preparation::PcrCatalogDigest { repair, .. } = &mut source.preparation {
        *repair = None;
    }
    let ar = a::Recipe {
        name: "overlap chain receipt".into(),
        fragments: vec![source],
        overlaps: vec![],
        circular: false,
        options: Default::default(),
    };
    let report = a::preview(&inputs, &target("ACGACCCCGGGG"), &ar).unwrap();
    assert!(report.preparation_complete);
    assert_eq!(report.prepared[0].sequence.as_str(), "ACGACCCCGGGG");
    report.validate().unwrap();
}
#[test]
fn wrong_target_internal_cut_and_fragment_index_do_not_select_a_matching_reference() {
    let inputs = vec![input("ACGACCCCGGGGTAAC")];
    let wrong = l::preview(
        &inputs,
        &target("AAAAAAAAAAAA"),
        &recipe(),
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(wrong.verification.status, Status::Invalid);
    assert_eq!(wrong.prepared[0].sequence.as_str(), "ACGACCCCGGGG");
    let mut r = recipe();
    if let Preparation::PcrCatalogDigest { fragment_index, .. } =
        &mut r.fragments[0].source.preparation
    {
        *fragment_index = 99;
    }
    let absent = l::preview(
        &inputs,
        &target("ACGACCCCGGGG"),
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert!(!absent.preparation_complete);
    assert_eq!(absent.verification.products.len(), 1);
    let inside = l::preview(
        &[input("ACGACCCCGGTCTCGGGGTAAC")],
        &target("ACGACCCCGGTCTCGGGG"),
        &recipe(),
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_ne!(
        inside.verification.status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn tampered_catalog_stale_original_and_aggregate_amplicon_work_fail_preflight() {
    let inputs = vec![input("ACGACCCCGGGGTAAC")];
    let t = target("ACGACCCCGGGG");
    let mut r = recipe();
    if let Preparation::PcrCatalogDigest { selection, .. } = &mut r.fragments[0].source.preparation
    {
        selection.catalog_content_hash = "0".repeat(64);
    }
    assert!(l::source_preflight(&[("s", &inputs[0].molecule)], &t, &r).is_err());
    let mut stale = inputs.clone();
    stale[0].revision = "0".repeat(64);
    assert!(
        l::compile(
            &stale,
            &t,
            &recipe(),
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .is_err()
    );
    let huge = input(&"A".repeat(100_000));
    let mut r = recipe();
    if let Preparation::PcrCatalogDigest {
        forward, reverse, ..
    } = &mut r.fragments[0].source.preparation
    {
        *forward = Primer::new("F", &"A".repeat(1000), "").unwrap();
        *reverse = Primer::new("R", &"T".repeat(1000), "").unwrap();
    }
    assert!(matches!(
        l::source_preflight(&[("s", &huge.molecule)], &t, &r),
        Err(collagene_model::Error::TooLarge)
    ));
}
