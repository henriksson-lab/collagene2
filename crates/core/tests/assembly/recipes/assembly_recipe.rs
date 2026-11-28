use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::{assembly::recipe as a, overlap::profiles::ProfileOptions},
    planning::plan::{Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Feature, Location, Molecule, Topology};
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
fn row(role: &str, id: &str, preparation: a::Preparation, reverse: bool) -> a::FragmentSpec {
    a::FragmentSpec {
        generation: None,
        role: role.into(),
        input_id: id.into(),
        preparation,
        reverse,
    }
}
fn recipe(rows: Vec<a::FragmentSpec>) -> a::Recipe {
    a::Recipe {
        name: "explicit original-source recipe".into(),
        overlaps: vec![15; rows.len() - 1],
        fragments: rows,
        circular: false,
        options: ProfileOptions::default(),
    }
}
fn target() -> Molecule {
    molecule("independently selected target", "TTGCAATCGTCGATGCTAGCGGTAC")
}
#[test]
fn original_whole_sources_with_explicit_reverse_are_verified_and_annotations_retained() {
    let mut left = molecule("original left", "TTGCAATCGTCGATGCTAGC");
    left.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("left mark".into()))],
        needs_review: false,
    });
    let right = molecule("actual reverse original", "GTACCGCTAGCATCGACGAT");
    let inputs = vec![input("doc-left", left), input("doc-right", right)];
    let recipe = recipe(vec![
        row("left", "doc-left", a::Preparation::WholeLinear, false),
        row("right", "doc-right", a::Preparation::WholeLinear, true),
    ]);
    let preview = a::preview(&inputs, &target(), &recipe).unwrap();
    assert!(preview.preparation_complete);
    assert_eq!(
        preview.prepared[1].sequence.as_str(),
        "GTACCGCTAGCATCGACGAT"
    );
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    let final_product = preview.verification.products.last().unwrap();
    assert_eq!(final_product.sequence, target().sequence);
    assert!(
        final_product
            .features
            .iter()
            .any(|f| f.label() == "left mark")
    );
    assert!(preview.profile.as_ref().unwrap().geometry_valid);
    assert!(
        preview
            .compiled
            .rows
            .iter()
            .all(|r| r.preparation_step.is_none())
    );
    assert_eq!(
        serde_json::to_value(&preview.compiled.plan.inputs[0].molecule).unwrap(),
        serde_json::to_value(&inputs[0].molecule).unwrap()
    );
    preview
        .validate_binding(&a::preview_revision(&inputs, &target(), &recipe).unwrap())
        .unwrap();
    let mut forged = preview.clone();
    forged.profile.as_mut().unwrap().source_revisions.swap(0, 1);
    assert!(forged.validate().is_err());
}
#[test]
fn actual_pcr_and_literal_ecorv_fragment_compile_without_source_extraction() {
    let inputs = vec![
        input("template", molecule("PCR template", "TTGCAATCGTCGATGCTAGC")),
        input(
            "cut-source",
            molecule("complete digest source", "GATATCGTCGATGCTAGCGGTAC"),
        ),
    ];
    let recipe = recipe(vec![
        row(
            "amplicon",
            "template",
            a::Preparation::ExactPcr {
                forward: Primer::new("ordered forward", "TTGCAATC", "").unwrap(),
                reverse: Primer::new("ordered reverse", "GCTAGCAT", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            false,
        ),
        row(
            "cut-arm",
            "cut-source",
            a::Preparation::Digest {
                enzymes: vec!["EcoRV".into()],
                fragment_index: 1,
                repair: None,
            },
            false,
        ),
    ]);
    let preview = a::preview(&inputs, &target(), &recipe).unwrap();
    assert_eq!(
        preview.prepared[0].sequence.as_str(),
        "TTGCAATCGTCGATGCTAGC"
    );
    assert_eq!(
        preview.prepared[1].sequence.as_str(),
        "ATCGTCGATGCTAGCGGTAC"
    );
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert!(
        matches!(&preview.compiled.plan.steps[0],Operation::Pcr{template:Source::Input(id),..} if id=="template")
    );
    assert!(
        matches!(&preview.compiled.plan.steps[1],Operation::Digest{template:Source::Input(id),fragment_index:1,..} if id=="cut-source")
    );
    assert_eq!(preview.compiled.assembly_step, 2);
    assert_eq!(
        inputs[1].molecule.sequence.as_str(),
        "GATATCGTCGATGCTAGCGGTAC"
    );
}
#[test]
fn repeated_source_roles_keep_order_and_explicit_repair_is_a_separate_step() {
    let inputs = vec![input(
        "one",
        molecule("unchanged source", "GAATTCTTGCAATCGTCGATGCTAGC"),
    )];
    let recipe = recipe(vec![
        row("first-copy", "one", a::Preparation::WholeLinear, false),
        row(
            "prepared-copy",
            "one",
            a::Preparation::Digest {
                enzymes: vec!["EcoRI".into()],
                fragment_index: 1,
                repair: Some(EndRepairMethod::Blunt),
            },
            true,
        ),
    ]);
    let compiled = a::compile(&inputs, &target(), &recipe).unwrap();
    assert_eq!(compiled.plan.inputs.len(), 1);
    assert_eq!(compiled.rows.len(), 2);
    assert_eq!(compiled.rows[1].preparation_step, Some(0));
    assert_eq!(compiled.rows[1].repair_step, Some(1));
    assert!(matches!(compiled.rows[1].prepared_source, Source::Step(1)));
    assert!(matches!(
        compiled.plan.steps[1],
        Operation::RepairEnds {
            template: Source::Step(0),
            method: EndRepairMethod::Blunt
        }
    ));
    assert_eq!(compiled.assembly_step, 2);
}
#[test]
fn circles_stale_bindings_and_failed_preparation_never_become_manufactured_linear_rows() {
    let mut source = molecule("raw circle", "TTGCAATCGTCGATGCTAGC");
    source.topology = Topology::Circular;
    let inputs = vec![input("circle", source)];
    let r = recipe(vec![row(
        "vector",
        "circle",
        a::Preparation::WholeLinear,
        false,
    )]);
    assert!(a::compile(&inputs, &target(), &r).is_err());
    let mut inputs = vec![input("left", molecule("linear", "TTGCAATCGTCGATGCTAGC"))];
    let r = recipe(vec![
        row(
            "failed-pcr",
            "left",
            a::Preparation::ExactPcr {
                forward: Primer::new("absent", "AAAAAAAA", "").unwrap(),
                reverse: Primer::new("absent", "CCCCCCCC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            false,
        ),
        row("later-whole", "left", a::Preparation::WholeLinear, false),
    ]);
    let preview = a::preview(&inputs, &target(), &r).unwrap();
    assert!(!preview.preparation_complete);
    assert!(preview.prepared.is_empty());
    assert!(preview.profile.is_none());
    assert_eq!(preview.verification.status, Status::Invalid);
    inputs[0].revision = "0".repeat(64);
    assert!(a::compile(&inputs, &target(), &r).is_err());
}
