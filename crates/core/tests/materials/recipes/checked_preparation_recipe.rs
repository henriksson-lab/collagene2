//! Literal source-generation declarations change real preparative cleavage, never target intent.
use collagene_core::{
    analysis::Primer,
    assembly_recipe::{CleavagePolicy, FragmentSpec, Preparation},
    catalog_digest, ligation_recipe as l,
    methylation::MethylationState as M,
    modified_duplex::GenerationAssumptions,
    verification::{self, Operation, PlanInput},
};
use collagene_model::{Molecule, Topology};
fn generation(state: M) -> GenerationAssumptions {
    GenerationAssumptions {
        dam: state,
        dcm: state,
        cpg: state,
    }
}
fn selection(name: &str) -> catalog_digest::CatalogSelection {
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == name)
        .unwrap()
        .id
        .clone();
    catalog_digest::select(&c, &[id]).unwrap()
}
fn input(m: Molecule) -> Vec<PlanInput> {
    vec![PlanInput {
        id: "source".into(),
        revision: m.revision(),
        molecule: m,
    }]
}
fn recipe(preparation: Preparation) -> l::Recipe {
    l::Recipe {
        name: "Explicit original CpG state".into(),
        fragments: vec![l::Fragment {
            source: FragmentSpec {
                role: "fragment".into(),
                input_id: "source".into(),
                preparation,
                reverse: false,
                generation: Some(generation(M::FullyModified)),
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    }
}
#[test]
fn blocked_original_site_cannot_be_cut_by_a_checked_preparation() {
    let source = Molecule::dna(
        "Original explicit fully modified CpG model",
        "AAACCCGGGTTT",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent requested right cut fragment",
        "GGGTTT",
        Topology::Linear,
    )
    .unwrap();
    let inputs = input(source);
    let r = recipe(Preparation::CatalogDigest {
        selection: selection("SmaI"),
        fragment_index: 1,
        repair: None,
    });
    let checked = l::compile(&inputs, &target, &r, CleavagePolicy::PositionalChecked).unwrap();
    assert!(matches!(
        &checked.plan.steps[0],
        Operation::DeclareSourceModifications { .. }
    ));
    assert!(
        matches!(&checked.plan.steps[1],Operation::CheckedDigest{selection,..} if selection.enzymes[0].enzyme.name=="SmaI")
    );
    let result = verification::verify(&checked.plan, &target);
    assert_eq!(result.products.len(), 1, "{:?}", result.findings);
    assert!(
        result.findings.iter().any(|s| s.contains("Blocked")),
        "{:?}",
        result.findings
    );
    let assumed = l::compile(&inputs, &target, &r, CleavagePolicy::AssumedComplete).unwrap();
    assert!(matches!(
        &assumed.plan.steps[1],
        Operation::CatalogDigest { .. }
    ));
    assert_ne!(
        verification::plan_revision(&checked.plan),
        verification::plan_revision(&assumed.plan)
    );
}
#[test]
fn original_modified_template_then_new_pcr_is_unmodified_before_checked_cutting() {
    let source = Molecule::dna(
        "Original fully modified CpG hypothesis",
        "AGTCCGATCCCGGGGCTACGTA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna("Independent cut DNA", "GGGGCTACGTA", Topology::Linear).unwrap();
    let inputs = input(source);
    let r = recipe(Preparation::PcrCatalogDigest {
        forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("R", "TACGTAGC", "").unwrap(),
        selection: selection("SmaI"),
        fragment_index: 1,
        repair: None,
        template_strand: collagene_core::chemistry::Strand::Top,
    });
    let compiled = l::compile(&inputs, &target, &r, CleavagePolicy::PositionalChecked).unwrap();
    assert!(matches!(
        compiled.plan.steps[0],
        Operation::DeclareSourceModifications { .. }
    ));
    assert!(matches!(compiled.plan.steps[1], Operation::Pcr { .. }));
    assert!(matches!(
        compiled.plan.steps[2],
        Operation::CheckedDigest { .. }
    ));
    let result = verification::verify(&compiled.plan, &target);
    assert_eq!(
        result.products.get(2).map(|m| m.sequence.as_str()),
        Some("GGGGCTACGTA"),
        "{:?}",
        result.findings
    );
}
#[test]
fn current_recipe_requires_explicit_generation_and_rejects_state_resetting_preparation() {
    let source = Molecule::dna("Original", "AAACCCGGGTTT", Topology::Linear).unwrap();
    let target = source.clone();
    let inputs = input(source);
    let mut r = recipe(Preparation::WholeLinear);
    let mut json = serde_json::to_value(&r).unwrap();
    json["fragments"][0]["source"]
        .as_object_mut()
        .unwrap()
        .remove("generation");
    assert!(serde_json::from_value::<l::Recipe>(json).is_err());
    r.fragments[0].source.preparation = Preparation::DeclaredDuplex {
        left: collagene_core::chemistry::Terminus::CLEAVED,
        right: collagene_core::chemistry::Terminus::CLEAVED,
    };
    assert!(l::source_preflight(&[(&inputs[0].id, &inputs[0].molecule)], &target, &r).is_err());
}
