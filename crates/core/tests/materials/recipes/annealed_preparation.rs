//! Literal adaptor construction: the target is declared independently of annealing/digest output.
use collagene_core::{
    chemistry::{AnnealOptions, EndTreatment, Terminus},
    cloning::{
        assembly::recipe::{self as assembly_recipe, FragmentSpec, Preparation},
        ligation::recipe as lig,
        restriction::catalog as catalog_digest,
    },
    planning::plan::{Operation, PlanInput, Source},
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
const TARGET: &str = "AATTCTTACGGGAATTGCTACG";
fn input(id: &str, letters: &str, topology: Topology) -> PlanInput {
    let molecule = Molecule::dna(id, letters, topology).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture() -> (Vec<PlanInput>, Molecule, lig::Recipe) {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let inputs = vec![
        input("backbone", "GGAATTCTTACG", Topology::Circular),
        input("top", "AATTGCTACG", Topology::Linear),
        input("bottom", "AATTCGTAGC", Topology::Linear),
    ];
    let target = Molecule::dna(
        "Independently specified22bp circle",
        TARGET,
        Topology::Circular,
    )
    .unwrap();
    let recipe = lig::Recipe {
        standard: None,
        name: "Original backbone and two separately ordered oligos".into(),
        circular: true,
        allow_nicked: false,
        fragments: vec![
            lig::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: "vector".into(),
                    input_id: "backbone".into(),
                    preparation: Preparation::CatalogDigest {
                        conditions: None,
                        selection: catalog_digest::select(&catalog, &[id]).unwrap(),
                        fragment_index: 0,
                        repair: None,
                    },
                    reverse: false,
                },
                treatment: None,
            },
            lig::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: "adaptor".into(),
                    input_id: "top".into(),
                    preparation: Preparation::AnnealedOligos {
                        bottom_input_id: "bottom".into(),
                        options: AnnealOptions {
                            min_paired_bases: 6,
                            ..Default::default()
                        },
                        candidate_index: 0,
                    },
                    reverse: false,
                },
                treatment: Some(EndTreatment::Phosphorylate),
            },
        ],
    };
    (inputs, target, recipe)
}
#[test]
fn original_catalog_backbone_and_actual_annealed_pair_ligate_to_literal_intent() {
    let (inputs, target, recipe) = fixture();
    let before = serde_json::to_vec(&inputs).unwrap();
    let preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    preview
        .validate_binding(
            &lig::preview_revision(
                &inputs,
                &target,
                &recipe,
                collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
            )
            .unwrap(),
        )
        .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        TARGET
    );
    assert!(
        preview
            .product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    assert_eq!(preview.junctions.len(), 2);
    assert_eq!(preview.compiled.plan.inputs.len(), 3);
    assert!(
        matches!(&preview.compiled.plan.steps[1], Operation::AnnealOligos { top: Source::Input(top), bottom: Source::Input(bottom), candidate_index: 0, .. } if top == "top" && bottom == "bottom")
    );
    assert!(matches!(
        &preview.compiled.plan.steps[2],
        Operation::TreatEnds {
            template: Source::Step(1),
            treatment: EndTreatment::Phosphorylate
        }
    ));
    assert_eq!(before, serde_json::to_vec(&inputs).unwrap());
    let decoded: lig::Preview =
        serde_json::from_slice(&serde_json::to_vec(&preview).unwrap()).unwrap();
    decoded.validate_binding(&preview.context_revision).unwrap();
}
#[test]
fn stale_secondary_binding_and_row_source_substitution_are_rejected() {
    let (mut inputs, target, recipe) = fixture();
    let mut preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    inputs[2].molecule.sequence = collagene_model::Sequence::dna("AATTCGTAGG").unwrap();
    assert!(
        lig::compile(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .is_err()
    );
    assert!(
        lig::preview_revision(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .is_err()
    );
    preview.compiled.plan.inputs[2].molecule = inputs[2].molecule.clone();
    assert!(preview.validate().is_err());
    let (inputs, target, recipe) = fixture();
    let mut preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    if let Operation::AnnealOligos { top, .. } = &mut preview.compiled.plan.steps[1] {
        *top = Source::Input("bottom".into());
    }
    preview.verification.plan_revision =
        collagene_core::planning::plan::DesignPlan::revision(&preview.compiled.plan);
    assert!(preview.validate().is_err());
}
#[test]
fn candidate_and_chemistry_changes_require_fresh_search_and_never_reuse_stock() {
    let (inputs, target, mut recipe) = fixture();
    let original = lig::preview_revision(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    if let Preparation::AnnealedOligos {
        candidate_index, ..
    } = &mut recipe.fragments[1].source.preparation
    {
        *candidate_index = 1;
    }
    assert_ne!(
        original,
        lig::preview_revision(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .unwrap()
    );
    let unavailable = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert!(!unavailable.preparation_complete);
    assert!(unavailable.product_physical.is_none());
    if let Preparation::AnnealedOligos {
        candidate_index,
        options,
        ..
    } = &mut recipe.fragments[1].source.preparation
    {
        *candidate_index = 0;
        options.top_chemistry = Terminus::UNKNOWN;
        options.bottom_chemistry = Terminus::UNKNOWN;
    }
    let unknown = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_ne!(
        unknown.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(unknown.product_physical.is_none());
    let (inputs, mut wrong, recipe) = fixture();
    wrong.sequence = collagene_model::Sequence::dna("AATTCTTACGGGAATTGCTACC").unwrap();
    assert_eq!(
        lig::preview(
            &inputs,
            &wrong,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .unwrap()
        .verification
        .status,
        Status::Invalid
    );
}
#[test]
fn borrowed_guards_reject_missing_duplex_ambiguous_and_unbounded_oligos() {
    for change in 0..5 {
        let (mut inputs, target, mut recipe) = fixture();
        match change {
            0 => {
                inputs.pop();
            }
            1 => inputs[2].molecule.topology = Topology::Circular,
            2 => {
                inputs[2].molecule.sequence = collagene_model::Sequence::dna("AATTCGTAGN").unwrap()
            }
            3 => {
                inputs[2]
                    .molecule
                    .metadata
                    .insert("genbank.molecule_type".into(), "ds-DNA".into());
            }
            _ => {
                if let Preparation::AnnealedOligos { options, .. } =
                    &mut recipe.fragments[1].source.preparation
                {
                    options.max_candidates = usize::MAX;
                }
            }
        }
        let borrowed = inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>();
        assert!(
            lig::source_preflight(&borrowed, &target, &recipe).is_err(),
            "case {change}"
        );
    }
}
#[test]
fn shared_overlap_compiler_retains_the_real_secondary_input_and_operation() {
    let (inputs, target, recipe) = fixture();
    let rows = recipe
        .fragments
        .iter()
        .map(|f| f.source.clone())
        .collect::<Vec<_>>();
    let prepared = assembly_recipe::compile_preparations(
        &inputs,
        &target,
        "Shared native preparation",
        &rows,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(prepared.inputs.len(), 3);
    assert_eq!(prepared.rows.len(), 2);
    assert!(matches!(
        prepared.steps[1],
        Operation::AnnealOligos {
            candidate_index: 0,
            ..
        }
    ));
    assert!(matches!(prepared.rows[1].prepared_source, Source::Step(1)));
}

#[test]
fn diagnostics_locate_two_original_annealing_without_accepting_forged_sources() {
    use collagene_core::cloning::assembly::diagnostics::{self, RowState, Stage};
    let (inputs, target, mut recipe) = fixture();
    let preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    let located = diagnostics::ligation(&preview).unwrap();
    assert_eq!(located.rows.len(), 2);
    assert!(
        located
            .rows
            .iter()
            .all(|r| r.state == RowState::Reconstructed)
    );
    assert_eq!(located.rows[1].preparation_step, Some(1));
    assert_eq!(located.rows[1].total_steps, 2); // actual anneal and explicit phosphorylation
    assert!(located.failures.is_empty());
    for forged in [
        (
            Source::Input("backbone".into()),
            Source::Input("bottom".into()),
        ),
        (Source::Input("top".into()), Source::Input("top".into())),
        (Source::Input("top".into()), Source::Input("unbound".into())),
        (Source::Input("top".into()), Source::Step(0)),
    ] {
        let mut changed = preview.clone();
        if let Operation::AnnealOligos { top, bottom, .. } = &mut changed.compiled.plan.steps[1] {
            *top = forged.0;
            *bottom = forged.1;
        }
        changed.verification.plan_revision =
            collagene_core::planning::plan::DesignPlan::revision(&changed.compiled.plan);
        assert!(diagnostics::ligation(&changed).is_err());
    }
    if let Preparation::AnnealedOligos {
        candidate_index, ..
    } = &mut recipe.fragments[1].source.preparation
    {
        *candidate_index = 1;
    }
    let missing = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    let located = diagnostics::ligation(&missing).unwrap();
    let failure = located
        .failures
        .into_iter()
        .next()
        .expect("one documented failure");
    assert_eq!(
        (failure.step, failure.row, failure.stage),
        (1, Some(1), Stage::Annealing)
    );
    assert_eq!(located.rows[1].state, RowState::Rejected);
}
