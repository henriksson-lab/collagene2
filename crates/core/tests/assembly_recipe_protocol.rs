//! Fresh preparation/profile preview uses original source roles and an independently supplied target.
use collagene_core::{
    Molecule, Topology,
    assembly_recipe::{self, FragmentSpec, Preparation, Recipe},
    overlap_profiles::ProfileOptions,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::PlanInput,
};
fn fixture() -> (Vec<PlanInput>, Molecule, Recipe) {
    let inputs = [
        ("left", "TTGCAACGTCGATGCTAGCA"),
        ("right", "ACGTCGATGCTAGCAGGTAC"),
    ]
    .into_iter()
    .map(|(id, dna)| {
        let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
        PlanInput {
            id: id.into(),
            revision: molecule.revision(),
            molecule,
        }
    })
    .collect();
    let target = Molecule::dna(
        "independent target",
        "TTGCAACGTCGATGCTAGCAGGTAC",
        Topology::Linear,
    )
    .unwrap();
    let recipe = Recipe {
        name: "Literal ordered overlap".into(),
        fragments: ["left", "right"]
            .into_iter()
            .map(|id| FragmentSpec {
                generation: None,
                role: id.into(),
                input_id: id.into(),
                preparation: Preparation::WholeLinear,
                reverse: false,
            })
            .collect(),
        overlaps: vec![15],
        circular: false,
        options: ProfileOptions::default(),
    };
    (inputs, target, recipe)
}
fn request(inputs: Vec<PlanInput>, target: Molecule, recipe: Recipe) -> Vec<u8> {
    let revision = assembly_recipe::preview_revision(&inputs, &target, &recipe).unwrap();
    protocol::encode_request(&Request {
        id: 314,
        revision,
        job: AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        },
    })
    .unwrap()
}
#[test]
fn schema_thirty_one_keeps_original_roles_target_and_exact_junction_evidence() {
    let (inputs, target, recipe) = fixture();
    let expected = assembly_recipe::preview_revision(&inputs, &target, &recipe).unwrap();
    let bytes = request(inputs.clone(), target.clone(), recipe.clone());
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 92);
    let output = protocol::process_request(&bytes).unwrap();
    let response = protocol::decode_response(&output).unwrap();
    assert_eq!(response.revision, expected);
    let AnalysisResult::OverlapRecipePreview(preview) = response.result else {
        panic!("expected prepared overlap preview")
    };
    assert_eq!(preview.context_revision, expected);
    assert!(preview.preparation_complete);
    assert_eq!(preview.prepared.len(), 2);
    assert_eq!(preview.compiled.plan.target_revision, target.revision());
    assert_eq!(preview.compiled.plan.inputs[0].revision, inputs[0].revision);
    let profile = preview.profile.unwrap();
    assert!(profile.geometry_valid && profile.complete);
    assert_eq!(profile.junctions[0].homology, "ACGTCGATGCTAGCA");
    for version in [1u16, 29, 30] {
        let mut old = bytes.clone();
        old[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(protocol::process_request(&old).is_err());
        let mut old = output.clone();
        old[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(protocol::decode_response(&old).is_err());
    }
}
#[test]
fn target_or_preparation_changes_cannot_reuse_an_old_context() {
    let (inputs, target, recipe) = fixture();
    let old = assembly_recipe::preview_revision(&inputs, &target, &recipe).unwrap();
    let changed_target = target.replace(0, 1, "C").unwrap();
    let mut changed_recipe = recipe.clone();
    changed_recipe.fragments[0].reverse = true;
    for (target, recipe) in [(changed_target, recipe), (target, changed_recipe)] {
        let bytes = protocol::encode_request(&Request {
            id: 314,
            revision: old.clone(),
            job: AnalysisJob::PreviewOverlapRecipe {
                inputs: inputs.clone(),
                target,
                recipe,
            },
        })
        .unwrap();
        assert!(protocol::process_request(&bytes).is_err());
    }
}
