//! Three independently specified PCR fragments, explicit phosphorylation and target.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation},
    chemistry::EndTreatment,
    ligation_recipe::{self, Fragment, Recipe},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::{PlanInput, Status},
};
fn fixture() -> (Vec<PlanInput>, Molecule, Recipe) {
    let inputs = [("a", "AAACCC"), ("b", "GGGTTT"), ("c", "ATGCGT")]
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
    let fragments = [
        ("a", "AAA", "GGG"),
        ("b", "GGG", "AAA"),
        ("c", "ATG", "ACG"),
    ]
    .into_iter()
    .map(|(id, f, r)| Fragment {
        source: FragmentSpec {
            generation: None,
            role: id.into(),
            input_id: id.into(),
            preparation: Preparation::ExactPcr {
                forward: Primer::new(format!("{id}-forward"), f, "").unwrap(),
                reverse: Primer::new(format!("{id}-reverse"), r, "").unwrap(),
            },
            reverse: false,
        },
        treatment: Some(EndTreatment::Phosphorylate),
    })
    .collect();
    let target =
        Molecule::dna("independent target", "AAACCCGGGTTTATGCGT", Topology::Linear).unwrap();
    (
        inputs,
        target,
        Recipe {
            name: "three-piece ligation".into(),
            fragments,
            circular: false,
            allow_nicked: false,
        },
    )
}
fn request(inputs: Vec<PlanInput>, target: Molecule, recipe: Recipe) -> Vec<u8> {
    let revision = ligation_recipe::preview_revision(
        &inputs,
        &target,
        &recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    protocol::encode_request(&Request {
        id: 955,
        revision,
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            inputs,
            target,
            recipe,
        },
    })
    .unwrap()
}
#[test]
fn fresh_three_piece_ligation_has_sealed_physical_state_and_schema_thirty_two() {
    let (inputs, target, recipe) = fixture();
    let bytes = request(inputs, target, recipe);
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 93);
    let response = protocol::process_request(&bytes).unwrap();
    let decoded = protocol::decode_response(&response).unwrap();
    let AnalysisResult::LigationRecipePreview(preview) = decoded.result else {
        panic!("expected ligation preview")
    };
    assert!(preview.preparation_complete);
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(preview.prepared.len(), 3);
    let actual = preview.product_physical.unwrap();
    assert_eq!(actual.top_length, 18);
    assert!(actual.nicks.is_empty());
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "AAACCCGGGTTTATGCGT"
    );
    for old in [1u16, 31] {
        let mut downgraded = bytes.clone();
        downgraded[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(protocol::process_request(&downgraded).is_err());
        let mut downgraded = response.clone();
        downgraded[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(protocol::decode_response(&downgraded).is_err());
    }
}
#[test]
fn treatment_and_target_are_part_of_worker_binding_and_retained_nicks_are_explicit() {
    let (inputs, target, mut recipe) = fixture();
    let prior = ligation_recipe::preview_revision(
        &inputs,
        &target,
        &recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    recipe.fragments[0].treatment = None;
    recipe.allow_nicked = true;
    let stale = protocol::encode_request(&Request {
        id: 955,
        revision: prior,
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            inputs: inputs.clone(),
            target: target.clone(),
            recipe: recipe.clone(),
        },
    })
    .unwrap();
    assert!(protocol::process_request(&stale).is_err());
    let fresh = protocol::decode_response(
        &protocol::process_request(&request(inputs.clone(), target, recipe.clone())).unwrap(),
    )
    .unwrap();
    let AnalysisResult::LigationRecipePreview(preview) = fresh.result else {
        panic!("expected ligation preview")
    };
    assert_eq!(preview.product_physical.unwrap().nicks.len(), 1);
    let wrong = Molecule::dna(
        "independent wrong target",
        "AAACCCGGGTTTATGCGA",
        Topology::Linear,
    )
    .unwrap();
    let fresh = protocol::decode_response(
        &protocol::process_request(&request(inputs, wrong, recipe)).unwrap(),
    )
    .unwrap();
    let AnalysisResult::LigationRecipePreview(preview) = fresh.result else {
        panic!("expected ligation preview")
    };
    assert_eq!(preview.verification.status, Status::Invalid);
}
