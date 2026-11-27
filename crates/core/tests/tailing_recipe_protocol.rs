//! Source-bound TA/GC recipe previews use one appended worker selector, not a new plan format.
use collagene_core::{
    Molecule, Sequence, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation},
    chemistry::EndTreatment,
    ligation_recipe,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    tailing::TailingBase,
    tailing_recipe::{self, Fragment, Recipe},
    verification::{self, Operation, PlanInput, Status},
    workspace,
};
fn fixture(gc: bool) -> (Vec<PlanInput>, Molecule, Recipe) {
    let inputs = [("vector", "AAACCC"), ("insert", "GGGTTT")]
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
        (
            "vector",
            "AAA",
            "GGG",
            if gc { TailingBase::C } else { TailingBase::T },
        ),
        (
            "insert",
            "GGG",
            "AAA",
            if gc { TailingBase::G } else { TailingBase::A },
        ),
    ]
    .into_iter()
    .map(|(id, f, r, base)| Fragment {
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
        action: collagene_core::tailing_recipe::Action::AddTail { base },
        treatment: Some(EndTreatment::Phosphorylate),
    })
    .collect();
    let target = Molecule::dna(
        "Independent literal tailed14bp circle",
        if gc {
            "AAACCCCGGGTTTG"
        } else {
            "AAACCCTGGGTTTA"
        },
        Topology::Circular,
    )
    .unwrap();
    (
        inputs,
        target,
        Recipe {
            profile: if gc {
                collagene_core::tailing_recipe::Profile::Gc
            } else {
                collagene_core::tailing_recipe::Profile::Ta
            },

            name: "Explicit original-source two-row tailing".into(),
            fragments,
            circular: true,
            allow_nicked: false,
        },
    )
}
fn request(inputs: Vec<PlanInput>, target: Molecule, recipe: Recipe) -> Request {
    Request {
        id: 156,
        revision: tailing_recipe::preview_revision(&inputs, &target, &recipe).unwrap(),
        job: AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        },
    }
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    bytes
}
fn preview(request: &Request) -> (Vec<u8>, tailing_recipe::Preview) {
    let bytes = protocol::process_request(&protocol::encode_request(request).unwrap()).unwrap();
    let response = protocol::decode_response(&bytes).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    let AnalysisResult::TailingRecipePreview(value) = response.result else {
        panic!("expected actual tailing recipe preview148")
    };
    (bytes, value)
}
fn circle_eq(a: &str, b: &str) -> bool {
    let reverse = Sequence::dna(a).unwrap().reverse_complement().unwrap();
    a.len() == b.len()
        && (format!("{a}{a}").contains(b)
            || format!("{}{}", reverse.as_str(), reverse.as_str()).contains(b))
}
#[test]
fn tailing_recipe_reconstructs_literal_ta_and_gc_with_current_operations() {
    assert!(protocol::WORKER_SCHEMA_VERSION >= 70);
    assert!(protocol::CANDIDATE_SCHEMA_VERSION >= 16);
    assert!(workspace::PLAN_SCHEMA_VERSION >= 14);
    for gc in [false, true] {
        let (inputs, target, recipe) = fixture(gc);
        let expected = target.sequence.as_str().to_owned();
        let request = request(inputs, target, recipe);
        let encoded = protocol::encode_request(&request).unwrap();
        let ((_, _, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
        assert_eq!(tag, 154);
        let (bytes, value) = preview(&request);
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 146);
        value.validate_binding(&request.revision).unwrap();
        assert_eq!(
            value.verification.engine_version,
            verification::ENGINE_VERSION
        );
        assert!(value.preparation_complete);
        assert_eq!(
            value.verification.status,
            Status::ConsistentUnderAssumptions
        );
        let physical = value.product_physical.as_ref().unwrap();
        assert_eq!(physical.top_length, 14);
        assert_eq!(physical.bottom_length, 14);
        assert!(physical.nicks.is_empty());
        assert_eq!(physical.topology, Topology::Circular);
        assert!(circle_eq(
            value
                .verification
                .products
                .last()
                .unwrap()
                .sequence
                .as_str(),
            &expected
        ));
        assert_eq!(
            value
                .compiled
                .plan
                .steps
                .iter()
                .filter(|op| matches!(op, Operation::Tail3Prime { .. }))
                .count(),
            2
        );
        let mut old = encoded;
        old[4..6].copy_from_slice(&69u16.to_le_bytes());
        assert!(protocol::process_request(&old).is_err());
        let mut old = bytes;
        old[4..6].copy_from_slice(&69u16.to_le_bytes());
        assert!(protocol::decode_response(&old).is_err());
    }
}
#[test]
fn tail_base_treatment_and_independent_target_are_bound_and_unsealed_nicks_are_retained() {
    let (inputs, target, mut recipe) = fixture(false);
    let prior = tailing_recipe::preview_revision(&inputs, &target, &recipe).unwrap();
    recipe.fragments[1].treatment = None;
    recipe.allow_nicked = true;
    let mut stale = request(inputs.clone(), target.clone(), recipe.clone());
    stale.revision = prior;
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let req = request(inputs.clone(), target, recipe.clone());
    let (_, value) = preview(&req);
    assert_eq!(value.product_physical.as_ref().unwrap().nicks.len(), 2);
    assert!(value.preparation_complete);
    let mut tampered = value.clone();
    tampered.recipe.fragments[0].action = tailing_recipe::Action::AddTail {
        base: TailingBase::A,
    };
    let forged = Response {
        id: req.id,
        revision: req.revision.clone(),
        result: AnalysisResult::TailingRecipePreview(tampered),
    };
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &forged))
            .is_err()
    );
    let wrong = Molecule::dna(
        "Independent wrong14bp target",
        "AAAAAAAAAAAAAA",
        Topology::Circular,
    )
    .unwrap();
    let (_, value) = preview(&request(inputs, wrong, recipe));
    assert_eq!(value.verification.status, Status::Invalid);
}
