//! Appended preparation capability is gated separately from old generic annealing plans.
use collagene_core::{
    chemistry::{AnnealOptions, Terminus},
    cloning::{
        assembly::recipe::{self as assembly_recipe, FragmentSpec, Preparation},
        ligation::{recipe as ligation_recipe, tailing::TailingBase, tailing_recipe},
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification,
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Vec<PlanInput>, Molecule, ligation_recipe::Recipe) {
    let inputs = [
        ("vector", "AAAAGAATTCTTTT", Topology::Circular),
        ("top", "AATTCCGATGCA", Topology::Linear),
        ("bottom", "AATTTGCATCGG", Topology::Linear),
    ]
    .into_iter()
    .map(|(id, dna, topology)| {
        let molecule = Molecule::dna(id, dna, topology).unwrap();
        PlanInput {
            id: id.into(),
            revision: molecule.revision(),
            molecule,
        }
    })
    .collect();
    let mut options = AnnealOptions::default();
    options.min_paired_bases = 8;
    options.top_chemistry = Terminus::CLEAVED;
    options.bottom_chemistry = Terminus::CLEAVED;
    let fragments = vec![
        ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "vector".into(),
                input_id: "vector".into(),
                preparation: Preparation::Digest {
                    enzymes: vec!["EcoRI".into()],
                    fragment_index: 0,
                    repair: None,
                },
                reverse: false,
            },
            treatment: None,
        },
        ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "top".into(),
                preparation: Preparation::AnnealedOligos {
                    bottom_input_id: "bottom".into(),
                    options,
                    candidate_index: 0,
                },
                reverse: false,
            },
            treatment: None,
        },
    ];
    (
        inputs,
        Molecule::dna(
            "independent26bp circle",
            "AATTCTTTTAAAAGAATTCCGATGCA",
            Topology::Circular,
        )
        .unwrap(),
        ligation_recipe::Recipe {
            standard: None,
            name: "Actual originals".into(),
            fragments,
            circular: true,
            allow_nicked: false,
        },
    )
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    b
}

// Explicit historical operation4 fixture; current recipe compilers intentionally emit operation29.

#[test]
fn all_three_recipe_requests_reject_downgraded_preparation_and_preview_rows_require71() {
    assert!(protocol::schema::WORKER_SCHEMA_VERSION >= 71);
    assert!(protocol::schema::CANDIDATE_SCHEMA_VERSION >= 16);
    let (inputs, target, recipe) = fixture();
    let request = Request {
        id: 1,
        revision: ligation_recipe::preview_revision(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
            inputs: inputs.clone(),
            target: target.clone(),
            recipe: recipe.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    let bytes = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    let response = protocol::response::Response::decode(&bytes).unwrap();
    let AnalysisResult::LigationRecipePreview(p) = &response.result else {
        panic!("ligation preview")
    };
    assert_eq!(
        p.verification.status,
        verification::evidence::Status::ConsistentUnderAssumptions
    );
    assert_eq!(p.compiled.plan.inputs.len(), 3);
    assert_eq!(p.compiled.rows.len(), 2);
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 70, &response)).is_err());
    // New compiler output requires outcomes capability; explicit historical
    // generic operations retain their original worker70 contract.
    let verify = Request {
        id: 2,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: p.compiled.plan.clone(),
            target: target.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 70, &verify)).is_err());
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 91, &verify)).is_err());
    let overlap = assembly_recipe::Recipe {
        standard: None,
        name: recipe.name.clone(),
        fragments: recipe.fragments.iter().map(|f| f.source.clone()).collect(),
        overlaps: vec![4, 4],
        circular: true,
        options: Default::default(),
    };
    let request = Request {
        id: 3,
        revision: assembly_recipe::preview_revision(&inputs, &target, &overlap).unwrap(),
        job: AnalysisJob::PreviewOverlapRecipe {
            inputs: inputs.clone(),
            target: target.clone(),
            recipe: overlap,
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    let tail = tailing_recipe::Recipe {
        standard: None,
        profile: collagene_core::cloning::ligation::tailing_recipe::Profile::Ta,

        name: recipe.name,
        fragments: recipe
            .fragments
            .into_iter()
            .map(|f| tailing_recipe::Fragment {
                source: f.source,
                action: collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                    base: TailingBase::A,
                },
                treatment: f.treatment,
            })
            .collect(),
        circular: true,
        allow_nicked: false,
    };
    let request = Request {
        id: 4,
        revision: tailing_recipe::preview_revision(&inputs, &target, &tail).unwrap(),
        job: AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe: tail,
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    let bytes = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    let response = protocol::response::Response::decode(&bytes).unwrap();
    assert!(matches!(
        response.result,
        AnalysisResult::TailingRecipePreview(_)
    ));
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 70, &response)).is_err());
}
