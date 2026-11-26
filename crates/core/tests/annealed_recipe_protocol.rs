//! Appended preparation capability is gated separately from old generic annealing plans.
use collagene_core::{
    Molecule, Topology,
    assembly_recipe::{self, FragmentSpec, Preparation},
    chemistry::{AnnealOptions, Terminus},
    ligation_recipe,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    tailing::TailingBase,
    tailing_recipe,
    verification::{self, PlanInput},
};
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
#[test]
fn all_three_recipe_requests_reject_downgraded_preparation_and_preview_rows_require71() {
    assert_eq!(protocol::WORKER_SCHEMA_VERSION, 71);
    assert_eq!(protocol::CANDIDATE_SCHEMA_VERSION, 16);
    let (inputs, target, recipe) = fixture();
    let request = Request {
        id: 1,
        revision: ligation_recipe::preview_revision(&inputs, &target, &recipe).unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            inputs: inputs.clone(),
            target: target.clone(),
            recipe: recipe.clone(),
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    let bytes = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    let response = protocol::decode_response(&bytes).unwrap();
    let AnalysisResult::LigationRecipePreview(p) = &response.result else {
        panic!("ligation preview")
    };
    assert_eq!(
        p.verification.status,
        verification::Status::ConsistentUnderAssumptions
    );
    assert_eq!(p.compiled.plan.inputs.len(), 3);
    assert_eq!(p.compiled.rows.len(), 2);
    assert!(protocol::decode_response(&wire(b"CGR\0", 70, &response)).is_err());
    // The same compiled generic plan still consists of established operations;
    // raising the recipe capability never raises candidate or Verify semantics.
    let verify = Request {
        id: 2,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: p.compiled.plan.clone(),
            target: target.clone(),
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &verify)).is_ok());
    let overlap = assembly_recipe::Recipe {
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
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    let tail = tailing_recipe::Recipe {
        name: recipe.name,
        fragments: recipe
            .fragments
            .into_iter()
            .map(|f| tailing_recipe::Fragment {
                source: f.source,
                base: TailingBase::A,
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
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    let bytes = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    let response = protocol::decode_response(&bytes).unwrap();
    assert!(matches!(
        response.result,
        AnalysisResult::TailingRecipePreview(_)
    ));
    assert!(protocol::decode_response(&wire(b"CGR\0", 70, &response)).is_err());
}
#[test]
fn ordinary_existing_recipe_and_generic_annealing_plan_keep_old70_compatibility() {
    let (inputs, target, mut recipe) = fixture();
    recipe.fragments[1].source.preparation = Preparation::WholeLinear;
    let request = Request {
        id: 10,
        revision: ligation_recipe::preview_revision(&inputs, &target, &recipe).unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &request)).is_ok());
    let bytes = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    let response: Response = protocol::decode_response(&bytes).unwrap();
    assert!(protocol::decode_response(&wire(b"CGR\0", 70, &response)).is_ok());
}
#[test]
fn paired_editable_draft_attachment_cannot_claim_old70_workspace_transport() {
    use collagene_core::workspace::{Workspace, WorkspaceAttachment};
    let mut workspace = Workspace::from_documents(
        "opaque editable draft",
        vec![Molecule::dna("source", "ACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let payload = br#"{"version":4}"#.to_vec();
    workspace.attachments.insert(
        blake3::hash(&payload).to_hex().to_string(),
        WorkspaceAttachment {
            kind: "assembly-draft-v4".into(),
            name: "editable paired draft".into(),
            media_type: "application/json".into(),
            bytes: payload,
        },
    );
    let request = Request {
        id: 20,
        revision: workspace
            .head(&workspace.active_document)
            .unwrap()
            .revision(),
        job: AnalysisJob::SaveWorkspace {
            workspace: workspace.clone(),
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &request)).is_err());
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_ok());
    // Previous editable attachment kinds remain transport-compatible. This is
    // an opaque-history gate, never a claim to validate the web draft schema.
    workspace.attachments.values_mut().next().unwrap().kind = "assembly-draft-v3".into();
    let request = Request {
        id: 21,
        revision: workspace
            .head(&workspace.active_document)
            .unwrap()
            .revision(),
        job: AnalysisJob::SaveWorkspace { workspace },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 70, &request)).is_ok());
}
