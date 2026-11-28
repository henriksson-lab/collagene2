//! A pinned catalog operation cannot hide in historical plan/recipe/archive headers.
use collagene_core::{
    cloning::restriction::catalog::{self as catalog_digest, CatalogSelection},
    documents::{
        transfer::{self, Candidate},
        workspace::{self, Workspace},
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification,
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (DesignPlan, Molecule, CatalogSelection) {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let selection = catalog_digest::select(&catalog, &["nhei".into()]).unwrap();
    let molecule =
        Molecule::dna("original NheI circle", "AAAAGCTAGCCCC", Topology::Circular).unwrap();
    let target = Molecule::dna(
        "independent cut reference",
        "CTAGCCCCAAAAG",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "pinned NheI cut".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: molecule.revision(),
            molecule,
        }],
        steps: vec![Operation::CatalogDigest {
            template: Source::Input("source".into()),
            selection: selection.clone(),
            fragment_index: 0,
        }],
    };
    (plan, target, selection)
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
#[test]
fn catalog_opcode_and_fresh_report_require_their_own_semantics_version() {
    let (plan, target, _) = fixture();
    assert_eq!(postcard::to_allocvec(&plan.steps[0]).unwrap()[0], 17);
    let candidate = Candidate::new(plan.clone());
    let bytes = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    transfer::load_candidate(&bytes).unwrap();
    assert!(transfer::load_candidate(&wire(b"CGP\0", 10, &candidate)).is_err());
    let request = Request {
        id: 170,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan,
            target: target.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 36, &request)).is_err());
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let AnalysisResult::Verification(report) = response.result else {
        panic!("expected reconstruction")
    };
    assert_eq!(report.engine_version, verification::ENGINE_VERSION);
    assert_eq!(
        report.status,
        verification::evidence::Status::ConsistentUnderAssumptions
    );
    assert_eq!(report.products[0].sequence.as_str(), "CTAGCCCCAAAAG");
    let response = Response {
        id: 170,
        revision: target.revision(),
        result: AnalysisResult::Verification(report),
    };
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 36, &response)).is_err());
}
#[test]
fn catalog_preparation_in_existing_recipe_jobs_cannot_claim_old_headers() {
    use collagene_core::cloning::assembly::recipe::{
        self as assembly_recipe, FragmentSpec, Preparation, Recipe,
    };
    let (plan, target, selection) = fixture();
    let recipe = Recipe {
        name: "one cut overlap input".into(),
        fragments: vec![FragmentSpec {
            generation: None,
            role: "backbone".into(),
            input_id: "source".into(),
            preparation: Preparation::CatalogDigest {
                selection,
                fragment_index: 0,
                repair: None,
            },
            reverse: false,
        }],
        overlaps: vec![],
        circular: false,
        options: Default::default(),
    };
    assert_eq!(
        postcard::to_allocvec(&recipe.fragments[0].preparation).unwrap()[0],
        4
    );
    let revision = assembly_recipe::preview_revision(&plan.inputs, &target, &recipe).unwrap();
    let request = Request {
        id: 171,
        revision,
        job: AnalysisJob::PreviewOverlapRecipe {
            inputs: plan.inputs,
            target,
            recipe,
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 36, &request)).is_err());
    let response = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    protocol::response::Response::decode(&response).unwrap();
}
#[test]
fn archives_and_nested_workspace_jobs_keep_catalog_operation_version_guards() {
    use std::io::{Cursor, Read, Write};
    let (plan, target, _) = fixture();
    let mut frame = Workspace::from_documents(
        "pinned cloning",
        vec![plan.inputs[0].molecule.clone(), target],
        1,
    )
    .unwrap();
    frame.candidates.push(Candidate::new(plan));
    let archive = workspace::save_workspace(&frame).unwrap();
    let request = Request {
        id: 172,
        revision: blake3::hash(&archive).to_hex().to_string(),
        job: AnalysisJob::LoadWorkspace {
            bytes: archive.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 36, &request)).is_err());
    let mut input = zip::ZipArchive::new(Cursor::new(archive)).unwrap();
    let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for i in 0..input.len() {
        let mut entry = input.by_index(i).unwrap();
        let name = entry.name().to_owned();
        let mut bytes = Vec::new();
        entry.read_to_end(&mut bytes).unwrap();
        if name == "manifest.json" {
            let mut value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(value["plan_schema"], workspace::PLAN_SCHEMA_VERSION);
            value["plan_schema"] = 8.into();
            bytes = serde_json::to_vec(&value).unwrap();
        }
        output
            .start_file(name, zip::write::SimpleFileOptions::default())
            .unwrap();
        output.write_all(&bytes).unwrap();
    }
    assert!(workspace::load_workspace(&output.finish().unwrap().into_inner()).is_err());
}
