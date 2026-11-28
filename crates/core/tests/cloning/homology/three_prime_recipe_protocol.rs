//! New physical homology material stays explicit across worker, archive and saved-run boundaries.
use collagene_core::{
    chemistry::Terminus,
    cloning::{
        assembly::{
            junctions::{self as assembly_junctions, Method},
            recipe::{self as assembly_recipe, FragmentSpec, Preparation},
        },
        overlap::{annealing::Options, resection::Cuts},
    },
    documents::{
        transfer::{self, Candidate},
        workspace::{self, Workspace},
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
const A: &str = "AGTCCGATGCTA";
const B: &str = "GCTATTCGAATTCAGTC";
const CIRCLE: &str = "AGTCCGATGCTATTCGAATTC";
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn policy() -> Options {
    Options {
        planned_rescue: None,
        overlaps: vec![4, 4],
        topology: Topology::Circular,
        cuts: vec![Cuts { left: 4, right: 4 }; 2],
        min_paired_bases: 4,
        max_comparisons: 100_000,
        max_candidates: 100,
    }
}
fn fixture() -> (Vec<PlanInput>, Molecule, assembly_recipe::Recipe) {
    (
        vec![input("a", A), input("b", B)],
        Molecule::dna("Independent sealed circle", CIRCLE, Topology::Circular).unwrap(),
        assembly_recipe::Recipe {
            standard: None,
            name: "Explicit physical homology closure".into(),
            fragments: ["a", "b"]
                .into_iter()
                .map(|id| FragmentSpec {
                    generation: None,
                    role: id.into(),
                    input_id: id.into(),
                    preparation: Preparation::DeclaredDuplex {
                        left: Terminus::PCR,
                        right: Terminus::PCR,
                    },
                    reverse: false,
                })
                .collect(),
            overlaps: vec![4, 4],
            circular: true,
            options: Default::default(),
        },
    )
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&payload).as_bytes());
    out.extend(payload);
    out
}
#[test]
fn worker88_job168_operation26_evidence5_and_portable_schemas_bind_annealing_without_fill() {
    let (inputs, target, recipe) = fixture();
    let context =
        assembly_recipe::annealing_preview_revision(&inputs, &target, &recipe, &policy()).unwrap();
    assert_ne!(
        context,
        assembly_recipe::preview_revision(&inputs, &target, &recipe).unwrap()
    );
    let request = Request {
        id: 168,
        revision: context.clone(),
        job: AnalysisJob::PreviewAnnealingRecipe {
            inputs: inputs.clone(),
            target: target.clone(),
            recipe: recipe.clone(),
            policy: policy(),
        },
    };
    let encoded = protocol::request::Request::encode(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 163);
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 87, &request)).is_err());
    let bytes = protocol::execution::process_request(&encoded).unwrap();
    let response = protocol::response::Response::decode(&bytes).unwrap();
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 87, &response)).is_err());
    let AnalysisResult::OverlapRecipePreview(preview) = response.result else {
        panic!("physical homology preview");
    };
    preview.validate_binding(&context).unwrap();
    assert_eq!(
        preview.verification.status,
        Status::Incomplete,
        "{:?}",
        preview.verification.findings
    );
    assert!(
        preview
            .verification
            .engine_version
            .eq(collagene_core::verification::ENGINE_VERSION)
    );
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        CIRCLE
    );
    let op = preview.compiled.plan.steps.last().unwrap();
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(op).unwrap())
            .unwrap()
            .0,
        25
    );
    let evidence = &preview
        .verification
        .physical_processing_reports
        .last()
        .unwrap()
        .evidence;
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(evidence).unwrap())
            .unwrap()
            .0,
        5
    );
    let candidate = Candidate::new(preview.compiled.plan.clone());
    let current = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &current[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&wire(b"CGP\0", 20, &candidate)).is_err());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 19, &candidate)).is_err());
    let mut w = Workspace::from_documents(
        "Portable physical candidate",
        vec![
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target.clone(),
        ],
        2,
    )
    .unwrap();
    w.candidates.push(candidate);
    let loaded = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    assert_eq!(loaded.candidates.len(), 1);
    let junctions =
        assembly_junctions::inspect_plan(&preview.compiled.plan, &target, &Default::default())
            .unwrap();
    assert!(matches!(junctions.method, Some(Method::HomologyAnnealing)));
    assert_eq!(junctions.junctions.len(), 2);
    let mut mismatched = policy();
    mismatched.overlaps[0] = 3;
    assert!(
        assembly_recipe::annealing_source_preflight(
            &inputs
                .iter()
                .map(|input| (input.id.as_str(), &input.molecule))
                .collect::<Vec<_>>(),
            &target,
            &recipe,
            &mismatched
        )
        .is_err()
    );
    let wrong = Request {
        revision: target.revision(),
        ..request
    };
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
}

#[test]
fn standalone_resection_uses_operation25_evidence4_and_never_hides_under_old_frames() {
    use collagene_core::planning::plan::{DesignPlan, Operation, Source};
    let target = Molecule::dna("Independent shortened top", "AGTCCGATG", Topology::Linear).unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Declared accessible source then shorten actual strands".into(),
        target_revision: target.revision(),
        inputs: vec![input("a", A)],
        steps: vec![
            Operation::DeclareSourceTermini {
                template: Source::Input("a".into()),
                left: Terminus::PCR,
                right: Terminus::PCR,
            },
            Operation::ResectThreePrime {
                source: Source::Step(0),
                cuts: Cuts { left: 2, right: 3 },
            },
        ],
    };
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&plan.steps[1]).unwrap())
            .unwrap()
            .0,
        24
    );
    let candidate = Candidate::new(plan.clone());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 19, &candidate)).is_err());
    let request = Request {
        id: 25,
        revision: target.revision(),
        job: AnalysisJob::Verify { plan, target },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 87, &request)).is_err());
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 87, &response)).is_err());
    let AnalysisResult::Verification(report) = response.result else {
        panic!("resection verification")
    };
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AGTCCGATG"
    );
    let evidence = &report.physical_processing_reports.last().unwrap().evidence;
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(evidence).unwrap())
            .unwrap()
            .0,
        4
    );
    let collagene_core::verification::evidence::PhysicalProcessingEvidence::ThreePrimeResection(
        evidence,
    ) = evidence
    else {
        panic!("actual resection")
    };
    assert_eq!(evidence.product.bottom.as_str(), "TAGCATCGGA");
    assert_eq!(evidence.product.bottom_start, 2);
}
