//! A branched physical preparation crosses existing recipe jobs without changing operations.
use collagene_core::{
    cloning::{
        assembly::recipe::{
            self as assembly_recipe, FragmentSpec, Preparation, Recipe, RunPreparation,
        },
        overlap::extension::{
            self as overlap_extension,
            pcr::{PrimerPair, Request as ExtensionRequest},
        },
    },
    documents::{
        transfer::{self, Candidate},
        workspace::Workspace,
    },
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::Selection,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const MUTANT: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const VECTOR: &str = "GATCGTCACAGTACGGTACCTAGC";
const TARGET: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCACAGTACGGTACCTAGC";
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("Actual forward", f, ft).unwrap(),
        reverse: Primer::new("Actual reverse", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn extension() -> ExtensionRequest {
    ExtensionRequest {
        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("AGTCCGAT", "", "TGACGATC", ""),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn recipe(preparation: Preparation) -> Recipe {
    Recipe {
        standard: None,
        name: "Mutate an insert before multipart assembly".into(),
        fragments: vec![
            FragmentSpec {
                generation: None,
                role: "edited_insert".into(),
                input_id: "insert".into(),
                preparation,
                reverse: false,
            },
            FragmentSpec {
                generation: None,
                role: "vector".into(),
                input_id: "vector".into(),
                preparation: Preparation::ExactPcr {
                    forward: Primer::new("Vector F", "GATCGTCA", "").unwrap(),
                    reverse: Primer::new("Vector R", "GCTAGGTA", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                reverse: false,
            },
        ],
        overlaps: vec![8],
        circular: false,
        options: Default::default(),
    }
}
fn fixture() -> (Vec<PlanInput>, Molecule, Recipe) {
    (
        vec![
            input(
                "insert",
                Molecule::dna("Original stock", ORIGINAL, Topology::Linear).unwrap(),
            ),
            input(
                "vector",
                Molecule::dna("Vector stock", VECTOR, Topology::Linear).unwrap(),
            ),
        ],
        Molecule::dna(
            "Independent edited multipart requirement",
            TARGET,
            Topology::Linear,
        )
        .unwrap(),
        recipe(Preparation::OverlapExtensionPcr {
            request: extension(),
        }),
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
fn request(inputs: Vec<PlanInput>, target: Molecule, recipe: Recipe) -> Request {
    Request {
        id: 82,
        revision: assembly_recipe::preview_revision(&inputs, &target, &recipe).unwrap(),
        job: AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        },
    }
}
#[test]
fn branched_extension_preparation_roundtrips_current_plan_and_independent_target() {
    let (inputs, target, recipe) = fixture();
    let sources_before = serde_json::to_value(&inputs).unwrap();
    assert_eq!(
        postcard::to_allocvec(&recipe.fragments[0].preparation).unwrap()[0],
        10
    );
    let req = request(inputs.clone(), target.clone(), recipe);
    let error = protocol::execution::process_request(&wire(b"CGQ\0", 81, &req)).unwrap_err();
    assert!(
        matches!(error, collagene_model::Error::Unsupported(_)),
        "{error}"
    );
    let bytes =
        protocol::execution::process_request(&protocol::request::Request::encode(&req).unwrap())
            .unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::response::Response::decode(&bytes).unwrap();
    let AnalysisResult::OverlapRecipePreview(preview) = &response.result else {
        panic!("overlap preview");
    };
    preview.validate().unwrap();
    assert!(preview.preparation_complete);
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.plan.steps.len(), 6);
    assert_eq!(preview.compiled.assembly_step, 5);
    assert_eq!(
        preview.prepared[0]
            .as_ref()
            .expect("prepared row must reconstruct")
            .sequence
            .as_str(),
        MUTANT
    );
    assert_eq!(
        preview.prepared[1]
            .as_ref()
            .expect("prepared row must reconstruct")
            .sequence
            .as_str(),
        VECTOR
    );
    assert_eq!(
        preview
            .verification
            .steps
            .iter()
            .map(|step| step
                .product
                .as_ref()
                .expect("each declared step must reconstruct"))
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![26, 30, 38, 38, 24, 54]
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
    assert_eq!(serde_json::to_value(inputs).unwrap(), sources_before);
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 81, &response)).is_err());
    assert!(
        preview
            .verification
            .engine_version
            .eq(collagene_core::verification::ENGINE_VERSION)
    );
    let candidate = Candidate::new(preview.compiled.plan.clone());
    let encoded_candidate = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &encoded_candidate[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&encoded_candidate).is_ok());
    let verify = Request {
        id: 81,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: candidate.plan.clone(),
            target: target.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 81, &verify)).is_err());
    for job in [
        AnalysisJob::ExportProtocolHtml {
            plan: candidate.plan.clone(),
            target: target.clone(),
        },
        AnalysisJob::ExportProtocolPdf {
            plan: candidate.plan,
            target,
        },
    ] {
        let req = Request {
            id: 82,
            revision: verify.revision.clone(),
            job,
        };
        let result = protocol::response::Response::decode(
            &protocol::execution::process_request(
                &protocol::request::Request::encode(&req).unwrap(),
            )
            .unwrap(),
        )
        .unwrap()
        .result;
        match result {
            AnalysisResult::ProtocolHtml(html) => {
                assert!(html.contains("Reciprocal extension") || html.contains("reanneal"))
            }
            AnalysisResult::ProtocolPdf(pdf) => assert!(pdf.bytes.starts_with(b"%PDF-")),
            other => panic!("unexpected protocol result {other:?}"),
        }
    }
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn flat_saved_run_preparation_variant7_cannot_hide_under_old_from_run_envelope() {
    let (mut inputs, target, _) = fixture();
    let source = inputs[0].molecule.clone();
    let original_target = Molecule::dna(
        "Original independent stock requirement",
        ORIGINAL,
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Prepare original stock".into(),
        target_revision: original_target.revision(),
        inputs: vec![input("insert", source.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("insert".into()),
            forward: Primer::new("Original F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Original R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let workspace =
        Workspace::from_documents("Exact stock events", vec![source, original_target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &workspace,
        &plan,
        binding(&workspace, 1),
        vec![InputEventBinding {
            input_id: "insert".into(),
            event: binding(&workspace, 0),
        }],
        &cloning_outputs::Options {
            name: "Prepared original stock".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selected = Selection {
        run_id: commit.receipt.run.run_id,
        step_index: 0,
        event: commit.receipt.construct,
    };
    let origin = collagene_core::material::capsule::capsule(&commit.workspace, &selected).unwrap();
    inputs[0] = input(
        "insert",
        commit.workspace.revisions[&selected.event.revision].clone(),
    );
    let child = RunPreparation::OverlapExtensionPcr {
        request: extension(),
    };
    assert_eq!(postcard::to_allocvec(&child).unwrap()[0], 7);
    let req = request(
        inputs,
        target,
        recipe(Preparation::FromRun {
            origin,
            preparation: child,
        }),
    );
    for version in [74, 81] {
        let error =
            protocol::execution::process_request(&wire(b"CGQ\0", version, &req)).unwrap_err();
        assert!(
            matches!(error, collagene_model::Error::Unsupported(_)),
            "{error}"
        );
    }
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(&protocol::request::Request::encode(&req).unwrap())
            .unwrap(),
    )
    .unwrap();
    let AnalysisResult::OverlapRecipePreview(preview) = response.result else {
        panic!("saved material preview");
    };
    preview.validate().unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.plan.steps.len(), 7);
    assert!(matches!(
        preview.compiled.plan.steps[0],
        Operation::ReconstructRunProduct { .. }
    ));
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        TARGET
    );
}
