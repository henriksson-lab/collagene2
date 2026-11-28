//! New physical homology material stays explicit across worker, archive and saved-run boundaries.
use collagene_core::{
    chemistry::Terminus,
    cloning::{
        assembly::{
            junctions::{self as assembly_junctions, Method},
            recipe::{self as assembly_recipe, FragmentSpec, Preparation, RunPreparation},
        },
        ligation::recipe as ligation_recipe,
        overlap::assembly::{Options, OuterEndPolicy, Resection},
        restriction::catalog as catalog_digest,
    },
    documents::{
        transfer::{self, Candidate},
        workspace::{self, Workspace},
    },
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::Selection,
    },
    planning::plan::{Operation, PlanInput, Source},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const A: &str = "AGTCCGATGCTA";
const B: &str = "GCTATTCGAATTCAGTC";
const CIRCLE: &str = "AGTCCGATGCTATTCGAATTC";
const CUT: &str = "AATTCAGTCCGATGCTATTCG";
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
        overlaps: vec![4, 4],
        topology: Topology::Circular,
        resection: Resection::ExactHomology,
        outer_end_policy: OuterEndPolicy::JoiningEndsOnly,
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
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved() -> cloning_outputs::Commit {
    let (inputs, target, recipe) = fixture();
    let compiled = assembly_recipe::compile_homology(&inputs, &target, &recipe, &policy()).unwrap();
    let w = Workspace::from_documents(
        "Native homology origins",
        vec![
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target,
        ],
        2,
    )
    .unwrap();
    cloning_outputs::finalize(
        &w,
        &compiled.plan,
        binding(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: binding(&w, 0),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: binding(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Saved actual sealed circle".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
}
#[test]
fn worker86_job167_operation24_evidence3_and_portable_schemas_keep_explicit_policy_binding() {
    let (inputs, target, recipe) = fixture();
    let context =
        assembly_recipe::homology_preview_revision(&inputs, &target, &recipe, &policy()).unwrap();
    assert_ne!(
        context,
        assembly_recipe::preview_revision(&inputs, &target, &recipe).unwrap()
    );
    let request = Request {
        id: 167,
        revision: context.clone(),
        job: AnalysisJob::PreviewHomologyRecipe {
            inputs: inputs.clone(),
            target: target.clone(),
            recipe: recipe.clone(),
            policy: policy(),
        },
    };
    let encoded = protocol::request::Request::encode(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 162);
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 85, &request)).is_err());
    let bytes = protocol::execution::process_request(&encoded).unwrap();
    let response = protocol::response::Response::decode(&bytes).unwrap();
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 85, &response)).is_err());
    let AnalysisResult::OverlapRecipePreview(preview) = response.result else {
        panic!("physical homology preview");
    };
    preview.validate_binding(&context).unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
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
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        CIRCLE
    );
    let op = preview.compiled.plan.steps.last().unwrap();
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(op).unwrap())
            .unwrap()
            .0,
        23
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
        3
    );
    let candidate = Candidate::new(preview.compiled.plan.clone());
    let current = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &current[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&wire(b"CGP\0", 19, &candidate)).is_err());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 18, &candidate)).is_err());
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
    assert!(matches!(junctions.method, Some(Method::HomologyAssembly)));
    assert_eq!(junctions.junctions.len(), 2);
    let mut mismatched = policy();
    mismatched.overlaps[0] = 3;
    assert!(
        assembly_recipe::homology_source_preflight(
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
fn saved_homology_replays_actual_material_for_pcr_and_digest_and_retains_historical_evidence() {
    let commit = saved();
    let record = run_lineage::load_run(&commit.workspace, &commit.receipt.run.run_id).unwrap();
    assert_eq!(record.plan_schema, workspace::PLAN_SCHEMA_VERSION);
    let mut old = record.clone();
    old.plan_schema = 16;
    assert!(run_lineage::save_run_record(&old).is_err());
    let roundtrip =
        workspace::load_workspace(&workspace::save_workspace(&commit.workspace).unwrap()).unwrap();
    assert!(run_lineage::load_run(&roundtrip, &record.run_id).is_ok());
    let mut missing = roundtrip.clone();
    missing.attachments.retain(|_, attachment| {
        attachment.kind != collagene_core::material::physical_evidence::KIND
    });
    assert!(missing.validate().is_err());
    let selection = Selection {
        run_id: record.run_id,
        step_index: 2,
        event: commit.receipt.construct.clone(),
    };
    let origin = collagene_core::material::capsule::capsule(&roundtrip, &selection).unwrap();
    let source = roundtrip.revisions[&selection.event.revision].clone();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|enzyme| enzyme.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let pin = catalog_digest::select(&catalog, &[id]).unwrap();
    for (preparation, expected) in [
        (
            RunPreparation::ExactPcr {
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "GAATTCGA", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            CIRCLE,
        ),
        (
            RunPreparation::CatalogDigest {
                selection: pin,
                fragment_index: 0,
                repair: None,
            },
            CUT,
        ),
    ] {
        let inputs = vec![PlanInput {
            id: "saved".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }];
        let target = Molecule::dna(
            "Independent material continuation",
            expected,
            Topology::Linear,
        )
        .unwrap();
        let recipe = ligation_recipe::Recipe {
            name: "Actual homology product reuse".into(),
            fragments: vec![ligation_recipe::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: "saved".into(),
                    input_id: "saved".into(),
                    preparation: Preparation::FromRun {
                        origin: origin.clone(),
                        preparation,
                    },
                    reverse: false,
                },
                treatment: None,
            }],
            circular: false,
            allow_nicked: false,
        };
        let request = Request {
            id: 86,
            revision: ligation_recipe::preview_revision(
                &inputs,
                &target,
                &recipe,
                collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
            )
            .unwrap(),
            job: AnalysisJob::PreviewLigationRecipe {
                policy: collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
                inputs,
                target,
                recipe,
            },
        };
        assert!(protocol::execution::process_request(&wire(b"CGQ\0", 85, &request)).is_err());
        let result = protocol::response::Response::decode(
            &protocol::execution::process_request(
                &protocol::request::Request::encode(&request).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
        let AnalysisResult::LigationRecipePreview(preview) = result.result else {
            panic!("material continuation");
        };
        assert_eq!(
            preview.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            preview.verification.findings
        );
        assert_eq!(
            preview
                .verification
                .products
                .last()
                .unwrap()
                .sequence
                .as_str(),
            expected
        );
    }
    // The prior sequence-only profile path still cannot bootstrap PCR material.
    let (inputs, target, recipe) = fixture();
    let mut plan = assembly_recipe::compile(&inputs, &target, &recipe)
        .unwrap()
        .plan;
    let linear = Molecule::dna("Independent linear PCR", CIRCLE, Topology::Linear).unwrap();
    plan.target_revision = linear.revision();
    plan.steps.push(Operation::Pcr {
        template: Source::Step(2),
        forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("R", "GAATTCGA", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    });
    let report = verification::engine::verify(&plan, &linear);
    assert_ne!(report.status, Status::ConsistentUnderAssumptions);
    assert_eq!(report.products.len(), 3);
}
