//! Shared mismatch preparation is new; its existing physical operation stays compatible.
use collagene_core::{
    analysis::Primer,
    assembly_recipe::{self, FragmentSpec, Preparation, RunPreparation},
    cloning_outputs, ligation_recipe,
    mutagenesis::PcrOptions,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    tailing::TailingBase,
    tailing_recipe,
    transfer::{self, Candidate},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
const MUTATED: &str = "AGACCGATGCTACGTAGGACTTCAGACCTGGA";
fn mismatch() -> (Primer, Primer, PcrOptions) {
    let mut o = PcrOptions::default();
    o.binding.max_mismatches = 1;
    o.binding.min_three_prime_exact = 2;
    (
        Primer::new("F mismatch", "AGACCGAT", "").unwrap(),
        Primer::new("R mismatch", "TCCAGGTC", "").unwrap(),
        o,
    )
}
fn preparation() -> Preparation {
    let (forward, reverse, options) = mismatch();
    Preparation::MismatchPcr {
        forward,
        reverse,
        options,
    }
}
fn input(molecule: Molecule) -> PlanInput {
    PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fragment(preparation: Preparation) -> FragmentSpec {
    FragmentSpec {
        generation: None,
        role: "source".into(),
        input_id: "source".into(),
        preparation,
        reverse: false,
    }
}
fn wire(magic: &[u8; 4], v: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(v.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
fn old_rejected(r: &Request) {
    let error = protocol::execution::process_request(&wire(b"CGQ\0", 75, r)).unwrap_err();
    assert!(
        matches!(error, collagene_model::Error::Unsupported(_)),
        "{error}"
    );
}

// Explicit historical operation4 fixture; current recipe compilers intentionally emit operation29.

#[test]
fn mismatch_preparations_roundtrip_current_plans_and_reject_unsupported_versions() {
    assert!(protocol::schema::WORKER_SCHEMA_VERSION >= 76);
    assert_eq!(postcard::to_allocvec(&preparation()).unwrap()[0], 8);
    let inputs = vec![input(
        Molecule::dna("original", ORIGINAL, Topology::Linear).unwrap(),
    )];
    let target = Molecule::dna(
        "Independent mismatched requirement",
        MUTATED,
        Topology::Linear,
    )
    .unwrap();
    let fragments = vec![ligation_recipe::Fragment {
        source: fragment(preparation()),
        treatment: None,
    }];
    let recipe = ligation_recipe::Recipe {
        name: "Explicit mismatched source preparation".into(),
        fragments,
        circular: false,
        allow_nicked: false,
    };
    let request = Request {
        id: 1,
        revision: ligation_recipe::preview_revision(
            &inputs,
            &target,
            &recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            inputs: inputs.clone(),
            target: target.clone(),
            recipe,
        },
    };
    old_rejected(&request);
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let AnalysisResult::LigationRecipePreview(p) = &response.result else {
        panic!("ligation preview")
    };
    assert_eq!(p.verification.status, Status::ConsistentUnderAssumptions);
    assert_eq!(
        p.verification.products.last().unwrap().sequence.as_str(),
        MUTATED
    );
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 75, &response)).is_err());
    let current_candidate = Candidate::new(p.compiled.plan.clone());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 3, &current_candidate)).is_err());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 22, &current_candidate)).is_err());
    let current_verify = Request {
        id: 20,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: current_candidate.plan,
            target: target.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 75, &current_verify)).is_err());
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 91, &current_verify)).is_err());
    let candidate = Candidate::new(p.compiled.plan.clone());
    assert!(transfer::load_candidate(&wire(b"CGP\0", 3, &candidate)).is_err());
    let verify = Request {
        id: 2,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: candidate.plan,
            target: target.clone(),
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 75, &verify)).is_err());
    let overlap = assembly_recipe::Recipe {
        name: "Mismatch overlap preparation".into(),
        fragments: vec![fragment(preparation())],
        overlaps: vec![],
        circular: false,
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
    old_rejected(&request);
    let tail = tailing_recipe::Recipe {
        profile: collagene_core::tailing_recipe::Profile::Ta,

        name: "Mismatch before explicit tailing".into(),
        fragments: vec![tailing_recipe::Fragment {
            source: fragment(preparation()),
            action: collagene_core::tailing_recipe::Action::AddTail {
                base: TailingBase::A,
            },
            treatment: None,
        }],
        circular: false,
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
    old_rejected(&request);
}
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn saved_run_flat_mismatch_variant5_cannot_hide_under_existing_from_run_floor74_or75() {
    let source = Molecule::dna("Original exact source", ORIGINAL, Topology::Linear).unwrap();
    let target =
        Molecule::dna("Independent exact requirement", ORIGINAL, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Saved exact PCR".into(),
        target_revision: target.revision(),
        inputs: vec![input(source.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("exact F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("exact R", "TACAGGTC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Exact source events", vec![source, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved exact PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
    let origin = run_material::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let inputs = vec![input(commit.workspace.revisions[&event.revision].clone())];
    let target = Molecule::dna(
        "Independent mutation after saved PCR",
        MUTATED,
        Topology::Linear,
    )
    .unwrap();
    let (forward, reverse, options) = mismatch();
    let run = RunPreparation::MismatchPcr {
        forward,
        reverse,
        options,
    };
    assert_eq!(postcard::to_allocvec(&run).unwrap()[0], 5);
    let recipe = ligation_recipe::Recipe {
        name: "Saved source explicit mismatched PCR".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: fragment(Preparation::FromRun {
                origin,
                preparation: run,
            }),
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    };
    let request = Request {
        id: 5,
        revision: ligation_recipe::preview_revision(
            &inputs,
            &target,
            &recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            inputs,
            target,
            recipe,
        },
    };
    old_rejected(&request);
    let error = protocol::execution::process_request(&wire(b"CGQ\0", 74, &request)).unwrap_err();
    assert!(matches!(error, collagene_model::Error::Unsupported(_)));
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let AnalysisResult::LigationRecipePreview(p) = response.result else {
        panic!("saved source preview")
    };
    assert_eq!(p.verification.status, Status::ConsistentUnderAssumptions);
    assert!(
        p.compiled
            .plan
            .steps
            .iter()
            .any(|s| matches!(s, Operation::ReconstructRunProduct { .. }))
    );
    assert!(
        p.compiled
            .plan
            .steps
            .iter()
            .any(|s| matches!(s, Operation::PcrWithMismatches { .. }))
    );
}
