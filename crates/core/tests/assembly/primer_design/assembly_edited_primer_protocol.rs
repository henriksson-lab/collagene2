//! Whole edited assembly primer proposals bind the independent intent before science.
use collagene_core::{
    cloning::{
        assembly::{
            edit_intent::{self as assembly_edit_intent, Edit},
            edited_primers::{self as assembly_edited_primer_design, Request as DesignRequest},
            intent::{self as assembly_intent, Contribution, Selector},
            junction_edit as assembly_junction_edit,
            recipe::Preparation,
            substitution_intent as assembly_substitution_intent,
        },
        overlap::extension::primers::Options,
        pcr::primers::DesignOptions,
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
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const VECTOR: &str = "GCTTACCGATGACCTAGTACGGCATTCGAAGC";
const TARGET: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAATATGCTTACCGATGACCTAGTACGGCATTCGAAGC";
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture() -> (Vec<PlanInput>, DesignRequest) {
    (
        vec![
            input(
                "original",
                Molecule::dna("Original source", ORIGINAL, Topology::Linear).unwrap(),
            ),
            input(
                "vector",
                Molecule::dna("Vector source", VECTOR, Topology::Linear).unwrap(),
            ),
        ],
        DesignRequest {
            intent: assembly_edit_intent::Request {
                base: assembly_substitution_intent::Request {
                    base: assembly_intent::Request {
                        name: "Independent mutation and literal junction".into(),
                        fragments: vec![
                            Contribution {
                                role: "insert".into(),
                                input_id: "original".into(),
                                selector: Selector::ReferenceRange {
                                    start: 0,
                                    end: 40,
                                    reverse: false,
                                },
                            },
                            Contribution {
                                role: "vector".into(),
                                input_id: "vector".into(),
                                selector: Selector::ReferenceRange {
                                    start: 0,
                                    end: 32,
                                    reverse: false,
                                },
                            },
                        ],
                        circular: false,
                        material_bindings: vec![],
                    },
                    junction_edits: vec![assembly_junction_edit::Edit {
                        junction: 0,
                        remove_left: 0,
                        remove_right: 0,
                        inserted: "ATAT".into(),
                    }],
                    substitutions: vec![],
                },
                edits: vec![Edit {
                    role: "insert".into(),
                    start: 16,
                    end: 20,
                    expected: "GGAC".into(),
                    replacement: "TT".into(),
                }],
            },
            overlap: 8,
            options: Options {
                homology_flank: 8,
                min_paired_bases: 8,
                primer_options: DesignOptions {
                    min_length: 8,
                    max_length: 8,
                    min_gc_fraction: 0.,
                    max_gc_fraction: 1.,
                    max_tm_difference_c: 100.,
                    max_homopolymer: 200,
                    max_three_prime_complementarity: 200,
                    ..Default::default()
                },
                ..Default::default()
            },
            profile_options: Default::default(),
            finishing: collagene_core::cloning::assembly::recipe::OverlapFinishing::ReferenceOnly,
        },
    )
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
fn execute(
    inputs: Vec<PlanInput>,
    request: DesignRequest,
    workspace: Option<Workspace>,
) -> protocol::response::Response {
    let revision =
        assembly_edited_primer_design::design_revision(&inputs, &request, workspace.as_ref())
            .unwrap();
    let r = Request {
        id: 166,
        revision,
        job: AnalysisJob::DesignEditedOverlapPrimers {
            inputs,
            request,
            workspace,
        },
    };
    protocol::response::Response::decode(
        &protocol::execution::process_request(&protocol::request::Request::encode(&r).unwrap())
            .unwrap(),
    )
    .unwrap()
}
#[test]
fn appended166_158_whole_design_preserves_literal_junction_and_internal_edit_against_independent_target()
 {
    let (inputs, request) = fixture();
    let before = serde_json::to_value(&inputs).unwrap();
    let context = assembly_edited_primer_design::design_revision(&inputs, &request, None).unwrap();
    let r = Request {
        id: 166,
        revision: context.clone(),
        job: AnalysisJob::DesignEditedOverlapPrimers {
            inputs: inputs.clone(),
            request: request.clone(),
            workspace: None,
        },
    };
    let bytes = protocol::request::Request::encode(&r).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 161);
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 84, &r)).is_err());
    let reply = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
    assert_eq!(tag, 153);
    let response = protocol::response::Response::decode(&reply).unwrap();
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 84, &response)).is_err());
    let AnalysisResult::EditedOverlapPrimersDesigned(report) = response.result else {
        panic!("whole edited primer design");
    };
    report.validate_binding(&context).unwrap();
    assert_eq!(report.intent.target.sequence.as_str(), TARGET);
    assert_eq!(report.rows.len(), 2);
    assert_eq!(report.rows[0].designs.len(), 3);
    assert_eq!(report.rows[1].designs.len(), 1);
    let Preparation::OverlapExtensionPcr { request: prepared } =
        &report.recipe.fragments[0].preparation
    else {
        panic!("edited row needs physical OE");
    };
    assert_eq!(
        prepared.outer.reverse.tail, "CGGTAAGCATAT",
        "outgoing tail is literalATAT plus actual next-core homology, in reverse oligo orientation"
    );
    assert!(matches!(
        report.recipe.fragments[1].preparation,
        Preparation::ExactPcr { .. }
    ));
    assert_eq!(report.recipe.overlaps, vec![8]);
    assert_eq!(
        report.preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.preview.verification.findings
    );
    assert_eq!(
        report
            .preview
            .verification
            .steps
            .iter()
            .map(|step| step
                .product
                .as_ref()
                .expect("each declared step must reconstruct"))
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![26, 30, 38, 50, 32, 74]
    );
    assert_eq!(
        report
            .preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        TARGET
    );
    let candidate = Candidate::new(report.preview.compiled.plan.clone());
    let encoded_candidate = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &encoded_candidate[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&encoded_candidate).is_ok());
    assert_eq!(
        report.preview.verification.engine_version,
        verification::ENGINE_VERSION
    );
    assert_eq!(
        verification::engine::verify(&candidate.plan, &report.intent.target).status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
}
#[test]
fn whole_context_binds_options_sources_and_intent_and_reverse_preparation_stays_in_original_orientation()
 {
    let (inputs, request) = fixture();
    let context = assembly_edited_primer_design::design_revision(&inputs, &request, None).unwrap();
    let mut changed_source = inputs.clone();
    changed_source[0]
        .molecule
        .name
        .push_str(" new immutable source");
    changed_source[0].revision = changed_source[0].molecule.revision();
    let mut changed_options = request.clone();
    changed_options.options.max_candidates = 99;
    let mut changed_intent = request.clone();
    changed_intent.intent.base.junction_edits[0].inserted = "CGCG".into();
    for (inputs, request) in [
        (changed_source, request.clone()),
        (inputs.clone(), changed_options),
        (inputs.clone(), changed_intent),
    ] {
        let r = Request {
            id: 166,
            revision: context.clone(),
            job: AnalysisJob::DesignEditedOverlapPrimers {
                inputs,
                request,
                workspace: None,
            },
        };
        assert!(
            protocol::execution::process_request(&protocol::request::Request::encode(&r).unwrap())
                .is_err()
        );
    }
    let AnalysisResult::EditedOverlapPrimersDesigned(mut report) =
        execute(inputs.clone(), request.clone(), None).result
    else {
        panic!("design");
    };
    report.recipe.overlaps[0] = 7;
    assert!(report.validate_binding(&context).is_err());
    let mut reverse_inputs = inputs;
    reverse_inputs[0] = input(
        "original",
        Molecule::dna(
            "Original reverse stock",
            "TGACGATCTACAGGTCTGAAGTCCTACGTAGCATCGGACT",
            Topology::Linear,
        )
        .unwrap(),
    );
    let mut reverse = request;
    if let Selector::ReferenceRange { reverse, .. } =
        &mut reverse.intent.base.base.fragments[0].selector
    {
        *reverse = true;
    } else {
        panic!("reference selector required by fixture");
    }
    let AnalysisResult::EditedOverlapPrimersDesigned(report) =
        execute(reverse_inputs, reverse, None).result
    else {
        panic!("reverse whole design");
    };
    assert_eq!(report.intent.target.sequence.as_str(), TARGET);
    assert_eq!(
        report.preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.preview.verification.findings
    );
    assert!(report.recipe.fragments[0].reverse);
    let Preparation::OverlapExtensionPcr { request } = &report.recipe.fragments[0].preparation
    else {
        panic!("reverse OE");
    };
    assert_eq!(request.outer.forward.tail, "CGGTAAGCATAT");
    assert!(request.outer.reverse.tail.is_empty());
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
fn per_role_saved_origin_is_retained_and_replayed_without_turning_its_metadata_into_stock_authority()
 {
    let (mut inputs, mut request) = fixture();
    let source = inputs[0].molecule.clone();
    let target = Molecule::dna(
        "Independent saved original requirement",
        ORIGINAL,
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Saved actual original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![input("original", source.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Saved role source", vec![source, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved PCR source".into(),
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
        "original",
        commit.workspace.revisions[&selected.event.revision].clone(),
    );
    request
        .intent
        .base
        .base
        .material_bindings
        .push(assembly_intent::MaterialBinding {
            input_id: "original".into(),
            event: origin.selection.event.clone(),
            origin: collagene_core::material::run::MaterialOrigin::Run(origin.clone()),
        });
    let AnalysisResult::EditedOverlapPrimersDesigned(report) =
        execute(inputs, request, Some(commit.workspace)).result
    else {
        panic!("saved role whole design");
    };
    report.validate_binding(&report.context_revision).unwrap();
    let Preparation::FromRun { origin: actual, .. } = &report.recipe.fragments[0].preparation
    else {
        panic!("saved role remains explicit");
    };
    assert_eq!(actual.selection, origin.selection);
    assert_eq!(report.intent.target.sequence.as_str(), TARGET);
    assert_eq!(
        report.preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.preview.verification.findings
    );
    assert!(matches!(
        report.preview.compiled.plan.steps[0],
        Operation::ReconstructRunProduct { .. }
    ));
}
