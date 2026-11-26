//! The physical OE product, not the original stock, supplies catalog digest cuts.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{self, FragmentSpec, Preparation, RunPreparation},
    catalog_digest,
    chemistry::EndRepairMethod,
    cloning_outputs, ligation_recipe,
    overlap_extension::{self},
    overlap_extension_digest_inspection::{self, Request as InspectionRequest},
    overlap_extension_pcr::{PrimerPair, Request as Extension},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    transfer::{self, Candidate},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{Workspace, WorkspaceAttachment},
};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const AMPLICON: &str = "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT";
const FILLED: &str = "AATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATT";
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("F", f, ft).unwrap(),
        reverse: Primer::new("R", r, rt).unwrap(),
    }
}
fn extension() -> Extension {
    Extension {
        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("AGTCCGAT", "ATGCGAATTC", "TGACGATC", "ATGCGAATTC"),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn selection() -> catalog_digest::CatalogSelection {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|enzyme| enzyme.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    catalog_digest::select(&catalog, &[id]).unwrap()
}
fn source() -> Molecule {
    Molecule::dna("Original unchanged source", ORIGINAL, Topology::Linear).unwrap()
}
fn input(source: Molecule) -> PlanInput {
    PlanInput {
        id: "source".into(),
        revision: source.revision(),
        molecule: source,
    }
}
fn recipe(preparation: Preparation) -> ligation_recipe::Recipe {
    ligation_recipe::Recipe {
        name: "Edit, actually amplify, digest and fill".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                role: "insert".into(),
                input_id: "source".into(),
                preparation,
                reverse: false,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    }
}
fn compound() -> Preparation {
    Preparation::OverlapExtensionPcrCatalogDigest {
        request: extension(),
        selection: selection(),
        fragment_index: 1,
        repair: Some(EndRepairMethod::FillFivePrimeOverhangs),
    }
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
fn execute(request: &Request) -> protocol::Response {
    protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(request).unwrap()).unwrap(),
    )
    .unwrap()
}
#[test]
fn compound_preparation11_requires84_but_compiled_operations_remain_candidate18_plan16_science22() {
    assert_eq!(postcard::to_allocvec(&compound()).unwrap()[0], 11);
    let source = source();
    let before = source.revision();
    let inputs = vec![input(source.clone())];
    let target = Molecule::dna(
        "Independent filled digest fragment",
        FILLED,
        Topology::Linear,
    )
    .unwrap();
    let recipe = recipe(compound());
    let request = Request {
        id: 84,
        revision: ligation_recipe::preview_revision(&inputs, &target, &recipe).unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            inputs,
            target: target.clone(),
            recipe,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 83, &request)).is_err());
    let response = execute(&request);
    assert!(protocol::decode_response(&wire(b"CGR\0", 83, &response)).is_err());
    let AnalysisResult::LigationRecipePreview(preview) = response.result else {
        panic!("recipe preview");
    };
    preview.validate().unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.plan.steps.len(), 7);
    assert_eq!(preview.verification.products[3].sequence.as_str(), AMPLICON);
    assert_eq!(preview.verification.products[4].sequence.len(), 44);
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        FILLED
    );
    assert_eq!(preview.compiled.rows[0].preparation.repair_step, Some(5));
    assert!(matches!(
        preview.compiled.rows[0].preparation.prepared_source,
        Source::Step(5)
    ));
    assert!(preview.verification.engine_version.ends_with("/science-22"));
    let candidate = Candidate::new(preview.compiled.plan);
    let encoded_candidate = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &encoded_candidate[4..6],
        &protocol::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&encoded_candidate).is_ok());
    // These existing operations still retain their historical candidate18 floor.
    let mut historical_candidate = encoded_candidate;
    historical_candidate[4..6].copy_from_slice(&18u16.to_le_bytes());
    assert!(transfer::load_candidate(&historical_candidate).is_ok());
    let verify = Request {
        id: 81,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: candidate.plan,
            target,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 81, &verify)).is_ok());
    assert_eq!(source.revision(), before);
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved() -> (Molecule, run_material::Capsule) {
    let source = source();
    let target = Molecule::dna(
        "Original independent saved stock",
        ORIGINAL,
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual saved original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![input(source.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let w =
        Workspace::from_documents("Actual saved source events", vec![source, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Original prepared source".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct;
    let selected = Selection {
        run_id: commit.receipt.run.run_id,
        step_index: 0,
        event: event.clone(),
    };
    let capsule = run_material::capsule(&commit.workspace, &selected).unwrap();
    (commit.workspace.revisions[&event.revision].clone(), capsule)
}
#[test]
fn flat_saved_child8_and_actual_origin_digest_inspection_keep_closed_event_bindings() {
    let (source, origin) = saved();
    let child = RunPreparation::OverlapExtensionPcrCatalogDigest {
        request: extension(),
        selection: selection(),
        fragment_index: 1,
        repair: Some(EndRepairMethod::FillFivePrimeOverhangs),
    };
    assert_eq!(postcard::to_allocvec(&child).unwrap()[0], 8);
    let target = Molecule::dna(
        "Independent selected filled fragment",
        FILLED,
        Topology::Linear,
    )
    .unwrap();
    let inputs = vec![input(source.clone())];
    let recipe = recipe(Preparation::FromRun {
        origin: origin.clone(),
        preparation: child,
    });
    let request = Request {
        id: 84,
        revision: ligation_recipe::preview_revision(&inputs, &target, &recipe).unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 83, &request)).is_err());
    let AnalysisResult::LigationRecipePreview(preview) = execute(&request).result else {
        panic!("saved compound");
    };
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.plan.steps.len(), 8);
    let request = InspectionRequest {
        source_role: "insert".into(),
        preparation: extension(),
        selection: selection(),
        options: Default::default(),
    };
    assert!(
        overlap_extension_digest_inspection::source_preflight(&source, None, &request).is_err()
    );
    let preview =
        overlap_extension_digest_inspection::inspect(&source, Some(&origin), &request).unwrap();
    preview
        .validate_binding(
            &overlap_extension_digest_inspection::inspection_revision(
                &source,
                Some(&origin),
                &request,
            )
            .unwrap(),
        )
        .unwrap();
    assert!(preview.selectable());
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        AMPLICON
    );
    assert_eq!(
        preview.origin_evidence.as_ref().unwrap().selection,
        origin.selection
    );
}
#[test]
fn target_free165_157_inspection_scans_actual58bp_product_and_retains_incomplete_terminal_search() {
    let source = source();
    let request = InspectionRequest {
        source_role: "insert".into(),
        preparation: extension(),
        selection: selection(),
        options: Default::default(),
    };
    let context =
        overlap_extension_digest_inspection::inspection_revision(&source, None, &request).unwrap();
    let job = Request {
        id: 165,
        revision: context.clone(),
        job: AnalysisJob::InspectOverlapExtensionDigest {
            source: source.clone(),
            origin: None,
            request: request.clone(),
        },
    };
    let encoded = protocol::encode_request(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 165);
    assert!(protocol::process_request(&wire(b"CGQ\0", 83, &job)).is_err());
    let reply = protocol::process_request(&encoded).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
    assert_eq!(tag, 157);
    let response = protocol::decode_response(&reply).unwrap();
    assert!(protocol::decode_response(&wire(b"CGR\0", 83, &response)).is_err());
    let AnalysisResult::OverlapExtensionDigestInspection(preview) = response.result else {
        panic!("actual OE digest");
    };
    preview.validate_binding(&context).unwrap();
    assert!(!preview.incomplete);
    assert!(preview.selectable());
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        AMPLICON
    );
    assert_eq!(
        preview
            .pcr_stages
            .iter()
            .map(|s| s.product_length.unwrap())
            .collect::<Vec<_>>(),
        vec![26, 30, 58]
    );
    assert_eq!(
        preview
            .inspection
            .as_ref()
            .unwrap()
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 44, 9]
    );
    let mut changed = job.clone();
    changed.revision = source.revision();
    assert!(protocol::process_request(&protocol::encode_request(&changed).unwrap()).is_err());
    let mut bounded = request;
    bounded.preparation.extension.max_comparisons = 18;
    let conditional =
        overlap_extension_digest_inspection::inspect(&source, None, &bounded).unwrap();
    assert!(conditional.incomplete);
    assert!(
        conditional.selectable(),
        "the declared candidate remains selectable with explicit uncertainty"
    );
    assert!(
        !conditional
            .extension
            .as_ref()
            .unwrap()
            .terminal_pairings
            .complete
    );
}
#[test]
fn ordinary_draft9_payload_requires84_without_changing_old_draft8_floor() {
    let mut w = Workspace::from_documents("Ordinary draft", vec![source()], 0).unwrap();
    let bytes = b"{}".to_vec();
    w.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        WorkspaceAttachment {
            kind: "assembly-draft-v9".into(),
            name: "ordinary draft".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let request = Request {
        id: 84,
        revision: w.head(&w.active_document).unwrap().revision(),
        job: AnalysisJob::SaveWorkspace { workspace: w },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 83, &request)).is_err());
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_ok());
}
