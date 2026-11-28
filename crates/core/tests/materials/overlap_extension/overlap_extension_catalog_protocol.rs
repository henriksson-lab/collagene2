//! The physical OE product, not the original stock, supplies catalog digest cuts.
use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation, RunPreparation},
        ligation::recipe as ligation_recipe,
        overlap::extension::{
            self as overlap_extension,
            digest::{self as overlap_extension_digest_inspection, Request as InspectionRequest},
            pcr::{PrimerPair, Request as Extension},
        },
        restriction::catalog as catalog_digest,
    },
    documents::{
        transfer::{self, Candidate},
        workspace::Workspace,
    },
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, Selection},
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
const AMPLICON: &str = "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT";
const FILLED: &str = "AATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATT";
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("F", f, ft).unwrap(),
        reverse: Primer::new("R", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
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
                generation: None,
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
fn execute(request: &Request) -> protocol::response::Response {
    protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(request).unwrap(),
        )
        .unwrap(),
    )
    .unwrap()
}
#[test]
fn extension_digest_preparation_preserves_physical_steps_and_current_encoding() {
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
            target: target.clone(),
            recipe,
        },
    };
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 83, &request)).is_err());
    let response = execute(&request);
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 83, &response)).is_err());
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
    assert!(
        preview
            .verification
            .engine_version
            .eq(collagene_core::verification::ENGINE_VERSION)
    );
    let candidate = Candidate::new(preview.compiled.plan);
    let encoded_candidate = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &encoded_candidate[4..6],
        &protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(transfer::load_candidate(&encoded_candidate).is_ok());
    // Unsupported candidate envelopes cannot relabel this current plan.
    let mut old_header = encoded_candidate;
    old_header[4..6].copy_from_slice(&18u16.to_le_bytes());
    assert!(transfer::load_candidate(&old_header).is_err());
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
fn saved() -> (
    Molecule,
    collagene_core::material::capsule::Capsule,
    Workspace,
) {
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
            template_strand: collagene_core::chemistry::Strand::Top,
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
    let capsule = collagene_core::material::capsule::capsule(&commit.workspace, &selected).unwrap();
    (
        commit.workspace.revisions[&event.revision].clone(),
        capsule,
        commit.workspace,
    )
}
#[test]
fn flat_saved_child8_and_actual_origin_digest_inspection_keep_closed_event_bindings() {
    let (source, origin, workspace) = saved();
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
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 83, &request)).is_err());
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
        overlap_extension_digest_inspection::source_preflight(&source, None, &request, None)
            .is_err()
    );
    let material = run_material::MaterialOrigin::Run(origin.clone());
    let preview = overlap_extension_digest_inspection::inspect(
        &source,
        Some(&material),
        &request,
        Some(&workspace),
    )
    .unwrap();
    preview
        .validate_binding(
            &overlap_extension_digest_inspection::inspection_revision(
                &source,
                Some(&material),
                &request,
                Some(&workspace),
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
        material.selection()
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
        overlap_extension_digest_inspection::inspection_revision(&source, None, &request, None)
            .unwrap();
    let job = Request {
        id: 165,
        revision: context.clone(),
        job: AnalysisJob::InspectOverlapExtensionDigest {
            source: source.clone(),
            origin: None,
            workspace: None,
            request: request.clone(),
        },
    };
    let encoded = protocol::request::Request::encode(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 160);
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 83, &job)).is_err());
    let reply = protocol::execution::process_request(&encoded).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
    assert_eq!(tag, 152);
    let response = protocol::response::Response::decode(&reply).unwrap();
    assert!(protocol::response::Response::decode(&wire(b"CGR\0", 83, &response)).is_err());
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
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&changed).unwrap()
        )
        .is_err()
    );
    let mut bounded = request;
    bounded.preparation.extension.max_comparisons = 18;
    let conditional =
        overlap_extension_digest_inspection::inspect(&source, None, &bounded, None).unwrap();
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
