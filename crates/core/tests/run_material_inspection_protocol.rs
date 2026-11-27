//! Origin discovery is documentary; explicit preparation binds one exact retained run.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation, RunPreparation},
    cloning_outputs::{self, Commit},
    ligation_recipe,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{Workspace, WorkspaceAttachment},
};
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn original() -> Commit {
    let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent requirement",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("R", "GTTACCCC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Recorded PCR", vec![source, target], 1).unwrap();
    cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
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

// Explicit historical operation4 fixture; current recipe compilers intentionally emit operation29.

#[test]
fn discovery_keeps74_tags_but_current_retained_workspace_requires92_and_exact_context() {
    assert!(protocol::WORKER_SCHEMA_VERSION >= 74);
    let commit = original();
    let binding = commit.receipt.construct.clone();
    let selection = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: binding.clone(),
    };
    for selected in [None, Some(selection.clone())] {
        let r = Request {
            id: 158,
            revision: run_material::inspection_revision(&commit.workspace, &binding, &selected)
                .unwrap(),
            job: AnalysisJob::InspectRunMaterial {
                workspace: commit.workspace.clone(),
                event: binding.clone(),
                selection: selected.clone(),
            },
        };
        let bytes = protocol::encode_request(&r).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 153);
        assert!(protocol::process_request(&wire(b"CGQ\0", 73, &r)).is_err());
        // The whole retained workspace was produced now: its native run plan
        // schema21/science26 is not a historical worker74 fixture.
        assert!(protocol::process_request(&wire(b"CGQ\0", 74, &r)).is_err());
        assert!(protocol::process_request(&wire(b"CGQ\0", 91, &r)).is_err());
        let reply = protocol::process_request(&bytes).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(tag, 145);
        let response = protocol::decode_response(&reply).unwrap();
        let AnalysisResult::RunMaterialInspection(report) = &response.result else {
            panic!("inspection")
        };
        report.validate_binding(&r.revision).unwrap();
        assert_eq!(report.event, binding);
        assert_eq!(report.origins, vec![selection.clone()]);
        assert_eq!(report.selected.is_some(), selected.is_some());
        assert!(protocol::decode_response(&wire(b"CGR\0", 73, &response)).is_err());
        let mut stale = r;
        stale.revision = "0".repeat(64);
        assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
        let mut stale = response;
        stale.revision = "0".repeat(64);
        assert!(
            protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &stale))
                .is_err()
        );
    }
}
#[test]
fn saved_run_recipe_preserves_exact_reconstruction_in_current_encoding() {
    let commit = original();
    let binding = commit.receipt.construct.clone();
    let selection = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: binding.clone(),
    };
    let origin = run_material::capsule(&commit.workspace, &selection).unwrap();
    let molecule = commit.workspace.revisions[&binding.revision].clone();
    let inputs = vec![PlanInput {
        id: "saved".into(),
        revision: binding.revision,
        molecule,
    }];
    let target = Molecule::dna(
        "Independent continuation",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let recipe = ligation_recipe::Recipe {
        name: "Explicit retained run source".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "saved".into(),
                input_id: "saved".into(),
                preparation: Preparation::FromRun {
                    origin,
                    preparation: RunPreparation::WholeLinear,
                },
                reverse: false,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    };
    let r = Request {
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
            inputs,
            target: target.clone(),
            recipe,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 73, &r)).is_err());
    let response = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&r).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::LigationRecipePreview(p) = &response.result else {
        panic!("preview")
    };
    assert!(protocol::decode_response(&wire(b"CGR\0", 73, &response)).is_err());
    let verify = Request {
        id: 2,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: p.compiled.plan.clone(),
            target,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 73, &verify)).is_err());
    assert!(protocol::process_request(&wire(b"CGQ\0", 91, &verify)).is_err());
}
#[test]
fn saved_run_actual_pcr_digest_uses_origin_bound_job159_not_plain_reference_inspector() {
    use collagene_core::{
        catalog_digest,
        pcr_digest_inspection::{self, Request as DigestRequest},
        restriction_inspection::InspectionOptions,
    };
    let commit = original();
    let selected = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let origin = run_material::capsule(&commit.workspace, &selected).unwrap();
    let molecule = commit.workspace.revisions[&selected.event.revision].clone();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let request = DigestRequest {
        source_role: "savedPCR".into(),
        forward: Primer::new("F", "ACGACCCC", "ATGCGAATTC").unwrap(),
        reverse: Primer::new("R", "GTTACCCC", "ATGCGAATTC").unwrap(),
        selection: catalog_digest::select(&catalog, &[id]).unwrap(),
        options: InspectionOptions::default(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    assert!(pcr_digest_inspection::source_preflight(&molecule, &request).is_err());
    let revision = pcr_digest_inspection::inspection_revision_from_material(
        &molecule,
        &run_material::MaterialOrigin::Run(origin.clone()),
        &request,
        &commit.workspace,
    )
    .unwrap();
    let r = Request {
        id: 159,
        revision: revision.clone(),
        job: AnalysisJob::InspectMaterialPcrDigest {
            molecule,
            origin: run_material::MaterialOrigin::Run(origin),
            request,
            workspace: commit.workspace.clone(),
        },
    };
    let bytes = protocol::encode_request(&r).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 154);
    assert!(protocol::process_request(&wire(b"CGQ\0", 73, &r)).is_err());
    let reply = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
    assert_eq!(tag, 146);
    let response = protocol::decode_response(&reply).unwrap();
    let AnalysisResult::MaterialPcrDigestInspection(report) = &response.result else {
        panic!("origin-bound PCR digest")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(report.origin_revision.len(), 64);
    assert_eq!(
        report.preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "ATGCGAATTCACGACCCCGGGGTAACGAATTCGCAT"
    );
    let inspection = report.preview.inspection.as_ref().unwrap();
    assert!(inspection.products_determined);
    assert_eq!(
        inspection
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 22, 9]
    );
    assert!(protocol::decode_response(&wire(b"CGR\0", 73, &response)).is_err());
    let mut stale = r;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let mut changed = response;
    changed.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &changed))
            .is_err()
    );
}
