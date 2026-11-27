//! Direct saved material uses actual replayed strands rather than source metadata.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::FivePrime,
    cloning_outputs::{self, Commit},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    restriction_inspection::InspectionOptions,
    run_catalog_digest,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
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
    let source = Molecule::dna("Original", "ACGAGAATTCCGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent requirement",
        "ACGAGAATTCCGTAAC",
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
            forward: Primer::new("F", "ACGAGAAT", "").unwrap(),
            reverse: Primer::new("R", "GTTACGGA", "").unwrap(),
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
#[test]
fn actual_saved_pcr_catalog_inspection_roundtrips160_152_and_rejects_old_headers_and_stale_context()
{
    assert!(protocol::WORKER_SCHEMA_VERSION >= 75);
    let commit = original();
    let event = commit.receipt.construct.clone();
    let origin = run_material::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let molecule = commit.workspace.revisions[&event.revision].clone();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection = catalog_digest::select(&catalog, &[id]).unwrap();
    let options = InspectionOptions::default();
    let revision = run_catalog_digest::inspection_revision(
        &catalog, &selection, "saved", &molecule, &origin, &options,
    )
    .unwrap();
    let r = Request {
        id: 160,
        revision: revision.clone(),
        job: AnalysisJob::InspectRunCatalogRestriction {
            source_role: "saved".into(),
            molecule,
            origin,
            selection,
            options,
        },
    };
    let bytes = protocol::encode_request(&r).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 155);
    assert!(protocol::process_request(&wire(b"CGQ\0", 74, &r)).is_err());
    let reply = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
    assert_eq!(tag, 147);
    let response = protocol::decode_response(&reply).unwrap();
    let AnalysisResult::RunCatalogRestrictionInspection(report) = &response.result else {
        panic!("saved catalog inspection")
    };
    report.validate_binding(&revision).unwrap();
    let inner = &report.inspection.inspection;
    assert!(inner.scan_complete && inner.products_determined);
    assert_eq!(
        inner
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 11]
    );
    assert_eq!(report.fragment_modifications.len(), 2);
    assert_eq!(
        inner.fragments[0]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        FivePrime::Hydroxyl
    );
    assert_eq!(
        inner.fragments[1]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        FivePrime::Phosphate
    );
    assert!(protocol::decode_response(&wire(b"CGR\0", 74, &response)).is_err());
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
