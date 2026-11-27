//! One target-free OE derivation through the current binary worker and archive.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    material_derivation as native, overlap_extension,
    overlap_extension_pcr::{PrimerPair, Request as Preparation},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    run_lineage::EventBinding,
    workspace::{self, Workspace},
};
const STOCK: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const PRODUCT: &str = "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT";
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("New forward", f, ft).unwrap(),
        reverse: Primer::new("New reverse", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn preparation() -> Preparation {
    Preparation {
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
fn fixture() -> (Workspace, EventBinding) {
    let w = Workspace::from_documents(
        "Target-free source",
        vec![Molecule::dna("Stock", STOCK, Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let d = &w.documents[0];
    let event = EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    };
    (w, event)
}
fn tag(v: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(v).unwrap())
        .unwrap()
        .0
}
fn roundtrip(job: AnalysisJob, context: String) -> Response {
    let q = Request {
        id: 77,
        revision: context.clone(),
        job,
    };
    let request = protocol::encode_request(&q).unwrap();
    assert_eq!(
        &request[4..6],
        &protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::AnalysisSession::default()
        .process_request(&request)
        .unwrap();
    let r = protocol::decode_response(&response).unwrap();
    assert_eq!(r.id, 77);
    assert_eq!(r.revision, context);
    r
}
fn wire(r: &Response) -> Vec<u8> {
    let payload = postcard::to_allocvec(r).unwrap();
    let mut out = b"CGR\0".to_vec();
    out.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&payload).as_bytes());
    out.extend(payload);
    out
}
#[test]
fn literal_outer_product_installs_without_a_target_then_replays_from_worker_archive() {
    let (w, event) = fixture();
    let original = workspace::save_workspace(&w).unwrap();
    let origin = native::SourceOrigin::OriginalReference;
    let prep = preparation();
    let context = native::overlap_extension_preview_revision(&w, &event, &origin, &prep).unwrap();
    let job = AnalysisJob::InspectMaterialOverlapExtension {
        workspace: w.clone(),
        source: event.clone(),
        origin: origin.clone(),
        preparation: prep.clone(),
    };
    assert_eq!(tag(&job), 188);
    let r = roundtrip(job, context.clone());
    assert_eq!(tag(&r.result), 178);
    let AnalysisResult::MaterialOverlapExtensionPreview(preview) = r.result else {
        panic!("OE preview")
    };
    preview.validate_binding(&context).unwrap();
    assert!(preview.reusable());
    assert_eq!(preview.product.as_ref().unwrap().sequence.as_str(), PRODUCT);
    let request = native::Request {
        source: event,
        origin,
        derivation: native::Derivation::OverlapExtensionPcr {
            preparation: prep,
            selection: preview.selection().unwrap(),
        },
        name: "Actual edited outer product".into(),
        save_primers: true,
    };
    let expected = native::revision(&w, &request).unwrap();
    let r = roundtrip(
        AnalysisJob::InstallMaterialDerivation {
            workspace: w.clone(),
            request,
        },
        expected,
    );
    let AnalysisResult::MaterialDerivationInstalled(installed) = r.result else {
        panic!("install")
    };
    let archive = workspace::load_workspace(&installed.archive).unwrap();
    installed.receipt.validate_archive(&archive).unwrap();
    let fresh = native::reconstruct(&archive, &installed.receipt.product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), PRODUCT);
    assert_eq!(fresh.tracked().duplex.top.as_str(), PRODUCT);
    assert!(!fresh.incomplete());
    assert!(installed.receipt.primer_attachment.is_some());
    assert_eq!(workspace::save_workspace(&w).unwrap(), original);
}
#[test]
fn stale_context_foreign_event_and_forged_reply_are_refused_before_material_installation() {
    let (w, event) = fixture();
    let origin = native::SourceOrigin::OriginalReference;
    let prep = preparation();
    let context = native::overlap_extension_preview_revision(&w, &event, &origin, &prep).unwrap();
    let mut changed = prep.clone();
    changed.outer.forward.tail = "AAAA".into();
    let q = Request {
        id: 1,
        revision: context.clone(),
        job: AnalysisJob::InspectMaterialOverlapExtension {
            workspace: w.clone(),
            source: event.clone(),
            origin: origin.clone(),
            preparation: changed,
        },
    };
    assert!(protocol::process_request(&protocol::encode_request(&q).unwrap()).is_err());
    let mut alias = event.clone();
    alias.event_id = "not-an-owned-event".into();
    let q = Request {
        id: 2,
        revision: context.clone(),
        job: AnalysisJob::InspectMaterialOverlapExtension {
            workspace: w.clone(),
            source: alias,
            origin: origin.clone(),
            preparation: prep.clone(),
        },
    };
    assert!(protocol::encode_request(&q).is_err());
    let preview = native::preview_overlap_extension(&w, &event, &origin, &prep).unwrap();
    let mut response = Response {
        id: 3,
        revision: context,
        result: AnalysisResult::MaterialOverlapExtensionPreview(preview),
    };
    let AnalysisResult::MaterialOverlapExtensionPreview(report) = &mut response.result else {
        unreachable!()
    };
    report.context_revision = "0".repeat(64);
    assert!(protocol::decode_response(&wire(&response)).is_err());
    let preview = native::preview_overlap_extension(&w, &event, &origin, &prep).unwrap();
    let mut selection = preview.selection().unwrap();
    selection.expected_product_revision = "0".repeat(64);
    let request = native::Request {
        source: event,
        origin,
        derivation: native::Derivation::OverlapExtensionPcr {
            preparation: prep,
            selection,
        },
        name: "Withheld altered selection".into(),
        save_primers: true,
    };
    assert!(native::prepare_install(&w, &request).is_err());
    assert_eq!(w.documents.len(), 1);
    assert!(w.attachments.is_empty());
}
