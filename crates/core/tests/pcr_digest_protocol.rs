//! Actual primer-added sites are inspected without fabricating an independent target.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    pcr_digest_inspection::{self, Request as InspectionRequest},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    restriction_inspection::InspectionOptions,
};
fn fixture() -> (Molecule, InspectionRequest) {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    (
        Molecule::dna("Original, no EcoRI site", "AAACCCGGGTTT", Topology::Linear).unwrap(),
        InspectionRequest {
            source_role: "original".into(),
            forward: Primer::new("forward added site", "AAA", "ATGCGAATTC").unwrap(),
            reverse: Primer::new("reverse added site", "AAA", "ATGCGAATTC").unwrap(),
            selection: catalog_digest::select(&catalog, &[id]).unwrap(),
            options: InspectionOptions::default(),
        },
    )
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    b
}
#[test]
fn unique_actual_pcr_digest_roundtrips_new_tags_and_rejects_downgrade_and_stale_context() {
    assert!(protocol::WORKER_SCHEMA_VERSION >= 72);
    assert!(protocol::CANDIDATE_SCHEMA_VERSION >= 16);
    let (molecule, inspection) = fixture();
    let revision = pcr_digest_inspection::inspection_revision(&molecule, &inspection).unwrap();
    let request = Request {
        id: 157,
        revision: revision.clone(),
        job: AnalysisJob::InspectPcrDigest {
            molecule,
            request: inspection,
        },
    };
    let encoded = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 155);
    assert!(protocol::process_request(&wire(b"CGQ\0", 71, &request)).is_err());
    let response_bytes = protocol::process_request(&encoded).unwrap();
    let ((_, _, tag), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&response_bytes[42..]).unwrap();
    assert_eq!(tag, 147);
    let response = protocol::decode_response(&response_bytes).unwrap();
    let AnalysisResult::PcrDigestInspection(preview) = &response.result else {
        panic!("PCR digest inspection")
    };
    preview.validate_binding(&revision).unwrap();
    assert_eq!(preview.pcr_geometries.len(), 1);
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "ATGCGAATTCAAACCCGGGTTTGAATTCGCAT"
    );
    let inspection = preview.inspection.as_ref().unwrap();
    assert!(inspection.scan_complete && inspection.products_determined);
    assert_eq!(inspection.fragments.len(), 3);
    assert_eq!(
        inspection
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 18, 9]
    );
    assert!(protocol::decode_response(&wire(b"CGR\0", 71, &response)).is_err());
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let mut stale_reply: Response = response;
    stale_reply.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&wire(
            b"CGR\0",
            protocol::WORKER_SCHEMA_VERSION,
            &stale_reply
        ))
        .is_err()
    );
}
#[test]
fn ambiguous_actual_pcr_geometry_has_no_digest_products() {
    let (_, request) = fixture();
    let molecule = Molecule::dna(
        "Repeated primer sites",
        "AAACCCGGGTTTAAACCCGGGTTT",
        Topology::Linear,
    )
    .unwrap();
    let revision = pcr_digest_inspection::inspection_revision(&molecule, &request).unwrap();
    let r = Request {
        id: 2,
        revision,
        job: AnalysisJob::InspectPcrDigest { molecule, request },
    };
    let response = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&r).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::PcrDigestInspection(p) = response.result else {
        panic!("inspection")
    };
    assert!(p.pcr_geometries.len() > 1);
    assert!(p.amplicon.is_none());
    assert!(p.inspection.is_none());
}
