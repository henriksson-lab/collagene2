//! Exact role/source/geometry/options binding across the native and WASM worker seam.
use collagene_core::{
    Error, Molecule, Topology,
    analysis::{self, Enzyme},
    protocol::{self, AnalysisJob as J, AnalysisResult as R, Request, Response},
    restriction_inspection::*,
};
fn fixture() -> (Molecule, Vec<Enzyme>, InspectionOptions) {
    let m = Molecule::dna("digest", "AAAAGAATTCCCCCGGATCCAAAA", Topology::Linear).unwrap();
    let enzymes = analysis::starter_enzymes()
        .into_iter()
        .filter(|e| matches!(e.name.as_str(), "EcoRI" | "BamHI"))
        .collect();
    (m, enzymes, InspectionOptions::default())
}
fn raw(magic: &[u8; 4], version: u16, body: &[u8]) -> Vec<u8> {
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((body.len() as u32).to_le_bytes());
    out.extend(blake3::hash(body).as_bytes());
    out.extend(body);
    out
}
fn run(revision: String, job: J) -> Response {
    let bytes = protocol::encode_request(&Request {
        id: 901,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    assert_eq!(&bytes[4..6], &21u16.to_le_bytes());
    let response = protocol::decode_response(
        &protocol::AnalysisSession::default()
            .process_request(&bytes)
            .unwrap(),
    )
    .unwrap();
    assert_eq!(response.id, 901);
    assert_eq!(response.revision, revision);
    response
}
#[test]
fn inspection_and_fresh_product_preserve_exact_evidence_and_bind_every_input() {
    let (m, es, o) = fixture();
    let role = "source-event";
    let revision = batch_revision(role, &m, &es, &o).unwrap();
    let R::RestrictionInspection(report) = run(
        revision.clone(),
        J::InspectRestriction {
            source_role: role.into(),
            molecule: m.clone(),
            enzymes: es.clone(),
            options: o.clone(),
        },
    )
    .result
    else {
        panic!()
    };
    report.validate_bindings_for(role, &m, &es, &o).unwrap();
    assert!(report.scan_complete && report.products_determined);
    assert_eq!(
        report
            .sites
            .iter()
            .map(|s| s.recognition_start)
            .collect::<Vec<_>>(),
        vec![4, 14]
    );
    assert_eq!(
        report
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 10, 9]
    );
    let revision = product_revision(role, &m, &es, &o, 1).unwrap();
    let R::DigestProduct(p) = run(
        revision,
        J::MaterializeDigestFragment {
            source_role: role.into(),
            molecule: m.clone(),
            enzymes: es.clone(),
            options: o.clone(),
            fragment_index: 1,
        },
    )
    .result
    else {
        panic!()
    };
    p.validate_bindings_for(role, &m, &es, &o, 1).unwrap();
    assert_eq!(p.molecule.sequence.as_str(), "AATTCCCCCG");
    assert_eq!(p.summary.top_start, 5);
    assert_ne!(
        product_revision(role, &m, &es, &o, 0).unwrap(),
        p.batch_revision
    );
    assert_ne!(
        batch_revision("other-event", &m, &es, &o).unwrap(),
        report.batch_revision
    );
    let mut changed = es.clone();
    changed[0].top_cut += 1;
    assert!(
        report
            .validate_bindings_for(role, &m, &changed, &o)
            .is_err()
    );
    let mut changed = o.clone();
    changed.max_comparisons -= 1;
    assert!(
        report
            .validate_bindings_for(role, &m, &es, &changed)
            .is_err()
    );
    let mut reversed = es.clone();
    reversed.reverse();
    assert_ne!(
        batch_revision(role, &m, &reversed, &o).unwrap(),
        report.batch_revision
    );
    let request = Request {
        id: 901,
        revision: report.batch_revision.clone(),
        job: J::MaterializeDigestFragment {
            source_role: role.into(),
            molecule: m.clone(),
            enzymes: es.clone(),
            options: o.clone(),
            fragment_index: 1,
        },
    };
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
    let forged = Response {
        id: 901,
        revision: "0".repeat(64),
        result: R::RestrictionInspection(report),
    };
    assert!(
        protocol::decode_response(&raw(b"CGR\0", 20, &postcard::to_allocvec(&forged).unwrap()))
            .is_err()
    );
    let forged = Response {
        id: 901,
        revision: "0".repeat(64),
        result: R::DigestProduct(p),
    };
    assert!(
        protocol::decode_response(&raw(b"CGR\0", 20, &postcard::to_allocvec(&forged).unwrap()))
            .is_err()
    );
}
#[test]
fn possible_edge_and_limited_scans_keep_evidence_without_product_authority() {
    let (_, es, mut o) = fixture();
    for (dna, enzyme) in [
        ("AAAAGAANTCAAAA", es[0].clone()),
        (
            "GGTCTCA",
            analysis::starter_enzymes()
                .into_iter()
                .find(|e| e.name == "BsaI")
                .unwrap(),
        ),
    ] {
        let m = Molecule::dna("uncertain or substrate edge", dna, Topology::Linear).unwrap();
        let enzymes = vec![enzyme];
        let R::RestrictionInspection(r) = run(
            batch_revision("s", &m, &enzymes, &o).unwrap(),
            J::InspectRestriction {
                source_role: "s".into(),
                molecule: m.clone(),
                enzymes: enzymes.clone(),
                options: o.clone(),
            },
        )
        .result
        else {
            panic!()
        };
        assert!(r.scan_complete);
        assert!(!r.products_determined);
        assert!(!r.sites.is_empty());
        assert!(!r.sites[0].diagnostics.is_empty());
        let result = run(
            product_revision("s", &m, &enzymes, &o, 0).unwrap(),
            J::MaterializeDigestFragment {
                source_role: "s".into(),
                molecule: m,
                enzymes,
                options: o.clone(),
                fragment_index: 0,
            },
        )
        .result;
        assert!(matches!(result, R::Error(_)));
    }
    let (m, es, _) = fixture();
    o.max_comparisons = 1;
    let R::RestrictionInspection(r) = run(
        batch_revision("s", &m, &es, &o).unwrap(),
        J::InspectRestriction {
            source_role: "s".into(),
            molecule: m,
            enzymes: es,
            options: o,
        },
    )
    .result
    else {
        panic!()
    };
    assert!(!r.scan_complete);
    assert!(!r.products_determined);
    assert!(r.limits.contains(&InspectionLimit::Comparisons));
}
#[test]
fn appended_tags_require_twenty_before_body_decode_and_inspection_keeps_ordinary_cap() {
    for (magic, tags) in [
        (b"CGQ\0", vec![67u32, 68, 69]),
        (b"CGR\0", vec![60u32, 61, 62]),
    ] {
        for tag in tags {
            let body = postcard::to_allocvec(&(902u32, "binding", tag)).unwrap();
            for version in 1..20 {
                let bytes = raw(magic, version, &body);
                let e = if magic == b"CGQ\0" {
                    protocol::process_request(&bytes).unwrap_err()
                } else {
                    protocol::decode_response(&bytes).unwrap_err()
                };
                assert!(e.to_string().contains("schema 20"), "{e}");
            }
            if matches!(tag, 67 | 68 | 60 | 61) {
                let mut large = body;
                large.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
                let bytes = raw(magic, 20, &large);
                let e = if magic == b"CGQ\0" {
                    protocol::process_request(&bytes).unwrap_err()
                } else {
                    protocol::decode_response(&bytes).unwrap_err()
                };
                assert!(matches!(e, Error::TooLarge));
            }
        }
    }
    let (mut m, es, o) = fixture();
    m.metadata
        .insert("oversized".into(), "x".repeat(32 * 1024 * 1024));
    assert!(matches!(
        batch_revision("s", &m, &es, &o),
        Err(Error::TooLarge)
    ));
}

#[test]
fn native_installation_reconstructs_product_and_preserves_exact_source_event() {
    use collagene_core::{
        restriction_history as history,
        run_lineage::EventBinding,
        workspace::{Workspace, load_workspace},
    };
    let (m, es, o) = fixture();
    let w = Workspace::from_documents("source workspace", vec![m], 0).unwrap();
    let source = EventBinding {
        document_id: w.documents[0].id.clone(),
        event_id: w.documents[0].head_event.clone().unwrap(),
        revision: w.documents[0].head.clone(),
    };
    let revision = history::installation_revision(&w, &source, &es, &o, 1).unwrap();
    let R::DigestInstalled(installed) = run(
        revision,
        J::InstallDigestFragment {
            workspace: w.clone(),
            source: source.clone(),
            enzymes: es.clone(),
            options: o.clone(),
            fragment_index: 1,
        },
    )
    .result
    else {
        panic!()
    };
    installed.validate().unwrap();
    let restored = load_workspace(&installed.archive).unwrap();
    installed.validate_archive(&restored).unwrap();
    assert_eq!(restored.documents.len(), 2);
    assert_eq!(restored.active_document, installed.document_id);
    assert_eq!(
        restored
            .head(&installed.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "AATTCCCCCG"
    );
    let event = restored
        .history
        .iter()
        .find(|e| e.id == installed.event_id)
        .unwrap();
    assert_eq!(event.parents, vec![source.event_id.clone()]);
    assert_eq!(
        restored.head(&source.document_id).unwrap().revision(),
        source.revision
    );
    let mut wrong = source;
    wrong.event_id = format!("event-{}", "0".repeat(64));
    assert!(history::installation_revision(&w, &wrong, &es, &o, 1).is_err());
    let forged = Response {
        id: 901,
        revision: "0".repeat(64),
        result: R::DigestInstalled(installed),
    };
    assert!(
        protocol::decode_response(&raw(b"CGR\0", 20, &postcard::to_allocvec(&forged).unwrap()))
            .is_err()
    );
}
