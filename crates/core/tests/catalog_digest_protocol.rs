//! Appended catalog wire jobs bind actual edition/IDs/geometry while legacy layouts stay readable.
use collagene_core::{
    Error, Molecule, Topology, analysis,
    catalog_digest::*,
    protocol::{self, AnalysisJob as J, AnalysisResult as R, Request, Response},
    restriction_history::SourceBinding,
    restriction_inspection::InspectionOptions,
    workspace::Workspace,
};
fn fixture() -> (
    catalog::Catalog,
    CatalogSelection,
    Molecule,
    InspectionOptions,
) {
    let c = builtin_catalog().unwrap();
    let s = select(&c, &["haeiii".into()]).unwrap();
    let m = Molecule::dna("literal HaeIII", "AAGGCCTT", Topology::Linear).unwrap();
    (c, s, m, InspectionOptions::default())
}
fn raw(magic: &[u8; 4], version: u16, p: &[u8]) -> Vec<u8> {
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(p).as_bytes());
    b.extend(p);
    b
}
fn run(revision: String, job: J) -> Response {
    let r = Request {
        id: 1001,
        revision: revision.clone(),
        job,
    };
    let b = protocol::encode_request(&r).unwrap();
    assert_eq!(
        &b[4..6],
        &collagene_core::protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let reply = protocol::decode_response(
        &protocol::AnalysisSession::default()
            .process_request(&b)
            .unwrap(),
    )
    .unwrap();
    assert_eq!(reply.id, r.id);
    assert_eq!(reply.revision, revision);
    reply
}
#[test]
fn exact_catalog_geometry_source_and_options_are_bound_before_fresh_products() {
    let (c, s, m, o) = fixture();
    let role = "event-role";
    let rev = batch_revision(&c, &s, role, &m, &o).unwrap();
    let R::CatalogRestrictionInspection(r) = run(
        rev,
        J::InspectCatalogRestriction {
            selection: s.clone(),
            source_role: role.into(),
            molecule: m.clone(),
            options: o.clone(),
        },
    )
    .result
    else {
        panic!()
    };
    r.validate_bindings_for(&c, &s, role, &m, &o).unwrap();
    assert_eq!(r.selection.enzymes[0].id, "haeiii");
    assert!(r.inspection.products_determined);
    assert_eq!(r.inspection.sites[0].recognition_start, 2);
    assert_eq!(r.inspection.sites[0].top_cut, Some(4));
    assert_eq!(
        r.inspection
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![4, 4]
    );
    assert_ne!(r.batch_revision, r.inspection.batch_revision);
    let R::CatalogDigestProduct(p) = run(
        product_revision(&c, &s, role, &m, &o, 0).unwrap(),
        J::MaterializeCatalogDigestFragment {
            selection: s.clone(),
            source_role: role.into(),
            molecule: m.clone(),
            options: o.clone(),
            fragment_index: 0,
        },
    )
    .result
    else {
        panic!()
    };
    p.validate_bindings_for(&c, &s, role, &m, &o, 0).unwrap();
    assert_eq!(p.product.molecule.sequence.as_str(), "AAGG");
    assert_eq!(
        p.product.molecule.metadata[PRODUCT_REVISION_KEY],
        p.batch_revision
    );
    assert_eq!(
        unverified_selection(&p.product.molecule).unwrap(),
        Some(s.clone())
    );
    assert_ne!(
        product_revision(&c, &s, role, &m, &o, 1).unwrap(),
        p.batch_revision
    );
    assert!(
        r.validate_bindings_for(&c, &s, "other-role", &m, &o)
            .is_err()
    );
    let mut changed = s.clone();
    changed.catalog_content_hash = "0".repeat(64);
    assert!(batch_revision(&c, &changed, role, &m, &o).is_err());
    changed = s.clone();
    changed.enzymes[0].enzyme.top_cut = 1;
    assert!(batch_revision(&c, &changed, role, &m, &o).is_err());
    let mut options = o.clone();
    options.max_sites -= 1;
    assert!(r.validate_bindings_for(&c, &s, role, &m, &options).is_err());
    let pair = select(&c, &["haeiii".into(), "ecori".into()]).unwrap();
    let reverse = select(&c, &["ecori".into(), "haeiii".into()]).unwrap();
    assert_ne!(
        batch_revision(&c, &pair, role, &m, &o).unwrap(),
        batch_revision(&c, &reverse, role, &m, &o).unwrap()
    );
    let forged = Response {
        id: 1001,
        revision: "0".repeat(64),
        result: R::CatalogRestrictionInspection(r),
    };
    assert!(
        protocol::decode_response(&raw(
            b"CGR\0",
            protocol::WORKER_SCHEMA_VERSION,
            &postcard::to_allocvec(&forged).unwrap()
        ))
        .is_err()
    );
    let forged = Response {
        id: 1001,
        revision: "0".repeat(64),
        result: R::CatalogDigestProduct(p),
    };
    assert!(
        protocol::decode_response(&raw(
            b"CGR\0",
            protocol::WORKER_SCHEMA_VERSION,
            &postcard::to_allocvec(&forged).unwrap()
        ))
        .is_err()
    );
}
#[test]
fn installed_catalog_fragment_has_exact_event_parent_and_native_recipe_pin() {
    let (c, s, m, o) = fixture();
    let w = Workspace::from_documents("catalog native", vec![m], 0).unwrap();
    let d = &w.documents[0];
    let b = SourceBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    };
    let rev = installation_revision(&c, &s, &w, &b, &o, 1).unwrap();
    let R::CatalogDigestInstalled(p) = run(
        rev,
        J::InstallCatalogDigestFragment {
            selection: s.clone(),
            workspace: w.clone(),
            source: b.clone(),
            options: o.clone(),
            fragment_index: 1,
        },
    )
    .result
    else {
        panic!()
    };
    let decoded = collagene_core::workspace::load_workspace(&p.installed.archive).unwrap();
    p.validate_archive_bindings_for(&c, &s, &w, &b, &o, 1, &decoded)
        .unwrap();
    p.validate_archive_for(&c, &s, &w, &b, &o, 1, &decoded)
        .unwrap();
    assert_eq!(
        decoded
            .head(&p.installed.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "CCTT"
    );
    assert_eq!(
        decoded
            .history
            .iter()
            .find(|e| e.id == p.installed.event_id)
            .unwrap()
            .parents,
        vec![b.event_id.clone()]
    );
    assert_eq!(decoded.documents.len(), 2);
    assert_eq!(decoded.head(&b.document_id).unwrap().revision(), b.revision);
    assert_eq!(
        unverified_selection(decoded.head(&p.installed.document_id).unwrap()).unwrap(),
        Some(s)
    );
    let forged = Response {
        id: 1001,
        revision: "0".repeat(64),
        result: R::CatalogDigestInstalled(p),
    };
    assert!(
        protocol::decode_response(&raw(
            b"CGR\0",
            protocol::WORKER_SCHEMA_VERSION,
            &postcard::to_allocvec(&forged).unwrap()
        ))
        .is_err()
    );
}
#[test]
fn catalog_tags_require_twenty_one_and_only_installation_extends_workspace_cap() {
    for (magic, tags) in [
        (b"CGQ\0", vec![70u32, 71, 72]),
        (b"CGR\0", vec![63u32, 64, 65]),
    ] {
        for tag in tags {
            let p = postcard::to_allocvec(&(1002u32, "binding", tag)).unwrap();
            for version in 1..21 {
                let b = raw(magic, version, &p);
                let e = if magic == b"CGQ\0" {
                    protocol::process_request(&b).unwrap_err()
                } else {
                    protocol::decode_response(&b).unwrap_err()
                };
                assert!(matches!(e, collagene_core::Error::Unsupported(_)), "{e}");
            }
            let mut large = p;
            large.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
            let b = raw(magic, protocol::WORKER_SCHEMA_VERSION, &large);
            let e = if magic == b"CGQ\0" {
                protocol::process_request(&b).unwrap_err()
            } else {
                protocol::decode_response(&b).unwrap_err()
            };
            if matches!(tag, 70 | 71 | 63 | 64) {
                assert!(matches!(e, Error::TooLarge));
            } else {
                assert!(
                    !matches!(e, Error::TooLarge),
                    "exact workspace tag remains whitelisted: {e}"
                );
            }
        }
    }
}
