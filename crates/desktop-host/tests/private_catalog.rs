//! Real offline IPC/repository restart for private archives, never an HTTP path.
use collagene_core::{
    Molecule, Topology,
    workspace::{self, Workspace, WorkspaceAttachment},
};
use collagene_desktop_host::{DesktopHost, WindowAccess};
use collagene_desktop_protocol as wire;
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
use wire::library;
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-private-catalog-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )))
    }
}
impl Drop for Disk {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn call(
    host: &mut DesktopHost,
    access: &WindowAccess,
    id: u64,
    request: library::Request,
) -> wire::Outcome {
    let request = wire::RequestEnvelope {
        request_id: id,
        operation_id: wire::OperationId::new(format!("op-{id:032x}")).unwrap(),
        command: wire::Command::Library(request),
    };
    let bytes = wire::encode_request(&request).unwrap();
    let reply = host.dispatch_raw(access, &bytes).unwrap();
    wire::decode_response_for(&reply, &request).unwrap().outcome
}
fn saved(out: wire::Outcome) -> library::Document {
    match out {
        wire::Outcome::Success(wire::Success::Library(library::Success::Imported(d))) => d,
        other => panic!("{other:?}"),
    }
}
fn frame() -> Workspace {
    let raw=br#"{"@context":{"sbol":"http://sbols.org/v3#"},"@graph":[{"@id":"https://private.invalid/unrelated","@type":"https://private.invalid/Activity","https://private.invalid/note":"PRESERVE-EXACT-ORIGINAL-LOCALLY"}]}"#.to_vec();
    let hash = library::Digest::of(&raw).as_str().to_string();
    let mut stock = Molecule::dna(
        "Local plasmid",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    for (k, v) in [
        ("sbol:source_hash", hash.as_str()),
        ("sbol:component_iri", "https://private.invalid/c"),
        ("sbol:sequence_iri", "https://private.invalid/s"),
        (
            "sbol:exchange_profile",
            "collagene-sbol-3.1-structural-jsonld-1",
        ),
    ] {
        stock.metadata.insert(k.into(), v.into());
    }
    let mut w = Workspace::from_documents(
        "Full local snapshot",
        vec![
            stock,
            Molecule::dna("Private second plasmid", "CCCGGGTTTAAA", Topology::Circular).unwrap(),
        ],
        0,
    )
    .unwrap();
    w.attachments.insert(
        hash,
        WorkspaceAttachment {
            kind: "sbol-design-exchange-v1".into(),
            name: "Exact source.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw,
        },
    );
    w
}
fn payload(w: &Workspace) -> library::Payload {
    let bytes = workspace::save_workspace(w).unwrap();
    let d = &w.documents[0];
    library::Payload {
        kind: library::PayloadKind::WorkspaceArchive,
        archive_hash: library::Digest::of(&bytes),
        bytes: library::Blob::new(bytes).unwrap(),
        selection: Some(library::NativeSelection {
            document_id: library::NativeId::new(&d.id).unwrap(),
            revision: library::Digest::new(&d.head).unwrap(),
            head_event: d
                .head_event
                .as_ref()
                .map(|e| library::NativeId::new(e).unwrap()),
        }),
        sharing: Some(library::WorkspaceSharing::ExplicitWholeWorkspace),
    }
}
#[test]
fn private_local_catalog_keeps_original_graph_native_history_and_unrelated_local_branches_after_restart()
 {
    let disk = Disk::new();
    let w = frame();
    let original = payload(&w);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let import = library::Request::ImportDocument {
        folder: None,
        name: library::Name::new("Private offline plasmid").unwrap(),
        payload: original.clone(),
    };
    let doc = saved(call(&mut host, &access, 1, import.clone()));
    drop(host);
    let mut host = DesktopHost::open(&disk.0).unwrap();
    let access = host.authorize_window().unwrap();
    let replay = saved(call(&mut host, &access, 1, import));
    assert_eq!(replay, doc);
    match call(
        &mut host,
        &access,
        2,
        library::Request::OpenDocument {
            id: doc.id.clone(),
            revision: None,
        },
    ) {
        wire::Outcome::Success(wire::Success::Library(library::Success::Opened {
            payload,
            ..
        })) => {
            assert_eq!(payload, original);
            let reopened = workspace::load_workspace(payload.bytes.as_slice()).unwrap();
            assert_eq!(
                serde_json::to_value(&reopened).unwrap(),
                serde_json::to_value(&w).unwrap()
            );
        }
        other => panic!("{other:?}"),
    }
    // The archive and selected native root remain exact in immutable history too.
    match call(
        &mut host,
        &access,
        3,
        library::Request::ListRevisions {
            id: doc.id,
            page: library::PageRequest::default(),
        },
    ) {
        wire::Outcome::Success(wire::Success::Library(library::Success::Revisions(page))) => {
            assert_eq!(page.revisions.len(), 1);
            assert_eq!(page.revisions[0].selection, original.selection);
            assert_eq!(page.revisions[0].archive_hash, original.archive_hash);
        }
        other => panic!("{other:?}"),
    }
}
