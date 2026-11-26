use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    protocol,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material,
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
};
use collagene_desktop_mcp::{NativeSelection, create};
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_sharing::Selection;
use collagene_repository as repository;
use rmcp::{ServiceExt, model::CallToolRequestParams};
use serde_json::{Value, json};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Root(PathBuf);
impl Root {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-desktop-mcp-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )))
    }
}
impl Drop for Root {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn event(frame: &Workspace, index: usize) -> EventBinding {
    let d = &frame.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        event_id: d.head_event.clone().unwrap(),
    }
}
async fn call(
    client: &rmcp::service::RunningService<rmcp::RoleClient, ()>,
    name: &str,
    args: Value,
) -> Value {
    let reply = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
#[test]
fn selection_parser_requires_exact_native_catalog_pins_and_empty_session_never_opens_root() {
    let pin = format!("doc-{}@{}", "1".repeat(32), "a".repeat(64));
    NativeSelection::parse(&pin).unwrap();
    for invalid in [
        "",
        "document-1",
        "doc-1@head",
        &format!("{pin}@extra"),
        &format!(" {pin}"),
    ] {
        assert!(NativeSelection::parse(invalid).is_err());
    }
    let root = Root::new();
    assert!(create(&root.0, &[]).unwrap().preloaded.is_empty());
    assert!(
        !root.0.exists(),
        "empty MCP must not open even the native repository"
    );
    let selections = vec![NativeSelection::parse(&pin).unwrap(); 33];
    assert!(create(&root.0, &selections).is_err());
    assert!(!root.0.exists());
}
#[tokio::test]
async fn selected_saved_pcr_closure_replays_fromrun_without_private_documents_or_native_writes() {
    let root = Root::new();
    let _desktop = collagene_desktop_host::DesktopHost::open(&root.0).unwrap();
    let mut repo =
        collagene_repository::Repository::open(root.0.join("repository.sqlite")).unwrap();
    let source = Molecule::dna(
        "Original source",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent target",
        source.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let secret = Molecule::dna(
        "Private unpublished plasmid",
        "GGGGTTTTAAAACCCCGGGGTTTT",
        Topology::Circular,
    )
    .unwrap();
    let secret_revision = secret.revision();
    let frame = Workspace::from_documents(
        "Private local working archive",
        vec![source.clone(), target.clone(), secret],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Real native PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let saved = run_lineage::record_run(
        &frame,
        &plan,
        event(&frame, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&frame, 0),
        }],
        false,
    )
    .unwrap();
    let product = &saved.receipt.products[0];
    let selected = Selection {
        document_id: product.document_id.clone(),
        revision: product.revision.clone(),
        head_event: Some(product.event_id.clone()),
    };
    let origin = run_material::Selection {
        run_id: saved.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            revision: product.revision.clone(),
            event_id: product.event_id.clone(),
        },
    };
    let payload = repository::private_workspace(&saved.workspace, &selected).unwrap();
    let original_bytes = payload.bytes.clone();
    let document = match repo
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Reusable selected stock".into(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        _ => panic!(),
    };
    let checkpoint = repo.checkpoint_version();
    let local_context = repo.local_context().unwrap();
    let pin = NativeSelection {
        document: document.id.clone(),
        revision: document.head.id.clone(),
    };
    let working_path = root.0.join("working.sqlite");
    let working_before = std::fs::read(&working_path).unwrap();
    let session = create(&root.0, &[pin]).unwrap();
    assert_eq!(
        std::fs::read(&working_path).unwrap(),
        working_before,
        "read-only MCP preload must preserve the live desktop working store"
    );
    let inventory = session.preloaded[0].inventory.as_ref().unwrap();
    assert!(!inventory.revisions.contains(&secret_revision));
    assert!(inventory.omitted.documents >= 1);
    let workspace = session.preloaded[0].receipt["workspace"].clone();
    let reference = json!({"id":workspace["id"],"revision":workspace["revision"]});
    let (a, b) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        session
            .into_server()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    let client = ().serve(b).await.unwrap();
    let listed = call(&client, "list_documents", json!({})).await;
    assert!(!listed.to_string().contains(&secret_revision));
    assert!(!listed.to_string().contains("Private unpublished plasmid"));
    let replay = call(
        &client,
        "replay_workspace_run",
        json!({"workspace":reference,"run_id":saved.receipt.run_id}),
    )
    .await;
    assert_eq!(replay["current_status"], "ConsistentUnderAssumptions");
    assert_eq!(replay["products_match"], true);
    let material = call(
        &client,
        "prepare_run_material",
        json!({"workspace":reference,"selection":origin}),
    )
    .await;
    assert_eq!(material["capsule_available"], true);
    // Session-created artifacts stay in MCP memory, never the native repository.
    repo.reload().unwrap();
    assert_eq!(repo.checkpoint_version(), checkpoint);
    assert_eq!(repo.local_context().unwrap(), local_context);
    assert_eq!(
        repo.open_document(&document.id, Some(&document.head.id))
            .unwrap()
            .bytes,
        original_bytes
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn bare_native_snapshot_is_literal_and_failed_exact_selection_never_serves_partial_session() {
    let root = Root::new();
    let desktop = collagene_desktop_host::DesktopHost::open(&root.0).unwrap();
    let mut repo =
        collagene_repository::Repository::open(root.0.join("repository.sqlite")).unwrap();
    let molecule =
        Molecule::dna("Literal native file", "AACCGGTTAACCGGTT", Topology::Linear).unwrap();
    let bytes = protocol::save_document(&molecule).unwrap();
    let document = match repo
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Literal".into(),
            payload: store::Payload {
                kind: store::PayloadKind::NativeDocument,
                archive_hash: store::Digest::parse(blake3::hash(&bytes).to_hex().to_string())
                    .unwrap(),
                bytes,
                selection: None,
                sharing: None,
            },
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        _ => panic!(),
    };
    let pin = NativeSelection {
        document: document.id.clone(),
        revision: document.head.id.clone(),
    };
    // A GUI working database can contain unsaved private work. The MCP factory
    // must not even open it (the sentinel is intentionally not SQLite).
    drop(desktop); // Close this fixture's real working store before replacing it.
    let working_path = root.0.join("working.sqlite");
    let mut options = std::fs::OpenOptions::new();
    options.write(true).truncate(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.mode(0o600);
    }
    use std::io::Write;
    let mut working = options.open(&working_path).unwrap();
    working
        .write_all(b"private unsaved GUI recovery sentinel")
        .unwrap();
    drop(working);
    let session = create(&root.0, &[pin.clone()]).unwrap();
    assert_eq!(
        std::fs::read(&working_path).unwrap(),
        b"private unsaved GUI recovery sentinel"
    );
    assert!(session.preloaded[0].inventory.is_none());
    let bad = NativeSelection {
        document: document.id.clone(),
        revision: store::RevisionId::parse("0".repeat(64)).unwrap(),
    };
    assert!(create(&root.0, &[pin, bad]).is_err());
    let (a, b) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        session
            .into_server()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    let client = ().serve(b).await.unwrap();
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 1);
    assert_eq!(
        call(&client, "list_workspaces", json!({})).await["workspaces"],
        json!([])
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}
