use super::*;
use crate::{AccessKey, MemberCredential};
use collagene_lab_library::collaboration::{ActorId, Author, Role};
use std::{path::PathBuf, sync::atomic::AtomicU64};

struct Fixture {
    directory: PathBuf,
    hosted: HostedMcp,
    alice: String,
    bob: String,
}
impl Fixture {
    fn new(limits: HostedMcpLimits) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let directory = std::env::temp_dir().join(format!(
            "collagene-hosted-owner-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&directory).unwrap();
        let (alice, a) = Self::credential(&directory, "alice");
        let (bob, b) = Self::credential(&directory, "bob");
        let members = MemberRegistry::new(vec![a, b]).unwrap();
        let root = directory.join("managed");
        std::fs::create_dir(&root).unwrap();
        let hosted = HostedMcp::new(
            StorageRoot::existing(root).unwrap(),
            limits,
            members,
            "lab_fixture".into(),
            CanonicalOrigin::parse("https://lab.example").unwrap(),
        )
        .unwrap();
        Self {
            directory,
            hosted,
            alice,
            bob,
        }
    }
    fn credential(directory: &std::path::Path, actor: &str) -> (String, MemberCredential) {
        let key = directory.join(format!("{actor}.key"));
        AccessKey::generate_file(&key).unwrap();
        let authorization = format!("Bearer {}", std::fs::read_to_string(&key).unwrap().trim());
        let member = MemberCredential::new(
            AccessKey::read(&key).unwrap(),
            Author::new(ActorId::parse(actor).unwrap(), actor).unwrap(),
            Role::Contributor,
        )
        .unwrap();
        (authorization, member)
    }
    fn key(&self, actor: &str, label: &str) -> WorkspaceKey {
        WorkspaceKey::new(&self.hosted.0.lab, actor, label).unwrap()
    }
    fn lease(&self, actor: &str, label: &str) -> std::result::Result<Lease, StatusCode> {
        let global = self
            .hosted
            .0
            .permits
            .clone()
            .try_acquire_owned()
            .map_err(|_| StatusCode::TOO_MANY_REQUESTS)?;
        self.hosted.0.lease(self.key(actor, label), true, global)
    }
    fn request(&self, actor: &str, label: &str, method: Method, body: Body) -> Request<Body> {
        Request::builder()
            .method(method)
            .uri(format!("/api/mcp/{label}"))
            .header(header::HOST, "lab.example")
            .header(
                header::AUTHORIZATION,
                if actor == "alice" {
                    &self.alice
                } else {
                    &self.bob
                },
            )
            .header(crate::EXPECTED_LAB_HEADER, &self.hosted.0.lab)
            .header(header::CONTENT_TYPE, "application/json")
            .header(header::ACCEPT, "application/json, text/event-stream")
            .header("mcp-protocol-version", "2026-07-28")
            .body(body)
            .unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        self.hosted.shutdown();
        let _ = std::fs::remove_dir_all(&self.directory);
    }
}

#[test]
fn every_initialize_uses_session_quota_admission_regardless_of_named_protocol() {
    for version in ["2024-11-05", "2025-06-18", "2026-07-28", "2099-01-01"] {
        let message = serde_json::json!({
            "jsonrpc": "2.0", "id": 1, "method": "initialize",
            "params": {"protocolVersion": version, "capabilities": {},
                "clientInfo": {"name": "quota-regression", "version": "1"}}
        });
        assert!(Workspace::starts_protocol_session(&message), "{version}");
    }
    assert!(Workspace::starts_protocol_session(
        &serde_json::json!({"method":"initialize"})
    ));
    for method in [
        "tools/list",
        "tools/call",
        "notifications/initialized",
        "ping",
    ] {
        assert!(!Workspace::starts_protocol_session(&serde_json::json!({
            "method": method, "params": {"protocolVersion":"2026-07-28"}
        })));
    }
}

#[test]
fn scope_keys_reject_paths_and_bind_member_and_startup_lab() {
    for invalid in [
        "",
        "../escape",
        ".",
        "a/b",
        "a\\b",
        "a%2Fb",
        "x?token",
        "space x",
        "é",
        "new\nline",
    ] {
        assert!(WorkspaceKey::new("lab", "alice", invalid).is_err());
    }
    assert!(WorkspaceKey::new("lab", "alice", &"a".repeat(65)).is_err());
    let a = WorkspaceKey::new("lab-a", "alice", "design_1").unwrap();
    assert!(a != WorkspaceKey::new("lab-a", "bob", "design_1").unwrap());
    assert!(a != WorkspaceKey::new("lab-b", "alice", "design_1").unwrap());
    assert!(a != WorkspaceKey::new("lab-a", "alice", "design_2").unwrap());
    let mut limits = HostedMcpLimits::default();
    limits.max_body_bytes = 16 * 1024 * 1024 + 1;
    assert!(limits.validate().is_err());
}

#[test]
fn member_workspace_global_and_live_request_quotas_share_one_pool() {
    let mut limits = HostedMcpLimits::default();
    limits.max_workspaces = 3;
    limits.max_workspaces_per_member = 2;
    let fixture = Fixture::new(limits);
    let a1 = fixture.lease("alice", "one").unwrap();
    let a2 = fixture.lease("alice", "one").unwrap();
    assert_eq!(
        fixture.lease("alice", "one").err(),
        Some(StatusCode::TOO_MANY_REQUESTS)
    );
    let a3 = fixture.lease("alice", "two").unwrap();
    assert_eq!(
        fixture.lease("alice", "three").err(),
        Some(StatusCode::TOO_MANY_REQUESTS)
    );
    let b1 = fixture.lease("bob", "one").unwrap();
    assert!(!Arc::ptr_eq(&a1.workspace, &b1.workspace));
    assert!(!Arc::ptr_eq(&a1.workspace, &a3.workspace));
    assert_eq!(
        fixture.lease("bob", "two").err(),
        Some(StatusCode::TOO_MANY_REQUESTS)
    );
    fixture
        .hosted
        .0
        .sweep(Instant::now() + Duration::from_secs(24 * 60 * 60));
    assert_eq!(
        fixture.hosted.0.workspaces.lock().unwrap().len(),
        3,
        "active scopes cannot idle-expire"
    );
    drop((a1, a2, a3, b1));
    fixture
        .hosted
        .0
        .sweep(Instant::now() + Duration::from_secs(24 * 60 * 60));
    assert_eq!(fixture.hosted.0.workspaces.lock().unwrap().len(), 0);
}

#[test]
fn retired_scope_holds_quota_until_leases_finish_and_replacement_is_fresh() {
    let mut limits = HostedMcpLimits::default();
    limits.max_workspaces = 1;
    limits.max_workspaces_per_member = 1;
    let fixture = Fixture::new(limits);
    let lease = fixture.lease("alice", "design").unwrap();
    let old = lease.workspace.clone();
    fixture
        .hosted
        .0
        .remove(&fixture.key("alice", "design"), &old)
        .unwrap();
    assert!(old.service.config.cancellation_token.is_cancelled());
    assert_eq!(
        fixture.lease("alice", "design").err(),
        Some(StatusCode::TOO_MANY_REQUESTS)
    );
    fixture.hosted.0.sweep(Instant::now());
    assert_eq!(fixture.hosted.0.workspaces.lock().unwrap().len(), 1);
    drop(lease);
    fixture.hosted.0.sweep(Instant::now());
    let fresh = fixture.lease("alice", "design").unwrap();
    assert!(!Arc::ptr_eq(&old, &fresh.workspace));
    assert!(
        !fresh
            .workspace
            .service
            .config
            .cancellation_token
            .is_cancelled()
    );
    assert!(!fresh.workspace.engine.has_pending_compute());
}

#[test]
fn delayed_delete_cannot_retire_a_replacement_under_the_same_workspace_key() {
    let fixture = Fixture::new(HostedMcpLimits::default());
    let key = fixture.key("alice", "design");
    let first_delete = fixture.lease("alice", "design").unwrap();
    let delayed_delete = fixture.lease("alice", "design").unwrap();
    assert!(Arc::ptr_eq(
        &first_delete.workspace,
        &delayed_delete.workspace
    ));

    fixture
        .hosted
        .0
        .remove(&key, &first_delete.workspace)
        .unwrap();
    let replacement = fixture.lease("alice", "design").unwrap();
    assert!(!Arc::ptr_eq(
        &delayed_delete.workspace,
        &replacement.workspace
    ));
    assert_eq!(
        fixture.hosted.0.remove(&key, &delayed_delete.workspace),
        Err(StatusCode::CONFLICT),
        "a request admitted before deletion cannot act on its replacement"
    );
    let pool = fixture.hosted.0.workspaces.lock().unwrap();
    assert!(Arc::ptr_eq(
        pool.current.get(&key).unwrap(),
        &replacement.workspace
    ));
    assert!(
        !replacement
            .workspace
            .service
            .config
            .cancellation_token
            .is_cancelled()
    );
    drop(pool);
    drop((first_delete, delayed_delete));
    fixture.hosted.0.sweep(Instant::now());
    assert_eq!(fixture.hosted.0.workspaces.lock().unwrap().len(), 1);
}

#[test]
fn revocation_cancels_active_scope_and_shutdown_forbids_new_admission() {
    let fixture = Fixture::new(HostedMcpLimits::default());
    let lease = fixture.lease("alice", "design").unwrap();
    let (_, bob) = Fixture::credential(&fixture.directory, "replacement");
    fixture.hosted.0.members.replace(vec![bob]).unwrap();
    fixture.hosted.0.sweep(Instant::now());
    assert!(
        lease
            .workspace
            .service
            .config
            .cancellation_token
            .is_cancelled()
    );
    assert_eq!(fixture.hosted.0.workspaces.lock().unwrap().retired.len(), 1);
    drop(lease);
    fixture.hosted.0.sweep(Instant::now());
    assert_eq!(fixture.hosted.0.workspaces.lock().unwrap().len(), 0);
    fixture.hosted.shutdown();
    assert_eq!(
        fixture.lease("replacement", "new").err(),
        Some(StatusCode::SERVICE_UNAVAILABLE)
    );
}

#[tokio::test]
#[ignore = "application HTTP acceptance is deferred until all in-scope implementation is complete"]
async fn http_auth_origin_lab_and_member_isolation_use_actual_scientific_tools() {
    let fixture = Fixture::new(HostedMcpLimits::default());
    let rpc = |name: &str, arguments: serde_json::Value| {
        serde_json::json!({
            "jsonrpc":"2.0", "id":1, "method":"tools/call", "params":{
                "name":name, "arguments":arguments,
                "_meta":{
                    "io.modelcontextprotocol/protocolVersion":"2026-07-28",
                    "io.modelcontextprotocol/clientInfo":{"name":"collagene-fixture","version":"1"},
                    "io.modelcontextprotocol/clientCapabilities":{}
                }
            }
        })
        .to_string()
    };
    let mut absent = fixture.request("alice", "one", Method::POST, Body::from("not-json"));
    absent.headers_mut().remove(header::AUTHORIZATION);
    assert_eq!(
        fixture.hosted.serve(absent).await.status(),
        StatusCode::UNAUTHORIZED
    );
    assert_eq!(fixture.hosted.0.workspaces.lock().unwrap().len(), 0);
    let mut wrong_origin = fixture.request("alice", "one", Method::POST, Body::empty());
    wrong_origin.headers_mut().insert(
        header::ORIGIN,
        HeaderValue::from_static("https://foreign.example"),
    );
    assert_eq!(
        fixture.hosted.serve(wrong_origin).await.status(),
        StatusCode::FORBIDDEN
    );
    let mut wrong_lab = fixture.request("alice", "one", Method::POST, Body::empty());
    wrong_lab.headers_mut().insert(
        crate::EXPECTED_LAB_HEADER,
        HeaderValue::from_static("other-lab"),
    );
    assert_eq!(
        fixture.hosted.serve(wrong_lab).await.status(),
        StatusCode::PRECONDITION_FAILED
    );
    let imported = fixture
        .hosted
        .serve(fixture.request(
            "alice",
            "one",
            Method::POST,
            Body::from(rpc(
                "import_document",
                serde_json::json!({"text":">retained\nACGTTGCA"}),
            )),
        ))
        .await;
    assert_eq!(imported.status(), StatusCode::OK);
    let import_bytes = to_bytes(imported.into_body(), 1024 * 1024).await.unwrap();
    let report: serde_json::Value = serde_json::from_slice(&import_bytes).unwrap();
    assert!(
        report.get("error").is_none(),
        "actual native import must succeed"
    );
    for (member, workspace, expected) in
        [("alice", "one", 1), ("alice", "two", 0), ("bob", "one", 0)]
    {
        let reply = fixture
            .hosted
            .serve(fixture.request(
                member,
                workspace,
                Method::POST,
                Body::from(rpc("list_documents", serde_json::json!({}))),
            ))
            .await;
        assert_eq!(reply.status(), StatusCode::OK);
        let bytes = to_bytes(reply.into_body(), 1024 * 1024).await.unwrap();
        let value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        let result = &value["result"];
        let data = result.get("structuredContent").cloned().unwrap_or_else(|| {
            serde_json::from_str(result["content"][0]["text"].as_str().unwrap()).unwrap()
        });
        assert_eq!(data["total"], expected);
    }
    let deleted = fixture
        .hosted
        .serve(fixture.request("alice", "one", Method::DELETE, Body::empty()))
        .await;
    assert_eq!(deleted.status(), StatusCode::NO_CONTENT);
}

#[test]
fn operator_storage_cannot_overlap_assets_or_change_after_state_is_shared() {
    fn state(fixture: &Fixture, assets: &std::path::Path) -> ServerState {
        let mut database = collagene_lab_library::Library::open_in_memory().unwrap();
        let lab = database.create_lab("Fixture lab").unwrap();
        let access = database.trusted_host_access(&lab.id).unwrap();
        ServerState::new(
            database,
            access,
            fixture.hosted.0.members.clone(),
            "https://lab.example",
            crate::AssetSnapshot::load(assets).unwrap(),
            crate::Limits::default(),
        )
        .unwrap()
    }
    let fixture = Fixture::new(HostedMcpLimits::default());
    let assets = fixture.directory.join("assets");
    std::fs::create_dir(&assets).unwrap();
    std::fs::write(
        assets.join("index.html"),
        "<html><head></head><body>fixture</body></html>",
    )
    .unwrap();
    let child = assets.join("nested");
    std::fs::create_dir(&child).unwrap();
    for root in [&assets, &fixture.directory, &child] {
        assert!(
            state(&fixture, &assets)
                .with_hosted_mcp(
                    StorageRoot::existing(root).unwrap(),
                    HostedMcpLimits::default()
                )
                .is_err()
        );
    }
    let original = state(&fixture, &assets);
    let held = original.clone();
    assert!(
        original
            .with_hosted_mcp(fixture.hosted.0.root.clone(), HostedMcpLimits::default())
            .is_err()
    );
    assert!(held.0.hosted_mcp.is_none());
    assert!(
        state(&fixture, &assets)
            .with_hosted_mcp(fixture.hosted.0.root.clone(), HostedMcpLimits::default())
            .is_ok()
    );
}

#[tokio::test]
#[ignore = "application HTTP acceptance is deferred until all in-scope implementation is complete"]
async fn legacy_protocol_session_ids_do_not_cross_member_or_workspace_scopes() {
    let fixture = Fixture::new(HostedMcpLimits::default());
    let initialize = serde_json::json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{
        "protocolVersion":"2025-11-25","capabilities":{},"clientInfo":{"name":"fixture","version":"1"}
    }});
    let mut request = fixture.request(
        "alice",
        "one",
        Method::POST,
        Body::from(initialize.to_string()),
    );
    request.headers_mut().insert(
        "mcp-protocol-version",
        HeaderValue::from_static("2025-11-25"),
    );
    let response = fixture.hosted.serve(request).await;
    assert_eq!(response.status(), StatusCode::OK);
    let session = response.headers().get(SESSION_HEADER).unwrap().clone();
    let _ = to_bytes(response.into_body(), 1024 * 1024).await.unwrap();
    // Establish the other scopes directly, so refusal is the SDK's actual
    // independent session manager rather than just an absent workspace label.
    drop(fixture.lease("bob", "one").unwrap());
    drop(fixture.lease("alice", "two").unwrap());
    for (member, workspace) in [("bob", "one"), ("alice", "two")] {
        let mut request = fixture.request(member, workspace, Method::GET, Body::empty());
        request
            .headers_mut()
            .insert(SESSION_HEADER, session.clone());
        request.headers_mut().insert(
            "mcp-protocol-version",
            HeaderValue::from_static("2025-11-25"),
        );
        assert_eq!(
            fixture.hosted.serve(request).await.status(),
            StatusCode::NOT_FOUND
        );
    }
}

#[tokio::test]
#[ignore = "application HTTP acceptance is deferred until all in-scope implementation is complete"]
async fn authentication_precedes_body_poll_and_body_byte_and_time_limits_are_enforced() {
    use axum::body::Bytes;
    let mut limits = HostedMcpLimits::default();
    limits.max_body_bytes = 8;
    limits.body_timeout = Duration::from_millis(10);
    let fixture = Fixture::new(limits);
    let polled = Arc::new(AtomicBool::new(false));
    let observe = polled.clone();
    let body = Body::from_stream(futures_util::stream::once(async move {
        observe.store(true, Ordering::Release);
        Ok::<_, std::convert::Infallible>(Bytes::from_static(b"not-json"))
    }));
    let mut absent = fixture.request("alice", "one", Method::POST, body);
    absent.headers_mut().remove(header::AUTHORIZATION);
    assert_eq!(
        fixture.hosted.serve(absent).await.status(),
        StatusCode::UNAUTHORIZED
    );
    assert!(!polled.load(Ordering::Acquire));
    let too_large = fixture.request("alice", "one", Method::POST, Body::from(vec![b' '; 9]));
    assert_eq!(
        fixture.hosted.serve(too_large).await.status(),
        StatusCode::PAYLOAD_TOO_LARGE
    );
    let pending = Body::from_stream(futures_util::stream::pending::<
        std::result::Result<Bytes, std::convert::Infallible>,
    >());
    assert_eq!(
        fixture
            .hosted
            .serve(fixture.request("alice", "one", Method::POST, pending))
            .await
            .status(),
        StatusCode::REQUEST_TIMEOUT
    );
}
