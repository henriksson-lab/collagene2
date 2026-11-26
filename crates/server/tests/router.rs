//! In-process Router fixtures only: no listener, deployment or application HTTP probe.
use axum::{
    Router,
    body::{Body, Bytes, to_bytes},
    http::{Request, StatusCode, header},
};
use collagene_lab_library::Library;
use collagene_library_protocol as wire;
use collagene_server::{AccessKey, AssetSnapshot, Limits, ServerState};
use std::{
    convert::Infallible,
    path::PathBuf,
    pin::Pin,
    sync::{
        Arc,
        atomic::{AtomicU64, AtomicUsize, Ordering},
    },
    task::{Context, Poll},
    time::Duration,
};
use tower::ServiceExt;
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "collagene-server-fixture-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
struct TestSession {
    authorization: String,
    lab: wire::LabId,
}
fn fixture(limits: Limits) -> (Router, TestSession) {
    source_fixture(limits, collagene_source_gateway::NativeBackend)
}
fn source_fixture<B: collagene_source_gateway::Backend>(
    limits: Limits,
    backend: B,
) -> (Router, TestSession) {
    let temp = Temp::new();
    let assets = temp.0.join("assets");
    std::fs::create_dir(&assets).unwrap();
    std::fs::write(
        assets.join("index.html"),
        b"<html><head><title>Snapshot</title></head><body>snapshot</body></html>",
    )
    .unwrap();
    std::fs::write(
        assets.join("worker-deadbeef.wasm"),
        b"original hashed worker",
    )
    .unwrap();
    let keyfile = temp.0.join("access.key");
    AccessKey::generate_file(&keyfile).unwrap();
    let authorization = format!(
        "Bearer {}",
        std::fs::read_to_string(&keyfile).unwrap().trim()
    );
    let key = AccessKey::read(&keyfile).unwrap();
    let snapshot = AssetSnapshot::load(&assets).unwrap();
    let mut db = Library::open_in_memory().unwrap();
    let lab = db.create_lab("Trusted startup lab").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let state = ServerState::with_source_backend(
        db,
        access,
        key,
        "https://lab.example",
        snapshot,
        limits,
        backend,
    )
    .unwrap();
    // Startup files may disappear: all HTTP static reads are immutable snapshots.
    drop(temp);
    (
        collagene_server::router(state),
        TestSession {
            authorization,
            lab: wire::LabId::new(lab.id.as_str()).unwrap(),
        },
    )
}
fn envelope(request: wire::Request) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: 17,
        request,
    }
}
fn browse() -> wire::RequestEnvelope {
    envelope(wire::Request::Browse {
        parent: None,
        page: wire::PageRequest {
            offset: 0,
            limit: 20,
            expected_catalog_generation: None,
        },
    })
}
fn request(auth: &TestSession, body: Body) -> Request<Body> {
    Request::builder()
        .method("POST")
        .uri("/api/library")
        .header(header::AUTHORIZATION, auth.authorization.as_str())
        .header(collagene_server::EXPECTED_LAB_HEADER, auth.lab.as_str())
        .header(header::HOST, "lab.example")
        .header(header::ORIGIN, "https://lab.example")
        .header(header::CONTENT_TYPE, wire::MEDIA_TYPE)
        .body(body)
        .unwrap()
}
async fn call(
    app: &Router,
    auth: &TestSession,
    request_frame: &wire::RequestEnvelope,
) -> wire::ResponseEnvelope {
    let response = app
        .clone()
        .oneshot(request(
            auth,
            Body::from(wire::encode_request(request_frame).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(response.headers()[header::CONTENT_TYPE], wire::MEDIA_TYPE);
    assert_eq!(response.headers()[header::CACHE_CONTROL], "no-store");
    assert!(
        !response
            .headers()
            .contains_key(header::ACCESS_CONTROL_ALLOW_ORIGIN)
    );
    let bytes = to_bytes(
        response.into_body(),
        wire::MAX_MESSAGE_BYTES + wire::HEADER_BYTES,
    )
    .await
    .unwrap();
    wire::decode_response_for(&bytes, request_frame).unwrap()
}
struct PendingBody(Arc<AtomicUsize>);
impl http_body::Body for PendingBody {
    type Data = Bytes;
    type Error = Infallible;
    fn poll_frame(
        self: Pin<&mut Self>,
        _: &mut Context<'_>,
    ) -> Poll<Option<std::result::Result<http_body::Frame<Bytes>, Infallible>>> {
        self.0.fetch_add(1, Ordering::SeqCst);
        Poll::Pending
    }
}
#[tokio::test]
async fn authentication_origin_and_size_checks_precede_body_polling() {
    let (app, auth) = fixture(Limits::default());
    let polled = Arc::new(AtomicUsize::new(0));
    for (token, origin, host, content_type, length, expected) in [
        (
            "Bearer wrong",
            "https://lab.example",
            "lab.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::UNAUTHORIZED,
        ),
        (
            auth.authorization.as_str(),
            "https://evil.example",
            "lab.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::FORBIDDEN,
        ),
        (
            auth.authorization.as_str(),
            "null",
            "lab.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::FORBIDDEN,
        ),
        (
            auth.authorization.as_str(),
            "https://lab.example",
            "evil.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::FORBIDDEN,
        ),
        (
            auth.authorization.as_str(),
            "https://lab.example",
            "lab.example",
            "text/plain",
            None,
            StatusCode::UNSUPPORTED_MEDIA_TYPE,
        ),
        (
            auth.authorization.as_str(),
            "https://lab.example",
            "lab.example",
            wire::MEDIA_TYPE,
            Some("999999999"),
            StatusCode::PAYLOAD_TOO_LARGE,
        ),
    ] {
        let mut req = request(&auth, Body::new(PendingBody(polled.clone())));
        req.headers_mut()
            .insert(header::AUTHORIZATION, token.parse().unwrap());
        req.headers_mut()
            .insert(header::ORIGIN, origin.parse().unwrap());
        req.headers_mut()
            .insert(header::HOST, host.parse().unwrap());
        req.headers_mut()
            .insert(header::CONTENT_TYPE, content_type.parse().unwrap());
        if let Some(length) = length {
            req.headers_mut()
                .insert(header::CONTENT_LENGTH, length.parse().unwrap());
        }
        let response = app.clone().oneshot(req).await.unwrap();
        assert_eq!(response.status(), expected);
    }
    assert_eq!(polled.load(Ordering::SeqCst), 0);
    let mut req = request(&auth, Body::from(wire::encode_request(&browse()).unwrap()));
    req.headers_mut().remove(header::ORIGIN);
    assert_eq!(
        app.clone().oneshot(req).await.unwrap().status(),
        StatusCode::OK
    );
}
#[tokio::test]
async fn binary_dispatch_cas_and_native_payloads_are_bound_to_startup_lab() {
    let (app, auth) = fixture(Limits::default());
    let empty = call(&app, &auth, &browse()).await;
    assert_eq!(
        empty.context.as_ref().unwrap().name.as_str(),
        "Trusted startup lab"
    );
    let molecule = collagene_core::Molecule::dna(
        "pExample",
        "AGTCCGATGCTA",
        collagene_core::Topology::Circular,
    )
    .unwrap();
    let bytes = collagene_core::protocol::save_document(&molecule).unwrap();
    let payload = wire::Payload {
        kind: wire::PayloadKind::NativeDocument,
        archive_hash: wire::Digest::new(blake3::hash(&bytes).to_hex().as_str()).unwrap(),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: None,
        sharing: None,
    };
    let imported = call(
        &app,
        &auth,
        &envelope(wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("pExample").unwrap(),
            payload: payload.clone(),
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Imported(document)) = imported.outcome else {
        panic!("import failed")
    };
    let saved = call(
        &app,
        &auth,
        &envelope(wire::Request::SaveRevision {
            id: document.id.clone(),
            base_revision: document.head.id.clone(),
            payload: payload.clone(),
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Saved {
        document: updated, ..
    }) = saved.outcome
    else {
        panic!("save failed")
    };
    assert_eq!(updated.head.parent.as_ref(), Some(&document.head.id));
    let stale = call(
        &app,
        &auth,
        &envelope(wire::Request::SaveRevision {
            id: document.id.clone(),
            base_revision: document.head.id.clone(),
            payload,
        }),
    )
    .await;
    let wire::Outcome::Failure(error) = stale.outcome else {
        panic!("stale write accepted")
    };
    assert_eq!(error.code, wire::ErrorCode::Conflict);
    let opened = call(
        &app,
        &auth,
        &envelope(wire::Request::OpenDocument {
            id: document.id,
            revision: Some(document.head.id),
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Opened { payload, .. }) = opened.outcome else {
        panic!("original unavailable")
    };
    assert_eq!(
        collagene_core::protocol::load_document(payload.bytes.as_slice())
            .unwrap()
            .revision(),
        molecule.revision()
    );
}
#[tokio::test]
async fn workspace_saves_require_independently_approved_document_closure() {
    let (app, auth) = fixture(Limits::default());
    let molecule = collagene_core::Molecule::dna(
        "Selected",
        "AGTCCGATGCTA",
        collagene_core::Topology::Circular,
    )
    .unwrap();
    let frame = collagene_core::workspace::Workspace::from_documents(
        "One selected closure",
        vec![molecule.clone()],
        0,
    )
    .unwrap();
    let selected = &frame.documents[0];
    let selection = wire::NativeSelection {
        document_id: wire::Text::new(&selected.id).unwrap(),
        revision: wire::Digest::new(&selected.head).unwrap(),
        head_event: selected
            .head_event
            .as_ref()
            .map(|event| wire::Text::new(event).unwrap()),
    };
    let bytes = collagene_core::workspace::save_workspace(&frame).unwrap();
    let payload = wire::Payload {
        kind: wire::PayloadKind::WorkspaceArchive,
        archive_hash: wire::Digest::new(blake3::hash(&bytes).to_hex().as_str()).unwrap(),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: Some(selection.clone()),
        sharing: Some(wire::WorkspaceSharing::DocumentClosure),
    };
    let imported = call(
        &app,
        &auth,
        &envelope(wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("Selected").unwrap(),
            payload: payload.clone(),
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Imported(document)) = imported.outcome else {
        panic!("minimal native closure rejected: {:?}", imported.outcome)
    };
    let opened = call(
        &app,
        &auth,
        &envelope(wire::Request::OpenDocument {
            id: document.id,
            revision: None,
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Opened {
        payload: retained, ..
    }) = opened.outcome
    else {
        panic!("archive unavailable")
    };
    assert_eq!(retained.bytes, payload.bytes);
    assert_eq!(retained.selection, payload.selection);
    let mut whole = payload;
    whole.sharing = Some(wire::WorkspaceSharing::ExplicitWholeWorkspace);
    let forbidden = call(
        &app,
        &auth,
        &envelope(wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("Unapproved whole workspace").unwrap(),
            payload: whole,
        }),
    )
    .await;
    assert!(matches!(
        forbidden.outcome,
        wire::Outcome::Failure(wire::Failure {
            code: wire::ErrorCode::Forbidden,
            ..
        })
    ));
    let unrelated = collagene_core::Molecule::dna(
        "Unrelated private plasmid",
        "TTCGAATTCAGTC",
        collagene_core::Topology::Linear,
    )
    .unwrap();
    let frame = collagene_core::workspace::Workspace::from_documents(
        "Unrelated browser workspace",
        vec![molecule, unrelated],
        0,
    )
    .unwrap();
    let document = &frame.documents[0];
    let bytes = collagene_core::workspace::save_workspace(&frame).unwrap();
    let extra = wire::Payload {
        kind: wire::PayloadKind::WorkspaceArchive,
        archive_hash: wire::Digest::new(blake3::hash(&bytes).to_hex().as_str()).unwrap(),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: Some(wire::NativeSelection {
            document_id: wire::Text::new(&document.id).unwrap(),
            revision: wire::Digest::new(&document.head).unwrap(),
            head_event: document
                .head_event
                .as_ref()
                .map(|event| wire::Text::new(event).unwrap()),
        }),
        sharing: Some(wire::WorkspaceSharing::DocumentClosure),
    };
    let rejected = call(
        &app,
        &auth,
        &envelope(wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("Do not publish unrelated").unwrap(),
            payload: extra,
        }),
    )
    .await;
    assert!(matches!(
        rejected.outcome,
        wire::Outcome::Failure(wire::Failure {
            code: wire::ErrorCode::Forbidden,
            ..
        })
    ));
}
#[tokio::test]
async fn request_timeouts_and_active_body_concurrency_are_bounded() {
    let limits = Limits {
        max_concurrent_library_requests: 1,
        body_timeout: Duration::from_millis(100),
        ..Limits::default()
    };
    let (app, auth) = fixture(limits);
    let polled = Arc::new(AtomicUsize::new(0));
    let req = request(&auth, Body::new(PendingBody(polled.clone())));
    let pending = tokio::spawn(app.clone().oneshot(req));
    while polled.load(Ordering::SeqCst) == 0 {
        tokio::task::yield_now().await;
    }
    let response = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::TOO_MANY_REQUESTS);
    assert_eq!(
        pending.await.unwrap().unwrap().status(),
        StatusCode::REQUEST_TIMEOUT
    );
    assert!(matches!(
        call(&app, &auth, &browse()).await.outcome,
        wire::Outcome::Success(_)
    ));
    let (app, auth) = fixture(Limits {
        max_body_bytes: 8,
        ..Limits::default()
    });
    let response = app
        .oneshot(request(&auth, Body::from(vec![0u8; 9])))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::PAYLOAD_TOO_LARGE);
}
#[tokio::test]
async fn static_snapshot_preserves_hash_urls_without_filesystem_or_directory_access() {
    let (app, _) = fixture(Limits::default());
    for path in ["/", "/worker-deadbeef.wasm"] {
        let response = app
            .clone()
            .oneshot(Request::builder().uri(path).body(Body::empty()).unwrap())
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::OK);
        let etag = response.headers()[header::ETAG].clone();
        let cached = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri(path)
                    .header(header::IF_NONE_MATCH, etag)
                    .body(Body::empty())
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(cached.status(), StatusCode::NOT_MODIFIED);
    }
    for path in [
        "/%2e%2e/access.key",
        "/.hidden",
        "/missing/",
        "/api/library/",
        "/%5csecret",
    ] {
        assert_eq!(
            app.clone()
                .oneshot(Request::builder().uri(path).body(Body::empty()).unwrap())
                .await
                .unwrap()
                .status(),
            StatusCode::NOT_FOUND
        );
    }
    assert_eq!(
        app.oneshot(
            Request::builder()
                .method("POST")
                .uri("/index.html")
                .body(Body::empty())
                .unwrap()
        )
        .await
        .unwrap()
        .status(),
        StatusCode::METHOD_NOT_ALLOWED
    );
}
#[test]
fn key_creation_is_exclusive_private_and_not_an_implicit_startup_action() {
    let temp = Temp::new();
    let path = temp.0.join("key");
    assert!(AccessKey::read(&path).is_err());
    AccessKey::generate_file(&path).unwrap();
    let first = std::fs::read(&path).unwrap();
    assert_eq!(first.len(), 65);
    assert!(AccessKey::generate_file(&path).is_err());
    assert_eq!(std::fs::read(&path).unwrap(), first);
    AccessKey::read(&path).unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::{PermissionsExt, symlink};
        assert_eq!(
            std::fs::metadata(&path).unwrap().permissions().mode() & 0o777,
            0o600
        );
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o644)).unwrap();
        assert!(AccessKey::read(&path).is_err());
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o600)).unwrap();
        let link = temp.0.join("link");
        symlink(&path, &link).unwrap();
        assert!(AccessKey::read(&link).is_err());
    }
}
#[cfg(unix)]
#[test]
fn static_symlinks_are_rejected_before_snapshotting() {
    use std::os::unix::fs::symlink;
    let temp = Temp::new();
    std::fs::write(temp.0.join("index.html"), "static").unwrap();
    symlink("index.html", temp.0.join("alias.html")).unwrap();
    assert!(AssetSnapshot::load(&temp.0).is_err());
}

#[tokio::test]
async fn server_capability_marker_is_in_memory_exact_and_does_not_mutate_integrity_assets() {
    const MARKER: &str = "<meta name=\"collagene-library\" content=\"enabled\">";
    let temp = Temp::new();
    let root = temp.0.join("assets");
    std::fs::create_dir(&root).unwrap();
    let original = r#"<!doctype html><html><head><!-- <meta name="collagene-library" content="disabled"> --><script>const inert = '<meta name="collagene-library" content="disabled">';</script><script src="app-aabb.js" integrity="sha384-do-not-change" type="module"></script><meta name="collagene-worker" content="/worker-deadbeef.wasm"></head><body></body></html>"#;
    std::fs::write(root.join("index.html"), original).unwrap();
    let js = b"export const unchanged = true;";
    let wasm = b"unchanged hashed worker bytes";
    std::fs::write(root.join("app-aabb.js"), js).unwrap();
    std::fs::write(root.join("worker-deadbeef.wasm"), wasm).unwrap();
    let snapshot = AssetSnapshot::load(&root).unwrap();
    assert_eq!(
        std::fs::read(root.join("index.html")).unwrap(),
        original.as_bytes()
    );
    let key_path = temp.0.join("key");
    AccessKey::generate_file(&key_path).unwrap();
    let key = AccessKey::read(&key_path).unwrap();
    let mut db = Library::open_in_memory().unwrap();
    let lab = db.create_lab("Lab").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let app = collagene_server::router(
        ServerState::new(
            db,
            access,
            key,
            "https://lab.example",
            snapshot,
            Limits::default(),
        )
        .unwrap(),
    );
    let source_marker = "<meta name=\"collagene-source-gateway\" content=\"enabled\">";
    let expected = original.replacen(
        "<head>",
        &format!("<head>\n{source_marker}\n\n{MARKER}\n"),
        1,
    );
    let response = app
        .clone()
        .oneshot(Request::builder().uri("/").body(Body::empty()).unwrap())
        .await
        .unwrap();
    assert_eq!(
        response.headers()[header::ETAG],
        format!("\"{}\"", blake3::hash(expected.as_bytes()).to_hex())
    );
    assert_eq!(
        to_bytes(response.into_body(), 1 << 20)
            .await
            .unwrap()
            .as_ref(),
        expected.as_bytes()
    );
    for (path, bytes) in [
        ("/app-aabb.js", js.as_slice()),
        ("/worker-deadbeef.wasm", wasm.as_slice()),
    ] {
        let response = app
            .clone()
            .oneshot(Request::builder().uri(path).body(Body::empty()).unwrap())
            .await
            .unwrap();
        assert_eq!(
            to_bytes(response.into_body(), 1 << 20)
                .await
                .unwrap()
                .as_ref(),
            bytes
        );
    }
    let valid = format!("<html><head>{MARKER}{source_marker}</head><body></body></html>");
    std::fs::write(root.join("index.html"), &valid).unwrap();
    AssetSnapshot::load(&root).unwrap();
    AssetSnapshot::load(&root).unwrap();
    assert_eq!(
        std::fs::read_to_string(root.join("index.html")).unwrap(),
        valid
    );
    for invalid in [
        format!("<html><head>{MARKER}{MARKER}</head></html>"),
        "<html><head><meta name='collagene-library' name='other' content='enabled'></head></html>".into(),
        "<html><head><meta name='collagene-library' content='enabled' content='disabled'></head></html>".into(),
        "<html><head><meta name='collagene-library' content='disabled'></head></html>".into(),
        "<html><head><meta name='collagene-library'></head></html>".into(),
        format!("<html><head></head><body>{MARKER}</body></html>"),
        "<html><body>No explicit head</body></html>".into(),
        "<html><head><meta name='collagene-source-gateway' content='disabled'></head></html>".into(),
        format!("<html><head>{source_marker}{source_marker}</head></html>"),
        "<html><head><meta name='collagene-source-gateway' name='other' content='enabled'></head></html>".into(),
        "<html><head><meta name='collagene-source-gateway' content='enabled' content='disabled'></head></html>".into(),
        format!("<html><head></head><body>{source_marker}</body></html>"),
    ] {
        std::fs::write(root.join("index.html"),&invalid).unwrap();assert!(AssetSnapshot::load(&root).is_err(),"accepted {invalid}");assert_eq!(std::fs::read_to_string(root.join("index.html")).unwrap(),invalid);
    }
    // Rejected replacement files cannot affect the already-installed snapshot.
    let response = app
        .oneshot(
            Request::builder()
                .uri("/index.html")
                .body(Body::empty())
                .unwrap(),
        )
        .await
        .unwrap();
    assert_eq!(
        to_bytes(response.into_body(), 1 << 20)
            .await
            .unwrap()
            .as_ref(),
        expected.as_bytes()
    );
}

#[tokio::test]
async fn outgoing_frame_clones_and_slices_keep_the_native_response_admission_lease() {
    let (app, auth) = fixture(Limits {
        max_concurrent_library_requests: 1,
        ..Limits::default()
    });
    let molecule = collagene_core::Molecule::dna(
        "Large retained reply",
        &"ACGT".repeat(256 * 1024),
        collagene_core::Topology::Circular,
    )
    .unwrap();
    let bytes = collagene_core::protocol::save_document(&molecule).unwrap();
    let payload = wire::Payload {
        kind: wire::PayloadKind::NativeDocument,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: None,
        sharing: None,
    };
    let imported = call(
        &app,
        &auth,
        &envelope(wire::Request::ImportDocument {
            folder: None,
            name: wire::Text::new("Large retained reply").unwrap(),
            payload,
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Imported(document)) = imported.outcome else {
        panic!("fixture import failed")
    };
    let open = envelope(wire::Request::OpenDocument {
        id: document.id,
        revision: None,
    });
    let response = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&open).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    let mut body = response.into_body();
    let frame = std::future::poll_fn(|cx| http_body::Body::poll_frame(Pin::new(&mut body), cx))
        .await
        .unwrap()
        .unwrap();
    let data = frame.into_data().expect("binary data frame");
    assert!(data.len() > 1024 * 1024);
    drop(body); // Transport has already obtained the frame; Body no longer owns it.
    let retained_clone = data.clone();
    let retained_slice = retained_clone.slice(0..1);
    drop(data);
    let blocked = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(blocked.status(), StatusCode::TOO_MANY_REQUESTS);
    drop(retained_clone); // Even a one-byte slice retains the large backing owner.
    let blocked = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(blocked.status(), StatusCode::TOO_MANY_REQUESTS);
    drop(retained_slice);
    assert!(matches!(
        call(&app, &auth, &browse()).await.outcome,
        wire::Outcome::Success(_)
    ));
}

fn sync_request(auth: &TestSession, body: Body) -> Request<Body> {
    let mut request = request(auth, body);
    *request.uri_mut() = "/api/library/sync".parse().unwrap();
    request.headers_mut().insert(
        header::CONTENT_TYPE,
        collagene_library_sync::MEDIA_TYPE.parse().unwrap(),
    );
    request
}
fn sync_mutation(
    n: u128,
    operation: collagene_library_sync::MutationOperation,
) -> collagene_library_sync::Mutation {
    collagene_library_sync::Mutation {
        id: collagene_library_sync::MutationId::new(format!("mut-{n:032x}")).unwrap(),
        operation,
    }
}
fn sync_envelope(
    auth: &TestSession,
    request: collagene_library_sync::Request,
) -> collagene_library_sync::RequestEnvelope {
    collagene_library_sync::RequestEnvelope {
        request_id: 23,
        expected_lab: auth.lab.clone(),
        request,
    }
}
async fn sync_call(
    app: &Router,
    auth: &TestSession,
    frame: &collagene_library_sync::RequestEnvelope,
) -> collagene_library_sync::ResponseEnvelope {
    use collagene_library_sync as sync;
    let response = app
        .clone()
        .oneshot(sync_request(
            auth,
            Body::from(sync::encode_request(frame).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(response.headers()[header::CONTENT_TYPE], sync::MEDIA_TYPE);
    assert_eq!(response.headers()[header::CACHE_CONTROL], "no-store");
    let bytes = to_bytes(
        response.into_body(),
        sync::MAX_MESSAGE_BYTES + sync::HEADER_BYTES,
    )
    .await
    .unwrap();
    sync::decode_response_for(&bytes, frame).unwrap()
}
fn native_sync_payload(sequence: &str) -> wire::Payload {
    let molecule =
        collagene_core::Molecule::dna("pReliable", sequence, collagene_core::Topology::Circular)
            .unwrap();
    let bytes = collagene_core::protocol::save_document(&molecule).unwrap();
    wire::Payload {
        kind: wire::PayloadKind::NativeDocument,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: None,
        sharing: None,
    }
}
#[tokio::test]
async fn sync_route_authenticates_before_polling_and_separates_media_and_binary_namespaces() {
    use collagene_library_sync as sync;
    let (app, auth) = fixture(Limits {
        max_concurrent_library_requests: 1,
        ..Limits::default()
    });
    let polled = Arc::new(AtomicUsize::new(0));
    let apply = sync_envelope(
        &auth,
        sync::Request::Apply {
            mutation: sync_mutation(
                1,
                sync::MutationOperation::CreateFolder {
                    parent: None,
                    name: wire::Text::new("Sync folder").unwrap(),
                },
            ),
        },
    );
    let mut unauthorized = sync_request(&auth, Body::new(PendingBody(polled.clone())));
    unauthorized
        .headers_mut()
        .insert(header::AUTHORIZATION, "Bearer wrong".parse().unwrap());
    assert_eq!(
        app.clone().oneshot(unauthorized).await.unwrap().status(),
        StatusCode::UNAUTHORIZED
    );
    let mut foreign = sync_request(&auth, Body::new(PendingBody(polled.clone())));
    foreign
        .headers_mut()
        .insert(header::ORIGIN, "https://other.example".parse().unwrap());
    assert_eq!(
        app.clone().oneshot(foreign).await.unwrap().status(),
        StatusCode::FORBIDDEN
    );
    let mut wrong_media = sync_request(&auth, Body::new(PendingBody(polled.clone())));
    wrong_media
        .headers_mut()
        .insert(header::CONTENT_TYPE, wire::MEDIA_TYPE.parse().unwrap());
    assert_eq!(
        app.clone().oneshot(wrong_media).await.unwrap().status(),
        StatusCode::UNSUPPORTED_MEDIA_TYPE
    );
    assert_eq!(polled.load(Ordering::SeqCst), 0);
    let v1_as_sync = sync_request(&auth, Body::from(wire::encode_request(&browse()).unwrap()));
    assert_eq!(
        app.clone().oneshot(v1_as_sync).await.unwrap().status(),
        StatusCode::BAD_REQUEST
    );
    let sync_as_v1 = request(&auth, Body::from(sync::encode_request(&apply).unwrap()));
    assert_eq!(
        app.clone().oneshot(sync_as_v1).await.unwrap().status(),
        StatusCode::BAD_REQUEST
    );
    // Both routes share admission, including an outstanding encoded response.
    let response = app
        .clone()
        .oneshot(sync_request(
            &auth,
            Body::from(sync::encode_request(&apply).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    let blocked = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(blocked.status(), StatusCode::TOO_MANY_REQUESTS);
    drop(response);
    let listed = call(&app, &auth, &browse()).await;
    assert!(
        matches!(listed.outcome,wire::Outcome::Success(wire::Success::Browsed(page))if page.entries.len()==1)
    );
}
#[tokio::test]
async fn lost_sync_import_reply_is_lookup_replayed_without_a_second_document() {
    use collagene_library_sync as sync;
    let (app, auth) = fixture(Limits::default());
    let mutation = sync_mutation(
        10,
        sync::MutationOperation::ImportDocument {
            folder: None,
            name: wire::Text::new("One durable plasmid").unwrap(),
            payload: native_sync_payload("AGTCCGATGCTA"),
        },
    );
    let apply = sync_envelope(
        &auth,
        sync::Request::Apply {
            mutation: mutation.clone(),
        },
    );
    let lost = app
        .clone()
        .oneshot(sync_request(
            &auth,
            Body::from(sync::encode_request(&apply).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(lost.status(), StatusCode::OK);
    drop(lost); // Commit succeeded; client never read its receipt.
    let lookup = sync_envelope(
        &auth,
        sync::Request::Lookup {
            mutation_id: mutation.id.clone(),
            expected_intent: mutation.intent_digest().unwrap(),
        },
    );
    let found = sync_call(&app, &auth, &lookup).await;
    let sync::ResponseOutcome::Receipt(first) = found.outcome else {
        panic!("lost receipt unavailable")
    };
    let repeated = sync_call(&app, &auth, &apply).await;
    assert_eq!(
        repeated.outcome,
        sync::ResponseOutcome::Receipt(first.clone())
    );
    let sync::MutationOutcome::Applied(sync::Applied::Imported(document)) = first.outcome else {
        panic!("import failed")
    };
    let listing = call(&app, &auth, &browse()).await;
    let wire::Outcome::Success(wire::Success::Browsed(page)) = listing.outcome else {
        panic!("legacy browse failed")
    };
    assert_eq!(page.entries.len(), 1);
    let opened = call(
        &app,
        &auth,
        &envelope(wire::Request::OpenDocument {
            id: document.id,
            revision: None,
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Opened {
        revision, payload, ..
    }) = opened.outcome
    else {
        panic!("legacy open failed")
    };
    assert_eq!(revision.number, 1);
    assert_eq!(
        collagene_core::protocol::load_document(&payload.bytes)
            .unwrap()
            .sequence
            .as_str(),
        "AGTCCGATGCTA"
    );
    let changed = sync::Mutation {
        id: mutation.id,
        operation: sync::MutationOperation::CreateFolder {
            parent: None,
            name: wire::Text::new("Different intent").unwrap(),
        },
    };
    let reply = sync_call(
        &app,
        &auth,
        &sync_envelope(&auth, sync::Request::Apply { mutation: changed }),
    )
    .await;
    assert!(matches!(
        reply.outcome,
        sync::ResponseOutcome::IntentMismatch { .. }
    ));
}
#[tokio::test]
async fn sync_cas_rejections_are_durable_original_snapshots_and_legacy_heads_stay_intact() {
    use collagene_library_sync as sync;
    let (app, auth) = fixture(Limits::default());
    let imported = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Apply {
                mutation: sync_mutation(
                    20,
                    sync::MutationOperation::ImportDocument {
                        folder: None,
                        name: wire::Text::new("pCAS").unwrap(),
                        payload: native_sync_payload("AGTCCGATGCTA"),
                    },
                ),
            },
        ),
    )
    .await;
    let sync::ResponseOutcome::Receipt(sync::Receipt {
        outcome: sync::MutationOutcome::Applied(sync::Applied::Imported(original)),
        ..
    }) = imported.outcome
    else {
        panic!("import failed")
    };
    let saved = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Apply {
                mutation: sync_mutation(
                    21,
                    sync::MutationOperation::SaveRevision {
                        id: original.id.clone(),
                        base_revision: original.head.id.clone(),
                        payload: native_sync_payload("AGTCCGATGCTATT"),
                    },
                ),
            },
        ),
    )
    .await;
    let sync::ResponseOutcome::Receipt(sync::Receipt {
        outcome:
            sync::MutationOutcome::Applied(sync::Applied::Saved {
                document: second, ..
            }),
        ..
    }) = saved.outcome
    else {
        panic!("save failed")
    };
    let stale = sync_envelope(
        &auth,
        sync::Request::Apply {
            mutation: sync_mutation(
                22,
                sync::MutationOperation::SaveRevision {
                    id: original.id.clone(),
                    base_revision: original.head.id.clone(),
                    payload: native_sync_payload("AGTCCGATGCTAGG"),
                },
            ),
        },
    );
    let rejected = sync_call(&app, &auth, &stale).await;
    assert!(
        matches!(&rejected.outcome,sync::ResponseOutcome::Receipt(sync::Receipt{outcome:sync::MutationOutcome::Rejected(sync::Rejection{code:sync::RejectionCode::Conflict,conflict:Some(sync::MutationConflict::Document{current}),..}),..})if current.head.id==second.head.id)
    );
    let latest = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Apply {
                mutation: sync_mutation(
                    23,
                    sync::MutationOperation::SaveRevision {
                        id: original.id.clone(),
                        base_revision: second.head.id,
                        payload: native_sync_payload("AGTCCGATGCTACC"),
                    },
                ),
            },
        ),
    )
    .await;
    assert!(matches!(
        latest.outcome,
        sync::ResponseOutcome::Receipt(sync::Receipt {
            outcome: sync::MutationOutcome::Applied(sync::Applied::Saved { .. }),
            ..
        })
    ));
    assert_eq!(
        sync_call(&app, &auth, &stale).await.outcome,
        rejected.outcome
    );
    let original_bytes = call(
        &app,
        &auth,
        &envelope(wire::Request::OpenDocument {
            id: original.id.clone(),
            revision: Some(original.head.id),
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Opened { payload, .. }) = original_bytes.outcome
    else {
        panic!("original revision lost")
    };
    assert_eq!(
        collagene_core::protocol::load_document(&payload.bytes)
            .unwrap()
            .sequence
            .as_str(),
        "AGTCCGATGCTA"
    );
    let head = call(
        &app,
        &auth,
        &envelope(wire::Request::OpenDocument {
            id: original.id,
            revision: None,
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Opened {
        revision, payload, ..
    }) = head.outcome
    else {
        panic!("current unavailable")
    };
    assert_eq!(revision.number, 3);
    assert_eq!(
        collagene_core::protocol::load_document(&payload.bytes)
            .unwrap()
            .sequence
            .as_str(),
        "AGTCCGATGCTACC"
    );
}
#[tokio::test]
async fn sync_payload_validation_rejects_bad_native_bytes_without_creating_a_receipt() {
    use collagene_library_sync as sync;
    let (app, auth) = fixture(Limits::default());
    let bytes = b"not a native document".to_vec();
    let payload = wire::Payload {
        kind: wire::PayloadKind::NativeDocument,
        archive_hash: wire::Digest::of(&bytes),
        bytes: wire::Blob::new(bytes).unwrap(),
        selection: None,
        sharing: None,
    };
    let mutation = sync_mutation(
        30,
        sync::MutationOperation::ImportDocument {
            folder: None,
            name: wire::Text::new("Invalid").unwrap(),
            payload,
        },
    );
    let reply = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Apply {
                mutation: mutation.clone(),
            },
        ),
    )
    .await;
    assert!(matches!(
        reply.outcome,
        sync::ResponseOutcome::Failure(wire::Failure {
            code: wire::ErrorCode::Invalid,
            ..
        })
    ));
    let lookup = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Lookup {
                mutation_id: mutation.id.clone(),
                expected_intent: mutation.intent_digest().unwrap(),
            },
        ),
    )
    .await;
    assert!(matches!(
        lookup.outcome,
        sync::ResponseOutcome::Unknown { .. }
    ));
    let listed = call(&app, &auth, &browse()).await;
    assert!(
        matches!(listed.outcome,wire::Outcome::Success(wire::Success::Browsed(page))if page.entries.is_empty())
    );
}

#[tokio::test]
async fn sync_workspace_writes_retain_exact_history_and_reject_unrelated_document_sharing() {
    use collagene_library_sync as sync;
    fn archived(documents: Vec<collagene_core::Molecule>) -> wire::Payload {
        let frame = collagene_core::workspace::Workspace::from_documents(
            "Selected publication",
            documents,
            0,
        )
        .unwrap();
        let selected = &frame.documents[0];
        let selection = wire::NativeSelection {
            document_id: wire::Text::new(&selected.id).unwrap(),
            revision: wire::Digest::new(&selected.head).unwrap(),
            head_event: selected
                .head_event
                .as_ref()
                .map(|event| wire::Text::new(event).unwrap()),
        };
        let bytes = collagene_core::workspace::save_workspace(&frame).unwrap();
        wire::Payload {
            kind: wire::PayloadKind::WorkspaceArchive,
            archive_hash: wire::Digest::of(&bytes),
            bytes: wire::Blob::new(bytes).unwrap(),
            selection: Some(selection),
            sharing: Some(wire::WorkspaceSharing::DocumentClosure),
        }
    }
    let (app, auth) = fixture(Limits::default());
    let selected = collagene_core::Molecule::dna(
        "Selected",
        "AGTCCGATGCTA",
        collagene_core::Topology::Circular,
    )
    .unwrap();
    let valid = archived(vec![selected.clone()]);
    let imported = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Apply {
                mutation: sync_mutation(
                    40,
                    sync::MutationOperation::ImportDocument {
                        folder: None,
                        name: wire::Text::new("Selected native history").unwrap(),
                        payload: valid.clone(),
                    },
                ),
            },
        ),
    )
    .await;
    let sync::ResponseOutcome::Receipt(sync::Receipt {
        outcome: sync::MutationOutcome::Applied(sync::Applied::Imported(document)),
        ..
    }) = imported.outcome
    else {
        panic!("valid closure rejected")
    };
    let opened = call(
        &app,
        &auth,
        &envelope(wire::Request::OpenDocument {
            id: document.id,
            revision: None,
        }),
    )
    .await;
    let wire::Outcome::Success(wire::Success::Opened { payload, .. }) = opened.outcome else {
        panic!("saved archive unavailable")
    };
    assert_eq!(payload.bytes, valid.bytes);
    assert_eq!(payload.selection, valid.selection);
    let restored = collagene_core::workspace::load_workspace(&payload.bytes).unwrap();
    assert!(!restored.history.is_empty());
    assert_eq!(
        restored.documents[0].head_event.as_deref(),
        payload
            .selection
            .as_ref()
            .unwrap()
            .head_event
            .as_ref()
            .map(|event| event.as_str())
    );
    let private = collagene_core::Molecule::dna(
        "Unrelated private document",
        "TTCGAATTCAGTC",
        collagene_core::Topology::Linear,
    )
    .unwrap();
    let rejected = sync_mutation(
        41,
        sync::MutationOperation::ImportDocument {
            folder: None,
            name: wire::Text::new("Do not leak unrelated work").unwrap(),
            payload: archived(vec![selected, private]),
        },
    );
    let intent = rejected.intent_digest().unwrap();
    let id = rejected.id.clone();
    let reply = sync_call(
        &app,
        &auth,
        &sync_envelope(&auth, sync::Request::Apply { mutation: rejected }),
    )
    .await;
    assert!(matches!(
        reply.outcome,
        sync::ResponseOutcome::Failure(wire::Failure {
            code: wire::ErrorCode::Forbidden,
            ..
        })
    ));
    let lookup = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Lookup {
                mutation_id: id,
                expected_intent: intent,
            },
        ),
    )
    .await;
    assert!(matches!(
        lookup.outcome,
        sync::ResponseOutcome::Unknown { .. }
    ));
}

#[tokio::test]
async fn expected_lab_is_a_prewrite_precondition_even_when_an_access_key_is_reused() {
    use collagene_library_sync as sync;
    let (app, auth) = fixture(Limits::default());
    let header = collagene_server::EXPECTED_LAB_HEADER;
    let wrong = wire::LabId::new(
        if auth.lab.as_str() == "lab-ffffffffffffffffffffffffffffffff" {
            "lab-00000000000000000000000000000000"
        } else {
            "lab-ffffffffffffffffffffffffffffffff"
        },
    )
    .unwrap();
    let polled = Arc::new(AtomicUsize::new(0));
    for sync_route in [false, true] {
        for mode in ["wrong", "malformed", "duplicate"] {
            let mut req = if sync_route {
                sync_request(&auth, Body::new(PendingBody(polled.clone())))
            } else {
                request(&auth, Body::new(PendingBody(polled.clone())))
            };
            match mode {
                "wrong" => {
                    req.headers_mut()
                        .insert(header, wrong.as_str().parse().unwrap());
                }
                "malformed" => {
                    req.headers_mut()
                        .insert(header, "not-a-lab-id".parse().unwrap());
                }
                _ => {
                    req.headers_mut()
                        .append(header, wrong.as_str().parse().unwrap());
                }
            }
            assert_eq!(
                app.clone().oneshot(req).await.unwrap().status(),
                if mode == "wrong" {
                    StatusCode::PRECONDITION_FAILED
                } else {
                    StatusCode::BAD_REQUEST
                }
            );
        }
    }
    let mut absent = sync_request(&auth, Body::new(PendingBody(polled.clone())));
    absent.headers_mut().remove(header);
    assert_eq!(
        app.clone().oneshot(absent).await.unwrap().status(),
        StatusCode::PRECONDITION_REQUIRED
    );
    assert_eq!(polled.load(Ordering::SeqCst), 0);
    // Legacy discovery remains possible before the client knows the lab ID.
    let mut bootstrap = request(&auth, Body::from(wire::encode_request(&browse()).unwrap()));
    bootstrap.headers_mut().remove(header);
    assert_eq!(
        app.clone().oneshot(bootstrap).await.unwrap().status(),
        StatusCode::OK
    );
    let folder = wire::FolderId::new("fld-00000000000000000000000000000001").unwrap();
    let doc = wire::DocumentId::new("doc-00000000000000000000000000000001").unwrap();
    let revision = wire::RevisionId::new("a".repeat(64)).unwrap();
    let name = wire::Text::new("Must not mutate without a pin").unwrap();
    let payload = native_sync_payload("AGTCCGATGCTA");
    for operation in [
        wire::Request::CreateFolder {
            parent: None,
            name: name.clone(),
        },
        wire::Request::RenameFolder {
            id: folder.clone(),
            expected_generation: 1,
            name: name.clone(),
        },
        wire::Request::MoveFolder {
            id: folder.clone(),
            expected_generation: 1,
            parent: None,
        },
        wire::Request::DeleteEmptyFolder {
            id: folder,
            expected_generation: 1,
        },
        wire::Request::ImportDocument {
            folder: None,
            name: name.clone(),
            payload: payload.clone(),
        },
        wire::Request::SaveRevision {
            id: doc.clone(),
            base_revision: revision,
            payload,
        },
        wire::Request::RenameDocument {
            id: doc.clone(),
            expected_generation: 1,
            name,
        },
        wire::Request::MoveDocument {
            id: doc,
            expected_generation: 1,
            folder: None,
        },
    ] {
        let mut req = request(
            &auth,
            Body::from(wire::encode_request(&envelope(operation)).unwrap()),
        );
        req.headers_mut().remove(header);
        assert_eq!(
            app.clone().oneshot(req).await.unwrap().status(),
            StatusCode::PRECONDITION_REQUIRED
        );
    }
    let mutation = sync_mutation(
        90,
        sync::MutationOperation::CreateFolder {
            parent: None,
            name: wire::Text::new("Pinned sync folder").unwrap(),
        },
    );
    let intent = mutation.intent_digest().unwrap();
    let mut contradictory = sync_envelope(
        &auth,
        sync::Request::Apply {
            mutation: mutation.clone(),
        },
    );
    contradictory.expected_lab = wrong;
    assert_eq!(
        app.clone()
            .oneshot(sync_request(
                &auth,
                Body::from(sync::encode_request(&contradictory).unwrap())
            ))
            .await
            .unwrap()
            .status(),
        StatusCode::PRECONDITION_FAILED
    );
    let unknown = sync_call(
        &app,
        &auth,
        &sync_envelope(
            &auth,
            sync::Request::Lookup {
                mutation_id: mutation.id.clone(),
                expected_intent: intent,
            },
        ),
    )
    .await;
    assert!(matches!(
        unknown.outcome,
        sync::ResponseOutcome::Unknown { .. }
    ));
    let listed = call(&app, &auth, &browse()).await;
    assert!(
        matches!(listed.outcome,wire::Outcome::Success(wire::Success::Browsed(page))if page.entries.is_empty())
    );
    assert!(matches!(
        call(
            &app,
            &auth,
            &envelope(wire::Request::CreateFolder {
                parent: None,
                name: wire::Text::new("Pinned v1 folder").unwrap()
            })
        )
        .await
        .outcome,
        wire::Outcome::Success(wire::Success::CreatedFolder(_))
    ));
    assert!(matches!(
        sync_call(
            &app,
            &auth,
            &sync_envelope(&auth, sync::Request::Apply { mutation })
        )
        .await
        .outcome,
        sync::ResponseOutcome::Receipt(sync::Receipt {
            outcome: sync::MutationOutcome::Applied(sync::Applied::CreatedFolder(_)),
            ..
        })
    ));
}

struct FixtureSource {
    calls: Arc<AtomicUsize>,
    pending: bool,
    started: Arc<tokio::sync::Notify>,
    dropped: Arc<tokio::sync::Notify>,
}
struct SourceDrop(Arc<tokio::sync::Notify>);
impl Drop for SourceDrop {
    fn drop(&mut self) {
        self.0.notify_one();
    }
}
impl collagene_source_gateway::Backend for FixtureSource {
    async fn import(
        &self,
        _: String,
        _: Option<String>,
    ) -> Result<collagene_core::formats::ImportReport, collagene_sources::SourceError> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        let _drop = SourceDrop(self.dropped.clone());
        self.started.notify_one();
        if self.pending {
            std::future::pending::<()>().await;
        }
        Ok(collagene_core::formats::ImportReport {
            molecules: vec![
                collagene_core::Molecule::dna(
                    "Fixture only",
                    &"ACGT".repeat(300_000),
                    collagene_core::Topology::Circular,
                )
                .unwrap(),
            ],
            warnings: vec![],
        })
    }
    async fn catalog(
        &self,
        _: String,
        _: Option<String>,
    ) -> Result<collagene_sources::AddgeneCatalog, collagene_sources::SourceError> {
        panic!("unexpected catalog backend call")
    }
    async fn selection(
        &self,
        _: Arc<collagene_sources::AddgeneCatalog>,
        _: u64,
        _: Option<String>,
    ) -> Result<collagene_core::formats::ImportReport, collagene_sources::SourceError> {
        panic!("unexpected selection backend call")
    }
}
fn fake_source(pending: bool) -> FixtureSource {
    FixtureSource {
        calls: Arc::new(AtomicUsize::new(0)),
        pending,
        started: Arc::new(tokio::sync::Notify::new()),
        dropped: Arc::new(tokio::sync::Notify::new()),
    }
}
fn source_request(auth: &TestSession, body: Body) -> Request<Body> {
    let mut req = request(auth, body);
    *req.uri_mut() = collagene_source_gateway::PATH.parse().unwrap();
    req.headers_mut().insert(
        header::CONTENT_TYPE,
        collagene_source_gateway::CONTENT_TYPE.parse().unwrap(),
    );
    req
}
fn source_frame() -> Vec<u8> {
    use collagene_source_gateway::wire as source;
    source::encode_request(&source::Request {
        id: 51,
        client: [17; 32],
        operation: source::Operation::Import {
            source: "https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1".into(),
            addgene_token: None,
        },
    })
    .unwrap()
}
#[tokio::test]
async fn source_route_auth_lab_origin_media_and_bounds_reject_before_body_poll() {
    let backend = fake_source(false);
    let calls = backend.calls.clone();
    let (app, auth) = source_fixture(Limits::default(), backend);
    let polled = Arc::new(AtomicUsize::new(0));
    for case in 0..8 {
        let mut req = source_request(&auth, Body::new(PendingBody(polled.clone())));
        let expected = match case {
            0 => {
                req.headers_mut().remove(header::AUTHORIZATION);
                StatusCode::UNAUTHORIZED
            }
            1 => {
                req.headers_mut()
                    .insert(header::HOST, "evil.example".parse().unwrap());
                StatusCode::FORBIDDEN
            }
            2 => {
                req.headers_mut()
                    .insert(header::ORIGIN, "https://evil.example".parse().unwrap());
                StatusCode::FORBIDDEN
            }
            3 => {
                req.headers_mut()
                    .remove(collagene_server::EXPECTED_LAB_HEADER);
                StatusCode::PRECONDITION_REQUIRED
            }
            4 => {
                req.headers_mut().insert(
                    collagene_server::EXPECTED_LAB_HEADER,
                    "lab-00000000000000000000000000000000".parse().unwrap(),
                );
                StatusCode::PRECONDITION_FAILED
            }
            5 => {
                req.headers_mut().append(
                    collagene_server::EXPECTED_LAB_HEADER,
                    auth.lab.as_str().parse().unwrap(),
                );
                StatusCode::BAD_REQUEST
            }
            6 => {
                req.headers_mut()
                    .insert(header::CONTENT_TYPE, wire::MEDIA_TYPE.parse().unwrap());
                StatusCode::UNSUPPORTED_MEDIA_TYPE
            }
            _ => {
                req.headers_mut()
                    .insert(header::CONTENT_LENGTH, "16385".parse().unwrap());
                StatusCode::PAYLOAD_TOO_LARGE
            }
        };
        assert_eq!(app.clone().oneshot(req).await.unwrap().status(), expected);
    }
    assert_eq!(polled.load(Ordering::SeqCst), 0);
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    let oversized = app
        .clone()
        .oneshot(source_request(&auth, Body::from(vec![0; 16385])))
        .await
        .unwrap();
    assert_eq!(oversized.status(), StatusCode::PAYLOAD_TOO_LARGE);
    let wrong_namespace = app
        .oneshot(source_request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(wrong_namespace.status(), StatusCode::BAD_REQUEST);
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}
#[tokio::test]
async fn source_reply_uses_shared_frame_lease_and_keeps_v1_available_after_release() {
    use collagene_source_gateway::wire as source;
    let backend = fake_source(false);
    let calls = backend.calls.clone();
    let (app, auth) = source_fixture(
        Limits {
            max_concurrent_library_requests: 1,
            ..Limits::default()
        },
        backend,
    );
    let response = app
        .clone()
        .oneshot(source_request(&auth, Body::from(source_frame())))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(
        response.headers()[header::CONTENT_TYPE],
        collagene_source_gateway::CONTENT_TYPE
    );
    let mut body = response.into_body();
    let frame = std::future::poll_fn(|cx| http_body::Body::poll_frame(Pin::new(&mut body), cx))
        .await
        .unwrap()
        .unwrap()
        .into_data()
        .unwrap();
    drop(body);
    assert!(frame.len() > 1024 * 1024);
    let decoded = source::decode_response(&frame).unwrap();
    assert_eq!(decoded.id, 51);
    assert!(matches!(decoded.result, Ok(source::Output::Imported(_))));
    assert_eq!(calls.load(Ordering::SeqCst), 1);
    let slice = frame.slice(0..1);
    drop(frame);
    let blocked = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(blocked.status(), StatusCode::TOO_MANY_REQUESTS);
    drop(slice);
    assert!(matches!(
        call(&app, &auth, &browse()).await.outcome,
        wire::Outcome::Success(_)
    ));
    // An authenticated generic URL cannot invoke even the injected backend.
    let invalid = source::encode_request(&source::Request {
        id: 52,
        client: [17; 32],
        operation: source::Operation::Import {
            source: "https://evil.example/arbitrary.gb".into(),
            addgene_token: None,
        },
    })
    .unwrap();
    let response = app
        .oneshot(source_request(&auth, Body::from(invalid)))
        .await
        .unwrap();
    let bytes = to_bytes(response.into_body(), source::MAX_RESPONSE)
        .await
        .unwrap();
    assert!(source::decode_response(&bytes).unwrap().result.is_err());
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}
#[tokio::test]
async fn dropping_source_caller_cleans_up_provider_future_and_active_nonce() {
    let backend = fake_source(true);
    let started = backend.started.clone();
    let dropped = backend.dropped.clone();
    let calls = backend.calls.clone();
    let (app, auth) = source_fixture(
        Limits {
            max_concurrent_library_requests: 1,
            ..Limits::default()
        },
        backend,
    );
    for _ in 0..2 {
        let req = source_request(&auth, Body::from(source_frame()));
        let app = app.clone();
        let task = tokio::spawn(async move { app.oneshot(req).await });
        tokio::time::timeout(Duration::from_secs(2), started.notified())
            .await
            .unwrap();
        task.abort();
        let _ = task.await;
        tokio::time::timeout(Duration::from_secs(2), dropped.notified())
            .await
            .unwrap();
    }
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert!(matches!(
        call(&app, &auth, &browse()).await.outcome,
        wire::Outcome::Success(_)
    ));
}

#[tokio::test]
async fn source_cancel_bypasses_occupied_heavy_lane_without_releasing_reply_leases() {
    use collagene_source_gateway::wire as source;
    let backend = fake_source(true);
    let started = backend.started.clone();
    let dropped = backend.dropped.clone();
    let (app, auth) = source_fixture(
        Limits {
            max_concurrent_library_requests: 1,
            ..Limits::default()
        },
        backend,
    );
    let req = source_request(&auth, Body::from(source_frame()));
    let provider_app = app.clone();
    let task = tokio::spawn(async move { provider_app.oneshot(req).await });
    tokio::time::timeout(Duration::from_secs(2), started.notified())
        .await
        .unwrap();
    let blocked = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(blocked.status(), StatusCode::TOO_MANY_REQUESTS);
    let cancel = source::encode_request(&source::Request {
        id: 52,
        client: [17; 32],
        operation: source::Operation::Cancel { request_id: 51 },
    })
    .unwrap();
    let cancel_reply = tokio::time::timeout(
        Duration::from_secs(2),
        app.clone()
            .oneshot(source_request(&auth, Body::from(cancel))),
    )
    .await
    .unwrap()
    .unwrap();
    assert_eq!(cancel_reply.status(), StatusCode::OK);
    let cancel_bytes = to_bytes(cancel_reply.into_body(), 4096).await.unwrap();
    assert!(matches!(
        source::decode_response(&cancel_bytes).unwrap().result,
        Ok(source::Output::Canceled { found: true })
    ));
    // Cancellation is delivered through the gateway, without aborting the caller.
    tokio::time::timeout(Duration::from_secs(2), dropped.notified())
        .await
        .unwrap();
    let reply = tokio::time::timeout(Duration::from_secs(2), task)
        .await
        .unwrap()
        .unwrap()
        .unwrap();
    let mut body = reply.into_body();
    let frame = std::future::poll_fn(|cx| http_body::Body::poll_frame(Pin::new(&mut body), cx))
        .await
        .unwrap()
        .unwrap()
        .into_data()
        .unwrap();
    assert!(matches!(
        source::decode_response(&frame).unwrap().result,
        Err(source::Failure {
            code: source::ErrorCode::Canceled,
            ..
        })
    ));
    let slice = frame.slice(0..1);
    drop(frame);
    drop(body);
    let blocked = app
        .clone()
        .oneshot(request(
            &auth,
            Body::from(wire::encode_request(&browse()).unwrap()),
        ))
        .await
        .unwrap();
    assert_eq!(blocked.status(), StatusCode::TOO_MANY_REQUESTS);
    // Dispose also works while the large-response admission is retained. Its
    // bounded reply owns a separate ingress slot, not the heavy lease.
    let dispose = source::encode_request(&source::Request {
        id: 53,
        client: [17; 32],
        operation: source::Operation::DisposeCatalog { handle: [0; 32] },
    })
    .unwrap();
    let response = app
        .clone()
        .oneshot(source_request(&auth, Body::from(dispose)))
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    drop(response);
    drop(slice);
    assert!(matches!(
        call(&app, &auth, &browse()).await.outcome,
        wire::Outcome::Success(_)
    ));
}
