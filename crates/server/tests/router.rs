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
fn fixture(limits: Limits) -> (Router, String) {
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
    let state = ServerState::new(db, access, key, "https://lab.example", snapshot, limits).unwrap();
    // Startup files may disappear: all HTTP static reads are immutable snapshots.
    drop(temp);
    (collagene_server::router(state), authorization)
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
fn request(auth: &str, body: Body) -> Request<Body> {
    Request::builder()
        .method("POST")
        .uri("/api/library")
        .header(header::AUTHORIZATION, auth)
        .header(header::HOST, "lab.example")
        .header(header::ORIGIN, "https://lab.example")
        .header(header::CONTENT_TYPE, wire::MEDIA_TYPE)
        .body(body)
        .unwrap()
}
async fn call(
    app: &Router,
    auth: &str,
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
            &auth,
            "https://evil.example",
            "lab.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::FORBIDDEN,
        ),
        (
            &auth,
            "null",
            "lab.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::FORBIDDEN,
        ),
        (
            &auth,
            "https://lab.example",
            "evil.example",
            wire::MEDIA_TYPE,
            None,
            StatusCode::FORBIDDEN,
        ),
        (
            &auth,
            "https://lab.example",
            "lab.example",
            "text/plain",
            None,
            StatusCode::UNSUPPORTED_MEDIA_TYPE,
        ),
        (
            &auth,
            "https://lab.example",
            "lab.example",
            wire::MEDIA_TYPE,
            Some("999999999"),
            StatusCode::PAYLOAD_TOO_LARGE,
        ),
    ] {
        let mut req = request(token, Body::new(PendingBody(polled.clone())));
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
    let expected = original.replacen("<head>", &format!("<head>\n{MARKER}\n"), 1);
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
    let valid = format!("<html><head>{MARKER}</head><body></body></html>");
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
