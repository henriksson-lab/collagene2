//! Authenticated same-origin HTTP adapter. No client identity/lab fields grant
//! permission; operator-configured individual credentials share one startup lab.
mod assets;
mod auth;
mod hosted_mcp;
mod members;
pub use assets::AssetSnapshot;
pub use auth::AccessKey;
use axum::{
    Router,
    body::{Body, Bytes, to_bytes},
    extract::State,
    http::{HeaderMap, HeaderValue, Method, Request, StatusCode, header},
    middleware::{self, Next},
    response::{IntoResponse, Response},
    routing::{any, post},
};
use collagene_lab_library::{LabAccess, Library};
use collagene_library_host as host;
use collagene_library_protocol as wire;
pub use hosted_mcp::HostedMcpLimits;
pub use members::{AuthenticatedMember, MemberCredential, MemberRegistry};
use wire::discussion;

use collagene_library_sync as sync;
use collagene_source_gateway::{self as sources, wire as source_wire};
use std::{
    future::Future,
    pin::Pin,
    sync::{Arc, Mutex},
    time::Duration,
};
use tokio::sync::{OwnedSemaphorePermit, Semaphore};

pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid server configuration: {0}")]
    Configuration(&'static str),
    #[error("operating-system random generation failed")]
    Entropy,
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Storage(#[from] collagene_lab_library::Error),
    #[error(transparent)]
    Backup(#[from] collagene_lab_library::backup::BackupError),
}
#[derive(Debug, Clone)]
pub struct Limits {
    pub max_concurrent_library_requests: usize,
    pub max_body_bytes: usize,
    pub body_timeout: Duration,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_concurrent_library_requests: 2,
            max_body_bytes: wire::MAX_MESSAGE_BYTES + wire::HEADER_BYTES,
            body_timeout: Duration::from_secs(30),
        }
    }
}
impl Limits {
    fn validate(&self) -> Result<()> {
        if self.max_concurrent_library_requests == 0
            || self.max_concurrent_library_requests > 16
            || self.max_body_bytes == 0
            || self.max_body_bytes > wire::MAX_MESSAGE_BYTES + wire::HEADER_BYTES
            || self.body_timeout.is_zero()
            || self.body_timeout > Duration::from_secs(120)
        {
            return Err(Error::Configuration("request limits"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct CanonicalOrigin(String);
impl CanonicalOrigin {
    fn parse(value: &str) -> Result<Self> {
        if value.trim() != value || value.chars().any(char::is_control) {
            return Err(Error::Configuration("origin"));
        }
        let url = url::Url::parse(value).map_err(|_| Error::Configuration("origin"))?;
        if !matches!(url.scheme(), "http" | "https")
            || url.host_str().is_none()
            || !url.username().is_empty()
            || url.password().is_some()
            || url.query().is_some()
            || url.fragment().is_some()
            || url.path() != "/"
        {
            return Err(Error::Configuration("origin"));
        }
        Ok(Self(url.origin().ascii_serialization()))
    }
    fn accepts(&self, headers: &HeaderMap) -> bool {
        let mut hosts = headers.get_all(header::HOST).iter();
        let Some(host) = hosts.next().and_then(|h| h.to_str().ok()) else {
            return false;
        };
        if hosts.next().is_some()
            || host.contains(['/', '\\', '@', '?', '#'])
            || host.chars().any(char::is_whitespace)
        {
            return false;
        }
        let scheme = self.0.split_once("://").map(|p| p.0).unwrap_or("");
        let Ok(host_origin) = Self::parse(&format!("{scheme}://{host}")) else {
            return false;
        };
        if host_origin.0 != self.0 {
            return false;
        }
        let mut origins = headers.get_all(header::ORIGIN).iter();
        if let Some(origin) = origins.next() {
            if origins.next().is_some() {
                return false;
            }
            let Ok(origin) = origin.to_str() else {
                return false;
            };
            let Ok(origin) = Self::parse(origin) else {
                return false;
            };
            if origin.0 != self.0 {
                return false;
            }
        }
        true
    }
}
type SourceFuture = Pin<Box<dyn Future<Output = source_wire::Response> + Send>>;
type SourceExecute =
    dyn Fn(source_wire::Request, sources::Cancellation) -> SourceFuture + Send + Sync;
struct Inner {
    source: Box<SourceExecute>,
    source_ingress: Arc<Semaphore>,
    database: Mutex<Library>,
    access: LabAccess,
    members: MemberRegistry,
    origin: CanonicalOrigin,
    assets: AssetSnapshot,
    permits: Arc<Semaphore>,
    limits: Limits,
    hosted_mcp: Option<hosted_mcp::HostedMcp>,
}
#[derive(Clone)]
pub struct ServerState(Arc<Inner>);
impl ServerState {
    /// All parameters are trusted operator configuration. Authentication on every
    /// API request grants this fixed lab only; archives never choose another lab.
    pub fn new(
        database: Library,
        access: LabAccess,
        members: MemberRegistry,
        origin: &str,
        assets: AssetSnapshot,
        limits: Limits,
    ) -> Result<Self> {
        Self::with_source_backend(
            database,
            access,
            members,
            origin,
            assets,
            limits,
            sources::NativeBackend,
        )
    }
    /// Injection is a trusted native host decision, never a request parameter.
    pub fn with_source_backend<B: sources::Backend>(
        database: Library,
        access: LabAccess,
        members: MemberRegistry,
        origin: &str,
        assets: AssetSnapshot,
        limits: Limits,
        backend: B,
    ) -> Result<Self> {
        limits.validate()?;
        database.lab(&access)?;
        // Only execute() is exposed: this config never creates another listener
        // or authorizes requests; the owning server enforces its lab policy.
        let gateway = sources::Gateway::with_backend(sources::Config::default(), backend);
        Ok(Self(Arc::new(Inner {
            source_ingress: Arc::new(Semaphore::new(4)),
            source: Box::new(move |request, cancellation| {
                let gateway = gateway.clone();
                Box::pin(async move { gateway.execute(request, cancellation).await })
            }),
            database: Mutex::new(database),
            access,
            members,
            origin: CanonicalOrigin::parse(origin)?,
            assets,
            permits: Arc::new(Semaphore::new(limits.max_concurrent_library_requests)),
            limits,
            hosted_mcp: None,
        })))
    }

    /// Trusted opt-in: hosted scientific workspaces have no shared-library or
    /// arbitrary filesystem authority. Configure before cloning this state.
    pub fn with_hosted_mcp(
        mut self,
        root: collagene_local_storage::StorageRoot,
        limits: HostedMcpLimits,
    ) -> Result<Self> {
        let inner = Arc::get_mut(&mut self.0).ok_or(Error::Configuration(
            "configure MCP before sharing server state",
        ))?;
        if inner.assets.overlaps_root(root.canonical_path()) {
            return Err(Error::Configuration(
                "MCP storage and asset roots must be separate",
            ));
        }
        inner.hosted_mcp = Some(hosted_mcp::HostedMcp::new(
            root,
            limits,
            inner.members.clone(),
            inner.access.lab_id().as_str().to_owned(),
            inner.origin.clone(),
        )?);
        Ok(self)
    }

    /// Close scientific workspace transports before graceful listener shutdown.
    pub fn shutdown_hosted_mcp(&self) {
        if let Some(hosted) = &self.0.hosted_mcp {
            hosted.shutdown();
        }
    }
}
pub fn router(state: ServerState) -> Router {
    Router::new()
        .route("/api/library", post(library))
        .route("/api/library/sync", post(reliable_library))
        .route(sources::PATH, post(source_import))
        .route(discussion::ENDPOINT, post(discussion_request))
        .route("/api/mcp/{workspace}", any(hosted_mcp::request))
        .fallback(static_asset)
        .layer(middleware::from_fn(response_headers))
        .with_state(state)
}
async fn response_headers(request: Request<Body>, next: Next) -> Response {
    let mut response = next.run(request).await;
    response.headers_mut().insert(
        header::X_CONTENT_TYPE_OPTIONS,
        HeaderValue::from_static("nosniff"),
    );
    response
        .headers_mut()
        .insert(header::X_FRAME_OPTIONS, HeaderValue::from_static("DENY"));
    response.headers_mut().insert(
        header::REFERRER_POLICY,
        HeaderValue::from_static("no-referrer"),
    );
    response
}
fn rejection(status: StatusCode) -> Response {
    let mut response = (
        status,
        status.canonical_reason().unwrap_or("Request rejected"),
    )
        .into_response();
    response
        .headers_mut()
        .insert(header::CACHE_CONTROL, HeaderValue::from_static("no-store"));
    if status == StatusCode::UNAUTHORIZED {
        response
            .headers_mut()
            .insert(header::WWW_AUTHENTICATE, HeaderValue::from_static("Bearer"));
    }
    response
}
#[derive(Clone, Copy)]
enum Endpoint {
    Library,
    Reliable,
    Discussion,
    Source,
}
impl Endpoint {
    fn media_type(self) -> &'static str {
        match self {
            Self::Library => wire::MEDIA_TYPE,
            Self::Reliable => sync::MEDIA_TYPE,
            Self::Source => sources::CONTENT_TYPE,
            Self::Discussion => discussion::MEDIA_TYPE,
        }
    }
}
async fn library(State(state): State<ServerState>, request: Request<Body>) -> Response {
    binary_library(state, request, Endpoint::Library).await
}
async fn reliable_library(State(state): State<ServerState>, request: Request<Body>) -> Response {
    binary_library(state, request, Endpoint::Reliable).await
}
async fn discussion_request(State(state): State<ServerState>, request: Request<Body>) -> Response {
    binary_library(state, request, Endpoint::Discussion).await
}
async fn source_import(State(state): State<ServerState>, request: Request<Body>) -> Response {
    binary_library(state, request, Endpoint::Source).await
}
pub const EXPECTED_LAB_HEADER: &str = "x-collagene-expected-lab";
fn expected_lab_pin(
    headers: &HeaderMap,
    expected: &str,
    required: bool,
) -> std::result::Result<bool, StatusCode> {
    let mut values = headers.get_all(EXPECTED_LAB_HEADER).iter();
    let Some(value) = values.next() else {
        return if required {
            Err(StatusCode::PRECONDITION_REQUIRED)
        } else {
            Ok(false)
        };
    };
    if values.next().is_some() {
        return Err(StatusCode::BAD_REQUEST);
    }
    let pin = value
        .to_str()
        .ok()
        .and_then(|value| wire::LabId::new(value).ok())
        .ok_or(StatusCode::BAD_REQUEST)?;
    if pin.as_str() != expected {
        return Err(StatusCode::PRECONDITION_FAILED);
    }
    Ok(true)
}
async fn binary_library(
    state: ServerState,
    request: Request<Body>,
    endpoint: Endpoint,
) -> Response {
    // Authenticate before reading, buffering, decoding or dispatching any body.
    let headers = request.headers();
    let mut authorization = headers.get_all(header::AUTHORIZATION).iter();
    let Some(credential) = authorization.next().cloned() else {
        return rejection(StatusCode::UNAUTHORIZED);
    };
    if authorization.next().is_some()
        || state
            .0
            .members
            .authenticate(credential.as_bytes())
            .is_none()
    {
        return rejection(StatusCode::UNAUTHORIZED);
    }
    if !state.0.origin.accepts(headers) {
        return rejection(StatusCode::FORBIDDEN);
    }
    let lab_pinned = match expected_lab_pin(
        headers,
        state.0.access.lab_id().as_str(),
        matches!(
            endpoint,
            Endpoint::Reliable | Endpoint::Source | Endpoint::Discussion
        ),
    ) {
        Ok(present) => present,
        Err(status) => return rejection(status),
    };
    if request.uri().query().is_some()
        || request.uri().authority().is_some()
        || request.uri().scheme().is_some()
    {
        return rejection(StatusCode::BAD_REQUEST);
    }
    if headers
        .get(header::CONTENT_TYPE)
        .is_none_or(|value| value.as_bytes() != endpoint.media_type().as_bytes())
        || headers.get_all(header::CONTENT_TYPE).iter().count() != 1
    {
        return rejection(StatusCode::UNSUPPORTED_MEDIA_TYPE);
    }
    if headers
        .get(header::CONTENT_ENCODING)
        .is_some_and(|value| value.as_bytes() != b"identity")
    {
        return rejection(StatusCode::UNSUPPORTED_MEDIA_TYPE);
    }
    if headers.get_all(header::CONTENT_LENGTH).iter().count() > 1
        || headers.get_all(header::CONTENT_ENCODING).iter().count() > 1
    {
        return rejection(StatusCode::BAD_REQUEST);
    }
    let max_body = if matches!(endpoint, Endpoint::Source) {
        state.0.limits.max_body_bytes.min(source_wire::MAX_REQUEST)
    } else if matches!(endpoint, Endpoint::Discussion) {
        state
            .0
            .limits
            .max_body_bytes
            .min(discussion::MAX_REQUEST_BYTES)
    } else {
        state.0.limits.max_body_bytes
    };
    if let Some(length) = headers.get(header::CONTENT_LENGTH) {
        let Ok(length) = length.to_str().unwrap_or("").parse::<u64>() else {
            return rejection(StatusCode::BAD_REQUEST);
        };
        if length > max_body as u64 {
            return rejection(StatusCode::PAYLOAD_TOO_LARGE);
        }
    }
    // Tiny source frames are admitted separately so Cancel/Dispose can reach
    // an occupied gateway. Ordinary source work still needs the heavy lease.
    let admission = if matches!(endpoint, Endpoint::Source) {
        &state.0.source_ingress
    } else {
        &state.0.permits
    };
    let Ok(permit) = admission.clone().try_acquire_owned() else {
        return rejection(StatusCode::TOO_MANY_REQUESTS);
    };
    let permit = Arc::new(permit);
    let bytes = match tokio::time::timeout(
        state.0.limits.body_timeout,
        to_bytes(request.into_body(), max_body),
    )
    .await
    {
        Ok(Ok(bytes)) => bytes,
        Ok(Err(_)) => return rejection(StatusCode::PAYLOAD_TOO_LARGE),
        Err(_) => return rejection(StatusCode::REQUEST_TIMEOUT),
    };
    if state
        .0
        .members
        .authenticate(credential.as_bytes())
        .is_none()
    {
        return rejection(StatusCode::UNAUTHORIZED);
    }
    if matches!(endpoint, Endpoint::Source) {
        return execute_source(state, bytes, permit).await;
    }
    // Cancellation cannot release admission while native work is still running.
    // The same lease also follows encoded response bytes through transport frames.
    let work_permit = permit.clone();
    let result = tokio::task::spawn_blocking(move || {
        let _permit = work_permit;
        enum Decoded {
            Library(wire::RequestEnvelope),
            Reliable(sync::RequestEnvelope),
            Discussion(discussion::RequestEnvelope),
        }
        let decoded = match endpoint {
            Endpoint::Source => return Err(StatusCode::BAD_REQUEST),
            Endpoint::Discussion => Decoded::Discussion(
                discussion::decode_request(&bytes).map_err(|_| StatusCode::BAD_REQUEST)?,
            ),
            Endpoint::Library => {
                let request = wire::decode_request(&bytes).map_err(|_| StatusCode::BAD_REQUEST)?;
                if !lab_pinned && request.request.is_mutation() {
                    return Err(StatusCode::PRECONDITION_REQUIRED);
                }
                Decoded::Library(request)
            }
            Endpoint::Reliable => {
                // Validate the intended database in the borrowed frame prefix
                // before allocating a decoded archive or touching SQLite.
                let prefix = sync::request_prefix(&bytes).map_err(|_| StatusCode::BAD_REQUEST)?;
                if prefix.expected_lab != state.0.access.lab_id().as_str() {
                    return Err(StatusCode::PRECONDITION_FAILED);
                }
                let request = sync::decode_request(&bytes).map_err(|_| StatusCode::BAD_REQUEST)?;
                Decoded::Reliable(request)
            }
        };
        let mut db = state
            .0
            .database
            .lock()
            .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
        // Resolve identity after waiting for the database lock: membership may
        // have changed while this request's body was decoded or work was queued.
        let member = state
            .0
            .members
            .authenticate(credential.as_bytes())
            .ok_or(StatusCode::UNAUTHORIZED)?;
        match &decoded {
            Decoded::Library(request) if !member.permits_library(&request.request) => {
                return Err(StatusCode::FORBIDDEN);
            }
            Decoded::Reliable(request) if !member.permits_sync(&request.request) => {
                return Err(StatusCode::FORBIDDEN);
            }
            _ => (),
        }
        match decoded {
            Decoded::Discussion(request) => {
                let access = db
                    .trusted_member_access(&state.0.access, member.author, member.role)
                    .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
                let response = host::discussion::dispatch(&mut db, &access, &request)
                    .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
                discussion::encode_response(&response).map_err(|_| StatusCode::SERVICE_UNAVAILABLE)
            }
            Decoded::Reliable(request) => {
                // The native adapter performs exact receipt lookup before a retry,
                // and fresh native/publication validation before any new mutation.
                // It exposes no checkpoint, SQL, local-scope or administration API.
                let response = sync::native::dispatch(&mut db, &state.0.access, &request).map_err(
                    |error| match error {
                        sync::native::DispatchError::WrongLab => StatusCode::PRECONDITION_FAILED,
                        _ => StatusCode::SERVICE_UNAVAILABLE,
                    },
                )?;
                sync::encode_response(&response).map_err(|_| StatusCode::SERVICE_UNAVAILABLE)
            }
            Decoded::Library(request) => {
                let payload = match &request.request {
                    wire::Request::ImportDocument { payload, .. }
                    | wire::Request::SaveRevision { payload, .. } => Some(payload),
                    _ => None,
                };
                let scope = payload
                    .map(|payload| host::document_publication_scope(&state.0.access, payload))
                    .transpose();
                // A failed publication review is dispatched without a scope so the host
                // returns its bounded typed failure, retaining request correlation.
                let scope = scope.ok().flatten().flatten();
                let response = host::dispatch(&mut db, &state.0.access, &request, scope.as_ref())
                    .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
                wire::encode_response(&response).map_err(|_| StatusCode::SERVICE_UNAVAILABLE)
            }
        }
    })
    .await;
    match result {
        Ok(Ok(bytes)) => leased_response(bytes, permit, endpoint.media_type()),
        Ok(Err(status)) => rejection(status),
        Err(_) => rejection(StatusCode::SERVICE_UNAVAILABLE),
    }
}
// Dropping an HTTP caller cancels provider work and its active nonce entry.
struct CancelSource(sources::Cancellation);
impl Drop for CancelSource {
    fn drop(&mut self) {
        self.0.cancel();
    }
}
async fn execute_source(
    state: ServerState,
    bytes: Bytes,
    permit: Arc<OwnedSemaphorePermit>,
) -> Response {
    let request = match source_wire::decode_request(&bytes) {
        Ok(request) => request,
        Err(_) => return rejection(StatusCode::BAD_REQUEST),
    };
    let control = matches!(
        &request.operation,
        source_wire::Operation::Cancel { .. } | source_wire::Operation::DisposeCatalog { .. }
    );
    let permit = if control {
        // Control replies retain the bounded small ingress slot, never a large
        // response slot. Four slow bodies/replies cannot grow this lane.
        permit
    } else {
        let Ok(heavy) = state.0.permits.clone().try_acquire_owned() else {
            return rejection(StatusCode::TOO_MANY_REQUESTS);
        };
        drop(permit);
        Arc::new(heavy)
    };
    let cancellation = CancelSource(sources::Cancellation::default());
    // The gateway also enforces its shorter provider deadline and retains its
    // own native admission until an already-running bounded parser finishes.
    let response = match tokio::time::timeout(
        Duration::from_secs(90),
        (state.0.source)(request, cancellation.0.clone()),
    )
    .await
    {
        Ok(response) => response,
        Err(_) => return rejection(StatusCode::GATEWAY_TIMEOUT),
    };
    let lease = permit.clone();
    let encoded = tokio::task::spawn_blocking(move || {
        let _lease = lease;
        let bytes = source_wire::encode_response(&response)?;
        if control && bytes.len() > 4096 {
            return Err(source_wire::WireError::TooLarge);
        }
        Ok(bytes)
    })
    .await;
    match encoded {
        Ok(Ok(bytes)) => leased_response(bytes, permit, sources::CONTENT_TYPE),
        _ => rejection(StatusCode::SERVICE_UNAVAILABLE),
    }
}
fn leased_response(
    bytes: Vec<u8>,
    permit: Arc<OwnedSemaphorePermit>,
    media: &'static str,
) -> Response {
    let leased = Bytes::from_owner(LeaseBytes {
        data: bytes,
        _lease: permit,
    });
    let mut response = Response::new(Body::from(leased));
    response
        .headers_mut()
        .insert(header::CONTENT_TYPE, HeaderValue::from_static(media));
    response
        .headers_mut()
        .insert(header::CACHE_CONTROL, HeaderValue::from_static("no-store"));
    response
}
// The lease belongs to the Bytes backing allocation, rather than just Body:
// Hyper can retain a data frame after polling/dropping its response Body. Clones
// and slices keep the owner alive until the final reference is released.
struct LeaseBytes {
    data: Vec<u8>,
    _lease: Arc<OwnedSemaphorePermit>,
}
impl AsRef<[u8]> for LeaseBytes {
    fn as_ref(&self) -> &[u8] {
        &self.data
    }
}
async fn static_asset(State(state): State<ServerState>, request: Request<Body>) -> Response {
    if request.method() != Method::GET && request.method() != Method::HEAD {
        return rejection(StatusCode::METHOD_NOT_ALLOWED);
    }
    let Some(path) = assets::asset_path(request.uri().path()) else {
        return rejection(StatusCode::NOT_FOUND);
    };
    let Some(asset) = state.0.assets.files.get(&path) else {
        return rejection(StatusCode::NOT_FOUND);
    };
    let unchanged = request
        .headers()
        .get(header::IF_NONE_MATCH)
        .is_some_and(|value| value.as_bytes() == asset.etag.as_bytes());
    let mut response = Response::new(if unchanged || request.method() == Method::HEAD {
        Body::empty()
    } else {
        Body::from(asset.bytes.clone())
    });
    if unchanged {
        *response.status_mut() = StatusCode::NOT_MODIFIED;
    } else {
        response.headers_mut().insert(
            header::CONTENT_LENGTH,
            HeaderValue::from(asset.bytes.len() as u64),
        );
    }
    response
        .headers_mut()
        .insert(header::CONTENT_TYPE, HeaderValue::from_static(asset.mime));
    if let Ok(etag) = HeaderValue::from_str(&asset.etag) {
        response.headers_mut().insert(header::ETAG, etag);
    }
    response
        .headers_mut()
        .insert(header::CACHE_CONTROL, HeaderValue::from_static("no-cache"));
    response
}
