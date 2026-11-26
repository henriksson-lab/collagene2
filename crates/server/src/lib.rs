//! Authenticated same-origin HTTP adapter. No client identity/lab fields grant
//! permission; one trusted startup lab is shared by holders of its access key.
mod assets;
mod auth;
pub use assets::AssetSnapshot;
pub use auth::AccessKey;
use axum::{
    Router,
    body::{Body, Bytes, to_bytes},
    extract::State,
    http::{HeaderMap, HeaderValue, Method, Request, StatusCode, header},
    middleware::{self, Next},
    response::{IntoResponse, Response},
    routing::post,
};
use collagene_lab_library::{LabAccess, Library};
use collagene_library_host as host;
use collagene_library_protocol as wire;
use std::{
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
struct Inner {
    database: Mutex<Library>,
    access: LabAccess,
    key: AccessKey,
    origin: CanonicalOrigin,
    assets: AssetSnapshot,
    permits: Arc<Semaphore>,
    limits: Limits,
}
#[derive(Clone)]
pub struct ServerState(Arc<Inner>);
impl ServerState {
    /// All parameters are trusted operator configuration. Authentication on every
    /// API request grants this fixed lab only; archives never choose another lab.
    pub fn new(
        database: Library,
        access: LabAccess,
        key: AccessKey,
        origin: &str,
        assets: AssetSnapshot,
        limits: Limits,
    ) -> Result<Self> {
        limits.validate()?;
        database.lab(&access)?;
        Ok(Self(Arc::new(Inner {
            database: Mutex::new(database),
            access,
            key,
            origin: CanonicalOrigin::parse(origin)?,
            assets,
            permits: Arc::new(Semaphore::new(limits.max_concurrent_library_requests)),
            limits,
        })))
    }
}
pub fn router(state: ServerState) -> Router {
    Router::new()
        .route("/api/library", post(library))
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
async fn library(State(state): State<ServerState>, request: Request<Body>) -> Response {
    // Authenticate before reading, buffering, decoding or dispatching any body.
    let headers = request.headers();
    let mut authorization = headers.get_all(header::AUTHORIZATION).iter();
    let valid = authorization
        .next()
        .is_some_and(|value| state.0.key.authenticates(value.as_bytes()))
        && authorization.next().is_none();
    if !valid {
        return rejection(StatusCode::UNAUTHORIZED);
    }
    if !state.0.origin.accepts(headers) {
        return rejection(StatusCode::FORBIDDEN);
    }
    if request.uri().query().is_some() {
        return rejection(StatusCode::BAD_REQUEST);
    }
    if headers
        .get(header::CONTENT_TYPE)
        .is_none_or(|value| value.as_bytes() != wire::MEDIA_TYPE.as_bytes())
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
    if headers.get_all(header::CONTENT_LENGTH).iter().count() > 1 {
        return rejection(StatusCode::BAD_REQUEST);
    }
    if let Some(length) = headers.get(header::CONTENT_LENGTH) {
        let Ok(length) = length.to_str().unwrap_or("").parse::<u64>() else {
            return rejection(StatusCode::BAD_REQUEST);
        };
        if length > state.0.limits.max_body_bytes as u64 {
            return rejection(StatusCode::PAYLOAD_TOO_LARGE);
        }
    }
    let Ok(permit) = state.0.permits.clone().try_acquire_owned() else {
        return rejection(StatusCode::TOO_MANY_REQUESTS);
    };
    let permit = Arc::new(permit);
    let bytes = match tokio::time::timeout(
        state.0.limits.body_timeout,
        to_bytes(request.into_body(), state.0.limits.max_body_bytes),
    )
    .await
    {
        Ok(Ok(bytes)) => bytes,
        Ok(Err(_)) => return rejection(StatusCode::PAYLOAD_TOO_LARGE),
        Err(_) => return rejection(StatusCode::REQUEST_TIMEOUT),
    };
    // Cancellation cannot release admission while native work is still running.
    // The same lease also follows encoded response bytes through transport frames.
    let work_permit = permit.clone();
    let result = tokio::task::spawn_blocking(move || {
        let _permit = work_permit;
        let request = wire::decode_request(&bytes).map_err(|_| StatusCode::BAD_REQUEST)?;
        let mut db = state
            .0
            .database
            .lock()
            .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
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
    })
    .await;
    match result {
        Ok(Ok(bytes)) => {
            let leased = Bytes::from_owner(LeaseBytes {
                data: bytes,
                _lease: permit,
            });
            let mut response = Response::new(Body::from(leased));
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static(wire::MEDIA_TYPE),
            );
            response
                .headers_mut()
                .insert(header::CACHE_CONTROL, HeaderValue::from_static("no-store"));
            response
        }
        Ok(Err(status)) => rejection(status),
        Err(_) => rejection(StatusCode::SERVICE_UNAVAILABLE),
    }
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
