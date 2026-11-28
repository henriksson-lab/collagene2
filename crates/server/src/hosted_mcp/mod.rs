//! Member-scoped scientific workspaces over the SDK's Streamable HTTP transport.
//! Workspace labels select memory state, never filesystem paths or shared-library rights.
use crate::{
    CanonicalOrigin, Error, MemberRegistry, Result, ServerState, expected_lab_pin, rejection,
};
use axum::{
    body::{Body, to_bytes},
    extract::State,
    http::{HeaderMap, HeaderValue, Method, Request, StatusCode, header},
    response::{IntoResponse, Response},
};
use collagene_local_storage::StorageRoot;
use collagene_mcp::server::CollageneServer;
use futures_util::StreamExt;
use rmcp::transport::streamable_http_server::{
    StreamableHttpServerConfig, StreamableHttpService, session::local::LocalSessionManager,
};
use std::{
    collections::BTreeMap,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, AtomicUsize, Ordering},
    },
    time::{Duration, Instant},
};
use tokio::sync::{OwnedSemaphorePermit, Semaphore};

const SESSION_HEADER: &str = "mcp-session-id";
const MAX_WORKSPACE_LABEL: usize = 64;
const SWEEP_INTERVAL: Duration = Duration::from_secs(30);

#[derive(Clone, Debug)]
pub struct HostedMcpLimits {
    pub max_workspaces: usize,
    pub max_workspaces_per_member: usize,
    pub max_concurrent_requests: usize,
    pub max_concurrent_workspace_requests: usize,
    pub max_legacy_sessions_per_workspace: usize,
    pub max_body_bytes: usize,
    pub body_timeout: Duration,
    pub idle_timeout: Duration,
    pub response_lifetime: Duration,
}
impl Default for HostedMcpLimits {
    fn default() -> Self {
        Self {
            max_workspaces: 16,
            max_workspaces_per_member: 4,
            max_concurrent_requests: 16,
            max_concurrent_workspace_requests: 2,
            max_legacy_sessions_per_workspace: 4,
            max_body_bytes: 12 * 1024 * 1024,
            body_timeout: Duration::from_secs(30),
            idle_timeout: Duration::from_secs(30 * 60),
            response_lifetime: Duration::from_secs(30 * 60),
        }
    }
}
impl HostedMcpLimits {
    fn validate(&self) -> Result<()> {
        if self.max_workspaces == 0
            || self.max_workspaces > 64
            || self.max_workspaces_per_member == 0
            || self.max_workspaces_per_member > self.max_workspaces
            || self.max_concurrent_requests == 0
            || self.max_concurrent_requests > 64
            || self.max_concurrent_workspace_requests == 0
            || self.max_concurrent_workspace_requests > self.max_concurrent_requests
            || self.max_legacy_sessions_per_workspace == 0
            || self.max_legacy_sessions_per_workspace > 8
            || self.max_body_bytes == 0
            || self.max_body_bytes > 16 * 1024 * 1024
            || self.body_timeout.is_zero()
            || self.body_timeout > Duration::from_secs(120)
            || self.idle_timeout.is_zero()
            || self.idle_timeout > Duration::from_secs(24 * 60 * 60)
            || self.response_lifetime.is_zero()
            || self.response_lifetime > Duration::from_secs(60 * 60)
        {
            return Err(Error::Configuration("hosted MCP limits"));
        }
        Ok(())
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct WorkspaceKey {
    lab: String,
    member: String,
    label: String,
}
impl WorkspaceKey {
    fn new(lab: &str, member: &str, label: &str) -> std::result::Result<Self, StatusCode> {
        if label.is_empty()
            || label.len() > MAX_WORKSPACE_LABEL
            || !label
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_'))
        {
            return Err(StatusCode::BAD_REQUEST);
        }
        Ok(Self {
            lab: lab.into(),
            member: member.into(),
            label: label.into(),
        })
    }
}

type HttpService = StreamableHttpService<CollageneServer, LocalSessionManager>;
struct Workspace {
    engine: CollageneServer,
    service: HttpService,
    manager: Arc<LocalSessionManager>,
    initialization: tokio::sync::Mutex<()>,
    permits: Arc<Semaphore>,
    active: AtomicUsize,
    last_used: Mutex<Instant>,
}
impl Workspace {
    fn starts_protocol_session(message: &serde_json::Value) -> bool {
        // The SDK routes every initialize through its legacy session lifecycle,
        // even when params names the modern stateless protocol version.
        message.get("method").and_then(|value| value.as_str()) == Some("initialize")
    }
    fn new(root: &StorageRoot, limits: &HostedMcpLimits, origin: &CanonicalOrigin) -> Result<Self> {
        let engine = CollageneServer::with_storage_policy(root.clone())?;
        let server = engine.clone();
        let mut manager = LocalSessionManager::default();
        manager.session_config.channel_capacity = 8;
        manager.session_config.init_timeout = Some(Duration::from_secs(30));
        let manager = Arc::new(manager);
        let mut config = StreamableHttpServerConfig::default();
        let authority = origin
            .0
            .split_once("://")
            .ok_or(Error::Configuration("MCP origin"))?
            .1;
        config.allowed_hosts = vec![authority.into()];
        config.allowed_origins = vec![origin.0.clone()];
        config.max_request_body_bytes = limits.max_body_bytes;
        config.json_response = true;
        config.stateless_protocol_metadata_required = true;
        config.sse_keep_alive = Some(Duration::from_secs(15));
        Ok(Self {
            engine,
            service: HttpService::new(move || Ok(server.clone()), manager.clone(), config),
            manager,
            initialization: tokio::sync::Mutex::new(()),
            permits: Arc::new(Semaphore::new(limits.max_concurrent_workspace_requests)),
            active: AtomicUsize::new(0),
            last_used: Mutex::new(Instant::now()),
        })
    }
    fn close(&self) {
        self.engine.close_compute();
        self.service.config.cancellation_token.cancel();
    }
    fn is_quiescent(&self) -> bool {
        self.active.load(Ordering::Acquire) == 0 && !self.engine.has_pending_compute()
    }
    fn is_idle(&self, now: Instant, idle: Duration) -> bool {
        self.is_quiescent()
            && self
                .last_used
                .lock()
                .map_or(true, |last| now.saturating_duration_since(*last) >= idle)
    }
}
impl Drop for Workspace {
    fn drop(&mut self) {
        self.close();
    }
}

struct Lease {
    workspace: Arc<Workspace>,
    _global_permit: OwnedSemaphorePermit,
    _workspace_permit: OwnedSemaphorePermit,
}
impl Drop for Lease {
    fn drop(&mut self) {
        if let Ok(mut last) = self.workspace.last_used.lock() {
            *last = Instant::now();
        }
        self.workspace.active.fetch_sub(1, Ordering::Release);
    }
}

struct RetiredWorkspace {
    key: WorkspaceKey,
    workspace: Arc<Workspace>,
}
#[derive(Default)]
struct WorkspacePool {
    current: BTreeMap<WorkspaceKey, Arc<Workspace>>,
    retired: Vec<RetiredWorkspace>,
}
impl WorkspacePool {
    fn len(&self) -> usize {
        self.current.len() + self.retired.len()
    }
    fn member_count(&self, member: &str) -> usize {
        self.current
            .keys()
            .filter(|key| key.member == member)
            .count()
            + self
                .retired
                .iter()
                .filter(|entry| entry.key.member == member)
                .count()
    }
    fn retire(
        &mut self,
        key: &WorkspaceKey,
        expected: &Arc<Workspace>,
    ) -> std::result::Result<(), StatusCode> {
        let current = self.current.get(key).ok_or(StatusCode::NOT_FOUND)?;
        if !Arc::ptr_eq(current, expected) {
            return Err(StatusCode::CONFLICT);
        }
        let workspace = self.current.remove(key).ok_or(StatusCode::NOT_FOUND)?;
        workspace.close();
        self.retired.push(RetiredWorkspace {
            key: key.clone(),
            workspace,
        });
        Ok(())
    }
    fn sweep(
        &mut self,
        now: Instant,
        limits: &HostedMcpLimits,
        members: &MemberRegistry,
        shutdown: bool,
    ) {
        let retired = &mut self.retired;
        self.current.retain(|key, workspace| {
            let remove = shutdown
                || !members.contains_actor(&key.member)
                || workspace.is_idle(now, limits.idle_timeout);
            if remove {
                workspace.close();
                retired.push(RetiredWorkspace {
                    key: key.clone(),
                    workspace: workspace.clone(),
                });
            }
            !remove
        });
        // Dropping an HTTP/RPC future does not stop spawn_blocking. Retired
        // engines keep both quota and storage until that native permit returns.
        self.retired.retain(|entry| !entry.workspace.is_quiescent());
    }
}

struct Inner {
    root: StorageRoot,
    limits: HostedMcpLimits,
    members: MemberRegistry,
    lab: String,
    origin: CanonicalOrigin,
    workspaces: Mutex<WorkspacePool>,
    permits: Arc<Semaphore>,
    sweeping: AtomicBool,
    closed: AtomicBool,
}
impl Inner {
    fn sweep(&self, now: Instant) {
        if let Ok(mut workspaces) = self.workspaces.lock() {
            workspaces.sweep(
                now,
                &self.limits,
                &self.members,
                self.closed.load(Ordering::Acquire),
            );
        }
    }
    fn shutdown(&self) {
        self.closed.store(true, Ordering::Release);
        self.sweep(Instant::now());
    }
    fn lease(
        &self,
        key: WorkspaceKey,
        create: bool,
        global: OwnedSemaphorePermit,
    ) -> std::result::Result<Lease, StatusCode> {
        let mut workspaces = self
            .workspaces
            .lock()
            .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
        if self.closed.load(Ordering::Acquire) {
            return Err(StatusCode::SERVICE_UNAVAILABLE);
        }
        if !workspaces.current.contains_key(&key) {
            if !create {
                return Err(StatusCode::NOT_FOUND);
            }
            if workspaces.len() >= self.limits.max_workspaces
                || workspaces.member_count(&key.member) >= self.limits.max_workspaces_per_member
            {
                return Err(StatusCode::TOO_MANY_REQUESTS);
            }
            let workspace = Workspace::new(&self.root, &self.limits, &self.origin)
                .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
            workspaces.current.insert(key.clone(), Arc::new(workspace));
        }
        let workspace = workspaces
            .current
            .get(&key)
            .ok_or(StatusCode::NOT_FOUND)?
            .clone();
        let permit = workspace
            .permits
            .clone()
            .try_acquire_owned()
            .map_err(|_| StatusCode::TOO_MANY_REQUESTS)?;
        workspace.active.fetch_add(1, Ordering::AcqRel);
        Ok(Lease {
            workspace,
            _global_permit: global,
            _workspace_permit: permit,
        })
    }
    fn remove(
        &self,
        key: &WorkspaceKey,
        expected: &Arc<Workspace>,
    ) -> std::result::Result<(), StatusCode> {
        let mut workspaces = self
            .workspaces
            .lock()
            .map_err(|_| StatusCode::SERVICE_UNAVAILABLE)?;
        workspaces.retire(key, expected)
    }
}
impl Drop for Inner {
    fn drop(&mut self) {
        self.shutdown();
    }
}

#[derive(Clone)]
pub(crate) struct HostedMcp(Arc<Inner>);
impl HostedMcp {
    pub(crate) fn new(
        root: StorageRoot,
        limits: HostedMcpLimits,
        members: MemberRegistry,
        lab: String,
        origin: CanonicalOrigin,
    ) -> Result<Self> {
        limits.validate()?;
        Ok(Self(Arc::new(Inner {
            root,
            permits: Arc::new(Semaphore::new(limits.max_concurrent_requests)),
            limits,
            members,
            lab,
            origin,
            workspaces: Mutex::new(WorkspacePool::default()),
            sweeping: AtomicBool::new(false),
            closed: AtomicBool::new(false),
        })))
    }
    pub(crate) fn shutdown(&self) {
        self.0.shutdown();
    }
    fn start_sweeper(&self) {
        if self.0.sweeping.swap(true, Ordering::AcqRel) {
            return;
        }
        let weak = Arc::downgrade(&self.0);
        tokio::spawn(async move {
            loop {
                tokio::time::sleep(SWEEP_INTERVAL).await;
                let Some(inner) = weak.upgrade() else {
                    break;
                };
                inner.sweep(Instant::now());
                if inner.closed.load(Ordering::Acquire) {
                    break;
                }
            }
        });
    }
    fn admit(
        &self,
        request: &Request<Body>,
    ) -> std::result::Result<(HeaderValue, String, WorkspaceKey), StatusCode> {
        let headers = request.headers();
        let mut authorization = headers.get_all(header::AUTHORIZATION).iter();
        let credential = authorization
            .next()
            .ok_or(StatusCode::UNAUTHORIZED)?
            .clone();
        if authorization.next().is_some() {
            return Err(StatusCode::UNAUTHORIZED);
        }
        let member = self
            .0
            .members
            .authenticate(credential.as_bytes())
            .ok_or(StatusCode::UNAUTHORIZED)?;
        if !self.0.origin.accepts(headers) {
            return Err(StatusCode::FORBIDDEN);
        }
        expected_lab_pin(headers, &self.0.lab, true)?;
        if request.uri().query().is_some()
            || request.uri().authority().is_some()
            || request.uri().scheme().is_some()
        {
            return Err(StatusCode::BAD_REQUEST);
        }
        let label = request
            .uri()
            .path()
            .strip_prefix("/api/mcp/")
            .ok_or(StatusCode::BAD_REQUEST)?;
        let key = WorkspaceKey::new(&self.0.lab, member.author.id.as_str(), label)?;
        if !matches!(
            *request.method(),
            Method::GET | Method::POST | Method::DELETE
        ) {
            return Err(StatusCode::METHOD_NOT_ALLOWED);
        }
        if headers.get_all(SESSION_HEADER).iter().count() > 1
            || headers.get(SESSION_HEADER).is_some_and(|v| {
                v.is_empty()
                    || v.len() > 128
                    || !v.as_bytes().iter().all(|b| (0x21..=0x7e).contains(b))
            })
        {
            return Err(StatusCode::BAD_REQUEST);
        }
        if headers.get_all(header::CONTENT_ENCODING).iter().count() > 1
            || headers
                .get(header::CONTENT_ENCODING)
                .is_some_and(|v| v.as_bytes() != b"identity")
        {
            return Err(StatusCode::UNSUPPORTED_MEDIA_TYPE);
        }
        if request.method() == Method::POST && !json_content_type(headers) {
            return Err(StatusCode::UNSUPPORTED_MEDIA_TYPE);
        }
        Ok((credential, member.author.id.as_str().into(), key))
    }
    async fn serve(&self, request: Request<Body>) -> Response {
        let (credential, actor, key) = match self.admit(&request) {
            Ok(v) => v,
            Err(s) => {
                let mut response = rejection(s);
                if s == StatusCode::METHOD_NOT_ALLOWED {
                    response
                        .headers_mut()
                        .insert(header::ALLOW, HeaderValue::from_static("GET, POST, DELETE"));
                }
                return response;
            }
        };
        self.start_sweeper();
        self.0.sweep(Instant::now());
        let global = match self.0.permits.clone().try_acquire_owned() {
            Ok(p) => p,
            Err(_) => return rejection(StatusCode::TOO_MANY_REQUESTS),
        };
        let (mut parts, body) = request.into_parts();
        let has_session = parts.headers.contains_key(SESSION_HEADER);
        // Acquire scope and global leases before buffering any body. Four scopes
        // with two requests each also bound a member's aggregate ingress.
        let lease = match self.0.lease(
            key.clone(),
            parts.method == Method::POST && !has_session,
            global,
        ) {
            Ok(l) => l,
            Err(s) => return rejection(s),
        };
        let request_deadline = tokio::time::Instant::now() + self.0.limits.response_lifetime;
        let limit = if parts.method == Method::POST {
            self.0.limits.max_body_bytes
        } else {
            0
        };
        let body_deadline =
            request_deadline.min(tokio::time::Instant::now() + self.0.limits.body_timeout);
        let bytes = match tokio::time::timeout_at(body_deadline, to_bytes(body, limit)).await {
            Err(_) => return rejection(StatusCode::REQUEST_TIMEOUT),
            Ok(Err(_)) => return rejection(StatusCode::PAYLOAD_TOO_LARGE),
            Ok(Ok(bytes)) => bytes,
        };
        // A queued or slow body cannot retain revoked membership. Never forward
        // a bearer credential into SDK request metadata or scientific handlers.
        if self
            .0
            .members
            .authenticate(credential.as_bytes())
            .is_none_or(|m| m.author.id.as_str() != actor)
        {
            return rejection(StatusCode::UNAUTHORIZED);
        }
        parts.headers.remove(header::AUTHORIZATION);
        if parts.method == Method::DELETE && !has_session {
            return match self.0.remove(&key, &lease.workspace) {
                Ok(()) => StatusCode::NO_CONTENT.into_response(),
                Err(s) => rejection(s),
            };
        }
        let legacy_initialize = if parts.method == Method::POST {
            let Ok(value) = serde_json::from_slice::<serde_json::Value>(&bytes) else {
                return rejection(StatusCode::BAD_REQUEST);
            };
            Workspace::starts_protocol_session(&value)
        } else {
            false
        };
        let workspace = lease.workspace.clone();
        // Serialize legacy handshakes before session-count admission. Modern
        // stateless requests use the same explicitly named scientific workspace.
        let initializing = if legacy_initialize && !has_session {
            match tokio::time::timeout_at(request_deadline, workspace.initialization.lock()).await {
                Ok(guard) => Some(guard),
                Err(_) => return rejection(StatusCode::GATEWAY_TIMEOUT),
            }
        } else {
            None
        };
        if initializing.is_some()
            && workspace.manager.sessions.read().await.len()
                >= self.0.limits.max_legacy_sessions_per_workspace
        {
            return rejection(StatusCode::TOO_MANY_REQUESTS);
        }
        if self
            .0
            .members
            .authenticate(credential.as_bytes())
            .is_none_or(|m| m.author.id.as_str() != actor)
            || workspace.service.config.cancellation_token.is_cancelled()
        {
            return rejection(StatusCode::UNAUTHORIZED);
        }
        drop(credential);
        let response = match tokio::time::timeout_at(
            request_deadline,
            workspace
                .service
                .handle(Request::from_parts(parts, Body::from(bytes))),
        )
        .await
        {
            Ok(response) => response,
            Err(_) => return rejection(StatusCode::GATEWAY_TIMEOUT),
        };
        drop(initializing);
        let (mut parts, body) = response.into_parts();
        parts
            .headers
            .insert(header::CACHE_CONTROL, HeaderValue::from_static("no-store"));
        let deadline = Box::pin(tokio::time::sleep_until(request_deadline));
        let stream = futures_util::stream::unfold(
            (Body::new(body).into_data_stream(), lease, deadline),
            |(mut stream, lease, mut deadline)| async move {
                let cancellation = lease.workspace.service.config.cancellation_token.clone();
                tokio::select! {
                    _ = &mut deadline => None,
                    _ = cancellation.cancelled() => None,
                    frame = stream.next() => frame.map(|frame| (frame, (stream, lease, deadline))),
                }
            },
        );
        Response::from_parts(parts, Body::from_stream(stream))
    }
}
fn json_content_type(headers: &HeaderMap) -> bool {
    headers.get_all(header::CONTENT_TYPE).iter().count() == 1
        && headers
            .get(header::CONTENT_TYPE)
            .and_then(|v| v.to_str().ok())
            .is_some_and(|v| {
                v.split(';')
                    .next()
                    .is_some_and(|v| v.trim().eq_ignore_ascii_case("application/json"))
            })
}
pub(crate) async fn request(State(state): State<ServerState>, request: Request<Body>) -> Response {
    let Some(hosted) = &state.0.hosted_mcp else {
        return rejection(StatusCode::NOT_FOUND);
    };
    hosted.serve(request).await
}

#[cfg(test)]
mod tests;
