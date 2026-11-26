use crate::wire::{self, ErrorCode, Failure, Operation, Output};
use collagene_core::formats::ImportReport;
use collagene_sources::{
    AddgeneCatalog, Config as SourceConfig, Source, SourceClient, SourceError, parse_source,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    future::Future,
    net::{IpAddr, Ipv4Addr, SocketAddr},
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};
use tokio::sync::{Semaphore, watch};
use url::Url;

#[derive(Clone)]
pub struct Cancellation {
    state: watch::Sender<bool>,
}
impl Default for Cancellation {
    fn default() -> Self {
        Self {
            state: watch::channel(false).0,
        }
    }
}
impl Cancellation {
    pub fn cancel(&self) {
        self.state.send_replace(true);
    }
    pub fn is_cancelled(&self) -> bool {
        *self.state.borrow()
    }
    pub async fn cancelled(&self) {
        let mut state = self.state.subscribe();
        loop {
            if *state.borrow_and_update() {
                return;
            }
            if state.changed().await.is_err() {
                return;
            }
        }
    }
}
#[derive(Clone, Debug)]
pub struct Limits {
    pub concurrent_imports: usize,
    pub connections: usize,
    pub retained_catalogs: usize,
    pub request_timeout: Duration,
    pub catalog_lifetime: Duration,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            concurrent_imports: 2,
            connections: 4,
            retained_catalogs: 2,
            request_timeout: Duration::from_secs(70),
            catalog_lifetime: Duration::from_secs(120),
        }
    }
}
#[derive(Debug, thiserror::Error)]
#[error("Invalid source-gateway host policy or resource limits")]
pub struct ConfigError;
/// No Debug implementation: the optional gateway capability is secret.
#[derive(Clone)]
pub struct Config {
    pub(crate) bind: SocketAddr,
    pub(crate) origins: BTreeSet<String>,
    pub(crate) authorities: BTreeSet<String>,
    pub(crate) access_digest: Option<[u8; 32]>,
    pub(crate) limits: Limits,
}
impl Default for Config {
    fn default() -> Self {
        Self::loopback(SocketAddr::new(IpAddr::V4(Ipv4Addr::LOCALHOST), 8084))
            .expect("fixed loopback policy")
    }
}
impl Config {
    pub fn loopback(bind: SocketAddr) -> Result<Self, ConfigError> {
        if !bind.ip().is_loopback() || bind.port() == 0 {
            return Err(ConfigError);
        }
        let origin = format!("http://{bind}");
        Self::configured(bind, vec![origin], None, false, Limits::default())
    }
    /// Configure exact frontend origins for a loopback host or reverse proxy.
    /// A public owning host must instead use remote(), even when it mounts the
    /// route in-process instead of invoking this crate's standalone listener.
    pub fn loopback_origins(
        bind: SocketAddr,
        origins: Vec<String>,
        capability: Option<String>,
    ) -> Result<Self, ConfigError> {
        if !bind.ip().is_loopback() {
            return Err(ConfigError);
        }
        Self::configured(bind, origins, capability, false, Limits::default())
    }
    /// Explicit non-loopback host policy. TLS termination belongs to the serving
    /// host/proxy; permitted browser origins here must use HTTPS for remote bind.
    /// No forwarded-origin header is trusted by this crate.
    pub fn remote(
        bind: SocketAddr,
        origins: Vec<String>,
        capability: String,
    ) -> Result<Self, ConfigError> {
        Self::configured(bind, origins, Some(capability), true, Limits::default())
    }
    pub fn with_limits(mut self, limits: Limits) -> Result<Self, ConfigError> {
        validate_limits(&limits)?;
        self.limits = limits;
        Ok(self)
    }
    pub fn bind_address(&self) -> SocketAddr {
        self.bind
    }
    pub fn allowed_origins(&self) -> &BTreeSet<String> {
        &self.origins
    }
    fn configured(
        bind: SocketAddr,
        origins: Vec<String>,
        capability: Option<String>,
        remote: bool,
        limits: Limits,
    ) -> Result<Self, ConfigError> {
        validate_limits(&limits)?;
        if bind.port() == 0
            || (!bind.ip().is_loopback() && !remote)
            || origins.is_empty()
            || origins.len() > 16
        {
            return Err(ConfigError);
        }
        let digest = match capability {
            Some(value)
                if (32..=256).contains(&value.len())
                    && value.bytes().all(|b| b.is_ascii_graphic()) =>
            {
                Some(*blake3::hash(value.as_bytes()).as_bytes())
            }
            Some(_) => return Err(ConfigError),
            None => None,
        };
        if remote && digest.is_none() {
            return Err(ConfigError);
        }
        let mut origins_set = BTreeSet::new();
        let mut authorities = BTreeSet::new();
        for origin in origins {
            let (canonical, authority) = origin_parts(&origin).ok_or(ConfigError)?;
            if canonical != origin || (remote && !origin.starts_with("https://")) {
                return Err(ConfigError);
            }
            origins_set.insert(canonical);
            authorities.insert(authority);
        }
        Ok(Self {
            bind,
            origins: origins_set,
            authorities,
            access_digest: digest,
            limits,
        })
    }
}
fn validate_limits(l: &Limits) -> Result<(), ConfigError> {
    if !(1..=8).contains(&l.concurrent_imports)
        || !(1..=16).contains(&l.connections)
        || !(1..=4).contains(&l.retained_catalogs)
        || l.request_timeout < Duration::from_secs(1)
        || l.request_timeout > Duration::from_secs(90)
        || l.catalog_lifetime < Duration::from_secs(1)
        || l.catalog_lifetime > Duration::from_secs(300)
    {
        Err(ConfigError)
    } else {
        Ok(())
    }
}
pub(crate) fn origin_parts(value: &str) -> Option<(String, String)> {
    if value.len() > 2048 {
        return None;
    }
    let url = Url::parse(value).ok()?;
    if !matches!(url.scheme(), "http" | "https")
        || url.host_str().is_none()
        || !url.username().is_empty()
        || url.password().is_some()
        || url.query().is_some()
        || url.fragment().is_some()
        || url.path() != "/"
    {
        return None;
    }
    Some((
        url.origin().ascii_serialization(),
        url[url::Position::BeforeHost..url::Position::AfterPort].to_ascii_lowercase(),
    ))
}
/// Injectable only at native host construction; browser requests cannot replace
/// this implementation. NativeBackend delegates all DNS/redirect/parser policy.
pub trait Backend: Send + Sync + 'static {
    fn import(
        &self,
        source: String,
        token: Option<String>,
    ) -> impl Future<Output = Result<ImportReport, SourceError>> + Send;
    fn catalog(
        &self,
        source: String,
        token: Option<String>,
    ) -> impl Future<Output = Result<AddgeneCatalog, SourceError>> + Send;
    fn selection(
        &self,
        catalog: Arc<AddgeneCatalog>,
        id: u64,
        token: Option<String>,
    ) -> impl Future<Output = Result<ImportReport, SourceError>> + Send;
}
#[derive(Default)]
pub struct NativeBackend;
fn client(token: Option<String>) -> Result<SourceClient, SourceError> {
    let config = match token {
        Some(token) => SourceConfig::default().with_addgene_api_token(token)?,
        None => SourceConfig::default(),
    };
    SourceClient::new(config)
}
impl Backend for NativeBackend {
    async fn import(
        &self,
        source: String,
        token: Option<String>,
    ) -> Result<ImportReport, SourceError> {
        client(token)?.import_url(&source).await
    }
    async fn catalog(
        &self,
        source: String,
        token: Option<String>,
    ) -> Result<AddgeneCatalog, SourceError> {
        client(token)?.addgene_catalog(&source).await
    }
    async fn selection(
        &self,
        catalog: Arc<AddgeneCatalog>,
        id: u64,
        token: Option<String>,
    ) -> Result<ImportReport, SourceError> {
        client(token)?.import_addgene_selection(&catalog, id).await
    }
}
type Owner = [u8; 32];
type Handle = [u8; 32];
struct CatalogEntry {
    owner: Owner,
    expires: Instant,
    catalog: Arc<AddgeneCatalog>,
}
#[derive(Default)]
struct State {
    active: BTreeMap<(Owner, u64), Cancellation>,
    catalogs: BTreeMap<Handle, CatalogEntry>,
    cleanup_running: bool,
}
// One weakly-owned timer per gateway, not one task per catalog. Newly inserted
// catalogs have the same lifetime and therefore cannot precede the oldest one.
fn start_catalog_cleanup(state: &Arc<Mutex<State>>) {
    let state = Arc::downgrade(state);
    tokio::spawn(async move {
        loop {
            let deadline = {
                let Some(state) = state.upgrade() else {
                    return;
                };
                let Ok(mut state) = state.lock() else {
                    return;
                };
                state
                    .catalogs
                    .retain(|_, entry| entry.expires > Instant::now());
                match state.catalogs.values().map(|entry| entry.expires).min() {
                    Some(deadline) => deadline,
                    None => {
                        state.cleanup_running = false;
                        return;
                    }
                }
            };
            tokio::time::sleep_until(tokio::time::Instant::from_std(deadline)).await;
        }
    });
}
struct Inner<B> {
    pub(crate) config: Config,
    backend: Arc<B>,
    permits: Arc<Semaphore>,
    state: Arc<Mutex<State>>,
}
pub struct Gateway<B: Backend = NativeBackend> {
    inner: Arc<Inner<B>>,
}
impl<B: Backend> Clone for Gateway<B> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}
struct Active {
    state: Arc<Mutex<State>>,
    key: (Owner, u64),
}
impl Drop for Active {
    fn drop(&mut self) {
        if let Ok(mut state) = self.state.lock() {
            state.active.remove(&self.key);
        }
    }
}
struct AbortOnDrop<T>(tokio::task::JoinHandle<T>);
impl<T> Drop for AbortOnDrop<T> {
    fn drop(&mut self) {
        self.0.abort();
    }
}
enum Work {
    Import {
        source: String,
        token: Option<String>,
    },
    Catalog {
        source: String,
        token: Option<String>,
    },
    Select {
        catalog: Arc<AddgeneCatalog>,
        id: u64,
        token: Option<String>,
    },
}
enum Finished {
    Import(ImportReport),
    Catalog(AddgeneCatalog),
}
fn invalid() -> Failure {
    Failure::new(
        ErrorCode::InvalidRequest,
        "Unsupported source-import request.",
    )
}
fn internal() -> Failure {
    Failure::new(
        ErrorCode::Internal,
        "The source-import gateway could not complete this request.",
    )
}
fn provider(error: SourceError) -> Failure {
    let retry = match &error {
        SourceError::RateLimited {
            retry_after_seconds,
        }
        | SourceError::HttpStatus {
            retry_after_seconds,
            ..
        } => *retry_after_seconds,
        _ => None,
    };
    let code = match &error {
        SourceError::TooLarge => ErrorCode::TooLarge,
        SourceError::Timeout => ErrorCode::Timeout,
        _ => ErrorCode::Provider,
    };
    // SourceError has intentionally sanitized messages; it never embeds a URL,
    // token, network exception or provider response body. Still bound the preview.
    Failure {
        code,
        message: error.to_string().chars().take(512).collect(),
        retry_after_seconds: retry,
    }
}
fn approved(source: &str, catalog: bool) -> Result<(), Failure> {
    if source.len() > 8192 {
        return Err(invalid());
    }
    match parse_source(source).map_err(provider)? {
        Source::Addgene { .. } if catalog => Ok(()),
        Source::Ncbi { .. } | Source::Ena { .. } | Source::Ddbj { .. } if !catalog => Ok(()),
        // Generic Direct redirects are deliberately not admitted: only the
        // source client can enforce its complete redirect boundary.
        _ => Err(invalid()),
    }
}
fn token_guard(token: &Option<String>) -> Result<(), Failure> {
    if token
        .as_ref()
        .is_some_and(|s| s.is_empty() || s.len() > 4096 || !s.bytes().all(|b| b.is_ascii_graphic()))
    {
        Err(invalid())
    } else {
        Ok(())
    }
}
impl Gateway<NativeBackend> {
    pub fn native(config: Config) -> Self {
        Self::with_backend(config, NativeBackend)
    }
}
impl<B: Backend> Gateway<B> {
    pub fn with_backend(config: Config, backend: B) -> Self {
        let permits = Arc::new(Semaphore::new(config.limits.concurrent_imports));
        Self {
            inner: Arc::new(Inner {
                config,
                backend: Arc::new(backend),
                permits,
                state: Arc::new(Mutex::new(State::default())),
            }),
        }
    }
    pub fn config(&self) -> &Config {
        &self.inner.config
    }
    /// Trusted native/Tauri host entry point. HTTP callers must use handle_http,
    /// which checks host, origin, media type and gateway access capability first.
    pub async fn execute(&self, request: wire::Request, external: Cancellation) -> wire::Response {
        let id = request.id;
        let result = self.execute_inner(request, external).await;
        wire::Response { id, result }
    }
    async fn execute_inner(
        &self,
        request: wire::Request,
        external: Cancellation,
    ) -> Result<Output, Failure> {
        if request.client == [0; 32] {
            return Err(invalid());
        }
        if external.is_cancelled() {
            return Err(Failure::new(
                ErrorCode::Canceled,
                "Source request canceled.",
            ));
        }
        let owner = *blake3::hash(&request.client).as_bytes();
        let key = (owner, request.id);
        let work = match request.operation {
            Operation::Cancel { request_id } => {
                let state = self.inner.state.lock().map_err(|_| internal())?;
                let found = state.active.get(&(owner, request_id));
                if let Some(cancellation) = found {
                    cancellation.cancel();
                }
                return Ok(Output::Canceled {
                    found: found.is_some(),
                });
            }
            Operation::DisposeCatalog { handle } => {
                let mut state = self.inner.state.lock().map_err(|_| internal())?;
                let found = state
                    .catalogs
                    .get(&handle)
                    .is_some_and(|e| e.owner == owner);
                if found {
                    state.catalogs.remove(&handle);
                }
                return Ok(Output::CatalogDisposed { found });
            }
            Operation::Import {
                source,
                addgene_token,
            } => {
                approved(&source, false)?;
                token_guard(&addgene_token)?;
                if addgene_token.is_some() {
                    return Err(invalid());
                }
                Work::Import {
                    source,
                    token: None,
                }
            }
            Operation::Catalog {
                source,
                addgene_token,
            } => {
                approved(&source, true)?;
                token_guard(&addgene_token)?;
                Work::Catalog {
                    source,
                    token: addgene_token,
                }
            }
            Operation::Select {
                handle,
                response_blake3,
                sequence_id,
                addgene_token,
            } => {
                token_guard(&addgene_token)?;
                if sequence_id == 0 || response_blake3.len() != 64 {
                    return Err(invalid());
                }
                let mut state = self.inner.state.lock().map_err(|_| internal())?;
                state.catalogs.retain(|_, e| e.expires > Instant::now());
                let entry = state
                    .catalogs
                    .get(&handle)
                    .filter(|e| e.owner == owner && e.catalog.response_blake3() == response_blake3)
                    .ok_or_else(|| {
                        Failure::new(
                            ErrorCode::ExpiredCatalog,
                            "The reviewed catalog is unavailable. Refresh and select again.",
                        )
                    })?;
                // No user-provided signed download URL or fabricated catalog can
                // reach the native selection client.
                entry.catalog.download_url(sequence_id).map_err(provider)?;
                Work::Select {
                    catalog: entry.catalog.clone(),
                    id: sequence_id,
                    token: addgene_token,
                }
            }
        };
        let runtime = tokio::runtime::Handle::try_current().map_err(|_| internal())?;
        let permit = self
            .inner
            .permits
            .clone()
            .try_acquire_owned()
            .map_err(|_| {
                Failure::new(
                    ErrorCode::Busy,
                    "Source-import capacity is in use. Retry after an active request finishes.",
                )
            })?;
        let cancellation = Cancellation::default();
        {
            let mut state = self.inner.state.lock().map_err(|_| internal())?;
            if state.active.contains_key(&key) {
                return Err(Failure::new(
                    ErrorCode::Busy,
                    "This client request ID is already active.",
                ));
            }
            state.active.insert(key, cancellation.clone());
        }
        let _active = Active {
            state: self.inner.state.clone(),
            key,
        };
        let backend = self.inner.backend.clone();
        let started = Instant::now();
        let mut task = AbortOnDrop(runtime.spawn(async move {
            // Cancellation drops network/rate-wait futures. A native parser that
            // is already executing synchronously finishes within its own bounded
            // input limits; its permit is not released early by HTTP cancellation.
            let _permit = permit;
            match work {
                Work::Import { source, token } => {
                    backend.import(source, token).await.map(Finished::Import)
                }
                Work::Catalog { source, token } => {
                    backend.catalog(source, token).await.map(Finished::Catalog)
                }
                Work::Select { catalog, id, token } => backend
                    .selection(catalog, id, token)
                    .await
                    .map(Finished::Import),
            }
        }));
        let finished = tokio::select! {
            result=&mut task.0=>result.map_err(|_|internal())?.map_err(provider)?,
            _=cancellation.cancelled()=>return Err(Failure::new(ErrorCode::Canceled,"Source request canceled.")),
            _=external.cancelled()=>return Err(Failure::new(ErrorCode::Canceled,"Source request canceled.")),
            _=tokio::time::sleep(self.inner.config.limits.request_timeout)=>return Err(Failure::new(ErrorCode::Timeout,"Source request exceeded the gateway deadline.")),
        };
        // Tokio deadlines cannot preempt synchronous parsing. Even if both the
        // native result and timer became ready before this task was scheduled,
        // a late result must not be reported as an on-time successful import.
        if started.elapsed() > self.inner.config.limits.request_timeout {
            return Err(Failure::new(
                ErrorCode::Timeout,
                "Source request exceeded the gateway deadline.",
            ));
        }
        if cancellation.is_cancelled() || external.is_cancelled() {
            return Err(Failure::new(
                ErrorCode::Canceled,
                "Source request canceled.",
            ));
        }
        match finished {
            Finished::Import(report) => {
                if report.molecules.len() != 1
                    || report.warnings.len() > 2048
                    || report.warnings.iter().any(|w| w.len() > 16384)
                {
                    return Err(Failure::new(
                        ErrorCode::TooLarge,
                        "Source import result exceeds gateway limits.",
                    ));
                }
                collagene_core::model::bounds::model_preflight(
                    &report.molecules.iter().collect::<Vec<_>>(),
                )
                .map_err(|_| {
                    Failure::new(
                        ErrorCode::TooLarge,
                        "Source import result exceeds gateway limits.",
                    )
                })?;
                report.molecules[0].validate().map_err(|_| {
                    Failure::new(
                        ErrorCode::Provider,
                        "The imported source document is invalid.",
                    )
                })?;
                Ok(Output::Imported(report))
            }
            Finished::Catalog(catalog) => {
                if catalog.choices().len() > 1000 {
                    return Err(Failure::new(
                        ErrorCode::TooLarge,
                        "Source catalog exceeds gateway limits.",
                    ));
                }
                let mut handle = [0; 32];
                getrandom::fill(&mut handle).map_err(|_| internal())?;
                let summary = wire::Catalog {
                    handle,
                    response_blake3: catalog.response_blake3().into(),
                    plasmid_id: catalog.plasmid_id(),
                    name: catalog.name().into(),
                    choices: catalog
                        .choices()
                        .iter()
                        .map(|c| wire::Choice {
                            sequence_id: c.sequence_id,
                            description: c.description.clone(),
                            origin: match c.origin {
                                collagene_sources::AddgeneSequenceOrigin::Depositor => {
                                    wire::Origin::Depositor
                                }
                                collagene_sources::AddgeneSequenceOrigin::Addgene => {
                                    wire::Origin::Addgene
                                }
                            },
                            completeness: match c.completeness {
                                collagene_sources::AddgeneCompleteness::Full => {
                                    wire::Completeness::Full
                                }
                                collagene_sources::AddgeneCompleteness::Partial => {
                                    wire::Completeness::Partial
                                }
                            },
                            length: c.length,
                            catalog_sequence_blake3: c.catalog_sequence_blake3.clone(),
                        })
                        .collect(),
                    expires_in_seconds: self.inner.config.limits.catalog_lifetime.as_secs(),
                };
                let mut state = self.inner.state.lock().map_err(|_| internal())?;
                state.catalogs.retain(|_, e| e.expires > Instant::now());
                if state.catalogs.len() >= self.inner.config.limits.retained_catalogs
                    || state.catalogs.contains_key(&handle)
                {
                    return Err(Failure::new(
                        ErrorCode::Busy,
                        "Catalog capacity is in use. Dispose an old catalog or wait for expiry.",
                    ));
                }
                let start_cleanup = !state.cleanup_running;
                state.cleanup_running = true;
                state.catalogs.insert(
                    handle,
                    CatalogEntry {
                        owner,
                        expires: Instant::now() + self.inner.config.limits.catalog_lifetime,
                        catalog: Arc::new(catalog),
                    },
                );
                drop(state);
                if start_cleanup {
                    start_catalog_cleanup(&self.inner.state);
                }
                Ok(Output::Catalog(summary))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Molecule, Topology};
    use std::sync::atomic::{AtomicUsize, Ordering};
    const CATALOG: &str = r#"{"id":12345,"name":"Independent gateway catalog","sequences":{"public_user_full_sequences":[{"sequence_id":41,"sequence_description":"original reviewed record","length":16,"sequence":"aaaaccccggggtttt","genbank_api_url":"https://files.example.org/41.gb?signature=TRANSIENT_LINK_SECRET"}],"public_addgene_full_sequences":[],"public_user_partial_sequences":[],"public_addgene_partial_sequences":[]}}"#;
    const GB: &str = "LOCUS       GATEWAY_TEST              16 bp    DNA     circular SYN 01-JAN-2026\nDEFINITION  Independent source gateway fixture.\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaaccccgg ggtttt\n//\n";
    #[derive(Default)]
    struct FixtureBackend {
        catalog_calls: AtomicUsize,
        selection_calls: AtomicUsize,
    }
    impl Backend for Arc<FixtureBackend> {
        async fn import(
            &self,
            _source: String,
            _token: Option<String>,
        ) -> Result<ImportReport, SourceError> {
            Ok(ImportReport {
                molecules: vec![
                    Molecule::dna("Imported original", "ACGT", Topology::Linear).unwrap(),
                ],
                warnings: vec![],
            })
        }
        async fn catalog(
            &self,
            _source: String,
            token: Option<String>,
        ) -> Result<AddgeneCatalog, SourceError> {
            assert_eq!(token.as_deref(), Some("PRIVATE_SOURCE_TOKEN"));
            self.catalog_calls.fetch_add(1, Ordering::Relaxed);
            collagene_sources::parse_addgene_catalog(12345, CATALOG.as_bytes())
        }
        async fn selection(
            &self,
            catalog: Arc<AddgeneCatalog>,
            id: u64,
            token: Option<String>,
        ) -> Result<ImportReport, SourceError> {
            assert_eq!(token.as_deref(), Some("PRIVATE_SOURCE_TOKEN"));
            self.selection_calls.fetch_add(1, Ordering::Relaxed);
            assert!(catalog.download_url(id)?.contains("TRANSIENT_LINK_SECRET"));
            catalog.import_download(id, GB.as_bytes(), "2026-11-26T00:00:00Z")
        }
    }
    fn request(client: u8, id: u64, operation: Operation) -> wire::Request {
        wire::Request {
            id,
            client: [client; 32],
            operation,
        }
    }
    #[tokio::test]
    async fn catalog_selection_keeps_original_native_binding_and_client_scope() {
        let backend = Arc::new(FixtureBackend::default());
        let gateway = Gateway::with_backend(Config::default(), backend.clone());
        let first = gateway
            .execute(
                request(
                    7,
                    1,
                    Operation::Catalog {
                        source: "addgene:12345".into(),
                        addgene_token: Some("PRIVATE_SOURCE_TOKEN".into()),
                    },
                ),
                Cancellation::default(),
            )
            .await;
        let wire = wire::encode_response(&first).unwrap();
        assert!(
            !wire
                .windows("PRIVATE_SOURCE_TOKEN".len())
                .any(|s| s == b"PRIVATE_SOURCE_TOKEN")
        );
        assert!(
            !wire
                .windows("TRANSIENT_LINK_SECRET".len())
                .any(|s| s == b"TRANSIENT_LINK_SECRET")
        );
        let Output::Catalog(catalog) = first.result.unwrap() else {
            panic!("catalog");
        };
        assert_eq!(catalog.choices[0].sequence_id, 41);
        for (client, hash) in [(8, catalog.response_blake3.clone()), (7, "0".repeat(64))] {
            let result = gateway
                .execute(
                    request(
                        client,
                        2,
                        Operation::Select {
                            handle: catalog.handle,
                            response_blake3: hash,
                            sequence_id: 41,
                            addgene_token: Some("PRIVATE_SOURCE_TOKEN".into()),
                        },
                    ),
                    Cancellation::default(),
                )
                .await;
            assert_eq!(result.result.unwrap_err().code, ErrorCode::ExpiredCatalog);
        }
        let selected = gateway
            .execute(
                request(
                    7,
                    3,
                    Operation::Select {
                        handle: catalog.handle,
                        response_blake3: catalog.response_blake3.clone(),
                        sequence_id: 41,
                        addgene_token: Some("PRIVATE_SOURCE_TOKEN".into()),
                    },
                ),
                Cancellation::default(),
            )
            .await;
        let Output::Imported(imported) = selected.result.unwrap() else {
            panic!("import");
        };
        assert_eq!(imported.molecules[0].sequence.as_str(), "AAAACCCCGGGGTTTT");
        assert_eq!(imported.molecules[0].metadata["source:identifier"], "12345");
        assert_eq!(backend.catalog_calls.load(Ordering::Relaxed), 1);
        assert_eq!(backend.selection_calls.load(Ordering::Relaxed), 1);
        let disposed = gateway
            .execute(
                request(
                    7,
                    4,
                    Operation::DisposeCatalog {
                        handle: catalog.handle,
                    },
                ),
                Cancellation::default(),
            )
            .await;
        assert!(matches!(
            disposed.result,
            Ok(Output::CatalogDisposed { found: true })
        ));
        let selected = gateway
            .execute(
                request(
                    7,
                    5,
                    Operation::Select {
                        handle: catalog.handle,
                        response_blake3: catalog.response_blake3,
                        sequence_id: 41,
                        addgene_token: None,
                    },
                ),
                Cancellation::default(),
            )
            .await;
        assert_eq!(selected.result.unwrap_err().code, ErrorCode::ExpiredCatalog);
    }
    struct Blocking {
        started: tokio::sync::Notify,
        active: Arc<AtomicUsize>,
    }
    struct Activity(Arc<AtomicUsize>);
    impl Drop for Activity {
        fn drop(&mut self) {
            self.0.fetch_sub(1, Ordering::Relaxed);
        }
    }
    impl Backend for Arc<Blocking> {
        async fn import(
            &self,
            _source: String,
            _token: Option<String>,
        ) -> Result<ImportReport, SourceError> {
            self.active.fetch_add(1, Ordering::Relaxed);
            let _active = Activity(self.active.clone());
            self.started.notify_one();
            std::future::pending().await
        }
        async fn catalog(
            &self,
            _source: String,
            _token: Option<String>,
        ) -> Result<AddgeneCatalog, SourceError> {
            Err(SourceError::InvalidSource)
        }
        async fn selection(
            &self,
            _catalog: Arc<AddgeneCatalog>,
            _id: u64,
            _token: Option<String>,
        ) -> Result<ImportReport, SourceError> {
            Err(SourceError::InvalidSource)
        }
    }
    #[tokio::test]
    async fn cancellation_is_scoped_and_busy_requests_do_not_queue_native_work() {
        let backend = Arc::new(Blocking {
            started: tokio::sync::Notify::new(),
            active: Arc::new(AtomicUsize::new(0)),
        });
        let config = Config::default()
            .with_limits(Limits {
                concurrent_imports: 1,
                ..Limits::default()
            })
            .unwrap();
        let gateway = Gateway::with_backend(config, backend.clone());
        let running = gateway.clone();
        let task = tokio::spawn(async move {
            running
                .execute(
                    request(
                        7,
                        1,
                        Operation::Import {
                            source: "AB000263.1".into(),
                            addgene_token: None,
                        },
                    ),
                    Cancellation::default(),
                )
                .await
        });
        backend.started.notified().await;
        let busy = gateway
            .execute(
                request(
                    7,
                    2,
                    Operation::Import {
                        source: "AB000263.1".into(),
                        addgene_token: None,
                    },
                ),
                Cancellation::default(),
            )
            .await;
        assert_eq!(busy.result.unwrap_err().code, ErrorCode::Busy);
        let wrong = gateway
            .execute(
                request(8, 3, Operation::Cancel { request_id: 1 }),
                Cancellation::default(),
            )
            .await;
        assert!(matches!(
            wrong.result,
            Ok(Output::Canceled { found: false })
        ));
        let canceled = gateway
            .execute(
                request(7, 4, Operation::Cancel { request_id: 1 }),
                Cancellation::default(),
            )
            .await;
        assert!(matches!(
            canceled.result,
            Ok(Output::Canceled { found: true })
        ));
        let response = tokio::time::timeout(Duration::from_secs(1), task)
            .await
            .unwrap()
            .unwrap();
        assert_eq!(response.result.unwrap_err().code, ErrorCode::Canceled);
        tokio::time::timeout(Duration::from_secs(1), async {
            while backend.active.load(Ordering::Relaxed) != 0 {
                tokio::task::yield_now().await;
            }
        })
        .await
        .unwrap();
    }
    #[test]
    fn remote_listener_requires_explicit_token_and_https_host_policy() {
        let remote = "0.0.0.0:8084".parse().unwrap();
        assert!(Config::loopback(remote).is_err());
        assert!(
            Config::remote(
                remote,
                vec!["https://collagene.example".into()],
                "short".into()
            )
            .is_err()
        );
        assert!(
            Config::remote(
                remote,
                vec!["http://collagene.example".into()],
                "x".repeat(32)
            )
            .is_err()
        );
        assert!(
            Config::remote(
                remote,
                vec!["https://collagene.example".into()],
                "x".repeat(32)
            )
            .is_ok()
        );
    }
}
