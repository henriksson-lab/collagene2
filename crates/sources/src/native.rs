use super::*;
use reqwest::{
    Client,
    dns::{Addrs, Name, Resolve, Resolving},
    redirect::Policy,
};
use std::{
    future::Future,
    io,
    net::SocketAddr,
    sync::Arc,
    time::{Duration, Instant},
};
use tokio::sync::Mutex;

const REQUEST_TIMEOUT: Duration = Duration::from_secs(30);
const MAX_REDIRECTS: usize = 5;
const NCBI_INTERVAL: Duration = Duration::from_millis(350);
struct RateState {
    last: Option<Instant>,
    cooldown: Option<(Instant, u64)>,
}
static NCBI_RATE: Mutex<RateState> = Mutex::const_new(RateState {
    last: None,
    cooldown: None,
});
static ADDGENE_RATE: Mutex<RateState> = Mutex::const_new(RateState {
    last: None,
    cooldown: None,
});

/// Process-wide conservative E-utilities throttle, including separate server sessions.
async fn throttle_ncbi() -> Result<(), SourceError> {
    throttle(&NCBI_RATE, NCBI_INTERVAL).await
}
async fn throttle(rate: &Mutex<RateState>, interval: Duration) -> Result<(), SourceError> {
    let mut state = rate.lock().await;
    if let Some((started, seconds)) = state.cooldown {
        let elapsed = started.elapsed().as_secs();
        if elapsed < seconds {
            return Err(SourceError::RateLimited {
                retry_after_seconds: Some(seconds - elapsed),
            });
        }
        state.cooldown = None;
    }
    if let Some(previous) = state.last {
        let next = previous + interval;
        if next > Instant::now() {
            tokio::time::sleep_until(tokio::time::Instant::from_std(next)).await;
        }
    }
    state.last = Some(Instant::now());
    Ok(())
}
async fn note_cooldown(seconds: u64) {
    record_cooldown(&NCBI_RATE, seconds).await
}
async fn record_cooldown(rate: &Mutex<RateState>, seconds: u64) {
    let mut state = rate.lock().await;
    let remaining = state.cooldown.map_or(0, |(started, old)| {
        old.saturating_sub(started.elapsed().as_secs())
    });
    state.cooldown = Some((Instant::now(), seconds.max(remaining)));
}
fn public_addresses(
    addresses: Vec<SocketAddr>,
) -> Result<Addrs, Box<dyn std::error::Error + Send + Sync>> {
    if addresses.is_empty() || addresses.iter().any(|address| !is_public_ip(address.ip())) {
        return Err(io::Error::other("source DNS address is not public").into());
    }
    Ok(Box::new(addresses.into_iter()))
}
struct PublicResolver;
impl Resolve for PublicResolver {
    fn resolve(&self, name: Name) -> Resolving {
        Box::pin(async move {
            let addresses = tokio::time::timeout(
                Duration::from_secs(8),
                tokio::net::lookup_host((name.as_str(), 0)),
            )
            .await
            .map_err(|_| io::Error::new(io::ErrorKind::TimedOut, "source DNS timeout"))??;
            public_addresses(addresses.take(32).collect())
        })
    }
}
fn download_error(error: reqwest::Error) -> SourceError {
    if error.is_timeout() {
        SourceError::Timeout
    } else {
        SourceError::DownloadFailed
    }
}
#[derive(Default)]
struct FetchResponse {
    status: u16,
    location: Option<String>,
    content_type: Option<String>,
    content_length: Option<u64>,
    retry_after: Option<String>,
    body: Vec<u8>,
}
trait Transport: Send + Sync {
    fn get(
        &self,
        url: Url,
        addgene_token: Option<&str>,
    ) -> impl Future<Output = Result<FetchResponse, SourceError>> + Send;
}
struct NativeTransport {
    client: Client,
}
impl Transport for NativeTransport {
    async fn get(
        &self,
        url: Url,
        addgene_token: Option<&str>,
    ) -> Result<FetchResponse, SourceError> {
        validate_public_url(&url)?;
        if url.host_str() == Some("eutils.ncbi.nlm.nih.gov") {
            throttle_ncbi().await?;
        }
        if url.host_str() == Some(addgene::API_HOST) {
            // A conservative client policy, not a claim about the provider's contractual limit.
            throttle(&ADDGENE_RATE, Duration::from_secs(1)).await?;
        } else if addgene_token.is_some() {
            return Err(SourceError::InvalidConfig);
        }
        let mut request = self
            .client
            .get(url)
            .header("Accept-Encoding", "identity")
            .header(
                "Accept",
                "application/json, text/plain;q=0.9, application/octet-stream;q=0.9",
            );
        if let Some(token) = addgene_token {
            let mut value = reqwest::header::HeaderValue::from_str(&format!("Token {token}"))
                .map_err(|_| SourceError::InvalidConfig)?;
            value.set_sensitive(true);
            request = request.header(reqwest::header::AUTHORIZATION, value);
        }
        let mut response = request.send().await.map_err(download_error)?;
        let status = response.status().as_u16();
        let header = |name: &str| {
            response
                .headers()
                .get(name)
                .and_then(|v| v.to_str().ok())
                .map(str::to_owned)
        };
        let mut result = FetchResponse {
            status,
            location: header("location"),
            content_type: header("content-type"),
            content_length: response.content_length(),
            retry_after: header("retry-after"),
            body: vec![],
        };
        // Error/redirect bodies are never parsed or echoed, and need not be downloaded.
        if !(200..300).contains(&status) {
            if matches!(status, 429 | 503)
                && response.url().host_str() == Some("eutils.ncbi.nlm.nih.gov")
            {
                if let Some(seconds) = retry_after(result.retry_after.as_deref()) {
                    note_cooldown(seconds).await;
                }
            }
            if matches!(status, 429 | 503) && response.url().host_str() == Some(addgene::API_HOST) {
                if let Some(seconds) = retry_after(result.retry_after.as_deref()) {
                    record_cooldown(&ADDGENE_RATE, seconds).await;
                }
            }
            return Ok(result);
        }
        if result
            .content_length
            .is_some_and(|length| length > MAX_IMPORT_BYTES as u64)
        {
            return Err(SourceError::TooLarge);
        }
        if header("content-encoding").is_some_and(|v| !v.eq_ignore_ascii_case("identity")) {
            return Err(SourceError::UnexpectedContent);
        }
        while let Some(chunk) = response.chunk().await.map_err(download_error)? {
            if chunk.len() > MAX_IMPORT_BYTES.saturating_sub(result.body.len()) {
                return Err(SourceError::TooLarge);
            }
            result.body.extend_from_slice(&chunk);
        }
        if result
            .content_length
            .is_some_and(|length| length != result.body.len() as u64)
        {
            return Err(SourceError::Truncated);
        }
        Ok(result)
    }
}
fn retry_after(value: Option<&str>) -> Option<u64> {
    let value = value?;
    if let Ok(seconds) = value.parse::<u64>() {
        return Some(seconds);
    }
    let date = chrono::DateTime::parse_from_rfc2822(value)
        .ok()?
        .with_timezone(&chrono::Utc);
    Some((date - chrono::Utc::now()).num_seconds().max(0) as u64)
}
async fn retrieve<T: Transport>(
    transport: &T,
    source: &Source,
    config: &Config,
) -> Result<Vec<u8>, SourceError> {
    let url = source.download_url_with_config(config)?;
    retrieve_bytes(
        transport,
        url,
        config,
        matches!(
            source,
            Source::Ncbi { .. } | Source::Ena { .. } | Source::Ddbj { .. }
        ),
        false,
        false,
    )
    .await
}
async fn retrieve_bytes<T: Transport>(
    transport: &T,
    mut url: Url,
    config: &Config,
    fixed_origin: bool,
    addgene_auth: bool,
    json_response: bool,
) -> Result<Vec<u8>, SourceError> {
    let initial_origin = url.origin();
    for hop in 0..=MAX_REDIRECTS {
        validate_public_url(&url)?;
        let token = if addgene_auth && url.host_str() == Some(addgene::API_HOST) {
            config.addgene_api_token.as_deref()
        } else {
            None
        };
        let response = transport.get(url.clone(), token).await?;
        if matches!(response.status, 301 | 302 | 303 | 307 | 308) {
            if hop == MAX_REDIRECTS {
                return Err(SourceError::RedirectRejected);
            }
            let target = url
                .join(
                    response
                        .location
                        .as_deref()
                        .ok_or(SourceError::RedirectRejected)?,
                )
                .map_err(|_| SourceError::RedirectRejected)?;
            validate_public_url(&target).map_err(|_| SourceError::RedirectRejected)?;
            // API key/contact only belong to the fixed NCBI origin; never follow it elsewhere.
            if fixed_origin && target.origin() != initial_origin {
                return Err(SourceError::RedirectRejected);
            }
            url = target;
            continue;
        }
        let retry = retry_after(response.retry_after.as_deref());
        if response.status == 429 {
            return Err(SourceError::RateLimited {
                retry_after_seconds: retry,
            });
        }
        if !(200..300).contains(&response.status) {
            return Err(SourceError::HttpStatus {
                status: response.status,
                retry_after_seconds: retry,
            });
        }
        if response.body.len() > MAX_IMPORT_BYTES
            || response
                .content_length
                .is_some_and(|n| n > MAX_IMPORT_BYTES as u64)
        {
            return Err(SourceError::TooLarge);
        }
        if response
            .content_length
            .is_some_and(|n| n != response.body.len() as u64)
        {
            return Err(SourceError::Truncated);
        }
        if response.content_type.as_deref().is_some_and(|v| {
            v.split(';')
                .next()
                .unwrap_or_default()
                .trim()
                .eq_ignore_ascii_case("text/html")
        }) {
            return Err(SourceError::UnexpectedContent);
        }
        if !json_response {
            detect_format(&response.body)?;
        }
        return Ok(response.body);
    }
    Err(SourceError::RedirectRejected)
}

#[derive(Clone)]
pub struct SourceClient {
    transport: Arc<NativeTransport>,
    config: Config,
}
async fn discover_addgene<T: Transport>(
    transport: &T,
    plasmid_id: u64,
    config: &Config,
) -> Result<AddgeneCatalog, SourceError> {
    if config.addgene_api_token.is_none() {
        return Err(SourceError::AddgeneAccessRequired { plasmid_id });
    }
    let url =
        Url::parse(&addgene_catalog_url(plasmid_id)?).map_err(|_| SourceError::InvalidSource)?;
    let bytes = tokio::time::timeout(
        Duration::from_secs(60),
        retrieve_bytes(transport, url, config, true, true, true),
    )
    .await
    .map_err(|_| SourceError::Timeout)??;
    parse_addgene_catalog(plasmid_id, &bytes)
}
async fn retrieve_addgene_selection<T: Transport>(
    transport: &T,
    catalog: &AddgeneCatalog,
    sequence_id: u64,
    config: &Config,
) -> Result<ImportReport, SourceError> {
    let url =
        Url::parse(&catalog.download_url(sequence_id)?).map_err(|_| SourceError::AddgeneSchema)?;
    if url.host_str() == Some(addgene::API_HOST) && config.addgene_api_token.is_none() {
        return Err(SourceError::AddgeneAccessRequired {
            plasmid_id: catalog.plasmid_id,
        });
    }
    let bytes = tokio::time::timeout(
        Duration::from_secs(60),
        retrieve_bytes(transport, url, config, false, true, false),
    )
    .await
    .map_err(|_| SourceError::Timeout)??;
    let now = chrono::Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Nanos, true);
    catalog.import_download(sequence_id, &bytes, &now)
}
impl SourceClient {
    pub fn new(config: Config) -> Result<Self, SourceError> {
        let client = Client::builder()
            .redirect(Policy::none())
            .no_proxy()
            .dns_resolver(Arc::new(PublicResolver))
            .timeout(REQUEST_TIMEOUT)
            .connect_timeout(Duration::from_secs(10))
            .user_agent(concat!("Collagene2/", env!("CARGO_PKG_VERSION")))
            .build()
            .map_err(|_| SourceError::InvalidConfig)?;
        Ok(Self {
            transport: Arc::new(NativeTransport { client }),
            config,
        })
    }
    pub async fn import_url(&self, input: &str) -> Result<ImportReport, SourceError> {
        self.import_url_selected(input, None).await
    }
    pub async fn addgene_catalog(&self, input: &str) -> Result<AddgeneCatalog, SourceError> {
        let Source::Addgene { plasmid_id } = parse_source(input)? else {
            return Err(SourceError::InvalidSource);
        };
        discover_addgene(self.transport.as_ref(), plasmid_id, &self.config).await
    }
    /// Addgene always requires an explicit catalog sequence choice, including a single choice.
    /// The catalog is re-fetched so a selected ID cannot silently use stale/expired signed links.
    pub async fn import_url_selected(
        &self,
        input: &str,
        sequence_id: Option<u64>,
    ) -> Result<ImportReport, SourceError> {
        let source = parse_source(input)?;
        if matches!(source, Source::Addgene { .. }) {
            let catalog = self.addgene_catalog(input).await?;
            let id = sequence_id.ok_or_else(|| {
                if catalog.choices.is_empty() {
                    SourceError::AddgeneNoSequences {
                        plasmid_id: catalog.plasmid_id,
                    }
                } else {
                    SourceError::AddgeneSelectionRequired {
                        plasmid_id: catalog.plasmid_id,
                        sequences: catalog.choices.len(),
                    }
                }
            })?;
            return self.import_addgene_selection(&catalog, id).await;
        }
        if sequence_id.is_some() {
            return Err(SourceError::UnsupportedSelection);
        }
        // Also bounds the whole redirect chain and time spent awaiting the shared rate gate.
        let bytes = tokio::time::timeout(
            Duration::from_secs(60),
            retrieve(self.transport.as_ref(), &source, &self.config),
        )
        .await
        .map_err(|_| SourceError::Timeout)??;
        let now = chrono::Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Nanos, true);
        import_download(&source, &bytes, &now)
    }
    /// Download a selection from an already-discovered immutable catalog, without re-fetching it.
    /// An expired link fails explicitly; the caller can refresh and review choices again.
    pub async fn import_addgene_selection(
        &self,
        catalog: &AddgeneCatalog,
        sequence_id: u64,
    ) -> Result<ImportReport, SourceError> {
        retrieve_addgene_selection(self.transport.as_ref(), catalog, sequence_id, &self.config)
            .await
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;
    struct Mock {
        responses: Mutex<VecDeque<FetchResponse>>,
        requests: Mutex<Vec<Url>>,
        authorizations: Mutex<Vec<Option<String>>>,
    }
    impl Mock {
        fn new(responses: Vec<FetchResponse>) -> Self {
            Self {
                responses: Mutex::new(responses.into()),
                requests: Mutex::new(vec![]),
                authorizations: Mutex::new(vec![]),
            }
        }
    }
    impl Transport for Mock {
        async fn get(
            &self,
            url: Url,
            addgene_token: Option<&str>,
        ) -> Result<FetchResponse, SourceError> {
            self.requests.lock().await.push(url);
            self.authorizations
                .lock()
                .await
                .push(addgene_token.map(str::to_owned));
            self.responses
                .lock()
                .await
                .pop_front()
                .ok_or(SourceError::DownloadFailed)
        }
    }
    fn fasta() -> FetchResponse {
        FetchResponse {
            status: 200,
            content_type: Some("text/plain".into()),
            body: b">fixture\nACGT\n".to_vec(),
            ..Default::default()
        }
    }
    #[tokio::test]
    async fn redirects_are_validated_and_no_credentials_follow_ncbi_off_origin() {
        let source = parse_source("https://sequence.example.org/old.fa?signature=SECRET").unwrap();
        let mock = Mock::new(vec![
            FetchResponse {
                status: 302,
                location: Some("https://files.example.net/new.fa?signature=OTHER".into()),
                ..Default::default()
            },
            fasta(),
        ]);
        assert!(
            retrieve(&mock, &source, &Config::default())
                .await
                .unwrap()
                .starts_with(b">fixture")
        );
        let mock = Mock::new(vec![FetchResponse {
            status: 302,
            location: Some("http://127.0.0.1/internal".into()),
            ..Default::default()
        }]);
        assert_eq!(
            retrieve(&mock, &source, &Config::default())
                .await
                .unwrap_err(),
            SourceError::RedirectRejected
        );
        let ncbi = parse_source("NC_001422.1").unwrap();
        let config = Config::default().with_ncbi_api_key("SECRET").unwrap();
        let mock = Mock::new(vec![FetchResponse {
            status: 302,
            location: Some("https://files.example.net/steal".into()),
            ..Default::default()
        }]);
        assert_eq!(
            retrieve(&mock, &ncbi, &config).await.unwrap_err(),
            SourceError::RedirectRejected
        );
        assert_eq!(mock.requests.lock().await.len(), 1);
    }
    #[tokio::test]
    async fn expired_links_html_partial_responses_and_retry_hints_are_explicit() {
        let source = parse_source("https://sequence.example.org/file.gb?token=SECRET").unwrap();
        for (response, expected) in [
            (
                FetchResponse {
                    status: 403,
                    ..Default::default()
                },
                SourceError::HttpStatus {
                    status: 403,
                    retry_after_seconds: None,
                },
            ),
            (
                FetchResponse {
                    status: 429,
                    retry_after: Some("120".into()),
                    ..Default::default()
                },
                SourceError::RateLimited {
                    retry_after_seconds: Some(120),
                },
            ),
            (
                FetchResponse {
                    content_type: Some("text/html".into()),
                    ..fasta()
                },
                SourceError::UnexpectedContent,
            ),
            (
                FetchResponse {
                    content_length: Some(100),
                    ..fasta()
                },
                SourceError::Truncated,
            ),
            (
                FetchResponse {
                    content_length: Some(MAX_IMPORT_BYTES as u64 + 1),
                    ..fasta()
                },
                SourceError::TooLarge,
            ),
        ] {
            let mock = Mock::new(vec![response]);
            assert_eq!(
                retrieve(&mock, &source, &Config::default())
                    .await
                    .unwrap_err(),
                expected
            );
        }
        let mock = Mock::new(
            (0..6)
                .map(|_| FetchResponse {
                    status: 302,
                    location: Some("/again.fa".into()),
                    ..Default::default()
                })
                .collect(),
        );
        assert_eq!(
            retrieve(&mock, &source, &Config::default())
                .await
                .unwrap_err(),
            SourceError::RedirectRejected
        );
    }
    #[test]
    fn dns_answer_filter_rejects_mixed_public_private_and_transition_addresses() {
        assert!(public_addresses(vec!["8.8.8.8:443".parse().unwrap()]).is_ok());
        for blocked in [
            "127.0.0.1:443",
            "10.0.0.1:443",
            "169.254.169.254:443",
            "[::ffff:127.0.0.1]:443",
            "[2002:7f00:1::]:443",
        ] {
            assert!(
                public_addresses(vec![
                    "8.8.8.8:443".parse().unwrap(),
                    blocked.parse().unwrap()
                ])
                .is_err()
            );
        }
    }
    #[tokio::test]
    async fn ncbi_requests_share_a_conservative_rate_gate() {
        throttle_ncbi().await.unwrap();
        let start = Instant::now();
        throttle_ncbi().await.unwrap();
        assert!(start.elapsed() >= Duration::from_millis(330));
    }
    #[tokio::test]
    async fn ena_ddbj_transport_uses_annotated_routes_and_never_forwards_other_provider_credentials()
     {
        let config = Config::default()
            .with_ncbi_api_key("NCBI_SECRET")
            .unwrap()
            .with_contact_email("fixture@example.org")
            .unwrap()
            .with_addgene_api_token("ADDGENE_SECRET")
            .unwrap();
        for input in ["ena:AB123456.2", "ddbj:AB123456.2"] {
            let source = parse_source(input).unwrap();
            let mock = Mock::new(vec![fasta()]);
            retrieve(&mock, &source, &config).await.unwrap();
            let requests = mock.requests.lock().await;
            assert_eq!(requests[0].as_str(), source.download_url().unwrap());
            assert!(!requests[0].as_str().contains("SECRET"));
            assert!(!requests[0].as_str().contains("email"));
            assert_eq!(*mock.authorizations.lock().await, vec![None]);
            drop(requests);
            let mock = Mock::new(vec![FetchResponse {
                status: 302,
                location: Some("https://files.example.org/elsewhere.gb".into()),
                ..Default::default()
            }]);
            assert_eq!(
                retrieve(&mock, &source, &config).await.unwrap_err(),
                SourceError::RedirectRejected
            );
            assert_eq!(mock.requests.lock().await.len(), 1);
        }
    }
    fn catalog_response() -> FetchResponse {
        FetchResponse {
            status:200, content_type:Some("application/json".into()),
            body:serde_json::to_vec(&serde_json::json!({"id":12345,"name":"Fixture","sequences":{
                "public_user_full_sequences":[],"public_user_partial_sequences":[],"public_addgene_partial_sequences":[],
                "public_addgene_full_sequences":[{"sequence_id":42,"sequence_description":"Fixture full record","length":4,"sequence":"ACGT","genbank_api_url":"https://api.developers.addgene.org/download/genbank/42/","genbank_url":"https://api.developers.addgene.org/download/genbank/42/"}]
            }})).unwrap(),..Default::default()
        }
    }
    #[tokio::test]
    async fn official_catalog_and_selected_download_keep_token_at_provider_origin() {
        let config = Config::default()
            .with_addgene_api_token("FIXTURE_TOKEN")
            .unwrap();
        let mock = Mock::new(vec![catalog_response(),FetchResponse {
            status:302,location:Some("https://files.example.org/42.gb?signature=TRANSIENT_SECRET".into()),..Default::default()
        },FetchResponse {
            status:200,content_type:Some("text/plain".into()),body:b"LOCUS       fixture                    4 bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 acgt\n//\n".to_vec(),..Default::default()
        }]);
        let catalog = discover_addgene(&mock, 12345, &config).await.unwrap();
        assert_eq!(catalog.choices.len(), 1);
        let report = retrieve_addgene_selection(&mock, &catalog, 42, &config)
            .await
            .unwrap();
        assert_eq!(report.molecules[0].metadata["source:sequence_id"], "42");
        let requests = mock.requests.lock().await;
        assert_eq!(
            requests[0].as_str(),
            "https://api.developers.addgene.org/catalog/plasmid-with-sequences/12345/"
        );
        assert_eq!(
            requests[1].as_str(),
            "https://api.developers.addgene.org/download/genbank/42/"
        );
        assert_eq!(
            *mock.authorizations.lock().await,
            vec![
                Some("FIXTURE_TOKEN".into()),
                Some("FIXTURE_TOKEN".into()),
                None
            ]
        );
        assert!(
            !serde_json::to_string(&report.molecules)
                .unwrap()
                .contains("TRANSIENT_SECRET")
        );
    }
    #[tokio::test]
    async fn addgene_catalog_redirects_missing_access_expired_links_and_rate_hints_are_explicit() {
        let config = Config::default()
            .with_addgene_api_token("FIXTURE_TOKEN")
            .unwrap();
        let mock = Mock::new(vec![]);
        assert_eq!(
            discover_addgene(&mock, 12345, &Config::default())
                .await
                .unwrap_err(),
            SourceError::AddgeneAccessRequired { plasmid_id: 12345 }
        );
        assert!(mock.requests.lock().await.is_empty());
        let mock = Mock::new(vec![FetchResponse {
            status: 302,
            location: Some("https://files.example.org/leak".into()),
            ..Default::default()
        }]);
        assert_eq!(
            discover_addgene(&mock, 12345, &config).await.unwrap_err(),
            SourceError::RedirectRejected
        );
        assert_eq!(mock.requests.lock().await.len(), 1);
        for (status, expected) in [
            (
                401,
                SourceError::HttpStatus {
                    status: 401,
                    retry_after_seconds: Some(120),
                },
            ),
            (
                403,
                SourceError::HttpStatus {
                    status: 403,
                    retry_after_seconds: Some(120),
                },
            ),
            (
                429,
                SourceError::RateLimited {
                    retry_after_seconds: Some(120),
                },
            ),
        ] {
            let mock = Mock::new(vec![FetchResponse {
                status,
                retry_after: Some("120".into()),
                ..Default::default()
            }]);
            assert_eq!(
                discover_addgene(&mock, 12345, &config).await.unwrap_err(),
                expected
            );
        }
        let native = SourceClient::new(config).unwrap();
        assert!(matches!(
            native
                .transport
                .get(
                    Url::parse("https://files.example.org/file.gb").unwrap(),
                    Some("FIXTURE_TOKEN")
                )
                .await,
            Err(SourceError::InvalidConfig)
        ));
    }
}
