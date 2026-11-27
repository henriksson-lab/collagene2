//! Real native HTTP transport. Authentication material is injected, memory-only
//! and outside Debug/serde/repository metadata. No redirect, proxy or HTTP retry
//! can replay a credential-bearing mutation to another destination.
use crate::{Channel, RemoteConfig, RemoteTransport, TransportError};
use collagene_library_protocol as library;
use collagene_library_sync as sync;
use reqwest::{
    blocking::{Client, Request},
    header::{self, HeaderMap, HeaderValue},
};
use std::{
    collections::BTreeMap,
    io::Read,
    time::{Duration, Instant},
};

/// Database scope pin checked by the HTTP host before reading a request body.
pub const EXPECTED_LAB_HEADER: &str = "x-collagene-expected-lab";

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum HttpError {
    #[error("invalid HTTP transport configuration")]
    Configuration,
    #[error("remote must be a trusted HTTPS origin (or explicitly enabled literal loopback HTTP)")]
    Origin,
    #[error("invalid or mismatched binary request frame")]
    RequestFrame,
    #[error("credential provider unavailable")]
    Credential,
    #[error("HTTP request failed or was interrupted")]
    Network,
    #[error("HTTP request timed out")]
    Timeout,
    #[error("HTTP server returned status {0}")]
    Status(u16),
    #[error("unexpected HTTP response media type or content encoding")]
    MediaType,
    #[error("HTTP response exceeds the binary message bound")]
    TooLarge,
    #[error("HTTP content length does not match its body")]
    Length,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("credential provider unavailable")]
pub struct CredentialError;
/// Deliberately no Debug, Display or serde implementation. The underlying
/// sensitive HeaderValue is redacted in standard reqwest request diagnostics.
#[derive(Clone)]
pub struct BearerToken {
    header: HeaderValue,
}
impl BearerToken {
    pub fn new(token: &str) -> Result<Self, CredentialError> {
        if token.is_empty()
            || token.len() > 4096
            || !token.bytes().all(|b| (0x21..=0x7e).contains(&b))
        {
            return Err(CredentialError);
        }
        let mut bytes = Vec::with_capacity(7 + token.len());
        bytes.extend_from_slice(b"Bearer ");
        bytes.extend_from_slice(token.as_bytes());
        let mut header = HeaderValue::from_bytes(&bytes).map_err(|_| CredentialError)?;
        header.set_sensitive(true);
        Ok(Self { header })
    }
}
/// Implement this using an OS credential store or an existing authenticated host.
/// Its errors must not include raw credentials. It is never serialized by us.
pub trait CredentialProvider {
    fn bearer_for(&mut self, remote: &RemoteConfig) -> Result<BearerToken, CredentialError>;
}
/// Optional process-memory provider for a native session; not persisted anywhere.
#[derive(Default)]
pub struct MemoryCredentials {
    tokens: BTreeMap<(String, String, String), BearerToken>,
}
fn credential_scope(remote: &RemoteConfig) -> Result<(String, String, String), CredentialError> {
    remote.validate().map_err(|_| CredentialError)?;
    let url = url::Url::parse(&remote.endpoint).map_err(|_| CredentialError)?;
    if url.path() != "/" {
        return Err(CredentialError);
    }
    Ok((
        remote.id.clone(),
        url.origin().ascii_serialization(),
        remote.lab.as_str().into(),
    ))
}
impl MemoryCredentials {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn insert(
        &mut self,
        remote: &RemoteConfig,
        token: BearerToken,
    ) -> Result<(), CredentialError> {
        self.tokens.insert(credential_scope(remote)?, token);
        Ok(())
    }
    pub fn remove(&mut self, remote_id: &str) {
        self.tokens.retain(|(id, _, _), _| id != remote_id);
    }
    pub fn clear(&mut self) {
        self.tokens.clear();
    }
}
impl CredentialProvider for MemoryCredentials {
    fn bearer_for(&mut self, remote: &RemoteConfig) -> Result<BearerToken, CredentialError> {
        self.tokens
            .get(&credential_scope(remote)?)
            .cloned()
            .ok_or(CredentialError)
    }
}
#[derive(Debug, Clone, Copy)]
pub struct HttpOptions {
    /// Applies only to literal 127.0.0.0/8 and ::1. Hostname/DNS guesses do not
    /// qualify, including `localhost`; use an explicit loopback address instead.
    pub allow_loopback_http: bool,
    pub connect_timeout: Duration,
    /// Request/per-read timeout and the between-chunk overall budget. A read
    /// already blocked at budget expiry can consume one additional read timeout.
    pub total_timeout: Duration,
}
impl Default for HttpOptions {
    fn default() -> Self {
        Self {
            allow_loopback_http: false,
            connect_timeout: Duration::from_secs(5),
            total_timeout: Duration::from_secs(60),
        }
    }
}
impl HttpOptions {
    fn validate(self) -> Result<(), HttpError> {
        if self.connect_timeout < Duration::from_millis(10)
            || self.connect_timeout > Duration::from_secs(60)
            || self.total_timeout < self.connect_timeout
            || self.total_timeout > Duration::from_secs(300)
        {
            return Err(HttpError::Configuration);
        }
        Ok(())
    }
}
/// Holds no credential in its client default headers. Provider tokens are looked
/// up only after origin/frame validation and attached to that one fixed request.
pub struct HttpRemoteTransport<P> {
    provider: P,
    client: Client,
    options: HttpOptions,
    last_error: Option<HttpError>,
}
impl<P: CredentialProvider> HttpRemoteTransport<P> {
    pub fn new(provider: P, options: HttpOptions) -> Result<Self, HttpError> {
        options.validate()?;
        let client = Client::builder()
            .redirect(reqwest::redirect::Policy::none())
            .retry(reqwest::retry::never())
            .no_proxy()
            .use_rustls_tls()
            .no_gzip()
            .no_brotli()
            .no_deflate()
            .no_zstd()
            .https_only(!options.allow_loopback_http)
            .connect_timeout(options.connect_timeout)
            .timeout(options.total_timeout)
            .build()
            .map_err(|_| HttpError::Configuration)?;
        Ok(Self {
            provider,
            client,
            options,
            last_error: None,
        })
    }
    /// Returns only a static classified failure, never a URL, header or body.
    pub fn last_error(&self) -> Option<HttpError> {
        self.last_error
    }
    /// For rotating/removing credentials in the trusted native session. Nothing
    /// here is copied into RemoteConfig, checkpoints or biological archives.
    pub fn credentials_mut(&mut self) -> &mut P {
        &mut self.provider
    }
    pub fn exchange_checked(
        &mut self,
        remote: &RemoteConfig,
        channel: Channel,
        bytes: &[u8],
    ) -> Result<Vec<u8>, HttpError> {
        let client = self.client.clone();
        let result = exchange_with(
            &mut self.provider,
            &self.client,
            self.options,
            remote,
            channel,
            bytes,
            move |request| {
                let response = client.execute(request).map_err(|e| {
                    if e.is_timeout() {
                        HttpError::Timeout
                    } else {
                        HttpError::Network
                    }
                })?;
                Ok(ResponseBody {
                    status: response.status().as_u16(),
                    headers: response.headers().clone(),
                    body: response,
                })
            },
        );
        self.last_error = result.as_ref().err().copied();
        result
    }
}
impl<P: CredentialProvider> RemoteTransport for HttpRemoteTransport<P> {
    fn exchange(
        &mut self,
        remote: &RemoteConfig,
        channel: Channel,
        request: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        self.exchange_checked(remote, channel, request)
            .map_err(|_| TransportError::Unavailable)
    }
}
fn media(channel: Channel) -> &'static str {
    match channel {
        Channel::Library => library::MEDIA_TYPE,
        Channel::Sync => sync::MEDIA_TYPE,
        Channel::Discussion => library::discussion::MEDIA_TYPE,
    }
}
fn response_limit(channel: Channel) -> usize {
    match channel {
        Channel::Library => library::MAX_MESSAGE_BYTES + library::HEADER_BYTES,
        Channel::Sync => sync::MAX_RESPONSE_BYTES + sync::HEADER_BYTES,
        Channel::Discussion => library::discussion::MAX_RESPONSE_BYTES,
    }
}
fn destination(
    remote: &RemoteConfig,
    channel: Channel,
    options: HttpOptions,
) -> Result<url::Url, HttpError> {
    remote.validate().map_err(|_| HttpError::Origin)?;
    let mut url = url::Url::parse(&remote.endpoint).map_err(|_| HttpError::Origin)?;
    if url.path() != "/" {
        return Err(HttpError::Origin);
    }
    if url.scheme() != "https" {
        let loopback = match url.host() {
            Some(url::Host::Ipv4(ip)) => ip.is_loopback(),
            Some(url::Host::Ipv6(ip)) => ip.is_loopback(),
            _ => false,
        };
        if url.scheme() != "http" || !options.allow_loopback_http || !loopback {
            return Err(HttpError::Origin);
        }
    }
    url.set_path(match channel {
        Channel::Library => "/api/library",
        Channel::Sync => sync::ENDPOINT,
        Channel::Discussion => library::discussion::ENDPOINT,
    });
    Ok(url)
}
/// Validate the borrowed frame without duplicating/deserializing a 64MiB native
/// archive. The server still performs typed decoding; repository validates each
/// typed response's correlation/context and full mutation receipt separately.
fn request_frame(
    channel: Channel,
    bytes: &[u8],
    expected_lab: &library::LabId,
) -> Result<(), HttpError> {
    if channel == Channel::Discussion {
        return library::discussion::decode_request(bytes)
            .map(|_| ())
            .map_err(|_| HttpError::RequestFrame);
    }
    if channel == Channel::Sync {
        let prefix = sync::request_prefix(bytes).map_err(|_| HttpError::RequestFrame)?;
        return if prefix.expected_lab == expected_lab.as_str() {
            Ok(())
        } else {
            Err(HttpError::RequestFrame)
        };
    }
    let (magic, version, cap) = (
        b"CLQ\0",
        library::SCHEMA_VERSION,
        library::MAX_MESSAGE_BYTES,
    );
    if bytes.len() < 42
        || bytes.get(..4) != Some(magic.as_slice())
        || u16::from_le_bytes([bytes[4], bytes[5]]) != version
    {
        return Err(HttpError::RequestFrame);
    }
    let size = u32::from_le_bytes(
        bytes[6..10]
            .try_into()
            .map_err(|_| HttpError::RequestFrame)?,
    ) as usize;
    if size == 0 || size > cap || size.checked_add(42) != Some(bytes.len()) {
        return Err(HttpError::RequestFrame);
    }
    let mut declared = String::with_capacity(64);
    for b in &bytes[10..42] {
        use std::fmt::Write;
        let _ = write!(&mut declared, "{b:02x}");
    }
    if library::Digest::of(&bytes[42..]).as_str() != declared {
        return Err(HttpError::RequestFrame);
    }
    Ok(())
}
struct ResponseBody<R> {
    status: u16,
    headers: HeaderMap,
    body: R,
}
fn exchange_with<P, F, R>(
    provider: &mut P,
    client: &Client,
    options: HttpOptions,
    remote: &RemoteConfig,
    channel: Channel,
    bytes: &[u8],
    execute: F,
) -> Result<Vec<u8>, HttpError>
where
    P: CredentialProvider,
    F: FnOnce(Request) -> Result<ResponseBody<R>, HttpError>,
    R: Read,
{
    options.validate()?;
    let url = destination(remote, channel, options)?;
    request_frame(channel, bytes, &remote.lab)?;
    let credential = provider
        .bearer_for(remote)
        .map_err(|_| HttpError::Credential)?;
    let deadline = Instant::now()
        .checked_add(options.total_timeout)
        .ok_or(HttpError::Configuration)?;
    let request = client
        .post(url)
        .header(header::CONTENT_TYPE, media(channel))
        .header(header::ACCEPT, media(channel))
        .header(header::ACCEPT_ENCODING, "identity")
        .header(EXPECTED_LAB_HEADER, remote.lab.as_str())
        .header(header::AUTHORIZATION, credential.header)
        .body(bytes.to_vec())
        .build()
        .map_err(|_| HttpError::Configuration)?;
    let mut response = execute(request)?;
    if response.status != 200 {
        return Err(HttpError::Status(response.status));
    }
    if !single_header(&response.headers, header::CONTENT_TYPE)?
        .is_some_and(|s| s.trim().eq_ignore_ascii_case(media(channel)))
        || single_header(&response.headers, header::CONTENT_ENCODING)?
            .is_some_and(|s| !s.trim().eq_ignore_ascii_case("identity"))
    {
        return Err(HttpError::MediaType);
    }
    let length = single_header(&response.headers, header::CONTENT_LENGTH)?
        .map(|s| s.parse::<u64>().map_err(|_| HttpError::Length))
        .transpose()?;
    bounded_body(
        &mut response.body,
        length,
        response_limit(channel),
        deadline,
    )
}
fn single_header(headers: &HeaderMap, name: header::HeaderName) -> Result<Option<&str>, HttpError> {
    let mut values = headers.get_all(name).iter();
    let first = values.next();
    if values.next().is_some() {
        return Err(HttpError::MediaType);
    }
    first
        .map(|v| v.to_str().map_err(|_| HttpError::MediaType))
        .transpose()
}
/// A missing/lying Content-Length cannot turn into an unbounded read. At most
/// cap+1 bytes are requested; overflow bytes are never appended to the result.
fn bounded_body<R: Read>(
    body: &mut R,
    length: Option<u64>,
    cap: usize,
    deadline: Instant,
) -> Result<Vec<u8>, HttpError> {
    if length.is_some_and(|n| n > cap as u64) {
        return Err(HttpError::TooLarge);
    }
    let mut result = Vec::with_capacity(length.unwrap_or(8192).min(65536).min(cap as u64) as usize);
    let mut buffer = [0u8; 16384];
    loop {
        if Instant::now() >= deadline {
            return Err(HttpError::Timeout);
        }
        let take = (cap - result.len()).saturating_add(1).min(buffer.len());
        let n = body.read(&mut buffer[..take]).map_err(|e| {
            if e.kind() == std::io::ErrorKind::TimedOut {
                HttpError::Timeout
            } else {
                HttpError::Network
            }
        })?;
        if Instant::now() >= deadline {
            return Err(HttpError::Timeout);
        }
        if n == 0 {
            break;
        }
        if n > cap - result.len() {
            return Err(HttpError::TooLarge);
        }
        result.extend_from_slice(&buffer[..n]);
    }
    if length.is_some_and(|n| n != result.len() as u64) {
        return Err(HttpError::Length);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        cell::Cell,
        io::{self, Cursor},
    };
    fn remote(endpoint: &str) -> RemoteConfig {
        RemoteConfig {
            id: "rem-0123456789abcdef0123456789abcdef".into(),
            name: "Configured lab".into(),
            endpoint: endpoint.into(),
            lab: library::LabId::new("lab-0123456789abcdef0123456789abcdef").unwrap(),
        }
    }
    fn frame(channel: Channel) -> Vec<u8> {
        match channel {
            Channel::Discussion => {
                library::discussion::encode_request(&library::discussion::RequestEnvelope {
                    request_id: 19,
                    request: library::discussion::Request::Identify,
                })
                .unwrap()
            }
            Channel::Library => library::encode_request(&library::RequestEnvelope {
                request_id: 17,
                request: library::Request::Browse {
                    parent: None,
                    page: Default::default(),
                },
            })
            .unwrap(),
            Channel::Sync => sync::encode_request(&sync::RequestEnvelope {
                request_id: 18,
                expected_lab: library::LabId::new("lab-0123456789abcdef0123456789abcdef").unwrap(),
                request: sync::Request::Lookup {
                    mutation_id: sync::MutationId::new("mut-0123456789abcdef0123456789abcdef")
                        .unwrap(),
                    expected_intent: library::Digest::of(b"queued intent"),
                },
            })
            .unwrap(),
        }
    }
    fn headers(channel: Channel) -> HeaderMap {
        let mut h = HeaderMap::new();
        h.insert(
            header::CONTENT_TYPE,
            HeaderValue::from_static(media(channel)),
        );
        h
    }
    fn transport() -> HttpRemoteTransport<MemoryCredentials> {
        HttpRemoteTransport::new(MemoryCredentials::new(), HttpOptions::default()).unwrap()
    }
    #[test]
    fn origin_policy_does_not_turn_names_paths_or_external_http_into_loopback_authority() {
        let secure = HttpOptions::default();
        let local = HttpOptions {
            allow_loopback_http: true,
            ..secure
        };
        for endpoint in [
            "https://lab.example/path",
            "https://name:secret@lab.example",
            "https://lab.example?token=secret",
            "https://lab.example#secret",
            "file:///tmp/plasmids",
            "http://lab.example",
            "http://localhost:8085",
            "http://127.0.0.1.evil.example",
        ] {
            assert!(
                destination(&remote(endpoint), Channel::Library, local).is_err(),
                "{endpoint}"
            );
        }
        assert!(destination(&remote("http://127.0.0.1:8085"), Channel::Library, secure).is_err());
        assert_eq!(
            destination(&remote("http://127.0.0.1:8085"), Channel::Library, local)
                .unwrap()
                .as_str(),
            "http://127.0.0.1:8085/api/library"
        );
        assert_eq!(
            destination(&remote("http://[::1]:8085"), Channel::Sync, local)
                .unwrap()
                .as_str(),
            "http://[::1]:8085/api/library/sync"
        );
        assert_eq!(
            destination(&remote("https://lab.example"), Channel::Sync, secure)
                .unwrap()
                .as_str(),
            "https://lab.example/api/library/sync"
        );
        assert!(
            HttpRemoteTransport::new(
                MemoryCredentials::new(),
                HttpOptions {
                    connect_timeout: Duration::ZERO,
                    ..secure
                }
            )
            .is_err()
        );
        assert!(
            HttpRemoteTransport::new(
                MemoryCredentials::new(),
                HttpOptions {
                    total_timeout: Duration::from_secs(301),
                    ..secure
                }
            )
            .is_err()
        );
    }
    #[test]
    fn fixed_requests_attach_sensitive_token_only_to_chosen_origin_and_channel() {
        let transport = transport();
        let origin = remote("https://lab.example:9443");
        let token = "offline-fixture-secret";
        let mut credentials = MemoryCredentials::new();
        credentials
            .insert(&origin, BearerToken::new(token).unwrap())
            .unwrap();
        for channel in [Channel::Library, Channel::Sync, Channel::Discussion] {
            let bytes = frame(channel);
            let called = Cell::new(0);
            let reply = exchange_with(
                &mut credentials,
                &transport.client,
                transport.options,
                &origin,
                channel,
                &bytes,
                |request| {
                    called.set(called.get() + 1);
                    assert_eq!(request.method(), reqwest::Method::POST);
                    assert_eq!(
                        request.url().origin(),
                        url::Url::parse(&origin.endpoint).unwrap().origin()
                    );
                    assert_eq!(
                        request.url().path(),
                        match channel {
                            Channel::Library => "/api/library",
                            Channel::Sync => sync::ENDPOINT,
                            Channel::Discussion => library::discussion::ENDPOINT,
                        }
                    );
                    assert_eq!(request.headers()[header::CONTENT_TYPE], media(channel));
                    assert_eq!(request.headers()[header::ACCEPT], media(channel));
                    assert_eq!(request.headers()[header::ACCEPT_ENCODING], "identity");
                    assert_eq!(request.headers()[EXPECTED_LAB_HEADER], origin.lab.as_str());
                    let auth = &request.headers()[header::AUTHORIZATION];
                    assert!(auth.is_sensitive());
                    assert_eq!(auth.as_bytes(), format!("Bearer {token}").as_bytes());
                    assert!(!format!("{request:?}").contains(token));
                    Ok(ResponseBody {
                        status: 200,
                        headers: headers(channel),
                        body: Cursor::new(b"bounded opaque reply".to_vec()),
                    })
                },
            )
            .unwrap();
            assert_eq!(called.get(), 1);
            assert_eq!(reply, b"bounded opaque reply");
        }
        let mut changed = origin.clone();
        changed.endpoint = "https://other.example".into();
        assert!(credentials.bearer_for(&changed).is_err());
        let mut changed = origin.clone();
        changed.lab = library::LabId::new("lab-ffffffffffffffffffffffffffffffff").unwrap();
        assert!(credentials.bearer_for(&changed).is_err());
        let mut canonical = origin.clone();
        canonical.endpoint.push('/');
        assert!(credentials.bearer_for(&canonical).is_ok());
        credentials.remove(&origin.id);
        assert!(credentials.bearer_for(&origin).is_err());
        for token in ["", "header\r\ninjection", "has whitespace", "nonascii-é"] {
            assert!(BearerToken::new(token).is_err());
        }
        assert!(BearerToken::new(&"a".repeat(4097)).is_err());
    }
    struct SpyProvider<'a>(&'a Cell<usize>);
    impl CredentialProvider for SpyProvider<'_> {
        fn bearer_for(&mut self, _: &RemoteConfig) -> Result<BearerToken, CredentialError> {
            self.0.set(self.0.get() + 1);
            BearerToken::new("fixture")
        }
    }
    #[test]
    fn invalid_frames_and_origins_do_not_request_credentials_or_invoke_executor() {
        let transport = transport();
        let calls = Cell::new(0);
        let mut provider = SpyProvider(&calls);
        let origin = remote("https://lab.example");
        let good = frame(Channel::Library);
        let mut corrupt = good.clone();
        *corrupt.last_mut().unwrap() ^= 1;
        let mut future = good.clone();
        future[4] = 99;
        let mut trailing = good.clone();
        trailing.push(0);
        let mut old_discussion = frame(Channel::Discussion);
        old_discussion[4..6].copy_from_slice(&1u16.to_le_bytes());
        for (channel, bytes) in [
            (Channel::Discussion, good.clone()),
            (Channel::Discussion, old_discussion),
            (
                Channel::Discussion,
                vec![0; library::discussion::MAX_REQUEST_BYTES + 1],
            ),
            (Channel::Sync, good.clone()),
            (Channel::Library, corrupt),
            (Channel::Library, future),
            (Channel::Library, trailing),
            (Channel::Library, vec![]),
        ] {
            let result = exchange_with(
                &mut provider,
                &transport.client,
                transport.options,
                &origin,
                channel,
                &bytes,
                |_| -> Result<ResponseBody<Cursor<Vec<u8>>>, HttpError> {
                    panic!("invalid request reached executor")
                },
            );
            assert_eq!(result.unwrap_err(), HttpError::RequestFrame);
        }
        let mut wrong_pin = sync::decode_request(&frame(Channel::Sync)).unwrap();
        wrong_pin.expected_lab =
            library::LabId::new("lab-ffffffffffffffffffffffffffffffff").unwrap();
        let wrong_pin = sync::encode_request(&wrong_pin).unwrap();
        let result = exchange_with(
            &mut provider,
            &transport.client,
            transport.options,
            &origin,
            Channel::Sync,
            &wrong_pin,
            |_| -> Result<ResponseBody<Cursor<Vec<u8>>>, HttpError> {
                panic!("mismatched database pin reached executor")
            },
        );
        assert_eq!(result.unwrap_err(), HttpError::RequestFrame);
        let result = exchange_with(
            &mut provider,
            &transport.client,
            transport.options,
            &remote("http://external.example"),
            Channel::Library,
            &good,
            |_| -> Result<ResponseBody<Cursor<Vec<u8>>>, HttpError> {
                panic!("untrusted origin reached executor")
            },
        );
        assert_eq!(result.unwrap_err(), HttpError::Origin);
        assert_eq!(calls.get(), 0);
    }
    struct MustNotRead;
    impl Read for MustNotRead {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            panic!("rejected HTTP response body was read")
        }
    }
    #[test]
    fn redirects_http_errors_wrong_media_encoding_and_advertised_overflow_stop_before_body() {
        let transport = transport();
        let origin = remote("https://lab.example");
        let mut provider = MemoryCredentials::new();
        provider
            .insert(&origin, BearerToken::new("fixture").unwrap())
            .unwrap();
        let bytes = frame(Channel::Library);
        for status in [301, 302, 307, 308, 401, 403, 429, 500] {
            let called = Cell::new(0);
            let result = exchange_with(
                &mut provider,
                &transport.client,
                transport.options,
                &origin,
                Channel::Library,
                &bytes,
                |_| {
                    called.set(called.get() + 1);
                    let mut h = headers(Channel::Library);
                    h.insert(
                        header::LOCATION,
                        HeaderValue::from_static("https://different.example/steal"),
                    );
                    Ok(ResponseBody {
                        status,
                        headers: h,
                        body: MustNotRead,
                    })
                },
            );
            assert_eq!(result.unwrap_err(), HttpError::Status(status));
            assert_eq!(called.get(), 1);
        }
        for h in [
            HeaderMap::new(),
            {
                let mut h = headers(Channel::Library);
                h.insert(header::CONTENT_TYPE, HeaderValue::from_static("text/html"));
                h
            },
            {
                let mut h = headers(Channel::Library);
                h.insert(header::CONTENT_ENCODING, HeaderValue::from_static("gzip"));
                h
            },
            {
                let mut h = headers(Channel::Library);
                h.append(
                    header::CONTENT_TYPE,
                    HeaderValue::from_static(library::MEDIA_TYPE),
                );
                h
            },
        ] {
            let result = exchange_with(
                &mut provider,
                &transport.client,
                transport.options,
                &origin,
                Channel::Library,
                &bytes,
                |_| {
                    Ok(ResponseBody {
                        status: 200,
                        headers: h,
                        body: MustNotRead,
                    })
                },
            );
            assert_eq!(result.unwrap_err(), HttpError::MediaType);
        }
        let mut h = headers(Channel::Library);
        h.insert(
            header::CONTENT_LENGTH,
            HeaderValue::from_str(&(response_limit(Channel::Library) + 1).to_string()).unwrap(),
        );
        let result = exchange_with(
            &mut provider,
            &transport.client,
            transport.options,
            &origin,
            Channel::Library,
            &bytes,
            |_| {
                Ok(ResponseBody {
                    status: 200,
                    headers: h,
                    body: MustNotRead,
                })
            },
        );
        assert_eq!(result.unwrap_err(), HttpError::TooLarge);
        assert!(response_limit(Channel::Sync) < response_limit(Channel::Library));
    }
    struct Counting<'a> {
        left: usize,
        read: &'a Cell<usize>,
    }
    impl Read for Counting<'_> {
        fn read(&mut self, b: &mut [u8]) -> io::Result<usize> {
            let n = b.len().min(self.left);
            b[..n].fill(0);
            self.left -= n;
            self.read.set(self.read.get() + n);
            Ok(n)
        }
    }
    #[test]
    fn bounded_reader_checks_missing_lying_and_truncated_lengths_with_one_overflow_probe() {
        let deadline = || Instant::now() + Duration::from_secs(3600);
        let count = Cell::new(0);
        let mut huge = Counting {
            left: usize::MAX,
            read: &count,
        };
        assert_eq!(
            bounded_body(&mut huge, None, 32, deadline()).unwrap_err(),
            HttpError::TooLarge
        );
        assert_eq!(count.get(), 33);
        let count = Cell::new(0);
        let mut huge = Counting {
            left: usize::MAX,
            read: &count,
        };
        assert_eq!(
            bounded_body(&mut huge, Some(999), 32, deadline()).unwrap_err(),
            HttpError::TooLarge
        );
        assert_eq!(count.get(), 0);
        assert_eq!(
            bounded_body(&mut Cursor::new(vec![0; 10]), Some(9), 32, deadline()).unwrap_err(),
            HttpError::Length
        );
        assert_eq!(
            bounded_body(&mut Cursor::new(vec![0; 8]), Some(9), 32, deadline()).unwrap_err(),
            HttpError::Length
        );
        assert_eq!(
            bounded_body(&mut Cursor::new(vec![1; 32]), Some(32), 32, deadline()).unwrap(),
            vec![1; 32]
        );
        assert_eq!(
            bounded_body(&mut Cursor::new(vec![1; 8]), None, 32, deadline()).unwrap(),
            vec![1; 8]
        );
        assert_eq!(
            bounded_body(&mut MustNotRead, None, 32, Instant::now()).unwrap_err(),
            HttpError::Timeout
        );
        struct Broken;
        impl Read for Broken {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::other(
                    "untrusted backend message must not escape",
                ))
            }
        }
        assert_eq!(
            bounded_body(&mut Broken, None, 32, deadline()).unwrap_err(),
            HttpError::Network
        );
    }
}
