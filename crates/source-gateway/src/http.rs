//! Same-origin route and optional HTTP/1 serving adapter. This module serves no
//! static files, reads no filesystem paths and launches no external process.
use crate::{
    service::{Backend, Cancellation, Config, Gateway, origin_parts},
    wire::{self, ErrorCode, Failure},
};
use ::http::{HeaderMap, HeaderValue, Method, Request, Response, StatusCode, header};
use bytes::Bytes;
use http_body_util::{BodyExt, Full};
use hyper::{body::Incoming, server::conn::http1, service::service_fn};
use hyper_util::rt::{TokioIo, TokioTimer};
use std::{convert::Infallible, io, sync::Arc, time::Duration};
use tokio::{net::TcpListener, sync::Semaphore, task::JoinSet};
pub const PATH: &str = "/api/source-import/v1";
pub const CONTENT_TYPE: &str = "application/vnd.collagene.source+postcard";

fn forbidden() -> Failure {
    Failure::new(
        ErrorCode::Forbidden,
        "This source-import request is not permitted by the serving host.",
    )
}
fn authorize(
    config: &Config,
    method: &Method,
    uri: &::http::Uri,
    headers: &HeaderMap,
) -> Result<(), Failure> {
    if method != Method::POST
        || uri.path() != PATH
        || uri.query().is_some()
        || uri.authority().is_some()
        || uri.scheme().is_some()
    {
        return Err(Failure::new(
            ErrorCode::InvalidRequest,
            "Unknown source-import route or method.",
        ));
    }
    if headers.len() > 32
        || headers
            .iter()
            .map(|(name, value)| name.as_str().len() + value.as_bytes().len())
            .sum::<usize>()
            > 8192
    {
        return Err(forbidden());
    }
    for key in [
        header::HOST,
        header::ORIGIN,
        header::AUTHORIZATION,
        header::CONTENT_TYPE,
        header::CONTENT_LENGTH,
        header::CONTENT_ENCODING,
    ] {
        if headers.get_all(&key).iter().count() > 1 {
            return Err(forbidden());
        }
    }
    let host = headers
        .get(header::HOST)
        .and_then(|h| h.to_str().ok())
        .ok_or_else(forbidden)?
        .to_ascii_lowercase();
    if !config.authorities.contains(&host) {
        return Err(forbidden());
    }
    if let Some(origin) = headers.get(header::ORIGIN) {
        let value = origin.to_str().map_err(|_| forbidden())?;
        let (canonical, authority) = origin_parts(value).ok_or_else(forbidden)?;
        if canonical != value || !config.origins.contains(&canonical) || authority != host {
            return Err(forbidden());
        }
    }
    if headers.get_all("sec-fetch-site").iter().count() > 1
        || headers
            .get("sec-fetch-site")
            .is_some_and(|v| !matches!(v.to_str(), Ok("same-origin" | "none")))
    {
        return Err(forbidden());
    }
    if headers
        .get(header::CONTENT_TYPE)
        .and_then(|h| h.to_str().ok())
        != Some(CONTENT_TYPE)
        || headers
            .get(header::CONTENT_ENCODING)
            .is_some_and(|h| h.as_bytes() != b"identity")
    {
        return Err(Failure::new(
            ErrorCode::InvalidRequest,
            "Use the source-import binary content type without content encoding.",
        ));
    }
    if let Some(length) = headers.get(header::CONTENT_LENGTH) {
        let length = length
            .to_str()
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .ok_or_else(forbidden)?;
        if length > wire::MAX_REQUEST {
            return Err(Failure::new(
                ErrorCode::TooLarge,
                "Source-import request body exceeds16KiB.",
            ));
        }
    }
    match (&config.access_digest, headers.get(header::AUTHORIZATION)) {
        (None, None) => {}
        (Some(expected), Some(value)) => {
            let value = value.to_str().map_err(|_| forbidden())?;
            let (scheme, token) = value.split_once(' ').ok_or_else(forbidden)?;
            if !scheme.eq_ignore_ascii_case("Bearer")
                || !(32..=256).contains(&token.len())
                || !token.bytes().all(|b| b.is_ascii_graphic())
            {
                return Err(forbidden());
            }
            let actual = blake3::hash(token.as_bytes());
            let mut difference = 0u8;
            for (a, b) in actual.as_bytes().iter().zip(expected) {
                difference |= *a ^ *b;
            }
            if difference != 0 {
                return Err(forbidden());
            }
        }
        _ => return Err(forbidden()),
    }
    Ok(())
}
fn status(failure: &Failure) -> StatusCode {
    match failure.code {
        ErrorCode::Forbidden => StatusCode::FORBIDDEN,
        ErrorCode::Busy => StatusCode::SERVICE_UNAVAILABLE,
        ErrorCode::Timeout => StatusCode::GATEWAY_TIMEOUT,
        ErrorCode::Canceled => StatusCode::REQUEST_TIMEOUT,
        ErrorCode::ExpiredCatalog => StatusCode::GONE,
        ErrorCode::Provider => StatusCode::UNPROCESSABLE_ENTITY,
        ErrorCode::TooLarge => StatusCode::PAYLOAD_TOO_LARGE,
        ErrorCode::Internal => StatusCode::INTERNAL_SERVER_ERROR,
        ErrorCode::InvalidRequest => StatusCode::BAD_REQUEST,
    }
}
fn response(packet: wire::Response) -> Response<Full<Bytes>> {
    let mut status_code = packet.result.as_ref().err().map_or(StatusCode::OK, status);
    let body = match wire::encode_response(&packet) {
        Ok(bytes) => bytes,
        Err(_) => {
            status_code = StatusCode::PAYLOAD_TOO_LARGE;
            wire::encode_response(&wire::Response {
                id: packet.id,
                result: Err(Failure::new(
                    ErrorCode::TooLarge,
                    "Source-import reply exceeds64MiB; no partial document was returned.",
                )),
            })
            .unwrap_or_default()
        }
    };
    let mut result = Response::new(Full::new(Bytes::from(body)));
    *result.status_mut() = status_code;
    let headers = result.headers_mut();
    headers.insert(header::CONTENT_TYPE, HeaderValue::from_static(CONTENT_TYPE));
    headers.insert(header::CACHE_CONTROL, HeaderValue::from_static("no-store"));
    headers.insert(
        "x-content-type-options",
        HeaderValue::from_static("nosniff"),
    );
    headers.insert(
        "cross-origin-resource-policy",
        HeaderValue::from_static("same-origin"),
    );
    headers.insert("referrer-policy", HeaderValue::from_static("no-referrer"));
    // No Access-Control-Allow-Origin or cookie is ever emitted.
    result
}
impl<B: Backend> Gateway<B> {
    /// Reusable complete-body router for native hosts and in-process fixtures.
    /// The network adapter below also enforces bounds while reading, before this
    /// function is entered. Host authorization is never delegated to CORS alone.
    pub async fn handle_http(
        &self,
        request: Request<Bytes>,
        cancel: Cancellation,
    ) -> Response<Full<Bytes>> {
        if let Err(error) = authorize(
            self.config(),
            request.method(),
            request.uri(),
            request.headers(),
        ) {
            return response(wire::Response {
                id: 0,
                result: Err(error),
            });
        }
        let (_, body) = request.into_parts();
        let decoded = wire::decode_request(&body);
        drop(body);
        let packet = match decoded {
            Ok(request) => self.execute(request, cancel).await,
            Err(error) => wire::Response {
                id: 0,
                result: Err(Failure::new(
                    match error {
                        wire::WireError::TooLarge => ErrorCode::TooLarge,
                        _ => ErrorCode::InvalidRequest,
                    },
                    "Invalid bounded source-import message.",
                )),
            },
        };
        response(packet)
    }
}
async fn read_request(request: Request<Incoming>) -> Result<Request<Bytes>, Failure> {
    let (parts, mut incoming) = request.into_parts();
    let expected = parts
        .headers
        .get(header::CONTENT_LENGTH)
        .and_then(|h| h.to_str().ok())
        .and_then(|h| h.parse::<usize>().ok());
    let mut body = Vec::new();
    while let Some(frame) = incoming.frame().await {
        let frame = frame.map_err(|_| {
            Failure::new(
                ErrorCode::InvalidRequest,
                "The source-import request body was incomplete.",
            )
        })?;
        if let Ok(bytes) = frame.into_data() {
            if bytes.len() > wire::MAX_REQUEST.saturating_sub(body.len()) {
                return Err(Failure::new(
                    ErrorCode::TooLarge,
                    "Source-import request body exceeds16KiB.",
                ));
            }
            body.try_reserve(bytes.len()).map_err(|_| {
                Failure::new(
                    ErrorCode::TooLarge,
                    "The source-import body could not be buffered.",
                )
            })?;
            body.extend_from_slice(&bytes);
        }
    }
    if expected.is_some_and(|length| length != body.len()) {
        return Err(Failure::new(
            ErrorCode::InvalidRequest,
            "The source-import request length disagrees with its body.",
        ));
    }
    Ok(Request::from_parts(parts, Bytes::from(body)))
}
/// Optional standalone native HTTP/1 adapter. The embedding host may instead
/// mount handle_http at PATH beside its immutable frontend. Bind is loopback by
/// default; remote serving requires Config::remote and explicit HTTPS origins and
/// a capability. Terminate TLS in the owning host/proxy; do not expose this plain
/// HTTP listener as a credential-bearing public frontend.
///
/// Shutdown and dropped connections abort native I/O jobs. The native source
/// parser's bounded synchronous phase cannot be forcibly preempted; its permit
/// remains owned until that phase finishes rather than leaking unlimited work.
pub async fn serve<B: Backend>(gateway: Gateway<B>, shutdown: Cancellation) -> io::Result<()> {
    let listener = TcpListener::bind(gateway.config().bind_address()).await?;
    let connections = Arc::new(Semaphore::new(gateway.config().limits.connections));
    let mut tasks = JoinSet::new();
    loop {
        while tasks.try_join_next().is_some() {}
        tokio::select! {
            _=shutdown.cancelled()=>break,
            accepted=listener.accept()=>{
                let (stream,_peer)=accepted?;
                let Ok(permit)=connections.clone().try_acquire_owned() else {drop(stream);continue;};
                let gateway=gateway.clone();
                tasks.spawn(async move {
                    let _permit=permit;
                    let deadline=gateway.config().limits.request_timeout+Duration::from_secs(45);
                    let service=service_fn(move|request:Request<Incoming>|{
                        let gateway=gateway.clone();
                        async move {
                            if let Err(error)=authorize(gateway.config(),request.method(),request.uri(),request.headers()){return Ok::<_,Infallible>(response(wire::Response{id:0,result:Err(error)}));}
                            let result=tokio::time::timeout(Duration::from_secs(15),read_request(request)).await;
                            let result=match result {
                                Ok(Ok(request))=>gateway.handle_http(request,Cancellation::default()).await,
                                Ok(Err(error))=>response(wire::Response{id:0,result:Err(error)}),
                                Err(_)=>response(wire::Response{id:0,result:Err(Failure::new(ErrorCode::Timeout,"Source-import request body timed out."))}),
                            };
                            Ok::<_,Infallible>(result)
                        }
                    });
                    let mut builder=http1::Builder::new();
                    builder.keep_alive(false).max_headers(32).max_buf_size(16*1024).timer(TokioTimer::new()).header_read_timeout(Duration::from_secs(10));
                    let _=tokio::time::timeout(deadline,builder.serve_connection(TokioIo::new(stream),service)).await;
                });
            }
            _=tasks.join_next(),if !tasks.is_empty()=>{},
        }
    }
    tasks.abort_all();
    while tasks.join_next().await.is_some() {}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::formats::sequence::ImportReport;
    use collagene_model::{Molecule, Topology};
    use collagene_sources::{AddgeneCatalog, SourceError};
    use std::sync::atomic::{AtomicUsize, Ordering};
    #[derive(Default)]
    struct Fake {
        calls: AtomicUsize,
    }
    impl Backend for Arc<Fake> {
        async fn import(
            &self,
            _source: String,
            _token: Option<String>,
        ) -> Result<ImportReport, SourceError> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            Ok(ImportReport {
                molecules: vec![
                    Molecule::dna("Provider record", "ACGT", Topology::Linear).unwrap(),
                ],
                warnings: vec![],
            })
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
    fn request(origin: &str, host: &str, source: &str) -> Request<Bytes> {
        let packet = wire::Request {
            id: 17,
            client: [7; 32],
            operation: wire::Operation::Import {
                source: source.into(),
                addgene_token: None,
            },
        };
        Request::builder()
            .method(Method::POST)
            .uri(PATH)
            .header(header::HOST, host)
            .header(header::ORIGIN, origin)
            .header(header::CONTENT_TYPE, CONTENT_TYPE)
            .body(Bytes::from(wire::encode_request(&packet).unwrap()))
            .unwrap()
    }
    #[tokio::test]
    async fn same_origin_binary_route_reuses_native_document_dto_without_network() {
        let backend = Arc::new(Fake::default());
        let gateway = Gateway::with_backend(Config::default(), backend.clone());
        let response = gateway
            .handle_http(
                request("http://127.0.0.1:8084", "127.0.0.1:8084", "ena:AB000263.1"),
                Cancellation::default(),
            )
            .await;
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(response.headers()[header::CACHE_CONTROL], "no-store");
        assert!(
            response
                .headers()
                .get("access-control-allow-origin")
                .is_none()
        );
        let body = response.into_body().collect().await.unwrap().to_bytes();
        let result = wire::decode_response(&body).unwrap();
        assert_eq!(result.id, 17);
        let wire::Output::Imported(report) = result.result.unwrap() else {
            panic!("document");
        };
        assert_eq!(report.molecules[0].sequence.as_str(), "ACGT");
        assert_eq!(backend.calls.load(Ordering::Relaxed), 1);
    }
    #[tokio::test]
    async fn cross_origin_wrong_host_direct_fetch_and_unbounded_body_do_not_reach_backend() {
        let backend = Arc::new(Fake::default());
        let gateway = Gateway::with_backend(Config::default(), backend.clone());
        for req in [
            request("https://other.example", "127.0.0.1:8084", "AB000263.1"),
            request("http://127.0.0.1:8084", "other.example", "AB000263.1"),
            request(
                "http://127.0.0.1:8084",
                "127.0.0.1:8084",
                "https://files.example.org/arbitrary.gb",
            ),
            request(
                "http://127.0.0.1:8084",
                "127.0.0.1:8084",
                "file:///etc/passwd",
            ),
        ] {
            assert!(
                !gateway
                    .handle_http(req, Cancellation::default())
                    .await
                    .status()
                    .is_success()
            );
        }
        let mut too_big = request("http://127.0.0.1:8084", "127.0.0.1:8084", "AB000263.1");
        *too_big.body_mut() = Bytes::from(vec![0; wire::MAX_REQUEST + 1]);
        assert_eq!(
            gateway
                .handle_http(too_big, Cancellation::default())
                .await
                .status(),
            StatusCode::PAYLOAD_TOO_LARGE
        );
        assert_eq!(backend.calls.load(Ordering::Relaxed), 0);
    }
    #[tokio::test]
    async fn configured_access_capability_is_checked_before_decoding_or_fetching() {
        let backend = Arc::new(Fake::default());
        let config = Config::remote(
            "0.0.0.0:8084".parse().unwrap(),
            vec!["https://collagene.example".into()],
            "x".repeat(32),
        )
        .unwrap();
        let gateway = Gateway::with_backend(config, backend.clone());
        let request_without = request(
            "https://collagene.example",
            "collagene.example",
            "AB000263.1",
        );
        assert_eq!(
            gateway
                .handle_http(request_without, Cancellation::default())
                .await
                .status(),
            StatusCode::FORBIDDEN
        );
        let mut request_with = request(
            "https://collagene.example",
            "collagene.example",
            "AB000263.1",
        );
        request_with.headers_mut().insert(
            header::AUTHORIZATION,
            HeaderValue::from_str(&format!("Bearer {}", "x".repeat(32))).unwrap(),
        );
        assert_eq!(
            gateway
                .handle_http(request_with, Cancellation::default())
                .await
                .status(),
            StatusCode::OK
        );
        assert_eq!(backend.calls.load(Ordering::Relaxed), 1);
    }
}
