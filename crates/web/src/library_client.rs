//! Same-origin shared-library transport. Credentials live only in this session;
//! signing out aborts pending requests and invalidates their eventual replies.
use collagene_library_protocol as wire;
use std::{cell::RefCell, collections::BTreeMap, rc::Rc};
use std::{future::Future, pin::Pin};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::JsFuture;

const ENDPOINT: &str = "/api/library";
const REQUEST_TIMEOUT_MS: i32 = 60_000;
const MAX_FRAME_BYTES: usize = wire::HEADER_BYTES + wire::MAX_MESSAGE_BYTES;

/// The folder UI consumes this boundary. Browser requests use the configured
/// server; Tauri dispatches to its embedded library without a network login.
pub trait LibraryClient {
    fn is_local(&self) -> bool {
        false
    }
    fn remote_key(&self) -> String;
    fn context(&self) -> Option<wire::LabContext>;
    /// Changes whenever this login is invalidated, including a reconnect to the same lab.
    fn publication_epoch(&self) -> u64;
    fn connect(&self, key: String) -> Result<(), String>;
    fn disconnect(&self);
    fn request(
        &self,
        operation: wire::Request,
    ) -> Pin<Box<dyn Future<Output = Result<wire::ResponseEnvelope, String>> + '_>>;
}

pub type Client = Rc<dyn LibraryClient>;
thread_local! {
    // One app-scoped login; source import borrows this session without ever
    // receiving its bearer value. Library panes created by browser() share it.
    static BROWSER_SESSION: Rc<Session> = Rc::new(Session::default());
}
pub(crate) fn browser_session() -> Rc<Session> {
    BROWSER_SESSION.with(Clone::clone)
}
pub fn browser() -> Client {
    if crate::desktop_rpc::is_desktop() {
        Rc::new(DesktopSession::default())
    } else {
        browser_session()
    }
}

#[derive(Clone, Default)]
pub struct Session(Rc<RefCell<State>>);

#[derive(Default)]
struct State {
    epoch: u64,
    next_request: u64,
    key: String,
    context: Option<wire::LabContext>,
    pending: BTreeMap<u64, web_sys::AbortController>,
}

impl LibraryClient for Session {
    fn remote_key(&self) -> String {
        web_sys::window()
            .and_then(|window| window.location().origin().ok())
            .unwrap_or_default()
    }
    fn context(&self) -> Option<wire::LabContext> {
        Session::context(self)
    }
    fn publication_epoch(&self) -> u64 {
        self.0.borrow().epoch
    }
    fn connect(&self, key: String) -> Result<(), String> {
        Session::connect(self, key)
    }
    fn disconnect(&self) {
        Session::disconnect(self);
    }
    fn request(
        &self,
        operation: wire::Request,
    ) -> Pin<Box<dyn Future<Output = Result<wire::ResponseEnvelope, String>> + '_>> {
        Box::pin(Session::request(self, operation))
    }
}

impl Session {
    /// A connection is confirmed only by a successfully validated server reply.
    pub fn connect(&self, key: String) -> Result<(), String> {
        let key = key.trim();
        if key.len() != 64 || !key.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err("Enter the 64-character lab access key.".into());
        }
        self.disconnect();
        self.0.borrow_mut().key = key.to_owned();
        Ok(())
    }

    pub fn disconnect(&self) {
        let mut state = self.0.borrow_mut();
        state.epoch = state.epoch.wrapping_add(1);
        state.key.clear();
        state.context = None;
        for (_, controller) in std::mem::take(&mut state.pending) {
            controller.abort();
        }
    }

    pub fn context(&self) -> Option<wire::LabContext> {
        self.0.borrow().context.clone()
    }

    /// Authenticated same-origin source gateway. The caller owns a dialog abort
    /// signal; signing out also aborts this fetch through the shared pending map.
    pub async fn source_request(
        &self,
        bytes: &[u8],
        external: &web_sys::AbortSignal,
    ) -> Result<Vec<u8>, String> {
        use collagene_source_gateway::wire as source;
        const MEDIA: &str = "application/vnd.collagene.source+postcard";
        if bytes.len() > source::MAX_REQUEST || bytes.get(..4) != Some(b"CGSQ".as_slice()) {
            return Err("Invalid bounded source request.".into());
        }
        if external.aborted() {
            return Err("Source request cancelled.".into());
        }
        let controller =
            web_sys::AbortController::new().map_err(|_| "Cannot start source request")?;
        let (id, epoch, key, context) = {
            let mut state = self.0.borrow_mut();
            let context = state
                .context
                .clone()
                .ok_or("Connect to the lab library first, then retry source import.")?;
            if state.key.is_empty() {
                return Err("Connect to the lab library first, then retry source import.".into());
            }
            state.next_request = state
                .next_request
                .checked_add(1)
                .ok_or("Library request IDs exhausted; reload the app.")?;
            let id = state.next_request;
            state.pending.insert(id, controller.clone());
            (id, state.epoch, state.key.clone(), context)
        };
        let mut lifetime = RequestLifetime {
            state: self.0.clone(),
            request_id: id,
            controller: controller.clone(),
            timeout: None,
        };
        let abort = controller.clone();
        let listener = Closure::<dyn FnMut()>::new(move || abort.abort());
        external
            .add_event_listener_with_callback("abort", listener.as_ref().unchecked_ref())
            .map_err(|_| "Cannot attach source cancellation")?;
        let _external = ExternalAbort {
            signal: external.clone(),
            listener,
        };
        if external.aborted() {
            controller.abort();
        }
        let abort = controller.clone();
        let timeout = Closure::<dyn FnMut()>::new(move || abort.abort());
        let window = web_sys::window().ok_or("No browser window")?;
        let timer = window
            .set_timeout_with_callback_and_timeout_and_arguments_0(
                timeout.as_ref().unchecked_ref(),
                90_000,
            )
            .map_err(|_| "Cannot bound source request duration")?;
        lifetime.timeout = Some((timer, timeout));
        let init = web_sys::RequestInit::new();
        init.set_method("POST");
        init.set_mode(web_sys::RequestMode::SameOrigin);
        init.set_credentials(web_sys::RequestCredentials::Omit);
        init.set_redirect(web_sys::RequestRedirect::Error);
        init.set_signal(Some(&controller.signal()));
        init.set_body(&js_sys::Uint8Array::from(bytes));
        let http = web_sys::Request::new_with_str_and_init("/api/source-import/v1", &init)
            .map_err(|_| "Cannot prepare source gateway request")?;
        http.headers()
            .set("Content-Type", MEDIA)
            .map_err(|_| "Cannot set source content type")?;
        http.headers()
            .set("Accept", MEDIA)
            .map_err(|_| "Cannot set source response type")?;
        http.headers()
            .set("Authorization", &format!("Bearer {key}"))
            .map_err(|_| "Cannot authorize source request")?;
        http.headers()
            .set(wire::EXPECTED_LAB_HEADER, context.lab_id.as_str())
            .map_err(|_| "Cannot pin source gateway lab")?;
        drop(key);
        let response = JsFuture::from(window.fetch_with_request(&http))
            .await
            .map_err(|_| "Source request interrupted or timed out.")?
            .dyn_into::<web_sys::Response>()
            .map_err(|_| "Invalid source gateway response")?;
        if self.0.borrow().epoch != epoch {
            return Err("The library login changed; source response discarded.".into());
        }
        match response.status() {
            401 | 412 => {
                self.disconnect();
                return Err(
                    "The lab login or database changed. Reconnect to the lab library first.".into(),
                );
            }
            403 => return Err("Source import is not permitted by this lab server.".into()),
            404 | 405 => return Err("Source import is not enabled on this server.".into()),
            413 => return Err("Source gateway request exceeds its limit.".into()),
            429 | 503 => return Err("Source gateway is busy; wait briefly before retrying.".into()),
            _ if !response.ok() => {
                return Err("Source gateway request failed; no document was imported.".into());
            }
            _ => {}
        }
        if response
            .headers()
            .get("content-type")
            .ok()
            .flatten()
            .is_none_or(|s| s.split(';').next().unwrap_or("").trim() != MEDIA)
        {
            return Err("Unexpected source gateway response type.".into());
        }
        let bytes = read_bounded_limit(&response, source::MAX_RESPONSE)
            .await
            .map_err(|_| "Source response was interrupted or exceeded its 64 MiB bound.")?;
        if self.0.borrow().epoch != epoch || external.aborted() {
            return Err("Source response discarded after cancellation or login change.".into());
        }
        Ok(bytes)
    }

    pub async fn request(
        &self,
        operation: wire::Request,
    ) -> Result<wire::ResponseEnvelope, String> {
        crate::library_capture::ensure_public_document_request(&operation)?;
        let controller =
            web_sys::AbortController::new().map_err(|_| "Cannot start a library request")?;
        let (request_id, epoch, key, context) = {
            let mut state = self.0.borrow_mut();
            if state.key.is_empty() {
                return Err("Connect to the lab library first.".into());
            }
            let mutation = matches!(
                &operation,
                wire::Request::CreateFolder { .. }
                    | wire::Request::RenameFolder { .. }
                    | wire::Request::MoveFolder { .. }
                    | wire::Request::DeleteEmptyFolder { .. }
                    | wire::Request::ImportDocument { .. }
                    | wire::Request::SaveRevision { .. }
                    | wire::Request::RenameDocument { .. }
                    | wire::Request::MoveDocument { .. }
            );
            if mutation && state.context.is_none() {
                return Err("Load and review the selected lab before sending changes.".into());
            }
            state.next_request = state.next_request.checked_add(1).ok_or(
                "This library session has exhausted its request identifiers; reload the app.",
            )?;
            let id = state.next_request;
            state.pending.insert(id, controller.clone());
            (id, state.epoch, state.key.clone(), state.context.clone())
        };
        let mut lifetime = RequestLifetime {
            state: self.0.clone(),
            request_id,
            controller: controller.clone(),
            timeout: None,
        };
        let abort = controller.clone();
        let on_timeout = Closure::<dyn FnMut()>::new(move || abort.abort());
        let window = web_sys::window().ok_or("No browser window")?;
        let timeout = window
            .set_timeout_with_callback_and_timeout_and_arguments_0(
                on_timeout.as_ref().unchecked_ref(),
                REQUEST_TIMEOUT_MS,
            )
            .map_err(|_| "Cannot bound the library request duration")?;
        lifetime.timeout = Some((timeout, on_timeout));
        let request = wire::RequestEnvelope {
            request_id,
            request: operation,
        };
        let bytes = wire::encode_request(&request).map_err(|error| error.to_string())?;
        let body = js_sys::Uint8Array::from(bytes.as_slice());
        let init = web_sys::RequestInit::new();
        init.set_method("POST");
        init.set_mode(web_sys::RequestMode::SameOrigin);
        init.set_credentials(web_sys::RequestCredentials::Omit);
        init.set_redirect(web_sys::RequestRedirect::Error);
        init.set_signal(Some(&controller.signal()));
        init.set_body(&body);
        let http = web_sys::Request::new_with_str_and_init(ENDPOINT, &init)
            .map_err(|_| "Cannot prepare the library request")?;
        http.headers()
            .set("Content-Type", wire::MEDIA_TYPE)
            .map_err(|_| "Cannot prepare the library content type")?;
        http.headers()
            .set("Accept", wire::MEDIA_TYPE)
            .map_err(|_| "Cannot prepare the library response type")?;
        http.headers()
            .set("Authorization", &format!("Bearer {key}"))
            .map_err(|_| "Cannot prepare lab authorization")?;
        if let Some(context) = &context {
            http.headers()
                .set(wire::EXPECTED_LAB_HEADER, context.lab_id.as_str())
                .map_err(|_| "Cannot pin the requested lab database")?;
        }
        drop(key);
        let response = JsFuture::from(window.fetch_with_request(&http))
            .await
            .map_err(|_| "Library request interrupted or timed out; local work is unchanged.")?
            .dyn_into::<web_sys::Response>()
            .map_err(|_| "Invalid library response")?;
        if self.0.borrow().epoch != epoch {
            return Err("The library session changed; this reply was discarded.".into());
        }
        match response.status() {
            401 => {
                self.disconnect();
                return Err("The lab access key was not accepted. Connect again.".into());
            }
            412 => {
                self.disconnect();
                return Err("The server's lab changed. Reconnect and review the selected database before saving.".into());
            }
            428 => {
                return Err(
                    "The server requires a database identity pin. Reload this app before saving."
                        .into(),
                );
            }
            404 | 405 => return Err("Shared library access is not enabled on this server.".into()),
            413 => return Err("The document exceeds this server's library size limit.".into()),
            429 | 503 => {
                return Err(
                    "The library server is busy. Retry without changing your local work.".into(),
                );
            }
            _ if !response.ok() => {
                return Err(format!(
                    "Library request failed (HTTP {}); local work is unchanged.",
                    response.status()
                ));
            }
            _ => {}
        }
        if response
            .headers()
            .get("content-type")
            .ok()
            .flatten()
            .is_none_or(|value| value.split(';').next().unwrap_or("").trim() != wire::MEDIA_TYPE)
        {
            return Err("The server returned an unexpected library response type.".into());
        }
        let bytes = read_bounded(&response).await?;
        let reply = wire::decode_response_for(&bytes, &request)
            .map_err(|error| format!("Library reply was rejected: {error}"))?;
        if let Some(context) = &context {
            reply
                .validate_for_in_context(&request, &context.lab_id)
                .map_err(|_| "The server's lab changed; reconnect before continuing.")?;
        }
        let mut state = self.0.borrow_mut();
        if state.epoch != epoch {
            return Err("The library session changed; this reply was discarded.".into());
        }
        if let Some(previous) = &state.context {
            if reply
                .context
                .as_ref()
                .is_none_or(|context| context.lab_id != previous.lab_id)
            {
                return Err("The server's lab changed; reconnect before continuing.".into());
            }
        } else if let Some(context) = &reply.context {
            state.context = Some(context.clone());
        }
        Ok(reply)
    }
}

struct RequestLifetime {
    state: Rc<RefCell<State>>,
    request_id: u64,
    controller: web_sys::AbortController,
    timeout: Option<(i32, Closure<dyn FnMut()>)>,
}
struct ExternalAbort {
    signal: web_sys::AbortSignal,
    listener: Closure<dyn FnMut()>,
}
impl Drop for ExternalAbort {
    fn drop(&mut self) {
        let _ = self
            .signal
            .remove_event_listener_with_callback("abort", self.listener.as_ref().unchecked_ref());
    }
}
impl Drop for RequestLifetime {
    fn drop(&mut self) {
        self.controller.abort();
        self.state.borrow_mut().pending.remove(&self.request_id);
        if let Some((id, _)) = self.timeout.take() {
            if let Some(window) = web_sys::window() {
                window.clear_timeout_with_handle(id);
            }
        }
    }
}

async fn read_bounded(response: &web_sys::Response) -> Result<Vec<u8>, String> {
    read_bounded_limit(response, MAX_FRAME_BYTES).await
}
async fn read_bounded_limit(response: &web_sys::Response, limit: usize) -> Result<Vec<u8>, String> {
    if response
        .headers()
        .get("content-length")
        .ok()
        .flatten()
        .and_then(|value| value.parse::<usize>().ok())
        .is_some_and(|size| size > limit)
    {
        return Err("The library response exceeds its size bound.".into());
    }
    let reader = response
        .body()
        .ok_or("Library response has no body")?
        .get_reader()
        .dyn_into::<web_sys::ReadableStreamDefaultReader>()
        .map_err(|_| "Cannot read the library response")?;
    let mut bytes = Vec::new();
    loop {
        let result = JsFuture::from(reader.read())
            .await
            .map_err(|_| "Library response was interrupted; local work is unchanged.")?;
        if js_sys::Reflect::get(&result, &"done".into())
            .ok()
            .and_then(|value| value.as_bool())
            .unwrap_or(false)
        {
            break;
        }
        let value = js_sys::Reflect::get(&result, &"value".into())
            .map_err(|_| "Invalid library response chunk")?;
        let chunk = js_sys::Uint8Array::new(&value);
        if chunk.length() as usize > limit.saturating_sub(bytes.len()) {
            let _ = reader.cancel();
            return Err("The library response exceeds its size bound.".into());
        }
        bytes.extend(chunk.to_vec());
    }
    reader.release_lock();
    Ok(bytes)
}

/// Offline local catalog. The native shell supplies repository authority; no
/// lab access key or network endpoint is accepted by this adapter.
#[derive(Default)]
struct DesktopSession {
    context: RefCell<Option<wire::LabContext>>,
    epoch: std::cell::Cell<u64>,
}
impl LibraryClient for DesktopSession {
    fn is_local(&self) -> bool {
        true
    }
    fn remote_key(&self) -> String {
        "collagene-desktop://local".into()
    }
    fn context(&self) -> Option<wire::LabContext> {
        self.context.borrow().clone()
    }
    fn publication_epoch(&self) -> u64 {
        self.epoch.get()
    }
    fn connect(&self, _key: String) -> Result<(), String> {
        Ok(())
    }
    fn disconnect(&self) {
        self.epoch.set(self.epoch.get().wrapping_add(1));
    }
    fn request(
        &self,
        operation: wire::Request,
    ) -> Pin<Box<dyn Future<Output = Result<wire::ResponseEnvelope, String>> + '_>> {
        Box::pin(async move {
            use collagene_desktop_protocol as desktop;
            let epoch = self.epoch.get();
            if self.context.borrow().is_none() {
                let bootstrap = crate::desktop_rpc::call(desktop::Command::Bootstrap).await?;
                if epoch != self.epoch.get() {
                    return Err("Local library view closed".into());
                }
                match bootstrap {
                    desktop::Success::Bootstrap(b) => {
                        *self.context.borrow_mut() = Some(b.library_context.ok_or(
                            "Saved library is paused; inspect native recovery status. Working-editor autosave remains available."
                        )?)
                    }
                    _ => return Err("Unexpected local library bootstrap".into()),
                }
            }
            let request = crate::desktop_rpc::request(desktop::Command::Library(operation))?;
            let outcome = crate::desktop_rpc::dispatch(&request).await?;
            if epoch != self.epoch.get() {
                return Err("Local library view closed".into());
            }
            let outcome = match outcome {
                desktop::Outcome::Success(desktop::Success::Library(success)) => {
                    wire::Outcome::Success(success)
                }
                desktop::Outcome::Failure(desktop::Failure::Library(failure)) => {
                    wire::Outcome::Failure(failure)
                }
                desktop::Outcome::Failure(failure) => {
                    return Err(crate::desktop_rpc::failure(&failure));
                }
                _ => return Err("Unexpected local library response".into()),
            };
            Ok(wire::ResponseEnvelope {
                request_id: request.request_id,
                context: self.context(),
                outcome,
            })
        })
    }
}
