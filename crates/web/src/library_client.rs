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

/// The folder UI consumes this boundary. A future offline Tauri adapter can
/// dispatch the same messages to its embedded library without a network login.
pub trait LibraryClient {
    fn remote_key(&self) -> String;
    fn context(&self) -> Option<wire::LabContext>;
    fn connect(&self, key: String) -> Result<(), String>;
    fn disconnect(&self);
    fn request(
        &self,
        operation: wire::Request,
    ) -> Pin<Box<dyn Future<Output = Result<wire::ResponseEnvelope, String>> + '_>>;
}

pub type Client = Rc<dyn LibraryClient>;
pub fn browser() -> Client {
    Rc::new(Session::default())
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

    pub async fn request(
        &self,
        operation: wire::Request,
    ) -> Result<wire::ResponseEnvelope, String> {
        let controller =
            web_sys::AbortController::new().map_err(|_| "Cannot start a library request")?;
        let (request_id, epoch, key, context) = {
            let mut state = self.0.borrow_mut();
            if state.key.is_empty() {
                return Err("Connect to the lab library first.".into());
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
    if response
        .headers()
        .get("content-length")
        .ok()
        .flatten()
        .and_then(|value| value.parse::<usize>().ok())
        .is_some_and(|size| size > MAX_FRAME_BYTES)
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
        if chunk.length() as usize > MAX_FRAME_BYTES.saturating_sub(bytes.len()) {
            let _ = reader.cancel();
            return Err("The library response exceeds its size bound.".into());
        }
        bytes.extend(chunk.to_vec());
    }
    reader.release_lock();
    Ok(bytes)
}
