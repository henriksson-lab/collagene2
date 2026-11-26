//! Rust/WASM adapter for Tauri's public raw IPC API. No network or JS helpers.
use collagene_desktop_protocol as wire;
use std::cell::Cell;
use wasm_bindgen::{JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
thread_local! { static NEXT: Cell<u64> = const { Cell::new(0) }; }
thread_local! {
    static LANE: std::rc::Rc<futures_util::lock::Mutex<()>> =
        std::rc::Rc::new(futures_util::lock::Mutex::new(()));
}

fn invoke() -> Result<(JsValue, js_sys::Function), String> {
    let window = web_sys::window().ok_or("No application window")?;
    let tauri = js_sys::Reflect::get(&window, &"__TAURI__".into())
        .map_err(|_| "Desktop IPC unavailable")?;
    let core =
        js_sys::Reflect::get(&tauri, &"core".into()).map_err(|_| "Desktop IPC unavailable")?;
    let function = js_sys::Reflect::get(&core, &"invoke".into())
        .map_err(|_| "Desktop IPC unavailable")?
        .dyn_into::<js_sys::Function>()
        .map_err(|_| "Desktop IPC unavailable")?;
    Ok((core, function))
}
pub fn available() -> bool {
    invoke().is_ok()
}
/// A native window must fail closed when its IPC initialization is broken;
/// falling back to an ephemeral WebView database would lose work on restart.
pub fn is_desktop() -> bool {
    if available() {
        return true;
    }
    let Some(window) = web_sys::window() else {
        return false;
    };
    let location = window.location();
    matches!(
        (
            location.protocol().ok().as_deref(),
            location.hostname().ok().as_deref()
        ),
        (Some("tauri:"), Some("localhost")) | (Some("https:"), Some("tauri.localhost"))
    )
}

pub fn request(command: wire::Command) -> Result<wire::RequestEnvelope, String> {
    let request_id = NEXT.with(|next| {
        let id = next
            .get()
            .checked_add(1)
            .ok_or("Desktop request IDs exhausted")?;
        next.set(id);
        Ok::<_, &str>(id)
    })?;
    // Random operation identities survive overlapping views and distinguish a
    // reopened WebView from a previous session's request counter.
    let random = js_sys::Uint8Array::new_with_length(16);
    let window = web_sys::window().ok_or("No application window")?;
    let crypto =
        js_sys::Reflect::get(&window, &"crypto".into()).map_err(|_| "No secure random source")?;
    let get = js_sys::Reflect::get(&crypto, &"getRandomValues".into())
        .map_err(|_| "No secure random source")?
        .dyn_into::<js_sys::Function>()
        .map_err(|_| "No secure random source")?;
    get.call1(&crypto, &random)
        .map_err(|_| "Cannot create desktop operation identity")?;
    let mut id = String::from("op-");
    use std::fmt::Write;
    for byte in random.to_vec() {
        let _ = write!(id, "{byte:02x}");
    }
    Ok(wire::RequestEnvelope {
        request_id,
        operation_id: wire::OperationId::new(id).map_err(|e| e.to_string())?,
        command,
    })
}

pub async fn dispatch(request: &wire::RequestEnvelope) -> Result<wire::Outcome, String> {
    dispatch_guarded(request, || Ok(())).await
}

pub async fn dispatch_guarded(
    request: &wire::RequestEnvelope,
    guard: impl Fn() -> Result<(), String>,
) -> Result<wire::Outcome, String> {
    // All native catalog, working-file and preference operations share the
    // host's single admission lane. Cancellation must bypass it so a pending
    // replacement can be stopped while another request occupies the host.
    let lane = LANE.with(Clone::clone);
    let _guard = if matches!(
        request.command,
        wire::Command::Cancel { .. } | wire::Command::OperationStatus { .. }
    ) {
        None
    } else {
        Some(lane.lock().await)
    };
    guard()?;
    let bytes = wire::encode_request(request).map_err(|e| e.to_string())?;
    let (core, function) = invoke()?;
    let promise = function
        .call2(
            &core,
            &"collagene_rpc".into(),
            &js_sys::Uint8Array::from(bytes.as_slice()),
        )
        .map_err(|_| "Cannot dispatch desktop operation")?
        .dyn_into::<js_sys::Promise>()
        .map_err(|_| "Desktop IPC did not return a promise")?;
    let result = JsFuture::from(promise)
        .await
        .map_err(|_| "The local repository reply was interrupted. Preserve unsaved changes before reloading.")?;
    if !result.is_instance_of::<js_sys::ArrayBuffer>() {
        return Err("Desktop IPC did not return binary data".into());
    }
    let bytes = js_sys::Uint8Array::new(&result);
    if bytes.length() as usize > wire::HEADER_BYTES + wire::MAX_RESPONSE_BYTES {
        return Err("Desktop response exceeds its limit".into());
    }
    Ok(wire::decode_response_for(&bytes.to_vec(), request)
        .map_err(|e| e.to_string())?
        .outcome)
}
pub async fn call(command: wire::Command) -> Result<wire::Success, String> {
    match dispatch(&request(command)?).await? {
        wire::Outcome::Success(success) => Ok(success),
        wire::Outcome::Failure(error) => Err(failure(&error)),
    }
}
pub fn failure(error: &wire::Failure) -> String {
    match error {
        wire::Failure::Invalid { message } => message.to_string(),
        wire::Failure::Conflict { .. } => "Local files changed in another window. Reload before saving; export your current work to preserve it.".into(),
        wire::Failure::Cancelled => "Desktop operation cancelled".into(),
        wire::Failure::Unavailable => "Local repository is unavailable".into(),
        wire::Failure::Library(failure) => failure.message.to_string(),
    }
}
