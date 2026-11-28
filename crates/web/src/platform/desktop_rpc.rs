//! Rust/WASM adapter for Tauri's public raw IPC API. No network or JS helpers.
use collagene_desktop_protocol as wire;
use std::cell::Cell;
use wasm_bindgen::{JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
thread_local! { static NEXT: Cell<u64> = const { Cell::new(0) }; }
thread_local! {
    static EXPORT_ACTIVE: Cell<bool> = const { Cell::new(false) };
    static EXPORT_LANE: std::rc::Rc<futures_util::lock::Mutex<()>> =
        std::rc::Rc::new(futures_util::lock::Mutex::new(()));
    static SOURCE_LANE: std::rc::Rc<futures_util::lock::Mutex<()>> =
        std::rc::Rc::new(futures_util::lock::Mutex::new(()));
    static SYNC_LANE: std::rc::Rc<futures_util::lock::Mutex<()>> =
        std::rc::Rc::new(futures_util::lock::Mutex::new(()));
    static LANE: std::rc::Rc<futures_util::lock::Mutex<()>> =
        std::rc::Rc::new(futures_util::lock::Mutex::new(()));
}

/// Reserve before copying an artifact: at most one large export is resident.
pub struct ExportReservation {
    _private: (),
}
impl Drop for ExportReservation {
    fn drop(&mut self) {
        EXPORT_ACTIVE.with(|active| active.set(false));
    }
}
pub fn reserve_export() -> Result<ExportReservation, String> {
    EXPORT_ACTIVE.with(|active| {
        if active.replace(true) {
            Err(
                "Another export is still being saved. Wait for its receipt before exporting again."
                    .into(),
            )
        } else {
            Ok(ExportReservation { _private: () })
        }
    })
}
pub async fn export_reserved(
    _reservation: ExportReservation,
    name: String,
    bytes: Vec<u8>,
    guard: impl Fn() -> Result<(), String>,
) -> Result<String, String> {
    let request = request(wire::Command::ExportFile {
        suggested_name: wire::library::Text::new(name)
            .map_err(|_| "Export filename exceeds its bound")?,
        bytes: wire::library::Blob::new(bytes).map_err(|_| "Export exceeds 64 MiB")?,
    })?;
    let lane = EXPORT_LANE.with(Clone::clone);
    let _lane = lane.lock().await;
    guard()?;
    let bytes = wire::encode_request(&request).map_err(|e| e.to_string())?;
    let bytes = invoke_frame(&bytes, wire::HEADER_BYTES + wire::MAX_RESPONSE_BYTES).await?;
    match wire::decode_response_for(&bytes, &request)
        .map_err(|e| e.to_string())?
        .outcome
    {
        wire::Outcome::Success(wire::Success::Exported { name, .. }) => Ok(name.into_string()),
        wire::Outcome::Failure(error) => Err(failure(&error)),
        _ => Err("Native export returned an unexpected receipt.".into()),
    }
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
    // Local catalog, working-file and preference operations serialize. Remote
    // HTTP sync uses a separate lane so it cannot delay local autosaves.
    // Cancellation bypasses both lanes.
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
    let bytes = invoke_frame(&bytes, wire::HEADER_BYTES + wire::MAX_RESPONSE_BYTES).await?;
    Ok(wire::decode_response_for(&bytes, request)
        .map_err(|e| e.to_string())?
        .outcome)
}
async fn invoke_frame(bytes: &[u8], response_limit: usize) -> Result<Vec<u8>, String> {
    let (core, function) = invoke()?;
    let promise = function
        .call2(
            &core,
            &"collagene_rpc".into(),
            &js_sys::Uint8Array::from(bytes),
        )
        .map_err(|_| "Cannot dispatch desktop operation")?
        .dyn_into::<js_sys::Promise>()
        .map_err(|_| "Desktop IPC did not return a promise")?;
    let result = JsFuture::from(promise).await.map_err(|_| "The local repository reply was interrupted. Preserve unsaved changes before reloading.")?;
    if !result.is_instance_of::<js_sys::ArrayBuffer>() {
        return Err("Desktop IPC did not return binary data".into());
    }
    let bytes = js_sys::Uint8Array::new(&result);
    if bytes.length() as usize > response_limit {
        return Err("Desktop response exceeds its limit".into());
    }
    Ok(bytes.to_vec())
}
/// Metadata-only remote operations have their own bounded native lane. They can
/// wait for HTTPS while local working-file and catalog operations keep running.
pub fn sync_client() -> crate::platform::desktop_sync_client::Client {
    use collagene_desktop_sync_protocol as sync;
    crate::platform::desktop_sync_client::Client::new(
        std::rc::Rc::new(|command| {
            let pin = request(wire::Command::Bootstrap)?;
            Ok(sync::RequestEnvelope {
                request_id: pin.request_id,
                operation_id: pin.operation_id,
                command,
            })
        }),
        std::rc::Rc::new(|bytes| {
            Box::pin(async move {
                if !is_desktop()
                    || bytes.len() > sync::MAX_REQUEST_BYTES
                    || bytes.get(..4) != Some(b"DSQ\0".as_slice())
                {
                    return Err("Invalid native database request".into());
                }
                let lane = SYNC_LANE.with(Clone::clone);
                let _guard = lane.lock().await;
                invoke_frame(&bytes, sync::MAX_RESPONSE_BYTES).await
            })
        }),
    )
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

/// Source controls must be able to cancel an active fetch while ordinary source
/// requests serialize. An aborted queued request is never sent to the native host.
pub async fn source_raw(
    bytes: Vec<u8>,
    before_dispatch: Option<web_sys::AbortSignal>,
) -> Result<Vec<u8>, String> {
    use collagene_source_gateway::wire as source;
    if !is_desktop() || bytes.len() > source::MAX_REQUEST {
        return Err("Native source import is unavailable".into());
    }
    let request = source::decode_request(&bytes).map_err(|_| "Invalid source request")?;
    let control = matches!(
        request.operation,
        source::Operation::Cancel { .. } | source::Operation::DisposeCatalog { .. }
    );
    let lane = SOURCE_LANE.with(Clone::clone);
    let _guard = if control {
        None
    } else {
        Some(lane.lock().await)
    };
    if before_dispatch
        .as_ref()
        .is_some_and(web_sys::AbortSignal::aborted)
    {
        return Err("Source request cancelled before dispatch".into());
    }
    invoke_frame(&bytes, source::MAX_RESPONSE).await
}
