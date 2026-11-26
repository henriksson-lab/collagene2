#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
#[cfg(any(target_os = "android", target_os = "ios"))]
compile_error!("Collagene desktop supports Windows, macOS and Linux; mobile is deferred.");
mod policy;
mod settings;
use collagene_desktop_host::{
    DesktopHost, HostControl, JobStart, MAX_REQUEST_BYTES, MAX_RESPONSE_BYTES, OperationState,
    WindowAccess,
};
use collagene_desktop_sync_host::{ProductionAdapter, SyncHost};
use collagene_desktop_sync_protocol as sync;
use collagene_source_gateway::{self as sources, wire as source_wire};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicUsize, Ordering},
};
use tauri::Manager;

#[derive(Clone)]
struct NativeState {
    host: Arc<Mutex<DesktopHost>>,
    sync: Arc<Mutex<SyncHost<ProductionAdapter>>>,
    control: HostControl,
    access: WindowAccess,
    local_active: Arc<AtomicUsize>,
    sync_active: Arc<AtomicUsize>,
    source_active: Arc<AtomicUsize>,
    sources: Arc<sources::Gateway<sources::NativeBackend>>,
    source_scope: [u8; 32],
}
struct JobLease(Arc<AtomicUsize>);
impl JobLease {
    fn acquire(active: &Arc<AtomicUsize>) -> Result<Self, String> {
        active
            .compare_exchange(0, 1, Ordering::AcqRel, Ordering::Acquire)
            .map_err(|_| {
                "Native operation is busy; retain its operation ID before retrying.".to_string()
            })?;
        Ok(Self(active.clone()))
    }
}
impl Drop for JobLease {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::Release);
    }
}
fn bounded_response(bytes: Vec<u8>, limit: usize) -> Result<tauri::ipc::Response, String> {
    if bytes.len() > limit {
        return Err("Native reply exceeded its binary protocol limit.".into());
    }
    Ok(tauri::ipc::Response::new(bytes))
}
fn sync_failure(request: &sync::RequestEnvelope, cancelled: bool) -> sync::ResponseEnvelope {
    sync::ResponseEnvelope {
        request_id: request.request_id,
        operation_id: request.operation_id.clone(),
        context: None,
        outcome: sync::Outcome::Failure(sync::Failure::Error {
            code: if cancelled {
                sync::ErrorCode::Cancelled
            } else {
                sync::ErrorCode::Unavailable
            },
            message: sync::library::Text::new(if cancelled {
                "Database operation cancelled before acceptance."
            } else {
                "Database operation interrupted; review its durable job before retrying."
            })
            .expect("bounded static message"),
        }),
    }
}
fn dispatch_sync(state: &NativeState, bytes: &[u8]) -> Result<Vec<u8>, String> {
    let request =
        sync::decode_request(bytes).map_err(|_| "Invalid native database frame.".to_owned())?;
    // Hash the whole validated frame, including its DSQ namespace. DTQ uses the
    // same registry; identical payload tags from different protocols must not alias.
    let digest = sync::library::Digest::of(bytes);
    let job = match state
        .control
        .begin_job(&state.access, &request.operation_id, digest)
    {
        Ok(JobStart::Cached(reply)) => {
            sync::decode_response_for(&reply, &request)
                .map_err(|_| "Cached database reply did not match its request.".to_owned())?;
            return Ok(reply);
        }
        Ok(JobStart::Run(job)) => job,
        Err(error) => {
            let cancelled = matches!(error, collagene_desktop_host::HostError::Cancelled);
            return sync::encode_response(&sync_failure(&request, cancelled))
                .map_err(|_| "Cannot encode native database failure.".into());
        }
    };
    let result = state
        .sync
        .lock()
        .map_err(|_| "Native database service is unavailable.".to_owned())?
        .dispatch_raw(bytes, &mut || job.accept().is_ok());
    let mut response = match result {
        Ok(reply) => sync::decode_response_for(&reply, &request)
            .map_err(|_| "Native database reply did not match its request.".to_owned())?,
        Err(_) => sync_failure(&request, false),
    };
    if job.cancelled_before_acceptance() {
        response = sync_failure(&request, true);
    }
    response
        .validate_for(&request)
        .map_err(|_| "Native database response binding failed.".to_owned())?;
    let terminal = match &response.outcome {
        sync::Outcome::Success(_) => OperationState::Completed,
        sync::Outcome::Failure(sync::Failure::Error {
            code: sync::ErrorCode::Cancelled,
            ..
        }) => OperationState::Cancelled,
        _ => OperationState::Failed,
    };
    let reply = sync::encode_response(&response)
        .map_err(|_| "Cannot encode native database response.".to_owned())?;
    job.finish(terminal, &reply).map_err(|_| {
        "Native database status could not be retained; review its durable job.".to_owned()
    })?;
    Ok(reply)
}
struct CancelSourceOnDrop(sources::Cancellation);
impl Drop for CancelSourceOnDrop {
    fn drop(&mut self) {
        self.0.cancel();
    }
}
async fn dispatch_source(
    state: &NativeState,
    bytes: &[u8],
) -> Result<tauri::ipc::Response, String> {
    state
        .control
        .check_window(&state.access)
        .map_err(|_| "Native source window access is invalid".to_owned())?;
    let mut request =
        source_wire::decode_request(bytes).map_err(|_| "Invalid native source frame".to_owned())?;
    if request.client == [0; 32] {
        return Err("A private source session is required".into());
    }
    let control = matches!(
        request.operation,
        source_wire::Operation::Cancel { .. } | source_wire::Operation::DisposeCatalog { .. }
    );
    let _lease = if control {
        None
    } else {
        Some(JobLease::acquire(&state.source_active)?)
    };
    // A client nonce is correlation, never authority. Bind it to this native
    // host instance so catalogs cannot be replayed into a different instance.
    request.client = *blake3::keyed_hash(&state.source_scope, &request.client).as_bytes();
    let cancellation = sources::Cancellation::default();
    let _cancel_on_drop = CancelSourceOnDrop(cancellation.clone());
    let response = state.sources.execute(request, cancellation).await;
    let frame =
        tauri::async_runtime::spawn_blocking(move || source_wire::encode_response(&response))
            .await
            .map_err(|_| "Native source response was interrupted".to_owned())?
            .map_err(|_| "Native source response exceeds its format limit".to_owned())?;
    bounded_response(frame, source_wire::MAX_RESPONSE)
}
#[tauri::command]
async fn collagene_rpc(
    request: tauri::ipc::Request<'_>,
    window: tauri::WebviewWindow,
    state: tauri::State<'_, NativeState>,
) -> Result<tauri::ipc::Response, String> {
    if window.label() != "main" || !window.url().ok().as_ref().is_some_and(policy::trusted_url) {
        return Err("Native repository access is restricted to the embedded main window.".into());
    }
    let tauri::ipc::InvokeBody::Raw(bytes) = request.body() else {
        return Err("Native requests must use the binary repository protocol.".into());
    };
    if bytes.get(..4) == Some(b"CGSQ".as_slice()) {
        if bytes.len() > source_wire::MAX_REQUEST {
            return Err("Native source request exceeds its limit".into());
        }
        return dispatch_source(state.inner(), bytes).await;
    }
    let syncing = match bytes.get(..4) {
        Some(b"DSQ\0") => true,
        Some(b"DTQ\0") => false,
        _ => return Err("Unknown native binary protocol.".into()),
    };
    let limit = if syncing {
        sync::MAX_REQUEST_BYTES
    } else {
        MAX_REQUEST_BYTES
    };
    let response_limit = if syncing {
        sync::MAX_RESPONSE_BYTES
    } else {
        MAX_RESPONSE_BYTES
    };
    if bytes.len() > limit {
        return Err("Native request exceeds its binary protocol limit.".into());
    }
    if !syncing {
        if let Some(response) = state
            .control
            .dispatch_raw(&state.access, bytes)
            .map_err(|_| "Native control request was rejected.".to_string())?
        {
            return bounded_response(response, response_limit);
        }
    }
    // Remote HTTPS and local working-file writes have separate admission lanes
    // and separate host mutexes. Their shared catalog uses immutable-head and
    // checkpoint CAS, so network latency cannot stop local working-file saves.
    let lease = JobLease::acquire(if syncing {
        &state.sync_active
    } else {
        &state.local_active
    })?;
    let bytes = bytes.clone();
    let state = state.inner().clone();
    let response = tauri::async_runtime::spawn_blocking(move || {
        let _lease = lease;
        if syncing {
            return dispatch_sync(&state, &bytes);
        }
        let mut host = state
            .host
            .lock()
            .map_err(|_| "Native repository is unavailable; retain local work.".to_string())?;
        host.dispatch_raw(&state.access, &bytes).map_err(|_| {
            "Native request failed; retain its ID and review status before retrying.".to_string()
        })
    })
    .await
    .map_err(|_| "Native task could not complete; review its durable status.".to_string())??;
    bounded_response(response, response_limit)
}
fn run(options: settings::GuiOptions) -> Result<(), Box<dyn std::error::Error>> {
    let allow_loopback_http = options.allow_loopback_http;
    let storage_override = options.storage_root;
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![collagene_rpc])
        .setup(move |app| {
            let configuration = app.path().app_config_dir()?;
            let default_root = app.path().app_data_dir()?.join("repository");
            let root = settings::select(&configuration, &default_root, storage_override)?;
            policy::separate_assets(&root)?;
            let mut host = DesktopHost::open(&root)?;
            let access = host.authorize_window()?;
            let control = host.control();
            let sync = SyncHost::production(
                host.open_repository_connection()?,
                collagene_repository::http::HttpOptions {
                    allow_loopback_http,
                    ..Default::default()
                },
            )?;
            let mut source_scope = [0; 32];
            getrandom::fill(&mut source_scope)
                .map_err(|_| "Cannot initialize private source scope")?;
            settings::persist(&configuration, &root)?;
            app.manage(NativeState {
                host: Arc::new(Mutex::new(host)),
                sync: Arc::new(Mutex::new(sync)),
                control,
                access,
                local_active: Arc::new(AtomicUsize::new(0)),
                sync_active: Arc::new(AtomicUsize::new(0)),
                source_active: Arc::new(AtomicUsize::new(0)),
                sources: Arc::new(sources::Gateway::native(sources::Config::default())),
                source_scope,
            });
            tauri::WebviewWindowBuilder::new(
                app,
                "main",
                tauri::WebviewUrl::App("index.html".into()),
            )
            .title("Collagene")
            .inner_size(1440.0, 960.0)
            .min_inner_size(1024.0, 700.0)
            .use_https_scheme(true)
            .incognito(true)
            .devtools(false)
            .on_navigation(policy::trusted_url)
            .on_new_window(|_, _| tauri::webview::NewWindowResponse::Deny)
            .build()?;
            Ok(())
        })
        .run(tauri::generate_context!())?;
    Ok(())
}
fn run_mcp(options: settings::McpOptions) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let selections = options
        .documents
        .iter()
        .map(|value| collagene_desktop_mcp::NativeSelection::parse(value))
        .collect::<Result<Vec<_>, _>>()?;
    let session = match options.storage_root {
        Some(root) => collagene_desktop_mcp::create_with_storage(&root, &selections, Some(&root))?,
        None if selections.is_empty() => collagene_desktop_mcp::Session::empty(),
        None => return Err("Saved MCP selections require an explicit storage root".into()),
    };
    // No Tauri application/event loop, native window, vault or child process is
    // created in stdio mode. stdout belongs exclusively to the MCP transport.
    tauri::async_runtime::block_on(collagene_desktop_mcp::serve(session))
}
fn main() {
    let launch = match settings::command_line() {
        Ok(launch) => launch,
        Err(error) => {
            eprintln!("Invalid Collagene arguments: {error}");
            std::process::exit(2);
        }
    };
    match launch {
        settings::Launch::Help => println!("{}", settings::HELP),
        settings::Launch::Version => println!("Collagene {}", env!("CARGO_PKG_VERSION")),
        settings::Launch::Mcp(options) => {
            if let Err(error) = run_mcp(options) {
                eprintln!("Collagene MCP could not start: {error}");
                std::process::exit(1);
            }
        }
        settings::Launch::Gui(options) => {
            if let Err(error) = run(options) {
                eprintln!("Collagene desktop could not start: {error}");
                std::process::exit(1);
            }
        }
    }
}
