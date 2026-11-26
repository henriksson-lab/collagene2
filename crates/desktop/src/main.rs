#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
#[cfg(any(target_os = "android", target_os = "ios"))]
compile_error!("Collagene desktop supports Windows, macOS and Linux; mobile is deferred.");
mod policy;
mod settings;
use collagene_desktop_host::{
    DesktopHost, HostControl, MAX_REQUEST_BYTES, MAX_RESPONSE_BYTES, WindowAccess,
};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicUsize, Ordering},
};
use tauri::Manager;

#[derive(Clone)]
struct NativeState {
    host: Arc<Mutex<DesktopHost>>,
    control: HostControl,
    access: WindowAccess,
    active: Arc<AtomicUsize>,
}
struct JobLease(Arc<AtomicUsize>);
impl JobLease {
    fn acquire(active: &Arc<AtomicUsize>) -> Result<Self, String> {
        // Repository mutations serialize; control/cancellation bypasses this
        // admission gate and never waits on the repository mutex.
        active
            .compare_exchange(0, 1, Ordering::AcqRel, Ordering::Acquire)
            .map_err(|_| "Native repository is busy; retry the same operation ID.".to_string())?;
        Ok(Self(active.clone()))
    }
}
impl Drop for JobLease {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::Release);
    }
}
fn bounded_response(bytes: Vec<u8>) -> Result<tauri::ipc::Response, String> {
    if bytes.len() > MAX_RESPONSE_BYTES {
        return Err("Native reply exceeded its binary protocol limit.".into());
    }
    Ok(tauri::ipc::Response::new(bytes))
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
    if bytes.len() > MAX_REQUEST_BYTES {
        return Err("Native request exceeds its binary protocol limit.".into());
    }
    // The host recognizes only bounded control messages here; archive/native
    // decoding and all repository work happen on the blocking worker below.
    if let Some(response) = state
        .control
        .dispatch_raw(&state.access, bytes)
        .map_err(|_| "Native control request was rejected.".to_string())?
    {
        return bounded_response(response);
    }
    let lease = JobLease::acquire(&state.active)?;
    let bytes = bytes.clone();
    let state = state.inner().clone();
    let response = tauri::async_runtime::spawn_blocking(move || {
        // A dropped IPC promise does not release the job slot while its native
        // transaction is still running. Durable IDs/status resolve lost replies.
        let _lease = lease;
        let mut host = state
            .host
            .lock()
            .map_err(|_| "Native repository is unavailable; retain local work.".to_string())?;
        host.dispatch_raw(&state.access, &bytes).map_err(|_| {
            "Native request failed; retain its ID and check status before retrying.".to_string()
        })
    })
    .await
    .map_err(|_| {
        "Native request task could not complete; check its durable status.".to_string()
    })??;
    bounded_response(response)
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let storage_override = settings::command_line()?;
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
            settings::persist(&configuration, &root)?;
            app.manage(NativeState {
                host: Arc::new(Mutex::new(host)),
                control,
                access,
                active: Arc::new(AtomicUsize::new(0)),
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
fn main() {
    if std::env::args_os().skip(1).collect::<Vec<_>>().as_slice()
        == [std::ffi::OsString::from("--help")]
    {
        println!(
            "Usage: collagene-desktop [--storage-root DIRECTORY]\nThe native repository root is persisted privately; frontend requests cannot change it."
        );
        return;
    }
    if let Err(error) = run() {
        eprintln!("Collagene desktop could not start: {error}");
        std::process::exit(1);
    }
}
