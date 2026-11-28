use collagene_core::verification;

use leptos::prelude::*;
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::spawn_local;

thread_local! { static EXPORT_NOTICE: std::cell::RefCell<Option<RwSignal<String>>> = const { std::cell::RefCell::new(None) }; }
/// An owned, bounded export snapshot. Native admission precedes its byte copy.
pub(crate) struct ExportArtifact {
    name: String,
    mime: String,
    bytes: Vec<u8>,
    reservation: Option<crate::platform::desktop_rpc::ExportReservation>,
    notice: Option<RwSignal<String>>,
}
pub(crate) fn prepare_download(
    name: &str,
    mime: &str,
    bytes: &[u8],
) -> Result<ExportArtifact, String> {
    if bytes.len() > 64 * 1024 * 1024 {
        return Err("Export exceeds the 64 MiB file limit.".into());
    }
    let reservation = if crate::platform::desktop_rpc::is_desktop() {
        Some(crate::platform::desktop_rpc::reserve_export()?)
    } else {
        None
    };
    let name = if reservation.is_some() {
        let mut safe: String = name
            .chars()
            .map(|c| {
                if c.is_control() || matches!(c, '/' | '\\' | ':') {
                    '_'
                } else {
                    c
                }
            })
            .collect();
        while safe.len() > 256 {
            safe.pop();
        }
        if safe.trim().is_empty() {
            safe = "export.bin".into();
        }
        safe
    } else {
        name.into()
    };
    Ok(ExportArtifact {
        name,
        mime: mime.into(),
        bytes: bytes.to_vec(),
        reservation,
        notice: EXPORT_NOTICE.with(|notice| *notice.borrow()),
    })
}
impl ExportArtifact {
    pub async fn save(self, guard: impl Fn() -> Result<(), String>) -> Result<(), String> {
        if let Some(reservation) = self.reservation {
            let name = crate::platform::desktop_rpc::export_reserved(
                reservation,
                self.name,
                self.bytes,
                guard,
            )
            .await?;
            if let Some(notice) = self.notice {
                let _ = notice.try_set(format!("Saved export: exports/{name}"));
            }
            Ok(())
        } else {
            guard()?;
            browser_download(&self.name, &self.mime, &self.bytes)
        }
    }
    pub fn complete(self, done: impl FnOnce(Result<(), String>) + 'static) {
        if self.reservation.is_none() {
            done(browser_download(&self.name, &self.mime, &self.bytes));
        } else {
            wasm_bindgen_futures::spawn_local(async move { done(self.save(|| Ok(())).await) });
        }
    }
    pub fn complete_guarded(
        self,
        guard: impl Fn() -> Result<(), String> + 'static,
        done: impl FnOnce(Result<(), String>) + 'static,
    ) {
        if self.reservation.is_none() {
            done(guard().and_then(|()| browser_download(&self.name, &self.mime, &self.bytes)));
        } else {
            wasm_bindgen_futures::spawn_local(async move { done(self.save(guard).await) });
        }
    }
}
pub(crate) async fn download(name: &str, mime: &str, bytes: &[u8]) -> Result<(), String> {
    prepare_download(name, mime, bytes)?.save(|| Ok(())).await
}
pub(crate) async fn download_guarded(
    name: &str,
    mime: &str,
    bytes: &[u8],
    guard: impl Fn() -> Result<(), String>,
) -> Result<(), String> {
    prepare_download(name, mime, bytes)?.save(guard).await
}
/// Completion runs after the native durable receipt or browser dispatch, never
/// merely after queueing. Capturing owned bytes keeps later edits out of the file.
pub(crate) fn download_with_completion(
    name: &str,
    mime: &str,
    bytes: &[u8],
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    match prepare_download(name, mime, bytes) {
        Ok(artifact) => artifact.complete(done),
        Err(error) => done(Err(error)),
    }
}
pub(crate) fn download_with_completion_guarded(
    name: &str,
    mime: &str,
    bytes: &[u8],
    guard: impl Fn() -> Result<(), String> + 'static,
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    complete_download_guarded(prepare_download(name, mime, bytes), guard, done);
}
pub(crate) fn complete_download(
    result: Result<ExportArtifact, String>,
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    match result {
        Ok(artifact) => artifact.complete(done),
        Err(error) => done(Err(error)),
    }
}
pub(crate) fn complete_download_guarded(
    result: Result<ExportArtifact, String>,
    guard: impl Fn() -> Result<(), String> + 'static,
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    match result {
        Ok(artifact) => artifact.complete_guarded(guard, done),
        Err(error) => done(Err(error)),
    }
}

fn browser_download(name: &str, mime: &str, bytes: &[u8]) -> Result<(), String> {
    let array = js_sys::Uint8Array::from(bytes);
    let parts = js_sys::Array::new();
    parts.push(&array);
    let options = web_sys::BlobPropertyBag::new();
    options.set_type(mime);
    let blob = web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &options)
        .map_err(|_| "Cannot prepare download")?;
    let url =
        web_sys::Url::create_object_url_with_blob(&blob).map_err(|_| "Cannot create download")?;
    let window = web_sys::window().ok_or("No browser window")?;
    let document = window.document().ok_or("No browser document")?;
    let anchor = document
        .create_element("a")
        .map_err(|_| "Cannot create download link")?;
    anchor
        .set_attribute("href", &url)
        .map_err(|_| "Cannot set download URL")?;
    anchor
        .set_attribute("download", name)
        .map_err(|_| "Cannot set filename")?;
    let anchor = anchor.unchecked_into::<web_sys::HtmlElement>();
    if let Some(body) = document.body() {
        let _ = body.append_child(&anchor);
    }
    anchor.click();
    anchor.remove();
    // The click queues navigation; defer revocation until the browser consumes the URL.
    let cleanup = wasm_bindgen::closure::Closure::once_into_js(move || {
        let _ = web_sys::Url::revoke_object_url(&url);
    });
    let _ =
        window.set_timeout_with_callback_and_timeout_and_arguments_0(cleanup.unchecked_ref(), 1000);
    Ok(())
}

/// Owns the notification lifetime for bounded browser/native file exports.
#[derive(Clone, Copy)]
pub(crate) struct ExportService {
    pub(crate) notice: RwSignal<String>,
}
impl ExportService {
    pub(crate) fn install() -> Self {
        let native_export_notice = RwSignal::new(String::new());
        EXPORT_NOTICE.with(|notice| *notice.borrow_mut() = Some(native_export_notice));
        on_cleanup(move || {
            EXPORT_NOTICE.with(|notice| {
                if *notice.borrow() == Some(native_export_notice) {
                    *notice.borrow_mut() = None;
                }
            })
        });
        Self {
            notice: native_export_notice,
        }
    }
}
