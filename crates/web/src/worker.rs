use collagene_core::protocol::{self, request::Request, response::Response};
use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::{JsCast, closure::Closure};

// The build pins the loader, JavaScript glue and WASM together. A fixed loader
// URL can otherwise combine cached files from different application releases.
fn worker_url() -> Result<String, String> {
    let missing = || "Analysis worker assets are missing from this build.".to_string();
    let document = web_sys::window()
        .and_then(|window| window.document())
        .ok_or_else(missing)?;
    let url = document
        .query_selector("meta[name=collagene-worker]")
        .map_err(|_| missing())?
        .and_then(|element| element.get_attribute("content"))
        .ok_or_else(missing)?;
    let bundle = url
        .rsplit_once("/workers/")
        .map(|(_, tail)| tail)
        .or_else(|| url.strip_prefix("workers/"));
    let valid = bundle.is_some_and(|tail| {
        tail.split_once('/').is_some_and(|(hash, file)| {
            hash.len() == 64
                && hash
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                && file == "analysis-worker_loader.js"
        })
    });
    if !valid || url.len() > 4096 || url.chars().any(char::is_control) {
        return Err(missing());
    }
    Ok(url)
}

/// One isolated worker per explicit job. Termination enforces the job timeout.
/// A retained worker pool and chunked document cache can replace this transport later.
pub async fn analyze(request: Request) -> Result<Response, String> {
    analyze_inner(request, None).await
}
/// Stops an isolated computation promptly when its owning view is disposed.
pub async fn analyze_cancellable(
    request: Request,
    alive: std::sync::Arc<std::sync::atomic::AtomicBool>,
) -> Result<Response, String> {
    analyze_inner(request, Some(alive)).await
}
async fn analyze_inner(
    request: Request,
    alive: Option<std::sync::Arc<std::sync::atomic::AtomicBool>>,
) -> Result<Response, String> {
    let bytes = protocol::request::Request::encode(&request).map_err(|e| e.to_string())?;
    let worker = web_sys::Worker::new(&worker_url()?)
        .map_err(|_| "Cannot start analysis worker".to_string())?;
    let (sender, receiver) = futures_channel::oneshot::channel();
    let sender = Rc::new(RefCell::new(Some(sender)));
    let on_message = {
        let sender = sender.clone();
        let worker = worker.clone();
        Closure::<dyn FnMut(web_sys::MessageEvent)>::new(move |event: web_sys::MessageEvent| {
            if let Some(control) = event.data().as_string() {
                if control == "ready" {
                    let array = js_sys::Uint8Array::from(bytes.as_slice());
                    let transfer = js_sys::Array::new();
                    transfer.push(&array.buffer());
                    if worker
                        .post_message_with_transfer(&array, &transfer)
                        .is_err()
                    {
                        if let Some(sender) = sender.borrow_mut().take() {
                            let _ = sender.send(Err("Cannot send analysis job".into()));
                        }
                    }
                } else if let Some(sender) = sender.borrow_mut().take() {
                    let _ = sender.send(Err(control));
                }
            } else if let Some(sender) = sender.borrow_mut().take() {
                let array = js_sys::Uint8Array::new(&event.data());
                let result = if array.length() as usize
                    > collagene_core::protocol::envelope::MAX_WORKSPACE_MESSAGE_BYTES + 42
                {
                    Err("Analysis result is too large".into())
                } else {
                    protocol::response::Response::decode(&array.to_vec()).map_err(|e| e.to_string())
                };
                let _ = sender.send(result);
            }
        })
    };
    let on_error = {
        let sender = sender.clone();
        Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            if let Some(sender) = sender.borrow_mut().take() {
                let _ = sender.send(Err("Analysis worker failed to load or run".into()));
            }
        })
    };
    worker.set_onmessage(Some(on_message.as_ref().unchecked_ref()));
    worker.set_onerror(Some(on_error.as_ref().unchecked_ref()));
    let timeout = {
        let sender = sender.clone();
        Closure::<dyn FnMut()>::new(move || {
            if let Some(sender) = sender.borrow_mut().take() {
                let _=sender.send(Err("Analysis exceeded its 30 second limit; use a smaller input or more specific query".into()));
            }
        })
    };
    let window = web_sys::window().ok_or("No browser window")?;
    let timer = window
        .set_timeout_with_callback_and_timeout_and_arguments_0(
            timeout.as_ref().unchecked_ref(),
            30_000,
        )
        .map_err(|_| "Cannot schedule analysis timeout")?;
    let cancellation = alive.map(|alive| {
        let sender = sender.clone();
        Closure::<dyn FnMut()>::new(move || {
            if !alive.load(std::sync::atomic::Ordering::SeqCst) {
                if let Some(sender) = sender.borrow_mut().take() {
                    let _ = sender.send(Err("Analysis was cancelled".into()));
                }
            }
        })
    });
    let cancellation_timer = cancellation.as_ref().and_then(|callback| {
        window
            .set_interval_with_callback_and_timeout_and_arguments_0(
                callback.as_ref().unchecked_ref(),
                50,
            )
            .ok()
    });
    let result = receiver
        .await
        .map_err(|_| "Analysis was interrupted".to_string());
    window.clear_timeout_with_handle(timer);
    if let Some(timer) = cancellation_timer {
        window.clear_interval_with_handle(timer);
    }
    worker.set_onmessage(None);
    worker.set_onerror(None);
    worker.terminate();
    let response = result??;
    if response.id != request.id || response.revision != request.revision {
        return Err("Analysis worker response does not match its request".into());
    }
    Ok(response)
}

/// A retained worker for one immutable sequence-view source. Source installation is
/// transferred once; viewport requests contain only its revision/options/window.
pub struct CachedWorker {
    worker: web_sys::Worker,
    state: Rc<RefCell<CachedState>>,
    _message: Closure<dyn FnMut(web_sys::MessageEvent)>,
    _error: Closure<dyn FnMut(web_sys::Event)>,
}
struct CachedState {
    ready: bool,
    closed: bool,
    pending: Option<futures_channel::oneshot::Sender<Result<Response, String>>>,
    queued: Option<Vec<u8>>,
}
fn post_bytes(worker: &web_sys::Worker, bytes: &[u8]) -> Result<(), String> {
    let array = js_sys::Uint8Array::from(bytes);
    let transfer = js_sys::Array::new();
    transfer.push(&array.buffer());
    worker
        .post_message_with_transfer(&array, &transfer)
        .map_err(|_| "Cannot send cached analysis job".into())
}
impl CachedWorker {
    pub fn new() -> Result<Self, String> {
        let worker = web_sys::Worker::new(&worker_url()?)
            .map_err(|_| "Cannot start sequence track worker".to_string())?;
        let state = Rc::new(RefCell::new(CachedState {
            ready: false,
            closed: false,
            pending: None,
            queued: None,
        }));
        let on_message = {
            let state = state.clone();
            let worker = worker.clone();
            Closure::<dyn FnMut(web_sys::MessageEvent)>::new(move |event: web_sys::MessageEvent| {
                if let Some(control) = event.data().as_string() {
                    if control == "ready" {
                        let queued = {
                            let mut s = state.borrow_mut();
                            s.ready = true;
                            s.queued.take()
                        };
                        if let Some(bytes) = queued {
                            if let Err(error) = post_bytes(&worker, &bytes) {
                                if let Some(tx) = state.borrow_mut().pending.take() {
                                    let _ = tx.send(Err(error));
                                }
                            }
                        }
                        return;
                    }
                    if let Some(tx) = state.borrow_mut().pending.take() {
                        let _ = tx.send(Err(control));
                    }
                } else {
                    let array = js_sys::Uint8Array::new(&event.data());
                    let result = if array.length() as usize
                        > collagene_core::formats::MAX_IMPORT_BYTES + 42
                    {
                        Err("Track report is too large".into())
                    } else {
                        protocol::response::Response::decode(&array.to_vec())
                            .map_err(|e| e.to_string())
                    };
                    if let Some(tx) = state.borrow_mut().pending.take() {
                        let _ = tx.send(result);
                    }
                }
            })
        };
        let on_error = {
            let state = state.clone();
            Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
                let mut s = state.borrow_mut();
                s.closed = true;
                s.queued = None;
                if let Some(tx) = s.pending.take() {
                    let _ = tx.send(Err("Sequence track worker failed".into()));
                }
            })
        };
        worker.set_onmessage(Some(on_message.as_ref().unchecked_ref()));
        worker.set_onerror(Some(on_error.as_ref().unchecked_ref()));
        Ok(Self {
            worker,
            state,
            _message: on_message,
            _error: on_error,
        })
    }
    pub fn close(&self) {
        self.worker.terminate();
        let mut state = self.state.borrow_mut();
        state.closed = true;
        state.queued = None;
        if let Some(sender) = state.pending.take() {
            let _ = sender.send(Err("Sequence track worker replaced or closed".into()));
        }
    }
    pub async fn ask(&self, request: Request) -> Result<Response, String> {
        let id = request.id;
        let revision = request.revision.clone();
        let bytes = protocol::request::Request::encode(&request).map_err(|e| e.to_string())?;
        let (tx, rx) = futures_channel::oneshot::channel();
        let ready = {
            let mut state = self.state.borrow_mut();
            if state.closed {
                return Err(
                    "Track worker is closed; change tracks or reopen the document to retry.".into(),
                );
            }
            if state.pending.is_some() {
                return Err("Track worker already has a request in flight".into());
            }
            state.pending = Some(tx);
            state.ready
        };
        if ready {
            if let Err(error) = post_bytes(&self.worker, &bytes) {
                self.state.borrow_mut().pending = None;
                return Err(error);
            }
        } else {
            self.state.borrow_mut().queued = Some(bytes);
        }
        let timeout = {
            let state = self.state.clone();
            let worker = self.worker.clone();
            Closure::<dyn FnMut()>::new(move || {
                let mut s = state.borrow_mut();
                s.closed = true;
                s.queued = None;
                worker.terminate();
                if let Some(tx) = s.pending.take() {
                    let _ = tx.send(Err(
                        "Sequence tracks exceeded the 30 second job limit".into()
                    ));
                }
            })
        };
        let window = web_sys::window().ok_or("No browser window")?;
        let timer = window
            .set_timeout_with_callback_and_timeout_and_arguments_0(
                timeout.as_ref().unchecked_ref(),
                30_000,
            )
            .map_err(|_| "Cannot schedule track job timeout")?;
        let result = rx.await.map_err(|_| "Track job interrupted".to_string());
        window.clear_timeout_with_handle(timer);
        let response = result??;
        if response.id != id || response.revision != revision {
            return Err("Track worker response does not match its request".into());
        }
        Ok(response)
    }
}
impl Drop for CachedWorker {
    fn drop(&mut self) {
        self.worker.set_onmessage(None);
        self.worker.set_onerror(None);
        self.worker.terminate();
        let mut s = self.state.borrow_mut();
        s.closed = true;
        s.pending = None;
        s.queued = None;
    }
}
