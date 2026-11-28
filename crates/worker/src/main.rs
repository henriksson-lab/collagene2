#[cfg(target_arch = "wasm32")]
fn main() {
    use wasm_bindgen::{JsCast, closure::Closure};
    let scope: web_sys::DedicatedWorkerGlobalScope = js_sys::global().unchecked_into();
    let sender = scope.clone();
    let mut session = collagene_core::protocol::session::AnalysisSession::default();
    let handler =
        Closure::<dyn FnMut(web_sys::MessageEvent)>::new(move |event: web_sys::MessageEvent| {
            let array = js_sys::Uint8Array::new(&event.data());
            if array.length() as usize
                > collagene_core::protocol::envelope::MAX_WORKSPACE_MESSAGE_BYTES + 42
            {
                let _ = sender.post_message(&"Worker input is too large".into());
                return;
            }
            match session.process_request(&array.to_vec()) {
                Ok(bytes) => {
                    let response = js_sys::Uint8Array::from(bytes.as_slice());
                    let transfer = js_sys::Array::new();
                    transfer.push(&response.buffer());
                    let _ = sender.post_message_with_transfer(&response, &transfer);
                }
                Err(error) => {
                    let _ = sender.post_message(&error.to_string().into());
                }
            }
        });
    scope.set_onmessage(Some(handler.as_ref().unchecked_ref()));
    handler.forget();
    let _ = scope.post_message(&"ready".into());
}

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    println!("This binary runs inside a browser Web Worker.");
}
