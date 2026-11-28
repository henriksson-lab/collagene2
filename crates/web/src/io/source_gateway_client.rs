//! Per-dialog native/hosted source-import client. Catalog responses contain only
//! native handles and choices; no provider download URL crosses this boundary.
use collagene_source_gateway::wire;
use std::{
    cell::{Cell, RefCell},
    collections::BTreeMap,
    rc::Rc,
};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::spawn_local;

pub fn enabled() -> bool {
    crate::platform::desktop_rpc::is_desktop()
        || web_sys::window()
            .and_then(|w| w.document())
            .and_then(|d| {
                d.query_selector("meta[name='collagene-source-gateway'][content='enabled']")
                    .ok()
                    .flatten()
            })
            .is_some()
}
struct State {
    nonce: [u8; 32],
    next: Cell<u64>,
    closed: Cell<bool>,
    active: RefCell<BTreeMap<u64, web_sys::AbortController>>,
    catalogs: RefCell<Vec<[u8; 32]>>,
}
#[derive(Clone)]
pub struct Client(Rc<State>);
#[derive(Clone, Copy)]
enum Expected {
    Import,
    Catalog(u64),
    Cancel,
    Dispose,
}
impl Client {
    pub fn new() -> Result<Self, String> {
        if !enabled() {
            return Err("Source gateway is unavailable in this host.".into());
        }
        let random = js_sys::Uint8Array::new_with_length(32);
        let window = web_sys::window().ok_or("No application window")?;
        let crypto = js_sys::Reflect::get(&window, &"crypto".into())
            .map_err(|_| "Secure source session unavailable")?;
        let get = js_sys::Reflect::get(&crypto, &"getRandomValues".into())
            .map_err(|_| "Secure source session unavailable")?
            .dyn_into::<js_sys::Function>()
            .map_err(|_| "Secure source session unavailable")?;
        get.call1(&crypto, &random)
            .map_err(|_| "Cannot create a private source session")?;
        let mut nonce = [0; 32];
        random.copy_to(&mut nonce);
        Ok(Self(Rc::new(State {
            nonce,
            next: Cell::new(0),
            closed: Cell::new(false),
            active: RefCell::new(BTreeMap::new()),
            catalogs: RefCell::new(vec![]),
        })))
    }
    fn next(&self) -> Result<u64, String> {
        let next = self
            .0
            .next
            .get()
            .checked_add(1)
            .ok_or("Source session IDs exhausted; reopen this dialog.")?;
        self.0.next.set(next);
        Ok(next)
    }
    pub async fn request(&self, operation: wire::Operation) -> Result<wire::Output, String> {
        let expected = match &operation {
            wire::Operation::Import { .. } => Expected::Import,
            wire::Operation::Catalog { source, .. } => {
                match collagene_sources::parse_source(source) {
                    Ok(collagene_sources::Source::Addgene { plasmid_id }) => {
                        Expected::Catalog(plasmid_id)
                    }
                    _ => {
                        return Err(
                            "Choose an Addgene plasmid before loading its sequence catalog.".into(),
                        );
                    }
                }
            }
            wire::Operation::Select { .. } => Expected::Import,
            wire::Operation::Cancel { .. } => Expected::Cancel,
            wire::Operation::DisposeCatalog { .. } => Expected::Dispose,
        };
        let control = matches!(expected, Expected::Cancel | Expected::Dispose);
        if self.0.closed.get() && !control {
            return Err("Source dialog closed.".into());
        }
        let id = self.next()?;
        let request = wire::Request {
            id,
            client: self.0.nonce,
            operation,
        };
        let bytes = wire::encode_request(&request)
            .map_err(|_| "Source request exceeds its supported bound.")?;
        drop(request); // release the provider-token request before awaiting I/O
        let controller =
            web_sys::AbortController::new().map_err(|_| "Cannot prepare source cancellation")?;
        if !control {
            self.0.active.borrow_mut().insert(id, controller.clone());
        }
        let _guard = Active {
            state: self.0.clone(),
            id,
        };
        let bytes = if crate::platform::desktop_rpc::is_desktop() {
            crate::platform::desktop_rpc::source_raw(
                bytes,
                if control {
                    None
                } else {
                    Some(controller.signal())
                },
            )
            .await
            .map_err(|_| "Native source request interrupted or cancelled.".to_owned())?
        } else {
            crate::library::library_client::browser_session()
                .source_request(&bytes, &controller.signal())
                .await?
        };
        if bytes.len() > wire::MAX_RESPONSE {
            return Err("Source gateway response exceeds its size bound.".into());
        }
        let response =
            wire::decode_response(&bytes).map_err(|_| "Invalid source gateway response.")?;
        if response.id != id {
            return Err("Source gateway response belongs to another request.".into());
        }
        let output = response.result.map_err(|f| safe_failure(&f))?;
        let matches = matches!(
            (expected, &output),
            (Expected::Import, wire::Output::Imported(_))
                | (Expected::Cancel, wire::Output::Canceled { .. })
                | (Expected::Dispose, wire::Output::CatalogDisposed { .. })
        ) || matches!((expected,&output),(Expected::Catalog(id),wire::Output::Catalog(c)) if id==c.plasmid_id);
        if !matches {
            return Err("Source gateway returned a different operation result.".into());
        }
        if let wire::Output::Catalog(c) = &output {
            validate_catalog(c)?;
            self.0.catalogs.borrow_mut().push(c.handle);
        }
        if !control && (self.0.closed.get() || controller.signal().aborted()) {
            if let wire::Output::Catalog(c) = &output {
                self.dispose_catalog(c.handle);
            }
            return Err("Source response discarded after cancellation.".into());
        }
        Ok(output)
    }
    pub fn cancel_active(&self) {
        let active = std::mem::take(&mut *self.0.active.borrow_mut());
        for (request_id, controller) in active {
            controller.abort();
            let client = self.clone();
            spawn_local(async move {
                let _ = client.request(wire::Operation::Cancel { request_id }).await;
            });
        }
    }
    pub fn dispose_catalog(&self, handle: [u8; 32]) {
        self.0.catalogs.borrow_mut().retain(|h| h != &handle);
        let client = self.clone();
        spawn_local(async move {
            let _ = client
                .request(wire::Operation::DisposeCatalog { handle })
                .await;
        });
    }
    pub fn reset(&self) {
        self.cancel_active();
        let handles = std::mem::take(&mut *self.0.catalogs.borrow_mut());
        for handle in handles {
            self.dispose_catalog(handle);
        }
    }
    pub fn close(&self) {
        self.0.closed.set(true);
        self.reset();
    }
}
struct Active {
    state: Rc<State>,
    id: u64,
}
impl Drop for Active {
    fn drop(&mut self) {
        self.state.active.borrow_mut().remove(&self.id);
    }
}
fn validate_catalog(c: &wire::Catalog) -> Result<(), String> {
    let mut ids = std::collections::BTreeSet::new();
    if c.plasmid_id == 0
        || c.name.len() > 8192
        || c.choices.len() > 1000
        || c.response_blake3.len() != 64
        || !c
            .response_blake3
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        || c.expires_in_seconds == 0
        || c.choices
            .iter()
            .any(|v| v.sequence_id == 0 || !ids.insert(v.sequence_id) || v.description.len() > 8192)
    {
        return Err("Invalid bounded source catalog.".into());
    }
    Ok(())
}
fn safe_failure(f: &wire::Failure) -> String {
    // Never surface provider/HTTP error strings that could contain a signed URL,
    // token or internal endpoint. The typed result is enough to guide a retry.
    use wire::ErrorCode::*;
    match f.code {
        InvalidRequest=>"The source request is invalid. Check the accession or URL.",
        Forbidden=>"This source request is not permitted. Check the lab login or provider access.",
        Busy=>"Source import is busy. Wait briefly and retry.",Timeout=>"Source retrieval timed out. Retry or import a downloaded sequence file.",
        Canceled=>"Source retrieval cancelled.",ExpiredCatalog=>"The sequence catalog expired. Refresh it and choose a sequence again.",
        Provider=>"The provider could not supply this record. Check its identifier and authorized access, or import a downloaded sequence file.",
        TooLarge=>"The source exceeds its supported download or catalog size.",Internal=>"The source gateway could not complete this request.",
    }.into()
}
