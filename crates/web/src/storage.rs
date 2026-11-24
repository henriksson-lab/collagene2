use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::{JsCast, JsValue, closure::Closure};
use web_sys::{Event, IdbDatabase, IdbRequest, IdbTransaction, IdbTransactionMode};

async fn request(request: &IdbRequest) -> Result<JsValue, String> {
    let (tx, rx) = futures_channel::oneshot::channel();
    let tx = Rc::new(RefCell::new(Some(tx)));
    let success = {
        let tx = tx.clone();
        let request = request.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(tx) = tx.borrow_mut().take() {
                let _ = tx.send(
                    request
                        .result()
                        .map_err(|_| "Browser database request failed".into()),
                );
            }
        })
    };
    let error = {
        let tx = tx.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(tx) = tx.borrow_mut().take() {
                let _ = tx.send(Err("Browser database request failed".into()));
            }
        })
    };
    request.set_onsuccess(Some(success.as_ref().unchecked_ref()));
    request.set_onerror(Some(error.as_ref().unchecked_ref()));
    let result = rx
        .await
        .map_err(|_| "Browser database request was interrupted".to_string());
    request.set_onsuccess(None);
    request.set_onerror(None);
    result?
}
pub async fn open() -> Result<IdbDatabase, String> {
    let factory = web_sys::window()
        .ok_or("No browser window")?
        .indexed_db()
        .map_err(|_| "IndexedDB is unavailable")?
        .ok_or("IndexedDB is unavailable")?;
    let open = factory
        .open_with_u32("collagene", 1)
        .map_err(|_| "Cannot open browser database")?;
    let upgrade = {
        let open = open.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Ok(db) = open.result().and_then(|v| v.dyn_into::<IdbDatabase>()) {
                let _ = db.create_object_store("workspace");
            }
        })
    };
    open.set_onupgradeneeded(Some(upgrade.as_ref().unchecked_ref()));
    let result = request(open.unchecked_ref()).await;
    open.set_onupgradeneeded(None);
    result?
        .dyn_into::<IdbDatabase>()
        .map_err(|_| "Invalid database connection".into())
}
pub async fn load(db: &IdbDatabase) -> Result<Option<Vec<u8>>, String> {
    load_key(db, "current").await
}
pub async fn load_previous(db: &IdbDatabase) -> Result<Option<Vec<u8>>, String> {
    load_key(db, "previous").await
}
async fn load_key(db: &IdbDatabase, key: &str) -> Result<Option<Vec<u8>>, String> {
    let tx = db
        .transaction_with_str("workspace")
        .map_err(|_| "Cannot read autosave")?;
    let pending = tx
        .object_store("workspace")
        .map_err(|_| "Cannot read autosave store")?
        .get(&JsValue::from_str(key))
        .map_err(|_| "Cannot read autosave record")?;
    let value = request(&pending).await?;
    if value.is_undefined() {
        return Ok(None);
    }
    Ok(Some(js_sys::Uint8Array::new(&value).to_vec()))
}
async fn committed(tx: &IdbTransaction) -> Result<(), String> {
    let (sender, receiver) = futures_channel::oneshot::channel();
    let sender = Rc::new(RefCell::new(Some(sender)));
    let complete = {
        let sender = sender.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(sender) = sender.borrow_mut().take() {
                let _ = sender.send(Ok(()));
            }
        })
    };
    let abort = {
        let sender = sender.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(sender) = sender.borrow_mut().take() {
                let _ = sender.send(Err(
                    "Autosave transaction failed; export a file to keep your work".into(),
                ));
            }
        })
    };
    tx.set_oncomplete(Some(complete.as_ref().unchecked_ref()));
    tx.set_onabort(Some(abort.as_ref().unchecked_ref()));
    tx.set_onerror(Some(abort.as_ref().unchecked_ref()));
    let result = receiver
        .await
        .map_err(|_| "Autosave interrupted".to_string());
    tx.set_oncomplete(None);
    tx.set_onabort(None);
    tx.set_onerror(None);
    result?
}
pub async fn save(db: &IdbDatabase, bytes: &[u8]) -> Result<(), String> {
    let tx = db
        .transaction_with_str_and_mode("workspace", IdbTransactionMode::Readwrite)
        .map_err(|_| "Cannot start autosave")?;
    let value = js_sys::Uint8Array::from(bytes);
    tx.object_store("workspace")
        .map_err(|_| "Cannot access autosave store")?
        .put_with_key(&value, &JsValue::from_str("current"))
        .map_err(|_| "Cannot write autosave")?;
    committed(&tx).await
}

/// Replacing a workspace and retaining its recovery copy is one transaction.
/// Quota/abort failure leaves both prior slots intact.
pub async fn replace_workspace(
    db: &IdbDatabase,
    current: &[u8],
    previous: &[u8],
) -> Result<(), String> {
    replace_workspace_cancellable(db, current, previous, Rc::new(RefCell::new(None))).await
}

/// A caller-owned handle permits aborting an uncommitted atomic replacement on
/// cancellation/unmount. Successful completion is the commit point.
pub async fn replace_workspace_cancellable(
    db: &IdbDatabase,
    current: &[u8],
    previous: &[u8],
    pending: Rc<RefCell<Option<IdbTransaction>>>,
) -> Result<(), String> {
    let tx = db
        .transaction_with_str_and_mode("workspace", IdbTransactionMode::Readwrite)
        .map_err(|_| "Cannot start workspace recovery transaction")?;
    let store = tx
        .object_store("workspace")
        .map_err(|_| "Cannot access workspace recovery store")?;
    *pending.borrow_mut() = Some(tx.clone());
    // A synchronous second-request failure must also roll back the first put.
    // Dropping an IndexedDB transaction does not abort its queued requests.
    let queued = (|| {
        store
            .put_with_key(
                &js_sys::Uint8Array::from(previous),
                &JsValue::from_str("previous"),
            )
            .map_err(|_| "Cannot retain previous workspace")?;
        store
            .put_with_key(
                &js_sys::Uint8Array::from(current),
                &JsValue::from_str("current"),
            )
            .map_err(|_| "Cannot restore workspace")?;
        Ok::<_, &str>(())
    })();
    if let Err(error) = queued {
        let _ = tx.abort();
        pending.borrow_mut().take();
        return Err(error.into());
    }
    let result = committed(&tx).await;
    pending.borrow_mut().take();
    result
}
