//! URL import selects the native/authenticated host gateway when available.
//! The direct CORS importer remains only for the static preview.
use collagene_core::formats::sequence::{ImportReport, MAX_IMPORT_BYTES};
use collagene_sources::{
    AddgeneCatalog, Config, Source, SourceError, addgene_catalog_url, import_download,
    parse_addgene_catalog, parse_source,
};
use leptos::prelude::*;
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::{JsFuture, spawn_local};

async fn fetch_bounded(
    url: &str,
    signal: &web_sys::AbortSignal,
    addgene_token: Option<&str>,
) -> Result<Vec<u8>, String> {
    let init = web_sys::RequestInit::new();
    init.set_method("GET");
    init.set_mode(web_sys::RequestMode::Cors);
    init.set_credentials(web_sys::RequestCredentials::Omit);
    init.set_signal(Some(signal));
    if addgene_token.is_some() {
        let parsed = web_sys::Url::new(url).map_err(|_| "Invalid provider request")?;
        if parsed.origin() != "https://api.developers.addgene.org" {
            return Err("Addgene credentials cannot be sent to another origin.".into());
        }
        // Browser APIs cannot safely inspect cross-origin redirect locations. Native retrieval
        // handles the official redirect route; browser credential-bearing requests reject redirects.
        init.set_redirect(web_sys::RequestRedirect::Error);
    }
    let request = web_sys::Request::new_with_str_and_init(url, &init)
        .map_err(|_| "Could not prepare source request")?;
    if let Some(token) = addgene_token {
        request
            .headers()
            .set("Authorization", &format!("Token {token}"))
            .map_err(|_| "Cannot prepare provider authorization")?;
    }
    let window = web_sys::window().ok_or("No browser window")?;
    let response=JsFuture::from(window.fetch_with_request(&request)).await.map_err(|_|"Download failed or was blocked by authentication/CORS. Use an authorized native/MCP import or open the downloaded file.")?.dyn_into::<web_sys::Response>().map_err(|_|"Invalid source response")?;
    if !response.ok() {
        return Err(format!(
            "Sequence source returned HTTP {}. Expired download links need to be refreshed.",
            response.status()
        ));
    }
    if response
        .headers()
        .get("content-type")
        .ok()
        .flatten()
        .is_some_and(|s| s.to_ascii_lowercase().contains("text/html"))
    {
        return Err("The URL returned an HTML/login page, not a sequence file.".into());
    }
    if response
        .headers()
        .get("content-length")
        .ok()
        .flatten()
        .and_then(|v| v.parse::<u64>().ok())
        .is_some_and(|n| n > MAX_IMPORT_BYTES as u64)
    {
        return Err("Source exceeds the 32 MiB download limit.".into());
    }
    let stream = response
        .body()
        .ok_or("Source response has no readable body")?;
    let reader = stream
        .get_reader()
        .dyn_into::<web_sys::ReadableStreamDefaultReader>()
        .map_err(|_| "Cannot read source stream")?;
    let mut bytes = Vec::new();
    loop {
        let result = JsFuture::from(reader.read())
            .await
            .map_err(|_| "Source stream was interrupted; no documents were imported")?;
        if js_sys::Reflect::get(&result, &"done".into())
            .ok()
            .and_then(|v| v.as_bool())
            .unwrap_or(false)
        {
            break;
        }
        let value =
            js_sys::Reflect::get(&result, &"value".into()).map_err(|_| "Invalid source chunk")?;
        let chunk = js_sys::Uint8Array::new(&value);
        if chunk.length() as usize > MAX_IMPORT_BYTES.saturating_sub(bytes.len()) {
            let _ = reader.cancel();
            return Err("Source exceeds the 32 MiB download limit.".into());
        }
        bytes.extend(chunk.to_vec());
    }
    reader.release_lock();
    if bytes.is_empty() {
        return Err("Sequence source returned an empty file.".into());
    }
    Ok(bytes)
}

fn provider_description(input: &str) -> String {
    match parse_source(input) {
        Ok(Source::Ncbi { accession }) => format!("NCBI nuccore · {accession} · annotated GenBank record. The returned accession/version must match the request."),
        Ok(Source::Ena { accession }) => format!("ENA · {accession} · annotated EMBL record. The returned accession/version must match the request."),
        Ok(Source::Ddbj { accession }) => format!("DDBJ · {accession} · annotated GenBank record. The returned accession/version must match the request."),
        Ok(Source::Addgene { plasmid_id }) => format!("Addgene plasmid {plasmid_id} · approved API access and an explicit sequence-record choice are required."),
        Ok(Source::Direct { .. }) => "Direct HTTPS sequence file · provider accession checks require a recognized database URL or explicit provider identifier.".into(),
        Err(_) => "Use ena: or ddbj: to choose those providers. Bare nucleotide accessions use NCBI; bare plasmid numbers use Addgene.".into(),
    }
}

#[component]
pub fn DirectSourceImport(on_import: Callback<ImportReport, Result<(), String>>) -> impl IntoView {
    let input = RwSignal::new(String::new());
    let token = RwSignal::new(String::new());
    let catalog = RwSignal::new(None::<AddgeneCatalog>);
    let selected = RwSignal::new(String::new());
    let generation = RwSignal::new(0u64);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let alive = Rc::new(Cell::new(true));
    let controller = Rc::new(RefCell::new(None::<web_sys::AbortController>));
    let last_request = Rc::new(Cell::new(0.0f64));
    let input_reset = controller.clone();
    let invalidate = StoredValue::new_local(move || {
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(abort) = input_reset.borrow_mut().take() {
            abort.abort();
        }
        catalog.set(None);
        selected.set(String::new());
        busy.set(false);
        status.set(String::new());
    });
    let selection_reset = StoredValue::new_local(controller.clone());
    let cleanup_state = StoredValue::new_local((alive.clone(), controller.clone()));
    on_cleanup(move || {
        cleanup_state.with_value(|(alive, controller)| {
            alive.set(false);
            if let Some(controller) = controller.borrow_mut().take() {
                controller.abort();
            }
        })
    });
    let import = move |_| {
        if busy.get_untracked() {
            return;
        }
        let source = match parse_source(&input.get_untracked()) {
            Ok(s) => s,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        enum Mode {
            Direct(Source, String),
            Catalog(u64, String),
            Selection(AddgeneCatalog, u64, String),
        }
        let password = token.get_untracked();
        let mode = match source {
            Source::Addgene { plasmid_id } => {
                if password.is_empty() {
                    status.set(SourceError::AddgeneAccessRequired { plasmid_id }.to_string());
                    return;
                }
                if let Err(e) = Config::default().with_addgene_api_token(&password) {
                    status.set(e.to_string());
                    return;
                }
                if let Some(catalog) = catalog.get_untracked() {
                    let id = match selected.get_untracked().parse::<u64>() {
                        Ok(id) => id,
                        Err(_) => {
                            status
                                .set("Select an Addgene sequence record before importing.".into());
                            return;
                        }
                    };
                    match catalog.download_url(id) {
                        Ok(url) => Mode::Selection(catalog, id, url),
                        Err(e) => {
                            status.set(e.to_string());
                            return;
                        }
                    }
                } else {
                    match addgene_catalog_url(plasmid_id) {
                        Ok(url) => Mode::Catalog(plasmid_id, url),
                        Err(e) => {
                            status.set(e.to_string());
                            return;
                        }
                    }
                }
            }
            source => match source.download_url() {
                Ok(url) => Mode::Direct(source, url),
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            },
        };
        let now = js_sys::Date::now();
        let interval = if matches!(mode, Mode::Direct(..)) {
            350.0
        } else {
            1000.0
        };
        if now - last_request.get() < interval {
            status.set("Please wait briefly before the next source request.".into());
            return;
        }
        last_request.set(now);
        let abort = match web_sys::AbortController::new() {
            Ok(a) => a,
            Err(_) => {
                status.set("Cannot initialize cancellable source request".into());
                return;
            }
        };
        *controller.borrow_mut() = Some(abort.clone());
        busy.set(true);
        status.set("Downloading and validating source…".into());
        let alive = alive.clone();
        let controller = controller.clone();
        let expected_generation = generation.get_untracked();
        spawn_local(async move {
            let window = web_sys::window().unwrap();
            let timeout_abort = abort.clone();
            let timeout = Closure::<dyn FnMut()>::new(move || timeout_abort.abort());
            let timer = window
                .set_timeout_with_callback_and_timeout_and_arguments_0(
                    timeout.as_ref().unchecked_ref(),
                    45_000,
                )
                .ok();
            enum Outcome {
                Catalog(AddgeneCatalog),
                Import(ImportReport),
            }
            let url = match &mode {
                Mode::Direct(_, url) | Mode::Catalog(_, url) | Mode::Selection(_, _, url) => url,
            };
            let auth = if !matches!(mode, Mode::Direct(..))
                && web_sys::Url::new(url)
                    .ok()
                    .is_some_and(|u| u.origin() == "https://api.developers.addgene.org")
            {
                Some(password.as_str())
            } else {
                None
            };
            let result = fetch_bounded(url, &abort.signal(), auth)
                .await
                .and_then(|bytes| {
                    let retrieved = js_sys::Date::new_0()
                        .to_iso_string()
                        .as_string()
                        .unwrap_or_default();
                    match &mode {
                        Mode::Catalog(id, _) => parse_addgene_catalog(*id, &bytes)
                            .map(Outcome::Catalog)
                            .map_err(|e| e.to_string()),
                        Mode::Selection(catalog, id, _) => catalog
                            .import_download(*id, &bytes, &retrieved)
                            .map(Outcome::Import)
                            .map_err(|e| e.to_string()),
                        Mode::Direct(source, _) => import_download(source, &bytes, &retrieved)
                            .map(Outcome::Import)
                            .map_err(|e| e.to_string()),
                    }
                });
            if let Some(timer) = timer {
                window.clear_timeout_with_handle(timer);
            }
            abort.abort();
            if !alive.get() || generation.get_untracked() != expected_generation {
                return;
            }
            controller.borrow_mut().take();
            busy.set(false);
            match result {
                Ok(Outcome::Catalog(found)) => {
                    let count = found.choices().len();
                    status.set(if count==0 {"Addgene returned no available sequence records for this plasmid.".into()}else{format!("Found {count} Addgene sequence records. Choose a record and click Import source again.")});
                    catalog.set(Some(found));
                    selected.set(String::new());
                }
                Ok(Outcome::Import(report)) => {
                    if let Err(e) = on_import.run(report) {
                        status.set(e);
                    }
                }
                Err(error) => status.set(error),
            }
        });
    };
    view! {
        <p>"Import an annotated nucleotide record from NCBI, ENA or DDBJ, an Addgene plasmid sequence, or a direct HTTPS sequence file. A new document retains the source identity and original record evidence."</p>
        <label>"Sequence URL or provider identifier"<input class="source-url-input" aria-label="Sequence source URL or accession" type="text" prop:value=move||input.get() on:input=move|e|{invalidate.with_value(|f|f());input.set(event_target_value(&e));} placeholder="ena:AB000263.1 or ddbj:AB000263.1"/></label>
        <p class="helper source-provider-status">{move||provider_description(&input.get())}</p>
        <details class="source-provider-examples"><summary>"Provider examples and record scope"</summary><p class="helper">"NCBI: L09137.2 or https://www.ncbi.nlm.nih.gov/nuccore/L09137.2. ENA: ena:AB000263.1 or https://www.ebi.ac.uk/ena/browser/view/AB000263.1. DDBJ: ddbj:AB000263.1 or https://getentry.ddbj.nig.ac.jp/getentry/na/AB000263.1. Addgene: addgene:12345 or a plasmid page URL."</p><p class="helper">"ENA and DDBJ links retrieve annotated records. Specify accession.version to pin a version; an unversioned request records the actual returned version. Search results, ranges and multi-record responses are rejected for database record imports. Retrieval does not verify biological completeness or annotation correctness."</p></details>
        <Show when=move||matches!(parse_source(&input.get()),Ok(Source::Addgene { .. }))>
            <label>"Addgene API token"<input aria-label="Addgene API token" type="password" autocomplete="off" prop:value=move||token.get() on:input=move|e|{invalidate.with_value(|f|f());token.set(event_target_value(&e));}/></label>
            <p class="helper">"Requires approved catalog:retrieve-with-sequences access. The token stays in this dialog's memory and is sent only to Addgene's API. A token alone does not resolve browser CORS: the observed localhost authorization preflight fails. Use authorized native/MCP import or open a file downloaded through your Addgene account."</p>
        </Show>
        {move||catalog.get().map(|found|view!{<div class="addgene-choices"><h3>{found.name().to_owned()}</h3><label>"Sequence record"<select aria-label="Addgene sequence record" prop:value=move||selected.get() on:change=move|e|{generation.update(|n|*n=n.wrapping_add(1));selection_reset.with_value(|c|if let Some(abort)=c.borrow_mut().take(){abort.abort();});busy.set(false);status.set(String::new());selected.set(event_target_value(&e));}><option value="">"Choose a sequence…"</option>{found.choices().iter().cloned().map(|choice|view!{<option value=choice.sequence_id.to_string()>{format!("{} · {:?} · {:?} · {} bases · {}",choice.sequence_id,choice.origin,choice.completeness,choice.length.map(|n|n.to_string()).unwrap_or_else(||"unknown length".into()),choice.description)}</option>}).collect_view()}</select></label><p class="helper">"Partial records are not complete plasmids. Provider categories are retained as source evidence."</p></div>})}
        <p class="helper">"Downloads are limited to 32 MiB and 45 seconds. Browser CORS restrictions may require a native/MCP import. No website login cookies are sent."</p>
        <p class="source-import-status" role="status">{move||status.get()}</p>
        <button class="button primary" disabled=move||busy.get() on:click=import>"Import source"</button>
    }
}

/// Do not fall back to direct fetch after a hosted/native failure. A broken native
/// bridge or missing lab login must stay explicit rather than change authority.
#[component]
pub fn SourceImport(on_import: Callback<ImportReport, Result<(), String>>) -> impl IntoView {
    if crate::io::source_gateway_client::enabled() {
        view! { <HostedSourceImport on_import=on_import/> }.into_any()
    } else {
        view! { <DirectSourceImport on_import=on_import/> }.into_any()
    }
}

#[component]
fn HostedSourceImport(on_import: Callback<ImportReport, Result<(), String>>) -> impl IntoView {
    use collagene_source_gateway::wire;
    let client = StoredValue::new_local(crate::io::source_gateway_client::Client::new());
    let input = RwSignal::new(String::new());
    let token = RwSignal::new(String::new());
    let catalog = RwSignal::new(None::<wire::Catalog>);
    let selected = RwSignal::new(String::new());
    let generation = RwSignal::new(0u64);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let alive = StoredValue::new_local(Rc::new(Cell::new(true)));
    let invalidate = Callback::new(move |()| {
        generation.update(|n| *n = n.wrapping_add(1));
        client.with_value(|c| {
            if let Ok(c) = c {
                c.reset();
            }
        });
        catalog.set(None);
        selected.set(String::new());
        busy.set(false);
        status.set(String::new());
    });
    on_cleanup(move || {
        alive.with_value(|a| a.set(false));
        client.with_value(|c| {
            if let Ok(c) = c {
                c.close();
            }
        });
    });
    let import = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let client = match client.get_value() {
            Ok(c) => c,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let source = match parse_source(&input.get_untracked()) {
            Ok(s) => s,
            Err(_) => {
                status.set("Enter a supported provider identifier or HTTPS sequence URL.".into());
                return;
            }
        };
        let operation = match source {
            Source::Addgene { plasmid_id: _ } => {
                let password = token.get_untracked();
                if password.is_empty() {
                    status.set("Enter an approved Addgene API token, or import a sequence file downloaded through your Addgene account.".into());
                    return;
                }
                if !password.is_empty()
                    && Config::default().with_addgene_api_token(&password).is_err()
                {
                    status.set("Enter a valid Addgene API token.".into());
                    return;
                }
                let addgene_token = if password.is_empty() {
                    None
                } else {
                    Some(password)
                };
                if let Some(found) = catalog.get_untracked() {
                    let id = match selected.get_untracked().parse::<u64>() {
                        Ok(id) if found.choices.iter().any(|v| v.sequence_id == id) => id,
                        _ => {
                            status
                                .set("Choose the exact Addgene sequence record to import.".into());
                            return;
                        }
                    };
                    wire::Operation::Select {
                        handle: found.handle,
                        response_blake3: found.response_blake3,
                        sequence_id: id,
                        addgene_token,
                    }
                } else {
                    wire::Operation::Catalog {
                        source: input.get_untracked(),
                        addgene_token,
                    }
                }
            }
            _ => wire::Operation::Import {
                source: input.get_untracked(),
                addgene_token: None,
            },
        };
        busy.set(true);
        status.set("Retrieving and validating source…".into());
        let expected_generation = generation.get_untracked();
        let alive = alive.get_value();
        spawn_local(async move {
            let result = client.request(operation).await;
            if !alive.get() || generation.get_untracked() != expected_generation {
                if let Ok(wire::Output::Catalog(found)) = result {
                    client.dispose_catalog(found.handle);
                }
                return;
            }
            busy.set(false);
            match result {
                Ok(wire::Output::Catalog(found)) => {
                    status.set(if found.choices.is_empty(){"No available sequence records were returned. Try an authorized downloaded sequence file.".into()}else{format!("Found {} records. Review their scope and explicitly choose one to import.",found.choices.len())});
                    catalog.set(Some(found));
                    selected.set(String::new());
                }
                Ok(wire::Output::Imported(report)) => {
                    // App's original epoch-capturing callback remains the final
                    // workspace commit guard. The host response grants no bypass.
                    client.reset();
                    catalog.set(None);
                    selected.set(String::new());
                    if let Err(error) = on_import.run(report) {
                        status.set(error);
                    }
                }
                Ok(_) => {
                    status.set("Source gateway returned an unexpected operation result.".into())
                }
                Err(error) => status.set(error),
            }
        });
    });
    view! {
        <p>"Import an annotated nucleotide record from NCBI, ENA or DDBJ, an Addgene sequence, or a direct HTTPS sequence file. The host retrieves and validates the original source evidence."</p>
        <label>"Sequence URL or provider identifier"<input class="source-url-input" type="text" maxlength="8192" aria-label="Sequence source URL or accession" prop:value=move||input.get() on:input=move|e|{invalidate.run(());input.set(event_target_value(&e));} placeholder="ena:AB000263.1 or addgene:12345"/></label>
        <p class="helper source-provider-status">{move||provider_description(&input.get())}</p>
        <p class="helper">{if crate::platform::desktop_rpc::is_desktop(){"Retrieved by the native application; no lab login is needed."}else{"Uses the same authenticated lab session as the Lab library. Connect there first. Your lab key is never shared with a sequence provider."}}</p>
        <Show when=move||matches!(parse_source(&input.get()),Ok(Source::Addgene{..}))>
            <label>"Addgene API token (dialog memory only)"<input aria-label="Addgene API token" type="password" autocomplete="off" maxlength="4096" prop:value=move||token.get() on:input=move|e|{invalidate.run(());token.set(event_target_value(&e));}/></label>
            <p class="helper">"Requires approved Addgene catalog access. This token is sent to the host for this request and only used for Addgene; it is never saved in the library or native vault. You can also open a sequence file downloaded through your Addgene account."</p>
        </Show>
        {move||catalog.get().map(|found|view!{
            <div class="addgene-choices"><h3>{found.name}</h3>
                <label>"Sequence record"<select aria-label="Addgene sequence record" prop:value=move||selected.get() on:change=move|e|{
                    generation.update(|n|*n=n.wrapping_add(1));
                    client.with_value(|c|if let Ok(c)=c{c.cancel_active();});
                    busy.set(false);status.set(String::new());selected.set(event_target_value(&e));
                }><option value="">"Choose a sequence…"</option>
                    {found.choices.into_iter().map(|choice|view!{<option value=choice.sequence_id.to_string()>{format!("{} · {:?} · {:?} · {} bases · {}",choice.sequence_id,choice.origin,choice.completeness,choice.length.map(|n|n.to_string()).unwrap_or_else(||"unknown length".into()),choice.description)}</option>}).collect_view()}
                </select></label>
                <p class="helper">"Partial records are not complete plasmids. Provider categories and original source evidence are retained. The catalog expires; refresh if it is no longer available."</p>
            </div>
        })}
        <p class="helper">"Downloads remain bounded to 32 MiB; the host enforces provider and timeout limits. Changing the source cancels retrieval and disposes its catalog."</p>
        <p class="source-import-status" role="status" aria-live="polite">{move||status.get()}</p>
        <button class="button secondary" disabled=move||!busy.get() on:click=move|_|{invalidate.run(());status.set("Source retrieval cancelled.".into());}>"Cancel retrieval"</button>
        <button class="button secondary" disabled=move||busy.get()||catalog.with(Option::is_none) on:click=move|_|{invalidate.run(());import.run(());}>"Refresh catalog"</button>
        <button class="button primary" disabled=move||busy.get() on:click=move|_|import.run(())>{move||if catalog.with(Option::is_some){"Import selected sequence"}else{"Import source"}}</button>
    }
}
