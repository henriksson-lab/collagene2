//! Read-only exports of exact released events. Protocols use fresh native checks;
//! portable sequence-only artifacts remain explicitly historical evidence.
use crate::worker::CachedWorker;
use collagene_core::{
    combinatorial_export::{self, Archive, Options},
    combinatorial_release::{self, Manifest},
    protocol::{AnalysisJob, AnalysisResult, Request},
    workspace::Workspace,
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone)]
struct Loaded {
    hash: String,
    epoch: u64,
    manifest: Arc<Manifest>,
}

#[component]
pub fn AssemblyArchiveView(
    snapshot: Callback<(), Result<Workspace, String>>,
    epoch: RwSignal<u64>,
    selected: RwSignal<Option<String>>,
    busy: RwSignal<bool>,
    attachment_available: Callback<String, bool>,
) -> impl IntoView {
    let loaded = RwSignal::new(None::<Arc<Loaded>>);
    let releases = RwSignal::new(Vec::<combinatorial_release::Summary>::new());
    let defaults = Options::default();
    let genbank = RwSignal::new(defaults.genbank);
    let fasta = RwSignal::new(defaults.fasta);
    let csv = RwSignal::new(defaults.primer_csv);
    let tsv = RwSignal::new(defaults.primer_tsv);
    let html = RwSignal::new(defaults.protocol_html);
    let pdf = RwSignal::new(defaults.protocol_pdf);
    let subset = RwSignal::new(false);
    let chosen = RwSignal::new(std::collections::BTreeSet::<String>::new());
    let row_page = RwSignal::new(0usize);
    let release_page = RwSignal::new(0usize);
    let prepared = RwSignal::new(None::<Arc<Archive>>);
    let exporting = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let cancel = Callback::new(move |()| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        exporting.set(false);
    });
    Effect::new(move || {
        let hash = selected.get();
        let now = epoch.get();
        if busy.get() {
            return;
        }
        if loaded.with_untracked(|current| {
            current
                .as_ref()
                .is_some_and(|value| Some(&value.hash) == hash.as_ref() && value.epoch == now)
        }) {
            return;
        }
        cancel.run(());
        prepared.set(None);
        row_page.set(0);
        let result = (|| -> Result<Option<Arc<Loaded>>, String> {
            let Some(hash) = hash else {
                return Ok(None);
            };
            let frame = snapshot.run(())?;
            let manifest = combinatorial_release::load(&frame, &hash).map_err(|e| e.to_string())?;
            chosen.set(
                manifest
                    .rows
                    .iter()
                    .map(|row| row.construct_id.clone())
                    .collect(),
            );
            Ok(Some(Arc::new(Loaded {
                hash,
                epoch: now,
                manifest: Arc::new(manifest),
            })))
        })();
        match result {
            Ok(value) => {
                loaded.set(value);
                status.set(String::new());
            }
            Err(error) => {
                loaded.set(None);
                status.set(error);
            }
        }
    });
    Effect::new(move || {
        let _ = (
            genbank.get(),
            fasta.get(),
            csv.get(),
            tsv.get(),
            html.get(),
            pdf.get(),
            subset.get(),
            chosen.get(),
        );
        cancel.run(());
        prepared.set(None);
    });
    let options = Callback::new(move |()| Options {
        genbank: genbank.get_untracked(),
        fasta: fasta.get_untracked(),
        primer_csv: csv.get_untracked(),
        primer_tsv: tsv.get_untracked(),
        protocol_html: html.get_untracked(),
        protocol_pdf: pdf.get_untracked(),
        construct_ids: if subset.get_untracked() {
            Some(chosen.with_untracked(|ids| ids.iter().cloned().collect()))
        } else {
            None
        },
    });
    let guard = Callback::new(
        move |(saved, fields): (Arc<Loaded>, Vec<u8>)| -> Result<(), String> {
            if epoch.get_untracked() != saved.epoch
                || selected.get_untracked().as_deref() != Some(saved.hash.as_str())
                || !attachment_available.run(saved.hash.clone())
            {
                return Err(
                    "The released workspace or manifest changed; prepare the archive again.".into(),
                );
            }
            if serde_json::to_vec(&options.run(())).map_err(|e| e.to_string())? != fields {
                return Err("Export options changed; prepare the archive again.".into());
            }
            Ok(())
        },
    );
    let create = {
        let alive = alive.clone();
        move |_| {
            if exporting.get_untracked() || busy.get_untracked() {
                return;
            }
            let input = (|| -> Result<_, String> {
                let saved = loaded
                    .get_untracked()
                    .ok_or("Choose a saved release first.")?;
                let config = options.run(());
                let fields = serde_json::to_vec(&config).map_err(|e| e.to_string())?;
                guard.run((saved.clone(), fields.clone()))?;
                let frame = snapshot.run(())?;
                let revision = combinatorial_export::revision(&frame, &saved.hash, &config)
                    .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((saved, config, fields, frame, revision, w))
            })();
            let (saved, config, fields, frame, revision, w) = match input {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            worker.set_value(Some(w.clone()));
            exporting.set(true);
            prepared.set(None);
            status.set("Preparing the complete selected archive in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::ExportCombinatorialOutputs {
                            workspace: frame,
                            manifest_hash: saved.hash.clone(),
                            options: config,
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                exporting.set(false);
                worker.set_value(None);
                let result = (|| -> Result<Archive, String> {
                    guard.run((saved.clone(), fields))?;
                    let response = response?;
                    if response.id != token || response.revision != revision {
                        return Err("Archive response differs from its request.".into());
                    }
                    let archive = match response.result {
                        AnalysisResult::CombinatorialOutputsArchive(archive) => archive,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected bulk archive response.".into()),
                    };
                    archive.validate().map_err(|e| e.to_string())?;
                    if archive.context_revision != revision || archive.manifest_hash != saved.hash {
                        return Err("Archive manifest or options differ from the request.".into());
                    }
                    Ok(archive)
                })();
                match result {
                    Ok(archive) => {
                        prepared.set(Some(Arc::new(archive)));
                        status.set("Archive prepared. Review format losses and evidence scope, then download explicitly.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let download = move |_| {
        let result = (|| {
            let saved = loaded
                .get_untracked()
                .ok_or("Choose the released manifest again.")?;
            let fields = serde_json::to_vec(&options.run(())).map_err(|e| e.to_string())?;
            guard.run((saved.clone(), fields.clone()))?;
            let archive = prepared
                .get_untracked()
                .ok_or("Prepare the archive first.")?;
            let artifact = crate::app::prepare_download(
                "collagene-constructs.zip",
                &archive.mime,
                &archive.bytes,
            )?;
            Ok::<_, String>((artifact, saved, fields, generation.get_untracked()))
        })();
        match result {
            Ok((artifact, saved, fields, ticket)) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    if generation.get_untracked() != ticket {
                        return Err("Archive preparation changed.".into());
                    }
                    guard.run((saved.clone(), fields.clone()))
                };
                let completion = check.clone();
                artifact.complete_guarded(check,move|result|{if completion().is_err(){return;}status.set(match result{Ok(())=>"Bulk archive saved. Sequence/primer exchange and current protocol checks retain their distinct evidence scope.".into(),Err(error)=>error});});
            }
            Err(error) => status.set(error),
        }
    };
    let refresh = move |_| match snapshot
        .run(())
        .and_then(|frame| combinatorial_release::list(&frame).map_err(|e| e.to_string()))
    {
        Ok(all) => {
            releases.set(all);
            release_page.set(0);
        }
        Err(error) => status.set(error),
    };
    view!{<details class="assembly-archive"><summary>"Released constructs and bulk export"</summary>
        <p class="helper">"Exports use exact historical released events, including their original source and target links. Sequence and primer exchange is historical data; requesting HTML/PDF protocols performs a current native check. Export never edits documents or restores experimental authority."</p>
        <button class="text-tool" disabled=move||busy.get() on:click=refresh>"Load saved releases"</button>
        {move||releases.with(|all|{
            let pages=all.len().div_ceil(10).max(1);let page=release_page.get().min(pages-1);
            let choices=all.iter().skip(page*10).take(10).map(|item|{let hash=item.manifest_hash.clone();view!{<article><p>{format!("{} released constructs · recorded {}",item.row_count,item.recorded_engine)}</p><button class="text-tool" disabled=move||busy.get() on:click=move|_|selected.set(Some(hash.clone()))>"Choose release"</button></article>}}).collect_view();
            view!{<div>{choices}<p>{format!("Release page {} of {pages}",page+1)}</p><button class="text-tool" disabled=page==0 on:click=move|_|release_page.set(page.saturating_sub(1))>"Previous releases"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|release_page.set(page+1)>"Next releases"</button></div>}
        })}
        {move||loaded.get().map(|saved|{
            let all=&saved.manifest.rows;let pages=all.len().div_ceil(20).max(1);let page=row_page.get().min(pages-1);
            let rows=all.iter().skip(page*20).take(20).map(|row|{
                let id=row.construct_id.clone();let check=id.clone();
                view!{<label><input type="checkbox" disabled=move||!subset.get()||exporting.get() prop:checked=move||chosen.with(|all|all.contains(&check)) on:change=move|event|{let selected=event_target_checked(&event);chosen.update(|all|{if selected{all.insert(id.clone());}else{all.remove(&id);}});}/>{format!("{} · {} · event {}",row.options.name,row.construct_id,row.output.construct.event_id)}</label>}
            }).collect_view();
            view!{<section><p>{format!("Release manifest {} · {} constructs",saved.hash,all.len())}</p>
                <label><input type="checkbox" prop:checked=move||subset.get() on:change=move|event|subset.set(event_target_checked(&event))/>"Export a selected subset"</label>
                {rows}<p>{format!("Construct page {} of {pages}",page+1)}</p><button class="text-tool" disabled=page==0 on:click=move|_|row_page.set(page.saturating_sub(1))>"Previous constructs"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|row_page.set(page+1)>"Next constructs"</button>
            </section>}
        })}
        <div class="assembly-export-formats">{[("GenBank",genbank),("FASTA",fasta),("Primer CSV",csv),("Primer TSV",tsv),("HTML protocols · current check",html),("PDF protocols · current check",pdf)].into_iter().map(|(label,value)|view!{<label><input type="checkbox" prop:checked=move||value.get() on:change=move|event|value.set(event_target_checked(&event))/>{label}</label>}).collect_view()}</div>
        <button class="button secondary" disabled=move||busy.get()||exporting.get()||loaded.get().is_none() on:click=create>"Prepare selected bulk archive"</button>
        <button class="text-tool" disabled=move||!exporting.get() on:click=move|_|{cancel.run(());status.set("Archive generation cancelled. No download was created.".into());}>"Cancel archive"</button>
        <ArchiveFiles prepared/>
        <button class="button primary" disabled=move||prepared.get().is_none() on:click=download>"Download prepared ZIP"</button>
        <p role="status" class="helper">{move||status.get()}</p>
    </details>}.into_any()
}

#[component]
fn ArchiveFiles(prepared: RwSignal<Option<Arc<Archive>>>) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let selected = RwSignal::new(None::<usize>);
    let note_page = RwSignal::new(0usize);
    let warnings_page = RwSignal::new(0usize);
    Effect::new(move || {
        if prepared.get().is_none() {
            page.set(0);
            selected.set(None);
            note_page.set(0);
            warnings_page.set(0);
        }
    });
    view!{<section>{move||prepared.get().map(|archive|{
        let pages=archive.files.len().div_ceil(40).max(1);let current=page.get().min(pages-1);
        let files=archive.files.iter().enumerate().skip(current*40).take(40).map(|(index,file)|view!{<article><p>{format!("{} · {} bytes · {}",file.path,file.bytes,file.mime)}</p><button class="text-tool" on:click=move|_|{selected.set(Some(index));note_page.set(0);}>"Review file details and losses"</button></article>}).collect_view();
        let notes=archive.warnings.len().div_ceil(40).max(1);let note=warnings_page.get().min(notes-1);
        let warnings=archive.warnings.iter().skip(note*40).take(40).map(|warning|view!{<li>{warning.clone()}</li>}).collect_view();
        view!{<section><h4>"Prepared archive review"</h4><p>{format!("{} bytes · {} files · ZIP hash {}",archive.bytes.len(),archive.files.len(),archive.blake3)}</p>
            <ul>{warnings}</ul><p>{format!("Archive note page {} of {notes}",note+1)}</p><button class="text-tool" disabled=note==0 on:click=move|_|warnings_page.set(note.saturating_sub(1))>"Previous archive notes"</button><button class="text-tool" disabled={note+1>=notes} on:click=move|_|warnings_page.set(note+1)>"Next archive notes"</button>
            {files}<p>{format!("File page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous files"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next files"</button>
        </section>}
    })}
    {move||selected.get().and_then(|index|prepared.with(|archive|archive.as_ref().and_then(|archive|archive.files.get(index).cloned()))).map(|file|{
        let pages=file.warnings.len().div_ceil(40).max(1);let current=note_page.get().min(pages-1);
        let warnings=file.warnings.iter().skip(current*40).take(40).map(|warning|view!{<li>{warning.clone()}</li>}).collect_view();
        view!{<section><h4>{file.path}</h4><p>{format!("File hash {} · construct {}",file.blake3,file.construct_id.as_deref().unwrap_or("archive metadata"))}</p><ul>{warnings}</ul><p>{format!("File loss-note page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|note_page.set(current.saturating_sub(1))>"Previous loss notes"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|note_page.set(current+1)>"Next loss notes"</button></section>}
    })}</section>}.into_any()
}
