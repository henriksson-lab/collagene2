//! User-owned feature libraries and reviewable automatic annotation proposals.
use crate::worker;
use collagene_core::{
    Molecule,
    feature_library::{self, FeatureLibrary, ScanOptions, ScanReport},
    protocol::{AnalysisJob, AnalysisResult, Request},
};
use leptos::{html, prelude::*};
use std::{
    collections::BTreeSet,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::{JsFuture, spawn_local};

#[derive(Clone)]
pub struct LibraryRevision {
    pub revision: String,
    pub library: Arc<FeatureLibrary>,
}
impl LibraryRevision {
    pub fn new(library: FeatureLibrary) -> Self {
        Self {
            revision: library.revision(),
            library: Arc::new(library),
        }
    }
}

async fn analyze(job: AnalysisJob, revision: String, id: u32) -> Result<AnalysisResult, String> {
    let response = worker::analyze(Request {
        id,
        revision: revision.clone(),
        job,
    })
    .await?;
    if response.id != id || response.revision != revision {
        return Err("Annotation worker returned a different document binding.".into());
    }
    match response.result {
        AnalysisResult::Error(error) => Err(error),
        result => Ok(result),
    }
}
#[component]
pub fn FeatureLibraryPanel(
    molecule: Memo<Molecule>,
    selected_feature: RwSignal<Option<usize>>,
    libraries: RwSignal<Vec<LibraryRevision>>,
    on_keep: Callback<FeatureLibrary, Result<(), String>>,
    on_commit: Callback<(String, Molecule), Result<(), String>>,
    selection: RwSignal<Option<(usize, usize)>>,
) -> impl IntoView {
    let name = RwSignal::new("My features".to_string());
    let version = RwSignal::new("1".to_string());
    let current = RwSignal::new(String::new());
    let json = RwSignal::new(String::new());
    let input = NodeRef::<html::Input>::new();
    let status = RwSignal::new(String::new());
    let mismatches = RwSignal::new("0".to_string());
    let comparisons = RwSignal::new("20000000".to_string());
    let maximum = RwSignal::new("1000".to_string());
    let both = RwSignal::new(true);
    let report = RwSignal::new(None::<Arc<ScanReport>>);
    let selected = RwSignal::new(BTreeSet::<String>::new());
    let page = RwSignal::new(0usize);
    let entry_page = RwSignal::new(0usize);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    let library = Memo::new_owning(move |_| {
        (
            libraries.with(|items| {
                items
                    .iter()
                    .find(|item| item.revision == current.get())
                    .map(|item| item.library.clone())
            }),
            true,
        )
    });
    Effect::new(move |_| {
        let _ = (
            molecule.with(Molecule::revision),
            current.get(),
            mismatches.get(),
            comparisons.get(),
            maximum.get(),
            both.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        selected.set(BTreeSet::new());
        page.set(0);
        entry_page.set(0);
        busy.set(false);
        status.set(String::new());
    });
    let keep = Callback::new(move |lib: FeatureLibrary| {
        let revision = lib.revision();
        on_keep.run(lib)?;
        current.set(revision);
        Ok::<_, String>(())
    });
    let create = move |_| match FeatureLibrary::new(name.get_untracked(), version.get_untracked())
        .map_err(|e| e.to_string())
        .and_then(|lib| keep.run(lib))
    {
        Ok(()) => status.set(
            "Empty library revision retained. Select a source annotation to add a template.".into(),
        ),
        Err(error) => status.set(error),
    };
    let add = move |_| {
        let Some(mut lib) = library.get_untracked().map(|l| l.as_ref().clone()) else {
            status.set("Create or choose a library first.".into());
            return;
        };
        let Some(index) = selected_feature.get_untracked() else {
            status.set("Select a source annotation in Features or the map first.".into());
            return;
        };
        match lib
            .add_feature(&molecule.get_untracked(), index)
            .map_err(|e| e.to_string())
            .and_then(|_| keep.run(lib))
        {
            Ok(()) => {
                status.set("Feature template retained as a new immutable library revision.".into())
            }
            Err(error) => status.set(error),
        }
    };
    let add_all = move |_| {
        let Some(mut lib) = library.get_untracked().map(|l| l.as_ref().clone()) else {
            status.set("Choose a library first.".into());
            return;
        };
        let document = molecule.get_untracked();
        let indices = document
            .features
            .iter()
            .enumerate()
            .filter(|(_, feature)| feature.kind != "source")
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        if indices.is_empty() {
            status.set("This document has no non-source annotations to add.".into());
            return;
        }
        for index in &indices {
            if let Err(error) = lib.add_feature(&document, *index) {
                status.set(format!(
                    "No library changes saved. Source annotation {}: {error}",
                    index + 1
                ));
                return;
            }
        }
        match keep.run(lib) {
            Ok(()) => status.set(format!(
                "{} source annotations retained in one new library revision.",
                indices.len()
            )),
            Err(error) => status.set(error),
        }
    };
    let import = move |_| match json
        .with_untracked(|text| feature_library::load_library(text.as_bytes()))
        .map_err(|e| e.to_string())
        .and_then(|lib| keep.run(lib))
    {
        Ok(()) => {
            json.set(String::new());
            status.set("Feature library imported. Source identifiers are recorded provenance, not independently authenticated evidence.".into());
        }
        Err(error) => status.set(error),
    };
    let file_alive = alive.clone();
    let read = move |event: web_sys::Event| {
        let target = event_target::<web_sys::HtmlInputElement>(&event);
        let file = target.files().and_then(|files| files.get(0));
        target.set_value("");
        let Some(file) = file else {
            return;
        };
        if file.size() > 8. * 1024. * 1024. {
            status.set("Feature library file exceeds the 8 MiB web import limit.".into());
            return;
        }
        let alive = file_alive.clone();
        spawn_local(async move {
            let result = JsFuture::from(file.array_buffer()).await;
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            match result {
                Ok(value) => match String::from_utf8(js_sys::Uint8Array::new(&value).to_vec()) {
                    Ok(text) => json.set(text),
                    Err(_) => status.set("Feature library file is not valid UTF-8.".into()),
                },
                Err(_) => status.set("Feature library file could not be read.".into()),
            }
        });
    };
    let scan_alive = alive.clone();
    let scan = move |_| {
        let Some(lib) = library.get_untracked() else {
            status.set("Choose a library first.".into());
            return;
        };
        let options = (|| {
            Ok::<_, String>(ScanOptions {
                both_strands: both.get_untracked(),
                max_mismatches: mismatches
                    .get_untracked()
                    .parse()
                    .map_err(|_| "Invalid mismatch limit")?,
                max_comparisons: comparisons
                    .get_untracked()
                    .parse()
                    .map_err(|_| "Invalid comparison budget")?,
                max_hits: maximum
                    .get_untracked()
                    .parse()
                    .map_err(|_| "Invalid hit limit")?,
            })
        })();
        let options = match options {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let source = molecule.get_untracked();
        let revision = source.revision();
        let lib_revision = lib.revision();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        report.set(None);
        selected.set(BTreeSet::new());
        status.set("Searching template placements…".into());
        let alive = scan_alive.clone();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ScanFeatureLibrary {
                    molecule: source,
                    library: lib.as_ref().clone(),
                    options,
                },
                revision.clone(),
                id,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::FeatureLibraryScan(value))
                    if value.source_revision == revision
                        && value.library_revision == lib_revision =>
                {
                    status.set(if value.complete {
                        "Search completed under the recorded options. Review proposed annotations."
                            .into()
                    } else {
                        "Search stopped at its recorded limit; results are incomplete.".into()
                    });
                    report.set(Some(Arc::new(value)));
                }
                Ok(_) => status.set("Unexpected feature-library response.".into()),
                Err(error) => status.set(error),
            }
        });
    };
    let apply_alive = alive;
    let apply = Callback::new(move |(): ()| {
        let (Some(lib), Some(computed)) = (library.get_untracked(), report.get_untracked()) else {
            return;
        };
        let ids = selected.get_untracked().into_iter().collect::<Vec<_>>();
        if ids.is_empty() {
            return;
        }
        let source = molecule.get_untracked();
        let revision = source.revision();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        status.set("Rechecking selected proposals…".into());
        let alive = apply_alive.clone();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ApplyFeatureLibrary {
                    molecule: source,
                    library: lib.as_ref().clone(),
                    report: computed.as_ref().clone(),
                    selected: ids,
                },
                revision.clone(),
                id,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AnnotatedDocument(edited)) => {
                    match on_commit.run((revision, edited)) {
                        Ok(()) => {
                            report.set(None);
                            selected.set(BTreeSet::new());
                            status.set("Selected annotations applied. Inferred biological annotations remain marked for review.".into());
                        }
                        Err(error) => status.set(error),
                    }
                }
                Ok(_) => status.set("Unexpected annotation-apply response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    view! {<section class="feature-library-panel"><h2>"Feature libraries"</h2><p class="helper">"Build templates from your own annotations, then search another document on both strands and across circular origins. Proposals show identity and ambiguity; adding them never establishes a cloning verification result."</p>
        <label>"Library revision"<select aria-label="Feature library revision" prop:value=move||current.get() on:change=move|e|current.set(event_target_value(&e))><option value="">"Choose a library"</option>{move||libraries.get().into_iter().map(|lib|{let revision=lib.revision;let lib=lib.library;view!{<option value=revision.clone()>{format!("{} · v{} · {} entries · {}",lib.name,lib.version,lib.entries.len(),revision.chars().take(10).collect::<String>())}</option>}}).collect_view()}</select></label>
        <details><summary>"Create a library"</summary><label>"Library name"<input aria-label="Feature library name" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label><label>"Library version"<input aria-label="Feature library version" prop:value=move||version.get() on:input=move|e|version.set(event_target_value(&e))/></label><button class="button secondary full" on:click=create>"Create feature library"</button></details>
        <button class="button secondary full" disabled=move||library.get().is_none()||selected_feature.get().is_none() on:click=add>"Add selected annotation to library"</button>
        <button class="button secondary full" disabled=move||library.get().is_none() on:click=add_all>"Add all source document annotations"</button>
        {move||library.get().map(|lib|{let count=lib.entries.len();let start=entry_page.get().min(count.saturating_sub(1)/100)*100;view!{<p class="helper">{format!("{} retained entries. Library revision {}",count,current.get())}</p><ul class="template-list">{lib.entries.iter().skip(start).take(100).map(|entry|view!{<li>{format!("{} · {} bases · {}",entry.name,entry.sequence.len(),entry.kind)}</li>}).collect_view()}</ul><Show when=move||{count>100}><div class="annotation-pages"><button class="button secondary" disabled=move||entry_page.get()==0 on:click=move|_|entry_page.update(|page|*page=page.saturating_sub(1))>"Previous library entries"</button><button class="button secondary" disabled=move||{(entry_page.get()+1)*100>=count} on:click=move|_|entry_page.update(|page|*page+=1)>"Next library entries"</button></div></Show>}})}
        <details><summary>"Import or export library JSON"</summary><button class="button secondary" on:click=move|_|if let Some(element)=input.get(){element.click();}>"Choose feature library file"</button><input node_ref=input type="file" accept=".json,.cg2features" style="display:none" on:change=read/><textarea aria-label="Feature library JSON" rows="6" prop:value=move||json.get() on:input=move|e|json.set(event_target_value(&e))></textarea><button class="button secondary" on:click=import>"Import feature library JSON"</button><button class="button secondary" disabled=move||library.get().is_none() on:click=move|_|if let Some(lib)=library.get_untracked(){match feature_library::save_library(lib.as_ref()){Ok(bytes)=>crate::app::download_with_completion("feature-library.cg2features","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}}>"Download feature library"</button></details>
        <h3>"Automatic annotation"</h3><p class="helper">"Searches whole contiguous entry templates. This does not search partial homology or splice alignments. Degenerate-compatible bases are distinguished from concrete matches."</p><label>"Allowed mismatches"<input aria-label="Feature library mismatches" type="number" min="0" prop:value=move||mismatches.get() on:input=move|e|mismatches.set(event_target_value(&e))/></label><label>"Comparison budget"<input aria-label="Feature library comparison budget" type="number" prop:value=move||comparisons.get() on:input=move|e|comparisons.set(event_target_value(&e))/></label><label>"Maximum proposals"<input aria-label="Feature library maximum proposals" type="number" prop:value=move||maximum.get() on:input=move|e|maximum.set(event_target_value(&e))/></label><label class="checkbox-label"><input aria-label="Feature library both orientations" type="checkbox" prop:checked=move||both.get() on:change=move|e|both.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Both orientations"</label><button class="button primary full" disabled=move||busy.get()||library.get().is_none() on:click=scan>"Search feature library"</button><p class="feature-library-status" role="status">{move||status.get()}</p>
        {move||report.get().map(|computed|{let for_json=computed.clone();let count=computed.hits.len();let effective=page.get().min(count.saturating_sub(1)/100);view!{<section class="feature-library-report"><h3>{format!("{} annotation proposals",count)}</h3><p>{format!("Complete: {} · {} comparisons",computed.complete,computed.comparisons)}</p><details class="annotation-evidence"><summary>"Search model and limits"</summary><ul>{computed.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul></details>{computed.hits.iter().skip(effective*100).take(100).cloned().map(|hit|{let id=hit.id;let checkbox=id.clone();let loc=collagene_core::editing::location_text(&hit.location).unwrap_or_else(|_|"unresolved".into());let start=hit.start;let end=(hit.start+hit.length).min(molecule.with_untracked(|m|m.sequence.len()));view!{<article class="library-hit"><label class="checkbox-label"><input type="checkbox" aria-label=format!("Select library proposal {id}") disabled=hit.already_annotated prop:checked=move||selected.with(|ids|ids.contains(&checkbox)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();selected.update(|ids|if checked{ids.insert(id.clone());}else{ids.remove(&id);});}/><strong>{hit.feature.label().to_string()}</strong></label><button class="text-tool" on:click=move|_|{selected_feature.set(None);selection.set(Some((start,end)));}>{loc}</button><p>{format!("Identity {:.2}–{:.2}% · coverage {:.2}% · {} mismatches · {} ambiguous-compatible bases{}",hit.identity_lower*100.,hit.identity_upper*100.,hit.coverage*100.,hit.mismatches,hit.compatible_ambiguous,if hit.already_annotated{" · already annotated"}else{""})}</p><p class="helper">"Proposed annotation · review required"</p><details class="annotation-evidence"><summary>"Match evidence"</summary><ul>{hit.warnings.into_iter().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul></details></article>}}).collect_view()}<div class="annotation-pages"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous library proposals"</button><button class="button secondary" disabled=move||{(page.get()+1)*100>=count} on:click=move|_|page.update(|p|*p+=1)>"Next library proposals"</button></div><button class="button primary full" disabled=move||busy.get()||selected.with(BTreeSet::is_empty) on:click=move|_|apply.run(())>"Apply selected library annotations"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(for_json.as_ref()){Ok(bytes)=>crate::app::download_with_completion("feature-library-search.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}>"Download library search report"</button></section>}})}
    </section>}
}

#[component]
pub fn AnnotationFilesPanel(
    molecule: Memo<Molecule>,
    source: RwSignal<String>,
    format: RwSignal<collagene_core::annotation_formats::AnnotationFormat>,
    on_commit: Callback<(String, Molecule), Result<(), String>>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
) -> impl IntoView {
    use collagene_core::annotation_formats::{AnnotationFormat, AnnotationImportReport};
    let seqid = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let input = NodeRef::<html::Input>::new();
    let report = RwSignal::new(None::<Arc<AnnotationImportReport>>);
    let selected = RwSignal::new(BTreeSet::<usize>::new());
    let page = RwSignal::new(0usize);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    Effect::new(move |_| {
        let _ = (
            molecule.with(Molecule::revision),
            source.with(|_| ()),
            format.get(),
            seqid.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        selected.set(BTreeSet::new());
        page.set(0);
        busy.set(false);
        status.set(String::new());
    });
    let retained_sources = Memo::new(move |_| {
        molecule.with(|document| {
            document
                .metadata
                .keys()
                .filter_map(|key| {
                    key.strip_prefix("annotation_import:")
                        .and_then(|tail| tail.strip_suffix(":original_base64"))
                        .map(str::to_string)
                })
                .collect::<Vec<_>>()
        })
    });
    let file_alive = alive.clone();
    let read = move |event: web_sys::Event| {
        let target = event_target::<web_sys::HtmlInputElement>(&event);
        let file = target.files().and_then(|files| files.get(0));
        target.set_value("");
        let Some(file) = file else {
            return;
        };
        if file.size() > 8. * 1024. * 1024. {
            status.set("Annotation file exceeds the 8 MiB web import limit.".into());
            return;
        }
        let alive = file_alive.clone();
        let bound = molecule.with_untracked(Molecule::revision);
        let epoch = generation.get_untracked();
        spawn_local(async move {
            let result = JsFuture::from(file.array_buffer()).await;
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            if molecule.with_untracked(Molecule::revision) != bound
                || generation.get_untracked() != epoch
            {
                status.set("Annotation source or reference changed while reading the file; choose it again.".into());
                return;
            }
            match result {
                Ok(value) => match String::from_utf8(js_sys::Uint8Array::new(&value).to_vec()) {
                    Ok(text) => source.set(text),
                    Err(_) => status.set(
                        "Annotation file is not valid UTF-8; source evidence was not changed."
                            .into(),
                    ),
                },
                Err(_) => status.set("Annotation file could not be read.".into()),
            }
        });
    };
    let import_alive = alive.clone();
    let import = move |_| {
        if source.with_untracked(|text| text.len())
            > collagene_core::annotation_formats::MAX_ANNOTATION_BYTES
        {
            status.set("Annotation source exceeds the 8 MiB import limit.".into());
            return;
        }
        let source_mol = molecule.get_untracked();
        let revision = source_mol.revision();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        report.set(None);
        selected.set(BTreeSet::new());
        status.set("Checking annotation coordinates and reference mapping…".into());
        let alive = import_alive.clone();
        let text = source.get_untracked();
        let selected_format = format.get_untracked();
        let reference_id = seqid.get_untracked();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ImportAnnotations {
                    molecule: source_mol,
                    text,
                    format: selected_format,
                    seqid: reference_id,
                },
                revision.clone(),
                id,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AnnotationImport(value))
                    if value.reference_revision == revision =>
                {
                    status.set("Annotation proposals parsed. Review the exact reference mapping and selected features before applying.".into());
                    report.set(Some(Arc::new(value)));
                }
                Ok(_) => status.set("Unexpected annotation-import response.".into()),
                Err(error) => status.set(error),
            }
        });
    };
    let apply_alive = alive.clone();
    let apply = Callback::new(move |(): ()| {
        let Some(computed) = report.get_untracked() else {
            return;
        };
        let indices = selected.get_untracked().into_iter().collect::<Vec<_>>();
        if indices.is_empty() {
            return;
        }
        let source_mol = molecule.get_untracked();
        let revision = source_mol.revision();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        status.set("Rechecking source evidence and selected annotations…".into());
        let alive = apply_alive.clone();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ApplyAnnotations {
                    molecule: source_mol,
                    report: computed.as_ref().clone(),
                    selected: indices,
                },
                revision.clone(),
                id,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AnnotatedDocument(edited)) => {
                    match on_commit.run((revision, edited)) {
                        Ok(()) => {
                            report.set(None);
                            selected.set(BTreeSet::new());
                            status.set("Selected annotations applied with original source evidence. Undo restores the previous document.".into());
                        }
                        Err(error) => status.set(error),
                    }
                }
                Ok(_) => status.set("Unexpected annotation-apply response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let export = move |_| {
        let source_mol = molecule.get_untracked();
        let revision = source_mol.revision();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        status.set("Preparing annotation export…".into());
        let alive = alive.clone();
        let selected_format = format.get_untracked();
        let reference_id = seqid.get_untracked();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ExportAnnotations {
                    molecule: source_mol,
                    format: selected_format,
                    seqid: reference_id,
                },
                revision,
                id,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AnnotationExport(value)) => {
                    let (name, mime) = match selected_format {
                        AnnotationFormat::Gff3 => ("annotations.gff3", "text/plain"),
                        AnnotationFormat::Gtf => ("annotations.gtf", "text/plain"),
                        AnnotationFormat::Bed => ("annotations.bed", "text/plain"),
                    };
                    crate::app::download_with_completion_guarded(
                        name,
                        mime,
                        value.text.as_bytes(),
                        move || {
                            if generation.try_get_untracked() == Some(id) {
                                Ok(())
                            } else {
                                Err("The annotation source changed before export acceptance."
                                    .into())
                            }
                        },
                        move |result| {
                            if status.try_get_untracked().is_none()
                                || generation.try_get_untracked() != Some(id)
                            {
                                return;
                            }
                            match result {
                                Ok(()) => status.set(if value.warnings.is_empty() {
                                    "Annotation export prepared.".into()
                                } else {
                                    value.warnings.join(" ")
                                }),
                                Err(error) => status.set(error),
                            }
                        },
                    )
                }
                Ok(_) => status.set("Unexpected annotation-export response.".into()),
                Err(error) => status.set(error),
            }
        });
    };
    view! {<section class="annotation-files-panel"><h2>"Annotation files"</h2><p class="helper">"Import GFF3, GTF or BED annotations onto the active reference. Enter its exact sequence ID from the file; matching names are never guessed. Other reference IDs are reported separately. Coordinates are reviewed before anything is added."</p>
        <label>"Annotation format"<select aria-label="Annotation file format" prop:value=move||match format.get(){AnnotationFormat::Gff3=>"gff3",AnnotationFormat::Gtf=>"gtf",AnnotationFormat::Bed=>"bed"} on:change=move|e|format.set(match event_target_value(&e).as_str(){"bed"=>AnnotationFormat::Bed,"gtf"=>AnnotationFormat::Gtf,_=>AnnotationFormat::Gff3})><option value="gff3">"GFF3"</option><option value="gtf">"GTF"</option><option value="bed">"BED"</option></select></label>
        <label>"Reference sequence ID in file"<input aria-label="Annotation reference ID" prop:value=move||seqid.get() on:input=move|e|seqid.set(event_target_value(&e)) placeholder="Exact GFF3 seqid / GTF seqname / BED chrom"/></label>
        <button class="button secondary" on:click=move|_|if let Some(element)=input.get(){element.click();}>"Choose annotation file"</button><input node_ref=input type="file" accept=".gff3,.gff,.gtf,.bed,.txt" style="display:none" on:change=read/>
        <textarea aria-label="Annotation source text" rows="8" prop:value=move||source.get() on:input=move|e|source.set(event_target_value(&e)) placeholder="Paste GFF3, GTF or BED rows"></textarea>
        <button class="button primary full" disabled=move||busy.get()||seqid.with(|s|s.trim().is_empty()) on:click=import>"Review annotation import"</button><button class="button secondary full" disabled=move||busy.get()||seqid.with(|s|s.trim().is_empty()) on:click=export>"Export current annotations"</button>
        <p class="annotation-file-status" role="status">{move||status.get()}</p>
        <Show when=move||retained_sources.with(|sources|!sources.is_empty())><details class="retained-annotation-sources"><summary>"Retained original annotation sources"</summary><p class="helper">"These are original input files. Their coordinates may describe an earlier document revision."</p>{move||retained_sources.get().into_iter().map(|hash|{let label=format!("Download retained annotation source {}",hash.chars().take(10).collect::<String>());view!{<button class="button secondary full" on:click=move|_|match collagene_core::annotation_formats::original_annotation_source(&molecule.get_untracked(),&hash){Ok(bytes)=>crate::app::download_with_completion("retained-annotations.txt","text/plain",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}>{label}</button>}}).collect_view()}</details></Show>
        {move||report.get().map(|computed|{let original=computed.clone();let for_json=computed.clone();let count=computed.proposals.len();let start=page.get().min(count.saturating_sub(1)/100)*100;view!{<section class="annotation-file-report"><h3>{format!("{} imported annotation proposals",count)}</h3><p>{format!("Reference ID: {} · {}",computed.seqid,computed.reference_revision)}</p><p>{format!("Other sequence IDs: {}",computed.skipped_seqids.join(", "))}</p><ul>{computed.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul><button class="button secondary full" on:click=move|_|selected.update(|indices|indices.extend(start..(start+100).min(count)))>"Select imported annotations on page"</button>{computed.proposals.iter().enumerate().skip(start).take(100).map(|(index,proposal)|{let label=proposal.feature.label().to_string();let location=proposal.feature.location.clone();let location_text=collagene_core::editing::location_text(&location).unwrap_or_else(|_|"unresolved".into());let feature=proposal.feature.clone();view!{<article class="annotation-file-proposal"><label class="checkbox-label"><input type="checkbox" aria-label=format!("Select imported annotation {index}") prop:checked=move||selected.with(|indices|indices.contains(&index)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();selected.update(|indices|if checked{indices.insert(index);}else{indices.remove(&index);});}/><strong>{label}</strong></label><button class="text-tool" on:click=move|_|{selected_feature.set(None);selection.set(feature.location.spans().first().map(|(a,b,_)|(*a,*b)));}>{location_text}</button><p>{format!("{} · source rows {:?}{}",proposal.feature.kind,proposal.source_rows,if proposal.feature.needs_review{" · review required"}else{""})}</p><ul>{proposal.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul></article>}}).collect_view()}<div class="annotation-pages"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous imported annotations"</button><button class="button secondary" disabled=move||{(page.get()+1)*100>=count} on:click=move|_|page.update(|p|*p+=1)>"Next imported annotations"</button></div><button class="button primary full" disabled=move||busy.get()||selected.with(BTreeSet::is_empty) on:click=move|_|apply.run(())>"Apply selected imported annotations"</button><button class="button secondary full" on:click=move|_|match original.original_bytes(){Ok(bytes)=>crate::app::download_with_completion("original-annotations.txt","text/plain",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}>"Download original annotation source"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(for_json.as_ref()){Ok(bytes)=>crate::app::download_with_completion("annotation-import.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error)}}),Err(error)=>status.set(error.to_string())}>"Download annotation import report"</button></section>}})}
    </section>}
}
