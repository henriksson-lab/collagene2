//! Immutable VCF/BCF browsing and explicit literal-allele reference checking.
//! A source genotype never chooses or installs an allele automatically.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    Alphabet, Molecule,
    protocol::{AnalysisJob, AnalysisResult, Request},
    variant_reference::{self as reference, Outcome},
    variants_session::{self as session, ExportFormat, FileFormat},
};
use leptos::prelude::*;
use serde::Serialize;
use std::{
    io::{self, Write},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

const PAGE: usize = 100;
const TEXT_LIMIT: usize = 2000;
fn clip(value: &str) -> String {
    let prefix = value.chars().take(TEXT_LIMIT).collect::<String>();
    if prefix.len() < value.len() {
        format!(
            "{prefix}\n… preview clipped; original/normalized source downloads retain complete fields."
        )
    } else {
        prefix
    }
}
fn field(row: &session::Row, index: usize, limit: usize) -> String {
    row.vcf
        .trim_end()
        .split('\t')
        .nth(index)
        .unwrap_or(".")
        .chars()
        .take(limit)
        .collect()
}
fn format_name(format: FileFormat) -> &'static str {
    match format {
        FileFormat::Vcf => "VCF",
        FileFormat::Bcf => "BCF",
    }
}
fn export_name(format: ExportFormat) -> &'static str {
    match format {
        ExportFormat::Original => "original",
        ExportFormat::Vcf => "VCF",
        ExportFormat::Bcf => "BCF",
    }
}
fn coordinates(report: &reference::Report) -> Vec<(usize, usize)> {
    report
        .start
        .zip(report.end)
        .filter(|(a, b)| a < b)
        .into_iter()
        .collect()
}
struct JsonBytes(Vec<u8>);
impl Write for JsonBytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > session::MAX_EXPORT_BYTES)
        {
            return Err(io::Error::other(
                "Variant evidence download exceeds byte limit",
            ));
        }
        self.0
            .try_reserve(bytes.len())
            .map_err(|_| io::Error::other("Variant evidence allocation"))?;
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
#[derive(Clone)]
struct Mapping {
    document: usize,
    event: String,
    index: usize,
    chromosome: String,
    alternate: usize,
    serial: u32,
}
struct FreshCheck {
    mapping: Mapping,
    options: reference::Options,
    value: Arc<session::Checked>,
}
struct FreshPreview {
    checked: Arc<FreshCheck>,
    value: Arc<session::Previewed>,
}
async fn ask(
    worker: Rc<CachedWorker>,
    id: u32,
    revision: String,
    job: AnalysisJob,
) -> Result<AnalysisResult, String> {
    let response = worker
        .ask(Request {
            id,
            revision: revision.clone(),
            job,
        })
        .await?;
    if response.id != id || response.revision != revision {
        return Err("Variant response belongs to different inputs.".into());
    }
    match response.result {
        AnalysisResult::Error(error) => Err(error),
        value => Ok(value),
    }
}

/// The app remounts this component for a new immutable file. Document tuple item3
/// is the live document/event/epoch identity; names or equal DNA are insufficient.
#[component]
pub fn VariantPanel(
    source: Arc<Vec<u8>>,
    format: FileFormat,
    name: String,
    documents: Memo<Vec<(usize, String, String)>>,
    get_document: Callback<usize, Option<Molecule>>,
    on_navigate: Callback<(usize, Vec<(usize, usize)>)>,
    on_preview: Callback<(Molecule, String), Result<(), String>>,
) -> impl IntoView {
    let bounded = !source.is_empty() && source.len() <= session::MAX_SOURCE_BYTES;
    let source_hash = StoredValue::new(if bounded {
        blake3::hash(source.as_slice()).to_hex().to_string()
    } else {
        String::new()
    });
    let prepare_binding = StoredValue::new(if bounded {
        session::prepare_revision(source.as_slice(), format).map_err(|e| e.to_string())
    } else {
        Err("VCF/BCF source must be nonempty and at most16MiB. No partial collection was installed.".into())
    });
    let source = StoredValue::new(source);
    let filename = StoredValue::new(name.clone());
    let prepared = RwSignal::new(None::<Arc<session::Prepared>>);
    let page = RwSignal::new(None::<Arc<session::Page>>);
    let selected = RwSignal::new(None::<Arc<session::Row>>);
    let document = RwSignal::new(None::<usize>);
    let chromosome = RwSignal::new(String::new());
    let alternate = RwSignal::new(String::new());
    let checked = RwSignal::new(None::<Arc<FreshCheck>>);
    let preview = RwSignal::new(None::<Arc<FreshPreview>>);
    let converted = RwSignal::new(None::<Arc<session::Exported>>);
    let output_format = RwSignal::new(ExportFormat::Vcf);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let mapping_serial = RwSignal::new(0u32);
    let request_id = StoredValue::new(0u32);
    let client = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let context = Memo::new(move |_| {
        documents.with(|docs| {
            document.get().and_then(|chosen| {
                docs.iter()
                    .find(|(index, _, _)| *index == chosen)
                    .map(|(index, _, event)| (*index, event.clone()))
            })
        })
    });
    Effect::new(move |_| {
        let _ = context.get();
        let _ = chromosome.get();
        let _ = alternate.get();
        let _ = selected.with(|row| row.as_ref().map(|row| row.index));
        mapping_serial.update(|serial| *serial = serial.wrapping_add(1));
        checked.set(None);
        preview.set(None);
    });
    let current_mapping = move |mapping: &Mapping| {
        context.get_untracked() == Some((mapping.document, mapping.event.clone()))
            && selected
                .with_untracked(|row| row.as_ref().is_some_and(|row| row.index == mapping.index))
            && chromosome.get_untracked() == mapping.chromosome
            && alternate.get_untracked().parse::<usize>().ok() == Some(mapping.alternate)
            && mapping_serial.get_untracked() == mapping.serial
    };
    let current_check = move || {
        checked
            .get()
            .filter(|value| current_mapping(&value.mapping))
    };
    let current_preview = move || {
        preview
            .get()
            .filter(|value| current_mapping(&value.checked.mapping))
    };
    let close = Callback::new(move |reason: String| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = client.get_value() {
            worker.close();
        }
        client.set_value(None);
        busy.set(false);
        prepared.set(None);
        page.set(None);
        selected.set(None);
        checked.set(None);
        preview.set(None);
        converted.set(None);
        status.set(reason);
    });
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = client.get_value() {
            worker.close();
        }
        client.set_value(None);
    });
    let next_id = move || {
        let id = request_id.get_value().wrapping_add(1);
        request_id.set_value(id);
        id
    };
    let page_alive = alive.clone();
    let load_page = Callback::new(move |start: usize| {
        if busy.get_untracked() {
            return;
        }
        let (Some(worker), Some(summary)) = (client.get_value(), prepared.get_untracked()) else {
            return;
        };
        let hash = source_hash.get_value();
        let revision = match session::page_revision(&hash, start, PAGE) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = page_alive.clone();
        busy.set(true);
        selected.set(None);
        checked.set(None);
        preview.set(None);
        status.set("Loading reported variant rows…".into());
        spawn_local(async move {
            let result = ask(
                worker,
                id,
                revision.clone(),
                AnalysisJob::VariantPage {
                    source_hash: hash.clone(),
                    start,
                    count: PAGE,
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::VariantPage(value)) => {
                    let valid = value.validate().map_err(|e| e.to_string()).and_then(|_| {
                        if value.revision == revision
                            && value.source_hash == hash
                            && value.start == start
                            && value.count == PAGE
                            && value.total == summary.summary.records
                        {
                            Ok(())
                        } else {
                            Err("Variant page source/offset binding changed.".into())
                        }
                    });
                    match valid {
                        Ok(()) => {
                            page.set(Some(Arc::new(value)));
                            status.set(String::new());
                        }
                        Err(error) => status.set(error),
                    }
                }
                Ok(_) => status.set("Unexpected variant page response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let prepare_alive = alive.clone();
    let prepare = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let revision = match prepare_binding.get_value() {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        close.run(String::new());
        chromosome.set(String::new());
        alternate.set(String::new());
        let worker = match CachedWorker::new() {
            Ok(value) => Rc::new(value),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        client.set_value(Some(worker.clone()));
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = prepare_alive.clone();
        let hash = source_hash.get_value();
        let bytes = source.with_value(|source| source.as_ref().clone());
        busy.set(true);
        status.set("Parsing the entire bounded variant file…".into());
        spawn_local(async move {
            let result = ask(
                worker,
                id,
                revision.clone(),
                AnalysisJob::PrepareVariants { bytes, format },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::VariantsPrepared(value)) => {
                    let valid = value.validate().map_err(|e| e.to_string()).and_then(|_| {
                        if value.revision == revision && value.summary.source_hash == hash {
                            Ok(())
                        } else {
                            Err("Prepared variants have a different source binding.".into())
                        }
                    });
                    match valid {
                        Ok(()) => {
                            prepared.set(Some(Arc::new(value)));
                            status.set("Complete file parsed; calls and genotypes remain source assertions.".into());
                            load_page.run(0);
                        }
                        Err(error) => close.run(error),
                    }
                }
                Ok(_) => close.run("Unexpected variant preparation response.".into()),
                Err(error) => close.run(error),
            }
        });
    });
    let check_alive = alive.clone();
    let check = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let (Some(worker), Some(row), Some((document, event))) = (
            client.get_value(),
            selected.get_untracked(),
            context.get_untracked(),
        ) else {
            status.set("Choose a record and actual DNA document.".into());
            return;
        };
        let name = chromosome.get_untracked();
        if name.is_empty() || name != row.reference_name {
            status.set(
                "Explicitly map the selected row's chromosome label to the chosen DNA document."
                    .into(),
            );
            return;
        }
        let alt = match alternate.get_untracked().parse::<usize>() {
            Ok(index) if (1..=64).contains(&index) => index,
            _ => {
                status.set(
                    "Enter the one-based alternate index,1–64; genotype fields never choose it."
                        .into(),
                );
                return;
            }
        };
        let Some(reference) = get_document.run(document) else {
            status.set("Selected reference is unavailable or exceeds the reference limits.".into());
            return;
        };
        if reference.sequence.alphabet() != Alphabet::Dna {
            status.set("Choose an actual DNA document.".into());
            return;
        }
        // Borrowed native-model guard precedes canonical hashing; the app getter
        // applies its borrowed policy before making the owned source copy too.
        if let Err(error) = collagene_core::model::bounds::model_preflight(&[&reference]) {
            status.set(error.to_string());
            return;
        }
        let reference_revision = reference.revision();
        let hash = source_hash.get_value();
        let options = reference::Options {
            source_blake3: hash.clone(),
            record_index: row.index,
            reference_name: name.clone(),
            expected_reference_revision: reference_revision,
            alternate_index: alt,
        };
        if let Err(error) = session::reference_preflight(&reference, &options) {
            status.set(error.to_string());
            return;
        }
        let revision = match session::check_revision_for_revision(
            &hash,
            &options.expected_reference_revision,
            &options,
        ) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let mapping = Mapping {
            document,
            event,
            index: row.index,
            chromosome: name,
            alternate: alt,
            serial: mapping_serial.get_untracked(),
        };
        let expected = options.clone();
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = check_alive.clone();
        busy.set(true);
        checked.set(None);
        preview.set(None);
        status.set("Checking this explicit REF/ALT against the actual DNA…".into());
        spawn_local(async move {
            let result = ask(
                worker,
                id,
                revision.clone(),
                AnalysisJob::CheckVariantReference {
                    source_hash: hash.clone(),
                    reference,
                    options,
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            if !current_mapping(&mapping) {
                status.set("Mapping changed; old REF evidence was discarded.".into());
                return;
            }
            match result {
                Ok(AnalysisResult::VariantReferenceChecked(value)) => {
                    match value.validate_bindings(&hash, &expected).and_then(|()| {
                        if value.revision == revision {
                            Ok(())
                        } else {
                            Err(collagene_core::Error::Unsupported(
                                "Variant response revision changed".into(),
                            ))
                        }
                    }) {
                        Ok(()) => {
                            checked.set(Some(Arc::new(FreshCheck {
                                mapping,
                                options: expected,
                                value: Arc::new(value),
                            })));
                            status.set(String::new());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Unexpected variant reference-check response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let preview_alive = alive.clone();
    let make_preview = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let (Some(worker), Some(fresh)) = (client.get_value(), current_check()) else {
            return;
        };
        if fresh.value.report.outcome != Outcome::ReferenceAlleleMatches {
            return;
        }
        let Some(reference) = get_document.run(fresh.mapping.document) else {
            status.set("Reference is no longer available.".into());
            return;
        };
        if let Err(error) = session::reference_preflight(&reference, &fresh.options) {
            status.set(error.to_string());
            return;
        }
        let revision =
            match session::preview_revision(&source_hash.get_value(), &reference, &fresh.options) {
                Ok(value) => value,
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = preview_alive.clone();
        busy.set(true);
        preview.set(None);
        status.set(
            "Freshly checking REF again before constructing a literal alternate design…".into(),
        );
        spawn_local(async move {
            let result = ask(
                worker,
                id,
                revision.clone(),
                AnalysisJob::PreviewVariantAllele {
                    source_hash: source_hash.get_value(),
                    reference,
                    options: fresh.options.clone(),
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            if !current_mapping(&fresh.mapping) {
                status.set("Reference event or allele choice changed; preview discarded.".into());
                return;
            }
            match result {
                Ok(AnalysisResult::VariantAllelePreview(value)) => {
                    match value
                        .validate_bindings(&source_hash.get_value(), &fresh.options)
                        .and_then(|()| {
                            if value.revision == revision {
                                Ok(())
                            } else {
                                Err(collagene_core::Error::Unsupported(
                                    "Variant preview revision changed".into(),
                                ))
                            }
                        }) {
                        Ok(()) => {
                            preview.set(Some(Arc::new(FreshPreview {
                                checked: fresh,
                                value: Arc::new(value),
                            })));
                            status.set(String::new());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Unexpected variant preview response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let export_alive = alive.clone();
    let convert = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let Some(worker) = client.get_value() else {
            return;
        };
        let hash = source_hash.get_value();
        let format = output_format.get_untracked();
        let cap = session::MAX_EXPORT_BYTES;
        let revision = match session::export_revision(&hash, format, cap) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = export_alive.clone();
        busy.set(true);
        converted.set(None);
        status.set("Preparing normalized variant exchange…".into());
        spawn_local(async move {
            let result = ask(
                worker,
                id,
                revision.clone(),
                AnalysisJob::ExportVariants {
                    source_hash: hash.clone(),
                    format,
                    max_output_bytes: cap,
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            if output_format.get_untracked() != format {
                status.set("Export choice changed; previous conversion discarded.".into());
                return;
            }
            match result {
                Ok(AnalysisResult::VariantsExported(value)) => {
                    let valid = value.validate().map_err(|e| e.to_string()).and_then(|_| {
                        if value.revision == revision
                            && value.source_hash == hash
                            && value.format == format
                            && value.max_output_bytes == cap
                        {
                            Ok(())
                        } else {
                            Err("Variant export source/format binding changed.".into())
                        }
                    });
                    match valid {
                        Ok(()) => {
                            converted.set(Some(Arc::new(value)));
                            status.set(
                                "Normalized export ready; review its notices before downloading."
                                    .into(),
                            );
                        }
                        Err(error) => status.set(error),
                    }
                }
                Ok(_) => status.set("Unexpected variant export response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let dispose_alive = alive.clone();
    let dispose = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let Some(worker) = client.get_value() else {
            return;
        };
        let hash = source_hash.get_value();
        let revision = match session::dispose_revision(&hash) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = dispose_alive.clone();
        busy.set(true);
        spawn_local(async move {
            let result = ask(
                worker,
                id,
                revision.clone(),
                AnalysisJob::DisposeVariants {
                    source_hash: hash.clone(),
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            match result {
                Ok(AnalysisResult::VariantsDisposed(value))
                    if value.validate().is_ok()
                        && value.revision == revision
                        && value.source_hash == hash =>
                {
                    close.run("Parsed variant records released; original source retained.".into())
                }
                Ok(_) => close.run("Variant worker closed.".into()),
                Err(error) => close.run(error),
            }
        });
    });
    let download_evidence = Callback::new(move |()| {
        let Some(fresh) = current_check() else {
            return;
        };
        let design = current_preview();
        #[derive(Serialize)]
        struct Evidence<'a> {
            interpretation: &'static str,
            document_event: &'a str,
            checked: &'a session::Checked,
            preview: Option<&'a session::Previewed>,
        }
        let value = Evidence {
            interpretation: "Saved source checking and literal-design evidence only. Recompute from the original file and actual current DNA/event. No allele presence, genotype, clinical meaning, build identity or inherited chemistry is established.",
            document_event: &fresh.mapping.event,
            checked: fresh.value.as_ref(),
            preview: design.as_ref().map(|p| p.value.as_ref()),
        };
        let mut bytes = JsonBytes(vec![]);
        match serde_json::to_writer(&mut bytes, &value) {
            Ok(()) => {
                if let Err(error) = download(
                    "variant-reference-evidence.json",
                    "application/json",
                    &bytes.0,
                ) {
                    status.set(error)
                }
            }
            Err(error) => status.set(error.to_string()),
        }
    });
    prepare.run(());
    view! {<section class="variants-panel" aria-label="VCF/BCF variant evidence">
        <h2>{format!("{} · {} variant evidence",name,format_name(format))}</h2>
        <p class="helper">"Imported calls, quality, FILTER, INFO, FORMAT and genotypes are declarations. Names and REF agreement do not establish genome build, allele presence, ploidy, function or clinical meaning. Only an explicitly selected literal ALT can become a new reference-only design preview."</p>
        <div class="variant-source-actions"><button class="button secondary" disabled=!bounded on:click=move|_|{if let Err(error)=source.with_value(|bytes|download(&filename.get_value(),"application/octet-stream",bytes.as_slice())){status.set(error);}}>"Download original VCF/BCF"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|prepare.run(())>"Reopen variant source"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|close.run("Task cancelled and variant worker closed.".into())>"Cancel and close variant session"</button><button class="text-tool" disabled=move||busy.get()||prepared.with(Option::is_none) on:click=move|_|dispose.run(())>"Release parsed variant records"</button></div>
        <p class="variant-status" role="status">{move||status.get()}</p>
        {move||prepared.get().map(|value|{let summary=&value.summary;view!{<div class="variant-summary"><p>{format!("Complete source parsed: {} bytes · {} records · retained-buffer charge {} bytes (not process RSS) · hash {}",summary.source_bytes,summary.records,summary.retained_charge,summary.source_hash)}</p><details><summary>{format!("Reported header: {} samples / {} contigs (not verified identities)",summary.samples.len(),summary.contigs.len())}</summary><p>{format!("Sample-column order{}: {}",if summary.samples.len()>16{" (first16)"}else{""},summary.samples.iter().take(16).map(|sample|sample.chars().take(80).collect::<String>()).collect::<Vec<_>>().join(", "))}</p>{summary.contigs.iter().take(100).map(|contig|view!{<p>{format!("{} · declared length {:?}",contig.name,contig.length)}</p>}).collect_view()}<p class="helper">"Header preview is bounded; source and normalized downloads retain all accepted metadata."</p></details></div>}})}
        <div class="variant-pages"><button class="text-tool" disabled=move||busy.get()||page.with(|p|p.as_ref().is_none_or(|p|p.start==0)) on:click=move|_|if let Some(value)=page.get_untracked(){load_page.run(value.start.saturating_sub(PAGE));}>"Previous 100 variants"</button><span>{move||page.with(|p|p.as_ref().map(|p|format!("{}–{} of {}",if p.rows.is_empty(){0}else{p.start+1},p.start+p.rows.len(),p.total)).unwrap_or_default())}</span><button class="text-tool" disabled=move||busy.get()||page.with(|p|p.as_ref().is_none_or(|p|p.start+p.rows.len()>=p.total)) on:click=move|_|if let Some(value)=page.get_untracked(){load_page.run(value.start+PAGE);}>"Next 100 variants"</button></div>
        <div class="variant-table-wrap"><table class="variant-rows"><thead><tr><th>"Record"</th><th>"Reported chromosome"</th><th>"POS (1-based)"</th><th>"Declared REF"</th><th>"Declared ALT"</th><th>"QUAL / FILTER"</th></tr></thead><tbody><For each=move||page.with(|p|p.as_ref().map(|p|p.rows.iter().cloned().map(Arc::new).collect::<Vec<_>>()).unwrap_or_default()) key=|row|row.index children={move|row:Arc<session::Row>|{let index=row.index;let chosen=row.clone();view!{<tr data-variant-index=index class:active=move||selected.with(|row|row.as_ref().is_some_and(|row|row.index==index))><td><button class="text-tool" disabled=move||busy.get() on:click=move|_|{selected.set(Some(chosen.clone()));alternate.set(String::new());checked.set(None);preview.set(None);}>{format!("Select {}",index+1)}</button></td><td>{row.reference_name.clone()}</td><td>{row.position_1based.map(|p|p.to_string()).unwrap_or_else(||"unavailable".into())}</td><td>{field(&row,3,80)}</td><td>{field(&row,4,80)}</td><td>{format!("{} / {}",field(&row,5,40),field(&row,6,80))}</td></tr>}}}/></tbody></table></div>
        {move||selected.get().map(|row|view!{<details open class="variant-record"><summary>{format!("Canonical VCF record {} · bounded text preview",row.index+1)}</summary><pre>{clip(&row.vcf)}</pre><p class="helper">"Columns are CHROM POS ID REF ALT QUAL FILTER INFO FORMAT followed by header-ordered sample/genotype columns. Preview clips long fields; this is normalized syntax, not the original lexical byte artifact. Genotypes never select an alternate here."</p><p>{format!("INFO (prefix): {}",field(&row,7,500))}</p><p>{format!("FORMAT (prefix): {}",field(&row,8,500))}</p><pre>{clip(&row.vcf.trim_end().split('\t').skip(9).take(8).enumerate().map(|(i,value)|format!("Sample column {}: {}\n",i+1,value.chars().take(180).collect::<String>())).collect::<String>())}</pre></details>})}
        <div class="variant-mapping"><label>"Actual DNA reference"<select aria-label="Variant DNA reference" prop:value=move||document.get().map(|index|index.to_string()).unwrap_or_default() on:change=move|event|document.set(event_target_value(&event).parse().ok())><option value="">"Choose an actual DNA document…"</option>{move||documents.get().into_iter().map(|(index,label,_)|view!{<option value=index.to_string()>{label}</option>}).collect_view()}</select></label><label>"Explicit reported chromosome mapping"<input aria-label="Variant chromosome mapping" maxlength="1024" placeholder="Enter the selected record's exact chromosome label" prop:value=move||chromosome.get() on:input=move|event|chromosome.set(event_target_value(&event))/></label><label>"Selected ALT index (one based)"<input type="number" min="1" max="64" aria-label="Variant alternate index" placeholder="1 selects first ALT, independently of genotype" prop:value=move||alternate.get() on:input=move|event|alternate.set(event_target_value(&event))/></label><button class="button primary" disabled=move||busy.get()||selected.with(Option::is_none)||context.with(Option::is_none)||chromosome.with(String::is_empty)||alternate.with(String::is_empty) on:click=move|_|check.run(())>"Check selected REF against DNA"</button></div>
        {move||current_check().map(|fresh|{let report=&fresh.value.report;let matches=report.outcome==Outcome::ReferenceAlleleMatches;let nav=fresh.clone();let parts=coordinates(report);view!{<section class="variant-reference-check"><h3>"Fresh actual-reference check"</h3><p>{format!("{:?}",report.outcome)}</p><p>{format!("Checked forward reference interval {:?}..{:?} (zero based); explicitly selected ALT {}",report.start,report.end,fresh.options.alternate_index)}</p><details><summary>"Reference binding"</summary><code>{fresh.mapping.event.clone()}</code></details><p>"Declared REF"</p><pre>{clip(&report.declared_reference)}</pre><p>"Actual mapped DNA"</p><pre>{report.actual_reference.as_deref().map(clip).unwrap_or_else(||"Unavailable; no usable forward span inferred.".into())}</pre><p>"Declared selected ALT"</p><pre>{clip(&report.declared_alternative)}</pre><p class="helper">{report.interpretation.clone()}</p><button class="text-tool" disabled=parts.is_empty() on:click=move|_|if current_mapping(&nav.mapping){on_navigate.run((nav.mapping.document,parts.clone()));}>"Show checked reference interval"</button><button class="button secondary" disabled=move||busy.get()||!matches on:click=move|_|make_preview.run(())>"Construct fresh alternate design preview"</button><button class="text-tool" on:click=move|_|download_evidence.run(())>"Download variant checking evidence"</button></section>}})}
        {move||current_preview().map(|design|{let open=design.clone();let molecule=&design.value.preview.molecule;view!{<section class="variant-design-preview"><h3>"Literal alternate design preview"</h3><p>{format!("{} bases · {:?} · {} annotations",molecule.sequence.len(),molecule.topology,molecule.features.len())}</p><p class="helper">"This is a newly derived reference-only design. It is not prepared DNA, an observed genotype, a confirmed variant or inherited chemical state. Opening it never changes the source reference document."</p><button class="button secondary" disabled=move||busy.get() on:click=move|_|{if !current_mapping(&open.checked.mapping){return;}let molecule=&open.value.preview.molecule;if molecule.metadata.get("source:legacy:chemistry_unmapped").map(String::as_str)!=Some("true")||molecule.metadata.get("variant:reference_only").map(String::as_str)!=Some("true"){status.set("Preview lost its reference-only chemistry guards.".into());return;}if let Err(error)=on_preview.run((molecule.clone(),open.checked.mapping.event.clone())){status.set(error);}}>"Open alternate design preview"</button></section>}})}
        <div class="variant-export"><label>"Normalized exchange format"<select aria-label="Variant export format" prop:value=move||export_name(output_format.get()) on:change=move|event|{output_format.set(if event_target_value(&event)=="BCF"{ExportFormat::Bcf}else{ExportFormat::Vcf});converted.set(None);}><option value="VCF">"VCF"</option><option value="BCF">"BCF"</option></select></label><button class="button secondary" disabled=move||busy.get()||prepared.with(Option::is_none) on:click=move|_|convert.run(())>"Prepare normalized variant export"</button>{move||converted.get().map(|value|{let artifact=value.clone();view!{<div class="variant-export-notices"><p>"Normalized exchange may change lexical formatting. Original is the byte-preserving artifact; no calls or genotypes are added to force conversion."</p>{value.notices.iter().take(20).map(|notice|view!{<p>{notice.clone()}</p>}).collect_view()}<button class="text-tool" on:click=move|_|{if artifact.source_hash!=source_hash.get_value(){return;}let filename=format!("variants.normalized.{}",if artifact.format==ExportFormat::Bcf{"bcf"}else{"vcf"});if let Err(error)=download(&filename,"application/octet-stream",&artifact.bytes){status.set(error);}}>"Download normalized variant file"</button></div>}})}</div>
    </section>}
}
