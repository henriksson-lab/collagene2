//! Immutable SAM/BAM source browsing and explicit source-bound reference checks.
//! Imported CIGAR/MAPQ/flags are assertions, never coverage or variant evidence.
use crate::platform::worker::CachedWorker;
use collagene_core::{
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    review::reads::{
        reference::{self, Operation, Relation},
        session::{self, FileFormat},
    },
};
use collagene_model::{Alphabet, Molecule};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

const PAGE: usize = 100;
const PREVIEW_COLUMNS: usize = 400;
const PREVIEW_TEXT: usize = 2000;

fn format_name(format: FileFormat) -> &'static str {
    match format {
        FileFormat::Sam => "SAM",
        FileFormat::Bam => "BAM",
    }
}
fn extension(format: FileFormat) -> &'static str {
    match format {
        FileFormat::Sam => "sam",
        FileFormat::Bam => "bam",
    }
}
fn clipped(value: &str) -> String {
    value.chars().take(PREVIEW_TEXT).collect()
}
fn cigar(record: &session::RecordPreview) -> String {
    let text = record
        .cigar
        .iter()
        .map(|op| format!("{}{}", op.length, op.operation))
        .collect::<String>();
    if record.truncated.cigar {
        format!("{text} … ({} operations total)", record.cigar_count)
    } else if text.is_empty() {
        "*".into()
    } else {
        text
    }
}
fn record_details(record: &session::RecordPreview) -> String {
    let mut text = format!(
        "Original record BLAKE3: {}\nStored calls (reported order): {}\n",
        record.original_blake3,
        record.sequence.as_deref().unwrap_or("* (sequence absent)")
    );
    if record.truncated.sequence {
        text.push_str(&format!(
            "Call preview clipped; {} reported calls total.\n",
            record.sequence_length.unwrap_or(0)
        ));
    }
    match &record.qualities {
        Some(values) => text.push_str(&format!(
            "Raw Phred scores{}: {:?}\n",
            if record.truncated.qualities {
                " (clipped)"
            } else {
                ""
            },
            values
        )),
        None => text.push_str("Reported qualities absent.\n"),
    }
    text.push_str("Raw scores are uncalibrated source evidence.\n");
    for tag in &record.tags {
        text.push_str(&format!(
            "{}{} [{}]: {}{}\n",
            tag.tag[0] as char,
            tag.tag[1] as char,
            tag.wire_type as char,
            tag.value_summary,
            if tag.truncated { " …" } else { "" }
        ));
    }
    if record.truncated.tags {
        text.push_str(&format!(
            "Tag preview clipped; {} tags total.\n",
            record.tag_count
        ));
    }
    clipped(&text)
}
fn call_preview(report: &reference::Report) -> String {
    let (mut ref_row, mut marker, mut query_row) = (String::new(), String::new(), String::new());
    let mut emitted = 0usize;
    'runs: for run in &report.runs {
        if emitted == PREVIEW_COLUMNS {
            break;
        }
        if matches!(run.operation, Operation::Deletion | Operation::Skip) {
            let count = (run.reference_end - run.reference_start)
                .min((PREVIEW_COLUMNS - emitted) as u64) as usize;
            for _ in 0..count {
                ref_row.push('?');
                query_row.push(if run.operation == Operation::Skip {
                    '~'
                } else {
                    '-'
                });
                marker.push(' ');
            }
            emitted += count;
        } else {
            for column in report
                .columns
                .get(run.column_start..run.column_end)
                .unwrap_or(&[])
            {
                if emitted == PREVIEW_COLUMNS {
                    break 'runs;
                }
                ref_row.push(column.reference_symbol.unwrap_or(
                    if run.operation == Operation::Insertion {
                        '-'
                    } else {
                        '.'
                    },
                ));
                query_row.push(column.query_symbol.unwrap_or('?'));
                marker.push(match column.relation {
                    Relation::CertainAgreement => '|',
                    Relation::CertainDisagreement => '!',
                    Relation::PossibleAgreement => ':',
                    Relation::ReferenceEqualResolved => '=',
                    Relation::InsertionCall => '+',
                    Relation::SoftClippedCall => 's',
                    Relation::Unresolved => '?',
                });
                emitted += 1;
            }
        }
    }
    clipped(&format!(
        "REF {ref_row}\n    {marker}\nQRY {query_row}\n\nPreview capped at400 columns. | certain agreement; ! disagreement; : possible; = reference-derived equality. D/N reference symbols are not copied into this report and show ?; - deletion, ~ skip. Stored call/reference indexes remain in the evidence below."
    ))
}
fn coordinate_preview(report: &reference::Report) -> String {
    let mut text = String::from(
        "Zero-based stored query / implied original query / reference coordinates; raw reported Phred.\n",
    );
    for call in report.columns.iter().take(18) {
        text.push_str(&format!(
            "q{} original{:?} ref{:?} unwrapped{:?} {}→{:?} Q{:?} {:?}{}\n",
            call.query_index,
            call.implied_original_query,
            call.reference_position,
            call.unwrapped_reference_position,
            call.reported_symbol,
            call.query_symbol,
            call.raw_phred,
            call.relation,
            if call.contradicts_cigar {
                " CIGAR contradiction"
            } else {
                ""
            }
        ));
    }
    if report.columns.len() > 18 {
        text.push_str("Coordinate preview clipped; download checked evidence for every column.\n");
    }
    clipped(&text)
}
fn run_preview(report: &reference::Report) -> String {
    clipped(
        &report
            .runs
            .iter()
            .take(20)
            .map(|run| {
                format!(
                    "{:?}: query {}..{}, reference {}..{}, original {:?}\n",
                    run.operation,
                    run.query_start,
                    run.query_end,
                    run.reference_start,
                    run.reference_end,
                    run.implied_original_query_span
                )
            })
            .collect::<String>(),
    )
}
fn reference_parts(report: &reference::Report) -> Vec<(usize, usize)> {
    let mut positions = report
        .columns
        .iter()
        .filter_map(|column| column.reference_position)
        .filter(|&p| p < report.reference_length)
        .collect::<Vec<_>>();
    positions.sort_unstable();
    positions.dedup();
    let mut parts: Vec<(usize, usize)> = vec![];
    for pos in positions {
        if let Some(last) = parts.last_mut().filter(|last| last.1 == pos) {
            last.1 = pos + 1;
        } else {
            parts.push((pos, pos + 1));
        }
    }
    parts
}
#[derive(Clone)]
struct FreshCheck {
    value: Arc<session::Checked>,
    document: usize,
    context: String,
    reference_name: String,
}
impl PartialEq for FreshCheck {
    fn eq(&self, other: &Self) -> bool {
        self.document == other.document
            && self.context == other.context
            && self.reference_name == other.reference_name
            && Arc::ptr_eq(&self.value, &other.value)
    }
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
        return Err("Read-alignment response belongs to different inputs".into());
    }
    match response.result {
        AnalysisResult::Error(error) => Err(error),
        value => Ok(value),
    }
}

/// One retained worker per immutable imported file. The app must remount for a
/// different source. The third document tuple item is document/event/epoch
/// identity; names or equal DNA never substitute for an explicit mapping.
#[component]
pub fn ReadAlignmentPanel(
    source: Arc<Vec<u8>>,
    format: FileFormat,
    name: String,
    documents: Memo<Vec<(usize, String, String)>>,
    get_document: Callback<usize, Option<Molecule>>,
    on_navigate: Callback<(usize, Vec<(usize, usize)>)>,
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
        Err("SAM/BAM imports must be nonempty and at most16MiB. No records were installed.".into())
    });
    let source = StoredValue::new(source);
    let filename = StoredValue::new(name.clone());
    let prepared = RwSignal::new(None::<Arc<session::Prepared>>);
    let page = RwSignal::new(None::<Arc<session::Page>>);
    let selected = RwSignal::new(None::<session::RecordPreview>);
    let selected_document = RwSignal::new(None::<usize>);
    let header_reference = RwSignal::new(String::new());
    let checked = RwSignal::new(None::<FreshCheck>);
    let exported = RwSignal::new(None::<Arc<session::Exported>>);
    let export_format = RwSignal::new(match format {
        FileFormat::Sam => FileFormat::Bam,
        FileFormat::Bam => FileFormat::Sam,
    });
    let busy = RwSignal::new(false);
    let message = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let request_id = StoredValue::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let selected_context = Memo::new(move |_| {
        documents.with(|docs| {
            selected_document.get().and_then(|selected| {
                docs.iter()
                    .find(|(index, _, _)| *index == selected)
                    .map(|(index, _, identity)| (*index, identity.clone()))
            })
        })
    });
    let close_session = Callback::new(move |reason: String| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(client) = worker.get_value() {
            client.close();
        }
        worker.set_value(None);
        busy.set(false);
        prepared.set(None);
        page.set(None);
        selected.set(None);
        checked.set(None);
        exported.set(None);
        message.set(reason);
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(client) = worker.get_value() {
            client.close();
        }
        worker.set_value(None);
    });
    // Inputs changing while a check is in flight revoke the displayed/future
    // check. The worker operation may finish, but its result cannot be published.
    Effect::new(move |_| {
        let _ = selected_context.get();
        let _ = header_reference.get();
        let _ = selected.with(|record| {
            record
                .as_ref()
                .map(|record| (record.index, record.original_blake3.clone()))
        });
        checked.set(None);
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
        let (Some(client), Some(summary)) = (worker.get_value(), prepared.get_untracked()) else {
            return;
        };
        let hash = source_hash.get_value();
        let revision = match session::page_revision(&hash, start, PAGE) {
            Ok(value) => value,
            Err(error) => {
                message.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = page_alive.clone();
        busy.set(true);
        selected.set(None);
        checked.set(None);
        message.set("Loading source records…".into());
        spawn_local(async move {
            let result = ask(
                client,
                id,
                revision.clone(),
                AnalysisJob::ReadAlignmentPage {
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
                Ok(AnalysisResult::ReadAlignmentPage(value)) => {
                    let valid = value.validate().map_err(|e| e.to_string()).and_then(|_| {
                        if value.revision == revision
                            && value.source_hash == hash
                            && value.start == start
                            && value.count == PAGE
                            && value.total == summary.summary.records
                        {
                            Ok(())
                        } else {
                            Err("Read page has different source/offset binding".into())
                        }
                    });
                    match valid {
                        Ok(()) => {
                            page.set(Some(Arc::new(value)));
                            message.set(String::new());
                        }
                        Err(error) => message.set(error),
                    }
                }
                Ok(_) => message.set("Unexpected record page response".into()),
                Err(error) => message.set(error),
            }
        });
    });
    let import_alive = alive.clone();
    let prepare = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let revision = match prepare_binding.get_value() {
            Ok(value) => value,
            Err(error) => {
                message.set(error);
                return;
            }
        };
        close_session.run(String::new());
        header_reference.set(String::new());
        let client = match CachedWorker::new() {
            Ok(value) => Rc::new(value),
            Err(error) => {
                message.set(error);
                return;
            }
        };
        worker.set_value(Some(client.clone()));
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = import_alive.clone();
        let hash = source_hash.get_value();
        // The immutable source has passed the fixed byte bound before this copy.
        let bytes = source.with_value(|source| source.as_ref().clone());
        busy.set(true);
        message.set("Parsing the complete SAM/BAM source…".into());
        spawn_local(async move {
            let result = ask(
                client,
                id,
                revision.clone(),
                AnalysisJob::PrepareReadAlignments { bytes, format },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::ReadAlignmentsPrepared(value)) => {
                    let valid = value.validate().map_err(|e| e.to_string()).and_then(|_| {
                        if value.revision == revision && value.summary.source_blake3 == hash {
                            Ok(())
                        } else {
                            Err("Prepared records have a different source binding".into())
                        }
                    });
                    match valid {
                        Ok(()) => {
                            prepared.set(Some(Arc::new(value)));
                            message.set("Complete file parsed; imported alignment claims remain unverified.".into());
                            load_page.run(0);
                        }
                        Err(error) => close_session.run(error),
                    }
                }
                Ok(_) => close_session.run("Unexpected source preparation response".into()),
                Err(error) => close_session.run(error),
            }
        });
    });
    let check_alive = alive.clone();
    let check_record = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let (Some(client), Some(record), Some((document, context))) = (
            worker.get_value(),
            selected.get_untracked(),
            selected_context.get_untracked(),
        ) else {
            message.set("Choose a source record and an actual DNA document.".into());
            return;
        };
        let reference_name = header_reference.get_untracked();
        if reference_name.is_empty()
            || !prepared.with_untracked(|value| {
                value.as_ref().is_some_and(|value| {
                    value
                        .summary
                        .references
                        .iter()
                        .any(|r| r.name == reference_name)
                })
            })
        {
            message.set("Explicitly choose the reported header reference to map.".into());
            return;
        }
        let Some(reference) = get_document.run(document) else {
            message.set(
                "The chosen document is unavailable or exceeds reference-check limits.".into(),
            );
            return;
        };
        if reference.sequence.alphabet() != Alphabet::Dna {
            message.set("Choose an actual DNA document for this reference mapping.".into());
            return;
        }
        if let Err(error) = session::reference_preflight(&reference) {
            message.set(error.to_string());
            return;
        }
        let reference_revision = reference.revision();
        let options = reference::Options {
            reference_name: reference_name.clone(),
            expected_reference_revision: reference_revision.clone(),
            ..Default::default()
        };
        let hash = source_hash.get_value();
        let revision = match session::check_revision_for_revision(
            &hash,
            record.index,
            &reference_revision,
            &options,
        ) {
            Ok(value) => value,
            Err(error) => {
                message.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = check_alive.clone();
        let expected_options = options.clone();
        busy.set(true);
        checked.set(None);
        message
            .set("Checking reported calls against the explicitly selected DNA reference…".into());
        spawn_local(async move {
            let result = ask(
                client,
                id,
                revision.clone(),
                AnalysisJob::CheckReadAlignment {
                    source_hash: hash.clone(),
                    index: record.index,
                    reference,
                    options,
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            busy.set(false);
            let current = selected_context.get_untracked() == Some((document, context.clone()))
                && header_reference.get_untracked() == reference_name
                && selected.with_untracked(|current| {
                    current.as_ref().is_some_and(|current| {
                        current.index == record.index
                            && current.original_blake3 == record.original_blake3
                    })
                });
            if !current {
                message.set("Reference/record context changed; the old check was discarded. Check the current mapping again.".into());
                return;
            }
            match result {
                Ok(AnalysisResult::ReadAlignmentChecked(value)) => match value.validate_bindings(
                    &revision,
                    &hash,
                    record.index,
                    &record.original_blake3,
                    &reference_revision,
                    &expected_options,
                ) {
                    Ok(()) => {
                        checked.set(Some(FreshCheck {
                            value: Arc::new(value),
                            document,
                            context,
                            reference_name,
                        }));
                        message.set(String::new());
                    }
                    Err(error) => message.set(error.to_string()),
                },
                Ok(_) => message.set("Unexpected reference-check response".into()),
                Err(error) => message.set(error),
            }
        });
    });
    let export_alive = alive.clone();
    let convert = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let Some(client) = worker.get_value() else {
            return;
        };
        let hash = source_hash.get_value();
        let format = export_format.get_untracked();
        let cap = session::MAX_EXPORT_BYTES;
        let revision = match session::export_revision(&hash, format, cap) {
            Ok(value) => value,
            Err(error) => {
                message.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = export_alive.clone();
        busy.set(true);
        exported.set(None);
        message.set("Converting the sealed source records…".into());
        spawn_local(async move {
            let result = ask(
                client,
                id,
                revision.clone(),
                AnalysisJob::ExportReadAlignments {
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
            if export_format.get_untracked() != format {
                message.set("Export format changed; the previous conversion was discarded.".into());
                return;
            }
            match result {
                Ok(AnalysisResult::ReadAlignmentsExported(value)) => {
                    let valid = value.validate().map_err(|e| e.to_string()).and_then(|_| {
                        if value.revision == revision
                            && value.source_hash == hash
                            && value.format == format
                            && value.bytes.len() <= cap
                        {
                            Ok(())
                        } else {
                            Err("Converted export has a different source/format binding".into())
                        }
                    });
                    match valid {
                        Ok(()) => {
                            exported.set(Some(Arc::new(value)));
                            message.set(
                                "Conversion prepared. Review its notices before downloading."
                                    .into(),
                            );
                        }
                        Err(error) => message.set(error),
                    }
                }
                Ok(_) => message.set("Unexpected converted export response".into()),
                Err(error) => message.set(error),
            }
        });
    });
    let dispose_alive = alive.clone();
    let dispose = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let Some(client) = worker.get_value() else {
            return;
        };
        let hash = source_hash.get_value();
        let revision = match session::dispose_revision(&hash) {
            Ok(value) => value,
            Err(error) => {
                message.set(error.to_string());
                return;
            }
        };
        let ticket = generation.get_untracked();
        let id = next_id();
        let alive = dispose_alive.clone();
        busy.set(true);
        spawn_local(async move {
            let result = ask(
                client,
                id,
                revision.clone(),
                AnalysisJob::DisposeReadAlignments {
                    source_hash: hash.clone(),
                },
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != ticket {
                return;
            }
            match result {
                Ok(AnalysisResult::ReadAlignmentsDisposed(value))
                    if value.validate().is_ok()
                        && value.revision == revision
                        && value.source_hash == hash =>
                {
                    close_session
                        .run("Parsed records released; original bytes remain available.".into())
                }
                Ok(_) => close_session.run("The parsed-record session was closed.".into()),
                Err(error) => close_session.run(error),
            }
        });
    });
    let current_check = Memo::new(move |_| {
        checked.get().filter(|check| {
            selected_context.get() == Some((check.document, check.context.clone()))
                && header_reference.get() == check.reference_name
                && selected.with(|record| {
                    record.as_ref().is_some_and(|record| {
                        record.index == check.value.check.record_index
                            && record.original_blake3 == check.value.check.record_blake3
                    })
                })
        })
    });
    prepare.run(());
    view! {
        <section class="read-alignments-panel" aria-label="SAM/BAM read alignments">
            <h2>{format!("{} · {} read alignments",name,format_name(format))}</h2>
            <p class="helper">"These are imported alignment assertions. Reference names, flags, MAPQ and CIGAR do not establish unique mapping, coverage, consensus or variants. Choose an actual DNA document and an explicit header-name mapping to check a record."</p>
            <div class="read-alignment-actions">
                <button class="button secondary" disabled=!bounded on:click=move|_|{let result=source.with_value(|source|crate::io::download::ExportArtifact::prepare(&filename.get_value(),"application/octet-stream",source.as_slice()));crate::io::download::complete_download(result,move|result|{if message.try_get_untracked().is_some(){if let Err(error)=result{message.set(error);}}});}>"Download original SAM/BAM"</button>
                <button class="button secondary" disabled=move||busy.get() on:click=move|_|prepare.run(())>"Reopen source records"</button>
                <button class="button secondary" disabled=move||!busy.get() on:click=move|_|close_session.run("Task cancelled and worker closed. Reopen the source to continue.".into())>"Cancel and close session"</button>
                <button class="text-tool" disabled=move||busy.get()||prepared.with(Option::is_none) on:click=move|_|dispose.run(())>"Release parsed records"</button>
            </div>
            <p class="read-alignment-status" role="status">{move||message.get()}</p>
            {move||prepared.get().map(|value|{let s=&value.summary;view!{
                <div class="read-alignment-summary"><p>{format!("Complete {} source: {} bytes · {} records · BLAKE3 {}",s.format,s.source_bytes,s.records,s.source_blake3)}</p><p>{format!("Reported flags: {} unmapped · {} secondary · {} supplementary. These flags are source claims.",s.flag_unmapped,s.flag_secondary,s.flag_supplementary)}</p><p>{s.interpretation.clone()}</p>{s.warnings.iter().take(20).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
                    <details><summary>{format!("Reported reference dictionary · {} entries (not verified build identity)",s.references.len())}</summary>{s.references.iter().take(100).map(|reference|view!{<p>{format!("{} · {} bases · {} (header claim)",reference.name,reference.length,if reference.circular{"circular"}else{"linear"})}</p>}).collect_view()}{(s.references.len()>100).then(||view!{<p>"Header preview shows the first100 references; the mapping menu retains the accepted dictionary."</p>})}</details>
                </div>
            }})}
            <div class="read-alignment-pages"><button class="text-tool" disabled=move||busy.get()||page.with(|p|p.as_ref().is_none_or(|p|p.start==0)) on:click=move|_|if let Some(p)=page.get_untracked(){load_page.run(p.start.saturating_sub(PAGE));}>"Previous100 records"</button><span>{move||page.with(|p|p.as_ref().map(|p|format!("{}–{} of {} reported records",if p.records.is_empty(){0}else{p.start+1},p.start+p.records.len(),p.total)).unwrap_or_default())}</span><button class="text-tool" disabled=move||busy.get()||page.with(|p|p.as_ref().is_none_or(|p|p.start+p.records.len()>=p.total)) on:click=move|_|if let Some(p)=page.get_untracked(){load_page.run(p.start+PAGE);}>"Next100 records"</button></div>
            <div class="read-alignment-table-wrap"><table class="read-alignment-records"><thead><tr><th>"Record"</th><th>"Reported name"</th><th>"Reference / position"</th><th>"Flags / MAPQ"</th><th>"Reported CIGAR"</th></tr></thead><tbody>
                <For each=move||page.with(|p|p.as_ref().map(|p|p.records.clone()).unwrap_or_default()) key=|record|record.index children=move|record:session::RecordPreview|{let index=record.index;let chosen=record.clone();view!{<tr data-record-index=index class:active=move||selected.with(|r|r.as_ref().is_some_and(|r|r.index==index))><td><button class="text-tool" disabled=move||busy.get() on:click=move|_|{selected.set(Some(chosen.clone()));checked.set(None);}>{format!("Select {}",index+1)}</button></td><td>{format!("{}{}",record.name.as_deref().unwrap_or("*"),if record.truncated.name{" …"}else{""})}</td><td>{format!("{}{} · {}",record.reference.as_deref().unwrap_or("*"),if record.truncated.reference{" …"}else{""},record.start.map(|p|(p as u64+1).to_string()).unwrap_or_else(||"unavailable".into()))}</td><td>{format!("0x{:04x} / {}",record.flags,record.mapping_quality.map(|q|q.to_string()).unwrap_or_else(||"unavailable".into()))}</td><td>{cigar(&record)}</td></tr>}}/>
            </tbody></table></div>
            {move||selected.get().map(|record|view!{<details open><summary>{format!("Reported record {} · bounded source preview",record.index+1)}</summary><pre>{record_details(&record)}</pre></details>})}
            <div class="read-alignment-reference-mapping">
                <label>"Actual DNA document"<select aria-label="Read alignment DNA reference" prop:value=move||selected_document.get().map(|index|index.to_string()).unwrap_or_default() on:change=move|event|selected_document.set(event_target_value(&event).parse::<usize>().ok())><option value="">"Choose an actual document…"</option>{move||documents.get().into_iter().map(|(index,name,_)|view!{<option value=index.to_string()>{name}</option>}).collect_view()}</select></label>
                <label>"Reported header reference"<select aria-label="Read alignment header reference" prop:value=move||header_reference.get() on:change=move|event|header_reference.set(event_target_value(&event))><option value="">"Choose a header name explicitly…"</option>{move||prepared.with(|value|value.as_ref().map(|value|value.summary.references.iter().map(|reference|view!{<option value=reference.name.clone()>{format!("{} · {} bases (header)",reference.name,reference.length)}</option>}).collect_view()).into_view())}</select></label>
                <p class="helper">"A matching name does not authorize this mapping. Changing the document event, header mapping or selected record revokes the current check."</p>
                <button class="button primary" disabled=move||busy.get()||selected.with(Option::is_none)||selected_context.with(Option::is_none)||header_reference.with(String::is_empty)||prepared.with(Option::is_none) on:click=move|_|check_record.run(())>"Check record against DNA reference"</button>
            </div>
            {move||current_check.get().map(|fresh|{let report=&fresh.value.check.report;let stats=&report.summary;let preview=call_preview(report);let coords=coordinate_preview(report);let runs=run_preview(report);let bound=report.reference_bound;let parts=reference_parts(report);let nav=fresh.clone();let evidence=fresh.clone();view!{
                <div class="read-alignment-check"><h3>"Fresh source-bound reference check"</h3><p>{format!("Reference mapped: {} · complete: {} · orientation: {} · source record {}",bound,report.complete,if report.reverse{"reverse (reported orientation)"}else{"forward (reported orientation)"},fresh.value.check.record_index+1)}</p><code>{fresh.context.clone()}</code><p>{format!("Certain observed agreements {} · certain disagreements {} · possible agreements {} · reference-derived '=' calls {} · unresolved {} · CIGAR contradictions {}",stats.certain_agreements,stats.certain_disagreements,stats.possible_agreements,stats.reference_equal_resolved,stats.unresolved_calls,stats.cigar_contradictions)}</p><p>{format!("Insertions {} · deletions {} · skipped reference bases {} · soft-clipped calls {} · hard-clipped absent bases {}",stats.inserted_calls,stats.deleted_reference_bases,stats.skipped_reference_bases,stats.soft_clipped_calls,stats.hard_clipped_absent_bases)}</p><p>{format!("Incomplete reasons: {:?}",report.incomplete)}</p>{report.warnings.iter().take(20).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}<pre class="read-alignment-gapped-preview">{preview}</pre><details><summary>"Stored/original call coordinates and reported operations"</summary><pre>{coords}</pre><pre>{runs}</pre></details>
                    <button class="text-tool" disabled=parts.is_empty()||!bound on:click=move|_|{if current_check.with_untracked(|current|current.as_ref().is_some_and(|current|current.context==nav.context&&current.value.revision==nav.value.revision)){on_navigate.run((nav.document,parts.clone()));}}>"Show checked call coordinates"</button>
                    <button class="text-tool" on:click=move|_|{if !current_check.with_untracked(|current|current.as_ref().is_some_and(|current|current.context==evidence.context&&current.value.revision==evidence.value.revision)){return;}let value=serde_json::json!({"interpretation":"Saved source-bound checking evidence only; recompute for changed source/reference/event. No coverage, variants, consensus or unique-alignment claim.","document_context":evidence.context,"checked":evidence.value.as_ref()});match serde_json::to_vec(&value){Ok(bytes)if bytes.len()<=session::MAX_EXPORT_BYTES=>crate::io::download::download_with_completion("read-alignment-check.json","application/json",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(error)=result {message.set(error)}}),Ok(_)=>message.set("Checked evidence exceeds the download byte limit.".into()),Err(error)=>message.set(error.to_string())}}>"Download checked source evidence"</button>
                </div>
            }})}
            <div class="read-alignment-export"><label>"Converted exchange format"<select aria-label="Read alignment export format" prop:value=move||format_name(export_format.get()) on:change=move|event|{export_format.set(if event_target_value(&event)=="SAM"{FileFormat::Sam}else{FileFormat::Bam});exported.set(None);}><option value="SAM">"SAM"</option><option value="BAM">"BAM / BGZF"</option></select></label><button class="button secondary" disabled=move||busy.get()||prepared.with(Option::is_none) on:click=move|_|convert.run(())>"Prepare converted export"</button>{move||exported.get().map(|value|{let download_value=value.clone();view!{<div class="read-alignment-export-notices"><p>"Converted output is not the original byte artifact. Unsupported values or unavailable calls/qualities are never invented to force conversion."</p>{value.notices.iter().map(|notice|view!{<p>{notice.clone()}</p>}).collect_view()}<button class="text-tool" on:click=move|_|{if download_value.source_hash!=source_hash.get_value(){return;}crate::io::download::download_with_completion(&format!("read-alignments.converted.{}",extension(download_value.format)),"application/octet-stream",&download_value.bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(error)=result {message.set(error);}})}>"Download converted file"</button></div>}})}</div>
        </section>
    }
}
