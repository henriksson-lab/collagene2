//! Reference-free read assembly with exact source-event binding and fresh replay before use.
use crate::platform::worker::CachedWorker;
use crate::viewer::map::renderer::ResizeListener;
use collagene_core::{
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    review::contigs::{
        self as contig, ConsensusState, Contig, ContigOptions, ContigReport, ContigTrimOptions,
        OverlapScoring, ReadInput,
    },
};
use collagene_model::Molecule;
use leptos::{html, prelude::*};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadInfo {
    pub document_id: String,
    pub event_id: String,
    pub revision: String,
    pub name: String,
    pub length: usize,
    pub eligible: bool,
    pub quality_label: String,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadBinding {
    pub document_id: String,
    pub event_id: String,
    pub revision: String,
}
#[derive(Clone)]
pub struct ArchiveRequest {
    pub report: Arc<ContigReport>,
    pub bindings: Vec<ReadBinding>,
    pub contig_id: Option<usize>,
    pub allow_uncertain: bool,
    pub epoch: u64,
}
#[derive(Clone)]
pub struct ReadNavigation {
    pub binding: ReadBinding,
    pub call_index: usize,
    pub epoch: u64,
}
#[derive(Clone)]
struct BoundAssembly {
    report: Arc<ContigReport>,
    bindings: Vec<ReadBinding>,
    epoch: u64,
}
fn whole<T: std::str::FromStr>(s: &str, label: &str) -> Result<T, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {label}."))
}
fn percent(s: &str) -> Result<u16, String> {
    let value = s
        .trim()
        .parse::<f64>()
        .map_err(|_| "Enter an error percentage in 0–100.".to_string())?;
    if !value.is_finite()
        || !(0.0..=100.0).contains(&value)
        || (value * 100.0 - (value * 100.0).round()).abs() > 0.000001
    {
        return Err("Error percentage accepts at most two decimal places in 0–100.".into());
    }
    Ok((value * 100.0).round() as u16)
}
fn selected_bindings(infos: &[ReadInfo], selected: &[String]) -> Result<Vec<ReadBinding>, String> {
    selected
        .iter()
        .map(|id| {
            let info = infos
                .iter()
                .find(|i| &i.document_id == id)
                .ok_or("A selected source document is no longer open.")?;
            if !info.eligible {
                return Err(format!("{} is not a linear DNA read.", info.name));
            }
            Ok(ReadBinding {
                document_id: info.document_id.clone(),
                event_id: info.event_id.clone(),
                revision: info.revision.clone(),
            })
        })
        .collect()
}
fn inputs(
    get_read: Callback<String, Option<Molecule>>,
    bindings: &[ReadBinding],
) -> Result<Vec<ReadInput>, String> {
    bindings
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let molecule = get_read
                .run(b.document_id.clone())
                .ok_or("A selected source document is no longer open.")?;
            if molecule.revision() != b.revision {
                return Err("A selected source changed. Assemble the current reads again.".into());
            }
            Ok(ReadInput {
                id: format!("read_{}", i + 1),
                molecule,
            })
        })
        .collect()
}
struct Json(Vec<u8>);
impl std::io::Write for Json {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > 32 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "Contig JSON exceeds the browser's 32 MiB export bound",
            ));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn state_label(state: &ConsensusState) -> &'static str {
    match state {
        ConsensusState::SupportedBase => "Supported base",
        ConsensusState::SupportedGap => "Supported deletion/gap",
        ConsensusState::ConflictingCalls => "Conflicting calls",
        ConsensusState::GapConflict => "Base/gap conflict",
        ConsensusState::InsufficientSupport => "Insufficient support",
        ConsensusState::LowQualityOrAmbiguousOnly => "Low quality or ambiguity only",
    }
}
fn current_contig(data: &BoundAssembly, index: usize) -> Option<&Contig> {
    data.report.contigs.get(index)
}
fn needs_ack(data: &BoundAssembly, c: &Contig) -> bool {
    !data.report.complete
        || !c.complete
        || !c.layout_unique
        || c.sequence.bytes().any(|b| !b"ACGT".contains(&b))
}

#[component]
fn LayoutCanvas(
    data: RwSignal<Option<Arc<BoundAssembly>>>,
    contig_index: RwSignal<usize>,
    start: RwSignal<usize>,
    selected_column: RwSignal<Option<usize>>,
    row_start: RwSignal<usize>,
    expanded: RwSignal<bool>,
) -> impl IntoView {
    let canvas = NodeRef::<html::Canvas>::new();
    let container = NodeRef::<html::Div>::new();
    let resize = StoredValue::new_local(None::<ResizeListener>);
    let resized = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, Ordering::Relaxed));
    Effect::new(move |_| {
        if let Some(element) = container.get() {
            if resize.with_value(|r| r.is_none()) {
                resize.set_value(
                    ResizeListener::new(element.unchecked_ref(), move || {
                        let _ = resized.try_update(|v| *v = v.wrapping_add(1));
                    })
                    .ok(),
                );
            }
        }
    });
    Effect::new(move |_| {
        let _ = (resized.get(), expanded.get());
        let first = start.get();
        let row_first = row_start.get();
        let selected = selected_column.get();
        let chosen = contig_index.get();
        let bound = data.get();
        let Some(c) = canvas.get() else {
            return;
        };
        let Some(bound) = bound else {
            return;
        };
        let Some(contig) = current_contig(&bound, chosen) else {
            return;
        };
        let width = c.get_bounding_client_rect().width().max(240.);
        let rows = contig
            .placements
            .iter()
            .skip(row_first)
            .take(12)
            .collect::<Vec<_>>();
        let height = 100. + rows.len() as f64 * 22.;
        let ratio = web_sys::window()
            .map(|w| w.device_pixel_ratio())
            .unwrap_or(1.)
            .clamp(1., 2.);
        c.set_width((width * ratio) as u32);
        c.set_height((height * ratio) as u32);
        let raw: &web_sys::HtmlCanvasElement = &c;
        let _ = web_sys::HtmlElement::style(raw.unchecked_ref::<web_sys::HtmlElement>())
            .set_property("height", &format!("{height}px"));
        let Ok(Some(ctx)) = c.get_context("2d") else {
            return;
        };
        let Ok(ctx) = ctx.dyn_into::<web_sys::CanvasRenderingContext2d>() else {
            return;
        };
        let _ = ctx.set_transform(ratio, 0., 0., ratio, 0., 0.);
        ctx.set_fill_style_str("#ffffff");
        ctx.fill_rect(0., 0., width, height);
        ctx.set_font("12px monospace");
        ctx.set_text_baseline("middle");
        let visible = ((width - 116.) / 9.).floor().max(1.) as usize;
        let last = (first + visible.min(160)).min(contig.columns.len());
        ctx.set_fill_style_str("#33594c");
        let _ = ctx.fill_text("Consensus", 4., 43.);
        for (j, column) in contig.columns.iter().enumerate().take(last).skip(first) {
            let x = 116. + (j - first) as f64 * 9.;
            if selected == Some(j) {
                ctx.set_fill_style_str("#e0eecd");
                ctx.fill_rect(x - 1., 22., 10., height - 26.);
            }
            if j == first || j % 10 == 0 {
                ctx.set_fill_style_str("#75887d");
                let _ = ctx.fill_text(&(j + 1).to_string(), x, 12.);
            }
            ctx.set_fill_style_str(if matches!(column.state, ConsensusState::SupportedBase) {
                "#286b50"
            } else {
                "#ad6d26"
            });
            let _ = ctx.fill_text(&column.base.unwrap_or('-').to_string(), x, 43.);
            ctx.set_fill_style_str("#809087");
            ctx.set_font("9px monospace");
            let _ = ctx.fill_text(&column.observations.len().to_string(), x, 61.);
            ctx.set_font("12px monospace");
        }
        for (row, p) in rows.iter().enumerate() {
            let y = 86. + row as f64 * 22.;
            ctx.set_fill_style_str("#365549");
            let name = bound
                .report
                .reads
                .get(p.read_index)
                .map(|r| r.id.as_str())
                .unwrap_or("read");
            let _ = ctx.fill_text(
                &format!("{} {}", name, if p.reverse { "←" } else { "→" }),
                4.,
                y,
            );
            for column in contig.columns.iter().take(last).skip(first) {
                if let Some(call) = column
                    .observations
                    .iter()
                    .find(|o| o.read_index == p.read_index)
                {
                    let x = 116. + (column.layout_index - first) as f64 * 9.;
                    ctx.set_fill_style_str(if call.eligible { "#2d6352" } else { "#b59470" });
                    let _ = ctx.fill_text(&call.base.to_string(), x, y);
                } else if column.gap_reads.contains(&p.read_index) {
                    let x = 116. + (column.layout_index - first) as f64 * 9.;
                    ctx.set_fill_style_str("#8d7963");
                    let _ = ctx.fill_text("-", x, y);
                }
            }
        }
    });
    view! {<div node_ref=container class="contig-layout-container"><canvas node_ref=canvas class="contig-layout-canvas" tabindex="0" aria-label="Contig layout. Click a column to inspect exact calls. Arrow keys move the selected column." on:click=move|e|{let Some(c)=canvas.get()else{return;};let x=e.client_x()as f64-c.get_bounding_client_rect().left();if x<116.{return;}let index=start.get_untracked()+((x-116.)/9.).floor()as usize;if data.with_untracked(|b|b.as_ref().and_then(|b|current_contig(b,contig_index.get_untracked())).is_some_and(|c|index<c.columns.len())){selected_column.set(Some(index));}} on:keydown=move|e|{let count=data.with_untracked(|b|b.as_ref().and_then(|b|current_contig(b,contig_index.get_untracked())).map_or(0,|c|c.columns.len()));if count==0{return;}let at=selected_column.get_untracked().unwrap_or(start.get_untracked());let next=match e.key().as_str(){"ArrowRight"=>Some((at+1).min(count-1)),"ArrowLeft"=>Some(at.saturating_sub(1)),"Home"=>Some(0),"End"=>Some(count-1),_=>None};if let Some(next)=next{e.prevent_default();selected_column.set(Some(next));if next<start.get_untracked()||next>=start.get_untracked()+10{start.set(next.saturating_sub(5));}}}/></div>}
}

#[component]
pub fn ContigPanel(
    read_options: Memo<Vec<ReadInfo>>,
    get_read: Callback<String, Option<Molecule>>,
    epoch: RwSignal<u64>,
    on_record: Callback<ArchiveRequest, Result<(), String>>,
    on_navigate: Callback<ReadNavigation>,
) -> impl IntoView {
    let selected = RwSignal::new(Vec::<String>::new());
    let source_page = RwSignal::new(0usize);
    let trim_windows = RwSignal::new(false);
    let threshold = RwSignal::new("20".to_string());
    let trim_window = RwSignal::new("5".to_string());
    let passing = RwSignal::new("4".to_string());
    let minimum = RwSignal::new("20".to_string());
    let unscored = RwSignal::new(false);
    let min_overlap = RwSignal::new("20".to_string());
    let error = RwSignal::new("10".to_string());
    let both = RwSignal::new(true);
    let support = RwSignal::new("1".to_string());
    let both_support = RwSignal::new(false);
    let match_score = RwSignal::new("2".to_string());
    let ambiguous_score = RwSignal::new("0".to_string());
    let mismatch_score = RwSignal::new("-3".to_string());
    let gap_open = RwSignal::new("-4".to_string());
    let gap_extend = RwSignal::new("-1".to_string());
    let pair_cells = RwSignal::new("4000000".to_string());
    let total_cells = RwSignal::new("128000000".to_string());
    let max_edges = RwSignal::new("4096".to_string());
    let max_columns = RwSignal::new("100000".to_string());
    let options = Memo::new(move |_| -> Result<ContigOptions, String> {
        Ok(ContigOptions {
            trim: ContigTrimOptions {
                quality_windows: trim_windows.get(),
                raw_score_threshold: whole(&threshold.get(), "raw confidence threshold")?,
                window_size: whole(&trim_window.get(), "trim window")?,
                minimum_passing: whole(&passing.get(), "passing calls per window")?,
                min_retained_bases: whole(&minimum.get(), "minimum retained bases")?,
                allow_unscored: unscored.get(),
            },
            scoring: OverlapScoring {
                match_score: whole(&match_score.get(), "match score")?,
                ambiguous_score: whole(&ambiguous_score.get(), "ambiguity score")?,
                mismatch_score: whole(&mismatch_score.get(), "mismatch score")?,
                gap_open: whole(&gap_open.get(), "gap opening score")?,
                gap_extend: whole(&gap_extend.get(), "gap extension score")?,
            },
            min_overlap: whole(&min_overlap.get(), "minimum overlap")?,
            max_error_bp: percent(&error.get())?,
            both_strands: both.get(),
            min_support: whole(&support.get(), "minimum call support")?,
            require_both_strands: both_support.get(),
            max_pair_cells: whole(&pair_cells.get(), "pair cell budget")?,
            max_total_cells: whole(&total_cells.get(), "total cell budget")?,
            max_edges: whole(&max_edges.get(), "edge limit")?,
            max_columns: whole(&max_columns.get(), "column limit")?,
        })
    });
    let bindings = Memo::new(move |_| selected_bindings(&read_options.get(), &selected.get()));
    let data = RwSignal::new(None::<Arc<BoundAssembly>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let contig_index = RwSignal::new(0usize);
    let column_start = RwSignal::new(0usize);
    let column_page = RwSignal::new(0usize);
    let row_start = RwSignal::new(0usize);
    let selected_column = RwSignal::new(None::<usize>);
    let edge_page = RwSignal::new(0usize);
    let selected_edge = RwSignal::new(None::<usize>);
    let expanded = RwSignal::new(false);
    let acknowledge = RwSignal::new(false);
    Effect::new(move |_| {
        let _ = (bindings.get(), options.get(), epoch.get());
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        data.set(None);
        busy.set(false);
        acknowledge.set(false);
        contig_index.set(0);
        column_start.set(0);
        column_page.set(0);
        row_start.set(0);
        selected_column.set(None);
        selected_edge.set(None);
        status.set(String::new());
    });
    let cancel = Callback::new(move |()| {
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        status.set("Assembly cancelled. No consensus was installed.".into());
    });
    let guard = Callback::new(move |bound: Arc<BoundAssembly>| -> Result<(), String> {
        if epoch.get_untracked() != bound.epoch
            || bindings.get_untracked() != Ok(bound.bindings.clone())
            || options.get_untracked() != Ok(bound.report.options.clone())
        {
            return Err("Read identities, source events or options changed. Assemble the current reads again.".into());
        }
        Ok(())
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            let opts = match options.get_untracked() {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let refs = match bindings.get_untracked() {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if refs.is_empty() {
                status.set("Select linear DNA read documents before assembling.".into());
                return;
            }
            let reads = match inputs(get_read, &refs) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let batch = match contig::batch_revision(&reads, &opts) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            if let Some(w) = worker.get_value() {
                w.close();
            }
            let w = match CachedWorker::new() {
                Ok(v) => Rc::new(v),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            generation.update(|n| *n = n.wrapping_add(1));
            let token = generation.get_untracked();
            let expected_epoch = epoch.get_untracked();
            data.set(None);
            busy.set(true);
            status.set("Comparing read ends and reconstructing source-supported layouts…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: batch,
                        job: AnalysisJob::AssembleContigs {
                            reads: reads.clone(),
                            options: opts.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || epoch.get_untracked() != expected_epoch
                    || bindings.get_untracked() != Ok(refs.clone())
                    || options.get_untracked() != Ok(opts.clone())
                {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                match response {
                    Ok(response) => match response.result {
                        AnalysisResult::Contigs(report) => {
                            match report.validate_bindings_for(&reads, &opts) {
                                Ok(()) => {
                                    status.set(if report.complete{"Declared overlap search completed. Consensus is read evidence, not experimental identity."}else{"Assembly reached a declared limit. Missing edges or contigs do not establish absence."}.into());
                                    data.set(Some(Arc::new(BoundAssembly {
                                        report: Arc::new(report),
                                        bindings: refs,
                                        epoch: expected_epoch,
                                    })));
                                }
                                Err(e) => status.set(e.to_string()),
                            }
                        }
                        AnalysisResult::Error(e) => status.set(e),
                        _ => status.set("Unexpected assembly response.".into()),
                    },
                    Err(e) => status.set(e),
                }
            });
        }
    };
    // Every released DNA document comes from a second, source-bound reconstruction.
    let release = {
        let alive = alive.clone();
        Callback::new(move |kind: &'static str| {
            let Some(bound) = data.get_untracked() else {
                return;
            };
            if let Err(e) = guard.run(bound.clone()) {
                status.set(e);
                return;
            }
            let chosen = contig_index.get_untracked();
            let Some(old) = bound.report.contigs.get(chosen) else {
                return;
            };
            let Some(old_molecule) = old.molecule.as_ref() else {
                status
                    .set("This layout has no independently releasable consensus document.".into());
                return;
            };
            if needs_ack(&bound, old) && !acknowledge.get_untracked() {
                status.set("Acknowledge the unresolved or incomplete evidence before using this consensus.".into());
                return;
            }
            if kind == "open" {
                match on_record.run(ArchiveRequest { report:bound.report.clone(),bindings:bound.bindings.clone(),contig_id:Some(old.id),allow_uncertain:acknowledge.get_untracked(),epoch:bound.epoch }) {
                    Ok(())=>status.set("Review the source-bound archival transaction before recording this consensus.".into()),
                    Err(error)=>status.set(error),
                }
                return;
            }
            let expected = old_molecule.revision();
            let expected_sequence = old.sequence.clone();
            let old_read_indices = old
                .placements
                .iter()
                .map(|p| p.read_index)
                .collect::<Vec<_>>();
            let contig_id = old.id;
            let reads = match inputs(get_read, &bound.bindings) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let opts = bound.report.options.clone();
            let revision = match contig::batch_revision(&reads, &opts) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(v) => Rc::new(v),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            generation.update(|n| *n = n.wrapping_add(1));
            let token = generation.get_untracked();
            busy.set(true);
            status.set(
                "Replaying this assembly against its exact current source events before release…"
                    .into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision,
                        job: AnalysisJob::AssembleContigs {
                            reads: reads.clone(),
                            options: opts.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                if let Err(e) = guard.run(bound.clone()) {
                    status.set(e);
                    return;
                }
                let result:Result<(),String> = async {
                    let response = result?;
                    let AnalysisResult::Contigs(fresh) = response.result else {
                        return Err(
                            "Fresh assembly replay failed or returned another result.".into()
                        );
                    };
                    fresh
                        .validate_bindings_for(&reads, &opts)
                        .map_err(|e| e.to_string())?;
                    let c = fresh
                        .contigs
                        .iter()
                        .find(|c| c.id == contig_id)
                        .ok_or("The selected layout changed during replay.")?;
                    let m = c
                        .molecule
                        .clone()
                        .ok_or("The fresh layout has no releasable consensus.")?;
                    if m.revision() != expected
                        || c.sequence != expected_sequence
                        || c.placements
                            .iter()
                            .map(|p| p.read_index)
                            .collect::<Vec<_>>()
                            != old_read_indices
                    {
                        return Err("The fresh source reconstruction differs from the displayed consensus. Assemble again before release.".into());
                    }
                    let export_guard=||guard.run(bound.clone());
                    if kind == "fasta" {
                        crate::io::download::download_guarded(
                            "contig-consensus.fasta",
                            "text/plain",
                            collagene_core::formats::sequence::export_fasta(&m).as_bytes(),&export_guard,
                        ).await?;
                    } else {
                        let text = collagene_core::formats::sequence::export_genbank(&m)
                            .map_err(|e| e.to_string())?;
                        crate::io::download::download_guarded("contig-consensus.gb", "text/plain", text.as_bytes(),&export_guard).await?;
                    }
                    Ok(())
                }.await;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || guard.run(bound.clone()).is_err()
                {
                    return;
                }
                match result{Ok(())=>status.set("Freshly reconstructed consensus released with its source bindings. Saved evidence never restores current verification.".into()),Err(e)=>status.set(e)}
            });
        })
    };
    let keep = Callback::new(move |()| {
        let Some(bound) = data.get_untracked() else {
            return;
        };
        if let Err(error) = guard.run(bound.clone()) {
            status.set(error);
            return;
        }
        if let Err(error) = on_record.run(ArchiveRequest {
            report: bound.report.clone(),
            bindings: bound.bindings.clone(),
            contig_id: None,
            allow_uncertain: true,
            epoch: bound.epoch,
        }) {
            status.set(error);
        }
    });
    let can_release = Memo::new(move |_| {
        data.with(|b| {
            b.as_ref()
                .and_then(|b| {
                    current_contig(b, contig_index.get())
                        .map(|c| c.molecule.is_some() && (!needs_ack(b, c) || acknowledge.get()))
                })
                .unwrap_or(false)
        })
    });
    let export_report = Callback::new(move |()| {
        let Some(bound) = data.get_untracked() else {
            return;
        };
        if let Err(e) = guard.run(bound.clone()) {
            status.set(e);
            return;
        }
        let mut bytes = Json(vec![]);
        match serde_json::to_writer_pretty(&mut bytes, &*bound.report) {
            Ok(()) => crate::io::download::download_with_completion(
                "contig-assembly-report.json",
                "application/json",
                &bytes.0,
                move |result| {
                    if status.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        status.set(e)
                    }
                },
            ),
            Err(e) => status.set(e.to_string()),
        }
    });
    view! {<section class="contig-panel" class:contig-expanded=move||expanded.get()>
        <div class="contig-heading"><h2>"Assemble reads"</h2><Show when=move||expanded.get()><button class="button secondary" on:click=move|_|expanded.set(false)>"Close expanded contigs"</button></Show></div>
        <p class="helper">"Build contigs from the selected read records without an external reference. Branches, tied overlaps, gaps and uncertain calls remain visible. Native trace scores are not automatically calibrated error probabilities."</p>
        <details class="contig-read-picker" open><summary>"Read documents"</summary><p>{move||format!("{} reads selected in explicit order",selected.with(Vec::len))}</p><button class="button secondary" on:click=move|_|selected.set(read_options.with(|values|values.iter().filter(|v|v.eligible).take(32).map(|v|v.document_id.clone()).collect()))>"Select linear DNA reads"</button><button class="button secondary" on:click=move|_|selected.set(vec![])>"Clear selected reads"</button>
            {move||read_options.with(|values|values.iter().skip(source_page.get()*25).take(25).map(|info|{let id=info.document_id.clone();let check_id=id.clone();let select_id=id.clone();view!{<label class="contig-read-choice"><input type="checkbox" aria-label=format!("Include contig read {}",info.document_id) disabled=!info.eligible prop:checked=move||selected.with(|s|s.contains(&check_id)) on:change=move|e|{let checked=event_target_checked(&e);selected.update(|values|{if checked{if !values.contains(&select_id){values.push(select_id.clone());}}else{values.retain(|id|id!=&select_id);}});}/><span><strong>{info.name.clone()}</strong><small>{format!("{} bases · {} · {}",info.length,info.quality_label,info.document_id)}</small></span></label>}}).collect_view())}
            <button class="button secondary" disabled=move||source_page.get()==0 on:click=move|_|source_page.update(|p|*p=p.saturating_sub(1))>"Previous source reads"</button><button class="button secondary" disabled=move||{(source_page.get()+1)*25>=read_options.with(Vec::len)} on:click=move|_|source_page.update(|p|*p+=1)>"Next source reads"</button>
            <ol class="contig-selected-order">{move||selected.get().iter().enumerate().map(|(i,id)|{let up=id.clone();let down=id.clone();let remove=id.clone();let label=read_options.with(|rows|rows.iter().find(|r|&r.document_id==id).map(|r|r.name.clone()).unwrap_or_else(||"Missing read".into()));view!{<li>{format!("read_{} · {label}",i+1)}<button class="text-tool" disabled=i==0 on:click=move|_|selected.update(|values|{if let Some(i)=values.iter().position(|v|v==&up){if i>0{values.swap(i,i-1);}}})>"Move read up"</button><button class="text-tool" on:click=move|_|selected.update(|values|{if let Some(i)=values.iter().position(|v|v==&down){if i+1<values.len(){values.swap(i,i+1);}}})>"Move read down"</button><button class="text-tool" on:click=move|_|selected.update(|values|values.retain(|v|v!=&remove))>"Remove selected read"</button></li>}}).collect_view()}</ol>
        </details>
        <details><summary>"Trimming and overlap policy"</summary>
            <label class="checkbox-label"><input type="checkbox" aria-label="Contig quality-window trimming" prop:checked=move||trim_windows.get() on:change=move|e|trim_windows.set(event_target_checked(&e))/>"Trim by native-score windows"</label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Contig allow unscored reads" prop:checked=move||unscored.get() on:change=move|e|unscored.set(event_target_checked(&e))/>"Explicitly allow unscored read calls"</label>
            <label>"Raw confidence threshold"<input aria-label="Contig raw score threshold" prop:value=move||threshold.get() on:input=move|e|threshold.set(event_target_value(&e))/></label>
            <label>"Trim window"<input aria-label="Contig trim window" prop:value=move||trim_window.get() on:input=move|e|trim_window.set(event_target_value(&e))/></label>
            <label>"Passing calls per window"<input aria-label="Contig passing calls" prop:value=move||passing.get() on:input=move|e|passing.set(event_target_value(&e))/></label>
            <label>"Minimum retained bases"<input aria-label="Contig minimum retained bases" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Minimum overlap"<input aria-label="Contig minimum overlap" prop:value=move||min_overlap.get() on:input=move|e|min_overlap.set(event_target_value(&e))/></label>
            <label>"Maximum non-certain overlap columns (%)"<input aria-label="Contig maximum error percent" prop:value=move||error.get() on:input=move|e|error.set(event_target_value(&e))/></label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Contig both strands" prop:checked=move||both.get() on:change=move|e|both.set(event_target_checked(&e))/>"Search both read orientations"</label>
            <label>"Minimum eligible call support"<input aria-label="Contig minimum support" prop:value=move||support.get() on:input=move|e|support.set(event_target_value(&e))/></label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Contig require both strands" prop:checked=move||both_support.get() on:change=move|e|both_support.set(event_target_checked(&e))/>"Require calls from both original strands"</label>
        </details>
        <details><summary>"Scoring and resource bounds"</summary>
            <label>"Match score"<input aria-label="Contig match score" prop:value=move||match_score.get() on:input=move|e|match_score.set(event_target_value(&e))/></label><label>"Ambiguity score"<input aria-label="Contig ambiguity score" prop:value=move||ambiguous_score.get() on:input=move|e|ambiguous_score.set(event_target_value(&e))/></label><label>"Mismatch score"<input aria-label="Contig mismatch score" prop:value=move||mismatch_score.get() on:input=move|e|mismatch_score.set(event_target_value(&e))/></label><label>"Gap opening score"<input aria-label="Contig gap open" prop:value=move||gap_open.get() on:input=move|e|gap_open.set(event_target_value(&e))/></label><label>"Gap extension score"<input aria-label="Contig gap extend" prop:value=move||gap_extend.get() on:input=move|e|gap_extend.set(event_target_value(&e))/></label>
            <label>"Per-pair cells"<input aria-label="Contig pair cell budget" prop:value=move||pair_cells.get() on:input=move|e|pair_cells.set(event_target_value(&e))/></label><label>"Total cells"<input aria-label="Contig total cell budget" prop:value=move||total_cells.get() on:input=move|e|total_cells.set(event_target_value(&e))/></label><label>"Overlap edge limit"<input aria-label="Contig edge limit" prop:value=move||max_edges.get() on:input=move|e|max_edges.set(event_target_value(&e))/></label><label>"Layout column limit"<input aria-label="Contig column limit" prop:value=move||max_columns.get() on:input=move|e|max_columns.set(event_target_value(&e))/></label>
        </details>
        <div class="contig-actions"><button class="button primary" disabled=move||busy.get() on:click=run>"Assemble contigs"</button><button class="button secondary" disabled=move||!busy.get() on:click=move|_|cancel.run(())>"Cancel contig assembly"</button><button class="button secondary" on:click=move|_|expanded.update(|v|*v=!*v)>"Expand contigs"</button></div>
        <p class="contig-status" role="status">{move||status.get()}</p>
        <Show when=move||data.with(|d|d.is_some())><button class="button secondary" disabled=move||busy.get() on:click=move|_|keep.run(())>"Keep assembly evidence"</button></Show>
        <ContigEvidence data=data busy=busy state=EvidenceState{contig_index,column_start,column_page,row_start,selected_column,edge_page,selected_edge,acknowledge} expanded=expanded on_navigate=on_navigate release=release export_report=export_report can_release=can_release/>

    </section>}
}

#[derive(Clone, Copy)]
struct EvidenceState {
    contig_index: RwSignal<usize>,
    column_start: RwSignal<usize>,
    column_page: RwSignal<usize>,
    row_start: RwSignal<usize>,
    selected_column: RwSignal<Option<usize>>,
    edge_page: RwSignal<usize>,
    selected_edge: RwSignal<Option<usize>>,
    acknowledge: RwSignal<bool>,
}
fn evidence_state() -> EvidenceState {
    EvidenceState {
        contig_index: RwSignal::new(0),
        column_start: RwSignal::new(0),
        column_page: RwSignal::new(0),
        row_start: RwSignal::new(0),
        selected_column: RwSignal::new(None),
        edge_page: RwSignal::new(0),
        selected_edge: RwSignal::new(None),
        acknowledge: RwSignal::new(false),
    }
}
#[component]
fn ContigEvidence(
    data: RwSignal<Option<Arc<BoundAssembly>>>,
    busy: RwSignal<bool>,
    state: EvidenceState,
    expanded: RwSignal<bool>,
    on_navigate: Callback<ReadNavigation>,
    release: Callback<&'static str>,
    export_report: Callback<()>,
    can_release: Memo<bool>,
) -> impl IntoView {
    let EvidenceState {
        contig_index,
        column_start,
        column_page,
        row_start,
        selected_column,
        edge_page,
        selected_edge,
        acknowledge,
    } = state;
    view! {        {move||data.get().map(|bound|view!{<div class="contig-result" data-complete=bound.report.complete.to_string()><p>{format!("{} contigs · {} overlap candidates · {} cells · search complete: {}",bound.report.contigs.len(),bound.report.overlaps.len(),bound.report.cells_used,bound.report.complete)}</p><button class="button secondary" on:click=move|_|export_report.run(())>"Download contig report"</button>
        <details><summary>"Read trimming and evidence"</summary>{bound.report.reads.iter().enumerate().map(|(i,r)|view!{<article class="contig-read-report"><strong>{format!("{} · {}",r.id,r.name)}</strong><p>{format!("Current calls {} · retained {}–{} · {} eligible · {}",r.source_length,r.trim_start+1,r.trim_end,r.eligible_calls,r.quality_source)}</p><p>{format!("Original calls {} · source retained {}–{} · source reversed {}",r.source_call_length,r.retained_start+1,r.retained_end,r.source_reverse)}</p>{r.duplicate_group.clone().map(|group|view!{<p class="helper">{format!("Duplicate evidence group: {group}; independent physical observations are not inferred.")}</p>})}{r.warnings.iter().cloned().map(|w|view!{<p class="helper">{w}</p>}).collect_view()}{r.error.clone().map(|e|view!{<p role="alert">{e}</p>})}<button class="button secondary" disabled={r.trim_start>=r.trim_end} on:click={let binding=bound.bindings.get(i).cloned();let at=r.trim_start;let source_epoch=bound.epoch;move|_|if let Some(binding)=binding.clone(){on_navigate.run(ReadNavigation{binding,call_index:at,epoch:source_epoch});}}>"Open retained read call"</button></article>}).collect_view()}</details>
        <details class="contig-overlaps" open><summary>"Overlap graph and alternatives"</summary><p class="helper">"Every retained directed overlap is available in pages. Selected edges describe this model's layout; other edges and tied paths are not erased."</p>
            <svg class="contig-overlap-graph" viewBox="0 0 400 260" role="img" aria-label="Current page of directed read overlaps">{bound.report.overlaps.iter().skip(edge_page.get()*25).take(25).enumerate().map(|(i,e)|{let y=15+i as i32*9;view!{<g><text x="2" y=y>{format!("{}{}",e.from+1,if e.from_reverse{"−"}else{"+"})}</text><line x1="35" y1=y-3 x2="350" y2=y-3 stroke=if e.selected{"#287c5e"}else if e.accepted{"#829d8d"}else{"#bb9a74"}/><text x="355" y=y>{format!("{}{}",e.to+1,if e.to_reverse{"−"}else{"+"})}</text></g>}}).collect_view()}</svg>
            {bound.report.overlaps.iter().skip(edge_page.get()*25).take(25).map(|e|{let id=e.id;view!{<button class="button secondary contig-overlap-entry" data-overlap-id=e.id on:click=move|_|selected_edge.set(Some(id))>{format!("{}{} → {}{} · score {} · {} paired · {} errors · paths {} · {}",e.from+1,if e.from_reverse{"−"}else{"+"},e.to+1,if e.to_reverse{"−"}else{"+"},e.score,e.paired_bases,e.error_columns,e.optimal_paths_saturated,if e.selected{"selected"}else if e.accepted{"alternative"}else{"excluded"})}</button>}}).collect_view()}
            <button class="button secondary" disabled=move||edge_page.get()==0 on:click=move|_|edge_page.update(|p|*p=p.saturating_sub(1))>"Previous overlaps"</button><button class="button secondary" disabled=move||{(edge_page.get()+1)*25>=data.with(|b|b.as_ref().map_or(0,|b|b.report.overlaps.len()))} on:click=move|_|edge_page.update(|p|*p+=1)>"Next overlaps"</button>
            {move||data.get().and_then(|b|selected_edge.get().and_then(|id|b.report.overlaps.iter().find(|e|e.id==id).cloned())).map(|e|view!{<div class="contig-overlap-detail"><p>{format!("{:?}; accepted {}; selected {}; oriented-trimmed positions {}–{} / {}–{}",e.kind,e.accepted,e.selected,e.from_start,e.from_end,e.to_start,e.to_end)}</p>{e.reason.map(|r|view!{<p>{r}</p>})}<pre>{e.reference_aligned}</pre><pre>{e.query_aligned}</pre></div>})}
        </details>
        <h3>"Contig layouts"</h3><select aria-label="Contig layout" prop:value=move||contig_index.get().to_string() on:change=move|e|{contig_index.set(event_target_value(&e).parse().unwrap_or(0));column_start.set(0);column_page.set(0);selected_column.set(None);row_start.set(0);acknowledge.set(false);}>{bound.report.contigs.iter().enumerate().map(|(i,c)|view!{<option value=i.to_string()>{format!("{} · {} bases · {} reads",c.name,c.sequence.len(),c.placements.len())}</option>}).collect_view()}</select>
        {move||data.get().and_then(|b|current_contig(&b,contig_index.get()).map(|c|(c.sequence.len(),c.columns.len(),c.placements.len(),c.complete,c.layout_unique,c.warnings.clone()))).map(|(len,columns,reads,complete,unique,warnings)|view!{<p>{format!("{len} consensus bases · {columns} layout columns · {reads} read placements · complete {complete} · layout unique under this model {unique}")}</p>{warnings.into_iter().map(|w|view!{<p class="helper">{w}</p>}).collect_view()}})}
        <details><summary>"Read placement coordinates"</summary>{move||data.get().map(|bound|current_contig(&bound,contig_index.get()).map(|c|c.placements.iter().map(|p|{let source=bound.report.reads.get(p.read_index);let binding=bound.bindings.get(p.read_index).cloned();let call=p.calls.first().map(|c|c.call_index);let source_epoch=bound.epoch;view!{<article class="contig-placement"><p>{format!("{} · layout columns {}–{} · {} mapped calls · current read reversed {} · original attachment reversed {}",source.map(|r|r.name.as_str()).unwrap_or("read"),p.start+1,p.end,p.calls.len(),p.reverse,source.is_some_and(|r|r.source_reverse)^p.reverse)}</p><button class="button secondary" disabled=call.is_none() on:click=move|_|if let(Some(binding),Some(call_index))=(binding.clone(),call){on_navigate.run(ReadNavigation{binding,call_index,epoch:source_epoch});}>"Open first placed call"</button></article>}}).collect_view()).into_view())}</details>
        <label>"First layout column (1-based)"<input aria-label="Contig first column" prop:value=move||(column_start.get()+1).to_string() on:change=move|e|{if let Ok(value)=event_target_value(&e).parse::<usize>(){let max=data.with_untracked(|b|b.as_ref().and_then(|b|current_contig(b,contig_index.get_untracked())).map_or(1,|c|c.columns.len()).max(1));column_start.set(value.saturating_sub(1).min(max-1));}}/></label>
        <LayoutCanvas data=data contig_index=contig_index start=column_start selected_column=selected_column row_start=row_start expanded=expanded/>
        <div class="contig-actions"><button class="button secondary" on:click=move|_|column_start.update(|p|*p=p.saturating_sub(50))>"Previous contig columns"</button><button class="button secondary" on:click=move|_|{let max=data.with_untracked(|b|b.as_ref().and_then(|b|current_contig(b,contig_index.get_untracked())).map_or(1,|c|c.columns.len()).max(1));column_start.update(|p|*p=(*p+50).min(max-1));}>"Next contig columns"</button><button class="button secondary" on:click=move|_|row_start.update(|p|*p=p.saturating_sub(12))>"Previous layout reads"</button><button class="button secondary" disabled=move||{row_start.get()+12>=data.with(|b|b.as_ref().and_then(|b|current_contig(b,contig_index.get())).map_or(0,|c|c.placements.len()))} on:click=move|_|row_start.update(|p|*p+=12)>"Next layout reads"</button></div>
        <p class="helper">"The canvas paints at most 160 columns and 12 read rows. Exact column evidence remains available below; no consensus base is invented from an external reference."</p>
        {move||data.get().and_then(|b|selected_column.get().and_then(|at|current_contig(&b,contig_index.get()).and_then(|c|c.columns.get(at)).map(|c|(b.clone(),c.clone())))).map(|(b,c)|view!{<article class="contig-column-detail" data-layout-column=c.layout_index><h4>{format!("Layout column {} · {}",c.layout_index+1,state_label(&c.state))}</h4><p>{format!("Consensus position {:?} · A/C/G/T counts {:?} · gaps {:?} · eligible gaps {:?}",c.consensus_position.map(|v|v+1),c.counts,c.gap_reads,c.eligible_gap_reads)}</p>{c.observations.into_iter().map(|o|{let source=b.report.reads.get(o.read_index);let label=source.map(|r|r.name.clone()).unwrap_or_else(||"read".into());let binding=b.bindings.get(o.read_index).cloned();let source_epoch=b.epoch;view!{<div class="contig-call"><p>{format!("{}: {} · current call {} · original call {} · raw score {:?} · eligible {} · current/original trace sample indices {:?}/{:?}",label,o.base,o.call_index+1,o.stored_call_index+1,o.raw_score,o.eligible,o.sample_index,o.stored_sample_index)}</p><button class="button secondary" on:click=move|_|if let Some(binding)=binding.clone(){on_navigate.run(ReadNavigation{binding,call_index:o.call_index,epoch:source_epoch});}>"Open source call"</button></div>}}).collect_view()}</article>})}
        <details><summary>"Exact layout columns"</summary>{move||data.get().map(|b|current_contig(&b,contig_index.get()).map(|c|c.columns.iter().skip(column_page.get()*25).take(25).map(|c|{let at=c.layout_index;view!{<button class="button secondary contig-column-entry" data-layout-column=c.layout_index on:click=move|_|{selected_column.set(Some(at));column_start.set(at.saturating_sub(5));}>{format!("Column {} · {} · {} · {} calls",c.layout_index+1,c.base.unwrap_or('-'),state_label(&c.state),c.observations.len())}</button>}}).collect_view()).into_view())}<button class="button secondary" disabled=move||column_page.get()==0 on:click=move|_|column_page.update(|p|*p=p.saturating_sub(1))>"Previous evidence columns"</button><button class="button secondary" disabled=move||{(column_page.get()+1)*25>=data.with(|b|b.as_ref().and_then(|b|current_contig(b,contig_index.get())).map_or(0,|c|c.columns.len()))} on:click=move|_|column_page.update(|p|*p+=1)>"Next evidence columns"</button></details>
        <label class="checkbox-label"><input type="checkbox" aria-label="Acknowledge uncertain contig consensus" prop:checked=move||acknowledge.get() on:change=move|e|acknowledge.set(event_target_checked(&e))/>"I understand unresolved calls, incomplete searches or ambiguous layouts do not establish a unique molecule."</label>
        <p class="helper">"FASTA omits quality, call and source-event evidence. GenBank comments retain producer summaries, not the full per-call/gap/history model. Keep the native workspace and JSON report for reconstruction."</p>
        <div class="contig-actions"><button class="button primary" disabled=move||busy.get()||!can_release.get() on:click=move|_|release.run("open")>"Open reconstructed contig"</button><button class="button secondary" disabled=move||busy.get()||!can_release.get() on:click=move|_|release.run("fasta")>"Download contig FASTA"</button><button class="button secondary" disabled=move||busy.get()||!can_release.get() on:click=move|_|release.run("genbank")>"Download contig GenBank"</button></div>
        <h4>"Unassembled reads"</h4>{bound.report.unassembled.iter().map(|r|view!{<p class="helper">{format!("read {}: {}",r.read_index+1,r.reason)}</p>}).collect_view()}{bound.report.limits.iter().chain(&bound.report.warnings).map(|w|view!{<p class="helper">{w.clone()}</p>}).collect_view()}
    </div>})}}
}

#[component]
pub fn ContigArchivePanel(
    request: ArchiveRequest,
    snapshot: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    validate: Callback<collagene_core::documents::workspace::Workspace, Result<(), String>>,
    on_restore: Callback<collagene_core::documents::workspace::Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let request = StoredValue::new(request);
    let status = RwSignal::new(String::new());
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, Ordering::Relaxed));
    let record = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(database) = db.get_value() else {
                status.set(
                    "Browser storage is unavailable; the assembly has not been recorded.".into(),
                );
                return;
            };
            let input = request.get_value();
            if epoch.get_untracked() != input.epoch {
                status.set(
                    "The workspace changed. Reopen the source-bound assembly before recording."
                        .into(),
                );
                return;
            }
            busy.set(true);
            save_generation.update(|n| *n = n.wrapping_add(1));
            let prepared = (|| -> Result<_, String> {
                let frame = snapshot.run(())?;
                let fingerprint = collagene_core::material::lineage::workspace_fingerprint(&frame)
                    .map_err(|e| e.to_string())?;
                if input.bindings.len() != input.report.reads.len() {
                    return Err(
                        "Read-event binding count does not match the assembly report.".into(),
                    );
                }
                let bindings = input
                    .bindings
                    .iter()
                    .zip(&input.report.reads)
                    .map(
                        |(b, r)| collagene_core::material::lineage::InputEventBinding {
                            input_id: r.id.clone(),
                            event: collagene_core::material::lineage::EventBinding {
                                document_id: b.document_id.clone(),
                                event_id: b.event_id.clone(),
                                revision: b.revision.clone(),
                            },
                        },
                    )
                    .collect::<Vec<_>>();
                let saved = collagene_core::review::contigs::history::SavedAssembly::new(
                    input.report.as_ref().clone(),
                    bindings.clone(),
                )
                .map_err(|e| e.to_string())?;
                let reads =
                    collagene_core::review::contigs::history::resolve_inputs(&frame, &saved)
                        .map_err(|e| e.to_string())?;
                let revision = contig::batch_revision(&reads, &input.report.options)
                    .map_err(|e| e.to_string())?;
                Ok((frame, fingerprint, bindings, reads, revision))
            })();
            let (frame, fingerprint, bindings, reads, revision) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    busy.set(false);
                    return;
                }
            };
            let alive = alive.clone();
            status.set(
                "Reconstructing the exact read events before an atomic workspace update…".into(),
            );
            spawn_local(async move {
                let result=async{
                let response=crate::platform::worker::analyze_cancellable(Request{id:u32::MAX-30,revision,job:AnalysisJob::AssembleContigs{reads:reads.clone(),options:input.report.options.clone()}},alive.clone()).await?;
                if response.id!=u32::MAX-30{return Err("Assembly response identity mismatch.".into());}
                let fresh=match response.result{AnalysisResult::Contigs(v)=>v,AnalysisResult::Error(e)=>return Err(e),_=>return Err("Unexpected assembly replay response.".into())};
                fresh.validate_bindings_for(&reads,&input.report.options).map_err(|e|e.to_string())?;
                let guard=||->Result<(),String>{if !alive.load(Ordering::Relaxed)||epoch.get_untracked()!=input.epoch{return Err("The dialog closed or the workspace changed. Nothing was recorded.".into());}let current=snapshot.run(())?;if collagene_core::material::lineage::workspace_fingerprint(&current).map_err(|e|e.to_string())?!=fingerprint{return Err("The workspace changed during reconstruction. Nothing was recorded.".into());}Ok(())};
                guard()?;
                let saved=collagene_core::review::contigs::history::SavedAssembly::new(fresh.clone(),bindings.clone()).map_err(|e|e.to_string())?;
                let(hash,attachment)=collagene_core::review::contigs::history::make_attachment(&saved).map_err(|e|e.to_string())?;
                let mut next=frame.clone();next.attachments.insert(hash.clone(),attachment);
                if let Some(id)=input.contig_id{
                    let old=input.report.contigs.iter().find(|c|c.id==id).ok_or("The selected historical contig is unavailable.")?;
                    let current=fresh.contigs.iter().find(|c|c.id==id).ok_or("The selected layout changed on replay. Inspect the fresh assembly before opening it.")?;
                    let mut product=current.molecule.clone().ok_or("The fresh layout has no releasable consensus molecule.")?;
                    if old.sequence!=current.sequence||old.molecule.as_ref().map(Molecule::revision)!=Some(product.revision())||old.placements.iter().map(|p|p.read_index).collect::<Vec<_>>()!=current.placements.iter().map(|p|p.read_index).collect::<Vec<_>>(){return Err("The fresh reconstruction differs from the displayed layout. Reassemble before opening it.".into());}
                    if (!fresh.complete||!current.complete||!current.layout_unique||current.sequence.bytes().any(|b|!b"ACGT".contains(&b)))&&!input.allow_uncertain{return Err("Acknowledge incomplete or uncertain consensus evidence before opening this layout.".into());}
                    let scientific_product_revision=product.revision();
                    product.metadata.insert("contig:scientific_product_revision".into(),scientific_product_revision);
                    product.metadata.insert("contig:evidence_hash".into(),hash.clone());
                    let histories=collagene_core::documents::history::DocumentHistory::restore_all(&next).map_err(|e|e.to_string())?;
                    let source_bindings=current.placements.iter().map(|p|bindings.get(p.read_index).ok_or("A contributing source binding is missing.")).collect::<Result<Vec<_>,_>>()?;
                    let sources=source_bindings.iter().map(|b|histories.iter().find(|h|h.id==b.event.document_id).map(|h|(h,b.event.event_id.as_str())).ok_or("A contributing read history is unavailable.")).collect::<Result<Vec<_>,_>>()?;
                    let base=format!("assembly-{}-{}",&hash[..24],id);let mut document_id=base.clone();let mut suffix=1usize;while next.documents.iter().any(|d|d.id==document_id){document_id=format!("{base}-{suffix}");suffix+=1;}
                    let derived=collagene_core::documents::history::DocumentHistory::derive_from_events(&sources,document_id.clone(),product,format!("Read contig assembly {hash} contig {id}; computational evidence, not experimental provenance")).map_err(|e|e.to_string())?;
                    derived.contribute(&mut next).map_err(|e|e.to_string())?;next.active_document=document_id;
                }
                next.validate().map_err(|e|e.to_string())?;validate.run(next.clone())?;guard()?;
                let current=crate::workspace::session::encode_workspace(next.clone(),u32::MAX-31).await?;let previous=crate::workspace::session::encode_workspace(frame,u32::MAX-32).await?;guard()?;
                crate::platform::storage::replace_workspace(&database,&current,&previous).await?;
                Ok::<_,String>(next)
            }.await;
                match result {
                    Ok(frame) => {
                        loaded.set(true);
                        message.set("Assembly evidence recorded with exact source events. Historical records do not restore verification; consensus use requires fresh reconstruction.".into());
                        on_restore.run(frame);
                    }
                    Err(e) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("Assembly was not recorded. {e}"));
                        }
                    }
                }
                busy.set(false);
            });
        }
    };
    view! {<section class="contig-archive-panel"><p>"Recording replays the assembly against its exact retained read events, then updates the current and recovery workspaces in one transaction."</p><p class="helper">"A stored report is historical producer evidence, not a verified physical molecule. Read heads and undo/redo remain unchanged. Opening a consensus is an explicit action and gives its document every contributing source-event parent."</p><p>{move||request.with_value(|r|format!("{} source events · {} recorded contigs · {}",r.bindings.len(),r.report.contigs.len(),if r.contig_id.is_some(){"Record evidence and open the selected consensus"}else{"Record evidence without opening a consensus"}))}</p><button class="button primary" disabled=move||busy.get() on:click=record>"Record read assembly"</button><p role="status">{move||status.get()}</p></section>}
}

#[component]
pub fn SavedAssemblyPanel(
    attachment_hash: String,
    snapshot: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    epoch: RwSignal<u64>,
    on_record: Callback<ArchiveRequest, Result<(), String>>,
    on_navigate: Callback<ReadNavigation>,
) -> impl IntoView {
    let initial = (|| -> Result<_, String> {
        let frame = snapshot.run(())?;
        frame.validate().map_err(|e| e.to_string())?;
        let a = frame
            .attachments
            .get(&attachment_hash)
            .ok_or("The saved assembly attachment is unavailable.")?;
        if a.kind != collagene_core::review::contigs::history::WORKSPACE_ATTACHMENT_KIND {
            return Err("This attachment is not read-assembly evidence.".into());
        }
        let saved = collagene_core::review::contigs::history::load_record(&a.bytes)
            .map_err(|e| e.to_string())?;
        let _ = collagene_core::review::contigs::history::resolve_inputs(&frame, &saved)
            .map_err(|e| e.to_string())?;
        let bindings = saved
            .inputs
            .iter()
            .map(|b| ReadBinding {
                document_id: b.event.document_id.clone(),
                event_id: b.event.event_id.clone(),
                revision: b.event.revision.clone(),
            })
            .collect::<Vec<_>>();
        Ok((saved, bindings))
    })();
    let error = initial.as_ref().err().cloned().unwrap_or_default();
    let data = RwSignal::new(initial.as_ref().ok().map(|(saved, bindings)| {
        Arc::new(BoundAssembly {
            report: Arc::new(saved.report.clone()),
            bindings: bindings.clone(),
            epoch: epoch.get_untracked(),
        })
    }));
    let historical = StoredValue::new(initial.ok().map(|(record, _)| Arc::new(record)));
    let fresh = RwSignal::new(false);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(error);
    let state = evidence_state();
    let expanded = RwSignal::new(false);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    Effect::new(move |_| {
        let _ = epoch.get();
        fresh.set(false);
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
    });
    let replay = {
        let alive = alive.clone();
        move |_| {
            let Some(saved) = historical.get_value() else {
                return;
            };
            let prepared = (|| -> Result<_, String> {
                let frame = snapshot.run(())?;
                let reads =
                    collagene_core::review::contigs::history::resolve_inputs(&frame, &saved)
                        .map_err(|e| e.to_string())?;
                let revision = contig::batch_revision(&reads, &saved.report.options)
                    .map_err(|e| e.to_string())?;
                Ok((reads, revision))
            })();
            let (reads, revision) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            generation.update(|n| *n = n.wrapping_add(1));
            let token = generation.get_untracked();
            let expected_epoch = epoch.get_untracked();
            busy.set(true);
            fresh.set(false);
            status
                .set("Replaying the saved exact source events with the current assembler…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision,
                        job: AnalysisJob::AssembleContigs {
                            reads: reads.clone(),
                            options: saved.report.options.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || epoch.get_untracked() != expected_epoch
                {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                let result = (|| -> Result<_, String> {
                    let response = response?;
                    let report = match response.result {
                        AnalysisResult::Contigs(v) => v,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected saved-assembly replay response.".into()),
                    };
                    report
                        .validate_bindings_for(&reads, &saved.report.options)
                        .map_err(|e| e.to_string())?;
                    let frame = snapshot.run(())?;
                    let current =
                        collagene_core::review::contigs::history::resolve_inputs(&frame, &saved)
                            .map_err(|e| e.to_string())?;
                    if contig::batch_revision(&current, &saved.report.options)
                        .map_err(|e| e.to_string())?
                        != report.batch_revision
                    {
                        return Err("The retained source bindings changed during replay.".into());
                    }
                    let bindings = saved
                        .inputs
                        .iter()
                        .map(|b| ReadBinding {
                            document_id: b.event.document_id.clone(),
                            event_id: b.event.event_id.clone(),
                            revision: b.event.revision.clone(),
                        })
                        .collect();
                    Ok((report, bindings))
                })();
                match result {
                    Ok((report, bindings)) => {
                        let old = saved
                            .report
                            .contigs
                            .iter()
                            .map(|c| {
                                (
                                    &c.sequence,
                                    c.placements
                                        .iter()
                                        .map(|p| p.read_index)
                                        .collect::<Vec<_>>(),
                                )
                            })
                            .collect::<Vec<_>>();
                        let new = report
                            .contigs
                            .iter()
                            .map(|c| {
                                (
                                    &c.sequence,
                                    c.placements
                                        .iter()
                                        .map(|p| p.read_index)
                                        .collect::<Vec<_>>(),
                                )
                            })
                            .collect::<Vec<_>>();
                        status.set(format!("Fresh replay completed. Engine matches recorded: {}. Consensus/layout membership matches recorded: {}. This is session evidence, not experimental identity.",report.engine_version==saved.report.engine_version,old==new));
                        data.set(Some(Arc::new(BoundAssembly {
                            report: Arc::new(report),
                            bindings,
                            epoch: expected_epoch,
                        })));
                        fresh.set(true);
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let release = {
        let alive = alive.clone();
        Callback::new(move |kind: &'static str| {
            if !alive.load(Ordering::Relaxed) || !fresh.get_untracked() {
                status.set("Replay this historical assembly before using a consensus.".into());
                return;
            }
            let Some(bound) = data.get_untracked() else {
                return;
            };
            if bound.epoch != epoch.get_untracked() {
                fresh.set(false);
                status.set("The workspace changed. Replay the retained evidence again.".into());
                return;
            }
            let Some(c) = current_contig(&bound, state.contig_index.get_untracked()) else {
                return;
            };
            if needs_ack(&bound, c) && !state.acknowledge.get_untracked() {
                status.set(
                    "Acknowledge incomplete or uncertain consensus evidence before use.".into(),
                );
                return;
            }
            let Some(molecule) = c.molecule.as_ref() else {
                status.set("This layout has no releasable consensus document.".into());
                return;
            };
            let saved = match collagene_core::review::contigs::history::SavedAssembly::new(
                bound.report.as_ref().clone(),
                bound
                    .bindings
                    .iter()
                    .zip(&bound.report.reads)
                    .map(
                        |(b, r)| collagene_core::material::lineage::InputEventBinding {
                            input_id: r.id.clone(),
                            event: collagene_core::material::lineage::EventBinding {
                                document_id: b.document_id.clone(),
                                event_id: b.event_id.clone(),
                                revision: b.revision.clone(),
                            },
                        },
                    )
                    .collect(),
            ) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let frame = match snapshot.run(()) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Err(e) = collagene_core::review::contigs::history::resolve_inputs(&frame, &saved)
            {
                status.set(e.to_string());
                return;
            }
            if kind == "open" {
                if let Err(e) = on_record.run(ArchiveRequest {
                    report: bound.report.clone(),
                    bindings: bound.bindings.clone(),
                    contig_id: Some(c.id),
                    allow_uncertain: state.acknowledge.get_untracked(),
                    epoch: bound.epoch,
                }) {
                    status.set(e);
                }
            } else {
                let text = if kind == "fasta" {
                    Ok(collagene_core::formats::sequence::export_fasta(molecule))
                } else {
                    collagene_core::formats::sequence::export_genbank(molecule)
                };
                match text {
                    Ok(text) => crate::io::download::download_with_completion(
                        if kind == "fasta" {
                            "contig-consensus.fasta"
                        } else {
                            "contig-consensus.gb"
                        },
                        "text/plain",
                        text.as_bytes(),
                        move |result| {
                            if status.try_get_untracked().is_none() {
                                return;
                            }
                            if let Err(e) = result {
                                status.set(e)
                            }
                        },
                    ),
                    Err(e) => status.set(e.to_string()),
                }
            }
        })
    };
    let can_release = Memo::new(move |_| {
        fresh.get()
            && data.with(|b| {
                b.as_ref()
                    .and_then(|b| {
                        current_contig(b, state.contig_index.get()).map(|c| {
                            c.molecule.is_some() && (!needs_ack(b, c) || state.acknowledge.get())
                        })
                    })
                    .unwrap_or(false)
            })
    });
    let export_report = Callback::new(move |()| {
        let Some(bound) = data.get_untracked() else {
            return;
        };
        let mut bytes = Json(vec![]);
        match serde_json::to_writer_pretty(&mut bytes, &*bound.report) {
            Ok(()) => crate::io::download::download_with_completion(
                "contig-assembly-report.json",
                "application/json",
                &bytes.0,
                move |result| {
                    if status.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        status.set(e)
                    }
                },
            ),
            Err(e) => status.set(e.to_string()),
        }
    });
    view! {<section class="saved-contig-panel" class:contig-expanded=move||expanded.get()><h2>"Saved read assembly"</h2><p class="saved-contig-authority">{move||if fresh.get(){"Fresh replay in this session"}else{"Historical producer evidence · unverified on opening"}}</p><p class="helper">"Saved statuses do not restore current verification. Replay resolves the exact retained input events even when current source heads have advanced. Opening an old source call may fork its retained revision without changing that source's current head."</p><button class="button primary" disabled=move||busy.get() on:click=replay>"Replay saved read assembly"</button><button class="button secondary" on:click=move|_|expanded.update(|v|*v=!*v)>{move||if expanded.get(){"Close expanded contigs"}else{"Expand contigs"}}</button><p role="status">{move||status.get()}</p><ContigEvidence data=data busy=busy state=state expanded=expanded on_navigate=on_navigate release=release export_report=export_report can_release=can_release/></section>}
}
