//! Virtualized Canvas2D chromatogram. Derived views always come from validated original evidence.
use crate::renderer::ResizeListener;
use collagene_core::traces::{self, QualityEvidence, Trace, TraceChannels};
use collagene_model::Molecule;
use leptos::{html, prelude::*};
use std::sync::Arc;
use wasm_bindgen::{JsCast, closure::Closure};
use web_sys::{CanvasRenderingContext2d as Ctx, HtmlCanvasElement};

const PAD: f64 = 24.;
const BASE_BIN: usize = 64;
const COLORS: [&str; 4] = ["#2f8f61", "#367cbc", "#b88820", "#c45569"];
const STYLE: &str = ".trace-view{display:flex;flex-direction:column;flex:1;min-height:0;background:white;overflow:hidden}.trace-controls{display:flex;align-items:center;flex-wrap:wrap;gap:8px;padding:10px 14px;border-bottom:1px solid #e4e9e1;font-size:12px}.trace-controls button{border:1px solid #d6dfd8;background:white;border-radius:4px;min-width:32px;padding:5px 8px;color:#365147}.trace-export{position:relative}.trace-export>summary{cursor:pointer;border:1px solid #d6dfd8;border-radius:4px;padding:5px 8px;color:#365147}.trace-export-content{position:absolute;z-index:10;right:0;top:100%;width:min(510px,80vw);max-height:65vh;overflow:auto;background:white;border:1px solid #d6dfd8;border-radius:5px;box-shadow:0 8px 28px #163b2226;padding:12px;line-height:1.5;overflow-wrap:anywhere}.trace-export-content p{margin:0 0 10px}.trace-export-content select{max-width:100%;padding:4px}.trace-export-content ul{padding-left:18px}.trace-export-actions{display:flex;gap:6px;margin:10px 0}.trace-export-review{border-top:1px solid #e4e9e1;margin-top:10px;padding-top:10px}.trace-export-review>label{display:block;margin-bottom:10px}.trace-controls button:hover{background:#edf4ef}.trace-legend{display:flex;gap:10px;margin-right:auto}.trace-legend span{font-weight:650}.trace-note,.trace-selection{font-size:12px;line-height:1.5;padding:8px 14px;margin:0;color:#65796e}.trace-scroller{position:relative;overflow:auto;flex:1;min-height:220px;outline-offset:-3px}.trace-spacer{position:relative;height:100%;min-height:220px}.trace-canvas{position:sticky;left:0;top:0;display:block;cursor:crosshair;outline-offset:-3px}.trace-error{padding:18px;color:#8c4c43;line-height:1.6;font-size:14px;overflow-wrap:anywhere}.trace-access{padding:8px 14px;font-size:12px;max-height:180px;overflow:auto;border-top:1px solid #e4e9e1}.trace-access summary{cursor:pointer;color:#365147}.trace-table-wrap{overflow:auto}.trace-access table{border-collapse:collapse;text-align:left;white-space:nowrap;font-size:11px}.trace-access th,.trace-access td{padding:5px 8px;border-bottom:1px solid #edf1ed}.trace-access caption{text-align:left;padding:8px 0;white-space:normal;color:#65796e}";

struct SignalPyramid {
    levels: Vec<Vec<(i32, i32)>>,
}
impl SignalPyramid {
    fn new(values: &[i32]) -> Self {
        let mut levels = vec![
            values
                .chunks(BASE_BIN)
                .map(|chunk| {
                    let low = *chunk.iter().min().unwrap();
                    let high = *chunk.iter().max().unwrap();
                    (low, high)
                })
                .collect::<Vec<_>>(),
        ];
        while levels.last().unwrap().len() > 1 {
            let next = levels
                .last()
                .unwrap()
                .chunks(2)
                .map(|pair| {
                    pair.iter()
                        .copied()
                        .fold((i32::MAX, i32::MIN), |(low, high), (a, b)| {
                            (low.min(a), high.max(b))
                        })
                })
                .collect();
            levels.push(next);
        }
        Self { levels }
    }
}
struct PlotData {
    trace: Arc<Trace>,
    pyramids: [SignalPyramid; 4],
    minimum: i32,
    maximum: i32,
}
impl PartialEq for PlotData {
    fn eq(&self, other: &Self) -> bool {
        self.trace.source_hash == other.trace.source_hash
            && self.trace.orientation == other.trace.orientation
            && self.trace.name == other.trace.name
    }
}
fn channels(c: &TraceChannels) -> [&[i32]; 4] {
    [&c.a, &c.c, &c.g, &c.t]
}
impl PlotData {
    fn new(trace: Trace) -> Self {
        let values = channels(&trace.channels);
        let minimum = values
            .iter()
            .flat_map(|v| v.iter())
            .copied()
            .min()
            .unwrap_or(0)
            .min(0);
        let maximum = values
            .iter()
            .flat_map(|v| v.iter())
            .copied()
            .max()
            .unwrap_or(1)
            .max(0);
        let pyramids = std::array::from_fn(|i| SignalPyramid::new(values[i]));
        Self {
            trace: Arc::new(trace),
            pyramids,
            minimum,
            maximum,
        }
    }
    fn max_zoom(&self) -> f64 {
        (16_000_000. / self.trace.channels.len() as f64)
            .min(128.)
            .max(0.125)
    }
}
fn evidence(trace: &Trace, index: usize) -> String {
    match &trace.qualities {
        Some(QualityEvidence::AbifBasecallerScores(values)) => values
            .get(index)
            .map(|q| format!("raw ABIF score {q}"))
            .unwrap_or_else(|| "unavailable".into()),
        Some(QualityEvidence::ScfAccuracy(q)) => format!(
            "SCF A:{}, C:{}, G:{}, T:{}",
            q.a[index], q.c[index], q.g[index], q.t[index]
        ),
        None => "unavailable".into(),
    }
}
fn quality_label(trace: &Trace) -> &'static str {
    match trace.qualities {
        Some(QualityEvidence::AbifBasecallerScores(_)) => {
            "Raw ABIF caller scores; scale not asserted Phred."
        }
        Some(QualityEvidence::ScfAccuracy(_)) => {
            "Four stored SCF accuracy estimates; no Phred conversion."
        }
        None => "No stored quality evidence; no scores inferred.",
    }
}
fn virtual_width(data: &PlotData, scale: f64, width: f64) -> f64 {
    (data.trace.channels.len() as f64 * scale + 2. * PAD).max(width)
}
fn visible_calls(trace: &Trace, scroll: f64, scale: f64, width: f64) -> (usize, usize) {
    let start = ((scroll - PAD) / scale).max(0.);
    let end = ((scroll + width - PAD) / scale).max(0.);
    let first = trace
        .peak_positions
        .partition_point(|&p| (p as f64) < start);
    let last = trace.peak_positions.partition_point(|&p| (p as f64) <= end);
    (first, last)
}
fn schedule_scroll(container: NodeRef<html::Div>, scroll: RwSignal<f64>, value: f64) {
    let Some(element) = container.get_untracked() else {
        return;
    };
    let callback = Closure::once_into_js(move || {
        element.set_scroll_left(value.max(0.).round() as i32);
        let _ = scroll.try_set(element.scroll_left() as f64);
    });
    if let Some(window) = web_sys::window() {
        let _ = window.request_animation_frame(callback.unchecked_ref());
    }
}
fn change_zoom(
    container: NodeRef<html::Div>,
    zoom: RwSignal<f64>,
    scroll: RwSignal<f64>,
    next: f64,
    focus: f64,
) {
    let width = container
        .get_untracked()
        .map(|e| e.client_width() as f64)
        .unwrap_or(500.);
    zoom.set(next);
    schedule_scroll(container, scroll, PAD + focus * next - width / 2.);
}
fn focus_sample(
    data: &PlotData,
    selection: Option<(usize, usize)>,
    scroll: f64,
    scale: f64,
    width: f64,
) -> f64 {
    selection
        .and_then(|(a, _)| data.trace.peak_positions.get(a))
        .map(|&p| p as f64)
        .unwrap_or(
            ((scroll + width / 2. - PAD) / scale)
                .clamp(0., data.trace.channels.len().saturating_sub(1) as f64),
        )
}
fn paint(
    canvas: &HtmlCanvasElement,
    data: &PlotData,
    width: f64,
    height: f64,
    scale: f64,
    scroll: f64,
    selection: Option<(usize, usize)>,
) -> Result<(), String> {
    let dpr = web_sys::window()
        .map(|w| w.device_pixel_ratio())
        .unwrap_or(1.)
        .clamp(1., 2.);
    let pixels_w = (width * dpr).round() as u32;
    let pixels_h = (height * dpr).round() as u32;
    if canvas.width() != pixels_w {
        canvas.set_width(pixels_w);
    }
    if canvas.height() != pixels_h {
        canvas.set_height(pixels_h);
    }
    canvas
        .unchecked_ref::<web_sys::HtmlElement>()
        .style()
        .set_property("width", &format!("{width}px"))
        .map_err(|_| "Cannot size chromatogram canvas")?;
    canvas
        .unchecked_ref::<web_sys::HtmlElement>()
        .style()
        .set_property("height", &format!("{height}px"))
        .map_err(|_| "Cannot size chromatogram canvas")?;
    let ctx = canvas
        .get_context("2d")
        .map_err(|_| "Cannot acquire chromatogram context")?
        .ok_or("Canvas2D is unavailable")?
        .dyn_into::<Ctx>()
        .map_err(|_| "Canvas2D is unavailable")?;
    ctx.set_transform(dpr, 0., 0., dpr, 0., 0.)
        .map_err(|_| "Cannot scale chromatogram")?;
    ctx.set_fill_style_str("#ffffff");
    ctx.fill_rect(0., 0., width, height);
    let trace = &data.trace;
    let n = trace.channels.len();
    let x = |sample: f64| PAD + sample * scale - scroll;
    let plot_top = 78.;
    let plot_bottom = (height - 35.).max(plot_top + 50.);
    let amplitude = (f64::from(data.maximum) - f64::from(data.minimum)).max(1.);
    let y = |value: i32| {
        plot_bottom
            - (f64::from(value) - f64::from(data.minimum)) / amplitude * (plot_bottom - plot_top)
    };
    if let Some((a, b)) = selection.filter(|&(a, b)| a < b && a < trace.sequence.len()) {
        let last = (b - 1).min(trace.sequence.len() - 1);
        let left = if a == 0 {
            0.
        } else {
            (f64::from(trace.peak_positions[a - 1]) + f64::from(trace.peak_positions[a])) / 2.
        };
        let right = if last + 1 >= trace.sequence.len() {
            n as f64
        } else {
            (f64::from(trace.peak_positions[last]) + f64::from(trace.peak_positions[last + 1])) / 2.
        };
        ctx.set_fill_style_str("#e3f2eb");
        let left = x(left).clamp(0., width);
        let right = x(right).clamp(0., width);
        ctx.fill_rect(left, 0., (right - left).max(2.), height);
    }
    ctx.set_stroke_style_str("#dfe7e1");
    ctx.set_line_width(1.);
    ctx.begin_path();
    ctx.move_to(0., y(0));
    ctx.line_to(width, y(0));
    ctx.stroke();
    let first = (((scroll - PAD) / scale).floor().max(0.) as usize)
        .saturating_sub(1)
        .min(n - 1);
    let last = (((scroll + width - PAD) / scale).ceil().max(0.) as usize)
        .saturating_add(2)
        .min(n);
    let raw = channels(&trace.channels);
    for channel in 0..4 {
        ctx.set_stroke_style_str(COLORS[channel]);
        ctx.set_line_width(1.3);
        ctx.begin_path();
        let samples_per_pixel = 1. / scale;
        if samples_per_pixel < BASE_BIN as f64 {
            for (position, &value) in raw[channel][first..last].iter().enumerate() {
                let px = x((first + position) as f64);
                if position == 0 {
                    ctx.move_to(px, y(value));
                } else {
                    ctx.line_to(px, y(value));
                }
            }
            if last - first == 1 {
                ctx.line_to(x(first as f64) + 1., y(raw[channel][first]));
            }
        } else {
            let mut level = 0;
            let mut bin = BASE_BIN;
            while level + 1 < data.pyramids[channel].levels.len()
                && (bin * 2) as f64 <= samples_per_pixel
            {
                level += 1;
                bin *= 2;
            }
            let values = &data.pyramids[channel].levels[level];
            let first_bin = (first / bin).min(values.len() - 1);
            let last_bin = last.div_ceil(bin).min(values.len());
            for (offset, &(low, high)) in values[first_bin..last_bin].iter().enumerate() {
                let center = ((first_bin + offset) * bin + bin / 2).min(n - 1);
                let px = x(center as f64);
                ctx.move_to(px, y(low));
                ctx.line_to(px, y(high));
            }
        }
        ctx.stroke();
    }
    ctx.set_font("11px system-ui");
    ctx.set_fill_style_str("#788b80");
    ctx.set_text_align("left");
    let _ = ctx.fill_text(&format!("{}", data.maximum), 5., plot_top - 4.);
    let _ = ctx.fill_text(&format!("{}", data.minimum), 5., plot_bottom + 13.);
    let tick_step = ((last.saturating_sub(first) as f64) / (width / 95.).max(1.))
        .ceil()
        .max(1.) as usize;
    let mut tick = first.div_ceil(tick_step) * tick_step;
    while tick < last {
        let px = x(tick as f64);
        ctx.set_fill_style_str("#788b80");
        let _ = ctx.fill_text(&format!("{tick}"), px, height - 7.);
        tick = tick.saturating_add(tick_step);
    }
    let (first_call, last_call) = visible_calls(trace, scroll, scale, width);
    let max_labels = (width / 16.).floor().max(1.) as usize + 2;
    let stride = (last_call - first_call).div_ceil(max_labels).max(1);
    let mut previous_x = f64::NEG_INFINITY;
    for i in (first_call..last_call).step_by(stride) {
        let px = x(trace.peak_positions[i] as f64);
        if px - previous_x < 14. {
            continue;
        }
        previous_x = px;
        ctx.set_text_align("center");
        ctx.set_font("10px system-ui");
        ctx.set_fill_style_str("#7c8b82");
        let _ = ctx.fill_text(&format!("{}", i + 1), px, 17.);
        let base = trace.sequence.as_str().as_bytes()[i];
        let color = match base {
            b'A' => COLORS[0],
            b'C' => COLORS[1],
            b'G' => COLORS[2],
            b'T' => COLORS[3],
            _ => "#52695d",
        };
        ctx.set_font("bold 15px monospace");
        ctx.set_fill_style_str(color);
        let _ = ctx.fill_text(&(base as char).to_string(), px, 37.);
        ctx.set_font("9px system-ui");
        ctx.set_fill_style_str("#7c8b82");
        if let Some(QualityEvidence::AbifBasecallerScores(values)) = &trace.qualities {
            let _ = ctx.fill_text(&values[i].to_string(), px, 54.);
        }
        ctx.set_stroke_style_str("#e6ece7");
        ctx.set_line_width(0.6);
        ctx.begin_path();
        ctx.move_to(px, 62.);
        ctx.line_to(px, plot_bottom);
        ctx.stroke();
    }
    Ok(())
}

#[component]
pub fn TraceView(
    molecule: Memo<Molecule>,
    source_scope: Memo<String>,
    selection: RwSignal<Option<(usize, usize)>>,
) -> impl IntoView {
    let data = Memo::new(move |_| {
        traces::trace_from_molecule(&molecule.get())
            .map(|trace| Arc::new(PlotData::new(trace)))
            .map_err(|e| e.to_string())
    });
    let export_trace = Memo::new(move |_| data.get().ok().map(|d| d.trace.clone()));
    // Host scope includes workspace epoch and exact document/event identity;
    // equal sequence attachments in separate documents must not share reviews.
    let export_revision =
        Memo::new(move |_| format!("{}:{}", source_scope.get(), molecule.with(|m| m.revision())));
    let canvas = NodeRef::<html::Canvas>::new();
    let container = NodeRef::<html::Div>::new();
    let zoom = RwSignal::new(2f64);
    let scroll = RwSignal::new(0f64);
    let viewport_width = RwSignal::new(500f64);
    let resize_epoch = RwSignal::new(0u32);
    let render_error = RwSignal::new(String::new());
    let resize = StoredValue::new_local(None::<ResizeListener>);
    let observed = StoredValue::new_local(None::<web_sys::Element>);
    Effect::new(move |_| {
        if let Some(element) = container.get() {
            let target: web_sys::Element = element.unchecked_ref::<web_sys::Element>().clone();
            if observed.with_value(|old| old.as_ref() == Some(&target)) {
                return;
            }
            resize.set_value(
                ResizeListener::new(&target, move || {
                    let _ = resize_epoch.try_update(|n| *n = n.wrapping_add(1));
                })
                .ok(),
            );
            observed.set_value(Some(target));
        }
    });
    Effect::new(move |_| {
        let Ok(data) = data.get() else {
            return;
        };
        let width = container
            .get_untracked()
            .map(|e| e.client_width() as f64)
            .unwrap_or(500.)
            .max(200.);
        let peaks = &data.trace.peak_positions;
        let spacing = if peaks.len() > 1 {
            (peaks[peaks.len() - 1] - peaks[0]) as f64 / (peaks.len() - 1) as f64
        } else {
            10.
        };
        let start = (20. / spacing.max(1.))
            .max((width - 2. * PAD) / data.trace.channels.len().min(200) as f64)
            .clamp(0.125, data.max_zoom());
        zoom.set(start);
        scroll.set(0.);
        schedule_scroll(container, scroll, 0.);
    });
    Effect::new(move |_| {
        let Ok(data) = data.get() else {
            return;
        };
        let selected = selection.get();
        let scale = zoom.get();
        let left = scroll.get();
        let _ = resize_epoch.get();
        let (Some(canvas), Some(element)) = (canvas.get(), container.get()) else {
            return;
        };
        let width = (element.client_width() as f64).max(200.);
        let height = (element.client_height() as f64).clamp(220., 640.);
        if viewport_width.get_untracked() != width {
            viewport_width.set(width);
        }
        let result = paint(
            canvas.unchecked_ref(),
            &data,
            width,
            height,
            scale,
            left,
            selected,
        );
        render_error.set(result.err().unwrap_or_default());
    });
    Effect::new(move |_| {
        let Ok(data) = data.get() else {
            return;
        };
        let selected = selection.get();
        let scale = zoom.get();
        if let (Some((a, _)), Some(element)) = (selected, container.get()) {
            if let Some(&peak) = data.trace.peak_positions.get(a) {
                let x = PAD + peak as f64 * scale;
                let left = element.scroll_left() as f64;
                let width = element.client_width() as f64;
                if x < left + PAD || x > left + width - PAD {
                    schedule_scroll(
                        container,
                        scroll,
                        (x - width / 2.)
                            .clamp(0., (virtual_width(&data, scale, width) - width).max(0.)),
                    );
                }
            }
        }
    });
    let zoom_by = move |factor: f64| {
        let Ok(data) = data.get_untracked() else {
            return;
        };
        let width = viewport_width.get_untracked();
        let current = zoom.get_untracked();
        let fit = ((width - 2. * PAD) / data.trace.channels.len() as f64).max(0.00001);
        let next = (current * factor).clamp(fit.min(0.125), data.max_zoom());
        let focus = focus_sample(
            &data,
            selection.get_untracked(),
            scroll.get_untracked(),
            current,
            width,
        );
        change_zoom(container, zoom, scroll, next, focus);
    };
    let fit = move |_| {
        if let Ok(data) = data.get_untracked() {
            let width = viewport_width.get_untracked();
            let scale = ((width - 2. * PAD) / data.trace.channels.len() as f64)
                .clamp(0.00001, data.max_zoom());
            change_zoom(
                container,
                zoom,
                scroll,
                scale,
                data.trace.channels.len() as f64 / 2.,
            );
        }
    };
    let rows = Memo::new(move |_| {
        let Ok(data) = data.get() else {
            return (Vec::<usize>::new(), 0usize);
        };
        let (first, last) =
            visible_calls(&data.trace, scroll.get(), zoom.get(), viewport_width.get());
        (
            (first..last.min(first + 200)).collect::<Vec<_>>(),
            last - first,
        )
    });
    let selected_text = move || {
        let Ok(data) = data.get() else {
            return String::new();
        };
        let Some((a, b)) = selection.get() else {
            return "Select a called peak to inspect its stored signal and quality evidence."
                .into();
        };
        let trace = &data.trace;
        if a >= trace.sequence.len() || b <= a {
            return "Selection is outside the trace calls.".into();
        }
        let b = b.min(trace.sequence.len());
        let sample = trace.peak_positions[a] as usize;
        format!(
            "Selected bases {}–{}; first call {} at display sample {} (zero-based). {}; signals A {}, C {}, G {}, T {}. No calibrated quality is inferred.",
            a + 1,
            b,
            trace.sequence.as_str().as_bytes()[a] as char,
            sample,
            evidence(trace, a),
            trace.channels.a[sample],
            trace.channels.c[sample],
            trace.channels.g[sample],
            trace.channels.t[sample]
        )
    };
    view! {
        <div class="trace-view"><style>{STYLE}</style>
            <Show when=move||data.get().is_err()><p class="trace-error" role="alert">{move||data.get().err().unwrap_or_default()}" Original evidence remains attached; a stale current sequence is never silently substituted into the chromatogram."</p></Show>
            <Show when=move||data.get().is_ok()>
                <div class="trace-controls"><div class="trace-legend" aria-label="Trace channel legend"><span style:color=COLORS[0]>"A"</span><span style:color=COLORS[1]>"C"</span><span style:color=COLORS[2]>"G"</span><span style:color=COLORS[3]>"T"</span></div><crate::scf_export_view::ScfExportControls document_revision=export_revision trace=export_trace/><button title="Zoom trace out" aria-label="Zoom trace out" on:click=move|_|zoom_by(1./1.5)>"−"</button><button title="Fit entire trace" on:click=fit>"Fit trace"</button><button title="Zoom trace in" aria-label="Zoom trace in" on:click=move|_|zoom_by(1.5)>"+"</button></div>
                <p class="trace-note" id="trace-help">{move||data.get().ok().map(|d|format!("{} calls · {} samples/channel. {} Signals share one amplitude scale; negatives are retained. Scroll to pan; arrow keys select calls. Overview traces use min/max envelopes.",d.trace.sequence.len(),d.trace.channels.len(),quality_label(&d.trace))).unwrap_or_default()}</p>
                <Show when=move||!render_error.get().is_empty()><p class="trace-error" role="alert">{move||render_error.get()}</p></Show>
                <div class="trace-scroller" node_ref=container tabindex="0" role="region" aria-label="Horizontally scrollable chromatogram" on:scroll=move|_|if let Some(element)=container.get(){scroll.set(element.scroll_left() as f64);}>
                    <div class="trace-spacer" style:width=move||data.get().ok().map(|d|format!("{}px",virtual_width(&d,zoom.get(),viewport_width.get()))).unwrap_or_else(||"100%".into())>
                        <canvas node_ref=canvas class="trace-canvas" tabindex="0" aria-label="Chromatogram. Click a called peak or use left/right, Home and End keys to select calls." aria-describedby="trace-help" on:keydown=move|event|{
                            let Ok(data)=data.get_untracked()else{return;};let n=data.trace.sequence.len();let current=selection.get_untracked().map(|(a,_)|a.min(n-1)).unwrap_or(0);
                            let next=match event.key().as_str(){"ArrowRight"=>Some((current+1).min(n-1)),"ArrowLeft"=>Some(current.saturating_sub(1)),"Home"=>Some(0),"End"=>Some(n-1),_=>None};if let Some(next)=next{event.prevent_default();selection.set(Some((next,next+1)));}
                        } on:click=move|event|{
                            let Ok(data)=data.get_untracked()else{return;};let Some(canvas)=canvas.get()else{return;};let rect=canvas.get_bounding_client_rect();let sample=((event.client_x() as f64-rect.left()+scroll.get_untracked()-PAD)/zoom.get_untracked()).clamp(0.,data.trace.channels.len().saturating_sub(1) as f64);
                            let peaks=&data.trace.peak_positions;let next=peaks.partition_point(|&p|(p as f64)<sample);let index=if next==0{0}else if next>=peaks.len(){peaks.len()-1}else if sample-(peaks[next-1] as f64)<(peaks[next] as f64)-sample{next-1}else{next};
                            let range=if event.shift_key(){selection.get_untracked().map(|(a,_)|(a.min(index),a.max(index)+1)).unwrap_or((index,index+1))}else{(index,index+1)};selection.set(Some(range));
                        }></canvas>
                    </div>
                </div>
                <p class="trace-selection" aria-live="polite">{selected_text}</p>
                <details class="trace-access"><summary>"Visible called bases and raw evidence"</summary><div class="trace-table-wrap"><table><caption>{move||format!("{} visible calls; showing at most 200. Calls are one-based, sample indexes zero-based. Stored scores retain source semantics.",rows.get().1)}</caption><thead><tr><th>"Call"</th><th>"Base"</th><th>"Sample"</th><th>"Stored quality evidence"</th></tr></thead><tbody>{move||data.get().ok().map(|data|rows.get().0.into_iter().map(|i|{let trace=&data.trace;view!{<tr><td><button on:click=move|_|selection.set(Some((i,i+1)))>{i+1}</button></td><td>{(trace.sequence.as_str().as_bytes()[i] as char).to_string()}</td><td>{trace.peak_positions[i]}</td><td>{evidence(trace,i)}</td></tr>}}).collect_view())}</tbody></table></div></details>
                <div class="sr-only" aria-live="polite">{move||selection.get().and_then(|(a,b)|data.get().ok().and_then(|d|d.trace.sequence.as_str().get(a..b).map(|s|format!("Selected bases {} to {}: {}",a+1,b,s)))).unwrap_or_default()}</div>
            </Show>
        </div>
    }
}
