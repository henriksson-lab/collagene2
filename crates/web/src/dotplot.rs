//! Revision-bound window dot plots with retained WebGL anchors and exact hit inspection.
#[path = "dotplot_projection.rs"]
mod projection;
use crate::{renderer::ResizeListener, worker::CachedWorker};
use collagene_core::{
    dotplot::{self, DotHit, DotPlotOptions, DotPlotReport, Span},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule, complement};
use leptos::{html, prelude::*};
use std::{
    cell::Cell,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::spawn_local;
use web_sys::{
    CanvasRenderingContext2d as Ctx, HtmlCanvasElement, WebGl2RenderingContext as Gl, WebGlBuffer,
    WebGlProgram,
};

#[derive(Clone, PartialEq, Eq)]
pub struct QueryInfo {
    pub identity: String,
    pub revision: String,
    pub name: String,
    pub length: usize,
    pub alphabet: Alphabet,
}
#[derive(Clone)]
struct BoundPlot {
    report: Arc<DotPlotReport>,
    reference: Arc<Molecule>,
    query: Arc<Molecule>,
    reference_identity: String,
    query_identity: String,
    epoch: u64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
struct Viewport {
    x: f64,
    y: f64,
    width: f64,
    height: f64,
}
impl Viewport {
    fn all(r: &DotPlotReport) -> Self {
        Self {
            x: 0.,
            y: 0.,
            width: r.reference_length as f64,
            height: r.query_length as f64,
        }
    }
    fn clamp(self, r: &DotPlotReport) -> Self {
        let width = self.width.clamp(1., r.reference_length as f64);
        let height = self.height.clamp(1., r.query_length as f64);
        Self {
            x: self.x.clamp(0., r.reference_length as f64 - width),
            y: self.y.clamp(0., r.query_length as f64 - height),
            width,
            height,
        }
    }
    fn zoom(self, r: &DotPlotReport, factor: f64, anchor: (f64, f64)) -> Self {
        let width = (self.width * factor).clamp(1., r.reference_length as f64);
        let height = (self.height * factor).clamp(1., r.query_length as f64);
        Self {
            x: self.x + (self.width - width) * anchor.0,
            y: self.y + (self.height - height) * anchor.1,
            width,
            height,
        }
        .clamp(r)
    }
}
const LEFT: f64 = 52.;
const TOP: f64 = 18.;
const RIGHT: f64 = 14.;
const BOTTOM: f64 = 44.;
fn unit(a: Alphabet) -> &'static str {
    match a {
        Alphabet::Protein => "residues",
        Alphabet::Dna | Alphabet::Rna => "nt",
    }
}
fn spans(values: &[Span]) -> String {
    values
        .iter()
        .map(|s| format!("{}–{}", s.start + 1, s.end))
        .collect::<Vec<_>>()
        .join(" + ")
}
fn hit_label(h: &DotHit) -> String {
    format!(
        "Window {} · R {} · Q {} · {} · {}/{} compatible, {}/{} certain",
        h.id,
        spans(&h.reference_spans),
        spans(&h.query_spans),
        if h.reverse {
            "reverse complement"
        } else {
            "forward"
        },
        h.compatible_matches,
        h.window,
        h.certain_matches,
        h.window
    )
}
fn preview(m: &Molecule, start: usize, window: usize, reverse: bool) -> String {
    let bases = m.sequence.as_str().as_bytes();
    let mut text = String::new();
    for i in 0..window.min(200) {
        let at = if reverse {
            (start + bases.len() - i % bases.len()) % bases.len()
        } else {
            (start + i) % bases.len()
        };
        let b = bases[at];
        text.push(if reverse {
            complement(b, m.sequence.alphabet())
        } else {
            b
        } as char);
    }
    if window > 200 {
        text.push_str("… (first 200 symbols)");
    }
    text
}
fn whole<T: std::str::FromStr>(s: &str, label: &str) -> Result<T, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {label}."))
}
struct Json(Vec<u8>);
impl std::io::Write for Json {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > 32 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "Dot-plot JSON exceeds this browser's 32 MiB output bound",
            ));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
struct Gpu {
    gl: Gl,
    program: WebGlProgram,
    buffer: WebGlBuffer,
    count: i32,
}
impl Gpu {
    fn new(canvas: &HtmlCanvasElement, r: &DotPlotReport) -> Result<Self, String> {
        let gl = canvas
            .get_context("webgl2")
            .map_err(|_| "WebGL2 unavailable")?
            .ok_or("WebGL2 unavailable")?
            .dyn_into::<Gl>()
            .map_err(|_| "WebGL2 unavailable")?;
        let shader = |kind, source: &str| -> Result<_, String> {
            let s = gl.create_shader(kind).ok_or("Cannot create plot shader")?;
            gl.shader_source(&s, source);
            gl.compile_shader(&s);
            if !gl
                .get_shader_parameter(&s, Gl::COMPILE_STATUS)
                .as_bool()
                .unwrap_or(false)
            {
                return Err(gl
                    .get_shader_info_log(&s)
                    .unwrap_or_else(|| "Cannot compile plot shader".into()));
            }
            Ok(s)
        };
        let vertex = shader(Gl::VERTEX_SHADER, projection::VERTEX_SHADER)?;
        let fragment = shader(
            Gl::FRAGMENT_SHADER,
            "#version 300 es\nprecision mediump float;in vec3 v_color;out vec4 out_color;void main(){out_color=vec4(v_color,0.85);}",
        )?;
        let program = gl.create_program().ok_or("Cannot create plot program")?;
        gl.attach_shader(&program, &vertex);
        gl.attach_shader(&program, &fragment);
        gl.link_program(&program);
        gl.delete_shader(Some(&vertex));
        gl.delete_shader(Some(&fragment));
        if !gl
            .get_program_parameter(&program, Gl::LINK_STATUS)
            .as_bool()
            .unwrap_or(false)
        {
            return Err("Cannot link plot program".into());
        }
        let buffer = gl
            .create_buffer()
            .ok_or("Cannot create plot vertex buffer")?;
        let mut data = Vec::with_capacity(r.hits.len() * 3);
        for h in &r.hits {
            data.extend_from_slice(&[
                h.reference_start as f32,
                h.query_start as f32,
                match (h.reverse, h.certain) {
                    (false, true) => 0.,
                    (true, true) => 1.,
                    (false, false) => 2.,
                    (true, false) => 3.,
                },
            ]);
        }
        gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&buffer));
        gl.buffer_data_with_array_buffer_view(
            Gl::ARRAY_BUFFER,
            &js_sys::Float32Array::from(data.as_slice()),
            Gl::STATIC_DRAW,
        );
        Ok(Self {
            gl,
            program,
            buffer,
            count: r.hits.len() as i32,
        })
    }
    fn draw(
        &self,
        canvas: &HtmlCanvasElement,
        _r: &DotPlotReport,
        v: Viewport,
        width: f64,
        height: f64,
        dpr: f64,
    ) {
        let gl = &self.gl;
        gl.viewport(0, 0, canvas.width() as i32, canvas.height() as i32);
        gl.disable(Gl::SCISSOR_TEST);
        gl.clear_color(1., 1., 1., 1.);
        gl.clear(Gl::COLOR_BUFFER_BIT);
        gl.use_program(Some(&self.program));
        gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&self.buffer));
        let position = gl.get_attrib_location(&self.program, "a_anchor") as u32;
        gl.enable_vertex_attrib_array(position);
        gl.vertex_attrib_pointer_with_i32(position, 3, Gl::FLOAT, false, 0, 0);
        let camera = projection::camera(v.x, v.y, v.width, v.height);
        gl.uniform2f(
            gl.get_uniform_location(&self.program, "u_origin_integer")
                .as_ref(),
            camera.integer[0],
            camera.integer[1],
        );
        gl.uniform2f(
            gl.get_uniform_location(&self.program, "u_origin_fraction")
                .as_ref(),
            camera.fraction[0],
            camera.fraction[1],
        );
        gl.uniform2f(
            gl.get_uniform_location(&self.program, "u_span").as_ref(),
            camera.span[0],
            camera.span[1],
        );
        gl.uniform2f(
            gl.get_uniform_location(&self.program, "u_size").as_ref(),
            width as f32,
            height as f32,
        );
        gl.uniform2f(
            gl.get_uniform_location(&self.program, "u_plot").as_ref(),
            (width - LEFT - RIGHT) as f32,
            (height - TOP - BOTTOM) as f32,
        );
        gl.uniform1f(
            gl.get_uniform_location(&self.program, "u_dpr").as_ref(),
            dpr as f32,
        );
        gl.enable(Gl::SCISSOR_TEST);
        gl.scissor(
            (LEFT * dpr) as i32,
            (BOTTOM * dpr) as i32,
            ((width - LEFT - RIGHT) * dpr) as i32,
            ((height - TOP - BOTTOM) * dpr) as i32,
        );
        gl.draw_arrays(Gl::POINTS, 0, self.count);
        gl.disable(Gl::SCISSOR_TEST);
    }
}
impl Drop for Gpu {
    fn drop(&mut self) {
        self.gl.delete_buffer(Some(&self.buffer));
        self.gl.delete_program(Some(&self.program));
    }
}
fn colors(reverse: bool, certain: bool) -> &'static str {
    match (reverse, certain) {
        (false, true) => "#2c7b5c",
        (true, true) => "#995087",
        (false, false) => "#ba8529",
        (true, false) => "#d1662e",
    }
}
fn paint_axes(
    ctx: &Ctx,
    r: &DotPlotReport,
    v: Viewport,
    width: f64,
    height: f64,
    selected: Option<usize>,
    fallback: bool,
) {
    let w = width - LEFT - RIGHT;
    let h = height - TOP - BOTTOM;
    let x = |p: f64| LEFT + (p - v.x) / v.width * w;
    let y = |p: f64| TOP + (p - v.y) / v.height * h;
    ctx.clear_rect(0., 0., width, height);
    if fallback {
        ctx.set_fill_style_str("white");
        ctx.fill_rect(0., 0., width, height);
        let mut bins = vec![0u8; 64 * 64];
        for hit in &r.hits {
            let px = (hit.reference_start as f64 + 0.5 - v.x) / v.width;
            let py = (hit.query_start as f64 + 0.5 - v.y) / v.height;
            if (0.0..1.0).contains(&px) && (0.0..1.0).contains(&py) {
                let category = match (hit.reverse, hit.certain) {
                    (false, true) => 0,
                    (true, true) => 1,
                    (false, false) => 2,
                    (true, false) => 3,
                };
                bins[(py * 64.) as usize * 64 + (px * 64.) as usize] |= 1 << category;
            }
        }
        for (i, bits) in bins.into_iter().enumerate() {
            if bits != 0 {
                for category in 0..4 {
                    if bits & (1 << category) != 0 {
                        ctx.set_fill_style_str(colors(category % 2 == 1, category < 2));
                        ctx.fill_rect(
                            LEFT + (i % 64) as f64 * w / 64.,
                            TOP + (i / 64) as f64 * h / 64.,
                            (w / 64.).max(1.),
                            (h / 64.).max(1.),
                        );
                    }
                }
            }
        }
    }
    ctx.set_stroke_style_str("#d8e4dc");
    ctx.set_line_width(1.);
    ctx.stroke_rect(LEFT, TOP, w, h);
    ctx.set_font("11px system-ui");
    ctx.set_fill_style_str("#597366");
    for i in 0..=4 {
        let t = i as f64 / 4.;
        let a = v.x + t * v.width;
        let b = v.y + t * v.height;
        ctx.set_text_align("center");
        let _ = ctx.fill_text(
            &format!("{}", (a.floor() as usize).min(r.reference_length - 1) + 1),
            LEFT + t * w,
            height - 23.,
        );
        ctx.set_text_align("right");
        let _ = ctx.fill_text(
            &format!("{}", (b.floor() as usize).min(r.query_length - 1) + 1),
            LEFT - 7.,
            TOP + t * h + 4.,
        );
    }
    ctx.set_text_align("center");
    let _ = ctx.fill_text(
        &format!("Reference ({})", unit(r.alphabet)),
        LEFT + w / 2.,
        height - 5.,
    );
    ctx.save();
    let _ = ctx.translate(12., TOP + h / 2.);
    let _ = ctx.rotate(-std::f64::consts::FRAC_PI_2);
    let _ = ctx.fill_text(&format!("Query ({})", unit(r.alphabet)), 0., 0.);
    ctx.restore();
    if let Some(hit) = selected.and_then(|i| r.hits.get(i)) {
        let px = x(hit.reference_start as f64 + 0.5);
        let py = y(hit.query_start as f64 + 0.5);
        ctx.save();
        ctx.begin_path();
        ctx.rect(LEFT, TOP, w, h);
        ctx.clip();
        ctx.set_stroke_style_str("#1860a5");
        ctx.set_line_width(1.5);
        ctx.begin_path();
        let _ = ctx.arc(px, py, 6., 0., std::f64::consts::TAU);
        ctx.stroke();
        // Split a selected window at each circular-origin crossing, never joining
        // opposite plot edges. This is a window trajectory, not a gapped alignment.
        let mut cuts = vec![0, hit.window];
        let rx = r.reference_length - hit.reference_start;
        let qx = if hit.reverse {
            hit.query_start + 1
        } else {
            r.query_length - hit.query_start
        };
        for cut in [rx, qx] {
            if cut > 0 && cut < hit.window {
                cuts.push(cut);
            }
        }
        cuts.sort_unstable();
        cuts.dedup();
        for pair in cuts.windows(2) {
            let a = pair[0];
            let b = pair[1] - 1;
            let coord = |i: usize| {
                (
                    (hit.reference_start + i) % r.reference_length,
                    if hit.reverse {
                        (hit.query_start + r.query_length - i % r.query_length) % r.query_length
                    } else {
                        (hit.query_start + i) % r.query_length
                    },
                )
            };
            let (a0, a1) = coord(a);
            let (b0, b1) = coord(b);
            ctx.begin_path();
            ctx.move_to(x(a0 as f64 + 0.5), y(a1 as f64 + 0.5));
            ctx.line_to(x(b0 as f64 + 0.5), y(b1 as f64 + 0.5));
            ctx.stroke();
        }
        ctx.restore();
    }
}

#[component]
fn PlotCanvas(
    data: RwSignal<Option<Arc<BoundPlot>>>,
    selected: RwSignal<Option<usize>>,
    viewport: RwSignal<Viewport>,
    expanded: RwSignal<bool>,
) -> impl IntoView {
    let canvas = NodeRef::<html::Canvas>::new();
    let overlay = NodeRef::<html::Canvas>::new();
    let container = NodeRef::<html::Div>::new();
    let resize = StoredValue::new_local(None::<ResizeListener>);
    let resize_epoch = RwSignal::new(0u32);
    let gpu = StoredValue::new_local(None::<Gpu>);
    let gpu_revision = StoredValue::new(String::new());
    let fallback = RwSignal::new(false);
    let queued = StoredValue::new_local(Rc::new(Cell::new(false)));
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, Ordering::Relaxed));
    Effect::new(move |_| {
        if let Some(element) = container.get() {
            if resize.with_value(|r| r.is_none()) {
                resize.set_value(
                    ResizeListener::new(element.unchecked_ref(), move || {
                        let _ = resize_epoch.try_update(|n| *n = n.wrapping_add(1));
                    })
                    .ok(),
                );
            }
        }
    });
    Effect::new(move |_| {
        let _ = (
            viewport.get(),
            selected.get(),
            resize_epoch.get(),
            expanded.get(),
        );
        let bound = data.get();
        let key = bound
            .as_ref()
            .map(|b| {
                format!(
                    "{}:{}:{}",
                    b.report.reference_revision, b.report.query_revision, b.report.options_revision
                )
            })
            .unwrap_or_default();
        if key != gpu_revision.get_value()
            || (canvas.get().is_some()
                && bound.is_some()
                && gpu.with_value(|g| g.is_none())
                && !fallback.get_untracked())
        {
            gpu.set_value(None);
            gpu_revision.set_value(key);
            if let (Some(c), Some(b)) = (canvas.get(), bound.as_ref()) {
                match Gpu::new(&c, &b.report) {
                    Ok(value) => {
                        gpu.set_value(Some(value));
                        fallback.set(false);
                    }
                    Err(_) => fallback.set(true),
                }
            }
        }
        if queued.with_value(|q| q.replace(true)) {
            return;
        }
        let alive = alive.clone();
        let frame = Closure::once_into_js(move || {
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            queued.with_value(|q| q.set(false));
            let (Some(c), Some(a), Some(wrapper), Some(bound)) = (
                canvas.get_untracked(),
                overlay.get_untracked(),
                container.get_untracked(),
                data.get_untracked(),
            ) else {
                return;
            };
            let width = f64::from(wrapper.client_width()).max(160.);
            let height = f64::from(wrapper.client_height()).max(220.);
            let dpr = web_sys::window()
                .map(|w| w.device_pixel_ratio())
                .unwrap_or(1.)
                .clamp(1., 2.);
            for layer in [
                c.unchecked_ref::<HtmlCanvasElement>(),
                a.unchecked_ref::<HtmlCanvasElement>(),
            ] {
                let w = (width * dpr).round() as u32;
                let h = (height * dpr).round() as u32;
                if layer.width() != w {
                    layer.set_width(w);
                }
                if layer.height() != h {
                    layer.set_height(h);
                }
            }
            let v = viewport.get_untracked().clamp(&bound.report);
            gpu.with_value(|g| {
                if let Some(g) = g {
                    g.draw(&c, &bound.report, v, width, height, dpr);
                }
            });
            if let Ok(Some(ctx)) = a.get_context("2d") {
                let ctx: Ctx = ctx.unchecked_into();
                let _ = ctx.set_transform(dpr, 0., 0., dpr, 0., 0.);
                paint_axes(
                    &ctx,
                    &bound.report,
                    v,
                    width,
                    height,
                    selected.get_untracked(),
                    fallback.get_untracked(),
                );
            }
            let _ = c.set_attribute("data-painted", "true");
            let _ = c.set_attribute("data-hit-count", &bound.report.hits.len().to_string());
            let _ = c.set_attribute("data-view-x", &v.x.to_string());
            let _ = c.set_attribute("data-view-y", &v.y.to_string());
            let _ = c.set_attribute("data-view-width", &v.width.to_string());
            let _ = c.set_attribute(
                "data-renderer",
                if fallback.get_untracked() {
                    "canvas-aggregate"
                } else {
                    "webgl2"
                },
            );
        });
        if let Some(w) = web_sys::window() {
            let _ = w.request_animation_frame(frame.unchecked_ref());
        }
    });
    let drag = StoredValue::new_local(None::<(f64, f64, Viewport)>);
    let moved = StoredValue::new_local(false);
    let pick = Callback::new(move |(client_x, client_y): (f64, f64)| {
        let (Some(c), Some(b)) = (canvas.get_untracked(), data.get_untracked()) else {
            return;
        };
        let rect = c.get_bounding_client_rect();
        let v = viewport.get_untracked();
        let (px, py) = (client_x - rect.left(), client_y - rect.top());
        let (w, h) = (rect.width() - LEFT - RIGHT, rect.height() - TOP - BOTTOM);
        let mut nearest = None;
        let mut best = 100.;
        for (i, hit) in b.report.hits.iter().enumerate() {
            let x = LEFT + (hit.reference_start as f64 + 0.5 - v.x) / v.width * w;
            let y = TOP + (hit.query_start as f64 + 0.5 - v.y) / v.height * h;
            let distance = (x - px).powi(2) + (y - py).powi(2);
            if distance < best {
                nearest = Some(i);
                best = distance;
            }
        }
        selected.set(nearest);
    });
    view! {<div class="dotplot-plot" node_ref=container><canvas class="dotplot-canvas" node_ref=canvas tabindex="0" aria-label="Window dot plot. Drag to pan, scroll to zoom around the pointer, click a window anchor to inspect it. Arrow keys pan; plus and minus zoom; Home fits both sequences."
    on:pointerdown=move|e|{if let Some(c)=canvas.get(){let _=c.set_pointer_capture(e.pointer_id());}drag.set_value(Some((e.client_x()as f64,e.client_y()as f64,viewport.get_untracked())));moved.set_value(false);}
    on:pointermove=move|e|{if let(Some((x,y,v)),Some(c),Some(b))=(drag.get_value(),canvas.get(),data.get_untracked()){let dx=e.client_x()as f64-x;let dy=e.client_y()as f64-y;if dx.abs()+dy.abs()>3.{moved.set_value(true);}let rect=c.get_bounding_client_rect();viewport.set(Viewport{x:v.x-dx/(rect.width()-LEFT-RIGHT)*v.width,y:v.y-dy/(rect.height()-TOP-BOTTOM)*v.height,..v}.clamp(&b.report));}}
    on:pointerup=move|e|{let was=drag.get_value().is_some();drag.set_value(None);if was&&!moved.get_value(){pick.run((e.client_x()as f64,e.client_y()as f64));}}
    on:pointercancel=move|_|drag.set_value(None)
    on:wheel=move|e|{e.prevent_default();if let(Some(c),Some(b))=(canvas.get(),data.get_untracked()){let r=c.get_bounding_client_rect();let anchor=(((e.client_x()as f64-r.left()-LEFT)/(r.width()-LEFT-RIGHT)).clamp(0.,1.),((e.client_y()as f64-r.top()-TOP)/(r.height()-TOP-BOTTOM)).clamp(0.,1.));viewport.update(|v|*v=v.zoom(&b.report,if e.delta_y()>0.{1.2}else{1./1.2},anchor));}}
    on:keydown=move|e|{let Some(b)=data.get_untracked()else{return;};let v=viewport.get_untracked();let next=match e.key().as_str(){"ArrowRight"=>Some(Viewport{x:v.x+v.width*0.1,..v}.clamp(&b.report)),"ArrowLeft"=>Some(Viewport{x:v.x-v.width*0.1,..v}.clamp(&b.report)),"ArrowDown"=>Some(Viewport{y:v.y+v.height*0.1,..v}.clamp(&b.report)),"ArrowUp"=>Some(Viewport{y:v.y-v.height*0.1,..v}.clamp(&b.report)),"+"|"="=>Some(v.zoom(&b.report,0.5,(0.5,0.5))),"-"=>Some(v.zoom(&b.report,2.,(0.5,0.5))),"Home"=>Some(Viewport::all(&b.report)),_=>None};if let Some(v)=next{e.prevent_default();viewport.set(v);}}></canvas><canvas class="dotplot-axes" node_ref=overlay aria-hidden="true"></canvas></div><Show when=move||fallback.get()><p class="helper">"WebGL2 is unavailable. This drawing combines all visible retained anchors into a 64×64 grid; the hit list and exports keep exact coordinates."</p></Show>}
}

#[component]
pub fn DotPlotPanel(
    molecule: Memo<Molecule>,
    source_identity: Memo<String>,
    active_index: Memo<usize>,
    query_options: Memo<Vec<(usize, String)>>,
    query_info: Callback<usize, Option<QueryInfo>>,
    get_query: Callback<usize, Option<Molecule>>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<(bool, usize, Vec<(usize, usize)>)>,
) -> impl IntoView {
    let query_index = RwSignal::new(active_index.get_untracked());
    let info = Memo::new(move |_| query_info.run(query_index.get()));
    let window = RwSignal::new("15".to_string());
    let mismatches = RwSignal::new("0".to_string());
    let stride = RwSignal::new("1".to_string());
    let both = RwSignal::new(true);
    let circular = RwSignal::new(true);
    let work = RwSignal::new("20000000".to_string());
    let hits = RwSignal::new("100000".to_string());
    let index = RwSignal::new("1000000".to_string());
    let protein = Memo::new(move |_| molecule.with(|m| m.sequence.alphabet() == Alphabet::Protein));
    Effect::new(move |_| {
        let _ = (source_identity.get(), epoch.get(), active_index.get());
        let alphabet = molecule.with(|m| m.sequence.alphabet());
        if info.get_untracked().is_none_or(|q| q.alphabet != alphabet) {
            query_index.set(active_index.get_untracked());
        }
    });
    let pair = Memo::new(move |_| (source_identity.get(), info.get().map(|i| i.identity)));
    Effect::new(move |_| {
        let _ = pair.get();
        let length = molecule
            .with(|m| m.sequence.len())
            .min(info.get().map_or(15, |i| i.length));
        window.set(length.min(15).max(1).to_string());
        mismatches.set("0".into());
    });
    let options = Memo::new(move |_| -> Result<DotPlotOptions, String> {
        let w: usize = whole(&window.get(), "window length")?;
        let miss: usize = whole(&mismatches.get(), "allowed mismatches")?;
        if miss > w {
            return Err("Allowed mismatches cannot exceed window length.".into());
        }
        let o = DotPlotOptions {
            window: w,
            min_matches: w - miss,
            stride: whole(&stride.get(), "sampling stride")?,
            both_strands: both.get() && !protein.get(),
            wrap_circular: circular.get(),
            max_work: whole(&work.get(), "work budget")?,
            max_hits: whole(&hits.get(), "hit budget")?,
            max_index_entries: whole(&index.get(), "index budget")?,
        };
        Ok(o)
    });
    let data = RwSignal::new(None::<Arc<BoundPlot>>);
    let selected = RwSignal::new(None::<usize>);
    let viewport = RwSignal::new(Viewport {
        x: 0.,
        y: 0.,
        width: 1.,
        height: 1.,
    });
    let expanded = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = active.get_value() {
            w.close();
        }
    });
    Effect::new(move |_| {
        let _ = (
            source_identity.get(),
            info.get(),
            options.get(),
            epoch.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
        busy.set(false);
        data.set(None);
        selected.set(None);
        page.set(0);
        status.set(String::new());
    });
    let guard = Callback::new(move |()| -> Result<Arc<BoundPlot>, String> {
        let b = data.get_untracked().ok_or("Calculate a dot plot first.")?;
        let q = query_info
            .run(query_index.get_untracked())
            .ok_or("The query document is unavailable.")?;
        if source_identity.get_untracked() != b.reference_identity
            || q.identity != b.query_identity
            || epoch.get_untracked() != b.epoch
            || options.get_untracked() != Ok(b.report.options.clone())
            || q.revision != b.report.query_revision
            || molecule.with_untracked(Molecule::revision) != b.report.reference_revision
        {
            return Err("Sources or options changed. Calculate this comparison again.".into());
        }
        Ok(b)
    });
    let start_alive = alive.clone();
    let run = Callback::new(move |()| {
        let opts = match options.get_untracked() {
            Ok(o) => o,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let Some(query) = get_query.run(query_index.get_untracked()) else {
            status.set("Choose an open query document.".into());
            return;
        };
        let Some(qinfo) = query_info.run(query_index.get_untracked()) else {
            return;
        };
        let source = molecule.get_untracked();
        let identity = source_identity.get_untracked();
        let source_epoch = epoch.get_untracked();
        let revision = source.revision();
        if let Some(w) = active.get_value() {
            w.close();
        }
        let worker = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        active.set_value(Some(worker.clone()));
        generation.update(|g| *g = g.wrapping_add(1));
        let token = generation.get_untracked();
        data.set(None);
        selected.set(None);
        page.set(0);
        busy.set(true);
        status.set("Comparing sequence windows in the background…".into());
        let alive = start_alive.clone();
        spawn_local(async move {
            let result = worker
                .ask(Request {
                    id: token,
                    revision,
                    job: AnalysisJob::DotPlot {
                        molecule: source.clone(),
                        query: query.clone(),
                        options: opts.clone(),
                    },
                })
                .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || source_identity.get_untracked() != identity
                || query_info.run(query_index.get_untracked()).as_ref() != Some(&qinfo)
                || epoch.get_untracked() != source_epoch
                || options.get_untracked() != Ok(opts.clone())
            {
                return;
            }
            busy.set(false);
            active.set_value(None);
            let result = result.and_then(|r| match r.result {
                AnalysisResult::DotPlot(report) => {
                    report
                        .validate_bindings_for(&source, &query, &opts)
                        .map_err(|e| e.to_string())?;
                    Ok(report)
                }
                AnalysisResult::Error(e) => Err(e),
                _ => Err("Unexpected dot-plot response.".into()),
            });
            match result {
                Ok(report) => {
                    viewport.set(Viewport::all(&report));
                    status.set(if !report.grid_complete{"Search stopped at a declared limit. Empty regions do not establish absence of similarity."}else if report.sampled{"Sampled window grid completed. Unsampled regions were not checked."}else{"Window comparison completed under the declared model."}.into());
                    data.set(Some(Arc::new(BoundPlot {
                        report: Arc::new(report),
                        reference: Arc::new(source),
                        query: Arc::new(query),
                        reference_identity: identity,
                        query_identity: qinfo.identity,
                        epoch: source_epoch,
                    })));
                }
                Err(e) => status.set(e),
            }
        });
    });
    let export_json = Callback::new(move |()| match guard.run(()) {
        Ok(b) => {
            let mut bytes = Json(vec![]);
            match serde_json::to_writer_pretty(&mut bytes, b.report.as_ref()) {
                Ok(()) => crate::app::download_with_completion(
                    "sequence-dotplot.json",
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
        }
        Err(e) => status.set(e),
    });
    let export_svg = Callback::new(move |()| match guard.run(()) {
        Ok(b) => match dotplot::export_svg(&b.report) {
            Ok(svg) => crate::app::download_with_completion(
                "sequence-dotplot.svg",
                "image/svg+xml",
                svg.as_bytes(),
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
        },
        Err(e) => status.set(e),
    });
    let focus = Callback::new(move |index: usize| {
        selected.set(Some(index));
        if let Some(b) = data.get_untracked() {
            if let Some(h) = b.report.hits.get(index) {
                let w = (h.window as f64 * 8.)
                    .max(10.)
                    .min(b.report.reference_length as f64);
                let height = (h.window as f64 * 8.)
                    .max(10.)
                    .min(b.report.query_length as f64);
                viewport.set(
                    Viewport {
                        x: h.reference_start as f64 - w / 2.,
                        y: h.query_start as f64 - height / 2.,
                        width: w,
                        height,
                    }
                    .clamp(&b.report),
                );
            }
        }
    });
    let navigate = Callback::new(move |query: bool| {
        if let Ok(b) = guard.run(()) {
            if let Some(h) = selected.get_untracked().and_then(|i| b.report.hits.get(i)) {
                let spans = if query {
                    &h.query_spans
                } else {
                    &h.reference_spans
                };
                on_navigate.run((
                    query,
                    query_index.get_untracked(),
                    spans.iter().map(|s| (s.start, s.end)).collect(),
                ));
            }
        }
    });
    view! {<section class="dotplot-panel" class:dotplot-expanded=move||expanded.get() role="region" aria-label="Sequence dot plot" on:keydown=move|e|if e.key()=="Escape"{expanded.set(false);}>
        <div class="dotplot-heading"><h2>"Dot plot"</h2><button class="button secondary" on:click=move|_|expanded.update(|v|*v=!*v)>{move||if expanded.get(){"Close expanded plot"}else{"Expand plot"}}</button></div>
        <p class="helper">"Compare the active reference with an open query. Each dot is the first aligned position of a retained qualifying window, using original source coordinates. This is a window comparison, not a gapped alignment or assembly verification."</p><details class="dotplot-bindings"><summary>"Source identities"</summary><p class="helper">{move||format!("Reference event: {}",source_identity.get())}</p><p class="helper">{move||info.get().map(|i|format!("Query event: {}",i.identity)).unwrap_or_else(||"No query".into())}</p></details>
        <label>"Query document"<select aria-label="Dot plot query" prop:value=move||query_index.get().to_string() on:change=move|e|query_index.set(event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
        <div class="dotplot-options"><label>"Window length"<input aria-label="Dot plot window" type="number" min="1" prop:value=move||window.get() on:input=move|e|window.set(event_target_value(&e))/></label><label>"Allowed incompatible positions"<input aria-label="Dot plot mismatches" type="number" min="0" prop:value=move||mismatches.get() on:input=move|e|mismatches.set(event_target_value(&e))/></label><label>"Window-start stride"<input aria-label="Dot plot stride" type="number" min="1" prop:value=move||stride.get() on:input=move|e|stride.set(event_target_value(&e))/></label></div>
        <label class="checkbox-label"><input aria-label="Dot plot both strands" type="checkbox" prop:checked=move||both.get()&&!protein.get() disabled=move||protein.get() on:change=move|e|both.set(event_target_checked(&e))/>"Both nucleotide strands"</label><Show when=move||protein.get()><p class="helper">"Protein windows use forward residue coordinates; reverse complement is unavailable."</p></Show><label class="checkbox-label"><input aria-label="Dot plot circular origins" type="checkbox" prop:checked=move||circular.get() on:change=move|e|circular.set(event_target_checked(&e))/>"Include origins of documents declared circular"</label>
        <details class="dotplot-budgets"><summary>"Search limits"</summary><label>"Work budget"<input aria-label="Dot plot work limit" type="number" min="1" prop:value=move||work.get() on:input=move|e|work.set(event_target_value(&e))/></label><label>"Retained windows"<input aria-label="Dot plot hit limit" type="number" min="1" prop:value=move||hits.get() on:input=move|e|hits.set(event_target_value(&e))/></label><label>"Index entries"<input aria-label="Dot plot index limit" type="number" min="1" prop:value=move||index.get() on:input=move|e|index.set(event_target_value(&e))/></label></details>
        {move||options.get().err().map(|e|view!{<p class="dialog-error">{e}</p>})}
        <button class="button primary full" disabled=move||busy.get()||options.get().is_err()||info.get().is_none() on:click=move|_|run.run(())>"Calculate dot plot"</button><Show when=move||busy.get()><button class="button secondary full" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);status.set("Dot-plot calculation cancelled.".into());}>"Cancel dot plot"</button></Show><p class="dotplot-status" role="status">{move||status.get()}</p>
        <Show when=move||data.with(|d|d.is_some())><div class="dotplot-result" data-complete=move||data.with(|d|d.as_ref().map(|b|b.report.complete).unwrap_or(false).to_string()) data-sampled=move||data.with(|d|d.as_ref().map(|b|b.report.sampled).unwrap_or(false).to_string())>
            <p>{move||data.with(|d|d.as_ref().map(|b|format!("{} retained windows · checked {} of {} scheduled windows · {}",b.report.hits.len(),b.report.examined_windows,b.report.total_windows,b.report.strategy)).unwrap_or_default())}</p><p class="helper">{move||data.with(|d|d.as_ref().map(|b|format!("Scheduled grid complete: {} · full-domain complete: {} · sampled: {}",b.report.grid_complete,b.report.complete,b.report.sampled)).unwrap_or_default())}</p>
            <div class="dotplot-controls"><button class="button secondary" on:click=move|_|if let Some(b)=data.get_untracked(){viewport.set(Viewport::all(&b.report));}>"Fit dot plot"</button><button class="button secondary" on:click=move|_|if let Some(b)=data.get_untracked(){viewport.update(|v|*v=v.zoom(&b.report,0.5,(0.5,0.5)));}>"Zoom in dot plot"</button><button class="button secondary" on:click=move|_|if let Some(b)=data.get_untracked(){viewport.update(|v|*v=v.zoom(&b.report,2.,(0.5,0.5)));}>"Zoom out dot plot"</button></div>
            <PlotCanvas data=data selected=selected viewport=viewport expanded=expanded/>
            <p class="dotplot-legend"><span>"● Forward"</span><span>"● Reverse complement"</span><span>"● Ambiguity-compatible"</span></p><p class="helper">"Drag to pan; wheel to zoom. Arrow keys pan, +/− zoom and Home fits. Coordinates below are one-based inclusive; exported spans are zero-based half-open."</p>
            {move||data.get().and_then(|b|selected.get().and_then(|i|b.report.hits.get(i).cloned().map(|h|(b,h)))).map(|(b,h)|view!{<section class="dotplot-inspection" data-hit-id=h.id data-reference-start=h.reference_start data-query-start=h.query_start data-reverse=h.reverse.to_string()><h3>"Selected window"</h3><p>{hit_label(&h)}</p><p>{format!("Reference spans {} · query spans {}",spans(&h.reference_spans),spans(&h.query_spans))}</p><pre>{format!("Reference: {}\nQuery:     {}",preview(&b.reference,h.reference_start,h.window,false),preview(&b.query,h.query_start,h.window,h.reverse))}</pre><button class="button secondary" on:click=move|_|navigate.run(false)>"Select reference window"</button><button class="button secondary" on:click=move|_|navigate.run(true)>"Open query window"</button></section>})}
            {move||data.get().map(|b|view!{<div class="dotplot-windows"><h3>"Retained windows"</h3><p class="helper">"Every retained window remains accessible here and in the exports, including overlapping dots."</p>{b.report.hits.iter().enumerate().skip(page.get()*25).take(25).map(|(i,h)|{let label=hit_label(h);view!{<button class="dotplot-hit" class:active=move||selected.get()==Some(i) data-hit-id=h.id on:click=move|_|focus.run(i)>{label}</button>}}).collect_view()}<button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous dot-plot windows"</button><button class="button secondary" disabled=move||{(page.get()+1)*25>=data.with(|d|d.as_ref().map_or(0,|b|b.report.hits.len()))} on:click=move|_|page.update(|p|*p+=1)>"Next dot-plot windows"</button>{b.report.limits.iter().chain(&b.report.warnings).take(100).map(|text|view!{<p class="helper">{text.clone()}</p>}).collect_view()}</div>})}
            <button class="button secondary full" on:click=move|_|export_json.run(())>"Download dot-plot JSON"</button><button class="button secondary full" on:click=move|_|export_svg.run(())>"Download dot-plot SVG"</button>
        </div></Show>
    </section>}
}
