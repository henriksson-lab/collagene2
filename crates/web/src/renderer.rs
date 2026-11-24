use crate::{
    annotation_geometry::{self, PreparedGeometry},
    map_labels::{self, LabelStats, View},
    restriction_map,
};
use collagene_core::{Molecule, Topology, model::Feature};
use js_sys::Float32Array;
use std::{
    cell::RefCell,
    f64::consts::{PI, TAU},
    rc::Rc,
    sync::Arc,
};
use wasm_bindgen::{JsCast, closure::Closure};
use web_sys::{
    CanvasRenderingContext2d as Ctx, HtmlCanvasElement, WebGl2RenderingContext as Gl, WebGlProgram,
    WebGlShader,
};

pub fn color(feature: &Feature) -> [f32; 4] {
    match feature.kind.as_str() {
        "CDS" => [0.17, 0.53, 0.46, 1.0],
        "promoter" => [0.84, 0.57, 0.23, 1.0],
        "rep_origin" => [0.40, 0.42, 0.70, 1.0],
        "primer_bind" => [0.31, 0.55, 0.75, 1.0],
        "terminator" => [0.76, 0.39, 0.43, 1.0],
        _ => [0.48, 0.64, 0.59, 1.0],
    }
}
pub fn css_color(feature: &Feature) -> String {
    let c = color(feature);
    format!(
        "rgb({},{},{})",
        (c[0] * 255.) as u8,
        (c[1] * 255.) as u8,
        (c[2] * 255.) as u8
    )
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MapHit {
    pub feature_index: Option<usize>,
    pub restriction_index: Option<usize>,
    pub span: (usize, usize),
    /// All selected reference parts; span remains the first/nearest part.
    pub spans: Vec<(usize, usize)>,
}

#[derive(Clone)]
pub struct MapRenderer(Rc<RefCell<State>>);
pub struct ResizeListener {
    observer: web_sys::ResizeObserver,
    _callback: Closure<dyn FnMut(js_sys::Array, web_sys::ResizeObserver)>,
    active: Rc<std::cell::Cell<bool>>,
}
impl ResizeListener {
    pub fn new(element: &web_sys::Element, callback: impl Fn() + 'static) -> Result<Self, String> {
        let active = Rc::new(std::cell::Cell::new(true));
        let active_for_callback = active.clone();
        let observed = element.clone();
        let previous = std::cell::Cell::new(None::<(i32, i32, u64, u64)>);
        let action = Rc::new(callback);
        let callback =
            Closure::<dyn FnMut(js_sys::Array, web_sys::ResizeObserver)>::new(move |_, _| {
                if !active_for_callback.get() {
                    return;
                }
                // Scrollbars can change the usable viewport while its border
                // rectangle stays fixed. Match the dimensions the canvas uses.
                let rect = observed.get_bounding_client_rect();
                let size = (
                    observed.client_width(),
                    observed.client_height(),
                    rect.width().to_bits(),
                    rect.height().to_bits(),
                );
                if previous.get() == Some(size) {
                    return;
                }
                previous.set(Some(size));
                // Never resize a canvas during the ResizeObserver delivery itself.
                let action = action.clone();
                let active = active_for_callback.clone();
                let frame = Closure::once_into_js(move || {
                    if active.get() {
                        action();
                    }
                });
                if let Some(window) = web_sys::window() {
                    let _ = window.request_animation_frame(frame.unchecked_ref());
                }
            });
        let observer = web_sys::ResizeObserver::new(callback.as_ref().unchecked_ref())
            .map_err(|_| "Cannot observe canvas size")?;
        observer.observe(element);
        Ok(Self {
            observer,
            _callback: callback,
            active,
        })
    }
}
impl Drop for ResizeListener {
    fn drop(&mut self) {
        self.active.set(false);
        self.observer.disconnect();
    }
}
#[derive(Clone, Copy)]
struct MapDrag {
    last_angle: f64,
    initial_x: f64,
    initial_y: f64,
    moved: bool,
    start_base: usize,
    end_base: usize,
    start_position: f64,
    sweep: f64,
}

/// A drag selects inclusive base cells, not just the final shortest arc. Keeping
/// a fractional source start distinguishes a near-full turn from a tiny move
/// within the same base, including one-base and origin-crossing references.
fn circular_selection(n: usize, start_position: f64, sweep: f64) -> Vec<(usize, usize)> {
    if n == 0 || !start_position.is_finite() || !sweep.is_finite() {
        return vec![];
    }
    let start = start_position.floor() as i64;
    let end = (start_position + sweep / TAU * n as f64).floor() as i64;
    if sweep.abs() >= TAU || start.abs_diff(end).saturating_add(1) >= n as u64 {
        return vec![(0, n)];
    }
    let first = start.rem_euclid(n as i64) as usize;
    let last = end.rem_euclid(n as i64) as usize;
    if sweep >= 0. {
        if first <= last {
            vec![(first, last + 1)]
        } else {
            vec![(first, n), (0, last + 1)]
        }
    } else if last <= first {
        vec![(last, first + 1)]
    } else {
        vec![(0, first + 1), (last, n)]
    }
}
fn shortest_angle(from: f64, to: f64) -> f64 {
    let raw = to - from;
    let delta = (raw + PI).rem_euclid(TAU) - PI;
    if delta == -PI && raw > 0. { PI } else { delta }
}
struct State {
    canvas: HtmlCanvasElement,
    overlay: HtmlCanvasElement,
    gl: Gl,
    program: WebGlProgram,
    buffer: web_sys::WebGlBuffer,
    vertices: i32,
    molecule: Molecule,
    revision: String,
    interaction_context: String,
    geometry: PreparedGeometry,
    geometry_valid: bool,
    terminals: Vec<Option<(usize, usize, bool)>>,
    annotation_total: usize,
    labels: Option<map_labels::Cache>,
    label_hits: Vec<map_labels::Hit>,
    restriction_sites: Arc<Vec<collagene_core::restriction_inspection::InspectionSite>>,
    restriction_cache: Option<restriction_map::Cache>,
    restriction_hits: Vec<restriction_map::Hit>,
    selected_restriction: Option<usize>,
    label_stats: Option<LabelStats>,
    on_labels: Rc<dyn Fn(LabelStats)>,
    overlay_ctx: Ctx,
    selected_feature: Option<usize>,
    rendered_geometry_count: usize,
    geometry_limited: bool,
    display_title: String,
    display_subtitle: String,
    hit_dataset: String,
    rotation: f64,
    zoom: f64,
    drag: Option<MapDrag>,
    selection_mode: bool,
    range_preview: Vec<(usize, usize)>,
    queued: bool,
    _resize: Option<ResizeListener>,
}

fn shader(gl: &Gl, kind: u32, source: &str) -> Result<WebGlShader, String> {
    let shader = gl.create_shader(kind).ok_or("Cannot create shader")?;
    gl.shader_source(&shader, source);
    gl.compile_shader(&shader);
    if gl
        .get_shader_parameter(&shader, Gl::COMPILE_STATUS)
        .as_bool()
        == Some(true)
    {
        Ok(shader)
    } else {
        Err(gl.get_shader_info_log(&shader).unwrap_or_default())
    }
}
fn arc_segments(start: f64, end: f64) -> usize {
    (((end - start).abs() * 32.).ceil() as usize).clamp(1, 128)
}
fn arc(out: &mut Vec<f32>, start: f64, end: f64, inner: f64, outer: f64, color: [f32; 4]) {
    let segments = arc_segments(start, end);
    let mut vertex = |angle: f64, radius: f64| {
        out.extend_from_slice(&[
            (angle.cos() * radius) as f32,
            (-angle.sin() * radius) as f32,
            color[0],
            color[1],
            color[2],
            color[3],
        ]);
    };
    for i in 0..segments {
        let a = start + (end - start) * i as f64 / segments as f64;
        let b = start + (end - start) * (i + 1) as f64 / segments as f64;
        vertex(a, inner);
        vertex(a, outer);
        vertex(b, outer);
        vertex(a, inner);
        vertex(b, outer);
        vertex(b, inner);
    }
}

fn quad(out: &mut Vec<f32>, a: f64, b: f64, c: f64, d: f64, color: [f32; 4]) {
    for (x, y) in [(a, c), (a, d), (b, d), (a, c), (b, d), (b, c)] {
        out.extend_from_slice(&[x as f32, y as f32, color[0], color[1], color[2], color[3]]);
    }
}

impl MapRenderer {
    pub fn new(
        canvas: HtmlCanvasElement,
        overlay: HtmlCanvasElement,
        molecule: Molecule,
        on_labels: Rc<dyn Fn(LabelStats)>,
    ) -> Result<Self, String> {
        let overlay_ctx = overlay
            .get_context("2d")
            .map_err(|_| "Cannot access map label canvas")?
            .ok_or("Map label canvas is unavailable")?
            .dyn_into::<Ctx>()
            .map_err(|_| "Invalid map label context")?;
        let prepared = annotation_geometry::prepare(&molecule);
        let geometry_valid = prepared.is_ok();
        let geometry = prepared.unwrap_or_default();
        let terminals = prepare_terminals(&molecule, &geometry);
        let annotation_total = molecule.features.len();
        // A failed borrowed geometry guard must not be followed by a rich-model hash.
        let revision = if geometry_valid {
            molecule.revision()
        } else {
            String::new()
        };
        let labels =
            geometry_valid.then(|| map_labels::Cache::prepare(&molecule, &geometry, &overlay_ctx));
        let display_title = molecule.name.chars().take(28).collect();
        let display_subtitle = subtitle(&molecule);
        let gl: Gl = canvas
            .get_context("webgl2")
            .map_err(|_| "Cannot access WebGL2")?
            .ok_or("WebGL2 is unavailable; sequence editing remains available")?
            .dyn_into()
            .map_err(|_| "Invalid WebGL2 context")?;
        let vertex = shader(
            &gl,
            Gl::VERTEX_SHADER,
            "#version 300 es\nin vec2 a_pos;in vec4 a_color;uniform vec2 u_aspect;uniform float u_rotation;uniform float u_zoom;out vec4 v_color;void main(){float c=cos(u_rotation),s=sin(u_rotation);vec2 p=vec2(a_pos.x*c+a_pos.y*s,a_pos.y*c-a_pos.x*s);gl_Position=vec4(p*u_aspect*u_zoom,0.,1.);v_color=a_color;}",
        )?;
        let fragment = shader(
            &gl,
            Gl::FRAGMENT_SHADER,
            "#version 300 es\nprecision mediump float;in vec4 v_color;out vec4 out_color;void main(){out_color=v_color;}",
        )?;
        let program = gl.create_program().ok_or("Cannot create WebGL program")?;
        gl.attach_shader(&program, &vertex);
        gl.attach_shader(&program, &fragment);
        gl.link_program(&program);
        gl.delete_shader(Some(&vertex));
        gl.delete_shader(Some(&fragment));
        if gl
            .get_program_parameter(&program, Gl::LINK_STATUS)
            .as_bool()
            != Some(true)
        {
            return Err(gl.get_program_info_log(&program).unwrap_or_default());
        }
        let buffer = gl.create_buffer().ok_or("Cannot create GPU buffer")?;
        let renderer = Self(Rc::new(RefCell::new(State {
            canvas,
            overlay,
            gl,
            program,
            buffer,
            vertices: 0,
            revision,
            interaction_context: String::new(),
            molecule,
            geometry,
            geometry_valid,
            terminals,
            annotation_total,
            labels,
            label_hits: Vec::new(),
            restriction_sites: Arc::new(Vec::new()),
            restriction_cache: None,
            restriction_hits: Vec::new(),
            selected_restriction: None,
            label_stats: None,
            on_labels,
            overlay_ctx,
            selected_feature: None,
            rendered_geometry_count: 0,
            geometry_limited: !geometry_valid,
            display_title,
            display_subtitle,
            hit_dataset: String::new(),
            rotation: 0.,
            zoom: 1.,
            drag: None,
            selection_mode: false,
            range_preview: Vec::new(),
            queued: false,
            _resize: None,
        })));
        let weak = Rc::downgrade(&renderer.0);
        let resize = ResizeListener::new(renderer.0.borrow().canvas.unchecked_ref(), move || {
            if let Some(state) = weak.upgrade() {
                MapRenderer(state).queue_draw();
            }
        })?;
        renderer.0.borrow_mut()._resize = Some(resize);
        renderer.update(&[]);
        Ok(renderer)
    }
    pub fn set_document_spans(
        &self,
        molecule: &Molecule,
        revision: &str,
        selection: &[(usize, usize)],
        selected_feature: Option<usize>,
    ) {
        {
            let mut state = self.0.borrow_mut();
            if state.revision != revision {
                // A gesture started on another source cannot select this document.
                state.drag = None;
                state.range_preview.clear();
                let prepared = annotation_geometry::prepare(molecule);
                state.geometry_valid = prepared.is_ok();
                state.geometry = prepared.unwrap_or_default();
                state.terminals = prepare_terminals(molecule, &state.geometry);
                state.annotation_total = molecule.features.len();
                // Reject the annotation cache without recursively cloning an invalid tree.
                // This private fallback retains only reference symbols for baseline display.
                state.molecule = if state.geometry_valid {
                    molecule.clone()
                } else {
                    Molecule {
                        name: molecule.name.chars().take(28).collect(),
                        sequence: molecule.sequence.clone(),
                        topology: molecule.topology,
                        features: vec![],
                        metadata: Default::default(),
                        original_genbank: None,
                    }
                };
                state.revision = revision.into();
                state.labels = state.geometry_valid.then(|| {
                    map_labels::Cache::prepare(molecule, &state.geometry, &state.overlay_ctx)
                });
                state.label_hits.clear();
                state.restriction_sites = Arc::new(Vec::new());
                state.restriction_cache = None;
                state.restriction_hits.clear();
                state.selected_restriction = None;
                let _ = state.canvas.remove_attribute("data-restriction-source");
                let _ = state.canvas.set_attribute("data-restriction-hits", "[]");
                state.hit_dataset.clear();
                let _ = state.canvas.remove_attribute("data-label-source");
                let _ = state.canvas.remove_attribute("title");
                state.display_title = molecule.name.chars().take(28).collect();
                state.display_subtitle = subtitle(molecule);
            }
            state.selected_feature = selected_feature;
        }
        self.update(selection);
    }

    /// Gestures belong to a document event, including equal-content documents.
    pub fn set_interaction_context(&self, context: &str) {
        let mut state = self.0.borrow_mut();
        if state.interaction_context != context {
            state.drag = None;
            state.range_preview.clear();
            state.interaction_context = context.into();
            drop(state);
            self.queue_draw();
        }
    }

    /// False preserves circular rotation; true turns a drag into a reference range.
    pub fn set_selection_mode(&self, selected: bool) {
        let mut state = self.0.borrow_mut();
        if state.selection_mode == selected {
            return;
        }
        state.selection_mode = selected;
        state.drag = None;
        state.range_preview.clear();
        drop(state);
        self.queue_draw();
    }

    pub fn pointer_cancel(&self) {
        let mut state = self.0.borrow_mut();
        let changed = state.drag.take().is_some() || !state.range_preview.is_empty();
        state.range_preview.clear();
        drop(state);
        if changed {
            self.queue_draw();
        }
    }

    /// Bind only a freshly inspected report for this source revision. The host
    /// also clears on document/event identity changes, including equal DNA.
    pub fn set_restriction_sites(
        &self,
        revision: &str,
        sites: Arc<Vec<collagene_core::restriction_inspection::InspectionSite>>,
        selected: Option<usize>,
    ) -> bool {
        let mut state = self.0.borrow_mut();
        if revision != state.revision {
            return false;
        }
        let changed = !Arc::ptr_eq(&state.restriction_sites, &sites);
        if !changed && state.selected_restriction == selected {
            return true;
        }
        if changed {
            let cache = restriction_map::Cache::prepare(
                &sites,
                state.molecule.sequence.len(),
                state.molecule.topology == Topology::Circular,
                &state.overlay_ctx,
            );
            state.restriction_sites = sites;
            state.restriction_cache = Some(cache);
            state.restriction_hits.clear();
            let _ = state.canvas.remove_attribute("data-restriction-source");
            let _ = state.canvas.set_attribute("data-restriction-hits", "[]");
            let _ = state.canvas.remove_attribute("title");
        }
        state.selected_restriction = selected;
        drop(state);
        self.queue_draw();
        true
    }
    pub fn zoom(&self, factor: f64) {
        if !factor.is_finite() || factor <= 0. {
            return;
        }
        let mut s = self.0.borrow_mut();
        s.drag = None;
        s.range_preview.clear();
        s.zoom = (s.zoom * factor).clamp(0.6, 2.5);
        drop(s);
        self.queue_draw();
    }
    pub fn reset(&self) {
        let mut s = self.0.borrow_mut();
        s.drag = None;
        s.range_preview.clear();
        s.rotation = 0.;
        s.zoom = 1.;
        drop(s);
        self.queue_draw();
    }
    pub fn pointer_down(&self, x: f64, y: f64) {
        self.pointer_cancel();
        if !x.is_finite() || !y.is_finite() {
            return;
        }
        let mut s = self.0.borrow_mut();
        let rect = s.canvas.get_bounding_client_rect();
        if rect.width() <= 0. || rect.height() <= 0. || s.molecule.sequence.is_empty() {
            return;
        }
        let local_x = x - rect.left();
        let local_y = y - rect.top();
        let Some(base) =
            reference_position(&s, local_x, local_y, rect.width(), rect.height(), true)
        else {
            return;
        };
        let angle = (local_y - rect.height() / 2.).atan2(local_x - rect.width() / 2.);
        let start_position =
            (angle - s.rotation + PI / 2.).rem_euclid(TAU) / TAU * s.molecule.sequence.len() as f64;
        s.drag = Some(MapDrag {
            last_angle: angle,
            initial_x: x,
            initial_y: y,
            moved: false,
            start_base: base,
            end_base: base,
            start_position,
            sweep: 0.,
        });
    }
    pub fn pointer_move(&self, x: f64, y: f64) {
        if !x.is_finite() || !y.is_finite() {
            self.pointer_cancel();
            return;
        }
        let mut s = self.0.borrow_mut();
        let rect = s.canvas.get_bounding_client_rect();
        let hit = s.label_hits.iter().rev().find(|hit| {
            x - rect.left() >= hit.rect.x
                && x - rect.left() < hit.rect.right()
                && y - rect.top() >= hit.rect.y
                && y - rect.top() < hit.rect.bottom()
        });
        let title = hit
            .and_then(|hit| s.molecule.features.get(hit.feature_index))
            .map(Feature::label)
            .or_else(|| {
                s.restriction_hits
                    .iter()
                    .rev()
                    .find(|hit| hit.contains(x - rect.left(), y - rect.top()))
                    .and_then(|hit| s.restriction_cache.as_ref()?.tooltip(hit.index))
            });
        if let Some(title) = title {
            if s.canvas.get_attribute("title").as_deref() != Some(title) {
                let _ = s.canvas.set_attribute("title", title);
            }
        } else {
            let _ = s.canvas.remove_attribute("title");
        }
        if s.drag.is_some() {
            advance_gesture(&mut s, x, y);
            drop(s);
            self.queue_draw();
        }
    }
    pub fn pointer_up(&self, x: f64, y: f64) -> Option<MapHit> {
        if !x.is_finite() || !y.is_finite() {
            self.pointer_cancel();
            return None;
        }
        let mut s = self.0.borrow_mut();
        // Include the release position even when no final move was dispatched.
        advance_gesture(&mut s, x, y);
        let drag = s.drag.take()?;
        let spans = std::mem::take(&mut s.range_preview);
        if drag.moved {
            let selected = s.selection_mode;
            drop(s);
            self.queue_draw();
            return if selected {
                spans.first().copied().map(|span| MapHit {
                    feature_index: None,
                    restriction_index: None,
                    span,
                    spans,
                })
            } else {
                None
            };
        }
        if s.molecule.sequence.is_empty() {
            return None;
        }
        let rect = s.canvas.get_bounding_client_rect();
        if rect.width() <= 0. || rect.height() <= 0. {
            return None;
        }
        let local_x = x - rect.left();
        let local_y = y - rect.top();
        // Measured text boxes take precedence over the source-coordinate fallback.
        if let Some(hit) = s.label_hits.iter().rev().find(|hit| {
            local_x >= hit.rect.x
                && local_x < hit.rect.right()
                && local_y >= hit.rect.y
                && local_y < hit.rect.bottom()
        }) {
            let pos = reference_position(&s, local_x, local_y, rect.width(), rect.height(), true)?;
            return nearest_span(&s, hit.feature_index, pos).map(|span| MapHit {
                feature_index: Some(hit.feature_index),
                restriction_index: None,
                span,
                spans: feature_spans(&s, hit.feature_index, span),
            });
        }
        if let Some(hit) = s
            .restriction_hits
            .iter()
            .rev()
            .find(|hit| hit.contains(local_x, local_y))
        {
            return Some(MapHit {
                feature_index: None,
                restriction_index: Some(hit.index),
                span: hit.span,
                spans: vec![hit.span],
            });
        }
        let pos = reference_position(&s, local_x, local_y, rect.width(), rect.height(), false)?;
        for located in s.geometry.located.iter().take(s.rendered_geometry_count) {
            if let Some(&(a, b, _)) = located.spans.iter().find(|&&(a, b, _)| pos >= a && pos < b) {
                return Some(MapHit {
                    feature_index: Some(located.index),
                    restriction_index: None,
                    span: (a, b),
                    spans: feature_spans(&s, located.index, (a, b)),
                });
            }
        }
        Some(MapHit {
            feature_index: None,
            restriction_index: None,
            span: (pos, pos + 1),
            spans: vec![(pos, pos + 1)],
        })
    }
    fn update(&self, selection: &[(usize, usize)]) {
        let mut s = self.0.borrow_mut();
        const MAX_FLOATS: usize = 3_000_000;
        // Reserve only the bounded numeric payload estimate, so an ordinary
        // plasmid does not allocate the full 12 MB geometry ceiling on selection.
        let mut vertices = Vec::with_capacity(geometry_capacity(&s, selection, MAX_FLOATS));
        let mut limited = !s.geometry_valid;
        let mut rendered = 0usize;
        let n = s.molecule.sequence.len();
        if n != 0 {
            let nf = n as f64;
            if s.molecule.topology == Topology::Linear {
                quad(
                    &mut vertices,
                    -0.8,
                    0.8,
                    -0.005,
                    0.005,
                    [0.78, 0.83, 0.81, 1.],
                );
            } else {
                arc(&mut vertices, 0., TAU, 0.65, 0.66, [0.78, 0.83, 0.81, 1.]);
            }
            for (ordinal, located) in s.geometry.located.iter().enumerate() {
                let feature = &s.molecule.features[located.index];
                let terminal = s.terminals[ordinal];
                // Admit the complete feature, preserving the prepared prefix used
                // by the label cache. No half-rendered joined item claims a label.
                let needed = located
                    .spans
                    .iter()
                    .try_fold(0usize, |sum, &(a, b, reverse)| {
                        let floats = if s.molecule.topology == Topology::Linear {
                            36
                        } else {
                            let start = a as f64 / nf * TAU - PI / 2.;
                            let end = b as f64 / nf * TAU - PI / 2.;
                            let width = end - start;
                            let segments = arc_segments(start, end);
                            segments * 36
                                + usize::from(terminal == Some((a, b, reverse)) && width > 0.035)
                                    * 18
                        };
                        sum.checked_add(floats)
                    });
                if needed.is_none_or(|needed| needed > MAX_FLOATS - vertices.len()) {
                    limited = true;
                    break;
                }
                let rgba = color(feature);
                for &(a, b, reverse) in &located.spans {
                    if s.molecule.topology == Topology::Linear {
                        let y = 0.08 + (ordinal % 3) as f64 * 0.12;
                        quad(
                            &mut vertices,
                            a as f64 / nf * 1.6 - 0.8,
                            b as f64 / nf * 1.6 - 0.8,
                            y,
                            y + 0.045,
                            rgba,
                        );
                    } else {
                        let radius = 0.67 + (ordinal % 3) as f64 * 0.057;
                        let start = a as f64 / nf * TAU - PI / 2.;
                        let end = b as f64 / nf * TAU - PI / 2.;
                        arc(&mut vertices, start, end, radius, radius + 0.04, rgba);
                        if terminal == Some((a, b, reverse)) && end - start > 0.035 {
                            let head = if reverse { start } else { end };
                            let base = if reverse { start + 0.04 } else { end - 0.04 };
                            for (angle, r) in [
                                (head, radius + 0.02),
                                (base, radius - 0.015),
                                (base, radius + 0.055),
                            ] {
                                vertices.extend_from_slice(&[
                                    (angle.cos() * r) as f32,
                                    (-angle.sin() * r) as f32,
                                    rgba[0],
                                    rgba[1],
                                    rgba[2],
                                    rgba[3],
                                ]);
                            }
                        }
                    }
                }
                rendered += 1;
            }
            for &(a, b) in selection {
                if a >= b || b > n {
                    continue;
                }
                let needed = if s.molecule.topology == Topology::Linear {
                    36
                } else {
                    arc_segments(a as f64 / nf * TAU - PI / 2., b as f64 / nf * TAU - PI / 2.) * 36
                };
                if needed > MAX_FLOATS - vertices.len() {
                    limited = true;
                    break;
                }
                if s.molecule.topology == Topology::Linear {
                    quad(
                        &mut vertices,
                        a as f64 / nf * 1.6 - 0.8,
                        b as f64 / nf * 1.6 - 0.8,
                        -0.08,
                        -0.04,
                        [0.09, 0.30, 0.27, 1.],
                    );
                } else {
                    arc(
                        &mut vertices,
                        a as f64 / nf * TAU - PI / 2.,
                        b as f64 / nf * TAU - PI / 2.,
                        0.59,
                        0.63,
                        [0.09, 0.30, 0.27, 1.],
                    );
                }
            }
        }
        s.rendered_geometry_count = rendered;
        s.geometry_limited = limited;
        let _ = s
            .canvas
            .set_attribute("data-rendered-feature-count", &rendered.to_string());
        let _ = s.canvas.set_attribute(
            "data-geometry-limited",
            if limited { "true" } else { "false" },
        );
        s.gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&s.buffer));
        let array = Float32Array::from(vertices.as_slice());
        s.gl.buffer_data_with_array_buffer_view(Gl::ARRAY_BUFFER, &array, Gl::STATIC_DRAW);
        s.vertices = (vertices.len() / 6) as i32;
        drop(s);
        self.queue_draw();
    }
    pub fn queue_draw(&self) {
        {
            let mut s = self.0.borrow_mut();
            if s.queued {
                return;
            }
            s.queued = true;
        }
        let weak = Rc::downgrade(&self.0);
        let callback = Closure::once_into_js(move || {
            if let Some(s) = weak.upgrade() {
                let mut s = s.borrow_mut();
                s.queued = false;
                let changed = draw(&mut s);
                let on_labels = s.on_labels.clone();
                drop(s);
                if let Some(stats) = changed {
                    on_labels(stats);
                }
            }
        });
        if let Some(window) = web_sys::window() {
            let _ = window.request_animation_frame(callback.unchecked_ref());
        }
    }
}

impl Drop for State {
    fn drop(&mut self) {
        self.gl.delete_buffer(Some(&self.buffer));
        self.gl.delete_program(Some(&self.program));
    }
}

fn prepare_terminals(
    molecule: &Molecule,
    geometry: &PreparedGeometry,
) -> Vec<Option<(usize, usize, bool)>> {
    geometry
        .located
        .iter()
        .map(|feature| {
            if molecule.sequence.alphabet() == collagene_core::Alphabet::Protein {
                None
            } else {
                collagene_core::tracks::terminal_span(&molecule.features[feature.index].location)
            }
        })
        .collect()
}

fn geometry_capacity(s: &State, selection: &[(usize, usize)], cap: usize) -> usize {
    let n = s.molecule.sequence.len();
    if n == 0 {
        return 0;
    }
    let nf = n as f64;
    let circular = s.molecule.topology == Topology::Circular;
    let mut total = if circular {
        arc_segments(0., TAU) * 36
    } else {
        36
    };
    for (ordinal, feature) in s.geometry.located.iter().enumerate() {
        for &(a, b, reverse) in &feature.spans {
            let amount = if circular {
                let start = a as f64 / nf * TAU - PI / 2.;
                let end = b as f64 / nf * TAU - PI / 2.;
                arc_segments(start, end) * 36
                    + usize::from(
                        s.terminals[ordinal] == Some((a, b, reverse)) && end - start > 0.035,
                    ) * 18
            } else {
                36
            };
            total = total.saturating_add(amount);
            if total >= cap {
                return cap;
            }
        }
    }
    for &(a, b) in selection {
        if a >= b || b > n {
            continue;
        }
        total = total.saturating_add(if circular {
            arc_segments(a as f64 / nf * TAU - PI / 2., b as f64 / nf * TAU - PI / 2.) * 36
        } else {
            36
        });
        if total >= cap {
            return cap;
        }
    }
    total
}

fn subtitle(m: &Molecule) -> String {
    match m.topology {
        Topology::Linear => format!(
            "{} {} · linear",
            m.sequence.len(),
            m.sequence.alphabet().map_unit()
        ),
        Topology::Circular => format!(
            "{} {}",
            m.sequence.len(),
            m.sequence.alphabet().short_unit()
        ),
    }
}

fn reference_position(
    s: &State,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
    clamp_linear: bool,
) -> Option<usize> {
    let n = s.molecule.sequence.len();
    if n == 0 || width <= 0. || height <= 0. {
        return None;
    }
    let normalized = if s.molecule.topology == Topology::Circular {
        ((y - height / 2.).atan2(x - width / 2.) - s.rotation + PI / 2.).rem_euclid(TAU) / TAU
    } else {
        let radius = map_labels::radius(View {
            width,
            height,
            rotation: s.rotation,
            zoom: s.zoom,
            topology: s.molecule.topology,
        });
        if radius <= 0. {
            return None;
        }
        let p = ((x - width / 2.) / radius + 0.8) / 1.6;
        if !clamp_linear && !(0.0..=1.0).contains(&p) {
            return None;
        }
        p.clamp(0., 1.)
    };
    Some(((normalized * n as f64).floor() as usize).min(n - 1))
}

fn advance_gesture(s: &mut State, x: f64, y: f64) {
    let Some(mut drag) = s.drag.take() else {
        return;
    };
    let rect = s.canvas.get_bounding_client_rect();
    if !x.is_finite()
        || !y.is_finite()
        || rect.width() <= 0.
        || rect.height() <= 0.
        || s.molecule.sequence.is_empty()
    {
        s.range_preview.clear();
        return;
    }
    let local_x = x - rect.left();
    let local_y = y - rect.top();
    let angle = (local_y - rect.height() / 2.).atan2(local_x - rect.width() / 2.);
    let delta = shortest_angle(drag.last_angle, angle);
    drag.last_angle = angle;
    drag.moved |= (x - drag.initial_x).hypot(y - drag.initial_y) > 4.;
    if s.selection_mode {
        drag.sweep = (drag.sweep + delta).clamp(-4. * TAU, 4. * TAU);
        if let Some(base) =
            reference_position(s, local_x, local_y, rect.width(), rect.height(), true)
        {
            drag.end_base = base;
        }
        if drag.moved {
            s.range_preview = if s.molecule.topology == Topology::Circular {
                circular_selection(s.molecule.sequence.len(), drag.start_position, drag.sweep)
            } else {
                vec![(
                    drag.start_base.min(drag.end_base),
                    drag.start_base.max(drag.end_base) + 1,
                )]
            };
        }
    } else if drag.moved && s.molecule.topology == Topology::Circular {
        s.rotation = (s.rotation + delta).rem_euclid(TAU);
    }
    s.drag = Some(drag);
}
fn feature_spans(s: &State, index: usize, nearest: (usize, usize)) -> Vec<(usize, usize)> {
    let Some(located) = s
        .geometry
        .located
        .iter()
        .find(|located| located.index == index)
    else {
        return vec![nearest];
    };
    let mut spans = located
        .spans
        .iter()
        .map(|&(a, b, _)| (a, b))
        .collect::<Vec<_>>();
    if let Some(first) = spans.iter().position(|span| *span == nearest) {
        spans.rotate_left(first);
    }
    if spans.is_empty() {
        spans.push(nearest);
    }
    spans
}

fn nearest_span(s: &State, index: usize, pos: usize) -> Option<(usize, usize)> {
    let located = s.geometry.located.get(
        s.geometry
            .located
            .binary_search_by_key(&index, |f| f.index)
            .ok()?,
    )?;
    let n = s.molecule.sequence.len();
    located
        .spans
        .iter()
        .min_by_key(|&&(a, b, _)| {
            if pos >= a && pos < b {
                return 0;
            }
            let direct = if pos < a { a - pos } else { pos - (b - 1) };
            if s.molecule.topology == Topology::Circular {
                direct.min((pos + n - a) % n).min((b - 1 + n - pos) % n)
            } else {
                direct
            }
        })
        .map(|&(a, b, _)| (a, b))
}

fn draw(s: &mut State) -> Option<LabelStats> {
    let rect = s.canvas.get_bounding_client_rect();
    let w = rect.width().max(0.);
    let h = rect.height().max(0.);
    if w <= 0. || h <= 0. || s.molecule.sequence.is_empty() {
        s.drag = None;
        s.range_preview.clear();
    }
    let view = View {
        width: w,
        height: h,
        rotation: s.rotation,
        zoom: s.zoom,
        topology: s.molecule.topology,
    };
    let radius = map_labels::radius(view);
    let dpr = web_sys::window()
        .map(|w| w.device_pixel_ratio())
        .filter(|v| v.is_finite() && *v > 0.)
        .unwrap_or(1.)
        .min(2.);
    let pw = (w * dpr).ceil().max(1.) as u32;
    let ph = (h * dpr).ceil().max(1.) as u32;
    if s.canvas.width() != pw || s.canvas.height() != ph {
        s.canvas.set_width(pw);
        s.canvas.set_height(ph);
    }
    if s.overlay.width() != pw || s.overlay.height() != ph {
        s.overlay.set_width(pw);
        s.overlay.set_height(ph);
    }
    let gl = &s.gl;
    gl.viewport(0, 0, pw as i32, ph as i32);
    gl.clear_color(0.985, 0.990, 0.987, 1.);
    gl.clear(Gl::COLOR_BUFFER_BIT);
    gl.use_program(Some(&s.program));
    gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&s.buffer));
    let pos = gl.get_attrib_location(&s.program, "a_pos") as u32;
    let color = gl.get_attrib_location(&s.program, "a_color") as u32;
    gl.enable_vertex_attrib_array(pos);
    gl.vertex_attrib_pointer_with_i32(pos, 2, Gl::FLOAT, false, 24, 0);
    gl.enable_vertex_attrib_array(color);
    gl.vertex_attrib_pointer_with_i32(color, 4, Gl::FLOAT, false, 24, 8);
    gl.uniform2f(
        gl.get_uniform_location(&s.program, "u_aspect").as_ref(),
        if w > 0. { (radius * 2. / w) as f32 } else { 0. },
        if h > 0. { (radius * 2. / h) as f32 } else { 0. },
    );
    gl.uniform1f(
        gl.get_uniform_location(&s.program, "u_rotation").as_ref(),
        if s.molecule.topology == Topology::Circular {
            s.rotation as f32
        } else {
            0.
        },
    );
    gl.uniform1f(gl.get_uniform_location(&s.program, "u_zoom").as_ref(), 1.);
    gl.draw_arrays(Gl::TRIANGLES, 0, s.vertices);
    let ctx = &s.overlay_ctx;
    let _ = ctx.set_transform(
        if w > 0. { pw as f64 / w } else { 1. },
        0.,
        0.,
        if h > 0. { ph as f64 / h } else { 1. },
        0.,
        0.,
    );
    ctx.clear_rect(0., 0., w, h);
    let _ = s.canvas.set_attribute(
        "data-map-selection-mode",
        if s.selection_mode { "select" } else { "rotate" },
    );
    let _ = s.canvas.set_attribute(
        "data-range-preview",
        &serde_json::to_string(&s.range_preview).unwrap_or_else(|_| "[]".into()),
    );
    if w > 0. && h > 0. && !s.range_preview.is_empty() && !s.molecule.sequence.is_empty() {
        ctx.save();
        ctx.set_fill_style_str("rgba(25,128,104,0.25)");
        ctx.set_stroke_style_str("rgba(25,128,104,0.4)");
        let n = s.molecule.sequence.len() as f64;
        for &(a, b) in &s.range_preview {
            if a >= b || b > s.molecule.sequence.len() {
                continue;
            }
            if s.molecule.topology == Topology::Linear {
                let x = w / 2. + (a as f64 / n * 1.6 - 0.8) * radius;
                ctx.fill_rect(x, h / 2. - 9., (b - a) as f64 / n * 1.6 * radius, 18.);
            } else {
                ctx.set_line_width((radius * 0.08).max(4.));
                ctx.begin_path();
                let _ = ctx.arc(
                    w / 2.,
                    h / 2.,
                    radius * 0.615,
                    a as f64 / n * TAU - PI / 2. + s.rotation,
                    b as f64 / n * TAU - PI / 2. + s.rotation,
                );
                ctx.stroke();
            }
        }
        ctx.restore();
    }
    if w > 0. && h > 0. {
        ctx.set_text_align("center");
        ctx.set_font("600 19px system-ui");
        ctx.set_fill_style_str("#254b44");
        let linear = s.molecule.topology == Topology::Linear;
        let _ = ctx.fill_text(
            &s.display_title,
            w / 2.,
            h / 2. + if linear { 85. } else { -6. },
        );
        ctx.set_font("14px system-ui");
        ctx.set_fill_style_str("#72817c");
        let _ = ctx.fill_text(
            &s.display_subtitle,
            w / 2.,
            h / 2. + if linear { 112. } else { 20. },
        );
        ctx.set_font("12px system-ui");
        ctx.set_fill_style_str("#87938e");
        let n = s.molecule.sequence.len();
        if n > 0 {
            if linear {
                for i in 0..=4 {
                    let x = w / 2. + (i as f64 / 4. * 1.6 - 0.8) * radius;
                    let _ = ctx.fill_text(&(i * n / 4).max(1).to_string(), x, h / 2. + 28.);
                }
            } else {
                for i in 0..8 {
                    let a = i as f64 / 8. * TAU - PI / 2. + s.rotation;
                    let _ = ctx.fill_text(
                        &(i * n / 8 + 1).to_string(),
                        w / 2. + a.cos() * radius * 0.56,
                        h / 2. + a.sin() * radius * 0.56 + 4.,
                    );
                }
            }
        }
    }
    let restrictions = s
        .restriction_cache
        .as_mut()
        .map(|cache| cache.draw(ctx, view, s.selected_restriction));
    let (total, markers, labels, limited) = if let Some(frame) = restrictions {
        s.restriction_hits = frame.hits;
        (frame.total, frame.markers, frame.labels, frame.limited)
    } else {
        s.restriction_hits.clear();
        (0, 0, 0, false)
    };
    for (key, value) in [
        ("total", total),
        ("displayed", markers),
        ("omitted", total.saturating_sub(markers)),
        ("labels-displayed", labels),
        ("labels-omitted", total.saturating_sub(labels)),
    ] {
        let value = value.to_string();
        let name = format!("data-restriction-{key}");
        if s.canvas.get_attribute(&name).as_deref() != Some(value.as_str()) {
            let _ = s.canvas.set_attribute(&name, &value);
        }
    }
    let value = if limited { "true" } else { "false" };
    if s.canvas
        .get_attribute("data-restriction-limited")
        .as_deref()
        != Some(value)
    {
        let _ = s.canvas.set_attribute("data-restriction-limited", value);
    }
    let dataset=serde_json::to_string(&s.restriction_hits.iter().filter(|hit|hit.label).map(|hit|serde_json::json!({
        "index":hit.index,"x":hit.rect.x,"y":hit.rect.y,"width":hit.rect.width,"height":hit.rect.height,
    })).collect::<Vec<_>>()).unwrap_or_else(|_|"[]".into());
    if s.canvas.get_attribute("data-restriction-hits").as_deref() != Some(dataset.as_str()) {
        let _ = s.canvas.set_attribute("data-restriction-hits", &dataset);
    }
    let _ = s
        .canvas
        .set_attribute("data-restriction-source", &s.revision);
    let frame = s
        .labels
        .as_mut()
        .map(|cache| cache.draw(ctx, &view, s.selected_feature, s.rendered_geometry_count));
    let (hits, mut stats) = if let Some(frame) = frame {
        (frame.hits, frame.stats)
    } else {
        (
            Vec::new(),
            LabelStats {
                total: s.annotation_total,
                omitted: s.annotation_total,
                limited: true,
                ..LabelStats::default()
            },
        )
    };
    stats.limited |= s.geometry_limited;
    let hits_changed = s.label_hits.len() != hits.len()
        || s.label_hits
            .iter()
            .zip(&hits)
            .any(|(a, b)| a.feature_index != b.feature_index || a.rect != b.rect);
    s.label_hits = hits;
    // Numeric QA hit rectangles are bounded by the 128-label display budget.
    // Serialize only when displayed rectangles/identities actually changed.
    if hits_changed || s.hit_dataset.is_empty() {
        let dataset = serde_json::to_string(&s.label_hits.iter().map(|hit|serde_json::json!({
            "index":hit.feature_index,"x":hit.rect.x,"y":hit.rect.y,"width":hit.rect.width,"height":hit.rect.height,
        })).collect::<Vec<_>>()).unwrap_or_else(|_|"[]".into());
        let _ = s.canvas.set_attribute("data-label-hits", &dataset);
        s.hit_dataset = dataset;
    }
    // Publish the source binding only after the matching frame/hits are ready.
    // A new source with zero labels must replace the previous source's hit data.
    if s.canvas.get_attribute("data-label-source").as_deref() != Some(s.revision.as_str()) {
        let _ = s.canvas.set_attribute("data-label-source", &s.revision);
    }
    if s.label_stats.as_ref() != Some(&stats) {
        for (key, value) in [
            ("total", stats.total),
            ("displayed", stats.displayed),
            ("omitted", stats.omitted),
            ("abbreviated", stats.abbreviated),
            ("unplaced", stats.unplaced),
        ] {
            let _ = s
                .canvas
                .set_attribute(&format!("data-label-{key}"), &value.to_string());
        }
        let _ = s.canvas.set_attribute(
            "data-label-limited",
            if stats.limited { "true" } else { "false" },
        );
        s.label_stats = Some(stats.clone());
        Some(stats)
    } else {
        None
    }
}
