use collagene_core::{Molecule, Topology, model::Feature};
use js_sys::Float32Array;
use std::{
    cell::RefCell,
    f64::consts::{PI, TAU},
    rc::Rc,
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
struct State {
    canvas: HtmlCanvasElement,
    overlay: HtmlCanvasElement,
    gl: Gl,
    program: WebGlProgram,
    buffer: web_sys::WebGlBuffer,
    vertices: i32,
    molecule: Molecule,
    revision: String,
    rotation: f64,
    zoom: f64,
    drag: Option<(f64, f64, f64, bool)>,
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
fn arc(out: &mut Vec<f32>, start: f64, end: f64, inner: f64, outer: f64, color: [f32; 4]) {
    let segments = (((end - start).abs() * 32.).ceil() as usize).clamp(1, 128);
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
    ) -> Result<Self, String> {
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
            revision: molecule.revision(),
            molecule,
            rotation: 0.,
            zoom: 1.,
            drag: None,
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
    ) {
        {
            let mut state = self.0.borrow_mut();
            if state.revision != revision {
                state.molecule = molecule.clone();
                state.revision = revision.into();
            }
        }
        self.update(selection);
    }
    pub fn zoom(&self, factor: f64) {
        let mut s = self.0.borrow_mut();
        s.zoom = (s.zoom * factor).clamp(0.6, 2.5);
        drop(s);
        self.queue_draw();
    }
    pub fn reset(&self) {
        let mut s = self.0.borrow_mut();
        s.rotation = 0.;
        s.zoom = 1.;
        drop(s);
        self.queue_draw();
    }
    pub fn pointer_down(&self, x: f64, y: f64) {
        let mut s = self.0.borrow_mut();
        let rect = s.canvas.get_bounding_client_rect();
        let angle =
            (y - rect.top() - rect.height() / 2.).atan2(x - rect.left() - rect.width() / 2.);
        s.drag = Some((angle, x, y, false));
    }
    pub fn pointer_move(&self, x: f64, y: f64) {
        let mut s = self.0.borrow_mut();
        if let Some((last, initial_x, initial_y, moved)) = s.drag {
            let rect = s.canvas.get_bounding_client_rect();
            let angle =
                (y - rect.top() - rect.height() / 2.).atan2(x - rect.left() - rect.width() / 2.);
            let moved = moved || (x - initial_x).hypot(y - initial_y) > 4.;
            if moved && s.molecule.topology == Topology::Circular {
                s.rotation += angle - last;
            }
            s.drag = Some((angle, initial_x, initial_y, moved));
            drop(s);
            self.queue_draw();
        }
    }
    pub fn pointer_up(&self, x: f64, y: f64) -> Option<(usize, usize)> {
        let mut s = self.0.borrow_mut();
        let moved = s.drag.take().map(|v| v.3).unwrap_or(true);
        if moved {
            return None;
        }
        let rect = s.canvas.get_bounding_client_rect();
        let angle =
            (y - rect.top() - rect.height() / 2.).atan2(x - rect.left() - rect.width() / 2.);
        let n = s.molecule.sequence.len();
        let pos = if s.molecule.topology == Topology::Circular {
            (((angle - s.rotation + PI / 2.).rem_euclid(TAU)) / TAU * n as f64).floor() as usize
        } else {
            let radius = rect.width().min(rect.height()) / 2. * s.zoom;
            let normalized = ((x - rect.left() - rect.width() / 2.) / radius + 0.8) / 1.6;
            if !(0.0..=1.0).contains(&normalized) {
                return None;
            }
            (normalized * n as f64).floor() as usize
        };
        s.molecule
            .features
            .iter()
            .filter(|f| f.kind != "source")
            .flat_map(|f| f.location.spans())
            .find(|(a, b, _)| pos >= *a && pos < *b)
            .map(|(a, b, _)| (a, b))
            .or(Some((pos.min(n - 1), (pos + 1).min(n))))
    }
    fn update(&self, selection: &[(usize, usize)]) {
        let mut s = self.0.borrow_mut();
        let mut vertices = vec![];
        // At most 12 MiB of static geometry. Unlike a fixed feature prefix, this
        // includes every ordinary dense-fixture feature and reports simplification.
        const MAX_FLOATS: usize = 3_000_000;
        let mut limited = false;
        let mut rendered = 0usize;
        let n = s.molecule.sequence.len() as f64;
        if s.molecule.topology == Topology::Linear {
            quad(
                &mut vertices,
                -0.8,
                0.8,
                -0.005,
                0.005,
                [0.78, 0.83, 0.81, 1.],
            );
            'linear_features: for (i, f) in s
                .molecule
                .features
                .iter()
                .filter(|f| f.kind != "source")
                .enumerate()
            {
                for (a, b, _) in f.location.spans() {
                    if vertices.len().saturating_add(4626) > MAX_FLOATS {
                        limited = true;
                        break 'linear_features;
                    }
                    let y = 0.08 + (i % 3) as f64 * 0.12;
                    quad(
                        &mut vertices,
                        a as f64 / n * 1.6 - 0.8,
                        b as f64 / n * 1.6 - 0.8,
                        y,
                        y + 0.045,
                        color(f),
                    );
                }
                rendered += 1;
            }
            for &(a, b) in selection {
                if vertices.len().saturating_add(4626) > MAX_FLOATS {
                    limited = true;
                    break;
                }
                quad(
                    &mut vertices,
                    a as f64 / n * 1.6 - 0.8,
                    b as f64 / n * 1.6 - 0.8,
                    -0.08,
                    -0.04,
                    [0.09, 0.30, 0.27, 1.],
                );
            }
        } else {
            arc(&mut vertices, 0., TAU, 0.65, 0.66, [0.78, 0.83, 0.81, 1.]);
            let n = s.molecule.sequence.len() as f64;
            'circular_features: for (i, f) in s
                .molecule
                .features
                .iter()
                .filter(|f| f.kind != "source")
                .enumerate()
            {
                let terminal =
                    if s.molecule.sequence.alphabet() == collagene_core::Alphabet::Protein {
                        None
                    } else {
                        collagene_core::tracks::terminal_span(&f.location)
                    };
                for (a, b, reverse) in f.location.spans() {
                    if vertices.len().saturating_add(4626) > MAX_FLOATS {
                        limited = true;
                        break 'circular_features;
                    }
                    let radius = 0.67 + (i % 3) as f64 * 0.057;
                    let start = a as f64 / n * TAU - PI / 2.;
                    let end = b as f64 / n * TAU - PI / 2.;
                    arc(&mut vertices, start, end, radius, radius + 0.04, color(f));
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
                                color(f)[0],
                                color(f)[1],
                                color(f)[2],
                                1.,
                            ]);
                        }
                    }
                }
                rendered += 1;
            }
            for &(a, b) in selection {
                if vertices.len().saturating_add(4626) > MAX_FLOATS {
                    limited = true;
                    break;
                }
                arc(
                    &mut vertices,
                    a as f64 / n * TAU - PI / 2.,
                    b as f64 / n * TAU - PI / 2.,
                    0.59,
                    0.63,
                    [0.09, 0.30, 0.27, 1.],
                );
            }
        }
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
                draw(&mut s);
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

fn draw(s: &mut State) {
    let rect = s.canvas.get_bounding_client_rect();
    let w = rect.width().max(1.);
    let h = rect.height().max(1.);
    let dpr = web_sys::window()
        .map(|w| w.device_pixel_ratio())
        .unwrap_or(1.)
        .min(2.);
    let pw = (w * dpr) as u32;
    let ph = (h * dpr) as u32;
    if s.canvas.width() != pw || s.canvas.height() != ph {
        s.canvas.set_width(pw);
        s.canvas.set_height(ph);
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
        (w.min(h) / w) as f32,
        (w.min(h) / h) as f32,
    );
    gl.uniform1f(
        gl.get_uniform_location(&s.program, "u_rotation").as_ref(),
        if s.molecule.topology == Topology::Circular {
            s.rotation as f32
        } else {
            0.
        },
    );
    gl.uniform1f(
        gl.get_uniform_location(&s.program, "u_zoom").as_ref(),
        s.zoom as f32,
    );
    gl.draw_arrays(Gl::TRIANGLES, 0, s.vertices);
    let Some(ctx) = s
        .overlay
        .get_context("2d")
        .ok()
        .flatten()
        .and_then(|c| c.dyn_into::<Ctx>().ok())
    else {
        return;
    };
    let _ = ctx.set_transform(dpr, 0., 0., dpr, 0., 0.);
    ctx.clear_rect(0., 0., w, h);
    if s.molecule.topology == Topology::Linear {
        let radius = w.min(h) / 2. * s.zoom;
        let n = s.molecule.sequence.len();
        ctx.set_text_align("center");
        ctx.set_font("600 19px system-ui");
        ctx.set_fill_style_str("#254b44");
        let _ = ctx.fill_text(
            &s.molecule.name.chars().take(28).collect::<String>(),
            w / 2.,
            h / 2. + 85.,
        );
        ctx.set_font("14px system-ui");
        ctx.set_fill_style_str("#72817c");
        let _ = ctx.fill_text(
            &format!("{n} {} · linear", s.molecule.sequence.alphabet().map_unit()),
            w / 2.,
            h / 2. + 112.,
        );
        ctx.set_font("12px system-ui");
        for i in 0..=4 {
            let x = w / 2. + (i as f64 / 4. * 1.6 - 0.8) * radius;
            let _ = ctx.fill_text(&format!("{}", (i * n / 4).max(1)), x, h / 2. + 28.);
        }
        for (i, f) in s
            .molecule
            .features
            .iter()
            .filter(|f| f.kind != "source")
            .take(12)
            .enumerate()
        {
            if let Some((a, b, _)) = f.location.spans().first() {
                let x = w / 2. + ((*a + *b) as f64 / 2. / n as f64 * 1.6 - 0.8) * radius;
                let y = h / 2. - (0.08 + (i % 3) as f64 * 0.12 + 0.06) * radius;
                ctx.set_font("14px system-ui");
                ctx.set_fill_style_str(&css_color(f));
                let _ = ctx.fill_text(&f.label().chars().take(20).collect::<String>(), x, y);
            }
        }
        return;
    }
    ctx.set_text_align("center");
    ctx.set_fill_style_str("#254b44");
    ctx.set_font("600 19px system-ui");
    let title: String = s.molecule.name.chars().take(28).collect();
    let _ = ctx.fill_text(&title, w / 2., h / 2. - 6.);
    ctx.set_font("14px system-ui");
    ctx.set_fill_style_str("#72817c");
    let _ = ctx.fill_text(
        &format!(
            "{} {}",
            s.molecule.sequence.len(),
            s.molecule.sequence.alphabet().short_unit()
        ),
        w / 2.,
        h / 2. + 20.,
    );
    let radius = w.min(h) / 2. * s.zoom;
    ctx.set_font("12px system-ui");
    ctx.set_fill_style_str("#87938e");
    let n = s.molecule.sequence.len();
    for i in 0..8 {
        let a = i as f64 / 8. * TAU - PI / 2. + s.rotation;
        let _ = ctx.fill_text(
            &format!("{}", i * n / 8 + 1),
            w / 2. + a.cos() * radius * 0.56,
            h / 2. + a.sin() * radius * 0.56 + 4.,
        );
    }
    ctx.set_font("14px system-ui");
    let mut occupied: Vec<(f64, f64, f64, f64)> = vec![];
    for f in s
        .molecule
        .features
        .iter()
        .filter(|f| f.kind != "source")
        .take(12)
    {
        if let Some((a, b, _)) = f.location.spans().first() {
            let angle = (*a + *b) as f64 / 2. / n as f64 * TAU - PI / 2. + s.rotation;
            let mut label: String = f.label().chars().take(28).collect();
            let measure = |text: &str| {
                ctx.measure_text(text)
                    .map(|m| m.width())
                    .unwrap_or(text.len() as f64 * 7.5)
            };
            let mut shortened = false;
            while measure(&label) > w * 0.43 && label.chars().count() > 3 {
                label.pop();
                shortened = true;
            }
            if shortened {
                label.push('…');
            }
            let text_width = measure(&label);
            let x = (w / 2. + angle.cos() * radius * 0.94).clamp(
                12. + text_width / 2.,
                (w - 12. - text_width / 2.).max(12. + text_width / 2.),
            );
            let initial_y = h / 2. + angle.sin() * radius * 0.94 + 4.;
            let placement = [0., -20., 20., -40., 40., -60., 60.]
                .into_iter()
                .map(|offset| initial_y + offset)
                .find(|y| {
                    *y >= 20.
                        && *y <= h - 14.
                        && !occupied.iter().any(|(a, b, c, d)| {
                            x + text_width / 2. + 5. > *a
                                && x - text_width / 2. - 5. < *c
                                && *y + 5. > *b
                                && *y - 16. < *d
                        })
                });
            let Some(y) = placement else {
                continue;
            };
            occupied.push((x - text_width / 2., y - 16., x + text_width / 2., y + 5.));
            ctx.set_stroke_style_str("#c8d5cd");
            ctx.set_line_width(0.8);
            ctx.begin_path();
            ctx.move_to(
                w / 2. + angle.cos() * radius * 0.82,
                h / 2. + angle.sin() * radius * 0.82,
            );
            ctx.line_to(x, y - 6.);
            ctx.stroke();
            ctx.set_fill_style_str(&css_color(f));
            let _ = ctx.fill_text(&label, x, y);
        }
    }
}
