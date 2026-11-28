//! Bounded raster text atlas on a dedicated transparent WebGL2 surface.
//! Positions are Canvas2D alphabetic baselines in CSS pixels. This module does
//! no sequence interpretation; unsupported input requests a complete 2D fallback.
use js_sys::Float32Array;
use wasm_bindgen::JsCast;
use web_sys::{
    CanvasRenderingContext2d as Ctx, HtmlCanvasElement, WebGl2RenderingContext as Gl, WebGlBuffer,
    WebGlProgram, WebGlShader, WebGlTexture, WebGlUniformLocation,
};

pub const MAX_GLYPHS: usize = 65_536;
const ATLAS_SIDE: u32 = 2048;
const MAX_SURFACE_SIDE: f64 = 8192.;
const MAX_SURFACE_PIXELS: f64 = 16_777_216.;
const FIRST: u32 = 32;
const LAST: u32 = 126;
const FLOATS_PER_GLYPH: usize = 12;

#[derive(Clone, Copy, Debug)]
pub struct Glyph {
    pub character: char,
    pub x: f64,
    pub y: f64,
    pub color: [f32; 4],
}
#[derive(Clone, Copy, Default)]
struct Cell {
    offset_x: f64,
    offset_y: f64,
    width: u32,
    height: u32,
    atlas_x: u32,
    atlas_y: u32,
}
struct Atlas {
    cells: [Cell; (LAST - FIRST + 1) as usize],
    width: u32,
    height: u32,
    dpr: f64,
}

pub struct GlyphRenderer {
    canvas: HtmlCanvasElement,
    gl: Gl,
    program: WebGlProgram,
    buffer: WebGlBuffer,
    texture: WebGlTexture,
    viewport: WebGlUniformLocation,
    sampler: WebGlUniformLocation,
    font: String,
    atlas: Option<Atlas>,
    vertices: Vec<f32>,
}
const VERTEX: &str = r#"#version 300 es
precision highp float;
layout(location=0) in vec4 rect;
layout(location=1) in vec4 uv_rect;
layout(location=2) in vec4 color;
uniform vec2 viewport;
out vec2 uv;
out vec4 ink;
void main() {
    const vec2 corners[6] = vec2[6](vec2(0,0),vec2(1,0),vec2(0,1),vec2(0,1),vec2(1,0),vec2(1,1));
    vec2 p = corners[gl_VertexID];
    vec2 xy = rect.xy + p * rect.zw;
    gl_Position = vec4(xy.x / viewport.x * 2.0 - 1.0, 1.0 - xy.y / viewport.y * 2.0, 0, 1);
    uv = uv_rect.xy + p * uv_rect.zw;
    ink = color;
}"#;
const FRAGMENT: &str = r#"#version 300 es
precision highp float;
uniform sampler2D atlas;
in vec2 uv;
in vec4 ink;
out vec4 output_color;
void main() { float alpha = ink.a * texture(atlas, uv).a; output_color = vec4(ink.rgb * alpha, alpha); }
"#;
fn shader(gl: &Gl, kind: u32, source: &str) -> Result<WebGlShader, String> {
    let shader = gl
        .create_shader(kind)
        .ok_or("Cannot allocate glyph shader")?;
    gl.shader_source(&shader, source);
    gl.compile_shader(&shader);
    if gl
        .get_shader_parameter(&shader, Gl::COMPILE_STATUS)
        .as_bool()
        == Some(true)
    {
        Ok(shader)
    } else {
        let error = gl
            .get_shader_info_log(&shader)
            .unwrap_or_else(|| "Glyph shader compilation failed".into());
        gl.delete_shader(Some(&shader));
        Err(error)
    }
}
fn program(gl: &Gl) -> Result<WebGlProgram, String> {
    let vertex = shader(gl, Gl::VERTEX_SHADER, VERTEX)?;
    let fragment = match shader(gl, Gl::FRAGMENT_SHADER, FRAGMENT) {
        Ok(value) => value,
        Err(error) => {
            gl.delete_shader(Some(&vertex));
            return Err(error);
        }
    };
    let Some(program) = gl.create_program() else {
        gl.delete_shader(Some(&vertex));
        gl.delete_shader(Some(&fragment));
        return Err("Cannot allocate glyph program".into());
    };
    gl.attach_shader(&program, &vertex);
    gl.attach_shader(&program, &fragment);
    gl.link_program(&program);
    let linked = gl
        .get_program_parameter(&program, Gl::LINK_STATUS)
        .as_bool()
        == Some(true);
    gl.detach_shader(&program, &vertex);
    gl.detach_shader(&program, &fragment);
    gl.delete_shader(Some(&vertex));
    gl.delete_shader(Some(&fragment));
    if linked {
        Ok(program)
    } else {
        let error = gl
            .get_program_info_log(&program)
            .unwrap_or_else(|| "Glyph shader linking failed".into());
        gl.delete_program(Some(&program));
        Err(error)
    }
}
fn canvas() -> Result<HtmlCanvasElement, String> {
    web_sys::window()
        .and_then(|w| w.document())
        .ok_or("No canvas document")?
        .create_element("canvas")
        .map_err(|_| "Cannot allocate glyph atlas")?
        .dyn_into()
        .map_err(|_| "Cannot create glyph atlas canvas".into())
}
fn context(canvas: &HtmlCanvasElement) -> Result<Ctx, String> {
    canvas
        .get_context("2d")
        .map_err(|_| "Cannot initialize glyph atlas")?
        .ok_or("No glyph atlas context")?
        .dyn_into()
        .map_err(|_| "Invalid glyph atlas context".into())
}
impl GlyphRenderer {
    /// The supplied canvas is exclusively WebGL; it never receives a 2D context.
    /// A lost context requires constructing a new renderer after restoration.
    pub fn new(canvas: &HtmlCanvasElement, font: &str) -> Result<Self, String> {
        if font.is_empty() || font.len() > 256 {
            return Err("Glyph font exceeds bounds".into());
        }
        canvas
            .style()
            .set_property("pointer-events", "none")
            .map_err(|_| "Cannot style glyph canvas")?;
        canvas
            .style()
            .set_property("background", "transparent")
            .map_err(|_| "Cannot style glyph canvas")?;
        let options = js_sys::Object::new();
        for (key, value) in [
            ("alpha", true),
            ("premultipliedAlpha", true),
            ("antialias", false),
            ("depth", false),
            ("stencil", false),
        ] {
            js_sys::Reflect::set(&options, &key.into(), &value.into())
                .map_err(|_| "Cannot configure glyph context")?;
        }
        let gl: Gl = canvas
            .get_context_with_context_options("webgl2", &options)
            .map_err(|_| "Cannot initialize glyph WebGL2")?
            .ok_or("Glyph WebGL2 unavailable")?
            .dyn_into()
            .map_err(|_| "Invalid glyph WebGL2 context")?;
        if gl.is_context_lost() {
            return Err("Glyph WebGL2 context lost".into());
        }
        let program = program(&gl)?;
        let Some(buffer) = gl.create_buffer() else {
            gl.delete_program(Some(&program));
            return Err("Cannot allocate glyph buffer".into());
        };
        let Some(texture) = gl.create_texture() else {
            gl.delete_buffer(Some(&buffer));
            gl.delete_program(Some(&program));
            return Err("Cannot allocate glyph atlas texture".into());
        };
        let locations = gl
            .get_uniform_location(&program, "viewport")
            .zip(gl.get_uniform_location(&program, "atlas"));
        let Some((viewport, sampler)) = locations else {
            gl.delete_texture(Some(&texture));
            gl.delete_buffer(Some(&buffer));
            gl.delete_program(Some(&program));
            return Err("Glyph shader uniforms unavailable".into());
        };
        Ok(Self {
            canvas: canvas.clone(),
            gl,
            program,
            buffer,
            texture,
            viewport,
            sampler,
            font: font.into(),
            atlas: None,
            vertices: Vec::new(),
        })
    }
    fn prepare_atlas(&mut self, dpr: f64) -> Result<(), String> {
        if self.atlas.as_ref().is_some_and(|a| a.dpr == dpr) {
            return Ok(());
        }
        let canvas = canvas()?;
        canvas.set_width(1);
        canvas.set_height(1);
        let ctx = context(&canvas)?;
        ctx.set_font("1px serif");
        ctx.set_font(&self.font);
        if ctx.font() == "1px serif" && self.font != "1px serif" {
            return Err("Invalid glyph font".into());
        }
        ctx.set_text_align("left");
        ctx.set_text_baseline("alphabetic");
        let mut cells = [Cell::default(); (LAST - FIRST + 1) as usize];
        let (mut x, mut y, mut row_height, mut used_width) = (0u32, 0u32, 0u32, 0u32);
        for (index, cell) in cells.iter_mut().enumerate() {
            let text = char::from_u32(FIRST + index as u32).unwrap().to_string();
            let m = ctx
                .measure_text(&text)
                .map_err(|_| "Cannot measure glyph")?;
            let (left, right, top, bottom) = (
                -m.actual_bounding_box_left(),
                m.actual_bounding_box_right(),
                -m.actual_bounding_box_ascent(),
                m.actual_bounding_box_descent(),
            );
            if [left, right, top, bottom].iter().any(|v| !v.is_finite()) {
                return Err("Nonfinite glyph metrics".into());
            }
            let left = (left * dpr).floor() - 2.;
            let top = (top * dpr).floor() - 2.;
            let width = ((right * dpr).ceil() + 2. - left).max(1.);
            let height = ((bottom * dpr).ceil() + 2. - top).max(1.);
            if width > ATLAS_SIDE as f64
                || height > ATLAS_SIDE as f64
                || left.abs() > ATLAS_SIDE as f64
                || top.abs() > ATLAS_SIDE as f64
            {
                return Err("Glyph metrics exceed atlas bounds".into());
            }
            let (width, height) = (width as u32, height as u32);
            if x + width > ATLAS_SIDE {
                x = 0;
                y += row_height;
                row_height = 0;
            }
            if y + height > ATLAS_SIDE {
                return Err("Glyph atlas exceeds 2048 square bound".into());
            }
            *cell = Cell {
                offset_x: left / dpr,
                offset_y: top / dpr,
                width,
                height,
                atlas_x: x,
                atlas_y: y,
            };
            x += width;
            used_width = used_width.max(x);
            row_height = row_height.max(height);
        }
        let height = (y + row_height).max(1);
        let width = used_width.max(1);
        canvas.set_width(width);
        canvas.set_height(height);
        // Resizing resets all Canvas2D state; restore the exact baseline/font.
        ctx.set_transform(dpr, 0., 0., dpr, 0., 0.)
            .map_err(|_| "Cannot scale glyph atlas")?;
        ctx.set_font(&self.font);
        ctx.set_text_align("left");
        ctx.set_text_baseline("alphabetic");
        ctx.set_fill_style_str("#fff");
        for (index, cell) in cells.iter().enumerate() {
            let text = char::from_u32(FIRST + index as u32).unwrap().to_string();
            ctx.fill_text(
                &text,
                cell.atlas_x as f64 / dpr - cell.offset_x,
                cell.atlas_y as f64 / dpr - cell.offset_y,
            )
            .map_err(|_| "Cannot rasterize glyph atlas")?;
        }
        let gl = &self.gl;
        gl.active_texture(Gl::TEXTURE0);
        gl.bind_texture(Gl::TEXTURE_2D, Some(&self.texture));
        gl.pixel_storei(Gl::UNPACK_FLIP_Y_WEBGL, 0);
        gl.pixel_storei(Gl::UNPACK_PREMULTIPLY_ALPHA_WEBGL, 0);
        gl.tex_parameteri(Gl::TEXTURE_2D, Gl::TEXTURE_MIN_FILTER, Gl::LINEAR as i32);
        gl.tex_parameteri(Gl::TEXTURE_2D, Gl::TEXTURE_MAG_FILTER, Gl::LINEAR as i32);
        gl.tex_parameteri(Gl::TEXTURE_2D, Gl::TEXTURE_WRAP_S, Gl::CLAMP_TO_EDGE as i32);
        gl.tex_parameteri(Gl::TEXTURE_2D, Gl::TEXTURE_WRAP_T, Gl::CLAMP_TO_EDGE as i32);
        gl.tex_image_2d_with_u32_and_u32_and_html_canvas_element(
            Gl::TEXTURE_2D,
            0,
            Gl::RGBA8 as i32,
            Gl::RGBA,
            Gl::UNSIGNED_BYTE,
            &canvas,
        )
        .map_err(|_| "Cannot upload glyph atlas")?;
        if gl.is_context_lost() || gl.get_error() != Gl::NO_ERROR {
            return Err("Glyph atlas GPU upload failed".into());
        }
        self.atlas = Some(Atlas {
            cells,
            width,
            height,
            dpr,
        });
        Ok(())
    }
    /// All validation precedes frame mutation. An error requires hiding this
    /// surface and drawing the complete text frame with the caller's 2D path.
    pub fn paint(
        &mut self,
        width: f64,
        height: f64,
        dpr: f64,
        glyphs: &[Glyph],
    ) -> Result<(), String> {
        let result = self.paint_frame(width, height, dpr, glyphs);
        if result.is_err() {
            self.clear();
        }
        result
    }
    fn paint_frame(
        &mut self,
        width: f64,
        height: f64,
        dpr: f64,
        glyphs: &[Glyph],
    ) -> Result<(), String> {
        if self.gl.is_context_lost() {
            return Err("Glyph WebGL2 context lost".into());
        }
        if !width.is_finite()
            || !height.is_finite()
            || width < 0.
            || height < 0.
            || !dpr.is_finite()
            || !(0.25..=4.).contains(&dpr)
            || glyphs.len() > MAX_GLYPHS
        {
            return Err("Glyph frame exceeds numeric/instance bounds".into());
        }
        let pw = (width * dpr).ceil();
        let ph = (height * dpr).ceil();
        if pw > MAX_SURFACE_SIDE || ph > MAX_SURFACE_SIDE || pw * ph > MAX_SURFACE_PIXELS {
            return Err("Glyph surface exceeds bounded backing size".into());
        }
        for glyph in glyphs {
            if !(FIRST..=LAST).contains(&(glyph.character as u32))
                || !glyph.x.is_finite()
                || !glyph.y.is_finite()
                || glyph.x.abs() > 1_000_000.
                || glyph.y.abs() > 1_000_000.
                || glyph
                    .color
                    .iter()
                    .any(|c| !c.is_finite() || !(0.0..=1.0).contains(c))
            {
                return Err("Unsupported glyph or invalid glyph geometry/color".into());
            }
        }
        if width == 0. || height == 0. {
            self.clear();
            return Ok(());
        }
        self.prepare_atlas(dpr)?;
        let atlas = self.atlas.as_ref().ok_or("No glyph atlas")?;
        self.vertices.clear();
        let count = glyphs.len() * FLOATS_PER_GLYPH;
        self.vertices
            .try_reserve_exact(count)
            .map_err(|_| "Cannot allocate bounded glyph instances")?;
        for glyph in glyphs {
            let cell = atlas.cells[(glyph.character as u32 - FIRST) as usize];
            self.vertices.extend_from_slice(&[
                (glyph.x + cell.offset_x) as f32,
                (glyph.y + cell.offset_y) as f32,
                (cell.width as f64 / dpr) as f32,
                (cell.height as f64 / dpr) as f32,
                cell.atlas_x as f32 / atlas.width as f32,
                cell.atlas_y as f32 / atlas.height as f32,
                cell.width as f32 / atlas.width as f32,
                cell.height as f32 / atlas.height as f32,
                glyph.color[0],
                glyph.color[1],
                glyph.color[2],
                glyph.color[3],
            ]);
        }
        if self.canvas.width() != pw as u32 {
            self.canvas.set_width(pw as u32);
        }
        if self.canvas.height() != ph as u32 {
            self.canvas.set_height(ph as u32);
        }
        let gl = &self.gl;
        gl.viewport(0, 0, pw as i32, ph as i32);
        gl.disable(Gl::DEPTH_TEST);
        gl.disable(Gl::SCISSOR_TEST);
        gl.disable(Gl::CULL_FACE);
        gl.clear_color(0., 0., 0., 0.);
        gl.clear(Gl::COLOR_BUFFER_BIT);
        gl.enable(Gl::BLEND);
        gl.blend_func_separate(
            Gl::ONE,
            Gl::ONE_MINUS_SRC_ALPHA,
            Gl::ONE,
            Gl::ONE_MINUS_SRC_ALPHA,
        );
        gl.use_program(Some(&self.program));
        gl.uniform2f(Some(&self.viewport), width as f32, height as f32);
        gl.uniform1i(Some(&self.sampler), 0);
        gl.active_texture(Gl::TEXTURE0);
        gl.bind_texture(Gl::TEXTURE_2D, Some(&self.texture));
        gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&self.buffer));
        let array = Float32Array::from(self.vertices.as_slice());
        gl.buffer_data_with_array_buffer_view(Gl::ARRAY_BUFFER, &array, Gl::DYNAMIC_DRAW);
        for index in 0..3 {
            gl.enable_vertex_attrib_array(index);
            gl.vertex_attrib_pointer_with_i32(index, 4, Gl::FLOAT, false, 48, (index * 16) as i32);
            gl.vertex_attrib_divisor(index, 1);
        }
        gl.draw_arrays_instanced(Gl::TRIANGLES, 0, 6, glyphs.len() as i32);
        if gl.is_context_lost() || gl.get_error() != Gl::NO_ERROR {
            return Err("Glyph GPU frame failed".into());
        }
        Ok(())
    }
    pub fn clear(&mut self) {
        if !self.gl.is_context_lost() {
            self.gl.disable(Gl::SCISSOR_TEST);
            self.gl.clear_color(0., 0., 0., 0.);
            self.gl.clear(Gl::COLOR_BUFFER_BIT);
        }
    }
}
impl Drop for GlyphRenderer {
    fn drop(&mut self) {
        self.gl.delete_texture(Some(&self.texture));
        self.gl.delete_buffer(Some(&self.buffer));
        self.gl.delete_program(Some(&self.program));
    }
}
