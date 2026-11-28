use crate::biology::structure::structure_camera::Camera;
use collagene_core::protein::structure::coordinates::scene::Scene;
use wasm_bindgen::JsCast;
use web_sys::{HtmlCanvasElement, WebGl2RenderingContext as Gl, WebGlBuffer, WebGlProgram};
pub struct Renderer {
    gl: Gl,
    program: WebGlProgram,
    buffer: WebGlBuffer,
    points: i32,
    lines: i32,
}
impl Renderer {
    pub fn new(canvas: &HtmlCanvasElement, scene: &Scene) -> Result<Self, String> {
        if scene.points.len() > collagene_core::protein::structure::coordinates::MAX_ATOMS
            || scene.edges.len() > collagene_core::protein::structure::coordinates::MAX_CONNECTIONS
        {
            return Err("Structure scene exceeds its geometry bounds".into());
        }
        let mut vertices = Vec::with_capacity((scene.points.len() + scene.edges.len() * 2) * 7);
        let mut positions = std::collections::BTreeMap::new();
        for p in &scene.points {
            vertices.extend_from_slice(&p.position);
            vertices.extend_from_slice(&p.color);
            vertices.push(p.atom_index as f32);
            positions.insert(p.atom_index, p.position);
        }
        for e in &scene.edges {
            for i in e.atoms {
                vertices.extend_from_slice(
                    positions
                        .get(&i)
                        .ok_or("Connection endpoint missing from scene")?,
                );
                vertices.extend_from_slice(&[0.48, 0.53, 0.59]);
                vertices.push(-1.0);
            }
        }
        let gl = canvas
            .get_context("webgl2")
            .map_err(|_| "WebGL2 unavailable")?
            .ok_or("WebGL2 unavailable")?
            .dyn_into::<Gl>()
            .map_err(|_| "Invalid WebGL2 context")?;
        let vert = crate::viewer::map::renderer::shader(
            &gl,
            Gl::VERTEX_SHADER,
            r#"#version 300 es
layout(location=0) in vec3 position;
layout(location=1) in vec3 color;
layout(location=2) in float atom_id;
uniform float selected;
uniform vec2 angles;
uniform vec2 pan;
uniform float scale;
uniform float aspect;
uniform float point_size;
out vec3 tint;
void main(){float cy=cos(angles.x),sy=sin(angles.x),cp=cos(angles.y),sp=sin(angles.y);vec3 r=vec3(cy*position.x+sy*position.z,position.y,-sy*position.x+cy*position.z);vec3 p=vec3(r.x,cp*r.y-sp*r.z,sp*r.y+cp*r.z);gl_Position=vec4(p.x*scale/aspect+pan.x,p.y*scale+pan.y,-p.z*0.45,1.0);gl_PointSize=point_size;tint=(selected>=0.0 && abs(atom_id-selected)<0.1)?vec3(1.0,0.55,0.10):color;}
"#,
        )?;
        let frag = match crate::viewer::map::renderer::shader(
            &gl,
            Gl::FRAGMENT_SHADER,
            r#"#version 300 es
precision highp float;
in vec3 tint;
uniform bool atoms;
out vec4 output_color;
void main(){float shade=1.0;if(atoms){vec2 p=gl_PointCoord*2.0-1.0;float d=dot(p,p);if(d>1.0)discard;shade=0.45+0.55*sqrt(1.0-d);}output_color=vec4(tint*shade,1.0);}
"#,
        ) {
            Ok(s) => s,
            Err(e) => {
                gl.delete_shader(Some(&vert));
                return Err(e);
            }
        };
        let program = match gl.create_program() {
            Some(p) => p,
            None => {
                gl.delete_shader(Some(&vert));
                gl.delete_shader(Some(&frag));
                return Err("Cannot create structure shader program".into());
            }
        };
        gl.attach_shader(&program, &vert);
        gl.attach_shader(&program, &frag);
        gl.link_program(&program);
        gl.delete_shader(Some(&vert));
        gl.delete_shader(Some(&frag));
        if gl
            .get_program_parameter(&program, Gl::LINK_STATUS)
            .as_bool()
            != Some(true)
        {
            let e = gl.get_program_info_log(&program).unwrap_or_default();
            gl.delete_program(Some(&program));
            return Err(e);
        }
        let buffer = match gl.create_buffer() {
            Some(b) => b,
            None => {
                gl.delete_program(Some(&program));
                return Err("Cannot create structure geometry buffer".into());
            }
        };
        gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&buffer));
        let array = js_sys::Float32Array::from(vertices.as_slice());
        gl.buffer_data_with_array_buffer_view(Gl::ARRAY_BUFFER, &array, Gl::STATIC_DRAW);
        Ok(Self {
            gl,
            program,
            buffer,
            points: scene.points.len() as i32,
            lines: (scene.edges.len() * 2) as i32,
        })
    }
    pub fn draw(
        &self,
        canvas: &HtmlCanvasElement,
        camera: Camera,
        picked: Option<u32>,
    ) -> Result<(), String> {
        if self.gl.is_context_lost() {
            return Err("WebGL context lost; reopen the structure view to rebuild it.".into());
        }
        let rect = canvas.get_bounding_client_rect();
        let mut dpr = web_sys::window()
            .map(|w| w.device_pixel_ratio())
            .unwrap_or(1.0)
            .clamp(1.0, 3.0);
        dpr = dpr
            .min(4096.0 / rect.width().max(1.0))
            .min(4096.0 / rect.height().max(1.0));
        let width = (rect.width() * dpr).round().clamp(1.0, 4096.0) as u32;
        let height = (rect.height() * dpr).round().clamp(1.0, 4096.0) as u32;
        if canvas.width() != width {
            canvas.set_width(width);
        }
        if canvas.height() != height {
            canvas.set_height(height);
        }
        let gl = &self.gl;
        gl.viewport(0, 0, width as i32, height as i32);
        gl.clear_color(0.97, 0.98, 0.99, 1.0);
        gl.enable(Gl::DEPTH_TEST);
        gl.clear(Gl::COLOR_BUFFER_BIT | Gl::DEPTH_BUFFER_BIT);
        gl.use_program(Some(&self.program));
        gl.bind_buffer(Gl::ARRAY_BUFFER, Some(&self.buffer));
        for location in [0, 1] {
            gl.enable_vertex_attrib_array(location);
            gl.vertex_attrib_pointer_with_i32(
                location,
                3,
                Gl::FLOAT,
                false,
                28,
                (location * 12) as i32,
            );
        }
        gl.enable_vertex_attrib_array(2);
        gl.vertex_attrib_pointer_with_i32(2, 1, Gl::FLOAT, false, 28, 24);
        let uniform = |name: &str| gl.get_uniform_location(&self.program, name);
        gl.uniform2f(
            uniform("angles").as_ref(),
            camera.yaw as f32,
            camera.pitch as f32,
        );
        gl.uniform2f(
            uniform("pan").as_ref(),
            camera.pan[0] as f32,
            camera.pan[1] as f32,
        );
        gl.uniform1f(uniform("scale").as_ref(), (0.82 * camera.zoom) as f32);
        gl.uniform1f(uniform("aspect").as_ref(), width as f32 / height as f32);
        gl.uniform1f(uniform("point_size").as_ref(), (12.0 * dpr) as f32);
        gl.uniform1f(
            uniform("selected").as_ref(),
            picked.map(|i| i as f32).unwrap_or(-1.0),
        );
        gl.uniform1i(uniform("atoms").as_ref(), 1);
        gl.draw_arrays(Gl::POINTS, 0, self.points);
        gl.uniform1i(uniform("atoms").as_ref(), 0);
        gl.draw_arrays(Gl::LINES, self.points, self.lines);
        Ok(())
    }
}
impl Drop for Renderer {
    fn drop(&mut self) {
        self.gl.delete_buffer(Some(&self.buffer));
        self.gl.delete_program(Some(&self.program));
    }
}
