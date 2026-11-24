//! Camera-relative integer anchors preserve base-level precision on long references.
//! Source indices are at most 10,000,000, hence exactly representable in f32.
#[derive(Debug, Clone, Copy)]
pub struct Camera {
    pub integer: [f32; 2],
    pub fraction: [f32; 2],
    pub span: [f32; 2],
}
pub fn camera(x: f64, y: f64, width: f64, height: f64) -> Camera {
    Camera {
        integer: [x.floor() as f32, y.floor() as f32],
        fraction: [(x - x.floor()) as f32, (y - y.floor()) as f32],
        span: [width as f32, height as f32],
    }
}
pub const VERTEX_SHADER: &str = "#version 300 es\nprecision highp float;in vec3 a_anchor;uniform vec2 u_origin_integer;uniform vec2 u_origin_fraction;uniform vec2 u_span;uniform vec2 u_size;uniform vec2 u_plot;uniform float u_dpr;out vec3 v_color;void main(){vec2 local=(a_anchor.xy-u_origin_integer)+(vec2(0.5)-u_origin_fraction);vec2 p=vec2(52.0,18.0)+local/u_span*u_plot;gl_Position=vec4(p.x/u_size.x*2.0-1.0,1.0-p.y/u_size.y*2.0,0.0,1.0);gl_PointSize=3.0*u_dpr;v_color=a_anchor.z<0.5?vec3(0.17,0.48,0.36):a_anchor.z<1.5?vec3(0.60,0.31,0.53):a_anchor.z<2.5?vec3(0.73,0.52,0.16):vec3(0.82,0.40,0.18);}";
