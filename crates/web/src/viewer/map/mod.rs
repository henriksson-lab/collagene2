//! viewer / map ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod glyph_renderer;
#[cfg(target_arch = "wasm32")]
pub(crate) mod map_labels;
#[cfg(target_arch = "wasm32")]
pub(crate) mod renderer;
