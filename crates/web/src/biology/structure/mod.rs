//! biology / structure ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod structure_camera;
#[cfg(target_arch = "wasm32")]
pub(crate) mod structure_renderer;
#[cfg(target_arch = "wasm32")]
pub(crate) mod structure_view;
