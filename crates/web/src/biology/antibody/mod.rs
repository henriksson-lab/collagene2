//! biology / antibody ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod alignment_antibody_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod alignment_antibody_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod antibody_model;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod antibody_overlay;
#[cfg(target_arch = "wasm32")]
pub(crate) mod antibody_overlay_draw;
#[cfg(target_arch = "wasm32")]
pub(crate) mod antibody_view;
