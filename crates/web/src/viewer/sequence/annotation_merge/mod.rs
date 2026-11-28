//! Exact-source annotation merge controls and their non-scientific draft guards.
#[cfg(target_arch = "wasm32")]
pub(super) mod attribution;
#[cfg(target_arch = "wasm32")]
pub(super) mod donor;
mod model;
#[cfg(target_arch = "wasm32")]
pub(super) mod view;
