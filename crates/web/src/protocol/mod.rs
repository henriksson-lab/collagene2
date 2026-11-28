//! protocol ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod plan_editing_host;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod plan_editing_sources;
#[cfg(target_arch = "wasm32")]
pub(crate) mod plan_editing_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod plan_processing_model;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod protocol_layout_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod protocol_layout_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod protocol_quantities_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod protocol_quantities_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod step_description;
