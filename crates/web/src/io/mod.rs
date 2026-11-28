//! io ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod download;
#[cfg(target_arch = "wasm32")]
pub(crate) mod map_export;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod map_export_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod oligo_import;
#[cfg(target_arch = "wasm32")]
pub(crate) mod sbol_exchange;
#[cfg(target_arch = "wasm32")]
pub(crate) mod sequence_compat_export;
#[cfg(target_arch = "wasm32")]
pub(crate) mod sequence_compat_import;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod sequence_export_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod source_gateway_client;
#[cfg(target_arch = "wasm32")]
pub(crate) mod source_import;
