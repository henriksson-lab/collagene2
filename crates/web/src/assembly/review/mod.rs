//! assembly / review ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_archive_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_finalize;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_output_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_outputs;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_product_export_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_product_exports;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_product_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_product_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_release_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_review;
