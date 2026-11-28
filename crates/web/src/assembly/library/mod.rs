//! assembly / library ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_batch_projection;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_batch_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_comparison_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_comparison_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_library_add_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_library_commit;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_library_history_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_library_topology;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_library_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_saved_designs_commit;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_saved_designs_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_saved_designs_view;
