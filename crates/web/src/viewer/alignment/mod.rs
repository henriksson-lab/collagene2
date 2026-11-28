//! viewer / alignment ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod alignment_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod compare;
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod dotplot;
#[cfg(target_arch = "wasm32")]
pub(crate) mod pairwise_context_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod pairwise_export_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod pairwise_pages;
#[cfg(target_arch = "wasm32")]
pub(crate) mod pairwise_result_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod pairwise_search_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod pairwise_source_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod sequence_similarity_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod sequence_similarity_view;

#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod design_consensus;
