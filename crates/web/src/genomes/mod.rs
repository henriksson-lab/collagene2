//! genomes ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod genome_store;
#[cfg(target_arch = "wasm32")]
pub(crate) mod genome_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod genomes;
#[cfg(target_arch = "wasm32")]
pub(crate) mod variant_source_library;
#[cfg(target_arch = "wasm32")]
pub(crate) mod variants_view;
