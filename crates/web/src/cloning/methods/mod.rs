//! cloning / methods ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod checked_cloning;
#[cfg(target_arch = "wasm32")]
pub(crate) mod crispr;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod domestication_form;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod domestication_review;
#[cfg(target_arch = "wasm32")]
pub(crate) mod domestication_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod multisite;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod multisite_form;
#[cfg(target_arch = "wasm32")]
pub(crate) mod mutagenesis;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod mutation_review;
#[cfg(target_arch = "wasm32")]
pub(crate) mod recombination;
#[cfg(target_arch = "wasm32")]
pub(crate) mod topo_actions_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod topo_graph_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod topo_graph_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod topo_preparation;
#[cfg(target_arch = "wasm32")]
pub(crate) mod topo_view;
