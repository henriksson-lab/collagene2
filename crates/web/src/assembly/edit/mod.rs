//! assembly / edit ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_edit_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_edit_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_fusion_design;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_fusion_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_junction_navigation_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_junctions_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_substitution_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_substitution_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod literal_intent_view;
