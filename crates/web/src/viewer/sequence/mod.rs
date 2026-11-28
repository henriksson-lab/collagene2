//! viewer / sequence ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod annotation_controller;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod annotation_geometry;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod annotation_layout;
#[cfg(target_arch = "wasm32")]
pub(crate) mod annotation_table;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod annotation_visibility;
#[cfg(target_arch = "wasm32")]
pub(crate) mod editing;
#[cfg(target_arch = "wasm32")]
pub(crate) mod editors;
#[cfg(target_arch = "wasm32")]
pub(crate) mod feature_tools;
#[cfg(target_arch = "wasm32")]
pub(crate) mod search_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod selection;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod selection_toolbar;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod sequence_edit_transition;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod sequence_geometry;
#[cfg(target_arch = "wasm32")]
pub(crate) mod sequence_view;

#[cfg(any(target_arch = "wasm32", test))]
mod annotation_merge;
