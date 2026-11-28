//! assembly / draft ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_alternatives;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_alternatives_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_draft_actions;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_multisite_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_pool_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_profile_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_role_ids;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_strategy_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_workspace;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod model;
