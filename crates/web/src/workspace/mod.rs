//! workspace ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod autosave;
#[cfg(target_arch = "wasm32")]
pub(crate) mod clipboard;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod clipboard_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod document_explorer;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod document_tabs_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod documents;
#[cfg(target_arch = "wasm32")]
pub(crate) mod reference_controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod reference_transfer_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod session;
#[cfg(target_arch = "wasm32")]
pub(crate) mod shell_state;
#[cfg(target_arch = "wasm32")]
pub(crate) mod tabs;
#[cfg(target_arch = "wasm32")]
pub(crate) mod workspace_commit;
#[cfg(target_arch = "wasm32")]
pub(crate) mod workspace_layout;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod workspace_preferences;
#[cfg(target_arch = "wasm32")]
pub(crate) mod workspace_view;
