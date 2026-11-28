//! library ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod catalog_backup_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod catalog_backup_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod library_bindings;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod library_capture;
#[cfg(target_arch = "wasm32")]
pub(crate) mod library_client;
#[cfg(target_arch = "wasm32")]
pub(crate) mod library_discussion_client;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod library_discussion_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod library_discussion_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod library_documents;
#[cfg(target_arch = "wasm32")]
pub(crate) mod library_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod remote_library_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod remote_removal_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod remote_removal_view;
