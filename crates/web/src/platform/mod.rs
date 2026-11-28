//! platform ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod commit_lease;
#[cfg(target_arch = "wasm32")]
pub(crate) mod desktop_rpc;
#[cfg(target_arch = "wasm32")]
pub(crate) mod desktop_shell;
#[cfg(target_arch = "wasm32")]
pub(crate) mod desktop_sync_client;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod native_publication_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod native_publication_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod storage;
#[cfg(target_arch = "wasm32")]
pub(crate) mod transient_menus;
#[cfg(target_arch = "wasm32")]
pub(crate) mod worker;
